// Earth Engine Code Editor. See data/covariate_transfer/README.md.
// The AlphaEarth Foundations Satellite Embedding dataset is produced by
// Google and Google DeepMind.
var settings = {
  aoiAsset: 'projects/northstarlabs/assets/aoi_bc',
  localAsset: 'projects/northstarlabs/assets/blue_carbon_cores_1794411291716',
  janousekAsset: 'projects/northstarlabs/assets/janousek_core_locations',
  year: 2024,
  driveFolder: 'blue_carbon_covariate_transfer',
  // Both the local cores and the raster use this same 10 m grid.
  localCrs: 'EPSG:32610',
  rasterBufferM: 100,
  // A global metre-based sampling grid for Janousek cores across UTM zones.
  referenceCrs: 'EPSG:3857',
  localFields: {
    id: 'core_id', latitude: 'latitude', longitude: 'longitude',
    stocks: ['stock_30', 'stock_50', 'stock_100'], stockMultiplier: 1
  },
  referenceFields: {
    id: 'core_id', latitude: 'latitude', longitude: 'longitude',
    stocks: ['carbon_stock_30cm', 'carbon_stock_50cm', 'carbon_stock_100cm'],
    stockMultiplier: 0.1, depth: 'core_depth',
    // The bundled table uses core IDs that can repeat between studies.
    study: 'study_id'
  }
};

if (!settings.aoiAsset || !settings.localAsset || !settings.janousekAsset) {
  throw new Error('Set the three Earth Engine asset IDs at the top of this script.');
}

var noData = -9999;
var transform = [10, 0, 0, 0, -10, 0];
var depths = [30, 50, 100];
var stockColumns = ['0-30', '0-50', '0-100'];
var nativeBands = [];
var bands = [];
var zeroes = [];
for (var i = 0; i < 64; i++) {
  nativeBands.push('A' + (i < 10 ? '0' : '') + i);
  bands.push('var' + (i + 1));
  zeroes.push(0);
}
var selectors = ['Core ID', 'Latitude', 'Longitude'].concat(
  stockColumns, bands, ['embedding_year']
);

function asText(value) {
  var present = ee.Algorithms.If(ee.Algorithms.IsEqual(value, null), '', value);
  var isString = ee.Algorithms.IsEqual(ee.Algorithms.ObjectType(present), 'String');
  var numeric = ee.Number(ee.Algorithms.If(isString, 0, present));
  return ee.String(ee.Algorithms.If(isString, present, numeric.format()))
    .replace('^\\s+|\\s+$', '', 'g');
}

function readNumber(value, prefixPattern) {
  var string = asText(value);
  if (prefixPattern) string = string.replace(prefixPattern, '');
  var pattern = '^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$';
  return ee.Number.parse(ee.String(ee.Algorithms.If(string.match(pattern).size(), string, '-9999')));
}

function sourceCoordinates(row, fields) {
  var longitude = readNumber(row.get(fields.longitude));
  var latitude = readNumber(row.get(fields.latitude));
  var validColumns = longitude.gte(-180).and(longitude.lte(180))
    .and(latitude.gte(-90)).and(latitude.lte(90));

  // CSV assets can store the imported X/Y values in the point geometry.
  var geometry = row.geometry();
  var safeGeometry = ee.Geometry(ee.Algorithms.If(
    ee.Algorithms.IsEqual(geometry, null), ee.Geometry.MultiPoint([]), geometry
  ));
  var geometryCoordinates = ee.Algorithms.If(
    ee.Algorithms.IsEqual(safeGeometry.type(), 'Point'),
    safeGeometry.transform('EPSG:4326', 1).coordinates(), [noData, noData]
  );
  return ee.List(ee.Algorithms.If(validColumns, [longitude, latitude], geometryCoordinates));
}

function standardize(collection, fields) {
  return collection.map(function(row) {
    var id = asText(row.get(fields.id));
    if (fields.study) {
      var study = asText(row.get(fields.study));
      id = ee.String(ee.Algorithms.If(study.length().gt(0).and(id.length().gt(0)),
        study.cat('/').cat(id), id));
    }
    var coordinates = sourceCoordinates(row, fields);
    var record = ee.Feature(null, {
      'Core ID': id,
      'Latitude': coordinates.get(1),
      'Longitude': coordinates.get(0),
      'embedding_year': settings.year
    });
    var measuredDepth = fields.depth ? readNumber(row.get(fields.depth), '^[>~]\\s*') : null;
    stockColumns.forEach(function(column, index) {
      var value = noData;
      if (fields.stocks) {
        var stock = readNumber(row.get(fields.stocks[index]));
        var usable = measuredDepth ? stock.gte(0).and(measuredDepth.gte(depths[index])) : stock.gte(0);
        value = ee.Algorithms.If(usable, stock.multiply(fields.stockMultiplier), noData);
      }
      record = record.set(column, value);
    });
    return record;
  }).filter(ee.Filter.neq('Core ID', ''))
    .filter(ee.Filter.gte('Latitude', -90)).filter(ee.Filter.lte('Latitude', 90))
    .filter(ee.Filter.gte('Longitude', -180)).filter(ee.Filter.lte('Longitude', 180))
    .map(function(row) {
      return row.setGeometry(ee.Geometry.Point([row.get('Longitude'), row.get('Latitude')]));
    });
}

var boundary = ee.FeatureCollection(settings.aoiAsset);
var explicitAoi = boundary.filter(ee.Filter.eq('feature_type', 'aoi'));
var aoi = ee.FeatureCollection(ee.Algorithms.If(explicitAoi.size().gt(0), explicitAoi, boundary)).geometry(1);
var localInput = ee.FeatureCollection(settings.localAsset);
var referenceInput = ee.FeatureCollection(settings.janousekAsset);
var localPoints = standardize(localInput, settings.localFields);
var referencePoints = standardize(referenceInput, settings.referenceFields);
var embeddings = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL');

function embeddingsForYear(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  return embeddings.filterDate(start, start.advance(1, 'year')).select(nativeBands, bands);
}

function mosaic(collection, region) {
  // Supply the band schema even where the collection has no coverage.
  var empty = ee.Image.constant(zeroes).rename(bands).toFloat().updateMask(ee.Image.constant(0));
  var tiles = collection.filterBounds(region).sort('system:index').map(function(tile) {
    return tile.toFloat();
  });
  return ee.ImageCollection.fromImages([empty]).merge(tiles).mosaic();
}

function samplePoints(image, points, crs) {
  // Retain original point properties and rows when the embedding is masked.
  return image.unmask(noData, false).reduceRegions({
    collection: points, reducer: ee.Reducer.first(), crs: crs,
    crsTransform: transform, tileScale: 4
  });
}

var exportPrefix = 'alphaearth_' + settings.year + '_';

function exportTable(table, suffix) {
  var name = exportPrefix + suffix;
  Export.table.toDrive({
    collection: table, description: name, folder: settings.driveFolder,
    fileNamePrefix: name, fileFormat: 'CSV', selectors: selectors
  });
}

var annual = embeddingsForYear(settings.year);
var localImage = mosaic(annual, aoi).clip(aoi);
var referenceImage = mosaic(annual, referencePoints.geometry(1));
var localSamples = samplePoints(localImage, localPoints, settings.localCrs);
var referenceSamples = samplePoints(referenceImage, referencePoints, settings.referenceCrs);
exportTable(localSamples, 'local_cores');
exportTable(referenceSamples, 'janousek_cores');

var rasterName = exportPrefix + 'aoi';
// Without a buffer projection, Earth Engine interprets the distance in metres.
var rasterAoi = settings.rasterBufferM === 0 ? aoi : aoi.buffer({distance: settings.rasterBufferM});
var rasterImage = mosaic(annual, rasterAoi).clip(rasterAoi);
Export.image.toDrive({
  image: rasterImage.unmask(noData, false), description: rasterName,
  folder: settings.driveFolder, fileNamePrefix: rasterName,
  region: rasterAoi, crs: settings.localCrs, crsTransform: transform,
  maxPixels: 1e10, fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true, noData: noData}
});

print('Embedding year', settings.year);
print('Local input example', localInput.first());
print('Reference input example', referenceInput.first());
print('Local rows, rows with coordinates and IDs', localInput.size(), localPoints.size());
print('Reference rows, rows with coordinates and IDs', referenceInput.size(), referencePoints.size());
print('Annual embedding tiles intersecting the AOI', annual.filterBounds(aoi).size());
print('GeoTIFF buffer (m)', settings.rasterBufferM);
print('Local sample', localSamples.first());
var referencePreviewPoints = referencePoints.limit(1);
var referencePreviewImage = mosaic(annual, referencePreviewPoints.geometry(1));
print('Reference sample', samplePoints(referencePreviewImage, referencePreviewPoints, settings.referenceCrs).first());

var comparisonYear = settings.year - 1;
var previousImage = mosaic(embeddingsForYear(comparisonYear), aoi).clip(aoi);
// A01, A16 and A09 are renamed var2, var17 and var10 in the mosaics.
var rgb = {min: -0.3, max: 0.3, bands: ['var2', 'var17', 'var10']};
var similarity = localImage.multiply(previousImage).reduce(ee.Reducer.sum());
var years = comparisonYear + '–' + settings.year;
Map.setOptions('SATELLITE');
Map.centerObject(aoi);
Map.addLayer(previousImage, rgb, comparisonYear + ' embeddings', false);
Map.addLayer(localImage, rgb, settings.year + ' embeddings');
Map.addLayer(similarity, {min: 0, max: 1, palette: ['white', 'black']},
  years + ' similarity (brighter = less similar)', false);
Map.addLayer(ee.FeatureCollection([ee.Feature(aoi)]).style({
  color: 'f4bc56', fillColor: '00000000', width: 2
}), {}, 'AOI');
Map.addLayer(referencePoints.style({
  color: '202020', fillColor: 'ff79c6', pointSize: 5, width: 1
}), {}, 'Janousek cores');
Map.addLayer(localPoints.style({
  color: '202020', fillColor: '63e6ff', pointSize: 8, width: 1
}), {}, 'Local cores');
