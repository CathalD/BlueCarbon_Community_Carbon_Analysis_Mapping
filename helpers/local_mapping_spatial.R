# Shared AOI, core locations and raster grid for the local and Bayesian maps.

read_mapping_inputs <- function(aoi_file, mapping_file, analysis_crs, cell_size_m) {
  if (!file.exists(aoi_file)) stop("AOI file not found: ", aoi_file)
  if (!file.exists(mapping_file)) {
    stop("Mapping data not found: ", mapping_file, ". Run script_1.R first.")
  }
  if (length(cell_size_m) != 1L || !is.numeric(cell_size_m) ||
      !is.finite(cell_size_m) || cell_size_m <= 0) {
    stop("cell_size_m must be a positive number of metres.")
  }
  crs <- sf::st_crs(analysis_crs)
  if (is.na(crs) || !grepl("^(PROJCRS|PROJCS)\\[", crs$wkt) ||
      !tolower(crs$units_gdal) %in% c("metre", "meter", "m")) {
    stop("analysis_crs must be a projected CRS in metres.")
  }

  boundary <- sf::st_read(aoi_file, quiet = TRUE)
  if (!nrow(boundary) || is.na(sf::st_crs(boundary)) ||
      any(!sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("The AOI file must contain polygon features with a coordinate reference system.")
  }
  boundary <- sf::st_make_valid(sf::st_transform(boundary, crs))
  geometry <- sf::st_geometry(boundary)
  feature_area <- as.numeric(sf::st_area(geometry))
  if (any(!is.finite(feature_area) | feature_area <= 0)) {
    stop("Every AOI file feature must have positive polygon area.")
  }

  feature_type <- if ("feature_type" %in% names(boundary)) {
    tolower(trimws(as.character(boundary$feature_type)))
  } else rep(NA_character_, nrow(boundary))
  explicit_aoi <- feature_type %in% "aoi"
  aoi_geometry <- sf::st_union(geometry[if (any(explicit_aoi)) explicit_aoi else rep(TRUE, nrow(boundary))])
  aoi <- sf::st_sf(geometry = aoi_geometry)
  aoi_area <- as.numeric(sf::st_area(aoi))

  zone_labels <- if ("zone" %in% names(boundary)) {
    trimws(as.character(boundary$zone))
  } else rep(NA_character_, nrow(boundary))
  has_label <- !is.na(zone_labels) & nzchar(zone_labels)
  if (any(feature_type %in% "zone" & !has_label)) {
    stop("Every feature marked as a zone needs a zone name.")
  }
  is_zone <- has_label & !explicit_aoi
  zone_names <- sort(unique(zone_labels[is_zone]))
  zone_geometry <- sf::st_sfc(lapply(zone_names, function(name) {
    same_zone <- is_zone & zone_labels %in% name
    clipped <- sf::st_intersection(sf::st_union(geometry[same_zone]), aoi_geometry)
    if (!length(clipped) || as.numeric(sf::st_area(clipped)) <= 0) {
      stop("Zone is outside the AOI: ", name)
    }
    clipped[[1]]
  }), crs = crs)
  neighbours <- sf::st_intersects(zone_geometry)
  # Projecting a shared edge with different vertex spacing can produce slivers.
  overlap_tolerance_m2 <- min(1, 1e-8 * aoi_area)
  for (i in seq_along(neighbours)) {
    for (j in neighbours[[i]][neighbours[[i]] > i]) {
      overlap <- sf::st_intersection(zone_geometry[i], zone_geometry[j])
      overlap_area <- sum(as.numeric(sf::st_area(overlap)))
      if (overlap_area > overlap_tolerance_m2) {
        stop("Zones overlap inside the AOI: ", zone_names[i], " and ", zone_names[j])
      }
      if (overlap_area > 0) {
        zone_geometry[j] <- sf::st_difference(zone_geometry[j], zone_geometry[i])
      }
    }
  }
  zones <- sf::st_sf(zone = zone_names, geometry = zone_geometry)
  zones$area_m2 <- as.numeric(sf::st_area(zones))
  if (any(zones$area_m2 <= 0)) stop("Each zone needs positive area inside the AOI.")

  records <- utils::read.csv(mapping_file, check.names = FALSE,
                             colClasses = c("Core ID" = "character", "Plot ID" = "character"))
  required <- c("Core ID", "Plot ID", "Latitude", "Longitude")
  missing <- setdiff(required, names(records))
  if (length(missing) || anyDuplicated(names(records))) {
    stop("Mapping data needs unique columns including ", paste(required, collapse = ", "), ".")
  }
  if (!nrow(records)) stop("Mapping data has no cores.")
  for (name in c("Core ID", "Plot ID")) {
    records[[name]] <- trimws(records[[name]])
    if (anyNA(records[[name]]) || any(!nzchar(records[[name]]))) {
      stop("Mapping data contains missing values in ", name, ".")
    }
  }
  if (anyDuplicated(records[["Core ID"]])) stop("Mapping data contains duplicate Core IDs.")
  records$Latitude <- suppressWarnings(as.numeric(records$Latitude))
  records$Longitude <- suppressWarnings(as.numeric(records$Longitude))
  invalid <- !is.finite(records$Latitude) | abs(records$Latitude) > 90 |
    !is.finite(records$Longitude) | abs(records$Longitude) > 180
  if (any(invalid)) {
    stop("Invalid WGS84 coordinates for cores: ", paste(records[["Core ID"]][invalid], collapse = ", "))
  }
  cores <- sf::st_as_sf(records, coords = c("Longitude", "Latitude"),
                        crs = 4326, remove = FALSE)
  cores <- sf::st_transform(cores, crs)
  outside <- lengths(sf::st_intersects(cores, aoi)) == 0L
  if (any(outside)) {
    stop("Core locations outside the AOI: ", paste(cores[["Core ID"]][outside], collapse = ", "),
         ". Check that the AOI matches the mapping data.")
  }

  cores$zone_assigned <- rep(NA_character_, nrow(cores))
  if (nrow(zones)) {
    recorded <- if ("zone" %in% names(cores)) trimws(as.character(cores$zone)) else {
      rep(NA_character_, nrow(cores))
    }
    matches <- sf::st_intersects(cores, zones)
    for (i in seq_along(matches)) {
      hit <- zones$zone[matches[[i]]]
      declared <- !is.na(recorded[i]) && nzchar(recorded[i])
      if (declared && !recorded[i] %in% hit) {
        stop("Recorded zone ", recorded[i], " disagrees with geometry for core ",
             cores[["Core ID"]][i], ".")
      }
      if (length(hit) > 1L && !declared) {
        stop("Core ", cores[["Core ID"]][i],
             " lies on a shared zone boundary. Supply its zone in the mapping data.")
      }
      if (length(hit)) cores$zone_assigned[i] <- if (declared) recorded[i] else hit[1L]
    }
  }

  # Keep a common origin and use raster cell indices even when edge cells are clipped.
  bounds <- sf::st_bbox(aoi)
  lower <- floor(bounds[c("xmin", "ymin")] / cell_size_m) * cell_size_m
  upper <- ceiling(bounds[c("xmax", "ymax")] / cell_size_m) * cell_size_m
  dimensions <- as.integer(round((upper - lower) / cell_size_m))
  template <- terra::rast(ncols = dimensions[1L], nrows = dimensions[2L],
                          xmin = lower[1L], xmax = upper[1L],
                          ymin = lower[2L], ymax = upper[2L], crs = crs$wkt)
  grid_geometry <- sf::st_make_grid(aoi, cellsize = cell_size_m,
                                   offset = unname(lower), n = dimensions)
  candidates <- grid_geometry[lengths(sf::st_intersects(grid_geometry, aoi_geometry)) > 0L]
  centers <- sf::st_centroid(candidates)
  cell_id <- terra::cellFromXY(template, sf::st_coordinates(centers))
  grid <- sf::st_sf(cell_id = cell_id, geometry = candidates, agr = "constant")
  cells <- sf::st_intersection(grid, aoi_geometry)
  cells$area_m2 <- as.numeric(sf::st_area(cells))
  cells <- cells[cells$area_m2 > 0, ]
  if (any(sf::st_geometry_type(cells) == "GEOMETRYCOLLECTION")) {
    cells <- suppressWarnings(sf::st_collection_extract(cells, "POLYGON"))
  }
  cells <- cells[order(cells$cell_id), ]
  if (!nrow(cells) || anyDuplicated(cells$cell_id) ||
      abs(sum(cells$area_m2) - aoi_area) > max(1e-4, aoi_area * 1e-8)) {
    stop("The clipped grid does not reproduce the AOI area.")
  }

  centers <- centers[match(cells$cell_id, cell_id)]
  off_aoi <- lengths(sf::st_intersects(centers, aoi_geometry)) == 0L
  if (any(off_aoi)) {
    centers[off_aoi] <- sf::st_point_on_surface(sf::st_geometry(cells)[off_aoi])
  }
  points <- sf::st_sf(cell_id = cells$cell_id, geometry = centers)
  list(aoi = aoi, zones = zones, cores = cores, cells = cells, points = points,
       template = template, area_m2 = aoi_area)
}

export_mapping_raster <- function(inputs, values, path) {
  if (!is.data.frame(values) || !ncol(values) || nrow(values) != nrow(inputs$cells) ||
      anyNA(names(values)) || any(!nzchar(names(values))) || anyDuplicated(names(values)) ||
      !all(vapply(values, is.numeric, logical(1)))) {
    stop("Raster values must have one row per grid cell and a name for each numeric band.")
  }
  raster <- terra::rast(inputs$template, nlyrs = ncol(values))
  pixels <- matrix(NA_real_, nrow = terra::ncell(raster), ncol = ncol(values))
  pixels[inputs$cells$cell_id, ] <- as.matrix(values)
  terra::values(raster) <- pixels
  names(raster) <- names(values)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(terra::writeRaster(raster, path, datatype = "FLT8S", overwrite = TRUE))
}
