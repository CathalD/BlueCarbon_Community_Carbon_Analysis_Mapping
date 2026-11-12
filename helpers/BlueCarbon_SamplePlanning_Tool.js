// =================================================================================
// BLUE CARBON SAMPLING DESIGN TOOL — 2026
// WWF-Canada · Coastal Blue Carbon Hub
//
// Companion to Part 2 (Project Planning) of the Blue Carbon Eelgrass Workshop.
// Every number this tool reports also appears in the workshop's Sample Allocation
// Calculator; Section 4 below proves they agree.
//
//   SECTION 0   Configuration and resource links
//   SECTION 1   Statistics core          (BCStats)
//   SECTION 2   Geometry and feasibility (BCGeom)
//   SECTION 3   Test harness             (BCTest)
//   SECTION 4   Run the tests
//
// To run the tests: set RUN_SELF_TEST to true and press Run. Results print to
// the Console panel. Set it to false for the released app.
// =================================================================================

var RUN_SELF_TEST = true;

// Works in the GEE Code Editor (print) and in plain JavaScript (console.log).
var LOG = (typeof print === 'function') ? print : console.log;


// =================================================================================
// === SECTION 0 — CONFIGURATION ===================================================
// =================================================================================

var CONFIG = {

  VERSION: '2026.1',

  // --- Defaults, matched to the workshop text -----------------------------------
  PLOT_M2:            100,    // 10 x 10 m plot per core (Appendix A3)
  CONFIDENCE:         0.90,   // Step 4
  MARGIN_OF_ERROR:    0.20,   // Step 4
  MIN_CORES_PER_ZONE: 5,      // Appendix A7
  ALLOCATION:         'proportional',
  DEFAULT_ECOSYSTEM:  'Eelgrass',
  DEFAULT_DEPTH_CM:   30,
  SEED:               42,

  // --- Resource links shown to the user -----------------------------------------
  // TODO: replace REPO with the real path once the workshop repo is public.
  REPO: 'https://github.com/WWF-Canada-SKI/BlueCarbon_Eelgrass_Workshop',

  LINKS: {
    calculator: '/raw/main/02_Project_Planning/BlueCarbon_SampleAllocation_2026.xlsx',
    planningGuide: '/blob/main/02_Project_Planning/README.md',
    appendixA: '/blob/main/02_Project_Planning/README.md#appendix-a--a-brief-lesson-in-sampling-logic',
    fieldGuide: '/blob/main/Coastal-Blue-Carbon-Field-Guide-FINAL.pdf',
    // TODO: point these at the real files once uploaded.
    datasheets: '/blob/main/03_Field_Methods/',
    workedExample: '/blob/main/Worked_Example/02_Project_Planning.md'
  },

  linkTo: function (key) { return this.REPO + this.LINKS[key]; }
};


// =================================================================================
// === SECTION 1 — STATISTICS CORE =================================================
// ===
// === Every reported number originates here. Verified against
// === BlueCarbon_SampleAllocation_2026.xlsx by the tests in Section 3.
// =================================================================================

var BCStats = {

  VERSION: '2026.1',

  // --- Inverse normal -----------------------------------------------------------
  // Replaces the cubic polynomial used previously, which under-estimated z by 10%
  // at 98% confidence and 47% at 99.9%. Since n scales with z squared, those were
  // 20% and 118% errors in the core count.
  // Acklam's algorithm; relative error below 1.15e-9.

  invNorm: function (p) {
    if (!(p > 0 && p < 1)) return NaN;

    var a = [-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00];
    var b = [-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01];
    var c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00];
    var d = [ 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
              3.754408661907416e+00];

    var pLow = 0.02425, pHigh = 1 - pLow, q, r;

    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    }
    if (p <= pHigh) {
      q = p - 0.5; r = q * q;
      return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
             (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
    }
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
  },

  z: function (confidence) { return this.invNorm(1 - (1 - confidence) / 2); },

  // --- Population size ----------------------------------------------------------
  // Each core represents a plot, not a point (Appendix A3).

  populationSize: function (areaM2, plotM2) {
    if (!areaM2 || !plotM2 || plotM2 <= 0) return NaN;
    return Math.floor(areaM2 / plotM2);
  },

  // --- Sample size, one uniform area --------------------------------------------
  //   n = z^2 N CV^2 / ( (N-1) E^2 + z^2 CV^2 )       (Appendix A3)

  srsSampleSize: function (o) {
    var z  = this.z(o.confidence);
    var N  = this.populationSize(o.areaM2, o.plotM2);
    var cv = o.sd / o.mean;
    var E  = o.marginOfError;

    if (!isFinite(N) || N < 2 || !isFinite(cv) || cv <= 0) {
      return { ok: false, reason: 'Need a valid area, plot size, mean and standard deviation.' };
    }

    var n = (z*z * N * cv*cv) / ((N - 1) * E*E + z*z * cv*cv);

    return { ok: true, n: Math.ceil(n), nExact: n, z: z, N: N, cv: cv,
             samplingFraction: Math.ceil(n) / N, inputs: o };
  },

  // --- Sample size, stratified --------------------------------------------------
  // The allocation choice drives BOTH the formula and the split. Mixing them is
  // the error in the earlier tools: they sized with the Neyman formula and then
  // allocated by area, which under-samples because
  //   (sum Ni si)^2  <=  N sum(Ni si^2)      [Cauchy-Schwarz]
  //
  //   proportional   n = z^2 N V / ((N-1) E^2 + z^2 V),  V = pooled var / mean^2
  //                  split proportional to area          (Appendix A7)
  //   neyman         n = (sum Ni si)^2 / ((N Eabs / z)^2 + sum Ni si^2)
  //                  split proportional to Ni si
  //
  // strata: [{ name, areaM2, mean, sd }, ...]

  stratifiedSampleSize: function (o) {
    var strata = [], i;
    for (i = 0; i < o.strata.length; i++) {
      var s0 = o.strata[i];
      if (s0 && s0.areaM2 > 0 && s0.mean > 0 && s0.sd >= 0) strata.push(s0);
    }
    if (strata.length === 0) {
      return { ok: false, reason: 'No zone has an area, mean and standard deviation.' };
    }

    var allocation = o.allocation || 'proportional';
    var minPer     = (o.minPerStratum === undefined) ? CONFIG.MIN_CORES_PER_ZONE : o.minPerStratum;
    var z          = this.z(o.confidence);
    var E          = o.marginOfError;
    var plot       = o.plotM2;
    var j;

    var G = 0;
    for (j = 0; j < strata.length; j++) G += strata[j].areaM2;
    var N = this.populationSize(G, plot);

    var overallMean = 0, pooledVar = 0, sumNiSi = 0, sumNiSi2 = 0;
    for (j = 0; j < strata.length; j++) {
      var st = strata[j];
      var w  = st.areaM2 / G;
      var Ni = this.populationSize(st.areaM2, plot);
      st._Ni = Ni;
      overallMean += w * st.mean;
      pooledVar   += w * st.sd * st.sd;
      sumNiSi     += Ni * st.sd;
      sumNiSi2    += Ni * st.sd * st.sd;
    }

    var Eabs = overallMean * E;
    var V    = pooledVar / (overallMean * overallMean);
    var nExact;

    if (allocation === 'neyman') {
      nExact = (sumNiSi * sumNiSi) / (Math.pow(N * Eabs / z, 2) + sumNiSi2);
    } else {
      nExact = (z*z * N * V) / ((N - 1) * E*E + z*z * V);
    }

    var n = Math.ceil(nExact);
    var allocated = [], total = 0;

    for (j = 0; j < strata.length; j++) {
      var stx   = strata[j];
      var share = (allocation === 'neyman')
                    ? (stx._Ni * stx.sd) / sumNiSi
                    : stx.areaM2 / G;
      var raw   = n * share;
      var cores = Math.max(Math.ceil(raw), minPer);
      total += cores;
      allocated.push({ name: stx.name, code: stx.code, areaM2: stx.areaM2, mean: stx.mean, sd: stx.sd,
                       cv: stx.sd / stx.mean, Ni: stx._Ni, share: share,
                       rawCores: raw, cores: cores,
                       flooredToMinimum: Math.ceil(raw) < minPer });
    }

    return { ok: true, allocation: allocation, n: n, nExact: nExact, nAllocated: total,
             z: z, N: N, totalAreaM2: G, overallMean: overallMean,
             pooledVariance: pooledVar, V: V, marginOfErrorAbsolute: Eabs,
             strata: allocated };
  },

  // --- Achieved precision, after fieldwork (Appendix A8) ------------------------
  // z throughout, matching the workbook and the UNFCCC A6.4 tool.

  achievedPrecisionSRS: function (o) {
    var z   = this.z(o.confidence);
    var se  = Math.sqrt((1 - o.n / o.N) * (o.sampleSd * o.sampleSd) / o.n);
    var rme = z * se / o.sampleMean;
    return { ok: true, z: z, se: se, rme: rme, target: o.target, pass: rme <= o.target };
  },

  achievedPrecisionStratified: function (o) {
    var G = 0, j;
    for (j = 0; j < o.strata.length; j++) G += o.strata[j].areaM2;

    var z = this.z(o.confidence), overallMean = 0, varSum = 0;
    for (j = 0; j < o.strata.length; j++) {
      var s = o.strata[j];
      var w = s.areaM2 / G;
      var Nh = this.populationSize(s.areaM2, o.plotM2);
      overallMean += w * s.sampleMean;
      varSum += w * w * (1 - s.cores / Nh) * (s.sampleSd * s.sampleSd) / s.cores;
    }

    var se = Math.sqrt(varSum), rme = z * se / overallMean;
    return { ok: true, z: z, se: se, rme: rme, overallMean: overallMean,
             target: o.target, pass: rme <= o.target };
  },

  // --- Priors, Janousek et al. (2025) -------------------------------------------
  // One figure per ecosystem, pooled across the Pacific Northwest: Alaska,
  // British Columbia, Washington and Oregon. California and south are excluded —
  // stocks there behave differently and the region is not comparable.
  //
  // Computed from core-level values in the published dataset, not by averaging
  // published regional means, so the SD is the real spread between cores.
  //
  // Mangrove is absent: it does not occur north of California.
  //
  // Janousek, C.N., Krause, J.R., Drexler, J.Z., Buffington, K.J., Poppe, K.L.,
  // Peck, E., et al. (2025). Blue carbon stocks along the Pacific coast of North
  // America are mainly driven by local rather than regional factors.
  // Global Biogeochemical Cycles, 39, e2024GB008239.

  //                    0             1      2     3      4      5     6
  //                 ecosystem      m30   sd30  n30   m100  sd100  n100
  PNW_PRIORS: [
    ['Eelgrass',                    24.8, 16.8, 175,  86.5,  41.7,  42],
    ['Salt marsh',                  88.4, 36.9, 351, 229.4,  89.4, 157],
    ['Tidal swamp',                112.5, 39.3,  62, 356.9, 104.9,  39],
    ['Tideflat',                    26.1, 13.8,  55,  null,  null,   9]
  ],

  // Priors above this CV get a warning: the campaign will be unusually large.
  HIGH_CV: 0.80,

  // Below this the formula's answer is not usable in practice. One core has no
  // standard deviation, so the achieved precision in Appendix A8 cannot be
  // computed at all, and the estimate can never be checked. Surfaced to the user
  // rather than applied silently — the earlier tools floored at 10 with no notice.
  MIN_USABLE_CORES: 5,

  // Fewer cores than this behind a prior makes it indicative only.
  THIN_EVIDENCE: 20,

  regionalPrior: function (ecosystem, depthCm) {
    var deep = (depthCm === 100);
    for (var i = 0; i < this.PNW_PRIORS.length; i++) {
      var p = this.PNW_PRIORS[i];
      if (p[0] !== ecosystem) continue;
      var mean  = deep ? p[4] : p[1];
      var sd    = deep ? p[5] : p[2];
      var cores = deep ? p[6] : p[3];
      if (mean === null || sd === null || cores < 3) {
        return { ok: false, ecosystem: ecosystem, depthCm: depthCm, cores: cores,
                 reason: 'Only ' + cores + ' published cores reach ' + depthCm +
                         ' cm in ' + ecosystem.toLowerCase() + ' here — not enough to ' +
                         'size a campaign. Use the top 30 cm, or your own pilot data.' };
      }
      return { ok: true, ecosystem: ecosystem, depthCm: depthCm,
               mean: mean, sd: sd, cv: sd / mean, cores: cores,
               source: 'Pacific Northwest average, Janousek et al. (2025)',
               indicative: cores < this.THIN_EVIDENCE,
               highVariability: (sd / mean) > this.HIGH_CV };
    }
    return { ok: false, reason: 'No published data for ' + ecosystem + '.' };
  },

  // Only ecosystems that occur north of California.
  ecosystems: function () {
    var out = [];
    for (var i = 0; i < this.PNW_PRIORS.length; i++) out.push(this.PNW_PRIORS[i][0]);
    return out;
  }
};


// =================================================================================
// === SECTION 2 — GEOMETRY AND FEASIBILITY ========================================
// ===
// === Fixes the failure that made the earlier tools unusable on real eelgrass
// === meadows: analysis scale was fixed at 250 m, so a 5 ha site resolved to zero
// === pixels and stratum areas came back empty.
// =================================================================================

var BCGeom = {

  VERSION: '2026.1',

  SCALE_LADDER: [10, 20, 30, 50, 100, 250],  // never finer than Sentinel-2
  TARGET_PIXELS: 5000,
  BUFFER_LADDER: [50, 25, 10, 5, 0],
  MIN_AREA_RETAINED: 0.60,

  // Square packing is a theoretical ceiling. Placing cores at random with a
  // minimum separation stalls well below it.
  FILL_WARN:  0.25,
  FILL_LIMIT: 0.50,

  // --- Analysis scale -----------------------------------------------------------

  analysisScale: function (areaM2) {
    if (!(areaM2 > 0)) return { ok: false, reason: 'Area must be greater than zero.' };

    var ideal = Math.sqrt(areaM2 / this.TARGET_PIXELS);
    var scale = this.SCALE_LADDER[0];
    for (var i = 0; i < this.SCALE_LADDER.length; i++) {
      if (this.SCALE_LADDER[i] <= ideal) scale = this.SCALE_LADDER[i];
    }
    var pixels = Math.floor(areaM2 / (scale * scale));

    return { ok: true, scale: scale, pixels: pixels, areaHa: areaM2 / 10000,
             coarse: pixels < 100,
             note: 'Analysis at ' + scale + ' m gives about ' + pixels +
                   ' pixels across the site.' };
  },

  // --- Which stratification methods can resolve this site -----------------------

  SOURCES: [
    { id: 'draw',       label: 'Draw them by hand',                    nativeScale: null, minPixels: 0 },
    { id: 'upload',     label: 'Upload my own boundaries',             nativeScale: null, minPixels: 0 },
    { id: 'dynamic',    label: 'Land cover (Dynamic World, 10 m)',     nativeScale: 10,  minPixels: 200 },
    { id: 'copernicus', label: 'Land cover (Copernicus, 100 m)',       nativeScale: 100, minPixels: 200 },
    { id: 'covariates', label: 'Satellite imagery grouping (30 m)',    nativeScale: 30,  minPixels: 300 },
    { id: 'embeddings', label: 'Satellite Embeddings grouping (10 m)', nativeScale: 10,  minPixels: 300 }
  ],

  stratificationOptions: function (areaM2) {
    var out = [];
    for (var i = 0; i < this.SOURCES.length; i++) {
      var s = this.SOURCES[i];
      var entry = { id: s.id, label: s.label, available: true, reason: '' };
      if (s.nativeScale) {
        var px = Math.floor(areaM2 / (s.nativeScale * s.nativeScale));
        entry.pixels = px;
        if (px < s.minPixels) {
          entry.available = false;
          entry.reason = 'Only ' + px + ' pixels fit inside this site at ' +
                         s.nativeScale + ' m. Too coarse to split it up.';
        }
      }
      out.push(entry);
    }
    return out;
  },

  // --- Core spacing -------------------------------------------------------------
  // Two cores 3 m apart both claim a 100 m2 plot, so they are one observation
  // counted twice. Cores are kept at least one plot width apart.

  minSpacing: function (plotM2) { return Math.sqrt(plotM2); },

  capacity: function (areaM2, plotM2) {
    var s = this.minSpacing(plotM2);
    return Math.floor(areaM2 / (s * s));
  },

  practicalCapacity: function (areaM2, plotM2) {
    return Math.floor(this.capacity(areaM2, plotM2) * this.FILL_LIMIT);
  },

  // --- Inward buffer ------------------------------------------------------------
  // Client-side estimate assuming a compact shape. The Earth Engine layer verifies
  // against real geometry and steps down again for long thin strata.

  retainedAfterBuffer: function (areaM2, bufferM) {
    var side = Math.sqrt(areaM2) - 2 * bufferM;
    return side <= 0 ? 0 : side * side;
  },

  chooseBuffer: function (areaM2, coresNeeded, plotM2) {
    for (var i = 0; i < this.BUFFER_LADDER.length; i++) {
      var b = this.BUFFER_LADDER[i];
      var retained = this.retainedAfterBuffer(areaM2, b);
      if (retained / areaM2 < this.MIN_AREA_RETAINED) continue;
      if (this.practicalCapacity(retained, plotM2) < coresNeeded) continue;
      return { ok: true, buffer: b, retainedM2: retained,
               retainedFraction: retained / areaM2,
               note: b === 0
                 ? 'No edge buffer — the zone is too small to give any up.'
                 : b + ' m edge buffer, keeping ' +
                   Math.round(100 * retained / areaM2) + '% of the zone.' };
    }
    return { ok: false, buffer: 0, retainedM2: areaM2, retainedFraction: 1,
             note: 'No buffer possible. This zone is barely larger than the cores it must hold.' };
  },

  // --- Feasibility --------------------------------------------------------------
  // Run before any point generation. Every message names the fix.

  feasibility: function (o) {
    var plot   = o.plotM2;
    var minPer = (o.minPerStratum === undefined) ? CONFIG.MIN_CORES_PER_ZONE : o.minPerStratum;
    var strata = o.strata || [{ name: 'Whole site', areaM2: o.areaM2, cores: o.cores }];
    var problems = [], warnings = [], j;

    var total = 0;
    for (j = 0; j < strata.length; j++) total += strata[j].areaM2;
    var N = Math.floor(total / plot);

    if (N < 1) {
      problems.push('The site (' + Math.round(total) + ' m2) is smaller than one ' + plot +
                    ' m2 plot. Check the boundary, or use a smaller plot size.');
      return { ok: false, problems: problems, warnings: warnings, N: N };
    }

    if (strata.length * minPer > N) {
      problems.push(strata.length + ' zones at a minimum of ' + minPer + ' cores each needs ' +
                    (strata.length * minPer) + ' plots, but the site only holds ' + N +
                    '. Use fewer zones or a smaller plot size.');
    }

    for (j = 0; j < strata.length; j++) {
      var s    = strata[j];
      var name = s.name || ('Zone ' + (j + 1));
      var Nh   = Math.floor(s.areaM2 / plot);

      if (Nh < minPer) {
        problems.push(name + ' holds only ' + Nh + ' plots but needs at least ' + minPer +
                      '. Merge it into a neighbouring zone, or drop it.');
        continue;
      }

      var cores = s.cores || minPer;
      var buf   = this.chooseBuffer(s.areaM2, cores, plot);
      var slots = this.capacity(buf.retainedM2, plot);
      var fill  = slots > 0 ? cores / slots : Infinity;

      if (fill > this.FILL_LIMIT) {
        problems.push(name + ' cannot hold ' + cores + ' cores kept ' + this.minSpacing(plot) +
                      ' m apart. About ' + Math.floor(slots * this.FILL_LIMIT) +
                      ' is the practical limit for this zone. Reduce the core count, ' +
                      'enlarge the zone, or merge it.');
      } else if (fill > this.FILL_WARN) {
        warnings.push(name + ' will be densely sampled (' + cores + ' cores in about ' +
                      slots + ' possible positions). Placement may take several attempts.');
      }

      if (!buf.ok) {
        warnings.push(name + ' is too small for an edge buffer, so cores may land near its boundary.');
      }
    }

    var sc = this.analysisScale(total);
    if (sc.ok && sc.coarse) {
      warnings.push('This site is small (' + sc.areaHa.toFixed(1) + ' ha). Areas are computed at ' +
                    sc.scale + ' m, the finest available, giving about ' + sc.pixels +
                    ' pixels. Expect the area to be approximate.');
    }

    return { ok: problems.length === 0, problems: problems, warnings: warnings,
             N: N, scale: sc.scale, totalAreaM2: total };
  }
};


// =================================================================================
// === SECTION 3 — TEST HARNESS ====================================================
// ===
// === Checks the two sections above against BlueCarbon_SampleAllocation_2026.xlsx
// === and against hand calculations. Pure JavaScript; no Earth Engine calls, so it
// === runs instantly and needs no assets.
// =================================================================================

var BCTest = {

  lines: [], passed: 0, failed: 0,

  // ---- helpers -----------------------------------------------------------------
  say:  function (t) { this.lines.push(t); },
  head: function (t) { this.say(''); this.say('=== ' + t + ' ==='); },
  pad:  function (s, n, right) {
    s = String(s);
    while (s.length < n) { s = right ? (' ' + s) : (s + ' '); }
    return s;
  },
  near: function (label, got, want, tol) {
    tol = (tol === undefined) ? 1e-6 : tol;
    var ok = isFinite(got) && Math.abs(got - want) <= tol;
    ok ? this.passed++ : this.failed++;
    this.say('  ' + (ok ? 'PASS  ' : 'FAIL  ') + this.pad(label, 34) +
             (ok ? '' : '  got ' + got + ', want ' + want));
    return ok;
  },
  eq: function (label, got, want) {
    var ok = (got === want);
    ok ? this.passed++ : this.failed++;
    this.say('  ' + (ok ? 'PASS  ' : 'FAIL  ') + this.pad(label, 34) +
             (ok ? '' : '  got ' + got + ', want ' + want));
    return ok;
  },

  // ---- 3.1 inverse normal ------------------------------------------------------
  testInverseNormal: function () {
    this.head('3.1  Inverse normal, against published z values');
    this.near('z at 80%',   BCStats.z(0.80),  1.2815516);
    this.near('z at 90%',   BCStats.z(0.90),  1.6448536);
    this.near('z at 95%',   BCStats.z(0.95),  1.9599640);
    this.near('z at 98%',   BCStats.z(0.98),  2.3263479);
    this.near('z at 99%',   BCStats.z(0.99),  2.5758293);
    this.near('z at 99.9%', BCStats.z(0.999), 3.2905267);

    this.say('');
    this.say('  For comparison, the polynomial the old tools used:');
    var truth = [[0.90, 1.6448536], [0.95, 1.9599640], [0.99, 2.5758293], [0.999, 3.2905267]];
    for (var i = 0; i < truth.length; i++) {
      var conf = truth[i][0], t = truth[i][1], al = 1 - conf;
      var old = 2.41 - 10.9*al + 37.7*al*al - 57.9*al*al*al;
      this.say('    ' + this.pad((conf*100).toFixed(1) + '%', 7, true) +
               '   polynomial ' + old.toFixed(4) + '   true ' + t.toFixed(4) +
               '   core count off by ' + (((old*old)/(t*t) - 1) * 100).toFixed(1) + '%');
    }
  },

  // ---- 3.2 one uniform area ----------------------------------------------------
  testSRS: function () {
    this.head('3.2  One uniform area, against workbook sheet 1');
    this.say('  Workbook case: CV 0.5, 5 ha inlet, 90% confidence, +/-20%');
    var r = BCStats.srsSampleSize({ areaM2: 50000, plotM2: 100, mean: 100, sd: 50,
                                    confidence: 0.90, marginOfError: 0.20 });
    this.eq  ('possible plot locations N', r.N, 500);
    this.near('coefficient of variation',  r.cv, 0.50);
    this.eq  ('cores required',            r.n, 17);

    // and the same site with the prior the tool now offers
    var eg = BCStats.regionalPrior('Eelgrass', 30);
    var live = BCStats.srsSampleSize({ areaM2: 50000, plotM2: 100, mean: eg.mean, sd: eg.sd,
                                       confidence: 0.90, marginOfError: 0.20 });
    this.say('  Same site on the PNW eelgrass prior (CV ' + eg.cv.toFixed(2) + '): ' +
             live.n + ' cores');
  },

  // ---- 3.3 sensitivity grid ----------------------------------------------------
  testSensitivityGrid: function () {
    this.head('3.3  Sensitivity grid, against workbook sheet 5 (48 cells)');
    var Es  = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30];
    var CVs = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    var want = [[82,24,11,7,4,3], [129,40,19,11,7,5], [176,60,29,17,11,8],
                [220,82,40,24,16,11], [258,105,53,32,21,15], [291,129,67,40,27,19],
                [319,153,82,50,33,24], [343,176,98,60,40,29]];
    var bad = 0, row, col, line;

    this.say('     CV        5%    10%    15%    20%    25%    30%');
    for (row = 0; row < CVs.length; row++) {
      line = '    ' + this.pad(CVs[row].toFixed(2), 6);
      for (col = 0; col < Es.length; col++) {
        var got = BCStats.srsSampleSize({ areaM2: 50000, plotM2: 100, mean: 100,
                                          sd: 100 * CVs[row], confidence: 0.90,
                                          marginOfError: Es[col] }).n;
        if (got !== want[row][col]) bad++;
        line += this.pad(got, 7, true);
      }
      this.say(line);
    }
    this.eq('all 48 grid cells match', bad, 0);
  },

  // ---- 3.4 stratified ----------------------------------------------------------
  testStratified: function () {
    this.head('3.4  Stratified, against workbook sheet 2');
    var strata = [{ name: 'Dense meadow',  areaM2: 30000, mean: 20.6, sd: 11.9 },
                  { name: 'Sparse fringe', areaM2: 20000, mean: 17.1, sd:  6.9 }];

    var p = BCStats.stratifiedSampleSize({ strata: strata, plotM2: 100, confidence: 0.90,
                                           marginOfError: 0.20, allocation: 'proportional',
                                           minPerStratum: 5 });
    this.near('area-weighted overall mean', p.overallMean, 19.2, 1e-9);
    this.near('pooled variance',            p.pooledVariance, 104.01, 1e-9);
    this.near('V (variance / mean^2)',      p.V, 0.2821452);
    this.eq  ('cores before allocation',    p.n, 19);
    this.eq  ('dense meadow cores',         p.strata[0].cores, 12);
    this.eq  ('sparse fringe cores',        p.strata[1].cores, 8);
    this.eq  ('total after the minimum',    p.nAllocated, 20);

    this.say('');
    this.say('  Same site under Neyman allocation:');
    var ny = BCStats.stratifiedSampleSize({ strata: strata, plotM2: 100, confidence: 0.90,
                                            marginOfError: 0.20, allocation: 'neyman',
                                            minPerStratum: 5 });
    this.say('    proportional  ' + p.n + ' cores -> ' +
             p.strata[0].cores + ' / ' + p.strata[1].cores +
             '  = ' + p.nAllocated + ' collected');
    this.say('    neyman        ' + ny.n + ' cores -> ' +
             ny.strata[0].cores + ' / ' + ny.strata[1].cores +
             '  = ' + ny.nAllocated + ' collected');
    this.eq('neyman n never exceeds proportional', ny.n <= p.n, true);

    this.say('');
    this.say('  The bug in the earlier tools: they computed ' + ny.n +
             ' with the Neyman formula,');
    this.say('  then split it by area — a design that needs ' + p.n + '. Under-sampled.');
  },

  // ---- 3.5 achieved precision --------------------------------------------------
  testAchievedPrecision: function () {
    this.head('3.5  Achieved precision, against workbook sheet 4');
    var a = BCStats.achievedPrecisionSRS({ N: 500, n: 22, sampleMean: 21.4, sampleSd: 12.8,
                                           confidence: 0.90, target: 0.20 });
    this.near('standard error',            a.se, 2.6682565);
    this.near('achieved margin of error',  a.rme, 0.2050884);
    this.eq  ('verdict is a miss',         a.pass, false);
    this.say('  20.5% achieved against a 20% target — see Appendix A8.');
  },

  // ---- 3.6 priors --------------------------------------------------------------
  testPriors: function () {
    this.head('3.6  Priors — Pacific Northwest pooled, Janousek et al. (2025)');
    this.say('  Alaska + British Columbia + Washington + Oregon. Mangrove excluded:');
    this.say('  it does not occur north of California.');
    this.say('');
    this.say('    ecosystem        depth      mean      SD      CV    cores');

    var ecos = BCStats.ecosystems(), bad = 0, i, d, depths = [30, 100];
    for (i = 0; i < ecos.length; i++) {
      for (d = 0; d < depths.length; d++) {
        var p = BCStats.regionalPrior(ecos[i], depths[d]);
        if (p.ok) {
          this.say('    ' + this.pad(ecos[i], 16) + this.pad(depths[d] + ' cm', 10) +
                   this.pad(p.mean.toFixed(1), 9, true) + this.pad(p.sd.toFixed(1), 8, true) +
                   this.pad(p.cv.toFixed(2), 8, true) + this.pad(p.cores, 8, true) +
                   (p.indicative ? '  [indicative]' : '') +
                   (p.highVariability ? '  [high variability]' : ''));
          // a shifted column or a swapped field shows up here
          if (!(p.cv > 0.10 && p.cv < 2.0) || !(p.mean > 5) ||
              !(p.cores >= 3) || !(p.sd < p.mean * 2)) bad++;
        } else {
          this.say('    ' + this.pad(ecos[i], 16) + this.pad(depths[d] + ' cm', 10) +
                   '  not enough published cores (' + p.cores + ')');
        }
      }
    }
    this.eq('every prior has a plausible mean, SD, CV and core count', bad, 0);

    var eg = BCStats.regionalPrior('Eelgrass', 30);
    this.eq('eelgrass mean is 24.8',  eg.mean, 24.8);
    this.eq('eelgrass SD is 16.8',    eg.sd,   16.8);
    this.eq('eelgrass CV is 0.68',    eg.cv.toFixed(2), '0.68');
    this.eq('eelgrass pools 175 cores', eg.cores, 175);

    var deep = BCStats.regionalPrior('Eelgrass', 100);
    this.eq('eelgrass at 100 cm, mean is 86.5', deep.mean, 86.5);
    this.eq('eelgrass at 100 cm, SD is 41.7',   deep.sd,   41.7);
    this.eq('mean rises with depth, as it must', deep.mean > eg.mean, true);

    var flat = BCStats.regionalPrior('Tideflat', 100);
    this.eq('tideflat at 100 cm refuses rather than guessing', flat.ok, false);
    this.say('    ' + flat.reason);

    this.eq('mangrove is not offered',
            BCStats.ecosystems().join(',').indexOf('Mangrove'), -1);
    this.eq('four ecosystems offered', BCStats.ecosystems().length, 4);
  },

  // ---- 3.7 analysis scale ------------------------------------------------------
  testScale: function () {
    this.head('3.7  Analysis scale (the fixed-250 m failure)');
    var sites = [[50000, '5 ha inlet'], [500000, '50 ha'], [5e6, '500 ha'],
                 [5e7, '5,000 ha'], [5e8, '50,000 ha'], [5e9, '500,000 ha']];
    for (var i = 0; i < sites.length; i++) {
      var r = BCGeom.analysisScale(sites[i][0]);
      this.say('    ' + this.pad(sites[i][1], 14) + ' -> ' + this.pad(r.scale, 4, true) +
               ' m, ' + this.pad(r.pixels, 8, true) + ' px' +
               (r.coarse ? '   [flagged as coarse]' : ''));
    }
    this.say('');
    this.say('    At the old fixed 250 m, the 5 ha inlet resolved to ' +
             Math.floor(50000 / (250 * 250)) + ' pixels.');
    this.eq('5 ha site now uses 10 m',       BCGeom.analysisScale(50000).scale, 10);
    this.eq('500,000 ha site capped at 250 m', BCGeom.analysisScale(5e9).scale, 250);
  },

  // ---- 3.8 stratification availability -----------------------------------------
  testStratificationOptions: function () {
    this.head('3.8  Which stratification methods suit a 5 ha meadow');
    var opts = BCGeom.stratificationOptions(50000);
    for (var i = 0; i < opts.length; i++) {
      this.say('    ' + (opts[i].available ? 'offer  ' : 'hide   ') +
               this.pad(opts[i].label, 44) +
               (opts[i].reason ? opts[i].reason :
                 (opts[i].pixels !== undefined ? opts[i].pixels + ' px' : '')));
    }
    function byId(list, id) {
      for (var k = 0; k < list.length; k++) { if (list[k].id === id) return list[k]; }
    }
    this.eq('Copernicus 100 m hidden at 5 ha',   byId(opts, 'copernicus').available, false);
    this.eq('30 m covariates hidden at 5 ha',    byId(opts, 'covariates').available, false);
    this.eq('Satellite Embeddings offered',      byId(opts, 'embeddings').available, true);
    this.eq('Copernicus returns at 500 ha',
            byId(BCGeom.stratificationOptions(5e6), 'copernicus').available, true);
  },

  // ---- 3.9 edge buffer ---------------------------------------------------------
  testBuffer: function () {
    this.head('3.9  Edge buffer (the fixed-50 m failure)');
    var cases = [[50000, 22, 'whole 5 ha inlet'], [30000, 12, 'dense meadow, 3 ha'],
                 [20000,  8, 'sparse fringe, 2 ha'], [2000,  5, 'small zone, 0.2 ha']];
    for (var i = 0; i < cases.length; i++) {
      var b = BCGeom.chooseBuffer(cases[i][0], cases[i][1], 100);
      this.say('    ' + this.pad(cases[i][2], 22) + ' -> ' + b.note);
    }
    this.say('');
    this.say('    A fixed 50 m buffer would leave ' +
             Math.round(100 * BCGeom.retainedAfterBuffer(30000, 50) / 30000) +
             '% of the dense meadow.');
    this.eq('dense meadow steps down to 10 m', BCGeom.chooseBuffer(30000, 12, 100).buffer, 10);
    this.eq('5 ha inlet uses 25 m',            BCGeom.chooseBuffer(50000, 22, 100).buffer, 25);
  },

  // ---- 3.10 feasibility --------------------------------------------------------
  testFeasibility: function () {
    this.head('3.10  Feasibility checks');

    var strata = [{ name: 'Dense meadow',  areaM2: 30000, mean: 20.6, sd: 11.9 },
                  { name: 'Sparse fringe', areaM2: 20000, mean: 17.1, sd:  6.9 }];
    var alloc = BCStats.stratifiedSampleSize({ strata: strata, plotM2: 100, confidence: 0.90,
                                               marginOfError: 0.20, allocation: 'proportional',
                                               minPerStratum: 5 });
    var forGeom = [], i;
    for (i = 0; i < alloc.strata.length; i++) {
      forGeom.push({ name: alloc.strata[i].name, areaM2: alloc.strata[i].areaM2,
                     cores: alloc.strata[i].cores });
    }
    var f = BCGeom.feasibility({ plotM2: 100, minPerStratum: 5, strata: forGeom });
    this.say('    Tsawwassen design, ' + alloc.nAllocated + ' cores: ' +
             (f.ok ? 'feasible' : 'not feasible') +
             '   (N = ' + f.N + ', analysis at ' + f.scale + ' m)');
    for (i = 0; i < f.warnings.length; i++) this.say('      warning: ' + f.warnings[i]);
    this.eq('workshop design is feasible', f.ok, true);

    this.say('');
    this.say('    Designs that should be rejected:');

    var sliver = BCGeom.feasibility({ plotM2: 100, minPerStratum: 5,
                                      strata: [{ name: 'Sliver', areaM2: 300, cores: 5 }] });
    for (i = 0; i < sliver.problems.length; i++) this.say('      - ' + sliver.problems[i]);
    this.eq('0.03 ha zone rejected', sliver.ok, false);

    var packed = BCGeom.feasibility({ plotM2: 100, minPerStratum: 5,
                                      strata: [{ name: 'Packed', areaM2: 8000, cores: 70 }] });
    for (i = 0; i < packed.problems.length; i++) this.say('      - ' + packed.problems[i]);
    this.eq('70 cores in 0.8 ha rejected', packed.ok, false);

    var sub = BCGeom.feasibility({ plotM2: 100, minPerStratum: 5, areaM2: 60, cores: 5 });
    for (i = 0; i < sub.problems.length; i++) this.say('      - ' + sub.problems[i]);
    this.eq('site smaller than one plot rejected', sub.ok, false);
  },

  // ---- 3.11 layout geometry ----------------------------------------------------
  testLayout: function () {
    this.head('3.11  Layout geometry');

    // 5 ha, near-square: the shape of the worked example
    var square = [[-123.09153,49.00399],[-123.08847,49.00399],
                  [-123.08847,49.00600],[-123.09153,49.00600],[-123.09153,49.00399]];
    // 5 ha, 400 m by 125 m: the shape of a shore-fringing meadow
    var wide   = [[-123.0927,49.00494],[-123.0873,49.00494],
                  [-123.0873,49.00606],[-123.0927,49.00606],[-123.0927,49.00494]];

    var c = BCLayout.centroidOf(square);
    var back = BCLayout.toLonLat(
      BCLayout.toLocal(-123.0900, 49.0050, c).x,
      BCLayout.toLocal(-123.0900, 49.0050, c).y, c);
    this.eq('metric frame round-trips',
            BCLayout.metresBetween({ lon: -123.0900, lat: 49.0050 }, back) < 0.01, true);

    // cramped on purpose: 22 cores into a 60 m box must collide
    var rnd = BCLayout.rng(7), cand = [], i;
    for (i = 0; i < 3000; i++) {
      cand.push(BCLayout.toLonLat((rnd() - 0.5) * 60, (rnd() - 0.5) * 60, c));
    }
    var kept   = BCLayout.thin(cand, 10, 22);
    var before = BCLayout.minSeparation(cand.slice(0, 22));
    var after  = BCLayout.minSeparation(kept);
    this.say('    untouched, closest pair ' + before.toFixed(1) +
             ' m; after thinning ' + after.toFixed(1) + ' m');
    this.eq('all 22 cores placed', kept.length, 22);
    this.eq('minimum spacing respected', after >= 10, true);
    this.eq('thinning changed the result', before < 10, true);

    var g = BCLayout.lattice(square, 22, 50000, CONFIG.SEED);
    this.near('grid spacing = sqrt(area/n)', g.spacingM, Math.sqrt(50000 / 22), 0.01);
    this.eq('grid offers at least the cores asked for', g.points.length >= 22, true);

    // A large site needing few cores gives a step wider than the site. The earlier
    // version started at a random offset up to a full step in and could walk past
    // the box entirely, returning nothing and failing with "Invalid geometry".
    var big = [[-123.12, 48.98], [-123.06, 48.98], [-123.06, 49.03],
               [-123.12, 49.03], [-123.12, 48.98]];
    var sparse = BCLayout.lattice(big, 1, 11480000, CONFIG.SEED);
    this.say('    1 core over 1,148 ha: step ' + Math.round(sparse.spacingM) +
             ' m, ' + sparse.points.length + ' position(s)');
    this.eq('never returns an empty lattice', sparse.points.length >= 1, true);
    var many = 0, sd;
    for (sd = 1; sd <= 40; sd++) {
      if (BCLayout.lattice(big, 1, 11480000, sd).points.length < 1) many++;
    }
    this.eq('holds across 40 seeds', many, 0);

    var axW = BCLayout.principalAxis(wide);
    var axS = BCLayout.principalAxis(square);
    this.say('    400x125 m zone: elongation ' + axW.elongation.toFixed(2) +
             ', axis ' + axW.angleDeg.toFixed(0) + ' deg');
    this.say('    near-square zone: elongation ' + axS.elongation.toFixed(2) +
             ' — no usable long axis');
    this.eq('long axis found on an elongated zone', Math.abs(axW.angleDeg) < 5, true);
    this.eq('elongated zone marked reliable', axW.reliable, true);
    this.eq('near-square zone marked unreliable', axS.reliable, false);

    var t = BCLayout.transects(wide, 22, 4, CONFIG.SEED);
    this.say('    ' + t.lines + ' transects at ' + t.bearingDeg.toFixed(0) + ' deg, ' +
             Math.round(t.spacingAcrossM) + ' m apart');
    this.eq('transect count honoured', t.lines, 4);
    this.eq('transects follow the long axis', Math.abs(t.bearingDeg - 90) < 5, true);

    var tOver = BCLayout.transects(wide, 22, 4, CONFIG.SEED, 30);
    this.eq('supplied bearing overrides the axis', Math.abs(tOver.bearingDeg - 30) < 1, true);
    this.eq('near-square zone warns instead of guessing',
            BCLayout.transects(square, 22, 4, CONFIG.SEED).note.length > 0, true);

    var subs = BCLayout.compositeSubsamples({ lon: -123.09, lat: 49.005 }, 5, 8, CONFIG.SEED);
    var far = 0;
    for (i = 0; i < subs.length; i++) {
      far = Math.max(far, BCLayout.metresBetween({ lon: -123.09, lat: 49.005 }, subs[i]));
    }
    this.eq('5 subsamples generated', subs.length, 5);
    this.eq('subsamples stay inside the radius', far <= 8.01, true);

    this.eq('same seed reproduces the design',
      JSON.stringify(BCLayout.lattice(square, 22, 50000, 42)) ===
      JSON.stringify(BCLayout.lattice(square, 22, 50000, 42)), true);
    this.eq('a different seed does not',
      JSON.stringify(BCLayout.lattice(square, 22, 50000, 42)) !==
      JSON.stringify(BCLayout.lattice(square, 22, 50000, 43)), true);
  },

  // ---- runner ------------------------------------------------------------------
  run: function () {
    this.lines = []; this.passed = 0; this.failed = 0;

    this.say('BLUE CARBON SAMPLING TOOL — SELF TEST');
    this.say('version ' + CONFIG.VERSION +
             '   ·   reference: BlueCarbon_SampleAllocation_2026.xlsx');

    this.testInverseNormal();
    this.testSRS();
    this.testSensitivityGrid();
    this.testStratified();
    this.testAchievedPrecision();
    this.testPriors();
    this.testScale();
    this.testStratificationOptions();
    this.testBuffer();
    this.testFeasibility();
    this.testLayout();

    this.say('');
    this.say('-----------------------------------------------------------');
    this.say(this.failed === 0
      ? 'ALL ' + this.passed + ' CHECKS PASSED'
      : this.passed + ' passed, ' + this.failed + ' FAILED');
    this.say('-----------------------------------------------------------');

    LOG(this.lines.join('\n'));
    return { passed: this.passed, failed: this.failed };
  }
};


// =================================================================================
// === SECTION 4 — RUN =============================================================
// =================================================================================

// Tests run at the end of the file, once every section is defined.




// =================================================================================
// === SECTION 5 — EARTH ENGINE LAYER ==============================================
// ===
// === Everything that touches Earth Engine. Holds no statistics and no thresholds
// === of its own — it asks Sections 1 and 2 for those, then does the spatial work.
// ===
// === Earth Engine is asynchronous: results arrive in callbacks, not return
// === values. Every method here takes a callback(result, error).
// =================================================================================

var BCEarth = {

  MAX_PIXELS: 1e10,
  MAX_ERROR: 1,
  TILE_SCALE: 4,
  EMBEDDING_YEAR: 2024,
  EMBEDDING_SCALE: 10,

  // A ~5 ha rectangle off Tsawwassen, BC. A fixture for testing the geometry
  // path at the size the workshop actually targets — not a real meadow boundary.
  TEST_AOI: ee.Geometry.Polygon([[
    [-123.09153, 49.00399], [-123.08847, 49.00399],
    [-123.08847, 49.00600], [-123.09153, 49.00600],
    [-123.09153, 49.00399]
  ]]),

  // --- Area ---------------------------------------------------------------------
  // Geodesic area from the geometry itself, not by counting pixels. Exact at any
  // size, and immune to the scale problem that broke the earlier tools.

  areaOf: function (geometry, callback) {
    geometry.area({ maxError: this.MAX_ERROR }).evaluate(function (m2, err) {
      if (err || !m2) { callback(null, err || 'Could not measure that boundary.'); return; }
      callback({ areaM2: m2, areaHa: m2 / 10000, scale: BCGeom.analysisScale(m2) }, null);
    });
  },

  // --- Stratification: land cover -----------------------------------------------

  LANDCOVER: {
    dynamic: {
      collection: 'GOOGLE/DYNAMICWORLD/V1', band: 'label', scale: 10,
      labels: { 0: 'Water', 1: 'Trees', 2: 'Grass', 3: 'Flooded vegetation', 4: 'Crops',
                5: 'Shrub and scrub', 6: 'Built area', 7: 'Bare ground', 8: 'Snow and ice' }
    },
    copernicus: {
      image: 'COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019',
      band: 'discrete_classification', scale: 100,
      labels: { 20: 'Shrubland', 30: 'Herbaceous vegetation', 50: 'Urban', 60: 'Bare or sparse',
                80: 'Permanent water', 90: 'Herbaceous wetland', 200: 'Ocean' }
    }
  },

  landcoverClasses: function (aoi, sourceId, callback) {
    var cfg = this.LANDCOVER[sourceId];
    if (!cfg) { callback(null, 'Unknown land cover source: ' + sourceId); return; }

    var img = cfg.collection
      ? ee.ImageCollection(cfg.collection).filterBounds(aoi)
          .filterDate('2022-01-01', '2023-01-01').select(cfg.band).mode()
      : ee.Image(cfg.image).select(cfg.band);

    img.rename('class').clip(aoi).reduceRegion({
      reducer: ee.Reducer.frequencyHistogram(), geometry: aoi,
      scale: cfg.scale, maxPixels: this.MAX_PIXELS, tileScale: this.TILE_SCALE
    }).evaluate(function (res, err) {
      if (err || !res || !res['class']) {
        callback(null, err || 'No land cover classes found inside this boundary.');
        return;
      }
      var out = [], code;
      for (code in res['class']) {
        if (res['class'].hasOwnProperty(code)) {
          out.push({ code: parseInt(code, 10), pixels: res['class'][code],
                     name: cfg.labels[code] || ('Class ' + code) });
        }
      }
      out.sort(function (a, b) { return b.pixels - a.pixels; });
      callback({ classes: out, image: img, scale: cfg.scale }, null);
    });
  },

  landcoverStrata: function (aoi, sourceId, selected, renames) {
    var cfg = this.LANDCOVER[sourceId];
    var img = cfg.collection
      ? ee.ImageCollection(cfg.collection).filterBounds(aoi)
          .filterDate('2022-01-01', '2023-01-01').select(cfg.band).mode()
      : ee.Image(cfg.image).select(cfg.band);

    var from = [], to = [], zones = [], i;
    for (i = 0; i < selected.length; i++) {
      from.push(selected[i].code);
      to.push(i);
      zones.push({ code: i, sourceCode: selected[i].code,
                   name: (renames && renames[selected[i].code]) || selected[i].name });
    }

    var remapped = img.rename('c').remap(from, to, -999);
    return { image: remapped.updateMask(remapped.neq(-999)).rename('zone').clip(aoi),
             zones: zones };
  },

  // --- Stratification: unsupervised grouping ------------------------------------
  // Satellite Embeddings are the default for small sites: 64 bands at 10 m, so a
  // 5 ha meadow still holds ~500 pixels. Tiles are served per UTM zone and a
  // mosaic inherits a 1-degree default projection, so the projection is set
  // explicitly — without this the clusterer trains on effectively one pixel.

  embeddingImage: function (aoi) {
    var y = this.EMBEDDING_YEAR;
    return ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
      .filterDate(y + '-01-01', (y + 1) + '-01-01')
      .filterBounds(aoi)
      .mosaic()
      .setDefaultProjection('EPSG:4326', null, this.EMBEDDING_SCALE)
      .clip(aoi);
  },

  covariateStack: function (aoi) {
    var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterBounds(aoi).filterDate('2022-01-01', '2024-01-01')
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)).median()
      .select(['B2','B3','B4','B8','B11','B12'],
              ['blue','green','red','nir','swir1','swir2']);

    var ndvi = s2.normalizedDifference(['nir','red']).rename('ndvi');
    var ndwi = s2.normalizedDifference(['green','nir']).rename('ndwi');
    var dem  = ee.Image('USGS/SRTMGL1_003').rename('elevation');

    return s2.addBands([ndvi, ndwi, dem])
             .setDefaultProjection('EPSG:4326', null, 30).clip(aoi);
  },

  clusterStrata: function (aoi, sourceId, k, areaM2, callback) {
    var img   = (sourceId === 'embeddings') ? this.embeddingImage(aoi) : this.covariateStack(aoi);
    var scale = (sourceId === 'embeddings') ? this.EMBEDDING_SCALE : 30;

    // Training size follows the site, so a small meadow is not asked for more
    // pixels than it contains.
    var available = Math.floor(areaM2 / (scale * scale));
    var training  = Math.max(200, Math.min(5000, Math.floor(available * 0.5)));

    if (available < k * 20) {
      callback(null, 'This site holds about ' + available + ' pixels at ' + scale +
                     ' m — too few to separate ' + k + ' zones. Use fewer zones, or draw them.');
      return;
    }

    var sample = img.sample({ region: aoi, scale: scale, numPixels: training,
                              seed: CONFIG.SEED, geometries: false });

    sample.size().evaluate(function (n, err) {
      if (err || !n || n < k * 10) {
        callback(null, 'Only ' + (n || 0) + ' usable pixels were found. Try drawing the zones instead.');
        return;
      }
      var clusterer = ee.Clusterer.wekaKMeans({ nClusters: k, seed: CONFIG.SEED }).train(sample);
      var zones = [];
      for (var i = 0; i < k; i++) zones.push({ code: i, sourceCode: i, name: 'Zone ' + (i + 1) });
      callback({ image: img.cluster(clusterer).rename('zone').clip(aoi),
                 zones: zones, trainingPixels: n, scale: scale }, null);
    });
  },

  // --- Zone areas ---------------------------------------------------------------
  // Uses the adaptive scale from Section 2 rather than a fixed 250 m.

  zoneAreas: function (zoneImage, zones, aoi, scale, callback) {
    ee.Image.pixelArea().addBands(zoneImage.rename('zone')).reduceRegion({
      reducer: ee.Reducer.sum().group({ groupField: 1, groupName: 'zone' }),
      geometry: aoi, scale: scale, maxPixels: this.MAX_PIXELS, tileScale: this.TILE_SCALE
    }).evaluate(function (res, err) {
      if (err || !res || !res.groups || res.groups.length === 0) {
        callback(null, err || 'No zones were found inside the boundary.');
        return;
      }
      var out = [], i, j;
      for (i = 0; i < res.groups.length; i++) {
        var g = res.groups[i], name = 'Zone ' + g.zone;
        for (j = 0; j < zones.length; j++) { if (zones[j].code === g.zone) name = zones[j].name; }
        if (g.sum > 0) out.push({ code: g.zone, name: name, areaM2: g.sum });
      }
      out.sort(function (a, b) { return b.areaM2 - a.areaM2; });
      callback({ zones: out, scale: scale }, null);
    });
  },

  // --- Zones the user supplies, drawn or uploaded --------------------------------

  vectorZoneAreas: function (features, callback) {
    var fc = ee.FeatureCollection(features).map(function (f) {
      return f.set('area_m2', f.geometry().area({ maxError: 1 }));
    });
    ee.Dictionary({ names: fc.aggregate_array('zone'),
                    areas: fc.aggregate_array('area_m2') })
      .evaluate(function (r, err) {
        if (err || !r || !r.names || r.names.length === 0) {
          callback(null, err || 'Those zones have no area.'); return;
        }
        var out = [], i;
        for (i = 0; i < r.names.length; i++) {
          out.push({ code: i, name: r.names[i], areaM2: r.areas[i] });
        }
        callback({ zones: out }, null);
      });
  },

  uploadedZones: function (assetId, field, aoi, callback) {
    var fc;
    try { fc = ee.FeatureCollection(assetId.trim()); }
    catch (e) { callback(null, 'Could not open that asset.'); return; }

    var tagged = fc.filterBounds(aoi).map(function (f) {
      return f.set('zone', ee.String(f.get(field)))
              .set('area_m2', f.geometry().area({ maxError: 1 }));
    });

    ee.Dictionary({ names: tagged.aggregate_array('zone'),
                    areas: tagged.aggregate_array('area_m2') })
      .evaluate(function (r, err) {
        if (err || !r || !r.names || r.names.length === 0) {
          callback(null, 'No features with a "' + field + '" value fell inside your site.');
          return;
        }
        // several polygons may share a zone name; add their areas together
        var byName = {}, order = [], i;
        for (i = 0; i < r.names.length; i++) {
          var n = String(r.names[i]);
          if (byName[n] === undefined) { byName[n] = 0; order.push(n); }
          byName[n] += r.areas[i];
        }
        var out = [];
        for (i = 0; i < order.length; i++) {
          out.push({ code: i, name: order[i], areaM2: byName[order[i]] });
        }
        callback({ zones: out, collection: tagged }, null);
      });
  },

  // A raster zone has to become a polygon before cores can be placed inside it.
  zoneGeometry: function (zoneImage, code, aoi, scale, callback) {
    var vec = zoneImage.eq(code).selfMask().reduceToVectors({
      geometry: aoi, scale: scale, geometryType: 'polygon',
      eightConnected: false, maxPixels: this.MAX_PIXELS, tileScale: this.TILE_SCALE
    });
    vec.size().evaluate(function (n, err) {
      if (err || !n) { callback(null, err || 'That zone has no mapped area.'); return; }
      callback(vec.union(1).geometry(), null);
    });
  },

  // --- Buffer, verified against real geometry -----------------------------------
  // Section 2 estimates the buffer assuming a compact shape. A long thin fringe
  // loses far more than that estimate, so the ladder is walked again here against
  // the actual polygon.

  verifyBuffer: function (geometry, coresNeeded, plotM2, callback) {
    var self = this;
    geometry.area({ maxError: this.MAX_ERROR }).evaluate(function (original, err) {
      if (err || !original) { callback(null, err || 'Could not measure that zone.'); return; }

      var ladder = BCGeom.BUFFER_LADDER;

      function attempt(i) {
        if (i >= ladder.length) {
          callback({ buffer: 0, geometry: geometry, areaM2: original,
                     retainedFraction: 1, verified: true,
                     note: 'No edge buffer — this zone is too small to give any up.' }, null);
          return;
        }
        var b = ladder[i];
        if (b === 0) { attempt(ladder.length); return; }

        var shrunk = geometry.buffer(-b, self.MAX_ERROR);
        shrunk.area({ maxError: self.MAX_ERROR }).evaluate(function (kept, e2) {
          var frac = (kept || 0) / original;
          if (e2 || !kept ||
              frac < BCGeom.MIN_AREA_RETAINED ||
              BCGeom.practicalCapacity(kept, plotM2) < coresNeeded) {
            attempt(i + 1);
            return;
          }
          callback({ buffer: b, geometry: shrunk, areaM2: kept,
                     retainedFraction: frac, verified: true,
                     note: b + ' m edge buffer, keeping ' + Math.round(frac * 100) +
                           '% of the zone.' }, null);
        });
      }
      attempt(0);
    });
  }
};


// =================================================================================
// === SECTION 5T — EARTH ENGINE TESTS =============================================
// ===
// === Needs a live connection, so these run separately from Section 3 and print
// === as they return. Set RUN_EE_TEST to true to run them.
// =================================================================================

var RUN_EE_TEST = false;

var BCEarthTest = {

  run: function (aoi) {
    aoi = aoi || BCEarth.TEST_AOI;
    LOG('EARTH ENGINE LAYER — TEST  (results arrive one at a time)');

    BCEarth.areaOf(aoi, function (r, err) {
      if (err) { LOG('  area: FAILED — ' + err); return; }

      LOG('  area          ' + Math.round(r.areaM2) + ' m2  (' + r.areaHa.toFixed(2) + ' ha)');
      LOG('  scale chosen  ' + r.scale.scale + ' m, about ' + r.scale.pixels + ' pixels');
      LOG('  expected      about 50000 m2 for the test fixture');

      var opts = BCGeom.stratificationOptions(r.areaM2), i;
      LOG('  stratification methods offered at this size:');
      for (i = 0; i < opts.length; i++) {
        LOG('    ' + (opts[i].available ? 'offer  ' : 'hide   ') + opts[i].label);
      }

      BCEarth.clusterStrata(aoi, 'embeddings', 2, r.areaM2, function (c, cErr) {
        if (cErr) { LOG('  embeddings: FAILED — ' + cErr); return; }
        LOG('  embeddings    trained on ' + c.trainingPixels + ' pixels at ' + c.scale + ' m');

        BCEarth.zoneAreas(c.image, c.zones, aoi, r.scale.scale, function (z, zErr) {
          if (zErr) { LOG('  zone areas: FAILED — ' + zErr); return; }

          var total = 0, k;
          for (k = 0; k < z.zones.length; k++) {
            LOG('    ' + z.zones[k].name + ': ' + Math.round(z.zones[k].areaM2) + ' m2');
            total += z.zones[k].areaM2;
          }
          LOG('  zone areas sum to ' + Math.round(total) + ' m2, against a boundary of ' +
              Math.round(r.areaM2) + ' m2  (' +
              (100 * total / r.areaM2).toFixed(1) + '% — should be near 100)');
        });
      });

      BCEarth.verifyBuffer(aoi, 22, CONFIG.PLOT_M2, function (b, bErr) {
        if (bErr) { LOG('  buffer: FAILED — ' + bErr); return; }
        LOG('  buffer        ' + b.note);
        LOG('  estimated     ' + BCGeom.chooseBuffer(50000, 22, 100).note);
      });
    });
  }
};

if (RUN_EE_TEST) { BCEarthTest.run(); }


// =================================================================================
// === SECTION 6 — LAYOUT ENGINE ===================================================
// ===
// === Where cores go inside a zone. Stratification (Section 5) and layout (here)
// === are independent: any zone definition combines with any layout, which is why
// === this replaces three separate tools.
// ===
// === 6A is pure geometry and is covered by the tests in Section 3.
// === 6B wraps it in Earth Engine calls.
// =================================================================================

// --- 6A · GEOMETRY, NO EARTH ENGINE ---------------------------------------------

var BCLayout = {

  LAYOUTS: [
    { id: 'random',    label: 'Random',
      blurb: 'Cores scattered at random. The default when a zone looks uniform.' },
    { id: 'grid',      label: 'Even grid',
      blurb: 'Cores on a regular lattice. Guarantees even coverage.' },
    { id: 'transect',  label: 'Shore-parallel transects',
      blurb: 'Lines following the long axis of the zone. Recommended for eelgrass.' },
    { id: 'composite', label: 'Composite plots',
      blurb: 'Clusters of subsamples combined into one sample. Fewer lab analyses.' }
  ],

  // --- local metric frame -------------------------------------------------------
  // Coordinates near a small site are converted to metres about a centroid, the
  // geometry is done flat, then converted back. Error over a few kilometres is
  // far below the precision anything here needs.

  M_PER_DEG_LAT: 111320,

  toLocal: function (lon, lat, c) {
    return { x: (lon - c.lon) * this.M_PER_DEG_LAT * Math.cos(c.lat * Math.PI / 180),
             y: (lat - c.lat) * this.M_PER_DEG_LAT };
  },

  toLonLat: function (x, y, c) {
    return { lon: c.lon + x / (this.M_PER_DEG_LAT * Math.cos(c.lat * Math.PI / 180)),
             lat: c.lat + y / this.M_PER_DEG_LAT };
  },

  centroidOf: function (ring) {
    var sx = 0, sy = 0;
    for (var i = 0; i < ring.length; i++) { sx += ring[i][0]; sy += ring[i][1]; }
    return { lon: sx / ring.length, lat: sy / ring.length };
  },

  metresBetween: function (a, b) {
    var R = 6371000;
    var p1 = a.lat * Math.PI / 180, p2 = b.lat * Math.PI / 180;
    var dp = (b.lat - a.lat) * Math.PI / 180, dl = (b.lon - a.lon) * Math.PI / 180;
    var h = Math.sin(dp/2)*Math.sin(dp/2) +
            Math.cos(p1)*Math.cos(p2)*Math.sin(dl/2)*Math.sin(dl/2);
    return 2 * R * Math.asin(Math.min(1, Math.sqrt(h)));
  },

  // --- minimum separation -------------------------------------------------------
  // Earth Engine's randomPoints has no spacing control, so two cores can land 3 m
  // apart while each claims a 100 m2 plot — one observation counted twice. Points
  // are oversampled, then thinned greedily in the order returned, which keeps the
  // result reproducible for a given seed.

  thin: function (candidates, spacingM, wanted) {
    var kept = [], i, j, ok;
    for (i = 0; i < candidates.length && kept.length < wanted; i++) {
      ok = true;
      for (j = 0; j < kept.length; j++) {
        if (this.metresBetween(candidates[i], kept[j]) < spacingM) { ok = false; break; }
      }
      if (ok) kept.push(candidates[i]);
    }
    return kept;
  },

  minSeparation: function (points) {
    var lowest = Infinity, i, j;
    for (i = 0; i < points.length; i++) {
      for (j = i + 1; j < points.length; j++) {
        var d = this.metresBetween(points[i], points[j]);
        if (d < lowest) lowest = d;
      }
    }
    return points.length < 2 ? Infinity : lowest;
  },

  // --- deterministic jitter -----------------------------------------------------
  // A tiny generator so grid offsets and transect jitter reproduce exactly from
  // the seed, rather than depending on when a callback happened to return.

  rng: function (seed) {
    var s = seed || 1;
    return function () { s = (s * 1103515245 + 12345) % 2147483648; return s / 2147483648; };
  },

  // --- even grid ----------------------------------------------------------------

  lattice: function (ring, wanted, areaM2, seed) {
    var c = this.centroidOf(ring), i;
    var spacing = Math.sqrt(areaM2 / wanted);

    var xs = [], ys = [];
    for (i = 0; i < ring.length; i++) {
      var p = this.toLocal(ring[i][0], ring[i][1], c);
      xs.push(p.x); ys.push(p.y);
    }
    var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs);
    var minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);

    // Centre the lattice in the bounding box and jitter by less than half a step.
    // The previous version started at a random offset up to a full step in, which
    // on a large site with few cores could step straight past the box and return
    // nothing at all.
    var width = maxX - minX, height = maxY - minY;
    var cols = Math.max(1, Math.floor(width / spacing) + 1);
    var rows = Math.max(1, Math.floor(height / spacing) + 1);

    var rand = this.rng(seed);
    var startX = minX + (width  - (cols - 1) * spacing) / 2 + (rand() - 0.5) * spacing * 0.4;
    var startY = minY + (height - (rows - 1) * spacing) / 2 + (rand() - 0.5) * spacing * 0.4;

    var out = [], i, j;
    for (i = 0; i < cols; i++) {
      for (j = 0; j < rows; j++) {
        out.push(this.toLonLat(startX + i * spacing, startY + j * spacing, c));
      }
    }
    return { points: out, spacingM: spacing, cols: cols, rows: rows };
  },

  // --- principal axis -----------------------------------------------------------
  // The long axis of an eelgrass meadow generally follows the shore, so transects
  // laid along it are shore-parallel in practice. This is a proxy for a real
  // coastline, and the tool says so rather than implying otherwise.

  // A near-square zone has no meaningful long axis — the answer is arbitrary and
  // transects laid along it would be too. Elongation is reported so the caller can
  // say so and ask for a bearing instead of quietly picking one.
  MIN_ELONGATION: 1.20,

  principalAxis: function (ring) {
    var c = this.centroidOf(ring), pts = [], i;

    // drop the closing vertex, which would otherwise weight one corner twice
    var n = ring.length;
    if (n > 2 && ring[0][0] === ring[n-1][0] && ring[0][1] === ring[n-1][1]) n -= 1;
    for (i = 0; i < n; i++) pts.push(this.toLocal(ring[i][0], ring[i][1], c));

    var mx = 0, my = 0;
    for (i = 0; i < pts.length; i++) { mx += pts[i].x; my += pts[i].y; }
    mx /= pts.length; my /= pts.length;

    var sxx = 0, syy = 0, sxy = 0;
    for (i = 0; i < pts.length; i++) {
      var dx = pts[i].x - mx, dy = pts[i].y - my;
      sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }

    var angle = 0.5 * Math.atan2(2 * sxy, sxx - syy);   // radians, 0 = east-west

    // eigenvalues of the 2x2 covariance give the long:short ratio
    var mid  = (sxx + syy) / 2;
    var half = Math.sqrt(Math.pow((sxx - syy) / 2, 2) + sxy * sxy);
    var big  = mid + half, small = Math.max(mid - half, 1e-9);
    var elongation = Math.sqrt(big / small);

    return { angleRad: angle,
             angleDeg: angle * 180 / Math.PI,
             elongation: elongation,
             reliable: elongation >= this.MIN_ELONGATION,
             centre: c,
             centreLocal: { x: mx, y: my } };
  },

  // --- shore-parallel transects -------------------------------------------------

  // bearingDegOverride: set it when the zone is too round for an axis to mean
  // anything, or when the operator knows the shoreline runs another way.
  transects: function (ring, cores, transectCount, seed, bearingDegOverride) {
    var axis = this.principalAxis(ring);
    var c = axis.centre;

    var useRad = (bearingDegOverride === undefined || bearingDegOverride === null)
                   ? axis.angleRad
                   : (90 - bearingDegOverride) * Math.PI / 180;

    var ca = Math.cos(useRad), sa = Math.sin(useRad);
    var i;

    // rotate the ring into axis-aligned coordinates
    var u = [], v = [];
    for (i = 0; i < ring.length; i++) {
      var p = this.toLocal(ring[i][0], ring[i][1], c);
      u.push(p.x * ca + p.y * sa);          // along the axis
      v.push(-p.x * sa + p.y * ca);         // across it
    }
    var minU = Math.min.apply(null, u), maxU = Math.max.apply(null, u);
    var minV = Math.min.apply(null, v), maxV = Math.max.apply(null, v);

    var lines = Math.max(1, transectCount || Math.max(2, Math.round(Math.sqrt(cores / 2))));
    var perLine = Math.ceil(cores / lines);
    var rand = this.rng(seed);

    var out = [], line, k;
    for (line = 0; line < lines; line++) {
      var vPos = minV + (maxV - minV) * (line + 0.5) / lines;
      for (k = 0; k < perLine; k++) {
        var uPos = minU + (maxU - minU) * (k + 0.5) / perLine;
        uPos += (rand() - 0.5) * (maxU - minU) / (perLine * 4);   // small jitter
        var x = uPos * ca - vPos * sa;
        var y = uPos * sa + vPos * ca;
        out.push(this.toLonLat(x, y, c));
      }
    }
    return { points: out, lines: lines, perLine: perLine,
             bearingDeg: (90 - (useRad * 180 / Math.PI) + 360) % 180,
             elongation: axis.elongation,
             axisReliable: axis.reliable,
             bearingSource: (bearingDegOverride === undefined || bearingDegOverride === null)
                              ? 'long axis of the zone' : 'bearing you supplied',
             note: axis.reliable
                     ? ''
                     : 'This zone is nearly as wide as it is long, so it has no clear ' +
                       'long axis. Set the shoreline bearing yourself, or use the even grid.',
             spacingAcrossM: (maxV - minV) / lines,
             spacingAlongM: (maxU - minU) / perLine };
  },

  // --- composite plots ----------------------------------------------------------
  // One composite = several subsamples pooled into a single lab analysis. Cuts
  // analysis cost, at the price of losing within-plot variability.

  compositeSubsamples: function (centre, subsampleCount, radiusM, seed) {
    var rand = this.rng(seed), out = [], i;
    for (i = 0; i < subsampleCount; i++) {
      var ang = 2 * Math.PI * (i + rand() * 0.5) / subsampleCount;
      var rad = radiusM * Math.sqrt(0.25 + 0.75 * rand());
      var p = this.toLocal(centre.lon, centre.lat, centre);
      out.push(this.toLonLat(p.x + rad * Math.cos(ang), p.y + rad * Math.sin(ang), centre));
    }
    return out;
  }
};


// --- 6B · EARTH ENGINE WRAPPERS -------------------------------------------------

var BCPlace = {

  OVERSAMPLE: 10,

  // Candidate points are generated server-side, brought back once, then thinned
  // in Section 6A. Core counts are in the tens, so this is a single small round
  // trip rather than a per-point conversation.

  place: function (o, callback) {
    var self    = this;
    var spacing = BCGeom.minSpacing(o.plotM2);

    o.geometry.bounds(BCEarth.MAX_ERROR).coordinates().evaluate(function (ring, err) {
      if (err || !ring || !ring[0]) { callback(null, err || 'Could not read that zone.'); return; }
      var box = ring[0];

      if (o.layout === 'grid' || o.layout === 'transect') {
        var built = (o.layout === 'grid')
          ? BCLayout.lattice(box, o.cores, o.areaM2, CONFIG.SEED + (o.zoneIndex || 0))
          : BCLayout.transects(box, o.cores, o.transectCount,
                               CONFIG.SEED + (o.zoneIndex || 0), o.bearingDeg);

        if (!built.points || built.points.length === 0) {
          callback(null, 'That layout produced no positions for this zone. ' +
                         'Try the random layout.');
          return;
        }

        // keep only what falls inside the real zone, not the bounding box
        var fc = ee.FeatureCollection(built.points.map(function (p) {
          return ee.Feature(ee.Geometry.Point([p.lon, p.lat]));
        })).filterBounds(o.geometry);

        fc.geometry().coordinates().evaluate(function (inside, e2) {
          if (e2 || !inside || inside.length === 0) {
            callback(null, 'No positions fell inside this zone. Try the random layout.');
            return;
          }
          var cand = inside.map(function (p) { return { lon: p[0], lat: p[1] }; });
          var kept = BCLayout.thin(cand, spacing, o.cores);
          callback({ points: kept, layout: o.layout, requested: o.cores,
                     placed: kept.length, meta: built,
                     minSeparationM: BCLayout.minSeparation(kept) }, null);
        });
        return;
      }

      // random, and the centres for composite
      var wanted = (o.layout === 'composite') ? o.compositeCount : o.cores;
      var over   = Math.min(wanted * self.OVERSAMPLE, 3000);

      ee.FeatureCollection.randomPoints({
        region: o.geometry, points: over,
        seed: CONFIG.SEED + (o.zoneIndex || 0), maxError: BCEarth.MAX_ERROR
      }).geometry().coordinates().evaluate(function (raw, e3) {
        if (e3 || !raw || raw.length === 0) {
          callback(null, e3 || 'Could not place points inside this zone.');
          return;
        }
        var cand = raw.map(function (p) { return { lon: p[0], lat: p[1] }; });
        var kept = BCLayout.thin(cand, spacing, wanted);

        if (kept.length < wanted) {
          callback({ points: kept, layout: o.layout, requested: wanted, placed: kept.length,
                     minSeparationM: BCLayout.minSeparation(kept),
                     shortfall: 'Only ' + kept.length + ' of ' + wanted +
                                ' cores fit while staying ' + spacing +
                                ' m apart. The zone is close to full.' }, null);
          return;
        }

        if (o.layout === 'composite') {
          var all = [], i, j;
          for (i = 0; i < kept.length; i++) {
            var subs = BCLayout.compositeSubsamples(kept[i], o.subsamples,
                                                    o.compositeRadiusM,
                                                    CONFIG.SEED + i);
            for (j = 0; j < subs.length; j++) {
              all.push({ lon: subs[j].lon, lat: subs[j].lat,
                         composite: i + 1, subsample: j + 1 });
            }
          }
          callback({ points: kept, subsamples: all, layout: 'composite',
                     requested: wanted, placed: kept.length,
                     minSeparationM: BCLayout.minSeparation(kept) }, null);
          return;
        }

        callback({ points: kept, layout: 'random', requested: wanted, placed: kept.length,
                   minSeparationM: BCLayout.minSeparation(kept) }, null);
      });
    });
  },

  // Turns placed points into an exportable collection carrying enough metadata
  // that a reviewer can reconstruct the design from the file alone.

  toFeatures: function (points, zoneName, layout, plotM2) {
    // An empty list would build an empty collection, and asking Earth Engine for
    // its geometry raises "Invalid geometry" rather than anything informative.
    if (!points || points.length === 0) return null;
    return ee.FeatureCollection(points.map(function (p, i) {
      return ee.Feature(ee.Geometry.Point([p.lon, p.lat]), {
        core_id: 'BC_' + (i + 1 < 10 ? '0' : '') + (i + 1),
        zone: zoneName, layout: layout, plot_m2: plotM2,
        lon: p.lon, lat: p.lat,
        composite: p.composite || null, subsample: p.subsample || null
      });
    }));
  }
};



// =================================================================================
// === SECTION 7 — USER INTERFACE ==================================================
// ===
// === Six steps matching Part 2 of the workshop. Plain language on the surface;
// === anything statistical sits behind a "Show the working" toggle.
// ===
// === The interface holds no arithmetic. It collects answers, calls Sections 1,
// === 2, 5 and 6, and displays what comes back.
// =================================================================================

var UI = {
  TITLE:   { fontSize: '22px', fontWeight: 'bold', color: '#0F4C5C', margin: '8px 8px 0 8px' },
  SUB:     { fontSize: '13px', color: '#5B6B70', margin: '0 8px 10px 8px' },
  STEP:    { fontSize: '15px', fontWeight: 'bold', color: '#0F4C5C', margin: '16px 8px 2px 8px' },
  ASK:     { fontSize: '13px', color: '#26343A', margin: '2px 8px 6px 8px' },
  HINT:    { fontSize: '11px', color: '#7A8B90', margin: '0 8px 6px 8px' },
  RESULT:  { fontSize: '20px', fontWeight: 'bold', color: '#0F4C5C', margin: '6px 8px' },
  NOTE:    { fontSize: '12px', color: '#26343A', margin: '2px 8px' },
  WARN:    { fontSize: '12px', color: '#B26A00', margin: '4px 8px' },
  ERROR:   { fontSize: '12px', color: '#B3261E', margin: '4px 8px' },
  OK:      { fontSize: '12px', color: '#1A7A8C', margin: '4px 8px' },
  LINK:    { fontSize: '12px', color: '#1A7A8C', margin: '2px 8px' },
  PANEL:   { width: '430px', border: '1px solid #D8E2E5' },
  BTN:     { stretch: 'horizontal', margin: '8px' },
  WIDE:    { stretch: 'horizontal', margin: '0 8px 4px 8px' },
  ZONE_COLOURS: ['#1A7A8C', '#C8763C', '#5B8C5A', '#8C5B8C', '#4C6E8C', '#A8843C']
};

var State = {
  aoi: null, areaM2: null, scale: null,
  ecosystem: CONFIG.DEFAULT_ECOSYSTEM, depth: CONFIG.DEFAULT_DEPTH_CM, prior: null,
  zoneMode: 'none', zoneImage: null, zoneAreas: null, zoneRevision: 0,
  zoneKeep: {}, zoneGeoms: {}, drawnZones: [], draft: null, classWidgets: [],
  design: null, cores: null, layout: 'random', transectCount: 4, bearing: null,
  compositeCount: 10, subsamples: 5, compositeRadius: 8,
  placed: null, features: null,

  reset: function () {
    this.aoi = null; this.areaM2 = null; this.scale = null;
    this.prior = null; this.zoneMode = 'none'; this.design = null;
    this.drawnZones = []; this.draft = null; this.classWidgets = [];
    this.placed = null; this.features = null;
  }
};

ui.root.clear();
var map   = ui.Map();
var panel = ui.Panel({ style: UI.PANEL });
ui.root.add(ui.SplitPanel(panel, map, 'horizontal', false));
map.setCenter(-123.09, 49.005, 11);
map.setOptions('SATELLITE');

function label(t, s)  { return ui.Label(t, s); }
function clearPanel(p) { p.clear(); }

panel.add(label('Blue Carbon Sampling Design', UI.TITLE));
panel.add(label('Turn a boundary and a precision target into a list of core locations. ' +
                'Companion to Part 2 of the eelgrass workshop.', UI.SUB));

// --- STEP 1 · boundary ----------------------------------------------------------

panel.add(label('Step 1 · Where are you working?', UI.STEP));
panel.add(label('Draw your site on the map, or point at a boundary you have already saved.', UI.ASK));

var sourceSelect = ui.Select({
  items: ['Draw it on the map', 'Use a saved boundary'], value: 'Draw it on the map',
  style: UI.WIDE,
  onChange: function (v) { assetRow.style().set('shown', v === 'Use a saved boundary'); }
});
var assetBox = ui.Textbox({ placeholder: 'users/you/your_site', style: UI.WIDE });
var assetRow = ui.Panel([assetBox], null, { shown: false });
var areaOut  = ui.Panel();

panel.add(sourceSelect); panel.add(assetRow);
panel.add(ui.Button({ label: 'Measure this site', style: UI.BTN, onClick: measureSite }));
panel.add(areaOut);

// --- STEP 2 · zones -------------------------------------------------------------

panel.add(label('Step 2 · Is the site all one thing?', UI.STEP));
panel.add(label('If part of the site is denser, deeper, or under different management, ' +
                'split it into zones now. Splitting by something that actually drives ' +
                'carbon gives a tighter answer for the same number of cores.', UI.ASK));

var zoneSelect = ui.Select({ items: ['Measure the site first'], style: UI.WIDE,
                             onChange: chooseZoneMode });

// how many groups, for the automatic methods
var zoneCount = ui.Slider({ min: 2, max: 6, value: 3, step: 1, style: UI.WIDE });
var zoneCountRow = ui.Panel([label('How many groups should it look for?', UI.HINT), zoneCount],
                            null, { shown: false });

// drawing zones by hand
var zoneNameBox = ui.Textbox({ placeholder: 'Name this zone, e.g. Dense meadow', style: UI.WIDE });
var drawnList   = ui.Panel();
var drawRow = ui.Panel([
  label('Name the zone, draw it on the map, then add it. Repeat for each zone.', UI.HINT),
  zoneNameBox,
  ui.Button({ label: 'Draw this zone on the map', style: UI.BTN, onClick: startZoneDrawing }),
  ui.Button({ label: 'Add the zone I just drew',  style: UI.BTN, onClick: addDrawnZone }),
  drawnList
], null, { shown: false });

// zones from an uploaded asset
var zoneAssetBox = ui.Textbox({ placeholder: 'users/you/your_zones', style: UI.WIDE });
var zoneFieldBox = ui.Textbox({ placeholder: 'Which column holds the zone name?', style: UI.WIDE });
var uploadRow = ui.Panel([
  label('An asset with one polygon per zone, and a column naming each.', UI.HINT),
  zoneAssetBox, zoneFieldBox
], null, { shown: false });

// land cover classes, chosen after they are read from the map
var classRow = ui.Panel(null, null, { shown: false });

var zoneOut     = ui.Panel();
var zoneChoices = ui.Panel();

panel.add(zoneSelect);
panel.add(zoneCountRow); panel.add(drawRow); panel.add(uploadRow); panel.add(classRow);
panel.add(ui.Button({ label: 'Build the zones', style: UI.BTN, onClick: buildZones }));
panel.add(zoneOut);
panel.add(zoneChoices);

// --- STEP 3 · what is being measured, and how variable it is --------------------

panel.add(label('Step 3 · What are you measuring?', UI.STEP));
panel.add(label('Sediment carbon, in the site and zones you just defined. How variable ' +
                'that carbon is drives the number of cores more than anything else you ' +
                'choose here.', UI.ASK));

var ecoSelect = ui.Select({
  items: BCStats.ecosystems(), value: CONFIG.DEFAULT_ECOSYSTEM, style: UI.WIDE,
  onChange: function (v) { State.ecosystem = v; refreshPriors(); }
});
var depthSelect = ui.Select({
  items: [{ label: 'Top 30 cm', value: 30 }, { label: 'Full metre', value: 100 }],
  value: 30, style: UI.WIDE,
  onChange: function (v) { State.depth = v; refreshPriors(); }
});
var priorSelect = ui.Select({
  items: [{ label: 'Use the published regional average', value: 'regional' },
          { label: 'I have my own data from this site',  value: 'own' }],
  value: 'regional', style: UI.WIDE, onChange: choosePrior
});
var priorOut = ui.Panel();

panel.add(ecoSelect);
panel.add(label('How deep are you coring?', UI.HINT));
panel.add(depthSelect);
panel.add(label('How much carbon should we expect, and how patchy?', UI.HINT));
panel.add(priorSelect); panel.add(priorOut);

// --- STEP 4 · precision ---------------------------------------------------------

panel.add(label('Step 4 · How precise do you need to be?', UI.STEP));
panel.add(label('Decide this before you see the answer, not after.', UI.ASK));

var moeSelect = ui.Select({
  items: [{ label: 'Within 10% — demanding', value: 0.10 },
          { label: 'Within 20% — usual for coastal MMRV', value: 0.20 },
          { label: 'Within 30% — a first look', value: 0.30 }],
  value: CONFIG.MARGIN_OF_ERROR, style: UI.WIDE, onChange: recompute
});
var confSelect = ui.Select({
  items: [{ label: '80% sure', value: 0.80 },
          { label: '90% sure — usual', value: 0.90 },
          { label: '95% sure', value: 0.95 }],
  value: CONFIG.CONFIDENCE, style: UI.WIDE, onChange: recompute
});
var allocSelect = ui.Select({
  items: [{ label: 'Split cores by zone area', value: 'proportional' },
          { label: 'Send more cores to patchier zones', value: 'neyman' }],
  value: 'proportional', style: UI.WIDE, onChange: recompute
});
var designOut = ui.Panel();
var mathPanel = ui.Panel(null, null, { shown: false });
var mathToggle = ui.Checkbox({
  label: 'Show the working', value: false, style: { margin: '4px 8px' },
  onChange: function (v) { mathPanel.style().set('shown', v); }
});

panel.add(moeSelect); panel.add(confSelect);
panel.add(label('When the site has zones:', UI.HINT)); panel.add(allocSelect);
panel.add(designOut); panel.add(mathToggle); panel.add(mathPanel);

// --- STEP 5 · layout ------------------------------------------------------------

panel.add(label('Step 5 · Where do the cores go?', UI.STEP));

var layoutSelect = ui.Select({
  items: BCLayout.LAYOUTS.map(function (l) { return { label: l.label, value: l.id }; }),
  value: 'random', style: UI.WIDE,
  onChange: function (v) {
    State.layout = v;
    transectRow.style().set('shown', v === 'transect');
    compositeRow.style().set('shown', v === 'composite');
    var blurb = '';
    for (var i = 0; i < BCLayout.LAYOUTS.length; i++) {
      if (BCLayout.LAYOUTS[i].id === v) blurb = BCLayout.LAYOUTS[i].blurb;
    }
    layoutHint.setValue(blurb);
  }
});
var layoutHint  = label(BCLayout.LAYOUTS[0].blurb, UI.HINT);
var transectBox = ui.Slider({ min: 2, max: 10, value: 4, step: 1, style: UI.WIDE });
var bearingBox  = ui.Textbox({ placeholder: 'Shoreline bearing in degrees (optional)', style: UI.WIDE });
var transectRow = ui.Panel([label('How many transect lines?', UI.HINT), transectBox, bearingBox],
                           null, { shown: false });
var subsampleBox = ui.Slider({ min: 3, max: 10, value: 5, step: 1, style: UI.WIDE });
var compositeRow = ui.Panel([label('Subsamples per composite:', UI.HINT), subsampleBox],
                            null, { shown: false });
var placeOut = ui.Panel();

panel.add(layoutSelect); panel.add(layoutHint);
panel.add(transectRow); panel.add(compositeRow);
panel.add(ui.Button({ label: 'Place the cores', style: UI.BTN, onClick: placeCores }));
panel.add(placeOut);

// --- STEP 6 · export ------------------------------------------------------------

panel.add(label('Step 6 · Take it to the field', UI.STEP));
var formatSelect = ui.Select({ items: ['CSV', 'GeoJSON', 'KML', 'SHP'], value: 'CSV', style: UI.WIDE });
var exportOut = ui.Panel();
var boundaryExportOut = ui.Panel();
var boundaryExportRequest = 0;
var boundaryExportSelect = ui.Select({
  items: [
    { label: 'Entire AOI', value: 'aoi' },
    { label: 'AOI with all zones', value: 'aoi_with_zones' },
    { label: 'Selected zones only', value: 'selected_zones' }
  ],
  value: 'aoi', style: UI.WIDE, onChange: clearBoundaryExport
});
var methodsOut = ui.Panel();

panel.add(formatSelect);
panel.add(ui.Button({ label: 'Download core locations', style: UI.BTN, onClick: exportCores }));
panel.add(exportOut);
panel.add(boundaryExportSelect);
panel.add(label('Selected zones use the checkboxes in Step 2.', UI.HINT));
panel.add(ui.Button({ label: 'Download boundary GeoJSON', style: UI.BTN, onClick: exportBoundary }));
panel.add(boundaryExportOut);
panel.add(label('Methods paragraph — paste this into your report:', UI.HINT));
panel.add(methodsOut);

// --- resources ------------------------------------------------------------------

panel.add(label('Workshop materials', UI.STEP));
[['Planning guide (Part 2)', 'planningGuide'],
 ['Sample allocation calculator', 'calculator'],
 ['The maths behind this (Appendix A)', 'appendixA'],
 ['Field guide', 'fieldGuide'],
 ['Field datasheets', 'datasheets'],
 ['Worked example', 'workedExample']].forEach(function (r) {
  var l = ui.Label(r[0], UI.LINK);
  l.setUrl(CONFIG.linkTo(r[1]));
  panel.add(l);
});

panel.add(ui.Button({
  label: 'Start over', style: { stretch: 'horizontal', margin: '16px 8px' },
  onClick: function () {
    State.reset(); map.layers().reset(); map.drawingTools().clear();
    clearBuiltZones(); clearPanel(drawnList);
    boundaryExportSelect.setValue('aoi', false);
    [areaOut, priorOut, zoneOut, designOut, mathPanel, placeOut, exportOut, methodsOut]
      .forEach(clearPanel);
    zoneSelect.items().reset(['Measure the site first']);
  }
}));

var tools = map.drawingTools();
tools.setShown(true); tools.setDrawModes(['polygon', 'rectangle']); tools.setShape('polygon');


// =================================================================================
// === SECTION 7B — BEHAVIOUR ======================================================
// =================================================================================

function getBoundary() {
  if (sourceSelect.getValue() === 'Draw it on the map') {
    var layers = map.drawingTools().layers();
    if (layers.length() === 0 || layers.get(0).geometries().length() === 0) return null;
    return layers.get(0).toGeometry();
  }
  var id = assetBox.getValue();
  if (!id) return null;
  return ee.FeatureCollection(id.trim()).geometry();
}

function measureSite() {
  clearPanel(areaOut);
  clearBuiltZones(); clearPanel(zoneOut);
  State.areaM2 = null; State.scale = null; State.design = null;
  var aoi = getBoundary();
  State.aoi = aoi;
  if (!aoi) {
    areaOut.add(label('Draw a polygon on the map first, or enter a saved boundary.', UI.ERROR));
    return;
  }
  areaOut.add(label('Measuring…', UI.NOTE));
  map.layers().reset();
  map.centerObject(State.aoi, 15);
  map.addLayer(State.aoi, { color: 'FFFFFF' }, 'Site boundary');

  BCEarth.areaOf(aoi, function (r, err) {
    if (State.aoi !== aoi) return;
    clearPanel(areaOut);
    if (err) { areaOut.add(label('Could not measure that: ' + err, UI.ERROR)); return; }

    State.areaM2 = r.areaM2;
    State.scale  = r.scale.scale;
    areaOut.add(label(r.areaHa.toFixed(2) + ' hectares  ·  ' +
                      BCStats.populationSize(r.areaM2, CONFIG.PLOT_M2) +
                      ' possible core positions', UI.OK));

    var opts = BCGeom.stratificationOptions(r.areaM2),
        items = [{ label: 'No — treat it as one area', value: 'none' }], i;
    for (i = 0; i < opts.length; i++) {
      if (opts[i].available) items.push({ label: opts[i].label, value: opts[i].id });
    }
    zoneSelect.items().reset(items);
    zoneSelect.setValue('none');

    var hidden = [];
    for (i = 0; i < opts.length; i++) { if (!opts[i].available) hidden.push(opts[i].label); }
    if (hidden.length > 0) {
      areaOut.add(label('Too coarse for a site this size, so not offered: ' +
                        hidden.join('; ') + '.', UI.HINT));
    }
    refreshPriors();
  });
}

var ownMean = ui.Textbox({ placeholder: 'Mean carbon, Mg C per hectare', style: UI.WIDE });
var ownSd   = ui.Textbox({ placeholder: 'Standard deviation between cores', style: UI.WIDE });

function refreshPriors() { choosePrior(priorSelect.getValue() || 'regional'); }

function choosePrior(mode) {
  clearPanel(priorOut);

  if (mode === 'own') {
    priorOut.add(label('From a pilot, an earlier survey, or nearby cores. This beats any ' +
                       'regional figure — local variability is what actually sets the ' +
                       'core count.', UI.HINT));
    priorOut.add(ownMean); priorOut.add(ownSd);
    priorOut.add(ui.Button({
      label: 'Use these numbers', style: UI.BTN,
      onClick: function () {
        var m = parseFloat(ownMean.getValue()), sd = parseFloat(ownSd.getValue());
        if (!(m > 0) || !(sd >= 0)) {
          priorOut.add(label('Enter a mean above zero and a standard deviation.', UI.ERROR));
          return;
        }
        State.prior = { mean: m, sd: sd, cv: sd / m, cores: null,
                        source: 'your own data from this site' };
        priorOut.add(label('Using your figures: carbon varies by ' +
                           Math.round(100 * sd / m) + '% between cores.', UI.OK));
        if (sd / m < 0.15) {
          priorOut.add(label('That is unusually even for coastal sediment. Check the ' +
                             'standard deviation is between cores, not a standard error.', UI.WARN));
        }
        recompute();
      }
    }));
    return;
  }

  var p = BCStats.regionalPrior(State.ecosystem, State.depth);
  if (!p.ok) {
    State.prior = null;
    priorOut.add(label(p.reason, UI.ERROR));
    return;
  }

  State.prior = p;
  priorOut.add(label('Expect about ' + p.mean.toFixed(0) + ' Mg C per hectare, varying by ' +
                     Math.round(p.cv * 100) + '% between cores.', UI.NOTE));
  priorOut.add(label(p.source + ', from ' + p.cores + ' cores across Alaska, British ' +
                     'Columbia, Washington and Oregon.', UI.HINT));
  if (p.highVariability) {
    priorOut.add(label('Carbon here varies a lot between cores, so this site will need an ' +
                       'unusually large campaign. A short pilot would likely save effort.', UI.WARN));
  }
  if (p.indicative) {
    priorOut.add(label('Based on only ' + p.cores + ' cores — treat it as indicative.', UI.WARN));
  }
  priorOut.add(label('A regional average understates how patchy any single site is. If your ' +
                     'own cores come back more variable than this, that is expected.', UI.HINT));
  recompute();
}

function chooseZoneMode(v) {
  State.zoneMode = v;
  clearBuiltZones();
  clearPanel(zoneOut); clearPanel(classRow);
  zoneCountRow.style().set('shown', v === 'embeddings' || v === 'covariates');
  drawRow.style().set('shown',      v === 'draw');
  uploadRow.style().set('shown',    v === 'upload');
  classRow.style().set('shown',     v === 'dynamic' || v === 'copernicus');
  if (v === 'dynamic' || v === 'copernicus') loadLandcoverClasses();
  recompute();
}

// --- drawing zones by hand ------------------------------------------------------

function zoneDraft() {
  var dt = map.drawingTools();
  if (!State.draft) {
    State.draft = ui.Map.GeometryLayer({ geometries: [], name: 'New zone', color: 'C8763C' });
    dt.layers().add(State.draft);
  }
  dt.setSelected(State.draft);
  return State.draft;
}

function startZoneDrawing() {
  clearPanel(zoneOut);
  if (!zoneNameBox.getValue()) {
    zoneOut.add(label('Give the zone a name first, so the export can identify it.', UI.ERROR));
    return;
  }
  var dt = map.drawingTools();
  zoneDraft().geometries().reset([]);
  dt.setShape('polygon');
  dt.draw();
  zoneOut.add(label('Drawing. Click to place corners, double-click to close the shape, ' +
                    'then press "Add the zone I just drew".', UI.NOTE));
}

function addDrawnZone() {
  clearPanel(zoneOut);
  var name = zoneNameBox.getValue();
  if (!name) { zoneOut.add(label('Give the zone a name first.', UI.ERROR)); return; }

  var draft = zoneDraft();
  if (draft.geometries().length() === 0) {
    zoneOut.add(label('Nothing drawn yet. Press "Draw this zone on the map" first.', UI.ERROR));
    return;
  }
  var geom = draft.geometries().get(0);
  State.drawnZones = State.drawnZones || [];
  State.drawnZones.push(ee.Feature(geom, { zone: name }));

  var colour = UI.ZONE_COLOURS[(State.drawnZones.length - 1) % UI.ZONE_COLOURS.length];
  map.addLayer(ee.FeatureCollection([ee.Feature(geom)])
                 .style({ color: colour, fillColor: colour + '55', width: 2 }),
               {}, 'Zone: ' + name);

  draft.geometries().reset([]);
  map.drawingTools().stop();
  zoneNameBox.setValue('');
  listDrawnZones();
  zoneOut.add(label('Added "' + name + '". Name and draw the next one, or press ' +
                    '"Build the zones".', UI.OK));
}

function listDrawnZones() {
  clearPanel(drawnList);
  var n = (State.drawnZones || []).length;
  if (n === 0) { drawnList.add(label('No zones drawn yet.', UI.HINT)); return; }
  drawnList.add(label(n + ' zone' + (n === 1 ? '' : 's') + ' drawn.', UI.OK));
  drawnList.add(ui.Button({
    label: 'Clear the drawn zones', style: { stretch: 'horizontal', margin: '4px 8px' },
    onClick: function () {
      State.drawnZones = []; clearBuiltZones();
      map.layers().reset();
      if (State.aoi) map.addLayer(State.aoi, { color: 'FFFFFF' }, 'Site boundary');
      State.draft = null;
      listDrawnZones(); recompute();
    }
  }));
}

// --- land cover classes ---------------------------------------------------------

function loadLandcoverClasses() {
  clearPanel(classRow);
  if (!State.aoi) { classRow.add(label('Measure the site first.', UI.ERROR)); return; }
  classRow.add(label('Reading the land cover map…', UI.NOTE));

  BCEarth.landcoverClasses(State.aoi, State.zoneMode, function (r, err) {
    clearPanel(classRow);
    if (err) { classRow.add(label(err, UI.ERROR)); return; }

    State.classWidgets = [];
    classRow.add(label('Found ' + r.classes.length + ' cover types inside your site. ' +
                       'Tick the ones to keep as zones, and rename them if you like.', UI.HINT));
    for (var i = 0; i < r.classes.length; i++) {
      var cls = r.classes[i];
      var box = ui.Checkbox({ label: cls.name + '  (' + cls.pixels + ' px)', value: true,
                              style: { margin: '2px 8px' } });
      var ren = ui.Textbox({ placeholder: 'Rename (optional)',
                             style: { stretch: 'horizontal', margin: '0 8px 4px 8px' } });
      classRow.add(box); classRow.add(ren);
      State.classWidgets.push({ box: box, rename: ren, cls: cls });
    }
  });
}

// --- build ----------------------------------------------------------------------

function clearBuiltZones() {
  State.zoneRevision++;
  State.zoneAreas = null; State.zoneImage = null; State.zoneGeoms = {};
  State.zoneKeep = {};
  clearPanel(zoneChoices);
  clearBoundaryExport();
}

function buildZones() {
  clearPanel(zoneOut);
  if (!State.aoi || !State.areaM2) {
    zoneOut.add(label('Measure the site first.', UI.ERROR)); return;
  }
  clearBuiltZones();
  var revision = State.zoneRevision;
  recompute();

  if (State.zoneMode === 'none') {
    zoneOut.add(label('Treating the whole site as one area.', UI.OK));
    return;
  }

  if (State.zoneMode === 'draw') {
    if (!State.drawnZones || State.drawnZones.length === 0) {
      zoneOut.add(label('Draw at least one zone first.', UI.ERROR)); return;
    }
    zoneOut.add(label('Measuring your zones…', UI.NOTE));
    var drawn = State.drawnZones;
    BCEarth.vectorZoneAreas(drawn, function (r, err) {
      if (State.zoneRevision !== revision) return;
      clearPanel(zoneOut);
      if (err) { zoneOut.add(label(err, UI.ERROR)); return; }
      for (var i = 0; i < r.zones.length; i++) {
        State.zoneGeoms[r.zones[i].code] = ee.Feature(drawn[i]).geometry();
      }
      finishZones(r.zones);
    });
    return;
  }

  if (State.zoneMode === 'upload') {
    var asset = zoneAssetBox.getValue(), field = zoneFieldBox.getValue();
    if (!asset || !field) {
      zoneOut.add(label('Give both the asset path and the column that names each zone.',
                        UI.ERROR));
      return;
    }
    zoneOut.add(label('Loading your zones…', UI.NOTE));
    BCEarth.uploadedZones(asset, field, State.aoi, function (r, err) {
      if (State.zoneRevision !== revision) return;
      clearPanel(zoneOut);
      if (err) { zoneOut.add(label(err, UI.ERROR)); return; }
      for (var i = 0; i < r.zones.length; i++) {
        State.zoneGeoms[r.zones[i].code] =
          r.collection.filter(ee.Filter.eq('zone', r.zones[i].name)).geometry();
      }
      map.addLayer(r.collection.style({ color: 'C8763C', fillColor: 'C8763C55', width: 2 }),
                   {}, 'Zones');
      finishZones(r.zones);
    });
    return;
  }

  if (State.zoneMode === 'dynamic' || State.zoneMode === 'copernicus') {
    var picked = [], renames = {}, w;
    for (var j = 0; j < (State.classWidgets || []).length; j++) {
      w = State.classWidgets[j];
      if (w.box.getValue()) {
        picked.push(w.cls);
        if (w.rename.getValue()) renames[w.cls.code] = w.rename.getValue();
      }
    }
    if (picked.length === 0) {
      zoneOut.add(label('Tick at least one cover type to use as a zone.', UI.ERROR)); return;
    }
    var built = BCEarth.landcoverStrata(State.aoi, State.zoneMode, picked, renames);
    State.zoneImage = built.image;
    map.addLayer(built.image, { min: 0, max: Math.max(1, built.zones.length - 1),
                                palette: UI.ZONE_COLOURS }, 'Zones');
    zoneOut.add(label('Measuring each cover type…', UI.NOTE));
    BCEarth.zoneAreas(built.image, built.zones, State.aoi, State.scale, function (z, err) {
      if (State.zoneRevision !== revision) return;
      clearPanel(zoneOut);
      if (err) { zoneOut.add(label(err, UI.ERROR)); return; }
      finishZones(z.zones);
    });
    return;
  }

  // embeddings or satellite covariates
  zoneOut.add(label('Grouping the site…', UI.NOTE));
  BCEarth.clusterStrata(State.aoi, State.zoneMode, zoneCount.getValue(), State.areaM2,
    function (c, err) {
      if (State.zoneRevision !== revision) return;
      clearPanel(zoneOut);
      if (err) { zoneOut.add(label(err, UI.ERROR)); return; }
      State.zoneImage = c.image;
      map.addLayer(c.image, { min: 0, max: c.zones.length - 1, palette: UI.ZONE_COLOURS },
                   'Zones');
      BCEarth.zoneAreas(c.image, c.zones, State.aoi, State.scale, function (z, e2) {
        if (State.zoneRevision !== revision) return;
        if (e2) { zoneOut.add(label(e2, UI.ERROR)); return; }
        finishZones(z.zones);
      });
    });
}

// --- what came back, and which of it to sample ----------------------------------

function finishZones(zones) {
  State.zoneAreas = zones;
  clearBoundaryExport();

  var sum = 0, i;
  for (i = 0; i < zones.length; i++) sum += zones[i].areaM2;
  if (sum / State.areaM2 < 0.95) {
    zoneOut.add(label('Zones cover ' + Math.round(100 * sum / State.areaM2) +
                      '% of the boundary. The rest was not classified and will not ' +
                      'be sampled.', UI.WARN));
  }
  renderZoneChoices();
  recompute();
}

function renderZoneChoices() {
  clearPanel(zoneChoices);
  if (!State.zoneAreas) return;

  zoneChoices.add(label('Which zones do you want to sample?', UI.ASK));
  zoneChoices.add(label('Untick anything you are not sampling — deep water, bare flat, ' +
                        'land. Cores are only shared out among the zones you keep, and ' +
                        'the area of the rest drops out of the calculation.', UI.HINT));

  var i;
  for (i = 0; i < State.zoneAreas.length; i++) {
    (function (z) {
      State.zoneKeep[z.code] = true;
      zoneChoices.add(ui.Checkbox({
        label: z.name + ' — ' + (z.areaM2 / 10000).toFixed(2) + ' ha',
        value: true, style: { margin: '2px 8px' },
        onChange: function (v) {
          State.zoneKeep[z.code] = v;
          clearBoundaryExport();
          recompute();
        }
      }));
    })(State.zoneAreas[i]);
  }
}

function keptZones() {
  if (!State.zoneAreas) return null;
  var out = [], i;
  for (i = 0; i < State.zoneAreas.length; i++) {
    if (State.zoneKeep[State.zoneAreas[i].code]) out.push(State.zoneAreas[i]);
  }
  return out;
}

function recompute() {
  clearPanel(designOut); clearPanel(mathPanel);
  State.design = null;
  if (!State.areaM2 || !State.prior) return;

  var E = moeSelect.getValue(), conf = confSelect.getValue();
  var result, i;

  var kept = keptZones();
  if (kept && kept.length === 0) {
    designOut.add(label('Every zone is unticked, so there is nothing to sample.', UI.ERROR));
    return;
  }

  if (kept && kept.length > 1) {
    var strata = [];
    for (i = 0; i < kept.length; i++) {
      strata.push({ name: kept[i].name, code: kept[i].code, areaM2: kept[i].areaM2,
                    mean: State.prior.mean, sd: State.prior.sd });
    }
    result = BCStats.stratifiedSampleSize({
      strata: strata, plotM2: CONFIG.PLOT_M2, confidence: conf, marginOfError: E,
      allocation: allocSelect.getValue(), minPerStratum: CONFIG.MIN_CORES_PER_ZONE });
  } else {
    // one zone kept, or none defined: sample that area alone
    var area = (kept && kept.length === 1) ? kept[0].areaM2 : State.areaM2;
    result = BCStats.srsSampleSize({
      areaM2: area, plotM2: CONFIG.PLOT_M2, mean: State.prior.mean,
      sd: State.prior.sd, confidence: conf, marginOfError: E });
  }

  if (!result.ok) { designOut.add(label(result.reason, UI.ERROR)); return; }
  State.design = result;

  var formula = result.nAllocated || result.n;
  var total   = Math.max(formula, BCStats.MIN_USABLE_CORES);
  State.cores = total;

  designOut.add(label(total + ' cores', UI.RESULT));
  designOut.add(label('to know the mean within ' + Math.round(E * 100) + '%, ' +
                      Math.round(conf * 100) + '% of the time.', UI.NOTE));

  if (formula < BCStats.MIN_USABLE_CORES) {
    designOut.add(label('The formula asks for only ' + formula + '. Collect at least ' +
                        BCStats.MIN_USABLE_CORES + ' regardless — with fewer you cannot ' +
                        'measure how variable the site is, so the precision you actually ' +
                        'achieved can never be checked.', UI.WARN));
    if (State.prior && State.prior.cv < 0.15) {
      designOut.add(label('The variability behind this is very low (CV ' +
                          State.prior.cv.toFixed(2) + '). Check you picked the right ' +
                          'ecosystem and depth — a figure this low is unusual for ' +
                          'coastal sediment.', UI.WARN));
    }
  }

  if (result.strata) {
    for (i = 0; i < result.strata.length; i++) {
      designOut.add(label('   ' + result.strata[i].name + ': ' + result.strata[i].cores +
                          ' cores' + (result.strata[i].flooredToMinimum
                            ? '  (raised to the ' + CONFIG.MIN_CORES_PER_ZONE + '-core minimum)' : ''),
                          UI.NOTE));
    }
  }

  // feasibility
  var forGeom = result.strata
    ? result.strata.map(function (s) { return { name: s.name, areaM2: s.areaM2, cores: s.cores }; })
    : [{ name: 'Whole site', areaM2: State.areaM2, cores: result.n }];
  var f = BCGeom.feasibility({ plotM2: CONFIG.PLOT_M2,
                               minPerStratum: CONFIG.MIN_CORES_PER_ZONE, strata: forGeom });
  for (i = 0; i < f.problems.length; i++) designOut.add(label(f.problems[i], UI.ERROR));
  for (i = 0; i < f.warnings.length; i++) designOut.add(label(f.warnings[i], UI.WARN));

  // what a different precision target would cost
  var alt = [0.10, 0.20, 0.30], line = [];
  for (i = 0; i < alt.length; i++) {
    var a = BCStats.srsSampleSize({ areaM2: State.areaM2, plotM2: CONFIG.PLOT_M2,
              mean: State.prior.mean, sd: State.prior.sd, confidence: conf, marginOfError: alt[i] });
    if (a.ok) line.push(Math.round(alt[i] * 100) + '% needs ' + a.n);
  }
  designOut.add(label('For comparison: ' + line.join(', ') + '.', UI.HINT));

  // the working
  mathPanel.add(label('Cochran, with the finite-population correction:', UI.HINT));
  mathPanel.add(label('n = z² N CV² / ((N−1) E² + z² CV²)', UI.NOTE));
  mathPanel.add(label('z = ' + result.z.toFixed(4) + '   N = ' + result.N +
                      '   CV = ' + (result.cv || Math.sqrt(result.V)).toFixed(3) +
                      '   E = ' + E, UI.NOTE));
  mathPanel.add(label('Each core stands for a ' + CONFIG.PLOT_M2 + ' m² plot, so N is the ' +
                      'site area divided by that.', UI.HINT));
  if (result.allocation) {
    mathPanel.add(label('Allocation: ' + result.allocation + '. Cores are rounded up per zone ' +
                        'and floored at ' + CONFIG.MIN_CORES_PER_ZONE +
                        ', so the total runs slightly above the formula.', UI.HINT));
  }
  var a1 = ui.Label('Full derivation — Appendix A', UI.LINK);
  a1.setUrl(CONFIG.linkTo('appendixA'));
  mathPanel.add(a1);
}

function placeCores() {
  clearPanel(placeOut);
  if (!State.design) { placeOut.add(label('Work through Steps 1 to 4 first.', UI.ERROR)); return; }

  var bearing = parseFloat(bearingBox.getValue());
  var common = {
    plotM2: CONFIG.PLOT_M2, layout: State.layout,
    transectCount: transectBox.getValue(),
    bearingDeg: isNaN(bearing) ? null : bearing,
    subsamples: subsampleBox.getValue(), compositeRadiusM: State.compositeRadius
  };

  var kept = keptZones();
  var zoned = State.design.strata && kept && kept.length > 1;

  // whole site, one zone
  if (!zoned) {
    var cores = State.cores || State.design.n;
    var area  = (kept && kept.length === 1) ? kept[0].areaM2 : State.areaM2;
    placeOut.add(label('Placing ' + cores + ' cores…', UI.NOTE));

    withZoneGeometry(kept && kept.length === 1 ? kept[0] : null, State.aoi, function (g, gErr) {
      if (gErr) { clearPanel(placeOut); placeOut.add(label(gErr, UI.ERROR)); return; }
      var o = { geometry: g, areaM2: area, cores: cores, zoneIndex: 0,
                compositeCount: Math.max(5, Math.round(cores / 2)) };
      for (var k in common) { if (common.hasOwnProperty(k)) o[k] = common[k]; }
      BCPlace.place(o, function (r, err) {
        clearPanel(placeOut);
        if (err) { placeOut.add(label(err, UI.ERROR)); return; }
        finishPlacement([{ zone: 'Whole site', result: r }]);
      });
    });
    return;
  }

  // several zones: place inside each, one after the other
  placeOut.add(label('Placing cores zone by zone…', UI.NOTE));
  var results = [], idx = 0;
  var strata = State.design.strata;
  var revision = State.zoneRevision;

  function next() {
    if (State.zoneRevision !== revision) return;
    if (idx >= strata.length) { clearPanel(placeOut); finishPlacement(results); return; }
    var st = strata[idx];
    var zoneRecord = null, i;
    for (i = 0; i < kept.length; i++) { if (kept[i].code === st.code) zoneRecord = kept[i]; }

    withZoneGeometry(zoneRecord, State.aoi, function (g, gErr) {
      if (gErr) {
        results.push({ zone: st.name, error: gErr });
        idx++; next(); return;
      }
      var o = { geometry: g, areaM2: st.areaM2, cores: st.cores, zoneIndex: idx,
                compositeCount: Math.max(5, Math.round(st.cores / 2)) };
      for (var k in common) { if (common.hasOwnProperty(k)) o[k] = common[k]; }
      BCPlace.place(o, function (r, err) {
        results.push(err ? { zone: st.name, error: err } : { zone: st.name, result: r });
        idx++; next();
      });
    });
  }
  next();
}

// Drawn and uploaded zones already have a polygon. A zone that came out of a
// raster is converted to a polygon on demand for placement or download.
function withZoneGeometry(zoneRecord, fallback, callback) {
  if (!zoneRecord) { callback(fallback, null); return; }
  if (State.zoneGeoms[zoneRecord.code]) {
    callback(State.zoneGeoms[zoneRecord.code], null); return;
  }
  if (!State.zoneImage) {
    callback(null, 'Could not find the boundary for ' + zoneRecord.name + '.'); return;
  }

  var revision = State.zoneRevision;
  BCEarth.zoneGeometry(State.zoneImage, zoneRecord.code, State.aoi, State.scale,
    function (g, err) {
      if (State.zoneRevision !== revision) { callback(null, 'The zones changed.'); return; }
      if (err || !g) { callback(null, 'Could not outline ' + zoneRecord.name + ': ' + err); return; }
      State.zoneGeoms[zoneRecord.code] = g;
      callback(g, null);
    });
}

function finishPlacement(results) {
  clearPanel(placeOut);

  var all = [], collections = [], totalPlaced = 0, worst = Infinity, i, j;

  for (i = 0; i < results.length; i++) {
    var r = results[i];
    if (r.error) { placeOut.add(label(r.zone + ': ' + r.error, UI.ERROR)); continue; }

    var fc = BCPlace.toFeatures(r.result.points, r.zone, r.result.layout, CONFIG.PLOT_M2);
    if (!fc) { placeOut.add(label(r.zone + ': no positions were produced.', UI.ERROR)); continue; }

    collections.push(fc);
    totalPlaced += r.result.placed;
    if (r.result.minSeparationM < worst) worst = r.result.minSeparationM;
    for (j = 0; j < r.result.points.length; j++) all.push(r.result.points[j]);

    if (results.length > 1) {
      placeOut.add(label(r.zone + ': ' + r.result.placed + ' cores', UI.NOTE));
    }
    if (r.result.shortfall) placeOut.add(label(r.zone + ': ' + r.result.shortfall, UI.WARN));
    if (r.result.meta && r.result.meta.note) {
      placeOut.add(label(r.zone + ': ' + r.result.meta.note, UI.WARN));
    }
  }

  if (collections.length === 0) {
    placeOut.add(label('No cores could be placed.', UI.ERROR));
    return;
  }

  var merged = collections[0];
  for (i = 1; i < collections.length; i++) merged = merged.merge(collections[i]);

  State.features = merged;
  State.placed = { placed: totalPlaced, points: all, layout: State.layout,
                   minSeparationM: worst };

  map.addLayer(merged, { color: 'FFCC00' }, 'Core locations');
  placeOut.add(label(totalPlaced + ' cores placed.', UI.OK));
  placeOut.add(label('Closest pair ' + worst.toFixed(0) + ' m apart; no two cores share ' +
                     'a plot.', UI.HINT));
  writeMethods();
}

function writeMethods() {
  clearPanel(methodsOut);
  if (!State.design || !State.placed) return;

  var d = State.design, E = moeSelect.getValue(), conf = confSelect.getValue();
  var cv = d.cv || Math.sqrt(d.V);
  var layoutName = '';
  for (var i = 0; i < BCLayout.LAYOUTS.length; i++) {
    if (BCLayout.LAYOUTS[i].id === State.layout) layoutName = BCLayout.LAYOUTS[i].label.toLowerCase();
  }

  var text =
    'Sediment cores (n = ' + State.placed.placed + ') were located across ' +
    ((d.totalAreaM2 || State.areaM2) / 10000).toFixed(1) + ' ha of ' +
    State.ecosystem.toLowerCase() +
    ' using a ' + layoutName + ' layout' +
    (d.strata ? ', allocated across ' + d.strata.length + ' zones by ' + d.allocation +
                ' allocation' : '') +
    '. Sample size was calculated by Cochran\'s formula with a finite-population ' +
    'correction, for a relative margin of error of ' + Math.round(E * 100) + '% at ' +
    Math.round(conf * 100) + '% confidence, taking each core to represent a ' +
    CONFIG.PLOT_M2 + ' m\u00B2 plot. Expected variability (CV = ' + cv.toFixed(2) +
    ') came from ' + (State.prior.source || 'the site') +
    '. Cores were placed at least ' + BCGeom.minSpacing(CONFIG.PLOT_M2) +
    ' m apart. Achieved precision should be recalculated from the collected cores ' +
    'before the estimate is reported.';

  methodsOut.add(ui.Label(text, { fontSize: '11px', color: '#26343A', margin: '2px 8px' }));
}

function exportCores() {
  clearPanel(exportOut);
  if (!State.features) { exportOut.add(label('Place the cores first.', UI.ERROR)); return; }
  var fmt = formatSelect.getValue();
  var url = State.features.getDownloadURL({
    format: fmt, filename: 'blue_carbon_cores_' + (new Date()).getTime() });
  var l = ui.Label('Download ' + fmt, UI.LINK);
  l.setUrl(url);
  exportOut.add(l);
  exportOut.add(label('Includes core_id, zone, layout, plot size and coordinates.', UI.HINT));
}

function clearBoundaryExport() {
  boundaryExportRequest++;
  clearPanel(boundaryExportOut);
}

function exportBoundary() {
  clearBoundaryExport();
  if (!State.aoi || !State.areaM2) {
    boundaryExportOut.add(label('Measure the site first.', UI.ERROR)); return;
  }

  var scope = boundaryExportSelect.getValue();
  var zones = [];
  if (scope !== 'aoi') {
    if (!State.zoneAreas || !State.zoneAreas.length) {
      boundaryExportOut.add(label('Build the zones in Step 2 first.', UI.ERROR)); return;
    }
    zones = scope === 'selected_zones' ? keptZones() : State.zoneAreas;
    if (!zones.length) {
      boundaryExportOut.add(label('Tick at least one zone in Step 2.', UI.ERROR)); return;
    }
  }

  var features = [];
  if (scope !== 'selected_zones') {
    features.push(ee.Feature(State.aoi, { feature_type: 'aoi' }));
  }
  var request = boundaryExportRequest;
  var filename = scope === 'aoi' ? 'aoi' : 'blue_carbon_' + scope;
  boundaryExportOut.add(label('Preparing GeoJSON…', UI.NOTE));

  function showDownload(url, err) {
    if (request !== boundaryExportRequest) return;
    clearPanel(boundaryExportOut);
    if (!url) {
      boundaryExportOut.add(label(err ? 'Could not prepare the GeoJSON. ' + err :
                                 'Could not prepare the GeoJSON.', UI.ERROR));
      return;
    }
    var link = label('Download GeoJSON', UI.LINK);
    link.setUrl(url);
    boundaryExportOut.add(link);
  }

  function addZone(index) {
    if (request !== boundaryExportRequest) return;
    if (index === zones.length) {
      try {
        ee.FeatureCollection(features).getDownloadURL({
          format: 'GeoJSON', filename: filename, callback: showDownload
        });
      } catch (err) {
        showDownload(null, err.message);
      }
      return;
    }
    var zone = zones[index];
    withZoneGeometry(zone, null, function (geometry, err) {
      if (request !== boundaryExportRequest) return;
      if (err) { showDownload(null, err); return; }
      features.push(ee.Feature(geometry, { feature_type: 'zone', zone: zone.name,
                                          zone_code: zone.code }));
      addZone(index + 1);
    });
  }

  addZone(0);
}


if (typeof module !== 'undefined') {
  module.exports = { CONFIG: CONFIG, BCStats: BCStats, BCGeom: BCGeom, BCTest: BCTest,
                     BCEarth: BCEarth, BCLayout: BCLayout, BCPlace: BCPlace };
}


// =================================================================================
// === RUN THE SELF TEST ===========================================================
// === Last, so that every section above it is defined.
// =================================================================================

if (RUN_SELF_TEST) { BCTest.run(); }
