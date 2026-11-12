# AlphaEarth-weighted Janousek priors and a local Matérn update.
# Run from the workshop folder. See data/covariate_transfer/README.md.

settings <- list(
  input_folder = "data/covariate_transfer",
  output_folder = "outputs/script3_covariate_transfer",
  aoi_file = "data/boundary/aoi.geojson",
  mapping_file = "outputs/script1/mapping_data.csv",
  planning_cores_file = "data/field/blue_carbon_cores_1794411291716.csv",
  analysis_crs = "EPSG:32610",
  embedding_year = 2024,
  cell_size_m = 250,
  core_match_tolerance_m = 2,
  bandwidth = 0.10,
  local_model_file = "outputs/script2_local/mapping_run.rds",
  default_range_m = 5000,
  default_smoothness = 0.5,
  default_observation_log_sd = 0.1
)

source("helpers/local_mapping_spatial.R")
source("helpers/local_mapping_models.R")

depths <- c(30, 50, 100)
stock_cols <- paste0("0-", depths)
bands <- paste0("var", seq_len(64))
point_cols <- c("Core ID", "Latitude", "Longitude", stock_cols, bands, "embedding_year")

write_table <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "")
}

numbers <- function(x, label) {
  value <- trimws(as.character(x))
  result <- suppressWarnings(as.numeric(value))
  blank <- is.na(value) | value %in% c("", "NA", "NaN", "-9999")
  if (any(!blank & !is.finite(result))) stop("Invalid numeric values in ", label, ".")
  result[blank | result %in% -9999] <- NA_real_
  result
}

check_coordinates <- function(data, label) {
  data$Latitude <- numbers(data$Latitude, paste(label, "Latitude"))
  data$Longitude <- numbers(data$Longitude, paste(label, "Longitude"))
  if (anyNA(data$Latitude) || anyNA(data$Longitude) ||
      any(abs(data$Latitude) > 90 | abs(data$Longitude) > 180)) {
    stop("Invalid WGS84 coordinates in ", label, ".")
  }
  data
}

match_core_rows <- function(cores, points, label) {
  points <- check_coordinates(points, label)
  points_sf <- sf::st_as_sf(points, coords = c("Longitude", "Latitude"), crs = 4326)
  points_xy <- sf::st_coordinates(sf::st_transform(points_sf, sf::st_crs(cores)))
  distances <- distance_matrix(sf::st_coordinates(cores), points_xy)
  ids <- cores[["Core ID"]]
  matches <- vapply(seq_along(ids), function(i) {
    nearby <- which(distances[i, ] <= settings$core_match_tolerance_m)
    same_id <- nearby[points[["Core ID"]][nearby] %in% ids[i]]
    if (length(same_id) == 1L) return(same_id)
    if (length(nearby) != 1L) {
      stop("Expected one matching location in ", label, " for ", ids[i],
           " within ", settings$core_match_tolerance_m, " m. Found ", length(nearby), ".")
    }
    nearby
  }, integer(1))
  if (anyDuplicated(matches)) stop("More than one local core matched a single row in ", label, ".")
  points[matches, , drop = FALSE]
}

core_stocks <- function(cores) {
  columns <- paste0("C Stock 0-", depths, "cm")
  if (!all(c(columns, "Stock coverage (cm)") %in% names(cores))) {
    stop("Run script_1.R to generate the depth-specific stocks and stock coverage.")
  }
  coverage <- numbers(cores[["Stock coverage (cm)"]], "Script 1 stock coverage")
  stocks <- setNames(lapply(columns, function(column) numbers(cores[[column]], column)), stock_cols)
  # A short core can supply a valid stock at every target depth.
  for (i in seq_along(depths)) {
    stock <- stocks[[i]]
    stock[is.na(coverage) | coverage <= 0 | is.na(stock) | stock < 0] <- NA_real_
    stocks[[i]] <- stock
  }
  as.data.frame(stocks, check.names = FALSE)
}

prepare_uploads <- function(inputs, stocks) {
  cores <- inputs$cores
  coordinates <- sf::st_drop_geometry(cores)[, c("Latitude", "Longitude")]
  if (!is.null(settings$planning_cores_file)) {
    planned <- utils::read.csv(settings$planning_cores_file, check.names = FALSE)
    if (!all(c("lat", "lon") %in% names(planned)) || !nrow(planned)) {
      stop("The planning export must contain points with lat and lon columns.")
    }
    coordinates <- match_core_rows(cores, data.frame(Latitude = planned$lat, Longitude = planned$lon),
                                   settings$planning_cores_file)
  }
  upload <- file.path(settings$input_folder, "upload")
  dir.create(upload, recursive = TRUE, showWarnings = FALSE)
  scratch <- tempfile("aoi_shapefile_")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  aoi <- sf::st_transform(inputs$aoi, 4326)
  aoi$type <- "aoi"
  sf::st_write(aoi, file.path(scratch, "aoi.shp"), quiet = TRUE)
  zip_path <- file.path(scratch, "aoi.zip")
  if (utils::zip(zip_path, list.files(scratch, full.names = TRUE), flags = "-j") != 0 ||
      !file.exists(zip_path) || !file.copy(zip_path, file.path(upload, "aoi.zip"), overwrite = TRUE)) {
    stop("Could not create the Earth Engine Shapefile upload ZIP.")
  }
  names(stocks) <- paste0("stock_", depths)
  write_table(data.frame(core_id = cores[["Core ID"]], latitude = coordinates$Latitude,
                         longitude = coordinates$Longitude, stocks),
              file.path(upload, "local_points_upload.csv"))
  message("Prepared the AOI and ", nrow(cores), " matched cores in ", upload, ".")
}

embedding_stems <- function(kind, legacy_prefix) {
  year <- settings$embedding_year
  c(paste0("alphaearth_", year, "_", kind), paste0(legacy_prefix, "_", year))
}

read_embedding_table <- function(kind, legacy_prefix) {
  stems <- file.path(settings$input_folder, embedding_stems(kind, legacy_prefix))
  candidates <- paste0(rep(stems, each = 2L), c(".csv", ".xlsx"))
  path <- candidates[file.exists(candidates)]
  if (length(path) != 1L) {
    stop("Expected one Earth Engine export for ", kind, ". Keep one of: ",
         paste(basename(candidates), collapse = ", "), ".")
  }
  table <- if (endsWith(path, ".xlsx")) {
    as.data.frame(readxl::read_excel(path, col_types = "text", .name_repair = "minimal"))
  } else {
    utils::read.csv(path, check.names = FALSE, colClasses = "character")
  }
  missing <- setdiff(point_cols, names(table))
  if (!nrow(table) || anyDuplicated(names(table)) || length(missing)) {
    stop("Invalid or empty embedding table ", path, ". Missing columns: ",
         paste(missing, collapse = ", "), ".")
  }
  table <- check_coordinates(table[, point_cols], path)
  table[["Core ID"]] <- trimws(table[["Core ID"]])
  if (anyNA(table[["Core ID"]]) || any(!nzchar(table[["Core ID"]]))) {
    stop("Missing Core ID in ", path, ".")
  }
  for (column in c(stock_cols, bands, "embedding_year")) {
    table[[column]] <- numbers(table[[column]], paste(path, column))
  }
  if (!all(table$embedding_year %in% settings$embedding_year)) {
    stop("Embedding year does not match settings in ", path, ".")
  }
  if (any(as.matrix(table[, stock_cols]) < 0, na.rm = TRUE)) {
    stop("Negative stocks in ", path, ". Expected kg C/m².")
  }
  table
}

read_embedding_raster <- function() {
  stems <- embedding_stems("aoi", "aoi_embeddings")
  candidates <- lapply(stems, function(stem) {
    list.files(settings$input_folder, paste0("^", stem, "(-[^.]+)?[.]tif$"),
               full.names = TRUE)
  })
  present <- lengths(candidates) > 0L
  if (!any(present)) stop("Missing ", stems[1L], ".tif in ", settings$input_folder, ".")
  if (sum(present) > 1L) stop("Keep one AOI raster export per embedding year.")
  files <- candidates[[which(present)]]
  if (length(files) > 1L && paste0(stems[present], ".tif") %in% basename(files)) {
    stop("Keep either the single AOI GeoTIFF or its tiles in the input folder.")
  }
  # Earth Engine preserves band order even if the GeoTIFF has generic names.
  for (file in files) {
    tile <- terra::rast(file)
    if (terra::nlyr(tile) != length(bands) || !nzchar(terra::crs(tile)) ||
        (any(names(tile) %in% bands) && !identical(names(tile), bands))) {
      stop("Expected 64 ordered AlphaEarth bands and a CRS in ", file, ".")
    }
  }
  raster <- if (length(files) == 1L) terra::rast(files) else terra::vrt(files)
  if (terra::nlyr(raster) != length(bands)) stop("The AOI tiles did not produce a 64-band mosaic.")
  names(raster) <- bands
  raster
}

normalize_embeddings <- function(values) {
  values <- as.matrix(values)
  values[!is.finite(values) | abs(values) > 1 + 1e-6] <- NA_real_
  norm <- sqrt(rowSums(values^2))
  valid <- is.finite(norm) & norm > 0
  values[valid, ] <- values[valid, , drop = FALSE] / norm[valid]
  values[!valid, ] <- NA_real_
  values
}

transfer_prior <- function(queries, donors, stocks) {
  result <- data.frame(prior_mean_kg_c_m2 = rep(NA_real_, nrow(queries)),
                       prior_sd_kg_c_m2 = NA_real_, effective_donors = NA_real_,
                       max_cosine_similarity = NA_real_)
  if (length(stocks) < 2L) return(result)
  valid <- which(!is.na(queries[, 1]))
  if (!length(valid)) return(result)
  center <- mean(stocks)
  deviation <- stocks - center
  for (start in seq.int(1L, length(valid), by = 256L)) {
    rows <- valid[start:min(length(valid), start + 255L)]
    similarity <- tcrossprod(queries[rows, , drop = FALSE], donors)
    similarity[] <- pmax(-1, pmin(1, similarity))
    best <- apply(similarity, 1L, max)
    weights <- exp((similarity - best) / settings$bandwidth)
    weights <- weights / rowSums(weights)
    squared_weights <- rowSums(weights^2)
    mean_deviation <- drop(weights %*% deviation)
    variance <- pmax(0, drop(weights %*% deviation^2) - mean_deviation^2)
    correction <- 1 - squared_weights
    variance <- ifelse(correction > .Machine$double.eps, variance / correction, NA_real_)
    result[rows, ] <- list(center + mean_deviation, sqrt(variance), 1 / squared_weights, best)
  }
  result
}

log_prior <- function(prior) {
  mean <- prior$prior_mean_kg_c_m2
  sd <- prior$prior_sd_kg_c_m2
  valid <- is.finite(mean) & mean > 0 & is.finite(sd) & sd > 0
  variance <- rep(NA_real_, length(mean))
  variance[valid] <- log1p((sd[valid] / mean[valid])^2)
  list(mean = log(ifelse(valid, mean, NA_real_)) - variance / 2, sd = sqrt(variance))
}

mapping_parameters <- function(inputs) {
  defaults <- data.frame(depth_cm = depths, range_m = settings$default_range_m,
                         smoothness = settings$default_smoothness,
                         observation_log_sd = settings$default_observation_log_sd,
                         parameter_source = "Script 3 defaults")
  path <- settings$local_model_file
  if (is.null(path) || !file.exists(path)) return(defaults)
  run <- readRDS(path)
  checksums <- unname(tools::md5sum(c(settings$aoi_file, settings$mapping_file)))
  if (!identical(run$input_md5, checksums) ||
      !isTRUE(sf::st_crs(run$settings$analysis_crs) == sf::st_crs(inputs$aoi))) {
    stop("Script 2 used a different AOI, core table or CRS. Rerun it or set local_model_file = NULL.")
  }
  columns <- c("range_m", "smoothness", "observation_log_sd")
  if (!all(c("depth_cm", columns) %in% names(run$models)) ||
      anyDuplicated(run$models$depth_cm)) {
    stop("Missing or duplicate Matérn model settings in ", path, ".")
  }
  matched <- match(depths, run$models$depth_cm)
  have_fit <- !is.na(matched)
  defaults[have_fit, columns] <- run$models[matched[have_fit], columns]
  defaults$parameter_source[have_fit] <- "Script 2 selected Matérn"
  defaults
}

update_prior <- function(local_prior, grid_prior, local_stock, core_xy, grid_xy, model) {
  observed <- log_prior(local_prior)
  target <- log_prior(grid_prior)
  used <- is.finite(observed$mean) & is.finite(observed$sd) &
    is.finite(local_stock) & local_stock > 0
  predicted <- is.finite(target$mean) & is.finite(target$sd)
  posterior <- data.frame(posterior_mean_kg_c_m2 = rep(NA_real_, nrow(grid_prior)),
                          posterior_sd_kg_c_m2 = NA_real_)
  if (any(used) && any(predicted)) {
    locations <- core_xy[used, , drop = FALSE]
    correlation <- matern_correlation(distance_matrix(locations, locations),
                                      model$range_m, model$smoothness)
    covariance <- tcrossprod(observed$sd[used]) * correlation
    diag(covariance) <- diag(covariance) + model$observation_log_sd^2
    lower <- tryCatch(t(chol(covariance)), error = function(e) {
      stop("Local covariance could not be factored at ", model$depth_cm,
           " cm. Check coincident cores and the observation error.")
    })
    residual <- forwardsolve(lower, log(local_stock[used]) - observed$mean[used])
    target_rows <- which(predicted)
    for (start in seq.int(1L, length(target_rows), by = 256L)) {
      rows <- target_rows[start:min(length(target_rows), start + 255L)]
      distances <- distance_matrix(locations, grid_xy[rows, , drop = FALSE])
      cross <- outer(observed$sd[used], target$sd[rows]) *
        matern_correlation(distances, model$range_m, model$smoothness)
      conditioned <- forwardsolve(lower, cross)
      moments <- log_stock_moments(
        target$mean[rows] + drop(crossprod(residual, conditioned)),
        pmax(0, target$sd[rows]^2 - colSums(conditioned^2))
      )
      posterior[rows, ] <- unname(as.list(moments))
    }
  }
  list(posterior = posterior, used = used)
}

run_transfer <- function(inputs, stocks) {
  reference <- read_embedding_table("janousek_cores", "janousek_embeddings")
  if (anyDuplicated(reference[, c("Core ID", "Latitude", "Longitude")])) {
    stop("The Janousek export has duplicate cores at the same coordinates.")
  }
  local <- match_core_rows(inputs$cores, read_embedding_table("local_cores", "local_embeddings"),
                           "local embeddings")
  local_export <- as.matrix(local[, stock_cols])
  stock_values <- as.matrix(stocks)
  supplied <- is.finite(local_export)
  if (any(supplied & (is.na(stock_values) |
      abs(local_export - stock_values) > 1e-6 * pmax(1, stock_values)))) {
    stop("Local export stocks disagree with Script 1. Regenerate the upload and export.")
  }
  raster <- read_embedding_raster()
  xy_wgs84 <- sf::st_coordinates(sf::st_transform(inputs$points, 4326))
  raster_points <- sf::st_transform(inputs$points, terra::crs(raster))
  values <- terra::extract(raster, terra::vect(raster_points), ID = FALSE)
  grid_covariates <- normalize_embeddings(values)
  values[is.na(grid_covariates[, 1]), ] <- NA_real_
  grid_ids <- paste0("cell_", inputs$points$cell_id)
  pred_points <- data.frame("Core ID" = grid_ids, Latitude = xy_wgs84[, 2],
                            Longitude = xy_wgs84[, 1], "0-30" = NA_real_,
                            "0-50" = NA_real_, "0-100" = NA_real_, values,
                            embedding_year = settings$embedding_year, check.names = FALSE)
  local_covariates <- normalize_embeddings(local[, bands])
  donor_covariates <- normalize_embeddings(reference[, bands])
  core_xy <- sf::st_coordinates(inputs$cores)
  grid_xy <- sf::st_coordinates(inputs$points)
  models <- mapping_parameters(inputs)
  if (any(!is.finite(models$range_m) | models$range_m <= 0 |
          !is.finite(models$smoothness) | models$smoothness <= 0 |
          !is.finite(models$observation_log_sd) | models$observation_log_sd < 0)) {
    stop("Matérn parameters must be finite, with positive range and smoothness and nonnegative error.")
  }
  folder <- settings$output_folder
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  summaries <- vector("list", length(depths))
  local_priors <- summaries
  predictions <- summaries
  area <- inputs$cells$area_m2
  data_source <- inputs$cores$data_source
  if (is.null(data_source)) data_source <- rep(NA_character_, nrow(inputs$cores))
  source_label <- paste(unique(stats::na.omit(data_source)), collapse = ", ")
  for (i in seq_along(depths)) {
    depth <- depths[i]
    message("Estimating transfer model for 0–", depth, " cm.")
    donor_stock <- reference[[stock_cols[i]]]
    donors <- is.finite(donor_stock) & !is.na(donor_covariates[, 1])
    donor_vectors <- donor_covariates[donors, , drop = FALSE]
    local_prior <- transfer_prior(local_covariates, donor_vectors, donor_stock[donors])
    grid_prior <- transfer_prior(grid_covariates, donor_vectors, donor_stock[donors])
    update <- update_prior(local_prior, grid_prior, stocks[[i]], core_xy, grid_xy, models[i, ])
    layers <- cbind(grid_prior, update$posterior)
    export_mapping_raster(inputs, layers, file.path(folder, paste0("transfer_", depth, "cm.tif")))
    uncertainty <- layers[, c("prior_sd_kg_c_m2", "posterior_sd_kg_c_m2")]
    export_mapping_raster(inputs, uncertainty,
                          file.path(folder, paste0("transfer_uncertainty_", depth, "cm.tif")))

    local_priors[[i]] <- data.frame("Core ID" = inputs$cores[["Core ID"]],
                                    "Plot ID" = inputs$cores[["Plot ID"]], depth_cm = depth,
                                    Latitude = inputs$cores$Latitude, Longitude = inputs$cores$Longitude,
                                    data_source = data_source, observed_kg_c_m2 = stocks[[i]],
                                    has_embedding = !is.na(local_covariates[, 1]),
                                    used_in_update = update$used, local_prior, check.names = FALSE)
    predictions[[i]] <- data.frame("Core ID" = grid_ids, cell_id = inputs$cells$cell_id,
                                   depth_cm = depth, Latitude = pred_points$Latitude,
                                   Longitude = pred_points$Longitude, area_m2 = area,
                                   layers, check.names = FALSE)
    summaries[[i]] <- do.call(rbind, lapply(c("prior", "posterior"), function(surface) {
      density <- layers[[paste0(surface, "_mean_kg_c_m2")]]
      covered <- is.finite(density)
      covered_area <- sum(area[covered])
      total <- if (any(covered)) sum(density[covered] * area[covered]) / 1000 else NA_real_
      data.frame(models[i, ], surface = surface, embedding_year = settings$embedding_year,
                 data_source = source_label, bandwidth = settings$bandwidth,
                 cell_size_m = settings$cell_size_m,
                 aoi_ha = inputs$area_m2 / 10000, covered_ha = covered_area / 10000,
                 coverage_fraction = min(1, covered_area / inputs$area_m2),
                 n_reference_cores = sum(donors), n_local_cores_used = sum(update$used),
                 mean_covered_kg_c_m2 = if (any(covered)) 1000 * total / covered_area else NA_real_,
                 covered_stock_Mg_C = total,
                 aoi_stock_Mg_C = if (all(covered)) total else NA_real_)
    }))
  }
  write_table(pred_points, file.path(settings$input_folder, "pred_points.csv"))
  write_table(do.call(rbind, predictions), file.path(folder, "transfer_predictions.csv"))
  write_table(do.call(rbind, local_priors), file.path(folder, "transfer_local_priors.csv"))
  summary <- do.call(rbind, summaries)
  write_table(summary, file.path(folder, "transfer_area_summary.csv"))
  missing_posterior <- summary$surface == "posterior" & summary$coverage_fraction == 0
  if (any(missing_posterior)) {
    warning("No posterior coverage at depth(s) ", paste(summary$depth_cm[missing_posterior], collapse = ", "),
            " cm. Check embedding coverage, reference stocks and usable local cores.", call. = FALSE)
  }
  message("Wrote prediction inputs and transfer outputs to ", settings$input_folder, " and ", folder, ".")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--prepare")) stop("Usage: Rscript Script_3_covariate_transfer.R [--prepare]")
stopifnot(is.finite(settings$core_match_tolerance_m), settings$core_match_tolerance_m >= 0,
          is.finite(settings$bandwidth), settings$bandwidth > 0)
inputs <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                              settings$analysis_crs, settings$cell_size_m)
stocks <- core_stocks(inputs$cores)
if (identical(args, "--prepare")) prepare_uploads(inputs, stocks) else run_transfer(inputs, stocks)
