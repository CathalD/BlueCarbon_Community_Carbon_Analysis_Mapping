# Script 2 | Local sediment carbon mapping
# Run from the workshop folder. Requires sf, terra and ggplot2.
# Edit these settings for a new site, then run Script_2.2_Bayesian.R if desired.

settings <- list(
  aoi_file = "data/boundary/aoi.geojson",
  mapping_file = "outputs/script1/mapping_data.csv",
  analysis_crs = "EPSG:32610",   # WGS84 / UTM zone 10N for this coastal AOI
  cell_size_m = 30,
  depths_cm = c(30, 50, 100),
  # Candidate levels from the November workshop. Inf uses all training cores.
  idw_powers = c(0, 0.5, 1, 2, 3, 4),
  idw_neighbours = c(2, 3, 4, Inf),
  matern_scales_m = c(5000, 10000, 20000, 40000, 80000),
  matern_smoothness = c(0.5, 1.5, 2.5, 10),
  matern_nugget_fraction = 0.05  # Fraction of log-stock variance attributed to core-level noise
)
output_folder <- "outputs/script2_local"

source("helpers/local_mapping_spatial.R")
source("helpers/local_mapping_models.R")
library(ggplot2)

stopifnot(
  length(settings$depths_cm) > 0, all(is.finite(settings$depths_cm)),
  all(settings$depths_cm > 0), !anyDuplicated(settings$depths_cm),
  length(settings$idw_powers) > 0, all(is.finite(settings$idw_powers) & settings$idw_powers >= 0),
  length(settings$idw_neighbours) > 0,
  all(settings$idw_neighbours >= 1 & settings$idw_neighbours == floor(settings$idw_neighbours)),
  length(settings$matern_scales_m) > 0,
  all(is.finite(settings$matern_scales_m) & settings$matern_scales_m > 0),
  length(settings$matern_smoothness) > 0,
  all(is.finite(settings$matern_smoothness) & settings$matern_smoothness > 0),
  is.finite(settings$matern_nugget_fraction), settings$matern_nugget_fraction > 0,
  settings$matern_nugget_fraction < 1
)
paths <- c(settings$aoi_file, settings$mapping_file)
input_md5 <- unname(tools::md5sum(paths))
inputs <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                              settings$analysis_crs, settings$cell_size_m)
cores <- inputs$cores
core_xy <- sf::st_coordinates(cores)
target_xy <- sf::st_coordinates(inputs$points)
plot_id <- cores[["Plot ID"]]
data_sources <- unique(na.omit(cores$data_source))
map_subtitle <- if (length(data_sources)) paste(data_sources, collapse = ", ") else NULL
columns <- paste0("C Stock 0-", settings$depths_cm, "cm")
if (!all(columns %in% names(cores))) {
  stop("Script 1 stock column(s) missing: ", paste(setdiff(columns, names(cores)), collapse = ", "))
}
for (id in unique(plot_id)) {
  plot_zones <- unique(na.omit(cores$zone_assigned[plot_id == id]))
  if (length(plot_zones) > 1L) stop("Plot ", id, " has cores in more than one zone.")
}

# 1. Area estimates. Average replicate cores within a plot, then weight
# the plots equally. Uncertainty here is the SE of a total under independent
# sampled plots. It is different from uncertainty at an unsampled location.
area_result <- function(stock, selected, area_m2, depth_cm, method, zone) {
  plot_stocks <- tapply(stock[selected], plot_id[selected], mean)
  n <- length(plot_stocks)
  mean_stock <- if (n) mean(plot_stocks) else NA_real_
  stock_sd <- if (n > 1L) sd(plot_stocks) else NA_real_
  se <- stock_sd / sqrt(n)
  data.frame(
    depth_cm = depth_cm, method = method, zone = zone,
    n_cores = sum(selected), n_plots = n, area_m2 = area_m2,
    covered_area_m2 = if (n) area_m2 else 0,
    mean_kg_c_m2 = mean_stock, plot_sd_kg_c_m2 = stock_sd,
    total_mg_c = mean_stock * area_m2 / 1000, total_se_mg_c = se * area_m2 / 1000
  )
}

area_rows <- list()
depth_data <- list()
for (i in seq_along(settings$depths_cm)) {
  depth_cm <- settings$depths_cm[i]
  stock <- cores[[columns[i]]]
  if (!is.numeric(stock) && !all(is.na(stock))) stop("Stock values must be numeric in ", columns[i])
  eligible <- is.finite(stock) & stock >= 0
  if ("Stock coverage (cm)" %in% names(cores)) {
    coverage <- cores[["Stock coverage (cm)"]]
    # Script 1 has already limited the stock to the measured core depth.
    eligible <- eligible & is.finite(coverage) & coverage > 0
  }
  depth_data[[i]] <- list(stock = stock, eligible = eligible)

  pooled <- area_result(stock, eligible, inputs$area_m2, depth_cm, "Pooled", "Entire AOI")
  area_rows[[length(area_rows) + 1L]] <- pooled
  if (nrow(inputs$zones) == 0L) next

  # 1.2. Separate zone means, then add the zone totals. A missing zone or
  # part of the AOI outside all zones leaves the full stratified total unknown.
  zone_rows <- lapply(seq_len(nrow(inputs$zones)), function(j) {
    zone <- inputs$zones$zone[j]
    in_zone <- !is.na(cores$zone_assigned) & cores$zone_assigned == zone
    area_result(stock, eligible & in_zone, inputs$zones$area_m2[j], depth_cm, "By zone", zone)
  })
  zone_results <- do.call(rbind, zone_rows)
  area_tolerance <- max(1e-4, min(1, inputs$area_m2 * 1e-8))
  unzoned_area <- max(0, inputs$area_m2 - sum(zone_results$area_m2))
  if (unzoned_area > area_tolerance) {
    unzoned <- area_result(stock, rep(FALSE, nrow(cores)), unzoned_area, depth_cm,
                           "By zone", "Outside named zones")
    zone_results <- rbind(zone_results, unzoned)
  }
  stratified <- pooled
  stratified$method <- "Stratified"
  stratified$n_cores <- sum(zone_results$n_cores)
  stratified$n_plots <- sum(zone_results$n_plots)
  stratified$covered_area_m2 <- sum(zone_results$covered_area_m2)
  full_coverage <- inputs$area_m2 - stratified$covered_area_m2 <= area_tolerance
  stratified$total_mg_c <- if (full_coverage) sum(zone_results$total_mg_c, na.rm = TRUE) else NA_real_
  stratified$mean_kg_c_m2 <- stratified$total_mg_c * 1000 / inputs$area_m2
  stratified$plot_sd_kg_c_m2 <- NA_real_
  stratified$total_se_mg_c <- if (full_coverage) sqrt(sum(zone_results$total_se_mg_c^2)) else NA_real_
  area_rows[[length(area_rows) + 1L]] <- rbind(zone_results, stratified)
}
area_summary <- do.call(rbind, area_rows)

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
manifest <- file.path(output_folder, "mapping_run.rds")
previous_products <- list.files(
  output_folder, full.names = TRUE,
  pattern = paste0("^(mapping_run\\.rds|area_by_zone\\.gpkg|",
                   "(held_out_predictions|validation_summary|map_totals|",
                   "idw_rankings|matern_rankings|selected_models)\\.csv|",
                   "(area_based|idw|matern)_[0-9.]+cm\\.tif|local_comparison_[0-9.]+cm\\.png)$")
)
if (unlink(previous_products) != 0L) stop("Could not replace previous Script 2 outputs.")
write.csv(area_summary, file.path(output_folder, "area_summary.csv"), row.names = FALSE)
print(area_summary[, c("depth_cm", "method", "zone", "n_plots", "area_m2", "total_mg_c")])

zone_results <- area_summary[area_summary$method == "By zone" &
                               area_summary$zone %in% inputs$zones$zone, , drop = FALSE]
if (nrow(zone_results)) {
  zone_map <- inputs$zones[match(zone_results$zone, inputs$zones$zone), "zone", drop = FALSE]
  zone_map$depth_cm <- zone_results$depth_cm
  zone_map$n_plots <- zone_results$n_plots
  zone_map$mean_kg_c_m2 <- zone_results$mean_kg_c_m2
  zone_map$total_mg_c <- zone_results$total_mg_c
  sf::st_write(zone_map, file.path(output_folder, "area_by_zone.gpkg"),
               layer = "area_by_zone", delete_layer = TRUE, quiet = TRUE)
}

models <- data.frame(depth_cm = numeric(), range_m = numeric(), smoothness = numeric(),
                     nugget_fraction = numeric(), observation_log_sd = numeric(),
                     raster_md5 = character())
validation_rows <- list()
map_totals <- list()
depth_status <- list()
idw_rankings <- list()
matern_rankings <- list()
selected_models <- list()
target_distances <- distance_matrix(core_xy, target_xy)
core_distances <- distance_matrix(core_xy, core_xy)

# Predict every held-out core using only the other plots. Rank by log RMSE,
# preserving candidate order for exact ties. Return the winning predictions
# so validation and final map fitting use the same parameter choice.
rank_candidates <- function(candidates, log_stock, plots, predict_fold) {
  predictions <- lapply(seq_len(nrow(candidates)), function(j) {
    predicted <- numeric(length(log_stock))
    for (plot in unique(plots)) {
      test <- plots == plot
      predicted[test] <- predict_fold(candidates[j, ], !test, test)
    }
    predicted
  })
  errors <- do.call(cbind, predictions) - log_stock
  candidates$rmse_log <- sqrt(colMeans(errors^2))
  candidates$mae_log <- colMeans(abs(errors))
  candidates$bias_log <- colMeans(errors)
  ranked <- order(candidates$rmse_log)
  candidates <- candidates[ranked, , drop = FALSE]
  candidates$rank <- seq_len(nrow(candidates))
  candidates$selected <- candidates$rank == 1L
  list(table = candidates, predicted = predictions[[ranked[1L]]])
}

map_panel <- function(geometry, stock, method) {
  sf::st_sf(mean_kg_c_m2 = stock, method = method, geometry = sf::st_geometry(geometry))
}

# 2 and 3. Work on the log scale for IDW and ordinary Matérn kriging.
# Fit and compare each depth with the same leave-one-plot-out folds.
for (i in seq_along(settings$depths_cm)) {
  depth_cm <- settings$depths_cm[i]
  stock <- depth_data[[i]]$stock
  eligible <- depth_data[[i]]$eligible
  pooled <- area_summary$mean_kg_c_m2[area_summary$depth_cm == depth_cm &
                                      area_summary$method == "Pooled"]
  if (is.finite(pooled)) {
    area_map <- data.frame(mean_kg_c_m2 = rep(pooled, nrow(inputs$cells)))
    export_mapping_raster(inputs, area_map,
                          file.path(output_folder, paste0("area_based_", depth_cm, "cm.tif")))
  }
  positive <- which(eligible & stock > 0)
  folds <- unique(plot_id[positive])
  fit_gp <- length(folds) >= 3L
  depth_status[[i]] <- data.frame(
    depth_cm = depth_cm, mapped_cores = length(positive), mapped_plots = length(folds),
    zero_stocks_area_only = sum(eligible & stock == 0, na.rm = TRUE),
    incomplete_or_invalid_cores = sum(!eligible), idw = length(positive) > 0L,
    matern = fit_gp
  )
  if (!length(positive)) next
  y <- log(stock[positive])
  observed_distances <- core_distances[positive, positive, drop = FALSE]
  to_grid <- target_distances[positive, , drop = FALSE]

  # 2. IDW. Tune power and maximum neighbours on held-out plots.
  # With one plot, use the usual power of 2 and all cores for the map.
  idw_held_out <- rep(NA_real_, length(y))
  idw_choice <- data.frame(power = 2, neighbours = Inf, rmse_log = NA_real_)
  if (length(folds) >= 2L) {
    candidates <- expand.grid(power = settings$idw_powers, neighbours = settings$idw_neighbours,
                              KEEP.OUT.ATTRS = FALSE)
    tuning <- rank_candidates(candidates, y, plot_id[positive], function(candidate, train, test) {
      idw_predict(observed_distances[train, test, drop = FALSE], y[train],
                  candidate$power, candidate$neighbours)
    })
    idw_choice <- tuning$table[1L, ]
    idw_held_out <- tuning$predicted
    idw_rankings[[i]] <- cbind(depth_cm = depth_cm, tuning$table)
  }
  selected_models[[length(selected_models) + 1L]] <- data.frame(
    depth_cm = depth_cm, method = "IDW", power = idw_choice$power,
    neighbours = idw_choice$neighbours, scale_m = NA_real_, smoothness = NA_real_,
    nugget_fraction = NA_real_, rmse_log = idw_choice$rmse_log,
    selection = if (length(folds) >= 2L) "Lowest held-out log RMSE" else "One plot; fixed defaults"
  )

  # IDW has no spatial covariance model. The held-out log residual variance
  # provides an approximate point SD and lognormal mean correction.
  idw_log_mean <- idw_predict(to_grid, y, idw_choice$power, idw_choice$neighbours)
  idw_log_variance <- if (all(is.finite(idw_held_out))) var(y - idw_held_out) else NA_real_
  idw <- log_stock_moments(idw_log_mean, idw_log_variance)
  if (!is.finite(idw_log_variance)) idw$mean_kg_c_m2 <- exp(idw_log_mean)
  maps <- list(IDW = idw)
  total_sd <- c(IDW = NA_real_)
  predictions <- list(IDW = idw_held_out)
  export_mapping_raster(inputs, idw, file.path(output_folder, paste0("idw_", depth_cm, "cm.tif")))

  # 3. Local Matérn GP. Tune distance scale and smoothness, with the same
  # nugget fraction for every candidate. Refit the mean and variance per fold.
  if (fit_gp) {
    candidates <- expand.grid(scale_m = settings$matern_scales_m, smoothness = settings$matern_smoothness,
                              KEEP.OUT.ATTRS = FALSE)
    tuning <- rank_candidates(candidates, y, plot_id[positive], function(candidate, train, test) {
      fit <- matern_predict(observed_distances[train, train, drop = FALSE],
                            observed_distances[train, test, drop = FALSE], y[train],
                            candidate$scale_m, candidate$smoothness,
                            settings$matern_nugget_fraction)
      fit$log_mean
    })
    matern_choice <- tuning$table[1L, ]
    matern_rankings[[i]] <- cbind(depth_cm = depth_cm,
                                  nugget_fraction = settings$matern_nugget_fraction, tuning$table)
    selected_models[[length(selected_models) + 1L]] <- data.frame(
      depth_cm = depth_cm, method = "Matérn GP", power = NA_real_, neighbours = NA_real_,
      scale_m = matern_choice$scale_m, smoothness = matern_choice$smoothness,
      nugget_fraction = settings$matern_nugget_fraction, rmse_log = matern_choice$rmse_log,
      selection = "Lowest held-out log RMSE"
    )
    predictions[["Matérn GP"]] <- tuning$predicted
    gp <- matern_predict(observed_distances, to_grid, y, matern_choice$scale_m,
                         matern_choice$smoothness, settings$matern_nugget_fraction)
    mapped <- log_stock_moments(gp$log_mean, gp$log_variance)
    maps[["Matérn GP"]] <- mapped
    total_sd[["Matérn GP"]] <- lognormal_total_sd(
      mapped$mean_kg_c_m2, inputs$cells$area_m2, target_xy, matern_choice$scale_m,
      matern_choice$smoothness, gp$process_log_variance, gp$conditioned, gp$mean_term
    )
    raster_path <- file.path(output_folder, paste0("matern_", depth_cm, "cm.tif"))
    export_mapping_raster(inputs, cbind(mapped, log_mean = gp$log_mean,
                                        log_variance = gp$log_variance), raster_path)
    models[nrow(models) + 1L, ] <- list(
      depth_cm, matern_choice$scale_m, matern_choice$smoothness,
      settings$matern_nugget_fraction, gp$observation_log_sd, unname(tools::md5sum(raster_path))
    )
  }

  for (method in names(predictions)) {
    predicted <- predictions[[method]]
    validation_rows[[length(validation_rows) + 1L]] <- data.frame(
      depth_cm = depth_cm, method = method, core_id = cores[["Core ID"]][positive],
      plot_id = plot_id[positive], observed_log_stock = y, predicted_log_stock = predicted,
      error_log_stock = predicted - y
    )
  }

  panels <- list(map_panel(inputs$cells, pooled, "Pooled area mean"))
  if (nrow(zone_results)) {
    zoned <- zone_results[zone_results$depth_cm == depth_cm, ]
    zones <- inputs$zones[match(zoned$zone, inputs$zones$zone), ]
    panels[[length(panels) + 1L]] <- map_panel(zones, zoned$mean_kg_c_m2, "Zone area means")
  }
  for (method in names(maps)) {
    prediction <- maps[[method]]$mean_kg_c_m2
    panels[[length(panels) + 1L]] <- map_panel(inputs$cells, prediction, method)
    area <- inputs$cells$area_m2
    total <- sum(prediction * area) / 1000
    map_totals[[length(map_totals) + 1L]] <- data.frame(
      depth_cm = depth_cm, method = method, area_m2 = sum(area),
      mean_kg_c_m2 = total * 1000 / sum(area), total_mg_c = total,
      total_sd_mg_c = unname(total_sd[method]),
      uncertainty = if (method == "IDW") "Joint uncertainty unavailable" else
        if (is.finite(total_sd[method])) "Conditional spatial SD; covariance settings fixed" else
          "Joint calculation skipped above 10,000 grid cells"
    )
  }
  plot_data <- do.call(rbind, panels)
  plot_data$method <- factor(plot_data$method,
                            levels = c("Pooled area mean", "Zone area means", "IDW", "Matérn GP"))
  map <- ggplot() +
    geom_sf(data = plot_data, aes(fill = mean_kg_c_m2), colour = NA) +
    geom_sf(data = inputs$aoi, fill = NA, colour = "grey45", linewidth = 0.15) +
    geom_sf(data = cores[positive, ], shape = 21, fill = "white", size = 1.5, stroke = 0.3) +
    facet_wrap(~ method, ncol = 2) +
    scale_fill_viridis_c(name = "kg C/m²") +
    labs(title = paste0("Local sediment carbon stock, 0–", depth_cm, " cm"),
         subtitle = map_subtitle) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.spacing.x = grid::unit(0.6, "in"))
  ggsave(file.path(output_folder, paste0("local_comparison_", depth_cm, "cm.png")),
         map, width = 9, height = 7, dpi = 180)
}

if (length(validation_rows)) {
  held_out <- do.call(rbind, validation_rows)
  write.csv(held_out, file.path(output_folder, "held_out_predictions.csv"), row.names = FALSE)
  groups <- split(held_out, list(held_out$method, held_out$depth_cm), drop = TRUE)
  validation <- do.call(rbind, lapply(groups, function(group) {
    error <- group$error_log_stock[is.finite(group$error_log_stock)]
    data.frame(depth_cm = group$depth_cm[1], method = group$method[1],
               held_out_cores = length(error), rmse_log = sqrt(mean(error^2)),
               mae_log = mean(abs(error)), bias_log = mean(error))
  }))
  write.csv(validation, file.path(output_folder, "validation_summary.csv"), row.names = FALSE)
}
if (length(map_totals)) {
  write.csv(do.call(rbind, map_totals), file.path(output_folder, "map_totals.csv"), row.names = FALSE)
}
if (length(idw_rankings)) {
  write.csv(do.call(rbind, idw_rankings), file.path(output_folder, "idw_rankings.csv"), row.names = FALSE)
}
if (length(matern_rankings)) {
  write.csv(do.call(rbind, matern_rankings), file.path(output_folder, "matern_rankings.csv"), row.names = FALSE)
}
if (length(selected_models)) {
  selected <- do.call(rbind, selected_models)
  write.csv(selected, file.path(output_folder, "selected_models.csv"), row.names = FALSE)
  print(selected, row.names = FALSE)
}
write.csv(do.call(rbind, depth_status), file.path(output_folder, "depth_status.csv"), row.names = FALSE)
write.csv(models, file.path(output_folder, "matern_settings.csv"), row.names = FALSE)
if (!identical(input_md5, unname(tools::md5sum(paths)))) stop("An input changed during mapping. Run Script 2 again.")
saveRDS(list(settings = settings, models = models, input_md5 = input_md5,
             cell_id = inputs$cells$cell_id, data_sources = data_sources), manifest)
message("Local maps and tables saved to ", output_folder)
