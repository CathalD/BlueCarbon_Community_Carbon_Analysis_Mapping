# Part 2.2. Update a Janousek ecosystem prior with the local sediment cores.
# Run Script_2_local_mapping.R first, then run this script from the workshop folder.

prior_ecosystem <- "Salt marsh"
observation_log_sd <- NULL  # Use the local model's observation SD by depth. Set a number to override.
reference_file <- "data/reference/janousek_core_locations.csv.gz"
local_folder <- "outputs/script2_local"
output_folder <- "outputs/script2_bayesian"

library(sf)
library(terra)
library(ggplot2)
source("helpers/local_mapping_spatial.R")
source("helpers/local_mapping_models.R")

# 1. Load the same AOI, observations, grid and Matérn settings as Script 2.
run_file <- file.path(local_folder, "mapping_run.rds")
if (!file.exists(run_file)) {
  stop("Run Script_2_local_mapping.R to create ", run_file, ".")
}
local_run <- readRDS(run_file)
settings <- local_run$settings
input_md5 <- unname(tools::md5sum(c(settings$aoi_file, settings$mapping_file)))
if (anyNA(input_md5) || !identical(input_md5, local_run$input_md5)) {
  stop("The AOI or mapping data changed. Rerun Script_2_local_mapping.R first.")
}
if (nrow(local_run$models) == 0) {
  stop("The local script did not fit a Matérn model for any depth.")
}
if (!is.null(observation_log_sd) &&
    (length(observation_log_sd) != 1L || !is.numeric(observation_log_sd) ||
     !is.finite(observation_log_sd) || observation_log_sd < 0)) {
  stop("observation_log_sd must be NULL or a finite, nonnegative number.")
}
if (is.null(observation_log_sd) &&
    (is.null(local_run$models$observation_log_sd) ||
     any(!is.finite(local_run$models$observation_log_sd) | local_run$models$observation_log_sd < 0))) {
  stop("The local run is missing observation SDs. Rerun Script_2_local_mapping.R first.")
}

inputs <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                              settings$analysis_crs, settings$cell_size_m)
if (!identical(inputs$cells$cell_id, local_run$cell_id)) {
  stop("The local run has a different set of grid cells. Rerun Script_2_local_mapping.R first.")
}
grid_xy <- st_coordinates(inputs$points)[, c("X", "Y"), drop = FALSE]
cell_area <- inputs$cells$area_m2

# 2. Derive the prior from complete Janousek profiles in the same Pacific
# Northwest region used in Script 1. Source stocks are in Mg C/ha.
habitats <- c(EM = "Salt marsh", SG = "Seagrass", FL = "Tideflat", TS = "Tidal swamp")
if (length(prior_ecosystem) != 1L || !is.character(prior_ecosystem) ||
    is.na(prior_ecosystem)) {
  stop("Choose one Janousek ecosystem for prior_ecosystem.")
}
prior_choice <- tolower(trimws(prior_ecosystem))
ecosystem_index <- match(prior_choice, tolower(unname(habitats)))
if (is.na(ecosystem_index)) ecosystem_index <- match(toupper(prior_choice), names(habitats))
if (is.na(ecosystem_index)) {
  stop("Choose a Janousek ecosystem: ", paste(habitats, collapse = ", "), ".")
}
ecosystem_code <- names(habitats)[ecosystem_index]
ecosystem_name <- unname(habitats[ecosystem_index])
pnw <- c("AK", "BC", "WA", "OR")
janousek <- read.csv(reference_file, colClasses = "character")
janousek <- janousek[janousek$state %in% pnw & janousek$ecosystem %in% ecosystem_code, ]
if (anyDuplicated(janousek$core_id)) stop("Check duplicate Janousek core IDs.")
number <- function(x) suppressWarnings(as.numeric(x))
sampled_cm <- number(sub("^[>~]", "", janousek$core_depth))

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
old_outputs <- list.files(output_folder, full.names = TRUE,
  pattern = "^(posterior_[0-9]+cm\\.tif|comparison_[0-9]+cm\\.png|prior_parameters\\.csv|aoi_totals\\.csv|run_settings\\.csv)$")
if (length(old_outputs) && !all(file.remove(old_outputs))) {
  stop("Could not replace previous Bayesian outputs in ", output_folder, ".")
}
prior_rows <- list()
total_rows <- list()
for (i in seq_len(nrow(local_run$models))) {
  model <- local_run$models[i, ]
  depth <- model$depth_cm
  suffix <- paste0(depth, "cm")
  depth_observation_sd <- if (is.null(observation_log_sd)) model$observation_log_sd else observation_log_sd
  reference_column <- paste0("carbon_stock_", suffix)
  if (!reference_column %in% names(janousek)) {
    stop("Janousek data are missing the stock column for ", suffix, ".")
  }
  reference_stock <- number(janousek[[reference_column]]) / 10
  eligible <- is.finite(sampled_cm) & sampled_cm >= depth &
    is.finite(reference_stock) & reference_stock >= 0
  reference_stock <- reference_stock[eligible]
  n_reference <- length(reference_stock)
  prior_mean <- if (n_reference > 0) mean(reference_stock) else NA_real_
  prior_sd <- if (n_reference > 1) sd(reference_stock) else NA_real_
  estuaries <- janousek$estuary_id[eligible]
  prior_row <- data.frame(
    depth_cm = depth, ecosystem = ecosystem_name, ecosystem_code = ecosystem_code,
    reference_cores = n_reference,
    reference_estuaries = length(unique(estuaries[!is.na(estuaries) & nzchar(estuaries)])),
    prior_mean_kg_c_m2 = prior_mean, prior_sd_kg_c_m2 = prior_sd,
    prior_log_mean = NA_real_, prior_log_sd = NA_real_,
    local_matern_scale_m = model$range_m, matern_smoothness = model$smoothness,
    local_nugget_fraction = model$nugget_fraction,
    observation_log_sd = depth_observation_sd, field_cores = NA_integer_,
    status = "Fitted"
  )
  if (!is.finite(prior_sd) || prior_sd <= 0 || prior_mean <= 0) {
    prior_row$status <- "Insufficient reference stocks to estimate a positive mean and SD"
    prior_rows[[suffix]] <- prior_row
    warning(ecosystem_name, " at ", suffix, ": ", prior_row$status, ".", call. = FALSE)
    next
  }

  # Match the Janousek arithmetic mean and SD with a lognormal distribution.
  # This supplies a constant prior mean and variance on the log scale.
  prior_log_variance <- log1p((prior_sd / prior_mean)^2)
  prior_log_mean <- log(prior_mean) - prior_log_variance / 2
  prior_row$prior_log_mean <- prior_log_mean
  prior_row$prior_log_sd <- sqrt(prior_log_variance)

  # 3. Condition the prior on the local log stocks with the local Matérn
  # correlation. Reuse the local model's observation SD by default to make
  # the influence of the ecosystem prior easier to see.
  stock <- inputs$cores[[paste0("C Stock 0-", suffix)]]
  keep <- is.finite(stock) & stock > 0
  if ("Stock coverage (cm)" %in% names(inputs$cores)) {
    covered <- inputs$cores[["Stock coverage (cm)"]]
    # Local stocks already stop at the measured core depth in Script 1.
    keep <- keep & is.finite(covered) & covered > 0
  }
  if (!any(keep)) stop("No eligible local observations at ", suffix, ".")
  observed <- inputs$cores[keep, ]
  observed_xy <- st_coordinates(observed)[, c("X", "Y"), drop = FALSE]
  prior_row$field_cores <- nrow(observed)
  prior_rows[[suffix]] <- prior_row

  core_correlation <- matern_correlation(distance_matrix(observed_xy, observed_xy),
                                        model$range_m, model$smoothness)
  grid_correlation <- matern_correlation(distance_matrix(observed_xy, grid_xy),
                                        model$range_m, model$smoothness)
  observed_covariance <- prior_log_variance * core_correlation
  diag(observed_covariance) <- diag(observed_covariance) + depth_observation_sd^2
  core_to_grid_covariance <- prior_log_variance * grid_correlation

  # Cholesky solves implement K_go V^-1 (y - prior) and its variance correction.
  covariance_factor <- tryCatch(chol(observed_covariance), error = function(e) {
    stop("Cannot condition the prior at ", suffix,
         ". Check coincident cores and observation_log_sd. ", conditionMessage(e))
  })
  conditioned <- forwardsolve(t(covariance_factor), core_to_grid_covariance)
  residual <- log(stock[keep]) - prior_log_mean
  whitened_residual <- forwardsolve(t(covariance_factor), residual)
  updated_log_mean <- prior_log_mean + as.vector(crossprod(conditioned, whitened_residual))
  updated_log_variance <- prior_log_variance - colSums(conditioned^2)
  if (any(updated_log_variance < -1e-8 * prior_log_variance)) {
    stop("The posterior log variance is negative at ", suffix, ".")
  }
  updated_log_variance <- pmax(updated_log_variance, 0)
  posterior <- log_stock_moments(updated_log_mean, updated_log_variance)

  # 4. Compare predictions at the same cells and use clipped cell areas for totals.
  local_file <- file.path(local_folder, paste0("matern_", suffix, ".tif"))
  if (!file.exists(local_file) ||
      !identical(unname(tools::md5sum(local_file)), model$raster_md5)) {
    stop("The local Matérn raster is missing or has changed: ", local_file,
         ". Rerun Script_2_local_mapping.R first.")
  }
  local_mean_raster <- rast(local_file)[[1]]
  if (!compareGeom(inputs$template, local_mean_raster, stopOnError = FALSE)) {
    stop("The local Matérn raster does not match the prediction grid: ", local_file)
  }
  local_mean <- values(local_mean_raster, mat = FALSE)[inputs$cells$cell_id]
  if (any(!is.finite(local_mean))) {
    stop("The local Matérn raster is missing predictions inside the AOI: ", local_file)
  }
  posterior_mean <- posterior$mean_kg_c_m2
  bands <- data.frame(
    posterior_mean_kg_c_m2 = posterior_mean,
    posterior_sd_kg_c_m2 = posterior$sd_kg_c_m2,
    prior_mean_kg_c_m2 = rep(prior_mean, length(cell_area)),
    prior_sd_kg_c_m2 = rep(prior_sd, length(cell_area)),
    local_matern_mean_kg_c_m2 = local_mean,
    posterior_minus_local_kg_c_m2 = posterior_mean - local_mean,
    posterior_log_mean = updated_log_mean,
    posterior_log_variance = updated_log_variance
  )
  export_mapping_raster(inputs, bands, file.path(output_folder, paste0("posterior_", suffix, ".tif")))

  posterior_total_sd <- lognormal_total_sd(posterior_mean, cell_area, grid_xy,
    model$range_m, model$smoothness, prior_log_variance, conditioned)
  total_rows[[suffix]] <- data.frame(
    depth_cm = depth, method = c("Janousek prior", "Local Matérn", "Bayesian posterior"),
    area_ha = sum(cell_area) / 10000,
    mean_kg_c_m2 = c(prior_mean, weighted.mean(local_mean, cell_area),
                     weighted.mean(posterior_mean, cell_area)),
    total_mg_c = c(prior_mean * sum(cell_area), sum(local_mean * cell_area),
                   sum(posterior_mean * cell_area)) / 1000,
    posterior_total_sd_mg_c = c(NA_real_, NA_real_, posterior_total_sd)
  )

  labels <- c(paste0("Prior (", ecosystem_name, ")"), "Local Matérn", "Bayesian posterior")
  map_data <- rbind(
    st_sf(stock = bands$prior_mean_kg_c_m2, panel = labels[1], geometry = st_geometry(inputs$cells)),
    st_sf(stock = local_mean, panel = labels[2], geometry = st_geometry(inputs$cells)),
    st_sf(stock = posterior_mean, panel = labels[3], geometry = st_geometry(inputs$cells))
  )
  map_data$panel <- factor(map_data$panel, levels = labels)
  observation_source <- unique(na.omit(observed$data_source))
  core_caption <- if (length(observation_source) == 1L) {
    paste("Points show", tolower(observation_source))
  } else "Points show core locations"
  comparison <- ggplot(map_data) +
    geom_sf(aes(fill = stock), colour = NA) +
    geom_sf(data = observed, shape = 21, fill = "white", colour = "#1b1b1b",
            size = 1.4, stroke = 0.35) +
    facet_wrap(vars(panel), nrow = 1) +
    scale_fill_viridis_c(name = "kg C/m²") +
    labs(title = paste0("Sediment organic carbon stock, 0–", depth, " cm"),
         caption = core_caption) +
    theme_void() +
    theme(legend.position = "bottom", panel.spacing = grid::unit(1, "lines"))
  ggsave(file.path(output_folder, paste0("comparison_", suffix, ".png")),
         comparison, width = 10, height = 4.2, dpi = 180, bg = "white")
}

prior_parameters <- do.call(rbind, prior_rows)
write.csv(prior_parameters, file.path(output_folder, "prior_parameters.csv"), row.names = FALSE)
if (!length(total_rows)) stop("No Bayesian maps were fitted. See prior_parameters.csv.")
aoi_totals <- do.call(rbind, total_rows)
write.csv(aoi_totals, file.path(output_folder, "aoi_totals.csv"), row.names = FALSE)

run_settings <- data.frame(
  setting = c("prior_ecosystem", "reference_file", "reference_md5", "reference_region",
              "reference_doi", "observation_log_sd_setting", "local_run", "aoi_file",
              "aoi_md5", "mapping_file", "mapping_md5", "cell_size_m",
              "posterior_total_sd"),
  value = c(ecosystem_name, reference_file, unname(tools::md5sum(reference_file)),
            paste(pnw, collapse = ", "), "10.1029/2024GB008239",
            if (is.null(observation_log_sd)) "From local Matérn model by depth" else observation_log_sd,
            run_file, settings$aoi_file, input_md5[1], settings$mapping_file,
            input_md5[2], settings$cell_size_m,
            "Conditional on fixed prior, Matérn and observation settings; computed with covariance among AOI cells when within the grid size limit.")
)
write.csv(run_settings, file.path(output_folder, "run_settings.csv"), row.names = FALSE)
print(aoi_totals, row.names = FALSE)
