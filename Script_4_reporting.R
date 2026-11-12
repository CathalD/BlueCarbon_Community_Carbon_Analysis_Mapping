# Part 4 | Summarize the cores and the maps.
# Run from the workshop folder after Script 1 and the local mapping script.
# Run the Bayesian and transfer scripts to include their maps and summaries.
# Add --render-report when running from Rscript to build Blue_Carbon_Report.html.
# Packages: sf, terra, dplyr, tidyr, readr, ggplot2, leaflet, htmlwidgets, base64enc.

field_folder <- "outputs/script1"
local_folder <- "outputs/script2_local"
bayesian_folder <- "outputs/script2_bayesian"
transfer_folder <- "outputs/script3_covariate_transfer"
transfer_inputs_folder <- "data/covariate_transfer"
output_folder <- "outputs/script4_reporting"
reference_file <- "data/reference/janousek_core_locations.csv.gz"
estuary_file <- "data/reference/janousek_estuaries.csv"
bc_estuary_habitat <- "Salt marsh"

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
source("helpers/local_mapping_spatial.R")
source("helpers/view_outputs.R")

read_output <- function(folder, filename, ...) {
  path <- file.path(folder, filename)
  if (!file.exists(path)) stop("Missing ", path, ". Run the script that produces it first.")
  read_csv(path, show_col_types = FALSE, ...)
}

run <- readRDS(file.path(local_folder, "mapping_run.rds"))
settings <- run$settings
input_md5 <- unname(tools::md5sum(c(settings$aoi_file, settings$mapping_file)))
if (anyNA(input_md5) || !identical(input_md5, run$input_md5)) {
  stop("The boundary or mapping data have changed. Rerun the local mapping script first.")
}
spatial <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                               settings$analysis_crs, settings$cell_size_m)
stopifnot(identical(spatial$cells$cell_id, run$cell_id))
cores <- read_output(field_folder, "core_summary.csv", col_types = cols(
  `Core ID` = col_character(), `Plot ID` = col_character()
))
field_summary <- read_output(field_folder, "field_stock_summary.csv")
core_context <- read_output(field_folder, "janousek_comparison.csv")
sections <- read_output(field_folder, "sediment_sections.csv")
depth_status <- read_output(local_folder, "depth_status.csv")
area <- read_output(local_folder, "area_summary.csv")
read_optional <- function(folder, filename) {
  if (file.exists(file.path(folder, filename))) read_output(folder, filename)
}
model_settings <- read_optional(local_folder, "selected_models.csv")
mapped <- read_optional(local_folder, "map_totals.csv")
prior <- read_optional(bayesian_folder, "prior_parameters.csv")
bayesian_run <- read_optional(bayesian_folder, "run_settings.csv")
bayesian <- NULL
if (!is.null(bayesian_run)) {
  values <- setNames(bayesian_run$value, bayesian_run$setting)
  if (!identical(unname(values[c("aoi_md5", "mapping_md5")]), input_md5)) {
    stop("The Bayesian maps use different inputs. Rerun the Bayesian script.")
  }
  bayesian <- read_output(bayesian_folder, "aoi_totals.csv")
}

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
save_table <- function(data, name) write_csv(data, file.path(output_folder, paste0(name, ".csv")))
save_figure <- function(plot, name, width = 10, height = 5) {
  ggsave(file.path(output_folder, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 180, bg = "white")
}
depth_label <- function(x) paste0("0–", x, " cm")
source_label <- paste(unique(na.omit(cores$data_source)), collapse = ", ")
report_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# 1. Core coverage, measured stocks and data checks.
save_table(tibble(
  aoi_file = settings$aoi_file, area_ha = spatial$area_m2 / 10000,
  zones = nrow(spatial$zones), logged_cores = nrow(cores),
  plots = n_distinct(cores$`Plot ID`, na.rm = TRUE), mapped_cores = nrow(spatial$cores),
  cores_to_review = sum(is.na(cores$core_check) | cores$core_check != "OK"),
  sections = nrow(sections), grid_spacing_m = settings$cell_size_m, data_source = source_label,
  reference_file = reference_file, bc_estuary_habitat = bc_estuary_habitat
), "site_summary")
save_table(cores, "core_stocks")
save_table(field_summary, "field_stock_summary")
save_table(core_context, "core_reference_comparison")
save_table(count(sections, carbon_method, section_check, lab_review, name = "sections"),
           "section_checks")
save_table(depth_status, "mapping_coverage")

stocks <- cores |>
  select(`Core ID`, zone, matches("^C Stock 0-[0-9]+cm$")) |>
  pivot_longer(matches("^C Stock 0-[0-9]+cm$"), names_to = "depth_cm",
               names_pattern = "C Stock 0-(\\d+)cm", names_transform = list(depth_cm = as.integer),
               values_to = "stock_kg_m2") |>
  filter(is.finite(stock_kg_m2)) |>
  mutate(depth = factor(depth_label(depth_cm), levels = depth_label(sort(unique(depth_cm)))))
if (!nrow(stocks)) stop("No complete fixed-depth stocks found in the Script 1 core summary.")
save_figure(
  ggplot(stocks, aes(zone, stock_kg_m2, colour = zone)) +
    geom_boxplot(outlier.shape = NA, width = 0.52, alpha = 0.12) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 42), size = 2) +
    facet_wrap(vars(depth), scales = "free_y") +
    labs(title = "Carbon stock in the study cores", subtitle = source_label,
         x = NULL, y = "Carbon stock (kg C/m²)") +
    report_theme + guides(colour = "none"),
  "core_stocks"
)
profile_figures <- c("sediment_profiles.png", "cumulative_carbon_stock.png", "janousek_profiles.png")
available_figures <- file.path(field_folder, profile_figures)
available_figures <- available_figures[file.exists(available_figures)]
invisible(file.copy(available_figures, output_folder, overwrite = TRUE))

# 2. Compare stocks over the same measured depth. The B.C. estuary figure
# follows the NRE comparison and uses names from the published estuary lookup.
habitats <- c(EM = "Salt marsh", SG = "Seagrass", FL = "Tideflat", TS = "Tidal swamp")
if (!bc_estuary_habitat %in% habitats) {
  stop("Choose a B.C. estuary habitat from ", paste(habitats, collapse = ", "), ".")
}
depth_columns <- paste0("carbon_stock_", sort(unique(stocks$depth_cm)), "cm")
locations <- read_csv(reference_file, col_types = cols(.default = col_character()))
estuary_names <- read_csv(estuary_file, col_types = cols(.default = col_character())) |>
  mutate(estuary_region = state) |>
  separate_longer_delim(state, delim = "/")
if (!all(depth_columns %in% names(locations))) stop("Janousek stocks are unavailable for a reported depth.")
number <- function(x) suppressWarnings(as.numeric(x))
reference <- locations |>
  filter(state %in% c("AK", "BC", "WA", "OR"), ecosystem %in% names(habitats)) |>
  mutate(sampled_cm = number(sub("^[>~]", "", core_depth)), habitat = unname(habitats[ecosystem])) |>
  pivot_longer(all_of(depth_columns), names_to = "depth_cm",
               names_pattern = "carbon_stock_(\\d+)cm", names_transform = list(depth_cm = as.integer),
               values_transform = list(value = number)) |>
  filter(sampled_cm >= depth_cm, is.finite(value), value >= 0) |>
  mutate(stock_kg_m2 = value / 10, depth = factor(depth_label(depth_cm), levels = levels(stocks$depth)))

stock_distribution <- function(data, groups) {
  data |>
    group_by(across(all_of(groups))) |>
    summarise(cores = n_distinct(core_id), estuaries = n_distinct(estuary_id, na.rm = TRUE),
              mean_kg_m2 = mean(stock_kg_m2), sd_kg_m2 = sd(stock_kg_m2),
              median_kg_m2 = median(stock_kg_m2), q25_kg_m2 = quantile(stock_kg_m2, .25),
              q75_kg_m2 = quantile(stock_kg_m2, .75), .groups = "drop")
}
ecosystems <- bind_rows(
  mutate(reference, region = "Pacific Northwest"),
  reference |> filter(state == "BC") |> mutate(region = "British Columbia")
)
bc <- reference |> filter(state == "BC", !is.na(estuary_id), nzchar(estuary_id))
save_table(stock_distribution(ecosystems, c("region", "habitat", "depth_cm")), "janousek_ecosystems")
save_table(stock_distribution(bc, c("habitat", "estuary_id", "depth_cm")), "janousek_bc_estuaries")
named_estuaries <- reference |>
  filter(!is.na(estuary_id), nzchar(estuary_id)) |>
  left_join(estuary_names, by = c("estuary_id", "state"), relationship = "many-to-one") |>
  mutate(estuary_name = coalesce(estuary_name, estuary_id), estuary_region = coalesce(estuary_region, state))
save_table(stock_distribution(named_estuaries,
                             c("habitat", "estuary_id", "estuary_region", "estuary_name", "depth_cm")),
           "janousek_estuary_comparison")

local_medians <- stocks |> group_by(depth) |> summarise(stock_kg_m2 = median(stock_kg_m2))
context_figure <- function(data, group, title, subtitle, row_facet = FALSE) {
  graph <- ggplot(data, aes(x = stock_kg_m2, y = .data[[group]])) +
    geom_boxplot(fill = "#edf1f1", colour = "#52626a", outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitter(width = 0, height = 0.12, seed = 42),
               colour = "#52626a", alpha = 0.35, size = 1.1) +
    geom_vline(data = local_medians, aes(xintercept = stock_kg_m2), colour = "#ab3b26") +
    labs(title = title, subtitle = subtitle, x = "Carbon stock (kg C/m²)", y = NULL,
         caption = paste("Each point is a Janousek core. The red line is the study median.", source_label)) +
    report_theme
  if (row_facet) graph + facet_grid(region ~ depth, scales = "free_x") else
    graph + facet_wrap(vars(depth), scales = "free_x")
}
save_figure(context_figure(ecosystems, "habitat", "Carbon stocks by ecosystem",
                          "Janousek cores sampled to the reported depth", TRUE),
            "janousek_ecosystem_context", height = 6)
bc_selected <- named_estuaries |> filter(state == "BC", habitat == bc_estuary_habitat)
if (nrow(bc_selected)) {
  save_figure(context_figure(bc_selected, "estuary_name", "Carbon stocks across B.C. estuaries",
                            paste(bc_estuary_habitat, "reference cores grouped by estuary")),
              "janousek_bc_estuary_context", height = max(4, 2 + 0.45 * n_distinct(bc_selected$estuary_id)))
}

# 3. Report total stock and distinguish sampling SE from model SD.
by_zone <- area |> filter(method == "By zone") |> mutate(area_ha = area_m2 / 10000)
save_table(by_zone, "area_by_zone")
totals <- area |>
  filter(method %in% c("Pooled", "Stratified")) |>
  transmute(depth_cm, method = recode(method, Pooled = "Pooled area mean", Stratified = "Zone-weighted area mean"),
            area_ha = area_m2 / 10000, coverage_pct = 100 * covered_area_m2 / area_m2,
            mean_kg_c_m2, total_mg_c, uncertainty_mg_c = total_se_mg_c,
            uncertainty_type = if_else(is.finite(total_se_mg_c), "Sampling SE", NA_character_))
if (!is.null(mapped)) {
  totals <- bind_rows(totals, mapped |>
    transmute(depth_cm, method, area_ha = area_m2 / 10000, coverage_pct = 100,
              mean_kg_c_m2, total_mg_c, uncertainty_mg_c = total_sd_mg_c,
              uncertainty_type = if_else(is.finite(total_sd_mg_c), "Model SD", NA_character_)))
}
if (!is.null(bayesian)) {
  totals <- bind_rows(totals, bayesian |>
    filter(method %in% c("Janousek prior", "Bayesian posterior")) |>
    transmute(depth_cm, method, area_ha, coverage_pct = 100, mean_kg_c_m2, total_mg_c,
              uncertainty_mg_c = posterior_total_sd_mg_c,
              uncertainty_type = if_else(is.finite(posterior_total_sd_mg_c), "Posterior SD", NA_character_)))
  save_table(prior, "bayesian_prior_parameters")
}
save_table(arrange(totals, depth_cm), "mapping_totals")
if (!is.null(model_settings)) save_table(model_settings, "mapping_settings")
totals <- totals |> mutate(depth = factor(depth_label(depth_cm),
                                          levels = depth_label(sort(unique(depth_cm)))))
save_figure(
  ggplot(filter(totals, is.finite(total_mg_c)), aes(total_mg_c / 1000, method)) +
    geom_errorbar(data = filter(totals, is.finite(total_mg_c) & is.finite(uncertainty_mg_c)),
                  aes(xmin = (total_mg_c - uncertainty_mg_c) / 1000,
                      xmax = (total_mg_c + uncertainty_mg_c) / 1000, colour = uncertainty_type),
                  orientation = "y", width = 0.22) +
    geom_point(size = 2.5) + facet_wrap(vars(depth), scales = "free_x") +
    scale_y_discrete(labels = c("Pooled area mean" = "Area average",
                                "Zone-weighted area mean" = "Area average by zone",
                                "IDW" = "Distance-weighted average",
                                "Matérn GP" = "Local Matérn",
                                "Janousek prior" = "Janousek prior",
                                "Bayesian posterior" = "Bayesian estimate")) +
    scale_colour_discrete(labels = c("Sampling SE" = "Sampling standard error",
                                     "Model SD" = "Model standard deviation",
                                     "Posterior SD" = "Posterior standard deviation")) +
    labs(title = "Total carbon stock", subtitle = source_label,
         x = "Carbon stock (thousand tonnes)", y = "Method", colour = "Error bars") +
    report_theme,
  "mapping_totals"
)

# 4. Summarize transfer learning with the embeddings kernel.
transfer <- read_optional(transfer_folder, "transfer_area_summary.csv")
if (!is.null(transfer)) {
  transfer_predictions <- read_output(transfer_folder, "transfer_predictions.csv")
  transfer_stats <- transfer_predictions |>
    pivot_longer(cols = matches("^(prior|posterior)_(mean|sd)_kg_c_m2$"),
                 names_to = c("surface", ".value"),
                 names_pattern = "(prior|posterior)_(mean|sd)_kg_c_m2") |>
    filter(is.finite(mean), is.finite(sd)) |>
    mutate(cv = .stock_cv(mean, sd)) |>
    group_by(depth_cm, surface) |>
    summarise(prediction_cells_with_sd = n(),
              area_weighted_point_sd_kg_c_m2 = weighted.mean(sd, area_m2),
              min_point_sd_kg_c_m2 = min(sd), max_point_sd_kg_c_m2 = max(sd),
              cells_with_cv = sum(is.finite(cv)),
              area_weighted_cv = if (any(is.finite(cv))) weighted.mean(cv, area_m2, na.rm = TRUE) else NA_real_,
              min_cv = if (any(is.finite(cv))) min(cv, na.rm = TRUE) else NA_real_,
              max_cv = if (any(is.finite(cv))) max(cv, na.rm = TRUE) else NA_real_,
              median_effective_donors = median(effective_donors), .groups = "drop")
  transfer_summary <- transfer |>
    left_join(transfer_stats, by = c("depth_cm", "surface")) |>
    mutate(prediction_cells_with_sd = coalesce(prediction_cells_with_sd, 0L))
  save_table(transfer_summary, "transfer_summary")
  save_table(read_output(transfer_folder, "transfer_local_priors.csv", col_types = cols(
    `Core ID` = col_character(), `Plot ID` = col_character()
  )), "transfer_local_priors")

  transfer_view <- .transfer_view_inputs(transfer, spatial, settings, transfer_folder, transfer_inputs_folder,
                                         read_grid = read_mapping_inputs)
  transfer_spatial <- transfer_view$grid
  embedding_map <- transfer_view$embedding_map
}

# 5. Save a Leaflet page and a Quarto widget per depth for each set of models.
# Core points and the AOI outline stay available as separate overlays.
previous_maps <- list.files(output_folder, pattern = "^map_(transfer_)?(cv_)?[0-9.]+cm\\.(html|rds|tif)$",
                            full.names = TRUE)
unlink(c(previous_maps, file.path(output_folder, "map_libraries")), recursive = TRUE)
pooled <- filter(area, method == "Pooled")
write_maps <- function(grid, depths, transfer_only = FALSE) {
  map_index <- list()
  surface_stats <- list()
  for (depth_cm in sort(unique(depths))) {
    result <- .build_output_map(grid, depth_cm, source_label, local_folder, bayesian_folder,
                                transfer_folder, by_zone, prior, bayesian,
                                if (transfer_only) embedding_map else NULL, transfer_only, pooled = pooled)
    if (is.null(result)) next
    suffix <- paste0(depth_cm, "cm")
    prefix <- if (transfer_only) "map_transfer_" else "map_"
    map_file <- paste0(prefix, suffix, ".html")
    widget_file <- paste0(prefix, suffix, ".rds")
    htmlwidgets::saveWidget(result$widget, file.path(output_folder, map_file),
                            selfcontained = FALSE, libdir = "map_libraries", title = depth_label(depth_cm))
    saveRDS(result$widget, file.path(output_folder, widget_file))
    if (length(result$cv_rasters)) {
      cv <- do.call(c, unname(result$cv_rasters))
      names(cv) <- names(result$cv_rasters)
      terra::writeRaster(cv, file.path(output_folder, paste0(prefix, "cv_", suffix, ".tif")),
                         overwrite = TRUE)
    }
    map_index[[suffix]] <- tibble(depth_cm = depth_cm, map_file = map_file, widget_file = widget_file,
                                  default_surface = result$default_surface, surfaces = length(result$surfaces))
    surface_stats[[suffix]] <- result$stats
  }
  list(index = bind_rows(map_index), stats = bind_rows(surface_stats))
}

local_maps <- write_maps(spatial, settings$depths_cm)
save_table(local_maps$index, "interactive_maps")
save_table(local_maps$stats, "map_surface_summary")
if (!is.null(transfer)) {
  transfer_maps <- write_maps(transfer_spatial, transfer$depth_cm, transfer_only = TRUE)
  save_table(transfer_maps$index, "transfer_interactive_maps")
  save_table(transfer_maps$stats, "transfer_map_surface_summary")
}
message("Reporting tables, figures and interactive maps saved to ", output_folder)

if ("--render-report" %in% commandArgs(trailingOnly = TRUE)) {
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) stop("Quarto is required to render Blue_Carbon_Report.qmd.")
  if (system2(quarto, c("render", "Blue_Carbon_Report.qmd")) != 0L) {
    stop("Quarto could not render Blue_Carbon_Report.qmd.")
  }
}
