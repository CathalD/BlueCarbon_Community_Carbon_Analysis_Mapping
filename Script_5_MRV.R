# Script 5 | Baseline MRV and repeat-sampling plan
# Run from the workshop folder after Scripts 1 and 2, or via Run_and_Visualize.R.
# Planning values are site-specific. Leave unknown inputs as NA.
local({
  plan <- list(
    reporting_depth_cm = 30,
    meaningful_change_mg_c_ha = NA_real_,
    expected_change_mg_c_ha_yr = NA_real_, # Magnitude expected for the reporting estimand
    paired_change_sd_mg_c_ha = NA_real_,    # SD of paired plot differences from comparable repeat cores
    planned_plot_pairs = NA_integer_,
    baseline_date = NA_character_,          # Confirmed baseline campaign end date, YYYY-MM-DD
    confidence = 0.95,
    power = 0.80,
    focus_areas_per_zone = 2L
  )
  workbook <- "data/field/blue_carbon_field_cores.xlsm"
  field <- "outputs/script1"
  local_models <- "outputs/script2_local"
  out <- "outputs/script5_mrv"

  library(dplyr)
  library(ggplot2)
  source("helpers/local_mapping_spatial.R", local = TRUE)

  optional_positive <- c(plan$meaningful_change_mg_c_ha, plan$expected_change_mg_c_ha_yr,
                         plan$paired_change_sd_mg_c_ha)
  stopifnot(all(is.na(optional_positive) | (is.finite(optional_positive) & optional_positive > 0)),
            is.na(plan$planned_plot_pairs) || (is.finite(plan$planned_plot_pairs) &&
              plan$planned_plot_pairs >= 2 && plan$planned_plot_pairs %% 1 == 0),
            plan$confidence > 0.5, plan$confidence < 1, plan$power > 0.5, plan$power < 1,
            is.finite(plan$focus_areas_per_zone), plan$focus_areas_per_zone >= 1,
            plan$focus_areas_per_zone %% 1 == 0)
  baseline_date <- as.Date(plan$baseline_date, format = "%Y-%m-%d")
  if (!is.na(plan$baseline_date) && is.na(baseline_date)) {
    stop("Enter baseline_date as YYYY-MM-DD or leave it as NA.")
  }

  manifest <- file.path(local_models, "mapping_run.rds")
  if (!file.exists(manifest)) stop("Run Script_2_local_mapping.R to completion before Script 5.")
  run <- readRDS(manifest)
  settings <- run$settings
  depths <- settings$depths_cm
  if (!plan$reporting_depth_cm %in% depths) stop("Set reporting_depth_cm to a depth in Script 2.")
  if (!identical(run$input_md5, unname(tools::md5sum(c(settings$aoi_file, settings$mapping_file))))) {
    stop("The boundary or mapped core data changed. Rerun Script 2.")
  }
  core_summary <- file.path(field, "core_summary.csv")
  if (!file.exists(workbook) || !file.exists(core_summary) ||
      file.mtime(workbook) > min(file.mtime(core_summary), file.mtime(settings$mapping_file))) {
    stop("Run Scripts 1 and 2 with the current field workbook before Script 5.")
  }

  read_csv <- function(path, ...) readr::read_csv(path, show_col_types = FALSE, ...)
  id_types <- readr::cols(`Core ID` = readr::col_character(), `Plot ID` = readr::col_character())
  records <- read_csv(core_summary, col_types = id_types)
  raw <- readxl::read_excel(workbook, sheet = "2. Core Log", skip = 3, col_types = "list") |>
    filter(!vapply(`Core ID`, function(x) all(is.na(x)), logical(1)))
  read_id <- function(x) vapply(x, as.character, character(1))
  raw_ids <- read_id(raw$`Core ID`)
  raw_plot <- read_id(raw$`Plot ID`)
  if (anyDuplicated(raw_ids) || !setequal(raw_ids, records$`Core ID`)) {
    stop("Core IDs in the raw workbook and Script 1 summary differ. Rerun Scripts 1 and 2.")
  }
  raw_order <- match(records$`Core ID`, raw_ids)
  if (!identical(raw_plot[raw_order], records$`Plot ID`)) {
    stop("Plot IDs in the raw workbook and Script 1 summary differ. Rerun Scripts 1 and 2.")
  }
  date_text <- function(x) {
    if (all(is.na(x))) return(NA_character_)
    if (inherits(x, "POSIXt")) return(format(x, "%Y-%m-%d", tz = "UTC"))
    if (is.numeric(x)) return(as.character(as.Date(x, origin = "1899-12-30")))
    if (grepl("^\\d{4}-\\d{2}-\\d{2}$", x)) return(as.character(as.Date(x, "%Y-%m-%d")))
    NA_character_
  }
  records$coring_date <- vapply(raw$`Coring date`[raw_order], date_text, character(1))

  inputs <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                                settings$analysis_crs, settings$cell_size_m)
  if (!identical(inputs$cells$cell_id, run$cell_id)) stop("Grid differs from Script 2. Rerun Script 2.")
  cores <- inputs$cores
  core_ids <- cores$`Core ID`
  grid <- inputs$points
  cell_area <- inputs$cells$area_m2
  xy <- terra::xyFromCell(inputs$template, grid$cell_id)
  core_index <- match(core_ids, records$`Core ID`)
  if (anyNA(core_index) || !identical(cores$`Plot ID`, records$`Plot ID`[core_index])) {
    stop("Mapped core records differ from Script 1. Rerun Scripts 1 and 2.")
  }

  zones <- inputs$zones[c("zone", "geometry")]
  if (!nrow(zones)) {
    zones <- mutate(inputs$aoi, zone = "Study area")
    grid$zone_index <- rep(1L, nrow(grid))
    core_zone <- rep(1L, nrow(cores))
  } else {
    hits <- sf::st_intersects(grid, zones)
    grid$zone_index <- vapply(hits, function(hit) c(hit, NA_integer_)[1L], integer(1))
    core_zone <- match(cores$zone_assigned, zones$zone)
    remainder <- sf::st_difference(sf::st_geometry(inputs$aoi), sf::st_union(sf::st_geometry(zones)))
    if (sum(as.numeric(sf::st_area(remainder))) > max(1e-4, min(1, inputs$area_m2 * 1e-8))) {
      zones <- rbind(zones, sf::st_sf(zone = "Outside named zones", geometry = remainder))
      grid$zone_index[is.na(grid$zone_index)] <- nrow(zones)
      core_zone[is.na(core_zone)] <- nrow(zones)
    }
  }

  selected_file <- file.path(local_models, "selected_models.csv")
  validation_file <- file.path(local_models, "held_out_predictions.csv")
  raster_files <- c("IDW" = "idw", "Matérn GP" = "matern")
  empty_selection <- tibble(depth_cm = numeric(), method = character(), rmse_log = numeric())
  selected <- if (file.exists(selected_file)) read_csv(selected_file) else empty_selection
  held <- if (file.exists(validation_file)) {
    read_csv(validation_file, col_types = readr::cols(core_id = readr::col_character(),
                                                      plot_id = readr::col_character()))
  } else {
    tibble(depth_cm = numeric(), method = character(), core_id = character(), plot_id = character(),
           observed_log_stock = numeric(), predicted_log_stock = numeric())
  }
  model_totals_file <- file.path(local_models, "map_totals.csv")
  model_totals <- if (file.exists(model_totals_file)) read_csv(model_totals_file) else tibble(
    depth_cm = numeric(), method = character(), total_sd_mg_c = numeric(), uncertainty = character()
  )
  stopifnot(all(selected$method %in% names(raster_files)),
            !anyDuplicated(selected[c("depth_cm", "method")]),
            !anyDuplicated(held[c("depth_cm", "method", "core_id")]))
  nearest_distance <- function(points, locations) {
    if (!nrow(locations)) return(rep(NA_real_, nrow(points)))
    nearest <- sf::st_nearest_feature(points, locations)
    as.numeric(sf::st_distance(points, locations[nearest, ], by_element = TRUE))
  }

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  depth_results <- lapply(depths, function(depth) {
    column <- paste0("C Stock 0-", depth, "cm")
    stock <- cores[[column]]
    logged_stock <- records[[column]]
    if (!isTRUE(all.equal(stock, logged_stock[core_index], check.attributes = FALSE)) ||
        !isTRUE(all.equal(cores$`Stock coverage (cm)`, records$`Stock coverage (cm)`[core_index]))) {
      stop("The mapped stock values differ from Script 1. Rerun Scripts 1 and 2.")
    }
    available <- is.finite(stock) & stock >= 0
    full <- available & is.finite(cores$`Stock coverage (cm)`) &
      cores$`Stock coverage (cm)` >= depth - 1e-8
    review <- is.na(cores$core_check) | cores$core_check != "OK"
    support <- cores[full & !review, ]
    summary_rows <- selected[selected$depth_cm == depth, ]
    scores <- held[held$depth_cm == depth & is.finite(held$predicted_log_stock) &
                     is.finite(held$observed_log_stock), ]
    validation_index <- match(scores$core_id, core_ids)
    if (anyNA(validation_index) || any(scores$plot_id != cores$`Plot ID`[validation_index]) ||
        !isTRUE(all.equal(scores$observed_log_stock, log(stock[validation_index])))) {
      stop("Held-out records do not match the mapped stocks. Rerun Script 2.")
    }
    comparable <- intersect(summary_rows$method, unique(scores$method))
    common_ids <- if (length(comparable)) {
      Reduce(intersect, lapply(comparable, function(m) scores$core_id[scores$method == m]))
    } else character()
    ranking <- scores |>
      filter(method %in% comparable, core_id %in% common_ids) |>
      mutate(error_log = predicted_log_stock - observed_log_stock,
             error_mg_c_ha = 10 * (exp(predicted_log_stock) - exp(observed_log_stock))) |>
      group_by(method) |>
      summarise(held_out_cores = n(), held_out_plots = n_distinct(plot_id),
                rmse_log = sqrt(mean(error_log^2)), rmse_mg_c_ha = sqrt(mean(error_mg_c_ha^2)),
                mae_mg_c_ha = mean(abs(error_mg_c_ha)), bias_mg_c_ha = mean(error_mg_c_ha),
                .groups = "drop") |>
      arrange(rmse_log, method) |>
      mutate(depth_cm = depth, selected = row_number() == 1L, .before = 1)
    if (length(comparable) && !nrow(ranking)) stop("No common held-out cores at ", depth, " cm.")
    method <- if (nrow(ranking)) ranking$method[1L] else if (nrow(summary_rows) == 1L) {
      summary_rows$method[1L]
    } else NA_character_
    if (nrow(summary_rows) > 1L && !nrow(ranking)) {
      stop("Multiple fitted models have no comparable held-out predictions at ", depth, " cm.")
    }

    mean_stock <- point_sd <- rep(NA_real_, nrow(grid))
    if (!is.na(method)) {
      path <- file.path(local_models, paste0(raster_files[[method]], "_", depth, "cm.tif"))
      raster <- terra::rast(path)
      if (!isTRUE(terra::compareGeom(raster, inputs$template, stopOnError = FALSE)) ||
          !all(c("mean_kg_c_m2", "sd_kg_c_m2") %in% names(raster))) {
        stop("Unexpected model raster grid or bands at ", depth, " cm. Rerun Script 2.")
      }
      mapped <- terra::values(raster, mat = TRUE)[grid$cell_id, , drop = FALSE]
      mean_stock <- mapped[, "mean_kg_c_m2"] * 10
      point_sd <- mapped[, "sd_kg_c_m2"] * 10
    }
    predicted <- is.finite(mean_stock)
    mapped_area <- sum(cell_area[predicted])
    total <- if (any(predicted)) sum(mean_stock[predicted] * cell_area[predicted] / 10000) else NA_real_
    complete_area <- all(predicted)
    total_sd <- model_totals[model_totals$depth_cm == depth & model_totals$method %in% method, ]
    joint_sd <- if (complete_area && nrow(total_sd)) total_sd$total_sd_mg_c[1L] else NA_real_
    joint_status <- if (nrow(total_sd)) total_sd$uncertainty[1L] else "Unavailable"
    nearest <- nearest_distance(grid, support)
    mean_stock[!predicted] <- NA_real_
    point_sd[!is.finite(point_sd)] <- NA_real_

    # Each focus area is a grid point around which to randomize new plots.
    # Spread multiple focus areas by updating distance after each selection.
    zone_rows <- vector("list", nrow(zones))
    focus_rows <- vector("list", nrow(zones))
    for (j in seq_len(nrow(zones))) {
      z <- which(grid$zone_index %in% j)
      in_zone <- core_zone %in% j
      usable <- in_zone & full & !review
      plot_stocks <- tapply(10 * stock[usable], cores$`Plot ID`[usable], mean)
      n_plots <- length(plot_stocks)
      spatial_gap <- if (length(z) && any(is.finite(nearest[z]))) max(nearest[z], na.rm = TRUE) else NA_real_
      zone_rows[[j]] <- tibble(
        depth_cm = depth, zone = zones$zone[j], area_ha = as.numeric(sf::st_area(zones[j, ])) / 10000,
        full_depth_qc_plots = n_plots, cores_logged_on_map = sum(in_zone),
        short_or_invalid_cores = sum(in_zone & !full), flagged_cores = sum(in_zone & review),
        max_gap_m = spatial_gap,
        observed_plot_mean_mg_c_ha = if (n_plots) mean(plot_stocks) else NA_real_,
        observed_plot_sd_mg_c_ha = if (n_plots > 1L) sd(plot_stocks) else NA_real_,
        median_point_sd_mg_c_ha = if (any(is.finite(point_sd[z]))) {
          median(point_sd[z], na.rm = TRUE)
        } else NA_real_,
        action = case_when(
          n_plots == 0 ~ "Establish randomly selected full-depth plots after checking habitat and access.",
          n_plots == 1 ~ "Add independent plots before estimating within-zone variation.",
          any(in_zone & (!full | review)) ~ "Resolve flagged or short cores, then add plots in coverage gaps.",
          TRUE ~ "Add randomized plots in coverage gaps and retain a paired repeat panel."
        )
      )
      if (!length(z)) next
      spacing <- nearest[z]
      if (!any(is.finite(spacing))) {
        spacing <- nearest_distance(grid[z, ], cores)
      }
      chosen <- integer()
      for (k in seq_len(min(plan$focus_areas_per_zone, length(z)))) {
        i <- order(-spacing, grid$cell_id[z], na.last = TRUE)[1L]
        chosen <- c(chosen, z[i])
        gap <- nearest_distance(grid[z, ], grid[z[i], ])
        spacing <- ifelse(is.finite(spacing), pmin(spacing, gap), gap)
        spacing[z %in% chosen] <- -Inf
      }
      focus_rows[[j]] <- sf::st_sf(
        depth_cm = depth, zone = zones$zone[j], order_in_zone = seq_along(chosen),
        gap_m = nearest[chosen], point_sd_mg_c_ha = point_sd[chosen],
        geometry = sf::st_geometry(grid[chosen, ])
      )
    }
    zone_table <- bind_rows(zone_rows) |>
      arrange(full_depth_qc_plots, desc(max_gap_m), desc(area_ha)) |>
      mutate(priority_order = row_number(), .after = zone)
    focus <- do.call(rbind, Filter(Negate(is.null), focus_rows))
    if (is.null(focus)) {
      focus <- sf::st_sf(depth_cm = numeric(), zone = character(), order_in_zone = integer(),
                         gap_m = numeric(), point_sd_mg_c_ha = numeric(),
                         geometry = sf::st_sfc(crs = sf::st_crs(grid)))
    }
    qa <- tibble(
      depth_cm = depth, core_id = records$`Core ID`, plot_id = records$`Plot ID`,
      zone = records$zone, coring_date = records$coring_date, source = records$data_source,
      latitude = records$Latitude, longitude = records$Longitude, stock_mg_c_ha = logged_stock * 10,
      coverage_cm = records$`Stock coverage (cm)`, qc = records$core_check,
      mapped_location = core_id %in% core_ids,
      full_depth = is.finite(stock_mg_c_ha) & stock_mg_c_ha >= 0 &
        is.finite(coverage_cm) & coverage_cm >= depth - 1e-8,
      used_to_fit = mapped_location & is.finite(stock_mg_c_ha) & stock_mg_c_ha > 0 &
        is.finite(coverage_cm) & coverage_cm > 0 & !is.na(method)
    )
    winner_scores <- scores[scores$method %in% method & scores$core_id %in% common_ids, ]
    qa$held_out_prediction_mg_c_ha <- 10 * exp(winner_scores$predicted_log_stock[
      match(qa$core_id, winner_scores$core_id)
    ])
    qa$held_out_error_mg_c_ha <- qa$held_out_prediction_mg_c_ha - qa$stock_mg_c_ha
    summary <- tibble(
      depth_cm = depth, model = method,
      model_selection = if (nrow(ranking)) "Lowest comparable held-out log RMSE" else
        if (!is.na(method)) "Only fitted model; no held-out score" else "No fitted model",
      logged_cores = nrow(qa), full_depth_cores = sum(qa$full_depth),
      qc_flagged_cores = sum(is.na(qa$qc) | qa$qc != "OK"),
      full_depth_qc_mapped_plots = n_distinct(support$`Plot ID`),
      short_cores_in_fit = sum(qa$used_to_fit & !qa$full_depth),
      flagged_cores_in_fit = sum(qa$used_to_fit & (is.na(qa$qc) | qa$qc != "OK")),
      aoi_ha = inputs$area_m2 / 10000, covered_ha = mapped_area / 10000,
      coverage_pct = 100 * mapped_area / inputs$area_m2,
      covered_mean_mg_c_ha = if (mapped_area > 0) total / (mapped_area / 10000) else NA_real_,
      covered_stock_mg_c = total, aoi_stock_mg_c = if (complete_area) total else NA_real_,
      aoi_stock_sd_mg_c = joint_sd,
      aoi_relative_sd_pct = if (is.finite(joint_sd) && total > 0) 100 * joint_sd / total else NA_real_,
      aoi_uncertainty = joint_status,
      held_out_rmse_log = if (nrow(ranking)) ranking$rmse_log[1L] else NA_real_,
      change_status = "Baseline only; change and crediting deduction unestimated"
    )
    bands <- data.frame(stock_mg_c_ha = mean_stock, point_sd_mg_c_ha = point_sd,
                        gap_to_full_depth_core_m = nearest)
    export_mapping_raster(inputs, bands, file.path(out, paste0("mrv_", depth, "cm.tif")))

    # Render each quantity on its own scale and use the existing raster grid.
    map_layer <- function(value, label, palette, mark_focus = FALSE) {
      pixels <- data.frame(x = xy[, 1], y = xy[, 2], value = value)
      subtitle <- if (any(is.finite(value))) NULL else "Unavailable"
      p <- ggplot() +
        geom_tile(data = pixels, aes(x, y, fill = value),
                  width = settings$cell_size_m, height = settings$cell_size_m) +
        geom_sf(data = inputs$aoi, fill = NA, colour = "#364950", linewidth = 0.15) +
        geom_sf(data = support, size = 1.2, shape = 21, fill = "white", stroke = 0.4) +
        scale_fill_viridis_c(option = palette, na.value = "transparent") +
        labs(title = label, subtitle = subtitle, fill = NULL) +
        theme_void(base_size = 11) +
        theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
              plot.margin = margin(8, 12, 4, 12))
      if (mark_focus && nrow(focus)) {
        p <- p + geom_sf(data = focus, shape = 23, fill = "#fff0a6",
                        colour = "#2b2321", size = 2.5, stroke = 0.6)
      }
      suppressMessages(p + coord_sf(crs = sf::st_crs(grid), datum = NA))
    }
    figure <- file.path(out, paste0("mrv_", depth, "cm.png"))
    grDevices::png(figure, width = 1600, height = 640, res = 145)
    tryCatch({
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 3)))
      print(map_layer(mean_stock, "Mapped stock (Mg C/ha)", "viridis"),
            vp = grid::viewport(layout.pos.col = 1))
      print(map_layer(point_sd, "Point SD (Mg C/ha)", "magma"),
            vp = grid::viewport(layout.pos.col = 2))
      print(map_layer(nearest / 1000, "Distance to usable cores (km)", "cividis", TRUE),
            vp = grid::viewport(layout.pos.col = 3))
    }, finally = grDevices::dev.off())
    list(summary = summary, ranking = ranking, qa = qa, zones = zone_table, focus = focus,
         zone_shapes = left_join(zones, zone_table, by = "zone"), figure = figure)
  })

  collect <- function(name) bind_rows(lapply(depth_results, "[[", name))
  summary <- collect("summary")
  validation <- collect("ranking")
  core_qa <- collect("qa")
  zone_table <- collect("zones")
  focus <- do.call(rbind, lapply(depth_results, "[[", "focus"))
  zone_shapes <- do.call(rbind, lapply(depth_results, "[[", "zone_shapes"))
  report_depth <- plan$reporting_depth_cm
  current <- summary[summary$depth_cm == report_depth, ]
  source_label <- paste(unique(stats::na.omit(core_qa$source)), collapse = ", ")
  simulated <- any(grepl("simulat", core_qa$source, ignore.case = TRUE))
  dates <- sort(unique(stats::na.omit(core_qa$coring_date)))
  timing <- if (length(dates)) paste("Coring dates recorded:", paste(range(dates), collapse = " to ")) else
    "Coring dates are missing from the workbook."

  # Normal-approximation planning for an independent panel of paired plots.
  # The input SD is for change at the plot support, in the reporting depth.
  # Use a design-specific variance/power calculation before a formal survey.
  critical <- qnorm((1 + plan$confidence) / 2) + qnorm(plan$power)
  sd_change <- plan$paired_change_sd_mg_c_ha
  n_pairs <- plan$planned_plot_pairs
  target <- plan$meaningful_change_mg_c_ha
  rate <- plan$expected_change_mg_c_ha_yr
  mdd <- critical * sd_change / sqrt(n_pairs)
  required_pairs <- if (is.finite(sd_change) && is.finite(target)) {
    max(2, ceiling((critical * sd_change / target)^2))
  } else NA_real_
  years <- if (is.finite(mdd) && is.finite(rate)) max(1, ceiling(mdd / rate)) else NA_real_
  repeat_year <- if (is.finite(years) && !is.na(baseline_date)) {
    as.integer(format(baseline_date, "%Y")) + years
  } else NA_real_
  actions <- c(
    if (simulated) "The supplied cores are simulated. Establish a measured baseline before making an MRV claim.",
    if (current$short_cores_in_fit > 0 || current$flagged_cores_in_fit > 0)
      "The selected map uses short or flagged cores at the reporting depth. Resolve the flagged records, fill the depth gaps, and refit on comparable cores before using the map for MRV.",
    if (current$qc_flagged_cores > 0)
      "Review the core QA table with the field and lab teams before choosing the repeat panel.",
    if (!is.finite(current$aoi_stock_sd_mg_c))
      "A spatially aggregated model standard deviation is unavailable. Point SD cannot supply a project interval.",
    if (is.na(current$model) || current$coverage_pct < 99.9999)
      "Check unmapped areas and habitat coverage before reporting a full-area stock.",
    if (!is.finite(current$held_out_rmse_log))
      "Collect an independent set of plots to assess prediction error.",
    "Confirm habitat, eligibility and safe access in the mapped focus areas. Use them to define or adjust strata, then select new plots randomly with known inclusion probabilities. Retain plot IDs and draw reserve locations before fieldwork.",
    "For the next campaign, compare the probability-sample estimate with a model-assisted estimate using a map frozen before those samples are collected. Compare them at the same depth and reporting area.",
    "At repeat visits, sample close to the original plots without recoring the original holes. Match season, tidal conditions, depth and lab methods. Agree on a sediment elevation or soil-mass reference when relevant to the selected methodology.",
    if (!is.finite(years)) paste(
      "Run a paired repeat pilot to estimate variability in plot-level stock change.",
      "Set a meaningful change, expected annual change, and planned number of independent plot pairs",
      "before using the timing calculator. Use the monitoring frequency required by the selected methodology."
    ),
    if (is.finite(required_pairs) && is.finite(n_pairs) && n_pairs < required_pairs)
      sprintf("At the chosen power, plan for about %s independent plot pairs to detect %.2f Mg C/ha under the stated assumptions.",
              format(required_pairs, scientific = FALSE), target),
    if (is.finite(years))
      sprintf("The expected change reaches the planning detection limit after about %s whole year(s)%s. Confirm this timing with a design-specific power calculation and the applicable monitoring schedule.",
              format(years, scientific = FALSE), if (is.finite(repeat_year))
                paste0(" (same season around ", repeat_year, ")") else " after baseline"),
    "For a management-effect claim, define the counterfactual and repeat-sampling design before estimating net change. Use the applicable blue carbon methodology to calculate reportable removals and any deduction."
  )
  monitoring <- tibble(
    depth_cm = report_depth, baseline_date = baseline_date, confidence = plan$confidence,
    power = plan$power, meaningful_change_mg_c_ha = target,
    expected_change_mg_c_ha_yr = rate, paired_change_sd_mg_c_ha = sd_change,
    planned_plot_pairs = n_pairs, minimum_detectable_change_mg_c_ha = mdd,
    approximate_pairs_for_target = required_pairs, approximate_years_after_baseline = years,
    approximate_repeat_year = repeat_year,
    status = if (is.finite(years)) "Planning estimate only" else "Collect or enter planning inputs"
  )
  write_csv <- function(data, filename) readr::write_csv(data, file.path(out, filename), na = "")
  write_csv(summary, "mrv_statistics.csv")
  write_csv(validation, "model_validation.csv")
  write_csv(core_qa, "core_qa.csv")
  write_csv(zone_table, "zone_guidance.csv")
  write_csv(monitoring, "monitoring_plan.csv")
  write_csv(tibble(order = seq_along(actions), action = actions), "next_steps.csv")
  focus_ll <- sf::st_transform(focus, 4326)
  coords <- sf::st_coordinates(focus_ll)
  write_csv(mutate(sf::st_drop_geometry(focus_ll), longitude = coords[, 1], latitude = coords[, 2]),
            "sampling_focus_areas.csv")
  geopackage <- file.path(out, "sampling_guidance.gpkg")
  sf::st_write(zone_shapes, geopackage, layer = "zone_guidance", delete_layer = TRUE, quiet = TRUE)
  sf::st_write(focus, geopackage, layer = "sampling_focus_areas", delete_layer = TRUE, quiet = TRUE)

  fmt <- function(x, digits = 1) {
    ifelse(is.finite(x), formatC(x, digits = digits, format = "f", big.mark = ","), "Unavailable")
  }
  tags <- htmltools::tags
  html_table <- function(data) {
    tags$table(
      tags$thead(tags$tr(lapply(names(data), tags$th))),
      tags$tbody(lapply(seq_len(nrow(data)), function(i) {
        tags$tr(lapply(data, function(column) tags$td(as.character(column[i]))))
      }))
    )
  }
  summary_table <- summary |>
    transmute(Depth = paste0("0–", depth_cm, " cm"), Model = coalesce(model, "Unavailable"),
              `Cores reaching depth` = full_depth_cores, `Usable mapped plots` = full_depth_qc_mapped_plots,
              `Coverage (%)` = fmt(coverage_pct, 2), `Mapped mean (Mg C/ha)` = fmt(covered_mean_mg_c_ha),
              `AOI stock (Mg C)` = fmt(aoi_stock_mg_c, 0),
              `AOI stock SD (Mg C)` = fmt(aoi_stock_sd_mg_c, 0))
  focus_table <- zone_table |>
    filter(depth_cm == report_depth) |>
    transmute(Zone = zone, `Area (ha)` = fmt(area_ha, 0), `Usable plots` = full_depth_qc_plots,
              `Largest gap (m)` = fmt(max_gap_m, 0),
              `Observed plot mean (Mg C/ha)` = fmt(observed_plot_mean_mg_c_ha), Action = action)
  validation_table <- validation |>
    transmute(Depth = paste0("0–", depth_cm, " cm"), Model = method, Selected = if_else(selected, "Yes", ""),
              `Held-out plots` = held_out_plots, `Log RMSE` = fmt(rmse_log, 3),
              `RMSE (Mg C/ha)` = fmt(rmse_mg_c_ha), `Bias (Mg C/ha)` = fmt(bias_mg_c_ha))
  main_figure <- depth_results[[match(report_depth, depths)]]$figure
  image_tag <- function(path, depth) {
    tags$img(src = base64enc::dataURI(file = path, mime = "image/png"), width = "100%",
             alt = paste0("Mapped stock, point uncertainty and distance to usable cores at 0–", depth, " cm"))
  }
  htmltools::save_html(tags$html(
    lang = "en", tags$head(tags$meta(charset = "utf-8"),
      tags$title("Blue carbon MRV and sampling plan"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$style("
        body { color:#172d36; background:#f6f8f7; font:16px/1.55 system-ui,sans-serif; margin:0 }
        main { max-width:1080px; margin:auto; padding:2.5rem clamp(1rem,4vw,3rem); background:white }
        h1,h2,h3 { line-height:1.2; color:#154d4c } h2 { margin-top:2.4rem }
        table { width:100%; border-collapse:collapse; font-size:.92rem }
        th,td { padding:.55rem .65rem; text-align:left; border-bottom:1px solid #d9e4e1 }
        th { background:#eef4f1 } li { margin-bottom:.55rem }
        .table { overflow-x:auto } .note { padding:.85rem 1rem; background:#eef4f1 }
        img { display:block; height:auto } details { margin-top:1rem }"
      )
    ),
    tags$body(tags$main(
      tags$h1("Blue carbon MRV and sampling plan"),
      tags$p(class = "note", paste0(source_label, ". Baseline stock assessment. Stock change, ",
                                   "net removals and credit deductions are unestimated.")),
      tags$p(timing),
      tags$h2("Baseline stock"),
      tags$div(class = "table", html_table(summary_table)),
      tags$p("The AOI total weights each estimate by the area of the grid cell inside the boundary. A full-area total requires complete map coverage. Usable plots have a geolocated core reaching the target depth with no Script 1 QC flag."),
      tags$h2(paste0("Where to sample next at 0–", report_depth, " cm")),
      image_tag(main_figure, report_depth),
      tags$p("White circles mark usable cores. Diamonds identify coverage gaps for field reconnaissance. They are screening locations. Use a probability design to choose the actual plots. IDW point SD is a residual-based proxy and tends to rise with mapped stock. Treat the gap map as the primary coverage guide."),
      tags$div(class = "table", html_table(focus_table)),
      tags$p("Zone plot means describe the sampled plots. Repeated cores at one plot are averaged before comparing plots."),
      tags$h2("Repeat sampling"),
      tags$p(sprintf("Planning uses a %.0f%% two-sided confidence level and %.0f%% power. It assumes independent paired plots, a common SD of paired differences and linear change. It estimates detectability, which is separate from a crediting rule.",
                     100 * plan$confidence, 100 * plan$power)),
      if (is.finite(mdd)) tags$p(paste("Approximate minimum detectable change:", fmt(mdd, 2), "Mg C/ha.")),
      tags$ol(lapply(actions, tags$li)),
      tags$h2("Model checks"),
      if (nrow(validation_table)) tags$div(class = "table", html_table(validation_table)),
      tags$p("Model selection compares IDW and Matérn using common held-out plots and log-stock error. These folds also tuned the candidates, so the scores are exploratory. Bias and RMSE in stock units use back-transformed log predictions. The regional Bayesian and AlphaEarth transfer surfaces lack matching held-out scores in this workflow."),
      if (any(summary$short_cores_in_fit > 0 | summary$flagged_cores_in_fit > 0))
        tags$p("The workshop fit used short or flagged cores at one or more depths. See the core QA table for the affected records."),
      lapply(which(depths != report_depth), function(i) tags$details(
        tags$summary(paste0("Map at 0–", depths[i], " cm")),
        image_tag(depth_results[[i]]$figure, depths[i])
      ))
    ))
  ), file = file.path(out, "mrv_report.html"))
  message("MRV report, tables and spatial outputs saved to ", out, ".")
})
