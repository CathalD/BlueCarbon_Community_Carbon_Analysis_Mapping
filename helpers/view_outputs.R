# Build the same map widgets for Script 4 and for viewing results directly in R.
.build_output_map <- function(grid, depth_cm, source_label, local_folder, bayesian_folder,
                              transfer_folder, by_zone = NULL, prior = NULL, bayesian = NULL,
                              embedding_map = NULL, transfer_only = FALSE, layer = NULL) {
  surface_stats <- list()
  basemap <- "Esri.WorldImagery"
  suffix <- paste0(depth_cm, "cm")
  layers <- list()
  layer_values <- list()
  layer_kind <- character()
  add_band <- function(raster, band, label, kind) {
    if (!band %in% names(raster)) stop("Missing band ", band, " for ", label, ".")
    layer <- raster[[band]]
    if (!terra::compareGeom(layer, grid$template, stopOnError = FALSE)) {
      stop("The raster grid differs from the expected mapping grid for ", label, ".")
    }
    z <- terra::values(layer, mat = FALSE)[grid$cells$cell_id]
    if (!any(is.finite(z))) return(invisible(NULL))
    layers[[label]] <<- layer
    layer_values[[label]] <<- z
    layer_kind[label] <<- kind
    use <- is.finite(z)
    surface_stats[[length(surface_stats) + 1L]] <<- tibble::tibble(
      depth_cm = depth_cm, layer = label, value_type = kind, cells = sum(use),
      coverage_pct = 100 * sum(grid$cells$area_m2[use]) / grid$area_m2,
      min_kg_c_m2 = min(z[use]), max_kg_c_m2 = max(z[use]),
      area_weighted_mean_kg_c_m2 = weighted.mean(z[use], grid$cells$area_m2[use])
    )
  }
  if (transfer_only) {
    estimate <- terra::rast(file.path(transfer_folder, paste0("transfer_", suffix, ".tif")))
    sd_path <- file.path(transfer_folder, paste0("transfer_uncertainty_", suffix, ".tif"))
    uncertainty <- if (file.exists(sd_path)) terra::rast(sd_path) else estimate
    correction <- estimate[["posterior_mean_kg_c_m2"]] - estimate[["prior_mean_kg_c_m2"]]
    names(correction) <- "correction"
    add_band(estimate, "prior_mean_kg_c_m2", "Janousek prior", "Estimate")
    add_band(correction, "correction", "Local Matérn correction", "Stock change")
    add_band(estimate, "posterior_mean_kg_c_m2", "Updated estimate", "Estimate")
    add_band(uncertainty, "prior_sd_kg_c_m2", "Prior uncertainty (SD)", "Point SD")
    add_band(uncertainty, "posterior_sd_kg_c_m2", "Updated uncertainty (SD)", "Point SD")
  } else {
    for (method in c("area_based", "idw", "matern")) {
      path <- file.path(local_folder, paste0(method, "_", suffix, ".tif"))
      if (!file.exists(path)) next
      raster <- terra::rast(path)
      label <- c(area_based = "Pooled area mean", idw = "IDW", matern = "Local Matérn")[method]
      add_band(raster, "mean_kg_c_m2", paste(label, "estimate", sep = " · "), "Estimate")
      if ("sd_kg_c_m2" %in% names(raster)) {
        add_band(raster, "sd_kg_c_m2", paste(label, "SD", sep = " · "), "Point SD")
      }
    }
    if (!is.null(bayesian) && any(prior$depth_cm == depth_cm & prior$status == "Fitted")) {
      raster <- terra::rast(file.path(bayesian_folder, paste0("posterior_", suffix, ".tif")))
      local_band <- layers[["Local Matérn · estimate"]]
      if (is.null(local_band) || !isTRUE(all.equal(
        terra::values(raster[["local_matern_mean_kg_c_m2"]], mat = FALSE),
        terra::values(local_band, mat = FALSE), tolerance = 1e-8))) {
        stop("The Bayesian output differs from the local Matérn run. Rerun the Bayesian script.")
      }
      add_band(raster, "prior_mean_kg_c_m2", "Janousek prior · estimate", "Estimate")
      add_band(raster, "prior_sd_kg_c_m2", "Janousek prior · SD", "Point SD")
      add_band(raster, "posterior_mean_kg_c_m2", "Bayesian posterior · estimate", "Estimate")
      add_band(raster, "posterior_sd_kg_c_m2", "Bayesian posterior · SD", "Point SD")
    }
  }
  if (!length(layers)) return(NULL)

  zones <- if (!transfer_only) {
    grid$zones |> dplyr::left_join(dplyr::filter(by_zone, depth_cm == .env$depth_cm) |>
      dplyr::select(zone, n_plots, mean_kg_c_m2, total_mg_c, total_se_mg_c), by = "zone")
  }
  has_zone_stock <- any(is.finite(zones$mean_kg_c_m2))
  values_for <- function(kind) {
    unlist(layer_values[layer_kind == kind], use.names = FALSE)
  }
  make_scale <- function(values, palette) {
    highest <- max(values, na.rm = TRUE)
    breaks <- pretty(c(0, max(highest, 0.01)), n = 5)
    list(breaks = breaks, colours = hcl.colors(length(breaks) - 1L, palette))
  }
  scales <- list(Estimate = make_scale(c(values_for("Estimate"), zones$mean_kg_c_m2), "Viridis"))
  if (any(layer_kind == "Point SD")) scales[["Point SD"]] <- make_scale(values_for("Point SD"), "YlOrRd")
  if (any(layer_kind == "Stock change")) {
    limit <- max(0.01, abs(values_for("Stock change")), na.rm = TRUE)
    positive_breaks <- pretty(c(0, limit), n = 4)
    breaks <- c(-rev(positive_breaks[-1L]), positive_breaks)
    scales[["Stock change"]] <- list(breaks = breaks, colours = hcl.colors(length(breaks) - 1L, "Blue-Red 3"))
  }

  map <- mapview::mapview(grid$aoi, layer.name = "AOI boundary", color = "#263746",
                          alpha.regions = 0, map.types = basemap, legend = FALSE, homebutton = FALSE,
                          popup = FALSE)
  for (label in names(layers)) {
    scale <- scales[[layer_kind[label]]]
    map <- map + mapview::mapview(layers[[label]], layer.name = label, at = scale$breaks,
                                  col.regions = scale$colours, na.color = "transparent",
                                  method = "ngb", map.types = basemap, legend = FALSE, homebutton = FALSE,
                                  query.type = "click", query.prefix = "kg C/m²")
  }
  surface_names <- names(layers)
  if (has_zone_stock) {
    zone_label <- "Zone area means · estimate"
    zone_popup <- zones |> dplyr::filter(is.finite(mean_kg_c_m2)) |>
      dplyr::select(Zone = zone, Plots = n_plots, `Stock (kg C/m²)` = mean_kg_c_m2,
             `Total (Mg C)` = total_mg_c, `Sampling SE (Mg C)` = total_se_mg_c)
    zone_stocks <- unique(zone_popup$`Stock (kg C/m²)`)
    zone_palette <- if (length(zone_stocks) == 1L) {
      scales$Estimate$colours[findInterval(zone_stocks, scales$Estimate$breaks, all.inside = TRUE)]
    } else function(n) hcl.colors(n, "Viridis")
    map <- map + mapview::mapview(zone_popup, zcol = "Stock (kg C/m²)", layer.name = zone_label,
                                  at = scales$Estimate$breaks, col.regions = zone_palette,
                                  color = "#263746", na.color = "transparent", alpha.regions = 0.85,
                                  map.types = basemap, legend = FALSE, homebutton = FALSE)
    surface_names <- c(surface_names, zone_label)
  }
  stock_column <- paste0("C Stock 0-", suffix)
  points <- grid$cores |>
    dplyr::select(`Core ID`, `Plot ID`, Zone = zone_assigned, `Stock (kg C/m²)` = dplyr::all_of(stock_column))
  map <- map + mapview::mapview(points, layer.name = "Core locations", label = points$`Core ID`,
                                color = "#263746", col.regions = "white", alpha.regions = 1,
                                cex = 5, map.types = basemap, legend = FALSE, homebutton = FALSE)
  widget <- map@map
  if (transfer_only) {
    widget <- leaflet::addRasterImage(widget, embedding_map, method = "ngb", group = "AlphaEarth embeddings")
    surface_names <- c("AlphaEarth embeddings", surface_names)
  }
  default_surface <- if (transfer_only) "Updated estimate" else "Local Matérn · estimate"
  if (!default_surface %in% surface_names) default_surface <- surface_names[1L]
  if (!is.null(layer)) {
    if (!is.character(layer) || length(layer) != 1L || !layer %in% surface_names) {
      stop("Choose a layer from ", paste(surface_names, collapse = ", "), ".", call. = FALSE)
    }
    default_surface <- layer
  }
  widget <- widget |>
    leaflet::removeLayersControl() |>
    leaflet::hideGroup(setdiff(surface_names, default_surface)) |>
    leaflet::addLayersControl(baseGroups = if (transfer_only) character() else surface_names,
                               overlayGroups = c("AOI boundary", "Core locations"),
                               options = leaflet::layersControlOptions(collapsed = FALSE))
  for (kind in names(scales)) {
    scale <- scales[[kind]]
    labels <- paste(head(scale$breaks, -1), tail(scale$breaks, -1), sep = "–")
    widget <- leaflet::addLegend(widget, position = if (kind == "Estimate") "bottomleft" else "bottomright",
                                 colors = scale$colours, labels = labels, opacity = 1,
                                 className = paste("info legend", paste0("map-legend-", match(kind, names(scales)))),
                                 title = paste0(kind, " (kg C/m²)"))
  }
  layer_menu <- NULL
  if (transfer_only) {
    layer_menu <- htmltools::tags$label("Map layer ",
      htmltools::tags$select(class = "transfer-map-layer",
        lapply(surface_names, function(label) {
          htmltools::tags$option(label, value = label,
                                 selected = if (label == default_surface) "selected" else NULL)
        })
      )
    )
    widget <- htmlwidgets::onRender(widget, '
      function(el, x, data) {
        const map = this;
        const menu = el.parentElement.querySelector(".transfer-map-layer");
        const showSelectedLayer = () => {
          const selected = menu.value;
          const visibleLegend = data.legends[data.groups.indexOf(selected)];
          data.groups.forEach(name => {
            const layer = map.layerManager.getLayerGroup(name, false);
            if (name === selected && !map.hasLayer(layer)) map.addLayer(layer);
            if (name !== selected && map.hasLayer(layer)) map.removeLayer(layer);
          });
          el.querySelectorAll(".legend").forEach(legend => {
            legend.hidden = !legend.classList.contains(visibleLegend);
          });
        };
        menu.addEventListener("change", showSelectedLayer);
        showSelectedLayer();
      }',
      data = list(groups = surface_names,
                  legends = paste0("map-legend-", match(unname(layer_kind[surface_names]), names(scales)))))
  }
  widget <- htmlwidgets::prependContent(widget, htmltools::tags$div(
    style = "padding: 10px 18px; font: 15px system-ui, sans-serif; color: #263746;",
    htmltools::tags$strong(paste(if (transfer_only) "Transfer Learning with embeddings kernel" else "Sediment carbon",
                                paste0("0–", depth_cm, " cm"))),
    htmltools::tags$p(paste0(source_label, ". Click a stock layer to see its value.")),
    layer_menu
  ))
  list(widget = widget, surfaces = surface_names, default_surface = default_surface,
       stats = dplyr::bind_rows(surface_stats))
}

.transfer_view_inputs <- function(transfer, spatial, settings, transfer_folder, transfer_inputs_folder,
                                  read_grid) {
  transfer_spacing <- unique(transfer$cell_size_m)
  if (length(transfer_spacing) != 1L || !is.finite(transfer_spacing)) {
    stop("The transfer summary must specify one grid spacing. Rerun Script 3.")
  }
  first_transfer <- terra::rast(file.path(transfer_folder, paste0("transfer_", min(transfer$depth_cm), "cm.tif")))
  transfer_spatial <- spatial
  if (!terra::compareGeom(first_transfer[[1L]], spatial$template, stopOnError = FALSE)) {
    transfer_spatial <- read_grid(settings$aoi_file, settings$mapping_file,
                                  terra::crs(first_transfer), transfer_spacing)
  }
  if (!all(abs(transfer$aoi_ha - transfer_spatial$area_m2 / 10000) < 1e-5)) {
    stop("The transfer summary uses a different AOI. Rerun Script 3 with the reporting boundary.")
  }

  # Use the same embedding values as the model, on its prediction grid.
  rgb_columns <- c("var2", "var17", "var10")
  embeddings <- .read_view_output(transfer_inputs_folder, "pred_points.csv",
                                  col_select = dplyr::all_of(c("Core ID", "embedding_year", rgb_columns)))
  embedding_cells <- match(paste0("cell_", transfer_spatial$cells$cell_id), embeddings[["Core ID"]])
  if (anyNA(embedding_cells) || anyDuplicated(embeddings[["Core ID"]]) ||
      !all(embeddings$embedding_year[embedding_cells] %in% unique(transfer$embedding_year))) {
    stop("Prediction point embeddings differ from the transfer run. Rerun Script 3.")
  }
  rgb <- as.matrix(embeddings[embedding_cells, rgb_columns])
  rgb[] <- round(pmin(1, pmax(0, (rgb + 0.3) / 0.6)) * 255)
  embedding_map <- terra::rast(transfer_spatial$template, nlyrs = length(rgb_columns))
  rgb_cells <- matrix(NA_real_, terra::ncell(embedding_map), length(rgb_columns))
  rgb_cells[transfer_spatial$cells$cell_id, ] <- rgb
  terra::values(embedding_map) <- rgb_cells
  terra::RGB(embedding_map) <- 1:3
  list(grid = transfer_spatial, embedding_map = embedding_map)
}

.read_view_output <- function(folder, filename, ...) {
  path <- file.path(folder, filename)
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run the analysis scripts first.", call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, ...)
}

# Open a Mapview map in the RStudio Viewer. Both map sets include a layer menu.
view_map <- function(depth_cm = 30, model = c("local", "transfer"), layer = NULL,
                     project_dir = ".") {
  model <- match.arg(model)
  if (length(depth_cm) != 1L || !is.numeric(depth_cm) || !is.finite(depth_cm) || depth_cm <= 0) {
    stop("Supply one depth in centimetres, for example 30.", call. = FALSE)
  }
  local_folder <- file.path(project_dir, "outputs", "script2_local")
  bayesian_folder <- file.path(project_dir, "outputs", "script2_bayesian")
  transfer_folder <- file.path(project_dir, "outputs", "script3_covariate_transfer")
  run_path <- file.path(local_folder, "mapping_run.rds")
  if (!file.exists(run_path)) stop("Run Script_2_local_mapping.R first.", call. = FALSE)
  run <- readRDS(run_path)
  settings <- run$settings
  project_path <- function(path) {
    if (grepl("^(/|~|[A-Za-z]:)", path)) path.expand(path) else file.path(project_dir, path)
  }
  settings$aoi_file <- project_path(settings$aoi_file)
  settings$mapping_file <- project_path(settings$mapping_file)
  input_md5 <- unname(tools::md5sum(c(settings$aoi_file, settings$mapping_file)))
  if (anyNA(input_md5) || !identical(input_md5, run$input_md5)) {
    stop("The boundary or mapping data have changed. Rerun the local mapping script.", call. = FALSE)
  }
  sys.source(file.path(project_dir, "helpers", "local_mapping_spatial.R"), envir = environment())
  grid <- read_mapping_inputs(settings$aoi_file, settings$mapping_file,
                              settings$analysis_crs, settings$cell_size_m)
  source_label <- paste(unique(stats::na.omit(grid$cores$data_source)), collapse = ", ")
  if (!nzchar(source_label)) source_label <- "Carbon stock estimates"

  if (model == "transfer") {
    transfer <- .read_view_output(transfer_folder, "transfer_area_summary.csv")
    available_depths <- sort(unique(transfer$depth_cm))
    if (!depth_cm %in% available_depths) {
      stop("Available transfer depths (cm): ", paste(available_depths, collapse = ", "), ".", call. = FALSE)
    }
    transfer_view <- .transfer_view_inputs(transfer, grid, settings, transfer_folder,
                                           file.path(project_dir, "data", "covariate_transfer"),
                                           read_grid = read_mapping_inputs)
    result <- .build_output_map(transfer_view$grid, depth_cm, source_label, local_folder,
                                bayesian_folder, transfer_folder, embedding_map = transfer_view$embedding_map,
                                transfer_only = TRUE, layer = layer)
  } else {
    available_depths <- sort(unique(settings$depths_cm))
    if (!depth_cm %in% available_depths) {
      stop("Available local depths (cm): ", paste(available_depths, collapse = ", "), ".", call. = FALSE)
    }
    by_zone <- .read_view_output(local_folder, "area_summary.csv") |>
      dplyr::filter(method == "By zone")
    prior <- bayesian <- NULL
    if (file.exists(file.path(bayesian_folder, "run_settings.csv"))) {
      bayesian_run <- .read_view_output(bayesian_folder, "run_settings.csv")
      values <- stats::setNames(bayesian_run$value, bayesian_run$setting)
      if (!identical(unname(values[c("aoi_md5", "mapping_md5")]), input_md5)) {
        stop("The Bayesian maps use different inputs. Rerun the Bayesian script.", call. = FALSE)
      }
      prior <- .read_view_output(bayesian_folder, "prior_parameters.csv")
      bayesian <- .read_view_output(bayesian_folder, "aoi_totals.csv")
    }
    result <- .build_output_map(grid, depth_cm, source_label, local_folder, bayesian_folder,
                                transfer_folder, by_zone, prior, bayesian, layer = layer)
  }
  if (is.null(result)) stop("No map layers were generated for that depth.", call. = FALSE)
  result$widget
}

# Return a table for printing in the console or opening with View().
view_table <- function(output = c("totals", "cores", "transfer", "coverage", "site", "zones",
                                  "local_layers", "transfer_layers"),
                       depth_cm = NULL, project_dir = ".") {
  output <- match.arg(output)
  file <- c(totals = "mapping_totals", cores = "core_stocks", transfer = "transfer_summary",
            coverage = "mapping_coverage", site = "site_summary", zones = "area_by_zone",
            local_layers = "map_surface_summary", transfer_layers = "transfer_map_surface_summary")
  result <- .read_view_output(file.path(project_dir, "outputs", "script4_reporting"),
                              paste0(file[[output]], ".csv"))
  if (!is.null(depth_cm)) {
    if (!"depth_cm" %in% names(result)) {
      stop("This table does not have a depth column. Read it without depth_cm.", call. = FALSE)
    }
    available_depths <- sort(unique(result$depth_cm))
    if (!depth_cm %in% available_depths) {
      stop("Available depths (cm): ", paste(available_depths, collapse = ", "), ".", call. = FALSE)
    }
    result <- result[result$depth_cm == depth_cm, ]
  }
  result
}

# Return a ggplot for the main stock results. Add ggplot layers as usual.
view_plot <- function(output = c("totals", "cores", "transfer"), depth_cm = NULL, project_dir = ".") {
  output <- match.arg(output)
  data <- view_table(output, if (output == "cores") NULL else depth_cm, project_dir)
  depth <- function(x) factor(paste0("0–", x, " cm"), levels = paste0("0–", sort(unique(x)), " cm"))
  theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  if (output == "cores") {
    data <- tidyr::pivot_longer(data, cols = dplyr::matches("^C Stock 0-[0-9]+cm$"),
                                names_to = "depth_cm", names_pattern = "C Stock 0-(\\d+)cm",
                                names_transform = list(depth_cm = as.integer), values_to = "stock")
    if (!is.null(depth_cm)) {
      if (!depth_cm %in% data$depth_cm) stop("That core depth is not available.", call. = FALSE)
      data <- data[data$depth_cm == depth_cm, ]
    }
    data$depth <- depth(data$depth_cm)
    return(ggplot2::ggplot(data[is.finite(data$stock), ], ggplot2::aes(zone, stock, colour = zone)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 42)) +
      ggplot2::facet_wrap(ggplot2::vars(depth), scales = "free_y") +
      ggplot2::labs(title = "Carbon stock in the study cores", x = "Zone", y = "Carbon stock (kg C/m²)") +
      ggplot2::guides(colour = "none") + theme)
  }

  data$depth <- depth(data$depth_cm)
  if (output == "transfer") {
    data$estimate <- factor(data$surface, levels = c("prior", "posterior"),
                            labels = c("Janousek prior", "Updated estimate"))
    return(ggplot2::ggplot(data, ggplot2::aes(estimate, mean_covered_kg_c_m2)) +
      ggplot2::geom_point(size = 3, colour = "#176A61") +
      ggplot2::facet_wrap(ggplot2::vars(depth)) +
      ggplot2::labs(title = "Carbon stock with and without the local update",
                    x = "Estimate", y = "Mean carbon stock (kg C/m²)") + theme)
  }

  ggplot2::ggplot(data[is.finite(data$total_mg_c), ], ggplot2::aes(total_mg_c / 1000, method)) +
    ggplot2::geom_errorbar(data = data[is.finite(data$total_mg_c) & is.finite(data$uncertainty_mg_c), ],
                          ggplot2::aes(xmin = (total_mg_c - uncertainty_mg_c) / 1000,
                                       xmax = (total_mg_c + uncertainty_mg_c) / 1000,
                                       colour = uncertainty_type),
                          orientation = "y", width = 0.22) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(ggplot2::vars(depth), scales = "free_x") +
    ggplot2::labs(title = "Carbon stock in the study area", x = "Carbon stock (thousand tonnes)",
                  y = "Method", colour = "Error bars") + theme
}
