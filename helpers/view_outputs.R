# Relative uncertainty is defined only for a positive estimate and a valid SD.
.stock_cv <- function(stock, sd) {
  if (inherits(stock, "SpatRaster")) {
    return(terra::ifel(is.finite(stock) & stock > 0 & is.finite(sd) & sd >= 0,
                       sd / stock, NA_real_))
  }
  valid <- is.finite(stock) & stock > 0 & is.finite(sd) & sd >= 0
  cv <- rep(NA_real_, length(stock))
  cv[valid] <- sd[valid] / stock[valid]
  cv
}

# Build the same Leaflet map for R, standalone HTML and the Quarto report.
.build_output_map <- function(grid, depth_cm, source_label, local_folder, bayesian_folder,
                              transfer_folder, by_zone = NULL, prior = NULL, bayesian = NULL,
                              embedding_map = NULL, transfer_only = FALSE, layer = NULL,
                              pooled = NULL) {
  suffix <- paste0(depth_cm, "cm")
  layers <- list()
  layer_values <- list()
  layer_kind <- character()
  surface_stats <- list()
  add_raster <- function(raster, label, kind) {
    if (!terra::compareGeom(raster, grid$template, stopOnError = FALSE)) {
      stop("The raster grid differs from the expected mapping grid for ", label, ".")
    }
    z <- terra::values(raster, mat = FALSE)[grid$cells$cell_id]
    use <- is.finite(z)
    if (!any(use)) return(invisible(NULL))
    layers[[label]] <<- raster
    layer_values[[label]] <<- z
    layer_kind[label] <<- kind
    values <- z[use]
    area <- grid$cells$area_m2[use]
    summary <- c(min(values), max(values), weighted.mean(values, area))
    stock_summary <- if (kind == "CV") rep(NA_real_, 3) else summary
    surface_stats[[label]] <<- tibble::tibble(
      depth_cm = depth_cm, layer = label, value_type = kind,
      units = if (kind == "CV") "SD / mean" else "kg C/m²",
      cells = sum(use), coverage_pct = 100 * sum(area) / grid$area_m2,
      min_value = summary[1], max_value = summary[2], area_weighted_mean = summary[3],
      min_kg_c_m2 = stock_summary[1], max_kg_c_m2 = stock_summary[2],
      area_weighted_mean_kg_c_m2 = stock_summary[3]
    )
  }
  add_stock <- function(stock, sd, estimate_label, sd_label, cv_label) {
    add_raster(stock, estimate_label, "Estimate")
    if (is.null(sd)) return(invisible(NULL))
    add_raster(sd, sd_label, "Point SD")
    add_raster(.stock_cv(stock, sd), cv_label, "CV")
  }

  if (transfer_only) {
    estimate <- terra::rast(file.path(transfer_folder, paste0("transfer_", suffix, ".tif")))
    sd_path <- file.path(transfer_folder, paste0("transfer_uncertainty_", suffix, ".tif"))
    uncertainty <- if (file.exists(sd_path)) terra::rast(sd_path) else estimate
    for (surface in c("prior", "posterior")) {
      label <- if (surface == "prior") "Prior" else "Updated"
      add_stock(estimate[[paste0(surface, "_mean_kg_c_m2")]],
                uncertainty[[paste0(surface, "_sd_kg_c_m2")]],
                if (surface == "prior") "Janousek prior" else "Updated estimate",
                paste0(label, " uncertainty (SD)"), paste0(label, " uncertainty (CV)"))
    }
    add_raster(estimate[["posterior_mean_kg_c_m2"]] - estimate[["prior_mean_kg_c_m2"]],
               "Local Matérn correction", "Stock change")
  } else {
    for (method in c("area_based", "idw", "matern")) {
      path <- file.path(local_folder, paste0(method, "_", suffix, ".tif"))
      if (!file.exists(path)) next
      raster <- terra::rast(path)
      label <- c(area_based = "Pooled area mean", idw = "IDW", matern = "Local Matérn")[[method]]
      stock <- raster[["mean_kg_c_m2"]]
      sd <- if ("sd_kg_c_m2" %in% names(raster)) raster[["sd_kg_c_m2"]]
      if (method == "area_based" && "plot_sd_kg_c_m2" %in% names(pooled)) {
        plot_sd <- pooled$plot_sd_kg_c_m2[match(depth_cm, pooled$depth_cm)]
        if (is.finite(plot_sd)) sd <- terra::ifel(is.finite(stock), plot_sd, NA_real_)
      }
      add_stock(stock, sd, paste(label, "estimate", sep = " · "),
                paste(label, "SD", sep = " · "), paste(label, "CV", sep = " · "))
    }
    if (!is.null(bayesian) && any(prior$depth_cm == depth_cm & prior$status == "Fitted")) {
      raster <- terra::rast(file.path(bayesian_folder, paste0("posterior_", suffix, ".tif")))
      local_band <- layers[["Local Matérn · estimate"]]
      if (is.null(local_band) || !isTRUE(all.equal(
        terra::values(raster[["local_matern_mean_kg_c_m2"]], mat = FALSE),
        terra::values(local_band, mat = FALSE), tolerance = 1e-8))) {
        stop("The Bayesian output differs from the local Matérn run. Rerun the Bayesian script.")
      }
      for (surface in c("prior", "posterior")) {
        label <- if (surface == "prior") "Janousek prior" else "Bayesian posterior"
        add_stock(raster[[paste0(surface, "_mean_kg_c_m2")]], raster[[paste0(surface, "_sd_kg_c_m2")]],
                  paste(label, "estimate", sep = " · "), paste(label, "SD", sep = " · "),
                  paste(label, "CV", sep = " · "))
      }
    }
  }
  if (!length(layers)) return(NULL)

  zones <- if (!transfer_only && nrow(grid$zones) && !is.null(by_zone)) {
    grid$zones |>
      dplyr::left_join(dplyr::filter(by_zone, depth_cm == .env$depth_cm) |>
        dplyr::select(zone, n_plots, mean_kg_c_m2, dplyr::any_of("plot_sd_kg_c_m2"),
                      total_mg_c, total_se_mg_c), by = "zone") |>
      dplyr::filter(is.finite(mean_kg_c_m2)) |>
      sf::st_transform(4326)
  }
  if (!is.null(zones) && nrow(zones)) {
    if (!"plot_sd_kg_c_m2" %in% names(zones)) zones$plot_sd_kg_c_m2 <- NA_real_
    zones$cv <- .stock_cv(zones$mean_kg_c_m2, zones$plot_sd_kg_c_m2)
  }
  palettes <- c(Estimate = "Viridis", "Point SD" = "YlOrRd", CV = "YlOrBr", "Stock change" = "Blue-Red 3")
  titles <- c(Estimate = "Stock estimate (kg C/m²)", "Point SD" = "Standard deviation (kg C/m²)",
              CV = "Uncertainty CV (SD / mean)", "Stock change" = "Stock change (kg C/m²)")
  scales <- list()
  for (kind in unique(layer_kind)) {
    values <- c(unlist(layer_values[layer_kind == kind], use.names = FALSE),
                switch(kind, Estimate = zones$mean_kg_c_m2, "Point SD" = zones$plot_sd_kg_c_m2, CV = zones$cv))
    limit <- max(pretty(c(0, max(0.01, abs(values), na.rm = TRUE) * 1.000001), n = 5))
    domain <- c(if (kind == "Stock change") -limit else 0, limit)
    scales[[kind]] <- list(domain = domain,
      palette = leaflet::colorNumeric(grDevices::hcl.colors(256, palettes[[kind]]),
                                      domain = domain, na.color = "transparent"))
  }

  aoi <- sf::st_transform(grid$aoi, 4326)
  bbox <- sf::st_bbox(aoi)
  widget <- leaflet::leaflet(height = 560) |>
    leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery) |>
    leaflet::addPolygons(data = aoi, group = "AOI boundary", color = "#253746",
                         fill = FALSE, weight = 2, options = leaflet::pathOptions(interactive = FALSE)) |>
    leaflet::fitBounds(unname(bbox[["xmin"]]), unname(bbox[["ymin"]]),
                       unname(bbox[["xmax"]]), unname(bbox[["ymax"]]))
  click_values <- list()
  for (label in names(layers)) {
    raster <- leaflet::projectRasterForLeaflet(layers[[label]], "ngb")
    widget <- leaflet::addRasterImage(widget, raster, group = label, project = FALSE,
                                      colors = scales[[layer_kind[[label]]]]$palette, opacity = 0.88)
    # Store click values in zlib format and decode only the selected layer in the browser.
    values <- writeBin(terra::values(raster, mat = FALSE), raw(), size = 4L, endian = "little")
    click_values[[label]] <- list(
      extent = unname(as.vector(terra::ext(raster))), nrow = terra::nrow(raster), ncol = terra::ncol(raster),
      data = base64enc::base64encode(memCompress(values, type = "gzip"))
    )
  }

  display_number <- function(x) ifelse(is.finite(x), formatC(x, digits = 2, format = "f"), "Unavailable")
  if (!is.null(zones) && nrow(zones)) {
    zone_popup <- sprintf("<strong>%s</strong><br>Plots: %s<br>Stock: %s kg C/m²<br>Plot SD: %s kg C/m²<br>CV: %s",
                          htmltools::htmlEscape(zones$zone), zones$n_plots,
                          display_number(zones$mean_kg_c_m2), display_number(zones$plot_sd_kg_c_m2),
                          display_number(zones$cv))
    for (kind in intersect(c("Estimate", "Point SD", "CV"), names(scales))) {
      values <- switch(kind, Estimate = zones$mean_kg_c_m2, "Point SD" = zones$plot_sd_kg_c_m2, CV = zones$cv)
      if (!any(is.finite(values))) next
      label <- paste("Zone area means", switch(kind, Estimate = "estimate", "Point SD" = "SD", CV = "CV"), sep = " · ")
      widget <- leaflet::addPolygons(widget, data = zones[is.finite(values), ], group = label,
        color = "#253746", weight = 1, fillColor = scales[[kind]]$palette(values[is.finite(values)]),
        fillOpacity = 0.88, popup = zone_popup[is.finite(values)])
      layer_kind[label] <- kind
    }
  }
  points <- sf::st_transform(grid$cores, 4326)
  stock_column <- paste0("C Stock 0-", suffix)
  core_popup <- sprintf("<strong>%s</strong><br>Plot: %s<br>Zone: %s<br>Stock: %s kg C/m²",
    htmltools::htmlEscape(points[["Core ID"]]), htmltools::htmlEscape(points[["Plot ID"]]),
    htmltools::htmlEscape(points$zone_assigned), display_number(points[[stock_column]]))
  widget <- leaflet::addCircleMarkers(widget, data = points, group = "Core locations",
    radius = 4, color = "#253746", weight = 1, fillColor = "white", fillOpacity = 1, popup = core_popup)

  surface_names <- names(layer_kind)
  if (!is.null(embedding_map)) {
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
    leaflet::hideGroup(setdiff(surface_names, default_surface)) |>
    leaflet::addLayersControl(overlayGroups = c("AOI boundary", "Core locations"),
                               options = leaflet::layersControlOptions(collapsed = TRUE))
  legend_ids <- stats::setNames(paste0("stock-legend-", seq_along(scales)), names(scales))
  for (kind in names(scales)) {
    scale <- scales[[kind]]
    widget <- leaflet::addLegend(widget, position = "bottomright", pal = scale$palette,
      values = scale$domain, title = titles[[kind]], opacity = 1,
      labFormat = leaflet::labelFormat(digits = 2),
      className = paste("info legend", legend_ids[[kind]]))
  }

  tags <- htmltools::tags
  widget <- htmlwidgets::prependContent(widget, tags$div(class = "stock-map-controls",
    style = "padding: .4rem 0 .65rem; font: 15px system-ui, sans-serif;",
    tags$label("Map layer ",
      tags$select(class = "stock-map-layer", lapply(surface_names, function(label) {
        tags$option(label, value = label, selected = if (label == default_surface) "selected")
      }))
    ),
    tags$p(style = "margin: .3rem 0 0;", paste0(source_label, ". Click the map to read the selected value."))
  ))
  widget <- htmlwidgets::onRender(widget, '
    function(el, x, data) {
      const binding = HTMLWidgets.find("#" + el.id);
      const map = binding.getMap();
      const menu = el.previousElementSibling.querySelector(".stock-map-layer");
      const selectLayer = () => {
        if (!map.layerManager) return;
        const index = data.groups.indexOf(menu.value);
        data.groups.forEach((name, i) => {
          const group = map.layerManager.getLayerGroup(name, false);
          if (i === index && !map.hasLayer(group)) map.addLayer(group);
          if (i !== index && map.hasLayer(group)) map.removeLayer(group);
        });
        el.querySelectorAll(".legend").forEach(legend => {
          legend.style.display = data.legends[index] && legend.classList.contains(data.legends[index]) ? "" : "none";
        });
        map.closePopup();
      };
      menu.addEventListener("change", selectLayer);
      map.on("click", async event => {
        const selected = menu.value;
        const raster = data.values[selected];
        if (!raster) return;
        const point = map.options.crs.project(event.latlng);
        const [xmin, xmax, ymin, ymax] = raster.extent;
        const column = Math.floor((point.x - xmin) / (xmax - xmin) * raster.ncol);
        const row = Math.floor((ymax - point.y) / (ymax - ymin) * raster.nrow);
        if (column < 0 || column >= raster.ncol || row < 0 || row >= raster.nrow) return;
        if (!raster.bytes) {
          const bytes = Uint8Array.from(atob(raster.data), c => c.charCodeAt(0));
          const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate"));
          raster.bytes = new Response(stream).arrayBuffer().then(buffer => new DataView(buffer));
        }
        const bytes = await raster.bytes;
        const value = bytes.getFloat32(4 * (row * raster.ncol + column), true);
        if (!Number.isFinite(value) || menu.value !== selected) return;
        const content = document.createElement("span");
        content.textContent = selected + ": " + value.toFixed(2) +
          (data.units[data.groups.indexOf(selected)] || "");
        map.openPopup(content, event.latlng);
      });
      // Leaflet defers drawing a map in a hidden depth panel until it has a size.
      let initialized = false;
      new ResizeObserver(() => {
        if (!el.clientWidth || !el.clientHeight) return;
        binding.resize(el.clientWidth, el.clientHeight);
        if (!initialized) {
          selectLayer();
          initialized = true;
        }
      }).observe(el);
    }', data = list(groups = surface_names, legends = unname(legend_ids[layer_kind[surface_names]]),
                   units = ifelse(layer_kind[surface_names] == "CV", " (SD / mean)", " kg C/m²"),
                   values = click_values))

  cv_rasters <- layers[names(layer_kind)[layer_kind == "CV"]]
  cv_rasters <- Filter(Negate(is.null), cv_rasters)
  if (!is.null(zones) && any(is.finite(zones$cv))) {
    zone_shapes <- terra::project(terra::vect(zones), terra::crs(grid$template))
    zone_cv <- terra::rasterize(zone_shapes, grid$template, field = "cv")
    cv_rasters[["Zone area means · CV"]] <- zone_cv
  }
  list(widget = widget, surfaces = surface_names, default_surface = default_surface,
       stats = dplyr::bind_rows(surface_stats), cv_rasters = cv_rasters)
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

# Open a Leaflet map in the RStudio Viewer. Both map sets include a layer menu.
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
    area <- .read_view_output(local_folder, "area_summary.csv")
    by_zone <- dplyr::filter(area, method == "By zone")
    pooled <- dplyr::filter(area, method == "Pooled")
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
                                transfer_folder, by_zone, prior, bayesian, layer = layer, pooled = pooled)
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
