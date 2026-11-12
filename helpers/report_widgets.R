# Switch a set of already embedded HTML panels without needing a Shiny server.
report_map_picker <- function(index, folder, title) {
  tags <- htmltools::tags
  if (is.null(index) || !nrow(index)) return(tags$p(paste("No", title, "maps were exported for this run.")))
  if (!"widget_file" %in% names(index)) {
    return(tags$p("Rerun Script_4_reporting.R to embed the updated maps in this report."))
  }
  paths <- file.path(folder, index$widget_file)
  if (any(!file.exists(paths))) stop("Some map widgets are missing. Rerun Script_4_reporting.R.")
  labels <- paste0("0–", index$depth_cm, " cm")
  htmltools::tagList(
    tags$div(class = "report-maps",
      tags$label("Map depth ",
        tags$select(lapply(seq_along(labels), function(i) tags$option(labels[i], value = i)))
      ),
      lapply(seq_along(paths), function(i) {
        tags$div("data-map-depth" = i, hidden = if (i != 1L) "hidden",
                 htmltools::as.tags(readRDS(paths[i])))
      })
    ),
    tags$script(htmltools::HTML('
      (() => {
        const container = document.currentScript.previousElementSibling;
        const menu = container.querySelector("select");
        menu.addEventListener("change", () => {
          container.querySelectorAll("[data-map-depth]").forEach(panel => {
            panel.hidden = panel.dataset.mapDepth !== menu.value;
          });
        });
      })();'))
  )
}

report_estuary_picker <- function(regions, estuaries, study, cores, habitat) {
  tags <- htmltools::tags
  fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 2, format = "f"), "—")
  study <- study |> dplyr::filter(complete > 0)
  study <- dplyr::bind_rows(lapply(seq_len(nrow(study)), function(i) {
    depth <- study$depth_cm[i]
    values <- cores[[paste0("C Stock 0-", depth, "cm")]]
    values <- values[is.finite(values)]
    dplyr::tibble(depth_cm = depth, cores = study$complete[i], mean_kg_m2 = study$mean_kg_m2[i],
                 median_kg_m2 = study$median_kg_m2[i], sd_kg_m2 = study$sd_kg_m2[i],
                 q25_kg_m2 = if (length(values)) unname(stats::quantile(values, .25)) else NA_real_,
                 q75_kg_m2 = if (length(values)) unname(stats::quantile(values, .75)) else NA_real_,
                 group = "Study cores", estuary_key = NA_character_)
  }))
  region_rows <- regions |>
    dplyr::filter(.data$habitat == .env$habitat) |>
    dplyr::mutate(group = region, estuary_key = NA_character_)
  estuary_rows <- estuaries |>
    dplyr::filter(.data$habitat == .env$habitat) |>
    dplyr::mutate(group = paste0(estuary_name, " (", estuary_region, ")"),
                  estuary_key = paste(estuary_region, estuary_id, sep = "-"))
  options <- estuary_rows |>
    dplyr::distinct(estuary_key, group) |>
    dplyr::arrange(group)
  data <- dplyr::bind_rows(study, region_rows, estuary_rows) |>
    dplyr::arrange(depth_cm)
  if (!nrow(data)) return(tags$p("No comparable stocks are available."))

  row_chart <- function(data, i, max_x) {
    x <- function(value) paste0(100 * value / max_x, "%")
    estuary <- data$estuary_key[i]
    style <- if (identical(data$group[i], "Study cores")) "comparison-row comparison-study" else "comparison-row"
    tags$div(class = style, "data-estuary" = if (!is.na(estuary)) estuary,
      hidden = if (!is.na(estuary)) "hidden",
      tags$span(data$group[i]),
      tags$div(class = "comparison-track",
        tags$span(class = "comparison-iqr", style = paste0("left:", x(data$q25_kg_m2[i]),
          ";width:", x(data$q75_kg_m2[i] - data$q25_kg_m2[i]))),
        tags$span(class = "comparison-median", style = paste0("left:", x(data$median_kg_m2[i]))),
        tags$span(class = "comparison-median-value", fmt(data$median_kg_m2[i]))
      )
    )
  }
  chart <- lapply(split(data, data$depth_cm), function(rows) {
    maximum <- max(0.1, rows$q75_kg_m2, rows$median_kg_m2, na.rm = TRUE)
    max_x <- max(pretty(c(0, maximum * 1.05)))
    tags$div(class = "comparison-depth",
      tags$strong(paste0("0–", rows$depth_cm[1], " cm")),
      lapply(seq_len(nrow(rows)), function(i) row_chart(rows, i, max_x)),
      tags$div(class = "comparison-axis", tags$span(),
        tags$div(tags$span("0"), tags$span(paste0(format(max_x, trim = TRUE), " kg C/m²"))))
    )
  })

  headings <- c("Depth", "Group", "Cores", "Mean", "SD", "Median", "Middle 50%")
  htmltools::tagList(
    tags$style(htmltools::HTML('
      .estuary-comparison .comparison-depth { margin: 1rem 0; }
      .estuary-comparison .comparison-row, .estuary-comparison .comparison-axis {
        display: grid; grid-template-columns: minmax(8rem, 13rem) minmax(6rem, 1fr);
        gap: .75rem; align-items: center; padding: .15rem 0; font-size: .9rem;
      }
      .estuary-comparison [hidden], .report-maps [hidden] { display: none !important; }
      .comparison-track { position: relative; min-height: 1.6rem; background: #edf3f3; }
      .comparison-iqr, .comparison-median { position: absolute; top: .35rem; height: .9rem; background: #82a9ad; }
      .comparison-median { width: .2rem; background: #154d4c; }
      .comparison-study .comparison-iqr { background: #dda08f; }
      .comparison-study .comparison-median { background: #ab3b26; }
      .comparison-median-value { position: absolute; right: .35rem; top: .12rem; font-weight: 600; color: #142d36; background: #edf3f3; }
      .comparison-axis > div { display: flex; justify-content: space-between; }
      .stock-map-layer { max-width: 100%; }
    ')),
    tags$div(class = "estuary-comparison",
      tags$label("Compare an estuary ",
        tags$select(tags$option("Regional comparisons only", value = ""),
          lapply(seq_len(nrow(options)), function(i) tags$option(options$group[i], value = options$estuary_key[i]))
        )
      ),
      tags$p("Bars show the middle 50% of stocks. The marker and printed value show the median."),
      tags$div("aria-hidden" = "true", chart),
      tags$div(class = "table-responsive",
        tags$table(class = "table table-sm",
          tags$caption(paste0(habitat, " comparisons by depth. Stock statistics are in kg C/m².")),
          tags$thead(tags$tr(lapply(headings, function(name) tags$th(scope = "col", name)))),
          tags$tbody(lapply(seq_len(nrow(data)), function(i) {
            estuary <- data$estuary_key[i]
            tags$tr("data-estuary" = if (!is.na(estuary)) estuary,
              hidden = if (!is.na(estuary)) "hidden",
              tags$td(paste0("0–", data$depth_cm[i], " cm")), tags$th(scope = "row", data$group[i]),
              tags$td(data$cores[i]), tags$td(fmt(data$mean_kg_m2[i])), tags$td(fmt(data$sd_kg_m2[i])),
              tags$td(fmt(data$median_kg_m2[i])),
              tags$td(paste(fmt(data$q25_kg_m2[i]), "to", fmt(data$q75_kg_m2[i])))
            )
          }))
        )
      )
    ),
    tags$script(htmltools::HTML('
      (() => {
        const container = document.currentScript.previousElementSibling;
        const menu = container.querySelector("select");
        menu.addEventListener("change", () => {
          container.querySelectorAll("[data-estuary]").forEach(row => {
            row.hidden = row.dataset.estuary !== menu.value;
          });
        });
      })();'))
  )
}
