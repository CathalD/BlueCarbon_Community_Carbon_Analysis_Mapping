# Part 1. Summarize sediment cores and compare with Pacific Northwest blue carbon.
# Run from the workshop folder. Install these packages once if needed:
# install.packages(c("readxl", "dplyr", "readr", "tidyr", "ggplot2"))

library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

workbook <- "data/field/blue_carbon_field_cores.xlsm"
output_folder <- "outputs/script1"
depths_cm <- c(30, 50, 100)
pnw <- c("AK", "BC", "WA", "OR")
habitats <- c(EM = "Salt marsh", SG = "Seagrass", FL = "Tideflat", TS = "Tidal swamp")
tol <- 1e-8

# 1. Read the workbook inputs. Row 4 contains the column names.
read_tab <- function(sheet) read_excel(workbook, sheet = sheet, skip = 3, col_types = "text")
number <- function(x) suppressWarnings(as.numeric(x))
valid_pct <- function(x) is.finite(x) & x >= 0 & x <= 100

settings <- read_tab("7. Settings")
setting <- function(name) number(settings$Value[match(name, settings$Parameter)])
default_diameter <- setting("CORER_DIAMETER_CM")
loi_intercept <- setting("LOI_OC_INTERCEPT")
loi_slope <- setting("LOI_OC_SLOPE")
workbook_note <- read_excel(workbook, sheet = "0. Instructions", range = "B2",
                            col_names = FALSE, .name_repair = "minimal")[[1]]
simulated <- any(grepl("simulat", workbook_note, ignore.case = TRUE))
field_label <- if (simulated) "Simulated cores" else "Field cores"

plots <- read_tab("1. Plot & Site Log") |>
  filter(!is.na(`Plot ID`)) |>
  transmute(
    `Plot ID`, `Site ID`, plot_latitude = number(Latitude), plot_longitude = number(Longitude),
    plot_zone = Location, `Observed habitat`, `Wetland type`
  )

cores <- read_tab("2. Core Log") |>
  filter(!is.na(`Core ID`)) |>
  transmute(
    `Core ID`, `Plot ID`, Latitude = number(Latitude), Longitude = number(Longitude),
    zone = `Sampling zone`, inserted_cm = number(`Inserted length (cm)`),
    recovered_cm = number(`Recovered length (cm)`),
    compaction_factor = case_when(
      is.finite(inserted_cm) & is.finite(recovered_cm) &
        recovered_cm > 0 & inserted_cm >= recovered_cm ~ inserted_cm / recovered_cm,
      is.na(`Inserted length (cm)`) & is.na(`Recovered length (cm)`) &
        tolower(`Compaction assumption`) == "assume none" ~ 1
    ),
    diameter_cm = if_else(is.na(`Internal diameter override (cm)`),
                          default_diameter, number(`Internal diameter override (cm)`))
  )

if (anyDuplicated(plots$`Plot ID`) || anyDuplicated(cores$`Core ID`)) {
  stop("Plot IDs and Core IDs must be unique in their logs.")
}

# Match explicit habitat labels or Janousek codes. Unknown habitats stay unclassified.
habitat_codes <- c(
  em = "EM", "salt marsh" = "EM", saltmarsh = "EM", "tidal marsh" = "EM",
  "emergent marsh" = "EM", sg = "SG", seagrass = "SG", eelgrass = "SG",
  "seagrass meadow" = "SG", "eelgrass meadow" = "SG",
  fl = "FL", tideflat = "FL", "tidal flat" = "FL", mudflat = "FL", "mud flat" = "FL",
  ts = "TS", "tidal swamp" = "TS", "tidal forest" = "TS", "tidal forested wetland" = "TS"
)
cores <- cores |>
  left_join(plots, by = "Plot ID") |>
  mutate(
    Latitude = coalesce(Latitude, plot_latitude), Longitude = coalesce(Longitude, plot_longitude),
    zone = coalesce(zone, plot_zone, "Unspecified"),
    ecosystem = unname(habitat_codes[tolower(coalesce(`Observed habitat`, `Wetland type`))]),
    habitat = coalesce(unname(habitats[ecosystem]), "Unclassified")
  ) |>
  select(-plot_latitude, -plot_longitude, -plot_zone)

# 2. Convert recovered depths to in-situ depths. The correction assumes full
# recovery and uniform compression. Dry mass and carbon mass are conserved.
sections <- read_tab("3. Sediment Data") |>
  transmute(
    workbook_row = row_number() + 4L, `Plot ID`, `Core ID`, `Section ID`,
    recovered_top = `Recovered top (cm)`, recovered_bottom = `Recovered bottom (cm)`,
    dry_mass = `Whole-slice dry mass (g)`, direct_oc = `Organic carbon (%)`,
    total_c = `Total carbon (%)`, inorganic_c = `Inorganic carbon (%)`, loi = `LOI₅₅₀ (%)`
  ) |>
  filter(if_any(-workbook_row, ~ !is.na(.x))) |>
  left_join(
    cores |> select(`Core ID`, core_plot = `Plot ID`, recovered_cm, compaction_factor, diameter_cm),
    by = "Core ID"
  ) |>
  mutate(
    ids_ok = !is.na(core_plot) & !is.na(`Plot ID`) & `Plot ID` == core_plot &
      `Plot ID` %in% plots$`Plot ID`,
    measured_oc = number(direct_oc), tc = number(total_c), ic = number(inorganic_c),
    loi_pct = number(loi), calibrated_oc = loi_intercept + loi_slope * loi_pct,
    carbon_method = case_when(
      valid_pct(measured_oc) ~ "Direct organic carbon",
      is.na(direct_oc) & valid_pct(tc) & valid_pct(ic) & tc >= ic ~ "Total minus inorganic carbon",
      is.na(direct_oc) & is.na(total_c) & is.na(inorganic_c) &
        valid_pct(loi_pct) & valid_pct(calibrated_oc) ~ "Locally calibrated LOI"
    ),
    oc_pct = case_when(
      carbon_method == "Direct organic carbon" ~ measured_oc,
      carbon_method == "Total minus inorganic carbon" ~ tc - ic,
      carbon_method == "Locally calibrated LOI" ~ calibrated_oc
    ),
    recovered_top_cm = number(recovered_top), recovered_bottom_cm = number(recovered_bottom),
    recovered_thickness_cm = recovered_bottom_cm - recovered_top_cm,
    dry_mass_g = number(dry_mass),
    geometry_ok = is.finite(compaction_factor) & is.finite(recovered_top_cm) &
      is.finite(recovered_bottom_cm) & recovered_top_cm >= 0 & recovered_thickness_cm > 0 &
      (is.na(recovered_cm) | recovered_bottom_cm <= recovered_cm + tol),
    top_cm = if_else(geometry_ok, recovered_top_cm * compaction_factor, NA_real_),
    bottom_cm = if_else(geometry_ok, recovered_bottom_cm * compaction_factor, NA_real_),
    thickness_cm = bottom_cm - top_cm, midpoint_cm = (top_cm + bottom_cm) / 2,
    area_cm2 = pi * (diameter_cm / 2)^2,
    mass_ok = is.finite(dry_mass_g) & dry_mass_g > 0 & is.finite(diameter_cm) & diameter_cm > 0,
    recovered_bd = if_else(geometry_ok & mass_ok,
                           dry_mass_g / (area_cm2 * recovered_thickness_cm), NA_real_),
    in_situ_bd = recovered_bd / compaction_factor,
    carbon_density = 10 * in_situ_bd * oc_pct,  # mg C/cm³
    section_stock = if_else(geometry_ok & mass_ok,
                            10 * (dry_mass_g * oc_pct / 100) / area_cm2, NA_real_)  # kg C/m²
  ) |>
  arrange(`Core ID`, recovered_top_cm, recovered_bottom_cm) |>
  group_by(`Core ID`) |>
  mutate(
    duplicate_section = duplicated(`Section ID`) | duplicated(`Section ID`, fromLast = TRUE),
    previous_bottom = lag(cummax(if_else(geometry_ok, bottom_cm, -Inf)), default = -Inf),
    next_top = lead(recovered_top_cm * compaction_factor, default = Inf),
    overlap = geometry_ok & (top_cm < previous_bottom - tol |
      coalesce(bottom_cm > next_top + tol, FALSE)),
    touches_previous = coalesce(abs(top_cm - lag(bottom_cm, default = 0)) <= tol, FALSE),
    unplaced_section = any(!is.finite(recovered_top_cm)),
    in_stock_profile = cumall(ids_ok & !is.na(`Section ID`) & !duplicate_section & !overlap &
      touches_previous & is.finite(section_stock) & !unplaced_section),
    cumulative_stock = if_else(in_stock_profile, cumsum(section_stock), NA_real_),
    section_check = case_when(
      !ids_ok ~ "Core or plot ID does not match the logs",
      is.na(`Section ID`) ~ "Missing Section ID",
      duplicate_section ~ "Duplicate Section ID within core",
      !geometry_ok ~ "Check recovered depths or compaction inputs",
      overlap ~ "Overlapping slices",
      !mass_ok ~ "Check dry mass or internal diameter",
      is.na(oc_pct) ~ "No valid organic carbon result",
      unplaced_section ~ "A slice in this core has no usable starting depth",
      !touches_previous ~ "Gap from surface or previous slice",
      !in_stock_profile ~ "Below a break in the stock profile",
      TRUE ~ "OK"
    ),
    lab_review = case_when(
      recovered_bd < setting("QC_BD_MIN") | recovered_bd > setting("QC_BD_MAX") ~ "Review recovered bulk density",
      oc_pct > setting("QC_CARBON_PCT_MAX") ~ "Review organic carbon concentration"
    )
  ) |>
  ungroup()

if (!nrow(cores) || !nrow(sections)) stop("Enter core and sediment rows before running the analysis.")

# 3. Sum stock to the target depth or the sampled bottom, whichever comes first.
# A break in the profile still prevents estimating stock below that break.
stock_to_depth <- function(depth, top, bottom, stock, reaches_sampled_bottom) {
  if (!length(bottom)) return(NA_real_)
  sampled_bottom <- max(bottom)
  if (!reaches_sampled_bottom && sampled_bottom < depth - tol) return(NA_real_)
  end <- min(depth, sampled_bottom)
  overlap_cm <- pmax(0, pmin(bottom, end) - top)
  carbon_per_cm <- stock / (bottom - top)
  sum(carbon_per_cm * overlap_cm)
}

stock_columns <- paste0("C Stock 0-", depths_cm, "cm")
profile <- sections |>
  group_by(`Core ID`) |>
  mutate(reaches_sampled_bottom = all(in_stock_profile %in% TRUE)) |>
  ungroup() |>
  filter(in_stock_profile)
core_depths <- profile |>
  group_by(`Core ID`) |>
  summarise(
    `Stock coverage (cm)` = max(bottom_cm),
    `Total Carbon Stock` = sum(section_stock),
    `Avg Bulk Density` = weighted.mean(in_situ_bd, thickness_cm),
    `Avg Carbon (%)` = weighted.mean(oc_pct, thickness_cm),
    .groups = "drop"
  )
for (i in seq_along(depths_cm)) {
  stocks <- profile |>
    group_by(`Core ID`) |>
    summarise(stock = stock_to_depth(depths_cm[i], top_cm, bottom_cm, section_stock,
                                    first(reaches_sampled_bottom)),
              .groups = "drop") |>
    rename(!!stock_columns[i] := stock)
  core_depths <- left_join(core_depths, stocks, by = "Core ID")
}

section_counts <- sections |>
  filter(!is.na(`Core ID`)) |>
  group_by(`Core ID`) |>
  summarise(
    `Total Sections` = n(),
    `Max Depth (cm)` = if (any(is.finite(bottom_cm))) max(bottom_cm, na.rm = TRUE) else NA_real_,
    `Sections flagged` = sum(section_check != "OK" | !is.na(lab_review)),
    .groups = "drop"
  )

core_summary <- cores |>
  left_join(section_counts, by = "Core ID") |>
  left_join(core_depths, by = "Core ID") |>
  mutate(
    `Total Sections` = coalesce(`Total Sections`, 0L),
    `Sections flagged` = coalesce(`Sections flagged`, 0L),
    data_source = field_label,
    core_check = case_when(
      is.na(`Site ID`) ~ "Check Plot ID or Site ID",
      !is.finite(compaction_factor) ~ "Check compaction inputs",
      compaction_factor > setting("QC_COMPACTION_FACTOR_MAX") ~ "Review high compaction factor",
      !is.finite(diameter_cm) | diameter_cm <= 0 ~ "Check internal diameter",
      `Total Sections` == 0 ~ "No sediment sections",
      `Sections flagged` > 0 ~ "Review sediment sections",
      TRUE ~ "OK"
    )
  ) |>
  relocate(`Core ID`, `Plot ID`, Latitude, Longitude, zone, habitat,
           `Total Sections`, `Max Depth (cm)`, `Stock coverage (cm)`, all_of(stock_columns),
           `Total Carbon Stock`, core_check)

field_stocks <- core_summary |>
  select(`Core ID`, zone, ecosystem, habitat, data_source,
         stock_coverage_cm = `Stock coverage (cm)`, all_of(stock_columns)) |>
  pivot_longer(all_of(stock_columns), names_to = "depth_cm",
               names_pattern = "C Stock 0-(\\d+)cm", names_transform = list(depth_cm = as.integer),
               values_to = "stock_kg_m2") |>
  mutate(short_core = is.finite(stock_kg_m2) & stock_coverage_cm < depth_cm - tol)

depth_summary <- bind_rows(mutate(field_stocks, zone = "All cores"), field_stocks) |>
  group_by(zone, depth_cm, data_source) |>
  summarise(
    cores = n(), complete = sum(is.finite(stock_kg_m2)), short_cores = sum(short_core),
    mean_kg_m2 = if (complete > 0) mean(stock_kg_m2, na.rm = TRUE) else NA_real_,
    sd_kg_m2 = sd(stock_kg_m2, na.rm = TRUE),
    median_kg_m2 = median(stock_kg_m2, na.rm = TRUE),
    min_kg_m2 = if (complete > 0) min(stock_kg_m2, na.rm = TRUE) else NA_real_,
    max_kg_m2 = if (complete > 0) max(stock_kg_m2, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# 4. Janousek et al. (2025), doi:10.1029/2024GB008239.
# Compare published stocks only where recorded sampling reaches the target depth.
# Published stocks in Mg C/ha are divided by 10 to obtain kg C/m².
locations <- read_csv("data/reference/janousek_core_locations.csv.gz",
                      col_types = cols(.default = col_character())) |>
  filter(state %in% pnw, ecosystem %in% names(habitats)) |>
  mutate(sampled_cm = number(sub("^[>~]", "", core_depth)))
if (anyDuplicated(locations$core_id)) stop("Check duplicate Janousek Core IDs before joining layers.")

reference_stocks <- locations |>
  select(core_id, state, ecosystem, estuary_id, sampled_cm, starts_with("carbon_stock_")) |>
  pivot_longer(starts_with("carbon_stock_"), names_to = "depth_cm",
               names_pattern = "carbon_stock_(\\d+)cm", names_transform = list(depth_cm = as.integer),
               values_transform = list(value = number), values_to = "value") |>
  filter(sampled_cm >= depth_cm, is.finite(value), value >= 0) |>
  mutate(stock_kg_m2 = value / 10, habitat = unname(habitats[ecosystem]))

reference_summary <- reference_stocks |>
  group_by(ecosystem, depth_cm) |>
  summarise(
    habitat = first(habitat), ref_cores = n(), ref_estuaries = n_distinct(estuary_id, na.rm = TRUE),
    ref_q25 = quantile(stock_kg_m2, 0.25), ref_median = median(stock_kg_m2),
    ref_q75 = quantile(stock_kg_m2, 0.75), reference_values = list(stock_kg_m2),
    .groups = "drop"
  )

comparison <- field_stocks |>
  left_join(select(reference_summary, -habitat), by = c("ecosystem", "depth_cm")) |>
  rowwise() |>
  mutate(
    percentile = if (is.finite(stock_kg_m2) && length(reference_values)) {
      100 * mean(reference_values <= stock_kg_m2)
    } else NA_real_,
    comparison = case_when(
      !is.finite(stock_kg_m2) ~ "Incomplete field coverage",
      is.na(ecosystem) ~ "Habitat unclassified",
      is.na(ref_cores) ~ "No Janousek cores with matching habitat and depth",
      short_core ~ "Shorter study core compared with target-depth references",
      TRUE ~ "Matched by habitat and depth"
    )
  ) |>
  ungroup() |>
  select(-reference_values)

# 5. Plots and tables. Each line or point represents a core unless labeled otherwise.
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
write_csv(core_summary, file.path(output_folder, "core_summary.csv"))
write_csv(depth_summary, file.path(output_folder, "field_stock_summary.csv"))
write_csv(
  sections |> select(workbook_row, `Core ID`, `Plot ID`, `Section ID`, recovered_top_cm,
                     recovered_bottom_cm, top_cm, bottom_cm, dry_mass_g, diameter_cm,
                     compaction_factor, recovered_bd, in_situ_bd, oc_pct, carbon_method,
                     carbon_density, section_stock, cumulative_stock, section_check, lab_review),
  file.path(output_folder, "sediment_sections.csv")
)
write_csv(comparison, file.path(output_folder, "janousek_comparison.csv"))
write_csv(select(reference_summary, -reference_values),
          file.path(output_folder, "janousek_reference_summary.csv"))
write_csv(
  core_summary |> filter(is.finite(Latitude), is.finite(Longitude),
                         abs(Latitude) <= 90, abs(Longitude) <= 180, !is.na(`Plot ID`)),
  file.path(output_folder, "mapping_data.csv")
)

metrics <- c(oc_pct = "Organic carbon (%)", in_situ_bd = "Dry bulk density (g/cm³)",
             carbon_density = "Carbon density (mg C/cm³)")
field_profiles <- sections |>
  left_join(select(cores, `Core ID`, zone), by = "Core ID") |>
  group_by(`Core ID`) |>
  mutate(profile_run = cumsum(!touches_previous)) |>
  ungroup() |>
  pivot_longer(all_of(names(metrics)), names_to = "metric", values_to = "value") |>
  mutate(value = if_else(ids_ok & !duplicate_section & !overlap, value, NA_real_),
         metric = factor(metric, levels = names(metrics), labels = metrics))

field_plot <- ggplot(field_profiles, aes(value, midpoint_cm, colour = zone)) +
  geom_path(aes(group = interaction(`Core ID`, profile_run)), alpha = 0.65, na.rm = TRUE) +
  geom_point(size = 1.2, alpha = 0.8, na.rm = TRUE) +
  scale_y_reverse() +
  facet_wrap(~ metric, scales = "free_x") +
  labs(title = "Sediment by depth", subtitle = field_label,
       caption = "Depth and dry bulk density account for core compression.",
       x = NULL, y = "Depth (cm)", colour = "Zone") +
  theme_bw() +
  theme(panel.spacing.x = grid::unit(0.55, "cm"))

cumulative <- bind_rows(
  profile |> distinct(`Core ID`) |> mutate(bottom_cm = 0, cumulative_stock = 0),
  profile |> select(`Core ID`, bottom_cm, cumulative_stock)
) |>
  left_join(select(cores, `Core ID`, zone), by = "Core ID")
cumulative_plot <- ggplot(cumulative, aes(cumulative_stock, bottom_cm, group = `Core ID`, colour = zone)) +
  geom_hline(yintercept = depths_cm, linetype = "dotted", colour = "grey65") +
  geom_path() + geom_point(size = 1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Cumulative carbon stock", subtitle = field_label,
       x = "Carbon stock (kg C/m²)", y = "Depth (cm)", colour = "Zone") +
  theme_bw()

depth_label <- function(x) factor(paste0("0–", x, " cm"), levels = paste0("0–", depths_cm, " cm"))
zone_stocks <- depth_summary |>
  filter(zone != "All cores", complete > 0) |>
  mutate(depth = depth_label(depth_cm), lower = mean_kg_m2 - sd_kg_m2,
         upper = mean_kg_m2 + sd_kg_m2)
zone_sd <- zone_stocks |>
  filter(is.finite(sd_kg_m2)) |>
  pivot_longer(c(lower, upper), values_to = "bound")

stock_context_plot <- function(reference, title, group_label) {
  reference <- mutate(reference, depth = depth_label(depth_cm))
  counts <- count(reference, depth, reference_group)

  ggplot(reference, aes(stock_kg_m2, reference_group)) +
    geom_boxplot(fill = "grey93", outlier.shape = NA, width = 0.58) +
    geom_point(position = position_jitter(width = 0, height = 0.12, seed = 42),
               colour = "grey30", alpha = 0.28, size = 0.9) +
    geom_blank(data = zone_sd, aes(x = bound), inherit.aes = FALSE) +
    geom_vline(data = zone_sd, aes(xintercept = bound, colour = zone),
               linewidth = 0.65, linetype = "dashed", alpha = 0.45, show.legend = FALSE) +
    geom_vline(data = zone_stocks, aes(xintercept = mean_kg_m2, colour = zone), linewidth = 0.8) +
    geom_text(data = counts, aes(x = Inf, label = paste0(n, if_else(n == 1L, " core", " cores"))),
              hjust = 1.12, size = 3) +
    facet_wrap(~ depth, scales = "free_x", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.32))) +
    labs(
      title = title,
      subtitle = field_label,
      caption = paste0("Grey points show Janousek cores. Solid lines show zone averages.\n",
                       "Dashed lines show one standard deviation above and below the average."),
      x = "Carbon stock (kg C/m²)", y = group_label, colour = "Zone"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

stock_plot <- stock_context_plot(
  mutate(reference_stocks, reference_group = habitat),
  "Carbon stock by ecosystem",
  "Ecosystem"
)

# The BC estuary rows are subsets of the BC and Pacific Northwest marsh rows.
marsh_stocks <- filter(reference_stocks, ecosystem == "EM")
bc_marsh <- filter(marsh_stocks, state == "BC")
bc_estuaries <- bc_marsh |>
  filter(!is.na(estuary_id)) |>
  mutate(reference_group = paste("BC estuary", estuary_id))
marsh_groups <- c(sort(unique(bc_estuaries$reference_group)), "British Columbia", "Pacific Northwest")
marsh_reference <- bind_rows(
  bc_estuaries,
  mutate(bc_marsh, reference_group = "British Columbia"),
  mutate(marsh_stocks, reference_group = "Pacific Northwest")
) |>
  mutate(reference_group = factor(reference_group, levels = rev(marsh_groups)))
marsh_plot <- stock_context_plot(marsh_reference, "Salt marsh carbon stock by estuary and region",
                                 "Estuary or region")
marsh_summary <- marsh_reference |>
  group_by(reference_group, depth_cm) |>
  summarise(ref_cores = n(), ref_median = median(stock_kg_m2),
            ref_q25 = quantile(stock_kg_m2, 0.25), ref_q75 = quantile(stock_kg_m2, 0.75),
            .groups = "drop")
write_csv(marsh_summary, file.path(output_folder, "janousek_bc_marsh_summary.csv"))

# Give each core one value per 10-cm depth bin before calculating group profiles.
bin_midpoint <- function(depth) 10 * floor(depth / 10) + 5

zone_profiles <- field_profiles |>
  filter(is.finite(value), is.finite(midpoint_cm)) |>
  mutate(bin_cm = bin_midpoint(midpoint_cm)) |>
  group_by(zone, metric, bin_cm, `Core ID`) |>
  summarise(value = mean(value), .groups = "drop_last") |>
  summarise(value = mean(value), .groups = "drop_last") |>
  arrange(bin_cm, .by_group = TRUE) |>
  mutate(profile_run = cumsum(bin_cm - lag(bin_cm, default = first(bin_cm)) > 10)) |>
  ungroup()

reference_layers <- read_csv(
  "data/reference/janousek_layers_v1.csv.gz",
  col_types = cols_only(profile_id = col_character(), upper_depth = col_double(),
                        lower_depth = col_double(), BDOD = col_double(), OrgC_pct = col_double(),
                        source_interval_reused = col_integer())
) |>
  inner_join(select(locations, core_id, state, ecosystem, sampled_cm), by = c("profile_id" = "core_id")) |>
  filter(source_interval_reused == 0, upper_depth >= 0, lower_depth > upper_depth,
         lower_depth <= sampled_cm, (upper_depth + lower_depth) / 2 < max(depths_cm)) |>
  transmute(
    profile_id, state, ecosystem, bin_cm = bin_midpoint((upper_depth + lower_depth) / 2),
    oc_pct = if_else(valid_pct(OrgC_pct), OrgC_pct, NA_real_),
    in_situ_bd = if_else(is.finite(BDOD) & BDOD > 0, BDOD, NA_real_),
    carbon_density = 10 * in_situ_bd * oc_pct
  ) |>
  pivot_longer(all_of(names(metrics)), names_to = "metric", values_to = "value") |>
  filter(is.finite(value)) |>
  group_by(profile_id, state, ecosystem, bin_cm, metric) |>
  summarise(value = median(value), .groups = "drop")

reference_profiles <- bind_rows(
  mutate(reference_layers, reference_group = "Pacific Northwest"),
  reference_layers |> filter(state == "BC", ecosystem == "EM") |>
    mutate(reference_group = "BC salt marsh")
) |>
  group_by(reference_group, metric, bin_cm) |>
  summarise(cores = n(), q10 = quantile(value, 0.10), q25 = quantile(value, 0.25),
            median = median(value), q75 = quantile(value, 0.75), q90 = quantile(value, 0.90),
            .groups = "drop_last") |>
  filter(cores >= 10) |>
  arrange(bin_cm, .by_group = TRUE) |>
  mutate(profile_run = cumsum(bin_cm - lag(bin_cm, default = first(bin_cm)) > 10)) |>
  ungroup() |>
  mutate(metric = factor(metric, levels = names(metrics), labels = metrics),
         reference_group = factor(reference_group, levels = c("Pacific Northwest", "BC salt marsh")))

reference_plot <- ggplot(reference_profiles, aes(median, bin_cm, group = profile_run)) +
  geom_ribbon(aes(xmin = q10, xmax = q90), orientation = "y", fill = "grey89") +
  geom_ribbon(aes(xmin = q25, xmax = q75), orientation = "y", fill = "grey72") +
  geom_path(colour = "grey35", linewidth = 0.65) +
  geom_point(colour = "grey35", size = 0.85) +
  geom_path(data = zone_profiles, aes(value, bin_cm, colour = zone,
                                     group = interaction(zone, profile_run)), linewidth = 0.65) +
  geom_point(data = zone_profiles, aes(value, bin_cm, colour = zone), size = 1.5) +
  scale_y_reverse() +
  facet_grid(reference_group ~ metric, scales = "free_x") +
  labs(title = "Sediment by depth compared with Janousek cores",
       subtitle = paste(field_label, "averaged by zone every 10 cm."),
       caption = paste0("The grey line shows the Janousek median. Dark shading shows the middle 50% of values.\n",
                        "Light shading shows the middle 80%. Each range uses at least 10 cores."),
       x = NULL, y = "Depth (cm)", colour = "Zone") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_folder, "sediment_profiles.png"), field_plot, width = 11, height = 5, dpi = 200)
ggsave(file.path(output_folder, "cumulative_carbon_stock.png"), cumulative_plot, width = 7, height = 5, dpi = 200)
ggsave(file.path(output_folder, "carbon_stock_context.png"), stock_plot, width = 11, height = 5, dpi = 200)
ggsave(file.path(output_folder, "janousek_profiles.png"), reference_plot, width = 11, height = 7, dpi = 200)
ggsave(file.path(output_folder, "janousek_bc_marsh_stocks.png"), marsh_plot, width = 11, height = 6, dpi = 200)

print(depth_summary, n = Inf)
message(nrow(cores), " cores and ", nrow(sections), " sections. Results in ", output_folder, ".")
