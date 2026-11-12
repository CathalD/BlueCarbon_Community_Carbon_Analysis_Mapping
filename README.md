# Community Blue Carbon Analysis and Mapping

This workshop uses simulated sediment cores to compare methods for mapping blue carbon stocks and uncertainty.

## Download the workshop from R

Copy and paste this block into the R or RStudio console.

```r
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis", repos = "https://cloud.r-project.org")
}
usethis::use_course("CathalD/BlueCarbon_Community_Carbon_Analysis_Mapping")
```

Confirm the download location when prompted. Once the workshop folder opens, follow **Run the workshop** below.

## Run the workshop

Open R or RStudio with this folder as the working directory. Install Quarto and the R packages once, then run the scripts in order through the workshop runner.

```r
install.packages(c("readxl", "dplyr", "readr", "tidyr", "ggplot2", "sf", "terra",
                   "mapview", "leaflet", "htmlwidgets", "htmltools", "knitr",
                   "rmarkdown", "base64enc", "quarto"))
source("Run_and_Visualize.R")
```

The runner reads the supplied workshop inputs in `data/`, writes results to `outputs/` and builds the HTML report. The AlphaEarth raster is supplied as lossless GeoTIFF tiles. Keep the tiles together in `data/covariate_transfer/`. Script 3 reads them as a single raster. See the [Earth Engine workflow](data/covariate_transfer/README.md) to prepare inputs for another site.

Local stock is summed to the target depth or the sampled bottom, whichever comes first. When a section crosses the target, the calculation uses its carbon per centimetre for the part inside the target interval. A valid 20 cm core contributes the same measured stock at the 30, 50 and 100 cm targets.

## View results in R

```r
view_map(30)                                             # Local and Bayesian maps
view_map(30, "transfer")                                 # Transfer map with a layer menu
view_map(30, "transfer", layer = "Local Matérn correction")
view_plot("totals")                                      # Stock by method and depth
view_plot("cores")                                       # Stock by core
view_plot("transfer")                                    # Transfer stock estimates
View(view_table("transfer", depth_cm = 30))               # Transfer summary in RStudio
view_table("totals")                                     # Mapping totals in the console
```

Map depths are 30, 50 and 100 cm. The transfer map also offers the AlphaEarth embeddings, Janousek prior, updated estimate and both uncertainty layers. To load the viewing functions after restarting R with existing outputs, run `source("helpers/view_outputs.R")`.

## Build the report

To rebuild [the report](Blue_Carbon_Report.qmd) from existing results, run Quarto from the project folder.

```sh
quarto render Blue_Carbon_Report.qmd
```

The rendered report, prediction tables, figures and saved maps are generated locally. The repository contains the source files and workshop inputs needed to recreate them.

## Data sources

The reference cores come from Janousek et al. (2025). [Reference dataset](https://doi.org/10.25573/SERC.28127486) · [Paper](https://doi.org/10.1029/2024GB008239)

The AlphaEarth Foundations Satellite Embedding dataset is produced by Google and Google DeepMind. The included 2024 exports are clipped and tiled for the workshop and are licensed under CC BY 4.0. [Dataset and attribution](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL)
