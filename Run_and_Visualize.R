source("script_1.R")
source("Script_2_local_mapping.R")
source("Script_2.2_Bayesian.R")
source("Script_3_covariate_transfer.R")
source("Script_4_reporting.R")

source("helpers/view_outputs.R")

# Run any of these lines in the R console to explore the results.
# view_map(30)
# view_map(30, "transfer")
# view_map(30, "transfer", layer = "Local Matérn correction")
# view_plot("totals")
# view_plot("cores")
# view_plot("transfer")
# View(view_table("transfer", depth_cm = 30))
# Make report:
quarto::quarto_render("Blue_Carbon_Report.qmd")
