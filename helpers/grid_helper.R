# One grid and one mask for the final mapping and Bayesian products.
# Values describe predictions at cell centres, not computed block averages.

create_analysis_grid <- function() {
  terra::projNetwork(FALSE)
  dir.create(grid_folder, recursive = TRUE, showWarnings = FALSE)
  boundary <- terra::project(terra::vect(aoi_file), analysis_crs)
  template <- terra::rast(
    xmin = floor(terra::xmin(boundary) / cell_size_m) * cell_size_m,
    xmax = ceiling(terra::xmax(boundary) / cell_size_m) * cell_size_m,
    ymin = floor(terra::ymin(boundary) / cell_size_m) * cell_size_m,
    ymax = ceiling(terra::ymax(boundary) / cell_size_m) * cell_size_m,
    resolution = cell_size_m, crs = analysis_crs
  )
  terra::values(template) <- 1
  names(template) <- "template"
  analysis_mask <- terra::mask(template, boundary, touches = TRUE)
  names(analysis_mask) <- "aoi_mask"
  cell_ids <- terra::mask(terra::init(template, "cell"), analysis_mask)
  names(cell_ids) <- "grid_id"

  terra::writeRaster(template, grid_file, overwrite = TRUE)
  terra::writeRaster(analysis_mask, mask_file, overwrite = TRUE)
  terra::writeRaster(cell_ids, cell_id_file, datatype = "INT4S", overwrite = TRUE)
  terra::writeVector(boundary, projected_aoi_file, overwrite = TRUE)
  utils::write.csv(data.frame(
    crs = analysis_crs, cell_size_m = cell_size_m,
    origin_x = 0, origin_y = 0,
    mask_rule = "Cells intersecting the AOI; clip cell areas for totals",
    boundary_md5 = unname(tools::md5sum(aoi_file))
  ), file.path(grid_folder, "grid_settings.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(grid_folder, "session_info.txt"))
  cat("\nShared grid:", nrow(template), "rows x",
      ncol(template), "columns;", cell_size_m, "m cells.\n")
  invisible(load_analysis_grid())
}

load_analysis_grid <- function() {
  terra::projNetwork(FALSE)
  template <- terra::rast(grid_file)
  analysis_mask <- terra::rast(mask_file)
  cell_ids <- terra::rast(cell_id_file)
  terra::compareGeom(template, analysis_mask, cell_ids)
  list(template = template, mask = analysis_mask, cell_ids = cell_ids,
       aoi = terra::vect(projected_aoi_file))
}
