# One cached context source, extent and projection for all maps of an AOI.
.context_cache <- new.env(parent = emptyenv())
.context_data <- function() {
  if (!exists("data", .context_cache, inherits = FALSE)) {
    path <- "data/context/north_america.gpkg"
    if (!file.exists(path)) stop("Offline map context is missing: ", path, ". Restore the course data/context folder.")
    .context_cache$data <- list(land = sf::st_read(path, layer="countries", quiet=TRUE),
      lakes = sf::st_read(path, layer="lakes", quiet=TRUE))
  }
  .context_cache$data
}
.context_extent <- function(aoi) {
  b <- sf::st_bbox(aoi); dx <- b["xmax"]-b["xmin"]; dy <- b["ymax"]-b["ymin"]
  pad <- max(dx, dy) * .22
  c(xmin=unname(b["xmin"]-pad), xmax=unname(b["xmax"]+pad),
    ymin=unname(b["ymin"]-pad), ymax=unname(b["ymax"]+pad))
}
.context_base <- function(aoi) {
  z <- .context_data(); b <- .context_extent(aoi)
  # Transform before cropping so local geometry and analytical coordinates agree.
  land <- suppressWarnings(sf::st_crop(sf::st_transform(z$land, sf::st_crs(aoi)), b))
  water <- suppressWarnings(sf::st_crop(sf::st_transform(z$lakes, sf::st_crs(aoi)), b))
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data=land, fill="#f2f0e9", colour="#a6aaa6", linewidth=.25) +
    ggplot2::geom_sf(data=water, fill="#deedf3", colour=NA)
  if (nrow(land)) {
    labels <- suppressWarnings(sf::st_point_on_surface(land))
    p <- p + ggplot2::geom_sf_text(data=labels, ggplot2::aes(label=ADMIN),
      size=2.5, colour="grey40", check_overlap=TRUE)
  }
  p
}
.context_finish <- function(p, aoi) {
  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(file=NULL)
    on.exit(grDevices::dev.off(),add=TRUE)
  }
  b <- .context_extent(aoi); dx <- b["xmax"]-b["xmin"]; dy <- b["ymax"]-b["ymin"]
  # Scale is projected metres; LAEA preserves area, not all distances.
  bar <- 10^floor(log10(dx*.18)); if (bar < dx*.08) bar <- bar*2
  label <- if (bar>=1000) paste(format(bar/1000,trim=TRUE),"km") else paste(format(bar,trim=TRUE),"m")
  x <- b["xmin"]+.05*dx; y <- b["ymin"]+.05*dy
  z <- .context_data()
  loc <- sf::st_transform(z$land, 3978)
  centre <- sf::st_transform(sf::st_centroid(sf::st_geometry(aoi)), 3978)
  inset <- ggplot2::ggplot() + ggplot2::geom_sf(data=loc, fill="grey92", colour="grey70", linewidth=.1) +
    ggplot2::geom_sf(data=centre, colour="#b2182b", size=2) +
    ggplot2::coord_sf(crs=3978, xlim=c(-3000000,3300000), ylim=c(-800000,4400000), expand=FALSE) +
    ggplot2::theme_void() + ggplot2::theme(panel.background=ggplot2::element_rect(fill="white",colour="grey60"))
  p + ggplot2::coord_sf(crs=sf::st_crs(aoi), xlim=b[c("xmin","xmax")],
      ylim=b[c("ymin","ymax")], expand=FALSE, datum=sf::st_crs(4326)) +
    ggplot2::annotate("segment", x=x, xend=x+bar, y=y, yend=y, linewidth=1) +
    ggplot2::annotate("text", x=x+bar/2, y=y+.035*dy, label=label, size=3) +
    ggplot2::annotation_custom(ggplot2::ggplotGrob(inset),
      xmin=b["xmax"]-.25*dx, xmax=b["xmax"]-.01*dx,
      ymin=b["ymin"]+.01*dy, ymax=b["ymin"]+.23*dy) +
    ggplot2::labs(caption="Made with Natural Earth • 1:50m context • approximate scale", x=NULL,y=NULL) +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="#deedf3",colour=NA),
      panel.grid.major=ggplot2::element_line(colour="grey85",linewidth=.2),
      axis.text=ggplot2::element_text(size=7), legend.title=ggplot2::element_text(size=10),
      plot.caption=ggplot2::element_text(size=7), strip.text=ggplot2::element_text(size=11))
}

regional_location_map <- function(aoi) {
  z <- .context_data(); crs <- sf::st_crs(aoi)
  b <- sf::st_bbox(aoi); pad <- max(150000, max(b["xmax"]-b["xmin"],b["ymax"]-b["ymin"])*.3)
  region <- sf::st_as_sfc(sf::st_bbox(c(xmin=unname(b["xmin"]-pad),ymin=unname(b["ymin"]-pad),
    xmax=unname(b["xmax"]+pad),ymax=unname(b["ymax"]+pad)),crs=crs))
  region <- sf::st_sf(geometry=region)
  extent <- .context_extent(region)
  p <- .context_base(region) + ggplot2::geom_sf(data=aoi,fill="#b2182b",colour="#b2182b",alpha=.65)
  for(layer in c("marine","places")) {
    labels <- sf::st_read("data/context/north_america.gpkg",layer=layer,quiet=TRUE)
    labels <- suppressWarnings(sf::st_crop(sf::st_transform(labels,crs),extent))
    if(nrow(labels)) {
      if(layer=="marine") labels <- suppressWarnings(sf::st_point_on_surface(labels))
      p <- p + ggplot2::geom_sf_text(data=labels,ggplot2::aes(label=name),size=3,
        colour=if(layer=="marine") "#386679" else "grey25",check_overlap=TRUE)
    }
  }
  p <- p + ggplot2::geom_sf_text(data=sf::st_sf(geometry=sf::st_centroid(sf::st_geometry(aoi))),
    label="AOI",nudge_y=pad*.15,size=3)
  .context_finish(p,region)
}
