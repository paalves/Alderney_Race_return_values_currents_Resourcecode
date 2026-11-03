required_packages <- c(
  "ggplot2", "sf", "patchwork", "rnaturalearth", "rnaturalearthdata",
  "resourcecodedata"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install the missing packages first: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

site_lon <- -2.010
site_lat <- 49.712
site_node <- 158348L
show_mesh <- TRUE

regional_extent <- c(xmin = -5.0, xmax = 10.0, ymin = 41.0, ymax = 52.0)
inset_extent <- c(
  xmin = site_lon - 0.40, xmax = site_lon + 0.40,
  ymin = site_lat - 0.40, ymax = site_lat + 0.40
)

inset_position <- c(left = 0.657, bottom = 0.510, right = 0.985, top = 0.985)
connector_end <- c(
  longitude = unname(
    regional_extent[["xmin"]] +
      inset_position[["left"]] * diff(regional_extent[c("xmin", "xmax")])
  ),
  latitude = unname(
    regional_extent[["ymin"]] +
      inset_position[["bottom"]] * diff(regional_extent[c("ymin", "ymax")])
  )
)

colours <- c(
  sea = "#EAF4F8",
  land = "#F5F1E8",
  coastline = "#25343B",
  border = "#87939A",
  grid = "#B9CDD5",
  mesh = "#4CA6A8",
  accent = "#D1495B"
)

utils::data(
  list = c("rscd_field", "rscd_triangles", "rscd_coastline", "rscd_islands"),
  package = "resourcecodedata",
  envir = environment()
)

field <- rscd_field
triangles <- as.matrix(rscd_triangles)

if (nrow(triangles) != 3L && ncol(triangles) == 3L) {
  triangles <- t(triangles)
}
if (nrow(triangles) != 3L) {
  stop("rscd_triangles must have three rows (one node index per vertex).", call. = FALSE)
}

distance2 <-
  ((field$longitude - site_lon) * cos(site_lat * pi / 180))^2 +
  (field$latitude - site_lat)^2
nearest_row <- which.min(distance2)
nearest_node <- field$node[nearest_row]

if (is.finite(site_node) && site_node %in% field$node) {
  documented_row <- match(site_node, field$node)
  documented_distance <- sqrt(
    ((field$longitude[documented_row] - site_lon) * cos(site_lat * pi / 180))^2 +
      (field$latitude[documented_row] - site_lat)^2
  )
  if (documented_distance > 0.03) {
    warning(
      "FIELD node ", site_node, " is not the closest node to the requested coordinates; ",
      "nearest node is ", nearest_node, "."
    )
  }
}

site <- data.frame(longitude = site_lon, latitude = site_lat)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
regional_land <- suppressWarnings(sf::st_crop(world, regional_extent))


inset_box <- sf::st_as_sfc(sf::st_bbox(
  c(
    xmin = inset_extent[["xmin"]], ymin = inset_extent[["ymin"]],
    xmax = inset_extent[["xmax"]], ymax = inset_extent[["ymax"]]
  ),
  crs = sf::st_crs(4326)
))

resource_ocean <- sf::st_sfc(
  sf::st_polygon(list(as.matrix(
    rscd_coastline[, c("longitude", "latitude")]
  ))),
  crs = sf::st_crs(4326)
)
resource_ocean <- sf::st_make_valid(resource_ocean)

local_crs <- sf::st_crs(32630)
inset_box_utm <- sf::st_transform(inset_box, local_crs)
resource_ocean_utm <- sf::st_transform(resource_ocean, local_crs)
inset_ocean_utm <- suppressWarnings(
  sf::st_intersection(resource_ocean_utm, inset_box_utm)
)
inset_mainland <- suppressWarnings(
  sf::st_difference(inset_box_utm, sf::st_union(inset_ocean_utm))
)
inset_mainland <- sf::st_transform(inset_mainland, 4326)

island_groups <- split(rscd_islands, rscd_islands$ID)
island_overlaps_inset <- vapply(island_groups, function(island) {
  max(island$longitude) >= inset_extent[["xmin"]] &&
    min(island$longitude) <= inset_extent[["xmax"]] &&
    max(island$latitude) >= inset_extent[["ymin"]] &&
    min(island$latitude) <= inset_extent[["ymax"]]
}, logical(1))
island_groups <- island_groups[island_overlaps_inset]

inset_islands <- sf::st_sf(
  ID = as.integer(names(island_groups)),
  geometry = sf::st_sfc(
    lapply(island_groups, function(island) {
      sf::st_polygon(list(as.matrix(
        island[, c("longitude", "latitude")]
      )))
    }),
    crs = sf::st_crs(4326)
  )
)
inset_islands <- sf::st_make_valid(inset_islands)

inset_water_utm <- inset_ocean_utm
if (nrow(inset_islands)) {
  inset_islands_utm <- sf::st_transform(inset_islands, local_crs)
  inset_water_utm <- suppressWarnings(
    sf::st_difference(inset_ocean_utm, sf::st_union(inset_islands_utm))
  )
}

#' Edges of the hindcast mesh triangles lying within the map extent.
mesh_segments <- function(field, triangles, extent, padding = 0.03) {
  inside <- with(
    field,
    longitude >= extent[["xmin"]] - padding &
      longitude <= extent[["xmax"]] + padding &
      latitude >= extent[["ymin"]] - padding &
      latitude <= extent[["ymax"]] + padding
  )

  if (identical(as.integer(field$node), seq_len(nrow(field)))) {
    triangle_rows <- triangles
  } else {
    triangle_rows <- matrix(
      match(as.vector(triangles), field$node),
      nrow = 3L
    )
  }

  valid <- colSums(is.na(triangle_rows)) == 0L
  triangle_rows <- triangle_rows[, valid, drop = FALSE]
  local <- colSums(
    matrix(inside[triangle_rows], nrow = 3L)
  ) == 3L
  triangle_rows <- t(triangle_rows[, local, drop = FALSE])

  if (!nrow(triangle_rows)) {
    return(data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
  }

  edge_from <- c(triangle_rows[, 1], triangle_rows[, 2], triangle_rows[, 3])
  edge_to <- c(triangle_rows[, 2], triangle_rows[, 3], triangle_rows[, 1])
  edges <- unique(data.frame(
    from = pmin(edge_from, edge_to),
    to = pmax(edge_from, edge_to)
  ))

  data.frame(
    x = field$longitude[edges$from],
    y = field$latitude[edges$from],
    xend = field$longitude[edges$to],
    yend = field$latitude[edges$to]
  )
}

mesh <- if (show_mesh) {
  mesh_segments(field, triangles, inset_extent)
} else {
  data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
}

mesh_clipped <- if (nrow(mesh)) {
  mesh_lines <- lapply(seq_len(nrow(mesh)), function(i) {
    matrix(
      c(mesh$x[i], mesh$y[i], mesh$xend[i], mesh$yend[i]),
      nrow = 2L, byrow = TRUE
    )
  })
  mesh_geometry <- sf::st_sfc(
    sf::st_multilinestring(mesh_lines),
    crs = sf::st_crs(4326)
  )
  mesh_geometry <- sf::st_transform(mesh_geometry, local_crs)
  mesh_geometry <- suppressWarnings(
    sf::st_intersection(mesh_geometry, inset_water_utm)
  )
  sf::st_sf(geometry = sf::st_transform(mesh_geometry, 4326))
} else {
  sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(4326)))
}


#' Longitude tick labels in degrees with an E/W suffix.
longitude_labels <- function(x) {
  suffix <- ifelse(x < 0, "W", ifelse(x > 0, "E", ""))
  paste0(format(abs(x), trim = TRUE, scientific = FALSE), "\u00b0", suffix)
}

#' Latitude tick labels in degrees with an N/S suffix.
latitude_labels <- function(x) {
  suffix <- ifelse(x < 0, "S", ifelse(x > 0, "N", ""))
  paste0(format(abs(x), trim = TRUE, scientific = FALSE), "\u00b0", suffix)
}

map_theme <- ggplot2::theme_minimal(base_family = "sans", base_size = 11) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = colours[["sea"]], colour = NA),
    panel.grid.major = ggplot2::element_line(
      colour = colours[["grid"]], linewidth = 0.3, linetype = "dashed"
    ),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(colour = "#52636B", size = 9),
    axis.ticks = ggplot2::element_line(colour = colours[["grid"]]),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(7, 7, 7, 7)
  )


regional_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = regional_land,
    fill = colours[["land"]],
    colour = colours[["border"]],
    linewidth = 0.3
  ) +
  ggplot2::annotate(
    "segment",
    x = inset_extent[["xmax"]], y = inset_extent[["ymin"]],
    xend = connector_end[["longitude"]], yend = connector_end[["latitude"]],
    colour = colours[["accent"]], linewidth = 1.15, lineend = "round"
  ) +
  ggplot2::annotate(
    "rect",
    xmin = inset_extent[["xmin"]], xmax = inset_extent[["xmax"]],
    ymin = inset_extent[["ymin"]], ymax = inset_extent[["ymax"]],
    fill = NA, colour = colours[["accent"]], linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = site,
    ggplot2::aes(x = longitude, y = latitude),
    shape = 21, size = 3.2, stroke = 1,
    fill = colours[["accent"]], colour = "white"
  ) +
  ggplot2::annotate(
    "segment",
    x = -4.35, xend = -4.35, y = 50.85, yend = 51.35,
    colour = colours[["coastline"]], linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed")
  ) +
  ggplot2::annotate(
    "text", x = -4.35, y = 51.53, label = "N",
    colour = colours[["coastline"]], fontface = "bold", size = 3.4
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(-5, 10, by = 5), labels = longitude_labels
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(42, 52, by = 2), labels = latitude_labels
  ) +
  ggplot2::coord_sf(
    crs = sf::st_crs(4326),
    xlim = regional_extent[c("xmin", "xmax")],
    ylim = regional_extent[c("ymin", "ymax")],
    expand = FALSE,
    clip = "on"
  ) +
  map_theme +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(
      fill = NA, colour = "#6F7E84", linewidth = 0.5
    ),
    axis.text = ggplot2::element_text(
      size = 12, colour = "#465860"
    ),
    axis.ticks = ggplot2::element_line(
      colour = "#73878F", linewidth = 0.45
    ),
    axis.ticks.length = grid::unit(0.12, "cm")
  )


inset_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = mesh_clipped,
    colour = colours[["mesh"]], alpha = 0.20, linewidth = 0.14
  ) +
  ggplot2::geom_sf(
    data = inset_mainland,
    fill = colours[["land"]], colour = NA
  ) +
  ggplot2::geom_sf(
    data = inset_islands,
    fill = colours[["land"]],
    colour = colours[["coastline"]],
    linewidth = 0.45
  ) +
  ggplot2::geom_path(
    data = rscd_coastline,
    ggplot2::aes(x = longitude, y = latitude),
    colour = colours[["coastline"]], linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = site,
    ggplot2::aes(x = longitude, y = latitude),
    shape = 21, size = 4.0, stroke = 1.1,
    fill = colours[["accent"]], colour = "white"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(-2.4, -2.0, -1.6), labels = longitude_labels
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(49.4, 49.8, 50.1), labels = latitude_labels
  ) +
  ggplot2::coord_sf(
    crs = sf::st_crs(4326),
    xlim = inset_extent[c("xmin", "xmax")],
    ylim = inset_extent[c("ymin", "ymax")],
    expand = FALSE,
    clip = "on"
  ) +
  map_theme +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(
      fill = NA, colour = colours[["accent"]], linewidth = 1.15
    ),
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

location_map <- regional_map +
  patchwork::inset_element(
    inset_map,
    left = inset_position[["left"]], bottom = inset_position[["bottom"]],
    right = inset_position[["right"]], top = inset_position[["top"]],
    align_to = "panel"
  )


output_dir <- figure_dir("site_description")
png_file <- file.path(output_dir, "alderney_race_location_map.png")
pdf_file <- file.path(output_dir, "alderney_race_location_map.pdf")

ggplot2::ggsave(
  png_file, location_map,
  width = 8.4, height = 7.1, units = "in", dpi = 450, bg = "white"
)
ggplot2::ggsave(
  pdf_file, location_map,
  width = 8.4, height = 7.1, units = "in", device = grDevices::cairo_pdf,
  bg = "white"
)

message("Saved: ", normalizePath(png_file, winslash = "/", mustWork = FALSE))
message("Saved: ", normalizePath(pdf_file, winslash = "/", mustWork = FALSE))

print(location_map)
