library(ggplot2)
library(patchwork)

diform3d_tangent_color <- "#7B2CBF"
diform3d_obs_color <- "#111111"
diform3d_palette <- viridisLite::viridis(256)

#' Add transparency to colours.
diform3d_alpha <- function(col, alpha = 1) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

#' Map values to palette colours between vmin and vmax.
diform3d_color_values <- function(z, vmin, vmax, alpha = 1) {
  z <- pmax(vmin, pmin(vmax, z))
  idx <- floor((z - vmin) / (vmax - vmin) * (length(diform3d_palette) - 1)) + 1
  idx[!is.finite(idx)] <- 1
  diform3d_alpha(diform3d_palette[idx], alpha)
}

#' Colours of surface cells from the mean of their four corner values.
diform3d_surface_colors <- function(z, vmin, vmax, alpha = 0.74) {
  z_cell <- 0.25 * (z[-nrow(z), -ncol(z), drop = FALSE] + z[-1, -ncol(z), drop = FALSE] +
                      z[-nrow(z), -1, drop = FALSE] + z[-1, -1, drop = FALSE])
  diform3d_color_values(z_cell, vmin, vmax, alpha)
}

#' Open a 2400 x 1150 PNG device.
diform3d_open_png <- function(path) {
  ragg::agg_png(path, width = 2400, height = 1150, res = 160, background = "white")
}

#' Half-space normals and offsets of the contour in unnormalized working-space coordinates.
diform3d_constraints <- function(diform_result, rp) {
  rv <- diform_result$rv_matrix[, paste0("T", rp)]
  dirs <- diform_result$directions
  med <- diform_result$norm_info$med
  std <- diform_result$norm_info$std
  normal <- sweep(dirs, 2, std, "/")
  offset <- -as.numeric(normal %*% med) - rv
  list(normal = normal, offset = offset, rv = rv, dirs = dirs)
}

#' Convert working-space coordinates to current speed, angle and response.
diform3d_physical <- function(cspd, hl, ht) {
  list(Cspd = cspd, angularDifference = (atan2(ht, hl) * 180 / pi) %% 360, response = hl^2 + ht^2)
}

#' Observations in working-space coordinates.
diform3d_obs_working <- function(obs) {
  a <- obs$angularDifference * pi / 180
  h <- sqrt(obs$response)
  cbind(Cspd = obs$Cspd, H_L = h * cos(a), H_T = h * sin(a))
}

#' Axis limits of the physical and working-space plots.
diform3d_limits <- function(all_vertices, obs) {
  ow <- diform3d_obs_working(obs)
  a <- max(abs(c(all_vertices[, 2], all_vertices[, 3], ow[, "H_L"], ow[, "H_T"]))) * 1.08
  list(
    cspd = c(0, max(c(all_vertices[, 1], ow[, "Cspd"])) * 1.08),
    angle = c(0, 360),
    response = c(0, max(c(all_vertices[, 2]^2 + all_vertices[, 3]^2, obs$response)) * 1.08),
    working_y = c(-a, a),
    working_z = c(-a, a)
  )
}

#' Contour response on a grid of current speed and angle, from the half-space constraints.
diform3d_surface_physical <- function(constraints, c_grid, angle_grid) {
  C <- matrix(c_grid, nrow = length(c_grid), ncol = length(angle_grid))
  A <- matrix(angle_grid, nrow = length(c_grid), ncol = length(angle_grid), byrow = TRUE)
  n <- constraints$normal
  b <- constraints$offset
  ca <- cos(A * pi / 180)
  sa <- sin(A * pi / 180)
  lower <- matrix(0, nrow = nrow(C), ncol = ncol(C))
  upper <- matrix(Inf, nrow = nrow(C), ncol = ncol(C))
  valid <- matrix(TRUE, nrow = nrow(C), ncol = ncol(C))
  for (i in seq_len(nrow(n))) {
    alpha <- n[i, 2] * ca + n[i, 3] * sa
    beta <- n[i, 1] * C + b[i]
    rhs <- -beta
    pos <- alpha > 1e-12
    neg <- alpha < -1e-12
    flat <- abs(alpha) <= 1e-12
    upper[pos] <- pmin(upper[pos], rhs[pos] / alpha[pos])
    lower[neg] <- pmax(lower[neg], rhs[neg] / alpha[neg])
    valid[flat] <- valid[flat] & (beta[flat] <= 1e-10)
  }
  lower0 <- matrix(pmax(0, lower), nrow = nrow(lower), ncol = ncol(lower))
  valid <- valid & is.finite(upper) & upper >= lower0 & upper > 0
  hs <- upper^2
  hs[!valid] <- NA_real_
  hs
}

#' Normalize a vector to unit length.
diform3d_unit <- function(x) {
  x / sqrt(sum(x^2))
}

#' Cross product of two three-dimensional vectors.
diform3d_cross <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
}

#' Two orthonormal vectors spanning the plane orthogonal to a normal.
diform3d_plane_basis <- function(normal) {
  n <- diform3d_unit(normal)
  a <- if (abs(n[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  e1 <- diform3d_unit(a - sum(a * n) * n)
  list(e1 = e1, e2 = diform3d_unit(diform3d_cross(n, e1)))
}

#' Text keys identifying directions, used for matching.
diform3d_direction_key <- function(x) {
  apply(round(as.matrix(x), 10), 1, paste, collapse = ":")
}

#' Directions and triangles of the reflected simplex mesh for subdivision m.
diform3d_direction_mesh_for_m <- function(m) {
  base <- do.call(rbind, lapply(0:m, function(i) {
    do.call(rbind, lapply(0:(m - i), function(j) c(i, j, m - i - j)))
  }))
  base <- base / m
  base <- base[rowSums(base^2) > 1e-12, , drop = FALSE]
  signs <- as.matrix(expand.grid(c(1, -1), c(1, -1), c(1, -1)))
  dirs_all <- do.call(rbind, lapply(seq_len(nrow(signs)), function(i) sweep(base, 2, signs[i, ], "*")))
  dirs <- unique(round(dirs_all / sqrt(rowSums(dirs_all^2)), 10))
  colnames(dirs) <- c("d1", "d2", "d3")
  key_to_index <- setNames(seq_len(nrow(dirs)), diform3d_direction_key(dirs))

  #' Key of the normalized direction for simplex indices i, j and a sign pattern.
  point_key <- function(i, j, sign) {
    w <- sign * c(i, j, m - i - j) / m
    paste(round(w / sqrt(sum(w^2)), 10), collapse = ":")
  }

  tris <- list()
  for (s in seq_len(nrow(signs))) {
    sign <- signs[s, ]
    for (i in 0:(m - 1)) {
      for (j in 0:(m - i - 1)) {
        tri1 <- c(key_to_index[[point_key(i, j, sign)]], key_to_index[[point_key(i + 1, j, sign)]],
                  key_to_index[[point_key(i, j + 1, sign)]])
        if (length(unique(tri1)) == 3) tris[[length(tris) + 1L]] <- tri1
        if (i + j <= m - 2) {
          tri2 <- c(key_to_index[[point_key(i + 1, j, sign)]], key_to_index[[point_key(i + 1, j + 1, sign)]],
                    key_to_index[[point_key(i, j + 1, sign)]])
          if (length(unique(tri2)) == 3) tris[[length(tris) + 1L]] <- tri2
        }
      }
    }
  }
  list(dirs = dirs, triangles = unique(t(apply(do.call(rbind, tris), 1, sort))))
}

#' Triangulation of the direction set, recovered by matching it to the simplex mesh of the same size.
diform3d_direction_triangles <- function(dirs) {
  input_keys <- diform3d_direction_key(dirs)
  m <- 1L
  repeat {
    mesh <- diform3d_direction_mesh_for_m(m)
    if (nrow(mesh$dirs) == nrow(dirs)) break
    m <- m + 1L
  }
  matrix(match(diform3d_direction_key(mesh$dirs), input_keys)[mesh$triangles], ncol = 3)
}

#' Supporting plane of the contour at the vertex with the largest response, from its adjacent triangular facets.
diform3d_support_info <- function(vertices, triangles) {
  i_star <- which.max(vertices[, 2]^2 + vertices[, 3]^2)
  x_star <- vertices[i_star, ]
  interior <- colMeans(vertices)
  grad <- diform3d_unit(c(0, 2 * x_star[2], 2 * x_star[3]))
  facets <- triangles[apply(triangles == i_star, 1, any), , drop = FALSE]
  best_score <- -Inf
  for (j in seq_len(nrow(facets))) {
    pts <- vertices[facets[j, ], , drop = FALSE]
    normal <- diform3d_cross(pts[2, ] - pts[1, ], pts[3, ] - pts[1, ])
    offset <- -sum(normal * pts[1, ])
    if (sum(normal * interior) + offset > 0) {
      normal <- -normal
      offset <- -offset
    }
    score <- sum(diform3d_unit(normal) * grad)
    if (score > best_score) {
      best_score <- score
      best_normal <- normal
      best_offset <- offset
    }
  }
  list(x_star = as.numeric(x_star), normal = diform3d_unit(best_normal),
       offset = best_offset / sqrt(sum(best_normal^2)))
}

#' Grid of the supporting-plane boundary mapped to physical coordinates around the maximum.
diform3d_halfspace_physical_boundary_local <- function(support_info, limits, n_cspd = 70, n_ang = 90,
                                                       cspd_fraction = 0.32, angle_half_width = 34) {
  star <- diform3d_physical(support_info$x_star[1], support_info$x_star[2], support_info$x_star[3])
  c_half <- cspd_fraction * diff(limits$cspd)
  cg <- seq(max(limits$cspd[1], star$Cspd - c_half), min(limits$cspd[2], star$Cspd + c_half), length.out = n_cspd)
  ag <- seq(star$angularDifference - angle_half_width, star$angularDifference + angle_half_width, length.out = n_ang)
  ag <- ag[ag >= limits$angle[1] & ag <= limits$angle[2]]
  C <- matrix(cg, nrow = length(cg), ncol = length(ag))
  A <- matrix(ag, nrow = length(cg), ncol = length(ag), byrow = TRUE)
  n <- support_info$normal
  alpha <- n[2] * cos(A * pi / 180) + n[3] * sin(A * pi / 180)
  r <- -(n[1] * C + support_info$offset) / alpha
  h <- r^2
  h[abs(alpha) <= 1e-12 | r < 0 | h < limits$response[1] | h > limits$response[2]] <- NA_real_
  list(X = C, Y = A, Z = h)
}

#' Draw a grid of quadrilaterals on a perspective plot, skipping cells with missing values or a large vertical span.
diform3d_project_grid <- function(pm, X, Y, Z, col, border, max_z_span) {
  for (i in seq_len(nrow(X) - 1)) {
    for (j in seq_len(ncol(X) - 1)) {
      zz <- c(Z[i, j], Z[i + 1, j], Z[i + 1, j + 1], Z[i, j + 1])
      if (any(!is.finite(zz)) || diff(range(zz)) > max_z_span) next
      xy <- grDevices::trans3d(c(X[i, j], X[i + 1, j], X[i + 1, j + 1], X[i, j + 1]),
                               c(Y[i, j], Y[i + 1, j], Y[i + 1, j + 1], Y[i, j + 1]), zz, pm)
      graphics::polygon(xy$x, xy$y, col = col, border = border)
    }
  }
}

#' Polygon where the supporting plane cuts the plotting box.
diform3d_halfspace_box_plane <- function(support_info, limits) {
  verts <- as.matrix(expand.grid(x = limits[, 1], y = limits[, 2], z = limits[, 3]))
  edges <- combn(seq_len(8), 2)
  edge_diff <- abs(verts[edges[1, ], , drop = FALSE] - verts[edges[2, ], , drop = FALSE])
  edges <- edges[, rowSums(edge_diff > 1e-12) == 1, drop = FALSE]
  s <- as.numeric(verts %*% support_info$normal + support_info$offset)
  pts <- verts[s >= -1e-10, , drop = FALSE]
  for (e in seq_len(ncol(edges))) {
    i <- edges[1, e]
    j <- edges[2, e]
    if ((s[i] < -1e-10 && s[j] > 1e-10) || (s[i] > 1e-10 && s[j] < -1e-10)) {
      pts <- rbind(pts, verts[i, ] + s[i] / (s[i] - s[j]) * (verts[j, ] - verts[i, ]))
    }
  }
  pts <- unique(round(pts, 10))
  idx <- which(abs(as.numeric(pts %*% support_info$normal + support_info$offset)) < 1e-8)
  basis <- diform3d_plane_basis(support_info$normal)
  centered <- pts[idx, , drop = FALSE] - matrix(colMeans(pts[idx, , drop = FALSE]), length(idx), 3, byrow = TRUE)
  h <- grDevices::chull(as.numeric(centered %*% basis$e1), as.numeric(centered %*% basis$e2))
  pts[idx[h], , drop = FALSE]
}

#' Draw a vertical colour bar.
diform3d_draw_colorbar <- function(vmax, label) {
  old <- par(mar = c(5.4, 0.4, 4.6, 4.2), xaxs = "i", yaxs = "i")
  on.exit(par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, vmax))
  y <- seq(0, vmax, length.out = length(diform3d_palette) + 1L)
  graphics::rect(0.18, y[-length(y)], 0.72, y[-1], col = diform3d_palette, border = NA)
  graphics::rect(0.18, 0, 0.72, vmax, border = "black", lwd = 0.9)
  ticks <- pretty(c(0, vmax), n = 6)
  ticks <- ticks[ticks >= 0 & ticks <= vmax]
  graphics::axis(4, at = ticks, pos = 0.78, las = 1, cex.axis = 1.20)
  graphics::mtext(label, side = 4, line = 2.7, cex = 1.35)
}

#' Draw a framed block of text lines at the top left of the plot.
diform3d_add_annotation <- function(labels, cex = 1.25) {
  usr <- par("usr")
  x <- usr[1] - 0.090 * diff(usr[1:2])
  y <- usr[4] + 0.025 * diff(usr[3:4])
  line_h <- strheight("M", cex = cex) * 1.35
  pad_x <- diff(usr[1:2]) * 0.012
  pad_y <- diff(usr[3:4]) * 0.014
  h <- line_h * length(labels) + 2 * pad_y
  w <- max(strwidth(labels, cex = cex)) + 2 * pad_x
  graphics::rect(x - pad_x, y - h + pad_y, x + w - pad_x, y + pad_y,
                 col = diform3d_alpha("white", 0.90), border = "lightgray")
  graphics::text(x, y - seq_along(labels) * line_h + line_h * 0.45, labels = labels, adj = c(0, 0.5), cex = cex)
}

#' Draw the contour surface in physical coordinates.
diform3d_plot_physical <- function(vertices, constraints, support_info, all_vertices, obs, obs_bg,
                                   response_info, elev, azim, show_support_plane, vmax) {
  limits <- diform3d_limits(all_vertices, obs)
  c_grid <- seq(min(vertices[, 1]), max(vertices[, 1]), length.out = 80)
  a_grid <- seq(0, 360, length.out = 180)
  hs <- diform3d_surface_physical(constraints, c_grid, a_grid)

  old <- par(mar = c(5.2, 5.4, 5.0, 1.1), xpd = NA)
  on.exit(par(old), add = TRUE)

  pm <- graphics::persp(
    x = c_grid, y = a_grid, z = hs, xlim = limits$cspd, ylim = limits$angle, zlim = limits$response,
    col = diform3d_surface_colors(hs, 0, vmax, alpha = 0.76), border = NA, shade = NA, ticktype = "detailed",
    theta = azim, phi = elev, r = 7, expand = 0.88,
    xlab = "Current speed [m s⁻¹]", ylab = "Angular difference [°]", zlab = response_info$label,
    main = "(b) Physical space",
    cex.lab = 1.34, cex.axis = 1.06, cex.main = 1.42, font.main = 2
  )

  if (show_support_plane) {
    boundary <- diform3d_halfspace_physical_boundary_local(support_info, limits, n_cspd = 64, n_ang = 99)
    diform3d_project_grid(pm, boundary$X, boundary$Y, boundary$Z,
                          col = diform3d_alpha(diform3d_tangent_color, 0.34),
                          border = diform3d_alpha(diform3d_tangent_color, 0.24),
                          max_z_span = diff(limits$response) * 0.25)
    star <- diform3d_physical(support_info$x_star[1], support_info$x_star[2], support_info$x_star[3])
    graphics::points(grDevices::trans3d(star$Cspd, star$angularDifference, star$response, pm),
                     pch = 18, col = diform3d_tangent_color, cex = 1.08)
  }

  graphics::points(grDevices::trans3d(obs_bg$Cspd, obs_bg$angularDifference, obs_bg$response, pm),
                   pch = 16, cex = 0.17, col = diform3d_alpha(diform3d_obs_color, 0.040))

  phys <- diform3d_physical(vertices[, 1], vertices[, 2], vertices[, 3])
  i <- which.max(phys$response)
  diform3d_add_annotation(c(
    paste0("Max ", response_info$short, " on contour"),
    sprintf("%s = %.3f %s", response_info$short, phys$response[i], response_info$unit),
    sprintf("Current speed = %.2f m s⁻¹", phys$Cspd[i]),
    sprintf("Angular difference = %.0f°", phys$angularDifference[i])
  ))
}

#' Draw the contour in working-space coordinates.
diform3d_plot_working <- function(vertices, triangles, support_info, all_vertices, obs, obs_bg,
                                  response_info, elev, azim, show_support_plane, vmax) {
  limits <- diform3d_limits(all_vertices, obs)
  comp <- response_info$components

  old <- par(mar = c(4.8, 4.9, 4.8, 1.0), xpd = NA)
  on.exit(par(old), add = TRUE)

  pm <- graphics::persp(
    x = limits$cspd, y = limits$working_y, z = matrix(limits$working_z[1], 2, 2), zlim = limits$working_z,
    col = NA, border = NA, ticktype = "detailed", theta = azim, phi = elev, r = 7, expand = 0.88,
    xlab = "Current speed [m s⁻¹]", ylab = comp[1], zlab = comp[2],
    main = "(a) Transformed space",
    cex.lab = 1.28, cex.axis = 1.04, cex.main = 1.34, font.main = 2
  )

  polygons <- lapply(seq_len(nrow(triangles)), function(i) vertices[triangles[i, ], , drop = FALSE])
  depth <- vapply(polygons, function(p) mean(p[, 1] + p[, 2] + p[, 3]), numeric(1))
  for (p in polygons[order(depth, decreasing = TRUE)]) {
    xy <- grDevices::trans3d(p[, 1], p[, 2], p[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = diform3d_color_values(mean(p[, 2]^2 + p[, 3]^2), 0, vmax, 0.58), border = NA)
  }

  if (show_support_plane) {
    lim <- rbind(c(limits$cspd[1], limits$working_y[1], limits$working_z[1]),
                 c(limits$cspd[2], limits$working_y[2], limits$working_z[2]))
    plane <- diform3d_halfspace_box_plane(support_info, lim)
    xy <- grDevices::trans3d(plane[, 1], plane[, 2], plane[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = diform3d_alpha(diform3d_tangent_color, 0.32),
                      border = diform3d_alpha(diform3d_tangent_color, 0.88), lwd = 0.8)
    graphics::points(grDevices::trans3d(support_info$x_star[1], support_info$x_star[2], support_info$x_star[3], pm),
                     pch = 18, col = diform3d_tangent_color, cex = 1.08)
  }

  bg_pts <- diform3d_obs_working(obs_bg)
  graphics::points(grDevices::trans3d(bg_pts[, 1], bg_pts[, 2], bg_pts[, 3], pm),
                   pch = 16, cex = 0.17, col = diform3d_alpha(diform3d_obs_color, 0.040))
}

#' Draw the working-space and physical contours side by side for every return period.
plot_diform_3d_r <- function(diform_result, data, response_info, out_dir, elev = 45, azim = 40,
                             show_support_plane = FALSE) {
  obs <- data.frame(Cspd = data$Cspd, angularDifference = data$angleDiff, response = data[[diform_result$response]])
  obs_bg <- obs[unique(round(seq(1, nrow(obs), length.out = 20000))), , drop = FALSE]
  rps <- diform_result$return_periods
  all_vertices <- do.call(rbind, lapply(diform_result$contours_3d, `[[`, "vertices_3d"))
  triangles <- diform3d_direction_triangles(diform_result$directions)
  vmax <- max(all_vertices[, 2]^2 + all_vertices[, 3]^2) * 1.03
  suffix <- if (show_support_plane) "_with_support_plane" else ""

  for (rp in rps) {
    vertices <- diform_result$contours_3d[[as.character(rp)]]$vertices_3d
    constraints <- diform3d_constraints(diform_result, rp)
    support_info <- diform3d_support_info(vertices, triangles)
    diform3d_open_png(file.path(out_dir, sprintf("diform_3d_%syr%s.png", format(rp, trim = TRUE), suffix)))
    layout(matrix(1:3, 1, 3), widths = c(8.3, 8.3, 1.3))
    par(cex = 1)
    diform3d_plot_working(vertices, triangles, support_info, all_vertices, obs, obs_bg,
                          response_info, elev, azim, show_support_plane, vmax)
    diform3d_plot_physical(vertices, constraints, support_info, all_vertices, obs, obs_bg,
                           response_info, elev, azim, show_support_plane, vmax)
    diform3d_draw_colorbar(vmax, response_info$label)
    grDevices::dev.off()
  }
}

diform2d_theme <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 11), axis.title = element_text(size = 13), plot.margin = margin(6, 14, 6, 6))

#' Contour response on a grid of current speed and angle, as a long table.
diform2d_surface <- function(diform_result, rp, c_grid, a_grid) {
  z <- diform3d_surface_physical(diform3d_constraints(diform_result, rp), c_grid, a_grid)
  data.table(Cspd = rep(c_grid, times = length(a_grid)), angle = rep(a_grid, each = length(c_grid)),
             value = as.vector(z), rp = rp)
}

#' Top view of the physical space coloured by the response, with the cut bands and the maximum.
diform2d_map <- function(s, cuts, pk, vmax, dc, da) {
  s <- s[!is.na(value)]
  cmax <- max(s$Cspd)
  bands <- rbindlist(lapply(cuts, function(x) {
    if (x$along == "angle") data.table(xmin = 0, xmax = 360, ymin = x$at - dc, ymax = x$at + dc,
                                       x0 = 0, x1 = 360, y0 = x$at, y1 = x$at, lx = 350, ly = x$at, tag = x$tag)
    else data.table(xmin = x$at - da, xmax = x$at + da, ymin = 0, ymax = cmax,
                    x0 = x$at, x1 = x$at, y0 = 0, y1 = cmax, lx = x$at, ly = cmax - 0.25, tag = x$tag)
  }))
  ggplot(s, aes(angle, Cspd)) +
    geom_raster(aes(fill = value), interpolate = TRUE) +
    geom_contour(aes(z = value), colour = "white", linewidth = 0.3, alpha = 0.6, bins = 10) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
              fill = "white", alpha = 0.3, colour = NA) +
    geom_segment(data = bands, aes(x = x0, xend = x1, y = y0, yend = y1), inherit.aes = FALSE,
                 colour = "grey15", linewidth = 0.35, linetype = "dashed", alpha = 0.8) +
    geom_label(data = bands, aes(lx, ly, label = tag), inherit.aes = FALSE, size = 4.2, fontface = "bold",
               colour = "grey15", fill = alpha("white", 0.75), linewidth = 0, label.padding = unit(0.1, "lines")) +
    geom_point(data = pk, shape = 4, size = 4, stroke = 1.4, colour = "black") +
    scale_fill_viridis_c(limits = c(0, vmax), guide = "none") +
    scale_x_continuous(breaks = seq(0, 360, 45), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "(c) Top view of the physical space", x = "Angular difference [°]", y = "Current speed [m s⁻¹]") +
    diform2d_theme
}

#' Cut of every return-period contour at fixed current speed or fixed angle, over the observations in the band.
diform2d_cut <- function(diform_result, data, response_info, cut, c_grid, a_grid, dc, da) {
  rps <- diform_result$return_periods
  y <- data[[diform_result$response]]
  if (cut$along == "angle") {
    sl <- rbindlist(lapply(rps, function(rp) diform2d_surface(diform_result, rp, cut$at, a_grid)))[, x := angle]
    keep <- abs(data$Cspd - cut$at) < dc
    ob <- data.table(x = data$angleDiff[keep], y = y[keep])
    xs <- scale_x_continuous(breaks = seq(0, 360, 90), limits = c(0, 360), expand = c(0, 0))
    lb <- labs(title = sprintf("%s Current speed = %.2f m s⁻¹", cut$tag, cut$at), x = "Angular difference [°]")
  } else {
    sl <- rbindlist(lapply(rps, function(rp) diform2d_surface(diform_result, rp, c_grid, cut$at)))[, x := Cspd]
    d <- abs(data$angleDiff - cut$at)
    keep <- pmin(d, 360 - d) < da
    ob <- data.table(x = data$Cspd[keep], y = y[keep])
    xs <- scale_x_continuous(limits = c(0, ceiling(max(c_grid))), expand = c(0, 0))
    lb <- labs(title = sprintf("%s Angular difference = %.0f°", cut$tag, cut$at), x = "Current speed [m s⁻¹]")
  }
  ggplot() +
    geom_point(data = ob, aes(x, y), colour = "grey60", size = 0.4, alpha = 0.25, shape = 16) +
    geom_vline(xintercept = cut$mark, linetype = "dashed", linewidth = 0.4, colour = "grey30") +
    geom_line(data = sl[!is.na(value)], aes(x, value, colour = factor(rp)), linewidth = 0.7) +
    scale_colour_manual(values = setNames(viridisLite::viridis(length(rps), end = 0.9), rps), name = "Return period [yr]") +
    xs + lb + labs(y = response_info$label) + diform2d_theme + theme(plot.title = element_text(size = 13))
}

#' Composite figure for every return period: 3D panels, top view of the physical space and cuts through the maximum.
plot_diform_composite <- function(diform_result, data, response_info, out_dir, second_ridge = FALSE, dc = 0.1, da = 2.5) {
  panel_dir <- tempfile("diform3d_")
  dir.create(panel_dir)
  on.exit(unlink(panel_dir, recursive = TRUE), add = TRUE)
  plot_diform_3d_r(diform_result, data, response_info, panel_dir)
  all_vertices <- do.call(rbind, lapply(diform_result$contours_3d, `[[`, "vertices_3d"))
  vmax <- max(all_vertices[, 2]^2 + all_vertices[, 3]^2) * 1.03
  c_grid <- seq(0, max(all_vertices[, 1]), by = 0.01)
  a_grid <- seq(0, 360, by = 0.5)

  for (rp in diform_result$return_periods) {
    tag <- format(rp, trim = TRUE)
    s <- diform2d_surface(diform_result, rp, c_grid, a_grid)
    bt <- diform_result$contours_3d[[as.character(rp)]]$vertices_bt
    i <- which.max(bt[[diform_result$response]])
    pk <- data.table(Cspd = bt$Cspd[i], angle = bt$angleDiff[i])
    cuts <- list(list(along = "angle", at = pk$Cspd, mark = pk$angle), list(along = "cspd", at = pk$angle, mark = pk$Cspd))
    if (second_ridge) {
      other <- s[!is.na(value) & (angle < 180) != (pk$angle < 180)]
      r2 <- other[which.max(value)]
      cuts[[1]]$mark <- c(pk$angle, r2$angle)
      cuts[[3]] <- list(along = "cspd", at = r2$angle, mark = r2$Cspd)
    }
    for (k in seq_along(cuts)) cuts[[k]]$tag <- sprintf("(%s)", letters[k + 3])

    img <- png::readPNG(file.path(panel_dir, sprintf("diform_3d_%syr.png", tag)))
    h_top <- 12 * nrow(img) / ncol(img)
    h_bottom <- 2.6 + 1.8 * length(cuts)
    right <- wrap_plots(lapply(cuts, function(x) diform2d_cut(diform_result, data, response_info, x, c_grid, a_grid, dc, da)), ncol = 1)
    bottom <- (diform2d_map(s, cuts, pk, vmax, dc, da) | right) + plot_layout(widths = c(1.9, 1), guides = "collect")
    p <- wrap_plots(wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE)), bottom, ncol = 1, heights = c(h_top, h_bottom))
    ggsave(file.path(out_dir, sprintf("diform_composite_%syr.png", tag)), p, width = 12, height = h_top + h_bottom,
           dpi = 250, device = ragg::agg_png)
  }
}
