diform3d_alpha <- function(col, alpha = 1) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

diform3d_tangent_color <- "#7B2CBF"
diform3d_halfspace_color <- "#7A7A7A"
diform3d_obs_color <- "#111111"
diform3d_orange_point <- "#D97706"

diform3d_sample_label <- function(diform_result) {
  n_years <- diform_result$n_years
  n_years <- if (is.finite(n_years) && abs(n_years - round(n_years)) < 1e-6) {
    round(n_years)
  } else {
    round(n_years, 2)
  }
  source <- if (isFALSE(diform_result$bootstrap)) "observed data" else "simulation"
  sprintf("%s yr %s", format(n_years, big.mark = ",", scientific = FALSE), source)
}

diform3d_palette <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "#30123B", "#4145AB", "#4675ED", "#39A2FC",
    "#1BCFD4", "#24E36F", "#A4FC3C", "#F1E51D",
    "#FEA02F", "#F04D25", "#BA0C2F", "#7A0403"
  ))(n)
}

diform3d_hull_boundary_2d <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(data.frame(x = x, y = y))
  i <- tryCatch(c(grDevices::chull(x, y), grDevices::chull(x, y)[1]),
                error = function(e) seq_along(x))
  data.frame(x = x[i], y = y[i])
}

diform3d_polar_envelope <- function(H_L, H_T, n = 361) {
  ok <- is.finite(H_L) & is.finite(H_T)
  H_L <- H_L[ok]
  H_T <- H_T[ok]
  if (length(H_L) < 3) return(data.frame(angDiff = numeric(), resp = numeric()))

  h <- grDevices::chull(H_L, H_T)
  p <- rbind(cbind(H_L[h], H_T[h]), c(H_L[h[1]], H_T[h[1]]))
  a <- seq(0, 360, length.out = n)
  r <- numeric(n)

  for (i in seq_along(a)) {
    u <- c(cospi(a[i] / 180), sinpi(a[i] / 180))
    ri <- -Inf

    for (j in seq_len(nrow(p) - 1)) {
      v <- p[j + 1, ] - p[j, ]
      M <- cbind(u, -v)
      if (abs(det(M)) < 1e-10) next
      z <- solve(M, p[j, ])
      if (z[1] >= 0 && z[2] >= 0 && z[2] <= 1) ri <- max(ri, z[1])
    }

    r[i] <- if (is.finite(ri)) ri else NA_real_
  }

  data.frame(angDiff = a, resp = r^2)
}

diform3d_contour_slice_hs <- function(vertices, constraints, hs_value) {
  polygons <- diform3d_working_polygons(vertices, constraints)
  if (!length(polygons)) return(list())

  edge_intersections <- function(p1, p2) {
    d <- p2 - p1
    a <- d[2]^2 + d[3]^2
    b <- 2 * (p1[2] * d[2] + p1[3] * d[3])
    c <- p1[2]^2 + p1[3]^2 - hs_value

    if (abs(a) < 1e-12) {
      if (abs(b) < 1e-12) return(matrix(numeric(), ncol = 3))
      t <- -c / b
      if (is.finite(t) && t >= -1e-8 && t <= 1 + 1e-8) {
        return(matrix(p1 + pmin(1, pmax(0, t)) * d, nrow = 1))
      }
      return(matrix(numeric(), ncol = 3))
    }

    disc <- b^2 - 4 * a * c
    if (!is.finite(disc) || disc < -1e-10) return(matrix(numeric(), ncol = 3))
    disc <- max(0, disc)
    roots <- c((-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a))
    roots <- unique(round(roots[is.finite(roots) & roots >= -1e-8 & roots <= 1 + 1e-8], 12))
    if (!length(roots)) return(matrix(numeric(), ncol = 3))
    t(vapply(roots, function(t) p1 + pmin(1, pmax(0, t)) * d, numeric(3)))
  }

  out <- list()
  k <- 0L
  for (poly in polygons) {
    pts <- matrix(numeric(), ncol = 3)
    n <- nrow(poly)
    if (n < 3) next
    for (i in seq_len(n)) {
      j <- if (i == n) 1L else i + 1L
      pts <- rbind(pts, edge_intersections(poly[i, ], poly[j, ]))
    }
    if (nrow(pts) < 2) next
    pts <- unique(round(pts, 10))
    if (nrow(pts) < 2) next

    if (nrow(pts) > 2) {
      center <- colMeans(pts)
      pts <- pts[order(rowSums((pts - matrix(center, nrow(pts), 3, byrow = TRUE))^2),
                       decreasing = TRUE), , drop = FALSE]
    }

    p <- diform3d_physical(pts[1:2, 1], pts[1:2, 2], pts[1:2, 3])
    if (abs(diff(p$angularDifference)) > 180) next
    k <- k + 1L
    out[[k]] <- data.frame(Cspd = p$Cspd, angularDifference = p$angularDifference)
  }

  out
}

diform3d_color_values <- function(z, vmin, vmax, palette, alpha = 1) {
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) vmax <- vmin + 1
  z <- pmax(vmin, pmin(vmax, z))
  idx <- floor((z - vmin) / (vmax - vmin) * (length(palette) - 1)) + 1
  idx[!is.finite(idx)] <- 1
  diform3d_alpha(palette[idx], alpha)
}

diform3d_surface_colors <- function(z, vmin, vmax, palette, alpha = 0.74) {
  z_cell <- 0.25 * (
    z[-nrow(z), -ncol(z), drop = FALSE] +
      z[-1, -ncol(z), drop = FALSE] +
      z[-nrow(z), -1, drop = FALSE] +
      z[-1, -1, drop = FALSE]
  )
  diform3d_color_values(z_cell, vmin, vmax, palette, alpha)
}

diform3d_open_png <- function(path, width = 1520, height = 1200, res = 160) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    grDevices::png(path, width = width, height = height, res = res,
                   bg = "white", type = "cairo-png")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    grDevices::png(path, width = width, height = height, res = res, bg = "white")
  }
}

diform3d_axis_max <- function(value, default) {
  if (is.null(value)) return(default)
  value <- as.numeric(value)[1]
  if (!is.finite(value)) return(default)
  value
}

diform3d_get_rps <- function(diform_result, rp = "all") {
  available <- as.numeric(names(diform_result$contours_3d))
  available <- available[is.finite(available)]
  if (identical(rp, "all")) return(sort(available))
  sort(intersect(as.numeric(rp), available))
}

diform3d_warn_if_nonmonotone_returns <- function(diform_result, tol = 1e-8) {
  rv <- diform_result$rv_matrix
  rps <- as.numeric(diform_result$return_periods)
  if (is.null(rv) || is.null(rps) || ncol(rv) < 2) return(invisible(FALSE))

  o <- order(rps)
  rps <- rps[o]
  rv <- rv[, o, drop = FALSE]
  d <- rv[, -1, drop = FALSE] - rv[, -ncol(rv), drop = FALSE]
  bad <- which(is.finite(d) & d < -tol, arr.ind = TRUE)
  if (!nrow(bad)) return(invisible(FALSE))

  b <- bad[1, ]
  warning(
    sprintf(
      paste0(
        "D-IFORM directional return values are not monotone in %d direction/period ",
        "interval(s). Contours may not be nested. First issue: direction %d, ",
        "T=%s -> T=%s, %.4g -> %.4g."
      ),
      nrow(bad),
      b[1],
      format(rps[b[2]], trim = TRUE),
      format(rps[b[2] + 1], trim = TRUE),
      rv[b[1], b[2]],
      rv[b[1], b[2] + 1]
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

diform3d_rv_for_rp <- function(diform_result, rp) {
  cn <- paste0("T", rp)
  if (cn %in% colnames(diform_result$rv_matrix)) return(diform_result$rv_matrix[, cn])
  k <- match(rp, diform_result$return_periods)
  if (is.na(k)) stop("Return period not found in D-IFORM result: ", rp)
  diform_result$rv_matrix[, k]
}

diform3d_constraints <- function(diform_result, rp) {
  rv <- diform3d_rv_for_rp(diform_result, rp)
  dirs <- as.matrix(diform_result$directions)
  med <- as.numeric(diform_result$norm_info$med)
  std <- as.numeric(diform_result$norm_info$std)
  normal <- sweep(dirs, 2, std, "/")
  offset <- -as.numeric(normal %*% med) - rv
  list(normal = normal, offset = offset, rv = rv, med = med, std = std, dirs = dirs)
}

diform3d_vertices <- function(diform_result, rp) {
  x <- diform_result$contours_3d[[as.character(rp)]]$vertices_3d
  nm <- diform_result$component_names
  if (is.null(nm) || length(nm) != 3) nm <- c("Cspd", "H_L", "H_T")
  colnames(x) <- nm
  as.matrix(x)
}

diform3d_polar_meta <- function(diform_result = NULL, diform_spec = NULL) {
  if (is.null(diform_spec) && !is.null(diform_result)) diform_spec <- diform_result$diform_spec
  diform_spec <- normalise_diform_spec(diform_spec)
  linear <- diform_spec$components[[1]]
  circular <- diform_spec$components[[2]]
  list(
    linear = linear$variable,
    angle = circular$angle,
    response = diform_spec$response,
    h_names = circular$names,
    angle_unit = if (is.null(circular$angle_unit)) "degrees" else circular$angle_unit,
    period = if (is.null(circular$period)) 360 else circular$period,
    origin = if (is.null(circular$origin)) 0 else circular$origin,
    diform_spec = diform_spec
  )
}

diform3d_spec_field <- function(diform_spec, field, variable) {
  if (is.null(diform_spec) || is.null(diform_spec[[field]])) return(NULL)
  x <- diform_spec[[field]]
  if (!is.null(x[[variable]])) return(x[[variable]])
  NULL
}

diform3d_pretty_name <- function(variable) {
  if (identical(variable, "Cspd")) return("Current speed")
  if (identical(variable, "angularDifference") || identical(variable, "angleDiff")) {
    return("Angular difference")
  }
  variable
}

diform3d_default_unit <- function(variable, angle_unit = "degrees", for_value = FALSE) {
  if (identical(variable, "Cspd")) return("m.s^-1")
  if (identical(variable, "angularDifference") || identical(variable, "angleDiff")) {
    if (identical(angle_unit, "radians")) return("rad")
    return(if (isTRUE(for_value)) "\u00B0" else "deg")
  }
  NULL
}

diform3d_normalise_value_unit <- function(unit) {
  if (is.null(unit) || !nzchar(unit)) return(NULL)
  if (unit %in% c("deg", "degree", "degrees")) return("\u00B0")
  unit
}

diform3d_label_for <- function(variable, diform_spec = NULL, angle_unit = "degrees") {
  label <- diform3d_spec_field(diform_spec, "variable_labels", variable)
  unit <- diform3d_spec_field(diform_spec, "variable_units", variable)
  if (is.null(unit) || !nzchar(unit)) {
    unit <- diform3d_default_unit(variable, angle_unit = angle_unit, for_value = FALSE)
  }
  name <- if (!is.null(label)) label else diform3d_pretty_name(variable)
  if (is.null(unit) || !nzchar(unit)) return(name)
  paste0(name, " [", unit, "]")
}

diform3d_display_name_for <- function(variable, diform_spec = NULL) {
  label <- diform3d_spec_field(diform_spec, "variable_labels", variable)
  if (!is.null(label)) label else diform3d_pretty_name(variable)
}

diform3d_working_labels <- function(meta) {
  if (is.null(meta$h_names) || length(meta$h_names) != 2) {
    meta$h_names <- c("H_L", "H_T")
  }
  list(
    x = diform3d_label_for(meta$linear, meta$diform_spec),
    y = diform3d_label_for(meta$h_names[1], meta$diform_spec),
    z = diform3d_label_for(meta$h_names[2], meta$diform_spec),
    title = paste(
      diform3d_display_name_for(meta$linear, meta$diform_spec),
      diform3d_display_name_for(meta$h_names[1], meta$diform_spec),
      diform3d_display_name_for(meta$h_names[2], meta$diform_spec),
      sep = ", "
    )
  )
}

diform3d_component_labels <- function(comp_names, diform_spec = NULL,
                                      response_label = NULL,
                                      response_col = NULL) {
  out <- vapply(comp_names, diform3d_label_for, character(1), diform_spec = diform_spec)
  if (!is.null(response_col) && !is.null(response_label) &&
      length(out) >= response_col && nzchar(response_label)) {
    response_component <- NULL
    if (!is.null(diform_spec)) {
      response_component <- diform_spec$response_component
      if (is.null(response_component)) response_component <- diform_spec$response
    }
    if (is.null(response_component) || identical(comp_names[response_col], response_component)) {
      out[response_col] <- response_label
    }
  }
  out
}

diform3d_box_label <- function(variable, value, digits = NULL, angle_unit = "degrees",
                               diform_spec = NULL) {
  unit <- diform3d_spec_field(diform_spec, "variable_units", variable)
  if (is.null(unit) || !nzchar(unit)) {
    unit <- diform3d_default_unit(variable, angle_unit = angle_unit, for_value = TRUE)
  }
  unit <- diform3d_normalise_value_unit(unit)
  label <- diform3d_spec_field(diform_spec, "variable_labels", variable)
  name <- if (!is.null(label)) label else diform3d_pretty_name(variable)
  if (is.null(digits)) digits <- if (identical(unit, "\u00B0")) 0 else 2
  value_fmt <- sprintf(paste0("%.", digits, "f"), value)
  if (is.null(unit) || !nzchar(unit)) return(sprintf("%s = %s", name, value_fmt))
  sep <- if (identical(unit, "\u00B0")) "" else " "
  sprintf("%s = %s%s%s", name, value_fmt, sep, unit)
}

diform3d_physical <- function(cspd, hl, ht, meta = NULL) {
  response <- hl^2 + ht^2
  angle <- atan2(ht, hl) * 180 / pi
  if (!is.null(meta) && identical(meta$angle_unit, "radians")) {
    angle <- ((atan2(ht, hl) - meta$origin) %% meta$period) + meta$origin
  } else if (!is.null(meta)) {
    angle <- ((angle - meta$origin) %% meta$period) + meta$origin
  } else {
    angle <- angle %% 360
  }
  list(Cspd = cspd, angularDifference = angle, response = response)
}

diform3d_obs <- function(obs_data, excess_df = NULL, response = "hs", diform_spec = NULL) {
  if (is.null(obs_data) || !nrow(obs_data)) return(NULL)
  meta <- diform3d_polar_meta(diform_spec = diform_spec)
  obs <- data.frame(
    Cspd = obs_data[[meta$linear]],
    angularDifference = ((obs_data[[meta$angle]] - meta$origin) %% meta$period) + meta$origin,
    response = obs_data[[meta$response]]
  )
  obs$is_peak <- FALSE
  if (!is.null(excess_df) && nrow(excess_df) &&
      "time" %in% names(obs_data) && "time" %in% names(excess_df)) {
    obs$is_peak <- obs_data$time %in% excess_df$time
  }
  obs <- obs[complete.cases(obs[, c("Cspd", "angularDifference", "response")]), ]
  rownames(obs) <- NULL
  obs
}

diform3d_obs_working <- function(obs) {
  a <- obs$angularDifference * pi / 180
  h <- sqrt(pmax(obs$response, 0))
  cbind(Cspd = obs$Cspd, H_L = h * cos(a), H_T = h * sin(a))
}

diform3d_signed <- function(pts, constraints) {
  as.matrix(pts) %*% t(constraints$normal) +
    matrix(constraints$offset, nrow = nrow(pts),
           ncol = length(constraints$offset), byrow = TRUE)
}

diform3d_outside <- function(obs, constraints, tol = 1e-10) {
  if (is.null(obs) || !nrow(obs)) return(logical())
  pts <- diform3d_obs_working(obs)
  apply(diform3d_signed(pts, constraints) > tol, 1, any)
}

diform3d_halfspace_mask <- function(obs, support_info, tol = 1e-10) {
  if (is.null(obs) || !nrow(obs)) return(logical())
  pts <- diform3d_obs_working(obs)
  as.numeric(pts %*% support_info$normal + support_info$offset) > tol
}

diform3d_highlight_masks <- function(obs, constraints, support_info) {
  if (is.null(obs) || !nrow(obs)) return(list(orange = logical(), green = logical()))
  outside <- diform3d_outside(obs, constraints)
  halfspace <- diform3d_halfspace_mask(obs, support_info)
  peaks <- isTRUE(any(obs$is_peak)) & obs$is_peak
  list(
    orange = peaks & outside & !halfspace,
    green = peaks & outside & halfspace
  )
}

diform3d_peak_mask <- function(masks) {
  orange <- masks$orange
  green <- masks$green
  if (is.null(orange)) orange <- logical()
  if (is.null(green)) green <- logical()
  if (!length(orange)) return(green)
  if (!length(green)) return(orange)
  orange | green
}

diform3d_limits <- function(all_vertices, obs = NULL, z_max = NULL, working_z_max = NULL) {
  c_all <- all_vertices[, 1]
  hl_all <- all_vertices[, 2]
  ht_all <- all_vertices[, 3]
  hs_all <- hl_all^2 + ht_all^2

  if (!is.null(obs) && nrow(obs)) {
    ow <- diform3d_obs_working(obs)
    c_all <- c(c_all, ow[, "Cspd"])
    hl_all <- c(hl_all, ow[, "H_L"])
    ht_all <- c(ht_all, ow[, "H_T"])
    hs_all <- c(hs_all, obs$response)
  }

  a <- max(abs(c(hl_all, ht_all)), na.rm = TRUE) * 1.08
  response_upper <- diform3d_axis_max(z_max, max(hs_all, na.rm = TRUE) * 1.08)
  working_y <- c(-a, a)
  working_z <- c(-a, diform3d_axis_max(working_z_max, a))

  list(
    cspd = c(0, max(c_all, na.rm = TRUE) * 1.08),
    angle = c(0, 360),
    response = c(0, response_upper),
    working = working_y,
    working_y = working_y,
    working_z = working_z
  )
}

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

diform3d_unit <- function(x) {
  n <- sqrt(sum(x^2))
  if (!is.finite(n) || n == 0) return(x)
  x / n
}

diform3d_cross <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

diform3d_plane_basis <- function(normal) {
  n <- diform3d_unit(normal)
  a <- if (abs(n[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  e1 <- a - sum(a * n) * n
  e1 <- diform3d_unit(e1)
  e2 <- diform3d_unit(diform3d_cross(n, e1))
  list(e1 = e1, e2 = e2)
}

diform3d_direction_key <- function(x) {
  apply(round(as.matrix(x), 10), 1, paste, collapse = ":")
}

diform3d_direction_mesh_for_m <- function(m) {
  base <- do.call(rbind, lapply(0:m, function(i) {
    do.call(rbind, lapply(0:(m - i), function(j) c(i, j, m - i - j)))
  }))
  base <- base / m
  base <- base[rowSums(base^2) > 1e-12, , drop = FALSE]
  signs <- as.matrix(expand.grid(c(1, -1), c(1, -1), c(1, -1)))
  dirs_all <- do.call(rbind, lapply(seq_len(nrow(signs)), function(i) {
    sweep(base, 2, signs[i, ], "*")
  }))
  dirs_all <- round(dirs_all / sqrt(rowSums(dirs_all^2)), 10)
  dirs <- unique(dirs_all)
  colnames(dirs) <- c("d1", "d2", "d3")
  key_to_index <- setNames(seq_len(nrow(dirs)), diform3d_direction_key(dirs))

  point_key <- function(i, j, sign) {
    w <- c(i, j, m - i - j) / m
    w <- sign * w
    w <- round(w / sqrt(sum(w^2)), 10)
    paste(w, collapse = ":")
  }

  tris <- list()
  k <- 0L
  for (s in seq_len(nrow(signs))) {
    sign <- signs[s, ]
    for (i in 0:(m - 1)) {
      for (j in 0:(m - i - 1)) {
        tri1 <- c(
          key_to_index[[point_key(i, j, sign)]],
          key_to_index[[point_key(i + 1, j, sign)]],
          key_to_index[[point_key(i, j + 1, sign)]]
        )
        if (length(unique(tri1)) == 3) {
          k <- k + 1L
          tris[[k]] <- tri1
        }
        if (i + j <= m - 2) {
          tri2 <- c(
            key_to_index[[point_key(i + 1, j, sign)]],
            key_to_index[[point_key(i + 1, j + 1, sign)]],
            key_to_index[[point_key(i, j + 1, sign)]]
          )
          if (length(unique(tri2)) == 3) {
            k <- k + 1L
            tris[[k]] <- tri2
          }
        }
      }
    }
  }

  triangles <- do.call(rbind, tris)
  triangles <- unique(t(apply(triangles, 1, sort)))
  list(dirs = dirs, triangles = triangles)
}

diform3d_direction_triangles <- function(dirs, max_m = 80) {
  dirs <- round(as.matrix(dirs), 10)
  input_keys <- diform3d_direction_key(dirs)
  for (m in seq_len(max_m)) {
    mesh <- diform3d_direction_mesh_for_m(m)
    if (nrow(mesh$dirs) != nrow(dirs)) next
    mesh_keys <- diform3d_direction_key(mesh$dirs)
    if (!setequal(input_keys, mesh_keys)) next
    remap <- match(mesh_keys, input_keys)
    triangles <- matrix(remap[mesh$triangles], ncol = 3)
    triangles <- triangles[complete.cases(triangles), , drop = FALSE]
    if (nrow(triangles)) return(triangles)
  }
  NULL
}

diform3d_support_info <- function(vertices, constraints) {
  hs <- vertices[, 2]^2 + vertices[, 3]^2
  i_star <- which.max(hs)
  x_star <- vertices[i_star, ]
  interior <- colMeans(vertices)
  grad <- diform3d_unit(c(0, 2 * x_star[2], 2 * x_star[3]))

  triangles <- NULL
  if (!is.null(constraints$dirs) && nrow(vertices) == nrow(constraints$dirs)) {
    triangles <- diform3d_direction_triangles(constraints$dirs)
  }

  if (!is.null(triangles) && nrow(triangles)) {
    facets <- triangles[apply(triangles == i_star, 1, any), , drop = FALSE]
    best_score <- -Inf
    best_normal <- NULL
    best_offset <- NULL

    for (j in seq_len(nrow(facets))) {
      pts <- vertices[facets[j, ], , drop = FALSE]
      normal <- diform3d_cross(pts[2, ] - pts[1, ], pts[3, ] - pts[1, ])
      if (sqrt(sum(normal^2)) <= 1e-12) next
      offset <- -sum(normal * pts[1, ])
      if (sum(normal * interior) + offset > 0) {
        normal <- -normal
        offset <- -offset
      }
      score <- sum(diform3d_unit(normal) * grad)
      if (is.finite(score) && score > best_score) {
        best_score <- score
        best_normal <- normal
        best_offset <- offset
      }
    }

    if (!is.null(best_normal)) {
      phys <- diform3d_physical(x_star[1], x_star[2], x_star[3])
      return(list(
        x_star = as.numeric(x_star),
        hs_star = phys$response,
        angle_star = phys$angularDifference,
        normal = diform3d_unit(best_normal),
        offset = best_offset / sqrt(sum(best_normal^2)),
        source = "surface_facet"
      ))
    }
  }

  signed <- as.numeric(x_star %*% t(constraints$normal) + constraints$offset)
  scale <- pmax(1, abs(constraints$rv))
  active <- which(abs(signed) <= 1e-6 * scale)
  if (!length(active)) active <- which.max(signed)
  score <- apply(constraints$normal[active, , drop = FALSE], 1, function(z) {
    sum(diform3d_unit(z) * grad)
  })
  i <- active[which.max(score)]
  phys <- diform3d_physical(x_star[1], x_star[2], x_star[3])
  list(
    x_star = as.numeric(x_star),
    hs_star = phys$response,
    angle_star = phys$angularDifference,
    normal = as.numeric(constraints$normal[i, ]),
    offset = constraints$offset[i],
    source = "active_constraint"
  )
}

diform3d_plane_patch <- function(support_info, vertices, scale = 0.32, n = 16) {
  ranges <- apply(vertices, 2, range)
  radius <- scale * max(ranges[2, ] - ranges[1, ])
  basis <- diform3d_plane_basis(support_info$normal)
  u <- seq(-radius, radius, length.out = n)
  v <- seq(-radius, radius, length.out = n)
  grid <- expand.grid(u = u, v = v)
  pts <- t(apply(grid, 1, function(g) {
    support_info$x_star + g[1] * basis$e1 + g[2] * basis$e2
  }))
  list(
    X = matrix(pts[, 1], n, n),
    Y = matrix(pts[, 2], n, n),
    Z = matrix(pts[, 3], n, n)
  )
}

diform3d_halfspace_physical_layers <- function(support_info, limits,
                                               n_cspd = 80, n_ang = 180,
                                               n_layers = 7) {
  cg <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  ag <- seq(limits$angle[1], limits$angle[2], length.out = n_ang)
  C <- matrix(cg, nrow = length(cg), ncol = length(ag))
  A <- matrix(ag, nrow = length(cg), ncol = length(ag), byrow = TRUE)
  n <- support_info$normal
  alpha <- n[2] * cos(A * pi / 180) + n[3] * sin(A * pi / 180)
  s0 <- n[1] * C + support_info$offset
  rmax <- sqrt(limits$response[2])
  r_low <- matrix(NA_real_, nrow(C), ncol(C))
  r_up <- matrix(NA_real_, nrow(C), ncol(C))

  pos <- alpha > 1e-12
  neg <- alpha < -1e-12
  flat <- abs(alpha) <= 1e-12
  r_low[pos] <- pmax(0, -s0[pos] / alpha[pos])
  r_up[pos] <- rmax
  r_low[neg] <- 0
  r_up[neg] <- -s0[neg] / alpha[neg]
  r_low[flat] <- 0
  r_up[flat] <- ifelse(s0[flat] >= 0, rmax, NA_real_)

  r_low <- matrix(pmax(0, pmin(r_low, rmax)), nrow = nrow(C), ncol = ncol(C))
  r_up <- matrix(pmax(0, pmin(r_up, rmax)), nrow = nrow(C), ncol = ncol(C))
  valid <- is.finite(r_low) & is.finite(r_up) & r_up >= r_low + 1e-10

  lapply(seq(0.10, 0.95, length.out = n_layers), function(t) {
    R <- r_low + t * (r_up - r_low)
    H <- R^2
    H[!valid] <- NA_real_
    list(X = C, Y = A, Z = H)
  })
}

diform3d_halfspace_physical_boundary <- function(support_info, limits,
                                                 n_cspd = 80, n_ang = 180) {
  cg <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  ag <- seq(limits$angle[1], limits$angle[2], length.out = n_ang)
  C <- matrix(cg, nrow = length(cg), ncol = length(ag))
  A <- matrix(ag, nrow = length(cg), ncol = length(ag), byrow = TRUE)
  n <- support_info$normal

  alpha <- n[2] * cos(A * pi / 180) + n[3] * sin(A * pi / 180)
  beta <- n[1] * C + support_info$offset
  r <- -beta / alpha
  h <- r^2
  h[abs(alpha) <= 1e-12 | r < 0 | h < limits$response[1] | h > limits$response[2]] <- NA_real_

  list(X = C, Y = A, Z = h)
}

diform3d_halfspace_physical_boundary_local <- function(support_info, limits,
                                                       n_cspd = 70, n_ang = 90,
                                                       cspd_fraction = 0.32,
                                                       angle_half_width = 34) {
  star <- diform3d_physical(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3]
  )
  c_half <- cspd_fraction * diff(limits$cspd)
  cg <- seq(
    max(limits$cspd[1], star$Cspd - c_half),
    min(limits$cspd[2], star$Cspd + c_half),
    length.out = n_cspd
  )

  ag <- seq(star$angularDifference - angle_half_width,
            star$angularDifference + angle_half_width,
            length.out = n_ang)
  ag <- ag[ag >= limits$angle[1] & ag <= limits$angle[2]]
  if (length(ag) < 4) {
    ag <- seq(limits$angle[1], limits$angle[2], length.out = n_ang)
  }

  C <- matrix(cg, nrow = length(cg), ncol = length(ag))
  A <- matrix(ag, nrow = length(cg), ncol = length(ag), byrow = TRUE)
  n <- support_info$normal

  alpha <- n[2] * cos(A * pi / 180) + n[3] * sin(A * pi / 180)
  beta <- n[1] * C + support_info$offset
  r <- -beta / alpha
  h <- r^2
  h[abs(alpha) <= 1e-12 | r < 0 | h < limits$response[1] | h > limits$response[2]] <- NA_real_

  list(X = C, Y = A, Z = h)
}

diform3d_project_grid <- function(pm, X, Y, Z, col, border = NA, zlim = NULL,
                                  max_z_span = NULL) {
  nr <- nrow(X)
  nc <- ncol(X)
  for (i in seq_len(nr - 1)) {
    for (j in seq_len(nc - 1)) {
      zz <- c(Z[i, j], Z[i + 1, j], Z[i + 1, j + 1], Z[i, j + 1])
      if (any(!is.finite(zz))) next
      if (!is.null(zlim) && (max(zz) < zlim[1] || min(zz) > zlim[2])) next
      if (!is.null(max_z_span) && diff(range(zz)) > max_z_span) next
      xy <- grDevices::trans3d(
        c(X[i, j], X[i + 1, j], X[i + 1, j + 1], X[i, j + 1]),
        c(Y[i, j], Y[i + 1, j], Y[i + 1, j + 1], Y[i, j + 1]),
        zz,
        pm
      )
      graphics::polygon(xy$x, xy$y, col = col, border = border)
    }
  }
}

diform3d_draw_colorbar <- function(vmin, vmax, palette, label) {
  old <- par(mar = c(5.4, 0.4, 4.6, 4.2), xaxs = "i", yaxs = "i")
  on.exit(par(old), add = TRUE)
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) {
    vmax <- vmin + 1
  }
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(vmin, vmax))
  y <- seq(vmin, vmax, length.out = length(palette) + 1L)
  graphics::rect(0.18, y[-length(y)], 0.72, y[-1],
                 col = palette, border = NA)
  graphics::rect(0.18, vmin, 0.72, vmax, border = "black", lwd = 0.9)
  ticks <- pretty(c(vmin, vmax), n = 6)
  ticks <- ticks[is.finite(ticks) & ticks >= vmin & ticks <= vmax]
  graphics::axis(4, at = ticks, pos = 0.78, las = 1, cex.axis = 1.20)
  graphics::mtext(label, side = 4, line = 2.7, cex = 1.35)
}

diform3d_polar_projector <- function(r_max, z_max, elev = 24, azim = -42) {
  elev_rad <- elev * pi / 180
  azim_rad <- azim * pi / 180
  z_scale <- if (is.finite(z_max) && z_max > 0) 1.30 * r_max / z_max else 1

  function(x, y, z) {
    xr <- x * cos(azim_rad) - y * sin(azim_rad)
    yr <- x * sin(azim_rad) + y * cos(azim_rad)
    zr <- z * z_scale
    list(
      x = xr,
      y = yr * cos(elev_rad) + zr * sin(elev_rad),
      depth = yr * sin(elev_rad) - zr * cos(elev_rad)
    )
  }
}

diform3d_polar_segments <- function(project, x, y, z, col, lwd = 1) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 2) return(invisible(NULL))
  if (length(col) == 1L) col <- rep(col, length(x) - 1L)

  p <- project(x, y, z)
  ok_seg <- ok[-length(ok)] & ok[-1]
  if (!any(ok_seg)) return(invisible(NULL))

  graphics::segments(
    p$x[-length(p$x)][ok_seg],
    p$y[-length(p$y)][ok_seg],
    p$x[-1][ok_seg],
    p$y[-1][ok_seg],
    col = col[ok_seg],
    lwd = lwd
  )
  invisible(NULL)
}

diform3d_polar_angle_rad <- function(angle_deg) {
  (90 - angle_deg) * pi / 180
}

diform3d_polar_project_grid <- function(project, X, Y, Z, col, border = NA,
                                        zlim = NULL, max_z_span = NULL) {
  nr <- nrow(X)
  nc <- ncol(X)
  if (is.null(dim(col))) {
    col <- matrix(col, nrow = nr - 1, ncol = nc - 1)
  }
  cells <- vector("list", max(0, (nr - 1) * (nc - 1)))
  k <- 0L

  for (i in seq_len(nr - 1)) {
    for (j in seq_len(nc - 1)) {
      zz <- c(Z[i, j], Z[i + 1, j], Z[i + 1, j + 1], Z[i, j + 1])
      if (any(!is.finite(zz))) next
      if (!is.null(zlim) && (max(zz) < zlim[1] || min(zz) > zlim[2])) next
      if (!is.null(max_z_span) && diff(range(zz)) > max_z_span) next

      p <- project(
        c(X[i, j], X[i + 1, j], X[i + 1, j + 1], X[i, j + 1]),
        c(Y[i, j], Y[i + 1, j], Y[i + 1, j + 1], Y[i, j + 1]),
        zz
      )

      k <- k + 1L
      cells[[k]] <- list(
        x = p$x,
        y = p$y,
        depth = mean(p$depth, na.rm = TRUE),
        col = col[i, j]
      )
    }
  }

  if (!k) return(invisible(NULL))
  cells <- cells[seq_len(k)]
  depth <- vapply(cells, `[[`, numeric(1), "depth")
  cells <- cells[order(depth)]

  for (cell in cells) {
    graphics::polygon(cell$x, cell$y, col = cell$col, border = border)
  }
  invisible(NULL)
}

diform3d_text_box <- function(x, y, labels, cex = 0.82) {
  usr <- par("usr")
  line_h <- strheight("M", cex = cex) * 1.35
  widths <- strwidth(labels, cex = cex)
  pad_x <- diff(usr[1:2]) * 0.012
  pad_y <- diff(usr[3:4]) * 0.014
  h <- line_h * length(labels) + 2 * pad_y
  w <- max(widths) + 2 * pad_x
  graphics::rect(x - pad_x, y - h + pad_y, x + w - pad_x, y + pad_y,
                 col = diform3d_alpha("white", 0.90), border = "lightgray")
  yy <- y - seq_along(labels) * line_h + line_h * 0.45
  graphics::text(x, yy, labels = labels, adj = c(0, 0.5), cex = cex)
}

diform3d_add_annotation <- function(lines) {
  usr <- par("usr")
  x <- usr[1] - 0.090 * diff(usr[1:2])
  y <- usr[4] + 0.025 * diff(usr[3:4])
  diform3d_text_box(x, y, lines, cex = 1.25)
}

diform3d_subsample_obs <- function(obs, max_background_points = 0, seed = 123) {
  if (is.null(obs) || !nrow(obs)) return(obs)
  if (is.null(max_background_points) || max_background_points <= 0 ||
      nrow(obs) <= max_background_points) {
    return(obs)
  }
  set.seed(seed)
  obs[sort(sample.int(nrow(obs), max_background_points)), , drop = FALSE]
}

diform3d_draw_points_physical <- function(pm, obs, masks, limits, max_labels = 6,
                                          show_peak_highlights = TRUE,
                                          obs_background = obs) {
  if (is.null(obs) || !nrow(obs)) return(invisible(NULL))
  if (!is.null(obs_background) && nrow(obs_background)) {
    xyz <- grDevices::trans3d(obs_background$Cspd, obs_background$angularDifference,
                              obs_background$response, pm)
    graphics::points(xyz, pch = 16, cex = 0.17,
                     col = diform3d_alpha(diform3d_obs_color, 0.040))
  }
  if (!isTRUE(show_peak_highlights)) return(invisible(NULL))

  draw_mask <- function(mask, bg, cex) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    p <- grDevices::trans3d(obs$Cspd[mask], obs$angularDifference[mask],
                            obs$response[mask], pm)
    graphics::points(p, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }
  peak_mask <- diform3d_peak_mask(masks)
  draw_mask(peak_mask, diform3d_orange_point, 0.90)

  label_mask <- function(mask, col) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    idx <- which(mask)
    idx <- idx[order(obs$response[idx], decreasing = TRUE)]
    idx <- idx[seq_len(min(max_labels, length(idx)))]
    p <- grDevices::trans3d(obs$Cspd[idx], obs$angularDifference[idx],
                            obs$response[idx] + 0.012 * limits$response[2], pm)
    graphics::text(p$x, p$y, sprintf("Hs=%.2f", obs$response[idx]),
                   cex = 0.70, col = col, pos = 3)
  }
  label_mask(peak_mask, "darkorange")
  invisible(NULL)
}

diform3d_draw_points_working <- function(pm, obs, masks, limits, max_labels = 6,
                                         show_peak_highlights = TRUE,
                                         obs_background = obs) {
  if (is.null(obs) || !nrow(obs)) return(invisible(NULL))
  pts <- diform3d_obs_working(obs)
  if (!is.null(obs_background) && nrow(obs_background)) {
    bg_pts <- diform3d_obs_working(obs_background)
    xyz <- grDevices::trans3d(bg_pts[, 1], bg_pts[, 2], bg_pts[, 3], pm)
    graphics::points(xyz, pch = 16, cex = 0.17,
                     col = diform3d_alpha(diform3d_obs_color, 0.040))
  }
  if (!isTRUE(show_peak_highlights)) return(invisible(NULL))

  draw_mask <- function(mask, bg, cex) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    p <- grDevices::trans3d(pts[mask, 1], pts[mask, 2], pts[mask, 3], pm)
    graphics::points(p, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }
  peak_mask <- diform3d_peak_mask(masks)
  draw_mask(peak_mask, diform3d_orange_point, 0.90)

  label_mask <- function(mask, col) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    idx <- which(mask)
    hs <- pts[, 2]^2 + pts[, 3]^2
    idx <- idx[order(hs[idx], decreasing = TRUE)]
    idx <- idx[seq_len(min(max_labels, length(idx)))]
    p <- grDevices::trans3d(
      pts[idx, 1], pts[idx, 2],
      pts[idx, 3] + 0.018 * diff(limits$working_z), pm
    )
    graphics::text(p$x, p$y, sprintf("Hs=%.2f", hs[idx]),
                   cex = 0.70, col = col, pos = 3)
  }
  label_mask(peak_mask, "darkorange")
  invisible(NULL)
}

diform3d_add_legend <- function(show_support_plane = TRUE) {
  legend <- c("Original data", "Peaks outside contour")
  pch <- c(16, 21)
  pt_bg <- c(NA, diform3d_orange_point)
  cols <- c(diform3d_alpha(diform3d_obs_color, 0.62), "black")

  if (isTRUE(show_support_plane)) {
    legend <- c(legend, "Tangent half-space boundary")
    pch <- c(pch, 15)
    pt_bg <- c(pt_bg, diform3d_alpha(diform3d_tangent_color, 0.32))
    cols <- c(cols, diform3d_alpha(diform3d_tangent_color, 0.90))
  }

  graphics::legend(
    "bottomright",
    inset = c(0.015, 0.020),
    legend = legend,
    pch = pch,
    pt.bg = pt_bg,
    col = cols,
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.94,
    y.intersp = 1.14
  )
}

diform3d_plot_physical <- function(path, vertices, constraints, support_info,
                                   all_vertices, obs, masks,
                                   response_label, response_short, response_unit,
                                   n_years, sample_label, rp, elev, azim,
                                   n_ang, n_cspd, show_points,
                                   show_peak_highlights, show_support_plane, plane_scale,
                                   max_background_points, max_labels,
                                   z_max = NULL, meta = NULL) {
  if (is.null(meta)) meta <- list(linear = "Cspd", angle = "angularDifference")
  palette <- diform3d_palette(256)
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  vmax <- max(all_vertices[, 2]^2 + all_vertices[, 3]^2, na.rm = TRUE) * 1.03
  c_grid <- seq(min(vertices[, 1]), max(vertices[, 1]), length.out = n_cspd)
  a_grid <- seq(0, 360, length.out = n_ang)
  hs <- diform3d_surface_physical(constraints, c_grid, a_grid)
  cols <- diform3d_surface_colors(hs, 0, vmax, palette, alpha = 0.76)

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(8.3, 1.3))
  old <- par(mar = c(5.2, 5.4, 5.0, 1.1), xpd = NA)
  on.exit(par(old), add = TRUE)

  pm <- graphics::persp(
    x = c_grid, y = a_grid, z = hs,
    xlim = limits$cspd,
    ylim = limits$angle,
    zlim = limits$response,
    col = cols,
    border = NA,
    shade = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = diform3d_label_for(meta$linear, meta$diform_spec),
    ylab = diform3d_label_for(meta$angle, meta$diform_spec,
                              angle_unit = meta$angle_unit),
    zlab = response_label,
    main = sprintf("D-IFORM %s-yr contour\n%s",
                   format(rp, trim = TRUE),
                   sample_label),
    cex.lab = 1.34,
    cex.axis = 1.06,
    cex.main = 1.42,
    font.main = 2
  )

  if (isTRUE(show_support_plane)) {
    boundary <- diform3d_halfspace_physical_boundary_local(
      support_info,
      limits,
      n_cspd = max(45, floor(n_cspd * 0.80)),
      n_ang = max(55, floor(n_ang * 0.55))
    )
    diform3d_project_grid(
      pm,
      boundary$X, boundary$Y, boundary$Z,
      col = diform3d_alpha(diform3d_tangent_color, 0.34),
      border = diform3d_alpha(diform3d_tangent_color, 0.24),
      zlim = limits$response,
      max_z_span = diff(limits$response) * 0.25
    )
    star <- diform3d_physical(support_info$x_star[1], support_info$x_star[2],
                              support_info$x_star[3], meta = meta)
    p <- grDevices::trans3d(star$Cspd, star$angularDifference, star$response, pm)
    graphics::points(p, pch = 18, col = diform3d_tangent_color, cex = 1.08)
  }

  obs_bg <- diform3d_subsample_obs(obs, max_background_points)
  if (isTRUE(show_points)) {
    diform3d_draw_points_physical(pm, obs, masks, limits, max_labels,
                                  show_peak_highlights = show_peak_highlights,
                                  obs_background = obs_bg)
  }

  phys_star <- diform3d_physical(vertices[, 1], vertices[, 2], vertices[, 3], meta = meta)
  i <- which.max(phys_star$response)
  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""
  diform3d_add_annotation(c(
    paste0("Max ", response_short, " on contour"),
    sprintf("%s = %.3f%s", response_short, phys_star$response[i], unit),
    diform3d_box_label(meta$linear, phys_star$Cspd[i],
                       diform_spec = meta$diform_spec),
    diform3d_box_label(meta$angle, phys_star$angularDifference[i],
                       angle_unit = meta$angle_unit,
                       diform_spec = meta$diform_spec)
  ))

  if (isTRUE(show_points) && isTRUE(show_peak_highlights) && !is.null(obs) && nrow(obs)) {
    diform3d_add_legend(show_support_plane)
  }
  diform3d_draw_colorbar(0, vmax, palette, response_label)
  invisible(path)
}

diform3d_plot_polar_physical <- function(path, vertices, constraints, support_info,
                                         all_vertices, obs, masks,
                                         response_label, response_short, response_unit,
                                         n_years, rp, elev, azim,
                                         n_ang, n_cspd, show_points,
                                         max_background_points, max_labels,
                                         z_max = NULL) {
  palette <- diform3d_palette(256)
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  if (!is.null(obs) && nrow(obs) && any(is.finite(obs$Cspd))) {
    limits$cspd <- c(0, max(obs$Cspd, na.rm = TRUE))
  }
  vmax <- max(all_vertices[, "H_L"]^2 + all_vertices[, "H_T"]^2, na.rm = TRUE) * 1.03
  if (!is.finite(vmax) || vmax <= 0) vmax <- limits$response[2]

  c_grid <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  a_grid <- seq(0, 360, length.out = n_ang)
  hs <- diform3d_surface_physical(constraints, c_grid, a_grid)
  hs[hs < limits$response[1] | hs > limits$response[2]] <- NA_real_

  C <- matrix(c_grid, nrow = length(c_grid), ncol = length(a_grid))
  A <- matrix(a_grid, nrow = length(c_grid), ncol = length(a_grid), byrow = TRUE)
  A_display <- diform3d_polar_angle_rad(A)
  X <- C * cos(A_display)
  Y <- C * sin(A_display)
  Z <- hs

  r_max <- max(limits$cspd[2], 1)
  z_top <- limits$response[2]
  project <- diform3d_polar_projector(r_max, z_top, elev = elev, azim = azim)

  ring_ticks <- pretty(c(0, limits$cspd[2]), n = 6)
  ring_ticks <- ring_ticks[is.finite(ring_ticks) & ring_ticks >= 0 & ring_ticks <= limits$cspd[2]]
  ring_ticks <- sort(unique(c(0, ring_ticks, limits$cspd[2])))
  angle_ticks <- seq(0, 315, by = 45)
  ring_theta <- seq(0, 2 * pi, length.out = 241)

  label_r <- r_max * 1.12

  base_x <- c(
    as.vector(outer(ring_ticks, cos(ring_theta))),
    label_r * cos(diform3d_polar_angle_rad(angle_ticks)),
    as.vector(X)
  )
  base_y <- c(
    as.vector(outer(ring_ticks, sin(ring_theta))),
    label_r * sin(diform3d_polar_angle_rad(angle_ticks)),
    as.vector(Y)
  )
  base_z <- c(
    rep(0, length(ring_ticks) * length(ring_theta)),
    rep(0, length(angle_ticks)),
    as.vector(Z)
  )

  if (!is.null(obs) && nrow(obs)) {
    a_obs <- diform3d_polar_angle_rad(obs$angularDifference)
    base_x <- c(base_x, obs$Cspd * cos(a_obs))
    base_y <- c(base_y, obs$Cspd * sin(a_obs))
    base_z <- c(base_z, obs$response)
  }

  p_all <- project(base_x, base_y, base_z)
  xlim <- range(p_all$x[is.finite(p_all$x)], na.rm = TRUE)
  ylim <- range(p_all$y[is.finite(p_all$y)], na.rm = TRUE)
  xlim <- xlim + c(-1, 1) * 0.10 * diff(xlim)
  ylim <- ylim + c(-1, 1) * 0.10 * diff(ylim)

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(8.3, 1.3))
  old <- par(mar = c(5.2, 5.4, 5.0, 1.1), xpd = NA)
  on.exit(par(old), add = TRUE)

  graphics::plot(
    0, 0,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = xlim,
    ylim = ylim,
    asp = 1,
    main = sprintf("D-IFORM %s-yr contour in polar physical space\nbase radius = Cspd, height = %s",
                   format(rp, trim = TRUE), response_short),
    cex.main = 1.12,
    font.main = 2
  )

  for (cs in ring_ticks) {
    p <- project(cs * cos(ring_theta), cs * sin(ring_theta), rep(0, length(ring_theta)))
    graphics::lines(p$x, p$y, col = "grey78", lty = 3)
  }

  for (ang in angle_ticks) {
    rr <- seq(0, r_max, length.out = 80)
    aa <- diform3d_polar_angle_rad(ang)
    p <- project(rr * cos(aa), rr * sin(aa), rep(0, length(rr)))
    graphics::lines(p$x, p$y, col = "grey83", lty = 3)
  }

  label_ticks <- ring_ticks[ring_ticks > 0]
  label_ticks <- label_ticks[seq_len(min(length(label_ticks), 5))]

  cols <- diform3d_surface_colors(Z, 0, vmax, palette, alpha = 0.82)
  diform3d_polar_project_grid(
    project,
    X, Y, Z,
    col = cols,
    border = NA,
    zlim = limits$response
  )

  # Sparse guide lines, closer to the classic contour figures than a full wireframe.
  surface_col <- function(v, alpha = 0.70) {
    diform3d_color_values(v, 0, vmax, palette, alpha = alpha)
  }
  guide_angle <- seq(0, 315, by = 45)
  guide_idx <- unique(pmax(1, pmin(ncol(X), round(guide_angle / 360 * (ncol(X) - 1)) + 1)))
  for (j in guide_idx) {
    seg_hs <- rowMeans(cbind(Z[-nrow(Z), j], Z[-1, j]), na.rm = FALSE)
    diform3d_polar_segments(
      project,
      X[, j], Y[, j], Z[, j],
      col = diform3d_alpha("grey20", 0.22),
      lwd = 0.8
    )
  }

  rim_idx <- which.max(c_grid)
  seg_hs <- colMeans(rbind(Z[rim_idx, -ncol(Z)], Z[rim_idx, -1]), na.rm = FALSE)
  diform3d_polar_segments(
    project,
    X[rim_idx, ], Y[rim_idx, ], Z[rim_idx, ],
    col = surface_col(seg_hs, alpha = 0.92),
    lwd = 1.25
  )

  for (ang in angle_ticks) {
    aa <- diform3d_polar_angle_rad(ang)
    txt <- project(label_r * cos(aa), label_r * sin(aa), 0)
    graphics::text(txt$x, txt$y, paste0(ifelse(ang == 0, "0/360", ang), " deg"),
                   cex = 0.68, col = "grey25")
  }

  for (cs in label_ticks) {
    aa <- 20 * pi / 180
    txt <- project(cs * cos(aa), cs * sin(aa), 0)
    label <- if (abs(cs - limits$cspd[2]) < 1e-8 && abs(cs - round(cs)) > 1e-8) {
      sprintf("%.2f m.s^-1", cs)
    } else {
      sprintf("%.1f m.s^-1", cs)
    }
    graphics::text(txt$x, txt$y, label, cex = 0.62, col = "grey35")
  }

  obs_bg <- diform3d_subsample_obs(obs, max_background_points)
  if (isTRUE(show_points) && !is.null(obs_bg) && nrow(obs_bg)) {
    a <- diform3d_polar_angle_rad(obs_bg$angularDifference)
    p <- project(obs_bg$Cspd * cos(a), obs_bg$Cspd * sin(a), obs_bg$response)
    graphics::points(p$x, p$y, pch = 16, cex = 0.15,
                     col = diform3d_alpha(diform3d_obs_color, 0.10))
  }

  draw_mask <- function(mask, bg, cex) {
    if (is.null(obs) || !nrow(obs) || !length(mask) || !any(mask)) return(invisible(NULL))
    a <- diform3d_polar_angle_rad(obs$angularDifference[mask])
    p <- project(obs$Cspd[mask] * cos(a), obs$Cspd[mask] * sin(a), obs$response[mask])
    graphics::points(p$x, p$y, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }

  if (isTRUE(show_points)) {
    draw_mask(diform3d_peak_mask(masks), diform3d_orange_point, 0.88)
  }

  phys_star <- diform3d_physical(vertices[, "Cspd"], vertices[, "H_L"], vertices[, "H_T"])
  i <- which.max(phys_star$response)
  star_a <- diform3d_polar_angle_rad(phys_star$angularDifference[i])
  star <- project(
    phys_star$Cspd[i] * cos(star_a),
    phys_star$Cspd[i] * sin(star_a),
    phys_star$response[i]
  )
  graphics::points(star$x, star$y, pch = 18, cex = 1.75,
                   col = diform3d_tangent_color)
  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""
  diform3d_add_annotation(c(
    paste0("Max ", response_short, " on contour"),
    sprintf("%s = %.3f%s", response_short, phys_star$response[i], unit),
    diform3d_box_label("Cspd", phys_star$Cspd[i]),
    diform3d_box_label("angularDifference", phys_star$angularDifference[i])
  ))

  graphics::legend(
    "topright",
    inset = 0.01,
    legend = c(
      "Original data",
      "Contour surface, coloured by Hs",
      paste0("Max ", response_short, " on contour"),
      "Peaks outside contour"
    ),
    pch = c(16, NA, 18, 21),
    lty = c(NA, 1, NA, NA),
    pt.bg = c(NA, NA, NA, diform3d_orange_point),
    col = c(diform3d_alpha(diform3d_obs_color, 0.62),
            palette[210], diform3d_tangent_color, "black"),
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.70,
    y.intersp = 1.10
  )

  usr <- par("usr")
  graphics::text(
    usr[1] + 0.01 * diff(usr[1:2]),
    usr[3] + 0.02 * diff(usr[3:4]),
    labels = "Physical polar view: 0 deg at top; radius = Cspd; colour = Hs.",
    adj = c(0, 0),
    cex = 0.76,
    col = "grey25"
  )

  diform3d_draw_colorbar(0, vmax, palette, response_label)
  invisible(path)
}

diform3d_plot_polar_2d <- function(path, vertices, constraints,
                                   all_vertices, obs, masks,
                                   response_label, response_short, response_unit,
                                   n_years, rp, n_ang, n_cspd,
                                   show_points, max_background_points,
                                   z_max = NULL) {
  palette <- diform3d_palette(256)
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  if (!is.null(obs) && nrow(obs) && any(is.finite(obs$Cspd))) {
    limits$cspd <- c(0, max(obs$Cspd, na.rm = TRUE))
  }

  vmax <- max(all_vertices[, "H_L"]^2 + all_vertices[, "H_T"]^2, na.rm = TRUE) * 1.03
  if (!is.finite(vmax) || vmax <= 0) vmax <- limits$response[2]

  c_grid <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  a_grid <- seq(0, 360, length.out = n_ang)
  hs <- diform3d_surface_physical(constraints, c_grid, a_grid)
  hs[hs < limits$response[1] | hs > limits$response[2]] <- NA_real_

  C <- matrix(c_grid, nrow = length(c_grid), ncol = length(a_grid))
  A <- matrix(a_grid, nrow = length(c_grid), ncol = length(a_grid), byrow = TRUE)
  A_display <- diform3d_polar_angle_rad(A)
  X <- C * cos(A_display)
  Y <- C * sin(A_display)

  cols <- diform3d_surface_colors(hs, 0, vmax, palette, alpha = 0.90)
  if (is.null(dim(cols))) cols <- matrix(cols, nrow = nrow(X) - 1, ncol = ncol(X) - 1)

  r_max <- max(limits$cspd[2], 1)
  ring_ticks <- pretty(c(0, limits$cspd[2]), n = 6)
  ring_ticks <- ring_ticks[is.finite(ring_ticks) & ring_ticks >= 0 & ring_ticks <= limits$cspd[2]]
  ring_ticks <- sort(unique(c(0, ring_ticks, limits$cspd[2])))
  angle_ticks <- seq(0, 315, by = 45)
  ring_theta <- seq(0, 2 * pi, length.out = 361)

  diform3d_open_png(path, width = 1500, height = 1200, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(8.3, 1.3))
  old <- par(mar = c(4.8, 4.8, 4.8, 1.0), xpd = NA)
  on.exit(par(old), add = TRUE)

  graphics::plot(
    0, 0,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = c(-1.17, 1.17) * r_max,
    ylim = c(-1.17, 1.17) * r_max,
    asp = 1,
    main = sprintf("D-IFORM %s-yr contour in polar physical view\nradius = Cspd, angle = direction, colour = %s",
                   format(rp, trim = TRUE), response_short),
    cex.main = 1.10,
    font.main = 2
  )

  for (i in seq_len(nrow(X) - 1)) {
    for (j in seq_len(ncol(X) - 1)) {
      zz <- c(hs[i, j], hs[i + 1, j], hs[i + 1, j + 1], hs[i, j + 1])
      if (any(!is.finite(zz))) next
      graphics::polygon(
        c(X[i, j], X[i + 1, j], X[i + 1, j + 1], X[i, j + 1]),
        c(Y[i, j], Y[i + 1, j], Y[i + 1, j + 1], Y[i, j + 1]),
        col = cols[i, j],
        border = NA
      )
    }
  }

  for (cs in ring_ticks) {
    graphics::lines(cs * cos(ring_theta), cs * sin(ring_theta),
                    col = diform3d_alpha("grey55", 0.38), lty = 3)
  }

  for (ang in angle_ticks) {
    aa <- diform3d_polar_angle_rad(ang)
    graphics::lines(c(0, r_max * cos(aa)), c(0, r_max * sin(aa)),
                    col = diform3d_alpha("grey35", 0.28), lwd = 0.8)
    graphics::text(
      1.09 * r_max * cos(aa),
      1.09 * r_max * sin(aa),
      paste0(ifelse(ang == 0, "0/360", ang), " deg"),
      cex = 0.72,
      col = "grey25"
    )
  }

  label_ticks <- ring_ticks[ring_ticks > 0]
  label_ticks <- label_ticks[seq_len(min(length(label_ticks), 5))]
  for (cs in label_ticks) {
    aa <- diform3d_polar_angle_rad(35)
    label <- if (abs(cs - limits$cspd[2]) < 1e-8 && abs(cs - round(cs)) > 1e-8) {
      sprintf("%.2f m.s^-1", cs)
    } else {
      sprintf("%.1f m.s^-1", cs)
    }
    graphics::text(cs * cos(aa), cs * sin(aa), label, cex = 0.66, col = "grey25")
  }

  obs_bg <- diform3d_subsample_obs(obs, max_background_points)
  if (isTRUE(show_points) && !is.null(obs_bg) && nrow(obs_bg)) {
    a <- diform3d_polar_angle_rad(obs_bg$angularDifference)
    graphics::points(obs_bg$Cspd * cos(a), obs_bg$Cspd * sin(a),
                     pch = 16, cex = 0.12,
                     col = diform3d_alpha(diform3d_obs_color, 0.10))
  }

  draw_mask <- function(mask, bg, cex) {
    if (is.null(obs) || !nrow(obs) || !length(mask) || !any(mask)) return(invisible(NULL))
    a <- diform3d_polar_angle_rad(obs$angularDifference[mask])
    graphics::points(obs$Cspd[mask] * cos(a), obs$Cspd[mask] * sin(a),
                     pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }

  if (isTRUE(show_points)) {
    draw_mask(diform3d_peak_mask(masks), diform3d_orange_point, 0.90)
  }

  phys_star <- diform3d_physical(vertices[, "Cspd"], vertices[, "H_L"], vertices[, "H_T"])
  i <- which.max(phys_star$response)
  star_a <- diform3d_polar_angle_rad(phys_star$angularDifference[i])
  graphics::points(phys_star$Cspd[i] * cos(star_a), phys_star$Cspd[i] * sin(star_a),
                   pch = 18, cex = 1.65, col = diform3d_tangent_color)

  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""
  diform3d_add_annotation(c(
    paste0("Max ", response_short, " on contour"),
    sprintf("%s = %.3f%s", response_short, phys_star$response[i], unit),
    diform3d_box_label("Cspd", phys_star$Cspd[i]),
    diform3d_box_label("angularDifference", phys_star$angularDifference[i])
  ))

  graphics::legend(
    "topright",
    inset = 0.01,
    legend = c(
      "Original data",
      "Contour surface, coloured by Hs",
      paste0("Max ", response_short, " on contour"),
      "Peaks outside contour"
    ),
    pch = c(16, NA, 18, 21),
    lty = c(NA, 1, NA, NA),
    pt.bg = c(NA, NA, NA, diform3d_orange_point),
    col = c(diform3d_alpha(diform3d_obs_color, 0.62),
            palette[210], diform3d_tangent_color, "black"),
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.70,
    y.intersp = 1.10
  )

  usr <- par("usr")
  graphics::text(
    usr[1] + 0.01 * diff(usr[1:2]),
    usr[3] + 0.02 * diff(usr[3:4]),
    labels = "2D polar view: 0 deg at top; radius = Cspd; colour = Hs.",
    adj = c(0, 0),
    cex = 0.76,
    col = "grey25"
  )

  diform3d_draw_colorbar(0, vmax, palette, response_label)
  invisible(path)
}

diform3d_angle_in_sector <- function(angle, start, end) {
  angle <- angle %% 360
  start <- start %% 360
  end <- end %% 360
  if (start < end) {
    angle >= start & angle < end
  } else {
    angle >= start | angle < end
  }
}

diform3d_sector_label <- function(start, end) {
  start_txt <- format(round(start %% 360), trim = TRUE)
  end_txt <- format(round(end %% 360), trim = TRUE)
  if ((start %% 360) > (end %% 360)) {
    paste0(start_txt, "-360 / 0-", end_txt, " deg")
  } else {
    paste0(start_txt, "-", end_txt, " deg")
  }
}

diform3d_finite_range <- function(...) {
  x <- unlist(list(...), use.names = FALSE)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  range(x)
}

diform3d_plot_sector_hs_cspd <- function(path, vertices, constraints,
                                         all_vertices, obs, masks,
                                         response_label, response_short,
                                         response_unit, rp, n_ang, n_cspd,
                                         show_points, max_background_points,
                                         z_max = NULL, sector_width = 30) {
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  if (!is.null(obs) && nrow(obs) && any(is.finite(obs$Cspd))) {
    limits$cspd <- c(0, max(obs$Cspd, na.rm = TRUE))
  }

  c_grid <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  a_grid <- seq(0, 360, length.out = max(n_ang, 361))
  hs <- diform3d_surface_physical(constraints, c_grid, a_grid)
  hs[hs < limits$response[1] | hs > limits$response[2]] <- NA_real_

  centers <- seq(0, 330, by = sector_width)
  starts <- (centers - sector_width / 2) %% 360
  ends <- (centers + sector_width / 2) %% 360

  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x)) fun(x) else NA_real_
  }

  sectors <- vector("list", length(centers))
  for (i in seq_along(centers)) {
    in_grid <- diform3d_angle_in_sector(a_grid, starts[i], ends[i])
    vals <- hs[, in_grid, drop = FALSE]
    sectors[[i]] <- list(
      center = centers[i],
      start = starts[i],
      end = ends[i],
      label = diform3d_sector_label(starts[i], ends[i]),
      hi = apply(vals, 1, finite_summary, fun = max)
    )
  }

  obs_plot <- diform3d_subsample_obs(obs, max_background_points)
  obs_response <- if (!is.null(obs_plot) && nrow(obs_plot)) obs_plot$response else numeric()
  hi_all <- unlist(lapply(sectors, `[[`, "hi"), use.names = FALSE)
  y_lim <- diform3d_finite_range(0, hi_all, obs_response, limits$response)
  y_lim[1] <- 0
  y_lim[2] <- max(y_lim[2] * 1.04, 1)
  x_lim <- limits$cspd
  unit_text <- if (nzchar(response_unit)) paste0(" [", response_unit, "]") else ""

  diform3d_open_png(path, width = 1850, height = 1350, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(
    mfrow = c(3, 4),
    mar = c(3.0, 3.4, 2.4, 0.9),
    oma = c(3.8, 3.2, 4.0, 0.8),
    xpd = FALSE
  )
  on.exit(par(old), add = TRUE)

  draw_mask <- function(mask, in_sector, bg, cex = 0.70) {
    if (!isTRUE(show_points) || is.null(obs) || !nrow(obs) ||
        !length(mask) || !any(mask)) return(invisible(NULL))
    idx <- mask & in_sector
    if (!any(idx, na.rm = TRUE)) return(invisible(NULL))
    graphics::points(obs$Cspd[idx], obs$response[idx],
                     pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }

  for (sector in sectors) {
    graphics::plot(
      NA_real_, NA_real_,
      type = "n",
      xlim = x_lim,
      ylim = y_lim,
      xlab = "",
      ylab = "",
      main = sector$label,
      cex.main = 0.95,
      axes = FALSE
    )
    graphics::axis(1, cex.axis = 0.75)
    graphics::axis(2, cex.axis = 0.75, las = 1)
    graphics::box(col = "grey35")
    graphics::grid(col = diform3d_alpha("grey65", 0.35), lty = 3)

    in_obs <- if (!is.null(obs_plot) && nrow(obs_plot)) {
      diform3d_angle_in_sector(obs_plot$angularDifference, sector$start, sector$end)
    } else {
      logical()
    }

    if (isTRUE(show_points) && !is.null(obs_plot) && nrow(obs_plot) && any(in_obs)) {
      graphics::points(obs_plot$Cspd[in_obs], obs_plot$response[in_obs],
                       pch = 16, cex = 0.16,
                       col = diform3d_alpha(diform3d_obs_color, 0.16))
    }

    ok_line <- is.finite(sector$hi)
    if (sum(ok_line) >= 2) {
      graphics::lines(c_grid[ok_line], sector$hi[ok_line],
                      col = "#B42318", lwd = 1.8)
    }

    in_obs_full <- if (!is.null(obs) && nrow(obs)) {
      diform3d_angle_in_sector(obs$angularDifference, sector$start, sector$end)
    } else {
      logical()
    }
    draw_mask(diform3d_peak_mask(masks), in_obs_full, diform3d_orange_point, 0.86)
  }

  graphics::mtext(
    sprintf("D-IFORM %s-yr projected %s-Cspd envelopes by angular sector",
            format(rp, trim = TRUE), response_short),
    outer = TRUE,
    side = 3,
    line = 2.0,
    font = 2,
    cex = 1.05
  )
  graphics::mtext("Cspd [m.s^-1]", outer = TRUE, side = 1, line = 1.3, cex = 0.90)
  graphics::mtext(response_label, outer = TRUE, side = 2, line = 1.9, cex = 0.90)
  graphics::mtext(
    paste0(
      "Red line: sectoral max ", response_short, unit_text,
      " from the 3D contour."
    ),
    outer = TRUE,
    side = 1,
    line = 2.5,
    cex = 0.72,
    col = "grey25"
  )
  invisible(path)
}

diform3d_face_polygons <- function(vertices, constraints, tol = 1e-6) {
  signed <- diform3d_signed(vertices, constraints)
  scale <- pmax(1, abs(constraints$rv))
  polys <- list()
  k <- 0L

  for (i in seq_len(ncol(signed))) {
    idx <- which(abs(signed[, i]) <= tol * scale[i])
    if (length(idx) < 3) next
    pts <- vertices[idx, , drop = FALSE]
    n <- constraints$normal[i, ]
    basis <- diform3d_plane_basis(n)
    center <- colMeans(pts)
    centered <- pts - matrix(center, nrow(pts), 3, byrow = TRUE)
    uv <- cbind(
      as.numeric(centered %*% basis$e1),
      as.numeric(centered %*% basis$e2)
    )
    hull <- tryCatch(grDevices::chull(uv[, 1], uv[, 2]), error = function(e) integer())
    if (length(hull) < 3) next
    k <- k + 1L
    polys[[k]] <- pts[hull, , drop = FALSE]
  }

  polys
}

diform3d_working_polygons <- function(vertices, constraints) {
  if (!is.null(constraints$dirs) && nrow(vertices) == nrow(constraints$dirs)) {
    triangles <- diform3d_direction_triangles(constraints$dirs)
    if (!is.null(triangles) && nrow(triangles)) {
      return(lapply(seq_len(nrow(triangles)), function(i) {
        vertices[triangles[i, ], , drop = FALSE]
      }))
    }
  }
  diform3d_face_polygons(vertices, constraints)
}

diform3d_draw_polygons <- function(pm, polygons, palette, vmin, vmax, alpha = 0.58) {
  if (!length(polygons)) return(invisible(NULL))
  depth <- vapply(polygons, function(p) mean(p[, 1] + p[, 2] + p[, 3]), numeric(1))
  polygons <- polygons[order(depth, decreasing = TRUE)]
  for (p in polygons) {
    hs <- p[, 2]^2 + p[, 3]^2
    col <- diform3d_color_values(mean(hs), vmin, vmax, palette, alpha)
    xy <- grDevices::trans3d(p[, 1], p[, 2], p[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = col, border = NA)
  }
}

diform3d_halfspace_working_box <- function(support_info, limits) {
  xlim <- limits$cspd
  ylim <- limits$working_y
  zlim <- limits$working_z
  verts <- as.matrix(expand.grid(x = xlim, y = ylim, z = zlim))
  edges <- combn(seq_len(8), 2)
  edge_diff <- abs(verts[edges[1, ], , drop = FALSE] - verts[edges[2, ], , drop = FALSE])
  edges <- edges[, rowSums(edge_diff > 1e-12) == 1,
                 drop = FALSE]

  s <- as.numeric(verts %*% support_info$normal + support_info$offset)
  pts <- verts[s >= -1e-10, , drop = FALSE]

  for (e in seq_len(ncol(edges))) {
    i <- edges[1, e]
    j <- edges[2, e]
    if ((s[i] < -1e-10 && s[j] > 1e-10) ||
        (s[i] > 1e-10 && s[j] < -1e-10)) {
      t <- s[i] / (s[i] - s[j])
      pts <- rbind(pts, verts[i, ] + t * (verts[j, ] - verts[i, ]))
    }
  }

  if (nrow(pts) < 4) return(list())
  pts <- unique(round(pts, 10))
  n <- support_info$normal

  face_from_axis <- function(axis, value) {
    idx <- which(abs(pts[, axis] - value) < 1e-8)
    if (length(idx) < 3) return(NULL)
    other <- setdiff(seq_len(3), axis)
    h <- grDevices::chull(pts[idx, other[1]], pts[idx, other[2]])
    pts[idx[h], , drop = FALSE]
  }

  plane_face <- function() {
    idx <- which(abs(as.numeric(pts %*% n + support_info$offset)) < 1e-8)
    if (length(idx) < 3) return(NULL)
    basis <- diform3d_plane_basis(n)
    center <- colMeans(pts[idx, , drop = FALSE])
    centered <- pts[idx, , drop = FALSE] - matrix(center, length(idx), 3, byrow = TRUE)
    uv <- cbind(
      as.numeric(centered %*% basis$e1),
      as.numeric(centered %*% basis$e2)
    )
    h <- grDevices::chull(uv[, 1], uv[, 2])
    pts[idx[h], , drop = FALSE]
  }

  faces <- list(
    face_from_axis(1, xlim[1]), face_from_axis(1, xlim[2]),
    face_from_axis(2, ylim[1]), face_from_axis(2, ylim[2]),
    face_from_axis(3, zlim[1]), face_from_axis(3, zlim[2]),
    plane_face()
  )
  Filter(Negate(is.null), faces)
}

diform3d_plot_working <- function(path, vertices, constraints, support_info,
                                  all_vertices, obs, masks,
                                  response_label, response_short, response_unit,
                                  n_years, sample_label, rp, elev, azim,
                                  show_points, show_support_plane,
                                  show_peak_highlights, plane_scale,
                                  max_background_points, max_labels,
                                  working_z_max = NULL,
                                  meta = NULL) {
  if (is.null(meta)) {
    meta <- list(linear = "Cspd", angle = "angularDifference",
                 h_names = c("H_L", "H_T"), diform_spec = NULL)
  }
  working_labels <- diform3d_working_labels(meta)
  palette <- diform3d_palette(256)
  limits <- diform3d_limits(all_vertices, obs, working_z_max = working_z_max)
  vmax <- max(all_vertices[, "H_L"]^2 + all_vertices[, "H_T"]^2, na.rm = TRUE) * 1.03

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(8.3, 1.3))
  old <- par(mar = c(4.8, 4.9, 4.8, 1.0), xpd = NA)
  on.exit(par(old), add = TRUE)

  empty_z <- matrix(limits$working_z[1], 2, 2)
  pm <- graphics::persp(
    x = limits$cspd,
    y = limits$working_y,
    z = empty_z,
    zlim = limits$working_z,
    col = NA,
    border = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = working_labels$x,
    ylab = working_labels$y,
    zlab = working_labels$z,
    main = sprintf("D-IFORM %s-yr contour in (%s) space\n%s",
                   format(rp, trim = TRUE),
                   working_labels$title,
                   sample_label),
    cex.lab = 1.28,
    cex.axis = 1.04,
    cex.main = 1.34,
    font.main = 2
  )

  polygons <- diform3d_working_polygons(vertices, constraints)
  diform3d_draw_polygons(pm, polygons, palette, 0, vmax, alpha = 0.58)

  if (isTRUE(show_support_plane)) {
    lim <- rbind(
      c(limits$cspd[1], limits$working_y[1], limits$working_z[1]),
      c(limits$cspd[2], limits$working_y[2], limits$working_z[2])
    )
    colnames(lim) <- c("Cspd", "H_L", "H_T")
    parts <- diform3d_halfspace_box_parts_generic(support_info, lim)
    if (!is.null(parts$plane)) {
      diform3d_draw_face_list(
        pm,
        list(parts$plane),
        col = diform3d_alpha(diform3d_tangent_color, 0.32),
        border = diform3d_alpha(diform3d_tangent_color, 0.88)
      )
    }
    p <- grDevices::trans3d(support_info$x_star[1], support_info$x_star[2],
                            support_info$x_star[3], pm)
    graphics::points(p, pch = 18, col = diform3d_tangent_color, cex = 1.08)
  }

  obs_bg <- diform3d_subsample_obs(obs, max_background_points)
  if (isTRUE(show_points)) {
    diform3d_draw_points_working(pm, obs, masks, limits, max_labels,
                                 show_peak_highlights = show_peak_highlights,
                                 obs_background = obs_bg)
  }

  if (isTRUE(show_points) && isTRUE(show_peak_highlights) && !is.null(obs) && nrow(obs)) {
    diform3d_add_legend(show_support_plane)
  }
  diform3d_draw_colorbar(0, vmax, palette, response_label)
  invisible(path)
}

diform3d_response_component_index <- function(diform_result, vertices) {
  comp <- NULL
  if (!is.null(diform_result$diform_spec$response_component)) {
    comp <- diform_result$diform_spec$response_component
  }
  if (is.null(comp)) comp <- diform_result$response
  idx <- match(comp, colnames(vertices))
  if (is.na(idx)) idx <- 3L
  idx
}

diform3d_obs_generic <- function(diform_result, obs_data = NULL, excess_df = NULL) {
  if (is.null(obs_data) || !nrow(obs_data)) return(NULL)
  x <- tryCatch(
    transform_to_diform_space(obs_data, diform_result$response, diform_result$diform_spec),
    error = function(e) NULL
  )
  if (is.null(x)) return(NULL)

  nm <- diform_result$component_names
  if (is.null(nm) || length(nm) != 3) nm <- colnames(x)
  colnames(x) <- nm

  obs <- as.data.frame(x)
  obs$is_peak <- FALSE
  if (!is.null(excess_df) && nrow(excess_df) &&
      "time" %in% names(obs_data) && "time" %in% names(excess_df)) {
    obs$is_peak <- obs_data$time %in% excess_df$time
  }
  obs <- obs[complete.cases(obs[, nm, drop = FALSE]), , drop = FALSE]
  rownames(obs) <- NULL
  obs
}

diform3d_limits_generic <- function(all_vertices, obs = NULL, response_col = 3L,
                                    response_axis_max = NULL) {
  lim <- matrix(NA_real_, nrow = 2, ncol = 3)
  colnames(lim) <- colnames(all_vertices)
  for (j in seq_len(3)) {
    vals <- all_vertices[, j]
    if (!is.null(obs) && nrow(obs)) vals <- c(vals, obs[[j]])
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      r <- c(-1, 1)
    } else {
      r <- range(vals)
      span <- diff(r)
      if (!is.finite(span) || span <= 0) span <- max(1, abs(r[1]) * 0.1)
      r <- r + c(-1, 1) * span * 0.08
    }
    lim[, j] <- r
  }
  if (!is.null(response_axis_max)) {
    mx <- diform3d_axis_max(response_axis_max, lim[2, response_col])
    lim[2, response_col] <- max(mx, lim[1, response_col] + 1e-8)
  }
  lim
}

diform3d_support_info_generic <- function(vertices, constraints, response_col = 3L) {
  values <- vertices[, response_col]
  i_star <- which.max(values)
  x_star <- vertices[i_star, ]
  interior <- colMeans(vertices)
  grad <- rep(0, 3)
  grad[response_col] <- 1

  triangles <- NULL
  if (!is.null(constraints$dirs) && nrow(vertices) == nrow(constraints$dirs)) {
    triangles <- diform3d_direction_triangles(constraints$dirs)
  }

  if (!is.null(triangles) && nrow(triangles)) {
    facets <- triangles[apply(triangles == i_star, 1, any), , drop = FALSE]
    best_score <- -Inf
    best_normal <- NULL
    best_offset <- NULL

    for (j in seq_len(nrow(facets))) {
      pts <- vertices[facets[j, ], , drop = FALSE]
      normal <- diform3d_cross(pts[2, ] - pts[1, ], pts[3, ] - pts[1, ])
      if (sqrt(sum(normal^2)) <= 1e-12) next
      offset <- -sum(normal * pts[1, ])
      if (sum(normal * interior) + offset > 0) {
        normal <- -normal
        offset <- -offset
      }
      score <- sum(diform3d_unit(normal) * grad)
      if (is.finite(score) && score > best_score) {
        best_score <- score
        best_normal <- normal
        best_offset <- offset
      }
    }

    if (!is.null(best_normal)) {
      return(list(
        x_star = as.numeric(x_star),
        target_value = values[i_star],
        response_col = response_col,
        normal = diform3d_unit(best_normal),
        offset = best_offset / sqrt(sum(best_normal^2)),
        source = "surface_facet"
      ))
    }
  }

  signed <- as.numeric(x_star %*% t(constraints$normal) + constraints$offset)
  scale <- pmax(1, abs(constraints$rv))
  active <- which(abs(signed) <= 1e-6 * scale)
  if (!length(active)) active <- which.max(signed)
  score <- apply(constraints$normal[active, , drop = FALSE], 1, function(z) {
    sum(diform3d_unit(z) * grad)
  })
  i <- active[which.max(score)]
  list(
    x_star = as.numeric(x_star),
    target_value = values[i_star],
    response_col = response_col,
    normal = as.numeric(constraints$normal[i, ]),
    offset = constraints$offset[i],
    source = "active_constraint"
  )
}

diform3d_halfspace_box_generic <- function(support_info, limits) {
  verts <- as.matrix(expand.grid(x = limits[, 1], y = limits[, 2], z = limits[, 3]))
  edges <- combn(seq_len(8), 2)
  edge_diff <- abs(verts[edges[1, ], , drop = FALSE] - verts[edges[2, ], , drop = FALSE])
  edges <- edges[, rowSums(edge_diff > 1e-12) == 1,
                 drop = FALSE]

  s <- as.numeric(verts %*% support_info$normal + support_info$offset)
  pts <- verts[s >= -1e-10, , drop = FALSE]

  for (e in seq_len(ncol(edges))) {
    i <- edges[1, e]
    j <- edges[2, e]
    if ((s[i] < -1e-10 && s[j] > 1e-10) ||
        (s[i] > 1e-10 && s[j] < -1e-10)) {
      t <- s[i] / (s[i] - s[j])
      pts <- rbind(pts, verts[i, ] + t * (verts[j, ] - verts[i, ]))
    }
  }

  if (nrow(pts) < 4) return(list())
  pts <- unique(round(pts, 10))
  n <- support_info$normal

  face_from_axis <- function(axis, value) {
    idx <- which(abs(pts[, axis] - value) < 1e-8)
    if (length(idx) < 3) return(NULL)
    other <- setdiff(seq_len(3), axis)
    h <- grDevices::chull(pts[idx, other[1]], pts[idx, other[2]])
    pts[idx[h], , drop = FALSE]
  }

  plane_face <- function() {
    idx <- which(abs(as.numeric(pts %*% n + support_info$offset)) < 1e-8)
    if (length(idx) < 3) return(NULL)
    basis <- diform3d_plane_basis(n)
    center <- colMeans(pts[idx, , drop = FALSE])
    centered <- pts[idx, , drop = FALSE] - matrix(center, length(idx), 3, byrow = TRUE)
    uv <- cbind(
      as.numeric(centered %*% basis$e1),
      as.numeric(centered %*% basis$e2)
    )
    h <- grDevices::chull(uv[, 1], uv[, 2])
    pts[idx[h], , drop = FALSE]
  }

  faces <- list(
    face_from_axis(1, limits[1, 1]), face_from_axis(1, limits[2, 1]),
    face_from_axis(2, limits[1, 2]), face_from_axis(2, limits[2, 2]),
    face_from_axis(3, limits[1, 3]), face_from_axis(3, limits[2, 3]),
    plane_face()
  )
  Filter(Negate(is.null), faces)
}

diform3d_halfspace_box_parts_generic <- function(support_info, limits) {
  verts <- as.matrix(expand.grid(x = limits[, 1], y = limits[, 2], z = limits[, 3]))
  edges <- combn(seq_len(8), 2)
  edge_diff <- abs(verts[edges[1, ], , drop = FALSE] - verts[edges[2, ], , drop = FALSE])
  edges <- edges[, rowSums(edge_diff > 1e-12) == 1,
                 drop = FALSE]

  s <- as.numeric(verts %*% support_info$normal + support_info$offset)
  pts <- verts[s >= -1e-10, , drop = FALSE]

  for (e in seq_len(ncol(edges))) {
    i <- edges[1, e]
    j <- edges[2, e]
    if ((s[i] < -1e-10 && s[j] > 1e-10) ||
        (s[i] > 1e-10 && s[j] < -1e-10)) {
      t <- s[i] / (s[i] - s[j])
      pts <- rbind(pts, verts[i, ] + t * (verts[j, ] - verts[i, ]))
    }
  }

  if (nrow(pts) < 4) return(list(boundary = list(), plane = NULL))
  pts <- unique(round(pts, 10))
  n <- support_info$normal

  face_from_axis <- function(axis, value) {
    idx <- which(abs(pts[, axis] - value) < 1e-8)
    if (length(idx) < 3) return(NULL)
    other <- setdiff(seq_len(3), axis)
    h <- grDevices::chull(pts[idx, other[1]], pts[idx, other[2]])
    pts[idx[h], , drop = FALSE]
  }

  plane_face <- function() {
    idx <- which(abs(as.numeric(pts %*% n + support_info$offset)) < 1e-8)
    if (length(idx) < 3) return(NULL)
    basis <- diform3d_plane_basis(n)
    center <- colMeans(pts[idx, , drop = FALSE])
    centered <- pts[idx, , drop = FALSE] - matrix(center, length(idx), 3, byrow = TRUE)
    uv <- cbind(
      as.numeric(centered %*% basis$e1),
      as.numeric(centered %*% basis$e2)
    )
    h <- grDevices::chull(uv[, 1], uv[, 2])
    pts[idx[h], , drop = FALSE]
  }

  boundary <- list(
    face_from_axis(1, limits[1, 1]), face_from_axis(1, limits[2, 1]),
    face_from_axis(2, limits[1, 2]), face_from_axis(2, limits[2, 2]),
    face_from_axis(3, limits[1, 3]), face_from_axis(3, limits[2, 3])
  )

  list(
    boundary = Filter(Negate(is.null), boundary),
    plane = plane_face()
  )
}

diform3d_halfspace_local_layers_generic <- function(support_info, vertices,
                                                    scale = 0.34,
                                                    n = 16,
                                                    n_layers = 6) {
  ranges <- apply(vertices, 2, range)
  radius <- scale * max(ranges[2, ] - ranges[1, ])
  basis <- diform3d_plane_basis(support_info$normal)
  normal <- diform3d_unit(support_info$normal)
  u <- seq(-radius, radius, length.out = n)
  v <- seq(-radius, radius, length.out = n)
  grid <- expand.grid(u = u, v = v)
  offsets <- seq(0.06, 0.70, length.out = n_layers) * radius

  lapply(offsets, function(d) {
    center <- support_info$x_star + d * normal
    pts <- t(apply(grid, 1, function(g) {
      center + g[1] * basis$e1 + g[2] * basis$e2
    }))
    list(
      X = matrix(pts[, 1], n, n),
      Y = matrix(pts[, 2], n, n),
      Z = matrix(pts[, 3], n, n)
    )
  })
}

diform3d_draw_polygons_generic <- function(pm, polygons, response_col,
                                           palette, vmin, vmax, alpha = 0.58) {
  if (!length(polygons)) return(invisible(NULL))
  depth <- vapply(polygons, function(p) mean(p[, 1] + p[, 2] + p[, 3]), numeric(1))
  polygons <- polygons[order(depth, decreasing = TRUE)]
  for (p in polygons) {
    col <- diform3d_color_values(mean(p[, response_col]), vmin, vmax, palette, alpha)
    xy <- grDevices::trans3d(p[, 1], p[, 2], p[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = col, border = NA)
  }
}

diform3d_highlight_masks_generic <- function(obs, constraints, support_info, tol = 1e-10) {
  if (is.null(obs) || !nrow(obs)) return(list(orange = logical(), green = logical()))
  pts <- as.matrix(obs[, seq_len(3), drop = FALSE])
  outside <- apply(diform3d_signed(pts, constraints) > tol, 1, any)
  halfspace <- as.numeric(pts %*% support_info$normal + support_info$offset) > tol
  peaks <- isTRUE(any(obs$is_peak)) & obs$is_peak
  list(
    orange = peaks & outside & !halfspace,
    green = peaks & outside & halfspace
  )
}

diform3d_draw_points_generic <- function(pm, obs, masks, limits, response_col,
                                         max_labels = 6,
                                         show_peak_highlights = TRUE,
                                         obs_background = obs) {
  if (is.null(obs) || !nrow(obs)) return(invisible(NULL))
  pts <- as.matrix(obs[, seq_len(3), drop = FALSE])
  if (!is.null(obs_background) && nrow(obs_background)) {
    bg_pts <- as.matrix(obs_background[, seq_len(3), drop = FALSE])
    xyz <- grDevices::trans3d(bg_pts[, 1], bg_pts[, 2], bg_pts[, 3], pm)
    graphics::points(xyz, pch = 16, cex = 0.17,
                     col = diform3d_alpha(diform3d_obs_color, 0.040))
  }
  if (!isTRUE(show_peak_highlights)) return(invisible(NULL))

  draw_mask <- function(mask, bg, cex) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    p <- grDevices::trans3d(pts[mask, 1], pts[mask, 2], pts[mask, 3], pm)
    graphics::points(p, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }
  peak_mask <- diform3d_peak_mask(masks)
  draw_mask(peak_mask, diform3d_orange_point, 0.90)

  label_mask <- function(mask, col) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    idx <- which(mask)
    idx <- idx[order(pts[idx, response_col], decreasing = TRUE)]
    idx <- idx[seq_len(min(max_labels, length(idx)))]
    p <- grDevices::trans3d(
      pts[idx, 1], pts[idx, 2],
      pts[idx, 3] + 0.018 * diff(limits[, 3]), pm
    )
    graphics::text(p$x, p$y, sprintf("%.2f", pts[idx, response_col]),
                   cex = 0.70, col = col, pos = 3)
  }
  label_mask(peak_mask, "darkorange")
  invisible(NULL)
}

diform3d_add_legend_generic <- function(response_short, show_support_plane = TRUE,
                                        show_exterior_halfspace = FALSE) {
  legend <- c("Original data", "Peaks outside contour")
  pch <- c(16, 21)
  pt_bg <- c(NA, diform3d_orange_point)
  cols <- c(diform3d_alpha(diform3d_obs_color, 0.62), "black")

  if (isTRUE(show_support_plane) && isTRUE(show_exterior_halfspace)) {
    legend <- c(legend, "Exterior half-space (local)")
    pch <- c(pch, 15)
    pt_bg <- c(pt_bg, diform3d_alpha(diform3d_halfspace_color, 0.48))
    cols <- c(cols, diform3d_alpha(diform3d_halfspace_color, 0.85))
  }

  if (isTRUE(show_support_plane)) {
    legend <- c(legend, paste0("Tangent plane at max ", response_short))
    pch <- c(pch, 15)
    pt_bg <- c(pt_bg, diform3d_alpha(diform3d_tangent_color, 0.34))
    cols <- c(cols, diform3d_alpha(diform3d_tangent_color, 0.88))
  }

  graphics::legend(
    "bottomright",
    inset = c(0.015, 0.020),
    legend = legend,
    pch = pch,
    pt.bg = pt_bg,
    col = cols,
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.94,
    y.intersp = 1.14
  )
}

diform3d_plot_generic <- function(path, vertices, constraints, support_info,
                                  all_vertices, obs, masks,
                                  response_label, response_short, response_unit,
                                  n_years, sample_label, rp, elev, azim,
                                  show_points, show_support_plane,
                                  show_peak_highlights, plane_scale,
                                  max_background_points, max_labels,
                                  z_max = NULL,
                                  space_label = "selected coordinate space",
                                  show_exterior_halfspace = FALSE,
                                  diform_spec = NULL) {
  palette <- diform3d_palette(256)
  response_col <- support_info$response_col
  limits <- diform3d_limits_generic(all_vertices, obs, response_col, z_max)
  vals <- all_vertices[, response_col]
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) {
    vmin <- 0
    vmax <- 1
  }
  comp_names <- colnames(vertices)
  comp_labels <- diform3d_component_labels(
    comp_names,
    diform_spec = diform_spec,
    response_label = response_label,
    response_col = response_col
  )

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2), 1, 2), widths = c(8.3, 1.3))
  old <- par(mar = c(5.2, 5.4, 5.0, 1.1), xpd = NA)
  on.exit(par(old), add = TRUE)

  empty_z <- matrix(limits[1, 3], 2, 2)
  pm <- graphics::persp(
    x = limits[, 1],
    y = limits[, 2],
    z = empty_z,
    zlim = limits[, 3],
    col = NA,
    border = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = comp_labels[1],
    ylab = comp_labels[2],
    zlab = comp_labels[3],
    main = sprintf("D-IFORM %s-yr contour\n%s",
                   format(rp, trim = TRUE),
                   sample_label),
    cex.lab = 1.28,
    cex.axis = 1.04,
    cex.main = 1.34,
    font.main = 2
  )

  if (isTRUE(show_support_plane) && isTRUE(show_exterior_halfspace)) {
    layers <- diform3d_halfspace_local_layers_generic(
      support_info, vertices,
      scale = max(plane_scale * 1.15, 0.28),
      n = 16,
      n_layers = 6
    )
    for (layer in layers) {
      diform3d_project_grid(pm, layer$X, layer$Y, layer$Z,
                            col = diform3d_alpha(diform3d_halfspace_color, 0.075),
                            border = NA, zlim = limits[, 3])
    }
  }

  polygons <- diform3d_working_polygons(vertices, constraints)
  diform3d_draw_polygons_generic(pm, polygons, response_col, palette, vmin, vmax, alpha = 0.58)

  if (isTRUE(show_support_plane)) {
    parts <- diform3d_halfspace_box_parts_generic(support_info, limits)
    if (!is.null(parts$plane)) {
      diform3d_draw_face_list(
        pm,
        list(parts$plane),
        col = diform3d_alpha(diform3d_tangent_color, 0.32),
        border = diform3d_alpha(diform3d_tangent_color, 0.88)
      )
    }
    p <- grDevices::trans3d(support_info$x_star[1], support_info$x_star[2],
                            support_info$x_star[3], pm)
    graphics::points(p, pch = 18, col = diform3d_tangent_color, cex = 1.08)
  }

  obs_bg <- diform3d_subsample_obs(obs, max_background_points)
  if (isTRUE(show_points)) {
    diform3d_draw_points_generic(pm, obs, masks, limits, response_col, max_labels,
                                 show_peak_highlights = show_peak_highlights,
                                 obs_background = obs_bg)
  }

  if (isTRUE(show_points) && isTRUE(show_peak_highlights) && !is.null(obs) && nrow(obs)) {
    diform3d_add_legend_generic(response_short, show_support_plane,
                                show_exterior_halfspace)
  }
  diform3d_draw_colorbar(vmin, vmax, palette, response_label)
  invisible(path)
}

diform3d_draw_face_list <- function(pm, faces, col, border) {
  if (!length(faces)) return(invisible(NULL))
  depth <- vapply(faces, function(p) mean(p[, 1] + p[, 2] + p[, 3]), numeric(1))
  faces <- faces[order(depth, decreasing = TRUE)]
  for (p in faces) {
    xy <- grDevices::trans3d(p[, 1], p[, 2], p[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = col, border = border, lwd = 0.8)
  }
  invisible(NULL)
}

diform3d_plot_halfspace_generic <- function(path, vertices, support_info,
                                            all_vertices, obs,
                                            response_label, response_short,
                                            response_unit, n_years, rp,
                                            elev, azim, z_max = NULL,
                                            space_label = "selected coordinate space") {
  response_col <- support_info$response_col
  limits <- diform3d_limits_generic(all_vertices, obs, response_col, z_max)
  comp_names <- colnames(vertices)
  parts <- diform3d_halfspace_box_parts_generic(support_info, limits)

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mar = c(4.8, 4.9, 5.1, 1.2), xpd = NA)
  on.exit(par(old), add = TRUE)

  empty_z <- matrix(limits[1, 3], 2, 2)
  pm <- graphics::persp(
    x = limits[, 1],
    y = limits[, 2],
    z = empty_z,
    zlim = limits[, 3],
    col = NA,
    border = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = comp_names[1],
    ylab = comp_names[2],
    zlab = comp_names[3],
    main = sprintf("Exterior half-space for D-IFORM %s-yr contour\n%s",
                   format(rp, trim = TRUE), space_label),
    cex.lab = 0.98,
    cex.axis = 0.82,
    cex.main = 1.14,
    font.main = 2
  )

  diform3d_draw_face_list(
    pm,
    parts$boundary,
    col = diform3d_alpha("#4C78A8", 0.22),
    border = diform3d_alpha("#2F4B7C", 0.55)
  )

  if (!is.null(parts$plane)) {
    diform3d_draw_face_list(
      pm,
      list(parts$plane),
      col = diform3d_alpha(diform3d_tangent_color, 0.42),
      border = diform3d_alpha(diform3d_tangent_color, 0.88)
    )
  }

  star <- grDevices::trans3d(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3],
    pm
  )
  graphics::points(star, pch = 18, col = diform3d_tangent_color, cex = 1.25)

  graphics::legend(
    "topright",
    inset = 0.01,
    legend = c("Exterior half-space", "Tangent-plane boundary", paste0("Max ", response_short, " point")),
    pch = c(15, 15, 18),
    pt.bg = c(diform3d_alpha("#4C78A8", 0.38),
              diform3d_alpha(diform3d_tangent_color, 0.46),
              NA),
    col = c(diform3d_alpha("#2F4B7C", 0.70),
            diform3d_alpha(diform3d_tangent_color, 0.90),
            diform3d_tangent_color),
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.76,
    y.intersp = 1.12
  )

  invisible(path)
}

diform3d_plot_halfspace_working <- function(path, vertices, support_info,
                                            all_vertices, obs,
                                            response_short, response_unit,
                                            n_years, rp, elev, azim,
                                            working_z_max = NULL,
                                            meta = NULL) {
  if (is.null(meta)) {
    meta <- list(linear = "Cspd", angle = "angularDifference",
                 h_names = c("H_L", "H_T"), diform_spec = NULL)
  }
  working_labels <- diform3d_working_labels(meta)
  limits <- diform3d_limits(all_vertices, obs, working_z_max = working_z_max)
  lim <- rbind(
    c(limits$cspd[1], limits$working_y[1], limits$working_z[1]),
    c(limits$cspd[2], limits$working_y[2], limits$working_z[2])
  )
  colnames(lim) <- c("Cspd", "H_L", "H_T")
  parts <- diform3d_halfspace_box_parts_generic(support_info, lim)

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mar = c(4.8, 4.9, 5.1, 1.2), xpd = NA)
  on.exit(par(old), add = TRUE)

  empty_z <- matrix(lim[1, 3], 2, 2)
  pm <- graphics::persp(
    x = lim[, 1],
    y = lim[, 2],
    z = empty_z,
    zlim = lim[, 3],
    col = NA,
    border = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = working_labels$x,
    ylab = working_labels$y,
    zlab = working_labels$z,
    main = sprintf("Exterior half-space for D-IFORM %s-yr contour\n(%s) space",
                   format(rp, trim = TRUE),
                   working_labels$title),
    cex.lab = 0.95,
    cex.axis = 0.80,
    cex.main = 1.14,
    font.main = 2
  )

  diform3d_draw_face_list(
    pm,
    parts$boundary,
    col = diform3d_alpha("#4C78A8", 0.22),
    border = diform3d_alpha("#2F4B7C", 0.55)
  )

  if (!is.null(parts$plane)) {
    diform3d_draw_face_list(
      pm,
      list(parts$plane),
      col = diform3d_alpha(diform3d_tangent_color, 0.42),
      border = diform3d_alpha(diform3d_tangent_color, 0.88)
    )
  }

  star <- grDevices::trans3d(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3],
    pm
  )
  graphics::points(star, pch = 18, col = diform3d_tangent_color, cex = 1.25)

  graphics::legend(
    "topright",
    inset = 0.01,
    legend = c("Exterior half-space", "Tangent-plane boundary", paste0("Max ", response_short, " point")),
    pch = c(15, 15, 18),
    pt.bg = c(diform3d_alpha("#4C78A8", 0.38),
              diform3d_alpha(diform3d_tangent_color, 0.46),
              NA),
    col = c(diform3d_alpha("#2F4B7C", 0.70),
            diform3d_alpha(diform3d_tangent_color, 0.90),
            diform3d_tangent_color),
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.76,
    y.intersp = 1.12
  )

  invisible(path)
}

diform3d_plot_halfspace_physical <- function(path, support_info,
                                             all_vertices, obs,
                                             response_label, response_short,
                                             response_unit, n_years, rp,
                                             elev, azim, n_ang, n_cspd,
                                             z_max = NULL) {
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  layers <- diform3d_halfspace_physical_layers(
    support_info,
    limits,
    n_cspd = n_cspd,
    n_ang = n_ang,
    n_layers = 9
  )
  boundary <- diform3d_halfspace_physical_boundary(
    support_info,
    limits,
    n_cspd = n_cspd,
    n_ang = n_ang
  )

  diform3d_open_png(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mar = c(4.8, 4.9, 5.1, 1.2), xpd = NA)
  on.exit(par(old), add = TRUE)

  empty_z <- matrix(limits$response[1], 2, 2)
  pm <- graphics::persp(
    x = limits$cspd,
    y = limits$angle,
    z = empty_z,
    zlim = limits$response,
    col = NA,
    border = NA,
    ticktype = "detailed",
    theta = azim,
    phi = elev,
    r = 7,
    expand = 0.88,
    xlab = "Current speed [m.s^-1]",
    ylab = "Angular difference [deg]",
    zlab = response_label,
    main = sprintf("Exterior half-space mapped to physical space\nD-IFORM %s-yr contour",
                   format(rp, trim = TRUE)),
    cex.lab = 1.02,
    cex.axis = 0.82,
    cex.main = 1.14,
    font.main = 2
  )

  max_cell_span <- diff(limits$response) * 0.25
  for (layer in layers) {
    diform3d_project_grid(
      pm,
      layer$X, layer$Y, layer$Z,
      col = diform3d_alpha("#4C78A8", 0.055),
      border = NA,
      zlim = limits$response,
      max_z_span = max_cell_span
    )
  }

  diform3d_project_grid(
    pm,
    boundary$X, boundary$Y, boundary$Z,
    col = diform3d_alpha(diform3d_tangent_color, 0.55),
    border = diform3d_alpha(diform3d_tangent_color, 0.25),
    zlim = limits$response,
    max_z_span = max_cell_span
  )

  star <- diform3d_physical(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3]
  )
  p <- grDevices::trans3d(star$Cspd, star$angularDifference, star$response, pm)
  graphics::points(p, pch = 18, col = diform3d_tangent_color, cex = 1.25)

  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""
  diform3d_add_annotation(c(
    "Blue sheets: samples inside n.x + b > 0",
    "Purple surface: boundary n.x + b = 0",
    sprintf("%s = %.3f%s", response_short, star$response, unit),
    diform3d_box_label("angularDifference", star$angularDifference)
  ))

  graphics::legend(
    "topright",
    inset = 0.01,
    legend = c("Mapped exterior half-space", "Mapped tangent boundary", paste0("Max ", response_short, " point")),
    pch = c(15, 15, 18),
    pt.bg = c(diform3d_alpha("#4C78A8", 0.38),
              diform3d_alpha(diform3d_tangent_color, 0.46),
              NA),
    col = c(diform3d_alpha("#2F4B7C", 0.70),
            diform3d_alpha(diform3d_tangent_color, 0.90),
            diform3d_tangent_color),
    bty = "o",
    bg = diform3d_alpha("white", 0.92),
    box.col = "lightgray",
    cex = 0.76,
    y.intersp = 1.12
  )

  invisible(path)
}

diform3d_plot_halfspace_physical_slice <- function(path, support_info,
                                                   all_vertices, obs,
                                                   response_label, response_short,
                                                   response_unit, n_years, rp,
                                                   z_max = NULL,
                                                   n_ang = 720,
                                                   n_h = 320) {
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  c0 <- support_info$x_star[1]
  ag <- seq(limits$angle[1], limits$angle[2], length.out = n_ang)
  hg <- seq(limits$response[1], limits$response[2], length.out = n_h)
  n <- support_info$normal

  signed <- outer(ag, hg, function(a, h) {
    r <- sqrt(pmax(h, 0))
    n[1] * c0 +
      n[2] * r * cos(a * pi / 180) +
      n[3] * r * sin(a * pi / 180) +
      support_info$offset
  })

  star <- diform3d_physical(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3]
  )
  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""

  diform3d_open_png(path, width = 1450, height = 1050, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mar = c(5.0, 5.1, 4.5, 2.0), xpd = NA)
  on.exit(par(old), add = TRUE)

  graphics::image(
    x = ag,
    y = hg,
    z = signed > 0,
    col = c(diform3d_alpha("white", 1), diform3d_alpha("#4C78A8", 0.36)),
    xlab = "Angular difference [deg]",
    ylab = response_label,
    main = sprintf("Exterior half-space slice at Cspd = %.2f m.s^-1\nD-IFORM %s-yr contour",
                   c0, format(rp, trim = TRUE)),
    useRaster = TRUE,
    axes = FALSE
  )
  graphics::contour(
    x = ag,
    y = hg,
    z = signed,
    levels = 0,
    add = TRUE,
    drawlabels = FALSE,
    col = diform3d_tangent_color,
    lwd = 2.2
  )
  graphics::points(star$angularDifference, star$response,
                   pch = 18, cex = 1.4, col = diform3d_tangent_color)
  graphics::axis(1, at = seq(0, 360, by = 45))
  graphics::axis(2, las = 1)
  graphics::box()
  graphics::legend(
    "topright",
    legend = c("n.x + b > 0", "n.x + b = 0", paste0("Max ", response_short)),
    pch = c(15, NA, 18),
    lty = c(NA, 1, NA),
    pt.bg = c(diform3d_alpha("#4C78A8", 0.36), NA, NA),
    col = c(diform3d_alpha("#2F4B7C", 0.85), diform3d_tangent_color, diform3d_tangent_color),
    lwd = c(NA, 2.2, NA),
    bty = "o",
    bg = diform3d_alpha("white", 0.94),
    box.col = "lightgray",
    cex = 0.82
  )
  graphics::mtext(sprintf("%s = %.3f%s, Angular difference = %.0f\u00B0",
                          response_short, star$response, unit, star$angularDifference),
                  side = 3, line = 0.3, cex = 0.82)

  invisible(path)
}

diform3d_draw_halfspace_slice_panel <- function(x, y, signed, xlab, ylab, main,
                                                point_x, point_y,
                                                obs_x = NULL,
                                                obs_y = NULL,
                                                contour_x = NULL,
                                                contour_y = NULL,
                                                contour_segments = NULL,
                                                x_at = NULL,
                                                y_at = NULL) {
  graphics::image(
    x = x,
    y = y,
    z = signed > 0,
    col = c(diform3d_alpha("white", 1), diform3d_alpha("#4C78A8", 0.34)),
    xlab = xlab,
    ylab = ylab,
    main = main,
    useRaster = TRUE,
    axes = FALSE
  )
  if (!is.null(obs_x) && !is.null(obs_y)) {
    ok <- is.finite(obs_x) & is.finite(obs_y)
    if (any(ok)) {
      graphics::points(obs_x[ok], obs_y[ok],
                       pch = 16, cex = 0.22,
                       col = diform3d_alpha("#111827", 0.48))
    }
  }
  if (!is.null(contour_x) && !is.null(contour_y)) {
    ok <- is.finite(contour_x) & is.finite(contour_y)
    if (sum(ok) >= 2) {
      graphics::lines(contour_x[ok], contour_y[ok],
                      col = "#111827", lwd = 1.7)
    }
  }
  if (!is.null(contour_segments) && length(contour_segments)) {
    for (seg in contour_segments) {
      ok <- is.finite(seg[, 1]) & is.finite(seg[, 2])
      if (sum(ok) >= 2) {
        graphics::lines(seg[ok, 1], seg[ok, 2],
                        col = "#111827", lwd = 1.7)
      }
    }
  }
  if (min(signed, na.rm = TRUE) <= 0 && max(signed, na.rm = TRUE) >= 0) {
    graphics::contour(
      x = x,
      y = y,
      z = signed,
      levels = 0,
      add = TRUE,
      drawlabels = FALSE,
      col = diform3d_tangent_color,
      lwd = 2.2
    )
  }
  graphics::points(point_x, point_y, pch = 18, cex = 1.35, col = diform3d_tangent_color)
  if (is.null(x_at)) graphics::axis(1) else graphics::axis(1, at = x_at)
  if (is.null(y_at)) graphics::axis(2, las = 1) else graphics::axis(2, at = y_at, las = 1)
  graphics::box()
}

diform3d_plot_halfspace_physical_slices <- function(path, support_info,
                                                    vertices,
                                                    constraints,
                                                    all_vertices, obs,
                                                    response_label, response_short,
                                                    response_unit, n_years, rp,
                                                    z_max = NULL,
                                                    n_grid = 360) {
  limits <- diform3d_limits(all_vertices, obs, z_max = z_max)
  star <- diform3d_physical(
    support_info$x_star[1],
    support_info$x_star[2],
    support_info$x_star[3]
  )
  n <- support_info$normal
  b <- support_info$offset

  cg <- seq(limits$cspd[1], limits$cspd[2], length.out = n_grid)
  ag <- seq(limits$angle[1], limits$angle[2], length.out = n_grid)
  hg <- seq(limits$response[1], limits$response[2], length.out = n_grid)

  signed_c_fixed <- outer(ag, hg, function(a, h) {
    r <- sqrt(pmax(h, 0))
    n[1] * star$Cspd +
      n[2] * r * cos(a * pi / 180) +
      n[3] * r * sin(a * pi / 180) +
      b
  })

  ca_star <- cos(star$angularDifference * pi / 180)
  sa_star <- sin(star$angularDifference * pi / 180)
  signed_a_fixed <- outer(cg, hg, function(c, h) {
    r <- sqrt(pmax(h, 0))
    n[1] * c + n[2] * r * ca_star + n[3] * r * sa_star + b
  })

  r_star <- sqrt(pmax(star$response, 0))
  signed_h_fixed <- outer(cg, ag, function(c, a) {
    n[1] * c +
      n[2] * r_star * cos(a * pi / 180) +
      n[3] * r_star * sin(a * pi / 180) +
      b
  })

  obs_plot <- diform3d_subsample_obs(obs, max_background_points = 20000)
  phys_v <- diform3d_physical(vertices[, "Cspd"], vertices[, "H_L"], vertices[, "H_T"])
  contour_ang_h <- diform3d_polar_envelope(vertices[, "H_L"], vertices[, "H_T"], n = 361)
  contour_c_h <- diform3d_hull_boundary_2d(phys_v$Cspd, phys_v$response)
  contour_c_ang <- diform3d_contour_slice_hs(vertices, constraints, star$response)

  unit <- if (nzchar(response_unit)) paste0(" ", response_unit) else ""

  diform3d_open_png(path, width = 2200, height = 850, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mfrow = c(1, 3), mar = c(4.8, 4.8, 4.2, 1.0), oma = c(0, 0, 3.2, 0))
  on.exit(par(old), add = TRUE)

  diform3d_draw_halfspace_slice_panel(
    x = ag,
    y = hg,
    signed = signed_c_fixed,
    xlab = "Angular difference [deg]",
    ylab = response_label,
    main = sprintf("Cspd = %.2f m.s^-1", star$Cspd),
    point_x = star$angularDifference,
    point_y = star$response,
    obs_x = obs_plot$angularDifference,
    obs_y = obs_plot$response,
    contour_x = contour_ang_h$angDiff,
    contour_y = contour_ang_h$resp,
    x_at = seq(0, 360, by = 90)
  )

  diform3d_draw_halfspace_slice_panel(
    x = cg,
    y = hg,
    signed = signed_a_fixed,
    xlab = "Current speed [m.s^-1]",
    ylab = response_label,
    main = sprintf("Angular difference = %.0f\u00B0", star$angularDifference),
    point_x = star$Cspd,
    point_y = star$response,
    obs_x = obs_plot$Cspd,
    obs_y = obs_plot$response,
    contour_x = contour_c_h$x,
    contour_y = contour_c_h$y
  )

  diform3d_draw_halfspace_slice_panel(
    x = cg,
    y = ag,
    signed = signed_h_fixed,
    xlab = "Current speed [m.s^-1]",
    ylab = "Angular difference [deg]",
    main = sprintf("%s = %.3f%s", response_short, star$response, unit),
    point_x = star$Cspd,
    point_y = star$angularDifference,
    obs_x = obs_plot$Cspd,
    obs_y = obs_plot$angularDifference,
    contour_segments = contour_c_ang,
    y_at = seq(0, 360, by = 90)
  )

  graphics::mtext(
    sprintf("Mapped exterior half-space slices through max %s point - D-IFORM %s-yr contour",
            response_short, format(rp, trim = TRUE)),
    outer = TRUE,
    side = 3,
    line = 1.0,
    cex = 1.18,
    font = 2
  )
  graphics::mtext(
    "Blue: n.x + b > 0; violet: tangent half-space boundary; black: projected D-IFORM contour; grey: observations",
    outer = TRUE,
    side = 3,
    line = -0.7,
    cex = 0.82
  )

  invisible(path)
}

plot_diform_3d_generic_r <- function(diform_result, obs_data = NULL, excess_df = NULL,
                                     response_label = "Response",
                                     response_short = "Response",
                                     response_unit = "",
                                     out_dir = ".",
                                     rp = "all",
                                     elev = 22,
                                     azim = -35,
                                     show_points = TRUE,
                                     show_peak_highlights = TRUE,
                                     show_support_plane = TRUE,
                                     plane_scale = 0.28,
                                     max_background_points = 20000,
                                     max_labels = 6,
                                     z_max = NULL,
                                     show_exterior_halfspace = FALSE,
                                     make_halfspace_plot = TRUE) {
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  obs <- diform3d_obs_generic(diform_result, obs_data, excess_df)
  rps <- diform3d_get_rps(diform_result, rp)
  if (!length(rps)) stop("No D-IFORM return period available for rp = ", paste(rp, collapse = ", "))

  all_rps <- diform3d_get_rps(diform_result, "all")
  all_vertices <- do.call(rbind, lapply(all_rps, function(x) diform3d_vertices(diform_result, x)))
  response_col <- diform3d_response_component_index(diform_result, all_vertices)
  all_linear <- !is.null(diform_result$diform_spec$components) &&
    all(vapply(diform_result$diform_spec$components, function(x) identical(x$type, "linear"), logical(1)))
  space_label <- if (isTRUE(all_linear)) "original coordinate space" else "selected D-IFORM coordinates"
  file_prefix <- if (isTRUE(all_linear)) "diform_3d" else "diform_3d_generic"
  sample_label <- diform3d_sample_label(diform_result)
  files <- character()

  for (rp_i in rps) {
    vertices <- diform3d_vertices(diform_result, rp_i)
    constraints <- diform3d_constraints(diform_result, rp_i)
    support_info <- diform3d_support_info_generic(vertices, constraints, response_col)
    masks <- diform3d_highlight_masks_generic(obs, constraints, support_info)

    path <- file.path(out_dir, sprintf("%s_%syr.png",
                                       file_prefix,
                                       format(rp_i, trim = TRUE)))
    diform3d_plot_generic(
      path = path,
      vertices = vertices,
      constraints = constraints,
      support_info = support_info,
      all_vertices = all_vertices,
      obs = obs,
      masks = masks,
      response_label = response_label,
      response_short = response_short,
      response_unit = response_unit,
      n_years = diform_result$n_years,
      sample_label = sample_label,
      rp = rp_i,
      elev = elev,
      azim = azim,
      show_points = show_points,
      show_support_plane = show_support_plane,
      show_peak_highlights = show_peak_highlights,
      plane_scale = plane_scale,
      max_background_points = max_background_points,
      max_labels = max_labels,
      z_max = z_max,
      space_label = space_label,
      show_exterior_halfspace = show_exterior_halfspace,
      diform_spec = diform_result$diform_spec
    )
    files <- c(files, path)
  }

  invisible(files)
}

plot_diform_3d_r <- function(diform_result, obs_data = NULL, excess_df = NULL,
                             response_label = "Response",
                             response_short = "Response",
                             response_unit = "",
                             out_dir = ".",
                             rp = "all",
                             elev = 22,
                             azim = -35,
                             n_ang = 180,
                             n_cspd = 80,
                             show_points = TRUE,
                             show_peak_highlights = TRUE,
                             show_support_plane = TRUE,
                             plane_scale = 0.28,
                             max_background_points = 20000,
                             max_labels = 6,
                             z_max = NULL,
                             working_z_max = NULL,
                             show_exterior_halfspace = FALSE,
                             make_halfspace_plot = TRUE) {
  if (!identical(diform_result$plot_type, "polar_response")) {
    return(plot_diform_3d_generic_r(
      diform_result = diform_result,
      obs_data = obs_data,
      excess_df = excess_df,
      response_label = response_label,
      response_short = response_short,
      response_unit = response_unit,
      out_dir = out_dir,
      rp = rp,
      elev = elev,
      azim = azim,
      show_points = show_points,
      show_peak_highlights = show_peak_highlights,
      show_support_plane = show_support_plane,
      plane_scale = plane_scale,
      max_background_points = max_background_points,
      max_labels = max_labels,
      z_max = z_max,
      show_exterior_halfspace = show_exterior_halfspace,
      make_halfspace_plot = make_halfspace_plot
    ))
  }

  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  response <- diform_result$response
  obs <- diform3d_obs(obs_data, excess_df, response, diform_result$diform_spec)
  meta <- diform3d_polar_meta(diform_result)
  rps <- diform3d_get_rps(diform_result, rp)
  if (!length(rps)) stop("No D-IFORM return period available for rp = ", paste(rp, collapse = ", "))

  all_rps <- diform3d_get_rps(diform_result, "all")
  all_vertices <- do.call(rbind, lapply(all_rps, function(x) diform3d_vertices(diform_result, x)))
  files <- character()
  sample_label <- diform3d_sample_label(diform_result)
  diform3d_warn_if_nonmonotone_returns(diform_result)

  for (rp_i in rps) {
    vertices <- diform3d_vertices(diform_result, rp_i)
    constraints <- diform3d_constraints(diform_result, rp_i)
    support_info <- diform3d_support_info(vertices, constraints)
    masks <- diform3d_highlight_masks(obs, constraints, support_info)

    physical_path <- file.path(out_dir, sprintf("diform_3d_%syr.png", format(rp_i, trim = TRUE)))
    diform3d_plot_physical(
      path = physical_path,
      vertices = vertices,
      constraints = constraints,
      support_info = support_info,
      all_vertices = all_vertices,
      obs = obs,
      masks = masks,
      response_label = response_label,
      response_short = response_short,
      response_unit = response_unit,
      n_years = diform_result$n_years,
      sample_label = sample_label,
      rp = rp_i,
      elev = elev,
      azim = azim,
      n_ang = n_ang,
      n_cspd = n_cspd,
      show_points = show_points,
      show_peak_highlights = show_peak_highlights,
      show_support_plane = show_support_plane,
      plane_scale = plane_scale,
      max_background_points = max_background_points,
      max_labels = max_labels,
      z_max = z_max,
      meta = meta
    )
    files <- c(files, physical_path)

    suffix <- if (isTRUE(show_support_plane)) "_with_support_plane" else "_working"
    working_path <- file.path(out_dir, sprintf("diform_3d_%syr%s.png",
                                               format(rp_i, trim = TRUE), suffix))
    diform3d_plot_working(
      path = working_path,
      vertices = vertices,
      constraints = constraints,
      support_info = support_info,
      all_vertices = all_vertices,
      obs = obs,
      masks = masks,
      response_label = response_label,
      response_short = response_short,
      response_unit = response_unit,
      n_years = diform_result$n_years,
      sample_label = sample_label,
      rp = rp_i,
      elev = elev,
      azim = azim,
      show_points = show_points,
      show_support_plane = show_support_plane,
      show_peak_highlights = show_peak_highlights,
      plane_scale = plane_scale,
      max_background_points = max_background_points,
      max_labels = max_labels,
      working_z_max = working_z_max,
      meta = meta
    )
    files <- c(files, working_path)
  }

  invisible(files)
}

diform3d_compare_layout <- function(n) {
  if (n <= 3) return(c(1, n))
  nc <- ceiling(sqrt(n))
  nr <- ceiling(n / nc)
  c(nr, nc)
}

diform3d_compare_limits <- function(vertices_list, working_z_max = NULL) {
  all_vertices <- do.call(rbind, vertices_list)
  lim <- matrix(NA_real_, nrow = 2, ncol = 3)
  colnames(lim) <- colnames(all_vertices)

  for (j in seq_len(3)) {
    vals <- all_vertices[, j]
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      r <- c(-1, 1)
    } else {
      r <- range(vals)
      span <- diff(r)
      if (!is.finite(span) || span <= 0) span <- max(1, abs(r[1]) * 0.1)
      r <- r + c(-1, 1) * span * 0.08
    }
    lim[, j] <- r
  }

  if (!is.null(working_z_max)) {
    lim[2, 3] <- max(diform3d_axis_max(working_z_max, lim[2, 3]), lim[1, 3] + 1e-8)
  }

  lim
}

diform3d_draw_compare_polygons <- function(pm, polygons, fill, border) {
  if (!length(polygons)) return(invisible(NULL))
  depth <- vapply(polygons, function(p) mean(p[, 1] + p[, 2] + p[, 3]), numeric(1))
  polygons <- polygons[order(depth, decreasing = TRUE)]

  for (p in polygons) {
    xy <- grDevices::trans3d(p[, 1], p[, 2], p[, 3], pm)
    graphics::polygon(xy$x, xy$y, col = fill, border = border, lwd = 0.6)
  }
}

diform3d_compare_title <- function(diform_result, name) {
  label <- if (!is.null(diform_result$scenario_label)) diform_result$scenario_label else name
  sprintf("%s\n%s", label, diform3d_sample_label(diform_result))
}

plot_diform_3d_compare_r <- function(diform_results,
                                     out_dir = ".",
                                     rp,
                                     response_label = "Response",
                                     response_short = "Response",
                                     elev = 30,
                                     azim = 60,
                                     working_z_max = NULL,
                                     file_prefix = "diform_3d_compare") {
  if (is.null(names(diform_results))) names(diform_results) <- rep("", length(diform_results))
  missing_names <- names(diform_results) == ""
  names(diform_results)[missing_names] <- paste0("scenario_", which(missing_names))

  has_rp <- vapply(diform_results, function(x) length(diform3d_get_rps(x, rp)) > 0, logical(1))
  diform_results <- diform_results[has_rp]

  vertices_list <- lapply(diform_results, diform3d_vertices, rp = rp)
  lim <- diform3d_compare_limits(vertices_list, working_z_max)
  comp_names <- colnames(vertices_list[[1]])
  if (is.null(comp_names) || length(comp_names) != 3) comp_names <- c("x1", "x2", "x3")
  response_col <- diform3d_response_component_index(diform_results[[1]], vertices_list[[1]])
  comp_labels <- diform3d_component_labels(
    comp_names,
    diform_spec = diform_results[[1]]$diform_spec,
    response_label = response_label,
    response_col = response_col
  )

  layout_dim <- diform3d_compare_layout(length(diform_results))
  cols <- grDevices::hcl.colors(length(diform_results), "Dark 3")
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  path <- file.path(out_dir, sprintf("%s_%syr.png", file_prefix, format(rp, trim = TRUE)))

  diform3d_open_png(
    path,
    width = max(1500, 760 * layout_dim[2]),
    height = max(1100, 680 * layout_dim[1]),
    res = 150
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mfrow = layout_dim, mar = c(4.2, 4.2, 4.2, 1.0), oma = c(0, 0, 3.0, 0), xpd = NA)
  on.exit(par(old), add = TRUE)

  for (i in seq_along(diform_results)) {
    diform_result <- diform_results[[i]]
    vertices <- vertices_list[[i]]
    constraints <- diform3d_constraints(diform_result, rp)
    polygons <- diform3d_working_polygons(vertices, constraints)
    empty_z <- matrix(lim[1, 3], 2, 2)

    pm <- graphics::persp(
      x = lim[, 1],
      y = lim[, 2],
      z = empty_z,
      zlim = lim[, 3],
      col = NA,
      border = NA,
      ticktype = "detailed",
      theta = azim,
      phi = elev,
      r = 7,
      expand = 0.88,
      xlab = comp_labels[1],
      ylab = comp_labels[2],
      zlab = comp_labels[3],
      main = diform3d_compare_title(diform_result, names(diform_results)[i]),
      cex.lab = 0.86,
      cex.axis = 0.70,
      cex.main = 0.88,
      font.main = 2
    )

    diform3d_draw_compare_polygons(
      pm,
      polygons,
      fill = diform3d_alpha(cols[i], 0.30),
      border = diform3d_alpha(cols[i], 0.80)
    )

    response_col <- diform3d_response_component_index(diform_result, vertices)
    if (identical(diform_result$plot_type, "polar_response") &&
        all(c("H_L", "H_T") %in% colnames(vertices))) {
      response_values <- vertices[, "H_L"]^2 + vertices[, "H_T"]^2
      response_col <- which.max(response_values)
      star <- vertices[response_col, ]
    } else {
      star <- vertices[which.max(vertices[, response_col]), ]
    }
    p <- grDevices::trans3d(star[1], star[2], star[3], pm)
    graphics::points(p, pch = 18, col = cols[i], cex = 1.05)
  }

  graphics::mtext(
    sprintf("D-IFORM %s-yr contour comparison in transformed 3D space", format(rp, trim = TRUE)),
    outer = TRUE,
    side = 3,
    line = 1.0,
    cex = 1.1,
    font = 2
  )
  graphics::mtext(
    sprintf("Common axes; marker shows the maximum %s point on each contour.", response_short),
    outer = TRUE,
    side = 3,
    line = -0.8,
    cex = 0.78
  )

  invisible(path)
}

plot_diform_3d_compare_physical_r <- function(diform_results,
                                              out_dir = ".",
                                              rp,
                                              response_label = "Response",
                                              response_short = "Response",
                                              elev = 30,
                                              azim = 60,
                                              n_ang = 180,
                                              n_cspd = 80,
                                              z_max = NULL,
                                              file_prefix = "diform_3d_compare_physical") {
  if (is.null(names(diform_results))) names(diform_results) <- rep("", length(diform_results))
  missing_names <- names(diform_results) == ""
  names(diform_results)[missing_names] <- paste0("scenario_", which(missing_names))

  is_physical <- vapply(diform_results, function(x) {
    identical(x$plot_type, "polar_response")
  }, logical(1))
  diform_results <- diform_results[is_physical]
  if (!length(diform_results)) return(invisible(NULL))

  has_rp <- vapply(diform_results, function(x) length(diform3d_get_rps(x, rp)) > 0, logical(1))
  diform_results <- diform_results[has_rp]

  vertices_list <- lapply(diform_results, diform3d_vertices, rp = rp)
  all_vertices <- do.call(rbind, vertices_list)
  meta <- diform3d_polar_meta(diform_results[[1]])
  limits <- diform3d_limits(all_vertices, obs = NULL, z_max = z_max)
  c_grid <- seq(limits$cspd[1], limits$cspd[2], length.out = n_cspd)
  a_grid <- seq(limits$angle[1], limits$angle[2], length.out = n_ang)

  layout_dim <- diform3d_compare_layout(length(diform_results))
  cols <- grDevices::hcl.colors(length(diform_results), "Dark 3")
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  path <- file.path(out_dir, sprintf("%s_%syr.png", file_prefix, format(rp, trim = TRUE)))

  diform3d_open_png(
    path,
    width = max(1500, 760 * layout_dim[2]),
    height = max(900, 650 * layout_dim[1]),
    res = 150
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mfrow = layout_dim, mar = c(4.2, 4.4, 4.2, 1.0), oma = c(0, 0, 3.0, 0), xpd = NA)
  on.exit(par(old), add = TRUE)

  for (i in seq_along(diform_results)) {
    diform_result <- diform_results[[i]]
    vertices <- vertices_list[[i]]
    constraints <- diform3d_constraints(diform_result, rp)
    hs <- diform3d_surface_physical(constraints, c_grid, a_grid)

    pm <- graphics::persp(
      x = c_grid,
      y = a_grid,
      z = hs,
      xlim = limits$cspd,
      ylim = limits$angle,
      zlim = limits$response,
      col = diform3d_alpha(cols[i], 0.38),
      border = NA,
      shade = NA,
      ticktype = "detailed",
      theta = azim,
      phi = elev,
      r = 7,
      expand = 0.88,
      xlab = diform3d_label_for(meta$linear, meta$diform_spec),
      ylab = diform3d_label_for(meta$angle, meta$diform_spec,
                                angle_unit = meta$angle_unit),
      zlab = response_label,
      main = diform3d_compare_title(diform_result, names(diform_results)[i]),
      cex.lab = 0.84,
      cex.axis = 0.68,
      cex.main = 0.88,
      font.main = 2
    )

    phys <- diform3d_physical(vertices[, 1], vertices[, 2], vertices[, 3], meta = meta)
    i_star <- which.max(phys$response)
    p <- grDevices::trans3d(
      phys$Cspd[i_star],
      phys$angularDifference[i_star],
      phys$response[i_star],
      pm
    )
    graphics::points(p, pch = 18, col = cols[i], cex = 1.05)
  }

  graphics::mtext(
    sprintf("D-IFORM %s-yr contour comparison in physical space", format(rp, trim = TRUE)),
    outer = TRUE,
    side = 3,
    line = 1.0,
    cex = 1.1,
    font = 2
  )
  graphics::mtext(
    sprintf("Common axes; marker shows the maximum %s point on each contour.", response_short),
    outer = TRUE,
    side = 3,
    line = -0.8,
    cex = 0.78
  )

  invisible(path)
}
