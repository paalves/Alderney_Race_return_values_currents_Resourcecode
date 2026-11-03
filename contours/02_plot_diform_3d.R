#' Add transparency to colours.
diform3d_alpha <- function(col, alpha = 1) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

diform3d_tangent_color <- "#7B2CBF"
diform3d_halfspace_color <- "#7A7A7A"
diform3d_obs_color <- "#111111"
diform3d_orange_point <- "#D97706"

#' Subtitle giving the length of the observation record.
diform3d_data_label <- function(diform_result) {
  n_years <- diform_result$n_years
  n_years <- if (is.finite(n_years) && abs(n_years - round(n_years)) < 1e-6) {
    round(n_years)
  } else {
    round(n_years, 2)
  }
  sprintf("%s yr observed data", format(n_years, big.mark = ",", scientific = FALSE))
}

#' Colour palette for the response values.
diform3d_palette <- function(n = 256) {
  grDevices::hcl.colors(n, "Plasma")
}

#' Map values to palette colours between vmin and vmax.
diform3d_color_values <- function(z, vmin, vmax, palette, alpha = 1) {
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) vmax <- vmin + 1
  z <- pmax(vmin, pmin(vmax, z))
  idx <- floor((z - vmin) / (vmax - vmin) * (length(palette) - 1)) + 1
  idx[!is.finite(idx)] <- 1
  diform3d_alpha(palette[idx], alpha)
}

#' Colours of surface cells from the mean of their four corner values.
diform3d_surface_colors <- function(z, vmin, vmax, palette, alpha = 0.74) {
  z_cell <- 0.25 * (
    z[-nrow(z), -ncol(z), drop = FALSE] +
      z[-1, -ncol(z), drop = FALSE] +
      z[-nrow(z), -1, drop = FALSE] +
      z[-1, -1, drop = FALSE]
  )
  diform3d_color_values(z_cell, vmin, vmax, palette, alpha)
}

#' Open a PNG device, falling back to the base device if ragg is unavailable.
diform3d_open_png <- function(path, width = 1520, height = 1200, res = 160) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    ragg::agg_png(path, width = width, height = height, res = res, background = "white")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    grDevices::png(path, width = width, height = height, res = res, bg = "white")
  }
}

#' Return a valid numeric axis limit or the default.
diform3d_axis_max <- function(value, default) {
  if (is.null(value)) return(default)
  value <- as.numeric(value)[1]
  if (!is.finite(value)) return(default)
  value
}

#' Return periods available in a D-IFORM result, optionally restricted to a selection.
diform3d_get_rps <- function(diform_result, rp = "all") {
  available <- as.numeric(names(diform_result$contours_3d))
  available <- available[is.finite(available)]
  if (identical(rp, "all")) return(sort(available))
  sort(intersect(as.numeric(rp), available))
}

#' Warn if directional return values decrease with the return period.
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

#' Directional return values for one return period.
diform3d_rv_for_rp <- function(diform_result, rp) {
  cn <- paste0("T", rp)
  if (cn %in% colnames(diform_result$rv_matrix)) return(diform_result$rv_matrix[, cn])
  k <- match(rp, diform_result$return_periods)
  if (is.na(k)) stop("Return period not found in D-IFORM result: ", rp)
  diform_result$rv_matrix[, k]
}

#' Half-space normals and offsets of the contour in unnormalized working-space coordinates.
diform3d_constraints <- function(diform_result, rp) {
  rv <- diform3d_rv_for_rp(diform_result, rp)
  dirs <- as.matrix(diform_result$directions)
  med <- as.numeric(diform_result$norm_info$med)
  std <- as.numeric(diform_result$norm_info$std)
  normal <- sweep(dirs, 2, std, "/")
  offset <- -as.numeric(normal %*% med) - rv
  list(normal = normal, offset = offset, rv = rv, med = med, std = std, dirs = dirs)
}

#' Contour vertices in working-space coordinates for one return period.
diform3d_vertices <- function(diform_result, rp) {
  x <- diform_result$contours_3d[[as.character(rp)]]$vertices_3d
  nm <- diform_result$component_names
  if (is.null(nm) || length(nm) != 3) nm <- c("Cspd", "H_L", "H_T")
  colnames(x) <- nm
  as.matrix(x)
}

#' Variable names and angular settings of the linear and circular components.
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

#' Value of a per-variable field of the specification, or NULL.
diform3d_spec_field <- function(diform_spec, field, variable) {
  if (is.null(diform_spec) || is.null(diform_spec[[field]])) return(NULL)
  x <- diform_spec[[field]]
  if (!is.null(x[[variable]])) return(x[[variable]])
  NULL
}

#' Readable name of a variable.
diform3d_pretty_name <- function(variable) {
  if (identical(variable, "Cspd")) return("Current speed")
  if (identical(variable, "angularDifference") || identical(variable, "angleDiff")) {
    return("Angular difference")
  }
  variable
}

#' Default unit of a variable for axis labels or annotated values.
diform3d_default_unit <- function(variable, angle_unit = "degrees", for_value = FALSE) {
  if (identical(variable, "Cspd")) return("m s\u207b\u00b9")
  if (identical(variable, "angularDifference") || identical(variable, "angleDiff")) {
    if (identical(angle_unit, "radians")) return("rad")
    return(if (isTRUE(for_value)) "\u00B0" else "deg")
  }
  NULL
}

#' Replace degree spellings by the degree symbol.
diform3d_normalise_value_unit <- function(unit) {
  if (is.null(unit) || !nzchar(unit)) return(NULL)
  if (unit %in% c("deg", "degree", "degrees")) return("\u00B0")
  unit
}

#' Axis label of a variable with its unit.
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

#' Display name of a variable.
diform3d_display_name_for <- function(variable, diform_spec = NULL) {
  label <- diform3d_spec_field(diform_spec, "variable_labels", variable)
  if (!is.null(label)) label else diform3d_pretty_name(variable)
}

#' Axis labels and title of the working-space plot.
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

#' Override the names of the two Cartesian response components.
diform3d_set_working_component_labels <- function(meta, labels = NULL) {
  if (is.null(labels)) return(meta)
  if (length(labels) != 2 || anyNA(labels) || any(!nzchar(labels))) {
    stop("working_component_labels must contain two non-empty labels.", call. = FALSE)
  }
  if (is.null(meta$h_names) || length(meta$h_names) != 2) {
    meta$h_names <- c("H_L", "H_T")
  }
  if (is.null(meta$diform_spec$variable_labels)) {
    meta$diform_spec$variable_labels <- list()
  }
  meta$diform_spec$variable_labels[meta$h_names] <- as.list(as.character(labels))
  meta
}

#' Formatted 'name = value unit' line for the annotation box.
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

#' Convert working-space coordinates to current speed, angle and response.
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

#' Observations in physical coordinates, flagged as declustered peaks when available.
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

#' Observations in working-space coordinates.
diform3d_obs_working <- function(obs) {
  a <- obs$angularDifference * pi / 180
  h <- sqrt(pmax(obs$response, 0))
  cbind(Cspd = obs$Cspd, H_L = h * cos(a), H_T = h * sin(a))
}

#' Signed values of every half-space constraint at the given points (positive outside).
diform3d_signed <- function(pts, constraints) {
  as.matrix(pts) %*% t(constraints$normal) +
    matrix(constraints$offset, nrow = nrow(pts),
           ncol = length(constraints$offset), byrow = TRUE)
}

#' Flag observations lying outside the contour.
diform3d_outside <- function(obs, constraints, tol = 1e-10) {
  if (is.null(obs) || !nrow(obs)) return(logical())
  pts <- diform3d_obs_working(obs)
  apply(diform3d_signed(pts, constraints) > tol, 1, any)
}

#' Flag observations lying beyond the supporting plane.
diform3d_halfspace_mask <- function(obs, support_info, tol = 1e-10) {
  if (is.null(obs) || !nrow(obs)) return(logical())
  pts <- diform3d_obs_working(obs)
  as.numeric(pts %*% support_info$normal + support_info$offset) > tol
}

#' Masks of the peaks outside the contour and beyond the supporting plane.
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

#' Union of the highlight masks.
diform3d_peak_mask <- function(masks) {
  orange <- masks$orange
  green <- masks$green
  if (is.null(orange)) orange <- logical()
  if (is.null(green)) green <- logical()
  if (!length(orange)) return(green)
  if (!length(green)) return(orange)
  orange | green
}

#' Axis limits of the physical and working-space plots.
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
  n <- sqrt(sum(x^2))
  if (!is.finite(n) || n == 0) return(x)
  x / n
}

#' Cross product of two three-dimensional vectors.
diform3d_cross <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' Two orthonormal vectors spanning the plane orthogonal to a normal.
diform3d_plane_basis <- function(normal) {
  n <- diform3d_unit(normal)
  a <- if (abs(n[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  e1 <- a - sum(a * n) * n
  e1 <- diform3d_unit(e1)
  e2 <- diform3d_unit(diform3d_cross(n, e1))
  list(e1 = e1, e2 = e2)
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
  dirs_all <- do.call(rbind, lapply(seq_len(nrow(signs)), function(i) {
    sweep(base, 2, signs[i, ], "*")
  }))
  dirs_all <- round(dirs_all / sqrt(rowSums(dirs_all^2)), 10)
  dirs <- unique(dirs_all)
  colnames(dirs) <- c("d1", "d2", "d3")
  key_to_index <- setNames(seq_len(nrow(dirs)), diform3d_direction_key(dirs))

  #' Key of the normalized direction for simplex indices i, j and a sign pattern.
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

#' Triangulation of a direction set, recovered by matching it to a simplex mesh.
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

#' Supporting plane of the contour at the vertex with the largest response.
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

#' Grid of the supporting-plane boundary mapped to physical coordinates around the maximum.
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

#' Draw a grid of quadrilaterals on a perspective plot.
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

#' Draw a vertical colour bar.
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

#' Draw a framed block of text lines.
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

#' Place the annotation box at the top left of the plot.
diform3d_add_annotation <- function(lines) {
  usr <- par("usr")
  x <- usr[1] - 0.090 * diff(usr[1:2])
  y <- usr[4] + 0.025 * diff(usr[3:4])
  diform3d_text_box(x, y, lines, cex = 1.25)
}

#' Randomly subsample the background observations.
diform3d_limit_observed_points <- function(obs, max_background_points = 0) {
  if (is.null(obs) || !nrow(obs)) return(obs)
  if (is.null(max_background_points) || max_background_points <= 0 ||
      nrow(obs) <= max_background_points) {
    return(obs)
  }
  indices <- unique(round(seq(1, nrow(obs), length.out = max_background_points)))
  obs[indices, , drop = FALSE]
}

#' Draw observations and highlighted peaks on the physical-space plot.
diform3d_draw_points_physical <- function(pm, obs, masks, limits, max_labels = 6,
                                          show_peak_highlights = TRUE,
                                          obs_background = obs,
                                          response_short = "Response") {
  if (is.null(obs) || !nrow(obs)) return(invisible(NULL))
  if (!is.null(obs_background) && nrow(obs_background)) {
    xyz <- grDevices::trans3d(obs_background$Cspd, obs_background$angularDifference,
                              obs_background$response, pm)
    graphics::points(xyz, pch = 16, cex = 0.17,
                     col = diform3d_alpha(diform3d_obs_color, 0.040))
  }
  if (!isTRUE(show_peak_highlights)) return(invisible(NULL))

  #' Draw the observations selected by a mask.
  draw_mask <- function(mask, bg, cex) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    p <- grDevices::trans3d(obs$Cspd[mask], obs$angularDifference[mask],
                            obs$response[mask], pm)
    graphics::points(p, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }
  peak_mask <- diform3d_peak_mask(masks)
  draw_mask(peak_mask, diform3d_orange_point, 0.90)

  #' Label the largest observations selected by a mask.
  label_mask <- function(mask, col) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    idx <- which(mask)
    idx <- idx[order(obs$response[idx], decreasing = TRUE)]
    idx <- idx[seq_len(min(max_labels, length(idx)))]
    p <- grDevices::trans3d(obs$Cspd[idx], obs$angularDifference[idx],
                            obs$response[idx] + 0.012 * limits$response[2], pm)
    graphics::text(p$x, p$y, sprintf("%s=%.2f", response_short, obs$response[idx]),
                   cex = 0.70, col = col, pos = 3)
  }
  label_mask(peak_mask, "darkorange")
  invisible(NULL)
}

#' Draw observations and highlighted peaks on the working-space plot.
diform3d_draw_points_working <- function(pm, obs, masks, limits, max_labels = 6,
                                         show_peak_highlights = TRUE,
                                         obs_background = obs,
                                         response_short = "Response") {
  if (is.null(obs) || !nrow(obs)) return(invisible(NULL))
  pts <- diform3d_obs_working(obs)
  if (!is.null(obs_background) && nrow(obs_background)) {
    bg_pts <- diform3d_obs_working(obs_background)
    xyz <- grDevices::trans3d(bg_pts[, 1], bg_pts[, 2], bg_pts[, 3], pm)
    graphics::points(xyz, pch = 16, cex = 0.17,
                     col = diform3d_alpha(diform3d_obs_color, 0.040))
  }
  if (!isTRUE(show_peak_highlights)) return(invisible(NULL))

  #' Draw the observations selected by a mask.
  draw_mask <- function(mask, bg, cex) {
    if (!length(mask) || !any(mask)) return(invisible(NULL))
    p <- grDevices::trans3d(pts[mask, 1], pts[mask, 2], pts[mask, 3], pm)
    graphics::points(p, pch = 21, cex = cex, col = "black", bg = bg, lwd = 0.45)
  }
  peak_mask <- diform3d_peak_mask(masks)
  draw_mask(peak_mask, diform3d_orange_point, 0.90)

  #' Label the largest observations selected by a mask.
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
    graphics::text(p$x, p$y, sprintf("%s=%.2f", response_short, hs[idx]),
                   cex = 0.70, col = col, pos = 3)
  }
  label_mask(peak_mask, "darkorange")
  invisible(NULL)
}

#' Draw the legend of the point symbols.
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

#' Draw the contour surface in physical coordinates and save it as a PNG file.
diform3d_plot_physical <- function(path, vertices, constraints, support_info,
                                   all_vertices, obs, masks,
                                   response_label, response_short, response_unit,
                                   n_years, data_label, rp, elev, azim,
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
                   data_label),
    cex.lab = 1.34,
    cex.axis = 1.06,
    cex.main = 1.42,
    font.main = 2
  )
  graphics::mtext("(b)", side = 3, adj = 0, line = 2.6, font = 2, cex = 1.6)

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

  obs_bg <- diform3d_limit_observed_points(obs, max_background_points)
  if (isTRUE(show_points)) {
    diform3d_draw_points_physical(pm, obs, masks, limits, max_labels,
                                  show_peak_highlights = show_peak_highlights,
                                  obs_background = obs_bg,
                                  response_short = response_short)
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

#' Contour faces as polygons of the vertices lying on each active half-space.
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

#' Polygons of the working-space contour surface, from the direction triangulation when available.
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

#' Draw polygons back to front, coloured by response.
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

#' Draw the contour in working-space coordinates and save it as a PNG file.
diform3d_plot_working <- function(path, vertices, constraints, support_info,
                                  all_vertices, obs, masks,
                                  response_label, response_short, response_unit,
                                  n_years, data_label, rp, elev, azim,
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
                   data_label),
    cex.lab = 1.28,
    cex.axis = 1.04,
    cex.main = 1.34,
    font.main = 2
  )
  graphics::mtext("(a)", side = 3, adj = 0, line = 2.6, font = 2, cex = 1.6)

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

  obs_bg <- diform3d_limit_observed_points(obs, max_background_points)
  if (isTRUE(show_points)) {
    diform3d_draw_points_working(pm, obs, masks, limits, max_labels,
                                 show_peak_highlights = show_peak_highlights,
                                 obs_background = obs_bg,
                                 response_short = response_short)
  }

  if (isTRUE(show_points) && isTRUE(show_peak_highlights) && !is.null(obs) && nrow(obs)) {
    diform3d_add_legend(show_support_plane)
  }
  diform3d_draw_colorbar(0, vmax, palette, response_label)
  invisible(path)
}

#' Polygon where the supporting plane cuts the plotting box, and the box faces on its outer side.
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

  #' Part of the box face at a given axis value, as a convex polygon.
  face_from_axis <- function(axis, value) {
    idx <- which(abs(pts[, axis] - value) < 1e-8)
    if (length(idx) < 3) return(NULL)
    other <- setdiff(seq_len(3), axis)
    h <- grDevices::chull(pts[idx, other[1]], pts[idx, other[2]])
    pts[idx[h], , drop = FALSE]
  }

  #' Polygon where the supporting plane cuts the box.
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

#' Draw a list of faces back to front.
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

#' Draw the physical and working-space contour figures for the selected return periods.
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
                             working_component_labels = NULL,
                             show_exterior_halfspace = FALSE,
                             make_halfspace_plot = TRUE,
                             make_working_plot = TRUE) {
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  response <- diform_result$response
  obs <- diform3d_obs(obs_data, excess_df, response, diform_result$diform_spec)
  meta <- diform3d_polar_meta(diform_result)
  meta <- diform3d_set_working_component_labels(meta, working_component_labels)
  rps <- diform3d_get_rps(diform_result, rp)
  if (!length(rps)) stop("No D-IFORM return period available for rp = ", paste(rp, collapse = ", "))

  all_rps <- diform3d_get_rps(diform_result, "all")
  all_vertices <- do.call(rbind, lapply(all_rps, function(x) diform3d_vertices(diform_result, x)))
  files <- character()
  data_label <- diform3d_data_label(diform_result)
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
      data_label = data_label,
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

    if (isTRUE(make_working_plot)) {
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
        data_label = data_label,
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
  }

  invisible(files)
}

