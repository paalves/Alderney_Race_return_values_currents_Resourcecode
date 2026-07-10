library(data.table)
library(ggplot2)
library(viridisLite)
library(Rcpp)


sourceCpp(file.path(root, "pot_decluster.cpp"))

default_diform_spec <- function(response = "hs") {
  list(
    components = list(
      list(type = "linear", variable = "Cspd", name = "Cspd"),
      list(type = "circular_radius", radius = response, angle = "angularDifference",
           transform = "sqrt", names = c("H_L", "H_T"),
           angle_unit = "degrees", period = 360, origin = 0)
    ),
    response = response,
    response_component = NULL,
    plot_type = "polar_response"
  )
}

directional_diform_spec <- function(linear, angle, response = "hs",
                                    angle_unit = "degrees",
                                    period = NULL,
                                    origin = NULL,
                                    transform = "sqrt") {
  list(
    components = list(
      list(type = "linear", variable = linear, name = linear),
      list(
        type = "circular_radius",
        radius = response,
        angle = angle,
        transform = transform,
        names = c("H_L", "H_T"),
        angle_unit = angle_unit,
        period = period,
        origin = origin
      )
    ),
    response = response,
    response_component = NULL,
    plot_type = "polar_response"
  )
}

linear_diform_spec <- function(variables, response = "hs", response_component = response) {
  list(
    components = lapply(variables, function(v) list(type = "linear", variable = v, name = v)),
    response = response,
    response_component = response_component,
    plot_type = "linear"
  )
}

normalise_diform_spec <- function(diform_spec = NULL, response = "hs") {
  if (is.null(diform_spec)) diform_spec <- default_diform_spec(response)
  if (is.character(diform_spec)) diform_spec <- linear_diform_spec(diform_spec, response)
  if (is.null(diform_spec$response)) diform_spec$response <- response
  diform_spec
}

diform_component_names <- function(diform_spec) {
  unlist(lapply(diform_spec$components, function(comp) {
    if (identical(comp$type, "circular_radius")) return(comp$names)
    if (!is.null(comp$name)) comp$name else comp$variable
  }), use.names = FALSE)
}

transform_radius <- function(x, transform = "identity") {
  if (identical(transform, "sqrt")) return(sqrt(pmax(x, 0)))
  if (identical(transform, "log")) return(log(x))
  x
}

inverse_radius <- function(x, transform = "identity") {
  if (identical(transform, "sqrt")) return(x^2)
  if (identical(transform, "log")) return(exp(x))
  x
}

angle_to_radians <- function(x, unit = "degrees") {
  if (identical(unit, "radians")) return(x)
  x * pi / 180
}

angle_from_radians <- function(x, unit = "degrees", period = NULL, origin = NULL) {
  out <- if (identical(unit, "radians")) x else x * 180 / pi
  if (is.null(period)) period <- if (identical(unit, "radians")) 2 * pi else 360
  if (is.null(origin)) origin <- 0
  ((out - origin) %% period) + origin
}

transform_to_diform_space <- function(dt, response = "hs", diform_spec = NULL) {
  diform_spec <- normalise_diform_spec(diform_spec, response)
  cols <- list()

  for (comp in diform_spec$components) {
    if (identical(comp$type, "linear")) {
      if (!comp$variable %in% names(dt)) {
        stop("D-IFORM variable not found in data: ", comp$variable,
             ". Check cfg$diform_spec against cfg$covariates.", call. = FALSE)
      }
      nm <- if (!is.null(comp$name)) comp$name else comp$variable
      cols[[nm]] <- dt[[comp$variable]]
    } else if (identical(comp$type, "circular_radius")) {
      missing <- setdiff(c(comp$angle, comp$radius), names(dt))
      if (length(missing)) {
        stop("D-IFORM variable(s) not found in data: ", paste(missing, collapse = ", "),
             ". Check cfg$diform_spec against cfg$covariates.", call. = FALSE)
      }
      a <- angle_to_radians(dt[[comp$angle]], comp$angle_unit)
      r <- transform_radius(dt[[comp$radius]], comp$transform)
      cols[[comp$names[1]]] <- r * cos(a)
      cols[[comp$names[2]]] <- r * sin(a)
    } else {
      stop("Unsupported D-IFORM component type: ", comp$type)
    }
  }

  X <- as.matrix(as.data.frame(cols))
  storage.mode(X) <- "double"
  if (ncol(X) != 3) stop("This D-IFORM implementation currently requires exactly 3 transformed components.")
  X
}

normalise_diform <- function(X) {
  med <- apply(X, 2, median)
  std <- apply(X, 2, sd)
  list(Y = sweep(sweep(X, 2, med), 2, std, "/"), med = med, std = std)
}

unnormalise_diform <- function(Y, norm_info) sweep(sweep(Y, 2, norm_info$std, "*"), 2, norm_info$med, "+")

backtransform_to_original <- function(pts, response = "hs", diform_spec = NULL) {
  diform_spec <- normalise_diform_spec(diform_spec, response)
  out <- list()
  j <- 1L

  for (comp in diform_spec$components) {
    if (identical(comp$type, "linear")) {
      nm <- if (!is.null(comp$name)) comp$name else comp$variable
      out[[comp$variable]] <- pts[, j]
      if (!identical(nm, comp$variable)) out[[nm]] <- pts[, j]
      j <- j + 1L
    } else if (identical(comp$type, "circular_radius")) {
      x <- pts[, j]
      y <- pts[, j + 1L]
      out[[comp$angle]] <- angle_from_radians(
        atan2(y, x),
        unit = comp$angle_unit,
        period = comp$period,
        origin = comp$origin
      )
      out[[comp$radius]] <- inverse_radius(sqrt(x^2 + y^2), comp$transform)
      j <- j + 2L
    }
  }

  as.data.frame(out)
}

generate_direction_vectors_3d <- function(simplex_resolution = 20) {
  w <- (0:simplex_resolution) / simplex_resolution
  x <- do.call(rbind, lapply(w, function(w1) do.call(rbind, lapply(w, function(w2) {
    w3 <- 1 - w1 - w2
    if (w3 >= -1e-10) c(w1, w2, max(0, w3))
  }))))
  s <- as.matrix(expand.grid(c(1, -1), c(1, -1), c(1, -1)))
  x <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) sweep(x, 2, s[i, ], "*")))
  x <- x[rowSums(x^2) > 1e-10, ]
  x <- unique(round(x / sqrt(rowSums(x^2)), 10))
  colnames(x) <- c("d1", "d2", "d3")
  x
}

fit_gp_constrained <- function(exc, max_obs = NULL) {
  exc <- exc[exc > 0]
  n <- length(exc)
  if (n < 3) return(list(sigma = mean(exc), xi = -0.1))
  m <- mean(exc)
  xi0 <- max(-0.99, min(-0.01, 0.5 * (1 - m^2 / var(exc))))
  sig0 <- max(0.01, m * (1 - xi0))
  nll <- function(par) {
    sig <- par[1]
    xi <- par[2]
    z <- 1 + xi * exc / sig
    if (sig <= 0 || xi >= 0 || xi < -1 || any(z <= 0) || (!is.null(max_obs) && max_obs > -sig / xi)) 1e12
    else n * log(sig) + (1 + 1 / xi) * sum(log(z))
  }
  o <- optim(c(sig0, xi0), nll, method = "Nelder-Mead", control = list(maxit = 2000))
  if (o$convergence) list(sigma = sig0, xi = xi0) else list(sigma = o$par[1], xi = o$par[2])
}

gp_return_value <- function(sigma, xi, lambda, T_years) {
  if (abs(xi) < 1e-8) sigma * log(lambda * T_years) else sigma / xi * ((lambda * T_years)^xi - 1)
}

compute_return_values_pot <- function(Y, dirs, return_periods, n_years, h_hours = 72, zeta = 0.1, n_threads = 0) {
  dec <- pot_decluster_all(Y, dirs, h_hours, zeta, n_threads)
  rv <- matrix(NA_real_, nrow(dirs), length(return_periods), dimnames = list(NULL, paste0("T", return_periods)))
  fits <- data.frame(mu = dec$mu, sigma = NA_real_, xi = NA_real_, n_peaks = dec$n_peaks,
                     n_exc = lengths(dec$exceedances), lambda = lengths(dec$exceedances) / n_years,
                     method = NA_character_)

  for (j in seq_len(nrow(dirs))) {
    peaks <- dec$peaks[[j]]
    exc <- dec$exceedances[[j]]
    if (length(exc) < 5 || fits$lambda[j] < 0.01) {
      fits$method[j] <- "empirical"
      lam <- length(peaks) / n_years
      rv[j, ] <- sapply(return_periods, function(T) as.numeric(quantile(peaks, pmin(1, pmax(0, 1 - 1 / (lam * T))), type = 7)))
    } else {
      gp <- fit_gp_constrained(exc, max(peaks) - fits$mu[j])
      fits$sigma[j] <- gp$sigma
      fits$xi[j] <- gp$xi
      fits$method[j] <- "gpd"
      rv[j, ] <- fits$mu[j] + sapply(return_periods, gp_return_value, sigma = gp$sigma, xi = gp$xi, lambda = fits$lambda[j])
    }
  }
  list(rv = rv, fits = fits, peaks = dec$peaks, exceedances = dec$exceedances)
}

build_diform_contours_3d <- function(dirs, rv_matrix, norm_info, return_periods,
                                     response = "hs", diform_spec = NULL) {
  dots <- dirs %*% t(dirs)
  out <- list()
  comp_names <- diform_component_names(normalise_diform_spec(diform_spec, response))

  for (k in seq_along(return_periods)) {
    rv <- rv_matrix[, k]
    pts <- t(vapply(seq_len(nrow(dirs)), function(j) {
      pos <- dots[, j] > 1e-10
      min(rv[pos] / dots[pos, j]) * dirs[j, ]
    }, numeric(3)))
    pts <- pts[complete.cases(pts), , drop = FALSE]
    if (nrow(pts) >= 4) {
      orig <- unnormalise_diform(pts, norm_info)
      colnames(orig) <- comp_names
      out[[as.character(return_periods[k])]] <- list(vertices_norm = pts, vertices_3d = orig,
                                                     vertices_bt = backtransform_to_original(orig, response, diform_spec),
                                                     rp = return_periods[k])
    }
  }
  out
}

run_diform_contours <- function(sim, response = "hs", n_years,
                                return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                simplex_resolution = 20, h_hours = 72, zeta = 0.1, n_threads = 0,
                                diform_spec = NULL) {
  diform_spec <- normalise_diform_spec(diform_spec, response)
  X <- transform_to_diform_space(sim, response, diform_spec)
  norm <- normalise_diform(X)
  dirs <- generate_direction_vectors_3d(simplex_resolution)
  pot <- compute_return_values_pot(norm$Y, dirs, return_periods, n_years, h_hours, zeta, n_threads)
  contours <- build_diform_contours_3d(dirs, pot$rv, norm, return_periods, response, diform_spec)
  list(contours_3d = contours, directions = dirs, rv_matrix = pot$rv, fits = pot$fits,
       peaks = pot$peaks, exceedances = pot$exceedances, norm_info = norm, X_transformed = X,
       summary = diform_contour_summary(contours, response), response = response,
       n_years = n_years, return_periods = return_periods,
       diform_spec = diform_spec,
       component_names = colnames(X),
       plot_type = diform_spec$plot_type)
}

diform_contour_summary <- function(contours_3d, response = "hs") {
  rbindlist(lapply(contours_3d, function(cc) {
    bt <- cc$vertices_bt
    if (response %in% names(bt)) {
      i <- which.max(bt[[response]])
      out <- data.table(return_period = cc$rp, n_vertices = nrow(bt),
                        max_resp = round(max(bt[[response]], na.rm = TRUE), 3))
      for (nm in names(bt)) out[[paste0(nm, "_at_max")]] <- round(bt[[nm]][i], 3)
      out
    } else {
      data.table(return_period = cc$rp, n_vertices = nrow(bt))
    }
  }))
}

diform_sample_label <- function(diform_result) {
  n_years <- diform_result$n_years
  n_years <- if (is.finite(n_years) && abs(n_years - round(n_years)) < 1e-6) {
    round(n_years)
  } else {
    round(n_years, 2)
  }
  source <- if (isFALSE(diform_result$bootstrap)) "observed data" else "simulation"
  sprintf("%s yr %s", format(n_years, big.mark = ",", scientific = FALSE), source)
}

hull_boundary_2d <- function(x, y) {
  i <- c(chull(x, y), chull(x, y)[1])
  data.frame(x = x[i], y = y[i])
}

diform_plot_spec_field <- function(diform_spec, field, variable) {
  if (is.null(diform_spec) || is.null(diform_spec[[field]])) return(NULL)
  x <- diform_spec[[field]]
  if (!is.null(x[[variable]])) x[[variable]] else NULL
}

diform_plot_pretty_name <- function(variable) {
  if (identical(variable, "Cspd")) return("Current speed")
  if (identical(variable, "angularDifference") || identical(variable, "angleDiff")) {
    return("Angular difference")
  }
  variable
}

diform_plot_label_for <- function(variable, diform_spec = NULL,
                                  is_angle = FALSE,
                                  angle_unit = "degrees") {
  label <- diform_plot_spec_field(diform_spec, "variable_labels", variable)
  unit <- diform_plot_spec_field(diform_spec, "variable_units", variable)
  if ((is.null(unit) || !nzchar(unit)) && isTRUE(is_angle)) {
    unit <- if (identical(angle_unit, "radians")) "rad" else "deg"
  }
  if ((is.null(unit) || !nzchar(unit)) && identical(variable, "Cspd")) {
    unit <- "m.s^-1"
  }

  name <- if (!is.null(label)) label else diform_plot_pretty_name(variable)
  if (is.null(unit) || !nzchar(unit)) return(name)
  paste0(name, " [", unit, "]")
}

diform_polar_meta <- function(diform_result, response = NULL) {
  response <- if (is.null(response)) diform_result$response else response
  diform_spec <- normalise_diform_spec(diform_result$diform_spec, response)
  is_linear <- vapply(diform_spec$components, function(comp) {
    identical(comp$type, "linear")
  }, logical(1))
  is_circular <- vapply(diform_spec$components, function(comp) {
    identical(comp$type, "circular_radius")
  }, logical(1))

  if (sum(is_linear) != 1 || sum(is_circular) != 1) {
    stop("Polar D-IFORM projection needs one linear component and one circular-radius component.",
         call. = FALSE)
  }

  linear <- diform_spec$components[[which(is_linear)]]
  circular <- diform_spec$components[[which(is_circular)]]
  angle_unit <- if (is.null(circular$angle_unit)) "degrees" else circular$angle_unit
  period <- if (is.null(circular$period)) {
    if (identical(angle_unit, "radians")) 2 * pi else 360
  } else {
    circular$period
  }
  origin <- if (is.null(circular$origin)) 0 else circular$origin
  transform <- if (is.null(circular$transform)) "identity" else circular$transform

  list(
    linear = linear$variable,
    linear_name = if (!is.null(linear$name)) linear$name else linear$variable,
    angle = circular$angle,
    response = circular$radius,
    angle_unit = angle_unit,
    period = period,
    origin = origin,
    transform = transform,
    h_names = circular$names,
    diform_spec = diform_spec
  )
}

diform_observed_polar <- function(obs_data, meta) {
  if (is.null(obs_data) || !nrow(obs_data)) return(NULL)

  out <- data.frame(
    linear = obs_data[[meta$linear]],
    angle = ((obs_data[[meta$angle]] - meta$origin) %% meta$period) + meta$origin,
    resp = obs_data[[meta$response]]
  )
  out <- out[complete.cases(out), , drop = FALSE]
  rownames(out) <- NULL
  out
}

diform_angle_breaks <- function(meta) {
  if (identical(meta$angle_unit, "radians")) {
    return(seq(meta$origin, meta$origin + meta$period, length.out = 9))
  }
  by <- if (abs(meta$period - 360) < 1e-8) 45 else meta$period / 8
  seq(meta$origin, meta$origin + meta$period, by = by)
}

polar_envelope <- function(H_L, H_T, n = 361,
                           angle_unit = "degrees",
                           period = NULL,
                           origin = NULL,
                           transform = "sqrt") {
  h <- chull(H_L, H_T)
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

  angle <- angle_from_radians(a * pi / 180, unit = angle_unit,
                              period = period, origin = origin)
  out <- data.frame(angDiff = angle, resp = inverse_radius(r, transform))
  out[order(out$angDiff), , drop = FALSE]
}

plot_diform_projections <- function(diform_result,
                                    obs_data = NULL,
                                    response = "hs",
                                    resp_label = "Hs [m]") {
  meta <- diform_polar_meta(diform_result, response)
  obs <- diform_observed_polar(obs_data, meta)
  contours <- diform_result$contours_3d
  rp_vals <- diform_result$return_periods
  rp_colors <- viridisLite::turbo(length(rp_vals))
  names(rp_colors) <- paste0(rp_vals, " yr")
  linear_label <- diform_plot_label_for(meta$linear, meta$diform_spec)
  angle_label <- diform_plot_label_for(meta$angle, meta$diform_spec,
                                       is_angle = TRUE,
                                       angle_unit = meta$angle_unit)
  angle_breaks <- diform_angle_breaks(meta)

  theme_diform <- theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "right",
          panel.grid.minor = element_blank())

  cdf1 <- do.call(rbind, lapply(names(contours), function(rk) {
    v <- contours[[rk]]$vertices_3d
    resp_v <- inverse_radius(sqrt(v[, 2]^2 + v[, 3]^2), meta$transform)
    hb <- hull_boundary_2d(v[, 1], resp_v)
    data.frame(linear = hb$x, resp = hb$y, rp = paste0(rk, " yr"))
  }))
  cdf1$rp <- factor(cdf1$rp, levels = paste0(rp_vals, " yr"))

  p1 <- ggplot()
  if (!is.null(obs)) {
    p1 <- p1 + geom_point(
      data = obs,
      aes(linear, resp), colour = "black", size = 0.8, alpha = 1)
  }
  p1 <- p1 +
    geom_path(data = cdf1, aes(linear, resp, colour = rp), linewidth = 0.9) +
    scale_colour_manual(values = rp_colors, name = "Return period") +
    labs(title = sprintf("D-IFORM: %s vs %s (%s)",
                         response, diform_plot_pretty_name(meta$linear),
                         diform_sample_label(diform_result)),
         x = linear_label,
         y = resp_label) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_diform

  cdf2 <- do.call(rbind, lapply(names(contours), function(rk) {
    v <- contours[[rk]]$vertices_3d
    hb <- hull_boundary_2d(v[, 2], v[, 3])
    data.frame(H_L = hb$x, H_T = hb$y, rp = paste0(rk, " yr"))
  }))
  cdf2$rp <- factor(cdf2$rp, levels = paste0(rp_vals, " yr"))

  p2 <- ggplot()
  if (!is.null(obs)) {
    ang <- angle_to_radians(obs$angle, meta$angle_unit)
    radius <- transform_radius(obs$resp, meta$transform)
    p2 <- p2 + geom_point(
      data = data.frame(H_L = radius * cos(ang), H_T = radius * sin(ang)),
      aes(H_L, H_T), colour = "black", size = 0.8, alpha = 1)
  }
  p2 <- p2 +
    geom_path(data = cdf2, aes(H_L, H_T, colour = rp), linewidth = 0.9) +
    scale_colour_manual(values = rp_colors, name = "Return period") +
    labs(title = sprintf("D-IFORM: transformed %s decomposition by %s",
                         response, diform_plot_pretty_name(meta$angle)),
         x = diform_plot_label_for(meta$h_names[1], meta$diform_spec),
         y = diform_plot_label_for(meta$h_names[2], meta$diform_spec)) +
    coord_equal() +
    theme_diform

  cdf3 <- do.call(rbind, lapply(names(contours), function(rk) {
    v <- contours[[rk]]$vertices_3d
    df <- polar_envelope(v[, 2], v[, 3],
                         angle_unit = meta$angle_unit,
                         period = meta$period,
                         origin = meta$origin,
                         transform = meta$transform)
    df$rp <- paste0(rk, " yr")
    df
  }))
  cdf3$rp <- factor(cdf3$rp, levels = paste0(rp_vals, " yr"))

  p3 <- ggplot()
  if (!is.null(obs)) {
    p3 <- p3 + geom_point(
      data = data.frame(angDiff = obs$angle, resp = obs$resp),
      aes(angDiff, resp), colour = "black", size = 0.8, alpha = 1)
  }
  p3 <- p3 +
    geom_path(data = cdf3, aes(angDiff, resp, colour = rp),
              linewidth = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = rp_colors, name = "Return period") +
    scale_x_continuous(breaks = angle_breaks) +
    labs(title = sprintf("D-IFORM: %s vs angular difference (%s)",
                         response, diform_sample_label(diform_result)),
         x = angle_label,
         y = resp_label) +
    coord_cartesian(xlim = c(meta$origin, meta$origin + meta$period), ylim = c(0, NA)) +
    theme_diform

  max_rp_key <- as.character(max(rp_vals))
  p4 <- NULL

  if (!is.null(contours[[max_rp_key]])) {
    v <- contours[[max_rp_key]]$vertices_3d
    df4 <- na.omit(polar_envelope(v[, 2], v[, 3],
                                  angle_unit = meta$angle_unit,
                                  period = meta$period,
                                  origin = meta$origin,
                                  transform = meta$transform))
    max_pt <- df4[which.max(df4$resp), ]

    p4 <- ggplot()
    if (!is.null(obs)) {
      p4 <- p4 + geom_point(
        data = data.frame(angDiff = obs$angle, resp = obs$resp),
        aes(angDiff, resp), colour = "black", size = 0.8, alpha = 1)
    }
    p4 <- p4 +
      geom_path(data = df4, aes(angDiff, resp),
                colour = rp_colors[length(rp_colors)], linewidth = 1.2) +
      geom_point(data = max_pt, aes(angDiff, resp),
                 colour = "red", shape = 18, size = 5) +
      geom_text(data = max_pt, aes(
        angDiff, resp,
        label = sprintf("Max %s = %.3f\nangDiff = %.0f deg",
                        resp_label, resp, angDiff %% 360)),
        vjust = -0.5, size = 5, colour = "black") +
      scale_x_continuous(breaks = angle_breaks) +
      coord_cartesian(xlim = c(meta$origin, meta$origin + meta$period), ylim = c(0, NA)) +
      labs(title = sprintf("D-IFORM: %s-yr %s envelope vs angular difference",
                           max_rp_key, resp_label),
           x = angle_label,
           y = resp_label) +
      theme_diform
  }

  list(resp_cspd = p1,
       hl_ht = p2,
       resp_angdiff = p3,
       resp_angdiff_max = p4)
}

plot_diform_generic_projections <- function(diform_result,
                                            obs_data = NULL,
                                            response = "hs",
                                            resp_label = "Response") {
  contours <- diform_result$contours_3d
  contour_rps <- as.numeric(names(contours))
  contour_rps <- contour_rps[is.finite(contour_rps)]
  comp_names <- diform_result$component_names
  if (is.null(comp_names) || length(comp_names) != 3) {
    comp_names <- diform_component_names(diform_result$diform_spec)
  }
  comp_labels <- vapply(comp_names, diform_plot_label_for, character(1),
                        diform_spec = diform_result$diform_spec)
  response_component <- diform_result$diform_spec$response_component
  if (is.null(response_component)) response_component <- response
  response_idx <- which(comp_names == response_component)
  if (length(response_idx) && !is.null(resp_label) && nzchar(resp_label)) {
    comp_labels[response_idx] <- resp_label
  }

  rp_colors <- viridisLite::turbo(length(contour_rps))
  names(rp_colors) <- paste0(contour_rps, " yr")

  theme_diform <- theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "right",
          panel.grid.minor = element_blank())

  obs_x <- NULL
  if (!is.null(obs_data) && nrow(obs_data)) {
    obs_x <- tryCatch(
      transform_to_diform_space(obs_data, response, diform_result$diform_spec),
      error = function(e) NULL
    )
    if (!is.null(obs_x)) colnames(obs_x) <- comp_names
  }

  make_plot <- function(i, j) {
    cdf <- do.call(rbind, lapply(names(contours), function(rk) {
      v <- contours[[rk]]$vertices_3d
      hb <- hull_boundary_2d(v[, i], v[, j])
      data.frame(x = hb$x, y = hb$y, rp = paste0(rk, " yr"))
    }))
    cdf$rp <- factor(cdf$rp, levels = paste0(contour_rps, " yr"))

    p <- ggplot()
    if (!is.null(obs_x)) {
      p <- p + geom_point(
        data = data.frame(x = obs_x[, i], y = obs_x[, j]),
        aes(x, y), colour = "black", size = 0.8, alpha = 1
      )
    }
    p +
      geom_path(data = cdf, aes(x, y, colour = rp), linewidth = 0.9) +
      scale_colour_manual(values = rp_colors, name = "Return period") +
      labs(
        title = sprintf("D-IFORM: %s vs %s (%s)",
                        comp_labels[j], comp_labels[i],
                        diform_sample_label(diform_result)),
        x = comp_labels[i],
        y = comp_labels[j]
      ) +
      theme_diform
  }

  pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  out <- lapply(pairs, function(idx) make_plot(idx[1], idx[2]))
  names(out) <- vapply(pairs, function(idx) {
    gsub("[^A-Za-z0-9_]+", "_", paste(comp_names[idx], collapse = "_"))
  }, character(1))
  out
}

plot_diform_all <- function(diform_result, obs_data = NULL, response = "hs", resp_label = "Hs [m]") {
  if (identical(diform_result$plot_type, "polar_response")) {
    return(plot_diform_projections(diform_result, obs_data, response, resp_label))
  }
  plot_diform_generic_projections(diform_result, obs_data, response, resp_label)
}
