library(data.table)
library(ggplot2)
library(viridisLite)
library(Rcpp)


sourceCpp(file.path(root, "pot_decluster.cpp"))

#' Default specification: current speed plus the response and wave-current angle as two Cartesian components.
default_diform_spec <- function(response = "hs") {
  list(
    components = list(
      list(type = "linear", variable = "Cspd", name = "Cspd"),
      list(type = "circular_radius", radius = response, angle = "angleDiff",
           transform = "sqrt", names = c("H_L", "H_T"),
           angle_unit = "degrees", period = 360, origin = 0)
    ),
    response = response,
    response_component = NULL,
    plot_type = "polar_response"
  )
}

#' Specification combining a linear variable with a radius-angle pair mapped to Cartesian components.
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

#' Specification in which every variable is used directly as a Cartesian coordinate.
linear_diform_spec <- function(variables, response = "hs", response_component = response) {
  list(
    components = lapply(variables, function(v) list(type = "linear", variable = v, name = v)),
    response = response,
    response_component = response_component,
    plot_type = "linear"
  )
}

#' Return a complete specification, building the default one when none is given.
normalise_diform_spec <- function(diform_spec = NULL, response = "hs") {
  if (is.null(diform_spec)) diform_spec <- default_diform_spec(response)
  if (is.character(diform_spec)) diform_spec <- linear_diform_spec(diform_spec, response)
  if (is.null(diform_spec$response)) diform_spec$response <- response
  diform_spec
}

#' Names of the working-space coordinates defined by a specification.
diform_component_names <- function(diform_spec) {
  unlist(lapply(diform_spec$components, function(comp) {
    if (identical(comp$type, "circular_radius")) return(comp$names)
    if (!is.null(comp$name)) comp$name else comp$variable
  }), use.names = FALSE)
}

#' Apply the radial transform (identity, square root or log).
transform_radius <- function(x, transform = "identity") {
  if (identical(transform, "sqrt")) return(sqrt(pmax(x, 0)))
  if (identical(transform, "log")) return(log(x))
  x
}

#' Invert the radial transform.
inverse_radius <- function(x, transform = "identity") {
  if (identical(transform, "sqrt")) return(x^2)
  if (identical(transform, "log")) return(exp(x))
  x
}

#' Convert angles to radians.
angle_to_radians <- function(x, unit = "degrees") {
  if (identical(unit, "radians")) return(x)
  x * pi / 180
}

#' Convert radians to the requested angular unit and wrap them to one period.
angle_from_radians <- function(x, unit = "degrees", period = NULL, origin = NULL) {
  out <- if (identical(unit, "radians")) x else x * 180 / pi
  if (is.null(period)) period <- if (identical(unit, "radians")) 2 * pi else 360
  if (is.null(origin)) origin <- 0
  ((out - origin) %% period) + origin
}

#' Map observations to the working space defined by the specification.
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

#' Centre each working coordinate by its median and scale it by its standard deviation.
normalise_diform <- function(X) {
  med <- apply(X, 2, median)
  std <- apply(X, 2, sd)
  list(Y = sweep(sweep(X, 2, med), 2, std, "/"), med = med, std = std)
}

#' Map normalized points back to working-space coordinates.
unnormalise_diform <- function(Y, norm_info) sweep(sweep(Y, 2, norm_info$std, "*"), 2, norm_info$med, "+")

#' Map working-space points back to the original physical variables.
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

#' Unit direction vectors on the sphere from a reflected simplex grid.
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

#' Fit a generalized Pareto distribution to positive excesses by maximum likelihood.
fit_gp <- function(exc) {
  exc <- exc[exc > 0]
  n <- length(exc)
  if (n < 3) return(list(sigma = mean(exc), xi = 0))
  m <- mean(exc)
  #' Negative log-likelihood of the generalized Pareto distribution (log-scale, shape).
  nll <- function(par) {
    sig <- exp(par[1])
    xi <- par[2]
    if (!is.finite(sig) || !is.finite(xi)) return(1e12)
    if (abs(xi) < 1e-7) return(n * log(sig) + sum(exc) / sig)
    z <- 1 + xi * exc / sig
    if (any(z <= 0)) 1e12
    else n * log(sig) + (1 + 1 / xi) * sum(log(z))
  }
  start <- c(log(max(m, .Machine$double.eps)), 0)
  o <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 2000))
  if (o$convergence) list(sigma = exp(start[1]), xi = start[2]) else list(sigma = exp(o$par[1]), xi = o$par[2])
}

#' Return level above the threshold for a GPD with yearly exceedance rate lambda.
gp_return_value <- function(sigma, xi, lambda, T_years) {
  if (abs(xi) < 1e-8) sigma * log(lambda * T_years) else sigma / xi * ((lambda * T_years)^xi - 1)
}

#' Decluster each directional projection and estimate its return values with a GPD (empirical quantiles when exceedances are too few).
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
      gp <- fit_gp(exc)
      fits$sigma[j] <- gp$sigma
      fits$xi[j] <- gp$xi
      fits$method[j] <- "gpd"
      rv[j, ] <- fits$mu[j] + sapply(return_periods, gp_return_value, sigma = gp$sigma, xi = gp$xi, lambda = fits$lambda[j])
    }
  }
  list(rv = rv, fits = fits, peaks = dec$peaks, exceedances = dec$exceedances)
}

#' Contour vertices for each return period from the half-space intersection along every direction.
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

#' Run the D-IFORM analysis: transformation, normalization, directional return values and contours.
run_diform_contours <- function(data, response = "hs", n_years,
                                return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                simplex_resolution = 20, h_hours = 72, zeta = 0.1, n_threads = 0,
                                diform_spec = NULL) {
  diform_spec <- normalise_diform_spec(diform_spec, response)
  X <- transform_to_diform_space(data, response, diform_spec)
  norm <- normalise_diform(X)
  dirs <- generate_direction_vectors_3d(simplex_resolution)
  pot <- compute_return_values_pot(norm$Y, dirs, return_periods, n_years, h_hours, zeta, n_threads)
  contours <- build_diform_contours_3d(dirs, pot$rv, norm, return_periods, response, diform_spec)
  list(data = data, contours_3d = contours, directions = dirs, rv_matrix = pot$rv, fits = pot$fits,
       peaks = pot$peaks, exceedances = pot$exceedances, norm_info = norm, X_transformed = X,
       summary = diform_contour_summary(contours, response), response = response,
       n_years = n_years, return_periods = return_periods,
       diform_spec = diform_spec,
       component_names = colnames(X),
       plot_type = diform_spec$plot_type)
}

#' Maximum response on each contour and the covariate values where it occurs.
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

