library(data.table)
library(Rcpp)

sourceCpp("contours/pot_decluster.cpp")

#' Working-space coordinates (Cspd, H_L, H_T): current speed and the square root of the response along the wave-current angle.
transform_to_diform_space <- function(data, response) {
  a <- data$angleDiff * pi / 180
  r <- sqrt(data[[response]])
  cbind(Cspd = data$Cspd, H_L = r * cos(a), H_T = r * sin(a))
}

#' Centre each working coordinate by its median and scale it by its standard deviation.
normalise_diform <- function(X) {
  med <- apply(X, 2, median)
  std <- apply(X, 2, sd)
  list(Y = sweep(sweep(X, 2, med), 2, std, "/"), med = med, std = std)
}

#' Map normalized points back to working-space coordinates.
unnormalise_diform <- function(Y, norm_info) sweep(sweep(Y, 2, norm_info$std, "*"), 2, norm_info$med, "+")

#' Map working-space points back to current speed, wave-current angle and response.
backtransform_to_original <- function(pts, response) {
  out <- data.frame(Cspd = pts[, 1], angleDiff = (atan2(pts[, 3], pts[, 2]) * 180 / pi) %% 360,
                    r = sqrt(pts[, 2]^2 + pts[, 3]^2)^2)
  names(out)[3] <- response
  out
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

#' Fit a generalized Pareto distribution to excesses by maximum likelihood (exponential start values if the optimizer does not converge).
fit_gp <- function(exc) {
  n <- length(exc)
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
  start <- c(log(mean(exc)), 0)
  o <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 2000))
  if (o$convergence) list(sigma = exp(start[1]), xi = start[2]) else list(sigma = exp(o$par[1]), xi = o$par[2])
}

#' Return level above the threshold for a GPD with yearly exceedance rate lambda.
gp_return_value <- function(sigma, xi, lambda, T_years) {
  if (abs(xi) < 1e-8) sigma * log(lambda * T_years) else sigma / xi * ((lambda * T_years)^xi - 1)
}

#' Decluster each directional projection and estimate its return values with a GPD.
compute_return_values_pot <- function(Y, dirs, return_periods, n_years, h_hours = 72, zeta = 0.1, n_threads = 0) {
  dec <- pot_decluster_all(Y, dirs, h_hours, zeta, n_threads)
  rv <- matrix(NA_real_, nrow(dirs), length(return_periods), dimnames = list(NULL, paste0("T", return_periods)))
  fits <- data.frame(mu = dec$mu, sigma = NA_real_, xi = NA_real_, n_peaks = dec$n_peaks,
                     n_exc = lengths(dec$exceedances), lambda = lengths(dec$exceedances) / n_years)
  for (j in seq_len(nrow(dirs))) {
    gp <- fit_gp(dec$exceedances[[j]])
    fits$sigma[j] <- gp$sigma
    fits$xi[j] <- gp$xi
    rv[j, ] <- fits$mu[j] + sapply(return_periods, gp_return_value, sigma = gp$sigma, xi = gp$xi, lambda = fits$lambda[j])
  }
  list(rv = rv, fits = fits, exceedances = dec$exceedances)
}

#' Contour vertices for each return period from the half-space intersection along every direction.
build_diform_contours_3d <- function(dirs, rv_matrix, norm_info, return_periods, response) {
  dots <- dirs %*% t(dirs)
  out <- list()
  for (k in seq_along(return_periods)) {
    rv <- rv_matrix[, k]
    pts <- t(vapply(seq_len(nrow(dirs)), function(j) {
      pos <- dots[, j] > 1e-10
      min(rv[pos] / dots[pos, j]) * dirs[j, ]
    }, numeric(3)))
    orig <- unnormalise_diform(pts, norm_info)
    colnames(orig) <- c("Cspd", "H_L", "H_T")
    out[[as.character(return_periods[k])]] <- list(vertices_3d = orig,
                                                   vertices_bt = backtransform_to_original(orig, response),
                                                   rp = return_periods[k])
  }
  out
}

#' Maximum response on each contour and the covariate values where it occurs.
diform_contour_summary <- function(contours_3d, response) {
  rbindlist(lapply(contours_3d, function(cc) {
    bt <- cc$vertices_bt
    i <- which.max(bt[[response]])
    out <- data.table(return_period = cc$rp, n_vertices = nrow(bt), max_resp = round(max(bt[[response]]), 3))
    for (nm in names(bt)) out[[paste0(nm, "_at_max")]] <- round(bt[[nm]][i], 3)
    out
  }))
}

#' Run the D-IFORM analysis: transformation, normalization, directional return values and contours.
run_diform_contours <- function(data, response, n_years, return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                simplex_resolution = 20, h_hours = 72, zeta = 0.05, n_threads = 0) {
  norm <- normalise_diform(transform_to_diform_space(data, response))
  dirs <- generate_direction_vectors_3d(simplex_resolution)
  pot <- compute_return_values_pot(norm$Y, dirs, return_periods, n_years, h_hours, zeta, n_threads)
  contours <- build_diform_contours_3d(dirs, pot$rv, norm, return_periods, response)
  list(contours_3d = contours, directions = dirs, rv_matrix = pot$rv, fits = pot$fits,
       exceedances = pot$exceedances, norm_info = norm, summary = diform_contour_summary(contours, response),
       response = response, n_years = n_years, return_periods = return_periods)
}
