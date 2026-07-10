library(data.table)
library(ggplot2)
library(viridisLite)
library(Rcpp)

if (!exists("pot_decluster_all")) {
  sourceCpp(file.path(root, "pot_decluster.cpp"))
}

linear_diform_spec_2d <- function(variables, response, response_component = response) {
  list(
    components = lapply(variables, function(v) list(type = "linear", variable = v, name = v)),
    response = response,
    response_component = response_component,
    plot_type = "linear"
  )
}

diform2d_normalise <- function(X) {
  med <- apply(X, 2, median)
  std <- apply(X, 2, sd)
  std[!is.finite(std) | std <= 0] <- 1
  list(Y = sweep(sweep(X, 2, med), 2, std, "/"), med = med, std = std)
}

diform2d_unnormalise <- function(Y, norm_info) {
  sweep(sweep(Y, 2, norm_info$std, "*"), 2, norm_info$med, "+")
}

generate_direction_vectors_2d <- function(n = 120) {
  a <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  out <- cbind(cos(a), sin(a))
  colnames(out) <- c("d1", "d2")
  out
}

fit_gp_constrained_2d <- function(exc, max_obs = NULL) {
  exc <- exc[exc > 0]
  n <- length(exc)
  if (n < 3) return(list(sigma = mean(exc), xi = -0.1))

  m <- mean(exc)
  xi0 <- max(-0.99, min(-0.01, 0.5 * (1 - m^2 / var(exc))))
  sig0 <- max(0.01, m * (1 - xi0))
  if (!is.null(max_obs)) sig0 <- max(sig0, -xi0 * max_obs * 1.001)

  nll <- function(par) {
    sig <- par[1]
    xi <- par[2]
    z <- 1 + xi * exc / sig
    if (sig <= 0 || xi >= 0 || xi < -1 || any(z <= 0) ||
        (!is.null(max_obs) && max_obs > -sig / xi)) {
      1e12
    } else {
      n * log(sig) + (1 + 1 / xi) * sum(log(z))
    }
  }

  o <- optim(c(sig0, xi0), nll, method = "Nelder-Mead", control = list(maxit = 2000))
  if (o$convergence) list(sigma = sig0, xi = xi0) else list(sigma = o$par[1], xi = o$par[2])
}

gp_return_value_2d <- function(sigma, xi, lambda, T_years) {
  if (abs(xi) < 1e-8) sigma * log(lambda * T_years) else sigma / xi * ((lambda * T_years)^xi - 1)
}

compute_return_values_pot_2d <- function(Y, dirs, return_periods, n_years,
                                         h_hours = 72, zeta = 0.1, n_threads = 0) {
  dec <- pot_decluster_all(cbind(Y, 0), cbind(dirs, 0), h_hours, zeta, n_threads)
  rv <- matrix(NA_real_, nrow(dirs), length(return_periods),
               dimnames = list(NULL, paste0("T", return_periods)))
  fits <- data.frame(mu = dec$mu, sigma = NA_real_, xi = NA_real_,
                     n_peaks = dec$n_peaks, n_exc = lengths(dec$exceedances),
                     lambda = lengths(dec$exceedances) / n_years,
                     method = NA_character_)

  for (j in seq_len(nrow(dirs))) {
    peaks <- dec$peaks[[j]]
    exc <- dec$exceedances[[j]]
    if (length(exc) < 5 || fits$lambda[j] < 0.01) {
      fits$method[j] <- "empirical"
      lam <- length(peaks) / n_years
      p <- pmin(1, pmax(0, 1 - 1 / (lam * return_periods)))
      rv[j, ] <- as.numeric(quantile(peaks, p, type = 7))
    } else {
      gp <- fit_gp_constrained_2d(exc, max(peaks) - fits$mu[j])
      fits$sigma[j] <- gp$sigma
      fits$xi[j] <- gp$xi
      fits$method[j] <- "gpd"
      rv[j, ] <- fits$mu[j] + gp_return_value_2d(gp$sigma, gp$xi, fits$lambda[j], return_periods)
    }
  }

  list(rv = rv, fits = fits, peaks = dec$peaks, exceedances = dec$exceedances)
}

halfspace_polygon_2d <- function(dirs, R, tol = 1e-8) {
  pts <- list()
  k <- 0L
  n <- nrow(dirs)
  scale <- pmax(1, abs(R))

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      A <- rbind(dirs[i, ], dirs[j, ])
      if (abs(det(A)) < 1e-10) next
      x <- as.numeric(solve(A, c(R[i], R[j])))
      if (all(as.numeric(dirs %*% x) <= R + tol * scale)) {
        k <- k + 1L
        pts[[k]] <- x
      }
    }
  }

  if (!length(pts)) return(matrix(numeric(), ncol = 2))
  P <- unique(round(do.call(rbind, pts), 10))
  if (nrow(P) < 3) return(matrix(numeric(), ncol = 2))
  ctr <- colMeans(P)
  P <- P[order(atan2(P[, 2] - ctr[2], P[, 1] - ctr[1])), , drop = FALSE]
  colnames(P) <- c("d1", "d2")
  P
}

build_diform_contours_2d <- function(dirs, rv_matrix, norm_info, return_periods,
                                     variables = c("tp", "hs")) {
  out <- list()

  for (k in seq_along(return_periods)) {
    P_norm <- halfspace_polygon_2d(dirs, rv_matrix[, k])
    if (nrow(P_norm) >= 3) {
      P <- diform2d_unnormalise(P_norm, norm_info)
      P <- pmax(P, 0)
      colnames(P) <- variables
      out[[as.character(return_periods[k])]] <- list(
        vertices_norm = P_norm,
        vertices_2d = P,
        rp = return_periods[k]
      )
    }
  }

  out
}

diform_contour_summary_2d <- function(contours_2d, response = "hs", covariate = "tp") {
  rbindlist(lapply(contours_2d, function(cc) {
    v <- as.data.frame(cc$vertices_2d)
    i <- which.max(v[[response]])
    data.table(
      return_period = cc$rp,
      n_vertices = nrow(v),
      max_resp = round(v[[response]][i], 3),
      covariate_at_max = round(v[[covariate]][i], 3)
    )
  }))
}

run_diform_contours_2d <- function(sim, response = "hs", covariate = "tp", n_years,
                                   return_periods = c(1, 2, 5, 10, 20, 50, 100),
                                   n_directions = 120, h_hours = 72, zeta = 0.1,
                                   n_threads = 0) {
  vars <- c(covariate, response)
  diform_spec <- linear_diform_spec_2d(vars, response, response)
  x <- as.data.frame(sim)[, vars, drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  X <- as.matrix(x)
  storage.mode(X) <- "double"

  norm <- diform2d_normalise(X)
  dirs <- generate_direction_vectors_2d(n_directions)
  pot <- compute_return_values_pot_2d(norm$Y, dirs, return_periods, n_years,
                                      h_hours, zeta, n_threads)
  contours <- build_diform_contours_2d(dirs, pot$rv, norm, return_periods, vars)

  list(
    contours_2d = contours,
    directions = dirs,
    rv_matrix = pot$rv,
    fits = pot$fits,
    peaks = pot$peaks,
    exceedances = pot$exceedances,
    norm_info = norm,
    X = X,
    X_transformed = X,
    response = response,
    covariate = covariate,
    variables = vars,
    component_names = vars,
    diform_spec = diform_spec,
    plot_type = diform_spec$plot_type,
    n_years = n_years,
    return_periods = return_periods,
    summary = diform_contour_summary_2d(contours, response, covariate)
  )
}

diform2d_sample_label <- function(diform_result) {
  n_years <- diform_result$n_years
  n_years <- if (is.finite(n_years) && abs(n_years - round(n_years)) < 1e-6) {
    round(n_years)
  } else {
    round(n_years, 2)
  }
  source <- if (isFALSE(diform_result$bootstrap)) "observed data" else "simulation"
  sprintf("%s yr %s", format(n_years, big.mark = ",", scientific = FALSE), source)
}

diform2d_contour_df <- function(diform_result, rp = "all") {
  rps <- as.numeric(names(diform_result$contours_2d))
  rps <- rps[is.finite(rps)]
  if (!identical(rp, "all")) rps <- intersect(rps, as.numeric(rp))

  rbindlist(lapply(rps, function(x) {
    v <- as.data.table(diform_result$contours_2d[[as.character(x)]]$vertices_2d)
    if (nrow(v)) v <- rbind(v, v[1])
    v[, return_period := x]
    v
  }), fill = TRUE)
}

plot_diform_2d_all <- function(diform_result, obs_data = NULL,
                               covariate_label = "Tp [s]",
                               response_label = "Hs [m]") {
  covariate <- diform_result$covariate
  response <- diform_result$response
  x_var <- response
  y_var <- covariate
  cdf <- diform2d_contour_df(diform_result, "all")
  cdf[, rp_label := paste0(return_period, " yr")]
  cdf[, rp_label := factor(rp_label, levels = paste0(diform_result$return_periods, " yr"))]
  rp_vals <- sort(unique(cdf$return_period))
  rp_colors <- viridisLite::turbo(length(rp_vals))
  names(rp_colors) <- paste0(rp_vals, " yr")

  p <- ggplot()
  if (!is.null(obs_data) && nrow(obs_data)) {
    p <- p + geom_point(
      data = data.frame(x = obs_data[[x_var]], y = obs_data[[y_var]]),
      aes(x, y), colour = "black", alpha = 0.55, size = 0.65
    )
  }

  p +
    geom_path(data = cdf, aes(.data[[x_var]], .data[[y_var]], colour = rp_label),
              linewidth = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = rp_colors, name = "Return period") +
    labs(
      title = sprintf("2D D-IFORM: %s vs %s (%s)",
                      response, covariate, diform2d_sample_label(diform_result)),
      x = response_label,
      y = covariate_label
    ) +
    coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    theme_minimal(base_size = 17) +
    theme(panel.grid.minor = element_blank(),
           plot.title = element_text(face = "bold", size = 21),
           axis.title = element_text(size = 18),
           legend.title = element_text(size = 16),
           legend.text = element_text(size = 15))
}

save_diform_2d_all_r <- function(diform_result, obs_data = NULL, out_dir = ".",
                                  covariate_label = "Tp [s]",
                                  response_label = "Hs [m]",
                                  width = 10, height = 9) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "diform_2d_all_return_periods.png")
  p_contour <- plot_diform_2d_all(diform_result, obs_data, covariate_label, response_label)
  grDevices::png(path, width = round(width * 180), height = round(height * 180),
                 res = 180, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p_contour)
  invisible(path)
}

diform2d_compare_layout <- function(n) {
  if (n <= 3) return(c(1, n))
  nc <- ceiling(sqrt(n))
  nr <- ceiling(n / nc)
  c(nr, nc)
}

diform2d_compare_title <- function(diform_result, name) {
  label <- if (!is.null(diform_result$scenario_label)) diform_result$scenario_label else name
  sprintf("%s\n%s", label, diform2d_sample_label(diform_result))
}

diform2d_limits <- function(vertices_list, obs = NULL, x_var = "hs", y_var = "tp") {
  all_vertices <- rbindlist(lapply(vertices_list, as.data.table), fill = TRUE)
  xv <- all_vertices[[x_var]]
  yv <- all_vertices[[y_var]]
  if (!is.null(obs) && nrow(obs)) {
    xv <- c(xv, obs[[x_var]])
    yv <- c(yv, obs[[y_var]])
  }

  pad <- function(x, lower_zero = FALSE) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(0, 1))
    r <- range(x)
    span <- diff(r)
    if (!is.finite(span) || span <= 0) span <- max(1, abs(r[1]) * 0.1)
    r <- r + c(-1, 1) * span * 0.06
    if (lower_zero) r[1] <- max(0, r[1])
    r
  }

  list(x = pad(xv, lower_zero = TRUE), y = pad(yv, lower_zero = TRUE))
}

subsample_obs_2d <- function(obs, max_points = Inf) {
  if (is.null(obs) || !nrow(obs) || nrow(obs) <= max_points) return(obs)
  obs[unique(round(seq(1, nrow(obs), length.out = max_points))), , drop = FALSE]
}

plot_diform_2d_compare_r <- function(diform_results, obs_data = NULL, out_dir = ".",
                                     rp, covariate_label = "Tp [s]",
                                     response_label = "Hs [m]",
                                     max_background_points = Inf,
                                     file_prefix = "diform_2d_compare") {
  if (is.null(names(diform_results))) names(diform_results) <- rep("", length(diform_results))
  missing_names <- names(diform_results) == ""
  names(diform_results)[missing_names] <- paste0("scenario_", which(missing_names))

  has_rp <- vapply(diform_results, function(x) as.character(rp) %in% names(x$contours_2d), logical(1))
  diform_results <- diform_results[has_rp]

  covariate <- diform_results[[1]]$covariate
  response <- diform_results[[1]]$response
  x_var <- response
  y_var <- covariate
  vertices_list <- lapply(diform_results, function(x) {
    as.data.frame(x$contours_2d[[as.character(rp)]]$vertices_2d)
  })

  obs <- NULL
  if (!is.null(obs_data) && nrow(obs_data)) {
    obs <- as.data.frame(obs_data)[, c(covariate, response), drop = FALSE]
    obs <- obs[complete.cases(obs), , drop = FALSE]
    obs <- subsample_obs_2d(obs, max_background_points)
  }

  lim <- diform2d_limits(vertices_list, obs, x_var, y_var)
  cols <- grDevices::hcl.colors(length(diform_results), "Dark 3")
  is_raw <- vapply(seq_along(diform_results), function(i) {
    identical(diform_results[[i]]$model, "raw") ||
      grepl("raw|no GAM", paste(names(diform_results)[i], diform_results[[i]]$scenario_label), ignore.case = TRUE)
  }, logical(1))
  line_order <- c(which(!is_raw), which(is_raw))
  ltys <- ifelse(is_raw, 2, 1)
  lwds <- ifelse(is_raw, 3.0, 2.2)
  pchs <- ifelse(is_raw, 17, 18)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  path <- file.path(out_dir, sprintf("%s_%syr.png", file_prefix, format(rp, trim = TRUE)))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1650, height = 1250, res = 160, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(mar = c(6.2, 6.4, 5.6, 1.4), xpd = NA)
  on.exit(par(old), add = TRUE)

  plot(NA, NA, xlim = lim$x, ylim = lim$y, xlab = response_label,
       ylab = covariate_label,
       main = sprintf("2D D-IFORM %s-yr contour comparison", format(rp, trim = TRUE)),
       cex.lab = 1.70, cex.axis = 1.35, cex.main = 1.72)

  legend_labels <- character(length(diform_results))
  for (i in seq_along(diform_results)) {
    P <- vertices_list[[i]]
    polygon(P[[x_var]], P[[y_var]],
            col = grDevices::adjustcolor(cols[i], 0.06),
            border = NA)
    legend_labels[i] <- if (!is.null(diform_results[[i]]$scenario_label)) {
      diform_results[[i]]$scenario_label
    } else {
      names(diform_results)[i]
    }
  }

  if (!is.null(obs) && nrow(obs)) {
    points(obs[[x_var]], obs[[y_var]], pch = 16,
           cex = 0.38, col = grDevices::adjustcolor("black", 0.45))
  }

  for (i in line_order) {
    P <- vertices_list[[i]]
    x <- c(P[[x_var]], P[[x_var]][1])
    y <- c(P[[y_var]], P[[y_var]][1])
    lines(x, y, col = cols[i], lwd = lwds[i], lty = ltys[i])
    j <- which.max(P[[response]])
    points(P[[x_var]][j], P[[y_var]][j],
           pch = pchs[i], cex = 1.2, col = cols[i], bg = cols[i])
  }

  legend("topright",
         legend = legend_labels,
         col = cols,
         lwd = lwds,
         lty = ltys,
         pch = pchs,
         pt.cex = 1.20,
         bg = grDevices::adjustcolor("white", 0.92),
         box.col = "gray80",
         cex = 1.20,
         y.intersp = 1.12)

  invisible(path)
}
