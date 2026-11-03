library(ggplot2)

#' Quantile function of the generalized Pareto distribution.
diform_gp_quantile <- function(p, sigma, xi) {
  p <- pmin(pmax(p, 0), 1 - .Machine$double.eps)
  if (abs(xi) < 1e-8) {
    -sigma * log1p(-p)
  } else {
    sigma / xi * ((1 - p)^(-xi) - 1)
  }
}

#' Wald confidence intervals for the GPD scale (log scale) and shape from the observed information.
diform_gp_parameter_ci <- function(excesses, sigma, xi, level = 0.95) {
  excesses <- excesses[is.finite(excesses) & excesses > 0]
  n <- length(excesses)
  missing_ci <- c(
    sigma_lower = NA_real_, sigma_upper = NA_real_,
    xi_lower = NA_real_, xi_upper = NA_real_
  )
  if (n < 3 || !is.finite(sigma) || sigma <= 0 || !is.finite(xi)) return(missing_ci)

  #' Negative log-likelihood of the generalized Pareto distribution (log-scale, shape).
  nll <- function(par) {
    scale <- exp(par[1])
    shape <- par[2]
    if (!is.finite(scale) || !is.finite(shape)) return(1e12)
    if (abs(shape) < 1e-7) return(n * log(scale) + sum(excesses) / scale)
    support <- 1 + shape * excesses / scale
    if (any(support <= 0)) return(1e12)
    n * log(scale) + (1 + 1 / shape) * sum(log(support))
  }

  hessian <- tryCatch(
    stats::optimHess(c(log(sigma), xi), nll),
    error = function(e) NULL
  )
  if (is.null(hessian) || any(!is.finite(hessian))) return(missing_ci)
  hessian <- (hessian + t(hessian)) / 2
  if (any(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values <= 0)) return(missing_ci)

  covariance <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(covariance) || any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    return(missing_ci)
  }

  critical <- stats::qnorm((1 + level) / 2)
  se_log_sigma <- sqrt(covariance[1, 1])
  se_xi <- sqrt(covariance[2, 2])
  c(
    sigma_lower = exp(log(sigma) - critical * se_log_sigma),
    sigma_upper = exp(log(sigma) + critical * se_log_sigma),
    xi_lower = xi - critical * se_xi,
    xi_upper = xi + critical * se_xi
  )
}

#' Short text label of a direction vector.
diform_direction_label <- function(u, digits = 2) {
  u[abs(u) < 0.5 * 10^(-digits)] <- 0
  values <- formatC(u, format = "f", digits = digits)
  values <- sub("\\.?0+$", "", values)
  paste0("u = (", paste(values, collapse = ", "), ")")
}

#' Coordinate directions plus a few random oblique directions for the QQ diagnostics.
select_diform_qq_directions <- function(diform_result, n_random = 3L, seed = 2023L) {
  dirs <- as.matrix(diform_result$directions)
  fits <- diform_result$fits
  valid <- which(
    fits$method == "gpd" &
      is.finite(fits$sigma) & fits$sigma > 0 &
      is.finite(fits$xi) &
      lengths(diform_result$exceedances) >= 5
  )

  axes <- diag(ncol(dirs))
  axis_index <- vapply(seq_len(nrow(axes)), function(i) {
    candidate <- valid[which.max(drop(dirs[valid, , drop = FALSE] %*% axes[i, ]))]
    if (sum(dirs[candidate, ] * axes[i, ]) < 1 - 1e-8) {
      stop("No fitted GPD is available for axis direction ", i, ".", call. = FALSE)
    }
    candidate
  }, integer(1))

  random_pool <- setdiff(valid, axis_index)
  random_pool <- random_pool[apply(abs(dirs[random_pool, , drop = FALSE]), 1, max) < 0.95]
  if (length(random_pool) < n_random) {
    stop("Not enough oblique directions with a fitted GPD for the diagnostic plot.", call. = FALSE)
  }

  set.seed(seed)
  random_index <- sample(random_pool, n_random)

  data.frame(
    index = c(axis_index, random_index),
    order = seq_len(length(axis_index) + n_random),
    stringsAsFactors = FALSE
  )
}

#' QQ points, pointwise envelopes and parameter summaries for the selected directions.
build_diform_gpd_qq_data <- function(diform_result, selection, level = 0.95) {
  dirs <- as.matrix(diform_result$directions)
  alpha <- 1 - level

  qq <- lapply(seq_len(nrow(selection)), function(k) {
    j <- selection$index[k]
    fit <- diform_result$fits[j, ]
    observed <- sort(as.numeric(diform_result$exceedances[[j]]))
    n <- length(observed)
    rank <- seq_len(n)
    probability <- (rank - 0.5) / n
    panel <- paste0("(", letters[k], ") ", diform_direction_label(dirs[j, ]))

    data.frame(
      panel = panel,
      theoretical = diform_gp_quantile(probability, fit$sigma, fit$xi),
      observed = observed,
      lower = diform_gp_quantile(qbeta(alpha / 2, rank, n - rank + 1), fit$sigma, fit$xi),
      upper = diform_gp_quantile(qbeta(1 - alpha / 2, rank, n - rank + 1), fit$sigma, fit$xi),
      stringsAsFactors = FALSE
    )
  })

  qq <- do.call(rbind, qq)
  panel_levels <- vapply(seq_len(nrow(selection)), function(k) {
    j <- selection$index[k]
    paste0("(", letters[k], ") ", diform_direction_label(dirs[j, ]))
  }, character(1))
  qq$panel <- factor(qq$panel, levels = panel_levels)

  stats <- lapply(seq_len(nrow(selection)), function(k) {
    j <- selection$index[k]
    fit <- diform_result$fits[j, ]
    intervals <- diform_gp_parameter_ci(
      diform_result$exceedances[[j]], fit$sigma, fit$xi, level
    )
    sigma_line <- if (all(is.finite(intervals[c("sigma_lower", "sigma_upper")]))) {
      sprintf("\u03c3 = %.2f [%.2f ; %.2f]", fit$sigma,
              intervals["sigma_lower"], intervals["sigma_upper"])
    } else {
      sprintf("\u03c3 = %.2f", fit$sigma)
    }
    xi_line <- if (all(is.finite(intervals[c("xi_lower", "xi_upper")]))) {
      sprintf("\u03be = %.2f [%.2f ; %.2f]", fit$xi,
              intervals["xi_lower"], intervals["xi_upper"])
    } else {
      sprintf("\u03be = %.2f", fit$xi)
    }
    data.frame(
      panel = factor(panel_levels[k], levels = panel_levels),
      label = paste(sprintf("n = %d", fit$n_exc), sigma_line, xi_line, sep = "\n"),
      stringsAsFactors = FALSE
    )
  })

  list(points = qq, stats = do.call(rbind, stats))
}

#' Faceted QQ plots of the directional GPD fits, optionally saved to file.
plot_diform_gpd_qq <- function(diform_result,
                               n_random = 3L,
                               seed = 2023L,
                               level = 0.95,
                               out_file = NULL,
                               width = 13.5,
                               height = 8.6,
                               dpi = 300) {
  selection <- select_diform_qq_directions(diform_result, n_random, seed)
  diagnostic <- build_diform_gpd_qq_data(diform_result, selection, level)
  qq <- diagnostic$points
  panel_stats <- diagnostic$stats

  limit <- max(qq$observed, qq$theoretical, qq$upper, na.rm = TRUE) * 1.04
  panel_stats$x <- 0.035 * limit
  panel_stats$y <- 0.965 * limit

  p <- ggplot(qq, aes(x = theoretical, y = observed)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#8ECAE6",
      alpha = 0.34,
      linewidth = 0
    ) +
    geom_abline(intercept = 0, slope = 1, color = "#D1495B", linewidth = 0.8) +
    geom_point(color = "#173F5F", fill = "white", shape = 21,
               size = 1.9, stroke = 0.55, alpha = 0.88) +
    geom_text(
      data = panel_stats,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 4.8,
      lineheight = 1.05,
      color = "#59636B"
    ) +
    facet_wrap(~panel, ncol = 3, drop = FALSE) +
    coord_cartesian(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    labs(
      title = "Directional GPD fit diagnostics",
      x = "Theoretical GPD quantiles (normalised projection)",
      y = "Observed excesses (normalised projection)"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold", size = 25, color = "#17242C"),
      axis.title = element_text(size = 17, color = "#263238"),
      axis.text = element_text(size = 14, color = "#53616A"),
      strip.text = element_text(face = "bold", size = 14.5, color = "#21313A",
                                lineheight = 1.15, margin = margin(10, 5, 10, 5)),
      strip.background = element_rect(fill = "#F1F5F7", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E3E9EC", linewidth = 0.45),
      panel.spacing = grid::unit(1.05, "lines"),
      plot.margin = margin(14, 18, 12, 14),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved: ", out_file)
  }

  attr(p, "direction_selection") <- selection
  p
}
