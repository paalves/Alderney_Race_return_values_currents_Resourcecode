library(ggplot2)

#' Quantile function of the generalized Pareto distribution.
diform_gp_quantile <- function(p, sigma, xi) {
  if (abs(xi) < 1e-8) -sigma * log1p(-p) else sigma / xi * ((1 - p)^(-xi) - 1)
}

#' Wald confidence intervals for the GPD scale (log scale) and shape from the observed information.
diform_gp_parameter_ci <- function(excesses, sigma, xi, level = 0.95) {
  n <- length(excesses)
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
  hessian <- stats::optimHess(c(log(sigma), xi), nll)
  covariance <- solve((hessian + t(hessian)) / 2)
  critical <- stats::qnorm((1 + level) / 2)
  se_log_sigma <- sqrt(covariance[1, 1])
  se_xi <- sqrt(covariance[2, 2])
  c(sigma_lower = exp(log(sigma) - critical * se_log_sigma), sigma_upper = exp(log(sigma) + critical * se_log_sigma),
    xi_lower = xi - critical * se_xi, xi_upper = xi + critical * se_xi)
}

#' Short text label of a direction vector.
diform_direction_label <- function(u, digits = 2) {
  u[abs(u) < 0.5 * 10^(-digits)] <- 0
  values <- sub("[.]?0+$", "", formatC(u, format = "f", digits = digits))
  paste0("u = (", paste(values, collapse = ", "), ")")
}

#' Indices of the three coordinate directions plus n_random random oblique directions.
select_diform_qq_directions <- function(dirs, n_random = 3L, seed = 2023L) {
  axis_index <- vapply(1:3, function(i) which.max(dirs[, i]), integer(1))
  random_pool <- setdiff(seq_len(nrow(dirs)), axis_index)
  random_pool <- random_pool[apply(abs(dirs[random_pool, , drop = FALSE]), 1, max) < 0.95]
  set.seed(seed)
  c(axis_index, sample(random_pool, n_random))
}

#' Faceted QQ plots of the directional GPD fits with 95 % pointwise envelopes, saved to file.
plot_diform_gpd_qq <- function(diform_result, out_file, level = 0.95) {
  dirs <- diform_result$directions
  selection <- select_diform_qq_directions(dirs)
  panel_levels <- paste0("(", letters[seq_along(selection)], ") ",
                         vapply(selection, function(j) diform_direction_label(dirs[j, ]), character(1)))
  alpha <- 1 - level

  qq <- do.call(rbind, lapply(seq_along(selection), function(k) {
    fit <- diform_result$fits[selection[k], ]
    observed <- sort(diform_result$exceedances[[selection[k]]])
    n <- length(observed)
    rank <- seq_len(n)
    data.frame(
      panel = panel_levels[k],
      theoretical = diform_gp_quantile((rank - 0.5) / n, fit$sigma, fit$xi),
      observed = observed,
      lower = diform_gp_quantile(qbeta(alpha / 2, rank, n - rank + 1), fit$sigma, fit$xi),
      upper = diform_gp_quantile(qbeta(1 - alpha / 2, rank, n - rank + 1), fit$sigma, fit$xi)
    )
  }))
  qq$panel <- factor(qq$panel, levels = panel_levels)

  panel_stats <- do.call(rbind, lapply(seq_along(selection), function(k) {
    fit <- diform_result$fits[selection[k], ]
    ci <- diform_gp_parameter_ci(diform_result$exceedances[[selection[k]]], fit$sigma, fit$xi, level)
    data.frame(
      panel = factor(panel_levels[k], levels = panel_levels),
      label = paste(sprintf("n = %d", fit$n_exc),
                    sprintf("σ = %.2f [%.2f ; %.2f]", fit$sigma, ci["sigma_lower"], ci["sigma_upper"]),
                    sprintf("ξ = %.2f [%.2f ; %.2f]", fit$xi, ci["xi_lower"], ci["xi_upper"]), sep = "\n")
    )
  }))
  limit <- max(qq$observed, qq$theoretical, qq$upper) * 1.04
  panel_stats$x <- 0.035 * limit
  panel_stats$y <- 0.965 * limit

  p <- ggplot(qq, aes(x = theoretical, y = observed)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#8ECAE6", alpha = 0.34, linewidth = 0) +
    geom_abline(intercept = 0, slope = 1, color = "#D1495B", linewidth = 0.8) +
    geom_point(color = "#173F5F", fill = "white", shape = 21, size = 1.9, stroke = 0.55, alpha = 0.88) +
    geom_text(data = panel_stats, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 4.8, lineheight = 1.05, color = "#59636B") +
    facet_wrap(~panel, ncol = 3, drop = FALSE) +
    coord_cartesian(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    labs(title = "Directional GPD fit diagnostics",
         x = "Theoretical GPD quantiles (normalised projection)",
         y = "Observed excesses (normalised projection)") +
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
  ggsave(out_file, p, width = 13.5, height = 8.6, dpi = 300, bg = "white")
}
