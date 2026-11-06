run_length <- 72
return_periods <- c(1, 2, 5, 10, 20, 50, 100, 200)
obs_per_year <- 8766

cluster_maxima <- decluster_maxima(
  x = response_values,
  threshold = chosen_threshold,
  run_length = run_length
)
excesses <- cluster_maxima$excess

fit_declust <- fit_stationary_gpd(excesses, threshold = 0)

gpd_scale <- unname(fit_declust$param["scale"])
gpd_shape <- unname(fit_declust$param["shape"])
gpd_se <- c(scale = unname(fit_declust$std.err["scale"]), shape = unname(fit_declust$std.err["shape"]))

fit_summary <- data.frame(
  parameter = c("scale", "shape"),
  estimate = c(gpd_scale, gpd_shape),
  se = unname(gpd_se),
  lower = c(gpd_scale, gpd_shape) - 1.96 * unname(gpd_se),
  upper = c(gpd_scale, gpd_shape) + 1.96 * unname(gpd_se)
)

n_years <- length(response_values) / obs_per_year
n_clusters <- nrow(cluster_maxima)
lambda_declust <- n_clusters / n_years
var_lambda <- lambda_declust / n_years

rv_declust <- data.frame(
  period = return_periods,
  level = NA_real_,
  lower = NA_real_,
  upper = NA_real_,
  method = "GPD declustered"
)

#' POT return level for a return period, exceedance rate, scale and shape.
gpd_return_level <- function(period, threshold, lambda, scale, shape) {
  threshold + (scale / shape) * ((lambda * period)^shape - 1)
}

for (i in seq_along(return_periods)) {
  period <- return_periods[i]
  m <- lambda_declust * period

  rv_declust$level[i] <- gpd_return_level(
    period = period,
    threshold = chosen_threshold,
    lambda = lambda_declust,
    scale = gpd_scale,
    shape = gpd_shape
  )

  vcov_gpd <- fit_declust$var.cov
  theta_gradient <- c(
    (m^gpd_shape - 1) / gpd_shape,
    -gpd_scale * (m^gpd_shape - 1) / gpd_shape^2 + gpd_scale * m^gpd_shape * log(m) / gpd_shape
  )
  lambda_gradient <- gpd_scale * period * m^(gpd_shape - 1)

  variance_theta <- as.numeric(
    t(theta_gradient) %*% vcov_gpd %*% theta_gradient
  )
  variance_lambda <- lambda_gradient^2 * var_lambda
  variance_level <- variance_theta + variance_lambda

  se_level <- sqrt(variance_level)
  rv_declust$lower[i] <- rv_declust$level[i] - 1.96 * se_level
  rv_declust$upper[i] <- rv_declust$level[i] + 1.96 * se_level
}

return_values <- rv_declust

stationary_analysis <- list(
  response = response_info[c("name", "tag", "label", "unit")],
  threshold = chosen_threshold,
  run_length = run_length,
  n_years = n_years,
  n_clusters = n_clusters,
  cluster_rate = lambda_declust,
  cluster_maxima = cluster_maxima,
  excesses = excesses,
  fit = fit_declust,
  fit_summary = fit_summary,
  return_values = return_values
)
