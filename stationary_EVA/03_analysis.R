library(data.table)
library(evd)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
source(file.path(project_root(), "stationary_EVA", "00_stationary_helpers.R"))
data <- read_rscd_2026()

response_variable <- get0("response_variable", ifnotfound = "hs", inherits = TRUE)
response_info <- stationary_response_info(data, response_variable)
response_values <- response_info$values

chosen_threshold <- as.numeric(quantile(response_values, 0.96))
run_length <- 72
return_periods <- c(1, 2, 5, 10, 20, 50, 100, 200)
obs_per_year <- 8766

#' Runs declustering returning the cluster maxima.
decluster_maxima <- function(x, threshold, run_length) {
  exceedance_index <- which(is.finite(x) & x > threshold)

  if (!length(exceedance_index)) {
    return(data.frame(
      index = integer(),
      value = numeric(),
      excess = numeric()
    ))
  }

  cluster_id <- cumsum(c(TRUE, diff(exceedance_index) > run_length))
  clusters <- split(exceedance_index, cluster_id)
  peak_index <- vapply(
    clusters,
    function(index) index[which.max(x[index])],
    integer(1)
  )

  data.frame(
    index = peak_index,
    value = x[peak_index],
    excess = x[peak_index] - threshold
  )
}

cluster_maxima <- decluster_maxima(
  x = response_values,
  threshold = chosen_threshold,
  run_length = run_length
)
excesses <- cluster_maxima$excess
excesses <- excesses[is.finite(excesses) & excesses > 0]

fit_declust <- fit_stationary_gpd(excesses, threshold = 0)

gpd_scale <- unname(fit_declust$param["scale"])
gpd_shape <- unname(fit_declust$param["shape"])
gpd_has_se <- !is.null(fit_declust$std.err) &&
  all(c("scale", "shape") %in% names(fit_declust$std.err)) &&
  !is.null(fit_declust$var.cov)

gpd_se <- if (gpd_has_se) {
  c(
    scale = unname(fit_declust$std.err["scale"]),
    shape = unname(fit_declust$std.err["shape"])
  )
} else {
  c(scale = NA_real_, shape = NA_real_)
}

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
  m <- lambda * period

  if (abs(shape) < 1e-6) {
    threshold + scale * log(m)
  } else {
    threshold + (scale / shape) * (m^shape - 1)
  }
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

  if (!gpd_has_se || !is.finite(m) || m <= 0) next

  vcov_gpd <- fit_declust$var.cov
  theta_gradient <- if (abs(gpd_shape) < 1e-6) {
    c(log(m), 0)
  } else {
    c(
      (m^gpd_shape - 1) / gpd_shape,
      -gpd_scale * (m^gpd_shape - 1) / gpd_shape^2 +
        gpd_scale * m^gpd_shape * log(m) / gpd_shape
    )
  }

  lambda_gradient <- if (abs(gpd_shape) < 1e-6) {
    gpd_scale * period / m
  } else {
    gpd_scale * period * m^(gpd_shape - 1)
  }

  variance_theta <- as.numeric(
    t(theta_gradient) %*% vcov_gpd %*% theta_gradient
  )
  variance_lambda <- lambda_gradient^2 * var_lambda
  variance_level <- variance_theta + variance_lambda

  if (is.finite(variance_level) && variance_level >= 0) {
    se_level <- sqrt(variance_level)
    rv_declust$lower[i] <- rv_declust$level[i] - 1.96 * se_level
    rv_declust$upper[i] <- rv_declust$level[i] + 1.96 * se_level
  }
}

return_values <- rv_declust

stationary_fitted_response <- response_info$name
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
