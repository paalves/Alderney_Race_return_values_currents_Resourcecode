library(data.table)

predict_threshold <- function(fit_obj, newdata) {
  as.numeric(predict(fit_obj$fit_ald, newdata = newdata, type = "response")[, 1])
}

predict_gpd_params <- function(fit_obj, newdata) {
  p <- predict(fit_obj$fit_gpd, newdata = newdata, type = "response")
  list(scale = as.numeric(p[, 1]), shape = as.numeric(p[, 2]))
}

qgpd_excess <- function(u, scale, shape) {
  pmax(ifelse(abs(shape) < 1e-6,
              -scale * log1p(-u),
              scale / shape * ((1 - u)^(-shape) - 1)), 0)
}

auto_bandwidth <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || !is.finite(stats::sd(x)) || stats::sd(x) == 0) return(0)
  1.06 * stats::sd(x) * length(x)^(-1 / 5)
}

jitter_linear <- function(x, sd, lower = -Inf, upper = Inf) {
  if (!is.finite(sd) || sd <= 0) return(x)
  out <- x + stats::rnorm(length(x), 0, sd)
  bad <- which(out < lower | out > upper | !is.finite(out))
  n_iter <- 0L
  while (length(bad) && n_iter < 200L) {
    out[bad] <- x[bad] + stats::rnorm(length(bad), 0, sd)
    bad <- which(out < lower | out > upper | !is.finite(out))
    n_iter <- n_iter + 1L
  }
  out[bad] <- pmin(upper, pmax(lower, x[bad]))
  out
}

jitter_circular <- function(x, sd, period = 360, origin = 0) {
  if (!is.finite(sd) || sd <= 0) return(x)
  ((x + stats::rnorm(length(x), 0, sd) - origin) %% period) + origin
}

jitter_covariate_values <- function(x0, cv, data, peaks) {
  if (!isTRUE(cv$jitter)) return(x0)
  if (isTRUE(cv$circular)) {
    sd <- if (is.null(cv$jitter_sd)) 10 else cv$jitter_sd
    period <- if (is.null(cv$period)) 360 else cv$period
    origin <- if (is.null(cv$origin)) 0 else cv$origin
    return(jitter_circular(x0, sd = sd, period = period, origin = origin))
  }

  sd <- if (is.null(cv$jitter_sd)) auto_bandwidth(peaks[[cv$name]]) else cv$jitter_sd
  lower <- if (is.null(cv$lower)) min(data[[cv$name]], na.rm = TRUE) else cv$lower
  upper <- if (is.null(cv$upper)) max(data[[cv$name]], na.rm = TRUE) else cv$upper
  jitter_linear(x0, sd = sd, lower = lower, upper = upper)
}

jitter_peak_covariates <- function(base, fit_obj) {
  out <- copy(base)
  peaks <- fit_obj$excess_df
  for (cv in fit_obj$covariates) {
    if (!cv$name %in% names(out)) next
    out[[cv$name]] <- jitter_covariate_values(
      x0 = base[[cv$name]],
      cv = cv,
      data = fit_obj$data,
      peaks = peaks
    )
  }
  out
}

run_simulation_prep <- function(fit_obj) {
  list(
    exceedance_rate = mean(fit_obj$data$excess > 0),
    covariates = fit_obj$covariates,
    response = fit_obj$response
  )
}

infer_observation_years <- function(data, time_col = "time") {
  if (time_col %in% names(data)) {
    tt <- as.POSIXct(data[[time_col]])
    tt <- sort(unique(tt[!is.na(tt)]))
    if (length(tt) >= 2) {
      step_days <- stats::median(as.numeric(diff(tt), units = "days"), na.rm = TRUE)
      span_days <- as.numeric(max(tt) - min(tt), units = "days")
      if (is.finite(step_days) && step_days > 0) span_days <- span_days + step_days
      return(span_days / 365.25)
    }
  }

  if ("obs_year" %in% names(data)) {
    years <- unique(data$obs_year[!is.na(data$obs_year)])
    if (length(years)) return(length(years))
  }

  NA_real_
}

use_observed_data <- function(fit_obj, n_years = NULL) {
  sim <- copy(fit_obj$data)
  response <- fit_obj$response

  if (is.null(n_years)) n_years <- infer_observation_years(sim)
  if (!is.finite(n_years) || n_years <= 0) {
    stop("Could not infer a positive observation duration. Provide n_years explicitly.")
  }

  if ("time" %in% names(sim)) {
    sim[, sim_year := as.integer(format(as.POSIXct(time), "%Y"))]
  } else {
    sim[, sim_year := NA_integer_]
  }

  sim[, extreme := excess > 0]
  peak_idx <- which(sim$extreme)

  if (!length(peak_idx)) {
    return(list(sim = sim, peak_log = data.table(), n_years = n_years, bootstrap = FALSE))
  }

  peak_rows <- copy(sim[peak_idx])
  par_obs <- predict_gpd_params(fit_obj, peak_rows)
  set(sim, peak_idx, "gpd_scale", par_obs$scale)
  set(sim, peak_idx, "gpd_shape", par_obs$shape)

  peak_log <- data.table(
    sim_row = peak_idx,
    sim_year = sim$sim_year[peak_idx],
    response_observed = peak_rows[[response]],
    response_est_jitter = peak_rows[[response]],
    response_est_obscov = peak_rows[[response]],
    thresh_obs = peak_rows$threshold,
    thresh_jit = peak_rows$threshold,
    scale_obs = par_obs$scale,
    scale_jit = par_obs$scale,
    shape_obs = par_obs$shape,
    shape_jit = par_obs$shape,
    u_draw = NA_real_
  )

  for (cv in fit_obj$covariates) {
    if (cv$name %in% names(peak_rows)) {
      peak_log[, (paste0(cv$name, "_obs")) := peak_rows[[cv$name]]]
      peak_log[, (paste0(cv$name, "_jit")) := peak_rows[[cv$name]]]
    }
  }

  list(sim = sim, peak_log = peak_log, n_years = n_years, bootstrap = FALSE)
}

simulate_block_resample <- function(fit_obj, pred_grids = NULL, n_years, seed = 42,
                                    h_cspd = NULL, sigma_ang_deg = NULL) {
  set.seed(seed)
  data <- copy(fit_obj$data)
  response <- fit_obj$response

  data[, obs_year := as.integer(format(as.POSIXct(time), "%Y"))]
  available_years <- unique(data$obs_year)

  sim <- rbindlist(lapply(seq_len(n_years), function(i) {
    x <- data[obs_year == sample(available_years, 1)]
    x[, sim_year := i]
    x
  }))

  if (!is.null(h_cspd)) {
    i <- match("Cspd", vapply(fit_obj$covariates, `[[`, character(1), "name"))
    if (!is.na(i)) fit_obj$covariates[[i]]$jitter_sd <- h_cspd
  }
  if (!is.null(sigma_ang_deg)) {
    i <- match("angularDifference", vapply(fit_obj$covariates, `[[`, character(1), "name"))
    if (!is.na(i)) fit_obj$covariates[[i]]$jitter_sd <- sigma_ang_deg
  }

  peak_idx <- which(sim$excess > 0)
  n_peaks <- length(peak_idx)
  sim[, extreme := FALSE]

  if (!n_peaks) {
    return(list(sim = sim, peak_log = data.table(), n_years = n_years, bootstrap = TRUE))
  }

  obs_rows <- copy(sim[peak_idx])
  star_rows <- jitter_peak_covariates(obs_rows, fit_obj)

  thresh_obs <- predict_threshold(fit_obj, obs_rows)
  par_obs <- predict_gpd_params(fit_obj, obs_rows)
  thresh_star <- predict_threshold(fit_obj, star_rows)
  par_star <- predict_gpd_params(fit_obj, star_rows)

  u <- stats::runif(n_peaks)
  excess_star <- qgpd_excess(u, par_star$scale, par_star$shape)
  excess_obs <- qgpd_excess(u, par_obs$scale, par_obs$shape)

  response_star <- thresh_star + excess_star
  response_obs_cov <- thresh_obs + excess_obs

  for (cv in fit_obj$covariates) {
    if (cv$name %in% names(star_rows)) {
      set(sim, peak_idx, cv$name, star_rows[[cv$name]])
    }
  }
  set(sim, peak_idx, "threshold", thresh_star)
  set(sim, peak_idx, response, response_star)
  set(sim, peak_idx, "extreme", TRUE)
  set(sim, peak_idx, "gpd_scale", par_star$scale)
  set(sim, peak_idx, "gpd_shape", par_star$shape)

  peak_log <- data.table(
    sim_row = peak_idx,
    sim_year = sim$sim_year[peak_idx],
    response_observed = obs_rows[[response]],
    response_est_jitter = response_star,
    response_est_obscov = response_obs_cov,
    thresh_obs = thresh_obs,
    thresh_jit = thresh_star,
    scale_obs = par_obs$scale,
    scale_jit = par_star$scale,
    shape_obs = par_obs$shape,
    shape_jit = par_star$shape,
    u_draw = u
  )

  for (cv in fit_obj$covariates) {
    if (cv$name %in% names(obs_rows)) {
      peak_log[, (paste0(cv$name, "_obs")) := obs_rows[[cv$name]]]
      peak_log[, (paste0(cv$name, "_jit")) := star_rows[[cv$name]]]
    }
  }

  list(sim = sim, peak_log = peak_log, n_years = n_years, bootstrap = TRUE)
}

quick_sim <- function(fit_obj, n_years = 10, seed = 1) {
  prep <- run_simulation_prep(fit_obj)
  simulate_block_resample(fit_obj, prep, n_years, seed)
}
