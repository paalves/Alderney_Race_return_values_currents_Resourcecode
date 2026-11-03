library(data.table)
library(evgam)
library(extRemes)
library(parallel)

EVGAM_TAU <- 0.96
EVGAM_DECLUSTER_R <- 72
EVGAM_N_FOLDS <- 8L
EVGAM_LOG_SCORE_CAP <- 30
EVGAM_ALD <- "s(Cspd, bs = 'tp', k = 5) + s(angleDiff, bs = 'cp', k = 5, m = 1)"
EVGAM_SCALE <- "s(Cspd, bs = 'ts', k = 10) + s(angleDiff, bs = 'cc', k = 10)"
EVGAM_CANDIDATES <- list(
  retained = c(paste("excess ~", EVGAM_SCALE), "~ 1"),
  constant = c("excess ~ 1", "~ 1"),
  interaction = c(paste("excess ~", EVGAM_SCALE, "+ ti(Cspd, angleDiff, bs = c('tp', 'cc'), k = c(5, 5), m = c(1, 1))"), "~ 1"),
  smooth_shape = c(paste("excess ~", EVGAM_SCALE), "~ s(Cspd, bs = 'ts', k = 5)")
)

#' Add the year (April-March) and the cross-validation fold to the data.
evgam_folds <- function(data, n_folds = EVGAM_N_FOLDS) {
  d <- as.data.table(data)[, .(time, hs, Os, Cspd, angleDiff)]
  d[, year := as.integer(format(time - 90 * 86400, "%Y"))]
  d[, fold := (match(year, sort(unique(year))) - 1L) %% n_folds + 1L]
  d[]
}

#' Fit the ALD quantile-regression threshold on the selected rows.
evgam_threshold <- function(d, response, rows, tau = EVGAM_TAU) {
  evgam(list(as.formula(paste(response, "~", EVGAM_ALD)), as.formula(paste("~", EVGAM_ALD))),
        d[rows], family = "ald", ald.args = list(tau = tau))
}

#' Decluster the exceedances of the fitted threshold and return the cluster maxima.
evgam_peaks <- function(d, response, fit_ald, r = EVGAM_DECLUSTER_R) {
  u <- predict(fit_ald, newdata = d, type = "response")$location
  ex <- pmax(as.numeric(decluster(d[[response]], threshold = u, r = r)) - u, 0)
  i <- which(ex > 0)
  data.table(d[i, .(time, Cspd, angleDiff, year, fold)], y = d[[response]][i], threshold = u[i], excess = ex[i])
}

#' Fit a GPD with evgam to the threshold excesses.
evgam_fit_gpd <- function(forms, peaks) {
  evgam(lapply(forms, as.formula), data = peaks, family = "gpd")
}

#' Logarithmic score (negative GPD log-density) of excesses, capped beyond the upper endpoint.
gpd_log_score <- function(excess, scale, shape, cap = EVGAM_LOG_SCORE_CAP) {
  z <- 1 + shape * excess / scale
  ld <- ifelse(abs(shape) < 1e-8, -log(scale) - excess / scale,
               ifelse(z > 0, -log(scale) - (1 + 1 / shape) * log(pmax(z, 1e-300)), -Inf))
  -pmax(ld, -cap)
}

#' Refit the threshold and candidate GPDs without fold k and score the held-out cluster maxima.
evgam_cv_fold <- function(d, response, k) {
  pk <- evgam_peaks(d, response, evgam_threshold(d, response, d$fold != k))
  train <- pk[fold != k]; test <- pk[fold == k]
  rbindlist(lapply(names(EVGAM_CANDIDATES), function(nm) {
    p <- predict(evgam_fit_gpd(EVGAM_CANDIDATES[[nm]], train), newdata = test, type = "response")
    data.table(response, fold = k, model = nm, log_score = mean(gpd_log_score(test$excess, p$scale, p$shape)))
  }))
}

#' Block cross-validation of the candidate models and full-data fit of the retained model.
select_evgam <- function(data, response, n_cores = 6L) {
  d <- evgam_folds(data)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, { library(data.table); library(evgam); library(extRemes) }))
  clusterExport(cl, c("EVGAM_TAU", "EVGAM_DECLUSTER_R", "EVGAM_LOG_SCORE_CAP", "EVGAM_ALD", "EVGAM_SCALE",
                      "EVGAM_CANDIDATES", "evgam_threshold", "evgam_peaks", "evgam_fit_gpd", "gpd_log_score",
                      "evgam_cv_fold"))
  cv <- rbindlist(parLapply(cl, seq_len(EVGAM_N_FOLDS), function(k, d, response) evgam_cv_fold(d, response, k),
                            d = d, response = response))
  summary <- cv[, .(log_score = mean(log_score), se = sd(log_score) / sqrt(.N)), by = model]
  fit_ald <- evgam_threshold(d, response, rep(TRUE, nrow(d)))
  peaks <- evgam_peaks(d, response, fit_ald)
  list(best_fit = list(response = response, fit_ald = fit_ald, fit_gpd = evgam_fit_gpd(EVGAM_CANDIDATES$retained, peaks),
                       excess_df = peaks, formulas = list(ald = EVGAM_ALD, gpd = EVGAM_CANDIDATES$retained)),
       cv = cv, summary = summary)
}
