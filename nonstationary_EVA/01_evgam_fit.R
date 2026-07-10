library(data.table)
library(evgam)
library(extRemes)

covariate_spec <- function(name,
                           source = NULL,
                           circular = FALSE,
                           period = NULL,
                           origin = NULL,
                           ald_bs = NULL,
                           ald_k = NULL,
                           ald_m = NULL,
                           gpd_bs = NULL,
                           gpd_k = NULL,
                           gpd_m = NULL,
                           shape_bs = NULL,
                           shape_k = NULL,
                           shape_m = NULL,
                           scale_bs = NULL,
                           scale_k = NULL,
                           scale_m = NULL,
                           jitter = TRUE,
                           jitter_sd = NULL,
                           lower = NULL,
                           upper = NULL) {
  Filter(Negate(is.null), list(
    name = name,
    source = source,
    circular = circular,
    period = period,
    origin = origin,
    ald_bs = ald_bs,
    ald_k = ald_k,
    ald_m = ald_m,
    gpd_bs = gpd_bs,
    gpd_k = gpd_k,
    gpd_m = gpd_m,
    shape_bs = shape_bs,
    shape_k = shape_k,
    shape_m = shape_m,
    scale_bs = scale_bs,
    scale_k = scale_k,
    scale_m = scale_m,
    jitter = jitter,
    jitter_sd = jitter_sd,
    lower = lower,
    upper = upper
  ))
}

default_covariate_specs <- function() {
  list(
    covariate_spec(
      name = "Cspd",
      circular = FALSE,
      ald_bs = "tp", ald_k = 5,
      gpd_bs = "cr", gpd_k = 15, gpd_m = 1,
      scale_bs = "cr", scale_k = 3,
      jitter = TRUE,
      lower = 0
    ),
    covariate_spec(
      name = "angularDifference",
      source = "angleDiff",
      circular = TRUE,
      period = 360,
      ald_bs = "cp", ald_k = 5, ald_m = 1,
      gpd_bs = "cc", gpd_k = 15, gpd_m = 1,
      jitter = TRUE, jitter_sd = 10
    )
  )
}

normalise_covariate_specs <- function(covariates = NULL) {
  if (is.null(covariates)) covariates <- default_covariate_specs()
  lapply(covariates, function(x) {
    if (is.character(x)) x <- list(name = x)
    x$circular <- isTRUE(x$circular)
    if (is.null(x$period) && x$circular) x$period <- 360
    if (is.null(x$jitter)) x$jitter <- TRUE
    x
  })
}

prepare_eva_data <- function(d, covariates) {
  for (cv in covariates) {
    if (!cv$name %in% names(d) && !is.null(cv$source) && cv$source %in% names(d)) {
      d[, (cv$name) := get(cv$source)]
    }
    if (isTRUE(cv$circular) && cv$name %in% names(d)) {
      period <- if (is.null(cv$period)) 360 else cv$period
      origin <- if (is.null(cv$origin)) 0 else cv$origin
      d[, (cv$name) := ((get(cv$name) - origin) %% period) + origin]
    }
  }
  d
}

fmt_value <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 1) return(as.character(x))
  paste0("c(", paste(x, collapse = ", "), ")")
}

smooth_term <- function(cv, prefix = "gpd") {
  bs <- cv[[paste0(prefix, "_bs")]]
  k <- cv[[paste0(prefix, "_k")]]
  m <- cv[[paste0(prefix, "_m")]]
  if (identical(prefix, "scale")) {
    if (is.null(bs)) bs <- cv$shape_bs
    if (is.null(k)) k <- cv$shape_k
    if (is.null(m)) m <- cv$shape_m
  }

  if (is.null(bs)) bs <- if (isTRUE(cv$circular)) "cc" else "cr"
  if (is.null(k)) k <- if (identical(prefix, "scale")) 3 else 10

  args <- c(
    cv$name,
    sprintf("bs = '%s'", bs),
    sprintf("k = %s", fmt_value(k))
  )
  if (!is.null(m)) args <- c(args, sprintf("m = %s", fmt_value(m)))
  sprintf("s(%s)", paste(args, collapse = ", "))
}

tensor_term <- function(covariates, prefix = "gpd", k = NULL) {
  if (length(covariates) < 2) return(character())
  bs <- vapply(covariates, function(cv) {
    x <- cv[[paste0(prefix, "_bs")]]
    if (is.null(x)) x <- if (isTRUE(cv$circular)) "cc" else "cr"
    x
  }, character(1))
  m <- lapply(covariates, function(cv) cv[[paste0(prefix, "_m")]])
  names <- vapply(covariates, `[[`, character(1), "name")
  if (is.null(k)) k <- rep(5, length(covariates))

  args <- c(
    names,
    sprintf("bs = c(%s)", paste(sprintf("'%s'", bs), collapse = ", ")),
    sprintf("k = c(%s)", paste(k, collapse = ", "))
  )
  if (all(vapply(m, Negate(is.null), logical(1)))) {
    args <- c(args, sprintf("m = c(%s)", paste(unlist(m), collapse = ", ")))
  }
  sprintf("ti(%s)", paste(args, collapse = ", "))
}

make_location_formula <- function(response, covariates, prefix,
                                  include_interaction = FALSE,
                                  interaction_k = NULL,
                                  lhs_override = NULL) {
  terms <- vapply(covariates, smooth_term, character(1), prefix = prefix)
  if (isTRUE(include_interaction) && length(covariates) >= 2) {
    terms <- c(terms, tensor_term(covariates, prefix = prefix, k = interaction_k))
  }
  lhs <- if (is.null(lhs_override)) response else lhs_override
  as.formula(sprintf("%s ~ %s", lhs, paste(terms, collapse = " + ")))
}

make_scale_formula <- function(covariates, scale_covariates = NULL) {
  if (is.null(scale_covariates)) {
    if ("Cspd" %in% vapply(covariates, `[[`, character(1), "name")) {
      scale_covariates <- "Cspd"
    } else {
      return(~ 1)
    }
  }
  if (identical(scale_covariates, character()) || identical(scale_covariates, FALSE)) {
    return(~ 1)
  }
  keep <- vapply(covariates, function(cv) cv$name %in% scale_covariates, logical(1))
  terms <- vapply(covariates[keep], smooth_term, character(1), prefix = "scale")
  if (!length(terms)) return(~ 1)
  as.formula(sprintf("~ %s", paste(terms, collapse = " + ")))
}

run_evgam_fit <- function(csv_path,
                          tau = 0.96,
                          decluster_r = 72,
                          response = "hs",
                          covariates = NULL,
                          ald_interaction = FALSE,
                          ald_interaction_k = NULL,
                          gpd_interaction = TRUE,
                          gpd_interaction_k = NULL,
                          gpd_shape_covariates = NULL,
                          gpd_scale_covariates = NULL,
                          make_plots = TRUE,
                          verbose = TRUE) {
  covariates <- normalise_covariate_specs(covariates)
  d <- prepare_eva_data(fread(csv_path), covariates)

  f_ald <- make_location_formula(
    response = response,
    covariates = covariates,
    prefix = "ald",
    include_interaction = ald_interaction,
    interaction_k = ald_interaction_k
  )

  if (isTRUE(verbose)) message("Fitting ALD model...")
  fit_ald <- evgam(list(f_ald, ~ 1), d, family = "ald", ald.args = list(tau = tau))
  p_ald <- NULL
  if (isTRUE(make_plots)) {
    plot(fit_ald, scheme = 2)
    p_ald <- recordPlot()
  }

  d[, threshold := predict(fit_ald, newdata = d, type = "response")[, 1]]
  d[, excess := pmax(decluster(.SD[[response]], threshold = threshold, r = decluster_r) - threshold, 0)]
  exc <- d[excess > 0]

  f_gpd_loc <- make_location_formula(
    response = response,
    covariates = covariates,
    prefix = "gpd",
    include_interaction = gpd_interaction,
    interaction_k = gpd_interaction_k,
    lhs_override = "excess"
  )
  if (is.null(gpd_shape_covariates)) gpd_shape_covariates <- gpd_scale_covariates
  f_gpd_scale <- make_scale_formula(covariates, gpd_shape_covariates)

  if (isTRUE(verbose)) message("Fitting GPDs...")
  fit_gpd <- evgam(
    list(f_gpd_loc, f_gpd_scale),
    data = exc,
    family = "gpd"
  )

  qq <- predict(fit_gpd, newdata = exc, type = "qqplot")
  p_gpd <- NULL
  if (isTRUE(make_plots)) {
    plot(fit_gpd, scheme = 2, shade = FALSE)
    p_gpd <- recordPlot()
    par(mfrow = c(1, 1))
  }

  list(
    data = d,
    fit_ald = fit_ald,
    fit_gpd = fit_gpd,
    excess_df = exc,
    response = response,
    covariates = covariates,
    formulas = list(
      ald = f_ald,
      gpd_location = f_gpd_loc,
      gpd_shape = f_gpd_scale,
      gpd_scale = f_gpd_scale
    ),
    plots = list(ald_splines = p_ald, gpd_splines = p_gpd, qq_gpd = qq)
  )
}
