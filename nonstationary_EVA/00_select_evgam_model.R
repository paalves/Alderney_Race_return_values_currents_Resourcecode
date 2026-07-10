if (!exists("root")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                        error = function(e) NA_character_)
  root <- if (is.na(this_file)) getwd() else dirname(this_file)
}

source(file.path(root, "01_evgam_fit.R"))

selection_safe_name <- function(x) {
  gsub("_+$", "", gsub("^_+", "", gsub("[^A-Za-z0-9_]+", "_", x)))
}

selection_value <- function(x, name, default = NULL) {
  if (!is.null(x[[name]])) x[[name]] else default
}

selection_r_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (identical(x, TRUE)) return("TRUE")
  if (identical(x, FALSE)) return("FALSE")
  if (is.character(x)) {
    if (length(x) == 1) return(sprintf("\"%s\"", x))
    return(sprintf("c(%s)", paste(sprintf("\"%s\"", x), collapse = ", ")))
  }
  if (is.numeric(x) || is.integer(x)) {
    if (length(x) == 1) return(format(x, trim = TRUE, scientific = FALSE))
    return(sprintf("c(%s)", paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")))
  }
  deparse(x)
}

selection_formula_text <- function(x) {
  paste(deparse(x), collapse = " ")
}

as_covariate_names <- function(covariate1, covariate2 = NULL) {
  covariates <- c(covariate1, covariate2)
  covariates <- covariates[!is.na(covariates) & nzchar(covariates)]
  if (!length(covariates)) stop("Provide at least covariate1.")
  covariates
}

as_options <- function(x) {
  if (is.null(x)) return(list(NULL))
  if (is.list(x) && !is.data.frame(x)) return(x)
  as.list(x)
}

param_by_covariate <- function(x, n_cov, name) {
  if (is.null(x)) return(rep(list(NULL), n_cov))
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) == 1 && n_cov > 1) x <- rep(x, n_cov)
    if (length(x) != n_cov) {
      stop(name, " must have length 1 or one element per covariate.")
    }
    return(x)
  }
  if (length(x) == 1) return(rep(list(x), n_cov))
  if (length(x) == n_cov) return(as.list(x))
  stop(name, " must have length 1 or one value per covariate.")
}

fixed_param_by_covariate <- function(x, n_cov, name) {
  out <- param_by_covariate(x, n_cov, name)
  lapply(out, function(z) {
    if (is.null(z)) return(NULL)
    if (length(z) != 1) stop(name, " is fixed for ALD: use one value per covariate.")
    z[[1]]
  })
}

one_cov_grid <- function(bs, k, m) {
  out <- list()
  i <- 0L
  for (bb in as_options(bs)) {
    for (kk in as_options(k)) {
      for (mm in as_options(m)) {
        i <- i + 1L
        out[[i]] <- list(bs = bb, k = kk, m = mm)
      }
    }
  }
  out
}

combine_cov_grids <- function(per_cov) {
  out <- list(list())
  for (opts in per_cov) {
    out <- unlist(lapply(out, function(prefix) {
      lapply(opts, function(opt) c(prefix, list(opt)))
    }), recursive = FALSE)
  }
  lapply(out, function(combo) {
    list(
      bs = lapply(combo, `[[`, "bs"),
      k = lapply(combo, `[[`, "k"),
      m = lapply(combo, `[[`, "m")
    )
  })
}

candidate_param_grid <- function(bs, k, m, n_cov, prefix) {
  bs <- param_by_covariate(bs, n_cov, paste0(prefix, "_bs"))
  k <- param_by_covariate(k, n_cov, paste0(prefix, "_k"))
  m <- param_by_covariate(m, n_cov, paste0(prefix, "_m"))
  per_cov <- lapply(seq_len(n_cov), function(i) one_cov_grid(bs[[i]], k[[i]], m[[i]]))
  combine_cov_grids(per_cov)
}

make_covariate_specs <- function(covariates, ald, gpd_location, gpd_shape) {
  lapply(seq_along(covariates), function(i) {
    covariate_spec(
      name = covariates[i],
      ald_bs = ald$bs[[i]],
      ald_k = ald$k[[i]],
      ald_m = ald$m[[i]],
      gpd_bs = gpd_location$bs[[i]],
      gpd_k = gpd_location$k[[i]],
      gpd_m = gpd_location$m[[i]],
      shape_bs = gpd_shape$bs[[i]],
      shape_k = gpd_shape$k[[i]],
      shape_m = gpd_shape$m[[i]],
      lower = 0
    )
  })
}

shape_covariate_options <- function(x, covariates) {
  if (identical(x, "all")) return(list(covariates))
  if (identical(x, "none")) return(list(FALSE))
  if (is.list(x) && !is.data.frame(x)) return(x)
  list(x)
}

covariate_spec_to_code <- function(cv, indent = "    ") {
  fields <- c(
    "name",
    "source",
    "circular",
    "period",
    "origin",
    "ald_bs", "ald_k", "ald_m",
    "gpd_bs", "gpd_k", "gpd_m",
    "shape_bs", "shape_k", "shape_m",
    "jitter_sd",
    "upper",
    "lower"
  )
  fields <- fields[fields %in% names(cv)]
  lines <- c(
    paste0(indent, "covariate_spec("),
    paste0(
      indent, "  ",
      fields, " = ", vapply(fields, function(nm) selection_r_value(cv[[nm]]), character(1)),
      c(rep(",", length(fields) - 1), "")
    ),
    paste0(indent, ")")
  )
  lines
}

candidate_config_code <- function(candidate, object_name = "covariates_hs_selected") {
  covariates <- candidate$covariates
  cov_lines <- unlist(lapply(seq_along(covariates), function(i) {
    x <- covariate_spec_to_code(covariates[[i]])
    if (i < length(covariates)) x[length(x)] <- paste0(x[length(x)], ",")
    x
  }), use.names = FALSE)

  paste(c(
    paste0(object_name, " <- list("),
    cov_lines,
    ")",
    "",
    sprintf("gpd_interaction = %s,", selection_r_value(selection_value(candidate, "gpd_interaction", TRUE))),
    sprintf("gpd_interaction_k = %s,", selection_r_value(selection_value(candidate, "gpd_interaction_k", NULL))),
    sprintf("gpd_shape_covariates = %s", selection_r_value(selection_value(candidate, "gpd_shape_covariates", NULL)))
  ), collapse = "\n")
}

make_evgam_candidates <- function(covariate1,
                                  covariate2 = NULL,
                                  ald_bs = "tp",
                                  ald_k = 5,
                                  ald_m = NULL,
                                  gpd_location_bs = "cr",
                                  gpd_location_k = 10,
                                  gpd_location_m = 1,
                                  gpd_shape_bs = "cr",
                                  gpd_shape_k = 3,
                                  gpd_shape_m = NULL,
                                  gpd_interaction = c(FALSE, TRUE),
                                  gpd_interaction_k = list(NULL, c(5, 5)),
                                  gpd_shape_covariates = "all") {
  covariates <- as_covariate_names(covariate1, covariate2)
  n_cov <- length(covariates)

  ald <- list(
    bs = fixed_param_by_covariate(ald_bs, n_cov, "ald_bs"),
    k = fixed_param_by_covariate(ald_k, n_cov, "ald_k"),
    m = fixed_param_by_covariate(ald_m, n_cov, "ald_m")
  )
  gpd_location_grid <- candidate_param_grid(
    gpd_location_bs, gpd_location_k, gpd_location_m, n_cov, "gpd_location"
  )
  gpd_shape_grid <- candidate_param_grid(
    gpd_shape_bs, gpd_shape_k, gpd_shape_m, n_cov, "gpd_shape"
  )
  shape_options <- shape_covariate_options(gpd_shape_covariates, covariates)
  interaction_k_options <- as_options(gpd_interaction_k)

  candidates <- list()
  for (loc in gpd_location_grid) {
    for (shape in gpd_shape_grid) {
      for (inter in gpd_interaction) {
        ik_opts <- if (isTRUE(inter)) interaction_k_options else list(NULL)
        for (ik in ik_opts) {
          for (shape_vars in shape_options) {
            shape_label <- if (identical(shape_vars, FALSE) ||
                               identical(shape_vars, character()) ||
                               is.null(shape_vars)) {
              "constshape"
            } else {
              paste(shape_vars, collapse = "_")
            }
            label <- paste(
              paste0("ald_", paste(covariates, unlist(ald$bs), unlist(ald$k), sep = "", collapse = "_")),
              paste0("gpd_", paste(unlist(loc$bs), unlist(loc$k), sep = "", collapse = "_")),
              paste0("shape_", paste(unlist(shape$bs), unlist(shape$k), sep = "", collapse = "_"), "_", shape_label),
              paste0("inter_", isTRUE(inter), "_ik_", if (is.null(ik)) "none" else paste(ik, collapse = "_")),
              sep = "__"
            )
            candidates[[selection_safe_name(label)]] <- list(
              label = label,
              covariates = make_covariate_specs(covariates, ald, loc, shape),
              gpd_interaction = isTRUE(inter),
              gpd_interaction_k = ik,
              gpd_shape_covariates = shape_vars
            )
          }
        }
      }
    }
  }

  candidates
}

make_cspd_tp_evgam_candidates <- make_evgam_candidates
make_cspd_tp_gpd_candidates <- make_evgam_candidates

gpd_pit <- function(excess, scale, shape, eps = 1e-10) {
  out <- rep(NA_real_, length(excess))
  ok <- is.finite(excess) & is.finite(scale) & is.finite(shape) & scale > 0
  near_zero <- ok & abs(shape) < 1e-8
  out[near_zero] <- 1 - exp(-excess[near_zero] / scale[near_zero])

  curved <- ok & !near_zero
  curved_idx <- which(curved)
  z <- 1 + shape[curved_idx] * excess[curved_idx] / scale[curved_idx]
  tmp <- rep(NA_real_, length(curved_idx))
  valid <- z > 0
  tmp[valid] <- 1 - z[valid]^(-1 / shape[curved_idx][valid])
  out[curved_idx] <- tmp

  pmin(1 - eps, pmax(eps, out))
}

score_gpd_qq <- function(fit_obj,
                         response = fit_obj$response) {
  exc <- fit_obj$excess_df
  if (is.null(exc) || !nrow(exc)) stop("No GPD exceedances available for QQ scoring.")

  pars <- predict(fit_obj$fit_gpd, newdata = exc, type = "response")
  pit <- gpd_pit(
    excess = exc$excess,
    scale = as.numeric(pars[, 1]),
    shape = as.numeric(pars[, 2])
  )
  resid <- -log1p(-pit)

  ok <- is.finite(resid)
  resid <- resid[ok]
  if (length(resid) < 5) stop("Too few finite QQ residuals to score the model.")

  o <- order(resid)
  resid <- resid[o]
  expected <- stats::qexp(stats::ppoints(length(resid)))

  err <- resid - expected
  score <- sum(abs(err))

  list(
    score = score,
    n = length(resid),
    qq = data.table(
      expected = expected,
      observed = resid,
      error = err
    )
  )
}

plot_gpd_qq <- function(score_obj, main = "GPD QQ diagnostic") {
  qq <- score_obj$qq
  lim <- range(c(qq$expected, qq$observed), finite = TRUE)

  plot(qq$expected, qq$observed,
       pch = 16, cex = 0.8,
       col = grDevices::adjustcolor("black", 0.45),
       xlim = lim, ylim = lim,
       xlab = "Theoretical",
       ylab = "Empirical",
       main = main,
       cex.lab = 1.35,
       cex.axis = 1.15,
       cex.main = 1.15)
  abline(0, 1, col = "red3", lwd = 2)
  invisible(score_obj)
}

fit_fixed_ald_for_selection <- function(csv_path,
                                        response = "hs",
                                        tau = 0.96,
                                        decluster_r = 72,
                                        ald_covariates,
                                        ald_interaction = FALSE,
                                        ald_interaction_k = NULL,
                                        make_plot = FALSE) {
  ald_covariates <- normalise_covariate_specs(ald_covariates)
  d <- prepare_eva_data(fread(csv_path), ald_covariates)

  f_ald <- make_location_formula(
    response = response,
    covariates = ald_covariates,
    prefix = "ald",
    include_interaction = ald_interaction,
    interaction_k = ald_interaction_k
  )

  fit_ald <- evgam(list(f_ald, ~ 1), d, family = "ald", ald.args = list(tau = tau))
  p_ald <- NULL
  if (isTRUE(make_plot)) {
    plot(fit_ald, scheme = 2)
    p_ald <- recordPlot()
  }

  d[, threshold := predict(fit_ald, newdata = d, type = "response")[, 1]]
  d[, excess := pmax(decluster(.SD[[response]], threshold = threshold, r = decluster_r) - threshold, 0)]
  exc <- d[excess > 0]

  list(
    data = d,
    fit_ald = fit_ald,
    excess_df = exc,
    response = response,
    covariates = ald_covariates,
    formulas = list(ald = f_ald),
    plots = list(ald_splines = p_ald),
    tau = tau,
    decluster_r = decluster_r
  )
}

fit_one_gpd_candidate <- function(ald_fit, candidate) {
  response <- ald_fit$response
  covariates <- normalise_covariate_specs(candidate$covariates)
  d <- prepare_eva_data(copy(ald_fit$data), covariates)
  exc <- d[excess > 0]

  f_gpd_loc <- make_location_formula(
    response = response,
    covariates = covariates,
    prefix = "gpd",
    include_interaction = selection_value(candidate, "gpd_interaction", TRUE),
    interaction_k = selection_value(candidate, "gpd_interaction_k", NULL),
    lhs_override = "excess"
  )
  f_gpd_shape <- make_scale_formula(
    covariates,
    selection_value(candidate, "gpd_shape_covariates", NULL)
  )

  fit_gpd <- evgam(
    list(f_gpd_loc, f_gpd_shape),
    data = exc,
    family = "gpd"
  )

  qq <- predict(fit_gpd, newdata = exc, type = "qqplot")
  fit <- list(
    data = d,
    fit_ald = ald_fit$fit_ald,
    fit_gpd = fit_gpd,
    excess_df = exc,
    response = response,
    covariates = covariates,
    formulas = list(
      ald = ald_fit$formulas$ald,
      gpd_location = f_gpd_loc,
      gpd_shape = f_gpd_shape,
      gpd_scale = f_gpd_shape
    ),
    plots = list(
      ald_splines = ald_fit$plots$ald_splines,
      gpd_splines = NULL,
      qq_gpd = qq
    )
  )

  qq_score <- score_gpd_qq(fit_obj = fit)

  list(fit = fit, qq_score = qq_score)
}

select_evgam_by_qq <- function(csv_path,
                               candidates,
                               response = "hs",
                               tau = 0.96,
                               decluster_r = 72,
                               ald_covariates = NULL,
                               ald_interaction = FALSE,
                               ald_interaction_k = NULL,
                               plot_ald = FALSE,
                               plot_best = TRUE,
                               verbose = TRUE) {
  if (is.null(names(candidates))) names(candidates) <- paste0("model_", seq_along(candidates))
  if (!length(candidates)) stop("At least one GPD candidate is required.")
  if (is.null(ald_covariates)) ald_covariates <- candidates[[1]]$covariates

  if (isTRUE(verbose)) message("Fitting ALD model...")
  fixed_ald <- fit_fixed_ald_for_selection(
    csv_path = csv_path,
    response = response,
    tau = tau,
    decluster_r = decluster_r,
    ald_covariates = ald_covariates,
    ald_interaction = ald_interaction,
    ald_interaction_k = ald_interaction_k,
    make_plot = plot_ald
  )

  summary_rows <- vector("list", length(candidates))
  names(summary_rows) <- names(candidates)
  best <- NULL
  best_score_value <- Inf

  if (isTRUE(verbose)) message("Fitting GPDs...")
  for (i in seq_along(candidates)) {
    nm <- names(candidates)[i]
    candidate <- candidates[[nm]]
    label <- selection_value(candidate, "label", nm)

    summary_rows[[i]] <- tryCatch({
      fitted <- fit_one_gpd_candidate(
        ald_fit = fixed_ald,
        candidate = candidate
      )

      row <- data.table(
        name = nm,
        label = label,
        ok = TRUE,
        score = fitted$qq_score$score,
        n_excess = fitted$qq_score$n,
        ald = selection_formula_text(fitted$fit$formulas$ald),
        gpd_location = selection_formula_text(fitted$fit$formulas$gpd_location),
        gpd_shape = selection_formula_text(fitted$fit$formulas$gpd_shape),
        error = NA_character_
      )

      if (is.finite(fitted$qq_score$score) && fitted$qq_score$score < best_score_value) {
        if (!is.null(best)) {
          best$fit <- NULL
          best$qq_score <- NULL
          best$candidate <- NULL
          invisible(gc(verbose = FALSE))
        }
        best_score_value <- fitted$qq_score$score
        best <- list(
          name = nm,
          label = label,
          candidate = candidate,
          fit = fitted$fit,
          qq_score = fitted$qq_score
        )
      }

      rm(fitted)
      invisible(gc(verbose = FALSE))
      row
    }, error = function(e) {
      data.table(
        name = nm,
        label = label,
        ok = FALSE,
        score = NA_real_,
        n_excess = NA_integer_,
        ald = NA_character_,
        gpd_location = NA_character_,
        gpd_shape = NA_character_,
        error = conditionMessage(e)
      )
    })
  }

  summary <- rbindlist(summary_rows, fill = TRUE)

  setorder(summary, ok, score, na.last = TRUE)
  summary[, rank := ifelse(ok, frank(score, ties.method = "first", na.last = "keep"), NA_integer_)]
  setorder(summary, rank, na.last = TRUE)

  if (is.null(best)) stop("No EVGAM candidate fitted successfully.")
  best_name <- best$name

  best_config_code <- candidate_config_code(best$candidate)

  if (isTRUE(plot_best)) {
    plot_gpd_qq(
      best$qq_score,
      main = paste0("Best GPD QQ: ", best$label)
    )
    plot(best$fit$fit_gpd, scheme = 2, shade = FALSE)
    par(mfrow = c(1, 1))
  }

  invisible(list(
    best_name = best_name,
    best_label = best$label,
    best_fit = best$fit,
    best_candidate = best$candidate,
    best_config_code = best_config_code,
    best_score = best$qq_score,
    summary = summary
  ))
}

if (FALSE) {
  candidates <- make_evgam_candidates(
    covariate1 = "Wspd",
    covariate2 = "tp",
    ald_bs = c("tp", "tp"),
    ald_k = c(5, 5),
    ald_m = list(NULL, NULL),
    gpd_location_bs = list(c("cr", "tp"), c("cr", "tp")),
    gpd_location_k = list(15, c(10, 15), c(15,15)),
    gpd_location_m = list(1, 1),
    gpd_shape_covariates = "Wspd",
    gpd_shape_bs = list("cr", NULL),
    gpd_shape_k = list(5, NULL),
    gpd_shape_m = list(NULL, NULL),
    gpd_interaction = c(FALSE, TRUE),
    gpd_interaction_k = list(NULL, c(5, 5))
  )

  selection <- select_evgam_by_qq(
    csv_path = data_file(),
    candidates = candidates,
    response = "hs",
    tau = 0.96,
    decluster_r = 72,
    plot_best = TRUE
  )

  selection$best_fit
  selection$best_config_code
}
