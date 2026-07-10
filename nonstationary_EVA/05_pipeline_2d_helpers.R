cfg_value_2d <- function(x, name, default = NULL) {
  if (!is.null(x[[name]])) x[[name]] else default
}

scenario_value_2d <- function(scenario, cfg, name, default = NULL) {
  if (!is.null(scenario[[name]])) scenario[[name]] else cfg_value_2d(cfg, name, default)
}

safe_name_2d <- function(x) {
  gsub("_+$", "", gsub("^_+", "", gsub("[^A-Za-z0-9_]+", "_", x)))
}

normalise_diform_2d_scenarios <- function(cfg) {
  if (is.null(cfg$scenarios)) {
    return(list(default = list(
      model = "gam",
      bootstrap = cfg_value_2d(cfg, "bootstrap", TRUE),
      N_YEARS = cfg_value_2d(cfg, "N_YEARS"),
      observed_n_years = cfg_value_2d(cfg, "observed_n_years"),
      seed = cfg_value_2d(cfg, "seed", 42)
    )))
  }

  scenarios <- cfg$scenarios
  if (is.null(names(scenarios))) names(scenarios) <- rep("", length(scenarios))
  missing_names <- names(scenarios) == ""
  names(scenarios)[missing_names] <- paste0("scenario_", which(missing_names))
  scenarios
}

scenario_label_2d <- function(name, scenario, do_bootstrap, n_years) {
  if (!is.null(scenario$label)) return(scenario$label)
  model <- if (!is.null(scenario$model)) scenario$model else "gam"
  if (isTRUE(do_bootstrap)) {
    return(sprintf("%s: %s, bootstrap N=%s", name, model, format(n_years, big.mark = ",")))
  }
  sprintf("%s: %s, observed data", name, model)
}

compare_2d_periods <- function(cfg, scenario_results) {
  rp <- cfg$compare_2d_rp
  if (is.null(rp)) return(numeric())

  if (identical(rp, "all")) {
    available <- lapply(scenario_results, function(x) {
      as.numeric(names(x$diform$contours_2d))
    })
    return(sort(unique(unlist(available, use.names = FALSE))))
  }

  as.numeric(rp)
}

prepare_raw_2d_data <- function(csv_path, covariates) {
  prepare_eva_data(fread(csv_path), covariates)
}

use_raw_observed_data_2d <- function(raw_data, n_years = NULL) {
  sim <- copy(raw_data)
  if (is.null(n_years)) n_years <- infer_observation_years(sim)

  if ("time" %in% names(sim)) {
    sim[, sim_year := as.integer(format(as.POSIXct(time), "%Y"))]
  } else {
    sim[, sim_year := NA_integer_]
  }
  sim[, extreme := FALSE]

  list(sim = sim, peak_log = data.table(), n_years = n_years, bootstrap = FALSE)
}

raw_block_resample_2d <- function(raw_data, n_years, seed = 42) {
  set.seed(seed)
  data <- copy(raw_data)
  data[, obs_year := as.integer(format(as.POSIXct(time), "%Y"))]
  available_years <- unique(data$obs_year[!is.na(data$obs_year)])

  sim <- rbindlist(lapply(seq_len(n_years), function(i) {
    x <- data[obs_year == sample(available_years, 1)]
    x[, sim_year := i]
    x
  }))
  sim[, extreme := FALSE]

  list(sim = sim, peak_log = data.table(), n_years = n_years, bootstrap = TRUE)
}

jitter_raw_covariates_2d <- function(sim, raw_data, covariates) {
  out <- copy(sim)
  for (cv in covariates) {
    if (!cv$name %in% names(out)) next
    out[[cv$name]] <- jitter_covariate_values(
      x0 = out[[cv$name]],
      cv = cv,
      data = raw_data,
      peaks = raw_data
    )
  }
  out
}

run_diform_2d_scenario <- function(fit, raw_data, cfg, scenario, scenario_name, tag, multi_scenario) {
  model <- scenario_value_2d(scenario, cfg, "model", "gam")
  do_bootstrap <- !isFALSE(scenario_value_2d(scenario, cfg, "bootstrap", TRUE))

  if (identical(model, "raw")) {
    prep <- NULL
    if (isTRUE(do_bootstrap)) {
      sim_out <- raw_block_resample_2d(
        raw_data = raw_data,
        n_years = scenario_value_2d(scenario, cfg, "N_YEARS"),
        seed = scenario_value_2d(scenario, cfg, "seed", 42)
      )
      if (isTRUE(scenario_value_2d(scenario, cfg, "jitter_covariates", FALSE))) {
        sim_out$sim <- jitter_raw_covariates_2d(
          sim = sim_out$sim,
          raw_data = raw_data,
          covariates = cfg$covariates
        )
      }
    } else {
      sim_out <- use_raw_observed_data_2d(
        raw_data = raw_data,
        n_years = scenario_value_2d(scenario, cfg, "observed_n_years")
      )
    }
  } else if (identical(model, "gam") && do_bootstrap) {
    prep <- run_simulation_prep(fit)
    n_years <- scenario_value_2d(scenario, cfg, "N_YEARS")

    sim_out <- simulate_block_resample(
      fit_obj = fit,
      pred_grids = prep,
      n_years = n_years,
      seed = scenario_value_2d(scenario, cfg, "seed", 42)
    )
  } else if (identical(model, "gam")) {
    prep <- NULL
    sim_out <- use_observed_data(
      fit_obj = fit,
      n_years = scenario_value_2d(scenario, cfg, "observed_n_years")
    )
  }

  analysis_years <- sim_out$n_years
  label <- scenario_label_2d(scenario_name, scenario, do_bootstrap, analysis_years)

  diform <- run_diform_contours_2d(
    sim = sim_out$sim,
    response = cfg$response,
    covariate = cfg$covariate,
    n_years = analysis_years,
    return_periods = cfg$return_periods,
    n_directions = cfg$n_directions,
    h_hours = cfg$h_hours,
    zeta = cfg$zeta,
    n_threads = cfg$n_threads
  )
  diform$bootstrap <- do_bootstrap
  diform$model <- model
  diform$scenario <- scenario_name
  diform$scenario_label <- label
  out_dir <- if (isTRUE(multi_scenario)) {
    file.path(cfg$plot_dir, tag, safe_name_2d(scenario_name))
  } else {
    file.path(cfg$plot_dir, tag)
  }

  individual_plot <- NULL
  if (isTRUE(scenario_value_2d(scenario, cfg, "make_individual_2d_plots", TRUE))) {
    individual_plot <- save_diform_2d_all_r(
      diform_result = diform,
      obs_data = raw_data,
      out_dir = out_dir,
      covariate_label = cfg$covariate_label,
      response_label = cfg$response_label
    )
  }

  top20 <- if ("extreme" %in% names(sim_out$sim) && any(sim_out$sim$extreme)) {
    sim_out$sim[extreme == TRUE][order(-get(cfg$response))]
  } else {
    sim_out$sim[order(-get(cfg$response))]
  }
  top20 <- head(top20, 20)
  if ("threshold" %in% names(top20)) {
    top20[, gpd_excess := get(cfg$response) - threshold]
  }
  top_cols <- c(cfg$response, "threshold", "gpd_excess", "gpd_scale", "gpd_shape", cfg$covariate)
  top_cols <- top_cols[top_cols %in% names(top20)]
  top_print <- top20[, top_cols, with = FALSE]
  setnames(top_print, cfg$response, "value")

  list(
    sim_prep = prep,
    sim = sim_out$sim,
    peak_log = sim_out$peak_log,
    bootstrap = do_bootstrap,
    model = model,
    n_years = analysis_years,
    diform = diform,
    top_events = top_print,
    plots = list(
      ald_splines = if (is.null(fit)) NULL else fit$plots$ald_splines,
      gpd_splines = if (is.null(fit)) NULL else fit$plots$gpd_splines,
      diform_2d = individual_plot
    )
  )
}

run_diform_2d_pipeline <- function(csv_path, tag, cfg) {
  scenarios <- normalise_diform_2d_scenarios(cfg)
  scenario_models <- vapply(scenarios, function(x) {
    scenario_value_2d(x, cfg, "model", "gam")
  }, character(1))

  raw_data <- prepare_raw_2d_data(csv_path, cfg$covariates)
  raw_data <- raw_data[
    is.finite(get(cfg$response)) & is.finite(get(cfg$covariate)) &
      get(cfg$response) > 0 & get(cfg$covariate) > 0
  ]

  fit <- NULL
  if (any(scenario_models == "gam")) {
    fit <- run_evgam_fit(
      csv_path = csv_path,
      tau = cfg$tau,
      decluster_r = cfg$decluster_r,
      response = cfg$response,
      covariates = cfg$covariates,
      ald_interaction = cfg$ald_interaction,
      ald_interaction_k = cfg$ald_interaction_k,
      gpd_interaction = cfg$gpd_interaction,
      gpd_interaction_k = cfg$gpd_interaction_k,
      gpd_shape_covariates = cfg$gpd_shape_covariates,
      gpd_scale_covariates = cfg$gpd_scale_covariates
    )
  }

  multi_scenario <- length(scenarios) > 1 || !is.null(cfg$scenarios)
  scenario_results <- lapply(names(scenarios), function(scenario_name) {
    run_diform_2d_scenario(
      fit = fit,
      raw_data = raw_data,
      cfg = cfg,
      scenario = scenarios[[scenario_name]],
      scenario_name = scenario_name,
      tag = tag,
      multi_scenario = multi_scenario
    )
  })
  names(scenario_results) <- names(scenarios)

  comparison_2d <- NULL
  rp_compare <- compare_2d_periods(cfg, scenario_results)
  if (multi_scenario && length(rp_compare)) {
    diform_results <- lapply(scenario_results, `[[`, "diform")
    comparison_2d <- lapply(rp_compare, function(rp) {
      plot_diform_2d_compare_r(
        diform_results = diform_results,
        obs_data = raw_data,
        out_dir = file.path(cfg$plot_dir, tag),
        rp = rp,
        covariate_label = cfg$covariate_label,
        response_label = cfg$response_label
      )
    })
    names(comparison_2d) <- paste0("T", rp_compare)
  }

  primary <- scenario_results[[1]]
  list(
    tag = tag,
    fit_obj = fit,
    raw_data = raw_data,
    scenarios = scenario_results,
    return_periods = cfg$return_periods,
    plot_dir = cfg$plot_dir,
    response_label = cfg$response_label,
    covariate_label = cfg$covariate_label,
    comparison_2d = comparison_2d,
    sim_prep = primary$sim_prep,
    sim = primary$sim,
    peak_log = primary$peak_log,
    bootstrap = primary$bootstrap,
    n_years = primary$n_years,
    diform = primary$diform,
    top_events = primary$top_events,
    plots = primary$plots
  )
}
