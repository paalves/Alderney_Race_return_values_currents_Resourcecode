safe_name <- function(x) {
  gsub("_+$", "", gsub("^_+", "", gsub("[^A-Za-z0-9_]+", "_", x)))
}

observed_simulation <- function(observed_n_years = NULL, label = NULL) {
  list(bootstrap = FALSE, observed_n_years = observed_n_years, label = label)
}

bootstrap_simulation <- function(n_years, seed = 42, label = NULL) {
  list(bootstrap = TRUE, n_years = n_years, N_YEARS = n_years, seed = seed, label = label)
}

scenario_value <- function(scenario, name, default = NULL) {
  if (!is.null(scenario[[name]])) scenario[[name]] else default
}

normalise_triplet_scenarios <- function(scenarios) {
  if (is.null(scenarios)) scenarios <- list(default = bootstrap_simulation(100))
  if (is.null(names(scenarios))) names(scenarios) <- rep("", length(scenarios))
  missing_names <- names(scenarios) == ""
  names(scenarios)[missing_names] <- paste0("scenario_", which(missing_names))
  scenarios
}

scenario_label <- function(name, scenario, do_bootstrap, n_years) {
  if (!is.null(scenario$label)) return(scenario$label)
  if (isTRUE(do_bootstrap)) {
    return(sprintf("%s: bootstrap, N=%s", name, format(n_years, big.mark = ",")))
  }
  sprintf("%s: observed data", name)
}

triplet_covariate_units <- function(triplet, cov1_unit = NULL, cov2_unit = NULL) {
  units <- list(cov1_unit, cov2_unit)
  names(units) <- triplet$covariates
  units <- units[!vapply(units, function(x) {
    is.null(x) || !length(x) || anyNA(x) || !nzchar(x)
  }, logical(1))]
  units
}

triplet_covariate_labels <- function(triplet, cov1_label = NULL, cov2_label = NULL) {
  labels <- list(cov1_label, cov2_label)
  names(labels) <- triplet$covariates
  labels <- labels[!vapply(labels, function(x) {
    is.null(x) || !length(x) || anyNA(x) || !nzchar(x)
  }, logical(1))]
  labels
}

add_triplet_units_to_diform_spec <- function(diform_spec, triplet,
                                             cov1_unit = NULL, cov2_unit = NULL) {
  units <- triplet_covariate_units(triplet, cov1_unit, cov2_unit)
  if (!length(units)) return(diform_spec)
  diform_spec$variable_units <- units
  diform_spec
}

add_triplet_labels_to_diform_spec <- function(diform_spec, triplet,
                                              cov1_label = NULL, cov2_label = NULL) {
  labels <- triplet_covariate_labels(triplet, cov1_label, cov2_label)
  if (!length(labels)) return(diform_spec)
  if (is.null(diform_spec$variable_labels)) diform_spec$variable_labels <- list()
  diform_spec$variable_labels[names(labels)] <- labels
  diform_spec
}

diform_response_component_labels <- function(response_short, response_component_labels = NULL) {
  if (!is.null(response_component_labels)) return(response_component_labels)
  prefix <- gsub("[^A-Za-z]", "", response_short)
  prefix <- if (nzchar(prefix)) toupper(substr(prefix, 1, 1)) else "R"
  paste0(prefix, c("L", "T"))
}

add_response_component_labels_to_diform_spec <- function(diform_spec,
                                                         response_short,
                                                         response_component_labels = NULL) {
  if (is.null(diform_spec$components)) return(diform_spec)
  is_circular <- vapply(diform_spec$components, function(comp) {
    identical(comp$type, "circular_radius")
  }, logical(1))
  if (!any(is_circular)) return(diform_spec)

  circular <- diform_spec$components[[which(is_circular)[1]]]
  if (is.null(circular$names) || length(circular$names) != 2) return(diform_spec)

  labels <- diform_response_component_labels(response_short, response_component_labels)
  if (is.null(diform_spec$variable_labels)) diform_spec$variable_labels <- list()
  diform_spec$variable_labels[circular$names] <- as.list(labels)
  diform_spec
}

plot3d_spec <- function(make = TRUE,
                        rp = "all",
                        elev = 50,
                        azim = 50,
                        n_ang = 180,
                        n_cspd = 80,
                        observed_points = TRUE,
                        peak_highlights = FALSE,
                        tangent_plane = FALSE,
                        exterior_halfspace = FALSE,
                        max_background_points = 0,
                        max_labels = 6,
                        plane_scale = 0.28,
                        z_max = NULL,
                        working_z_max = NULL,
                        make_halfspace_plot = FALSE) {
  list(
    make = make,
    rp = rp,
    elev = elev,
    azim = azim,
    n_ang = n_ang,
    n_cspd = n_cspd,
    observed_points = observed_points,
    peak_highlights = peak_highlights,
    tangent_plane = tangent_plane,
    exterior_halfspace = exterior_halfspace,
    max_background_points = max_background_points,
    max_labels = max_labels,
    plane_scale = plane_scale,
    z_max = z_max,
    working_z_max = working_z_max,
    make_halfspace_plot = make_halfspace_plot
  )
}

run_triplet_diform_scenario <- function(fit,
                                        triplet,
                                        diform_spec,
                                        scenario,
                                        scenario_name,
                                        return_periods,
                                        simplex_resolution = 20,
                                        h_hours = 72,
                                        zeta = 0.05,
                                        n_threads = 0,
                                        top_events = 20) {
  do_bootstrap <- !isFALSE(scenario_value(scenario, "bootstrap", TRUE))

  if (do_bootstrap) {
    prep <- run_simulation_prep(fit)
    n_years <- scenario_value(scenario, "n_years", scenario_value(scenario, "N_YEARS"))
    sim_out <- simulate_block_resample(
      fit_obj = fit,
      pred_grids = prep,
      n_years = n_years,
      seed = scenario_value(scenario, "seed", 42)
    )
  } else {
    prep <- NULL
    sim_out <- use_observed_data(
      fit_obj = fit,
      n_years = scenario_value(scenario, "observed_n_years")
    )
  }

  analysis_years <- sim_out$n_years
  label <- scenario_label(scenario_name, scenario, do_bootstrap, analysis_years)

  diform <- run_diform_contours(
    sim = sim_out$sim,
    response = triplet$response,
    n_years = analysis_years,
    return_periods = return_periods,
    simplex_resolution = simplex_resolution,
    h_hours = h_hours,
    zeta = zeta,
    n_threads = n_threads,
    diform_spec = diform_spec
  )
  diform$bootstrap <- do_bootstrap
  diform$scenario <- scenario_name
  diform$scenario_label <- label

  top_events_table <- NULL
  if (is.finite(top_events) && top_events > 0 && "extreme" %in% names(sim_out$sim)) {
    top <- sim_out$sim[extreme == TRUE][order(-get(triplet$response))]
    top <- head(top, top_events)
    if (nrow(top)) {
      top[, gpd_excess := get(triplet$response) - threshold]
      top_cols <- unique(c(
        triplet$response,
        "threshold",
        "gpd_excess",
        "gpd_scale",
        "gpd_shape",
        triplet$covariates
      ))
      top_cols <- top_cols[top_cols %in% names(top)]
      top_print <- top[, top_cols, with = FALSE]
      setnames(top_print, triplet$response, "value")
      top_events_table <- top_print
    }
  }

  list(
    sim_prep = prep,
    sim = sim_out$sim,
    peak_log = sim_out$peak_log,
    bootstrap = do_bootstrap,
    n_years = analysis_years,
    diform = diform,
    top_events = top_events_table
  )
}

run_triplet_diform_pipeline <- function(triplet,
                                        fit,
                                        scenarios,
                                        return_periods,
                                        tag = "triplet",
                                        plot_dir = file.path(root, "plots"),
                                        response_label = paste0(triplet$response, " [m]"),
                                        response_short = triplet$response,
                                        response_unit = "",
                                        cov1_unit = NULL,
                                        cov2_unit = NULL,
                                        cov1_label = NULL,
                                        cov2_label = NULL,
                                        response_component_labels = NULL,
                                        simplex_resolution = 20,
                                        h_hours = 72,
                                        zeta = 0.05,
                                        n_threads = 0,
                                        top_events = 20) {
  scenarios <- normalise_triplet_scenarios(scenarios)
  diform_spec <- contour_from_triplet(triplet)
  diform_spec <- add_triplet_units_to_diform_spec(
    diform_spec,
    triplet,
    cov1_unit = cov1_unit,
    cov2_unit = cov2_unit
  )
  diform_spec <- add_triplet_labels_to_diform_spec(
    diform_spec,
    triplet,
    cov1_label = cov1_label,
    cov2_label = cov2_label
  )
  diform_spec <- add_response_component_labels_to_diform_spec(
    diform_spec,
    response_short = response_short,
    response_component_labels = response_component_labels
  )
  scenario_results <- lapply(names(scenarios), function(scenario_name) {
    run_triplet_diform_scenario(
      fit = fit,
      triplet = triplet,
      diform_spec = diform_spec,
      scenario = scenarios[[scenario_name]],
      scenario_name = scenario_name,
      return_periods = return_periods,
      simplex_resolution = simplex_resolution,
      h_hours = h_hours,
      zeta = zeta,
      n_threads = n_threads,
      top_events = top_events
    )
  })
  names(scenario_results) <- names(scenarios)

  primary <- scenario_results[[1]]
  list(
    tag = tag,
    triplet = triplet,
    fit_obj = fit,
    diform_spec = diform_spec,
    scenarios = scenario_results,
    return_periods = return_periods,
    plot_dir = plot_dir,
    response_label = response_label,
    response_short = response_short,
    response_unit = response_unit,
    covariate_units = triplet_covariate_units(triplet, cov1_unit, cov2_unit),
    covariate_labels = triplet_covariate_labels(triplet, cov1_label, cov2_label),
    response_component_labels = diform_response_component_labels(response_short, response_component_labels),
    sim_prep = primary$sim_prep,
    sim = primary$sim,
    peak_log = primary$peak_log,
    bootstrap = primary$bootstrap,
    n_years = primary$n_years,
    diform = primary$diform,
    top_events = primary$top_events
  )
}

plot_diform_projection_runs <- function(runs) {
  plots <- lapply(runs$scenarios, function(res) {
    plot_diform_all(
      diform_result = res$diform,
      obs_data = runs$fit_obj$data,
      response = runs$triplet$response,
      resp_label = runs$response_label
    )
  })
  invisible(plots)
}

plot_diform_3d_runs <- function(runs, plot_spec = plot3d_spec()) {
  if (!isTRUE(plot_spec$make)) return(invisible(NULL))
  multi_scenario <- length(runs$scenarios) > 1

  files <- lapply(names(runs$scenarios), function(scenario_name) {
    res <- runs$scenarios[[scenario_name]]
    out_dir <- if (multi_scenario) {
      file.path(runs$plot_dir, runs$tag, safe_name(scenario_name))
    } else {
      file.path(runs$plot_dir, runs$tag)
    }

    plot_diform_3d_r(
      diform_result = res$diform,
      obs_data = runs$fit_obj$data,
      excess_df = runs$fit_obj$excess_df,
      response_label = runs$response_label,
      response_short = runs$response_short,
      response_unit = runs$response_unit,
      out_dir = out_dir,
      rp = plot_spec$rp,
      elev = plot_spec$elev,
      azim = plot_spec$azim,
      n_ang = plot_spec$n_ang,
      n_cspd = plot_spec$n_cspd,
      show_points = isTRUE(plot_spec$observed_points),
      show_peak_highlights = isTRUE(plot_spec$peak_highlights),
      show_support_plane = isTRUE(plot_spec$tangent_plane),
      show_exterior_halfspace = isTRUE(plot_spec$exterior_halfspace),
      plane_scale = plot_spec$plane_scale,
      max_background_points = plot_spec$max_background_points,
      max_labels = plot_spec$max_labels,
      z_max = plot_spec$z_max,
      working_z_max = plot_spec$working_z_max,
      make_halfspace_plot = isTRUE(plot_spec$make_halfspace_plot)
    )
  })
  names(files) <- names(runs$scenarios)
  invisible(files)
}
