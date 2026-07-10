read_triplet_data <- function(csv_path, n_max = 200000) {
  d <- data.table::fread(csv_path, nrows = n_max)
  as.data.frame(d)
}

triplet_source_value <- function(sources, name) {
  if (is.null(sources)) return(NULL)
  if (!is.null(sources[[name]])) return(sources[[name]])
  NULL
}

prepare_triplet_detection_data <- function(data, variables, sources = NULL) {
  out <- data
  for (nm in variables) {
    src <- triplet_source_value(sources, nm)
    if (!nm %in% names(out) && !is.null(src) && src %in% names(out)) {
      out[[nm]] <- out[[src]]
    }
  }
  out
}

detect_angle_support <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(NULL)
  xr <- range(x)
  span <- diff(xr)
  tol_deg <- 5
  tol_rad <- 0.08

  if (xr[1] >= -tol_deg && xr[2] <= 360 + tol_deg && (xr[2] > 300 || span > 300)) {
    return(list(unit = "degrees", period = 360, origin = 0, support = "0_360"))
  }
  if (xr[1] >= -180 - tol_deg && xr[2] <= 180 + tol_deg &&
      (xr[1] < -120 && xr[2] > 120 || span > 300)) {
    return(list(unit = "degrees", period = 360, origin = -180, support = "minus180_180"))
  }
  if (xr[1] >= -tol_rad && xr[2] <= 2 * pi + tol_rad && (xr[2] > 5.5 || span > 5.5)) {
    return(list(unit = "radians", period = 2 * pi, origin = 0, support = "0_2pi"))
  }
  if (xr[1] >= -pi - tol_rad && xr[2] <= pi + tol_rad &&
      (xr[1] < -2.6 && xr[2] > 2.6 || span > 5.5)) {
    return(list(unit = "radians", period = 2 * pi, origin = -pi, support = "minuspi_pi"))
  }

  NULL
}

detect_directional_variable <- function(data, covariates, direction = "auto") {
  if (identical(direction, "none") || identical(direction, FALSE)) return(NULL)

  if (!identical(direction, "auto")) {
    if (!direction %in% covariates) {
      stop("Forced directional variable is not in the triplet covariates: ", direction,
           call. = FALSE)
    }
    support <- detect_angle_support(data[[direction]])
    if (is.null(support)) {
      support <- list(unit = "degrees", period = 360, origin = 0, support = "forced")
    }
    return(c(list(variable = direction), support))
  }

  detected <- lapply(covariates, function(v) {
    support <- detect_angle_support(data[[v]])
    if (is.null(support)) return(NULL)
    c(list(variable = v), support)
  })
  detected <- Filter(Negate(is.null), detected)

  if (length(detected) > 1) {
    stop(
      "More than one directional covariate was detected: ",
      paste(vapply(detected, `[[`, character(1), "variable"), collapse = ", "),
      ". Force direction = 'none' or direction = '<variable>'.",
      call. = FALSE
    )
  }
  if (!length(detected)) return(NULL)
  detected[[1]]
}

triplet_spec <- function(csv_path,
                         x,
                         y,
                         z,
                         sources = NULL,
                         direction = "auto",
                         data = NULL) {
  if (is.null(data)) data <- read_triplet_data(csv_path)
  covariates <- c(x, y)
  data <- prepare_triplet_detection_data(data, c(covariates, z), sources)
  missing <- setdiff(c(covariates, z), names(data))
  if (length(missing)) {
    stop("Triplet variable(s) missing from data: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  directional <- detect_directional_variable(data, covariates, direction = direction)
  linear_covariates <- setdiff(covariates, if (is.null(directional)) character() else directional$variable)

  list(
    type = "triplet",
    csv_path = csv_path,
    x = x,
    y = y,
    z = z,
    response = z,
    covariates = covariates,
    sources = sources,
    direction = directional,
    linear_covariates = linear_covariates,
    variable_ranges = lapply(c(covariates, z), function(v) range(data[[v]], finite = TRUE))
  )
}

triplet_covariate_spec <- function(triplet, name,
                                   ald_bs = NULL, ald_k = NULL, ald_m = NULL,
                                   gpd_bs = NULL, gpd_k = NULL, gpd_m = NULL,
                                   shape_bs = NULL, shape_k = NULL, shape_m = NULL) {
  is_directional <- !is.null(triplet$direction) && identical(name, triplet$direction$variable)
  source <- triplet_source_value(triplet$sources, name)
  rng <- triplet$variable_ranges[[match(name, c(triplet$covariates, triplet$response))]]
  lower <- if (!is_directional && length(rng) == 2 && is.finite(rng[1]) && rng[1] >= 0) 0 else NULL

  covariate_spec(
    name = name,
    source = source,
    circular = is_directional,
    period = if (is_directional) triplet$direction$period else NULL,
    origin = if (is_directional) triplet$direction$origin else NULL,
    ald_bs = ald_bs,
    ald_k = ald_k,
    ald_m = ald_m,
    gpd_bs = gpd_bs,
    gpd_k = gpd_k,
    gpd_m = gpd_m,
    shape_bs = shape_bs,
    shape_k = shape_k,
    shape_m = shape_m,
    jitter_sd = if (is_directional) {
      if (identical(triplet$direction$unit, "radians")) 10 * pi / 180 else 10
    } else {
      NULL
    },
    lower = lower
  )
}

merge_triplet_covariate_metadata <- function(cv, triplet) {
  base <- triplet_covariate_spec(
    triplet = triplet,
    name = cv$name,
    ald_bs = cv$ald_bs,
    ald_k = cv$ald_k,
    ald_m = cv$ald_m,
    gpd_bs = cv$gpd_bs,
    gpd_k = cv$gpd_k,
    gpd_m = cv$gpd_m,
    shape_bs = cv$shape_bs,
    shape_k = cv$shape_k,
    shape_m = cv$shape_m
  )
  base
}

merge_triplet_candidate_metadata <- function(candidate, triplet) {
  candidate$covariates <- lapply(candidate$covariates, merge_triplet_covariate_metadata,
                                 triplet = triplet)
  candidate
}

triplet_role_defaults <- function(triplet,
                                  numeric_bs = c("cr", "tp"),
                                  directional_bs = c("cc", "cp")) {
  lapply(triplet$covariates, function(v) {
    is_directional <- !is.null(triplet$direction) && identical(v, triplet$direction$variable)
    if (is_directional) {
      list(
        ald_bs = "cp",
        ald_k = 5,
        ald_m = 1,
        gpd_location_bs = directional_bs,
        gpd_location_k = c(10, 15),
        gpd_location_m = 1,
        gpd_shape_bs = directional_bs,
        gpd_shape_k = 5,
        gpd_shape_m = 1
      )
    } else {
      list(
        ald_bs = "tp",
        ald_k = 5,
        ald_m = NULL,
        gpd_location_bs = numeric_bs,
        gpd_location_k = 15,
        gpd_location_m = 1,
        gpd_shape_bs = numeric_bs,
        gpd_shape_k = c(3, 5),
        gpd_shape_m = NULL
      )
    }
  })
}

triplet_default_param <- function(defaults, name) {
  lapply(defaults, `[[`, name)
}

triplet_search_param <- function(search, defaults, name) {
  value <- search[[name]]
  if (is.null(value)) triplet_default_param(defaults, name) else value
}

triplet_reference_basis <- function(triplet, covariate, numeric_bs, directional_bs) {
  is_directional <- !is.null(triplet$direction) && identical(covariate, triplet$direction$variable)
  if (is_directional) directional_bs else numeric_bs
}

shape_covariate_is_used <- function(candidate, covariate) {
  shape_vars <- selection_value(candidate, "gpd_shape_covariates", NULL)
  if (identical(shape_vars, FALSE) || identical(shape_vars, character()) || is.null(shape_vars)) {
    return(FALSE)
  }
  covariate %in% shape_vars
}

candidate_basis_deviations <- function(candidate,
                                       triplet,
                                       numeric_bs = "cr",
                                       directional_bs = "cc") {
  sum(vapply(candidate$covariates, function(cv) {
    ref <- triplet_reference_basis(triplet, cv$name, numeric_bs, directional_bs)
    n <- as.integer(!identical(cv$gpd_bs, ref))
    if (shape_covariate_is_used(candidate, cv$name)) {
      n <- n + as.integer(!identical(cv$shape_bs, ref))
    }
    n
  }, integer(1)))
}

thin_evgam_candidates_by_basis <- function(candidates,
                                           triplet,
                                           max_basis_deviations = 1,
                                           numeric_bs = "cr",
                                           directional_bs = "cc") {
  if (is.null(max_basis_deviations) || !is.finite(max_basis_deviations)) return(candidates)
  keep <- vapply(candidates, function(candidate) {
    candidate_basis_deviations(candidate, triplet, numeric_bs, directional_bs) <= max_basis_deviations
  }, logical(1))
  candidates[keep]
}

evgam_candidate_signature <- function(candidate, response) {
  covariates <- normalise_covariate_specs(candidate$covariates)
  gpd_location <- make_location_formula(
    response = response,
    covariates = covariates,
    prefix = "gpd",
    include_interaction = selection_value(candidate, "gpd_interaction", TRUE),
    interaction_k = selection_value(candidate, "gpd_interaction_k", NULL),
    lhs_override = "excess"
  )
  gpd_shape <- make_scale_formula(
    covariates,
    selection_value(candidate, "gpd_shape_covariates", NULL)
  )
  paste(
    paste(deparse(gpd_location), collapse = " "),
    paste(deparse(gpd_shape), collapse = " "),
    sep = " || "
  )
}

deduplicate_evgam_candidates <- function(candidates, response) {
  signatures <- vapply(candidates, evgam_candidate_signature, character(1), response = response)
  candidates[!duplicated(signatures)]
}

evgam_search_spec <- function(gpd_interaction = c(FALSE, TRUE),
                              gpd_interaction_k = list(c(5, 5)),
                              gpd_shape_covariates = "all",
                              numeric_bs = c("cr", "tp"),
                              directional_bs = c("cc", "cp"),
                              gpd_location_bs = NULL,
                              gpd_location_k = NULL,
                              gpd_shape_bs = NULL,
                              basis_grid = c("sparse", "full"),
                              reference_numeric_bs = "cr",
                              reference_directional_bs = "cc",
                              max_basis_deviations = NULL,
                              deduplicate = TRUE,
                              plot_best = TRUE,
                              verbose = TRUE) {
  basis_grid <- match.arg(basis_grid)
  if (is.null(max_basis_deviations)) {
    max_basis_deviations <- if (identical(basis_grid, "sparse")) 1 else Inf
  }

  list(
    gpd_interaction = gpd_interaction,
    gpd_interaction_k = gpd_interaction_k,
    gpd_shape_covariates = gpd_shape_covariates,
    numeric_bs = numeric_bs,
    directional_bs = directional_bs,
    gpd_location_bs = gpd_location_bs,
    gpd_location_k = gpd_location_k,
    gpd_shape_bs = gpd_shape_bs,
    basis_grid = basis_grid,
    reference_numeric_bs = reference_numeric_bs,
    reference_directional_bs = reference_directional_bs,
    max_basis_deviations = max_basis_deviations,
    deduplicate = deduplicate,
    plot_best = plot_best,
    verbose = verbose
  )
}

make_evgam_candidates_for_triplet <- function(triplet, search = evgam_search_spec()) {
  defaults <- triplet_role_defaults(
    triplet,
    numeric_bs = search$numeric_bs,
    directional_bs = search$directional_bs
  )
  candidates <- make_evgam_candidates(
    covariate1 = triplet$covariates[1],
    covariate2 = triplet$covariates[2],
    ald_bs = unlist(triplet_default_param(defaults, "ald_bs"), use.names = FALSE),
    ald_k = unlist(triplet_default_param(defaults, "ald_k"), use.names = FALSE),
    ald_m = triplet_default_param(defaults, "ald_m"),
    gpd_location_bs = triplet_search_param(search, defaults, "gpd_location_bs"),
    gpd_location_k = triplet_search_param(search, defaults, "gpd_location_k"),
    gpd_location_m = triplet_default_param(defaults, "gpd_location_m"),
    gpd_shape_bs = triplet_search_param(search, defaults, "gpd_shape_bs"),
    gpd_shape_k = triplet_default_param(defaults, "gpd_shape_k"),
    gpd_shape_m = triplet_default_param(defaults, "gpd_shape_m"),
    gpd_interaction = search$gpd_interaction,
    gpd_interaction_k = search$gpd_interaction_k,
    gpd_shape_covariates = search$gpd_shape_covariates
  )
  candidates <- lapply(candidates, merge_triplet_candidate_metadata, triplet = triplet)
  candidates <- thin_evgam_candidates_by_basis(
    candidates,
    triplet = triplet,
    max_basis_deviations = search$max_basis_deviations,
    numeric_bs = search$reference_numeric_bs,
    directional_bs = search$reference_directional_bs
  )
  if (isTRUE(search$deduplicate)) {
    candidates <- deduplicate_evgam_candidates(candidates, response = triplet$response)
  }
  candidates
}

select_best_evgam_for_triplet <- function(triplet,
                                          tau = 0.96,
                                          decluster_r = 72,
                                          search = evgam_search_spec(),
                                          ald_interaction = FALSE,
                                          ald_interaction_k = NULL) {
  candidates <- make_evgam_candidates_for_triplet(triplet, search)
  select_evgam_by_qq(
    csv_path = triplet$csv_path,
    candidates = candidates,
    response = triplet$response,
    tau = tau,
    decluster_r = decluster_r,
    ald_covariates = candidates[[1]]$covariates,
    ald_interaction = ald_interaction,
    ald_interaction_k = ald_interaction_k,
    plot_best = isTRUE(search$plot_best),
    verbose = isTRUE(search$verbose)
  )
}

contour_from_triplet <- function(triplet) {
  if (is.null(triplet$direction)) {
    return(linear_diform_spec(
      variables = c(triplet$covariates, triplet$response),
      response = triplet$response,
      response_component = triplet$response
    ))
  }

  if (length(triplet$linear_covariates) != 1) {
    stop("A 3D directional D-IFORM triplet needs exactly one non-directional covariate.",
         call. = FALSE)
  }

  directional_diform_spec(
    linear = triplet$linear_covariates[1],
    angle = triplet$direction$variable,
    response = triplet$response,
    angle_unit = triplet$direction$unit,
    period = triplet$direction$period,
    origin = triplet$direction$origin
  )
}
