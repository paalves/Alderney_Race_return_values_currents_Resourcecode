.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

root <- project_file("nonstationary_EVA")

source(file.path(root, "01_evgam_fit.R"))
source(file.path(root, "00_select_evgam_model.R"))
source(file.path(root, "02_simulation.R"))
source(file.path(root, "03_diform_contours.R"))
source(file.path(root, "plot_diform_3d.R"))
source(file.path(root, "04b_pipeline_triplet.R"))
source(file.path(root, "04_pipeline_helpers.R"))

#figures 18 and 19. to switch from Os to Hs, change args in "triplet_spec" and "run_triplet_diform_pipeline" functions below. also in plot_3d_spec further down, change z_max argument<
#to adapt figure scale. use elev/azim args to get a different view of the 3d figures. 

dataset_tag <- "Alderney_Race"
csv_path <- data_file()
selection_cache <- project_file("outputs", paste0("evgam_selection_", dataset_tag, ".rds"))

RUN_EVGAM_SELECTION <- TRUE #you can set it to FALSE after running it once. re-enable it if you want to change the covariates or the response. 
RUN_DIFORM_PIPELINE <- TRUE

triplet <- triplet_spec(
  csv_path = csv_path,
  x = "Cspd",
  y = "angleDiff",
  z = "Hs",
  sources = list(angularDifference = "Angular difference"),
  direction = "auto"
)

if (isTRUE(RUN_EVGAM_SELECTION)) {
  selection <- select_best_evgam_for_triplet(
    triplet = triplet,
    tau = 0.96,
    decluster_r = 72,
    search = evgam_search_spec(
      gpd_interaction = c(FALSE, TRUE),
      gpd_interaction_k = list(c(5, 5)),
      gpd_shape_covariates = "all",
      numeric_bs = c("cr", "tp"),
      directional_bs = c("cc", "cp"),
      basis_grid = "sparse",
      plot_best = TRUE,
      verbose = TRUE
    )
  )
  dir.create(dirname(selection_cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(selection, selection_cache)
} else {
  if (exists("selection", inherits = TRUE)) {
    selection <- get("selection", inherits = TRUE)
  } else {
    selection <- readRDS(selection_cache)
  }
}

best_evgam_fit <- selection$best_fit

scenarios <- list(
  boot_1000 = bootstrap_simulation(n_years = 1000, seed = 70)
)

if (isTRUE(RUN_DIFORM_PIPELINE)) {
  runs <- run_triplet_diform_pipeline(
    triplet = triplet,
    fit = best_evgam_fit,
    scenarios = scenarios,
    return_periods = c(1, 2, 5, 10, 20, 50, 100),
    tag = dataset_tag,
    plot_dir = file.path(root, "plots"),
    response_label = "Hs [m]",
    response_short = "Hs",
    response_unit = "m",
    cov1_unit = "m.s^-1",
    cov2_unit = "°",
    cov1_label = "Current speed",
    cov2_label = "Angular difference",
    response_component_labels = c("HL", "HT"),
    simplex_resolution = 20, #produces 1600 directions. computationally intensive. you can lower this to 10 reduces it to 400 directions.
    h_hours = 72,
    zeta = 0.05,
    n_threads = 0,
    top_events = 20
  )

  projection_plots <- plot_diform_projection_runs(runs)

  plot3d_files <- plot_diform_3d_runs(
    runs,
    plot3d_spec(
      rp = "all",
      elev = 45,
      azim = 40,
      observed_points = TRUE,
      peak_highlights = FALSE,
      tangent_plane = FALSE,
      max_background_points = 0,
      z_max = 8,
      working_z_max = NULL,
      make_halfspace_plot = FALSE
    )
  )

  result <- c(
    runs,
    list(
      selection = selection,
      projection_plots = projection_plots,
      plot3d_files = plot3d_files
    )
  )
} else {
  result <- list(
    triplet = triplet,
    selection = selection,
    best_evgam_fit = best_evgam_fit
  )
}
