.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

root <- project_file("nonstationary_EVA")

source(file.path(root, "01_evgam_fit.R"))
source(file.path(root, "02_simulation.R"))
source(file.path(root, "03b_diform_contours_2d.R"))
source(file.path(root, "05_pipeline_2d_helpers.R"))

#be careful, this code has heavier computational cost compared to the 3D pipeline because of the 3000-year simulation. 
#this code is for making figure 17

dataset_tag <- "Raz_Blanchard"
csv_path <- data_file()

covariates_hs_tp <- list(
  covariate_spec(
    name = "tp",
    source = "Tp",
    circular = FALSE,
    ald_bs = "tp", ald_k = 6,
    gpd_bs = "cr", gpd_k = 12, gpd_m = 1,
    lower = 0
  )
)

cfg <- list(
  response = "hs",
  covariate = "tp",
  tau = 0.96,
  decluster_r = 72,
  covariates = covariates_hs_tp,
  ald_interaction = FALSE,
  ald_interaction_k = NULL,
  gpd_interaction = FALSE,
  gpd_interaction_k = NULL,
  gpd_scale_covariates = NULL,
  N_YEARS = 100,
  bootstrap = FALSE,
  observed_n_years = NULL,
  scenarios = list(
    raw_observed = list(
      model = "raw",
      bootstrap = FALSE,
      label = "no GAM, observed data"
    ),
    raw_boot_1000 = list(
      model = "raw",
      bootstrap = TRUE,
      N_YEARS = 1000,
      jitter_covariates = TRUE,
      label = "no GAM, 1000-yr resampling"
    ),
    gam_observed = list(
      model = "gam",
      bootstrap = FALSE,
      label = "GAM, observed data"
    ),
    gam_boot_1000 = list(
      model = "gam",
      bootstrap = TRUE,
      N_YEARS = 1000,
      label = "GAM, 1000-yr resampling"
    ),
    gam_boot_3000 = list(
      model = "gam",
      bootstrap = TRUE,
      N_YEARS = 3000,
      label = "GAM, 3000-yr resampling"
    )
  ),
  seed = 70,
  return_periods = c(1, 2, 5, 10, 20, 50, 100),
  n_directions = 30,
  h_hours = 72,
  zeta = 0.05,
  n_threads = 0,
  plot_dir = file.path(root, "plots_2d"),
  covariate_label = "Tp [s]",
  response_label = "Hs [m]",
  make_individual_2d_plots = FALSE,
  compare_2d_rp = "all"
)

result_2d <- run_diform_2d_pipeline(
  csv_path = csv_path,
  tag = dataset_tag,
  cfg = cfg
)
