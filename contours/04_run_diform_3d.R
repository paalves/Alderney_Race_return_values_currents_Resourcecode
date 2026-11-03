response_variable <- "Os"

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

root <- project_file("contours")
source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
source(file.path(root, "01_diform_contours_3d.R"))
source(file.path(root, "02_plot_diform_3d.R"))
source(file.path(root, "03_plot_diform_diagnostics.R"))

dataset_tag <- "rscd_2026"
response_tag <- tolower(response_variable)
response_info <- list(
  hs = list(label = "Hs [m]", short = "Hs", unit = "m", components = c("H_L", "H_T")),
  os = list(label = "Os [m s\u207b\u00b9]", short = "Os", unit = "m s\u207b\u00b9", components = c("O_L", "O_T"))
)[[response_tag]]
data <- read_rscd_2026()
data <- data[
  is.finite(Cspd) & is.finite(angleDiff) & is.finite(Os) &
    Cspd >= 0 & Os > 0
]

analysis_years <- observation_years(data)
if (!is.finite(analysis_years) || analysis_years <= 0) {
  stop("The observation period could not be determined from data_rscd_alderney_2026.csv.",
       call. = FALSE)
}

diform_spec <- directional_diform_spec(
  linear = "Cspd",
  angle = "angleDiff",
  response = response_variable,
  angle_unit = "degrees",
  period = 360,
  origin = 0
)
diform_spec$variable_labels <- list(
  Cspd = "Current speed",
  angleDiff = "Angular difference"
)
diform_spec$variable_units <- list(
  Cspd = "m s\u207b\u00b9",
  angleDiff = "\u00b0"
)
result_3d <- run_diform_contours(
  data = data,
  response = response_variable,
  n_years = analysis_years,
  return_periods = c(1, 2, 5, 10, 20, 50, 100),
  simplex_resolution = 20,
  h_hours = 72,
  zeta = 0.05,
  n_threads = 0,
  diform_spec = diform_spec
)

result_3d$data_label <- diform3d_data_label(result_3d)

plot_dir <- figure_dir("contours", dataset_tag, response_tag)
plot3d_files <- plot_diform_3d_r(
  diform_result = result_3d,
  obs_data = data,
  response_label = response_info$label,
  response_short = response_info$short,
  response_unit = response_info$unit,
  out_dir = plot_dir,
  rp = "all",
  elev = 45,
  azim = 40,
  show_points = TRUE,
  show_peak_highlights = FALSE,
  show_support_plane = FALSE,
  max_background_points = 20000,
  z_max = NULL,
  working_z_max = NULL,
  working_component_labels = response_info$components,
  make_halfspace_plot = FALSE,
  make_working_plot = TRUE
)

result_3d$plot3d_files <- plot3d_files

gpd_qq_file <- file.path(plot_dir, "diform_gpd_qq_diagnostics.png")
gpd_qq_plot <- plot_diform_gpd_qq(
  diform_result = result_3d,
  n_random = 3,
  seed = 2023,
  level = 0.95,
  out_file = gpd_qq_file
)

result_3d$gpd_qq_file <- gpd_qq_file

result_3d
