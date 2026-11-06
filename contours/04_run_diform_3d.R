response_variable <- "Os"

source("paths.R")
source("importation/02_read_rscd_2026.R")
source("contours/01_diform_contours_3d.R")
source("contours/02_plot_diform_3d.R")
source("contours/03_plot_diform_diagnostics.R")

response_tag <- tolower(response_variable)
response_info <- list(
  hs = list(label = "Hs [m]", short = "Hs", unit = "m", components = c("H_L", "H_T")),
  os = list(label = "Os [m s⁻¹]", short = "Os", unit = "m s⁻¹", components = c("O_L", "O_T"))
)[[response_tag]]

data <- read_rscd_2026()
result_3d <- run_diform_contours(data, response_variable, observation_years(data))

plot_dir <- figure_dir("contours", "rscd_2026", response_tag)
plot_diform_3d_r(result_3d, data, response_info, plot_dir, show_support_plane = FALSE)
plot_diform_gpd_qq(result_3d, file.path(plot_dir, "diform_gpd_qq_diagnostics.png"))
