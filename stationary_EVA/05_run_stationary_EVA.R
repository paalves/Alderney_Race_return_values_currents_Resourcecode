response_variable <- "Os"

source("paths.R")
source("importation/02_read_rscd_2026.R")
source("stationary_EVA/00_stationary_helpers.R")

data <- read_rscd_2026()
response_info <- stationary_response_info(data, response_variable)
response_values <- response_info$values
chosen_threshold <- as.numeric(quantile(response_values, 0.96))

for (script in c("01_threshold_selection.R", "02_declustering_sensitivity.R", "03_analysis.R", "04_EVA_plots.R")) {
  source(file.path("stationary_EVA", script), local = environment())
}

stationary_eva_result <- list(
  response = response_info[c("name", "tag", "label", "unit")],
  data_summary = data_summary,
  threshold_selection = list(
    chosen = threshold_summary,
    alternative = alternative_threshold,
    diagnostics = param_data
  ),
  declustering_sensitivity = cluster_results,
  model = stationary_analysis,
  figures = list(
    threshold_selection = threshold_figure_file,
    declustering_sensitivity = declustering_figure_file,
    gpd_model_fit = gpd_figure_file
  )
)

dir.create("outputs", showWarnings = FALSE)
saveRDS(stationary_eva_result, file.path("outputs", paste0("stationary_EVA_", response_info$tag, ".rds")))

print(threshold_summary)
print(fit_summary)
print(return_values)
