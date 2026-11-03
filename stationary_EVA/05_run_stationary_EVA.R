response_variable <- "Os"  

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

stationary_root <- project_file("stationary_EVA")
stationary_scripts <- c(
  "01_threshold_selection.R",
  "02_declustering_sensitivity.R",
  "03_analysis.R",
  "04_EVA_plots.R"
)

for (script in stationary_scripts) {
  message("Running ", script, " for response ", response_variable, "...")
  source(file.path(stationary_root, script), local = environment())
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

result_file <- project_file(
  "outputs",
  paste0("stationary_EVA_", response_info$tag, ".rds")
)
dir.create(dirname(result_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(stationary_eva_result, result_file)

message("Stationary EVA complete for ", response_info$label, ".")
message("Results: ", result_file)
message("Figures: ", dirname(gpd_figure_file))

print(threshold_summary)
print(fit_summary)
print(return_values)

invisible(stationary_eva_result)
