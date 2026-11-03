.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
source(project_file("nonstationary_EVA", "01_select_evgam.R"))

dataset_tag <- "rscd_2026"
data <- read_rscd_2026(data_file("data_rscd_alderney_2026.csv"))
RUN_MODEL_SELECTION <- TRUE

result_evgam <- lapply(c(hs = "hs", os = "Os"), function(response) {
  cache <- project_file("outputs", sprintf("evgam_selection_%s_%s.rds", dataset_tag, tolower(response)))
  if (isTRUE(RUN_MODEL_SELECTION)) {
    selection <- select_evgam(data, response)
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    saveRDS(selection, cache)
  } else {
    selection <- readRDS(cache)
  }
  message("Cross-validated log score, ", response, ":")
  print(selection$summary)
  selection
})
