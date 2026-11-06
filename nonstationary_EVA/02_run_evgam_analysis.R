source("importation/02_read_rscd_2026.R")
source("nonstationary_EVA/01_select_evgam.R")

dataset_tag <- "rscd_2026"
data <- read_rscd_2026()
RUN_MODEL_SELECTION <- TRUE

result_evgam <- lapply(c(hs = "hs", os = "Os"), function(response) {
  cache <- file.path("outputs", sprintf("evgam_selection_%s_%s.rds", dataset_tag, tolower(response)))
  if (RUN_MODEL_SELECTION) {
    selection <- select_evgam(data, response)
    dir.create("outputs", showWarnings = FALSE)
    saveRDS(selection, cache)
  } else {
    selection <- readRDS(cache)
  }
  message("Cross-validated log score, ", response, ":")
  print(selection$summary)
  selection
})
