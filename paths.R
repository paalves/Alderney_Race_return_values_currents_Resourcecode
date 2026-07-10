project_root <- function() {
  root <- Sys.getenv("PAPER_CODE_ROOT", unset = "")
  if (nzchar(root)) {
    return(normalizePath(root, winslash = "/", mustWork = TRUE))
  }

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  path_file <- ofiles[basename(ofiles) == "paths.R"]

  if (length(path_file)) {
    return(dirname(normalizePath(tail(path_file, 1), winslash = "/", mustWork = TRUE)))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (basename(wd) %in% c("stationary_EVA", "site_description", "nonstationary_EVA")) {
    wd <- dirname(wd)
  }
  wd
}

project_file <- function(...) {
  file.path(project_root(), ...)
}

data_file <- function(filename = "data_alderney_rscd_2024.csv") {
  env_path <- Sys.getenv("DATA_ALDERNEY_RSCD_2024_PATH", unset = "")
  if (nzchar(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = TRUE))
  }

  candidates <- c(
    project_file("data", filename),
    project_file(filename),
    file.path(dirname(project_root()), filename)
  )
  hit <- candidates[file.exists(candidates)][1]

  if (is.na(hit)) {
    stop(
      "Could not find ", filename,
      ". Set DATA_ALDERNEY_RSCD_2024_PATH or put the file in data/.",
      call. = FALSE
    )
  }

  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

figure_dir <- function(...) {
  path <- project_file("figures", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

source_project_paths <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(start, dirname(start), dirname(dirname(start))))
  helper <- file.path(candidates, "paths.R")
  helper <- helper[file.exists(helper)][1]

  if (is.na(helper)) {
    stop("Could not find paths.R from: ", start, call. = FALSE)
  }

  source(helper)
}
