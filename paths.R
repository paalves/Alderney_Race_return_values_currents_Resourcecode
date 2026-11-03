#' Project root directory, from PAPER_CODE_ROOT, the location of paths.R or the working directory.
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
  while (basename(wd) %in% c("stationary_EVA", "site_description", "nonstationary_EVA", "contours", "importation")) {
    wd <- dirname(wd)
  }
  wd
}

#' Path to a file inside the project.
project_file <- function(...) {
  file.path(project_root(), ...)
}

#' Path to the hindcast data file, overridable by DATA_RSCD_2026_PATH.
data_file <- function(filename = "data_rscd_alderney_2026.csv") {
  env_path <- Sys.getenv("DATA_RSCD_2026_PATH", unset = "")
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
      ". Set DATA_RSCD_2026_PATH or put the file in data/.",
      call. = FALSE
    )
  }

  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

#' Create and return a figure directory inside the project.
figure_dir <- function(...) {
  path <- project_file("figures", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' Find and source paths.R from the working directory or its parents.
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
