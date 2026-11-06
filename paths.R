#' Create and return a figure directory.
figure_dir <- function(...) {
  path <- file.path("figures", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
