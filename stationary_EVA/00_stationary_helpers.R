#' Label, unit, display precision and values of the response (hs or Os).
stationary_response_info <- function(data, response_variable) {
  info <- if (response_variable == "hs") {
    list(label = "Hs", unit = "m", digits = 2L)
  } else {
    list(label = "Os", unit = "m s⁻¹", digits = 3L)
  }
  c(list(name = response_variable, tag = tolower(response_variable), values = data[[response_variable]]), info)
}

#' Axis label with the unit in parentheses.
stationary_axis_label <- function(prefix, response_info) {
  sprintf("%s (%s)", prefix, response_info$unit)
}

#' Label giving the threshold value and its unit.
stationary_threshold_label <- function(value, response_info, prefix = "Q96") {
  paste0(prefix, " = ", formatC(value, format = "f", digits = response_info$digits), " ", response_info$unit)
}

#' Output path of a stationary-analysis figure for the response.
stationary_figure_file <- function(response_info, filename) {
  file.path(figure_dir("stationary_EVA", response_info$tag), filename)
}

#' Fit a stationary GPD above the threshold by maximum likelihood.
fit_stationary_gpd <- function(x, threshold) {
  evd::fpot(x, threshold = threshold, model = "gpd")
}

#' Runs declustering: maxima of the clusters of exceedances separated by more than run_length steps.
decluster_maxima <- function(x, threshold, run_length) {
  exceedance_index <- which(x > threshold)
  clusters <- split(exceedance_index, cumsum(c(TRUE, diff(exceedance_index) > run_length)))
  peak_index <- vapply(clusters, function(index) index[which.max(x[index])], integer(1))
  data.frame(index = peak_index, value = x[peak_index], excess = x[peak_index] - threshold)
}
