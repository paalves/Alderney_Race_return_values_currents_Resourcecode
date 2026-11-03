#' Find the response column and its label, unit and display precision.
stationary_response_info <- function(data, response_variable = "hs") {
  response_variable <- trimws(as.character(response_variable)[1])
  column_index <- match(tolower(response_variable), tolower(names(data)))

  if (is.na(column_index)) {
    stop(
      "Unknown stationary-EVA response '", response_variable,
      "'. Available columns include: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  response_name <- names(data)[column_index]
  response_key <- tolower(response_name)
  metadata <- switch(
    response_key,
    hs = list(label = "Hs", unit = "m", digits = 2L),
    os = list(label = "Os", unit = "m s\u207b\u00b9", digits = 3L),
    list(label = response_name, unit = "", digits = 3L)
  )

  values <- as.numeric(data[[response_name]])
  values <- values[is.finite(values)]
  if (length(values) < 100L) {
    stop("Too few finite observations for response '", response_name, "'.", call. = FALSE)
  }

  c(
    list(
      name = response_name,
      tag = response_key,
      values = values
    ),
    metadata
  )
}

#' Unit suffix appended to labels.
stationary_unit_suffix <- function(response_info) {
  if (nzchar(response_info$unit)) paste0(" ", response_info$unit) else ""
}

#' Axis label with the unit in parentheses.
stationary_axis_label <- function(prefix, response_info) {
  if (nzchar(response_info$unit)) {
    sprintf("%s (%s)", prefix, response_info$unit)
  } else {
    prefix
  }
}

#' Label giving the threshold value and its unit.
stationary_threshold_label <- function(value, response_info, prefix = "Q96") {
  paste0(
    prefix, " = ",
    formatC(value, format = "f", digits = response_info$digits),
    stationary_unit_suffix(response_info)
  )
}

#' Output path of a stationary-analysis figure for the response.
stationary_figure_file <- function(response_info, filename) {
  file.path(
    figure_dir("stationary_EVA", response_info$tag),
    filename
  )
}

#' Fit a stationary GPD above the threshold, without standard errors if the information matrix is singular.
fit_stationary_gpd <- function(x, threshold) {
  tryCatch(
    evd::fpot(x, threshold = threshold, model = "gpd"),
    error = function(e) {
      if (!grepl(
        "observed information matrix is singular",
        conditionMessage(e),
        fixed = TRUE
      )) {
        stop(e)
      }

      message("GPD point estimates retained; standard errors unavailable.")
      evd::fpot(x, threshold = threshold, model = "gpd", std.err = FALSE)
    }
  )
}
