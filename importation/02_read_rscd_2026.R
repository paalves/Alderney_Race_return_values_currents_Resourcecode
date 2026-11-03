#' Start of the common analysis window.
analysis_start_time <- function() {
  as.POSIXct("1994-04-01 00:00:00", tz = "UTC")
}

#' Last day of the common analysis window.
analysis_end_date <- function() {
  as.Date("2026-03-31")
}

#' Upper time bound of the analysis window.
analysis_end_time <- function() {
  as.POSIXct(analysis_end_date() + 1, tz = "UTC")
}

#' Keep the rows inside the analysis window.
filter_rscd_analysis_period <- function(data, time_col = "time") {
  if (!time_col %in% names(data)) return(data)

  time <- as.POSIXct(data[[time_col]], tz = "UTC")
  keep <- !is.na(time) &
    time >= analysis_start_time() &
    time < analysis_end_time()
  data[which(keep), ]
}

#' Format times, derive the wave-current angle and tidal phase when missing and restrict to the analysis window.
prepare_rscd_2026 <- function(data, apply_analysis_period = TRUE) {
  data <- data.table::as.data.table(data)

  if ("time" %in% names(data)) {
    data[, time := as.POSIXct(time, tz = "UTC")]
  }

  if (!"angleDiff" %in% names(data) &&
      all(c("dp", "Cdirdeg") %in% names(data))) {
    wave_to_deg <- (data$dp + 180) %% 360
    data[, angleDiff := (wave_to_deg - Cdirdeg) %% 360]
  }

  if (!"fje" %in% names(data) && all(c("Cspd", "Cdirdeg") %in% names(data))) {
    data[, fje := data.table::fifelse(Cspd < 1, "étale", data.table::fifelse(Cdirdeg %/% 180 == 1, "jusant", "flot"))]
  }

  if ("fje" %in% names(data)) {
    data[, fje := factor(fje)]
  }

  if (isTRUE(apply_analysis_period)) {
    data <- filter_rscd_analysis_period(data)
  }

  data
}

#' Read and prepare the ResourceCode hindcast CSV file.
read_rscd_2026 <- function(path = NULL) {
  if (is.null(path)) path <- data_file("data_rscd_alderney_2026.csv")
  prepare_rscd_2026(data.table::fread(path))
}

#' Length of the record in years, from the time span and sampling step.
observation_years <- function(data, time_col = "time") {
  if (!time_col %in% names(data)) return(NA_real_)
  tt <- sort(unique(as.POSIXct(data[[time_col]], tz = "UTC")))
  tt <- tt[!is.na(tt)]
  if (length(tt) < 2) return(NA_real_)
  step_days <- stats::median(as.numeric(diff(tt), units = "days"), na.rm = TRUE)
  span_days <- as.numeric(max(tt) - min(tt), units = "days")
  if (is.finite(step_days) && step_days > 0) span_days <- span_days + step_days
  span_days / 365.25
}
