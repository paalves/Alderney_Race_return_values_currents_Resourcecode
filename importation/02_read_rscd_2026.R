#' Read the hindcast CSV and keep the analysis window (1994-04-01 to 2026-03-31 UTC).
read_rscd_2026 <- function(path = "data_rscd_alderney_2026.csv") {
  data <- data.table::fread(path)
  data[, time := as.POSIXct(time, tz = "UTC")]
  data[, fje := factor(fje)]
  data[time >= as.POSIXct("1994-04-01", tz = "UTC") & time < as.POSIXct("2026-04-01", tz = "UTC")]
}

#' Length of the hourly record in years.
observation_years <- function(data) {
  (as.numeric(difftime(max(data$time), min(data$time), units = "days")) + 1 / 24) / 365.25
}
