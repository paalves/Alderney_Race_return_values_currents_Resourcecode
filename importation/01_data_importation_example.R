node_id <- 158348L
start_time <- as.POSIXct("1994-01-01 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2020-12-31 23:00:00", tz = "UTC")

parameters <- c(
  "hs",
  "dp",
  "tp",
  "wlv",
  "ucur", "vcur",
  "uwnd", "vwnd",
  "dpt"
)

data <- resourcecode::get_parameters(
  node = node_id,
  start = start_time,
  end = end_time,
  parameters = parameters
)

data$time <- as.POSIXct(data$time, tz = "UTC")
data <- data[order(data$time), ]

data$Cspd <- sqrt(data$ucur^2 + data$vcur^2)
data$Cdir <- atan2(data$vcur, data$ucur)

data$Cdirdeg <- (90 - data$Cdir * 180 / pi) %% 360
data$Cdirdeg[!is.finite(data$Cspd) | data$Cspd == 0] <- NA_real_
data$fje <- factor(ifelse(data$Cspd < 1, "étale", ifelse(data$Cdirdeg %/% 180 == 1, "jusant", "flot")))

data$wave_to_deg <- (data$dp + 180) %% 360
data$angleDiff <- (data$wave_to_deg - data$Cdirdeg) %% 360

data$Wspd <- sqrt(data$uwnd^2 + data$vwnd^2)
data$Wdir <- atan2(data$vwnd, data$uwnd)
data$Wdirdeg <- (90 - data$Wdir * 180 / pi) %% 360

data <- data[-1, ]

#' Download 1D spectra for a ResourceCode point and compute near-bed orbital speed at height z.
compute_orbital_speed <- function(point = "W001995N49729", start = "1994-01-01", end = "1994-02-01", z = 10) {
  spectra_1d <- resourcecode::get_1d_spectrum(point = point, start = format(start), end = format(end))
  orbital_speeds <- resourcecode::compute_orbital_speeds(
    spec = t(spectra_1d$ef),
    freq = spectra_1d$freq,
    depth = spectra_1d$forcings$dpt,
    z = z
  )
  orbital_data <- data.frame(
    time = as.POSIXct(spectra_1d$forcings$time, tz = "UTC"),
    Ou = orbital_speeds[, 1],
    Ov = orbital_speeds[, 2]
  )
  orbital_data$Os <- sqrt(orbital_data$Ou^2 + orbital_data$Ov^2)
  orbital_data
}
