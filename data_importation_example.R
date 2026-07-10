#example of data importation from ResourceCode using the data point for the paper.
#exact localization of the data is available at https://resourcecode.ifremer.fr/explore?pointId=158438

#the data I use in the paper includes the 2021-2024 time period, which isn't available through the 'resourcecode' and 'resourcecodedata' R packages yet.

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

required_columns <- c("time", parameters)

data$time <- as.POSIXct(data$time, tz = "UTC")
data <- data[order(data$time), ]

data$Cspd <- sqrt(data$ucur^2 + data$vcur^2)
data$Cdir <- atan2(data$vcur, data$ucur)

data$Cdirdeg <- (90 - data$Cdir * 180 / pi) %% 360
data$Cdirdeg[!is.finite(data$Cspd) | data$Cspd == 0] <- NA_real_

data$wave_to_deg <- (data$dp + 180) %% 360
data$angularDifference <- (data$wave_to_deg - data$Cdirdeg) %% 360

data$Wspd <- sqrt(data$uwnd^2 + data$vwnd^2)
data$Wdir <- atan2(data$vwnd, data$uwnd)
data$Wdirdeg <- (90 - data$Wdir * 180 / pi) %% 360


head(data)

data2024 <- read.csv("data_alderney_rscd_2024.csv")
data2024 <- data2024[as.POSIXct(data2024$time, tz = "UTC") < as.POSIXct("2021-01-01", tz = "UTC"), , drop = FALSE]
data <- data[-1,] #1994-01-01 00:00:00 is missing in dataset data_alderney_rscd_2024.csv

plot(data$hs, data2024$hs)
plot(data$tp, data2024$tp)
plot(data$angularDifference, data2024$angleDiff)