library(data.table)
library(ggplot2)
library(gridExtra)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

data <- fread(data_file())

data[, TideCondition := fifelse(
  fje == "flot", "Flood",
  fifelse(fje == "jusant", "Ebb", "Slack")
)]
data[, TideCondition := factor(TideCondition, levels = c("Slack", "Flood", "Ebb"))]

data[, dir_bin := as.numeric(as.character(cut(
  Cdirdeg %% 360,
  breaks = seq(0, 360, by = 10),
  include.lowest = TRUE,
  labels = seq(5, 355, by = 10)
)))]

theme_dist <- theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

theme_rose <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "gray20", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "gray50", size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

tide_colors <- c("Slack" = "#66BB6A", "Flood" = "#5C6BC0", "Ebb" = "#EF5350")

make_hist <- function(var, fill_col, title, xlab, adjust = 1.5) {
  brk <- seq(0, max(data[[var]], na.rm = TRUE), length.out = 50)
  ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      breaks = brk,
      fill = fill_col,
      color = "black",
      linewidth = 0.2
    ) +
    geom_density(linewidth = 0.75, color = "red", adjust = adjust) +
    labs(title = title, x = xlab, y = "Density") +
    theme_dist
}

p1 <- make_hist("hs", "lightblue", "Significant wave height distribution", expression(H[s]~"[m]"))
p2 <- make_hist("tp", "lightgreen", "Wave peak period distribution", expression(T[p]~"[s]"))
p3 <- make_hist("Cspd", "lavender", "Current speed distribution", "Current speed [m/s]", adjust = 2)
p4 <- make_hist("OrbSpeed", "orange", "Orbital speed distribution", "Orbital speed [m/s]")

rose_scale_x <- scale_x_continuous(
  limits = c(0, 360),
  breaks = seq(0, 315, 45),
  labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
)

p_rose_cur <- ggplot(data, aes(x = dir_bin)) +
  geom_bar(aes(fill = after_stat(count)), width = 10) +
  coord_polar(start = 0) +
  rose_scale_x +
  scale_fill_gradient(low = "#E0F7FA", high = "purple") +
  labs(title = "Current direction rose") +
  theme_rose

p_rose_wave <- ggplot(data, aes(x = dp)) +
  geom_histogram(bins = 16, center = 11.25, fill = "darkblue", color = "white", linewidth = 0.1) +
  coord_polar(start = 0) +
  rose_scale_x +
  labs(title = "Wave direction rose") +
  theme_rose

p_hist_diff <- ggplot(data, aes(x = angleDiff, fill = TideCondition, color = TideCondition)) +
  geom_histogram(
    breaks = seq(0, 360, by = 5),
    linewidth = 0.1,
    alpha = 0.45,
    aes(y = after_stat(density))
  ) +
  geom_density(linewidth = 0.8, fill = NA, adjust = 1.5, key_glyph = "path") +
  scale_fill_manual(values = tide_colors) +
  scale_color_manual(values = tide_colors) +
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45)) +
  labs(
    title = "Angular difference between wave and current directions",
    subtitle = "180 deg means waves oppose currents. 0/360 deg means waves follow currents. 90/270 deg means waves are perpendicular to currents.",
    x = "Angular difference [deg]",
    y = "Density",
    fill = "Tide condition",
    color = "Tide condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "gray20"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80"),
    legend.position = "top",
    legend.box.background = element_rect(color = "gray90", fill = "white")
  )

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

grid.arrange(
  arrangeGrob(p_rose_cur, p_rose_wave, ncol = 2),
  p_hist_diff,
  nrow = 2,
  heights = c(1, 0.8)
)

data[, time := as.POSIXct(time)]
data_sub <- data[261450:261550]

tide_rle <- rle(as.character(data_sub$TideCondition))
end_idx <- cumsum(tide_rle$lengths)
start_idx <- end_idx - tide_rle$lengths + 1

tide_bg <- data.table(
  condition = tide_rle$values,
  xmin = data_sub$time[start_idx],
  xmax = data_sub$time[end_idx]
)

p_timeseries <- ggplot() +
  geom_rect(
    data = tide_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
    alpha = 0.45,
    inherit.aes = FALSE
  ) +
  geom_line(data = data_sub, aes(x = time, y = wlv, color = "Water level [m]"), linewidth = 0.7) +
  geom_line(data = data_sub, aes(x = time, y = Cspd, color = "Current speed [m.s-1]"), linewidth = 0.7) +
  scale_fill_manual(
    values = c("Slack" = "white", "Flood" = "#C8E6C9", "Ebb" = "#FFF9C4"),
    name = "Tidal condition"
  ) +
  scale_color_manual(
    values = c("Water level [m]" = "#2196F3", "Current speed [m.s-1]" = "#F44336"),
    name = NULL
  ) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%d %b.\n%Hh") +
  labs(
    title = "Water levels and current speed over four days in the Alderney Race",
    subtitle = "Oct 29 - Nov 3, 2023",
    x = NULL,
    y = NULL,
    caption = "Flood/ebb distinction based on La maree oceanique cotiere, Bernard Simon"
  ) +
  guides(
    fill = guide_legend(order = 2, override.aes = list(alpha = 0.45)),
    color = guide_legend(order = 1)
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, color = "gray40", size = 10)
  )

p_timeseries
