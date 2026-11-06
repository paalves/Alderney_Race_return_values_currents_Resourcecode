library(ggplot2)
library(patchwork)
library(data.table)

source("paths.R")

source("importation/02_read_rscd_2026.R")
data <- read_rscd_2026()

#' Hexagonal-bin density panel of y against x on a log colour scale.
hex_panel <- function(x, y, xlab, ylab, title) {
  plot_data <- data.frame(x = x, y = y)

  ggplot(plot_data, aes(x = x, y = y)) +
    geom_hex(bins = 80, aes(fill = after_stat(count))) +
    scale_fill_viridis_c(
      option = "viridis",
      trans = "log10",
      limits = c(1, 1000),
      oob = scales::squish,
      name = "Observations",
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1,000+"),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(4.2, "cm"),
        barheight = grid::unit(0.45, "cm")
      )
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = grid::unit(4.2, "cm"),
      legend.key.height = grid::unit(0.5, "cm")
    )
}

hs_tp <- hex_panel(
  data$hs, data$tp,
  xlab = expression(H[s] ~ "[m]"),
  ylab = expression(T[p] ~ "[s]"),
  title = "Wave peak period vs.\nsignificant wave height"
)
hs_current <- hex_panel(
  data$Cspd, data$hs,
  xlab = expression("Current speed [m"~s^{-1} * "]"),
  ylab = expression(H[s] ~ "[m]"),
  title = "Significant wave height vs.\ncurrent speed"
)
hs_angle <- hex_panel(
  data$angleDiff, data$hs,
  xlab = "Angular difference [°]",
  ylab = expression(H[s] ~ "[m]"),
  title = "Significant wave height vs.\nangular difference"
)

os_tp <- hex_panel(
  data$Os, data$tp,
  xlab = expression(O[s] ~ "[m"~s^{-1} * "]"),
  ylab = expression(T[p] ~ "[s]"),
  title = "Wave peak period vs.\norbital speed"
)
os_current <- hex_panel(
  data$Cspd, data$Os,
  xlab = expression("Current speed [m"~s^{-1} * "]"),
  ylab = expression(O[s] ~ "[m"~s^{-1} * "]"),
  title = "Orbital speed vs.\ncurrent speed"
)
os_angle <- hex_panel(
  data$angleDiff, data$Os,
  xlab = "Angular difference [°]",
  ylab = expression(O[s] ~ "[m"~s^{-1} * "]"),
  title = "Orbital speed vs.\nangular difference"
)

combined <- wrap_plots(
  hs_tp, hs_current, hs_angle,
  os_tp, os_current, os_angle,
  ncol = 3,
  guides = "collect"
) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir("site_description"), "hs_scatter_hexbin.png"),
  combined,
  width = 14,
  height = 10.5,
  dpi = 300
)

print(combined)
