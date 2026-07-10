library(ggplot2)
library(patchwork)
library(data.table)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

data <- fread(data_file())

hex_panel <- function(x, y, xlab, title,
                      ylab = expression(H[s]~"(m)"),
                      breaks = c(1, 10, 100, 1000)) {
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_hex(bins = 80, aes(fill = after_stat(count))) +
    scale_fill_viridis_c(
      option = "turbo",
      trans  = "log10",
      name   = "Observations",
      breaks = breaks,
      labels = scales::comma(breaks),
      guide  = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(3.4, "cm"),
        barheight = grid::unit(0.45, "cm")
      )
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 18, face = "bold"),
      axis.title       = element_text(size = 21),
      axis.text        = element_text(size = 18),
      legend.position  = "bottom",
      legend.title     = element_text(size = 19, face = "bold"),
      legend.text      = element_text(size = 17),
      legend.key.width = grid::unit(3.4, "cm"),
      legend.key.height = grid::unit(0.55, "cm")
    )
}

pa <- hex_panel(data$hs,        data$tp, expression(H[s]~"(m)"),
                "Wave peak period\nvs. significant wave height",
                ylab = expression(T[p]~"(s)"))
pb <- hex_panel(data$Cspd,      data$hs, "Current speed (m/s)",
                "Significant wave height\nvs. current speed",
                breaks = c(1, 10, 100))
pc <- hex_panel(data$angleDiff, data$hs, "Angular difference (deg)",
                "Significant wave height\nvs. angular difference",
                breaks = c(1, 10, 100, 1000))

combined <- pa + pb + pc + plot_layout(ncol = 3)

ggsave(file.path(figure_dir("site_description"), "hs_scatter_hexbin.png"),
       combined, width = 14, height = 5.5, dpi = 300)
