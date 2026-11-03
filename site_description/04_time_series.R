library(data.table)
library(ggplot2)
library(patchwork)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
data <- read_rscd_2026()

required_variables <- c("time", "hs", "Os")
missing_variables <- setdiff(required_variables, names(data))
if (length(missing_variables)) {
  stop(
    "Missing variables for the time-series figure: ",
    paste(missing_variables, collapse = ", "),
    call. = FALSE
  )
}

setorder(data, time)

#' Time-series panel of one variable.
time_series_panel <- function(variable, colour, title, ylab, show_x = TRUE) {
  plot_data <- data.frame(
    time = data$time,
    value = data[[variable]]
  )
  plot_data <- plot_data[
    !is.na(plot_data$time) & is.finite(plot_data$value),
  ]

  p <- ggplot(plot_data, aes(x = time, y = value)) +
    geom_line(
      colour = colour,
      linewidth = 0.28,
      alpha = 0.85
    ) +
    scale_x_datetime(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.45),
      panel.grid.minor = element_blank()
    )

  if (!show_x) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )
  }

  p
}

hs_time_series <- time_series_panel(
  variable = "hs",
  colour = "#2A819F",
  title = "Significant wave height",
  ylab = expression(H[s] ~ "[m]"),
  show_x = FALSE
)

os_time_series <- time_series_panel(
  variable = "Os",
  colour = "#D06C34",
  title = "Orbital speed",
  ylab = expression(O[s] ~ "[m"~s^{-1} * "]"),
  show_x = TRUE
)

combined_time_series <- hs_time_series / os_time_series +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave(
  file.path(figure_dir("site_description"), "hs_time_series.png"),
  combined_time_series,
  width = 14,
  height = 9.5,
  dpi = 300
)

print(combined_time_series)
