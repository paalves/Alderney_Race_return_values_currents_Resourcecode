library(ggplot2)
library(patchwork)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "stationary_EVA", "00_stationary_helpers.R"))

response_variable <- get0("response_variable", ifnotfound = "hs", inherits = TRUE)
fit_is_missing_or_stale <-
  !exists("fit_declust", inherits = FALSE) ||
  !exists("stationary_fitted_response", inherits = FALSE) ||
  !identical(tolower(stationary_fitted_response), tolower(response_variable))

if (fit_is_missing_or_stale) {
  source(file.path(project_root(), "stationary_EVA", "03_analysis.R"))
}

#' Quantile function of the generalized Pareto distribution.
gpd_quantile <- function(probability, scale, shape) {
  if (abs(shape) < 1e-6) {
    -scale * log(1 - probability)
  } else {
    (scale / shape) * ((1 - probability)^(-shape) - 1)
  }
}

#' Density of the generalized Pareto distribution.
gpd_density <- function(x, scale, shape) {
  if (abs(shape) < 1e-6) {
    return(exp(-x / scale) / scale)
  }

  support <- 1 + shape * x / scale > 0
  density <- rep(NA_real_, length(x))
  density[support] <- (1 / scale) *
    (1 + shape * x[support] / scale)^(-1 / shape - 1)
  density
}

scale_est <- unname(fit_declust$param["scale"])
shape_est <- unname(fit_declust$param["shape"])

if (all(is.finite(fit_summary$lower[fit_summary$parameter == "shape"]))) {
  gpd_subtitle <- bquote(
    xi == .(sprintf("%.3f", shape_est)) ~ "[" ~
      .(sprintf("%.3f", fit_summary$lower[fit_summary$parameter == "shape"])) ~ "," ~
      .(sprintf("%.3f", fit_summary$upper[fit_summary$parameter == "shape"])) ~ "]" ~
      "   " ~ sigma == .(sprintf("%.3f", scale_est)) ~ "[" ~
      .(sprintf("%.3f", fit_summary$lower[fit_summary$parameter == "scale"])) ~ "," ~
      .(sprintf("%.3f", fit_summary$upper[fit_summary$parameter == "scale"])) ~ "]"
  )
} else {
  gpd_subtitle <- bquote(
    xi == .(sprintf("%.3f", shape_est)) ~ "   " ~
      sigma == .(sprintf("%.3f", scale_est)) ~
      "   (standard errors unavailable)"
  )
}

plot_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 19, hjust = 0.5),
    plot.subtitle = element_text(
      face = "bold",
      size = 16,
      colour = "grey35",
      hjust = 0.5
    ),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

qq_data <- data.frame(
  theoretical = gpd_quantile(
    (seq_along(excesses) - 0.5) / length(excesses),
    scale = scale_est,
    shape = shape_est
  ),
  empirical = sort(excesses)
)

qq_plot <- ggplot(qq_data, aes(x = theoretical, y = empirical)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey65",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    colour = "#4C8BC2",
    size = 1.5,
    alpha = 0.65
  ) +
  labs(
    title = "QQ plot",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  plot_theme

density_data <- data.frame(
  excess = seq(0, max(excesses), length.out = 400)
)
density_data$density <- gpd_density(
  density_data$excess,
  scale = scale_est,
  shape = shape_est
)

density_plot <- ggplot(data.frame(excess = excesses), aes(x = excess)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 28,
    fill = "#9FC4DF",
    colour = "white",
    linewidth = 0.2
  ) +
  geom_line(
    data = density_data,
    aes(x = excess, y = density),
    colour = "#2F75B5",
    linewidth = 1.1
  ) +
  labs(
    title = "Density plot",
    subtitle = gpd_subtitle,
    x = stationary_axis_label("Excess over threshold", response_info),
    y = "Density"
  ) +
  plot_theme

return_level_plot <- ggplot(rv_declust, aes(x = period, y = level)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#6F9FCB",
    alpha = 0.25,
    na.rm = TRUE
  ) +
  geom_line(
    colour = "#2F75B5",
    linewidth = 1
  ) +
  geom_point(
    shape = 21,
    size = 2.8,
    stroke = 1,
    colour = "#2F75B5",
    fill = "white"
  ) +
  scale_x_log10(
    breaks = return_periods,
    labels = return_periods,
    expand = expansion(mult = c(0.04, 0.05))
  ) +
  labs(
    title = "Return level plot",
    x = "Return period (years)",
    y = paste0(response_info$label, " [", response_info$unit, "]")
  ) +
  plot_theme

gpd_model_fit <- (qq_plot | density_plot | return_level_plot) +
  plot_annotation(
    title = paste("GPD model fit"),
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 24,
        hjust = 0
      )
    )
  ) +
  plot_layout(widths = c(1, 1, 1))

gpd_figure_file <- stationary_figure_file(
  response_info,
  "gpd_model_fit.png"
)

ggsave(
  gpd_figure_file,
  gpd_model_fit,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

gpd_model_fit
