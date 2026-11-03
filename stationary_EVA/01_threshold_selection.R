library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
source(file.path(project_root(), "stationary_EVA", "00_stationary_helpers.R"))
data <- read_rscd_2026()

response_variable <- get0("response_variable", ifnotfound = "hs", inherits = TRUE)
response_info <- stationary_response_info(data, response_variable)
response_values <- response_info$values

chosen_threshold <- as.numeric(quantile(response_values, 0.96))
thresholds <- as.numeric(quantile(
  response_values,
  probs = seq(0.90, 0.995, by = 0.005)
))

theme_diag <- theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

#' Vertical line and label marking the chosen threshold.
vline_threshold <- function(y_pos, label_text) {
  list(
    geom_vline(xintercept = chosen_threshold, linetype = "dashed",
               color = "blue", linewidth = 1),
    annotate("text", x = chosen_threshold, y = y_pos,
             label = label_text,
             hjust = -0.08, vjust = 1, size = 4.5, color = "blue",
             fontface = "bold")
  )
}

mrl_data <- do.call(rbind, lapply(thresholds, function(u) {
  excesses <- response_values[response_values > u] - u
  n  <- length(excesses)
  me <- mean(excesses)
  se <- sd(excesses) / sqrt(n)
  data.frame(threshold = u, mean_excess = me,
             lower = me - 1.96 * se, upper = me + 1.96 * se)
}))

p1 <- ggplot(mrl_data, aes(x = threshold, y = mean_excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  vline_threshold(max(mrl_data$upper),
                  paste(
                    "Chosen threshold",
                    stationary_threshold_label(chosen_threshold, response_info),
                    sep = "\n"
                  )) +
  labs(
    title    = "Mean residual life plot",
    x        = stationary_axis_label("Threshold u", response_info),
    y        = stationary_axis_label("Mean excess", response_info)
  ) +
  theme_diag

param_data <- do.call(rbind, lapply(thresholds, function(u) {
  excesses <- response_values[response_values > u] - u
  if (length(excesses) <= 50) return(NULL)

  fit <- tryCatch(
    fit_stationary_gpd(response_values, threshold = u),
    error = function(e) {
      if (!grepl("observed information matrix is singular",
                 conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }

      message(sprintf(
        "Threshold %.3f: standard errors unavailable; keeping point estimates.",
        u
      ))
      evd::fpot(response_values, threshold = u, model = "gpd", std.err = FALSE)
    }
  )

  has_se <- !is.null(fit$std.err) &&
    all(c("scale", "shape") %in% names(fit$std.err))
  scale_se <- if (has_se) unname(fit$std.err["scale"]) else NA_real_
  shape_se <- if (has_se) unname(fit$std.err["shape"]) else NA_real_
  scale_est <- unname(fit$param["scale"])
  shape_est <- unname(fit$param["shape"])

  data.frame(
    threshold     = u,
    scale         = scale_est,
    shape         = shape_est,
    scale_lower   = scale_est - 1.96 * scale_se,
    scale_upper   = scale_est + 1.96 * scale_se,
    shape_lower   = shape_est - 1.96 * shape_se,
    shape_upper   = shape_est + 1.96 * shape_se,
    n_exceedances = length(excesses)
  )
}))

p2 <- ggplot(param_data, aes(x = threshold, y = shape)) +
  geom_ribbon(aes(ymin = shape_lower, ymax = shape_upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  vline_threshold(max(param_data$shape_upper, na.rm = TRUE),
                  stationary_threshold_label(chosen_threshold, response_info)) +
  labs(
    title    = "Shape parameter stability",
    x        = stationary_axis_label("Threshold u", response_info),
    y        = expression(paste("Shape parameter ", xi))
  ) +
  theme_diag

p3 <- ggplot(param_data, aes(x = threshold, y = scale)) +
  geom_ribbon(aes(ymin = scale_lower, ymax = scale_upper), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  vline_threshold(max(param_data$scale_upper, na.rm = TRUE),
                  stationary_threshold_label(chosen_threshold, response_info)) +
  labs(
    title    = "Modified scale parameter stability",
    x        = stationary_axis_label("Threshold u", response_info),
    y        = expression(paste("Modified scale ", sigma^"*"))
  ) +
  theme_diag

p4 <- ggplot(param_data, aes(x = threshold, y = n_exceedances)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "darkred", alpha = 0.7) +
  annotate("text", x = min(param_data$threshold), y = 115,
           label = "Minimum recommended", hjust = 0, size = 4, color = "darkred") +
  vline_threshold(max(param_data$n_exceedances),
                  stationary_threshold_label(chosen_threshold, response_info)) +
  scale_y_log10(labels = comma) +
  labs(
    title    = "Number of threshold exceedances",
    x        = stationary_axis_label("Threshold u", response_info),
    y        = "Number of exceedances"
  ) +
  theme_diag

p1 <- p1 + labs(tag = "(a)")
p2 <- p2 + labs(tag = "(b)")
p3 <- p3 + labs(tag = "(c)")

combined_plot <- grid.arrange(
  p1, p2, p3,
  nrow = 1,
  top  = textGrob(
    paste("GPD threshold selection diagnostics for", response_info$label),
    gp = gpar(fontsize = 19, fontface = "bold")
  )
)

threshold_figure_file <- stationary_figure_file(
  response_info,
  "threshold_selection_diagnostics.png"
)

ggsave(
  threshold_figure_file,
  combined_plot,
  width = 16,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

combined_plot

fit_chosen <- fit_stationary_gpd(response_values, threshold = chosen_threshold)
fit_chosen_has_se <- !is.null(fit_chosen$std.err) &&
  all(c("scale", "shape") %in% names(fit_chosen$std.err))

threshold_summary <- data.frame(
  response = response_info$name,
  threshold = chosen_threshold,
  n_exceedances = sum(response_values > chosen_threshold),
  percent_exceedances = 100 * mean(response_values > chosen_threshold),
  shape = fit_chosen$param["shape"],
  shape_ci_half_width = if (fit_chosen_has_se) {
    1.96 * fit_chosen$std.err["shape"]
  } else {
    NA_real_
  },
  scale = fit_chosen$param["scale"],
  scale_ci_half_width = if (fit_chosen_has_se) {
    1.96 * fit_chosen$std.err["scale"]
  } else {
    NA_real_
  }
)

param_data$shape_ci_width <- param_data$shape_upper - param_data$shape_lower
stable <- subset(param_data, n_exceedances >= 100)

if (nrow(stable) > 0) {
  best <- stable[which.min(stable$shape_ci_width), ]
  alternative_threshold <- data.frame(
    response = response_info$name,
    threshold = best$threshold,
    n_exceedances = best$n_exceedances,
    percent_exceedances = 100 * mean(response_values > best$threshold),
    shape = best$shape,
    shape_ci_half_width = best$shape_upper - best$shape,
    scale = best$scale,
    scale_ci_half_width = best$scale_upper - best$scale
  )
} else {
  alternative_threshold <- NULL
}

data_summary <- data.frame(
  response = response_info$name,
  n_observations = length(response_values),
  mean = mean(response_values),
  maximum = max(response_values),
  q95 = as.numeric(quantile(response_values, 0.95)),
  q96 = as.numeric(quantile(response_values, 0.96)),
  q99 = as.numeric(quantile(response_values, 0.99))
)
