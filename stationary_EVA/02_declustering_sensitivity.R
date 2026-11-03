library(data.table)
library(POT)
library(evd)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(viridis)

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

#' Runs declustering: number of clusters and their indices.
decluster_runs <- function(x, threshold, run_length) {
  exceed <- x > threshold
  if (sum(exceed) == 0) return(list(n_clusters = 0, cluster_indices = integer(0)))
  
  exceed_indices <- which(exceed)
  gaps <- diff(exceed_indices)
  new_cluster <- c(TRUE, gaps > run_length)
  cluster_id <- cumsum(new_cluster)
  
  cluster_indices <- sapply(unique(cluster_id), function(i) {
    members <- exceed_indices[cluster_id == i]
    members[which.max(x[members])]
  })
  
  list(n_clusters = length(cluster_indices), cluster_indices = cluster_indices)
}

thresholds <- as.numeric(quantile(
  response_values,
  probs = seq(0.90, 0.995, by = 0.01)
))
run_lengths <- c(12, 24, 48, 72, 96, 120)
cluster_results <- data.frame()

for (u in thresholds) {
  n_exceedances <- sum(response_values > u)
  
  for (r in run_lengths) {
    clusters <- decluster_runs(response_values, threshold = u, run_length = r)
    n_events <- clusters$n_clusters
    
    if (n_events >= 50 && n_exceedances > 0) {
      cluster_maxima <- response_values[clusters$cluster_indices]
      
      fit <- fit_stationary_gpd(cluster_maxima, threshold = u)

      has_se <- !is.null(fit$std.err) &&
        all(c("shape", "scale") %in% names(fit$std.err))

      cluster_results <- rbind(cluster_results, data.frame(
        threshold = u, run_length = r, n_exceedances = n_exceedances,
        n_clusters = n_events, cluster_ratio = n_events / n_exceedances,
        shape = fit$param["shape"], scale = fit$param["scale"],
        shape_se = if (has_se) fit$std.err["shape"] else NA_real_,
        scale_se = if (has_se) fit$std.err["scale"] else NA_real_
      ))
    }
  }
}

plot_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p1 <- ggplot(cluster_results, aes(x = threshold, y = n_clusters, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$n_clusters), 
           label = stationary_threshold_label(chosen_threshold, response_info),
           hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Number of independent storm events",
       x = stationary_axis_label("Threshold u", response_info),
       y = "Number of independent events") +
  plot_theme

p2 <- ggplot(cluster_results, aes(x = threshold, y = cluster_ratio, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, y = max(cluster_results$cluster_ratio), 
           label = stationary_threshold_label(chosen_threshold, response_info),
           hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Extremal index change",
       x = stationary_axis_label("Threshold u", response_info),
       y = "Extremal index (theta)") +
  plot_theme

p3 <- ggplot(cluster_results, aes(x = threshold, y = shape, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_ribbon(aes(ymin = shape - 1.96 * shape_se, ymax = shape + 1.96 * shape_se, 
                  fill = as.factor(run_length)), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$shape + 1.96 * cluster_results$shape_se, na.rm = TRUE), 
           label = stationary_threshold_label(chosen_threshold, response_info),
           hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Shape parameter stability",
       x = stationary_axis_label("Threshold u", response_info),
       y = expression(paste("Shape parameter ", xi))) +
  plot_theme

p4 <- ggplot(cluster_results, aes(x = threshold, y = scale, 
                                   color = as.factor(run_length), group = run_length)) +
  geom_ribbon(aes(ymin = scale - 1.96 * scale_se, ymax = scale + 1.96 * scale_se, 
                  fill = as.factor(run_length)), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = chosen_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  annotate("text", x = chosen_threshold, 
           y = max(cluster_results$scale + 1.96 * cluster_results$scale_se, na.rm = TRUE), 
           label = stationary_threshold_label(chosen_threshold, response_info),
           hjust = -0.1, size = 6, fontface="bold", color = "blue") +
  scale_color_viridis_d(name = "Run length\n(hours)") +
  scale_fill_viridis_d(name = "Run length\n(hours)") +
  labs(title = "Scale parameter stability",
       x = stationary_axis_label("Threshold u", response_info),
       y = expression(paste("Scale parameter ", sigma))) +
  plot_theme

p1 <- p1 + labs(tag = "(a)")
p2 <- p2 + labs(tag = "(b)")
p3 <- p3 + labs(tag = "(c)")
p4 <- p4 + labs(tag = "(d)")

declustering_sensitivity_plot <- arrangeGrob(
  p1, p2, p3, p4,
  nrow = 2,
  ncol = 2,
  top = textGrob(
    paste("Declustering sensitivity analysis for", response_info$label),
    gp = gpar(fontsize = 20, fontface = "bold")
  )
)

declustering_figure_file <- stationary_figure_file(
  response_info,
  "declustering_sensitivity.png"
)

ggsave(
  declustering_figure_file,
  declustering_sensitivity_plot,
  width = 13,
  height = 9,
  dpi = 300,
  bg = "white"
)

grid.newpage()
grid.draw(declustering_sensitivity_plot)
