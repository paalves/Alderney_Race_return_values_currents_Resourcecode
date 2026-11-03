library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(data.table)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

source(file.path(project_root(), "importation", "02_read_rscd_2026.R"))
data <- read_rscd_2026()

required_variables <- c("time", "hs", "dp", "tp", "Cspd", "Os")
missing_variables <- setdiff(required_variables, names(data))
if (length(missing_variables)) {
  stop(
    "Missing variables for the storm-evolution figure: ",
    paste(missing_variables, collapse = ", "),
    call. = FALSE
  )
}

HS_THRESHOLD <- 5.5
WINDOW       <- 50

#' Indices of storm peaks above the threshold, at least min_sep hours apart.
find_storm_peaks <- function(hs, threshold = HS_THRESHOLD, min_sep = WINDOW) {
  n      <- length(hs)
  peaks  <- c()
  in_storm <- FALSE
  best_idx <- NA
  best_val <- -Inf

  for (i in seq_len(n)) {
    if (hs[i] > threshold) {
      in_storm <- TRUE
      if (hs[i] > best_val) { best_val <- hs[i]; best_idx <- i }
    } else {
      if (in_storm) { peaks <- c(peaks, best_idx); best_val <- -Inf; best_idx <- NA }
      in_storm <- FALSE
    }
  }
  if (in_storm && !is.na(best_idx)) peaks <- c(peaks, best_idx)

  if (length(peaks) <= 1) return(peaks)
  kept <- peaks[1]
  for (p in peaks[-1]) {
    if ((p - tail(kept, 1)) >= min_sep) kept <- c(kept, p)
    else if (hs[p] > hs[tail(kept, 1)]) kept[length(kept)] <- p
  }
  kept
}

peak_idx <- find_storm_peaks(data$hs)

storm_list <- lapply(seq_along(peak_idx), function(s) {
  pi   <- peak_idx[s]
  lags <- seq(-WINDOW, WINDOW)
  rows <- pi + lags
  keep <- rows >= 1 & rows <= nrow(data)

  df <- data.frame(
    storm_id  = s,
    peak_date = format(data$time[pi], "%d/%m/%Y"),
    peak_hs   = data$hs[pi],
    lag       = lags[keep],
    hs        = data$hs[rows[keep]],
    dp        = data$dp[rows[keep]],
    tp        = data$tp[rows[keep]],
    Cspd      = data$Cspd[rows[keep]],
    Os        = data$Os[rows[keep]]
  )
  df
})

storm_df <- bind_rows(storm_list)
storm_summary <- storm_df %>%
  distinct(storm_id, peak_date, peak_hs) %>%
  filter(peak_hs > HS_THRESHOLD) %>%
  arrange(storm_id) %>%
  mutate(label = sprintf("%s (Hs: %.1f m)", peak_date, peak_hs))

storm_df <- storm_df %>%
  select(-any_of("label")) %>%
  left_join(storm_summary %>% select(storm_id, label), by = "storm_id")

ordered_labels <- storm_summary$label

storm_summary$label <- factor(storm_summary$label, levels = ordered_labels)
storm_df$label <- factor(storm_df$label, levels = ordered_labels)

n_storms <- nrow(storm_summary)
storm_colors <- setNames(
  hue_pal()(n_storms),
  ordered_labels
)

#' Mean and mean plus or minus one standard deviation across storms at each lag.
make_ribbon <- function(df, var) {
  df %>%
    group_by(lag) %>%
    summarise(
      mn  = mean(.data[[var]], na.rm = TRUE),
      sdd = sd(.data[[var]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(lo = mn - sdd, hi = mn + sdd)
}

rib_hs   <- make_ribbon(storm_df, "hs")
rib_dp   <- make_ribbon(storm_df, "dp")
rib_tp   <- make_ribbon(storm_df, "tp")
rib_cspd <- make_ribbon(storm_df, "Cspd")
rib_os   <- make_ribbon(storm_df, "Os")

BASE <- 16

#' Panel of the storm trajectories with their mean and spread.
storm_panel <- function(var, ribbon_df, ylab, title,
                        extra_layers = NULL, xlab = NULL) {
  p <- ggplot() +
    geom_ribbon(data = ribbon_df,
                aes(x = lag, ymin = lo, ymax = hi),
                fill = "grey60", alpha = 0.9) +
    geom_line(data = storm_df,
              aes(x = lag, y = .data[[var]], colour = label),
              linewidth = 0.3, alpha = 0.85) +
    geom_line(data = ribbon_df,
              aes(x = lag, y = mn),
              colour = "black", linewidth = 1.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30",
               linewidth = 0.7) +
    scale_colour_manual(values = storm_colors, name = "Storm peak (date & max Hs)") +
    labs(title = title,
         x = xlab, y = ylab) +
    theme_minimal(base_size = BASE) +
    theme(
      plot.title       = element_text(face = "bold", size = BASE + 2),
      axis.title       = element_text(face = "bold", size = BASE + 1),
      axis.text        = element_text(size = BASE),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    )

  if (!is.null(extra_layers)) p <- p + extra_layers
  p
}

pa <- storm_panel(
  "hs", rib_hs,
  ylab     = expression(H[s]~"[m]"),
  title    = "Significant wave height evolution",
  extra_layers = geom_hline(yintercept = HS_THRESHOLD, linetype = "dotted",
                            colour = "grey50", linewidth = 0.6)
)

pb <- storm_panel(
  "dp", rib_dp,
  ylab     =  expression(D[p]~"[°]"),
  title    = "Peak wave direction evolution"
)

pc <- storm_panel(
  "tp", rib_tp,
  ylab     = expression(T[p]~"[s]"),
  title    = "Peak wave period evolution",
  xlab     = "Observations relative to Hs peak"
)

pd <- storm_panel(
  "Cspd", rib_cspd,
  ylab     = expression("Current speed [m"~s^{-1}*"]"),
  title    = "Current speed evolution"
)

pe <- storm_panel(
  "Os", rib_os,
  ylab     = expression(O[s]~"[m"~s^{-1}*"]"),
  title    = "Orbital speed evolution"
)

pa <- pa + labs(title = paste("(a)", pa$labels$title))
pe <- pe + labs(title = paste("(b)", pe$labels$title))
pb <- pb + labs(title = "(c) Wave direction evolution")
pc <- pc + labs(title = paste("(d)", pc$labels$title))
pd <- pd + labs(title = paste("(e)", pd$labels$title))

library(cowplot)

legend_panel <- ggplot(storm_summary,
                       aes(x = storm_id, y = peak_hs, colour = label)) +
  geom_point(size = 2) +
  scale_colour_manual(values = storm_colors,
                      name   = "Storm peak (date & max Hs)") +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = BASE - 1),
    legend.title     = element_text(face = "bold", size = BASE),
    legend.key.width = unit(1.1, "lines"),
    legend.key.height = unit(1.1, "lines"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(t = 8)
  )
shared_legend <- cowplot::get_plot_component(legend_panel, "guide-box-bottom")
if (is.null(shared_legend)) {
  shared_legend <- cowplot::get_plot_component(legend_panel, "^guide-box$")
}

grid_top <- plot_grid(pa, pe, ncol = 2, align = "hv")
grid_bottom <- plot_grid(pb, pc, pd, ncol = 3, align = "hv")
grid <- plot_grid(
  grid_top,
  grid_bottom,
  ncol = 1,
  rel_heights = c(1, 1)
)

title_gg <- ggdraw() +
  draw_label("Severe storms: evolution around Hs peak",
             fontface = "bold", size = BASE + 6, x = 0.5, hjust = 0.5)

final <- plot_grid(
  title_gg,
  grid,
  shared_legend,
  ncol    = 1,
  rel_heights = c(0.06, 1, 0.28)
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(figure_dir("site_description"), "storms.png"),
       final, width = 14, height = 13, dpi = 300, bg = "white")
