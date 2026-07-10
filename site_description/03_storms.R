library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(data.table)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

data <- fread(data_file())

HS_THRESHOLD <- 5.5
WINDOW       <- 50

find_storm_peaks <- function(hs, threshold = HS_THRESHOLD, min_sep = WINDOW) {
  n      <- length(hs)
  peaks  <- c()
  in_storm <- FALSE
  best_idx <- NA
  best_val <- -Inf

  for (i in seq_len(n)) {
    if (hs[i] >= threshold) {
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
    peak_hs   = round(data$hs[pi], 1),
    lag       = lags[keep],
    hs        = data$hs[rows[keep]],
    dp        = data$dp[rows[keep]],
    tp        = data$tp[rows[keep]],
    Cspd      = data$Cspd[rows[keep]]
  )
  df
})

storm_df <- bind_rows(storm_list)
storm_df$label <- paste0(storm_df$peak_date, " (Hs: ", storm_df$peak_hs, "m)")

ordered_labels <- storm_df %>%
  distinct(storm_id, label) %>%
  arrange(storm_id) %>%
  pull(label)

storm_df$label <- factor(storm_df$label, levels = ordered_labels)

n_storms <- length(peak_idx)
storm_colors <- setNames(
  hue_pal()(n_storms),
  ordered_labels
)

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

BASE <- 12

storm_panel <- function(var, ribbon_df, ylab, title, subtitle,
                        extra_layers = NULL) {
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
    labs(title = title, subtitle = subtitle,
         x = "Observations relative to Hs peak", y = ylab) +
    theme_minimal(base_size = BASE) +
    theme(
      plot.title       = element_text(face = "bold", size = BASE + 1),
      plot.subtitle    = element_text(colour = "grey45", size = BASE - 2),
      axis.title       = element_text(face = "bold", size = BASE),
      axis.text        = element_text(size = BASE - 2),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    )

  if (!is.null(extra_layers)) p <- p + extra_layers
  p
}

pa <- storm_panel(
  "hs", rib_hs,
  ylab     = expression(H[s]~"[m]"),
  title    = "Significant wave height (Hs) evolution\nand mean curve + sd",
  subtitle = "Time relative to peak Hs (t=0 at maximum wave height)",
  extra_layers = geom_hline(yintercept = HS_THRESHOLD, linetype = "dotted",
                            colour = "grey50", linewidth = 0.6)
)

pb <- storm_panel(
  "dp", rib_dp,
  ylab     =  expression(D[p]~"[deg]"),
  title    = "Peak wave direction (Dp) evolution\nand mean curve + sd",
  subtitle = "0 deg = North, 90 deg = East, 180 deg = South, 270 deg = West"
)

pc <- storm_panel(
  "tp", rib_tp,
  ylab     = expression(T[p]~"[s]"),
  title    = "Peak wave period (Tp) evolution\nand mean curve + sd",
  subtitle = "Wave period relative to storm peak"
)

pd <- storm_panel(
  "Cspd", rib_cspd,
  ylab     = expression("Current speed [m.s"^{-1}*"]"),
  title    = "Current speed (Cspd) evolution\nand mean curve + sd",
  subtitle = "Current velocity relative to storm peak"
)

legend_panel <- ggplot(storm_df,
                       aes(x = lag, y = hs, colour = label)) +
  geom_line() +
  scale_colour_manual(values = storm_colors,
                      name   = "Storm peak (date & max Hs)") +
  theme_minimal(base_size = BASE - 1) +
theme(legend.position   = "bottom",
      legend.text       = element_text(size = BASE - 4),  
      legend.title      = element_text(face = "bold", size = BASE - 3),
      legend.key.width  = unit(0.5, "lines"),   
      legend.key.height = unit(0.5, "lines"),   
      legend.spacing.x  = unit(0.1, "cm"),      
      legend.spacing.y  = unit(0.05, "cm"),      
      legend.margin     = margin(t = 8),
      legend.byrow      = TRUE)
shared_legend <- cowplot::get_legend(legend_panel)

library(cowplot)

grid <- plot_grid(pa, pb, pc, pd, ncol = 2, align = "hv")

title_gg <- ggdraw() +
  draw_label("Severe storms: evolution around Hs peak",
             fontface = "bold", size = BASE + 4, x = 0.5, hjust = 0.5)

final <- plot_grid(
  title_gg,
  grid,
  shared_legend,
  ncol    = 1,
  rel_heights = c(0.04, 1, 0.12)
)

ggsave(file.path(figure_dir("site_description"), "severe_storms_evolution.png"),
       final, width = 12, height = 9, dpi = 300)
