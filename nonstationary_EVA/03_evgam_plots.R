source("paths.R")

source("importation/02_read_rscd_2026.R")
library(data.table)
library(evgam)
library(ggplot2)
library(patchwork)

dataset_tag <- "rscd_2026"
out_dir <- figure_dir("nonstationary_EVA")
files <- list(hs = c(effects = "hs_gpd.png", qq = "qq_hs.png", ald = "hs_ald.png"),
              os = c(effects = "os_gpd.png", qq = "qq_os.png", ald = "os_ald.png"))
theme_p <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))
#' Axis label of a covariate.
x_label <- function(v) if (v == "Cspd") expression(C[s]~"[m s"^-1*"]") else expression(theta[wc]~"[°]")
tags <- plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
#' Axis breaks of a covariate.
x_breaks <- function(v) if (v == "Cspd") 0:4 else seq(0, 360, 90)

#' Centred partial effect of one covariate with 95 % pointwise intervals from posterior simulation.
partial_effect <- function(fit, v, fixed, par = 1, nsim = 1000) {
  gr <- if (v == "Cspd") seq(0, 4.75, length.out = 100) else seq(0, 359, length.out = 100)
  nd <- data.table(Cspd = rep(fixed["Cspd"], 100), angleDiff = rep(fixed["angleDiff"], 100))
  set(nd, j = v, value = gr)
  sm <- simulate(fit, nsim = nsim, newdata = nd, type = "link")[[par]]
  sm <- sweep(sm, 2, colMeans(sm))
  data.table(x = gr, fit = rowMeans(sm), lo = apply(sm, 1, quantile, 0.025), hi = apply(sm, 1, quantile, 0.975))
}

#' Partial-effect panel with its interval band and a rug of the data.
effect_panel <- function(pe, rug, v, title, ylab, rug_alpha) {
  ggplot(pe, aes(x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1f4e79", alpha = 0.2) +
    geom_line(aes(y = fit), colour = "#1f4e79", linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_rug(data = data.table(x = rug), aes(x), alpha = rug_alpha, length = unit(0.03, "npc")) +
    scale_x_continuous(breaks = x_breaks(v)) +
    labs(title = title, x = x_label(v), y = ylab) +
    theme_p
}

fits <- lapply(c(hs = "hs", os = "os"), function(r)
  readRDS(file.path("outputs", sprintf("evgam_selection_%s_%s.rds", dataset_tag, r)))$best_fit)

set.seed(2026)
for (r in c("hs", "os")) {
  bf <- fits[[r]]
  pk <- bf$excess_df
  fixed <- c(Cspd = median(pk$Cspd), angleDiff = 100)
  lab <- if (r == "hs") "H[s]" else "O[s]"

  panels <- lapply(c("Cspd", "angleDiff"), function(v)
    effect_panel(partial_effect(bf$fit_gpd, v, fixed), pk[[v]], v,
                 parse(text = sprintf("log~sigma~','~%s", lab)), "Partial effect [-]", 0.4))
  ggsave(file.path(out_dir, files[[r]][["effects"]]), (panels[[1]] | panels[[2]]) + tags,
         width = 10, height = 4, dpi = 300, bg = "white")

  p <- predict(bf$fit_gpd, newdata = pk, type = "response")
  res <- sort(-log((1 + p$shape * pk$excess / p$scale)^(-1 / p$shape)))
  qq <- ggplot(data.table(th = qexp(ppoints(length(res))), obs = res), aes(th, obs)) +
    geom_abline(linetype = 2, colour = "grey50") +
    geom_point(size = 1) +
    labs(title = parse(text = lab), x = "Theoretical exponential quantile [-]", y = "Empirical residual [-]") +
    theme_p
  ggsave(file.path(out_dir, files[[r]][["qq"]]), qq, width = 4.5, height = 4.3, dpi = 300, bg = "white")
}

hourly <- read_rscd_2026()
for (r in c("hs", "os")) {
  bf <- fits[[r]]
  rug <- hourly[sample(.N, 3000)]
  fixed <- c(Cspd = median(hourly$Cspd), angleDiff = 100)
  lab <- if (r == "hs") "H[s]" else "O[s]"
  unit <- if (r == "hs") "Partial effect [m]" else expression("Partial effect [m s"^-1*"]")
  panels <- list()
  for (par in 1:2) for (v in c("Cspd", "angleDiff")) {
    title <- parse(text = if (par == 1) sprintf("u~','~%s", lab) else sprintf("log~sigma[ALD]~','~%s", lab))
    panels[[length(panels) + 1]] <- effect_panel(partial_effect(bf$fit_ald, v, fixed, par), rug[[v]], v,
                                                 title, if (par == 1) unit else "Partial effect [-]", 0.15)
  }
  ggsave(file.path(out_dir, files[[r]][["ald"]]), (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) + tags,
         width = 10, height = 7.5, dpi = 300, bg = "white")
}
