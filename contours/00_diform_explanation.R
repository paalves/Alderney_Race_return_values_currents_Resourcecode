source("paths.R")
source("importation/02_read_rscd_2026.R")
source("contours/01_diform_contours_3d.R")

RETURN_PERIODS <- 100
N_DIRECTIONS <- 100
H_HOURS <- 72
ZETA <- 0.1
N_THREADS <- 0
OUT_FILE <- file.path(figure_dir("contours"), "diform_Hs_Tp_geometry.png")

#' Return n unit vectors evenly spaced on the circle.
directions_2d <- function(n) {
  a <- seq(0, 2 * pi, length.out = n + 1)[-1]
  cbind(cos(a), sin(a))
}

#' Vertices of the polygon formed by the intersection of the directional half-planes.
halfspace_polygon_2d <- function(dirs, R) {
  pts <- do.call(rbind, lapply(seq_len(nrow(dirs) - 1), function(i) {
    do.call(rbind, lapply((i + 1):nrow(dirs), function(j) {
      A <- rbind(dirs[i, ], dirs[j, ])
      if (abs(det(A)) < 1e-10) return(NULL)
      x <- solve(A, c(R[i], R[j]))
      if (all(dirs %*% x <= R + 1e-8)) x
    }))
  }))
  pts <- unique(round(pts, 10))
  ctr <- colMeans(pts)
  pts[order(atan2(pts[, 2] - ctr[2], pts[, 1] - ctr[1])), ]
}

#' Contour point along direction j and index of its active constraint.
radial_boundary_point <- function(j, dirs, R) {
  d <- drop(dirs %*% dirs[j, ])
  ok <- d > 1e-10
  a <- which(ok)[which.min(R[ok] / d[ok])]
  r <- min(R[ok] / d[ok])
  list(j = j, active = a, r = r, p = r * dirs[j, ])
}

#' Pick a direction that is its own active constraint, close to the preferred one.
choose_highlight_direction <- function(dirs, R, preferred = 2L) {
  a <- sapply(seq_len(nrow(dirs)), \(j) radial_boundary_point(j, dirs, R)$active)
  s <- which(a == seq_along(a))
  if (preferred %in% s) preferred else if (length(s)) s[which.min(abs(s - preferred))] else preferred
}

#' Draw the boundary line of a directional half-plane within the plot limits.
draw_support_line <- function(d, R, xlim, ylim, col = "gray72", lty = 2, lwd = 1) {
  x0 <- d * R / sum(d^2)
  v <- c(-d[2], d[1])
  t <- seq(-50, 50, length.out = 1500)
  x <- x0[1] + t * v[1]
  y <- x0[2] + t * v[2]
  keep <- x >= xlim[1] - 1 & x <= xlim[2] + 1 & y >= ylim[1] - 1 & y <= ylim[2] + 1
  lines(x[keep], y[keep], col = col, lty = lty, lwd = lwd)
}

#' Draw the two-panel illustration of the D-IFORM construction in normalized (Hs, Tp) space.
plot_diform_geometry_hstp <- function(Y, dirs, R, T, out_file) {
  P <- halfspace_polygon_2d(dirs, R)
  j <- choose_highlight_direction(dirs, R)
  rb <- radial_boundary_point(j, dirs, R)
  p <- rb$p

  xlim <- range(c(Y[, 1], P[, 1])) + c(-0.6, 0.6)
  ylim <- range(c(Y[, 2], P[, 2])) + c(-0.6, 0.6)
  dx <- diff(xlim)
  dy <- diff(ylim)
  al <- 0.075 * min(dx, dy)

  png(out_file, width = 1800, height = 900, res = 180)
  par(mfrow = c(1, 2), mar = c(3.8, 3.4, 3.2, 1.2),
      xaxs = "i", yaxs = "i", cex = 0.82, cex.main = 0.98)

  plot(NA, NA, xlim = xlim, ylim = ylim, asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = "D-IFORM contour construction")
  mtext("(a)", side = 3, adj = 0, line = 1.2, font = 2)
  points(Y, pch = 16, cex = 0.5, col = adjustcolor("black", 0.75))
  for (i in seq_len(nrow(dirs))) draw_support_line(dirs[i, ], R[i], xlim, ylim)
  polygon(P, col = adjustcolor("#9ecae1", 0.42), border = "#08519c", lwd = 2)
  points(0, 0, pch = 19, cex = 0.75)
  text(0.015 * dx, -0.035 * dy, "0", cex = 0.76)
  for (i in seq_len(nrow(dirs)))
    arrows(0, 0, al * dirs[i, 1], al * dirs[i, 2],
           length = 0.055, lwd = 0.9, col = "#cb181d")
  for (i in unique(round(seq(1, nrow(dirs), length.out = 6))))
    text(1.25 * al * dirs[i, 1], 1.25 * al * dirs[i, 2],
         bquote(d[.(i)]), col = "#cb181d", cex = 0.76)
  text(4, 3.25, bquote(C(T)), cex = 1, font = 2, col = "#08306b")

  plot(NA, NA, xlim = xlim, ylim = ylim, asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = "D-IFORM interpretation for one direction")
  mtext("(b)", side = 3, adj = 0, line = 1.2, font = 2)
  points(Y, pch = 16, cex = 0.5, col = adjustcolor("black", 0.75))
  polygon(P, col = adjustcolor("#9ecae1", 0.42), border = "#08519c", lwd = 2)
  points(0, 0, pch = 19, cex = 0.75)
  text(0.015 * dx, -0.05 * dy, "0", cex = 0.76)

  for (i in seq_len(nrow(dirs)))
    arrows(0, 0, 0.8 * al * dirs[i, 1], 0.8 * al * dirs[i, 2],
           length = 0.045, lwd = 0.75, col = adjustcolor("#cb181d", 0.28))

  draw_support_line(dirs[rb$active, ], R[rb$active], xlim, ylim,
                    col = "#238b45", lty = 1, lwd = 1.55)
  segments(0, 0, p[1], p[2], lwd = 1.8, col = "#e31a1c")
  arrows(0, 0, 1.08 * p[1], 1.08 * p[2], length = 0.065, lwd = 1.8, col = "#e31a1c")
  points(p[1], p[2], pch = 19, cex = 1, col = "#e31a1c")
  text(1.13 * p[1], 1.13 * p[2] + 0.025 * dy, bquote(d[.(j)]),
       col = "#d95f0e", cex = 1)
  text(p[1] + 0.025 * dx, p[2] + 0.005 * dy,
       bquote(p[.(j)] == r[.(j)]^"*" * d[.(j)]),
       col = "#e31a1c", cex = 1, pos = 4)
  text(0.5 * p[1], 0.5 * p[2] - 0.055 * dy,
       bquote(r[.(j)]^"*"), col = "#e31a1c", cex = 1)
  text(xlim[1] + 0.5 * dx, ylim[2] - 0.08 * dy,
       bquote("active constraint: " ~ d[.(rb$active)]^T %.% x == R[.(rb$active)](T)),
       col = "#238b45", cex = 1, adj = 0)

  invisible(dev.off())
}

dt <- read_rscd_2026()
Y <- normalise_diform(as.matrix(dt[, .(hs, tp)]))$Y
dirs <- directions_2d(N_DIRECTIONS)
pot <- compute_return_values_pot(cbind(Y, 0), cbind(dirs, 0), RETURN_PERIODS, observation_years(dt), H_HOURS, ZETA, N_THREADS)
plot_diform_geometry_hstp(Y, dirs, pot$rv[, 1], RETURN_PERIODS, OUT_FILE)
