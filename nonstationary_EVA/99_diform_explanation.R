library(data.table)
library(Rcpp)

.path_helper <- c("paths.R", "../paths.R", "../../paths.R")
source(.path_helper[file.exists(.path_helper)][1])
rm(.path_helper)

sourceCpp(project_file("nonstationary_EVA", "pot_decluster.cpp"))

DATA_PATH <- data_file()
N_YEARS <- 31
RETURN_PERIODS <- 100
N_DIRECTIONS <- 100
H_HOURS <- 72
ZETA <- 0.1
N_THREADS <- 0
OUT_FILE <- file.path(figure_dir("nonstationary_EVA"), "diform_Hs_Tp_geometry.png")

normalise_matrix <- function(X) {
  m <- apply(X, 2, median)
  s <- apply(X, 2, sd)
  list(Y = sweep(sweep(X, 2, m, "-"), 2, s, "/"), med = m, sd = s)
}

directions_2d <- function(n) {
  a <- seq(0, 2 * pi, length.out = n + 1)[-1]
  cbind(cos(a), sin(a))
}

fit_gp_constrained <- function(exc, max_obs = NULL) {
  exc <- exc[exc > 0]
  n <- length(exc)

  if (n < 3) return(list(sigma = mean(exc), xi = -0.1))

  m <- mean(exc)
  v <- var(exc)

  xi0 <- max(-0.99, min(-0.01, 0.5 * (1 - m^2 / v)))
  sig0 <- max(0.01, m * (1 - xi0))

  if (!is.null(max_obs))
    sig0 <- max(sig0, -xi0 * max_obs * 1.001)

  nll <- function(par) {
    sig <- par[1]
    xi <- par[2]
    z <- 1 + xi * exc / sig

    if (sig <= 0 || xi >= 0 || xi < -1 || any(z <= 0) ||
        (!is.null(max_obs) && max_obs > -sig / xi)) {
      1e12
    } else {
      n * log(sig) + (1 + 1 / xi) * sum(log(z))
    }
  }

  o <- optim(c(sig0, xi0), nll, method = "Nelder-Mead",
             control = list(maxit = 2000))

  list(sigma = o$par[1], xi = o$par[2])
}

gp_return_value <- function(sigma, xi, lambda, T) {
  if (abs(xi) < 1e-8) sigma * log(lambda * T)
  else sigma / xi * ((lambda * T)^xi - 1)
}

compute_return_values_pot_2d <- function(Y, dirs, return_periods, n_years,
                                         h_hours = 72, zeta = 0.1, n_threads = 0) {
  dec <- pot_decluster_all(cbind(Y, 0), cbind(dirs, 0), h_hours, zeta, n_threads)
  rv <- matrix(NA_real_, nrow(dirs), length(return_periods),
               dimnames = list(NULL, paste0("T", return_periods)))
  fits <- data.frame(mu = dec$mu, sigma = NA_real_, xi = NA_real_,
                     n_peaks = dec$n_peaks, n_exc = NA_integer_,
                     lambda = NA_real_, method = NA_character_)

  for (j in seq_len(nrow(dirs))) {
    peaks <- dec$peaks[[j]]
    exc <- dec$exceedances[[j]]
    fits$n_exc[j] <- length(exc)
    fits$lambda[j] <- length(exc) / n_years

    if (length(exc) < 5 || fits$lambda[j] < 0.01) {
      fits$method[j] <- "empirical"
      p <- pmin(1, pmax(0, 1 - 1 / ((length(peaks) / n_years) * return_periods)))
      rv[j, ] <- as.numeric(quantile(peaks, p, type = 7))
    } else {
      gp <- fit_gp_constrained(exc, max(peaks) - dec$mu[j])
      fits$sigma[j] <- gp$sigma
      fits$xi[j] <- gp$xi
      fits$method[j] <- "gpd"
      rv[j, ] <- dec$mu[j] + gp_return_value(gp$sigma, gp$xi, fits$lambda[j], return_periods)
    }
  }

  list(rv = rv, fits = fits, peaks = dec$peaks, exceedances = dec$exceedances)
}

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

containment_fraction <- function(Y, dirs, R) {
  mean(apply(dirs %*% t(Y) <= R + 1e-8, 2, all))
}

choose_return_period <- function(Y, dirs, rv, return_periods) {
  f <- sapply(seq_along(return_periods), \(k) containment_fraction(Y, dirs, rv[, k]))
  i <- which(f >= 1)[1]
  if (is.na(i)) i <- length(return_periods)
  list(i = i, T = return_periods[i], R = rv[, i],
       containment = f[i],
       table = data.frame(return_period = return_periods, containment_fraction = f))
}

radial_boundary_point <- function(j, dirs, R) {
  d <- drop(dirs %*% dirs[j, ])
  ok <- d > 1e-10
  a <- which(ok)[which.min(R[ok] / d[ok])]
  r <- min(R[ok] / d[ok])
  list(j = j, active = a, r = r, p = r * dirs[j, ])
}

choose_highlight_direction <- function(dirs, R, preferred = 2L) {
  a <- sapply(seq_len(nrow(dirs)), \(j) radial_boundary_point(j, dirs, R)$active)
  s <- which(a == seq_along(a))
  if (preferred %in% s) preferred else if (length(s)) s[which.min(abs(s - preferred))] else preferred
}

draw_support_line <- function(d, R, xlim, ylim, col = "gray72", lty = 2, lwd = 1) {
  x0 <- d * R / sum(d^2)
  v <- c(-d[2], d[1])
  t <- seq(-50, 50, length.out = 1500)
  x <- x0[1] + t * v[1]
  y <- x0[2] + t * v[2]
  keep <- x >= xlim[1] - 1 & x <= xlim[2] + 1 & y >= ylim[1] - 1 & y <= ylim[2] + 1
  lines(x[keep], y[keep], col = col, lty = lty, lwd = lwd)
}

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

  dev.off()
  list(polygon = P, highlighted_direction = j,
       active_constraint = rb$active, boundary_point = p)
}

dt <- fread(DATA_PATH)[is.finite(hs) & is.finite(tp) & hs > 0 & tp > 0]
setorderv(dt, intersect(names(dt), c("time", "Time", "datetime", "Datetime", "date", "Date", "timestamp"))[1])

Y <- normalise_matrix(as.matrix(dt[, .(hs, tp)]))$Y
dirs <- directions_2d(N_DIRECTIONS)

pot <- compute_return_values_pot_2d(Y, dirs, RETURN_PERIODS, N_YEARS, H_HOURS, ZETA, N_THREADS)
choice <- choose_return_period(Y, dirs, pot$rv, RETURN_PERIODS)

fig <- plot_diform_geometry_hstp(Y, dirs, choice$R, choice$T, OUT_FILE)
