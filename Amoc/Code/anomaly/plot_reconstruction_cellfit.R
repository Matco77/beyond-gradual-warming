# ============================================================
# PER-CELL FIT QUALITY of the tasmin/tasmax reconstruction
# (EC-Earth3, Europe window; mirrors the per-cell per-calendar-month
# OLS of 1.reconstruct_tasmin_tasmax_ECEarth3.py, incl. its variance
# guard var <= 1e-6 K^2 -> slope fallback 1).
#
# Five map pages answer "is the fit good in THIS cell?":
#  (1) R^2 per cell: MEDIAN over the 12 calendar months, and the
#      WORST month. R^2 = share of interannual variance explained.
#  (2) In-sample residual RMSE [K]: median and worst month. This is
#      the error scale that actually enters the reconstruction --
#      low R^2 with tiny RMSE (quiet maritime cells) is harmless.
#  (3) OUT-OF-SAMPLE check: fit on piControl years 1..250, score on
#      the remaining years (temporal split -> also tests the
#      stationarity assumption). Maps: out-of-sample RMSE and its
#      RATIO to in-sample (ratio ~ 1 = no overfitting/drift).
#  (4) Fitted slopes b (median over months) for tasmin and tasmax --
#      diverging scale centred on 1.
#  (5) EXTRAPOLATION EXPOSURE: % of hosing months whose tas falls
#      OUTSIDE the cell's piControl min-max for that calendar month
#      (per experiment) -- where the line is evaluated on trust.
# Console: area-weighted summary quantiles of every metric.
# ============================================================

library(ncdf4)
library(fields)
have_maps <- requireNamespace("maps", quietly = TRUE)

.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: pic_root, read_europe_cube, publish_file.

pic_dir <- file.path(pic_root, "EC-Earth3")
hos_dir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/NAHosMIP/EC-Earth3"
outdir  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"

files_by_year <- function(v) {
  f <- list.files(pic_dir, pattern = sprintf("^%s_Amon_EC-Earth3_piControl_r1i1p1f1_gr_", v), full.names = TRUE)
  names(f) <- sub(".*_gr_(\\d{4})01-.*", "\\1", f); f
}
Ft <- files_by_year("tas"); Fx <- files_by_year("tasmax"); Fn <- files_by_year("tasmin")
yy <- sort(intersect(intersect(names(Ft), names(Fx)), names(Fn)))
cat("piControl years:", length(yy), sprintf("(%s-%s)\n", head(yy, 1), tail(yy, 1)))

## build cell x time matrices (Europe crop; float64)
read_year_mat <- function(f) { cb <- read_europe_cube(f); list(cb = cb, M = matrix(as.numeric(cb$v), ncol = 12)) }
grab <- function(F) {
  cols <- vector("list", length(yy))
  for (i in seq_along(yy)) { cols[[i]] <- read_year_mat(F[[yy[i]]])$M
    if (i %% 100 == 0) cat("  ...", i, "years\n") }
  do.call(cbind, cols)
}
cb0 <- read_europe_cube(Ft[[yy[1]]]); lon <- cb0$lon; lat <- cb0$lat
nl <- length(lon); na <- length(lat)
cat("Europe window:", nl, "x", na, "cells\n")
cat("reading tas ...\n");    Xall <- grab(Ft)
cat("reading tasmax ...\n"); Yx   <- grab(Fx)
cat("reading tasmin ...\n"); Yn   <- grab(Fn)
nyr <- length(yy); mo <- rep(1:12, nyr)

GUARD <- 1e-6   # same threshold as the Python script

## per-month OLS moments on a chosen set of years; returns list of 12 x cells stats
fit_stats <- function(Y, yrs_fit, yrs_val) {
  r2  <- rmse_in <- rmse_out <- b12 <- matrix(NA_real_, nrow(Y), 12)
  for (m in 1:12) {
    cf <- (yrs_fit - 1) * 12 + m; cv <- (yrs_val - 1) * 12 + m
    X <- Xall[, cf, drop = FALSE]; Yy <- Y[, cf, drop = FALSE]
    mx <- rowMeans(X); my <- rowMeans(Yy)
    xa <- X - mx; ya <- Yy - my
    vx <- rowMeans(xa^2); vy <- rowMeans(ya^2); cv2 <- rowMeans(xa * ya)
    b  <- ifelse(vx <= GUARD, 1.0, cv2 / vx)           # guard mirrors the .py
    a  <- my - b * mx
    r2[, m]      <- ifelse(vy <= GUARD, NA, (cv2^2) / (vx * vy))
    rmse_in[, m] <- sqrt(pmax(vy - 2 * b * cv2 + b^2 * vx, 0))
    b12[, m]     <- b
    if (length(cv)) {
      Xv <- Xall[, cv, drop = FALSE]; Yv <- Y[, cv, drop = FALSE]
      rmse_out[, m] <- sqrt(rowMeans((Yv - a - b * Xv)^2))
    }
  }
  list(r2 = r2, rmse_in = rmse_in, rmse_out = rmse_out, b = b12)
}

all_yrs <- seq_len(nyr); half <- seq_len(nyr %/% 2)
Sx  <- fit_stats(Yx, all_yrs, integer(0))              # full fit (as production)
Sn  <- fit_stats(Yn, all_yrs, integer(0))
Vx  <- fit_stats(Yx, half, setdiff(all_yrs, half))     # split-sample validation
Vn  <- fit_stats(Yn, half, setdiff(all_yrs, half))

## extrapolation exposure of the hosing tas, per experiment
hos_files <- c(g01 = file.path(hos_dir, "tas_Amon_EC-Earth3_hos-g01-hos_185001-189912.nc"),
               u03 = file.path(hos_dir, "tas_Amon_EC-Earth3_hos-u03-hos_185001-194912.nc"))
extrap <- lapply(hos_files, function(f) {
  cb <- read_europe_cube(f); H <- matrix(as.numeric(cb$v), nrow(Xall))
  hmo <- rep(1:12, ncol(H) / 12); out <- matrix(NA_real_, nrow(H), 12)
  for (m in 1:12) {
    cp <- (all_yrs - 1) * 12 + m
    lo <- apply(Xall[, cp], 1, min); hi <- apply(Xall[, cp], 1, max)
    Hm <- H[, hmo == m, drop = FALSE]
    out[, m] <- rowMeans(Hm < lo | Hm > hi)
  }
  100 * rowMeans(out)                                   # % of all hosing months outside support
})

## ---- summaries (area-weighted quantiles) ----
W <- as.vector(matrix(cos(lat * pi / 180), nl, na, byrow = TRUE))
wq <- function(x, p = c(.05, .5, .95)) { ok <- is.finite(x)
  o <- order(x[ok]); xx <- x[ok][o]; ww <- W[ok][o]
  approx((cumsum(ww) - 0.5 * ww) / sum(ww), xx, xout = p, rule = 2)$y }
summ <- function(lab, x) cat(sprintf("  %-34s q05 %6.3f   median %6.3f   q95 %6.3f\n", lab, wq(x)[1], wq(x)[2], wq(x)[3]))
cat("\n=== per-cell fit quality (area-weighted quantiles over the window) ===\n")
summ("tasmin R2 (median month)",   apply(Sn$r2, 1, median, na.rm = TRUE))
summ("tasmin R2 (worst month)",    apply(Sn$r2, 1, min,    na.rm = TRUE))
summ("tasmax R2 (median month)",   apply(Sx$r2, 1, median, na.rm = TRUE))
summ("tasmax R2 (worst month)",    apply(Sx$r2, 1, min,    na.rm = TRUE))
summ("tasmin RMSE in-sample [K]",  apply(Sn$rmse_in, 1, median))
summ("tasmax RMSE in-sample [K]",  apply(Sx$rmse_in, 1, median))
summ("tasmin RMSE out-of-sample",  apply(Vn$rmse_out, 1, median))
summ("tasmax RMSE out-of-sample",  apply(Vx$rmse_out, 1, median))
summ("tasmin slope b (median)",    apply(Sn$b, 1, median))
summ("tasmax slope b (median)",    apply(Sx$b, 1, median))
summ("extrapolation % g01",        extrap$g01)
summ("extrapolation % u03",        extrap$u03)
cat(sprintf("  guarded cells (any month, b==1 fallback): tasmin %d, tasmax %d\n",
            sum(apply(Sn$b == 1, 1, any)), sum(apply(Sx$b == 1, 1, any))))

## ---- map pages ----
fld <- function(x) matrix(x, nl, na)
draw <- function(z, ttl, pal, zlim, unit) {
  image.plot(lon, lat, fld(z), col = pal, zlim = zlim, xlab = "lon", ylab = "lat",
             main = ttl, cex.main = 0.9, legend.lab = unit, legend.line = 2.3)
  if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
}
pal_r2   <- hcl.colors(64, "viridis")
pal_rmse <- hcl.colors(64, "Reds 3", rev = TRUE)
pal_b    <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(64)

pdf_final <- file.path(outdir, "reconstruction_cellfit.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 11, height = 8.5)
par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(3.2, 3.2, 2.2, 4.4), mgp = c(2, 0.6, 0))

draw(apply(Sn$r2, 1, median, na.rm = TRUE), "tasmin R2, median month", pal_r2, c(0, 1), "R2")
draw(apply(Sn$r2, 1, min,    na.rm = TRUE), "tasmin R2, worst month",  pal_r2, c(0, 1), "R2")
draw(apply(Sx$r2, 1, median, na.rm = TRUE), "tasmax R2, median month", pal_r2, c(0, 1), "R2")
draw(apply(Sx$r2, 1, min,    na.rm = TRUE), "tasmax R2, worst month",  pal_r2, c(0, 1), "R2")
mtext("Per-cell per-month fit: share of interannual variance explained (piControl, full record)", outer = TRUE, cex = 0.95, font = 2)

zr <- c(0, as.numeric(quantile(c(apply(Sn$rmse_in, 1, max), apply(Sx$rmse_in, 1, max)), 0.98, na.rm = TRUE)))
draw(apply(Sn$rmse_in, 1, median), "tasmin RMSE, median month [K]", pal_rmse, zr, "K")
draw(apply(Sn$rmse_in, 1, max),    "tasmin RMSE, worst month [K]",  pal_rmse, zr, "K")
draw(apply(Sx$rmse_in, 1, median), "tasmax RMSE, median month [K]", pal_rmse, zr, "K")
draw(apply(Sx$rmse_in, 1, max),    "tasmax RMSE, worst month [K]",  pal_rmse, zr, "K")
mtext("In-sample residual RMSE: the reconstruction's error scale in kelvin (low R2 + low RMSE = harmless quiet cell)", outer = TRUE, cex = 0.95, font = 2)

ratio_n <- apply(Vn$rmse_out, 1, median) / pmax(apply(Vn$rmse_in, 1, median), 1e-9)
ratio_x <- apply(Vx$rmse_out, 1, median) / pmax(apply(Vx$rmse_in, 1, median), 1e-9)
draw(apply(Vn$rmse_out, 1, median), "tasmin OUT-of-sample RMSE [K]", pal_rmse, zr, "K")
draw(pmin(ratio_n, 2), "tasmin out/in RMSE ratio (1 = no drift)", pal_b, c(0.5, 1.5), "ratio")
draw(apply(Vx$rmse_out, 1, median), "tasmax OUT-of-sample RMSE [K]", pal_rmse, zr, "K")
draw(pmin(ratio_x, 2), "tasmax out/in RMSE ratio (1 = no drift)", pal_b, c(0.5, 1.5), "ratio")
mtext("Temporal split-sample validation: fit on piControl years 1-250, scored on the rest (tests fit AND stationarity)", outer = TRUE, cex = 0.95, font = 2)

draw(apply(Sn$b, 1, median), "tasmin slope b, median month", pal_b, c(0, 2), "b")
draw(apply(Sn$b, 1, mad),    "tasmin slope spread across months (MAD)", pal_rmse, c(0, 0.5), "K/K")
draw(apply(Sx$b, 1, median), "tasmax slope b, median month", pal_b, c(0, 2), "b")
draw(apply(Sx$b, 1, mad),    "tasmax slope spread across months (MAD)", pal_rmse, c(0, 0.5), "K/K")
mtext("Fitted slopes (diverging scale centred on b = 1) and their seasonal spread", outer = TRUE, cex = 0.95, font = 2)

par(mfrow = c(1, 2))
draw(extrap$g01, "g01: % hosing months outside piControl range", pal_rmse, c(0, 100), "%")
draw(extrap$u03, "u03: % hosing months outside piControl range", pal_rmse, c(0, 100), "%")
mtext("Extrapolation exposure: where the fitted line is evaluated beyond its training support", outer = TRUE, cex = 0.95, font = 2)

dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
