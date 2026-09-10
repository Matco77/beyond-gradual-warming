# ============================================================
# RECONSTRUCTION DIAGNOSTIC: piControl tas vs tasmin/tasmax scatter
# with fitted lines. PURE piControl -- no hosing data on any axis.
#
# 2 x 2 page (EC-Earth3 only -- HadGEM has native extremes):
#   rows    = tasmin, tasmax
#   columns = (1) RAW monthly pairs: the full seasonal range, pooled
#                 OLS line, dotted identity;
#             (2) MONTH-CENTRED pairs: (tas - its month mean) vs
#                 (extreme - its month mean) -- the INTERANNUAL
#                 anomalies the per-calendar-month regression of
#                 1.reconstruct_tasmin_tasmax_ECEarth3.py actually
#                 fits. The raw column's tightness is mostly the
#                 seasonal cycle; this column shows the signal that
#                 determines the fitted slopes.
# Dots coloured by calendar month (Jan -> Dec). Data = Europe-mean
# months of piControl years 2259-2757 (inner join of the three
# variables, mirroring xr.align "inner" in the Python script).
#
# HONESTY: the real reconstruction is fitted PER CELL and PER MONTH;
# this figure is the Europe-mean ILLUSTRATION of that relationship.
# It is a diagnostic, not the fit itself.
# ============================================================

library(ncdf4)

## rename-proof: source amoc_common.R from THIS script's own folder
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: pic_root, europe_monthly_means, publish_file.

pic_dir <- file.path(pic_root, "EC-Earth3")
hos_dir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/NAHosMIP/EC-Earth3"
outdir  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"

files_by_year <- function(v) {
  f <- list.files(pic_dir, pattern = sprintf("^%s_Amon_EC-Earth3_piControl_r1i1p1f1_gr_", v), full.names = TRUE)
  names(f) <- sub(".*_gr_(\\d{4})01-.*", "\\1", f); f
}
Ft <- files_by_year("tas"); Fx <- files_by_year("tasmax"); Fn <- files_by_year("tasmin")
yy <- sort(intersect(intersect(names(Ft), names(Fx)), names(Fn)))
cat("piControl years common to tas/tasmax/tasmin:", length(yy),
    sprintf("(%s-%s)\n", head(yy, 1), tail(yy, 1)))

## Europe-mean monthly series, year by year (each file = 12 months)
Tm <- Xm <- Nm <- numeric(0)
for (i in seq_along(yy)) {
  y  <- yy[i]
  Tm <- c(Tm, europe_monthly_means(Ft[[y]]))
  Xm <- c(Xm, europe_monthly_means(Fx[[y]]))
  Nm <- c(Nm, europe_monthly_means(Fn[[y]]))
  if (i %% 100 == 0) cat("  ...", i, "years\n")
}
mo <- rep(1:12, length(yy))

mcol <- hcl.colors(12, "viridis")

panel_raw <- function(ext, extlab) {
  fit <- lm(ext ~ Tm); r2 <- summary(fit)$r.squared
  plot(Tm, ext, col = mcol[mo], pch = 16, cex = 0.3,
       xlab = "piControl tas, Europe mean [K]", ylab = sprintf("piControl %s, Europe mean [K]", extlab),
       main = sprintf("%s ~ tas  (raw months)", extlab), cex.main = 0.95)
  abline(0, 1, lty = 3, col = "grey55")
  abline(fit, lwd = 2)
  mtext(sprintf("pooled OLS: b=%.3f, R2=%.3f, n=%d  (range mostly = seasonal cycle)",
                coef(fit)[2], r2, length(ext)), side = 3, line = 0.1, cex = 0.62, col = "grey30")
  cat(sprintf("  %s raw     | pooled b=%.3f R2=%.3f\n", extlab, coef(fit)[2], r2))
}

panel_centred <- function(ext, extlab) {
  xa <- Tm - ave(Tm, mo); ya <- ext - ave(ext, mo)      # per-calendar-month anomalies
  fit <- lm(ya ~ xa); r2 <- summary(fit)$r.squared
  bm  <- vapply(1:12, function(m) unname(coef(lm(ya[mo == m] ~ xa[mo == m]))[2]), numeric(1))
  plot(xa, ya, col = mcol[mo], pch = 16, cex = 0.3,
       xlab = "tas anomaly vs its month mean [K]", ylab = sprintf("%s anomaly vs its month mean [K]", extlab),
       main = sprintf("%s ~ tas  (month-centred = what the fit learns)", extlab), cex.main = 0.95)
  abline(h = 0, v = 0, lty = 3, col = "grey55")
  abline(fit, lwd = 2)
  mtext(sprintf("pooled b=%.3f, R2=%.3f   |   per-month slopes %.2f..%.2f",
                coef(fit)[2], r2, min(bm), max(bm)), side = 3, line = 0.1, cex = 0.62, col = "grey30")
  cat(sprintf("  %s centred | pooled b=%.3f R2=%.3f | monthly b in [%.2f, %.2f]\n",
              extlab, coef(fit)[2], r2, min(bm), max(bm)))
}

pdf_final <- file.path(outdir, "reconstruction_scatter.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 10.5, height = 9)
par(mfrow = c(2, 2), oma = c(0, 0, 4, 0), mar = c(3.6, 3.6, 2.6, 1), mgp = c(2.1, 0.6, 0))
panel_raw(Nm, "tasmin"); panel_centred(Nm, "tasmin")
panel_raw(Xm, "tasmax"); panel_centred(Xm, "tasmax")
mtext("EC-Earth3 piControl ONLY: tas vs tasmin/tasmax, Europe-mean months (coloured Jan→Dec)",
      outer = TRUE, line = 2.4, cex = 1.0, font = 2)
mtext(paste("left = raw monthly pairs (seasonal range), pooled OLS + dotted identity;",
            "right = month-centred interannual anomalies: the relationship the per-cell per-month regression actually fits"),
      outer = TRUE, line = 1.1, cex = 0.72, col = "grey30")
dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
