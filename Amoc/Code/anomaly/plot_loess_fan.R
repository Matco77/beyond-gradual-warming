# ============================================================
# LOESS + CONFIDENCE BAND on the trajectory, over a SPATIAL FAN
# that uses the full grid's variability WITHOUT pooling to one mean.
#
# Per variable (tas, pr), one page, one panel per model x protocol:
#   * grey fan  = AREA-WEIGHTED SPATIAL QUANTILES of the ANNUAL-mean
#     anomaly field across all Europe-window cells: 5-95% (light) and
#     25-75% (dark). This is spatial HETEROGENEITY (how differently
#     the grid responds), NOT sampling error of a mean.
#   * black     = area-weighted spatial MEDIAN trajectory.
#   * red       = loess(median ~ year, degree 1) with a pointwise
#     ±1.96*se band from predict(se = TRUE).
#
# Honesty notes (also in the caption):
#   - the loess band is POINTWISE and assumes iid residuals; annual
#     means reduce, but do not remove, interannual autocorrelation ->
#     the band is indicative, not a formal 95% envelope. The formal
#     effect test remains the matched-window ratio in the contrast CSV.
#   - fan and loess band answer different questions: spread of the
#     response across space vs uncertainty of the central trajectory.
# ============================================================

library(ncdf4)

## rename-proof: source amoc_common.R from THIS script's own folder
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: read_europe_cube, panel_grid, publish_file.

runs <- list(
  list(model = "EC-Earth3",       proto = "g01", folder = "ECHearth3_anomaly"),
  list(model = "EC-Earth3",       proto = "u03", folder = "ECHearth3_anomaly"),
  list(model = "HadGEM3-GC31-LL", proto = "g01", folder = "LL_anomaly"),
  list(model = "HadGEM3-GC31-LL", proto = "u03", folder = "LL_anomaly"),
  list(model = "HadGEM3-GC31-MM", proto = "g01", folder = "MM_anomaly"),
  list(model = "HadGEM3-GC31-MM", proto = "u03", folder = "MM_anomaly"),
  list(model = "IPSL-CM6A-LR",    proto = "u03", folder = "IPSL_anomaly")   # u03 only (no g01)
)
base   <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/anomaly_output"
outdir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"
anomaly_path <- function(r, vn) file.path(base, r$folder,
  if (r$folder %in% c("ECHearth3_anomaly", "IPSL_anomaly")) sprintf("%s_Amon_%s_hos-%s-hos_anomaly.nc", vn, r$model, r$proto)
  else sprintf("%s_anomaly_%s-hos_minus_piControl_1850-1949.nc", vn, r$proto))

annual <- function(v) {
  d <- dim(v); ny <- d[3] %/% 12
  apply(array(v[, , 1:(12 * ny), drop = FALSE], c(d[1], d[2], 12, ny)), c(1, 2, 4), mean, na.rm = TRUE)
}

## area-weighted quantiles of a field (weights = cos lat, cells NA-safe)
wquant <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w); x <- x[ok]; w <- w[ok]
  o  <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = p, rule = 2)$y
}

QP <- c(0.05, 0.25, 0.50, 0.75, 0.95)

pdf_final <- file.path(outdir, "loess_fan_grid.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 11, height = 8)

for (vn in c("tas", "pr")) {
  unit <- if (vn == "pr") "mm/day" else "K"

  FAN <- Filter(Negate(is.null), lapply(runs, function(r) {
    p <- anomaly_path(r, vn); if (!file.exists(p)) return(NULL)
    cb <- read_europe_cube(p); A <- annual(cb$v); ny <- dim(A)[3]
    W  <- matrix(cos(cb$lat * pi / 180), dim(A)[1], dim(A)[2], byrow = TRUE)
    Q  <- t(vapply(seq_len(ny), function(y) wquant(A[, , y], W, QP), numeric(length(QP))))  # ny x 5
    x  <- seq_len(ny)
    sp <- max(0.30, min(0.75, 35 / ny))
    lo <- loess(Q[, 3] ~ x, span = sp, degree = 1)
    pr <- predict(lo, data.frame(x = x), se = TRUE)
    cat(sprintf("  %-16s %s %s | median(final yr) %+.2f | loess span %.2f\n",
                r$model, r$proto, vn, Q[ny, 3], sp))
    list(r = r, x = x, Q = Q, fit = pr$fit, ci = 1.96 * pr$se.fit)
  }))
  if (!length(FAN)) next

  yl <- range(c(unlist(lapply(FAN, function(z) z$Q)), 0), na.rm = TRUE)
  par(mfrow = panel_grid(length(FAN)), oma = c(0, 0, 4, 0), mar = c(3.4, 3.6, 2, 1), mgp = c(2, 0.6, 0))
  for (z in FAN) {
    plot(NA, xlim = range(z$x), ylim = yl, xlab = "Years from start of hosing",
         ylab = sprintf("anomaly [%s]", unit), main = sprintf("%s %s", z$r$model, z$r$proto), cex.main = 0.9)
    polygon(c(z$x, rev(z$x)), c(z$Q[, 1], rev(z$Q[, 5])), col = "grey88", border = NA)  # 5-95%
    polygon(c(z$x, rev(z$x)), c(z$Q[, 2], rev(z$Q[, 4])), col = "grey72", border = NA)  # 25-75%
    abline(h = 0, lty = 3)
    lines(z$x, z$Q[, 3], col = "black", lwd = 1.4)                                      # spatial median
    polygon(c(z$x, rev(z$x)), c(z$fit - z$ci, rev(z$fit + z$ci)),
            col = adjustcolor("#B2182B", alpha.f = 0.25), border = NA)                  # loess CI
    lines(z$x, z$fit, col = "#B2182B", lwd = 2.4)                                       # loess fit
  }
  mtext(sprintf("%s: spatial fan of the full grid + LOESS of the spatial median", vn),
        outer = TRUE, line = 2.2, cex = 1.0, font = 2)
  mtext(paste("grey fan = area-weighted SPATIAL quantiles across cells (5-95% light, 25-75% dark) = heterogeneity, not sampling error;",
              "red = loess(median) ±1.96·se, pointwise, iid assumption (autocorrelation not adjusted) — indicative"),
        outer = TRUE, line = 0.9, cex = 0.72, col = "grey30")
}
dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
