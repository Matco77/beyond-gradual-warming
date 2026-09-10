# ============================================================
# PDFs (kernel densities) of the tas / tasmin / tasmax ANOMALIES.
# READ-ONLY on the anomaly outputs -- nothing is modified.
#
# One page, one panel per model x protocol. Per panel, three curves:
#   tas (black), tasmin (blue), tasmax (red) --
# the distribution of LATE-RUN (final-third) monthly anomalies over
# ALL Europe-window cells, area-weighted (cos lat). Dotted line at 0.
#
# What to look for ("come cambiano"):
#   * the whole distribution shifting left = the cooling;
#   * tasmin shifting MORE than tasmax (fitted slopes 1.045 vs 0.968)
#     -> the diurnal temperature range widens slightly under cooling;
#     the panel subtitle prints the three means and d(DTR) =
#     mean(tasmax) - mean(tasmin) anomaly difference.
# Densities are computed on the full data; only the DISPLAY x-range is
# trimmed to the pooled 0.2-99.8% quantiles so sea-ice tails do not
# flatten every curve.
# ============================================================

library(ncdf4)

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

VCOL <- c(tas = "grey20", tasmin = "#2166AC", tasmax = "#B2182B")

## late-run (final-third) pooled values + matching cos-lat weights
late_vals <- function(r, vn) {
  p <- anomaly_path(r, vn); if (!file.exists(p)) return(NULL)
  cb <- read_europe_cube(p); nt <- dim(cb$v)[3]; pl <- max(12, round(nt / 3))
  v  <- cb$v[, , (nt - pl + 1):nt, drop = FALSE]
  w  <- rep(as.vector(matrix(cos(cb$lat * pi / 180), dim(v)[1], dim(v)[2], byrow = TRUE)), dim(v)[3])
  list(x = as.vector(v), w = w / sum(w))
}

D <- lapply(runs, function(r) {
  out <- lapply(c("tas", "tasmin", "tasmax"), function(vn) late_vals(r, vn))
  names(out) <- c("tas", "tasmin", "tasmax")
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  mns <- vapply(out, function(z) sum(z$x * z$w), numeric(1))
  cat(sprintf("  %-16s %s | mean tas %+.2f  tasmin %+.2f  tasmax %+.2f  -> d(DTR) %+.3f K\n",
              r$model, r$proto, mns["tas"], mns["tasmin"], mns["tasmax"],
              mns["tasmax"] - mns["tasmin"]))
  list(r = r, d = lapply(out, function(z) density(z$x, weights = z$w, n = 1024)),
       mns = mns)
})
D <- Filter(Negate(is.null), D)

## shared display range: pooled 0.2-99.8% quantiles (display only)
allx <- unlist(lapply(D, function(z) lapply(z$d, function(dd) dd$x)))
xr   <- as.numeric(quantile(allx, c(0.002, 0.998)))

pdf_final <- file.path(outdir, "pdf_tas_triple.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 11, height = 8)
par(mfrow = panel_grid(length(D)), oma = c(0, 0, 4, 0), mar = c(3.4, 3.6, 2.6, 1), mgp = c(2, 0.6, 0))
for (z in D) {
  ymax <- max(vapply(z$d, function(dd) max(dd$y), numeric(1)))
  plot(NA, xlim = xr, ylim = c(0, ymax * 1.04), xlab = "anomaly [K]", ylab = "density",
       main = sprintf("%s %s", z$r$model, z$r$proto), cex.main = 0.95)
  abline(v = 0, lty = 3, col = "grey55")
  for (vn in names(VCOL)) lines(z$d[[vn]], col = VCOL[vn], lwd = 2)
  mtext(sprintf("means: tas %+.2f  min %+.2f  max %+.2f   |   d(DTR) %+.2f K",
                z$mns["tas"], z$mns["tasmin"], z$mns["tasmax"],
                z$mns["tasmax"] - z$mns["tasmin"]),
        side = 3, line = 0.1, cex = 0.62, col = "grey30")
}
for (k in seq_along(VCOL)) mtext(names(VCOL)[k], col = VCOL[k], side = 3, line = 0.6,
                                 at = seq(0.3, 0.7, length.out = 3)[k], outer = TRUE, cex = 0.9, font = 2)
mtext("PDFs of late-run (final-third) monthly anomalies -- all Europe-window cell-months, area-weighted",
      outer = TRUE, line = 2.2, cex = 1.0, font = 2)
mtext("dotted = 0 (no change); display x-range trimmed to pooled 0.2-99.8% quantiles (densities computed on all data)",
      outer = TRUE, line = 1.0, cex = 0.72, col = "grey30")
dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
