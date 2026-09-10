# ============================================================
# PDFs of the HOSING runs' ABSOLUTE tas / tasmin / tasmax [K].
# No piControl anywhere. READ-ONLY.
#
# One page, one panel per model x protocol; per panel three kernel
# densities over ALL Europe-window cell-months of the FULL run,
# area-weighted: tas (dark grey), tasmin (blue), tasmax (red).
# Subtitle: weighted mean ± sd of each + mean diurnal-range proxy
# DTR = mean(tasmax) - mean(tasmin).
#
# EC-Earth3 tasmin/tasmax are the RECONSTRUCTED files: identical
# shape family to tas with a clean ±3 K offset = the QC pass;
# spikes/clipping/order violations would be the failure signs.
# ============================================================

library(ncdf4)

.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: read_europe_cube, panel_grid, publish_file.

DS     <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
outdir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"

runs <- list(
  list(model = "EC-Earth3",       proto = "g01"),
  list(model = "EC-Earth3",       proto = "u03"),
  list(model = "HadGEM3-GC31-LL", proto = "g01"),
  list(model = "HadGEM3-GC31-LL", proto = "u03"),
  list(model = "HadGEM3-GC31-MM", proto = "g01"),
  list(model = "HadGEM3-GC31-MM", proto = "u03"),
  list(model = "IPSL-CM6A-LR",    proto = "u03")   # u03 only (no g01)
)
hos_files <- function(model, vn, proto) {
  dir <- file.path(DS, "NAHosMIP", model)
  if (model == "EC-Earth3") {
    if (vn == "tas") list.files(dir, pattern = sprintf("^tas_Amon_EC-Earth3_hos-%s-hos_\\d", proto), full.names = TRUE)
    else            list.files(dir, pattern = sprintf("^%s_Amon_EC-Earth3_hos-%s-hos_reconstructed", vn, proto), full.names = TRUE)
  } else if (model == "IPSL-CM6A-LR") {
    # IPSL: tas/pr native (DRS name _u03-hos_), tasmin/tasmax reconstructed (like EC)
    if (vn %in% c("tas", "pr")) list.files(dir, pattern = sprintf("^%s_Amon_IPSL-CM6A-LR_%s-hos_", vn, proto), full.names = TRUE)
    else                        list.files(dir, pattern = sprintf("^%s_Amon_IPSL-CM6A-LR_hos-%s-hos_reconstructed", vn, proto), full.names = TRUE)
  } else list.files(dir, pattern = sprintf("^%s_Amon_%s_%s-hos_", vn, model, proto), full.names = TRUE)
}

pool <- function(files) {
  xs <- list(); ws <- list()
  for (f in files) {
    cb <- read_europe_cube(f)
    xs[[f]] <- as.vector(cb$v)
    ws[[f]] <- rep(as.vector(matrix(cos(cb$lat * pi / 180), dim(cb$v)[1], dim(cb$v)[2], byrow = TRUE)), dim(cb$v)[3])
  }
  x <- unlist(xs, use.names = FALSE); w <- unlist(ws, use.names = FALSE)
  ok <- is.finite(x); x <- x[ok]; w <- w[ok] / sum(w[ok])
  list(x = x, w = w, mean = sum(x * w), sd = sqrt(sum(w * (x - sum(x * w))^2)))
}

VCOL <- c(tas = "grey20", tasmin = "#2166AC", tasmax = "#B2182B")

pdf_final <- file.path(outdir, "pdf_hosing_triple.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 11, height = 8)
par(mfrow = panel_grid(length(runs)), oma = c(0, 0, 4, 0), mar = c(3.2, 3.4, 2.6, 0.8), mgp = c(2, 0.6, 0))

for (r in runs) {
  sets <- lapply(c(tas = "tas", tasmin = "tasmin", tasmax = "tasmax"),
                 function(vn) pool(hos_files(r$model, vn, r$proto)))
  dens <- lapply(sets, function(s) density(s$x, weights = s$w, n = 1024))
  xr <- range(unlist(lapply(dens, `[[`, "x"))); ymax <- max(unlist(lapply(dens, function(d) max(d$y))))
  plot(NA, xlim = xr, ylim = c(0, ymax * 1.05), xlab = "K", ylab = "density",
       main = sprintf("%s %s%s", r$model, r$proto,
                      if (r$model == "EC-Earth3") "  (min/max reconstructed)" else ""), cex.main = 0.9)
  for (k in names(VCOL)) lines(dens[[k]], col = VCOL[k], lwd = 2)
  mtext(sprintf("tas %.1f±%.2f | min %.1f±%.2f | max %.1f±%.2f | DTR %.2f K",
                sets$tas$mean, sets$tas$sd, sets$tasmin$mean, sets$tasmin$sd,
                sets$tasmax$mean, sets$tasmax$sd, sets$tasmax$mean - sets$tasmin$mean),
        side = 3, line = 0.1, cex = 0.55, col = "grey30")
  cat(sprintf("  %-16s %s | tas %7.2f ± %5.2f | min %7.2f ± %5.2f | max %7.2f ± %5.2f | DTR %5.2f\n",
              r$model, r$proto, sets$tas$mean, sets$tas$sd, sets$tasmin$mean, sets$tasmin$sd,
              sets$tasmax$mean, sets$tasmax$sd, sets$tasmax$mean - sets$tasmin$mean))
}
for (k in seq_along(VCOL)) mtext(names(VCOL)[k], col = VCOL[k], side = 3, line = 0.6,
                                 at = seq(0.35, 0.65, length.out = 3)[k], outer = TRUE, cex = 0.9, font = 2)
mtext("HOSING runs only: absolute tas / tasmin / tasmax [K], full run, all Europe cell-months (area-weighted)",
      outer = TRUE, line = 2.2, cex = 1.0, font = 2)
mtext("QC: three curves must be same shape family, ordered min < tas < max (~3 K offsets); spikes/clipping = weirdness",
      outer = TRUE, line = 1.0, cex = 0.72, col = "grey30")
dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
