# ============================================================
# QC PDFs of ABSOLUTE tas / tasmin / tasmax [K] -- piControl vs the
# FULL hosing runs (no windowing). READ-ONLY on all datasets.
#
# 3 x 3 page: rows = models, columns = tas / tasmin / tasmax.
# Per panel three kernel densities over all Europe-window cell-months
# (area-weighted): piControl reference (grey), hos g01 (orange),
# hos u03 (red). Subtitle prints weighted mean ± sd of each curve.
#
# Purpose: sanity check. Weirdness would show as spikes, clipping,
# unphysical bumps, or a variance collapse -- especially in the
# EC-Earth3 tasmin/tasmax, which are RECONSTRUCTED files.
# The pooled distribution mixes seasons and locations, so it is wide
# and multi-modal by nature; compare curves, not absolute shape.
#
# piControl reference = the anomaly baseline of each model:
#   EC-Earth3 full control (yearly files); HadGEM 1850-1949 chunks.
# ============================================================

library(ncdf4)

.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: pic_root, read_europe_cube, publish_file.

DS     <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
outdir <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"

MODELS <- c("EC-Earth3", "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "IPSL-CM6A-LR")
VARS   <- c("tas", "tasmin", "tasmax")

## which files feed each curve
pic_files <- function(model, vn) {
  f <- sort(list.files(file.path(pic_root, model), pattern = paste0("^", vn, "_Amon_"), full.names = TRUE))
  if (model == "EC-Earth3") f else {
    tok <- sub(".*_(\\d{6})-(\\d{6})\\.nc$", "\\1", f)
    f[tok >= "185001" & tok <= "194901"]           # chunks starting inside 1850-1949
  }
}
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

## pool values + cos-lat weights from a set of files
pool <- function(files) {
  xs <- list(); ws <- list()
  for (f in files) {
    cb <- read_europe_cube(f)
    xs[[f]] <- as.vector(cb$v)
    ws[[f]] <- rep(as.vector(matrix(cos(cb$lat * pi / 180), dim(cb$v)[1], dim(cb$v)[2], byrow = TRUE)), dim(cb$v)[3])
  }
  x <- unlist(xs, use.names = FALSE); w <- unlist(ws, use.names = FALSE)
  ok <- is.finite(x); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(list(x = numeric(0), w = numeric(0), mean = NA_real_, sd = NA_real_))  # e.g. IPSL has no g01
  w <- w / sum(w); mu <- sum(x * w)
  list(x = x, w = w, mean = mu, sd = sqrt(sum(w * (x - mu)^2)))
}

CCOL <- c(pic = "grey35", g01 = "#E69F00", u03 = "#D55E00")

pdf_final <- file.path(outdir, "pdf_tas_absolute.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 12, height = 9)
par(mfrow = c(4, 3), oma = c(0, 0, 4, 0), mar = c(3.2, 3.4, 2.6, 0.8), mgp = c(2, 0.6, 0))

for (model in MODELS) for (vn in VARS) {
  sets <- list(pic = pool(pic_files(model, vn)),
               g01 = pool(hos_files(model, vn, "g01")),
               u03 = pool(hos_files(model, vn, "u03")))
  dens <- lapply(sets, function(s) if (length(s$x) > 1) density(s$x, weights = s$w, n = 1024) else NULL)
  have <- !vapply(dens, is.null, logical(1))            # IPSL has no g01 -> skip that curve
  xr <- range(unlist(lapply(dens[have], `[[`, "x"))); ymax <- max(unlist(lapply(dens[have], function(d) max(d$y))))
  plot(NA, xlim = xr, ylim = c(0, ymax * 1.05), xlab = "K", ylab = "density",
       main = sprintf("%s  %s%s", model, vn,
                      if (model %in% c("EC-Earth3", "IPSL-CM6A-LR") && vn != "tas") "  (reconstructed)" else ""), cex.main = 0.9)
  for (k in names(sets)) if (!is.null(dens[[k]])) lines(dens[[k]], col = CCOL[k], lwd = 2)
  mtext(sprintf("pic %.1f±%.2f | g01 %.1f±%.2f | u03 %.1f±%.2f",
                sets$pic$mean, sets$pic$sd, sets$g01$mean, sets$g01$sd, sets$u03$mean, sets$u03$sd),
        side = 3, line = 0.1, cex = 0.55, col = "grey30")
  cat(sprintf("  %-16s %-7s | pic %7.2f ± %5.2f | g01 %7.2f ± %5.2f | u03 %7.2f ± %5.2f\n",
              model, vn, sets$pic$mean, sets$pic$sd, sets$g01$mean, sets$g01$sd, sets$u03$mean, sets$u03$sd))
}
for (k in seq_along(CCOL)) mtext(names(CCOL)[k], col = CCOL[k], side = 3, line = 0.6,
                                 at = seq(0.35, 0.65, length.out = 3)[k], outer = TRUE, cex = 0.9, font = 2)
mtext("QC: absolute tas / tasmin / tasmax [K] -- piControl baseline vs FULL hosing runs (all Europe cell-months, area-weighted)",
      outer = TRUE, line = 2.2, cex = 1.0, font = 2)
mtext("look for: mean shift (cooling), sd change (variance), spikes/clipping/bumps = weirdness; EC tasmin/tasmax are the reconstructed files",
      outer = TRUE, line = 1.0, cex = 0.72, col = "grey30")
dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
