# ============================================================
# ANIMATION FRAMES: year-by-year warming/cooling + drying/wetting.
#
# For each model x protocol run, renders one PNG per model year:
#   left  = ANNUAL-mean tas anomaly [K]      (RdBu, blue = cooling)
#   right = ANNUAL-mean pr  anomaly [mm/day] (BrBG, brown = drying)
# Colour scales are FIXED across the whole run (98th pct of |value|
# pooled over all years, clipped at the ends), so the animation shows
# real growth of the signal, not a rescaling artefact. Title carries
# the year counter and the Europe-mean values of that year.
#
# Frames go to the directory given as the first CLI argument
# (one subdir naming pattern <model>_<proto>_%03d.png); the caller
# assembles them into GIFs (ImageMagick), e.g.:
#   magick -delay 20 -loop 0 frames/EC-Earth3_u03_*.png anim.gif
# ============================================================

library(ncdf4)
library(fields)
have_maps <- requireNamespace("maps", quietly = TRUE)

## rename-proof: source amoc_common.R from THIS script's own folder
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: read_europe_cube, pal_for, wmean, CLIP_Q.

frame_dir <- commandArgs(TRUE)[1]
stopifnot(!is.na(frame_dir)); dir.create(frame_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list(
  list(model = "EC-Earth3",       proto = "g01", folder = "ECHearth3_anomaly"),
  list(model = "EC-Earth3",       proto = "u03", folder = "ECHearth3_anomaly"),
  list(model = "HadGEM3-GC31-LL", proto = "g01", folder = "LL_anomaly"),
  list(model = "HadGEM3-GC31-LL", proto = "u03", folder = "LL_anomaly"),
  list(model = "HadGEM3-GC31-MM", proto = "g01", folder = "MM_anomaly"),
  list(model = "HadGEM3-GC31-MM", proto = "u03", folder = "MM_anomaly"),
  list(model = "IPSL-CM6A-LR",    proto = "u03", folder = "IPSL_anomaly")   # u03 only (no g01)
)
base <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/anomaly_output"
anomaly_path <- function(r, vn) file.path(base, r$folder,
  if (r$folder %in% c("ECHearth3_anomaly", "IPSL_anomaly")) sprintf("%s_Amon_%s_hos-%s-hos_anomaly.nc", vn, r$model, r$proto)
  else sprintf("%s_anomaly_%s-hos_minus_piControl_1850-1949.nc", vn, r$proto))

annual <- function(v) {
  d <- dim(v); ny <- d[3] %/% 12
  apply(array(v[, , 1:(12 * ny), drop = FALSE], c(d[1], d[2], 12, ny)), c(1, 2, 4), mean, na.rm = TRUE)
}

for (r in runs) {
  pt <- anomaly_path(r, "tas"); pp <- anomaly_path(r, "pr")
  if (!file.exists(pt) || !file.exists(pp)) { cat("[skip]", r$model, r$proto, "\n"); next }
  ct <- read_europe_cube(pt); cp <- read_europe_cube(pp)
  At <- annual(ct$v); Ap <- annual(cp$v); ny <- dim(At)[3]
  zt <- as.numeric(quantile(abs(At), CLIP_Q, na.rm = TRUE))
  zp <- as.numeric(quantile(abs(Ap), CLIP_Q, na.rm = TRUE))
  key <- sprintf("%s_%s", r$model, r$proto)
  cat(sprintf("[frames] %s: %d years, scales ±%.2f K / ±%.2f mm/day\n", key, ny, zt, zp))
  for (y in seq_len(ny)) {
    png(file.path(frame_dir, sprintf("%s_%03d.png", key, y)), width = 900, height = 430, res = 100)
    par(mfrow = c(1, 2), oma = c(0, 0, 3.2, 0), mar = c(3, 3, 1.6, 4.4), mgp = c(1.8, 0.5, 0))
    image.plot(ct$lon, ct$lat, pmin(pmax(At[, , y], -zt), zt), zlim = c(-zt, zt),
               col = pal_for("tas"), xlab = "lon", ylab = "lat",
               main = "tas anomaly [K]  (blue = cooling)", cex.main = 0.85,
               legend.lab = "K", legend.line = 2.2)
    if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
    image.plot(cp$lon, cp$lat, pmin(pmax(Ap[, , y], -zp), zp), zlim = c(-zp, zp),
               col = pal_for("pr"), xlab = "lon", ylab = "lat",
               main = "pr anomaly [mm/day]  (brown = drying)", cex.main = 0.85,
               legend.lab = "mm/day", legend.line = 2.2)
    if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
    mtext(sprintf("%s %s — year %d of %d    (Europe mean: %+.2f K, %+.2f mm/day)",
                  r$model, r$proto, y, ny,
                  wmean(At[, , y], ct$lon, ct$lat), wmean(Ap[, , y], cp$lon, cp$lat)),
          outer = TRUE, line = 1.4, font = 2, cex = 1.0)
    mtext(sprintf("annual means; scales fixed for the whole run, clipped at the %.0fth pct: ±%.2f K / ±%.2f mm/day",
                  CLIP_Q * 100, zt, zp), outer = TRUE, line = 0.2, cex = 0.7, col = "grey30")
    dev.off()
    if (y %% 25 == 0) cat("   ...", y, "frames\n")
  }
}
cat("FRAMES ->", frame_dir, "\n")
