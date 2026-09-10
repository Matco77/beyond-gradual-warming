# ============================================================
# QUADRANT view (wet/dry x hot/cold) + TIME-SPACE evolution.
#
# Six pages, one panel per model x protocol, Europe window, native
# grids (no interpolation):
#  (1-3) MAPS of the FIRST / SECOND / FINAL third of the run: each cell
#      classified by the SIGN of its mean tas and pr anomaly over that
#      (year-based) third ->
#         cold+dry / cold+wet / hot+dry / hot+wet
#      (hot: dtas >= 0; wet: dpr >= 0). Panel title reports the year
#      window and the area-weighted % of each class. Sign-only view:
#      magnitude and robustness are NOT encoded here (per-cell effect
#      sizes live in seasonal_fields_significance). The final third is
#      computed from whole years, so it matches the month-based
#      delta_*_lastthird exports to within a few months.
#  (4) EVOLUTION of those class shares through the run: annual-mean
#      anomalies classified per year -> stacked area summing to 100%.
#  (5) HOVMOELLER, tas: zonal-mean (over the window's longitudes)
#      ANNUAL anomaly, latitude x time -- the space-time fingerprint of
#      the AMOC response deepening/spreading by latitude.
#  (6) Same Hovmoeller for pr (mm/day).
#
# Colours (CVD-validated; hue = temperature, darkness = wetness):
#   cold+dry #4393C3   cold+wet #2166AC   hot+dry #D6604D   hot+wet #B2182B
# Inputs: the full monthly anomaly files (all pages). pr arrives in
# mm/day via amoc_common's reader.
# ============================================================

library(ncdf4)
library(fields)
have_maps <- requireNamespace("maps", quietly = TRUE)

## rename-proof: source amoc_common.R from THIS script's own folder
## (keeps working if the folder is moved/renamed; run via `Rscript "<this file>"`).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)   # Rscript encodes spaces in --file= as ~+~
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: europe, CLIP_Q, publish_file, panel_grid, pal_for,
## read_europe_cube.

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
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

anomaly_path <- function(r, vn) file.path(base, r$folder,
  if (r$folder %in% c("ECHearth3_anomaly", "IPSL_anomaly")) sprintf("%s_Amon_%s_hos-%s-hos_anomaly.nc", vn, r$model, r$proto)
  else sprintf("%s_anomaly_%s-hos_minus_piControl_1850-1949.nc", vn, r$proto))

QCOL <- c("#4393C3", "#2166AC", "#D6604D", "#B2182B")     # validated, see header
QLAB <- c("cold+dry", "cold+wet", "hot+dry", "hot+wet")

classify <- function(t, p) {              # 1 cold+dry, 2 cold+wet, 3 hot+dry, 4 hot+wet
  cl <- ifelse(t >= 0, 3, 1) + ifelse(p >= 0, 1, 0)
  cl[!is.finite(t) | !is.finite(p)] <- NA
  cl
}
shares <- function(cl, lon, lat) {        # area-weighted % of window per class
  W  <- matrix(cos(lat * pi / 180), length(lon), length(lat), byrow = TRUE)
  ok <- !is.na(cl); tot <- sum(W[ok])
  vapply(1:4, function(k) 100 * sum(W[ok & cl == k]) / tot, numeric(1))
}
qlegend <- function() for (k in 1:4)
  mtext(QLAB[k], col = QCOL[k], side = 3, line = 0.6,
        at = seq(0.14, 0.86, length.out = 4)[k], outer = TRUE, cex = 0.8, font = 2)

annual <- function(v) {                   # monthly cube -> annual-mean cube
  d <- dim(v); ny <- d[3] %/% 12
  apply(array(v[, , 1:(12 * ny), drop = FALSE], c(d[1], d[2], 12, ny)), c(1, 2, 4), mean, na.rm = TRUE)
}

pdf_final <- file.path(outdir, "quadrants_evolution.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")   # render locally, publish to OneDrive at the end
pdf(pdf_tmp, width = 11, height = 8)

## ---------- load annual anomaly cubes once (all pages) ----------
AR <- Filter(Negate(is.null), lapply(runs, function(r) {
  pt <- anomaly_path(r, "tas"); pp <- anomaly_path(r, "pr")
  if (!file.exists(pt) || !file.exists(pp)) return(NULL)
  ct <- read_europe_cube(pt); cp <- read_europe_cube(pp)
  cat("  [cubes]", r$model, r$proto, "annualised\n")
  list(r = r, lon = ct$lon, lat = ct$lat, At = annual(ct$v), Ap = annual(cp$v))
}))

## ---------- PAGES 1-3: quadrant maps, FIRST / SECOND / FINAL third ----------
third_idx <- function(ny, t) { k <- ny %/% 3; if (t < 3) ((t - 1) * k + 1):(t * k) else (2 * k + 1):ny }
TLAB <- c("First-third (yr 1..n/3)", "Second-third", "Final-third")
for (t3 in 1:3) {
  P <- lapply(AR, function(d) {
    idx <- third_idx(dim(d$At)[3], t3)
    Tm  <- apply(d$At[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    Pm  <- apply(d$Ap[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    cl  <- classify(Tm, Pm); s <- shares(cl, d$lon, d$lat)
    cat(sprintf("  %-16s %s %-14s | CW %2.0f%%  CD %2.0f%%  HD %2.0f%%  HW %2.0f%%\n",
                d$r$model, d$r$proto, TLAB[t3], s[2], s[1], s[3], s[4]))
    list(r = d$r, lon = d$lon, lat = d$lat, cl = cl, s = s, yrs = range(idx))
  })
  par(mfrow = panel_grid(length(P)), oma = c(0, 0, 4, 0), mar = c(3.2, 3.2, 3.0, 1), mgp = c(2, 0.6, 0))
  for (d in P) {
    image(d$lon, d$lat, matrix(0, length(d$lon), length(d$lat)), col = "grey85", zlim = c(0, 1),
          xlab = "lon", ylab = "lat", cex.main = 0.85,
          main = sprintf("%s %s  (yr %d–%d)\nCW %.0f%%  CD %.0f%%  HD %.0f%%  HW %.0f%%",
                         d$r$model, d$r$proto, d$yrs[1], d$yrs[2], d$s[2], d$s[1], d$s[3], d$s[4]))
    image(d$lon, d$lat, d$cl, col = QCOL, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5), add = TRUE)
    if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
  }
  qlegend()
  mtext(sprintf("%s state, SIGN of tas & pr delta per cell  |  %% = area-weighted share of the Europe window",
                TLAB[t3]),
        outer = TRUE, line = 2.2, cex = 0.95, font = 2)
}

## ---------- PAGE 4: evolution of the class shares ----------
par(mfrow = panel_grid(length(AR)), oma = c(0, 0, 4, 0), mar = c(3.2, 3.4, 2, 1), mgp = c(2, 0.6, 0))
for (d in AR) {
  ny <- dim(d$At)[3]; x <- seq_len(ny)
  SH  <- t(vapply(x, function(y) shares(classify(d$At[, , y], d$Ap[, , y]), d$lon, d$lat), numeric(4)))
  cum <- cbind(0, t(apply(SH, 1, cumsum)))          # ny x 5 cumulative bounds
  plot(NA, xlim = c(1, ny), ylim = c(0, 100), xlab = "Years from start of hosing",
       ylab = "% of window area", main = sprintf("%s %s", d$r$model, d$r$proto), cex.main = 0.9)
  for (k in 1:4) polygon(c(x, rev(x)), c(cum[, k], rev(cum[, k + 1])),
                         col = QCOL[k], border = "white", lwd = 0.6)
}
qlegend()
mtext("Evolution of the quadrant shares (annual-mean anomalies, sign classification, area-weighted)",
      outer = TRUE, line = 2.2, cex = 0.95, font = 2)

## ---------- PAGES 5-6: Hovmoeller latitude x time ----------
hov_page <- function(vn, unit) {
  H <- lapply(AR, function(d) {
    A <- if (vn == "tas") d$At else d$Ap
    list(r = d$r, lat = d$lat, h = apply(A, c(2, 3), mean, na.rm = TRUE))   # lat x year
  })
  zmax <- as.numeric(quantile(abs(unlist(lapply(H, function(z) as.vector(z$h)))), CLIP_Q, na.rm = TRUE))
  if (!is.finite(zmax) || zmax == 0) zmax <- 1
  par(mfrow = panel_grid(length(H)), oma = c(0, 0, 3, 0), mar = c(3.2, 3.4, 2, 4.2), mgp = c(2, 0.6, 0))
  for (z in H) {
    ny <- ncol(z$h)
    image.plot(seq_len(ny), z$lat, t(pmin(pmax(z$h, -zmax), zmax)), zlim = c(-zmax, zmax),
               col = pal_for(vn), xlab = "Years from start of hosing", ylab = "lat",
               main = sprintf("%s %s", z$r$model, z$r$proto), cex.main = 0.9,
               legend.lab = unit, legend.line = 2.3)
  }
  mtext(sprintf("Hovmoeller: zonal-mean (%g..%gE) ANNUAL %s anomaly, latitude x time  |  scale clipped at %.0fth pct = ±%.2f %s",
                europe$lon[1], europe$lon[2], vn, CLIP_Q * 100, zmax, unit),
        outer = TRUE, cex = 0.92, font = 2)
}
hov_page("tas", "K")
hov_page("pr", "mm/day")

dev.off()
publish_file(pdf_tmp, pdf_final)
cat("PLOTTED ->", pdf_final, "\n")
