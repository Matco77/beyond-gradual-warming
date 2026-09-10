# ============================================================
# SEASONAL delta fields + per-cell SIGNIFICANCE, NATIVE grids.
#
# For DJF and JJA, per model x variable, on EACH MODEL'S OWN grid
# (Europe window cropped, NO interpolation):
#   * delta  = late-run (last third) seasonal mean anomaly, per cell;
#   * snr    = delta / (that cell's control spread of multidecadal
#              seasonal means) = a per-cell EFFECT SIZE. |snr| >= 2
#              outlines where the seasonal cooling is robust against
#              LOCAL natural variability. This is an effect size, not
#              a per-cell p-value: it avoids field-significance / FDR
#              and the short-control tail problem.
#
# Significance is computed on PER-YEAR seasonal means (within-season
# weather averaged out) and matched to the delta's averaging length.
# CAVEAT: the control sd uses ALL archived control years (MAX_PIC_FILES = Inf):
# EC-Earth3 501 yr and HadGEM-LL 2000 yr pin sigma well; HadGEM-MM 500 yr gives
# ~15 independent 33-yr windows per cell -- adequate, contour still indicative.
#
# CALENDAR: monthly files start in JANUARY (CMIP Amon standard).
# Reuses the v8 Europe crop. Native grids only. No nested indexing
# beyond clean per-file / per-row maps.
# ============================================================

library(ncdf4)
library(fields)
have_maps <- requireNamespace("maps", quietly = TRUE)

## rename-proof: source amoc_common.R from THIS script's own folder
## (keeps working if the folder is moved/renamed; run via `Rscript "<this file>"`).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)   # Rscript encodes spaces in --file= as ~+~
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: europe, MAX_PIC_FILES, EFFECT_K, CLIP_Q, publish_file,
## panel_grid, pal_for, pic_dir_for, and read_europe_cube (used below as the
## native-grid reader).

## Files verified to start in January (CMIP Amon). NOTE: "DJF" is the
## CALENDAR-year variant (Jan, Feb, Dec of the same year); it is applied
## identically to hosing and control, so the effect size stays consistent.
season_months <- list(DJF = c(12, 1, 2), JJA = c(6, 7, 8))

## ---- native-grid Europe cube: identical to amoc_common's read_europe_cube ----
## (despite the _raw name it inherits the pr kg m-2 s-1 -> mm/day conversion
## done once in amoc_common.R, so the mm/day labels below are correct)
read_europe_cube_raw <- read_europe_cube

## ---- monthly cube -> per-year seasonal-mean cube (nlon x nlat x ny) ----
yearly_season <- function(cube, season) {
  d <- dim(cube); ny <- d[3] %/% 12
  arr <- array(cube[, , 1:(12 * ny), drop = FALSE], c(d[1], d[2], 12, ny))
  apply(arr[, , season_months[[season]], , drop = FALSE], c(1, 2, 4), mean, na.rm = TRUE)
}

## ---- per-row (per-cell) linear detrend ----
detrend_rows <- function(M) {
  T <- ncol(M); t <- seq_len(T); tc <- matrix(t - mean(t), ncol = 1)
  b <- as.vector((M %*% tc) / sum(tc^2))
  M - (outer(b, t) + matrix(rowMeans(M) - b * mean(t), nrow(M), T))
}

## ---- per-row sd of length-w running means ----
roll_mean_sd <- function(M, w) {
  T <- ncol(M); if (T < w) return(rep(NA_real_, nrow(M)))
  cs <- cbind(0, t(apply(M, 1, cumsum)))
  wm <- (cs[, (w + 1):(T + 1), drop = FALSE] - cs[, 1:(T - w + 1), drop = FALSE]) / w
  sqrt(pmax(rowMeans(wm^2) - rowMeans(wm)^2, 0))
}

## ---- piControl dir + control per-year seasonal cube (native grid) ----
## Reshape each control file to (cells x years) and bind along time, so it
## works whether a file holds 1 year or many (HadGEM files bundle years).
control_season_cube <- function(p, vn, season) {
  pf  <- sort(list.files(pic_dir_for(p), pattern = paste0("^", vn, "_Amon_"), full.names = TRUE))
  pf  <- head(pf, MAX_PIC_FILES)
  cb1 <- read_europe_cube_raw(pf[1])
  cols <- lapply(pf, function(f) {
    ys <- yearly_season(read_europe_cube_raw(f)$v, season)   # nlon x nlat x ny_file
    matrix(ys, dim(ys)[1] * dim(ys)[2], dim(ys)[3])          # cells x ny_file
  })
  list(lon = cb1$lon, lat = cb1$lat, MC = do.call(cbind, cols))   # cells x total_years
}

## ---- runs + paths ----
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

vars     <- c("tas", "pr", "tasmax", "tasmin")
seas_dir <- file.path(dirname(base), "delta_fields", "seasonal")
outdir   <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly"
dir.create(seas_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir,   showWarnings = FALSE, recursive = TRUE)

write_seasonal_nc <- function(out, lon, lat, delta, snr, vn) {
  ud <- if (vn == "pr") "mm/day" else "K"
  dl <- ncdim_def("lon", "degrees_east", lon); da <- ncdim_def("lat", "degrees_north", lat)
  vd <- ncvar_def(vn, ud, list(dl, da), -9999, paste(vn, "seasonal late-run mean anomaly"))
  vs <- ncvar_def(paste0(vn, "_snr"), "ratio", list(dl, da), -9999, "delta / control multidecadal seasonal-mean sd")
  tmp <- tempfile(fileext = ".nc")                 # build off OneDrive, then publish
  nc <- nc_create(tmp, list(vd, vs))
  d2 <- delta; d2[is.na(d2)] <- -9999; s2 <- snr; s2[is.na(s2)] <- -9999
  ncvar_put(nc, vd, d2); ncvar_put(nc, vs, s2); nc_close(nc)
  publish_file(tmp, out)
}

## ---- compute one (run, var, season): delta + snr on native grid ----
compute_seasonal <- function(r, vn, season) {
  p <- anomaly_path(r, vn); if (!file.exists(p)) return(NULL)
  ac  <- read_europe_cube_raw(p)
  A   <- yearly_season(ac$v, season); nyA <- dim(A)[3]; plY <- max(5, round(nyA / 3))
  MA  <- matrix(A, dim(A)[1] * dim(A)[2], nyA)
  dvec <- rowMeans(MA[, (nyA - plY + 1):nyA, drop = FALSE], na.rm = TRUE)

  cc  <- control_season_cube(p, vn, season)
  stopifnot(length(cc$lon) == length(ac$lon), length(cc$lat) == length(ac$lat),
            nrow(cc$MC) == length(dvec))
  svec <- roll_mean_sd(detrend_rows(cc$MC), plY)

  delta <- matrix(dvec, length(ac$lon), length(ac$lat)); delta[is.nan(delta)] <- NA
  snr   <- matrix(dvec / svec, length(ac$lon), length(ac$lat)); snr[is.nan(snr)] <- NA

  out <- file.path(seas_dir, sprintf("seasonal_%s_%s_%s_%s_lastthird.nc", r$model, r$proto, vn, season))
  write_seasonal_nc(out, ac$lon, ac$lat, delta, snr, vn)

  ## robust fraction AREA-WEIGHTED (cos-lat) so it really is "% of Europe",
  ## not "% of grid cells" (equal counting overweights high latitudes ~2.7x)
  W  <- matrix(cos(ac$lat * pi / 180), length(ac$lon), length(ac$lat), byrow = TRUE)
  ok <- is.finite(snr)
  robust <- if (any(ok)) sum(W[ok] * (abs(snr[ok]) >= EFFECT_K)) / sum(W[ok]) else NA_real_
  cat(sprintf("  %-16s %s %-6s %s | mean delta %+.2f | robust cells %4.0f%%\n",
              r$model, r$proto, vn, season,
              sum(delta * W, na.rm = TRUE) / sum(W * !is.na(delta)), 100 * robust))
  list(r = r, lon = ac$lon, lat = ac$lat, delta = delta, snr = snr, robust = robust)
}

# ============================================================
# RUN: produce NetCDFs for both seasons, then panel them.
# ============================================================

## ---- render the robust-% summary as a table page ----
render_robust_table <- function(S) {
  S$pct <- round(100 * S$robust)
  w <- reshape(S[, c("model", "proto", "variable", "season", "pct")],
               idvar = c("model", "proto", "variable"), timevar = "season", direction = "wide")
  ord <- order(match(w$variable, vars),
               match(paste(w$model, w$proto), vapply(runs, function(r) paste(r$model, r$proto), character(1))))
  w <- w[ord, ]
  plot.new(); plot.window(c(0, 1), c(0, 1))
  title("Robust-cell fraction: % of Europe (area-weighted) with |effect| >= 2 per cell")
  cx <- c(0.03, 0.42, 0.64, 0.80); y0 <- 0.93; dy <- min(0.045, 0.84 / (nrow(w) + 1))
  text(cx, y0, c("Model / protocol", "Variable", "DJF %", "JJA %"), font = 2, adj = 0, cex = 0.8)
  for (i in seq_len(nrow(w))) {
    yy <- y0 - i * dy
    text(cx[1], yy, sprintf("%s %s", w$model[i], w$proto[i]), adj = 0, cex = 0.75)
    text(cx[2], yy, w$variable[i], adj = 0, cex = 0.75)
    text(cx[3], yy, w[i, "pct.DJF"], adj = 0, cex = 0.75)
    text(cx[4], yy, w[i, "pct.JJA"], adj = 0, cex = 0.75)
  }
}

SUMMARY <- data.frame()

pdf_final <- file.path(outdir, "seasonal_fields_significance.pdf")
pdf_tmp   <- tempfile(fileext = ".pdf")        # render locally, publish to OneDrive at the end
pdf(pdf_tmp, width = 11, height = 8)
for (vn in vars) for (season in c("DJF", "JJA")) {
  cat(sprintf("== %s %s ==\n", vn, season))
  F <- Filter(Negate(is.null), lapply(runs, compute_seasonal, vn = vn, season = season))
  if (!length(F)) next
  for (d in F) SUMMARY <- rbind(SUMMARY, data.frame(model = d$r$model, proto = d$r$proto,
                                                    variable = vn, season = season, robust = d$robust))

  zmax <- as.numeric(quantile(abs(unlist(lapply(F, function(d) as.vector(d$delta)))), CLIP_Q, na.rm = TRUE))
  if (!is.finite(zmax) || zmax == 0) zmax <- 1
  pal <- pal_for(vn); unit <- if (vn == "pr") "mm/day" else "K"

  par(mfrow = panel_grid(length(F)), oma = c(0, 0, 3, 0), mar = c(3.2, 3.2, 2.4, 4.2), mgp = c(2, 0.6, 0))
  for (d in F) {
    zr <- pmin(pmax(d$delta, -zmax), zmax)
    zr[abs(d$snr) < EFFECT_K | is.na(d$snr)] <- NA          # keep ONLY robust cells in colour
    image(d$lon, d$lat, matrix(0, length(d$lon), length(d$lat)), col = "grey85", zlim = c(0, 1),
          xlab = "lon", ylab = "lat", main = sprintf("%s %s", d$r$model, d$r$proto), cex.main = 0.92)
    image.plot(d$lon, d$lat, zr, zlim = c(-zmax, zmax), col = pal, add = TRUE,
               legend.lab = unit, legend.line = 2.3)
    if (have_maps) maps::map("world", add = TRUE, interior = FALSE, lwd = 0.5)
  }
  mtext(sprintf("%s seasonal delta - %s [%s]   |   scale +/-%.2f (98th pct), coloured = robust (|effect|>=2), grey = within noise",
                season, vn, unit, zmax), outer = TRUE, cex = 0.92, font = 2)
}
if (nrow(SUMMARY)) render_robust_table(SUMMARY)
dev.off()
publish_file(pdf_tmp, pdf_final)

if (nrow(SUMMARY)) {
  SUMMARY$robust_pct <- round(100 * SUMMARY$robust)
  wide <- reshape(SUMMARY[, c("model", "proto", "variable", "season", "robust_pct")],
                  idvar = c("model", "proto", "variable"), timevar = "season", direction = "wide")
  names(wide) <- sub("robust_pct\\.", "robust_pct_", names(wide))
  csv_path <- file.path(seas_dir, "seasonal_robust_summary.csv")
  csv_tmp  <- tempfile(fileext = ".csv")
  write.csv(wide, csv_tmp, row.names = FALSE)
  publish_file(csv_tmp, csv_path)
  cat("Robust-% summary ->", csv_path, "\n")
}
cat("\nSeasonal fields ->", seas_dir, "\nSeasonal PDF   ->", pdf_final, "\nDONE\n")
