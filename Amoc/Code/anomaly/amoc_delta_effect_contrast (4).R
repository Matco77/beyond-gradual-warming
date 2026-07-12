# ============================================================
# AMOC-hosing: DELTA fields + AMOC-EFFECT size + MODEL CONTRAST
#
# Three outputs, nothing more:
#  (1) EFFECT?  plateau / sd(control plateau-length means).
#               Reported as an effect SIZE (ratio), not a tail
#               p-value, so the short HadGEM control (100 yr) is
#               NOT a degrees-of-freedom problem. ratio >= EFFECT_K
#               => the imposed AMOC weakening has a European effect;
#               ratio < EFFECT_K => inside control spread (no effect
#               yet, e.g. EC g01: weak 0.1 Sv, 50 yr, not settled).
#  (2) FINGERPRINT: NW-minus-Med gradient. AMOC cooling is
#               NW-amplified, so a real effect has NW more negative
#               than the Mediterranean (grad < 0). A uniform offset
#               (drift) is not the AMOC. The gradient also differs
#               across models -> the model-difference story.
#  (3) DELTA:   per-cell LATE-RUN (plateau-window) mean anomaly
#               field, written to NetCDF per model x variable on the
#               Europe window = the matured weak-AMOC climate state
#               to drive the European impact variables.
#
# Headline is the LEVEL (plateau), never a growing slope: constant-
# step hosing -> NEW EQUILIBRIUM (ramp-then-plateau), so a flat
# anomaly in the late run is the CORRECT, expected result. A model
# still cooling at the end (HadGEM-MM) -> its delta is a LOWER bound.
#
# Reuses the v8 readers verbatim (same Europe window, same cos-lat
# weighting, same de-drifted piControl). One pass, ONE top-level map
# over the SAME 24 anomaly files (EC-Earth3 g01+u03, HadGEM-LL g01+u03,
# HadGEM-MM g01+u03). No nested loops.
# ============================================================

library(ncdf4)   # Pierce, D. (2023). ncdf4: Interface to Unidata netCDF.
library(fields)  # Nychka et al. (2021). fields: Tools for Spatial Data.

## rename-proof: source amoc_common.R from THIS script's own folder
## (keeps working if the folder is moved/renamed; run via `Rscript "<this file>"`).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.self <- gsub("~+~", " ", .self, fixed = TRUE)   # Rscript encodes spaces in --file= as ~+~
source(file.path(if (length(.self)) dirname(normalizePath(.self)) else getwd(), "amoc_common.R"))
## amoc_common.R provides: europe, EFFECT_K, nw_box, med_box, publish_file,
## pic_dir_for, read_europe_cube, europe_monthly_means, mean_null, box_mean, and
## get_control_series (now deseasonalised + de-drifted -> same ratio as the plots).

## ---- write a 2D delta field to NetCDF (built locally, then published) ----
write_delta_nc <- function(out_path, lon, lat, fld, vn) {
  units    <- if (vn == "pr") "mm/day" else "K"
  longname <- paste0(vn, " late-run (plateau-window) mean anomaly")
  dlon <- ncdim_def("lon", "degrees_east",  lon)
  dlat <- ncdim_def("lat", "degrees_north", lat)
  vdef <- ncvar_def(vn, units, list(dlon, dlat), missval = -9999, longname = longname)
  tmp  <- tempfile(fileext = ".nc")          # build off OneDrive to dodge sync locks
  nc   <- nc_create(tmp, vdef)
  fld2 <- fld; fld2[is.na(fld2)] <- -9999
  ncvar_put(nc, vdef, fld2)
  nc_close(nc)
  publish_file(tmp, out_path)
}

## ---- ONE anomaly file -> stats row (+ writes the delta nc) ----
process_one <- function(path, delta_dir) {
  model <- switch(basename(dirname(path)),
                  "ECHearth3_anomaly" = "EC-Earth3",
                  "LL_anomaly"        = "HadGEM3-GC31-LL",
                  "MM_anomaly"        = "HadGEM3-GC31-MM")
  proto <- if (grepl("u03", basename(path))) "u03" else "g01"
  vn    <- strsplit(basename(path), "_")[[1]][1]

  ## Europe-mean monthly anomaly -> settled PLATEAU (last third)
  ts      <- europe_monthly_means(path)
  n       <- length(ts); pl <- max(12, round(n / 3))
  plateau <- mean(tail(ts, pl))
  run_yr  <- round(n / 12)

  ## EFFECT SIZE = plateau / control plateau-length-mean sd (matched timescale)
  ctrl   <- get_control_series(path); ctrl_a <- ctrl - mean(ctrl)
  sd_wm  <- sd(mean_null(ctrl_a, pl))
  ratio  <- abs(plateau) / sd_wm
  effect <- if (is.finite(ratio) && ratio >= EFFECT_K) "YES" else "no"

  ## per-cell LATE-RUN mean -> the delta field (matured weak-AMOC state)
  cb   <- read_europe_cube(path)
  dv   <- dim(cb$v); M <- matrix(cb$v, dv[1] * dv[2], dv[3])
  late <- rowMeans(M[, (n - pl + 1):n, drop = FALSE], na.rm = TRUE)
  fld  <- matrix(late, dv[1], dv[2]); fld[is.nan(fld)] <- NA

  ## FINGERPRINT: NW minus Mediterranean (AMOC cooling => grad < 0)
  grad <- box_mean(fld, cb$lon, cb$lat, nw_box) - box_mean(fld, cb$lon, cb$lat, med_box)

  out_nc <- file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, proto, vn))
  write_delta_nc(out_nc, cb$lon, cb$lat, fld, vn)

  cat(sprintf("  %-16s %s %-6s | run=%3dyr  plateau=%+.3f  sd_win=%.3f  ratio=%5.2f  effect=%-3s  NW-Med=%+.3f\n",
              model, proto, vn, run_yr, plateau, sd_wm, ratio, effect, grad))

  data.frame(model = model, protocol = proto, variable = vn, run_years = run_yr,
             plateau = round(plateau, 3), sd_window = round(sd_wm, 3),
             effect_ratio = round(ratio, 2), effect = effect,
             NW_minus_Med = round(grad, 3), delta_file = basename(out_nc),
             stringsAsFactors = FALSE)
}

# ============================================================
# RUN: the 24 anomaly files; delta NetCDFs + one contrast CSV.
# (Any anomaly not built yet -- e.g. a u03-hos run still missing for a model --
#  is skipped with a notice instead of aborting the whole contrast.)
# ============================================================
base      <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/anomaly_output"
delta_dir <- file.path(dirname(base), "delta_fields")   # .../datasets/delta_fields
dir.create(delta_dir, showWarnings = FALSE, recursive = TRUE)

files <- c(
  file.path(base, "ECHearth3_anomaly/pr_Amon_EC-Earth3_hos-g01-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/pr_Amon_EC-Earth3_hos-u03-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tas_Amon_EC-Earth3_hos-g01-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tas_Amon_EC-Earth3_hos-u03-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tasmax_Amon_EC-Earth3_hos-g01-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tasmax_Amon_EC-Earth3_hos-u03-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tasmin_Amon_EC-Earth3_hos-g01-hos_anomaly.nc"),
  file.path(base, "ECHearth3_anomaly/tasmin_Amon_EC-Earth3_hos-u03-hos_anomaly.nc"),
  file.path(base, "LL_anomaly/pr_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tas_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tasmax_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tasmin_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/pr_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tas_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tasmax_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "LL_anomaly/tasmin_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/pr_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tas_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tasmax_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tasmin_anomaly_g01-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/pr_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tas_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tasmax_anomaly_u03-hos_minus_piControl_1850-1949.nc"),
  file.path(base, "MM_anomaly/tasmin_anomaly_u03-hos_minus_piControl_1850-1949.nc")
)

rows    <- lapply(files, function(p) {
  if (!file.exists(p)) { cat("MISSING (skipped):", basename(p), "\n"); return(NULL) }
  cat("processing:", basename(p), "\n"); process_one(p, delta_dir)
})
rows    <- Filter(Negate(is.null), rows)          # drop not-yet-built anomalies
contrast <- do.call(rbind, rows)
contrast <- contrast[order(contrast$variable, contrast$protocol, contrast$model), ]

csv_path <- file.path(delta_dir, "amoc_effect_model_contrast.csv")
csv_tmp  <- tempfile(fileext = ".csv")
write.csv(contrast, csv_tmp, row.names = FALSE)
publish_file(csv_tmp, csv_path)

cat("\n================ AMOC EFFECT + MODEL CONTRAST ================\n")
print(contrast, row.names = FALSE)
cat("\n  effect = YES when |plateau| >=", EFFECT_K, "x control plateau-length-mean sd\n")
cat("  NW_minus_Med < 0  => NW-amplified cooling = AMOC fingerprint (not a uniform offset)\n")
cat("  NOTE: HadGEM-MM is still cooling at year 100 -> its delta field is a LOWER bound.\n")
cat("\nDelta fields ->", delta_dir, "\nContrast CSV ->", csv_path, "\nALL DONE\n")

## ---- to eyeball any delta field afterwards (two lines) ----
# d <- nc_open(file.path(delta_dir, "delta_EC-Earth3_u03_tas_lastthird.nc"))
# image.plot(ncvar_get(d,"lon"), ncvar_get(d,"lat"), ncvar_get(d,"tas")); nc_close(d)
