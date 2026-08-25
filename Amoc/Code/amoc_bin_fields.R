# Per-bin AMOC effect fields, MONTHLY, on the native model grid cropped to the E-OBS box.
# Author: Marco Bova
#
# WHY THIS EXISTS. amoc_effect_bins_<MODEL>.nc has dims (bin, lat, lon): the month axis is already
# collapsed, it carries only delta_tas and delta_pr, and delta_pr is an ADDITIVE mm/day anomaly.
# The scenario replay needs all three of tasmin/tasmax/tas (heat_daily reads tx, frost_daily reads
# tn, within_day needs the diurnal amplitude (tx-tn)/2) and a MULTIPLICATIVE precip factor, per
# CALENDAR MONTH. So the bin fields are rebuilt here from the monthly anomaly cubes that
# anomaly_output/ already holds, which have every year and every month.
#
# BINNING, reproduced from the existing files and verified against them. delta_Sv of hosing year y
# = M26 hos(y) - mean(M26 con) (u03-hos; M26 carries no Sv series for g01, so g01 cannot be binned).
# Bins are 1 Sv wide by floor(), and a bin is kept only with >= 3 years. The M26 hos series is
# longer than the anomaly cube for some models (LL: 145 valid years vs 100 in the cube), so only
# the first min(valid, cube) years are usable. With that rule this script reproduces the n_years
# vector of all four amoc_effect_bins_* files EXACTLY - which is what pins the year -> bin mapping.
#
# Output: amoc_bin_fields_<MODEL>_u03.nc, dims (lon, lat, month, bin)
#   delta_tas / delta_tasmin / delta_tasmax   K, additive, bin mean per calendar month
#   pr_ratio                                  dimensionless, hosing/piControl, bin mean per month
# The consumer applies the [0.1, 10] clamp (scenario_engine.R), so the same clamp lands on the bin
# means and on the per-year fields of the uncertainty band.
suppressMessages({library(ncdf4)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

MODELS    <- c("IPSL-CM6A-LR", "EC-Earth3", "HadGEM3-GC3-1LL", "HadGEM3-GC3-1MM")
VARS      <- c("tas", "tasmin", "tasmax", "pr")     # pr is read from the RATIO file, not the anomaly
MIN_YEARS <- 3L                                     # same rule as the existing bins files
# crop box = E-OBS extent (lon -40.375..75.375, lat 25.375..75.375) + 3 deg, enough margin for the
# bilinear interpolation onto E-OBS cell centres to always have four model neighbours.
LON0 <- -44; LON1 <- 79; LAT0 <- 22; LAT1 <- 79

apath <- function(model, var) {
  b <- file.path(d, "anomaly_output"); kind <- if (var == "pr") "ratio" else "anomaly"
  switch(model,
    "IPSL-CM6A-LR"    = file.path(b, "IPSL_anomaly",      sprintf("%s_Amon_IPSL-CM6A-LR_hos-u03-hos_%s.nc", var, kind)),
    "EC-Earth3"       = file.path(b, "ECHearth3_anomaly", sprintf("%s_Amon_EC-Earth3_hos-u03-hos_%s.nc",    var, kind)),
    "HadGEM3-GC3-1LL" = file.path(b, "LL_anomaly", if (var == "pr") "pr_ratio_u03-hos_over_piControl_1850-1949.nc"
                                                   else sprintf("%s_anomaly_u03-hos_minus_piControl_1850-1949.nc", var)),
    "HadGEM3-GC3-1MM" = file.path(b, "MM_anomaly", if (var == "pr") "pr_ratio_u03-hos_over_piControl_1850-1949.nc"
                                                   else sprintf("%s_anomaly_u03-hos_minus_piControl_1850-1949.nc", var)))
}

# (hosing year, calendar month) of every time step, decoded from the time VALUES - never from the
# position. Two calendars are in play (IPSL/EC gregorian, LL/MM 360_day, where a naive days-since
# decode is wrong by two centuries) and, more to the point, MM's tasmin/tasmax carry 1197 steps
# against 1200 for tas with the SAME first and last value: three months are missing mid-series, so
# index (y-1)*12+m would silently shift every month after the gap into the wrong month and bin.
time_map <- function(nc) {
  t <- nc$dim$time$vals; cal <- nc$dim$time$calendar
  if (!is.null(cal) && grepl("360", cal)) { k <- floor(t / 30); mo <- (k %% 12) + 1L; yr <- k %/% 12 }
  else { dt <- as.Date(t, origin = sub(" .*", "", sub(".*since *", "", nc$dim$time$units)))
         mo <- as.integer(format(dt, "%m")); yr <- as.integer(format(dt, "%Y")) }
  # hosing year index is 1-based from the first year on the axis; first/last value agree across the
  # variables of a model, so the indices align even when a variable is missing interior months.
  data.frame(i = seq_along(t), yi = as.integer(yr - min(yr) + 1L), mo = as.integer(mo))
}

## ---- delta_Sv per hosing year, and the bins ------------------------------------------------
sv_bins <- function(model, n_cube_years) {
  nc <- nc_open(file.path(d, sprintf("M26/M26_%s.nc", model)))
  msk <- function(x) { x[!is.finite(x) | abs(x) > 1e10] <- NA; x }   # M26 fill value ~9.97e36, no attribute
  hos <- msk(ncvar_get(nc, paste0("hos_", model))); con <- msk(ncvar_get(nc, paste0("con_", model)))
  nc_close(nc)
  nv <- sum(!is.na(hos)); stopifnot(identical(which(!is.na(hos)), seq_len(nv)))   # valid block is a prefix
  n  <- min(nv, n_cube_years)                                  # years BOTH M26 and the cube have
  dsv <- (hos - mean(con, na.rm = TRUE))[seq_len(n)]
  lo  <- floor(dsv); keep <- names(which(table(lo) >= MIN_YEARS))
  bins <- lapply(sort(as.integer(keep)), function(b)
    list(lo = b, hi = b + 1, years = which(lo == b), dsv = mean(dsv[lo == b])))
  list(bins = bins, n_used = n, dsv = dsv)
}

## ---- one model ------------------------------------------------------------------------------
build <- function(model) {
  nc0 <- nc_open(apath(model, "tas"))
  TM0 <- time_map(nc0); ncube <- max(TM0$yi)
  lon <- as.numeric(ncvar_get(nc0, "lon")); lat <- as.numeric(ncvar_get(nc0, "lat"))
  nc_close(nc0)
  B <- sv_bins(model, ncube); bins <- B$bins

  # lat is increasing in all four cubes -> the Europe band is one contiguous ncdf4 slab
  stopifnot(all(diff(lat) > 0))
  li <- which(lat >= LAT0 & lat <= LAT1); lat_e <- lat[li]
  # lon comes as 0..360; Europe straddles the seam, so convert to -180..180 and reorder in memory
  lon180 <- ((lon + 180) %% 360) - 180; ox <- order(lon180)
  oi <- ox[lon180[ox] >= LON0 & lon180[ox] <= LON1]; lon_e <- lon180[oi]
  stopifnot(all(diff(lon_e) > 0))                              # interp.surface needs increasing x

  nx <- length(lon_e); ny <- length(lat_e); nb <- length(bins)
  cat(sprintf("%-16s grid %dx%d -> Europe %dx%d | years used %d | bins %d (%s)\n",
              model, length(lon), length(lat), nx, ny, B$n_used, nb,
              paste(sapply(bins, function(b) sprintf("%d:%d", b$lo, length(b$years))), collapse = " ")))

  out <- list()
  for (v in VARS) {
    nc <- nc_open(apath(model, v)); TM <- time_map(nc)
    # whole Europe slab at once (<= ~150 MB per variable), then index by (year, month) in memory
    A <- ncvar_get(nc, v, start = c(1, li[1], 1), count = c(-1, ny, -1))[oi, , , drop = FALSE]
    nc_close(nc)
    m <- array(NA_real_, c(nx, ny, 12, nb)); nonfin <- 0L; gaps <- 0L
    for (bi in seq_len(nb)) for (mo in 1:12) {
      idx <- TM$i[TM$yi %in% bins[[bi]]$years & TM$mo == mo]
      gaps <- gaps + (length(bins[[bi]]$years) - length(idx))   # months this variable simply lacks
      if (!length(idx)) next
      S <- A[, , idx, drop = FALSE]
      f <- is.finite(S); nonfin <- nonfin + sum(!f); S[!f] <- 0 # Inf turns up in the pr ratio on dry cells
      cn <- rowSums(f, dims = 2); sm <- rowSums(S, dims = 2)
      m[, , mo, bi] <- ifelse(cn > 0, sm / cn, NA_real_)        # mean over the bin's years, that month
    }
    out[[if (v == "pr") "pr_ratio" else paste0("delta_", v)]] <- m
    cat(sprintf("   %-7s done | missing (year,month) slices %d | non-finite values %d (%.4f%%)\n",
                v, gaps, nonfin, 100 * nonfin / max(1, length(A))))
    rm(A); flush.console()
  }

  ## ---- write ---------------------------------------------------------------------------------
  dx <- ncdim_def("lon", "degrees_east",  lon_e); dy <- ncdim_def("lat", "degrees_north", lat_e)
  dm <- ncdim_def("month", "1", 1:12);            db <- ncdim_def("bin", "Sv", sapply(bins, function(b) b$lo + 0.5))
  un <- c(delta_tas = "K", delta_tasmin = "K", delta_tasmax = "K", pr_ratio = "1")
  vars <- lapply(names(out), function(k) ncvar_def(k, un[[k]], list(dx, dy, dm, db), NA, prec = "double"))
  meta <- list(ncvar_def("bin_lo", "Sv", db, NA, prec = "double"),
               ncvar_def("bin_hi", "Sv", db, NA, prec = "double"),
               ncvar_def("n_years", "1", db, NA, prec = "double"),
               ncvar_def("delta_sv_mean", "Sv", db, NA, prec = "double"))
  f <- file.path(d, sprintf("amoc_bin_fields_%s_u03.nc", model))
  nc <- nc_create(f, c(vars, meta))
  for (k in names(out)) ncvar_put(nc, k, out[[k]])
  ncvar_put(nc, "bin_lo",  sapply(bins, `[[`, "lo"))
  ncvar_put(nc, "bin_hi",  sapply(bins, `[[`, "hi"))
  ncvar_put(nc, "n_years", sapply(bins, function(b) length(b$years)))
  ncvar_put(nc, "delta_sv_mean", sapply(bins, `[[`, "dsv"))
  ncatt_put(nc, 0, "model", model)
  ncatt_put(nc, 0, "experiment", "u03-hos (NAHosMIP)")
  ncatt_put(nc, 0, "delta_sv_definition", "M26 hos minus mean(con)")
  ncatt_put(nc, 0, "binning", sprintf("floor to %g Sv, bins with >= %d hosing years kept", 1, MIN_YEARS))
  ncatt_put(nc, 0, "years_used", sprintf("first %d hosing years (min of M26 valid and cube length)", B$n_used))
  ncatt_put(nc, 0, "pr_note", "hosing/piControl ratio, bin mean per calendar month; [0.1,10] clamp applied by the consumer")
  ncatt_put(nc, 0, "grid", "native model grid cropped to the E-OBS box + 3 deg, lon rewrapped to -180..180")
  ncatt_put(nc, 0, "created", format(Sys.time(), "%Y-%m-%d %H:%M"))
  nc_close(nc)
  cat("   wrote", basename(f), "\n\n")
  sapply(bins, function(b) length(b$years))
}

## ---- check against the existing files: the n_years vectors must match exactly ---------------
# This is the whole validation of the year -> bin mapping. The old files were built by a script
# that is not in the repo, so reproducing their n_years from M26 is the only way to prove this
# rebuild bins the same years - and therefore that the monthly fields describe the same bins.
EXPECTED <- list("IPSL-CM6A-LR"    = c(5, 8, 19, 11, 13, 11, 8, 9, 6, 8),
                 "EC-Earth3"       = c(4, 15, 20, 21, 8, 8, 9, 5, 6, 3),
                 "HadGEM3-GC3-1LL" = c(12, 19, 24, 10, 6, 7, 6, 8, 6),
                 "HadGEM3-GC3-1MM" = c(4, 8, 13, 21, 10, 7, 7, 8, 4, 4, 4))
got <- sapply(MODELS, build, simplify = FALSE)
ok  <- sapply(MODELS, function(m) identical(as.integer(got[[m]]), as.integer(EXPECTED[[m]])))
for (m in MODELS) cat(sprintf("%-16s n_years matches amoc_effect_bins_%s.nc: %s\n", m, m, ok[[m]]))
if (!all(ok)) stop("n_years does NOT reproduce the existing bins files - the year -> bin mapping is wrong")
cat("\nall four models reproduce the existing bin membership\n")
