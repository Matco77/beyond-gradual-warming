# Shared engine for every delta-method scenario replay (delta_fields runs and AMOC bin runs alike).
# Author: Marco Bova
#
# One copy of: the E-OBS plumbing, the perturbation, the indicator construction, the aggregation to
# the estimation geography, and the zero-delta self-check. The runners on top of it differ only in
# where the per-cell monthly delta comes from and what they write out. Same reason
# weather_indicators.R exists: if the regressor is not built identically in estimation and in
# replay, beta is not transportable, and the surest way to break that is a second copy of the code.
#
# The perturbation is applied to the DAILY field, BEFORE any threshold - the delta shifts the whole
# distribution and the hinges are then evaluated on the shifted days. Applying a delta to an
# already-computed degree-day total instead would assume the threshold never binds differently,
# which is exactly what a cooling scenario violates.
#
# Interpolation onto the E-OBS grid: BILINEAR on cell centres, via fields::interp.surface, for both
# temperature and precipitation. That is the mechanism scenario_replay_crop.R already used for the
# precip ratio; it is used here for temperature too because terra rejects the native model grids
# ("lat not regularly spaced"), so one path serves all fields instead of two.
suppressMessages({library(terra); library(ncdf4); library(Matrix); library(data.table); library(fields)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

PR_CLAMP <- c(0.1, 10)      # multiplicative precip factor bounds (ymondiv spikes in dry cells)

## ---- E-OBS axes, read once ------------------------------------------------------------------
.ax <- local({
  nc <- nc_open(ncf("tg")); lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude")
  tv <- as.Date(ncvar_get(nc, "time"), origin = "1950-01-01"); nc_close(nc)
  list(lon = lon, lat = lat, nlon = length(lon), nlat = length(lat), tvec = tv,
       yr = as.integer(format(tv, "%Y")), mo = as.integer(format(tv, "%m")))
})

## ---- branch setup ---------------------------------------------------------------------------
# "energy": population-weighted country aggregate, HDD/CDD via the within-day integration.
# "crop"  : area-weighted NUTS3 aggregate, daily indicators, monthly then summed over the window.
eobs_setup <- function(branch) {
  tmpl <- rast(ncf("tg"))[[1]]
  if (branch == "energy") {
    w <- fread(file.path(d, "12.eobs_population_weights.csv"))
    cells <- sort(unique(w$eobs_cell_id)); grp <- sort(unique(w$cntr))
    W <- sparseMatrix(i = match(w$cntr, grp), j = match(w$eobs_cell_id, cells),
                      x = w$w_pop, dims = c(length(grp), length(cells)))
    YRS <- YRS_ENERGY
    # the Oct-Mar season labelled YRS[1] needs Oct-Dec of the year before it
    read_yrs <- c(min(YRS) - 1L, YRS); vars <- c("tg", "tx", "tn")
  } else {
    w <- fread(file.path(d, "6.eobs_to_nuts_crosswalk.csv"))
    cells <- sort(unique(w$eobs_cell_id)); grp <- sort(unique(w$NUTS_ID))
    W <- sparseMatrix(i = match(w$NUTS_ID, grp), j = match(w$eobs_cell_id, cells),
                      x = w$nuts_weight, dims = c(length(grp), length(cells)))
    YRS <- YRS_CROP; read_yrs <- YRS; vars <- c("tg", "tx", "tn", "rr")
  }
  xy <- xyFromCell(tmpl, cells)
  loni <- match(round(xy[, 1], 3), round(.ax$lon, 3)); lati <- match(round(xy[, 2], 3), round(.ax$lat, 3))
  stopifnot(!anyNA(loni), !anyNA(lati))
  list(branch = branch, cells = cells, xy = xy, grp = grp, W = W, YRS = YRS,
       read_yrs = read_yrs, vars = vars, ncdf_row = (lati - 1) * .ax$nlon + loni)
}

## ---- model field -> E-OBS cells (bilinear on cell centres) ----------------------------------
# Z is [lon, lat] on the model grid; lon/lat must be increasing (the bin files are written that
# way). Cells outside the model domain get `outside`: 0 for an additive delta, 1 for a ratio.
to_eobs <- function(S, lon, lat, Z, outside = 0) {
  v <- interp.surface(list(x = lon, y = lat, z = Z), S$xy)
  v[!is.finite(v)] <- outside
  v
}

## ---- weighted aggregation, NA-safe (renormalised over the cells with data that day) ---------
.wmean <- function(S, v) { p <- !is.na(v); v[!p] <- 0
  r <- as.matrix(S$W %*% v) / as.matrix(S$W %*% (p * 1)); r[!is.finite(r)] <- NA; r }
.ind_month <- function(S, mat, mm, fun) {
  v <- fun(mat); pres <- !is.na(v); v[!pres] <- 0
  reg <- as.matrix(S$W %*% v) / as.matrix(S$W %*% (pres * 1)); reg[!is.finite(reg)] <- NA
  t(rowsum(t(reg), mm))                                     # group x month
}

## ---- run a set of scenarios ------------------------------------------------------------------
# scen: named list; each element is list(dtg, dtx, dtn, Rr), matrices [length(S$cells) x 12] of
# per-calendar-month deltas already on the E-OBS cells (Rr only used by the crop branch).
# Scenarios run in chunks so E-OBS is read once per chunk instead of once per scenario, while peak
# memory stays bounded; the read is ~3 s/year against ~5 s/year/scenario of indicator work, so the
# chunk size trades a few minutes of re-reading against RAM and is not worth tuning further.
run_scenarios <- function(S, scen, chunk = 8L, verbose = TRUE) {
  ids <- names(scen); res <- list()
  for (g in split(seq_along(ids), ceiling(seq_along(ids) / chunk))) {
    ncs <- lapply(S$vars, function(v) nc_open(ncf(v))); names(ncs) <- S$vars
    rd <- function(v, ti) { a <- ncvar_get(ncs[[v]], v, start = c(1, 1, ti[1]), count = c(-1, -1, length(ti)))
      dim(a) <- c(.ax$nlon * .ax$nlat, length(ti)); a[S$ncdf_row, , drop = FALSE] }
    acc <- vector("list", length(g)); names(acc) <- ids[g]
    for (k in seq_along(acc)) acc[[k]] <- list()
    for (y in S$read_yrs) {
      ti <- which(.ax$yr == y); mm <- .ax$mo[ti]
      RAW <- lapply(S$vars, function(v) rd(v, ti)); names(RAW) <- S$vars
      for (k in seq_along(g)) {
        sc <- scen[[ids[g[k]]]]
        TG <- RAW$tg + sc$dtg[, mm]; TX <- RAW$tx + sc$dtx[, mm]; TN <- RAW$tn + sc$dtn[, mm]
        if (S$branch == "energy") {
          acc[[k]][[length(acc[[k]]) + 1L]] <- data.table(
            cntr = rep(S$grp, length(ti)), date = rep(.ax$tvec[ti], each = length(S$grp)),
            hdd = as.vector(.wmean(S, within_day(TG, TX, TN, hdd_hinge))),
            cdd = as.vector(.wmean(S, within_day(TG, TX, TN, cdd_hinge))))
        } else {
          RR <- RAW$rr * sc$Rr[, mm]                        # temp additive, precip multiplicative
          # NUTS3 x month first, rounded to 3dp, then summed over the window: the exact two-step of
          # cropweather_eobs_nuts3.R (:140 rounds the monthly file, :152 sums the window).
          Mo <- list(gdd    = round(.ind_month(S, TG, mm, gdd_daily), 3),
                     heat   = round(.ind_month(S, TX, mm, heat_daily), 3),
                     frost  = round(.ind_month(S, TN, mm, frost_daily), 3),
                     precip = round(.ind_month(S, RR, mm, identity), 3))
          acc[[k]][[length(acc[[k]]) + 1L]] <- rbindlist(lapply(names(CROP_WINDOWS), function(lab) {
            keep <- as.integer(colnames(Mo$gdd)) %in% CROP_WINDOWS[[lab]]
            ws <- function(A) round(rowSums(A[, keep, drop = FALSE]), 3)   # NA month -> NA season
            data.table(NUTS_ID = S$grp, year = y, window = lab, gdd = ws(Mo$gdd), heat = ws(Mo$heat),
                       frost = ws(Mo$frost), precip = ws(Mo$precip)) }))
        }
      }
      if (verbose) { cat("."); flush.console() }
    }
    for (nc in ncs) nc_close(nc)
    for (k in seq_along(g)) res[[ids[g[k]]]] <- .reduce(S, rbindlist(acc[[k]]))[, id := ids[g[k]]]
    if (verbose) cat(sprintf(" [%s]\n", paste(ids[g], collapse = ", ")))
  }
  rbindlist(res, use.names = TRUE)
}

CROP_WINDOWS <- list("Mar-Jul" = 3:7, "Apr-Aug" = 4:8)   # main window + the robustness window

# Reduction to the estimation panel geometry. The energy side reproduces
# energyweather_eobs_country.R:114-117 line for line, the >=150-day rule included; the crop side is
# already reduced per year above and only needs the read window trimmed.
.reduce <- function(S, dt) {
  if (S$branch == "crop") return(dt[year %in% S$YRS][order(NUTS_ID, year, window)])
  dt[, `:=`(y = year(date), m = month(date))]
  cal <- dt[y %in% S$YRS, .(hdd_calendar = sum(hdd)), by = .(cntr, year = y)]
  oct <- dt[m %in% c(10, 11, 12, 1, 2, 3)][, .(hdd_octmar = sum(hdd), nd = .N),
            by = .(cntr, year = y + (m >= 10L))][nd >= 150 & year %in% S$YRS, .(cntr, year, hdd_octmar)]
  jja <- dt[m %in% 6:8 & y %in% S$YRS, .(cdd_jja = sum(cdd)), by = .(cntr, year = y)]
  out <- Reduce(function(a, b) merge(a, b, by = c("cntr", "year"), all = TRUE),
                list(cal, oct, jja))[order(cntr, year)]
  out[, c("hdd_calendar", "hdd_octmar", "cdd_jja") := lapply(.SD, round, 2),
      .SDcols = c("hdd_calendar", "hdd_octmar", "cdd_jja")][]   # same rounding as 13.
}

## ---- zero-delta self-check -------------------------------------------------------------------
# With delta = 0 and ratio = 1 the engine is the estimation pipeline with a "+0" and a "*1" in it,
# so any difference from 8. / 13. is a defect of the engine (wrong cells, wrong month mapping,
# wrong window rule), not a climate signal. Every number downstream is a DIFFERENCE of these
# indicators, so a failure here biases all of them silently -> hard stop, no fallback.
zero_scenario <- function(S) list(dtg = matrix(0, length(S$cells), 12), dtx = matrix(0, length(S$cells), 12),
                                  dtn = matrix(0, length(S$cells), 12), Rr = matrix(1, length(S$cells), 12))
selfcheck <- function(S, verbose = TRUE) {
  z <- run_scenarios(S, list("ZERO-DELTA" = zero_scenario(S)), verbose = verbose)
  if (S$branch == "energy") {
    key <- c("cntr", "year"); vars <- c("hdd_calendar", "hdd_octmar", "cdd_jja")
    h <- fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))[year %in% S$YRS]
  } else {
    key <- c("NUTS_ID", "year", "window"); vars <- c("gdd", "heat", "frost", "precip")
    h <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[year %in% S$YRS]
  }
  cmp <- merge(z, h[, c(key, vars), with = FALSE], by = key, suffixes = c("_scen", "_hist"), all = TRUE)
  ok <- TRUE
  for (v in vars) {
    a <- cmp[[paste0(v, "_scen")]]; b <- cmp[[paste0(v, "_hist")]]
    nmis <- sum(xor(is.na(a), is.na(b))); dmax <- max(abs(a - b), na.rm = TRUE)
    cat(sprintf("  %-13s max|scen-hist| %.3g | NA mismatch %d\n", v, dmax, nmis))
    ok <- ok && nmis == 0L && dmax == 0
  }
  cat(sprintf("  rows %d | groups %d | years %d-%d\n", nrow(cmp), uniqueN(cmp[[key[1]]]),
              min(S$YRS), max(S$YRS)))
  if (!ok) stop(sprintf("ZERO-DELTA SELF-CHECK FAILED (%s): the engine does not reproduce the historical indicators", S$branch))
  cat("  zero-delta self-check PASSED\n")
  invisible(TRUE)
}
