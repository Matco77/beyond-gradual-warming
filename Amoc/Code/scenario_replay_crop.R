# Delta-method scenario replay, delta_fields runs - CROP side.
# Perturb the daily E-OBS record by an AMOC delta, then recompute the NUTS3 crop-window indicators
# with the SAME construction used in estimation. That machinery lives in scenario_engine.R, shared
# with the energy runner and the AMOC bin runner.
#
# Temperature: ADDITIVE delta (delta_fields/delta_<MODEL>_<EXP>_{tas,tasmax,tasmin}_lastthird.nc).
# Precip:      MULTIPLICATIVE ratio R = hosing/piControl (anomaly_output .../pr_*ratio.nc, last-third
#              mean per cell, clamped to [0.1,10]) -> rr_scen = rr * R. A ratio keeps precip >= 0 and
#              scales with local climate; an additive precip delta would go negative in dry cells.
#
# These are the "last third" deltas, the end-state of each hosing run. For the impact-vs-Sv curve
# see scenario_replay_bins.R, which resolves the delta by calendar month and by delta_Sv bin.
# Both estimation windows are produced (Mar-Jul main, Apr-Aug robustness = the _mj / _aa regressors
# of regression.R), over YRS_CROP so scenario and historical share the estimation sample.
# Remaining TODO: these deltas and ratios are ANNUAL, applied to every month.
# Author: Marco Bova
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))

delta_dir <- file.path(d, "delta_fields")
RUNS <- lapply(list.files(delta_dir, "^delta_.*_tas_lastthird\\.nc$"), function(f) {
  s <- sub("^delta_(.*)_tas_lastthird\\.nc$", "\\1", f)
  list(model = sub("_[^_]*$", "", s), exp = sub(".*_", "", s))
})
cat("runs:", paste(sapply(RUNS, function(r) paste(r$model, r$exp)), collapse = " | "), "\n\n")

tmpl <- rast(ncf("tg"))[[1]]
month_delta <- function(S, model, exp, var) {               # additive, terra bilinear (regular grid)
  ann <- values(resample(rast(file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, exp, var))),
                         tmpl, method = "bilinear"))[, 1]
  matrix(ann[S$cells], length(S$cells), 12)
}
ratio_path <- function(model, exp) {
  base <- file.path(d, "anomaly_output")
  if (model == "EC-Earth3")         file.path(base, "ECHearth3_anomaly", sprintf("pr_Amon_EC-Earth3_hos-%s-hos_ratio.nc", exp))
  else if (model == "IPSL-CM6A-LR") file.path(base, "IPSL_anomaly",     sprintf("pr_Amon_IPSL-CM6A-LR_hos-%s-hos_ratio.nc", exp))
  else file.path(base, if (grepl("LL", model)) "LL_anomaly" else "MM_anomaly",
                 sprintf("pr_ratio_%s-hos_over_piControl_1850-1949.nc", exp))
}
month_ratio <- function(S, model, exp) {                    # multiplicative precip factor at E-OBS cells
  # Read via ncdf4 (terra rejects the native model grid: "lat not regularly spaced"). The ratio is
  # DIMENSIONLESS -> do NOT apply read_europe_cube's pr x86400. The clamp is applied on the model
  # grid, BEFORE interpolation, which is where these runs were built (the bin runner clamps after
  # interpolation instead; the two differ only next to a clamped cell, and no clamped cell is
  # within one model cell of any E-OBS cell these branches use).
  nc <- nc_open(ratio_path(model, exp))
  lon <- as.numeric(ncvar_get(nc, "lon")); lat <- as.numeric(ncvar_get(nc, "lat"))
  R3 <- ncvar_get(nc, "pr"); nc_close(nc)                   # [lon, lat, time], hosing/piControl ratio
  nt <- dim(R3)[3]; Rm <- apply(R3[, , (nt - nt %/% 3 + 1):nt, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  lon <- ((lon + 180) %% 360) - 180; ox <- order(lon); lon <- lon[ox]; Rm <- Rm[ox, ]   # interp needs increasing x,y
  if (lat[1] > lat[length(lat)]) { oy <- order(lat); lat <- lat[oy]; Rm <- Rm[, oy] }
  Rm[!is.finite(Rm)] <- 1; Rm <- pmin(pmax(Rm, PR_CLAMP[1]), PR_CLAMP[2])
  v <- interp.surface(list(x = lon, y = lat, z = Rm), S$xy) # bilinear to E-OBS cell centres
  v[!is.finite(v)] <- 1                                     # cells outside the model domain: no precip change
  matrix(v, length(S$cells), 12)                            # annual ratio for every month
}

S <- eobs_setup("crop")
cat("-- zero-delta self-check --\n"); selfcheck(S)

scen <- lapply(RUNS, function(r) list(dtg = month_delta(S, r$model, r$exp, "tas"),
                                      dtx = month_delta(S, r$model, r$exp, "tasmax"),
                                      dtn = month_delta(S, r$model, r$exp, "tasmin"),
                                      Rr  = month_ratio(S, r$model, r$exp)))
names(scen) <- sapply(RUNS, function(r) paste(r$model, r$exp, sep = "|"))
res <- run_scenarios(S, scen)
res[, c("model", "exp") := tstrsplit(id, "|", fixed = TRUE)][, id := NULL]

VAL <- c("gdd", "heat", "frost", "precip")
setcolorder(res, c("model", "exp", "NUTS_ID", "year", "window", VAL))
fwrite(res, file.path(d, "scenario_crop_window.csv"))
cat("wrote scenario_crop_window.csv | rows", nrow(res), "| runs", uniqueN(res[, .(model, exp)]),
    "| windows", paste(unique(res$window), collapse = ","), "\n")

## ---- what the perturbation did, per run x window ---------------------------------------------
hist <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[year %in% S$YRS]
for (r in RUNS) for (lab in names(CROP_WINDOWS)) {
  z <- res[model == r$model & exp == r$exp & window == lab]; h <- hist[window == lab]
  cat(sprintf("=== CROP %s %s vs historical (%d-%d), NUTS3-year %s means ===\n",
              r$model, r$exp, min(S$YRS), max(S$YRS), lab))
  for (v in VAL) { a <- mean(z[[v]], na.rm = TRUE); b <- mean(h[[v]], na.rm = TRUE)
    cat(sprintf("  %-6s %8.2f -> %8.2f (%+.1f%%)\n", v, b, a, 100 * (a / b - 1))) }
  cat("\n")
}
