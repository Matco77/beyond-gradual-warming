# Delta-method scenario replay, delta_fields runs - ENERGY side.
# Perturb the daily E-OBS record by a climate delta, then recompute the country HDD/CDD with the
# SAME construction used in estimation. All of that now lives in scenario_engine.R, shared with the
# crop runner and with the AMOC bin runner, so there is one machine and one self-check rather than
# a copy per scenario source.
#
# These are the "last third" deltas: the end-state of each hosing run, one level shift per cell.
# They are NOT the per-bin fields - for the impact-vs-Sv curve see scenario_replay_bins.R, which
# resolves the delta by calendar month and by delta_Sv bin. This runner is kept as the reference
# point those bins are read against.
#
# Runs are AUTO-DISCOVERED from delta_fields/ (every model x protocol that has a tas delta), so
# IPSL-CM6A-LR (u03 only) and any new model are picked up with no code edit. Replayed over the
# energy estimation window (YRS_ENERGY) so scenario and historical share the estimation sample.
# Remaining TODO: these deltas are ANNUAL and applied to every month (the seasonal DJF/JJA fields
# exist but carry no CRS, so terra cannot grid them). The bin fields do not have that limitation.
# Author: Marco Bova
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))

delta_dir <- file.path(d, "delta_fields")
# auto-discover (model, exp) from the tas delta files: delta_<MODEL>_<EXP>_tas_lastthird.nc
RUNS <- lapply(list.files(delta_dir, "^delta_.*_tas_lastthird\\.nc$"), function(f) {
  s <- sub("^delta_(.*)_tas_lastthird\\.nc$", "\\1", f)        # <MODEL>_<EXP>
  list(model = sub("_[^_]*$", "", s), exp = sub(".*_", "", s)) # model may carry hyphens, exp is the last _token
})
cat("runs discovered:", paste(sapply(RUNS, function(r) paste(r$model, r$exp)), collapse = " | "), "\n\n")

# delta_fields sit on a regular grid that terra can read, so they keep the terra bilinear resample
# they were built with. (The bin fields go through fields::interp.surface instead, because terra
# rejects the native model grids: "lat not regularly spaced". Both are bilinear on cell centres.)
tmpl <- rast(ncf("tg"))[[1]]
month_delta <- function(S, model, exp, var) {
  ann <- values(resample(rast(file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, exp, var))),
                         tmpl, method = "bilinear"))[, 1]
  matrix(ann[S$cells], length(S$cells), 12)                    # annual delta, every month
}

S <- eobs_setup("energy")
cat("-- zero-delta self-check --\n"); selfcheck(S)

scen <- lapply(RUNS, function(r) list(dtg = month_delta(S, r$model, r$exp, "tas"),
                                      dtx = month_delta(S, r$model, r$exp, "tasmax"),
                                      dtn = month_delta(S, r$model, r$exp, "tasmin"),
                                      Rr  = matrix(1, length(S$cells), 12)))   # unused on this branch
names(scen) <- sapply(RUNS, function(r) paste(r$model, r$exp, sep = "|"))
res <- run_scenarios(S, scen)
res[, c("model", "exp") := tstrsplit(id, "|", fixed = TRUE)][, id := NULL]

VAL <- c("hdd_calendar", "hdd_octmar", "cdd_jja")
setcolorder(res, c("model", "exp", "cntr", "year", VAL))
fwrite(res, file.path(d, "scenario_energy_country.csv"))
cat("\nwrote scenario_energy_country.csv | rows", nrow(res), "| runs", uniqueN(res[, .(model, exp)]), "\n")

## ---- what the perturbation did, per run ------------------------------------------------------
hist <- fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))[year %in% S$YRS]
M <- function(x) mean(x, na.rm = TRUE)          # MT (Malta) is off-domain in E-OBS -> NA hdd, skip it
for (r in RUNS) {
  z <- res[model == r$model & exp == r$exp]
  cat("\n=== AMOC", r$model, r$exp, "vs historical (", min(S$YRS), "-", max(S$YRS),
      "), pop-weighted country-year means ===\n")
  for (v in VAL) cat(sprintf("  %-13s hist %6.0f -> scen %6.0f (%+.1f%%)\n",
                             v, M(hist[[v]]), M(z[[v]]), 100 * (M(z[[v]]) / M(hist[[v]]) - 1)))
}
