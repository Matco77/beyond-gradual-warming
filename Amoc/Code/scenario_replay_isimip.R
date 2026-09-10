# ISIMIP branch: one perturbed-climate scenario per model, replayed through the SAME engine as the
# AMOC bin branch (scenario_replay_bins.R). Only the delta source differs.
# Author: Marco Bova
#
# ISIMIP has daily data, so nothing here is forced to use the delta method for lack of alternative -
# it is a deliberate choice to keep the two branches procedurally identical: the AMOC branch must
# perturb E-OBS because NAHosMIP is monthly-only, and running ISIMIP directly on its own daily
# levels while AMOC goes through a delta would mean the two effects differ partly by METHOD, not
# only by forcing. Instead both branches perturb the SAME E-OBS baseline with a monthly delta, so
# a comparison between them isolates the climate signal and climatological biases cancel on each
# side. See Appendix G (methodology) for the argument in full and its cost: a monthly delta cannot
# express a change in the SHAPE of the daily distribution, in neither branch.
#
# There are no delta_Sv bins on this branch - ISIMIP is one emissions trajectory, not an ensemble
# indexed by AMOC state - so this is one scenario per model, not forty.
#
# EMISSIONS SCENARIO. From ISIMIP_SCEN (default ssp126). ssp126 reads/writes the bare file names
# (unchanged); any other scenario gets a _<scen> suffix on both the delta-field input and the
# scenario outputs, so ssp370 runs alongside ssp126 without disturbing it.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))

# ISIMIP's own model list, kept separate from scenario_replay_bins.R's AMOC-branch MODELS: the two
# never need to be the same set, and coupling them would make an AMOC-only model addition silently
# require an ISIMIP field that does not exist.
MODELS <- c("ipsl-cm6a-lr", "ec-earth3")
SCEN   <- Sys.getenv("ISIMIP_SCEN", "ssp126")
SFX    <- if (SCEN == "ssp126") "" else paste0("_", SCEN)
field_file <- function(m) file.path(d, sprintf("isimip_delta_fields_%s%s.nc", m, SFX))

## ---- per-model field -> E-OBS cells (bilinear on cell centres, same as bin_scenarios()) --------
isimip_scenario <- function(S, model) {
  nc  <- nc_open(field_file(model))
  lon <- as.numeric(ncvar_get(nc, "lon")); lat <- as.numeric(ncvar_get(nc, "lat"))
  V   <- list(dtg = "delta_tas", dtx = "delta_tasmax", dtn = "delta_tasmin", Rr = "pr_ratio")
  A   <- lapply(V, function(v) ncvar_get(nc, v))          # [lon, lat, month]
  nc_close(nc)
  f <- lapply(names(V), function(nm) {
    M <- sapply(1:12, function(m) to_eobs(S, lon, lat, A[[nm]][, , m], outside = if (nm == "Rr") 1 else 0))
    if (nm == "Rr") pmin(pmax(M, PR_CLAMP[1]), PR_CLAMP[2]) else M
  })
  names(f) <- names(V); f
}

## ---- run one branch over the two models ---------------------------------------------------------
run_branch <- function(branch, out_file) {
  S <- eobs_setup(branch)
  cat(sprintf("\n########## ISIMIP | %s branch | %d E-OBS cells | %d groups | %d-%d ##########\n",
              toupper(branch), length(S$cells), length(S$grp), min(S$YRS), max(S$YRS)))
  cat("-- zero-delta self-check (gate: must be exact before any scenario is trusted) --\n")
  selfcheck(S)

  scen <- setNames(lapply(MODELS, function(m) isimip_scenario(S, m)), MODELS)
  if (branch != "energy") for (m in MODELS) { cat(sprintf(" %s:\n", m)); order_check(S, scen[[m]]) }

  t0  <- Sys.time()
  res <- run_scenarios(S, scen)
  cat(sprintf("-- %d scenarios in %.1f min --\n", length(scen),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  out <- res[, model := id][, id := NULL]
  key <- switch(branch, energy = c("cntr", "year"), crop = c("NUTS_ID", "year", "window"),
                crop_gddwin = c("NUTS_ID", "year"))
  val <- switch(branch, energy = c("hdd_calendar", "hdd_octmar", "cdd_jja"),
                crop = c("gdd", "heat", "frost", "precip"),
                crop_gddwin = c("gdd", "heat", "frost", "precip", "n_day", "closed"))
  setcolorder(out, c("model", key, val)); setorderv(out, c("model", key))
  fwrite(out, file.path(d, out_file))
  cat(sprintf("wrote %s | rows %d | models %d\n", out_file, nrow(out), uniqueN(out$model)))

  h <- switch(branch,
    energy      = fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))[year %in% S$YRS],
    crop        = fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[year %in% S$YRS],
    crop_gddwin = fread(file.path(d, "11.crop_weather_gdd_window.csv"))[year %in% S$YRS][, closed := 1L])
  base <- sapply(val, function(v) mean(h[[v]], na.rm = TRUE))
  cat("\n-- mean indicator by model (historical baseline:",
      paste(sprintf("%s %.1f", val, base), collapse = ", "), ") --\n")
  smry <- out[, lapply(.SD, mean, na.rm = TRUE), by = model, .SDcols = val]
  for (i in seq_len(nrow(smry))) cat(sprintf(" %-14s %s\n", smry$model[i],
    paste(sprintf("%s %8.1f (%+6.1f%%)", val, unlist(smry[i, val, with = FALSE]),
                  100 * (unlist(smry[i, val, with = FALSE]) / base - 1)), collapse = "  ")))
  invisible(out)
}

OUT <- c(energy      = sprintf("scenario_isimip_energy_country%s.csv", SFX),
         crop        = sprintf("scenario_isimip_crop_window%s.csv", SFX),
         crop_gddwin = sprintf("scenario_isimip_crop_gddwin%s.csv", SFX))
BRANCHES <- if (length(commandArgs(TRUE))) commandArgs(TRUE) else names(OUT)
for (b in BRANCHES) run_branch(b, OUT[[b]])
cat(sprintf("\nISIMIP branch done (%s).\n", SCEN))
