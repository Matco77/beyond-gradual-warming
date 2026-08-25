# AMOC branch: one perturbed-climate scenario per delta_Sv bin, per model.
# Author: Marco Bova
#
# Reads the MONTHLY per-bin fields (amoc_bin_fields.R), interpolates them onto the E-OBS grid, and
# replays the daily E-OBS record through scenario_engine.R. The delta is applied per CALENDAR MONTH
# to every day of that month, and BEFORE the thresholds - the hinges are evaluated on the shifted
# days, which is the whole point: under cooling, days cross the 5/28 C crop thresholds and the
# 15/24 C energy dead-bands differently, and a delta applied to a finished degree-day total cannot
# see that.
#
# Output (regenerable, so kept out of git - the crop table is ~4.2 M rows, past GitHub's limit):
#   scenario_bins_energy_country.csv.gz   schema of 13. + model, bin_id, bin_lo, bin_hi, delta_sv, n_years
#   scenario_bins_crop_window.csv.gz      schema of  8. + the same columns
# fread() reads .csv.gz transparently, so a panel builder still only changes a path.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))

MODELS <- c("IPSL-CM6A-LR", "EC-Earth3", "HadGEM3-GC3-1LL", "HadGEM3-GC3-1MM")
bin_file <- function(m) file.path(d, sprintf("amoc_bin_fields_%s_u03.nc", m))

## ---- per-bin fields -> E-OBS cells ------------------------------------------------------------
# Bilinear onto E-OBS cell centres (scenario_engine.R::to_eobs), one call per (variable, month,
# bin). Temperature deltas fall back to 0 outside the model domain, the precip factor to 1 - i.e.
# "no change", never a fabricated one. The [0.1, 10] clamp lands here, after interpolation, so the
# bin means and the per-year fields of the Phase-2 band get exactly the same treatment.
bin_scenarios <- function(S, model) {
  nc  <- nc_open(bin_file(model))
  lon <- as.numeric(ncvar_get(nc, "lon")); lat <- as.numeric(ncvar_get(nc, "lat"))
  V   <- list(dtg = "delta_tas", dtx = "delta_tasmax", dtn = "delta_tasmin", Rr = "pr_ratio")
  A   <- lapply(V, function(v) ncvar_get(nc, v))                    # [lon, lat, month, bin]
  lo  <- as.numeric(ncvar_get(nc, "bin_lo"));  hi  <- as.numeric(ncvar_get(nc, "bin_hi"))
  dsv <- as.numeric(ncvar_get(nc, "delta_sv_mean")); ny <- as.integer(ncvar_get(nc, "n_years"))
  nc_close(nc)

  scen <- lapply(seq_along(lo), function(k) {
    f <- lapply(names(V), function(nm) {
      M <- sapply(1:12, function(m) to_eobs(S, lon, lat, A[[nm]][, , m, k], outside = if (nm == "Rr") 1 else 0))
      if (nm == "Rr") pmin(pmax(M, PR_CLAMP[1]), PR_CLAMP[2]) else M
    })
    names(f) <- names(V); f
  })
  names(scen) <- sprintf("%s|%g", model, lo)
  attr(scen, "meta") <- data.table(id = names(scen), model = model, bin_id = (lo + hi) / 2,
                                   bin_lo = lo, bin_hi = hi, delta_sv = round(dsv, 4), n_years = ny)
  scen
}

## ---- run one branch over every model x bin ----------------------------------------------------
run_branch <- function(branch, out_file) {
  S <- eobs_setup(branch)
  cat(sprintf("\n########## %s branch | %d E-OBS cells | %d groups | %d-%d ##########\n",
              toupper(branch), length(S$cells), length(S$grp), min(S$YRS), max(S$YRS)))
  cat("-- zero-delta self-check (gate: must be exact before any scenario is trusted) --\n")
  selfcheck(S)

  scen <- list(); meta <- list()
  for (m in MODELS) { s <- bin_scenarios(S, m); meta[[m]] <- attr(s, "meta"); scen <- c(scen, s) }
  meta <- rbindlist(meta)
  cat(sprintf("-- %d scenarios (%s) --\n", length(scen),
              paste(sprintf("%s:%d", MODELS, sapply(MODELS, function(m) sum(meta$model == m))), collapse = " ")))

  t0  <- Sys.time()
  res <- run_scenarios(S, scen)
  cat(sprintf("-- %d scenarios in %.1f min --\n", length(scen),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  out <- merge(res, meta, by = "id", sort = FALSE)[, id := NULL]
  key <- switch(branch, energy = c("cntr", "year"), crop = c("NUTS_ID", "year", "window"),
                crop_gddwin = c("NUTS_ID", "year"))
  val <- switch(branch, energy = c("hdd_calendar", "hdd_octmar", "cdd_jja"),
                crop = c("gdd", "heat", "frost", "precip"),
                crop_gddwin = c("gdd", "heat", "frost", "precip", "n_day", "closed"))
  setcolorder(out, c("model", "bin_id", "bin_lo", "bin_hi", "delta_sv", "n_years", key, val))
  setorderv(out, c("model", "bin_id", key))
  fwrite(out, file.path(d, out_file))
  cat(sprintf("wrote %s | rows %d | models %d | bins %d\n",
              out_file, nrow(out), uniqueN(out$model), uniqueN(out[, .(model, bin_id)])))

  # what the perturbation actually did, per model x bin, against the unperturbed record
  h <- switch(branch,
    energy      = fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))[year %in% S$YRS],
    crop        = fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[year %in% S$YRS],
    crop_gddwin = fread(file.path(d, "11.crop_weather_gdd_window.csv"))[year %in% S$YRS][, closed := 1L])
  base <- sapply(val, function(v) mean(h[[v]], na.rm = TRUE))
  smry <- out[, c(.(n = n_years[1]), lapply(.SD, function(x) mean(x, na.rm = TRUE))),
              by = .(model, bin_id), .SDcols = val]
  cat("\n-- mean indicator by bin (historical baseline:",
      paste(sprintf("%s %.1f", val, base), collapse = ", "), ") --\n")
  for (m in MODELS) {
    cat(sprintf(" %s\n", m))
    z <- smry[model == m][order(bin_id)]
    for (i in seq_len(nrow(z))) cat(sprintf("   bin %6.1f Sv (n=%2d) %s\n", z$bin_id[i], z$n[i],
      paste(sprintf("%s %8.1f (%+6.1f%%)", val, unlist(z[i, val, with = FALSE]),
                    100 * (unlist(z[i, val, with = FALSE]) / base - 1)), collapse = "  ")))
  }
  invisible(out)
}

OUT <- c(energy      = "scenario_bins_energy_country.csv.gz",
         crop        = "scenario_bins_crop_window.csv.gz",
         crop_gddwin = "scenario_bins_crop_gddwin.csv.gz")
# branches selectable from the command line so a single one can be (re)run without redoing the rest
BRANCHES <- if (length(commandArgs(TRUE))) commandArgs(TRUE) else names(OUT)
for (b in BRANCHES) run_branch(b, OUT[[b]])
cat("\nAMOC bin branch done.\n")
