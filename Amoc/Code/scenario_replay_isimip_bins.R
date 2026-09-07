# ISIMIP branch, bin-resolved: one perturbed-climate scenario per delta_Sv bin, per model, matched
# level-for-level to the NAHosMIP bins (isimip_bin_fields.R). Same engine, same interpolation, same
# self-check as the NAHosMIP bin branch (scenario_replay_bins.R) - only the field source differs, so
# this file duplicates its bin_scenarios()/run_branch() glue rather than parametrising that script,
# to keep the already-validated NAHosMIP path untouched (Appendix G).
# Author: Marco Bova
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))

MODELS <- c("ipsl-cm6a-lr", "ec-earth3")           # file-id list; display name comes from the .nc attribute
# emissions scenario, from ISIMIP_SCEN (default ssp126). ssp126 reads/writes bare names (unchanged);
# any other scenario gets a _<scen> suffix on the bin-field input and the scenario outputs.
SCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
SFX  <- if (SCEN == "ssp126") "" else paste0("_", SCEN)
bin_file <- function(m) file.path(d, sprintf("isimip_bin_fields_%s%s.nc", m, SFX))

## ---- per-bin fields -> E-OBS cells (identical to scenario_replay_bins.R::bin_scenarios) --------
bin_scenarios <- function(S, model) {
  nc  <- nc_open(bin_file(model))
  disp <- ncatt_get(nc, 0, "model")$value
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
  names(scen) <- sprintf("%s|%g", disp, lo)
  attr(scen, "meta") <- data.table(id = names(scen), model = disp, bin_id = (lo + hi) / 2,
                                   bin_lo = lo, bin_hi = hi, delta_sv = round(dsv, 4), n_years = ny)
  scen
}

## ---- run one branch over both models' bins ------------------------------------------------------
run_branch <- function(branch, out_file) {
  S <- eobs_setup(branch)
  cat(sprintf("\n########## ISIMIP-BINS | %s branch | %d E-OBS cells | %d groups | %d-%d ##########\n",
              toupper(branch), length(S$cells), length(S$grp), min(S$YRS), max(S$YRS)))
  cat("-- zero-delta self-check (gate: must be exact before any scenario is trusted) --\n")
  selfcheck(S)

  scen <- list(); meta <- list()
  for (m in MODELS) { s <- bin_scenarios(S, m); meta[[m]] <- attr(s, "meta"); scen <- c(scen, s) }
  meta <- rbindlist(meta)
  cat(sprintf("-- %d scenarios (%s) --\n", length(scen),
              paste(sprintf("%s:%d", unique(meta$model), sapply(unique(meta$model), function(m) sum(meta$model == m))), collapse = " ")))

  if (branch != "energy") for (m in unique(meta$model)) {
    k <- meta[model == m][which.min(bin_id)]$id
    cat(sprintf(" %s, bin %s Sv:\n", m, meta[id == k]$bin_id))
    order_check(S, scen[[k]])
  }
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
  invisible(out)
}

OUT <- c(energy      = sprintf("scenario_isimip_bins_energy_country%s.csv.gz", SFX),
         crop        = sprintf("scenario_isimip_bins_crop_window%s.csv.gz", SFX),
         crop_gddwin = sprintf("scenario_isimip_bins_crop_gddwin%s.csv.gz", SFX))
BRANCHES <- if (length(commandArgs(TRUE))) commandArgs(TRUE) else names(OUT)
for (b in BRANCHES) run_branch(b, OUT[[b]])
cat(sprintf("\nISIMIP bin branch done (%s).\n", SCEN))
