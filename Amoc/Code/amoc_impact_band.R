# Phase 2, uncertainty band: run the WHOLE chain on each hosing year of a bin separately, and use
# the spread of the FINAL results as the band.
# Author: Marco Bova
#
# Why not propagate the sd of the delta field analytically: the indicators sit behind thresholds
# (5/28 C for crop, the 15/24 C dead-bands for energy), so the mean of the effects is not the effect
# of the mean field. A cell-wise sd of the delta cannot be pushed through a hinge. Running the chain
# per year and dispersing at the end is the only version of "uncertainty" that keeps the
# non-linearity, which is exactly where the AMOC signal lives.
#
# Scope: IPSL-CM6A-LR and EC-Earth3 (98 + 99 = 197 hosing years). HadGEM LL/MM keep the central
# estimate only.
#
# Nothing about the chain is re-implemented here: the fields come from the same builder
# (amoc_bin_fields.R, mode = "year"), the replay from the same engine, and the beta application from
# the same crop_effect()/energy_effect() as the central estimate - only the grouping key changes.
# So the band cannot silently disagree with the number it is a band around.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
cat("beta and weights loaded from amoc_impact.R (central estimate recomputed silently)\n")

BAND_MODELS <- c("IPSL-CM6A-LR", "EC-Earth3")
KEY <- c(BINKEY, "hos_year")
BND <- grab(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_replay_bins.R"),
            c("bin_file", "bin_scenarios"))
BND$bin_file <- function(m) file.path(d, sprintf("amoc_year_fields_%s_u03.nc", m))   # per-year fields

year_scen <- function(S, m) {
  sc <- BND$bin_scenarios(S, m)
  nc <- nc_open(BND$bin_file(m)); hy <- as.integer(ncvar_get(nc, "hos_year")); nc_close(nc)
  # bin_scenarios names a scenario by its bin edge. Once every year is its own scenario that name is
  # NOT unique, and run_scenarios looks scenarios up by name - so rename before anything reads them.
  names(sc) <- sprintf("%s|y%03d", m, hy)
  mt <- attr(sc, "meta"); mt[, `:=`(id = names(sc), hos_year = hy)]
  attr(sc, "meta") <- mt
  sc
}

## ---- one model, one branch: replay every year, then apply beta straight away ------------------
# CHECKPOINTED per (branch, model). The energy branch alone is ~5 h and this job has already been
# killed once mid-run; a cache file per branch x model means an interruption costs at most one
# model, not the whole chain. Delete the cache file to force a recompute.
# Reduced to EU-aggregated effects before caching: the crop branch is ~10M scenario rows per model,
# and the effects are four orders of magnitude smaller.
cache <- function(branch, m) file.path(d, sprintf("amoc_band_%s_%s.csv", branch, gsub("[^A-Za-z0-9]", "", m)))

band_one <- function(branch, m) {
  f <- cache(branch, m)
  if (file.exists(f)) { cat(sprintf("-- %-6s | %-14s | cached, skipped\n", toupper(branch), m)); return(fread(f)) }
  S  <- eobs_setup(branch)
  sc <- year_scen(S, m); meta <- attr(sc, "meta")
  cat(sprintf("\n-- %s | %s | %d hosing years --\n", toupper(branch), m, length(sc)))
  t0  <- Sys.time()
  res <- merge(run_scenarios(S, sc), meta, by = "id", sort = FALSE)[, id := NULL]
  cat(sprintf("   replay %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  eu <- if (branch == "crop") {
    merge(crop_effect(res, KEY), w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "component")][, branch := "crop"][]
  } else {
    merge(energy_effect(res, KEY), w_pop, by = c("cntr", "country_id"))[
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "fuel", "component")][
      , branch := paste("energy", fuel)][]
  }
  fwrite(eu, f); eu
}

# crop first: ~7x cheaper per scenario, so a failure surfaces in minutes rather than hours
BAND <- rbindlist(lapply(c("crop", "energy"),
                         function(b) rbindlist(lapply(BAND_MODELS, function(m) band_one(b, m)), fill = TRUE)),
                  fill = TRUE)

## ---- disperse across the years of each bin -----------------------------------------------------
# Also reported: the central estimate (chain run on the bin-MEAN field) next to the mean of the
# per-year effects. The gap between them IS the non-linearity the band exists to capture - if the
# two agreed, the band would not have been worth the compute.
central <- rbind(
  fread(file.path(d, "amoc_impact_crop_eu.csv"))[, .(model, bin_id, component, central = dln, branch = "crop")],
  fread(file.path(d, "amoc_impact_energy_eu.csv"))[, .(model, bin_id, component, central = dln,
                                                       branch = paste("energy", fuel))])

band <- BAND[, .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)),
             by = .(branch, model, bin_id, component)]
band <- merge(band, central, by = c("branch", "model", "bin_id", "component"), all.x = TRUE)
fwrite(band[order(branch, model, bin_id, component)], file.path(d, "amoc_band_summary.csv"))

tot  <- BAND[, .(dln = sum(dln)), by = .(branch, model, bin_id, hos_year)][   # total = sum of components
  , .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)), by = .(branch, model, bin_id)]
ctot <- central[, .(central = sum(central)), by = .(branch, model, bin_id)]
tot  <- merge(tot, ctot, by = c("branch", "model", "bin_id"), all.x = TRUE)
fwrite(tot[order(branch, model, bin_id)], file.path(d, "amoc_band_total.csv"))

pc <- function(x) 100 * (exp(x) - 1)
for (br in sort(unique(tot$branch))) for (m in BAND_MODELS) {
  z <- tot[branch == br & model == m][order(bin_id)]
  if (!nrow(z)) next
  cat(sprintf("\n=== %s | %s | total effect, spread across the bin's hosing years ===\n", br, m))
  cat(sprintf("  %7s %4s %10s %10s %9s %10s %10s\n", "bin Sv", "n", "central%", "mean%", "sd(log)%", "min%", "max%"))
  for (i in seq_len(nrow(z)))
    cat(sprintf("  %7.1f %4d %9.2f%% %9.2f%% %8.2f%% %9.2f%% %9.2f%%\n", z$bin_id[i], z$n[i],
                pc(z$central[i]), pc(z$mean[i]), 100 * z$sd[i], pc(z$lo[i]), pc(z$hi[i])))
}
cat("\nwrote amoc_band_{summary,total}.csv\n")
