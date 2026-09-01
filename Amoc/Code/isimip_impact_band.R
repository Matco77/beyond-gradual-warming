# ISIMIP branch, uncertainty band: run the WHOLE chain on each individual ssp126 year of a bin
# separately, and use the spread of the FINAL results as the band - same construction as
# amoc_impact_band.R, applied to isimip_bin_fields.R's bins instead of NAHosMIP's.
# Author: Marco Bova
#
# Why this exists. Each ISIMIP bin already averages over several ssp126 years (5 to 17) - the exact
# situation amoc_impact_band.R was built for on the NAHosMIP side. Appendix J.6 flagged the absence
# of this band as an open gap, not as established to be unnecessary; this fills it.
#
# THREE branches here, not two: energy and crop (spec A) mirror amoc_impact_band.R exactly; crop
# spec C (thermal-time window, GW_X/gw_fit from amoc_impact_gddwin.R) is a NEW branch with no
# NAHosMIP-side counterpart - the NAHosMIP band was never built for spec C. Kept separate and
# labelled "crop_gddwin" so the asymmetry is visible in the output, not hidden.
#
# Nothing about the chain is reimplemented: the fields come from isimip_bin_fields.R's mode = "year",
# the replay from the same scenario_engine.R, and the beta application from the same
# crop_effect()/energy_effect()/gw_fit as isimip_impact_bins.R's central estimate - only the
# grouping key changes (bin_id -> bin_id + ssp_year), so the band cannot silently disagree with the
# number it is a band around.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_gddwin.R"))))
cat("beta, weights and both crop specifications loaded (AMOC central estimate recomputed silently)\n")

ISIMIP_MODELS <- c("ipsl-cm6a-lr", "ec-earth3")     # file ids, as scenario_replay_isimip_bins.R uses
KEY <- c(BINKEY, "ssp_year")
BND <- grab(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_replay_isimip_bins.R"),
            c("bin_file", "bin_scenarios"))
BND$bin_file <- function(m) file.path(d, sprintf("isimip_year_fields_%s.nc", m))   # per-year fields

year_scen <- function(S, m) {
  sc <- BND$bin_scenarios(S, m)
  nc <- nc_open(BND$bin_file(m)); yr <- as.integer(ncvar_get(nc, "ssp_year")); nc_close(nc)
  # bin_scenarios() names a scenario by "<display model>|<bin lo>". Once every ssp126 year is its
  # own scenario that name is NOT unique (several years share a bin) - rename before anything reads
  # them, same fix amoc_impact_band.R applies for the NAHosMIP hosing years.
  names(sc) <- sprintf("%s|y%d", m, yr)
  mt <- attr(sc, "meta"); mt[, `:=`(id = names(sc), ssp_year = yr)]
  attr(sc, "meta") <- mt
  sc
}

## ---- one model, one branch: replay every ssp126 year, then apply beta straight away -------------
# CHECKPOINTED per (branch, model): a 48h job with no checkpoint is a 48h job that restarts from
# zero on any interruption. Delete the cache file to force a recompute.
cache <- function(branch, m) file.path(d, sprintf("isimip_band_%s_%s.csv", branch, gsub("[^A-Za-z0-9]", "", m)))

band_one <- function(branch, m) {
  f <- cache(branch, m)
  if (file.exists(f)) { cat(sprintf("-- %-11s | %-14s | cached, skipped\n", toupper(branch), m)); return(fread(f)) }
  S  <- eobs_setup(branch)
  sc <- year_scen(S, m); meta <- attr(sc, "meta")
  cat(sprintf("\n-- %s | %s | %d ssp126 years --\n", toupper(branch), m, length(sc)))
  t0  <- Sys.time()
  res <- merge(run_scenarios(S, sc), meta, by = "id", sort = FALSE)[, id := NULL]
  cat(sprintf("   replay %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  eu <- if (branch == "crop") {
    merge(crop_effect(res, KEY), w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "component")][, branch := "crop"][]
  } else if (branch == "crop_gddwin") {
    h <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip)]
    sc2 <- copy(res); setnames(sc2, c("gdd", "heat", "frost", "precip"), paste0(c("gdd", "heat", "frost", "precip"), "_s"))
    z <- merge(sc2, h, by = c("NUTS_ID", "year"))
    z[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
             d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2)]
    ce <- rbindlist(lapply(names(gw_fit), function(k) {
      b <- coef(gw_fit[[k]]); y <- merge(z, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
      for (v in GW_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
      y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
        by = c(KEY, "cntr", "year"), .SDcols = paste0("e_", GW_X)][
        , c(lapply(.SD, mean), .(crop = k)), by = c(KEY, "cntr"), .SDcols = paste0("e_", GW_X)]
    }))
    setnames(ce, paste0("e_", GW_X), GW_X)
    ce <- melt(ce, id.vars = c(KEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
    merge(ce, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "component")][, branch := "crop_gddwin"][]
  } else {
    merge(energy_effect(res, KEY), w_pop, by = c("cntr", "country_id"))[
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "fuel", "component")][
      , branch := paste("energy", fuel)][]
  }
  fwrite(eu, f); eu
}

# crop_gddwin first: cheapest indicator set (no precip-window robustness pass), so a construction
# error surfaces in the shortest possible run before the expensive branches are attempted.
BAND <- rbindlist(lapply(c("crop_gddwin", "crop", "energy"),
                         function(b) rbindlist(lapply(ISIMIP_MODELS, function(m) band_one(b, m)), fill = TRUE)),
                  fill = TRUE)

## ---- disperse across the ssp126 years of each bin ------------------------------------------------
central <- rbind(
  fread(file.path(d, "isimip_bin_impact_crop_eu.csv"))[, .(model, bin_id, component, central = dln, branch = "crop")],
  fread(file.path(d, "isimip_bin_impact_cropgw_eu.csv"))[, .(model, bin_id, component, central = dln, branch = "crop_gddwin")],
  fread(file.path(d, "isimip_bin_impact_energy_eu.csv"))[, .(model, bin_id, component, central = dln,
                                                              branch = paste("energy", fuel))])

band <- BAND[, .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)),
             by = .(branch, model, bin_id, component)]
band <- merge(band, central, by = c("branch", "model", "bin_id", "component"), all.x = TRUE)
fwrite(band[order(branch, model, bin_id, component)], file.path(d, "isimip_band_summary.csv"))

tot  <- BAND[, .(dln = sum(dln)), by = .(branch, model, bin_id, ssp_year)][   # total = sum of components
  , .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)), by = .(branch, model, bin_id)]
ctot <- central[, .(central = sum(central)), by = .(branch, model, bin_id)]
tot  <- merge(tot, ctot, by = c("branch", "model", "bin_id"), all.x = TRUE)
fwrite(tot[order(branch, model, bin_id)], file.path(d, "isimip_band_total.csv"))

pc <- function(x) 100 * (exp(x) - 1)
for (br in sort(unique(tot$branch))) for (m in unique(tot[branch == br]$model)) {
  z <- tot[branch == br & model == m][order(bin_id)]
  if (!nrow(z)) next
  cat(sprintf("\n=== %s | %s | total effect, spread across the bin's ssp126 years ===\n", br, m))
  cat(sprintf("  %7s %4s %10s %10s %9s %10s %10s\n", "bin Sv", "n", "central%", "mean%", "sd(log)%", "min%", "max%"))
  for (i in seq_len(nrow(z)))
    cat(sprintf("  %7.1f %4d %9.2f%% %9.2f%% %8.2f%% %9.2f%% %9.2f%%\n", z$bin_id[i], z$n[i],
                pc(z$central[i]), pc(z$mean[i]), 100 * z$sd[i], pc(z$lo[i]), pc(z$hi[i])))
}
cat("\nwrote isimip_band_{summary,total}.csv\n")
