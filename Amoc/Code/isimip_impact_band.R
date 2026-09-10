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
# crop_effect()/energy_effect()/gw_fit as isimip_bin_impact_bins.R's central estimate - only the
# grouping key changes (bin_id -> bin_id + ssp_year), so the band cannot silently disagree with the
# number it is a band around.
#
# The two crop branches are kept PER CROP (Soft/Durum wheat, Spring/Winter barley) all the way to
# the output; energy has no crop split. The all-crop number of the central estimate is the
# area-weighted mean of the per-crop rows, so nothing is lost by carrying the split.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_engine.R"))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_gddwin.R"))))
cat("beta, weights and both crop specifications loaded (AMOC central estimate recomputed silently)\n")

ISIMIP_MODELS <- c("ipsl-cm6a-lr", "ec-earth3")     # file ids, as scenario_replay_isimip_bins.R uses
KEY <- c(BINKEY, "ssp_year")
# emissions scenario, from ISIMIP_SCEN (default ssp126). ssp126 reads/writes bare names (unchanged);
# any other scenario carries a _<scen> suffix on the per-year field input, the per-crop central
# inputs, the checkpoint caches, and the band outputs.
SCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
SFX  <- if (SCEN == "ssp126") "" else paste0("_", SCEN)
BND <- grab(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/scenario_replay_isimip_bins.R"),
            c("bin_file", "bin_scenarios"))
BND$bin_file <- function(m) file.path(d, sprintf("isimip_year_fields_%s%s.nc", m, SFX))   # per-year fields

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
cache <- function(branch, m) file.path(d, sprintf("isimip_band_%s_%s%s.csv", branch, gsub("[^A-Za-z0-9]", "", m), SFX))

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
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "component", "crop")][, branch := "crop"][]
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
      , .(dln = weighted.mean(dln, w)), by = c(KEY, "component", "crop")][, branch := "crop_gddwin"][]
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

# ALL-CROP band, reconstructed from the per-crop BAND above - no second replay needed. The EU
# aggregation inside band_one() is a weighted mean over (cntr, crop) cells; a weighted mean is
# linear, so summing per-crop EU numbers with weight W_k = total area of crop k (sum over cntr) is
# EXACTLY the same all-crop EU number crop_effect()/amoc_impact.R would produce directly - checked
# algebraically, not assumed. crop = NA marks "no split", same convention energy already used before
# the split existed.
# crop = "ALL" (not NA) marks the all-crop row: NA does not survive the fwrite/fread round-trip
# (comes back as ""), so a downstream is.na() filter would silently break.
Wk        <- w_crop[, .(W = sum(w)), by = crop]
all_crop  <- merge(BAND[!is.na(crop)], Wk, by = "crop")[
  , .(dln = weighted.mean(dln, W)), by = .(branch, model, bin_id, delta_sv, n_years, ssp_year, component)][
  , crop := "ALL"]
BAND <- rbindlist(list(BAND, all_crop), use.names = TRUE, fill = TRUE)
BAND[is.na(crop), crop := "ALL"]   # energy rows: never split, so they already ARE the all-crop row

## ---- disperse across the ssp126 years of each bin ------------------------------------------------
# central estimate, kept PER CROP for the two crop branches (energy has no crop split), PLUS the
# all-crop row (crop = NA) read straight from the *_eu.csv files - same reason as the BAND
# reconstruction above, but here the all-crop file already exists, no need to re-derive it.
wc <- w_crop[, .(w = sum(w)), by = .(cntr, crop)]
cpc <- function(file, br) merge(fread(file.path(d, file)), wc, by = c("cntr", "crop"))[
  , .(central = weighted.mean(dln, w)), by = .(model, bin_id, crop, component)][, branch := br][]
allcrop <- function(file, br) fread(file.path(d, file))[
  , .(model, bin_id, crop = "ALL", component, central = dln, branch = br)]
central <- rbind(
  cpc(sprintf("isimip_bin_impact_crop_country%s.csv", SFX),   "crop"),
  cpc(sprintf("isimip_bin_impact_cropgw_country%s.csv", SFX), "crop_gddwin"),
  allcrop(sprintf("isimip_bin_impact_crop_eu%s.csv", SFX),   "crop"),
  allcrop(sprintf("isimip_bin_impact_cropgw_eu%s.csv", SFX), "crop_gddwin"),
  fread(file.path(d, sprintf("isimip_bin_impact_energy_eu%s.csv", SFX)))[
    , .(model, bin_id, crop = "ALL", component, central = dln, branch = paste("energy", fuel))],
  use.names = TRUE)

band <- BAND[, .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)),
             by = .(branch, model, bin_id, component, crop)]
band <- merge(band, central, by = c("branch", "model", "bin_id", "component", "crop"), all.x = TRUE)
fwrite(band[order(branch, model, bin_id, crop, component)], file.path(d, sprintf("isimip_band_summary%s.csv", SFX)))

tot  <- BAND[, .(dln = sum(dln)), by = .(branch, model, bin_id, crop, ssp_year)][   # total = sum of components
  , .(n = .N, mean = mean(dln), sd = sd(dln), lo = min(dln), hi = max(dln)), by = .(branch, model, bin_id, crop)]
ctot <- central[, .(central = sum(central)), by = .(branch, model, bin_id, crop)]
tot  <- merge(tot, ctot, by = c("branch", "model", "bin_id", "crop"), all.x = TRUE)
fwrite(tot[order(branch, model, bin_id, crop)], file.path(d, sprintf("isimip_band_total%s.csv", SFX)))

pc <- function(x) 100 * (exp(x) - 1)
for (br in sort(unique(tot$branch))) for (m in unique(tot[branch == br]$model))
  for (cr in unique(tot[branch == br & model == m]$crop)) {
  z <- tot[branch == br & model == m & crop == cr][order(bin_id)]
  if (!nrow(z)) next
  cat(sprintf("\n=== %s | %s | %s | total effect, spread across the bin's ssp126 years ===\n",
              br, m, if (cr == "ALL") "all crops" else cr))
  cat(sprintf("  %7s %4s %10s %10s %9s %10s %10s\n", "bin Sv", "n", "central%", "mean%", "sd(log)%", "min%", "max%"))
  for (i in seq_len(nrow(z)))
    cat(sprintf("  %7.1f %4d %9.2f%% %9.2f%% %8.2f%% %9.2f%% %9.2f%%\n", z$bin_id[i], z$n[i],
                pc(z$central[i]), pc(z$mean[i]), 100 * z$sd[i], pc(z$lo[i]), pc(z$hi[i])))
}
cat("\nwrote isimip_band_{summary,total}.csv (per crop for the crop branches)\n")
