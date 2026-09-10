# ISIMIP branch, bin-resolved: apply the estimated response functions to the ISIMIP delta_Sv-bin
# scenarios (isimip_bin_fields.R + scenario_replay_isimip_bins.R).
# Author: Marco Bova
#
# Identical logic to isimip_impact.R (single-delta) and to the AMOC branch (Appendix H): effect =
# sum_k beta_k * (X_scen - X_hist), climate regressors only, beta never re-estimated. amoc_impact.R
# and amoc_impact_gddwin.R are sourced (silently) for the fitted coefficients and the fixed
# aggregation weights - crop_effect()/energy_effect() already default to BINKEY = (model, bin_id,
# delta_sv, n_years), which is exactly the scenario file's grouping here, so unlike isimip_impact.R
# (KEY = "model") no key override is needed.
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_gddwin.R"))))
cat("beta, weights and both crop specifications loaded (AMOC central estimate recomputed silently)\n")

pct <- function(x) 100 * (exp(x) - 1)

# emissions scenario, from ISIMIP_SCEN (default ssp126). ssp126 reads/writes bare names (unchanged);
# any other scenario carries a _<scen> suffix on the scenario inputs and the impact outputs.
SCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
SFX  <- if (SCEN == "ssp126") "" else paste0("_", SCEN)

## ---- fixed-window crop (spec A) and energy, at BINKEY resolution --------------------------------
ce_A <- crop_effect(sprintf("scenario_isimip_bins_crop_window%s.csv.gz", SFX))
ee   <- energy_effect(sprintf("scenario_isimip_bins_energy_country%s.csv.gz", SFX))

ce_A_eu <- merge(ce_A, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]
ee_eu <- merge(ee, w_pop, by = c("cntr", "country_id"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "fuel", "component")]

## ---- thermal-window crop (spec C), same construction as isimip_impact.R -------------------------
sc_C <- fread(file.path(d, sprintf("scenario_isimip_bins_crop_gddwin%s.csv.gz", SFX)))
h_C  <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip)]
setnames(sc_C, c("gdd", "heat", "frost", "precip"), paste0(c("gdd", "heat", "frost", "precip"), "_s"))
z_C  <- merge(sc_C, h_C, by = c("NUTS_ID", "year"))
z_C[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
           d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2)]

ce_C <- rbindlist(lapply(names(gw_fit), function(k) {
  b <- coef(gw_fit[[k]]); y <- merge(z_C, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
  for (v in GW_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
  y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
    by = c(BINKEY, "cntr", "year"), .SDcols = paste0("e_", GW_X)][
    , c(lapply(.SD, mean), .(crop = k)), by = c(BINKEY, "cntr"), .SDcols = paste0("e_", GW_X)]
}))
setnames(ce_C, paste0("e_", GW_X), GW_X)
ce_C <- melt(ce_C, id.vars = c(BINKEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
ce_C_eu <- merge(ce_C, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]

fwrite(ce_A,    file.path(d, sprintf("isimip_bin_impact_crop_country%s.csv", SFX)))
fwrite(ce_A_eu, file.path(d, sprintf("isimip_bin_impact_crop_eu%s.csv", SFX)))
fwrite(ce_C,    file.path(d, sprintf("isimip_bin_impact_cropgw_country%s.csv", SFX)))
fwrite(ce_C_eu, file.path(d, sprintf("isimip_bin_impact_cropgw_eu%s.csv", SFX)))
fwrite(ee,      file.path(d, sprintf("isimip_bin_impact_energy_country%s.csv", SFX)))
fwrite(ee_eu,   file.path(d, sprintf("isimip_bin_impact_energy_eu%s.csv", SFX)))

## ---- report: both crop specs, per model x bin, next to the NAHosMIP curve at the same levels ----
cat("\n================ ISIMIP bins: crop, effect on ln(yield), area-weighted Europe ================\n")
cat("A = fixed Mar-Jul (benchmark).  C = thermal-time window. Same caveat as Appendix I: neither\n")
cat("total below is a defensible standalone estimate.\n")
wA <- dcast(ce_A_eu, model + bin_id + delta_sv + n_years ~ component, value.var = "dln")
wC <- dcast(ce_C_eu, model + bin_id + delta_sv + n_years ~ component, value.var = "dln")
compA <- setdiff(names(wA), c(BINKEY)); compC <- setdiff(names(wC), c(BINKEY))
wA[, total := rowSums(.SD), .SDcols = compA]; wC[, total := rowSums(.SD), .SDcols = compC]
for (m in unique(wA$model)) {
  cat(sprintf("\n %s\n", m))
  a <- wA[model == m][order(bin_id)]; c_ <- wC[model == m][order(bin_id)]
  for (i in seq_len(nrow(a))) {
    cat(sprintf("  bin %5.1f Sv (dSv=%+.2f, n=%d)  A TOTAL %+6.2f%%  C TOTAL %+6.2f%%\n",
                a$bin_id[i], a$delta_sv[i], a$n_years[i], pct(a$total[i]), pct(c_$total[i])))
  }
}

cat("\n================ ISIMIP bins: energy, effect on ln(energy per capita), pop-weighted Europe ================\n")
for (f in unique(ee_eu$fuel)) {
  z <- dcast(ee_eu[fuel == f], model + bin_id + delta_sv + n_years ~ component, value.var = "dln")
  comp <- setdiff(names(z), BINKEY); z[, total := rowSums(.SD), .SDcols = comp]
  cat(sprintf("\n %s\n", f))
  for (m in unique(z$model)) {
    zz <- z[model == m][order(bin_id)]
    for (i in seq_len(nrow(zz)))
      cat(sprintf("  %-14s bin %5.1f Sv (n=%d)  TOTAL %+6.2f%%\n", m, zz$bin_id[i], zz$n_years[i], pct(zz$total[i])))
  }
}

cat("\nwrote isimip_bin_impact_{crop,cropgw,energy}_{country,eu}.csv\n")
