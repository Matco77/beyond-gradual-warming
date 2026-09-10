# Phase 3: apply the estimated response functions to the ISIMIP scenario.
# Author: Marco Bova
#
# Identical logic to the AMOC branch (Appendix H): effect = sum_k beta_k * (X_scen - X_hist), climate
# regressors only, beta never re-estimated. Nothing is reimplemented here - amoc_impact.R and
# amoc_impact_gddwin.R are sourced (silently) for the fitted coefficients, the fixed aggregation
# weights, and crop_effect()/energy_effect(), which already accept an arbitrary grouping `key`
# (Appendix H, H.2; used the same way by the uncertainty band, Appendix I). Only the key changes,
# from the AMOC bin key to a single "model" column, because ISIMIP carries one scenario per model,
# not forty bins.
#
# No per-year band on this branch. Phase 2's band disperses the chain across the individual hosing
# years that make up a bin; ISIMIP has no such ensemble - each scenario is a single 30-year
# climatological delta (2071-2100 vs 1985-2014), not a set of years to disperse over. A band here
# would require inter-annual variability WITHIN each 30-year window, which is a different quantity
# and is not computed.
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_gddwin.R"))))
cat("beta, weights and both crop specifications loaded (AMOC central estimate recomputed silently)\n")

KEY <- "model"
pct <- function(x) 100 * (exp(x) - 1)

# emissions scenario, from ISIMIP_SCEN (default ssp126). ssp126 reads/writes bare names (unchanged);
# any other scenario carries a _<scen> suffix on both the scenario inputs and the impact outputs.
SCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
SFX  <- if (SCEN == "ssp126") "" else paste0("_", SCEN)

## ---- fixed-window crop (spec A) and energy -----------------------------------------------------
ce_A <- crop_effect(sprintf("scenario_isimip_crop_window%s.csv", SFX), KEY)
ee   <- energy_effect(sprintf("scenario_isimip_energy_country%s.csv", SFX), KEY)

ce_A_eu <- merge(ce_A, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(KEY, "component")]
ee_eu <- merge(ee, w_pop, by = c("cntr", "country_id"))[
  , .(dln = weighted.mean(dln, w)), by = c(KEY, "fuel", "component")]

## ---- thermal-window crop (spec C), same construction as amoc_impact_gddwin.R --------------------
sc_C <- fread(file.path(d, sprintf("scenario_isimip_crop_gddwin%s.csv", SFX)))
h_C  <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip)]
setnames(sc_C, c("gdd", "heat", "frost", "precip"), paste0(c("gdd", "heat", "frost", "precip"), "_s"))
z_C  <- merge(sc_C, h_C, by = c("NUTS_ID", "year"))
z_C[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
           d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2)]

ce_C <- rbindlist(lapply(names(gw_fit), function(k) {
  b <- coef(gw_fit[[k]]); y <- merge(z_C, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
  for (v in GW_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
  y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
    by = c(KEY, "cntr", "year"), .SDcols = paste0("e_", GW_X)][
    , c(lapply(.SD, mean), .(crop = k)), by = c(KEY, "cntr"), .SDcols = paste0("e_", GW_X)]
}))
setnames(ce_C, paste0("e_", GW_X), GW_X)
ce_C <- melt(ce_C, id.vars = c(KEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
ce_C_eu <- merge(ce_C, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(KEY, "component")]

fwrite(ce_A,    file.path(d, sprintf("isimip_impact_crop_country%s.csv", SFX)))
fwrite(ce_A_eu, file.path(d, sprintf("isimip_impact_crop_eu%s.csv", SFX)))
fwrite(ce_C,    file.path(d, sprintf("isimip_impact_cropgw_country%s.csv", SFX)))
fwrite(ce_C_eu, file.path(d, sprintf("isimip_impact_cropgw_eu%s.csv", SFX)))
fwrite(ee,      file.path(d, sprintf("isimip_impact_energy_country%s.csv", SFX)))
fwrite(ee_eu,   file.path(d, sprintf("isimip_impact_energy_eu%s.csv", SFX)))

## ---- coverage, same restriction as the AMOC branch and reported for the same reason -------------
.scn <- fread(file.path(d, sprintf("scenario_isimip_crop_window%s.csv", SFX)), select = "NUTS_ID")
cat(sprintf("\ncoverage | crop: scenario %d NUTS3 -> effect on %d NUTS3 in %d countries (%.0f%%)\n",
            uniqueN(.scn$NUTS_ID), uniqueN(w_crop$NUTS_ID), uniqueN(ce_A$cntr),
            100 * uniqueN(w_crop$NUTS_ID) / uniqueN(.scn$NUTS_ID)))
cat(sprintf("coverage | energy: effect on %d countries (electricity), %d (gas)\n",
            uniqueN(ee[fuel == "Electricity"]$cntr), uniqueN(ee[fuel == "Natural gas"]$cntr)))

## ---- report: both crop specs side by side, as Appendix I requires -------------------------------
cat("\n================ ISIMIP: crop, effect on ln(yield), area-weighted Europe ================\n")
cat("A = fixed Mar-Jul (benchmark).  C = thermal-time window.\n")
cat("!! Same caveat as the AMOC branch (Appendix I): A is contaminated by phenological\n")
cat("!! misalignment, C cannot identify the thermal dose by construction. Neither total below\n")
cat("!! is a defensible standalone estimate; both are reported so the reader sees the same\n")
cat("!! specification sensitivity here as on the AMOC branch.\n")
wA <- dcast(ce_A_eu, model ~ component, value.var = "dln")
wC <- dcast(ce_C_eu, model ~ component, value.var = "dln")
compA <- setdiff(names(wA), "model"); compC <- setdiff(names(wC), "model")
wA[, total := rowSums(.SD), .SDcols = compA]; wC[, total := rowSums(.SD), .SDcols = compC]
for (m in wA$model) {
  a <- wA[model == m]; c_ <- wC[model == m]
  cat(sprintf("\n %s\n", m))
  cat(sprintf("  A: %s  TOTAL %+.2f%%\n",
              paste(sprintf("%s %+.2f%%", compA, pct(unlist(a[, compA, with = FALSE]))), collapse = "  "),
              pct(a$total)))
  cat(sprintf("  C: %s  TOTAL %+.2f%%\n",
              paste(sprintf("%s %+.2f%%", compC, pct(unlist(c_[, compC, with = FALSE]))), collapse = "  "),
              pct(c_$total)))
}

cat("\n================ ISIMIP: energy, effect on ln(energy per capita), pop-weighted Europe ================\n")
for (f in unique(ee_eu$fuel)) {
  z <- dcast(ee_eu[fuel == f], model ~ component, value.var = "dln")
  comp <- setdiff(names(z), "model"); z[, total := rowSums(.SD), .SDcols = comp]
  cat(sprintf("\n %s\n", f))
  for (i in seq_len(nrow(z))) cat(sprintf("  %-14s %s  TOTAL %+.2f%%\n", z$model[i],
    paste(sprintf("%s %+.2f%%", comp, pct(unlist(z[i, comp, with = FALSE]))), collapse = "  "),
    pct(z$total[i])))
}

cat("\nwrote isimip_impact_{crop,cropgw,energy}_{country,eu}.csv\n")
