# ISIMIP crop impact under spec D (thermal window + window length), symmetric with Appendix J's A/C.
# Author: Marco Bova
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_specD.R"))

sc_D_iso <- fread(file.path(d, "scenario_isimip_crop_gddwin.csv"))
h_D_iso  <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip, n_day)]
setnames(sc_D_iso, c("gdd","heat","frost","precip","n_day"), paste0(c("gdd","heat","frost","precip","n_day"), "_s"))
zDi <- merge(sc_D_iso, h_D_iso, by = c("NUTS_ID", "year"))
zDi[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
           d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2, d_n_day = n_day_s - n_day)]

ceD_iso <- rbindlist(lapply(names(d_fit), function(k) {
  b <- coef(d_fit[[k]]); y <- merge(zDi, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
  for (v in D_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
  y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
    by = c("model", "cntr", "year"), .SDcols = paste0("e_", D_X)][
    , c(lapply(.SD, mean), .(crop = k)), by = c("model", "cntr"), .SDcols = paste0("e_", D_X)]
}))
setnames(ceD_iso, paste0("e_", D_X), D_X)
ceD_iso <- melt(ceD_iso, id.vars = c("model", "cntr", "crop"), variable.name = "component", value.name = "dln")
ceD_iso_eu <- merge(ceD_iso, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c("model", "component")]
fwrite(ceD_iso,    file.path(d, "isimip_impact_specD_country.csv"))
fwrite(ceD_iso_eu, file.path(d, "isimip_impact_specD_eu.csv"))

## ---- report: A, C, D side by side for ISIMIP ----------------------------------------------------
A_iso <- fread(file.path(d, "isimip_impact_crop_eu.csv"))[, .(dlnA = sum(dln)), by = model]
C_iso <- fread(file.path(d, "isimip_impact_cropgw_eu.csv"))[, .(dlnC = sum(dln)), by = model]
D_iso <- ceD_iso_eu[, .(dlnD = sum(dln)), by = model]
cmp_iso <- merge(merge(A_iso, C_iso, by = "model"), D_iso, by = "model")
cat("\n============ ISIMIP: crop, spec A vs spec C vs spec D, Europe ============\n")
for (i in seq_len(nrow(cmp_iso)))
  cat(sprintf("  %-14s A %+7.2f%%   C %+7.2f%%   D %+7.2f%%\n", cmp_iso$model[i],
              pct(cmp_iso$dlnA[i]), pct(cmp_iso$dlnC[i]), pct(cmp_iso$dlnD[i])))
cat("\nwrote isimip_impact_specD_{country,eu}.csv\n")
