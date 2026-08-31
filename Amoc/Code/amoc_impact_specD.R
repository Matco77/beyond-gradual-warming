# Crop impact, spec D: the thermal-time window plus window LENGTH as a regressor.
# Author: Marco Bova
#
# Appendix I, §I.6-I.7 identified why neither existing specification identifies a crop effect:
# spec A (fixed Mar-Jul) confounds the thermal dose with phenological misalignment (a warm year
# samples a later growth stage inside a window that does not move); spec C (thermal-time window)
# cannot identify the dose at all, because the window is DEFINED by accumulating a fixed amount of
# thermal time, so gdd inside it is ~constant by construction (within-sd of gdd collapses to 0.0% of
# its raw variance under spec C - Appendix I, §I.6) and the fixed effects absorb what little is
# left. Appendix I, §I.7 found where the AMOC signal actually goes once the window is spec C's: not
# into the thermal dose, but into how LONG the season took to accumulate it - n_day, which the
# thermal-window construction already computes (57-322 days historically, real variation) but which
# neither A nor C carries as a regressor.
#
# Spec D is spec C's five regressors plus n_day. It is not offered as the final answer - it is the
# specification Appendix I said was missing, estimated for the first time here.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact_gddwin.R"))

D_X <- c(GW_X, "n_day")                                     # gdd, heat, frost, precip, precip2, n_day
f_D <- as.formula(paste("ln_yield ~", paste(D_X, collapse = " + "), "|", FE))

## ---- beta on the historical panel, IDENTICAL sample construction as spec C ---------------------
d_fit <- lapply(setNames(R$crop_set, R$crop_set), function(k) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  z <- merge(s[, .(NUTS_ID, cntr, year, ln_yield, t, t2)], hist_gw, by = c("NUTS_ID", "year"))
  feols(f_D, R$common_obs(z, list(f_D), ~cntr), cluster = ~cntr, notes = FALSE)
})

## ---- scenario minus historical, per component (n_day: plain difference, not a squared term) ----
sc_D <- fread(file.path(d, "scenario_bins_crop_gddwin.csv.gz"))
h_D  <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip, n_day)]
setnames(sc_D, c("gdd", "heat", "frost", "precip", "n_day"), paste0(c("gdd", "heat", "frost", "precip", "n_day"), "_s"))
zD <- merge(sc_D, h_D, by = c("NUTS_ID", "year"))
zD[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
          d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2, d_n_day = n_day_s - n_day)]

ceD <- rbindlist(lapply(names(d_fit), function(k) {
  b <- coef(d_fit[[k]]); y <- merge(zD, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
  for (v in D_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
  y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
    by = c(BINKEY, "cntr", "year"), .SDcols = paste0("e_", D_X)][
    , c(lapply(.SD, mean), .(crop = k)), by = c(BINKEY, "cntr"), .SDcols = paste0("e_", D_X)]
}))
setnames(ceD, paste0("e_", D_X), D_X)
ceD <- melt(ceD, id.vars = c(BINKEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
ceD_eu <- merge(ceD, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]
fwrite(ceD,    file.path(d, "amoc_impact_specD_country.csv"))
fwrite(ceD_eu, file.path(d, "amoc_impact_specD_eu.csv"))

## ---- report: A, C, D side by side --------------------------------------------------------------
D <- ceD_eu[, .(dlnD = sum(dln)), by = .(model, bin_id, n_years)]
cmpD <- merge(cmp, D, by = c("model", "bin_id", "n_years"))

cat("\n============ CROP: spec A vs spec C vs spec D (thermal window + n_day), Europe ============\n")
cat("A = fixed Mar-Jul (benchmark).  C = thermal window.  D = C + window LENGTH as a regressor.\n")
for (m in unique(cmpD$model)) {
  cat(sprintf("\n %s\n", m))
  zz <- cmpD[model == m][order(bin_id)]
  cat(sprintf("  %7s %4s %10s %10s %10s %9s\n", "bin Sv", "n", "A total%", "C total%", "D total%", "win days"))
  for (i in seq_len(nrow(zz)))
    cat(sprintf("  %7.1f %4d %9.2f%% %9.2f%% %9.2f%% %9.0f\n", zz$bin_id[i], zz$n_years[i],
                pct(zz$dlnA[i]), pct(zz$dlnC[i]), pct(zz$dlnD[i]), zz$n_day[i]))
}

cat("\n-- beta, spec D (historical panel) --\n")
for (k in names(d_fit)) cat(sprintf("  %-14s N=%6d  %s\n", k, nobs(d_fit[[k]]),
  paste(sprintf("%s=%+.3e(t%+.2f)", D_X, coef(d_fit[[k]])[D_X],
                coef(d_fit[[k]])[D_X] / se(d_fit[[k]])[D_X]), collapse = "  ")))

## ---- identifying variation check, same test as Appendix I, §I.6 on spec C ----------------------
cat("\n-- within-sd of n_day, spec D (share of raw variance surviving the fixed effects) --\n")
for (k in R$crop_set) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  z <- R$common_obs(merge(s[, .(NUTS_ID, cntr, year, ln_yield, t, t2)], hist_gw, by = c("NUTS_ID", "year")),
                    list(f_D), ~cntr)
  w <- resid(feols(as.formula(paste("n_day ~ 1 |", FE)), z, notes = FALSE))
  cat(sprintf("  %-14s within sd %5.1f days (%.1f%% of raw variance) | mean n_day %.0f\n",
              k, sd(w), 100 * var(w) / var(z$n_day), mean(z$n_day)))
}
cat("\nwrote amoc_impact_specD_{country,eu}.csv\n")
