# Crop impact under the GDD-WINDOW specification, reported next to the fixed Mar-Jul one.
# Author: Marco Bova
#
# Why both. beta_gdd < 0 - which is what makes the fixed-window result say that a weaker AMOC RAISES
# European yields - does not survive moving the window from the calendar to accumulated thermal
# time: soft wheat loses significance, spring barley changes sign, winter barley collapses to zero,
# only durum survives. The 28 C cap is inert (verified in crop_gdd_diagnostics.R); the window is
# what does the work. A fixed calendar window in a warm year samples a later phenological stage, so
# "more GDD in Mar-Jul" partly measures window misalignment rather than a thermal dose-response.
# Neither specification is obviously the truth, so both are reported and the reader is told the sign
# is a specification choice, not a finding.
#
# Anchors do NOT move with the scenario: g_open/g_close stand for the crop's thermal requirement,
# which does not change because the climate does. Under strong cooling the close anchor can become
# unreachable - the crop never completes its cycle inside the year. That is reported as a result
# (share of region-years that fail to close, per bin), not patched away.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))

GW_X <- c("gdd", "heat", "frost", "precip", "precip2")
FE   <- "NUTS_ID[t, t2] + year"
f_gw <- as.formula(paste("ln_yield ~", paste(GW_X, collapse = " + "), "|", FE))

## ---- beta on the historical GDD-window panel ---------------------------------------------------
hist_gw <- fread(file.path(d, "11.crop_weather_gdd_window.csv"))[, precip2 := precip^2][]
gw_fit <- lapply(setNames(R$crop_set, R$crop_set), function(k) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)          # same yield sample as the benchmark
  z <- merge(s[, .(NUTS_ID, cntr, year, ln_yield, t, t2)], hist_gw, by = c("NUTS_ID", "year"))
  feols(f_gw, R$common_obs(z, list(f_gw), ~cntr), cluster = ~cntr, notes = FALSE)
})

## ---- scenario minus historical, per component --------------------------------------------------
sc <- fread(file.path(d, "scenario_bins_crop_gddwin.csv.gz"))
h  <- hist_gw[, .(NUTS_ID, year, gdd, heat, frost, precip)]
setnames(sc, c("gdd", "heat", "frost", "precip"), paste0(c("gdd", "heat", "frost", "precip"), "_s"))
z <- merge(sc, h, by = c("NUTS_ID", "year"))
z[, `:=`(d_gdd = gdd_s - gdd, d_heat = heat_s - heat, d_frost = frost_s - frost,
         d_precip = precip_s - precip, d_precip2 = precip_s^2 - precip^2)]   # precip^2 recomputed

ce <- rbindlist(lapply(names(gw_fit), function(k) {
  b <- coef(gw_fit[[k]]); y <- merge(z, w_crop[crop == k], by = "NUTS_ID", allow.cartesian = TRUE)
  for (v in GW_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
  y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
    by = c(BINKEY, "cntr", "year"), .SDcols = paste0("e_", GW_X)][
    , c(lapply(.SD, mean), .(crop = k)), by = c(BINKEY, "cntr"), .SDcols = paste0("e_", GW_X)]
}))
setnames(ce, paste0("e_", GW_X), GW_X)
ce <- melt(ce, id.vars = c(BINKEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
ce_eu <- merge(ce, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]
fwrite(ce,    file.path(d, "amoc_impact_cropgw_country.csv"))
fwrite(ce_eu, file.path(d, "amoc_impact_cropgw_eu.csv"))

## ---- cycles that never close -------------------------------------------------------------------
# Share of region-years whose accumulated thermal time never reaches the close anchor. In the
# historical record this is ~0 by construction (the anchors are its own 1990-2010 mean); wherever it
# rises under a scenario, the crop is not completing its cycle inside the calendar year - a bigger
# statement than any coefficient in the table above.
nocl <- sc[, .(pct_open = 100 * mean(closed == 0L), n_day = mean(n_day)), by = .(model, bin_id, delta_sv)]

## ---- report ------------------------------------------------------------------------------------
pct <- function(x) 100 * (exp(x) - 1)
A <- fread(file.path(d, "amoc_impact_crop_eu.csv"))[, .(dlnA = sum(dln)), by = .(model, bin_id, n_years)]
C <- ce_eu[, .(dlnC = sum(dln)), by = .(model, bin_id, n_years)]
cmp <- merge(merge(A, C, by = c("model", "bin_id", "n_years")), nocl, by = c("model", "bin_id"))

cat("\n============ CROP: fixed Mar-Jul vs GDD-window, Europe, area-weighted ============\n")
cat("A = benchmark (capped GDD, fixed Mar-Jul).  C = uncapped GDD on the thermal-time window.\n")
cat("cycle-open % = region-years whose thermal time never reaches the close anchor.\n")
for (m in unique(cmp$model)) {
  cat(sprintf("\n %s\n", m))
  zz <- cmp[model == m][order(bin_id)]
  cat(sprintf("  %7s %4s %11s %11s %11s %9s\n", "bin Sv", "n", "A total%", "C total%", "cycle-open%", "win days"))
  for (i in seq_len(nrow(zz)))
    cat(sprintf("  %7.1f %4d %10.2f%% %10.2f%% %10.1f%% %9.0f\n", zz$bin_id[i], zz$n_years[i],
                pct(zz$dlnA[i]), pct(zz$dlnC[i]), zz$pct_open[i], zz$n_day[i]))
}

cat("\n-- beta, GDD-window specification (historical panel) --\n")
for (k in names(gw_fit)) cat(sprintf("  %-14s N=%6d  %s\n", k, nobs(gw_fit[[k]]),
  paste(sprintf("%s=%+.3e(t%+.2f)", GW_X, coef(gw_fit[[k]])[GW_X],
                coef(gw_fit[[k]])[GW_X] / se(gw_fit[[k]])[GW_X]), collapse = "  ")))
cat("\nwrote amoc_impact_cropgw_{country,eu}.csv\n")
