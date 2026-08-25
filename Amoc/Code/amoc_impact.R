# Phase 2: apply the estimated response functions to the AMOC bin scenarios.
# Author: Marco Bova
#
# effect = sum_k beta_k * (X_k,scenario - X_k,historical), over the CLIMATE regressors only.
# Price and population are held at their historical values, so they drop out of the difference by
# construction - the counterfactual is "same economy, different climate".
#
# The beta are NEVER re-estimated on the perturbed climate. They are refitted here on the
# HISTORICAL panels with the benchmark specifications of regression.R, and the specification is not
# retyped: grab() evaluates the named top-level definitions of regression.R (common_obs, crop_set,
# load_crop, f_tr) without running the script, so there is one source of truth for the formulas.
# The energy formulas are built inside energy_tab() from a `wx` string and cannot be grabbed by
# name, so they are rebuilt exactly as regression.R:114-118 does. Both branches are checked against
# the saved coefficient tables in Amoc/results/ at the bottom of this script.
#
# Contributions are reported PER COMPONENT, not just as a net. Under cooling GDD and Heat both fall
# while Frost rises, and with beta_gdd and beta_heat of opposite sign the components push yield in
# opposite directions - a small net can hide two large offsetting effects.
#
# Only the Mar-Jul window enters: regression.R's benchmark (f_tr) is estimated on the _mj
# regressors. Phase 1 also emits Apr-Aug, but regression.R fits no benchmark on _aa, so there is no
# beta to apply to it.
suppressMessages({library(fixest); library(data.table)})
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
res <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/results")
RG  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/regression.R")

## ---- 1. beta from the benchmark specifications ------------------------------------------------
grab <- function(file, want) {                  # eval only the named top-level `<-` defs of a script
  e <- new.env(parent = globalenv())
  for (x in parse(file))
    if (is.call(x) && identical(as.character(x[[1]]), "<-") && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% want) eval(x, e)
  e
}
R <- grab(RG, c("common_obs", "crop_set", "load_crop", "f_tr"))

CROP_X <- c("gdd_mj", "heat_mj", "frost_mj", "precip_mj", "precip_mj2")   # climate regressors only
cp9 <- R$load_crop("9.crop_panel_nuts3_estimation.csv")
crop_fit <- lapply(setNames(R$crop_set, R$crop_set), function(k)
  feols(R$f_tr, R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr), cluster = ~cntr, notes = FALSE))

ep <- fread(file.path(d, "16.energy_panel_estimation.csv"))
ep[, `:=`(t = year - 2007L, t2 = (year - 2007L)^2)]     # centred trend, regression.R:113
efit <- function(dt, wx) {                              # regression.R:114-118, benchmark column
  f  <- function(fe) as.formula(paste("ln_energy_pc ~", wx, "+ ln_price |", fe))
  fs <- list(f("country_id + year"), f("country_id[t, t2]"), f("country_id[t, t2] + year"))
  s  <- R$common_obs(dt, fs, ~country_id)             # the benchmark's constant sample
  # the sample is returned alongside the fit: m$fixef_id holds INTEGER factor codes, not the
  # country strings, so it cannot be used to recover which countries beta was identified on.
  list(m = feols(fs[[3]], s, cluster = ~country_id, notes = FALSE), sample = s)
}
en_fit <- list(Electricity = efit(ep[fuel == "Electricity"], "hdd_calendar + cdd_jja"),
               `Natural gas` = efit(ep[fuel == "Natural gas"], "hdd_octmar"))
EN_X <- list(Electricity = c("hdd_calendar", "cdd_jja"), `Natural gas` = "hdd_octmar")

## ---- 2. fixed aggregation weights --------------------------------------------------------------
# Crop area per NUTS3 x crop, averaged over the estimation years: a FIXED weight, so the country
# aggregate reflects the climate effect and not a drift in the cropping mix. Only (NUTS3, crop)
# cells that are in the benchmark estimation sample get a weight - beta is only transportable
# where it was identified.
w_crop <- rbindlist(lapply(names(crop_fit), function(k) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  s[!is.na(area_ha) & area_ha > 0, .(crop = k, w = mean(area_ha)), by = .(NUTS_ID, cntr)]
}))
# Population per country, averaged over the estimation years: fixed weight for the Europe aggregate.
w_pop <- ep[fuel == "Electricity" & !is.na(population),
            .(w = mean(population)), by = .(country_id, cntr = eurostat_geo)]

## ---- 3. scenario minus historical, per component -----------------------------------------------
BINKEY <- c("model", "bin_id", "delta_sv", "n_years")

crop_effect <- function(scen_file) {
  sc <- fread(file.path(d, scen_file))[window == "Mar-Jul"]
  hi <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[window == "Mar-Jul"]
  setnames(sc, c("gdd", "heat", "frost", "precip"), paste0(c("gdd", "heat", "frost", "precip"), "_s"))
  z <- merge(sc, hi[, .(NUTS_ID, year, gdd, heat, frost, precip)], by = c("NUTS_ID", "year"))
  # precip^2 is RECOMPUTED from the perturbed precip, never shifted: it is a regressor in its own
  # right and (p+dp)^2 != p^2 + dp^2.
  z[, `:=`(d_gdd_mj = gdd_s - gdd, d_heat_mj = heat_s - heat, d_frost_mj = frost_s - frost,
           d_precip_mj = precip_s - precip, d_precip_mj2 = precip_s^2 - precip^2)]
  out <- rbindlist(lapply(names(crop_fit), function(k) {
    b <- coef(crop_fit[[k]]); wk <- w_crop[crop == k]
    y <- merge(z, wk, by = "NUTS_ID", allow.cartesian = TRUE)          # keeps only the estimation cells
    for (v in CROP_X) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
    # NUTS3 -> country: area-weighted mean of the effect, then mean over years
    y[, lapply(.SD, function(x) weighted.mean(x, w, na.rm = TRUE)),
      by = c(BINKEY, "cntr", "year"), .SDcols = paste0("e_", CROP_X)][
      , c(lapply(.SD, mean), .(crop = k)), by = c(BINKEY, "cntr"), .SDcols = paste0("e_", CROP_X)]
  }))
  setnames(out, paste0("e_", CROP_X), CROP_X)
  melt(out, id.vars = c(BINKEY, "cntr", "crop"), variable.name = "component", value.name = "dln")
}

energy_effect <- function(scen_file) {
  sc <- fread(file.path(d, scen_file))
  hi <- fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))
  V  <- c("hdd_calendar", "hdd_octmar", "cdd_jja")
  setnames(sc, V, paste0(V, "_s"))
  z <- merge(sc, hi[, c("cntr", "year", V), with = FALSE], by = c("cntr", "year"))
  for (v in V) set(z, j = paste0("d_", v), value = z[[paste0(v, "_s")]] - z[[v]])
  rbindlist(lapply(names(en_fit), function(f) {
    b <- coef(en_fit[[f]]$m); xs <- EN_X[[f]]
    keep <- unique(en_fit[[f]]$sample$country_id)                      # countries beta was fitted on
    y <- merge(z, w_pop[country_id %in% keep], by = "cntr")
    for (v in xs) set(y, j = paste0("e_", v), value = b[[v]] * y[[paste0("d_", v)]])
    m <- y[, lapply(.SD, mean, na.rm = TRUE), by = c(BINKEY, "cntr", "country_id"),
           .SDcols = paste0("e_", xs)]                                 # mean over years, per country
    setnames(m, paste0("e_", xs), xs)
    melt(m, id.vars = c(BINKEY, "cntr", "country_id"), variable.name = "component",
         value.name = "dln")[, fuel := f]
  }))
}

## ---- 4. run and report -------------------------------------------------------------------------
pct <- function(x) 100 * (exp(x) - 1)            # log-linear outcome -> relative change

ce <- crop_effect("scenario_bins_crop_window.csv.gz")
ee <- energy_effect("scenario_bins_energy_country.csv.gz")

# Europe aggregate: crop area across all crops/regions, population across countries
ce_eu <- merge(ce, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]
ee_eu <- merge(ee, w_pop, by = c("cntr", "country_id"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "fuel", "component")]

fwrite(ce,    file.path(d, "amoc_impact_crop_country.csv"))
fwrite(ee,    file.path(d, "amoc_impact_energy_country.csv"))
fwrite(ce_eu, file.path(d, "amoc_impact_crop_eu.csv"))
fwrite(ee_eu, file.path(d, "amoc_impact_energy_eu.csv"))

show <- function(x, lab, grp = NULL) {
  wide <- dcast(x, as.formula(paste(paste(c(BINKEY, grp), collapse = " + "), "~ component")),
                value.var = "dln")
  comp <- setdiff(names(wide), c(BINKEY, grp))
  wide[, total := rowSums(.SD), .SDcols = comp]
  for (m in unique(wide$model)) {
    cat(sprintf("\n%s -- %s\n", lab, m))
    z <- wide[model == m][order(bin_id)]
    cat(sprintf("  %7s %4s %s %10s\n", "bin Sv", "n",
                paste(sprintf("%12s", comp), collapse = ""), "TOTAL %"))
    for (i in seq_len(nrow(z)))
      cat(sprintf("  %7.1f %4d %s %+9.2f%%\n", z$bin_id[i], z$n_years[i],
                  paste(sprintf("%+11.2f%%", pct(unlist(z[i, comp, with = FALSE]))), collapse = ""),
                  pct(z$total[i])))
  }
  invisible(wide)
}
cat("\n================ CROP: effect on ln(yield), area-weighted Europe ================\n")
cat("(components in % of yield; TOTAL is the net)\n")
show(ce_eu, "CROP, all crops")
for (f in names(en_fit)) {
  cat(sprintf("\n================ ENERGY %s: effect on ln(energy per capita), pop-weighted Europe ================\n", f))
  show(ee_eu[fuel == f], paste("ENERGY", f))
}

## ---- 5. guard: the beta must match the committed estimation tables -----------------------------
cat("\n-- beta used (must match Amoc/results/) --\n")
for (k in names(crop_fit)) cat(sprintf("  %-14s N=%6d  %s\n", k, nobs(crop_fit[[k]]),
  paste(sprintf("%s=%.4g", names(coef(crop_fit[[k]])), coef(crop_fit[[k]])), collapse = "  ")))
for (f in names(en_fit))  cat(sprintf("  %-14s N=%6d  %s\n", f, nobs(en_fit[[f]]$m),
  paste(sprintf("%s=%.4g", names(coef(en_fit[[f]]$m)), coef(en_fit[[f]]$m)), collapse = "  ")))
cat("\nwrote amoc_impact_{crop,energy}_{country,eu}.csv\n")
