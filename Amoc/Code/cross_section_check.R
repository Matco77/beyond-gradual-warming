# Cross-sectional check for the extrapolation flagged in Appendix I, section I.4.
# Author: Marco Bova
#
# The benchmark specification (unit FE + year FE + unit-specific quadratic trends) is deliberately
# confounded on the cross-section: it identifies beta from WITHIN-unit year-to-year wobble only,
# which Appendix I showed is a narrow range - the AMOC scenario moves the regressor up to 9.4
# within-sd beyond it. This script asks whether that same slope also holds BETWEEN units, using the
# classic between estimator: collapse each unit to its time-mean and regress the mean outcome on
# the mean regressors, no fixed effects. If the between slope resembles the within slope, the
# within-slope's extrapolation to an out-of-sample regressor level has an empirical basis - two
# largely independent sources of variation (year-to-year wobble vs. permanent cross-unit
# differences) point the same way. If they disagree, the linear form may not be transportable to
# the levels the AMOC scenario reaches, regardless of how precisely the within slope is estimated.
#
# This is diagnostic, not a replacement estimator: the beta applied to the AMOC/ISIMIP scenarios
# are and remain the within (benchmark) coefficients - re-reading this file changes nothing
# upstream. It also does not re-estimate anything used elsewhere in the pipeline.
suppressMessages({library(fixest); library(data.table)})
d  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
RG <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/regression.R")
grab <- function(file, want) {                  # same mechanism as amoc_impact.R: read regression.R's
  e <- new.env(parent = globalenv())            # named definitions without executing the script
  for (x in parse(file))
    if (is.call(x) && identical(as.character(x[[1]]), "<-") && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% want) eval(x, e)
  e
}
R <- grab(RG, c("common_obs", "crop_set", "load_crop", "f_tr"))

## ---- CROP: between-NUTS3, on the benchmark's own estimation sample ---------------------------
CROP_X <- c("gdd_mj", "heat_mj", "frost_mj", "precip_mj", "precip_mj2")
cp9 <- R$load_crop("9.crop_panel_nuts3_estimation.csv")
f_between_crop <- as.formula(paste("ln_yield ~", paste(CROP_X, collapse = " + ")))

cat("=== CROP: within (benchmark) vs between (cross-NUTS3) slope on gdd_mj ===\n")
crop_cmp <- rbindlist(lapply(R$crop_set, function(k) {
  s  <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)     # IDENTICAL sample the benchmark uses
  w  <- feols(R$f_tr, s, cluster = ~cntr, notes = FALSE)      # within (already fitted in amoc_impact.R;
                                                               # refit here so this script is self-contained)
  m  <- s[, lapply(.SD, mean), by = NUTS_ID, .SDcols = c("ln_yield", CROP_X)]   # between estimator:
  b  <- feols(f_between_crop, m, notes = FALSE)                                # one row per unit, no FE
  data.table(crop = k, n_units = nrow(m),
             within_gdd = coef(w)["gdd_mj"], within_se = se(w)["gdd_mj"],
             between_gdd = coef(b)["gdd_mj"], between_se = se(b)["gdd_mj"], between_r2 = r2(b, "r2"))
}))
crop_cmp[, ratio := between_gdd / within_gdd]
print(crop_cmp[, .(crop, n_units, within_gdd = signif(within_gdd,3), within_se = signif(within_se,3),
                   between_gdd = signif(between_gdd,3), between_se = signif(between_se,3),
                   between_r2 = round(between_r2,3), ratio = round(ratio,2))])

## ---- ENERGY: between-country, on the benchmark's own estimation sample -----------------------
ep <- fread(file.path(d, "16.energy_panel_estimation.csv"))
ep[, `:=`(t = year - 2007L, t2 = (year - 2007L)^2)]
FE <- "country_id[t, t2] + year"
efit <- function(dt, wx) {
  f  <- function(fe) as.formula(paste("ln_energy_pc ~", wx, "+ ln_price |", fe))
  fs <- list(f("country_id + year"), f("country_id[t, t2]"), f("country_id[t, t2] + year"))
  s  <- R$common_obs(dt, fs, ~country_id)
  list(w = feols(fs[[3]], s, cluster = ~country_id, notes = FALSE), sample = s)
}
en_cmp <- rbindlist(lapply(list(list("Electricity", "hdd_calendar + cdd_jja", c("hdd_calendar","cdd_jja")),
                                list("Natural gas", "hdd_octmar", "hdd_octmar")), function(z) {
  r  <- efit(ep[fuel == z[[1]]], z[[2]])
  m  <- r$sample[, lapply(.SD, mean), by = country_id, .SDcols = c("ln_energy_pc", "ln_price", z[[3]])]
  fb <- as.formula(paste("ln_energy_pc ~", z[[2]], "+ ln_price"))
  b  <- feols(fb, m, notes = FALSE)
  rbindlist(lapply(z[[3]], function(v) data.table(
    fuel = z[[1]], var = v, n_units = nrow(m),
    within = coef(r$w)[v], within_se = se(r$w)[v],
    between = coef(b)[v], between_se = se(b)[v], between_r2 = r2(b, "r2"))))
}))
en_cmp[, ratio := between / within]
cat("\n=== ENERGY: within (benchmark) vs between (cross-country) slope ===\n")
print(en_cmp[, .(fuel, var, n_units, within = signif(within,3), within_se = signif(within_se,3),
                 between = signif(between,3), between_se = signif(between_se,3),
                 between_r2 = round(between_r2,3), ratio = round(ratio,2))])

fwrite(crop_cmp, file.path(d, "crossection_check_crop.csv"))
fwrite(en_cmp,   file.path(d, "crossection_check_energy.csv"))
cat("\nwrote crossection_check_{crop,energy}.csv\n")
