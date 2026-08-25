# Is beta_gdd < 0 a phenological effect, or is GDD absorbing something else?
# Author: Marco Bova
#
# The benchmark gives beta_gdd < 0 for all four crops, which drives the whole "AMOC weakening raises
# European yields" result. The wild-cluster bootstrap says that coefficient is PRECISELY estimated;
# it says nothing about whether it is the causal effect of growing degree days. A negative GDD
# coefficient on C3 cereals is against prior: the 28 C cap exists precisely to separate useful
# warmth from heat stress, and heat stress is already carried by heat_mj.
#
# Three candidate explanations, tested here:
#   (1) COLLINEARITY. gdd and heat both rise with temperature. With heat_mj insignificant, the model
#       may be loading the whole thermal signal onto gdd. Tested by the WITHIN correlation (after
#       partialling out the benchmark's fixed effects - the raw correlation is not the relevant one,
#       identification is within-region year-to-year) and by dropping heat from the spec.
#   (2) HETEROGENEITY BY CLIMATE. NUTS3 fixed effects absorb the level of gdd, so latitude cannot
#       confound the intercept - but a single pooled beta is still a weighted average of region-
#       specific slopes. If warmth helps in Denmark and hurts in Sicily, the pooled number can come
#       out negative without describing anywhere. Tested by refitting within terciles of regional
#       mean gdd.
#   (3) FIXED WINDOW. Mar-Jul does not follow phenology, and phenology shifts with temperature -
#       so in a warm year the window covers a different crop stage. Tested with the GDD-window of
#       gdd_window_daily.R (11.crop_weather_gdd_window.csv), whose edges are set by cumulative GDD
#       and which ALSO drops the 28 C cap: pmax(tg-5,0) instead of pmax(pmin(tg,28)-5,0).
#
# Nothing in regression.R is modified: its definitions are read with the same grab() used in
# amoc_impact.R, and every alternative is fitted here.
suppressMessages({library(fixest); library(data.table)})
d  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
RG <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/regression.R")
grab <- function(file, want) {
  e <- new.env(parent = globalenv())
  for (x in parse(file))
    if (is.call(x) && identical(as.character(x[[1]]), "<-") && is.name(x[[2]]) &&
        as.character(x[[2]]) %in% want) eval(x, e)
  e
}
R   <- grab(RG, c("common_obs", "crop_set", "load_crop", "f_tr"))
cp9 <- R$load_crop("9.crop_panel_nuts3_estimation.csv")
FE  <- "NUTS_ID[t, t2] + year"

## ---- (1a) how collinear are gdd and heat, WITHIN the benchmark's fixed effects? ---------------
cat("=== gdd_mj vs heat_mj correlation ===\n")
cat("  raw = across all region-years; within = after partialling out", FE, "\n")
for (k in R$crop_set) {
  s  <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  rr <- function(v) resid(feols(as.formula(sprintf("%s ~ 1 | %s", v, FE)), s, notes = FALSE))
  g <- rr("gdd_mj"); h <- rr("heat_mj")
  # variance inflation of gdd from the other four regressors, in the within transform
  X <- sapply(c("heat_mj", "frost_mj", "precip_mj", "precip_mj2"), rr)
  r2 <- summary(lm(g ~ X))$r.squared
  cat(sprintf("  %-14s raw %+.3f | within %+.3f | VIF(gdd|others) %.2f\n",
              k, cor(s$gdd_mj, s$heat_mj), cor(g, h), 1 / (1 - r2)))
}

## ---- (1b) does beta_gdd move when heat is dropped? --------------------------------------------
cat("\n=== beta_gdd with and without heat_mj in the spec ===\n")
f_noheat <- as.formula(paste("ln_yield ~ gdd_mj + frost_mj + precip_mj + precip_mj2 |", FE))
for (k in R$crop_set) {
  s  <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  b1 <- coef(feols(R$f_tr,   s, cluster = ~cntr, notes = FALSE))["gdd_mj"]
  b2 <- coef(feols(f_noheat, s, cluster = ~cntr, notes = FALSE))["gdd_mj"]
  cat(sprintf("  %-14s benchmark %+.3e | heat dropped %+.3e | change %+.1f%%\n",
              k, b1, b2, 100 * (b2 / b1 - 1)))
}

## ---- (2) is the pooled beta hiding opposite slopes across climates? ---------------------------
cat("\n=== beta_gdd by tercile of regional mean gdd_mj (cold / mid / warm) ===\n")
for (k in R$crop_set) {
  s <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  s[, gbar := mean(gdd_mj), by = NUTS_ID]
  s[, terc := cut(gbar, quantile(unique(data.table(NUTS_ID, gbar))$gbar, 0:3/3),
                  labels = c("cold", "mid", "warm"), include.lowest = TRUE)]
  out <- sapply(c("cold", "mid", "warm"), function(g) {
    z <- s[terc == g]; if (uniqueN(z$cntr) < 2) return(NA_real_)
    coef(feols(R$f_tr, z, cluster = ~cntr, notes = FALSE))["gdd_mj"] })
  cat(sprintf("  %-14s cold %+.3e | mid %+.3e | warm %+.3e  (mean gdd %.0f / %.0f / %.0f)\n", k,
              out[1], out[2], out[3],
              s[terc == "cold", mean(gdd_mj)], s[terc == "mid", mean(gdd_mj)], s[terc == "warm", mean(gdd_mj)]))
}

## ---- (3) uncapped GDD on the phenological window ----------------------------------------------
# 11.crop_weather_gdd_window.csv: gdd = pmax(tg-5, 0) with NO 28 C cap, summed over a window whose
# edges are set by cumulative GDD anchors, so it moves with temperature instead of sitting at
# Mar-Jul. Both of the remaining candidate explanations change at once here; if the sign flips, the
# next step is to separate them.
cat("\n=== benchmark (capped GDD, fixed Mar-Jul) vs uncapped GDD on the GDD-window ===\n")
gw <- fread(file.path(d, "11.crop_weather_gdd_window.csv"))
gw[, precip2 := precip^2]
f_gw <- as.formula(paste("ln_yield ~ gdd + heat + frost + precip + precip2 |", FE))
for (k in R$crop_set) {
  s  <- R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr)
  m1 <- feols(R$f_tr, s, cluster = ~cntr, notes = FALSE)
  z  <- merge(s[, .(NUTS_ID, cntr, year, ln_yield, t, t2)], gw, by = c("NUTS_ID", "year"))
  z  <- R$common_obs(z, list(f_gw), ~cntr)
  m2 <- feols(f_gw, z, cluster = ~cntr, notes = FALSE)
  cat(sprintf("\n  %-14s                gdd          heat         frost        precip           N\n", k))
  cat(sprintf("    benchmark   %+.3e   %+.3e   %+.3e   %+.3e   %6d\n",
              coef(m1)["gdd_mj"], coef(m1)["heat_mj"], coef(m1)["frost_mj"], coef(m1)["precip_mj"], nobs(m1)))
  cat(sprintf("    GDD-window  %+.3e   %+.3e   %+.3e   %+.3e   %6d   [t(gdd)=%+.2f]\n",
              coef(m2)["gdd"], coef(m2)["heat"], coef(m2)["frost"], coef(m2)["precip"], nobs(m2),
              coef(m2)["gdd"] / se(m2)["gdd"]))
}
cat("\nNOTE: the GDD-window column changes BOTH the cap and the window at once.\n")
