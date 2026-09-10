# Robustness BOUND: how much of the crop impact rests on beta_gdd alone.
# Author: Marco Bova
#
# NOT an alternative preferred specification. Dropping gdd_mj is mis-specification, not a fix:
# gdd_mj and heat_mj correlate +0.52..+0.57 WITHIN the benchmark's fixed effects
# (crop_gdd_diagnostics.R), so removing gdd does not delete the thermal signal, it reassigns it to
# heat/frost/precip - moving the contamination from a coefficient that is already flagged as
# uninterpretable into three that are currently readable. The mirror experiment is already on
# record: dropping heat_mj instead moves beta_gdd by -13% to +33%.
#
# What it IS for: a bound. The whole "AMOC weakening raises European yields" sign comes from
# beta_gdd < 0, which does not survive spec A -> C (soft wheat loses significance, spring barley
# changes sign, winter barley collapses to zero; the 28 C cap is inert, the window does the work).
# This quantifies the exposure - how far the crop number moves if that one coefficient is removed
# entirely - so the reader sees the magnitude of what the caveat is about instead of taking it on
# trust.
#
# SAMPLE. Both columns are anchored to f_tr's common_obs sample, i.e. the BENCHMARK's rows, not
# each spec's own. Otherwise the comparison would mix a coefficient change with a sample change -
# the same reason regression.R's build-up tables anchor to the most saturated spec.
invisible(capture.output(
  source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/amoc_impact.R"))))
cat("beta, weights and crop_effect() loaded from amoc_impact.R\n")

ISCEN <- Sys.getenv("ISIMIP_SCEN", "ssp126")
ISFX  <- if (ISCEN == "ssp126") "" else paste0("_", ISCEN)

NOGDD_X <- setdiff(CROP_X, "gdd_mj")
f_ng    <- as.formula(paste("ln_yield ~", paste(NOGDD_X, collapse = " + "), "| NUTS_ID[t, t2] + year"))
nogdd_fit <- lapply(setNames(R$crop_set, R$crop_set), function(k)
  feols(f_ng, R$common_obs(cp9[crop == k], list(R$f_tr), ~cntr), cluster = ~cntr, notes = FALSE))

## ---- omitted-variable bias on the coefficients that stay ---------------------------------------
cat("\n=== beta with and without gdd_mj (same rows) - what the other coefficients absorb ===\n")
for (k in R$crop_set) {
  b1 <- coef(crop_fit[[k]]); b2 <- coef(nogdd_fit[[k]])
  cat(sprintf("\n  %-14s %14s %14s %10s\n", k, "benchmark", "drop gdd", "change"))
  for (v in NOGDD_X)
    cat(sprintf("    %-12s %+13.3e %+14.3e %+9.1f%%\n", v, b1[[v]], b2[[v]], 100 * (b2[[v]] / b1[[v]] - 1)))
  cat(sprintf("    %-12s %+13.3e %14s\n", "gdd_mj", b1[["gdd_mj"]], "(dropped)"))
}

## ---- re-apply beta with the reduced spec -------------------------------------------------------
# ponytail: crop_effect() reads crop_fit/CROP_X from the enclosing env, so swapping the two globals
# reuses the whole beta-application and aggregation path unchanged rather than cloning it. Restored
# below. Upgrade path if a third spec ever needs this: make them arguments of crop_effect().
eu <- function(ce) merge(ce, w_crop[, .(w = sum(w)), by = .(cntr, crop)], by = c("cntr", "crop"))[
  , .(dln = weighted.mean(dln, w)), by = c(BINKEY, "component")]

bench_fit <- crop_fit; bench_X <- CROP_X
run <- function(file) {
  crop_fit <<- bench_fit; CROP_X <<- bench_X;  a <- eu(crop_effect(file))[, spec := "benchmark"]
  crop_fit <<- nogdd_fit; CROP_X <<- NOGDD_X;  b <- eu(crop_effect(file))[, spec := "drop_gdd"]
  crop_fit <<- bench_fit; CROP_X <<- bench_X
  rbindlist(list(a, b))
}
CE <- rbindlist(list(
  run("scenario_bins_crop_window.csv.gz")[, branch := "NAHosMIP (A)"],
  run(sprintf("scenario_isimip_bins_crop_window%s.csv.gz", ISFX))[, branch := sprintf("ISIMIP %s (T)", ISCEN)]
), fill = TRUE)

tot <- CE[, .(dln = sum(dln)), by = .(branch, spec, model, bin_id, delta_sv, n_years)]
tot[, pct := pct(dln)]
W <- dcast(tot, branch + model + bin_id + delta_sv + n_years ~ spec, value.var = "pct")
W[, `:=`(gdd_share = benchmark - drop_gdd)]      # percentage points the gdd term contributes
setorder(W, branch, model, -bin_id)
fwrite(CE, file.path(d, sprintf("crop_dropgdd_components%s.csv", ISFX)))
fwrite(W,  file.path(d, sprintf("crop_dropgdd_bound%s.csv", ISFX)))

stopifnot(nrow(W) > 0, !anyNA(W$benchmark), !anyNA(W$drop_gdd))

cat(sprintf("\n=== CROP total, Europe: benchmark vs drop-gdd (BOUND, not a preferred spec) ===\n"))
cat("gdd pp = percentage points the gdd_mj term contributes to the benchmark total\n")
for (br in unique(W$branch)) {
  cat(sprintf("\n %s\n", br))
  z <- W[branch == br]
  cat(sprintf("  %-14s %6s %4s %11s %11s %10s\n", "model", "bin", "n", "benchmark%", "drop gdd%", "gdd pp"))
  for (i in seq_len(nrow(z)))
    cat(sprintf("  %-14s %6.1f %4d %+10.2f%% %+10.2f%% %+9.2f\n", z$model[i], z$bin_id[i],
                z$n_years[i], z$benchmark[i], z$drop_gdd[i], z$gdd_share[i]))
}
cat("\n!! Read as EXPOSURE, not as a corrected estimate: the drop-gdd column is mis-specified\n")
cat("!! (gdd_mj and heat_mj correlate ~0.55 within the FE, so heat/frost/precip absorb the\n")
cat("!! thermal signal - see the beta table above). It bounds how much of the crop result is\n")
cat("!! carried by the one coefficient that spec A -> C shows is not interpretable.\n")
cat(sprintf("\nwrote crop_dropgdd_bound%s.csv and crop_dropgdd_components%s.csv\n", ISFX, ISFX))
