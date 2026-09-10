# Shared weather-indicator construction. Sourced by BOTH the historical weather scripts and the
# scenario replay so estimation and replay use the IDENTICAL formulas - the econometric
# transportability requirement (a response function is only valid to replay on a counterfactual
# climate if the regressor is built the same way in both). Do NOT duplicate these formulas.
# Author: Marco Bova

## --- within-day integration -----------------------------------------------------------------
# Integrate a daily hinge f over the diurnal cycle, approximated by a sine centred on the daily
# mean tg with amplitude (tx-tn)/2. This recovers threshold-straddling days: a summer day with
# mean < 24C but a hot afternoon still generates cooling demand, which f(daily mean) reads as 0.
# Matters where the threshold sits inside the daily swing (CDD@24 strongly; HDD@15 mildly; crop
# GDD@5/28 barely). Reduces EXACTLY to the daily-mean value when tx==tn. Matrices in, matrix out.
.PH24 <- 2 * pi * (0:23) / 24                         # 24 sub-daily samples (converges <0.02 dd)
within_day <- function(tg, tx, tn, hinge, ph = .PH24) {
  A <- (tx - tn) / 2
  s <- tg * 0
  for (h in ph) s <- s + hinge(tg + A * sin(h))
  s / length(ph)
}

## --- energy hinges (proposal slide 19) ------------------------------------------------------
# dead-band conventions: HDD counts (18-T) only below 15C; CDD counts (T-21) only at/above 24C.
hdd_hinge <- function(t) ifelse(t < 15,  18 - t, 0)
cdd_hinge <- function(t) ifelse(t >= 24, t - 21, 0)

## --- crop daily indicators (proposal slide 14) ----------------------------------------------
gdd_daily   <- function(tg) pmax(pmin(tg, 28) - 5, 0)  # moderate warmth, CAPPED at 28 (Schlenker-Roberts)
heat_daily  <- function(tx) pmax(tx - 28, 0)           # extreme heat above 28 (daily max)
frost_daily <- function(tn) (tn < 0) * 1               # frost day (daily min < 0)

## --- estimation sample windows --------------------------------------------------------------
# The replay must be run on the SAME years the beta were estimated on, otherwise the scenario-minus-
# historical difference mixes the climate perturbation with a change of sample. These are the year
# ranges of the two estimation panels (9.crop_panel_nuts3_estimation.csv,
# 16.energy_panel_estimation.csv); they live here because scenario_replay.R and
# scenario_replay_crop.R both need them and must not drift apart.
YRS_CROP   <- 1989:2023        # crop panel (CropStatHarm yield coverage)
YRS_ENERGY <- 1990:2024        # energy panel (Eurostat nrg_bal_c coverage)

## --- self-check -----------------------------------------------------------------------------
if (sys.nframe() == 0) {                               # runs only when executed directly
  m <- matrix(c(10, 22, 26), 1)                        # daily means; treat as tg
  # tx==tn -> within_day == daily-mean hinge
  stopifnot(all(abs(within_day(m, m, m, cdd_hinge) - cdd_hinge(m)) < 1e-9))
  # straddling day: mean 22 (<24) but max 30 -> daily-mean CDD 0, within-day CDD > 0
  stopifnot(cdd_hinge(matrix(22,1)) == 0,
            within_day(matrix(22,1), matrix(30,1), matrix(14,1), cdd_hinge) > 0)
  cat("weather_indicators.R self-check OK\n")
}
