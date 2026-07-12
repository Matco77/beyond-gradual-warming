# NUTS3 crop estimation panel: CropStatHarm subnational yield  x  E-OBS NUTS3 weather.
# Unit = region (NUTS3 2016) x crop x year. Ready for the slide-16 response function:
#   ln(Yield_{r,k,t}) = b1k GDD + b2k Heat + b3k Frost + b4k Precip + b5k Precip^2
#                       + alpha_{r,k} + lambda_t + eps
# estimated at NUTS3 (within-country weather variation), not aggregated to country.
# Author: Marco Bova
suppressMessages({library(data.table)})

d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

# --- weather: NUTS3 x year, crop-window sums; widen so both windows sit on one row
w <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))
w[, win := fifelse(window == "Mar-Jul", "mj", "aa")]     # mj = main (Mar-Jul), aa = robustness (Apr-Aug)
wide <- dcast(w, NUTS_ID + year ~ win, value.var = c("gdd", "heat", "frost", "precip"))

# --- yield: CropStatHarm, NUTS3 only, positive observed yield (t/ha)
# Trusted as-is except ONE broken cell: ES300 durum 2012 reports 4 ha vs 14426 t -> 3606 t/ha
# (the region yields 1.9-5.2 t/ha every other year) - a source area error with 1336x leverage.
# Every other suspect cell is kept by choice.
csh <- fread(file.path(d, "6.CropStatHarm_prepared.csv"))
y <- csh[nuts_level == 3 & !is.na(yield_t_ha) & yield_t_ha > 0 &
           !(NUTS_ID == "ES300" & crop == "Durum wheat" & year == 2012),
         .(NUTS_ID, cntr = substr(NUTS_ID, 1, 2), crop, year,
           area_ha, production_t, yield_t_ha, ln_yield = log(yield_t_ha))]

# --- direct join on region + year (weather is crop-agnostic; crop response via b_k)
panel <- wide[y, on = .(NUTS_ID, year), nomatch = 0]
setcolorder(panel, c("NUTS_ID", "cntr", "crop", "year",
                     "yield_t_ha", "ln_yield", "area_ha", "production_t"))
setorder(panel, NUTS_ID, crop, year)

fwrite(panel, file.path(d, "9.crop_panel_nuts3_estimation.csv"))
cat("panel rows:", nrow(panel), "| regions:", uniqueN(panel$NUTS_ID),
    "| countries:", uniqueN(panel$cntr), "| years:", paste(range(panel$year), collapse = "-"),
    "| crops:", uniqueN(panel$crop), "\n")

# self-check: join integrity, no non-finite outcomes/regressors, weather attached
stopifnot(
  nrow(panel) == nrow(y[NUTS_ID %in% wide$NUTS_ID]),        # every in-domain yield row kept
  all(is.finite(panel$ln_yield)),
  panel[is.na(gdd_mj) | is.na(precip_mj), .N] == 0 |        # main-window weather present
    panel[is.na(gdd_mj), .N] < nrow(panel)                  # (early-year gaps may leave few NA)
)

# NOTE for the regression:
#  - Crops: {Total wheat, Total barley} are aggregates of {Soft,Durum} / {Spring,Winter}.
#    Use ONE set per regression, never totals + components together (double counting).
#  - FE: region-crop (NUTS_ID^crop) + year. Cluster SE by country (cntr) or region.
#  - Main weather = *_mj (Mar-Jul); robustness = *_aa (Apr-Aug). Precip^2 added in-model.
#  - Country-level exposure/valuation comes later by aggregating PREDICTED d ln(yield)
#    up to country with crop-area weights - estimation stays at NUTS3.
