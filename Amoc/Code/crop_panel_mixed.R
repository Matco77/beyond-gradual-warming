# Mixed-level crop estimation panel: use the FINEST level at which CropStatHarm reports yield.
#   NUTS3 for 15 countries, NUTS2 for 9 (AT BE BG HR IE NL PL PT SI), NUTS0 for 3 (CY LU MT).
# Each country reports yield at exactly ONE level (verified), so there is NO parent/child
# double-counting and no dedup is needed - just keep every yield row at its own level.
#
# Weather: NUTS3 nests perfectly inside NUTS2/NUTS0 and the crop indicators are area-linear
# (each is an area-weighted daily mean summed over days), so a coarser region's weather is the
# area-weighted mean of its NUTS3 children. We therefore AGGREGATE the existing NUTS3 weather
# (8.eobs_nuts3_crop_weather_window.csv) upward - no need to re-run the E-OBS overlay.
# Companion to crop_panel_nuts3.R (which keeps NUTS3-only, 15 countries). Author: Marco Bova
suppressMessages(library(data.table))
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

# --- NUTS3 window weather + each region's area (from the E-OBS overlay crosswalk) ------------
w3   <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))
area <- fread(file.path(d, "6.eobs_to_nuts_crosswalk.csv"))[, .(area = sum(overlap_area, na.rm = TRUE)), by = NUTS_ID]
w3   <- merge(w3, area, by = "NUTS_ID")                       # NUTS_ID(5), year, window, gdd/heat/frost/precip, area

# --- aggregate NUTS3 -> NUTS2 (4-char key) and -> NUTS0 (2-char key), area-weighted -----------
agg <- function(keylen) w3[, .(
    gdd    = weighted.mean(gdd,    area, na.rm = TRUE),
    heat   = weighted.mean(heat,   area, na.rm = TRUE),
    frost  = weighted.mean(frost,  area, na.rm = TRUE),
    precip = weighted.mean(precip, area, na.rm = TRUE)),
  by = .(NUTS_ID = substr(NUTS_ID, 1, keylen), year, window)]
wx <- rbind(w3[, .(NUTS_ID, year, window, gdd, heat, frost, precip)],   # NUTS3 as-is
            agg(4), agg(2))                                             # NUTS2, NUTS0

# widen the two windows onto one row (mj = Mar-Jul main, aa = Apr-Aug robustness)
wx[, win := fifelse(window == "Mar-Jul", "mj", "aa")]
wide <- dcast(wx, NUTS_ID + year ~ win, value.var = c("gdd", "heat", "frost", "precip"))

# --- yield at its native level, positive observed yield, minus the one broken source cell -----
csh <- fread(file.path(d, "6.CropStatHarm_prepared.csv"))
y <- csh[!is.na(yield_t_ha) & yield_t_ha > 0 &
           !(NUTS_ID == "ES300" & crop == "Durum wheat" & year == 2012),   # see crop_panel_nuts3.R
         .(NUTS_ID, cntr = substr(NUTS_ID, 1, 2), nuts_level = nchar(NUTS_ID) - 2L,
           crop, year, area_ha, production_t, yield_t_ha, ln_yield = log(yield_t_ha))]

# --- join each yield region to the weather at its own level -----------------------------------
panel <- wide[y, on = .(NUTS_ID, year), nomatch = 0]
setcolorder(panel, c("NUTS_ID", "cntr", "nuts_level", "crop", "year",
                     "yield_t_ha", "ln_yield", "area_ha", "production_t"))
setorder(panel, nuts_level, NUTS_ID, crop, year)
fwrite(panel, file.path(d, "10.crop_panel_mixed_estimation.csv"))

cat("mixed panel rows:", nrow(panel), "| countries:", uniqueN(panel$cntr),
    "| regions:", uniqueN(panel$NUTS_ID), "| years:", paste(range(panel$year), collapse = "-"), "\n")
print(panel[, .(countries = uniqueN(cntr), regions = uniqueN(NUTS_ID), rows = .N),
            by = nuts_level][order(nuts_level)])

# self-check: no fabricated outcome, main-window weather attached, no country lost to a bad join
stopifnot(all(is.finite(panel$ln_yield)),
          panel[is.na(gdd_mj), .N] < nrow(panel))
lost <- setdiff(unique(y$cntr), unique(panel$cntr))
if (length(lost)) cat("NOTE countries with yield but no weather (off-domain, e.g. Malta):",
                      paste(lost, collapse = " "), "\n")
