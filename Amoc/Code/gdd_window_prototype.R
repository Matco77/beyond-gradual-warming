# PROTOTYPE - data-driven thermal-time (GDD) growing window. NO agronomic thresholds needed.
#
# Motivation: crops develop on accumulated heat (thermal time), not calendar dates, so the
# sensitive period falls earlier in warm regions/years and later in cool ones. A fixed Mar-Jul
# window mistimes it across space and over time, and cannot follow a shifted-climate scenario.
#
# Method (self-calibrating, threshold-free): for each region, read off how much GDD has
# accumulated by the OPEN (end of Feb) and CLOSE (end of Jul) of the baseline Mar-Jul window,
# averaged over a baseline period -> region anchors (g_open, g_close). For EVERY year, the
# growing window is the set of months whose cumulative GDD lies in (g_open, g_close]. By
# construction it equals Mar-Jul in the baseline; in a warm year cumulative GDD crosses the
# anchors earlier, so the window shifts earlier - automatically, and it will shift under an
# AMOC-cooled scenario too. Monthly resolution here (uses file 7, no E-OBS re-read); a daily
# version is the production upgrade. Author: Marco Bova
suppressMessages(library(data.table))
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

mo <- fread(file.path(d, "7.eobs_nuts3_crop_weather_monthly.csv"))   # NUTS_ID, year, month, gdd/heat/frost/precip
setorder(mo, NUTS_ID, year, month)
mo[, cumgdd := cumsum(gdd), by = .(NUTS_ID, year)]                   # heat accumulated so far this year

BASE <- 1990:2010                                                   # calibration baseline
FIX  <- 3:7                                                         # window we anchor to (Mar-Jul)
# region anchors = cumulative GDD at the open (end of month before FIX) and close (last FIX month)
anch <- mo[year %in% BASE, .(g_open  = mean(cumgdd[month == min(FIX) - 1L]),
                             g_close = mean(cumgdd[month == max(FIX)])), by = NUTS_ID]

# thermal window per region-year: months whose end-of-month cumGDD sits in (g_open, g_close]
gw <- merge(mo, anch, by = "NUTS_ID")[cumgdd > g_open & cumgdd <= g_close,
       .(gdd = sum(gdd), heat = sum(heat), frost = sum(frost), precip = sum(precip),
         m_open = min(month), m_close = max(month), n_month = .N),
       by = .(NUTS_ID, year)]
fwrite(gw, file.path(d, "11.crop_weather_gdd_window.csv"))

## ---- demonstration ------------------------------------------------------------------------
cat("1) baseline validation - thermal window should reproduce Mar-Jul (open~3, close~7):\n")
print(gw[year %in% BASE, .(mean_open = round(mean(m_open), 2), mean_close = round(mean(m_close), 2),
                           mean_len = round(mean(n_month), 2))])

cat("\n2) warming shift - window opening month, cool 1990s vs warm 2015-2024 (lower = earlier):\n")
sh <- gw[year %in% c(1990:1999, 2015:2024)][, era := fifelse(year < 2000, "1990s", "2015-24")]
print(sh[, .(mean_open_month = round(mean(m_open), 2), mean_close_month = round(mean(m_close), 2)), by = era])

cat("\n3) secular trend - regression of window-open month on year (negative = shifting earlier):\n")
tr <- gw[, .(open = mean(m_open)), by = year]
print(round(coef(lm(open ~ year, tr)), 4))

cat("\n4) vs fixed Mar-Jul: how much does the thermal window's GDD differ (share of region-years\n",
    "   whose window is NOT exactly months 3-7):\n")
cat("   ", round(100 * gw[!(m_open == 3 & m_close == 7), .N] / nrow(gw), 1), "% of region-years\n")
cat("\nwrote 11.crop_weather_gdd_window.csv | region-years:", nrow(gw), "\n")
