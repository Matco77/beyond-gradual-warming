# Covariate mappings: showcase maps of every crop (NUTS3) and energy (country) covariate.
# Author: Marco Bova
#
# Reads the estimation panels and renders one value per region = mean over 2000-2019.
#   8.eobs_nuts3_crop_weather_window.csv   -> crop weather (GDD/Heat/Frost/Precip), Mar-Jul
#   9.crop_panel_nuts3_estimation.csv      -> crop yield by crop (t/ha)
#   16.energy_panel_estimation.csv         -> HDD/CDD, per-capita use, price, population
# Sequential viridis ramps (perceptually uniform, colourblind-safe), equal-area EPSG:3035.
# Writes four PNG figures to Amoc/figures/.
suppressMessages({library(sf); library(ggplot2); library(data.table); library(grid)})
sf_use_s2(FALSE)
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
out <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
YRS <- 1989:2024   # full crop-module period (yield 1989-2023, weather/energy through 2024)
EUR <- coord_sf(xlim = c(2.5e6, 6.1e6), ylim = c(1.4e6, 5.4e6), expand = FALSE)  # continental Europe, LAEA

# --- geometries (EPSG:3035) --------------------------------------------------
nuts <- st_read(file.path(d, "ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip"), quiet = TRUE)
n3 <- st_transform(nuts[nuts$LEVL_CODE == 3, "NUTS_ID"], 3035)
n0 <- st_transform(nuts[nuts$LEVL_CODE == 0, "NUTS_ID"], 3035)

# --- one choropleth ----------------------------------------------------------
mk <- function(sfd, fill, title, unit, opt, dir = 1, base = n0) {
  ggplot() +
    geom_sf(data = base, fill = "grey92", color = NA) +
    geom_sf(data = sfd, aes(fill = .data[[fill]]), color = "white", linewidth = 0.04) +
    scale_fill_viridis_c(option = opt, direction = dir, na.value = "grey85",
                         name = unit, guide = guide_colorbar(barwidth = 0.6, barheight = 4)) +
    EUR +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10.5, hjust = 0.01, margin = margin(b = 2)),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          plot.margin = margin(4, 4, 4, 4), plot.background = element_rect(fill = "white", color = NA))
}
# arrange independent-legend maps into one PNG (base grid, no extra package)
arrange_png <- function(plots, ncol, file, w, h) {
  png(file, width = w, height = h, res = 150); on.exit(dev.off())
  grid.newpage(); nr <- ceiling(length(plots) / ncol)
  pushViewport(viewport(layout = grid.layout(nr, ncol)))
  for (i in seq_along(plots)) print(plots[[i]], vp = viewport(
    layout.pos.row = ((i - 1) %/% ncol) + 1, layout.pos.col = ((i - 1) %% ncol) + 1))
}

## ============================ CROP =========================================
# weather: main window Mar-Jul, mean over YRS, all NUTS3
w <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[window == "Mar-Jul" & year %in% YRS]
wm <- w[, .(gdd = mean(gdd, na.rm = TRUE), heat = mean(heat, na.rm = TRUE),
            frost = mean(frost, na.rm = TRUE), precip = mean(precip, na.rm = TRUE)), by = NUTS_ID]
cw <- merge(n3, wm, by = "NUTS_ID")
p_gdd   <- mk(cw, "gdd",   "GDD  (Mar-Jul growing degree-days, base 5C)", "C.days", "D")
p_heat  <- mk(cw, "heat",  "Heat  (Mar-Jul degree-days over 28C)",       "C.days", "B")
p_frost <- mk(cw, "frost", "Frost  (Mar-Jul frost-days, tn<0)",          "days",   "G", dir = -1)
p_prec  <- mk(cw, "precip","Precipitation  (Mar-Jul total)",             "mm",     "G")
arrange_png(list(p_gdd, p_heat, p_frost, p_prec), 2,
            file.path(out, "map_crop_weather.png"), 1500, 1500)
cat("map_crop_weather.png\n")

# yield: mean t/ha over YRS by NUTS3 x crop (same unit -> shared scale, facet)
y <- fread(file.path(d, "9.crop_panel_nuts3_estimation.csv"))[year %in% YRS,
        .(yield = mean(yield_t_ha, na.rm = TRUE)), by = .(NUTS_ID, crop)]
cy <- merge(n3, y, by = "NUTS_ID")
gy <- ggplot() +
  geom_sf(data = n0, fill = "grey92", color = NA) +
  geom_sf(data = cy, aes(fill = yield), color = NA) +
  facet_wrap(~crop, ncol = 3) +
  # binned scale, breaks dense below 5 t/ha (where most regions sit) - more low-end contrast,
  # and the top open bin absorbs the few high outliers instead of stretching the whole ramp
  scale_fill_viridis_b(option = "D", breaks = c(1, 2, 3, 4, 5, 6, 8, 10),
                       na.value = "grey85", name = "t/ha") +
  EUR + labs(title = "Mean crop yield by crop, NUTS3  (1989-2023)") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, margin = margin(b = 4)),
        strip.text = element_text(face = "bold", size = 9), legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(6, 6, 6, 6))
ggsave(file.path(out, "map_crop_yield.png"), gy, width = 10, height = 7.5, dpi = 150)
cat("map_crop_yield.png\n")

## ============================ ENERGY =======================================
ep <- fread(file.path(d, "16.energy_panel_estimation.csv"))[year %in% YRS]
# climate + population are per country-year (same across fuels) -> dedup
clim <- unique(ep[, .(eurostat_geo, year, hdd_calendar, hdd_octmar, cdd_jja, population)])
cm <- clim[, .(hdd_calendar = mean(hdd_calendar, na.rm = TRUE), hdd_octmar = mean(hdd_octmar, na.rm = TRUE),
               cdd_jja = mean(cdd_jja, na.rm = TRUE), population = mean(population, na.rm = TRUE)),
           by = eurostat_geo]
ec <- merge(n0, cm, by.x = "NUTS_ID", by.y = "eurostat_geo")
p_hddc <- mk(ec, "hdd_calendar", "HDD  (calendar year, base 18/15C)",  "C.days", "G", dir = -1)
p_hddo <- mk(ec, "hdd_octmar",   "HDD  (Oct-Mar heating season)",      "C.days", "G", dir = -1)
p_cdd  <- mk(ec, "cdd_jja",      "CDD  (Jun-Aug cooling, base 21/24C)","C.days", "B")
arrange_png(list(p_hddc, p_hddo, p_cdd), 3,
            file.path(out, "map_energy_climate.png"), 2100, 850)
cat("map_energy_climate.png\n")

# per-capita household use by fuel (toe/person), price by fuel, population
fuse <- ep[fuel %in% c("Electricity", "Natural gas"),
           .(pc = mean(energy_ktoe * 1e3 / population, na.rm = TRUE)), by = .(eurostat_geo, fuel)]
fpr  <- ep[fuel %in% c("Electricity", "Natural gas") & !is.na(price),
           .(price = mean(price, na.rm = TRUE)), by = .(eurostat_geo, fuel)]
me <- function(dt, fl) merge(n0, dt[fuel == fl], by.x = "NUTS_ID", by.y = "eurostat_geo")
p_elu <- mk(me(fuse, "Electricity"), "pc",    "Electricity use per capita",   "toe/person", "D")
p_gau <- mk(me(fuse, "Natural gas"), "pc",    "Gas use per capita",           "toe/person", "C")
p_elp <- mk(me(fpr,  "Electricity"), "price", "Electricity price (incl. tax)","EUR/kWh",    "F")
p_gap <- mk(me(fpr,  "Natural gas"), "price", "Gas price (incl. tax)",        "EUR/GJ",     "F")
p_pop <- mk(ec, "population", "Population", "people", "D")
arrange_png(list(p_elu, p_gau, p_elp, p_gap, p_pop), 3,
            file.path(out, "map_energy_use_price_pop.png"), 2100, 1500)
cat("map_energy_use_price_pop.png\n")
cat("DONE ->", out, "\n")
