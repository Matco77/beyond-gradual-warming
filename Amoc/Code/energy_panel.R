# 0.country_crosswalk.csv (slide 9) + energy estimation panel (slides 17,19-21).
# Merges Eurostat household energy (country names) with population-weighted HDD/CDD
# (2-letter codes) on one canonical country_id. Scope: EU27 + NO + UK + CH.
# Author: Marco Bova
suppressMessages({library(data.table)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")

# --- country crosswalk: country_id (ISO3) | country_name (Eurostat English) |
#     faostat_area (FAOSTAT spelling) | eurostat_geo (2-letter NUTS0 / E-OBS cntr)
cw <- fread(text = "
country_id,country_name,faostat_area,eurostat_geo
AUT,Austria,Austria,AT
BEL,Belgium,Belgium,BE
BGR,Bulgaria,Bulgaria,BG
HRV,Croatia,Croatia,HR
CYP,Cyprus,Cyprus,CY
CZE,Czechia,Czechia,CZ
DNK,Denmark,Denmark,DK
EST,Estonia,Estonia,EE
FIN,Finland,Finland,FI
FRA,France,France,FR
DEU,Germany,Germany,DE
GRC,Greece,Greece,EL
HUN,Hungary,Hungary,HU
IRL,Ireland,Ireland,IE
ITA,Italy,Italy,IT
LVA,Latvia,Latvia,LV
LTU,Lithuania,Lithuania,LT
LUX,Luxembourg,Luxembourg,LU
MLT,Malta,Malta,MT
NLD,Netherlands,Netherlands (Kingdom of the),NL
POL,Poland,Poland,PL
PRT,Portugal,Portugal,PT
ROU,Romania,Romania,RO
SVK,Slovakia,Slovakia,SK
SVN,Slovenia,Slovenia,SI
ESP,Spain,Spain,ES
SWE,Sweden,Sweden,SE
NOR,Norway,Norway,NO
GBR,United Kingdom,United Kingdom of Great Britain and Northern Ireland,UK
CHE,Switzerland,Switzerland,CH
")
fwrite(cw, file.path(d, "0.country_crosswalk.csv"))
cat("crosswalk: countries", nrow(cw), "\n")

# --- energy outcome (country name) -> country_id
en <- fread(file.path(d, "2.nrg_bal_c_prepared.csv"))          # Country, Year, Fuel, Energy_Use_Household_ktoe
en <- merge(en, cw[, .(Country = country_name, country_id, eurostat_geo)], by = "Country")

# --- HDD/CDD (2-letter code) -> country_id
w <- fread(file.path(d, "13.eobs_country_energy_weather_weighted.csv"))   # cntr, year, hdd_calendar, hdd_octmar, cdd_jja
w <- merge(w, cw[, .(cntr = eurostat_geo, country_id)], by = "cntr")

# --- annual population control (Eurostat DEMO_GIND) -> 15.population_country_year
pp <- fread(file.path(d, "totalpopulation.csv"))
pp[, geo := fifelse(geo == "Germany including former GDR", "Germany", geo)]  # -> crosswalk name
pop <- merge(pp[indic_de == "Population on 1 January - total"],
             cw[, .(geo = country_name, country_id)], by = "geo")
pop <- pop[, .(country_id, year = TIME_PERIOD, population = OBS_VALUE,
               ln_population = log(OBS_VALUE))][!is.na(population)]
fwrite(pop[order(country_id, year)], file.path(d, "15.population_country_year.csv"))

# --- panel: country_id x fuel x year, energy outcome + weather + population
panel <- merge(
  en[Energy_Use_Household_ktoe > 0,
     .(country_id, eurostat_geo, country_name = Country, fuel = Fuel, year = Year,
       energy_ktoe = Energy_Use_Household_ktoe, ln_energy = log(Energy_Use_Household_ktoe))],
  w[, .(country_id, year, hdd_octmar, hdd_calendar, cdd_jja)],
  by = c("country_id", "year"))
panel <- merge(panel, pop[, .(country_id, year, population, ln_population)],
               by = c("country_id", "year"), all.x = TRUE)
panel[, ln_energy_pc := ln_energy - ln_population]          # per-capita outcome (slide 17)

# fuel price control (built by energy_prices.R): household gas + electricity, annual incl. taxes,
# spliced 1990-2024 (partial country coverage - e.g. no NO/CH gas price).
pf <- file.path(d, "14.energy_prices_prepared.csv")
if (file.exists(pf)) {
  panel <- merge(panel, fread(pf)[, .(country_id, fuel, year, price, ln_price)],
                 by = c("country_id", "fuel", "year"), all.x = TRUE)
} else panel[, `:=`(price = NA_real_, ln_price = NA_real_)]

setcolorder(panel, c("country_id","eurostat_geo","country_name","fuel","year",
                     "energy_ktoe","ln_energy","population","ln_population","ln_energy_pc",
                     "price","ln_price","hdd_octmar","hdd_calendar","cdd_jja"))
setorder(panel, country_id, fuel, year)
fwrite(panel, file.path(d, "16.energy_panel_estimation.csv"))
cat("population: rows", nrow(pop), "| missing pop:", panel[is.na(population), .N],
    "| price coverage:", panel[!is.na(ln_price), .N], "of", nrow(panel),
    "rows (years", paste(range(panel[!is.na(ln_price), year]), collapse = "-"), ")\n")

cat("panel rows:", nrow(panel), "| countries:", uniqueN(panel$country_id),
    "| fuels:", paste(unique(panel$fuel), collapse = ","),
    "| years:", paste(range(panel$year), collapse = "-"), "\n")

# checks: energy scope countries all matched; report weather coverage
cat("energy scope countries matched:", uniqueN(en$country_id),
    "| in crosswalk but no energy (expect CH):",
    paste(setdiff(cw$country_id, en$country_id), collapse = ","), "\n")
stopifnot(!anyNA(panel$ln_energy))                       # every outcome present
na_w <- panel[is.na(hdd_calendar), unique(country_id)]
cat("rows missing weather:", panel[is.na(hdd_calendar), .N],
    "| only in:", paste(na_w, collapse = ","),
    "(Malta = Sicily proxy, sporadic E-OBS coastal gaps; regression drops these)\n")

# NOTE for the regression (slides 20-21):
#  Gas:   ln(energy) ~ hdd_octmar (or hdd_calendar) + ln(price) | country_id + year
#  Elec:  ln(energy) ~ hdd_calendar + cdd_jja       + ln(price) | country_id + year
#  fuel == "Natural gas" / "Electricity"; "Total" is a robustness outcome.
#  Controls attached: fuel prices (14.energy_prices_prepared, gas+elec 1990-2024, incl. taxes,
#  partial country coverage) and annual population (15.population_country_year -> ln_energy_pc).
