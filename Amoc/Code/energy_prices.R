# 14.energy_prices_prepared.csv - household gas & electricity prices, annual, incl. taxes.
# Author: Marco Bova
#
# Sources (Eurostat, TOTAL price incl. all taxes I_TAX, EUR, standard medium household band;
# bi-annual S1/S2 averaged to annual):
#   Electricity: nrg_pc_204_h (<=2007, band Dc 4161150) spliced with nrg_pc_204 (2007+,
#                band DC KWH2500-4999). FULL 1985-2025.
#   Gas:         nrg_pc_202_h (<=2007, band D2 4141100) spliced with nrg_pc_202 (2007+,
#                band GJ20-199), unit GJ_GCV. FULL 1985-2025.
suppressMessages({library(data.table)})
d  <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
pr <- file.path(d, "PricesEectricitiesandgas")
cw <- fread(file.path(d, "0.country_crosswalk.csv"))

# generic bi-annual reader: keep one (band, unit, I_TAX, EUR) slice, average semesters -> annual
read_price <- function(file, dimnames, keep, fuel) {
  x <- fread(file.path(pr, file), sep = "\t", header = TRUE)
  setnames(x, 1, "dims"); x[, (dimnames) := tstrsplit(dims, ",")]
  for (k in names(keep)) x <- x[get(k) == keep[[k]]]
  scols <- grep("-S[12]", names(x), value = TRUE)
  m <- melt(x, id.vars = "geo", measure.vars = scols, variable.name = "period", value.name = "val")
  m[, val := as.numeric(gsub("[^0-9.-]", "", val))]
  m <- m[!is.na(val) & val > 0][, year := as.integer(substr(trimws(period), 1, 4))]
  m[, .(price = mean(val)), by = .(geo, year)][, fuel := fuel][]
}

old_dims <- c("freq","product","consom","unit","tax","currency","geo")
new_dims <- c("freq","siec","nrg_cons","unit","tax","currency","geo")

# electricity: historical (<2007) + plain biannual (>=2007), spliced at 2007
el_h <- read_price("estat_nrg_pc_204_h.tsv", old_dims,
                   list(consom="4161150", unit="KWH", tax="I_TAX", currency="EUR"), "Electricity")
el_n <- read_price("estat_nrg_pc_204.tsv", new_dims,
                   list(nrg_cons="KWH2500-4999", unit="KWH", tax="I_TAX", currency="EUR"), "Electricity")
el <- rbind(el_h[year < 2007], el_n[year >= 2007])

# gas: historical (<2007) + plain biannual (>=2007), spliced at 2007
ga_h <- read_price("estat_nrg_pc_202_h.tsv", old_dims,
                   list(consom="4141100", unit="GJ_GCV", tax="I_TAX", currency="EUR"), "Natural gas")
ga_n <- read_price("estat_nrg_pc_202.tsv", new_dims,
                   list(nrg_cons="GJ20-199", unit="GJ_GCV", tax="I_TAX", currency="EUR"), "Natural gas")
ga <- rbind(ga_h[year < 2007], ga_n[year >= 2007])

prices <- merge(rbindlist(list(el, ga)), cw[, .(geo = eurostat_geo, country_id)], by = "geo")
prices <- prices[, .(country_id, fuel, year, price, ln_price = log(price))][order(country_id, fuel, year)]
fwrite(prices, file.path(d, "14.energy_prices_prepared.csv"))
prices[, .(years = paste(range(year), collapse = "-"), n = .N, ncntry = uniqueN(country_id)), by = fuel]
