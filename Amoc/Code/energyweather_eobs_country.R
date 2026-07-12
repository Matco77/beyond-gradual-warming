# Energy weather: population-weighted country HDD/CDD from E-OBS (proposal slides 18-19).
# Author: Marco Bova
#
# WHY COUNTRY, NOT NUTS3: the energy outcome (Eurostat nrg_bal_c household gas /
# electricity, 2.nrg_bal_c_prepared.csv) exists ONLY at country level - there is no
# subnational household-energy series. Unlike the crop module (CropStatHarm gave
# subnational yield), NUTS3 HDD/CDD would have no subnational outcome to regress on,
# and population-weighting to NUTS3 then up to country equals weighting grid->country
# directly. So the finest thing that matters is the 1km population grid used for the
# weights; the aggregation target is the country.
#
# Population weights (slide 18): w^pop_{g,c} = pop(cell g ∩ country c) / pop(country c).
# Daily indicators (slide 19): HDD_{g,d}=1(T<15)(18-T),  CDD_{g,d}=1(T>=24)(T-21), integrated over
# the WITHIN-DAY temperature distribution (diurnal sine on tg with amplitude (tx-tn)/2), not from
# the daily mean alone - the mean misses cooling on days that straddle 24C (34.5% of summer days).
# Country aggregates: HDD^Calendar (Jan-Dec), HDD^OctMar (Oct(t-1)-Mar t), CDD^JJA (Jun-Aug).
suppressMessages({library(terra); library(sf); library(ncdf4); library(Matrix); library(data.table)})
sf_use_s2(FALSE)
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

## --- Step 1: E-OBS population weights (1km Census-2021 pop, EPSG:3035) ---------
pop    <- rast(file.path(d, "GISCO_population_grid/ESTAT_OBS-VALUE-T_2021_V2.tiff"))
tmpl   <- rast(ncf("tg"))[[1]]
cellid <- init(tmpl, "cell")
cell_1km <- project(cellid, pop, method = "near")          # each 1km pop pixel -> its E-OBS cell
nuts0 <- st_read(file.path(d, "ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip"), quiet = TRUE)
nuts0 <- st_transform(nuts0[nuts0$LEVL_CODE == 0, "NUTS_ID"], 3035)
nuts0$idx <- seq_len(nrow(nuts0))
cntr_1km <- rasterize(vect(nuts0), pop, field = "idx")
z <- data.table(pop = values(pop)[,1], cell = values(cell_1km)[,1], idx = values(cntr_1km)[,1])
z <- z[!is.na(pop) & pop > 0 & !is.na(cell) & !is.na(idx)][, cntr := nuts0$NUTS_ID[idx]]

# snap population that landed on a data-less E-OBS cell (near-resample can pick a sea
# cell for coastal / tiny-island pixels, e.g. Malta) to the nearest reliable land cell.
# "reliable" = data on EVERY day of a recent fully-covered sample; this rejects flickering
# coastal sea cells that hold data only on rare days.
rr <- rast(ncf("tg"))
recent <- which(as.integer(format(terra::time(rr), "%Y")) >= 2015)
samp <- recent[round(seq(1, length(recent), length.out = 12))]
present <- app(rr[[samp]], function(x) sum(!is.na(x)))
valid_cells <- which(values(present)[, 1] == length(samp))
bad <- setdiff(unique(z$cell), valid_cells)
if (length(bad)) {
  # equal-area coords (EPSG:3035) so nearest = true ground distance, matching the crop crosswalk
  to3035 <- function(cells) st_coordinates(st_transform(
    st_as_sf(as.data.frame(xyFromCell(tmpl, cells)), coords = 1:2, crs = 4326), 3035))
  vx <- to3035(valid_cells); bx <- to3035(bad)
  nn <- integer(length(bad)); dist_km <- numeric(length(bad))
  for (i in seq_along(bad)) {
    d2 <- (vx[, 1] - bx[i, 1])^2 + (vx[, 2] - bx[i, 2])^2
    j <- which.min(d2); nn[i] <- valid_cells[j]; dist_km[i] <- sqrt(d2[j]) / 1000
  }
  # A data-less cell still holds REAL coastal/island population: snap it to the nearest cell
  # that has weather so those people stay in the country mean (the weight IS population -
  # dropping reachable people would bias the mean inland). But cap the reach: beyond max_snap_km
  # the pixel is off the E-OBS domain (Azores/Madeira sit ~1500 km from any land cell) and a
  # mainland proxy would be fabricated weather - drop it, as the crop module drops off-domain
  # regions. The cap sits in the empty gap between the near-continental islands kept
  # (Malta->Sicily ~110 km, farthest Greek/Italian islets ~220 km) and the Atlantic islands
  # dropped (>1500 km); their energy stays in the national outcome, only the weather is absent.
  max_snap_km <- 300
  off <- bad[dist_km > max_snap_km]
  cat("snap: kept", length(bad) - length(off), "cells (",
      round(z[cell %in% bad & !cell %in% off, sum(pop)] / 1e6, 3), "M ppl, max",
      round(max(dist_km[dist_km <= max_snap_km]), 0), "km) | off-domain dropped", length(off),
      "cells (", round(z[cell %in% off, sum(pop)] / 1e6, 3), "M ppl, up to",
      round(max(c(0, dist_km[dist_km > max_snap_km])), 0), "km)\n")
  z <- z[!cell %in% off]
  z[cell %in% bad, cell := nn[match(cell, bad)]]
}
wp <- z[, .(pop = sum(pop)), by = .(eobs_cell_id = cell, cntr)][, w_pop := pop / sum(pop), by = cntr][]
fwrite(wp[order(cntr, -w_pop)], file.path(d, "12.eobs_population_weights.csv"))
chk <- wp[, sum(w_pop), by = cntr]
cat("pop weights: countries", uniqueN(wp$cntr), "| max|sum(w)-1|", signif(max(abs(chk$V1 - 1)), 3),
    "| pop assigned", round(sum(wp$pop)/1e6, 1), "M\n")

## --- Step 2: population-weighted country HDD/CDD ------------------------------
cells <- sort(unique(wp$eobs_cell_id)); col <- match(wp$eobs_cell_id, cells)
cntrs <- sort(unique(wp$cntr));         row <- match(wp$cntr, cntrs)
Wp <- sparseMatrix(i = row, j = col, x = wp$w_pop, dims = c(length(cntrs), length(cells)))

xy <- xyFromCell(rast(ncf("tg"))[[1]], cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0, "longitude"); nc_lat <- ncvar_get(nc0, "latitude")
tvec <- as.Date(ncvar_get(nc0, "time"), origin = "1950-01-01"); nc_close(nc0)
nlon <- length(nc_lon)
loni <- match(round(xy[,1], 3), round(nc_lon, 3)); lati <- match(round(xy[,2], 3), round(nc_lat, 3))
stopifnot(!anyNA(loni), !anyNA(lati))
ncdf_row <- (lati - 1) * nlon + loni
yr <- as.integer(format(tvec, "%Y")); years <- sort(unique(yr))

ncs <- list(tg = nc_open(ncf("tg")), tx = nc_open(ncf("tx")), tn = nc_open(ncf("tn")))
wmean <- function(v) { p <- !is.na(v); v[!p] <- 0                  # pop-weighted country mean, NA-safe
  r <- as.matrix(Wp %*% v) / as.matrix(Wp %*% (p * 1)); r[!is.finite(r)] <- NA; r }
rd <- function(nc, v, ti) { a <- ncvar_get(nc, v, start = c(1,1,ti[1]), count = c(-1,-1,length(ti)))
  dim(a) <- c(nlon * length(nc_lat), length(ti)); a[ncdf_row, , drop = FALSE] }
daily <- vector("list", length(years))
for (k in seq_along(years)) {
  ti <- which(yr == years[k])
  TG <- rd(ncs$tg,"tg",ti); TX <- rd(ncs$tx,"tx",ti); TN <- rd(ncs$tn,"tn",ti)   # within-day needs tx/tn
  daily[[k]] <- data.table(cntr = rep(cntrs, length(ti)), date = rep(tvec[ti], each = length(cntrs)),
                           hdd = as.vector(wmean(within_day(TG, TX, TN, hdd_hinge))),
                           cdd = as.vector(wmean(within_day(TG, TX, TN, cdd_hinge))))
  cat(years[k], "done\n"); flush.console()
}
for (nc in ncs) nc_close(nc)
dt <- rbindlist(daily); dt[, `:=`(y = year(date), m = month(date))]

# daily population-weighted country HDD/CDD (kept for flexible re-windowing)
fwrite(dt[, .(cntr, date, hdd = round(hdd, 3), cdd = round(cdd, 3))],
       file.path(d, "13.eobs_country_hdd_cdd_daily.csv"))

cal <- dt[, .(hdd_calendar = sum(hdd)), by = .(cntr, year = y)]
oct <- dt[m %in% c(10,11,12,1,2,3)][, .(hdd_octmar = sum(hdd), nd = .N),
          by = .(cntr, year = y + (m >= 10L))][nd >= 150, .(cntr, year, hdd_octmar)]  # drop incomplete cross-year seasons
jja <- dt[m %in% 6:8, .(cdd_jja = sum(cdd)), by = .(cntr, year = y)]
res <- Reduce(function(a, b) merge(a, b, by = c("cntr","year"), all = TRUE), list(cal, oct, jja))[order(cntr, year)]
res[, c("hdd_calendar","hdd_octmar","cdd_jja") := lapply(.SD, round, 2),
    .SDcols = c("hdd_calendar","hdd_octmar","cdd_jja")]
fwrite(res, file.path(d, "13.eobs_country_energy_weather_weighted.csv"))
cat("country HDD/CDD: rows", nrow(res), "| countries", uniqueN(res$cntr),
    "| years", paste(range(res$year), collapse = "-"), "\n")

# NEXT (slide 20-21, not here): merge 13 with 2.nrg_bal_c (household gas/electricity),
# fuel prices, and population -> country x fuel x year energy estimation panel.
