# Crop weather from E-OBS at NUTS3 (2016), monthly.
# Implements the crop-module weather pipeline of the thesis proposal, slides 13-14.
# Author: Marco Bova
#
# Slide 14 daily grid indicators (per E-OBS cell g, day d) - never average temp first:
#   GDD_{g,d}   = max(min(tg_{g,d}, 28) - 5, 0)   # moderate warmth, CAPPED at 28C: beyond 28
#                                                 # heat stops helping C3 cereals (photorespiration,
#                                                 # accelerated grain-fill) - the harmful part is
#                                                 # captured separately by Heat (Schlenker-Roberts).
#   Heat_{g,d}  = max(tx_{g,d} - 28, 0)
#   Frost_{g,d} = 1(tn_{g,d} < 0)
#   Precip_{g,d}= rr_{g,d}
# Slide 13 area-overlap crosswalk (EPSG:3035, slide 9):
#   omega_{g,r} = Area(g ∩ r) / sum_h Area(h ∩ r),   sum_g omega_{g,r} = 1
# Aggregation to NUTS3 x month (monthly = flexible superset of any crop window;
# sum the window's months later for March-July main / April-August robustness):
#   Indicator_{r,month} = sum_{d in month} sum_g omega_{g,r} Indicator_{g,d}
#
# Outputs:
#   6.eobs_to_nuts_crosswalk.csv         eobs_cell_id, NUTS_ID, CNTR_CODE, overlap_area, nuts_weight, quality_flag
#   7.eobs_nuts3_crop_weather_monthly.csv NUTS_ID, year, month, gdd, heat, frost, precip
suppressMessages({library(terra); library(sf); library(ncdf4); library(Matrix); library(dplyr); library(data.table)})
sf_use_s2(FALSE)
# daily indicator formulas live in ONE place, shared with the scenario replay: a response
# function is only transportable if the regressor is built identically in estimation and replay.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))

d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

## ------------------------------------------------------------------ ##
## Step A+B: E-OBS grid -> NUTS3 2016 area-overlap crosswalk (slide 13)
## ------------------------------------------------------------------ ##
# land mask = union of cells ever observed (E-OBS daily coverage drifts in time -
# e.g. Sicily is NA on 1950-01-01 but present later); per-day gaps handled in C+D.
r     <- rast(ncf("tg"))
samp  <- unique(round(seq(1, nlyr(r), length.out = 80)))
valid <- app(r[[samp]], fun = function(x) as.integer(any(!is.na(x))))
valid[valid == 0] <- NA
rid    <- mask(init(r[[1]], "cell"), valid)                 # cell index on union-valid cells
cellsf <- st_as_sf(as.polygons(rid, aggregate = FALSE))
names(cellsf)[1] <- "eobs_cell_id"
cellsf <- st_transform(cellsf, 3035)

nuts3 <- st_read(file.path(d, "ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip"), quiet = TRUE) |>
  filter(LEVL_CODE == 3) |>
  select(NUTS_ID, CNTR_CODE) |>
  st_transform(3035) |>
  st_make_valid()

inter <- st_intersection(nuts3, cellsf)
inter$overlap_area <- as.numeric(st_area(inter))            # m^2, EPSG:3035
xwalk <- inter |>
  st_drop_geometry() |>
  filter(overlap_area > 0) |>
  group_by(NUTS_ID) |>
  mutate(nuts_weight = overlap_area / sum(overlap_area)) |>
  ungroup() |>
  transmute(eobs_cell_id, NUTS_ID, CNTR_CODE, overlap_area, nuts_weight, quality_flag = "overlap")

# fallback for regions smaller than / between grid cells: nearest land cell within
# 40 km (~1.5 cells). Farther regions are outside the E-OBS domain (Canaries, Azores,
# French overseas, Malta) and stay unassigned by design.
miss <- setdiff(nuts3$NUTS_ID, xwalk$NUTS_ID)
if (length(miss)) {
  mr <- nuts3[nuts3$NUTS_ID %in% miss, ]
  cc <- st_centroid(cellsf)
  ni <- st_nearest_feature(st_centroid(st_geometry(mr)), cc)
  dist_km <- as.numeric(st_distance(st_centroid(st_geometry(mr)), cc[ni, ], by_element = TRUE)) / 1000
  fb <- tibble(eobs_cell_id = cellsf$eobs_cell_id[ni], NUTS_ID = mr$NUTS_ID,
               CNTR_CODE = mr$CNTR_CODE, overlap_area = NA_real_, nuts_weight = 1,
               quality_flag = "nearest_fallback", dist_km = dist_km) |>
    filter(dist_km <= 40) |> select(-dist_km)
  xwalk <- bind_rows(xwalk, fb)
}

chk <- xwalk |> group_by(NUTS_ID) |> summarise(s = sum(nuts_weight), .groups = "drop")
cat("crosswalk: NUTS3 covered", nrow(chk), "of", nrow(nuts3),
    "| fallback", sum(xwalk$quality_flag == "nearest_fallback"),
    "| max|sum(w)-1|", signif(max(abs(chk$s - 1)), 3),
    "| off-domain unassigned", length(setdiff(nuts3$NUTS_ID, xwalk$NUTS_ID)), "\n")
write.csv(xwalk, file.path(d, "6.eobs_to_nuts_crosswalk.csv"), row.names = FALSE)

## ------------------------------------------------------------------ ##
## Step C+D: daily grid indicators -> NUTS3-month weather (slide 14)
## ------------------------------------------------------------------ ##
# fixed sparse weight matrix W: [n_nuts x n_cells]
cells    <- sort(unique(xwalk$eobs_cell_id)); col <- match(xwalk$eobs_cell_id, cells)
nuts_ids <- sort(unique(xwalk$NUTS_ID));      row <- match(xwalk$NUTS_ID, nuts_ids)
W <- sparseMatrix(i = row, j = col, x = xwalk$nuts_weight,
                  dims = c(length(nuts_ids), length(cells)))

# map each participating terra cell -> row of the ncdf4 [lon,lat] grid (lon fastest).
# ncdf4 hyperslabs are ~400x faster than terra per-cell indexing on this file.
xy  <- xyFromCell(r[[1]], cells)
nc0 <- nc_open(ncf("tg"))
nc_lon <- ncvar_get(nc0, "longitude"); nc_lat <- ncvar_get(nc0, "latitude")
tvec   <- as.Date(ncvar_get(nc0, "time"), origin = "1950-01-01")   # E-OBS time origin
nc_close(nc0)
nlon <- length(nc_lon)
loni <- match(round(xy[, 1], 3), round(nc_lon, 3))
lati <- match(round(xy[, 2], 3), round(nc_lat, 3))
stopifnot(!anyNA(loni), !anyNA(lati))
ncdf_row <- (lati - 1) * nlon + loni

yr <- as.integer(format(tvec, "%Y")); mo <- as.integer(format(tvec, "%m"))
years <- sort(unique(yr))
ncs <- list(tg = nc_open(ncf("tg")), tx = nc_open(ncf("tx")),
            tn = nc_open(ncf("tn")), rr = nc_open(ncf("rr")))

read_year <- function(nc, v, ti) {
  a <- ncvar_get(nc, v, start = c(1, 1, ti[1]), count = c(-1, -1, length(ti)))
  dim(a) <- c(length(nc_lon) * length(nc_lat), length(ti))
  a[ncdf_row, , drop = FALSE]
}
# area-weighted regional mean per day (renormalised over cells with data that day),
# then summed over the month; NA month only if a day is fully unobserved region-wide.
ind_month <- function(mat, mm, fun) {
  v <- fun(mat); pres <- !is.na(v); v[!pres] <- 0
  reg <- as.matrix(W %*% v) / as.matrix(W %*% (pres * 1))
  reg[!is.finite(reg)] <- NA
  t(rowsum(t(reg), mm))
}

out <- vector("list", length(years))
for (k in seq_along(years)) {
  y <- years[k]; ti <- which(yr == y); mm <- mo[ti]
  gdd   <- ind_month(read_year(ncs$tg, "tg", ti), mm, gdd_daily)     # weather_indicators.R
  heat  <- ind_month(read_year(ncs$tx, "tx", ti), mm, heat_daily)
  frost <- ind_month(read_year(ncs$tn, "tn", ti), mm, frost_daily)
  prec  <- ind_month(read_year(ncs$rr, "rr", ti), mm, identity)
  months <- sort(unique(mm))
  out[[k]] <- data.table(
    NUTS_ID = rep(nuts_ids, times = length(months)),
    year = y, month = rep(months, each = length(nuts_ids)),
    gdd = as.vector(gdd), heat = as.vector(heat),
    frost = as.vector(frost), precip = as.vector(prec))
  cat(y, "done\n"); flush.console()
}
for (nc in ncs) nc_close(nc)

res <- rbindlist(out)[order(NUTS_ID, year, month)]
res[, c("gdd","heat","frost","precip") := lapply(.SD, round, 3), .SDcols = c("gdd","heat","frost","precip")]
fwrite(res, file.path(d, "7.eobs_nuts3_crop_weather_monthly.csv"))
cat("monthly weather: rows", nrow(res), "| regions", uniqueN(res$NUTS_ID),
    "| years", paste(range(res$year), collapse = "-"), "\n")

## ------------------------------------------------------------------ ##
## Step E: sum months -> crop window S_k (slide 14)
## ------------------------------------------------------------------ ##
# GDD_{r,k,t} = sum_{d in S_k} sum_g omega_{g,r} GDD_{g,d} = sum of the window's
# monthly values. Same window for wheat and barley (crop k enters later, at the
# country crop-area weighting). A season is NA if any of its months is NA.
# Main: March-July (3:7). Robustness: April-August (4:8).
win <- function(months, lab) res[month %in% months,
  .(window = lab, gdd = sum(gdd), heat = sum(heat), frost = sum(frost), precip = sum(precip)),
  by = .(NUTS_ID, year)]
resw <- rbindlist(list(win(3:7, "Mar-Jul"), win(4:8, "Apr-Aug")))[order(NUTS_ID, year, window)]
fwrite(resw, file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))
cat("window weather: rows", nrow(resw), "| windows", paste(unique(resw$window), collapse = ","),
    "| incomplete-season NA", resw[is.na(gdd), .N], "\n")

# NEXT (slide 15, not done here): crop-area-weight NUTS3 -> country using the
# fixed weights w^crop_{r,c,k} to get country-crop-year GDD/Heat/Frost/Precip.
