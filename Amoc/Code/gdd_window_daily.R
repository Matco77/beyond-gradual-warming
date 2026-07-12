# Daily thermal-time (GDD) growing window - production version of gdd_window_prototype.R.
# Same self-calibrating idea (no agronomic thresholds), but the window is resolved to the DAY:
#   1. daily NUTS3 GDD/Heat/Frost/Precip from E-OBS (same overlay as cropweather_eobs_nuts3.R);
#   2. cumulative GDD within each year;
#   3. region anchors = mean baseline (1990-2010) cumGDD at the Mar-Jul window edges
#      (end of Feb = open, end of Jul = close);
#   4. window per region-year = the DAYS whose cumulative GDD lies in (g_open, g_close];
#   5. sum the daily indicators over those days.
# Warm years cross the anchors earlier -> window opens/closes earlier, smoothly. Writes
# 11.crop_weather_gdd_window.csv, then runs the crop response fixed-window vs GDD-window.
# Author: Marco Bova
suppressMessages({library(terra); library(sf); library(ncdf4); library(Matrix); library(data.table)})
sf_use_s2(FALSE)
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

## --- weight matrix W [n_nuts x n_cells] from the crop crosswalk (area weights) ---------------
xw <- fread(file.path(d, "6.eobs_to_nuts_crosswalk.csv"))
cells <- sort(unique(xw$eobs_cell_id)); col <- match(xw$eobs_cell_id, cells)
nuts_ids <- sort(unique(xw$NUTS_ID));   row <- match(xw$NUTS_ID, nuts_ids)
W <- sparseMatrix(i = row, j = col, x = xw$nuts_weight, dims = c(length(nuts_ids), length(cells)))

## --- E-OBS cell -> ncdf row mapping (same trick as cropweather) ------------------------------
r  <- rast(ncf("tg")); xy <- xyFromCell(r[[1]], cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0, "longitude"); nc_lat <- ncvar_get(nc0, "latitude")
tvec <- as.Date(ncvar_get(nc0, "time"), origin = "1950-01-01"); nc_close(nc0); nlon <- length(nc_lon)
loni <- match(round(xy[,1],3), round(nc_lon,3)); lati <- match(round(xy[,2],3), round(nc_lat,3))
stopifnot(!anyNA(loni), !anyNA(lati)); ncdf_row <- (lati-1)*nlon + loni
yr <- as.integer(format(tvec,"%Y")); mo <- as.integer(format(tvec,"%m"))

read_year <- function(nc, v, ti) { a <- ncvar_get(nc, v, start = c(1,1,ti[1]), count = c(-1,-1,length(ti)))
  dim(a) <- c(nlon*length(nc_lat), length(ti)); a[ncdf_row, , drop = FALSE] }
# area-weighted NUTS3 daily mean, NA-safe -> [n_nuts x n_days]
ind <- function(mat, fun) { v <- fun(mat); p <- !is.na(v); v[!p] <- 0
  m <- as.matrix(W %*% v) / as.matrix(W %*% (p*1)); m[!is.finite(m)] <- NA; m }
rowcumsum <- function(x) { x[is.na(x)] <- 0; t(apply(x, 1, cumsum)) }   # NA day contributes 0 heat

## --- PASS A: baseline anchors (cumGDD at end-of-Feb open and end-of-Jul close) ---------------
BASE <- 1990:2010
nc_tg <- nc_open(ncf("tg")); g_open <- g_close <- numeric(length(nuts_ids)); nb <- 0L
for (y in BASE) {
  ti <- which(yr == y); mm <- mo[ti]
  cg <- rowcumsum(ind(read_year(nc_tg, "tg", ti), function(x) pmax(x - 5, 0)))
  g_open  <- g_open  + cg[, max(which(mm == 2))]     # cumGDD by end of February
  g_close <- g_close + cg[, max(which(mm == 7))]     # cumGDD by end of July
  nb <- nb + 1L
}
nc_close(nc_tg); g_open <- g_open/nb; g_close <- g_close/nb
cat("anchors calibrated on", nb, "baseline years\n")

## --- PASS B: apply the thermal window to every estimation year, sum indicators ---------------
YRS <- 1989:2024
ncs <- list(tg = nc_open(ncf("tg")), tx = nc_open(ncf("tx")), tn = nc_open(ncf("tn")), rr = nc_open(ncf("rr")))
out <- vector("list", length(YRS))
for (k in seq_along(YRS)) {
  y <- YRS[k]; ti <- which(yr == y); doy <- as.integer(format(tvec[ti], "%j"))
  gdd  <- ind(read_year(ncs$tg, "tg", ti), function(x) pmax(x - 5,  0))
  heat <- ind(read_year(ncs$tx, "tx", ti), function(x) pmax(x - 28, 0))
  frost<- ind(read_year(ncs$tn, "tn", ti), function(x) (x < 0) * 1)
  prec <- ind(read_year(ncs$rr, "rr", ti), function(x) x)
  cg   <- rowcumsum(gdd)
  inwin <- (cg > g_open) & (cg <= g_close)            # [n_nuts x n_days] window membership
  sw <- function(m) rowSums(m * inwin, na.rm = TRUE)   # sum indicator over window days
  first_day <- apply(inwin, 1, function(z) if (any(z)) doy[which.max(z)] else NA_integer_)
  last_day  <- apply(inwin, 1, function(z) if (any(z)) doy[length(z) - which.max(rev(z)) + 1L] else NA_integer_)
  out[[k]] <- data.table(NUTS_ID = nuts_ids, year = y,
                         gdd = sw(gdd), heat = sw(heat), frost = sw(frost), precip = sw(prec),
                         open_doy = first_day, close_doy = last_day, n_day = rowSums(inwin))
  cat(y, "done\n"); flush.console()
}
for (nc in ncs) nc_close(nc)
res <- rbindlist(out)[n_day > 0]
res[, c("gdd","heat","frost","precip") := lapply(.SD, round, 3), .SDcols = c("gdd","heat","frost","precip")]
fwrite(res, file.path(d, "11.crop_weather_gdd_window.csv"))
cat("wrote 11.crop_weather_gdd_window.csv | region-years", nrow(res),
    "| mean window length", round(mean(res$n_day)), "days | mean open DOY", round(mean(res$open_doy, na.rm=TRUE)),
    "close DOY", round(mean(res$close_doy, na.rm=TRUE)), "\n")

## --- RESPONSE FUNCTION: fixed Mar-Jul vs GDD-window, side by side -----------------------------
suppressMessages(library(fixest))
cp <- fread(file.path(d, "9.crop_panel_nuts3_estimation.csv"))      # yield + Mar-Jul weather
gw <- merge(cp[, .(NUTS_ID, cntr, crop, year, ln_yield)], res, by = c("NUTS_ID", "year"))
cp[, precip_mj2 := precip_mj^2]; gw[, precip2 := precip^2]
crop_set <- c("Soft wheat", "Durum wheat", "Spring barley", "Winter barley")
for (kc in crop_set) {
  m_fix <- feols(ln_yield ~ gdd_mj + heat_mj + frost_mj + precip_mj + precip_mj2 | NUTS_ID + year,
                 cluster = ~cntr, data = cp[crop == kc])
  m_gdd <- feols(ln_yield ~ gdd + heat + frost + precip + precip2 | NUTS_ID + year,
                 cluster = ~cntr, data = gw[crop == kc])
  cat("\n===== CROP:", kc, "  (fixed Mar-Jul vs GDD-window) =====\n")
  print(etable(m_fix, m_gdd, headers = c("Mar-Jul", "GDD-window"), fitstat = ~ n + r2 + wr2))
}
