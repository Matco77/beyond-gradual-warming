# Delta-method scenario replay (SKETCH) - energy side.
# Perturb the daily E-OBS record by a climate delta (AMOC hosing or ISIMIP), then recompute the
# country HDD/CDD with the SAME within-day construction used in estimation (sources
# weather_indicators.R). This guarantees the response-function regressor is built identically in
# estimation and replay -> beta is transportable to the counterfactual climate.
#
# Why delta method: the AMOC signal (delta_fields/) is coarse (annual + DJF/JJA seasonal, ~2deg
# grid), so it can only supply a LEVEL shift per season, not sub-seasonal variability. The daily
# and within-day variability come from observed E-OBS; the model supplies only the change. Both
# scenarios (A=ISIMIP, B=AMOC) should run through this same module for a clean head-to-head.
#
# Crop side is the identical pattern on the NUTS3 crosswalk (6.eobs_to_nuts_crosswalk) with the
# crop indicators from weather_indicators.R; left as the documented extension.
# Author: Marco Bova
suppressMessages({library(terra); library(ncdf4); library(Matrix); library(data.table)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

## ---- scenario spec ------------------------------------------------------------------------
MODEL <- "HadGEM3-GC31-LL"; EXP <- "g01"    # AMOC hosing (g01 strong, u03 weaker). "" = no delta.
YRS   <- 2010:2019                          # baseline years to perturb (demo window; use full record for production)

## ---- 1. load AMOC deltas, regrid coarse model grid -> E-OBS grid ---------------------------
tmpl <- rast(ncf("tg"))[[1]]
# per-cell delta by calendar month. Sketch uses the ANNUAL delta for every month (single field,
# georeferences cleanly). SEASONAL refinement (DJF for Dec-Feb, JJA for Jun-Aug -> the winter-
# dominant AMOC cooling) needs delta_fields/seasonal/*.nc re-exported with a CRS: they currently
# carry no GDAL geotransform, so terra can't grid them. Plug the DJF/JJA fields in here once fixed.
month_delta <- function(var) {
  if (EXP == "") return(matrix(0, ncell(tmpl), 12))
  ann <- values(resample(rast(file.path(d, sprintf("delta_fields/delta_%s_%s_%s_lastthird.nc", MODEL, EXP, var))),
                          tmpl, method = "bilinear"))[, 1]
  matrix(ann, ncell(tmpl), 12)
}
DTG <- month_delta("tas"); DTX <- month_delta("tasmax"); DTN <- month_delta("tasmin")

## ---- 2. pop weights + ncdf cell mapping (same as energyweather) ----------------------------
wp <- fread(file.path(d, "12.eobs_population_weights.csv"))
cells <- sort(unique(wp$eobs_cell_id)); col <- match(wp$eobs_cell_id, cells)
cntrs <- sort(unique(wp$cntr));         row <- match(wp$cntr, cntrs)
Wp <- sparseMatrix(i = row, j = col, x = wp$w_pop, dims = c(length(cntrs), length(cells)))
xy <- xyFromCell(tmpl, cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0,"longitude"); nc_lat <- ncvar_get(nc0,"latitude")
tvec <- as.Date(ncvar_get(nc0,"time"), origin="1950-01-01"); nc_close(nc0); nlon <- length(nc_lon)
loni <- match(round(xy[,1],3), round(nc_lon,3)); lati <- match(round(xy[,2],3), round(nc_lat,3))
ncdf_row <- (lati-1)*nlon + loni; yr <- as.integer(format(tvec,"%Y")); mo <- as.integer(format(tvec,"%m"))
# delta per participating cell (rows align to `cells`)
dtg <- DTG[cells,]; dtx <- DTX[cells,]; dtn <- DTN[cells,]
wmean <- function(v){p<-!is.na(v);v[!p]<-0; r<-as.matrix(Wp%*%v)/as.matrix(Wp%*%(p*1));r[!is.finite(r)]<-NA;r}

## ---- 3. perturb daily E-OBS + recompute country HDD/CDD (identical construction) -----------
ncs <- list(tg=nc_open(ncf("tg")), tx=nc_open(ncf("tx")), tn=nc_open(ncf("tn")))
rd <- function(nc,v,ti){a<-ncvar_get(nc,v,start=c(1,1,ti[1]),count=c(-1,-1,length(ti)));dim(a)<-c(nlon*length(nc_lat),length(ti));a[ncdf_row,,drop=FALSE]}
out <- vector("list", length(YRS))
for (k in seq_along(YRS)) { ti <- which(yr==YRS[k]); mm <- mo[ti]
  TG<-rd(ncs$tg,"tg",ti)+dtg[,mm]; TX<-rd(ncs$tx,"tx",ti)+dtx[,mm]; TN<-rd(ncs$tn,"tn",ti)+dtn[,mm]  # delta method
  hdd<-wmean(within_day(TG,TX,TN,hdd_hinge)); cdd<-wmean(within_day(TG,TX,TN,cdd_hinge))
  out[[k]] <- data.table(cntr=rep(cntrs,length(ti)), date=rep(tvec[ti],each=length(cntrs)),
                         hdd=as.vector(hdd), cdd=as.vector(cdd)); cat(YRS[k],"done\n"); flush.console() }
for (nc in ncs) nc_close(nc)
dt <- rbindlist(out); dt[,`:=`(y=year(date),m=month(date))]
scen <- dt[, .(hdd_calendar=sum(hdd), cdd_jja=sum(cdd[m %in% 6:8])), by=.(cntr,year=y)]

## ---- 4. compare scenario vs historical (same years) ---------------------------------------
hist <- fread(file.path(d,"13.eobs_country_energy_weather_weighted.csv"))[year %in% YRS, .(cntr,year,hdd_calendar,cdd_jja)]
cmp <- merge(scen, hist, by=c("cntr","year"), suffixes=c("_scen","_hist"))
cat("\n=== AMOC", MODEL, EXP, "vs historical (", min(YRS),"-",max(YRS),"), pop-weighted country-year means ===\n")
cat("  HDD calendar: hist", round(mean(cmp$hdd_calendar_hist)), "-> scen", round(mean(cmp$hdd_calendar_scen)),
    "(", sprintf("%+.1f%%", 100*(mean(cmp$hdd_calendar_scen)/mean(cmp$hdd_calendar_hist)-1)), ")\n")
cat("  CDD JJA:      hist", round(mean(cmp$cdd_jja_hist)), "-> scen", round(mean(cmp$cdd_jja_scen)),
    "(", sprintf("%+.1f%%", 100*(mean(cmp$cdd_jja_scen)/mean(cmp$cdd_jja_hist)-1)), ")\n")
cat("  hottest-CDD countries, hist -> scen:\n")
print(cmp[year==max(YRS)][order(-cdd_jja_hist)][1:5, .(cntr, cdd_hist=round(cdd_jja_hist), cdd_scen=round(cdd_jja_scen))])

## ---- self-check: EXP=="" (zero delta) must reproduce historical exactly --------------------
# run this file with EXP<-"" to verify the replay machinery adds nothing when the delta is zero.
