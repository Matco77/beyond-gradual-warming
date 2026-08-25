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
# Runs are AUTO-DISCOVERED from delta_fields/ (every model x protocol that has a tas delta),
# so IPSL-CM6A-LR (u03 only) and any new model are picked up with no code edit.
# Replayed over the ENERGY estimation window (YRS_ENERGY, weather_indicators.R) so the scenario and
# the historical indicators share the sample the beta were estimated on. Remaining TODO: the
# delta_fields/ deltas are ANNUAL and applied to every month (seasonal DJF/JJA fields exist but
# carry no CRS, so terra cannot grid them).
#
# Crop side: see scenario_replay_crop.R (NUTS3 crosswalk + crop indicators, precip multiplicative).
# Author: Marco Bova
suppressMessages({library(terra); library(ncdf4); library(Matrix); library(data.table)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

## ---- scenario spec ------------------------------------------------------------------------
YRS     <- YRS_ENERGY       # baseline years to perturb = energy estimation panel (weather_indicators.R)
YRS_EXT <- c(min(YRS) - 1L, YRS)   # + the preceding Oct-Dec, which the first Oct-Mar season needs
delta_dir <- file.path(d, "delta_fields")
# auto-discover (model, exp) from the tas delta files: delta_<MODEL>_<EXP>_tas_lastthird.nc
RUNS <- lapply(list.files(delta_dir, "^delta_.*_tas_lastthird\\.nc$"), function(f) {
  s <- sub("^delta_(.*)_tas_lastthird\\.nc$", "\\1", f)        # <MODEL>_<EXP>
  list(model = sub("_[^_]*$", "", s), exp = sub(".*_", "", s)) # model may carry hyphens, exp is the last _token
})
cat("runs discovered:", paste(sapply(RUNS, function(r) paste(r$model, r$exp)), collapse = " | "), "\n\n")

## ---- model-INDEPENDENT setup (once) -------------------------------------------------------
tmpl <- rast(ncf("tg"))[[1]]
# pop weights + ncdf cell mapping (same as energyweather)
wp <- fread(file.path(d, "12.eobs_population_weights.csv"))
cells <- sort(unique(wp$eobs_cell_id)); col <- match(wp$eobs_cell_id, cells)
cntrs <- sort(unique(wp$cntr));         row <- match(wp$cntr, cntrs)
Wp <- sparseMatrix(i = row, j = col, x = wp$w_pop, dims = c(length(cntrs), length(cells)))
xy <- xyFromCell(tmpl, cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0,"longitude"); nc_lat <- ncvar_get(nc0,"latitude")
tvec <- as.Date(ncvar_get(nc0,"time"), origin="1950-01-01"); nc_close(nc0); nlon <- length(nc_lon)
loni <- match(round(xy[,1],3), round(nc_lon,3)); lati <- match(round(xy[,2],3), round(nc_lat,3))
ncdf_row <- (lati-1)*nlon + loni; yr <- as.integer(format(tvec,"%Y")); mo <- as.integer(format(tvec,"%m"))
wmean <- function(v){p<-!is.na(v);v[!p]<-0; r<-as.matrix(Wp%*%v)/as.matrix(Wp%*%(p*1));r[!is.finite(r)]<-NA;r}
hist <- fread(file.path(d,"13.eobs_country_energy_weather_weighted.csv"))

## ---- per-cell ANNUAL delta -> E-OBS grid (sketch: annual field for every month) ------------
# EXP == "" -> zero delta (the self-check path: must reproduce historical exactly).
month_delta <- function(model, exp, var) {
  if (exp == "") return(matrix(0, ncell(tmpl), 12))
  ann <- values(resample(rast(file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, exp, var))),
                          tmpl, method = "bilinear"))[, 1]
  matrix(ann, ncell(tmpl), 12)
}

## ---- one run: perturb daily E-OBS + recompute country HDD/CDD, compare to historical -------
run_one <- function(MODEL, EXP) {
  DTG <- month_delta(MODEL, EXP, "tas"); DTX <- month_delta(MODEL, EXP, "tasmax"); DTN <- month_delta(MODEL, EXP, "tasmin")
  dtg <- DTG[cells,]; dtx <- DTX[cells,]; dtn <- DTN[cells,]     # rows align to `cells`
  ncs <- list(tg=nc_open(ncf("tg")), tx=nc_open(ncf("tx")), tn=nc_open(ncf("tn")))
  rd <- function(nc,v,ti){a<-ncvar_get(nc,v,start=c(1,1,ti[1]),count=c(-1,-1,length(ti)));dim(a)<-c(nlon*length(nc_lat),length(ti));a[ncdf_row,,drop=FALSE]}
  out <- vector("list", length(YRS_EXT))
  for (k in seq_along(YRS_EXT)) { ti <- which(yr==YRS_EXT[k]); mm <- mo[ti]
    TG<-rd(ncs$tg,"tg",ti)+dtg[,mm]; TX<-rd(ncs$tx,"tx",ti)+dtx[,mm]; TN<-rd(ncs$tn,"tn",ti)+dtn[,mm]  # delta method
    hdd<-wmean(within_day(TG,TX,TN,hdd_hinge)); cdd<-wmean(within_day(TG,TX,TN,cdd_hinge))
    out[[k]] <- data.table(cntr=rep(cntrs,length(ti)), date=rep(tvec[ti],each=length(cntrs)),
                           hdd=as.vector(hdd), cdd=as.vector(cdd)) }
  for (nc in ncs) nc_close(nc)
  dt <- rbindlist(out); dt[,`:=`(y=year(date),m=month(date))]
  # SAME three aggregates as energyweather_eobs_country.R:114-117, including the >=150-day rule that
  # drops an incomplete cross-year Oct-Mar season. YRS[1]-1 is read only to complete the first one;
  # the last year's Oct-Dec would label a season YRS[n]+1, which is outside YRS and dropped here
  # exactly as the 150-day rule drops it in the historical file.
  cal <- dt[y %in% YRS, .(hdd_calendar = sum(hdd)), by = .(cntr, year = y)]
  oct <- dt[m %in% c(10,11,12,1,2,3)][, .(hdd_octmar = sum(hdd), nd = .N),
            by = .(cntr, year = y + (m >= 10L))][nd >= 150 & year %in% YRS, .(cntr, year, hdd_octmar)]
  jja <- dt[m %in% 6:8 & y %in% YRS, .(cdd_jja = sum(cdd)), by = .(cntr, year = y)]
  scen <- Reduce(function(a,b) merge(a,b,by=c("cntr","year"),all=TRUE), list(cal,oct,jja))[order(cntr,year)]
  scen[, c("hdd_calendar","hdd_octmar","cdd_jja") := lapply(.SD, round, 2),
       .SDcols = c("hdd_calendar","hdd_octmar","cdd_jja")]        # same rounding as 13. -> exact comparability

  h <- hist[year %in% YRS, .(cntr,year,hdd_calendar,hdd_octmar,cdd_jja)]
  cmp <- merge(scen, h, by=c("cntr","year"), suffixes=c("_scen","_hist"), all=TRUE)
  cmp[, `:=`(model=MODEL, exp=EXP)]                               # after the merge: never NA on an outer join
  M <- function(x) mean(x, na.rm = TRUE)      # MT (Malta) is off-domain in E-OBS -> NA hdd, skip it
  cat("\n=== AMOC", MODEL, EXP, "vs historical (", min(YRS),"-",max(YRS),"), pop-weighted country-year means ===\n")
  for (v in c("hdd_calendar","hdd_octmar","cdd_jja")) {
    a <- M(cmp[[paste0(v,"_scen")]]); b <- M(cmp[[paste0(v,"_hist")]])
    cat(sprintf("  %-13s hist %6.0f -> scen %6.0f (%+.1f%%)\n", v, b, a, 100*(a/b-1)))
  }
  cmp
}

## ---- self-check: EXP=="" (zero delta) MUST reproduce 13. exactly ---------------------------
# With a zero delta the replay is the estimation pipeline with one extra "+ 0" in it, so any
# difference is a defect of the replay machinery (wrong cells, wrong month mapping, wrong window
# rule), not of the climate signal. Everything downstream is a difference of these indicators, so
# a failure here silently biases every scenario number -> hard stop, no fallback.
selfcheck <- function() {
  cmp <- run_one("ZERO-DELTA", "")
  ok <- TRUE
  for (v in c("hdd_calendar","hdd_octmar","cdd_jja")) {
    a <- cmp[[paste0(v,"_scen")]]; b <- cmp[[paste0(v,"_hist")]]
    nmis <- sum(xor(is.na(a), is.na(b))); dmax <- max(abs(a-b), na.rm = TRUE)
    cat(sprintf("  %-13s max|scen-hist| %.3g | NA mismatch %d\n", v, dmax, nmis))
    ok <- ok && nmis == 0L && dmax == 0
  }
  cat(sprintf("  rows compared %d | countries %d | years %d-%d\n",
              nrow(cmp), uniqueN(cmp$cntr), min(YRS), max(YRS)))
  if (!ok) stop("ZERO-DELTA SELF-CHECK FAILED: replay does not reproduce 13. exactly - do not use the scenarios")
  cat("  zero-delta self-check PASSED\n")
}
selfcheck()

allscen <- rbindlist(lapply(RUNS, function(r) {
  cmp <- run_one(r$model, r$exp)
  cmp[, .(model, exp, cntr, year, hdd_calendar = hdd_calendar_scen,
          hdd_octmar = hdd_octmar_scen, cdd_jja = cdd_jja_scen)]
}))
fwrite(allscen, file.path(d, "scenario_energy_country.csv"))
cat("\nwrote scenario_energy_country.csv | rows", nrow(allscen), "| runs", uniqueN(allscen[,.(model,exp)]), "\n")
