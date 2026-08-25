# Delta-method scenario replay (SKETCH) - CROP side (the documented extension of scenario_replay.R).
# Perturb the daily E-OBS record by an AMOC delta, then recompute the NUTS3 crop-window indicators
# (GDD/Heat/Frost/Precip) with the IDENTICAL construction used in estimation
# (cropweather_eobs_nuts3.R / weather_indicators.R) -> the response-function regressor is built the
# same way in estimation and replay, so beta is transportable to the counterfactual climate.
#
# Temperature: ADDITIVE delta (delta_fields/delta_<MODEL>_<EXP>_{tas,tasmax,tasmin}_lastthird.nc).
# Precip:      MULTIPLICATIVE ratio R = hosing/piControl (anomaly_output .../pr_*ratio.nc, last-third
#              mean per cell, clamped to [0.1,10]) -> rr_scen = rr * R. A ratio keeps precip >= 0 and
#              scales with local climate; an additive precip delta would go negative in dry cells.
#
# Runs AUTO-DISCOVERED from delta_fields/ (IPSL-CM6A-LR u03 and any new model picked up, no edit).
# Replayed over the CROP estimation window (YRS_CROP, weather_indicators.R) so scenario and
# historical indicators share the sample the beta were estimated on, and over BOTH estimation
# windows (Mar-Jul main, Apr-Aug robustness - the _mj / _aa regressors of regression.R).
# Remaining TODO: the delta_fields/ delta and ratio are ANNUAL, applied to every month.
# Crop indicators use DAILY values directly (heat<-tx, frost<-tn: the diurnal range), matching the
# estimation script (no within-day integration on the crop side).
# Author: Marco Bova
suppressMessages({library(terra); library(ncdf4); library(Matrix); library(data.table); library(fields)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/weather_indicators.R"))
ncf <- function(v) file.path(d, sprintf("EOBS/%s_ens_mean_0.25deg_reg_v31.0e.nc", v))

YRS <- YRS_CROP                                # = crop estimation panel (weather_indicators.R)
WINDOWS <- list("Mar-Jul" = 3:7, "Apr-Aug" = 4:8)   # main window (slide 14) + the robustness window
delta_dir <- file.path(d, "delta_fields")
RUNS <- lapply(list.files(delta_dir, "^delta_.*_tas_lastthird\\.nc$"), function(f) {
  s <- sub("^delta_(.*)_tas_lastthird\\.nc$", "\\1", f)
  list(model = sub("_[^_]*$", "", s), exp = sub(".*_", "", s))
})
cat("runs:", paste(sapply(RUNS, function(r) paste(r$model, r$exp)), collapse = " | "), "\n\n")

## ---- model-INDEPENDENT setup (once): NUTS3 crosswalk W + ncdf cell mapping ------------------
tmpl <- rast(ncf("tg"))[[1]]
xwalk <- fread(file.path(d, "6.eobs_to_nuts_crosswalk.csv"))
cells    <- sort(unique(xwalk$eobs_cell_id)); col <- match(xwalk$eobs_cell_id, cells)
nuts_ids <- sort(unique(xwalk$NUTS_ID));      row <- match(xwalk$NUTS_ID, nuts_ids)
W <- sparseMatrix(i = row, j = col, x = xwalk$nuts_weight, dims = c(length(nuts_ids), length(cells)))
xy  <- xyFromCell(tmpl, cells)
nc0 <- nc_open(ncf("tg")); nc_lon <- ncvar_get(nc0,"longitude"); nc_lat <- ncvar_get(nc0,"latitude")
tvec <- as.Date(ncvar_get(nc0,"time"), origin="1950-01-01"); nc_close(nc0); nlon <- length(nc_lon)
loni <- match(round(xy[,1],3), round(nc_lon,3)); lati <- match(round(xy[,2],3), round(nc_lat,3))
stopifnot(!anyNA(loni), !anyNA(lati))
ncdf_row <- (lati-1)*nlon + loni; yr <- as.integer(format(tvec,"%Y")); mo <- as.integer(format(tvec,"%m"))

# same area-weight-per-day-then-sum-to-month as cropweather_eobs_nuts3.R (renormalise over present cells)
ind_month <- function(mat, mm, fun) {
  v <- fun(mat); pres <- !is.na(v); v[!pres] <- 0
  reg <- as.matrix(W %*% v) / as.matrix(W %*% (pres * 1)); reg[!is.finite(reg)] <- NA
  t(rowsum(t(reg), mm))                                     # NUTS x month
}

## ---- per-cell ANNUAL delta / ratio on the E-OBS grid (rows aligned to `cells`) --------------
month_delta <- function(model, exp, var) {                  # additive, mm/day-irrelevant (temps in K)
  if (exp == "") return(matrix(0, length(cells), 12))       # self-check path: no perturbation
  ann <- values(resample(rast(file.path(delta_dir, sprintf("delta_%s_%s_%s_lastthird.nc", model, exp, var))),
                          tmpl, method = "bilinear"))[, 1]
  matrix(ann[cells], length(cells), 12)
}
ratio_path <- function(model, exp) {
  base <- file.path(d, "anomaly_output")
  if (model == "EC-Earth3")        file.path(base, "ECHearth3_anomaly", sprintf("pr_Amon_EC-Earth3_hos-%s-hos_ratio.nc", exp))
  else if (model == "IPSL-CM6A-LR") file.path(base, "IPSL_anomaly",     sprintf("pr_Amon_IPSL-CM6A-LR_hos-%s-hos_ratio.nc", exp))
  else file.path(base, if (grepl("LL", model)) "LL_anomaly" else "MM_anomaly",
                 sprintf("pr_ratio_%s-hos_over_piControl_1850-1949.nc", exp))
}
month_ratio <- function(model, exp) {                       # multiplicative precip factor at E-OBS cells
  if (exp == "") return(matrix(1, length(cells), 12))       # self-check path: no perturbation
  # Read via ncdf4 (terra rejects the native model grid: "lat not regularly spaced"). The ratio is
  # DIMENSIONLESS -> do NOT apply read_europe_cube's pr x86400.
  nc <- nc_open(ratio_path(model, exp))
  lon <- as.numeric(ncvar_get(nc, "lon")); lat <- as.numeric(ncvar_get(nc, "lat"))
  R3 <- ncvar_get(nc, "pr"); nc_close(nc)                   # [lon, lat, time], hosing/piControl ratio
  nt <- dim(R3)[3]; Rm <- apply(R3[, , (nt - nt %/% 3 + 1):nt, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  lon <- ((lon + 180) %% 360) - 180; ox <- order(lon); lon <- lon[ox]; Rm <- Rm[ox, ]   # interp.surface needs increasing x,y
  if (lat[1] > lat[length(lat)]) { oy <- order(lat); lat <- lat[oy]; Rm <- Rm[, oy] }
  Rm[!is.finite(Rm)] <- 1; Rm <- pmin(pmax(Rm, 0.1), 10)    # clamp (ymondiv spikes in dry cells)
  v <- interp.surface(list(x = lon, y = lat, z = Rm), xy)   # bilinear to E-OBS cell centres (rows = cells)
  v[!is.finite(v)] <- 1                                     # cells outside the model domain: no precip change
  matrix(v, length(cells), 12)                              # annual ratio for every month
}

## ---- one run: perturb daily E-OBS + recompute NUTS3 window indicators, compare to historical
VARS <- c("gdd","heat","frost","precip")
hist <- fread(file.path(d, "8.eobs_nuts3_crop_weather_window.csv"))[year %in% YRS]
run_one <- function(MODEL, EXP) {
  dtg <- month_delta(MODEL, EXP, "tas"); dtx <- month_delta(MODEL, EXP, "tasmax"); dtn <- month_delta(MODEL, EXP, "tasmin")
  Rr  <- month_ratio(MODEL, EXP)
  ncs <- list(tg=nc_open(ncf("tg")), tx=nc_open(ncf("tx")), tn=nc_open(ncf("tn")), rr=nc_open(ncf("rr")))
  rd <- function(nc,v,ti){a<-ncvar_get(nc,v,start=c(1,1,ti[1]),count=c(-1,-1,length(ti)));dim(a)<-c(nlon*length(nc_lat),length(ti));a[ncdf_row,,drop=FALSE]}
  out <- vector("list", length(YRS))
  for (k in seq_along(YRS)) { ti <- which(yr==YRS[k]); mm <- mo[ti]
    TG<-rd(ncs$tg,"tg",ti)+dtg[,mm]; TX<-rd(ncs$tx,"tx",ti)+dtx[,mm]
    TN<-rd(ncs$tn,"tn",ti)+dtn[,mm]; RR<-rd(ncs$rr,"rr",ti)*Rr[,mm]        # temp additive, precip multiplicative
    # NUTS3 x month first, ROUNDED to 3dp, then summed over the window - the exact two-step of
    # cropweather_eobs_nuts3.R (:140 rounds the monthly file, :152 sums the window), so a zero
    # delta lands on the historical value bit for bit rather than 1e-13 away from it.
    Mo <- lapply(list(gdd=list(TG,gdd_daily), heat=list(TX,heat_daily),
                      frost=list(TN,frost_daily), precip=list(RR,identity)),
                 function(a) round(ind_month(a[[1]], mm, a[[2]]), 3))
    out[[k]] <- rbindlist(lapply(names(WINDOWS), function(lab) {
      keep <- as.integer(colnames(Mo$gdd)) %in% WINDOWS[[lab]]
      ws   <- function(A) round(rowSums(A[, keep, drop = FALSE]), 3)      # NA month -> NA season
      data.table(NUTS_ID=nuts_ids, year=YRS[k], window=lab,
                 gdd=ws(Mo$gdd), heat=ws(Mo$heat), frost=ws(Mo$frost), precip=ws(Mo$precip)) }))
  }
  for (nc in ncs) nc_close(nc)
  cmp <- merge(rbindlist(out), hist, by=c("NUTS_ID","year","window"),
               suffixes=c("_scen","_hist"), all=TRUE)
  cmp[, `:=`(model=MODEL, exp=EXP)]                                       # after the merge: never NA
  for (lab in names(WINDOWS)) {
    z <- cmp[window == lab]
    cat(sprintf("=== CROP %s %s vs historical (%d-%d), NUTS3-year %s means ===\n",
                MODEL, EXP, min(YRS), max(YRS), lab))
    for (v in VARS) { a <- mean(z[[paste0(v,"_scen")]], na.rm=TRUE); b <- mean(z[[paste0(v,"_hist")]], na.rm=TRUE)
      cat(sprintf("  %-6s %8.2f -> %8.2f (%+.1f%%)\n", v, b, a, 100*(a/b-1))) }
    cat("\n")
  }
  cmp
}

## ---- self-check: zero delta / unit ratio MUST reproduce 8. exactly --------------------------
# With dtg=dtx=dtn=0 and Rr=1 the replay is the estimation pipeline with a "+0" and a "*1" in it,
# so any difference is a defect of the replay machinery (wrong cells, wrong month mapping, wrong
# window sum), not of the climate signal. Every scenario number downstream is a difference of
# these indicators, so a failure here silently biases all of them -> hard stop, no fallback.
selfcheck <- function() {
  cmp <- run_one("ZERO-DELTA", "")
  ok <- TRUE
  for (lab in names(WINDOWS)) { z <- cmp[window == lab]
    for (v in VARS) {
      a <- z[[paste0(v,"_scen")]]; b <- z[[paste0(v,"_hist")]]
      nmis <- sum(xor(is.na(a), is.na(b))); dmax <- max(abs(a-b), na.rm = TRUE)
      cat(sprintf("  %-8s %-6s max|scen-hist| %.3g | NA mismatch %d\n", lab, v, dmax, nmis))
      ok <- ok && nmis == 0L && dmax == 0 } }
  cat(sprintf("  rows compared %d | NUTS3 %d | years %d-%d\n",
              nrow(cmp), uniqueN(cmp$NUTS_ID), min(YRS), max(YRS)))
  if (!ok) stop("ZERO-DELTA SELF-CHECK FAILED: replay does not reproduce 8. exactly - do not use the scenarios")
  cat("  zero-delta self-check PASSED\n\n")
}
selfcheck()

allscen <- rbindlist(lapply(RUNS, function(r) {
  cmp <- run_one(r$model, r$exp)
  setnames(cmp[, c("model","exp","NUTS_ID","year","window", paste0(VARS,"_scen")), with=FALSE],
           paste0(VARS,"_scen"), VARS)
}))
fwrite(allscen, file.path(d, "scenario_crop_window.csv"))
cat("wrote scenario_crop_window.csv | rows", nrow(allscen), "| runs", uniqueN(allscen[,.(model,exp)]),
    "| windows", paste(unique(allscen$window), collapse=","), "\n")
