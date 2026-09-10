# ============================================================
# amoc_common.R  --  shared helpers for the AMOC-hosing R scripts.
#
# Sourced by: amoc_delta_effect_contrast, plot_delta_timeseries,
#             seasonal_fields_significance, plot_delta_fields_clipped.
#
# ONE definition of each shared piece so the four scripts cannot drift.
# In particular get_control_series() lives ONLY here and DESEASONALISES
# then de-drifts, so the effect ratio reported in the contrast CSV and in
# the time-series plot legends are now the SAME number (previously the
# contrast script skipped the deseasonalise step -> a different sd -> a
# different ratio for the identical run).
# ============================================================

library(ncdf4)

## ---- shared options ----
europe        <- list(lon = c(-15, 40), lat = c(34, 72))
MAX_PIC_FILES <- Inf      # all control years = most reliable noise estimate
EFFECT_K      <- 2        # |plateau| >= K * control sd  => an effect
CLIP_Q        <- 0.98     # map colour scale = this quantile of |delta|
nw_box  <- list(lon = c(-15,  5), lat = c(50, 62))   # cold-blob edge
med_box <- list(lon = c(  0, 25), lat = c(34, 45))   # Mediterranean
pic_root <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/CMIP6_piControl"

## ---- OneDrive-safe publish: build to a local temp, then move/copy into
## ---- place with a short retry (the CloudStorage File Provider intermittently
## ---- locks files mid-sync, so a direct write can fail "Permission denied"). ----
publish_file <- function(tmp, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  for (k in 1:5) {
    if (file.exists(dest)) try(Sys.chmod(dest, "0644"), silent = TRUE)
    if (suppressWarnings(file.rename(tmp, dest))) return(invisible(dest))                                    # same device
    if (suppressWarnings(file.copy(tmp, dest, overwrite = TRUE))) { unlink(tmp); return(invisible(dest)) }   # cross-device
    Sys.sleep(k)
  }
  warning("could not publish to ", dest, " (OneDrive lock?); kept temp at ", tmp, call. = FALSE)
  invisible(tmp)
}

## ---- piControl directory for an anomaly path ----
pic_dir_for <- function(p) file.path(pic_root, switch(basename(dirname(p)),
  "ECHearth3_anomaly" = "EC-Earth3", "LL_anomaly" = "HadGEM3-GC31-LL",
  "MM_anomaly" = "HadGEM3-GC31-MM", "IPSL_anomaly" = "IPSL-CM6A-LR",
  stop("unknown model folder: ", basename(dirname(p)))))

## ---- one nc -> Europe CUBE (lon, lat, time) on the native grid (no interp) ----
read_europe_cube <- function(path, region = europe) {
  nc  <- nc_open(path); vn <- strsplit(basename(path), "_")[[1]][1]
  lon <- as.numeric(ncvar_get(nc, "lon")); lat <- as.numeric(ncvar_get(nc, "lat"))
  v   <- ncvar_get(nc, vn); nc_close(nc)
  lon <- ((lon + 180) %% 360) - 180
  ox  <- order(lon); lon <- lon[ox]; v <- v[ox, , ]
  if (lat[1] > lat[length(lat)]) { oy <- order(lat); lat <- lat[oy]; v <- v[, oy, ] }
  ix  <- which(lon >= region$lon[1] & lon <= region$lon[2])
  iy  <- which(lat >= region$lat[1] & lat <= region$lat[2])
  # pr is stored in the anomaly files as kg m-2 s-1 (CMIP6 native); convert ONCE
  # here to mm/day (x 86400 s/day). Every downstream "mm/day" label -- plots,
  # CSVs, and the exported delta_*_lastthird.nc -- relies on this single point.
  if (vn == "pr") v <- v * 86400
  list(lon = lon[ix], lat = lat[iy], v = v[ix, iy, , drop = FALSE], vn = vn)
}

## ---- area-weighted (cos-lat) mean series over selected flattened cells ----
## sel indexes cells of the (lon x lat) plane flattened lon-fastest (matrix()).
masked_mean_series <- function(v, lat, sel) {
  nl <- dim(v)[1]; na <- dim(v)[2]; nt <- dim(v)[3]
  if (!length(sel)) return(rep(NA_real_, nt))
  M  <- matrix(v, nl * na, nt)[sel, , drop = FALSE]
  ws <- as.vector(matrix(cos(lat * pi / 180), nl, na, byrow = TRUE))[sel]
  apply(M, 2, function(col) { ok <- is.finite(col)
                              if (!any(ok)) NA_real_ else sum(col[ok] * ws[ok]) / sum(ws[ok]) })
}

## ---- one nc -> all-cell area-weighted Europe MONTHLY mean series ----
europe_monthly_means <- function(path, region = europe) {
  cb <- read_europe_cube(path, region)
  masked_mean_series(cb$v, cb$lat, seq_len(dim(cb$v)[1] * dim(cb$v)[2]))
}

## ---- control plateau-length-mean null (overlapping sliding windows) ----
mean_null <- function(series, w) {
  n <- length(series); if (n < w) return(numeric(0))
  vapply(1:(n - w + 1), function(s) mean(series[s:(s + w - 1)]), numeric(1))
}

## ---- cached control-noise series: DESEASONALISE, then DE-DRIFT. ----
## The anomaly (hosing - piControl monthly climatology) already has its
## seasonal cycle removed, so the control it is judged against must be
## deseasonalised too: otherwise window-means of length pl (pl is in months
## and rarely a multiple of 12) carry a seasonal variance the numerator does
## not have, inflating the control sd and shrinking the effect ratio. Defined
## once here so contrast CSV and plot legends report the identical ratio.
.ctrl_cache <- new.env()
get_control_series <- function(p) {
  vn  <- strsplit(basename(p), "_")[[1]][1]
  key <- paste(basename(dirname(p)), vn)
  if (!is.null(.ctrl_cache[[key]])) return(.ctrl_cache[[key]])
  pf  <- sort(list.files(pic_dir_for(p), pattern = paste0("^", vn, "_Amon_"), full.names = TRUE))
  pf  <- head(pf, MAX_PIC_FILES)
  s   <- unlist(lapply(pf, europe_monthly_means))                 # contiguous MONTHLY series
  mo  <- ((seq_along(s) - 1) %% 12) + 1
  s   <- s - as.numeric(tapply(s, mo, mean)[as.character(mo)])     # deseasonalise
  s   <- residuals(lm(s ~ seq_along(s))) + mean(s)                 # de-drift (spin-up)
  .ctrl_cache[[key]] <- s
  cat("  control months for", key, "=", length(s), "(", round(length(s) / 12), "yr, deseasonalised + de-drifted )\n")
  s
}

## ---- area-weighted mean of a 2D field (whole field, or a lon/lat box) ----
wmean <- function(f, lon, lat) {
  W <- matrix(cos(lat * pi / 180), length(lon), length(lat), byrow = TRUE)
  sum(f * W, na.rm = TRUE) / sum(W * !is.na(f))
}
box_mean <- function(fld, lon, lat, box) {
  ix <- which(lon >= box$lon[1] & lon <= box$lon[2])
  iy <- which(lat >= box$lat[1] & lat <= box$lat[2])
  if (!length(ix) || !length(iy)) return(NA_real_)
  wmean(fld[ix, iy, drop = FALSE], lon[ix], lat[iy])
}

## ---- country name per cell centre (NA = sea), cached per model grid ----
## Used by plot_delta_timeseries (country panels) and amoc_delta_effect_contrast
## (land-only plateau). Caller must guard with requireNamespace("maps").
.country_cache <- new.env()
cell_country <- function(model_key, lon, lat) {
  if (!is.null(.country_cache[[model_key]])) return(.country_cache[[model_key]])
  g  <- expand.grid(lon = lon, lat = lat)             # lon fastest -> matches matrix() flatten
  cc <- sub(":.*", "", maps::map.where("world", g$lon, g$lat))  # strip subregion suffix
  .country_cache[[model_key]] <- cc; cc
}

## ---- annual-mean piControl climatology FIELD (Europe crop), cached ----
## Used for the relative (%) precipitation page: % = 100 * delta / clim.
## Annual mean = mean of the monthly means (equal month weights); pr arrives
## in mm/day via read_europe_cube, matching the delta's units.
.pic_field_cache <- new.env()
pic_annual_field <- function(model, vn) {
  key <- paste(model, vn)
  if (!is.null(.pic_field_cache[[key]])) return(.pic_field_cache[[key]])
  pf <- sort(list.files(file.path(pic_root, model), pattern = paste0("^", vn, "_Amon_"), full.names = TRUE))
  stopifnot(length(pf) > 0)
  cat("  [clim field]", model, vn, "from", length(pf), "control file(s) ...\n")
  acc <- NULL; nt <- 0; cb <- NULL
  for (f in pf) {
    cb  <- read_europe_cube(f)
    s   <- apply(cb$v, c(1, 2), sum)
    acc <- if (is.null(acc)) s else acc + s
    nt  <- nt + dim(cb$v)[3]
  }
  .pic_field_cache[[key]] <- list(lon = cb$lon, lat = cb$lat, clim = acc / nt)
  .pic_field_cache[[key]]
}

## ---- plotting helpers shared by the two map scripts ----
panel_grid <- function(n) { nc <- min(3, max(1, ceiling(sqrt(n)))); c(ceiling(n / nc), nc) }
pal_for <- function(vn) {
  if (vn == "pr")
    colorRampPalette(c("#8C510A","#BF812D","#DFC27D","#F6E8C3","#F5F5F5",
                       "#C7EAE5","#80CDC1","#35978F","#01665E"))(64)
  else
    colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
                       "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(64)
}

## ---- self-check (run demo_amoc_common(); NOT called on source) ----
## Confirms the fix: deseasonalise+de-drift must collapse a strong seasonal
## cycle + drift down to the underlying noise, so the control sd (and thus the
## effect-ratio denominator) reflects weather noise, not the seasonal swing.
demo_amoc_common <- function() {
  set.seed(1)
  seas <- rep(sin(2 * pi * (1:12) / 12) * 5, 100)                 # seasonal amp ~5, 100 yr
  raw  <- seas + rnorm(1200, 0, 0.1) + 0.01 * (1:1200)            # + noise + linear drift
  mo   <- ((seq_along(raw) - 1) %% 12) + 1
  ds   <- raw - as.numeric(tapply(raw, mo, mean)[as.character(mo)])
  ds   <- residuals(lm(ds ~ seq_along(ds))) + mean(ds)
  stopifnot(sd(ds) < 0.5)                                         # seasonal + drift removed
  stopifnot(sd(mean_null(ds - mean(ds), 120)) < sd(mean_null(raw - mean(raw), 120)))
  cat("demo_amoc_common: OK\n")
}
