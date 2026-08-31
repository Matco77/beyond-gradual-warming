# ISIMIP branch, resolved into delta_Sv bins matched to the NAHosMIP bins (amoc_bin_fields.R), so
# the two branches can be plotted level-for-level instead of AMOC-bins-vs-one-static-ISIMIP-point.
# Author: Marco Bova
#
# WHY THIS EXISTS. isimip_delta_fields.R computes ONE climatological delta per model (mean of
# ssp126 2071-2100 minus mean of historical 1985-2014). That collapses the AMOC information in the
# ssp126 trajectory to nothing - but plot_amoc_sv_ssp126.R proves that trajectory is NOT a single
# level: it is a year-by-year AMOC-at-26N series (Terhaar msftyz for IPSL, own vo reconstruction for
# EC-Earth3), covering exactly the years the ISIMIP3b download holds ssp126 climate data for
# (2071-2100 decade files). So each ssp126 year has both an AMOC value and a climate day-record,
# same as a NAHosMIP hosing year - the only reason ISIMIP had been treated as unbinned was that nothing
# had gone and looked at the AMOC series behind it.
#
# BIN DEFINITION, deliberately mirrored on amoc_bin_fields.R's sv_bins():
#   delta_Sv of ssp126 year y = amoc_ssp126(y) - mean(amoc_historical), per model, floor()'d to 1 Sv
#   bins, kept with >= MIN_YEARS. NAHosMIP's delta_Sv is hos(y) - mean(piControl): an unforced,
#   stationary reference. ISIMIP has no piControl counterpart in this pairing, so the reference used
#   here is the mean of the model's OWN full historical AMOC series (1850-2014) - the closest
#   available analogue (both are "the model's own baseline circulation, not the perturbed one"), but
#   NOT the same physical quantity: the historical series is not unforced (it carries the actual
#   1850-2014 forcing history), so a bin's delta_Sv value is not on identically the same footing as
#   the NAHosMIP bin it is plotted against. Reported as data, not smoothed over.
#
# LEVEL-FOR-LEVEL MATCHING, per explicit instruction: an ISIMIP bin is written only if the SAME
# model's amoc_bin_fields_<MODEL>_u03.nc also kept that bin_lo. NAHosMIP's hosing bins run far
# beyond what 30 years of ssp126 ever reaches (bin_lo down to -10 Sv vs ISIMIP's -4/-5), so this
# intersection is what "stop once NAHosMIP exceeds ISIMIP" means in practice: bins beyond ISIMIP's
# reach are simply absent from this file, nothing is extrapolated to fill them.
#
# Only IPSL-CM6A-LR and EC-Earth3 have a ssp126 AMOC reconstruction (plot_amoc_sv_ssp126.R) and
# ISIMIP3b download (isimip_delta_fields.R) - same two-model restriction as the Phase 2 uncertainty
# band, for an unrelated reason each time.
#
# Output: isimip_bin_fields_<model>.nc, SAME schema as amoc_bin_fields_<MODEL>_u03.nc (dims lon, lat,
# month, bin; delta_tas/delta_tasmin/delta_tasmax K additive, pr_ratio dimensionless; bin_lo/bin_hi/
# n_years/delta_sv_mean) so scenario_replay_isimip_bins.R can reuse bin_scenarios() unchanged.
suppressMessages({library(ncdf4)})
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/isimip_common.R"))

VARS      <- c("tas", "tasmin", "tasmax", "pr")
MIN_YEARS <- 3L
HIST_WIN  <- 1985:2014          # fixed historical climatology, same window as isimip_delta_fields.R
SSP_YEARS <- 2071:2100          # the only ssp126 years ISIMIP3b climate data was downloaded for

# file-model-id -> display name used everywhere else (COL, MODELS in plot_impact_common.R)
DISPLAY <- c("ipsl-cm6a-lr" = "IPSL-CM6A-LR", "ec-earth3" = "EC-Earth3")

## ---- AMOC-at-26N annual series, same files/readers as plot_amoc_sv_ssp126.R -------------------
read_terhaar <- function(f) { nc <- nc_open(f)
  d <- list(t = as.numeric(ncvar_get(nc, "year")), a = as.numeric(ncvar_get(nc, "amoc"))); nc_close(nc); d }
read_ecearth <- function(f) { nc <- nc_open(f)
  d <- list(t = floor(as.numeric(ncvar_get(nc, "time"))), a = as.numeric(ncvar_get(nc, "amoc"))); nc_close(nc); d }

amoc_series <- function(model) {
  if (model == "ipsl-cm6a-lr") {
    hist <- read_terhaar(file.path(d, "terhaar_amoc/amoc/amoc/26.5N/historical/amoc_historical_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc"))
    ssp  <- read_terhaar(file.path(d, "terhaar_amoc/amoc/amoc/26.5N/ssp126/amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc"))
  } else {
    hist <- read_ecearth(file.path(d, "ecearth3_amoc26N_vo_below500m_historical_1850_2014.nc"))
    ssp  <- read_ecearth(file.path(d, "ecearth3_amoc26N_vo_below500m_ssp126_2015_2100.nc"))
  }
  list(hist = hist, ssp = ssp)
}

## ---- delta_Sv per ssp126 year (restricted to SSP_YEARS), floored to bins, intersected with the --
## NAHosMIP bin_lo set already kept for this model -------------------------------------------------
sv_bins <- function(model) {
  A <- amoc_series(model)
  base <- mean(A$hist$a, na.rm = TRUE)
  keep_t <- A$ssp$t %in% SSP_YEARS
  yr <- A$ssp$t[keep_t]; dsv <- A$ssp$a[keep_t] - base
  lo_all <- floor(dsv)

  disp <- DISPLAY[[model]]
  nc <- nc_open(file.path(d, sprintf("amoc_bin_fields_%s_u03.nc", disp)))
  nahos_lo <- as.integer(round(ncvar_get(nc, "bin_lo"))); nc_close(nc)

  tab <- table(lo_all)
  cand <- as.integer(names(tab))[tab >= MIN_YEARS]
  common <- sort(intersect(cand, nahos_lo))
  dropped <- setdiff(cand, common)
  if (length(dropped)) cat(sprintf("  %-14s dropped bin(s) %s: not in NAHosMIP's kept bin_lo for this model\n",
                                    disp, paste(dropped, collapse = ",")))
  bins <- lapply(common, function(b) list(lo = b, hi = b + 1L, years = yr[lo_all == b], dsv = mean(dsv[lo_all == b])))
  list(bins = bins, base = base, n_ssp_years = length(yr))
}

## ---- one model -------------------------------------------------------------------------------
build <- function(model) {
  disp <- DISPLAY[[model]]
  g <- isimip_grid(model); lon <- g$lon; lat <- g$lat; oy <- g$oy
  B <- sv_bins(model); bins <- B$bins; nb <- length(bins)
  cat(sprintf("%-14s ssp126 years available %d (%d-%d) | historical AMOC baseline %.3f Sv | bins kept %d\n",
              disp, B$n_ssp_years, min(SSP_YEARS), max(SSP_YEARS), B$base, nb))
  cat(sprintf("   bins: %s\n", paste(sapply(bins, function(b)
    sprintf("%d:%d(dSv=%.2f)", b$lo, length(b$years), b$dsv)), collapse = " ")))
  if (!nb) { cat("   no bin overlaps NAHosMIP's range for this model - nothing written\n\n"); return(invisible(NULL)) }

  ## fixed historical climatology (same for every bin, exactly as isimip_delta_fields.R computes it)
  hist_clim <- setNames(lapply(VARS, function(v) clim_monthly(path_for(model, "historical", v), HIST_WIN)[, oy, , drop = FALSE]), VARS)

  out <- list()
  for (v in VARS) {
    key <- if (v == "pr") "pr_ratio" else paste0("delta_", v)
    m <- array(NA_real_, c(length(lon), length(lat), 12, nb))
    fs <- path_for(model, "ssp126", v)
    for (bi in seq_len(nb)) {
      C <- clim_monthly(fs, bins[[bi]]$years)[, oy, , drop = FALSE]
      m[, , , bi] <- if (v == "pr") { r <- C / hist_clim[[v]]; r[!is.finite(r)] <- 1; r } else C - hist_clim[[v]]
    }
    out[[key]] <- m
    cat(sprintf("   %-7s done | non-finite %d\n", v, sum(!is.finite(m))))
  }

  dx <- ncdim_def("lon", "degrees_east", lon); dy <- ncdim_def("lat", "degrees_north", lat)
  dm <- ncdim_def("month", "1", 1:12); db <- ncdim_def("bin", "Sv", sapply(bins, function(b) b$lo + 0.5))
  un <- c(delta_tas = "K", delta_tasmin = "K", delta_tasmax = "K", pr_ratio = "1")
  vars <- lapply(names(out), function(k) ncvar_def(k, un[[k]], list(dx, dy, dm, db), NA, prec = "double"))
  meta <- list(ncvar_def("bin_lo", "Sv", db, NA, prec = "double"),
               ncvar_def("bin_hi", "Sv", db, NA, prec = "double"),
               ncvar_def("n_years", "1", db, NA, prec = "double"),
               ncvar_def("delta_sv_mean", "Sv", db, NA, prec = "double"))
  f <- file.path(d, sprintf("isimip_bin_fields_%s.nc", model))
  nc <- nc_create(f, c(vars, meta))
  for (k in names(out)) ncvar_put(nc, k, out[[k]])
  ncvar_put(nc, "bin_lo", sapply(bins, `[[`, "lo")); ncvar_put(nc, "bin_hi", sapply(bins, `[[`, "hi"))
  ncvar_put(nc, "n_years", sapply(bins, function(b) length(b$years)))
  ncvar_put(nc, "delta_sv_mean", sapply(bins, `[[`, "dsv"))
  ncatt_put(nc, 0, "model", disp)
  ncatt_put(nc, 0, "experiment", "ssp126 (ISIMIP3b w5e5, bias-adjusted)")
  ncatt_put(nc, 0, "delta_sv_definition", "AMOC(26N, ssp126 year) minus mean(AMOC(26N, historical 1850-2014)) - NOT a piControl baseline, see file header")
  ncatt_put(nc, 0, "binning", sprintf("floor to 1 Sv, bins with >= %d ssp126 years kept, INTERSECTED with NAHosMIP's kept bin_lo for this model", MIN_YEARS))
  ncatt_put(nc, 0, "historical_reference_window", paste(range(HIST_WIN), collapse = "-"))
  ncatt_put(nc, 0, "ssp126_years_available", paste(range(SSP_YEARS), collapse = "-"))
  ncatt_put(nc, 0, "pr_note", "ssp126-bin/historical ratio, bin mean per calendar month; [0.1,10] clamp applied by the consumer")
  ncatt_put(nc, 0, "grid", "ISIMIP3b native grid (already cropped to Europe at download), lat reordered increasing")
  ncatt_put(nc, 0, "created", format(Sys.time(), "%Y-%m-%d %H:%M"))
  nc_close(nc)
  cat("   wrote", basename(f), "\n\n")
  invisible(bins)
}

invisible(lapply(c("ipsl-cm6a-lr", "ec-earth3"), build))
cat("ISIMIP bin fields done.\n")
