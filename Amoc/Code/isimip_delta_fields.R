# Phase 3: monthly ISIMIP3b delta fields, on the same footing as the AMOC bin fields.
# Author: Marco Bova
#
# ISIMIP3b atmospheric variables (tas, tasmin, tasmax, pr) are archived DAILY, not monthly - the
# only ISIMIP3b "monthly" time step is ocean/marine forcing (chl, thetao, tos, ...), which is not
# what this branch needs (verified against the ISIMIP data API before downloading anything). The
# per-calendar-month climatology used here is therefore built by aggregating the daily record, not
# received pre-aggregated - but the number that comes out is the one a native monthly product would
# have given, since a mean of daily values IS the monthly mean.
#
# This mirrors amoc_bin_fields.R deliberately: same delta definition (additive for temperature,
# ratio for precipitation), same per-calendar-month resolution, same output shape (lon, lat, month),
# so scenario_replay_isimip.R can reuse bin_scenarios()'s interpolation code unchanged. The only
# structural difference is that ISIMIP has no ΔSv bins - one scenario per model, not forty - because
# ISIMIP is a single emissions trajectory (ssp126), not an ensemble indexed by AMOC state.
#
# Reference and future windows: historical 1985-2014, ssp126 2071-2100. Chosen to match the
# decade-file boundaries as closely as the ISIMIP archive allows (historical decades start in 1981,
# so 1981-1984 is read but excluded from the mean) and to give each window a full 30 years.
d <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/isimip_common.R"))

MODELS  <- c("ipsl-cm6a-lr", "ec-earth3")
VARS    <- c("tas", "tasmin", "tasmax", "pr")
WINDOWS <- list(historical = 1985:2014, ssp126 = 2071:2100)

## ---- one model -----------------------------------------------------------------------------
build <- function(model) {
  g <- isimip_grid(model); lon <- g$lon; lat <- g$lat; oy <- g$oy

  cat(sprintf("%-14s grid %dx%d | windows hist %d-%d (%d yr) ssp126 %d-%d (%d yr)\n", model,
              length(lon), length(lat), min(WINDOWS$historical), max(WINDOWS$historical), length(WINDOWS$historical),
              min(WINDOWS$ssp126), max(WINDOWS$ssp126), length(WINDOWS$ssp126)))

  clim <- list()
  for (scen in names(WINDOWS)) for (v in VARS) {
    fs <- path_for(model, scen, v)
    stopifnot(length(fs) > 0)
    C <- clim_monthly(fs, WINDOWS[[scen]])[, oy, , drop = FALSE]     # reorder to increasing lat
    clim[[paste(scen, v, sep = "_")]] <- C
    cat(sprintf("   %-10s %-7s done | files %d | non-finite %d\n", scen, v, length(fs), sum(!is.finite(C))))
  }

  ## delta: additive for temperature (K), ratio for precipitation (dimensionless, units cancel so
  ## kg m-2 s-1 needs no conversion). Clamp is NOT applied here - the engine applies [0.1, 10] after
  ## interpolation, exactly as for the AMOC bin fields, so both branches get identical treatment.
  delta_tas    <- clim$ssp126_tas    - clim$historical_tas
  delta_tasmin <- clim$ssp126_tasmin - clim$historical_tasmin
  delta_tasmax <- clim$ssp126_tasmax - clim$historical_tasmax
  pr_ratio     <- clim$ssp126_pr / clim$historical_pr
  pr_ratio[!is.finite(pr_ratio)] <- 1                # undefined where historical pr is 0: "no change"

  dx <- ncdim_def("lon", "degrees_east", lon); dy <- ncdim_def("lat", "degrees_north", lat)
  dm <- ncdim_def("month", "1", 1:12)
  vars <- list(ncvar_def("delta_tas", "K", list(dx, dy, dm), NA, prec = "double"),
               ncvar_def("delta_tasmin", "K", list(dx, dy, dm), NA, prec = "double"),
               ncvar_def("delta_tasmax", "K", list(dx, dy, dm), NA, prec = "double"),
               ncvar_def("pr_ratio", "1", list(dx, dy, dm), NA, prec = "double"))
  f <- file.path(d, sprintf("isimip_delta_fields_%s.nc", model))
  nc <- nc_create(f, vars)
  ncvar_put(nc, "delta_tas", delta_tas); ncvar_put(nc, "delta_tasmin", delta_tasmin)
  ncvar_put(nc, "delta_tasmax", delta_tasmax); ncvar_put(nc, "pr_ratio", pr_ratio)
  ncatt_put(nc, 0, "model", model)
  ncatt_put(nc, 0, "scenario", "ssp126")
  ncatt_put(nc, 0, "reference_window", paste(range(WINDOWS$historical), collapse = "-"))
  ncatt_put(nc, 0, "future_window", paste(range(WINDOWS$ssp126), collapse = "-"))
  ncatt_put(nc, 0, "delta_definition", "ssp126 monthly climatology minus historical monthly climatology (additive T, ratio pr)")
  ncatt_put(nc, 0, "note", "monthly climatology built by aggregating ISIMIP3b DAILY data - no native monthly atmospheric product exists")
  ncatt_put(nc, 0, "created", format(Sys.time(), "%Y-%m-%d %H:%M"))
  nc_close(nc)
  cat("   wrote", basename(f), "\n\n")
}

invisible(lapply(MODELS, build))
cat("ISIMIP delta fields done.\n")
