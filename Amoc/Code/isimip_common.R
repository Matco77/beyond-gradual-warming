# Shared ISIMIP3b file access, used by isimip_delta_fields.R (one climatology per model) and
# isimip_bin_fields.R (one climatology per AMOC-matched delta_Sv bin). Extracted so the daily-file
# reading and per-calendar-month aggregation exist in exactly one place - the two callers must
# aggregate identically or their deltas are not comparable to each other.
# Author: Marco Bova
suppressMessages({library(ncdf4)})
ID <- file.path(d, "Isimip3b")   # `d` (datasets root) is set by the caller before sourcing this

path_for <- function(model, scen, var) {
  sort(Sys.glob(file.path(ID, sprintf("%s_r1i1p1f1_w5e5_%s_%s_lon*.nc", model, scen, var))))
}

## ---- per-calendar-month climatology over a set of files, restricted to a year window ----------
# Accumulates sum and count PER CELL (not just per month), matching the NA-safe style of
# amoc_bin_fields.R, though ISIMIP is not expected to carry missing values inside the Europe box.
clim_monthly <- function(files, years) {
  acc <- NULL; cnt <- NULL; nx <- ny <- NULL
  for (f in files) {
    nc  <- nc_open(f)
    v   <- names(nc$var)[1]
    tv  <- nc$dim$time$vals
    org <- sub(" .*", "", sub(".*since *", "", nc$dim$time$units))
    dt  <- as.Date(tv, origin = org)
    yr  <- as.integer(format(dt, "%Y")); mo <- as.integer(format(dt, "%m"))
    keep <- which(yr %in% years)
    if (!length(keep)) { nc_close(nc); next }                 # decade file entirely outside the window
    A <- ncvar_get(nc, v, start = c(1, 1, min(keep)), count = c(-1, -1, max(keep) - min(keep) + 1))
    nc_close(nc)
    keep_local <- keep - min(keep) + 1L                       # positions within the just-read slab
    if (is.null(acc)) { nx <- dim(A)[1]; ny <- dim(A)[2]; acc <- array(0, c(nx, ny, 12)); cnt <- array(0L, c(nx, ny, 12)) }
    for (m in 1:12) {
      idx <- keep_local[mo[keep] == m]
      if (!length(idx)) next
      S <- A[, , idx, drop = FALSE]
      fin <- is.finite(S); S[!fin] <- 0
      acc[, , m] <- acc[, , m] + rowSums(S, dims = 2)
      cnt[, , m] <- cnt[, , m] + rowSums(fin, dims = 2)
    }
  }
  out <- acc / cnt; out[cnt == 0] <- NA_real_
  out
}

## ---- ISIMIP grid, read once from any historical tas file (all vars/scenarios share the grid) ----
# lat is DECREASING in the ISIMIP files (75.75 -> 25.25); interp.surface (scenario_engine.R) needs
# increasing y. Returns lon as-is (already -180..180, increasing - no rewrap needed) and the
# reordering index `oy` for lat, which the caller applies to every field it reads.
isimip_grid <- function(model) {
  nc0 <- nc_open(path_for(model, "historical", "tas")[1])
  lon <- as.numeric(ncvar_get(nc0, "lon")); lat <- as.numeric(ncvar_get(nc0, "lat")); nc_close(nc0)
  stopifnot(all(diff(lon) > 0))
  oy <- order(lat)
  list(lon = lon, lat = lat[oy], oy = oy)
}
