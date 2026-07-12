#!/usr/bin/env python3
# =====================================================================
# Reconstruct monthly tasmin and tasmax for EC-Earth3 NAHosMIP (hosing)
# from EC-Earth3 piControl  --  Annalisa's method.
#
#   Per cell, per CALENDAR month, learn in piControl:
#       tasmax = a_max + b_max * tas
#       tasmin = a_min + b_min * tas
#   then push the hosing tas through those lines.
#
#   slope b  -> how the day-night spread responds as the mean shifts
#   intercept a -> anchors the absolute level
#
# Years need NOT match: the fit is per calendar month (all Januaries
# together, etc.), so piControl 2259-2759 teaches hosing 1850-1949.
#
# Grounding: Weedon et al. (2010) WATCH Forcing Data; stationary
# mean<->extreme statistical link.  Everything MONTHLY, native gr cell.
# If hosing grid differs from piControl, the fitted COEFFICIENTS are
# regridded onto the hosing grid (nearest-neighbour) -- no data invented.
# =====================================================================

import os
import glob
import numpy as np
import xarray as xr

# ---------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------
PICONTROL_DIR = ("/Users/Bova/Library/CloudStorage/"
                 "OneDrive-UniversitàCommercialeLuigiBocconi/"
                 "1.Tesi/Amoc/datasets/CMIP6_piControl/EC-Earth3")

HOSING_DIR = ("/Users/Bova/Library/CloudStorage/"
              "OneDrive-UniversitàCommercialeLuigiBocconi/"
              "1.Tesi/Amoc/datasets/NAHosMIP/EC-Earth3")

OUTPUT_DIR = HOSING_DIR        # write reconstructed files beside the hosing tas

# piControl files that teach the regression (Amon = monthly)
PIC_TAS_GLOB    = os.path.join(PICONTROL_DIR, "tas_Amon_EC-Earth3_piControl_*_gr_*.nc")
PIC_TASMAX_GLOB = os.path.join(PICONTROL_DIR, "tasmax_Amon_EC-Earth3_piControl_*_gr_*.nc")
PIC_TASMIN_GLOB = os.path.join(PICONTROL_DIR, "tasmin_Amon_EC-Earth3_piControl_*_gr_*.nc")

# Hosing experiments to process: label -> glob for its monthly tas file(s)
HOSING_EXPERIMENTS = {
    "hos-g01-hos": os.path.join(HOSING_DIR, "tas_Amon_EC-Earth3_hos-g01-hos_*.nc"),
    "hos-u03-hos": os.path.join(HOSING_DIR, "tas_Amon_EC-Earth3_hos-u03-hos_*.nc"),
}

# ---------------------------------------------------------------------
# 2. HELPERS
# ---------------------------------------------------------------------
def open_many(glob_pattern, varname):
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(
            f"No files matched:\n  {glob_pattern}\n"
            "Check the folder path and the variable label in the filename."
        )
    print(f"  {varname}: {len(files)} file(s) "
          f"[{os.path.basename(files[0])} ... {os.path.basename(files[-1])}]")
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_mfdataset(files, combine="by_coords", chunks={"time": 120},
                           decode_times=time_coder, data_vars="minimal",
                           coords="minimal", compat="override")
    if varname not in ds:
        cands = [v for v in ds.data_vars if ds[v].ndim >= 3]
        varname = cands[0]
    return ds[varname]


def monthly_regression(mean_da, target_da):
    """Per cell, per calendar month, OLS fit target = a + b*mean across years.
       Vectorised over the grid (no python loop). Returns a(month,...), b(month,...)."""
    mbar = mean_da.groupby("time.month").mean("time")            # E[x]
    tbar = target_da.groupby("time.month").mean("time")          # E[y]
    xy   = (mean_da * target_da).groupby("time.month").mean("time")  # E[xy]
    xx   = (mean_da * mean_da).groupby("time.month").mean("time")    # E[x^2]
    cov = xy - mbar * tbar
    var = xx - mbar * mbar
    # Guard zero/near-zero variance cells (poles, masked, never-varying):
    # where piControl tas does not move across years there is no slope to
    # learn -> b=0 (extreme tracks the mean by the climatological offset a).
    # Avoids inf/nan leaking into the output. Defensible: stationary offset.
    b = xr.where(np.abs(var) > 1e-12, cov / var, 0.0).rename("slope")
    a = (tbar - b * mbar).rename("intercept")
    return a, b


def grids_match(da_a, da_b):
    """True if lat/lon coords are identical (same cells)."""
    for d in ("lat", "lon"):
        if d not in da_a.coords or d not in da_b.coords:
            return False
        if da_a[d].shape != da_b[d].shape:
            return False
        if not np.allclose(da_a[d].values, da_b[d].values):
            return False
    return True


def regrid_coeff_to(coeff, target_da):
    """Nearest-neighbour put a (month,lat,lon) coeff onto target's grid.
       Used only if the hosing grid differs from piControl."""
    return coeff.interp(lat=target_da["lat"], lon=target_da["lon"], method="nearest")


def apply_regression(a, b, hos_mean_da):
    months = hos_mean_da["time.month"]
    out = (a.sel(month=months) + b.sel(month=months) * hos_mean_da)
    return out.drop_vars("month", errors="ignore")


def write_var(template_da, data_da, varname, out_path):
    out = data_da.to_dataset(name=varname)
    longname = ("Daily Maximum Near-Surface Air Temperature" if varname == "tasmax"
                else "Daily Minimum Near-Surface Air Temperature")
    out[varname].attrs.update(standard_name="air_temperature",
                              long_name=longname,
                              units=template_da.attrs.get("units", "K"))
    out.attrs["reconstruction_method"] = (
        "Per-cell per-calendar-month OLS regression of extreme on mean fitted "
        "on EC-Earth3 piControl, applied to hosing tas (Weedon et al. 2010). "
        "Assumes stationary mean-extreme relationship.")
    enc = {v: {"zlib": True, "complevel": 4} for v in out.data_vars}
    out.to_netcdf(out_path, encoding=enc)
    print(f"    wrote {os.path.basename(out_path)}")


# ---------------------------------------------------------------------
# 3. MAIN
# ---------------------------------------------------------------------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Opening piControl mean / max / min ...")
    pic_tas    = open_many(PIC_TAS_GLOB,    "tas")
    pic_tasmax = open_many(PIC_TASMAX_GLOB, "tasmax")
    pic_tasmin = open_many(PIC_TASMIN_GLOB, "tasmin")

    pic_tas, pic_tasmax, pic_tasmin = xr.align(pic_tas, pic_tasmax, pic_tasmin,
                                               join="inner")
    print(f"  paired piControl steps: {pic_tas.sizes.get('time','NA')}")

    print("Fitting per-cell per-month regressions on piControl ...")
    a_max, b_max = monthly_regression(pic_tas, pic_tasmax)
    a_min, b_min = monthly_regression(pic_tas, pic_tasmin)

    for label, tas_glob in HOSING_EXPERIMENTS.items():
        print(f"\n=== Hosing experiment: {label} ===")
        try:
            hos_tas = open_many(tas_glob, "tas")
        except FileNotFoundError as e:
            print("  SKIP:", e)
            continue

        # grid check -- regrid coefficients only if needed
        aM, bM, aN, bN = a_max, b_max, a_min, b_min
        if grids_match(pic_tas, hos_tas):
            print("  grid: identical to piControl (no regridding)")
            # Same cells, but on-disk lat/lon labels can differ by float dust
            # (45.0000001 vs 45.0). Snap coeff labels onto hosing's EXACTLY so
            # downstream xr.where / anomaly subtraction align with join='exact'.
            snap = {}
            for d in ("lat", "lon"):
                if d in hos_tas.coords:
                    snap[d] = hos_tas[d]
            def _snap(x):
                return x.assign_coords({d: snap[d] for d in snap if d in x.coords})
            aM, bM = _snap(a_max), _snap(b_max)
            aN, bN = _snap(a_min), _snap(b_min)
        else:
            print("  grid: DIFFERS from piControl -> regridding coefficients (nearest)")
            aM = regrid_coeff_to(a_max, hos_tas); bM = regrid_coeff_to(b_max, hos_tas)
            aN = regrid_coeff_to(a_min, hos_tas); bN = regrid_coeff_to(b_min, hos_tas)

        print("  applying regression ...")
        hos_tasmax = apply_regression(aM, bM, hos_tas)
        hos_tasmin = apply_regression(aN, bN, hos_tas)

        # enforce tasmax >= tas >= tasmin where a fit overshoots.
        # Operate on raw arrays (same shape/order) then rewrap with hosing
        # coords -- never triggers xarray's strict-align on float-dust labels.
        hos_tasmax = xr.DataArray(
            np.maximum(hos_tasmax.data, hos_tas.data),
            dims=hos_tas.dims, coords=hos_tas.coords, name="tasmax",
            attrs=hos_tasmax.attrs)
        hos_tasmin = xr.DataArray(
            np.minimum(hos_tasmin.data, hos_tas.data),
            dims=hos_tas.dims, coords=hos_tas.coords, name="tasmin",
            attrs=hos_tasmin.attrs)

        # Absolute reconstructed fields ONLY. The hosing-minus-piControl anomaly is
        # built by the CDO step (1.ECHEarth3nahos_anomaly_cdo_explicit.sh, ymonsub),
        # so there is a SINGLE anomaly owner; this script emits just the absolutes
        # that step consumes.
        write_var(pic_tasmax, hos_tasmax, "tasmax",
                  os.path.join(OUTPUT_DIR, f"tasmax_Amon_EC-Earth3_{label}_reconstructed.nc"))
        write_var(pic_tasmin, hos_tasmin, "tasmin",
                  os.path.join(OUTPUT_DIR, f"tasmin_Amon_EC-Earth3_{label}_reconstructed.nc"))

    print("\nDone. Files written to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
