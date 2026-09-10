#!/usr/bin/env python3
# =====================================================================
# Reconstruct monthly tasmin and tasmax for IPSL-CM6A-LR NAHosMIP (hosing)
# from IPSL-CM6A-LR piControl  --  clone of 1.reconstruct_tasmin_tasmax_ECEarth3.py,
# same method, IPSL paths/names.
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
# together, etc.), so piControl 1850-2349 teaches hosing 1850-....
#
# Grounding: Weedon et al. (2010) WATCH Forcing Data; stationary
# mean<->extreme statistical link.  Everything MONTHLY, native gr cell.
# If hosing grid differs from piControl, the fitted COEFFICIENTS are
# regridded onto the hosing grid (nearest-neighbour) -- no data invented.
#
# Numerics: moments are computed on CENTRED anomalies in float64.
# The naive E[x^2]-E[x]^2 on ~300 K float32 data loses the whole
# interannual variance (< 1 K^2) to catastrophic cancellation;
# (x - xbar)*(y - ybar) products do not.
#
# Self-check (synthetic, no data needed):
#   python 1.reconstruct_tasmin_tasmax_IPSL.py --selfcheck
# =====================================================================

import os
import sys
import glob
import numpy as np
import xarray as xr

# ---------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------
PICONTROL_DIR = ("/Users/Bova/Library/CloudStorage/"
                 "OneDrive-UniversitàCommercialeLuigiBocconi/"
                 "1.Tesi/Amoc/datasets/CMIP6_piControl/IPSL-CM6A-LR")

HOSING_DIR = ("/Users/Bova/Library/CloudStorage/"
              "OneDrive-UniversitàCommercialeLuigiBocconi/"
              "1.Tesi/Amoc/datasets/NAHosMIP/IPSL-CM6A-LR")

OUTPUT_DIR = HOSING_DIR        # write reconstructed files beside the hosing tas

# piControl files that teach the regression (Amon = monthly).
# Variant pinned to r1i1p1f1: the folder holds exactly one member; if a
# second member ever lands here the glob must not silently blend them.
PIC_TAS_GLOB    = os.path.join(PICONTROL_DIR, "tas_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_*.nc")
PIC_TASMAX_GLOB = os.path.join(PICONTROL_DIR, "tasmax_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_*.nc")
PIC_TASMIN_GLOB = os.path.join(PICONTROL_DIR, "tasmin_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_*.nc")

# Hosing experiments to process: label -> glob for its monthly tas file(s).
# IPSL has u03-hos only (no g01 in the NAHosMIP atmospheric tar). The archive
# filename uses the DRS token "u03-hos"; the label key stays "hos-u03-hos" so the
# reconstructed output name matches the EC-Earth3 / downstream convention.
HOSING_EXPERIMENTS = {
    "hos-u03-hos": os.path.join(HOSING_DIR, "tas_Amon_IPSL-CM6A-LR_u03-hos_*.nc"),
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
        raise KeyError(
            f"'{varname}' not found in {os.path.basename(files[0])}; "
            f"data variables present: {list(ds.data_vars)}")
    da = ds[varname]
    units = da.attrs.get("units")
    if units != "K":
        raise ValueError(f"{varname}: expected units 'K', got {units!r}")
    return da


def monthly_regression(mean_da, target_da):
    """Per cell, per calendar month, OLS fit target = a + b*mean across years.
       Vectorised over the grid (no python loop). Returns a(month,...), b(month,...).
       Centred moments in float64 -- see header note on cancellation."""
    x = mean_da.astype("float64")
    y = target_da.astype("float64")
    mbar = x.groupby("time.month").mean("time")                   # E[x]
    tbar = y.groupby("time.month").mean("time")                   # E[y]
    xa = x.groupby("time.month") - mbar                           # x - E[x]
    ya = y.groupby("time.month") - tbar                           # y - E[y]
    cov = (xa * ya).groupby("time.month").mean("time")
    var = (xa * xa).groupby("time.month").mean("time")
    # Guard cells where piControl tas barely moves across years (var below
    # 1e-6 K^2 over ~500 years = effectively constant, e.g. masked/fill).
    # Fallback b=1, a = tbar - mbar: the extreme TRACKS the mean shifted by
    # the climatological offset. (b=0 would freeze the cell at piControl
    # climatology -- visibly wrong under a several-K hosing signal.)
    guarded = var <= 1e-6
    b = xr.where(guarded, 1.0, cov / var).rename("slope")
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
       Used only if the hosing grid differs from piControl.
       fill_value='extrapolate': target cells outside the source range take
       the nearest edge value instead of silently becoming NaN."""
    return coeff.interp(lat=target_da["lat"], lon=target_da["lon"],
                        method="nearest", kwargs={"fill_value": "extrapolate"})


def apply_regression(a, b, hos_mean_da):
    """Group-wise a[month] + b[month]*tas. The groupby arithmetic keeps the
       coefficients at 12 months (never materialises a per-timestep copy)."""
    out = (hos_mean_da.groupby("time.month") * b).groupby("time.month") + a
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
        "on IPSL-CM6A-LR piControl, applied to hosing tas (Weedon et al. 2010). "
        "Assumes stationary mean-extreme relationship.")
    out.attrs["fit_source"] = (
        f"IPSL-CM6A-LR piControl r1i1p1f1 gr, "
        f"{int(template_da['time.year'][0])}-{int(template_da['time.year'][-1])}, "
        f"{template_da.sizes['time']} months")
    enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32"}
           for v in out.data_vars}
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

    nt_before = (pic_tas.sizes["time"], pic_tasmax.sizes["time"], pic_tasmin.sizes["time"])
    grid_before = (pic_tas.sizes["lat"], pic_tas.sizes["lon"])
    pic_tas, pic_tasmax, pic_tasmin = xr.align(pic_tas, pic_tasmax, pic_tasmin,
                                               join="inner")
    print(f"  paired piControl steps: {pic_tas.sizes['time']} "
          f"(tas/tasmax/tasmin had {nt_before[0]}/{nt_before[1]}/{nt_before[2]})")
    if (pic_tas.sizes["lat"], pic_tas.sizes["lon"]) != grid_before:
        raise RuntimeError(
            "inner align dropped grid cells: piControl variables carry "
            "different lat/lon labels -- fix the inputs, do not fit on a hole-y grid")

    print("Fitting per-cell per-month regressions on piControl ...")
    a_max, b_max = monthly_regression(pic_tas, pic_tasmax)
    a_min, b_min = monthly_regression(pic_tas, pic_tasmin)
    # One pass over piControl for all four coefficient maps (25 MB in memory);
    # without this, every output file written below would re-read ~1000 files.
    coeffs = xr.Dataset({"a_max": a_max, "b_max": b_max,
                         "a_min": a_min, "b_min": b_min}).compute()
    a_max, b_max = coeffs["a_max"], coeffs["b_max"]
    a_min, b_min = coeffs["a_min"], coeffs["b_min"]
    for nm in ("b_max", "b_min"):
        bb = coeffs[nm]
        print(f"  {nm}: min {float(bb.min()):+.3f}, "
              f"median {float(bb.median()):+.3f}, max {float(bb.max()):+.3f}, "
              f"cells |b|>3: {int((np.abs(bb) > 3).sum())}, "
              f"guarded (b==1 fallback): {int((bb == 1.0).sum())}")

    for label, tas_glob in HOSING_EXPERIMENTS.items():
        print(f"\n=== Hosing experiment: {label} ===")
        try:
            hos_tas = open_many(tas_glob, "tas")
        except FileNotFoundError as e:
            print("  SKIP:", e)
            continue

        # grid check -- regrid coefficients only if needed
        if grids_match(pic_tas, hos_tas):
            print("  grid: identical to piControl (no regridding)")
            # Same cells, but on-disk lat/lon labels can differ by float dust
            # (45.0000001 vs 45.0). Snap coeff labels onto hosing's EXACTLY so
            # downstream groupby arithmetic aligns instead of inner-joining
            # away the whole grid.
            snap = {d: hos_tas[d] for d in ("lat", "lon") if d in hos_tas.coords}
            def _snap(x):
                return x.assign_coords({d: snap[d] for d in snap if d in x.coords})
            aM, bM = _snap(a_max), _snap(b_max)
            aN, bN = _snap(a_min), _snap(b_min)
        else:
            print("  grid: DIFFERS from piControl -> regridding coefficients (nearest)")
            aM = regrid_coeff_to(a_max, hos_tas); bM = regrid_coeff_to(b_max, hos_tas)
            aN = regrid_coeff_to(a_min, hos_tas); bN = regrid_coeff_to(b_min, hos_tas)

        print("  applying regression ...")
        # load() everything BEFORE to_netcdf: writing a live dask graph through
        # the netCDF4 backend deadlocks (dask thread lock vs HDF5 lock while the
        # input files are still open). Arrays are <= 0.7 GB, they fit in memory.
        hos_tas = hos_tas.load()
        pred_max = apply_regression(aM, bM, hos_tas).astype("float32").load()
        pred_min = apply_regression(aN, bN, hos_tas).astype("float32").load()

        # Audit: did the reconstruction invent NaNs (bad regrid / coefficient
        # hole), and how often does the ordering clamp below fire?
        nan_pred = int(np.isnan(pred_max.data).sum()) + int(np.isnan(pred_min.data).sum())
        nan_tas = int(np.isnan(hos_tas.data).sum())
        if nan_pred > 2 * nan_tas:
            raise RuntimeError(
                f"reconstruction introduced NaNs not present in hosing tas "
                f"(pred: {nan_pred}, tas: {nan_tas})")
        clamp_max = float((pred_max.data < hos_tas.data).mean())
        clamp_min = float((pred_min.data > hos_tas.data).mean())
        print(f"  clamp fired: tasmax {clamp_max:.3%}, tasmin {clamp_min:.3%} of points")

        # enforce tasmax >= tas >= tasmin where a fit overshoots.
        # Operate on raw arrays (same shape/order) then rewrap with hosing
        # coords -- never triggers xarray's strict-align on float-dust labels.
        hos_tasmax = xr.DataArray(
            np.maximum(pred_max.data, hos_tas.data),
            dims=hos_tas.dims, coords=hos_tas.coords, name="tasmax",
            attrs={"clamp_fraction": clamp_max})
        hos_tasmin = xr.DataArray(
            np.minimum(pred_min.data, hos_tas.data),
            dims=hos_tas.dims, coords=hos_tas.coords, name="tasmin",
            attrs={"clamp_fraction": clamp_min})

        # Absolute reconstructed fields ONLY. The hosing-minus-piControl anomaly is
        # built by the CDO step (the IPSL anomaly script, ymonsub), so there is a
        # SINGLE anomaly owner; this script emits just the absolutes it consumes.
        write_var(pic_tasmax, hos_tasmax, "tasmax",
                  os.path.join(OUTPUT_DIR, f"tasmax_Amon_IPSL-CM6A-LR_{label}_reconstructed.nc"))
        write_var(pic_tasmin, hos_tasmin, "tasmin",
                  os.path.join(OUTPUT_DIR, f"tasmin_Amon_IPSL-CM6A-LR_{label}_reconstructed.nc"))

    print("\nDone. Files written to:", OUTPUT_DIR)


# ---------------------------------------------------------------------
# 4. SELF-CHECK (synthetic; validates fit + apply + guard, no data files)
# ---------------------------------------------------------------------
def _selfcheck():
    rng = np.random.default_rng(0)
    time = xr.date_range("2000-01-01", periods=1200, freq="MS", use_cftime=True)
    coords = {"time": time, "lat": [0.0, 1.0], "lon": [0.0, 1.0]}
    months = np.array([t.month for t in time.values])

    # tropical-ocean regime: tiny interannual spread (0.1 K) around ~280 K --
    # exactly where the old uncentred float32 formula collapses.
    x = 280.0 + rng.normal(0.0, 0.1, size=(1200, 2, 2))
    x[:, 0, 0] = 280.0                                   # constant cell -> guard path
    b_true = 1.0 + 0.05 * months / 12.0
    a_true = 5.0 - 0.1 * months
    y = a_true[:, None, None] + b_true[:, None, None] * x

    X = xr.DataArray(x.astype("float32"), dims=("time", "lat", "lon"), coords=coords)
    Y = xr.DataArray(y.astype("float32"), dims=("time", "lat", "lon"), coords=coords)
    a, b = monthly_regression(X, Y)

    b_exp = xr.DataArray(1.0 + 0.05 * np.arange(1, 13) / 12.0, dims="month")
    a_exp = xr.DataArray(5.0 - 0.1 * np.arange(1, 13), dims="month")
    varying = b.isel(lat=[0, 1], lon=[1])                # skip the constant cell
    assert float(np.abs(varying - b_exp).max()) < 1e-2, "slope not recovered"
    assert float(np.abs(a.isel(lat=1, lon=1) - a_exp).max()) < 3.0, "intercept off"
    # guard cell: b must be exactly the fallback 1.0
    assert float(np.abs(b.isel(lat=0, lon=0) - 1.0).max()) == 0.0, "guard not applied"

    # apply on a shifted 'hosing' series: prediction must equal a + b*tas
    x2 = X - 8.0
    pred = apply_regression(a, b, x2)
    expected = a.sel(month=xr.DataArray(months, dims="time")) + \
               b.sel(month=xr.DataArray(months, dims="time")) * x2
    assert float(np.abs(pred - expected.drop_vars("month")).max()) < 1e-9, "apply mismatch"
    print("selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
