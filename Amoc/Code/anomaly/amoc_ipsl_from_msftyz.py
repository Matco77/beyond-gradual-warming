#!/usr/bin/env python3
# =====================================================================
# IPSL-CM6A-LR AMOC at ~26N from the raw msftyz field, for the scenarios
# Terhaar et al. did NOT publish (Terhaar has historical / ssp126 / ssp585
# only; this project needs ssp370 for the ISIMIP high-emissions branch).
#
# METHOD -- reproduce the number in
#   terhaar_amoc/.../26.5N/ssp126/amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc
# from raw msftyz (Omon, gn):
#     msftyz, Atlantic basin (3basin == 2)
#       -> row: nav_lat nearest 26.5N  (y = 228, lat 26.056N on the IPSL grid)
#       -> annual mean over the 12 months           (Le Bars "O1" order)
#       -> max over depth levels DEEPER THAN 500 m  (Terhaar / NAHosMIP convention)
#       -> / 1026 / 1e6                             (kg s-1 -> Sv)
#   verify_msftyz_conversion.py picked this O1 / >500 m / density combination:
#   against Terhaar ssp126 (2015-2100) it gives rmse 0.889 Sv, mean bias -0.20 Sv.
#
# BIAS CORRECTION. The ~-0.20 Sv bias is a near-constant offset (grid geometry,
# basin box vs model mask, the parameterised transport msftyz includes but this
# row-max approximation partly drops). It cancels in an anomaly taken against a
# baseline computed the SAME way -- but isimip_bin_fields.R takes the IPSL ssp370
# anomaly against Terhaar's HISTORICAL series, so the offset would NOT cancel.
# We therefore add  offset = mean(Terhaar_ssp126 - mymethod_ssp126)  over the
# 2015-2100 overlap, putting the ssp370 series on Terhaar's level to first order.
# The correction is one scalar, written into the file's attributes.
#
# OUTPUT  Amoc/datasets/amoc_ipsl_msftyz_26N_<scen>_2015_2100.nc
#   vars: year (int), amoc (Sv) -- the schema isimip_bin_fields.R::read_terhaar reads.
#
# USAGE   ./.amoc_venv/bin/python Amoc/Code/anomaly/amoc_ipsl_from_msftyz.py ssp370
# =====================================================================
import datetime
import os
import sys

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
D = os.path.join(ROOT, "Amoc", "datasets")

RAW = os.path.join(D, "msftyz_ipsl", "msftyz_Omon_IPSL-CM6A-LR_{scen}_r1i1p1f1_gn_201501-210012.nc")
TERHAAR_SSP126 = os.path.join(D, "terhaar_amoc", "amoc", "amoc", "26.5N", "ssp126",
                              "amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc")
LAT_SEL = 26.5
DEPTH_MIN = 500.0
RHO = 1026.0


def amoc_series(scen):
    """Annual AMOC-at-26N (Sv) from raw msftyz for one scenario: O1 order, max below 500 m."""
    ds = xr.open_dataset(RAW.format(scen=scen))
    basin = ds["3basin"].values
    atl = int(np.where(basin == 2)[0][0])
    nav_lat = ds["nav_lat"].isel(x=0).values
    yj = int(np.abs(nav_lat - LAT_SEL).argmin())
    col = ds["msftyz"].isel({"3basin": atl, "y": yj, "x": 0})   # (time, olevel)
    col = col.where(col != 1e20)
    deep = ds["olevel"].values > DEPTH_MIN
    ann = col.groupby(col.time.dt.year).mean("time")            # (year, olevel)
    amoc = ann.isel(olevel=deep).max("olevel").values / RHO / 1e6
    years = ann["year"].values.astype(int)
    return years, amoc, float(nav_lat[yj])


def offset_to_terhaar():
    """mean(Terhaar_ssp126 - mymethod_ssp126) over the 2015-2100 overlap."""
    y126, a126, _ = amoc_series("ssp126")
    ref = xr.open_dataset(TERHAAR_SSP126)
    ry = ref["year"].values.astype(int)
    ra = np.asarray(ref["amoc"].squeeze().values, dtype="float64")
    rmap = {int(y): float(v) for y, v in zip(ry, ra)}
    both = [(rmap[int(y)], a) for y, a in zip(y126, a126) if int(y) in rmap]
    r = np.array([t for t, _ in both]); m = np.array([x for _, x in both])
    return float(np.mean(r - m)), len(both)


def main():
    scen = sys.argv[1] if len(sys.argv) > 1 else "ssp370"
    if not os.path.exists(RAW.format(scen=scen)):
        sys.exit(f"raw msftyz not found: {RAW.format(scen=scen)}")

    off, n = offset_to_terhaar()
    years, amoc_raw, lat_used = amoc_series(scen)
    amoc = amoc_raw + off

    out = os.path.join(D, f"amoc_ipsl_msftyz_26N_{scen}_2015_2100.nc")
    ds = xr.Dataset(
        {"amoc": ("year", amoc), "amoc_uncorrected": ("year", amoc_raw)},
        coords={"year": ("year", years)})
    ds["amoc"].attrs = dict(
        long_name=f"IPSL-CM6A-LR AMOC volume transport at {lat_used:.3f}N, maximum below 500 m ({scen})",
        units="Sv",
        method=("raw msftyz (Omon gn), Atlantic basin 3basin==2, nav_lat row nearest 26.5N, "
                "annual mean of the 12 months, then max over levels > 500 m, divided by 1026 and 1e6"),
        bias_correction_Sv=f"{off:+.4f}",
        bias_correction_note=(f"added scalar = mean(Terhaar_ssp126 - this_method_ssp126) over {n} "
                              "overlapping years 2015-2100, to place the series on Terhaar's level "
                              "so the anomaly against Terhaar's historical series in isimip_bin_fields.R "
                              "is method-consistent to first order"),
        validation=("this O1 / >500 m / density combination reproduces Terhaar ssp126 to rmse "
                    "0.889 Sv (verify_msftyz_conversion.py); correlation 0.78, so it is an "
                    "approximation of Terhaar's section, adequate for a 1-Sv-bin anomaly"))
    ds["amoc_uncorrected"].attrs = dict(long_name="same series before the bias correction", units="Sv")
    ds["year"].attrs = dict(long_name="calendar year")
    ds.attrs = dict(
        title=f"IPSL-CM6A-LR AMOC at 26N from raw msftyz ({scen}), Terhaar-calibrated",
        source_file=os.path.basename(RAW.format(scen=scen)),
        source_variable="msftyz (Omon, gn), IPSL-CM6A-LR r1i1p1f1, ESGF v20190119",
        reference_series=os.path.relpath(TERHAAR_SSP126, D),
        generated_by=os.path.basename(__file__),
        created=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    ds.to_netcdf(out)
    print(f"wrote {os.path.basename(out)} | {years.min()}-{years.max()} | "
          f"offset {off:+.3f} Sv (n={n}) | mean {amoc.mean():.2f} Sv | "
          f"2015-20 {amoc[:6].mean():.2f} -> 2091-2100 {amoc[-10:].mean():.2f} Sv")


if __name__ == "__main__":
    main()
