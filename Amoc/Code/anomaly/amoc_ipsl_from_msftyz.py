#!/usr/bin/env python3
# =====================================================================
# IPSL-CM6A-LR AMOC at ~26N from the raw msftyz field: historical, ssp126, ssp370.
#
# WHY ALL THREE FROM ONE PIPELINE, instead of taking Terhaar's published series.
#   Terhaar et al. published historical / piControl / ssp126 / ssp585 but NOT ssp370,
#   so the high-emissions branch had to be derived here anyway. Deriving only ssp370
#   and keeping Terhaar for ssp126 would mean the two scenarios of the SAME model are
#   not the same diagnostic, and a ssp126-vs-ssp370 comparison would partly measure
#   the change of method. So every IPSL series this project uses is built here, with
#   one pipeline, from the standard ScenarioMIP r1i1p1f1 2015-2100 files.
#
# METHOD
#     msftyz, Atlantic basin (3basin == 2)
#       -> row: nav_lat nearest 26.5N  (y = 228, lat 26.056N on the IPSL grid)
#       -> annual mean of the 12 months
#       -> max over depth levels DEEPER THAN 500 m   (Terhaar / NAHosMIP M26 convention)
#       -> / 1026 / 1e6                              (kg s-1 -> Sv)
#
# VALIDATION GATE (--self-check, and run automatically before writing).
#   On the HISTORICAL period this pipeline reproduces Terhaar's own IPSL series to
#   corr 0.9971, bias -0.117 Sv, debiased rmse 0.081 Sv over 1850-2014. That is the
#   evidence the extraction is the same diagnostic Terhaar computes.
#
#   It does NOT reproduce Terhaar's ssp126 as closely (corr 0.78, debiased rmse 0.87 Sv).
#   That is a DATA difference, not a method one: Terhaar's ssp126 file runs to 2214 —
#   the extended ssp126 run — while the standard ScenarioMIP file used here ends in 2100,
#   and their 2071-2100 means differ by 0.64 Sv (9.43 vs 10.07). Checked, not assumed:
#   a lag scan (-3..+3 yr) puts the best match at lag 0, so it is not a year-axis
#   misalignment, and Terhaar's ssp126 is visibly less variable (sd of first differences
#   0.90 vs 1.18 Sv) while his historical matches ours (1.25 vs 1.23). Two different runs.
#
# ANOMALIES. Consumers take amoc(scenario year) - mean(amoc(historical)), both from the
# files this script writes, so the subtraction never crosses methods and the constant
# offset against Terhaar (-0.117 Sv) cancels inside it.
#
# OUTPUT  Amoc/datasets/amoc_ipsl_msftyz_26N_<exp>_<y0>_<y1>.nc
#   vars: year (int), amoc (Sv) — the schema isimip_bin_fields.R::read_terhaar reads.
#
# USAGE
#   ./.amoc_venv/bin/python Amoc/Code/anomaly/amoc_ipsl_from_msftyz.py            # all three
#   ./.amoc_venv/bin/python Amoc/Code/anomaly/amoc_ipsl_from_msftyz.py --self-check
# =====================================================================
import argparse
import datetime
import os
import sys

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
D = os.path.join(ROOT, "Amoc", "datasets")

RAW = {
    "historical": "msftyz_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc",
    "ssp126":     "msftyz_Omon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "ssp370":     "msftyz_Omon_IPSL-CM6A-LR_ssp370_r1i1p1f1_gn_201501-210012.nc",
}
YEARS = {"historical": (1850, 2014), "ssp126": (2015, 2100), "ssp370": (2015, 2100)}
TERHAAR_HIST = os.path.join(D, "terhaar_amoc", "amoc", "amoc", "26.5N", "historical",
                            "amoc_historical_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc")
LAT_SEL, DEPTH_MIN, RHO = 26.5, 500.0, 1026.0
GATE_CORR, GATE_RMSE = 0.99, 0.20        # historical agreement with Terhaar; observed 0.9971 / 0.081


def series(exp):
    """Annual AMOC-at-26N (Sv) for one experiment."""
    ds = xr.open_dataset(os.path.join(D, "msftyz_ipsl", RAW[exp]))
    atl = int(np.where(ds["3basin"].values == 2)[0][0])
    nav_lat = ds["nav_lat"].isel(x=0).values
    yj = int(np.abs(nav_lat - LAT_SEL).argmin())
    col = ds["msftyz"].isel({"3basin": atl, "y": yj, "x": 0})
    col = col.where(col != 1e20)
    deep = ds["olevel"].values > DEPTH_MIN
    ann = col.groupby(col.time.dt.year).mean("time")
    return (ann["year"].values.astype(int),
            ann.isel(olevel=deep).max("olevel").values / RHO / 1e6,
            float(nav_lat[yj]))


def self_check(verbose=True):
    """Gate: the historical series must match Terhaar's own IPSL historical."""
    y, a, _ = series("historical")
    ref = xr.open_dataset(TERHAAR_HIST)
    rmap = dict(zip(ref["year"].values.astype(int),
                    np.asarray(ref["amoc"].squeeze().values, dtype="float64")))
    k = [i for i, yy in enumerate(y) if int(yy) in rmap]
    x = a[k]; r = np.array([rmap[int(y[i])] for i in k])
    bias = float(np.mean(x - r))
    rmse_db = float(np.sqrt(np.mean((x - bias - r) ** 2)))
    corr = float(np.corrcoef(x, r)[0, 1])
    if verbose:
        print(f"self-check vs Terhaar historical (n={len(r)}): "
              f"corr {corr:.4f} | bias {bias:+.3f} Sv | debiased rmse {rmse_db:.3f} Sv")
    ok = corr >= GATE_CORR and rmse_db <= GATE_RMSE
    if not ok:
        raise SystemExit(f"FAILED gate: need corr >= {GATE_CORR} and rmse <= {GATE_RMSE}")
    if verbose:
        print("PASS")
    return corr, bias, rmse_db


def write(exp, y, a, lat_used, gate):
    y0, y1 = YEARS[exp]
    corr, bias, rmse_db = gate
    ds = xr.Dataset({"amoc": ("year", a)}, coords={"year": ("year", y)})
    ds["amoc"].attrs = dict(
        long_name=f"IPSL-CM6A-LR AMOC volume transport at {lat_used:.3f}N, maximum below 500 m ({exp})",
        units="Sv",
        method=("raw msftyz (Omon gn), Atlantic basin 3basin==2, nav_lat row nearest 26.5N, "
                "annual mean of the 12 months, then max over levels > 500 m, divided by 1026 and 1e6"),
        validation=(f"the historical series from this same pipeline matches Terhaar's published IPSL "
                    f"historical at corr {corr:.4f}, bias {bias:+.3f} Sv, debiased rmse {rmse_db:.3f} Sv "
                    f"over 1850-2014"),
        anomaly_note=("take anomalies against the historical file written by this same script, never "
                      "against Terhaar's, so the subtraction does not cross methods"),
        scenario_consistency=("historical, ssp126 and ssp370 are all built by this one pipeline from the "
                              "standard ScenarioMIP r1i1p1f1 2015-2100 files, so a ssp126-vs-ssp370 "
                              "comparison does not partly measure a change of method"))
    ds["year"].attrs = dict(long_name="calendar year")
    ds.attrs = dict(
        title=f"IPSL-CM6A-LR AMOC at 26N from raw msftyz ({exp})",
        source_file=RAW[exp],
        source_variable="msftyz (Omon, gn), IPSL-CM6A-LR r1i1p1f1",
        terhaar_ssp126_note=("Terhaar published no ssp370, and his ssp126 is a DIFFERENT run (extends to "
                             "2214, 2071-2100 mean 10.07 Sv vs 9.43 Sv here; best lag 0, so not a year "
                             "misalignment). His ssp126 is therefore not used by this project's IPSL branch."),
        generated_by=os.path.basename(__file__),
        created=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    p = os.path.join(D, f"amoc_ipsl_msftyz_26N_{exp}_{y0}_{y1}.nc")
    ds.to_netcdf(p)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-check", action="store_true", help="validate against Terhaar, no writes")
    a = ap.parse_args()
    gate = self_check()
    if a.self_check:
        return 0
    for exp in RAW:
        y, s, lat_used = series(exp)
        k = (y >= YEARS[exp][0]) & (y <= YEARS[exp][1])
        y, s = y[k], s[k]
        p = write(exp, y, s, lat_used, gate)
        tail = s[-30:].mean()
        print(f"{os.path.basename(p):48s} {y.min()}-{y.max()}  mean {s.mean():6.3f} Sv  "
              f"last-30yr {tail:6.3f} Sv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
