#!/usr/bin/env python3
# =====================================================================
# EC-Earth3 AMOC at ~26N from vo, maximum below 500 m
#
#   Regenerates:
#     Amoc/datasets/ecearth3_amoc26N_vo_below500m_{historical,ssp126}_*.nc
#     Amoc/datasets/ecearth3_amoc26p5N_vo_below500m_{historical,ssp126}_*.nc
#
# METHOD -- Compute_amoc_from_vo.py, D. Le Bars, github.com/dlebars/CMIP_SeaLevel
#
#       vo * zonal_length * vertical_length
#         -> sum over longitude
#         -> cumsum over lev (from the surface downward)
#         -> max over lev
#         -> / 1e6                                (m3/s -> Sv)
#
#   with ONE modification: the vertical maximum is taken over levels
#   with lev > 500 m instead of the whole water column, to match the
#   convention of Terhaar et al. ("msftyz maximum below 500m") and of
#   Jackson et al. NAHosMIP M26.  For EC-Earth3 the threshold turns out
#   to be inert -- the maximum sits at 697.3 m in every single year of
#   historical and ssp126 -- so it is a documentation fix, not a
#   numerical one.  It is applied anyway so the convention is explicit.
#
# GRID ROW.  Le Bars crops lon 280-355E / lat 20-40N, then picks the row
#   whose mean latitude is closest to 26N.  On the ORCA1 V-grid that is
#   j=187 (0-based), mean latitude 25.5990N.  The alternative j=188
#   (26.4955N) is the other face of the T cell that Terhaar labels
#   rlat=26.056N.  Both are produced.  j=187 is the RECOMMENDED one:
#   tested against msftyz at rlat=26.056N over 18 decadal historical
#   years it matches better on every metric (closure-corrected bias
#   -0.178 vs -0.436 Sv, rmse 0.202 vs 0.454, profile rmse 0-2000 m
#   0.989 vs 1.287).  j=188 is kept as a grid-row sensitivity bound
#   (~0.23 Sv on the level, ~0 on the anomaly).
#
# KNOWN BIAS.  The section does not close: Psi at the bottom is about
#   -1.55 Sv instead of 0 (grid-geometry approximation cos(lat)*RE*dlon
#   instead of the true e1v, partial bottom cells, longitude-box basin
#   instead of the model mask, and vo excluding the parameterized GM
#   bolus that msftyz includes).  That is ~70% of the -2.19 Sv offset
#   against Terhaar.  It is a CONSTANT bias -- 1850-1900 mean -1.642 Sv
#   vs 2091-2100 -1.697 Sv, t = -0.85, total drift over 250 years
#   -0.149 Sv against an AMOC decline of -4.0 Sv -- so it cancels in
#   anomalies.  It is written out as psi_bottom for auditing, and NOT
#   subtracted.  Starting from vmo instead of vo would remove most of
#   it; vmo is published for EC-Earth3 on ESGF.
#
# DATA.  vo Omon gn r1i1p1f1, ESGF v20200918, node esgf.ceda.ac.uk.
#   The full record is 45.7 GB.  This script never downloads it: it
#   pulls one grid row over the Atlantic crop via OPeNDAP
#   (vo[:, :, j, 208:283]), ~270 KB per file, ~50 GB -> ~50 MB.
#   Slices are cached under CACHE, so re-runs are free.
#
# USAGE
#   python3 amoc_from_vo_below500m.py --self-check   # reproduce Le Bars, no writes
#   python3 amoc_from_vo_below500m.py                # both rows -> Amoc/datasets/
#   python3 amoc_from_vo_below500m.py --row 187 --outdir /tmp
# =====================================================================

import argparse
import datetime
import glob
import os
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATASETS = os.path.join(ROOT, "Amoc", "datasets")
CACHE = os.environ.get("AMOC_VO_CACHE", os.path.join(DATASETS, "cache_vo_rows"))

ESGF_SEARCH = "https://esgf-data.dkrz.de/esg-search/search/"
ESGF_VERSION = "v20200918"
ESGF_NODE = "esgf.ceda.ac.uk"
DATASET_ID = {
    "historical": "CMIP6.CMIP.EC-Earth-Consortium.EC-Earth3.historical.r1i1p1f1.Omon.vo.gn." + ESGF_VERSION,
    "ssp126": "CMIP6.ScenarioMIP.EC-Earth-Consortium.EC-Earth3.ssp126.r1i1p1f1.Omon.vo.gn." + ESGF_VERSION,
    "ssp370": "CMIP6.ScenarioMIP.EC-Earth-Consortium.EC-Earth3.ssp370.r1i1p1f1.Omon.vo.gn." + ESGF_VERSION,
}
YEARS = {"historical": (1850, 2014), "ssp126": (2015, 2100), "ssp370": (2015, 2100)}
EXPERIMENTS = ("historical", "ssp126", "ssp370")

# Le Bars, Compute_amoc_from_vo.py
RE = 6.371e6                       # :133  radius of the Earth
LON_MIN, LON_MAX = 280, 355        # :32   Atlantic crop
LAT_MIN, LAT_MAX = 20, 40          # :33
LAT_SEL = 26                       # :29   latitude to select
DEPTH_MIN = 500.0                  # <-- the one modification

ROWS = {187: "26N", 188: "26p5N"}   # V-grid row (0-based) -> filename tag
I0, I1 = 208, 283                   # i-block of the Le Bars crop (0-based, end-exclusive)


# ---------------------------------------------------------------- ESGF

def esgf_urls(exp):
    """OPeNDAP endpoints for one experiment, sorted by start year."""
    import json
    import urllib.parse
    import urllib.request

    q = urllib.parse.urlencode({
        "type": "File", "project": "CMIP6", "source_id": "EC-Earth3",
        "experiment_id": exp, "variable_id": "vo", "table_id": "Omon",
        "variant_label": "r1i1p1f1", "format": "application/solr+json", "limit": 1000,
        "fields": "title,url",
    })
    with urllib.request.urlopen(ESGF_SEARCH + "?" + q, timeout=180) as r:
        docs = json.load(r)["response"]["docs"]
    urls = [u.split("|")[0][:-5] for d in docs for u in d["url"] if u.endswith("|OPENDAP")]
    if not urls:
        raise RuntimeError(f"ESGF returned no OPeNDAP endpoints for {exp}")
    return sorted(urls)


def grid():
    """Ocean grid (2D latitude/longitude, lev, lev_bnds), cached."""
    p = os.path.join(CACHE, "grid_ecearth3.nc")
    if not os.path.exists(p):
        os.makedirs(CACHE, exist_ok=True)
        url = esgf_urls("ssp126")[0]
        xr.open_dataset(url)[["latitude", "longitude", "lev", "lev_bnds"]].load().to_netcdf(p)
    return xr.open_dataset(p)


def fetch_rows(exp, j, workers=6):
    """Cache vo[:, :, j, I0:I1] for every file of one experiment."""
    out = os.path.join(CACHE, f"j{j}_{exp}")
    os.makedirs(out, exist_ok=True)
    urls = esgf_urls(exp)
    todo = [u for u in urls
            if not os.path.exists(os.path.join(out, u.split("/")[-1].replace(".nc", ".npz")))]
    if not todo:
        return out

    def grab(url):
        dst = os.path.join(out, url.split("/")[-1].replace(".nc", ".npz"))
        for att in range(5):
            try:
                with xr.open_dataset(url, decode_times=False) as ds:
                    v = ds["vo"].isel(j=j, i=slice(I0, I1)).load().values.astype("float32")
                np.savez_compressed(dst, vo=v)
                return
            except Exception:
                if att == 4:
                    raise
                time.sleep(3 * (att + 1))

    print(f"  fetching {len(todo)} files for {exp}, row j={j} "
          f"(~{0.27 * len(todo):.0f} MB) ...", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(grab, todo))
    print(f"  done in {time.time() - t0:.0f}s", flush=True)
    return out


# ------------------------------------------------------------ Le Bars

def section(g, j):
    """zonal_length and the section coordinates for V-grid row j."""
    lat = g.latitude.values[j, I0:I1]
    lon = g.longitude.values[j, I0:I1]
    dlon = lon[1:] - lon[:-1]
    dlon = np.append(dlon, dlon[-1])                                    # :136
    zonal_length = np.cos(np.radians(lat)) * RE * np.radians(dlon)      # :137
    return lat, lon, zonal_length


def default_row(g):
    """The row Le Bars' own code selects: crop, then nearest mean latitude to 26N."""
    lat, lon = g.latitude.values, g.longitude.values
    mask = ((lon > LON_MIN) & (lon < LON_MAX) & (lat >= LAT_MIN) & (lat <= LAT_MAX))
    jj = np.where(mask.any(axis=1))[0]
    ii = np.where(mask.any(axis=0))[0]
    avg_lat = lat[jj.min():jj.max() + 1, ii.min():ii.max() + 1].mean(axis=1)
    return int(jj.min() + np.abs(avg_lat - LAT_SEL).argmin()), int(ii.min()), int(ii.max() + 1)


def yearly(cache_dir):
    """Annual-mean vo section per year: unweighted mean of the 12 months, as in Le Bars."""
    years, out = [], []
    for f in sorted(glob.glob(os.path.join(cache_dir, "vo_Omon_EC-Earth3_*.npz"))):
        v = np.load(f)["vo"].astype("float64")
        v = np.where(np.abs(v) > 1e19, np.nan, v)
        out.append(np.nanmean(v, axis=0))
        years.append(int(os.path.basename(f).split("_gn_")[1][:4]))
    if not out:
        raise RuntimeError(f"no cached slices in {cache_dir}")
    return np.array(years), np.stack(out)


def streamfunction(vo_yr, zonal_length, vertical_length):
    """(year, lev, x) m/s -> (year, lev) Sv. Le Bars :148-154."""
    vol = vo_yr * zonal_length[None, None, :] * vertical_length[None, :, None]
    return np.nancumsum(np.nansum(vol, axis=2), axis=1) / 1e6           # :151, :154, :278


def amoc(psi, lev, depth_min=DEPTH_MIN):
    """Maximum below depth_min (None -> whole column, i.e. Le Bars verbatim :157)."""
    if depth_min is None:
        return np.nanmax(psi, axis=1), lev[np.nanargmax(psi, axis=1)]
    sel = lev > depth_min
    return np.nanmax(psi[:, sel], axis=1), lev[sel][np.nanargmax(psi[:, sel], axis=1)]


def series(exp, j, g):
    lat, lon, zl = section(g, j)
    lev = g.lev.values
    vl = g.lev_bnds.values[:, 1] - g.lev_bnds.values[:, 0]              # :145
    years, vo_yr = yearly(fetch_rows(exp, j))
    psi = streamfunction(vo_yr, zl, vl)
    deep, d_deep = amoc(psi, lev, DEPTH_MIN)
    full, d_full = amoc(psi, lev, None)
    return dict(years=years, deep=deep, full=full, depth=d_deep, depth_full=d_full,
                psi_bottom=psi[:, -1], lat=lat, lon=lon, lev=lev)


# ----------------------------------------------------------- validate

def self_check(g):
    """Reproduce Le Bars' published 26N historical series with no depth threshold.

    This is the gate: if the pipeline does not land on his numbers, nothing
    downstream is trustworthy.  Tolerance 1e-3 Sv; observed 1e-6.
    """
    ref = os.path.join(DATASETS, "lebars_amoc", "vo", "cmip6_amoc_vo",
                       "cmip6_amoc_vo_historical_EC-Earth3_1850_2014.nc")
    if not os.path.exists(ref):
        print(f"SKIP self-check: reference not found at {ref}")
        return True
    j, i0, i1 = default_row(g)
    assert (j, i0, i1) == (187, I0, I1), f"crop drifted: got j={j}, i={i0}:{i1}"
    s = series("historical", j, g)
    with xr.open_dataset(ref) as d:
        lb = d["amoc"].squeeze().sel(latitude=LAT_SEL).values
        assert (np.floor(d["time"].values).astype(int) == s["years"]).all(), "year axis mismatch"
    err = np.abs(s["full"] - lb).max()
    print(f"self-check: row j={j} ({s['lat'].mean():.4f}N), "
          f"max|mine - Le Bars| = {err:.6f} Sv over {len(lb)} years")
    assert err < 1e-3, f"FAILED: {err:.6f} Sv exceeds 1e-3"
    d_thr = np.abs(s["full"] - s["deep"]).max()
    print(f"            depth threshold shifts historical by {d_thr:.6f} Sv "
          f"(max depth {np.unique(np.round(s['depth_full'], 1))} m)")
    print("PASS")
    return True


# -------------------------------------------------------------- write

RATIONALE = {
    187: ("V-grid row selected by Le Bars' own criterion (Atlantic crop, then nearest mean "
          "latitude to 26N). RECOMMENDED. Tested against msftyz at rlat=26.056N over 18 decadal "
          "historical years it matches better than the alternative row j=188 on every metric: "
          "closure-corrected bias -0.178 vs -0.436 Sv, rmse 0.202 vs 0.454 Sv, profile rmse "
          "0-2000 m 0.989 vs 1.287 Sv."),
    188: ("V-grid row j=188, the northern face of the T cell that carries Terhaar's label "
          "rlat=26.056N; the companion file uses j=187, the southern face of the same T cell. "
          "Chosen a priori as the nearer row in latitude (0.4395 vs 0.4570 deg) and on the NEMO "
          "C-grid indexing argument. EMPIRICAL CHECK CONTRADICTS THAT CHOICE: comparing "
          "vo-derived streamfunction profiles against msftyz at rlat=26.056N over 18 decadal "
          "historical years, row j=187 matches better on every metric (closure-corrected bias "
          "-0.178 vs -0.436 Sv, rmse 0.202 vs 0.454 Sv, profile rmse 0-2000 m 0.989 vs 1.287 Sv). "
          "The companion file at j=187 is the better representation of Terhaar's section; this "
          "file is retained for comparison only."),
}

DESCRIPTION = (
    "AMOC volume transport computed by integrating the CMIP6 meridional velocity vo: "
    "vo * zonal_length * vertical_length, summed zonally over the Atlantic, cumulatively "
    "integrated in the vertical from the surface downward, then maximised over levels DEEPER "
    "THAN 500 m. Method follows Compute_amoc_from_vo.py from D. Le Bars, CMIP_SeaLevel "
    "(github.com/dlebars/CMIP_SeaLevel), with the single modification that the vertical maximum "
    "is restricted to lev > 500 m instead of the full water column, to match the convention of "
    "Terhaar et al. (\"msftyz maximum below 500m\") and of Jackson et al. NAHosMIP M26."
)


def write(exp, j, s, outdir):
    y0, y1 = YEARS[exp]
    tag = ROWS[j]
    other = ROWS[188 if j == 187 else 187]
    ds = xr.Dataset(
        {"amoc": ("time", s["deep"]),
         "amoc_full_column": ("time", s["full"]),
         "depth_of_max": ("time", s["depth"]),
         "psi_bottom": ("time", s["psi_bottom"])},
        coords={"time": ("time", s["years"].astype("float64") + 0.5)})
    ds["amoc"].attrs = dict(
        long_name=f"AMOC volume transport at {s['lat'].mean():.2f}N, maximum below 500 m",
        units="Sv",
        description=DESCRIPTION,
        latitude_used=(f"{s['lat'].mean():.4f} degrees_north (EC-Earth3 ORCA1 V-grid row "
                       f"j={j} 0-based / {j + 1} 1-based; row spans "
                       f"{s['lat'].min():.4f} to {s['lat'].max():.4f})"),
        grid_row_rationale=RATIONALE[j],
        depth_threshold="lev > 500 m (vertical maximum restricted to levels below 500 m)",
        atlantic_definition=(f"longitude box {s['lon'].min():.3f} to {s['lon'].max():.3f} "
                             "degrees_east (Le Bars crop 280-355E); not a model basin mask"),
        unit_conversion="m3 s-1 divided by 1e6; no density factor (vo is a velocity)")
    ds["amoc_full_column"].attrs = dict(
        long_name="AMOC volume transport, maximum over full water column", units="Sv",
        description=("Same pipeline without the depth threshold; at j=187 this reproduces Le Bars "
                     "cmip6_amoc_vo_*_EC-Earth3_*.nc at 26N to 1e-6 Sv. Kept for traceability."))
    ds["depth_of_max"].attrs = dict(long_name="Depth of the maximum used for amoc", units="m")
    ds["psi_bottom"].attrs = dict(
        long_name="Streamfunction at the deepest level (mass-closure residual)", units="Sv",
        description=(f"Psi at {s['lev'][-1]:.0f} m. Should be 0 for a closed section; departure "
                     "from 0 measures the mass imbalance of the vo integration (grid-geometry "
                     "approximation, partial cells, basin box, missing parameterized transport). "
                     "Constant in time (no significant drift), so it cancels in anomalies. "
                     "Diagnostic only, not subtracted from amoc."))
    ds["time"].attrs = dict(long_name="time", units="year",
                            description="year + 0.5 (mid-year), as in Le Bars")
    ds.attrs = dict(
        title=f"EC-Earth3 AMOC at {s['lat'].mean():.2f}N from vo, maximum below 500 m ({exp})",
        source_variable="vo (Omon, gn), EC-Earth3 r1i1p1f1",
        esgf_dataset=DATASET_ID[exp], esgf_version=ESGF_VERSION, esgf_data_node=ESGF_NODE,
        esgf_access=f"OPeNDAP subset, j={j + 1} (1-based), i={I0 + 1}:{I1} (1-based)",
        method_reference="Compute_amoc_from_vo.py, D. Le Bars, github.com/dlebars/CMIP_SeaLevel",
        method_modification="vertical maximum restricted to lev > 500 m",
        validation="full-column variant reproduces Le Bars 26N historical to 1e-6 Sv at row j=187",
        companion_file=f"ecearth3_amoc{other}_vo_below500m_{exp}_{y0}_{y1}.nc",
        generated_by=os.path.basename(__file__),
        created=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    if j == 188:
        ds.attrs["recommended_file"] = ("ecearth3_amoc26N_vo_below500m_*.nc (j=187) — better match "
                                        "to Terhaar; see amoc:grid_row_rationale")
    p = os.path.join(outdir, f"ecearth3_amoc{tag}_vo_below500m_{exp}_{y0}_{y1}.nc")
    ds.to_netcdf(p)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--row", type=int, choices=sorted(ROWS), action="append",
                    help="V-grid row (default: both)")
    ap.add_argument("--outdir", default=DATASETS)
    ap.add_argument("--self-check", action="store_true", help="validate against Le Bars, no writes")
    a = ap.parse_args()

    g = grid()
    if a.self_check:
        return 0 if self_check(g) else 1

    for j in (a.row or sorted(ROWS)):
        for exp in EXPERIMENTS:
            s = series(exp, j, g)
            print(f"{os.path.basename(write(exp, j, s, a.outdir))}  "
                  f"{s['lat'].mean():.4f}N  mean {s['deep'].mean():.3f} Sv  "
                  f"bottom {s['psi_bottom'].mean():+.3f} Sv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
