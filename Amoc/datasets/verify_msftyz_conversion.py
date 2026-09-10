import numpy as np, xarray as xr

base = "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
msf = xr.open_dataset(f"{base}/msftyz_ipsl/msftyz_Omon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gn_201501-210012.nc")
ref = xr.open_dataset(f"{base}/terhaar_amoc/amoc/amoc/26.5N/ssp126/amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc")

v = msf["msftyz"]                       # (time, 3basin, olevel, y, x)
basin = msf["3basin"].values            # [1,2,3]
atl = int(np.where(basin == 2)[0][0])   # Atlantic index

# nearest lat 26.5N on nav_lat(y,x); x=1
nav_lat = msf["nav_lat"].isel(x=0).values   # (y,)
yj = int(np.abs(nav_lat - 26.5).argmin())
lat_eff = float(nav_lat[yj])
# reference stored the model latitude it actually used:
yj_ref = int(np.abs(nav_lat - 26.95257).argmin())
print(f"[nearest 26.5 -> y={yj} lat={lat_eff:.5f}] "
      f"[ref-row y={yj_ref} lat={nav_lat[yj_ref]:.5f}] "
      f"neighbours={np.round(nav_lat[yj_ref-1:yj_ref+2],5)}")
import sys
if "--refrow" in sys.argv:
    yj, lat_eff = yj_ref, float(nav_lat[yj_ref])

col = v.isel({"3basin": atl, "y": yj, "x": 0})   # (time, olevel)
col = col.where(col != 1e20)

olevel = msf["olevel"].values
deep = olevel > 500

# ORDER 1: annual-mean over months, then max over depth
ann = col.groupby(col.time.dt.year).mean("time")   # (year, olevel)
raw_full = ann.max("olevel").values
raw_deep = ann.isel(olevel=deep).max("olevel").values
yrs = ann["year"].values

# ORDER 2: max over depth per month, then annual-mean
mmax_full = col.max("olevel").groupby(col.time.dt.year).mean("time").values
mmax_deep = col.isel(olevel=deep).max("olevel").groupby(col.time.dt.year).mean("time").values

def conv(x, mode):
    return x/1e6 if mode == "a" else x/1026.0/1e6

# reference series aligned to msftyz years (2015-2100)
ry = ref["year"].values
ramoc = ref["amoc"].isel(x=0).values
rmap = {int(y): float(a) for y, a in zip(ry, ramoc)}
rvec = np.array([rmap[int(y)] for y in yrs])

print(f"units(msftyz)       = {v.attrs.get('units')!r}")
print(f"standard_name       = {v.attrs.get('standard_name')!r}")
print(f"Atlantic 3basin idx = {atl} (value {basin[atl]:.0f})")
print(f"lat nearest 26.5N   = {lat_eff:.5f}  (ref nav_lat 26.95257)")
print(f"years overlap       = {yrs.min()}-{yrs.max()} ({len(yrs)})")
print(f"ref amoc mean/range = {rvec.mean():.2f} / {rvec.min():.2f}-{rvec.max():.2f} Sv\n")

series = [
    ("O1 annual-mean then FULL max", raw_full),
    ("O1 annual-mean then max>500m", raw_deep),
    ("O2 monthly FULL max then mean", mmax_full),
    ("O2 monthly max>500m then mean", mmax_deep),
]
for name, raw in series:
    print(f"--- {name} ---")
    for mode, tag in [("a", "a) /1e6        [m3/s]"), ("b", "b) /1026/1e6   [kg/s]")]:
        s = conv(raw, mode)
        r = np.corrcoef(s, rvec)[0, 1]
        md = float(np.mean(s - rvec))
        rmse = float(np.sqrt(np.mean((s - rvec)**2)))
        print(f"  {tag}: mean={s.mean():7.2f} Sv | corr={r:.5f} | mean_diff={md:+.3f} | rmse={rmse:.3f}")
    print()
