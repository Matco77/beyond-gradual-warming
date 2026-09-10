import numpy as np, xarray as xr

base = "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
lb_h = xr.open_dataset(f"{base}/lebars_amoc/vo/cmip6_amoc_vo/cmip6_amoc_vo_historical_EC-Earth3_1850_2014.nc")
lb_s = xr.open_dataset(f"{base}/lebars_amoc/vo/cmip6_amoc_vo/cmip6_amoc_vo_ssp126_EC-Earth3_2015_2100.nc")
th_h = xr.open_dataset(f"{base}/terhaar_amoc/amoc/amoc/26.5N/historical/amoc_historical_EC-Earth-Consortium_EC-Earth3_r1i1p1f1.nc")

print("Le Bars latitudes:", lb_h["latitude"].values, "| model:", lb_h["model"].values)
print("Terhaar lat:", float(th_h["lat"].values),
      "| sector:", bytes(th_h["sector"].values).decode(errors="replace").strip())

# pick 26N latitude in Le Bars
lats = lb_h["latitude"].values
li = int(np.abs(lats - 26).argmin())
print(f"Le Bars picked latitude index {li} = {lats[li]}N\n")

def lb_series(ds):
    a = ds["amoc"].isel(model=0, latitude=li).values
    yr = np.array([int(str(t)[:4]) if not np.issubdtype(ds['time'].dtype, np.floating)
                   else int(t) for t in ds["time"].values])
    # time may be decimal year or datetime; handle below
    return a, ds

lb_h_a = lb_h["amoc"].isel(model=0, latitude=li).values
lb_s_a = lb_s["amoc"].isel(model=0, latitude=li).values
# Le Bars time -> year
def to_year(ds):
    t = ds["time"].values
    if np.issubdtype(t.dtype, np.datetime64):
        return t.astype('datetime64[Y]').astype(int) + 1970
    return np.floor(t).astype(int)   # time is year+0.5 (mid-year)
lb_h_y = to_year(lb_h); lb_s_y = to_year(lb_s)

th_a = th_h["amoc"].values
th_y = th_h["year"].values.astype(int)

# align historical 1850-2014
lbm = {int(y): float(v) for y, v in zip(lb_h_y, lb_h_a)}
thm = {int(y): float(v) for y, v in zip(th_y, th_a)}
yrs = np.array(sorted(set(lbm) & set(thm)))
L = np.array([lbm[y] for y in yrs]); T = np.array([thm[y] for y in yrs])

def rm(x, w=5): return np.convolve(x, np.ones(w)/w, 'valid')
print("=== TEST 1: HISTORICAL 1850-2014 (Le Bars vs Terhaar) ===")
print(f"years {yrs.min()}-{yrs.max()} (n={len(yrs)})")
print(f"mean  Le Bars={L.mean():.2f} Sv   Terhaar={T.mean():.2f} Sv")
print(f"mean_diff (LB-TH) = {np.mean(L-T):+.3f} Sv")
print(f"rmse              = {np.sqrt(np.mean((L-T)**2)):.3f} Sv")
print(f"corr raw          = {np.corrcoef(L,T)[0,1]:.4f}")
print(f"corr 5yr-smoothed = {np.corrcoef(rm(L),rm(T))[0,1]:.4f}\n")

# TEST 3: continuity hist->ssp126
lbsm = {int(y): float(v) for y, v in zip(lb_s_y, lb_s_a)}
m0514 = np.mean([lbm[y] for y in range(2005,2015)])
m1524 = np.mean([lbsm[y] for y in range(2015,2025)])
print("=== TEST 3: CONTINUITY historical->ssp126 (Le Bars) ===")
print(f"mean 2005-2014 (hist) = {m0514:.3f} Sv")
print(f"mean 2015-2024 (ssp)  = {m1524:.3f} Sv")
print(f"jump                  = {m1524-m0514:+.3f} Sv\n")

# TEST 4: RAPID
print("=== TEST 4: PLAUSIBILITY vs RAPID ~17 Sv @26N ===")
print(f"Le Bars mean 2005-2014 = {m0514:.3f} Sv | RAPID~17 | diff = {m0514-17:+.3f} Sv")
print(f"Terhaar mean 2005-2014 = {np.mean([thm[y] for y in range(2005,2015)]):.3f} Sv")
