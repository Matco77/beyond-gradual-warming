import xarray as xr, numpy as np, glob, os, csv, re

BASE = "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
ABS = "assente"
rows = []

def attr(ds, var, key):
    v = ds[var].attrs.get(key)
    return v if v not in (None, "") else ABS

FILL = 9.969209968386869e+36   # NC_FILL_DOUBLE (reso come "_" da ncdump)
def stats(a):
    a = np.asarray(a, float)
    a = a[~np.isnan(a) & (a != FILL)]
    if a.size == 0: return ("assente","assente","assente","assente",0)
    return (round(float(a.mean()),4), round(float(a.std(ddof=1)),4),
            round(float(a.min()),4), round(float(a.max()),4), int(a.size))

def timeinfo(ds):
    tname = "time" if "time" in ds.variables else ("year" if "year" in ds.variables else None)
    if tname is None: return (ABS, "assente", "assente", "assente", ABS)
    t = ds[tname].values
    u = ds[tname].attrs.get("units", ABS) or ABS
    return (tname, int(t.size), round(float(t[0]),4), round(float(t[-1]),4), u)

def latinfo(ds):
    """Return 'value (varname)' for every stored latitude coordinate; else assente."""
    out = []
    for name in ds.variables:
        v = ds[name]
        sn = v.attrs.get("standard_name","")
        ln = v.attrs.get("long_name","")
        if name in ("lat","nav_lat","latitude") or sn in ("latitude","grid_latitude") \
           or "latitud" in ln.lower():
            vals = np.atleast_1d(v.values).astype(float)
            for x in vals:
                out.append(f"{x:.4f} ({name})")
    return "; ".join(out) if out else ABS

def add(source, model, exp, ds, var):
    tname,tn,t0,t1,tu = timeinfo(ds)
    mu,sd,mn,mx,nv = stats(ds[var].values)
    rows.append({
        "fonte": source, "modello": model, "esperimento": exp,
        "variabile": var,
        "units": attr(ds, var, "units"),
        "latitudine": latinfo(ds),
        "long_name_o_description": ds[var].attrs.get("long_name",
                                    ds[var].attrs.get("description", ABS)) or ABS,
        "time_var": tname, "n_time": tn, "time_primo": t0, "time_ultimo": t1,
        "time_units": tu,
        "n_valid": nv, "media": mu, "std": sd, "min": mn, "max": mx,
    })

# ---------- 1. M26 (all models, all internal vars) ----------
for f in sorted(glob.glob(f"{BASE}/nahosmip_m26/M26_*.nc")):
    model = re.sub(r"^M26_|\.nc$","", os.path.basename(f))
    ds = xr.open_dataset(f, decode_times=False)
    for var in ds.data_vars:
        if var == "time": continue
        m = re.match(r"^(con|hos|r20|r50|r70|r100)_", var)
        exp = m.group(1) if m else var
        add("M26", model, exp, ds, var)

# ---------- 2. Le Bars (EC-Earth3 historical + ssp126) ----------
for exp,fn in [("historical","cmip6_amoc_vo_historical_EC-Earth3_1850_2014.nc"),
               ("ssp126","cmip6_amoc_vo_ssp126_EC-Earth3_2015_2100.nc")]:
    ds = xr.open_dataset(f"{BASE}/lebars_amoc/vo/cmip6_amoc_vo/{fn}", decode_times=False)
    # one row per stored latitude value (dim latitude)
    lats = np.atleast_1d(ds["latitude"].values)
    latu = ds["latitude"].attrs.get("units", ABS) or ABS
    for i,lv in enumerate(lats):
        tname,tn,t0,t1,tu = timeinfo(ds)
        mu,sd,mn,mx,nv = stats(ds["amoc"].isel(model=0, latitude=i).values)
        rows.append({
            "fonte":"LeBars","modello":"EC-Earth3","esperimento":exp,"variabile":"amoc",
            "units": ds["amoc"].attrs.get("units", ABS) or ABS,
            "latitudine": f"{float(lv):.4f} (latitude)",
            "long_name_o_description": ds["amoc"].attrs.get("long_name",
                                        ds["amoc"].attrs.get("description",ABS)) or ABS,
            "time_var":tname,"n_time":tn,"time_primo":t0,"time_ultimo":t1,"time_units":tu,
            "n_valid":nv,"media":mu,"std":sd,"min":mn,"max":mx,
        })

# ---------- 3. Terhaar (EC-Earth3 + IPSL, all experiments present) ----------
tmap = {"EC-Earth3":"EC-Earth-Consortium_EC-Earth3", "IPSL-CM6A-LR":"IPSL_IPSL-CM6A-LR"}
for exp in ("historical","piControl","ssp126","ssp585"):
    for model,tag in tmap.items():
        f = f"{BASE}/terhaar_amoc/amoc/amoc/26.5N/{exp}/amoc_{exp}_{tag}_r1i1p1f1.nc"
        if not os.path.exists(f): continue
        ds = xr.open_dataset(f, decode_times=False)
        add("Terhaar", model, exp, ds, "amoc")

# order & write
order = {"M26":0,"LeBars":1,"Terhaar":2}
rows.sort(key=lambda r:(order[r["fonte"]], r["modello"], r["esperimento"], str(r["variabile"])))
cols = ["fonte","modello","esperimento","variabile","units","latitudine",
        "long_name_o_description","time_var","n_time","time_primo","time_ultimo",
        "time_units","n_valid","media","std","min","max"]
out = f"{BASE}/inventario_amoc.csv"
with open(out,"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
print("rows:", len(rows), "->", out)
# print compact
for r in rows:
    print("|".join(str(r[c]) for c in ["fonte","modello","esperimento","variabile","units","latitudine","time_units","n_time","n_valid","time_primo","time_ultimo","media","std","min","max"]))
