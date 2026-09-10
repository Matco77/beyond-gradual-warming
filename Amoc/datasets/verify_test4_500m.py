import numpy as np, xarray as xr, glob, os

def atlantic_idx(ds, sectorvar, basindim):
    sec = ds[sectorvar].values
    labels = []
    for row in sec:
        if isinstance(row, bytes):
            labels.append(row.decode(errors="replace").strip())
        else:
            labels.append(b"".join(x for x in row if isinstance(x, bytes)).decode(errors="replace").strip()
                          if row.dtype.kind == "S" else "".join(row.astype(str)).strip())
    return labels

# ---------- EC-Earth3 historical (byte-string sector) ----------
base = "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
files = sorted(glob.glob(f"{base}/nahosmip_m26/ecearth3_msftyz_hist/*.nc"))
d0 = xr.open_dataset(files[0])
sec = d0["sector"].values
labels = ["".join(c.decode() if isinstance(c, bytes) else c for c in row).strip()
          if sec.dtype.kind in ("S","U") and sec.ndim==2 else str(row) for row in sec]
# robust decode for char(basin,strlen)
raw = xr.open_dataset(files[0], decode_cf=False)["sector"].values
labels = [b"".join(raw[i]).decode(errors="replace").strip() if raw.dtype.kind=="S"
          else "".join(raw[i].astype(str)).strip() for i in range(raw.shape[0])]
print("EC-Earth3 sector labels (raw byte-strings decoded):", labels)
atl = [i for i,l in enumerate(labels) if "atlantic" in l.lower()][0]
rlat = d0["rlat"].values
rj = int(np.abs(rlat-26.5).argmin())
lev = d0["lev"].values
print(f"Atlantic basin idx={atl} ({labels[atl]}) | rlat nearest 26.5 = {rlat[rj]:.4f} | levels {lev.min():.1f}-{lev.max():.0f}m")

full=[]; deep=[]; dfull=[]; ddeep=[]
for f in files:
    ds=xr.open_dataset(f)
    v=ds["msftyz"].isel(basin=atl, rlat=rj)          # (time, lev)
    v=v.where(v!=v.attrs.get("_FillValue",1e20))
    prof=v.mean("time").values/1026/1e6              # annual mean profile -> Sv
    full.append(np.nanmax(prof)); dfull.append(lev[np.nanargmax(prof)])
    mask=lev>500
    deep.append(np.nanmax(prof[mask])); ddeep.append(lev[mask][np.nanargmax(prof[mask])])
full=np.array(full); deep=np.array(deep)
diff=full-deep
print("\n=== TEST 4a EC-Earth3 historical: full-column max vs max>500m ===")
print(f"years sampled: {len(files)}")
print(f"full-max  mean={full.mean():.3f} Sv | deep-max mean={deep.mean():.3f} Sv")
print(f"diff (full-deep): mean={diff.mean():+.4f} Sv | max={np.abs(diff).max():.4f} Sv")
print(f"depth of max  full: {np.unique(np.round(dfull)).astype(int)} m")
print(f"depth of max  >500: {np.unique(np.round(ddeep)).astype(int)} m")

# ---------- IPSL ssp126 control check (numeric 3basin) ----------
ip=xr.open_dataset(f"{base}/msftyz_ipsl/msftyz_Omon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gn_201501-210012.nc")
b=ip["3basin"].values; iatl=int(np.where(b==2)[0][0])
nav=ip["nav_lat"].isel(x=0).values; ij=int(np.abs(nav-26.5).argmin())
ol=ip["olevel"].values
col=ip["msftyz"].isel({"3basin":iatl,"y":ij,"x":0}); col=col.where(col!=1e20)
ann=col.groupby(col.time.dt.year).mean("time").values/1026/1e6   # (year, olevel)
fmax=np.nanmax(ann,axis=1); dmask=ol>500
dmax=np.nanmax(ann[:,dmask],axis=1)
dd=fmax-dmax
print("\n=== TEST 4b IPSL ssp126 (control): full vs >500m ===")
print(f"diff (full-deep): mean={dd.mean():+.4f} Sv | max={np.abs(dd).max():.4f} Sv")
print(f"depth of max full: {sorted(set(np.round(ol[np.nanargmax(ann,axis=1)]).astype(int)))} m")
