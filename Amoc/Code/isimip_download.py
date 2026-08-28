#!/usr/bin/env python3
"""Download ISIMIP3b daily Europe cutouts for the Phase-3 branch.

Author: Marco Bova

Run with the project venv, which has isimip-client:
    ./.amoc_venv/bin/python Amoc/Code/isimip_download.py

WHY THE SERVER-SIDE CUTOUT IS MANDATORY, not an optimisation. The files we need are
2 models x 2 experiments x 4 variables x 3-4 decades = 56 global files, 105.6 GB. Measured
throughput from files.isimip.org is ~0.5-0.7 MB/s, so the globals would take ~49 h against
~3 h for the Europe cutouts (~7.6 GB, the box is 7% of the global grid). And subsetting them
ourselves is not possible: the variables are stored gzip-compressed with chunk shape
(1, 360, 720) - one chunk is one whole global day - so reading any European cell requires
fetching and decompressing the entire global slice for that day. Byte-range requests, which the
server does support, buy nothing. Only a process sitting on the file can cut it.

API VERSION: the files API is v2 (https://files.isimip.org/api/v2). The v1 endpoint still
answers, but with a maintenance page (HTTP 502), which is easy to mistake for an outage.
isimip_client.ISIMIPClient defaults to v2, so let it own the URL rather than hardcoding one.

Resumable: a (model, scenario, variable) group whose cutouts are already in OUT_DIR is skipped,
so the script can be killed and restarted without losing work.
"""
import os, re, sys, time, urllib.parse, urllib.request, json
from isimip_client.client import ISIMIPClient

ROOT    = os.path.expanduser("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi")
OUT_DIR = os.path.join(ROOT, "Amoc/datasets/Isimip3b")
DATA_API = "https://data.isimip.org/api/v1/files/"

MODELS = ["ipsl-cm6a-lr", "ec-earth3"]
VARS   = ["tas", "tasmin", "tasmax", "pr"]
# reference vs future window. Files come in decade chunks, so historical pulls 1981-2020 and
# ssp126 pulls 2071-2100; the exact windows are applied when the monthly delta is built.
WINDOWS = {"historical": (1985, 2014), "ssp126": (2071, 2100)}
# Wider than the ISIMIP suggestion (-15..40, 34..72), which covers only 95.0% of the E-OBS cells
# the two branches use. This box covers 100% (Iceland included), so the ISIMIP and AMOC branches
# run on the identical geography and the comparison carries no procedural difference.
WEST, EAST, SOUTH, NORTH = -45.0, 46.0, 25.0, 76.0

POLL = 15


def log(m):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}", flush=True)


def wanted():
    out = []
    for fm in MODELS:
        for sc, (y0, y1) in WINDOWS.items():
            for v in VARS:
                q = urllib.parse.urlencode({
                    "simulation_round": "ISIMIP3b", "climate_forcing": fm,
                    "climate_scenario": sc, "climate_variable": v,
                    "time_step": "daily", "page_size": 1000})
                with urllib.request.urlopen(DATA_API + "?" + q, timeout=90) as r:
                    js = json.load(r)
                for f in js["results"]:
                    m = re.search(r"_(\d{4})_(\d{4})\.nc$", f["name"])
                    if not m:
                        continue
                    a, b = int(m.group(1)), int(m.group(2))
                    if b < y0 or a > y1:
                        continue
                    out.append({"model": fm, "scen": sc, "var": v,
                                "name": f["name"], "path": f["path"], "size": f["size"]})
    return out


def already(fm, sc, v, n_expected):
    """A group is done when OUT_DIR holds n_expected cutouts matching model+scenario+variable."""
    if not os.path.isdir(OUT_DIR):
        return False
    pat = re.compile(rf"^{re.escape(fm)}_.*_{re.escape(sc)}_{re.escape(v)}_.*\.nc$")
    return sum(bool(pat.match(f)) for f in os.listdir(OUT_DIR)) >= n_expected


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    files = wanted()
    tot = sum(f["size"] for f in files)
    log(f"{len(files)} global files | {tot/1e9:.1f} GB global | ~{tot*0.072/1e9:.1f} GB after cutout")

    groups = {}
    for f in files:
        groups.setdefault((f["model"], f["scen"], f["var"]), []).append(f)

    client = ISIMIPClient()
    todo = [(k, v) for k, v in sorted(groups.items()) if not already(*k, len(v))]
    log(f"{len(todo)} of {len(groups)} cutout jobs to run "
        f"({len(groups)-len(todo)} already present)")

    failed = []
    for i, ((fm, sc, v), fs) in enumerate(todo, 1):
        label = f"{fm}_{sc}_{v}"
        try:
            t0 = time.time()
            r = client.cutout_bbox([f["path"] for f in fs],
                                   west=WEST, east=EAST, south=SOUTH, north=NORTH)
            log(f"[{i}/{len(todo)}] {label}: {len(fs)} files, job {r.get('id','?')[:8]} {r.get('status')}")
            while r.get("status") in ("queued", "started"):
                time.sleep(POLL)
                r = client.get_job(r["job_url"])
            if r.get("status") != "finished":
                log(f"    FAILED status={r.get('status')}")
                failed.append(label); continue
            client.download(r["file_url"], path=OUT_DIR, extract=True)
            log(f"    done in {time.time()-t0:.0f}s")
        except Exception as e:
            log(f"    ERROR {e!r}")
            failed.append(label)

    have = len([f for f in os.listdir(OUT_DIR) if f.endswith(".nc")])
    log(f"finished | .nc in {OUT_DIR}: {have} | failed groups: {len(failed)}"
        + (f" -> {failed}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
