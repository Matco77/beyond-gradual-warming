#!/usr/bin/env python3
"""Download ISIMIP3b daily Europe cutouts for the Phase-3 branch.

Author: Marco Bova

WHY A WAITER AND NOT A PLAIN DOWNLOAD. The global files are 105.6 GB for what we need
(2 models x 2 experiments x 4 variables x 3-4 decades). Measured throughput to
files.isimip.org is ~0.5-0.7 MB/s, so pulling the globals would take ~49 h. The Europe box is
~7% of the global grid, which brings it to ~7.6 GB and ~3.5 h - so the server-side cutout is
not an optimisation here, it is the only viable route. At the time of writing the cutout API
(https://files.isimip.org/api/v1) returns HTTP 502 while the plain file server is up, so this
script waits for the service to come back rather than falling back to the 49-hour path.

Resumable: a cutout already present in OUT_DIR is skipped, so it can be killed and restarted.
Stop it with:  kill $(cat <scratch>/isimip_dl.pid)
"""
import json, os, re, sys, time, urllib.error, urllib.parse, urllib.request, zipfile

ROOT     = os.path.expanduser("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi")
OUT_DIR  = os.path.join(ROOT, "Amoc/datasets/Isimip3b")
DATA_API = "https://data.isimip.org/api/v1/files/"
PROC_API = "https://files.isimip.org/api/v1"

MODELS = ["ipsl-cm6a-lr", "ec-earth3"]
VARS   = ["tas", "tasmin", "tasmax", "pr"]
# reference window vs future window. The files come in decade chunks, so historical pulls
# 1981-2020 and ssp126 pulls 2071-2100; the exact windows are applied when the delta is built.
WINDOWS = {"historical": (1985, 2014), "ssp126": (2071, 2100)}
# bbox = [south, north, west, east]. Wider than the ISIMIP suggestion (-15..40, 34..72), which
# covers only 95.0% of the E-OBS cells the two branches use; this box covers 100% (Iceland
# included), so the ISIMIP and AMOC branches run on the identical geography.
BBOX = [25.0, 76.0, -45.0, 46.0]

POLL_SERVICE = 300     # seconds between availability probes while the API is down
POLL_JOB     = 20      # seconds between job-status polls
MAX_WAIT_H   = 72      # give up waiting for the service after this many hours


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def get_json(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)


def post_json(url, payload, timeout=180):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def wanted_files():
    """Ask the data API which global files cover the windows we need."""
    out = []
    for model in MODELS:
        for scen, (y0, y1) in WINDOWS.items():
            for var in VARS:
                q = urllib.parse.urlencode({
                    "simulation_round": "ISIMIP3b", "climate_forcing": model,
                    "climate_scenario": scen, "climate_variable": var,
                    "time_step": "daily", "page_size": 1000})
                for f in get_json(DATA_API + "?" + q)["results"]:
                    m = re.search(r"_(\d{4})_(\d{4})\.nc$", f["name"])
                    if not m:
                        continue
                    a, b = int(m.group(1)), int(m.group(2))
                    if b < y0 or a > y1:            # decade chunk outside the window
                        continue
                    out.append({"model": model, "scen": scen, "var": var,
                                "name": f["name"], "path": f["path"], "size": f["size"]})
    return out


def service_up():
    try:
        urllib.request.urlopen(PROC_API, timeout=30)
        return True
    except urllib.error.HTTPError as e:
        return e.code < 500                        # 4xx means it is answering, just not to a GET
    except Exception:
        return False


def wait_for_service():
    t0 = time.time()
    while not service_up():
        waited = (time.time() - t0) / 3600
        if waited > MAX_WAIT_H:
            log(f"cutout API still down after {MAX_WAIT_H} h - giving up")
            return False
        log(f"cutout API down (waited {waited:.1f} h), retrying in {POLL_SERVICE//60} min")
        time.sleep(POLL_SERVICE)
    log("cutout API is answering")
    return True


def run_job(paths, label):
    """Submit one cutout job and block until the zip is downloaded and extracted."""
    job = post_json(PROC_API, {"task": "cutout_bbox", "paths": paths, "bbox": BBOX})
    url = job.get("job_url") or job.get("url")
    log(f"  {label}: job {job.get('id')} status={job.get('status')}")
    while job.get("status") in ("queued", "started", "pending", None):
        time.sleep(POLL_JOB)
        if not url:
            log(f"  {label}: no job_url in response, raw={json.dumps(job)[:300]}")
            return False
        job = get_json(url)
    if job.get("status") != "finished":
        log(f"  {label}: job ended status={job.get('status')} raw={json.dumps(job)[:300]}")
        return False
    zurl = job.get("file_url")
    zpath = os.path.join(OUT_DIR, f"_{label}.zip")
    log(f"  {label}: downloading {zurl}")
    urllib.request.urlretrieve(zurl, zpath)
    with zipfile.ZipFile(zpath) as z:
        z.extractall(OUT_DIR)
    os.remove(zpath)
    log(f"  {label}: done")
    return True


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    files = wanted_files()
    tot = sum(f["size"] for f in files)
    log(f"{len(files)} global files cover the windows | global volume {tot/1e9:.1f} GB "
        f"| Europe cutout ~{tot*0.072/1e9:.1f} GB")

    # one job per (model, scenario, variable): 3-4 decade files each, 16 jobs in total
    groups = {}
    for f in files:
        groups.setdefault((f["model"], f["scen"], f["var"]), []).append(f)

    todo = []
    for key, fs in sorted(groups.items()):
        label = "_".join(key)
        # already extracted? the cutout keeps the stem and appends the bbox, so match on the stem
        stems = [os.path.splitext(f["name"])[0].replace("_global_", "_") for f in fs]
        have = os.listdir(OUT_DIR) if os.path.isdir(OUT_DIR) else []
        if all(any(s.split("_daily_")[-1] in h and key[2] in h and key[1] in h for h in have)
               for s in stems) and have:
            log(f"  {label}: already present, skipped")
            continue
        todo.append((label, [f["path"] for f in fs]))
    log(f"{len(todo)} cutout jobs to run")

    if not todo:
        log("nothing to do")
        return 0
    if not wait_for_service():
        return 1

    failed = []
    for label, paths in todo:
        try:
            if not run_job(paths, label):
                failed.append(label)
        except Exception as e:                     # one bad job must not kill the whole run
            log(f"  {label}: ERROR {e!r}")
            failed.append(label)
    log(f"finished | ok {len(todo)-len(failed)} | failed {len(failed)}"
        + (f" -> {failed}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
