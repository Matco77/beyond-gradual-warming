# AMOC-weakening pipeline (HadGEM + EC-Earth3)

Two steps: build the hosing-minus-piControl **anomalies** (CDO), then compute
**effect size / fingerprint / delta fields** and the model-contrast table (R).

Both steps now cover **two NAHosMIP hosing protocols** for the HadGEM models:

| protocol  | meaning      |
|-----------|--------------|
| `g01-hos` | weak hosing  |
| `u03-hos` | strong hosing (added) |

for all four variables `tas`, `tasmin`, `tasmax`, `pr`, on
`HadGEM3-GC31-LL` and `HadGEM3-GC31-MM`. (EC-Earth3 already had both protocols.)

## 1. Build anomalies — `nahos_anomaly_cdo.sh`

```bash
bash scripts/nahos_anomaly_cdo.sh        # needs CDO:  brew install cdo
```

For every `model × protocol × variable` it:

1. `mergetime` of **all** hosing monthly chunk files — auto-discovered by glob,
   so the (different) year-ranges of the `u03-hos` runs need not be hard-coded;
2. `mergetime` of the piControl chunk(s) overlapping the **1850–1949** baseline
   window — auto-selected from the chunk filename date-token, so it works for
   LL (one 100-yr file) and MM (five 20-yr files) alike;
3. `monmean` (a no-op guard on already-monthly `Amon` data) → `ymonmean`
   monthly climatology;
4. `ymonsub`: anomaly = hosing − piControl climatology.

A `model × protocol` with no hosing data is **skipped with a warning** (so it is
safe to run even before a `u03-hos` run has been downloaded). Outputs land in:

```
anomaly_output/LL_anomaly/<var>_anomaly_<proto>_minus_piControl_1850-1949.nc
anomaly_output/MM_anomaly/<var>_anomaly_<proto>_minus_piControl_1850-1949.nc
```

For **`pr` only**, a second, *multiplicative* field is also written — `R = hosing /
piControl climatology` — the quantity the ISIMIP stressing step applies as
`pr_AMOC = pr_ISIMIP · R`. Temperature stays additive (a ratio of interval-scale K
is meaningless), and the additive `pr` anomaly above is unchanged, so the effect /
delta / fingerprint diagnostics are untouched:

```
anomaly_output/{LL,MM}_anomaly/pr_ratio_<proto>_over_piControl_1850-1949.nc
anomaly_output/ECHearth3_anomaly/pr_Amon_EC-Earth3_hos-<proto>-hos_ratio.nc
```

(EC-Earth3 is built by `1.ECHEarth3nahos_anomaly_cdo_explicit.sh`, which divides by
its whole-run 2259–2759 climatology; the HadGEM builder divides by the 1850–1949
baseline. Both write `R = hosing / piControl climatology` for `pr` only.)

The baseline is protocol-independent, so `g01-hos` and `u03-hos` anomalies share
one common piControl reference and stay directly comparable. The `g01-hos`
filenames are reproduced exactly as before (backward compatible).

> Edit the `DATASETS` path (and, if needed, `PROTOCOLS=(u03-hos)` to rebuild only
> the strong-hosing anomalies) at the top of the script.

## 2. Effect / fingerprint / delta — `amoc_delta_effect_contrast.R`

```bash
Rscript scripts/amoc_delta_effect_contrast.R     # needs: ncdf4, fields
```

Consumes the 24 anomaly files (EC-Earth3, HadGEM-LL, HadGEM-MM × `g01`+`u03`)
and, per model × protocol × variable, writes:

- a **delta** NetCDF (late-run / plateau-window mean anomaly over Europe),
- a row in `delta_fields/amoc_effect_model_contrast.csv` with the plateau,
  the control-window effect ratio (effect = YES when `|plateau| ≥ K·sd`), and
  the NW-minus-Mediterranean **fingerprint** gradient.

Protocol (`u03` vs `g01`) and model are inferred from each anomaly's filename /
folder, so the new HadGEM `u03-hos` files are picked up automatically. Any
anomaly not built yet is skipped with a notice instead of aborting the run.

## Notes

- **OneDrive resilience.** Both scripts build every NetCDF/CSV in a *local*
  scratch dir and then publish it into the `CloudStorage/OneDrive…` tree with a
  short retry. OneDrive's macOS File Provider briefly locks files mid-sync, which
  otherwise makes a direct write fail with `cdi error (cdf__create): … Permission
  denied`. If publishing still can't get the lock after the retries, the file is
  kept in scratch and the run continues (it no longer aborts the whole batch).
  Override the scratch location with `AMOC_SCRATCH=/some/local/dir`. For a fully
  clean run you can also pause OneDrive syncing first.
- **All available hosing years are used.** The hosing glob merges *every* chunk
  present for a `model × protocol × variable`. For `HadGEM3-GC31-MM g01-hos` that
  is the full ~100-yr run (five 20-yr chunks), not the 40 yr the old explicit
  script hard-coded — which matches the R script's "MM still cooling at year 100"
  note. To restrict to a shorter window, narrow the glob or move the unwanted
  chunks out of the `NAHosMIP/<model>/` folder.
