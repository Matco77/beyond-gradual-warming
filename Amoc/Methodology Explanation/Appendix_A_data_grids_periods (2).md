# Appendix A — Data, Spatial Grids, and Reference Periods

This appendix documents the climate-model inputs used by the AMOC-weakening
anomaly pipeline, the grids on which the scripts operate, and the reference
periods used to build the monthly hosing-minus-piControl climatologies. The
statements below describe the behavior of the supplied scripts; where a claim
would require external metadata or a completed run log, it is identified as such.

## A.1 Input data

The workflow uses monthly (`Amon`) NetCDF fields for near-surface air temperature
(`tas`), daily-minimum and daily-maximum near-surface air temperature expressed
as monthly fields (`tasmin`, `tasmax`), and precipitation (`pr`). The scripts are
configured for three source models:

- EC-Earth3;
- HadGEM3-GC31-LL;
- HadGEM3-GC31-MM.

Two NAHosMIP hosing labels are processed where matching files exist:
`g01-hos` and `u03-hos`. For the HadGEM3 models, the anomaly builder loops over
both hosing labels and all four variables directly. For EC-Earth3, the anomaly
builder uses native hosing files for `tas` and `pr`, and expects reconstructed
absolute hosing files for `tasmin` and `tasmax`. Those reconstructed files are
produced by the EC-Earth3 reconstruction script, which fits per-cell,
per-calendar-month regressions of each temperature extreme on `tas` using
EC-Earth3 piControl data.

The code assumes the local directory structure used in the scripts. Paths are
absolute and must be edited if the project is moved.

## A.2 Spatial grids

The downstream R diagnostics read each anomaly file on its own grid, wrap
longitudes to `[-180, 180]`, sort longitudes and latitudes as needed, crop the
European analysis window lon `[-15, 40]`, lat `[34, 72]`, and do not interpolate
during that read/crop step. Precipitation is converted from flux units to
`mm/day` in the shared R reader by multiplying by `86400`.

The grid-handling behavior differs by stage and model:

- **HadGEM anomaly construction:** the HadGEM shell script assumes that hosing
  and piControl files are grid-compatible and applies CDO `ymonsub` directly. It
  does not perform a grid-size check, `setgrid`, or remapping.
- **EC-Earth3 anomaly construction:** the EC-Earth3 shell script compares the
  hosing and climatology grid sizes. If the sizes match, it copies the hosing
  grid description onto the climatology with `setgrid`, which relabels the grid
  but does not interpolate. If the sizes differ, it remaps the climatology to the
  hosing grid: `remapbil` for temperature variables and `remapcon` for
  precipitation.
- **EC-Earth3 `tasmin`/`tasmax` reconstruction:** the Python script checks
  whether piControl and hosing latitude/longitude coordinates match within
  `numpy.allclose`. If they match, coefficient coordinates are snapped to the
  hosing coordinates. If they differ, the fitted coefficients are placed on the
  hosing grid using nearest-neighbour interpolation.

The following grid dimensions are descriptive source-data metadata, not values
hard-coded or independently verified by the supplied scripts:

| Source                 | Native grid label | Grid (lon × lat) | Approximate resolution |
|------------------------|-------------------|------------------|------------------------|
| EC-Earth3              | `gr`              | 512 × 256        | about 0.703° × 0.703°  |
| HadGEM3-GC31-LL        | `gn`              | 192 × 144        | about 1.875° × 1.250°  |
| HadGEM3-GC31-MM        | `gn`              | 432 × 324        | about 0.833° × 0.556°  |

Claims such as exact coordinate equality, floating-point differences of a
specific magnitude, or verified grid dimensions require a separate grid-check
output, for example from `cdo griddes`, `ncdump`, or an archived run log.

## A.3 Reference periods and rationale

For each variable and grid cell, the additive anomaly is computed as the hosing
monthly field minus a piControl monthly climatology matched by calendar month:

```text
anomaly(cell, t) = V_hosing(cell, t) - climatology_piControl(cell, month(t)).
```

The climatology is a twelve-month mean seasonal cycle built with CDO `ymonmean`,
and the anomaly is formed with `ymonsub`. The hosing and piControl years do not
need to have the same calendar years because the subtraction matches only the
seasonal phase.

The reference construction is model-specific:

| Model | Hosing files used by code | piControl reference used by code | Notes |
|---|---|---|---|
| EC-Earth3 | `hos-g01-hos` and `hos-u03-hos`; native `tas`/`pr` plus reconstructed `tasmin`/`tasmax` | all EC-Earth3 piControl files matched by the configured glob; script comments describe these as model years 2259-2759 | no filename date filtering is applied in the EC-Earth3 anomaly builder |
| HadGEM3-GC31-LL | all matching hosing chunks for `g01-hos` and `u03-hos` | piControl chunks whose filename date token overlaps `185001-194912` | overlapping chunks are selected, but the script does not crop a chunk internally if it extends beyond the target window |
| HadGEM3-GC31-MM | all matching hosing chunks for `g01-hos` and `u03-hos` | piControl chunks whose filename date token overlaps `185001-194912` | same selection logic as HadGEM3-GC31-LL |

The HadGEM builder therefore produces an exact 1850-1949 climatology only when
the selected source chunks align with that period. The EC-Earth3 builder uses the
full set of files matched by its piControl glob. The script comments justify this
with the absence of a contemporaneous EC-Earth3 branch-window file and describe
the available control as post-spin-up, but the supplied code does not itself test
quasi-stationarity or quantify residual drift.

## A.4 Output

The anomaly builders write monthly additive anomaly files for `tas`, `tasmin`,
`tasmax`, and `pr`. For `pr`, they also write a multiplicative ratio field:

```text
R(cell, t) = pr_hosing(cell, t) / climatology_pr_piControl(cell, month(t)).
```

The R diagnostics described in Appendix C read the additive anomaly files, not
the precipitation ratio files.

Coverage depends on the way each builder finds hosing files. The HadGEM builder
merges all hosing chunks matched by its glob for each model, protocol, and
variable. The EC-Earth3 reconstruction script opens all matching hosing `tas`
files for reconstruction, but the EC-Earth3 anomaly shell script selects the
first matching native `tas` and `pr` file with `ls ... | head -1`; therefore the
native EC-Earth3 `tas`/`pr` anomaly period is the period contained in that first
matched file.

The output anomaly fields are the climate-model AMOC-hosing deltas used by later
R diagnostics and, where implemented elsewhere, by later coupling or scenario
replay steps. The supplied anomaly and diagnostic scripts do not themselves
perform the full ISIMIP coupling.

---

**References.**
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975-1995.
Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological forcing
dataset for land surface- and hydrological-models. WATCH Technical Report.
