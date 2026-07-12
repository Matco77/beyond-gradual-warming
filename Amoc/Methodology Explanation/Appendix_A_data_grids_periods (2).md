# Appendix A — Data, Spatial Grids, and Reference Periods

This appendix documents the climate model data used to construct the stylised
AMOC-weakening signal, the spatial grids on which the inputs and outputs are
defined, and the time periods chosen for each experiment together with their
justification.

## A.1 Input data

The analysis draws on two classes of simulation from the CMIP6 archive, at
monthly resolution (`Amon`), for the variables near-surface air temperature
(`tas`), daily minimum and maximum temperature (`tasmin`, `tasmax`), and
precipitation (`pr`):

1. **piControl** — the unforced pre-industrial control simulation, used as the
   baseline climate state.
2. **NAHosMIP hosing experiments** (Jackson et al., 2023) — idealised
   freshwater-forcing ("hosing") runs that mechanically weaken the AMOC. Two
   protocols are used: `g01-hos` (primary) and `u03-hos` (secondary, longer).

Three source models are employed: EC-Earth3, HadGEM3-GC31-LL, and
HadGEM3-GC31-MM. In the hosing runs, EC-Earth3 provides only `tas` and `pr`;
its `tasmin` and `tasmax` are therefore **reconstructed** from `tas` using a
per-cell, per-calendar-month linear regression of each extreme on the mean,
fitted on piControl (after Weedon et al., 2010). Both HadGEM3 configurations
provide all four variables natively.

## A.2 Spatial grids

All computations are performed on each model's **native grid**; no spatial
interpolation is applied within a model, because the hosing and piControl
fields of a given model share the same grid. The output anomaly fields inherit
this native grid.

| Source                | Native grid label | Grid (lon × lat) | Resolution (lon × lat) |
|-----------------------|-------------------|------------------|------------------------|
| EC-Earth3 (T255)      | `gr`              | 512 × 256        | 0.703° × 0.703°        |
| HadGEM3-GC31-LL (N96) | `gn`              | 192 × 144        | 1.875° × 1.250°        |
| HadGEM3-GC31-MM (N216)| `gn`              | 432 × 324        | 0.833° × 0.556°        |

For EC-Earth3 the piControl and hosing grids were verified to be identical:
longitudes coincide exactly and latitudes differ only by floating-point noise
(≈10⁻¹³ degrees). Coordinate labels are reconciled by an exact relabelling
(`setgrid`) rather than interpolation, so the subtraction is performed
cell-by-cell with no smoothing. Genuine regridding (bilinear for temperature,
first-order conservative for precipitation) is reserved for the subsequent
coupling step, in which the anomaly is mapped onto the 0.5° ISIMIP3b grid.

## A.3 Reference periods and rationale

The anomaly is defined, for each grid cell and calendar month *m*, as

    anomaly(cell, t) = V_hosing(cell, t) − V̄_piControl(cell, m),

where V̄_piControl is the piControl monthly climatology (the mean seasonal
cycle). The subtraction is matched by calendar month, so the absolute model
years of the hosing and piControl segments need not coincide; only the seasonal
phase is matched. The full hosing time series is retained, while piControl is
reduced to a twelve-month climatology.

The piControl reference window is, in principle, the **contemporaneous parallel
window** — the control years from which the hosing run branched — because this
choice cancels any slow drift in the control simulation. The availability of
that window differs by model:

| Model            | Hosing run (period)        | piControl baseline period | Basis |
|------------------|----------------------------|---------------------------|-------|
| EC-Earth3        | g01-hos 1850–1899; u03-hos 1850–1949 | full piControl 2259–2759 | parallel window unavailable (see below) |
| HadGEM3-GC31-LL  | g01-hos 2050–2149          | 1850–1949                 | parallel branch window |
| HadGEM3-GC31-MM  | g01-hos 2050–2149          | 1850–1949                 | parallel branch window |

For the two HadGEM3 configurations the hosing runs branch from piControl years
1850–1949, which are present in the archive; these years are used directly as
the baseline, so control drift is cancelled.

For EC-Earth3 the run metadata indicate that the hosing experiment branches
from the **piControl-spinup** phase (`branch_time_in_parent = 0`), an epoch not
retained in the published piControl output (which spans model years 2259–2759).
No contemporaneous parallel window is therefore available. The **entire
published piControl (2259–2759)** is used instead as the unforced reference.
This is justified because the EC-Earth3 piControl is quasi-stationary after
spin-up, so its long-term mean seasonal cycle is a representative control state.
The single consequence is that any residual long-term drift in the control is
not perfectly removed; given the quasi-stationarity of the post-spin-up control
this effect is small, and a formal drift correction is left to future work.

## A.4 Output

For each model and hosing experiment the procedure yields monthly anomaly
fields for `tas`, `tasmin`, `tasmax`, and `pr`, on the model's native grid,
spanning the full hosing period (50 or 100 years). These fields constitute the
stylised AMOC-weakening signal that is subsequently superimposed on the
observational-era ISIMIP3b baseline to assess impacts on European energy demand
and production.

---

**References.**
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975–1995.
Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological forcing
dataset for land surface- and hydrological-models. WATCH Technical Report.
