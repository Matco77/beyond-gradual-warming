# Appendix B — Reconstruction of Monthly Temperature Extremes for EC-Earth3

This appendix documents the EC-Earth3 reconstruction step used by the supplied
pipeline. In the EC-Earth3 anomaly builder, hosing `tas` and `pr` are read from
native EC-Earth3 hosing files, while hosing `tasmin` and `tasmax` are expected as
reconstructed absolute monthly fields. The reconstruction script creates those
absolute `tasmin` and `tasmax` files from hosing `tas` using relationships fitted
on EC-Earth3 piControl data. The HadGEM3 configurations do not use this
reconstruction script in the supplied workflow.

## B.1 Method

The script opens EC-Earth3 piControl `tas`, `tasmax`, and `tasmin` files matched
by configured filename globs, aligns the three series on their common time steps,
and fits independent ordinary least-squares relationships for each grid cell and
calendar month:

```text
tasmax = a_max(cell, month) + b_max(cell, month) * tas
tasmin = a_min(cell, month) + b_min(cell, month) * tas
```

The coefficients are computed from piControl monthly statistics:

```text
b = Cov(tas, extreme) / Var(tas)
a = E(extreme) - b * E(tas)
```

where the expectations are taken over all piControl time steps belonging to the
same calendar month. The fitted coefficients are then selected by the hosing
`time.month` coordinate and applied to hosing `tas`.

The fit is independent for every grid cell and calendar month. The code does not
apply spatial smoothing, spatial regularisation, or neighbour constraints. Each
slope measures the fitted sensitivity of the corresponding extreme (`tasmax` or
`tasmin`) to monthly mean temperature. The implied response of the diurnal range
would depend on the difference between the fitted `tasmax` and `tasmin`
sensitivities; the script does not directly fit diurnal temperature range.

Because the regression is matched by calendar month, piControl and hosing model
years do not need to overlap. The actual hosing period is whatever is contained
in the files matched by the configured hosing globs.

## B.2 Numerical treatment

The script guards cells with zero or near-zero piControl `tas` variance. Where
`abs(Var(tas)) <= 1e-12`, the slope is set to zero to avoid division by zero or
near-zero values. In that case, the unconstrained fitted value is the
calendar-month climatological extreme, because `a = E(extreme)` when `b = 0`.
The code does not report how many cells are affected by this safeguard.

After reconstruction, the script enforces physical ordering point by point:

```text
tasmax = max(reconstructed tasmax, hosing tas)
tasmin = min(reconstructed tasmin, hosing tas)
```

This guarantees `tasmax >= tas` and `tasmin <= tas` in the written fields. The
script does not report how often this ordering clamp changes values or whether
changes occur over Europe.

Grid handling is conditional. If the EC-Earth3 piControl and hosing latitude and
longitude coordinates match within `numpy.allclose`, the fitted coefficient
coordinates are snapped to the hosing coordinates and no interpolation is used.
If the grids differ, the fitted coefficients are interpolated to the hosing grid
with nearest-neighbour interpolation.

## B.3 Outputs and validation status

For each processed hosing experiment, the script writes two absolute monthly
NetCDF files beside the EC-Earth3 hosing `tas` files:

```text
tasmax_Amon_EC-Earth3_<label>_reconstructed.nc
tasmin_Amon_EC-Earth3_<label>_reconstructed.nc
```

The script does not write anomaly files. The EC-Earth3 CDO anomaly builder later
uses these reconstructed absolute fields as the hosing-side inputs and subtracts
the EC-Earth3 piControl monthly climatology with `ymonsub`.

No standalone validation script, synthetic test output, global range check, or
missing-value report is included in the supplied files. Claims about exact
synthetic recovery, physically plausible global ranges, or absence of missing
values must therefore be supported by a separate validation log before being
reported as results.

## B.4 Assumptions and limitations

The reconstruction assumes that the piControl-calibrated monthly relationship
between mean temperature and each extreme remains applicable under hosing. This
stationarity assumption is recorded in the metadata of the written NetCDF files.
The method is a monthly statistical reconstruction: it captures the fitted
response of monthly `tasmin` and `tasmax` to hosing `tas`, but it does not model
changes in sub-monthly temperature-distribution shape or daily variability.

The reconstructed fields should therefore be interpreted as first-order monthly
inputs for the anomaly construction step, not as a physical simulation of daily
temperature extremes.

---

**References.**
Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological forcing
dataset for land surface- and hydrological-models. WATCH Technical Report.
