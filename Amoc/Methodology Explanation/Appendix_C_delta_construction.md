# Appendix C — The AMOC Delta: Construction and Statistical Diagnosis

This appendix documents how the supplied R scripts reduce monthly hosing-minus-
piControl anomaly files to late-run delta fields and diagnostic effect-size
summaries. The diagnostics are computed on each model's own grid after cropping
to Europe. EC-Earth3 `tasmin` and `tasmax` anomalies are based on the
reconstructed monthly fields described in Appendix B.

The code reports effect-size diagnostics rather than calibrated p-values. It
labels a run as an effect when the late-run Europe-mean anomaly is large relative
to control variability at the same averaging length; it does not prove physical
equilibration or causality by itself.

## C.1 The anomaly field

For each model, hosing protocol, and variable, the monthly additive anomaly is
the hosing series minus a piControl monthly climatology matched by calendar
month:

```text
anomaly(cell, t) = V_hosing(cell, t) - climatology_piControl(cell, month(t)).
```

The climatology is built with CDO `ymonmean`, and the subtraction is performed
with `ymonsub`. Within each model, both hosing protocols are compared against the
same piControl reference used by that model's anomaly builder. The EC-Earth3 and
HadGEM builders use different reference sources, as described in Appendix A.

Grid handling is not identical for all models. The EC-Earth3 anomaly builder
snaps climatology grid labels with `setgrid` when hosing and climatology grid
sizes match, and remaps the climatology when grid sizes differ (`remapbil` for
temperature variables and `remapcon` for precipitation). The HadGEM anomaly
builder does not perform this grid-size check or remap; it assumes that hosing
and piControl files are already compatible and applies `ymonsub` directly.

For precipitation, the anomaly builders also write a multiplicative companion
field:

```text
R = hosing / piControl climatology.
```

The diagnostics in this appendix read the additive anomaly files. They do not
read the precipitation ratio files.

All Europe-mean and box-mean reductions use cosine-of-latitude area weights.
The shared R reader crops lon `[-15, 40]`, lat `[34, 72]` and converts `pr` to
`mm/day` by multiplying by `86400`.

## C.2 The delta as a late-run mean

The delta diagnostic is computed as a late-run mean, not as a fitted trend. For a
Europe-mean monthly anomaly series `a_1, ..., a_n`, the plateau-window length is

```text
pl = max(12, round(n / 3)).
```

The Europe-mean plateau is

```text
plateau = mean(a_{n-pl+1}, ..., a_n).
```

The delta field is the same final-window mean computed independently at each
grid cell:

```text
Delta(cell) = mean of the last pl monthly values of anomaly(cell, .).
```

One NetCDF delta file is written per processed model, protocol, and variable.
Averaging over the final window reduces month-to-month variability in the delta
estimate. The script does not test whether the run has equilibrated. If a model
is still evolving at the end of the available run, that interpretation must come
from an external trend or settling diagnostic; the code itself only reports the
late-run window mean.

## C.3 Detection: an effect size against the control

The effect ratio compares the late-run Europe-mean plateau with internal
variability estimated from the corresponding model and variable's piControl
series. The control Europe-mean series is processed in the shared R helper:

1. read all matched piControl files for the model and variable;
2. compute the cosine-weighted Europe-mean monthly series;
3. remove that control series' own 12-month climatology;
4. remove a linear trend and restore the mean.

The null distribution is the set of all overlapping control means with the same
length `pl` as the hosing late-run window:

```text
sigma_pl = sd({ mean(c_s, ..., c_{s+pl-1}) : s = 1, ..., N-pl+1 }).
```

The effect ratio is

```text
effect_ratio = abs(plateau) / sigma_pl.
```

The script labels the effect as `YES` when the ratio is finite and at least
`K = 2`; otherwise it labels it `no`. This is an effect-size threshold, not a
p-value or calibrated false-positive probability. It avoids assigning tail
probabilities from short control runs, but the denominator still depends on the
length, variability, and drift treatment of the available control data.

## C.4 The AMOC fingerprint

The spatial fingerprint diagnostic is a two-box contrast computed on the delta
field:

```text
NW_minus_Med = mean_NW(Delta) - mean_Med(Delta).
```

The boxes are:

- NW: lon `[-15, 5]`, lat `[50, 62]`;
- Mediterranean: lon `[0, 25]`, lat `[34, 45]`.

Both box means are cosine-of-latitude weighted. Negative values are interpreted
as stronger cooling in north-western Europe than in the Mediterranean under the
assumed AMOC-cooling fingerprint. The code computes this contrast; it does not
by itself prove that a pattern is caused by AMOC weakening.

## C.5 Per-cell seasonal effect fields

The seasonal script computes DJF and JJA fields separately for each processed
model, protocol, and variable. Monthly fields are reshaped into Jan-Dec blocks,
and seasonal means are computed by averaging selected months within each block:

```text
DJF = months 12, 1, and 2 within the same Jan-Dec block
JJA = months 6, 7, and 8
```

Thus DJF is implemented as a same-block month selection; the code does not shift
December into the following meteorological winter.

For each anomaly cell and season, the number of seasonal years in the late-run
window is

```text
plY = max(5, round(number_of_seasonal_years / 3)).
```

The seasonal delta is the mean of the final `plY` seasonal years. The
corresponding control seasonal series is computed cell by cell, linearly
detrended, converted to all overlapping running means of length `plY`, and
reduced to a standard deviation. The per-cell signal-to-noise field is

```text
snr(cell) = seasonal_delta(cell) / seasonal_control_running_mean_sd(cell).
```

Cells with `abs(snr) >= 2` are treated as robust in the seasonal maps. The
reported robust-cell fraction is the unweighted fraction of finite European grid
cells satisfying this threshold. The script reports an effect-size field; it does
not compute per-cell p-values and does not apply a field-significance or FDR
correction.

## C.6 Outputs

The delta/effect script writes:

- one delta NetCDF per processed model, protocol, and variable;
- `amoc_effect_model_contrast.csv`, with model, protocol, variable, run length,
  plateau, control-window standard deviation, effect ratio, effect label,
  `NW_minus_Med`, and the delta filename.

The seasonal script writes:

- one seasonal NetCDF per processed model, protocol, variable, and season; each
  file contains both the seasonal delta variable and its corresponding
  `<variable>_snr` field;
- `seasonal_robust_summary.csv`;
- `seasonal_fields_significance.pdf`.

The plotting scripts also produce a clipped delta-field PDF and a time-series
PDF. The time-series script can add country pages and a `country_plateaus.csv`
file when the R `maps` package is available; otherwise those country outputs are
skipped and the Europe-mean pages are still produced.

## C.7 Assumptions and limitations

The delta is a late-run mean over the available hosing anomaly period. It should
not be described as a proven equilibrium response unless an additional settling
diagnostic supports that interpretation. The effect ratio is an effect-size
diagnostic and should not be interpreted as a p-value. Control noise estimates
use deseasonalisation and linear de-drifting for Europe-mean diagnostics; any
nonlinear residual drift is not explicitly modeled.

The seasonal robust maps are thresholded signal-to-noise maps, not multiple-
testing-corrected significance maps. The robust-cell percentage is unweighted by
area. The diagnostics use additive anomalies for all variables; the
multiplicative precipitation ratio is written for later coupling workflows but
is not read by the delta, effect-ratio, fingerprint, or seasonal diagnostic
scripts described here.

---

**References.**
Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
scatterplots. *Journal of the American Statistical Association*, 74(368),
829-836.
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975-1995.
