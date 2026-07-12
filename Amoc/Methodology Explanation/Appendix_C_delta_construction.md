# Appendix C — The AMOC Delta: Construction and Statistical Diagnosis

This appendix documents how the AMOC-weakening **delta** is constructed from the
hosing and control simulations of Appendix A, and the statistical criteria by
which it is judged a genuine forced response rather than internal variability.
The construction is a calendar-matched anomaly reduced to its settled equilibrium
level; the diagnosis rests on **effect sizes** — signal-to-noise ratios
calibrated against the control simulation — rather than on tail p-values, a
choice dictated by the short length of the control runs. Three diagnostics are
reported: a detection ratio, a spatial fingerprint, and a per-cell seasonal
significance field. Throughout, computations are on each model's native grid
(Appendix A.2), and the EC-Earth3 extremes are the reconstructed fields of
Appendix B.

## C.1 The anomaly field

For each model, hosing protocol, and variable, the monthly anomaly is the hosing
series minus the piControl **monthly climatology**, matched by calendar month.
The climatology is the control's mean seasonal cycle, obtained by averaging every
control year within each calendar month (`ymonmean`); the anomaly is then formed
month by month (`ymonsub`),

    anomaly(cell, t) = V_hosing(cell, t) − V̄_piControl(cell, m(t)),     m(t) = calendar month of t,

so the seasonal cycle is removed by construction and only the seasonal phase, not
the absolute model year, need be matched (Appendix A.3). Both hosing protocols
(`g01`, weak; `u03`, strong) are differenced against the *same* control baseline,
so their deltas are directly comparable. The subtraction is strictly cell by
cell: where the hosing and control grids of a model coincide in size their
coordinate labels are snapped exactly (`setgrid`, no interpolation), and only a
true grid-size mismatch triggers a genuine remap (`remapbil` for temperature,
`remapcon` for precipitation). For precipitation a multiplicative companion field
`R = hosing / climatology` (`ymondiv`) is also written, for the coupling step
(Appendix E); the additive anomaly above is the one every diagnostic below reads.

All regional reductions use **cosine-of-latitude area weights**, so each cell
contributes in proportion to the area it subtends, and are taken over the
European window lon [−15°, 40°], lat [34°, 72°].

## C.2 The delta as an equilibrium level

The hosing forcing is a **constant freshwater flux held fixed** for the length of
the run. A system driven by such a step responds by **ramp-then-plateau** — it
approaches a new equilibrium and levels off — so the quantity of interest is the
*level* of the settled response, not a trend. The delta is therefore estimated as
a **late-run mean**. Writing the Europe-mean monthly anomaly series as
*a₁ … aₙ*, and *pl* = max(12, ⌊n/3⌋) for the plateau window,

    plateau = mean( a_{n−pl+1}, …, a_n )               (last third of the run, ≥ 1 yr)

Averaging over the plateau window is variance reduction: it suppresses residual
internal variability in the estimate of the equilibrium shift. The **delta field**
is the same last-third mean taken per grid cell,

    Δ(cell) = mean of the last pl months of anomaly(cell, ·),

written to NetCDF per model × protocol × variable. A flat late-run anomaly is
thus the correct, expected outcome, not a null result; where a model has not yet
levelled off (HadGEM3-GC31-MM is still cooling at year 100) its delta is a
**lower bound** on the true equilibrium shift.

## C.3 Detection: an effect size against the control

Whether a plateau is a genuine response is judged against the **internal
variability of the control**, at the *same averaging timescale*. The control
Europe-mean series is first **deseasonalised** (its own calendar-month
climatology removed) and **de-drifted** (a linear trend removed, the mean
restored), because the anomaly it is compared with already carries neither a
seasonal cycle nor spin-up drift; leaving either in would inflate the reference
spread and understate the effect. From the processed control series *c* the null
distribution of plateau-length means is formed over all overlapping windows of
length *pl*,

    σ_pl = sd( { mean(c_s, …, c_{s+pl−1}) : s = 1 … N−pl+1 } ),

and the delta is standardised by it,

    effect ratio = | plateau | / σ_pl,      effect declared when  ratio ≥ K,  with  K = 2.

This is a **signal-to-noise effect size**, not a significance probability. The
choice is deliberate: the HadGEM control is only about 100 years, so a tail
p-value would be governed by scarce degrees of freedom rather than by the
science; an effect size instead reports the magnitude of the shift relative to
natural variability, and is not a degrees-of-freedom problem. The K = 2 threshold
is the conventional two-standard-deviation bar. A run that sits inside the control
spread — EC-Earth3 `g01`, a 0.1 Sv forcing over 50 years, not yet settled —
correctly returns *no effect*.

## C.4 The AMOC fingerprint

An effect that clears C.3 could in principle still be a spatially uniform offset —
residual control drift — rather than the AMOC. The two are separated by **spatial
pattern**: AMOC weakening cools the north-west (the sub-polar "cold blob") far
more than the Mediterranean. A two-box contrast on the delta field measures this,

    fingerprint = mean_NW(Δ) − mean_Med(Δ),

with NW = lon [−15°, 5°], lat [50°, 62°] and Mediterranean = lon [0°, 25°],
lat [34°, 45°], each area-weighted. A genuine AMOC signal has the north-west more
negative than the Mediterranean, i.e. **fingerprint < 0**; a uniform drift gives
≈ 0. This is a minimal projection of the response onto its expected spatial
signature, and it also separates the models from one another.

## C.5 Per-cell seasonal significance

The detection logic of C.3 is repeated cell by cell, separately for winter (DJF)
and summer (JJA), to map *where* the cooling is robust. Each monthly field is
reduced to **per-year seasonal means** (the three season-months averaged within
each year), and the delta of a cell is the mean of its last-third seasonal years.
The per-cell noise is that cell's own control seasonal-mean series, **linearly
detrended**, reduced to the standard deviation of its plateau-length running
means. The ratio is a per-cell signal-to-noise field,

    snr(cell) = Δ_season(cell) / σ_season(cell),      robust when  | snr | ≥ 2,

and the summary statistic is the **fraction of European cells that are robust**.
Reported as an effect-size field, it avoids per-cell p-values and hence the
multiple-testing (field-significance / FDR) problem. The precision of σ_season is
control-length limited: EC-Earth3 (≈ 500-yr control) pins it well, whereas a
100-yr HadGEM control gives only about three independent windows per cell, so its
contour is indicative rather than exact.

## C.6 Outputs

The procedure yields, per model × protocol × variable: a delta NetCDF (C.2); a
row in a model-contrast table carrying the plateau, σ_pl, the effect ratio and
verdict (C.3), and the fingerprint gradient (C.4); and a seasonal NetCDF pair
(delta and snr) with the robust-cell fraction (C.5). These are accompanied by
diagnostic figures — the Europe-mean trajectory with its LOESS fit (Cleveland,
1979), settled plateau, and pooled control central-95 % band; the delta maps on a
shared, symmetric scale clipped at the 98th percentile of |Δ| (so that a few
domain-edge or cold-blob cells do not wash out the dominant signal); and the
seasonal robust-cell maps.

## C.7 Assumptions and limitations

The delta is an **equilibrium-level** estimate: for a model that has not settled
it is a lower bound (C.2). Detection is an **effect size**, chosen for honesty
about the short control; it reports magnitude relative to internal variability
rather than a calibrated false-positive rate. The control is de-drifted only
linearly; for EC-Earth3, whose contemporaneous parallel window is unavailable,
the whole quasi-stationary post-spin-up control (model years 2259–2759) is used
as the reference, so a residual long-term drift is not perfectly removed
(Appendix A.3). Finally, the diagnostics read the additive anomaly for every
variable; the multiplicative precipitation field enters only the coupling of
Appendix E.

---

**References.**
Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
scatterplots. *Journal of the American Statistical Association*, 74(368),
829–836.
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975–1995.
