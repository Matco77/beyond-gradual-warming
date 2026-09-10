# The ISIMIP Application Step — Implementation Proposal

*Operational specification of the delta coupling (Appendix E), with the data-driven
safeguards identified in the July 2026 audit of the anomaly pipeline.
Status: proposal — this step is specified here but not yet implemented in code.*

---

## 0. Purpose and positioning

Appendices A–C construct, for each climate model (EC-Earth3, HadGEM3-GC31-LL,
HadGEM3-GC31-MM) and each hosing protocol (g01, u03), monthly
hosing-minus-piControl anomaly fields on the model's native grid. Appendix E
states the *concept* for combining that signal with a bias-adjusted daily
baseline (ISIMIP3b, GFDL-ESM4 adjusted to W5E5, SSP3-7.0, 0.5°, Europe,
2015–2030): temperature is shifted additively, precipitation is scaled
multiplicatively, matched by calendar month.

This document is the missing layer between the two: the **step-by-step
operational recipe** — which epoch of the anomaly series is used, how the
precipitation ratio is actually constructed and sanitised, in what order
regridding and safeguards are applied, and how the result is validated. Each
choice is stated with its rationale and, where it was measured, the number that
motivates it. Nothing here changes the method of Appendix E; it makes it
executable and safe.

## 1. Inputs (all verified on disk, 2026-07-14)

| Input | Content | Role |
|---|---|---|
| `<var>_…_anomaly.nc` (×3 models ×2 protocols ×4 vars) | **full monthly series** of hosing − piControl-climatology (50 or 100 yr), native grid | source of the delta |
| piControl monthly climatologies (12 months, per model/var) | reference seasonal cycle (rebuilt in one CDO call where not cached) | denominator of the precipitation ratio |
| `pr_…_ratio.nc` | stored monthly ratio hosing/climatology, **unclipped** | diagnostics only — *not* applied directly (§4) |
| ISIMIP3b daily `tas, tasmin, tasmax, pr` (GFDL-ESM4 → W5E5, SSP3-7.0, 0.5°, Europe window, 2015–2030) | absolute daily baseline | the climate being perturbed |

Temperature anomalies are in K; precipitation fields in kg m⁻² s⁻¹ (the mm/day
conversion lives in the diagnostics, not in the files).

## 2. Design principles (recap, three lines)

1. The climate model contributes **only the change** (hosing − control); its
   absolute state and its biases never touch the impact models (Appendix E.1).
2. The baseline contributes the **absolute daily and within-day climate**; its
   variability is inherited unchanged — the delta shifts means only.
3. Estimation and replay build every weather covariate with the **identical
   code path**, so the estimated response functions transport to the perturbed
   climate (Appendix D.2 / E.5).

## 3. Step-by-step specification

**S1 — Temporal collapse: from anomaly series to a 12-month delta.**
The anomaly files retain the full hosing time axis; the coupling needs one
climatological delta per calendar month. For each model × protocol × variable,
average the anomaly over the **settled plateau** of the run (the LOESS-settled
window already defined for the delta diagnostics in Appendix C; operationally
the last third of the run is the fallback):

    Δ_m = mean over plateau years of  [X_hos(m, y) − clim_pic(m)] ,   m = 1…12

This yields twelve maps per variable. Averaging ~17–33 realisations of each
calendar month suppresses internal variability, and — as a side benefit —
collapsing to climatological months makes the models' different calendars
(proleptic-Gregorian vs 360-day) irrelevant. The result is an
**equilibrium-level** delta in the sense of Appendix C.7; a transient
(year-by-year) application remains possible later from the same files but is
not proposed here.

**S2 — Precipitation ratio: ratio of means, never mean of ratios.**
The multiplicative factor for month *m* is constructed from the collapsed
quantities:

    R_m = ( clim_pic(m) + Δpr_m ) / clim_pic(m)  =  P̄_hos(m) / clim_pic(m)

i.e. the ratio of the plateau-mean hosing precipitation to the control
climatology. It is **not** the time-average of the stored monthly ratio fields:
averaging ratios weights the noisiest (driest) months astronomically — the
stored fields contain values up to 3.7 × 10⁶ and, for EC-Earth3, negative
values of order 10¹⁸ where the ~500-year control climatology is a float-noise
negative (§4). The ratio-of-means form is algebraically the "1 + relative
change" of the additive delta and is far better conditioned. The stored ratio
files keep their role as per-month diagnostics.

**S3 — Regrid to the baseline grid.**
Regrid to the 0.5° ISIMIP grid following Appendix E.3: **bilinear** for the
three temperature deltas. For precipitation, regrid the *ingredients*, not the
ratio: remap `clim_pic(m)` and `Δpr_m` **conservatively** (preserving the
area-integrated water flux, which is meaningful for fluxes), then form `R_m` on
the target grid. Interpolating an already-formed ratio would (a) apply
flux-conserving weights to a dimensionless intensive quantity, where
"conservation" has no physical meaning, and (b) smear the extreme dry-cell
values of §4 into their neighbours before they can be masked. This keeps
Appendix E's promise that regridding is the single genuine interpolation of the
pipeline, confined to the coupling step.

**S4 — Sanitise the ratio (on the target grid).**
Two safeguards, in this order, with counts reported:

*Dry-cell mask.* Where the regridded control climatology is below a physical
floor — proposed: clim_pic(m) ≤ 0.01 mm/day (≈ 1.16 × 10⁻⁷ kg m⁻² s⁻¹) — set
`R_m := 1` (no perturbation). Rationale: where the control has essentially no
rain, the ratio estimates 0/0; no meaningful relative signal exists there, and
the threshold also absorbs the negative float-noise climatology cells found in
EC-Earth3. Over the European window the affected fraction is expected to be
very small (the measured global extremes concentrate in deserts and polar
cells); the actual fraction is a reported diagnostic, not an assumption.

*Clamp.* Bound the surviving ratios to a physical band, proposed
**R_m ∈ [0.1, 10]**: a tenfold monthly-mean change at cell scale already
exceeds any plausible forced response in this domain, so values beyond the band
are treated as denominator noise, not signal. The prior measurement on the
stored monthly fields (Europe: max ≈ 110, ~0.3 % of cell-months above 5)
suggests the plateau-mean `R_m` will rarely touch the bounds over land Europe —
the clamp is a guard rail, not a shaper. Sensitivity: rerun with [0.05, 20]
and report the change in final impact numbers (§6).

**S5 — Apply, day by calendar month (Appendix E.4 unchanged).**
For every baseline day *d* in month *m(d)*, per cell:

    tas'    = tas_ISIMIP    + Δtas_m      tasmin' = tasmin_ISIMIP + Δtasmin_m
    tasmax' = tasmax_ISIMIP + Δtasmax_m   pr'     = pr_ISIMIP · R_m

Monthly steps (rather than day-interpolated deltas) are proposed: the
downstream indicators are monthly or window aggregates of daily hinges, so a
smooth sub-monthly delta would change results negligibly while complicating the
audit trail. This can be revisited if a daily-threshold indicator proves
sensitive to month-boundary steps.

**S6 — Physical repair of the diurnal triple.**
The three temperature deltas are applied independently, so in months where
Δtasmin − Δtasmax exceeds a day's baseline diurnal margins the perturbed triple
can cross. Enforce, per day and cell:

    tasmin' := min(tasmin', tas')      tasmax' := max(tasmax', tas')

and **report the fraction of day-cells altered** (expected: rare and tiny in
magnitude; if it is not, that itself is a finding about the delta fields).
`pr' ≥ 0` holds by construction of the multiplicative form.

**S7 — Rebuild the covariates with the estimation code.**
From (`tas'`, `tasmin'`, `tasmax'`, `pr'`) recompute every impact covariate by
calling the same functions used in estimation (`weather_indicators.R`: capped
GDD, heat above 28 °C, frost, precipitation windows; within-day HDD/CDD hinges
with population weights). No re-implementation: identical code path, per the
transportability requirement.

**S8 — Scenario pair and ensemble.**
Scenario A = unperturbed ISIMIP3b baseline; Scenario B = the same days with the
delta applied (E.5). B − A isolates the AMOC signal with weather, bias and
scenario held fixed. The full exercise runs the six delta sets (3 model
configurations × 2 protocols) against the one common baseline, so the spread of
B − A across the six is the (small-sample) structural uncertainty of the AMOC
signal; results are reported per member and as the range.

**S9 — Diagnostics bundle (produced by the same run).**
(i) *Zero-delta identity:* running the module with Δ = 0, R = 1 must reproduce
Scenario A exactly — the machinery adds nothing. (ii) *Applied-mean check:* the
month-mean difference B − A equals Δ_m (and the pr ratio equals R_m) by
construction; verified numerically as a wiring test. (iii) Masked and clamped
cell fractions per month (S4), repair fraction (S6). (iv) Europe maps of the
applied Δ and R for visual comparison with the Appendix C delta maps.

## 4. Why the stored ratio files are not applied directly

Measured on the archived fields (2026-07-14): global maxima of the monthly
ratio are ≈ 3.6 × 10³ (HadGEM3-LL), ≈ 1.9 × 10³ (HadGEM3-MM), ≈ 3.7 × 10⁶
(EC-Earth3 u03); EC-Earth3 additionally contains **negative** ratios of order
−10¹⁸, because its long-term control climatology is a tiny negative number in a
few dry cells (float-level noise in regridded spectral-model precipitation).
These values are not physical signal; they are the 0/0 limit of the ratio
estimator. The stored fields were deliberately written unclipped (faithful
record); the consequence, adopted here, is that **all sanitisation is the
application step's responsibility** — mask first, clamp second, on the target
grid, with the affected fractions reported rather than hidden.

## 5. Deferred alternatives

| Choice made | Alternative kept open | Trigger to revisit |
|---|---|---|
| plateau-mean 12-month delta | transient year-by-year application | interest in the *path*, not the settled state |
| clamp [0.1, 10] | [0.05, 20] sensitivity run | clamped fraction over land Europe not negligible |
| monthly step deltas | daily interpolation of Δ | threshold indicator sensitive at month edges |
| ISIMIP3b carrier (E.2) | same module on E-OBS (prototype exists for energy) | head-to-head with estimation-grid replay |

## 6. Limitations (inherited and specific)

The delta shifts the monthly-mean seasonal cycle only: daily and within-day
variability, wet-day frequency, and the shape of extremes are inherited from
the baseline (E.6). The bias adjustment is assumed stationary under the
perturbation. One ensemble member per model; one baseline ESM. The clamp band
is a judgement call — defensible, stated, and tested by sensitivity, but not
derivable from first principles. HadGEM3's baseline is the first archived
control century rather than a verified parallel window (Appendix A.3; measured
sensitivity +0.11 K global, ≤ 6.1 K in isolated sea-ice cell-months).

## 7. Implementation note

One script (mirroring the structure of the existing energy-side prototype
`scenario_replay.R`, which already demonstrates S5/S7/S9-i for temperature on
the observational grid): CDO for S1–S3 (timmean over plateau, ymonmean, remap),
R for S4–S9. One practical prerequisite: the HadGEM anomaly script currently
deletes the piControl climatologies after use; keeping them (one-line change)
provides the S2 denominator without recomputation.

---

**References.**
Hempel, S., Frieler, K., Warszawski, L., Schewe, J., & Piontek, F. (2013). A
trend-preserving bias correction — the ISI-MIP approach. *Earth System
Dynamics*, 4, 219–236.
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975–1995.
Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling
with ISIMIP3BASD (v1.0). *Geoscientific Model Development*, 12, 3055–3070.
