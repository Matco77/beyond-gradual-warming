# Appendix C — Coupling the AMOC Signal onto the ISIMIP3b Baseline

This appendix documents how the stylised AMOC-weakening signal of Appendix A —
the monthly hosing-minus-piControl anomaly fields on each CMIP6 model's native
grid — is combined with a bias-adjusted daily baseline to construct the
perturbed ("AMOC-stressed") climate that drives the European energy and crop
impact models. The construction follows the *delta (change) method*: the climate
model supplies only the **change** it projects, while the **absolute, daily
climate** is taken from an observation-adjusted reference. This isolates the AMOC
response from the models' absolute-state biases, which the impact indicators —
being threshold-based — would otherwise transmit.

## C.1 Two archives, two roles

The perturbed climate is built from two independent products, each used only for
the quantity it estimates reliably:

| Archive   | Product                                          | Role in the coupling |
|-----------|--------------------------------------------------|----------------------|
| CMIP6     | NAHosMIP hosing − piControl anomaly (Appendix A) | the **delta**: the AMOC-attributable *change*, used as a difference only |
| ISIMIP3b  | bias-adjusted daily fields                       | the **baseline**: the absolute, daily climate onto which the delta is added |

CMIP6 output enters exclusively as the hosing-minus-control difference, in which
the shared model bias cancels (Appendix A.2); its absolute state, and its monthly
resolution, are never used to set the climate the impact models see. Conversely
the ISIMIP3b baseline supplies the absolute level and the daily and within-day
variability that the coarse, monthly delta cannot.

## C.2 The ISIMIP3b baseline

The baseline is the ISIMIP3b bias-adjusted product for a single Earth-system
model, **GFDL-ESM4** (`r1i1p1f1`), statistically bias-adjusted and downscaled to
the observational reference **W5E5** with the ISIMIP3BASD method (Lange, 2019),
under scenario **SSP3-7.0**. The four variables `tas`, `tasmin`, `tasmax`, `pr`
are used at **daily** resolution on the **0.5° × 0.5°** ISIMIP grid, over the
European window (lon −13.2° to 45.5°, lat 34.2° to 71.7°) and the years
2015–2030.

Two properties make this the appropriate baseline. First, the impact indicators
are **threshold-based** — cooling degree-days accumulate only above 24 °C,
growing degree-days are capped at 28 °C, frost is counted below 0 °C — so an
uncorrected absolute-temperature bias of a few kelvin would move days across
these thresholds and corrupt the covariate; the W5E5 bias adjustment removes that
bias. Second, the indicators integrate the **diurnal cycle**, which requires a
physically consistent daily (`tas`, `tasmin`, `tasmax`) triple; ISIMIP3b provides
`tasmin` and `tasmax` natively, so the reconstruction of Appendix B is required
only on the CMIP6 delta side (EC-Earth3), not on the baseline.

## C.3 Spatial coupling

The delta fields live on each model's native grid (Appendix A.2); the baseline is
on the 0.5° ISIMIP grid. The delta is therefore regridded onto the ISIMIP grid —
**bilinearly for temperature** and by **first-order conservative remapping for
precipitation** (which preserves the area-integrated water flux, a property
bilinear interpolation does not guarantee). This is the **single genuine
interpolation** in the entire pipeline: all anomaly construction in Appendix A is
deliberately interpolation-free on native grids, precisely so that regridding
never enters the hosing-minus-control subtraction and is confined to this
coupling step.

## C.4 Temporal coupling

The delta is a **monthly climatology** — twelve values per cell, one per calendar
month (Appendix A.3) — while the baseline is daily. The two are matched by
calendar month: every baseline day in month *m* receives that month's delta.
Writing *d* for a day falling in calendar month *m*(*d*), and dropping the cell
index,

    T_AMOC(d) = T_ISIMIP(d) + ΔT(m(d))          (temperature, additive)
    P_AMOC(d) = P_ISIMIP(d) · R(m(d))           (precipitation, multiplicative)

applied to `tas`, `tasmin`, and `tasmax` with their respective deltas.
Temperature is perturbed **additively**, because a shift of an interval-scale
quantity (K) is the meaningful operation, and it preserves the baseline diurnal
range and day-to-day sequence. Precipitation is perturbed **multiplicatively**,
through the ratio field R = hosing / piControl-climatology exported for `pr`
alongside the additive anomaly (Appendix A): a non-negative, ratio-scale variable
is scaled rather than offset, which keeps `P_AMOC ≥ 0` and rescales wet and dry
days proportionally.

## C.5 The scenario pair

The coupling defines a matched pair of daily climates on one common baseline,

    Scenario A (baseline) = ISIMIP3b, unperturbed
    Scenario B (AMOC)     = ISIMIP3b + delta      (the equations of C.4)

Because A and B share the identical bias-adjusted daily baseline, their
difference B − A isolates the imposed AMOC signal, with baseline bias and
internal day-to-day weather held fixed. Each scenario yields AMOC-consistent
daily `tas`, `tasmin`, `tasmax`, `pr` on the 0.5° grid over 2015–2030; these
drive the European energy-demand and crop-yield response functions, with every
weather covariate rebuilt by the same construction used in estimation, so that
the estimated response is transportable to the counterfactual climate.

## C.6 Assumptions and limitations

The delta is a shift of the **monthly-mean** seasonal cycle: like the
extreme-temperature reconstruction of Appendix B, it is a first-order
construction that transmits the change in the *mean* but not any change in the
*shape* of the sub-monthly distribution — the daily and within-day variability of
Scenario B is inherited unchanged from the ISIMIP3b baseline. The
additive-temperature and multiplicative-precipitation choices preserve,
respectively, the baseline diurnal range and the baseline wet-day structure.
Finally, the coupling assumes the ISIMIP3b bias adjustment, calibrated on the
observational reference, remains valid under the perturbed state (a stationarity
assumption on the bias, parallel to that of Appendix B.4). A single source model
(GFDL-ESM4) defines the baseline, so baseline structural uncertainty is not
sampled; the model spread of Appendix A enters through the delta alone.

---

**References.**
Cucchi, M., et al. (2020). WFDE5: bias-adjusted ERA5 reanalysis data for impact
studies. *Earth System Science Data*, 12, 2097–2120.
Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1
(GFDL-ESM4.1). *Journal of Advances in Modeling Earth Systems*, 12, e2019MS002015.
Jackson, L. C., et al. (2023). Understanding AMOC stability: the North Atlantic
Hosing Model Intercomparison Project. *Geoscientific Model Development*, 16,
1975–1995.
Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling
with ISIMIP3BASD (v1.0). *Geoscientific Model Development*, 12, 3055–3070.
