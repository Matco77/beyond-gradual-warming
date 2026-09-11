# Two identification strategies for the AMOC impact estimates

**Marco Bova — supervision memo, September 2026**

*Beyond Gradual Warming*: European crop yields and household energy demand under a stylised
AMOC weakening. Two response functions are estimated on historical observational data and then
transported to counterfactual climates. **Strategy A is estimated and the full replay is built.
Strategy B is a proposal, not yet estimated.** This memo sets both out and asks for a ruling on
five points.

---

## 1. What both strategies share

Neither of the following is in question.

**Two modules, each at the finest resolution its outcome exists at.** Crop: NUTS3 × crop × year,
yield from CropStatHarm (Ronchetti et al. 2024), 76,979 rows, 816 NUTS3, 15 countries, 1989–2023.
Energy: country × fuel × year, household gas and electricity final use from Eurostat `nrg_bal_c`,
2,924 rows, 29 countries, 1990–2024. The asymmetry is a data constraint — no sub-national
household-energy series exists.

**Weather always built on the E-OBS 0.25° daily grid, then aggregated up** — area-weighted to
NUTS3 for crop, population-weighted (fixed 2021 GISCO grid) to country for energy. Non-linear
transforms are applied per cell-day before any aggregation, since `E[f(T)] ≠ f(E[T])`.

**The identification claim for β.** Unit fixed effects, year fixed effects and unit-specific time
trends; weather enters as a deviation from the unit's own smooth path, net of common annual
shocks (Blanc & Schlenker 2017). Both strategies rest on this.

**The transport method.** Two monthly anomaly fields, each a difference within one model so the
additive component of model bias cancels:

- Δ^W = ISIMIP3b ssp126, 2071–2100 minus 1985–2014 (models: IPSL-CM6A-LR, EC-Earth3)
- Δ^A = NAHosMIP `u03-hos` minus matched piControl, grouped into 1 Sv bins of simulated
  overturning strength at 26°N (models: IPSL-CM6A-LR, EC-Earth3, HadGEM3-GC3-1LL, -1MM)

Temperature deltas are additive (K); precipitation is a multiplicative ratio. Both are applied to
the same observed E-OBS baseline, so climatological model bias cancels on each side and the two
arms differ by forcing rather than by method.

**The four arms.** I_W (warming alone), I_A (AMOC weakening alone, from the present climate),
I_{A|W} (AMOC weakening on a warmed background), I_tot. The conditional arm requires assuming the
two forcings superpose linearly and the AMOC response is not itself state-dependent; the
unconditional arms do not.

---

## 2. Strategy A — threshold exposures (estimated)

**Crop**, per crop *k* separately, on the March–July window (April–August as robustness):

```
ln(Yield_{r,t}) = β₁·GDD + β₂·Heat + β₃·Frost + β₄·Precip + β₅·Precip²
                + α_r + λ_t + δ_r·t + γ_r·t² + ε_{r,t}
```

with per cell-day `GDD = max(min(tg,28) − 5, 0)`, `Heat = max(tx − 28, 0)`, `Frost = 1(tn < 0)`,
`Precip = rr`. The 28 °C cap makes beneficial warmth and harmful extremes non-overlapping
(Schlenker–Roberts convention).

**Energy**, per fuel:

```
ln(E_{c,t}/Pop_{c,t}) = β₁·HDD + β₂·CDD + β₃·ln(Price) + μ_c + λ_t + δ_c·t + γ_c·t² + ε_{c,t}
```

with `HDD = 1(T<15)(18−T)` and `CDD = 1(T≥24)(T−21)` **integrated over the within-day temperature
path** — a sine centred on `tg` with amplitude `(tx−tn)/2` — rather than evaluated at the daily
mean. This matters for cooling: 34.5% of summer days straddle 24 °C, and evaluating at the mean
understates total CDD by roughly 1.7×. Electricity uses calendar-year HDD and JJA CDD; gas uses
October–March HDD and no CDD.

**Where the non-linearity is identified from.** The daily distribution. Three first-order
coefficients, each on a separately constructed exposure, rather than curvature in a single
variable.

**Transport.** The monthly delta is added to every E-OBS day of that calendar month and the
indicators are recomputed with the *same* construction module used in estimation. The hinges are
therefore evaluated on the shifted days — under cooling, days cross the 5/28 °C and 15/24 °C
thresholds differently, which a delta applied to a finished degree-day total cannot see. A
zero-delta self-check requires the engine to reproduce the historical indicators exactly before
any scenario is trusted.

**Inference.** Clustered on country (15 crop, 29 energy); two-way clustering; restricted wild
cluster bootstrap (Rademacher weights, null imposed; Cameron–Gelbach–Miller 2008, Roodman et al.
2019) on the fixed-effect-partialled model. Default standard errors never used.

**Estimated.** Electricity HDD ≈ 6.5×10⁻⁵, gas HDD ≈ 1.4×10⁻⁴, electricity CDD ≈ 2×10⁻⁴, all
correctly signed under the benchmark. Crop: see §4.

---

## 3. Strategy B — quadratics in seasonal means and second moments (proposed)

**Crop:**

```
ln Y_{r,t} = β₁·Tx + β₂·Tx² + β₃·Tn + β₄·Tn² + β₅·P + β₆·P²
           + α_r + λ_t + δ_r·t + ε_{r,t}
```

where Tx and Tn are growing-season means of daily maximum and minimum temperature, and Tx², Tn²
are **means of daily squares, not squares of means**, formed at cell-day level before aggregation.
Since `E[T²] = (E[T])² + Var(T)`, the regressor carries within-month daily variance, and part of
the identifying variation in β₂ comes from differences in daily dispersion at equal seasonal mean.

**Energy:**

```
ln(E/Pop) = θ₁·Tw + θ₂·Tw² + θ₃·Ts + θ₄·Ts² + θ₅·ln Price + μ_c + λ_t + δ_c·t + ν_{c,t}
```

with Tw the October(t−1)–March(t) mean and Ts the June–August mean; gas omits the summer terms.

Note: unit-specific trends are **linear** here, against quadratic in Strategy A.

**Where the non-linearity is identified from.** Curvature — the residual variation in Tx²
orthogonal to Tx, after unit FE, year FE and unit-specific trends.

**Transport, in closed form.** Shifting every day of month *m* in cell *g* by Δ gives

```
ΔTx  = Σ_m (n_m/D) Σ_g ω⁰ · Δ
ΔTx² = Σ_m (n_m/D) Σ_g ω⁰ · [2·Δ·B + Δ²]
```

from `E[(T+Δ)²] − E[T²] = 2ΔE[T] + Δ²`, which depends on the daily distribution only through its
mean. B is taken from the observational climatology, so no model temperature is ever squared. The
four arms then satisfy an exact identity, with the interaction term equal to `2β₂⟨Δ^W·Δ^A⟩` —
**proportional to β₂** — rather than obtained by subtraction.

**Inference.** Moving-block bootstrap over years (contiguous blocks of 3–5 years, all
cross-sectional units retained within a sampled year), preferred to country clustering on
few-cluster grounds; two-way clustering and spatial-HAC (Conley) reported alongside.

**Not estimated.** This would require rebuilding the E-OBS pipeline to produce monthly means and
monthly means-of-squares per cell, and re-estimating both response functions.

---

## 4. The problem neither strategy resolves: the growing-season window

This is the most consequential open item and it is independent of the A-vs-B choice.

Under Strategy A with the fixed March–July window, **β_gdd is negative** for the crops estimated
(soft wheat, durum wheat, spring barley, winter barley). The diagnosis in the code is that a fixed
calendar window confounds the thermal dose with phenological misalignment: a warm year samples a
later growth stage inside a window that does not move, so "more GDD in March–July" partly measures
misalignment rather than a dose–response.

Three specifications have been tried:

| Spec | Window | Result |
|---|---|---|
| **A** | Fixed Mar–Jul | β_gdd < 0. Sign does not survive moving the window. |
| **C** | Thermal time — window opens and closes at fixed accumulated GDD | β_gdd not identified **by construction**: the window is defined by accumulating a fixed thermal amount, so within-unit variation in GDD collapses to **0.0% of raw variance** and the fixed effects absorb what remains. Soft wheat loses significance, spring barley changes sign, winter barley collapses to zero. |
| **D** | Spec C plus window **length** (`n_day`, 57–322 days historically) as a regressor | The AMOC signal under a thermal window goes into how long the season takes, not into the dose. Estimated for the first time; not yet adopted. |

The per-year uncertainty band crosses zero in 16 of 20 ΔSv bins. The crop *components* (GDD,
Heat, Frost, Precip separately) remain informative under all three; the crop *net* currently does
not, and is flagged in the code as not reportable. The energy side has no equivalent problem.

Strategy B would inherit this: it fixes a common calendar window across crops and regions and
records that as a limitation.

---

## 5. Questions

**Q1 — The crop window.** Is spec D (thermal window + season length) the right resolution, or
does this call for a different design: crop- and region-specific phenological calendars, a
decomposition into two or three growth phases, or dropping GDD and identifying the crop response
off Heat, Frost and Precipitation alone? Related: if no specification yields a stable sign, is the
honest output the components rather than a net yield effect?

**Q2 — Where should the non-linearity be identified from?** Strategy A identifies it from the
daily distribution through hinge transforms; Strategy B from curvature in a seasonal mean after
three layers of temporal controls. Is there enough within-unit variation left to identify β₂ in B?
The stakes are asymmetric: under B the interaction term in the four-arm decomposition is
*proportional to β₂*, so a weakly identified curvature makes the state-dependence result — the
paper's central object — unidentified as well. A counter-consideration in B's favour: the hinge
thresholds in A (5, 28, 15, 24 °C) are literature conventions, not estimated, whereas a quadratic
imposes no such constants.

**Q3 — Trend order.** The estimated benchmark uses unit-specific **quadratic** trends: 1,632 trend
parameters across 816 NUTS3 over ~35 years. It is needed on the energy side — with year fixed
effects alone the spatially correlated common European weather is absorbed by λ_t, leaving HDD
insignificant and wrong-signed for gas. But two trend parameters per unit may also absorb
medium-frequency climate variation. Linear trends and restricted cubic splines have not been
tested. Which is the defensible benchmark, and should the choice be made on out-of-sample
predictive performance by held-out year?

**Q4 — Transport under a non-marginal shift.** β is identified from inter-annual anomalies; the
AMOC scenario moves the regressors several within-unit standard deviations beyond the estimation
support (measured and reported per bin). The envelope theorem licenses panel variation for a
*marginal* change. Two sub-questions: (i) how far does that argument stretch here; (ii) a
long-difference design (Burke & Emerick 2016) across the 816 NUTS3 is estimated — but should it be
read as a long-run *bound*, given that the declining-sensitivity literature concerns adaptation to
*heat*, and air conditioning has been found not to mitigate cold? The direction of the
short-run/long-run wedge in a **cooling** scenario appears unsigned.

**Q5 — Inference.** The wild cluster bootstrap addresses the few-cluster problem in the
cross-section. It does not address serial dependence in the outcome, nor spatial correlation in
the interpolated weather beyond what country clustering absorbs. Is that sufficient, or should
spatial-HAC (Conley) and a block bootstrap over years be added as in Strategy B?

**Q6 (secondary) — Which arm is the headline?** I_A requires no additivity assumption but answers
a counterfactual question about the present climate. I_{A|W} answers the policy-relevant question
but assumes the AMOC response measured on an unforced background transfers to a warmed one.
Reporting both is possible; which should carry the result?

---

**References.** Blanc & Schlenker (2017) *REEP* 11(2); Burke & Emerick (2016) *AEJ: Policy* 8(3);
Cameron, Gelbach & Miller (2008) *ReStat* 90(3); Cornes et al. (2018) *JGR-Atmospheres* 123;
Fisher, Hanemann, Roberts & Schlenker (2012) *AER* 102(7); Jackson et al. (2023) *GMD* 16;
Mérel & Gammans (2021) *AJAE* 103(4); Ronchetti et al. (2024) *ESSD* 16; Roodman et al. (2019)
*Stata Journal* 19(1); Schlenker & Roberts (2009) *PNAS* 106(37).
