# Appendix I — Uncertainty, Identification, and the Limits of the Estimates

This appendix reports what the impact estimates of Appendix H can and cannot
support. It covers three distinct questions, which are often conflated and are kept
separate here: how much the result moves with the particular hosing years that
happen to fall in a bin; whether the estimated coefficients identify a
dose–response relationship; and how far outside its estimation support each response
function is being evaluated.

The conclusion on the crop side is negative, and is stated as such.

---

## I.1 Why the uncertainty band is computed by re-running the chain

The natural shortcut — propagating the cell-wise standard deviation of the
perturbation field through the impact calculation — is unavailable, because the
indicators sit behind thresholds. A standard deviation cannot be pushed through a
hinge: the mean of the effects is not the effect of the mean field, and the
discrepancy is exactly the non-linearity the delta method exists to capture.

The band is therefore constructed by running the **entire chain** — perturbation,
indicators, aggregation, coefficients — separately on each individual hosing year of
a bin, and dispersing the final results. Scope: IPSL-CM6A-LR and EC-Earth3, 98 + 99
= 197 hosing years, both impact branches. The HadGEM configurations retain the
central estimate only.

The per-year fields come from the same builder, the replay from the same engine, and
the coefficients from the same application code as the central estimate, with only
the grouping key changed (Appendix F, §F.5). The band therefore cannot disagree with
the central estimate for reasons of implementation.

---

## I.2 The size of the non-linearity

Comparing the central estimate (the chain run on the bin-**mean** field) with the
mean of the per-year effects isolates the non-linearity directly. Across the 20 bins
of the two models:

| branch | central − mean of per-year effects | median |
|---|---|---|
| crop | +0.77 … +2.12 pp | +1.47 pp |
| energy, electricity | −0.41 … −0.11 pp | −0.18 pp |
| energy, natural gas | −0.03 … 0.00 pp | −0.01 pp |

The ordering follows the physics of the indicators. Heating degree days over
October–March are a near-linear integral of temperature, and the two quantities agree
to two decimals — for gas, the band **demonstrates that it was unnecessary**.
Electricity adds a hard 24 °C dead-band through the cooling term and shows a small
gap. The crop indicators carry two hard temperature thresholds, a step function at
0 °C for frost, and a quadratic precipitation term, and the gap reaches 2 pp — the
same order of magnitude as the effect itself in the weaker bins.

This is a useful result in both directions: it quantifies where the delta method's
non-linearity matters, and where it does not.

---

## I.3 Dispersion across hosing years

The more consequential finding is the width of the band relative to the signal.

| | mean | s.d. | range across the bin's years |
|---|---|---|---|
| crop, IPSL, −7.5 Sv | +4.59 % | 3.82 | **−4.03 % … +12.34 %** |
| crop, EC-Earth3, −3.5 Sv | +1.90 % | 7.42 | **−10.01 % … +10.94 %** |
| electricity, IPSL, −8.5 Sv | +6.63 % | 3.04 | +1.40 % … +10.28 % |
| gas, IPSL, −8.5 Sv | +10.17 % | 4.50 | +1.81 % … +16.82 % |

**Sixteen of the twenty crop bins have a range that crosses zero**, against six of
twenty for each energy branch. On the crop side the *sign* of the effect is not
stable with respect to which hosing years happen to populate a bin. The energy bins
remain positive throughout.

This is an independent line of evidence from the specification problem of §I.5, and
it points the same way: the fixed-window crop result does not support a claim about
the sign of the impact.

**A related check.** The IPSL central estimates are not monotone in ΔSv — bin −8.5
returns a larger effect than bin −9.5. Testing the per-year effects of the two bins
against each other gives a difference of 1.70 pp with Welch *t* = 1.66 and
*p* = 0.126. The non-monotonicity is sampling noise, not structure, and the small
occupancy of the extreme bins (5 and 8 years) is the reason.

---

## I.4 Extrapolation beyond the estimation support

A benchmark with unit fixed effects, year fixed effects and unit-specific quadratic
trends discards the cross-section deliberately — it is confounded — and identifies
the coefficient from year-to-year variation alone. The relevant question is then how
far the scenario moves the regressor **relative to that variation**.

| regressor | within s.d. | share of raw variance surviving the fixed effects | scenario shift at IPSL −8.5 Sv |
|---|---|---|---|
| HDD calendar | 115.9 | 1.3 % | **+9.4 s.d.** |
| HDD Oct–Mar | 105.6 | 3.2 % | +6.6 s.d. |
| CDD JJA | 30.7 | 4.9 % | −1.8 s.d. |
| gdd, fixed window | 62.3 | 3.9 % | −5.7 s.d. |

The electricity coefficient is identified on an annual wobble of about 116 HDD —
some 4 % of a mean of roughly 3 018 — and the scenario asks it to predict a shift of
about 1 100 HDD. The linear functional form is being evaluated some nine standard
deviations outside the range in which it was tested. Only the cooling term stays at
a defensible distance.

This does not invalidate the energy result, and the heating–consumption relationship
is physically closer to linear than most. But the honest description is *a declared
extrapolation*, not *a robust estimate*, and the earlier characterisation of the
energy branch as raising no doubts was wrong.

## I.4b The cross-sectional check

The check proposed above was run: does the **between-unit** relationship — the
classic between estimator, each unit collapsed to its time-mean, regressed with no
fixed effects — resemble the within slope the benchmark identifies? The between
estimator draws on a different source of variation (permanent differences across
units) than the within estimator (year-to-year wobble inside a unit), so agreement
between the two is not guaranteed by construction; if it holds, the within slope's
extrapolation to an out-of-sample level has support beyond the narrow range that
identifies it. Script: `Amoc/Code/cross_section_check.R`.

**Crop (gdd_mj, between-NUTS3, on the benchmark's own estimation sample):**

| crop | n (NUTS3) | within (t) | between (t) | ratio |
|---|---|---|---|---|
| Soft wheat | 787 | −3.06×10⁻⁴ (−2.28) | −1.20×10⁻³ (**−10.59**) | 3.93 |
| Durum wheat | 345 | −3.72×10⁻⁴ (−2.36) | −5.10×10⁻⁴ (−3.40) | 1.37 |
| Spring barley | 603 | −5.93×10⁻⁴ (−4.34) | −6.94×10⁻⁴ (−6.87) | 1.17 |
| Winter barley | 730 | −6.37×10⁻⁴ (−4.14) | −1.03×10⁻³ (−10.10) | 1.61 |

Same sign on all four crops, and the between estimate is *more* precisely
determined than the within one (787 NUTS3 identify a slope far more tightly than
the within estimator's year-to-year wobble does), 1.2 to 3.9 times larger in
magnitude. This is real corroboration of one specific thing: **the sign of the
gdd_mj coefficient is not an artefact of the within-estimator's narrow identifying
variation** — an independent source of variation, differences between regions
rather than between years, points the same way.

**It does not corroborate the causal interpretation.** The phenological-misalignment
problem of §I.5–§I.6 is a property of the fixed Mar–Jul window itself, and that
window is fixed in the between comparison exactly as it is in the within one: a
warm region, like a warm year, is sampled at a later phenological stage by a
calendar window that does not move. The between check cannot distinguish "colder
regions have a genuinely different thermal dose–response" from "colder regions
also have systematically different phenological timing relative to a fixed
calendar window" — both would produce the same negative between-NUTS3 slope. So
this result licenses extrapolating the *sign* of the fixed-window effect with more
confidence than before; it does not rehabilitate spec A as a dose–response
estimate, and does not bear on spec C at all.

**Energy (between-country, on each benchmark's own estimation sample):**

| regressor | n (countries) | within (t) | between (t) | ratio |
|---|---|---|---|---|
| HDD calendar | 29 | 6.52×10⁻⁵ (4.20) | 1.93×10⁻⁴ (1.06) | 2.95 |
| CDD JJA | 29 | 1.53×10⁻⁴ (2.73) | 3.25×10⁻⁴ (0.24) | 2.13 |
| HDD Oct–Mar | 25 | 1.45×10⁻⁴ (2.49) | −4.49×10⁻⁴ (−0.95) | −3.11 |

**Inconclusive**, and should be reported as such rather than as either confirming or
refuting the extrapolation. With only 25–29 countries the between estimator has
essentially no power: none of the three point estimates is distinguishable from
zero (|t| ≤ 1.06), and HDD Oct–Mar's point estimate **changes sign** relative to
the within slope, though that reversal is itself not statistically meaningful.
HDD calendar and CDD JJA at least have the right sign and a plausible order of
magnitude (2–3×, comparable to the crop ratios); HDD Oct–Mar does not even offer
that. The extrapolation caveat of §I.4 therefore stands where it was — this check
had the power to corroborate it for crop, and did; it did not have the power to do
either for energy.

---

## I.5 The crop specification problem

The fixed-window crop result is driven entirely by a negative coefficient on growing
degree days. Since a negative thermal coefficient on C3 cereals runs against prior —
the 28 °C cap exists precisely to separate useful warmth from heat stress, which is
carried separately by the heat term — three alternative explanations were tested.

**Collinearity between the thermal and heat terms — rejected.** The within
correlation, after partialling out the benchmark's fixed effects, is +0.52 to +0.57,
with a variance inflation factor of 1.67 to 1.78. Dropping the heat term moves the
thermal coefficient by −13 % to +33 % but never changes its sign or order of
magnitude.

**Heterogeneity across climates — rejected.** Refitting within terciles of regional
mean thermal accumulation gives a negative coefficient in **all twelve** cases (four
crops × cold/medium/warm). The pooled coefficient is not an average concealing
opposite slopes; it is negative even in the coldest regions.

**The fixed calendar window — confirmed.** Separating the two features of the
thermal-time construction:

| crop | A: capped, Mar–Jul | B: uncapped, Mar–Jul | C: uncapped, thermal window |
|---|---|---|---|
| Soft wheat | −3.06 × 10⁻⁴ (t −2.28) | −3.09 × 10⁻⁴ (t −2.31) | −3.26 × 10⁻⁴ (t **−1.05**) |
| Durum wheat | −3.72 × 10⁻⁴ (t −2.36) | −3.78 × 10⁻⁴ (t −2.40) | −9.80 × 10⁻⁴ (t −2.56) |
| Spring barley | −5.93 × 10⁻⁴ (t −4.34) | −5.89 × 10⁻⁴ (t −4.40) | **+4.17 × 10⁻⁴** (t +0.91) |
| Winter barley | −6.37 × 10⁻⁴ (t −4.14) | −6.41 × 10⁻⁴ (t −4.12) | −5.93 × 10⁻⁵ (t **−0.26**) |

The 28 °C cap is **inert**: columns A and B agree to three significant figures on all
four crops, which is expected, since days with a *mean* temperature above 28 °C are
rare in a European March–July. The window does all the work. Under the thermal-time
window soft wheat loses significance, spring barley changes sign, winter barley
collapses to zero, and only durum survives.

The interpretation is that a fixed calendar window samples a later phenological stage
in a warm year, so "more thermal accumulation in March–July" partly measures window
misalignment with the crop's actual cycle rather than a thermal dose–response.
Consistent with this, the heat term becomes properly identified once the window is
aligned (soft wheat *t* = −2.26, spring barley *t* = −6.96), whereas it was never
significant under the fixed window.

---

## I.6 Why the thermal-time specification is not the answer either

The thermal-time window cannot be adopted as the corrected specification, for a
reason that is mechanical rather than empirical.

| crop | within s.d. of the thermal term, fixed window | thermal window | ratio |
|---|---|---|---|
| Soft wheat | 62.3 (3.9 % of raw variance) | 4.8 (**0.0 %**) | 0.078 |
| Durum wheat | 64.0 (4.6 %) | 5.3 (**0.0 %**) | 0.083 |
| Spring barley | 58.3 (4.8 %) | 4.6 (**0.0 %**) | 0.079 |
| Winter barley | 57.5 (3.7 %) | 4.8 (**0.0 %**) | 0.083 |

The window is *defined* by accumulating a fixed quantity of thermal time, so the
thermal total inside it is very nearly constant by construction; the fixed effects
then absorb essentially all of what remains. The coefficient in column C of §I.5 is
estimated on effectively no identifying variation — it is not a different answer, it
is no answer. The scenario nevertheless displaces the regressor by up to 3.4 within
standard deviations, so the coefficient is both unidentified and extrapolated.

The correct statement is therefore **not** that the sign of the crop effect depends
on the specification. It is:

> The fixed-window specification is contaminated by phenological misalignment. The
> thermal-time specification cannot measure the thermal dose, by construction.
> **Neither identifies the dose–response relationship, and a third specification is
> required.**

Consistent with this, the component decomposition under the thermal-time window
attributes essentially nothing to the thermal term (−0.01 % to −0.06 %), and the
residual effect is precipitation level against precipitation curvature almost
cancelling.

---

## I.7 Where the crop signal actually goes

Under the thermal-time window the temperature signal does not disappear; it migrates
to the **extensive margin**, which neither specification carries as a regressor:

| | region-years failing to complete the cycle | mean window length (baseline 150 days) |
|---|---|---|
| IPSL, −8.5 Sv | **15.1 %** | 157 |
| IPSL, −9.5 Sv | 9.7 % | 152 |
| EC-Earth3, −8.5 Sv | 4.2 % | 140 |
| HadGEM3-GC3-1MM, −11.5 Sv | 31.2 % | 174 |
| HadGEM3-GC3-1MM, −13.5 Sv | **43.1 %** | 185 |

Under strong cooling the season lengthens and, in a growing share of regions,
accumulated thermal time never reaches the closing anchor within the calendar year.
This is a larger statement than any coefficient in §I.5, and it is not represented in
either estimated response function. The obvious next specification carries window
length among the regressors — it is estimated in §I.7b.

---

## I.7b Spec D: the thermal window plus window length

Spec D is spec C's five regressors (§I.6) with `n_day` — the number of days inside
the thermal window, already computed by the window construction (§G.6) but never
used as a regressor by either A or C — added:

$$\ln y = \beta_1\,\text{gdd} + \beta_2\,\text{heat} + \beta_3\,\text{frost} +
\beta_4\,\text{precip} + \beta_5\,\text{precip}^2 + \beta_6\,\text{n\_day} +
\text{NUTS}_3[t,t^2] + \text{year}$$

fitted on the identical sample as specs A and C. Script:
`Amoc/Code/amoc_impact_specD.R` (AMOC branch) and `isimip_impact_specD.R` (ISIMIP).

**`n_day` is identified where `gdd` under spec C was not.** Repeating the test of
§I.6:

| crop | within s.d. of n_day | share of raw variance surviving the FE | mean n_day |
|---|---|---|---|
| Soft wheat | 11.8 days | **35.6 %** | 148 |
| Durum wheat | 8.7 days | **37.3 %** | 152 |
| Spring barley | 12.1 days | **31.1 %** | 146 |
| Winter barley | 11.0 days | **32.1 %** | 149 |

Compare to gdd's 0.0 % under spec C (§I.6): window length is not absorbed by the
fixed effects the way accumulated thermal time inside the window is, because the
window's *length* is exactly the margin that responds to a shifted climate while its
*content* is pinned by construction. This is the identification spec C structurally
lacked.

**The coefficient on n_day itself is not uniformly significant.** Positive in three
of four crops, but only durum clears a conventional threshold:

| crop | n_day coefficient (t) | gdd (t) | heat (t) |
|---|---|---|---|
| Soft wheat | +8.04×10⁻⁴ (1.01) | −3.60×10⁻⁴ (−1.09) | −6.45×10⁻⁴ (**−2.07**) |
| Durum wheat | +2.89×10⁻³ (**2.01**) | −1.15×10⁻³ (**−2.72**) | +1.70×10⁻⁴ (0.83) |
| Spring barley | −1.51×10⁻⁴ (−0.17) | +4.24×10⁻⁴ (0.90) | −1.13×10⁻³ (**−6.81**) |
| Winter barley | +1.16×10⁻³ (1.30) | −1.17×10⁻⁴ (−0.45) | −8.33×10⁻⁴ (−0.54) |

A longer season is associated with higher yield for durum (clearly) and weakly for
soft wheat and winter barley; spring barley shows no relationship. Spec D is
therefore not itself a clean, uniformly significant finding — but it does not need
to be, to answer the question this section asks.

**The question it answers: does adding a real, identified regressor pull the
picture back toward spec A's positive result?** No.

| bin ΔSv | A (fixed window) | C (thermal window) | D (+ n_day) |
|---|---|---|---|
| IPSL, −8.5 Sv | **+8.77 %** | −0.72 % | −0.56 % |
| IPSL, −6.5 Sv | +4.69 % | −2.18 % | −2.65 % |
| EC-Earth3, −8.5 Sv | **+6.05 %** | −1.44 % | −1.90 % |
| HadGEM3-GC3-1MM, −13.5 Sv | **+10.17 %** | +2.00 % | +2.87 % |
| ISIMIP, IPSL-CM6A-LR | **−8.70 %** | −0.30 % | +0.30 % |
| ISIMIP, EC-Earth3 | **−10.30 %** | −0.47 % | −0.49 % |

D sits close to C in every case shown, on both branches, and neither is anywhere
near A's magnitude or, on most bins, its sign. Where D moves away from C it is
usually further from A, not back toward it (IPSL −8.5 Sv: C = −0.72 %, D = −0.56 %,
both far from A's +8.77 %; HadGEM3-GC3-1MM −13.5 Sv: C = +2.00 %, D = +2.87 %, both
far below A's +10.17 %). Adding the one regressor Appendix I identified as missing
does not rehabilitate spec A; it leaves the picture where spec C left it — small,
inconsistently signed across bins and models, an order of magnitude below spec A —
now with a specification whose window-length term is at least identified.

This still does not amount to a defensible point estimate of the crop effect. It
narrows what remains open: not "which of three specifications is right", but "the
effect, if any, is small and this family of specifications cannot pin down its sign
at this magnitude" — which is itself a finding, and a different one from where §I.6
left off.

---

## I.8 Summary of what is and is not supported

**Supported.** The energy results have a stable sign across hosing years and across
models: AMOC weakening raises household energy demand, with heating dominating
cooling by roughly eight to one. The magnitude is subject to the extrapolation of
§I.4 and should be reported with it.

**Not supported.** A point estimate of the crop yield effect, from any of the three
specifications estimated. The fixed-window figure (spec A) fails on two independent
grounds — an unidentified thermal coefficient (§I.5) and a band that crosses zero in
16 of 20 bins (§I.3). The thermal-time window (spec C) cannot identify the thermal
dose by construction (§I.6). Adding window length (spec D, §I.7b) is identified
where C was not, but does not rehabilitate spec A's magnitude or, on most bins, its
sign — it lands close to spec C on both the AMOC and the ISIMIP branch. The
cross-sectional check (§I.4b) corroborates the *sign* of the fixed-window gdd
coefficient from an independent source of variation, but cannot speak to the
phenological-misalignment problem, which is a property of the fixed window that the
between comparison shares with the within one.

**Supported, narrowly.** That the large positive crop effect of spec A is a
specification artefact and not a robust finding: three independent specifications
that address the phenological-misalignment problem (C, D) or draw on independent
variation (the between estimator, §I.4b, which corroborates only the *sign* of A's
own coefficient, not its causal reading) all place the effect far below spec A's
magnitude, mostly negative or near zero rather than the +2 % to +10 % spec A
reports.

**Partially examined, and mixed.** The energy extrapolation of §I.4 was checked
cross-sectionally (§I.4b) and found inconclusive rather than confirmed: with only
25–29 countries the between estimator has no power to distinguish its point
estimates from zero, and one of the three regressors (HDD Oct–Mar) flips sign in
the point estimate. The extrapolation caveat therefore stands as stated in §I.4,
neither strengthened nor weakened with any real precision.

**Not examined.** The selection induced by the 812-of-1 507 region coverage
(Appendix H, §H.5); the stability of the crop-area weights; and whether the ΔSv
definition inherited in Appendix F is the appropriate measure of AMOC state.
