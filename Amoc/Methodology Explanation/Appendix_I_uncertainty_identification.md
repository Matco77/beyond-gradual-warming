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

A natural check exists and has not been run: the **between-country** relationship
spans precisely the range into which the scenario moves. If the cross-sectional
slope resembles the within slope, the extrapolation has an empirical basis.

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
either estimated response function. The obvious next specification would carry window
length, or an indicator of completion, among the regressors; that has not been
estimated here.

---

## I.8 Summary of what is and is not supported

**Supported.** The energy results have a stable sign across hosing years and across
models: AMOC weakening raises household energy demand, with heating dominating
cooling by roughly eight to one. The magnitude is subject to the extrapolation of
§I.4 and should be reported with it.

**Not supported.** Any statement about the sign or magnitude of the crop yield
effect. The fixed-window figure fails on two independent grounds — an unidentified
thermal coefficient (§I.5) and a band that crosses zero in 16 of 20 bins (§I.3) —
and the thermal-time alternative is not identified either (§I.6).

**Not examined.** The cross-sectional check of §I.4; the selection induced by the
812-of-1 507 region coverage (Appendix H, §H.5); the stability of the crop-area
weights; and whether the ΔSv definition inherited in Appendix F is the appropriate
measure of AMOC state.
