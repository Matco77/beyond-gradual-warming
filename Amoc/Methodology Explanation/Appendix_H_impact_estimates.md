# Appendix H — Applying the Response Functions to the Perturbed Climate

This appendix documents the step that converts perturbed climate indicators into
impacts: which coefficients are used and how they are obtained without risk of
drift from the estimation code, how the effect is computed and decomposed, how it
is aggregated, and — importantly — over what fraction of the scenario geography it
is defined at all.

Implementation: `Amoc/Code/amoc_impact.R` (fixed-window specification) and
`Amoc/Code/amoc_impact_gddwin.R` (thermal-time specification).

---

## H.1 The estimator

For each unit and year, the impact is the change in the linear index of the
estimated response function, evaluated only on the climate regressors:

$$\widehat{\Delta \ln y} \;=\; \sum_k \hat\beta_k \left( X_k^{\text{scen}} - X_k^{\text{hist}} \right)$$

The non-climate regressors — fuel price for energy, and implicitly population
through the per-capita outcome — are held at their historical values, so they
cancel identically in the difference. The counterfactual is therefore *the same
economy under a different climate*, not a joint projection of climate and economy.
Fixed effects and unit-specific trends also cancel, since the unit and the year are
the same on both sides of the difference.

Because the index is linear in the regressors, aggregation and evaluation commute
for fixed weights: the weighted mean of the effects equals the effect of the
weighted mean regressors. This is not true across *years* once thresholds are in
play, which is the subject of Appendix I.

**The coefficients are never re-estimated on perturbed data.** They are the
historical benchmark coefficients, applied to a difference of indicators.

---

## H.2 Obtaining the coefficients without drift

The benchmark specifications live in the estimation script. Retyping them into the
impact script would create a second source of truth that could silently diverge —
the same hazard that motivates a single shared definition of the indicator formulas.

Instead, the impact script evaluates the *named top-level definitions* of the
estimation script without executing it: the sample-construction function, the crop
set, the panel loader, and the crop benchmark formula are read directly from
`regression.R`. The energy formulas are built inside a function from a regressor
string and cannot be extracted by name; they are reconstructed, and the
reconstruction is verified against the committed coefficient tables.

Verification, against the tables saved by the estimation run:

| Specification | N | Coefficients |
|---|---|---|
| Soft wheat | 18 550 | gdd −3.056 × 10⁻⁴, heat −5.444 × 10⁻⁴, frost −3.590 × 10⁻³, precip 1.091 × 10⁻³, precip² −1.704 × 10⁻⁶ |
| Durum wheat | 7 282 | gdd −3.723 × 10⁻⁴, heat 2.080 × 10⁻⁴, frost −7.552 × 10⁻³ |
| Spring barley | 13 397 | gdd −5.933 × 10⁻⁴, heat −5.974 × 10⁻⁴, frost −2.189 × 10⁻³ |
| Winter barley | 15 726 | gdd −6.370 × 10⁻⁴, heat 4.451 × 10⁻⁴, frost −4.116 × 10⁻³ |
| Electricity | 844 | HDD calendar 6.523 × 10⁻⁵, CDD JJA 1.527 × 10⁻⁴ |
| Natural gas | 577 | HDD Oct–Mar 1.447 × 10⁻⁴ |

Both the coefficients and the sample sizes reproduce the committed tables exactly.

**Only the main crop window is used.** The replay produces both March–July and
April–August, but the benchmark is estimated on the March–July regressors only;
there is no April–August benchmark and therefore no coefficient to apply to it. The
second window is retained as an input for future robustness work, not silently
averaged into the result.

---

## H.3 The quadratic precipitation term

`precip²` is a regressor in its own right and is **recomputed from the perturbed
precipitation**, not obtained by shifting the historical squared term:

$$\Delta(P^2) = \left(P^{\text{scen}}\right)^2 - \left(P^{\text{hist}}\right)^2 \neq \left(\Delta P\right)^2$$

Treating it otherwise would misstate the curvature term wherever the perturbation is
not small relative to the level, which is most of the ΔSv range.

---

## H.4 Aggregation and weights

Crop effects are computed per NUTS3 × crop × year, aggregated to countries with
**fixed crop-area weights** (the mean sown area of that crop in that region over the
estimation years), then averaged over years, then aggregated across crops and
countries with total crop area for the European figure. Energy effects are computed
per country × year and aggregated to Europe with **fixed population weights**.

Weights are fixed at historical means so that the aggregate reflects the climate
effect rather than a drift in the cropping mix or in population.

---

## H.5 Coverage: where the estimate is defined

This is a material restriction and is stated explicitly because the merge that
imposes it is otherwise silent.

The scenario fields cover **1 507 NUTS3 regions**. The crop response function is
identified on the 816 regions that have yield data in the estimation panel, of which
**812 also carry a usable crop-area weight** (four have zero or missing sown area and
therefore receive no weight). Effects are computed on those 812 regions, in **15
countries** — 54 % of the scenario geography. The remaining 695 regions receive a
perturbed climate but no impact estimate, because a coefficient is transportable only
where it was identified.

For energy the outcome is national and the restriction is milder: **29 countries**
for electricity and **25** for natural gas, out of 33 in the weather aggregate.

Two consequences follow. First, the "European" crop figure is an aggregate over 15
countries, not over Europe. Second, the 812 regions are not a random sample of the
1 507: regions with harmonised sub-national yield statistics may differ
systematically from those without, in ways plausibly correlated with agricultural
intensity. **This selection has not been characterised**, and the aggregate should be
read as conditional on it.

---

## H.6 Decomposition by component

Contributions are reported per regressor, not only as a net. This is not
presentational: under cooling the components move in opposite directions and a small
net can conceal large offsetting terms.

IPSL-CM6A-LR, fixed-window specification, area-weighted European aggregate, effect on
yield in per cent:

| bin ΔSv | n | gdd | heat | frost | precip | precip² | **net** |
|---|---|---|---|---|---|---|---|
| −0.5 | 8 | +0.19 | −0.04 | +0.31 | −1.36 | +0.96 | **+0.06** |
| −3.5 | 8 | +6.42 | +0.37 | −2.48 | −2.20 | +1.63 | **+3.54** |
| −8.5 | 8 | +15.38 | +1.09 | −5.88 | −4.14 | +3.37 | **+8.77** |

At −8.5 Sv the net of +8.8 % contains a +15.4 % thermal term and a −5.9 % frost term,
while the precipitation level and curvature terms very nearly cancel (−4.1 % against
+3.4 %). Reporting only the net would conceal all of this.

The weakest bin returning ≈ 0 % on every component is an internal consistency check
that was not imposed: a near-zero AMOC perturbation produces near-zero indicator
changes and hence a near-zero effect.

Energy, effect on per-capita consumption:

| bin ΔSv | electricity (IPSL) | of which HDD / CDD | gas (IPSL) |
|---|---|---|---|
| −0.5 | +0.38 % | +0.40 / −0.02 | +0.29 % |
| −3.5 | +2.47 % | +2.92 / −0.44 | +4.73 % |
| −8.5 | +6.47 % | +7.41 / −0.87 | +10.16 % |

The heating response dominates the cooling response by roughly a factor of eight:
additional winter demand is not offset by reduced summer demand.

---

## H.7 Status of the two crop specifications

The fixed-window result above is reported **alongside** the thermal-time result, and
neither is presented as the estimate. The sign of the crop effect is not stable
across the two window definitions, for reasons documented in Appendix I: it is a
specification choice, not a finding. The corresponding thermal-time figures are
−0.7 % at −8.5 Sv against +8.8 % here.

The energy results are not affected by the window question — the energy windows are
calendar constructions (the heating year, the cooling season) that correspond to
actual consumption cycles rather than to a biological process whose timing shifts
with temperature. They are, however, subject to the extrapolation issue of
Appendix I, §I.4.

---

## H.8 What this appendix does not establish

The crop-area weights are means over the estimation period. The stability of the
cropping mix over that period has not been examined, and the weights have not been
validated against an independent source.

The coverage restriction of §H.5 is reported but not diagnosed: no test has been run
for whether the 816 regions differ systematically from the excluded 691 in ways that
would bias the aggregate.

Nothing in this appendix addresses whether the coefficients are causally
interpretable. That question — whether the estimated response is a dose–response
relationship or a correlate of something else — is the subject of Appendix I and, on
the crop side, is answered unfavourably.
