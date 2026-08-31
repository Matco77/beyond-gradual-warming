# Appendix J — The ISIMIP Branch

This appendix documents the second forcing source run through the impact chain: the
ISIMIP3b bias-adjusted climate projections, used to benchmark the AMOC-weakening
effect against a standard warming trajectory. It covers what is actually available
from ISIMIP3b (which differs from what was originally assumed), how the monthly
delta is built from daily data, and the impact estimates, reported with the same
specification caveat as the AMOC branch.

Implementation: `Amoc/Code/isimip_delta_fields.R`, `scenario_replay_isimip.R`,
`isimip_impact.R`.

---

## J.1 What ISIMIP3b actually provides

The premise that ISIMIP3b atmospheric data is available natively at monthly
resolution does not hold. Querying the ISIMIP data API by `time_step` shows:

| `time_step` | variables |
|---|---|
| `daily` | tas, tasmax, tasmin, pr, hurs, huss, ps, rlds, rsds, sfcwind, prsn, fwi |
| `monthly` | ocean/marine forcing only (thetao, tos, so, uo, vo, siconc, chl, …) |
| `annual` | derived indices (temp-annual-max, pr-extreme-99-9, …) |

The four variables this branch needs are distributed **only** at daily resolution;
the monthly time step exists for a different community (ocean biogeochemistry) and
carries none of them. The monthly climatology used here is therefore built by
aggregating the daily archive rather than received pre-aggregated. This changes
nothing about the resulting number — a mean of daily values *is* the monthly mean —
but it means the delta-method symmetry with the AMOC branch (§J.2) is a design
choice enabled by aggregating available daily data, not a property the archive
handed over directly.

Two further points, both checked rather than assumed. First, EC-Earth3 **is**
present among the ISIMIP3b bias-adjusted forcings (14 climate models in total, not
the five sometimes cited as "the ISIMIP3b core set"), so both AMOC-branch models of
primary interest are available here too. Second, the required cutout — the
server-side operation that returns a regional subset instead of a 105.6 GB global
download — is a distinct service from the file server and the metadata API; it
returned HTTP 502 against the **retired** `v1` endpoint (easily mistaken for an
outage) and worked immediately against the live `v2` endpoint once identified from
the `isimip-client` package defaults.

---

## J.2 Symmetric construction with the AMOC branch

The ISIMIP branch is deliberately run through the **same delta method and the same
engine** as the AMOC branch (Appendix G), rather than evaluating indicators directly
on the ISIMIP daily levels. Both branches perturb the identical E-OBS baseline with
a monthly delta; only the delta's source differs. This is a methodological choice,
not a data constraint — ISIMIP's daily archive would support a levels-based
calculation — made so that a comparison between the two branches isolates the
climate signal. Running ISIMIP on its own daily levels while AMOC goes through a
delta would confound the comparison with a difference of *procedure*, not only of
forcing, and any climatological bias in either dataset would not cancel.

**Delta construction.** For each model and each calendar month, a climatology is
built by aggregating daily ISIMIP3b values over two 30-year windows — 1985–2014
(historical) and 2071–2100 (ssp126) — chosen to align with the archive's decade-file
boundaries while giving each window a full 30 years. The delta is additive for
`tas`, `tasmin`, `tasmax` (K) and a ratio for `pr` (dimensionless, ssp126 climatology
over historical climatology); the [0.1, 10] bound is applied by the engine after
interpolation, identically to the AMOC bin fields. Mean warming over the European
domain: **+1.60 K (IPSL-CM6A-LR)**, **+2.31 K (EC-Earth3)**, with a mild seasonal
cycle (1.3–1.9 K IPSL, 2.1–2.6 K EC-Earth3 across calendar months) — physically
reasonable for a low-emissions trajectory 30 years out. Non-finite precipitation
ratios: none inside the European box for either model.

**One scenario per model, not forty.** ISIMIP is a single emissions trajectory
(ssp126), not an ensemble indexed by AMOC state, so there is no bin structure here
— the per-model field plays the role a single bin plays on the AMOC branch. Grid:
0.5°, matching ISIMIP3b/W5E5's native resolution; latitude arrives **decreasing**
(75.75° → 25.25°) and is re-sorted increasing before interpolation, the same
requirement the AMOC bin builder handles for its own grids.

---

## J.3 The zero-delta check and the ordering diagnostic

Because the scenario runs through the identical engine, it inherits the identical
gate: with a zero delta the replay reproduces the historical files exactly. All
three branches passed with **zero** maximum difference — `energy` against `13.`
(1 155 rows, 33 countries), `crop` against `8.` (105 490 rows), `crop_gddwin`
against `11.` (52 742 rows, window length included).

The ordering diagnostic of Appendix G, §G.5b was also run on this branch:

| model | `tx < tg` | `tg < tn` | mean diurnal range |
|---|---|---|---|
| E-OBS baseline, unperturbed | 0.248 % | 0.499 % | 9.06 K |
| IPSL-CM6A-LR | 0.508 % | 0.940 % | 9.23 K |
| EC-Earth3 | 0.399 % | 0.684 % | 9.08 K |

The violation rate rises more here than on the strongest AMOC bins (Appendix G
reported a maximum of 0.258 %/0.514 %), because the ISIMIP scenario applies a single
delta drawn from a 30-year climatological difference rather than a bin-specific
mean, and does not benefit from the partial cancellation across bins. The rates
remain small in absolute terms and the mechanism by which they are harmless is
unchanged (§G.5b): `within_day` is invariant to the sign of the diurnal amplitude,
and `heat_daily`/`frost_daily` are evaluated on days where the indicator is zero in
any case.

---

## J.4 Impact estimates

Coefficients, weights and both crop specifications are read from the AMOC branch's
scripts (`amoc_impact.R`, `amoc_impact_gddwin.R`) without modification — the same
`crop_effect()`/`energy_effect()` functions are reused with the grouping key reduced
to `model`, since there are no bins here. Coverage is identical to the AMOC branch,
because it is a property of the estimation panel, not of the forcing: 812 of 1 507
NUTS3 regions, 15 countries, 54 %.

**Sign check.** ISIMIP warms where the AMOC branch cools, so a consistent chain
should return opposite-signed effects on every indicator that responded to AMOC
weakening. It does, in both branches:

| | AMOC (weaker, e.g. −8.5 Sv) | ISIMIP (ssp126, 2071–2100) |
|---|---|---|
| HDD calendar | rises | **falls** (−2.9 to −3.2 %) |
| CDD JJA | falls | **rises** (+1.3 to +1.4 %) |
| crop, spec A, gdd term | positive | **negative** (−8.7 to −10.2 %) |

This is a property of the arithmetic, not of the specification, and it holds
regardless of which crop specification is trusted — confirming the sign flip is
not an artefact of the ISIMIP construction.

**Crop, all three specifications, area-weighted Europe:**

| | spec A (fixed Mar–Jul) | spec C (thermal-time window) | spec D (C + window length) |
|---|---|---|---|
| EC-Earth3 | −10.30 % | −0.47 % | −0.49 % |
| IPSL-CM6A-LR | −8.70 % | −0.30 % | +0.30 % |

The same caveat as Appendix I applies without modification: spec A carries the
phenological-misalignment problem (fixed calendar window, negative `gdd`
coefficient), and spec C is mechanically unable to identify a thermal dose-response
(the window is defined by a fixed accumulation of thermal time, so the `gdd` term
inside it contributes essentially nothing — **+0.01 %** here, matching the near-zero
contribution found throughout Appendix I, §I.6). Spec D (Appendix I, §I.7b) adds
window length, which unlike `gdd` under spec C keeps real identifying variation
(31–37 % of its raw variance survives the fixed effects, against 0.0 % for `gdd`);
its effect here sits close to spec C's on both models, **not** back toward spec A's
magnitude. **None of the three totals is a defensible standalone estimate of the
crop effect**; all three are reported so the specification sensitivity is visible
on this branch exactly as it is on the AMOC branch, rather than presenting a single
confident number that the AMOC-branch analysis has already
shown not to be trustworthy.

**Energy, pop-weighted Europe:**

| | electricity | gas |
|---|---|---|
| EC-Earth3 | −1.79 % | −4.18 % |
| IPSL-CM6A-LR | −1.65 % | −4.36 % |

**Extrapolation.** Following Appendix I, §I.4, the ISIMIP scenario shift was checked
against the within-estimator standard deviation that identifies each coefficient:

| regressor | within s.d. | ISIMIP shift |
|---|---|---|
| HDD calendar | 115.9 | −4.2 to −4.4 s.d. |
| HDD Oct–Mar | 105.6 | −2.9 to −3.1 s.d. |
| CDD JJA | 30.7 | +2.7 to +2.8 s.d. |
| gdd, fixed window | 62.3 | +3.4 to +3.9 s.d. |

This is a real extrapolation — three to four standard deviations beyond the range
that identifies the coefficients — but smaller than the AMOC branch's most extreme
bins, which reached 9.4 s.d. on the same regressor. The ISIMIP branch is therefore
evaluated closer to, though still outside, its estimation support.

---

## J.4b What the Phase 4 markers actually compare

The impact-vs-Sv figures (Appendix G onward, `Amoc/Code/plot_impact_curve_*.R`) mark
two reference quantities on the AMOC curve for IPSL-CM6A-LR and EC-Earth3: a
vertical line at an AMOC weakening level, and a horizontal line at this branch's own
effect. Both come from **the same underlying CMIP6 simulation** — checked against
file naming and metadata on both sides, not assumed. The ocean-circulation
reconstruction behind the vertical line is keyed
`amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc` (Terhaar, msftyz) and carries the source
attribute `vo (Omon, gn), EC-Earth3 r1i1p1f1` for the other model; the ISIMIP3b
downloads of §J.1 are named `ipsl-cm6a-lr_r1i1p1f1_...` and
`ec-earth3_r1i1p1f1_...`. Same model, same `r1i1p1f1` ensemble member, same ssp126
experiment on both sides.

**What differs is the diagnostic, not the simulation.** ISIMIP3b bias-adjusts only
the near-surface atmosphere of that run (§J.1) and carries no ocean output at all,
so the AMOC weakening implied by the run has to be read off a separate ocean
circulation reconstruction of the identical CMIP6 experiment (Terhaar's for IPSL,
an own reconstruction from `vo` for EC-Earth3), not off anything ISIMIP3b provides.
The vertical line is therefore not sourced from ISIMIP3b — but it is not an
independent or unrelated quantity either: it is the ocean state of the exact
simulated future whose surface climate produces the horizontal line's effect. The
two lines are two windows onto one simulated world, which is what makes their
comparison on the figure meaningful in principle — but §J.4c qualifies exactly how
meaningful, because the two windows are not the same *moment* of that world.
Neither NAHosMIP nor its `hos` terminology has any part in this — ssp126
is a real greenhouse-forced scenario, not the artificial freshwater perturbation
the AMOC-branch bins (Appendix F) are built from; the two branches are compared
because they are run through the identical replay engine (Appendix G), not because
either forcing resembles the other.

---

## J.4c The two lines are not from the same point in time

§J.4b establishes that the vertical and horizontal markers come from the same
simulation. They do not come from the same **time window within it**, and this
matters more than the "same simulation" fact does, because it is what determines
whether the crossing point on the figure can be read as a same-time prediction
check.

The horizontal line (§J.2) is a climatological delta over 2071–2100 versus
1985–2014. The vertical line's `target_delta_sv` is defined (Appendix A/E,
predating this branch) as the single most negative **annual** ssp126 anomaly over
the *entire* available ssp126 run, not an average over 2071–2100. Checked directly
against the raw netCDF series (`Amoc/Code/anomaly/plot_amoc_sv_ssp126.R`'s source
files), not assumed:

| model | ssp126 series | year of the AMOC minimum | AMOC at that minimum | AMOC averaged 2071–2100 (the horizontal line's window) |
|---|---|---|---|---|
| IPSL-CM6A-LR | 2015–**2214** | **2173** | 8.42 Sv | 10.07 Sv |
| EC-Earth3 | 2015–2100 | **2054** | 11.10 Sv | 12.89 Sv |

For IPSL-CM6A-LR the ssp126 run extends to 2214 — a long-tail continuation past the
usual 2100 endpoint — and the single most negative year falls **a century after**
the 2071–2100 window the horizontal line is built from. For EC-Earth3 the minimum
falls at **mid-century**, and the AMOC has partially recovered by 2071–2100
relative to that minimum (12.89 Sv against 11.10 Sv at the low point). Neither
model's vertical marker describes the state of the ocean during the years the
horizontal line's surface-climate effect is computed over.

One number could not be independently reproduced and is flagged rather than
asserted: subtracting the 1850–1900 historical mean (12.61 Sv for IPSL, computed
from the same historical file) from the 8.42 Sv minimum gives −4.19 Sv, not the
−3.7049 Sv the k10 file reports. The gap is not resolved here — the script that
built that file predates this branch and is not in the repository, so the exact
definition (a specific smoothing, a different minimum criterion, a different
historical sub-window) cannot be checked. Reported as an open discrepancy, not
papered over with a matching number.

**What this means for the figure.** The crossing point is not "what the AMOC branch
predicts for the period the ISIMIP effect describes." It is two reference points on
the same simulated trajectory, read at the trajectory's most extreme ocean state and
at a fixed late-century climate window respectively — a looser, weaker comparison
than "same simulation" by itself would suggest, and the figure and its caption
should be read accordingly.

---

## J.4d The vertical marker was dropped

§J.4c already weakens the vertical marker to a same-trajectory-but-different-moment
comparison. A second, more basic problem removes what was left of its value, and the
figures (`Amoc/Code/plot_impact_curve_crop.R`, `plot_impact_curve_energy.R`) no
longer draw it.

**The AMOC bins and the ssp126 marker are not the same physical quantity even at a
shared moment.** The bin construction (Appendix C) is `hosing − piControl
climatology`: piControl carries no rising greenhouse forcing at all, so the bins
isolate a **pure hosing effect** against an otherwise static climate. ssp126's AMOC
weakening, by contrast, occurs **inside** a real, if modest, greenhouse-driven
warming — the ocean circulation change and the background warming are not
separable in that simulation the way they are by construction in the hosing
protocol. A vertical marker built from the ssp126 diagnostic would place a
mixed (circulation + warming) quantity on an axis that measures a pure
circulation effect. This is very likely the deeper reason the two branches disagree
in sign throughout Appendices H and J (e.g. crop spec A: positive under the AMOC
bins, −8.7 to −10.3 % under ISIMIP) — not only, and probably not mainly, the
time-window mismatch of §J.4c, which is a smaller effect on top of this one.

Given that, a crossing point between the AMOC curve and a marked ssp126 Sv level
would not have been a meaningful check under any choice of time window: fixing
§J.4c's mismatch (e.g. marking the 2071–2100 AMOC mean instead of the
whole-trajectory minimum) would not fix this second problem. The marker is
therefore dropped rather than patched.

---

## J.4e The horizontal marker was dropped too

§J.4d's reasoning is not specific to the vertical line. The AMOC bins and the
ISIMIP/ssp126 branch measure different physical quantities — pure hosing-driven
circulation change against an unforced baseline, versus circulation change mixed
with real background warming — regardless of which axis a reference to the other
branch is drawn on. A horizontal line at the ISIMIP effect, laid across an AMOC-bin
panel, still invites reading it against the curve underneath it: does the curve
reach that level, at what bin, is that "close" or "far" — questions the same
physical-quantity mismatch makes unanswerable. Keeping the vertical marker's problem
in mind while keeping the horizontal one would have been inconsistent, so it is
dropped too, and the Phase 4 figures (`Amoc/Code/plot_impact_curve_crop.R`,
`plot_impact_curve_energy.R`) now carry no ISIMIP reference of either kind.

This does not remove the ISIMIP branch's results from the thesis — only from this
one figure. The effect estimates of §J.4 remain the branch's standalone output,
reported in their own tables (`isimip_impact_*.csv`) and in the text of §J.4, on
exactly the same footing as the AMOC branch's tables in Appendix H: two results
from two different forcing experiments, each interpretable on its own, presented
side by side in prose rather than merged onto one axis that would imply a
comparability neither branch's construction supports.

---

## J.5 What this appendix does not establish

No uncertainty band is computed for this branch. Appendix I's band disperses the
chain across the individual hosing years that populate an AMOC bin; ISIMIP has no
analogous ensemble; each scenario is a single climatological delta between two
30-year windows. A comparable quantity would require inter-annual variability
*within* each 30-year window, which is a different question and is not addressed
here.

The between-country check proposed in Appendix I, §I.4 as a way to give the
extrapolation an empirical basis has not been run on this branch either.

Only ssp126 is used. ISIMIP3b makes other scenarios available (e.g. ssp370, the
scenario the pre-existing local GFDL-ESM4 files happened to use); no comparison
across emissions pathways is attempted.
