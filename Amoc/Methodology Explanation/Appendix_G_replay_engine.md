# Appendix G — The Scenario Replay Engine

This appendix documents the machinery that turns a monthly perturbation field into
the regressors of the estimated response functions. It specifies the perturbation
itself, the order in which it is applied relative to the indicator thresholds, the
interpolation onto the observational grid, the aggregation to the estimation
geography, and the check that guarantees the replay reproduces the historical
record when the perturbation is zero.

Implementation: `Amoc/Code/scenario_engine.R`, with three runners on top of it
(`scenario_replay.R`, `scenario_replay_crop.R`, `scenario_replay_bins.R`).

---

## G.1 What "replay" means

The impact indicators are non-linear functions of **daily** values: growing degree
days are bounded below at 5 °C and above at 28 °C, heat accumulates only above
28 °C, frost is an indicator on the daily minimum crossing 0 °C, and the energy
indicators use dead-bands at 15 °C and 24 °C. None of these can be recovered from a
monthly mean.

The forcing, by contrast, is available only monthly: the NAHosMIP hosing
experiments archive monthly means. The two are reconciled by the delta method. The
observed daily record supplies the daily and within-day structure; the model
supplies only the *change*. Concretely, for each cell and each day:

$$T^{\text{scen}}_{c,d} = T^{\text{obs}}_{c,d} + \delta_{c,\,m(d)}, \qquad
P^{\text{scen}}_{c,d} = P^{\text{obs}}_{c,d} \cdot R_{c,\,m(d)}$$

where `m(d)` is the calendar month of day *d*. The delta of a month is applied to
every day of that month. Temperature is perturbed additively, precipitation
multiplicatively with the factor bounded to [0.1, 10]; a multiplicative factor
cannot drive precipitation negative in dry cells and scales with local climate,
which an additive delta would not.

The observational baseline is E-OBS v31.0e on a 0.25° grid, the same record on
which the response functions of Appendix D were estimated.

---

## G.2 The perturbation precedes the thresholds

The delta is applied to the daily field and the indicator hinges are then evaluated
on the shifted days. The alternative — perturbing an already-computed seasonal
degree-day total — would assume the thresholds bind identically before and after,
which is exactly what a change of climate violates.

The magnitude of the difference is not subtle. Taking the IPSL-CM6A-LR bin centred
on −3.5 Sv and applying it to March–July 2003 over the 1.86 million cell-days of
the crop grid:

| | observed | perturbed |
|---|---|---|
| mean daily maximum | 17.17 °C | 15.74 °C (−1.43 K) |
| cell-days above the 28 °C threshold | 246 883 | 197 960 (**−20 %**) |
| accumulated Heat | 790 355 | 644 551 (**−18.4 %**) |

A mean shift of −1.4 K removes a fifth of the qualifying days and cuts accumulated
Heat by 18 %, because Heat counts only the tail above 28 °C: a day at 29 °C stops
contributing altogether rather than contributing 3 units less. The relative response
of the indicator is roughly thirteen times that of the mean temperature. No
transformation applied to a finished seasonal total can reproduce this.

---

## G.3 Interpolation onto the observational grid

Model fields are interpolated onto E-OBS **cell centres** by bilinear interpolation
(`fields::interp.surface`). The same routine is used for temperature and for
precipitation.

The choice of routine is dictated by the data: `terra::resample` rejects the native
model grids with *"lat not regularly spaced"*, so a single path through
`interp.surface` serves all fields rather than maintaining two.

One deliberate inconsistency is retained. The legacy `delta_fields` runs sit on a
regular grid that `terra` accepts, and keep the `terra` bilinear resample with which
they were originally built, so that previously published numbers remain
reproducible. Both routines are bilinear on cell centres and agree to interpolation
precision. Similarly, the precipitation bound is applied on the model grid *before*
interpolation in the legacy path and *after* interpolation in the bin path; the two
differ only in cells adjacent to a bounded cell, and no bounded cell lies within one
model grid cell of any E-OBS cell used by either branch (Appendix F, §F.3).

Cells falling outside the model domain receive a neutral value — zero for an
additive delta, one for a multiplicative factor — never an extrapolated one.

---

## G.4 The three branches

The engine exposes three reduction paths, differing in geography, indicator set and
window definition. All three share the E-OBS reading, the perturbation and the
weighting machinery.

**`energy`** — population-weighted country aggregate. Daily HDD and CDD are computed
with the within-day integration over the diurnal cycle, then aggregated to countries
using fixed population weights, and summed into three seasonal aggregates: HDD over
the calendar year, HDD over October–March, and CDD over June–August. The Oct–Mar
season is a cross-year construction, so the season labelled with the first year of
the estimation window requires the preceding October–December; the engine therefore
reads one extra year at the start. Seasons with fewer than 150 days are dropped, as
in the historical construction.

**`crop`** — area-weighted NUTS3 aggregate on a fixed calendar window. Daily
indicators are aggregated to NUTS3 and summed by month, **rounded to three decimals
at the monthly stage**, and then summed over the window. The intermediate rounding is
not cosmetic: it reproduces the exact two-step arithmetic of the historical
construction, without which a zero perturbation would land ~10⁻¹³ away from the
historical value rather than on it. Both estimation windows are produced,
March–July (main) and April–August (robustness).

**`crop_gddwin`** — the same geography, but with window edges set by accumulated
thermal time rather than by the calendar. See §G.6.

Weights are fixed at their historical means throughout, so that a scenario aggregate
reflects the change in climate and not a change in the weighting.

---

## G.5 The zero-delta check

With a zero temperature delta and a unit precipitation factor, the engine is the
historical construction with a `+ 0` and a `× 1` inserted into it. Any departure from
the historical files is therefore a defect of the replay — wrong cells, wrong month
mapping, wrong window rule, wrong rounding — and not a climate signal. Since every
number downstream is a *difference* between a scenario and the historical baseline, a
defect here would bias all of them silently and invisibly.

The check is executed before any scenario is run and halts the pipeline on failure.
It is exact, not tolerance-based:

| branch | compared against | result |
|---|---|---|
| `energy` | `13.eobs_country_energy_weather_weighted.csv` | max abs. difference **0** on HDD calendar, HDD Oct–Mar, CDD JJA; 1 155 rows, 33 countries |
| `crop` | `8.eobs_nuts3_crop_weather_window.csv` | max abs. difference **0** on all four indicators, both windows; 105 490 rows, 1 507 NUTS3 |
| `crop_gddwin` | `11.crop_weather_gdd_window.csv` | max abs. difference **0** on all four indicators **and on window length**; 52 742 rows |

No missing-value pattern differs between replay and historical in any branch.

The same check served as the regression test when the three runners were
consolidated onto the shared engine: both previously published scenario files were
regenerated and are **byte-identical** to the committed versions (identical MD5), so
the consolidation provably changed nothing.

---

## G.5b Ordering of the perturbed temperature fields

Nothing in the delta method enforces `tn ≤ tg ≤ tx` after perturbation: three
independent additive deltas are applied to three fields, and

$$TX - TG = (tx^{\text{obs}} - tg^{\text{obs}}) + (\delta_{tx} - \delta_{tg})$$

can turn negative wherever the differential delta exceeds the observed diurnal
half-range on a given day. The differential deltas are not small: across the bin
fields `δtx − δtg` spans −1.15 to +4.61 K (IPSL, strongest bin) and −2.94 to
+4.22 K (HadGEM3-GC3-1MM). The engine therefore reports the violation rate as a
standing diagnostic (`order_check`), on the strongest bin of each model.

Measured on 4.43 million cell-days:

| | `tx < tg` | `tg < tn` | mean diurnal range |
|---|---|---|---|
| **E-OBS baseline, unperturbed** | 0.248 % | 0.499 % | 9.06 K |
| IPSL-CM6A-LR, −9.5 Sv | 0.151 % | 0.443 % | 9.21 K |
| EC-Earth3, −9.5 Sv | 0.245 % | 0.548 % | 8.86 K |
| HadGEM3-GC3-1LL, −8.5 Sv | 0.258 % | 0.514 % | 9.27 K |
| HadGEM3-GC3-1MM, −14.5 Sv | 0.139 % | 0.278 % | 10.17 K |

Two things follow. First, **the observational baseline already violates the
ordering**, on roughly 0.25 % of cell-days for `tx < tg` and 0.5 % for `tg < tn`.
This is a property of E-OBS: the three fields are interpolated from station data
independently, by three separate kriging passes, and the ordering that holds at each
station is not preserved on the grid. Both branches inherit it; neither introduces
it. Second, **the perturbation does not make it worse** — the worst case across the
four models is EC-Earth3 raising `tg < tn` from 0.499 % to 0.548 %, a rise of five
hundredths of a percentage point, while HadGEM3-GC3-1MM roughly halves both rates.

The consequences for the indicators were checked rather than assumed. The within-day
integration is **invariant to the sign** of the amplitude `(tx − tn)/2`, because the
sine sampled over a full cycle is symmetric about zero, so HDD and CDD are unaffected
even where the ordering inverts. `heat_daily` and `frost_daily` read `tx` and `tn` on
the affected days, but those are days on which the artefact arises precisely because
the day is cold (for `tx < tg`) or mild (for `tg < tn`), and the indicator is zero
there in any case.

---

## G.6 The thermal-time window

The `crop_gddwin` branch exists because the choice of window is not innocuous, and
under a cooling scenario the window is precisely what moves.

The window is defined as in the thermal-time construction: uncapped growing degree
days, `max(tg − 5, 0)`, are accumulated within the year, and the window is the set of
days whose cumulative total lies between two anchors. The anchors are the mean
*observed* cumulative totals at the end of February (opening) and the end of July
(closing), averaged over 1990–2010, per region.

**The anchors do not move with the scenario.** They represent the crop's thermal
requirement, which does not change because the climate does. This is what makes the
window shift under cooling — and what makes it possible for the closing anchor to
become unreachable, so that the accumulated thermal time never completes the cycle
within the calendar year. That case is recorded as an outcome (the share of
region-years failing to close, reported per bin) rather than patched away.

---

## G.7 Computational organisation

Scenarios are processed in chunks so that the E-OBS record is read once per chunk
rather than once per scenario, with bounded peak memory. Measured cost per scenario
over the full estimation window: **1.60 min** for the energy branch (dominated by the
24-point within-day integration) and **0.22 min** for the crop branch. The complete
set of 40 bin scenarios takes 64 min and 9 min respectively.

---

## G.8 Limitations of the coupling

**A monthly delta shifts the mean of the daily distribution without changing its
shape.** Day-to-day variance, persistence, and the frequency of extremes relative to
the mean are all inherited unchanged from the observed record. Any change in the
*shape* of the daily distribution under AMOC weakening — more variable winters, a
changed frequency of blocking, altered wet-day statistics — is therefore not
represented. What the coupling does capture is the interaction between a shifted
mean and a fixed threshold, which §G.2 shows is itself substantial.

This limitation is a property of the forcing data, not of the method: NAHosMIP
archives monthly means only. It is shared by the ISIMIP branch by construction, since
that branch is deliberately run through the identical machinery so that the two
effects can be compared without a procedural difference between them.

The window question of §G.6 is a second limitation of a different kind, and is not
resolved by either window definition; it is taken up in Appendix I.
