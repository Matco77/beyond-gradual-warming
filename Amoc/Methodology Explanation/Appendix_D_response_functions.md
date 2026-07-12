# Appendix D — Historical Response Functions: Data and Estimation Design

This appendix documents the empirical response functions estimated from
observational-era data — one for crop yield, one for household energy demand —
together with the datasets on which they are built and the reasoning behind each
design choice. These functions are the impact models: the stylised
AMOC-weakening signal constructed in Appendices A and B is superimposed on the
ISIMIP3b baseline and then read *through* the coefficients estimated here to
translate a shifted climate into changes in European crop production and energy
demand. The present appendix concerns the estimation step only; the scenario
replay is documented separately.

## D.1 Design principle

Two modules are estimated, each at the finest spatial resolution at which its
**outcome** is observed:

| Module | Estimation unit | Outcome | Native resolution of the outcome |
|--------|-----------------|---------|----------------------------------|
| Crop   | NUTS3 × crop × year | Yield (t ha⁻¹), CropStatHarm | Sub-national (NUTS3) |
| Energy | Country × fuel × year | Household gas / electricity final use (ktoe), Eurostat | National only |

The asymmetry is dictated by data, not preference: the estimation unit can be no
finer than the outcome. Crop yield is reported sub-nationally, so the crop
response is identified at NUTS3, where within-country cross-sectional weather
variation is large; no sub-national household-energy series exists, so the
energy response can only be identified at the country level. Weather, by
contrast, is always constructed on the fine E-OBS 0.25° (~28 km) grid and then
aggregated **up** to the estimation unit, so spatial detail is preserved as far
as each outcome permits.

The aggregation weight differs by module, matching the physical support of each
channel. Crop weather is **area**-weighted, because crops occupy land; energy
weather is **population**-weighted, because demand follows people. The two
weighting schemes are documented in §D.3 and §D.4.

## D.2 Weather indicators

All indicators are computed **per grid cell, per day, before any spatial or
temporal aggregation**. This ordering is not cosmetic: the transforms are
nonlinear, so for a nonlinear indicator *f*, `E[f(T)] ≠ f(E[T])`, and averaging
temperature before applying a threshold would return a biased value. Let `tg`,
`tx`, `tn` denote the daily mean, maximum, and minimum temperature of a cell,
and `rr` its daily precipitation.

**Crop (growing-season indicators).**

    GDD   = max( min(tg, 28) − 5, 0 )      growing degree-days, warmth capped at 28 °C
    Heat  = max( tx − 28, 0 )              extreme-heat degree-days above 28 °C
    Frost = 1( tn < 0 )                    frost-day indicator
    Precip = rr                            daily precipitation

The 28 °C cap on GDD makes the beneficial-warmth term (GDD) and the
harmful-extreme term (Heat) **non-overlapping**, following the Schlenker–Roberts
degree-day convention: beyond roughly 28 °C additional warmth ceases to help C3
cereals (photorespiration rises, grain-filling shortens), and that range is
carried by Heat rather than double-counted in GDD. The cap binds rarely, because
the daily-mean `tg` seldom exceeds 28 °C, but it is the physiologically correct
and literature-consistent construction.

**Energy (heating and cooling demand).** With hinge functions

    HDD = 1( T < 15 ) · (18 − T)
    CDD = 1( T ≥ 24 ) · (T − 21)

the daily degree-day value is obtained by **integrating the hinge over the
within-day temperature path** — a sine cycle centred on `tg` with amplitude
`(tx − tn)/2` — rather than by evaluating the hinge at the daily mean. The
correction is material for cooling: because the 24 °C cooling threshold sits
*inside* the typical daily swing, evaluating at the mean systematically
undercounts (roughly 34.5 % of summer days straddle 24 °C, understating total
CDD by a factor of about 1.7). It is negligible for crop GDD and mild for HDD,
whose 15 °C threshold is rarely straddled in winter — consistent with the
`E[f(T)] ≠ f(E[T])` bias biting only where a threshold falls within the daily
range. The within-day construction is also what makes the estimated coefficient
**transportable**: the same integrator applies under a shifted scenario climate,
so the coefficient is not tied to the historical position of the temperature
distribution relative to the threshold.

Daily indicators are aggregated to the estimation unit and then summed over the
relevant window: crop over **March–July** (main) or April–August (robustness);
energy HDD over the calendar year and the October–March heating season, energy
CDD over June–August. A single shared construction module builds these variables
identically for estimation and for the later scenario replay, which is the
precondition for a response function to be applied to a counterfactual climate.

## D.3 Crop datasets

**Outcome.** The yield panel is built from CropStatHarm (Ronchetti et al., 2024),
a peer-reviewed sub-national crop dataset already harmonised to NUTS 2016
geometry, so region codes join one-to-one with the weather crosswalk without any
bridging. Area is gap-filled from a second, area-only source (EU-JRC
`cropareaRegional`) where CropStatHarm is missing; CropStatHarm always takes
precedence on overlap, and **no yields are ever fabricated** — only area is
filled. Six crops are retained (Soft, Durum, and Total wheat; Spring, Winter, and
Total barley); because the totals are the component crops summed, any single
regression uses **one set only**, never both, to avoid double-counting.

**Weather.** Each E-OBS cell is mapped to the NUTS3 regions it intersects by
area-overlap weights

    ω(g, r) = Area(g ∩ r) / Σ_h Area(h ∩ r),

computed by exact geometric intersection in the equal-area projection EPSG:3035,
so the weights are true areal shares. Regions smaller than one grid cell receive
a nearest-cell fallback (capped at 40 km); regions off the E-OBS domain
(Canaries, Azores) remain unassigned. Daily cell indicators are combined into an
area-weighted regional mean — renormalised over the cells reporting data on a
given day, so a partially observed day still yields a value — and then summed to
the month and to the crop window.

**Estimation panels.** Two panels are produced. The primary panel keeps
CropStatHarm rows reported at **NUTS3** with observed positive yield (76,979
rows; 15 countries; 816 NUTS3 regions; six crops; 1989–2023). One broken source
cell is dropped by name (ES300 durum wheat 2012, an upstream area error yielding
3,606 t ha⁻¹); no blanket yield cap is imposed. A second, **mixed-level** panel
recovers the eleven countries that report yield only at NUTS2 or NUTS0, by
aggregating the existing NUTS3 weather *up* to the coarser unit by area weights
(no re-reading of E-OBS is required, since NUTS3 nest exactly inside coarser
units and the indicators are area-linear); it spans 86,591 rows and 26 countries.
Each country reports yield at exactly one NUTS level, verified, so there is no
parent–child double-counting. The mixed panel buys wider country **coverage** for
EU-level valuation, but its added units are coarser and contribute less
within-country weather variation — coverage, not identifying power — and it is
therefore reported as a robustness sample rather than the headline.

## D.4 Energy datasets

**Outcome.** Household final use of gas and electricity (ktoe) comes from the
Eurostat energy balance (`nrg_bal_c`). The outcome is expressed **per capita**,
`ln(E / Pop)`, which removes the mechanical scale effect of population and strips
the within-country demographic trend that a country fixed effect — absorbing only
the level — would otherwise leave in the residual.

**Weather.** Population weights are built from the 1 km GISCO Census-2021 grid:
each populated pixel is assigned to its E-OBS cell and its country, giving
`w(g, c) = pop(g) / pop(c)`. Population landing on a data-less E-OBS cell
(coastal or island pixels over sea) is snapped to the nearest valid land cell,
capped at **300 km** — retaining reachable islands (Malta to Sicily, ~110 km)
while dropping off-domain Atlantic islands (Azores, Madeira, ~1,500 km), whose
energy remains in the national outcome but whose weather cannot be fabricated.
Daily HDD/CDD are computed per cell by the within-day construction (§D.2),
population-weighted to the country, and summed into calendar-year HDD, Oct–Mar
heating-season HDD, and JJA CDD. The population weight is **fixed at its 2021
geography** deliberately: it fixes *where people live* as a baseline, while
time-varying population enters the panel separately as an annual control, because
a fixed total would simply be absorbed by the country fixed effect.

**Prices.** Annual household electricity and gas prices (taxes included) splice a
historical Eurostat series (≤2007) to the current series (2007+) at a consistent
consumption band, unit, tax level, and currency, averaging the two semesters to
an annual figure. The splice is required because Eurostat redefined the
consumption bands in 2007.

**Estimation panel.** Merging outcome, weather, population, and prices on a
canonical country identifier — every merge routed through an explicit ISO3 ↔
Eurostat ↔ FAOSTAT crosswalk, never raw names, to avoid coding traps (Greece
`EL ≠ GR`, the United Kingdom `UK ≠ GB`, and similar) — yields the energy panel
(2,924 rows; 29 countries; 1990–2024). Coverage gaps are automatic consequences
of data availability, not choices: Switzerland has no household-energy outcome,
Norway (and Cyprus, Finland, Malta) have no household gas price and so drop from
the gas regression.

## D.5 Response-function specifications

Both modules share one identification logic (Blanc & Schlenker, 2017). With
region and year fixed effects, weather enters each regression as a *deviation
from the local mean*: it is plausibly random, exogenous, and uncorrelated both
with time-invariant confounders (soil, baseline climate) and with common shocks
(prices, technology). Region-specific quadratic time trends further absorb
smooth, region-specific drift (technological yield growth; national efficiency,
income, and fuel-switching paths). What remains identifying is the year-to-year
weather anomaly, which is what a causal weather→outcome coefficient requires.

**Crop response function** (estimated separately for each crop *k*):

    ln(Yield_{r,t}) = β₁·GDD_{r,t} + β₂·Heat_{r,t} + β₃·Frost_{r,t}
                    + β₄·Precip_{r,t} + β₅·Precip²_{r,t}
                    + α_r + λ_t + δ_r·t + γ_r·t² + ε_{r,t}

where `α_r` is the NUTS3 region fixed effect, `λ_t` the year fixed effect, and
`δ_r·t + γ_r·t²` the region-specific quadratic trend. The quadratic precipitation
term captures the concave rainfall response — rain helps, excess hurts. Standard
errors are clustered by country.

**Energy demand function** (estimated separately for each fuel):

    ln( E_{c,t} / Pop_{c,t} ) = β₁·HDD_{c,t} + β₂·CDD_{c,t} + β₃·ln(Price_{c,t})
                              + μ_c + λ_t + δ_c·t + γ_c·t² + ε_{c,t}

with `μ_c` the country fixed effect and `δ_c·t + γ_c·t²` the country-specific
quadratic trend. Electricity uses calendar-year HDD and June–August CDD; gas uses
the October–March heating-season HDD and no CDD. Standard errors are clustered by
country.

The two functions share the **same fixed-effects benchmark** — unit fixed effect,
year fixed effect, and unit-specific quadratic trend — deliberately, so that crop
and energy rest on one common specification. Each is presented as a build-up
(pooled OLS → two-way fixed effects → the benchmark with trends), on the NUTS3
and mixed samples for crop, and pooled OLS → year fixed effects → country trends
→ benchmark for energy. The build-up is reported in full rather than only the
preferred column because it makes visible *why* the benchmark is needed: with
year fixed effects alone, the spatially correlated common European weather (a
cold winter hits every country at once) is absorbed by `λ_t`, leaving HDD
unidentified — insignificant, and wrong-signed for gas. Adding country-specific
trends, while retaining year fixed effects for genuinely common shocks, recovers
a significant and correctly signed heating response (electricity HDD ≈ 6.5 × 10⁻⁵,
gas HDD ≈ 1.4 × 10⁻⁴) and a small but identified cooling response (electricity
CDD ≈ 2 × 10⁻⁴), at the highest fit.

Across the columns of any build-up table the estimation sample is held constant:
every column is refitted on the rows retained by the most saturated
fixed-effects specification (which drops not only rows with missing values but
also singleton fixed-effect groups), so differences across columns reflect the
specification alone and not a moving sample.

## D.6 Inference and robustness

Because the benchmark clusters on country — 15 clusters for crop, 29 for energy —
cluster-robust asymptotics are unreliable, and default (independent) standard
errors are never used. On the benchmark specification the estimation reports, in
addition to the clustered standard errors: (i) two-way clustering by country and
year, to allow for common time shocks; and (ii) a restricted wild-cluster
bootstrap (Rademacher weights, null imposed; Cameron, Gelbach & Miller, 2008;
Roodman et al., 2019), computed on the fixed-effect-partialled model so that each
resample is a matrix operation rather than a full re-fit. A measurement-error
check re-estimates the benchmark after dropping the units whose grid centroid was
snapped to a non-local E-OBS cell (the Ionian NUTS3 for crop; Malta for energy);
these units are few and the estimates are unchanged, so the check is reported as
a defence rather than a live concern.

For crop, a long-difference design (Burke & Emerick, 2016) is estimated as a
complementary **long-run bound**: the long change in mean log yield is regressed
on the long change in mean climate across the 816 NUTS3 (early window 1999–2005
versus late 2016–2022, country fixed effects, area weights), so that
cross-sectional variation in climate *trends* identifies a response that permits
partial adaptation. It is reported as a bound rather than a point estimate,
because a record of about three decades with a modest trend leaves it imprecise.
The same design is not applied to energy, whose 29 cross-sectional units provide
too little identifying variation.

## D.7 Assumptions and limitations

1. **Transient versus permanent climate (the central caveat for the replay).**
   Panel fixed-effects models identify the response to *transient* weather
   anomalies. A *permanent* shift, such as an AMOC-weakening scenario, may elicit
   a different — possibly opposite-signed — long-run response as agents adapt
   (rotations, varieties, capital). The envelope theorem licenses using panel
   variation for a *marginal* change; a large AMOC shift is plausibly
   non-marginal, so the estimates are best read as short-run responses, bounded
   above (for crop) by the long-difference design of §D.6.

2. **Direction of the adaptation gap is scenario-specific.** The sign of the
   short-run/long-run wedge cannot be imported wholesale from the warming
   literature. Auffhammer (2022) finds that panel estimates *understate* the
   long-run electricity impact of *warming* because they miss extensive-margin
   air-conditioning adoption; an AMOC scenario is one of *cooling*, where the
   relevant extensive margin is heating capital and the sign of the gap is a
   distinct argument. The gap is therefore stated in both directions rather than
   assumed.

3. **Precipitation is not temperature.** The panel recovers the long-run response
   for temperature but not for precipitation (Mérel & Gammans, 2021; shown for
   French cereals, present in this panel, by Gammans, Mérel & Ortiz-Bobea, 2017).
   The `Precip` and `Precip²` terms are therefore not given a long-run structural
   reading in the replay, even where the degree-day terms may be.

4. **Measurement-error attenuation.** Fixed effects absorb between-unit variation,
   so interpolation noise in the gridded weather attenuates coefficients toward
   zero, an effect the literature shows can halve estimated impacts (Fisher et
   al., 2012). The snapped-unit robustness check (§D.6) bounds the part of this
   error attributable to non-local grid assignment.

5. **Income confounding is handled by assumption, not by a control.** No GDP
   control is included. The unit-specific quadratic trends already absorb smooth
   national income paths, and identification runs off exogenous weather
   anomalies, which would bias the weather coefficient only if income shocks
   correlated with weather anomalies within a country-year — a second-order
   concern. A GDP control would additionally be a *bad control* for the jointly
   determined energy–income relationship (Auffhammer & Mansur, 2014).

6. **Cooling response is modest.** The identified electricity CDD coefficient is
   small in absolute terms, reflecting genuinely low historical European
   household cooling demand (limited air-conditioning). It matters little for an
   AMOC *cooling* scenario, in which HDD is the load-bearing covariate, but it is
   reported because it is now identified rather than absorbed.

---

**References.**
Auffhammer, M. (2022). Climate adaptive response estimation: Short and long run
impacts of climate change on residential electricity and natural gas
consumption. *Journal of Environmental Economics and Management*, 114, 102669.
Auffhammer, M., & Mansur, E. T. (2014). Measuring climatic impacts on energy
consumption: A review of the empirical literature. *Energy Economics*, 46,
522–530.
Blanc, E., & Schlenker, W. (2017). The use of panel models in assessments of
climate impacts on agriculture. *Review of Environmental Economics and Policy*,
11(2), 258–279.
Burke, M., & Emerick, K. (2016). Adaptation to climate change: Evidence from US
agriculture. *American Economic Journal: Economic Policy*, 8(3), 106–140.
Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). Bootstrap-based
improvements for inference with clustered errors. *Review of Economics and
Statistics*, 90(3), 414–427.
Cornes, R. C., et al. (2018). An ensemble version of the E-OBS temperature and
precipitation data sets. *Journal of Geophysical Research: Atmospheres*, 123,
9391–9409.
Fisher, A. C., Hanemann, W. M., Roberts, M. J., & Schlenker, W. (2012). The
economic impacts of climate change: Comment. *American Economic Review*, 102(7),
3749–3760.
Gammans, M., Mérel, P., & Ortiz-Bobea, A. (2017). Negative impacts of climate
change on cereal yields: Statistical evidence from France. *Environmental
Research Letters*, 12(5), 054007.
Kolstad, C. D., & Moore, F. C. (2020). Estimating the economic impacts of climate
change using weather observations. *Review of Environmental Economics and
Policy*, 14(1), 1–24.
Mérel, P., & Gammans, M. (2021). Climate econometrics: Can the panel approach
account for long-run adaptation? *American Journal of Agricultural Economics*,
103(4), 1207–1238.
Ronchetti, G., et al. (2024). Harmonized European Union subnational crop
statistics. *Earth System Science Data*, 16, 1623–1649.
Roodman, D., MacKinnon, J. G., Nielsen, M. Ø., & Webb, M. D. (2019). Fast and
wild: Bootstrap inference in Stata using boottest. *The Stata Journal*, 19(1),
4–60.
Schlenker, W., & Roberts, M. J. (2009). Nonlinear temperature effects indicate
severe damages to U.S. crop yields under climate change. *Proceedings of the
National Academy of Sciences*, 106(37), 15594–15598.
