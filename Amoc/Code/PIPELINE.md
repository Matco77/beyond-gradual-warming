# Historical Data Pipeline — Methodology

**Thesis:** *Beyond Gradual Warming.* Estimate weather→crop-yield and weather→household-energy
**response functions** from historical data, to be replayed later on ISIMIP and AMOC-stressed
climate scenarios. This document covers the **historical pipeline only** (estimation), not the
scenario replay.

All code is R in `Amoc/Code/`; all data in `Amoc/datasets/`. Outputs are numbered CSVs (`0`–`16`).
Run each with `Rscript <file>`.

---

## 1. Design in one page

Two modules, each estimating a response function at the finest resolution its **outcome** exists at:

| Module | Estimation unit | Outcome | Why this unit |
|--------|-----------------|---------|---------------|
| **Crop** | NUTS3 × crop × year | CropStatHarm yield (t/ha) | Yield is reported subnationally → NUTS3 gives large cross-sectional weather variation |
| **Energy** | Country × fuel × year | Eurostat household gas/electricity (ktoe) | No subnational household-energy series exists → country is the finest possible |

The asymmetry is deliberate and dictated by data, not preference: the estimation unit can be no
finer than the outcome. Weather is always built on the fine **E-OBS 0.25° (~28 km) grid** and then
aggregated up to the estimation unit.

**Aggregation weights differ by module, matching what each channel is about:**
- Crop weather → **area**-weighted (crops occupy land).
- Energy weather → **population**-weighted (energy demand follows people).

**Core econometric principle (Blanc & Schlenker 2017, *REEP*):** with region + year fixed effects,
weather enters as *deviations from the local mean* — random, exogenous, uncorrelated with
time-invariant confounders (soil, baseline climate) and with common shocks (prices, technology).
This is what identifies a **causal** weather→outcome relationship.

**A rule applied throughout:** nonlinear indicators (degree-days, frost) are computed **per grid
cell per day first, then aggregated**. Because the transforms are nonlinear, `E[f(T)] ≠ f(E[T])` —
averaging temperature before applying the threshold would give the wrong answer.

---

## 2. Pipeline map and data lineage

```
RAW INPUTS                          SCRIPT                              OUTPUTS
─────────────────────────────────────────────────────────────────────────────────────
CropStatHarm CSV ─┐
cropareaRegional ─┼──► agriculturegeografical.R ──► 6.CropStatHarm_prepared.csv
NUTS 2016 shp ────┘

E-OBS (tg,tx,tn,rr) ─┐
NUTS 2016 shp ───────┼► cropweather_eobs_nuts3.R ──► 6.eobs_to_nuts_crosswalk.csv
                     │                                7.eobs_nuts3_crop_weather_monthly.csv
                     │                                8.eobs_nuts3_crop_weather_window.csv
6.prepared + 8 ──────┼► crop_panel_nuts3.R ─────────► 9.crop_panel_nuts3_estimation.csv   [CROP PANEL, 15 countries]
6.prepared + 8 ──────┴► crop_panel_mixed.R ─────────► 10.crop_panel_mixed_estimation.csv  [CROP PANEL, 26 countries]

E-OBS + GISCO pop grid ─► energyweather_eobs_country.R ► 12.eobs_population_weights.csv
     (uses weather_indicators.R)                        13.eobs_country_energy_weather_weighted.csv
Eurostat price TSVs ────► energy_prices.R ─────────────► 14.energy_prices_prepared.csv
nrg_bal_c + 13/14 + pop ► energy_panel.R ──────────────► 0.country_crosswalk.csv
                                                         15.population_country_year.csv
                                                         16.energy_panel_estimation.csv   [ENERGY PANEL]

9,10,16 ──────────────► regression.R ──────────────────► Amoc/results/*.txt              [ESTIMATES]
9,10,16 ──────────────► covariate_mappings.R ──────────► Amoc/figures/*.png              [MAPS]
```

`weather_indicators.R` is a **shared module** (not run directly) holding the indicator formulas, so
estimation and the later scenario replay build identical variables.

### 2.1 Starting datasets (raw inputs)
Complete paths, relative to the project root `~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/`.
Every generated file below lives in `Amoc/datasets/`.

| Starting dataset (complete path) | Source | Contents |
|---|---|---|
| `Amoc/datasets/CropStatHarm_sub-national_16234acd-d17d-49de-bcb2-908a366a12be_2025.v01_39.csv` | Ronchetti et al. 2024, ESSD | Subnational crop area / production / yield, harmonised to NUTS 2016 |
| `Amoc/datasets/cropareaRegional_b2ba07b7-73da-4d7f-96f6-4928b2f2b26e_2022.01_31.csv` | EU-JRC | Supplementary regional crop **area** (gap-fill only) |
| `Amoc/datasets/ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip` | Eurostat GISCO | NUTS 2016 region polygons, all levels, WGS84 |
| `Amoc/datasets/EOBS/tg_ens_mean_0.25deg_reg_v31.0e.nc` | E-OBS v31.0e | Daily **mean** temperature, 0.25° grid |
| `Amoc/datasets/EOBS/tx_ens_mean_0.25deg_reg_v31.0e.nc` | E-OBS v31.0e | Daily **max** temperature |
| `Amoc/datasets/EOBS/tn_ens_mean_0.25deg_reg_v31.0e.nc` | E-OBS v31.0e | Daily **min** temperature |
| `Amoc/datasets/EOBS/rr_ens_mean_0.25deg_reg_v31.0e.nc` | E-OBS v31.0e | Daily **precipitation** |
| `Amoc/datasets/GISCO_population_grid/ESTAT_OBS-VALUE-T_2021_V2.tiff` | Eurostat GISCO | 1 km Census-2021 population grid |
| `Amoc/datasets/PricesEectricitiesandgas/estat_nrg_pc_204_h.tsv` | Eurostat `nrg_pc_204_h` | Household **electricity** price, historical (≤2007) |
| `Amoc/datasets/PricesEectricitiesandgas/estat_nrg_pc_204.tsv` | Eurostat `nrg_pc_204` | Household **electricity** price, current (2007+) |
| `Amoc/datasets/PricesEectricitiesandgas/estat_nrg_pc_202_h.tsv` | Eurostat `nrg_pc_202_h` | Household **gas** price, historical (≤2007) |
| `Amoc/datasets/PricesEectricitiesandgas/estat_nrg_pc_202.tsv` | Eurostat `nrg_pc_202` | Household **gas** price, current (2007+) |
| `Amoc/datasets/2.nrg_bal_c_prepared.csv` | Eurostat `nrg_bal_c` (pre-prepared upstream) | Household gas/electricity final use (ktoe) — no regenerator, treat as raw |
| `Amoc/datasets/totalpopulation.csv` | Eurostat DEMO_GIND | Annual country population |

(`estat_nrg_pc_203*/205*.tsv` also exist in the price folder but are **not** used — only 202 (gas)
and 204 (electricity).)

### 2.2 Derived datasets — inputs → outputs
Numbered outputs are all written to `Amoc/datasets/`; here referenced by filename.

| Script | Reads | Writes |
|---|---|---|
| `agriculturegeografical.R` | `CropStatHarm_sub-national_…_2025.v01_39.csv`, `cropareaRegional_…_2022.01_31.csv`, `ref-nuts-2016-03m.shp/NUTS_RG_03M_2016_4326.shp.zip` | `6.CropStatHarm_prepared.csv` |
| `cropweather_eobs_nuts3.R` | `EOBS/{tg,tx,tn,rr}_ens_mean_0.25deg_reg_v31.0e.nc`, `ref-nuts-2016-03m.shp/…shp.zip` | `6.eobs_to_nuts_crosswalk.csv`, `7.eobs_nuts3_crop_weather_monthly.csv`, `8.eobs_nuts3_crop_weather_window.csv` |
| `crop_panel_nuts3.R` | `6.CropStatHarm_prepared.csv`, `8.eobs_nuts3_crop_weather_window.csv` | `9.crop_panel_nuts3_estimation.csv` |
| `crop_panel_mixed.R` | `6.CropStatHarm_prepared.csv`, `8.eobs_nuts3_crop_weather_window.csv`, `6.eobs_to_nuts_crosswalk.csv` (region areas) | `10.crop_panel_mixed_estimation.csv` |
| `energyweather_eobs_country.R` | `EOBS/{tg,tx,tn}_ens_mean_0.25deg_reg_v31.0e.nc`, `GISCO_population_grid/ESTAT_OBS-VALUE-T_2021_V2.tiff`, `ref-nuts-2016-03m.shp/…shp.zip`, `weather_indicators.R` (sourced) | `12.eobs_population_weights.csv`, `13.eobs_country_energy_weather_weighted.csv`, `13.eobs_country_hdd_cdd_daily.csv` |
| `energy_prices.R` | `PricesEectricitiesandgas/estat_nrg_pc_{204_h,204,202_h,202}.tsv`, `0.country_crosswalk.csv` | `14.energy_prices_prepared.csv` |
| `energy_panel.R` | `2.nrg_bal_c_prepared.csv`, `13.eobs_country_energy_weather_weighted.csv`, `14.energy_prices_prepared.csv`, `totalpopulation.csv` | `0.country_crosswalk.csv`, `15.population_country_year.csv`, `16.energy_panel_estimation.csv` |
| `regression.R` | `9.crop_panel_nuts3_estimation.csv`, `10.crop_panel_mixed_estimation.csv`, `16.energy_panel_estimation.csv`, `weather_indicators.R` | `Amoc/results/*.txt` |

**Energy run order (one soft cycle):** `energy_panel.R` writes `0.country_crosswalk.csv`, which
`energy_prices.R` reads; `energy_panel.R` in turn reads `energy_prices`' output `14`. So run
`energy_panel.R` → `energy_prices.R` → `energy_panel.R` once more so the prices land in panel `16`.

---

## 3. Weather indicators (definitions)

Per grid cell *g*, day *d*. `tg`=daily mean, `tx`=daily max, `tn`=daily min, `rr`=precipitation.

**Crop (growing-season, proposal slide 14):**
| Indicator | Formula | Meaning |
|-----------|---------|---------|
| GDD | `max(min(tg,28) − 5, 0)` | Growing degree-days; warmth **capped at 28 °C** (see §7) |
| Heat | `max(tx − 28, 0)` | Extreme-heat degree-days above 28 °C |
| Frost | `1(tn < 0)` | Frost day |
| Precip | `rr` | Daily precipitation |

**Energy (proposal slide 19):** with hinges `HDD=1(T<15)(18−T)`, `CDD=1(T≥24)(T−21)`, integrated
over the **within-day** temperature path (§7), not evaluated at the daily mean.

Daily indicators are aggregated to the estimation unit, then summed over the relevant window:
crop = **March–July** (main) / April–August (robustness); energy HDD = calendar-year and Oct–Mar
heating season, CDD = June–August.

---

## 4. Crop module

### 4.1 `agriculturegeografical.R` — crop outcome
**What:** builds the region × crop × year yield panel.
**How:**
1. Read CropStatHarm sub-national CSV; the `UNIT` field encodes the variable (`ha`=area,
   `t`=production, `t/ha`=yield) — pivot these into three columns, one row per region-crop-year.
2. Read `cropareaRegional` (a second, area-only source) and use it to **gap-fill area** where
   CropStatHarm is missing; CropStatHarm always wins on overlap (it is the peer-reviewed product).
   Area only — **no fabricated yields**.
3. Join to the NUTS 2016 shapefile; assert every region has a polygon.

**Output:** `6.CropStatHarm_prepared.csv` (region, level, crop, year, area, production, yield).
**Why:** CropStatHarm is harmonised to NUTS 2016 geometry, so region codes join 1:1 to the weather
crosswalk with no bridging. Six crops are kept (Soft/Durum/Total wheat, Spring/Winter/Total barley);
totals are components summed — the regression must use **one set**, never both (double-counting).

### 4.2 `cropweather_eobs_nuts3.R` — crop weather
**What:** turns the raw E-OBS grid into NUTS3 growing-season weather.

**Step A/B — spatial crosswalk (once).** Build E-OBS cell → NUTS3 area-overlap weights
`ω_{g,r} = Area(g ∩ r) / Σ Area(h ∩ r)`, computed by `st_intersection` in **EPSG:3035**
(equal-area, so weights are true areal shares). Regions smaller than a grid cell get a
nearest-cell fallback (≤40 km); regions off the E-OBS domain (Canaries, Azores) stay unassigned.
→ `6.eobs_to_nuts_crosswalk.csv`.

**Step C/D — daily indicators → NUTS3-month.** For each day, compute the four cell indicators,
take the area-weighted regional mean (renormalised over cells with data that day, so a
partially-observed day still yields a value; a fully-unobserved day is NA), and sum over the month.
→ `7.eobs_nuts3_crop_weather_monthly.csv` (1,356,300 rows = 1,507 regions × 75 years × 12 months).

**Step E — window sums.** Sum monthly values into the crop windows (Mar–Jul, Apr–Aug).
→ `8.eobs_nuts3_crop_weather_window.csv`.

**Why these choices:**
- *E-OBS read via `ncdf4` hyperslabs*, not `terra::extract` — measured ~400× faster on this file.
- *Nonlinear-before-aggregation* — indicators are computed per cell-day, never on averaged temperature.
- *Monthly then window* — the monthly file is a flexible superset; windows are month sums, so
  alternative windows need no re-read.

### 4.3 `crop_panel_nuts3.R` — crop estimation panel (15 countries)
**What:** joins NUTS3 yield to NUTS3 weather.
**How:** keep CropStatHarm rows at **NUTS3 only**, with observed `yield > 0`; drop one broken source
cell (ES300 durum 2012: 4 ha vs 14,426 t → 3,606 t/ha, an upstream area error). Join weather on
region × year (weather is crop-agnostic; the crop-specific response comes from the regression).
**Output:** `9.crop_panel_nuts3_estimation.csv` — **76,979 rows, 15 countries, 816 NUTS3 regions,
6 crops, 1989–2023.**
**Why NUTS3-only here:** the weather crosswalk exists only at NUTS3, so coarser CropStatHarm rows
(NUTS2/NUTS0) have no weather to join and are dropped.

### 4.4 `crop_panel_mixed.R` — mixed-level panel (26 countries)
**What:** recovers the 11 countries that report yield only at NUTS2 or NUTS0.
**How:** each country reports yield at **exactly one** NUTS level (verified — so no parent/child
double-counting, no dedup needed). Aggregate the existing NUTS3 weather **up** to NUTS2/NUTS0 by
area weights — no E-OBS re-read, since NUTS3 nest perfectly inside coarser units and the indicators
are area-linear. Join each yield row to weather at its own level.
**Output:** `10.crop_panel_mixed_estimation.csv` — **86,591 rows, 26 countries** (NUTS3 block
identical to panel 9).
**Why:** wider country coverage for EU-wide impact valuation. Trade-off: the added units are coarser,
so they contribute less within-country weather variation — they buy **coverage, not power**.

---

## 5. Energy module

### 5.1 `energyweather_eobs_country.R` — energy weather
**What:** population-weighted country HDD/CDD from E-OBS.

**Step 1 — population weights.** From the 1 km GISCO Census-2021 grid, assign each populated pixel
to its E-OBS cell and its country, giving `w_pop_{g,c} = pop(g) / pop(c)`. Population that lands on a
data-less E-OBS cell (coastal/island pixels over sea) is **snapped to the nearest valid land cell**,
capped at **300 km**: this keeps reachable islands (Malta→Sicily ~110 km) but drops off-domain
Atlantic islands (Azores/Madeira ~1,500 km), whose energy stays in the national outcome but whose
weather cannot be fabricated. → `12.eobs_population_weights.csv`.

**Step 2 — country HDD/CDD.** Compute daily HDD/CDD per cell using the **within-day** construction
(§7), population-weight to country, and sum into calendar-year HDD, Oct–Mar heating-season HDD
(seasons with ≥150 days), and JJA CDD. → `13.eobs_country_energy_weather_weighted.csv`
(2,475 rows = 33 countries × 75 years) and the daily file `13.eobs_country_hdd_cdd_daily.csv`.

**Why fixed 2021 population weights:** the proposal mandates a fixed baseline geography (where
people live); time-varying population enters separately as a control (§5.3), because a fixed total
would be absorbed by the country fixed effect.

### 5.2 `energy_prices.R` — fuel prices
**What:** annual household electricity and gas prices, taxes included.
**How:** read four Eurostat bi-annual TSVs. Each fuel splices a **historical** series (`_h`, ≤2007)
with the current series (2007+) at a fixed consumption band, unit, tax level (`I_TAX`), and currency
(EUR); average the S1/S2 semesters to annual.
**Output:** `14.energy_prices_prepared.csv`, 1985–2025.
**Why the splice:** Eurostat changed band definitions in 2007; the `_h`↔current splice at a
consistent band gives one continuous series (a small 2007 level jump is possible and is checked).

### 5.3 `energy_panel.R` — energy estimation panel
**What:** assembles the country × fuel × year panel.
**How:**
1. Build `0.country_crosswalk.csv` mapping ISO3 ↔ Eurostat/E-OBS 2-letter ↔ FAOSTAT name (handles
   traps: Greece `EL≠GR`, UK `UK≠GB`, "Germany including former GDR", etc.) — **every merge uses the
   crosswalk, never raw names.**
2. Merge household energy (`nrg_bal_c`) + weather (`13`) + annual population (`15`, from
   `totalpopulation.csv`) + prices (`14`), all on the canonical `country_id`.
3. Form the per-capita outcome `ln_energy_pc = ln_energy − ln_population`.

**Output:** `16.energy_panel_estimation.csv` — **2,924 rows, 29 countries, 1990–2024.**
**Why per-capita:** it removes the mechanical scale effect of population and strips the
within-country demographic trend that country fixed effects (which absorb only the *level*) leave in.

**Coverage (by data availability, not choice):**
| | Norway | Switzerland | UK |
|--|:--:|:--:|:--:|
| Electricity | ✓ | ✗ (no outcome) | ✓ |
| Gas | ✗ (no gas price) | ✗ | ✓ |

Also CY/FI/MT drop from gas (no household gas price). These are automatic, from `NA` price rows.

---

## 6. `weather_indicators.R` — shared construction module
Holds the indicator formulas in one place so **estimation and scenario replay build identical
variables** — the requirement for a response function to be transportable to a counterfactual
climate. Contains: the `within_day` integrator, the HDD/CDD hinges, the crop daily indicators, and a
self-check (within-day reduces to the daily-mean value when `tx = tn`; a threshold-straddling day
gives positive CDD where the daily mean gives zero).

---

## 7. Two indicator-construction refinements (and their rationale)

**(a) Capped GDD (`min(tg, 28)`).** Beyond ~28 °C, warmth stops helping C3 cereals (photorespiration
rises, grain-filling shortens) — the harmful range is captured separately by the Heat term. Capping
makes GDD (beneficial warmth) and Heat (harmful extremes) **non-overlapping**, following the
Schlenker–Roberts degree-day convention. Empirically the cap binds rarely (daily-mean `tg` seldom
exceeds 28 °C), so it changes little, but it is the physiologically correct, literature-consistent
design.

**(b) Within-day HDD/CDD (energy).** A day with mean 22 °C but a 31 °C afternoon has real cooling
demand, yet `CDD(mean) = 0`. Because the CDD threshold (24 °C) sits **inside** the daily temperature
swing, the daily mean systematically undercounts — 34.5 % of summer days straddle 24 °C, and total
CDD is understated ~1.7×. We therefore integrate the hinge over the within-day cycle (a sine centred
on `tg` with amplitude `(tx−tn)/2`). This matters strongly for CDD, mildly for HDD (15 °C threshold
is rarely straddled in winter), and negligibly for crop GDD — consistent with `E[f(T)] ≠ f(E[T])`
biting only where a threshold falls within the daily range. It is also what makes the response
function **transportable**: the same construction applies under a shifted scenario climate, so the
estimated coefficient is not tied to the historical position of the temperature distribution
relative to 24 °C.

---

## 8. `regression.R` — estimation

### 8.1 Crop response function (per crop *k*, estimated separately)
$$
\ln(\text{Yield}_{r,t}) = \beta_{1}\text{GDD}_{r,t}+\beta_{2}\text{Heat}_{r,t}+\beta_{3}\text{Frost}_{r,t}
+\beta_{4}\text{Precip}_{r,t}+\beta_{5}\text{Precip}^2_{r,t}+\alpha_{r}+\lambda_{t}+\delta_{r}t+\gamma_{r}t^2+\varepsilon_{r,t}
$$
- Mar–Jul window; region FE `α_r`, year FE `λ_t`, **region-specific quadratic time trends** `δ_r t + γ_r t²`
  (absorb region-specific technology yield growth); SE clustered by country.
- Presented as a **build-up**: pooled OLS → two-way FE → +region trends (the Blanc & Schlenker
  benchmark), on **both** the NUTS3 (15-country) and mixed (26-country) samples.
- `Precip²` gives the concave rainfall response (rain helps, excess hurts).

### 8.2 Energy demand function (per fuel, per-capita)
$$
\ln\!\frac{E_{c,t}}{\text{Pop}_{c,t}} = \beta_1\text{HDD}_{c,t}+\beta_2\text{CDD}^{\text{JJA}}_{c,t}\,(\text{elec})+\beta_3\ln(\text{Price}_{c,t})+\mu_c+\lambda_t+\delta_c t+\gamma_c t^2+\varepsilon_{c,t}
$$
- Electricity uses calendar HDD + JJA CDD; gas uses Oct–Mar HDD (no CDD).
- **Benchmark spec (headline) = the SAME structure as the crop model:** country FE `μ_c` + year FE
  `λ_t` + **country-specific quadratic trends** `δ_c t + γ_c t²` — all three (Blanc & Schlenker,
  Appendix A). SE clustered by country.
- Build-up: pooled OLS → year FE → country trends → **benchmark (main)**.
- **Why the benchmark, not year-FE-only:** year FE alone absorbs the spatially-correlated common
  European weather (a cold winter hits every country), leaving HDD **unidentified** — insignificant,
  and wrong-signed for gas. Adding country-specific trends (absorbing smooth national drift:
  efficiency, income, fuel-switching) while keeping year FE for common shocks **recovers a
  significant, correctly-signed** HDD (electricity 6.5e-5\*\*\*, gas 1.4e-4\*) and CDD (electricity
  2e-4\*), and the highest fit. The full build-up is reported so the year-FE absorbed-variation
  problem is visible rather than hidden. This makes crop and energy use one common benchmark.

### 8.3 A table convention
`common_obs()` refits every column of a build-up table on the rows the most-saturated
fixed-effects spec retains (fixed effects drop not only `NA` rows but **singleton** groups). This
holds the estimation sample constant across columns, so differences reflect the specification, not a
moving sample.

---

## 9. Visualization (not estimation)
- `covariate_mappings.R` → choropleths of every crop and energy covariate (`Amoc/figures/`).
- `energy_climate_nuts3.R` → NUTS3 population-weighted HDD/CDD maps, illustrating within-country
  climate variation the country panel necessarily collapses. **Visualization only** — energy is not
  re-estimated at NUTS3 (no subnational outcome).

---

## 10. Methodological caveats (to state in the write-up)

1. **Short-run vs long-run / adaptation (Blanc & Schlenker; Kolstad & Moore 2020).** Panel models
   identify the response to *transient* weather shocks. A *permanent* climate shift (an AMOC
   scenario) may elicit a different, possibly opposite-sign long-run response as agents adapt. The
   envelope theorem licenses using panel variation for *marginal* change; a large AMOC shift is
   plausibly non-marginal. **This is the central caveat for the scenario replay.** Two refinements:
   - *Crop long-run bound.* A long-difference design (Burke & Emerick 2016) exploits cross-sectional
     variation in climate *trends* across the 816 NUTS3 to bound the longer-run (partly adapted)
     response. Reported as a *bound*, not a point estimate — a 35-year record with a modest trend
     leaves it imprecise. Not run for energy (only 29 cross-sectional units → uninformative).
   - *Direction of the short-run bias is scenario-specific.* Auffhammer (2022) shows panel estimates
     *understate* long-run electricity impacts of **warming** because they miss extensive-margin AC
     adoption. The AMOC scenario is **cooling**: the relevant extensive margin is heating capital, so
     the sign of the adaptation gap is a *different* argument and must not borrow the warming sign.
     State the gap in both directions.
2. **Measurement-error attenuation.** Fixed effects absorb the between-region variation, so
   interpolation noise in the weather data (E-OBS) attenuates coefficients toward zero; the
   literature shows estimated impacts can halve with noisier weather data (Fisher et al. 2012). Units
   whose centroid was *snapped* to a non-local E-OBS grid cell inject this error into specific
   observations; `regression.R` reports a with/without-snapped robustness column (crop: drop the
   Ionian `EL62*` NUTS3; energy: drop Malta `MLT`, snapped to Sicily). Expected to move nothing (3-4
   tiny crop islands; 1 energy country) — reported as defense, not a live threat.
3. **Cooling response is modest but identified.** Under the benchmark spec the electricity CDD
   coefficient is significant (≈2e-4\*), once (i) the within-day construction fixes the daily-mean
   undercount and (ii) country-specific trends stop year FE from absorbing the common weather. The
   *marginal effect* stays small — European household cooling demand is genuinely low historically
   (limited air-conditioning) — but it is no longer an unidentified/insignificant term. It also
   matters little for an AMOC *cooling* scenario, where HDD is the load-bearing covariate.
4. **Explored and parked:** a thermal-time (GDD-accumulation) growing window was tested; it fits the
   historical data slightly *worse* (summing GDD over a GDD-defined window is near-circular and
   introduces variable window length), so the fixed Mar–Jul window is retained for estimation. The
   thermal window's value is for scenario timing, not historical fit.
5. **Precipitation is not temperature (Mérel & Gammans 2021).** The panel recovers the long-run
   response for *temperature* but **not** for *precipitation* — shown for French cereals in Gammans,
   Mérel & Ortiz-Bobea (2017), a country contained in the crop panel. So in the scenario replay the
   `Precip`/`Precip²` terms must **not** be given a long-run structural reading, even where the
   degree-day terms may be. Do not treat all weather covariates symmetrically.
6. **Few-cluster inference.** The benchmark clusters on country (crop 15, energy 29). Cluster-robust
   asymptotics are unreliable at that count, so `regression.R` reports, on the benchmark spec, (a)
   two-way clustering (country + year) for common time shocks and (b) a restricted wild-cluster
   bootstrap (Rademacher weights, null imposed; Cameron, Gelbach & Miller 2008; Roodman et al. 2019),
   computed on the FE-partialled model. Default (iid) SEs are never used.
7. **Income/GDP confounding — assumption, not omission.** No GDP control is included, by design. The
   benchmark's country-specific quadratic trends (`δ_c t + γ_c t²`) already absorb *smooth* country
   income paths, and FE identification runs off exogenous weather *anomalies*, which bias the weather
   coefficient only if income shocks correlate with weather anomalies within country-year (second
   order). A GDP control is also a *bad control* for jointly-determined energy–GDP (Auffhammer &
   Mansur 2014). Optional robustness if a GDP series is added; not a default.

---

## 11. Deliberate choices — do **not** "fix" these
- Crop at NUTS3, energy at country (outcome resolution, §1).
- Crop weather area-weighted, energy population-weighted (§1).
- Fixed 2021 population grid for weights; time-varying population as a separate control (§5.1, §5.3).
- Switzerland absent from energy estimation (no outcome); Norway absent from gas (no price).
- Six crops kept in the panel for flexibility — a regression must select one set (§4.1).
- One broken source cell dropped by name (ES300 durum 2012), no blanket yield cap (§4.3).
