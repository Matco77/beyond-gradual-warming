# Review guide — historical data pipeline (crop + energy modules)

Purpose: hand this to a reviewing agent to hunt for **mistakes, double-counting, unit errors,
join/coverage problems, and silent NA**. This document says what each script does, which
choices are *deliberate* (do not "fix" them), and the *specific* things to verify.

Thesis: "Beyond Gradual Warming" — historical weather→crop-yield and weather→household-energy
response functions, later replayed on ISIMIP + AMOC-stressed weather. Spec in
`Presentations and Proposal/tesi_proposal 2.pdf`. Only the HISTORICAL pipeline is built; the
scenario replay and valuation are NOT built yet.

All scripts are R, in `Amoc/Code/`, reading/writing `Amoc/datasets/`. They read E-OBS via
`ncdf4` hyperslabs (terra per-cell indexing is ~400× slower here — do not "simplify" back to
terra::extract). Spatial overlays in EPSG:3035 (equal-area). Run each with `Rscript <file>`.

---

## Pipeline map (script → key inputs → outputs)

CROP module (estimation unit = NUTS3 × crop × year):
- `agriculturegeografical.R` → `6.CropStatHarm_prepared.csv` (region×crop×year: area/production/yield,
  CropStatHarm harmonised to NUTS 2016 + cropareaRegional area gap-fill) and
  `6.NUTS2016_cropregions.gpkg`.
- `cropweather_eobs_nuts3.R` → `6.eobs_to_nuts_crosswalk.csv` (E-OBS cell→NUTS3 area weights ω),
  `7.eobs_nuts3_crop_weather_monthly.csv` (GDD/Heat/Frost/Precip per NUTS3×year×month),
  `8.eobs_nuts3_crop_weather_window.csv` (summed to crop windows Mar-Jul main, Apr-Aug robustness).
- `crop_panel_nuts3.R` → `9.crop_panel_nuts3_estimation.csv` (yield ⋈ weather, the estimation panel).

ENERGY module (estimation unit = country × fuel × year):
- `energyweather_eobs_country.R` → `12.eobs_population_weights.csv` (E-OBS cell→country pop weights),
  `13.eobs_country_energy_weather_weighted.csv` (annual HDD calendar+OctMar, CDD JJA),
  `13.eobs_country_hdd_cdd_daily.csv` (daily HDD/CDD).
- `energy_prices.R` → `14.energy_prices_prepared.csv` (household gas + electricity price, annual).
- `energy_panel.R` → `0.country_crosswalk.csv`, `15.population_country_year.csv`,
  `16.energy_panel_estimation.csv` (energy outcome ⋈ weather ⋈ population ⋈ price).

Upstream (not written by us, pre-existing): `1.FAOSTAT_prepared.csv`, `2.nrg_bal_c_prepared.csv`,
`3.EOBS_prepared.csv`, `4.chdd_a_prepared.csv`, `5.chdd_m_prepared.csv`. `Dataset_Preparation.R`
is an older scratch/prep file — section 6 there is superseded by `agriculturegeografical.R`.

Formulas (proposal slides 13-14, 18-19), verify the code matches EXACTLY:
- Crop daily per cell: `GDD=max(min(tg,28)-5,0)` (CAPPED at 28C, deliberate deviation from the
  proposal's uncapped GDD - beyond 28 warmth stops helping C3 cereals (photorespiration,
  accelerated grain-fill) and is captured instead by Heat; Schlenker-Roberts convention, keeps
  GDD and Heat non-overlapping), `Heat=max(tx-28,0)`, `Frost=1(tn<0)`, `Precip=rr`.
  NUTS3 aggregation `Indicator_{r,m} = sum_days sum_cells ω_{g,r}·Indicator_{g,d}` (area-weighted
  mean per day, summed over days). **Temperature must NOT be averaged before the nonlinear step.**
- Energy daily per cell: `HDD=1(T<15)(18-T)`, `CDD=1(T>=24)(T-21)`, T=tg. Country aggregation
  population-weighted. HDD^Calendar (Jan-Dec), HDD^OctMar (Oct(t-1)-Mar t), CDD^JJA (Jun-Aug).

---

## Deliberate decisions — do NOT flag these as bugs

1. **Crop estimated at NUTS3, energy at country.** Intentional asymmetry: CropStatHarm gives
   subnational yield (→ NUTS3 power); Eurostat energy is country-only (no subnational outcome
   exists → NUTS3 would be pointless, weights compose to the same country number).
2. **Weather is crop-agnostic** (same window both crops). Crop-specific response comes from
   crop-specific coefficients in the regression, not from the weather file. Correct per slide 14.
3. **Monthly crop file (7) then window sums (8).** Monthly is a superset; windows are month sums.
4. **Population weights fixed at 2021** (Census grid). Proposal mandates fixed baseline geography.
   A time-varying *annual* population enters only as the control `15` (time-varying, so NOT
   absorbed by country FE — a fixed 2021 total WOULD be, which is why the grid can't be the control).
5. **Switzerland absent from energy panel** — Eurostat `nrg_bal_c` has no CH household energy.
   Not a merge bug.
6. **Malta HDD/CDD = Sicily proxy** (nearest reliable land cell); a few residual NA years. Malta
   has no E-OBS land cell. Deliberate, flagged in `energyweather_eobs_country.R`.
7. **`_c` price component tables intentionally unused** (summing ≠ published total; start 2017).
   Prices come from `_h` (≤2007) spliced with plain `nrg_pc_202/204` (2007+).

---

## CHECK THIS — specific risks to verify

### A. Double counting (highest priority)
1. **Crop totals vs components.** `9.crop_panel_nuts3_estimation.csv` contains all 6 crops:
   `Total wheat = Soft + Durum`, `Total barley = Spring + Winter`. A regression must use ONE set
   (totals OR components), never both. Confirm any downstream regression filters crop accordingly.
   The panel itself SHOULD keep all 6 (flexibility) — verify the header comment warns about this.
2. **NUTS level mixing / parent-child.** Weather crosswalk (`6.eobs_to_nuts_crosswalk`) is NUTS3
   only. Crop panel filters `nuts_level==3`, so NUTS0/NUTS2 CropStatHarm rows are dropped (no
   NUTS3 weather). VERIFIED during build: `9` is NUTS3-only (all NUTS_ID length 5, no NUTS2 parent
   rows) — no level-mixing double count. Note `6.CropStatHarm_regression_panel.csv` is written by
   `agriculturegeografical.R` but read by NOTHING (dead output; `crop_panel_nuts3.R` uses
   `6.CropStatHarm_prepared.csv`). Safe to delete, or wire it if its parent/child dedup was intended.
3. **cropareaRegional gap-fill.** `6.CropStatHarm_prepared.csv` merged a second area source
   (`area_source`, `area_ha_alt`). It added AREA only, not yield. Verify `crop_panel_nuts3.R`'s
   `yield_t_ha>0` filter keeps only genuine CropStatHarm yield rows (no fabricated yields).

### B. Weights and aggregation
4. **ω sums to 1 per NUTS3, w_pop sums to 1 per country.** Both scripts assert this; re-verify on
   the output CSVs (group by region/country, sum weight ≈ 1). Watch NA-renormalisation: daily
   aggregation divides by present-weight, so a partially-observed day still yields a value —
   confirm this is intended and that fully-unobserved days give NA, not 0.
5. **Area weighting in EPSG:3035.** Crop crosswalk uses `st_intersection` area in 3035. Confirm
   not accidentally computed in lon/lat.
6. **Nearest-cell fallbacks.** Crop crosswalk snaps ≤40 km island regions (quality_flag=
   `nearest_fallback`); energy snaps population off data-less cells to nearest reliable land cell
   (validity = data on every day of a recent sample). Verify these don't silently assign a far/
   wrong cell (check distances; Malta→Sicily is ~80 km and IS a known approximation).

### C. Joins, units, coverage
7. **Country key.** `0.country_crosswalk.csv` maps ISO3 ↔ Eurostat/E-OBS 2-letter ↔ FAOSTAT name.
   Traps: Greece EL≠GR, UK≠GB, Germany in DEMO_GIND = "Germany including former GDR", Netherlands
   FAOSTAT = "Netherlands (Kingdom of the)". Verify every merge uses the crosswalk, not raw names.
8. **Units per fuel.** Electricity price EUR/kWh, gas price EUR/GJ — different units, fine because
   regressions are per-fuel. Energy outcome ktoe. Verify no cross-fuel unit mixing.
9. **Price band consistency across the 2007 splice.** Electricity: old band Dc(3500 kWh, code
   4161150) spliced to new DC(KWH2500-4999); gas: old D2(4141100) to new GJ20-199. Small level
   jump at 2007 possible — check for a discontinuity. Both use I_TAX (incl. taxes), EUR.
10. **Year coverage.** Crop yield 1989-2023; crop weather 1950-2024; energy 1990-2024; prices
    1990-2024 (electricity full, gas full after the 202 splice); population 1960-2026. Estimation
    windows = overlaps. Verify no accidental year drops or off-by-one in the OctMar cross-year
    season (season year = calendar year + (month>=10); first/last incomplete seasons dropped at
    <150 days).

### D. Silent NA / sanity
11. Re-run each script; confirm exit 0 and the reported row counts match (crop panel 76,980;
    energy panel 2924; weather monthly 1,356,300). Spot-check a warm vs cold region/country for
    plausible GDD/frost/HDD/CDD (e.g. Sicily vs Lapland; Finland vs Italy).
12. Confirm `ln_*` columns have no `-Inf`/`NaN` (log of 0 or negative).

---

## Country coverage — is Norway / Switzerland / UK in every control?

Declared scope (proposal slide 6): EU27 + Norway + UK + Switzerland. Actual coverage differs
by dataset availability:

| control / output            | Norway | Switzerland | UK  | notes |
|-----------------------------|:------:|:-----------:|:---:|-------|
| Crop panel (9)              |   NO   |     NO      | NO  | CropStatHarm has no NUTS3 yield for them; crop = 15 EU countries only |
| Energy outcome (16)         |  yes   |   **NO**    | yes | CH absent from Eurostat nrg_bal_c |
| Energy weather HDD/CDD (13) |  yes   |     yes     | yes | 33 countries |
| Population (15)             |  yes   |     yes     | yes | |
| Electricity price (14)      |  yes   |   **NO**    | yes | |
| Gas price (14)              | **NO** |   **NO**    | yes | Norway reports little/no household gas price |

Consequences to verify downstream:
- **Switzerland**: appears only in intermediates (13, 15); it is NOT in the energy estimation
  panel (no outcome). Effectively excluded from energy estimation.
- **Norway**: in the electricity model, but DROPS from the gas model (no gas price, minimal gas use).
- **UK**: fully covered on the energy side.
- **All three absent from the crop side** — the crop estimation is only 15 EU countries:
  CZ DE DK EE EL ES FI FR HU IT LT LV RO SE SK.

If the thesis needs NO/CH/UK crop results, CropStatHarm cannot supply them — would require FAOSTAT
national yield for those countries (a different, country-level crop design).

---

## Known gaps (not bugs, just not built)
- Scenario replay (ISIMIP Scenario A, AMOC-stressed Scenario B) — slides 22-27.
- Valuation / country exposure ranking — slides 28-29.
- Numbering gaps vs proposal (no 10/11 etc.) — our numbering diverged; harmless.
