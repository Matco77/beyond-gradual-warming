# Appendix M — Panel Balance of the Estimation Samples

This appendix documents how far the samples on which the response functions are
estimated are from a balanced panel, where the missing observations are, and why
they are missing. It covers the crop panels of specification F (the main
specification) and specification A (the March–July comparison), the mixed-level
crop panel used in the spec-F comparison column, and the two household-energy
panels. For every panel it reports the first and last year of each country.

All figures were computed on 24 September 2026 directly from the data files in
`Amoc/datasets/`: `2.nrg_bal_c_prepared.csv`, `6.CropStatHarm_prepared.csv`,
`6.eobs_to_nuts_crosswalk.csv`, `6.NUTS2016_cropregions.gpkg`,
`8.eobs_nuts3_crop_weather_window.csv`, `9.crop_panel_nuts3_estimation.csv`,
`10.crop_panel_mixed_estimation.csv`, `11e.crop_precip_cellsq.csv`,
`11g.crop_weather_spec_f.csv`, `14.energy_prices_prepared.csv` and
`16.energy_panel_estimation.csv`. The estimation samples were rebuilt with
`regression.R`'s own definitions (`common_obs()`, `add_f()`, `f_F`, `f_tr` and the
formulas of `energy_tab()`), and their sizes reproduce the *Observations* line of
every table in `Amoc/results/` (`crop_<crop>_n3_vs_mixed.txt`,
`energy_electricity.txt`, `energy_gas.txt`).

---

## M.1 Definitions

**Samples.** The crop regressions are run crop by crop on a panel of NUTS3 regions ×
harvest years. Specifications F and A are estimated on **the same rows**: for each
crop `regression.R` builds one sample with
`common_obs(add_f(cp9[crop == k]), list(f_F, f_tr), ~cntr)`, which keeps only the
rows usable by both formulas after removing missing values and fixed-effect
singletons, and repeats the removal until nothing more drops. The OLS, FE and
"+ spring wheat N/B" columns of the crop tables are fitted on the same sample, so
one description covers all of them; the "Drop snapped EL62*" robustness column uses
this sample without the EL62 regions. The mixed-level panel (§M.3) is a separate
sample, used only for the "Mixed FE+Trends" column.

The energy regressions are run fuel by fuel on countries × years. The
price-controlled columns (Pooled OLS, Year FE, Country trends, Benchmark, No price
same rows) share the sample `common_obs(dt, list(f_yfe, f_ctr, f_bm), ~country_id)`,
called the **benchmark sample** below. The "No price, full sample" column uses
`common_obs(dt_all, list(f_np), ~country_id)`, called the **weather-only sample**.
The "Drop Malta (MLT)" robustness column is the benchmark sample without Malta. The
population-elasticity check (`energy_pop_elasticity.txt`), which does not go through
`common_obs()`, is not covered.

**Balance.** A panel is balanced when every unit is observed in every year of its
calendar, the years from the first to the last year present in that sample. Three
measures are used:

- **fill** = rows / (units × calendar years), which equals 1 for a balanced panel;
- **empty cells** = units × calendar years − rows, split into years *before* a
  unit's first observation (late entry), years *after* its last observation (early
  exit), and *internal gaps* (years missing between a unit's first and last
  observation);
- **years per unit**, the number of years in which each unit is observed.

---

## M.2 Crop panels of specifications F and A (NUTS3)

### M.2.1 Degree of imbalance

**Table M.1** — Size and balance of the crop estimation samples (spec F = spec A).

| Crop | Rows | NUTS3 regions | Countries | Harvest years | Cells (regions × years) | Fill | Regions observed every year | Years per region (min / median / max) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Soft wheat | 18,550 | 787 | 15 | 1989–2023 (35) | 27,545 | 0.67 | 93 | 2 / 24 / 35 |
| Durum wheat | 7,280 | 345 | 7 | 1989–2023 (35) | 12,075 | 0.60 | 37 | 2 / 24 / 35 |
| Spring barley | 13,397 | 603 | 12 | 1989–2023 (35) | 21,105 | 0.63 | 86 | 2 / 24 / 35 |
| Winter barley | 15,725 | 730 | 13 | 1989–2023 (35) | 25,550 | 0.62 | 91 | 2 / 23.5 / 35 |
| Grain maize | 9,974 | 401 | 11 | 1989–2023 (35) | 14,035 | 0.71 | 87 | 2 / 25 / 35 |
| Sugar beet | 10,349 | 581 | 13 | 1989–2022 (34) | 19,754 | 0.52 | 0 | 2 / 18 / 33 |
| Sunflower | 6,967 | 322 | 7 | 1989–2023 (35) | 11,270 | 0.62 | 74 | 2 / 24 / 35 |

None of the seven panels is balanced. The fill ranges from 0.52 (sugar beet) to
0.71 (grain maize), and at most 23% of the regions of a crop are observed in every
year (sunflower, 74 of 322). For sugar beet no region is: it has no observation in
2023, so its calendar ends in 2022, and no region covers all 34 years.

**Table M.2** — Where the empty cells are. The last two columns count regions whose
first year is later, or whose last year is earlier, than that of their country.

| Crop | Empty cells | Before region's first year | After region's last year | Internal gaps | Internal gaps, % of in-span years | Regions with ≥ 1 internal gap | Regions starting after their country | Regions ending before their country |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Soft wheat | 8,995 | 7,174 | 1,427 | 394 | 2.1% | 126 | 70 | 104 |
| Durum wheat | 4,795 | 3,538 | 532 | 725 | 9.1% | 110 | 109 | 40 |
| Spring barley | 7,708 | 5,331 | 1,715 | 662 | 4.7% | 165 | 43 | 172 |
| Winter barley | 9,825 | 7,403 | 1,898 | 524 | 3.2% | 213 | 45 | 129 |
| Grain maize | 4,061 | 3,386 | 569 | 106 | 1.1% | 49 | 19 | 30 |
| Sugar beet | 9,405 | 6,112 | 2,465 | 828 | 7.4% | 174 | 92 | 267 |
| Sunflower | 4,303 | 3,420 | 486 | 397 | 5.4% | 75 | 38 | 24 |

Most of the imbalance is staggered entry. Between 65% (sugar beet) and 83% (grain
maize) of the empty cells are years before a region's first observation, 11–26% are
years after its last observation, and 2.6–15.1% are internal gaps. Internal gaps are
rare relative to the years a region spends in the panel: 1.1% (grain maize) to 9.1%
(durum wheat) of its in-span years. Entry happens mostly by whole countries
(§M.2.2), but not only: inside a country, 19–109 regions per crop start later and
24–267 stop earlier than the country itself.

### M.2.2 Start and end by country

**Table M.3** — First and last harvest year in the estimation sample, by country and
crop (number of NUTS3 regions in brackets). "—" = no NUTS3 row for that crop.

| Country | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CZ Czechia | 1998–2022 (14) | — | 1998–2022 (14) | 1998–2022 (14) | 2005–2022 (14) | 1998–2022 (11) | 2005–2020 (14) |
| DE Germany | 1999–2022 (345) | — | 1999–2022 (350) | 1999–2022 (356) | — | 1999–2022 (302) | — |
| DK Denmark | 2006–2022 (11) | — | 2006–2022 (11) | 2006–2022 (11) | 2011–2022 (11) | 2006–2022 (10) | — |
| EE Estonia | 2004–2022 (5) | — | 2004–2022 (5) | 2004–2022 (5) | — | — | — |
| EL Greece | 1998–2021 (46) | 1998–2021 (46) | — | 1998–2020 (46) | 2009–2021 (46) | 2009–2021 (22) | 2009–2021 (34) |
| ES Spain | 1998–2021 (51) | 1998–2021 (42) | 1998–2021 (48) | 1998–2021 (46) | 1998–2021 (51) | 1998–2021 (30) | 1998–2021 (40) |
| FI Finland | 1998–2022 (18) | — | 1999–2022 (19) | — | — | 1998–2022 (13) | — |
| FR France | 1989–2023 (95) | 1989–2023 (91) | 1989–2023 (92) | 1989–2023 (94) | 1989–2023 (94) | 1989–2020 (54) | 1989–2023 (91) |
| HU Hungary | 2002–2022 (20) | 2000–2022 (20) | 1996–2022 (20) | 1996–2022 (20) | 1996–2022 (20) | 1996–2022 (20) | 1996–2022 (20) |
| IT Italy | 1995–2022 (100) | 1995–2022 (100) | — | 2006–2020 (108) | 1995–2022 (104) | 2006–2022 (54) | 2006–2022 (81) |
| LT Lithuania | 2000–2022 (10) | — | 2000–2022 (10) | 2000–2022 (10) | 2002–2022 (10) | 1998–2022 (10) | — |
| LV Latvia | 1997–2022 (5) | — | 1997–2022 (5) | 1997–2022 (5) | — | — | — |
| RO Romania | 1998–2022 (42) | 1998–2022 (38) | — | — | 1998–2022 (42) | 1990–2022 (42) | 1998–2022 (42) |
| SE Sweden | 1990–2020 (17) | — | 1990–2022 (21) | 1995–2022 (7) | 2007–2022 (1) | 1990–2022 (5) | — |
| SK Slovakia | 2017–2018 (8) | 2017–2022 (8) | 2007–2022 (8) | 2007–2022 (8) | 2007–2022 (8) | 1997–2022 (8) | — |

- For every crop × country pair, the first and last year in Table M.3 are exactly the
  first and last year in which CropStatHarm reports a positive NUTS3 yield for a
  region that has weather. The start and end dates are those of the source: the
  pipeline does not shorten any country's series.
- France is the only country observed from 1989, and the only one observed in 2023
  (for every crop except sugar beet, which has no 2023 observation and whose French
  series ends in 2020). The other 71 country series start between 1990 and 2017, 21
  of them in 1998, and end between 2018 and 2022.
- Germany reports soft wheat, spring barley, winter barley and sugar beet, always
  from 1999 to 2022 and with 302–356 regions. It supplies 39.9% (soft wheat), 49.1%
  (spring barley), 47.2% (winter barley) and 51.9% (sugar beet) of the rows. It is
  absent from the durum wheat, grain maize and sunflower samples; its maize and
  sunflower yields are reported at NUTS1 (§M.3).
- Italy enters in 1995 for soft wheat, durum wheat and grain maize, but only in 2006
  for winter barley, sugar beet and sunflower. Greece enters in 1998 for the wheats
  and winter barley and in 2009 for grain maize, sugar beet and sunflower. Spain
  covers 1998–2021 for all seven crops.
- Two country series are thin: Slovak soft wheat (8 regions observed only in 2017
  and 2018) and Swedish grain maize (a single region, 2007–2022).

### M.2.3 Balance within countries

**Table M.4** — Balance inside each country's own span. Each cell gives the fill,
rows / (regions × years from the country's first to its last year), followed by the
number of regions observed in every year of that span over the number of regions.
† = the whole country is missing in some year between its first and last year.

| Country | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CZ Czechia | 1.00 · 14/14 | — | 1.00 · 14/14 | 1.00 · 14/14 | 0.96 · 13/14 | 1.00 · 11/11 | 0.62 · 0/14 † |
| DE Germany | 0.89 · 203/345 | — | 0.78 · 117/350 | 0.87 · 203/356 | — | 0.74 · 104/302 | — |
| DK Denmark | 1.00 · 11/11 | — | 1.00 · 11/11 | 1.00 · 11/11 | 1.00 · 11/11 | 0.71 · 3/10 | — |
| EE Estonia | 1.00 · 5/5 | — | 0.96 · 3/5 | 0.71 · 0/5 | — | — | — |
| EL Greece | 0.92 · 37/46 | 0.90 · 35/46 | — | 0.94 · 39/46 | 0.95 · 40/46 | 0.80 · 6/22 | 0.74 · 16/34 |
| ES Spain | 0.91 · 39/51 | 0.75 · 18/42 | 0.85 · 33/48 | 0.82 · 29/46 | 0.95 · 31/51 | 0.66 · 10/30 | 0.94 · 31/40 |
| FI Finland | 0.77 · 10/18 | — | 0.80 · 14/19 | — | — | 0.54 · 3/13 | — |
| FR France | 0.99 · 93/95 | 0.66 · 37/91 | 0.98 · 86/92 | 0.99 · 91/94 | 0.98 · 87/94 | 0.72 · 27/54 | 0.91 · 74/91 |
| HU Hungary | 1.00 · 20/20 | 0.91 · 8/20 | 0.98 · 18/20 | 0.99 · 18/20 | 0.99 · 18/20 | 0.92 · 12/20 | 0.99 · 18/20 |
| IT Italy | 0.89 · 69/100 | 0.90 · 68/100 | — | 0.87 · 0/108 † | 0.94 · 77/104 | 0.61 · 12/54 | 0.87 · 51/81 |
| LT Lithuania | 1.00 · 9/10 | — | 1.00 · 10/10 | 1.00 · 10/10 | 0.99 · 8/10 | 0.69 · 4/10 | — |
| LV Latvia | 1.00 · 5/5 | — | 1.00 · 5/5 | 1.00 · 5/5 | — | — | — |
| RO Romania | 0.98 · 38/42 | 0.48 · 1/38 | — | — | 0.98 · 38/42 | 0.67 · 14/42 | 0.92 · 34/42 |
| SE Sweden | 0.83 · 10/17 | — | 0.97 · 15/21 | 0.60 · 1/7 | 0.75 · 0/1 † | 0.83 · 2/5 | — |
| SK Slovakia | 1.00 · 8/8 | 0.83 · 4/8 | 0.98 · 6/8 | 1.00 · 8/8 | 0.98 · 6/8 | 0.77 · 4/8 | — |

† Years in which the whole country is missing: Czechia, sunflower: 2008–2013; Italy, winter barley: 2015; Sweden, grain maize: 2008–2010, 2015.

Balance inside a country varies widely. Some national series are complete: every
Czech, Danish and Latvian region is observed in every year of its country's span for
soft wheat, spring barley and winter barley. Others are thin: German spring barley
(117 of 350 regions complete, fill 0.78) and sugar beet (104 of 302, fill 0.74),
French durum wheat (37 of 91, fill 0.66) and Romanian durum wheat (1 of 38, fill
0.48). No Italian winter-barley region and no Czech sunflower region can be complete,
because the country itself is missing in 2015 and in 2008–2013 respectively. Sugar
beet is the least balanced crop inside countries, with a fill from 0.54 (Finland) to
1.00 (Czechia).

### M.2.4 Composition over time

**Table M.5** — Rows per harvest year (number of countries in brackets).

| Year | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1989 | 94 (1) | 49 (1) | 91 (1) | 94 (1) | 93 (1) | 41 (1) | 81 (1) |
| 1990 | 105 (2) | 50 (1) | 109 (2) | 94 (1) | 93 (1) | 83 (3) | 82 (1) |
| 1991 | 105 (2) | 52 (1) | 109 (2) | 94 (1) | 93 (1) | 87 (3) | 80 (1) |
| 1992 | 105 (2) | 51 (1) | 109 (2) | 94 (1) | 93 (1) | 85 (3) | 81 (1) |
| 1993 | 104 (2) | 48 (1) | 104 (2) | 94 (1) | 93 (1) | 83 (3) | 80 (1) |
| 1994 | 104 (2) | 47 (1) | 103 (2) | 93 (1) | 93 (1) | 82 (3) | 79 (1) |
| 1995 | 184 (3) | 128 (2) | 106 (2) | 96 (2) | 186 (2) | 82 (3) | 81 (1) |
| 1996 | 185 (3) | 132 (2) | 124 (3) | 112 (3) | 207 (3) | 100 (4) | 101 (2) |
| 1997 | 195 (4) | 133 (2) | 131 (4) | 119 (4) | 211 (3) | 108 (5) | 102 (2) |
| 1998 | 350 (9) | 217 (5) | 185 (6) | 216 (7) | 300 (5) | 162 (9) | 177 (4) |
| 1999 | 671 (10) | 218 (5) | 525 (8) | 549 (8) | 302 (5) | 424 (10) | 177 (4) |
| 2000 | 687 (11) | 233 (6) | 531 (9) | 560 (9) | 306 (5) | 418 (10) | 183 (4) |
| 2001 | 684 (11) | 235 (6) | 535 (9) | 557 (9) | 305 (5) | 408 (10) | 180 (4) |
| 2002 | 704 (12) | 227 (6) | 527 (9) | 556 (9) | 305 (6) | 407 (10) | 179 (4) |
| 2003 | 705 (12) | 229 (6) | 535 (9) | 561 (9) | 312 (6) | 413 (10) | 175 (4) |
| 2004 | 711 (13) | 225 (6) | 531 (10) | 557 (10) | 313 (6) | 406 (10) | 176 (4) |
| 2005 | 704 (13) | 223 (6) | 525 (10) | 554 (10) | 325 (7) | 394 (10) | 188 (5) |
| 2006 | 730 (14) | 234 (6) | 552 (11) | 675 (12) | 324 (7) | 452 (12) | 262 (6) |
| 2007 | 732 (14) | 238 (6) | 556 (12) | 695 (13) | 333 (9) | 437 (12) | 262 (6) |
| 2008 | 739 (14) | 255 (6) | 549 (12) | 693 (13) | 333 (8) | 407 (12) | 254 (5) |
| 2009 | 731 (14) | 262 (6) | 538 (12) | 689 (13) | 378 (9) | 423 (13) | 273 (6) |
| 2010 | 738 (14) | 275 (6) | 546 (12) | 675 (13) | 376 (9) | 423 (13) | 275 (6) |
| 2011 | 732 (14) | 278 (6) | 537 (12) | 675 (13) | 391 (11) | 431 (13) | 277 (6) |
| 2012 | 732 (14) | 281 (6) | 545 (12) | 672 (13) | 389 (11) | 448 (13) | 258 (6) |
| 2013 | 732 (14) | 285 (6) | 521 (12) | 657 (13) | 392 (11) | 444 (13) | 272 (6) |
| 2014 | 739 (14) | 278 (6) | 528 (12) | 674 (13) | 390 (11) | 442 (13) | 290 (7) |
| 2015 | 713 (14) | 279 (6) | 522 (12) | 561 (12) | 388 (10) | 424 (13) | 291 (7) |
| 2016 | 705 (14) | 284 (6) | 479 (12) | 644 (13) | 388 (11) | 368 (13) | 288 (7) |
| 2017 | 717 (15) | 303 (7) | 427 (12) | 607 (13) | 378 (11) | 327 (13) | 287 (7) |
| 2018 | 710 (15) | 307 (7) | 452 (12) | 605 (13) | 382 (11) | 334 (13) | 295 (7) |
| 2019 | 702 (14) | 304 (7) | 434 (12) | 612 (13) | 379 (11) | 341 (13) | 294 (7) |
| 2020 | 667 (14) | 302 (7) | 434 (12) | 608 (13) | 371 (11) | 331 (13) | 294 (7) |
| 2021 | 664 (13) | 306 (7) | 412 (12) | 463 (11) | 371 (11) | 281 (12) | 281 (6) |
| 2022 | 576 (11) | 227 (5) | 393 (11) | 426 (10) | 289 (9) | 253 (10) | 222 (4) |
| 2023 | 94 (1) | 85 (1) | 92 (1) | 94 (1) | 92 (1) | — | 90 (1) |

- Up to 1997 each crop sample holds 1–5 countries and 41–211 rows a year. The large
  entries come in 1998 (Spain, and depending on the crop Czechia, Greece, Romania,
  Finland and Lithuania) and in 1999 (Germany). From 1999 to 2021 the number of
  countries moves within a narrow band (soft wheat 10–15, durum wheat 5–7, spring
  barley 8–12, winter barley 8–13, grain maize 5–11, sugar beet 10–13, sunflower
  4–7).
- The number of rows still drifts inside that band. Spring barley falls from 556 rows
  in 2007 to 412 in 2021, entirely through Germany (320 → 169 German rows; the other
  countries go from 236 to 243). Sugar beet falls from 452 rows in 2006 to 281 in
  2021 (Germany 257 → 150, the other countries 195 → 131).
- The samples shrink in 2022, because Spain and Greece end in 2021, and in 2023 only
  France is left.
- Only French regions are present in these years: soft wheat and spring barley 1989
  and 2023; durum wheat, winter barley and grain maize 1989–1994 and 2023; sunflower
  1989–1995 and 2023; sugar beet 1989. They hold 188 (soft wheat, 1.0% of its sample), 382 (durum wheat, 5.2% of its sample), 183 (spring barley, 1.4% of its sample), 657 (winter barley, 4.2% of its sample), 650 (grain maize, 6.5% of its sample), 41 (sugar beet, 0.4% of its sample) and 654 (sunflower, 9.4% of its sample) rows. In those years the
  year effect absorbs the French mean, so the weather coefficients are identified
  only from differences between French regions.
- Between 87.4% (grain maize) and 93.7% (winter barley) of the rows fall in
  1998–2022, the years in which most countries are present.

### M.2.5 Years per region

**Table M.6** — Regions by number of years observed (rows in brackets).

| Years observed | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2–3 | 16 (36) | 10 (23) | 9 (23) | 8 (23) | 3 (7) | 40 (97) | 6 (13) |
| 4–5 | 6 (27) | 18 (83) | 6 (25) | 9 (40) | 2 (10) | 23 (102) | 10 (46) |
| 6–10 | 24 (196) | 40 (314) | 32 (253) | 25 (205) | 7 (63) | 63 (504) | 36 (314) |
| 11–20 | 106 (1745) | 60 (942) | 178 (2957) | 227 (3428) | 82 (1176) | 211 (3370) | 96 (1504) |
| 21–30 | 532 (12981) | 174 (4426) | 273 (6510) | 369 (8813) | 217 (5573) | 197 (4758) | 98 (2433) |
| 31–35 | 103 (3565) | 43 (1492) | 105 (3629) | 92 (3216) | 90 (3145) | 47 (1518) | 76 (2657) |

The median region is observed for 18 (sugar beet) to 25 (grain maize) years.
Regions observed for 21 years or more supply 61% (sugar beet) to 89% (soft wheat) of
the rows. At the other end, 3–40 regions per crop are observed for only two or three
years (§M.2.7).

### M.2.6 Why observations are missing

**Table M.7** — From the CropStatHarm NUTS3 records to the estimation sample (rows).

| Crop | CropStatHarm NUTS3 records | No yield | Yield = 0 | Excluded by hand | Region without weather | Panel 9 | Missing spec-F regressor | Singleton regions | Estimation sample |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Soft wheat | 47,388 | 28,615 | 152 | 0 | 67 | 18,554 | 0 | 4 | 18,550 |
| Durum wheat | 42,062 | 32,249 | 2,468 | 1 | 48 | 7,296 | 2 | 14 | 7,280 |
| Spring barley | 40,359 | 24,099 | 2,855 | 0 | 7 | 13,398 | 0 | 1 | 13,397 |
| Winter barley | 43,517 | 27,223 | 504 | 0 | 61 | 15,729 | 1 | 3 | 15,725 |
| Grain maize | 10,989 | 19 | 943 | 0 | 53 | 9,974 | 0 | 0 | 9,974 |
| Sugar beet | 13,537 | 369 | 2,796 | 0 | 0 | 10,372 | 0 | 23 | 10,349 |
| Sunflower | 10,102 | 18 | 3,115 | 0 | 0 | 6,969 | 0 | 2 | 6,967 |

- **No yield.** For the wheats and barleys almost all of these are area-only records
  from the area gap-fill source (`area_source = cropareaRegional`): 28,506 of 28,615 (soft wheat), 32,218 of 32,249 (durum wheat), 23,743 of 24,099 (spring barley), 26,888 of 27,223 (winter barley); grain maize, sugar beet and sunflower have no area-only records. No
  yield is ever imputed for them.
- **Yield = 0.** In almost every zero-yield record the area is also zero or missing
  (all but 2 soft-wheat, 2 durum-wheat, 50 spring-barley and 1 sugar-beet records), and production is zero or missing in all of them. These records
  describe a crop that was not grown, not a failed harvest, and the log outcome could
  not use them anyway. No yield in the source is negative.
- **Excluded by hand.** ES300 durum wheat 2012, a source area error documented in
  `crop_panel_nuts3.R`.
- **Region without weather.** Five NUTS3 regions have no cell in the E-OBS crosswalk
  (`6.eobs_to_nuts_crosswalk.csv`): EL413 Chios; EL422 Andros, Thira, Kea, Milos, Mykonos, Naxos, Paros, Syros, Tinos; ES703 El Hierro; ES706 La Gomera; ES708 Lanzarote (names from the NUTS 2016 geometry
  file).
- **Missing spec-F regressor.** Overwinter precipitation (`precip_w`) is missing in
  `11g` for ITG18 and ITG19 in harvest year 2016. These rows (durum wheat ITG18 and
  ITG19, winter barley ITG19) have every spec-A regressor, but they leave the spec-A
  regression as well, because the two specifications share one sample.
- **Singleton regions.** These are regions left with a single observation, which
  their own fixed effect fits exactly. Every such region drops together with its one
  row: soft wheat 4 (DE402, DEA16, ITF33, ITF45); durum wheat 14 (EL622, ES113, ES130, ES211, ES704, ES705, FRD12, ITC33, ITC34, ITC41, ITH44, RO215, RO321, RO412); spring barley 1 (DE27A); winter barley 3 (DE942, DEA16, SE121); grain maize 0; sugar beet 23 (DE12A, DE13A, DE142, DE272, DE94H, DED42, EL301, EL411, EL412, EL543, EL641, ES424, FRJ23, ITC12, ITF33, ITF44, ITF45, ITF61, ITF62, ITF63, ITH44, ITI1A, SE232); sunflower 2 (EL432, ITH44).

**Table M.8** — Internal gaps (region-years missing between a region's first and last
observation) by reason.

| Reason | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| No CropStatHarm record | 2 | 3 | 186 | 1 | 53 | 651 | 212 |
| Area-only record (cropareaRegional), no yield | 336 | 452 | 354 | 339 | 0 | 0 | 0 |
| CropStatHarm record without yield | 18 | 22 | 22 | 123 | 10 | 77 | 9 |
| Yield recorded as 0 | 38 | 245 | 100 | 60 | 43 | 100 | 176 |
| Excluded by hand (ES300, 2012) | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| Spec-F regressor missing (ITG18/ITG19, 2016) | 0 | 2 | 0 | 1 | 0 | 0 | 0 |
| **Total** | 394 | 725 | 662 | 524 | 106 | 828 | 397 |

Of the 3,636 internal gaps, only 4 are created by the pipeline: the ES300
exclusion and the three rows without `precip_w`. Every other internal gap is a
region-year for which the source holds no positive yield. Among the internal gaps
with a zero yield, the area is zero or missing in every case except 26 spring barley region-years.

### M.2.7 Consequences for estimation

**Selection.** The fixed-effects estimator stays consistent on an unbalanced panel
as long as whether a region-year is observed is unrelated to the idiosyncratic yield
shock, given the regressors, the fixed effects and the region trends. Entry and exit
that depend on the region or the country (a country joining the database in 1999),
or on the weather regressors themselves, do not bias it. Two checks were run on the
benchmark sample (Table M.9):

1. **Variable-addition test** (Nijman and Verbeek 1992; Wooldridge 2010, ch. 19).
   Two indicators are added to the spec-F benchmark: s(t−1) and s(t+1) equal 1 when
   the region is observed in the previous or the next harvest year of the sample's
   calendar. Standard errors are clustered by country, as in the benchmark. Under
   ignorable selection both coefficients are zero. Rows in the first and last
   calendar year drop, because one of the indicators is undefined there.
2. **Does weather predict whether a year is observed?** This is a linear probability
   model of the observation indicator on the spec-F weather regressors, with region
   and year fixed effects, over each region's own first-to-last span. Here standard
   errors are clustered by **region**. A country-clustered covariance matrix has rank
   at most G − 1; for durum wheat that is 6 (G = 7 countries), below its 9 weather
   regressors, so the joint test would not be defined. Region clustering is used for
   every crop so that the rows are comparable. It ignores the correlation between
   regions of the same country, so these p-values are, if anything, too small.

**Table M.9** — Selection diagnostics on the spec-F benchmark sample. Coefficients of
s(t−1) and s(t+1) are in log points, with p-values in brackets. The last column is
the change in the probability of observation, in percentage points, for a one
within-region standard deviation of the most influential weather variable.

| Crop | Rows | Rows with s(t−1) = 0 / s(t+1) = 0 | s(t−1) | s(t+1) | Joint p | In-span years missing | LPM joint F (p) | Largest 1-SD effect, pp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Soft wheat | 18,362 | 867 / 867 | −0.017 (0.365) | −0.043 (0.095) | 0.153 | 2.1% | 2.14 (0.024) | +0.57 (precip_s) |
| Durum wheat | 7,146 | 524 / 490 | +0.035 (0.334) | +0.002 (0.926) | 0.122 | 9.1% | 3.45 (<0.001) | −1.65 (heat_s) |
| Spring barley | 13,214 | 798 / 796 | −0.047 (0.134) | −0.020 (0.584) | 0.209 | 4.7% | 2.70 (0.020) | +0.74 (precip_s2) |
| Winter barley | 15,537 | 921 / 923 | −0.014 (0.563) | +0.002 (0.900) | 0.840 | 3.2% | 8.06 (<0.001) | −1.26 (gdd_s) |
| Grain maize | 9,789 | 369 / 370 | +0.008 (0.859) | +0.061 (0.275) | 0.525 | 1.1% | 1.61 (0.156) | −0.58 (precip_s) |
| Sugar beet | 10,055 | 794 / 592 | +0.081 (0.071) | +0.010 (0.607) | 0.139 | 7.4% | 3.05 (0.010) | −1.24 (gdd_s) |
| Sunflower | 6,794 | 352 / 344 | +0.020 (0.643) | +0.065 (0.025) | 0.062 | 5.4% | 2.70 (0.021) | −2.79 (precip_s) |

The variable-addition test does not reject jointly for any crop at the 5% level
(joint p from 0.062 to 0.840). One of the fourteen single coefficients has
p < 0.05 (sunflower s(t+1), p = 0.025), which is about what chance alone produces
over fourteen tests; with 7–15 clusters the clustered p-values are, if anything, too
small. Yields are not systematically different in the years next to a missing year.

Weather does predict, weakly, whether a year is observed. The joint F is
significant at 5% for six of the seven crops (not grain maize, p = 0.156), but one
standard deviation of any weather variable moves the probability of observation by
at most 2.8 percentage points (sunflower, season precipitation). This is selection
on the regressors, and it does not bias the estimator. It would matter only if the
yield response in the missing years differed from the response in the observed ones,
and the data cannot test that.

**Regions with two or three years.** The benchmark gives each region a fixed effect
and its own linear and quadratic trend, which makes three parameters. A region
observed two or three times is therefore fitted exactly: across these regions, the
largest absolute residual is between 6 × 10⁻¹⁴ and 2 × 10⁻¹² for every crop. Such a
region adds nothing to the weather coefficients but still counts in N. There are
3–40 such regions per crop, with 7–97 rows, or 0.07–0.94% of the sample (first row
of Table M.6). Dropping them leaves every spec-F coefficient unchanged: the largest
absolute difference is 1.2e-13. Regions with four or five years leave one or two
degrees of freedom after their trend (10–102 rows per crop).

**A country that carries no information.** In soft wheat all eight Slovak regions are
observed only in 2017 and 2018, so the whole Slovak cluster is fitted exactly. The
soft-wheat regression reports 15 country clusters, and the wild-cluster bootstrap
counts G = 15, but only 14 carry identifying variation. Dropping Slovakia leaves the
coefficients unchanged (largest absolute difference 6.4e-18). It changes the
clustered standard errors only through the small-sample factor, by at most
0.21%. No other crop has such a country.

---

## M.3 Mixed-level crop panel (spec-F comparison column)

The mixed panel adds the countries that report yields only at NUTS2 or NUTS0, and it
runs from 1975 to 2024. A unit is a region at the level at which its country reports
(NUTS3, NUTS2 or NUTS0). The panel is less balanced than the NUTS3 one.

**Table M.10** — Size and balance of the mixed-level samples.

| Crop | Rows | Units | Countries | Harvest years | Fill | Units observed every year | Years per unit (min / median / max) | Before unit's first year | After unit's last year | Internal gaps |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Soft wheat | 20,338 | 856 | 26 | 1975–2024 (50) | 0.48 | 1 | 2 / 24 / 50 | 19,589 | 2,372 | 501 |
| Durum wheat | 7,763 | 369 | 13 | 1986–2024 (39) | 0.54 | 0 | 2 / 24 / 38 | 4,895 | 949 | 784 |
| Spring barley | 14,729 | 654 | 18 | 1975–2023 (49) | 0.46 | 1 | 2 / 24 / 49 | 14,883 | 1,772 | 662 |
| Winter barley | 17,229 | 791 | 22 | 1975–2024 (50) | 0.44 | 1 | 2 / 23 / 50 | 19,030 | 2,766 | 525 |
| Grain maize | 11,505 | 466 | 20 | 1975–2023 (49) | 0.50 | 0 | 2 / 25 / 48 | 10,569 | 650 | 110 |
| Sugar beet | 11,943 | 636 | 21 | 1975–2023 (49) | 0.38 | 0 | 2 / 18 / 48 | 15,156 | 3,176 | 889 |
| Sunflower | 8,039 | 368 | 15 | 1975–2023 (49) | 0.45 | 0 | 2 / 23.5 / 48 | 8,880 | 540 | 573 |

In this panel the NUTS3 countries have exactly the first and last years of
Table M.3 (checked pair by pair). The countries it adds are listed below.

**Table M.11** — First and last harvest year of the countries added by the mixed
panel (number of units in brackets; reporting level after the name).

| Country | Soft wheat | Durum wheat | Spring barley | Winter barley | Grain maize | Sugar beet | Sunflower |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AT Austria (NUTS2) | 1995–2022 (9) | 1995–2022 (9) | 1975–2022 (9) | 1975–2022 (9) | 1975–2022 (9) | 1975–2022 (7) | 1975–2022 (9) |
| BE Belgium (NUTS2) | 1975–2022 (11) | — | 2009–2022 (10) | 2009–2022 (11) | 2011–2022 (11) | 1975–2022 (11) | — |
| BG Bulgaria (NUTS2) | 1995–2020 (6) | 1998–2022 (6) | — | — | 1991–2021 (6) | — | 1991–2021 (6) |
| CY Cyprus (NUTS0) | 2012–2024 (1) | 1987–2024 (1) | — | 2004–2024 (1) | — | — | — |
| HR Croatia (NUTS2) | 2008–2022 (2) | 2008–2022 (2) | 2008–2022 (2) | 2008–2022 (2) | 2005–2022 (2) | 2005–2022 (1) | 2005–2022 (2) |
| IE Ireland (NUTS2) | 1990–2023 (3) | — | — | 2010–2023 (3) | — | — | — |
| LU Luxembourg (NUTS0) | 1975–2024 (1) | 2020–2024 (1) | 1975–2023 (1) | 1975–2024 (1) | 2000–2023 (1) | 1975–2023 (1) | 2020–2023 (1) |
| NL Netherlands (NUTS2) | 1994–2022 (12) | — | 1994–2022 (12) | 1994–2022 (12) | 2008–2022 (12) | 1994–2022 (12) | — |
| PL Poland (NUTS2) | 2003–2022 (17) | — | 2003–2022 (17) | 2003–2022 (17) | 2003–2022 (17) | 1999–2022 (17) | 1999–2023 (17) |
| PT Portugal (NUTS2) | 1986–2022 (5) | 1986–2022 (5) | — | 2000–2020 (5) | 1986–2022 (5) | 1994–2008 (4) | 1986–2022 (5) |
| SI Slovenia (NUTS2) | 2007–2021 (2) | — | — | — | 2007–2023 (2) | 2019–2023 (2) | 2014–2023 (2) |
| SK Slovakia (NUTS2) | — | — | — | — | — | — | 1995–2023 (4) |

- Slovakia reports sunflower at NUTS2 (4 units, 1995–2023) and every other crop at
  NUTS3, so Slovak sunflower appears only in the mixed panel.
- Germany reports grain maize (13 NUTS1 units, 2010–2023, 167 rows) and sunflower
  (12 NUTS1 units, 1994–2023, 229 rows) at NUTS1. The mixed panel carries weather
  only at NUTS3, NUTS2 and NUTS0, so these rows are in neither crop panel.
- Before 1989 only a few countries are observed. In 1975–1985 Austria is the only
  country in the grain maize and sunflower samples, and in 1986 Portugal is the only
  one in durum wheat.
- Luxembourg (one NUTS0 unit) is the only unit observed in every year, for soft
  wheat, spring barley and winter barley. No unit is observed in every year for the
  other crops.

---

## M.4 Energy panels (country × year)

### M.4.1 Degree of imbalance

**Table M.12** — Size and balance of the energy samples.

| Fuel | Sample | Rows | Countries | Years | Cells | Fill | Countries observed every year | Years per country (min / median / max) | Before first year | After last year | Internal gaps |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Electricity | Panel 16 | 1,010 | 29 | 1990–2024 (35) | 1,015 | 0.995 | 28 | 30 / 35 / 35 | 0 | 5 | 0 |
| Electricity | Benchmark sample | 844 | 29 | 1990–2024 (35) | 1,015 | 0.832 | 3 | 20 / 30 / 35 | 163 | 5 | 3 |
| Electricity | Weather-only sample | 1,008 | 29 | 1990–2024 (35) | 1,015 | 0.993 | 27 | 30 / 35 / 35 | 0 | 5 | 2 |
| Natural gas | Panel 16 | 904 | 27 | 1990–2024 (35) | 945 | 0.957 | 22 | 25 / 35 / 35 | 36 | 5 | 0 |
| Natural gas | Benchmark sample | 577 | 25 | 1991–2024 (34) | 850 | 0.679 | 1 | 13 / 21 / 34 | 267 | 6 | 0 |
| Natural gas | Weather-only sample | 904 | 27 | 1990–2024 (35) | 945 | 0.957 | 22 | 25 / 35 / 35 | 36 | 5 | 0 |

The panel on disk is close to balanced. For electricity only the United Kingdom's
2020–2024 years are missing. For gas, the 41 empty cells are the late start of four
countries' outcome series (36 cells) and the United Kingdom's exit (5 cells). The
benchmark samples are clearly unbalanced, with a fill of 0.832 (electricity) and
0.679 (gas). Almost all of their empty cells are years before a country's first
observation: 163 of 171 for electricity and 267 of 273 for gas. The weather-only
samples are close to balanced again (0.993 and 0.957).

### M.4.2 Start and end by country

**Table M.13** — Electricity, by country: years in panel 16; years with a household
price in `14.` (the price file runs from 1985 to 2025, beyond the outcome); benchmark
sample; years missing inside the benchmark span; weather-only sample. Number of years
in brackets.

| Country | Panel 16 | Price (file 14) | Benchmark sample | Missing inside benchmark span | Weather-only sample |
| --- | --- | --- | --- | --- | --- |
| AUT Austria | 1990–2024 (35) | 1996–2025 (30) | 1996–2024 (29) | none | 1990–2024 (35) |
| BEL Belgium | 1990–2024 (35) | 1985–2025 (41) | 1990–2024 (35) | none | 1990–2024 (35) |
| BGR Bulgaria | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| HRV Croatia | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| CYP Cyprus | 1990–2024 (35) | 1999–2025 (27) | 1999–2024 (26) | none | 1990–2024 (35) |
| CZE Czechia | 1990–2024 (35) | 2000–2025 (26) | 2000–2024 (25) | none | 1990–2024 (35) |
| DNK Denmark | 1990–2024 (35) | 1985–2025 (41) | 1990–2024 (35) | none | 1990–2024 (35) |
| EST Estonia | 1990–2024 (35) | 2002–2025 (24) | 2002–2024 (23) | none | 1990–2024 (35) |
| FIN Finland | 1990–2024 (35) | 1995–2025 (31) | 1995–2024 (30) | none | 1990–2024 (35) |
| FRA France | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| DEU Germany | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| GRC Greece | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| HUN Hungary | 1990–2024 (35) | 1992–2025 (34) | 1992–2024 (33) | none | 1990–2024 (35) |
| IRL Ireland | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| ITA Italy | 1990–2024 (35) | 1991–2025 (34) | 1991–2024 (33) | 2007 | 1990–2024 (35) |
| LVA Latvia | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| LTU Lithuania | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| LUX Luxembourg | 1990–2024 (35) | 1985–2025 (41) | 1990–2024 (35) | none | 1990–2024 (35) |
| MLT Malta | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (32) | 2014, 2016 | 1990–2024 (33) |
| NLD Netherlands | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| NOR Norway | 1990–2024 (35) | 1995–2025 (31) | 1995–2024 (30) | none | 1990–2024 (35) |
| POL Poland | 1990–2024 (35) | 2000–2025 (26) | 2000–2024 (25) | none | 1990–2024 (35) |
| PRT Portugal | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| ROU Romania | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| SVK Slovakia | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| SVN Slovenia | 1990–2024 (35) | 1992–2025 (34) | 1992–2024 (33) | none | 1990–2024 (35) |
| ESP Spain | 1990–2024 (35) | 1991–2025 (35) | 1991–2024 (34) | none | 1990–2024 (35) |
| SWE Sweden | 1990–2024 (35) | 1996–2025 (30) | 1996–2024 (29) | none | 1990–2024 (35) |
| CHE Switzerland | — | — | — | — | — |
| GBR United Kingdom | 1990–2019 (30) | 1991–2020 (30) | 1991–2019 (29) | none | 1990–2019 (30) |

**Table M.14** — Natural gas, same layout.

| Country | Panel 16 | Price (file 14) | Benchmark sample | Missing inside benchmark span | Weather-only sample |
| --- | --- | --- | --- | --- | --- |
| AUT Austria | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| BEL Belgium | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| BGR Bulgaria | 2000–2024 (25) | 2004–2025 (22) | 2004–2024 (21) | none | 2000–2024 (25) |
| HRV Croatia | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| CYP Cyprus | — | — | — | — | — |
| CZE Czechia | 1990–2024 (35) | 2000–2025 (26) | 2000–2024 (25) | none | 1990–2024 (35) |
| DNK Denmark | 1990–2024 (35) | 1999–2025 (27) | 1999–2024 (26) | none | 1990–2024 (35) |
| EST Estonia | 1990–2024 (35) | 2003–2025 (23) | 2003–2024 (22) | none | 1990–2024 (35) |
| FIN Finland | 1990–2024 (35) | — | — | — | 1990–2024 (35) |
| FRA France | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| DEU Germany | 1990–2024 (35) | 1985–2025 (41) | 1991–2024 (34) | none | 1990–2024 (35) |
| GRC Greece | 1999–2024 (26) | 2012–2025 (14) | 2012–2024 (13) | none | 1999–2024 (26) |
| HUN Hungary | 1990–2024 (35) | 1992–2025 (34) | 1992–2024 (33) | none | 1990–2024 (35) |
| IRL Ireland | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| ITA Italy | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| LVA Latvia | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| LTU Lithuania | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| LUX Luxembourg | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| MLT Malta | — | — | — | — | — |
| NLD Netherlands | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| NOR Norway | 2000–2024 (25) | — | — | — | 2000–2024 (25) |
| POL Poland | 1990–2024 (35) | 2000–2023 (24) | 2000–2023 (24) | none | 1990–2024 (35) |
| PRT Portugal | 1997–2024 (28) | 2002–2025 (24) | 2002–2024 (23) | none | 1997–2024 (28) |
| ROU Romania | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| SVK Slovakia | 1990–2024 (35) | 2004–2025 (22) | 2004–2024 (21) | none | 1990–2024 (35) |
| SVN Slovenia | 1990–2024 (35) | 1995–2025 (31) | 1995–2024 (30) | none | 1990–2024 (35) |
| ESP Spain | 1990–2024 (35) | 2005–2025 (21) | 2005–2024 (20) | none | 1990–2024 (35) |
| SWE Sweden | 1990–2024 (35) | 1996–2025 (30) | 1996–2024 (29) | none | 1990–2024 (35) |
| CHE Switzerland | — | — | — | — | — |
| GBR United Kingdom | 1990–2019 (30) | 1991–2020 (30) | 1991–2019 (29) | none | 1990–2019 (30) |

For every country and both fuels, the panel-16 years are exactly the years in which
the Eurostat energy balance (`2.nrg_bal_c_prepared.csv`) reports a positive household
quantity.

### M.4.3 Why observations are missing

**Table M.15** — From panel 16 to the estimation samples (rows).

| Fuel | Panel 16 | No price | Missing weather | Singleton | Benchmark sample | Weather-only sample |
| --- | --- | --- | --- | --- | --- | --- |
| Electricity | 1,010 | 164 | 2 | 0 | 844 | 1,008 |
| Natural gas | 904 | 326 | 0 | 1 | 577 | 904 |

- **No price** is the only large cause. Household price series start between 1985
  and 2005 for electricity and between 1985 and 2012 for gas (Table M.16), and two
  countries have no gas price at all (Finland, Norway). Italy's electricity price is
  missing in 2007, which makes the only internal gap caused by price. Poland's gas
  price ends in 2023.
- **Missing weather.** Malta's calendar-year HDD (`hdd_calendar`) is missing in 2014
  and 2016. These two rows are the only difference between the electricity panel and
  its weather-only sample.
- **Singleton.** Germany is the only country with a gas price in 1990, so the 1990
  year effect would rest on one observation. That row drops, and the gas benchmark
  starts in 1991.
- **Countries or years without an outcome.** Switzerland has no row in the Eurostat
  outcome file. Cyprus and Malta record 0 ktoe of household gas in all 35 years.
  Bulgaria (1990–1999), Greece (1990–1998), Norway (1990–1999) and Portugal
  (1990–1996) record 0 before their gas series start, and a zero cannot enter the log
  outcome. The United Kingdom's outcome series ends in 2019 for both fuels.

**Table M.16** — Year in which each country enters the benchmark sample.

| Entry year | Electricity | Natural gas |
| --- | --- | --- |
| 1990 | BEL, DNK, LUX | — |
| 1991 | DEU, ESP, FRA, GBR, GRC, IRL, ITA, MLT, NLD, PRT | DEU, GBR |
| 1992 | HUN, SVN | HUN |
| 1995 | FIN, NOR | SVN |
| 1996 | AUT, SWE | SWE |
| 1999 | CYP | DNK |
| 2000 | CZE, POL | CZE, POL |
| 2002 | EST | PRT |
| 2003 | — | EST |
| 2004 | BGR, LTU, LVA, SVK | AUT, BGR, FRA, IRL, ITA, LTU, LVA, NLD, SVK |
| 2005 | HRV, ROU | BEL, ESP, HRV, LUX, ROU |
| 2012 | — | GRC |

The electricity benchmark holds 3 countries in 1990 and 13 in 1991, and all 29 by
2005. The gas benchmark is mostly a post-2004 panel: 506 of its 577 rows (88%) are
from 2004–2024, when 19–25 countries are present, while before 2004 it holds 2–10
countries a year.

### M.4.4 Consequences for estimation

- **No country is short.** The benchmark observes every country for at least 20 years
  (electricity) or 13 years (gas: Greece, 2012–2024), far more than the three
  parameters of the country fixed effect and trends. All 29 (electricity) and 25
  (gas) country clusters carry identifying variation.
- **Selection.** The variable-addition test of §M.2.7 does not reject for either fuel
  (Table M.17). Its power is low here, though. s(t−1) = 0 occurs almost only in the
  first year of a country's price series: 26 of the 29 electricity rows (the
  exceptions follow an internal gap: ITA 2008, MLT 2015, MLT 2017) and all 23 gas
  rows. s(t+1) = 0 occurs in only 4 and 2 rows. In practice the test asks whether the
  first year of a price series is unusual.
- **Separating the effect of the imbalance.** The imbalance of the energy benchmark
  is a by-product of the price control. The "No price, same rows" and "No price,
  full sample" columns differ only in their rows: the second uses the close-to-balanced
  weather-only sample (1,008 of 1,015 cells for electricity, 904 of 945 for gas).
  Comparing the two columns therefore isolates the effect of the price-driven loss
  of rows.

**Table M.17** — Variable-addition test on the energy benchmark (country-clustered).
Coefficients are in log points, with p-values in brackets.

| Fuel | Rows | Rows with s(t−1) = 0 / s(t+1) = 0 | s(t−1) | s(t+1) | Joint p |
| --- | --- | --- | --- | --- | --- |
| Electricity | 813 | 29 / 4 | +0.011 (0.503) | +0.052 (0.173) | 0.368 |
| Natural gas | 552 | 23 / 2 | +0.035 (0.457) | −0.060 (0.088) | 0.200 |

---

## M.5 Summary

- **No estimation panel is balanced.** For the crops (spec F = spec A) the fill is
  0.52–0.71 and at most 23% of the regions are observed every year. The mixed panel's
  fill is 0.38–0.54. The energy benchmark's fill is 0.832 (electricity) and 0.679
  (gas); the energy weather-only samples reach 0.993 and 0.957.
- **The imbalance is mostly staggered entry.** For the crops, 65–83% of the empty
  cells come before a region's first year. The first years are those of the source:
  France from 1989, the other countries between 1990 and 2017, 21 of 71 country
  series in 1998. For energy, 163 of 171 (electricity) and 267 of 273 (gas) empty
  cells come before a country's first year, which is set by when its household price
  series starts.
- **The pipeline removes little.** For the crops it removes 4 internal region-years,
  five island regions without weather, and the single-observation regions. For energy
  it removes two Maltese rows without HDD and one singleton. Internal gaps cover
  1.1–9.1% of the in-span region-years, and almost all of them come from the source.
- **No evidence of selection on shocks.** The variable-addition tests show no sign
  that selection is related to yield or consumption shocks. Weather predicts crop
  observation only weakly (at most 2.8 percentage points per standard deviation), and
  that is selection on the regressors.
- **Points to keep in mind.**
  - Soft wheat's Slovak cluster carries no information, so 14 clusters are
    informative, not 15.
  - Regions with 2–3 years are fitted exactly (0.07–0.94% of the rows).
  - In the first and last years the crop cross-section is France alone.
  - German maize and sunflower (reported at NUTS1) are in neither crop panel.
  - The gas benchmark is mostly a 2004–2024 panel.

**References.** Nijman, T. and Verbeek, M. (1992), "Nonresponse in panel data: the
impact on estimates of a life cycle consumption function", *Journal of Applied
Econometrics* 7(3), 243–257. Wooldridge, J. M. (2010), *Econometric Analysis of Cross
Section and Panel Data*, 2nd ed., MIT Press.
