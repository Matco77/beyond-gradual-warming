# Appendix N — From NAHosMIP to E-OBS: the subtraction and the daily application

The hosing experiment reaches the response functions in two moves:

1. **Subtraction** (model world, monthly, native grid): what the hosing run changes relative to piControl, grouped by how much the AMOC has weakened.
2. **Application** (observed world, daily, E-OBS 0.25°): that monthly change added to every observed day, before any indicator is computed.

Everything below follows the code (file and line references in §N.6). The numbers are read from the data files.

---

## N.0 The chain in one picture

```
 MODEL WORLD  (NAHosMIP u03-hos + piControl · native grid · MONTHLY)
 ─────────────────────────────────────────────────────────────────────────────────────────────
   piControl months ──[mean of each calendar month]──────────────►  climatology      12 per cell
                                                                          │
   hosing months ─────[minus that month's climatology (T);               │
                        divided by it (P)]───────────────────────►  anomaly of EVERY hosing month
                                                                          │
   AMOC index ────────[hosing year − mean piControl AMOC]─►  ΔSv(year) ─[1-Sv bin, ≥3 years]
                                                                          │
                                        [mean of the anomalies over the bin's years, month by month]
                                                                          ▼
                                                                 δ_bin(month, cell)  12 per cell per bin
 ─────────────────────────────────────────────────────────────────────────────────────────────
 OBSERVED WORLD  (E-OBS · 0.25° cells · DAILY · 1989–2023 crops, 1990–2024 energy)
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                        [bilinear to E-OBS cell centres; ratio clamped to 0.1–10]
                                                                          ▼
   observed day (tg, tx, tn, rr) ──────────────────────────────►  (+) δ of the day's month
                                     EVERY day, EVERY year              (×) for precipitation
                                                                          ▼
                                                                 shifted day
                                                                          │ [daily indicator: the thresholds
                                                                          │  act on each shifted day]
                                                                          ▼
                                             GDD · heat · frost · precip  |  HDD · CDD   per cell, per day
                                                                          │ [sum over the days of the month]
                                                                          ▼
                                                                 per cell, per month
                                                                          │ [area mean → NUTS3 · population mean → country]
                                                                          │ [sum over the months of the crop block / energy season]
                                                                          ▼
                                                   X_scen(region, year)   ── vs ──   X_obs(region, year)
                                                                                   (same chain, δ = 0;
                                                                                    identical to the
                                                                                    estimation files)
                                                                          ▼
                                             ΔX = X_scen − X_obs, year by year  →  β·ΔX  →  mean over years
```

---

## N.1 Step by step

| # | goes in | operation | comes out |
|---|---|---|---|
| 1 | piControl monthly fields, $K$ years | mean of each calendar month (`cdo ymonmean`) | climatology: 12 values per model cell |
| 2 | hosing monthly fields | each month **minus** the climatology of the same month (T), **divided** by it (P) (`ymonsub` / `ymondiv`) | an anomaly for **every** hosing month |
| 3 | AMOC index, hosing and piControl | $\Delta$Sv = hosing year − mean of the whole piControl AMOC series; 1-Sv bins by floor; keep bins with ≥ 3 years | the bin of each hosing year |
| 4 | anomalies + bins | mean over the bin's years, month by month | $\delta_b$: 12 values per model cell per bin |
| 5 | $\delta_b$ | bilinear to the E-OBS cell centres; outside the model domain "no change"; precipitation ratio clamped to [0.1, 10] | $\tilde\delta_b$: 12 values per E-OBS cell per bin |
| 6 | observed E-OBS days + $\tilde\delta_b$ | $\tilde\delta$ of the day's month **added** (T) or **multiplied** (P), on **every** day of **every** replay year | shifted daily weather |
| 7 | shifted days | daily indicator, threshold evaluated day by day | one value per cell per day |
| 8 | daily values | Σ days → month (per cell); weighted mean → NUTS3 (area) or country (population); Σ months → crop block or energy season | $X^{\text{scen}}$(region, year) |
| 9 | the same chain with $\delta=0$ | — | $X^{\text{obs}}$(region, year) = the estimation files exactly (gate) |
| 10 | $X^{\text{scen}}$, $X^{\text{obs}}$ | difference year by year, times $\beta$, mean over the replay years | impact |

Three things the table makes explicit:

- **Nothing is paired in time.** Hosing months are compared with a *climatology*, not with a parallel control year.
- **The change is a bin average.** It does not come from one hosing year.
- **The observed record is kept whole.** Each bin replays all 35 observed years, with their own weather, shifted by the same 12 monthly values.

---

## N.2 A worked example: one cell, one month (illustrative numbers)

| stage | numbers |
|---|---|
| piControl July mean, one model cell (Step 1) | daily maximum (tasmax) 24.0 °C · precipitation 60 mm |
| three hosing years with $\Delta$Sv = −4.2, −4.6, −4.9 → floor = −5 → bin **−4.5** (Step 3) | their July means: tasmax 22.1, 22.6, 22.2 °C · 48, 57, 54 mm |
| monthly anomalies (Step 2) | −1.9, −1.4, −1.8 °C · ratios 0.80, 0.95, 0.90 |
| bin delta for July (Step 4) | $\delta_{\text{tasmax}}$ = **−1.7 °C** · $R$ = **0.88** |
| E-OBS cell inside that model cell (Step 5) | $\tilde\delta$ ≈ −1.7 °C, $\tilde R$ ≈ 0.88 |
| 10 July 2003: tx 35.0 °C, rr 10.0 mm (Step 6) | tx **33.3 °C**, rr **8.8 mm** · heat 7.0 → 5.3 degree-days |
| 10 July 1995: tx 29.0 °C (Step 6–7) | tx **27.3 °C** · heat 1.0 → **0** (the day falls below 28 °C) |
| July of that cell (Step 8) | sum of the 31 shifted days, for every year 1989–2023, then the regional and block sums |

The daily mean (tg, with $\delta_{\text{tas}}$) and the daily minimum (tn, with $\delta_{\text{tasmin}}$) go
through the same steps, each with its own delta.

The second shifted day shows why the change is applied **before** the indicators. On a warm July
day the same −1.7 °C removes part of the heat; on a milder one it removes all of it. A change applied
to a monthly total could not tell the two apart.

---

## N.3 The formulas

$c$ model cell, $g$ E-OBS cell, $h$ hosing year, $k$ piControl year, $y$ replay year, $m$ month, $d$ day ($m_d$ its month), $b$ bin.

| step | formula |
|---|---|
| 1 | $\bar X^{\mathrm{pi}}_v(m,c)=\frac1K\sum_{k}X^{\mathrm{pi}}_v(k,m,c)$, same for $\bar P^{\mathrm{pi}}$; $v\in\{$tas, tasmin, tasmax$\}$ |
| 2 | $a_v(h,m,c)=X^{\mathrm{hos}}_v(h,m,c)-\bar X^{\mathrm{pi}}_v(m,c)$, $\quad r(h,m,c)=P^{\mathrm{hos}}(h,m,c)/\bar P^{\mathrm{pi}}(m,c)$ |
| 3 | $\Delta S(h)=M^{\mathrm{hos}}(h)-\frac1{K_M}\sum_k M^{\mathrm{pi}}(k)$, $\quad b(h)=\lfloor\Delta S(h)\rfloor$, $\quad Y_b=\{h: b(h)=b\}$, kept if $\lvert Y_b\rvert\ge3$ |
| 4 | $\delta_{v,b}(m,c)=\frac1{\lvert Y_b\rvert}\sum_{h\in Y_b}a_v(h,m,c)$, $\quad R_b(m,c)=\frac1{\lvert Y_b\rvert}\sum_{h\in Y_b}r(h,m,c)$ (non-finite ratios dropped) |
| 5 | $\tilde\delta_{v,b}(m,g),\ \tilde R_b(m,g)$ = bilinear at $g$; outside: 0 and 1; $\tilde R_b\leftarrow\min(\max(\tilde R_b,0.1),10)$ |
| 6 | $tg'_{g,d}=tg_{g,d}+\tilde\delta_{\mathrm{tas},b}(m_d,g)$ (tx with tasmax, tn with tasmin), $\quad rr'_{g,d}=rr_{g,d}\,\tilde R_b(m_d,g)$ |
| 7 | GDD $\max(\min(tg',28)-5,0)$ · heat $\max(tx'-28,0)$ · frost $\mathbf 1[tn'<0]$ · precip $rr'$ · HDD $\tfrac1{24}\sum_s\mathbf 1[T_s<15](18-T_s)$ · CDD $\tfrac1{24}\sum_s\mathbf 1[T_s\ge24](T_s-21)$, with $T_s=tg'+\tfrac{tx'-tn'}{2}\sin\tfrac{2\pi s}{24}$ |
| 8 | $I_{g,m,y}=\sum_{d\in(m,y)}i_{g,d}$; crops: $X_{r,m,y}=\frac{\sum_g\omega_{rg}I_{g,m,y}}{\sum_g\omega_{rg}\mathbf 1[\text{finite}]}$ (3 dp), $X^B_{r,y}=\sum_{m\in B}X_{r,m,y}$, precip² cell-first $\sum_g\omega_{rg}(P^B_{g,y})^2/\sum_g\omega_{rg}$; energy: $X^S_{j,y}=\frac{\sum_g\pi_{jg}\sum_{m\in S}I_{g,m,y}}{\sum_g\pi_{jg}\mathbf 1[\text{finite}]}$ |
| 9–10 | gate: $\tilde\delta=0,\tilde R=1\Rightarrow X^{\text{scen}}\equiv X^{\text{obs}}$; $\quad\Delta X_{r,y}=X^{\text{scen}}_{r,y}-X^{\text{obs}}_{r,y}$, then $\beta'\Delta X$ averaged over $y$ |

---

## N.4 The numbers behind each step

**Table N.1 — references and bins, per model** (u03-hos)

| model | native grid | climate reference (step 1) | AMOC reference (step 3) | hosing years used | $\Delta$Sv range | bins kept | years per bin | years in bins |
|---|---|---|---|---|---|---|---|---|
| IPSL-CM6A-LR | 2.50° × 1.27° | piControl 1850–2349, 500 yr | 500 yr, mean 12.48 Sv | 100 | −10.6 … +0.2 | −0.5 … −9.5 (10) | 5–19 | 98 |
| EC-Earth3 | 0.70° × 0.70° | piControl 2259–2759, 501 yr (tasmax 499) | 150 yr, mean 17.24 Sv | 100 | −9.4 … 0.0 | −0.5 … −9.5 (10) | 3–21 | 99 |
| HadGEM3-GC3.1-LL | 1.88° × 1.25° | piControl 1850–1949, 100 yr | 150 yr, mean 15.37 Sv | 100 | −9.3 … +0.6 | −0.5 … −8.5 (9) | 6–24 | 98 |
| HadGEM3-GC3.1-MM | 0.83° × 0.56° | piControl 1850–1949, 100 yr | 150 yr, mean 16.33 Sv | 99 | −14.8 … +1.0 | −4.5 … −14.5 (11) | 4–21 | 90 |

For IPSL-CM6A-LR and EC-Earth3, tasmin and tasmax on the hosing side are reconstructed (per-cell,
per-month OLS fitted on piControl); HadGEM3 provides them natively. The anomalies live on the native
hosing grid, cropped to 44°W–79°E, 22–79°N. The climatology is remapped (bilinear T, conservative P)
only if its grid size differs; otherwise its coordinates are copied.

**Table N.2 — where and when the change is applied**

| branch | E-OBS cells | replay years | extra year read | variables |
|---|---|---|---|---|
| crop | 12,150 (NUTS3 crosswalk) | 1989–2023 | 1988 (first overwinter block) | tg, tx, tn, rr |
| energy | 9,309 (populated cells, 33 countries) | 1990–2024 | 1989 (Oct–Dec of the first Oct–Mar season) | tg, tx, tn |

**Table N.3 — months summed (step 8)**

| crops (spec F, $y$ = harvest year) | overwinter block | season block |
|---|---|---|
| soft wheat, winter barley (autumn-sown) | Oct$_{y-1}$ – Feb$_y$ | Mar – Aug |
| durum wheat | Nov$_{y-1}$ – Feb$_y$ | Mar – Aug |
| soft wheat (spring-sown, northern regions) | — | Mar – Aug |
| spring barley | — | Feb – Aug |
| grain maize | — | Mar – Nov |
| sugar beet, sunflower | — | Mar – Sep |

| energy | months summed |
|---|---|
| HDD, calendar year (electricity) | Jan – Dec |
| HDD, heating season (gas) | Oct$_{y-1}$ – Mar$_y$ |
| CDD (electricity) | Jun – Aug |

---

## N.5 What the construction implies

| feature | consequence |
|---|---|
| Reference = long piControl climatology, not a parallel control | $\delta$ = mean hosing state minus mean unforced state; any piControl drift between the two periods enters $\delta$ (not measured). |
| AMOC and climate references from different piControl series (Table N.1) | $\Delta$Sv and $\delta$ are each referenced to their own control mean. |
| Bins of 3–24 hosing years | A bin mean still carries internal variability, most in the 3–4-year bins (EC-Earth3 −0.5 and −9.5; HadGEM3-MM −4.5, −5.5, −6.5, −14.5). |
| Hosing years grouped by $\Delta$Sv | Everything the freshwater forcing changes is attributed to the AMOC. |
| One $\delta$ for all replay years | The baseline is the observed climate, trend included. At −4.5 Sv (mean of the 4 models) the shift on 1989–98 vs 2014–23: season GDD −185 vs −191 dd, frost 0 to +4 %, precipitation +2 %, HDD −3 % / 0 % — but overwinter GDD +17 %, CDD −21.7 vs −26.8 dd (+23 %), heat −14.6 vs −20.5 dd (+40 %). Threshold-at-the-warm-end indicators depend on the baseline. |

---

## N.6 Code map

| step | file : lines |
|---|---|
| 1–2 | `Code/anomaly/1.IPSLnahos_anomaly_cdo_explicit.sh` : 45–48, 66, 99, 110 · `1.ECHEarth3nahos_anomaly_cdo_explicit.sh` : 50–53, 72, 113, 129 · `nahos_anomaly_cdo (3).sh` : 47–48, 145, 148, 168 |
| 2 (remap rule) | `1.IPSLnahos_anomaly_cdo_explicit.sh` : 50–53, 79–94 |
| 2 (tasmin/tasmax, IPSL & EC-Earth3) | `Code/anomaly/1.reconstruct_tasmin_tasmax_{IPSL,ECEarth3}.py` |
| 3 | `Code/amoc_bin_fields.R` : 60–72 (crop box : 31) |
| 4 | `Code/amoc_bin_fields.R` : 111–121 |
| 5 | `Code/scenario_replay_bins.R` : 39–43 · `Code/scenario_engine.R` : 25, 151–155 |
| 6 | `Code/scenario_engine.R` : 192–197 (temperature), 241 (precipitation) |
| 7 | `Code/weather_indicators.R` : 14–19, 23–24, 27–29 |
| 8 | `Code/scenario_engine.R` : 169–175 (crops), 198–218 (energy) · `Code/weather_indicators.R` : 58–65, 87 · `Code/crop_spec_f_blocks.R` (cell-first square) |
| 9 | `Code/scenario_engine.R` : 346–348 |
