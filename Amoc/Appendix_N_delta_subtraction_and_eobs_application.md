# Appendix N — From NAHosMIP to E-OBS: the subtraction and the daily application

Two operations turn the hosing experiment into weather that the response functions can read:
a **subtraction** (hosing minus piControl, then averaged by level of AMOC weakening) and an
**application** (the resulting monthly change added to every observed E-OBS day before any
indicator is computed). This appendix states both as the code implements them. The file and line
references are in §N.4, and the numbers in the tables are read from the data files.

---

## N.0 Notation

| symbol | meaning |
|---|---|
| $c$ / $g$ | native model grid cell / E-OBS 0.25° cell |
| $h$, $k$, $y$ | hosing year, piControl year, E-OBS (replay) year |
| $m$, $d$, $m_d$ | calendar month, day, calendar month of day $d$ |
| $b$ | $\Delta$Sv bin (1 Sv wide) |
| $X_v$, $P$ | monthly mean of $v\in\{\text{tas},\text{tasmin},\text{tasmax}\}$ (K); monthly precipitation |
| $M$ | AMOC strength index (NAHosMIP M26 files: `hos_<model>`, `con_<model>`), Sv |

---

## N.1 The subtraction: hosing minus piControl

**Step 1 — piControl climatology** (`cdo ymonmean`): 12 values per native cell, one per calendar month,
averaged over the whole reference window of $K$ years.

$$\bar X^{\mathrm{pi}}_v(m,c)=\frac1K\sum_{k=1}^{K}X^{\mathrm{pi}}_v(k,m,c),\qquad
\bar P^{\mathrm{pi}}(m,c)=\frac1K\sum_{k=1}^{K}P^{\mathrm{pi}}(k,m,c)$$

**Step 2 — anomaly of every hosing month** (`cdo ymonsub`, `cdo ymondiv`). Each month of the hosing run
minus (or, for precipitation, divided by) the climatology of the same calendar month. There is no
pairing with a parallel control year.

$$a_v(h,m,c)=X^{\mathrm{hos}}_v(h,m,c)-\bar X^{\mathrm{pi}}_v(m,c),\qquad
r(h,m,c)=\frac{P^{\mathrm{hos}}(h,m,c)}{\bar P^{\mathrm{pi}}(m,c)}$$

The anomaly lives on the native hosing grid. The climatology is put on that grid by a coordinate copy
when the two grids have the same size, and by a real remap (bilinear for temperature, conservative
for precipitation) only if they differ.

**Step 3 — AMOC weakening of each hosing year, and its bin.**

$$\Delta S(h)=M^{\mathrm{hos}}(h)-\frac1{K_M}\sum_{k=1}^{K_M}M^{\mathrm{pi}}(k),\qquad b(h)=\lfloor\Delta S(h)\rfloor$$

- Hosing years used: $h=1,\dots,n$, with $n=\min(\text{valid M26 hosing years},\ \text{anomaly years})$.
- Bins are 1 Sv wide, labelled $b+0.5$, and kept only if $|Y_b|\ge3$, where $Y_b=\{h\le n:\ b(h)=b\}$.
- A bin's level is $\overline{\Delta S}_b=\tfrac1{|Y_b|}\sum_{h\in Y_b}\Delta S(h)$.

**Step 4 — bin delta fields**: the mean, over the bin's hosing years, of the monthly anomalies.

$$\delta_{v,b}(m,c)=\frac1{|Y_b|}\sum_{h\in Y_b}a_v(h,m,c),\qquad
R_b(m,c)=\frac1{|Y_b|}\sum_{h\in Y_b}r(h,m,c)$$

$R_b$ is a mean of ratios. Non-finite values, such as $r=\infty$ where control precipitation is zero,
are dropped from the mean. Output: `amoc_bin_fields_<MODEL>_u03.nc`, dims (lon, lat, month, bin),
native grid cropped to 44°W–79°E, 22–79°N.

**Table N.1 — reference windows and bins** (NAHosMIP u03-hos)

| model | native grid (lon × lat) | climate reference (Step 1) | AMOC reference $K_M$, mean | hosing years used | $\Delta S$ range (Sv) | bins kept | years per bin | years in bins |
|---|---|---|---|---|---|---|---|---|
| IPSL-CM6A-LR | 2.50° × 1.27° | piControl 1850–2349 (500 yr) | 500 yr, 12.48 Sv | 100 | −10.6 … +0.2 | −0.5 … −9.5 (10) | 5–19 | 98 |
| EC-Earth3 | 0.70° × 0.70° | piControl 2259–2759 (501 yr; tasmax 499) | 150 yr, 17.24 Sv | 100 | −9.4 … 0.0 | −0.5 … −9.5 (10) | 3–21 | 99 |
| HadGEM3-GC3.1-LL | 1.88° × 1.25° | piControl 1850–1949 (100 yr) | 150 yr, 15.37 Sv | 100 | −9.3 … +0.6 | −0.5 … −8.5 (9) | 6–24 | 98 |
| HadGEM3-GC3.1-MM | 0.83° × 0.56° | piControl 1850–1949 (100 yr) | 150 yr, 16.33 Sv | 99 | −14.8 … +1.0 | −4.5 … −14.5 (11) | 4–21 | 90 |

For IPSL-CM6A-LR and EC-Earth3, tasmin and tasmax on the hosing side are reconstructed fields
(per-cell, per-month OLS fitted on piControl). HadGEM3 provides them natively.

---

## N.2 The application: summing the change onto E-OBS days

**Step 5 — onto the E-OBS cells.** Each (variable, month, bin) field is interpolated bilinearly to the
E-OBS cell centres. Outside the model domain the default is "no change" ($\tilde\delta=0$, $\tilde R=1$).
The ratio is then clamped: $\tilde R_b\leftarrow\min(\max(\tilde R_b,0.1),10)$.

**Step 6 — every day of every replay year** (additive on temperature, multiplicative on precipitation):

$$tg'_{g,d}=tg_{g,d}+\tilde\delta_{\mathrm{tas},b}(m_d,g),\quad
tx'_{g,d}=tx_{g,d}+\tilde\delta_{\mathrm{tasmax},b}(m_d,g),\quad
tn'_{g,d}=tn_{g,d}+\tilde\delta_{\mathrm{tasmin},b}(m_d,g),\quad
rr'_{g,d}=rr_{g,d}\,\tilde R_b(m_d,g)$$

The same 12 monthly values per cell are added to **every** day of their month in **every** replay
year. No hosing year is matched to an E-OBS year: each bin replays the whole observed record, with
its own daily and interannual variability, shifted by the bin's mean change.

**Table N.2 — replay domains**

| branch | E-OBS cells | years replayed | extra year read | variables |
|---|---|---|---|---|
| crop | 12,150 (NUTS3 crosswalk) | 1989–2023 | 1988 (first overwinter block) | tg, tx, tn, rr |
| energy | 9,309 (populated cells, 33 countries) | 1990–2024 | 1989 (Oct–Dec of the first Oct–Mar season) | tg, tx, tn |

**Step 7 — daily indicators, computed on the shifted days** (thresholds act day by day):

| indicator | daily value $i_{g,d}$ |
|---|---|
| GDD | $\max(\min(tg',28)-5,\ 0)$ |
| heat | $\max(tx'-28,\ 0)$ |
| frost | $\mathbf 1[tn'<0]$ |
| precipitation | $rr'$ |
| HDD | $\tfrac1{24}\sum_{s=0}^{23}\mathbf 1[T_s<15]\,(18-T_s)$ |
| CDD | $\tfrac1{24}\sum_{s=0}^{23}\mathbf 1[T_s\ge24]\,(T_s-21)$ |

where $T_s=tg'+\tfrac{tx'-tn'}{2}\sin(2\pi s/24)$: a within-day cycle sampled 24 times.

**Step 8 — sums, in the order cell → month → area → block or season**

$$I_{g,m,y}=\sum_{d\in(m,y)}i_{g,d}\qquad(\text{NA if any day of the month is missing})$$

- **Crops**, for NUTS3 region $r$ with crosswalk weights $\omega_{rg}$: first the regional monthly value,
  renormalised over the cells with a complete month and rounded to 3 dp, then the sum over the months
  of block $B$.

$$X_{r,m,y}=\frac{\sum_g\omega_{rg}\,I_{g,m,y}}{\sum_g\omega_{rg}\,\mathbf 1[I_{g,m,y}\ \text{finite}]},\qquad
X^{B}_{r,y}=\sum_{m\in B}X_{r,m,y}$$

  The squared precipitation term is cell-first: $\sum_g\omega_{rg}\big(P^{B}_{g,y}\big)^2\big/\sum_g\omega_{rg}$,
  never $(X^{B}_{r,y})^2$.

- **Energy**, for country $j$ with population weights $\pi_{jg}$: first the season sum per cell,
  $S\in\{$Jan–Dec (HDD), Oct$_{y-1}$–Mar$_y$ (HDD), Jun–Aug (CDD)$\}$, then the population-weighted mean.

$$X^{S}_{j,y}=\frac{\sum_g\pi_{jg}\sum_{m\in S}I_{g,m,y}}{\sum_g\pi_{jg}\,\mathbf 1[\,\cdot\ \text{finite}]}$$

**Table N.3 — crop blocks (spec F, $y$ = harvest year)**

| crop | overwinter block $W$ | season block $S$ |
|---|---|---|
| soft wheat (autumn-sown) | Oct$_{y-1}$ – Feb$_y$ | Mar – Aug |
| soft wheat (spring-sown, northern regions) | — | Mar – Aug |
| winter barley | Oct$_{y-1}$ – Feb$_y$ | Mar – Aug |
| durum wheat | Nov$_{y-1}$ – Feb$_y$ | Mar – Aug |
| spring barley | — | Feb – Aug |
| grain maize | — | Mar – Nov |
| sugar beet, sunflower | — | Mar – Sep |

**Step 9 — gate, and what is differenced.** With $\tilde\delta=0$ and $\tilde R=1$ the replay reproduces
the estimation files exactly (max |difference| = 0), and this is checked before any scenario is run.
Downstream, each replay year is compared with the same observed year,
$\Delta X_{r,y}=X^{\text{scen}}_{r,y}-X^{\text{obs}}_{r,y}$, and $\beta'\Delta X_{r,y}$ is averaged over the replay years.

---

## N.3 What the construction implies

| feature of the code | consequence |
|---|---|
| Reference = long piControl climatology, not a parallel control segment | $\delta$ = mean hosing state minus mean unforced state. Any piControl drift between the two periods enters $\delta$; it is not measured. |
| Climate and AMOC references come from different piControl series and windows (Table N.1) | $\Delta S$ and $\delta$ are each referenced to their own control mean. |
| Bins hold 3–24 hosing years | A bin mean still carries internal variability, most in the 3–4-year bins (EC-Earth3 −0.5 and −9.5; HadGEM3-MM −4.5, −5.5, −6.5 and −14.5). |
| All hosing-induced change is grouped by $\Delta S$ | $\delta$ attributes to the AMOC everything the freshwater forcing changes. |
| One $\delta$ for all replay years | The baseline is the observed climate of 1989–2023 (energy 1990–2024), trend included. Near-linear indicators are insensitive to it, threshold indicators are not. At −4.5 Sv (mean of the 4 models) the shift applied to 1989–98 vs 2014–23 is −185 vs −191 dd for season GDD, 0 to +4% for frost, +2% for precipitation, −3%/0% for HDD — but −65 vs −76 dd for overwinter GDD (+17%), −21.7 vs −26.8 dd for CDD (+23%), and −14.6 vs −20.5 dd for heat (+40%). Crop figures are soft-wheat NUTS3 means; energy figures are country means. |

---

## N.4 Code map

| step | file : lines |
|---|---|
| 1–2 piControl climatology, anomaly, ratio | `Code/anomaly/1.IPSLnahos_anomaly_cdo_explicit.sh` : 45–48, 66, 99, 110 · `1.ECHEarth3nahos_anomaly_cdo_explicit.sh` : 50–53, 72, 113, 129 · `nahos_anomaly_cdo (3).sh` : 47–48, 145, 148, 168 |
| 2 remap only on grid-size mismatch | `1.IPSLnahos_anomaly_cdo_explicit.sh` : 50–53, 79–94 |
| tasmin/tasmax reconstruction (IPSL, EC-Earth3) | `Code/anomaly/1.reconstruct_tasmin_tasmax_{IPSL,ECEarth3}.py` |
| 3 $\Delta S$, bins | `Code/amoc_bin_fields.R` : 60–72 (box : 31) |
| 4 bin means | `Code/amoc_bin_fields.R` : 111–121 |
| 5 interpolation, defaults, clamp | `Code/scenario_replay_bins.R` : 39–43 · `Code/scenario_engine.R` : 25, 151–155 |
| 6 daily perturbation | `Code/scenario_engine.R` : 192–197 (temperature), 241 (precipitation) |
| 7 daily indicators | `Code/weather_indicators.R` : 14–19, 23–24, 27–29 |
| 8 aggregation | `Code/scenario_engine.R` : 169–175 (crops), 198–218 (energy) · `Code/weather_indicators.R` : 58–65, 87 (blocks) · `Code/crop_spec_f_blocks.R` (cell-first square) |
| 9 zero-delta gate | `Code/scenario_engine.R` : 346–348 |
