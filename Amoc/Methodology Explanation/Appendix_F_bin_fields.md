# Appendix F — Per-Bin Monthly AMOC Effect Fields

This appendix documents how the AMOC-weakening signal is resolved **by level of
weakening** rather than as a single end-state, and how the resulting fields are
constructed, validated, and bounded. The object built here is a set of monthly
perturbation fields indexed by AMOC strength: for each model and each 1 Sv bin of
AMOC weakening, a per-calendar-month, cell-by-cell delta for `tas`, `tasmin`,
`tasmax`, and a multiplicative factor for `pr`. These fields are the input to the
replay of Appendix G.

---

## F.1 Why a bin decomposition, and why the fields had to be rebuilt

The impact question is not "what happens if the AMOC collapses" but "how does the
impact vary with the degree of weakening". That requires the forcing to be indexed
by AMOC strength, not reduced to the terminal state of each hosing run.

A set of per-bin fields already existed (`amoc_effect_bins_<MODEL>.nc`). They could
not be used, for three reasons that are properties of those files and not of the
underlying simulations:

1. **The month axis had already been collapsed.** The arrays are dimensioned
   `(bin, lat, lon)`. A single annual number per cell cannot express the seasonal
   structure of the AMOC response, which is strongly winter-weighted.
2. **Only `tas` and `pr` were retained.** The indicator definitions require all
   three temperature fields: `heat_daily` reads the daily maximum, `frost_daily`
   the daily minimum, and the within-day integration of the energy indicators uses
   the diurnal amplitude `(tx − tn)/2`. Applying a single `tas` delta to all three
   would fix the diurnal range by construction and suppress precisely the
   asymmetry the indicators are built to detect.
3. **Precipitation was stored as an additive anomaly in mm day⁻¹.** The coupling
   uses a multiplicative factor, which cannot go negative in dry cells and scales
   with local climate.

The fields were therefore rebuilt from the monthly anomaly cubes in
`anomaly_output/`, which retain every year and every calendar month for all four
variables. Script: `Amoc/Code/amoc_bin_fields.R`.

---

## F.2 Definition of the bins

For each hosing year *y* the AMOC state is taken from the M26 series (AMOC
streamfunction maximum at 26° N):

$$\Delta \mathrm{Sv}(y) = \mathrm{hos}(y) - \overline{\mathrm{con}}$$

where `hos` is the hosing run and `con` the control, averaged over its full length.
Bins are 1 Sv wide, assigned by `floor(ΔSv)`, and a bin is retained only if it
contains at least three hosing years.

Two restrictions follow from the data rather than from choice:

- **Only the `u03-hos` protocol can be binned.** M26 archives no Sv series for
  `g01`, so no AMOC state can be attached to its years. The `g01` runs therefore
  appear only in the end-state comparison of Appendix G, never in the bin
  decomposition.
- **Only the first `min(valid M26 years, cube years)` hosing years are usable.**
  For HadGEM3-GC3-1LL the M26 series has 145 valid years while the anomaly cube has
  100, so 100 are used. For HadGEM3-GC3-1MM the binding constraint is M26, with 99
  valid years.

The resulting bin structure:

| Model | Bins | ΔSv range covered | Years per bin |
|---|---|---|---|
| IPSL-CM6A-LR | 10 | −10 … 0 | 5, 8, 19, 11, 13, 11, 8, 9, 6, 8 |
| EC-Earth3 | 10 | −10 … 0 | 4, 15, 20, 21, 8, 8, 9, 5, 6, 3 |
| HadGEM3-GC3-1LL | 9 | −9 … 0 | 12, 19, 24, 10, 6, 7, 6, 8, 6 |
| HadGEM3-GC3-1MM | 11 | −15 … −4 | 4, 8, 13, 21, 10, 7, 7, 8, 4, 4, 4 |

Two features of this table matter downstream. Bin occupancy is highly uneven —
from 3 to 24 years — so the extreme bins rest on few years and carry
correspondingly wide dispersion (Appendix I). And **HadGEM3-GC3-1MM never reaches
weak weakening**: its range stops at −4 Sv, so there is no common ΔSv axis across
all four models at the weak end, and cross-model comparison is only meaningful
where the ranges overlap.

---

## F.3 Construction of the fields

For each bin, each calendar month, and each cell, the field is the mean over the
bin's hosing years of the monthly anomaly for that month. For temperature the
anomaly is additive (K); for precipitation the ratio file
(hosing / piControl, computed by `ymondiv`) is averaged, giving a dimensionless
multiplicative factor.

**Time indexing is done from the time *values*, never from the position.** This is
not defensive coding but a response to a defect found in the source data: for
HadGEM3-GC3-1MM, `tasmin` and `tasmax` carry 1197 time steps against 1200 for
`tas`, with identical first and last values. Three months are missing from the
middle of the series — July, August and September of hosing year 38. Positional
indexing `(y−1)·12 + m` would have shifted every month after the gap into the wrong
calendar month and the wrong bin, silently. Hosing year 38 falls in bin
[−10, −9), so for those three months that bin averages six years instead of seven.
This is the only such gap across the four models.

Two calendar conventions are in play and are decoded accordingly: proleptic
Gregorian for IPSL-CM6A-LR and EC-Earth3, `360_day` for both HadGEM configurations.
A naive days-since decode applied to the 360-day files is wrong by roughly two
centuries.

**Grid and extent.** Fields are kept on each model's native grid, cropped to the
E-OBS extent plus 3° (longitude −44 … 79, latitude 22 … 79) and with longitude
rewrapped from 0…360 to −180…180 and re-sorted. The margin guarantees four native
neighbours for every E-OBS cell centre in the bilinear interpolation of Appendix G.
No regridding is performed at this stage: interpolation happens once, at the point
of use.

**Non-finite precipitation ratios.** The `ymondiv` ratio is undefined where control
precipitation is zero. Across the European box, IPSL-CM6A-LR has 1 764 non-finite
values over 100 years (0.044 % of the field), confined to 15 cells all lying
between 22° N and 29° N — the Saharan and Arabian margins. None of these cells lies
within one model grid cell of any E-OBS cell used by either impact branch, so their
treatment is immaterial to every number reported. The other three models have none.
Non-finite values are excluded from the bin mean; the [0.1, 10] bound is applied by
the consumer, after interpolation, so that bin-mean fields and the per-year fields
of the uncertainty band receive identical treatment.

---

## F.4 Validation

Because the script that produced the pre-existing annual fields is not in the
repository, the rebuild is validated against those files on two independent axes.

**Bin membership.** Reconstructing ΔSv from M26 and applying the binning rule
reproduces the `n_years` vector of all four `amoc_effect_bins_*.nc` files
**exactly**, including the drop of bins with fewer than three years. Since the bin
counts are a deterministic function of the year-to-bin assignment, an exact match on
all 40 bins across four models establishes that the rebuild bins the same years as
the original.

**Field content.** Averaging the new monthly fields over the twelve months
reproduces the annual field of the existing files to within **1 × 10⁻⁶ K** on every
bin of every model — the float32 storage precision of the older files. The rebuild
therefore adds the month axis, the two extra temperature fields and the
multiplicative precipitation factor **without changing the signal** that was already
there.

Total build time is 11 s for all four models.

---

## F.5 Per-year fields for the uncertainty band

The same builder, run in `mode = "year"`, emits one field per individual hosing
year instead of one per bin: 98 for IPSL-CM6A-LR and 99 for EC-Earth3, i.e. every
year belonging to a retained bin. Reading, time indexing, cropping and the
non-finite rule are the identical code path; only the grouping of years differs.
This is deliberate — it removes the possibility that the uncertainty band of
Appendix I disagrees with the central estimate for reasons of implementation rather
than of substance.

---

## F.6 What this appendix does not establish

The ΔSv definition, `hos(y) − mean(con)`, is **reproduced** from the existing
fields, not independently validated. The evidence in §F.4 shows that the rebuild
assigns years to bins exactly as the previous construction did; it says nothing
about whether that definition is the appropriate measure of AMOC state. If the
definition is inappropriate, this appendix has faithfully reproduced the problem.

The bin decomposition also treats ΔSv as a sufficient statistic for the climate
state: two years with the same ΔSv are assumed exchangeable. Nothing here tests
that assumption, and the within-bin dispersion reported in Appendix I is the
closest available evidence on it.
