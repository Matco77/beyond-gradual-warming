# Appendix C - AMOC-hosing anomalies, late-window deltas, precipitation ratios, and diagnostics

This appendix documents how the AMOC-hosing anomaly fields are created and how they are reduced to the delta metrics used in the analysis. The code does not dynamically detect an equilibrium state. The main delta used in the analysis is an operational late-window mean: the mean of the final third of the available monthly anomaly series, with a minimum window of 12 months. The word “plateau” appears in some variable names, comments, and plot labels, but in this appendix I refer to the quantity more precisely as a late-window mean response unless discussing the code labels directly.

The code separates three related but distinct steps:

1. construction of monthly hosing-minus-piControl anomaly fields;
2. construction of late-window delta fields from those additive anomaly fields;
3. construction of diagnostics, including Europe-mean effect ratios, NW-minus-Mediterranean spatial contrasts, and seasonal signal-to-noise fields.

For temperature variables, the anomaly and delta fields are additive changes in kelvin. For precipitation, the code creates both an additive anomaly and a separate multiplicative ratio. The additive precipitation anomaly is used for the delta maps and diagnostics. The multiplicative precipitation ratio is a separate monthly field intended for later precipitation stressing. The uploaded scripts prepare these anomaly, delta, ratio, and diagnostic products; they do not apply the anomalies or ratios to ISIMIP daily data.

## C.1 Monthly anomaly construction

For each model, hosing protocol, variable, grid cell, and monthly time step, the additive anomaly is defined as the hosing monthly value minus the corresponding piControl monthly climatology:

```text
A_v(cell, t) = H_v(cell, t) - C_v(cell, month(t))
```

where:

```text
A_v(cell, t)          = additive anomaly for variable v
H_v(cell, t)          = hosing value for variable v at time t
C_v(cell, month(t))   = piControl climatological mean for the same calendar month
month(t)              = calendar month of the hosing time step
```

Thus, hosing January values are compared with the piControl mean January, hosing February values with the piControl mean February, and so on. The model years themselves do not need to match, because the subtraction is calendar-month matched rather than year-matched.

For HadGEM3-GC31-LL and HadGEM3-GC31-MM, the HadGEM anomaly script loops over:

```text
models:     HadGEM3-GC31-LL, HadGEM3-GC31-MM
protocols:  g01-hos, u03-hos
variables:  tas, tasmin, tasmax, pr
```

For each model-protocol-variable combination, the script discovers all matching hosing monthly chunks, selects piControl files whose filename date ranges overlap the baseline window 185001-194912, merges the hosing files with `cdo mergetime`, merges the selected piControl files with `cdo mergetime`, computes a 12-month piControl climatology with `cdo ymonmean`, and subtracts that climatology from the hosing series with `cdo ymonsub`.

The HadGEM anomaly formula implemented by the script is therefore:

```text
A_v = hosing_monthly_v - piControl_monthly_climatology_v
```

The HadGEM script does not internally remap the piControl climatology to the hosing grid. It assumes that the selected hosing and piControl files are already grid-compatible.

For EC-Earth3, the anomaly script follows the same conceptual method:

```text
A_v = hosing_monthly_v - piControl_monthly_climatology_v
```

but the reference climatology differs. The EC-Earth3 comments state that the hosing experiments branch from a piControl spin-up epoch that is not among the piControl files available on disk. Therefore, the EC-Earth3 script uses the full available EC-Earth3 piControl archive as the reference climatology. The comments identify this archive as model years 2259-2759 for `tas`, `pr`, and `tasmin`, while the available `tasmax` archive ends in 2757. The script enforces this through the available filename patterns rather than through an internal date-slicing command.

For EC-Earth3, native hosing files provide `tas` and `pr`. The `tasmin` and `tasmax` hosing inputs are reconstructed absolute fields produced outside this script, as described in Appendix B. Appendix C does not itself reconstruct `tasmin` or `tasmax`; it uses the reconstructed absolute fields as the hosing-side inputs and then applies the same hosing-minus-piControl anomaly construction used for the native variables.

## C.2 Grid handling

The HadGEM anomaly script assumes that the selected hosing and piControl files for a given model, protocol, and variable are already on compatible grids. It does not perform an explicit grid-size check, `setgrid`, or remapping step before the anomaly is computed.

The EC-Earth3 anomaly script includes explicit grid handling. For each variable and hosing experiment, it compares the CDO `gridsize` of the hosing file with the `gridsize` of the piControl climatology.

If the grid sizes are equal, the script writes the hosing grid description and applies it to the piControl climatology using `cdo setgrid`. This changes the grid metadata or coordinate labels of the climatology so that they match the hosing file, but it does not interpolate the climatology values.

If the grid sizes differ, the script remaps the piControl climatology to the hosing grid. The remapping method depends on the variable:

```text
pr                 -> conservative remapping, remapcon
tas, tasmin, tasmax -> bilinear remapping, remapbil
```

This means that the final EC-Earth3 anomaly is constructed on the hosing grid. A methodological caveat is that equal grid size does not by itself prove identical grid geometry. The `setgrid` branch is a practical grid-label correction, not a geometric proof.

## C.3 Precipitation: additive anomaly versus multiplicative ratio

Precipitation is handled in two different ways, and these outputs must not be confused.

First, the scripts write an additive precipitation anomaly:

```text
A_pr(cell, t) = H_pr(cell, t) - C_pr(cell, month(t))
```

This additive `pr` anomaly is stored in the same way as the temperature anomalies. In the raw anomaly NetCDFs it remains in native CMIP6 precipitation units, `kg m-2 s-1`. When the R scripts read it through `read_europe_cube()`, it is converted to `mm/day` by multiplying by 86400.

Second, for precipitation only, the anomaly-builder scripts also write a multiplicative ratio:

```text
R_pr(cell, t) = H_pr(cell, t) / C_pr(cell, month(t))
```

This ratio is dimensionless and its NetCDF unit is set to `1`. This is the object intended for later multiplicative precipitation stressing, for example:

```text
pr_adjusted_daily = pr_ISIMIP_daily * R_pr
```

Temperature variables do not receive ratio fields, because a ratio of interval-scale temperature in kelvin or degrees Celsius anomalies is not physically meaningful for this workflow.

The late-window delta file for precipitation is not the ratio. The delta script reads the additive precipitation anomaly, converts it to `mm/day`, and averages the final third of the additive anomaly series. Therefore:

```text
delta_pr_lastthird = mean of additive pr anomaly over the late window
```

and its unit is:

```text
mm/day
```

It is not:

```text
mean(R_pr)
```

The plotting script may also show precipitation as a relative percentage of the control annual-mean precipitation:

```text
relative_pr_delta_percent = 100 * delta_pr_lastthird / piControl_annual_mean_pr
```

This percentage is a diagnostic display of the additive delta relative to the control climatology. It is not the same as the monthly multiplicative `R_pr` field.

## C.4 Regional reading, unit conversion, and area weighting

The shared R helper `read_europe_cube()` reads one NetCDF file and returns a Europe-cropped data cube. It infers the variable name from the first underscore-delimited part of the filename, reads `lon`, `lat`, and the variable field, converts longitude to the -180 to 180 convention, sorts longitudes, sorts latitudes if necessary, and crops the data to the European analysis window:

```text
longitude: -15 to 40 degrees
latitude:   34 to 72 degrees
```

The returned cube has the structure:

```text
lon vector
lat vector
v array with dimensions lon x lat x time
variable name
```

For precipitation, `read_europe_cube()` converts the additive anomaly from `kg m-2 s-1` to `mm/day` using:

```text
pr_mm_day = pr_kg_m2_s * 86400
```

This conversion is applied once when R reads the precipitation anomaly for plotting, delta construction, or diagnostics. Temperature variables are left in kelvin. For temperature anomalies, a change of 1 K is numerically the same size as a change of 1 degree Celsius.

Regional means are computed with cosine-of-latitude weights. The helper `masked_mean_series()` flattens the longitude-latitude grid, keeps the requested cell indices, assigns each selected cell a weight proportional to:

```text
cos(latitude)
```

and returns one area-weighted monthly mean anomaly time series. This is used for Europe-wide means, land-only means, and country-level means where applicable.

## C.5 Late-window mean response and delta-field construction

Let the Europe-mean monthly anomaly series be:

```text
a_1, a_2, ..., a_n
```

where `n` is the number of monthly time steps available for a given anomaly file. The late-window length is defined as:

```text
pl = max(12, round(n / 3))
```

The Europe-mean late-window response is then:

```text
late_window_mean = mean(a_(n - pl + 1), ..., a_n)
```

In the code this value is called `plateau` or `plate`, depending on the script. Scientifically, it should be interpreted as a late-window mean response unless an additional trend diagnostic shows that the late window is effectively stable.

The per-cell delta field is computed analogously, but independently at each grid cell. If `A_v(cell, t)` is the monthly additive anomaly field for variable `v`, then the delta field is:

```text
delta_v(cell) = mean(A_v(cell, n - pl + 1 : n))
```

The R implementation flattens the Europe-cropped anomaly cube into a matrix with structure:

```text
rows    = grid cells
columns = monthly time steps
```

It then averages the final `pl` columns for each cell and reshapes the result back to a two-dimensional `lon x lat` field. The result is written to NetCDF as:

```text
delta_<model>_<protocol>_<variable>_lastthird.nc
```

The delta NetCDF contains one variable:

```text
tas, tasmin, tasmax, or pr
```

with units:

```text
K       for tas, tasmin, tasmax
mm/day  for pr
```

Missing values are written as `-9999`.

The late-window delta is therefore a time-mean anomaly field over the final third of the available run. It is not a slope, not a final single-month value, and not an automatically detected equilibrium.

For a 50-year monthly hosing run, `n` is approximately 600 months, so `pl = round(600 / 3) = 200` months, or about 16.7 years. For a 100-year monthly hosing run, `n` is approximately 1200 months, so `pl = 400` months, or about 33.3 years.

The final-third rule is an analysis convention in this workflow. It is consistent with the broader climate-model practice of summarising forced responses using late-period means rather than single final time steps, but the exact final-third window is not prescribed by the cited AMOC-hosing literature. If the response is still evolving during the final third, the late-window mean may understate the eventual end-state response.

## C.6 Late-window trend diagnostic

The delta script also computes a trend within the same late-window period used for the Europe-mean response. It takes the final `pl` monthly Europe-mean anomaly values:

```text
win = last pl values of the Europe-mean monthly anomaly series
```

then fits a linear regression:

```text
win ~ time_in_years
```

The slope is multiplied by 10, producing a trend in units per decade:

```text
plateau_trend = 10 * slope
```

This value is a settledness diagnostic. A small late-window trend supports the interpretation that the late-window mean is close to a stabilized response. A large late-window trend indicates that the run is still evolving, so the late-window mean should be treated as a lower-bound or incomplete-response summary rather than as a fully equilibrated plateau.

## C.7 Effect-size diagnosis against piControl variability

The code computes an effect-size diagnostic by comparing the late-window hosing response with the variability of matched-length means in the piControl simulation.

The control series is obtained by `get_control_series()`. This helper identifies the relevant piControl directory from the anomaly file path, reads all matching piControl files because `MAX_PIC_FILES` is set to `Inf`, computes Europe-mean monthly values using the same Europe crop and cosine-latitude weighting, removes the control series’ own monthly climatology, and removes a linear drift.

The control series is therefore processed as:

```text
raw piControl monthly Europe mean
-> deseasonalised control anomaly
-> linearly de-drifted control anomaly
```

The code then constructs a matched-window null distribution using overlapping running means of length `pl`:

```text
null_means = mean(control_s, ..., control_(s + pl - 1))
```

for every possible starting index `s`.

The matched-window control spread is:

```text
sd_window = sd(null_means)
```

The Europe-mean effect ratio is:

```text
effect_ratio = abs(late_window_mean) / sd_window
```

The shared threshold is:

```text
EFFECT_K = 2
```

The script labels the effect as `YES` when:

```text
effect_ratio >= 2
```

This is an effect-size or signal-to-noise diagnostic. It is not a p-value and should not be interpreted as a formal significance probability.

## C.8 Spatial fingerprint: NW-minus-Mediterranean contrast

The shared helper file defines two regional boxes:

```text
NW box:            longitude -15 to 5,   latitude 50 to 62
Mediterranean box: longitude   0 to 25,  latitude 34 to 45
```

For each delta field, the code computes an area-weighted mean in each box and reports:

```text
NW_minus_Med = mean_NW(delta) - mean_Med(delta)
```

For temperature variables, a negative value means that the north-western European box cooled more, or warmed less, than the Mediterranean box. This is interpreted as an AMOC-like spatial fingerprint because AMOC weakening is expected to produce stronger cooling or reduced warming around the North Atlantic and north-western Europe than around the Mediterranean.

For precipitation, the script still computes and reports the same NW-minus-Mediterranean number, but the code comments explicitly treat the precipitation value as a reference diagnostic rather than the main AMOC temperature-fingerprint logic.

## C.9 Land-only and country-level summaries

The delta script optionally computes a land-only Europe-mean late-window response if the R `maps` package is available. It classifies grid cells by whether their cell centre falls inside a country polygon. It then recomputes the area-weighted monthly anomaly time series using only those land cells and applies the same final-third averaging rule.

The time-series plotting script also optionally computes country-level anomaly series. For each country, it selects cells whose centres fall inside the country polygon, computes the same cosine-latitude weighted monthly anomaly series, smooths the trajectory for plotting, and computes the same late-window mean. If country extraction succeeds, the script writes a `country_plateaus.csv` file containing country-level late-window summaries.

These country and land calculations use cell-centre classification. They are useful for approximate land-impact interpretation, but they may be sensitive to grid resolution, coastlines, islands, and small countries.

## C.10 Seasonal delta and per-cell signal-to-noise fields

The uploaded seasonal script computes DJF and JJA seasonal diagnostics on each model’s native grid after applying the same Europe crop. It does not interpolate all models to a common grid.

The seasonal script defines:

```text
DJF = c(12, 1, 2)
JJA = c(6, 7, 8)
```

The code comments state that the monthly files start in January and that this is a calendar-year variant of DJF. This means DJF uses December, January, and February from the same reshaped model year index. The definition is applied consistently to hosing and control data.

For each monthly anomaly cube, the seasonal script reshapes the data into:

```text
lon x lat x month x year
```

and averages the selected season months to produce a yearly seasonal anomaly cube:

```text
A_season(cell, year)
```

The seasonal late-window length is defined in years rather than months:

```text
plY = max(5, round(number_of_years / 3))
```

The seasonal delta at each cell is:

```text
seasonal_delta(cell) = mean of the last plY seasonal anomaly values
```

For the control comparison, the script reads all matching piControl files, converts them to yearly seasonal means on the same native Europe grid, flattens each seasonal control field into cells by years, binds all control years together, detrends each cell’s control seasonal series, and computes the standard deviation of running means of length `plY`.

The per-cell seasonal signal-to-noise ratio is:

```text
snr(cell) = seasonal_delta(cell) / sd_control_running_means(cell)
```

This is a per-cell effect-size measure, not a p-value. The script treats cells with:

```text
abs(snr) >= 2
```

as robust against local matched-timescale control variability. It also computes an area-weighted robust-cell fraction: the percentage of the European window, weighted by cosine latitude, whose absolute seasonal effect size is at least 2.

The seasonal script writes NetCDF files containing both the seasonal delta field and the corresponding SNR field.

## C.11 Plotting outputs and display-only transformations

The delta plotting script reads the exported `delta_*_lastthird.nc` files. It does not recompute the deltas. For precipitation, it assumes the exported delta is already in `mm/day`, because the conversion from native precipitation units occurred once inside `read_europe_cube()` before the delta NetCDF was written.

The plotted maps use a shared, symmetric color scale per variable. The scale is clipped at the 98th percentile of absolute delta values. This clipping affects only the display. It does not change the saved delta NetCDF files.

For precipitation, the plotting script produces an additional relative diagnostic page:

```text
relative_pr_delta_percent = 100 * delta_pr / piControl_annual_mean_pr
```

where both numerator and denominator are in `mm/day`. The script masks cells where control precipitation is very small. This percentage plot is a diagnostic representation of the additive precipitation delta. It is not the same object as the monthly multiplicative precipitation ratio.

## C.12 Outputs produced by the uploaded scripts

The HadGEM anomaly script writes additive anomaly files for each available model, protocol, and variable:

```text
anomaly_output/LL_anomaly/<variable>_anomaly_<protocol>_minus_piControl_1850-1949.nc
anomaly_output/MM_anomaly/<variable>_anomaly_<protocol>_minus_piControl_1850-1949.nc
```

For precipitation, it also writes dimensionless ratio files:

```text
pr_ratio_<protocol>_over_piControl_1850-1949.nc
```

The EC-Earth3 anomaly script writes additive anomaly files into:

```text
anomaly_output/ECHearth3_anomaly/
```

and, for precipitation, writes one dimensionless ratio file per hosing experiment.

The delta and contrast script reads the 24 expected additive anomaly files:

```text
3 models x 2 protocols x 4 variables = 24 files
```

It skips missing files with a printed notice rather than aborting. For each processed file, it writes one delta NetCDF:

```text
delta_<model>_<protocol>_<variable>_lastthird.nc
```

and also writes a contrast table:

```text
amoc_effect_model_contrast.csv
```

The contrast table includes the model, protocol, variable, approximate run length in years, Europe-mean late-window response, land-only late-window response where available, late-window trend, matched-window control spread, effect ratio, effect label, NW-minus-Mediterranean contrast, and the corresponding delta filename.

The time-series plotting script writes:

```text
delta_timeseries_overlay.pdf
```

and, when country extraction is available, also writes:

```text
country_plateaus.csv
```

The seasonal script writes seasonal NetCDFs and a seasonal diagnostic PDF. The seasonal NetCDFs contain both seasonal delta and seasonal SNR fields for DJF and JJA.

## C.13 Main assumptions, limitations, and risks

The late-window mean is not proof of equilibrium. It summarizes the final third of each hosing run. It should be called a late-window mean response unless the late-window trend diagnostic supports a settled interpretation.

The final-third window is an operational choice in this analysis. It reduces interannual variability relative to a shorter final-period mean, but it can dilute the end-state response if the anomaly is still strengthening or weakening during the final third.

For HadGEM, the piControl baseline selection is based on filename date ranges overlapping 185001-194912. The script does not internally crop months within an overlapping file. If a selected file extends beyond the target window, all months in that selected file enter the climatology.

For EC-Earth3, the piControl climatology is built from the files matched by the piControl filename patterns. The years are therefore enforced by file availability, not by an explicit date-slicing command in the script.

For EC-Earth3 native `tas` and `pr`, the script selects the first matching hosing file with `ls ... | head -1`. This is safe only if the matching pattern returns one intended file per variable and protocol, or if the first match is known to be the complete intended file.

Grid handling differs by model. EC-Earth3 explicitly snaps or remaps the piControl climatology to the hosing grid. The HadGEM script assumes hosing and piControl files are already grid-compatible.

Equal grid size does not prove equal grid geometry. In the EC-Earth3 script, equal grid size triggers `setgrid`, which changes metadata but does not interpolate values.

Precipitation ratios can become very large where climatological precipitation is near zero. The scripts intentionally store the ratio faithfully and unclipped. Any later clipping, masking, or bounding must occur in the ISIMIP-application step, not in the anomaly or delta scripts described here.

The additive precipitation delta and the multiplicative precipitation ratio answer different questions. The additive `pr` delta describes the late-window mean precipitation change in `mm/day`. The ratio describes the multiplicative hosing/control precipitation factor to be applied in a later precipitation-stressing workflow.

The uploaded scripts do not apply either temperature anomalies or precipitation ratios to ISIMIP daily data. They only prepare anomaly files, ratio files, delta fields, plots, and diagnostic summaries.

## C.14 Concise methodological interpretation

The workflow first converts AMOC-hosing simulations into monthly anomalies by subtracting a piControl monthly climatology matched by calendar month. It then summarizes the late-run response by averaging the final third of each monthly anomaly series. For each grid cell, this produces a delta field: an additive temperature change in kelvin for `tas`, `tasmin`, and `tasmax`, and an additive precipitation change in `mm/day` for `pr`.

For precipitation, the workflow also creates a separate dimensionless monthly ratio equal to hosing precipitation divided by piControl climatological precipitation. That ratio is the object intended for later multiplicative precipitation adjustment. It is not the same as the saved precipitation delta field.

The diagnostic statistics compare the Europe-mean late-window response with matched-timescale piControl variability, report a simple NW-minus-Mediterranean spatial contrast, and optionally compute land, country, and seasonal signal-to-noise summaries. These diagnostics support interpretation of the AMOC-hosing response, but they should be described as effect-size and late-window response measures rather than as formal p-values or automatically detected equilibria.
