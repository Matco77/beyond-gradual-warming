# Appendix B - Reconstruction of Monthly Means of Daily Temperature Extremes for EC-Earth3

The EC-Earth3 hosing experiments used in this work archive the monthly mean
near-surface air temperature (`tas`) but not the accompanying monthly means of
the daily maximum and minimum near-surface air temperatures (`tasmax`,
`tasmin`). These two fields are required downstream and are therefore
reconstructed statistically. The reconstruction is calibrated on the EC-Earth3
pre-industrial control simulation (piControl), in which all three variables are
available. It is implemented here as a local empirical transfer from monthly
mean temperature to the monthly means of daily extremes: for each grid cell and
calendar month, an ordinary least-squares relationship is fitted between
piControl `tas` and the corresponding piControl `tasmax` or `tasmin`, and the
resulting coefficients are then applied to hosing `tas`.

This should be interpreted as a local statistical reconstruction inspired by,
or in the spirit of, WATCH/WFD-style temperature adjustment and ISIMIP-style
bias-adjustment methods, which treat mean temperature, temperature range, and
temperature extremes as statistically linked quantities. It is not a direct
implementation of the WATCH Forcing Data equations. In particular, Weedon et
al. describe adjustment of 2 m temperature using monthly mean temperature and
diurnal temperature range (DTR), whereas the present script fits an explicit
per-cell, per-calendar-month OLS relationship between model `tas` and model
`tasmax`/`tasmin`. The procedure operates entirely on monthly CMIP `Amon` data
and produces, for each hosing experiment, two files of absolute reconstructed
fields. The hosing-minus-piControl anomalies used in the analysis are derived
from these files in a separate CDO step, so that a single processing stage owns
the anomaly definition.

## B.1 Method

The method rests on a locally calibrated linear relationship between monthly
mean temperature and the monthly means of daily maximum and minimum
temperature. For every grid cell and every calendar month `m` separately, two
ordinary least-squares (OLS) regressions are fitted across the piControl time
series:

$$\mathrm{tasmax} = a_{\max}(c, m) + b_{\max}(c, m)\ \mathrm{tas},$$
$$\mathrm{tasmin} = a_{\min}(c, m) + b_{\min}(c, m)\ \mathrm{tas},$$

where `c` indexes the grid cell. The coefficients are obtained in closed form
from the sample moments of that cell and calendar month,

$$b = \frac{\mathrm{Cov}(\mathrm{tas}, \mathrm{extreme})}{\mathrm{Var}(\mathrm{tas})},
\qquad a = \overline{\mathrm{extreme}} - b\ \overline{\mathrm{tas}},$$

with averages taken over all piControl years sharing calendar month `m`. The
slope `b` expresses how the monthly mean of the daily maximum or minimum
temperature co-varies with the monthly mean temperature in that cell and
season. The intercept `a` anchors the absolute level, so that the fitted line
passes through the piControl climatological state of the cell. The
reconstruction then consists of evaluating the fitted lines at the hosing
monthly means.

Because the fit is stratified by calendar month - all Januaries pooled
together, all Februaries pooled together, and so on - only the seasonal phase
is matched between the two simulations. The piControl and hosing years do not
need to overlap in absolute time. The code reads the piControl and hosing
periods from the files selected by the file patterns; in the EC-Earth3 archive
used here, the available piControl calibration files are the long piControl
period used as the unforced reference, and the hosing files retain their own
experiment-specific time axes. The output metadata records the actual first
and last piControl years and the number of calibration months loaded from the
input files.

**Input data and pairing.** The calibration uses EC-Earth3 piControl monthly (`Amon`) files for `tas`, `tasmax`, and `tasmin`. The file patterns are restricted to the ensemble member `r1i1p1f1` and the native `gr` grid, so files from another ensemble member or grid would not accidentally be included. For each variable, all matching files are opened together and combined along their coordinate-defined time axis. The time coordinate is decoded using CF calendar conventions, which allows model-calendar years to be handled correctly. At loading, the code checks that the expected variable is actually present in the NetCDF file and that its units are kelvin. The three piControl variables are then aligned to their shared coordinates using an inner join. This keeps only the months that are present in all three variables. The code reports the number of time steps before and after this alignment and aborts if the alignment removes spatial grid cells, because that would indicate inconsistent latitude or longitude coordinates. Therefore, the regression is fitted only on paired piControl months where `tas`, `tasmax`, and `tasmin` are all available on the same grid.


## B.2 Numerical treatment

**Well-conditioned moment computation.** The covariance and variance are not
computed with the algebraically equivalent shortcut
`E[xy] - E[x]E[y]`. The input fields are stored in single precision at
absolute temperatures of roughly 300 K, so `E[x^2]` is of order `9 x 10^4`
while the interannual variance of a monthly mean can lie well below 1 K^2;
subtracting two large, nearly equal numbers in single precision then leaves
mostly rounding error, degrading the fitted slopes precisely where natural
variability is smallest. The implementation therefore casts the data to double
precision and centres them within each calendar month before forming any
products,

$$\mathrm{Cov} = \overline{(x-\bar{x})(y-\bar{y})}, \qquad
\mathrm{Var} = \overline{(x-\bar{x})^2},$$

which is the same estimator in exact arithmetic but numerically more stable.
Converts the data to double precision and subtracts the calendar-month mean before computing covariance and variance. 
This preserves the small interannual variability that determines the regression slope.

**Degenerate-variance fallback.** In cells whose piControl monthly mean is
effectively constant across years, the slope is undefined. The criterion is
`Var(tas) <= 10^-6 K^2`, an engineering threshold corresponding to a standard
deviation below one millikelvin over the calibration period. In those cases with var(tas) close to 0 the code uses a fallback slope of b=1. This means the reconstructed tasmax or tasmin follows the hosing monthly mean temperature one-for-one, while preserving the usual piControl offset between the extreme and the mean. For example, if piControl tasmax is typically 5 K warmer than tas, the fallback reconstructs hosing tasmax as hosing tas + 5 K. This is preferable to b=0, which would keep the reconstructed extreme fixed at the piControl climatology and would ignore the hosing temperature signal.

**Grid consistency between calibration and application.** The regression coefficients are learned from the piControl files, but they must be applied to the hosing `tas` file. Therefore, before applying the regression, the code checks whether the piControl and hosing latitude/longitude grids describe the same cells. If the grids are effectively identical, the coefficient values are left unchanged, but their latitude and longitude labels are replaced with the exact labels from the hosing file. This avoids technical alignment problems caused by tiny floating-point differences in coordinate labels. If the grids are genuinely different, the code does not modify the hosing `tas` data. Instead, it transfers the fitted coefficient maps to the hosing grid using nearest-neighbour interpolation.

**Application to the hosing simulation.** For each hosing month, the code uses the coefficient maps for the corresponding calendar month. For example, all hosing Januaries use the January coefficients, all hosing Februaries use the February coefficients, and so on. At each grid cell and time step, the reconstructed value is calculated as (a + b \cdot tas). This produces one reconstructed `tasmax` field and one reconstructed `tasmin` field on the hosing time axis and grid. The calculation is completed in memory before writing, and the predictions are stored as single-precision output.

**Quality control and ordering constraint.** After reconstruction, the code checks whether the predictions introduced missing values beyond those already present in the hosing `tas`. Since two output fields are produced, the combined number of missing values in reconstructed `tasmax` and `tasmin` should not exceed twice the number of missing values in hosing `tas`. The code also enforces the physical ordering (tasmax \ge tas \ge tasmin). If a reconstructed maximum is lower than the mean temperature, it is raised to the mean. If a reconstructed minimum is higher than the mean, it is lowered to the mean. The fraction of corrected points is printed and saved as a `clamp_fraction` attribute.

**Output.** For each available hosing experiment, the code writes two compressed NetCDF files: one for reconstructed `tasmax` and one for reconstructed `tasmin`. These files contain absolute monthly reconstructed fields, not anomalies. They use the hosing time axis and grid and include metadata describing the OLS reconstruction method, the stationarity assumption, and the piControl source used for calibration.

## B.3 Validation

The implementation includes a built-in synthetic self-check that runs without
any input data and exercises the fit-and-apply chain. It constructs 1200
monthly values, corresponding to 100 years, in single precision with a small
interannual spread around 280 K - deliberately the regime in which a naively
computed variance can collapse numerically - from prescribed month-dependent
slopes and intercepts. The test asserts that: (i) the fitted slopes recover
the prescribed values to better than `10^-2` in non-constant cells; (ii) a
cell with strictly constant input receives exactly the fallback slope `b = 1`;
and (iii) applying the fitted coefficients to a shifted series reproduces
`a + b tas` to machine precision. In addition, every production run reports
its own diagnostics, as described above: pairing counts of the calibration
sample, slope summary statistics, fallback counts, NaN accounting, and clamp
incidence per experiment.

## B.4 Assumptions and limitations

The central assumption is stationarity: the mean-extreme relationship
estimated from the internal variability of the control climate is taken to
remain valid under the hosing perturbation. Related to this, the fitted lines
are evaluated at whatever monthly means the hosing simulation produces; no
restriction confines the application to the range of means sampled during
calibration, so the linear form is trusted in extrapolation where the hosing
signal exceeds piControl variability. The ordering constraint corrects fitted
values that violate `tasmax >= tas >= tasmin` by truncation toward the mean, a
one-sided adjustment whose frequency is quantified by the recorded clamp
fractions.

Finally, the reconstruction is monthly, inheriting the temporal resolution of
the input `tas`: it estimates the monthly mean of the daily extremes, not a
new daily time series, and captures the response of these monthly extreme
fields to shifts in the monthly mean rather than any change in the sub-monthly
temperature distribution itself. Although the script metadata string names
Weedon et al. (2010), the methodological interpretation used here is more
cautious: Weedon/WFD provides background for local monthly temperature and DTR
adjustment, while the exact per-cell, per-calendar-month OLS mean-to-extreme
transfer is the empirical reconstruction implemented in this workflow.

---

**References.**

Weedon, G. P., Gomes, S., Viterbo, P., Shuttleworth, W. J., Blyth, E., Osterle,
H., Adam, J. C., Bellouin, N., Boucher, O., and Best, M. (2011). Creation of
the WATCH Forcing Data and its use to assess global and regional reference
crop evaporation over land during the twentieth century. *Journal of
Hydrometeorology*, 12, 823-848. https://doi.org/10.1175/2011JHM1369.1

Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological forcing
dataset for land surface- and hydrological-models. WATCH Technical Report No.
22. Relevant background section: 2b, 2 m temperature.

Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling
with ISIMIP3BASD (v1.0). *Geoscientific Model Development*, 12, 3055-3070.
https://doi.org/10.5194/gmd-12-3055-2019

