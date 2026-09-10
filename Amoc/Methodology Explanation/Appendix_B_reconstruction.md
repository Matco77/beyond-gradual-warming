# Appendix B — Reconstruction of Daily Temperature Extremes for EC-Earth3

The EC-Earth3 hosing experiments used in this work archive the monthly mean
near-surface air temperature (`tas`) but not the accompanying monthly means of
the daily maximum and minimum temperatures (`tasmax`, `tasmin`). These two
fields are required downstream and are therefore reconstructed statistically.
The reconstruction is calibrated on the EC-Earth3 pre-industrial control
simulation (piControl), in which all three variables are available, and
follows the stationary mean–extreme relationship approach of Weedon et al.
(2010). The procedure operates entirely on monthly data and on the native
model grid, and produces, for each hosing experiment, two files of absolute
reconstructed fields; the hosing-minus-piControl anomalies used in the
analysis are derived from these in a separate, subsequent step, so that a
single processing stage owns the anomaly definition.

## B.1 Method

The method rests on a locally calibrated linear relationship between the
monthly mean temperature and the monthly means of its daily extremes. For
every grid cell and every calendar month $m$ separately, two ordinary
least-squares (OLS) regressions are fitted across the piControl years:

$$\mathrm{tasmax} = a_{\max}(c, m) + b_{\max}(c, m)\,\mathrm{tas},$$
$$\mathrm{tasmin} = a_{\min}(c, m) + b_{\min}(c, m)\,\mathrm{tas},$$

where $c$ indexes the grid cell. The coefficients are obtained in closed form
from the sample moments of that cell and calendar month,

$$b = \frac{\mathrm{Cov}(\mathrm{tas}, \mathrm{extreme})}{\mathrm{Var}(\mathrm{tas})},
\qquad a = \overline{\mathrm{extreme}} - b\,\overline{\mathrm{tas}},$$

with averages taken over all piControl years sharing calendar month $m$. The
slope $b$ expresses how the diurnal extremes respond when the monthly mean of
that cell and season shifts; the intercept $a$ anchors the absolute level, so
that the fitted line passes through the piControl climatological state of the
cell. The reconstruction then consists of evaluating the fitted lines at the
hosing monthly means. Because the fit is stratified by calendar month — all
Januaries pooled together, and so on — only the seasonal phase is matched
between the two simulations: the piControl calibration period (model years
2259–2759; the tasmax archive ends at 2757) and the hosing period (1850–1949)
need not, and do not, overlap in absolute time.

**Input data and pairing.** The piControl inputs are the monthly (`Amon`)
fields of a single ensemble member (`r1i1p1f1`) on the native regular grid
(`gr`); the member and grid identifiers are pinned in the file-selection
patterns so that a second ensemble member appearing in the same directory
could not be silently blended into the calibration. The yearly file chunks of
each variable are concatenated along the time axis by coordinate value, with
time decoded through the CF calendar conventions so that model years far
outside the conventional datetime range are represented correctly. Two
validity checks are applied at load time: the expected variable must be
present in the files by name — a differently named variable is never
substituted — and its units must be kelvin. The three piControl series are
then reduced to their common time steps by an inner join; the numbers of time
steps before and after the join are reported, and the procedure aborts if the
join removes grid cells rather than time steps, which would indicate
inconsistent spatial coordinates among the input variables. Only the paired
sample enters the regression, so all three variables are fitted on identical
years.

## B.2 Numerical treatment

**Well-conditioned moment computation.** The covariance and variance are not
computed with the algebraically equivalent shortcut
$\mathrm{E}[xy] - \mathrm{E}[x]\mathrm{E}[y]$. The input fields are stored in
single precision at absolute temperatures of roughly 300 K, so
$\mathrm{E}[x^2]$ is of order $9\times10^4$ while the interannual variance of
a monthly mean can lie well below 1 K²; subtracting two large, nearly equal
numbers in single precision then leaves mostly rounding error (catastrophic
cancellation), degrading the fitted slopes precisely where natural
variability is smallest. The implementation therefore casts the data to
double precision and centres them within each calendar month before forming
any products,

$$\mathrm{Cov} = \overline{(x-\bar{x})(y-\bar{y})}, \qquad
\mathrm{Var} = \overline{(x-\bar{x})^2},$$

which is the same estimator in exact arithmetic but numerically stable.

**Degenerate-variance fallback.** In cells whose piControl monthly mean is
effectively constant across years, the slope is undefined. The criterion is
$\mathrm{Var}(\mathrm{tas}) \le 10^{-6}\,\mathrm{K}^2$, an engineering
threshold corresponding to a standard deviation below one millikelvin over
the full calibration period. Such cells receive $b = 1$ and
$a = \overline{\mathrm{extreme}} - \overline{\mathrm{tas}}$, so that the
reconstructed extreme tracks the hosing mean displaced by the piControl
climatological offset. The alternative of setting $b = 0$ would hold the cell
at the piControl climatology regardless of the hosing signal, which is
untenable where that signal reaches several kelvin. The number of cells
handled by this fallback is reported at run time.

**Coefficient evaluation and diagnostics.** The four coefficient maps
(12 calendar months × grid, for slope and intercept of each extreme) are
evaluated once and held in memory, so that the full multi-century calibration
archive is read a single time rather than once per output file. Summary
statistics of the fitted slopes — minimum, median, maximum, the count of
cells with $|b| > 3$, and the count of fallback cells — are printed for every
fit, providing an immediate plausibility check on the calibration.

**Grid consistency between calibration and application.** Before the
coefficients are applied, the piControl and hosing grids are compared: they
are considered identical only if both latitude and longitude exist in the two
datasets with equal lengths and numerically indistinguishable values. When
the grids are identical — the case for the runs processed here — the
coefficient maps are relabelled with the hosing coordinate values exactly.
This step exists because coordinates that are equal in substance can differ
in their stored floating-point representation, and label-exact agreement is
required for the subsequent cell-by-cell arithmetic; the relabelling changes
no data values. If the grids genuinely differed, the coefficient maps (not
the data) would instead be transferred to the hosing grid by
nearest-neighbour lookup, with edge cells outside the source range taking the
nearest boundary value rather than becoming undefined.

**Application to the hosing simulation.** Each hosing time step is assigned
the coefficient pair of its calendar month, and the reconstruction is
evaluated as $a(c, m) + b(c, m)\,\mathrm{tas}(c, t)$ for every cell and time
step, grouping the hosing series by calendar month so that the coefficient
maps are never expanded to the length of the time axis. The predictions are
cast back to single precision, matching the precision of the inputs.

**Quality control and ordering constraint.** Two checks follow. First, a NaN
audit: each prediction inherits, at most, the missing values already present
in the hosing input, so the procedure aborts if the two reconstructed fields
together contain more NaN values than twice the count in the hosing `tas` —
the exact number expected when no invalid values have been introduced by the
regression or the grid handling. Second, the physical ordering
$\mathrm{tasmax} \ge \mathrm{tas} \ge \mathrm{tasmin}$ is enforced pointwise:
a reconstructed maximum falling below the mean is raised to the mean, and a
reconstructed minimum above it is lowered to the mean. This truncation is
one-sided, and its incidence is therefore both printed per experiment and
stored in each output file as a `clamp_fraction` attribute, making the
extent of the correction auditable after the fact.

**Output.** For each hosing experiment, the reconstructed `tasmax` and
`tasmin` are written as compressed single-precision NetCDF files on the
hosing time axis and grid, with CF-standard variable metadata, a global
attribute stating the reconstruction method and its stationarity assumption,
and a provenance attribute recording the calibration source, ensemble
member, grid, period, and sample size. A missing hosing experiment is
skipped with a message rather than aborting the run, so partial data
availability does not block the remaining experiment.

## B.3 Validation

The implementation includes a built-in synthetic self-check that runs
without any input data and exercises the full fit-and-apply chain. It
constructs a century of monthly values in single precision with a small
interannual spread (0.1 K) around 280 K — deliberately the regime in which a
naively computed variance collapses numerically — from prescribed
month-dependent slopes and intercepts, and asserts that (i) the fitted
slopes recover the prescribed values to better than $10^{-2}$; (ii) a cell
with strictly constant input receives exactly the fallback slope $b = 1$;
and (iii) applying the fitted coefficients to a shifted series reproduces
$a + b\,\mathrm{tas}$ to machine precision. In addition, every production
run reports its own diagnostics, as described above: pairing counts of the
calibration sample, slope summary statistics, fallback counts, NaN
accounting, and clamp incidence per experiment.

## B.4 Assumptions and limitations

The central assumption, stated in the output metadata itself, is
stationarity: the mean–extreme relationship estimated from the internal
variability of the control climate is taken to remain valid under the hosing
perturbation. Related to this, the fitted lines are evaluated at whatever
monthly means the hosing simulation produces; no restriction confines the
application to the range of means sampled during calibration, so the linear
form is trusted in extrapolation where the hosing signal exceeds piControl
variability. The ordering constraint corrects fitted values that violate
$\mathrm{tasmax} \ge \mathrm{tas} \ge \mathrm{tasmin}$ by truncation toward
the mean, a one-sided adjustment whose frequency is quantified by the
recorded clamp fractions. Finally, the reconstruction is monthly, inheriting
the temporal resolution of the input `tas`: it estimates the monthly mean of
the daily extremes, not daily values, and captures the response of these
extremes to shifts in the monthly mean rather than any change in the
sub-monthly temperature distribution itself.

---

**References.**
Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological
forcing dataset for land surface- and hydrological-models. WATCH Technical
Report.
