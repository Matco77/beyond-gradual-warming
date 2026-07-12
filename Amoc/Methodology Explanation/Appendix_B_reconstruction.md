# Appendix B — Reconstruction of Daily Temperature Extremes for EC-Earth3

In the EC-Earth3 hosing experiments only the monthly mean near-surface air
temperature (`tas`) and precipitation (`pr`) are available; the daily minimum
and maximum temperatures (`tasmin`, `tasmax`) are not archived. These two
fields are reconstructed from the monthly mean using a statistical relationship
calibrated on the piControl simulation, following the approach of Weedon et al.
(2010). The HadGEM3 configurations require no such step, as they provide all
four variables natively.

## B.1 Method

The reconstruction exploits the fact that, within a stable climate, the monthly
mean temperature and its diurnal extremes are linearly related, and that this
relationship varies smoothly in space and across the seasonal cycle. For every
grid cell and every calendar month *m*, an ordinary least-squares (OLS)
regression of each extreme on the mean is fitted across the piControl years:

    tasmax = a_max(cell, m) + b_max(cell, m) · tas
    tasmin = a_min(cell, m) + b_min(cell, m) · tas

The slope and intercept are obtained in closed form from the piControl
statistics of that cell and month,

    b = Cov(tas, extreme) / Var(tas),     a = E[extreme] − b · E[tas],

where the expectations are taken over all piControl years sharing calendar
month *m*. The fitted coefficients are then applied to the hosing monthly mean
to obtain the reconstructed extremes.

The slope *b* measures how the diurnal temperature range responds as the
monthly mean shifts, while the intercept *a* anchors the absolute level so that
the climatological seasonal cycle is preserved. Fitting separately for each
calendar month removes the seasonal cycle from the regression and, importantly,
makes the calibration independent of the absolute model years: the piControl
segment (model years 2259–2759) and the hosing segment (1850–1949) need not
overlap in time, since only the seasonal phase is matched.

## B.2 Numerical treatment

Two safeguards are applied. In grid cells where the monthly mean does not vary
across the piControl years (Var(tas) = 0; e.g. permanently masked cells), the
slope is undefined; there the slope is set to zero, so that the extreme follows
the mean by its fixed climatological offset. After reconstruction the physical
ordering tasmax ≥ tas ≥ tasmin is enforced, pinning any over-shooting fitted
value to the mean. Both safeguards affect only degenerate cells and have no
effect over Europe.

The reconstruction is performed on the native EC-Earth3 grid (`gr`,
512 × 256). The piControl and hosing grids are identical to within
floating-point precision, so no spatial interpolation is involved.

## B.3 Validation

The implementation was verified against a synthetic case with a known linear
relationship, which it recovered exactly, and the zero-variance safeguard was
confirmed to return a zero slope without numerical error. The reconstructed
fields lie within a physically plausible range (e.g. reconstructed `tasmax`
spans approximately 210–323 K globally, with no missing values).

## B.4 Assumptions and limitations

The method assumes that the mean–extreme relationship calibrated on piControl
remains valid under AMOC weakening (a stationarity assumption). It therefore
captures the response of the extremes to a shift in the monthly mean, but not
any change in the *shape* of the sub-monthly temperature distribution that the
hosing forcing might induce; in this sense it is a first-order statistical
reconstruction rather than a physical simulation of extremes. The reconstructed
fields are monthly, consistent with the resolution of the input `tas`.

---

**References.**
Weedon, G. P., et al. (2010). The WATCH Forcing Data: a meteorological forcing
dataset for land surface- and hydrological-models. WATCH Technical Report.
