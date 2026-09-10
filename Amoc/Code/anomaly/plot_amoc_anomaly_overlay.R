## AMOC anomaly overlay: SSP projections (ΔSv, calendar years) vs U03 hosing
## (ΔSv, years since hosing start) on a shared y-axis, one panel per model.
## Same base-R + ncdf4 style as plot_amoc_sv_ssp126.R / plot_amoc_sv_hosing.R.
## Each series anomalised against its OWN source reference (never cross datasets).
## Raw annual values, no smoothing.
##
## TWO emissions scenarios on the bottom calendar axis, hosing stays dashed on the top axis.
## Scenario is carried by SHADE, not by line type: both projections are solid (SSP1-2.6 in the
## model colour, SSP3-7.0 in a darker shade of it), because a dotdash SSP3-7.0 was not reliably
## distinguishable from the dashed hosing line at this line width - checked on the rendered PDF,
## not assumed. Dashed therefore means "hosing, top axis" and nothing else.
##
## SOURCES - one pipeline per model, covering the historical reference and BOTH scenarios, so a
## SSP1-2.6-vs-SSP3-7.0 difference cannot partly measure a change of diagnostic:
##   EC-Earth3 : own vo below-500m reconstruction (amoc_from_vo_below500m.py).
##   IPSL      : raw msftyz (amoc_ipsl_from_msftyz.py). Terhaar's published IPSL series is NOT used:
##               he has no ssp370, and his ssp126 is a different run (extends to 2214, 2071-2100 mean
##               10.07 Sv vs 9.43 here). The msftyz pipeline reproduces his HISTORICAL at corr 0.9971
##               / debiased rmse 0.081 Sv, which is the gate that licenses the substitution.
## The two MODELS use different diagnostics - EC-Earth3 publishes no msftyz for the scenarios,
## checked on ESGF - but each model's anomaly is internally consistent, which is what matters here.
library(ncdf4)

ds  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
out <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly/amoc_anomaly_overlay.pdf"

cols <- c("EC-Earth3" = "#1b9e77", "IPSL-CM6A-LR" = "#e7298a")

## ---- readers (same conventions as the other two scripts) ----
read_terhaar <- function(f) {                     # amoc(year, x), Sv, calendar years
  nc <- nc_open(f)
  d  <- list(t = as.numeric(ncvar_get(nc, "year")), a = as.numeric(ncvar_get(nc, "amoc")))
  nc_close(nc); d
}
read_ecearth <- function(f) {                     # amoc(time), row already selected; time = year+0.5
  nc <- nc_open(f)
  d  <- list(t = floor(as.numeric(ncvar_get(nc, "time"))),
             a = as.numeric(ncvar_get(nc, "amoc")))
  nc_close(nc); d
}
read_hos <- function(m) {                         # M26: hos/con(time), time = years since hosing start
  nc <- nc_open(file.path(ds, "M26", paste0("M26_", m, ".nc")))
  d  <- list(t = as.numeric(ncvar_get(nc, "time")),
             hos = as.numeric(ncvar_get(nc, paste0("hos_", m))),
             con = as.numeric(ncvar_get(nc, paste0("con_", m))))
  nc_close(nc)
  d$hos[abs(d$hos) > 100] <- NA                   # unflagged fill in tail
  d$con[abs(d$con) > 100] <- NA
  d
}

## ---- projection anomalies (recompute 1850-1900 reference from each own historical) ----
ref_mean <- function(d) mean(d$a[d$t >= 1850 & d$t <= 1900], na.rm = TRUE)

cut2100 <- function(d) { k <- d$t <= 2100; list(t = d$t[k], a = d$a[k]) }

# per model: one reader, the historical reference file, and the file for each scenario. Same
# AMOC_SRC layout as isimip_bin_fields.R, so the figure and the binning cannot drift apart.
SRC <- list(
  "EC-Earth3" = list(
    read = read_ecearth, hist = "ecearth3_amoc26N_vo_below500m_historical_1850_2014.nc",
    scen = c(ssp126 = "ecearth3_amoc26N_vo_below500m_ssp126_2015_2100.nc",
             ssp370 = "ecearth3_amoc26N_vo_below500m_ssp370_2015_2100.nc")),
  "IPSL-CM6A-LR" = list(
    read = read_terhaar,                        # msftyz files carry the same (year, amoc) schema
    hist = "amoc_ipsl_msftyz_26N_historical_1850_2014.nc",
    scen = c(ssp126 = "amoc_ipsl_msftyz_26N_ssp126_2015_2100.nc",
             ssp370 = "amoc_ipsl_msftyz_26N_ssp370_2015_2100.nc")))
REF  <- sapply(names(SRC), function(m) ref_mean(SRC[[m]]$read(file.path(ds, SRC[[m]]$hist))))
LAB  <- c(ssp126 = "SSP1-2.6", ssp370 = "SSP3-7.0")
SCEN <- names(LAB)                              # SRC is keyed by MODEL, so scenarios come from LAB
# a darker version of a colour: a*colour, i.e. mixed toward black. Opaque, no alpha, so it needs
# no PDF transparency support (same reason plot_impact_common.R::tint() avoids alpha).
shade    <- function(col, a) { x <- col2rgb(col) / 255; rgb(x[1]*a, x[2]*a, x[3]*a) }
scen_col <- function(col, s) if (s == "ssp370") shade(col, 0.55) else col

proj <- lapply(SCEN, function(s) {
  z <- lapply(names(cols), function(m) {
    S <- SRC[[m]]
    d <- cut2100(S$read(file.path(ds, S$scen[[s]])))
    list(t = d$t, d = d$a - REF[[m]])
  }); names(z) <- names(cols); z
}); names(proj) <- SCEN

## ---- hosing anomalies: hos - mean(con) of same model ----
hos <- lapply(names(cols), function(m) {
  h <- read_hos(m); list(t = h$t, d = h$hos - mean(h$con, na.rm = TRUE))
}); names(hos) <- names(cols)

## shared y across everything so the two panels are directly comparable
ylim <- range(unlist(lapply(names(cols), function(m)
  c(unlist(lapply(SCEN, function(s) proj[[s]][[m]]$d)), hos[[m]]$d))), na.rm = TRUE)

## ---- plot ----
pdf(out, width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 5, 3.2), oma = c(0, 1.5, 2.5, 0))
for (m in names(cols)) {
  col <- cols[m]
  ## bottom axis: projections, calendar years 2015-2100 (ssp126 solid, ssp370 dotdash)
  plot(NA, xlim = c(2015, 2100), ylim = ylim, axes = FALSE,
       xlab = "Year (SSP projections, bottom)", ylab = "")
  axis(1); axis(2); box()
  abline(h = 0, col = "grey70")
  for (s in SCEN) lines(proj[[s]][[m]]$t, proj[[s]][[m]]$d, col = scen_col(col, s), lwd = 2)
  title(main = m, line = 3.3)
  ## top axis: hosing, years since hosing start 0-150 (dashed)
  par(new = TRUE)
  plot(hos[[m]]$t, hos[[m]]$d, type = "l", col = col, lwd = 1.5, lty = 2,
       xlim = c(0, 150), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  axis(3)
  mtext("Years since hosing start (U03, top)", side = 3, line = 1.6, cex = 0.75, col = "grey30")
  if (m == "EC-Earth3")
    legend("bottomleft", bg = "white", box.col = NA, cex = 0.85,   # opaque: the hosing line runs behind it
           legend = c(sprintf("%s (bottom axis)", LAB[SCEN]), "U03 hosing (dashed, top axis)"),
           col = c(sapply(SCEN, function(s) scen_col(col, s)), col),
           lty = c(rep(1, length(SCEN)), 2), lwd = c(rep(2, length(SCEN)), 1.5))
}
mtext("AMOC anomaly: SSP1-2.6 / SSP3-7.0 projections vs U03 hosing", outer = TRUE, cex = 1.15, font = 2, line = 0.8)
mtext("AMOC anomaly (Sv, vs preindustrial reference)", side = 2, outer = TRUE, line = 0.1)
dev.off()
cat("wrote", out, "\n\n")

## ---- diagnostics, per model x scenario ----
for (m in names(cols)) {
  ho <- range(hos[[m]]$d, na.rm = TRUE)
  cat(sprintf("%-13s ΔSv_hos [% .2f, % .2f]\n", m, ho[1], ho[2]))
  for (s in SCEN) {
    pr <- range(proj[[s]][[m]]$d, na.rm = TRUE)
    inside <- ho[1] <= pr[1] && ho[2] >= pr[2]
    target <- pr[1]                                      # most negative projected anomaly
    hit <- which(hos[[m]]$d <= target)                   # first hosing year reaching it
    hit_yr <- if (length(hit)) hos[[m]]$t[hit[1]] else NA
    cat(sprintf("%-13s   %s ΔSv_proj [% .2f, % .2f] | inside hosing range: %-5s | hosing reaches %.2f Sv at year %s\n",
                "", LAB[[s]], pr[1], pr[2], inside, target,
                ifelse(is.na(hit_yr), "never", as.character(hit_yr))))
  }
  cat("\n")
}
