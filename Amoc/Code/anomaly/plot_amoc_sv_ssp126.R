## AMOC at 26N (Sv): SSP1-2.6 decline vs historical, IPSL (own msftyz pipeline) + EC-Earth3 (own
## reconstruction from vo). Same style as plot_amoc_sv_hosing.R; raw annual values, no smoothing.
## Single continuous 1850-2100 axis: dashed = historical, solid = ssp126.
##
## This is the trajectory isimip_bin_fields.R bins, so it MUST read the same files that script does.
## Terhaar's published IPSL series is not used: he has no ssp370, and his ssp126 is a different run
## (extends to 2214, 2071-2100 mean 10.07 Sv against 9.43 Sv here), so keeping it would have made
## IPSL's two scenarios different diagnostics. The msftyz pipeline reproduces his HISTORICAL at
## corr 0.9971 / debiased rmse 0.081 Sv - the gate that licenses the substitution.
library(ncdf4)

ds  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
out <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly/amoc_sv_ssp126.pdf"

models <- c("EC-Earth3", "IPSL-CM6A-LR")
cols   <- c("#1b9e77", "#d95f02")   # match original green/orange

## IPSL msftyz product: amoc(year), calendar years, already Sv. Same (year, amoc) schema Terhaar used.
read_terhaar <- function(f) {
  nc <- nc_open(f)
  d  <- list(t = as.numeric(ncvar_get(nc, "year")),
             a = as.numeric(ncvar_get(nc, "amoc")))
  nc_close(nc); d
}
## Recomputed EC-Earth3 vo (max below 500 m): amoc(time), row already selected; time = year+0.5.
read_ecearth <- function(f) {
  nc <- nc_open(f)
  d  <- list(t = floor(as.numeric(ncvar_get(nc, "time"))),
             a = as.numeric(ncvar_get(nc, "amoc")))
  nc_close(nc); d
}

dat <- list(
  "EC-Earth3" = list(
    hist = read_ecearth(file.path(ds, "ecearth3_amoc26N_vo_below500m_historical_1850_2014.nc")),
    ssp  = read_ecearth(file.path(ds, "ecearth3_amoc26N_vo_below500m_ssp126_2015_2100.nc"))),
  "IPSL-CM6A-LR" = list(
    hist = read_terhaar(file.path(ds, "amoc_ipsl_msftyz_26N_historical_1850_2014.nc")),
    ssp  = read_terhaar(file.path(ds, "amoc_ipsl_msftyz_26N_ssp126_2015_2100.nc")))
)
## both ssp126 series already end in 2100; the cut is kept as a guard, not a fix.
for (m in models) { s <- dat[[m]]$ssp; keep <- s$t <= 2100; dat[[m]]$ssp <- list(t = s$t[keep], a = s$a[keep]) }

xlim <- c(1850, 2100)
ylim <- range(unlist(lapply(dat, function(d) c(d$hist$a, d$ssp$a))), na.rm = TRUE)

pdf(out, width = 9, height = 6)
par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = xlim, ylim = ylim, xlab = "Year",
     ylab = "AMOC strength at 26N (Sv)", main = "")
title(main = "AMOC decline under SSP1-2.6 vs historical", line = 2)
mtext("CMIP6 projections - own reconstructions: msftyz (IPSL-CM6A-LR), vo below 500 m (EC-Earth3)",
      side = 3, line = 0.5, cex = 0.85, col = "grey30")
grid(col = "grey90")
abline(v = 2015, col = "grey70", lty = 3, lwd = 1)   # hist / ssp126 split
for (i in seq_along(models)) {
  d <- dat[[models[i]]]
  lines(d$hist$t, d$hist$a, col = cols[i], lty = 2, lwd = 1)   # historical
  lines(d$ssp$t,  d$ssp$a,  col = cols[i], lty = 1, lwd = 2)   # ssp126
}
legend("bottomleft", bty = "n",
       legend = c(models, "ssp126 (solid)", "historical (dashed)"),
       col = c(cols, "black", "black"),
       lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1))
dev.off()
cat("wrote", out, "\n")

## Report: ssp126 minimum + 1850-1900 historical mean per model.
for (m in models) {
  d <- dat[[m]]
  hb <- d$hist$a[d$hist$t >= 1850 & d$hist$t <= 1900]
  cat(sprintf("%-13s ssp126 min = %5.2f Sv (year %d) | 1850-1900 hist mean = %5.2f Sv\n",
              m, min(d$ssp$a, na.rm = TRUE), d$ssp$t[which.min(d$ssp$a)],
              mean(hb, na.rm = TRUE)))
}
