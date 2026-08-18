## AMOC anomaly overlay: SSP1-2.6 projection (ΔSv, calendar years) vs U03 hosing
## (ΔSv, years since hosing start) on a shared y-axis, one panel per model.
## Same base-R + ncdf4 style as plot_amoc_sv_ssp126.R / plot_amoc_sv_hosing.R.
## Each series anomalised against its OWN source reference (never cross datasets).
## Raw annual values, no smoothing.
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

ec_hist  <- read_ecearth(file.path(ds, "ecearth3_amoc26N_vo_below500m_historical_1850_2014.nc"))
ec_ssp   <- read_ecearth(file.path(ds, "ecearth3_amoc26N_vo_below500m_ssp126_2015_2100.nc"))
ip_hist  <- read_terhaar(file.path(ds, "terhaar_amoc/amoc/amoc/26.5N/historical/amoc_historical_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc"))
ip_ssp   <- read_terhaar(file.path(ds, "terhaar_amoc/amoc/amoc/26.5N/ssp126/amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc"))
keep <- ip_ssp$t <= 2100; ip_ssp <- list(t = ip_ssp$t[keep], a = ip_ssp$a[keep])   # cut IPSL 2214 -> 2100

proj <- list(
  "EC-Earth3"    = list(t = ec_ssp$t, d = ec_ssp$a - ref_mean(ec_hist)),
  "IPSL-CM6A-LR" = list(t = ip_ssp$t, d = ip_ssp$a - ref_mean(ip_hist))
)
## ---- hosing anomalies: hos - mean(con) of same model ----
hos <- lapply(names(cols), function(m) {
  h <- read_hos(m); list(t = h$t, d = h$hos - mean(h$con, na.rm = TRUE))
}); names(hos) <- names(cols)

## shared y across everything so the two panels are directly comparable
ylim <- range(unlist(lapply(names(cols), function(m) c(proj[[m]]$d, hos[[m]]$d))), na.rm = TRUE)

## ---- plot ----
pdf(out, width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 5, 3.2), oma = c(0, 1.5, 2.5, 0))
for (m in names(cols)) {
  col <- cols[m]
  ## bottom axis: projection, calendar years 2015-2100 (solid)
  plot(proj[[m]]$t, proj[[m]]$d, type = "l", col = col, lwd = 2, lty = 1,
       xlim = c(2015, 2100), ylim = ylim, axes = FALSE,
       xlab = "Year (SSP1-2.6, bottom)", ylab = "")
  axis(1); axis(2); box()
  abline(h = 0, col = "grey70")
  title(main = m, line = 3.3)
  ## top axis: hosing, years since hosing start 0-150 (dashed)
  par(new = TRUE)
  plot(hos[[m]]$t, hos[[m]]$d, type = "l", col = col, lwd = 1.5, lty = 2,
       xlim = c(0, 150), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  axis(3)
  mtext("Years since hosing start (U03, top)", side = 3, line = 1.6, cex = 0.75, col = "grey30")
  if (m == "EC-Earth3")
    legend("bottomleft", bty = "n",
           legend = c("SSP1-2.6 (solid, bottom axis)", "U03 hosing (dashed, top axis)"),
           col = col, lty = c(1, 2), lwd = c(2, 1.5))
}
mtext("AMOC anomaly: SSP1-2.6 projection vs U03 hosing", outer = TRUE, cex = 1.15, font = 2, line = 0.8)
mtext("AMOC anomaly (Sv, vs preindustrial reference)", side = 2, outer = TRUE, line = 0.1)
dev.off()
cat("wrote", out, "\n\n")

## ---- diagnostics ----
for (m in names(cols)) {
  pr <- range(proj[[m]]$d, na.rm = TRUE); ho <- range(hos[[m]]$d, na.rm = TRUE)
  inside <- ho[1] <= pr[1] && ho[2] >= pr[2]
  target <- pr[1]                                        # most negative projected anomaly
  hit <- which(hos[[m]]$d <= target)                     # first hosing year reaching it
  hit_yr <- if (length(hit)) hos[[m]]$t[hit[1]] else NA
  cat(sprintf("%-13s ΔSv_proj [% .2f, % .2f]  ΔSv_hos [% .2f, % .2f]\n", m, pr[1], pr[2], ho[1], ho[2]))
  cat(sprintf("%-13s projection range fully inside hosing range: %s\n", "", inside))
  cat(sprintf("%-13s hosing reaches most-negative proj (%.2f Sv) first at year %s since hosing start\n\n",
              "", target, ifelse(is.na(hit_yr), "never", as.character(hit_yr))))
}
