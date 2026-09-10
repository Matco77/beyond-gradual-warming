## AMOC at 26N (Sv): hosing decline vs control, EC-Earth3 + HadGEM3.
## Reads M26 netCDFs as-is; no derived variables.
library(ncdf4)

dir  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets/M26"
out  <- "/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly/amoc_sv_hosing.pdf"

models <- c("EC-Earth3", "HadGEM3-GC3-1LL", "HadGEM3-GC3-1MM", "IPSL-CM6A-LR")
cols   <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

read1 <- function(m) {
  nc <- nc_open(file.path(dir, paste0("M26_", m, ".nc")))
  d  <- list(t   = ncvar_get(nc, "time"),
             hos = ncvar_get(nc, paste0("hos_", m)),
             con = ncvar_get(nc, paste0("con_", m)))
  nc_close(nc)
  ## netCDF fill values in the tail aren't flagged; drop non-physical Sv.
  d$hos[abs(d$hos) > 100] <- NA
  d$con[abs(d$con) > 100] <- NA
  d
}
dat <- lapply(models, read1); names(dat) <- models

## IPSL con runs to 500 yr but its hos is valid only to ~100; cap x at the
## longest hos-valid series so the short runs aren't squashed into the left.
xlim <- range(sapply(dat, function(d) range(d$t[is.finite(d$hos)], na.rm = TRUE)))
ylim <- range(sapply(dat, function(d) range(c(d$hos, d$con), na.rm = TRUE)))

pdf(out, width = 9, height = 6)
par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = xlim, ylim = ylim, xlab = "Years since hosing start",
     ylab = "AMOC strength at 26N (Sv)", main = "")
title(main = "AMOC decline under hosing vs control (no hosing)", line = 2)
mtext("NAHosMIP — hosing experiment U03 (0.3 Sv North Atlantic freshwater)",
      side = 3, line = 0.5, cex = 0.85, col = "grey30")
grid(col = "grey90")
for (i in seq_along(models)) {
  d <- dat[[models[i]]]
  lines(d$t, d$con, col = cols[i], lty = 2, lwd = 1)   # control, no hosing
  lines(d$t, d$hos, col = cols[i], lty = 1, lwd = 2)   # hosing decline
}
legend("bottomleft", bg = "white", box.col = NA, cex = 0.9,   # opaque bg: readable over crossing curves
       legend = c(models, "hosing (solid)", "control / no hosing (dashed)"),
       col = c(cols, "black", "black"),
       lty = c(1, 1, 1, 1, 1, 2), lwd = c(2, 2, 2, 2, 2, 1))
dev.off()
cat("wrote", out, "\n")
