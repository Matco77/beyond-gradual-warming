# Phase 4, per-crop view: the ISIMIP-branch impact and its ssp126 per-year band, one panel per crop
# instead of the area-weighted all-crop aggregate of plot_impact_curve_crop.R.
# Author: Marco Bova
#
# WHY THIS EXISTS. plot_impact_curve_crop.R plots one crop response - the area-weighted mean over
# Soft wheat, Durum wheat, Spring barley, Winter barley - so a reader cannot see whether the four
# crops move together or cancel. The beta ARE per crop (regression.R crop_set, one feols each) and
# isimip_impact_band.R now carries the crop split all the way to isimip_band_{summary,total}.csv;
# this figure just stops aggregating it away.
#
# Two specifications, same as the all-crop figure: spec A (fixed Mar-Jul window, branch "crop",
# solid + triangle + band) and spec C (thermal-time window, branch "crop_gddwin", dashed line only -
# no NAHosMIP-side band exists for spec C, Appendix J.6). IPSL-CM6A-LR and EC-Earth3 only, the two
# models with a ssp126 AMOC reconstruction, at the 3 delta_Sv bins they share with NAHosMIP.
#
# Bands and points sit at the bin's delta_sv (bin-mean AMOC weakening), not the integer bin label -
# same convention as plot_impact_curve_crop.R.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

CROPS <- c("Soft wheat", "Durum wheat", "Spring barley", "Winter barley")   # regression.R crop_set order
MDLS  <- c("IPSL-CM6A-LR", "EC-Earth3")

# per-crop TOTAL (summed components) central + band, straight from isimip_impact_band.R's output.
# crop == "ALL" is the area-weighted aggregate (that is plot_impact_curve_crop.R); drop it here.
T <- fread(file.path(d, "isimip_band_total.csv"))[branch %in% c("crop", "crop_gddwin") & crop %in% CROPS]
# delta_sv per (model, bin_id): the bands/points are drawn there, not at the integer bin_id
DSV <- unique(fread(file.path(d, "isimip_bin_impact_crop_country.csv"))[, .(model, bin_id, delta_sv)])
T   <- merge(T, DSV, by = c("model", "bin_id"))

yr <- range(pct(c(T$central, T$lo, T$hi)), na.rm = TRUE); yr <- yr + c(-1, 1) * 0.08 * diff(yr)
xr <- range(T$delta_sv) + c(-0.4, 0.4)

panel <- function(cr, show_legend = FALSE) {
  plot(NA, xlim = xr, ylim = yr, xlab = expression(Delta*"Sv"), ylab = "yield effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)
  for (mdl in MDLS) {
    col <- COL[mdl]
    a <- T[crop == cr & model == mdl & branch == "crop"][order(delta_sv)]
    g <- T[crop == cr & model == mdl & branch == "crop_gddwin"][order(delta_sv)]
    if (nrow(a)) {
      for (i in seq_len(nrow(a)))
        segments(a$delta_sv[i], pct(a$lo[i]), a$delta_sv[i], pct(a$hi[i]), col = tint(col, 0.5), lwd = 4, lend = 1)
      lines(a$delta_sv, pct(a$central), col = col, lwd = 1.6, lty = ISIMIP_LTY)
      points(a$delta_sv, pct(a$central), col = col, pch = ISIMIP_PCH, cex = 1.2)
    }
    if (nrow(g)) lines(g$delta_sv, pct(g$central), col = col, lwd = 1.3, lty = 2)
  }
  title(main = cr, cex.main = 1)
  if (show_legend) legend("bottomleft", bg = "white", box.col = NA, cex = 0.62,
    legend = c(MDLS, "spec A central + band", "spec C central"),
    col = c(COL[MDLS], "grey30", "grey30"), lty = c(NA, NA, ISIMIP_LTY, 2),
    pch = c(15, 15, ISIMIP_PCH, NA), lwd = c(NA, NA, 1.6, 1.3), pt.cex = c(1.4, 1.4, 1.2, NA))
}

pdf(file.path(out, "impact_curve_crop_bycrop.pdf"), width = 13, height = 8)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
for (i in seq_along(CROPS)) panel(CROPS[i], show_legend = (i == 1))

plot.new()
text(0.5, 0.90, "Crop yield effect vs AMOC-weakening bin (delta_Sv), per crop", cex = 1.0, font = 2)
text(0.5, 0.78, "ISIMIP3b/ssp126 branch only (Appendix J.5-J.6), IPSL-CM6A-LR and EC-Earth3", cex = 0.7, col = "grey30")
text(0.5, 0.70, "at the delta_Sv levels shared with NAHosMIP", cex = 0.7, col = "grey30")
text(0.5, 0.56, "solid triangle + dotted line = spec A (fixed Mar-Jul window)", cex = 0.7, col = "grey30")
text(0.5, 0.49, "narrow bar = spec A range across the bin's individual ssp126 years", cex = 0.7, col = "grey30")
text(0.5, 0.42, "dashed line = spec C (thermal-time window), no band (Appendix J.6)", cex = 0.7, col = "grey30")
text(0.5, 0.26, "y-axis shared across the four crop panels so magnitudes compare", cex = 0.66, col = "grey40")
text(0.5, 0.19, "the all-crop area-weighted aggregate is plot_impact_curve_crop.pdf", cex = 0.66, col = "grey40")
dev.off()
cat("wrote", file.path(out, "impact_curve_crop_bycrop.pdf"), "\n")
