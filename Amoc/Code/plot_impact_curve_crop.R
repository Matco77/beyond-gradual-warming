# Phase 4: crop impact vs AMOC weakening (delta_Sv), both crop specifications side by side.
# Author: Marco Bova
#
# Both specifications are plotted on every panel - spec A (fixed Mar-Jul window, solid) and spec C
# (thermal-time window, dashed) - never spec A alone. Appendix I found that spec A's sign comes
# entirely from a coefficient that is not identified once the window follows phenology, and spec C
# cannot identify the thermal dose by construction; showing only one would misrepresent the crop
# result as settled when it is not. The per-hosing-year band (Appendix I) is drawn only for spec A,
# on IPSL-CM6A-LR and EC-Earth3 - the two models it was computed for - as the range across the
# individual years that populate each bin, not an analytic error bar.
#
# Two reference markers, from TWO DIFFERENT sources - do not conflate them. The vertical line is
# each model's OWN ssp126 CMIP6 run's real projected AMOC decline (ocean circulation diagnostic,
# not ISIMIP3b - see plot_impact_common.R::ssp126_target_sv). The horizontal line is the effect the
# ISIMIP3b branch (Appendix J, bias-adjusted surface fields) computes independently at that
# warming. Where the AMOC curve crosses the vertical line is the prediction the delta-method
# scenario makes for that Sv level; the horizontal line is what the ISIMIP branch says directly.
# Agreement between the two is not assumed by construction - it is the actual check this figure
# exists to show.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

A <- totalise(fread(file.path(d, "amoc_impact_crop_eu.csv")),   c("model", "bin_id", "delta_sv", "n_years"))
C <- totalise(fread(file.path(d, "amoc_impact_cropgw_eu.csv")), c("model", "bin_id", "delta_sv", "n_years"))
BND <- fread(file.path(d, "amoc_band_total.csv"))[branch == "crop"]         # spec A only, IPSL+EC
IsoA <- totalise(fread(file.path(d, "isimip_impact_crop_eu.csv")),   "model")
IsoC <- totalise(fread(file.path(d, "isimip_impact_cropgw_eu.csv")), "model")

panel <- function(mdl, ylim, show_legend = FALSE) {
  # NOTE: the loop/argument variable is deliberately NOT called `model` - inside data.table's `[`,
  # `A[model == model]` resolves the right-hand `model` to the COLUMN itself before the calling
  # scope, so it silently becomes `A[model_col == model_col]` (always TRUE, every row). Verified
  # this actually happens (not a hypothetical) before shipping the plot.
  a <- A[model == mdl][order(delta_sv)]; c_ <- C[model == mdl][order(delta_sv)]
  col <- COL[mdl]; rgbcol <- col2rgb(col) / 255       # also used by the legend swatch below
  # per-model x-range (not the shared range across models): HadGEM3-GC3-1MM's bins run to -14.3 Sv
  # while IPSL/EC-Earth3 stop near -9.3, so a shared x-axis would leave most panels mostly empty.
  # The y-axis (effect size) IS shared (ylim, passed in) so magnitudes stay comparable across panels.
  xlim <- range(a$delta_sv, c_$delta_sv); xlim <- xlim + c(-1, 1) * 0.08 * diff(xlim)
  plot(NA, xlim = xlim, ylim = ylim, xlab = expression(Delta*"Sv"), ylab = "yield effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)

  bnd <- BND[model == mdl]
  if (nrow(bnd)) {
    for (i in seq_len(nrow(bnd)))
      segments(bnd$bin_id[i], pct(bnd$lo[i]), bnd$bin_id[i], pct(bnd$hi[i]),
               col = rgb(rgbcol[1], rgbcol[2], rgbcol[3], 0.35), lwd = 5, lend = 1)
  }
  lines(a$delta_sv, a$pct, col = col, lwd = 2); points(a$delta_sv, a$pct, col = col, pch = 16, cex = 1.1)
  lines(c_$delta_sv, c_$pct, col = col, lwd = 1.4, lty = 2); points(c_$delta_sv, c_$pct, col = col, pch = 21, bg = "white", cex = 0.9)
  label_n(a$delta_sv, a$pct, a$n_years, col)

  key <- ISIMIP_KEY[mdl]
  if (!is.na(key) && key %in% IsoA$model) {
    tgt <- ssp126_target_sv(mdl)
    abline(v = tgt, col = "grey40", lty = 3)
    yA <- IsoA[model == key]$pct; yC <- IsoC[model == key]$pct
    abline(h = yA, col = "grey40", lty = 2); abline(h = yC, col = "grey40", lty = 3)
    text(tgt, ylim[2], sprintf("ssp126 %.1f Sv", tgt), col = "grey30", cex = 0.65, pos = 2, srt = 90, offset = 0.3)
  }
  title(main = mdl, cex.main = 1)
  if (show_legend) legend("topright", bg = "white", box.col = NA, cex = 0.65,
    legend = c("spec A (Mar-Jul)", "spec C (thermal window)", "spec A per-year range", "ssp126 AMOC target", "ISIMIP effect (A / C)"),
    col = c(col, col, rgb(rgbcol[1], rgbcol[2], rgbcol[3], 0.5), "grey40", "grey40"), lty = c(1,2,NA,3,2), pch = c(16,21,15,NA,NA),
    lwd = c(2,1.4,NA,1,1), pt.cex = c(1,0.9,1.6,1,1))
}

xr <- range(A$delta_sv); yr <- range(pct(c(A$dln, C$dln, BND$lo, BND$hi)), na.rm = TRUE)
yr <- yr + c(-1, 1) * 0.08 * diff(yr)

pdf(file.path(out, "impact_curve_crop.pdf"), width = 13, height = 8)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
for (i in seq_along(MODELS)) panel(MODELS[i], yr, show_legend = (i == 1))

## ---- comparison panel: all four models overlaid, both specs, no band (too busy) ----------------
plot(NA, xlim = xr, ylim = yr, xlab = expression(Delta*"Sv"), ylab = "yield effect (%)")
grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)
for (m in MODELS) {
  a <- A[model == m][order(delta_sv)]; c_ <- C[model == m][order(delta_sv)]
  lines(a$delta_sv, a$pct, col = COL[m], lwd = 2); points(a$delta_sv, a$pct, col = COL[m], pch = 16, cex = 0.8)
  lines(c_$delta_sv, c_$pct, col = COL[m], lwd = 1.2, lty = 2)
}
title(main = "All models - comparison", cex.main = 1)
legend("topright", bg = "white", box.col = NA, cex = 0.7, legend = MODELS, col = COL[MODELS], lwd = 2, pch = 16)

plot.new()
text(0.5, 0.9, "Crop yield effect vs AMOC-weakening bin (delta_Sv)", cex = 1.05, font = 2)
text(0.5, 0.78, "solid = spec A (fixed Mar-Jul window, benchmark)", cex = 0.8)
text(0.5, 0.70, "dashed = spec C (thermal-time window)", cex = 0.8)
text(0.5, 0.58, "neither total is a standalone estimate (Appendix I) -", cex = 0.75, col = "grey30")
text(0.5, 0.51, "reported together so the specification sensitivity is visible", cex = 0.75, col = "grey30")
text(0.5, 0.36, "shaded band: range across the bin's individual hosing years,", cex = 0.7, col = "grey30")
text(0.5, 0.29, "spec A only, IPSL-CM6A-LR and EC-Earth3 (Appendix I)", cex = 0.7, col = "grey30")
text(0.5, 0.14, "vertical: ocean circulation of that SAME run (ISIMIP has no ocean output);", cex = 0.7, col = "grey30")
text(0.5, 0.07, "horizontal: the ISIMIP3b branch's own effect at that warming (Appendix J)", cex = 0.7, col = "grey30")
dev.off()
cat("wrote", file.path(out, "impact_curve_crop.pdf"), "\n")
