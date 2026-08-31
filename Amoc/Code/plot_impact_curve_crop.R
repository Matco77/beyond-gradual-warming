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
# NO ISIMIP reference marker on this figure - neither a vertical Sv line nor a horizontal effect
# line. Both were tried and dropped (Appendix J, J.4d-J.4e): the AMOC bins isolate a PURE hosing
# effect against an unforced piControl baseline (Appendix C), while the ISIMIP/ssp126 branch's
# effect (Appendix J) reflects hosing-like circulation change MIXED with real background greenhouse
# warming. They are not the same physical quantity, so overlaying either one on this curve implies
# a comparability that does not hold - not as a vertical marker at a mismatched moment (J.4c-d), and
# not as a horizontal marker either, since a flat reference line across the whole panel still visually
# invites reading it against the curve. The ISIMIP branch's numbers are a real, standalone result;
# they are reported in Appendix J's own tables, not on this figure.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

A <- totalise(fread(file.path(d, "amoc_impact_crop_eu.csv")),   c("model", "bin_id", "delta_sv", "n_years"))
C <- totalise(fread(file.path(d, "amoc_impact_cropgw_eu.csv")), c("model", "bin_id", "delta_sv", "n_years"))
BND <- fread(file.path(d, "amoc_band_total.csv"))[branch == "crop"]         # spec A only, IPSL+EC

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

  title(main = mdl, cex.main = 1)
  if (show_legend) legend("topright", bg = "white", box.col = NA, cex = 0.65,
    legend = c("spec A (Mar-Jul)", "spec C (thermal window)", "spec A per-year range"),
    col = c(col, col, rgb(rgbcol[1], rgbcol[2], rgbcol[3], 0.5)), lty = c(1,2,NA), pch = c(16,21,15),
    lwd = c(2,1.4,NA), pt.cex = c(1,0.9,1.6))
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
text(0.5, 0.34, "shaded band: range across the bin's individual hosing years,", cex = 0.7, col = "grey30")
text(0.5, 0.27, "spec A only, IPSL-CM6A-LR and EC-Earth3 (Appendix I)", cex = 0.7, col = "grey30")
text(0.5, 0.12, "no ISIMIP marker on this figure - the AMOC bins and the ISIMIP", cex = 0.65, col = "grey30")
text(0.5, 0.06, "branch measure different physical quantities (Appendix J, J.4e)", cex = 0.65, col = "grey30")
dev.off()
cat("wrote", file.path(out, "impact_curve_crop.pdf"), "\n")
