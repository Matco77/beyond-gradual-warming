# Phase 4: energy impact vs AMOC weakening (delta_Sv), one page per fuel.
# Author: Marco Bova
#
# Single specification here - unlike the crop branch, energy carries no window-choice ambiguity
# (Appendix I, H.7): the windows are calendar constructions (the heating year, Oct-Mar, JJA) that
# correspond to actual consumption cycles, not to a biological process whose timing shifts with
# temperature. The per-hosing-year band (Appendix I) is drawn for IPSL-CM6A-LR and EC-Earth3, the
# two models it was computed for.
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
#
# NAMING NOTE, applies throughout this file: loop/argument variables are `mdl` and `fl`, never
# `model` or `fuel`. Inside data.table's `[`, a bare `E[model == model]` or `E[fuel == fuel]`
# resolves the right-hand name to the COLUMN itself before the calling scope, so it silently
# matches every row regardless of the argument's value. Verified this happens (not hypothetical,
# checked in isolation) for both `model` and `fuel` before shipping the plot.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

E    <- totalise(fread(file.path(d, "amoc_impact_energy_eu.csv")), c("model", "bin_id", "delta_sv", "n_years", "fuel"))
BND  <- fread(file.path(d, "amoc_band_total.csv"))         # branch %in% {"energy Electricity","energy Natural gas"}
FUELS <- c("Electricity", "Natural gas")
BND_KEY <- c(Electricity = "energy Electricity", `Natural gas` = "energy Natural gas")

panel <- function(mdl, fl, ylim) {
  z   <- E[model == mdl & fuel == fl][order(delta_sv)]
  col <- COL[mdl]
  xlim <- range(z$delta_sv); xlim <- xlim + c(-1, 1) * 0.08 * diff(xlim)
  plot(NA, xlim = xlim, ylim = ylim, xlab = expression(Delta*"Sv"), ylab = "demand effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)

  bnd <- BND[model == mdl & branch == BND_KEY[fl]]
  if (nrow(bnd))
    for (i in seq_len(nrow(bnd)))
      segments(bnd$bin_id[i], pct(bnd$lo[i]), bnd$bin_id[i], pct(bnd$hi[i]),
               col = tint(col, 0.35), lwd = 5, lend = 1)

  lines(z$delta_sv, z$pct, col = col, lwd = 2); points(z$delta_sv, z$pct, col = col, pch = 16, cex = 1.1)
  label_n(z$delta_sv, z$pct, z$n_years, col)
  title(main = mdl, cex.main = 1)
}

pdf(file.path(out, "impact_curve_energy.pdf"), width = 13, height = 8)
for (fl in FUELS) {
  sub <- E[fuel == fl]
  bl  <- BND[branch == BND_KEY[fl]]
  yr  <- range(pct(c(sub$dln, bl$lo, bl$hi)), na.rm = TRUE)
  yr  <- yr + c(-1, 1) * 0.08 * diff(yr)
  xr  <- range(sub$delta_sv)

  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
  for (mdl in MODELS) panel(mdl, fl, yr)

  plot(NA, xlim = xr, ylim = yr, xlab = expression(Delta*"Sv"), ylab = "demand effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)
  for (mdl in MODELS) {
    z <- sub[model == mdl][order(delta_sv)]
    lines(z$delta_sv, z$pct, col = COL[mdl], lwd = 2); points(z$delta_sv, z$pct, col = COL[mdl], pch = 16, cex = 0.8)
  }
  title(main = "All models - comparison", cex.main = 1)
  legend("topleft", bg = "white", box.col = NA, cex = 0.7, legend = MODELS, col = COL[MODELS], lwd = 2, pch = 16)

  plot.new()
  text(0.5, 0.9, sprintf("%s demand effect vs AMOC-weakening bin (delta_Sv)", fl), cex = 1.0, font = 2)
  text(0.5, 0.76, "single specification (Appendix H, H.7): the energy windows are calendar", cex = 0.72, col = "grey30")
  text(0.5, 0.70, "constructions matching real consumption cycles, no window ambiguity", cex = 0.72, col = "grey30")
  text(0.5, 0.52, "shaded band: range across the bin's individual hosing years,", cex = 0.7, col = "grey30")
  text(0.5, 0.45, "IPSL-CM6A-LR and EC-Earth3 only (Appendix I)", cex = 0.7, col = "grey30")
  text(0.5, 0.27, "no ISIMIP marker on this figure - the AMOC bins and the ISIMIP", cex = 0.65, col = "grey30")
  text(0.5, 0.21, "branch measure different physical quantities (Appendix J, J.4e)", cex = 0.65, col = "grey30")
  text(0.5, 0.08, "caveat (Appendix I, I.4): the AMOC bins extrapolate the linear response", cex = 0.65, col = "grey40")
  text(0.5, 0.02, "function up to ~9 s.d. beyond the range that identifies the coefficients", cex = 0.65, col = "grey40")
}
dev.off()
cat("wrote", file.path(out, "impact_curve_energy.pdf"), "\n")
