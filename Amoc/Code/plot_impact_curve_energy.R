# Phase 4: energy impact vs AMOC weakening (delta_Sv), one page per fuel.
# Author: Marco Bova
#
# Single specification here - unlike the crop branch, energy carries no window-choice ambiguity
# (Appendix I, H.7): the windows are calendar constructions (the heating year, Oct-Mar, JJA) that
# correspond to actual consumption cycles, not to a biological process whose timing shifts with
# temperature. The per-hosing-year band (Appendix I) is drawn for IPSL-CM6A-LR and EC-Earth3, the
# two models it was computed for.
#
# ISIMIP overlay (Appendix J.5, plot_impact_common.R): filled triangles, IPSL-CM6A-LR and
# EC-Earth3 panels only, at the delta_Sv bins isimip_bin_fields.R shares with NAHosMIP. Residual
# caveat: the two branches' delta_Sv are not on identical footing (piControl-relative vs
# historical-relative baseline) - see plot_impact_common.R header.
#
# ISIMIP per-year band (Appendix J.6, isimip_impact_band.R): a narrow bar (thinner than the NAHosMIP
# band) at the ISIMIP bin's own delta_sv, spanning the range of the final effect across that bin's
# individual ssp126 years - same construction as the hosing-year band, applied to
# isimip_bin_fields.R's mode = "year" fields instead of NAHosMIP's. Both bands sit at their own
# curve's delta_sv (bin-mean AMOC weakening), not the integer bin label.
#
# NAMING NOTE, applies throughout this file: loop/argument variables are `mdl` and `fl`, never
# `model` or `fuel`. Inside data.table's `[`, a bare `E[model == model]` or `E[fuel == fuel]`
# resolves the right-hand name to the COLUMN itself before the calling scope, so it silently
# matches every row regardless of the argument's value. Verified this happens (not hypothetical,
# checked in isolation) for both `model` and `fuel` before shipping the plot.
source(path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Code/plot_impact_common.R"))

E    <- totalise(fread(file.path(d, "amoc_impact_energy_eu.csv")), c("model", "bin_id", "delta_sv", "n_years", "fuel"))
BND  <- fread(file.path(d, "amoc_band_total.csv"))         # branch %in% {"energy Electricity","energy Natural gas"}
ISB  <- fread(file.path(d, "isimip_band_total.csv"))[crop == "ALL"]   # ssp126 per-year band (crop col is "ALL" for energy; the crop branches also carry per-crop rows)
IS   <- totalise(fread(file.path(d, "isimip_bin_impact_energy_eu.csv")), c("model", "bin_id", "delta_sv", "n_years", "fuel"))
FUELS <- c("Electricity", "Natural gas")
BND_KEY <- c(Electricity = "energy Electricity", `Natural gas` = "energy Natural gas")

panel <- function(mdl, fl, ylim) {
  z   <- E[model == mdl & fuel == fl][order(delta_sv)]
  col <- COL[mdl]
  xlim <- range(z$delta_sv); xlim <- xlim + c(-1, 1) * 0.08 * diff(xlim)
  plot(NA, xlim = xlim, ylim = ylim, xlab = expression(Delta*"Sv"), ylab = "demand effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)

  # every band is drawn at its own curve's delta_sv (bin-mean AMOC weakening), NOT the bin label
  # bin_id, so the bar sits under the point it belongs to.
  bnd <- merge(BND[model == mdl & branch == BND_KEY[fl]], z[, .(bin_id, delta_sv)], by = "bin_id")
  if (nrow(bnd))
    for (i in seq_len(nrow(bnd)))
      segments(bnd$delta_sv[i], pct(bnd$lo[i]), bnd$delta_sv[i], pct(bnd$hi[i]),
               col = tint(col, 0.35), lwd = 5, lend = 1)

  lines(z$delta_sv, z$pct, col = col, lwd = 2); points(z$delta_sv, z$pct, col = col, pch = 16, cex = 1.1)
  label_n(z$delta_sv, z$pct, z$n_years, col)

  is <- IS[model == mdl & fuel == fl][order(delta_sv)]
  isb <- merge(ISB[model == mdl & branch == BND_KEY[fl]], is[, .(bin_id, delta_sv)], by = "bin_id")
  if (nrow(isb))
    for (i in seq_len(nrow(isb)))
      segments(isb$delta_sv[i], pct(isb$lo[i]), isb$delta_sv[i], pct(isb$hi[i]),
               col = tint(col, 0.5), lwd = 3, lend = 1)

  if (nrow(is)) {
    lines(is$delta_sv, is$pct, col = col, lwd = 1.4, lty = ISIMIP_LTY)
    points(is$delta_sv, is$pct, col = col, pch = ISIMIP_PCH, cex = 1.1)
    label_n(is$delta_sv, is$pct, is$n_years, col)
  }
  title(main = mdl, cex.main = 1)
}

pdf(file.path(out, "impact_curve_energy.pdf"), width = 13, height = 8)
for (fl in FUELS) {
  sub <- E[fuel == fl]
  bl  <- BND[branch == BND_KEY[fl]]
  isb <- ISB[branch == BND_KEY[fl]]
  isf <- IS[fuel == fl]
  yr  <- range(pct(c(sub$dln, bl$lo, bl$hi, isf$dln, isb$lo, isb$hi)), na.rm = TRUE)
  yr  <- yr + c(-1, 1) * 0.08 * diff(yr)
  xr  <- range(sub$delta_sv)

  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
  for (mdl in MODELS) panel(mdl, fl, yr)

  plot(NA, xlim = xr, ylim = yr, xlab = expression(Delta*"Sv"), ylab = "demand effect (%)")
  grid(col = "grey90"); abline(h = 0, col = "grey60", lty = 3)
  for (mdl in MODELS) {
    z <- sub[model == mdl][order(delta_sv)]
    lines(z$delta_sv, z$pct, col = COL[mdl], lwd = 2); points(z$delta_sv, z$pct, col = COL[mdl], pch = 16, cex = 0.8)
    is <- isf[model == mdl][order(delta_sv)]
    if (nrow(is)) points(is$delta_sv, is$pct, col = COL[mdl], pch = ISIMIP_PCH, cex = 0.9)
  }
  title(main = "All models - comparison", cex.main = 1)
  legend("topleft", bg = "white", box.col = NA, cex = 0.7,
    legend = c(MODELS, "ISIMIP-bin (ssp126)"), col = c(COL[MODELS], "grey30"),
    lwd = c(rep(2, length(MODELS)), NA), pch = c(rep(16, length(MODELS)), ISIMIP_PCH))

  plot.new()
  text(0.5, 0.9, sprintf("%s demand effect vs AMOC-weakening bin (delta_Sv)", fl), cex = 1.0, font = 2)
  text(0.5, 0.76, "single specification (Appendix H, H.7): the energy windows are calendar", cex = 0.72, col = "grey30")
  text(0.5, 0.70, "constructions matching real consumption cycles, no window ambiguity", cex = 0.72, col = "grey30")
  text(0.5, 0.52, "shaded band: range across the bin's individual hosing years,", cex = 0.7, col = "grey30")
  text(0.5, 0.45, "IPSL-CM6A-LR and EC-Earth3 only (Appendix I)", cex = 0.7, col = "grey30")
  text(0.5, 0.39, "narrow bar (thinner): same range across the ISIMIP bin's ssp126 years (Appendix J.6)", cex = 0.62, col = "grey30")
  text(0.5, 0.33, "triangles: ISIMIP3b/ssp126, binned by year-level AMOC (Appendix J.5),", cex = 0.65, col = "grey30")
  text(0.5, 0.27, "at the delta_Sv levels it shares with NAHosMIP - IPSL and EC-Earth3 only", cex = 0.65, col = "grey30")
  text(0.5, 0.21, "caveat: delta_Sv baseline differs (piControl vs historical mean), see Appendix J.5", cex = 0.62, col = "grey40")
  text(0.5, 0.08, "caveat (Appendix I, I.4): the AMOC bins extrapolate the linear response", cex = 0.65, col = "grey40")
  text(0.5, 0.02, "function up to ~9 s.d. beyond the range that identifies the coefficients", cex = 0.65, col = "grey40")
}
dev.off()
cat("wrote", file.path(out, "impact_curve_energy.pdf"), "\n")
