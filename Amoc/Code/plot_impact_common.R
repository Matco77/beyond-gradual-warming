# Shared constants and helpers for the Phase 4 impact-vs-Sv curves (crop and energy).
# Author: Marco Bova
suppressMessages({library(data.table)})
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
out <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly")

# same 4 colours / model order as plot_amoc_sv_hosing.R, so a reader sees one palette across the thesis
MODELS <- c("IPSL-CM6A-LR", "EC-Earth3", "HadGEM3-GC3-1LL", "HadGEM3-GC3-1MM")
COL    <- c("IPSL-CM6A-LR" = "#e7298a", "EC-Earth3" = "#1b9e77",
           "HadGEM3-GC3-1LL" = "#d95f02", "HadGEM3-GC3-1MM" = "#7570b3")

pct <- function(x) 100 * (exp(x) - 1)          # log-effect -> % change, used throughout Phase 2/3

# These figures carry NO ISIMIP reference marker, vertical or horizontal (there was one of each;
# both dropped - Appendix J, J.4d-J.4e). The AMOC bins are a PURE hosing effect against an unforced
# piControl baseline (Appendix C); the ISIMIP/ssp126 branch's effect (Appendix J) reflects
# hosing-like circulation change MIXED with real background greenhouse warming. Different physical
# quantities, so no marker built from one and overlaid on the other's curve would have been a
# meaningful comparison, whichever way it was drawn.

# sum a component-level effect table (model, bin_id[, delta_sv, n_years], component, dln) to one
# TOTAL row per group, in % (pct()) - used for the AMOC bin curve
totalise <- function(dt, by) dt[, .(dln = sum(dln)), by = by][, pct := pct(dln)][]

# n_years labels: small, offset above each point, so they read as annotation, not as ticks
label_n <- function(x, y, n, col) text(x, y, labels = n, pos = 3, cex = 0.55, col = col, offset = 0.25)

# A tinted, FULLY OPAQUE version of a colour - alpha*colour + (1-alpha)*white, i.e. exactly what
# alpha-blending that colour at `a` would look like against a white background, but written as a
# plain solid fill. Used for the band instead of rgb(..., alpha=a): a PDF alpha fill goes through an
# ExtGState/transparency-group operator (confirmed present in the file: `grep -c '/ca'` > 0), and
# not every viewer/renderer honours it - a plain opaque colour has no such dependency and renders
# identically everywhere. Found by the reader reporting the band invisible in their own viewer while
# it rendered fine in this session's own ghostscript check - a real cross-viewer gap, not assumed.
tint <- function(col, a) { x <- col2rgb(col) / 255; rgb(x[1]*a + (1-a), x[2]*a + (1-a), x[3]*a + (1-a)) }
