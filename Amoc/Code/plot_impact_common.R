# Shared constants and helpers for the Phase 4 impact-vs-Sv curves (crop and energy).
# Author: Marco Bova
suppressMessages({library(data.table)})
d   <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets")
out <- path.expand("~/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/Plots_Anomaly")

# same 4 colours / model order as plot_amoc_sv_hosing.R, so a reader sees one palette across the thesis
MODELS <- c("IPSL-CM6A-LR", "EC-Earth3", "HadGEM3-GC3-1LL", "HadGEM3-GC3-1MM")
COL    <- c("IPSL-CM6A-LR" = "#e7298a", "EC-Earth3" = "#1b9e77",
           "HadGEM3-GC3-1LL" = "#d95f02", "HadGEM3-GC3-1MM" = "#7570b3")
# ISIMIP scenario files key models in lowercase (ISIMIP's own file-naming convention); only these two
# were run through the ISIMIP branch (Appendix J) - HadGEM has no ISIMIP-forcing counterpart here.
ISIMIP_KEY <- c("IPSL-CM6A-LR" = "ipsl-cm6a-lr", "EC-Earth3" = "ec-earth3")

pct <- function(x) 100 * (exp(x) - 1)          # log-effect -> % change, used throughout Phase 2/3

# No ssp126-AMOC-level marker here (there was one; dropped - Appendix J, J.4d). It would have paired
# a moment along the ssp126 trajectory (the single most negative year, decades away from the
# horizontal line's 2071-2100 window - J.4c) with a quantity the AMOC bins do not isolate the same
# way: the bins are a pure hosing effect against an unforced piControl baseline (Appendix C), while
# ssp126's AMOC weakening happens alongside real background greenhouse warming. Not the same
# physical quantity at two moments - two different quantities - so no crossing point built from it
# would have been a meaningful check.

# sum a component-level effect table (model, bin_id[, delta_sv, n_years], component, dln) to one
# TOTAL row per group, in % (pct()) - used for both the AMOC bin curve and the ISIMIP single point
totalise <- function(dt, by) dt[, .(dln = sum(dln)), by = by][, pct := pct(dln)][]

# n_years labels: small, offset above each point, so they read as annotation, not as ticks
label_n <- function(x, y, n, col) text(x, y, labels = n, pos = 3, cex = 0.55, col = col, offset = 0.25)
