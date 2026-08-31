# Shared constants and helpers for the Phase 4 impact-vs-Sv curves (crop and energy).
# Author: Marco Bova
suppressMessages({library(ncdf4); library(data.table)})
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

# The AMOC decline each model's OWN ssp126 run projects, read from the k10 target files rather
# than hardcoded, so a re-run of that pipeline cannot silently drift from the number plotted here.
#
# NOT computed from the ISIMIP3b files: ISIMIP3b (Appendix J) bias-adjusts only the SURFACE
# atmosphere (tas/tasmin/tasmax/pr) and carries no ocean circulation output, so this target comes
# from the ocean diagnostic instead - Terhaar's msftyz-based reconstruction for IPSL-CM6A-LR, an
# own vo-based reconstruction for EC-Earth3 (Amoc/Code/anomaly/plot_amoc_sv_ssp126.R), as
# min(ssp126 AMOC) - mean(1850-1900 historical AMOC).
#
# But it IS the same underlying simulation as the Phase 3 ISIMIP branch, not an unrelated one -
# checked against both file sets' naming, which encodes CMIP6 model+member+experiment: Terhaar's
# file is "amoc_ssp126_IPSL_IPSL-CM6A-LR_r1i1p1f1.nc", the EC-Earth3 reconstruction's source
# attribute reads "vo (Omon, gn), EC-Earth3 r1i1p1f1", and the ISIMIP3b downloads are
# "ipsl-cm6a-lr_r1i1p1f1_..." / "ec-earth3_r1i1p1f1_...". Same model, same r1i1p1f1 ensemble
# member, same ssp126 experiment on both sides. So the vertical line (ocean circulation) and the
# horizontal line (Appendix J's effect, from the bias-adjusted atmosphere) are two diagnostics of
# the SAME simulated future, not two independent ones - which is what makes the crossing point a
# meaningful check rather than a comparison of unrelated quantities. HadGEM has no target file (no
# target computed for it in that earlier pipeline stage).
ssp126_target_sv <- function(model) {
  f <- file.path(d, sprintf("amoc_effect_%s_u03_target_k10.nc", model))
  if (!file.exists(f)) return(NA_real_)
  nc <- nc_open(f); v <- ncatt_get(nc, 0, "target_delta_sv")$value; nc_close(nc)
  as.numeric(sub(" Sv$", "", v))
}

# sum a component-level effect table (model, bin_id[, delta_sv, n_years], component, dln) to one
# TOTAL row per group, in % (pct()) - used for both the AMOC bin curve and the ISIMIP single point
totalise <- function(dt, by) dt[, .(dln = sum(dln)), by = by][, pct := pct(dln)][]

# n_years labels: small, offset above each point, so they read as annotation, not as ticks
label_n <- function(x, y, n, col) text(x, y, labels = n, pos = 3, cex = 0.55, col = col, offset = 0.25)
