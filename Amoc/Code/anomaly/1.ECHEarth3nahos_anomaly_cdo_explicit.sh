#!/usr/bin/env bash
###############################################################################
# NAHosMIP hosing  MINUS  piControl  monthly anomaly  --  EC-Earth3
#
# Same method as the HadGEM3-LL / MM scripts:
#   hosing monthly  -  piControl monthly CLIMATOLOGY (ymonmean)  =  anomaly
#   matched by calendar month, per native grid cell (gr).
#
# BASELINE WINDOW NOTE (read once, write one sentence in thesis):
#   EC-Earth3 hosing branched from piControl-spinup at branch_time_in_parent=0
#   (days since 1850-01-01). That spinup epoch is NOT among the piControl files
#   on disk (which cover model years 2259-2759). There is therefore no
#   contemporaneous parallel chunk to select, as there was for HadGEM3 (2050-
#   2149). We use the FULL available piControl (2259-2759) as the unforced
#   reference climatology. After spinup EC-Earth3 piControl is quasi-stationary,
#   so its long-term mean seasonal cycle is a sound reference state.
#   (NAHosMIP protocol: Jackson et al., 2023, Geosci. Model Dev. 16:1975-1995.)
#
# Variables:
#   tas, pr        -> native EC-Earth3 hosing files
#   tasmin, tasmax -> RECONSTRUCTED files (per-cell per-month OLS on piControl;
#                     Weedon et al. 2010). Reconstructed absolute fields used as
#                     the hosing side so the anomaly is built identically to the
#                     other variables.
#
# Two hosing experiments: hos-g01-hos (50 yr) and hos-u03-hos (100 yr).
#
# Requires CDO:  brew install cdo
###############################################################################
set -euo pipefail

# Talk on failure: if any command fails, print which line and command died,
# instead of exiting silently (the thing that made debugging hard before).
trap 'echo "ERROR: script failed at line $LINENO running: $BASH_COMMAND" >&2' ERR

ROOT="/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
PIC_DIR="${ROOT}/CMIP6_piControl/EC-Earth3"
HOS_DIR="${ROOT}/NAHosMIP/EC-Earth3"
OUT_DIR="${ROOT}/anomaly_output/ECHearth3_anomaly"
TMP="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$TMP"

# ---- piControl monthly file globs (501 chunk files, one year each) ----------
PIC_TAS_GLOB="${PIC_DIR}/tas_Amon_EC-Earth3_piControl_r1i1p1f1_gr_*.nc"
PIC_PR_GLOB="${PIC_DIR}/pr_Amon_EC-Earth3_piControl_r1i1p1f1_gr_*.nc"
PIC_TASMIN_GLOB="${PIC_DIR}/tasmin_Amon_EC-Earth3_piControl_r1i1p1f1_gr_*.nc"
PIC_TASMAX_GLOB="${PIC_DIR}/tasmax_Amon_EC-Earth3_piControl_r1i1p1f1_gr_*.nc"

remap_method () {          # pr conserves, temperatures bilinear
  case "$1" in
    pr) echo "remapcon" ;;
    *)  echo "remapbil" ;;
  esac
}

# build the piControl climatology ONCE per variable (whole run -> 12 months)
# args: varname  pic_glob   -> writes ${TMP}/${var}_pic_clim.nc
build_pic_clim () {
  local V="$1" GLOB="$2"
  local merged="${TMP}/${V}_pic_merged.nc"
  local clim="${TMP}/${V}_pic_clim.nc"
  if [ -f "$clim" ]; then echo "  pic clim $V already built"; return; fi
  echo "  [pic] mergetime ${V} (this is the slow step, 501 files) ..."
  cdo -O mergetime ${GLOB} "$merged"
  echo "  [pic] ymonmean ${V} -> 12-month climatology"
  cdo -O ymonmean "$merged" "$clim"
  rm -f "$merged"
}

# args: varname  hosing_file  pic_clim  label
do_var () {
  local V="$1" HOS="$2" CLIM="$3" LABEL="$4"
  echo "=========== $V  ($LABEL) ==========="
  local clim_use="$CLIM"
  local anom="${OUT_DIR}/${V}_Amon_EC-Earth3_${LABEL}_anomaly.nc"

  # grid check by SIZE (gridsize = number of cells), not full griddes text.
  # Identical native grids can show float-dust differences in lat/lon labels
  # (45.0000001 vs 45.0) that fool a text diff into a needless remap -- which
  # for pr (conservative) can smear values. If sizes match it IS the same grid:
  # snap the climatology onto the hosing grid with setgrid (exact coordinate
  # copy, NO interpolation, NO error). Only a true size mismatch -> real remap.
  #
  # NOTE: setgrid needs a GRID DESCRIPTION FILE, not a data .nc. So we first
  # dump the hosing grid with `griddes` into a text file, then setgrid that.
  # The size capture is wrapped so a non-zero pipe status cannot abort the
  # script under `set -e`.
  local hos_n pic_n
  hos_n=$(cdo -s griddes "$HOS"  2>/dev/null | awk -F= '/gridsize/{gsub(/ /,"",$2);print $2;exit}' || true)
  pic_n=$(cdo -s griddes "$CLIM" 2>/dev/null | awk -F= '/gridsize/{gsub(/ /,"",$2);print $2;exit}' || true)
  echo "  hosing cells=${hos_n:-?}  pic cells=${pic_n:-?}"

  if [ -n "$hos_n" ] && [ "$hos_n" = "$pic_n" ]; then
    echo "  grids same size -> snap labels with setgrid (no interpolation)"
    local griddes_file="${TMP}/${V}_${LABEL}_hos.griddes"
    cdo -s griddes "$HOS" > "$griddes_file"
    clim_use="${TMP}/${V}_${LABEL}_pic_clim_snap.nc"
    cdo -O setgrid,"$griddes_file" "$CLIM" "$clim_use"
  else
    local M; M=$(remap_method "$V")
    echo "  grids DIFFER in size (${pic_n:-?} -> ${hos_n:-?}) -> real remap with $M"
    clim_use="${TMP}/${V}_${LABEL}_pic_clim_remap.nc"
    cdo -O "$M","$HOS" "$CLIM" "$clim_use"
  fi

  echo "  anomaly = hosing monthly - pic climatology (ymonsub)"
  cdo -O ymonsub "$HOS" "$clim_use" "$anom"
  echo "  -> $(basename "$anom")"

  # pr ALSO gets the multiplicative ratio R = hosing / climatology, on the SAME
  # grid-matched clim_use and 2259-2759 reference as the additive anomaly above.
  # R is what the ISIMIP step applies as pr_AMOC = pr_ISIMIP * R; temperature stays
  # additive only (a ratio of interval-scale K is meaningless). The additive pr
  # anomaly is unchanged, so the diagnostics keep reading it as before.
  # ponytail: ymondiv blows up where climatological precip ~ 0 (deserts / a few dry
  # cells). Stored faithful/unclipped; clamp R to a physical band (e.g. [0.1, 10]) or
  # mask tiny-climatology cells in the ISIMIP-application step (matches HadGEM builder).
  if [ "$V" = "pr" ]; then
    local ratio="${OUT_DIR}/${V}_Amon_EC-Earth3_${LABEL}_ratio.nc"
    echo "  pr ratio = hosing / pic climatology (ymondiv)  (for ISIMIP stressing)"
    cdo -O ymondiv "$HOS" "$clim_use" "$ratio"
    echo "  -> $(basename "$ratio")"
  fi
}

# ---------------------------------------------------------------------------
# Build the four piControl climatologies once
# ---------------------------------------------------------------------------
echo "### Building piControl climatologies (whole run 2259-2759) ###"
build_pic_clim tas    "$PIC_TAS_GLOB"
build_pic_clim pr     "$PIC_PR_GLOB"
build_pic_clim tasmin "$PIC_TASMIN_GLOB"
build_pic_clim tasmax "$PIC_TASMAX_GLOB"

# ---------------------------------------------------------------------------
# Per experiment, per variable: subtract
#   tas/pr  -> native hosing files
#   tasmin/tasmax -> reconstructed files
# ---------------------------------------------------------------------------
for LABEL in hos-g01-hos hos-u03-hos; do
  echo
  echo "############## EXPERIMENT: $LABEL ##############"

  # native tas/pr filenames carry the year span; resolve by glob
  HOS_TAS=$(ls "${HOS_DIR}/tas_Amon_EC-Earth3_${LABEL}_"*.nc | head -1)
  HOS_PR=$(ls  "${HOS_DIR}/pr_Amon_EC-Earth3_${LABEL}_"*.nc  | head -1)
  # reconstructed tasmin/tasmax (no year in name)
  HOS_TASMIN="${HOS_DIR}/tasmin_Amon_EC-Earth3_${LABEL}_reconstructed.nc"
  HOS_TASMAX="${HOS_DIR}/tasmax_Amon_EC-Earth3_${LABEL}_reconstructed.nc"

  do_var tas    "$HOS_TAS"    "${TMP}/tas_pic_clim.nc"    "$LABEL"
  do_var pr     "$HOS_PR"     "${TMP}/pr_pic_clim.nc"     "$LABEL"
  do_var tasmin "$HOS_TASMIN" "${TMP}/tasmin_pic_clim.nc" "$LABEL"
  do_var tasmax "$HOS_TASMAX" "${TMP}/tasmax_pic_clim.nc" "$LABEL"
done

echo
echo "Done. Per-month, per-native-cell anomalies in:"
echo "  $OUT_DIR"
echo "hosing monthly  -  piControl (2259-2759) 12-month climatology."