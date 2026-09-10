#!/usr/bin/env bash
###############################################################################
# NAHosMIP hosing  MINUS  piControl  monthly anomaly  --  IPSL-CM6A-LR
#
# Clone of 1.ECHEarth3nahos_anomaly_cdo_explicit.sh, same method, IPSL paths.
#   hosing monthly  -  piControl monthly CLIMATOLOGY (ymonmean)  =  anomaly
#   matched by calendar month, per native grid cell (gr).
#
# BASELINE WINDOW NOTE:
#   IPSL piControl on disk covers model years 1850-2349 (the NAHosMIP tar chunk
#   for tas/pr; tasmin/tasmax downloaded from ESGF, same 185001-234912 chunk).
#   We use that FULL 500-yr piControl as the unforced reference climatology
#   (ymonmean over 1850-2349), same reasoning as the EC-Earth3 script: after
#   spinup IPSL piControl is quasi-stationary, so its long-term mean seasonal
#   cycle is a sound reference state. u03-hos branched from piControl (see the
#   :name production path .../PROD/piControl/hos-u03-hos/...).
#   (NAHosMIP protocol: Jackson et al., 2023, Geosci. Model Dev. 16:1975-1995.)
#
# Variables:
#   tas, pr        -> native IPSL hosing files (DRS name: _u03-hos_..._gr_)
#   tasmin, tasmax -> RECONSTRUCTED files (per-cell per-month OLS on piControl;
#                     Weedon et al. 2010), name _hos-u03-hos_reconstructed.nc
#
# One hosing experiment: hos-u03-hos (IPSL has no g01 in the atmospheric tar).
# NOTE: IPSL u03-hos atmospheric data is 150 yr (1850-1999) despite the
#       185001-214912 filename token.
#
# Requires CDO:  brew install cdo
###############################################################################
set -euo pipefail

trap 'echo "ERROR: script failed at line $LINENO running: $BASH_COMMAND" >&2' ERR

ROOT="/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"
PIC_DIR="${ROOT}/CMIP6_piControl/IPSL-CM6A-LR"
HOS_DIR="${ROOT}/NAHosMIP/IPSL-CM6A-LR"
OUT_DIR="${ROOT}/anomaly_output/IPSL_anomaly"
TMP="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$TMP"

# ---- piControl monthly file globs -------------------------------------------
# Restricted to the 185001-234912 chunk so all four variables share ONE uniform
# 500-yr baseline (1850-2349). tasmin/tasmax have 5 ESGF chunks on disk; the
# 2350-3849 ones are not used here (the fit already used only this chunk).
PIC_TAS_GLOB="${PIC_DIR}/tas_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-234912.nc"
PIC_PR_GLOB="${PIC_DIR}/pr_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-234912.nc"
PIC_TASMIN_GLOB="${PIC_DIR}/tasmin_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-234912.nc"
PIC_TASMAX_GLOB="${PIC_DIR}/tasmax_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-234912.nc"

remap_method () {          # pr conserves, temperatures bilinear
  case "$1" in
    pr) echo "remapcon" ;;
    *)  echo "remapbil" ;;
  esac
}

# build the piControl climatology ONCE per variable (whole run -> 12 months)
build_pic_clim () {
  local V="$1" GLOB="$2"
  local merged="${TMP}/${V}_pic_merged.nc"
  local clim="${TMP}/${V}_pic_clim.nc"
  if [ -f "$clim" ]; then echo "  pic clim $V already built"; return; fi
  echo "  [pic] mergetime ${V} ..."
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
  local anom="${OUT_DIR}/${V}_Amon_IPSL-CM6A-LR_${LABEL}_anomaly.nc"

  # grid check by SIZE (gridsize = number of cells). Identical native grids can
  # show float-dust label differences; if sizes match snap with setgrid (exact
  # coordinate copy, NO interpolation). Only a true size mismatch -> real remap.
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

  # pr ALSO gets the multiplicative ratio R = hosing / climatology (for ISIMIP
  # stressing: pr_AMOC = pr_ISIMIP * R). Temperature stays additive only.
  # ponytail: ymondiv blows up where climatological precip ~ 0 (deserts). Stored
  # faithful/unclipped; clamp R to a physical band or mask tiny-clim cells in the
  # ISIMIP-application step (matches EC/HadGEM builders).
  if [ "$V" = "pr" ]; then
    local ratio="${OUT_DIR}/${V}_Amon_IPSL-CM6A-LR_${LABEL}_ratio.nc"
    echo "  pr ratio = hosing / pic climatology (ymondiv)  (for ISIMIP stressing)"
    cdo -O setunit,'1' -ymondiv "$HOS" "$clim_use" "$ratio"
    echo "  -> $(basename "$ratio")"
  fi
}

# ---------------------------------------------------------------------------
# Build the four piControl climatologies once
# ---------------------------------------------------------------------------
echo "### Building piControl climatologies (1850-2349) ###"
build_pic_clim tas    "$PIC_TAS_GLOB"
build_pic_clim pr     "$PIC_PR_GLOB"
build_pic_clim tasmin "$PIC_TASMIN_GLOB"
build_pic_clim tasmax "$PIC_TASMAX_GLOB"

# ---------------------------------------------------------------------------
# Per experiment, per variable: subtract
#   tas/pr  -> native hosing files (DRS name _u03-hos_..._gr_)
#   tasmin/tasmax -> reconstructed files (_hos-u03-hos_reconstructed.nc)
# ---------------------------------------------------------------------------
for LABEL in hos-u03-hos; do
  echo
  echo "############## EXPERIMENT: $LABEL ##############"

  # native tas/pr carry the DRS token (u03-hos, not hos-u03-hos) + year span
  HOS_TAS=$(ls "${HOS_DIR}/tas_Amon_IPSL-CM6A-LR_u03-hos_"*.nc | head -1)
  HOS_PR=$(ls  "${HOS_DIR}/pr_Amon_IPSL-CM6A-LR_u03-hos_"*.nc  | head -1)
  # reconstructed tasmin/tasmax (label carries the hos- prefix)
  HOS_TASMIN="${HOS_DIR}/tasmin_Amon_IPSL-CM6A-LR_${LABEL}_reconstructed.nc"
  HOS_TASMAX="${HOS_DIR}/tasmax_Amon_IPSL-CM6A-LR_${LABEL}_reconstructed.nc"

  do_var tas    "$HOS_TAS"    "${TMP}/tas_pic_clim.nc"    "$LABEL"
  do_var pr     "$HOS_PR"     "${TMP}/pr_pic_clim.nc"     "$LABEL"
  do_var tasmin "$HOS_TASMIN" "${TMP}/tasmin_pic_clim.nc" "$LABEL"
  do_var tasmax "$HOS_TASMAX" "${TMP}/tasmax_pic_clim.nc" "$LABEL"
done

echo
echo "Done. Per-month, per-native-cell anomalies in:"
echo "  $OUT_DIR"
echo "hosing monthly  -  piControl 1850-2349 12-month climatology."
