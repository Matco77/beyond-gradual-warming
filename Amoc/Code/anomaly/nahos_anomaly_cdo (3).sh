#!/usr/bin/env bash
###############################################################################
# NAHosMIP hosing minus piControl anomaly  --  HadGEM3-GC31-LL and -MM
#
# Generalised, file-discovery version of the original per-protocol scripts.
# It now builds anomalies for BOTH hosing protocols
#       g01-hos   (weak hosing)
#       u03-hos   (strong hosing)            <-- added
# for BOTH HadGEM resolutions (LL, MM) and all four variables
#       tas  tasmin  tasmax  pr
#
# Method per (model x protocol x variable) -- identical to the original g01
# pipeline, so the g01 outputs are reproduced byte-for-meaning:
#   1. mergetime ALL hosing monthly chunk files (auto-discovered by glob)
#   2. mergetime the piControl 1850-1949 baseline chunk(s) (auto-selected)
#   3. monmean (harmless no-op on Amon monthly) then ymonmean  -> climatology
#   4. anomaly = hosing monthly - piControl monthly climatology (ymonsub)
#
# Why glob instead of hard-coded names:
#   the u03-hos runs do not necessarily span the same years as g01-hos, and
#   the year-ranges differ between LL and MM. Discovering the chunk files by
#   pattern makes the script correct for whatever ranges are actually present,
#   and it simply *skips with a warning* any (model,protocol) that has no data.
#
# Baseline note: the unperturbed reference is piControl 1850-1949 (the first
# 100-yr "parallel window"). It is protocol-independent, so g01-hos and u03-hos
# anomalies share one common baseline and stay directly comparable.
#
# Output names match what amoc_delta_effect_contrast.R expects:
#   anomaly_output/LL_anomaly/<var>_anomaly_<proto>_minus_piControl_1850-1949.nc
#   anomaly_output/MM_anomaly/<var>_anomaly_<proto>_minus_piControl_1850-1949.nc
#
# Requires CDO:  brew install cdo
# CDO streams from disk, so large inputs do NOT blow memory.
###############################################################################

set -euo pipefail

# ---- top-level configuration ------------------------------------------------
DATASETS="/Users/Bova/Library/CloudStorage/OneDrive-UniversitàCommercialeLuigiBocconi/1.Tesi/Amoc/datasets"

MODELS=(HadGEM3-GC31-LL HadGEM3-GC31-MM)
PROTOCOLS=(g01-hos u03-hos)          # set to (u03-hos) to (re)build only u03
VARS=(tas tasmin tasmax pr)

# piControl baseline window (the "parallel" 100 yr), as YYYYMM bounds.
BASE_START=185001
BASE_END=194912

# Local scratch for CDO intermediates + staging, kept OUT of OneDrive. The
# OneDrive File Provider briefly locks files mid-sync, so writing NetCDF straight
# into the CloudStorage tree intermittently fails with "Permission denied"; a
# local scratch sidesteps that and avoids sync churn (also faster I/O).
SCRATCH_BASE="${AMOC_SCRATCH:-${TMPDIR:-/tmp}/amoc_anom}"

# full model name -> the anomaly sub-folder the R script keys on
anomaly_folder () {
  case "$1" in
    HadGEM3-GC31-LL) echo "LL_anomaly" ;;
    HadGEM3-GC31-MM) echo "MM_anomaly" ;;
    *)               echo "${1}_anomaly" ;;
  esac
}

# first 6 digits (YYYYMM) of a piControl/hosing date token like
#   185001-194912     (monthly Amon)   -> 185001 / 194912
#   18500101-19491230 (daily)          -> 185001 / 194912
yyyymm () { local s="$1"; echo "${s:0:6}"; }

# Publish a freshly built file into an output dir that may live on OneDrive.
# OneDrive can lock a file/dir mid-sync -> retry with backoff; if it still will
# not take, keep the local copy and carry on (never abort the whole batch).
publish () {
  local src="$1" dest="$2" k
  for k in 1 2 3 4 5; do
    if [[ -e "$dest" ]]; then chmod u+w "$dest" 2>/dev/null || true; fi
    if mv -f "$src" "$dest" 2>/dev/null; then
      echo "  -> $(basename "$dest")"
      return 0
    fi
    echo "  [retry $k] output busy (OneDrive sync?), waiting ${k}s..." >&2
    sleep "$k"
  done
  echo "  [warn] could not publish into $(dirname "$dest") (OneDrive lock/permission?)." >&2
  echo "         kept local copy: ${src}" >&2
  return 1
}

# ---- one (model, protocol, variable) -> one anomaly file --------------------
do_var () {
  local MODEL="$1" PROTO="$2" V="$3"
  local HOS_DIR="${DATASETS}/NAHosMIP/${MODEL}"
  local PIC_DIR="${DATASETS}/CMIP6_piControl/${MODEL}"
  local OUT_DIR="${DATASETS}/anomaly_output/$(anomaly_folder "$MODEL")"
  local TMP="${SCRATCH_BASE}/$(anomaly_folder "$MODEL")/${PROTO}"   # local, not OneDrive
  mkdir -p "$OUT_DIR" "$TMP"

  # 1) discover the hosing chunk files for this model/protocol/variable
  local hos_files=() f
  for f in "${HOS_DIR}/${V}_Amon_${MODEL}_${PROTO}_"*.nc; do
    [[ -e "$f" ]] || continue
    hos_files+=("$f")
  done
  if [[ ${#hos_files[@]} -eq 0 ]]; then
    echo "  [skip] ${MODEL} ${PROTO} ${V}: no hosing files in ${HOS_DIR}" >&2
    return 0
  fi

  # 2) select the piControl chunk(s) overlapping BASE_START..BASE_END
  local pic_files=() tok start end
  for f in "${PIC_DIR}/${V}_Amon_${MODEL}_piControl_"*.nc; do
    [[ -e "$f" ]] || continue
    tok="${f##*_}"; tok="${tok%.nc}"            # e.g. 185001-194912
    start="$(yyyymm "${tok%%-*}")"              # chunk start YYYYMM
    end="$(yyyymm "${tok##*-}")"                # chunk end   YYYYMM
    if (( 10#$start <= 10#$BASE_END && 10#$end >= 10#$BASE_START )); then
      pic_files+=("$f")
    fi
  done
  if [[ ${#pic_files[@]} -eq 0 ]]; then
    echo "  [skip] ${MODEL} ${PROTO} ${V}: no piControl ${BASE_START}-${BASE_END} files in ${PIC_DIR}" >&2
    return 0
  fi

  local hos_m="${TMP}/${V}_${PROTO}_monthly.nc"
  local pic_merged="${TMP}/${V}_pic_merged.nc"
  local pic_m="${TMP}/${V}_pic_monthly.nc"
  local pic_clim="${TMP}/${V}_pic_clim.nc"
  local fname="${V}_anomaly_${PROTO}_minus_piControl_1850-1949.nc"
  local anom_tmp="${TMP}/${fname}"          # built locally first
  local anom="${OUT_DIR}/${fname}"          # then published to (OneDrive) OUT_DIR

  echo "=========== ${MODEL} ${PROTO} ${V}  (${#hos_files[@]} hosing, ${#pic_files[@]} piControl chunk(s)) ==========="

  echo "[1] merge hosing monthly"
  cdo -O mergetime "${hos_files[@]}" "$hos_m"

  echo "[2] merge piControl ${BASE_START}-${BASE_END} -> full baseline series"
  cdo -O mergetime "${pic_files[@]}" "$pic_merged"

  echo "[3] piControl monthly (Amon); monmean = harmless no-op pass"
  cdo -O monmean "$pic_merged" "$pic_m"

  echo "[4] monthly climatology (ymonmean over the 100-yr baseline)"
  cdo -O ymonmean "$pic_m" "$pic_clim"

  echo "[5] anomaly = hosing - climatology"
  cdo -O ymonsub "$hos_m" "$pic_clim" "$anom_tmp"

  echo "[6] publish -> ${OUT_DIR}"
  publish "$anom_tmp" "$anom" || true        # on give-up: local copy kept, batch continues

  # Precipitation ALSO gets a multiplicative ratio field R = hosing / climatology
  # (same 1850-1949 denominator as the additive anomaly). R is the quantity the
  # ISIMIP stressing step applies as pr_AMOC = pr_ISIMIP * R. Temperature stays
  # additive only (a ratio of interval-scale K is meaningless). UNITS: the
  # additive pr anomaly stays in native kg m-2 s-1 (CMIP6); the R diagnostics
  # convert ONCE at read (amoc_common.R read_europe_cube, x86400 -> mm/day)
  # before labelling plots/CSVs mm/day. The ratio itself is dimensionless: its
  # units attribute is set to "1" below.
  # ponytail: ymondiv blows up where climatological precip ~ 0 (global deserts, and a
  # few dry Mediterranean-summer cells: measured Europe max ~110x, ~0.3% of cell-months
  # > 5x). Stored FAITHFUL/unclipped here; clamp R to a physical band (e.g. [0.1, 10])
  # or mask tiny-climatology cells in the ISIMIP-application step, not in this field.
  if [[ "$V" == "pr" ]]; then
    local rname="pr_ratio_${PROTO}_over_piControl_1850-1949.nc"
    echo "[5b] pr ratio = hosing / climatology  (for ISIMIP stressing)"
    cdo -O setunit,'1' -ymondiv "$hos_m" "$pic_clim" "${TMP}/${rname}"
    echo "[6b] publish -> ${OUT_DIR}"
    publish "${TMP}/${rname}" "${OUT_DIR}/${rname}" || true
  fi

  rm -f "$pic_merged" "$pic_m" "$pic_clim" "$hos_m"
}

# ---- run: model x protocol x variable ---------------------------------------
for MODEL in "${MODELS[@]}"; do
  for PROTO in "${PROTOCOLS[@]}"; do
    for V in "${VARS[@]}"; do
      do_var "$MODEL" "$PROTO" "$V"
    done
  done
done

echo
echo "Done. Anomaly files under: ${DATASETS}/anomaly_output/{LL_anomaly,MM_anomaly}"
echo "Each is hosing monthly minus piControl 1850-1949 climatology, per protocol."
