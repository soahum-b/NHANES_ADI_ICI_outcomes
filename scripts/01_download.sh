#!/usr/bin/env bash
# Reproducible download of NHANES 2011-2018 + USDA FPED per-person day-1 files.
# CDC public XPT files:  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{startyr}/DataFiles/{FILE}_{suf}.xpt
# USDA FPED SFX archives: https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_{yy}_sas.exe
set -euo pipefail
cd "$(dirname "$0")/.."
CDC="data/raw/cdc"; FPED="data/raw/fped"
mkdir -p "$CDC" "$FPED"
CDCBASE="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
USDABASE="https://www.ars.usda.gov/ARSUserFiles/80400530/apps"

# suffix:startyear:fpedyy
CYCLES="G:2011:1112 H:2013:1314 I:2015:1516 J:2017:1718"
# CDC file stems (per participant): demographics, day1+day2 total nutrients, body measures, smoking
CDCFILES="DEMO DR1TOT DR2TOT BMX SMQ DR1IFF"   # DR1IFF = individual foods (for hPDI/uPDI)

get() { # url dest
  if [ -s "$2" ]; then echo "  cached  $2"; return; fi
  curl -sS -L --retry 3 --max-time 300 -o "$2" "$1"
  echo "  got     $2  ($(wc -c < "$2") bytes)"
}

echo "=== CDC XPT files ==="
for c in $CYCLES; do
  suf=${c%%:*}; rest=${c#*:}; yr=${rest%%:*}
  for f in $CDCFILES; do
    get "$CDCBASE/$yr/DataFiles/${f}_${suf}.xpt" "$CDC/${f}_${suf}.xpt"
  done
done

echo "=== USDA FPED per-person day-1 (SFX .exe -> .sas7bdat) ==="
for c in $CYCLES; do
  suf=${c%%:*}; rest=${c#*:}; yy=${rest#*:}
  exe="$FPED/FPED_DR1TOT_${yy}_sas.exe"
  get "$USDABASE/FPED_DR1TOT_${yy}_sas.exe" "$exe"
  # extract the .sas7bdat payload from the self-extracting archive
  ( cd "$FPED" && bsdtar -xf "FPED_DR1TOT_${yy}_sas.exe" )
  echo "  extracted fped_dr1tot_${yy}.sas7bdat"
done

echo "=== USDA FPED food-level day-1 (FPED_DR1IFF, for hPDI/uPDI) ==="
for c in $CYCLES; do
  suf=${c%%:*}; rest=${c#*:}; yy=${rest#*:}
  exe="$FPED/FPED_DR1IFF_${yy}_sas.exe"
  get "$USDABASE/FPED_DR1IFF_${yy}_sas.exe" "$exe"
  ( cd "$FPED" && bsdtar -xf "FPED_DR1IFF_${yy}_sas.exe" )
  echo "  extracted fped_dr1iff_${yy}.sas7bdat"
done
echo "=== DONE ==="
