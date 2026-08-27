#!/usr/bin/env bash
# A3 fallback screen (contingency for reports/A3_FALLBACK_SCREEN.md).
# Panel-first §12.1 step-2 screen for the two runner-up candidates of
# reports/PHASE7_H1_SCREEN.md. Takes ONLY cell_feature_matrix.h5 (panel
# verification) plus small metadata files (donor counting). Nothing else:
# no transcripts, no morphology, no boundaries.
# Destination is data/raw_h1_candidates/ -- NOT data/raw/ or data/raw_h1/.
set -u
DEST="${1:-/workspace/data/raw_h1_candidates}"
LOG=/workspace/logs/fetch_h1_candidates.log
mkdir -p "$DEST" /workspace/logs
B=https://ftp.ncbi.nlm.nih.gov/geo
echo "=== START $(date -u +%FT%TZ) ===" >>"$LOG"
get () { # $1 = url, $2 = subdir
  local fn out; fn=$(basename "$1"); out="$DEST/$2"; mkdir -p "$out"
  if [ -s "$out/$fn" ]; then echo "skip $fn" >>"$LOG"; return 0; fi
  if curl -sS -L --max-time 1800 --retry 3 --retry-delay 5 "$1" -o "$out/$fn.part" && [ -s "$out/$fn.part" ]; then
    mv "$out/$fn.part" "$out/$fn"; echo "OK   $2/$fn ($(stat -c %s "$out/$fn") bytes)" >>"$LOG"
  else
    rm -f "$out/$fn.part"; echo "FAIL $2/$fn" >>"$LOG"; return 1
  fi
}
# --- GSE336890, human kidney biopsy (AIN/ATI/reference), 9 Xenium Regions ---
get $B/samples/GSM9844nnn/GSM9844156/suppl/GSM9844156_Region01_cell_feature_matrix.h5 GSE336890
get $B/samples/GSM9844nnn/GSM9844157/suppl/GSM9844157_Region02_cell_feature_matrix.h5 GSE336890
get $B/samples/GSM9844nnn/GSM9844163/suppl/GSM9844163_Region14_cell_feature_matrix.h5 GSE336890
get $B/samples/GSM9844nnn/GSM9844164/suppl/GSM9844164_Region15_cell_feature_matrix.h5 GSE336890
get $B/series/GSE336nnn/GSE336890/suppl/GSE336890_cells_stats_dir.tar.gz             GSE336890
# --- GSE335963 (SuperSeries) -> Xenium subseries GSE335962, bone marrow biopsy ---
get $B/samples/GSM9824nnn/GSM9824312/suppl/GSM9824312_CH02_cell_feature_matrix.h5    GSE335963
get $B/samples/GSM9824nnn/GSM9824313/suppl/GSM9824313_CH15_cell_feature_matrix.h5    GSE335963
get $B/samples/GSM9824nnn/GSM9824314/suppl/GSM9824314_NC03_cell_feature_matrix.h5    GSE335963
get $B/samples/GSM9824nnn/GSM9824316/suppl/GSM9824316_NC05_cell_feature_matrix.h5    GSE335963
get $B/samples/GSM9824nnn/GSM9824312/suppl/GSM9824312_CH02_Metadata.csv.gz           GSE335963
get $B/samples/GSM9824nnn/GSM9824314/suppl/GSM9824314_NC03_Metadata.csv.gz           GSE335963
echo "=== DONE $(date -u +%FT%TZ) ===" >>"$LOG"
