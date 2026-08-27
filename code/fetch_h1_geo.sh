#!/usr/bin/env bash
# Phase 7 §12.2 — pull the H1 (GSE326743) Xenium analysis files from GEO.
#
# Why this is NOT fetch_xenium_bundle.sh: GSE310392 deposited one .tar.gz per
# sample with an ~11.7 GB .rds and a 4-11 GB morphology.ome.tif inside, so that
# script had to stream the archive and tear the pipe down mid-transfer.
# GSE326743 deposits each Xenium output as its OWN GEO supplementary file, so
# the same download discipline is just a URL filter — no tar streaming needed.
#
# Discipline is identical to §12.2: take counts, coordinates and segmentation
# polygons; skip transcripts.parquet, morphology*.ome.tif and nucleus
# boundaries. 28 files, ~530 MB, vs many GB if we took everything.
set -u
LIST="${1:-/workspace/logs/h1_sample_list.txt}"
DEST="${2:-/workspace/data/raw_h1}"
LOG=/workspace/logs/fetch_h1.log
mkdir -p "$DEST" /workspace/logs
echo "=== START $(date -u +%FT%TZ) ===" >>"$LOG"
fail=0
while IFS=$'\t' read -r gsm url; do
  [ -z "${gsm:-}" ] && continue
  fn=$(basename "$url")
  # sample dir keyed on the SPLNxx tag so section names read like the mouse arm
  samp=$(echo "$fn" | sed -E 's/^GSM[0-9]+_(SPLN[0-9]+)_.*/\1/')
  out="$DEST/$samp"; mkdir -p "$out"
  if [ -s "$out/$fn" ]; then echo "skip $fn (present)" >>"$LOG"; continue; fi
  if curl -sL --max-time 900 --retry 3 --retry-delay 5 "$url" -o "$out/$fn.part" \
     && [ -s "$out/$fn.part" ]; then
    mv "$out/$fn.part" "$out/$fn"; echo "OK   $samp/$fn ($(stat -c %s "$out/$fn") bytes)" >>"$LOG"
  else
    rm -f "$out/$fn.part"; echo "FAIL $samp/$fn" >>"$LOG"; fail=$((fail+1))
  fi
done < "$LIST"
echo "=== DONE $(date -u +%FT%TZ) failures=$fail ===" >>"$LOG"
exit $fail
