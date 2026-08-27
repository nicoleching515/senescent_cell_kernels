#!/usr/bin/env bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Pull ONLY the three Xenium analysis files from a GEO sample tarball.
#
# Why not just `curl | tar -x`:
#   Member ordering differs across GSE310392. Some archives put the Xenium
#   bundle first (~130 MB in); others front-load an ~11.7 GB processed .rds.
#   Every archive ends with morphology.ome.tif (4-11 GB) + transcripts.parquet
#   (1.9 GB), which Section 17.2 says explicitly not to download.
#
#   So: stream, extract the three members, and tear the pipeline down the
#   moment cells.parquet stops growing. cells.parquet is the last of the three
#   in every ordering observed, so a stable cells.parquet means we have all
#   three. Bundle-first archives finish in seconds; rds-first in minutes.
#
#   The pipeline runs under setsid in its own process group so we can kill the
#   whole group by PGID. Do NOT use `pkill -f "$GSM"` here: the accession
#   appears in this script's own command line, so it kills the caller too.
#
# Usage: fetch_xenium_bundle.sh <GSM> <TARBALL_BASENAME> <DEST_ROOT>

set -u
GSM="$1"; FN="$2"; DEST="$3"
PREFIX="${GSM%???}nnn"          # GEO masks the last 3 digits
URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/${PREFIX}/${GSM}/suppl/${FN}"
SAMPLE="${FN#${GSM}_}"; SAMPLE="${SAMPLE%.tar.gz}"
OUT="$DEST/$SAMPLE"
LOG="/workspace/logs/fetch_${GSM}.log"
mkdir -p /workspace/logs

WANT=(cell_feature_matrix.h5 cells.parquet cell_boundaries.parquet)
have_all() { for w in "${WANT[@]}"; do [ -s "$OUT/$w" ] || return 1; done; return 0; }

if have_all; then echo "$(date -u +%FT%TZ) $GSM already present, skipping" >>"$LOG"; exit 0; fi

mkdir -p "$OUT"
echo "$(date -u +%FT%TZ) START $GSM ($SAMPLE)" >>"$LOG"

setsid bash -c "curl -sL '$URL' \
  | gzip -dc 2>/dev/null \
  | head -c 31000000000 \
  | tar -x -C '$DEST' --wildcards \
      '*/cell_feature_matrix.h5' '*/cells.parquet' '*/cell_boundaries.parquet'" \
  >>"$LOG" 2>&1 &
PIPE=$!
PGID=$(ps -o pgid= -p "$PIPE" 2>/dev/null | tr -d ' ')

		# Watch the COMBINED size of all three targets, not just cells.parquet.
# Member ordering is not guaranteed identical across samples, so keying the
# teardown on one file risks either cutting early or streaming to the ceiling.
prev=-1; stable=0
while kill -0 "$PIPE" 2>/dev/null; do
  sleep 5
  cur=0
  for w in "${WANT[@]}"; do
    sz=$(stat -c %s "$OUT/$w" 2>/dev/null || echo 0); cur=$((cur+sz))
  done
  if [ "$cur" -gt 0 ] && [ "$cur" -eq "$prev" ]; then
    stable=$((stable+1))
    if [ "$stable" -ge 3 ] && have_all; then
      echo "$(date -u +%FT%TZ) $GSM bundle complete, tearing down stream early" >>"$LOG"
      # Group-kill ONLY if setsid actually isolated the pipeline into its own
      # process group. Without this guard a failed setsid means PGID == our own
      # group and `kill -TERM -$PGID` takes down this script AND the calling
      # runner — which is exactly what happened on GSM9295276.
      MYPGID=$(ps -o pgid= -p $$ 2>/dev/null | tr -d ' ')
      if [ -n "$PGID" ] && [ "$PGID" != "$MYPGID" ]; then
        kill -TERM -"$PGID" 2>/dev/null
      else
        pkill -P "$PIPE" 2>/dev/null; kill -TERM "$PIPE" 2>/dev/null
      fi
      break
    fi
  else
    stable=0
  fi
  prev=$cur
done
wait "$PIPE" 2>/dev/null

if have_all; then
  echo "$(date -u +%FT%TZ) OK $GSM" >>"$LOG"; ls -l "$OUT" >>"$LOG"; exit 0
else
  echo "$(date -u +%FT%TZ) FAILED $GSM (incomplete)" >>"$LOG"; ls -l "$OUT" >>"$LOG" 2>/dev/null; exit 1
fi
