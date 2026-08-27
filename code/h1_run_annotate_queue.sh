#!/usr/bin/env bash
# Phase 9 — annotate the remaining H1 sections one at a time.
# Sequential by design: peak RSS is ~24 GB for a 196k-cell section and the cgroup ceiling is
# 57.7 GB, so two large sections at once is not safe.  Largest first so it fails fast.
set -u
cd /workspace
for s in "$@"; do
  if [ -s "data/processed_h1/celltypes_h1_${s}.csv" ]; then
    echo "skip $s (present)"; continue
  fi
  echo "=== $(date -u +%FT%TZ) annotate $s ==="
  python3 code/h1_annotate.py "$s" > "logs/phase9/annot_${s}.log" 2>&1 \
    && echo "OK $s" || echo "FAIL $s"
  echo "mem after $s: $(( $(cat /sys/fs/cgroup/memory.current)/1048576 )) MB"
done
echo "=== queue done $(date -u +%FT%TZ) ==="
