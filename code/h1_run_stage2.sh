#!/usr/bin/env bash
# Phase 9 stage 2 — callers, A6 covariate and the Phase-3-format cache, per section.
# Sequential: each step loads one full section matrix.
set -u
cd /workspace
for s in "$@"; do
  echo "=== $(date -u +%FT%TZ) stage2 $s ==="
  [ -s "data/processed_h1/senders_h1_${s}.csv" ] || \
    python3 code/h1_callers.py "$s" > "logs/phase9/callers_${s}.log" 2>&1 || echo "FAIL callers $s"
  [ -s "data/processed_h1/anatomy_h1_${s}.csv" ] || \
    python3 code/h1_a6_compartments.py "$s" > "logs/phase9/a6_${s}.log" 2>&1 || echo "FAIL a6 $s"
  [ -s "data/processed_h1/cache3_h1/${s}.npz" ] || \
    python3 code/h1_prep_cache.py "$s" > "logs/phase9/cache_${s}.log" 2>&1 || echo "FAIL cache $s"
  echo "  mem: $(( $(cat /sys/fs/cgroup/memory.current)/1048576 )) MB"
  tail -1 "logs/phase9/cache_${s}.log" 2>/dev/null
done
echo "=== stage2 done $(date -u +%FT%TZ) ==="
