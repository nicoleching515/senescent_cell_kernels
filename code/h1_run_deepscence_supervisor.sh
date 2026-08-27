#!/usr/bin/env bash
# Keep at most MAXP DeepScence processes alive until every H1 section has a score.
# 3 concurrent runs OOM-killed each other at this 57.7 GB ceiling; 2 is the maximum.
set -u
MAXP=2
cd /workspace
while true; do
  todo=()
  for s in SPLN44 SPLN43 SPLN14 SPLN07 SPLN21 SPLN24 SPLN30; do
    [ -s "data/processed_h1/deepscence_h1_${s}.csv" ] || todo+=("$s")
  done
  [ ${#todo[@]} -eq 0 ] && { echo "=== all sections scored $(date -u +%FT%TZ) ==="; break; }
  running=$(pgrep -fc "h1_deepscence.py" || true)
  for s in "${todo[@]}"; do
    pgrep -f "h1_deepscence.py $s" >/dev/null && continue
    [ "$running" -ge "$MAXP" ] && break
    anon=$(grep '^anon ' /sys/fs/cgroup/memory.stat | awk '{print int($2/1073741824)}')
    echo "=== $(date -u +%FT%TZ) launch $s (running=$running anon=${anon}GB) ==="
    OMP_NUM_THREADS=32 nohup python3 code/h1_deepscence.py "$s" \
      >> /workspace/logs/phase9/deepscence_h1.log 2>&1 &
    running=$((running+1))
    sleep 60
  done
  sleep 45
done
