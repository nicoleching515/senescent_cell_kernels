#!/usr/bin/env bash
# Phase 10 / PI decision: the H1 DeepScence score becomes a MULTI-SEED CONSENSUS,
# because the frozen single-seed configuration is not reproducible on this arm
# (full-section seed 0 vs 1: Pearson r 0.372, top-5 % Jaccard 0.211, against M1's
# floor of 0.9955 / 0.761 -- CS_PHASE9_H1_AUDIT.md s10.5 / H10).
#
# Seeds are the five ALREADY PRE-REGISTERED composition-matched seeds
# (PREREG_PHASE8.md s3.8): 20260901..05.  No new seed value is invented.
# Configuration is otherwise the frozen one: denoise=False, published CDKN1A anchor,
# >= 20 counts/cell, NATIVE human panel, full sections, no subsampling.
#
# SEED-MAJOR order: every stopping point leaves a COMPLETE 7-section panel at k seeds
# rather than a complete seed set on a few sections.
# At most 2 concurrent runs -- 3 OOM-kill each other in this 57.7 GB cgroup (H8).
set -u
MAXP=${MAXP:-2}
SEEDS=${SEEDS:-"20260901 20260902 20260903 20260904 20260905"}
SECS=${SECS:-"SPLN21 SPLN07 SPLN44 SPLN14 SPLN30 SPLN43 SPLN24"}   # ascending n_cells
LOG=/workspace/logs/phase10/deepscence_seeds.log
cd /workspace
mkdir -p /workspace/logs/phase10
for sd in $SEEDS; do
  for s in $SECS; do
    out="data/processed_h1/deepscence_h1_nodn_seed${sd}_${s}.csv"
    [ -s "$out" ] && continue
    while [ "$(pgrep -fc "python3 code/h1_deepscence_dca.py" || echo 0)" -ge "$MAXP" ]; do sleep 30; done
    anon=$(grep '^anon ' /sys/fs/cgroup/memory.stat | awk '{print int($2/1073741824)}')
    echo "=== $(date -u +%FT%TZ) launch seed=$sd $s (anon=${anon}GB) ===" >> "$LOG"
    OMP_NUM_THREADS=24 nohup python3 code/h1_deepscence_dca.py --denoise-false \
        --seed "$sd" "$s" >> "$LOG" 2>&1 &
    sleep 90
  done
done
while [ "$(pgrep -fc "python3 code/h1_deepscence_dca.py" || echo 0)" -gt 0 ]; do sleep 30; done
echo "=== ALL DONE $(date -u +%FT%TZ) ===" >> "$LOG"
