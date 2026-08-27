#!/usr/bin/env bash
# Phase 10 — the five-seed DeepScence panel on H1 (PI decision D-A).
#
# Seeds are the five ALREADY PRE-REGISTERED composition-matched seeds (PREREG §3.8),
# 20260901..05.  Configuration is otherwise the frozen one: denoise=False, published
# CDKN1A anchor, >= 20 counts/cell, NATIVE human panel, full sections, no subsampling.
#
# SEED-MAJOR: seed k is completed across all 7 sections before seed k+1 starts, so any
# stopping point leaves a COMPLETE 7-section panel at k seeds.
#
# TWO GUARDS, both learned the hard way.
#  (1) At most MAXP=2 concurrent runs.  Three OOM-kill each other in this 57.7 GB cgroup
#      (Phase 9 deviation H8).
#  (2) A MEMORY-HEADROOM guard against `anon` in memory.stat, not `memory.current` (which
#      counts reclaimable page cache) and not `free`.  The first version of this script
#      had only guard (1), and when a run was OOM-killed within its 90 s settling window
#      the process count fell and the loop launched the next one straight into the same
#      wall: FIVE runs of seed 20260901 (SPLN44, SPLN14, SPLN30, SPLN43, SPLN24) were
#      launched at 90 s intervals between 20:57 and 21:03 UTC and all five were killed.
#      That is recorded in reports/CS_PHASE10_TWO_ARM.md as deviation T11.
#  (3) The outer `while` re-derives the outstanding list every pass, so an OOM-killed run
#      is RETRIED rather than silently skipped -- which is what let the failure above go
#      unnoticed for two hours.
set -u
MAXP=${MAXP:-2}
MIN_FREE_GB=${MIN_FREE_GB:-17}          # a full H1 section costs ~13 GB anonymous
SEEDS=${SEEDS:-"20260901 20260902 20260903 20260904 20260905"}
SECS=${SECS:-"SPLN21 SPLN07 SPLN44 SPLN14 SPLN30 SPLN43 SPLN24"}   # ascending n_cells
LOG=/workspace/logs/phase10/deepscence_seeds.log
cd /workspace
mkdir -p /workspace/logs/phase10

nrun() { pgrep -fc "code/h1_deepscence_dca.py --denoise-false" 2>/dev/null || echo 0; }
freegb() {
  a=$(grep '^anon ' /sys/fs/cgroup/memory.stat | awk '{print $2}')
  m=$(cat /sys/fs/cgroup/memory.max)
  echo $(( (m - a) / 1073741824 ))
}

for sd in $SEEDS; do
  while true; do
    todo=()
    for s in $SECS; do
      [ -s "data/processed_h1/deepscence_h1_nodn_seed${sd}_${s}.csv" ] || todo+=("$s")
    done
    [ ${#todo[@]} -eq 0 ] && { echo "=== seed $sd COMPLETE $(date -u +%FT%TZ) ===" >> "$LOG"; break; }
    launched=0
    for s in "${todo[@]}"; do
      pgrep -f "h1_deepscence_dca.py --denoise-false --seed ${sd} ${s}$" >/dev/null && continue
      while [ "$(nrun)" -ge "$MAXP" ] || [ "$(freegb)" -lt "$MIN_FREE_GB" ]; do sleep 30; done
      echo "=== $(date -u +%FT%TZ) launch seed=$sd $s (running=$(nrun) free=$(freegb)GB) ===" >> "$LOG"
      OMP_NUM_THREADS=24 nohup python3 code/h1_deepscence_dca.py --denoise-false \
          --seed "$sd" "$s" >> "$LOG" 2>&1 &
      launched=1
      sleep 120
    done
    # nothing could be launched this pass (all outstanding runs are already alive): wait
    [ $launched -eq 0 ] && sleep 60
  done
done
echo "=== ALL SEEDS DONE $(date -u +%FT%TZ) ===" >> "$LOG"
