#!/usr/bin/env bash
# Phase 10 - the five-seed DeepScence panel on H1 (PI decision D-A).
#
# Seeds are the five ALREADY PRE-REGISTERED composition-matched seeds (PREREG s3.8),
# 20260901..05.  Configuration is otherwise the frozen one: denoise=False, published
# CDKN1A anchor, >= 20 counts/cell, NATIVE human panel, full sections, no subsampling.
#
# SEED-MAJOR PRIORITY.  The work list is (seed, section) in seed-major order and the
# EARLIEST outstanding pair is always launched, so seed k finishes before seed k+1 starts
# -- except that when the current seed's last section is already running the free slot
# takes the next seed's first section rather than idling.  Completion ORDER changes
# nothing; only which (seed, section) pairs exist does, and that is recorded per row in
# results/phase10_h1/deepscence_consensus_coverage.csv.
#
# THREE GUARDS, all of them learned the hard way in this phase (deviation T11):
#  (1) MAXP=2 concurrent runs.  Three OOM-kill each other in this 57.7 GB cgroup (H8).
#  (2) A MEMORY-HEADROOM guard on `anon` in memory.stat -- not `memory.current`, which
#      counts reclaimable page cache, and not `free`.
#  (3) The outstanding list is re-derived every pass, so an OOM-killed run is RETRIED.
# Two bugs that defeated (1) and (3) and are fixed here:
#  * `pgrep -fc PAT || echo 0` emits "0\n0" when there is no match (pgrep -c prints 0 AND
#    exits 1), and `[ "0\n0" -ge 2 ]` is a shell ERROR rather than false, so the
#    concurrency guard silently stopped guarding.  `| wc -l` always yields one integer.
#  * pgrep matches /proc/pid/cmdline with NUL separators mapped to spaces and that string
#    can carry a trailing separator, so a `$` end-anchor does not match and the same
#    (seed, section) was launched twice.  Section names are mutually non-prefixing, so a
#    plain substring match is exact here.
set -u
MAXP=${MAXP:-2}
MIN_FREE_GB=${MIN_FREE_GB:-22}    # ~13 GB for a typical section, ~19 GB for SPLN24
SEEDS=${SEEDS:-"20260901 20260902 20260903 20260904 20260905"}
SECS=${SECS:-"SPLN21 SPLN07 SPLN44 SPLN14 SPLN30 SPLN43 SPLN24"}   # ascending n_cells
LOG=/workspace/logs/phase10/deepscence_seeds.log
cd /workspace
mkdir -p /workspace/logs/phase10

nrun() { pgrep -f "code/h1_deepscence_dca.py --denoise-false" 2>/dev/null | wc -l; }
freegb() {
  a=$(grep '^anon ' /sys/fs/cgroup/memory.stat | awk '{print $2}')
  m=$(cat /sys/fs/cgroup/memory.max)
  echo $(( (m - a) / 1073741824 ))
}

while true; do
  todo=()
  for sd in $SEEDS; do
    for s in $SECS; do
      [ -s "data/processed_h1/deepscence_h1_nodn_seed${sd}_${s}.csv" ] && continue
      pgrep -f "h1_deepscence_dca.py --denoise-false --seed ${sd} ${s}" >/dev/null && continue
      todo+=("${sd}:${s}")
    done
  done
  if [ ${#todo[@]} -eq 0 ]; then
    [ "$(nrun)" -eq 0 ] && break
    sleep 60
    continue
  fi
  while [ "$(nrun)" -ge "$MAXP" ] || [ "$(freegb)" -lt "$MIN_FREE_GB" ]; do sleep 30; done
  pair=${todo[0]}
  sd=${pair%%:*}
  s=${pair##*:}
  echo "=== $(date -u +%FT%TZ) launch seed=$sd $s (running=$(nrun) free=$(freegb)GB outstanding=${#todo[@]}) ===" >> "$LOG"
  OMP_NUM_THREADS=24 nohup python3 code/h1_deepscence_dca.py --denoise-false \
      --seed "$sd" "$s" >> "$LOG" 2>&1 &
  sleep 120
done
echo "=== ALL SEEDS DONE $(date -u +%FT%TZ) ===" >> "$LOG"
