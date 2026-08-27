#!/bin/bash
# Phase 8 task 8.5 / C7-D2 -- normalisation-probe stream (one section at a time).
#
# Scope was retargeted once §6 path 1 (real DCA) came up working: the point of the path-3
# probes is no longer to substitute for DCA but to say what a rescaling normalisation does
# compared with what DCA actually does.  Five depth-spanning sections for `mor` (including
# BOTH sections of the preserved two-section base, 7250 and 7259), three for `lib`, two
# determinism controls.  Compute freed by dropping the other six `mor` sections goes to
# running DCA on more sections instead.
#
# One section at a time: `free` reports the host's 251 GB but the CONTAINER cgroup
# (/sys/fs/cgroup/memory.max) is 57.7 GiB, and DeepScence holds several dense
# n_cells x 4845 float32 copies at once.  Another agent's M1 re-run shares this cgroup.
set -u
cd /workspace
L=/workspace/logs/deepscence_d2
LIM=$(cat /sys/fs/cgroup/memory.max)
wait_for_room () {
  while :; do
    free_gb=$(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 ))
    [ "$free_gb" -ge 28 ] && break
    sleep 30
  done
}
run1 () {
  wait_for_room
  echo "--- $1 $2 $(date -u +%FT%TZ) (cgroup free $(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 )) GB)"
  OMP_NUM_THREADS=16 python3 /workspace/code/run_deepscence_denoise_probe.py --config "$1" "$2" \
    > $L/$1_$2.log 2>&1 && echo "OK $1 $2" || echo "FAILED $1 $2"
}
echo "=== PROBE STREAM START $(date -u +%FT%TZ) ==="
while pgrep -f "run_deepscence_denoise_probe.py" > /dev/null; do sleep 30; done
run1 mor 7259_liver_sbr_Male_26-U1
run1 raw 7239_liver_sbr_Male_52-U1
run1 raw 7259_liver_sbr_Male_26-U1
run1 lib 7259_liver_sbr_Male_26-U1
run1 lib 7352_liver_sham_Male_2-U1
run1 lib 7248_liver_sham_Male_26-U1
echo "=== PROBE STREAM DONE $(date -u +%FT%TZ) ==="
