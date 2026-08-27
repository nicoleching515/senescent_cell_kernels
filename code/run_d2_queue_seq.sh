#!/bin/bash
# Phase 8 task 8.5 / C7-D2, SEQUENTIAL runner.
#
# Replaces the 2-way parallel run_d2_queue.sh after it OOM-killed a section.  The reason
# is worth recording: `free` reports this host's 251 GB, but the CONTAINER's cgroup limit
# is /sys/fs/cgroup/memory.max = 57.7 GiB (61,999,996,928 bytes), and DeepScence holds several dense
# n_cells x 4845 float32 copies of a section at once.  Two large sections plus the DCA
# bridge exceed it.  One section at a time does not.  Another agent's M1 re-run shares
# this cgroup, so overshooting endangers its work too.
set -u
cd /workspace
L=/workspace/logs/deepscence_d2
LIM=$(cat /sys/fs/cgroup/memory.max)
wait_for_room () {   # block until the cgroup has at least 30 GB free
  while :; do
    cur=$(cat /sys/fs/cgroup/memory.current)
    free_gb=$(( (LIM - cur) / 1073741824 ))
    [ "$free_gb" -ge 30 ] && break
    sleep 30
  done
}
run1 () {  # $1 config, $2 section
  wait_for_room
  echo "--- $1 $2 $(date -u +%FT%TZ) (cgroup free $(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 )) GB)"
  OMP_NUM_THREADS=16 python3 /workspace/code/run_deepscence_denoise_probe.py --config "$1" "$2" \
    > $L/$1_$2.log 2>&1 && echo "OK $1 $2" || echo "FAILED $1 $2"
}
echo "=== SEQ START $(date -u +%FT%TZ) ==="
# wait out anything still running from the parallel attempt
while pgrep -f "run_deepscence_denoise_probe.py" > /dev/null; do sleep 30; done
for s in 7250_liver_sham_Male_26-U1 7260_liver_sbr_Male_26-U1 7352_liver_sham_Male_2-U1 \
         7361_liver_sbr_Male_2-U1 7001_liver_sham_Male_52-U1 7448_liver_sbr_Male_10-U1 \
         7435_liver_sham_Male_10-U1 7259_liver_sbr_Male_26-U1; do
  run1 mor "$s"
done
echo "--- mor done $(date -u +%FT%TZ) ---"
for s in 7239_liver_sbr_Male_52-U1 7259_liver_sbr_Male_26-U1; do run1 raw "$s"; done
echo "--- raw determinism controls done $(date -u +%FT%TZ) ---"
for s in 7259_liver_sbr_Male_26-U1 7352_liver_sham_Male_2-U1 7248_liver_sham_Male_26-U1; do run1 lib "$s"; done
echo "=== SEQ DONE $(date -u +%FT%TZ) ==="
