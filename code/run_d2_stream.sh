#!/bin/bash
# Phase 8 task 8.5 / C7-D2 -- ONE job at a time, all configurations in priority order.
#
# Two things forced this shape.
#  1. The container cgroup is 57.7 GiB (/sys/fs/cgroup/memory.max), not the 251 GB `free`
#     reports for the host, and DeepScence holds several dense n_cells x 4845 float32
#     copies of a section at once.
#  2. The M1 re-run in the next seat over runs ten run_phase3_n8.py workers in the SAME
#     cgroup, ~33 GB between them.  Two of my jobs plus theirs does not fit; one does.
# 16 GB is the guard, set empirically: an 83k-cell section was OOM-killed with 11 GB free,
# because DeepScence holds the input, `original`, the read_dataset copy, the normalize copy
# and a raw_counts layer -- five dense n_cells x 4845 float32 arrays -- at the same time.
# Sections above ~150k cells are dropped from this stream for the same reason -- 7248
# already covers the deep end of the depth range and was scored before the M1 re-run
# started.
#
# Order is by what the report cannot do without:
#   raw 7239   the determinism control, and the cheapest job here; without it no other
#              delta has a noise floor
#   dca 7259   a SECOND section at the published default; one section cannot separate
#              "denoising changes the caller" from "denoising changes it here"
#   mor 7259   the §6 estimator at the shallow end of the depth range
#   dca 7352   a third section at the published default, sham arm
#   raw 7259   determinism control on a score committed BEFORE two container rebuilds
#   lib 7259   full library-size equalisation, shallow end
#   lib 7352   full library-size equalisation, mid
#   mor 7352
set -u
cd /workspace
L=/workspace/logs/deepscence_d2
SP=/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt
export DCA_VENV_PYTHON=$SP/v38/bin/python
export DCA_BRIDGE_SCRATCH=$SP/scratch
export DCA_THREADS=16
export OMP_NUM_THREADS=12
LIM=$(cat /sys/fs/cgroup/memory.max)
room () {
  while :; do
    f=$(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 ))
    [ "$f" -ge 16 ] && break
    sleep 30
  done
}
run1 () {  # $1 config, $2 section
  room
  echo "--- $1 $2 $(date -u +%FT%TZ) (cgroup free $(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 )) GB)"
  if [ "$1" = dca ]; then
    python3 /workspace/code/run_deepscence_dca.py "$2" > $L/dca_$2.log 2>&1
  else
    python3 /workspace/code/run_deepscence_denoise_probe.py --config "$1" "$2" > $L/$1_$2.log 2>&1
  fi
  [ $? -eq 0 ] && echo "OK $1 $2" || echo "FAILED $1 $2"
}
echo "=== D2 STREAM START $(date -u +%FT%TZ) ==="
run1 raw 7239_liver_sbr_Male_52-U1
run1 dca 7259_liver_sbr_Male_26-U1
run1 mor 7259_liver_sbr_Male_26-U1
run1 raw 7259_liver_sbr_Male_26-U1
run1 dca 7352_liver_sham_Male_2-U1
run1 lib 7259_liver_sbr_Male_26-U1
run1 lib 7352_liver_sham_Male_2-U1
run1 mor 7352_liver_sham_Male_2-U1
echo "=== D2 STREAM DONE $(date -u +%FT%TZ) ==="
