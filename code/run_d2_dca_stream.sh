#!/bin/bash
# Phase 8 task 8.5 / C7-D2 -- DeepScence at the PUBLISHED default (denoise=True) through
# the isolated DCA environment.  One section at a time, alongside the probe stream; see
# run_d2_probe_stream.sh for why the concurrency is capped.
set -u
SP=/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt
export DCA_VENV_PYTHON=$SP/v38/bin/python
export DCA_BRIDGE_SCRATCH=$SP/scratch
export DCA_THREADS=24
export OMP_NUM_THREADS=16
cd /workspace
LIM=$(cat /sys/fs/cgroup/memory.max)
echo "=== DCA STREAM START $(date -u +%FT%TZ) ==="
while pgrep -f "run_deepscence_dca.py" > /dev/null; do sleep 30; done
for s in 7352_liver_sham_Male_2-U1 7248_liver_sham_Male_26-U1; do
  while :; do
    free_gb=$(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 ))
    [ "$free_gb" -ge 28 ] && break
    sleep 30
  done
  echo "--- dca $s $(date -u +%FT%TZ) (cgroup free $(( (LIM - $(cat /sys/fs/cgroup/memory.current)) / 1073741824 )) GB)"
  python3 /workspace/code/run_deepscence_dca.py "$s" > /workspace/logs/deepscence_d2/dca_$s.log 2>&1 \
    && echo "OK dca $s" || echo "FAILED dca $s"
done
echo "=== DCA STREAM DONE $(date -u +%FT%TZ) ==="
