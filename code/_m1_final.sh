#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 / 8.7 -- final stage: wait for every in-flight job, then summaries and
# the ONE figure pass.  Waits are on log markers and file mtimes, never on
# pgrep patterns (which match the supervising shells themselves).
cd /workspace/code
L=/workspace/logs
START=1787812466
newer () { [ -f "$1" ] && [ "$(stat -c %Y "$1")" -gt "$START" ]; }

# 1. perturbation nulls: N7 at 1000, then the per-module variant, then curves
until grep -q M1_STAGE2_DONE $L/m1_chain.log 2>/dev/null; do sleep 20; done
echo "STAGE2 COMPLETE $(date -u)"
# 2. the eight small Phase 3 scripts
until grep -q M1_STAGE3_DONE $L/m1_stage3.log 2>/dev/null; do sleep 20; done
echo "STAGE3 COMPLETE $(date -u)"
# 3. Phase 5
until grep -q M1_STAGE4_DONE $L/m1_stage4.log 2>/dev/null; do sleep 20; done
echo "STAGE4 COMPLETE $(date -u)"
# 4. A7 under the C6 sender sets
until newer /workspace/results/phase3/a7_control_probe_fits.csv; do sleep 20; done
echo "A7 COMPLETE $(date -u)"

bash /workspace/code/_m1_rerun_stage5.sh >> $L/m1_chain.log 2>&1
echo "STAGE5 COMPLETE $(date -u)"
