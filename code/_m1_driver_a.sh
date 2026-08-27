#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
cd /workspace/code
L=/workspace/logs
while pgrep -f "run_phase3_nulls.py --stage main" >/dev/null; do sleep 20; done
echo "MAIN DONE $(date -u)  $(stat -c '%y %s' /workspace/results/phase3/main_fits.csv)"
bash _m1_rerun_stage3.sh >> $L/m1_chain.log 2>&1
echo "STAGE3 DONE $(date -u)"
bash _m1_rerun_stage4.sh >> $L/m1_chain.log 2>&1
echo "STAGE4 DONE $(date -u)"
