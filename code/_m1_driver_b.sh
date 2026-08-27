#!/bin/bash
cd /workspace/code
L=/workspace/logs
while pgrep -f "run_phase3_nulls.py --stage perm --sections inband --calls tierA_p95" >/dev/null \
   || pgrep -f "run_phase3_nulls.py --stage perm_c1 --sections inband --calls tierA_p95" >/dev/null; do sleep 20; done
echo "PRIMARY PERM STAGES DONE $(date -u)"
stat -c '%y %s %n' /workspace/results/phase3/perm_nulls.csv /workspace/results/phase3/perm_nulls_c1.csv
bash _m1_rerun_stage2.sh >> $L/m1_chain.log 2>&1
echo "STAGE2 DONE $(date -u)"
