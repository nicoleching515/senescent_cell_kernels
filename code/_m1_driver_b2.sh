#!/bin/bash
cd /workspace/code
L=/workspace/logs
PID=$1
while [ -d /proc/$PID ]; do sleep 20; done
echo "PRIMARY perm_c1 DONE $(date -u)"
stat -c '%y %s %n' /workspace/results/phase3/perm_nulls.csv /workspace/results/phase3/perm_nulls_c1.csv
bash _m1_rerun_stage2.sh >> $L/m1_chain.log 2>&1
echo "STAGE2 DONE $(date -u)"
