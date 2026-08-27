#!/bin/bash
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=/workspace/logs
while pgrep -f "run_phase3_nulls.py --stage main" >/dev/null \
   || pgrep -f "run_phase3_nulls.py --stage perm --" >/dev/null \
   || pgrep -f "run_phase3_nulls.py --stage perm_c1 --" >/dev/null; do sleep 20; done
echo "PRIMARY STAGES DONE $(date -u)"
for f in main_fits.csv perm_nulls.csv perm_nulls_c1.csv; do
  echo "  $f $(stat -c '%y %s' /workspace/results/phase3/$f)"
done
bash _m1_rerun_stage2.sh >> $L/m1_chain.log 2>&1
echo "STAGE2 DONE $(date -u)"
bash _m1_rerun_stage3.sh >> $L/m1_chain.log 2>&1
echo "STAGE3 DONE $(date -u)"
bash _m1_rerun_stage4.sh >> $L/m1_chain.log 2>&1
echo "STAGE4 DONE $(date -u)"
