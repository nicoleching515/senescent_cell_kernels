#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 / 8.7 -- chain the remaining M1 re-run stages once the three
# already-running primary stages have written their outputs.
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=/workspace/logs
while [ "$(grep -c '^\[main\]'    $L/m1_main.log)"    -lt 99 ] \
   || [ "$(grep -c '^\[perm\]'    $L/m1_perm.log)"    -lt 6  ] \
   || [ "$(grep -c '^\[perm_c1\]' $L/m1_perm_c1.log)" -lt 6  ]; do sleep 20; done
echo "PRIMARY STAGES DONE $(date -u)"
bash _m1_rerun_stage2.sh >> $L/m1_chain.log 2>&1
echo "STAGE2 DONE $(date -u)"
bash _m1_rerun_stage3.sh >> $L/m1_chain.log 2>&1
echo "STAGE3 DONE $(date -u)"
bash _m1_rerun_stage4.sh >> $L/m1_chain.log 2>&1
echo "STAGE4 DONE $(date -u)"
