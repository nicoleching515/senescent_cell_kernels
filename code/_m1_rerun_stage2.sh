#!/bin/bash
# Phase 8 / task 8.7 -- the sender-definition (N7) axis and the second
# pre-registered Tier A variant, at 1,000 permutations (master plan Sec 24.3).
set -x
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=/workspace/logs
N7="tierA_p90,tierA_p99,cdkn1a_pos,senepy_p95,senepy_p99"

# --- wave 1: N7, corrected in-tissue nulls + published bounding-box nulls ----
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband \
    --calls $N7 --n-perm 1000 --n-jobs 30 > $L/m1_perm_c1_n7.log 2>&1 &
P1=$!
python3 -u run_phase3_nulls.py --stage perm --sections inband \
    --calls $N7 --n-perm 1000 --n-jobs 15 > $L/m1_perm_n7.log 2>&1 &
P2=$!
wait $P1 $P2
echo M1_WAVE1_DONE

# --- wave 2: the per-module Tier A sensitivity variant -----------------------
# the sender set is module-specific, so this fans out over (section x module):
# 6 sections x 7 modules = 42 jobs for one call.
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband \
    --calls tierApm_p95 --n-perm 1000 --n-jobs 42 --tag _pm \
    > $L/m1_perm_c1_pm.log 2>&1 &
P3=$!
python3 -u run_phase3_nulls.py --stage perm --sections inband \
    --calls tierApm_p95 --n-perm 1000 --n-jobs 12 --tag _pm \
    > $L/m1_perm_pm.log 2>&1 &
P4=$!
wait $P3 $P4
echo M1_WAVE2_DONE

python3 -u run_phase3_nulls.py --stage curves --sections inband \
    --calls tierA_p95 --n-jobs 6 > $L/m1_curves.log 2>&1
echo M1_STAGE2_DONE
