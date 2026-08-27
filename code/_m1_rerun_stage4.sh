#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 / task 8.7 -- Phase 5.
set -x
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=/workspace/logs
python3 -u run_phase5_kernels.py --stage section  --sections inband --n-jobs 6 > $L/m1_p5_kern_sec.log 2>&1
python3 -u run_phase5_kernels.py --stage heldout  --sections inband --n-jobs 6 > $L/m1_p5_kern_held.log 2>&1
python3 -u run_phase5_kernels.py --stage proxdown --sections inband --n-jobs 6 > $L/m1_p5_kern_prox.log 2>&1
python3 -u run_phase5_super.py --stage section --sections inband --n-jobs 6 > $L/m1_p5_super_sec.log 2>&1
python3 -u run_phase5_super.py --stage nulls   --sections inband --n-jobs 6 > $L/m1_p5_super_null.log 2>&1
python3 -u run_phase5_super.py --stage heldout --sections inband --n-jobs 6 > $L/m1_p5_super_held.log 2>&1
for c in cdkn1a_pos senepy_p95; do
  python3 -u run_phase5_super.py --stage section --sections inband --call $c --tag _$c --n-jobs 6 > $L/m1_p5_super_sec_$c.log 2>&1
  python3 -u run_phase5_super.py --stage nulls   --sections inband --call $c --tag _$c --n-jobs 6 > $L/m1_p5_super_null_$c.log 2>&1
  python3 -u run_phase5_super.py --stage heldout --sections inband --call $c --tag _$c --n-jobs 6 > $L/m1_p5_super_held_$c.log 2>&1
done
python3 -u run_phase5_wc.py --stage both --sections inband --n-jobs 6 > $L/m1_p5_wc.log 2>&1
python3 -u run_phase5_wc.py --stage crossfit --sections inband --frac 0.75 --tag _f75 --n-jobs 6 > $L/m1_p5_wc75.log 2>&1
python3 -u _se_ratio_phase5.py        > $L/m1_p5_se.log 2>&1
python3 -u _spline_window_check.py    > $L/m1_p5_spline.log 2>&1
echo M1_STAGE4_DONE
