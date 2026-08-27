#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
set -x
cd /workspace/code
export OMP_NUM_THREADS=1
until [ -f /workspace/results/phase3/main_fits.csv ]; do sleep 20; done
until [ -f /workspace/results/phase3/perm_nulls.csv ]; do sleep 20; done
python3 -u /workspace/code/run_phase3_strat.py > /workspace/logs/phase3_strat2.log 2>&1
python3 -u /workspace/code/run_phase3_combined.py --sections inband > /workspace/logs/phase3_combined2.log 2>&1
python3 -u /workspace/code/run_phase3_poisson.py > /workspace/logs/phase3_poisson.log 2>&1
python3 -u /workspace/code/run_phase3_lamscale.py > /workspace/logs/phase3_lamscale.log 2>&1
echo FINISH_STAGE_DONE
