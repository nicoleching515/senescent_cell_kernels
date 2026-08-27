#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 / D15 -- composition-matched rerun protocol at 5 seeds, M1 arm.
# Stage A: PRIMARY Tier A (A_SENDER_FINAL_strict) -- matched (comp, full) + adjusted.
# Stage B: per-module Tier A sensitivity sets (A_sender_for_*.txt) -- matched (comp) + adjusted.
# Merge:   one fit-level CSV and one per-seed + summary CSV.
set -e
cd /workspace/code
R=/workspace/results/phase3
python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierA_p95 \
    --variants comp,full,comp_adj,type_adj,typecomp_adj --out-tag _tierA \
    > /workspace/logs/compmatch_tierA.log 2>&1
python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierApm_p95 \
    --variants comp,comp_adj,type_adj,typecomp_adj --out-tag _tierApm \
    > /workspace/logs/compmatch_tierApm.log 2>&1
python3 -u run_phase8_compmatch.py \
    --merge $R/compmatch_fits_tierA.csv,$R/compmatch_fits_tierApm.csv \
    > /workspace/logs/compmatch_merge.log 2>&1
echo CHAIN_DONE
