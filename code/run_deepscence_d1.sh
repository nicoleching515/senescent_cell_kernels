#!/usr/bin/env bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 7 C7/D1: DeepScence on the nine M1 sections that have no scores.
# 5 concurrent x 20 threads on a 112-core box; each section densifies to
# <10 GB so peak RSS stays far under the 251 GB ceiling.
set -u
cd /workspace
SECTIONS="7001_liver_sham_Male_52-U1 7239_liver_sbr_Male_52-U1 7248_liver_sham_Male_26-U1
7260_liver_sbr_Male_26-U1 7352_liver_sham_Male_2-U1 7361_liver_sbr_Male_2-U1
7435_liver_sham_Male_10-U1 7448_liver_sbr_Male_10-U1 7450_liver_sbr_Male_10-U1"
mkdir -p logs/deepscence_d1
echo "=== START $(date -u +%FT%TZ) ==="
for s in $SECTIONS; do echo "$s"; done | xargs -n1 -P5 -I{} \
  env OMP_NUM_THREADS=20 MKL_NUM_THREADS=20 \
  bash -c 'python3 /workspace/code/run_deepscence_all.py "$1" > /workspace/logs/deepscence_d1/"$1".log 2>&1 && echo "OK $1" || echo "FAILED $1"' _ {}
echo "=== ALL DONE $(date -u +%FT%TZ) ==="
ls -l /workspace/data/processed/deepscence_*.csv
