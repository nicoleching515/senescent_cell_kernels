#!/bin/bash
# Phase 8 / 8.7 -- the eight small Phase 3 scripts.  N8 has no internal
# parallelism and 200 random gene sets per section, so it is fanned out over
# sections here instead of looping serially.
set -x
cd /workspace/code
L=/workspace/logs
SECS="7259_liver_sbr_Male_26-U1 7260_liver_sbr_Male_26-U1 7001_liver_sham_Male_52-U1 7248_liver_sham_Male_26-U1 7352_liver_sham_Male_2-U1 7435_liver_sham_Male_10-U1 7239_liver_sbr_Male_52-U1 7361_liver_sbr_Male_2-U1 7450_liver_sbr_Male_10-U1"
for s in $SECS; do
  OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
    python3 -u run_phase3_n8.py --sections $s > $L/m1_n8_$s.log 2>&1 &
done
wait
cat $L/m1_n8_*.log > $L/m1_n8.log
echo M1_N8_DONE
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 -u run_phase3_strat.py            > $L/m1_strat.log 2>&1
python3 -u run_phase3_attribution.py      > $L/m1_attr.log 2>&1
python3 -u run_phase3_combined.py --sections inband > $L/m1_combined.log 2>&1
python3 -u run_phase3_poisson.py          > $L/m1_poisson.log 2>&1
python3 -u run_phase3_lamscale.py         > $L/m1_lamscale.log 2>&1
python3 -u _ripley.py                     > $L/m1_ripley.log 2>&1
python3 -u _correlogram.py                > $L/m1_correlogram.log 2>&1
echo M1_STAGE3_DONE
