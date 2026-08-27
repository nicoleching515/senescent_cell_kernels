#!/usr/bin/env bash
# Phase 10 — the Phase-5 producers on H1, run one stage at a time so they do not
# oversubscribe the box against the perturbation-null stages already running.
set -u
cd /workspace
L=/workspace/logs/phase10
for spec in "kernels section tierAmg_p95" "kernels section tierA_p95" \
            "kernels proxdown tierAmg_p95" "super section tierAmg_p95" \
            "super section tierA_p95" "kernels heldout tierAmg_p95"; do
  set -- $spec
  echo "=== $(date -u +%FT%TZ) $1 $2 $3 ===" >> $L/phase5.log
  python3 -u code/h1_run_phase5.py --which "$1" --stage "$2" --call "$3" --n-jobs 4 \
      >> $L/phase5.log 2>&1 || echo "FAILED $1 $2 $3" >> $L/phase5.log
done
echo "=== PHASE5 CHAIN DONE $(date -u +%FT%TZ) ===" >> $L/phase5.log
