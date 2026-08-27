#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 task 8.5 / C7-D2 job queue.  Two concurrent DeepScence runs, never more:
# DeepScence densifies to two n_cells x 4845 float32 tensors plus AnnData copies and
# five concurrent sections OOM-killed this box during the D1 run (logs/deepscence_d1_2026-08-27.log).
set -u
cd /workspace
L=/workspace/logs/deepscence_d2
echo "=== START $(date -u +%FT%TZ) ==="
run_all () {   # $1 = config, rest = sections
  cfg=$1; shift
  printf '%s\n' "$@" | xargs -P 2 -I{} bash -c \
    'OMP_NUM_THREADS=32 python3 /workspace/code/run_deepscence_denoise_probe.py --config '"$cfg"' "$1" \
       > '"$L"'/'"$cfg"'_"$1".log 2>&1 && echo "OK '"$cfg"' $1" || echo "FAILED '"$cfg"' $1"' _ {}
}
# interleaved big/small so the two largest sections are never resident at the same time
run_all mor 7248_liver_sham_Male_26-U1 7450_liver_sbr_Male_10-U1 \
            7250_liver_sham_Male_26-U1 7239_liver_sbr_Male_52-U1 \
            7260_liver_sbr_Male_26-U1 7352_liver_sham_Male_2-U1 \
            7361_liver_sbr_Male_2-U1 7001_liver_sham_Male_52-U1 \
            7448_liver_sbr_Male_10-U1 7435_liver_sham_Male_10-U1 \
            7259_liver_sbr_Male_26-U1
echo "--- mor done $(date -u +%FT%TZ) ---"
run_all raw 7239_liver_sbr_Male_52-U1 7259_liver_sbr_Male_26-U1
echo "--- raw controls done $(date -u +%FT%TZ) ---"
run_all lib 7248_liver_sham_Male_26-U1 7239_liver_sbr_Male_52-U1 \
            7250_liver_sham_Male_26-U1 7352_liver_sham_Male_2-U1 \
            7259_liver_sbr_Male_26-U1
echo "--- lib done $(date -u +%FT%TZ) ---"
run_all ds10 7248_liver_sham_Male_26-U1 7352_liver_sham_Male_2-U1 \
             7250_liver_sham_Male_26-U1 7259_liver_sbr_Male_26-U1
echo "=== ALL DONE $(date -u +%FT%TZ) ==="
