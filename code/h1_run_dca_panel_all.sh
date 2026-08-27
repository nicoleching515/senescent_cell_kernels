#!/usr/bin/env bash
# Phase 10 — extend PREREG §8 prediction P-vi from 1 section to all 7.
#
# P-vi's falsifier is "Delta rho <= 0 in >= 5 of 7 sections" and Phase 9 ran 1 of 7, so the
# prediction was contradicted-but-not-falsified.  A full-section denoise=True run costs a DCA
# autoencoder pass over 200-400k cells on CPU (TensorFlow 2.4 cannot use this box's CUDA-12
# GPU) and is not affordable beside the rest of Phase 10.  This runs the SAME fixed
# 20,000-cell panel design the mouse D2 study used (subsampling seed 12345, independent of
# --seed, so every configuration sees the identical cells) on the six sections Phase 9 did
# not cover, with the denoise=False companion on the identical cells.
#
# DECLARED: this evaluates P-vi on the 20,000-cell panel at 7 of 7 sections, and on the FULL
# section at 1 of 7 (SPLN21, from Phase 9).  On SPLN21 the two agree in sign and roughly in
# magnitude (Delta rho -0.2104 full, -0.3669 panel), which is why the panel is used, but it
# is an approximation and is reported as one.
set -u
cd /workspace
export DCA_VENV_PYTHON=/tmp/dca_env/v38/bin/python
export DCA_BRIDGE_SCRATCH=/tmp/dca_scratch
export DCA_THREADS=12
export OMP_NUM_THREADS=12
L=/workspace/logs/phase10/dca_panel_all.log
for s in SPLN07 SPLN14 SPLN24 SPLN30 SPLN43 SPLN44; do
  echo "=== $(date -u +%FT%TZ) $s denoise=True panel ===" >> $L
  python3 -u code/h1_deepscence_dca.py --subsample 20000 --seed 0 "$s" >> $L 2>&1 \
      || echo "FAIL dca $s" >> $L
  echo "=== $(date -u +%FT%TZ) $s denoise=False panel ===" >> $L
  python3 -u code/h1_deepscence_dca.py --subsample 20000 --seed 0 --denoise-false "$s" >> $L 2>&1 \
      || echo "FAIL nodn $s" >> $L
done
echo "=== DCA PANEL ALL DONE $(date -u +%FT%TZ) ===" >> $L
