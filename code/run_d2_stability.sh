#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Phase 8 task 8.5 / C7-D2 -- seed-stability panel on a fixed 20,000-cell subsample.
#
# Purpose: put a noise floor under the denoise=True vs denoise=False comparison without
# needing a second full-section run.  Four runs on the SAME 20,000 cells (subsampling seed
# fixed at 12345, independent of --seed):
#   raw  seed 0 / seed 1   -> DeepScence's own run-to-run spread at denoise=False
#   dca  seed 0 / seed 1   -> the spread once DCA is in the loop (DeepScence passes its
#                             random_state straight into dca(), so one seed moves both)
# and raw-seed0 vs dca-seed0 gives the denoise effect on the identical cells.
#
# 20k cells is ~0.4 GB dense, so this fits alongside the M1 re-run that is currently
# holding ~30 GB of the 57.7 GiB cgroup, which the full-section jobs do not.
set -u
# The DCA arm runs in a SEPARATE CPython 3.8 venv built by code/setup_dca_env.sh, which
# takes DCA_ENV_ROOT as a required parameter.  These runners used to hardcode one
# session-scoped path under /tmp -- meaningless on any other machine and gone after a pod
# restart (AUDIT_REPRODUCIBILITY D3).  DCA_ENV_ROOT now wins; the original path is kept only
# as the last fallback, and a missing interpreter is fatal rather than silent.
SP="${DCA_ENV_ROOT:-/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt}"
export DCA_VENV_PYTHON="${DCA_VENV_PYTHON:-$SP/v38/bin/python}"
export DCA_BRIDGE_SCRATCH="${DCA_BRIDGE_SCRATCH:-$SP/scratch}"
if [ ! -x "$DCA_VENV_PYTHON" ]; then
  echo "FATAL: no DCA interpreter at $DCA_VENV_PYTHON" >&2
  echo "  build it:  DCA_ENV_ROOT=/some/scratch bash code/setup_dca_env.sh" >&2
  echo "  then:      export DCA_ENV_ROOT=/some/scratch   (or DCA_VENV_PYTHON directly)" >&2
  exit 1
fi
mkdir -p "$DCA_BRIDGE_SCRATCH"
export DCA_THREADS=12
export OMP_NUM_THREADS=8
cd /workspace
S=7239_liver_sbr_Male_52-U1
L=/workspace/logs/deepscence_d2
echo "=== STABILITY START $(date -u +%FT%TZ) ==="
for sd in 0 1; do
  python3 code/run_deepscence_denoise_probe.py --config raw --subsample 20000 --seed $sd $S \
    > $L/stab_raw_seed$sd.log 2>&1 && echo "OK raw seed$sd" || echo "FAILED raw seed$sd"
done
for sd in 0 1; do
  python3 code/run_deepscence_dca.py --subsample 20000 --seed $sd $S \
    > $L/stab_dca_seed$sd.log 2>&1 && echo "OK dca seed$sd" || echo "FAILED dca seed$sd"
done
echo "=== STABILITY DONE $(date -u +%FT%TZ) ==="
