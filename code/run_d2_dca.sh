#!/bin/bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
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
export DCA_THREADS=24
export OMP_NUM_THREADS=16
cd /workspace
echo "=== DCA START $(date -u +%FT%TZ) ==="
for s in 7239_liver_sbr_Male_52-U1 7259_liver_sbr_Male_26-U1; do
  echo "--- $s $(date -u +%FT%TZ)"
  python3 /workspace/code/run_deepscence_dca.py "$s" > /workspace/logs/deepscence_d2/dca_$s.log 2>&1     && echo "OK dca $s" || echo "FAILED dca $s"
done
echo "=== DCA DONE $(date -u +%FT%TZ) ==="
