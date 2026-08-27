#!/bin/bash
set -u
export DCA_VENV_PYTHON=/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt/v38/bin/python
export DCA_BRIDGE_SCRATCH=/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt/scratch
export DCA_THREADS=24
export OMP_NUM_THREADS=16
cd /workspace
echo "=== DCA START $(date -u +%FT%TZ) ==="
for s in 7239_liver_sbr_Male_52-U1 7259_liver_sbr_Male_26-U1; do
  echo "--- $s $(date -u +%FT%TZ)"
  python3 /workspace/code/run_deepscence_dca.py "$s" > /workspace/logs/deepscence_d2/dca_$s.log 2>&1     && echo "OK dca $s" || echo "FAILED dca $s"
done
echo "=== DCA DONE $(date -u +%FT%TZ) ==="
