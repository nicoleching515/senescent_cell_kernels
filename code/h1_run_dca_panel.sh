#!/usr/bin/env bash
set -u
cd /workspace
export DCA_VENV_PYTHON=/tmp/dca_env/v38/bin/python
export DCA_BRIDGE_SCRATCH=/tmp/dca_scratch
export DCA_THREADS=24
export OMP_NUM_THREADS=24
# P-vii: three seeds on ONE fixed 20,000-cell subsample, exactly the M1 design
for s in 0 1 2; do
  python3 code/h1_deepscence_dca.py --subsample 20000 --seed $s SPLN21 || echo "FAIL seed $s"
done
# denoise=False companion on the SAME 20,000 cells, the seed-to-seed floor
echo "=== P-vi: full-section denoise=True on SPLN21 ==="
python3 code/h1_deepscence_dca.py SPLN21 || echo "FAIL SPLN21 full"
echo "=== dca h1 done $(date -u +%FT%TZ) ==="
