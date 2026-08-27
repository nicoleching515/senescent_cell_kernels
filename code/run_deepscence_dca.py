#!/usr/bin/env python3
"""Phase 8, task 8.5 / C7-D2 -- DeepScence at its PUBLISHED default, denoise=True.

§6 path 1.  The denoising step runs as real DCA 0.3.4 under TensorFlow 2.4.4 in an
isolated Python 3.8 venv, reached through code/_shims_dca_bridge (see that file).  The
main pinned Python 3.11 stack is untouched: it never imports TensorFlow, and the two
interpreters exchange nothing but a scipy .npz of counts and a numpy .npy of denoised
values.

Everything except `denoise` is identical to run_deepscence_all.py, so a difference
between deepscence_<section>.csv and deepscence_dca_<section>.csv is attributable to the
denoising step alone.

Env:
  DCA_VENV_PYTHON    interpreter of the isolated venv (required)
  DCA_BRIDGE_SCRATCH scratch dir with room for one dense n_cells x n_genes float32 array
  DCA_THREADS        threads for TensorFlow

Usage: run_deepscence_dca.py [--subsample N] <section_dir_name> [...]
"""
import sys, os, time, json, argparse
import numpy as np, pandas as pd

# The bridge must be imported BEFORE anything imports DeepScence, because
# DeepScence/api.py does `from dca.api import dca` at module scope and code/_shims/dca
# (which raises by design) is on sys.path for every other entry point in this project.
# Importing it here puts the bridge in sys.modules first; the assert makes the ordering a
# hard failure rather than a silent fallback to the raising stub.
sys.path.insert(0, '/workspace/code/_shims_dca_bridge')
import dca.api as _dca_api                                          # noqa: E402
assert 'dca_bridge' in open(_dca_api.__file__).read(), \
    'the raising stub, not the bridge, is on sys.path -- refusing to run'

sys.path.insert(0, '/workspace/code')
from run_deepscence_denoise_probe import load_section, PROC, META   # noqa: E402
import torch                                                        # noqa: E402
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "16")))
import DeepScence.api as api                                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--subsample', type=int, default=0,
                    help='score only N cells (seeded); for timing probes only')
    ap.add_argument('--seed', type=int, default=0,
                    help='random_state for BOTH DeepScence and DCA, exactly as '
                         'DeepScence/api.py wires it (it passes its own random_state '
                         'straight into dca()).  0 = the published/committed setting.')
    ap.add_argument('sections', nargs='+')
    a = ap.parse_args()
    for section in a.sections:
        suffix = ('dca' if not a.subsample else 'dca_sub%d' % a.subsample) \
                 + ('' if a.seed == 0 else '_seed%d' % a.seed)
        out = PROC + 'deepscence_%s_%s.csv' % (suffix, section)
        if os.path.exists(out):
            print('%s [%s]: already present, skipping' % (section, suffix), flush=True)
            continue
        A = load_section(section)
        if a.subsample and a.subsample < A.n_obs:
            # fixed subsampling seed, independent of --seed, so every seed sees the SAME cells
            idx = np.sort(np.random.default_rng(12345).choice(A.n_obs, a.subsample, replace=False))
            A = A[idx].copy()
        mapped_counts = np.asarray(A.X.sum(1)).ravel().astype(np.int64)
        ids_kept = A.obs_names.to_numpy()
        print('%s [%s]: %d cells x %d ortholog-mapped genes'
              % (section, suffix, A.n_obs, A.n_vars), flush=True)
        t = time.time()
        res = api.DeepScence(A, denoise=True, verbose=False, random_state=a.seed)
        mins = (time.time() - t) / 60
        pd.DataFrame({'cell_id': ids_kept,
                      'deepscence_score': np.round(res.obs['ds'].to_numpy(), 5),
                      'mapped_counts': mapped_counts}).to_csv(out, index=False)
        log = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in dict(res.uns['log']).items()}
        json.dump(dict(section=section, config=suffix, n_cells=int(len(ids_kept)),
                       n_genes=int(res.n_vars), denoise=True, random_state=a.seed,
                       anchor='published_CDKN1A', panel='ortholog_mapped_MGI_1to1',
                       min_counts_per_cell=20, direction_log=log,
                       dca_bridge=res.uns.get('dca_bridge') or dict(_dca_api.LAST_RUN), minutes=round(mins, 2)),
                  open(META + 'runmeta_%s_%s.json' % (suffix, section), 'w'),
                  indent=1, default=str)
        print('wrote %s (%d cells, %.1f min)' % (os.path.basename(out), len(ids_kept), mins),
              flush=True)


if __name__ == '__main__':
    main()
