"""A `dca.api.dca` that satisfies DeepScence's import but runs DCA out-of-process.

DeepScence/api.py does `from dca.api import dca` at module import time and, when
denoise=True, calls `dca(adata, random_state=random_state)` and relies on it mutating
adata.X in place.  This shim keeps that contract while executing the real DCA 0.3.4 in the
isolated Python 3.8 / TensorFlow 2.4.4 venv (see code/dca_denoise_worker.py).  It is in a
SEPARATE shim directory from code/_shims/dca, which still raises, so nothing that does not
explicitly ask for the bridge can pick it up.

Set DCA_VENV_PYTHON to the venv interpreter; DCA_BRIDGE_SCRATCH to a directory with room
for one dense n_cells x n_genes float32 array.
"""
import os, sys, json, subprocess, tempfile
import numpy as np
import scipy.sparse as sp

# DeepScence/api.py takes `original = adata.copy()` BEFORE calling dca(), so anything the
# bridge writes into adata.uns is dropped from the returned object.  Park the worker's
# provenance here instead; run_deepscence_dca.py reads it after the call.
LAST_RUN = {}

WORKER = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'dca_denoise_worker.py')


def dca(adata, random_state=0, **kwargs):
    py = os.environ.get('DCA_VENV_PYTHON')
    if not py or not os.path.exists(py):
        raise RuntimeError('DCA_VENV_PYTHON is not set to an existing interpreter; '
                           'the DCA bridge cannot run.')
    scratch = os.environ.get('DCA_BRIDGE_SCRATCH', tempfile.gettempdir())
    os.makedirs(scratch, exist_ok=True)
    d = tempfile.mkdtemp(dir=scratch, prefix='dca_')
    cin, cout, cmeta = (os.path.join(d, f) for f in ('counts.npz', 'denoised.npy', 'meta.json'))
    X = adata.X
    sp.save_npz(cin, sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr())
    cmd = [py, WORKER, cin, cout, str(random_state),
           os.environ.get('DCA_THREADS', os.environ.get('OMP_NUM_THREADS', '16')), cmeta]
    print('[dca-bridge] ' + ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    Xd = np.load(cout)
    assert Xd.shape == adata.shape, 'DCA returned %s, expected %s' % (Xd.shape, adata.shape)
    adata.X = Xd
    meta = json.load(open(cmeta))
    adata.uns['dca_bridge'] = meta
    LAST_RUN.clear(); LAST_RUN.update(meta)
    for f in (cin, cout):
        os.remove(f)
    print('[dca-bridge] denoised in %s min under TF %s / keras %s'
          % (meta['minutes'], meta['tensorflow'], meta['keras']), flush=True)
    return None
