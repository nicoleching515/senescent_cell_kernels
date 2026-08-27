#!/usr/bin/env python3
"""DCA denoising worker -- RUNS ONLY INSIDE THE ISOLATED PYTHON 3.8 VENV.

Phase 7 addendum §6 path 1.  DCA 0.3.4 pins tensorflow>=2.0,<2.5 and keras>=2.4,<2.6;
those TensorFlow wheels exist only for cp36/cp37/cp38, and this project's pinned stack is
Python 3.11.  The two are therefore not co-installable, and DCA is run out-of-process in a
separate interpreter that the main environment never imports.

Handoff format is deliberately primitive -- scipy .npz in, numpy .npy out -- because the
h5ad on-disk format changed between anndata 0.7.8 (the newest anndata that installs beside
TF 2.4) and anndata 0.12.19 (the pinned version).  Passing raw arrays sidesteps that
entirely.

Usage: dca_denoise_worker.py <counts.npz> <out.npy> <seed> <threads> <meta.json>
"""
import sys, os, json, time
import numpy as np
import scipy.sparse as sp

counts_in, out_npy, seed, threads, meta_out = sys.argv[1:6]
seed = int(seed); threads = int(threads)
os.environ['PYTHONHASHSEED'] = '0'

import warnings; warnings.filterwarnings('ignore')
import anndata, scanpy as sc, tensorflow as tf, keras
from dca.api import dca

X = sp.load_npz(counts_in).tocsr()
A = anndata.AnnData(X.astype(np.float32))
A.var_names = ['g%d' % i for i in range(A.n_vars)]
A.obs_names = ['c%d' % i for i in range(A.n_obs)]
t = time.time()
# exactly the call DeepScence/api.py makes: dca(adata, random_state=random_state),
# every other argument left at the DCA default (mode='denoise', ae_type='nb-conddisp',
# hidden_size=(64,32,64), epochs=300, batch_size=32, early_stop=15, RMSprop).
dca(A, random_state=seed, threads=threads)
mins = (time.time() - t) / 60
Xd = A.X
if sp.issparse(Xd):
    Xd = Xd.toarray()
Xd = np.asarray(Xd, dtype=np.float32)
np.save(out_npy, Xd)
json.dump(dict(tensorflow=tf.__version__, keras=keras.__version__,
               anndata=anndata.__version__, scanpy=sc.__version__,
               numpy=np.__version__, python=sys.version.split()[0],
               n_cells=int(A.n_obs), n_genes=int(A.n_vars),
               minutes=round(mins, 2), threads=threads, seed=seed,
               out_min=float(Xd.min()), out_max=float(Xd.max()),
               out_mean=float(Xd.mean())), open(meta_out, 'w'), indent=1)
print('DCA denoise done: %d x %d in %.1f min' % (A.n_obs, A.n_vars, mins), flush=True)
