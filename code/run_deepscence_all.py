#!/usr/bin/env python3
"""DeepScence sender scoring across ALL M1 sections — Phase 7 §5 (C7/D1).

Why this file exists rather than an edit to run_deepscence.py:
  run_deepscence.py hardcodes SAMP={'sham':7250,'sbr':7259} and writes
  deepscence_{sham,sbr}.csv. Phase 7 §5 requires the ELEVEN-section values to be
  reported ALONGSIDE the two-section values, so the original two outputs must
  survive untouched. This script is parameterised by section directory name and
  writes deepscence_<section>.csv, leaving the 2-section base intact.

Settings are deliberately IDENTICAL to run_deepscence.py — denoise=False,
random_state=0, >=20 counts/cell, MGI 1:1 ortholog remap, published CDKN1A
anchor. D1 is a COVERAGE change only. Changing denoise (D2) or the anchor (D3)
at the same time would confound "did the headline move because of coverage?"
with "did it move because we changed the caller?". Those run separately.

Usage: run_deepscence_all.py <section_dir_name> [...]
"""
import sys, os, time, csv, numpy as np, pandas as pd, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix

# torch grabs every core by default; we run several sections concurrently, so
# each process is pinned to its share by the launcher via OMP_NUM_THREADS.
import torch
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))

sys.path.insert(0, '/workspace/code/_shims')
import DeepScence.api as api

RAW = '/workspace/data/raw/'
PROC = '/workspace/data/processed/'
ORTH = '/workspace/genesets/mouse_human_orthologs_MGI.csv'

orth = {r['mouse_symbol']: r['human_symbol'] for r in csv.DictReader(open(ORTH))}

for section in sys.argv[1:]:
    out = PROC + 'deepscence_%s.csv' % section
    if os.path.exists(out):
        print('%s: already present, skipping' % section, flush=True)
        continue
    d = RAW + section + '/'
    f = h5py.File(d + 'cell_feature_matrix.h5', 'r')
    ft = np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids = np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm = np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc = np.array([x.decode() for x in f['matrix/barcodes'][:]])
    keep = (ft == 'Gene Expression') & np.char.startswith(ids, 'ENSMUSG')
    M = csc_matrix((f['matrix/data'][:].astype(np.float32),
                    f['matrix/indices'][:].astype(np.int32),
                    f['matrix/indptr'][:].astype(np.int64)),
                   shape=tuple(f['matrix/shape'][:]))
    A = ad.AnnData(csr_matrix(M[np.where(keep)[0], :].T),
                   obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=nm[keep]))
    A = A[np.asarray(A.X.sum(1)).ravel() >= 20].copy()
    hs = np.array([orth.get(g, '') for g in A.var_names]); A = A[:, hs != ''].copy()
    A.var_names = hs[hs != '']; A.var_names_make_unique()
    ids_kept = A.obs_names.to_numpy()
    print('%s: %d cells x %d ortholog-mapped genes' % (section, A.n_obs, A.n_vars), flush=True)
    t = time.time()
    res = api.DeepScence(A, denoise=False, verbose=False, random_state=0)
    mins = (time.time() - t) / 60
    print('%s: DeepScence %.1f min' % (section, mins), flush=True)
    pd.DataFrame({'cell_id': ids_kept,
                  'deepscence_score': np.round(res.obs['ds'].to_numpy(), 5)}
                 ).to_csv(out, index=False)
    print('wrote %s (%d cells, %.1f min)' % (os.path.basename(out), A.n_obs, mins), flush=True)
