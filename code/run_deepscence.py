#!/usr/bin/env python3
"""DeepScence (Qu et al., Cell Genomics 2025) sender scoring, Sec 10 method 1.

TWO DEVIATIONS, both forced and both documented:
 1. DeepScence ships a HUMAN core gene set (CoreScence v2) and hardcodes 'CDKN1A' for direction
    fixing. Our data is mouse. We rename mouse symbols to their human 1:1 orthologs using the
    MGI HOM_MouseHumanSequence report -- the SAME resource the paper's own authors used
    (see their Zenodo notebook, which reads informatics.jax.org/.../HOM_MouseHumanSequence.rpt).
    Only unambiguous 1:1 pairs are used; 4,845 of 5,097 panel genes map.
 2. denoise=False. DeepScence's default denoising step requires the `dca` package, which depends
    on an obsolete TensorFlow stack and does not install here. Scores are therefore computed on
    undenoised counts. This is a real departure from the published default and must be stated.
"""
import sys, os, time, csv, numpy as np, pandas as pd, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix
sys.path.insert(0,'/workspace/code/_shims')
import DeepScence.api as api

RAW='/workspace/data/raw/'; PROC='/workspace/data/processed/'
SAMP={'sham':'7250_liver_sham_Male_26-U1','sbr':'7259_liver_sbr_Male_26-U1'}
orth={r['mouse_symbol']:r['human_symbol'] for r in
      csv.DictReader(open('/workspace/genesets/mouse_human_orthologs_MGI.csv'))}
for tag in sys.argv[1:]:
    d=RAW+SAMP[tag]+'/'
    f=h5py.File(d+'cell_feature_matrix.h5','r')
    ft=np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids=np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm=np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc=np.array([x.decode() for x in f['matrix/barcodes'][:]])
    keep=(ft=='Gene Expression')&np.char.startswith(ids,'ENSMUSG')
    M=csc_matrix((f['matrix/data'][:].astype(np.float32),f['matrix/indices'][:].astype(np.int32),
                  f['matrix/indptr'][:].astype(np.int64)),shape=tuple(f['matrix/shape'][:]))
    A=ad.AnnData(csr_matrix(M[np.where(keep)[0],:].T), obs=pd.DataFrame(index=bc),
                 var=pd.DataFrame(index=nm[keep]))
    A=A[np.asarray(A.X.sum(1)).ravel()>=20].copy()
    hs=np.array([orth.get(g,'') for g in A.var_names]); A=A[:,hs!=''].copy()
    A.var_names=hs[hs!='']; A.var_names_make_unique()
    ids_kept=A.obs_names.to_numpy()
    print('%s: %d cells x %d ortholog-mapped genes'%(tag,A.n_obs,A.n_vars), flush=True)
    t=time.time(); res=api.DeepScence(A, denoise=False, verbose=False, random_state=0)
    print('%s: DeepScence %.1f min'%(tag,(time.time()-t)/60), flush=True)
    pd.DataFrame({'cell_id':ids_kept,'deepscence_score':np.round(res.obs['ds'].to_numpy(),5)}
                 ).to_csv(PROC+'deepscence_%s.csv'%tag, index=False)
    print('wrote deepscence_%s.csv'%tag, flush=True)
