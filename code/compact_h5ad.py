#!/usr/bin/env python3
"""Shrink the interim h5ad: drop the scaled DENSE X (4.8 GB/sample) and .raw duplication.
Keep sparse log-normalised X, sparse counts layer, PCA, and obs. ~11 GB -> <1.5 GB."""
import sys, scanpy as sc, anndata as ad, numpy as np, os
from scipy.sparse import csr_matrix
for tag in sys.argv[1:]:
    p='/workspace/data/interim/%s.h5ad'%tag
    print(tag,'reading %.1f GB ...'%(os.path.getsize(p)/1e9), flush=True)
    A=sc.read_h5ad(p)
    R=A.raw.to_adata()                              # log1p CP-median normalised, sparse
    N=ad.AnnData(csr_matrix(R.X), obs=A.obs.copy(), var=R.var.copy())
    N.layers['counts']=csr_matrix(A.layers['counts'])
    N.obsm['X_pca']=A.obsm['X_pca']
    if 'neighbors' in A.uns: N.uns['neighbors']=A.uns['neighbors']
    for k in ['connectivities','distances']:
        if k in A.obsp: N.obsp[k]=A.obsp[k]
    tmp=p+'.tmp'; N.write(tmp); os.replace(tmp,p)
    print('  -> %.2f GB'%(os.path.getsize(p)/1e9), flush=True)
