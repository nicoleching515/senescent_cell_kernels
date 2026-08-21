"""Shared loader: read counts straight from the Xenium h5 (fast, ~100 MB) instead of a
multi-GB h5ad. /workspace is a NETWORK volume, so large intermediates are I/O-bound."""
import numpy as np, pandas as pd, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix
RAW='/workspace/data/raw/'
SAMP={'sham':'7250_liver_sham_Male_26-U1','sbr':'7259_liver_sbr_Male_26-U1'}
# Phase 3: allow any full sample directory name as its own tag, so the same
# loaders work for the sections whose annotation CSVs are keyed on the full name.
import os as _os
for _d in sorted(_os.listdir(RAW)):
    if _os.path.exists(_os.path.join(RAW,_d,'cells.parquet')):
        SAMP.setdefault(_d,_d)

MIN_COUNTS, MIN_GENES = 20, 5
def resolve(tag):
    """tag may be a short alias ('sham','sbr') or a full sample directory name."""
    import os
    if tag in SAMP: return SAMP[tag]
    if os.path.isdir(RAW+tag): return tag
    raise KeyError('unknown sample tag %r'%tag)
def load(tag, lognorm=True, with_cells=True):
    d=RAW+resolve(tag)+'/'
    f=h5py.File(d+'cell_feature_matrix.h5','r')
    ft=np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids=np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm=np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc=np.array([x.decode() for x in f['matrix/barcodes'][:]])
    keep=(ft=='Gene Expression')&np.char.startswith(ids,'ENSMUSG')  # drops 9 genotyping probes
    M=csc_matrix((f['matrix/data'][:].astype(np.float32),f['matrix/indices'][:].astype(np.int32),
                  f['matrix/indptr'][:].astype(np.int64)),shape=tuple(f['matrix/shape'][:]))
    A=ad.AnnData(csr_matrix(M[np.where(keep)[0],:].T), obs=pd.DataFrame(index=bc),
                 var=pd.DataFrame(index=nm[keep]))
    tot=np.asarray(A.X.sum(1)).ravel(); ng=np.asarray((A.X>0).sum(1)).ravel()
    A=A[(tot>=MIN_COUNTS)&(ng>=MIN_GENES)].copy()
    if with_cells:
        c=pd.read_parquet(d+'cells.parquet').set_index('cell_id')
        c=c.rename(columns={'x_centroid':'x_um','y_centroid':'y_um'})
        A.obs=A.obs.join(c)
    A.layers['counts']=A.X.copy()
    if lognorm:
        import scanpy as sc
        sc.pp.normalize_total(A); sc.pp.log1p(A)
    return A
def add_leiden(A, tag):
    ct=pd.read_csv('/workspace/data/processed/celltypes_%s.csv'%tag).set_index('cell_id')
    A.obs['leiden']=ct['leiden'].reindex(A.obs_names).astype(str).values
    return A
