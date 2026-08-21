#!/usr/bin/env python3
"""D-A: Leiden clustering + marker-based cell type annotation for GSE310392 mouse Xenium."""
import sys, os, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix
sys.path.insert(0,'/workspace/code')
from markers_mouse_liver import MARKERS
warnings.filterwarnings('ignore')
sc.settings.n_jobs = 48
sc.settings.verbosity = 2

SAMPLES = {'7250_liver_sham_Male_26-U1':'sham', '7259_liver_sbr_Male_26-U1':'sbr'}
RAW='/workspace/data/raw/'; INTERIM='/workspace/data/interim/'; PROC='/workspace/data/processed/'
MIN_COUNTS, MIN_GENES, RES = 20, 5, 1.0

def load(sample):
    p = RAW+sample+'/cell_feature_matrix.h5'
    f = h5py.File(p,'r')
    ft = np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids= np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm = np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc = np.array([x.decode() for x in f['matrix/barcodes'][:]])
    ge = (ft=='Gene Expression')
    # DROP the 9 genotyping probes (non-ENSMUSG ids among Gene Expression features)
    geno = ge & ~np.char.startswith(ids,'ENSMUSG')
    print('  dropping %d genotyping probes: %s'%(geno.sum(), ', '.join(nm[geno])))
    keep = ge & ~geno
    M = csc_matrix((f['matrix/data'][:].astype(np.float32),
                    f['matrix/indices'][:].astype(np.int32),
                    f['matrix/indptr'][:].astype(np.int64)), shape=tuple(f['matrix/shape'][:]))
    X = csr_matrix(M[np.where(keep)[0],:].T)     # cells x genes
    A = ad.AnnData(X, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=nm[keep]))
    A.var['ensembl']=ids[keep]
    cells = pd.read_parquet(RAW+sample+'/cells.parquet').set_index('cell_id')
    A.obs = A.obs.join(cells)
    A.obs.rename(columns={'x_centroid':'x_um','y_centroid':'y_um'}, inplace=True)
    return A

for sample, tag in SAMPLES.items():
    print('\n'+'='*90); print('SAMPLE', sample); print('='*90, flush=True)
    A = load(sample)
    print('  loaded', A.shape, flush=True)
    sc.pp.calculate_qc_metrics(A, inplace=True, percent_top=None, log1p=False)
    n0=A.n_obs
    A = A[(A.obs.total_counts>=MIN_COUNTS)&(A.obs.n_genes_by_counts>=MIN_GENES)].copy()
    print('  QC: %d -> %d cells (%.2f%% kept; counts>=%d & genes>=%d)'%(n0,A.n_obs,100*A.n_obs/n0,MIN_COUNTS,MIN_GENES), flush=True)
    A.layers['counts']=A.X.copy()
    sc.pp.normalize_total(A)          # target = median counts
    sc.pp.log1p(A)
    A.raw = A
    sc.pp.scale(A, max_value=10)
    sc.tl.pca(A, n_comps=50, svd_solver='arpack')
    print('  PCA done', flush=True)
    sc.pp.neighbors(A, n_neighbors=15, n_pcs=50)
    print('  neighbors done', flush=True)
    sc.tl.leiden(A, resolution=RES, key_added='leiden', flavor='igraph', n_iterations=2, directed=False)
    print('  leiden done: %d clusters'%A.obs.leiden.nunique(), flush=True)

    # ---- marker scoring on log-normalised values (A.raw), with matched control sets
    B = A.raw.to_adata(); B.obs['leiden']=A.obs['leiden'].values
    avail={}
    for ct,gs in MARKERS.items():
        on=[g for g in gs if g in B.var_names]
        if len(on)<2: print('  !! %s has <2 markers on panel, skipped'%ct); continue
        avail[ct]=on
        sc.tl.score_genes(B, on, score_name='sc_'+ct, ctrl_size=max(50,len(on)*10))
    S = B.obs[['sc_'+c for c in avail]].copy(); S.columns=list(avail)
    prof = S.groupby(A.obs.leiden.values, observed=True).mean()
    Z = (prof - prof.mean())/ (prof.std()+1e-9)          # z across clusters, per cell type
    top  = Z.idxmax(axis=1); tv = Z.max(axis=1)
    second = Z.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
    margin = tv - second
    lab = pd.Series(np.where((tv>=1.0)&(margin>=0.25), top, 'Unknown'), index=Z.index)
    conf = pd.Series(np.clip(margin/ (tv.abs()+1e-9),0,1), index=Z.index).round(3)
    conf[lab=='Unknown']=0.0
    A.obs['cell_type']=lab.reindex(A.obs.leiden.values).values
    A.obs['cell_type_confidence']=conf.reindex(A.obs.leiden.values).values

    print('\n  --- cluster -> cell type (top marker z, margin over 2nd) ---')
    cnt=A.obs.leiden.value_counts()
    for cl in Z.index:
        print('   cl %-3s n=%-7d %-28s z=%+.2f margin=%.2f  2nd=%s'%(
            cl,cnt.get(cl,0),lab[cl],tv[cl],margin[cl],Z.loc[cl].nlargest(2).index[1]))
    print('\n  --- cell type totals ---'); print(A.obs.cell_type.value_counts().to_string())

    sc.tl.rank_genes_groups(A, 'leiden', method='wilcoxon', n_genes=12)
    top_de = pd.DataFrame(A.uns['rank_genes_groups']['names'])
    top_de.to_csv(PROC+'cluster_top_markers_%s.csv'%tag, index=False)
    Z.round(3).to_csv(PROC+'cluster_celltype_zscores_%s.csv'%tag)

    out = A.obs[['leiden','cell_type','cell_type_confidence']].copy()
    out.index.name='cell_id'; out.reset_index().to_csv(PROC+'celltypes_%s.csv'%tag, index=False)
    print('  wrote celltypes_%s.csv'%tag, flush=True)
    del B, S
    A.write(INTERIM+'%s.h5ad'%tag)
    print('  wrote interim h5ad', flush=True)
    json.dump({k:v for k,v in avail.items()}, open(PROC+'markers_used_%s.json'%tag,'w'), indent=1)
print('\nALL DONE')
