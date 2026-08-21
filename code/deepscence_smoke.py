import sys, time, numpy as np, pandas as pd, scanpy as sc, anndata as ad, h5py, csv
from scipy.sparse import csc_matrix, csr_matrix
sys.path.insert(0,'/workspace/code/_shims')
import DeepScence.api as api

orth=dict((r['mouse_symbol'],r['human_symbol']) for r in csv.DictReader(open('/workspace/genesets/mouse_human_orthologs_MGI.csv')))
RAW='/workspace/data/raw/7250_liver_sham_Male_26-U1/'
f=h5py.File(RAW+'cell_feature_matrix.h5','r')
ft=np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
ids=np.array([x.decode() for x in f['matrix/features/id'][:]]); nm=np.array([x.decode() for x in f['matrix/features/name'][:]])
bc=np.array([x.decode() for x in f['matrix/barcodes'][:]])
keep=(ft=='Gene Expression')&np.char.startswith(ids,'ENSMUSG')
M=csc_matrix((f['matrix/data'][:].astype(np.float32),f['matrix/indices'][:].astype(np.int32),
              f['matrix/indptr'][:].astype(np.int64)),shape=tuple(f['matrix/shape'][:]))
X=csr_matrix(M[np.where(keep)[0],:].T); genes=nm[keep]
A=ad.AnnData(X, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=genes))
A=A[np.asarray(A.X.sum(1)).ravel()>=20].copy()
hs=np.array([orth.get(g,'') for g in A.var_names])
n_map=(hs!='').sum(); print('mouse genes mapped to 1:1 human ortholog: %d / %d'%(n_map,A.n_vars))
A=A[:, hs!=''].copy(); A.var_names=hs[hs!='']; A.var_names_make_unique()
core=pd.read_csv('/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv')
c5=core[core.occurrence>=5]
print('CoreScence(occurrence>=5): %d genes; on our ortholog-mapped panel: %d'
      %(len(c5), c5.gene_symbol.isin(A.var_names).sum()))
print('  up  :', c5[(c5.direction=="up")&(c5.gene_symbol.isin(A.var_names))].gene_symbol.tolist())
print('  down:', c5[(c5.direction=="down")&(c5.gene_symbol.isin(A.var_names))].gene_symbol.tolist())
sub=A[np.random.default_rng(0).choice(A.n_obs, 20000, replace=False)].copy()
t=time.time(); print('\nrunning DeepScence on 20,000 cells, denoise=False ...', flush=True)
res=api.DeepScence(sub, denoise=False, verbose=True, random_state=0)
print('elapsed %.1f s'%(time.time()-t))
print('obs cols:', [c for c in res.obs.columns])
print(res.obs.filter(like='ds').describe().to_string() if len(res.obs.filter(like='ds').columns) else res.obs.head())
