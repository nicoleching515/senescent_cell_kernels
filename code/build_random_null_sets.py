#!/usr/bin/env python3
"""Tier E, third component (Sec 9): 500 random gene sets per response module, matched on
SIZE and MEAN EXPRESSION. Gives an empirical null for effect size that is independent of the
housekeeping-gene problem (Phase 1 sec 3.3).

Matching: genes binned into 20 quantile bins of mean log-normalised expression computed on
GSM9295284 (sham, all QC-passed cells). For each module, each random set reproduces the module's
per-bin membership counts exactly, drawing from the same bin and excluding the module's own genes.
"""
import numpy as np, h5py, os, glob, csv
from scipy.sparse import csc_matrix

RAW='/workspace/data/raw/7250_liver_sham_Male_26-U1/'
OUT='/workspace/genesets/E3_random_matched/'
os.makedirs(OUT, exist_ok=True)
N_SETS, N_BINS, SEED = 500, 20, 20260820

f=h5py.File(RAW+'cell_feature_matrix.h5','r')
ft=np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
ids=np.array([x.decode() for x in f['matrix/features/id'][:]])
nm=np.array([x.decode() for x in f['matrix/features/name'][:]])
keep=(ft=='Gene Expression') & np.char.startswith(ids,'ENSMUSG')   # drop 9 genotyping probes
M=csc_matrix((f['matrix/data'][:].astype(np.float32),f['matrix/indices'][:].astype(np.int32),
              f['matrix/indptr'][:].astype(np.int64)),shape=tuple(f['matrix/shape'][:]))
X=M[np.where(keep)[0],:].tocsr(); genes=nm[keep]
tot=np.asarray(X.sum(axis=0)).ravel(); ok=tot>=20
X=X[:,ok]; tot=tot[ok]
Xn=X.multiply((np.median(tot)/tot)[None,:].astype(np.float32)).tocsr(); Xn.data=np.log1p(Xn.data)
mean_expr=np.asarray(Xn.mean(axis=1)).ravel()
print('genes %d, cells %d'%(len(genes), Xn.shape[1]))

order=np.argsort(mean_expr)
binof=np.empty(len(genes),int)
for b,chunk in enumerate(np.array_split(order,N_BINS)): binof[chunk]=b
gidx={g:i for i,g in enumerate(genes)}
rng=np.random.default_rng(SEED)

with open('/workspace/genesets/E3_random_matched/_expression_bins.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['gene','mean_lognorm_expr','bin'])
    for i,g in enumerate(genes): w.writerow([g, round(float(mean_expr[i]),6), int(binof[i])])

for path in sorted(glob.glob('/workspace/genesets/B_*.txt')):
    mod=os.path.basename(path)[2:-4]
    G=[l.strip() for l in open(path) if l.strip() and l.strip() in gidx]
    gi=np.array([gidx[g] for g in G])
    counts=np.bincount(binof[gi], minlength=N_BINS)
    pool={b:np.setdiff1d(np.where(binof==b)[0], gi) for b in range(N_BINS)}
    short=[b for b in range(N_BINS) if counts[b]>len(pool[b])]
    if short: print('  !! %s: bins %s have fewer candidates than needed'%(mod,short))
    out=[]
    for s in range(N_SETS):
        pick=[]
        for b in range(N_BINS):
            if counts[b]==0: continue
            pick += list(rng.choice(pool[b], counts[b], replace=False))
        out.append([genes[i] for i in pick])
    with open(OUT+'%s.tsv'%mod,'w') as fh:
        for s,gs in enumerate(out): fh.write('%s\t%d\t%s\n'%(mod,s,','.join(gs)))
    obs=mean_expr[gi].mean()
    nullm=np.array([mean_expr[[gidx[g] for g in gs]].mean() for gs in out])
    print('%-24s n=%3d  module mean_expr=%.4f  null mean=%.4f (sd %.4f)  -> matched'
          %(mod,len(G),obs,nullm.mean(),nullm.std()))
print('\nwrote', OUT)
