#!/usr/bin/env python3
"""D-A (v2): marker-evidence cell type assignment.

FIX over v1: v1 z-scored each cell type's AGGREGATE score across clusters, which normalises away
magnitude -- a 2-gene marker set that is barely expressed anywhere still produced z=+1.3 in some
cluster and won the argmax. That put 28% of a liver into 'Erythroid cells' and 6% into 'Mast'.
v2 instead:
  (1) z-scores EACH MARKER GENE across clusters, then averages over the set, so every gene
      contributes equally and a set cannot win on one noisy gene;
  (2) GATES on evidence: the marker set must actually be detected in that cluster
      (mean detection rate >= MIN_DET) and reach MIN_Z;
  (3) DROPS marker sets with < MIN_MARKERS genes on panel as untrustworthy, rather than
      letting them compete. Types dropped for this reason are reported, not silently omitted.
"""
import sys, json, numpy as np, pandas as pd, scanpy as sc, warnings
sys.path.insert(0,'/workspace/code'); warnings.filterwarnings('ignore')
from markers_mouse_liver import MARKERS
import sasp_io as _io
MIN_MARKERS, MIN_DET, MIN_Z, MIN_MARGIN = 4, 0.05, 0.50, 0.20
INTERIM='/workspace/data/interim/'; PROC='/workspace/data/processed/'

EXTRA = {   # authors' vocabulary; zonated endothelial subtypes need their own markers
 'Central venous LSECs':'Wnt2 Wnt9b Rspo3 Thbd Kit Fabp4'.split(),
 'Portal endothelial cells':'Vwf Gja5 Jam2 Ntn4 Efnb2 Sox17 Dll4 Epas1'.split(),
}
M = {**MARKERS, **EXTRA}

for tag in sys.argv[1:]:
    print('\n'+'='*95); print('SAMPLE',tag); print('='*95)
    B = _io.add_leiden(_io.load(tag), tag)
    cl = B.obs['leiden'].astype(str).values
    clusters = sorted(set(cl), key=lambda s:int(s))
    avail, dropped = {}, {}
    for ct,gs in M.items():
        on=[g for g in gs if g in B.var_names]
        (avail if len(on)>=MIN_MARKERS else dropped)[ct]=on
    if dropped:
        print('DROPPED (fewer than %d markers on panel -- NOT assignable, reported as absent):'%MIN_MARKERS)
        for k,v in dropped.items(): print('   %-24s %d on-panel: %s'%(k,len(v),','.join(v)))
    allg = sorted({g for v in avail.values() for g in v})
    X = B[:, allg].X
    X = np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
    dfm = pd.DataFrame(X, columns=allg); dfm['cl']=cl
    cmean = dfm.groupby('cl', observed=True).mean().loc[clusters]
    cdet  = pd.DataFrame((X>0).astype(np.float32), columns=allg).assign(cl=cl)\
              .groupby('cl', observed=True).mean().loc[clusters]
    Zg = (cmean - cmean.mean())/(cmean.std()+1e-9)     # z PER GENE across clusters
    S   = pd.DataFrame({ct: Zg[v].mean(axis=1)  for ct,v in avail.items()})
    DET = pd.DataFrame({ct: cdet[v].mean(axis=1) for ct,v in avail.items()})
    Sg = S.where(DET>=MIN_DET, other=-np.inf)
    top=Sg.idxmax(axis=1); tv=Sg.max(axis=1)
    second=Sg.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
    margin=tv-second
    lab=pd.Series(np.where((tv>=MIN_Z)&(margin>=MIN_MARGIN), top, 'Unknown'), index=S.index)
    conf=pd.Series(np.round(np.clip(margin/(np.abs(tv)+1e-9),0,1),3), index=S.index)
    conf[lab=='Unknown']=0.0
    n=pd.Series(cl).value_counts()
    print('\n%-5s %-8s %-28s %7s %7s %6s  %s'%('cl','n','cell_type','z','margin','det','runner-up'))
    for c in clusters:
        print('%-5s %-8d %-28s %+7.2f %7.2f %6.2f  %s'%(c,n[c],lab[c],tv[c] if np.isfinite(tv[c]) else -9,
              margin[c] if np.isfinite(margin[c]) else 0, DET.loc[c,top[c]] if top[c] in DET else 0,
              Sg.loc[c].nlargest(2).index[1]))
    B.obs['cell_type']=lab.reindex(cl).values
    B.obs['cell_type_confidence']=conf.reindex(cl).values
    tot=B.obs.cell_type.value_counts()
    print('\n--- composition ---')
    for k,v in tot.items(): print('   %-28s %7d  %5.2f%%'%(k,v,100*v/len(B)))
    out=B.obs[['leiden','cell_type','cell_type_confidence']].copy()
    out.index.name='cell_id'; out.reset_index().to_csv(PROC+'celltypes_%s.csv'%tag,index=False)
    S.round(3).to_csv(PROC+'cluster_celltype_zscores_%s.csv'%tag)
    json.dump({'assignable':avail,'dropped_thin_marker_sets':dropped},
              open(PROC+'markers_used_%s.json'%tag,'w'), indent=1)
    print('wrote celltypes_%s.csv'%tag)
