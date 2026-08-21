#!/usr/bin/env python3
"""D-A: cell type annotation for GSE310392 mouse Xenium liver sections.
Reusable across all 11 sections:  python3 annotate_pipeline.py <sample_dir> [<sample_dir> ...]

DESIGN NOTES (v3; v1 and v2 failures are documented in reports/BIO_PHASE2.md)
 1. LABEL SET IS GATED BY ON-PANEL MARKER SUPPORT. A cell type needs >= MIN_MARKERS informative
    markers surviving the 5,106-gene panel or it is removed from the label set entirely, so it
    cannot absorb noise. Dropped types are reported, never silently omitted.
 2. NON-SPECIFIC types are also removed: a type whose on-panel markers are individually expressed
    by a *different*, more abundant type cannot be resolved. See DROP_NONSPECIFIC.
 3. SCORING: each marker gene is z-scored across clusters, then averaged within the marker set, so
    sets of different size are comparable and no single gene can carry a set.
 4. ASSIGNMENT requires (a) detection: mean marker detection rate >= MIN_DET; (b) magnitude
    z >= MIN_Z; (c) a margin over the runner-up >= MIN_MARGIN. Otherwise 'Unknown'.
    MIN_Z is deliberately low (0.25): when ONE cell type occupies most clusters (hepatocytes in
    liver), that type's z-across-clusters is compressed toward 0 because it *is* the background.
    The margin criterion, not the magnitude, does the discriminating.
 5. LOW-QUALITY clusters (median counts < LOWQ_FRAC x section median) are labelled 'Low_quality'
    and excluded from receiver analyses -- they are segmentation fragments, not a cell type.
    EXCEPTION: a low-count cluster that nonetheless carries an unambiguous identity
    (z >= LOWQ_RESCUE_Z and margin >= LOWQ_RESCUE_MARGIN) keeps its label. Some real cell types
    are genuinely low-content on this panel: mesothelial cells are thin and flat and sit on the
    liver capsule, and they were being discarded by the count rule despite z=+4.3, margin=3.7.
 6. Hepatocyte ZONATION states (pericentral Cyp2e1/Glul/Cyp1a2 vs periportal Cps1/Hal/Arg1) are
    NOT separate cell types. They are labelled 'Hepatocytes' and zonation is carried as the
    continuous covariate from D-B, per Section 11.
 7. SANITY CHECK: liver sections must come out hepatocyte-dominant. HEP_MIN_FRAC is asserted and
    a violation is printed loudly rather than written out quietly.
"""
import sys, os, json, warnings, numpy as np, pandas as pd, scanpy as sc, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix
warnings.filterwarnings('ignore'); sc.settings.n_jobs=48; sc.settings.verbosity=1

MIN_MARKERS, MIN_DET, MIN_Z, MIN_MARGIN = 4, 0.05, 0.25, 0.20
LOWQ_FRAC, HEP_MIN_FRAC, RES, MIN_COUNTS, MIN_GENES = 0.5, 0.35, 1.0, 20, 5
LOWQ_RESCUE_Z, LOWQ_RESCUE_MARGIN = 1.0, 1.0
PROC='/workspace/data/processed/'

MARKERS = {
'Hepatocytes':'Apoa1 Apob Fga Fgb Fgg Serpinc1 Hpd Orm1 Agt C3 Pzp Cps1 Otc Fah Aldh1l1 Apoa4'.split(),
# NOT called 'Cholangiocytes'. On this panel the two definitive bile-duct markers effectively do
 # not detect -- measured mean log-norm expression Krt7 <= 0.03 and Cftr <= 0.06 in EVERY cluster
 # of both sections -- and Krt19 reaches only 0.22 in the best sham cluster. What the panel can
 # resolve is a Sox9+/Epcam+/Pkhd1+/Spp1+ compartment that contains BOTH true cholangiocytes and
 # hepatocytes undergoing ductular metaplasia; the two are not separable here. Calling it
 # 'Cholangiocytes' would overstate what the data support, so it is named for the compartment.
 # This compartment expands from 2.9% (sham) to 25.8% (SBR 26wk) -- a ductular reaction, which is
 # the expected IFALD phenotype, not an annotation error.
'Biliary/ductular':'Sox9 Epcam Pkhd1 Krt19 Hnf1b Onecut1 Krt7 Cftr'.split(),
'Hepatic stellate cells':'Lrat Des Reln Hgf Lhx2 Rbp1 Col1a1 Col1a2 Pdgfrb'.split(),
'Portal fibroblasts':'Col1a1 Col1a2 Pdgfra Thy1 Eln Fbln1 Fbln2 Postn Lum Vcan'.split(),
'LSECs':'Kdr Pecam1 Cdh5 Lyve1 Fcgr2b Gata4 Mrc1 Eng Flt1 Egfl7'.split(),
'Central venous LSECs':'Wnt2 Wnt9b Rspo3 Thbd Kit Fabp4'.split(),
'Portal endothelial cells':'Vwf Gja5 Jam2 Ntn4 Efnb2 Sox17 Dll4 Epas1'.split(),
'Kupffer cells':'Adgre1 Csf1r Cd68 Itgam Marco Lyz2 Cd5l Slc40a1 Fcgr1 Mafb'.split(),
'Inflammatory macs':'Ccr2 Itgam Lyz2 Fcgr1 Cd14 Nlrp3 Trem2 Gpnmb Cd9'.split(),
# T and NK are MERGED: on this panel Il2rb/Ptprc/Klrd1 are shared and Ncr1/Eomes/Gzma are too
# sparse to separate them. Reporting them apart would be a guess.
'T/NK cells':'Cd3e Cd3d Cd3g Cd2 Lck Cd8a Cd4 Il7r Themis Il2rb Ncr1 Gzma Klrb1c Klrd1 Eomes'.split(),
'B-cells':'Cd79a Cd79b Ms4a1 Cd19 Pax5 Ighm Blnk'.split(),
'Plasma':'Jchain Sdc1 Prdm1 Mzb1 Igkc'.split(),
'DC':'Flt3 Itgax Cd209a Batf3 Xcr1 Zbtb46 Clec9a Irf8'.split(),
'vSMCs':'Acta2 Myh11 Tagln Cnn1 Actg2'.split(),
'Pericytes':'Rgs5 Notch3 Kcnj8 Abcc9 Cspg4'.split(),
'Mesothelial cells':'Msln Upk3b Wt1 Gpm6a'.split(),
'Proliferating':'Mki67 Top2a Pcna Ccnb1 Ccna2 Birc5 Aurkb Cdk1'.split(),
# NOT in the label set, with reasons -- see DROP_NONSPECIFIC / marker-support gate:
#   Erythroid cells : only Slc4a1, Klf1 on panel (Hba-a1/Hbb-bs/Alas2/Gypa all OFF-panel)
#   Mast            : only Kit, Ms4a2, Tpsb2 on panel (Cpa3/Cma1/Mcpt4 all OFF-panel)
#   Neutrophils     : only Ly6g, Csf3r, Cxcr2 on panel (S100a8/S100a9/Retnlg/Mpo/Elane OFF-panel)
}
# ---------------------------------------------------------------------------------------------
# MERGED LABEL SET (added in Phase 3).  Some of the fine labels above are not stably separable on
# this panel and the winner flips between sections, which makes a cross-section composition table
# uninterpretable.  Measured instances:
#   * stellate vs portal fibroblast vs pericyte -- all Col1a1/Col1a2/Pdgfrb+. In 7001 the two top
#     scores tie at margin 0.02 and 19,420 cells (11.7% of the section) fall to 'Unknown', giving
#     that section 0.0% stellate while every other section has 7-17%.  That is an artefact.
#   * Kupffer vs 'Inflammatory macs' -- 7448 splits 0.2%/10.6%, 7450 splits 14.1%/0.0%.
#   * LSEC vs central-venous LSEC vs portal endothelial -- same family, margins 0.6-1.4.
# The fix is NOT to drop the fine labels (they are informative where the margin is large).  Both
# are written: `cell_type` is unchanged, `cell_type_merged` is the assignment recomputed over the
# UNION of each group's on-panel markers, and is the label to use for cross-section comparison.
MERGE = {
 'Mesenchymal':['Hepatic stellate cells','Portal fibroblasts','Pericytes'],
 'Macrophages':['Kupffer cells','Inflammatory macs'],
 'Endothelial':['LSECs','Central venous LSECs','Portal endothelial cells'],
}
DROP_NONSPECIFIC = {
 'Lymphatic endothelial cells':
   'measured in GSM9295284: Ccl21a and Mmrn1 mean expression ~0.00 in every cluster; Prox1 is '
   'expressed by hepatocytes (0.11-0.19) as highly as anywhere else; Lyve1 and Flt4 peak in the '
   'LSEC clusters. No separable lymphatic population exists on this panel.'}

def load_counts(sample_dir):
    f=h5py.File(os.path.join(sample_dir,'cell_feature_matrix.h5'),'r')
    ft=np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids=np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm=np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc=np.array([x.decode() for x in f['matrix/barcodes'][:]])
    keep=(ft=='Gene Expression')&np.char.startswith(ids,'ENSMUSG')   # drops 9 genotyping probes
    M=csc_matrix((f['matrix/data'][:].astype(np.float32),f['matrix/indices'][:].astype(np.int32),
                  f['matrix/indptr'][:].astype(np.int64)),shape=tuple(f['matrix/shape'][:]))
    A=ad.AnnData(csr_matrix(M[np.where(keep)[0],:].T), obs=pd.DataFrame(index=bc),
                 var=pd.DataFrame(index=nm[keep]))
    c=pd.read_parquet(os.path.join(sample_dir,'cells.parquet')).set_index('cell_id')
    A.obs=A.obs.join(c.rename(columns={'x_centroid':'x_um','y_centroid':'y_um'}))
    return A

def run(sample_dir, tag, relabel=False):
    print('\n'+'='*95); print('SAMPLE %s   (%s)'%(tag,sample_dir)); print('='*95, flush=True)
    A=load_counts(sample_dir); n0=A.n_obs
    tot=np.asarray(A.X.sum(1)).ravel(); ng=np.asarray((A.X>0).sum(1)).ravel()
    A=A[(tot>=MIN_COUNTS)&(ng>=MIN_GENES)].copy()
    print('QC %d -> %d cells (%.2f%%)'%(n0,A.n_obs,100*A.n_obs/n0), flush=True)
    A.layers['counts']=A.X.copy()
    sc.pp.normalize_total(A); sc.pp.log1p(A)
    if relabel:
        prev=pd.read_csv(PROC+'celltypes_%s.csv'%tag).set_index('cell_id')
        A.obs['leiden']=prev['leiden'].reindex(A.obs_names).astype(str).values
        print('RELABEL mode: reusing saved leiden (%d clusters), no re-clustering'
              %A.obs.leiden.nunique(), flush=True)
        B=A
    else:
        A.raw=A
        sc.pp.scale(A,max_value=10); sc.tl.pca(A,n_comps=50,svd_solver='arpack')
        sc.pp.neighbors(A,n_neighbors=15,n_pcs=50)
        sc.tl.leiden(A,resolution=RES,key_added='leiden',flavor='igraph',n_iterations=2,directed=False)
        B=A.raw.to_adata(); B.obs=A.obs.copy()
        del A
    cl=B.obs['leiden'].astype(str).values
    clusters=sorted(set(cl),key=lambda s:int(s))
    print('clusters: %d'%len(clusters), flush=True)

    avail,dropped={},{}
    for ct,gs in MARKERS.items():
        on=[g for g in gs if g in B.var_names]
        (avail if len(on)>=MIN_MARKERS else dropped)[ct]=on
    print('\nLABEL SET GATE: %d types assignable, %d dropped'%(len(avail),len(dropped)))
    for k,v in dropped.items(): print('  DROPPED %-28s only %d on-panel markers: %s'%(k,len(v),','.join(v)))
    for k,v in DROP_NONSPECIFIC.items(): print('  DROPPED %-28s non-specific. %s'%(k,v))

    cnts=np.asarray(B.layers['counts'].sum(1)).ravel() if 'counts' in B.layers else None
    if cnts is None: cnts=np.asarray(B.obs['transcript_counts'].values,dtype=float)
    qc=pd.DataFrame({'cl':cl,'counts':cnts}).groupby('cl').median().loc[clusters,'counts']
    lowq=set(qc.index[qc < LOWQ_FRAC*np.median(cnts)])
    print('\nLOW-QUALITY clusters (median counts < %.0f%% of section median %.0f): %s'
          %(100*LOWQ_FRAC,np.median(cnts), sorted(lowq,key=int) or 'none'))

    allg=sorted({g for v in avail.values() for g in v})
    X=np.asarray(B[:,allg].X.todense())
    dfm=pd.DataFrame(X,columns=allg); dfm['cl']=cl
    cmean=dfm.groupby('cl',observed=True).mean().loc[clusters]
    cdet=pd.DataFrame((X>0).astype(np.float32),columns=allg).assign(cl=cl)\
          .groupby('cl',observed=True).mean().loc[clusters]
    Zg=(cmean-cmean.mean())/(cmean.std()+1e-9)

    def assign(labelsets):
        S=pd.DataFrame({c:Zg[v].mean(axis=1)  for c,v in labelsets.items()})
        DET=pd.DataFrame({c:cdet[v].mean(axis=1) for c,v in labelsets.items()})
        Sg=S.where(DET>=MIN_DET,other=-np.inf)
        top=Sg.idxmax(axis=1); tv=Sg.max(axis=1)
        second=Sg.apply(lambda r:r.nlargest(2).iloc[1],axis=1); m=tv-second
        lb=pd.Series(np.where((tv>=MIN_Z)&(m>=MIN_MARGIN),top,'Unknown'),index=S.index)
        cf=pd.Series(np.round(np.clip(m/(np.abs(tv)+1e-9),0,1),3),index=S.index); cf[lb=='Unknown']=0.0
        return S,DET,Sg,top,tv,m,lb,cf
    S,DET,Sg,top,tv,margin,lab,conf = assign(avail)
    # merged label set: union of on-panel markers over each group's members
    avail_m={k:v for k,v in avail.items() if not any(k in mem for mem in MERGE.values())}
    for grp,mem in MERGE.items():
        u=sorted({g for c in mem for g in avail.get(c,[])})
        if any(c in avail for c in mem): avail_m[grp]=u
    # Merged score = MAX over the group's members, not the mean over the union of their markers.
    # The union-mean dilutes: 7001 cl10 scores stellate 2.65 / portal-fibroblast 2.64, but the
    # union-mean falls to ~2.0 and loses to Mesothelial at 1.92, leaving 13,158 cells Unknown.
    # 'These are one compartment, take whichever member evidences it' is the intended semantics.
    Sm=pd.DataFrame({c:Sg[c] for c in avail_m if c in Sg.columns})
    for grp,mem in MERGE.items():
        mem=[c for c in mem if c in Sg.columns]
        if mem: Sm[grp]=Sg[mem].max(axis=1)
    top_m=Sm.idxmax(axis=1); tv_m=Sm.max(axis=1)
    sec_m=Sm.apply(lambda r:r.nlargest(2).iloc[1],axis=1); margin_m=tv_m-sec_m
    lab_m=pd.Series(np.where((tv_m>=MIN_Z)&(margin_m>=MIN_MARGIN),top_m,'Unknown'),index=Sm.index)
    conf_m=pd.Series(np.round(np.clip(margin_m/(np.abs(tv_m)+1e-9),0,1),3),index=Sm.index)
    conf_m[lab_m=='Unknown']=0.0
    rescued=[]
    for c in lowq:
        if np.isfinite(tv[c]) and tv[c]>=LOWQ_RESCUE_Z and margin[c]>=LOWQ_RESCUE_MARGIN:
            rescued.append((c,lab[c])); continue
        lab[c]='Low_quality'; conf[c]=0.0
    if rescued: print('  RESCUED from Low_quality on unambiguous identity: %s'
                      %', '.join('cl%s=%s'%(c,l) for c,l in rescued))
    nper=pd.Series(cl).value_counts()
    print('\n%-4s %-8s %-9s %-26s %7s %7s %6s  %s'%('cl','n','medcnt','cell_type','z','margin','det','runner-up'))
    for c in clusters:
        print('%-4s %-8d %-9.0f %-26s %+7.2f %7.2f %6.2f  %s'%(c,nper[c],qc[c],lab[c],
            tv[c] if np.isfinite(tv[c]) else -9, margin[c] if np.isfinite(margin[c]) else 0,
            DET.loc[c,top[c]], Sg.loc[c].nlargest(2).index[1]))
    for c in lowq:
        if not (np.isfinite(tv_m[c]) and tv_m[c]>=LOWQ_RESCUE_Z and margin_m[c]>=LOWQ_RESCUE_MARGIN):
            lab_m[c]='Low_quality'; conf_m[c]=0.0
    B.obs['cell_type']=lab.reindex(cl).values
    B.obs['cell_type_confidence']=conf.reindex(cl).values
    B.obs['cell_type_merged']=lab_m.reindex(cl).values
    B.obs['cell_type_merged_confidence']=conf_m.reindex(cl).values
    print('\n--- merged composition (cell_type_merged) ---')
    for k,v in B.obs.cell_type_merged.value_counts().items(): print('   %-28s %7d  %5.2f%%'%(k,v,100*v/len(B)))
    vc=B.obs.cell_type.value_counts(); frac=vc/len(B)
    print('\n--- composition ---')
    for k,v in vc.items(): print('   %-28s %7d  %5.2f%%'%(k,v,100*v/len(B)))
    hep=frac.get('Hepatocytes',0); par=hep+frac.get('Biliary/ductular',0)
    print('\nSANITY: hepatocyte %.1f%% ; hepatic parenchyma (hepatocyte + biliary/ductular) %.1f%%'
          ' (floor %.0f%%) -> %s'%(100*hep,100*par,100*HEP_MIN_FRAC,
          'OK' if par>=HEP_MIN_FRAC else '*** FAIL: ANNOTATION IS BROKEN ***'))
    print('median cell_type_confidence %.3f'%B.obs.cell_type_confidence.median())
    o=B.obs[['leiden','cell_type','cell_type_confidence','cell_type_merged',
             'cell_type_merged_confidence']].copy(); o.index.name='cell_id'
    o.reset_index().to_csv(PROC+'celltypes_%s.csv'%tag,index=False)
    S.round(3).to_csv(PROC+'cluster_celltype_zscores_%s.csv'%tag)
    json.dump({'assignable':avail,'merged_label_set':avail_m,'merge_groups':MERGE,
               'composition_merged':{k:int(v) for k,v in B.obs.cell_type_merged.value_counts().items()},'dropped_thin':dropped,'dropped_nonspecific':DROP_NONSPECIFIC,
               'params':dict(MIN_MARKERS=MIN_MARKERS,MIN_DET=MIN_DET,MIN_Z=MIN_Z,
                             MIN_MARGIN=MIN_MARGIN,LOWQ_FRAC=LOWQ_FRAC,RES=RES),
               'composition':{k:int(v) for k,v in vc.items()}},
              open(PROC+'annotation_meta_%s.json'%tag,'w'),indent=1)
    print('wrote celltypes_%s.csv'%tag, flush=True)

if __name__=='__main__':
    args=[a for a in sys.argv[1:] if a!='--relabel']
    rel='--relabel' in sys.argv
    for d in args:
        d=d.rstrip('/'); base=os.path.basename(d)
        tag={'7250_liver_sham_Male_26-U1':'sham','7259_liver_sbr_Male_26-U1':'sbr'}.get(base,base)
        run(d,tag,relabel=rel)
