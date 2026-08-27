#!/usr/bin/env python3
"""Phase 2 D-B (anatomy), D-C (sender calls), D-D (module scores). Runs per sample tag."""
import sys, os, json, glob, warnings, numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0,'/workspace/code')
import sasp_io
from scipy.spatial import cKDTree
from scipy import ndimage
from sklearn.cluster import DBSCAN
warnings.filterwarnings('ignore'); sc.settings.n_jobs=48; sc.settings.verbosity=0
PROC='/workspace/data/processed/'; GS='/workspace/genesets/'
def gl(n): return [l.strip() for l in open(GS+n+'.txt') if l.strip()]
PC,PP = gl('D_zonation_pericentral'), gl('D_zonation_periportal')
A_STRICT = gl('A_SENDER_FINAL_strict')
BMODS = {os.path.basename(p)[2:-4]:[l.strip() for l in open(p) if l.strip()]
         for p in sorted(glob.glob(GS+'B_*.txt'))}
# Phase 8 / D1: the seven per-module Tier A sensitivity sender sets.  One set
# per Tier B module, scored and thresholded exactly as A_SENDER_FINAL_strict
# is, so `tierApm_pNN` differs from `tierA_pNN` only in which genes define the
# sender score.  The names must be the Tier B module list, or a downstream
# `sender_mask(call, module=...)` would silently look up a column that is not
# there.
AMODS = {os.path.basename(p)[len('A_sender_for_'):-4]:
         [l.strip() for l in open(p) if l.strip()]
         for p in sorted(glob.glob(GS+'A_sender_for_*.txt'))}
assert sorted(AMODS) == sorted(BMODS), (
    'A_sender_for_* sets do not match the Tier B module list: %s vs %s'
    % (sorted(AMODS), sorted(BMODS)))
GRID_UM=25.0
EXCLUDE_FROM_STRATA = {'Low_quality','Unknown'}   # not receiver cell types

for tag in sys.argv[1:]:
    print('\n'+'='*90); print('SAMPLE',tag); print('='*90, flush=True)
    B = sasp_io.load(tag)
    ctdf = pd.read_csv(PROC+'celltypes_%s.csv'%tag).set_index('cell_id')
    B = B[B.obs_names.isin(ctdf.index)].copy()
    B.obs['cell_type']=ctdf['cell_type'].reindex(B.obs_names).values
    ct=B.obs['cell_type'].astype(str); xy=B.obs[['x_um','y_um']].to_numpy(); n=B.n_obs
    print('cells %d ; types %d'%(n,ct.nunique()), flush=True)

    # ---------------- D-B anatomy ----------------
    sc.tl.score_genes(B,[g for g in PC if g in B.var_names],score_name='_pc',ctrl_size=200)
    sc.tl.score_genes(B,[g for g in PP if g in B.var_names],score_name='_pp',ctrl_size=200)
    zon=(B.obs['_pc']-B.obs['_pp']).to_numpy()
    hep=(ct=='Hepatocytes').to_numpy()
    zon_z=(zon-zon[hep].mean())/zon[hep].std()      # standardised on hepatocytes
    q=np.quantile(zon_z[hep],[1/3,2/3])
    comp=np.where(zon_z<=q[0],'periportal',np.where(zon_z>=q[1],'pericentral','midzonal'))
    comp=np.where(hep,comp,'non_hepatocyte')
    # distance to tissue boundary: occupancy grid -> close/fill/open -> EDT
    x0,y0=xy.min(0)-GRID_UM; nx=int((xy[:,0].max()+GRID_UM-x0)/GRID_UM)+1; ny=int((xy[:,1].max()+GRID_UM-y0)/GRID_UM)+1
    ix=((xy[:,0]-x0)/GRID_UM).astype(int); iy=((xy[:,1]-y0)/GRID_UM).astype(int)
    occ=np.zeros((nx,ny),bool); occ[ix,iy]=True
    occ=ndimage.binary_fill_holes(ndimage.binary_closing(occ,np.ones((3,3))))
    occ=ndimage.binary_opening(occ,np.ones((3,3)))
    dist_bound=(ndimage.distance_transform_edt(occ)*GRID_UM)[ix,iy]
    # distance to portal triad, proxied by bile ducts (cholangiocyte foci)
    chol=np.where(ct=='Biliary/ductular')[0]; dist_pt=np.full(n,np.nan); nf=0
    if len(chol)>=30:
        lab=DBSCAN(eps=30,min_samples=10).fit(xy[chol]).labels_
        if (lab>=0).any():
            cents=np.array([xy[chol][lab==l].mean(0) for l in np.unique(lab[lab>=0])]); nf=len(cents)
            dist_pt=cKDTree(cents).query(xy,k=1)[0]
    print('boundary dist median %.1f um | portal-triad foci %d | dist median %.1f um'
          %(np.median(dist_bound),nf,np.nanmedian(dist_pt) if nf else np.nan))
    # validation: does the zonation score behave?
    r=np.corrcoef(zon_z[hep],dist_pt[hep])[0,1] if nf else np.nan
    print('corr(zonation, dist_to_portal_triad) within hepatocytes = %+.3f  (expect POSITIVE:'
          ' pericentral = far from portal triad)'%r)
    pd.DataFrame({'cell_id':B.obs_names,'cell_type':ct.values,'zonation_score':np.round(zon_z,4),
                  'compartment_label':comp,'dist_to_boundary_um':np.round(dist_bound,2),
                  'dist_to_portal_triad_um':np.round(dist_pt,2)}).to_csv(PROC+'anatomy_%s.csv'%tag,index=False)

    # ---------------- D-C senders ----------------
    cdk=np.asarray(B.layers['counts'][:,B.var_names.get_loc('Cdkn1a')].todense()).ravel()
    sc.tl.score_genes(B,[g for g in A_STRICT if g in B.var_names],score_name='_tierA',ctrl_size=200)
    tierA=B.obs['_tierA'].to_numpy()
    import senepy
    H=senepy.load_hubs(species='Mouse'); HC=senepy.load_hubs(species='Mouse',sig_type='cell_type')
    HUBMAP={'Hepatocytes':('Liver','hepatocyte',1),
            'LSECs':('Liver','endothelial cell of hepatic sinusoid',0),
            'Central venous LSECs':('Liver','endothelial cell of hepatic sinusoid',0),
            'Kupffer cells':('Liver','Kupffer cell',0),
            'Inflammatory macs':('Myeloid','macrophage'),
            'T/NK cells':('Lymphoid','T cell'),'B-cells':('Lymphoid','B cell')}
    sp=np.full(n,np.nan); used={}
    for cell_type,key in HUBMAP.items():
        m=(ct==cell_type).to_numpy()
        if m.sum()<50: continue
        hub=H.hubs.get(key) if len(key)==3 else HC.hubs.get(key)
        if hub is None: continue
        onp=[(g,i) for g,i in hub if g in B.var_names]
        if len(onp)<10: print('  senepy %s: %d hub genes on panel, skipped'%(cell_type,len(onp))); continue
        used[cell_type]=(str(key),len(hub),len(onp))
        sp[m]=senepy.score_hub(B[m].copy(),onp,verbose=False)
    print('SenePy hubs used:')
    for k,v in used.items(): print('   %-24s %s hub=%d on-panel=%d'%(k,v[0],v[1],v[2]))
    sen=pd.DataFrame({'cell_id':B.obs_names,'cell_type':ct.values,'cdkn1a_counts':cdk,
                      'cdkn1a_pos':(cdk>0).astype(int),'tierA_score':np.round(tierA,5),
                      'senepy_score':np.round(sp,5)})
    for q_ in [90,95,99]:
        f=np.zeros(n,int)
        for c in sen.cell_type.unique():
            if c in EXCLUDE_FROM_STRATA: continue
            m=(sen.cell_type==c).to_numpy()
            if m.sum()<20: continue
            f[m]=(tierA[m]>np.percentile(tierA[m],q_)).astype(int)
        sen['sender_flag_p%d'%q_]=f
    # --- Phase 8 / D1: the per-module Tier A sensitivity sender sets --------
    # Appended AFTER the primary columns and BEFORE the DeepScence merge, so
    # neither the existing columns nor their order moves and no row alignment
    # can drift.  Same scorer, same ctrl_size, same within-cell-type strict
    # percentile rule, same >= 20-cell stratum floor as the primary block.
    for name in sorted(AMODS):
        on=[g for g in AMODS[name] if g in B.var_names]
        sc.tl.score_genes(B,on,score_name='_tierApm',ctrl_size=200)
        v=B.obs['_tierApm'].to_numpy()
        sen['tierA_%s_score'%name]=np.round(v,5)
        for q_ in [90,95,99]:
            f=np.zeros(n,int)
            for c in sen.cell_type.unique():
                if c in EXCLUDE_FROM_STRATA: continue
                m=(sen.cell_type==c).to_numpy()
                if m.sum()<20: continue
                f[m]=(v[m]>np.percentile(v[m],q_)).astype(int)
            sen['sender_flag_%s_p%d'%(name,q_)]=f
        print('  A_sender_for_%-22s %3d genes, %3d on panel, p95 senders %6d'
              %(name,len(AMODS[name]),len(on),
                int(sen['sender_flag_%s_p95'%name].sum())))
    ds=PROC+'deepscence_%s.csv'%tag
    if os.path.exists(ds):
        sen=sen.merge(pd.read_csv(ds),on='cell_id',how='left')
    sen.to_csv(PROC+'senders_%s.csv'%tag,index=False)

    # Test 3 prevalence table, per cell type
    print('\n--- Sec 8 Test 3: sender prevalence per cell type ---')
    print('%-26s %8s %9s %9s %9s'%('cell_type','n','Cdkn1a+%','p95 flag%','in 1-20%?'))
    rows=[]
    for c,gsub in sen.groupby('cell_type'):
        if c in EXCLUDE_FROM_STRATA: continue
        cp=100*gsub.cdkn1a_pos.mean(); p95=100*gsub.sender_flag_p95.mean()
        ok='YES' if 1<=cp<=20 else ('LOW' if cp<1 else 'HIGH')
        rows.append((c,len(gsub),cp,p95,ok))
        print('%-26s %8d %9.2f %9.2f %9s'%(c,len(gsub),cp,p95,ok))
    pd.DataFrame(rows,columns=['cell_type','n','cdkn1a_pos_pct','sender_p95_pct','cdkn1a_in_1_20_band']
                 ).to_csv(PROC+'test3_prevalence_%s.csv'%tag,index=False)

    # ---------------- D-D module scores ----------------
    mod=pd.DataFrame({'cell_id':B.obs_names})
    for name,genes in BMODS.items():
        on=[g for g in genes if g in B.var_names]
        sc.tl.score_genes(B,on,score_name='_m',ctrl_size=max(200,len(on)*5))
        mod[name]=np.round(B.obs['_m'].to_numpy(),5)
    mod.to_csv(PROC+'modules_%s.csv'%tag,index=False)
    print('\nwrote anatomy/senders/modules/test3 for %s'%tag, flush=True)
    del B
print('\nDONE')
