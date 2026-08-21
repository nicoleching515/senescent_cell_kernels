#!/usr/bin/env python3
"""T3 / Deliverable 5 -- Tier C ligand-receptor plausibility on the GSE310392 mouse 5K panel.

Three questions, answered with measurements rather than assertion:
  Q1 EXPRESSION   Which Tier C ligands are actually detected in SENDERS, and which cognate
                  receptors in candidate RECEIVERS, per cell type, in the SBR arm?
  Q2 ENRICHMENT   Is the ligand enriched in senders over non-senders of the SAME cell type?
                  (a ligand expressed equally everywhere carries no sender-specific signal)
  Q3 IDENTIFIABILITY  For each ligand, what range of lambda could even be resolved?  The
                  binding quantity is the spatial density of ligand+ cells: the distance-to-
                  nearest-ligand+-cell distribution sets the dynamic range of the regressor.
                  If Il1a+ cells are 30x rarer than Ccl2+ cells, Il1a's d-distribution is
                  shifted far right and any fitted lambda will be larger for DENSITY reasons,
                  in the OPPOSITE direction to the Section 9 biological prediction.
"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/code')
import sasp_io, warnings, scanpy as sc
from scipy.spatial import cKDTree
warnings.filterwarnings('ignore'); sc.settings.verbosity=0
PROC='/workspace/data/processed/'; GS='/workspace/genesets/'; OUT='/workspace/results/phase3/'
def gl(n): return [l.strip() for l in open(GS+n+'.txt') if l.strip()]
LIG, REC = gl('C_ligands'), gl('C_receptors')
# Section 9 Tier C pairs, mouse. Cxcl8 has no mouse ortholog; Cxcl1/2/5 are the analogues.
PAIRS = [('Il1a',['Il1r1','Il1rap'],'membrane-bound; predicted SHORTEST lambda'),
         ('Tgfb1',['Tgfbr1','Tgfbr2'],'latent/ECM-tethered; short'),
         ('Tnf',['Tnfrsf1a','Tnfrsf1b'],'membrane + shed soluble; short-intermediate'),
         ('Il6',['Il6ra','Il6st'],'freely secreted; intermediate'),
         ('Igfbp3',[],'secreted carrier; no cognate receptor on panel'),
         ('Gdf15',[],'secreted; GFRAL off-panel'),
         ('Timp1',['Cd47'],'secreted; Cd63 off-panel'),
         ('Mmp3',[],'secreted protease; no receptor'),
         ('Thbs1',['Cd47','Sdc1'],'matricellular; short'),
         ('Ccl2',['Ccr2','Ackr3'],'small diffusible chemokine; predicted LONG'),
         ('Cxcl12',['Cxcr4','Ackr3','Dpp4'],'small chemokine; Dpp4 cleaves it -> range-limiting'),
         ('Cxcl1',['Cxcr2'],'CXCL8 analogue; predicted LONGEST'),
         ('Cxcl2',['Cxcr2'],'CXCL8 analogue; predicted LONGEST'),
         ('Cxcl5',['Cxcr2'],'CXCL8 analogue; predicted LONGEST')]
EXCL={'Low_quality','Unknown'}

def run(tag, arm):
    B=sasp_io.load(tag)
    sen=pd.read_csv(PROC+'senders_%s.csv'%tag).set_index('cell_id')
    ana=pd.read_csv(PROC+'anatomy_%s.csv'%tag).set_index('cell_id')
    B=B[B.obs_names.isin(sen.index)].copy()
    ct=sen['cell_type'].reindex(B.obs_names).astype(str).values
    cdkpos=sen['cdkn1a_pos'].reindex(B.obs_names).to_numpy().astype(bool)
    p95=sen['sender_flag_p95'].reindex(B.obs_names).to_numpy().astype(bool)
    xy=B.obs[['x_um','y_um']].to_numpy()
    C=np.asarray(B.layers['counts'].todense()) if False else B.layers['counts'].tocsc()
    vn=list(B.var_names); n=B.n_obs
    def counts(g):
        if g not in vn: return None
        return np.asarray(C[:,vn.index(g)].todense()).ravel()
    keep=~np.isin(ct,list(EXCL))
    rows_expr=[]; rows_id=[]
    # ---- Q1/Q2 : per-cell-type detection, and sender enrichment ----
    for g in sorted(set(LIG+REC)):
        v=counts(g)
        if v is None:
            rows_expr.append(dict(section=tag,arm=arm,gene=g,role='ligand' if g in LIG else 'receptor',
                                  cell_type='(ANY)',n=0,pct_pos=np.nan,mean_counts=np.nan,
                                  pct_pos_senders=np.nan,pct_pos_nonsenders=np.nan,enrich=np.nan,
                                  on_panel=False)); continue
        pos=v>0
        for c in sorted(set(ct[keep])):
            m=(ct==c)&keep
            if m.sum()<200: continue
            ms=m&cdkpos; mn=m&~cdkpos
            ps=100*pos[ms].mean() if ms.sum()>=30 else np.nan
            pn=100*pos[mn].mean() if mn.sum()>=30 else np.nan
            rows_expr.append(dict(section=tag,arm=arm,gene=g,
                role=('ligand' if g in LIG else '')+('receptor' if g in REC else ''),
                cell_type=c,n=int(m.sum()),pct_pos=round(100*pos[m].mean(),3),
                mean_counts=round(float(v[m].mean()),4),
                pct_pos_senders=None if np.isnan(ps) else round(ps,3),
                pct_pos_nonsenders=None if np.isnan(pn) else round(pn,3),
                enrich=None if (np.isnan(ps) or np.isnan(pn) or pn==0) else round(ps/pn,3),
                on_panel=True))
    # ---- Q3 : identifiability of lambda per ligand, from ligand+ SENDER cell density ----
    area=(xy[:,0].max()-xy[:,0].min())*(xy[:,1].max()-xy[:,1].min())/1e6  # mm^2
    tree_all=cKDTree(xy)
    med_nn=np.median(tree_all.query(xy,k=2)[0][:,1])
    for g,recs,note in PAIRS:
        v=counts(g)
        if v is None:
            rows_id.append(dict(section=tag,arm=arm,ligand=g,on_panel=False,note=note)); continue
        # "ligand-expressing sender" = Cdkn1a+ cell that also expresses the ligand
        src=np.where((v>0)&cdkpos&keep)[0]
        srcL=np.where((v>0)&keep)[0]
        d=cKDTree(xy[src]).query(xy,k=1)[0] if len(src)>=5 else np.full(n,np.nan)
        dL=cKDTree(xy[srcL]).query(xy,k=1)[0] if len(srcL)>=5 else np.full(n,np.nan)
        rows_id.append(dict(section=tag,arm=arm,ligand=g,on_panel=True,
            n_ligand_pos=int(len(srcL)), pct_ligand_pos=round(100*len(srcL)/keep.sum(),3),
            n_ligand_pos_senders=int(len(src)),
            dens_lig_per_mm2=round(len(srcL)/area,1),
            dens_ligsender_per_mm2=round(len(src)/area,2),
            d_med_to_ligand_um=round(float(np.nanmedian(dL)),1),
            d_p10_to_ligand_um=round(float(np.nanpercentile(dL,10)),1),
            d_med_to_ligsender_um=round(float(np.nanmedian(d)),1) if len(src)>=5 else None,
            d_p10_to_ligsender_um=round(float(np.nanpercentile(d,10)),1) if len(src)>=5 else None,
            recept_on_panel=';'.join(recs), note=note))
    meta=dict(section=tag,arm=arm,n_cells=int(n),n_used=int(keep.sum()),
              area_mm2=round(area,2),median_nn_um=round(float(med_nn),2),
              n_cdkn1a_pos=int((cdkpos&keep).sum()))
    return pd.DataFrame(rows_expr), pd.DataFrame(rows_id), meta

if __name__=='__main__':
    E=[];I=[];M=[]
    for spec in sys.argv[1:]:
        tag,arm=spec.split('=')
        print('...',tag,flush=True)
        e,i,m=run(tag,arm); E.append(e);I.append(i);M.append(m)
    E=pd.concat(E);I=pd.concat(I)
    E.to_csv(OUT+'tierC_expression_by_celltype.csv',index=False)
    I.to_csv(OUT+'tierC_ligand_identifiability.csv',index=False)
    json.dump(M,open(OUT+'tierC_meta.json','w'),indent=1)
    print(json.dumps(M,indent=1))
    print('\n=== Q3 ligand identifiability (SBR) ===')
    print(I[I.on_panel==True][['section','ligand','pct_ligand_pos','dens_lig_per_mm2',
         'd_med_to_ligand_um','n_ligand_pos_senders','d_med_to_ligsender_um','note']].to_string(index=False))
