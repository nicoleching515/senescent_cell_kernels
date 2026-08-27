#!/usr/bin/env python3
"""T4 -- characterise WHAT each senescence caller is actually picking up.

Phase 2 established the callers agree at or below chance (Jaccard ratio 0.60-1.66x, all
|rho|<0.03).  Disagreement alone is not a finding; WHICH AXIS each caller loads on is.

Three candidate axes, all measurable:
  (a) TECHNICAL DEPTH   total transcript counts / genes detected / cell area.  A score that is
      really a depth score will call the biggest cells.
  (b) CELL TYPE         global top-5% composition vs background composition.
  (c) ZONATION          within hepatocytes, periportal / midzonal / pericentral.

Thresholding note: Tier A's stored sender_flag_p* are WITHIN-cell-type by construction and are
therefore flat across cell types by definition -- useless for (b).  Every caller is re-thresholded
here at a GLOBAL top-5% so the comparison is apples to apples; the within-type version is also
reported.
"""
import sys, os, json, numpy as np, pandas as pd, warnings
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')
PROC='/workspace/data/processed/'; RAW='/workspace/data/raw/'; OUT='/workspace/results/phase3/'
SAMP={'sham':'7250_liver_sham_Male_26-U1','sbr':'7259_liver_sbr_Male_26-U1'}
# inverse of SAMP: full section directory name -> the preserved short filename used by
# data/processed/deepscence_{sham,sbr}.csv.  caller_disagree_all.py:_load() needs this so a
# spec given by section name still finds the two DeepScence score files that were run under
# the short names.  Identity for any tag not listed.
DS_ALIAS={v:k for k,v in SAMP.items()}
EXCL={'Low_quality','Unknown'}

def run(tag, arm):
    sen=pd.read_csv(PROC+'senders_%s.csv'%tag).set_index('cell_id')
    ana=pd.read_csv(PROC+'anatomy_%s.csv'%tag).set_index('cell_id')
    cells=pd.read_parquet(RAW+SAMP.get(tag,tag)+'/cells.parquet').set_index('cell_id')
    ds=PROC+'deepscence_%s.csv'%tag
    if os.path.exists(ds): sen=sen.join(pd.read_csv(ds).set_index('cell_id'),rsuffix='_ds')
    df=sen.join(ana[['zonation_score','compartment_label']]).join(
        cells[['transcript_counts','cell_area','nucleus_area','segmentation_method']])
    df['n_genes']=np.nan
    df=df[~df.cell_type.isin(EXCL)].copy()
    df['log_counts']=np.log1p(df.transcript_counts)
    SCORES=[c for c in ['tierA_score','senepy_score','deepscence_score','cdkn1a_counts'] if c in df]
    res={'section':tag,'arm':arm,'n':int(len(df))}

    # ---- (a) technical-depth loading ----
    tech=['transcript_counts','cell_area','nucleus_area']
    rows=[]
    for s in SCORES:
        v=df[s]; ok=v.notna()
        r=dict(section=tag,arm=arm,score=s,n_scored=int(ok.sum()))
        for t in tech:
            r['rho_'+t]=round(float(spearmanr(v[ok],df[t][ok]).statistic),4)
        r['rho_zonation']=round(float(spearmanr(v[ok],df.zonation_score[ok]).statistic),4)
        rows.append(r)
    tech_tab=pd.DataFrame(rows)

    # ---- global top-5% calls ----
    flags={}
    for s in SCORES:
        v=df[s].to_numpy(float); ok=np.isfinite(v)
        thr=np.nanpercentile(v[ok],95)
        f=np.zeros(len(df),bool); f[ok]=v[ok]>thr
        flags['G_'+s]=f
        # within-cell-type
        fw=np.zeros(len(df),bool)
        for c in df.cell_type.unique():
            m=(df.cell_type==c).to_numpy()&ok
            if m.sum()<50: continue
            fw[m]=v[m]>np.nanpercentile(v[m],95)
        flags['W_'+s]=fw
    F=pd.DataFrame(flags,index=df.index)

    # ---- (b) cell-type composition of the calls ----
    bg=df.cell_type.value_counts(normalize=True)
    comp=[]
    for k in [c for c in F.columns if c.startswith('G_')]:
        sub=df.cell_type[F[k].to_numpy()].value_counts(normalize=True)
        for c in bg.index:
            comp.append(dict(section=tag,arm=arm,caller=k[2:],cell_type=c,
                             bg_pct=round(100*bg[c],3),call_pct=round(100*sub.get(c,0.0),3),
                             enrichment=round(sub.get(c,0.0)/bg[c],3),
                             n_called=int((F[k].to_numpy()&(df.cell_type==c).to_numpy()).sum())))
    comp=pd.DataFrame(comp)

    # ---- (c) depth strata + zonation strata ----
    df['depth_q']=pd.qcut(df.transcript_counts,5,labels=['Q1_low','Q2','Q3','Q4','Q5_high'],duplicates='drop')
    strat=[]
    for k in [c for c in F.columns if c.startswith('G_')]:
        f=F[k].to_numpy()
        for q,gsub in df.groupby('depth_q',observed=True):
            m=(df.depth_q==q).to_numpy()
            strat.append(dict(section=tag,arm=arm,caller=k[2:],stratum_kind='depth_quintile',
                stratum=str(q),bg_pct=round(100*m.mean(),2),
                call_pct=round(100*(f&m).sum()/max(f.sum(),1),2),
                enrichment=round(((f&m).sum()/max(f.sum(),1))/m.mean(),3)))
        hep=(df.cell_type=='Hepatocytes').to_numpy()
        for q in ['periportal','midzonal','pericentral']:
            m=hep&(df.compartment_label==q).to_numpy()
            fh=f&hep
            strat.append(dict(section=tag,arm=arm,caller=k[2:],stratum_kind='zonation_hepatocyte',
                stratum=q,bg_pct=round(100*m.sum()/max(hep.sum(),1),2),
                call_pct=round(100*(fh&m).sum()/max(fh.sum(),1),2),
                enrichment=round(((fh&m).sum()/max(fh.sum(),1))/max(m.sum()/max(hep.sum(),1),1e-9),3)))
    strat=pd.DataFrame(strat)

    # ---- pairwise agreement, global and within-type ----
    pair=[]
    for pre,lab in [('G_','global_top5'),('W_','within_type_top5')]:
        ks=[c for c in F.columns if c.startswith(pre)]
        for i in range(len(ks)):
            for j in range(i+1,len(ks)):
                a=F[ks[i]].to_numpy(); b=F[ks[j]].to_numpy()
                ok=np.isfinite(df[ks[i][2:]].to_numpy(float))&np.isfinite(df[ks[j][2:]].to_numpy(float))
                a=a[ok]; b=b[ok]
                inter=(a&b).sum(); uni=(a|b).sum()
                jac=inter/uni if uni else np.nan
                chance=a.mean()*b.mean()*len(a)/max(uni,1)
                pair.append(dict(section=tag,arm=arm,threshold=lab,A=ks[i][2:],B=ks[j][2:],
                    n=int(ok.sum()),n_A=int(a.sum()),n_B=int(b.sum()),n_both=int(inter),
                    jaccard=round(float(jac),5),chance_jaccard=round(float(chance),5),
                    ratio=round(float(jac/chance),3) if chance else None,
                    spearman=round(float(spearmanr(df[ks[i][2:]][ok],df[ks[j][2:]][ok]).statistic),4)))
    return tech_tab,comp,strat,pd.DataFrame(pair),res

if __name__=='__main__':
    T=[];C=[];S=[];P=[];R=[]
    for spec in sys.argv[1:]:
        tag,arm=spec.split('='); print('...',tag,flush=True)
        t,c,s,p,r=run(tag,arm); T.append(t);C.append(c);S.append(s);P.append(p);R.append(r)
    pd.concat(T).to_csv(OUT+'caller_technical_loading.csv',index=False)
    pd.concat(C).to_csv(OUT+'caller_celltype_composition.csv',index=False)
    pd.concat(S).to_csv(OUT+'caller_strata.csv',index=False)
    pd.concat(P).to_csv(OUT+'caller_pairwise_agreement.csv',index=False)
    print('\n=== (a) technical / zonation loading (Spearman rho) ===')
    print(pd.concat(T).to_string(index=False))
    print('\n=== (b) cell-type enrichment of global top-5% calls ===')
    cc=pd.concat(C)
    print(cc.pivot_table(index=['section','cell_type'],columns='caller',values='enrichment').to_string())
    print('\n=== (c) strata enrichment ===')
    ss=pd.concat(S)
    print(ss.pivot_table(index=['section','stratum_kind','stratum'],columns='caller',values='enrichment').to_string())
    print('\n=== pairwise agreement ===')
    print(pd.concat(P).to_string(index=False))
