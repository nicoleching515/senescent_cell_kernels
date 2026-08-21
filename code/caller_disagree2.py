#!/usr/bin/env python3
"""T4 part 2 -- mechanism of the caller disagreement.

Adds to caller_disagree.py:
 (1) a NAIVE CoreScence score computed here (mean z of up-genes minus mean z of down-genes,
     occurrence>=5, mouse orthologs), to separate "the gene set loads on depth" from
     "the DeepScence autoencoder loads on depth".  We run DeepScence with denoise=False
     (forced; the DCA stack will not install), and DCA denoising is precisely the step that
     would normalise depth -- so this distinction matters and must be reported.
 (2) DEPTH-STRATIFIED agreement: if the disagreement is only a depth artefact, callers should
     agree better inside a depth quintile.
 (3) partial Spearman between callers after regressing out log transcript counts.
"""
import sys, os, csv, json, numpy as np, pandas as pd, warnings
sys.path.insert(0,'/workspace/code')
import sasp_io, scanpy as sc
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings('ignore'); sc.settings.verbosity=0
PROC='/workspace/data/processed/'; RAW='/workspace/data/raw/'; OUT='/workspace/results/phase3/'
EXCL={'Low_quality','Unknown'}

core=pd.read_csv('/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv')
core['occurrence']=pd.to_numeric(core.occurrence,errors='coerce')
core=core[core.occurrence>=5]
orth={r['mouse_symbol']:r['human_symbol'] for r in
      csv.DictReader(open('/workspace/genesets/mouse_human_orthologs_MGI.csv'))}
h2m={}
for m,h in orth.items(): h2m.setdefault(h,m)
UP=[h2m[g] for g in core[core.direction=='up'].gene_symbol if g in h2m]
DN=[h2m[g] for g in core[core.direction=='down'].gene_symbol if g in h2m]

def partial_spearman(x,y,z):
    rx,ry,rz=(rankdata(v) for v in (x,y,z))
    def resid(a,b):
        b=np.c_[np.ones(len(b)),b]; return a-b@np.linalg.lstsq(b,a,rcond=None)[0]
    return float(np.corrcoef(resid(rx,rz),resid(ry,rz))[0,1])

def run(tag,arm):
    B=sasp_io.load(tag)
    sen=pd.read_csv(PROC+'senders_%s.csv'%tag).set_index('cell_id')
    ds=PROC+'deepscence_%s.csv'%tag
    if os.path.exists(ds): sen=sen.join(pd.read_csv(ds).set_index('cell_id'))
    B=B[B.obs_names.isin(sen.index)].copy()
    keep=~sen['cell_type'].reindex(B.obs_names).isin(EXCL).to_numpy()
    up=[g for g in UP if g in B.var_names]; dn=[g for g in DN if g in B.var_names]
    sc.tl.score_genes(B,up,score_name='_u',ctrl_size=200)
    sc.tl.score_genes(B,dn,score_name='_d',ctrl_size=200)
    naive=(B.obs['_u']-B.obs['_d']).to_numpy()
    df=pd.DataFrame(dict(cell_id=B.obs_names,
        cell_type=sen['cell_type'].reindex(B.obs_names).values,
        counts=B.obs['transcript_counts'].to_numpy(float),
        cell_area=B.obs['cell_area'].to_numpy(float),
        nucleus_area=B.obs['nucleus_area'].to_numpy(float),
        naive_corescence=naive,
        tierA_score=sen['tierA_score'].reindex(B.obs_names).to_numpy(float),
        senepy_score=sen['senepy_score'].reindex(B.obs_names).to_numpy(float),
        deepscence_score=(sen['deepscence_score'].reindex(B.obs_names).to_numpy(float)
                          if 'deepscence_score' in sen else np.nan),
        cdkn1a_counts=sen['cdkn1a_counts'].reindex(B.obs_names).to_numpy(float)))[keep]
    df['n_genes']=np.asarray((B.layers['counts'][keep]>0).sum(1)).ravel()
    df['lc']=np.log1p(df.counts)
    SC=[c for c in ['tierA_score','senepy_score','deepscence_score','naive_corescence','cdkn1a_counts']
        if df[c].notna().any()]
    print('CoreScence occ>=5 : %d up / %d down human ; %d up / %d down on mouse panel'
          %((core.direction=='up').sum(),(core.direction=='down').sum(),len(up),len(dn)))
    rows=[]
    for s in SC:
        ok=df[s].notna()&df.nucleus_area.notna()
        ok2=df[s].notna()
        rows.append(dict(section=tag,arm=arm,score=s,n=int(ok2.sum()),
          rho_counts=round(spearmanr(df[s][ok2],df.counts[ok2]).statistic,4),
          rho_n_genes=round(spearmanr(df[s][ok2],df.n_genes[ok2]).statistic,4),
          rho_cell_area=round(spearmanr(df[s][ok2],df.cell_area[ok2]).statistic,4),
          rho_nucleus_area=round(spearmanr(df[s][ok],df.nucleus_area[ok]).statistic,4)))
    T=pd.DataFrame(rows)
    # depth-stratified agreement
    df['q']=pd.qcut(df.counts,5,labels=['Q1_low','Q2','Q3','Q4','Q5_high'],duplicates='drop')
    st=[]
    for i in range(len(SC)):
        for j in range(i+1,len(SC)):
            a,b=SC[i],SC[j]
            ok=df[a].notna()&df[b].notna()
            st.append(dict(section=tag,arm=arm,A=a,B=b,stratum='ALL',
              spearman=round(spearmanr(df[a][ok],df[b][ok]).statistic,4),
              partial_spearman_given_logcounts=round(partial_spearman(
                  df[a][ok].values,df[b][ok].values,df.lc[ok].values),4),
              jaccard_ratio=jr(df,a,b,ok)))
            for q in df.q.cat.categories:
                m=ok&(df.q==q)
                st.append(dict(section=tag,arm=arm,A=a,B=b,stratum=str(q),
                  spearman=round(spearmanr(df[a][m],df[b][m]).statistic,4),
                  partial_spearman_given_logcounts=None,
                  jaccard_ratio=jr(df,a,b,m)))
    return T,pd.DataFrame(st),df

def jr(df,a,b,m):
    x=df[a][m].to_numpy(); y=df[b][m].to_numpy()
    if len(x)<200: return None
    fa=x>np.percentile(x,95); fb=y>np.percentile(y,95)
    uni=(fa|fb).sum()
    if uni==0: return None
    jac=(fa&fb).sum()/uni; ch=fa.mean()*fb.mean()*len(x)/uni
    return round(float(jac/ch),3) if ch>0 else None

if __name__=='__main__':
    T=[];S=[]
    for spec in sys.argv[1:]:
        tag,arm=spec.split('='); print('...',tag,flush=True)
        t,s,df=run(tag,arm); T.append(t);S.append(s)
    T=pd.concat(T); S=pd.concat(S)
    T.to_csv(OUT+'caller_technical_loading2.csv',index=False)
    S.to_csv(OUT+'caller_depth_stratified_agreement.csv',index=False)
    print('\n=== technical loading (Spearman rho) ===');print(T.to_string(index=False))
    print('\n=== depth-stratified agreement ===')
    print(S.pivot_table(index=['section','A','B'],columns='stratum',values='jaccard_ratio').to_string())
    print('\n=== spearman, all vs partial|logcounts ===')
    print(S[S.stratum=='ALL'].to_string(index=False))
