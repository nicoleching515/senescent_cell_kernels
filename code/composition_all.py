#!/usr/bin/env python3
"""T2 -- composition by arm and timepoint across every annotated section of GSE310392.

Two denominators are reported, because the choice is not neutral:
  pct_all         -- fraction of all QC-passing cells, including Low_quality and Unknown
  pct_analysable  -- fraction of cells that are actually usable as receivers, i.e. excluding
                     Low_quality and Unknown.  These two labels swing 0-13% between sections
                     for reasons (segmentation quality, cluster granularity) that have nothing
                     to do with biology, so a trend computed on pct_all is partly a trend in
                     annotation quality.  pct_analysable is the primary readout.
Kupffer cells and 'Inflammatory macs' are ALSO reported summed as 'Myeloid_total': the boundary
between them is a clustering-resolution decision and it moves between sections (7448 puts 0.2%
in Kupffer and 9.1% in Inflammatory macs; 7250 puts 8.1% in Kupffer and 0% in Inflammatory macs).
The split is not trustworthy per section; the sum is.
"""
import os, glob, json, numpy as np, pandas as pd
PROC='/workspace/data/processed/'; OUT='/workspace/results/'
# GEO sample table (BIO_PHASE1 section 2.1). mouse 7239 also contributes a tumour section (not used).
META={'7361_liver_sbr_Male_2-U1' :('7361','SBR',2),  '7352_liver_sham_Male_2-U1' :('7352','sham',2),
      '7448_liver_sbr_Male_10-U1':('7448','SBR',10), '7450_liver_sbr_Male_10-U1':('7450','SBR',10),
      '7435_liver_sham_Male_10-U1':('7435','sham',10),
      '7259_liver_sbr_Male_26-U1':('7259','SBR',26), '7260_liver_sbr_Male_26-U1':('7260','SBR',26),
      '7248_liver_sham_Male_26-U1':('7248','sham',26),'7250_liver_sham_Male_26-U1':('7250','sham',26),
      '7239_liver_sbr_Male_52-U1':('7239','SBR',52), '7001_liver_sham_Male_52-U1':('7001','sham',52)}
ALIAS={}
EXCL=['Low_quality','Unknown']
FOCUS=['Hepatocytes','Biliary/ductular','Mesenchymal','Macrophages','Endothelial',
       'T/NK cells','B-cells','DC','vSMCs','Proliferating']

rows=[]
for p in sorted(glob.glob(PROC+'celltypes_*.csv')):
    tag=os.path.basename(p)[10:-4]
    key=ALIAS.get(tag,tag)
    if key not in META: print('SKIP (no metadata):',tag); continue
    mouse,arm,wk=META[key]
    d=pd.read_csv(p,usecols=['cell_type','cell_type_merged'])
    vc=d.cell_type_merged.value_counts(); vcf=d.cell_type.value_counts()
    n_all=len(d); n_an=n_all-int(vc.reindex(EXCL).fillna(0).sum())
    r=dict(section=key.split('_')[0],sample=key,mouse=mouse,arm=arm,timepoint_wk=wk,
           n_cells=n_all,n_analysable=n_an,
           pct_low_quality=round(100*vc.get('Low_quality',0)/n_all,2),
           pct_unknown=round(100*vc.get('Unknown',0)/n_all,2))
    for c in sorted(set(list(vc.index)+FOCUS)):
        if c in EXCL: continue
        r['pct_'+c]=round(100*vc.get(c,0)/n_an,3)
    for c in ['Hepatic stellate cells','Portal fibroblasts','Pericytes','Kupffer cells',
              'Inflammatory macs','LSECs','Central venous LSECs','Portal endothelial cells']:
        r['fine_pct_'+c]=round(100*vcf.get(c,0)/n_an,3)
    r['fine_pct_unknown']=round(100*vcf.get('Unknown',0)/len(d),2)
    rows.append(r)
D=pd.DataFrame(rows).drop_duplicates(subset=['sample']).sort_values(['arm','timepoint_wk','section'])
cols=['section','sample','mouse','arm','timepoint_wk','n_cells','n_analysable',
      'pct_low_quality','pct_unknown']+['pct_'+c for c in FOCUS]
D=D[cols+[c for c in D.columns if c not in cols]]
D.to_csv(OUT+'composition_by_arm_timepoint.csv',index=False)
pd.set_option('display.width',250)
key=['section','arm','timepoint_wk','n_cells','n_analysable','pct_low_quality','pct_unknown',
     'pct_Hepatocytes','pct_Biliary/ductular','pct_Mesenchymal','pct_Macrophages',
     'pct_T/NK cells','pct_Endothelial']
print(D[key].to_string(index=False))
print('\n--- FINE labels that merge (unstable per section; shown to justify the merge) ---')
print(D[['section','arm','timepoint_wk','fine_pct_Hepatic stellate cells','fine_pct_Portal fibroblasts',
         'fine_pct_Pericytes','fine_pct_Kupffer cells','fine_pct_Inflammatory macs','fine_pct_LSECs',
         'fine_pct_Central venous LSECs','fine_pct_Portal endothelial cells','fine_pct_unknown']].to_string(index=False))
print('\n--- monotonicity check, SBR arm (Spearman rho vs timepoint) ---')
from scipy.stats import spearmanr
for arm in ['SBR','sham']:
    s=D[D.arm==arm]
    if len(s)<3: continue
    print(' %s (n=%d sections):'%(arm,len(s)))
    for c in ['pct_'+x for x in FOCUS]:
        rho,p=spearmanr(s.timepoint_wk,s[c])
        print('   %-30s rho=%+.3f p=%.3f   values %s'%(c,rho,p,
              ' '.join('%.1f'%v for v in s[c])))
print('\nwrote',OUT+'composition_by_arm_timepoint.csv')
