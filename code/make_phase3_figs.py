#!/usr/bin/env python3
"""Phase 3 figures. Palette = the project's pre-validated categorical theme (sasp_palette.SERIES,
the dataviz-skill reference instance already validated in Phase 1); hues assigned in fixed order,
never cycled. node is unavailable in this container so validate_palette.js could not be re-run --
the palette is reused unchanged rather than eyeballed."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'/workspace/code')
import sasp_palette as P
P.apply_style(matplotlib)
FIG='/workspace/figures/'; RES='/workspace/results/'
ARM={'SBR':P.SERIES[1],'sham':P.SERIES[0]}

# ---------------------------------------------------------------- Fig A: composition
D=pd.read_csv(RES+'composition_by_arm_timepoint.csv')
panels=['Hepatocytes','Biliary/ductular','Mesenchymal','Macrophages','Endothelial','T/NK cells']
fig,axes=plt.subplots(2,3,figsize=(10.2,5.6),sharex=True)
for ax,c in zip(axes.ravel(),panels):
    for arm in ['sham','SBR']:
        s=D[D.arm==arm].sort_values('timepoint_wk')
        ax.plot(s.timepoint_wk,s['pct_'+c],'o',color=ARM[arm],ms=6,mec=P.SURFACE,mew=1.4,zorder=3)
        m=s.groupby('timepoint_wk')['pct_'+c].mean()
        ax.plot(m.index,m.values,'-',color=ARM[arm],lw=2,alpha=.85,zorder=2,label=arm)
    ax.set_title(c); ax.set_xticks([2,10,26,52]); ax.set_xlim(-2,58)
    ax.set_ylim(bottom=0)
for ax in axes[1]: ax.set_xlabel('weeks post-surgery')
for ax in axes[:,0]: ax.set_ylabel('% of analysable cells')
axes[0,0].legend(loc='upper right')
fig.suptitle('Cell-type composition of all 11 GSE310392 liver sections, by arm and timepoint\n'
             'points = sections (one mouse each); line = arm mean; merged label set',
             fontsize=10,y=1.005)
fig.tight_layout()
for _e in ('png','pdf'):
    fig.savefig(FIG+'fig_phase3_composition.'+_e,bbox_inches='tight')
plt.close(fig)
# the backing CSV: exactly the six panels' series, nothing else.  Before
# 2026-08-27 this file was committed without a producer (1351ce8).
D[['arm','timepoint_wk']+['pct_'+c for c in panels]].to_csv(
    FIG+'fig_phase3_composition_data.csv',index=False)

# ---------------------------------------------------------------- Fig B: caller depth bias
# BASIS.  This producer used to read caller_within_type_depth_bias.csv, which is
# byte-identical to results/phase3_pre_c6/ and disagrees with the committed
# artefact on all 10 tierA_score values (max |d| 0.069): the artefact plots the
# C6 table.  Re-running the old producer regressed the figure.  It now reads
# _2sec_c6.csv, which reproduces the artefact 40/40 exactly.
#
# COVERAGE.  The two-section base (7250 sham / 7259 SBR, caller_disagree.py:22)
# is 2 of the 11 sections, while the frozen configuration is 11-section
# coverage and caller_within_type_depth_bias_11sections.csv had been sitting
# unused.  Panel 3 now plots it, so the coverage leg is on the figure instead
# of only in a file.
B=pd.read_csv('/workspace/results/phase3/caller_within_type_depth_bias_2sec_c6.csv')
B11=pd.read_csv('/workspace/results/phase3/caller_within_type_depth_bias_11sections.csv')
NAME={'tierA_score':'Tier A score','senepy_score':'SenePy','deepscence_score':'DeepScence',
      'cdkn1a_counts':'Cdkn1a > 0'}
order=['tierA_score','deepscence_score','senepy_score','cdkn1a_counts']
col=dict(zip(order,P.SERIES[:4]))
QS=['Q1','Q2','Q3','Q4','Q5']
NSEC11=B11.section.nunique()
fig,axes=plt.subplots(1,3,figsize=(12.4,3.9),sharey=True)
for ax,(sec,lab) in zip(axes,[('sham','sham (7250, 26 wk)'),('sbr','SBR (7259, 26 wk)')]):
    s=B[B.section==sec]
    ax.axhline(1,color=P.MUTED,lw=1,ls='--',zorder=1)
    for k in order:
        q=s[s.caller==k].set_index('within_type_depth_quintile').reindex(QS).dropna(subset=['enrichment'])
        if not len(q): continue
        ax.plot(range(1,len(q)+1),q.enrichment,'-o',color=col[k],ms=6,mec=P.SURFACE,mew=1.4,
                label=NAME[k],zorder=3)
    ax.set_title(lab)
ax=axes[2]
ax.axhline(1,color=P.MUTED,lw=1,ls='--',zorder=1)
for k in order:
    for _sec,g in B11[B11.caller==k].groupby('section'):
        g=g.set_index('within_type_depth_quintile').reindex(QS)
        ax.plot(range(1,6),g.enrichment,'-',color=col[k],lw=0.7,alpha=0.30,zorder=2)
    med=B11[B11.caller==k].groupby('within_type_depth_quintile').enrichment.median().reindex(QS)
    ax.plot(range(1,6),med.values,'-o',color=col[k],ms=6,mec=P.SURFACE,mew=1.4,zorder=3)
ax.set_title('all %d sections (thin = each section,\nmarkers = median)'%NSEC11,fontsize=9)
for ax in axes:
    ax.set_yscale('log')
    ax.set_xticks([1,2,3,4,5]); ax.set_xticklabels(['Q1\nlow','Q2','Q3','Q4','Q5\nhigh'])
    ax.set_xlabel('transcript-count quintile, within cell type')
axes[0].set_ylabel('enrichment of top-5% calls\n(1 = no depth bias)')
h,l=axes[1].get_legend_handles_labels()          # order the legend by the value at Q5 so it
last={NAME[k]:B[(B.section=='sbr')&(B.caller==k)].set_index(                   # reads as a
      'within_type_depth_quintile').reindex(QS).enrichment.iloc[-1] for k in order}
o=sorted(range(len(l)),key=lambda i:-last[l[i]])
axes[2].legend([h[i] for i in o],[l[i] for i in o],loc='center left',bbox_to_anchor=(1.02,.5))
fig.suptitle('What each senescence caller actually selects: sequencing depth, within cell type\n'
             'left/centre: the two-section base (C6 tables); right: all %d sections'%NSEC11,
             fontsize=10,y=1.05)
fig.tight_layout()
for _e in ('png','pdf'):
    fig.savefig(FIG+'fig_phase3_caller_depth.'+_e,bbox_inches='tight')
plt.close(fig)
_basis='A_SENDER_FINAL_strict C6 (%d genes)'%sum(
    1 for _l in open('/workspace/genesets/A_SENDER_FINAL_strict.txt') if _l.strip())
_d=pd.concat([B.assign(basis='2 sections (7250 sham, 7259 SBR)'),
              B11.assign(basis='11 sections')],ignore_index=True)
_d['tierA_basis']=_basis
_d.to_csv(FIG+'fig_phase3_caller_depth_data.csv',index=False)

# ---------------------------------------------------------------- Fig C: Tier C identifiability
I=pd.read_csv('/workspace/results/phase3/tierC_ligand_identifiability.csv')
_I0=I[(I.on_panel==True)&I.section.str.contains('sbr')]
I=_I0[_I0.d_med_to_ligand_um.gt(0)]
CTRL={'Il1a','Ccl2','Cxcl1','Cxcl2','Cxcl5'}
fig,ax=plt.subplots(figsize=(6.4,4.6))
x=np.log10(I.dens_lig_per_mm2); y=np.log10(I.d_med_to_ligand_um)
b=np.polyfit(x,y,1); r=np.corrcoef(x,y)[0,1]
xs=np.linspace(x.min()-.15,x.max()+.15,50)
ax.plot(10**xs,10**np.polyval(b,xs),'-',color=P.MUTED,lw=1.6,zorder=1,
        label='fit: slope %.2f, $r^2$=%.3f'%(b[0],r*r))
ax.plot(10**xs,10**(np.log10(np.sqrt(np.log(2)/np.pi))+3-0.5*xs),'--',color=P.INK2,lw=1.2,zorder=1,
        label='Poisson prediction, slope $-1/2$')
for lab,m,cl in [('Section 9 control ligands',I.ligand.isin(CTRL),P.SERIES[1]),
                 ('other Tier C ligands',~I.ligand.isin(CTRL),P.SERIES[0])]:
    ax.plot(I.dens_lig_per_mm2[m],I.d_med_to_ligand_um[m],'o',color=cl,ms=6,mec=P.SURFACE,mew=1.2,
            ls='none',label=lab,zorder=3)
g=I.groupby('ligand')[['dens_lig_per_mm2','d_med_to_ligand_um']].median().sort_values(
    'dens_lig_per_mm2')
for i,(lig,row) in enumerate(g.iterrows()):
    ax.annotate(lig,(row.dens_lig_per_mm2,row.d_med_to_ligand_um),textcoords='offset points',
                xytext=(7,7) if i%2 else (-30,-13),fontsize=7.5,
                color=P.SERIES[1] if lig in CTRL else P.INK2)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('density of ligand$^+$ cells (per mm$^2$)')
ax.set_ylabel('median distance to nearest ligand$^+$ cell (µm)')
# the title used to claim "6 SBR sections x 14 Tier C ligands" = 84 points for
# a panel that draws 81: three (section, ligand) cells have a zero median
# distance and cannot go on a log axis.  Counted at draw time now.
ax.set_title('The Tier C distance regressor is a detection-rate readout\n'
             '%d SBR sections x %d Tier C ligands; %d of %d points '
             '(%d have zero median distance)'
             %(I.section.nunique(),I.ligand.nunique(),len(I),len(_I0),
               len(_I0)-len(I)),fontsize=9.5)
ax.legend(loc='upper right',fontsize=7.5)
fig.tight_layout()
for _e in ('png','pdf'):
    fig.savefig(FIG+'fig_phase3_tierC_identifiability.'+_e,bbox_inches='tight')
plt.close(fig)
_t=I[['section','ligand','dens_lig_per_mm2','d_med_to_ligand_um']].copy()
_t['is_section9_control']=I.ligand.isin(CTRL)
# same homogeneous-Poisson nearest-neighbour median as the dashed line above
_t['poisson_prediction_um']=np.sqrt(np.log(2)/(np.pi*I.dens_lig_per_mm2*1e-6))
_t.to_csv(FIG+'fig_phase3_tierC_identifiability_data.csv',index=False)
print('wrote 3 figures (png+pdf) and 3 _data.csv to',FIG)
