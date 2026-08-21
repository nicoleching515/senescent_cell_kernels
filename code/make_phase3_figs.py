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
fig.tight_layout(); fig.savefig(FIG+'fig_phase3_composition.png',bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Fig B: caller depth bias
B=pd.read_csv('/workspace/results/phase3/caller_within_type_depth_bias.csv')
NAME={'tierA_score':'Tier A score','senepy_score':'SenePy','deepscence_score':'DeepScence',
      'cdkn1a_counts':'Cdkn1a > 0'}
order=['tierA_score','deepscence_score','senepy_score','cdkn1a_counts']
col=dict(zip(order,P.SERIES[:4]))
fig,axes=plt.subplots(1,2,figsize=(8.6,3.6),sharey=True)
for ax,(sec,lab) in zip(axes,[('sham','sham (7250, 26 wk)'),('sbr','SBR (7259, 26 wk)')]):
    s=B[B.section==sec]
    ax.axhline(1,color=P.MUTED,lw=1,ls='--',zorder=1)
    for k in order:
        q=s[s.caller==k].sort_values('within_type_depth_quintile')
        if not len(q): continue
        ax.plot(range(1,len(q)+1),q.enrichment,'-o',color=col[k],ms=6,mec=P.SURFACE,mew=1.4,
                label=NAME[k],zorder=3)
    ax.set_yscale('log'); ax.set_title(lab)
    ax.set_xticks([1,2,3,4,5]); ax.set_xticklabels(['Q1\nlow','Q2','Q3','Q4','Q5\nhigh'])
    ax.set_xlabel('transcript-count quintile, within cell type')
axes[0].set_ylabel('enrichment of top-5% calls\n(1 = no depth bias)')
h,l=axes[1].get_legend_handles_labels()          # order the legend by the value at Q5 so it
last={NAME[k]:B[(B.section=='sbr')&(B.caller==k)].sort_values(                # reads as a
      'within_type_depth_quintile').enrichment.iloc[-1] for k in order}       # direct label
o=sorted(range(len(l)),key=lambda i:-last[l[i]])
axes[1].legend([h[i] for i in o],[l[i] for i in o],loc='center left',bbox_to_anchor=(1.02,.5))
fig.suptitle('What each senescence caller actually selects: sequencing depth, within cell type',
             fontsize=10,y=1.02)
fig.tight_layout(); fig.savefig(FIG+'fig_phase3_caller_depth.png',bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Fig C: Tier C identifiability
I=pd.read_csv('/workspace/results/phase3/tierC_ligand_identifiability.csv')
I=I[(I.on_panel==True)&I.section.str.contains('sbr')&I.d_med_to_ligand_um.gt(0)]
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
ax.set_title('The Tier C distance regressor is a detection-rate readout\n'
             '6 SBR sections x 14 Tier C ligands',fontsize=9.5)
ax.legend(loc='upper right',fontsize=7.5)
fig.tight_layout(); fig.savefig(FIG+'fig_phase3_tierC_identifiability.png',bbox_inches='tight')
print('wrote 3 figures to',FIG)
