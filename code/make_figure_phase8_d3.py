#!/usr/bin/env python3
"""Phase 8, task 8.6 -- D3: what fixes DeepScence's sign, and what it costs.

(a) every anchor's DEPTH-PARTIALLED correlation with the D1 score, per section.
    A senescence score must run OPPOSITE to proliferation and to Lmnb1 and WITH
    Cdkn1a; the panel shows where each anchor actually points, and which section
    is the odd one out.
(b) depth- and type-matched caller agreement with DeepScence under the published
    anchor, the re-anchored score, and the sign-invariant |score| call.

Palette: sasp_palette.apply_style, as every other make_figure*.py.  Writes .png
and .pdf plus a *_data.csv with every plotted number.
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/workspace/code')
import sasp_palette as P
P.apply_style(matplotlib)
RES = '/workspace/results/phase3/'; FIG = '/workspace/figures/'
ARM = {'sham': P.SERIES[0], 'SBR': P.SERIES[1]}

D = pd.read_csv(RES + 'deepscence_anchor_decisions.csv')
G = pd.read_csv(RES + 'caller_agreement_matched_d3_11sections.csv')
D['sec'] = [s.split('_')[0] for s in D.section]
D = D.sort_values('sec')

ANCH = [('prho_ds_prolif', 'proliferation set, 8 genes disjoint from every Tier A / Tier B file', 'ρ < 0'),
        ('prho_ds_Lmnb1', 'Lmnb1 — but it IS in B_downstream_arrest and B_secondary_senescence', 'ρ < 0'),
        ('prho_ds_caller_consensus', 'consensus of the other three callers', 'ρ > 0'),
        ('prho_ds_Cdkn1a', 'Cdkn1a — the PUBLISHED anchor', 'ρ > 0')]

fig = plt.figure(figsize=(12.2, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=.28)
plotted = []

# ---------------------------------------------------------------- (a)
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(D))
w = 0.2
for i, (col, lab, want) in enumerate(ANCH):
    ax.bar(x + (i - 1.5) * w, D[col], width=w, color=P.SERIES[i], zorder=3,
           label='%s   [expect %s]' % (lab, want))
    for s, v in zip(D.sec, D[col]):
        plotted.append(dict(panel='a', section=s, series=col, value=v))
ax.axhline(0, color=P.INK2, lw=1.1)
ax.set_xticks(x)
ax.set_xticklabels([f'{s}\n{a}' for s, a in zip(D.sec, D.arm)], fontsize=7.6)
ax.set_ylabel('Spearman ρ with the D1 DeepScence score,\nafter removing log transcript counts')
ax.set_title('(a)  After removing depth the two caller-free anchors agree with each other in '
             '11/11 sections\n        and with the published CDKN1A anchor in 10/11.  '
             'The exception is 7250.\n        ▼ = 7248 and 7435, where the published '
             'anchor itself goes negative once depth is removed.',
             loc='left', fontsize=9.5)
ax.legend(loc='lower left', fontsize=6.9, ncol=1, handlelength=1.4,
          labelspacing=.28, borderpad=.4)
i7250 = int(np.where(D.sec.to_numpy() == '7250')[0][0])
ax.axvspan(i7250 - .5, i7250 + .5, color=P.STATUS['critical'], alpha=.10, zorder=0)
ax.annotate('7250: the only section where the proliferation and\nLmnb1 anchors both invert '
            'the published sign', xy=(i7250, D.prho_ds_prolif.iloc[i7250]),
            xytext=(i7250 - 3.1, 0.125), fontsize=7.4, color=P.STATUS['critical'],
            linespacing=1.4,
            arrowprops=dict(arrowstyle='->', color=P.STATUS['critical'], lw=1.1))
for _s in ('7248', '7435'):
    _i = int(np.where(D.sec.to_numpy() == _s)[0][0])
    ax.plot([_i + 1.5 * w], [0.004], marker='v', color=P.SERIES[3], ms=7,
            mec=P.INK, mew=0.8, zorder=5, clip_on=False)

# ---------------------------------------------------------------- (b)
ax = fig.add_subplot(gs[0, 1])
ORD = ['published', 'prolif', 'abs']
NAME = {'published': 'published\n(CDKN1A anchor)', 'prolif': 're-anchored\n(proliferation)',
        'abs': 'sign-invariant\n|score|'}
PAIRS = [('tierA_score', 'vs Tier A'), ('cdkn1a_counts', 'vs Cdkn1a⁺ (circular)'),
         ('senepy_score', 'vs SenePy')]
ax.axhline(1.0, color=P.INK2, lw=1.1, ls='--', zorder=1)
for j, (pair, plab) in enumerate(PAIRS):
    ys, zs = [], []
    for v in ORD:
        g = G[(G.ds_variant == v) & (G.A == pair)]
        pooled = g.n_both.sum() / g.exp_both_stratified.sum()
        ys.append(pooled)
        zs.append((g.n_both.sum() - g.exp_both_stratified.sum())
                  / np.sqrt((g.sd_both_stratified ** 2).sum()))
        plotted.append(dict(panel='b', section=v, series=pair, value=pooled))
    ax.plot(range(len(ORD)), ys, '-o', color=P.SERIES[j], ms=8, mec=P.SURFACE,
            mew=1.4, zorder=3, label=plab)
    for k, (yv, zv) in enumerate(zip(ys, zs)):
        ax.annotate('%.2f' % yv, (k, yv), textcoords='offset points',
                    xytext=(0, [11, -16, -16][j]), ha='center', fontsize=7.2,
                    color=P.SERIES[j])
ax.set_xticks(range(len(ORD)))
ax.set_xticklabels([NAME[v] for v in ORD], fontsize=8)
ax.set_ylabel('pooled agreement / chance\n(depth- and type-matched, 11 sections)')
ax.set_title('(b)  Re-anchoring halves the circularity;\n        the magnitude call keeps the '
             'Tier A signal', loc='left', fontsize=9.5)
ax.legend(loc='center left', fontsize=7.6)

fig.suptitle('Phase 8 / task 8.6 (D3): re-anchoring DeepScence off CDKN1A',
             fontsize=11, y=1.0)
pd.DataFrame(plotted).to_csv(RES + 'figure_phase8_d3_data.csv', index=False)
for e in ('png', 'pdf'):
    fig.savefig(FIG + 'figure_phase8_d3.' + e, bbox_inches='tight')
plt.close(fig)
print('wrote figures/figure_phase8_d3.{png,pdf}')
