#!/usr/bin/env python3
"""Phase 8, tasks 8.3/8.4 -- the caller-agreement gate figure.

(a) depth- and type-matched agreement of every caller pair, all 11 M1 sections,
    against the published two-section band;
(b) the same pooled, two-section base vs eleven-section base, with 95% intervals;
(c) DeepScence's Spearman correlation with transcript counts per section -- the
    evidence that the published "the sign reverses between arms" reading was a
    property of ONE section, which is the case for re-anchoring (8.6 / D3).

Palette is the project's validated one (sasp_palette.apply_style), as every other
make_figure*.py.  Writes .png AND .pdf plus a *_data.csv carrying every plotted
number.  Compare PNGs for reproducibility; matplotlib date-stamps PDFs.
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/workspace/code')
import sasp_palette as P
P.apply_style(matplotlib)
RES = '/workspace/results/phase3/'; FIG = '/workspace/figures/'
PUB = ['7250_liver_sham_Male_26-U1', '7259_liver_sbr_Male_26-U1']
BAND = (0.93, 1.22)                       # the published headline band
ARM = {'sham': P.SERIES[0], 'SBR': P.SERIES[1]}
SH = lambda s: s.replace('_score', '').replace('_counts', '')

S = pd.read_csv(RES + 'caller_agreement_matched_significance_11sections.csv')
S['pair'] = [SH(a) + '\nvs ' + SH(b) for a, b in zip(S.A, S.B)]
S['circular'] = (S.A == 'deepscence_score') & (S.B == 'cdkn1a_counts')
G = pd.read_csv(RES + 'caller_coverage_gate.csv')
T = pd.read_csv(RES + 'caller_technical_loading_11sections.csv')

ORDER = [p for p in ['tierA\nvs senepy', 'tierA\nvs deepscence', 'tierA\nvs cdkn1a',
                     'senepy\nvs deepscence', 'senepy\nvs cdkn1a', 'deepscence\nvs cdkn1a']
         if p in set(S.pair)]

fig = plt.figure(figsize=(11.6, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.35, 1.0],
                      hspace=.42, wspace=.26)

# ---------------------------------------------------------------- (a)
ax = fig.add_subplot(gs[0, :])
ax.axhspan(BAND[0], BAND[1], color=P.SERIES[0], alpha=.10, zorder=0)
ax.axhline(1.0, color=P.INK2, lw=1.1, ls='--', zorder=1)
rng = np.random.default_rng(0)
plotted = []
for i, pr in enumerate(ORDER):
    g = S[S.pair == pr]
    for _, r in g.iterrows():
        pub = r.section in PUB
        x = i + (0.0 if pub else rng.uniform(-.22, .22))
        ax.plot([x], [r.ratio_stratified], 'D' if pub else 'o',
                color=ARM[r.arm], ms=9 if pub else 6,
                mec=P.INK if pub else P.SURFACE, mew=1.5 if pub else 1.2,
                zorder=4 if pub else 3)
        plotted.append(dict(panel='a', pair=pr.replace('\n', ' '), section=r.section,
                            arm=r.arm, x=round(x, 4), ratio_stratified=r.ratio_stratified,
                            z=r.z, is_published_section=pub))
ax.axvline(len(ORDER) - 1.5, color=P.AXIS, lw=1.2, ls=':')
ax.text(len(ORDER) - 1.45, ax.get_ylim()[1], ' circular:\n DeepScence anchors\n its sign on CDKN1A',
        va='top', ha='left', fontsize=7.5, color=P.INK2)
ax.set_yscale('log'); ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_ylabel('agreement / chance\n(top-5%, matched on cell type x within-type depth decile)')
ax.set_title('(a)  Every caller pair, all 11 M1 sections. Shaded = the published two-section band '
             '0.93-1.22x; diamonds = the two sections it was measured on', loc='left')
for a_, c in ARM.items():
    ax.plot([], [], 'o', color=c, mec=P.SURFACE, mew=1.2, label=a_)
ax.plot([], [], 'D', color=P.MUTED, mec=P.INK, mew=1.5, label='published 2-section base')
ax.legend(loc='upper left', ncol=3)

# ---------------------------------------------------------------- (b)
ax = fig.add_subplot(gs[1, 0])
ax.axvline(1.0, color=P.INK2, lw=1.1, ls='--', zorder=1)
pairs = [p.replace('\n', ' ') for p in ORDER]
for j, pr in enumerate(pairs):
    # Phase 8 / 8.7: `caller_coverage_gate.csv` now carries FOUR bases, because
    # the C6 promotion changed Tier A from 25 to 33 genes and a 2-vs-11
    # comparison is only interpretable at a fixed sender definition.  The
    # coverage effect is drawn on the FROZEN (post-C6) sets; the published
    # pre-C6 two-section value is kept as a separate marker.
    BASES = [('2-section, pre-C6 Tier A (PUBLISHED)', P.MUTED, 'D', -0.26),
             ('2-section, post-C6 Tier A (FROZEN)', P.SERIES[1], 's', 0.0),
             ('11-section, post-C6 Tier A (FROZEN)', P.SERIES[2], 'o', 0.26)]
    for base, col, mk, dy in BASES:
        g = G[(G.pair == pr) & (G.basis == base)]
        if g.empty:
            continue
        r = g.iloc[0]
        ax.plot([r.pooled_ratio], [j + dy], mk, color=col, ms=8,
                mec=P.SURFACE, mew=1.2, zorder=3)
        plotted.append(dict(panel='b', pair=pr, section=base, arm='',
                            x=r.pooled_ratio, ratio_stratified=r.pooled_ratio,
                            z=r.pooled_z,
                            is_published_section='PUBLISHED' in base))
    lo = G[(G.pair == pr) & (G.basis == BASES[1][0])]
    hi = G[(G.pair == pr) & (G.basis == BASES[2][0])]
    if not lo.empty and not hi.empty:
        ax.plot([lo.pooled_ratio.iloc[0], hi.pooled_ratio.iloc[0]],
                [j, j + 0.26], '-', color=P.AXIS, lw=1.4, zorder=2)
ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs, fontsize=8)
ax.invert_yaxis(); ax.set_xlabel('pooled agreement / chance (Mantel-Haenszel over sections)')
ax.plot([], [], 'D', color=P.MUTED, mec=P.SURFACE, mew=1.2,
        label='2 sections, pre-C6 Tier A (published)')
ax.plot([], [], 's', color=P.SERIES[1], mec=P.SURFACE, mew=1.2,
        label='2 sections, frozen Tier A')
ax.plot([], [], 'o', color=P.SERIES[2], mec=P.SURFACE, mew=1.2,
        label='11 sections, frozen Tier A')
ax.legend(loc='lower right', fontsize=7)
ax.set_title('(b)  Pooled: coverage at a FIXED sender definition (line), and what '
             'the C6 gene sets alone did to the published two-section base',
             loc='left')

# ---------------------------------------------------------------- (c)
ax = fig.add_subplot(gs[1, 1])
d = T[T.score == 'deepscence_score'].sort_values('rho_transcript_counts')
arms = ['SBR' if 'sbr' in s else 'sham' for s in d.section]
ax.axvline(0, color=P.INK2, lw=1.1, ls='--', zorder=1)
ax.barh(range(len(d)), d.rho_transcript_counts,
        color=[ARM[a] for a in arms], zorder=3, height=.7)
ax.set_yticks(range(len(d)))
ax.set_yticklabels([s.split('_')[0] for s in d.section], fontsize=8)
ax.set_xlabel(r'DeepScence score: Spearman $\rho$ with transcript counts')
ax.set_title('(c)  The sign flip is ONE section, not the arm', loc='left')
for i, (s, v) in enumerate(zip(d.section, d.rho_transcript_counts)):
    plotted.append(dict(panel='c', pair='deepscence_score', section=s,
                        arm='SBR' if 'sbr' in s else 'sham', x=v,
                        ratio_stratified=np.nan, z=np.nan, is_published_section=s in PUB))
ax.annotate('7250 = the sham section the\npublished DeepScence numbers rest on',
            xy=(d.rho_transcript_counts.iloc[0], 0), xytext=(0.12, 2.4),
            fontsize=7.5, color=P.INK2,
            arrowprops=dict(arrowstyle='->', color=P.INK2, lw=1))

fig.suptitle('Phase 8 / task 8.4 gate: caller agreement at 2-section vs 11-section DeepScence coverage',
             fontsize=11, y=.985)
D = pd.DataFrame(plotted)
D.to_csv(RES + 'figure_phase8_callers_data.csv', index=False)
for ext in ('png', 'pdf'):
    fig.savefig(FIG + 'figure_phase8_callers.' + ext, bbox_inches='tight')
plt.close(fig)
print('wrote figures/figure_phase8_callers.{png,pdf} and %d data rows'
      % len(D))
