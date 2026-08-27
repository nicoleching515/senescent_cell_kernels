#!/usr/bin/env python3
"""Render the C7-D2 tables as markdown, straight from results/phase8_d2/*.csv.

Every number in reports/CS_PHASE8_D2_DENOISE.md is emitted by this script so none of them
is hand-transcribed.

Producer of results/phase8_d2/d2_tables.md, which had no committed producer before:
    python code/report_d2_tables.py --out results/phase8_d2/d2_tables.md

Section 10.11 of PREREG_PHASE8.md: the `raw` rows of `rho_signed_dz_vs_depth` (-0.4673,
-0.1583; quoted in the prereg as -0.47 and -0.16) may never be reported -- they are the
direction of numerical noise on a shift of 0.0002-0.001 z-units.  This script used to hand
the author those two values pre-formatted in table C.  They are now WITHHELD in the rendered
table, with the reason printed beside it.  `--emit-prohibited` restores them for a purely
internal diagnostic; the rendered file must never be written with that flag.
"""
import pandas as pd, numpy as np, json, glob, os

OUT = '/workspace/results/phase8_d2/'
EMIT_PROHIBITED = False
DEPTH = {  # median panel-wide transcript_counts, from cells.parquet (see §2)
}


def md(df, cols=None, index=False):
    """Minimal markdown table writer.

    pandas' own .to_markdown() needs `tabulate`, which is not in requirements.txt, and
    installing anything into the main pinned environment for a report table is exactly the
    kind of drift this project has already been bitten by twice.  So: 12 lines here.
    """
    d = df[cols] if cols else df
    if index:
        d = d.reset_index()
    hdr = [str(c) for c in d.columns]
    def fmt(v):
        if isinstance(v, float):
            if abs(v) >= 1e6:
                return '%g' % v
            # correlations sit at 0.999999 in the determinism control; rounding them to
            # "1" would hide exactly the digit the control exists to show
            if 0.99 < abs(v) < 1.0:
                return '%.8f' % v
            return ('%.4f' % v).rstrip('0').rstrip('.')
        return str(v)
    rows = [[fmt(v) for v in r] for r in d.itertuples(index=False)]
    w = [max(len(hdr[i]), *(len(r[i]) for r in rows)) if rows else len(hdr[i])
         for i in range(len(hdr))]
    out = ['| ' + ' | '.join(h.ljust(w[i]) for i, h in enumerate(hdr)) + ' |',
           '|' + '|'.join('-' * (w[i] + 2) for i in range(len(hdr))) + '|']
    out += ['| ' + ' | '.join(r[i].ljust(w[i]) for i in range(len(hdr))) + ' |' for r in rows]
    return '\n'.join(out)


def main():
    ag = pd.read_csv(OUT + 'd2_agreement.csv')
    dp = pd.read_csv(OUT + 'd2_depth.csv')
    de = pd.read_csv(OUT + 'd2_deciles.csv')
    st = pd.read_csv(OUT + 'd2_normalisation_strength.csv') if os.path.exists(
        OUT + 'd2_normalisation_strength.csv') else pd.DataFrame()
    ag = ag.merge(dp[['section', 'config', 'median_transcript_counts',
                      'rho_depth_committed', 'rho_depth_alt', 'delta_rho_depth',
                      'rho_signed_dz_vs_depth', 'rho_abs_dz_vs_depth']],
                  on=['section', 'config'])
    ag['sec'] = ag.section.str.split('_').str[0]
    ag = ag.sort_values(['config', 'median_transcript_counts'])

    print('### A. Score agreement, per section (committed `denoise=False` vs alternative)\n')
    print(md(ag, ['config', 'sec', 'median_transcript_counts', 'n_cells', 'pearson_r',
                  'spearman_rho', 'anchor_sign_flipped', 'sign_alignment_ambiguous',
                  'mean_abs_dz', 'p95_abs_dz']))
    print('\n### B. Sender-status change at the operative thresholds\n')
    print(md(ag, ['config', 'sec', 'global_top5_n_called', 'global_top5_jaccard',
                  'global_top5_n_changed', 'global_top5_pct_of_called_changed',
                  'global_top1_jaccard', 'within_type_top5_jaccard',
                  'matched_top5_jaccard']))
    print('\n### C. Depth dependence of the shift\n')
    agc = ag.copy()
    if not EMIT_PROHIBITED:
        # PREREG_PHASE8.md section 10.11
        mask = agc.config.astype(str).str.fullmatch('raw')
        agc['rho_signed_dz_vs_depth'] = agc.rho_signed_dz_vs_depth.astype(object)
        agc.loc[mask, 'rho_signed_dz_vs_depth'] = 'WITHHELD'
    print(md(agc, ['config', 'sec', 'median_transcript_counts', 'rho_depth_committed',
                   'rho_depth_alt', 'delta_rho_depth', 'rho_signed_dz_vs_depth',
                   'rho_abs_dz_vs_depth']))
    if not EMIT_PROHIBITED:
        print('\n`rho_signed_dz_vs_depth` is WITHHELD for the `raw` rows under '
              'PREREG_PHASE8.md section 10.11: on a shift of 0.0002-0.001 z-units it is the '
              'direction of numerical noise, not a depth effect, and must not be quoted. '
              '`rho_abs_dz_vs_depth` and the `delta_rho_depth` column are unaffected.')
    if len(st):
        st['sec'] = st.section.str.split('_').str[0]
        print('\n### D. How much depth variation each configuration actually removed\n')
        print(md(st.sort_values(['config', 'sec']),
                 ['config', 'sec', 'sd_log_depth_before', 'sd_log_depth_after',
                  'frac_log_depth_variance_removed', 'depth_p90_over_p10_before',
                  'depth_p90_over_p10_after']))

    print('\n### E. Per-configuration summary across sections\n')
    g = ag.groupby('config').agg(
        n_sections=('section', 'nunique'),
        cells=('n_cells', 'sum'),
        min_r=('pearson_r', 'min'), median_r=('pearson_r', 'median'),
        max_r=('pearson_r', 'max'),
        min_J5=('global_top5_jaccard', 'min'),
        median_J5=('global_top5_jaccard', 'median'),
        max_J5=('global_top5_jaccard', 'max'),
        median_drho_depth=('delta_rho_depth', 'median'),
        max_abs_drho_depth=('delta_rho_depth', lambda s: s.abs().max())).reset_index()
    print(md(g.round(4)))

    print('\n### F. Score shift by within-section depth decile (mean z-score change)\n')
    for cfg in de.config.unique():
        sub = de[de.config == cfg]
        piv = sub.pivot_table(index='depth_decile', columns='section', values='mean_dz')
        piv.columns = [c.split('_')[0] for c in piv.columns]
        print('\n**%s**\n' % cfg)
        print(md(piv.round(4), index=True))
    print('\n### G. Top-5% call rate by depth decile, committed vs alternative\n')
    for cfg in de.config.unique():
        sub = de[de.config == cfg]
        p = sub.groupby('depth_decile')[['pct_called_committed', 'pct_called_alt',
                                         'pct_flipped']].mean().round(3)
        print('\n**%s** (mean over %d section(s))\n' % (cfg, sub.section.nunique()))
        print(md(p.round(3), index=True))

    cp_path = OUT + 'd2_celltype_composition.csv'
    if os.path.exists(cp_path):
        cp = pd.read_csv(cp_path)
        for cfg in cp.config.unique():
            sub = cp[cp.config == cfg]
            g2 = sub.groupby('cell_type').agg(
                bg_pct=('bg_pct', 'mean'),
                enrich_committed=('enrich_committed', 'mean'),
                enrich_alt=('enrich_alt', 'mean')).reset_index()
            g2['delta'] = (g2.enrich_alt - g2.enrich_committed).round(3)
            g2 = g2.sort_values('bg_pct', ascending=False)
            print('\n### G2. Cell-type enrichment of the top-5%% calls, %s '
                  '(mean over %d section(s))\n' % (cfg, sub.section.nunique()))
            print(md(g2.round(3)))

    print('\n### H. Run provenance\n')
    rows = []
    metas = sorted(glob.glob(OUT + 'runmeta_*.json'))
    # This section's row count is a function of what is on disk, not of any fixed set, so it
    # is stated rather than left to be inferred: a rebuilt tree with more (or fewer) runs
    # renders a different table and nothing would otherwise say so.
    print('_%d `runmeta_*.json` files present in %s at render time._\n' % (len(metas), OUT))
    for f in metas:
        d = json.load(open(f))
        dl = d.get('direction_log', {})
        rows.append(dict(config=d['config'], section=d['section'].split('_')[0],
                         n_cells=d['n_cells'], n_genes=d['n_genes'],
                         denoise=d['denoise'], anchor=d['anchor'],
                         node=dl.get('node'), reverse=dl.get('reverse'),
                         minutes=d.get('minutes')))
    print(md(pd.DataFrame(rows).sort_values(['config', 'section'])))


if __name__ == '__main__':
    import argparse, contextlib, sys
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', help='write the markdown here instead of stdout')
    ap.add_argument('--emit-prohibited', action='store_true',
                    help='internal diagnostic only: also print the section-10.11 `raw` '
                         'rho_signed_dz_vs_depth values.  Never use for a written file.')
    a = ap.parse_args()
    EMIT_PROHIBITED = a.emit_prohibited
    if a.out:
        if EMIT_PROHIBITED:
            sys.exit('refusing to write a file with --emit-prohibited (section 10.11)')
        with open(a.out, 'w') as fh, contextlib.redirect_stdout(fh):
            main()
        print('wrote', a.out)
    else:
        main()
