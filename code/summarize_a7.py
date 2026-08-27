#!/usr/bin/env python3
"""Phase 8, task 8.5 -- verdict for audit test A7 on the mouse arm.

A7 (addendum §13): the kernel fitted against negative-control-probe counts vs
distance-to-nearest-sender MUST BE FLAT.  §13 states the test but not a numeric
pass rule, so this file states one and reports against it three ways:

  (i)   AMPLITUDE.  |beta| / sd_y -- the fitted amplitude in units of the
        response's own SD inside the fitted cell set -- against the same
        quantity for the seven biological Tier B modules from main_fits.csv,
        the same sections, sender calls, cell types and estimator.
  (ii)  DIRECTION.  A real spatial technical gradient has a CONSISTENT SIGN
        across sections and cell types; noise does not.  Reported as the
        fraction of fits with beta > 0 and a section-clustered t-test on the
        mean amplitude, clustering on section (11 clusters) because fits from
        one section are not independent.
  (iii) POWER.  The median half-width of the 400-replicate spatial block
        bootstrap CI on beta/sd_y, i.e. the smallest gradient a single A7 fit
        could have resolved.  Without this the word "flat" is unfalsifiable.

Reads  results/phase3/a7_control_probe_{fits,provenance}.csv, main_fits.csv
Writes results/phase3/a7_summary.csv, a7_verdict.txt
"""
import sys, numpy as np, pandas as pd
from scipy.stats import t as tdist
RES = '/workspace/results/phase3/'
DESIGNS = [('base', 'naive (intercept only)'),
           ('n6', '+N6 neighbour baseline'),
           ('n5', '+N5 technical covariates'),
           ('n6n5', '+N6+N5 (full nuisance design)'),
           ('n2', 'N2 matched-decoy contrast')]
CALLS = ['tierA_p95', 'cdkn1a_pos']


def add(df):
    for k, _ in DESIGNS:
        df['e_' + k] = df['beta_' + k] / df.sd_y
        df['sig_' + k] = (df['beta_%s_lo' % k] * df['beta_%s_hi' % k]) > 0
        df['hw_' + k] = (df['beta_%s_hi' % k] - df['beta_%s_lo' % k]) / 2 / df.sd_y
    return df


def clustered(df, col, cluster='section'):
    """Mean of `col` with a section-clustered 95% CI (t on the cluster means)."""
    m = df.groupby(cluster)[col].mean().to_numpy()
    m = m[np.isfinite(m)]
    if m.size < 3:
        return np.nan, np.nan, np.nan, np.nan
    mu = m.mean(); se = m.std(ddof=1) / np.sqrt(m.size)
    tc = tdist.ppf(.975, m.size - 1)
    return mu, mu - tc * se, mu + tc * se, 2 * tdist.sf(abs(mu / se), m.size - 1)


def main():
    A = add(pd.read_csv(RES + 'a7_control_probe_fits.csv'))
    A = A[A.get('skip').isna()] if 'skip' in A.columns else A
    M = pd.read_csv(RES + 'main_fits.csv')
    M = add(M[(M.stratum == 'all') & M.call.isin(CALLS)].copy())
    M['response'] = 'BIOLOGICAL MODULES (reference)'
    prov = pd.read_csv(RES + 'a7_control_probe_provenance.csv')

    rows = []
    for src in (A, M):
        for resp, g in src.groupby('response'):
            for k, dl in DESIGNS:
                mu, lo, hi, p = clustered(g, 'e_' + k)
                rows.append(dict(response=resp, design=k, design_label=dl,
                                 n_fits=len(g), n_sections=g.section.nunique(),
                                 median_abs_amplitude=round(g['e_' + k].abs().median(), 4),
                                 median_signed_amplitude=round(g['e_' + k].median(), 4),
                                 frac_positive=round(float((g['beta_' + k] > 0).mean()), 3),
                                 frac_CI_excludes_zero=round(float(g['sig_' + k].mean()), 3),
                                 clustered_mean=round(mu, 4),
                                 clustered_lo=round(lo, 4), clustered_hi=round(hi, 4),
                                 clustered_p=('%.3g' % p),
                                 median_CI_halfwidth=round(g['hw_' + k].median(), 4)))
    S = pd.DataFrame(rows)
    S.to_csv(RES + 'a7_summary.csv', index=False)

    ctrl = S[S.response != 'BIOLOGICAL MODULES (reference)']
    bio = S[S.response == 'BIOLOGICAL MODULES (reference)'].set_index('design')
    L = []
    L.append('A7 -- negative-control-probe kernel, MOUSE arm (M1), 11 sections x '
             '%d sender calls' % len(CALLS))
    L.append('%d control fits, %d biological-module fits, same estimator '
             '(WINDOW=100um, 40-point lambda grid, MIN_RECEIVERS=2000, '
             '400-replicate spatial block bootstrap).' % (len(A), len(M)))
    L.append('')
    L.append('Control-probe sparsity (why power is the binding constraint):')
    for r, g in prov.groupby('response'):
        L.append('  %-22s mean %.4f counts/cell, %.2f%% of cells non-zero'
                 % (r, g.mean_per_cell.mean(), 100 * g.frac_cells_nonzero.mean()))
    L.append('')
    for k, dl in DESIGNS:
        c = ctrl[ctrl.design == k]; b = bio.loc[k]
        L.append('%-30s control |beta|/sd med %.4f (sign +: %.2f-%.2f, CI excl 0: '
                 '%.2f-%.2f) | modules %.4f (sign +: %.2f, CI excl 0: %.2f)'
                 % (dl, c.median_abs_amplitude.median(),
                    c.frac_positive.min(), c.frac_positive.max(),
                    c.frac_CI_excludes_zero.min(), c.frac_CI_excludes_zero.max(),
                    b.median_abs_amplitude, b.frac_positive, b.frac_CI_excludes_zero))
    L.append('')
    L.append('Smallest amplitude one A7 fit could resolve (median CI half-width): '
             '%.4f SD naive, %.4f SD conditioned.' %
             (ctrl[ctrl.design == 'base'].median_CI_halfwidth.median(),
              ctrl[ctrl.design == 'n6n5'].median_CI_halfwidth.median()))
    L.append('Biological module amplitude to be ruled out: %.4f SD naive, '
             '%.4f SD conditioned.' % (bio.loc['base'].median_abs_amplitude,
                                       bio.loc['n6n5'].median_abs_amplitude))
    txt = '\n'.join(L)
    open(RES + 'a7_verdict.txt', 'w').write(txt + '\n')
    pd.set_option('display.width', 250)
    print(txt)
    print()
    print(S.to_string(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
