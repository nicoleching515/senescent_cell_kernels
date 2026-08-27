#!/usr/bin/env python3
"""Phase 8, task 8.4 -- the gate: did the caller-agreement headline move when
DeepScence coverage went from 2 of 11 M1 sections to 11 of 11?

Published headline (BIO_PHASE3.md §4.4, two-section base): the depth- and
type-matched top-5% calls of four caller pairs overlap at **0.93-1.22x chance**.
DeepScence vs Cdkn1a+ is held out of that band because DeepScence anchors the sign
of its bottleneck on CDKN1A, so the pair is CIRCULAR and is never pooled with the
others.  That exclusion is preserved here.

Reads   results/phase3/caller_agreement_matched_significance_{11sections,verify2sec}.csv
Writes  results/phase3/caller_coverage_gate.csv          (per pair, both bases)
        results/phase3/caller_coverage_gate_headline.csv (the single headline row)
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import norm, binomtest
RES = '/workspace/results/phase3/'
CIRCULAR = ('deepscence_score', 'cdkn1a_counts')     # BIO_PHASE3 §4.4: never pooled
# the four pairs whose 2-section values define the published 0.93-1.22x band
HEADLINE_PAIRS = [('tierA_score', 'senepy_score'), ('tierA_score', 'deepscence_score'),
                  ('tierA_score', 'cdkn1a_counts'), ('senepy_score', 'cdkn1a_counts')]
# the six Test-3-admissible ("in-band") sections Phase 3 fits on.  7239 is excluded
# there because 45% of its hepatocytes are Cdkn1a+, above §8's 20% ceiling, and 7250
# / 7361 / 7448 / 7450 are outside the Phase 3 fitting scope for other reasons.  The
# band is reported both ways so the widening cannot be blamed on a section Phase 3
# would not have fitted anyway.
INBAND = ['7001_liver_sham_Male_52-U1', '7248_liver_sham_Male_26-U1',
          '7259_liver_sbr_Male_26-U1', '7260_liver_sbr_Male_26-U1',
          '7352_liver_sham_Male_2-U1', '7435_liver_sham_Male_10-U1']


def short(a, b):
    f = lambda s: s.replace('_score', '').replace('_counts', '')
    return f(a) + ' vs ' + f(b)


def pooled(g):
    """Mantel-Haenszel pooling of the per-section stratified overlap tables."""
    obs = g.n_both.sum(); exp = g.exp_both_stratified.sum()
    sd = np.sqrt((g.sd_both_stratified ** 2).sum())
    z = (obs - exp) / sd
    return obs, exp, obs / exp, z, 2 * norm.sf(abs(z))


# ---------------------------------------------------------------------------
# Phase 8 / 8.7.  A gate row is only interpretable if BOTH of its bases were
# scored on the SAME sender definition.  The C6 promotion changed Tier A from
# 25 to 33 genes, so there are now two coverage comparisons, not one, and they
# answer different questions:
#
#   pre-C6  2-section vs 11-section  -> "did COVERAGE move the headline?" (8.4)
#   post-C6 2-section vs 11-section  -> the same question under the FROZEN sets
#
# Mixing a pre-C6 2-section base with a post-C6 11-section base measures the
# gene-set change and the coverage change at once.  Every row therefore carries
# `tierA_definition` and `tierA_n_genes`.
PRE = '/workspace/results/phase3_pre_c6/'
BASES = [
    # label, file, tierA definition, n genes, is the published record?
    ('2-section, pre-C6 Tier A (PUBLISHED)',
     PRE + 'caller_agreement_matched_significance_verify2sec.csv',
     'A_SENDER_FINAL_strict pre-C6', 25, True),
    ('11-section, pre-C6 Tier A (task 8.4)',
     PRE + 'caller_agreement_matched_significance_11sections.csv',
     'A_SENDER_FINAL_strict pre-C6', 25, False),
    ('2-section, post-C6 Tier A (FROZEN)',
     RES + 'caller_agreement_matched_significance_2sec_c6.csv',
     'A_SENDER_FINAL_strict C6', 33, False),
    ('11-section, post-C6 Tier A (FROZEN)',
     RES + 'caller_agreement_matched_significance_11sections.csv',
     'A_SENDER_FINAL_strict C6', 33, False),
]


def _load():
    out = []
    for lab, f, tdef, tn, pub in BASES:
        if not os.path.exists(f):
            print('MISSING, skipped:', f)
            continue
        out.append((lab, pd.read_csv(f), tdef, tn, pub))
    return out


def main():
    bases = _load()
    new = [d for l, d, _t, _n, _p in bases if l.startswith('11-section, post')][0]
    rows = []
    for A, B in [tuple(x) for x in new[['A', 'B']].drop_duplicates().to_numpy()]:
        for base, d, tdef, tn, _pub in bases:
            g = d[(d.A == A) & (d.B == B)]
            if g.empty:
                continue
            o, e, r, z, p = pooled(g)
            rows.append(dict(basis=base, tierA_definition=tdef,
                             tierA_n_genes=tn,
                             pair=short(A, B), A=A, B=B,
                             circular=(A, B) == CIRCULAR,
                             n_sections=len(g),
                             ratio_min=round(g.ratio_stratified.min(), 3),
                             ratio_median=round(g.ratio_stratified.median(), 3),
                             ratio_max=round(g.ratio_stratified.max(), 3),
                             n_sections_above_chance=int((g.ratio_stratified > 1).sum()),
                             n_sections_sig_above=int((g.z > 1.96).sum()),
                             n_sections_sig_below=int((g.z < -1.96).sum()),
                             obs_overlap=int(o), exp_overlap=round(e, 1),
                             pooled_ratio=round(r, 3), pooled_z=round(z, 2),
                             pooled_p='%.3g' % p))
    tab = pd.DataFrame(rows).sort_values(['circular', 'pair', 'tierA_n_genes',
                                          'n_sections'])
    tab.to_csv(RES + 'caller_coverage_gate.csv', index=False)

    # ---- the single headline: the band over the four non-circular headline pairs ----
    head = []
    heads = [(l, d, t, n) for l, d, t, n, _p in bases]
    for l, d, t, n, _p in bases:
        if l.startswith('11-section'):
            heads.append((l.replace('11-section', '6-section, in-band only'),
                          d[d.section.isin(INBAND)], t, n))
    for base, d, tdef, tn in heads:
        m = pd.Series(False, index=d.index)
        for A, B in HEADLINE_PAIRS:
            m |= (d.A == A) & (d.B == B)
        g = d[m]
        o, e, r, z, p = pooled(g)
        nsec = g.section.nunique()
        head.append(dict(basis=base, tierA_definition=tdef, tierA_n_genes=tn,
                         n_sections=nsec, n_pairs=len(HEADLINE_PAIRS),
                         n_values=len(g),
                         band_low=round(g.ratio_stratified.min(), 3),
                         band_high=round(g.ratio_stratified.max(), 3),
                         median=round(g.ratio_stratified.median(), 3),
                         pooled_ratio=round(r, 3), pooled_z=round(z, 2),
                         pooled_p='%.3g' % p,
                         n_above_chance=int((g.ratio_stratified > 1).sum()),
                         sign_test_p='%.3g' % binomtest(int((g.ratio_stratified > 1).sum()),
                                                        len(g), 0.5).pvalue))
    H = pd.DataFrame(head)
    H.to_csv(RES + 'caller_coverage_gate_headline.csv', index=False)

    pd.set_option('display.width', 220)
    print('=== HEADLINE: the published band, recomputed at full coverage ===')
    print(H.to_string(index=False))
    print('\n=== per pair, 2-section vs 11-section ===')
    print(tab.to_string(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
