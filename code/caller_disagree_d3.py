#!/usr/bin/env python3
"""Phase 8, task 8.6 -- what re-anchoring and a sign-invariant call do to the
caller-agreement statistic.

Reuses caller_disagree_all's matched-strata machinery verbatim (same top-5%
inside each cell type x within-type depth decile, same marginal `ratio`, same
exact stratified null) and swaps only the DeepScence column:

  published    the D1 score as shipped, CDKN1A anchor           (the 8.3/8.4 number)
  prolif       re-anchored on the depth-partialled proliferation set (D3 primary)
  lmnb1        re-anchored on depth-partialled Lmnb1                 (D3 secondary)
  consensus    re-anchored on the other three callers                (D3, circular)
  abs          |score| -- the SIGN-INVARIANT call.  A top-5% by magnitude cannot
               be moved by any anchor, so whatever agreement survives here is the
               part of the DeepScence result that does not depend on the sign at all.

Writes results/phase3/caller_agreement_matched_d3_11sections.csv
"""
import sys, numpy as np, pandas as pd, warnings
sys.path.insert(0, '/workspace/code')
import caller_disagree_all as CA
warnings.filterwarnings('ignore')

PROC = '/workspace/data/processed/'
RES = '/workspace/results/phase3/'
VARIANTS = {'published': 'ds_published', 'prolif': 'ds_prolif_anchor',
            'lmnb1': 'ds_lmnb1_anchor', 'consensus': 'ds_consensus_anchor',
            'abs': 'ds_abs'}
OTHERS = ['tierA_score', 'senepy_score', 'cdkn1a_counts']


def main():
    rows = []
    for sec, arm, _ in CA.SECTIONS:
        df = CA._load(sec)
        d3 = pd.read_csv(PROC + 'deepscence_d3_%s.csv' % sec).set_index('cell_id')
        df = df.join(d3)
        print('...', sec, len(df), flush=True)
        for vname, col in VARIANTS.items():
            for other in OTHERS:
                ok = np.isfinite(df[col].to_numpy(float)) & np.isfinite(df[other].to_numpy(float))
                if ok.sum() < 200:
                    continue
                sub = df[ok]
                fa, sa = CA._matched_flags(sub, other, True)
                fb, _ = CA._matched_flags(sub, col, True)
                uni = int((fa | fb).sum()); inter = int((fa & fb).sum())
                if uni == 0:
                    continue
                jac = inter / uni
                ch = fa.mean() * fb.mean() * len(sub) / uni
                e, sd, nst = CA._stratified_null(fa, fb, sa)
                rows.append(dict(section=sec, arm=arm, ds_variant=vname,
                                 A=other, B='deepscence_' + vname,
                                 circular=(other == 'cdkn1a_counts' and vname == 'published')
                                 or (vname == 'consensus'),
                                 n=int(len(sub)), n_strata=nst,
                                 n_A=int(fa.sum()), n_B=int(fb.sum()), n_both=inter,
                                 jaccard=round(float(jac), 5), chance=round(float(ch), 5),
                                 ratio_marginal=round(float(jac / ch), 3) if ch else None,
                                 exp_both_stratified=round(e, 2),
                                 sd_both_stratified=round(sd, 3),
                                 ratio_stratified=round(inter / e, 3) if e else None,
                                 z=round((inter - e) / sd, 2) if sd else None))
    D = pd.DataFrame(rows)
    D.to_csv(RES + 'caller_agreement_matched_d3_11sections.csv', index=False)
    pd.set_option('display.width', 240)
    print()
    print(D.pivot_table(index=['ds_variant', 'A'], values='ratio_stratified',
                        aggfunc=['min', 'median', 'max']).round(3).to_string())
    print()
    p = D.groupby(['ds_variant', 'A']).apply(
        lambda g: pd.Series(dict(pooled=g.n_both.sum() / g.exp_both_stratified.sum(),
                                 z=(g.n_both.sum() - g.exp_both_stratified.sum())
                                 / np.sqrt((g.sd_both_stratified ** 2).sum()),
                                 n_above=(g.ratio_stratified > 1).sum())),
        include_groups=False)
    print(p.round(3).to_string())


if __name__ == '__main__':
    main()
