#!/usr/bin/env python3
"""Phase 8, task 8.3 -- the five caller-agreement tables at ELEVEN-section coverage.

Why this file exists rather than an edit to caller_disagree.py:
  caller_disagree.py writes fixed output paths (`caller_*.csv`).  Phase 7 addendum
  §5 requires the TWO-section values to be reported ALONGSIDE the eleven-section
  values, so the committed two-section tables must survive untouched.  This script
  REUSES caller_disagree.run() verbatim for the four tables that script already
  produces, adds the two depth/type-matched tables that §5 also lists, and writes
  everything under a `_11sections` suffix.

  The only edit made to caller_disagree.py is DS_ALIAS, which maps the two section
  directory names onto the preserved `deepscence_{sham,sbr}.csv` filenames.  Called
  with tag='sham' it behaves exactly as before; this is verified by --verify.
  (That edit was lost once and re-applied on 2026-08-27; it is now committed, as
  `DS_ALIAS={v:k for k,v in SAMP.items()}` in caller_disagree.py.)

What changed between the two runs is COVERAGE ONLY.  Every DeepScence score comes
from run_deepscence_all.py at settings identical to run_deepscence.py
(denoise=False, random_state=0, published CDKN1A anchor, MGI 1:1 ortholog remap,
>=20 counts/cell), so any movement in the headline is attributable to coverage.

Reconstruction note: the producer of `caller_agreement_depth_and_type_matched.csv`
and `caller_within_type_depth_bias.csv` was not committed to code/ (only its two
outputs were, dated 2026-08-20 18:02).  It is reconstructed here and validated by
`--verify`, which re-runs the two original sections through this file and requires
every cell of all six published tables to come back EXACTLY as committed.

Usage:
  caller_disagree_all.py --verify         # reproduce the 2-section tables, compare
  caller_disagree_all.py --set 2sec_c6    # the post-C6 two-section tables (*_2sec_c6.csv)
  caller_disagree_all.py --all            # write the 11-section tables
  add --out-dir DIR to write anywhere but results/phase3/ (used to verify without
  overwriting the committed tables).
"""
import sys, os, numpy as np, pandas as pd, warnings
sys.path.insert(0, '/workspace/code')
import caller_disagree as CD
warnings.filterwarnings('ignore')

OUT = '/workspace/results/phase3/'
SCORES = ['tierA_score', 'senepy_score', 'deepscence_score', 'cdkn1a_counts']
MIN_STRATUM = 50          # same floor caller_disagree.py uses for its within-type calls
N_DEPTH_DECILES = 10

# all eleven M1 sections, in the order BIO_PHASE3 §5 tabulates them
SECTIONS = [
    ('7001_liver_sham_Male_52-U1', 'sham', 52),
    ('7239_liver_sbr_Male_52-U1',  'SBR',  52),
    ('7248_liver_sham_Male_26-U1', 'sham', 26),
    ('7250_liver_sham_Male_26-U1', 'sham', 26),
    ('7259_liver_sbr_Male_26-U1',  'SBR',  26),
    ('7260_liver_sbr_Male_26-U1',  'SBR',  26),
    ('7352_liver_sham_Male_2-U1',  'sham', 2),
    ('7361_liver_sbr_Male_2-U1',   'SBR',  2),
    ('7435_liver_sham_Male_10-U1', 'sham', 10),
    ('7448_liver_sbr_Male_10-U1',  'SBR',  10),
    ('7450_liver_sbr_Male_10-U1',  'SBR',  10),
]


def _load(tag):
    """Exactly caller_disagree.run()'s input assembly, up to the EXCL filter."""
    sen = pd.read_csv(CD.PROC + 'senders_%s.csv' % tag).set_index('cell_id')
    ana = pd.read_csv(CD.PROC + 'anatomy_%s.csv' % tag).set_index('cell_id')
    cells = pd.read_parquet(CD.RAW + CD.SAMP.get(tag, tag) + '/cells.parquet').set_index('cell_id')
    ds = CD.PROC + 'deepscence_%s.csv' % CD.DS_ALIAS.get(tag, tag)
    if os.path.exists(ds):
        sen = sen.join(pd.read_csv(ds).set_index('cell_id'), rsuffix='_ds')
    df = sen.join(ana[['zonation_score', 'compartment_label']]).join(
        cells[['transcript_counts', 'cell_area', 'nucleus_area', 'segmentation_method']])
    return df[~df.cell_type.isin(CD.EXCL)].copy()


def _within_type_flags(sub, s):
    """Top-5% of `s` recomputed inside each cell type."""
    v = sub[s].to_numpy(float); ok = np.isfinite(v)
    f = np.zeros(len(sub), bool)
    ct = sub.cell_type.to_numpy()
    for c in pd.unique(ct):
        m = (ct == c) & ok
        if m.sum() < MIN_STRATUM:
            continue
        f[m] = v[m] > np.nanpercentile(v[m], 95)
    return f


def _matched_flags(sub, s, want_strata=False):
    """Top-5% of `s` recomputed inside each (cell type x within-type depth decile).

    With want_strata=True also returns an integer stratum id per cell (-1 = not in
    any admissible stratum), which is what the stratified-exact null below needs.
    """
    v = sub[s].to_numpy(float); ok = np.isfinite(v)
    tc = sub.transcript_counts.to_numpy(float)
    ct = sub.cell_type.to_numpy()
    f = np.zeros(len(sub), bool)
    sid = np.full(len(sub), -1, np.int32); nxt = 0
    for c in pd.unique(ct):
        m0 = (ct == c) & ok
        if m0.sum() < MIN_STRATUM:
            continue
        idx = np.where(m0)[0]
        try:
            d = pd.qcut(tc[idx], N_DEPTH_DECILES, labels=False, duplicates='drop')
        except Exception:
            d = np.zeros(len(idx), int)
        for k in np.unique(d):
            ii = idx[d == k]
            if len(ii) < MIN_STRATUM:
                continue
            f[ii] = v[ii] > np.nanpercentile(v[ii], 95)
            sid[ii] = nxt; nxt += 1
    return (f, sid) if want_strata else f


def _stratified_null(fa, fb, sid):
    """Exact conditional (Mantel-Haenszel) null for the overlap of two top-5% call
    sets that were each thresholded INSIDE the same strata.

    The published `ratio` column is observed_overlap / (n_A * n_B / n), i.e. an
    observed/expected overlap under a MARGINAL independence null.  That null
    ignores the fact that both calls are stratum-balanced by construction.  Under
    the correct null -- B's calls permuted at random WITHIN each stratum, A held
    fixed -- the overlap is a sum of independent hypergeometrics, one per stratum,
    whose mean and variance are closed-form.  Returns (expected, sd, n_strata).
    """
    m = sid >= 0
    if not m.any():
        return np.nan, np.nan, 0
    sid = sid[m]; a = fa[m]; b = fb[m]
    ns = np.bincount(sid).astype(float)
    na = np.bincount(sid, weights=a.astype(float))
    nb = np.bincount(sid, weights=b.astype(float))
    keep = ns >= 2
    ns, na, nb = ns[keep], na[keep], nb[keep]
    exp = float((na * nb / ns).sum())
    var = float((na * nb * (ns - na) * (ns - nb) / (ns * ns * (ns - 1))).sum())
    return exp, float(np.sqrt(var)), int(keep.sum())


def _within_type_depth_quintile(df):
    q = pd.Series(index=df.index, dtype=object)
    for c in df.cell_type.unique():
        m = (df.cell_type == c)
        if m.sum() < MIN_STRATUM:
            continue
        q[m] = pd.qcut(df.transcript_counts[m], 5,
                       labels=['Q%d' % i for i in range(1, 6)], duplicates='drop')
    return q


def extra_tables(tag, arm, label=None):
    """The two tables §5 lists that caller_disagree.py does not write."""
    label = label or tag
    df = _load(tag)
    sig = []

    # ---- within-type depth bias of the within-type top-5% calls ----
    df['wq'] = _within_type_depth_quintile(df)
    bias = []
    for s in SCORES:
        if s not in df or not df[s].notna().any():
            continue
        f = _within_type_flags(df, s)
        for qq in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
            m = (df.wq == qq).to_numpy()
            bg = m.mean()
            bias.append(dict(section=label, arm=arm, caller=s, within_type_depth_quintile=qq,
                             bg_pct=round(100 * bg, 2),
                             enrichment=round(float((f & m).sum() / max(f.sum(), 1) / bg), 3)))

    # ---- pairwise agreement, matched on cell type AND within-type depth decile ----
    pair = []
    for i in range(len(SCORES)):
        for j in range(i + 1, len(SCORES)):
            a_s, b_s = SCORES[i], SCORES[j]
            if a_s not in df or b_s not in df:
                continue
            ok = np.isfinite(df[a_s].to_numpy(float)) & np.isfinite(df[b_s].to_numpy(float))
            if ok.sum() < 200:
                continue
            sub = df[ok]
            fa, sa = _matched_flags(sub, a_s, True); fb, sb = _matched_flags(sub, b_s, True)
            uni = int((fa | fb).sum()); inter = int((fa & fb).sum())
            if uni == 0:
                continue
            jac = inter / uni
            ch = fa.mean() * fb.mean() * len(sub) / uni
            row = dict(section=label, arm=arm, A=a_s, B=b_s, n=int(len(sub)),
                       n_A=int(fa.sum()), n_B=int(fb.sum()), n_both=inter,
                       jaccard=round(float(jac), 5), chance=round(float(ch), 5),
                       ratio=round(float(jac / ch), 3) if ch else None)
            pair.append(row)
            # sa and sb are built from the same (cell type x depth decile) partition
            # of the same rows, so they are identical; assert rather than assume.
            assert (sa == sb).all(), 'stratum ids diverge for %s/%s' % (a_s, b_s)
            e, sd, nst = _stratified_null(fa, fb, sa)
            sig.append(dict(section=label, arm=arm, A=a_s, B=b_s, n=int(len(sub)),
                            n_strata=nst, n_A=int(fa.sum()), n_B=int(fb.sum()),
                            n_both=inter,
                            exp_both_marginal=round(float(fa.sum()) * fb.sum() / len(sub), 2),
                            ratio_marginal=round(float(jac / ch), 3) if ch else None,
                            exp_both_stratified=round(e, 2), sd_both_stratified=round(sd, 3),
                            ratio_stratified=round(inter / e, 3) if e else None,
                            z=round((inter - e) / sd, 2) if sd else None))
    return pd.DataFrame(bias), pd.DataFrame(pair), pd.DataFrame(sig)


def run_set(specs, suffix, out=None, require_deepscence=True):
    out = out or OUT
    if require_deepscence:
        # D6 guard: data/processed/deepscence_*.csv is gitignored, and every caller table
        # below silently drops the DeepScence columns when it is absent (exit 0, different
        # caller set, plausible-looking output).  Refuse instead.
        miss = [t for t, _, _ in specs
                if not os.path.exists(CD.PROC + 'deepscence_%s.csv' % CD.DS_ALIAS.get(t, t))]
        if miss:
            raise SystemExit('no DeepScence scores for: %s\n'
                             '  expected %sdeepscence_<tag>.csv (gitignored; rebuild with '
                             'run_deepscence_all.py)\n'
                             '  re-run with --allow-missing-deepscence to accept a caller set '
                             'WITHOUT DeepScence' % (miss, CD.PROC))
    T = []; C = []; S = []; P = []; B = []; M = []; G = []
    for tag, arm, label in specs:
        print('...', label, flush=True)
        t, c, s, p, _ = CD.run(tag, arm)
        for d in (t, c, s, p):
            d['section'] = label
        b, m, g = extra_tables(tag, arm, label)
        T.append(t); C.append(c); S.append(s); P.append(p); B.append(b); M.append(m); G.append(g)
    outs = {
        'caller_technical_loading':               pd.concat(T),
        'caller_celltype_composition':            pd.concat(C),
        'caller_strata':                          pd.concat(S),
        'caller_pairwise_agreement':              pd.concat(P),
        'caller_within_type_depth_bias':          pd.concat(B),
        'caller_agreement_depth_and_type_matched': pd.concat(M),
    }
    # the stratified-null test is a Phase 8 addition and has no 2-section
    # counterpart to overwrite, so it is written unconditionally alongside
    pd.concat(G).to_csv(out + 'caller_agreement_matched_significance' + suffix + '.csv',
                        index=False)
    print('wrote caller_agreement_matched_significance' + suffix + '.csv', flush=True)
    for k, v in outs.items():
        v.to_csv(out + k + suffix + '.csv', index=False)
        print('wrote', k + suffix + '.csv', len(v), 'rows', flush=True)
    return outs


def verify():
    """Re-run the two original sections and require an EXACT match to the committed tables."""
    outs = run_set([('sham', 'sham', 'sham'), ('sbr', 'SBR', 'sbr')], '_verify2sec')
    bad = 0
    for k, v in outs.items():
        ref = pd.read_csv(OUT + k + '.csv')
        new = pd.read_csv(OUT + k + '_verify2sec.csv')
        if list(ref.columns) != list(new.columns) or len(ref) != len(new):
            print('MISMATCH shape/columns:', k, ref.shape, new.shape); bad += 1; continue
        eq = ((ref == new) | (ref.isna() & new.isna())).all().all()
        print(('OK   ' if eq else 'DIFF ') + k, ref.shape)
        if not eq:
            bad += 1
            d = (ref != new) & ~(ref.isna() & new.isna())
            print(ref[d.any(1)].to_string())
    print('VERIFY:', 'PASS' if bad == 0 else 'FAIL (%d tables differ)' % bad)
    return bad


# the two-section set, in the order the published tables carry it.  This is the spec that
# produced results/phase3/caller_*_2sec_c6.csv -- the post-C6 base of
# summarize_caller_coverage.py:67 -- which until now had no committed producer at all
# (logs/m1_callers_2sec_c6.log is the only record that it ran).
TWO_SEC = [('sham', 'sham', 'sham'), ('sbr', 'SBR', 'sbr')]

SETS = {
    '2sec_c6':    (TWO_SEC, '_2sec_c6'),
    '2sec':       (TWO_SEC, ''),
    '11sections': ([(s, a, s) for s, a, _ in SECTIONS], '_11sections'),
}

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--verify', action='store_true',
                    help='re-run the two original sections and require an exact match')
    ap.add_argument('--all', action='store_true', help="alias for --set 11sections")
    ap.add_argument('--set', choices=sorted(SETS), help='which section set to write')
    ap.add_argument('--suffix', help='override the output suffix for --set')
    ap.add_argument('--out-dir', default=OUT,
                    help='write elsewhere (trailing slash added); default %s' % OUT)
    ap.add_argument('--allow-missing-deepscence', action='store_true',
                    help='proceed with a caller set that has no DeepScence column')
    a = ap.parse_args()
    out = a.out_dir if a.out_dir.endswith('/') else a.out_dir + '/'
    if a.verify:
        sys.exit(verify())
    name = '11sections' if a.all else a.set
    if not name:
        ap.error('nothing to do: pass --verify, --all, or --set {%s}' % ','.join(sorted(SETS)))
    specs, suffix = SETS[name]
    if a.suffix is not None:
        suffix = a.suffix
    run_set(specs, suffix, out=out, require_deepscence=not a.allow_missing_deepscence)
