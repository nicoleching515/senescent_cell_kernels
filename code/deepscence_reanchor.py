#!/usr/bin/env python3
"""Phase 8, task 8.6 -- D3: re-anchor DeepScence off CDKN1A.

WHY.  DeepScence fixes the sign of its bottleneck node by correlating it with
CDKN1A (`DeepScence/io.py::fix_score_direction`: if CDKN1A's rank in the
score-vs-gene correlation table falls in the bottom half, negate the score).
That creates two problems the paper has to own:
  * CIRCULARITY with the `Cdkn1a`+ caller (BIO_PHASE3 §4.4);
  * a polarity that FLIPS between sections -- correlation of `ds` with its own
    gene set was -0.350 in sham 7250 and +0.318 in SBR 7259.

WHAT RE-ANCHORING ACTUALLY IS.  `fix_score_direction` picks the node by
`corr_df["correlation"].abs().mean()`, which is INVARIANT to negating the score,
and then applies a single global sign to that node.  So the anchor decides
exactly one bit per section and nothing else: re-anchoring is a per-section sign
choice applied to the score already computed.  No DeepScence re-run is needed,
and the D1 outputs are not touched (8.3/8.4 must stay a pure coverage comparison).

THE FOUR RULES REPORTED
  published   as shipped: DeepScence's own CDKN1A anchor (the sign in the D1 files)
  prolif      sign so that corr(ds, proliferation score) < 0.  Senescence is
              proliferation arrest, so a senescence score must run opposite to a
              proliferation score.  The 8 proliferation genes used are on the
              mouse panel and in NEITHER any Tier A file NOR any Tier B module.
  lmnb1       sign so that corr(ds, Lmnb1) < 0.  REPORTED BUT NOT PRIMARY:
              the §7 brief proposes `Lmnb1`, but `Lmnb1` is a member of
              B_downstream_arrest and B_secondary_senescence (and of the
              non-strict Tier A variants), so it fails the brief's own
              "in neither Tier A nor any Tier B module" condition.
  consensus   sign so that corr(ds, mean within-type rank of the other three
              callers) > 0.  REPORTED BUT NOT PRIMARY: anchoring on the other
              callers makes the caller-agreement statistic partly circular in
              the opposite direction, exactly the failure mode D3 exists to fix.

SIGN-INVARIANT SUMMARY.  |ds| is emitted as its own caller, so a top-5% call can
be made by magnitude and compared against the signed calls.  A sign-invariant
statistic is the only kind the anchor cannot move.

Writes  data/processed/deepscence_d3_<section>.csv   (never overwrites D1)
        results/phase3/deepscence_anchor_decisions.csv
"""
import sys, os, csv, glob, numpy as np, pandas as pd, h5py
from scipy.sparse import csc_matrix
from scipy.stats import pearsonr, spearmanr, rankdata

RAW = '/workspace/data/raw/'; PROC = '/workspace/data/processed/'
RES = '/workspace/results/phase3/'; GS = '/workspace/genesets/'
ORTH = GS + 'mouse_human_orthologs_MGI.csv'
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _pkgdata   # resolves DeepScence/senepy package data; see AUDIT_REPRODUCIBILITY D2
CORE = _pkgdata.core_gs(verbose=False)
EXCL = {'Low_quality', 'Unknown'}

# proliferation anchor: on the mouse panel, and absent from every A_*.txt and
# every B_*.txt in genesets/ (verified at run time by _check_disjoint).
PROLIF = ['Kif20a', 'Ncaph', 'Anln', 'Ect2', 'Gtse1', 'Uhrf1', 'Fen1', 'Clspn']

SECTIONS = ['7001_liver_sham_Male_52-U1', '7239_liver_sbr_Male_52-U1',
            '7248_liver_sham_Male_26-U1', '7250_liver_sham_Male_26-U1',
            '7259_liver_sbr_Male_26-U1', '7260_liver_sbr_Male_26-U1',
            '7352_liver_sham_Male_2-U1', '7361_liver_sbr_Male_2-U1',
            '7435_liver_sham_Male_10-U1', '7448_liver_sbr_Male_10-U1',
            '7450_liver_sbr_Male_10-U1']
DS_ALIAS = {'7250_liver_sham_Male_26-U1': 'sham', '7259_liver_sbr_Male_26-U1': 'sbr'}


def _check_disjoint():
    """genesets/ is read-only here; this only reads, and it fails loudly."""
    ab = set()
    for f in glob.glob(GS + 'A_*.txt') + glob.glob(GS + 'B_*.txt'):
        ab |= {l.strip() for l in open(f) if l.strip()}
    bad = [g for g in PROLIF if g in ab]
    assert not bad, 'proliferation anchor overlaps Tier A/B: %s' % bad
    return len(ab)


def core_mouse_genes():
    orth = {r['mouse_symbol']: r['human_symbol'] for r in csv.DictReader(open(ORTH))}
    h2m = {}
    for m, h in orth.items():
        h2m.setdefault(h, m)
    core = pd.read_csv(CORE)
    core['occurrence'] = pd.to_numeric(core.occurrence, errors='coerce')
    core = core[core.occurrence >= 5]
    return {h2m[g]: d for g, d in zip(core.gene_symbol, core.direction) if g in h2m}


def load_genes(section, genes):
    """Per-cell counts for `genes`, for every barcode in the h5, as a DataFrame."""
    f = h5py.File(RAW + section + '/cell_feature_matrix.h5', 'r')
    ft = np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    nm = np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc = np.array([x.decode() for x in f['matrix/barcodes'][:]])
    want = [g for g in genes if g in set(nm[ft == 'Gene Expression'])]
    rows = [int(np.where((nm == g) & (ft == 'Gene Expression'))[0][0]) for g in want]
    M = csc_matrix((f['matrix/data'][:].astype(np.float32),
                    f['matrix/indices'][:].astype(np.int32),
                    f['matrix/indptr'][:].astype(np.int64)),
                   shape=tuple(f['matrix/shape'][:]))
    X = np.asarray(M[rows, :].todense()).T
    tot = np.asarray(M[np.where(ft == 'Gene Expression')[0], :].sum(0)).ravel()
    return pd.DataFrame(X, index=bc, columns=want), pd.Series(tot, index=bc)


def partial_spearman(x, y, zz):
    """Spearman of x and y after linearly removing the RANK of zz from both.
    Every anchor candidate on this panel is a detection-rate readout as much as a
    biological one -- proliferation genes are detected in 0.2-1.6 % of cells and
    detection scales with depth -- so an anchor decided on the raw correlation is
    partly deciding on sequencing depth.  Signs are taken from THIS, not the raw
    correlation."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(zz)
    rx, ry, rz = (rankdata(v[m]) for v in (x, y, zz))

    def resid(a, b):
        B = np.c_[np.ones(len(b)), b]
        return a - B @ np.linalg.lstsq(B, a, rcond=None)[0]
    return float(np.corrcoef(resid(rx, rz), resid(ry, rz))[0, 1])


def z(v):
    v = np.asarray(v, float)
    s = v.std()
    return (v - v.mean()) / (s if s > 1e-12 else 1.0)


def main():
    n_ab = _check_disjoint()
    core = core_mouse_genes()
    rows = []
    for sec in SECTIONS:
        ds = pd.read_csv(PROC + 'deepscence_%s.csv' % DS_ALIAS.get(sec, sec)).set_index('cell_id')
        sen = pd.read_csv(PROC + 'senders_%s.csv' % sec).set_index('cell_id')
        need = sorted(set(PROLIF) | {'Lmnb1', 'Cdkn1a'} | set(core))
        X, tot = load_genes(sec, need)
        X = X.reindex(ds.index); tot = tot.reindex(ds.index)
        s = ds.deepscence_score.to_numpy(float)
        lcpm = np.log1p(X.to_numpy(float) / np.maximum(tot.to_numpy(float), 1)[:, None] * 1e4)
        L = pd.DataFrame(lcpm, index=X.index, columns=X.columns)

        prolif = np.mean([z(L[g]) for g in PROLIF if g in L], 0)
        lmnb1 = L['Lmnb1'].to_numpy(float) if 'Lmnb1' in L else np.full(len(L), np.nan)

        # caller consensus: mean of within-cell-type ranks of the other three callers
        cc = sen.reindex(ds.index)
        keep = ~cc.cell_type.isin(EXCL)
        cons = np.zeros(len(cc)); nk = 0
        for col in ['tierA_score', 'senepy_score', 'cdkn1a_counts']:
            v = cc[col].to_numpy(float)
            r = np.full(len(v), np.nan)
            for t in cc.cell_type.dropna().unique():
                m = (cc.cell_type == t).to_numpy() & np.isfinite(v) & keep.to_numpy()
                if m.sum() < 50:
                    continue
                r[m] = rankdata(v[m]) / m.sum()
            cons = cons + np.nan_to_num(r - 0.5); nk += 1
        cons = cons / nk

        # DeepScence's own rule, reproduced as a diagnostic: rank of Cdkn1a among
        # the coreGS(occ>=5) genes' correlations with the score, descending.
        cg = [g for g in core if g in L.columns]
        cr = np.array([pearsonr(L[g].to_numpy(float), s)[0] if L[g].std() > 0 else 0.0
                       for g in cg])
        order = np.argsort(-cr)
        rank_cdkn1a = int(np.where(np.array(cg)[order] == 'Cdkn1a')[0][0]) if 'Cdkn1a' in cg else -1

        def r_(v):
            m = np.isfinite(v) & np.isfinite(s)
            return float(spearmanr(v[m], s[m]).statistic) if m.sum() > 100 else np.nan

        lcv = np.log1p(tot.to_numpy(float))

        def fold_stability(v, k=20, seed=0):
            """Fraction of k random equal folds whose sign of rho matches the
            whole-section sign.  An anchor that decides the polarity RELIABLY
            decides it the same way in every fold; one that is riding on noise
            does not.  This is the diagnostic that separates "small but real"
            from "an arbitrary bit"."""
            m = np.isfinite(v) & np.isfinite(s)
            vv, ss = np.asarray(v)[m], s[m]
            if vv.size < 2000:
                return np.nan
            ll = lcv[m]
            g = np.sign(partial_spearman(vv, ss, ll))
            f = np.random.default_rng(seed).permutation(vv.size) % k
            sg = [np.sign(partial_spearman(vv[f == i], ss[f == i], ll[f == i]))
                  for i in range(k)]
            return float(np.mean(np.array(sg) == g))

        r_prolif, r_lmnb1, r_cons = r_(prolif), r_(lmnb1), r_(cons)
        cdk = cc.cdkn1a_counts.to_numpy(float)
        r_cdkn1a = r_(cdk)
        lc = np.log1p(tot.to_numpy(float))
        p_prolif = partial_spearman(prolif, s, lc)
        p_lmnb1 = partial_spearman(lmnb1, s, lc)
        p_cons = partial_spearman(cons, s, lc)
        p_cdkn1a = partial_spearman(cdk, s, lc)
        # each rule's sign for the score as it currently stands (+1 = keep, -1 = flip),
        # taken from the DEPTH-PARTIALLED correlation
        sgn = dict(published=1,
                   prolif=-1 if p_prolif > 0 else 1,
                   lmnb1=-1 if p_lmnb1 > 0 else 1,
                   consensus=-1 if p_cons < 0 else 1)
        out = pd.DataFrame({'cell_id': ds.index,
                            'ds_published': np.round(s, 5),
                            'ds_prolif_anchor': np.round(sgn['prolif'] * s, 5),
                            'ds_lmnb1_anchor': np.round(sgn['lmnb1'] * s, 5),
                            'ds_consensus_anchor': np.round(sgn['consensus'] * s, 5),
                            'ds_abs': np.round(np.abs(s), 5)})
        out.to_csv(PROC + 'deepscence_d3_%s.csv' % sec, index=False)
        rows.append(dict(section=sec, arm='SBR' if 'sbr' in sec else 'sham',
                         n_cells=len(ds),
                         rho_ds_prolif=round(r_prolif, 4),
                         rho_ds_Lmnb1=round(r_lmnb1, 4),
                         rho_ds_caller_consensus=round(r_cons, 4),
                         rho_ds_Cdkn1a=round(r_cdkn1a, 4),
                         prho_ds_prolif=round(p_prolif, 4),
                         prho_ds_Lmnb1=round(p_lmnb1, 4),
                         prho_ds_caller_consensus=round(p_cons, 4),
                         prho_ds_Cdkn1a=round(p_cdkn1a, 4),
                         stab_prolif=fold_stability(prolif),
                         stab_Lmnb1=fold_stability(lmnb1),
                         stab_consensus=fold_stability(cons),
                         stab_Cdkn1a=fold_stability(cdk),
                         det_pct_prolif_any=round(100 * float(
                             (X[[g for g in PROLIF if g in X]].to_numpy() > 0).any(1).mean()), 3),
                         det_pct_Lmnb1=round(100 * float((X['Lmnb1'] > 0).mean()), 3),
                         det_pct_Cdkn1a=round(100 * float((X['Cdkn1a'] > 0).mean()), 3),
                         rho_tierA_prolif=round(float(spearmanr(
                             prolif, np.nan_to_num(cc.tierA_score.to_numpy(float))).statistic), 4),
                         cdkn1a_rank_in_coreGS=rank_cdkn1a, n_coreGS_on_panel=len(cg),
                         sign_published=1, sign_prolif=sgn['prolif'],
                         sign_lmnb1=sgn['lmnb1'], sign_consensus=sgn['consensus'],
                         prolif_agrees_with_published=sgn['prolif'] == 1,
                         lmnb1_agrees_with_published=sgn['lmnb1'] == 1,
                         consensus_agrees_with_published=sgn['consensus'] == 1))
        print('%s  partial rho | prolif %+.4f  Lmnb1 %+.4f  consensus %+.4f  '
              'Cdkn1a %+.4f | Cdkn1a rank %d/%d | signs P%+d L%+d C%+d'
              % (sec, p_prolif, p_lmnb1, p_cons, p_cdkn1a, rank_cdkn1a, len(cg),
                 sgn['prolif'], sgn['lmnb1'], sgn['consensus']), flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(RES + 'deepscence_anchor_decisions.csv', index=False)
    pd.set_option('display.width', 250)
    print()
    print(D.to_string(index=False))
    print('\n(%d genes across all A_*/B_* files checked for anchor disjointness)' % n_ab)


if __name__ == '__main__':
    main()
