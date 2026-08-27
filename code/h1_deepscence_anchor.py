#!/usr/bin/env python3
"""Phase 9 — the DeepScence sign-anchor diagnostic on H1 (PREREG §8 prediction P-ii, §3.9 D3).

`code/deepscence_reanchor.py` is the mouse producer; it is hard-wired to the mouse sections,
mouse symbols and `data/processed/`.  This is the same diagnostic on H1, importing
`partial_spearman` from that file so the definition cannot drift: Spearman of two variables
after linearly removing the RANK of log transcript counts from both.  On a 5K panel every
anchor candidate is partly a detection-rate readout, so the polarity decision is taken on the
depth-partialled correlation, never the raw one.

Reported per section:
  rho_raw / rho_partial of the score against CDKN1A counts (the PUBLISHED anchor)
  the same against the 8-gene proliferation anchor (D3 primary alternative; all eight are
  human symbols already and are absent from every A_*.txt and B_*.txt -- asserted at run time)
  the same against LMNB1 (secondary only -- it is inside two Tier B modules, P12)
  fold sign stability at k = 20 random folds, the P-ii falsifier (< 0.90 in >= 1 of 7)
  the rank of CDKN1A among the on-panel CoreScence genes by correlation with the score,
  which is DeepScence's own rule reproduced as a diagnostic

Usage: python3 code/h1_deepscence_anchor.py
Writes results/phase9_h1/deepscence_anchor_h1.csv
"""
import sys, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr, rankdata
sys.path.insert(0, "/workspace/code")
import h1_common as H
from deepscence_reanchor import partial_spearman, CORE

PROLIF = ["KIF20A", "NCAPH", "ANLN", "ECT2", "GTSE1", "UHRF1", "FEN1", "CLSPN"]
EXCL = {"Low_quality", "Unknown"}


def _assert_disjoint():
    import glob
    ab = set()
    for p in glob.glob(H.GS_HUMAN + "/A_*.txt") + glob.glob(H.GS_HUMAN + "/B_*.txt"):
        ab |= {l.strip() for l in open(p) if l.strip()}
    bad = [g for g in PROLIF if g in ab]
    assert not bad, "proliferation anchor overlaps Tier A/B: %s" % bad
    return len(ab)


def z(v):
    v = np.asarray(v, float); s = v.std()
    return (v - v.mean()) / (s if s > 1e-12 else 1.0)


def main():
    _assert_disjoint()
    # DeepScence's own scoring set is coreGS_v2 filtered to occurrence >= 5 (io.get_geneset,
    # default n=5) -- 39 genes, 33 of them on the H1 panel.  Reading the file unfiltered
    # would give 1,303 genes and a different diagnostic.
    core = pd.read_csv(CORE, index_col=0)
    core = core[core["occurrence"] >= 5]
    core_genes = set(core["gene_symbol"].astype(str)) if "gene_symbol" in core.columns \
        else set(core.index.astype(str))
    rows = []
    for sec in H.ALL_SECTIONS:
        dsp = H.PROC + "/deepscence_h1_%s.csv" % sec
        if not os.path.exists(dsp):
            print("skip", sec, "(no DeepScence score)"); continue
        ds = pd.read_csv(dsp).set_index("cell_id")
        sen = pd.read_csv(H.PROC + "/senders_h1_%s.csv" % sec).set_index("cell_id")
        X, names, bc, _ = H.load_matrix(sec, "gene")
        want = sorted((set(PROLIF) | {"LMNB1", "CDKN1A"} | core_genes) & set(names))
        cols = [int(np.where(names == g)[0][0]) for g in want]
        sub = pd.DataFrame(np.asarray(X[:, cols].todense()), index=bc, columns=want)
        tot = pd.Series(np.asarray(X.sum(1)).ravel(), index=bc)
        sub = sub.reindex(ds.index); tot = tot.reindex(ds.index)
        s = ds.deepscence_score.to_numpy(float)
        L = np.log1p(sub.to_numpy(float) / np.maximum(tot.to_numpy(float), 1)[:, None] * 1e4)
        L = pd.DataFrame(L, index=sub.index, columns=want)
        lc = np.log1p(tot.to_numpy(float))
        cc = sen.reindex(ds.index)
        cdk = cc.cdkn1a_counts.to_numpy(float)
        prolif = np.mean([z(L[g]) for g in PROLIF if g in L], 0)
        lmnb1 = L["LMNB1"].to_numpy(float) if "LMNB1" in L else np.full(len(L), np.nan)

        keep = ~cc.cell_type.isin(EXCL)
        cons = np.zeros(len(cc)); nk = 0
        for col in ["tierA_score", "senepy_score", "cdkn1a_counts"]:
            v = cc[col].to_numpy(float); r = np.full(len(v), np.nan)
            for t in cc.cell_type.dropna().unique():
                m = (cc.cell_type == t).to_numpy() & np.isfinite(v) & keep.to_numpy()
                if m.sum() < 50:
                    continue
                r[m] = rankdata(v[m]) / m.sum()
            cons = cons + np.nan_to_num(r - 0.5); nk += 1
        cons = cons / nk

        cg = [g for g in want if g in core_genes]
        cr = np.array([pearsonr(L[g].to_numpy(float), s)[0] if L[g].std() > 0 else 0.0
                       for g in cg])
        order = np.argsort(-cr)
        rank_cdkn1a = (int(np.where(np.array(cg)[order] == "CDKN1A")[0][0])
                       if "CDKN1A" in cg else -1)

        def r_(v):
            m = np.isfinite(v) & np.isfinite(s)
            return float(spearmanr(v[m], s[m]).statistic) if m.sum() > 100 else np.nan

        def fold_stability(v, k=20, seed=0):
            m = np.isfinite(v) & np.isfinite(s)
            vv, ss, ll = np.asarray(v)[m], s[m], lc[m]
            if vv.size < 2000:
                return np.nan
            g = np.sign(partial_spearman(vv, ss, ll))
            f = np.random.default_rng(seed).permutation(vv.size) % k
            sg = [np.sign(partial_spearman(vv[f == i], ss[f == i], ll[f == i]))
                  for i in range(k)]
            return float(np.mean(np.array(sg) == g))

        row = dict(section=sec, n_cells=int(len(ds)),
                   rho_score_vs_counts=round(r_(tot.to_numpy(float)), 4),
                   rho_raw_cdkn1a=round(r_(cdk), 4),
                   rho_partial_cdkn1a=round(partial_spearman(cdk, s, lc), 4),
                   stab_cdkn1a=round(fold_stability(cdk), 3),
                   rho_partial_prolif=round(partial_spearman(prolif, s, lc), 4),
                   stab_prolif=round(fold_stability(prolif), 3),
                   rho_partial_lmnb1=round(partial_spearman(lmnb1, s, lc), 4),
                   stab_lmnb1=round(fold_stability(lmnb1), 3),
                   rho_partial_consensus=round(partial_spearman(cons, s, lc), 4),
                   stab_consensus=round(fold_stability(cons), 3),
                   n_core_on_panel=len(cg), rank_cdkn1a_in_core=rank_cdkn1a)
        rows.append(row); print(row, flush=True)
        del X, sub
    d = pd.DataFrame(rows)
    d.to_csv(H.RESULTS + "/deepscence_anchor_h1.csv", index=False)
    pd.set_option("display.width", 260)
    print(d.to_string(index=False))
    print("\nP-ii falsifier: stability >= 0.90 in ALL 7 sections would falsify.  "
          "CDKN1A stability: %s" % d.stab_cdkn1a.tolist())
    print("wrote", H.RESULTS + "/deepscence_anchor_h1.csv")


if __name__ == "__main__":
    main()
