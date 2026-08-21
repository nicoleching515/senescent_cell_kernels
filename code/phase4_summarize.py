"""Phase 4 aggregation — Figure 4's numbers.

Two quantities per (method, LR pair, null), and the whole mechanistic verdict
lives in the difference between them:

  SIG SURVIVAL   fraction of interactions called significant on REAL
                 coordinates that are still called significant on SHUFFLED
                 coordinates.  This is the CellWHISPER number.
  SCORE SF       median over replicates of (null score / real score), the
                 project's own surviving-fraction metric (CS_PHASE3 s0).

  SCORE SF ~ 1 and SIG SURVIVAL ~ 1  ->  the null is too weak: the shuffled
      tissue really does still contain the signal, so rejecting is correct
      behaviour on an uninformative null.  (This is how OUR estimator fails.)
  SCORE SF ~ 0 but SIG SURVIVAL ~ 1  ->  the null is fine and the method's own
      significance test is miscalibrated: it keeps calling interactions after
      the quantity it is testing has collapsed.
"""
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd

P4 = "/workspace/results/phase4"
PARTS = f"{P4}/parts"
ALPHA = 0.05
NULLS = ["N3_lig", "N4_lig", "N3_type", "N0_perm", "N0_type"]
METHOD_ORDER = ["COMMOT", "CellChat v2*", "SpaTalk*", "NCEM linear*"]
METHOD_MAP = {"commot": "COMMOT", "cellchat": "CellChat v2*",
              "spatalk": "SpaTalk*", "ncem": "NCEM linear*"}
PAIR_ORDER = ["Ccl2->Ccr2", "Tnf->Tnfrsf1a/1b", "Tgfb1->Tgfbr1/2", "Il1a->Il1r1"]


def load():
    fs = [f for f in glob.glob(f"{PARTS}/*.csv")
          if "ncem_sweep" not in f and os.path.getsize(f) > 0]
    ds = []
    for f in fs:
        try:
            d = pd.read_csv(f)
        except pd.errors.EmptyDataError:
            continue          # a job whose ligand+ count fell below the floor
        if len(d):
            ds.append(d)
    df = pd.concat(ds, ignore_index=True)
    df["sig"] = df["p"] < ALPHA
    df["method"] = df["method"].map(METHOD_MAP)
    return df


def build(df, drop_diagonal=False):
    if drop_diagonal:
        df = df[df.sender != df.receiver]
    key = ["tile", "section", "arm", "method", "pair", "sender", "receiver"]
    real = df[df.cond == "real"].set_index(key)[["score", "sig", "p"]]
    real.columns = ["real_score", "real_sig", "real_p"]
    out = []
    for null in NULLS:
        nd = df[df.cond == null]
        if nd.empty:
            continue
        g = nd.groupby(key).agg(
            n_rep=("rep", "size"),
            null_sig_frac=("sig", "mean"),
            null_score_med=("score", "median"),
            null_score_mean=("score", "mean"),
            null_score_q05=("score", lambda v: np.nanpercentile(v, 5)),
            null_score_q95=("score", lambda v: np.nanpercentile(v, 95)),
            null_p_med=("p", "median"))
        g["null"] = null
        out.append(g.join(real, how="inner"))
    R = pd.concat(out).reset_index()
    with np.errstate(divide="ignore", invalid="ignore"):
        R["score_sf"] = np.where(np.abs(R.real_score) > 0,
                                 R.null_score_med / R.real_score, np.nan)
    return R


def headline(R):
    rows = []
    for (m, null), g in R.groupby(["method", "null"]):
        for pair in ["ALL"] + PAIR_ORDER:
            gg = g if pair == "ALL" else g[g.pair == pair]
            if gg.empty:
                continue
            s = gg[gg.real_sig]
            rows.append(dict(
                method=m, null=null, pair=pair,
                n_interactions=len(gg), n_real_sig=len(s),
                real_sig_rate=float(gg.real_sig.mean()),
                sig_survival=float(s.null_sig_frac.mean()) if len(s) else np.nan,
                sig_survival_all=float(gg.null_sig_frac.mean()),
                score_sf_median=float(np.nanmedian(s.score_sf)) if len(s) else np.nan,
                score_sf_q25=float(np.nanpercentile(s.score_sf, 25)) if len(s) else np.nan,
                score_sf_q75=float(np.nanpercentile(s.score_sf, 75)) if len(s) else np.nan,
                n_rep_median=float(gg.n_rep.median())))
        # false-discovery-style rate: significance rate under the null over ALL
        # interactions, not just the real-significant ones
    return pd.DataFrame(rows)


def verdict(H):
    """Classify each (method, null) by the PAIR of numbers, which is the whole
    mechanistic point (report section 2.5)."""
    rows = []
    for _, r in H[H.pair == "ALL"].iterrows():
        sf, sv = r.score_sf_median, r.sig_survival
        if not np.isfinite(sf) or not np.isfinite(sv):
            v = "n/a"
        elif sv <= 0.20 and (not np.isfinite(sf) or sf <= 0.35):
            v = "this null bites: score and calls both collapse"
        elif sv <= 0.20:
            v = "this null bites: calls collapse"
        elif sf <= 0.35:
            v = "MISCALIBRATED TEST: the score collapsed, the calls did not"
        elif r.null == "N0_perm":
            v = ("STATISTIC IS NOT SPATIAL: survives complete destruction of "
                 "the spatial arrangement")
        elif r.null == "N0_type":
            v = ("SURVIVES CELLWHISPER'S OWN CONTROL: within-cell-type location "
                 "permutation, which keeps each cell type's spatial "
                 "organisation and destroys only ligand-receptor proximity")
        else:
            v = ("null does not bite: the shuffled tissue still produces the "
                 "score the method measures")
        rows.append(dict(method=r.method, null=r.null,
                         score_sf=round(float(sf), 3),
                         sig_survival=round(float(sv), 3),
                         null_sig_rate=round(float(r.sig_survival_all), 4),
                         real_sig_rate=round(float(r.real_sig_rate), 4),
                         n_interactions=int(r.n_interactions),
                         n_real_sig=int(r.n_real_sig),
                         n_rep=int(r.n_rep_median), verdict=v))
    return pd.DataFrame(rows)


def main():
    df = load()
    print("parts loaded:", df.method.value_counts().to_dict())
    R = build(df)
    R.to_csv(f"{P4}/interactions.csv.gz", index=False)
    H = headline(R)
    H.to_csv(f"{P4}/headline.csv", index=False)
    # by surgical arm, and by section, as a robustness check
    rows = []
    for keycol in ("arm", "section"):
        for (m, null, k), g in R.groupby(["method", "null", keycol]):
            sgood = g[g.real_sig]
            rows.append(dict(split=keycol, level=k, method=m, null=null,
                             n_interactions=len(g), n_real_sig=len(sgood),
                             real_sig_rate=float(g.real_sig.mean()),
                             sig_survival=float(sgood.null_sig_frac.mean())
                             if len(sgood) else np.nan,
                             score_sf_median=float(np.nanmedian(sgood.score_sf))
                             if len(sgood) else np.nan))
    pd.DataFrame(rows).to_csv(f"{P4}/headline_by_split.csv", index=False)
    Rd = build(df, drop_diagonal=True)
    headline(Rd).to_csv(f"{P4}/headline_offdiagonal.csv", index=False)
    V = verdict(H)
    V.to_csv(f"{P4}/verdict.csv", index=False)
    piv = H[H.pair == "ALL"].pivot(index="method", columns="null",
                                   values="sig_survival")
    print("\nSIG SURVIVAL (fraction of real-significant interactions still "
          "significant under the null)\n")
    print(piv.reindex(METHOD_ORDER)[NULLS].round(3).to_string())
    piv2 = H[H.pair == "ALL"].pivot(index="method", columns="null",
                                    values="score_sf_median")
    print("\nSCORE SURVIVING FRACTION (median null score / real score)\n")
    print(piv2.reindex(METHOD_ORDER)[NULLS].round(3).to_string())
    print("\nRANK CORRELATION between the REAL score and the SHUFFLED score,\n"
          "across interactions (Spearman).  ~1 means the shuffle changed the\n"
          "ORDERING of interactions not at all, i.e. the statistic is not\n"
          "reading the geometry that was destroyed.\n")
    from scipy.stats import spearmanr
    cr = []
    for (m, null), g in R.groupby(["method", "null"]):
        g = g[np.isfinite(g.real_score) & np.isfinite(g.null_score_med)]
        rho = spearmanr(g.real_score, g.null_score_med).statistic if len(g) > 10 else np.nan
        cr.append(dict(method=m, null=null, spearman=rho, n=len(g)))
    CR = pd.DataFrame(cr)
    CR.to_csv(f"{P4}/score_rank_correlation.csv", index=False)
    print(CR.pivot(index="method", columns="null", values="spearman")
          .reindex(METHOD_ORDER)[NULLS].round(3).to_string())
    print("\nSIGNIFICANCE RATE OVER ALL INTERACTIONS (real vs each null)\n")
    A = H[H.pair == "ALL"]
    tab = A.pivot(index="method", columns="null", values="sig_survival_all")
    tab.insert(0, "real", A.groupby("method").real_sig_rate.first())
    print(tab.reindex(METHOD_ORDER)[["real"] + NULLS].round(3).to_string())
    print("\nREAL-DATA SIGNIFICANCE RATE\n")
    print(H[H.pair == "ALL"].groupby("method").real_sig_rate.first()
          .reindex(METHOD_ORDER).round(4).to_string())
    print("\nVERDICT\n")
    print(V.to_string(index=False))


if __name__ == "__main__":
    main()
