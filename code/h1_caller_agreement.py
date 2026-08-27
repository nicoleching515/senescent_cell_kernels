#!/usr/bin/env python3
"""Phase 9 Job B step 4 — caller agreement on H1, conditioned on cell type and transcript-depth
decile, exactly as `reports/BIO_PHASE3.md` §4.4 does for mouse.

The estimator is not reimplemented: `_within_type_flags`, `_matched_flags`,
`_stratified_null` and `_within_type_depth_quintile` are imported from
`code/caller_disagree_all.py` and called on an H1 dataframe.  The primary rule is the
DEPTH- AND TYPE-MATCHED one (PREREG §3.7 b): top-5 % recomputed inside each
(cell type x within-type transcript-count decile), strata >= 50 cells, N_DEPTH_DECILES = 10.

Also computed, because the pre-registration names them as §8 predictions:
  P-i    Spearman rho of deepscence_score against transcript counts, per section
  P-iii  the SIGN-INVARIANT (|score|) matched agreement of Tier A with DeepScence
  P-v    within-type depth quintile enrichment Q5/Q1 for every caller

Usage: python3 code/h1_caller_agreement.py [SPLN07 ...]
Writes results/phase9_h1/caller_*.csv
"""
import sys, os, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
sys.path.insert(0, "/workspace/code")
import h1_common as H
from caller_disagree_all import (_within_type_flags, _matched_flags, _stratified_null,
                                 _within_type_depth_quintile, SCORES, MIN_STRATUM)

PROC = H.PROC + "/"
EXCL = set(H.EXCLUDE_TYPES)


def load(section):
    sen = pd.read_csv(PROC + "senders_h1_%s.csv" % section).set_index("cell_id")
    ds = PROC + "deepscence_h1_%s.csv" % section
    if os.path.exists(ds):
        sen = sen.join(pd.read_csv(ds).set_index("cell_id"))
    cells = H.cells_table(section).set_index("cell_id")
    df = sen.join(cells[["transcript_counts", "cell_area", "nucleus_area",
                         "segmentation_method"]])
    return df[~df.cell_type.isin(EXCL)].copy()


def tables(section):
    df = load(section)
    scores = [s for s in SCORES if s in df and df[s].notna().any()]
    bias, pair, sig, tech = [], [], [], []
    df["wq"] = _within_type_depth_quintile(df)
    for s in scores:
        f = _within_type_flags(df, s)
        e = {}
        for qq in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
            m = (df.wq == qq).to_numpy(); bg = m.mean()
            e[qq] = float((f & m).sum() / max(f.sum(), 1) / bg) if bg else np.nan
            bias.append(dict(section=section, caller=s, within_type_depth_quintile=qq,
                             bg_pct=round(100 * bg, 2), enrichment=round(e[qq], 3)))
        v = df[s].to_numpy(float); ok = np.isfinite(v)
        rho = spearmanr(v[ok], df.transcript_counts.to_numpy(float)[ok]).statistic
        tech.append(dict(section=section, caller=s, n_scored=int(ok.sum()),
                         frac_cells_scored=round(float(ok.mean()), 4),
                         spearman_vs_transcript_counts=round(float(rho), 4),
                         q5_over_q1=round(e["Q5"] / e["Q1"], 3) if e.get("Q1") else None))
    # add |deepscence_score| as a fifth, sign-invariant score (PREREG §3.9 / P-iii)
    if "deepscence_score" in df:
        df["abs_deepscence_score"] = df.deepscence_score.abs()
        scores = scores + ["abs_deepscence_score"]
    for i in range(len(scores)):
        for j in range(i + 1, len(scores)):
            a_s, b_s = scores[i], scores[j]
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
            pair.append(dict(section=section, A=a_s, B=b_s, n=int(len(sub)),
                             n_A=int(fa.sum()), n_B=int(fb.sum()), n_both=inter,
                             jaccard=round(float(jac), 5), chance=round(float(ch), 5),
                             ratio=round(float(jac / ch), 3) if ch else None))
            assert (sa == sb).all()
            e_, sd, nst = _stratified_null(fa, fb, sa)
            sig.append(dict(section=section, A=a_s, B=b_s, n=int(len(sub)), n_strata=nst,
                            n_A=int(fa.sum()), n_B=int(fb.sum()), n_both=inter,
                            exp_both_marginal=round(float(fa.sum()) * fb.sum() / len(sub), 2),
                            ratio_marginal=round(float(jac / ch), 3) if ch else None,
                            exp_both_stratified=round(e_, 2),
                            sd_both_stratified=round(sd, 3),
                            ratio_stratified=round(inter / e_, 3) if e_ else None,
                            z=round((inter - e_) / sd, 2) if sd else None))
    return pd.DataFrame(bias), pd.DataFrame(pair), pd.DataFrame(sig), pd.DataFrame(tech)


if __name__ == "__main__":
    secs = sys.argv[1:] or [s for s in H.ALL_SECTIONS
                            if os.path.exists(PROC + "senders_h1_%s.csv" % s)]
    B, P_, S, T = [], [], [], []
    for s in secs:
        print("...", s, flush=True)
        b, p, g, t = tables(s); B.append(b); P_.append(p); S.append(g); T.append(t)
    pd.concat(B).to_csv(H.RESULTS + "/caller_within_type_depth_bias.csv", index=False)
    pd.concat(P_).to_csv(H.RESULTS + "/caller_agreement_depth_and_type_matched.csv", index=False)
    sig = pd.concat(S); sig.to_csv(H.RESULTS + "/caller_agreement_matched_significance.csv",
                                   index=False)
    pd.concat(T).to_csv(H.RESULTS + "/caller_technical_loading.csv", index=False)
    # pooled, as PREREG P1 defines it: sum(n_both) / sum(exp_both_stratified),
    # z from the pooled sd.  The circular DeepScence-vs-CDKN1A pair is EXCLUDED from any
    # pooled number (BIO_PHASE3.md §4.4) and is reported on its own row.
    g = sig.groupby(["A", "B"]).agg(n_sections=("section", "nunique"),
                                    n_both=("n_both", "sum"),
                                    exp=("exp_both_stratified", "sum"),
                                    var=("sd_both_stratified", lambda v: float((v ** 2).sum())),
                                    above_chance=("ratio_stratified",
                                                  lambda v: int((v > 1).sum()))).reset_index()
    g["pooled_ratio"] = (g.n_both / g.exp).round(3)
    g["z"] = ((g.n_both - g.exp) / np.sqrt(g["var"])).round(2)
    g["circular"] = ((g.A == "deepscence_score") & (g.B == "cdkn1a_counts")) | \
                    ((g.B == "deepscence_score") & (g.A == "cdkn1a_counts"))
    g.to_csv(H.RESULTS + "/caller_agreement_pooled.csv", index=False)
    pd.set_option("display.width", 250)
    print(g.to_string(index=False))
    print("wrote", H.RESULTS + "/caller_*.csv")
