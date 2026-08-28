#!/usr/bin/env python3
"""Phase 10 — the H1 DeepScence FIVE-SEED CONSENSUS score (PI decision D-A).

**This specification was written down before the consensus was read**, as D-A requires
("a producer-level choice that must be written down before the consensus is read, not
after").  It is not tunable after the fact.

WHY.  The frozen single-seed configuration is not seed-reproducible on H1: full section,
`denoise=False`, SPLN21, seed 0 vs 1, Pearson r = 0.3719, top-5 % Jaccard = 0.2107
(`results/phase9_h1/d2_stability.csv`), against an M1 floor of r = 0.99553 /
Jaccard = 0.7606 and an M1 same-seed determinism control of r = 0.99999913.
**M1 keeps its frozen single-seed score, so the two arms use different DeepScence
estimators, and that asymmetry is itself reportable.**

SEEDS.  The five ALREADY-FROZEN composition-matched seeds of PREREG §3.8 --
20260901, 20260902, 20260903, 20260904, 20260905.  No new seed value is introduced.
Configuration is otherwise the frozen one (PREREG §3.9): `denoise=False`, published
`CDKN1A` anchor, >= 20 counts/cell, NATIVE human panel, full sections, no subsampling.

THE RULE, in order.

 0. **SIGN ALIGNMENT FIRST.**  DeepScence fixes its bottleneck sign by correlating with
    `CDKN1A`, and that polarity is known to flip -- it inverted between the two surgical
    arms on M1 (-0.350 sham / +0.318 SBR).  **Z-scoring does not fix a sign flip**, and a
    per-cell median over a mixture of two polarities is meaningless.  So each seed's own
    polarity is measured against the published anchor -- the depth-partialled Spearman of
    its score with `CDKN1A` counts, `deepscence_reanchor.partial_spearman`, imported so the
    definition cannot drift -- and a seed whose anchor correlation is NEGATIVE while the
    majority are positive is multiplied by -1.  **The number of seeds flipped is reported,
    per section, whatever it is.**  The raw pairwise Pearson matrix is written beside it as
    an independent check.  Phase 9 found the anchor stable in 20/20 folds in all seven
    sections, so zero flips is the expectation -- it is verified, not assumed.
 1. z-score each aligned seed score WITHIN ITS SECTION: (x - mean) / sd.
 2. consensus = the PER-CELL MEDIAN of the five z-scored seeds.  Median, not mean, for
    robustness to one divergent seed -- which is not hypothetical: on M1 one of three
    `denoise=True` seeds gave a top-5 % set PERFECTLY DISJOINT from the other two
    (Jaccard 0.000, `results/phase8_d2/d2_stability.csv`).
 3. threshold as normal: the frozen within-cell-type percentile rule
    (`h1_callers.pct_flags`, fine labels, >= 20 cells, strict >) at p90 / p95 / p99.

DISPERSION -- BOTH are reported, because they diverge.  On M1 the same comparison gave
r = 0.99553 and Jaccard = 0.7606: a score that looks stable can sit on a call-set that is
not, and the call-set is what downstream analysis consumes.
  * SCORE   : the PER-CELL INTERQUARTILE RANGE across the five z-scored seeds, summarised
              per section by its own median and quartiles.
  * CALL-SET: the MEAN PAIRWISE top-5 % JACCARD across the five seeds, WITH ITS RANGE.

PARTIAL SEED SETS.  If fewer than five seeds have completed for a section, the section is
pooled over the seeds that exist, the count `n_seeds` is written into every output row, and
the estimator is named `consensus_k<n>`.  A section with fewer than 3 seeds is skipped and
listed.  Nothing is silently averaged over an unstated number of runs.

Usage: python3 code/h1_deepscence_consensus.py [--min-seeds 3] [SPLN07 ...]
Writes data/processed_h1/deepscence_consensus_h1_<sec>.csv and
       results/phase10_h1/deepscence_consensus_{sign,dispersion,jaccard,coverage}.csv
"""
import argparse, os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H
from deepscence_reanchor import partial_spearman          # definition cannot drift
from h1_cache_extend import pct_flags                      # = h1_callers.pct_flags

SEEDS = (20260901, 20260902, 20260903, 20260904, 20260905)
RES10 = "/workspace/results/phase10_h1"
PROC = H.PROC + "/"


def seed_path(section, seed):
    return PROC + "deepscence_h1_nodn_seed%d_%s.csv" % (seed, section)


def available(section):
    return [s for s in SEEDS if os.path.exists(seed_path(section, s))
            and os.path.getsize(seed_path(section, s)) > 0]


def run(section, min_seeds=3):
    have = available(section)
    if len(have) < min_seeds:
        return None, dict(section=section, n_seeds=len(have),
                          seeds="|".join(str(s) for s in have), status="SKIPPED")
    sen = pd.read_csv(PROC + "senders_h1_%s.csv" % section).set_index("cell_id")
    sen.index = sen.index.astype(str)
    ct = pd.read_csv(PROC + "celltypes_h1_%s.csv" % section).set_index("cell_id")
    ct.index = ct.index.astype(str)

    cols, counts = {}, None
    for s in have:
        d = pd.read_csv(seed_path(section, s))
        d["cell_id"] = d.cell_id.astype(str)
        d = d.set_index("cell_id")
        cols[s] = d["deepscence_score"]
        counts = d["counts"] if counts is None else counts
    X = pd.DataFrame(cols)
    idx = X.index
    cdk = sen["cdkn1a_counts"].reindex(idx).to_numpy(float)
    dep = counts.reindex(idx).to_numpy(float)
    types = ct["cell_type"].reindex(idx).astype(str)

    # --- 0. sign alignment against the PUBLISHED anchor -----------------------
    sign_rows, flip = [], {}
    rho = {s: partial_spearman(X[s].to_numpy(float), cdk, dep) for s in have}
    n_pos = sum(1 for s in have if rho[s] > 0)
    majority_positive = n_pos * 2 >= len(have)
    for s in have:
        inverted = (rho[s] < 0) if majority_positive else (rho[s] > 0)
        flip[s] = -1.0 if inverted else 1.0
        sign_rows.append(dict(section=section, seed=s, n_cells=len(idx),
                              rho_partial_cdkn1a=round(rho[s], 5),
                              majority_sign=("+" if majority_positive else "-"),
                              flipped=bool(inverted)))
    n_flipped = int(sum(1 for s in have if flip[s] < 0))

    # independent check: the raw pairwise Pearson matrix, BEFORE any flip; and the same
    # matrix AFTER alignment, which is the score-side dispersion in the units Figure 6c
    # plots against the call-set Jaccard.
    pear = []
    for a, b in itertools.combinations(have, 2):
        r = float(np.corrcoef(X[a].to_numpy(float), X[b].to_numpy(float))[0, 1])
        ra = float(np.corrcoef(X[a].to_numpy(float) * flip[a],
                               X[b].to_numpy(float) * flip[b])[0, 1])
        pear.append(dict(section=section, seed_a=a, seed_b=b,
                         pearson_r_raw=round(r, 5), pearson_r_aligned=round(ra, 5)))

    # --- 1. z-score within section, per aligned seed --------------------------
    Z = pd.DataFrame({s: (lambda v: (v - v.mean()) / (v.std() if v.std() > 1e-12 else 1.0))(
        X[s].to_numpy(float) * flip[s]) for s in have}, index=idx)

    # --- 2. per-cell median = the consensus; per-cell IQR = its dispersion -----
    A = Z.to_numpy(float)
    cons = np.median(A, axis=1)
    iqr = np.quantile(A, 0.75, axis=1) - np.quantile(A, 0.25, axis=1)

    # --- 3. threshold as normal -----------------------------------------------
    out = pd.DataFrame({"cell_id": idx, "deepscence_consensus": np.round(cons, 5),
                        "deepscence_seed_iqr": np.round(iqr, 5)})
    for q in (90, 95, 99):
        out["consensus_flag_p%d" % q] = pct_flags(cons, types, q).astype(int)
    # The per-seed z-scored columns are deliberately NOT written: each is exactly
    # (x * flip - mean) / sd of a per-seed score file already on disk, and at 7 sections x
    # 5 seeds they cost ~60 MB of a workspace quota this box exhausted (see §12).
    # `results/phase10_h1/deepscence_consensus_sign.csv` carries every flip needed to
    # reconstruct them.
    out.to_csv(PROC + "deepscence_consensus_h1_%s.csv" % section, index=False)

    # --- call-set dispersion: mean pairwise top-5 % Jaccard, with its range ----
    f5 = {s: pct_flags(Z[s].to_numpy(float), types, 95).astype(bool) for s in have}
    jac = []
    for a, b in itertools.combinations(have, 2):
        u = int((f5[a] | f5[b]).sum()); i = int((f5[a] & f5[b]).sum())
        jac.append(dict(section=section, seed_a=a, seed_b=b, n_a=int(f5[a].sum()),
                        n_b=int(f5[b].sum()), n_both=i,
                        jaccard_top5=round(i / u, 5) if u else np.nan))
    jv = np.array([r["jaccard_top5"] for r in jac], float)

    disp = dict(section=section, n_cells=len(idx), n_seeds=len(have),
                seeds="|".join(str(s) for s in have),
                estimator="consensus_k%d" % len(have),
                n_seeds_sign_flipped=n_flipped,
                score_iqr_median=round(float(np.median(iqr)), 5),
                score_iqr_q25=round(float(np.quantile(iqr, .25)), 5),
                score_iqr_q75=round(float(np.quantile(iqr, .75)), 5),
                score_iqr_p90=round(float(np.quantile(iqr, .90)), 5),
                jaccard_top5_mean=round(float(jv.mean()), 5),
                jaccard_top5_min=round(float(jv.min()), 5),
                jaccard_top5_max=round(float(jv.max()), 5),
                pearson_raw_min=round(min(r["pearson_r_raw"] for r in pear), 5),
                pearson_raw_max=round(max(r["pearson_r_raw"] for r in pear), 5),
                pearson_aligned_mean=round(float(np.mean([r["pearson_r_aligned"]
                                                          for r in pear])), 5),
                pearson_aligned_min=round(min(r["pearson_r_aligned"] for r in pear), 5),
                pearson_aligned_max=round(max(r["pearson_r_aligned"] for r in pear), 5),
                n_senders_p95_consensus=int(out["consensus_flag_p95"].sum()),
                status="OK")
    print("%s k=%d flips=%d  score IQR med %.3f  top5 Jaccard mean %.3f [%.3f, %.3f]"
          % (section, len(have), n_flipped, disp["score_iqr_median"],
             disp["jaccard_top5_mean"], disp["jaccard_top5_min"],
             disp["jaccard_top5_max"]), flush=True)
    return dict(sign=sign_rows, pearson=pear, jaccard=jac), disp


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-seeds", type=int, default=3)
    ap.add_argument("sections", nargs="*", default=None)
    a = ap.parse_args()
    secs = a.sections or list(H.ALL_SECTIONS)
    os.makedirs(RES10, exist_ok=True)
    sign, pear, jac, disp = [], [], [], []
    for s in secs:
        parts, d = run(s, a.min_seeds)
        disp.append(d)
        if parts:
            sign += parts["sign"]; pear += parts["pearson"]; jac += parts["jaccard"]
    pd.DataFrame(disp).to_csv(RES10 + "/deepscence_consensus_coverage.csv", index=False)
    if sign:
        pd.DataFrame(sign).to_csv(RES10 + "/deepscence_consensus_sign.csv", index=False)
        pd.DataFrame(pear).to_csv(RES10 + "/deepscence_consensus_pearson.csv", index=False)
        pd.DataFrame(jac).to_csv(RES10 + "/deepscence_consensus_jaccard.csv", index=False)
    print(pd.DataFrame(disp)[["section", "n_seeds", "estimator",
                              "n_seeds_sign_flipped", "score_iqr_median",
                              "jaccard_top5_mean", "status"]].to_string(index=False))
