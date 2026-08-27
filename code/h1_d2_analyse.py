#!/usr/bin/env python3
"""Phase 9 — PREREG §8 predictions P-vi and P-vii on H1.

P-vi  Switching `denoise=True` RAISES DeepScence's depth loading rather than lowering it, as
      it does on 3 of 3 M1 sections (Delta rho +0.13 ... +0.25).  Falsified if Delta rho <= 0.
P-vii The `denoise=True` seed instability recurs: >= 1 of 3 seeds on an H1 section gives a
      top-5 % sender set with Jaccard < 0.30 against the others.  Falsified if all three
      agree at Jaccard >= 0.60, as `denoise=False` does on M1.

Design copied from the mouse D2 run (`reports/CS_PHASE8_D2_DENOISE.md`): ONE fixed
20,000-cell subsample (subsampling seed 12345, independent of `--seed`, so every run sees the
same cells), three DeepScence seeds, nothing else changed, plus the `denoise=False` companion
on the identical cells as the seed-to-seed floor.

Usage: python3 code/h1_d2_analyse.py
Writes results/phase9_h1/d2_depth.csv, d2_stability.csv
"""
import sys, os, glob, itertools
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, "/workspace/code")
import h1_common as H


def load(tag, section):
    p = H.PROC + "/deepscence_h1_%s_%s.csv" % (tag, section)
    return pd.read_csv(p).set_index("cell_id") if os.path.exists(p) else None


def top5(v):
    return v > np.nanpercentile(v, 95)


def main():
    rows_depth, rows_stab = [], []
    # ---- P-vi: depth loading, denoise=False vs denoise=True -------------------------
    for section in H.ALL_SECTIONS:
        full_false = H.PROC + "/deepscence_h1_%s.csv" % section
        dca_full = load("dca", section)
        if dca_full is None or not os.path.exists(full_false):
            continue
        f = pd.read_csv(full_false).set_index("cell_id")
        j = f.join(dca_full, how="inner", rsuffix="_dca")
        cells = H.cells_table(section).set_index("cell_id").reindex(j.index)
        tc = cells.transcript_counts.to_numpy(float)
        r0 = spearmanr(j.deepscence_score.to_numpy(float), tc).statistic
        r1 = spearmanr(j.deepscence_score_dca.to_numpy(float), tc).statistic
        a, b = top5(j.deepscence_score.to_numpy(float)), top5(j.deepscence_score_dca.to_numpy(float))
        rows_depth.append(dict(section=section, scope="full section", n=len(j),
                               rho_denoise_False=round(float(r0), 4),
                               rho_denoise_True=round(float(r1), 4),
                               delta_rho=round(float(r1 - r0), 4),
                               ratio=round(float(r1 / r0), 3) if r0 else None,
                               sender_jaccard_False_vs_True=round(
                                   float((a & b).sum() / max((a | b).sum(), 1)), 4)))
    # the 20k panel, denoise=False seed 0 vs denoise=True seed 0, identical cells
    for section in H.ALL_SECTIONS:
        n0, d0 = load("nodn_sub20000", section), load("dca_sub20000", section)
        if n0 is None or d0 is None:
            continue
        j = n0.join(d0, how="inner", rsuffix="_dca")
        tc = j.counts.to_numpy(float)
        r0 = spearmanr(j.deepscence_score.to_numpy(float), tc).statistic
        r1 = spearmanr(j.deepscence_score_dca.to_numpy(float), tc).statistic
        a, b = top5(j.deepscence_score.to_numpy(float)), top5(j.deepscence_score_dca.to_numpy(float))
        rows_depth.append(dict(section=section, scope="20,000-cell panel", n=len(j),
                               rho_denoise_False=round(float(r0), 4),
                               rho_denoise_True=round(float(r1), 4),
                               delta_rho=round(float(r1 - r0), 4),
                               ratio=round(float(r1 / r0), 3) if r0 else None,
                               sender_jaccard_False_vs_True=round(
                                   float((a & b).sum() / max((a | b).sum(), 1)), 4)))

    # ---- P-vii: seed stability -------------------------------------------------------
    for section in H.ALL_SECTIONS:
        for cfg, label in (("dca_sub20000", "denoise=True"), ("nodn_sub20000", "denoise=False")):
            runs = {}
            for s in (0, 1, 2):
                tag = cfg if s == 0 else cfg + "_seed%d" % s
                d = load(tag, section)
                if d is not None:
                    runs[s] = d.deepscence_score
            for i, k in itertools.combinations(sorted(runs), 2):
                v1 = runs[i].to_numpy(float); v2 = runs[k].reindex(runs[i].index).to_numpy(float)
                a, b = top5(v1), top5(v2)
                rows_stab.append(dict(
                    section=section, config=label, seed_a=i, seed_b=k, n=len(v1),
                    pearson_r=round(float(pearsonr(v1, v2)[0]), 4),
                    spearman_r=round(float(spearmanr(v1, v2).statistic), 4),
                    top5_jaccard=round(float((a & b).sum() / max((a | b).sum(), 1)), 4),
                    n_cells_changing_status=int((a ^ b).sum())))
    # ---- full-section seed check on the FROZEN PRIMARY configuration ------------------
    # The 20,000-cell panel above could be a small-sample artefact, so the committed
    # random_state=0 full-section score is compared against a random_state=1 re-run of the
    # same configuration on the same cells.
    for section in H.ALL_SECTIONS:
        p0 = H.PROC + "/deepscence_h1_%s.csv" % section
        s1 = load("nodn_seed1", section)
        if s1 is None or not os.path.exists(p0):
            continue
        j = pd.read_csv(p0).set_index("cell_id").join(s1, how="inner", rsuffix="_s1")
        v1 = j.deepscence_score.to_numpy(float); v2 = j.deepscence_score_s1.to_numpy(float)
        a, b = top5(v1), top5(v2)
        rows_stab.append(dict(
            section=section, config="denoise=False, FULL section", seed_a=0, seed_b=1,
            n=len(j), pearson_r=round(float(pearsonr(v1, v2)[0]), 4),
            spearman_r=round(float(spearmanr(v1, v2).statistic), 4),
            top5_jaccard=round(float((a & b).sum() / max((a | b).sum(), 1)), 4),
            n_cells_changing_status=int((a ^ b).sum())))
        cells = H.cells_table(section).set_index("cell_id").reindex(j.index)
        tc = cells.transcript_counts.to_numpy(float)
        rows_depth.append(dict(section=section, scope="full section, seed check", n=len(j),
                               rho_denoise_False=round(float(spearmanr(v1, tc).statistic), 4),
                               rho_denoise_True=None, delta_rho=None, ratio=None,
                               sender_jaccard_False_vs_True=None,
                               rho_seed1=round(float(spearmanr(v2, tc).statistic), 4)))

    d1 = pd.DataFrame(rows_depth); d2 = pd.DataFrame(rows_stab)
    d1.to_csv(H.RESULTS + "/d2_depth.csv", index=False)
    d2.to_csv(H.RESULTS + "/d2_stability.csv", index=False)
    pd.set_option("display.width", 240)
    print("=== P-vi: depth loading, denoise=False vs denoise=True ===")
    print(d1.to_string(index=False) if len(d1) else "(no denoise=True run yet)")
    print("\nM1: 0.3891->0.6404, 0.3176->0.5314, 0.4096->0.5419 (x1.32-1.67, 3 of 3 sections)")
    print("\n=== P-vii: seed stability on one fixed 20,000-cell subsample ===")
    print(d2.to_string(index=False) if len(d2) else "(no seed panel yet)")
    print("\nM1: denoise=False seeds 0/1 r=0.9955 Jaccard=0.761; denoise=True seeds 0/1 and 1/2 "
          "r=0.57 Jaccard=0.000 (2,000 of 2,000 cells changing status)")


if __name__ == "__main__":
    main()
