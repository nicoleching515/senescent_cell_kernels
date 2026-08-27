#!/usr/bin/env python3
"""Phase 9 — verdict for audit test A7 on H1.

`code/summarize_a7.py` transplanted: same three readouts (amplitude in response-SD units,
direction with a SECTION-CLUSTERED CI, and power as the median bootstrap CI half-width), the
same designs, and the same biological-module reference — except that the reference comes from
`results/phase9_h1/h1_module_fits.csv` (produced by `code/h1_module_fits.py`, the frozen
stage-2 driver on the H1 cache) instead of the mouse `main_fits.csv`.

PREREG §6 R5: A7 on H1 is only adequately powered POOLED across sections.  The clustered
mean over 7 sections is therefore the reported quantity, and the per-fit half-width is
reported beside it so "flat" is falsifiable.

Reads  results/phase9_h1/a7_control_probe_{fits,provenance}.csv, h1_module_fits.csv
Writes results/phase9_h1/a7_summary.csv, a7_verdict.txt
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import t as tdist
sys.path.insert(0, "/workspace/code")
import h1_common as H

RES = H.RESULTS + "/"
DESIGNS = [("base", "naive (intercept only)"),
           ("n6", "+N6 neighbour baseline"),
           ("n5", "+N5 technical covariates"),
           ("n6n5", "+N6+N5 (full nuisance design)"),
           ("n2", "N2 matched-decoy contrast")]
CALLS = ["tierA_p95", "cdkn1a_pos"]


SD_FLOOR = 1e-6   # a fit whose response has no variance at all inside its receiver set


def add(df):
    # beta/sd_y is undefined where the response is constant.  On H1 the 609 negative control
    # CODEWORDS are 62x sparser than on M1 (0.0007 vs 0.0428 counts/cell) and 10 of 98
    # codeword fits have sd_y == 0 exactly, which makes their amplitude ratio explode.  Those
    # fits are dropped and counted rather than left to dominate a mean.
    df = df.copy()
    df["degenerate_sd_y"] = df.sd_y < SD_FLOOR
    for k, _ in DESIGNS:
        df["e_" + k] = df["beta_" + k] / df.sd_y
        df["sig_" + k] = (df["beta_%s_lo" % k] * df["beta_%s_hi" % k]) > 0
        df["hw_" + k] = (df["beta_%s_hi" % k] - df["beta_%s_lo" % k]) / 2 / df.sd_y
    return df


def clustered(df, col, cluster="section"):
    m = df.groupby(cluster)[col].mean().to_numpy()
    m = m[np.isfinite(m)]
    if m.size < 3:
        return np.nan, np.nan, np.nan, np.nan
    mu = m.mean(); se = m.std(ddof=1) / np.sqrt(m.size)
    tc = tdist.ppf(.975, m.size - 1)
    return mu, mu - tc * se, mu + tc * se, 2 * tdist.sf(abs(mu / se), m.size - 1)


def main():
    A = pd.read_csv(RES + "a7_control_probe_fits.csv")
    A = A[A.get("skip").isna()] if "skip" in A.columns else A
    A = add(A)
    M = pd.read_csv(RES + "h1_module_fits.csv")
    M = add(M[M.call.isin(CALLS)].copy())
    M["response"] = "BIOLOGICAL MODULES (reference)"
    prov = pd.read_csv(RES + "a7_control_probe_provenance.csv")
    rows = []
    for src in (A, M):
        for resp, g0 in src.groupby("response"):
            n_drop = int(g0.degenerate_sd_y.sum())
            g = g0[~g0.degenerate_sd_y]
            for k, dl in DESIGNS:
                mu, lo, hi, p = clustered(g, "e_" + k)
                rows.append(dict(response=resp, design=k, design_label=dl,
                                 n_fits=len(g), n_sections=g.section.nunique(),
                                 median_abs_amplitude=round(g["e_" + k].abs().median(), 4),
                                 median_signed_amplitude=round(g["e_" + k].median(), 4),
                                 frac_positive=round(float((g["beta_" + k] > 0).mean()), 3),
                                 frac_CI_excludes_zero=round(float(g["sig_" + k].mean()), 3),
                                 clustered_mean=round(mu, 4), clustered_lo=round(lo, 4),
                                 clustered_hi=round(hi, 4), clustered_p=("%.3g" % p),
                                 median_CI_halfwidth=round(g["hw_" + k].median(), 4),
                                 n_dropped_degenerate=int(n_drop)))
    S = pd.DataFrame(rows)
    S.to_csv(RES + "a7_summary.csv", index=False)
    ctrl = S[S.response != "BIOLOGICAL MODULES (reference)"]
    bio = S[S.response == "BIOLOGICAL MODULES (reference)"].set_index("design")
    L = ["A7 -- negative-control kernel, HUMAN arm (H1, GSE326743), %d sections x %d sender calls"
         % (A.section.nunique(), len(CALLS)),
         "%d control fits, %d biological-module fits, same estimator "
         "(WINDOW=100um, 40-point lambda grid, MIN_RECEIVERS=2000, "
         "400-replicate spatial block bootstrap)." % (len(A), len(M)), "",
         "PRIMARY pre-registered technical null: the 40 NEGATIVE CONTROL PROBES "
         "(`neg_control_probe`).  `all_controls` is the POOLED control feature set and must "
         "never be called 'negative control probes'.", "",
         "Control-feature sparsity (why power is the binding constraint):"]
    for r, g in prov.groupby("response"):
        L.append("  %-22s mean %.4f counts/cell, %.2f%% of cells non-zero"
                 % (r, g.mean_per_cell.mean(), 100 * g.frac_cells_nonzero.mean()))
    L.append("")
    for k, dl in DESIGNS:
        c = ctrl[ctrl.design == k]; b = bio.loc[k]
        p = c[c.response == "neg_control_probe"].iloc[0]
        L.append("%-30s PRIMARY 40 probes: %+.4f [%+.4f, %+.4f] p=%s | pooled controls "
                 "|beta|/sd med %.4f | modules %.4f (CI excl 0: %.2f)"
                 % (dl, p.clustered_mean, p.clustered_lo, p.clustered_hi, p.clustered_p,
                    c.median_abs_amplitude.median(), b.median_abs_amplitude,
                    b.frac_CI_excludes_zero))
    L.append("")
    L.append("Smallest amplitude one A7 fit could resolve (median CI half-width): "
             "%.4f SD naive, %.4f SD conditioned."
             % (ctrl[ctrl.design == "base"].median_CI_halfwidth.median(),
                ctrl[ctrl.design == "n6n5"].median_CI_halfwidth.median()))
    L.append("Biological module amplitude to be ruled out: %.4f SD naive, %.4f SD conditioned."
             % (bio.loc["base"].median_abs_amplitude, bio.loc["n6n5"].median_abs_amplitude))
    pr = ctrl[(ctrl.response == "neg_control_probe")].set_index("design")
    flat = all(float(pr.loc[k].clustered_lo) <= 0 <= float(pr.loc[k].clustered_hi)
               for k in ("base", "n6n5"))
    L.append("")
    L.append("VERDICT on the pre-registered primary response (40 negative control probes): "
             "%s -- the section-clustered mean %s zero under both the naive and the full "
             "N6+N5 design." % ("FLAT" if flat else "NOT FLAT",
                                "includes" if flat else "excludes"))
    txt = "\n".join(L)
    open(RES + "a7_verdict.txt", "w").write(txt + "\n")
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
    print(txt); print(); print(S.to_string(index=False))


if __name__ == "__main__":
    main()
