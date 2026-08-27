#!/usr/bin/env python3
"""Two-panel summary for the §29 objection-9 Moran's I run.

(a) The Voyager reproduction: per-feature Moran's I by h5 feature class, all
    11 M1 sections pooled.  Controls sit on zero; genes do not.
(b) The comparison that answers the objection: aggregate Moran's I on exactly
    the A7 responses, against the A7 naive kernel amplitude on the same
    features and cells.

Writes results/moran/moran_summary.png / .pdf.  Nothing under figures/ is
touched, so `check_figures_guard.py` is unaffected.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MOR = "/workspace/results/moran/"
ORDER = ["Gene Expression", "Deprecated Codeword", "Genomic Control",
         "Negative Control Probe", "Negative Control Codeword",
         "Unassigned Codeword"]
COL = {"Gene Expression": "#1b6ca8", "Deprecated Codeword": "#9e9e9e",
       "Genomic Control": "#d1495b", "Negative Control Probe": "#e07a5f",
       "Negative Control Codeword": "#f2b134",
       "Unassigned Codeword": "#7a8b99"}
# label offsets in panel (b); the three near-zero control points collide
OFF = {"all_controls": (7, 5), "neg_control_codeword": (7, -11),
       "neg_control_probe": (7, 6), "genomic_control": (-46, -12),
       "neg_probe_rate": (7, -4), "downstream_arrest": (-8, -13),
       "secondary_senescence": (6, 6), "oxidative_stress": (7, -3)}


def main():
    F = pd.read_csv(MOR + "moran_per_feature.csv.gz")
    P = pd.read_csv(MOR + "moran_pooled.csv")
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0))

    a = ax[0]
    rng = np.random.default_rng(0)
    for i, ft in enumerate(ORDER):
        v = F.loc[F.feature_type == ft, "moran_I"].dropna().to_numpy()
        if not len(v):
            continue
        y = i + (rng.random(len(v)) - .5) * .6
        a.scatter(v, y, s=2, alpha=.20, color=COL[ft], rasterized=True,
                  linewidths=0)
        a.plot([np.median(v)] * 2, [i - .38, i + .38], color="k", lw=2, zorder=5)
        a.text(0.62, i + .30, "n=%d  med %+.5f" % (len(v), np.median(v)),
               fontsize=7.5, ha="right", va="center")
    a.axvline(0, color="k", lw=.8, ls=":")
    a.set_yticks(range(len(ORDER)))
    a.set_yticklabels([o.replace("Negative Control ", "Neg. control\n")
                       for o in ORDER], fontsize=8)
    a.set_xlabel("Moran's I, per feature (k = 6 NN, row-standardised)")
    a.set_xlim(-0.03, 0.68)
    a.set_title("(a) Voyager's test, reproduced on M1\n"
                "one point = one feature in one of 11 sections", fontsize=10)

    b = ax[1]
    sub = P[P.kind.isin(["control", "module"])]
    for _, r in sub.iterrows():
        c = "#d1495b" if r.kind == "control" else "#1b6ca8"
        b.errorbar(r.a7_base_mean, r.I_raw_mean,
                   xerr=[[r.a7_base_mean - r.a7_base_lo],
                         [r.a7_base_hi - r.a7_base_mean]],
                   yerr=[[r.I_raw_mean - r.I_raw_lo],
                         [r.I_raw_hi - r.I_raw_mean]],
                   fmt="o", ms=6, color=c, ecolor=c, elinewidth=1, capsize=2)
        b.annotate(r.field, (r.a7_base_mean, r.I_raw_mean),
                   textcoords="offset points", xytext=OFF.get(r.field, (6, 4)),
                   fontsize=7.5, color=c)
    b.axhline(0, color="k", lw=.8, ls=":")
    b.axvline(0, color="k", lw=.8, ls=":")
    b.set_xlabel("A7 naive kernel amplitude  β/sd$_y$  (section-clustered mean)")
    b.set_ylabel("Moran's I, aggregate field (k = 6 NN)")
    b.set_title("(b) the two statistics on the same features and cells\n"
                "red = A7 control responses, blue = Tier B modules", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(MOR, "moran_summary." + ext), dpi=200)
    print("wrote", MOR + "moran_summary.png")


if __name__ == "__main__":
    main()
