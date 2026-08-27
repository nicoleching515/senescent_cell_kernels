"""Supplementary to Figure 4 — NCEM's length scale is not identified.

NCEM's linear variant selects its interaction radius as the r maximising
variance explained.  If that criterion is flat, the reported length scale is an
artefact of the argmax, and if it is nearly as high on randomised coordinates
the model is not measuring space at all.  Both are true here.
"""
import glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sasp_palette as PAL
PAL.apply_style(matplotlib)

CONDS = ["real", "N3_lig", "N4_lig", "N3_type", "N0_perm"]
CLAB = {"real": "real coordinates", "N3_lig": "N3  torus shift, ligand$^+$",
        "N4_lig": "N4  rotation, ligand$^+$", "N3_type": "N3t per-cell-type shift",
        "N0_perm": "N0  full coordinate permutation"}
CCOL = {"real": PAL.INK, "N3_lig": "#256abf", "N4_lig": "#5598e7",
        "N3_type": "#9ec5f4", "N0_perm": "#b9b7b0"}
MK = {"real": "o", "N3_lig": "s", "N4_lig": "^", "N3_type": "D", "N0_perm": "v"}
LS = {"real": "-", "N3_lig": "--", "N4_lig": "-.", "N3_type": ":", "N0_perm": (0, (3, 1, 1, 1))}


def main():
    d = pd.concat([pd.read_csv(f) for f in
                   glob.glob("/workspace/results/phase4/parts/ncem_sweep_*.csv")],
                  ignore_index=True).drop_duplicates(["tile", "cond", "radius"])
    d.to_csv("/workspace/results/phase4/ncem_radius_sweep.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    ax = axes[0]
    for c in CONDS:
        g = d[d.cond == c].groupby("radius").r2.median()
        ax.plot(g.index, g.values, marker=MK[c], ls=LS[c], color=CCOL[c],
                label=CLAB[c], lw=1.6, ms=4.5)
    ax.set_xlabel("NCEM neighbourhood radius r (µm)")
    ax.set_ylabel("median variance explained $R^2$\n(18 tiles × 4 LR pairs)")
    ax.set_title("a   the model-selection criterion is flat", loc="left",
                 fontsize=9.5, fontweight="bold")
    ax.legend(fontsize=7.2, loc="lower left", ncol=1)
    # These three numbers were typed in (0.0104, 0.0136, "77%") against a
    # 2026-08-20 vintage of `parts/`.  The parts were regenerated on
    # 2026-08-21 and 575 of the 720 r2 values moved, so the annotation drifted
    # away from the data the same script was plotting -- and this figure's own
    # output CSV was the only surviving record of the old state.  Computed now.
    r2 = d.groupby("cond").r2.median()
    r2_real, r2_rand = float(r2["real"]), float(r2["N0_perm"])
    ax.annotate("", xy=(103, r2_real), xytext=(103, r2_rand),
                arrowprops=dict(arrowstyle="<->", color=PAL.INK2, lw=1.0))
    ax.text(108, 0.5 * (r2_real + r2_rand),
            "%.1f%% of the explained\nvariance survives full\n"
            "coordinate randomisation" % (100 * r2_rand / r2_real),
            fontsize=7, va="center", ha="left",
            color=PAL.INK2, linespacing=1.5)
    ax.set_xlim(5, 165)
    ax.set_ylim(0, 0.019)

    ax = axes[1]
    best = d.loc[d.groupby(["tile", "cond"]).r2.idxmax()]
    for i, c in enumerate(CONDS):
        v = best[best.cond == c].radius.values
        j = (np.arange(len(v)) % 5 - 2) * 0.045
        ax.scatter(np.full(len(v), i) + j, v, s=26, color=CCOL[c],
                   marker=MK[c], edgecolor=PAL.INK2, linewidth=0.4, zorder=3)
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=PAL.INK, lw=1.6,
                zorder=4)
    ax.set_xticks(range(len(CONDS)))
    ax.set_xticklabels(["real", "N3", "N4", "N3t", "N0"], fontsize=8.5)
    ax.set_yscale("log")
    ax.set_yticks([10, 20, 50, 100]); ax.set_yticklabels([10, 20, 50, 100])
    ax.set_ylabel("selected radius (µm), one point per tile")
    ax.set_title("b   so the selected length scale is arbitrary", loc="left",
                 fontsize=9.5, fontweight="bold")
    fig.suptitle("NCEM linear* — the reported interaction length scale is not "
                 "identified in this tissue", fontsize=11, y=1.0)
    fig.text(0.5, -0.06, "* reimplementation of the published linear variant, "
             "not the ncem package.", ha="center", fontsize=7.5, color=PAL.INK2)
    fig.tight_layout()
    fig.savefig("/workspace/figures/figure4_supp_ncem_lengthscale.png",
                dpi=200, bbox_inches="tight")
    fig.savefig("/workspace/figures/figure4_supp_ncem_lengthscale.pdf",
                bbox_inches="tight")
    print(best.groupby("cond").radius.describe()[["min", "50%", "max"]].to_string())
    print(d.groupby("cond").r2.median().round(4).to_string())


if __name__ == "__main__":
    main()
