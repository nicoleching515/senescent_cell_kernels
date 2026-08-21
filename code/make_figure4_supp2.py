"""Supplementary to Figure 4 — why COMMOT's answer does not move.

COMMOT's optimal transport is highly sensitive to geometry: permuting
coordinates replaces essentially the entire cell-to-cell communication network.
But the transport conserves total ligand mass, and the cluster-level summary
that COMMOT tests keeps ~0.4-0.95 rank correlation with its real-coordinate
value, so the calls barely change.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sasp_palette as PAL
PAL.apply_style(matplotlib)

PAIRS = ["Ccl2->Ccr2", "Tnf->Tnfrsf1a/1b", "Tgfb1->Tgfbr1/2", "Il1a->Il1r1"]
PLAB = {"Ccl2->Ccr2": "Ccl2→Ccr2", "Tnf->Tnfrsf1a/1b": "Tnf→Tnfrsf1a/1b",
        "Tgfb1->Tgfbr1/2": "Tgfb1→Tgfbr1/2", "Il1a->Il1r1": "Il1a→Il1r1"}
MK = ["o", "s", "^", "D"]


def main():
    d = pd.read_csv("/workspace/results/phase4/commot_mechanism.csv")
    d["mass_ratio"] = d.total_mass_perm / d.total_mass_real
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.7))

    ax = axes[0]
    for i, p in enumerate(PAIRS):
        g = d[d.pair == p]
        ax.scatter(np.full(len(g), i), g.edge_jaccard, s=34, marker=MK[i],
                   color="#256abf", edgecolor=PAL.INK2, linewidth=0.4, zorder=3)
        ax.plot([i - .3, i + .3], [g.edge_jaccard.median()] * 2, color=PAL.INK, lw=1.6)
    ax.set_xticks(range(4)); ax.set_xticklabels([PLAB[p] for p in PAIRS],
                                                fontsize=7.5, rotation=18, ha="right")
    ax.set_ylim(0, 0.06); ax.set_ylabel("Jaccard overlap of communicating\ncell pairs, real vs N0-permuted")
    ax.set_title("a   the cell-level network is destroyed", loc="left",
                 fontsize=9.5, fontweight="bold")

    ax = axes[1]
    for i, p in enumerate(PAIRS):
        g = d[d.pair == p]
        ax.scatter(np.full(len(g), i), g.mass_ratio, s=34, marker=MK[i],
                   color="#9ec5f4", edgecolor=PAL.INK2, linewidth=0.4, zorder=3)
    ax.axhline(1.0, color=PAL.INK, lw=1.2, ls=":")
    ax.set_xticks(range(4)); ax.set_xticklabels([PLAB[p] for p in PAIRS],
                                                fontsize=7.5, rotation=18, ha="right")
    ax.set_ylim(0.9, 1.25)
    ax.set_ylabel("total transported ligand mass,\nN0-permuted ÷ real")
    ax.set_title("b   but the transported mass is conserved", loc="left",
                 fontsize=9.5, fontweight="bold")

    ax = axes[2]
    for i, p in enumerate(PAIRS):
        g = d[d.pair == p]
        ax.scatter(np.full(len(g), i), g.cluster_spearman, s=34, marker=MK[i],
                   color=PAL.INK, edgecolor="w", linewidth=0.4, zorder=3)
        ax.plot([i - .3, i + .3], [g.cluster_spearman.median()] * 2,
                color=PAL.STATUS["critical"], lw=1.8)
    ax.set_xticks(range(4)); ax.set_xticklabels([PLAB[p] for p in PAIRS],
                                                fontsize=7.5, rotation=18, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Spearman ρ of the cluster-level\nscore, real vs N0-permuted")
    ax.set_title("c   so the cluster summary barely moves", loc="left",
                 fontsize=9.5, fontweight="bold")

    fig.suptitle("COMMOT — the optimal transport reads the geometry; the "
                 "cluster-level summary that COMMOT tests throws it away",
                 fontsize=11, y=1.02)
    axes[1].text(0.06, 1.19, "Ccr2 is detected in ~1–3% of cells, so for Ccl2 not\n"
                 "all ligand mass finds a receiver within 100 µm and the\n"
                 "conservation is only approximate for that pair.",
                 fontsize=6.6, color=PAL.INK2, va="top", linespacing=1.4)
    fig.text(0.5, -0.13, "One point per section (first tile of each of the six "
             "admissible sections).  N0 = every cell coordinate permuted, i.e. no "
             "spatial information left in the data.\nCOMMOT 0.0.3, published "
             "software.", ha="center", fontsize=7.5, color=PAL.INK2,
             linespacing=1.5)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"/workspace/figures/figure4_supp_commot_mechanism.{e}",
                    dpi=200, bbox_inches="tight")
    print(d.groupby("pair")[["edge_jaccard", "mass_ratio", "cluster_spearman"]]
          .median().round(4).to_string())
    print("\nall pairs: jaccard med %.4f  mass ratio med %.6f  cluster rho med %.3f"
          % (d.edge_jaccard.median(), d.mass_ratio.median(), d.cluster_spearman.median()))
    print("significant real / N0 / shared: %d / %d / %d"
          % (d.nsig_real.sum(), d.nsig_perm.sum(), d.nsig_shared.sum()))


if __name__ == "__main__":
    main()
