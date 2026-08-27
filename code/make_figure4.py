"""Figure 4 — existing spatial CCC tools under the same coordinate nulls.

Master Plan Section 25, Figure 4.  Greyscale-safe: the five conditions are
encoded by MONOTONE LUMINANCE plus distinct hatches, so the figure survives a
black-and-white print.  Reimplementations carry an asterisk everywhere.
"""
from __future__ import annotations
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sasp_palette as PAL

PAL.apply_style(matplotlib)
P4 = "/workspace/results/phase4"
FIG = "/workspace/figures"

CONDS = ["real", "N3_lig", "N4_lig", "N3_type", "N0_perm", "N0_type"]
CLAB = {"real": "real coordinates",
        "N3_lig": "N3  torus shift, ligand$^+$ cells",
        "N4_lig": "N4  rotation, ligand$^+$ cells",
        "N3_type": "N3t torus shift, per cell type",
        "N0_perm": "N0  full coordinate permutation",
        "N0_type": "N0t within-cell-type permutation\n     (CellWHISPER's own control)"}
CCOL = {"real": PAL.INK, "N3_lig": "#256abf", "N4_lig": "#5598e7",
        "N3_type": "#9ec5f4", "N0_perm": "#cde2fb", "N0_type": "#f0a860"}
CHATCH = {"real": "", "N3_lig": "///", "N4_lig": "\\\\\\",
          "N3_type": "xxx", "N0_perm": "...", "N0_type": "++"}
NULLS = CONDS[1:]
METHODS = ["COMMOT", "CellChat v2*", "SpaTalk*", "NCEM linear*"]
MSUB = {"COMMOT": "published software (v0.0.3)",
        "CellChat v2*": "reimplementation of the statistic",
        "SpaTalk*": "reimplementation of the statistic",
        "NCEM linear*": "reimplementation of the linear variant"}
PAIRS = ["Ccl2->Ccr2", "Tnf->Tnfrsf1a/1b", "Tgfb1->Tgfbr1/2", "Il1a->Il1r1"]
PLAB = {"Ccl2->Ccr2": "Ccl2→\nCcr2", "Tnf->Tnfrsf1a/1b": "Tnf→\nTnfrsf1a/b",
        "Tgfb1->Tgfbr1/2": "Tgfb1→\nTgfbr1/2", "Il1a->Il1r1": "Il1a→\nIl1r1"}

P3 = "/workspace/results/phase3"


def load_ours():
    """Our own estimator's numbers, READ FROM results/phase3 — never typed in.

    These were hard-coded literals until 2026-08-27 on the pre-C6 160-fit
    basis (N3 1.0001 / N4 0.9641), which silently regressed the figure on
    every re-run.  The reportable population is now 153, not 160
    (CS_PHASE8_TORUS_VAR.md).  Definition of "reportable" is
    summarize_phase3.py:85 verbatim: in-band sections, the primary call,
    stratum=all, naive beta > 0 AND its block-bootstrap CI excluding zero.
    """
    import sys
    sys.path.insert(0, "/workspace/code")
    import sasp_phase3 as P
    import run_phase3_nulls as RN

    mf = pd.read_csv(f"{P3}/main_fits.csv")
    pn = pd.read_csv(f"{P3}/perm_nulls.csv")
    key = ["section", "celltype", "module"]
    d = mf[mf.section.isin(P.IN_BAND) & (mf.call == RN.PRIMARY_CALL)
           & (mf.stratum == "all")].copy()
    p = pn[pn.section.isin(P.IN_BAND) & (pn.call == RN.PRIMARY_CALL)]
    cols = [c for c in p.columns if c.endswith("_sf") or c.endswith("_p")]
    cols += ["beta_obs", "N3_null_mean"]
    d = d.merge(p[key + cols].drop_duplicates(key), on=key, how="left")
    rep = d[(d.beta_naive > 0) & (d.beta_base_lo > 0)]
    if len(rep) == 0:
        raise SystemExit("figure4: no reportable fits — refusing to draw "
                         "panel c's own-estimator bars from nothing")

    # the variance-corrected (N3-var/N4-var) run is the FROZEN PRIMARY null.
    # sf_summary_var.csv is written by summarize_phase3_var.py; perm_nulls_var
    # is its input.  Absent -> we say so rather than silently dropping it.
    var = {}
    fv = f"{P3}/sf_summary_var.csv"
    if os.path.exists(fv):
        sv = pd.read_csv(fv)
        sv = sv[sv.subset.str.startswith("PRIMARY")].set_index("variant")
        for k, v in (("N3_var_sf", "N3_var"), ("N4_var_sf", "N4_var")):
            if v in sv.index:
                var[k] = float(sv.loc[v, "median"])
                var[k[:2] + "_var_n"] = int(sv.loc[v, "n"])
    else:
        print("figure4 WARNING: %s absent — the PRIMARY (variance-corrected) "
              "null cannot be drawn" % fv)

    o = dict(n_reportable=int(len(rep)),
             N3_sf=float(rep.N3_sf.median()),
             N4_sf=float(rep.N4_sf.median()),
             N1_sf=float(rep.N1_sf.median()),
             N3_reject=float((rep.N3_p < 0.05).mean()),
             N4_reject=float((rep.N4_p < 0.05).mean()),
             beta_obs=float(rep.beta_obs.median()),
             N3_null_abs=float(rep.N3_null_mean.abs().median()))
    o.update(var)
    return o


OURS = load_ours()


def load():
    R = pd.read_csv(f"{P4}/interactions.csv.gz")
    H = pd.read_csv(f"{P4}/headline.csv")
    return R, H


def sig_rate_table(R):
    """Fraction of sender x receiver interactions called significant, per
    method x pair x condition."""
    rows = []
    for (m, p), g in R.groupby(["method", "pair"]):
        base = g.drop_duplicates(["tile", "sender", "receiver"])
        rows.append(dict(method=m, pair=p, cond="real",
                         rate=float(base.real_sig.mean()),
                         n=len(base)))
        for null in NULLS:
            gg = g[g.null == null]
            if gg.empty:
                continue
            rows.append(dict(method=m, pair=p, cond=null,
                             rate=float(gg.null_sig_frac.mean()), n=len(gg)))
    return pd.DataFrame(rows)


def main():
    R, H = load()
    S = sig_rate_table(R)
    HA = H[H.pair == "ALL"].set_index(["method", "null"])
    # the CSV behind every bar in the figure, panels a, b and c
    rows = []
    for _, r in S.iterrows():
        rows.append(dict(panel="a", method=r.method, pair=r["pair"], cond=r.cond,
                         quantity="fraction of cell-type pairs called significant",
                         value=r.rate, n_interactions=r.n))
    for (m, null), r in HA.iterrows():
        pass
    for m in METHODS:
        for null in NULLS:
            if (m, null) not in HA.index:
                continue
            r = HA.loc[(m, null)]
            rows.append(dict(panel="b", method=m, pair="ALL", cond=null,
                             quantity="fraction of real-significant interactions "
                                      "still significant after shuffling",
                             value=float(r.sig_survival),
                             n_interactions=int(r.n_real_sig)))
            rows.append(dict(panel="c", method=m, pair="ALL", cond=null,
                             quantity="surviving fraction of the score "
                                      "(median null / real)",
                             value=float(r.score_sf_median),
                             n_interactions=int(r.n_real_sig)))
    for null, v in (("N3_lig", OURS["N3_sf"]), ("N4_lig", OURS["N4_sf"])):
        rows.append(dict(panel="c", method="our SASP kernel estimator (CS_PHASE3)",
                         pair="ALL", cond=null,
                         quantity="surviving fraction of the score "
                                  "(median null / real)", value=v,
                         n_interactions=OURS["n_reportable"]))
    pd.DataFrame(rows).to_csv(f"{FIG}/figure4_data.csv", index=False)

    fig = plt.figure(figsize=(12.2, 8.8))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.05], hspace=0.60,
                          wspace=0.32, left=0.068, right=0.985,
                          top=0.825, bottom=0.205)

    # ---- row 1: significance rate, real vs each null, per method ----------
    for k, m in enumerate(METHODS):
        ax = fig.add_subplot(gs[0, k])
        sub = S[S.method == m]
        w = 0.16
        for ci, c in enumerate(CONDS):
            v = [float(sub[(sub.pair == p) & (sub.cond == c)].rate.mean())
                 if len(sub[(sub.pair == p) & (sub.cond == c)]) else np.nan
                 for p in PAIRS]
            ax.bar(np.arange(4) + (ci - (len(CONDS) - 1) / 2) * w, v, w, color=CCOL[c],
                   hatch=CHATCH[c], edgecolor=PAL.INK2, linewidth=0.5)
        ax.set_xticks(range(4))
        ax.set_xticklabels([PLAB[p] for p in PAIRS], fontsize=6.2,
                           linespacing=1.3)
        ax.set_ylim(0, max(0.05, S[S.method == m].rate.max() * 1.18))
        ax.set_title(f"{m}\n{MSUB[m]}", fontsize=8.5, linespacing=1.5, pad=7)
        if k == 0:
            ax.set_ylabel("fraction of cell-type pairs\ncalled significant",
                          fontsize=8)
        ax.tick_params(axis="x", length=0)
    fig.text(0.068, 0.888, "a   Interactions called significant, real vs shuffled "
             "coordinates", fontsize=9.5, fontweight="bold")

    # ---- row 2 left: significance survival (the headline) -----------------
    ax = fig.add_subplot(gs[1, :2])
    w = 0.19
    for ni, null in enumerate(NULLS):
        v = [HA.loc[(m, null), "sig_survival"] if (m, null) in HA.index else np.nan
             for m in METHODS]
        ax.bar(np.arange(len(METHODS)) + (ni - (len(NULLS) - 1) / 2) * w, v, w, color=CCOL[null],
               hatch=CHATCH[null], edgecolor=PAL.INK2, linewidth=0.5)
    ax.axhline(0.90, color=PAL.STATUS["critical"], lw=1.4, ls="--", zorder=5)
    ax.text(-0.45, 1.32, "CellWHISPER's reported >90% FPR, under the control "
            "these bars now use (N0t, orange)", fontsize=7.5, ha="left",
            va="center", color=PAL.STATUS["critical"])
    ax.text(-0.45, 1.20, "Our bars are CALL SURVIVAL, not an FPR: survival "
            "equals an FPR only if the\nshuffle is a true null \u2014 which is "
            "the assumption under test (\u00a75).", fontsize=6.4, ha="left",
            va="center", linespacing=1.35, color=PAL.INK2)
    for ni, null in enumerate(NULLS):
        for mi, m in enumerate(METHODS):
            v = HA.loc[(m, null), "sig_survival"] if (m, null) in HA.index else np.nan
            if np.isfinite(v):
                ax.text(mi + (ni - (len(NULLS) - 1) / 2) * w, v + 0.015, f"{v:.2f}", ha="center",
                        va="bottom", fontsize=6.4, color=PAL.INK2, rotation=90)
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(METHODS, fontsize=8.5)
    ax.set_ylim(0, 1.42)
    ax.set_ylabel("fraction of real-significant interactions\n"
                  "still significant after shuffling", fontsize=8.5)
    ax.set_title("b   Does the coordinate shuffle change what is CALLED?",
                 fontsize=9.5, loc="left", pad=16, fontweight="bold")

    # ---- row 2 right: score surviving fraction + our own estimator --------
    ax2 = fig.add_subplot(gs[1, 2:])
    for ni, null in enumerate(NULLS):
        v = [HA.loc[(m, null), "score_sf_median"] if (m, null) in HA.index else np.nan
             for m in METHODS]
        ax2.bar(np.arange(len(METHODS)) + (ni - (len(NULLS) - 1) / 2) * w, v, w, color=CCOL[null],
                hatch=CHATCH[null], edgecolor=PAL.INK2, linewidth=0.5)
    xo = len(METHODS) + 0.35
    ax2.bar([xo - 0.5 * w, xo + 0.5 * w], [OURS["N3_sf"], OURS["N4_sf"]], w,
            color=[CCOL["N3_lig"], CCOL["N4_lig"]],
            hatch=[CHATCH["N3_lig"], CHATCH["N4_lig"]],
            edgecolor=PAL.INK2, linewidth=0.5)
    for ni, null in enumerate(NULLS):
        for mi, m in enumerate(METHODS):
            v = HA.loc[(m, null), "score_sf_median"] if (m, null) in HA.index else np.nan
            if np.isfinite(v):
                ax2.text(mi + (ni - (len(NULLS) - 1) / 2) * w, max(v, 0) + 0.02, f"{v:.2f}",
                         ha="center", va="bottom", fontsize=6.4,
                         color=PAL.INK2, rotation=90)
    ax2.axhline(1.0, color=PAL.MUTED, lw=1.0, ls=":")
    for xx, vv in ((xo - 0.5 * w, OURS["N3_sf"]), (xo + 0.5 * w, OURS["N4_sf"])):
        ax2.text(xx, vv + 0.02, f"{vv:.2f}", ha="center", va="bottom",
                 fontsize=6.4, color=PAL.INK2, rotation=90)
    ax2.axvline(len(METHODS) - 0.42, color=PAL.AXIS, lw=1.0)
    ax2.set_xticks(list(range(len(METHODS))) + [xo])
    ax2.set_xticklabels(METHODS + ["our SASP\nkernel estimator"], fontsize=8.5)
    ax2.set_ylabel("surviving fraction of the score\n(median null ÷ real)",
                   fontsize=8.5)
    ax2.set_title("c   Does the coordinate shuffle change the SCORE?",
                  fontsize=9.5, loc="left", pad=16, fontweight="bold")
    ymax = max(1.25, np.nanmax([HA["score_sf_median"].max(), 1.05]) * 1.15)
    ax2.set_ylim(min(0, np.nanmin(HA["score_sf_median"].min()) * 1.15), ymax)

    handles = [Patch(facecolor=CCOL[c], hatch=CHATCH[c], edgecolor=PAL.INK2,
                     linewidth=0.5, label=CLAB[c]) for c in CONDS]
    fig.legend(handles=handles, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, 0.085), fontsize=8, handlelength=2.2)

    fig.suptitle("Existing spatial cell–cell-communication tools do not "
                 "distinguish real tissue from coordinate-shuffled tissue",
                 fontsize=12.5, y=0.975)
    fig.text(0.5, 0.925,
             "6 admissible sections × 3 tiles × 4 senescence-relevant "
             "ligand–receptor pairs; sender × receiver cell-type interactions; "
             "significance at α = 0.05 (NCEM: BH-FDR q < 0.05)",
             fontsize=8.5, ha="center", color=PAL.INK2)
    fig.text(0.5, 0.022,
             "* reimplementation of the published statistic, NOT the published software. "
             "CellChat v2 and SpaTalk are R packages and NCEM needs Python ≤ 3.10 + "
             "TensorFlow; only COMMOT ran as released software.\n"
             "COMMOT is run on 1.2 mm tiles at native density because commot 0.0.3 "
             "materialises a dense n×n distance matrix (105 GB for a whole section).\n"
             "CellChat is run with type=\"mean\": at its default triMean the ligand term "
             "is exactly zero for 3 of the 4 pairs on this panel, and it calls nothing at all.",
             fontsize=7.2, ha="center", color=PAL.INK2, linespacing=1.6)

    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/figure4.png", dpi=200)
    fig.savefig(f"{FIG}/figure4.pdf")
    print("wrote", f"{FIG}/figure4.png", f"{FIG}/figure4.pdf",
          f"{FIG}/figure4_data.csv")


if __name__ == "__main__":
    main()
