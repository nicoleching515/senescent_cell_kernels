#!/usr/bin/env python3
"""Figure 6 (new, Phase7 §19) — DEEPSCENCE NATIVE vs REMAPPED.  This carries §8.

(a) Caller agreement at full coverage, both arms, BEFORE and AFTER conditioning on cell
    type and transcript-depth decile.
(b) The sign anchor: the depth-partialled correlation of the score with its own published
    `CDKN1A` anchor, per section, per arm, with the fold-split sign stability.  This is
    where pre-registered prediction P-ii was FALSIFIED.
(c) Whether the instability transfers: seed-to-seed agreement of the score and of the
    top-5 % call set, on each arm.

Palette = the project's pre-validated categorical theme (`sasp_palette`).  Every plotted
number is written to figures/figure6_data.csv; .png and .pdf both emitted (PREREG §11).

READ EVERY DEEPSCENCE NUMBER WITH ITS FIVE ATTRIBUTES (PREREG §3.9 as extended by D-A):
coverage, denoise state, anchor, panel, and seed configuration.  M1 = 11/11 sections,
denoise=False, published CDKN1A anchor, ORTHOLOG-REMAPPED panel (4,845 of 5,097),
random_state=0 single run.  H1 = 7/7 sections, denoise=False, published CDKN1A anchor,
NATIVE human panel (5,093, no remap); panels (a) and (b) are at random_state=0 and are
therefore on a NON-PRIMARY estimator under D-A -- the directions are seed-robust, the
magnitudes are not, and panel (c) is the measurement of exactly that.

Usage: python3 code/make_figure6.py
"""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/workspace/code")
import sasp_palette as PAL
PAL.apply_style(matplotlib)

FIG = "/workspace/figures/"
R3 = "/workspace/results/phase3"
R8 = "/workspace/results/phase8_d2"
R9 = "/workspace/results/phase9_h1"
RH = "/workspace/results/phase10_h1"
C_M1, C_H1 = PAL.SERIES[0], PAL.SERIES[1]
ROWS = []
PAIRS = [(("tierA_score", "deepscence_score"), "Tier A × DeepScence", False),
         (("tierA_score", "cdkn1a_counts"), "Tier A × CDKN1A$^+$", False),
         (("tierA_score", "senepy_score"), "Tier A × SenePy", False),
         (("senepy_score", "deepscence_score"), "SenePy × DeepScence", False),
         (("deepscence_score", "cdkn1a_counts"),
          "DeepScence × CDKN1A$^+$  (CIRCULAR)", True)]


def _pick(d, a, b):
    m = ((d.A == a) & (d.B == b)) | ((d.A == b) & (d.B == a))
    return d[m].ratio.dropna().to_numpy(float)


def main():
    fig = plt.figure(figsize=(13.6, 12.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.95, 0.80],
                          width_ratios=[1.25, 1.0], hspace=0.46, wspace=0.26)

    # ---------------- (a) agreement before / after conditioning -----------
    ax = fig.add_subplot(gs[0, :])
    srcs = {}
    srcs[("M1", "before")] = pd.read_csv(
        f"{R3}/caller_pairwise_agreement_11sections.csv").query("threshold=='global_top5'")
    srcs[("M1", "after")] = pd.read_csv(
        f"{R3}/caller_agreement_depth_and_type_matched_11sections.csv")
    p = f"{RH}/caller_pairwise_agreement_global_h1.csv"
    if os.path.exists(p):
        srcs[("H1", "before")] = pd.read_csv(p)
    srcs[("H1", "after")] = pd.read_csv(
        f"{R9}/caller_agreement_depth_and_type_matched.csv")
    y = np.arange(len(PAIRS))[::-1]
    style = {("M1", "before"): dict(c=C_M1, off=+0.27, mfc=PAL.SURFACE, m="o"),
             ("M1", "after"):  dict(c=C_M1, off=+0.09, mfc=C_M1, m="o"),
             ("H1", "before"): dict(c=C_H1, off=-0.09, mfc=PAL.SURFACE, m="s"),
             ("H1", "after"):  dict(c=C_H1, off=-0.27, mfc=C_H1, m="s")}
    ax.axvline(1.0, color=PAL.INK2, lw=1.2, zorder=1)
    for k, st in style.items():
        d = srcs.get(k)
        if d is None:
            continue
        for yi, ((a, b), lab, circ) in zip(y, PAIRS):
            v = _pick(d, a, b)
            if not v.size:
                continue
            ax.plot([v.min(), v.max()], [yi + st["off"]] * 2, "-", color=st["c"], lw=2.0,
                    alpha=.45, zorder=2)
            ax.plot([np.median(v)], [yi + st["off"]], st["m"], color=st["c"],
                    mfc=st["mfc"], ms=6.4, mew=1.4, zorder=3)
            ROWS.append(dict(panel="a", arm=k[0], conditioning=k[1], pair=lab,
                             n_sections=int(v.size), ratio_median=round(float(np.median(v)), 4),
                             ratio_min=round(float(v.min()), 4),
                             ratio_max=round(float(v.max()), 4), circular=circ))
        ax.plot([], [], st["m"], color=st["c"], mfc=st["mfc"], ms=6.4, mew=1.4,
                label="%s, %s conditioning" % (k[0], k[1]))
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l, _ in PAIRS], fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlabel("top-5 % call overlap, as a ratio of chance  (1 = chance; "
                  "point = median over sections, bar = range)")
    ax.set_title("(a)  Caller agreement at full coverage, both arms, before and after "
                 "conditioning on cell type × depth decile.\n"
                 "No pair is at chance on either arm — the caller-independence claim is "
                 "dead in both (PREREG §10.2).  The circular pair is shown and is excluded "
                 "from every pooled claim.", loc="left")
    ax.set_ylim(-0.75, len(PAIRS) - 0.25)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=7.5)

    # ---------------- (b) the sign anchor ---------------------------------
    ax = fig.add_subplot(gs[1, 0])
    m = pd.read_csv(f"{R3}/deepscence_anchor_decisions.csv")
    h = pd.read_csv(f"{R9}/deepscence_anchor_h1.csv")
    ax.axhline(0, color=PAL.INK2, lw=1.2, zorder=1)
    xs = np.arange(len(m))
    ax.plot(xs, m.prho_ds_Cdkn1a, "o", color=C_M1, ms=7, mec=PAL.SURFACE, mew=1.3,
            label="M1, ortholog-remapped panel (11 sections)", zorder=3)
    for i, r in enumerate(m.itertuples()):
        if r.stab_Cdkn1a < 0.90:
            ax.annotate("stab %.2f" % r.stab_Cdkn1a, (i, r.prho_ds_Cdkn1a),
                        textcoords="offset points", xytext=(0, 9), ha="center",
                        fontsize=6.6, color=PAL.STATUS["critical"])
        ROWS.append(dict(panel="b", arm="M1", section=r.section,
                         rho_partial_cdkn1a=float(r.prho_ds_Cdkn1a),
                         fold_split_sign_stability=float(r.stab_Cdkn1a)))
    xs2 = np.arange(len(m), len(m) + len(h))
    ax.plot(xs2, h.rho_partial_cdkn1a, "s", color=C_H1, ms=7, mec=PAL.SURFACE, mew=1.3,
            label="H1, NATIVE human panel (7 sections)", zorder=3)
    for i, r in zip(xs2, h.itertuples()):
        ROWS.append(dict(panel="b", arm="H1", section=r.section,
                         rho_partial_cdkn1a=float(r.rho_partial_cdkn1a),
                         fold_split_sign_stability=float(r.stab_cdkn1a)))
    ax.axvline(len(m) - 0.5, color=PAL.AXIS, lw=1)
    ax.set_xticks(list(xs) + list(xs2))
    ax.set_xticklabels([s.split("_")[0] for s in m.section] + list(h.section),
                       rotation=90, fontsize=6.5)
    ax.set_ylabel(r"depth-partialled $\rho$(score, $CDKN1A$)")
    ax.set_title("(b)  P-ii, FALSIFIED.  The published anchor is stable in 20/20 folds in "
                 "all 7 H1 sections\n(stability 1.00 everywhere) and positive everywhere; "
                 "on the remapped mouse panel it is not.", loc="left")
    ax.legend(loc="upper left", fontsize=7.5)

    # ---------------- (c) does the instability transfer? ------------------
    ax = fig.add_subplot(gs[1, 1])
    sm = pd.read_csv(f"{R8}/d2_stability.csv")
    sh = pd.read_csv(f"{R9}/d2_stability.csv")
    pts = []
    r = sm[sm.pair == "raw_seed0 vs raw_seed1"].iloc[0]
    pts.append(("M1 denoise=False\nseed 0 vs 1 (20k)", r.pearson_r, r.top5_jaccard, C_M1, "o"))
    r = sm[sm.pair == "dca_seed0 vs dca_seed1"].iloc[0]
    pts.append(("M1 denoise=True\nseed 0 vs 1 (20k)", r.pearson_r, r.top5_jaccard, C_M1, "^"))
    r = sh[(sh.config == "denoise=False, FULL section")].iloc[0]
    pts.append(("H1 denoise=False\nseed 0 vs 1 (FULL, 196k)", r.pearson_r, r.top5_jaccard,
                C_H1, "s"))
    r = sh[(sh.config == "denoise=False") & (sh.seed_a == 0) & (sh.seed_b == 1)].iloc[0]
    pts.append(("H1 denoise=False\nseed 0 vs 1 (20k)", r.pearson_r, r.top5_jaccard,
                C_H1, "D"))
    r = sh[(sh.config == "denoise=True") & (sh.seed_a == 0) & (sh.seed_b == 1)].iloc[0]
    pts.append(("H1 denoise=True\nseed 0 vs 1 (20k)", r.pearson_r, r.top5_jaccard,
                C_H1, "^"))
    p = f"{RH}/deepscence_consensus_coverage.csv"
    if os.path.exists(p):
        cc = pd.read_csv(p)
        cc = cc[cc.status == "OK"]
        for r in cc.itertuples():
            pts.append(("H1 %s\n%d-seed consensus" % (r.section, r.n_seeds),
                        np.nan, r.jaccard_top5_mean, PAL.SERIES[3], "*"))
    OFFS = {0: (-9, -22), 1: (9, 4), 2: (9, 10), 3: (9, -16), 4: (9, 2)}
    for _i, (lab, pr, jc, colr, mk) in enumerate(pts):
        ROWS.append(dict(panel="c", label=lab.replace("\n", " "),
                         pearson_r=None if pr != pr else round(float(pr), 5),
                         top5_jaccard=round(float(jc), 5)))
        if pr == pr:
            ax.plot([pr], [jc], mk, color=colr, ms=9, mec=PAL.SURFACE, mew=1.3, zorder=3)
            dx, dy = OFFS.get(_i, (6, 6))
            ax.annotate(lab, (pr, jc), textcoords="offset points", xytext=(dx, dy),
                        fontsize=6.6, color=PAL.INK2,
                        ha="right" if dx < 0 else "left")
    ax.set_xlim(0, 1.12); ax.set_ylim(-0.06, 0.88)
    ax.plot([0, 1], [0, 1], "--", color=PAL.MUTED, lw=1)
    ax.set_xlabel("Pearson $r$ of the score between two seeds")
    ax.set_ylabel("Jaccard of the top-5 % call set")
    ax.set_title("(c)  The instability transfers, and gets worse natively.\n"
                 "Score agreement and call-set agreement are NOT the same quantity.",
                 loc="left")

    # ---------------- (d) P-vi: what denoising does to the depth loading ----
    ax = fig.add_subplot(gs[2, :])
    dd = pd.read_csv(f"{R9}/d2_depth.csv")
    dd = dd[dd.scope == "20,000-cell panel"].sort_values("section")
    m1_delta = [(0.3891, 0.6404, "7239"), (0.3176, 0.5314, "7259"),
                (0.4096, 0.5419, "7352")]        # results/phase8_d2/d2_depth.csv, config=="dca"
    xs = np.arange(len(m1_delta))
    ax.axhline(0, color=PAL.INK2, lw=1.2, zorder=1)
    for i, (a, b, lab) in enumerate(m1_delta):
        ax.plot([i, i], [0, b - a], "-", color=C_M1, lw=6, alpha=.55, solid_capstyle="butt",
                zorder=2)
        ax.plot([i], [b - a], "o", color=C_M1, ms=7, mec=PAL.SURFACE, mew=1.3, zorder=3)
        ROWS.append(dict(panel="d", arm="M1", section=lab, scope="full section",
                         rho_denoise_False=a, rho_denoise_True=b, delta_rho=round(b - a, 4)))
    xs2 = np.arange(len(m1_delta), len(m1_delta) + len(dd))
    for i, r in zip(xs2, dd.itertuples()):
        ax.plot([i, i], [0, r.delta_rho], "-", color=C_H1, lw=6, alpha=.55,
                solid_capstyle="butt", zorder=2)
        ax.plot([i], [r.delta_rho], "s", color=C_H1, ms=7, mec=PAL.SURFACE, mew=1.3, zorder=3)
        ROWS.append(dict(panel="d", arm="H1", section=r.section,
                         scope="20,000-cell panel",
                         rho_denoise_False=float(r.rho_denoise_False),
                         rho_denoise_True=float(r.rho_denoise_True),
                         delta_rho=float(r.delta_rho)))
    ax.axvline(len(m1_delta) - 0.5, color=PAL.AXIS, lw=1)
    ax.set_xticks(list(xs) + list(xs2))
    ax.set_xticklabels([l for _, _, l in m1_delta] + list(dd.section), fontsize=7)
    ax.set_ylabel(r"$\Delta\rho$ = $\rho$(score, counts) with denoise=True"
                  "\n" r"minus without")
    ax.plot([], [], "o-", color=C_M1, lw=4, alpha=.6, label="M1, full sections (3), "
            "ortholog-remapped")
    ax.plot([], [], "s-", color=C_H1, lw=4, alpha=.6, label="H1, 20,000-cell panel (7), "
            "native")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("(d)  P-vi, FALSIFIED.  Registered: denoise=True RAISES the depth loading, "
                 "as it does on 3 of 3 mouse sections.  Falsifier: $\\Delta\\rho \\leq 0$ in "
                 "$\\geq$ 5 of 7.\n"
                 "Measured on H1: $\\Delta\\rho \\leq 0$ in 6 of 7, and the SIGN of the "
                 "loading inverts in 4 of 7.  The mouse direction does not transfer.",
                 loc="left")

    fig.suptitle("Figure 6 — DeepScence, native human panel against ortholog-remapped "
                 "mouse panel.  Every number carries coverage / denoise / anchor / panel / "
                 "seed configuration.\n"
                 "Panels (a) and (b) on H1 are at random_state = 0, which under PI decision "
                 "D-A is no longer H1's primary DeepScence estimator: directions are "
                 "seed-robust, magnitudes are not.", fontsize=9.6, y=1.005)
    for e in ("png", "pdf"):
        fig.savefig(FIG + "figure6." + e, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(ROWS).to_csv(FIG + "figure6_data.csv", index=False)
    print("wrote figures/figure6.png/.pdf and figure6_data.csv (%d rows)" % len(ROWS))


if __name__ == "__main__":
    main()
