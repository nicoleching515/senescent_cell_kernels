"""Figure 1 -- Identifiability regime map (synthetic).  Master Plan Section 25.

(a) Bias in lambda_hat as a function of sender clustering and baseline
    autocorrelation length.
(b) CI coverage across the same grid.
(c) Recovery with and without the matched-decoy control.

Design follows the dataviz method: signed bias uses the blue<->red diverging
pair with a NEUTRAL GRAY midpoint at zero bias; coverage uses a single-hue
sequential ramp; estimator identity uses categorical slots in fixed order,
never cycled, always with a legend and direct labels so identity is never
carried by colour alone.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullLocator

import sasp_palette as P

P.apply_style(matplotlib)

RES = "/workspace/results"
FIG = "/workspace/figures"
LAM_TRUE = 30.0
KAPPA = [0.0, 1.5, 3.0, 4.5]
RATIO = [0.25, 0.5, 1.0, 2.0, 4.0]

# estimator display order = fixed categorical slot order
EST = [("naive", "naive", P.SERIES[0], "o"),
       ("decoyS", "matched-decoy (N2)", P.SERIES[1], "s"),
       ("nuis", "nuisance-conditioned (N5)", P.SERIES[2], "^"),
       ("decoyS_nuis", "N2 + N5", P.SERIES[3], "D")]


def grid_of(df, col, agg="mean"):
    g = np.full((len(KAPPA), len(RATIO)), np.nan)
    for i, k in enumerate(KAPPA):
        for j, r in enumerate(RATIO):
            s = df[(df.clustering == k) & (np.isclose(df.ell_over_lambda, r))][col]
            if len(s):
                g[i, j] = s.mean() if agg == "mean" else s.median()
    return g


def logx(ax, ticks):
    """Log x-axis with only the swept values labelled -- matplotlib's minor
    decade labels otherwise collide with them."""
    ax.set_xscale("log")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])


def panel(ax, letter, text):
    ax.text(-0.24, 1.30, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", color=P.INK, va="bottom", ha="left")
    ax.text(-0.24, 1.30, "     " + text, transform=ax.transAxes,
            fontsize=10.5, fontweight="bold", color=P.INK2,
            va="bottom", ha="left")


def heat(ax, G, cmap, vmin, vmax, title, fmt="{:+.0%}", ref=None):
    im = ax.imshow(G, cmap=cmap, vmin=vmin, vmax=vmax, origin="lower",
                   aspect="auto")
    ax.set_xticks(range(len(RATIO)))
    ax.set_xticklabels([f"{r:g}" for r in RATIO])
    ax.set_yticks(range(len(KAPPA)))
    ax.set_yticklabels([f"{k:g}" for k in KAPPA])
    ax.set_title(title, pad=6, color=P.INK)
    ax.grid(False)
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            if not np.isfinite(G[i, j]):
                continue
            v = G[i, j]
            rel = (v - vmin) / (vmax - vmin + 1e-12)
            col = "#ffffff" if (rel > 0.78 or rel < 0.22) else P.INK
            ax.text(j, i, fmt.format(v), ha="center", va="center",
                    fontsize=7.5, color=col)
    return im


def main():
    df = pd.read_csv(f"{RES}/sweep_all.csv")
    df["sweep"] = df["sweep"].fillna("null")   # "null" is read as NaN
    main_df = df[df.sweep == "main"].copy()
    curves = pd.read_csv(f"{RES}/figure1c_curves.csv")

    fig = plt.figure(figsize=(13.6, 13.6))
    gs = GridSpec(4, 6, figure=fig, hspace=0.78, wspace=0.62,
                  height_ratios=[1.0, 1.0, 1.0, 0.95],
                  left=0.075, right=0.965, top=0.895, bottom=0.05)

    fig.text(0.045, 0.972,
             "Figure 1  ·  Identifiability regime map for the SASP spatial "
             "response kernel (synthetic tissue)", fontsize=13.5,
             fontweight="bold", color=P.INK)
    fig.text(0.045, 0.953,
             f"lambda_true = {LAM_TRUE:g} um, beta_true = 1.0, sender "
             "prevalence 5%, confounder on. 30 seeds per grid cell; "
             "~11,900 cells per section, median NN distance 11.7 um.",
             fontsize=9, color=P.INK2)

    # ---------------- (a) bias in lambda_hat ----------------------------
    cols_a = [("relbias_lam_naive", "naive"),
              ("relbias_lam_decoyS", "matched-decoy (N2)"),
              ("relbias_lam_nuis", "nuisance-conditioned (N5)")]
    vlim = 1.2
    for c, (col, name) in enumerate(cols_a):
        ax = fig.add_subplot(gs[0, 2 * c:2 * c + 2])
        G = grid_of(main_df, col)
        im = heat(ax, G, P.DIV, -vlim, vlim, name)
        if c == 0:
            ax.set_ylabel("sender clustering  $\\kappa$")
            panel(ax, "a", "Relative bias in $\\hat\\lambda$")
        ax.set_xlabel("baseline autocorrelation  $\\ell/\\lambda_{true}$")
    cb = fig.colorbar(im, ax=fig.axes[-1], fraction=0.035, pad=0.03)
    cb.set_label("bias in $\\hat\\lambda$", fontsize=8)
    cb.outline.set_visible(False)

    # ---------------- (b) CI coverage -----------------------------------
    cols_b = [("cover_lam_naive_iid", "naive, iid asymptotic CI"),
              ("cover_lam_naive_blk", "naive, spatial block bootstrap"),
              ("cover_lam_decoyS_blk", "matched-decoy, block bootstrap")]
    for c, (col, name) in enumerate(cols_b):
        ax = fig.add_subplot(gs[1, 2 * c:2 * c + 2])
        G = grid_of(main_df, col)
        im = heat(ax, G, P.SEQ, 0.0, 1.0, name, fmt="{:.2f}")
        if c == 0:
            ax.set_ylabel("sender clustering  $\\kappa$")
            panel(ax, "b", "95% CI coverage of $\\lambda$  (nominal 0.95)")
        ax.set_xlabel("baseline autocorrelation  $\\ell/\\lambda_{true}$")
    cb = fig.colorbar(im, ax=fig.axes[-1], fraction=0.035, pad=0.03)
    cb.set_label("coverage", fontsize=8)
    cb.outline.set_visible(False)

    # ---------------- (c) recovery with / without decoy control ---------

    # c1: binned response curves, easy vs hard
    ax = fig.add_subplot(gs[2, 0:2])
    regs = sorted(curves.regime.unique())
    easy = [r for r in regs if r.startswith("easy")][0]
    hard = [r for r in regs if r.startswith("hard ")][0]
    for reg, ls, tag in [(easy, "-", "easy"), (hard, "--", "hard")]:
        cc = curves[(curves.regime == reg) & (curves.bin_um <= 200)
                    & (curves.n_cells >= 40)]
        ax.plot(cc.bin_um, cc.sender_curve, ls, color=P.SERIES[0], lw=2,
                label=f"sender, {tag}")
        ax.plot(cc.bin_um, cc.decoy_curve, ls, color=P.SERIES[1], lw=2,
                label=f"matched decoy, {tag}")
    cc = curves[(curves.regime == easy) & (curves.bin_um <= 200)
                & (curves.n_cells >= 40)]
    ax.plot(cc.bin_um, cc.truth, ":", color=P.MUTED, lw=1.8,
            label="planted kernel")
    ax.set_xlim(0, 200)
    ax.set_xlabel("distance to nearest sender / decoy (um)")
    ax.set_ylabel("mean response (far-field centred)")
    ax.set_title("binned curves: sender vs matched decoy", pad=6, color=P.INK)
    ax.legend(fontsize=7, loc="upper right")
    panel(ax, "c", "Recovery with and without the matched-decoy control")

    # c2: lambda_hat vs ell/lambda per estimator (kappa pooled)
    ax = fig.add_subplot(gs[2, 2:4])
    for key, lab, col, mk in EST:
        m = main_df.groupby("ell_over_lambda")[f"lam_{key}"].mean()
        lo = main_df.groupby("ell_over_lambda")[f"lam_{key}"].quantile(0.25)
        hi = main_df.groupby("ell_over_lambda")[f"lam_{key}"].quantile(0.75)
        ax.plot(m.index, m.values, marker=mk, color=col, label=lab)
        ax.fill_between(m.index, lo.values, hi.values, color=col, alpha=0.13,
                        linewidth=0)
    ax.axhline(LAM_TRUE, color=P.MUTED, ls=":", lw=1.6)
    ax.annotate("$\\lambda_{true}$", xy=(0.27, LAM_TRUE * 1.06), fontsize=8,
                color=P.INK2)
    logx(ax, RATIO)
    ax.set_xlabel("baseline autocorrelation  $\\ell/\\lambda_{true}$")
    ax.set_ylabel("$\\hat\\lambda$ (um)")
    ax.set_title("$\\hat\\lambda$ by estimator (pooled over $\\kappa$)",
                 pad=6, color=P.INK)
    ax.legend(loc="upper left", fontsize=7.5)

    # c3: signal retained -- the Section 29 objection-6 number
    ax = fig.add_subplot(gs[2, 4:6])
    clean = df[(df.sweep == "clean")]
    groups = [("clean\n(no confounder)", clean[clean.regime == "easy"]),
              ("$\\ell/\\lambda$ = 1", main_df[np.isclose(main_df.ell_over_lambda, 1.0)]),
              ("$\\ell/\\lambda$ = 4", main_df[np.isclose(main_df.ell_over_lambda, 4.0)])]
    w = 0.2
    for e, (key, lab, col, mk) in enumerate(EST):
        xs = np.arange(len(groups)) + (e - 1.5) * w
        vals = [g[f"signal_retained_{key}"].mean() for _, g in groups]
        errs = [g[f"signal_retained_{key}"].std() for _, g in groups]
        ax.bar(xs, vals, width=w * 0.88, color=col, label=lab, zorder=3)
        ax.errorbar(xs, vals, yerr=errs, fmt="none", ecolor=P.INK2, lw=1,
                    capsize=2, zorder=4)
    ax.axhline(1.0, color=P.MUTED, ls=":", lw=1.6)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("$\\hat\\beta / \\beta_{true}$")
    ax.set_title("signal retained: $\\hat\\beta/\\beta_{true}$, 1.0 = unbiased",
                 pad=6, color=P.INK)
    ax.set_ylim(0, 3.1)
    ax.legend(loc="upper left", fontsize=7, ncol=2)

    # ---------------- (d) supporting sweeps ------------------------------

    # d1: prevalence
    ax = fig.add_subplot(gs[3, 0:2])
    pv = df[df.sweep == "prev"]
    for key, lab, col, mk in EST[:3]:
        m = pv.groupby("prevalence")[f"lam_{key}"].mean()
        ax.plot(m.index * 100, m.values, marker=mk, color=col, label=lab)
    ax.axhline(LAM_TRUE, color=P.MUTED, ls=":", lw=1.6)
    ax.axvspan(2, 10, color=P.STATUS["good"], alpha=0.09, zorder=0)
    ax.annotate("plan's sweet spot\n2-10%", xy=(4.4, 0.60),
                xycoords=("data", "axes fraction"), fontsize=7,
                color=P.INK2, ha="center")
    logx(ax, [0.5, 1, 2, 5, 10, 20, 30])
    ax.set_xlabel("sender prevalence (%)")
    ax.set_ylabel("$\\hat\\lambda$ (um)")
    ax.set_title("sender prevalence (Section 8 Test 3)", pad=6, color=P.INK)
    ax.legend(fontsize=7, loc="upper left")
    panel(ax, "d", "Supporting sweeps")

    # d2: confounder strength
    ax = fig.add_subplot(gs[3, 2:4])
    cf = df[df.sweep == "conf"]
    for key, lab, col, mk in EST[:3]:
        m = cf.groupby("conf_strength")[f"signal_retained_{key}"].mean()
        ax.plot(m.index, m.values, marker=mk, color=col, label=lab)
    ax.axhline(1.0, color=P.MUTED, ls=":", lw=1.6)
    ax.set_xlabel("confounder strength")
    ax.set_ylabel("$\\hat\\beta / \\beta_{true}$")
    ax.set_title("confounder strength ($\\kappa$=3, $\\ell/\\lambda$=2)",
                 pad=6, color=P.INK)
    ax.legend(fontsize=7)

    # d3: over-confidence ratio + SMD
    ax = fig.add_subplot(gs[3, 4:6])
    m = main_df.groupby("ell_over_lambda")["se_ratio_naive"].mean()
    ax.plot(m.index, m.values, marker="o", color=P.SERIES[0],
            label="block-bootstrap sd / iid SE")
    ax.axhline(1.0, color=P.MUTED, ls=":", lw=1.6)
    ax2 = ax.twinx()
    ms = main_df.groupby("ell_over_lambda")["max_smd_after"].mean()
    ax2.plot(ms.index, ms.values, marker="s", color=P.SERIES[1],
             label="max SMD after matching")
    ax2.axhline(0.1, color=P.SERIES[1], ls="--", lw=1.2, alpha=0.6)
    ax2.set_ylim(0, 0.16); ax2.grid(False)
    ax2.set_ylabel("max SMD after matching", color=P.SERIES[1], fontsize=8.5)
    ax2.tick_params(axis="y", colors=P.SERIES[1])
    logx(ax, RATIO)
    ax.set_xlabel("baseline autocorrelation  $\\ell/\\lambda_{true}$")
    ax.set_ylabel("SE understatement factor", color=P.SERIES[0], fontsize=8.5)
    ax.tick_params(axis="y", colors=P.SERIES[0])
    ax.set_title("over-confidence, and matching balance", pad=6, color=P.INK)
    ax.annotate("SMD < 0.1 bar\n(Section 8 Test 5)", xy=(0.55, 0.105),
                xycoords=("data", "axes fraction"), fontsize=7,
                color=P.SERIES[1])

    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/figure1.png", bbox_inches="tight")
    fig.savefig(f"{FIG}/figure1.pdf", bbox_inches="tight")
    print(f"wrote {FIG}/figure1.png / .pdf")

    # ---- the underlying tidy CSV ---------------------------------------
    keep = ["clustering", "ell_over_lambda", "ripley50"]
    metrics = []
    for key, lab, _, _ in EST:
        metrics += [f"lam_{key}", f"beta_{key}", f"relbias_lam_{key}",
                    f"signal_retained_{key}"]
    metrics += ["cover_lam_naive_iid", "cover_lam_naive_blk",
                "cover_lam_nuis_blk", "cover_lam_decoyS_blk",
                "cover_beta_naive_blk", "cover_beta_decoyS_blk",
                "se_ratio_naive", "max_smd_before", "max_smd_after",
                "smd_pass", "u_sender_mean", "u_decoy_mean", "u_all_mean",
                "n_cells", "n_senders", "median_nn_um", "med_d_sender"]
    metrics = [m for m in metrics if m in main_df.columns]
    tab = (main_df.groupby(["clustering", "ell_over_lambda"])[metrics]
           .agg(["mean", "std"]).round(4))
    tab.columns = ["_".join(c) for c in tab.columns]
    tab.reset_index().to_csv(f"{FIG}/figure1_data.csv", index=False)
    print(f"wrote {FIG}/figure1_data.csv")


if __name__ == "__main__":
    main()
