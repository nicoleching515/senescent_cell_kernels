"""Figure 2a -- the naive decay curve on real tissue (Master Plan Sec 25).

"Reproduce the standard analysis, deliberately without controls, so you have the
thing you are going to test."  Panels a/b are that curve.  Panels c/d are what
Phase 1 says you must show next to it: which kernel family you chose, how far the
data actually reach, and how badly the standard iid CI understates uncertainty.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import sasp_palette as P
import sasp_kernels as K

P.apply_style(matplotlib)
FIG, RES = "/workspace/figures", "/workspace/results"
FAM_ORDER = ["exponential", "gaussian", "powerlaw", "step", "spline"]
FAM_COL = {f: P.SERIES[i] for i, f in enumerate(FAM_ORDER)}
FAM_MK = dict(zip(FAM_ORDER, ["o", "s", "^", "D", "v"]))
XMAX_UM, XMAX_NN = 130, 15


def short(s):
    p = s.split("_")
    return f"{p[0]} {p[2]} {p[4].split('-')[0]}wk"


def panel(ax, letter, text, dx=-0.22, dy=1.30):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", color=P.INK, va="bottom")
    ax.text(dx, dy, "     " + text, transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", color=P.INK2, va="bottom")


def main():
    fits = pd.read_csv(f"{RES}/real_fits.csv")
    fits = fits[fits.celltype == "ALL"]
    curves = pd.read_csv(f"{RES}/real_curves.csv")
    curves = curves[(curves.celltype == "ALL") & (curves["n"] >= 50)]
    um = fits[fits.units == "um"]
    adj = fits[fits.units == "um_covadj"]

    mods = sorted(curves.module.unique())
    samples = sorted(curves["sample"].unique())
    scol = {s: P.SERIES[i % 8] for i, s in enumerate(samples)}
    smk = {s: ["o", "s", "^", "D", "v", "P"][i % 6] for i, s in enumerate(samples)}
    ncol, nrow = 4, int(np.ceil(len(mods) / 4))

    fig = plt.figure(figsize=(15.2, 2.4 + 2.55 * 2 * nrow + 6.6))
    gs = GridSpec(2 * nrow + 2, 4, figure=fig, hspace=0.55, wspace=0.34,
                  left=0.065, right=0.975, top=0.885, bottom=0.045,
                  height_ratios=[1] * (2 * nrow) + [1.3, 1.3])

    sender = um.sender_call.iloc[0]
    fig.text(0.04, 0.982, "Figure 2a  ·  The naive distance-to-nearest-sender "
             "gradient on real liver tissue", fontsize=15.5,
             fontweight="bold", color=P.INK)
    fig.text(0.04, 0.966,
             f"Xenium Prime Mouse 5K · {len(samples)} liver sections from "
             f"{len(samples)} mice · sender call {sender} (5 % prevalence) · "
             "Tier B module scores · NO controls applied.", fontsize=9.5,
             color=P.INK2)
    fig.text(0.04, 0.951,
             "Receiver cell-type stratification NOT applied "
             "(Bio celltypes_*.csv pending) — all fits are on a single "
             "unstratified stratum.", fontsize=9,
             color=P.STATUS["serious"], fontweight="bold")
    fig.text(0.04, 0.936,
             "99 % of cells lie within 72–90 µm of a sender, so the plan's "
             "300 µm binning window is unreachable at this prevalence; panels "
             "are drawn over the range the data actually cover.",
             fontsize=9, color=P.INK2)

    # ---- (a) microns, (b) NN units ---------------------------------------
    for row, (xcol, xmax, xlab, tag) in enumerate([
            ("bin_mean_d_um", XMAX_UM, "distance to nearest sender (µm)", "a"),
            ("bin_nn_units", XMAX_NN, "distance (local median NN units)", "b")]):
        for i, mod in enumerate(mods):
            ax = fig.add_subplot(gs[row * nrow + i // ncol, i % ncol])
            for s in samples:
                c = curves[(curves.module == mod) & (curves["sample"] == s)]
                c = c[c[xcol] <= xmax]
                if not len(c):
                    continue
                ax.plot(c[xcol], c["mean"], marker=smk[s], ms=3.2,
                        color=scol[s], lw=1.7, label=short(s))
                if row == 0:
                    ax.fill_between(c[xcol], c["mean"] - 1.96 * c["sem"],
                                    c["mean"] + 1.96 * c["sem"],
                                    color=scol[s], alpha=0.15, linewidth=0)
            ax.set_title(mod.replace("_", " "), pad=5, fontsize=9, color=P.INK)
            ax.set_xlim(0, xmax)
            if i % ncol == 0:
                ax.set_ylabel("module score")
            ax.set_xlabel(xlab, fontsize=8)
            if i == 0:
                ax.legend(fontsize=6.4, loc="lower left")
                panel(ax, tag, "binned response vs distance (microns)"
                      if tag == "a" else
                      "same curves, distance in local median-NN units")

    # ---- (c1) families on a representative interior fit -------------------
    r0 = um[(um.family == "exponential") & (um.lam_at_bound == 0)]
    rep = r0.iloc[int(np.argsort(r0.n_cells.values)[-1])]
    rep_s, rep_m = rep["sample"], rep["module"]
    ax = fig.add_subplot(gs[2 * nrow, 0:2])
    c = curves[(curves.module == rep_m) & (curves["sample"] == rep_s)]
    c = c[c.bin_mean_d_um <= XMAX_UM]
    x, ybin, w = c.bin_mean_d_um.values, c["mean"].values, c["n"].values
    ax.plot(x, ybin, "o", ms=5, color=P.MUTED, lw=0, label="binned data",
            zorder=5)
    dd = np.linspace(0.5, XMAX_UM, 400)
    f = um[(um.module == rep_m) & (um["sample"] == rep_s)]
    for fam in FAM_ORDER:
        rr = f[f.family == fam]
        if not len(rr) or fam == "spline":
            continue
        rr = rr.iloc[0]
        pr = dict(lam=rr["lam"])
        if fam == "powerlaw":
            pr["p"] = rr["pow_p"]
        # rescale the fitted shape to the binned curve by weighted LS, so the
        # comparison is about SHAPE, not about the nuisance intercept
        kx = K.kernel_curve(fam, pr, x)
        A = np.column_stack([np.ones_like(kx), kx])
        cf = np.linalg.lstsq(A * np.sqrt(w)[:, None], ybin * np.sqrt(w),
                             rcond=None)[0]
        ax.plot(dd, cf[0] + cf[1] * K.kernel_curve(fam, pr, dd),
                color=FAM_COL[fam], lw=2,
                label=f"{fam}: d½={rr['d_half']:.0f} µm")
    ax.set_xlabel("distance to nearest sender (µm)")
    ax.set_ylabel("module score")
    ax.set_title(f"{short(rep_s)} · {rep_m.replace('_',' ')} "
                 f"(n={int(rep.n_cells):,})", pad=6, color=P.INK)
    ax.legend(fontsize=7.5)
    panel(ax, "c", "the family you pick changes the answer", dx=-0.11, dy=1.14)

    # ---- (c2) AIC ---------------------------------------------------------
    ax = fig.add_subplot(gs[2 * nrow, 2])
    g = um.copy()
    g["daic"] = g.groupby(["sample", "module"])["aic"].transform(lambda v: v - v.min())
    mean_d = g.groupby("family")["daic"].mean().reindex(FAM_ORDER)
    wins = (g.loc[g.groupby(["sample", "module"])["aic"].idxmin(), "family"]
            .value_counts().reindex(FAM_ORDER).fillna(0).astype(int))
    xs = np.arange(len(FAM_ORDER))
    ax.bar(xs, mean_d.values, color=[FAM_COL[f] for f in FAM_ORDER], zorder=3)
    top = float(np.nanmax(mean_d.values))
    for x_, f_ in zip(xs, FAM_ORDER):
        ax.text(x_, mean_d.values[x_] + top * 0.03, f"{wins[f_]}", ha="center",
                fontsize=8, color=P.INK, fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels(FAM_ORDER, rotation=30, ha="right")
    ax.set_ylabel("mean ΔAIC vs best")
    ax.set_title("model comparison\n(number above bar = fits won, of 35)",
                 pad=6, color=P.INK, fontsize=9)

    # ---- (c3) covariate-adjusted preview of N5 ----------------------------
    ax = fig.add_subplot(gs[2 * nrow, 3])
    a = adj[adj.family == "exponential"].set_index(["sample", "module"])["beta"].abs()
    b = um[um.family == "exponential"].set_index(["sample", "module"])["beta"].abs()
    ratio = (a / b).dropna()
    order = ratio.groupby("module").median().sort_values()
    ax.barh(np.arange(len(order)), order.values, color=P.SERIES[1], zorder=3)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([m.replace("_", " ") for m in order.index], fontsize=7.5)
    ax.axvline(1.0, color=P.MUTED, ls=":", lw=1.6)
    ax.set_xlabel(r"$|\hat\beta_{adjusted}|\ /\ |\hat\beta_{naive}|$")
    ax.set_title("preview of null N5: adjust for counts,\narea and local "
                 "density", pad=6, color=P.INK, fontsize=9)

    # ---- (d1) d_half by family, at-bound flagged --------------------------
    ax = fig.add_subplot(gs[2 * nrow + 1, 0:2])
    for i, fam in enumerate(FAM_ORDER):
        sub = um[um.family == fam]
        ib = sub.lam_at_bound.fillna(0).astype(bool).values
        xj = np.full(len(sub), i) + np.linspace(-0.26, 0.26, max(len(sub), 1))
        lo = np.maximum(sub.d_half - sub.ci_dhalf_blk_lo, 0)
        hi = np.maximum(sub.ci_dhalf_blk_hi - sub.d_half, 0)
        ax.errorbar(xj[~ib], sub.d_half.values[~ib],
                    yerr=[lo.values[~ib], hi.values[~ib]], fmt=FAM_MK[fam],
                    ms=5, color=FAM_COL[fam], lw=0, capsize=2, elinewidth=1,
                    label="interior optimum" if i == 0 else None)
        ax.scatter(xj[ib], sub.d_half.values[ib], marker="x", s=42,
                   color=P.STATUS["critical"], zorder=5,
                   label="λ̂ at grid bound (unidentified)" if i == 0 else None)
    ax.axhspan(0, 10, color=P.STATUS["warning"], alpha=0.16, zorder=0)
    ax.text(0.02, 0.02, "shaded: below the 6.7–9.8 µm resolution floor",
            transform=ax.transAxes, fontsize=7.5, color=P.INK2)
    ax.set_xticks(range(len(FAM_ORDER)))
    ax.set_xticklabels(FAM_ORDER, rotation=25, ha="right")
    ax.set_ylabel("$d_{1/2}$ (µm)")
    ax.set_title("half-decay distance, every section × module "
                 "(error bars: spatial block bootstrap)", pad=6, color=P.INK)
    ax.legend(fontsize=7.5, loc="upper right")
    panel(ax, "d", "length scale and its honest uncertainty",
          dx=-0.11, dy=1.14)

    # ---- (d2) SE understatement, interior fits only -----------------------
    ax = fig.add_subplot(gs[2 * nrow + 1, 2:])
    sub = um[(um.family != "spline") & np.isfinite(um.se_ratio)
             & (um.lam_at_bound == 0)]
    fams = [f for f in FAM_ORDER if f != "spline"]
    for i, fam in enumerate(fams):
        v = sub[sub.family == fam].se_ratio.values
        if not len(v):
            continue
        ax.scatter(np.full(len(v), i) + np.linspace(-0.22, 0.22, len(v)), v,
                   color=FAM_COL[fam], s=30, marker=FAM_MK[fam], zorder=3)
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=P.INK, lw=2,
                zorder=4)
    ax.axhline(1.0, color=P.MUTED, ls=":", lw=1.6)
    ax.axhline(7.9, color=P.STATUS["serious"], ls="--", lw=1.4)
    ax.text(3.45, 8.6, "Phase 1 synthetic max 7.9×", fontsize=7.5,
            ha="right", color=P.STATUS["serious"])
    ax.set_yscale("log")
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels(fams, rotation=25, ha="right")
    ax.set_ylabel("block bootstrap sd ÷ iid SE")
    ax.set_title("SE understatement factor of the standard iid CI\n"
                 "(interior optima only; black bar = median)",
                 pad=6, color=P.INK)

    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/figure2a.png", bbox_inches="tight")
    fig.savefig(f"{FIG}/figure2a.pdf", bbox_inches="tight")
    curves.to_csv(f"{FIG}/figure2a_curves.csv", index=False)
    fits.to_csv(f"{FIG}/figure2a_fits.csv", index=False)
    print("wrote figure2a.png/.pdf + figure2a_curves.csv + figure2a_fits.csv")


if __name__ == "__main__":
    main()
