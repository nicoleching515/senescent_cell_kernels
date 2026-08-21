#!/usr/bin/env python3
"""Figures 2b, 2c (Master Plan Section 25) and 2d (Phase 3 addition).

2b  binned response vs distance-to-nearest-sender, with the matched-decoy curve
    and the torus-shift null band on the same axes.
2c  surviving fraction of beta under each null.
2d  log(median distance to nearest sender) vs log(sender density), against the
    homogeneous-Poisson line of slope -1/2 -- i.e. whether the model's
    independent variable is a measurement or a sender-calling rate.
"""
from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/workspace/code")
import sasp_palette as PAL
import sasp_phase3 as P
import run_phase3_nulls as RN

PAL.apply_style(matplotlib)
RES = P.RESULTS
FIG = "/workspace/figures"


def lab(s):
    m = P.R.parse_sample(s) if hasattr(P, "R") else None
    parts = s.split("_")
    return f"{parts[2]} {parts[-1].split('-')[0]} wk ({parts[0]})"


def fig2b(modules=("tnfa_nfkb_proximal", "secondary_senescence"),
          celltype="Hepatocytes"):
    cur = pd.read_csv(f"{RES}/curves.csv")
    per = pd.read_csv(f"{RES}/perm_curves.csv") if os.path.exists(
        f"{RES}/perm_curves.csv") else pd.DataFrame()
    order = [s for s in P.IN_BAND if s in set(cur.section)]
    fig, axes = plt.subplots(len(modules), len(order),
                             figsize=(2.45 * len(order), 2.6 * len(modules)),
                             sharex=True)
    axes = np.atleast_2d(axes)
    for r, mod in enumerate(modules):
        for c, sec in enumerate(order):
            ax = axes[r, c]
            q = cur[(cur.section == sec) & (cur.celltype == celltype)
                    & (cur.module == mod)]
            if q.empty:
                ax.set_axis_off()
                continue
            x = 0.5 * (q.bin_lo + q.bin_hi)
            b = per[(per.section == sec) & (per.celltype == celltype)
                    & (per.module == mod) & (per.null == "N3")] \
                if not per.empty else pd.DataFrame()
            if not b.empty:
                bx = 0.5 * (b.bin_lo + b.bin_hi)
                ax.fill_between(bx, b["lo"], b["hi"], color=PAL.MUTED,
                                alpha=0.35, lw=0,
                                label="torus-shift null (N3) 95 %")
            ax.plot(x, q.mean_decoy, color=PAL.SERIES[1], lw=1.6, ls="--",
                    label="matched decoy (N2)")
            ax.errorbar(x, q.mean_obs, yerr=q.sem_obs, color=PAL.SERIES[0],
                        lw=1.6, marker="o", ms=2.6, capsize=0,
                        label="observed senders")
            ax.axhline(0, color=PAL.AXIS, lw=0.6)
            if r == 0:
                ax.set_title(lab(sec), fontsize=8)
            if c == 0:
                ax.set_ylabel(f"{mod}\nmean module score", fontsize=7.5)
            if r == len(modules) - 1:
                ax.set_xlabel("distance to nearest sender (µm)", fontsize=8)
            ax.set_xlim(0, RN.WINDOW_UM)
    axes[0, 0].legend(loc="upper right", fontsize=6)
    fig.suptitle("Figure 2b — the naive gradient, the matched-decoy curve and "
                 "the torus-shift null band\n"
                 f"{celltype}, six Test-3-admissible sections, sender call "
                 f"{RN.PRIMARY_CALL}, 5 µm bins to {RN.WINDOW_UM:.0f} µm",
                 fontsize=9.5, y=1.01)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2b.{e}", bbox_inches="tight")
    print("wrote figure2b")


NULL_LABELS = [
    ("sf_n2", "N2  matched decoy"),
    ("N3_sf", "N3  torus shift"),
    ("N4_sf", "N4  rotation"),
    ("sf_zon", "zonation covariate only"),
    ("sf_n6", "N6  receiver baseline"),
    ("N1_sf", "N1  label permutation"),
    ("sf_n5", "N5  nuisance conditioning"),
    ("sf_n6n5", "N5 + N6"),
    ("sf_n2n5n6", "N2 + N5 + N6"),
    ("sf_n8", "N8  scrambled response set"),
]


def fig2c(call=RN.PRIMARY_CALL):
    mf = pd.read_csv(f"{RES}/main_fits.csv")
    pf = pd.read_csv(f"{RES}/perm_nulls.csv") if os.path.exists(
        f"{RES}/perm_nulls.csv") else pd.DataFrame()
    n8p = [f"{RES}/n8_scrambled_{s}.csv" for s in P.IN_BAND]
    n8p = [p for p in n8p if os.path.exists(p)]
    n8 = pd.concat([pd.read_csv(p) for p in n8p], ignore_index=True) \
        if n8p else pd.DataFrame()
    key = ["section", "celltype", "module"]
    d = mf[mf.section.isin(P.IN_BAND) & (mf.call == call)
           & (mf.stratum == "all")].copy()
    if not pf.empty:
        d = d.merge(pf[key + ["N1_sf", "N3_sf", "N4_sf"]].drop_duplicates(key),
                    on=key, how="left")
    if not n8.empty:
        d = d.merge(n8[key + ["sf_n8"]].drop_duplicates(key), on=key, how="left")
    rep = d[(d.beta_naive > 0) & (d.beta_base_lo > 0)]
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ys = np.arange(len(NULL_LABELS))[::-1]
    rows, yticks, ylabels = [], [], []
    rng = np.random.default_rng(1)
    for y, (col, l) in zip(ys, NULL_LABELS):
        if col not in rep:
            continue
        v = rep[col].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        jit = (rng.random(v.size) - .5) * .36
        ax.scatter(np.clip(v, -1.15, 2.15), y + jit, s=12, alpha=.4,
                   color=PAL.SERIES[0], lw=0)
        q = np.quantile(v, [.25, .5, .75])
        ax.plot(q[[0, 2]], [y, y], color=PAL.INK, lw=3.2, solid_capstyle="butt")
        ax.plot([q[1]], [y], marker="|", ms=17, color=PAL.SERIES[3], mew=3)
        rows.append(dict(null=l, n=int(v.size), q25=q[0], median=q[1],
                         q75=q[2], frac_le_0=float((v <= 0).mean()),
                         frac_gt_05=float((v > .5).mean())))
        yticks.append(y)
        ylabels.append(f"{l}   (med {q[1]:+.2f})")
    ax.axvline(1.0, color=PAL.MUTED, lw=1, ls=":")
    ax.axvline(0.0, color=PAL.STATUS["critical"], lw=1.2)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8.5)
    ax.set_xlabel("surviving fraction of β̂   "
                  "(1 = the null removes nothing, 0 = it removes all of it)")
    ax.set_xlim(-1.2, 2.2)
    ax.set_title("Figure 2c — surviving fraction under each null\n"
                 f"six Test-3-admissible sections, {call}, receiver type × "
                 f"Tier B module; n = {len(rep)} of {len(d)} fits with a "
                 "positive, bootstrap-nonzero naive amplitude", fontsize=9.5)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2c.{e}", bbox_inches="tight")
    pd.DataFrame(rows).to_csv(f"{RES}/figure2c_data.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print("wrote figure2c")


def fig2d():
    d = pd.read_csv(f"{RES}/poisson_density.csv")
    f = pd.read_csv(f"{RES}/poisson_fits.csv")
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    calls = sorted(d.call.unique())
    for i, c in enumerate(calls):
        g = d[d.call == c]
        ax.scatter(g.sender_density_per_um2, g.median_d_um, s=34,
                   color=PAL.SERIES[i % len(PAL.SERIES)], label=c, lw=0,
                   alpha=.85)
    x = np.logspace(np.log10(d.sender_density_per_um2.min() * .8),
                    np.log10(d.sender_density_per_um2.max() * 1.2), 50)
    ax.plot(x, np.sqrt(np.log(2) / (np.pi * x)), color=PAL.INK, lw=1.6,
            ls="--", label="homogeneous Poisson,  slope −1/2")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("sender density (called senders per µm²)")
    ax.set_ylabel("median distance to nearest sender (µm)")
    r = f[f.subset.str.startswith("ALL")].iloc[0]
    ax.set_title("Figure 2d — the kernel's independent variable is a\n"
                 "sender-calling rate\n"
                 f"{int(r.n)} section × sender-definition combinations:  "
                 f"slope {r.slope:+.3f},  r² = {r.r2:.3f}", fontsize=9.5)
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2d.{e}", bbox_inches="tight")
    print("wrote figure2d")


if __name__ == "__main__":
    if os.path.exists(f"{RES}/curves.csv"):
        fig2b()
    if os.path.exists(f"{RES}/main_fits.csv"):
        fig2c()
    if os.path.exists(f"{RES}/poisson_density.csv"):
        fig2d()
