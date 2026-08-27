#!/usr/bin/env python3
"""Figures 2b, 2c (Master Plan Section 25) and 2d (Phase 3 addition).

2b  binned response vs distance-to-nearest-sender, with the matched-decoy curve
    and the torus-shift null bands (bounding-box and in-tissue) on the same axes.
2c  surviving fraction of beta under each null, including every C1 in-tissue
    variant and the FROZEN PRIMARY variance-corrected pair.
2d  log(median distance to nearest sender) vs log(sender density), against the
    homogeneous-Poisson line of slope -1/2 -- i.e. whether the model's
    independent variable is a measurement or a sender-calling rate.

REPAIRS, 2026-08-27
-------------------
* SILENT DEGRADATION REMOVED.  This file used to read perm_curves.csv and
  perm_nulls.csv behind `if os.path.exists(...) else pd.DataFrame()`, and
  fig2b/fig2c then skipped the null band and every null row without a word.  A
  missing input produced a COMPLETE-LOOKING figure with the null absent.  Every
  input is now required through _need(); if it is missing the script exits.
* The committed producer drew ONE band in 2b and TEN rows in 2c; the committed
  artefacts show THREE bands and NINETEEN rows.  The in-tissue variants are
  restored here (from perm_curves_c1.csv / perm_nulls_c1.csv), so the producer
  emits what the figures show.
* N3-var / N4-var, the FROZEN PRIMARY null, are added to 2c.  They CANNOT be
  added to 2b: perm_curves*.csv contains no `_var` draws at all, so no per-bin
  band for the primary null exists anywhere in the repository.  2b says so on
  its face rather than substituting a different null.
* 2b and 2d now write their backing CSVs.  Nothing in code/ wrote
  figure2b_data.csv or figure2d_data.csv before (CS_PHASE8_M1_RERUN.md:128
  asserts otherwise).
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
import summarize_phase3_c1 as SC1

PAL.apply_style(matplotlib)
RES = P.RESULTS
FIG = "/workspace/figures"


def _need(name, why):
    """Required input.  A missing null file used to yield a complete-looking
    figure with the null silently absent; it now stops the run."""
    p = f"{RES}/{name}"
    if not os.path.exists(p):
        sys.exit("make_figure2bc: %s is MISSING.\n  It supplies %s.\n"
                 "  Refusing to draw a figure that would look complete with "
                 "that missing." % (p, why))
    return pd.read_csv(p)


def lab(s):
    parts = s.split("_")
    return f"{parts[2]} {parts[-1].split('-')[0]} wk ({parts[0]})"


# --------------------------------------------------------------------- 2b ---
# (series key, source, null name in that source, colour, style, legend label)
BANDS = [
    ("null_N3", "perm_curves.csv", "N3", PAL.MUTED, dict(alpha=0.35, lw=0),
     "N3 torus shift, bounding box 95 %"),
    ("null_N3_snap", "perm_curves_c1.csv", "N3_snap", PAL.SERIES[2],
     dict(alpha=0.0, lw=1.1, ls=":"), "N3 shift, snapped into tissue 95 %"),
    ("null_N3_occ15", "perm_curves_c1.csv", "N3_occ15", PAL.SERIES[3],
     dict(alpha=0.0, lw=1.1, ls=":"), "N3 shift, occupancy-screened 95 %"),
]


def fig2b(modules=("tnfa_nfkb_proximal", "secondary_senescence"),
          celltype="Hepatocytes"):
    cur = _need("curves.csv", "the observed and matched-decoy curves")
    src = {"perm_curves.csv":
           _need("perm_curves.csv", "the published bounding-box N3 null band"),
           "perm_curves_c1.csv":
           _need("perm_curves_c1.csv", "the in-tissue (C1) N3 null bands")}
    for key, f, null, *_ in BANDS:
        if null not in set(src[f].null):
            sys.exit("make_figure2bc: %s carries no null=%r; figure 2b would "
                     "lose a band without saying so." % (f, null))

    order = [s for s in P.IN_BAND if s in set(cur.section)]
    rows = []
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
            for key, f, null, col, kw, leglab in BANDS:
                per = src[f]
                b = per[(per.section == sec) & (per.celltype == celltype)
                        & (per.module == mod) & (per.null == null)]
                if b.empty:
                    continue
                bx = 0.5 * (b.bin_lo + b.bin_hi)
                if kw["alpha"] > 0:
                    ax.fill_between(bx, b["lo"], b["hi"], color=col,
                                    label=leglab, **kw)
                else:
                    ax.plot(bx, b["lo"], color=col, lw=kw["lw"], ls=kw["ls"])
                    ax.plot(bx, b["hi"], color=col, lw=kw["lw"], ls=kw["ls"],
                            label=leglab)
                for _, rr in b.iterrows():
                    rows.append(dict(panel="2b", section=sec, celltype=celltype,
                                     module=mod, series=key, bin_lo=rr.bin_lo,
                                     bin_hi=rr.bin_hi, n=np.nan,
                                     value=rr["mean"], lo=rr["lo"],
                                     hi=rr["hi"], sem=np.nan))
            ax.plot(x, q.mean_decoy, color=PAL.SERIES[1], lw=1.6, ls="--",
                    label="matched decoy (N2)")
            ax.errorbar(x, q.mean_obs, yerr=q.sem_obs, color=PAL.SERIES[0],
                        lw=1.6, marker="o", ms=2.6, capsize=0,
                        label="observed senders")
            for _, rr in q.iterrows():
                rows.append(dict(panel="2b", section=sec, celltype=celltype,
                                 module=mod, series="observed",
                                 bin_lo=rr.bin_lo, bin_hi=rr.bin_hi, n=rr.n,
                                 value=rr.mean_obs, lo=np.nan, hi=np.nan,
                                 sem=rr.sem_obs))
                rows.append(dict(panel="2b", section=sec, celltype=celltype,
                                 module=mod, series="matched_decoy_N2",
                                 bin_lo=rr.bin_lo, bin_hi=rr.bin_hi,
                                 n=rr.n_decoy, value=rr.mean_decoy,
                                 lo=np.nan, hi=np.nan, sem=np.nan))
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
                 "the torus-shift null band,\nbounding-box and in-tissue "
                 "(correction C1)\n"
                 f"{celltype}, six Test-3-admissible sections, sender call "
                 f"{RN.PRIMARY_CALL}, 5 µm bins to {RN.WINDOW_UM:.0f} µm",
                 fontsize=9.5, y=1.01)
    fig.text(0.5, -0.035,
             "The frozen PRIMARY null (N3-var / N4-var, variance-corrected) "
             "has NO band here: perm_curves*.csv contains no variance-"
             "corrected draws, so no per-bin\ninterval for it exists in the "
             "repository. Its surviving fraction is in figure 2c "
             "(N3-var 0.996, N4-var 0.985).",
             fontsize=7.2, ha="center", color=PAL.STATUS["warning"],
             linespacing=1.5)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2b.{e}", bbox_inches="tight")
    plt.close(fig)
    D = pd.DataFrame(rows, columns=["panel", "section", "celltype", "module",
                                    "series", "bin_lo", "bin_hi", "n", "value",
                                    "lo", "hi", "sem"])
    D.to_csv(f"{RES}/figure2b_data.csv", index=False)
    print("wrote figure2b and figure2b_data.csv (%d rows, %d series)"
          % (len(D), D.series.nunique()))


# --------------------------------------------------------------------- 2c ---
BASE_NULLS = [
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
# the C1 in-tissue variants, from perm_nulls_c1.csv
C1_NULLS = [
    ("N3_tile", "N3-tile   shift inside solid-tissue tiles"),
    ("N3_occ", "N3-occ    ≤5 % out of tissue  (degenerate)"),
    ("N3_occ15", "N3-occ15  ≤15 % out of tissue"),
    ("N3_swap", "N3-swap   random real cell positions  (= N1)"),
    ("N3_snap", "N3-snap   shift, snapped into tissue"),
    ("N4_tile", "N4-tile   rotation inside solid-tissue tiles"),
    ("N4_occ", "N4-occ    ≤5 % out of tissue  (degenerate)"),
    ("N4_occ15", "N4-occ15  ≤15 % out of tissue"),
    ("N4_swap", "N4-swap   rotation, snapped into tissue"),
]
# the FROZEN PRIMARY null, from perm_nulls_var.csv
VAR_NULLS = [
    ("N3_var", "N3-var    variance-corrected shift   ★ PRIMARY"),
    ("N4_var", "N4-var    variance-corrected rotation ★ PRIMARY"),
]


def fig2c(call=RN.PRIMARY_CALL):
    mf = _need("main_fits.csv", "the reportable population and the N2/N5/N6 SFs")
    pf = _need("perm_nulls.csv", "the published N1/N3/N4 surviving fractions")
    c1 = _need("perm_nulls_c1.csv", "the in-tissue (C1) N3/N4 variants")
    va = _need("perm_nulls_var.csv",
               "the FROZEN PRIMARY variance-corrected N3-var/N4-var")
    n8p = [f"{RES}/n8_scrambled_{s}.csv" for s in P.IN_BAND]
    n8p = [p for p in n8p if os.path.exists(p)]
    if not n8p:
        sys.exit("make_figure2bc: no n8_scrambled_*.csv for the in-band "
                 "sections; the N8 row would vanish silently.")
    n8 = pd.concat([pd.read_csv(p) for p in n8p], ignore_index=True)

    key = ["section", "celltype", "module"]
    d = mf[mf.section.isin(P.IN_BAND) & (mf.call == call)
           & (mf.stratum == "all")].copy()
    d = d.merge(pf[key + ["N1_sf", "N3_sf", "N4_sf"]].drop_duplicates(key),
                on=key, how="left")
    d = d.merge(n8[key + ["sf_n8"]].drop_duplicates(key), on=key, how="left")
    rep = d[(d.beta_naive > 0) & (d.beta_base_lo > 0)]

    # the C1 and var runs are keyed on the same (section, celltype, module).
    c1f = c1[c1.scope == "full"] if "scope" in c1 else c1
    c1t = c1[c1.scope == "tile"] if "scope" in c1 else c1
    series = {}
    for col, _l in BASE_NULLS:
        if col in rep:
            series[col] = rep[col].to_numpy(float)
    for col, _l in C1_NULLS:
        src = c1t if col.endswith("_tile") else c1f
        m = rep[key].merge(src[key + [f"{col}_sf"]].drop_duplicates(key),
                           on=key, how="inner")
        series[col] = m[f"{col}_sf"].to_numpy(float)
    for col, _l in VAR_NULLS:
        m = rep[key].merge(va[key + [f"{col}_sf"]].drop_duplicates(key),
                           on=key, how="inner")
        series[col] = m[f"{col}_sf"].to_numpy(float)

    ALL = ([(c, l, "base") for c, l in BASE_NULLS]
           + [(c, l, "c1") for c, l in C1_NULLS]
           + [(c, l, "var") for c, l in VAR_NULLS])
    missing = [c for c, _l, _g in ALL if c not in series
               or not np.isfinite(series[c]).any()]
    if missing:
        sys.exit("make_figure2bc: no finite values for %s; figure 2c would "
                 "silently drop those rows." % ", ".join(missing))

    # rows are drawn top-to-bottom in ALL order, with a separator line before
    # each new group
    n_rows = len(ALL) + 2  # two separator slots
    ys, y = {}, n_rows - 1
    seps = []
    prev = None
    for col, _l, grp in ALL:
        if prev is not None and grp != prev:
            seps.append((y, grp))
            y -= 1
        ys[col] = y
        y -= 1
        prev = grp

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    rows, yticks, ylabels = [], [], []
    rng = np.random.default_rng(1)
    for col, l, grp in ALL:
        v = series[col]
        v = v[np.isfinite(v)]
        yy = ys[col]
        jit = (rng.random(v.size) - .5) * .36
        ax.scatter(np.clip(v, -1.15, 2.15), yy + jit, s=12, alpha=.4,
                   color=PAL.SERIES[0], lw=0)
        q = np.quantile(v, [.25, .5, .75])
        ax.plot(q[[0, 2]], [yy, yy], color=PAL.INK, lw=3.2,
                solid_capstyle="butt")
        ax.plot([q[1]], [yy], marker="|", ms=17,
                color=PAL.SERIES[5] if grp == "var" else PAL.SERIES[3], mew=3)
        rows.append(dict(null=l, group=grp, n=int(v.size), q25=q[0],
                         median=q[1], q75=q[2],
                         frac_le_0=float((v <= 0).mean()),
                         frac_gt_05=float((v > .5).mean())))
        yticks.append(yy)
        ylabels.append(f"{l}   (med {q[1]:+.2f})")
    for yy, grp in seps:
        ax.axhline(yy, color=PAL.AXIS, lw=1.0)
        yticks.append(yy)
        ylabels.append("— corrected, in-tissue (C1) —" if grp == "c1"
                       else "— FROZEN PRIMARY (variance-corrected) —")
    ax.axvline(1.0, color=PAL.MUTED, lw=1, ls=":")
    ax.axvline(0.0, color=PAL.STATUS["critical"], lw=1.2)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8.5)
    for t, lab_ in zip(ax.get_yticklabels(), ylabels):
        if "PRIMARY" in lab_:
            t.set_color(PAL.SERIES[5])
    ax.set_ylim(-0.8, n_rows - 0.2)
    ax.set_xlabel("surviving fraction of β̂   "
                  "(1 = the null removes nothing, 0 = it removes all of it)")
    ax.set_xlim(-1.2, 2.2)
    ax.set_title("Figure 2c — surviving fraction under each null, with the "
                 "corrected in-tissue N3/N4\nand the frozen primary "
                 "variance-corrected null\n"
                 f"six Test-3-admissible sections, {call}, receiver type × "
                 f"Tier B module; n = {len(rep)} of {len(d)} fits with a "
                 "positive, bootstrap-nonzero naive amplitude", fontsize=9.5)
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2c.{e}", bbox_inches="tight")
    plt.close(fig)
    R = pd.DataFrame(rows)
    R.to_csv(f"{RES}/figure2c_data.csv", index=False)
    print(R.round(3).to_string(index=False))
    print("wrote figure2c and figure2c_data.csv (%d rows)" % len(R))


# --------------------------------------------------------------------- 2d ---

def fig2d():
    d = _need("poisson_density.csv", "the 77 section x sender-definition points")
    f = _need("poisson_fits.csv", "the fitted slope and r2 in the title")
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
    plt.close(fig)
    obs = d.copy()
    obs["panel"] = "2d"
    obs["poisson_prediction_um"] = np.sqrt(
        np.log(2) / (np.pi * obs.sender_density_per_um2))
    obs["series"] = "observed"
    fit = f.copy()
    fit["panel"] = "2d"
    fit["series"] = "fitted_slope"
    pd.concat([obs, fit], ignore_index=True).to_csv(
        f"{RES}/figure2d_data.csv", index=False)
    print("wrote figure2d and figure2d_data.csv")


if __name__ == "__main__":
    fig2b()
    fig2c()
    fig2d()
