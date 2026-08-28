#!/usr/bin/env python3
"""
Phase 5 figures.

Figure 2a (REGENERATED, Section 25) — the naive decay curve, now stratified by
receiver cell type on the corrected six-section admissible set.  The Phase 2
version was fitted unstratified on provisional module scores and CS Phase 3 §3.1
showed 66 % of what it displayed is receiver cell-type composition, so it cannot
go in the paper as it stands.  The new version puts the composition artefact and
its correction on the same axes.

Figure 3 (Section 25) — controlled kernel estimates:
  (a) lambda-hat with donor-bootstrap CIs per module x receiver cell type
  (b) kernel family comparison, naive vs under full control
  (c) lambda_proximal vs lambda_downstream (Section 6.4)
  (d) nearest-sender vs superposition (Section 6.3), against the synthetic
      calibration and the torus-shift / label-permutation nulls
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, "/workspace/code")
import sasp_palette as PAL
import sasp_phase3 as P
import run_phase3_nulls as RN
import phase5_common as C5

PAL.apply_style(matplotlib)
FIG = "/workspace/figures"
RES5 = C5.RES5
# non-uniform bins: the distance distribution is concentrated below ~30 um
# (CS Phase 3 §2: median 32.6 um, p99 98.3 um), so uniform 5 um bins spend
# three quarters of the axis on the empty tail and hide the only region with
# resolution.
BINS = np.array([0., 2., 4., 6., 8., 10., 12.5, 15., 17.5, 20., 25., 30.,
                 35., 40., 50., 60., 75., 100.])
MID = 0.5 * (BINS[:-1] + BINS[1:])
# Tier B numbering is build_genesets.py:105-118 (TIER_B_SRC + B7_CUR), which
# Master Plan §9 lines 455-461 and every report follow:
#   B1 tnfa_nfkb_proximal  B2 il6_jak_stat3  B3 interferon_response
#   B4 downstream_arrest   B5 emt_ecm        B6 oxidative_stress
#   B7 secondary_senescence
# Until 2026-08-27 this dict labelled secondary_senescence "B5", emt_ecm "B6"
# and oxidative_stress "B7" -- a three-way rotation, and figure 2a's legend was
# the ONLY place in the repository using that numbering.  The consequence was
# direct: the module a reader identified as "B7" was oxidative stress (31
# genes), while the re-sourced B7 (secondary_senescence, 108 genes -- the
# frozen configuration's headline gene-set rebuild) appeared as "B5".
SHORT = {"tnfa_nfkb_proximal": "TNFA/NF-kB (B1)",
         "il6_jak_stat3": "IL6/JAK/STAT3 (B2)",
         "interferon_response": "Interferon (B3)",
         "downstream_arrest": "Arrest (B4)",
         "emt_ecm": "EMT/ECM (B5)",
         "oxidative_stress": "Oxidative stress (B6)",
         "secondary_senescence": "2nd senescence (B7)"}
# plotting order == the B1..B7 numbering, so the legend reads in order
ORDER = ["tnfa_nfkb_proximal", "il6_jak_stat3", "interferon_response",
         "downstream_arrest", "emt_ecm", "oxidative_stress",
         "secondary_senescence"]


# ---------------------------------------------------------------------------
# Figure 2a data
# ---------------------------------------------------------------------------

def _ts(t):
    import time
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(t))


def _stale_inputs(cache):
    """Inputs newer than the figure2a cache, newest first.

    figure2a caches its own binned input in figures/, and that cache is TRACKED and is one
    of the 52 artefacts check_figures_guard.py watches.  So `git clone` + this script
    regenerates figure2a byte-identically from the committed cache no matter what happened
    upstream, and the guard then reports OK -- its pass on figure2a is vacuous
    (AUDIT_REPRODUCIBILITY D1).  The only mitigation used to be one `rm -f` line in
    _m1_rerun_stage5.sh.  Compare mtimes instead, and refuse rather than warn.
    """
    import glob
    if not os.path.exists(cache):
        return []
    t0 = os.path.getmtime(cache)
    ins = (glob.glob(f"{P.CACHE3}/*.npz") + glob.glob("/workspace/genesets/*.txt")
           + glob.glob("/workspace/data/processed/senders_*.csv"))
    out = [(os.path.getmtime(x), x) for x in ins if os.path.getmtime(x) > t0]
    return sorted(out, reverse=True)


def check_cache_fresh(cache, accept_stale=False):
    """Refuse to plot from a cache that is older than its inputs."""
    stale = _stale_inputs(cache)
    if stale and not accept_stale:
        raise SystemExit(
            "figure2a: the cached %s is OLDER than %d of its inputs, so the figure would be\n"
            "drawn from superseded data -- and because the cache is itself a tracked,\n"
            "guard-watched artefact, check_figures_guard.py would still report OK.\n"
            "  newest stale input: %s (%s)\n"
            "  cache mtime:        %s\n"
            "Pass --rebuild to recompute it, or --accept-stale if you really mean it."
            % (cache, len(stale), stale[0][1], _ts(stale[0][0]), _ts(os.path.getmtime(cache))))
    return True


def build_fig2a_data(sections, call=RN.PRIMARY_CALL):
    """Binned response vs distance, three ways:

      unstrat  : z-scored within SECTION only -- what the field plots
      composition-only : each cell's response replaced by the mean response of
                 its cell type in that section, then binned the same way; this
                 curve is what receiver cell-type composition alone produces
      stratified : z-scored within SECTION x CELL TYPE, so composition is
                 removed by construction
    """
    acc = {}

    def add(key, ix, y, d):
        s = acc.setdefault(key, [np.zeros(BINS.size - 1),
                                 np.zeros(BINS.size - 1),
                                 np.zeros(BINS.size - 1)])
        cnt = np.bincount(ix, minlength=BINS.size - 1)
        tot = np.bincount(ix, weights=y, minlength=BINS.size - 1)
        tot2 = np.bincount(ix, weights=y * y, minlength=BINS.size - 1)
        s[0] += cnt
        s[1] += tot
        s[2] += tot2

    for sample in sections:
        sec = P.Sec(sample)
        sender = sec.sender_mask(call)
        coords = sec.coords.astype(float)
        d = P.dist_to_senders(coords, sender)
        keep = ((~np.isin(sec.celltype, P.EXCLUDE_TYPES)) & (~sender)
                & np.isfinite(d) & (d <= RN.WINDOW_UM))
        ii = np.flatnonzero(keep)
        dd = d[ii]
        ix = np.clip(np.digitize(dd, BINS) - 1, 0, BINS.size - 2)
        ct = sec.celltype[ii]
        for j, mod in enumerate(P.MODULES):
            y = sec.module(mod)[ii]
            yz = (y - y.mean()) / (y.std() + 1e-12)
            add(("unstrat", mod, "all"), ix, yz, dd)
            # composition-only surrogate: cell-type mean within this section
            comp = np.empty_like(yz)
            ys = np.empty_like(yz)
            for c in np.unique(ct):
                m = ct == c
                comp[m] = yz[m].mean()
                ys[m] = (y[m] - y[m].mean()) / (y[m].std() + 1e-12)
            add(("compositiononly", mod, "all"), ix, comp, dd)
            add(("strat", mod, "all"), ix, ys, dd)
            for c in np.unique(ct):
                m = ct == c
                if m.sum() < RN.MIN_RECEIVERS:
                    continue
                add(("strat", mod, c), ix[m], ys[m], dd[m])
                add(("unstrat", mod, c), ix[m], yz[m], dd[m])
        print(f"[fig2a] {sample}", flush=True)

    rows = []
    for (kind, mod, ct), (cnt, tot, tot2) in acc.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            mu = tot / np.where(cnt > 0, cnt, np.nan)
            var = tot2 / np.where(cnt > 0, cnt, np.nan) - mu ** 2
            sem = np.sqrt(np.maximum(var, 0)) / np.sqrt(np.maximum(cnt, 1))
        for b in range(BINS.size - 1):
            rows.append(dict(kind=kind, module=mod, celltype=ct,
                             bin_lo=BINS[b], bin_hi=BINS[b + 1],
                             n=cnt[b], mean=mu[b], sem=sem[b]))
    return pd.DataFrame(rows)


def _amp(q):
    """Contact amplitude: the first bin retained (>= 30 cells) minus the
    40-100 um plateau, in within-stratum sd units."""
    b = q[q.bin_lo >= 40.0]
    if q.empty or b.empty:
        return np.nan
    w = b.n.to_numpy(float)
    return float(q["mean"].iloc[0] - np.average(b["mean"], weights=w))


def fig2a(df):
    cts = [c for c in sorted(set(df[df.kind == "strat"].celltype)) if c != "all"]
    panels = ["ALL RECEIVERS"] + cts
    ncol = 4
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.15 * ncol, 2.7 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    amps = {}
    for k, pan in enumerate(panels):
        ax = axes[k // ncol, k % ncol]
        kind = "unstrat" if pan == "ALL RECEIVERS" else "strat"
        sel = "all" if pan == "ALL RECEIVERS" else pan
        aa = []
        ntot = 0
        for m, mod in enumerate(ORDER):
            q = df[(df.kind == kind) & (df.module == mod)
                   & (df.celltype == sel)].sort_values("bin_lo")
            q = q[q.n >= 30]
            if q.empty:
                continue
            ntot = max(ntot, q.n.sum())
            x = 0.5 * (q.bin_lo + q.bin_hi)
            ax.plot(x, q["mean"], color=PAL.SERIES[m % 8], lw=1.5,
                    label=SHORT[mod])
            aa.append(_amp(q))
        amps[pan] = float(np.nanmedian(aa)) if aa else np.nan
        # how much of the panel sits in the first bin -- the "contact spike"
        q0 = df[(df.kind == kind) & (df.celltype == sel)
                & (df.module == ORDER[0])].sort_values("bin_lo")
        f0 = (float(q0[q0.bin_hi <= 10.0].n.sum() / max(q0.n.sum(), 1))
              if len(q0) else np.nan)
        ax.axvspan(0, 10, color=PAL.STATUS["warning"], alpha=0.12, lw=0)
        if pan == "ALL RECEIVERS":
            for m, mod in enumerate(ORDER):
                q = df[(df.kind == "compositiononly") & (df.module == mod)
                       & (df.celltype == "all")].sort_values("bin_lo")
                q = q[q.n >= 30]
                x = 0.5 * (q.bin_lo + q.bin_hi)
                ax.plot(x, q["mean"], color=PAL.SERIES[m % 8], lw=1.2, ls=":")
            ax.set_title("UNSTRATIFIED — what the field plots\n"
                         "dotted: receiver cell-type composition alone",
                         fontsize=8.5, color=PAL.INK)
        else:
            ax.set_title(f"{pan}   n={ntot:,.0f}", fontsize=8.5, color=PAL.INK)
        ax.text(0.97, 0.05,
                f"contact amplitude {amps[pan]:+.2f} sd\n"
                f"only {100*f0:.1f}% of receivers are < 10 um",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8,
                color=PAL.INK2)
        ax.axhline(0, color=PAL.AXIS, lw=0.7)
        ax.set_xscale("function", functions=(np.sqrt, np.square))
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 5, 10, 20, 40, 60, 100])
        ax.set_ylim(-0.75, 0.65)
    for k in range(len(panels), nrow * ncol):
        axes[k // ncol, k % ncol].set_axis_off()
    for r in range(nrow):
        axes[r, 0].set_ylabel("response (sd)")
    for c in range(ncol):
        axes[min(nrow - 1, (len(panels) - 1) // ncol), c].set_xlabel(
            "distance to nearest sender (um), sqrt scale")
    h = [Line2D([], [], color=PAL.SERIES[m % 8], lw=1.8, label=SHORT[mod])
         for m, mod in enumerate(ORDER)]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.012))
    fig.suptitle(
        "Figure 2a  Binned response vs distance to nearest sender, by "
        "receiver cell type\n"
        "six Test-3-admissible sections, six mice, sender call tierA_p95 "
        "(union-strict Tier A, within-cell-type p95); responses z-scored "
        "within section (unstratified panel) or within section x cell type "
        "(all others)", fontsize=10, y=1.008)
    fig.tight_layout(rect=(0, 0.04, 1, 0.985))
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2a.{ext}", bbox_inches="tight")
    plt.close(fig)
    pd.Series(amps, name="contact_amplitude_sd").to_csv(
        f"{FIG}/figure2a_amplitudes.csv")
    print("wrote figure2a")
    print(pd.Series(amps).round(3).to_string())


# ---------------------------------------------------------------------------
# Figure 3
# ---------------------------------------------------------------------------

def fig3():
    # PROVENANCE.  Nothing in code/ wrote figures/figure3_data.csv until
    # 2026-08-27 (CS_PHASE8_M1_RERUN.md:129,230 asserts otherwise), and the
    # committed CSV covered only the `call = tierA_p95` sources, so 6 of panel
    # (d)'s 12 real boxplots -- the Cdkn1a+ and SenePy p95 caller rows -- had no
    # backing rows at all, and panel (a) carried 63 rows for the 42 CIs it
    # draws (the 21 zonation rows are not in the figure).  `parts` collects
    # exactly the frames the panels draw, and fig3 writes them at the end.
    parts = []

    def keep(panel, source, frame):
        f = frame.copy()
        f.insert(0, "source", source)
        f.insert(0, "panel", panel)
        parts.append(f)
        return frame

    fig = plt.figure(figsize=(12.4, 10.2))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.40)

    # (a) lambda-hat with donor-bootstrap CIs -------------------------------
    ax = fig.add_subplot(gs[0, 0])
    cd = pd.read_csv("/workspace/results/phase3/combined_donor.csv")
    cd = keep("3a", "combined_donor.csv", cd[cd.zone == "all"].copy())
    cts = sorted(set(cd.celltype))
    y = 0
    ticks, labs = [], []
    span = 0
    for ct in cts:
        for mod in ORDER:
            q = cd[(cd.celltype == ct) & (cd.module == mod)]
            if q.empty:
                continue
            r = q.iloc[0]
            lo = r.get("lam_full_donor_lo", np.nan)
            hi = r.get("lam_full_donor_hi", np.nan)
            full = (np.isfinite(lo) and np.isfinite(hi)
                    and lo <= 7.3 and hi >= 49.0)
            span += int(full)
            ax.plot([lo, hi], [y, y], color=(PAL.STATUS["critical"] if full
                                             else PAL.MUTED), lw=2.2,
                    solid_capstyle="butt", alpha=0.85)
            ax.plot([r.lam_full_profiled], [y], "o", ms=3.4,
                    color=PAL.SERIES[ORDER.index(mod) % 8])
            y += 1
        ticks.append(y - 3.5)
        labs.append(ct)
        y += 1.5
    ax.axvline(7, color=PAL.INK2, ls="--", lw=0.9)
    ax.axvline(50, color=PAL.INK2, ls="--", lw=0.9)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labs, fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlim(6, 105)
    ax.set_xlabel("lambda-hat (um), donor bootstrap over 6 animals")
    ax.legend(handles=[
        Line2D([], [], color=PAL.STATUS["critical"], lw=2.4,
               label="CI spans the whole grid"),
        Line2D([], [], color=PAL.MUTED, lw=2.4, label="CI narrower than grid")],
        loc="center right", fontsize=7.5)
    ax.set_title(f"(a) No length constant is identified\n"
                 f"{span} of {len(cd)} donor-bootstrap CIs on lambda under "
                 f"N5+N6 span all of [7, 50] um", fontsize=9.5, color=PAL.INK)

    # (b) kernel family comparison -----------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    kf = keep("3b", "kernel_families.csv",
              pd.read_csv(f"{RES5}/kernel_families.csv"))
    fams = ["exponential", "gaussian", "powerlaw", "step", "spline"]
    w = 0.36
    for i, (des, col, lab) in enumerate(
            [("naive", PAL.SERIES[1], "naive"),
             ("ctrl", PAL.SERIES[0], "under N5+N6 control")]):
        vals = [kf[(kf.design == des) & (kf.family == f)].d_aic_vs_cov
                for f in fams]
        xs = np.arange(len(fams)) + (i - 0.5) * w
        med = [v.median() for v in vals]
        lo = [v.quantile(.25) for v in vals]
        hi = [v.quantile(.75) for v in vals]
        ax.bar(xs, med, width=w, color=col, label=lab, alpha=0.9)
        ax.vlines(xs, lo, hi, color=PAL.INK2, lw=1.0)
    ax.axhline(0, color=PAL.INK, lw=1.0)
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels(fams, fontsize=8, rotation=12)
    ax.set_ylabel("dAIC vs covariates-only  (<0 = kernel helps)")
    ax.legend(loc="upper left")
    wins = kf[kf.design == "ctrl"].groupby("family").apply(
        lambda g: (g.d_aic_vs_best == 0).mean())
    ax.set_title("(b) Under control no family earns its place\n"
                 "step wins %.0f%% of fits by AIC and still only beats "
                 "no kernel in %.0f%%"
                 % (100 * wins.get("step", np.nan),
                    100 * (kf[(kf.design == "ctrl") & (kf.family == "step")]
                           .d_aic_vs_cov < 0).mean()),
                 fontsize=9.5, color=PAL.INK)

    # (c) proximal vs downstream -------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    pd_ = keep("3c", "proximal_vs_downstream.csv",
               pd.read_csv(f"{RES5}/proximal_vs_downstream.csv"))
    yv = np.arange(len(pd_))
    ax.hlines(yv, pd_.ratio_lo, pd_.ratio_hi, color=PAL.MUTED, lw=2.4)
    ax.plot(pd_.ratio, yv, "o", color=PAL.SERIES[0], ms=6)
    ax.axvline(1.0, color=PAL.INK, lw=1.0)
    ax.axvline(50 / 7, color=PAL.INK2, ls="--", lw=0.9)
    ax.axvline(7 / 50, color=PAL.INK2, ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_xticks([0.14, 0.25, 0.5, 1.0, 2.0, 4.0, 7.14])
    ax.set_xticklabels(["0.14", "0.25", "0.5", "1", "2", "4", "7.1"],
                       fontsize=8)
    ax.minorticks_off()
    ax.set_yticks(yv)
    ax.set_yticklabels(pd_.celltype, fontsize=8)
    ax.set_xlabel("lambda_proximal / lambda_downstream  (B1 / B4)")
    ax.set_title("(c) Section 6.4 is not estimable here\n"
                 "dashed = the widest ratio the [7, 50] um grid allows; "
                 "%d/%d CIs reach it"
                 % (int(((pd_.ratio_lo <= 0.15) & (pd_.ratio_hi >= 7.0)).sum()),
                    len(pd_)), fontsize=9.5, color=PAL.INK)

    # (d) superposition vs nearest -----------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    # CS_PHASE5 Sec 10.4 item 5: the panel carries all three near-independent
    # sender callers, not tierA_p95 alone -- model selection between kernel
    # families is sender-definition dependent and the panel must show it.
    CALLERS = [("tierA_p95", ""),
               ("Cdkn1a+", "_cdkn1a_pos"),
               ("SenePy p95", "_senepy_p95")]
    ss = keep("3d", "super_section.csv",
              pd.read_csv(f"{RES5}/super_section.csv"))
    sn = keep("3d", "super_nulls.csv",
              pd.read_csv(f"{RES5}/super_nulls.csv"))
    ss["v"] = 1000 * ss.d_aic / ss.n
    sn["v"] = 1000 * sn.d_aic / sn.n
    groups = [("planted superposition\n(synthetic, confounded)", None,
               PAL.STATUS["good"]),
              ("planted nearest\n(synthetic, confounded)", None,
               PAL.STATUS["good"]),
              ("REAL naive \u00b7 tierA_p95", ss[ss.design == "naive"].v,
               PAL.SERIES[1])]
    for cname, suf in CALLERS:
        cs = pd.read_csv(f"{RES5}/super_section{suf}.csv")
        cn = pd.read_csv(f"{RES5}/super_nulls{suf}.csv")
        if suf:  # the tierA_p95 pair is already in `parts`, unsuffixed
            keep("3d", f"super_section{suf}.csv", cs)
            keep("3d", f"super_nulls{suf}.csv", cn)
        cs["v"] = 1000 * cs.d_aic / cs.n
        cn["v"] = 1000 * cn.d_aic / cn.n
        groups += [
            (f"REAL N5+N6 ctrl \u00b7 {cname}", cs[cs.design == "ctrl"].v,
             PAL.SERIES[0]),
            (f"N1 permuted \u00b7 {cname}", cn[cn.null == "N1"].v,
             PAL.MUTED),
            (f"N3 torus-shift \u00b7 {cname}", cn[cn.null == "N3"].v,
             PAL.MUTED)]
    syn = pd.read_csv("/workspace/results/phase1b/misspec.csv")
    syn = keep("3d", "misspec.csv", syn[syn.fit_family == "exponential"])
    pv = syn.pivot_table(index=["true_superposition", "regime", "rep"],
                         columns="fit_mode", values=["aic", "n"])
    sv = ((pv["aic"]["superposition"] - pv["aic"]["nearest"])
          / pv["n"]["nearest"] * 1000).reset_index()
    sv.columns = ["true_sup", "regime", "rep", "v"]
    data = []
    for i, (lab, v, col) in enumerate(groups):
        if v is None:
            v = sv[(sv.true_sup == (i == 0)) & (sv.regime == "confounded")].v
        data.append(np.asarray(v, float))
    pos = np.arange(len(groups))
    bp = ax.boxplot([d[np.isfinite(d)] for d in data], positions=pos,
                    vert=False, widths=0.6, showfliers=False,
                    patch_artist=True)
    for patch, (lab, v, col) in zip(bp["boxes"], groups):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
        patch.set_edgecolor(PAL.INK2)
    for el in ("medians", "whiskers", "caps"):
        for a in bp[el]:
            a.set_color(PAL.INK)
    ax.axvline(0, color=PAL.INK, lw=1.0)
    ax.set_yticks(pos)
    ax.set_yticklabels([g[0] for g in groups], fontsize=7.2)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlabel("AIC(superposition) - AIC(nearest), per 1,000 cells\n"
                  "<0 favours superposition (dose), >0 favours nearest "
                  "(threshold)")
    ax.set_title("(d) Section 6.3 is answerable in synthetic tissue and not on "
                 "this data:\nthe real spread is ~3 orders of magnitude smaller, "
                 "at every sender caller", fontsize=9.5, color=PAL.INK)

    fig.suptitle("Figure 3  Controlled kernel estimates — nothing is "
                 "identified, and each panel says so with a bound",
                 fontsize=11.5, y=0.985)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure3.{ext}", bbox_inches="tight")
    plt.close(fig)
    D = pd.concat(parts, ignore_index=True)
    D.to_csv(f"{FIG}/figure3_data.csv", index=False)
    print("wrote figure3 and figure3_data.csv")
    print(D.groupby(["panel", "source"]).size().to_string())


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="2a,3")
    ap.add_argument("--rebuild", action="store_true",
                    help="recompute figure2a_stratified_curves.csv instead of "
                         "reusing the cached copy")
    ap.add_argument("--accept-stale", action="store_true",
                    help="draw figure2a from the cache even when an input is newer "
                         "than it (default: refuse)")
    a = ap.parse_args()
    w = a.which.split(",")
    if "2a" in w:
        f = f"{FIG}/figure2a_stratified_curves.csv"
        if os.path.exists(f) and not a.rebuild:
            check_cache_fresh(f, a.accept_stale)
            print("figure2a: REUSING the cached %s (mtime %s). It is NOT "
                  "rebuilt from the senders/genesets on disk -- pass --rebuild "
                  "after any change to either." % (f, _ts(os.path.getmtime(f))))
            df = pd.read_csv(f)
        else:
            df = build_fig2a_data(P.IN_BAND)
            df.to_csv(f, index=False)
        fig2a(df)
    if "3" in w:
        fig3()
