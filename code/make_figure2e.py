#!/usr/bin/env python3
"""Figure 2e-h — the coordinate nulls audited.

WHY THIS FILE EXISTS.  Until 2026-08-27 figure2e.png/.pdf and
results/phase3/figure2e_data.csv were in the repository with NO PRODUCER
ANYWHERE.  `grep -rl figure2e --include=*.py` hit exactly one file, the figure
guard.  reports/CS_PHASE8_C1_CLOSEOUT.md:305 cites "the docstring of fig2e()";
no such function existed.  The figure could not be regenerated, re-audited or
corrected, and it had drifted: panel (f) drew a black reference line labelled
"median lambda-hat 15.7 um" -- the interior-only median that
reports/COMPLETED_TASKS.md:165-169 and CS_PHASE8_TORUS_VAR.md:7-9 record as
WITHDRAWN, unsourced and circular -- while the panel's written caption
(CS_PHASE8_C1_CLOSEOUT.md §4.1-4.2) called the same line "12.8 um".  One drawn
line, three numbers.

This producer rebuilds figure2e_data.csv from upstream results and redraws the
figure from that CSV.  Every row of the committed CSV was verified to
reproduce exactly from the sources below before this file was written, so the
reconstruction is faithful to what made the artefact; the DELIBERATE changes
are listed under CHANGES.

SOURCES (nothing is typed in):
  results/phase3/sf_summary_var.csv       per-variant SF median/IQR + geometry
  results/phase3/sf_summary.csv           the PUBLISHED N3/N4 SF (panel g diamond)
  results/phase3/null_destructiveness.csv frac_in_occupancy (bounding-box nulls)
  results/phase3/null_destructiveness_var.csv  ... for the variance-corrected pair
  results/phase3/main_fits.csv            lambda-hat: median, IQR, railing rate
  results/phase3/a7_summary.csv           panel h amplitudes (section-clustered)
  results/phase3/a7_control_probe_fits.csv    panel h SCOPE (sections, cell types)
  results/phase3/a7_control_probe_curves.csv  panel h curve

CHANGES vs the 2026-08-20 artefact this replaces:
  1. Panel f's reference line is lambda-hat = 14.7321 um -- the median lam_naive
     over the 315 primary in-band tierA_p95 stratum=all fits, which
     summary_phase3.txt §6 prints as `medlam 14.7` -- NOT 15.6821, and it is
     drawn WITH its IQR and railing rate, because 60 % of those fits are pinned
     to a grid bound and a bare point estimate misrepresents that distribution.
  2. N3-var / N4-var, the FROZEN PRIMARY null (variance-corrected shift,
     Mrkvicka et al. 2021), are drawn in every panel that reports a null.  They
     were the only 2 of the 15 variants in sf_summary_var.csv that no figure in
     the repository plotted.
  3. Panel h's scope box no longer applies one clause to two populations.  The
     CURVE is six in-band sections, hepatocytes; the beta/sd amplitudes beside
     it are 11 sections x 9 cell types (a7_control_probe_fits.csv), which the
     old box called "six in-band sections, hepatocytes".  Both scopes are now
     read from the files and printed separately.
  4. The per-control-family amplitudes the box's "40 probes + 609 codewords +
     21 genomic controls" implies are now in figure2e_data.csv and summarised
     in the box, instead of only the pooled all_controls value.
"""
from __future__ import annotations

import os
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

PAL.apply_style(matplotlib)
RES = P.RESULTS
FIG = "/workspace/figures"

# family -> colour/label.  "primary" is the frozen configuration's null.
FAM_COL = {"published bounding box": PAL.SERIES[7],
           "in tissue": PAL.SERIES[0],
           "degenerate": PAL.SERIES[3],
           "equals N1": PAL.SERIES[1],
           "PRIMARY (variance-corrected)": PAL.SERIES[5]}

# variant -> (row label, family, which sf_summary_var variant supplies the SF,
#             which sf_summary.csv row supplies the PUBLISHED value, note)
ROWS = [
    ("N3_orig", "N3-orig", "published bounding box", "N3_orig_rerun",
     "N3 torus shift", "PUBLISHED, leaves the tissue"),
    ("N3_tile", "N3-tile", "in tissue", "N3_tile", None,
     "in-tissue, clustering-preserving"),
    ("N3_occ", "N3-occ", "degenerate", "N3_occ", None,
     "degenerate: displacement below the fitting window"),
    ("N3_occ15", "N3-occ15", "in tissue", "N3_occ15", None,
     "in-tissue, clustering-preserving"),
    ("N3_swap", "N3-swap", "equals N1", "N3_swap", None,
     "equals N1 (label permutation), not a torus shift"),
    ("N3_snap", "N3-snap", "in tissue", "N3_snap", None,
     "in-tissue, clustering-preserving"),
    ("N3_var", "N3-var", "PRIMARY (variance-corrected)", "N3_var", None,
     "FROZEN PRIMARY null: variance-corrected shift, tissue window W"),
    ("N4_orig", "N4-orig", "published bounding box", "N4_orig_rerun",
     "N4 rotation", "PUBLISHED, leaves the tissue"),
    ("N4_tile", "N4-tile", "in tissue", "N4_tile", None,
     "in-tissue, clustering-preserving"),
    ("N4_occ", "N4-occ", "degenerate", "N4_occ", None,
     "degenerate: displacement below the fitting window"),
    ("N4_occ15", "N4-occ15", "in tissue", "N4_occ15", None,
     "in-tissue, clustering-preserving"),
    ("N4_swap", "N4-swap", "in tissue", "N4_swap", None,
     "in-tissue, clustering-preserving"),
    ("N4_var", "N4-var", "PRIMARY (variance-corrected)", "N4_var", None,
     "FROZEN PRIMARY null: variance-corrected rotation, tissue window W"),
]
# separators drawn between the N3 block, the N4 block and the primary pair
A7_DESIGNS = [("base", "naive (intercept only)"),
              ("n6", "+N6 neighbour baseline"),
              ("n5", "+N5 technical covariates"),
              ("n6n5", "+N6+N5 (full nuisance design)"),
              ("n2", "N2 matched-decoy contrast")]
# Xenium 5K panel design constants, documented at run_a7_control_probes.py:10.
PANEL_CONTROLS = "40 probes + 609 codewords + 21 genomic controls"


def _require(path):
    if not os.path.exists(path):
        raise SystemExit("figure2e: %s is missing. This figure reports a null "
                         "battery; a missing input would silently remove a "
                         "null and leave a complete-looking figure. Refusing "
                         "to draw." % path)
    return path


# ---------------------------------------------------------------- data ------

def lambda_hat():
    """The primary fits' lambda-hat: median, IQR, railing rate, n.

    summarize_phase3.py's §1/§6 population verbatim: in-band sections, the
    primary sender call, stratum=all -- ALL 315 of them, railed included.
    Dropping the railed fits is the interior-only definition that produced the
    withdrawn 15.68.
    """
    mf = pd.read_csv(_require(f"{RES}/main_fits.csv"))
    d = mf[mf.section.isin(P.IN_BAND) & (mf.call == RN.PRIMARY_CALL)
           & (mf.stratum == "all")]
    v = d.lam_naive.to_numpy(float)
    q25, med, q75 = np.quantile(v, [.25, .5, .75])
    return dict(median=float(med), q25=float(q25), q75=float(q75),
                railed=float(d.lam_railed.mean()), n=int(len(d)))


def build_data():
    sv = pd.read_csv(_require(f"{RES}/sf_summary_var.csv"))
    sv = sv[sv.subset.str.startswith("PRIMARY")].set_index("variant")
    ss = pd.read_csv(_require(f"{RES}/sf_summary.csv"))
    ss = ss[ss.subset.str.startswith("PRIMARY")].set_index("null")
    nd = pd.read_csv(_require(f"{RES}/null_destructiveness.csv"))
    ndv = pd.read_csv(_require(f"{RES}/null_destructiveness_var.csv"))
    occ = pd.concat([nd, ndv], ignore_index=True)
    occ = occ[occ.section.isin(P.IN_BAND) & (occ.call == RN.PRIMARY_CALL)]
    occ = occ.groupby("null").frac_in_occupancy.median()

    rows = []
    for var, _lab, fam, sfsrc, pubrow, note in ROWS:
        if var not in sv.index:
            raise SystemExit("figure2e: variant %s absent from "
                             "sf_summary_var.csv" % var)
        g, s = sv.loc[var], sv.loc[sfsrc]
        rows.append(dict(
            null=var, family=fam,
            frac_retaining_a_neighbour=g.frac_retaining_a_neighbour,
            real_median_nbrs=g.real_median_nbrs,
            null_median_nbrs=g.null_median_nbrs,
            frac_in_occupancy=occ.get(var, np.nan),
            median_displacement_um=g.median_displacement_um,
            sf_median=s["median"], sf_q25=s.q25, sf_q75=s.q75,
            sf_published=(ss.loc[pubrow, "median"] if pubrow else np.nan),
            note=note))

    lam = lambda_hat()
    rows.append(dict(null="reference_median_lambda_hat_um", family="",
                     median_displacement_um=lam["median"],
                     sf_q25=lam["q25"], sf_q75=lam["q75"],
                     note="panel f line: median lam_naive over the %d primary "
                          "in-band %s stratum=all fits; IQR [%.1f, %.1f]; "
                          "%.0f %% railed at a grid bound. NOT the "
                          "interior-only median (15.68, withdrawn)."
                          % (lam["n"], RN.PRIMARY_CALL, lam["q25"], lam["q75"],
                             100 * lam["railed"])))
    rows.append(dict(null="reference_lambda_hat_frac_railed", family="",
                     sf_median=lam["railed"], note="panel f line annotation"))
    rows.append(dict(null="reference_N1_sf", family="",
                     sf_median=ss.loc["N1 stratified label permutation",
                                      "median"],
                     note="panel g line"))

    # ---- panel h ---------------------------------------------------------
    a7 = pd.read_csv(_require(f"{RES}/a7_summary.csv")).set_index(
        ["response", "design"])
    fits = pd.read_csv(_require(f"{RES}/a7_control_probe_fits.csv"))
    scope = dict(n_sections=int(fits.section.nunique()),
                 n_celltypes=int(fits.celltype.nunique()))
    for resp, tag in (("all_controls", "A7_all_controls"),
                      ("neg_control_probe", "A7_neg_control_probe"),
                      ("neg_control_codeword", "A7_neg_control_codeword"),
                      ("genomic_control", "A7_genomic_control"),
                      ("BIOLOGICAL MODULES (reference)", "A7_BIOLOGICAL")):
        for d, _dl in A7_DESIGNS:
            if (resp, d) not in a7.index:
                continue
            r = a7.loc[(resp, d)]
            rows.append(dict(
                null=f"{tag}_{d}", family="M1 negative-control-probe kernel",
                sf_median=r.clustered_mean, sf_q25=r.clustered_lo,
                sf_q75=r.clustered_hi,
                note="panel h: beta/sd_y, section-clustered mean [95%% CI], "
                     "p=%g, %d fits over %d sections x %d cell types"
                     % (r.clustered_p, int(r.n_fits), int(r.n_sections),
                        scope["n_celltypes"])))

    cur = pd.read_csv(_require(f"{RES}/a7_control_probe_curves.csv"))
    cur = cur[(cur.response == "all_controls") & cur.section.isin(P.IN_BAND)
              & (cur.call == RN.PRIMARY_CALL)]
    scope["curve_celltypes"] = sorted(cur.celltype.unique())
    scope["curve_sections"] = int(cur.section.nunique())
    for b, s in cur.groupby("bin_mid"):
        rows.append(dict(
            null="A7_curve_bin_%gum" % b,
            family="M1 negative-control-probe kernel",
            median_displacement_um=b,
            sf_median=float(np.average(s.mean_raw, weights=s.n)),
            sf_published=float(np.average(s.mean_resid, weights=s.n)),
            note="panel h curve: sf_median=naive mean z, "
                 "sf_published=N6+N5-residualised mean z, n=%d; %d in-band "
                 "sections, %s" % (int(s.n.sum()), s.section.nunique(),
                                   ", ".join(scope["curve_celltypes"]))))
    rows.append(dict(null="negative_control_probe_kernel_H1", family="PENDING",
                     note="Section 13 test A7, HUMAN half — needs H1, behind "
                          "the freeze"))

    D = pd.DataFrame(rows, columns=[
        "null", "family", "frac_retaining_a_neighbour", "real_median_nbrs",
        "null_median_nbrs", "frac_in_occupancy", "median_displacement_um",
        "sf_median", "sf_q25", "sf_q75", "sf_published", "note"])
    D.to_csv(f"{RES}/figure2e_data.csv", index=False)
    return D, lam, scope


# ---------------------------------------------------------------- draw ------

def draw(D, lam, scope):
    idx = D.set_index("null")
    labs = [r[1] for r in ROWS]
    fams = [r[2] for r in ROWS]
    y = np.arange(len(ROWS))[::-1]
    seps = [i for i in range(1, len(ROWS))
            if ROWS[i][0] in ("N4_orig", "N3_var", "N4_var")]

    fig = plt.figure(figsize=(19.0, 13.6))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.20, left=0.075,
                          right=0.975, top=0.845, bottom=0.135)

    def rowstyle(ax):
        ax.set_yticks(y)
        ax.set_yticklabels(labs, fontsize=10.5)
        for t, f in zip(ax.get_yticklabels(), fams):
            if f == "published bounding box":
                t.set_color(PAL.SERIES[7]); t.set_fontweight("bold")
            elif f == "PRIMARY (variance-corrected)":
                t.set_color(PAL.SERIES[5]); t.set_fontweight("bold")
            else:
                t.set_color(PAL.INK2)
        for i in seps:
            ax.axhline(y[i] + 0.5, color=PAL.AXIS, lw=1.0)
        ax.set_ylim(y.min() - 0.7, y.max() + 0.7)
        ax.grid(axis="y", visible=False)

    cols = [FAM_COL[f] for f in fams]

    # ---- e: does the null leave the tissue? ------------------------------
    ax = fig.add_subplot(gs[0, 0])
    v = [float(idx.loc[r[0], "frac_retaining_a_neighbour"]) for r in ROWS]
    ax.barh(y, v, color=cols, height=0.62)
    for yy, vv in zip(y, v):
        ax.text(vv + 0.012, yy, f"{vv:.3f}", va="center", fontsize=9,
                color=PAL.INK2)
    ax.axvline(1.0, color=PAL.MUTED, lw=1.0, ls=":")
    ax.set_xlim(0, 1.16)
    ax.set_xlabel("fraction of shifted senders with a real cell within the "
                  f"{RN.WINDOW_UM:.0f} µm window")
    ax.set_title("e   Does the null leave the tissue?", fontsize=13,
                 loc="left", fontweight="bold", pad=14)
    rowstyle(ax)

    # ---- f: does the null move anything? ---------------------------------
    ax = fig.add_subplot(gs[0, 1])
    v = [float(idx.loc[r[0], "median_displacement_um"]) for r in ROWS]
    ax.barh(y, v, color=cols, height=0.62)
    for yy, vv in zip(y, v):
        ax.text(vv * 1.10, yy, f"{vv:,.0f}", va="center", fontsize=9,
                color=PAL.INK2)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlim(0, 3.2e4)
    ax.axvline(lam["median"], color=PAL.INK, lw=1.8)
    ax.axvline(RN.WINDOW_UM, color=PAL.STATUS["critical"], lw=1.6, ls="-.")
    ax.set_xlabel("median distance a sender is actually moved (µm, log)")
    ax.set_title("f   Does the null move anything?", fontsize=13, loc="left",
                 fontweight="bold", pad=34)
    ax.text(lam["median"], 1.005,
            "median λ̂ %.1f µm\nIQR [%.1f, %.1f], %.0f %% railed"
            % (lam["median"], lam["q25"], lam["q75"], 100 * lam["railed"]),
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=9, color=PAL.INK, linespacing=1.35)
    ax.text(RN.WINDOW_UM * 1.15, 1.005, f"{RN.WINDOW_UM:.0f} µm\nwindow",
            transform=ax.get_xaxis_transform(), ha="left", va="bottom",
            fontsize=9, color=PAL.STATUS["critical"], linespacing=1.35)
    ax.text(0.0, -0.135,
            "N3-occ and N4-occ move senders barely further than λ̂.\n"
            "They are in tissue because they are near-identity, not because "
            "they are corrected.",
            transform=ax.transAxes, fontsize=9.5, color=PAL.SERIES[3],
            va="top", linespacing=1.5)
    rowstyle(ax)

    # ---- g: old vs new surviving fraction --------------------------------
    ax = fig.add_subplot(gs[1, 0])
    med = np.array([float(idx.loc[r[0], "sf_median"]) for r in ROWS])
    lo = np.array([float(idx.loc[r[0], "sf_q25"]) for r in ROWS])
    hi = np.array([float(idx.loc[r[0], "sf_q75"]) for r in ROWS])
    ax.barh(y, med, color=cols, height=0.62)
    ax.errorbar(0.5 * (lo + hi), y, xerr=0.5 * (hi - lo), fmt="none",
                ecolor=PAL.INK, elinewidth=1.6, capsize=0)
    for yy, vv, hh in zip(y, med, hi):
        ax.text(max(vv, hh) + 0.03, yy, f"{vv:.3f}", va="center", fontsize=9,
                color=PAL.INK2)
    pub = [(yy, float(idx.loc[r[0], "sf_published"]))
           for yy, r in zip(y, ROWS)
           if np.isfinite(float(idx.loc[r[0], "sf_published"]))]
    ax.scatter([p[1] for p in pub], [p[0] for p in pub], marker="D", s=48,
               color=PAL.INK, zorder=6)
    n1 = float(idx.loc["reference_N1_sf", "sf_median"])
    ax.axvline(n1, color=PAL.SERIES[1], lw=1.5, ls="-.")
    ax.text(n1, 1.005, f"N1 label\npermutation {n1:.3f}",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=9, color=PAL.SERIES[1], linespacing=1.35)
    ax.axvline(1.0, color=PAL.MUTED, lw=1.0, ls=":")
    ax.set_xlim(0, 2.0)
    ax.set_xlabel("surviving fraction of β̂   (bar = median, line = IQR,\n"
                  "◆ = the published value)", linespacing=1.6)
    ax.set_title("g   Old vs new surviving fraction", fontsize=13, loc="left",
                 fontweight="bold", pad=34)
    ax.text(0.0, -0.20,
            "N3-swap ≡ N1 — it destroys sender clustering, so it is NOT a "
            "corrected torus shift.",
            transform=ax.transAxes, fontsize=9.5, color=PAL.SERIES[1],
            va="top")
    rowstyle(ax)

    # ---- h: A7 negative-control-probe kernel ------------------------------
    ax = fig.add_subplot(gs[1, 1])
    C = D[D.null.str.startswith("A7_curve_bin_")].copy()
    C = C.sort_values("median_displacement_um")
    ax.plot(C.median_displacement_um, C.sf_median, "o-", ms=4.5, lw=1.8,
            color=PAL.SERIES[7], label="naive (intercept only) — NOT flat")
    ax.plot(C.median_displacement_um, C.sf_published, "o-", ms=4.5, lw=1.8,
            color=PAL.SERIES[2],
            label="after the N6+N5 nuisance design — flat")
    ax.axhline(0, color=PAL.MUTED, lw=1.2, ls="--")
    ax.set_xlabel("distance to nearest sender (µm)")
    ax.set_ylabel("negative-control counts (z, mean per 5 µm bin)")
    ax.set_title("h   Negative-control-probe kernel — M1 DONE, H1 pending",
                 fontsize=13, loc="left", fontweight="bold",
                 color=PAL.STATUS["good"], pad=14)
    ax.legend(loc="upper left", fontsize=9.5)

    A = pd.read_csv(f"{RES}/a7_summary.csv").set_index(["response", "design"])

    def a7(resp, d):
        r = idx.loc[f"A7_{resp}_{d}"]
        return r.sf_median, r.sf_q25, r.sf_q75, A.loc[
            ({"all_controls": "all_controls"}.get(resp, resp), d),
            "clustered_p"]

    nsec_curve = scope["curve_sections"]
    m, l, h, p = a7("all_controls", "base")
    m2, l2, h2, p2 = a7("all_controls", "n6n5")
    m3, _, _, p3 = a7("all_controls", "n2")
    bio_b = idx.loc["A7_BIOLOGICAL_base", "sf_median"]
    bio_c = idx.loc["A7_BIOLOGICAL_n6n5", "sf_median"]
    fam = " / ".join("%s %+.4f" % (t, idx.loc[f"A7_{k}_base", "sf_median"])
                     for t, k in (("probe", "neg_control_probe"),
                                  ("codeword", "neg_control_codeword"),
                                  ("genomic", "genomic_control")))
    box = (
        f"A7, Xenium 5K panel controls ({PANEL_CONTROLS}).\n"
        f"CURVE above: {nsec_curve} in-band sections, "
        f"{', '.join(scope['curve_celltypes'])}, sender call "
        f"{RN.PRIMARY_CALL}.\n"
        f"AMPLITUDES below: {scope['n_sections']} sections × "
        f"{scope['n_celltypes']} cell types — a WIDER scope than the curve.\n"
        f"β/sd, section-clustered mean [95% CI]:\n"
        f"   naive      {m:+.4f} [{l:+.4f}, {h:+.4f}]  p={p:g}   NOT flat\n"
        f"   +N6+N5     {m2:+.4f} [{l2:+.4f}, {h2:+.4f}]  p={p2:g}    flat\n"
        f"   N2 decoy   {m3:+.4f}  p={p3:g}   — N2 does NOT remove it\n"
        f"   per family, naive: {fam}\n"
        f"Tier B modules, same fits: {bio_b:+.4f} → {bio_c:+.4f}.  "
        f"H1 half behind the freeze.")
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_ - 0.78 * (hi_ - lo_), hi_)
    ax.text(0.985, 0.02, box, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.6, color=PAL.INK, linespacing=1.5, family="DejaVu Sans",
            bbox=dict(boxstyle="round,pad=0.5", fc=PAL.SURFACE,
                      ec=PAL.STATUS["good"], lw=1.4))

    handles = [Line2D([], [], marker="s", ls="", ms=10, color=c, label=f)
               for f, c in FAM_COL.items()]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=10.5,
               bbox_to_anchor=(0.5, 0.012))

    fig.suptitle("Figure 2e–h — the coordinate nulls audited: the published "
                 "torus shift left the tissue, and the\nin-tissue "
                 "replacements do not change the result", fontsize=16,
                 y=0.985, linespacing=1.5)
    fig.text(0.5, 0.922,
             "six Test-3-admissible sections, sender call "
             f"{RN.PRIMARY_CALL}; e–f are medians over sections, 20 draws "
             f"each; g is over the {int(pd.read_csv(f'{RES}/sf_summary.csv').iloc[0].n)} "
             "reportable fits, 1,000 permutations.  "
             "N3-var / N4-var are the FROZEN PRIMARY null (green).",
             fontsize=11, ha="center", color=PAL.INK2)
    fig.text(0.5, 0.900, "red = the published bounding-box nulls;   "
             "green = the frozen primary variance-corrected null;   "
             "λ̂ is quoted with its IQR and railing rate, never alone",
             fontsize=10, ha="center", color=PAL.INK2)

    os.makedirs(FIG, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/figure2e.{e}")
    print("wrote", f"{FIG}/figure2e.png", f"{FIG}/figure2e.pdf",
          f"{RES}/figure2e_data.csv")


def main():
    D, lam, scope = build_data()
    draw(D, lam, scope)


if __name__ == "__main__":
    main()
