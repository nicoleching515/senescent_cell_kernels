#!/usr/bin/env python3
"""Figure 5 (new, Phase7 §19) — TWO-ARM REPLICATION.

(a) Surviving fraction by null, mouse and human side by side.
(b) The §17 table as a forest plot: the controlled amplitude of each Tier B module on each
    arm, against that arm's OWN 80 %-power detectable bound (the R2 / R4 evaluation).
(c) The geometric predictions (R3): the Poisson identity, and the lambda-hat railing rate,
    against their expected values in both arms.

Palette = the project's pre-validated categorical theme (`sasp_palette`, the dataviz-skill
reference instance already validated in Phase 1); hues assigned in fixed order, never cycled.
Every plotted number is written to figures/figure5_data.csv; both .png and .pdf are emitted
(PREREG §11).  Nothing here computes a statistic: every value is read from a results file.

M1 primary call: tierA_p95 (fine labels).  H1 primary call: tierAmg_p95 == PREREG D-B's
tierA_merged_p95; the frozen-literal tierA_p95 is drawn as an open marker beside it.

Usage: python3 code/make_figure5.py
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
RH = "/workspace/results/phase10_h1"
M1_INBAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
             "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
             "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
H1_SECS = ["SPLN07", "SPLN14", "SPLN21", "SPLN24", "SPLN30", "SPLN43", "SPLN44"]
C_M1, C_H1, C_H1B = PAL.SERIES[0], PAL.SERIES[1], PAL.SERIES[3]
ROWS = []


def rep(res, call, sections):
    d = pd.read_csv(os.path.join(res, "main_fits.csv"))
    d = d[d.section.isin(sections) & (d.call == call) & (d.stratum == "all")]
    return d, d[(d.beta_naive > 0) & (d.sf_base.notna()) & (d.beta_base_lo > 0)]


def perm_sf(res, fn, call, sections, key, col, scope=None):
    p = os.path.join(res, fn)
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p)
    d = d[d.section.isin(sections) & (d.call == call)]
    if scope and "scope" in d.columns:
        d = d[d.scope == scope]
    if col not in d.columns:
        return None
    d = d.merge(key, on=["section", "celltype", "module"])
    v = d[col].dropna().to_numpy(float)
    return v if v.size else None


NULLS = [
    ("N1 label permutation",            "perm", "N1_sf",      None),
    ("N3 torus, published bbox",        "perm", "N3_sf",      None),
    ("N4 rotation, published bbox",     "perm", "N4_sf",      None),
    ("N3-var (PRIMARY corrected)",      "var",  "N3_var_sf",  "full"),
    ("N4-var (PRIMARY corrected)",      "var",  "N4_var_sf",  "full"),
    ("N3-tile (in-tissue, tile scope)", "c1",   "N3_tile_sf", "tile"),
    ("N4-tile (in-tissue, tile scope)", "c1",   "N4_tile_sf", "tile"),
    ("N2 matched decoy",                "main", "sf_n2",      None),
    ("N6 receiver baseline",            "main", "sf_n6",      None),
    ("anatomical covariate alone",      "main", "sf_zon",     None),
    ("N5 nuisance covariates",          "main", "sf_n5",      None),
    ("N5 + N6",                         "main", "sf_n6n5",    None),
    ("N2 + N5 + N6  (PRIMARY)",         "main", "sf_n2n5n6",  None),
]
FILES = {"perm": "perm_nulls.csv", "c1": "perm_nulls_c1.csv", "var": "perm_nulls_var.csv"}


def collect(res, call, sections):
    d, r = rep(res, call, sections)
    key = r[["section", "celltype", "module"]].drop_duplicates()
    out = {}
    for lab, fam, col, scope in NULLS:
        if fam == "main":
            v = r[col].dropna().to_numpy(float) if col in r else None
        else:
            v = perm_sf(res, FILES[fam], call, sections, key, col, scope)
        if v is not None and v.size:
            out[lab] = v
    return d, r, out


def main():
    d_m1, r_m1, sf_m1 = collect(R3, "tierA_p95", M1_INBAND)
    have_h1 = os.path.exists(os.path.join(RH, "main_fits.csv"))
    if have_h1:
        d_h1, r_h1, sf_h1 = collect(RH, "tierAmg_p95", H1_SECS)
        d_h1b, r_h1b, sf_h1b = collect(RH, "tierA_p95", H1_SECS)
    else:
        raise SystemExit("results/phase10_h1/main_fits.csv missing; run stage main first")

    fig = plt.figure(figsize=(13.4, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.28, 1.0], width_ratios=[1.05, 1.0],
                          hspace=0.34, wspace=0.26)

    # ---------------- (a) SF by null, both arms ---------------------------
    ax = fig.add_subplot(gs[0, :])
    labs = [l for l, *_ in NULLS if l in sf_m1 or l in sf_h1]
    y = np.arange(len(labs))[::-1]
    ax.axvline(1.0, color=PAL.MUTED, lw=1, ls="--", zorder=1)
    ax.axvline(0.0, color=PAL.INK2, lw=1.1, zorder=1)
    for series, colr, off, name in ((sf_m1, C_M1, +0.20, "M1 mouse liver (tierA_p95)"),
                                    (sf_h1, C_H1, 0.00, "H1 human spleen (tierAmg_p95)"),
                                    (sf_h1b, C_H1B, -0.20,
                                     "H1, frozen-literal tierA_p95")):
        for yi, lab in zip(y, labs):
            v = series.get(lab)
            if v is None:
                continue
            q1, md, q3 = np.quantile(v, [.25, .5, .75])
            ax.plot([q1, q3], [yi + off] * 2, "-", color=colr, lw=2.6, alpha=.55,
                    solid_capstyle="butt", zorder=2)
            ax.plot([md], [yi + off], "o", color=colr, ms=6.5, mec=PAL.SURFACE, mew=1.3,
                    zorder=3)
            ROWS.append(dict(panel="a", arm=name, null=lab, n=int(v.size),
                             median=round(float(md), 5), q25=round(float(q1), 5),
                             q75=round(float(q3), 5)))
        ax.plot([], [], "o-", color=colr, label=name, ms=6.5)
    ax.set_yticks(y); ax.set_yticklabels(labs)
    ax.set_xlabel("surviving fraction of the naive distance amplitude "
                  "(1 = null removes nothing, 0 = null removes it all)")
    ax.set_xlim(-0.45, 1.35)
    ax.set_title("(a)  Surviving fraction by null, both arms.  Point = median over that "
                 "arm's reportable fits; bar = inter-quartile range across fits (NOT a "
                 "confidence interval)", loc="left")
    ax.legend(loc="lower left", ncol=3)

    # ---------------- (b) forest: controlled amplitude per module ---------
    ax = fig.add_subplot(gs[1, 0])
    mods = sorted(set(r_m1.module) | set(r_h1.module))
    y = np.arange(len(mods))[::-1]
    bnd = {}
    for r, colr, off, name in ((r_m1, C_M1, +0.18, "M1 (tierA_p95)"),
                               (r_h1, C_H1, -0.18, "H1 (tierAmg_p95)")):
        se = ((r.beta_n2n5n6_hi - r.beta_n2n5n6_lo) / (2 * 1.959964) / r.sd_y)
        bnd[name] = float(2.802 * se.median())
        for yi, m in zip(y, mods):
            g = r[r.module == m]
            if not len(g):
                continue
            a = (g.beta_n2n5n6 / g.sd_y)
            q1, md, q3 = a.quantile(.25), a.median(), a.quantile(.75)
            ax.plot([q1, q3], [yi + off] * 2, "-", color=colr, lw=2.6, alpha=.55, zorder=2)
            ax.plot([md], [yi + off], "o", color=colr, ms=6.5, mec=PAL.SURFACE, mew=1.3,
                    zorder=3)
            ROWS.append(dict(panel="b", arm=name, module=m, n=int(len(g)),
                             ctrl_amp_median=round(float(md), 5),
                             q25=round(float(q1), 5), q75=round(float(q3), 5),
                             arm_power80_bound=round(bnd[name], 5)))
        ax.axvline(bnd[name], color=colr, lw=1.4, ls=":", zorder=1)
        ax.plot([], [], "o-", color=colr, label=name, ms=6.5)
    ax.axvline(0, color=PAL.INK2, lw=1.1, zorder=1)
    ax.set_yticks(y); ax.set_yticklabels([m.replace("_", " ") for m in mods], fontsize=8)
    ax.set_xlabel(r"controlled amplitude $|\beta|/\mathrm{sd}(y)$ under N2+N5+N6")
    ax.set_title("(b)  §17 as a forest plot: each Tier B module against its own arm's\n"
                 "80 %-power detectable bound (dotted).  Nothing crosses it on either arm.",
                 loc="left")
    ax.legend(loc="lower right")

    # ---------------- (c) geometric predictions ---------------------------
    gsc = gs[1, 1].subgridspec(1, 2, wspace=0.42)
    ax = fig.add_subplot(gsc[0, 0])
    for res, colr, name in ((R3, C_M1, "M1"), (RH, C_H1, "H1")):
        p = os.path.join(res, "poisson_density.csv")
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p)
        x = np.log(d.sender_density_per_um2.to_numpy(float))
        yv = np.log(d.median_d_um.to_numpy(float))
        ax.plot(x, yv, "o", color=colr, ms=4.2, alpha=.75, mec=PAL.SURFACE, mew=.7,
                label=name)
        f = pd.read_csv(os.path.join(res, "poisson_fits.csv"))
        f = f[f.subset.str.contains("ALL sections", na=False)].iloc[0]
        xs = np.linspace(x.min(), x.max(), 2)
        ax.plot(xs, f.intercept + f.slope * xs, "-", color=colr, lw=1.6, alpha=.9)
        ROWS.append(dict(panel="c_poisson", arm=name, n=int(len(d)),
                         slope=round(float(f.slope), 4), r2=round(float(f.r2), 4),
                         intercept=round(float(f.intercept), 4)))
    xs = np.array(ax.get_xlim())
    ax.plot(xs, np.log(0.4697) - 0.5 * xs, "--", color=PAL.INK2, lw=1.2,
            label="homogeneous Poisson\n(slope $-1/2$)")
    ax.set_xlabel("log sender density (µm$^{-2}$)")
    ax.set_ylabel("log median distance to\nnearest sender (µm)")
    ax.set_title("(c1)  R3a Poisson identity", loc="left")
    ax.legend(loc="upper right", fontsize=7)

    ax = fig.add_subplot(gsc[0, 1])
    bars = []
    for d, colr, name in ((d_m1, C_M1, "M1\ntierA_p95"),
                          (d_h1, C_H1, "H1\ntierAmg_p95"),
                          (d_h1b, C_H1B, "H1\ntierA_p95")):
        lo = float((d.lam_naive <= d.lam_grid_lo + 1e-9).mean())
        hi = float((d.lam_naive >= d.lam_grid_hi - 1e-9).mean())
        bars.append((name, lo, hi, colr, float(d.lam_naive.median()), len(d)))
    xp = np.arange(len(bars))
    ax.bar(xp, [b[1] for b in bars], .62, color=[b[3] for b in bars], alpha=.95,
           label="railed at the 7 µm floor")
    ax.bar(xp, [b[2] for b in bars], .62, bottom=[b[1] for b in bars],
           color=[b[3] for b in bars], alpha=.42, hatch="///",
           edgecolor=PAL.SURFACE, label="railed at the 50 µm ceiling")
    for i, b in enumerate(bars):
        ax.text(i, b[1] + b[2] + .02, "%.0f %%\nmed λ̂ %.1f µm" % (100 * (b[1] + b[2]), b[4]),
                ha="center", fontsize=7.5, color=PAL.INK2)
        ROWS.append(dict(panel="c_railing", arm=b[0].replace("\n", " "), n_fits=b[5],
                         frac_railed_floor=round(b[1], 4), frac_railed_ceiling=round(b[2], 4),
                         median_lam_naive_um=round(b[4], 3)))
    ax.set_xticks(xp); ax.set_xticklabels([b[0] for b in bars], fontsize=7.5)
    ax.set_ylim(0, 1.05); ax.set_ylabel("fraction of fits with λ̂ at a grid bound")
    ax.set_title("(c2)  R3b λ̂ railing", loc="left")
    ax.legend(loc="lower left", fontsize=7)

    fig.suptitle("Figure 5 — Two-arm replication: M1 mouse liver (GSE310392, 6 admissible "
                 "sections) against H1 human spleen (GSE326743, 7 sections)\n"
                 "Species and tissue are confounded by design and no difference below may "
                 "be attributed to either (PREREG §10.5).", fontsize=10, y=0.995)
    for e in ("png", "pdf"):
        fig.savefig(FIG + "figure5." + e, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(ROWS).to_csv(FIG + "figure5_data.csv", index=False)
    print("wrote figures/figure5.png/.pdf and figure5_data.csv (%d rows)" % len(ROWS))


if __name__ == "__main__":
    main()
