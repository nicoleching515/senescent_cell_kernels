#!/usr/bin/env python3
"""
Methods figures for the frozen gene sets (Phase 7 section 11 and section 17).

  gs1  disjointness / intersection matrices, both arms
  gs2  cross-arm gene-set symmetry, pre- and post-C6, with ortholog-intersected counts
  gs3  CoreScence circularity, both arms, pre-C6 -> frozen C6
  gs4  SenePy spleen coverage

Every value is read from a CSV/JSON produced by the build/gate scripts -- nothing is recomputed
here and nothing is typed in. (This was violated once: gs3's mouse bar was the string literal
`24/35`. It is now produced by `code/corescence_circularity.py`. See AUDIT_PHASE8_FACTCHECK M1.)
Each figure writes figures/<name>.png, .pdf and <name>_data.csv, and the _data.csv carries every
number the figure draws, including the caption block (see AUDIT_PHASE8_FACTCHECK M5).

Run: python3 /workspace/code/make_figure_genesets.py
"""
from __future__ import annotations
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, "/workspace/code")
import sasp_palette as PAL

PAL.apply_style(matplotlib)
RES = "/workspace/results/phase7_jobA"
FIG = "/workspace/figures"
BLAB = ["B1", "B2", "B3", "B4", "B5", "B6", "B7"]
BFULL = ["B1 tnfa_nfkb", "B2 il6_jak_stat3", "B3 interferon", "B4 downstream_arrest",
         "B5 emt_ecm", "B6 oxidative", "B7 secondary_sen"]


def save(fig, name):
    for e in ("png", "pdf"):
        fig.savefig(f"{FIG}/{name}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG}/{name}.png/.pdf")


def heat(ax, M, rowlab, collab, title, vmax=None, zero_is_meaningful=True):
    """Sequential SEQ heat map. Exact zeros are painted as a flat neutral so that
    'zero' is unambiguously distinguishable from 'small but non-zero' -- the whole
    point of the A x B panel."""
    M = np.asarray(M, float)
    vmax = vmax if vmax is not None else max(1.0, np.nanmax(M))
    masked = np.ma.masked_where(M == 0, M) if zero_is_meaningful else M
    cmap = PAL.SEQ.copy()
    cmap.set_bad(PAL.SURFACE)
    im = ax.imshow(masked, cmap=cmap, norm=Normalize(vmin=0, vmax=vmax), aspect="auto")
    ax.set_xticks(range(len(collab)))
    ax.set_xticklabels(collab, rotation=90, fontsize=7)
    ax.set_yticks(range(len(rowlab)))
    ax.set_yticklabels(rowlab, fontsize=7)
    ax.set_title(title, fontsize=8.5, color=PAL.INK)
    ax.grid(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v == 0 and zero_is_meaningful:
                ax.text(j, i, "0", ha="center", va="center", fontsize=6.5,
                        color=PAL.MUTED, weight="bold")
            else:
                ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6.5,
                        color="white" if v > 0.55 * vmax else PAL.INK)
    for s in ax.spines.values():
        s.set_edgecolor(PAL.AXIS)
    return im


# ---------------------------------------------------------------- gs1
def gs1():
    out = []
    fig, axes = plt.subplots(2, 3, figsize=(13.6, 8.4),
                             gridspec_kw=dict(width_ratios=[1.05, 1.05, 0.95]))
    for r, (arm, ab_csv, bb_csv, bbpre_csv) in enumerate([
            ("M1 mouse liver", f"{RES}/intersection_matrix_mouse.csv",
             f"{RES}/intersection_matrix_mouse_BxB.csv",
             f"{RES}/intersection_matrix_mouse_BxB_preC6.csv"),
            ("H1 human spleen", f"{RES}/intersection_matrix_human.csv",
             f"{RES}/intersection_matrix_BxB.csv",
             f"{RES}/intersection_matrix_BxB_preC6.csv")]):
        ab = pd.read_csv(ab_csv)
        cols = [c for c in ab.columns if c[:2] in BLAB and not c.endswith("_genes")][:7]
        NAMECOL = "tier_A_set"
        pre = ab[ab[NAMECOL].str.startswith(("A0", "A_ported")) |
                 ab[NAMECOL].str.contains("FAILED")].copy()
        # stable order across arms: candidate pool first, then the section 10 list that failed
        pre["_o"] = np.where(pre[NAMECOL].str.contains("FAILED"), 1, 0)
        pre = pre.sort_values("_o")
        frozen = ab[ab[NAMECOL].str.contains("PRIMARY|sensitivity", regex=True)]

        Mpre = pre[cols].to_numpy()
        heat(axes[r, 0], Mpre,
             [f"{s.split(' (')[0]}  (n={n})" for s, n in zip(pre[NAMECOL], pre["n_A_on_panel"])],
             BLAB, f"{arm}\n(a) candidate Tier A definitions x Tier B, BEFORE removal",
             vmax=Mpre.max())

        Mf = frozen[cols].to_numpy()
        rl = [f"{s.split(' (')[0]}  (n={n})" for s, n in zip(frozen[NAMECOL], frozen["n_A_on_panel"])]
        axb = axes[r, 1]
        heat(axb, Mf, rl, BLAB,
             "(b) FROZEN sets x Tier B — GATE CELLS OUTLINED, all 0\n"
             "(per-module sets need disjointness only from their own module)",
             vmax=max(1.0, Mf.max()))
        # mark the cells the section 11 gate actually asserts: every cell of the PRIMARY row,
        # and the own-module cell of each per-module sensitivity row.
        for i, name in enumerate(frozen[NAMECOL]):
            if "PRIMARY" in name:
                js = range(len(cols))
            else:
                mod = name.split("A_sender_for_")[1].split(" (")[0]
                js = [k for k, c in enumerate(cols) if c[3:] == mod or c.endswith(mod)]
            for j in js:
                axb.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                            edgecolor=PAL.STATUS["good"], lw=1.6, zorder=5))
                assert Mf[i, j] == 0, (arm, name, cols[j], Mf[i, j])

        bb = pd.read_csv(bb_csv).set_index("module")
        bbp = pd.read_csv(bbpre_csv).set_index("module")
        Mbb = bb[BLAB].to_numpy().astype(float)
        Mbbp = bbp[BLAB].to_numpy().astype(float)
        off = Mbb.copy()
        np.fill_diagonal(off, 0)
        heat(axes[r, 2], off,
             [f"{m.split('_', 1)[0]} {m.split('_', 1)[1][:14]} (n={n})"
              for m, n in zip(bb.index, bb["n_on_panel"])],
             BLAB, "(c) Tier B x Tier B cross-talk (diagonal blanked)\nnot required to be 0",
             vmax=max(1, off.max()))
        for i, m in enumerate(bb.index):
            for j, l in enumerate(BLAB):
                if i == j:
                    continue
                out.append(dict(arm=arm, panel="BxB", row=m, col=l,
                                n_frozen=int(Mbb[i, j]), n_preC6=int(Mbbp[i, j])))
        for df, panel in ((pre, "A0xB"), (frozen, "frozenAxB")):
            for _, row in df.iterrows():
                for c in cols:
                    out.append(dict(arm=arm, panel=panel, row=row[NAMECOL], col=c[:2],
                                    n_frozen=int(row[c]), n_preC6=""))
    fig.suptitle("Tier A / Tier B disjointness — Phase 7 §11 gate, both arms (frozen 2026-08-27)",
                 fontsize=11, color=PAL.INK, y=0.995)
    fig.text(0.5, -0.075,
             "(b) is the §11 gate. Green-outlined cells are what the gate asserts: the PRIMARY "
             "sender set must be disjoint from all seven modules, each per-module sensitivity set "
             "from its own module only. Every outlined cell is 0 on both arms (asserted in code, "
             "not eyeballed). Unoutlined non-zeros are expected and are not gate failures. "
             "(a) shows the collisions that had to be removed, and the §10 Tier A that failed. "
             "(c) is where C6's cost sits: B7∩B1 rose 10→26 and B7∩B5 6→20 in human — the "
             "n_preC6 column of the data file carries the pre-C6 values.",
             ha="center", fontsize=7.5, color=PAL.INK2, wrap=True)
    fig.tight_layout()
    save(fig, "figure_gs1_intersection_matrix")
    pd.DataFrame(out).to_csv(f"{FIG}/figure_gs1_intersection_matrix_data.csv", index=False)


# ---------------------------------------------------------------- gs2
def gs2():
    t = pd.read_csv(f"{RES}/crossarm_geneset_table.csv")
    x = json.load(open(f"{RES}/crossarm_geneset_table.json"))
    keep = t[t["set"].str.startswith(("A0", "A_SENDER_FINAL_strict", "B_"))].copy()
    lab = [s.replace("A_SENDER_FINAL_strict (PRIMARY)", "Tier A strict (PRIMARY)")
            .replace("A0 candidate pool (pre-disjointness)", "Tier A pool A0")
            .replace("B_", "B ") for s in keep["set"]]
    n = len(keep)
    y = np.arange(n)[::-1]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.9),
                             gridspec_kw=dict(width_ratios=[1.5, 1.5, 1.0]))

    def num(c):
        return pd.to_numeric(keep[c], errors="coerce").to_numpy(float)

    for ax, arm, pre, post in ((axes[0], "M1 mouse liver", num("mouse_preC6"), num("mouse_C6")),
                               (axes[1], "H1 human spleen", num("human_preC6"), num("human_C6"))):
        ax.barh(y + 0.20, pre, height=0.38, color=PAL.MUTED, label="pre-C6")
        ax.barh(y - 0.20, post, height=0.38, color=PAL.SERIES[0], label="frozen (C6)")
        # SHARED-MEMBERSHIP OVERLAY. Equal bar lengths read as identity; they are not. For the
        # two headline sets, draw how many of THIS arm's genes the other arm actually also has
        # (after ortholog mapping) as a dark inner segment, against this arm's own total.
        # Shared count is the MAP-GAP-CORRECTED one: a gene both arms carry is shared even
        # when the pinned MGI report has no row for it (AUDIT R2). The uncorrected
        # pinned-map count is carried alongside in the data CSV.
        shared = {d["set"]: d["n_overlap_map_gap_corrected"] for d in x["asymmetry"]
                  if "n_overlap_map_gap_corrected" in d}
        ov = {"A_SENDER_FINAL_strict (PRIMARY)": shared.get("A_SENDER_FINAL_strict_FROZEN"),
              "B_secondary_senescence": shared.get("B7_FROZEN")}
        for yy, a, b, nm in zip(y, pre, post, keep["set"]):
            if np.isfinite(a):
                ax.text(a + 3, yy + 0.20, f"{int(a)}", va="center", fontsize=6.5, color=PAL.INK2)
            if np.isfinite(b):
                ax.text(b + 3, yy - 0.20, f"{int(b)}", va="center", fontsize=6.5,
                        color=PAL.SERIES[0], weight="bold")
            if nm in ov and ov[nm] is not None and np.isfinite(b):
                k = ov[nm]
                ax.barh(yy - 0.20, k, height=0.38, color=PAL.INK, alpha=0.9, zorder=3)
                ax.text(1.5, yy - 0.56, f"only {k} of {int(b)} shared with the other arm",
                        va="center", ha="left", fontsize=6.4,
                        color=PAL.STATUS["critical"], weight="bold", zorder=4)
        ax.axvline(30, color=PAL.STATUS["critical"], lw=1.0, ls="--")
        ax.text(30, n - 0.2, " §11 floor = 30", fontsize=6.5, color=PAL.STATUS["critical"])
        ax.set_yticks(y)
        ax.set_yticklabels(lab if ax is axes[0] else [""] * n, fontsize=7.5)
        ax.set_xlabel("genes on that arm's panel")
        ax.set_title(arm, fontsize=9)
        ax.legend(loc="center right", fontsize=7)
        ax.grid(axis="y", visible=False)

    om = pd.to_numeric(keep["orthologue_intersected_mouse_C6"], errors="coerce").to_numpy(float)
    oh = pd.to_numeric(keep["orthologue_intersected_human_C6"], errors="coerce").to_numpy(float)
    ax = axes[2]
    ax.barh(y + 0.20, om, height=0.38, color=PAL.SERIES[2], label="mouse")
    ax.barh(y - 0.20, oh, height=0.38, color=PAL.SERIES[1], label="human")
    ax.set_yticks(y)
    ax.set_yticklabels([""] * n)
    ax.set_xlabel(f"genes on the ortholog-intersected panel (n={x['ortholog_intersected_panel']})")
    ax.set_title("(c) test A8 denominator", fontsize=9)
    ax.legend(loc="lower right", fontsize=7)
    ax.grid(axis="y", visible=False)

    a = [d for d in x["asymmetry"] if d["set"] == "A_SENDER_FINAL_strict_FROZEN"][0]
    b7 = [d for d in x["asymmetry"] if d["set"] == "B7_FROZEN"][0]
    axes[0].set_title("M1 mouse liver\n(dark segment = genes H1 also has, after ortholog mapping)",
                      fontsize=9)
    axes[1].set_title("H1 human spleen\n(dark segment = genes M1 also has, after ortholog mapping)",
                      fontsize=9)
    fig.suptitle("Cross-arm gene-set symmetry, pre- and post-C6 — Phase 7 §17", fontsize=11,
                 color=PAL.INK, y=1.02)
    fig.text(0.5, -0.135,
             "EQUAL SIZE IS NOT IDENTITY. The two frozen Tier A sets are both n=33 but share only "
             f"{a['n_overlap_map_gap_corrected']} members "
             f"({a['n_overlap']} through the pinned 1:1 MGI map, plus {a['human_only_map_gap']}, "
             "which both arms carry and the map has no row for). "
             f"Mouse-only, complete: {a['mouse_only_complete']} "
             f"({len(a['mouse_only_complete'].split())} = 33 − "
             f"{a['n_overlap_map_gap_corrected']}; H2afx is real — human H2AFX is in B4 and is "
             "removed from human Tier A by disjointness).  "
             f"Human-only: {a['human_only_real']}.  "
             f"Frozen B7 is {b7['n_mouse']} vs {b7['n_human']} with "
             f"{b7['n_overlap_map_gap_corrected']} shared ({b7['n_overlap']} by the pinned map "
             f"+ {b7['human_only_map_gap']}) — the residual gap is driven by REACTOME_SASP having "
             "40 mouse members and 111 human ones. The dark segments plot the map-gap-corrected "
             "shared counts; panel (c) is a DIFFERENT quantity (set ∩ ortholog-intersected "
             "panel) that happens to coincide at 88 for human B7.",
             ha="center", fontsize=7.5, color=PAL.INK2, wrap=True)
    fig.tight_layout()
    save(fig, "figure_gs2_crossarm_symmetry")
    # AUDIT M5: the dark inner segments and the whole footnote block used to live only in
    # crossarm_geneset_table.json, so a reader auditing this CSV could not check them. They
    # are appended here as long-format rows, tier = PLOTTED / FOOTNOTE.
    keep = keep.copy()
    keep["dark_segment_shared_with_other_arm"] = [
        ov.get(nm, "") if nm in ov else "" for nm in keep["set"]]
    extra = []
    for d, lab in ((a, "Tier A strict (PRIMARY)"), (b7, "B_secondary_senescence")):
        for k in ("n_mouse", "n_human", "n_mouse_mapped", "n_overlap",
                  "n_overlap_map_gap_corrected", "human_only", "human_only_map_gap",
                  "human_only_real", "mouse_only", "mouse_unmapped", "mouse_only_map_gap",
                  "mouse_only_real", "mouse_only_complete"):
            if k in d:
                extra.append(dict(tier="PLOTTED" if k.startswith("n_overlap") else "FOOTNOTE",
                                  set=lab, metric=k, value=d[k]))
    extra.append(dict(tier="PLOTTED", set="ortholog-intersected panel", metric="n_genes",
                      value=x["ortholog_intersected_panel"]))
    pd.concat([keep, pd.DataFrame(extra)], ignore_index=True).to_csv(
        f"{FIG}/figure_gs2_crossarm_symmetry_data.csv", index=False)


# ---------------------------------------------------------------- gs3
def gs3():
    g = json.load(open(f"{RES}/gate_result_human.json"))["corescence"]
    npanel = g["n_on_panel"]
    # AUDIT M1: the mouse bar used to be the literal (24, 35). It is now derived from files by
    # code/corescence_circularity.py -> corescence_circularity_mouse.json. 35 is reproducible
    # under no mapping convention; the project's own committed n8_disjointness_*.csv says 33.
    mc = json.load(open(f"{RES}/corescence_circularity_mouse.json"))
    mconv = mc["convention_to_cite"]
    mpre = mc["configurations"]["pre_C6"][mconv]
    mc6 = mc["configurations"]["C6_promoted"][mconv]
    bars = [("M1 mouse\npre-C6 B7 (n=38)", mpre["n_in_any_B"], mpre["n_on_panel"], PAL.MUTED),
            ("M1 mouse\nFROZEN B7 (n=108)", mc6["n_in_any_B"], mc6["n_on_panel"],
             PAL.SERIES[2]),
            ("H1 human\nsuperseded B7 (n=35)", g["superseded_n_in_any_B"], npanel, PAL.SERIES[0]),
            ("H1 human\nFROZEN B7 (n=116)", g["frozen_n_in_any_B"], npanel,
             PAL.STATUS["critical"])]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.4),
                             gridspec_kw=dict(width_ratios=[1.15, 1.05]))
    ax = axes[0]
    for i, (lab, k, tot, col) in enumerate(bars):
        f = 100 * k / tot
        ax.bar(i, f, color=col, width=0.62)
        ax.text(i, f + 1.6, f"{f:.0f}%\n{k}/{tot}", ha="center", fontsize=8.5,
                color=PAL.INK, weight="bold")
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([b[0] for b in bars], fontsize=7.5)
    ax.set_ylabel("CoreScence genes in ≥1 Tier B module (%)")
    ax.set_ylim(0, 104)
    ax.set_title("(a) DeepScence's own gene set vs the response modules, both arms", fontsize=9)
    ax.grid(axis="x", visible=False)

    per = g["per_module_frozen"]
    order = sorted(per, key=lambda m: -len(per[m]))
    ax = axes[1]
    yy = np.arange(len(order))[::-1]
    cols = [PAL.STATUS["critical"] if m == "secondary_senescence" else PAL.SERIES[0]
            for m in order]
    ax.barh(yy, [len(per[m]) for m in order], color=cols, height=0.6)
    for y_, m in zip(yy, order):
        ax.text(len(per[m]) + 0.25, y_, f"{len(per[m])}/{npanel} = {len(per[m])/npanel:.2f}",
                va="center", fontsize=7.5, color=PAL.INK2)
    ax.set_yticks(yy)
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("CoreScence genes shared with that module")
    ax.set_title("(b) H1 human, frozen Tier B, per module — B7 alone contributes "
                 f"{len(per['secondary_senescence'])}/{npanel}\n"
                 "(the mouse per-module counts are in the _data.csv)", fontsize=9)
    ax.grid(axis="y", visible=False)

    fig.suptitle("CoreScence circularity — a disclosed cost of the adopted C6 decision",
                 fontsize=11, color=PAL.INK, y=1.02)
    fig.text(0.5, -0.125,
             "DeepScence is the sender caller for Job B, so overlap between its gene set and the "
             "readouts is circularity by construction. Re-sourcing B7 from SenMayo ∪ Reactome SASP "
             "fixed the sender/readout collision but raised this figure on BOTH arms: mouse "
             f"{100 * mpre['frac']:.0f}% → {100 * mc6['frac']:.0f}%, human 76% → 88%. CITE 88% "
             "— it is the configuration that will be run. The mouse bars are the "
             f"{mc['n_on_panel_fallback']} CoreScence genes reachable on the 5,097-gene mouse "
             f"panel ({mc['n_on_panel_strict']} through the pinned 1:1 MGI map plus "
             f"{' and '.join(mc['fallback_extra'])}, which the map has no row for); the same "
             "convention code/run_phase3_n8.py uses, and the denominator in the committed "
             "results/phase3/n8_disjointness_*.csv. An earlier version of this panel showed a "
             "typed-in 24/35 = 69% for the mouse arm; that denominator is reproducible under no "
             "convention and understated the pre-C6 mouse circularity by ~10 points. CoreScence "
             "is human, so on H1 it runs natively and none of the human numbers is "
             "ortholog-mapping loss. Strip-and-refit is in the frozen run order.",
             ha="center", fontsize=7.5, color=PAL.INK2, wrap=True)
    fig.tight_layout()
    save(fig, "figure_gs3_corescence_circularity")
    src = {0: "corescence_circularity_mouse.json (pre_C6/%s)" % mconv,
           1: "corescence_circularity_mouse.json (C6_promoted/%s)" % mconv,
           2: "gate_result_human.json (superseded)",
           3: "gate_result_human.json (frozen)"}
    pd.DataFrame([dict(configuration=b[0].replace("\n", " "), n_circular=b[1], n_on_panel=b[2],
                       pct=round(100 * b[1] / b[2], 1), source=src[i])
                  for i, b in enumerate(bars)] +
                 [dict(configuration=f"H1 human frozen per-module: {m}", n_circular=len(per[m]),
                       n_on_panel=npanel, pct=round(100 * len(per[m]) / npanel, 1),
                       source="gate_result_human.json (per_module_frozen)")
                  for m in order] +
                 [dict(configuration="M1 mouse frozen per-module: %s" % m,
                       n_circular=len(v), n_on_panel=mc6["n_on_panel"],
                       pct=round(100 * len(v) / mc6["n_on_panel"], 1),
                       source="corescence_circularity_mouse.json (C6_promoted/%s)" % mconv)
                  for m, v in mc6["per_module"].items()]
                 ).to_csv(f"{FIG}/figure_gs3_corescence_circularity_data.csv", index=False)


# ---------------------------------------------------------------- gs4
def gs4():
    d = pd.read_csv(f"{RES}/senepy_spleen_coverage.csv")
    j = json.load(open(f"{RES}/senepy_spleen_coverage.json"))
    d = d[d["in_label_set"].astype(str).str.lower() == "true"].copy()
    d["on_panel"] = pd.to_numeric(d["on_panel"], errors="coerce").fillna(0)
    d["hub"] = np.where(d["usable"] == "yes",
                        d["hub_tissue"].fillna("") + " / " + d["hub_cell"].fillna(""), "—")
    d = d.sort_values(["usable", "on_panel"], ascending=[True, False])
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.6),
                             gridspec_kw=dict(width_ratios=[1.55, 1.0]))
    yy = np.arange(len(d))[::-1]
    cols = [PAL.SERIES[0] if u == "yes" else PAL.STATUS["critical"] for u in d["usable"]]
    ax = axes[0]
    ax.barh(yy, d["on_panel"], color=cols, height=0.62)
    for y_, v, hub, u in zip(yy, d["on_panel"], d["hub"], d["usable"]):
        ax.text(max(v, 0) + 12, y_, f"{int(v)}   {hub}" if u == "yes" else "no hub in any tissue",
                va="center", fontsize=7, color=PAL.INK2 if u == "yes" else PAL.STATUS["critical"])
    ax.axvline(j["min_on_panel"], color=PAL.MUTED, ls="--", lw=1.0)
    ax.set_yticks(yy)
    ax.set_yticklabels(d["cell_type"], fontsize=7.5)
    ax.set_xlabel("hub genes on the human panel (mouse arm's threshold = 10, dashed)")
    ax.set_xlim(0, float(d["on_panel"].max()) * 1.55)
    ax.set_title("(a) the 22 spleen cell types vs the best SenePy hub available", fontsize=9)
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    grp = (d[d["usable"] == "yes"].groupby("hub")["cell_type"]
           .apply(lambda s: sorted(s)).sort_values(key=lambda s: s.map(len), ascending=True))
    yy2 = np.arange(len(grp))
    ax.barh(yy2, [len(v) for v in grp], color=PAL.SERIES[3], height=0.6)
    for y_, (hub, v) in zip(yy2, grp.items()):
        ax.text(len(v) + 0.06, y_, ", ".join(x[:32] for x in v), va="center", fontsize=6.5,
                color=PAL.INK2)
    ax.set_yticks(yy2)
    ax.set_yticklabels(grp.index, fontsize=7.5)
    ax.set_xlabel("spleen cell types served by that one hub")
    ax.set_xlim(0, 5.6)
    ax.set_title("(b) hub collapse — distinct hubs are far fewer than labels", fontsize=9)
    ax.grid(axis="y", visible=False)

    fig.suptitle("SenePy coverage of spleen — caller 2 is not the same estimator across arms",
                 fontsize=11, color=PAL.INK, y=1.0)
    fig.text(0.5, -0.045,
             f"SenePy 1.0.1 ships {j['n_hubs']} human hubs across {len(j['senepy_tissues'])} "
             f"tissues: {', '.join(j['senepy_tissues'])}. SPLEEN IS ABSENT; LIVER IS PRESENT, so "
             "M1 used tissue-matched hubs and H1 cannot. Tissue-matched: 0. Cross-tissue "
             f"surrogate: {j['n_usable_surrogate']}. No hub at all: {j['n_no_hub']} "
             "(cDC1, cDC2, pDC, lymphatic endothelium, erythroid, megakaryocyte, mesothelial) — "
             "those cells get no SenePy score, so audit test A3's sender counts are over a subset "
             "for this caller and over all cells for the others.",
             ha="center", fontsize=7.5, color=PAL.INK2, wrap=True)
    fig.tight_layout()
    save(fig, "figure_gs4_senepy_coverage")
    d.to_csv(f"{FIG}/figure_gs4_senepy_coverage_data.csv", index=False)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    for f in (gs1, gs2, gs3, gs4):
        print(f"{f.__name__} ...")
        f()
    print("done")
