#!/usr/bin/env python3
"""Phase 10 — re-read every headline number of reports/CS_PHASE10_TWO_ARM.md from its file.

The analogue of Phase 9's `code/h1_verify_report.py`, and it exists for the same reason: this
project has repeatedly found its own written numbers not to match the files behind them
(`reports/AUDIT_PHASE8_FACTCHECK.md`).  Each check names the quantity, the file, the filter,
the value written in the report, and the value re-derived now.  Exits non-zero on any drift.

  python3 code/phase10_verify_report.py
"""
import json, os, sys
import numpy as np, pandas as pd

R3 = "/workspace/results/phase3"
R5 = "/workspace/results/phase5"
R9 = "/workspace/results/phase9_h1"
RH = "/workspace/results/phase10_h1"
RI = "/workspace/results/phase10_h1_iso"
M1_INBAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
             "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
             "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
H1 = ["SPLN07", "SPLN14", "SPLN21", "SPLN24", "SPLN30", "SPLN43", "SPLN44"]
CHECKS, FAILED = [], 0


def chk(name, written, got, tol=5e-4, src=""):
    global FAILED
    ok = (written == got) if isinstance(written, str) else abs(written - got) <= tol
    FAILED += (not ok)
    CHECKS.append((name, written, got, "OK" if ok else "**DRIFT**", src))


def rep(path, call, sections, stratum="all"):
    d = pd.read_csv(path)
    d = d[d.section.isin(sections) & (d.call == call) & (d.stratum == stratum)]
    return d, d[(d.beta_naive > 0) & (d.sf_base.notna()) & (d.beta_base_lo > 0)]


def perm(path, col, call, sections, key, scope=None):
    d = pd.read_csv(path)
    d = d[d.section.isin(sections) & (d.call == call)]
    if scope and "scope" in d.columns:
        d = d[d.scope == scope]
    d = d.merge(key, on=["section", "celltype", "module"])
    return float(d[col].dropna().median())


def main():
    # ---- the primary outcome and the bound, both arms, both H1 calls ----
    for lab, path, call, secs, w in [
            ("M1", f"{R3}/main_fits.csv", "tierA_p95", M1_INBAND,
             dict(nfit=315, nrep=153, naive=0.3288, ctrl=0.0288, bound=0.1833,
                  sf=0.0885, railed=0.600, lam=14.73)),
            ("H1 tierAmg_p95", f"{RH}/main_fits.csv", "tierAmg_p95", H1,
             dict(nfit=343, nrep=58, naive=0.1843, ctrl=0.0307, bound=0.1094,
                  sf=0.1582, railed=0.700, lam=38.16)),
            ("H1 tierA_p95", f"{RH}/main_fits.csv", "tierA_p95", H1,
             dict(nfit=343, nrep=61, naive=0.1807, ctrl=0.0345, bound=0.1083,
                  sf=0.1630, railed=0.741, lam=40.49))]:
        d, r = rep(path, call, secs)
        se = ((r.beta_n2n5n6_hi - r.beta_n2n5n6_lo) / (2 * 1.959964) / r.sd_y).median()
        chk(f"{lab}: n fits", w["nfit"], len(d), 0, path)
        chk(f"{lab}: n reportable", w["nrep"], len(r), 0, path)
        chk(f"{lab}: naive amplitude", w["naive"], float((r.beta_naive / r.sd_y).median()),
            src=path)
        chk(f"{lab}: controlled amplitude N2+N5+N6", w["ctrl"],
            float((r.beta_n2n5n6 / r.sd_y).median()), src=path)
        chk(f"{lab}: 80%-power bound", w["bound"], float(2.802 * se), src=path)
        chk(f"{lab}: SF N2+N5+N6", w["sf"], float(r.sf_n2n5n6.median()), src=path)
        chk(f"{lab}: lambda railed fraction", w["railed"], float(d.lam_railed.mean()),
            tol=6e-4, src=path)
        chk(f"{lab}: median lambda-hat (um)", w["lam"], float(d.lam_naive.median()),
            tol=6e-3, src=path)
        n_over = 0
        for m, g in r.groupby("module"):
            if float((g.beta_n2n5n6 / g.sd_y).median()) > 2.802 * se:
                n_over += 1
        chk(f"{lab}: Tier B modules above the arm's own bound", 0, n_over, 0, path)

    # ---- the perturbation nulls ----
    _, rm = rep(f"{R3}/main_fits.csv", "tierA_p95", M1_INBAND)
    km = rm[["section", "celltype", "module"]].drop_duplicates()
    _, rh = rep(f"{RH}/main_fits.csv", "tierAmg_p95", H1)
    kh = rh[["section", "celltype", "module"]].drop_duplicates()
    for lab, res, call, secs, key, w in [
            ("M1", R3, "tierA_p95", M1_INBAND, km,
             dict(N1=0.7074, N3=0.9989, N4=0.9472, N3var=0.9960, N4var=0.9851,
                  N3tile=0.9706, N4tile=0.9243)),
            ("H1", RH, "tierAmg_p95", H1, kh,
             dict(N1=0.3319, N3=0.9929, N4=0.9561, N3var=1.0000, N4var=0.9740,
                  N3tile=0.9507, N4tile=0.9063))]:
        chk(f"{lab}: SF N1", w["N1"], perm(f"{res}/perm_nulls.csv", "N1_sf", call, secs, key),
            src=f"{res}/perm_nulls.csv")
        chk(f"{lab}: SF N3 (published bbox)", w["N3"],
            perm(f"{res}/perm_nulls.csv", "N3_sf", call, secs, key), src=f"{res}/perm_nulls.csv")
        chk(f"{lab}: SF N4 (published bbox)", w["N4"],
            perm(f"{res}/perm_nulls.csv", "N4_sf", call, secs, key), src=f"{res}/perm_nulls.csv")
        chk(f"{lab}: SF N3-var (PRIMARY corrected)", w["N3var"],
            perm(f"{res}/perm_nulls_var.csv", "N3_var_sf", call, secs, key, "full"),
            src=f"{res}/perm_nulls_var.csv")
        chk(f"{lab}: SF N4-var", w["N4var"],
            perm(f"{res}/perm_nulls_var.csv", "N4_var_sf", call, secs, key, "full"),
            src=f"{res}/perm_nulls_var.csv")
        chk(f"{lab}: SF N3-tile (tile scope)", w["N3tile"],
            perm(f"{res}/perm_nulls_c1.csv", "N3_tile_sf", call, secs, key, "tile"),
            src=f"{res}/perm_nulls_c1.csv")
        chk(f"{lab}: SF N4-tile (tile scope)", w["N4tile"],
            perm(f"{res}/perm_nulls_c1.csv", "N4_tile_sf", call, secs, key, "tile"),
            src=f"{res}/perm_nulls_c1.csv")

    # ---- composition-matched, both arms ----
    for lab, path, call, w in [
            ("M1", f"{R3}/compmatch_reruns.csv", "tierA_p95",
             dict(comp=0.9837, typecomp=0.1461, typeadj=0.3414)),
            ("H1", f"{RH}/compmatch_reruns_h1.csv", "tierAmg_p95",
             dict(comp=0.9988, typecomp=0.0698, typeadj=0.2958))]:
        d = pd.read_csv(path)
        d = d[(d.row_type == "summary") & (d.scope_kind == "pooled") & (d.call == call)]
        for v, k in [("comp", "comp"), ("typecomp_adj", "typecomp"), ("type_adj", "typeadj")]:
            chk(f"{lab}: composition-matched SF [{v}]", w[k],
                float(d[d.variant == v].median_sf_matched.iloc[0]), src=path)

    # ---- A7, both arms ----
    for lab, path, w in [("M1", f"{R3}/a7_summary.csv", dict(base=-0.0225, n6n5=0.0015)),
                         ("H1", f"{R9}/a7_summary.csv", dict(base=-0.0118, n6n5=-0.0028))]:
        d = pd.read_csv(path)
        for design in ("base", "n6n5"):
            r = d[(d.response == "neg_control_probe") & (d.design == design)].iloc[0]
            chk(f"{lab}: A7 primary (40 probes), {design}", w[design],
                float(r.clustered_mean), tol=1e-4, src=path)

    # ---- R3 geometric predictions ----
    for lab, path, w in [("M1", f"{R3}/poisson_fits.csv", 0.9843),
                         ("H1", f"{RH}/poisson_fits.csv", 0.9783)]:
        d = pd.read_csv(path)
        r = d[d.subset.str.contains("ALL sections", na=False)].iloc[0]
        chk(f"{lab}: Poisson identity r2 (R3a)", w, float(r.r2), src=path)

    # ---- §8 predictions ----
    a = pd.read_csv(f"{R9}/deepscence_anchor_h1.csv")
    chk("P-ii: H1 CDKN1A fold-split stability, min over 7 sections", 1.0,
        float(a.stab_cdkn1a.min()), 0, f"{R9}/deepscence_anchor_h1.csv")
    chk("P-ii: H1 depth-partialled rho, min over 7", 0.1911,
        float(a.rho_partial_cdkn1a.min()), src=f"{R9}/deepscence_anchor_h1.csv")
    g = json.load(open("/workspace/results/phase7_jobA/gate_result_human.json"))["corescence"]
    chk("P-iv: native CoreScence circularity", 0.8788, float(g["frozen_frac"]),
        src="results/phase7_jobA/gate_result_human.json")
    p = pd.read_csv(f"{R9}/caller_agreement_pooled.csv")
    chk("P-iii: pooled TierA x |DeepScence| ratio", 1.102,
        float(p[(p.A == "tierA_score") & (p.B == "abs_deepscence_score")].pooled_ratio.iloc[0]),
        src=f"{R9}/caller_agreement_pooled.csv")
    chk("D-C: DeepScence x CDKN1A+ pooled (circular)", 6.436,
        float(p[(p.A == "deepscence_score") & (p.B == "cdkn1a_counts")].pooled_ratio.iloc[0]),
        src=f"{R9}/caller_agreement_pooled.csv")
    dd = pd.read_csv(f"{R9}/d2_depth.csv")
    dd = dd[dd.scope == "20,000-cell panel"]
    chk("P-vi: sections with delta_rho <= 0, of 7", 6, int((dd.delta_rho <= 0).sum()), 0,
        f"{R9}/d2_depth.csv")
    st = pd.read_csv(f"{R9}/d2_stability.csv")
    r = st[st.config == "denoise=False, FULL section"].iloc[0]
    chk("H10: full-section seed 0 vs 1 Pearson r", 0.3719, float(r.pearson_r), tol=1e-4,
        src=f"{R9}/d2_stability.csv")
    chk("H10: full-section seed 0 vs 1 top-5% Jaccard", 0.2107, float(r.top5_jaccard),
        tol=1e-4, src=f"{R9}/d2_stability.csv")

    # ---- the five-seed DeepScence consensus (D-A) ----
    if os.path.exists(f"{RH}/deepscence_consensus_coverage.csv"):
        cc = pd.read_csv(f"{RH}/deepscence_consensus_coverage.csv")
        cc = cc[cc.status == "OK"]
        chk("D-A: sections with a complete 5-seed panel", 7,
            int((cc.n_seeds == 5).sum()), 0, f"{RH}/deepscence_consensus_coverage.csv")
        chk("D-A: seeds whose score was SIGN-INVERTED", 1,
            int(cc.n_seeds_sign_flipped.sum()), 0, f"{RH}/deepscence_consensus_coverage.csv")
        chk("D-A: min mean pairwise top-5% Jaccard across seeds", 0.157,
            float(cc.jaccard_top5_mean.min()), tol=1e-3,
            src=f"{RH}/deepscence_consensus_coverage.csv")
        chk("D-A: max mean pairwise top-5% Jaccard across seeds", 0.698,
            float(cc.jaccard_top5_mean.max()), tol=1e-3,
            src=f"{RH}/deepscence_consensus_coverage.csv")
        sg = pd.read_csv(f"{RH}/deepscence_consensus_sign.csv")
        f = sg[sg.flipped]
        chk("D-A: the inverted run's anchor rho (SPLN07, seed 20260903)", -0.1492,
            float(f.rho_partial_cdkn1a.iloc[0]), tol=1e-4,
            src=f"{RH}/deepscence_consensus_sign.csv")

    # ---- the ortholog-intersected panel ----
    if os.path.exists(f"{RI}/iso_vs_full_headlines.csv"):
        d = pd.read_csv(f"{RI}/iso_vs_full_headlines.csv")
        for arm, call, panel, w in [("m1", "tierA_p95", "intersected", 0.1054),
                                    ("h1", "tierAmg_p95", "intersected", 0.1785),
                                    ("h1", "tierA_p95", "intersected", 0.2236)]:
            r = d[(d.arm == arm) & (d.call == call) & (d.panel == panel)].iloc[0]
            chk(f"ISO: {arm} {call} SF N2+N5+N6 on the intersected panel", w,
                float(r["SF_N2+N5+N6"]), src=f"{RI}/iso_vs_full_headlines.csv")

    # ---- Phase 5 ----
    for lab, path in [("M1", f"{R5}/kernel_families.csv"), ("H1", f"{RH}/kernel_families.csv")]:
        d = pd.read_csv(path); d = d[d.design == "ctrl"]
        w = 0.898 if lab == "M1" else 0.956
        chk(f"{lab}: step wins the AIC, controlled design", w,
            float((d[d.family == "step"].family == d[d.family == "step"].best_family).mean()),
            tol=1e-3, src=path)
    for lab, path in [("M1", f"{R5}/super_section.csv"), ("H1", f"{RH}/super_section.csv")]:
        d = pd.read_csv(path); d = d[d.design == "ctrl"]
        chk(f"{lab}: paired block-bootstrap superposition win fraction", 0.730,
            float(d.boot_win_sup.median()), tol=1e-3, src=path)

    w = max(len(c[0]) for c in CHECKS)
    print("%-*s  %12s  %12s  %s" % (w, "quantity", "in report", "re-derived", "status"))
    for n, a_, b_, s, src in CHECKS:
        fa = a_ if isinstance(a_, str) else ("%.4f" % a_ if isinstance(a_, float) else a_)
        fb = b_ if isinstance(b_, str) else ("%.4f" % b_ if isinstance(b_, float) else b_)
        print("%-*s  %12s  %12s  %s" % (w, n, fa, fb, s))
    print("\n%d checks, %d failed" % (len(CHECKS), FAILED))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
