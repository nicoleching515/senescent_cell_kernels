#!/usr/bin/env python3
"""Phase 10 — assemble the §17 two-arm comparison table FROM FILES.

Every row carries the absolute path of the file it came from and the exact filter used, so
the table is auditable line by line and nothing in the write-up has to be typed by hand.
Rows whose producer has not run yet are emitted with value `NOT_RUN` and the reason, rather
than being silently omitted.

  M1 = mouse liver, GSE310392, Xenium Prime Mouse 5K + 100 custom, 11 sections / 11 mice,
       6 §8-Test-3-admissible.  PRIMARY sender call `tierA_p95` (fine labels).
  H1 = human spleen, GSE326743, Xenium Prime Human 5K + 100 addon, 7 sections / 7 donors,
       all 7 analysed.  PRIMARY sender call `tierAmg_p95` == PREREG D-B's
       `tierA_merged_p95`; the frozen-literal fine-label `tierA_p95` is the sensitivity.

Usage: python3 code/two_arm_table.py
Writes results/phase10_h1/two_arm_table.csv  (long form: quantity, arm, panel, value, ...)
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")

R3 = "/workspace/results/phase3"
R5 = "/workspace/results/phase5"
RH = "/workspace/results/phase10_h1"
R9 = "/workspace/results/phase9_h1"
M1_INBAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
             "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
             "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
H1_SECS = ["SPLN07", "SPLN14", "SPLN21", "SPLN24", "SPLN30", "SPLN43", "SPLN44"]
ROWS = []


def add(quantity, arm, value, source, filt, panel="full", note=""):
    ROWS.append(dict(quantity=quantity, arm=arm, panel=panel,
                     value=value, source=source, filter=filt, note=note))


def missing(quantity, arm, source, reason, panel="full"):
    add(quantity, arm, "NOT_RUN", source, "", panel=panel, note=reason)


def _rep(mf, call, sections, stratum="all"):
    d = mf[mf.section.isin(sections) & (mf.call == call) & (mf.stratum == stratum)]
    return d, d[(d.beta_naive > 0) & (d.sf_base.notna()) & (d.beta_base_lo > 0)]


def main():
    # ---------------- platform / acquisition -------------------------------
    add("Platform / panel", "M1", "Xenium Prime Mouse 5K + 100 custom; 5,097 panel genes",
        f"{R9}/a8_panel_arithmetic.csv", "row mouse_panel")
    add("Platform / panel", "H1", "Xenium Prime Human 5K + 100 addon; 5,093 panel genes",
        f"{R9}/a8_panel_arithmetic.csv", "row human_panel")

    c = pd.read_csv("/workspace/results/composition_by_arm_timepoint.csv")
    add("Sections / donors", "M1", "11 / 11 mice; 6 admissible (§8 Test 3)",
        "/workspace/results/composition_by_arm_timepoint.csv",
        "len(c); c.mouse.nunique(); sasp_phase3.IN_BAND")
    add("Cells, total", "M1", int(c.n_cells.sum()),
        "/workspace/results/composition_by_arm_timepoint.csv", "c.n_cells.sum()")
    add("Cells, admissible sections", "M1",
        int(c[c.section.astype(str).str.split('_').str[0].isin(
            [s.split('_')[0] for s in M1_INBAND])].n_cells.sum()),
        "/workspace/results/composition_by_arm_timepoint.csv",
        "c[c.section in IN_BAND].n_cells.sum()")

    a1 = pd.read_csv(f"{R9}/a1_sections.csv")
    add("Sections / donors", "H1", "7 / 7 donors, ages 17-59, 4M/3F; all 7 analysed",
        f"{R9}/a1_sections.csv", "len(a1)")
    add("Cells, total", "H1", int(a1.n_cells.sum()), f"{R9}/a1_sections.csv",
        "a1.n_cells.sum()")
    add("Cells, admissible sections", "H1", int(a1.n_qc_pass.sum()),
        f"{R9}/a1_sections.csv", "a1.n_qc_pass.sum() (QC-passed; H1 has no Test-3 rule)")

    w3 = pd.read_csv(f"{R3}/window.csv")
    mn = w3.drop_duplicates("section").median_nn_um
    add("Median NN distance (um)", "M1", "%.2f - %.2f" % (mn.min(), mn.max()),
        f"{R3}/window.csv", "unique(section) -> median_nn_um")
    add("Median NN distance (um)", "H1",
        "%.2f - %.2f (all cells); %.2f - %.2f (QC-passed)"
        % (a1.median_nn_um_all.min(), a1.median_nn_um_all.max(),
           a1.median_nn_um_qc.min(), a1.median_nn_um_qc.max()),
        f"{R9}/a1_sections.csv", "median_nn_um_all / median_nn_um_qc",
        note="BELOW the frozen 7.0 um lambda-grid floor in every section; PREREG 3.1 "
             "forbids patching it")

    add("Transcript assignment rate", "M1", "88.27 % (section 7259 only)",
        "/workspace/logs/assignment_rate.log",
        "stdout of code/assignment_rate.py; no CSV exists for M1")
    ar = pd.read_csv(f"{R9}/a1_assignment_rate.csv")
    add("Transcript assignment rate", "H1",
        "%.2f - %.2f %% (3 sections)" % (ar.assigned_pct.min(), ar.assigned_pct.max()),
        f"{R9}/a1_assignment_rate.csv", "assigned_pct")

    # ---------------- the fits ---------------------------------------------
    mf3 = pd.read_csv(f"{R3}/main_fits.csv")
    m_all, m_rep = _rep(mf3, "tierA_p95", M1_INBAND)
    arms = [("M1", "tierA_p95", m_all, m_rep)]
    if os.path.exists(f"{RH}/main_fits.csv"):
        mfh = pd.read_csv(f"{RH}/main_fits.csv")
        for call in ("tierAmg_p95", "tierA_p95"):
            a, r = _rep(mfh, call, H1_SECS)
            arms.append(("H1 (%s)" % call, call, a, r))
    else:
        missing("naive amplitude", "H1", f"{RH}/main_fits.csv", "stage main not run")

    for label, call, d, rep in arms:
        src = f"{R3}/main_fits.csv" if label == "M1" else f"{RH}/main_fits.csv"
        flt = ("call=='%s' & stratum=='all' & section in %s"
               % (call, "IN_BAND(6)" if label == "M1" else "the 7 H1 sections"))
        add("Sender prevalence, primary call (%)", label,
            "%.3f - %.3f" % (100 * d.prevalence.min(), 100 * d.prevalence.max()), src, flt)
        add("Fits / reportable", label, "%d / %d (%.1f %%)"
            % (len(d), len(rep), 100 * len(rep) / max(len(d), 1)), src,
            flt + " ; reportable = beta_naive>0 & sf_base.notna() & beta_base_lo>0")
        add("Naive amplitude, median |b|/sd(y)", label,
            round(float((rep.beta_naive / rep.sd_y).median()), 4), src,
            flt + " ; median over reportable of beta_naive/sd_y")
        ctrl = rep.beta_n2n5n6 / rep.sd_y
        add("Controlled amplitude (N2+N5+N6), response-SD", label,
            "%.4f  IQR [%.4f, %.4f]" % (ctrl.median(), ctrl.quantile(.25),
                                        ctrl.quantile(.75)), src, flt)
        se = ((rep.beta_n2n5n6_hi - rep.beta_n2n5n6_lo) / (2 * 1.959964) / rep.sd_y)
        add("80 %-power detectable bound, response-SD", label,
            round(float(2.802 * se.median()), 4), src,
            flt + " ; 2.802 * median SE  (code/m1_headlines.py:30-32)")
        add("Controlled amplitude vs own bound", label,
            "%.4f vs %.4f -> %s" % (ctrl.median(), 2.802 * se.median(),
                                    "BELOW" if ctrl.median() < 2.802 * se.median()
                                    else "ABOVE"), src, flt)
        for k, lab in [("sf_n2", "SF, N2 matched decoy"),
                       ("sf_n5", "SF, N5 alone"),
                       ("sf_n6", "SF, N6 alone"),
                       ("sf_zon", "SF, anatomical covariate alone"),
                       ("sf_n6n5", "SF, N5+N6"),
                       ("sf_n2n5n6", "SF, N2+N5+N6 (PRIMARY OUTCOME)")]:
            v = rep[k]
            add(lab, label, "%.4f  IQR [%.4f, %.4f]"
                % (v.median(), v.quantile(.25), v.quantile(.75)), src,
                flt + " ; np.quantile over reportable fits")
        add("lambda-hat railed at a grid bound", label,
            "%.1f %% (%d at the %g um floor, %d at the %g um ceiling); median lam-hat %.2f um"
            % (100 * d.lam_railed.mean(),
               int((d.lam_naive <= d.lam_grid_lo + 1e-9).sum()),
               float(d.lam_grid_lo.iloc[0]),
               int((d.lam_naive >= d.lam_grid_hi - 1e-9).sum()),
               float(d.lam_grid_hi.iloc[0]), d.lam_naive.median()), src, flt)

    # ---------------- perturbation nulls ------------------------------------
    def perm_rows(label, res, call, sections, fn, cols, prefix):
        p = os.path.join(res, fn)
        base = mf3 if label == "M1" else (pd.read_csv(f"{RH}/main_fits.csv")
                                          if os.path.exists(f"{RH}/main_fits.csv") else None)
        if not os.path.exists(p) or base is None:
            for c in cols:
                missing("SF, %s%s" % (prefix, c.replace("_sf", "")), label, p,
                        "producer not run")
            return
        _, rep = _rep(base, call, sections)
        key = rep[["section", "celltype", "module"]].drop_duplicates()
        d = pd.read_csv(p)
        d = d[d.section.isin(sections) & (d.call == call)]
        for scope in (["full", "tile"] if "scope" in d.columns else [None]):
            dd = d[d.scope == scope] if scope else d
            dd = dd.merge(key, on=["section", "celltype", "module"])
            for c in cols:
                if c not in dd.columns:
                    continue
                v = dd[c].dropna()
                if not len(v):
                    continue
                nm = "SF, %s%s" % (prefix, c.replace("_sf", ""))
                if scope == "tile":
                    nm += " (tile scope)"
                add(nm, label, "%.4f  IQR [%.4f, %.4f]  n=%d"
                    % (v.median(), v.quantile(.25), v.quantile(.75), len(v)), p,
                    "call=='%s'%s ; merged onto the reportable fit key"
                    % (call, (" & scope=='%s'" % scope) if scope else ""))

    perm_rows("M1", R3, "tierA_p95", M1_INBAND, "perm_nulls.csv",
              ["N1_sf", "N3_sf", "N4_sf"], "")
    perm_rows("M1", R3, "tierA_p95", M1_INBAND, "perm_nulls_c1.csv",
              ["N3_tile_sf", "N4_tile_sf"], "C1 ")
    perm_rows("M1", R3, "tierA_p95", M1_INBAND, "perm_nulls_var.csv",
              ["N3_var_sf", "N4_var_sf"], "")
    for call in ("tierAmg_p95", "tierA_p95"):
        lab = "H1 (%s)" % call
        perm_rows(lab, RH, call, H1_SECS, "perm_nulls.csv",
                  ["N1_sf", "N3_sf", "N4_sf"], "")
        perm_rows(lab, RH, call, H1_SECS, "perm_nulls_c1.csv",
                  ["N3_tile_sf", "N4_tile_sf"], "C1 ")
        perm_rows(lab, RH, call, H1_SECS, "perm_nulls_var.csv",
                  ["N3_var_sf", "N4_var_sf"], "")

    # ---------------- geometric predictions (R3) ----------------------------
    for label, res in [("M1", R3), ("H1", RH)]:
        p = f"{res}/poisson_fits.csv"
        if os.path.exists(p):
            d = pd.read_csv(p)
            r = d[d.subset.str.contains("ALL sections", na=False)].iloc[0]
            add("Poisson identity (R3a): slope, r2", label,
                "%+.4f, r2 = %.4f" % (r.slope, r.r2), p,
                "subset == 'ALL sections x ALL sender definitions'")
        else:
            missing("Poisson identity (R3a)", label, p, "producer not run")
        p = f"{res}/null_destructiveness.csv"
        if os.path.exists(p):
            d = pd.read_csv(p)
            if label == "M1":
                d = d[d.section.isin(M1_INBAND)]
            g = d.groupby("null").frac_retaining_a_neighbour.median()
            add("R3c null destructiveness: senders keeping a real neighbour <=100 um",
                label, "; ".join("%s %.3f" % (k, v) for k, v in g.items()), p,
                "groupby('null').frac_retaining_a_neighbour.median()")
        else:
            missing("R3c null destructiveness", label, p, "producer not run")

    # ---------------- A7, the negative-control kernel -----------------------
    for label, p in [("M1", f"{R3}/a7_summary.csv"), ("H1", f"{R9}/a7_summary.csv")]:
        d = pd.read_csv(p)
        for resp, nm in [("neg_control_probe", "A7 PRIMARY: 40 negative control probes"),
                         ("all_controls", "A7 pooled negative-control features")]:
            for design, dn in [("base", "naive"), ("n6n5", "N6+N5")]:
                r = d[(d.response == resp) & (d.design == design)]
                if len(r):
                    r = r.iloc[0]
                    add("%s (%s)" % (nm, dn), label,
                        "%+.4f [%+.4f, %+.4f] p = %.4g" % (r.clustered_mean, r.clustered_lo,
                                                           r.clustered_hi, r.clustered_p),
                        p, "response=='%s' & design=='%s'" % (resp, design))

    # ---------------- composition-matched -----------------------------------
    for label, p, call in [("M1", f"{R3}/compmatch_reruns.csv", "tierA_p95"),
                           ("H1 (tierAmg_p95)", f"{RH}/compmatch_reruns_h1.csv",
                            "tierAmg_p95"),
                           ("H1 (tierA_p95)", f"{RH}/compmatch_reruns_h1.csv",
                            "tierA_p95")]:
        if not os.path.exists(p):
            missing("Composition-matched SF", label, p, "producer not run")
            continue
        d = pd.read_csv(p)
        d = d[(d.row_type == "summary") & (d.scope_kind == "pooled") & (d.call == call)]
        for v in ("comp", "full", "comp_adj", "type_adj", "typecomp_adj"):
            r = d[d.variant == v]
            if not len(r):
                continue
            r = r.iloc[0]
            add("Composition-matched SF [%s]%s" % (v, " (PRIMARY of the pair)"
                                                   if v == "typecomp_adj" else ""),
                label, "%.4f  95%% CI [%.4f, %.4f]  share removed %.1f %%  n_rep %d"
                % (r.median_sf_matched, r.median_sf_matched_lo, r.median_sf_matched_hi,
                   100 * r.median_comp_share, r.n_reportable), p,
                "row_type=='summary' & scope_kind=='pooled' & call=='%s' & variant=='%s'"
                % (call, v))

    d = pd.DataFrame(ROWS)
    os.makedirs(RH, exist_ok=True)
    d.to_csv(f"{RH}/two_arm_table.csv", index=False)
    print(d.to_string(index=False, max_colwidth=70))
    print("\n->", f"{RH}/two_arm_table.csv", d.shape)


if __name__ == "__main__":
    main()
