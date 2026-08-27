#!/usr/bin/env python3
"""Phase 8 / 8.7 -- extract the headline vector the Section 17 two-arm table and
`reports/CORRECTIONS.md` quote, from ONE results directory, so the pre-C6 and
post-C6 trees can be compared by running the same code twice.

  python3 m1_headlines.py /workspace/results/phase3 /workspace/results/phase5
"""
import sys, os, json
import numpy as np, pandas as pd

R3 = sys.argv[1] if len(sys.argv) > 1 else "/workspace/results/phase3"
R5 = sys.argv[2] if len(sys.argv) > 2 else "/workspace/results/phase5"
IN_BAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
           "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
           "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
CALL = "tierA_p95"
out = {}

mf = pd.read_csv(f"{R3}/main_fits.csv")
mf = mf[mf.section.isin(IN_BAND) & (mf.call == CALL) & (mf.stratum == "all")]
rep = mf[(mf.beta_naive > 0) & (mf.sf_base.notna()) & (mf.beta_base_lo > 0)]
out["n_fits"] = int(len(mf))
out["n_reportable"] = int(len(rep))
out["sender_prevalence_min"] = float(mf.prevalence.min())
out["sender_prevalence_max"] = float(mf.prevalence.max())
out["naive_amp_med"] = float((rep.beta_naive / rep.sd_y).median())
ctrl = rep.beta_n2n5n6 / rep.sd_y
out["ctrl_amp_med"] = float(ctrl.median())
out["ctrl_amp_iqr"] = [float(ctrl.quantile(.25)), float(ctrl.quantile(.75))]
se = ((rep.beta_n2n5n6_hi - rep.beta_n2n5n6_lo) / (2 * 1.959964) / rep.sd_y)
out["ctrl_amp_se_med"] = float(se.median())
out["power80_bound"] = float(2.802 * se.median())   # 1.96 + 0.842
for k, lab in [("sf_n2", "N2"), ("sf_n5", "N5"), ("sf_n6", "N6"),
               ("sf_zon", "zon"), ("sf_n6n5", "N5+N6"),
               ("sf_n2n5n6", "N2+N5+N6")]:
    out[f"SF_{lab}"] = float(rep[k].median())
# the published bracket on this row is the IQR over reportable fits, not a CI
out["SF_N2+N5+N6_iqr"] = [float(rep.sf_n2n5n6.quantile(.25)),
                          float(rep.sf_n2n5n6.quantile(.75))]
out["lam_railed_frac"] = float(mf.lam_railed.mean())
out["ctrl_pos_and_sig"] = int(((rep.beta_n2n5n6 > 0) &
                               (rep.beta_n2n5n6_lo > 0)).sum())

def _perm(f, cols):
    p = os.path.join(R3, f)
    if not os.path.exists(p):
        return
    d = pd.read_csv(p)
    d = d[d.section.isin(IN_BAND) & (d.call == CALL)]
    if "scope" in d.columns:
        d = d[d.scope == "full"] if "full" in set(d.scope) else d
    key = rep[["section", "celltype", "module"]].drop_duplicates()
    d = d.merge(key, on=["section", "celltype", "module"])
    for c in cols:
        if c in d.columns:
            out[c.replace("_sf", "")] = float(d[c].median())
    return d

_perm("perm_nulls.csv", ["N1_sf", "N3_sf", "N4_sf", "N1_full_sf"])
p = os.path.join(R3, "perm_nulls_c1.csv")
if os.path.exists(p):
    d = pd.read_csv(p)
    d = d[d.section.isin(IN_BAND) & (d.call == CALL)]
    key = rep[["section", "celltype", "module"]].drop_duplicates()
    d = d.merge(key, on=["section", "celltype", "module"])
    for nm in ("N3_orig", "N3_tile", "N3_occ", "N3_occ15", "N3_swap",
               "N3_snap", "N4_orig", "N4_tile", "N4_occ", "N4_occ15",
               "N4_swap"):
        c = f"{nm}_sf"
        if c in d.columns:
            v = d[c].dropna()
            if len(v):
                out[f"c1_{nm}"] = float(v.median())
        c = f"{nm}_full_sf"
        if c in d.columns:
            v = d[c].dropna()
            if len(v):
                out[f"c1_{nm}_fullsf"] = float(v.median())

n8 = [f"{R3}/n8_scrambled_{s}.csv" for s in IN_BAND]
n8 = [f for f in n8 if os.path.exists(f)]
if n8:
    d = pd.concat([pd.read_csv(f) for f in n8], ignore_index=True)
    c = "sf_n8" if "sf_n8" in d.columns else None
    if c:
        out["SF_N8"] = float(d[c].median())

p = f"{R3}/poisson_fits.csv"
if os.path.exists(p):
    d = pd.read_csv(p)
    r = d[d.subset.str.contains("ALL sections", na=False)] if "subset" in d else d
    if len(r):
        out["poisson_slope"] = float(r.slope.iloc[0])
        out["poisson_r2"] = float(r.r2.iloc[0])

p = f"{R3}/attribution.csv"
if os.path.exists(p):
    d = pd.read_csv(p)
    d = d[(d.band == "in_band") & (d.call == CALL)]
    out["attr_sf_comp_med"] = float(d.sf_comp.median())
    out["attr_sf_tech_med"] = float(d.sf_tech.median())
    out["attr_sf_dens_med"] = float(d.sf_dens.median())
    out["attr_sf_anat_med"] = float(d.sf_anat.median())

p = f"{R5}/kernel_families.csv"
if os.path.exists(p):
    d = pd.read_csv(p)
    if "design" in d.columns:
        for dg in sorted(set(d.design.dropna())):
            s = d[d.design == dg]
            if "aic_win" in s.columns:
                out[f"p5_{dg}_step_win"] = float(
                    (s[s.family == "step"].aic_win.mean()
                     if "family" in s.columns else np.nan))

p = f"{R3}/ripley.csv"
if os.path.exists(p):
    d = pd.read_csv(p)
    d = d[d.call == CALL]
    out["ripley_ratio_med"] = float(d.ripley_ratio.median())

print(json.dumps(out, indent=1, sort_keys=True))
