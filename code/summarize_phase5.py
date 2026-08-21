#!/usr/bin/env python3
"""Every number quoted in reports/CS_PHASE5.md, printed from the CSVs."""
from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P

R5 = "/workspace/results/phase5"
R3 = "/workspace/results/phase3"
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


def med(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.nan


def q(v, p):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.quantile(v, p)) if v.size else np.nan


# ===========================================================================
say("=" * 78)
say("T1  SUPERPOSITION vs NEAREST-SENDER  (Section 6.3)")
say("=" * 78)
ss = pd.read_csv(f"{R5}/super_section.csv")
sn = pd.read_csv(f"{R5}/super_nulls.csv")
ss["v"] = 1000 * ss.d_aic / ss.n
sn["v"] = 1000 * sn.d_aic / sn.n
syn = pd.read_csv("/workspace/results/phase1b/misspec.csv")
syn = syn[syn.fit_family == "exponential"]
pv = syn.pivot_table(index=["true_superposition", "regime", "rep"],
                     columns="fit_mode", values=["aic", "n"])
sv = ((pv["aic"]["superposition"] - pv["aic"]["nearest"])
      / pv["n"]["nearest"] * 1000).reset_index()
sv.columns = ["true_sup", "regime", "v"] if sv.shape[1] == 3 else \
    ["true_sup", "regime", "rep", "v"]

say("\nsynthetic calibration, dAIC(sup - near) per 1,000 cells:")
for (ts, rg), g in sv.groupby(["true_sup", "regime"]):
    say(f"   planted={'superposition' if ts else 'nearest      '} "
        f"{rg:11s}: median {med(g.v):+10.2f}   "
        f"correct verdict {(g.v < 0).mean() if ts else (g.v > 0).mean():.2f} "
        f"({len(g)} runs)")

c = ss[ss.design == "ctrl"]
say("\nreal data, per (section x receiver type x module) fit:")
for des, g in ss.groupby("design"):
    say(f"   {des:8s} n={len(g):4d}  dAIC/1k med {med(g.v):+8.4f} "
        f"[{q(g.v,.25):+.3f}, {q(g.v,.75):+.3f}]  "
        f"superposition wins {(g.d_aic < 0).mean():.3f}  "
        f"block-bootstrap win frac (median) "
        f"{med(g.boot_win_sup) if 'boot_win_sup' in g else np.nan:.3f}")
say("\nbound from the PAIRED spatial block bootstrap (same resample of blocks "
    "for both bases):")
for des, g in ss.groupby("design"):
    if des == "ctrl+N2":
        continue
    w = np.abs(np.c_[g.d_aic_per1k_lo, g.d_aic_per1k_hi]).max(1)
    say(f"   {des:6s} per-fit 95% CI on dAIC/1k: median "
        f"[{med(g.d_aic_per1k_lo):+.3f}, {med(g.d_aic_per1k_hi):+.3f}]; "
        f"excludes 0 in {((g.d_aic_per1k_lo>0)|(g.d_aic_per1k_hi<0)).mean():.3f}"
        f" of fits; |bound| median {med(w):.2f}, p90 {q(w,.9):.2f}")
    say(f"          paired bootstrap win fraction for superposition: median "
        f"{med(g.boot_win_sup):.3f}; decisive (>=0.95) in "
        f"{(g.boot_win_sup>=0.95).mean():.3f}, decisive for nearest (<=0.05) in "
        f"{(g.boot_win_sup<=0.05).mean():.3f}")
say(f"\namplitude bound (response-sd units, controlled): superposition median "
    f"{med(c.beta_sup_sd):+.4f} IQR [{q(c.beta_sup_sd,.25):+.4f}, "
    f"{q(c.beta_sup_sd,.75):+.4f}] p90 {q(c.beta_sup_sd,.9):.4f}; "
    f"nearest median {med(c.beta_near_sd):+.4f}")
nv = ss[ss.design == "naive"]
say(f"                          (naive): superposition median "
    f"{med(nv.beta_sup_sd):+.4f}; nearest median {med(nv.beta_near_sd):+.4f}")

say("\nnull stability (5 draws per section):")
for nm, g in sn.groupby("null"):
    say(f"   {nm} n={len(g):5d}  dAIC/1k med {med(g.v):+8.4f}  "
        f"superposition wins {(g.d_aic < 0).mean():.3f}   "
        f"per-draw " + str([round(float((h.d_aic < 0).mean()), 3)
                            for _, h in g.groupby('draw')]))
say(f"\nregressor correlation at the two lambda-hats: median "
    f"{med(c.r_regressors):.3f} [{q(c.r_regressors,.25):.3f}, "
    f"{q(c.r_regressors,.75):.3f}]")
say(f"lambda railed: nearest {c.railed_near.mean():.3f}  "
    f"superposition {c.railed_sup.mean():.3f};  median lambda "
    f"near {med(c.lam_near):.1f} um, sup {med(c.lam_sup):.1f} um")
say("\ndoes EITHER kernel beat the covariates-only model? (AIC)")
for des, g in ss.groupby("design"):
    if des == "ctrl+N2":
        continue
    say(f"   {des:6s}: nearest  median dAIC vs cov "
        f"{med(g.d_aic_near_vs_cov):+8.2f}, improves in "
        f"{(g.d_aic_near_vs_cov < 0).mean():.2f} of fits")
    say(f"   {des:6s}: superpos median dAIC vs cov "
        f"{med(g.d_aic_sup_vs_cov):+8.2f}, improves in "
        f"{(g.d_aic_sup_vs_cov < 0).mean():.2f} of fits")
for nm, g in sn.groupby("null"):
    say(f"   {nm:6s}: nearest {med(g.d_aic_near_vs_cov):+8.2f} "
        f"({(g.d_aic_near_vs_cov<0).mean():.2f}) | superpos "
        f"{med(g.d_aic_sup_vs_cov):+8.2f} "
        f"({(g.d_aic_sup_vs_cov<0).mean():.2f})")

h = pd.read_csv(f"{R5}/super_heldout.csv")
say("\nheld-out log-likelihood, leave-one-section-out (Section 24.6):")
for des, g in h.groupby("design"):
    p_ = g.pivot_table(index=["celltype", "module", "held_out"],
                       columns="basis", values=["ll_per_cell", "ll0_per_cell"])
    d = p_["ll_per_cell"]["sup"] - p_["ll_per_cell"]["near"]
    say(f"   {des:6s} folds={len(d)}  dLL/cell (sup-near) median "
        f"{med(d):+.3e}  superposition wins {(d > 0).mean():.3f}")
    for b in ("near", "sup"):
        dd = p_["ll_per_cell"][b] - p_["ll0_per_cell"][b]
        say(f"          {b:5s} vs covariates-only: median dLL/cell "
            f"{med(dd):+.3e}, beats it in {(dd > 0).mean():.3f} of folds")

sr_p = f"{R5}/se_ratio.csv"
if os.path.exists(sr_p):
    sr = pd.read_csv(sr_p)
    say("\niid vs spatial block bootstrap SE of the amplitude at the profiled "
        "lambda:")
    for (dz, b), g in sr.groupby(["design", "basis"]):
        say(f"   {dz:6s} {b:14s} n={len(g)}  SE ratio block/iid median "
            f"{med(g.se_ratio):.2f} [{q(g.se_ratio,.25):.2f}, "
            f"{q(g.se_ratio,.75):.2f}] p90 {q(g.se_ratio,.9):.2f}  |  "
            f"CI excludes 0: iid "
            f"{((g.ci_iid_lo>0)|(g.ci_iid_hi<0)).mean():.3f}  block "
            f"{((g.ci_blk_lo>0)|(g.ci_blk_hi<0)).mean():.3f}")

# ===========================================================================
say("\n" + "=" * 78)
say("T2  WINNER'S CURSE")
say("=" * 78)
m = pd.read_csv(f"{R3}/main_fits.csv")
d3 = m[(m.call == "tierA_p95") & (m.stratum == "all")
       & (m.section.isin(P.IN_BAND)) & (m.n >= 2000)]
sel3 = (d3.beta_naive > 0) & (d3.beta_base_lo > 0)
say(f"Phase 3 reference: {int(sel3.sum())}/{len(d3)} fits selected, "
    f"median SF(N5) {med(d3.loc[sel3,'sf_n5']):.4f}, "
    f"SF(N2+N5+N6) {med(d3.loc[sel3,'sf_n2n5n6']):.4f}")

for tag, f in (("0.50", "wc_crossfit.csv"), ("0.75", "wc_crossfit_f75.csv")):
    d = pd.read_csv(f"{R5}/{f}")
    rows = []
    for s, g in d.groupby("split"):
        s_ = g[g.selected_A == 1]
        rows.append(dict(sfA=med(s_.sf_n5_A), sfB=med(s_.sf_n5_B),
                         sfA_all=med(g.sf_n5_A), sfB_all=med(g.sf_n5_B),
                         sfA2=med(s_.sf_n2n5n6_A), sfB2=med(s_.sf_n2n5n6_B),
                         rate=g.selected_A.mean()))
    r = pd.DataFrame(rows)
    gap = r.sfA.mean() - r.sfB.mean()
    pl = r.sfA_all.mean() - r.sfB_all.mean()
    say(f"\nselection fraction {tag}  (selection rate {r.rate.mean():.3f}, "
        f"{len(r)} random splits)")
    say(f"   SF(N5) | selected on A : in-sample(A) {r.sfA.mean():.4f} "
        f"(sd {r.sfA.std():.4f})   held-out(B) {r.sfB.mean():.4f} "
        f"(sd {r.sfB.std():.4f})   gap {gap:+.4f}")
    say(f"   SF(N5) | ALL fits (placebo, no selection): "
        f"in-sample {r.sfA_all.mean():.4f}   held-out {r.sfB_all.mean():.4f}"
        f"   gap {pl:+.4f}")
    say(f"   placebo-corrected winner's curse = {gap - pl:+.4f}")
    say(f"   SF(N2+N5+N6) | selected: in {r.sfA2.mean():.4f}  "
        f"out {r.sfB2.mean():.4f}  gap {r.sfA2.mean()-r.sfB2.mean():+.4f}")

d = pd.read_csv(f"{R5}/wc_crossfit.csv")
u = d.drop_duplicates(["section", "celltype", "module"])
say(f"\nSF(N5) depends on sample size (no selection, median over fits):")
say(f"   100 blocks (full)  {med(u.sf_n5_full):.4f}")
r50 = med([med(g.sf_n5_A) for _, g in d.groupby('split')])
say(f"    50 blocks (half)  {r50:.4f}")
d75 = pd.read_csv(f"{R5}/wc_crossfit_f75.csv")
say(f"    25 blocks (quarter) "
    f"{med([med(g.sf_n5_B) for _, g in d75.groupby('split')]):.4f}")

sy = pd.read_csv(f"{R5}/wc_synthetic.csv")
sy = sy[sy.sweep == "ALL"]
say("\nsynthetic replication of the selection rule (Phase 1 sweep):")
for _, r in sy.iterrows():
    say(f"   beta_true={r.beta_true:.0f}  {r.subset:9s} n={int(r.n):5d}  "
        f"SF(N5) {r.sf_n5:.4f}  SF(N2) {r.sf_n2:.4f}  "
        f"SF(N2+N5) {r.sf_n2n5:.4f}   (selection rate {r.sel_rate:.3f})")

# ===========================================================================
say("\n" + "=" * 78)
say("T3  FIGURE 2a, STRATIFIED")
say("=" * 78)
amp = pd.read_csv("/workspace/figures/figure2a_amplitudes.csv", index_col=0)
say(amp.round(3).to_string())

# ===========================================================================
say("\n" + "=" * 78)
say("T4  KERNEL FAMILIES  (Section 6.2)")
say("=" * 78)
kf = pd.read_csv(f"{R5}/kernel_families.csv")
for des in ("naive", "ctrl"):
    g = kf[kf.design == des]
    say(f"\ndesign={des}  ({g.groupby(['section','celltype','module']).ngroups}"
        f" fits, stratified by receiver cell type)")
    t = g.groupby("family").agg(
        AIC_win=("d_aic_vs_best", lambda v: (v == 0).mean()),
        med_lam=("lam", "median"), med_d_half=("d_half", "median"),
        railed=("railed", "mean"),
        med_dAIC_vs_cov=("d_aic_vs_cov", "median"),
        beats_no_kernel=("d_aic_vs_cov", lambda v: (v < 0).mean()))
    say(t.round(3).to_string())
p_ = kf[kf.design == "ctrl"].pivot_table(
    index=["section", "celltype", "module"], columns="family", values="d_half")
say(f"\nd_half spread across families within a fit (max/min): "
    f"ctrl median {med(p_.max(1)/p_.min(1)):.2f}x")
p2 = kf[kf.design == "naive"].pivot_table(
    index=["section", "celltype", "module"], columns="family", values="d_half")
say(f"                                                       "
    f"naive median {med(p2.max(1)/p2.min(1)):.2f}x")

kh = pd.read_csv(f"{R5}/kernel_heldout.csv")
say("\nheld-out log-likelihood per family (Section 24.6):")
for des in ("naive", "ctrl"):
    g = kh[kh.design == des]
    pp = g.pivot_table(index=["celltype", "module", "held_out"],
                       columns="family", values="ll_per_cell")
    win = pp.idxmax(axis=1).value_counts(normalize=True)
    t = g.groupby("family").agg(
        med_dLL_vs_cov=("dll_vs_cov", "median"),
        beats_no_kernel=("dll_vs_cov", lambda v: (v > 0).mean()))
    t["heldout_win"] = win
    say(f"\ndesign={des}  ({len(pp)} leave-one-section-out folds)")
    say(t.round(3).fillna(0).to_string())

sw = pd.read_csv(f"{R5}/spline_window_check.csv")
say("\nwhy 'the spline wins 34/35' (CS Phase 2 §8) does not reproduce:")
for (c, st), g in sw.groupby(["cfg", "stratified"]):
    w = (g[g.family == "spline"].family
         == g[g.family == "spline"].best_family).mean()
    ws = (g[g.family == "step"].family
          == g[g.family == "step"].best_family).mean()
    nf = g.groupby(["section", "celltype", "module"]).ngroups
    say(f"   {c:24s} stratified={st}  ({nf:3d} fits, naive design): "
        f"spline wins {w:.3f}, step wins {ws:.3f}")

# ===========================================================================
say("\n" + "=" * 78)
say("T5  lambda_proximal vs lambda_downstream  (Section 6.4)")
say("=" * 78)
pdn = pd.read_csv(f"{R5}/proximal_vs_downstream.csv")
cols = ["celltype", "lam_prox_ctrl", "lam_down_ctrl", "ratio", "ratio_lo",
        "ratio_hi", "frac_prox_gt_down", "prox_sf_ctrl", "down_sf_ctrl"]
say(pdn[cols].round(3).to_string(index=False))
say(f"\ngrid ratio bounds: [{7/50:.3f}, {50/7:.3f}]; CIs reaching both bounds: "
    f"{int(((pdn.ratio_lo <= 0.15) & (pdn.ratio_hi >= 7.0)).sum())}/{len(pdn)}")

with open(f"{R5}/summary_phase5.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
print("\nwrote", f"{R5}/summary_phase5.txt")
