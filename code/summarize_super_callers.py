#!/usr/bin/env python3
"""T1b — the Phase 5 T1 superposition-vs-nearest verdict at every sender call.

Reads the tagged outputs of `run_phase5_super.py --call C --tag _C` (plus the
untagged tierA_p95 originals) and writes one tidy table,
`results/phase5/super_by_caller.csv`, with the five things the brief asks for:
median dAIC/1k with IQR, the PAIRED block-bootstrap CI, the observed/N1/N3
superposition win fractions, whether either kernel beats covariates-only, and
the leave-one-section-out held-out log-likelihood.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

R5 = "/workspace/results/phase5"
CALLS = [("tierA_p95", ""), ("cdkn1a_pos", "_cdkn1a_pos"),
         ("senepy_p95", "_senepy_p95")]
OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


def med(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.nan


def frac(v, sign=-1):
    """Fraction of finite entries on the requested side of 0; NaN if none."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if not v.size:
        return np.nan
    return float((v < 0).mean() if sign < 0 else (v > 0).mean())


def q(v, p):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return float(np.quantile(v, p)) if v.size else np.nan


rows = []
for call, tag in CALLS:
    fs = f"{R5}/super_section{tag}.csv"
    if not os.path.exists(fs):
        say(f"!! missing {fs}"); continue
    ss = pd.read_csv(fs)
    ss["v"] = 1000 * ss.d_aic / ss.n
    for des, g in ss.groupby("design"):
        w = np.abs(np.c_[g.d_aic_per1k_lo, g.d_aic_per1k_hi]).max(1)
        r = dict(call=call, design=des, n_fits=len(g),
                 med_senders=med(g.n_senders), med_n=med(g.n),
                 daic1k_med=med(g.v), daic1k_q25=q(g.v, .25),
                 daic1k_q75=q(g.v, .75),
                 sup_win=float((g.d_aic < 0).mean()),
                 boot_win_med=med(g.boot_win_sup) if "boot_win_sup" in g else np.nan,
                 ci_lo_med=med(g.d_aic_per1k_lo), ci_hi_med=med(g.d_aic_per1k_hi),
                 ci_excl0=float(((g.d_aic_per1k_lo > 0) |
                                 (g.d_aic_per1k_hi < 0)).mean()),
                 bound_med=med(w), bound_p90=q(w, .9),
                 decisive_sup=float((g.boot_win_sup >= .95).mean())
                 if "boot_win_sup" in g else np.nan,
                 decisive_near=float((g.boot_win_sup <= .05).mean())
                 if "boot_win_sup" in g else np.nan,
                 near_vs_cov_med=med(g.d_aic_near_vs_cov),
                 near_vs_cov_frac=frac(g.d_aic_near_vs_cov),
                 sup_vs_cov_med=med(g.d_aic_sup_vs_cov),
                 sup_vs_cov_frac=frac(g.d_aic_sup_vs_cov),
                 beta_sup_sd=med(g.beta_sup_sd) if "beta_sup_sd" in g else np.nan,
                 beta_near_sd=med(g.beta_near_sd) if "beta_near_sd" in g else np.nan,
                 r_regressors=med(g.r_regressors) if "r_regressors" in g else np.nan,
                 railed_near=float(g.railed_near.mean()) if "railed_near" in g else np.nan,
                 railed_sup=float(g.railed_sup.mean()) if "railed_sup" in g else np.nan,
                 lam_near_med=med(g.lam_near), lam_sup_med=med(g.lam_sup))
        rows.append(r)
    fn = f"{R5}/super_nulls{tag}.csv"
    if os.path.exists(fn):
        sn = pd.read_csv(fn)
        sn["v"] = 1000 * sn.d_aic / sn.n
        for nm, g in sn.groupby("null"):
            rows.append(dict(
                call=call, design=nm, n_fits=len(g), med_n=med(g.n),
                daic1k_med=med(g.v), daic1k_q25=q(g.v, .25),
                daic1k_q75=q(g.v, .75),
                sup_win=float((g.d_aic < 0).mean()),
                near_vs_cov_med=med(g.d_aic_near_vs_cov),
                near_vs_cov_frac=frac(g.d_aic_near_vs_cov),
                sup_vs_cov_med=med(g.d_aic_sup_vs_cov),
                sup_vs_cov_frac=frac(g.d_aic_sup_vs_cov),
                per_draw=str([round(float((h.d_aic < 0).mean()), 3)
                              for _, h in g.groupby("draw")]),
                lam_near_med=med(g.lam_near), lam_sup_med=med(g.lam_sup)))
    fh = f"{R5}/super_heldout{tag}.csv"
    if os.path.exists(fh):
        h = pd.read_csv(fh)
        for des, g in h.groupby("design"):
            p_ = g.pivot_table(index=["celltype", "module", "held_out"],
                               columns="basis",
                               values=["ll_per_cell", "ll0_per_cell"])
            d = p_["ll_per_cell"]["sup"] - p_["ll_per_cell"]["near"]
            dn = p_["ll_per_cell"]["near"] - p_["ll0_per_cell"]["near"]
            dsu = p_["ll_per_cell"]["sup"] - p_["ll0_per_cell"]["sup"]
            rows.append(dict(call=call, design=f"heldout_{des}", n_fits=len(d),
                             dll_med=med(d), sup_win=float((d > 0).mean()),
                             near_vs_cov_med=med(dn),
                             near_vs_cov_frac=float((dn > 0).mean()),
                             sup_vs_cov_med=med(dsu),
                             sup_vs_cov_frac=float((dsu > 0).mean())))

df = pd.DataFrame(rows)
df.to_csv(f"{R5}/super_by_caller.csv", index=False)

# ---------------------------------------------------------------------------
# two diagnostics quoted in the report, persisted so they are reproducible:
#   (a) is the cdkn1a_pos effect gene-set circularity?  Cdkn1a is a member of
#       4 of the 7 Tier B response modules; split on that.
#   (b) where superposition beats covariates-only, is lambda_hat at the ceiling
#       (i.e. is it acting as a regional density covariate, not a kernel)?
# ---------------------------------------------------------------------------
WITH_CDKN1A = {"downstream_arrest", "interferon_response",
               "secondary_senescence", "tnfa_nfkb_proximal"}
drows = []
for call, tag in CALLS:
    f = f"{R5}/super_section{tag}.csv"
    if not os.path.exists(f):
        continue
    c = pd.read_csv(f)
    c = c[c.design == "ctrl"].copy()
    c["v"] = 1000 * c.d_aic / c.n
    for gname, g in c.groupby(np.where(c.module.isin(WITH_CDKN1A),
                                       "Cdkn1a_in_module", "no_Cdkn1a")):
        drows.append(dict(call=call, split=gname, n=len(g), daic1k_med=med(g.v),
                          sup_win=float((g.d_aic < 0).mean()),
                          sup_vs_cov_med=med(g.d_aic_sup_vs_cov),
                          sup_vs_cov_frac=float((g.d_aic_sup_vs_cov < 0).mean())))
    w = c[c.d_aic_sup_vs_cov < 0]
    drows.append(dict(call=call, split="ALL", n=len(c), daic1k_med=med(c.v),
                      sup_win=float((c.d_aic < 0).mean()),
                      sup_vs_cov_med=med(c.d_aic_sup_vs_cov),
                      sup_vs_cov_frac=float((c.d_aic_sup_vs_cov < 0).mean()),
                      lam_sup_ceiling=float((c.lam_sup > 49.9).mean()),
                      lam_sup_floor=float((c.lam_sup < 7.01).mean()),
                      lam_sup_ceiling_where_sup_beats_cov=float((w.lam_sup > 49.9).mean())
                      if len(w) else np.nan,
                      lam_sup_med_where_sup_beats_cov=med(w.lam_sup) if len(w) else np.nan))
pd.DataFrame(drows).to_csv(f"{R5}/super_by_caller_diag.csv", index=False)

say("=" * 78)
say("T1b  SUPERPOSITION vs NEAREST ACROSS SENDER DEFINITIONS")
say("=" * 78)
say("\nsynthetic calibration (unchanged): planted nearest +15.2/+17.0 per 1k, "
    "planted superposition -57.2/-250.9, 15/15 correct.\n")
for call, _ in CALLS:
    d = df[df.call == call]
    if not len(d):
        continue
    c = d[d.design == "ctrl"]
    if not len(c):
        continue
    c = c.iloc[0]
    say(f"--- {call}  ({int(c.n_fits)} fits, median {int(c.med_senders)} "
        f"senders/section, median {int(c.med_n)} receivers/fit)")
    for des in ("naive", "ctrl", "ctrl+N2"):
        g = d[d.design == des]
        if not len(g):
            continue
        g = g.iloc[0]
        say(f"   {des:8s} dAIC/1k med {g.daic1k_med:+8.4f} "
            f"[{g.daic1k_q25:+.3f}, {g.daic1k_q75:+.3f}]  sup wins "
            f"{g.sup_win:.3f}  paired boot win "
            f"{g.boot_win_med:.3f}" if np.isfinite(g.boot_win_med) else
            f"   {des:8s} dAIC/1k med {g.daic1k_med:+8.4f} "
            f"[{g.daic1k_q25:+.3f}, {g.daic1k_q75:+.3f}]  sup wins "
            f"{g.sup_win:.3f}")
    say(f"   PAIRED bootstrap CI on dAIC/1k (ctrl): median "
        f"[{c.ci_lo_med:+.3f}, {c.ci_hi_med:+.3f}]; excludes 0 in "
        f"{c.ci_excl0:.3f} of fits; |bound| med {c.bound_med:.2f} p90 "
        f"{c.bound_p90:.2f}; decisive for sup {c.decisive_sup:.3f}, "
        f"for near {c.decisive_near:.3f}")
    for nm in ("N1", "N3"):
        g = d[d.design == nm]
        if len(g):
            g = g.iloc[0]
            say(f"   {nm} null   dAIC/1k med {g.daic1k_med:+8.4f}  sup wins "
                f"{g.sup_win:.3f}  per-draw {g.per_draw}")
    say(f"   vs covariates-only (ctrl, AIC): nearest {c.near_vs_cov_med:+7.2f} "
        f"(improves {c.near_vs_cov_frac:.2f}) | superposition "
        f"{c.sup_vs_cov_med:+7.2f} (improves {c.sup_vs_cov_frac:.2f})")
    for nm in ("N1", "N3"):
        g = d[d.design == nm]
        if len(g):
            g = g.iloc[0]
            say(f"      {nm}: nearest {g.near_vs_cov_med:+7.2f} "
                f"({g.near_vs_cov_frac:.2f}) | superposition "
                f"{g.sup_vs_cov_med:+7.2f} ({g.sup_vs_cov_frac:.2f})")
    for des in ("ctrl", "naive"):
        g = d[d.design == f"heldout_{des}"]
        if len(g):
            g = g.iloc[0]
            say(f"   held-out LL {des:5s} folds={int(g.n_fits)}  dLL/cell "
                f"(sup-near) {g.dll_med:+.3e}  sup wins {g.sup_win:.3f}  |  "
                f"near beats cov {g.near_vs_cov_frac:.3f}, sup beats cov "
                f"{g.sup_vs_cov_frac:.3f}")
    say(f"   amplitude (ctrl, response-sd): sup {c.beta_sup_sd:+.4f}, near "
        f"{c.beta_near_sd:+.4f}; regressor r {c.r_regressors:.3f}; railed "
        f"near {c.railed_near:.3f} sup {c.railed_sup:.3f}; median lambda "
        f"near {c.lam_near_med:.1f} sup {c.lam_sup_med:.1f} um")
    say("")

with open(f"{R5}/summary_super_callers.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
say("gene-overlap and lambda-railing diagnostics: super_by_caller_diag.csv")
say(pd.DataFrame(drows).round(4).to_string(index=False))
with open(f"{R5}/summary_super_callers.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
print(f"\nwrote {R5}/super_by_caller.csv, super_by_caller_diag.csv and "
      f"summary_super_callers.txt")
