#!/usr/bin/env python3
"""
Phase 5 T2 — how much of SF(N5) = 0.084 is WINNER'S CURSE?

Phase 3's headline is computed on the 160 of 315 fits that pass a NAIVE
significance screen (beta_naive > 0 AND spatial block-bootstrap CI excluding
zero).  Selecting on the naive amplitude picks fits whose beta_naive is
upward-biased by noise, and the surviving fraction SF = beta_controlled /
beta_naive carries that inflated beta_naive in its DENOMINATOR.  Regression to
the mean therefore deflates SF for a purely statistical reason, with no control
doing any work.  Phase 3 argued this was not driving the result (N3 returns
1.000 and N2 0.943 on the identical selected set) but never quantified it.

Two quantifications, both prespecified here:

  (A) CROSS-FITTING on the same real fits.  The 100 spatial blocks of each fit
      are split at random into halves A and B.  The selection rule is applied on
      half A; SF is then estimated on the held-out half B, where the denominator
      beta_naive_B is statistically independent of the selection.  The
      like-for-like in-sample comparator is SF estimated on half A itself, so
      the half-sample size is held constant and the ONLY difference between the
      two numbers is whether the estimate saw the selection.
          winner's curse = median SF(A | selected on A) - median SF(B | selected on A)
      Repeated over R random splits with pinned seeds.

  (B) SYNTHETIC REPLICATION of the selection rule on the Phase 1 sweep, where
      the truth is known: apply the identical rule to runs with a planted effect
      (beta_true = 1) and to runs with none (beta_true = 0), and measure how far
      selection moves the surviving fraction when the answer is known.

Everything reuses the Phase 3 estimator core: a half-sample is simply a 0/1
multiplicity vector over the blocks of the existing `BlockProfiler`, so no fit
is recomputed from the cells.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import run_phase3_nulls as RN
import phase5_common as C5

RES = C5.RES5
N_SPLIT = 20            # random block splits per fit
SEL_FRAC = 0.5          # fraction of blocks used for SELECTION (half A)
N_BOOT = 200            # block-bootstrap replicates for the selection CI
MIN_N = RN.MIN_RECEIVERS


def _sf(num, den):
    return num / den if (den is not None and np.isfinite(den)
                         and abs(den) > 1e-12) else np.nan


def _one_half(prof1, prof2, m, pp):
    """Naive fit + the controlled betas on the sub-sample defined by `m`."""
    base = prof1.fit1(m, pp["base"])
    t = base["t"]
    b0 = base["beta"]
    out = dict(lam=base["lam"], t=t, beta_naive=b0, n=base["n"],
               railed=int(t == 0 or t == prof1.lam_s.size - 1))
    out["beta_n5"] = prof2.beta_at(m, pp["n5"], t)[0]
    out["beta_n6"] = prof1.beta_at(m, pp["n6"], t)[0]
    out["beta_n6n5"] = prof1.beta_at(m, pp["n6n5"], t)[0]
    out["beta_n2n5n6"] = prof1.beta2_at(m, pp["n6n5"], t)[0]
    for k in ("n5", "n6", "n6n5", "n2n5n6"):
        out[f"sf_{k}"] = _sf(out[f"beta_{k}"], b0)
    return out


def _cell_job(sf: RN.SectionFit, ct, j, seed, sel_frac=SEL_FRAC):
    idx = sf.receivers(ct)
    n = int(idx.sum())
    if n < MIN_N:
        return None
    ii = np.flatnonzero(idx)
    y = sf.Y[ii, j].astype(float)
    d_s, d_d = sf.d_obs[ii], sf.d_dec[ii]
    bid = sf.bid[ii]
    ub, bid = np.unique(bid, return_inverse=True)
    nb = ub.size
    X1, X2, X3, pp = RN._designs(sf, idx, j)
    prof1 = E.BlockProfiler(d_s, d_d, y, X1, bid, nb, sf.lam, sf.lam)
    prof2 = E.BlockProfiler(d_s, None, y, X2, bid, nb, sf.lam, sf.lam)
    rng = np.random.default_rng(seed)
    one = np.ones(nb)

    full = _one_half(prof1, prof2, one, pp)
    # full-sample block bootstrap -> the Phase 3 selection rule
    bt = np.full(N_BOOT, np.nan)
    for r in range(N_BOOT):
        m = rng.multinomial(nb, np.full(nb, 1.0 / nb)).astype(float)
        bt[r] = prof1.beta_at(m, pp["base"], full["t"])[0]
    v = bt[np.isfinite(bt)]
    lo = float(np.quantile(v, .025)) if v.size > 20 else np.nan
    hi = float(np.quantile(v, .975)) if v.size > 20 else np.nan
    base = dict(section=sf.sec.name, arm=sf.sec.meta["condition"],
                celltype=ct, module=P.MODULES[j], n=n, nb=nb,
                lam_full=full["lam"], beta_full=full["beta_naive"],
                beta_full_lo=lo, beta_full_hi=hi,
                selected_full=int(full["beta_naive"] > 0 and lo > 0),
                sf_n5_full=full["sf_n5"], sf_n6_full=full["sf_n6"],
                sf_n6n5_full=full["sf_n6n5"],
                sf_n2n5n6_full=full["sf_n2n5n6"])

    rows = []
    for s in range(N_SPLIT):
        perm = rng.permutation(nb)
        cut = max(2, min(nb - 2, int(round(sel_frac * nb))))
        A, B = perm[:cut], perm[cut:]
        mA = np.zeros(nb); mA[A] = 1.0
        mB = np.zeros(nb); mB[B] = 1.0
        rA = _one_half(prof1, prof2, mA, pp)
        rB = _one_half(prof1, prof2, mB, pp)
        # selection CI computed WITHIN half A only
        bt = np.full(N_BOOT, np.nan)
        for r in range(N_BOOT):
            cnt = rng.multinomial(A.size, np.full(A.size, 1.0 / A.size))
            m = np.zeros(nb); m[A] = cnt
            bt[r] = prof1.beta_at(m, pp["base"], rA["t"])[0]
        v = bt[np.isfinite(bt)]
        loA = float(np.quantile(v, .025)) if v.size > 20 else np.nan
        selA = int(rA["beta_naive"] > 0 and np.isfinite(loA) and loA > 0)
        # SF on B, evaluated (i) at B's own lambda and (ii) at A's lambda
        sfB_lamA = {k: _sf(prof1.beta_at(mB, pp[k], rA["t"])[0],
                           prof1.beta_at(mB, pp["base"], rA["t"])[0])
                    for k in ("n6", "n6n5")}
        sfB_lamA["n5"] = _sf(prof2.beta_at(mB, pp["n5"], rA["t"])[0],
                             prof1.beta_at(mB, pp["base"], rA["t"])[0])
        rows.append(dict(base, split=s, selected_A=selA,
                         beta_A=rA["beta_naive"], beta_B=rB["beta_naive"],
                         lam_A=rA["lam"], lam_B=rB["lam"],
                         sf_n5_A=rA["sf_n5"], sf_n5_B=rB["sf_n5"],
                         sf_n6_A=rA["sf_n6"], sf_n6_B=rB["sf_n6"],
                         sf_n6n5_A=rA["sf_n6n5"], sf_n6n5_B=rB["sf_n6n5"],
                         sf_n2n5n6_A=rA["sf_n2n5n6"],
                         sf_n2n5n6_B=rB["sf_n2n5n6"],
                         sf_n5_B_lamA=sfB_lamA["n5"],
                         sf_n6n5_B_lamA=sfB_lamA["n6n5"]))
    return rows


def _section_job(sample, call, seed, sel_frac=SEL_FRAC):
    t0 = time.time()
    sf = RN.SectionFit(sample, call, seed)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    for ct in types:
        for j in range(len(P.MODULES)):
            r = _cell_job(sf, ct, j, seed + 7 * j + 101 * hash(ct) % 9973,
                          sel_frac)
            if r:
                rows += r
    print(f"[T2] {sample} {len(rows)} rows {time.time()-t0:.0f}s", flush=True)
    return rows


# ---------------------------------------------------------------------------
# (B) synthetic replication of the selection rule
# ---------------------------------------------------------------------------

def synthetic_selection():
    d = pd.read_csv("/workspace/results/sweep_all.csv")
    d = d[np.isfinite(d.beta_naive) & np.isfinite(d.beta_nuis)].copy()
    d["sf_n5"] = d.beta_nuis / d.beta_naive
    d["sf_n2"] = d.beta_decoyS / d.beta_naive
    d["sf_n2n5"] = d.beta_decoyS_nuis / d.beta_naive
    d["sel"] = (d.beta_naive > 0) & (d.ci_beta_naive_blk_lo > 0)
    rows = []
    for (bt, sweep), g in d.groupby(["beta_true_cfg", "sweep"]):
        for lab, sub in (("all", g), ("selected", g[g.sel])):
            if len(sub) < 5:
                continue
            rows.append(dict(beta_true=bt, sweep=sweep, subset=lab, n=len(sub),
                             sel_rate=float(g.sel.mean()),
                             sf_n5=float(np.nanmedian(sub.sf_n5)),
                             sf_n2=float(np.nanmedian(sub.sf_n2)),
                             sf_n2n5=float(np.nanmedian(sub.sf_n2n5)),
                             beta_naive=float(np.nanmedian(sub.beta_naive)),
                             beta_nuis=float(np.nanmedian(sub.beta_nuis))))
    for bt, g in d.groupby("beta_true_cfg"):
        for lab, sub in (("all", g), ("selected", g[g.sel])):
            rows.append(dict(beta_true=bt, sweep="ALL", subset=lab, n=len(sub),
                             sel_rate=float(g.sel.mean()),
                             sf_n5=float(np.nanmedian(sub.sf_n5)),
                             sf_n2=float(np.nanmedian(sub.sf_n2)),
                             sf_n2n5=float(np.nanmedian(sub.sf_n2n5)),
                             beta_naive=float(np.nanmedian(sub.beta_naive)),
                             beta_nuis=float(np.nanmedian(sub.beta_nuis))))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="both")
    ap.add_argument("--sections", default="inband")
    ap.add_argument("--call", default=RN.PRIMARY_CALL)
    ap.add_argument("--n-jobs", type=int, default=6)
    ap.add_argument("--frac", type=float, default=0.5)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    secs = {"inband": P.IN_BAND, "all": P.ALL_SECTIONS}.get(
        a.sections, a.sections.split(","))
    if a.stage in ("both", "synth"):
        s = synthetic_selection()
        s.to_csv(f"{RES}/wc_synthetic.csv", index=False)
        print(s.to_string(index=False))
    if a.stage in ("both", "crossfit"):
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_section_job)(s, a.call, P.MASTER_SEED + 900 * i, a.frac)
            for i, s in enumerate(secs))
        df = pd.DataFrame([r for rs in out for r in rs])
        df.to_csv(f"{RES}/wc_crossfit{a.tag}.csv", index=False)
        print(df.shape)
