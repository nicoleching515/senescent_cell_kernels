#!/usr/bin/env python3
"""Why does 'the spline wins 34/35' (CS Phase 2 §8) not reproduce in Phase 5?

Phase 5 fits the five families inside the Phase 3 window (100 um, lambda grid
[7, 50] um, spline dmax = 100) and WITHIN receiver cell type.  Phase 2 fitted
them unstratified over a 300 um window with lambda on [3, 400].  This isolates
which of the two changes is responsible, on the same six admissible sections.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import sasp_kernels as K
import run_phase5_kernels as T4

CFG = {
    "phase5_win100_grid7_50": dict(win=100.0, lam=RN.lam_grid(), dmax=100.0),
    "phase2_win300_grid3_400": dict(win=300.0,
                                    lam=np.exp(np.linspace(np.log(3.0),
                                                           np.log(400.0), 64)),
                                    dmax=300.0),
}


def job(sample, cfgname, stratified):
    cfg = CFG[cfgname]
    sf = RN.SectionFit(sample, RN.PRIMARY_CALL, P.MASTER_SEED)
    rng = np.random.default_rng(1)
    rows = []
    groups = ([t for t in sorted(set(sf.sec.celltype))
               if t not in P.EXCLUDE_TYPES] if stratified else [None])
    for ct in groups:
        m = ((~np.isin(sf.sec.celltype, P.EXCLUDE_TYPES)) & (~sf.sender)
             & np.isfinite(sf.d_obs) & (sf.d_obs <= cfg["win"]))
        if ct is not None:
            m &= sf.sec.celltype == ct
        if m.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(m)
        d = sf.d_obs[ii]
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        X = np.ones((ii.size, 1))
        _, KN = K.spline_basis(d, n_knots=T4.N_KNOTS, dmax=cfg["dmax"])
        for j, mod in enumerate(P.MODULES):
            y = sf.Y[ii, j].astype(float)
            res = {}
            for fam in T4.FAMS:
                res[fam] = T4.fit_family_basis(d, y, X, bid, ub.size, fam,
                                               cfg["lam"], knots=KN, rng=rng)
            best = min(T4.FAMS, key=lambda f: res[f]["aic"])
            for fam in T4.FAMS:
                rows.append(dict(section=sample, cfg=cfgname,
                                 stratified=int(stratified),
                                 celltype=ct or "ALL", module=mod, family=fam,
                                 n=int(ii.size), aic=res[fam]["aic"],
                                 lam=res[fam]["lam"],
                                 d_half=res[fam].get("d_half", np.nan),
                                 best_family=best,
                                 d_aic_vs_cov=res[fam]["aic"] - res[fam]["aic0"]))
    print(f"  {sample} {cfgname} strat={stratified}: {len(rows)}", flush=True)
    return rows


if __name__ == "__main__":
    jobs = [(s, c, st) for s in P.IN_BAND for c in CFG
            for st in (False, True)]
    out = Parallel(n_jobs=12, prefer="processes")(
        delayed(job)(*j) for j in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv("/workspace/results/phase5/spline_window_check.csv", index=False)
    for (c, st), g in df.groupby(["cfg", "stratified"]):
        w = g.groupby("family").apply(
            lambda h: (h.family == h.best_family).mean())
        nf = g.groupby(["section", "celltype", "module"]).ngroups
        print(f"\n{c}  stratified={st}  ({nf} fits)")
        print(w.round(3).to_string())
