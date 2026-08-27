#!/usr/bin/env python3
"""Phase 10 — the primary design WITH and WITHOUT the anatomical covariate (§15).

A6's red/white-pulp axis occupies the frozen N5 block's `zonation` slot and its
`zonation_sq` term, and its independent validation is weak: the follicle-distance check
V1 has the WRONG SIGN in 2 of 7 sections and never exceeds +0.17 raw
(`results/phase9_h1/a6_summary.csv`; `reports/CS_PHASE9_H1_AUDIT.md` §7).  §15 therefore
requires every anatomy-conditioned H1 quantity to be reported both ways.

Two of the three comparisons need no new fit -- `fit_cell` already emits them, because
`run_phase3_nulls._designs` nests base -> +N6 -> +N6+zonation -> +N6+N5:

    sf_n6      N6 only, NO anatomy, no other N5 column
    sf_n6zon   N6 + the anatomical term (+ its square)
    sf_n6n5    N6 + the FULL N5 block, which contains the anatomical term
    sf_zon     the anatomical term alone

The one that does NOT exist in the frozen output is the primary design with the anatomical
term REMOVED and everything else kept, because the nesting puts `zonation` before the rest.
This script adds exactly that one design -- X = [1, N6, N5 \\ {zonation, zonation_sq}] --
using the same `BlockProfiler`, the same cells, the same lambda index t0 = the naive
lambda-hat, and the same SF definition (beta at t0 / beta_naive at t0).  No estimator is
reimplemented and no threshold is touched.  Point estimates only: the frozen 400-replicate
bootstrap is not re-run, so this file carries no CI and says so.

Usage: python3 code/h1_a6_sensitivity.py [--calls tierAmg_p95,tierA_p95] [--n-jobs 7]
Writes results/phase10_h1/a6_sensitivity_fits.csv
"""
import argparse, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_phase10                                   # noqa: F401  ARM BINDING first
from joblib import Parallel, delayed
import sasp_phase3 as P
import h1_common as H

RES = h1_phase10.RESULTS10


def _job(sample, call, seed):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                # noqa: F401
    import numpy as _np
    import sasp_estimators as E
    import sasp_phase3 as _P
    import run_phase3_nulls as RN
    sf = RN.SectionFit(sample, call, seed)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in _P.EXCLUDE_TYPES]
    for t in types:
        idx = sf.receivers(t)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        for j in range(len(_P.MODULES)):
            y = sf.Y[idx, j].astype(float)
            d_s = sf.d_obs[idx]
            d_d = sf.d_dec[idx]
            bid = sf.bid[idx]
            ub, bid = _np.unique(bid, return_inverse=True)
            nb = ub.size
            X1, X2, X3, pp = RN._designs(sf, idx, j)
            one = _np.ones(nb)
            prof1 = E.BlockProfiler(d_s, d_d, y, X1, bid, nb, sf.lam, sf.lam)
            base = prof1.fit1(one, pp["base"])
            t0 = base["t"]
            b0 = base["beta"]
            # X1 = [1, N6, zon, zon^2, rest];  drop columns 2 and 3
            keep = [0, 1] + list(range(4, X1.shape[1]))
            Xna = X1[:, keep]
            profA = E.BlockProfiler(d_s, d_d, y, Xna, bid, nb, sf.lam, sf.lam)
            b_noanat = profA.beta_at(one, Xna.shape[1], t0)[0]
            b_n2_noanat = profA.beta2_at(one, Xna.shape[1], t0)[0]
            rows.append(dict(section=sample, call=call, celltype=t,
                             module=_P.MODULES[j], n=int(idx.sum()),
                             lam_naive=base["lam"], beta_naive=b0,
                             sd_y=float(y.std()),
                             beta_n6n5_noanat=float(b_noanat),
                             sf_n6n5_noanat=float(b_noanat / b0) if b0 else _np.nan,
                             beta_n2n5n6_noanat=float(b_n2_noanat),
                             sf_n2n5n6_noanat=float(b_n2_noanat / b0) if b0 else _np.nan))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default="tierAmg_p95,tierA_p95")
    ap.add_argument("--n-jobs", type=int, default=7)
    a = ap.parse_args()
    calls = a.calls.split(",")
    jobs = [(s, c, P.MASTER_SEED + 1000 * i + j)
            for i, s in enumerate(H.ALL_SECTIONS) for j, c in enumerate(calls)]
    t0 = time.time()
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(_job)(s, c, sd) for s, c, sd in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(RES + "/a6_sensitivity_fits.csv", index=False)
    print(df.shape, "->", RES + "/a6_sensitivity_fits.csv", "%.1f min" % ((time.time()-t0)/60))
