#!/usr/bin/env python3
"""
Where does the naive gradient come from?  A three-step decomposition that
Phase 2 could not do, because Bio's receiver cell-type labels had not landed.

For every (section, module):
  (a) UNSTRATIFIED, intercept only        -- exactly the Phase 2 / field analysis
  (b) UNSTRATIFIED + receiver cell-type dummies  -- Section 6.1's mu_{c_i}
  (c) UNSTRATIFIED + cell type + the rest of N5

beta_b / beta_a is the fraction of the naive amplitude that is NOT receiver
cell-type composition.  It is reported next to the binned-curve monotonicity
(Spearman rho over bins), which is the statistic Phase 2 used for Figure 2a.
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import run_phase3_nulls as RN


def job(sample, call):
    sf = RN.SectionFit(sample, call, P.MASTER_SEED)
    sec = sf.sec
    rows = []
    idx = sf.receivers()                       # all receiver types pooled
    ii = np.flatnonzero(idx)
    ct = sec.celltype[ii]
    types = np.array(sorted(set(ct)))
    D = np.column_stack([(ct == t).astype(float) for t in types[1:]])
    Z5 = sf.blocks["N5"][ii]
    bid = sf.bid[ii]
    ub, bid = np.unique(bid, return_inverse=True)
    nb = ub.size
    lam = sf.lam
    for j, mod in enumerate(P.MODULES):
        y = sf.Y[ii, j].astype(float)
        X = np.column_stack([np.ones(ii.size), D, Z5])
        prof = E.BlockProfiler(sf.d_obs[ii], sf.d_dec[ii], y, X, bid, nb,
                               lam, lam)
        one = np.ones(nb)
        a = prof.fit1(one, 1)
        t0 = a["t"]
        b = prof.beta_at(one, 1 + D.shape[1], t0)[0]
        c = prof.beta_at(one, X.shape[1], t0)[0]
        mu, cnt = RN._binned(sf.d_obs[ii], y)
        ok = np.isfinite(mu) & (cnt > 30)
        x = 0.5 * (RN.BINS[:-1] + RN.BINS[1:])
        rho = spearmanr(x[ok], mu[ok]).statistic if ok.sum() > 4 else np.nan
        rows.append(dict(section=sample, arm=sec.meta["condition"],
                         week=sec.meta["week"], call=call, module=mod,
                         n=int(ii.size), lam_naive=a["lam"],
                         beta_unstrat=a["beta"], beta_plus_celltype=b,
                         beta_plus_n5=c, sd_y=float(y.std()),
                         sf_celltype=b / a["beta"] if a["beta"] else np.nan,
                         sf_n5_unstrat=c / a["beta"] if a["beta"] else np.nan,
                         spearman_bins=rho, n_bins=int(ok.sum())))
    print(f"[strat] {sample} {call}", flush=True)
    return rows


if __name__ == "__main__":
    secs = P.ALL_SECTIONS
    out = Parallel(n_jobs=6, prefer="processes")(
        delayed(job)(s, RN.PRIMARY_CALL) for s in secs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{P.RESULTS}/stratification.csv", index=False)
    print(df.groupby("arm")[["spearman_bins", "sf_celltype",
                             "sf_n5_unstrat"]].median().round(3))
