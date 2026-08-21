#!/usr/bin/env python3
"""SE understatement factor for the Phase 5 estimands.

Phase 1 measured up to 7.9x understatement for iid asymptotic SEs under spatial
confounding and Phases 2-3 measured 4-5x on real tissue.  Every Phase 5 interval
is a spatial block bootstrap; this script measures what the iid interval would
have claimed instead, for the two NEW estimands (the Section 6.3 superposition
amplitude and the Section 6.2 length scale per family), so the ratio is on the
record rather than assumed.
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
import phase5_common as C5
import run_phase5_super as T1

N_BOOT = 300


def job(sample, seed):
    sf = RN.SectionFit(sample, RN.PRIMARY_CALL, seed)
    S = T1.Sup(sf)
    rng = np.random.default_rng(seed)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    for ct in types:
        idx = sf.receivers(ct)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        nb = ub.size
        for j, mod in enumerate(P.MODULES):
            y = sf.Y[ii, j].astype(float)
            X1, _, _, pp = RN._designs(sf, idx, j)
            BM = rng.multinomial(nb, np.full(nb, 1.0 / nb),
                                 size=N_BOOT).astype(float)
            for dname, X in (("ctrl", X1[:, :pp["n6n5"]]),
                             ("naive", X1[:, :1])):
              Q, _ = np.linalg.qr(X)
              yt = y - Q @ (Q.T @ y)
              for nm, Kall in (("nearest", S.Ks_near), ("superposition",
                                                        S.Ks_sup)):
                  Kb = Kall[ii]
                  r = C5.profile_basis(Kb, X, y, bid, nb, sf.lam)
                  t = r["t"]
                  k = Kb[:, t]
                  kt = k - Q @ (Q.T @ k)
                  kk = float(kt @ kt)
                  if kk <= 1e-12:
                      continue
                  beta = float(kt @ yt) / kk
                  resid = yt - beta * kt
                  s2 = float(resid @ resid) / max(y.size - X.shape[1] - 2, 1)
                  se_iid = float(np.sqrt(s2 / kk))
                  # block bootstrap of beta at the SAME grid index
                  prof = r["prof"]
                  bb = np.full(N_BOOT, np.nan)
                  for b in range(N_BOOT):
                      m = BM[b]
                      XX = np.tensordot(m, prof.XX, 1)
                      XK = np.tensordot(m, prof.XK, 1)
                      Xy = np.tensordot(m, prof.Xy, 1)
                      Gi = np.linalg.pinv(XX)
                      A = Gi @ XK
                      kkb = float(np.tensordot(m, prof.KK, 1)[t]
                                  - XK[:, t] @ A[:, t])
                      kyb = float(np.tensordot(m, prof.Ky, 1)[t]
                                  - XK[:, t] @ (Gi @ Xy))
                      if kkb > 1e-12:
                          bb[b] = kyb / kkb
                  se_blk = float(np.nanstd(bb))
                  rows.append(dict(section=sample, celltype=ct, module=mod,
                                   design=dname, basis=nm,
                                   n=int(ii.size), beta=beta,
                                   se_iid=se_iid, se_blk=se_blk,
                                   se_ratio=se_blk / se_iid if se_iid > 0
                                   else np.nan,
                                   ci_iid_lo=beta - 1.96 * se_iid,
                                   ci_iid_hi=beta + 1.96 * se_iid,
                                   ci_blk_lo=float(np.nanquantile(bb, .025)),
                                   ci_blk_hi=float(np.nanquantile(bb, .975))))
    print(f"[se] {sample} {len(rows)}", flush=True)
    return rows


if __name__ == "__main__":
    out = Parallel(n_jobs=6, prefer="processes")(
        delayed(job)(s, P.MASTER_SEED + 77 * i)
        for i, s in enumerate(P.IN_BAND))
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv("/workspace/results/phase5/se_ratio.csv", index=False)
    for (dz, b), g in df.groupby(["design", "basis"]):
        sig_i = ((g.ci_iid_lo > 0) | (g.ci_iid_hi < 0)).mean()
        sig_b = ((g.ci_blk_lo > 0) | (g.ci_blk_hi < 0)).mean()
        print(f"{dz:6s} {b:14s} n={len(g)}  SE ratio (block/iid) median "
              f"{g.se_ratio.median():.2f} "
              f"[{g.se_ratio.quantile(.25):.2f}, {g.se_ratio.quantile(.75):.2f}]"
              f"  p90 {g.se_ratio.quantile(.9):.2f}  |  CI excludes 0: "
              f"iid {sig_i:.3f}  block {sig_b:.3f}")
