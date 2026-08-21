#!/usr/bin/env python3
"""Which part of the N5 block removes the naive amplitude?

Nested sub-blocks, each added to the intercept-only design in turn and then
cumulatively, all evaluated at the naive lambda_hat so the numbers are
comparable (CS Phase 2 §10).

  tech   log transcript counts, log genes detected, log cell area, log nucleus area
  seg    segmentation_method (3 levels, spatially patterned)
  dens   local density at 25/50/100 um and the 1-NN distance
  comp   k-NN (k=20) receiver cell-type composition
  anat   zonation score (linear + quadratic) and log distance to tissue boundary
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import run_phase3_nulls as RN

BLOCKS = ["tech", "seg", "dens", "comp", "anat"]


def _cols(sf):
    n5 = sf.blocks["n5_cols"]
    Z = sf.blocks["N5"]
    ix = {
        "tech": [n5.index(c) for c in ("log_counts", "log_genes", "log_area",
                                       "log_nucarea")],
        "seg": [i for i, c in enumerate(n5) if c.startswith("seg_")],
        "dens": [n5.index(c) for c in ("log_dens25", "log_dens50",
                                       "log_dens100", "nn1_um")],
        "comp": [i for i, c in enumerate(n5) if c.startswith("knn_frac_")],
        "anat": [n5.index(c) for c in ("zonation", "zonation_sq",
                                       "log_dist_boundary")],
    }
    return Z, ix


def job(sample, call):
    sf = RN.SectionFit(sample, call, P.MASTER_SEED)
    Z, ix = _cols(sf)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    for t in types:
        idx = sf.receivers(t)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        nb = ub.size
        one = np.ones(nb)
        for j, mod in enumerate(P.MODULES):
            y = sf.Y[ii, j].astype(float)
            base = np.ones((ii.size, 1))
            prof0 = E.BlockProfiler(sf.d_obs[ii], None, y, base, bid, nb,
                                    sf.lam, sf.lam)
            b0 = prof0.fit1(one, 1)
            t0 = b0["t"]
            if not np.isfinite(b0["beta"]) or b0["beta"] == 0:
                continue
            r = dict(section=sample, arm=sf.sec.meta["condition"],
                     band=sf.sec.meta["band"], call=call, celltype=t,
                     module=mod, n=int(ii.size), lam=b0["lam"],
                     beta_naive=b0["beta"], sd_y=float(y.std()))
            # each block alone
            cum = []
            for blk in BLOCKS:
                Xb = np.column_stack([base, Z[np.ix_(ii, ix[blk])]])
                pb = E.BlockProfiler(sf.d_obs[ii], None, y, Xb, bid, nb,
                                     sf.lam, sf.lam)
                r[f"sf_{blk}"] = pb.beta_at(one, Xb.shape[1], t0)[0] / b0["beta"]
                cum += ix[blk]
                Xc = np.column_stack([base, Z[np.ix_(ii, cum)]])
                pc = E.BlockProfiler(sf.d_obs[ii], None, y, Xc, bid, nb,
                                     sf.lam, sf.lam)
                r[f"sf_cum_{blk}"] = pc.beta_at(one, Xc.shape[1],
                                                t0)[0] / b0["beta"]
            rows.append(r)
    print(f"[attr] {sample} {call} {len(rows)}", flush=True)
    return rows


if __name__ == "__main__":
    jobs = [(s, RN.PRIMARY_CALL) for s in P.IN_BAND]
    out = Parallel(n_jobs=6, prefer="processes")(
        delayed(job)(s, c) for s, c in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{P.RESULTS}/attribution.csv", index=False)
    rep = df[(df.beta_naive > 0)]
    print("\nSurviving fraction with EACH block alone (median over "
          f"{len(rep)} fits):")
    for b in BLOCKS:
        print(f"  {b:6s} {rep['sf_'+b].median():+.3f}")
    print("\nCumulative (added in order tech -> seg -> dens -> comp -> anat):")
    for b in BLOCKS:
        print(f"  +{b:6s} {rep['sf_cum_'+b].median():+.3f}")
