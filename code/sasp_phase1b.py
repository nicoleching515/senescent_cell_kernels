"""Phase 1b: the synthetic work still owed from Phase 1.

1. KERNEL FAMILY MISSPECIFICATION -- plant each family in turn, fit all five,
   and ask (a) does AIC recover the truth, (b) how wrong is the reported length
   scale when the wrong family is fitted?  Phase 1 planted and fitted only the
   exponential, so misspecification was completely untested.
2. SUPERPOSITION vs NEAREST-SENDER (Section 6.3) -- plant each, fit both,
   and ask whether model comparison can tell them apart.
3. NULLS N1 / N3 / N4 -- size (beta_true = 0) and power (real kernel), with the
   full null distributions.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.spatial import cKDTree

from sasp_sim import TissueConfig, simulate_tissue, planted_kernel
from sasp_estimators import block_ids, design_matrix
import sasp_kernels as K
from sasp_nulls import run_nulls

MASTER_SEED = 20260820
OUT = "/workspace/results/phase1b"
LAM_TRUE = 30.0
FAMS = ["exponential", "gaussian", "powerlaw", "step", "spline"]

REGIMES = {
    "clean": dict(clustering=0.0, conf_strength=0.0, autocorr_len_um=7.5),
    "confounded": dict(clustering=3.0, conf_strength=1.0, autocorr_len_um=60.0),
}


def true_d_half(family, lam, p=2.0):
    dd = np.linspace(0, 600, 6001)
    k = planted_kernel(dd, family, lam, p)
    below = np.flatnonzero(k <= 0.5)
    return float(dd[below[0]]) if below.size else np.nan


# --------------------------------------------------------------------------
# 1 + 2: misspecification and superposition
# --------------------------------------------------------------------------


def one_misspec(true_fam: str, superpos: bool, regime: str, rep: int,
                n_boot: int = 120) -> pd.DataFrame:
    cfg = TissueConfig(lambda_true_um=LAM_TRUE, kernel_family=true_fam,
                       superposition=superpos, prevalence=0.05,
                       **REGIMES[regime])
    sim = simulate_tissue(cfg, [MASTER_SEED, 11, hash(true_fam) % 997,
                                int(superpos), hash(regime) % 997, rep])
    idx = np.flatnonzero(~sim["sender_mask"])
    y = sim["r"][idx]
    d = sim["d_sender"][idx]
    X0 = design_matrix(sim, idx, with_nuisance=False)
    bid = block_ids(sim["coords"][idx], sim["window_um"], 8)
    rng = np.random.default_rng([MASTER_SEED, rep, 7])

    rows = []
    base = dict(true_family=true_fam, true_superposition=superpos,
                regime=regime, rep=rep, n=idx.size,
                true_d_half=true_d_half(true_fam, LAM_TRUE, cfg.kernel_p))
    for fam in FAMS:
        r = K.fit_family(d, y, X0, fam, bid, 64, n_boot=n_boot, rng=rng)
        rows.append(dict(base, fit_family=fam, fit_mode="nearest",
                         lam=r["params"].get("lam", np.nan),
                         d_half=r["d_half"], aic=r["aic"],
                         beta=r.get("beta", np.nan), rss=r["rss"],
                         k_params=r["k_params"]))

    # --- superposition FIT: sum of exponential kernels over all senders ----
    # profiled over lambda exactly like the nearest-sender fit
    send_xy = sim["coords"][sim["sender_mask"]]
    tc = cKDTree(sim["coords"][idx], boxsize=sim["window_um"])
    ts = cKDTree(send_xy, boxsize=sim["window_um"])
    lam_grid = K.LAM_GRID
    trunc = 6.0 * lam_grid.max()
    D = tc.sparse_distance_matrix(ts, min(trunc, sim["window_um"] / 2 - 1),
                                  output_type="coo_matrix")
    Ksup = np.zeros((idx.size, lam_grid.size))
    for j, lam in enumerate(lam_grid):
        v = np.exp(-D.data / lam)
        np.add.at(Ksup[:, j], D.row, v)
    prof = K.BasisBlockProfiler(Ksup, X0, y, bid, 64)
    rss, beta, _ = prof.profile(np.ones(64))
    t = int(np.nanargmin(rss))
    kp = X0.shape[1] + 2
    rows.append(dict(base, fit_family="exponential", fit_mode="superposition",
                     lam=float(lam_grid[t]), d_half=np.nan,
                     aic=idx.size * np.log(rss[t] / idx.size) + 2 * (kp + 1),
                     beta=float(beta[t]), rss=float(rss[t]), k_params=kp))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 3: nulls
# --------------------------------------------------------------------------


def one_null(regime: str, beta_true: float, rep: int, n_perm: int):
    cfg = TissueConfig(lambda_true_um=LAM_TRUE, beta_true=beta_true,
                       prevalence=0.05, **REGIMES[regime])
    sim = simulate_tissue(cfg, [MASTER_SEED, 12, hash(regime) % 997,
                                int(beta_true * 10), rep])
    r = run_nulls(sim, [MASTER_SEED, 12, rep, 3], n_perm=n_perm)
    r.update(regime=regime, beta_true_cfg=beta_true, rep=rep)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="misspec,nulls")
    ap.add_argument("--reps", type=int, default=15)
    ap.add_argument("--null-reps", type=int, default=12)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-jobs", type=int, default=46)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    parts = a.part.split(",")

    if "misspec" in parts and not os.path.exists(f"{OUT}/misspec.csv"):
        jobs = [(f, sup, reg, rep)
                for f in ["exponential", "gaussian", "powerlaw", "step"]
                for sup in (False, True)
                for reg in REGIMES
                for rep in range(a.reps)]
        print(f"misspec: {len(jobs)} runs", flush=True)
        res = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(one_misspec)(*j) for j in jobs)
        pd.concat(res, ignore_index=True).to_csv(f"{OUT}/misspec.csv",
                                                 index=False)
        print("wrote misspec.csv", flush=True)

    if "nulls" in parts and not os.path.exists(f"{OUT}/nulls.csv"):
        jobs = [(reg, bt, rep) for reg in REGIMES for bt in (0.0, 1.0)
                for rep in range(a.null_reps)]
        print(f"nulls: {len(jobs)} runs x {a.n_perm} permutations", flush=True)
        res = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(one_null)(r, b, p, a.n_perm) for r, b, p in jobs)
        pd.DataFrame(res).to_csv(f"{OUT}/nulls.csv", index=False)
        print("wrote nulls.csv", flush=True)


if __name__ == "__main__":
    main()
