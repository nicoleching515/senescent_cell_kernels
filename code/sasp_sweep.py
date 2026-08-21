"""
Sweep runner for the synthetic identifiability study (Master Plan Sec 22 Step 1.4).

Sweeps
  main   sender clustering  x  baseline autocorrelation length / lambda_true
         (the Figure 1a/1b grid, Section 25)
  prev   sender prevalence, spanning both failure ends named in Section 8 Test 3
  conf   confounder strength
  clean  no latent confounder, random senders, short autocorrelation -- the
         reference regime for "how much true signal does the decoy control
         remove" (Section 29, objection 6)
  null   beta_true = 0, to check no estimator manufactures an effect
  nsize  number of cells N at fixed density

Engineering (Section 18)
  * BLAS pinned to one thread per worker; parallelism is across runs
    (joblib, n_jobs=-2), which is the efficient level for this workload.
  * Every (sweep, config) block is checkpointed to /workspace/results/sweep/
    as soon as it finishes, and existing blocks are skipped on restart, so a
    crash or preemption costs one block (Section 18.4).
  * ALL SEEDS PINNED: every tissue and every analysis is seeded from
    (MASTER_SEED, sweep_hash, cfg_id, rep), so the entire figure is
    reproducible from this file alone (Section 24.8, Section 26 Day 11).
"""

from __future__ import annotations

import os
# must precede numpy import in the worker processes
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import time
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from sasp_sim import TissueConfig, simulate_tissue, ripley_ratio
from sasp_estimators import analyze

MASTER_SEED = 20260820
RESULTS = "/workspace/results/sweep"
LAM_TRUE = 30.0
N_BOOT = 400

SWEEP_ID = {"main": 1, "prev": 2, "conf": 3, "clean": 4, "null": 5, "nsize": 6}

CLUSTERING_LEVELS = [0.0, 1.5, 3.0, 4.5]
ELL_RATIOS = [0.25, 0.5, 1.0, 2.0, 4.0]


def build_configs() -> Dict[str, List[Tuple[int, dict, int]]]:
    """(cfg_id, config kwargs, n_reps) per sweep.  cfg_id is stable, so seeds
    and checkpoint filenames never shift when a sweep is extended."""
    S: Dict[str, List[Tuple[int, dict, int]]] = {}

    # --- main grid: clustering x autocorrelation length -------------------
    main = []
    cid = 0
    for kap in CLUSTERING_LEVELS:
        for rat in ELL_RATIOS:
            main.append((cid, dict(clustering=kap,
                                   autocorr_len_um=rat * LAM_TRUE,
                                   conf_strength=1.0, prevalence=0.05,
                                   ell_ratio=rat), 30))
            cid += 1
    S["main"] = main

    # --- sender prevalence (Section 8 Test 3: <1% and >20% break it) ------
    S["prev"] = [(i, dict(prevalence=p, clustering=1.5, autocorr_len_um=30.0,
                          conf_strength=1.0), 30)
                 for i, p in enumerate([0.005, 0.01, 0.02, 0.05, 0.10,
                                        0.20, 0.30])]

    # --- confounder strength ----------------------------------------------
    S["conf"] = [(i, dict(conf_strength=c, clustering=3.0,
                          autocorr_len_um=60.0, prevalence=0.05), 30)
                 for i, c in enumerate([0.0, 0.25, 0.5, 1.0, 1.5, 2.0])]

    # --- clean reference regimes ------------------------------------------
    S["clean"] = [
        (0, dict(clustering=0.0, conf_strength=0.0, autocorr_len_um=7.5,
                 prevalence=0.05, regime="easy"), 40),
        (1, dict(clustering=0.0, conf_strength=0.0, autocorr_len_um=7.5,
                 prevalence=0.05, base_grf_amp=0.0, gamma_density=0.0,
                 gamma_counts=0.0, density_mod=0.0, regime="pure"), 40),
    ]

    # --- null: no planted kernel ------------------------------------------
    S["null"] = [
        (0, dict(beta_true=0.0, clustering=0.0, conf_strength=0.0,
                 autocorr_len_um=7.5), 30),
        (1, dict(beta_true=0.0, clustering=1.5, conf_strength=1.0,
                 autocorr_len_um=30.0), 30),
        (2, dict(beta_true=0.0, clustering=3.0, conf_strength=1.0,
                 autocorr_len_um=120.0), 30),
        (3, dict(beta_true=0.0, clustering=3.0, conf_strength=2.0,
                 autocorr_len_um=120.0), 30),
    ]

    # --- N (window size at fixed density) ---------------------------------
    S["nsize"] = [(i, dict(window_um=w, clustering=1.5, autocorr_len_um=30.0,
                           conf_strength=1.0), 20)
                  for i, w in enumerate([1000.0, 1414.0, 2000.0, 2828.0])]
    return S


# --------------------------------------------------------------------------


def run_one(sweep: str, cfg_id: int, kw: dict, rep: int) -> dict:
    meta = {k: kw[k] for k in ("ell_ratio", "regime") if k in kw}
    cfg_kw = {k: v for k, v in kw.items() if k not in ("ell_ratio", "regime")}
    cfg = TissueConfig(lambda_true_um=LAM_TRUE, **cfg_kw)

    sid = SWEEP_ID[sweep]
    sim = simulate_tissue(cfg, [MASTER_SEED, sid, cfg_id, rep])
    rec = analyze(sim, [MASTER_SEED, sid, cfg_id, rep, 1], n_boot=N_BOOT)
    rec = {k: v for k, v in rec.items() if not k.startswith("_")}

    rec["ripley50"] = ripley_ratio(sim["coords"], sim["sender_mask"], 50.0)
    rec.update(sweep=sweep, cfg_id=cfg_id, rep=rep,
               clustering=cfg.clustering, autocorr_len_um=cfg.autocorr_len_um,
               conf_strength=cfg.conf_strength, prevalence=cfg.prevalence,
               window_um=cfg.window_um, beta_true_cfg=cfg.beta_true,
               ell_over_lambda=cfg.autocorr_len_um / LAM_TRUE)
    rec.update(meta)
    return rec


def run_sweep(sweep: str, configs, n_jobs: int = -2, force: bool = False):
    os.makedirs(RESULTS, exist_ok=True)
    frames = []
    for cfg_id, kw, n_reps in configs:
        path = os.path.join(RESULTS, f"{sweep}_{cfg_id:03d}.csv")
        if os.path.exists(path) and not force:
            frames.append(pd.read_csv(path))
            print(f"  [skip] {sweep}/{cfg_id:03d} (checkpoint exists)", flush=True)
            continue
        t0 = time.time()
        recs = Parallel(n_jobs=n_jobs, prefer="processes", batch_size=2)(
            delayed(run_one)(sweep, cfg_id, kw, rep) for rep in range(n_reps))
        df = pd.DataFrame(recs)
        df.to_csv(path, index=False)          # checkpoint immediately
        frames.append(df)
        print(f"  [done] {sweep}/{cfg_id:03d} n={n_reps} "
              f"{time.time() - t0:5.1f}s  {json.dumps(kw)}", flush=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweeps", default="main,prev,conf,clean,null,nsize")
    ap.add_argument("--n-jobs", type=int, default=-2)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    S = build_configs()
    out = []
    for name in args.sweeps.split(","):
        name = name.strip()
        if not name:
            continue
        print(f"=== sweep: {name} ({len(S[name])} configs) ===", flush=True)
        t0 = time.time()
        out.append(run_sweep(name, S[name], args.n_jobs, args.force))
        print(f"=== {name} finished in {time.time() - t0:.1f}s ===", flush=True)

    alldf = pd.concat(out, ignore_index=True)
    alldf.to_csv("/workspace/results/sweep_all.csv", index=False)
    print(f"\nwrote /workspace/results/sweep_all.csv  "
          f"{alldf.shape[0]} runs x {alldf.shape[1]} columns")


if __name__ == "__main__":
    main()
