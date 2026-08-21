"""Representative binned response curves for Figure 1c (Master Plan Sec 25).

For a handful of named regimes, regenerate tissues with `keep_curves=True` and
average the binned sender curve and the matched-decoy curve over seeds.  These
are the curves the field plots (Section 22 Step 2) with the decoy control that
Section 23 N2 prescribes drawn on top of them.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from sasp_sim import TissueConfig, simulate_tissue
from sasp_estimators import analyze, BIN_MID

MASTER_SEED = 20260820
LAM_TRUE = 30.0
N_SEEDS = 24

REGIMES = {
    "easy  (kappa=0, ell/lam=0.25, conf=0)":
        dict(clustering=0.0, autocorr_len_um=7.5, conf_strength=0.0),
    "moderate (kappa=1.5, ell/lam=1, conf=1)":
        dict(clustering=1.5, autocorr_len_um=30.0, conf_strength=1.0),
    "hard  (kappa=3, ell/lam=4, conf=1)":
        dict(clustering=3.0, autocorr_len_um=120.0, conf_strength=1.0),
    "hard+ (kappa=4.5, ell/lam=4, conf=2)":
        dict(clustering=4.5, autocorr_len_um=120.0, conf_strength=2.0),
}


def one(name, kw, rep):
    cfg = TissueConfig(lambda_true_um=LAM_TRUE, prevalence=0.05, **kw)
    sim = simulate_tissue(cfg, [MASTER_SEED, 99, hash(name) % 1000, rep])
    r = analyze(sim, [MASTER_SEED, 99, hash(name) % 1000, rep, 1],
                n_boot=50, keep_curves=True)
    return name, r["_curve_sender"], r["_curve_decoy"], r["_curve_counts"], \
        r["_curve_d"], r["lam_naive"], r["lam_decoyS"], r["beta_naive"], \
        r["beta_decoyS"]


if __name__ == "__main__":
    jobs = [(n, kw, s) for n, kw in REGIMES.items() for s in range(N_SEEDS)]
    out = Parallel(n_jobs=46, prefer="processes")(
        delayed(one)(*j) for j in jobs)

    rows = []
    for name in REGIMES:
        sel = [o for o in out if o[0] == name]
        cs = np.nanmean([o[1] for o in sel], axis=0)
        cd = np.nanmean([o[2] for o in sel], axis=0)
        cn = np.nanmean([o[3] for o in sel], axis=0)
        # centre both curves on their far-field level so they are comparable
        far = (BIN_MID > 200)
        cs0 = cs - np.nanmean(cs[far])
        cd0 = cd - np.nanmean(cd[far])
        for i, b in enumerate(BIN_MID):
            rows.append(dict(regime=name, bin_um=b, n_cells=cn[i],
                             sender_curve=cs0[i], decoy_curve=cd0[i],
                             diff_curve=cs0[i] - cd0[i],
                             truth=np.exp(-b / LAM_TRUE)))
    df = pd.DataFrame(rows)
    os.makedirs("/workspace/results", exist_ok=True)
    df.to_csv("/workspace/results/figure1c_curves.csv", index=False)
    print(df.groupby("regime").head(2))
    print("wrote /workspace/results/figure1c_curves.csv")
