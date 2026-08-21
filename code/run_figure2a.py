"""
Section 22 Step 2 -- naive estimation on real tissue, and Figure 2a.

"Reproduce the standard analysis, deliberately without controls, so you have
the thing you are going to test."  Bin cells by distance to nearest sender
(10 um bins out to 300 um), plot mean response per bin per receiver cell type
per tissue, fit every Section 6.2 kernel family, report lambda_hat_naive and
model comparison.

Phase 1 lessons carried forward:
  * report the iid asymptotic CI AND the spatial block bootstrap CI side by
    side (Phase 1 measured an SE understatement factor up to 7.9x);
  * lead with lambda / d_half, not beta (Phase 1's beta_true = 0 null showed
    the naive beta is manufactured 90-100% of the time under confounding);
  * report distance in units of the LOCAL median nearest-neighbour distance as
    well as in raw microns, because these sections differ in packing by ~43%
    and Phase 1 showed local density mechanically confounds
    distance-to-nearest-sender with no latent field required.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import sasp_real as R
import sasp_kernels as K
from prepare_samples import Cached

RESULTS = "/workspace/results/real"
MASTER_SEED = 20260820
BIN_EDGES = np.arange(0.0, 310.0, 10.0)
BIN_MID = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
N_BLOCKS_SIDE = 10
P_GRID_REAL = np.array([0.5, 1.0, 2.0, 4.0, 8.0])


def block_ids(coords, n_side=N_BLOCKS_SIDE):
    out = np.zeros(coords.shape[0], dtype=np.int64)
    for k in (0, 1):
        v = coords[:, k]
        q = np.quantile(v, np.linspace(0, 1, n_side + 1)[1:-1])
        out = out * n_side + np.searchsorted(q, v)
    return out


def binned(d, y, edges=BIN_EDGES):
    b = np.floor(d / 10.0).astype(np.int64)
    ok = (b >= 0) & (d < edges[-1])
    nb = edges.size - 1
    cnt = np.bincount(b[ok], minlength=nb).astype(float)
    tot = np.bincount(b[ok], weights=y[ok], minlength=nb)
    sq = np.bincount(b[ok], weights=y[ok] ** 2, minlength=nb)
    dm = np.bincount(b[ok], weights=d[ok], minlength=nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = tot / cnt
        sd = np.sqrt(np.maximum(sq / cnt - m ** 2, 0))
        se = sd / np.sqrt(np.maximum(cnt, 1))
        dmm = dm / cnt
    return m, cnt, se, dmm


def analyze_one(sample: str, module: str, sender_key: str,
                n_boot: int = 300, families=K.FAMILIES,
                stratify_celltype: bool = True):
    """Naive kernel estimation for one (section, module, sender call).

    Stratified by receiver cell type when a cell-type annotation is available
    (Bio's `celltypes_*.csv`); when it is not, a single 'ALL' stratum is used
    and the fact is recorded in `celltype_source` so no figure can silently
    imply a stratification that did not happen.
    """
    s = Cached(sample)
    if module not in s.modules() or sender_key not in s.senders():
        return None, None
    smask = s.sender(sender_key)
    y_all = s.module(module)
    d_all = s.dist(sender_key)
    coords_all = s.coords
    ct_all = s.celltype

    valid = (~smask) & np.isfinite(y_all) & np.isfinite(d_all)
    strata = ["ALL"]
    if stratify_celltype:
        cats = pd.unique(ct_all[valid])
        if len(cats) > 1:
            strata += [c for c in sorted(cats)
                       if (valid & (ct_all == c)).sum() >= 2000]

    rows, curves = [], []
    for stratum in strata:
        keep = valid if stratum == "ALL" else (valid & (ct_all == stratum))
        idx = np.flatnonzero(keep)
        if idx.size < 5000:
            continue
        coords = coords_all[idx]
        y = y_all[idx].astype(float)
        d = d_all[idx].astype(float)
        ct = ct_all[idx]

        # receiver cell-type intercepts are part of the model (Section 6.1)
        cols = [np.ones(idx.size)]
        if stratum == "ALL":
            for c in sorted(pd.unique(ct))[1:]:
                cols.append((ct == c).astype(float))
        X0 = np.column_stack(cols)

        # --- Tier-D style nuisance covariates, for the counts/density-adjusted
        # preview fit.  Phase 1 showed local density mechanically confounds
        # distance-to-nearest-sender; on this panel the module scores also
        # correlate 0.22-0.67 with total transcript count, and the sender score
        # itself correlates 0.29, so both need to be visible here.
        zz = []
        for c in ("transcript_counts", "cell_area", "density_25um",
                  "density_50um", "density_100um", "nn1_um"):
            v = s.col(c)[idx].astype(float)
            v = np.log1p(v) if c in ("transcript_counts", "cell_area") else v
            zz.append((v - v.mean()) / (v.std() + 1e-12))
        X_adj = np.column_stack(cols + zz)

        bid = block_ids(coords)
        nb = N_BLOCKS_SIDE ** 2
        rng = np.random.default_rng(np.random.SeedSequence(
            [MASTER_SEED, abs(hash(sample)) % 10**6,
             abs(hash(module)) % 10**6, abs(hash(stratum)) % 10**6]))

        m, cnt, se, dmm = binned(d, y)
        curves.append(pd.DataFrame(dict(
            sample=sample, module=module, sender_call=sender_key,
            celltype=stratum, bin_um=BIN_MID, bin_mean_d_um=dmm,
            mean=m, n=cnt, sem=se,
            bin_nn_units=BIN_MID / s.median_nn_um)))

        base = dict(sample=sample, module=module, sender_call=sender_key,
                    celltype=stratum, n_cells=int(idx.size),
                    n_senders=int(smask.sum()),
                    prevalence=float(smask.mean()),
                    median_nn_um=s.median_nn_um,
                    med_d_sender_um=float(np.median(d)),
                    celltype_source=s.sources.get("celltype", "?"),
                    sender_source=s.sources.get("sender", "?"),
                    module_source=s.sources.get("module", "?"),
                    condition=s.meta.get("condition", "?"),
                    week=s.meta.get("week", "?"),
                    mouse=s.meta.get("mouse", "?"))

        def pack(r, fam, units):
            return dict(base, family=fam, units=units,
                        lam=r["params"].get("lam", np.nan),
                        pow_p=r["params"].get("p", np.nan),
                        beta=r.get("beta", np.nan),
                        d_half=r["d_half"], d_05=r["d_05"],
                        amp=r.get("amp", np.nan),
                        aic=r["aic"],
                        delta_aic_vs_null=r["delta_aic_vs_null"],
                        rss=r["rss"], k_params=r["k_params"],
                        se_lam_iid=r.get("se_lam_iid", np.nan),
                        ci_lam_iid_lo=r.get("ci_lam_iid", (np.nan,) * 2)[0],
                        ci_lam_iid_hi=r.get("ci_lam_iid", (np.nan,) * 2)[1],
                        ci_lam_blk_lo=r.get("ci_lam_blk", (np.nan,) * 2)[0],
                        ci_lam_blk_hi=r.get("ci_lam_blk", (np.nan,) * 2)[1],
                        ci_dhalf_blk_lo=r.get("ci_d_half_blk", (np.nan,) * 2)[0],
                        ci_dhalf_blk_hi=r.get("ci_d_half_blk", (np.nan,) * 2)[1],
                        se_ratio=r.get("se_ratio", np.nan),
                        lam_at_bound=float(
                            np.isfinite(r["params"].get("lam", np.nan))
                            and (r["params"].get("lam", np.nan) <= K.LAM_LO * 1.02
                                 or r["params"].get("lam", np.nan) >= K.LAM_HI * 0.98)))

        for fam in families:
            r = K.fit_family(d, y, X0, fam, bid, nb, n_boot=n_boot, rng=rng,
                             p_grid=P_GRID_REAL)
            rows.append(pack(r, fam, "um"))

        # preview of null N5: same fit with counts / density / area covariates
        for fam in ("exponential", "spline"):
            r = K.fit_family(d, y, X_adj, fam, bid, nb, n_boot=n_boot,
                             rng=rng, p_grid=P_GRID_REAL)
            rows.append(pack(r, fam, "um_covadj"))

        # sensitivity: distance in units of the LOCAL median NN distance.
        # These sections differ in packing by ~43%, and Phase 1 showed local
        # density mechanically confounds distance-to-nearest-sender.
        r = K.fit_family(d / s.median_nn_um, y, X0, "exponential", bid, nb,
                         n_boot=n_boot, rng=rng,
                         lam_grid=K.LAM_GRID / s.median_nn_um)
        rows.append(pack(r, "exponential", "median_NN"))

    if not rows:
        return None, None
    return pd.DataFrame(rows), pd.concat(curves, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="")
    ap.add_argument("--modules", default="")
    ap.add_argument("--sender", default="tierA_p95")
    ap.add_argument("--n-jobs", type=int, default=24)
    ap.add_argument("--n-boot", type=int, default=300)
    args = ap.parse_args()

    os.makedirs(RESULTS, exist_ok=True)
    samples = args.samples.split(",") if args.samples else R.list_samples()
    # only sections whose cache exists; more sections keep arriving, so
    # prepare_samples.py is the gate, not this script
    from prepare_samples import CACHE
    have = {f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz")}
    missing = [s for s in samples if s not in have]
    if missing:
        print(f"NOT CACHED (run prepare_samples.py): {missing}", flush=True)
    samples = [s for s in samples if s in have]
    if args.modules:
        modules = args.modules.split(",")
    else:
        modules = sorted(Cached(samples[0]).modules())
    print(f"samples={samples}\nmodules={modules}\nsender={args.sender}",
          flush=True)

    jobs = [(s, m) for s in samples for m in modules]
    todo = [(s, m) for s, m in jobs
            if not os.path.exists(f"{RESULTS}/fit_{s}__{m}__{args.sender}.csv")]
    print(f"{len(jobs)} jobs, {len(jobs)-len(todo)} already checkpointed",
          flush=True)

    def run(sm):
        s, m = sm
        fit, curve = analyze_one(s, m, args.sender, n_boot=args.n_boot)
        if fit is None:
            return None
        fit.to_csv(f"{RESULTS}/fit_{s}__{m}__{args.sender}.csv", index=False)
        curve.to_csv(f"{RESULTS}/curve_{s}__{m}__{args.sender}.csv", index=False)
        return f"{s}/{m}"

    done = Parallel(n_jobs=args.n_jobs, prefer="processes", verbose=5)(
        delayed(run)(sm) for sm in todo)
    print("finished:", [d for d in done if d], flush=True)

    fits = pd.concat([pd.read_csv(f"{RESULTS}/{f}") for f in os.listdir(RESULTS)
                      if f.startswith("fit_")], ignore_index=True)
    curves = pd.concat([pd.read_csv(f"{RESULTS}/{f}") for f in os.listdir(RESULTS)
                        if f.startswith("curve_")], ignore_index=True)
    fits.to_csv("/workspace/results/real_fits.csv", index=False)
    curves.to_csv("/workspace/results/real_curves.csv", index=False)
    print(f"wrote real_fits.csv {fits.shape}, real_curves.csv {curves.shape}")


if __name__ == "__main__":
    main()
