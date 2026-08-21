"""
Phase 3 estimator core: designs, matched decoys, and the null battery for one
(section, sender call, receiver cell type, module) cell.

Everything is built on `sasp_estimators.BlockProfiler`, the sufficient-statistic
block profiler from Phase 1: per-block cross products are computed ONCE and every
spatial-block-bootstrap replicate is a weighted sum of small matrices.  Because
`BlockProfiler._acc(m, p)` uses `X[:, :p]`, a single profiler gives every NESTED
design (base -> +N6 -> +N6+N5) for free, and `fit2_shared` gives the shared-lambda
matched-decoy contrast that Phase 1 showed is the usable form of N2.
"""
from __future__ import annotations

import numpy as np
from typing import Dict, List, Optional

import sasp_estimators as E
import sasp_phase3 as P


# ---------------------------------------------------------------------------
# matched decoys (N2) -- same greedy 1-1 propensity algorithm as
# sasp_estimators.match_decoys, generalised off the synthetic `sim` dict.
# ---------------------------------------------------------------------------

def greedy_ps_match(ps: np.ndarray, sender: np.ndarray, strata: np.ndarray,
                    rng: np.random.Generator, caliper_sd: float = 0.25):
    """1-1 nearest-neighbour matching on a scalar score, without replacement,
    inside each stratum, with a caliper.  O(n log n): matching happens in 1-D
    by sort + searchsorted, never an (n, n) matrix."""
    n = ps.size
    send_ids = np.flatnonzero(sender)
    decoy = np.full(send_ids.size, -1, np.int64)
    used = np.zeros(n, bool)
    for k in np.unique(strata):
        s_k = np.flatnonzero(sender & (strata == k))
        c_k = np.flatnonzero((~sender) & (strata == k))
        if s_k.size == 0 or c_k.size == 0:
            continue
        order = np.argsort(ps[c_k], kind="stable")
        c_sorted = c_k[order]
        v_sorted = ps[c_sorted]
        cal = caliper_sd * float(np.std(ps[strata == k]) + 1e-12)
        for si in rng.permutation(s_k):
            pos = int(np.searchsorted(v_sorted, ps[si]))
            best, best_d = -1, np.inf
            lo, hi = pos - 1, pos
            while True:
                progressed = False
                if hi < v_sorted.size:
                    dh = abs(v_sorted[hi] - ps[si])
                    if dh < best_d:
                        if not used[c_sorted[hi]]:
                            best, best_d = c_sorted[hi], dh
                        hi += 1
                        progressed = True
                    else:
                        hi = v_sorted.size
                if lo >= 0:
                    dl = abs(v_sorted[lo] - ps[si])
                    if dl < best_d:
                        if not used[c_sorted[lo]]:
                            best, best_d = c_sorted[lo], dl
                        lo -= 1
                        progressed = True
                    else:
                        lo = -1
                if not progressed:
                    break
            if best >= 0 and best_d <= cal:
                used[best] = True
                decoy[np.searchsorted(send_ids, si)] = best
    ok = decoy >= 0
    return send_ids[ok], decoy[ok], float(ok.mean())


def match_decoys_section(sec: P.Sec, sender: np.ndarray, Zmatch: np.ndarray,
                         seed: int, caliper_sd: float = 0.25) -> Dict:
    """Matched decoys WITHIN this section and WITHIN cell type (Bio Phase 2:
    architecture differs by arm, so no pooling across sections or arms).

    Matching covariates: local density at 50 um, k-NN cell-type composition,
    log transcript counts, and -- per Section 11 -- the zonation score.
    """
    from sklearn.linear_model import LogisticRegression
    rng = np.random.default_rng(seed)
    ct = sec.celltype
    eligible = ~np.isin(ct, P.EXCLUDE_TYPES)
    types, _ = np.unique(ct, return_counts=True)
    Xp = np.column_stack([Zmatch] + [(ct == t).astype(float) for t in types[1:]])
    lr = LogisticRegression(max_iter=500, C=1.0)
    lr.fit(Xp[eligible], sender[eligible].astype(int))
    ps = np.full(sec.n, np.nan)
    ps[eligible] = E._logit(lr.predict_proba(Xp[eligible])[:, 1])
    ps = np.nan_to_num(ps, nan=-99.0)

    strata = np.searchsorted(types, ct)
    strata = np.where(eligible, strata, -1)
    s_idx, d_idx, rate = greedy_ps_match(ps, sender & eligible, strata, rng,
                                         caliper_sd)
    ms = np.zeros(sec.n, bool); ms[s_idx] = True
    md = np.zeros(sec.n, bool); md[d_idx] = True
    return dict(sender_matched=s_idx, decoy_idx=d_idx, match_rate=rate,
                smd_before=E._smd(Zmatch, sender & eligible, (~sender) & eligible),
                smd_after=E._smd(Zmatch, ms, md),
                max_smd_after=float(np.nanmax(np.abs(E._smd(Zmatch, ms, md)))),
                max_smd_before=float(np.nanmax(np.abs(
                    E._smd(Zmatch, sender & eligible, (~sender) & eligible)))))


# ---------------------------------------------------------------------------
# covariate blocks
# ---------------------------------------------------------------------------

def knn_composition(sec: P.Sec, types: np.ndarray) -> np.ndarray:
    nb = sec.celltype[sec.knn_idx]
    return np.stack([(nb == t).mean(axis=1) for t in types], axis=1)


def build_blocks(sec: P.Sec, sender: np.ndarray,
                 types: np.ndarray = None) -> Dict[str, np.ndarray]:
    """Tier D / Section 23 N5 covariates, plus the pieces N2 matches on.

    `types` fixes the k-NN-composition columns; pass P.CANON_TYPES when
    sections are pooled so every section contributes the same design columns.
    """
    if types is None:
        types = np.array(sorted(set(sec.celltype[~np.isin(sec.celltype,
                                                          P.EXCLUDE_TYPES)])))
    else:
        types = np.asarray(types)
    comp = knn_composition(sec, types)
    tc = np.log1p(sec.transcript_counts.astype(float))
    gd = np.log1p(sec.genes_detected.astype(float))
    ar = np.log1p(sec.cell_area.astype(float))
    na = np.log1p(sec.nucleus_area.astype(float))
    d25 = np.log1p(sec.z["density_25um"].astype(float))
    d50 = np.log1p(sec.z["density_50um"].astype(float))
    d100 = np.log1p(sec.z["density_100um"].astype(float))
    nn1 = sec.nn1_um.astype(float)
    zon = np.nan_to_num(sec.zonation_score.astype(float))
    db = np.log1p(np.nan_to_num(sec.dist_to_boundary_um.astype(float)))
    seg = sec.z["seg_code"].astype(int)
    segd = np.column_stack([(seg == k).astype(float) for k in (1, 2)])

    n5_cols = ["log_counts", "log_genes", "log_area", "log_nucarea",
               "log_dens25", "log_dens50", "log_dens100", "nn1_um",
               "zonation", "zonation_sq", "log_dist_boundary"]
    N5 = np.column_stack([tc, gd, ar, na, d25, d50, d100, nn1,
                          zon, zon ** 2, db])
    N5 = P.standardize(N5)
    N5 = np.column_stack([N5, P.standardize(comp), segd])
    n5_cols += [f"knn_frac_{t}" for t in types] + \
               [f"seg_{i}" for i in range(segd.shape[1])]

    Zmatch = np.column_stack([P.standardize(np.column_stack([d50, tc, zon])),
                              P.standardize(comp)])
    zmatch_cols = ["log_dens50", "log_counts", "zonation"] + \
                  [f"knn_frac_{t}" for t in types]
    return dict(N5=N5, n5_cols=n5_cols, Zmatch=Zmatch,
                zmatch_cols=zmatch_cols, types=types, zon=zon)


def neighbour_baseline(sec: P.Sec, Y: np.ndarray, sender: np.ndarray) -> np.ndarray:
    """N6: spatially smoothed expected response from the k=20 neighbourhood,
    EXCLUDING senders and excluding self.  Shape (n, n_modules)."""
    idx = sec.knn_idx
    ok = (~sender)[idx].astype(float)                    # (n, k)
    cnt = ok.sum(1)
    out = np.empty((sec.n, Y.shape[1]))
    for j in range(Y.shape[1]):
        num = (Y[idx, j] * ok).sum(1)
        out[:, j] = np.where(cnt > 0, num / np.maximum(cnt, 1), Y[:, j].mean())
    return out
