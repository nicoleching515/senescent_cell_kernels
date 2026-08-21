"""
Phase 5 shared machinery — Sections 6.2, 6.3, 6.4 and 24.6.

Everything here is additive on top of Phase 3.  In particular it REUSES
`run_phase3_nulls.SectionFit` (which owns the receiver definition, the window,
the lambda grid, the N5/N6 covariate blocks and the N2 matched decoys) and
`sasp_kernels.BasisBlockProfiler` (per-block Gram matrices, so a spatial block
bootstrap replicate is a weighted sum of small matrices).  No fitter is
rewritten.

Two things are genuinely new:

  1. `superposition_basis` — the Section 6.3 aggregate-sender regressor
     S_i(lam) = sum_{j in S} exp(-||x_i - x_j|| / lam), built with a chunked
     cKDTree sparse distance matrix truncated at 6*lam_max (contribution beyond
     that is < 0.25 %).  An (n, n) matrix is never formed.

  2. `heldout_ll` — leave-one-SECTION-out Gaussian held-out log-likelihood,
     which Section 24.6 requires and Phase 3 did not have (AIC only).
"""
from __future__ import annotations

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import sasp_kernels as K

RES5 = "/workspace/results/phase5"
os.makedirs(RES5, exist_ok=True)

TRUNC_MULT = 6.0          # truncate the superposition sum at 6*lam_max
CHUNK = 20000


# ---------------------------------------------------------------------------
# Section 6.3 regressor
# ---------------------------------------------------------------------------

def superposition_basis(coords: np.ndarray, sender_xy: np.ndarray,
                        lam_grid: np.ndarray, chunk: int = CHUNK) -> np.ndarray:
    """S_i(lam) = sum_{j in S} exp(-d_ij / lam) for every lam on the grid.

    Returns (n, n_lam) float64.  cKDTree only; the pair list is built one
    receiver chunk at a time so peak memory is bounded by `chunk` x (mean
    senders within 6*lam_max).
    """
    n = coords.shape[0]
    out = np.zeros((n, lam_grid.size))
    if sender_xy.shape[0] == 0:
        return out
    R = TRUNC_MULT * float(lam_grid.max())
    ts = cKDTree(sender_xy)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        tc = cKDTree(coords[a:b])
        D = tc.sparse_distance_matrix(ts, R, output_type="coo_matrix")
        row, dat = D.row, D.data
        if row.size == 0:
            continue
        for g, lam in enumerate(lam_grid):
            out[a:b, g] = np.bincount(row, weights=np.exp(-dat / lam),
                                      minlength=b - a)
    return out


def nearest_basis(d: np.ndarray, lam_grid: np.ndarray) -> np.ndarray:
    """K_i(lam) = exp(-d_i / lam), the Section 6.1 nearest-sender regressor."""
    return np.exp(-d[:, None] / lam_grid[None, :])


# ---------------------------------------------------------------------------
# profiling one basis under a fixed design
# ---------------------------------------------------------------------------

def profile_basis(Kb: np.ndarray, X: np.ndarray, y: np.ndarray,
                  bid: np.ndarray, nb: int, lam_grid: np.ndarray,
                  boot_m: np.ndarray = None) -> dict:
    """Profile a precomputed (n, G) basis against design X.

    Returns lam_hat, beta_hat, rss, rss0 (covariates only), aic, and -- when
    `boot_m` (n_boot, nb) is given -- the block-bootstrap distribution of the
    profiled rss.

    `boot_m` is passed IN rather than drawn here on purpose.  Two candidate
    bases must be evaluated on the SAME resample of blocks or the difference
    between them is pure resampling noise and the paired win fraction collapses
    to 0.5 by construction.
    """
    prof = K.BasisBlockProfiler(Kb, X, y, bid, nb)
    one = np.ones(nb)
    rss, beta, n = prof.profile(one)
    t = int(np.nanargmin(rss))
    rss0, _ = prof.fit_linear(one)
    p = X.shape[1]
    k_par = p + 1 + 1                      # covariates + beta + lambda
    out = dict(lam=float(lam_grid[t]), t=t, beta=float(beta[t]),
               rss=float(rss[t]), rss0=float(rss0), n=float(n),
               k_params=k_par,
               aic=float(n * np.log(rss[t] / n) + 2 * (k_par + 1)),
               aic0=float(n * np.log(rss0 / n) + 2 * (p + 1)),
               railed=int(t == 0 or t == lam_grid.size - 1))
    out["prof"] = prof
    if boot_m is not None:
        n_boot = boot_m.shape[0]
        rr = np.full(n_boot, np.nan)
        ll = np.full(n_boot, np.nan)
        for r in range(n_boot):
            m = boot_m[r]
            rs, _, nn = prof.profile(m)
            if np.all(~np.isfinite(rs)):
                continue
            j = int(np.nanargmin(rs))
            rr[r] = rs[j] / nn            # per-cell RSS, comparable across reps
            ll[r] = lam_grid[j]
        out["boot_rss_per_cell"] = rr
        out["boot_lam"] = ll
    return out


# ---------------------------------------------------------------------------
# Section 24.6 — held-out log-likelihood on left-out sections
# ---------------------------------------------------------------------------

def _prof_coef(prof, m, g):
    """(theta, beta, rss, n) of  y ~ X theta + beta K[:, g]  under block
    multiplicities `m`, read straight off the precomputed block Grams."""
    XX = np.tensordot(m, prof.XX, 1)
    XK = np.tensordot(m, prof.XK, 1)
    Xy = np.tensordot(m, prof.Xy, 1)
    Gi = np.linalg.pinv(XX)
    GX = Gi @ Xy
    A = Gi @ XK
    yy = float(m @ prof.yy) - float(Xy @ GX)
    kk = float(np.tensordot(m, prof.KK, 1)[g] - XK[:, g] @ A[:, g])
    ky = float(np.tensordot(m, prof.Ky, 1)[g] - XK[:, g] @ GX)
    beta = ky / kk if kk > 1e-12 else 0.0
    theta = GX - A[:, g] * beta
    rss = yy - (ky ** 2 / kk if kk > 1e-12 else 0.0)
    return theta, beta, rss, float(m @ prof.nn)


def heldout_ll(prof, parts: list, lam_grid: np.ndarray, basis_key: str,
               prof0=None) -> list:
    """Leave-one-SECTION-out Gaussian held-out log-likelihood (Section 24.6).

    `prof` is a `BasisBlockProfiler` whose BLOCKS ARE SECTIONS, so a training
    fold is just an indicator multiplicity vector and no refit on stacked data
    is ever needed.  `parts` carries the per-section test arrays.

    Both candidate bases must be z-scored within section before the profiler is
    built.  That is not cosmetic: the superposition sum scales with the
    section's sender density, so an unstandardised amplitude is not
    transferable between sections and the held-out comparison would be a
    density comparison rather than a kernel comparison.
    """
    S = len(parts)
    rows = []
    for h in range(S):
        m = np.ones(S)
        m[h] = 0.0
        rss = np.full(lam_grid.size, np.inf)
        for g in range(lam_grid.size):
            rss[g] = _prof_coef(prof, m, g)[2]
        g = int(np.nanargmin(rss))
        theta, beta, rss_tr, ntr = _prof_coef(prof, m, g)
        p = theta.size
        s2 = rss_tr / max(ntr - p - 2, 1)
        te = parts[h]
        pred = te["X"] @ theta + beta * te[basis_key][:, g]
        ll = K.gaussian_heldout_ll(te["y"], pred, s2)
        r = dict(held_out=te["name"], lam=float(lam_grid[g]), beta=float(beta),
                 n_test=int(te["y"].size), ll=ll,
                 ll_per_cell=ll / te["y"].size, rss_tr_per_cell=rss_tr / ntr,
                 rss_te_per_cell=float(((te["y"] - pred) ** 2).mean()))
        if prof0 is not None:
            XX = np.tensordot(m, prof0.XX, 1)
            Xy = np.tensordot(m, prof0.Xy, 1)
            th0 = np.linalg.pinv(XX) @ Xy
            n0 = float(m @ prof0.nn)
            rss0 = float(m @ prof0.yy) - float(Xy @ th0)
            s20 = rss0 / max(n0 - th0.size, 1)
            ll0 = K.gaussian_heldout_ll(te["y"], te["X"] @ th0, s20)
            r["ll0"] = ll0
            r["ll0_per_cell"] = ll0 / te["y"].size
            r["dll_vs_cov"] = ll - ll0
        rows.append(r)
    return rows


def zsc(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    s = v.std()
    return (v - v.mean()) / (s if s > 1e-12 else 1.0)


def zsc_cols(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, float)
    mu = M.mean(0)
    sd = M.std(0)
    sd[sd < 1e-12] = 1.0
    return (M - mu) / sd


# ---------------------------------------------------------------------------
# nulls, reused for T1 stability and T4
# ---------------------------------------------------------------------------

def torus_shift(rng, pts, lo, hi):
    span = hi - lo
    return lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span


def permute_within_type(rng, sender, celltype, eligible):
    out = np.zeros_like(sender)
    for c in np.unique(celltype):
        ix = np.flatnonzero((celltype == c) & eligible)
        k = int(sender[ix].sum())
        if k:
            out[rng.choice(ix, size=k, replace=False)] = True
    return out
