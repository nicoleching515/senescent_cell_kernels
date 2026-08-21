"""
Nulls N1, N3, N4 on synthetic tissue (Master Plan Section 23).

N1 -- cell-type-stratified label permutation.  Reassign sender labels at random
      among cells of the SAME cell type.  Preserves composition and
      architecture; destroys sender-specific signal.
N3 -- torus shift.  Translate the sender coordinate set by a random vector with
      wraparound.  Preserves sender clustering AND receiver autocorrelation;
      destroys their alignment.  This is the null CellWHISPER showed leading
      methods fail, so it is the one that matters most.
N4 -- rotation.  Rotate sender coordinates about the tissue centroid.

Section 18.2: N1 changes no geometry that a re-index cannot express, but the
distance-to-nearest-sender DOES change when the sender set changes, so N1 still
requires a tree query per permutation -- what it does NOT require is rebuilding
the tree over all cells.  We build the all-cell tree once and re-query it.
N3 and N4 genuinely move the senders and are the expensive ones, so they get
the parallelism (Section 18.3).

What we measure, which the plan does not specify but Phase 1 makes essential:
  * SIZE   -- with beta_true = 0, what fraction of nulls reject at alpha = 0.05?
              A null that rejects far more than 5% is not a valid null.
  * POWER  -- with a real planted kernel, what fraction reject?
  * the null distribution of lambda_hat and beta_hat, not just a p-value
    (Section 24.3).
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from scipy.spatial import cKDTree

from sasp_estimators import BlockProfiler, block_ids, LAM_S_GRID, design_matrix


def _stat(coords, sender_mask, y_idx, y, X0, bid, nb, lam_grid):
    """beta_hat and lambda_hat of the exponential kernel for a given sender set."""
    if sender_mask.sum() < 5:
        return np.nan, np.nan
    d = cKDTree(coords[sender_mask]).query(coords[y_idx], k=1, workers=1)[0]
    K = np.exp(-d[:, None] / lam_grid[None, :])
    from sasp_kernels import BasisBlockProfiler
    prof = BasisBlockProfiler(K, X0, y, bid, nb)
    rss, beta, _ = prof.profile(np.ones(nb))
    t = int(np.nanargmin(rss))
    return float(beta[t]), float(lam_grid[t])


def permute_within_type(rng, sender_mask, cell_type):
    """N1: reassign sender labels at random within each cell type."""
    out = np.zeros_like(sender_mask)
    for c in np.unique(cell_type):
        ix = np.flatnonzero(cell_type == c)
        k = int(sender_mask[ix].sum())
        if k:
            out[rng.choice(ix, size=k, replace=False)] = True
    return out


def torus_shift(rng, sender_coords, window):
    """N3: translate the sender set with wraparound."""
    v = rng.uniform(0, window, size=2)
    return (sender_coords + v) % window


def rotate(rng, sender_coords, window):
    """N4: rotate the sender set about the tissue centroid, wrapped."""
    th = rng.uniform(0, 2 * np.pi)
    c = np.array([window / 2, window / 2])
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    return ((sender_coords - c) @ R.T + c) % window


def run_nulls(sim: Dict, seed, n_perm: int = 200,
              nulls=("N1", "N3", "N4")) -> Dict[str, object]:
    """Observed statistic plus the full null distribution for each null."""
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    coords = sim["coords"]
    W = sim["window_um"]
    smask = sim["sender_mask"]
    ct = sim["cell_type"]

    idx = np.flatnonzero(~smask)
    y = sim["r"][idx]
    X0 = design_matrix(sim, idx, with_nuisance=False)
    bid = block_ids(coords[idx], W, n_blocks_side=8)
    nb = 64
    lam_grid = LAM_S_GRID

    beta_obs, lam_obs = _stat(coords, smask, idx, y, X0, bid, nb, lam_grid)
    # Statistic 2, evaluated at the OBSERVED lambda_hat and held fixed.
    # Re-profiling lambda inside the null compares two DIFFERENT models: when a
    # null destroys the spatial structure the profiled lambda rails to the top
    # of the grid, the kernel column becomes nearly constant, and beta stops
    # being an amplitude at all.  Holding lambda at lambda_obs keeps the null
    # and the observed statistic on the same scale, which is what a permutation
    # test requires.
    fixed = np.array([lam_obs])
    out = dict(beta_obs=beta_obs, lam_obs=lam_obs,
               beta_true=sim["cfg"].beta_true,
               lam_true=sim["cfg"].lambda_true_um)

    beta_obs_fix, _ = _stat(coords, smask, idx, y, X0, bid, nb, fixed)
    out["beta_obs_fixedlam"] = beta_obs_fix

    send_xy = coords[smask]
    for null in nulls:
        bs = np.full(n_perm, np.nan)
        ls = np.full(n_perm, np.nan)
        bf = np.full(n_perm, np.nan)
        for t in range(n_perm):
            if null == "N1":
                m = permute_within_type(rng, smask, ct)
                b, l = _stat(coords, m, idx, y, X0, bid, nb, lam_grid)
                bf[t] = _stat(coords, m, idx, y, X0, bid, nb, fixed)[0]
            else:
                new_xy = (torus_shift(rng, send_xy, W) if null == "N3"
                          else rotate(rng, send_xy, W))
                d = cKDTree(new_xy).query(coords[idx], k=1, workers=1)[0]
                from sasp_kernels import BasisBlockProfiler
                Kk = np.exp(-d[:, None] / lam_grid[None, :])
                prof = BasisBlockProfiler(Kk, X0, y, bid, nb)
                rss, beta, _ = prof.profile(np.ones(nb))
                j = int(np.nanargmin(rss))
                b, l = float(beta[j]), float(lam_grid[j])
                kf = np.exp(-d[:, None] / fixed[None, :])
                pf = BasisBlockProfiler(kf, X0, y, bid, nb)
                _, bfix, _ = pf.profile(np.ones(nb))
                bf[t] = float(bfix[0])
            bs[t], ls[t] = b, l
        ok = np.isfinite(bs)
        # one-sided empirical p-value on beta (larger = more signal)
        p = (1.0 + np.sum(bs[ok] >= beta_obs)) / (1.0 + ok.sum())
        out[f"{null}_p_beta"] = float(p)
        out[f"{null}_null_beta_mean"] = float(np.nanmean(bs))
        out[f"{null}_null_beta_q95"] = float(np.nanquantile(bs[ok], 0.95))
        out[f"{null}_null_lam_mean"] = float(np.nanmean(ls))
        out[f"{null}_reject05"] = float(p <= 0.05)
        okf = np.isfinite(bf)
        pf_ = (1.0 + np.sum(bf[okf] >= beta_obs_fix)) / (1.0 + okf.sum())
        out[f"{null}_p_beta_fixedlam"] = float(pf_)
        out[f"{null}_reject05_fixedlam"] = float(pf_ <= 0.05)
        out[f"{null}_null_beta_fixedlam_mean"] = float(np.nanmean(bf))
        out[f"{null}_surviving_frac_fixedlam"] = (
            float((beta_obs_fix - np.nanmean(bf)) / beta_obs_fix)
            if np.isfinite(beta_obs_fix) and beta_obs_fix != 0 else np.nan)
        # surviving fraction (Section 6.5): observed minus the null mean
        out[f"{null}_surviving_beta"] = float(beta_obs - np.nanmean(bs))
        out[f"{null}_surviving_frac"] = (
            float((beta_obs - np.nanmean(bs)) / beta_obs)
            if beta_obs not in (0, np.nan) and np.isfinite(beta_obs) and beta_obs != 0
            else np.nan)
    return out
