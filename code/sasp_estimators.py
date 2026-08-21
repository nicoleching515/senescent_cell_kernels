"""
Estimators for the SASP spatial response kernel (Master Plan Sections 6, 22, 23).

Implemented here
----------------
1. NAIVE estimator (what the field currently does, Section 22 Step 2):
     bin cells by distance-to-nearest-sender in 10 um bins out to 300 um,
     fit r = a + beta * exp(-d/lambda), report lambda_hat, beta_hat.
   Also a cell-level version, which is what gives the (over-confident)
   asymptotic CI that the coverage panel of Figure 1 is about.

2. MATCHED-DECOY estimator (null N2, Section 23 -- "the single most important
   number in the paper"):  for each sender pick a non-sender matched on cell
   type, local density (kNN count within 50 um) and kNN cell-type composition,
   by propensity-score nearest-neighbour matching (coarsened exact matching is
   also provided).  Recompute the kernel using decoys as senders.  The
   corrected effect is the sender curve MINUS the decoy curve.

3. NUISANCE-CONDITIONED fitting (null N5): the same cell-level fit with the
   Tier-D style covariates available in the synthetic data.

4. Standardized mean difference (SMD) on the matching covariates before and
   after matching.  Section 8 Test 5 sets the bar at SMD < 0.1.

Uncertainty
-----------
* asymptotic Gauss-Newton CI assuming iid errors  -> the "naive" CI
* spatial block bootstrap CI                      -> the honest alternative
Both are reported so the coverage failure can be decomposed into bias versus
understated standard error.

Compute rules (Section 18): cKDTree for all neighbour work; the geometry
(d_sender, d_decoy, covariates) is computed ONCE per tissue and every
bootstrap replicate is a re-indexing, never a re-query.
"""

from __future__ import annotations

from typing import Dict, Any, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

LAM_LO, LAM_HI = 3.0, 400.0
LAM_GRID = np.exp(np.linspace(np.log(LAM_LO), np.log(LAM_HI), 240))


# --------------------------------------------------------------------------
# core: profiled exponential-kernel least squares
# --------------------------------------------------------------------------


def _wls(X: np.ndarray, y: np.ndarray, w: Optional[np.ndarray]):
    if w is None:
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
        return coef, float(resid @ resid)
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    resid = y - X @ coef
    return coef, float((w * resid ** 2).sum())


def _profile_rss(d: np.ndarray, y: np.ndarray, X0: np.ndarray,
                 w: Optional[np.ndarray], lam_grid: np.ndarray) -> np.ndarray:
    """RSS(lambda) for the whole grid at once, by projecting out the linear
    block X0 and then solving the remaining 1-D problem in closed form.

    This is the compute-small trick of Section 18 applied to the fit itself:
    two BLAS-3 matmuls instead of `len(lam_grid)` separate least-squares
    solves.  It turns a ~0.7 s fit into a ~5 ms fit, which is what makes a
    sweep with block bootstraps affordable.
    """
    n = d.shape[0]
    sw = None if w is None else np.sqrt(w)
    Xw = X0 if sw is None else X0 * sw[:, None]
    yw = y if sw is None else y * sw
    # projector onto X0 (weighted)
    G = Xw.T @ Xw
    Ginv = np.linalg.pinv(G)
    Py = Ginv @ (Xw.T @ yw)
    y_res = yw - Xw @ Py
    yy = float(y_res @ y_res)

    out = np.empty(lam_grid.size)
    # chunk the lambda grid so the (n, n_lam) block stays modest
    chunk = max(1, int(4e7 // max(n, 1)))
    for a in range(0, lam_grid.size, chunk):
        lg = lam_grid[a:a + chunk]
        K = np.exp(-d[:, None] / lg[None, :])
        Kw = K if sw is None else K * sw[:, None]
        K_res = Kw - Xw @ (Ginv @ (Xw.T @ Kw))
        kk = np.einsum("ij,ij->j", K_res, K_res)
        ky = K_res.T @ y_res
        with np.errstate(divide="ignore", invalid="ignore"):
            out[a:a + chunk] = yy - np.where(kk > 1e-12, ky ** 2 / kk, 0.0)
    return out


def fit_exp_kernel(d: np.ndarray, y: np.ndarray,
                   X0: Optional[np.ndarray] = None,
                   w: Optional[np.ndarray] = None,
                   lam_grid: np.ndarray = LAM_GRID,
                   refine: bool = True) -> Dict[str, Any]:
    """Fit y = X0 @ theta + beta * exp(-d/lambda) (+ error).

    lambda is profiled out on a log grid (the linear parameters are solved
    exactly at each lambda), then parabolically refined.  This is far more
    robust than handing the whole thing to a generic optimiser, which
    routinely finds local minima on confounded data.

    Returns lambda_hat, beta_hat, asymptotic SEs (iid errors), and flags.
    """
    n = d.shape[0]
    if X0 is None:
        X0 = np.ones((n, 1))
    p0 = X0.shape[1]

    rss = _profile_rss(d, y, X0, w, lam_grid)
    t_best = int(np.argmin(rss))
    lam_hat = float(lam_grid[t_best])

    # parabolic refinement in log-lambda
    if refine and 0 < t_best < lam_grid.size - 1:
        x1, x2, x3 = np.log(lam_grid[t_best - 1:t_best + 2])
        y1, y2, y3 = rss[t_best - 1:t_best + 2]
        denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
        if abs(denom) > 0:
            a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
            b = (x3 ** 2 * (y1 - y2) + x2 ** 2 * (y3 - y1)
                 + x1 ** 2 * (y2 - y3)) / denom
            if a > 0:
                xv = -b / (2 * a)
                if x1 <= xv <= x3:
                    lam_hat = float(np.exp(xv))

    k = np.exp(-d / lam_hat)
    X = np.column_stack([X0, k])
    coef, rss_hat = _wls(X, y, w)
    beta_hat = float(coef[-1])

    # Gauss-Newton asymptotic covariance, iid errors.
    # d/dlambda of beta*exp(-d/lam) = beta * d / lam^2 * exp(-d/lam)
    J = np.column_stack([X0, k, beta_hat * d / lam_hat ** 2 * k])
    dof = max(n - J.shape[1], 1)
    sigma2 = rss_hat / dof
    JtJ = J.T @ J if w is None else J.T @ (J * w[:, None])
    try:
        cov = sigma2 * np.linalg.pinv(JtJ)
        se_beta = float(np.sqrt(max(cov[p0, p0], 0.0)))
        se_lam = float(np.sqrt(max(cov[p0 + 1, p0 + 1], 0.0)))
    except np.linalg.LinAlgError:
        se_beta = se_lam = np.nan

    at_bound = bool(lam_hat <= LAM_LO * 1.02 or lam_hat >= LAM_HI * 0.98)
    return dict(lam=lam_hat, beta=beta_hat, se_lam=se_lam, se_beta=se_beta,
                rss=rss_hat, sigma2=sigma2, at_bound=at_bound, n=n,
                coef=coef)


# --------------------------------------------------------------------------
# binning (Section 22 Step 2: 10 um bins out to 300 um)
# --------------------------------------------------------------------------

BIN_EDGES = np.arange(0.0, 310.0, 10.0)
BIN_MID = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
N_BINS = BIN_MID.size


def bin_index(d: np.ndarray) -> np.ndarray:
    """Bin id per cell, -1 for cells beyond 300 um."""
    b = np.floor(d / 10.0).astype(np.int64)
    b[(d >= BIN_EDGES[-1]) | (d < 0)] = -1
    return b


def binned_means(bidx: np.ndarray, y: np.ndarray, d: np.ndarray,
                 sub: Optional[np.ndarray] = None):
    """Mean of y and mean of d per distance bin.

    `sub` is an index array (a bootstrap resample); passing it makes a
    bootstrap replicate a pure re-indexing of precomputed geometry rather than
    a re-query of the KD-tree (Section 18.2).

    The mean distance inside each bin is returned and used as the regressor
    instead of the bin midpoint.  Using midpoints introduces a purely
    discretisation-driven upward bias in lambda_hat (Jensen, plus the hard-core
    process leaves the 0-10 um bin nearly empty), and we do not want that
    artefact contaminating a figure about confounding bias.
    """
    if sub is not None:
        bidx = bidx[sub]
        y = y[sub]
        d = d[sub]
    ok = bidx >= 0
    bb = bidx[ok]
    cnt = np.bincount(bb, minlength=N_BINS).astype(float)
    tot = np.bincount(bb, weights=y[ok], minlength=N_BINS)
    dtot = np.bincount(bb, weights=d[ok], minlength=N_BINS)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = tot / cnt
        dm = dtot / cnt
    return m, cnt, dm


def fit_binned_curve(m: np.ndarray, cnt: np.ndarray,
                     dm: Optional[np.ndarray] = None) -> Dict[str, Any]:
    x = BIN_MID if dm is None else np.where(np.isfinite(dm), dm, BIN_MID)
    ok = np.isfinite(m) & (cnt >= 5)
    if ok.sum() < 6:
        return dict(lam=np.nan, beta=np.nan, se_lam=np.nan, se_beta=np.nan,
                    at_bound=True, n=int(ok.sum()))
    return fit_exp_kernel(x[ok], m[ok],
                          X0=np.ones((int(ok.sum()), 1)), w=cnt[ok])


# --------------------------------------------------------------------------
# covariate design (Tier D analogue)
# --------------------------------------------------------------------------


def design_matrix(sim: Dict[str, Any], idx: np.ndarray,
                  with_nuisance: bool) -> np.ndarray:
    """Intercept + cell-type dummies (always: mu_{c_i} is part of the model),
    plus, if requested, the observable nuisance covariates z_i."""
    ct = sim["cell_type"][idx]
    n_types = int(sim["cfg"].n_types)
    cols = [np.ones(idx.size)]
    for k in range(1, n_types):
        cols.append((ct == k).astype(float))
    if with_nuisance:
        dens = sim["dens50"][idx]
        cols.append((dens - dens.mean()) / (dens.std() + 1e-12))
        comp = sim["knn_comp"][idx]
        for k in range(n_types - 1):
            c = comp[:, k]
            cols.append((c - c.mean()) / (c.std() + 1e-12))
        cols.append(sim["log_counts_z"][idx])
        de = sim["d_edge"][idx]
        cols.append((de - de.mean()) / (de.std() + 1e-12))
    return np.column_stack(cols)


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ coef + float(y.mean())


# --------------------------------------------------------------------------
# matched decoys (null N2 / Section 8 Test 5)
# --------------------------------------------------------------------------


def _logit(p, eps=1e-6):
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def match_decoys(sim: Dict[str, Any], rng: np.random.Generator,
                 method: str = "propensity",
                 caliper_sd: float = 0.25) -> Dict[str, Any]:
    """For each sender, select a non-sender of the SAME cell type matched on
    local density (kNN count within 50 um) and kNN cell-type composition.

    method="propensity": logistic propensity score, greedy 1-1 nearest-
      neighbour matching without replacement inside each cell type, with a
      caliper.  Matching is done in 1-D propensity space by sorting +
      searchsorted, so it is O(n log n) and never forms an (n, n) matrix.
    method="cem": coarsened exact matching on quantile-binned covariates.
    """
    n = sim["n_cells"]
    s = sim["sender_mask"]
    ct = sim["cell_type"]
    n_types = int(sim["cfg"].n_types)

    dens = sim["dens50"]
    Z = np.column_stack([
        (dens - dens.mean()) / (dens.std() + 1e-12),
        sim["knn_comp"][:, :n_types - 1],
        sim["log_counts_z"],
    ])
    cov_names = (["dens50"]
                 + [f"knn_comp_t{k}" for k in range(n_types - 1)]
                 + ["log_counts"])

    if method == "propensity":
        from sklearn.linear_model import LogisticRegression
        Xp = np.column_stack([Z] + [(ct == k).astype(float)
                                    for k in range(1, n_types)])
        lr = LogisticRegression(max_iter=400, C=1.0)
        lr.fit(Xp, s.astype(int))
        ps = _logit(lr.predict_proba(Xp)[:, 1])
    else:
        # coarsened exact matching: 5 quantile bins per covariate
        strata = np.zeros(n, dtype=np.int64)
        for j in range(Z.shape[1]):
            q = np.quantile(Z[:, j], np.linspace(0, 1, 6)[1:-1])
            strata = strata * 5 + np.searchsorted(q, Z[:, j])
        ps = strata.astype(float)

    decoy = np.full(int(s.sum()), -1, dtype=np.int64)
    send_ids = np.flatnonzero(s)
    used = np.zeros(n, dtype=bool)

    for k in range(n_types):
        s_k = np.flatnonzero(s & (ct == k))
        c_k = np.flatnonzero((~s) & (ct == k))
        if s_k.size == 0 or c_k.size == 0:
            continue
        order = np.argsort(ps[c_k], kind="stable")
        c_sorted = c_k[order]
        v_sorted = ps[c_sorted]
        cal = caliper_sd * float(np.std(ps[ct == k]) + 1e-12)
        # process senders in random order to avoid systematic greedy bias
        for si in rng.permutation(s_k):
            pos = int(np.searchsorted(v_sorted, ps[si]))
            best, best_d = -1, np.inf
            lo, hi = pos - 1, pos
            # expand outward until both sides exceed the current best distance
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

    matched = decoy >= 0
    decoy_idx = decoy[matched]
    sender_matched = send_ids[matched]

    smd_before = _smd(Z, s, ~s)
    m_s = np.zeros(n, bool); m_s[sender_matched] = True
    m_d = np.zeros(n, bool); m_d[decoy_idx] = True
    smd_after = _smd(Z, m_s, m_d)

    return dict(decoy_idx=decoy_idx, sender_matched=sender_matched,
                match_rate=float(matched.mean()),
                smd_before=smd_before, smd_after=smd_after,
                cov_names=cov_names,
                max_smd_after=float(np.nanmax(np.abs(smd_after))),
                ps=ps)


def _smd(Z: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Standardized mean difference per covariate (Section 8 Test 5)."""
    ma, mb = Z[a].mean(axis=0), Z[b].mean(axis=0)
    sa, sb = Z[a].std(axis=0), Z[b].std(axis=0)
    pooled = np.sqrt(0.5 * (sa ** 2 + sb ** 2)) + 1e-12
    return (ma - mb) / pooled


# --------------------------------------------------------------------------
# block-Gram profiler: the whole estimator suite in sufficient-statistic form
# --------------------------------------------------------------------------
#
# Section 18.2 says a permutation/bootstrap should be a re-indexing, not a
# re-query.  We take that one step further.  The spatial block bootstrap
# resamples whole blocks, so a replicate is fully described by an integer
# multiplicity m_b per block.  Every quantity any of our fits needs is a
# cross-product, and a cross-product over a weighted union of blocks is just
# sum_b m_b * (cross-product within block b).
#
# So we compute per-block cross-products ONCE (cost: one pass over the data)
# and every bootstrap replicate becomes a weighted sum of 64 small matrices
# plus a vectorised 2-D grid solve.  A 400-replicate block bootstrap of a
# two-kernel profiled nonlinear fit then costs ~0.2 s instead of ~30 s.


class BlockProfiler:
    """Precomputed sufficient statistics for

        y = X[:, :p] theta + beta_s K(d_s; lam_s) [+ beta_d K(d_d; lam_d)] + eps

    over spatial blocks, supporting arbitrary block multiplicities.
    """

    def __init__(self, d_s, d_d, y, X, bid, n_blocks, lam_s, lam_d):
        self.lam_s, self.lam_d = lam_s, lam_d
        self.d_s = d_s
        self.have_d = d_d is not None and np.all(np.isfinite(d_d))
        ns, nd = lam_s.size, (lam_d.size if self.have_d else 1)
        self.ns, self.nd, self.B = ns, nd, n_blocks
        p = X.shape[1]
        self.p_full = p

        Ks = np.exp(-d_s[:, None] / lam_s[None, :])
        Kd = (np.exp(-d_d[:, None] / lam_d[None, :]) if self.have_d
              else np.zeros((d_s.size, 1)))
        self.Ks_full, self.Kd_full = Ks, Kd

        self.XX = np.zeros((n_blocks, p, p))
        self.XKs = np.zeros((n_blocks, p, ns))
        self.XKd = np.zeros((n_blocks, p, nd))
        self.Xy = np.zeros((n_blocks, p))
        self.KsKs = np.zeros((n_blocks, ns))
        self.KdKd = np.zeros((n_blocks, nd))
        self.KsKd = np.zeros((n_blocks, ns, nd))
        self.Ksy = np.zeros((n_blocks, ns))
        self.Kdy = np.zeros((n_blocks, nd))
        self.yy = np.zeros(n_blocks)
        self.nb_cells = np.zeros(n_blocks)

        order = np.argsort(bid, kind="stable")
        starts = np.searchsorted(bid[order], np.arange(n_blocks + 1))
        for b in range(n_blocks):
            sl = order[starts[b]:starts[b + 1]]
            if sl.size == 0:
                continue
            Xb, Ksb, Kdb, yb = X[sl], Ks[sl], Kd[sl], y[sl]
            self.XX[b] = Xb.T @ Xb
            self.XKs[b] = Xb.T @ Ksb
            self.XKd[b] = Xb.T @ Kdb
            self.Xy[b] = Xb.T @ yb
            self.KsKs[b] = np.einsum("ij,ij->j", Ksb, Ksb)
            self.KdKd[b] = np.einsum("ij,ij->j", Kdb, Kdb)
            self.KsKd[b] = Ksb.T @ Kdb
            self.Ksy[b] = Ksb.T @ yb
            self.Kdy[b] = Kdb.T @ yb
            self.yy[b] = float(yb @ yb)
            self.nb_cells[b] = sl.size

    # -- weighted accumulation -------------------------------------------
    def _acc(self, m, p):
        XX = np.tensordot(m, self.XX, 1)[:p, :p]
        XKs = np.tensordot(m, self.XKs, 1)[:p]
        XKd = np.tensordot(m, self.XKd, 1)[:p]
        Xy = np.tensordot(m, self.Xy, 1)[:p]
        Gi = np.linalg.pinv(XX)
        GX = Gi @ Xy
        yy = float(m @ self.yy) - float(Xy @ GX)
        A = Gi @ XKs                                    # (p, ns)
        kk = np.tensordot(m, self.KsKs, 1) - np.einsum("ij,ij->j", XKs, A)
        ky = np.tensordot(m, self.Ksy, 1) - XKs.T @ GX
        Ad = Gi @ XKd
        dd = np.tensordot(m, self.KdKd, 1) - np.einsum("ij,ij->j", XKd, Ad)
        dy = np.tensordot(m, self.Kdy, 1) - XKd.T @ GX
        sd = np.tensordot(m, self.KsKd, 1) - XKs.T @ Ad
        n = float(m @ self.nb_cells)
        return yy, kk, ky, dd, dy, sd, n

    # -- single kernel (naive / nuisance-conditioned) ---------------------
    def fit1(self, m, p, refine=True):
        yy, kk, ky, _, _, _, n = self._acc(m, p)
        with np.errstate(divide="ignore", invalid="ignore"):
            rss = yy - np.where(kk > 1e-12, ky ** 2 / kk, 0.0)
        t = int(np.nanargmin(rss))
        lam = _refine(self.lam_s, rss, t) if refine else float(self.lam_s[t])
        beta = float(ky[t] / kk[t]) if kk[t] > 1e-12 else np.nan
        return dict(lam=lam, beta=beta, rss=float(rss[t]), n=n,
                    sigma2=float(rss[t]) / max(n - p - 2, 1), t=t)

    # -- two kernels (decoy-adjusted, null N2) ----------------------------
    def fit2(self, m, p, refine=True):
        if not self.have_d:
            r = self.fit1(m, p, refine)
            r.update(lam_d=np.nan, beta_d=np.nan)
            return r
        yy, kk, ky, dd, dy, sd, n = self._acc(m, p)
        det = kk[:, None] * dd[None, :] - sd ** 2
        good = det > 1e-10 * np.abs(kk[:, None] * dd[None, :])
        with np.errstate(divide="ignore", invalid="ignore"):
            bs = (ky[:, None] * dd[None, :] - sd * dy[None, :]) / det
            bd = (kk[:, None] * dy[None, :] - sd * ky[:, None]) / det
            rss = yy - (bs * ky[:, None] + bd * dy[None, :])
        rss = np.where(good & np.isfinite(rss), rss, np.inf)
        i, j = np.unravel_index(int(np.argmin(rss)), rss.shape)
        lam = _refine(self.lam_s, rss[:, j], i) if refine else float(self.lam_s[i])
        return dict(lam=lam, beta=float(bs[i, j]), lam_d=float(self.lam_d[j]),
                    beta_d=float(bd[i, j]), rss=float(rss[i, j]), n=n,
                    sigma2=float(rss[i, j]) / max(n - p - 4, 1), t=i, td=j)


    # -- fixed-lambda readouts (Phase 3) ----------------------------------
    # Added in Phase 3, additive only: every null must be evaluated at the
    # OBSERVED lambda, because re-profiling inside a null compares different
    # models and beta stops being an amplitude (CS Phase 2 Sec 10).
    def beta_at(self, m, p, t):
        """beta of the single-kernel fit at grid index t (lambda held fixed)."""
        yy, kk, ky, _, _, _, n = self._acc(m, p)
        if not np.isfinite(kk[t]) or kk[t] <= 1e-12:
            return np.nan, np.nan
        b = float(ky[t] / kk[t])
        return b, float(yy - b * ky[t])

    def beta2_at(self, m, p, t):
        """beta_sender of the SHARED-lambda two-kernel (matched-decoy) fit at
        grid index t.  Returns (beta_sender, beta_decoy, rss)."""
        if not self.have_d or self.lam_d.size != self.lam_s.size:
            b, r = self.beta_at(m, p, t)
            return b, np.nan, r
        yy, kk, ky, dd, dy, sd, n = self._acc(m, p)
        sdg = float(sd[t, t]); k = float(kk[t]); d = float(dd[t])
        det = k * d - sdg ** 2
        if not np.isfinite(det) or abs(det) <= 1e-10 * abs(k * d):
            return np.nan, np.nan, np.nan
        bs = (ky[t] * d - sdg * dy[t]) / det
        bd = (k * dy[t] - sdg * ky[t]) / det
        return float(bs), float(bd), float(yy - bs * ky[t] - bd * dy[t])

    # -- two kernels with a SHARED length scale ---------------------------
    def fit2_shared(self, m, p, refine=True):
        """Decoy-adjusted fit forcing lambda_d == lambda_s.

        This is the apples-to-apples version of the plan's
        "beta_true - beta_decoy".  With lambda free in both terms the two fits
        can land on completely different length scales -- in a long-range
        confounded regime the decoy fit runs off to lambda_d ~ 300 um -- and
        subtracting a long-range amplitude from a short-range one is not a
        meaningful contrast.  Constraining the scales makes beta_s read as
        "sender effect over and above whatever appears around a matched
        non-sender at the same distance".
        """
        if not self.have_d or self.lam_d.size != self.lam_s.size:
            r = self.fit1(m, p, refine)
            r.update(lam_d=np.nan, beta_d=np.nan)
            return r
        yy, kk, ky, dd, dy, sd, n = self._acc(m, p)
        sdg = np.diagonal(sd).copy()
        det = kk * dd - sdg ** 2
        good = det > 1e-10 * np.abs(kk * dd)
        with np.errstate(divide="ignore", invalid="ignore"):
            bs = (ky * dd - sdg * dy) / det
            bd = (kk * dy - sdg * ky) / det
            rss = yy - (bs * ky + bd * dy)
        rss = np.where(good & np.isfinite(rss), rss, np.inf)
        t = int(np.argmin(rss))
        lam = _refine(self.lam_s, rss, t) if refine else float(self.lam_s[t])
        return dict(lam=lam, beta=float(bs[t]), lam_d=lam, beta_d=float(bd[t]),
                    rss=float(rss[t]), n=n,
                    sigma2=float(rss[t]) / max(n - p - 3, 1), t=t, td=t)


def _refine(grid, rss, t):
    """Parabolic refinement of the profile minimum in log-lambda."""
    lam = float(grid[t])
    if 0 < t < grid.size - 1:
        x1, x2, x3 = np.log(grid[t - 1:t + 2])
        y1, y2, y3 = rss[t - 1:t + 2]
        den = (x1 - x2) * (x1 - x3) * (x2 - x3)
        if np.isfinite(y1) and np.isfinite(y3) and abs(den) > 0:
            A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / den
            B = (x3 ** 2 * (y1 - y2) + x2 ** 2 * (y3 - y1)
                 + x1 ** 2 * (y2 - y3)) / den
            if A > 0:
                xv = -B / (2 * A)
                if x1 <= xv <= x3:
                    lam = float(np.exp(xv))
    return lam


def _se_lambda(X, d_s, lam, beta, sigma2, extra=None):
    """Gauss-Newton asymptotic SE of lambda_hat under iid errors -- this is the
    CI a standard nonlinear-least-squares package would report, and the one the
    coverage panel of Figure 1 is about."""
    k = np.exp(-d_s / lam)
    cols = [X, k[:, None]]
    if extra is not None:
        cols.append(extra[:, None])
    cols.append((beta * d_s / lam ** 2 * k)[:, None])
    J = np.column_stack(cols)
    try:
        cov = sigma2 * np.linalg.pinv(J.T @ J)
        return float(np.sqrt(max(cov[-1, -1], 0.0))), \
               float(np.sqrt(max(cov[X.shape[1], X.shape[1]], 0.0)))
    except np.linalg.LinAlgError:
        return np.nan, np.nan


# --------------------------------------------------------------------------
# spatial block bootstrap
# --------------------------------------------------------------------------


def block_ids(coords: np.ndarray, window_um: float, n_blocks_side: int = 8):
    h = window_um / n_blocks_side
    bx = np.clip((coords[:, 0] / h).astype(np.int64), 0, n_blocks_side - 1)
    by = np.clip((coords[:, 1] / h).astype(np.int64), 0, n_blocks_side - 1)
    return bx * n_blocks_side + by


def _ci(vals: np.ndarray) -> Tuple[float, float]:
    v = vals[np.isfinite(vals)]
    if v.size < 10:
        return (np.nan, np.nan)
    return (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)))


def _cover(ci, truth) -> float:
    if not np.all(np.isfinite(ci)):
        return np.nan
    return float(ci[0] <= truth <= ci[1])


# --------------------------------------------------------------------------
# top-level analysis of one simulated tissue
# --------------------------------------------------------------------------

LAM_S_GRID = np.exp(np.linspace(np.log(LAM_LO), np.log(LAM_HI), 96))
LAM_D_GRID = LAM_S_GRID          # shared grid enables fit2_shared


def analyze(sim: Dict[str, Any], seed, n_boot: int = 400,
            match_method: str = "propensity",
            keep_curves: bool = False) -> Dict[str, Any]:
    """Run every estimator on one synthetic tissue and return a flat record.

    Estimators
      naive       cell-level single-kernel fit, cell-type intercepts only
      naive_bin   the binned-mean version the field actually plots
      nuis        + Tier-D nuisance covariates                        (N5)
      decoy       cell-level decoy-adjusted two-kernel fit            (N2)
      decoy_nuis  N2 + N5 together
    """
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    cfg = sim["cfg"]
    lam_true, beta_true = cfg.lambda_true_um, cfg.beta_true
    W = sim["window_um"]

    out: Dict[str, Any] = dict(
        n_cells=sim["n_cells"], n_senders=sim["n_senders"],
        prevalence_real=sim["prevalence_real"],
        median_nn_um=sim["median_nn_um"],
        lam_true=lam_true, beta_true=beta_true,
        med_d_sender=float(np.median(sim["d_sender"])),
    )

    # ---- matched decoys (N2 / Section 8 Test 5) ---------------------------
    mt = match_decoys(sim, rng, method=match_method)
    out["match_rate"] = mt["match_rate"]
    out["max_smd_before"] = float(np.nanmax(np.abs(mt["smd_before"])))
    out["max_smd_after"] = mt["max_smd_after"]
    out["smd_pass"] = float(mt["max_smd_after"] < 0.1)

    decoy_mask = np.zeros(sim["n_cells"], dtype=bool)
    decoy_mask[mt["decoy_idx"]] = True
    have_decoy = mt["decoy_idx"].size >= 20
    if have_decoy:
        d_decoy, _ = cKDTree(sim["coords"][decoy_mask]).query(
            sim["coords"], k=1, workers=1)
        # audit (latent u is NOT available to any estimator)
        out["u_sender_mean"] = float(sim["u_latent"][sim["sender_mask"]].mean())
        out["u_decoy_mean"] = float(sim["u_latent"][decoy_mask].mean())
        out["u_all_mean"] = float(sim["u_latent"].mean())
    else:
        d_decoy = None

    # ---- analysis set: exclude senders and decoys symmetrically -----------
    keep = ~sim["sender_mask"] & ~decoy_mask
    idx = np.flatnonzero(keep)
    y = sim["r"][idx]
    ds = sim["d_sender"][idx]
    dd = d_decoy[idx] if have_decoy else None

    X_nuis = design_matrix(sim, idx, with_nuisance=True)
    p_base = int(sim["cfg"].n_types)          # intercept + (n_types-1) dummies
    X_base = X_nuis[:, :p_base]
    p_nuis = X_nuis.shape[1]

    # ---- block profiler ---------------------------------------------------
    nbs = 8
    nb = nbs * nbs
    bid = block_ids(sim["coords"][idx], W, n_blocks_side=nbs)
    bp = BlockProfiler(ds, dd, y, X_nuis, bid, nb, LAM_S_GRID, LAM_D_GRID)
    ones = np.ones(nb)

    f_naive = bp.fit1(ones, p_base)
    f_nuis = bp.fit1(ones, p_nuis)
    f_dec = bp.fit2(ones, p_base)
    f_decn = bp.fit2(ones, p_nuis)
    f_decs = bp.fit2_shared(ones, p_base)
    f_decsn = bp.fit2_shared(ones, p_nuis)

    se_lam_n, se_bet_n = _se_lambda(X_base, ds, f_naive["lam"], f_naive["beta"],
                                    f_naive["sigma2"])
    se_lam_u, _ = _se_lambda(X_nuis, ds, f_nuis["lam"], f_nuis["beta"],
                             f_nuis["sigma2"])
    kd_col = bp.Kd_full[:, f_dec.get("td", 0)] if have_decoy else None
    se_lam_d, _ = _se_lambda(X_base, ds, f_dec["lam"], f_dec["beta"],
                             f_dec["sigma2"], extra=kd_col)

    out.update(lam_naive=f_naive["lam"], beta_naive=f_naive["beta"],
               se_lam_naive=se_lam_n, se_beta_naive=se_bet_n,
               lam_nuis=f_nuis["lam"], beta_nuis=f_nuis["beta"],
               se_lam_nuis=se_lam_u,
               lam_decoy=f_dec["lam"], beta_decoy=f_dec["beta"],
               se_lam_decoy=se_lam_d,
               lam_decoy_scale=f_dec["lam_d"], beta_decoy_term=f_dec["beta_d"],
               lam_decoy_nuis=f_decn["lam"], beta_decoy_nuis=f_decn["beta"],
               lam_decoyS=f_decs["lam"], beta_decoyS=f_decs["beta"],
               beta_decoyS_term=f_decs["beta_d"],
               lam_decoyS_nuis=f_decsn["lam"], beta_decoyS_nuis=f_decsn["beta"])

    # ---- binned curves (what the field plots) + the plan's literal scalar --
    bs_idx = bin_index(ds)
    ms, cs, dms = binned_means(bs_idx, y, ds)
    fb = fit_binned_curve(ms, cs, dms)
    out.update(lam_naive_bin=fb["lam"], beta_naive_bin=fb["beta"])
    if have_decoy:
        md, cd, dmd = binned_means(bin_index(dd), y, dd)
        fdo = fit_binned_curve(md, cd, dmd)
    else:
        md = dmd = np.full(N_BINS, np.nan)
        fdo = dict(lam=np.nan, beta=np.nan)
    out.update(lam_decoy_only=fdo["lam"], beta_decoy_only=fdo["beta"])
    out["beta_true_minus_decoy"] = out["beta_naive_bin"] - out["beta_decoy_only"]

    # ---- spatial block bootstrap (cheap, thanks to BlockProfiler) ---------
    lam_n = np.empty(n_boot); bet_n = np.empty(n_boot)
    lam_u = np.empty(n_boot); bet_u = np.empty(n_boot)
    lam_d = np.empty(n_boot); bet_d = np.empty(n_boot)
    lam_s2 = np.empty(n_boot); bet_s2 = np.empty(n_boot)
    for t in range(n_boot):
        m = np.bincount(rng.integers(0, nb, size=nb), minlength=nb).astype(float)
        r1 = bp.fit1(m, p_base); lam_n[t], bet_n[t] = r1["lam"], r1["beta"]
        r2 = bp.fit1(m, p_nuis); lam_u[t], bet_u[t] = r2["lam"], r2["beta"]
        r3 = bp.fit2(m, p_base); lam_d[t], bet_d[t] = r3["lam"], r3["beta"]
        r4 = bp.fit2_shared(m, p_base); lam_s2[t], bet_s2[t] = r4["lam"], r4["beta"]

    ci_n_iid = (f_naive["lam"] - 1.96 * se_lam_n, f_naive["lam"] + 1.96 * se_lam_n)
    ci_u_iid = (f_nuis["lam"] - 1.96 * se_lam_u, f_nuis["lam"] + 1.96 * se_lam_u)
    ci_d_iid = (f_dec["lam"] - 1.96 * se_lam_d, f_dec["lam"] + 1.96 * se_lam_d)
    ci_n_blk, ci_u_blk, ci_d_blk = _ci(lam_n), _ci(lam_u), _ci(lam_d)
    ci_s_blk = _ci(lam_s2)
    cib_n, cib_u, cib_d, cib_s = _ci(bet_n), _ci(bet_u), _ci(bet_d), _ci(bet_s2)

    out.update(
        ci_lam_naive_iid_lo=ci_n_iid[0], ci_lam_naive_iid_hi=ci_n_iid[1],
        ci_lam_nuis_iid_lo=ci_u_iid[0], ci_lam_nuis_iid_hi=ci_u_iid[1],
        ci_lam_decoy_iid_lo=ci_d_iid[0], ci_lam_decoy_iid_hi=ci_d_iid[1],
        ci_lam_naive_blk_lo=ci_n_blk[0], ci_lam_naive_blk_hi=ci_n_blk[1],
        ci_lam_nuis_blk_lo=ci_u_blk[0], ci_lam_nuis_blk_hi=ci_u_blk[1],
        ci_lam_decoy_blk_lo=ci_d_blk[0], ci_lam_decoy_blk_hi=ci_d_blk[1],
        ci_beta_naive_blk_lo=cib_n[0], ci_beta_naive_blk_hi=cib_n[1],
        ci_beta_decoy_blk_lo=cib_d[0], ci_beta_decoy_blk_hi=cib_d[1],
        boot_sd_lam_naive=float(np.nanstd(lam_n)),
        boot_sd_lam_nuis=float(np.nanstd(lam_u)),
        boot_sd_lam_decoy=float(np.nanstd(lam_d)),
        ci_lam_decoyS_blk_lo=ci_s_blk[0], ci_lam_decoyS_blk_hi=ci_s_blk[1],
        ci_beta_decoyS_blk_lo=cib_s[0], ci_beta_decoyS_blk_hi=cib_s[1],
        se_ratio_naive=(float(np.nanstd(lam_n)) / se_lam_n
                        if se_lam_n and np.isfinite(se_lam_n) and se_lam_n > 0
                        else np.nan),
    )

    # ---- coverage (Section 24.7) ------------------------------------------
    out["cover_lam_naive_iid"] = _cover(ci_n_iid, lam_true)
    out["cover_lam_nuis_iid"] = _cover(ci_u_iid, lam_true)
    out["cover_lam_decoy_iid"] = _cover(ci_d_iid, lam_true)
    out["cover_lam_naive_blk"] = _cover(ci_n_blk, lam_true)
    out["cover_lam_nuis_blk"] = _cover(ci_u_blk, lam_true)
    out["cover_lam_decoy_blk"] = _cover(ci_d_blk, lam_true)
    out["cover_beta_naive_blk"] = _cover(cib_n, beta_true)
    out["cover_beta_decoy_blk"] = _cover(cib_d, beta_true)
    out["cover_lam_decoyS_blk"] = _cover(ci_s_blk, lam_true)
    out["cover_beta_decoyS_blk"] = _cover(cib_s, beta_true)
    out["cover_beta_nuis_blk"] = _cover(cib_u, beta_true)

    # ---- bias summaries ----------------------------------------------------
    # beta_true == 0 in the null sweep, so every beta ratio must be guarded
    def _ratio(num, den):
        return num / den if den not in (0, 0.0) else np.nan

    for tag in ["naive", "naive_bin", "nuis", "decoy", "decoy_nuis",
                "decoyS", "decoyS_nuis"]:
        out[f"bias_lam_{tag}"] = out[f"lam_{tag}"] - lam_true
        out[f"relbias_lam_{tag}"] = (out[f"lam_{tag}"] - lam_true) / lam_true
        out[f"relbias_beta_{tag}"] = _ratio(out[f"beta_{tag}"] - beta_true,
                                            beta_true)

    # signal retained (Section 29, objection 6)
    for tag in ["naive", "nuis", "decoy", "decoy_nuis", "decoyS",
                "decoyS_nuis"]:
        out[f"signal_retained_{tag}"] = _ratio(out[f"beta_{tag}"], beta_true)
    out["signal_retained_subtract"] = _ratio(out["beta_true_minus_decoy"],
                                             beta_true)
    out["decoy_cost_vs_naive"] = _ratio(out["beta_decoy"], out["beta_naive"])

    if keep_curves:
        out["_bins"] = BIN_MID
        out["_curve_sender"] = ms
        out["_curve_decoy"] = md
        out["_curve_counts"] = cs
        out["_curve_d"] = dms
        out["_curve_d_decoy"] = dmd
    return out
