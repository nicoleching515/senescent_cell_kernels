"""
Kernel families and model comparison for the SASP spatial response kernel
(Master Plan Section 6.2), with the uncertainty machinery Phase 1 showed is
mandatory.

Families
  exponential  K(d) = exp(-d/lambda)            degradation-limited diffusion
  gaussian     K(d) = exp(-d^2 / (2 lambda^2))
  powerlaw     K(d) = (1 + d/lambda)^(-p)       heavier tail; 2 nonlinear params
  step         K(d) = 1[d < lambda]             is the response graded at all?
  spline       natural cubic B-spline basis in d, no parametric form

Every family is fitted by the same scheme used in Phase 1: the NONLINEAR
parameters are profiled on a grid, and at each grid point the linear
parameters (intercept, receiver cell-type dummies, any nuisance covariates,
and the kernel amplitude beta) are solved exactly.  This is far more robust
than a generic optimiser, which finds local minima routinely on spatially
confounded data.

Comparability across families: `d_half` (distance at which the kernel falls to
half its value at d=0) and `d_05` (falls to 5%) are reported for every family,
including the spline, so families with different parameterisations can be put
on one axis.

Uncertainty: Phase 1 measured an SE understatement factor of up to 7.9x for the
iid asymptotic CI under spatial confounding, so **both** the iid CI and a
spatial block bootstrap CI are always returned.  The block bootstrap re-uses the
Phase 1 trick: per-block Gram matrices are computed once and a replicate is a
weighted sum of small matrices, never a refit on resampled rows.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# grids
# ---------------------------------------------------------------------------

LAM_LO, LAM_HI = 3.0, 400.0
N_LAM = 64
LAM_GRID = np.exp(np.linspace(np.log(LAM_LO), np.log(LAM_HI), N_LAM))
P_GRID = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])


def basis_exponential(d, lam):
    return np.exp(-d[:, None] / lam[None, :])


def basis_gaussian(d, lam):
    return np.exp(-0.5 * (d[:, None] / lam[None, :]) ** 2)


def basis_step(d, lam):
    return (d[:, None] < lam[None, :]).astype(float)


def basis_powerlaw(d, lam, p):
    return (1.0 + d[:, None] / lam[None, :]) ** (-p)


def spline_basis(d, knots=None, n_knots: int = 6, dmax: float = 300.0):
    """Cubic B-spline design matrix in d (intercept column dropped).

    Returns (B, knots).  The knot vector MUST be carried from the fit to any
    later evaluation -- knots are data-quantile based, so recomputing them on a
    different grid silently produces a different basis and a meaningless curve.
    """
    from scipy.interpolate import BSpline
    if knots is None:
        inner = np.quantile(np.clip(d, 0, dmax),
                            np.linspace(0, 1, n_knots)[1:-1])
        knots = np.r_[[0.0] * 4, inner, [dmax] * 4]
    B = np.asarray(BSpline.design_matrix(
        np.clip(d, 0, dmax), knots, 3, extrapolate=True).todense())
    return B[:, 1:], knots


FAMILIES = ("exponential", "gaussian", "powerlaw", "step", "spline")


# ---------------------------------------------------------------------------
# block-Gram profiler for an arbitrary precomputed basis
# ---------------------------------------------------------------------------


class BasisBlockProfiler:
    """Sufficient statistics for  y = X theta + beta * K[:, j] + eps  over
    spatial blocks, for every column j of a precomputed basis K.

    A bootstrap replicate is an integer multiplicity per block, so every
    cross-product is sum_b m_b * (within-block cross-product).  Per-block Grams
    are computed once; a replicate is then a weighted sum of small matrices.
    """

    def __init__(self, K: np.ndarray, X: np.ndarray, y: np.ndarray,
                 bid: np.ndarray, n_blocks: int):
        n, G = K.shape
        p = X.shape[1]
        self.G, self.p, self.B = G, p, n_blocks
        self.XX = np.zeros((n_blocks, p, p))
        self.XK = np.zeros((n_blocks, p, G))
        self.Xy = np.zeros((n_blocks, p))
        self.KK = np.zeros((n_blocks, G))
        self.Ky = np.zeros((n_blocks, G))
        self.yy = np.zeros(n_blocks)
        self.nn = np.zeros(n_blocks)
        order = np.argsort(bid, kind="stable")
        starts = np.searchsorted(bid[order], np.arange(n_blocks + 1))
        for b in range(n_blocks):
            sl = order[starts[b]:starts[b + 1]]
            if sl.size == 0:
                continue
            Xb, Kb, yb = X[sl], K[sl], y[sl]
            self.XX[b] = Xb.T @ Xb
            self.XK[b] = Xb.T @ Kb
            self.Xy[b] = Xb.T @ yb
            self.KK[b] = np.einsum("ij,ij->j", Kb, Kb)
            self.Ky[b] = Kb.T @ yb
            self.yy[b] = float(yb @ yb)
            self.nn[b] = sl.size

    def profile(self, m: np.ndarray):
        """Return (rss per grid column, beta per column, n)."""
        XX = np.tensordot(m, self.XX, 1)
        XK = np.tensordot(m, self.XK, 1)
        Xy = np.tensordot(m, self.Xy, 1)
        Gi = np.linalg.pinv(XX)
        GX = Gi @ Xy
        yy = float(m @ self.yy) - float(Xy @ GX)
        A = Gi @ XK
        kk = np.tensordot(m, self.KK, 1) - np.einsum("ij,ij->j", XK, A)
        ky = np.tensordot(m, self.Ky, 1) - XK.T @ GX
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(kk > 1e-12, ky / kk, np.nan)
            rss = yy - np.where(kk > 1e-12, ky ** 2 / kk, 0.0)
        return rss, beta, float(m @ self.nn)

    def fit_linear(self, m: np.ndarray):
        """RSS of the covariates-only model (no kernel term)."""
        XX = np.tensordot(m, self.XX, 1)
        Xy = np.tensordot(m, self.Xy, 1)
        Gi = np.linalg.pinv(XX)
        return float(m @ self.yy) - float(Xy @ (Gi @ Xy)), float(m @ self.nn)


def _refine(grid, rss, t):
    lam = float(grid[t])
    if 0 < t < grid.size - 1 and np.all(np.isfinite(rss[t - 1:t + 2])):
        x1, x2, x3 = np.log(grid[t - 1:t + 2])
        y1, y2, y3 = rss[t - 1:t + 2]
        den = (x1 - x2) * (x1 - x3) * (x2 - x3)
        if abs(den) > 0:
            A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / den
            B = (x3 ** 2 * (y1 - y2) + x2 ** 2 * (y3 - y1)
                 + x1 ** 2 * (y2 - y3)) / den
            if A > 0:
                xv = -B / (2 * A)
                if x1 <= xv <= x3:
                    lam = float(np.exp(xv))
    return lam


# ---------------------------------------------------------------------------
# kernel shape summaries, comparable across families
# ---------------------------------------------------------------------------


def kernel_curve(family: str, params: Dict[str, float], d: np.ndarray,
                 spline_coef=None, spline_kw=None) -> np.ndarray:
    if family == "exponential":
        return np.exp(-d / params["lam"])
    if family == "gaussian":
        return np.exp(-0.5 * (d / params["lam"]) ** 2)
    if family == "step":
        return (d < params["lam"]).astype(float)
    if family == "powerlaw":
        return (1.0 + d / params["lam"]) ** (-params["p"])
    if family == "spline":
        B, _ = spline_basis(d, **(spline_kw or {}))
        return B @ spline_coef
    raise ValueError(family)


def shape_summary(family: str, params: Dict[str, float],
                  spline_coef=None, spline_kw=None) -> Dict[str, float]:
    """d_half and d_05: distances at which the kernel falls to 50 % and 5 % of
    its value at d = 0. Comparable across families."""
    dmax = 300.0 if family == "spline" else 600.0
    dd = np.linspace(0, dmax, 3001)
    k = kernel_curve(family, params, dd, spline_coef, spline_kw)
    if not np.all(np.isfinite(k)):
        return dict(d_half=np.nan, d_05=np.nan, amp=np.nan)
    # normalise from the value at d=0 down to the far-field level, so a spline
    # (whose far field is not pinned at 0) is on the same footing as the
    # parametric families (whose far field is 0 by construction).
    # far-field level: median over the last decile of the range, not the single
    # endpoint -- a spline is extrapolating there and its endpoint is noisy.
    k0, kinf = float(k[0]), float(np.median(k[int(0.9 * k.size):]))
    amp = k0 - kinf
    if abs(amp) < 1e-12:
        return dict(d_half=np.nan, d_05=np.nan, amp=amp)
    kn = (k - kinf) / amp
    if amp < 0:               # increasing with distance: no decay length
        return dict(d_half=np.nan, d_05=np.nan, amp=amp)
    out = dict(amp=amp)
    for tag, frac in (("d_half", 0.5), ("d_05", 0.05)):
        below = np.flatnonzero(kn <= frac)
        out[tag] = float(dd[below[0]]) if below.size else np.nan
    return out


# ---------------------------------------------------------------------------
# fitting one family
# ---------------------------------------------------------------------------


def fit_family(d: np.ndarray, y: np.ndarray, X0: np.ndarray, family: str,
               bid: np.ndarray, n_blocks: int, n_boot: int = 300,
               rng: Optional[np.random.Generator] = None,
               lam_grid: np.ndarray = LAM_GRID,
               p_grid: np.ndarray = P_GRID,
               spline_kw: Optional[dict] = None) -> Dict[str, object]:
    """Fit one kernel family; return point estimates, AIC, iid CI and spatial
    block bootstrap CI."""
    rng = rng or np.random.default_rng(0)
    n, p0 = X0.shape
    ones = np.ones(n_blocks)
    res: Dict[str, object] = dict(family=family, n=n)

    if family == "spline":
        B, knots = spline_basis(d, **(spline_kw or {}))
        spline_kw = dict(spline_kw or {}); spline_kw["knots"] = knots
        res["spline_knots"] = knots
        X = np.column_stack([X0, B])
        prof = BasisBlockProfiler(np.zeros((n, 1)), X, y, bid, n_blocks)
        rss0, _ = prof.fit_linear(ones)
        k_par = X.shape[1]
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res.update(params={}, spline_coef=coef[p0:], beta=np.nan,
                   rss=rss0, k_params=k_par)
        res.update(shape_summary("spline", {}, coef[p0:], spline_kw))
        res["aic"] = n * np.log(rss0 / n) + 2 * (k_par + 1)
        prof0 = BasisBlockProfiler(np.zeros((n, 1)), X0, y, bid, n_blocks)
        res["rss0"], _ = prof0.fit_linear(ones)
        res["delta_aic_vs_null"] = (
            res["aic"] - (n * np.log(res["rss0"] / n) + 2 * (p0 + 1)))
        res["se_lam_iid"] = np.nan
        res["ci_lam_iid"] = (np.nan, np.nan)
        res["ci_lam_blk"] = (np.nan, np.nan)
        res["se_ratio"] = np.nan
        # bootstrap the shape summaries
        dh = np.full(n_boot, np.nan)
        for t in range(n_boot):
            m = np.bincount(rng.integers(0, n_blocks, size=n_blocks),
                            minlength=n_blocks).astype(float)
            XX = np.tensordot(m, prof.XX, 1)
            Xy = np.tensordot(m, prof.Xy, 1)
            c = np.linalg.pinv(XX) @ Xy
            dh[t] = shape_summary("spline", {}, c[p0:], spline_kw)["d_half"]
        if np.isfinite(dh).sum() >= 10:
            res["ci_d_half_blk"] = tuple(np.nanquantile(dh, [0.025, 0.975]))
            res["boot_sd_d_half"] = float(np.nanstd(dh))
        else:
            res["ci_d_half_blk"] = (np.nan, np.nan)
            res["boot_sd_d_half"] = np.nan
        return res

    # --- families with nonlinear parameters --------------------------------
    if family == "powerlaw":
        best = None
        for pw in p_grid:
            K = basis_powerlaw(d, lam_grid, pw)
            prof = BasisBlockProfiler(K, X0, y, bid, n_blocks)
            rss, beta, _ = prof.profile(ones)
            t = int(np.nanargmin(rss))
            if best is None or rss[t] < best[0]:
                best = (rss[t], pw, t, prof, float(beta[t]),
                        _refine(lam_grid, rss, t))
            del K
        rss_hat, pw, t, prof, beta_hat, lam_hat = best
        params = dict(lam=lam_hat, p=float(pw))
        k_nl = 2
    else:
        bfun = dict(exponential=basis_exponential, gaussian=basis_gaussian,
                    step=basis_step)[family]
        K = bfun(d, lam_grid)
        prof = BasisBlockProfiler(K, X0, y, bid, n_blocks)
        rss, beta, _ = prof.profile(ones)
        t = int(np.nanargmin(rss))
        rss_hat = float(rss[t])
        beta_hat = float(beta[t])
        lam_hat = float(lam_grid[t]) if family == "step" else _refine(lam_grid, rss, t)
        params = dict(lam=lam_hat)
        k_nl = 1
        del K

    k_par = p0 + 1 + k_nl
    res.update(params=params, beta=beta_hat, rss=rss_hat, k_params=k_par)
    res.update(shape_summary(family, params))

    # --- iid asymptotic CI on lambda (Gauss-Newton) ------------------------
    sigma2 = rss_hat / max(n - k_par, 1)
    eps = max(lam_hat * 1e-3, 1e-3)
    kc = kernel_curve(family, params, d)
    kp = (kernel_curve(family, {**params, "lam": lam_hat + eps}, d) - kc) / eps
    J = np.column_stack([X0, kc, beta_hat * kp])
    try:
        cov = sigma2 * np.linalg.pinv(J.T @ J)
        se_lam = float(np.sqrt(max(cov[-1, -1], 0.0)))
        se_beta = float(np.sqrt(max(cov[p0, p0], 0.0)))
    except np.linalg.LinAlgError:
        se_lam = se_beta = np.nan
    res.update(se_lam_iid=se_lam, se_beta_iid=se_beta,
               ci_lam_iid=(lam_hat - 1.96 * se_lam, lam_hat + 1.96 * se_lam))

    # --- spatial block bootstrap -------------------------------------------
    lam_b = np.full(n_boot, np.nan)
    bet_b = np.full(n_boot, np.nan)
    for tt in range(n_boot):
        m = np.bincount(rng.integers(0, n_blocks, size=n_blocks),
                        minlength=n_blocks).astype(float)
        r, bta, _ = prof.profile(m)
        if np.all(~np.isfinite(r)):
            continue
        j = int(np.nanargmin(r))
        lam_b[tt] = (float(lam_grid[j]) if family == "step"
                     else _refine(lam_grid, r, j))
        bet_b[tt] = bta[j]
    dh_b = np.array([shape_summary(family, {**params, "lam": L})["d_half"]
                     for L in lam_b if np.isfinite(L)])
    res["ci_d_half_blk"] = (tuple(np.nanquantile(dh_b, [0.025, 0.975]))
                            if dh_b.size >= 10 else (np.nan, np.nan))
    res.update(ci_lam_blk=tuple(np.nanquantile(lam_b, [0.025, 0.975])),
               ci_beta_blk=tuple(np.nanquantile(bet_b, [0.025, 0.975])),
               boot_sd_lam=float(np.nanstd(lam_b)),
               se_ratio=(float(np.nanstd(lam_b)) / se_lam
                         if se_lam and np.isfinite(se_lam) and se_lam > 0
                         else np.nan))

    # --- AIC (Gaussian errors) ---------------------------------------------
    res["aic"] = n * np.log(rss_hat / n) + 2 * (k_par + 1)
    res["rss0"], _ = prof.fit_linear(ones)      # covariates-only reference
    res["delta_aic_vs_null"] = (
        res["aic"] - (n * np.log(res["rss0"] / n) + 2 * (p0 + 1)))
    return res


def gaussian_heldout_ll(y_new: np.ndarray, pred: np.ndarray,
                        sigma2_train: float) -> float:
    r = y_new - pred
    return float(-0.5 * y_new.size * np.log(2 * np.pi * sigma2_train)
                 - 0.5 * (r @ r) / sigma2_train)
