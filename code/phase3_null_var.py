#!/usr/bin/env python3
"""Phase 8 — the VARIANCE-CORRECTED random shift null (N3-var / N4-var).

WHY THIS FILE EXISTS
--------------------
`phase3_null_geom.py` (correction C1) replaced the whole-section bounding-box
torus shift with three in-tissue variants (tile / occ / swap).  An external
literature review then established that this is not new statistics:

  Lotwick HW, Silverman BW (1982), *Methods for analysing spatial processes of
  several types of points*, JRSS-B 44(3):406-413 -- the toroidal shift test;
  the window must be a RECTANGLE.

  Mrkvicka T, Dvorak J, Gonzalez JA, Mateu J (2021), *Revisiting the random
  shift approach for testing in spatial statistics*, Spatial Statistics
  42:100430; arXiv:1911.00240 -- diagnoses the liberality of the torus
  correction and proposes the RANDOM SHIFT WITH VARIANCE CORRECTION, which is
  the accepted remedy for IRREGULAR windows.

This module implements the Mrkvicka et al. (2021) variance correction, §2.1.3
and §2.2 of the paper, verbatim.  Quoting the paper (arXiv:1911.00240v2, via
ar5iv):

  §2.1.3 "For the shift vectors v_1,...,v_N denote W_i = W \\cap (W + v_i),
  i = 1,...,N, the smaller window where both the information about Phi and
  (Psi + v_i) is available.  The first step consists of producing the simulated
  values T_i = T(Phi|_{W_i}, (Psi+v_i)|_{W_i}; W_i), i = 1,...,N. ...
  Thus in the second step, the observed value T_0 and the simulated values
  T_1,...,T_N are standardized to have zero mean and equal variance, i.e. we
  subtract the overall mean Tbar = 1/(N+1) sum_{i=0}^{N} T_i and divide by
  sqrt(var(T_i)):   S_i = (T_i - Tbar)/sqrt(var(T_i)),  i = 0,...,N."

  §2.1.4 "The approaches using minus correction and variance correction can be
  applied in case of general (compact) observation windows."

  §2.2 "for testing the independence assumption for the two random fields Phi
  and Psi, a natural choice of the test statistic is the sample covariance
  cov(Phi(X), Psi(X)). ...  Hence var(T_i) ~= C/n_i where n_i is the number of
  sampling locations in W_i, i = 0,1,...,N, and C is a constant; thus setting
  var(T_i) ~= 1/n_i stabilizes the variance of S_i.  We denote such variance
  correction approach RS_count."

  §2.2 (RS_ker) "var^(T_i) = sum_{k=0}^{N} (T_k - Tbar)^2 * w_ik", with
  Nadaraya-Watson weights w_ik = K(||v_i - v_k||/h) / sum_j K(||v_i - v_j||/h),
  K the Epanechnikov kernel and h a bandwidth.

  §5 (simulation study) and §6 (Barro Colorado data): "shift vectors with
  uniform distribution on a disk centered in the origin and having radius 1/2"
  for W = [0,1]^2, and "random shift vectors distributed uniformly on a disk
  with radius 250 m" for a 1000x500 m^2 window -- i.e. a disk of radius equal
  to HALF THE SHORTER SIDE of the window.  That convention is adopted here.

MAPPING ONTO THIS PROJECT
-------------------------
The Phase 3 estimand is the amplitude beta of exp(-d/lambda) in

    y_c = gamma' X_c + beta * exp(-d(x_c, S)/lambda) + eps_c

over receiver cells c, where S is the sender point set.  Write
Psi(x) = exp(-d(x, S)/lambda).  Then translating the SENDERS by v is EXACTLY
translating the field Psi by v, because

    d(x, S + v) = d(x - v, S)      =>    (Psi + v)(x) = Psi(x - v),

so this is the RANDOM FIELD case of the paper, with

    Phi  = the module score field, observed at the cell locations X
    Psi  = exp(-d(., S)/lambda), computable anywhere S is observed
    T    = the sample covariance cov(Phi(X), Psi(X))   (paper's §2.2 choice)

and therefore the RS_count variance correction var(T_i) ~= 1/n_i applies
without modification, with n_i the number of RECEIVER CELLS in W_i.

W is the tissue mask -- irregular, which is the whole point: §2.1.4 says the
variance correction is defined for general compact windows, whereas the torus
correction is not.  W is built here with EXACTLY the construction
`phase3_null_geom.Geom` uses (25 um occupancy grid, 5x5 binary closing, hole
fill), so N3-var sees the same window as N3-occ.

The cell x is retained iff x in W and x - v in W.  Because every real cell sits
in an occupied grid square, x in W always holds, so retention reduces to
x - v in W.  Psi(x - v) is then evaluated from the FULL observed sender set,
which is the correct reading of "(Psi+v)|_{W_i}": the field value at x in W_i
is Psi(x-v) and x-v lies in W, where Psi is observed.  Senders translated out
of W contribute to no retained cell -- that is the "drop" of the drop-and-
standardize, and it is SOFT: the shift vector is never rejected, only the data
that has been translated out of the window is discarded.  Contrast N3-occ,
which imposes a HARD >= 95 % retention constraint on the OFFSET and therefore
collapses to near-identity offsets (`reports/CS_PHASE7_C1.md` §3).

N4-var
------
Mrkvicka et al. define the construction for TRANSLATIONS only.  N4-var is our
rotation analogue of the same drop-and-standardize principle, stated as an
extension and not as the published method: for a rotation R_theta about the
section centroid, W_i = W \\cap R_theta(W), a cell x is retained iff
x in W and R_theta^{-1} x in W, and Psi_theta(x) = Psi(R_theta^{-1} x).  The
RS_count standardization is unchanged because it depends on the shifted window
only through n_i.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from phase3_null_geom import GRID_UM

RADIUS_FRAC = 0.5      # disk radius = RADIUS_FRAC * (shorter bbox side);
                       # Mrkvicka et al. §5-§6 convention (see module docstring)
MIN_N_RETAINED = 50    # a draw retaining fewer cells than this yields no T_i


class VarGeom:
    """Tissue window W and the two variance-corrected shift constructions."""

    def __init__(self, coords, sender, grid_um=GRID_UM,
                 radius_frac=RADIUS_FRAC):
        self.coords = np.asarray(coords, float)
        self.sender = np.asarray(sender, bool)
        self.grid_um = float(grid_um)
        self.lo = self.coords.min(0)
        self.hi = self.coords.max(0)
        self.span = self.hi - self.lo
        self.cen = self.coords.mean(0)
        self.n = self.coords.shape[0]
        self.k = int(self.sender.sum())
        self.send_pts = self.coords[self.sender]

        # ---- W: exactly phase3_null_geom.Geom's tissue mask ---------------
        g = np.floor((self.coords - self.lo) / self.grid_um).astype(np.int64)
        self.nx, self.ny = int(g[:, 0].max()) + 1, int(g[:, 1].max()) + 1
        occ = np.zeros((self.nx, self.ny), bool)
        occ[g[:, 0], g[:, 1]] = True
        self.occ = occ
        self.tissue = ndi.binary_fill_holes(
            ndi.binary_closing(occ, np.ones((5, 5), bool)))
        self.occ_frac = float(occ.mean())
        self.tissue_frac = float(self.tissue.mean())

        # ---- the shift distribution --------------------------------------
        self.radius = float(radius_frac * self.span.min())

    # ------------------------------------------------------------------ #
    def in_W(self, pts):
        """Boolean: which of `pts` lie in the tissue window W."""
        pts = np.asarray(pts, float)
        ij = np.floor((pts - self.lo) / self.grid_um).astype(np.int64)
        ok = ((ij[:, 0] >= 0) & (ij[:, 0] < self.nx)
              & (ij[:, 1] >= 0) & (ij[:, 1] < self.ny))
        out = np.zeros(pts.shape[0], bool)
        if ok.any():
            out[ok] = self.tissue[ij[ok, 0], ij[ok, 1]]
        return out

    # ---- the two random moves ----------------------------------------- #
    def draw_shift(self, rng):
        """v ~ Uniform(disk of radius `self.radius` centred at the origin),
        the shift distribution of Mrkvicka et al. §5/§6."""
        r = self.radius * np.sqrt(rng.uniform())
        th = rng.uniform(0.0, 2 * np.pi)
        return np.array([r * np.cos(th), r * np.sin(th)])

    def draw_angle(self, rng):
        return float(rng.uniform(0.0, 2 * np.pi))

    @staticmethod
    def _R(th):
        return np.array([[np.cos(th), -np.sin(th)],
                         [np.sin(th), np.cos(th)]])

    # ---- pull-back maps: x |-> the point where (Psi + move) reads Psi --- #
    def pullback_shift(self, pts, v):
        return np.asarray(pts, float) - v

    def pullback_rot(self, pts, th):
        return (np.asarray(pts, float) - self.cen) @ self._R(th) + self.cen

    # ---- forward maps: where the senders actually go ------------------- #
    def push_shift(self, pts, v):
        return np.asarray(pts, float) + v

    def push_rot(self, pts, th):
        return (np.asarray(pts, float) - self.cen) @ self._R(th).T + self.cen

    # ------------------------------------------------------------------ #
    def draw(self, null, rng):
        """One draw of the move parameter: a shift vector for N3-var, an angle
        for N4-var."""
        if null == "N3_var":
            return self.draw_shift(rng)
        if null == "N4_var":
            return self.draw_angle(rng)
        raise ValueError(null)

    def pullback(self, null, pts, param):
        """x |-> the point at which the moved field reads the observed one:
        (Psi + v)(x) = Psi(x - v);  Psi_theta(x) = Psi(R_theta^{-1} x)."""
        if null == "N3_var":
            return self.pullback_shift(pts, param)
        if null == "N4_var":
            return self.pullback_rot(pts, param)
        raise ValueError(null)

    def push(self, null, pts, param):
        """Where the SENDERS actually go under the move (diagnostics only)."""
        if null == "N3_var":
            return self.push_shift(pts, param)
        if null == "N4_var":
            return self.push_rot(pts, param)
        raise ValueError(null)

    def displacement(self, null, param):
        """Median distance a sender is moved by this draw."""
        if null == "N3_var":
            return float(np.hypot(*np.asarray(param, float)))
        moved = self.push_rot(self.send_pts, param)
        return float(np.median(np.hypot(*(moved - self.send_pts).T)))


# ---------------------------------------------------------------------- #
# the Mrkvicka et al. (2021) standardization and Monte Carlo test
# ---------------------------------------------------------------------- #

def rs_count(T, n):
    """RS_count (paper §2.2): S_i = (T_i - Tbar) / sqrt(1/n_i).

    T[0], n[0] are the OBSERVED statistic and its sample size on the full
    window W; T[1:], n[1:] are the N shifted values on W_1,...,W_N.
    Tbar is the overall mean over i = 0,...,N (paper §2.1.3)."""
    T = np.asarray(T, float)
    n = np.asarray(n, float)
    ok = np.isfinite(T) & np.isfinite(n) & (n > 0)
    S = np.full(T.shape, np.nan)
    Tbar = float(T[ok].mean())
    S[ok] = (T[ok] - Tbar) * np.sqrt(n[ok])
    return S, Tbar


def _epanechnikov(u):
    u = np.abs(np.asarray(u, float))
    return np.where(u <= 1.0, 0.75 * (1.0 - u ** 2), 0.0)


def rs_ker(T, V, h):
    """RS_ker (paper §2.2, eqs (1)-(2)): the Nadaraya-Watson estimate of
    var(T_i) with the shift vectors as explanatory variables.

    V[0] must be the origin (the zero shift, v_0 = o, paper §2.1.3)."""
    T = np.asarray(T, float)
    V = np.asarray(V, float)
    ok = np.isfinite(T)
    Tbar = float(T[ok].mean())
    D = np.linalg.norm(V[:, None, :] - V[None, :, :], axis=2)
    Kw = _epanechnikov(D / float(h))
    Kw[:, ~ok] = 0.0
    den = Kw.sum(1)
    var = np.full(T.shape, np.nan)
    good = den > 0
    var[good] = (Kw[good] * np.where(ok, (T - Tbar) ** 2, 0.0)).sum(1) / den[good]
    S = np.full(T.shape, np.nan)
    use = good & ok & (var > 0)
    S[use] = (T[use] - Tbar) / np.sqrt(var[use])
    return S, var


def mc_pvalue(S, two_sided=True):
    """Classical Monte Carlo p-value for S_0 against S_1,...,S_N."""
    S = np.asarray(S, float)
    s0, rest = S[0], S[1:]
    rest = rest[np.isfinite(rest)]
    if rest.size == 0 or not np.isfinite(s0):
        return np.nan
    if two_sided:
        hit = int((np.abs(rest) >= abs(s0)).sum())
    else:
        hit = int((rest >= s0).sum())
    return (1.0 + hit) / (1.0 + rest.size)


def sample_cov(a, b):
    """The paper's s_n = 1/(n-1) sum (a_i - abar)(b_i - bbar)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n = a.size
    if n < 2:
        return np.nan
    return float(((a - a.mean()) * (b - b.mean())).sum() / (n - 1))
