"""Phase 4 — the four baseline methods, run identically on real and on
coordinate-shuffled tissue.

  COMMOT        : PUBLISHED SOFTWARE (commot 0.0.3, PyPI).
  CellChat v2   : REIMPLEMENTATION of the published statistic (R package not run).
  SpaTalk       : REIMPLEMENTATION of the published statistic (R package not run).
  NCEM (linear) : REIMPLEMENTATION of the published linear variant (package
                  requires Python<=3.10 + TensorFlow).

Every method returns the same object: a (n_type x n_type) score matrix and a
matching p-value matrix over sender x receiver cell types, for one ligand-
receptor pair on one set of coordinates.  That common interface is what makes
"significant on real, still significant on shifted" a comparable number.

No (n,n) matrix is ever built here; COMMOT builds one internally, which is why
COMMOT is run on tiles (see phase4_run.py).
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree
from scipy import stats

# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

MAXQ = 400       # query cap per cell type in _nn_dist_matrix; the median over
                 # 400 cells is within a few percent of the full-set median and
                 # the cap is applied identically to real and shuffled data


def trimean(v):
    if v.size == 0:
        return 0.0
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    return float((q1 + 2 * q2 + q3) / 4.0)


# CellChat's default group summary is the triMean.  On this panel EVERY one of
# the four ligands is detected in <8 % of cells, so Q25 = Q50 = Q75 = 0 and the
# triMean is IDENTICALLY ZERO in every cell type -- CellChat's default statistic
# is 0 for all four pairs and calls nothing.  That is reported as a finding; the
# null comparison then uses `type = "mean"`, which is a documented CellChat
# option (`computeCommunProb(type = ...)`).
CC_SUMMARY = "mean"


def _group_summary(x, code, K, how=None):
    how = how or CC_SUMMARY
    if how == "trimean":
        return np.array([trimean(x[code == k]) for k in range(K)])
    if how == "mean":
        cnt = np.bincount(code, minlength=K).astype(float)
        return np.bincount(code, weights=x, minlength=K) / np.maximum(cnt, 1.0)
    if how == "truncmean":                      # CellChat truncatedMean, trim 0.1
        out = np.zeros(K)
        for k in range(K):
            v = x[code == k]
            if v.size == 0:
                continue
            lo, hi = np.percentile(v, [10, 90])
            m = (v >= lo) & (v <= hi)
            out[k] = v[m].mean() if m.any() else 0.0
        return out
    raise ValueError(how)


def _group_trimean(x, code, K):
    return _group_summary(x, code, K)


def _nn_dist_matrix(xy, code, K):
    """d[i,j] = median distance from a cell of type i to the NEAREST cell of
    type j.  K tree builds and K*K small queries -- never an (n,n) matrix.
    workers=1 deliberately: the outer loop is already parallel over tiles, and
    thread fan-out on an 11k-point query costs more than it saves."""
    idxs = [np.flatnonzero(code == k) for k in range(K)]
    trees = [cKDTree(xy[i]) if i.size else None for i in idxs]
    Dm = np.full((K, K), np.inf)
    for i in range(K):
        if idxs[i].size == 0:
            continue
        # a median over >=MAXQ cells of type i is already precise to <1%;
        # deterministic stride, not a random draw, so the result is reproducible
        step = max(1, idxs[i].size // MAXQ)
        q = xy[idxs[i][::step]]
        for j in range(K):
            if trees[j] is None:
                continue
            kq = 2 if (i == j and idxs[j].size > 1) else 1
            dd, _ = trees[j].query(q, k=kq, workers=1)
            dd = np.atleast_2d(dd.T).T
            Dm[i, j] = np.median(dd[:, -1])   # self-match dropped when i == j
    return Dm


# ---------------------------------------------------------------------------
# CellChat v2 -- REIMPLEMENTATION of the published statistic
# ---------------------------------------------------------------------------
# Published statistic (Jin et al. 2021 Nat Commun; v2 spatial extension):
#   P_ij = Hill(L_i * R_j)  x  spatial constraint(d_ij)
#   L_i, R_j  = triMean of ligand / receptor expression in cell group i / j
#   Hill(x)   = x / (Kh + x),  Kh = 0.5
#   spatial   = (d0 / d_ij) truncated at 1, zeroed beyond the interaction range
#   p-value   = permutation over CELL GROUP LABELS (CellChat default 100)
# This is NOT the CellChat R package; it is the statistic that package computes.

CC_KH = 0.5
CC_RANGE_UM = 250.0        # CellChat v2 default range for secreted signalling


def cellchat_stat(lig, rec, xy, code, K, d0, D=None, how=None):
    L = _group_summary(lig, code, K, how)
    R = _group_summary(rec, code, K, how)
    if D is None:
        D = _nn_dist_matrix(xy, code, K)
    prod = np.outer(L, R)
    P = prod / (CC_KH + prod)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.minimum(1.0, d0 / D)
    s[~np.isfinite(s)] = 0.0
    s[D > CC_RANGE_UM] = 0.0
    return P * s


def cellchat(lig, rec, xy, code, K, d0, rng, n_perm=100):
    obs = cellchat_stat(lig, rec, xy, code, K, d0)
    ge = np.zeros((K, K))
    n = code.size
    for _ in range(n_perm):
        cp = code[rng.permutation(n)]
        ge += (cellchat_stat(lig, rec, xy, cp, K, d0) >= obs)
    return obs, ge / n_perm


def cellchat_many(pairs, xy, code, K, d0, rng, n_perm=100):
    """All LR pairs share one permutation stream and one geometry recompute."""
    obs = [cellchat_stat(l, r, xy, code, K, d0) for l, r in pairs]
    ge = [np.zeros((K, K)) for _ in pairs]
    n = code.size
    for _ in range(n_perm):
        cp = code[rng.permutation(n)]
        Dp = _nn_dist_matrix(xy, cp, K)
        for t, (l, r) in enumerate(pairs):
            ge[t] += (cellchat_stat(l, r, xy, cp, K, d0, D=Dp) >= obs[t])
    return obs, [g / n_perm for g in ge]


# ---------------------------------------------------------------------------
# SpaTalk -- REIMPLEMENTATION of the published statistic
# ---------------------------------------------------------------------------
# Published statistic (Shao et al. 2022 Nat Commun, `dec_cci`):
#   a spatial neighbour graph over cells (k-nearest neighbours);
#   for sender type A and receiver type B,
#     LRscore = sqrt(Lbar) * sqrt(Rbar) / (sqrt(Lbar)*sqrt(Rbar) + 1)
#   with Lbar the mean ligand expression over A-cells that sit on an A-B edge
#   and Rbar the mean receptor expression over the B-cells they connect to;
#   significance by permutation of cell labels (SpaTalk default 1000).
# This is NOT the SpaTalk R package; it is the statistic that package computes.

SPATALK_K = 10


def spatalk_edges(xy, k=SPATALK_K, max_um=None):
    d, j = cKDTree(xy).query(xy, k=k + 1, workers=1)
    src = np.repeat(np.arange(xy.shape[0]), k)
    dst = j[:, 1:].ravel()
    dd = d[:, 1:].ravel()
    if max_um is not None:
        m = dd <= max_um
        src, dst = src[m], dst[m]
    # undirected: keep both orientations so A->B and B->A are both available
    return np.concatenate([src, dst]), np.concatenate([dst, src])


def _spatalk_stat(lig, rec, src, dst, key, cnt, K):
    ls = np.bincount(key, weights=lig[src], minlength=K * K)
    rs = np.bincount(key, weights=rec[dst], minlength=K * K)
    den = np.maximum(cnt, 1.0)
    g = np.sqrt(ls / den) * np.sqrt(rs / den)
    out = np.where(cnt > 0, g / (g + 1.0), 0.0)
    return out.reshape(K, K)


def spatalk_many(pairs, xy, code, K, rng, n_perm=200, edges=None):
    src, dst = edges if edges is not None else spatalk_edges(xy)

    def prep(c):
        k = c[src] * K + c[dst]
        return k, np.bincount(k, minlength=K * K).astype(float)

    key, cnt = prep(code)
    obs = [_spatalk_stat(l, r, src, dst, key, cnt, K) for l, r in pairs]
    ge = [np.zeros((K, K)) for _ in pairs]
    n = code.size
    for _ in range(n_perm):
        cp = code[rng.permutation(n)]
        kp, cp_cnt = prep(cp)
        for t, (l, r) in enumerate(pairs):
            ge[t] += (_spatalk_stat(l, r, src, dst, kp, cp_cnt, K) >= obs[t])
    return obs, [g / n_perm for g in ge], cnt.reshape(K, K)


# ---------------------------------------------------------------------------
# NCEM linear variant -- REIMPLEMENTATION
# ---------------------------------------------------------------------------
# Published model (Fischer et al. 2023 Nat Biotech, "linear NCEM"):
#   expression of gene g in cell i of type t is linear in the cell-type
#   composition of i's spatial neighbourhood within radius r:
#       y_ig = mu_tg + sum_s beta_tsg * x_is + eps
#   Significant type couplings s->t are the beta with FDR < 0.05.
#   The interaction radius is chosen by sweeping r and maximising variance
#   explained -- that sweep is the "length scale" the method reports.
# The ncem package needs Python<=3.10 and TensorFlow, so this is the linear
# model written out directly; for the linear variant that is exact up to the
# optimiser (OLS closed form vs. keras SGD).

def ncem_linear(y, xy, code, K, radius, min_cells=50):
    """Return (beta, p) matrices [sender s, receiver t] for one response gene."""
    tree = cKDTree(xy)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    n = xy.shape[0]
    # bincount, not np.add.at: the latter is unbuffered and dominates runtime
    flat = np.concatenate([pairs[:, 0] * K + code[pairs[:, 1]],
                           pairs[:, 1] * K + code[pairs[:, 0]]])
    X = np.bincount(flat, minlength=n * K).astype(float).reshape(n, K)
    B = np.full((K, K), np.nan)
    P = np.full((K, K), np.nan)
    r2 = []
    for t in range(K):
        m = code == t
        nt = int(m.sum())
        if nt < min_cells:
            continue
        Xt = np.column_stack([np.ones(nt), X[m]])
        keep = Xt.std(0) > 0
        keep[0] = True
        Xt = Xt[:, keep]
        yt = y[m]
        beta, *_ = np.linalg.lstsq(Xt, yt, rcond=None)
        res = yt - Xt @ beta
        dof = nt - Xt.shape[1]
        if dof <= 1:
            continue
        s2 = float(res @ res) / dof
        try:
            XtX = np.linalg.pinv(Xt.T @ Xt)
        except np.linalg.LinAlgError:
            continue
        se = np.sqrt(np.maximum(np.diag(XtX) * s2, 1e-300))
        tstat = beta / se
        pv = 2 * stats.t.sf(np.abs(tstat), dof)
        cols = np.flatnonzero(keep)[1:] - 1
        B[cols, t] = beta[1:]
        P[cols, t] = pv[1:]
        tot = float(((yt - yt.mean()) ** 2).sum())
        r2.append((nt, 1.0 - float(res @ res) / tot if tot > 0 else np.nan))
    if r2:
        w = np.array([a for a, _ in r2], float)
        v = np.array([b for _, b in r2], float)
        R2 = float(np.nansum(w * v) / w.sum())
    else:
        R2 = np.nan
    return B, P, R2


def bh(p):
    """Benjamini-Hochberg over the finite entries of an array."""
    q = np.full(p.shape, np.nan)
    f = np.isfinite(p)
    v = p[f]
    if v.size == 0:
        return q
    o = np.argsort(v)
    m = v.size
    adj = np.minimum.accumulate((v[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    q[f] = out
    return q
