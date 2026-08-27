#!/usr/bin/env python3
"""Phase 8 — does OUR implementation of the variance correction reproduce
Mrkvicka et al. (2021)'s own headline result?

The paper's claim is that the random shift with TORUS correction is LIBERAL for
autocorrelated random fields, and that the random shift with VARIANCE
CORRECTION (RS_count) restores an acceptable significance level.  Their §5
design: W = [0,1]^2, X = 100 uniform sampling points, Phi and Psi independent
stationary centred Gaussian fields with isotropic exponential correlation
c(r) = exp(-r/s), s from 0.001 to 0.5, N = 999 shifts uniform on a disk of
radius 1/2, test statistic the sample covariance, 1000 replications.

This script runs a reduced version of that design (fewer replications and
shifts, so it runs in minutes) and adds a second panel the paper does not have:
the SAME comparison on an IRREGULAR window, where the torus correction is not
even defined and is applied anyway on the bounding box -- which is precisely
what `run_phase3_nulls.py --stage perm` did.

Reproducing the paper's qualitative result is the check that our
`phase3_null_var.rs_count` is the paper's estimator and not something adjacent.

-> results/phase3/var_sim_calibration.csv
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import phase3_null_var as V

RES = "/workspace/results/phase3"
GRID = 256          # field is synthesised on a GRID x GRID lattice over
PAD = 3             # a PAD x PAD-times-larger domain and cropped, so the
                    # realisation is NOT periodic -- otherwise the torus
                    # correction would be valid by construction.


def gaussian_field(rng, s, n=GRID, pad=PAD):
    """Stationary centred Gaussian field with c(r) = exp(-r/s) on [0,1]^2,
    by spectral synthesis on a pad-times-larger torus, then cropped."""
    m = n * pad
    u = np.fft.fftfreq(m, d=1.0 / m)
    dx = np.minimum(np.abs(u), m - np.abs(u)) / n      # lag in [0,1] units
    R = np.exp(-np.hypot(dx[:, None], dx[None, :]) / s)
    S = np.fft.fft2(R).real
    S[S < 0] = 0.0
    z = (rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))
    f = np.fft.ifft2(np.sqrt(S) * z).real * m
    f = f[:n, :n]
    return (f - f.mean()) / f.std()


def _bilinear(F, xy):
    """Field value at continuous coordinates in [0,1)^2."""
    n = F.shape[0]
    p = np.clip(np.asarray(xy, float) * n, 0, n - 1 - 1e-9)
    i0 = np.floor(p).astype(int)
    fr = p - i0
    i1 = np.minimum(i0 + 1, n - 1)
    return (F[i0[:, 0], i0[:, 1]] * (1 - fr[:, 0]) * (1 - fr[:, 1])
            + F[i1[:, 0], i0[:, 1]] * fr[:, 0] * (1 - fr[:, 1])
            + F[i0[:, 0], i1[:, 1]] * (1 - fr[:, 0]) * fr[:, 1]
            + F[i1[:, 0], i1[:, 1]] * fr[:, 0] * fr[:, 1])


def blob_mask(n=GRID):
    """An irregular, non-convex window: three overlapping discs minus a bite.
    Tuned to 0.743 of the bounding box, inside the 0.658-0.858 tissue fraction
    the six real sections show (`null_destructiveness.csv`)."""
    g = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(g, g, indexing="ij")
    m = (((X - .38) ** 2 + (Y - .42) ** 2 < .40 ** 2)
         | ((X - .68) ** 2 + (Y - .62) ** 2 < .36 ** 2)
         | ((X - .38) ** 2 + (Y - .84) ** 2 < .24 ** 2))
    m &= ~((X - .80) ** 2 + (Y - .12) ** 2 < .22 ** 2)
    return m


def in_mask(mask, xy):
    n = mask.shape[0]
    p = np.asarray(xy, float)
    ok = np.all((p >= 0) & (p < 1), axis=1)
    out = np.zeros(p.shape[0], bool)
    if ok.any():
        ij = np.clip((p[ok] * n).astype(int), 0, n - 1)
        out[ok] = mask[ij[:, 0], ij[:, 1]]
    return out


def solid_tiles(mask, k):
    """The sub-squares of a k x k partition of the unit square that lie
    ENTIRELY inside the window -- the analogue of `phase3_null_geom`'s
    >=98%-tissue solid tiles."""
    n = mask.shape[0]
    step = n // k
    out = []
    for a in range(k):
        for b in range(k):
            sub = mask[a * step:(a + 1) * step, b * step:(b + 1) * step]
            if sub.mean() >= 0.98:
                out.append((a, b))
    return out


def tile_pullback(X, v, k):
    """N3-tile's construction: ONE shared offset v, wrapped inside each cell of
    a k x k partition.  Returned as the pull-back of X, so the field is read at
    the wrapped point -- identical to shifting the other component forward.

    This is the "finite union of (aligned) rectangles" case of Mrkvicka et al.
    §2.1.4, which predicts MORE cracks in the autocorrelation structure and
    therefore MORE liberality than a single-rectangle torus."""
    side = 1.0 / k
    org = np.floor(X * k) / k
    return org + ((X - org - v) % side)


def one_rep(rng, s, n_pts, n_shift, mask=None, radius=0.5, tile_ks=(4, 8)):
    F1 = gaussian_field(rng, s)
    F2 = gaussian_field(rng, s)
    if mask is None:
        X = rng.uniform(size=(n_pts, 2))
        tmask = {k: None for k in tile_ks}
    else:
        X = np.empty((0, 2))
        while X.shape[0] < n_pts:
            c = rng.uniform(size=(4 * n_pts, 2))
            c = c[in_mask(mask, c)]
            X = np.vstack([X, c])
        X = X[:n_pts]
        # N3-tile restricts BOTH senders and receivers to solid tiles
        tmask = {}
        for k in tile_ks:
            sol = set(solid_tiles(mask, k))
            ti = np.floor(X * k).astype(int)
            tmask[k] = np.array([(a, b) in sol for a, b in ti])
    phi = _bilinear(F1, X)
    T0 = V.sample_cov(phi, _bilinear(F2, X))

    Tt = np.empty(n_shift)               # torus correction, whole window
    Tk = {k: np.empty(n_shift) for k in tile_ks}   # torus inside k x k tiles
    T0k = {}
    for k in tile_ks:
        sel = slice(None) if tmask[k] is None else tmask[k]
        T0k[k] = (V.sample_cov(phi[sel], _bilinear(F2, X[sel]))
                  if (tmask[k] is None or tmask[k].sum() >= 5) else np.nan)
    Tv = np.full(n_shift + 1, np.nan)    # variance correction
    nv = np.full(n_shift + 1, np.nan)
    Vv = np.zeros((n_shift + 1, 2))
    Tv[0], nv[0] = T0, n_pts
    for i in range(n_shift):
        r = radius * np.sqrt(rng.uniform())
        th = rng.uniform(0, 2 * np.pi)
        v = np.array([r * np.cos(th), r * np.sin(th)])
        Vv[i + 1] = v
        Tt[i] = V.sample_cov(phi, _bilinear(F2, (X - v) % 1.0))
        for k in tile_ks:
            sel = slice(None) if tmask[k] is None else tmask[k]
            if tmask[k] is not None and tmask[k].sum() < 5:
                Tk[k][i] = np.nan
                continue
            bp = tile_pullback(X[sel], v, k)
            Tk[k][i] = V.sample_cov(phi[sel], _bilinear(F2, bp))
        back = X - v
        keep = (np.all((back >= 0) & (back < 1), axis=1) if mask is None
                else in_mask(mask, back))
        ni = int(keep.sum())
        nv[i + 1] = ni
        if ni >= 5:
            Tv[i + 1] = V.sample_cov(phi[keep], _bilinear(F2, back[keep]))
    p_torus = V.mc_pvalue(np.concatenate([[T0], Tt]))
    p_tile = [V.mc_pvalue(np.concatenate([[T0k[k]], Tk[k]])) for k in tile_ks]
    S, _ = V.rs_count(Tv, nv)
    p_count = V.mc_pvalue(S)
    S2, _ = V.rs_ker(Tv, Vv, 0.3 * radius)
    p_ker = V.mc_pvalue(S2)
    p_raw = V.mc_pvalue(Tv)              # drop, but do NOT standardize
    return (p_torus, *p_tile, p_count, p_ker, p_raw)


def _batch(seed, s, n_pts, n_shift, mask, k):
    rng = np.random.default_rng(seed)
    return [one_rep(rng, s, n_pts, n_shift, mask) for _ in range(k)]


def main(scales=(0.02, 0.05, 0.15, 0.3), n_rep=300, n_pts=100, n_shift=199,
         seed=20260827, n_jobs=8):
    from joblib import Parallel, delayed
    rows = []
    for lab, mask in (("rectangle W=[0,1]^2 (paper's design)", None),
                      ("irregular W (three-disc blob, 74% of the box)",
                       blob_mask())):
        for s in scales:
            per = max(1, n_rep // n_jobs)
            out = Parallel(n_jobs=n_jobs, prefer="processes")(
                delayed(_batch)(seed + int(1000 * s) + 7919 * b, s, n_pts,
                                n_shift, mask, per)
                for b in range(n_jobs))
            P = np.array([r for o in out for r in o])
            row = dict(window=lab, corr_scale=s, n_rep=int(P.shape[0]),
                       n_shift=n_shift, n_pts=n_pts)
            names = (["torus"] + [f"torus_tile{k}x{k}" for k in (4, 8)]
                     + ["rs_count", "rs_ker", "drop_no_standardize"])
            for k, nm in enumerate(names):
                p = P[:, k]
                p = p[np.isfinite(p)]
                row[f"{nm}_n"] = int(p.size)
                row[f"{nm}_reject_05"] = float((p <= 0.05).mean())
                row[f"{nm}_reject_10"] = float((p <= 0.10).mean())
            rows.append(row)
            print(f"{lab[:28]:28s} s={s:<5} "
                  f"torus={row['torus_reject_05']:.3f} "
                  f"tile4={row['torus_tile4x4_reject_05']:.3f} "
                  f"tile8={row['torus_tile8x8_reject_05']:.3f} "
                  f"rs_count={row['rs_count_reject_05']:.3f} "
                  f"rs_ker={row['rs_ker_reject_05']:.3f} "
                  f"drop_only={row['drop_no_standardize_reject_05']:.3f}",
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/var_sim_calibration.csv", index=False)
    print(f"\nwrote {RES}/var_sim_calibration.csv")
    return df


if __name__ == "__main__":
    main()
