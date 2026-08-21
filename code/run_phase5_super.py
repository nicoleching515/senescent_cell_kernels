#!/usr/bin/env python3
"""
Phase 5 T1 — SUPERPOSITION vs NEAREST-SENDER on real tissue (Section 6.3).

    nearest        r_i = mu_c + beta * exp(-d_i / lam)            + gamma'z_i
    superposition  r_i = mu_c + beta * sum_j exp(-||x_i-x_j||/lam) + gamma'z_i

Both are profiled over the SAME lambda grid, under the SAME design, on the SAME
cells, and carry the same number of free parameters, so
    dAIC = AIC(superposition) - AIC(nearest) = n * log(RSS_sup / RSS_near)
exactly, and the comparison is not a parameter-count artefact.

The comparison is run under FULL CONTROL (N5 nuisance block + N6 receiver
baseline, and an N2 matched-decoy variant), because model-selecting between two
descriptions of the same confounder is not a result.  It is also run under the
N3 torus shift and the N1 stratified label permutation: if superposition still
wins on shifted senders, the verdict is about tissue geometry, not signalling,
and the report must say so.

Stages
  section   within-section AIC + block-bootstrap win fraction (obs)
  nulls     the same verdict on N3-shifted and N1-permuted senders
  heldout   leave-one-section-out held-out log-likelihood (Section 24.6)
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import sasp_kernels as K
import phase5_common as C5

RES = C5.RES5
N_BOOT = 200
MIN_N = RN.MIN_RECEIVERS


# ---------------------------------------------------------------------------
# residualised shared-lambda two-basis fit (the N2 variant)
# ---------------------------------------------------------------------------

def fit2_shared_basis(Ks, Kd, X, y):
    """RSS profile of  y ~ X g + b_s Ks[:,j] + b_d Kd[:,j]  (shared lambda).

    Residualise once on X (thin QR), then every grid point is a 2x2 solve.
    """
    Q, _ = np.linalg.qr(np.asarray(X, float))
    yt = y - Q @ (Q.T @ y)
    Kst = Ks - Q @ (Q.T @ Ks)
    Kdt = Kd - Q @ (Q.T @ Kd)
    a = np.einsum("ij,ij->j", Kst, Kst)
    c = np.einsum("ij,ij->j", Kdt, Kdt)
    b = np.einsum("ij,ij->j", Kst, Kdt)
    u = Kst.T @ yt
    v = Kdt.T @ yt
    yy = float(yt @ yt)
    det = a * c - b ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        bs = (u * c - b * v) / det
        bd = (a * v - b * u) / det
        rss = yy - (bs * u + bd * v)
    ok = np.isfinite(rss) & (det > 1e-10 * np.abs(a * c))
    rss = np.where(ok, rss, np.inf)
    return rss, bs, bd


# ---------------------------------------------------------------------------
# one section
# ---------------------------------------------------------------------------

class Sup:
    """Observed + decoy superposition bases for one section, built once."""

    def __init__(self, sf: RN.SectionFit):
        self.sf = sf
        lam = sf.lam
        self.lam = lam
        co = sf.coords
        self.Ks_near = C5.nearest_basis(sf.d_obs, lam)
        self.Ks_sup = C5.superposition_basis(co, co[sf.sender], lam)
        if sf.decoy.sum() > 10:
            self.Kd_near = C5.nearest_basis(sf.d_dec, lam)
            self.Kd_sup = C5.superposition_basis(co, co[sf.decoy], lam)
        else:
            self.Kd_near = self.Kd_sup = None


def _designs_for(sf, idx, j):
    X1, X2, X3, pp = RN._designs(sf, idx, j)
    return {"naive": X1[:, :1], "ctrl": X1[:, :pp["n6n5"]]}


def _section_job(sample, call, seed):
    t0 = time.time()
    sf = RN.SectionFit(sample, call, seed)
    S = Sup(sf)
    rng = np.random.default_rng(seed)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    for ct in types:
        idx = sf.receivers(ct)
        if idx.sum() < MIN_N:
            continue
        ii = np.flatnonzero(idx)
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        nb = ub.size
        Kn = S.Ks_near[ii]
        Ku = S.Ks_sup[ii]
        for j, mod in enumerate(P.MODULES):
            y = sf.Y[ii, j].astype(float)
            Xs = _designs_for(sf, idx, j)
            base = dict(section=sample, arm=sf.sec.meta["condition"],
                        week=sf.sec.meta["week"], call=call, celltype=ct,
                        module=mod, n=int(ii.size), n_senders=int(sf.sender.sum()))
            # ONE set of block resamples, shared by both bases (paired)
            BM = rng.multinomial(nb, np.full(nb, 1.0 / nb),
                                 size=N_BOOT).astype(float)
            for dname, X in Xs.items():
                rn = C5.profile_basis(Kn, X, y, bid, nb, S.lam, boot_m=BM)
                ru = C5.profile_basis(Ku, X, y, bid, nb, S.lam, boot_m=BM)
                dw = rn["boot_rss_per_cell"] - ru["boot_rss_per_cell"]
                win = float(np.nanmean(dw > 0))
                with np.errstate(divide="ignore", invalid="ignore"):
                    dab = (ii.size * np.log(ru["boot_rss_per_cell"]
                                            / rn["boot_rss_per_cell"])
                           * 1000.0 / ii.size)
                dab = dab[np.isfinite(dab)]
                r = dict(base, design=dname,
                         lam_near=rn["lam"], beta_near=rn["beta"],
                         rss_near=rn["rss"], aic_near=rn["aic"],
                         railed_near=rn["railed"],
                         lam_sup=ru["lam"], beta_sup=ru["beta"],
                         rss_sup=ru["rss"], aic_sup=ru["aic"],
                         railed_sup=ru["railed"],
                         aic0=rn["aic0"], rss0=rn["rss0"],
                         d_aic=ru["aic"] - rn["aic"],
                         d_aic_near_vs_cov=rn["aic"] - rn["aic0"],
                         d_aic_sup_vs_cov=ru["aic"] - ru["aic0"],
                         boot_win_sup=win,
                         d_aic_per1k=1000.0 * (ru["aic"] - rn["aic"]) / ii.size,
                         d_aic_per1k_lo=(float(np.quantile(dab, .025))
                                         if dab.size > 20 else np.nan),
                         d_aic_per1k_hi=(float(np.quantile(dab, .975))
                                         if dab.size > 20 else np.nan),
                         r_regressors=float(np.corrcoef(
                             Kn[:, rn["t"]], Ku[:, ru["t"]])[0, 1]),
                         sd_y=float(y.std()))
                # amplitude in response-sd units, for the bound
                r["beta_sup_sd"] = ru["beta"] * float(Ku[:, ru["t"]].std()) / (y.std() + 1e-12)
                r["beta_near_sd"] = rn["beta"] / (y.std() + 1e-12)
                rows.append(r)
            # --- N2 matched-decoy variant, controlled design, point est ----
            if S.Kd_near is not None:
                X = Xs["ctrl"]
                rn2, bs_n, _ = fit2_shared_basis(Kn, S.Kd_near[ii], X, y)
                ru2, bs_u, _ = fit2_shared_basis(Ku, S.Kd_sup[ii], X, y)
                tn = int(np.argmin(rn2)); tu = int(np.argmin(ru2))
                n = float(ii.size)
                p = X.shape[1] + 3
                rows.append(dict(base, design="ctrl+N2",
                                 lam_near=float(S.lam[tn]), beta_near=float(bs_n[tn]),
                                 rss_near=float(rn2[tn]),
                                 aic_near=n * np.log(rn2[tn] / n) + 2 * (p + 1),
                                 lam_sup=float(S.lam[tu]), beta_sup=float(bs_u[tu]),
                                 rss_sup=float(ru2[tu]),
                                 aic_sup=n * np.log(ru2[tu] / n) + 2 * (p + 1),
                                 d_aic=n * np.log(ru2[tu] / rn2[tn]),
                                 sd_y=float(y.std())))
    print(f"[T1 sec] {sample} {len(rows)} rows {time.time()-t0:.0f}s", flush=True)
    return rows


# ---------------------------------------------------------------------------
# nulls: is the verdict about signalling or about geometry?
# ---------------------------------------------------------------------------

def _null_job(sample, call, seed, n_draw):
    t0 = time.time()
    sf = RN.SectionFit(sample, call, seed)
    rng = np.random.default_rng(seed + 99)
    co = sf.coords
    lo, hi = co.min(0), co.max(0)
    eligible = ~np.isin(sf.sec.celltype,
                        P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    cells = []
    for ct in types:
        idx = sf.receivers(ct)
        if idx.sum() < MIN_N:
            continue
        ii = np.flatnonzero(idx)
        b = sf.bid[ii]
        ub, b = np.unique(b, return_inverse=True)
        cells.append(dict(ct=ct, ii=ii, bid=b, nb=ub.size,
                          X=[_designs_for(sf, idx, j)["ctrl"]
                             for j in range(len(P.MODULES))]))
    rows = []
    for dr in range(n_draw):
        for nm in ("N3", "N1"):
            if nm == "N3":
                pts = C5.torus_shift(rng, co[sf.sender], lo, hi)
            else:
                m1 = C5.permute_within_type(rng, sf.sender, sf.sec.celltype,
                                            eligible)
                pts = co[m1]
            dnear = P.dist_to_points(co, pts)
            Kn_all = C5.nearest_basis(dnear, sf.lam)
            Ku_all = C5.superposition_basis(co, pts, sf.lam)
            for cc in cells:
                ii = cc["ii"]
                Kn, Ku = Kn_all[ii], Ku_all[ii]
                for j, mod in enumerate(P.MODULES):
                    y = sf.Y[ii, j].astype(float)
                    X = cc["X"][j]
                    rn = C5.profile_basis(Kn, X, y, cc["bid"], cc["nb"], sf.lam)
                    ru = C5.profile_basis(Ku, X, y, cc["bid"], cc["nb"], sf.lam)
                    rows.append(dict(section=sample, call=call, null=nm,
                                     draw=dr, celltype=cc["ct"], module=mod,
                                     n=int(ii.size), lam_near=rn["lam"],
                                     lam_sup=ru["lam"], aic_near=rn["aic"],
                                     aic_sup=ru["aic"],
                                     d_aic=ru["aic"] - rn["aic"],
                                     d_aic_near_vs_cov=rn["aic"] - rn["aic0"],
                                     d_aic_sup_vs_cov=ru["aic"] - ru["aic0"]))
        print(f"  [T1 null] {sample} draw {dr+1}/{n_draw} "
              f"{time.time()-t0:.0f}s", flush=True)
    return rows


# ---------------------------------------------------------------------------
# held-out log-likelihood on left-out sections (Section 24.6)
# ---------------------------------------------------------------------------

def stage_heldout(sections, call, n_jobs):
    """Pool the in-band sections; leave one section out at a time."""
    def prep(sample, si):
        sf = RN.SectionFit(sample, call, P.MASTER_SEED + 31 * si,
                           types=P.CANON_TYPES_MERGED)
        S = Sup(sf)
        return sf, S

    packs = [prep(s, i) for i, s in enumerate(sections)]
    types = sorted(set.intersection(*[set(sf.sec.celltype) for sf, _ in packs])
                   - set(P.EXCLUDE_TYPES))

    def job(ct):
        out = []
        for j, mod in enumerate(P.MODULES):
            parts = []
            for (sf, S) in packs:
                idx = sf.receivers(ct)
                if idx.sum() < MIN_N:
                    continue
                ii = np.flatnonzero(idx)
                X1, _, _, pp = RN._designs(sf, idx, j)
                X = X1[:, :pp["n6n5"]].copy()
                X[:, 1:] = C5.zsc_cols(X[:, 1:])
                parts.append(dict(
                    name=sf.sec.name, y=C5.zsc(sf.Y[ii, j].astype(float)),
                    X=X, Xnaive=X[:, :1], near=C5.zsc_cols(S.Ks_near[ii]),
                    sup=C5.zsc_cols(S.Ks_sup[ii])))
            if len(parts) < 3:
                continue
            bid = np.concatenate([np.full(p["y"].size, a)
                                  for a, p in enumerate(parts)])
            nb = len(parts)
            ys = np.concatenate([p["y"] for p in parts])
            lam = packs[0][0].lam
            for dname, xk in (("ctrl", "X"), ("naive", "Xnaive")):
                Xs = np.vstack([p[xk] for p in parts])
                pp2 = [dict(p, X=p[xk]) for p in parts]
                profs = {k: K.BasisBlockProfiler(
                    np.vstack([p[k] for p in parts]), Xs, ys, bid, nb)
                    for k in ("near", "sup")}
                prof0 = K.BasisBlockProfiler(np.zeros((ys.size, 1)), Xs, ys,
                                             bid, nb)
                for k in ("near", "sup"):
                    for r in C5.heldout_ll(profs[k], pp2, lam, k, prof0):
                        out.append(dict(r, celltype=ct, module=mod, basis=k,
                                        design=dname, n_sections=nb))
        print(f"[T1 heldout] {ct} {len(out)} rows", flush=True)
        return out

    res = Parallel(n_jobs=min(n_jobs, len(types)), prefer="threads")(
        delayed(job)(t) for t in types)
    return pd.DataFrame([r for rs in res for r in rs])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--sections", default="inband")
    ap.add_argument("--call", default=RN.PRIMARY_CALL)
    ap.add_argument("--n-jobs", type=int, default=6)
    ap.add_argument("--n-draw", type=int, default=5)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    secs = {"inband": P.IN_BAND, "all": P.ALL_SECTIONS}.get(
        a.sections, a.sections.split(","))
    if a.stage == "section":
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_section_job)(s, a.call, P.MASTER_SEED + 500 * i)
            for i, s in enumerate(secs))
        df = pd.DataFrame([r for rs in out for r in rs])
        df.to_csv(f"{RES}/super_section{a.tag}.csv", index=False)
        print(df.shape)
    elif a.stage == "nulls":
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_null_job)(s, a.call, P.MASTER_SEED + 700 * i, a.n_draw)
            for i, s in enumerate(secs))
        df = pd.DataFrame([r for rs in out for r in rs])
        df.to_csv(f"{RES}/super_nulls{a.tag}.csv", index=False)
        print(df.shape)
    elif a.stage == "heldout":
        df = stage_heldout(secs, a.call, a.n_jobs)
        df.to_csv(f"{RES}/super_heldout{a.tag}.csv", index=False)
        print(df.shape)
    else:
        raise SystemExit(a.stage)
