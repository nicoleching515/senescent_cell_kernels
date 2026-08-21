#!/usr/bin/env python3
"""
Phase 5 T4 — KERNEL FAMILY COMPARISON UNDER FULL CONTROL (Section 6.2), plus
T5 lambda_proximal vs lambda_downstream (Section 6.4).

Families: exponential, gaussian, power law, step/threshold, natural cubic
spline.  Every family is a precomputed basis fed to the same
`sasp_kernels.BasisBlockProfiler`, so the nonlinear parameters are profiled on
a grid and the linear parameters solved exactly, under the SAME design and on
the SAME cells for all five.

Two things Phase 3 did not do and Section 24.6 requires:
  * the comparison is run under the full N5 + N6 control, not naively, so a
    family is not being selected to describe the confounder; and
  * model selection is scored by HELD-OUT log-likelihood on left-out sections
    as well as by AIC.

Phase 2 measured, on synthetic tissue, that the spline wins 0.87-0.93 of the
time when the data are confounded and only 0.07-0.33 of the time when they are
clean -- so "the spline wins" is a confounding marker, and the question here is
whether it survives conditioning.

T5 note on sender definition: the primary call `tierA_p95` is thresholded on
`tierA_score`, which `phase2_downstream.py` computes from
`genesets/A_SENDER_FINAL_strict.txt` -- the union-strict 25-gene Tier A set.
The B1-vs-B4 comparison therefore already uses ONE sender definition for both
modules, which is the condition Section 6.4 requires; no per-module sender set
enters any fit in this file.
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
MIN_N = RN.MIN_RECEIVERS
N_BOOT = 200
FAMS = ["exponential", "gaussian", "powerlaw", "step", "spline"]
P_POW = np.array([0.5, 1.0, 2.0, 4.0])       # power-law exponent grid
N_KNOTS = 6


def family_basis(d: np.ndarray, fam: str, lam: np.ndarray,
                 knots=None) -> tuple:
    """(basis, labels, n_nonlinear_params, knots).

    The spline has no lambda: it contributes its whole basis as extra columns
    at once, so it is fitted as a single 'grid point' whose amplitude count is
    the number of spline columns.
    """
    if fam == "exponential":
        return K.basis_exponential(d, lam), lam, 1, None
    if fam == "gaussian":
        return K.basis_gaussian(d, lam), lam, 1, None
    if fam == "step":
        return K.basis_step(d, lam), lam, 1, None
    if fam == "powerlaw":
        cols, labs = [], []
        for pw in P_POW:
            cols.append(K.basis_powerlaw(d, lam, pw))
            labs += [(l, pw) for l in lam]
        return np.hstack(cols), np.array([l for l, _ in labs]), 2, None
    if fam == "spline":
        B, kn = K.spline_basis(d, knots=knots, n_knots=N_KNOTS,
                               dmax=RN.WINDOW_UM)
        return B, None, 0, kn
    raise ValueError(fam)


def _pow_p(labs_idx, t):
    return float(P_POW[t // labs_idx])


def fit_family_basis(d, y, X, bid, nb, fam, lam, knots=None, rng=None,
                     n_boot=0):
    """Profile one family under design X; AIC and (optionally) a spatial block
    bootstrap of the profiled length scale."""
    B, labs, n_nl, kn = family_basis(d, fam, lam, knots)
    n = float(y.size)
    p = X.shape[1]
    if fam == "spline":
        Xf = np.column_stack([X, B])
        prof = K.BasisBlockProfiler(np.zeros((y.size, 1)), Xf, y, bid, nb)
        one = np.ones(nb)
        rss, _ = prof.fit_linear(one)
        k_par = Xf.shape[1]
        prof0 = K.BasisBlockProfiler(np.zeros((y.size, 1)), X, y, bid, nb)
        rss0, _ = prof0.fit_linear(one)
        out = dict(family=fam, lam=np.nan, p_pow=np.nan, beta=np.nan,
                   rss=float(rss), k_params=k_par,
                   aic=float(n * np.log(rss / n) + 2 * (k_par + 1)),
                   aic0=float(n * np.log(rss0 / n) + 2 * (p + 1)),
                   railed=0, n=n)
        out["basis"] = B
        out["knots"] = kn
        return out
    prof = K.BasisBlockProfiler(B, X, y, bid, nb)
    one = np.ones(nb)
    rss, beta, _ = prof.profile(one)
    t = int(np.nanargmin(rss))
    rss0, _ = prof.fit_linear(one)
    k_par = p + 1 + n_nl
    lam_hat = float(labs[t])
    out = dict(family=fam, lam=lam_hat,
               p_pow=(float(P_POW[t // lam.size]) if fam == "powerlaw"
                      else np.nan),
               beta=float(beta[t]), rss=float(rss[t]), k_params=k_par,
               aic=float(n * np.log(rss[t] / n) + 2 * (k_par + 1)),
               aic0=float(n * np.log(rss0 / n) + 2 * (p + 1)),
               railed=int(lam_hat in (lam[0], lam[-1])), n=n)
    par = dict(lam=lam_hat)
    if fam == "powerlaw":
        par["p"] = out["p_pow"]
    out.update({k: v for k, v in K.shape_summary(fam, par).items()
                if k in ("d_half", "d_05")})
    if n_boot and rng is not None:
        lb = np.full(n_boot, np.nan)
        for r in range(n_boot):
            m = rng.multinomial(nb, np.full(nb, 1.0 / nb)).astype(float)
            rs, _, _ = prof.profile(m)
            if np.all(~np.isfinite(rs)):
                continue
            lb[r] = labs[int(np.nanargmin(rs))]
        v = lb[np.isfinite(lb)]
        if v.size > 20:
            out["lam_blk_lo"] = float(np.quantile(v, .025))
            out["lam_blk_hi"] = float(np.quantile(v, .975))
    out["basis"] = B
    out["t"] = t
    return out


# ---------------------------------------------------------------------------
# within-section: AIC across families, naive vs controlled
# ---------------------------------------------------------------------------

def _section_job(sample, call, seed):
    t0 = time.time()
    sf = RN.SectionFit(sample, call, seed)
    rng = np.random.default_rng(seed)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    for ct in types:
        idx = sf.receivers(ct)
        if idx.sum() < MIN_N:
            continue
        ii = np.flatnonzero(idx)
        d = sf.d_obs[ii]
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        nb = ub.size
        for j, mod in enumerate(P.MODULES):
            y = sf.Y[ii, j].astype(float)
            X1, _, _, pp = RN._designs(sf, idx, j)
            for dname, X in (("naive", X1[:, :1]),
                             ("ctrl", X1[:, :pp["n6n5"]])):
                res = {}
                for fam in FAMS:
                    r = fit_family_basis(d, y, X, bid, nb, fam, sf.lam,
                                         rng=rng, n_boot=N_BOOT)
                    res[fam] = r
                best = min(FAMS, key=lambda f: res[f]["aic"])
                for fam in FAMS:
                    r = res[fam]
                    rows.append(dict(
                        section=sample, arm=sf.sec.meta["condition"],
                        call=call, celltype=ct, module=mod, design=dname,
                        n=int(ii.size), family=fam, lam=r["lam"],
                        p_pow=r["p_pow"], beta=r["beta"], rss=r["rss"],
                        aic=r["aic"], aic0=r["aic0"],
                        d_aic_vs_cov=r["aic"] - r["aic0"],
                        d_half=r.get("d_half", np.nan),
                        d_05=r.get("d_05", np.nan),
                        railed=r["railed"], k_params=r["k_params"],
                        lam_blk_lo=r.get("lam_blk_lo", np.nan),
                        lam_blk_hi=r.get("lam_blk_hi", np.nan),
                        best_family=best,
                        d_aic_vs_best=r["aic"] - res[best]["aic"],
                        sd_y=float(y.std())))
    print(f"[T4 sec] {sample} {len(rows)} rows {time.time()-t0:.0f}s",
          flush=True)
    return rows


# ---------------------------------------------------------------------------
# held-out log-likelihood on left-out sections, per family (Section 24.6)
# ---------------------------------------------------------------------------

def stage_heldout(sections, call, n_jobs):
    packs = [RN.SectionFit(s, call, P.MASTER_SEED + 31 * i,
                           types=P.CANON_TYPES_MERGED)
             for i, s in enumerate(sections)]
    types = sorted(set.intersection(*[set(sf.sec.celltype) for sf in packs])
                   - set(P.EXCLUDE_TYPES))
    lam = packs[0].lam
    # one common spline knot vector across sections, else the basis differs
    dall = np.concatenate([sf.d_obs[np.isfinite(sf.d_obs)] for sf in packs])
    _, KNOTS = K.spline_basis(dall[:200000], n_knots=N_KNOTS,
                              dmax=RN.WINDOW_UM)
    del dall

    def job(ct):
        out = []
        for j, mod in enumerate(P.MODULES):
            parts = []
            for sf in packs:
                idx = sf.receivers(ct)
                if idx.sum() < MIN_N:
                    continue
                ii = np.flatnonzero(idx)
                X1, _, _, pp = RN._designs(sf, idx, j)
                X = X1[:, :pp["n6n5"]].copy()
                X[:, 1:] = C5.zsc_cols(X[:, 1:])
                parts.append(dict(name=sf.sec.name, ii=ii, sf=sf,
                                  y=C5.zsc(sf.Y[ii, j].astype(float)),
                                  X=X, Xnaive=X[:, :1],
                                  d=sf.d_obs[ii]))
            if len(parts) < 3:
                continue
            nb = len(parts)
            bid = np.concatenate([np.full(p["y"].size, a)
                                  for a, p in enumerate(parts)])
            ys = np.concatenate([p["y"] for p in parts])
            for fam in FAMS:
                for p_ in parts:
                    B, labs, n_nl, _ = family_basis(p_["d"], fam, lam, KNOTS)
                    p_["B"] = C5.zsc_cols(B)
                for dname, xk in (("ctrl", "X"), ("naive", "Xnaive")):
                    Xs = np.vstack([p_[xk] for p_ in parts])
                    pp2 = [dict(p_, X=p_[xk]) for p_ in parts]
                    if fam == "spline":
                        # spline enters as covariate columns, not a 1-df term
                        rows = _heldout_block(pp2, xk, "B", nb)
                    else:
                        prof = K.BasisBlockProfiler(
                            np.vstack([p_["B"] for p_ in parts]), Xs, ys,
                            bid, nb)
                        prof0 = K.BasisBlockProfiler(
                            np.zeros((ys.size, 1)), Xs, ys, bid, nb)
                        labs2 = (labs if labs is not None else lam)
                        rows = C5.heldout_ll(prof, pp2, labs2, "B", prof0)
                    for r in rows:
                        out.append(dict(r, celltype=ct, module=mod,
                                        family=fam, design=dname))
        print(f"[T4 heldout] {ct} {len(out)} rows", flush=True)
        return out

    res = Parallel(n_jobs=min(n_jobs, max(len(types), 1)), prefer="threads")(
        delayed(job)(t) for t in types)
    return pd.DataFrame([r for rs in res for r in rs])


def _heldout_block(parts, xk, bkey, nb):
    """Leave-one-section-out for a MULTI-column kernel block (the spline)."""
    rows = []
    for h in range(nb):
        tr = [p for i, p in enumerate(parts) if i != h]
        te = parts[h]
        Atr = np.vstack([np.column_stack([p["X"], p[bkey]]) for p in tr])
        ytr = np.concatenate([p["y"] for p in tr])
        c, *_ = np.linalg.lstsq(Atr, ytr, rcond=None)
        r = ytr - Atr @ c
        s2 = float(r @ r) / max(Atr.shape[0] - Atr.shape[1], 1)
        pred = np.column_stack([te["X"], te[bkey]]) @ c
        ll = K.gaussian_heldout_ll(te["y"], pred, s2)
        X0tr = np.vstack([p["X"] for p in tr])
        c0, *_ = np.linalg.lstsq(X0tr, ytr, rcond=None)
        r0 = ytr - X0tr @ c0
        s20 = float(r0 @ r0) / max(X0tr.shape[0] - X0tr.shape[1], 1)
        ll0 = K.gaussian_heldout_ll(te["y"], te["X"] @ c0, s20)
        rows.append(dict(held_out=te["name"], lam=np.nan, beta=np.nan,
                         n_test=int(te["y"].size), ll=ll,
                         ll_per_cell=ll / te["y"].size, ll0=ll0,
                         ll0_per_cell=ll0 / te["y"].size, dll_vs_cov=ll - ll0,
                         rss_te_per_cell=float(((te["y"] - pred) ** 2).mean())))
    return rows


# ---------------------------------------------------------------------------
# T5 — lambda_proximal (B1) vs lambda_downstream (B4), union-strict senders
# ---------------------------------------------------------------------------

PROX, DOWN = "tnfa_nfkb_proximal", "downstream_arrest"


def stage_proxdown(sections, call, n_jobs):
    """Pooled per receiver cell type, donor (= animal = section) bootstrap on
    the DIFFERENCE log(lam_prox) - log(lam_down), under the full control.

    The two modules are read out against the SAME sender call, so a difference
    cannot be a sender-definition difference (Section 6.4 requirement).
    """
    import sasp_estimators as E
    jp, jd = P.MODULES.index(PROX), P.MODULES.index(DOWN)
    packs = [RN.SectionFit(s, call, P.MASTER_SEED + 31 * i,
                           types=P.CANON_TYPES_MERGED)
             for i, s in enumerate(sections)]
    types = sorted(set.intersection(*[set(sf.sec.celltype) for sf in packs])
                   - set(P.EXCLUDE_TYPES))
    lam = packs[0].lam
    rows = []
    for ct in types:
        prof = {}
        ok = True
        for tag, j in ((PROX, jp), (DOWN, jd)):
            parts = []
            for sf in packs:
                idx = sf.receivers(ct)
                if idx.sum() < MIN_N:
                    continue
                ii = np.flatnonzero(idx)
                X1, _, _, pp = RN._designs(sf, idx, j)
                X = X1[:, :pp["n6n5"]].copy()
                X[:, 1:] = C5.zsc_cols(X[:, 1:])
                parts.append(dict(y=C5.zsc(sf.Y[ii, j].astype(float)), X=X,
                                  d_s=sf.d_obs[ii], d_d=sf.d_dec[ii],
                                  n=ii.size))
            if len(parts) < 3:
                ok = False
                break
            nb = len(parts)
            bid = np.concatenate([np.full(p["n"], a)
                                  for a, p in enumerate(parts)])
            prof[tag] = (E.BlockProfiler(
                np.concatenate([p["d_s"] for p in parts]),
                np.concatenate([p["d_d"] for p in parts]),
                np.concatenate([p["y"] for p in parts]),
                np.vstack([p["X"] for p in parts]), bid, nb, lam, lam),
                np.vstack([p["X"] for p in parts]).shape[1], nb)
        if not ok:
            continue
        nb = prof[PROX][2]
        one = np.ones(nb)
        res = {}
        for tag in (PROX, DOWN):
            pr, pfull, _ = prof[tag]
            naive = pr.fit1(one, 1)
            ctrl = pr.fit1(one, pfull)
            n2 = pr.fit2_shared(one, pfull)
            res[tag] = dict(lam_naive=naive["lam"], beta_naive=naive["beta"],
                            lam_ctrl=ctrl["lam"], beta_ctrl=ctrl["beta"],
                            lam_n2=n2["lam"], beta_n2=n2["beta"],
                            beta_ctrl_fixed=pr.beta_at(one, pfull,
                                                       naive["t"])[0],
                            beta_n2_fixed=pr.beta2_at(one, pfull,
                                                      naive["t"])[0],
                            railed_ctrl=int(ctrl["lam"] <= lam[0] * 1.001
                                            or ctrl["lam"] >= lam[-1] * .999))
        rng = np.random.default_rng(P.MASTER_SEED + 4242)
        dl = []
        lp, ld = [], []
        for _ in range(2000):
            m = rng.multinomial(nb, np.full(nb, 1.0 / nb)).astype(float)
            try:
                a = prof[PROX][0].fit1(m, prof[PROX][1])["lam"]
                b = prof[DOWN][0].fit1(m, prof[DOWN][1])["lam"]
            except Exception:
                continue
            lp.append(a); ld.append(b)
            dl.append(np.log(a) - np.log(b))
        dl = np.asarray(dl); lp = np.asarray(lp); ld = np.asarray(ld)
        r = dict(celltype=ct, n_donors=nb, call=call,
                 lam_prox_ctrl=res[PROX]["lam_ctrl"],
                 lam_down_ctrl=res[DOWN]["lam_ctrl"],
                 log_ratio=float(np.log(res[PROX]["lam_ctrl"]
                                        / res[DOWN]["lam_ctrl"])),
                 log_ratio_lo=float(np.quantile(dl, .025)),
                 log_ratio_hi=float(np.quantile(dl, .975)),
                 ratio=float(np.exp(np.log(res[PROX]["lam_ctrl"]
                                           / res[DOWN]["lam_ctrl"]))),
                 ratio_lo=float(np.exp(np.quantile(dl, .025))),
                 ratio_hi=float(np.exp(np.quantile(dl, .975))),
                 frac_prox_gt_down=float((dl > 0).mean()),
                 lam_prox_lo=float(np.quantile(lp, .025)),
                 lam_prox_hi=float(np.quantile(lp, .975)),
                 lam_down_lo=float(np.quantile(ld, .025)),
                 lam_down_hi=float(np.quantile(ld, .975)),
                 grid_lo=float(lam[0]), grid_hi=float(lam[-1]))
        for tag, pre in ((PROX, "prox"), (DOWN, "down")):
            for k, v in res[tag].items():
                r[f"{pre}_{k}"] = v
            r[f"{pre}_sf_ctrl"] = (res[tag]["beta_ctrl_fixed"]
                                   / res[tag]["beta_naive"])
            r[f"{pre}_sf_n2ctrl"] = (res[tag]["beta_n2_fixed"]
                                     / res[tag]["beta_naive"])
        rows.append(r)
        print(f"[T5] {ct}: lam_prox {r['lam_prox_ctrl']:.1f} "
              f"lam_down {r['lam_down_ctrl']:.1f} "
              f"ratio {r['ratio']:.2f} [{r['ratio_lo']:.2f}, "
              f"{r['ratio_hi']:.2f}]", flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--sections", default="inband")
    ap.add_argument("--call", default=RN.PRIMARY_CALL)
    ap.add_argument("--n-jobs", type=int, default=6)
    a = ap.parse_args()
    secs = {"inband": P.IN_BAND, "all": P.ALL_SECTIONS}.get(
        a.sections, a.sections.split(","))
    if a.stage == "section":
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_section_job)(s, a.call, P.MASTER_SEED + 300 * i)
            for i, s in enumerate(secs))
        df = pd.DataFrame([r for rs in out for r in rs])
        df.to_csv(f"{RES}/kernel_families.csv", index=False)
        print(df.shape)
    elif a.stage == "heldout":
        df = stage_heldout(secs, a.call, a.n_jobs)
        df.to_csv(f"{RES}/kernel_heldout.csv", index=False)
        print(df.shape)
    elif a.stage == "proxdown":
        df = stage_proxdown(secs, a.call, a.n_jobs)
        df.to_csv(f"{RES}/proximal_vs_downstream.csv", index=False)
        print(df.shape)
    else:
        raise SystemExit(a.stage)
