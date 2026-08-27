#!/usr/bin/env python3
"""
Phase 3 — the full null battery N1-N8 on real liver sections.

Reports SURVIVING FRACTIONS (Section 6.5), never permutation p-values
(CS Phase 2 §10: all three permutation nulls reject at 42-100 % under a true
null, so their p-values carry no information).

Two families of null, two definitions of the surviving fraction:

  conditioning nulls (N2, N5, N6, zonation-only)
      SF = beta_controlled / beta_naive, both evaluated at the SAME lambda
      (lambda_hat from the naive fit), on the same cells.

  perturbation nulls (N1, N3, N4, N8-scrambled)
      SF = (beta_obs - mean(beta_null)) / beta_obs, lambda held fixed.
      This is the Phase 2 convention and it is the only one that is
      interpretable when the null destroys the kernel's identifiability.

Both carry spatial block bootstrap CIs (10x10 quantile blocks, 400 replicates).

Stage `perm_c1` is `perm` with the Phase 7 / C1 in-tissue geometry
(`phase3_null_geom.Geom`) in place of the whole-section bounding-box wrap:
the nine whole-section variants are fitted on the section's own receivers and
the two tile variants on the solid-tissue tiles only, so each row carries a
`scope` of `full` or `tile`.

Usage
  python3 -u run_phase3_nulls.py --stage window
  python3 -u run_phase3_nulls.py --stage main    --n-jobs 12 --calls all9
  python3 -u run_phase3_nulls.py --stage perm    --n-jobs 12
  python3 -u run_phase3_nulls.py --stage perm_c1 --n-jobs 12
  python3 -u run_phase3_nulls.py --stage perm_c1 --calls tierApm_p95 --tag _pm
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import phase3_core as C
import phase3_null_geom as G

RES = P.RESULTS
os.makedirs(RES, exist_ok=True)

# --- the window -----------------------------------------------------------
# Set from the observed distance-to-nearest-sender distribution, NOT from the
# plan's 300 um (CS Phase 2 §2: 99 % of cells lie within 72-90 um of a sender,
# so 300 um is unreachable and any lambda above ~40 um is extrapolation).
# `--stage window` measures it.  Adopted value: 100 um, the 99th percentile of
# distance-to-nearest-sender pooled over the six Test-3-admissible sections at
# the primary sender call (99.1 % of receivers retained; 96.5 % at cdkn1a_pos).
WINDOW_UM = 100.0
LAM_LO_FLOOR = 7.0            # ~ the median nearest-neighbour distance; the
                              # resolution floor of Section 8 Test 1
N_LAM = 40
N_BLOCKS_SIDE = 10
N_BOOT = 400
MIN_RECEIVERS = 2000

PRIMARY_CALL = "tierA_p95"
N7_CALLS = ["tierA_p90", "tierA_p95", "tierA_p99", "cdkn1a_pos",
            "senepy_p95", "senepy_p99"]
# Phase 8 / D1: the second pre-registered Tier A variant.  The sender score is
# built from the PER-MODULE set `genesets/A_sender_for_<module>.txt`, so one
# call is really seven sender definitions and every stage that uses it has to
# fan out over `P.MODULES` -- that is what `is_permodule` / `_expand` do.
TIERA_PM_CALLS = ["tierApm_p90", "tierApm_p95", "tierApm_p99"]
ALL9_CALLS = N7_CALLS + TIERA_PM_CALLS


def lam_grid(dmax: float = WINDOW_UM, med_nn: float = 0.0) -> np.ndarray:
    """FIXED grid, identical for every section so that lambda-hats are
    comparable.  Lower bound 7 um is at or below the median nearest-neighbour
    distance of every section (6.7-10.6 um), i.e. the resolution floor of
    Section 8 Test 1: anything railing there is not a reportable length scale.
    Upper bound is WINDOW/2, because an exponential with lambda > dmax/2 is not
    distinguishable from a linear trend over [0, dmax]."""
    return np.exp(np.linspace(np.log(LAM_LO_FLOOR), np.log(dmax / 2.0), N_LAM))


# ---------------------------------------------------------------------------
# the per-module sender axis (Phase 8 / D1)
# ---------------------------------------------------------------------------
# A per-module call names seven sender definitions, one per Tier B module, and
# each one may only be used on its own module: fitting `A_sender_for_emt_ecm`
# senders against the interferon response would be a different experiment.  So
# a per-module call fans out over `P.MODULES` at the JOB level (`_expand`) and
# is restricted to a single module at the FIT level (`_js`).
#
# Seeds.  `_expand` reproduces the pre-Phase-8 rule `base + step_i*i +
# step_j*j` exactly for the six non-per-module calls at the same section and
# call ordering, so every number computed before Phase 8 is reproducible from
# it unchanged.  A per-module job is displaced by PM_SEED_OFFSET * (module
# index + 1), so the seven sender sets of one call never share a matching
# draw, a permutation stream or a bootstrap stream, and no per-module seed can
# collide with a non-per-module one.
PM_SEED_OFFSET = 300000


def is_permodule(call: str) -> bool:
    """Is `call` one of the per-module Tier A sensitivity calls?"""
    return str(call).startswith("tierApm_")


def _js(module=None):
    """Module indices a job fits: all seven, or just the sender set's own."""
    if module is None:
        return list(range(len(P.MODULES)))
    if module not in P.MODULES:
        raise ValueError(f"unknown module {module!r}")
    return [P.MODULES.index(module)]


def _expand(sections, calls, base, step_i, step_j, modules=None):
    """(section, call) -> [(section, call, seed, module)], fanning per-module
    calls out over `modules`.  The seed rule for non-per-module calls is
    exactly the pre-Phase-8 `base + step_i*i + step_j*j`."""
    modules = list(P.MODULES) if modules is None else list(modules)
    jobs = []
    for i, s in enumerate(sections):
        for j, c in enumerate(calls):
            sd = int(base) + step_i * i + step_j * j
            if is_permodule(c):
                for mi, m in enumerate(modules):
                    jobs.append((s, c, sd + PM_SEED_OFFSET * (mi + 1), m))
            else:
                jobs.append((s, c, sd, None))
    return jobs


def _resolve_calls(spec):
    """--calls all | all9 | tierA_pm | <comma-separated>"""
    if spec == "all":
        return list(N7_CALLS)
    if spec == "all9":
        return list(ALL9_CALLS)
    if spec == "tierA_pm":
        return list(TIERA_PM_CALLS)
    return spec.split(",")


# ---------------------------------------------------------------------------
# stage 1: the window
# ---------------------------------------------------------------------------

def stage_window(sections):
    """Every sender definition on the N7 axis, on every section.

    This stage does not take `--calls`: the window is a property of the whole
    sender axis, so it always measures all nine calls, per-module calls fanned
    out over the seven Tier B modules (`sender_module`, empty otherwise).
    """
    rows = []
    for s in sections:
        sec = P.Sec(s)
        for call in ALL9_CALLS:
            for module in (P.MODULES if is_permodule(call) else [None]):
                try:
                    snd = sec.sender_mask(call, module=module)
                except Exception as e:
                    print(f"  {s} {call} {module or ''}: {e}")
                    continue
                if snd.sum() < 50:
                    continue
                d = P.dist_to_senders(sec.coords.astype(float), snd)
                recv = (~np.isin(sec.celltype, P.EXCLUDE_TYPES)) & (~snd)
                dr = d[recv]
                rows.append(dict(section=s, arm=sec.meta["condition"],
                                 week=sec.meta["week"], call=call,
                                 sender_module=module,
                                 n=int(recv.sum()), n_senders=int(snd.sum()),
                                 prevalence=float(snd.mean()),
                                 median_nn_um=sec.median_nn_um,
                                 d_p50=float(np.percentile(dr, 50)),
                                 d_p75=float(np.percentile(dr, 75)),
                                 d_p90=float(np.percentile(dr, 90)),
                                 d_p95=float(np.percentile(dr, 95)),
                                 d_p99=float(np.percentile(dr, 99)),
                                 d_max=float(dr.max()),
                                 frac_gt_80=float((dr > 80).mean()),
                                 frac_gt_150=float((dr > 150).mean()),
                                 frac_gt_300=float((dr > 300).mean())))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/window.csv", index=False)
    print(df.to_string(index=False))
    return df


# ---------------------------------------------------------------------------
# per-section setup shared by every fit
# ---------------------------------------------------------------------------

class SectionFit:
    def __init__(self, sample: str, call: str, seed: int, types=None,
                 labels: str = None, module: str = None):
        self.sec = sec = P.Sec(sample, labels=labels)
        self.call = call
        # Phase 8 / D1: `module` names the per-module Tier A sender set this
        # fit uses.  It is None for every other call, and then the sender set
        # is the PRIMARY one.  It is recorded on every row it produces so no
        # output is ambiguous about which Tier A definition made it.
        self.module = module
        self.sender_set = (f"A_sender_for_{module}" if module
                           else "A_SENDER_FINAL_strict")
        self.sender = sec.sender_mask(call, module=module)
        self.coords = sec.coords.astype(float)
        self.d_obs = P.dist_to_senders(self.coords, self.sender)
        self.blocks = C.build_blocks(sec, self.sender, types)
        self.Y = np.column_stack([sec.module(m) for m in P.MODULES])
        self.NB = C.neighbour_baseline(sec, self.Y, self.sender)
        self.bid = P.block_ids(self.coords, N_BLOCKS_SIDE)
        self.nb = N_BLOCKS_SIDE ** 2
        self.lam = lam_grid()
        # N2 matched decoys, WITHIN this section and within cell type
        self.match = C.match_decoys_section(sec, self.sender,
                                            self.blocks["Zmatch"], seed)
        dec = np.zeros(sec.n, bool)
        dec[self.match["decoy_idx"]] = True
        self.decoy = dec
        self.d_dec = (P.dist_to_senders(self.coords, dec) if dec.sum() > 10
                      else np.full(sec.n, np.nan))

    def receivers(self, celltype=None, extra=None):
        m = (~np.isin(self.sec.celltype, P.EXCLUDE_TYPES)) & (~self.sender) \
            & np.isfinite(self.d_obs) & (self.d_obs <= WINDOW_UM)
        if celltype is not None:
            m &= (self.sec.celltype == celltype)
        if extra is not None:
            m &= extra
        return m


def _designs(sf: SectionFit, idx: np.ndarray, j: int):
    """Nested designs.  BlockProfiler._acc(m, p) uses X[:, :p], so one
    profiler yields base -> +N6 -> +N6+zonation -> +N6+N5 for free."""
    one = np.ones((idx.sum(), 1))
    nb = sf.NB[idx, j][:, None]
    zon = sf.blocks["zon"][idx][:, None]
    Z5 = sf.blocks["N5"][idx]
    n5cols = sf.blocks["n5_cols"]
    zi = [n5cols.index("zonation"), n5cols.index("zonation_sq")]
    rest = np.delete(Z5, zi, axis=1)
    X1 = np.column_stack([one, nb, zon, zon ** 2, rest])
    p_base, p_n6, p_n6zon, p_full = 1, 2, 4, X1.shape[1]
    X2 = np.column_stack([one, Z5])
    X3 = np.column_stack([one, zon, zon ** 2])
    return X1, X2, X3, dict(base=p_base, n6=p_n6, n6zon=p_n6zon,
                            n6n5=p_full, n5=X2.shape[1], zon=3)


def fit_cell(sf: SectionFit, celltype: str, j: int, seed: int,
             extra_mask=None, tag: str = "") -> dict:
    """Everything for one (section, sender call, receiver type, module)."""
    idx = sf.receivers(celltype, extra_mask)
    n = int(idx.sum())
    out = dict(section=sf.sec.name, arm=sf.sec.meta["condition"],
               week=sf.sec.meta["week"], call=sf.call, celltype=celltype,
               module=P.MODULES[j], stratum=tag, n=n,
               sender_set=sf.sender_set,
               n_senders=int(sf.sender.sum()),
               prevalence=float(sf.sender.mean()),
               match_rate=sf.match["match_rate"],
               max_smd_before=sf.match["max_smd_before"],
               max_smd_after=sf.match["max_smd_after"])
    if n < MIN_RECEIVERS:
        out["skip"] = "too_few_receivers"
        return out

    y = sf.Y[idx, j].astype(float)
    d_s = sf.d_obs[idx]
    d_d = sf.d_dec[idx]
    bid = sf.bid[idx]
    # re-index blocks so every block is populated
    ub, bid = np.unique(bid, return_inverse=True)
    nb = ub.size
    X1, X2, X3, pp = _designs(sf, idx, j)

    prof1 = E.BlockProfiler(d_s, d_d, y, X1, bid, nb, sf.lam, sf.lam)
    prof2 = E.BlockProfiler(d_s, None, y, X2, bid, nb, sf.lam, sf.lam)
    prof3 = E.BlockProfiler(d_s, d_d, y, X3, bid, nb, sf.lam, sf.lam)
    one = np.ones(nb)

    base = prof1.fit1(one, pp["base"])
    t0 = base["t"]
    out.update(lam_naive=base["lam"], beta_naive=base["beta"],
               lam_grid_lo=float(sf.lam[0]), lam_grid_hi=float(sf.lam[-1]),
               lam_railed=int(t0 == 0 or t0 == sf.lam.size - 1),
               sd_y=float(y.std()))

    # profiled lambda under each design (Section 11 wants the with/without
    # zonation comparison as a panel)
    for nm, (pr, p) in dict(
            n6=(prof1, pp["n6"]), n6zon=(prof1, pp["n6zon"]),
            n6n5=(prof1, pp["n6n5"]), n5=(prof2, pp["n5"]),
            zon=(prof3, pp["zon"])).items():
        r = pr.fit1(one, p)
        out[f"lam_{nm}_profiled"] = r["lam"]
        out[f"beta_{nm}_profiled"] = r["beta"]
    r2 = prof1.fit2_shared(one, pp["base"])
    out["lam_n2_profiled"] = r2["lam"]
    out["beta_n2_profiled"] = r2["beta"]

    # --- fixed-lambda betas: the surviving-fraction numerators -------------
    def _bp(pr, p):
        return pr.beta_at(one, p, t0)[0]

    b = dict(base=out["beta_naive"])
    b["n6"] = _bp(prof1, pp["n6"])
    b["n6zon"] = _bp(prof1, pp["n6zon"])
    b["n6n5"] = _bp(prof1, pp["n6n5"])
    b["n5"] = _bp(prof2, pp["n5"])
    b["zon"] = _bp(prof3, pp["zon"])
    b["n2"] = prof1.beta2_at(one, pp["base"], t0)[0]
    b["n2n5n6"] = prof1.beta2_at(one, pp["n6n5"], t0)[0]
    b["n2zon"] = prof3.beta2_at(one, pp["zon"], t0)[0]
    for k, v in b.items():
        out[f"beta_{k}"] = v
        out[f"sf_{k}"] = v / b["base"] if b["base"] not in (0, np.nan) else np.nan

    # --- spatial block bootstrap ------------------------------------------
    rng = np.random.default_rng(seed)
    keys = list(b)
    boot = {k: np.full(N_BOOT, np.nan) for k in keys}
    for r in range(N_BOOT):
        m = rng.multinomial(nb, np.full(nb, 1.0 / nb)).astype(float)
        try:
            bb = dict(base=prof1.beta_at(m, pp["base"], t0)[0],
                      n6=prof1.beta_at(m, pp["n6"], t0)[0],
                      n6zon=prof1.beta_at(m, pp["n6zon"], t0)[0],
                      n6n5=prof1.beta_at(m, pp["n6n5"], t0)[0],
                      n5=prof2.beta_at(m, pp["n5"], t0)[0],
                      zon=prof3.beta_at(m, pp["zon"], t0)[0],
                      n2=prof1.beta2_at(m, pp["base"], t0)[0],
                      n2n5n6=prof1.beta2_at(m, pp["n6n5"], t0)[0],
                      n2zon=prof3.beta2_at(m, pp["zon"], t0)[0])
        except Exception:
            continue
        for k in keys:
            boot[k][r] = bb[k]
    for k in keys:
        v = boot[k][np.isfinite(boot[k])]
        if v.size > 20:
            out[f"beta_{k}_lo"], out[f"beta_{k}_hi"] = (
                float(np.quantile(v, .025)), float(np.quantile(v, .975)))
    denom = boot["base"]
    for k in keys:
        with np.errstate(divide="ignore", invalid="ignore"):
            sfv = boot[k] / denom
        sfv = sfv[np.isfinite(sfv)]
        if sfv.size > 20:
            out[f"sf_{k}_lo"], out[f"sf_{k}_hi"] = (
                float(np.quantile(sfv, .025)), float(np.quantile(sfv, .975)))
    return out


# ---------------------------------------------------------------------------
# stage 2: the conditioning nulls, all sections x cell types x modules
# ---------------------------------------------------------------------------

def _section_job(sample, call, seed, strata_mode, module=None):
    t0 = time.time()
    sf = SectionFit(sample, call, seed, module=module)
    js = _js(module)
    rows = []
    types = [t for t in sorted(set(sf.sec.celltype))
             if t not in P.EXCLUDE_TYPES]
    for t in types:
        if (sf.receivers(t)).sum() < MIN_RECEIVERS:
            continue
        for j in js:
            rows.append(fit_cell(sf, t, j, seed + 7 * j, tag="all"))
    # pooled over receiver types (cell-type dummies are absorbed by the
    # per-type fits above; the pooled row is reported for the figure only)
    if strata_mode == "zonation":
        comp = sf.sec.compartment
        for zone in ("periportal", "midzonal", "pericentral"):
            mz = comp == zone
            for j in js:
                rows.append(fit_cell(sf, "Hepatocytes", j, seed + 13 * j,
                                     extra_mask=mz, tag=zone))
    print(f"[main] {sample} {call} {module or ''} {len(rows)} rows "
          f"{time.time()-t0:.0f}s", flush=True)
    return rows


def stage_main(sections, calls, n_jobs, strata_mode="zonation"):
    jobs = _expand(sections, calls, P.MASTER_SEED, 1000, 1)
    out = Parallel(n_jobs=n_jobs, prefer="processes", verbose=5)(
        delayed(_section_job)(s, c, sd, strata_mode, m) for s, c, sd, m in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{RES}/main_fits.csv", index=False)
    print(df.shape)
    return df



# ---------------------------------------------------------------------------
# stage 3: the perturbation nulls N1 / N3 / N4
# ---------------------------------------------------------------------------
#
# Vectorised per Section 18.2: one cKDTree query per permutation gives the
# permuted distance-to-nearest-sender for EVERY cell, and that one array then
# serves every (receiver cell type x module) fit.  lambda is HELD FIXED at the
# observed lambda_hat (CS Phase 2 Sec 10) so the null and the observation are
# the same model.

def permute_within_type(rng, sender, celltype, eligible):
    out = np.zeros_like(sender)
    for c in np.unique(celltype):
        ix = np.flatnonzero((celltype == c) & eligible)
        k = int(sender[ix].sum())
        if k:
            out[rng.choice(ix, size=k, replace=False)] = True
    return out


def torus_shift(rng, pts, lo, hi):
    span = hi - lo
    return lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span


def rotate_about_centroid(rng, pts, cen, lo, hi):
    th = rng.uniform(0, 2 * np.pi)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    span = hi - lo
    return lo + ((pts - cen) @ R.T + cen - lo) % span


BINS = np.arange(0.0, WINDOW_UM + 1e-9, 5.0)


def _binned(d, y):
    ix = np.clip(np.digitize(d, BINS) - 1, 0, BINS.size - 2)
    cnt = np.bincount(ix, minlength=BINS.size - 1).astype(float)
    tot = np.bincount(ix, weights=y, minlength=BINS.size - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return tot / np.where(cnt > 0, cnt, np.nan), cnt


def _perm_job(sample, call, seed, n_perm, do_full, curves_for, module=None):
    t0 = time.time()
    sf = SectionFit(sample, call, seed, module=module)
    sec = sf.sec
    rng = np.random.default_rng(seed)
    eligible = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    lo = sf.coords.min(0); hi = sf.coords.max(0); cen = sf.coords.mean(0)
    send_pts = sf.coords[sf.sender]

    types = [t for t in sorted(set(sec.celltype)) if t not in P.EXCLUDE_TYPES]
    cells = []
    for t in types:
        idx = sf.receivers(t)
        if idx.sum() < MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        for j in _js(module):
            y = sf.Y[ii, j].astype(float)
            X1, X2, X3, pp = _designs(sf, idx, j)
            lam, beta, rss, t0i = P.profile_lambda(sf.d_obs[ii],
                                                   X1[:, :1], y, sf.lam)
            rec = dict(t=t, j=j, ii=ii, lam=lam, yc=y - y.mean(),
                       beta_obs=beta, n=ii.size)
            if do_full:
                f = P.FixedLambdaFitter(X1, y[:, None])
                rec["fitter"] = f
                rec["beta_obs_full"] = float(
                    f.beta(np.exp(-sf.d_obs[ii] / lam))[0])
            cells.append(rec)

    nulls = {"N1": np.zeros((len(cells), n_perm)),
             "N3": np.zeros((len(cells), n_perm)),
             "N4": np.zeros((len(cells), n_perm))}
    nulls_full = {k: np.zeros((len(cells), n_perm)) for k in nulls} if do_full \
        else None
    curve_acc = {}
    for nm in ("N1", "N3", "N4"):
        for (t, j) in curves_for:
            curve_acc[(nm, t, j)] = []

    for r in range(n_perm):
        dsets = {}
        m1 = permute_within_type(rng, sf.sender, sec.celltype, eligible)
        dsets["N1"] = P.dist_to_senders(sf.coords, m1)
        dsets["N3"] = P.dist_to_points(sf.coords,
                                       torus_shift(rng, send_pts, lo, hi))
        dsets["N4"] = P.dist_to_points(
            sf.coords, rotate_about_centroid(rng, send_pts, cen, lo, hi))
        for ci, rec in enumerate(cells):
            ii, lam = rec["ii"], rec["lam"]
            for nm, dd in dsets.items():
                k = np.exp(-dd[ii] / lam)
                kc = k - k.mean()
                kk = float(kc @ kc)
                nulls[nm][ci, r] = (kc @ rec["yc"]) / kk if kk > 1e-12 else np.nan
                if do_full:
                    nulls_full[nm][ci, r] = rec["fitter"].beta(k)[0]
                if (nm, rec["t"], rec["j"]) in curve_acc:
                    mu, _ = _binned(dd[ii], sf.Y[ii, rec["j"]].astype(float))
                    curve_acc[(nm, rec["t"], rec["j"])].append(mu)

    rows = []
    for ci, rec in enumerate(cells):
        o = dict(section=sample, arm=sec.meta["condition"],
                 week=sec.meta["week"], call=call, celltype=rec["t"],
                 module=P.MODULES[rec["j"]], n=rec["n"], lam=rec["lam"],
                 n_perm=n_perm, beta_obs=rec["beta_obs"],
                 sender_set=sf.sender_set)
        for nm in ("N1", "N3", "N4"):
            v = nulls[nm][ci]
            v = v[np.isfinite(v)]
            o[f"{nm}_null_mean"] = float(v.mean())
            o[f"{nm}_null_sd"] = float(v.std())
            o[f"{nm}_null_lo"] = float(np.quantile(v, .025))
            o[f"{nm}_null_hi"] = float(np.quantile(v, .975))
            o[f"{nm}_sf"] = float((rec["beta_obs"] - v.mean()) / rec["beta_obs"])
            o[f"{nm}_p"] = float((np.abs(v) >= abs(rec["beta_obs"])).mean())
            if do_full:
                w = nulls_full[nm][ci]; w = w[np.isfinite(w)]
                o[f"{nm}_full_null_mean"] = float(w.mean())
                o[f"{nm}_full_sf"] = float(
                    (rec["beta_obs_full"] - w.mean()) / rec["beta_obs_full"])
        if do_full:
            o["beta_obs_full"] = rec["beta_obs_full"]
        rows.append(o)

    crows = []
    for (nm, t, j), lst in curve_acc.items():
        if not lst:
            continue
        A = np.vstack(lst)
        for bi in range(BINS.size - 1):
            crows.append(dict(section=sample, call=call, celltype=t,
                              module=P.MODULES[j], null=nm,
                              bin_lo=BINS[bi], bin_hi=BINS[bi + 1],
                              mean=float(np.nanmean(A[:, bi])),
                              lo=float(np.nanquantile(A[:, bi], .025)),
                              hi=float(np.nanquantile(A[:, bi], .975))))
    print(f"[perm] {sample} {call} {module or ''} cells={len(cells)} "
          f"perms={n_perm} {time.time()-t0:.0f}s", flush=True)
    return rows, crows


def _perm_tag(calls, tag):
    """Output suffix: explicit `--tag` wins, else the historical rule."""
    return tag if tag else ("" if calls == [PRIMARY_CALL] else "_n7")


def stage_perm(sections, calls, n_jobs, n_perm, do_full=True, tag=""):
    curves_for = [("Hepatocytes", j) for j in range(len(P.MODULES))]
    jobs = _expand(sections, calls, P.MASTER_SEED, 5000, 17)
    out = Parallel(n_jobs=n_jobs, prefer="processes", verbose=5)(
        delayed(_perm_job)(s, c, sd, n_perm, do_full, curves_for, m)
        for s, c, sd, m in jobs)
    df = pd.DataFrame([r for rs, _ in out for r in rs])
    tag = _perm_tag(calls, tag)
    df.to_csv(f"{RES}/perm_nulls{tag}.csv", index=False)
    dc = pd.DataFrame([r for _, cs in out for r in cs])
    dc.to_csv(f"{RES}/perm_curves{tag}.csv", index=False)
    print(df.shape, dc.shape)
    return df


# ---------------------------------------------------------------------------
# stage 3b: the SAME perturbation nulls under the corrected in-tissue geometry
# ---------------------------------------------------------------------------
#
# Phase 7 / correction C1.  `stage_perm` moves the senders over the whole-section
# BOUNDING BOX, and a liver section is not a rectangle, so a large fraction of
# the shifted senders lands in the void outside the tissue (CS_PHASE4 s2.4).
# `phase3_null_geom.Geom` supplies the corrected constructions; this stage runs
# them.  Nothing about the fit changes -- same receivers, same window, same
# lambda grid, same fixed-lambda amplitude -- only where the senders go.
#
# Two scopes, because two of the eleven variants are not defined on the whole
# section:
#   full  the nine whole-section variants, on the section's own receivers and
#         the section's own sender set.  Identical to `stage_perm`'s population,
#         so `beta_obs` here and in `perm_nulls.csv` are the same number.
#   tile  the two solid-tissue-tile variants.  Only senders inside a tile are
#         moved, so only they define the observed distance field, and only
#         receivers inside a tile are fitted.  That is a different (smaller)
#         fit population, and it gets its own rows, its own lambda-hat and its
#         own `beta_obs`.
# One `scope` column separates them; `summarize_phase3_c1.py` reads it.


def _c1_receivers(sf, geom, d, celltype, scope):
    """`SectionFit.receivers` against an arbitrary distance field, optionally
    restricted to the solid-tissue tiles."""
    m = ((~np.isin(sf.sec.celltype, P.EXCLUDE_TYPES)) & (~sf.sender)
         & np.isfinite(d) & (d <= WINDOW_UM))
    if scope == "tile":
        m &= geom.in_tile
    if celltype is not None:
        m &= (sf.sec.celltype == celltype)
    return m


def _perm_c1_job(sample, call, seed, n_perm, do_full, curves_for, module=None):
    t0 = time.time()
    sf = SectionFit(sample, call, seed, module=module)
    sec = sf.sec
    rng = np.random.default_rng(seed)
    eligible = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    geom = G.Geom(sf.coords, sf.sender, eligible)

    tsend = np.zeros(sec.n, bool)
    tsend[geom.tile_sender_idx()] = True
    d_tile = (P.dist_to_senders(sf.coords, tsend) if tsend.any()
              else np.full(sec.n, np.nan))
    dists = {"full": sf.d_obs, "tile": d_tile}
    nulls_of = {"full": G.FULL_NULLS, "tile": G.TILE_NULLS}

    types = [t for t in sorted(set(sec.celltype)) if t not in P.EXCLUDE_TYPES]
    cells = {"full": [], "tile": []}
    for scope in ("full", "tile"):
        d_o = dists[scope]
        for t in types:
            idx = _c1_receivers(sf, geom, d_o, t, scope)
            if idx.sum() < MIN_RECEIVERS:
                continue
            ii = np.flatnonzero(idx)
            for j in _js(module):
                y = sf.Y[ii, j].astype(float)
                X1, X2, X3, pp = _designs(sf, idx, j)
                lam, beta, rss, t0i = P.profile_lambda(d_o[ii], X1[:, :1], y,
                                                       sf.lam)
                rec = dict(t=t, j=j, ii=ii, lam=lam, yc=y - y.mean(),
                           beta_obs=beta, n=ii.size, scope=scope)
                if do_full:
                    f = P.FixedLambdaFitter(X1, y[:, None])
                    rec["fitter"] = f
                    rec["beta_obs_full"] = float(
                        f.beta(np.exp(-d_o[ii] / lam))[0])
                cells[scope].append(rec)

    nulls = {nm: np.zeros((len(cells[sc]), n_perm))
             for sc, nms in nulls_of.items() for nm in nms}
    nulls_full = ({nm: np.zeros((len(cells[sc]), n_perm))
                   for sc, nms in nulls_of.items() for nm in nms}
                  if do_full else None)
    curve_acc = {}
    for nm in G.FULL_NULLS:
        for (t, j) in curves_for:
            curve_acc[(nm, t, j)] = []

    for r in range(n_perm):
        for nm in G.ALL_NULLS:
            scope = "tile" if nm in G.TILE_NULLS else "full"
            if not cells[scope]:
                continue
            pts = geom.draw(nm, rng)
            dd = P.dist_to_points(sf.coords, pts)
            for ci, rec in enumerate(cells[scope]):
                ii, lam = rec["ii"], rec["lam"]
                k = np.exp(-dd[ii] / lam)
                kc = k - k.mean()
                kk = float(kc @ kc)
                nulls[nm][ci, r] = (kc @ rec["yc"]) / kk if kk > 1e-12 else np.nan
                if do_full:
                    nulls_full[nm][ci, r] = rec["fitter"].beta(k)[0]
                if (nm, rec["t"], rec["j"]) in curve_acc:
                    mu, _ = _binned(dd[ii], sf.Y[ii, rec["j"]].astype(float))
                    curve_acc[(nm, rec["t"], rec["j"])].append(mu)

    rows = []
    for scope in ("full", "tile"):
        for ci, rec in enumerate(cells[scope]):
            o = dict(section=sample, arm=sec.meta["condition"],
                     week=sec.meta["week"], call=call, scope=scope,
                     celltype=rec["t"], module=P.MODULES[rec["j"]],
                     n=rec["n"], lam=rec["lam"], n_perm=n_perm,
                     beta_obs=rec["beta_obs"], sender_set=sf.sender_set,
                     n_senders_used=int(tsend.sum() if scope == "tile"
                                        else sf.sender.sum()),
                     n_solid_tiles=len(geom.tiles),
                     tile_cell_coverage=round(geom.tile_cov_cells, 4))
            for nm in nulls_of[scope]:
                v = nulls[nm][ci]
                v = v[np.isfinite(v)]
                o[f"{nm}_null_mean"] = float(v.mean())
                o[f"{nm}_null_sd"] = float(v.std())
                o[f"{nm}_null_lo"] = float(np.quantile(v, .025))
                o[f"{nm}_null_hi"] = float(np.quantile(v, .975))
                o[f"{nm}_sf"] = float(
                    (rec["beta_obs"] - v.mean()) / rec["beta_obs"])
                o[f"{nm}_p"] = float(
                    (np.abs(v) >= abs(rec["beta_obs"])).mean())
                if do_full:
                    w = nulls_full[nm][ci]; w = w[np.isfinite(w)]
                    o[f"{nm}_full_null_mean"] = float(w.mean())
                    o[f"{nm}_full_sf"] = float(
                        (rec["beta_obs_full"] - w.mean()) / rec["beta_obs_full"])
            if do_full:
                o["beta_obs_full"] = rec["beta_obs_full"]
            rows.append(o)

    crows = []
    for (nm, t, j), lst in curve_acc.items():
        if not lst:
            continue
        A = np.vstack(lst)
        for bi in range(BINS.size - 1):
            crows.append(dict(section=sample, call=call, celltype=t,
                              module=P.MODULES[j], null=nm,
                              bin_lo=BINS[bi], bin_hi=BINS[bi + 1],
                              mean=float(np.nanmean(A[:, bi])),
                              lo=float(np.nanquantile(A[:, bi], .025)),
                              hi=float(np.nanquantile(A[:, bi], .975))))
    print(f"[perm_c1] {sample} {call} {module or ''} "
          f"full={len(cells['full'])} tile={len(cells['tile'])} "
          f"perms={n_perm} {time.time()-t0:.0f}s", flush=True)
    return rows, crows


def stage_perm_c1(sections, calls, n_jobs, n_perm, do_full=True, tag=""):
    curves_for = [("Hepatocytes", j) for j in range(len(P.MODULES))]
    jobs = _expand(sections, calls, P.MASTER_SEED, 5000, 17)
    out = Parallel(n_jobs=n_jobs, prefer="processes", verbose=5)(
        delayed(_perm_c1_job)(s, c, sd, n_perm, do_full, curves_for, m)
        for s, c, sd, m in jobs)
    df = pd.DataFrame([r for rs, _ in out for r in rs])
    tag = _perm_tag(calls, tag)
    df.to_csv(f"{RES}/perm_nulls_c1{tag}.csv", index=False)
    dc = pd.DataFrame([r for _, cs in out for r in cs])
    dc.to_csv(f"{RES}/perm_curves_c1{tag}.csv", index=False)
    print(df.shape, dc.shape)
    return df


# ---------------------------------------------------------------------------
# stage 4: observed and matched-decoy binned curves (Figure 2b)
# ---------------------------------------------------------------------------

def stage_curves(sections, call, n_jobs):
    def job(sample):
        sf = SectionFit(sample, call, P.MASTER_SEED)
        rows = []
        for t in sorted(set(sf.sec.celltype)):
            if t in P.EXCLUDE_TYPES:
                continue
            idx = sf.receivers(t)
            if idx.sum() < MIN_RECEIVERS:
                continue
            ii = np.flatnonzero(idx)
            for j, mod in enumerate(P.MODULES):
                y = sf.Y[ii, j].astype(float)
                mu, cnt = _binned(sf.d_obs[ii], y)
                mud, cntd = _binned(sf.d_dec[ii], y)
                sd = np.sqrt(np.array([
                    y[(np.clip(np.digitize(sf.d_obs[ii], BINS) - 1, 0,
                               BINS.size - 2) == b)].var()
                    if cnt[b] > 1 else np.nan for b in range(BINS.size - 1)]))
                for b in range(BINS.size - 1):
                    rows.append(dict(section=sample, arm=sf.sec.meta["condition"],
                                     week=sf.sec.meta["week"], call=call,
                                     celltype=t, module=mod,
                                     bin_lo=BINS[b], bin_hi=BINS[b + 1],
                                     n=cnt[b], mean_obs=mu[b],
                                     sem_obs=sd[b] / np.sqrt(max(cnt[b], 1)),
                                     n_decoy=cntd[b], mean_decoy=mud[b]))
        print(f"[curves] {sample} {len(rows)}", flush=True)
        return rows
    out = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(job)(s) for s in sections)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{RES}/curves.csv", index=False)
    print(df.shape)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["window", "main", "perm", "perm_c1", "curves"])
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--sections", default="sbr")
    ap.add_argument("--calls", default=PRIMARY_CALL,
                    help="all | all9 | tierA_pm | comma-separated call names")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--tag", default="",
                    help="output suffix, e.g. _pm; default is '' for the "
                         "primary call alone and _n7 otherwise")
    a = ap.parse_args()
    secs = {"sbr": P.SBR, "sham": P.SHAM, "all": P.ALL_SECTIONS,
            "inband": P.IN_BAND, "excluded": P.OVER_CEILING + P.BELOW_FLOOR
            }.get(a.sections, a.sections.split(","))
    secs = [s for s in secs if os.path.exists(
        os.path.join(P.CACHE3, f"{s}.npz"))]
    calls = _resolve_calls(a.calls)
    print("sections:", secs, "\ncalls:", calls, flush=True)
    if a.stage == "window":
        stage_window(secs)
    elif a.stage == "main":
        stage_main(secs, calls, a.n_jobs)
    elif a.stage == "perm":
        stage_perm(secs, calls, a.n_jobs, a.n_perm,
                   do_full=(calls == [PRIMARY_CALL]), tag=a.tag)
    elif a.stage == "perm_c1":
        stage_perm_c1(secs, calls, a.n_jobs, a.n_perm,
                      do_full=(calls == [PRIMARY_CALL]), tag=a.tag)
    elif a.stage == "curves":
        stage_curves(secs, calls[0], a.n_jobs)
    else:
        raise SystemExit(f"unknown stage {a.stage}")
