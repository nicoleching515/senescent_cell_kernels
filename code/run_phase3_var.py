#!/usr/bin/env python3
"""Phase 8 — run the VARIANCE-CORRECTED random shift null (N3-var / N4-var).

    python3 -u code/run_phase3_var.py --stage perm --n-perm 1000 --n-jobs 6

Nothing here overwrites anything.  It writes NEW files:

    results/phase3/perm_nulls_var.csv          per (section, celltype, module)
    results/phase3/perm_draws_var.csv          per (section, null, draw) geometry
    results/phase3/null_destructiveness_var.csv   the C1 diagnostic, for var

The method is Mrkvicka T, Dvorak J, Gonzalez JA, Mateu J (2021), *Revisiting
the random shift approach for testing in spatial statistics*, Spatial
Statistics 42:100430 (arXiv:1911.00240), §2.1.3 + §2.2 -- see the module
docstring of `phase3_null_var.py` for the verbatim quotes and for how the
Phase 3 estimand maps onto the paper's random-field case.

`run_phase3_nulls.py` is imported, never modified: `SectionFit`, `_designs`,
`_expand`, `WINDOW_UM`, `MIN_RECEIVERS` and `PRIMARY_CALL` all come from it, so
the fit population, the lambda grid, the receiver definition and the covariate
design are identical to `--stage perm_c1` by construction.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.spatial import cKDTree

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import phase3_null_var as V

RES = P.RESULTS
NULLS = ("N3_var", "N4_var")
KER_BW_FRAC = (0.1, 0.2, 0.3)     # RS_ker bandwidths, as a fraction of the
                                  # shift-disk radius (paper used 0.05/0.10/
                                  # 0.15 on a unit-square window, i.e. 0.1-0.3
                                  # of its radius-1/2 shift disk)


def _var_job(sample, call, seed, n_perm, do_full=True):
    t0 = time.time()
    sf = RN.SectionFit(sample, call, seed)
    sec = sf.sec
    rng = np.random.default_rng(seed)
    VG = V.VarGeom(sf.coords, sf.sender)
    stree = cKDTree(VG.send_pts)

    types = [t for t in sorted(set(sec.celltype)) if t not in P.EXCLUDE_TYPES]
    cells = []
    for t in types:
        idx = sf.receivers(t)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        for j in range(len(P.MODULES)):
            y = sf.Y[ii, j].astype(float)
            lam, beta, _rss, _t = P.profile_lambda(
                sf.d_obs[ii], np.ones((ii.size, 1)), y, sf.lam)
            k0 = np.exp(-sf.d_obs[ii] / lam)
            rec = dict(t=t, j=j, ii=ii, lam=lam, y=y, n=ii.size,
                       beta_obs=beta, cov_obs=V.sample_cov(k0, y))
            if do_full:
                X1, _X2, _X3, _pp = RN._designs(sf, idx, j)
                rec["X1"] = X1
                f = P.FixedLambdaFitter(X1, y[:, None])
                rec["beta_obs_full"] = float(f.beta(k0)[0])
            cells.append(rec)

    if not cells:
        return [], []

    q = np.unique(np.concatenate([c["ii"] for c in cells]))
    for c in cells:
        c["pos"] = np.searchsorted(q, c["ii"])
    Xq = sf.coords[q]
    d_obs_q = sf.d_obs[q]
    # cells sharing a receiver type share `ii`, hence share the retained mask
    # and the null distances; only lambda and y differ across the 7 modules.
    groups = {}
    for c in cells:
        groups.setdefault(c["t"], []).append(c)
    groups = [(v[0]["pos"], v) for v in groups.values()]

    # T[0] is the observed statistic on the FULL window W (paper: T_0 uses W,
    # T_1..T_N use W_1..W_N).  Slot 0 of every accumulator is therefore the
    # observed value and V[0] = o, the zero shift.
    acc = {}
    for c in cells:
        for nm in NULLS:
            acc[(nm, id(c))] = dict(
                T=np.full(n_perm + 1, np.nan),      # sample covariance
                B=np.full(n_perm + 1, np.nan),      # beta (cov / var(k))
                T0w=np.full(n_perm + 1, np.nan),    # observed cov on W_i
                B0w=np.full(n_perm + 1, np.nan),    # observed beta on W_i
                BF=np.full(n_perm + 1, np.nan),     # beta, full design
                BF0w=np.full(n_perm + 1, np.nan),
                n=np.full(n_perm + 1, np.nan))
            a = acc[(nm, id(c))]
            a["T"][0] = c["cov_obs"]
            a["B"][0] = c["beta_obs"]
            a["T0w"][0] = c["cov_obs"]
            a["B0w"][0] = c["beta_obs"]
            a["n"][0] = c["n"]
            if do_full:
                a["BF"][0] = c["beta_obs_full"]
                a["BF0w"][0] = c["beta_obs_full"]

    Vshift = {nm: np.zeros((n_perm + 1, 2)) for nm in NULLS}
    drows = []
    for nm in NULLS:
        for r in range(n_perm):
            param = VG.draw(nm, rng)
            back = VG.pullback(nm, Xq, param)
            keep = VG.in_W(back)
            # the shift vector each draw is indexed by, for RS_ker.  For a
            # rotation the natural "shift vector" is the displacement of the
            # window centroid's orbit; we use (median sender displacement, 0)
            # rotated onto the angle, so that ||v_i - v_k|| is a metric on the
            # move parameter.
            if nm == "N3_var":
                Vshift[nm][r + 1] = param
                disp = float(np.hypot(*param))
            else:
                dsp = VG.displacement(nm, param)
                Vshift[nm][r + 1] = (dsp * np.cos(param), dsp * np.sin(param))
                disp = dsp
            dd = np.full(q.size, np.nan)
            if keep.any():
                dd[keep] = stree.query(back[keep], k=1, workers=2)[0]
            drows.append(dict(section=sample, call=call, null=nm, draw=r,
                              displacement_um=disp,
                              frac_cells_retained=float(keep.mean())))
            for gpos, gcells in groups:
                m = keep[gpos]
                ni = int(m.sum())
                for c in gcells:
                    acc[(nm, id(c))]["n"][r + 1] = ni
                if ni < V.MIN_N_RETAINED:
                    continue
                posm = gpos[m]
                d_null = dd[posm]
                d_real = d_obs_q[posm]
                for c in gcells:
                    a = acc[(nm, id(c))]
                    yv = c["y"][m]
                    yc = yv - yv.mean()
                    k = np.exp(-d_null / c["lam"])
                    kc = k - k.mean()
                    ky = float(kc @ yc)
                    kk = float(kc @ kc)
                    k0 = np.exp(-d_real / c["lam"])
                    k0c = k0 - k0.mean()
                    k0y = float(k0c @ yc)
                    k0k = float(k0c @ k0c)
                    a["T"][r + 1] = ky / (ni - 1)
                    a["T0w"][r + 1] = k0y / (ni - 1)
                    a["B"][r + 1] = ky / kk if kk > 1e-12 else np.nan
                    a["B0w"][r + 1] = k0y / k0k if k0k > 1e-12 else np.nan
                    if do_full:
                        f = P.FixedLambdaFitter(c["X1"][m], yv[:, None])
                        a["BF"][r + 1] = float(f.beta(k)[0])
                        a["BF0w"][r + 1] = float(f.beta(k0)[0])

    rows = []
    for c in cells:
        base = dict(section=sample, arm=sec.meta["condition"],
                    week=sec.meta["week"], call=call, scope="full",
                    celltype=c["t"], module=P.MODULES[c["j"]], n=c["n"],
                    lam=c["lam"], n_perm=n_perm, beta_obs=c["beta_obs"],
                    cov_obs=c["cov_obs"],
                    sender_set="A_SENDER_FINAL_strict",
                    n_senders_used=int(sf.sender.sum()),
                    shift_radius_um=VG.radius,
                    tissue_frac_bbox=round(VG.tissue_frac, 4))
        if do_full:
            base["beta_obs_full"] = c["beta_obs_full"]
        o = dict(base)
        for nm in NULLS:
            a = acc[(nm, id(c))]
            nT, nB = a["T"][1:], a["B"][1:]
            ok = np.isfinite(nB)
            o[f"{nm}_n_valid"] = int(ok.sum())
            o[f"{nm}_frac_cells_retained"] = float(
                np.nanmean(a["n"][1:]) / c["n"])
            o[f"{nm}_null_mean"] = float(np.nanmean(nB))
            o[f"{nm}_null_sd"] = float(np.nanstd(nB))
            o[f"{nm}_null_lo"] = float(np.nanquantile(nB, .025))
            o[f"{nm}_null_hi"] = float(np.nanquantile(nB, .975))
            # SF exactly as every other variant defines it: how much of the
            # observed amplitude the null does NOT reproduce.
            o[f"{nm}_sf"] = float(
                (c["beta_obs"] - np.nanmean(nB)) / c["beta_obs"])
            # window-matched SF: the observed statistic RE-FITTED on the same
            # W_i, so the retained-window effect cancels.
            o[f"{nm}_sf_wm"] = float(np.nanmean(
                (a["B0w"][1:] - nB) / a["B0w"][1:]))
            o[f"{nm}_beta_obs_wm"] = float(np.nanmean(a["B0w"][1:]))
            # the repo's naive permutation p, for comparability only
            o[f"{nm}_p_naive"] = float(
                np.nanmean(np.abs(nB) >= abs(c["beta_obs"])))
            # --- the paper's test: RS_count standardization ---------------
            for tag, arr in (("cov", a["T"]), ("beta", a["B"])):
                S, _Tbar = V.rs_count(arr, a["n"])
                o[f"{nm}_{tag}_S0"] = float(S[0])
                o[f"{nm}_{tag}_p_rscount"] = V.mc_pvalue(S, two_sided=True)
                o[f"{nm}_{tag}_p_rscount_1s"] = V.mc_pvalue(S, two_sided=False)
            # --- RS_ker, the kernel variance estimate ---------------------
            hmax = float(np.linalg.norm(Vshift[nm], axis=1).max())
            for bf in KER_BW_FRAC:
                S, _var = V.rs_ker(a["T"], Vshift[nm], max(bf * hmax, 1e-9))
                o[f"{nm}_cov_p_rsker{int(bf*100):02d}"] = V.mc_pvalue(S)
            # --- uncorrected: the same drawn statistics with NO variance
            #     standardization (the "minus-correction-less" comparison) ---
            o[f"{nm}_cov_p_raw"] = V.mc_pvalue(
                np.concatenate([[a["T"][0]], nT]))
            if do_full:
                nF = a["BF"][1:]
                o[f"{nm}_full_null_mean"] = float(np.nanmean(nF))
                o[f"{nm}_full_sf"] = float(
                    (c["beta_obs_full"] - np.nanmean(nF)) / c["beta_obs_full"])
                o[f"{nm}_full_sf_wm"] = float(np.nanmean(
                    (a["BF0w"][1:] - nF) / a["BF0w"][1:]))
        rows.append(o)
    print(f"[perm_var] {sample} {call} cells={len(cells)} perms={n_perm} "
          f"{time.time()-t0:.0f}s", flush=True)
    return rows, drows


def stage_perm(sections, calls, n_jobs, n_perm, do_full=True, tag=""):
    jobs = RN._expand(sections, calls, P.MASTER_SEED, 5000, 17)
    out = Parallel(n_jobs=n_jobs, prefer="processes", verbose=5)(
        delayed(_var_job)(s, c, sd, n_perm, do_full) for s, c, sd, _m in jobs)
    df = pd.DataFrame([r for rs, _ in out for r in rs])
    df.to_csv(f"{RES}/perm_nulls_var{tag}.csv", index=False)
    dd = pd.DataFrame([r for _, ds in out for r in ds])
    dd.to_csv(f"{RES}/perm_draws_var{tag}.csv", index=False)
    print(df.shape, dd.shape)
    return df


# ---------------------------------------------------------------------- #
# the C1 destructiveness diagnostic, for the variance-corrected variants
# ---------------------------------------------------------------------- #

def _diag_job(sample, call, n_rep, seed):
    import phase3_null_geom as G
    sec = P.Sec(sample)
    xy = sec.coords.astype(float)
    snd = sec.sender_mask(call)
    elig = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    VG = V.VarGeom(xy, snd)
    geom = G.Geom(xy, snd, elig)      # only for occ-grid comparability
    tree = geom.tree
    import phase3_null_diag as _ND
    rng = np.random.default_rng(seed + _ND.section_offset(sample))
    src = xy[snd]
    base = np.asarray(tree.query_ball_point(src, RN.WINDOW_UM, workers=-1,
                                            return_length=True), float) - 1.0
    rows = []
    for nm in NULLS:
        keep_all, keep_ret, med, occf, disp = [], [], [], [], []
        fsr, fcr = [], []
        for _ in range(n_rep):
            param = VG.draw(nm, rng)
            pushed = VG.push(nm, src, param)
            inW = VG.in_W(pushed)
            cnt = np.asarray(tree.query_ball_point(
                pushed, RN.WINDOW_UM, workers=-1, return_length=True), float)
            keep_all.append(float((cnt > 0).mean()))
            keep_ret.append(float((cnt[inW] > 0).mean()) if inW.any() else np.nan)
            med.append(float(np.median(cnt[inW])) if inW.any() else np.nan)
            occf.append(geom._in_occ(pushed))
            disp.append(VG.displacement(nm, param))
            fsr.append(float(inW.mean()))
            back = VG.pullback(nm, xy, param)
            fcr.append(float(VG.in_W(back).mean()))
        rows.append(dict(
            section=sample, arm=sec.meta["condition"], week=sec.meta["week"],
            call=call, null=nm, scope="whole section (tissue window W)",
            n_senders=int(snd.sum()), n_cells=int(sec.n),
            occ_frac_bbox=round(VG.occ_frac, 4),
            tissue_frac_bbox=round(VG.tissue_frac, 4),
            shift_radius_um=round(VG.radius, 1),
            real_median_nbrs=float(np.median(base)),
            null_median_nbrs=float(np.nanmean(med)),
            frac_retaining_a_neighbour=float(np.nanmean(keep_ret)),
            frac_retaining_a_neighbour_all_shifted=float(np.mean(keep_all)),
            frac_senders_retained=float(np.mean(fsr)),
            frac_cells_retained=float(np.mean(fcr)),
            frac_in_occupancy=float(np.mean(occf)),
            median_displacement_um=float(np.mean(disp))))
        r = rows[-1]
        print(f"  {sample[:4]} {nm:8s} keep={r['frac_retaining_a_neighbour']:.3f} "
              f"nbrs {r['real_median_nbrs']:.1f}->{r['null_median_nbrs']:.1f} "
              f"retW={r['frac_senders_retained']:.3f} "
              f"cellsW={r['frac_cells_retained']:.3f} "
              f"disp={r['median_displacement_um']:.0f}um", flush=True)
    return rows


def stage_diag(sections, call, n_rep, n_jobs, seed=20260827):
    out = Parallel(n_jobs=n_jobs, prefer="processes", verbose=5)(
        delayed(_diag_job)(s, call, n_rep, seed) for s in sections)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{RES}/null_destructiveness_var.csv", index=False)
    print(df.groupby("null", sort=False)[
        ["frac_retaining_a_neighbour", "real_median_nbrs", "null_median_nbrs",
         "frac_senders_retained", "frac_cells_retained",
         "median_displacement_um"]].median().round(3).to_string())
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["perm", "diag"])
    ap.add_argument("--sections", default="inband")
    ap.add_argument("--calls", default=RN.PRIMARY_CALL)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--n-rep", type=int, default=20)
    ap.add_argument("--n-jobs", type=int, default=6)
    ap.add_argument("--tag", default="")
    ap.add_argument("--no-full", action="store_true")
    a = ap.parse_args()
    sections = (P.IN_BAND if a.sections == "inband"
                else a.sections.split(","))
    calls = a.calls.split(",")
    os.makedirs(RES, exist_ok=True)
    if a.stage == "perm":
        stage_perm(sections, calls, a.n_jobs, a.n_perm,
                   do_full=not a.no_full, tag=a.tag)
    else:
        stage_diag(sections, calls[0], a.n_rep, a.n_jobs)


if __name__ == "__main__":
    main()
