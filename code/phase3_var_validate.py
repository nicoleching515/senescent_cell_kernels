#!/usr/bin/env python3
"""Phase 8 — does RS_count's var(T_i) ~= 1/n_i actually hold on OUR data?

Mrkvicka et al. (2021) §2.2 justify the RS_count correction factor 1/n_i with
Theorem 1: for two independent stationary random fields with compactly
supported autocovariances, sampled at points with bounded local multiplicity,
var(s_n) = Theta(1/n).  Our Phi is a module score at Xenium cell centroids and
our Psi is exp(-d(.,S)/lambda) -- neither is Gaussian and neither is stationary
in the strict sense, so the assumption is worth testing rather than assuming.

This script re-draws the N3-var / N4-var moves for a subset of sections,
stores the PER-DRAW statistic T_i and retained count n_i, and reports

    slope of log sd(T | n-bin) against log n  -- should be -0.5 under 1/n
    Spearman rho(n_i, |T_i - Tbar|)           -- should be negative

-> results/phase3/var_variance_check.csv

Nothing is overwritten.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import phase3_null_var as V
import summarize_phase3_c1 as SC1

RES = P.RESULTS
NULLS = ("N3_var", "N4_var")
N_BIN = 8


def _job(sample, call, n_draw, seed, wanted):
    sf = RN.SectionFit(sample, call, seed)
    sec = sf.sec
    rng = np.random.default_rng(seed + 991)
    VG = V.VarGeom(sf.coords, sf.sender)
    stree = cKDTree(VG.send_pts)
    cells = []
    for t in sorted(set(sec.celltype)):
        if t in P.EXCLUDE_TYPES:
            continue
        idx = sf.receivers(t)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        for j, mod in enumerate(P.MODULES):
            if (sample, t, mod) not in wanted:
                continue
            y = sf.Y[ii, j].astype(float)
            lam, beta, _r, _t = P.profile_lambda(
                sf.d_obs[ii], np.ones((ii.size, 1)), y, sf.lam)
            cells.append(dict(t=t, mod=mod, ii=ii, lam=lam, y=y, n=ii.size))
    if not cells:
        return []
    q = np.unique(np.concatenate([c["ii"] for c in cells]))
    for c in cells:
        c["pos"] = np.searchsorted(q, c["ii"])
    Xq = sf.coords[q]
    rec = []
    for nm in NULLS:
        for r in range(n_draw):
            param = VG.draw(nm, rng)
            back = VG.pullback(nm, Xq, param)
            keep = VG.in_W(back)
            dd = np.full(q.size, np.nan)
            if keep.any():
                dd[keep] = stree.query(back[keep], k=1, workers=2)[0]
            for c in cells:
                m = keep[c["pos"]]
                ni = int(m.sum())
                if ni < V.MIN_N_RETAINED:
                    continue
                yv = c["y"][m]
                k = np.exp(-dd[c["pos"][m]] / c["lam"])
                rec.append((sample, nm, c["t"], c["mod"], r, ni,
                            V.sample_cov(k, yv)))
    print(f"[var-check] {sample} draws={n_draw} rows={len(rec)}", flush=True)
    return rec


def main(sections=None, call=RN.PRIMARY_CALL, n_draw=400, n_cells=4):
    rep = SC1.reportable(call)
    sections = sections or P.IN_BAND[:2]
    rows = []
    for i, s in enumerate(sections):
        sub = rep[rep.section == s].head(n_cells)
        wanted = {(r.section, r.celltype, r.module) for r in sub.itertuples()}
        rows += _job(s, call, n_draw, P.MASTER_SEED + 5000 * i, wanted)
    d = pd.DataFrame(rows, columns=["section", "null", "celltype", "module",
                                    "draw", "n_i", "T_i"])
    out = []
    for (sec, nm, ct, mod), g in d.groupby(["section", "null", "celltype",
                                            "module"]):
        g = g[np.isfinite(g.T_i)]
        if len(g) < 50:
            continue
        Tbar = g.T_i.mean()
        rho, prho = spearmanr(g.n_i, np.abs(g.T_i - Tbar))
        qs = np.quantile(g.n_i, np.linspace(0, 1, N_BIN + 1))
        qs[-1] += 1
        b = np.clip(np.digitize(g.n_i, qs[1:-1]), 0, N_BIN - 1)
        xs, ys = [], []
        for k in range(N_BIN):
            gg = g[b == k]
            if len(gg) < 15:
                continue
            xs.append(np.log(gg.n_i.mean()))
            ys.append(np.log(gg.T_i.std(ddof=1)))
        slope = (np.polyfit(xs, ys, 1)[0] if len(xs) >= 3 else np.nan)
        out.append(dict(section=sec, null=nm, celltype=ct, module=mod,
                        n_draws=int(len(g)),
                        n_i_min=int(g.n_i.min()), n_i_max=int(g.n_i.max()),
                        spearman_n_vs_absdev=float(rho), spearman_p=float(prho),
                        log_sd_vs_log_n_slope=float(slope)))
    df = pd.DataFrame(out)
    df.to_csv(f"{RES}/var_variance_check.csv", index=False)
    print(df.to_string(index=False))
    print("\nmedian slope of log sd(T) on log n (RS_count predicts -0.5): "
          f"{df.log_sd_vs_log_n_slope.median():.3f}")
    print("median Spearman rho(n_i, |T_i - Tbar|) (predicted negative): "
          f"{df.spearman_n_vs_absdev.median():.3f}")
    return df


if __name__ == "__main__":
    main()
