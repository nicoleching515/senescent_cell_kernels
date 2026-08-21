#!/usr/bin/env python3
"""
How much of the between-section spread in lambda_hat is just sender density?

Companion to `run_phase3_poisson.py`.  For every (section, sender call,
receiver type, module) the exponential kernel is fitted on three distance
scales, all with the same relative window and the same relative grid:

  raw          d in microns
  poisson      d / d_pois, where d_pois = sqrt(ln2/(pi rho)) is the median
               nearest-sender distance a homogeneous Poisson process of the
               observed sender density would give -- i.e. distance measured in
               units of "how far apart the called senders happen to be"
  nn           d / median nearest-neighbour distance of the section (packing)

If lambda in microns varies across sections but lambda in Poisson units does
not, then the microns were reporting the sender calling rate.
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import run_phase3_poisson as PO

CALLS = ["tierA_p95", "cdkn1a_pos", "senepy_p95"]
CELLTYPES = ["Hepatocytes", "LSECs", "Kupffer cells", "Hepatic stellate cells",
             "Biliary/ductular"]
REL_GRID = np.exp(np.linspace(np.log(0.07), np.log(0.50), 40))  # x window


def job(sample, call):
    sec = P.Sec(sample)
    co = sec.coords.astype(float)
    area = PO.tissue_area_um2(co)
    snd = sec.sender_mask(call)
    if snd.sum() < 30:
        return []
    rho = snd.sum() / area
    d_pois = float(np.sqrt(np.log(2.0) / (np.pi * rho)))
    med_nn = sec.median_nn_um
    d = P.dist_to_senders(co, snd)
    Y = np.column_stack([sec.module(m) for m in P.MODULES])
    rows = []
    for ct in CELLTYPES:
        idx = ((sec.celltype == ct) & (~snd) & np.isfinite(d)
               & (d <= RN.WINDOW_UM))
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        X = np.ones((ii.size, 1))
        for j, mod in enumerate(P.MODULES):
            y = Y[ii, j].astype(float)
            r = dict(section=sample, arm=sec.meta["condition"],
                     band=("in_band" if sample in P.IN_BAND else
                           "over_ceiling" if sample in P.OVER_CEILING
                           else "below_floor"),
                     call=call, celltype=ct, module=mod, n=int(ii.size),
                     sender_density_per_um2=rho, d_poisson_um=d_pois,
                     median_nn_um=med_nn, sd_y=float(y.std()))
            for scale, unit in (("raw", 1.0), ("poisson", d_pois),
                                ("nn", med_nn)):
                dd = d[ii] / unit
                grid = REL_GRID * (RN.WINDOW_UM / unit)
                lam, beta, rss, t = P.profile_lambda(dd, X, y, grid)
                r[f"lam_{scale}"] = lam
                r[f"beta_{scale}"] = beta
                r[f"railed_{scale}"] = int(t in (0, grid.size - 1))
            rows.append(r)
    print(f"[lamscale] {sample} {call} {len(rows)}", flush=True)
    return rows


def report(df):
    def spread(sub, col):
        v = np.log(sub[col].to_numpy(float))
        v = v[np.isfinite(v)]
        return float(v.std())

    print("\nBetween-section spread of log lambda (sd), by scale, "
          "over in-band sections:")
    print(f"{'call':12s} {'celltype':24s} {'module':22s} "
          f"{'raw':>7s} {'poisson':>8s} {'nn':>7s}")
    rows = []
    ib = df[df.band == "in_band"]
    for (call, ct, mod), g in ib.groupby(["call", "celltype", "module"]):
        if g.section.nunique() < 4:
            continue
        a, b, c = (spread(g, "lam_raw"), spread(g, "lam_poisson"),
                   spread(g, "lam_nn"))
        rows.append(dict(call=call, celltype=ct, module=mod,
                         n_sections=int(g.section.nunique()),
                         sd_log_lam_raw=a, sd_log_lam_poisson=b,
                         sd_log_lam_nn=c))
        print(f"{call:12s} {ct:24s} {mod:22s} {a:7.3f} {b:8.3f} {c:7.3f}")
    r = pd.DataFrame(rows)
    if not r.empty:
        print("\nMEDIAN over (call x celltype x module): "
              f"raw {r.sd_log_lam_raw.median():.3f}  "
              f"poisson {r.sd_log_lam_poisson.median():.3f}  "
              f"nn {r.sd_log_lam_nn.median():.3f}")
    # R^2 of log lambda_raw on log sender density
    from numpy.linalg import lstsq
    print("\nR^2 of log(lambda_raw) on log(sender density), pooled within "
          "celltype x module across ALL sections and calls:")
    out = []
    for (ct, mod), g in df.groupby(["celltype", "module"]):
        g = g[np.isfinite(g.lam_raw)]
        if len(g) < 8:
            continue
        x = np.log(g.sender_density_per_um2.to_numpy())
        for col in ("lam_raw", "lam_poisson"):
            y = np.log(g[col].to_numpy(float))
            A = np.column_stack([np.ones(x.size), x])
            b, *_ = lstsq(A, y, rcond=None)
            r2 = 1 - ((y - A @ b) ** 2).sum() / max(((y - y.mean()) ** 2).sum(),
                                                    1e-12)
            out.append(dict(celltype=ct, module=mod, scale=col, n=len(g),
                            slope=float(b[1]), r2=float(r2)))
    o = pd.DataFrame(out)
    piv = o.pivot_table(index=["celltype", "module"], columns="scale",
                        values=["slope", "r2"])
    print(piv.round(3).to_string())
    print("\nmedian r2: raw "
          f"{o[o.scale=='lam_raw'].r2.median():.3f}, poisson-normalised "
          f"{o[o.scale=='lam_poisson'].r2.median():.3f}")
    print("median slope on log density: raw "
          f"{o[o.scale=='lam_raw'].slope.median():+.3f} "
          "(a pure density readout would give -0.500), poisson-normalised "
          f"{o[o.scale=='lam_poisson'].slope.median():+.3f}")
    r.to_csv(f"{P.RESULTS}/lamscale_spread.csv", index=False)
    o.to_csv(f"{P.RESULTS}/lamscale_density_r2.csv", index=False)


if __name__ == "__main__":
    jobs = [(s, c) for s in P.ALL_SECTIONS for c in CALLS]
    out = Parallel(n_jobs=8, prefer="processes")(
        delayed(job)(s, c) for s, c in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(f"{P.RESULTS}/lamscale.csv", index=False)
    report(df)
