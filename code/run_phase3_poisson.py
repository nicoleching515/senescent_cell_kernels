#!/usr/bin/env python3
"""
Is distance-to-nearest-sender a measurement, or a readout of the sender calling
rate?

For a homogeneous Poisson process of intensity rho, the distance from an
arbitrary point to the nearest event has P(d > r) = exp(-rho pi r^2), so

    median d = sqrt(ln 2 / (pi rho)) = 0.4697 * rho^(-1/2)

i.e. log(median d) = const - 0.5 log(rho), EXACTLY.  If the observed
(section x sender-definition) points fall on that line, then the independent
variable of the whole kernel model is, to first order, a deterministic function
of how many cells were called senders -- which tracks sequencing depth -- and
not of where senescence is.

This script measures it three ways:
  1. log(median d) vs log(sender density), slope and r^2, against -0.50.
  2. the ratio observed / Poisson-predicted median d, which is a pure
     clustering diagnostic (1.0 = Poisson, >1 = clustered senders).
  3. how much of the between-section spread in lambda_hat is predicted by
     sender density alone, and what lambda looks like on a density-normalised
     distance scale.
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import run_phase3_nulls as RN

GRID_UM = 25.0
CALLS = ["tierA_p90", "tierA_p95", "tierA_p99", "cdkn1a_pos",
         "senepy_p90", "senepy_p95", "senepy_p99"]


def tissue_area_um2(coords):
    """Occupancy-grid area, the same construction Bio used for
    dist_to_boundary_um (25 um cells, closing + fill + opening)."""
    from scipy import ndimage
    x0, y0 = coords.min(0) - GRID_UM
    ix = ((coords[:, 0] - x0) / GRID_UM).astype(int)
    iy = ((coords[:, 1] - y0) / GRID_UM).astype(int)
    occ = np.zeros((ix.max() + 2, iy.max() + 2), bool)
    occ[ix, iy] = True
    occ = ndimage.binary_fill_holes(ndimage.binary_closing(occ, np.ones((3, 3))))
    occ = ndimage.binary_opening(occ, np.ones((3, 3)))
    return float(occ.sum()) * GRID_UM ** 2


def main():
    rows = []
    for s in P.ALL_SECTIONS:
        sec = P.Sec(s)
        co = sec.coords.astype(float)
        area = tissue_area_um2(co)
        med_depth = float(np.median(sec.transcript_counts))
        for call in CALLS:
            try:
                snd = sec.sender_mask(call)
            except Exception:
                continue
            if snd.sum() < 30:
                continue
            d = P.dist_to_senders(co, snd)
            recv = (~np.isin(sec.celltype, P.EXCLUDE_TYPES)) & (~snd)
            dr = d[recv]
            rho = snd.sum() / area
            med_pred = np.sqrt(np.log(2.0) / (np.pi * rho))
            rows.append(dict(
                section=s, arm=sec.meta["condition"], week=sec.meta["week"],
                band=("in_band" if s in P.IN_BAND else
                      "over_ceiling" if s in P.OVER_CEILING else "below_floor"),
                call=call, n_cells=sec.n, area_um2=area,
                median_depth=med_depth, median_nn_um=sec.median_nn_um,
                n_senders=int(snd.sum()),
                sender_prevalence=float(snd.mean()),
                sender_density_per_um2=rho,
                median_d_um=float(np.median(dr)),
                iqr_d_um=float(np.percentile(dr, 75) - np.percentile(dr, 25)),
                p90_d_um=float(np.percentile(dr, 90)),
                poisson_median_pred_um=med_pred,
                obs_over_poisson=float(np.median(dr)) / med_pred))
            print(f"{s[:4]} {call:11s} rho={rho:.5f} med_d={np.median(dr):6.2f} "
                  f"pred={med_pred:6.2f} ratio={np.median(dr)/med_pred:.3f}",
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{P.RESULTS}/poisson_density.csv", index=False)

    def fit(sub, label):
        x = np.log(sub.sender_density_per_um2.to_numpy())
        y = np.log(sub.median_d_um.to_numpy())
        A = np.column_stack([np.ones(x.size), x])
        b, *_ = np.linalg.lstsq(A, y, rcond=None)
        r2 = 1 - ((y - A @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        print(f"\n{label}: n={x.size}  slope={b[1]:+.4f} (Poisson -0.5000)  "
              f"r2={r2:.4f}  intercept={b[0]:+.3f} "
              f"(Poisson {np.log(0.4697):+.3f})")
        print(f"   observed/Poisson median-d ratio: median "
              f"{sub.obs_over_poisson.median():.3f}, range "
              f"{sub.obs_over_poisson.min():.3f}-{sub.obs_over_poisson.max():.3f}")
        return dict(subset=label, n=int(x.size), slope=float(b[1]),
                    intercept=float(b[0]), r2=float(r2),
                    ratio_median=float(sub.obs_over_poisson.median()))

    out = [fit(df, "ALL sections x ALL sender definitions"),
           fit(df[df.call == "cdkn1a_pos"], "cdkn1a_pos only (across sections)"),
           fit(df[df.call.str.startswith("tierA")], "tierA percentile calls"),
           fit(df[df.band == "in_band"], "in-band sections only")]
    pd.DataFrame(out).to_csv(f"{P.RESULTS}/poisson_fits.csv", index=False)

    # depth -> sender prevalence
    c = df[df.call == "cdkn1a_pos"]
    from scipy.stats import spearmanr, pearsonr
    print("\nsection-level depth vs Cdkn1a+ prevalence: "
          f"spearman {spearmanr(c.median_depth, c.sender_prevalence).statistic:+.3f} "
          f"(n={len(c)})")
    print("section-level depth vs median d (cdkn1a_pos): "
          f"spearman {spearmanr(c.median_depth, c.median_d_um).statistic:+.3f}")
    print("\nwrote poisson_density.csv, poisson_fits.csv")


if __name__ == "__main__":
    main()
