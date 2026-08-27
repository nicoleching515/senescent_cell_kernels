#!/usr/bin/env python3
"""Phase 7 / C1 — how destructive is each Phase 3 coordinate null, physically?

The Phase 4 diagnostic (`code/phase4_diag.py`), applied to the Phase 3 sections
and the corrected N3/N4 variants:

  * frac_retaining_a_neighbour -- fraction of SHIFTED senders that still have at
    least one real cell within the 100 um fitting window.  A shifted sender in
    the void has no receivers and contributes nothing, which is exactly the
    artifact `CS_PHASE4.md` §2.4 flagged and Phase 3 never fixed.
  * real_median_nbrs / null_median_nbrs -- median number of real cells within
    100 um of a sender, before and after the null.  Phase 4 reports 150.8 vs
    149.4 on its tiles.
  * frac_in_occupancy -- fraction landing in an occupied 25 um grid cell.
  * median_displacement_um -- how far the null actually MOVED each sender.  A
    null that keeps every sender in tissue by not moving it is not a null; this
    column is what separates the two.

-> /workspace/results/phase3/null_destructiveness.csv
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import phase3_null_geom as G

WINDOW_UM = 100.0
N_REP = 20
CALL = "tierA_p95"
SEED = 20260820


def section_rows(sample, call=CALL, n_rep=N_REP, seed=SEED):
    sec = P.Sec(sample)
    xy = sec.coords.astype(float)
    snd = sec.sender_mask(call)
    elig = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    g = G.Geom(xy, snd, elig)
    tree = g.tree
    rng = np.random.default_rng(seed + int(sample[:4]))

    base_all = np.asarray(tree.query_ball_point(xy[snd], WINDOW_UM, workers=-1,
                                                return_length=True), float) - 1.0
    ti = g.tile_sender_idx()
    base_tile = np.asarray(tree.query_ball_point(xy[ti], WINDOW_UM, workers=-1,
                                                 return_length=True), float) - 1.0
    rows = []
    for null in G.ALL_NULLS:
        tiled = null in G.TILE_NULLS
        src = xy[ti] if tiled else xy[snd]
        base = base_tile if tiled else base_all
        keep, med, occf, disp = [], [], [], []
        for _ in range(n_rep):
            pts = g.draw(null, rng)
            cnt = np.asarray(tree.query_ball_point(pts, WINDOW_UM, workers=-1,
                                                   return_length=True), float)
            keep.append(float((cnt > 0).mean()))
            med.append(float(np.median(cnt)))
            occf.append(g._in_occ(pts))
            disp.append(float(np.median(np.hypot(*(pts - src).T))))
        rows.append(dict(
            section=sample, arm=sec.meta["condition"], week=sec.meta["week"],
            call=call, null=null,
            scope="solid tiles" if tiled else "whole section",
            n_senders=int(len(src)), n_cells=int(sec.n),
            occ_frac_bbox=round(g.occ_frac, 4),
            tissue_frac_bbox=round(g.tissue_frac, 4),
            n_solid_tiles=len(g.tiles),
            tile_cell_coverage=round(g.tile_cov_cells, 4),
            tile_sender_coverage=round(g.tile_cov_senders, 4),
            real_median_nbrs=float(np.median(base)),
            null_median_nbrs=float(np.mean(med)),
            frac_retaining_a_neighbour=float(np.mean(keep)),
            frac_in_occupancy=float(np.mean(occf)),
            median_displacement_um=float(np.mean(disp)),
            n_admissible_moves=(len(g._acc_off(G.OCC_TOL)) if null == "N3_occ"
                                else len(g._acc_off(G.OCC_TOL_RELAXED)) if null == "N3_occ15"
                                else len(g._acc_ang(G.OCC_TOL)) if null == "N4_occ"
                                else len(g._acc_ang(G.OCC_TOL_RELAXED)) if null == "N4_occ15"
                                else np.nan),
            n_candidate_moves=(g.nx * g.ny if null in ("N3_occ", "N3_occ15")
                               else G.N_ANGLE if null in ("N4_occ", "N4_occ15")
                               else np.nan)))
        print(f"  {sample[:4]} {null:9s} keep={rows[-1]['frac_retaining_a_neighbour']:.3f} "
              f"nbrs {rows[-1]['real_median_nbrs']:.1f}->{rows[-1]['null_median_nbrs']:.1f} "
              f"occ={rows[-1]['frac_in_occupancy']:.3f} "
              f"disp={rows[-1]['median_displacement_um']:.0f}um", flush=True)
    return rows


def main(sections=None):
    sections = sections or P.IN_BAND
    rows = []
    for s in sections:
        print(f"[diag] {s}", flush=True)
        rows += section_rows(s)
    df = pd.DataFrame(rows)
    out = f"{P.RESULTS}/null_destructiveness.csv"
    df.to_csv(out, index=False)
    print("\nwrote", out)
    print(df.groupby("null", sort=False)[
        ["frac_retaining_a_neighbour", "real_median_nbrs", "null_median_nbrs",
         "frac_in_occupancy", "median_displacement_um"]].median().round(3).to_string())
    return df


if __name__ == "__main__":
    main()
