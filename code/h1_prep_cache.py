#!/usr/bin/env python3
"""Phase 9 — build the H1 analogue of the Phase-3 section cache.

`code/sasp_phase3.prep` writes one .npz per mouse section and `sasp_phase3.Sec` reads it.
Every downstream estimator (`run_phase3_nulls.SectionFit`, `phase3_core.build_blocks`,
`match_decoys_section`, `neighbour_baseline`, `run_a7_control_probes`) is written against
`Sec` and nothing else, so writing the same .npz for H1 lets the audit tests A5 and A7 reuse
the FROZEN estimator byte-for-byte instead of a reimplementation.

The cache goes to data/processed_h1/cache3_h1/ ; data/processed/cache3/ is never touched.
The mouse `zonation_score` / `compartment` / `dist_to_portal_triad_um` slots carry the H1
red/white-pulp axis, its tertile label, and distance-to-nearest-follicle (test A6).

Usage: python3 code/h1_prep_cache.py SPLN07 [...]
"""
import sys, os
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
sys.path.insert(0, "/workspace/code")
import h1_common as H

CACHE = H.PROC + "/cache3_h1"


def prep(section, force=False):
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, section + ".npz")
    if os.path.exists(out) and not force:
        return "[skip] " + section
    cells = H.cells_table(section)
    X, gene_names, bc, _ = H.load_matrix(section, "gene")
    assert np.array_equal(bc, cells.cell_id.to_numpy().astype(str))

    def _read(kind):
        return pd.read_csv(H.PROC + "/%s_h1_%s.csv" % (kind, section)).set_index("cell_id")

    ct, mods, sen, anat = (_read(k) for k in ("celltypes", "modules", "senders", "anatomy"))
    cid = cells.cell_id.to_numpy().astype(str)
    keep = pd.Index(cid).isin(ct.index)
    cells = cells.loc[keep].reset_index(drop=True)
    X = X[keep]; cid = cid[keep]
    ct, mods, sen, anat = (t.reindex(cid) for t in (ct, mods, sen, anat))

    coords = np.column_stack([cells.x_centroid.to_numpy(float),
                              cells.y_centroid.to_numpy(float)])
    tree = cKDTree(coords)
    d = {}
    for r in (25.0, 50.0, 100.0):
        d["density_%dum" % int(r)] = np.asarray(
            tree.query_ball_point(coords, r=r, return_length=True, workers=-1),
            np.float32) - 1.0
    nn, nn_idx = tree.query(coords, k=21, workers=-1)
    d["nn1_um"] = nn[:, 1].astype(np.float32)
    d["knn_idx"] = nn_idx[:, 1:].astype(np.int32)
    d["median_nn_um"] = np.float64(np.median(nn[:, 1]))
    genes_det = np.asarray((X > 0).sum(axis=1)).ravel().astype(np.float32)
    segm = cells.segmentation_method.to_numpy().astype(str)
    seg_levels = np.array(sorted(set(segm)))
    d["seg_levels"] = seg_levels
    d["seg_code"] = np.searchsorted(seg_levels, segm).astype(np.int8)
    # nucleus_area is NaN for the 0.5-1.3 % of cells with no segmented nucleus;
    # build_blocks does log1p on it, so NaN would poison the N5 block.  Median-fill,
    # and record how many cells that touched.
    na = cells.nucleus_area.to_numpy(np.float32)
    n_na = int(np.isnan(na).sum())
    na = np.where(np.isnan(na), np.nanmedian(na), na)

    d.update(cell_id=cid, coords=coords.astype(np.float32),
             celltype=ct["cell_type"].to_numpy().astype(str),
             celltype_merged=ct["cell_type_merged"].to_numpy().astype(str),
             median_depth=np.float64(np.median(cells.transcript_counts.to_numpy(float))),
             ct_conf=ct["cell_type_confidence"].to_numpy(np.float32),
             transcript_counts=cells.transcript_counts.to_numpy(np.float32),
             genes_detected=genes_det,
             cell_area=cells.cell_area.to_numpy(np.float32),
             nucleus_area=na,
             n_nucleus_area_imputed=np.int64(n_na),
             zonation_score=anat["zonation_score"].to_numpy(np.float32),
             compartment=anat["compartment_label"].to_numpy().astype(str),
             dist_to_boundary_um=anat["dist_to_boundary_um"].to_numpy(np.float32),
             dist_to_portal_triad_um=anat["dist_to_portal_triad_um"].to_numpy(np.float32),
             cdkn1a_counts=sen["cdkn1a_counts"].to_numpy(np.float32),
             tierA_score=sen["tierA_score"].to_numpy(np.float32),
             senepy_score=sen["senepy_score"].to_numpy(np.float32))
    for q in (90, 95, 99):
        d["flag_p%d" % q] = sen["sender_flag_p%d" % q].fillna(0).to_numpy().astype(bool)
    for m in H.MODULES:
        d["mod__" + m] = mods[m].to_numpy(np.float32)
    np.savez_compressed(out, **d)
    return ("[done] %s n=%d medNN=%.2f seg_levels=%d nucarea_imputed=%d"
            % (section, len(cid), float(d["median_nn_um"]), len(seg_levels), n_na))


if __name__ == "__main__":
    for s in sys.argv[1:]:
        print(prep(s), flush=True)
