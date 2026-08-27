#!/usr/bin/env python3
"""Phase 9 test A6 — the arm-specific anatomical covariate for H1: RED PULP vs WHITE PULP.

§13's A6 text is written for lung (airway-to-alveolar) and is void for this arm
(reports/PHASE7_H1_SCREEN.md; code/spec_a6_compartments_human.py).  The spleen analogue of
liver zonation is the red-pulp / white-pulp axis, and the gene side of it was frozen before
the tag as the five `genesets/human/D_spleen_*.txt` sets.  This script does the half that
needed expression and could not be done before the freeze: build the covariate and validate it.

CONSTRUCTION, mirroring `code/phase2_downstream.py`'s D-B block one-for-one:
  mouse:  zon   = score(D_zonation_pericentral) - score(D_zonation_periportal),
          z-scored on HEPATOCYTES, tertiles -> periportal / midzonal / pericentral
  H1:     pulp  = score(D_spleen_red_pulp) - mean(score(follicle), score(tzone)),
          z-scored on ALL analysis cells, tertiles -> white_pulp / intermediate / red_pulp
  Declared difference: the mouse score is standardised on the one parenchymal type that
  carries the axis (hepatocytes).  Spleen has no such single type -- red pulp and white pulp
  are made of different cells -- so the standardisation is over all analysis cells.  Using any
  one cell type would define the axis out of existence.

  score_genes ctrl_size = 200, identical to phase2_downstream.py.

ALSO BUILT (the other Tier D anatomical columns the N5 block needs):
  dist_to_boundary_um   occupancy grid -> close/fill/open -> EDT, GRID_UM = 25.0, verbatim
  dist_to_follicle_um   DBSCAN(eps=30, min_samples=10) on follicular + germinal-centre B
                        cells -> follicle centroids -> distance.  This is the structural
                        analogue of the mouse `dist_to_portal_triad_um` (DBSCAN on
                        biliary/ductular cells) and occupies the same cache slot.

VALIDATION (A6 is a "build the covariate before fitting anything" test; a covariate that does
not track the anatomy is not a covariate):
  V1  corr(pulp axis, dist_to_follicle) within all cells -- expect POSITIVE, red pulp is far
      from follicles.  This is the exact analogue of the mouse check
      corr(zonation, dist_to_portal_triad) > 0.
  V2  mean pulp axis by the DEPOSITORS' Level_3 annotation.  External, marker-free-for-us
      validation: red pulp macrophages / red pulp fibroblasts / sinusoidal endothelium must
      sit at the red end, B cells / T cells / reticular cells at the white end.
  V3  the same by our own Job B label.
  V4  spatial structure: Moran's I of the axis on the 20-NN graph.  A real anatomical axis is
      strongly spatially autocorrelated; a depth artefact need not be.
  V5  correlation with log transcript counts -- reported against interest, because on a 5K
      panel any score correlates with depth and the mouse E2 control failed on exactly this.

Usage: python3 code/h1_a6_compartments.py SPLN07 [...]
Writes data/processed_h1/anatomy_h1_<sec>.csv and results/phase9_h1/a6_validation*.csv
"""
import sys, os, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.spatial import cKDTree
from scipy import ndimage
from sklearn.cluster import DBSCAN
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 32; sc.settings.verbosity = 0
sys.path.insert(0, "/workspace/code")
import h1_common as H

PROC = H.PROC + "/"
GRID_UM = 25.0
COMPARTMENTS = ["D_spleen_red_pulp", "D_spleen_white_pulp_follicle",
                "D_spleen_white_pulp_tzone", "D_spleen_marginal_zone",
                "D_spleen_capsule_trabecula"]
FOLLICLE_TYPES = ("Follicular B cells", "Germinal centre B cells")


def run(section):
    print("\n" + "=" * 90); print("A6", section); print("=" * 90, flush=True)
    ctdf = pd.read_csv(PROC + "celltypes_h1_%s.csv" % section).set_index("cell_id")
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = pd.Index(bc).isin(ctdf.index)
    cells = H.cells_table(section).set_index("cell_id")
    B = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]), var=pd.DataFrame(index=names))
    B.layers["counts"] = B.X.copy()
    sc.pp.normalize_total(B); sc.pp.log1p(B)
    ctdf = ctdf.reindex(B.obs_names)
    ct = ctdf["cell_type"].astype(str).to_numpy()
    xy = cells.reindex(B.obs_names)[["x_centroid", "y_centroid"]].to_numpy(float)
    tc = cells.reindex(B.obs_names)["transcript_counts"].to_numpy(float)
    n = B.n_obs

    sco = {}
    for c in COMPARTMENTS:
        g = [x for x in H.gl(c) if x in B.var_names]
        sc.tl.score_genes(B, g, score_name="_" + c, ctrl_size=200)
        sco[c] = B.obs["_" + c].to_numpy()
        print("  %-32s %2d/%2d genes on panel" % (c, len(g), len(H.gl(c))))

    white = 0.5 * (sco["D_spleen_white_pulp_follicle"] + sco["D_spleen_white_pulp_tzone"])
    pulp = sco["D_spleen_red_pulp"] - white
    pulp_z = (pulp - pulp.mean()) / pulp.std()
    q = np.quantile(pulp_z, [1 / 3, 2 / 3])
    comp = np.where(pulp_z <= q[0], "white_pulp",
                    np.where(pulp_z >= q[1], "red_pulp", "intermediate"))

    # spatially smoothed sensitivity (20-NN mean, self excluded)
    tree = cKDTree(xy); _, nnidx = tree.query(xy, k=21, workers=-1)
    pulp_smooth = pulp_z[nnidx[:, 1:]].mean(1)

    # distance to tissue boundary
    x0, y0 = xy.min(0) - GRID_UM
    nx = int((xy[:, 0].max() + GRID_UM - x0) / GRID_UM) + 1
    ny = int((xy[:, 1].max() + GRID_UM - y0) / GRID_UM) + 1
    ix = ((xy[:, 0] - x0) / GRID_UM).astype(int); iy = ((xy[:, 1] - y0) / GRID_UM).astype(int)
    occ = np.zeros((nx, ny), bool); occ[ix, iy] = True
    occ = ndimage.binary_fill_holes(ndimage.binary_closing(occ, np.ones((3, 3))))
    occ = ndimage.binary_opening(occ, np.ones((3, 3)))
    dist_bound = (ndimage.distance_transform_edt(occ) * GRID_UM)[ix, iy]

    # distance to nearest follicle centre
    foll = np.where(np.isin(ct, FOLLICLE_TYPES))[0]
    dist_f = np.full(n, np.nan); nf = 0
    if len(foll) >= 30:
        lab = DBSCAN(eps=30, min_samples=10).fit(xy[foll]).labels_
        if (lab >= 0).any():
            cents = np.array([xy[foll][lab == l].mean(0) for l in np.unique(lab[lab >= 0])])
            nf = len(cents)
            dist_f = cKDTree(cents).query(xy, k=1)[0]

    # ---- validation -------------------------------------------------------
    v1 = float(np.corrcoef(pulp_z, dist_f)[0, 1]) if nf else np.nan
    v1s = float(np.corrcoef(pulp_smooth, dist_f)[0, 1]) if nf else np.nan
    # Moran's I on the 20-NN graph
    z = pulp_z - pulp_z.mean()
    moran = float(n * (z * z[nnidx[:, 1:]].sum(1)).sum() / (20 * n * (z * z).sum()) * n / n)
    moran = float((z * z[nnidx[:, 1:]].mean(1)).sum() / (z * z).sum())
    v5 = float(np.corrcoef(pulp_z, np.log1p(tc))[0, 1])
    print("  boundary dist median %.1f um | follicle foci %d | dist median %.1f um"
          % (np.median(dist_bound), nf, np.nanmedian(dist_f) if nf else np.nan))
    print("  V1 corr(pulp axis, dist_to_follicle) = %+.3f  (expect POSITIVE)" % v1)
    print("  V1s same, 20-NN smoothed axis        = %+.3f" % v1s)
    print("  V4 Moran's I (20-NN) of the axis     = %+.3f" % moran)
    print("  V5 corr(pulp axis, log counts)       = %+.3f  (reported against interest)" % v5)

    pd.DataFrame({"cell_id": B.obs_names, "cell_type": ct,
                  "zonation_score": np.round(pulp_z, 4),
                  "pulp_axis_smooth": np.round(pulp_smooth, 4),
                  "compartment_label": comp,
                  "dist_to_boundary_um": np.round(dist_bound, 2),
                  "dist_to_portal_triad_um": np.round(dist_f, 2),
                  **{"score_" + c: np.round(v, 5) for c, v in sco.items()}}
                 ).to_csv(PROC + "anatomy_h1_%s.csv" % section, index=False)

    # V2 / V3 tables
    ann = H.annotations(section).set_index("cell_id").reindex(B.obs_names)
    rows = []
    for lvl, lab in (("depositor_L3", ann["Level_3_Annotations"].astype(str).to_numpy()),
                     ("jobB_fine", ct)):
        for t in pd.unique(lab):
            m = lab == t
            if m.sum() < 200:
                continue
            rows.append(dict(section=section, label_source=lvl, label=t, n=int(m.sum()),
                             mean_pulp_axis=round(float(pulp_z[m].mean()), 3),
                             mean_pulp_axis_smooth=round(float(pulp_smooth[m].mean()), 3),
                             median_dist_to_follicle_um=round(float(np.nanmedian(dist_f[m])), 1)
                             if nf else None))
    val = pd.DataFrame(rows)
    val.to_csv(H.RESULTS + "/a6_validation_by_label_%s.csv" % section, index=False)
    summ = dict(section=section, n_cells=int(n), n_follicle_foci=int(nf),
                v1_corr_axis_distfollicle=round(v1, 4) if nf else None,
                v1s_corr_smooth_distfollicle=round(v1s, 4) if nf else None,
                v4_moran_i_20nn=round(moran, 4),
                v5_corr_axis_logcounts=round(v5, 4),
                median_dist_boundary_um=round(float(np.median(dist_bound)), 2),
                frac_red=round(float((comp == "red_pulp").mean()), 4),
                frac_white=round(float((comp == "white_pulp").mean()), 4))
    print("  wrote anatomy_h1_%s.csv" % section, flush=True)
    del B
    return summ


if __name__ == "__main__":
    out = [run(s) for s in sys.argv[1:]]
    p = H.RESULTS + "/a6_summary.csv"
    d = pd.DataFrame(out)
    if os.path.exists(p):
        old = pd.read_csv(p)
        d = pd.concat([old[~old.section.isin(d.section)], d])
    d.sort_values("section").to_csv(p, index=False)
    print(d.to_string(index=False)); print("wrote", p)
