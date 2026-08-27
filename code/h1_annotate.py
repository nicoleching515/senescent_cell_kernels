#!/usr/bin/env python3
"""Phase 9 Job B step 1 — cell-type annotation for the H1 human spleen arm.

This is `code/annotate_pipeline.py` transplanted to H1 with NO change to any threshold:
MIN_MARKERS/MIN_DET/MIN_Z/MIN_MARGIN/LOWQ_FRAC/RES/MIN_COUNTS/MIN_GENES/LOWQ_RESCUE_* are
imported from that file so they cannot drift.  Three things differ, each declared:

  1. Loader — human panel (5,093 ENSG Gene Expression features) instead of the mouse
     ENSMUSG filter, via `h1_common.load_matrix`.
  2. Label set — `markers_human_spleen.MARKERS` (the frozen 22) plus the pre-registered
     plasma-cell exception (PREREG_PHASE8.md P6 / PI decision D2b), applied in
     `h1_common.marker_set()` rather than by editing the frozen generated file.
  3. Sanity check — annotate_pipeline asserts a hepatocyte floor.  Spleen has no single
     dominant parenchymal type, so the assertion is replaced by a *reported* haematopoietic
     fraction and a cross-check against the depositors' own Level_1..Level_4 annotations
     (P21).  Nothing is asserted away.

Outputs go to data/processed_h1/ ; data/processed/ is never written.
Usage: python3 code/h1_annotate.py SPLN07 [SPLN14 ...]
"""
import sys, os, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 32; sc.settings.verbosity = 1
sys.path.insert(0, "/workspace/code")
import h1_common as H
from annotate_pipeline import (MIN_MARKERS, MIN_DET, MIN_Z, MIN_MARGIN,
                               LOWQ_FRAC, RES, MIN_COUNTS, MIN_GENES,
                               LOWQ_RESCUE_Z, LOWQ_RESCUE_MARGIN)

MARKERS, MERGE, MIN_MARKERS_EXCEPTION = H.marker_set()
PROC = H.PROC + "/"


def run(section):
    print("\n" + "=" * 95); print("SECTION %s" % section); print("=" * 95, flush=True)
    X, names, bc, _ = H.load_matrix(section, "gene")
    cells = H.cells_table(section).set_index("cell_id")
    A = ad.AnnData(X, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=names))
    A.obs = A.obs.join(cells.rename(columns={"x_centroid": "x_um", "y_centroid": "y_um"}))
    n0 = A.n_obs
    tot = np.asarray(A.X.sum(1)).ravel(); ng = np.asarray((A.X > 0).sum(1)).ravel()
    A = A[(tot >= MIN_COUNTS) & (ng >= MIN_GENES)].copy()
    print("QC %d -> %d cells (%.2f%%)" % (n0, A.n_obs, 100 * A.n_obs / n0), flush=True)
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A); sc.pp.log1p(A)
    A.raw = A
    sc.pp.scale(A, max_value=10)
    sc.tl.pca(A, n_comps=50, svd_solver="arpack")
    print("pca done", flush=True)
    sc.pp.neighbors(A, n_neighbors=15, n_pcs=50)
    print("neighbors done", flush=True)
    sc.tl.leiden(A, resolution=RES, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    B = A.raw.to_adata(); B.obs = A.obs.copy(); B.layers["counts"] = A.layers["counts"]
    del A
    cl = B.obs["leiden"].astype(str).values
    clusters = sorted(set(cl), key=lambda s: int(s))
    print("clusters: %d" % len(clusters), flush=True)

    avail, dropped = {}, {}
    for ct, gs in MARKERS.items():
        on = [g for g in gs if g in B.var_names]
        floor = MIN_MARKERS_EXCEPTION.get(ct, MIN_MARKERS)
        (avail if len(on) >= floor else dropped)[ct] = on
    print("\nLABEL SET GATE: %d types assignable, %d dropped" % (len(avail), len(dropped)))
    for k, v in dropped.items():
        print("  DROPPED %-30s only %d on-panel markers: %s" % (k, len(v), ",".join(v)))
    for k, f in MIN_MARKERS_EXCEPTION.items():
        if k in avail:
            print("  EXCEPTION %-28s admitted at %d markers (PREREG P6 / D2b): %s"
                  % (k, f, ",".join(avail[k])))

    cnts = np.asarray(B.layers["counts"].sum(1)).ravel()
    qc = pd.DataFrame({"cl": cl, "counts": cnts}).groupby("cl").median().loc[clusters, "counts"]
    lowq = set(qc.index[qc < LOWQ_FRAC * np.median(cnts)])
    print("\nLOW-QUALITY clusters (median counts < %.0f%% of section median %.0f): %s"
          % (100 * LOWQ_FRAC, np.median(cnts), sorted(lowq, key=int) or "none"))

    allg = sorted({g for v in avail.values() for g in v})
    Xd = np.asarray(B[:, allg].X.todense())
    dfm = pd.DataFrame(Xd, columns=allg); dfm["cl"] = cl
    cmean = dfm.groupby("cl", observed=True).mean().loc[clusters]
    cdet = pd.DataFrame((Xd > 0).astype(np.float32), columns=allg).assign(cl=cl) \
             .groupby("cl", observed=True).mean().loc[clusters]
    del Xd, dfm
    Zg = (cmean - cmean.mean()) / (cmean.std() + 1e-9)

    S = pd.DataFrame({c: Zg[v].mean(axis=1) for c, v in avail.items()})
    DET = pd.DataFrame({c: cdet[v].mean(axis=1) for c, v in avail.items()})
    Sg = S.where(DET >= MIN_DET, other=-np.inf)
    top = Sg.idxmax(axis=1); tv = Sg.max(axis=1)
    second = Sg.apply(lambda r: r.nlargest(2).iloc[1], axis=1); margin = tv - second
    lab = pd.Series(np.where((tv >= MIN_Z) & (margin >= MIN_MARGIN), top, "Unknown"), index=S.index)
    conf = pd.Series(np.round(np.clip(margin / (np.abs(tv) + 1e-9), 0, 1), 3), index=S.index)
    conf[lab == "Unknown"] = 0.0

    avail_m = {k: v for k, v in avail.items() if not any(k in mem for mem in MERGE.values())}
    for grp, mem in MERGE.items():
        u = sorted({g for c in mem for g in avail.get(c, [])})
        if any(c in avail for c in mem): avail_m[grp] = u
    Sm = pd.DataFrame({c: Sg[c] for c in avail_m if c in Sg.columns})
    for grp, mem in MERGE.items():
        mem = [c for c in mem if c in Sg.columns]
        if mem: Sm[grp] = Sg[mem].max(axis=1)
    top_m = Sm.idxmax(axis=1); tv_m = Sm.max(axis=1)
    sec_m = Sm.apply(lambda r: r.nlargest(2).iloc[1], axis=1); margin_m = tv_m - sec_m
    lab_m = pd.Series(np.where((tv_m >= MIN_Z) & (margin_m >= MIN_MARGIN), top_m, "Unknown"),
                      index=Sm.index)
    conf_m = pd.Series(np.round(np.clip(margin_m / (np.abs(tv_m) + 1e-9), 0, 1), 3), index=Sm.index)
    conf_m[lab_m == "Unknown"] = 0.0

    rescued = []
    for c in lowq:
        if np.isfinite(tv[c]) and tv[c] >= LOWQ_RESCUE_Z and margin[c] >= LOWQ_RESCUE_MARGIN:
            rescued.append((c, lab[c])); continue
        lab[c] = "Low_quality"; conf[c] = 0.0
    if rescued:
        print("  RESCUED from Low_quality on unambiguous identity: %s"
              % ", ".join("cl%s=%s" % (c, l) for c, l in rescued))
    for c in lowq:
        if not (np.isfinite(tv_m[c]) and tv_m[c] >= LOWQ_RESCUE_Z
                and margin_m[c] >= LOWQ_RESCUE_MARGIN):
            lab_m[c] = "Low_quality"; conf_m[c] = 0.0

    nper = pd.Series(cl).value_counts()
    print("\n%-4s %-8s %-9s %-30s %7s %7s %6s  %s"
          % ("cl", "n", "medcnt", "cell_type", "z", "margin", "det", "runner-up"))
    for c in clusters:
        print("%-4s %-8d %-9.0f %-30s %+7.2f %7.2f %6.2f  %s"
              % (c, nper[c], qc[c], lab[c],
                 tv[c] if np.isfinite(tv[c]) else -9,
                 margin[c] if np.isfinite(margin[c]) else 0,
                 DET.loc[c, top[c]], Sg.loc[c].nlargest(2).index[1]))

    B.obs["cell_type"] = lab.reindex(cl).values
    B.obs["cell_type_confidence"] = conf.reindex(cl).values
    B.obs["cell_type_merged"] = lab_m.reindex(cl).values
    B.obs["cell_type_merged_confidence"] = conf_m.reindex(cl).values

    print("\n--- merged composition (cell_type_merged) ---")
    for k, v in B.obs.cell_type_merged.value_counts().items():
        print("   %-30s %7d  %5.2f%%" % (k, v, 100 * v / len(B)))
    vc = B.obs.cell_type.value_counts()
    print("\n--- fine composition ---")
    for k, v in vc.items():
        print("   %-30s %7d  %5.2f%%" % (k, v, 100 * v / len(B)))
    frac = vc / len(B)
    haem = sum(frac.get(t, 0) for t in
               ("Red pulp macrophages", "Monocytes", "cDC1", "cDC2", "pDC",
                "Follicular B cells", "Marginal zone B cells", "Germinal centre B cells",
                "CD4 T cells", "CD8 T cells", "NK cells", "Erythroid cells",
                "Megakaryocytes", "Neutrophils", "Plasma cells"))
    print("\nSANITY (reported, NOT asserted): haematopoietic fraction %.1f%% ; "
          "Unknown %.1f%% ; Low_quality %.1f%%"
          % (100 * haem, 100 * frac.get("Unknown", 0), 100 * frac.get("Low_quality", 0)))
    print("median cell_type_confidence %.3f" % B.obs.cell_type_confidence.median())

    o = B.obs[["leiden", "cell_type", "cell_type_confidence", "cell_type_merged",
               "cell_type_merged_confidence"]].copy(); o.index.name = "cell_id"
    o.reset_index().to_csv(PROC + "celltypes_h1_%s.csv" % section, index=False)
    S.round(3).to_csv(PROC + "cluster_celltype_zscores_h1_%s.csv" % section)
    json.dump({"assignable": avail, "merged_label_set": avail_m, "merge_groups": MERGE,
               "dropped_thin": dropped, "min_markers_exception": MIN_MARKERS_EXCEPTION,
               "composition": {k: int(v) for k, v in vc.items()},
               "composition_merged": {k: int(v) for k, v in
                                      B.obs.cell_type_merged.value_counts().items()},
               "params": dict(MIN_MARKERS=MIN_MARKERS, MIN_DET=MIN_DET, MIN_Z=MIN_Z,
                              MIN_MARGIN=MIN_MARGIN, LOWQ_FRAC=LOWQ_FRAC, RES=RES,
                              MIN_COUNTS=MIN_COUNTS, MIN_GENES=MIN_GENES),
               "n_cells_qc": int(B.n_obs), "n_cells_raw": int(n0)},
              open(PROC + "annotation_meta_h1_%s.json" % section, "w"), indent=1)
    print("wrote celltypes_h1_%s.csv" % section, flush=True)
    del B


if __name__ == "__main__":
    for s in sys.argv[1:]:
        run(s)
