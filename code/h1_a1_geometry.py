#!/usr/bin/env python3
"""Phase 9 test A1 — resolution, segmentation and cell-level QC for the H1 human spleen arm,
plus the characterisation of the depositors' annotation QC filter (PREREG_PHASE8.md P21).

A1 asks for: single-cell resolution confirmed, median nearest-neighbour distance
(M1: 6.74-10.61 um over 11 sections) and the transcript assignment rate (M1: 88.27 % on one
admissible section).  The assignment rate needs `transcripts.parquet`, which §12.3 deliberately
did NOT download; `code/h1_a1_assignment.py` fetches it for a subset of sections and is a
declared deviation.  This script covers everything obtainable from the acquired files.

Writes results/phase9_h1/a1_sections.csv, a1_segmentation.csv, a1_depositor_annotation.csv.
Reads nothing under data/raw/ and writes nothing under results/phase3/.
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
sys.path.insert(0, "/workspace/code")
import h1_common as H

os.makedirs(H.RESULTS, exist_ok=True)
os.makedirs(H.PROC, exist_ok=True)

rows, segrows, annrows = [], [], []
for s in H.ALL_SECTIONS:
    gsm, age, sex = H.SECTIONS[s]
    c = H.cells_table(s)
    xy = np.column_stack([c.x_centroid.to_numpy(float), c.y_centroid.to_numpy(float)])
    tree = cKDTree(xy)
    nn, _ = tree.query(xy, k=2, workers=-1)
    med_nn_all = float(np.median(nn[:, 1]))

    # QC-passed set: the frozen annotate_pipeline rule, >=20 counts and >=5 genes.
    # genes-detected needs the matrix; read it once here.
    X, names, bc, _ = H.load_matrix(s, "gene")
    assert np.array_equal(bc, c.cell_id.to_numpy().astype(str)), s
    tot = np.asarray(X.sum(1)).ravel()
    ng = np.asarray((X > 0).sum(1)).ravel()
    qc = (tot >= H.MIN_COUNTS) & (ng >= H.MIN_GENES)
    tq = cKDTree(xy[qc]); nnq, _ = tq.query(xy[qc], k=2, workers=-1)
    med_nn_qc = float(np.median(nnq[:, 1]))

    ann = H.annotations(s)
    covered = c.cell_id.isin(set(ann.cell_id)).to_numpy()
    lowq = ann.Level_1_Annotations.eq("Low quality").to_numpy()

    rows.append(dict(
        section=s, gsm=gsm, age=age, sex=sex,
        n_cells=len(c),
        n_qc_pass=int(qc.sum()), qc_pass_pct=round(100 * qc.mean(), 2),
        median_nn_um_all=round(med_nn_all, 3),
        median_nn_um_qc=round(med_nn_qc, 3),
        p05_nn_um=round(float(np.percentile(nn[:, 1], 5)), 3),
        p95_nn_um=round(float(np.percentile(nn[:, 1], 95)), 3),
        median_transcripts=float(np.median(c.transcript_counts)),
        median_genes_detected=float(np.median(ng)),
        median_cell_area_um2=round(float(np.median(c.cell_area)), 2),
        median_nucleus_area_um2=round(float(np.median(c.nucleus_area)), 2),
        median_nucleus_count=float(np.median(c.nucleus_count)),
        pct_multinucleate=round(100 * float((c.nucleus_count > 1).mean()), 2),
        pct_zero_nucleus=round(100 * float((c.nucleus_count == 0).mean()), 2),
        x_range_um=round(float(xy[:, 0].max() - xy[:, 0].min()), 1),
        y_range_um=round(float(xy[:, 1].max() - xy[:, 1].min()), 1),
        # control-feature burden, the only in-file readout related to mis-assignment
        pct_counts_neg_probe=round(100 * c.control_probe_counts.sum() / c.total_counts.sum(), 4),
        pct_counts_neg_codeword=round(100 * c.control_codeword_counts.sum() / c.total_counts.sum(), 4),
        pct_counts_genomic_control=round(100 * c.genomic_control_counts.sum() / c.total_counts.sum(), 4),
        pct_counts_unassigned_codeword=round(100 * c.unassigned_codeword_counts.sum() / c.total_counts.sum(), 4),
        pct_counts_deprecated_codeword=round(100 * c.deprecated_codeword_counts.sum() / c.total_counts.sum(), 4),
    ))
    for m, k in c.segmentation_method.value_counts().items():
        segrows.append(dict(section=s, segmentation_method=m, n=int(k),
                            pct=round(100 * k / len(c), 2)))
    annrows.append(dict(
        section=s, n_cells_matrix=len(c), n_cells_annotated=len(ann),
        pct_annotated=round(100 * len(ann) / len(c), 2),
        n_annotated_lowquality=int(lowq.sum()),
        pct_lowquality_of_annotated=round(100 * lowq.mean(), 2),
        # what the depositors' filter drops, in our terms
        median_transcripts_annotated=float(np.median(c.transcript_counts[covered])),
        median_transcripts_dropped=float(np.median(c.transcript_counts[~covered])),
        median_area_annotated=round(float(np.median(c.cell_area[covered])), 2),
        median_area_dropped=round(float(np.median(c.cell_area[~covered])), 2),
        pct_dropped_failing_our_qc=round(100 * float(qc[~covered].mean() == 0 and 1 or
                                                     (~qc[~covered]).mean()), 2),
        pct_our_qc_pass_among_dropped=round(100 * float(qc[~covered].mean()), 2),
        pct_our_qc_pass_among_annotated=round(100 * float(qc[covered].mean()), 2),
        n_annotated_failing_our_qc=int((~qc[covered]).sum()),
        n_levels=4, n_level3_labels=int(ann.Level_3_Annotations.nunique()),
        n_level4_labels=int(ann.Level_4_Annotations.nunique()),
    ))
    print(rows[-1], flush=True)
    del X

pd.DataFrame(rows).to_csv(H.RESULTS + "/a1_sections.csv", index=False)
pd.DataFrame(segrows).to_csv(H.RESULTS + "/a1_segmentation.csv", index=False)
pd.DataFrame(annrows).to_csv(H.RESULTS + "/a1_depositor_annotation.csv", index=False)
d = pd.DataFrame(rows)
print("\n=== A1 SUMMARY ===")
print("sections %d  cells %d  QC-pass %d (%.2f%%)"
      % (len(d), d.n_cells.sum(), d.n_qc_pass.sum(), 100 * d.n_qc_pass.sum() / d.n_cells.sum()))
print("median NN (all cells)  %.2f - %.2f um" % (d.median_nn_um_all.min(), d.median_nn_um_all.max()))
print("median NN (QC-passed)  %.2f - %.2f um" % (d.median_nn_um_qc.min(), d.median_nn_um_qc.max()))
print("M1 reference           6.74 - 10.61 um ; lam grid floor is a literal 7.0 um")
print("wrote", H.RESULTS + "/a1_{sections,segmentation,depositor_annotation}.csv")
