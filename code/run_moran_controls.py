#!/usr/bin/env python3
"""Moran's I on the Xenium negative-control features, M1 (mouse) arm.

WHY THIS EXISTS
---------------
`SASP_Kernel_Master_Plan.md` §29 objection 9 promises "our own Moran's I on the
controls alongside the kernel amplitude".  Nothing in the repo computed it
(`grep -ril moran` over code/ and results/ returned nothing).  The objection is
the designated response to the one piece of prior art that falsifies a project
claim: the Voyager Xenium vignette (Pachter lab) and Ren et al., Nat Commun 16
(2025), doi:10.1038/s41467-025-64292-3, both already compute Moran's I on
negative control probes and read a near-zero value as "no technical artifact
spatial trend".

The project's counter-claim is that a GLOBAL autocorrelation statistic and a
DISTANCE-TO-SENDER KERNEL are different questions.  This script tests that
claim on this data by computing both quantities on the same features and the
same cells and putting them side by side.

WHAT IT COMPUTES
----------------
Per section (11 M1 sections), on the SAME analysis cell set the Phase-3 /
A7 estimator uses (`sasp_phase3.Sec`, i.e. QC-passing cells carrying a Bio
cell-type label):

1. PER-FEATURE Moran's I on every one of the 13,590 h5 features individually
   -- 5,106 Gene Expression, 40 Negative Control Probe, 609 Negative Control
   Codeword, 21 Genomic Control, 7,806 Unassigned Codeword, 8 Deprecated
   Codeword.  This is the Voyager reproduction: it is the statistic Voyager
   actually plots, per feature, and it is the thing a reviewer holding the
   prior art will point at.

2. PER-CLASS AGGREGATE Moran's I on exactly the five A7 responses
   (`neg_control_probe`, `neg_control_codeword`, `genomic_control`,
   `all_controls`, `neg_probe_rate`), so the comparison against the A7 kernel
   amplitude is apples to apples -- same feature definition, same cells.

3. Positive controls: the seven Tier B biological module scores, plus
   `transcript_counts`, `cell_area` and local density as technical references.

4. Everything in (2) and (3) is also computed after removing the receiver
   cell-type means (`_ctcentred`), because the raw module I is inflated by the
   spatial arrangement of cell types, and A7/Phase-3 fits are within cell type.

SPATIAL WEIGHTS
---------------
Primary: k = 6 nearest neighbours, row-standardised.  6 is ~the number of
first-shell neighbours of a cell in a packed tissue and gives a mean neighbour
radius of ~20 um on these sections (median 1-NN 9.5-10 um).
Sensitivity: k in {4, 10, 20} and binary distance bands of 30, 50 and 100 um,
row-standardised.  The 100 um band is matched to the A7 fitting window, so the
"same scale" objection is answered directly.

Inference: analytic z under the randomisation assumption (Cliff & Ord), which
is what esda reports as `z_rand`/`p_rand`, plus a 999-replicate conditional
permutation p for the key features under the primary weights.  n ~ 2.4e5 per
section, so SE(I) ~ 2e-3 and essentially any non-zero I is "significant" --
the effect size, not the p, is what carries the argument.  That is stated in
the report.

VALIDATION
----------
`--validate` recomputes I with esda.Moran + libpysal.weights.KNN on a 4,000
cell subsample and asserts agreement with this module's vectorised routine.

Usage:
  run_moran_controls.py --validate
  run_moran_controls.py [--sections all|inband] [--n-jobs 3] [--perms 999]
"""
from __future__ import annotations

import os
import sys
import time
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from scipy.spatial import cKDTree
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P  # noqa: E402

RAW = "/workspace/data/raw"
OUT = "/workspace/results/moran"

# A7's five responses, and how they are built from the h5 feature classes.
CLASS_OF = {
    "neg_control_probe": "Negative Control Probe",
    "neg_control_codeword": "Negative Control Codeword",
    "genomic_control": "Genomic Control",
}
KNN_PRIMARY = 6
KNN_SENS = (4, 10, 20)
BANDS_UM = (30.0, 50.0, 100.0)


# ---------------------------------------------------------------------------
# weights
# ---------------------------------------------------------------------------
def knn_weights(coords: np.ndarray, k: int) -> sp.csr_matrix:
    """Row-standardised k-nearest-neighbour weights, self excluded."""
    n = coords.shape[0]
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1, workers=-1)
    idx = idx[:, 1:]                                   # drop self
    rows = np.repeat(np.arange(n, dtype=np.int64), k)
    W = sp.csr_matrix((np.full(n * k, 1.0 / k), (rows, idx.ravel())),
                      shape=(n, n))
    return W


def band_weights(coords: np.ndarray, r: float) -> tuple:
    """Row-standardised binary distance-band weights.  Returns (W, n_islands)."""
    tree = cKDTree(coords)
    A = tree.sparse_distance_matrix(tree, r, output_type="coo_matrix")
    A = A.tocsr()
    A.setdiag(0.0)
    A.eliminate_zeros()
    A.data[:] = 1.0
    rs = np.asarray(A.sum(1)).ravel()
    islands = int((rs == 0).sum())
    inv = np.where(rs > 0, 1.0 / np.maximum(rs, 1.0), 0.0)
    W = sp.diags(inv) @ A
    return W.tocsr(), islands


def w_moments(W: sp.csr_matrix) -> tuple:
    """S0, S1, S2 (Cliff & Ord)."""
    S0 = float(W.sum())
    B = W + W.T
    S1 = 0.5 * float(B.multiply(B).sum())
    rs = np.asarray(W.sum(1)).ravel()
    cs = np.asarray(W.sum(0)).ravel()
    S2 = float(((rs + cs) ** 2).sum())
    return S0, S1, S2


# ---------------------------------------------------------------------------
# Moran's I, vectorised over columns
# ---------------------------------------------------------------------------
def moran_I(W: sp.csr_matrix, S0: float, Y: np.ndarray) -> np.ndarray:
    """I for every column of Y (n x p, float64).  Constant columns -> NaN."""
    n = Y.shape[0]
    Z = Y - Y.mean(0, keepdims=True)
    den = (Z * Z).sum(0)
    num = (Z * (W @ Z)).sum(0)
    out = np.full(Y.shape[1], np.nan)
    ok = den > 0
    out[ok] = (n / S0) * num[ok] / den[ok]
    return out


def moran_inference(W, S0, S1, S2, Y, I):
    """Analytic normality and randomisation z / two-sided p, per column."""
    n = Y.shape[0]
    EI = -1.0 / (n - 1)
    Z = Y - Y.mean(0, keepdims=True)
    m2 = (Z ** 2).sum(0) / n
    m4 = (Z ** 4).sum(0) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        b2 = m4 / (m2 ** 2)
    V_norm = (n ** 2 * S1 - n * S2 + 3 * S0 ** 2) / (S0 ** 2 * (n ** 2 - 1)) - EI ** 2
    A = n * ((n ** 2 - 3 * n + 3) * S1 - n * S2 + 3 * S0 ** 2)
    B = b2 * ((n ** 2 - n) * S1 - 2 * n * S2 + 6 * S0 ** 2)
    C = (n - 1) * (n - 2) * (n - 3) * S0 ** 2
    V_rand = (A - B) / C - EI ** 2
    from scipy.stats import norm
    with np.errstate(invalid="ignore"):
        z_norm = (I - EI) / np.sqrt(V_norm)
        z_rand = (I - EI) / np.sqrt(np.maximum(V_rand, 0.0))
    p_norm = 2 * norm.sf(np.abs(z_norm))
    p_rand = 2 * norm.sf(np.abs(z_rand))
    return EI, z_norm, p_norm, z_rand, p_rand


def moran_perm(W, S0, Y, I, n_perm, seed):
    """Conditional-permutation two-sided p, per column of Y."""
    rng = np.random.default_rng(seed)
    n, p = Y.shape
    Z = Y - Y.mean(0, keepdims=True)
    den = (Z * Z).sum(0)
    ge = np.zeros(p, dtype=np.int64)
    obs = np.abs(I - (-1.0 / (n - 1)))
    for _ in range(n_perm):
        Zp = np.empty_like(Z)
        for j in range(p):
            Zp[:, j] = Z[rng.permutation(n), j]
        Ip = (n / S0) * (Zp * (W @ Zp)).sum(0) / den
        ge += (np.abs(Ip - (-1.0 / (n - 1))) >= obs).astype(np.int64)
    return (ge + 1.0) / (n_perm + 1.0)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_h5_all(sample: str):
    """cells x features CSC->CSR for ALL feature types, plus names and types."""
    path = os.path.join(RAW, sample, "cell_feature_matrix.h5")
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        names = np.array([x.decode() for x in g["features"]["name"][:]])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        n_feat, n_cells = g["shape"][:]
        M = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=(int(n_feat), int(n_cells)))
        barcodes = np.array([x.decode() for x in g["barcodes"][:]])
    return M.T.tocsr(), names, ftype, barcodes


def centre_by_celltype(v: np.ndarray, ct: np.ndarray) -> np.ndarray:
    out = v.astype(float).copy()
    for c in np.unique(ct):
        m = ct == c
        out[m] -= out[m].mean()
    return out


def build_fields(sec, Xc, ftype, cells):
    """The A7 responses, the seven modules and the technical references."""
    cls_sum = {}
    for resp, cl in CLASS_OF.items():
        cls_sum[resp] = np.asarray(Xc[:, ftype == cl].sum(1)).ravel().astype(float)
    cls_sum["all_controls"] = sum(cls_sum[r] for r in CLASS_OF)
    tc = cells["transcript_counts"].to_numpy(float)
    cls_sum["neg_probe_rate"] = cls_sum["neg_control_probe"] / (tc + 1.0)

    fields, kind = {}, {}
    for k, v in cls_sum.items():
        fields[k] = v
        kind[k] = "control"
    for m in P.MODULES:
        fields[m] = sec.module(m)
        kind[m] = "module"
    fields["transcript_counts"] = tc
    fields["cell_area"] = cells["cell_area"].to_numpy(float)
    fields["density_50um"] = sec.density_50um.astype(float)
    fields["unassigned_codeword"] = np.asarray(
        Xc[:, ftype == "Unassigned Codeword"].sum(1)).ravel().astype(float)
    for k in ("transcript_counts", "cell_area", "density_50um",
              "unassigned_codeword"):
        kind[k] = "technical"
    return fields, kind


# ---------------------------------------------------------------------------
# one section
# ---------------------------------------------------------------------------
def run_section(sample: str, n_perm: int, seed: int):
    t0 = time.time()
    sec = P.Sec(sample)
    coords = sec.coords.astype(float)
    cid = sec.cell_id.astype(str)
    ct = sec.celltype
    n = coords.shape[0]

    X, names, ftype, barcodes = load_h5_all(sample)
    pos = pd.Index(barcodes).get_indexer(pd.Index(cid))
    assert (pos >= 0).all(), "cache cell_id missing from h5 barcodes"
    Xc = X[pos].tocsc()
    del X

    cells = (pd.read_parquet(os.path.join(RAW, sample, "cells.parquet"))
             .set_index("cell_id").reindex(pd.Index(cid)))
    assert cells["transcript_counts"].notna().all()

    # provenance: the h5 class sums must equal the cells.parquet tallies
    prov_ok = {}
    for resp, cl in CLASS_OF.items():
        s = np.asarray(Xc[:, ftype == cl].sum(1)).ravel()
        col = {"neg_control_probe": "control_probe_counts",
               "neg_control_codeword": "control_codeword_counts",
               "genomic_control": "genomic_control_counts"}[resp]
        prov_ok[resp] = bool(np.array_equal(s, cells[col].to_numpy()))
    assert all(prov_ok.values()), (sample, prov_ok)

    fields, kind = build_fields(sec, Xc, ftype, cells)
    fnames = list(fields)
    Yraw = np.column_stack([fields[k] for k in fnames]).astype(float)
    Ycc = np.column_stack([centre_by_celltype(fields[k], ct) for k in fnames])

    rows_agg, rows_feat = [], []

    # ---- weights variants ------------------------------------------------
    variants = [("knn%d" % KNN_PRIMARY, KNN_PRIMARY, None)]
    variants += [("knn%d" % k, k, None) for k in KNN_SENS]
    variants += [("band%dum" % int(r), None, r) for r in BANDS_UM]

    for label, k, r in variants:
        if k is not None:
            W = knn_weights(coords, k)
            islands = 0
            mean_nb = float(k)
        else:
            W, islands = band_weights(coords, r)
            mean_nb = float(np.asarray((W > 0).sum(1)).ravel().mean())
        S0, S1, S2 = w_moments(W)
        primary = (label == "knn%d" % KNN_PRIMARY)

        for cset, Y in (("raw", Yraw), ("ctcentred", Ycc)):
            I = moran_I(W, S0, Y)
            EI, zn, pn, zr, pr = moran_inference(W, S0, S1, S2, Y, I)
            psim = np.full(len(fnames), np.nan)
            if primary and cset == "raw" and n_perm > 0:
                psim = moran_perm(W, S0, Y, I, n_perm, seed)
            for j, nm in enumerate(fnames):
                x = Y[:, j]
                rows_agg.append(dict(
                    section=sample, arm=sec.meta["condition"],
                    band=sec.meta["band"], n_cells=n,
                    weights=label, mean_neighbours=mean_nb, islands=islands,
                    centring=cset, field=nm, kind=kind[nm],
                    moran_I=I[j], E_I=EI, z_rand=zr[j], p_rand=pr[j],
                    z_norm=zn[j], p_norm=pn[j], p_sim=psim[j],
                    mean=float(x.mean()), sd=float(x.std()),
                    frac_nonzero=float(np.mean(Yraw[:, j] > 0))))

        # ---- per-feature (Voyager reproduction), primary weights only ----
        if primary:
            step = 256
            tot = np.asarray(Xc.sum(0)).ravel()
            nzc = np.asarray((Xc > 0).sum(0)).ravel()
            for a in range(0, Xc.shape[1], step):
                b = min(a + step, Xc.shape[1])
                Yb = np.asarray(Xc[:, a:b].todense(), dtype=float)
                Ib = moran_I(W, S0, Yb)
                _, _, _, zrb, prb = moran_inference(W, S0, S1, S2, Yb, Ib)
                for j in range(b - a):
                    rows_feat.append(dict(
                        section=sample, arm=sec.meta["condition"],
                        weights=label, feature=names[a + j],
                        feature_type=ftype[a + j],
                        moran_I=Ib[j], z_rand=zrb[j], p_rand=prb[j],
                        total_counts=float(tot[a + j]),
                        frac_cells_nonzero=float(nzc[a + j]) / n))
                del Yb
        del W
    print("[moran] %s n=%d %.0fs" % (sample, n, time.time() - t0), flush=True)
    return pd.DataFrame(rows_agg), pd.DataFrame(rows_feat)


# ---------------------------------------------------------------------------
# validation against esda / libpysal
# ---------------------------------------------------------------------------
def validate(sample="7250_liver_sham_Male_26-U1", n_sub=4000, seed=7):
    import libpysal
    from esda.moran import Moran
    sec = P.Sec(sample)
    rng = np.random.default_rng(seed)
    sub = rng.choice(sec.n, size=n_sub, replace=False)
    co = sec.coords.astype(float)[sub]
    cells = (pd.read_parquet(os.path.join(RAW, sample, "cells.parquet"))
             .set_index("cell_id").reindex(pd.Index(sec.cell_id.astype(str))))
    ys = {"control_probe_counts": cells["control_probe_counts"].to_numpy(float)[sub],
          "transcript_counts": cells["transcript_counts"].to_numpy(float)[sub],
          "emt_ecm": sec.module("emt_ecm")[sub]}
    w = libpysal.weights.KNN.from_array(co, k=KNN_PRIMARY)
    w.transform = "r"
    W = knn_weights(co, KNN_PRIMARY)
    S0, S1, S2 = w_moments(W)
    print("libpysal S0=%.4f mine S0=%.4f  S1 %.4f/%.4f  S2 %.4f/%.4f"
          % (w.s0, S0, w.s1, S1, w.s2, S2))
    ok = True
    for nm, y in ys.items():
        m = Moran(y, w, permutations=999)
        Y = y[:, None]
        I = moran_I(W, S0, Y)[0]
        _, zn, pn, zr, pr = moran_inference(W, S0, S1, S2, Y, np.array([I]))
        d_I, d_z = abs(I - m.I), abs(zr[0] - m.z_rand)
        ok &= (d_I < 1e-10) and (d_z < 1e-8)
        print("  %-22s esda I=%+.6f mine %+.6f (d=%.2e) | z_rand esda %+.4f "
              "mine %+.4f (d=%.2e) | p_rand esda %.4g mine %.4g | esda p_sim %.4g"
              % (nm, m.I, I, d_I, m.z_rand, zr[0], d_z, m.p_rand, pr[0], m.p_sim))
    print("VALIDATION", "PASS" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", default="all")
    ap.add_argument("--n-jobs", type=int, default=3)
    ap.add_argument("--perms", type=int, default=999)
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()
    if a.validate:
        sys.exit(0 if validate() else 1)
    os.makedirs(OUT, exist_ok=True)
    secs = P.ALL_SECTIONS if a.sections == "all" else P.IN_BAND
    res = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(run_section)(s, a.perms, P.MASTER_SEED + 13 * i)
        for i, s in enumerate(secs))
    agg = pd.concat([r[0] for r in res], ignore_index=True)
    feat = pd.concat([r[1] for r in res], ignore_index=True)
    agg.to_csv(os.path.join(OUT, "moran_fields.csv"), index=False)
    feat.to_csv(os.path.join(OUT, "moran_per_feature.csv.gz"), index=False)
    print(agg.shape, "->", OUT + "/moran_fields.csv")
    print(feat.shape, "->", OUT + "/moran_per_feature.csv.gz")


if __name__ == "__main__":
    main()
