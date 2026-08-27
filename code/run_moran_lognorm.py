#!/usr/bin/env python3
"""Sensitivity for `run_moran_controls.py`: per-feature Moran's I on
CP10K + log1p normalised values instead of raw counts.

Voyager's Xenium vignette computes Moran's I on the SFE's `logcounts` assay.
Raw counts are the conservative choice for a technical-confound diagnostic --
per-cell normalisation divides out total counts, which is the very quantity the
A7 confound runs through -- so the primary run uses raw counts and this script
measures how much the choice matters.

Three sections (one below-floor sham, one in-band sham, one over-ceiling sbr),
primary weights only.  Writes results/moran/moran_per_feature_lognorm.csv.gz.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P                       # noqa: E402
import run_moran_controls as R                # noqa: E402

SECTIONS = ["7250_liver_sham_Male_26-U1",      # below floor
            "7352_liver_sham_Male_2-U1",       # in band
            "7239_liver_sbr_Male_52-U1"]       # over ceiling


def main():
    os.makedirs(R.OUT, exist_ok=True)
    out = []
    for s in SECTIONS:
        t0 = time.time()
        sec = P.Sec(s)
        coords = sec.coords.astype(float)
        cid = sec.cell_id.astype(str)
        n = coords.shape[0]
        X, names, ftype, barcodes = R.load_h5_all(s)
        pos = pd.Index(barcodes).get_indexer(pd.Index(cid))
        assert (pos >= 0).all()
        Xc = X[pos]
        del X
        # CP10K + log1p over ALL features present in the matrix, matching the
        # SFE convention of normalising by the cell's total counts.
        tot = np.asarray(Xc.sum(1)).ravel()
        tot[tot == 0] = 1.0
        Xc = (sp.diags(1e4 / tot) @ Xc).tocsc()
        Xc.data = np.log1p(Xc.data)
        W = R.knn_weights(coords, R.KNN_PRIMARY)
        S0, S1, S2 = R.w_moments(W)
        step = 256
        for a in range(0, Xc.shape[1], step):
            b = min(a + step, Xc.shape[1])
            Yb = np.asarray(Xc[:, a:b].todense(), dtype=float)
            Ib = R.moran_I(W, S0, Yb)
            _, _, _, zr, pr = R.moran_inference(W, S0, S1, S2, Yb, Ib)
            for j in range(b - a):
                out.append(dict(section=s, weights="knn%d" % R.KNN_PRIMARY,
                                normalisation="cp10k_log1p",
                                feature=names[a + j],
                                feature_type=ftype[a + j],
                                moran_I=Ib[j], z_rand=zr[j], p_rand=pr[j]))
            del Yb
        print("[lognorm] %s n=%d %.0fs" % (s, n, time.time() - t0), flush=True)
    D = pd.DataFrame(out)
    D.to_csv(os.path.join(R.OUT, "moran_per_feature_lognorm.csv.gz"),
             index=False)
    print(D.groupby("feature_type").moran_I.describe()[
        ["count", "mean", "50%", "max"]].to_string())


if __name__ == "__main__":
    main()
