"""Phase 9 — shared loaders for the H1 human spleen arm (GEO GSE326743).

New file. Nothing under `data/processed/`, `results/phase3/`, `figures/` or `genesets/`
is written by anything that imports this module: H1 outputs go to
  data/processed_h1/     per-cell tables and caches
  results/phase9_h1/     audit tables and summaries

The M1 mouse loaders (`sasp_io.py`, `sasp_real.py`, `sasp_phase3.py`) are left untouched;
this is the human analogue of `sasp_real.load_expression` + `sasp_phase3.prep`.
"""
from __future__ import annotations

import gzip, io, os
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix, csr_matrix

RAW = "/workspace/data/raw_h1"
PROC = "/workspace/data/processed_h1"
RESULTS = "/workspace/results/phase9_h1"
GS_HUMAN = "/workspace/genesets/human"
MASTER_SEED = 20260820          # sasp_phase3.MASTER_SEED, unchanged

# section -> (GSM, age, sex).  Source: reports/PHASE7_H1_SCREEN.md §3.
SECTIONS = {
    "SPLN07": ("GSM9638040", 17, "F"),
    "SPLN14": ("GSM9638041", 37, "F"),
    "SPLN21": ("GSM9638042", 32, "M"),
    "SPLN24": ("GSM9638043", 32, "F"),
    "SPLN30": ("GSM9638044", 57, "M"),
    "SPLN43": ("GSM9638045", 31, "M"),
    "SPLN44": ("GSM9638046", 59, "M"),
}
ALL_SECTIONS = list(SECTIONS)

# frozen QC / annotation constants, copied verbatim from code/annotate_pipeline.py
MIN_COUNTS, MIN_GENES = 20, 5

# frozen exclusions, copied verbatim from code/sasp_phase3.py
EXCLUDE_TYPES = ("Low_quality", "Unknown", "unknown")
EXCLUDE_FROM_SENDERS = ("Proliferating",)

MODULES = ["downstream_arrest", "emt_ecm", "il6_jak_stat3", "interferon_response",
           "oxidative_stress", "secondary_senescence", "tnfa_nfkb_proximal"]


def path(section: str, kind: str) -> str:
    gsm = SECTIONS[section][0]
    return os.path.join(RAW, section, f"{gsm}_{section}_{kind}")


def read_parquet_gz(p: str) -> pd.DataFrame:
    """GEO ships these gzipped; pyarrow cannot read a .parquet.gz directly."""
    with gzip.open(p, "rb") as fh:
        return pd.read_parquet(io.BytesIO(fh.read()))


def cells_table(section: str) -> pd.DataFrame:
    return read_parquet_gz(path(section, "cells.parquet.gz"))


def annotations(section: str) -> pd.DataFrame:
    return pd.read_csv(path(section, "annotations.csv.gz"))


def load_features(section: str):
    with h5py.File(path(section, "cell_feature_matrix.h5"), "r") as f:
        ft = np.array([x.decode() for x in f["matrix/features/feature_type"][:]])
        ids = np.array([x.decode() for x in f["matrix/features/id"][:]])
        nm = np.array([x.decode() for x in f["matrix/features/name"][:]])
    return ft, ids, nm


def load_matrix(section: str, which: str = "gene"):
    """Return (csr cells x features, feature_names, barcodes).

    which='gene'    -> the 5,093 ENSG Gene Expression features (the panel)
    which='control' -> the 40 negative control probes + 609 codewords + 21 genomic controls
    which='all'     -> everything in the file
    """
    p = path(section, "cell_feature_matrix.h5")
    with h5py.File(p, "r") as f:
        ft = np.array([x.decode() for x in f["matrix/features/feature_type"][:]])
        ids = np.array([x.decode() for x in f["matrix/features/id"][:]])
        nm = np.array([x.decode() for x in f["matrix/features/name"][:]])
        bc = np.array([x.decode() for x in f["matrix/barcodes"][:]])
        if which == "gene":
            keep = (ft == "Gene Expression") & np.char.startswith(ids, "ENSG")
        elif which == "control":
            keep = np.isin(ft, ["Negative Control Probe", "Negative Control Codeword",
                                "Genomic Control"])
        else:
            keep = np.ones(ft.size, bool)
        M = csc_matrix((f["matrix/data"][:].astype(np.float32),
                        f["matrix/indices"][:].astype(np.int32),
                        f["matrix/indptr"][:].astype(np.int64)),
                       shape=tuple(f["matrix/shape"][:]))
    X = csr_matrix(M[np.where(keep)[0], :].T)      # cells x features
    return X, nm[keep], bc, ft[keep]


def gl(name: str):
    """Read a frozen human gene-set file."""
    return [l.strip() for l in open(os.path.join(GS_HUMAN, name + ".txt")) if l.strip()]


def marker_set():
    """The frozen 22-label spleen marker set, plus the pre-registered plasma-cell
    exception (PREREG_PHASE8.md P6 / PI decision D2b: MIN_MARKERS=4 is waived for
    'Plasma cells', whose three surviving on-panel markers are JCHAIN/MZB1/XBP1).

    `code/markers_human_spleen.py` is a frozen generated file under
    genesets/.geneset_manifest.json and is NOT edited; the exception is applied here.
    """
    import sys
    sys.path.insert(0, "/workspace/code")
    from markers_human_spleen import MARKERS, MERGE
    M = dict(MARKERS)
    M["Plasma cells"] = ["JCHAIN", "MZB1", "XBP1"]     # P6 exception, 3 markers
    return M, dict(MERGE), {"Plasma cells": 3}
