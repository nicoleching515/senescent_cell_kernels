"""
Real-tissue loader and geometry for the SASP spatial response kernel
(Master Plan Section 22 Step 2).

Data: 10x Xenium Prime Mouse 5K liver sections. Coordinates are already in
MICRONS (`x_centroid`, `y_centroid` in cells.parquet) -- unit is in every
column name, per Section 8 Test 1.

Design rule for this module: **sender calls, cell-type labels, anatomy and
module scores are SWAPPABLE INPUTS, never hardcoded.** The Bio agent is
producing `celltypes_/anatomy_/senders_/modules_{sample}.csv` in
/workspace/data/processed/ concurrently; this module uses them when present and
falls back to a clearly-labelled provisional call otherwise, so the pipeline can
be built and validated without blocking on them.

Compute rule (Section 18.1): cKDTree for everything. A 238k-cell section is
2.8e10 pairs; no (n, n) matrix is built anywhere.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

RAW = "/workspace/data/raw"
PROC = "/workspace/data/processed"
GENESETS = "/workspace/genesets"

# Feature types to KEEP. Everything else -- Negative Control Probe/Codeword,
# Genomic Control, Unassigned Codeword, Deprecated Codeword -- is dropped.
KEEP_FEATURE_TYPE = "Gene Expression"

# Xenium Prime Mouse 5K carries a handful of genotyping probes among the
# Gene Expression features; drop by name pattern if present.
GENOTYPING_RE = re.compile(r"(^|[_-])(genotyp|gt)([_-]|$)", re.I)


def list_samples() -> List[str]:
    return sorted(d for d in os.listdir(RAW)
                  if os.path.isdir(os.path.join(RAW, d))
                  and os.path.exists(os.path.join(RAW, d, "cells.parquet")))


def parse_sample(sample: str) -> Dict[str, object]:
    """`7250_liver_sham_Male_26-U1` -> mouse/tissue/condition/sex/week/section."""
    m = re.match(r"(\d+)_(\w+?)_(\w+?)_(\w+?)_(\d+)-(\w+)$", sample)
    if not m:
        return dict(sample=sample, mouse=sample, tissue="?", condition="?",
                    sex="?", week=np.nan, section=sample)
    mouse, tissue, cond, sex, week, sec = m.groups()
    return dict(sample=sample, mouse=mouse, tissue=tissue, condition=cond,
                sex=sex, week=int(week), section=f"{mouse}-{sec}")


# --------------------------------------------------------------------------
# expression matrix
# --------------------------------------------------------------------------


def load_expression(sample: str):
    """Return (csr cells x genes, gene_names, barcodes) for Gene Expression
    features only."""
    path = os.path.join(RAW, sample, "cell_feature_matrix.h5")
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        names = np.array([x.decode() for x in g["features"]["name"][:]])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        n_genes, n_cells = g["shape"][:]
        # stored CSC with cells as columns
        M = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=(int(n_genes), int(n_cells)))
        barcodes = np.array([x.decode() for x in g["barcodes"][:]])

    keep = (ftype == KEEP_FEATURE_TYPE) & ~np.array(
        [bool(GENOTYPING_RE.search(n)) for n in names])
    M = M[keep, :]
    names = names[keep]
    return M.T.tocsr(), names, barcodes          # cells x genes


def read_geneset(name: str) -> List[str]:
    path = name if os.path.sep in name else os.path.join(GENESETS, name)
    if not path.endswith(".txt"):
        path += ".txt"
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------


def normalize_counts(M: sp.csr_matrix, target: float = 1e4):
    """CP10K + log1p. Returns a dense-per-gene-friendly csr."""
    tot = np.asarray(M.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    inv = sp.diags(target / tot)
    L = inv @ M
    L.data = np.log1p(L.data)
    return L


def module_score(L: sp.csr_matrix, gene_names: np.ndarray,
                 genes: Sequence[str], zcap: float = 10.0):
    """Mean of per-gene z-scores over the module.

    Simple, transparent, and swappable -- if the Bio agent's `modules_*.csv`
    is present it is used instead (see `load_sample`). Genes absent from the
    panel are dropped and reported.
    """
    idx = {g: i for i, g in enumerate(gene_names)}
    cols = [idx[g] for g in genes if g in idx]
    missing = [g for g in genes if g not in idx]
    if not cols:
        return np.full(L.shape[0], np.nan), [], missing
    sub = L[:, cols].toarray()
    mu = sub.mean(axis=0)
    sd = sub.std(axis=0)
    sd[sd < 1e-9] = 1.0
    z = np.clip((sub - mu) / sd, -zcap, zcap)
    return z.mean(axis=1), [gene_names[c] for c in cols], missing


# --------------------------------------------------------------------------
# geometry (cKDTree only)
# --------------------------------------------------------------------------


def geometry(coords: np.ndarray, radii=(25.0, 50.0, 100.0), knn: int = 20):
    """Local density at several radii, kNN distances, median NN distance."""
    tree = cKDTree(coords)
    out = {}
    for r in radii:
        out[f"density_{int(r)}um"] = np.asarray(
            tree.query_ball_point(coords, r=r, return_length=True,
                                  workers=-1)) - 1.0
    nn, nn_idx = tree.query(coords, k=min(knn + 1, coords.shape[0]), workers=-1)
    out["nn1_um"] = nn[:, 1]
    out["median_nn_um"] = float(np.median(nn[:, 1]))
    return out, tree, nn_idx


def distance_to_set(coords: np.ndarray, member_mask: np.ndarray) -> np.ndarray:
    """Distance from every cell to the nearest member of `member_mask`."""
    if member_mask.sum() == 0:
        return np.full(coords.shape[0], np.nan)
    d, _ = cKDTree(coords[member_mask]).query(coords, k=1, workers=-1)
    return d


def knn_composition(labels: np.ndarray, knn_idx: np.ndarray,
                    categories: Sequence) -> np.ndarray:
    """Fraction of each label among the k nearest neighbours (self excluded)."""
    nb = labels[knn_idx[:, 1:]]
    return np.stack([(nb == c).mean(axis=1) for c in categories], axis=1)


# --------------------------------------------------------------------------
# the sample bundle
# --------------------------------------------------------------------------


@dataclass
class Sample:
    name: str
    meta: Dict[str, object]
    cells: pd.DataFrame           # one row per cell, coords in um
    coords: np.ndarray
    gene_names: np.ndarray
    modules: Dict[str, np.ndarray]
    sender: Dict[str, np.ndarray]        # sender-call name -> boolean mask
    celltype: np.ndarray
    celltype_source: str
    sender_source: str
    module_source: str
    knn_idx: np.ndarray
    median_nn_um: float

    def n(self) -> int:
        return self.coords.shape[0]


def _proc(sample: str, kind: str) -> Optional[pd.DataFrame]:
    p = os.path.join(PROC, f"{kind}_{sample}.csv")
    if os.path.exists(p):
        try:
            return pd.read_csv(p)
        except Exception:
            return None
    return None


def load_sample(sample: str,
                module_sets: Optional[Dict[str, str]] = None,
                min_counts: int = 20,
                provisional_sender_gene: str = "Cdkn1a") -> Sample:
    """Load one section with geometry, module scores, sender calls and cell
    types. Bio-agent annotations are used when present; otherwise a clearly
    labelled provisional call is made so the pipeline is not blocked.
    """
    meta = parse_sample(sample)
    cells = pd.read_parquet(os.path.join(RAW, sample, "cells.parquet"))
    M, gene_names, barcodes = load_expression(sample)
    assert M.shape[0] == cells.shape[0]
    if not np.array_equal(barcodes, cells["cell_id"].to_numpy().astype(str)):
        order = pd.Series(np.arange(len(barcodes)), index=barcodes)
        cells = cells.set_index("cell_id").loc[barcodes].reset_index()

    # QC: drop near-empty cells (they carry no usable module score)
    keep = cells["transcript_counts"].to_numpy() >= min_counts
    cells = cells.loc[keep].reset_index(drop=True)
    M = M[keep]

    coords = np.column_stack([cells["x_centroid"].to_numpy(float),
                              cells["y_centroid"].to_numpy(float)])
    geo, tree, knn_idx = geometry(coords)
    med_nn = geo.pop("median_nn_um")
    for k, v in geo.items():
        cells[k] = v

    L = normalize_counts(M)

    # ---- modules: Bio's file if present, else score here ------------------
    modules: Dict[str, np.ndarray] = {}
    mod_df = _proc(sample, "modules")
    module_source = "bio:modules_*.csv"
    if mod_df is not None and "cell_id" in mod_df.columns:
        mod_df = mod_df.set_index("cell_id").reindex(cells["cell_id"])
        for c in mod_df.columns:
            modules[c] = mod_df[c].to_numpy(float)
    else:
        module_source = "provisional:z-mean of Tier B sets (Bio pending)"
        if module_sets is None:
            module_sets = {f[2:-4]: f for f in sorted(os.listdir(GENESETS))
                           if f.startswith("B_") and f.endswith(".txt")}
        for mod, fname in module_sets.items():
            s, used, miss = module_score(L, gene_names, read_geneset(fname))
            modules[mod] = s
            modules.setdefault("_panel_n", {})
        modules = {k: v for k, v in modules.items() if not k.startswith("_")}

    # ---- senders: Bio's file if present, else provisional ------------------
    sender: Dict[str, np.ndarray] = {}
    sen_df = _proc(sample, "senders")
    sender_source = "bio:senders_*.csv"
    if sen_df is not None and "cell_id" in sen_df.columns:
        sen_df = sen_df.set_index("cell_id").reindex(cells["cell_id"])
        for c in sen_df.columns:
            v = sen_df[c]
            if v.dropna().isin([0, 1, True, False]).all():
                sender[c] = v.fillna(0).to_numpy().astype(bool)
    if not sender:
        sender_source = (f"provisional:{provisional_sender_gene}>0 "
                         "(Bio senders_*.csv pending)")
        gi = np.flatnonzero(gene_names == provisional_sender_gene)
        pos = (np.asarray(M[:, gi[0]].todense()).ravel() > 0 if gi.size
               else np.zeros(cells.shape[0], bool))
        sender[f"{provisional_sender_gene}_pos"] = pos
        # a Tier-A score threshold variant, for the N7 sensitivity axis
        a_score, _, _ = module_score(L, gene_names,
                                     read_geneset("A_SENDER_FINAL_strict.txt"))
        cells["tierA_score"] = a_score
        for q in (90, 95, 99):
            sender[f"tierA_p{q}"] = a_score >= np.percentile(a_score, q)

    # ---- cell types: Bio's file if present --------------------------------
    ct_df = _proc(sample, "celltypes")
    if ct_df is not None and "cell_id" in ct_df.columns:
        col = [c for c in ct_df.columns if c != "cell_id"][0]
        ct = (ct_df.set_index("cell_id").reindex(cells["cell_id"])[col]
              .fillna("unknown").to_numpy().astype(str))
        celltype_source = f"bio:celltypes_*.csv[{col}]"
    else:
        ct = np.full(cells.shape[0], "all_cells")
        celltype_source = "provisional:unstratified (Bio celltypes_*.csv pending)"

    an_df = _proc(sample, "anatomy")
    if an_df is not None and "cell_id" in an_df.columns:
        an_df = an_df.set_index("cell_id").reindex(cells["cell_id"])
        for c in an_df.columns:
            cells[f"anat_{c}"] = an_df[c].to_numpy()

    return Sample(name=sample, meta=meta, cells=cells, coords=coords,
                  gene_names=gene_names, modules=modules, sender=sender,
                  celltype=ct, celltype_source=celltype_source,
                  sender_source=sender_source, module_source=module_source,
                  knn_idx=knn_idx, median_nn_um=med_nn)
