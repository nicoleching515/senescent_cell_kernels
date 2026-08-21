"""Phase 4 — data layer for the existing-tools benchmark (Figure 4).

Loads, for each admissible section, the analysis cell set used by the Phase 3
null battery (cache3) plus log-normalised expression of the nine Tier C genes
Bio certified as defensible on this panel (BIO_PHASE3 section 3.4).

Nothing here builds an (n,n) matrix; all neighbour work is cKDTree.
"""
from __future__ import annotations
import os, numpy as np, pandas as pd
from scipy.spatial import cKDTree
import sasp_real as R
import sasp_phase3 as P

P4 = "/workspace/results/phase4"
CACHE4 = "/tmp/p4/cache4"
MASTER_SEED = 20260820

SECTIONS = list(P.IN_BAND)              # the six Section-8 Test-3 admissible sections
EXCLUDE_TYPES = P.EXCLUDE_TYPES

# BIO_PHASE3 section 3.4: the only Tier C pairs defensible on this panel.
# Ackr3 <=1% of cells, Il6 ~3 cells/section, Cxcl2 ~8, mouse has no CXCL8.
LR_PAIRS = [
    dict(pair="Ccl2->Ccr2",              pathway="CCL2",  ligand="Ccl2",
         receptors=("Ccr2",)),
    dict(pair="Tnf->Tnfrsf1a/1b",        pathway="TNF",   ligand="Tnf",
         receptors=("Tnfrsf1a", "Tnfrsf1b")),
    dict(pair="Tgfb1->Tgfbr1/2",         pathway="TGFB",  ligand="Tgfb1",
         receptors=("Tgfbr1", "Tgfbr2")),
    dict(pair="Il1a->Il1r1",             pathway="IL1",   ligand="Il1a",
         receptors=("Il1r1",)),
]
LR_GENES = sorted({g for d in LR_PAIRS for g in (d["ligand"],) + d["receptors"]})

WINDOW_UM = 100.0     # CS_PHASE3 section 2: the fitting window adopted project-wide


def build(sample: str, force: bool = False) -> str:
    os.makedirs(CACHE4, exist_ok=True)
    out = os.path.join(CACHE4, f"{sample}.npz")
    if os.path.exists(out) and not force:
        return f"[skip] {sample}"
    sec = P.Sec(sample)                       # merged labels, Phase 3 cell set
    cid = sec.cell_id.astype(str)
    M, names, bc = R.load_expression(sample)
    bc = bc.astype(str)
    pos = pd.Index(bc).get_indexer(cid)
    assert (pos >= 0).all(), "cache3 cell not in h5"
    gi = pd.Index(names).get_indexer(LR_GENES)
    assert (gi >= 0).all(), [g for g, j in zip(LR_GENES, gi) if j < 0]
    counts = np.asarray(M[pos][:, gi].todense(), dtype=np.float32)
    tot = np.asarray(M[pos].sum(1)).ravel().astype(np.float32)
    med = float(np.median(tot))
    ln = np.log1p(counts / np.maximum(tot, 1)[:, None] * med)   # scanpy-style
    np.savez_compressed(out, cell_id=cid, coords=sec.coords.astype(np.float32),
                        celltype=sec.celltype.astype(str),
                        counts=counts, lognorm=ln.astype(np.float32),
                        total_counts=tot, genes=np.array(LR_GENES))
    return f"[ok] {sample} n={len(cid)}"


class Sec4:
    """Analysis view of one section: labelled cells only, LR expression, coords."""

    def __init__(self, sample: str):
        z = np.load(os.path.join(CACHE4, f"{sample}.npz"), allow_pickle=False)
        keep = ~np.isin(z["celltype"], EXCLUDE_TYPES)
        self.name = sample
        self.meta = R.parse_sample(sample)
        self.coords = z["coords"][keep].astype(float)
        self.celltype = z["celltype"][keep]
        self.genes = list(z["genes"])
        self.X = z["lognorm"][keep].astype(float)      # cells x 9 genes
        self.counts = z["counts"][keep]
        self.cell_id = z["cell_id"][keep]
        self.types = np.array(sorted(set(self.celltype)))
        self.tcode = pd.Index(self.types).get_indexer(self.celltype)
        self.n = self.coords.shape[0]
        self.lo = self.coords.min(0)
        self.hi = self.coords.max(0)
        self.cen = self.coords.mean(0)
        d, _ = cKDTree(self.coords).query(self.coords, k=2, workers=-1)
        self.med_nn = float(np.median(d[:, 1]))

    def g(self, name):
        return self.X[:, self.genes.index(name)]

    def gc(self, name):
        return self.counts[:, self.genes.index(name)]

    def receptor(self, receptors):
        """Alternative receptors combine as an OR: a cell responds through
        whichever it carries.  Multi-subunit complexes would be a geometric
        mean; Tnfrsf1a/1b and Tgfbr1/2 are used here as alternatives."""
        return np.max(np.column_stack([self.g(r) for r in receptors]), axis=1)


# --------------------------------------------------------------------------
# coordinate nulls -- torus_shift / rotate are copied verbatim from
# run_phase3_nulls.py so Figure 4 is on the same footing as Figure 2c.
# --------------------------------------------------------------------------

def torus_shift(rng, pts, lo, hi):
    span = hi - lo
    return lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span


def rotate_about_centroid(rng, pts, cen, lo, hi):
    th = rng.uniform(0, 2 * np.pi)
    Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    span = hi - lo
    return lo + ((pts - cen) @ Rm.T + cen - lo) % span


NULLS = ("N3_lig", "N4_lig", "N3_type", "N0_perm")
NULL_LABEL = {
    "real":    "real coordinates",
    "N3_lig":  "N3  torus shift of ligand+ cells",
    "N4_lig":  "N4  rotation of ligand+ cells",
    "N3_type": "N3t per-cell-type torus shift",
    "N0_perm": "N0  full coordinate permutation",
    "N0_type": "N0t within-cell-type location permutation (CellWHISPER's own)",
}


def null_coords(sec: Sec4, null: str, rng, lig_mask=None):
    """Return a new (n,2) coordinate array under the requested null."""
    xy = sec.coords
    if null == "real":
        return xy
    if null == "N0_perm":
        # NOT CellWHISPER's control -- see the N0_type branch below and D7 B9.
        # This permutes across ALL cells and so also destroys the cell-type
        # spatial architecture; it is strictly more destructive than theirs.
        return xy[rng.permutation(sec.n)]
    out = xy.copy()
    if null in ("N3_lig", "N4_lig"):
        m = lig_mask
        if m is None or m.sum() < 5:
            return None
        pts = xy[m]
        out[m] = (torus_shift(rng, pts, sec.lo, sec.hi) if null == "N3_lig"
                  else rotate_about_centroid(rng, pts, sec.cen, sec.lo, sec.hi))
        return out
    if null == "N0_type":
        # CellWHISPER's ACTUAL randomised control (D7 sec 2.1 B9): permute cell
        # locations WITHIN each cell type, "preserving cell-type-specific spatial
        # organization and ligand-receptor (LR) expression while destroying spatial
        # proximity between ligand- and receptor-expressing cells".  Strictly less
        # destructive than N0_perm, which permutes across all cells and so also
        # destroys the cell-type architecture.
        code = getattr(sec, "tcode", None)
        if code is None:
            code = sec.code
        for t in range(len(sec.types)):
            idx = np.flatnonzero(code == t)
            if idx.size:
                out[idx] = xy[rng.permutation(idx)]
        return out
    if null == "N3_type":                       # independent shift per cell type
        code = getattr(sec, "tcode", None)
        if code is None:
            code = sec.code
        for t in range(len(sec.types)):
            m = code == t
            out[m] = torus_shift(rng, xy[m], sec.lo, sec.hi)
        return out
    raise ValueError(null)


if __name__ == "__main__":
    for s in SECTIONS:
        print(build(s), flush=True)
