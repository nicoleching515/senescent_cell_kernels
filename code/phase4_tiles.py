"""Phase 4 — spatially contiguous tiles.

COMMOT v0.0.3 materialises a DENSE (n,n) distance matrix
(`tools/_spatial_communication.py:390`), which is 105 GB for the 114,721-cell
section 7259 and kills the process.  That is a property of the published
software, not of our data, and it is reported as such.  So every method is run
on the same set of spatially CONTIGUOUS tiles at native cell density -- not on
random subsamples, which would change local density and therefore change what a
distance-thresholded method sees.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import phase4_data as D

TILE_UM = 1200.0
N_TILES = 3


def tiles_for(sec: D.Sec4, n_tiles=N_TILES, side=TILE_UM):
    xy = sec.coords
    lo, hi = sec.lo, sec.hi
    step = side / 2.0
    gx = np.arange(lo[0], hi[0] - side + 1e-9, step)
    gy = np.arange(lo[1], hi[1] - side + 1e-9, step)
    cand = []
    for x0 in gx:
        for y0 in gy:
            m = ((xy[:, 0] >= x0) & (xy[:, 0] < x0 + side) &
                 (xy[:, 1] >= y0) & (xy[:, 1] < y0 + side))
            cand.append((int(m.sum()), x0, y0))
    cand.sort(reverse=True)
    chosen = []
    for n, x0, y0 in cand:
        if any(abs(x0 - a) < side and abs(y0 - b) < side for _, a, b in chosen):
            continue
        chosen.append((n, x0, y0))
        if len(chosen) == n_tiles:
            break
    return chosen


class Tile:
    def __init__(self, sec: D.Sec4, x0, y0, ti, side=TILE_UM):
        xy = sec.coords
        m = ((xy[:, 0] >= x0) & (xy[:, 0] < x0 + side) &
             (xy[:, 1] >= y0) & (xy[:, 1] < y0 + side))
        self.idx = np.flatnonzero(m)
        self.name = f"{sec.name.split('_')[0]}_t{ti}"
        self.section = sec.name
        self.meta = sec.meta
        self.coords = xy[m]
        self.X = sec.X[m]
        self.counts = sec.counts[m]
        self.genes = sec.genes
        ctype = sec.celltype[m]
        # keep only types with enough cells for a group statistic to exist
        vc = pd.Series(ctype).value_counts()
        self.types = np.array(sorted(vc[vc >= 30].index))
        keep = np.isin(ctype, self.types)
        for a in ("idx", "coords", "X", "counts"):
            setattr(self, a, getattr(self, a)[keep])
        ctype = ctype[keep]
        self.celltype = ctype
        self.K = len(self.types)
        self.code = pd.Index(self.types).get_indexer(ctype)
        self.n = self.coords.shape[0]
        self.lo = self.coords.min(0); self.hi = self.coords.max(0)
        self.cen = self.coords.mean(0)
        d, _ = cKDTree(self.coords).query(self.coords, k=2, workers=1)
        self.med_nn = float(np.median(d[:, 1]))

    def g(self, nm):
        return self.X[:, self.genes.index(nm)]

    def gc(self, nm):
        return self.counts[:, self.genes.index(nm)]

    def receptor(self, rs):
        return np.max(np.column_stack([self.g(r) for r in rs]), axis=1)


def all_tiles():
    out = []
    for s in D.SECTIONS:
        sec = D.Sec4(s)
        for ti, (n, x0, y0) in enumerate(tiles_for(sec)):
            out.append(Tile(sec, x0, y0, ti))
    return out


if __name__ == "__main__":
    rows = []
    for t in all_tiles():
        rows.append(dict(tile=t.name, section=t.section, arm=t.meta["condition"],
                         week=t.meta["week"], n_cells=t.n, n_types=t.K,
                         med_nn_um=round(t.med_nn, 2),
                         **{f"pct_{p['ligand']}": round(100 * float((t.gc(p['ligand']) > 0).mean()), 2)
                            for p in D.LR_PAIRS}))
    df = pd.DataFrame(rows)
    df.to_csv("/workspace/results/phase4/tiles.csv", index=False)
    print(df.to_string(index=False))
