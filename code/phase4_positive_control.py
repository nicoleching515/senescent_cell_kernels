"""Positive controls: are these implementations sensitive to anything at all?

The obvious objection to §3 is "your reimplementations are simply insensitive".
Two synthetic controls answer it, on a tile-sized synthetic tissue with a KNOWN
ground truth and no confounding whatsoever.

  C1  CELL-TYPE SENSITIVITY.  Ligand expressed only in cell type A, receptor
      only in type B, background elsewhere.  A -> B is the true interaction and
      nothing else is.  Every method must find it.  This is what the statistic
      is designed to detect, and it verifies the code.

  C2  SPATIAL SENSITIVITY.  Identical cells, identical cell-type labels,
      identical marginal ligand and receptor expression -- the ONLY difference
      between the two conditions is that in one the ligand+ cells sit next to
      the receptor+ cells and in the other they have been torus-shifted away.
      A method that measures spatial communication must separate them.

C1 passing and C2 failing is not a bug: it is the mechanism of §4 demonstrated
on data where the right answer is known.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import phase4_methods as M
import phase4_data as D

OUT = "/workspace/results/phase4/positive_controls.csv"
SEED = 20260820
N = 8000
SIDE = 1200.0
TYPES = np.array(["A", "B", "C", "D"])
K = 4


class Fake:
    """Minimal stand-in for a Tile."""
    def __init__(self, xy, code, lig, rec):
        self.coords, self.code, self.K = xy, code, K
        self.types, self.n = TYPES, xy.shape[0]
        self.lig, self.rec = lig, rec
        self.lo, self.hi = xy.min(0), xy.max(0)
        self.cen = xy.mean(0)
        d, _ = cKDTree(xy).query(xy, k=2, workers=1)
        self.med_nn = float(np.median(d[:, 1]))


def _score(f, rng):
    out = {}
    s, p = M.cellchat_many([(f.lig, f.rec)], f.coords, f.code, K, f.med_nn,
                           rng, n_perm=200)
    out["CellChat v2*"] = (s[0], p[0])
    s, p, _ = M.spatalk_many([(f.lig, f.rec)], f.coords, f.code, K, rng,
                             n_perm=500)
    out["SpaTalk*"] = (s[0], p[0])
    B, P, _ = M.ncem_linear(f.rec, f.coords, f.code, K, 30.0)
    out["NCEM linear*"] = (B, M.bh(P))
    return out


def _commot(f, seed):
    import anndata as ad, commot as ct
    X = np.column_stack([f.lig, f.rec]).astype(np.float64)
    A = ad.AnnData(X, obs=pd.DataFrame({"celltype": TYPES[f.code]},
                                       index=[str(i) for i in range(f.n)]),
                   var=pd.DataFrame(index=["Lig", "Rec"]))
    A.obsm["spatial"] = np.asarray(f.coords, float)
    ct.tl.spatial_communication(A, database_name="pc",
                                df_ligrec=pd.DataFrame([["Lig", "Rec", "PW"]]),
                                dis_thr=100.0, pathway_sum=True)
    ct.tl.cluster_communication(A, database_name="pc", pathway_name="PW",
                                clustering="celltype", n_permutations=200,
                                random_seed=seed)
    r = A.uns["commot_cluster-celltype-pc-PW"]
    return (r["communication_matrix"].reindex(index=TYPES, columns=TYPES).values,
            r["communication_pvalue"].reindex(index=TYPES, columns=TYPES).values)


def main():
    rng = np.random.default_rng(SEED)
    xy = rng.uniform(0, SIDE, size=(N, 2))
    code = np.repeat(np.arange(K), N // K)
    code = code[rng.permutation(N)]
    rows = []

    # ---- C1: cell-type sensitivity -------------------------------------
    lig = np.where(code == 0, rng.random(N) < 0.30, rng.random(N) < 0.02).astype(float)
    rec = np.where(code == 1, rng.random(N) < 0.40, rng.random(N) < 0.02).astype(float)
    f = Fake(xy, code, lig, rec)
    res = _score(f, np.random.default_rng(SEED + 1))
    res["COMMOT"] = _commot(f, 7)
    for m, (S, P) in res.items():
        rows.append(dict(control="C1 cell-type sensitivity", method=m,
                         truth="A->B", p_AB=float(P[0, 1]),
                         p_CD=float(P[2, 3]), p_BA=float(P[1, 0]),
                         n_sig=int(np.nansum(P < 0.05)),
                         detected=bool(P[0, 1] < 0.05)))

    # ---- C2: spatial sensitivity ---------------------------------------
    # identical labels and identical MARGINAL ligand / receptor rates in both
    # conditions; only the ligand+ positions differ.
    lig = (rng.random(N) < 0.10).astype(float)
    tree = cKDTree(xy[lig > 0])
    near = np.array(tree.query_ball_point(xy, 25.0, workers=1, return_length=True)) > 0
    rec = np.zeros(N)
    cand = np.flatnonzero(near)
    take = rng.choice(cand, size=min(int(0.20 * N), cand.size), replace=False)
    rec[take] = 1.0                       # receptor+ ONLY next to ligand+ cells
    f_real = Fake(xy, code, lig, rec)
    xy_shift = xy.copy()
    m = lig > 0
    xy_shift[m] = D.torus_shift(np.random.default_rng(SEED + 3), xy[m],
                                xy.min(0), xy.max(0))
    f_null = Fake(xy_shift, code, lig, rec)
    for label, ff, sd in (("real (ligand+ adjacent to receptor+)", f_real, 11),
                          ("N3 torus-shifted ligand+", f_null, 11)):
        res = _score(ff, np.random.default_rng(SEED + 2))
        res["COMMOT"] = _commot(ff, sd)
        for meth, (S, P) in res.items():
            rows.append(dict(control="C2 spatial sensitivity", method=meth,
                             truth=label,
                             score_AB=float(S[0, 1]) if np.isfinite(S[0, 1]) else np.nan,
                             p_AB=float(P[0, 1]),
                             n_sig=int(np.nansum(P < 0.05)),
                             detected=bool(P[0, 1] < 0.05)))

    # ---- C3: WITHIN-cell-type spatial coupling -------------------------
    # A expresses the ligand, B the receptor, at IDENTICAL group-level rates in
    # both conditions.  The only difference is whether the receptor+ B cells sit
    # next to the ligand+ A cells or are scattered at random among the B cells.
    # This is precisely the extra thing a SPATIAL method claims to see.
    r3 = np.random.default_rng(SEED + 21)
    ligA = np.zeros(N)
    iA = np.flatnonzero(code == 0)
    ligA[r3.choice(iA, size=int(0.30 * iA.size), replace=False)] = 1.0
    iB = np.flatnonzero(code == 1)
    nrec = int(0.25 * iB.size)
    tree = cKDTree(xy[ligA > 0])
    dB = tree.query(xy[iB], k=1, workers=1)[0]
    coupled = np.zeros(N); coupled[iB[np.argsort(dB)[:nrec]]] = 1.0   # nearest to a ligand+ A
    scattered = np.zeros(N); scattered[r3.choice(iB, size=nrec, replace=False)] = 1.0
    for label, recv in (("coupled (receptor+ B next to ligand+ A)", coupled),
                        ("scattered (same rates, random B cells)", scattered)):
        ff = Fake(xy, code, ligA, recv)
        res = _score(ff, np.random.default_rng(SEED + 22))
        res["COMMOT"] = _commot(ff, 23)
        for meth, (S, P) in res.items():
            rows.append(dict(control="C3 within-type spatial coupling",
                             method=meth, truth=label,
                             score_AB=float(S[0, 1]) if np.isfinite(S[0, 1]) else np.nan,
                             p_AB=float(P[0, 1]),
                             n_sig=int(np.nansum(P < 0.05)),
                             detected=bool(P[0, 1] < 0.05)))

    # ---- C4: BETWEEN-cell-type geometry --------------------------------
    # Same expression as C1.  In one condition A and B are interleaved, in the
    # other they are segregated to opposite ends of the tile, far beyond the
    # 100 um interaction range.  This IS expressible at the cell-type level and
    # is the sanity check that the methods see geometry when it is theirs to see.
    r4 = np.random.default_rng(SEED + 31)
    lig4 = np.where(code == 0, r4.random(N) < 0.30, r4.random(N) < 0.02).astype(float)
    rec4 = np.where(code == 1, r4.random(N) < 0.40, r4.random(N) < 0.02).astype(float)
    xy_seg = xy.copy()
    xy_seg[code == 0, 0] = r4.uniform(0, 0.3 * SIDE, size=(code == 0).sum())
    xy_seg[code == 1, 0] = r4.uniform(0.7 * SIDE, SIDE, size=(code == 1).sum())
    for label, coords in (("interleaved A and B", xy), ("segregated A and B", xy_seg)):
        ff = Fake(coords, code, lig4, rec4)
        res = _score(ff, np.random.default_rng(SEED + 32))
        res["COMMOT"] = _commot(ff, 33)
        for meth, (S, P) in res.items():
            rows.append(dict(control="C4 between-type geometry", method=meth,
                             truth=label,
                             score_AB=float(S[0, 1]) if np.isfinite(S[0, 1]) else np.nan,
                             p_AB=float(P[0, 1]),
                             n_sig=int(np.nansum(P < 0.05)),
                             detected=bool(P[0, 1] < 0.05)))

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
