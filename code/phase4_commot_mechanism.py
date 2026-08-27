"""Why COMMOT's answer does not move when the coordinates do.

COMMOT solves a collective optimal-transport problem: the total ligand mass
that gets transported is fixed by expression, and the geometry decides only
WHICH cell pairs it flows between.  Averaging that flow over a sender x receiver
cell-type block therefore recovers something close to (total mass x composition),
which carries no geometry at all -- and the cluster-level test then permutes
labels with the transport plan HELD FIXED, so it never tests the geometry
either.

This script measures that directly: real vs fully coordinate-permuted, at the
cell-pair level and at the cell-type level, on the first tile of every section.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import spearmanr
import _shims.np2_compat            # noqa: F401  (MUST precede commot: commot 0.0.3 has
                                    # `rho = np.Inf` as a default argument and raises at
                                    # import under the pinned numpy 2.4.6)
import commot as ct
import phase4_data as D
from phase4_tiles import all_tiles

OUT = "/workspace/results/phase4/commot_mechanism.csv"


def run(tile, xy, pdef, seed):
    A = ad.AnnData(np.asarray(tile.X, dtype=np.float64),
                   obs=pd.DataFrame({"celltype": tile.celltype},
                                    index=[str(i) for i in range(tile.n)]),
                   var=pd.DataFrame(index=tile.genes))
    A.obsm["spatial"] = np.asarray(xy, float)
    df = pd.DataFrame([[pdef["ligand"], r, pdef["pathway"]] for r in pdef["receptors"]])
    ct.tl.spatial_communication(A, database_name="s", df_ligrec=df, dis_thr=100.,
                                pathway_sum=True)
    ct.tl.cluster_communication(A, database_name="s", pathway_name=pdef["pathway"],
                                clustering="celltype", n_permutations=100,
                                random_seed=seed)
    S = A.obsp["commot-s-" + pdef["pathway"]].tocoo()
    r = A.uns["commot_cluster-celltype-s-" + pdef["pathway"]]
    return S, r["communication_matrix"].values, r["communication_pvalue"].values


def main():
    rows = []
    tiles = [t for t in all_tiles() if t.name.endswith("_t0")]
    for tile in tiles:
        rng = np.random.default_rng(20260820 + hash(tile.name) % 9973)
        xy1 = D.null_coords(tile, "N0_perm", rng)
        for pdef in D.LR_PAIRS:
            S0, C0, P0 = run(tile, tile.coords, pdef, 5)
            S1, C1, P1 = run(tile, xy1, pdef, 5)
            k0 = set(zip(S0.row.tolist(), S0.col.tolist()))
            k1 = set(zip(S1.row.tolist(), S1.col.tolist()))
            rows.append(dict(
                tile=tile.name, section=tile.section, pair=pdef["pair"],
                n_edges_real=len(k0), n_edges_perm=len(k1),
                edge_jaccard=round(len(k0 & k1) / max(1, len(k0 | k1)), 4),
                total_mass_real=float(S0.data.sum()),
                total_mass_perm=float(S1.data.sum()),
                cluster_spearman=round(float(spearmanr(C0.ravel(), C1.ravel()).statistic), 4),
                nsig_real=int((P0 < 0.05).sum()), nsig_perm=int((P1 < 0.05).sum()),
                nsig_shared=int(((P0 < 0.05) & (P1 < 0.05)).sum()),
                n_cells=tile.n, n_types=tile.K))
            print(rows[-1], flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(df[["edge_jaccard", "cluster_spearman"]].describe().round(4).to_string())
    print("mass ratio perm/real:",
          (df.total_mass_perm / df.total_mass_real).describe().round(6).to_string())


if __name__ == "__main__":
    main()
