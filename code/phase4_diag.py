"""Phase 4 diagnostics: how destructive is each coordinate null, physically?

A torus shift moves ligand+ cells into space the tissue does not occupy.  Cells
that land in a void send nothing, which makes the SHIFTED condition weaker than
it would be under an ideal within-tissue shift -- so every surviving-fraction we
report is CONSERVATIVE.  This quantifies that.
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import phase4_data as D
from phase4_tiles import all_tiles

def main(n_rep=20, seed=20260820):
    rows = []
    for tile in all_tiles():
        T = cKDTree(tile.coords)
        rng = np.random.default_rng(seed + hash(tile.name) % 10000)
        for p in D.LR_PAIRS:
            m = tile.gc(p["ligand"]) > 0
            base = np.array(T.query_ball_point(tile.coords[m], D.WINDOW_UM,
                                               workers=1, return_length=True))
            for null in ("N3_lig", "N4_lig"):
                keep, dcen = [], []
                for _ in range(n_rep):
                    xy = D.null_coords(tile, null, rng, lig_mask=m)
                    if xy is None:
                        break
                    cnt = np.array(T.query_ball_point(xy[m], D.WINDOW_UM,
                                                      workers=1, return_length=True))
                    keep.append(float((cnt > 0).mean()))
                    dcen.append(float(np.median(cnt)))
                if not keep:
                    continue
                rows.append(dict(tile=tile.name, section=tile.section,
                                 pair=p["pair"], null=null,
                                 n_ligand_pos=int(m.sum()),
                                 pct_ligand_pos=round(100 * float(m.mean()), 3),
                                 real_median_nbrs=float(np.median(base)),
                                 null_median_nbrs=float(np.mean(dcen)),
                                 frac_retaining_a_neighbour=float(np.mean(keep))))
    df = pd.DataFrame(rows)
    df.to_csv("/workspace/results/phase4/null_destructiveness.csv", index=False)
    print(df.groupby(["null"])[["frac_retaining_a_neighbour", "real_median_nbrs",
                                "null_median_nbrs"]].median().round(3).to_string())
    print(df.groupby("pair").n_ligand_pos.describe()[["min", "50%", "max"]].to_string())
    return df

if __name__ == "__main__":
    main()
