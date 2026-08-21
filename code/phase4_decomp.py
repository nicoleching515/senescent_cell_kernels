"""How much of an edge-averaged LR score is spatial at all?

Any statistic of the form "average over A-B neighbour edges of f(sender)*g(receiver)"
has expectation E[f|A]*E[g|B] under a random rewiring of the graph.  The observed
value differs from that only by the spatial covariance of f and g across edges.
This measures the size of that difference -- i.e. the fraction of the score that
carries any geometry -- for the SpaTalk*-style statistic, analytically rather
than by permutation.
"""
import numpy as np, pandas as pd
import phase4_data as D
import phase4_methods as M
from phase4_tiles import all_tiles


def main():
    rows = []
    for t in all_tiles():
        src, dst = M.spatalk_edges(t.coords)
        key = t.code[src] * t.K + t.code[dst]
        cnt = np.bincount(key, minlength=t.K ** 2).astype(float)
        for p in D.LR_PAIRS:
            lig, rec = t.g(p["ligand"]), t.receptor(p["receptors"])
            obs, _ = M._spatalk_stat(lig, rec, src, dst, key, cnt, t.K), None
            # composition-only prediction: cell-type means, no graph at all
            L = np.array([lig[t.code == k].mean() for k in range(t.K)])
            R = np.array([rec[t.code == k].mean() for k in range(t.K)])
            gpred = np.sqrt(L)[:, None] * np.sqrt(R)[None, :]
            pred = gpred / (gpred + 1.0)
            m = (cnt.reshape(t.K, t.K) > 0) & (obs > 0)
            if m.sum() < 4:
                continue
            rows.append(dict(
                tile=t.name, section=t.section, pair=p["pair"],
                spearman_obs_vs_composition=float(
                    pd.Series(obs[m]).corr(pd.Series(pred[m]), method="spearman")),
                median_abs_rel_gap=float(np.median(np.abs(obs[m] - pred[m]) /
                                                   np.maximum(obs[m], 1e-12))),
                n_blocks=int(m.sum())))
    df = pd.DataFrame(rows)
    df.to_csv("/workspace/results/phase4/spatial_content_of_edge_scores.csv",
              index=False)
    print(df.groupby("pair")[["spearman_obs_vs_composition",
                              "median_abs_rel_gap"]].median().round(4).to_string())
    print("\noverall:", df.spearman_obs_vs_composition.median().round(4),
          df.median_abs_rel_gap.median().round(4))


if __name__ == "__main__":
    main()
