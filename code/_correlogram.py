"""Where do the real sections sit on the Phase 1 regime map?

Phase 1's x-axis is the baseline autocorrelation length ell relative to
lambda_true.  Here we measure ell directly on the real sections: the spatial
correlogram of each Tier B module after removing the receiver-cell-type mean,
fitted with an exponential, on a 20k-cell subsample (cKDTree, never an (n,n)
matrix).  The same is done for the zonation score and the Tier A sender score.
"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
from scipy.spatial import cKDTree

BINS = np.arange(0, 200.1, 10.0)
rows = []
for s in P.ALL_SECTIONS:
    sec = P.Sec(s)
    co = sec.coords.astype(float)
    rng = np.random.default_rng(P.MASTER_SEED)
    sub = rng.choice(co.shape[0], size=min(20000, co.shape[0]), replace=False)
    tree = cKDTree(co)
    pairs = tree.query_ball_point(co[sub], r=BINS[-1], workers=-1)
    ii = np.repeat(sub, [len(p) for p in pairs])
    jj = np.concatenate([np.asarray(p, int) for p in pairs])
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    dd = np.linalg.norm(co[ii] - co[jj], axis=1)
    bidx = np.clip(np.digitize(dd, BINS) - 1, 0, BINS.size - 2)
    fields = {m: sec.module(m) for m in P.MODULES}
    fields["zonation_score"] = np.nan_to_num(sec.zonation_score.astype(float))
    fields["tierA_score"] = np.nan_to_num(sec.tierA_score.astype(float))
    ct = sec.celltype
    for name, v in fields.items():
        v = v.astype(float).copy()
        for c in np.unique(ct):                    # remove cell-type means
            m = ct == c
            v[m] -= v[m].mean()
        v /= (v.std() + 1e-12)
        num = np.bincount(bidx, weights=v[ii] * v[jj], minlength=BINS.size - 1)
        cnt = np.bincount(bidx, minlength=BINS.size - 1).astype(float)
        rho = num / np.maximum(cnt, 1)
        x = 0.5 * (BINS[:-1] + BINS[1:])
        ok = (cnt > 100) & (rho > 0)
        ell = np.nan
        if ok.sum() >= 3:
            A = np.column_stack([np.ones(ok.sum()), -x[ok]])
            b = np.linalg.lstsq(A, np.log(rho[ok]), rcond=None)[0]
            ell = 1.0 / b[1] if b[1] > 0 else np.inf
        rows.append(dict(section=s, arm=sec.meta["condition"], field=name,
                         rho_10um=rho[0], rho_30um=rho[2], rho_50um=rho[4],
                         rho_100um=rho[9], ell_um=ell))
        print(f"{s[:4]} {name:22s} rho10={rho[0]:+.3f} rho50={rho[4]:+.3f} "
              f"rho100={rho[9]:+.3f} ell={ell:7.1f} um", flush=True)
pd.DataFrame(rows).to_csv(f"{P.RESULTS}/correlogram.csv", index=False)
print("wrote correlogram.csv")
