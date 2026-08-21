"""Section 8 Test 4: sender clustering, measured as the Ripley-K ratio at 50 um
against the N1 (within-cell-type label permutation) null.  Phase 1's synthetic
clustering axis kappa = 0 -> 4.5 corresponded to ratios 1.1-3.0 -> 10.9-12.6."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P, run_phase3_nulls as RN
from scipy.spatial import cKDTree
rows = []
for s in P.ALL_SECTIONS:
    sec = P.Sec(s); co = sec.coords.astype(float)
    elig = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    rng = np.random.default_rng(P.MASTER_SEED)
    for call in ("tierA_p95", "cdkn1a_pos", "senepy_p95"):
        snd = sec.sender_mask(call)
        if snd.sum() < 100: continue
        obs = np.asarray(cKDTree(co[snd]).query_ball_point(
            co[snd], r=50, return_length=True, workers=-1)).mean() - 1
        nulls = []
        for _ in range(10):
            m = RN.permute_within_type(rng, snd, sec.celltype, elig)
            nulls.append(np.asarray(cKDTree(co[m]).query_ball_point(
                co[m], r=50, return_length=True, workers=-1)).mean() - 1)
        rows.append(dict(section=s, band=sec.meta["band"], call=call,
                         obs=obs, null=float(np.mean(nulls)),
                         ripley_ratio=obs / float(np.mean(nulls))))
        print(rows[-1], flush=True)
d = pd.DataFrame(rows); d.to_csv(f"{P.RESULTS}/ripley.csv", index=False)
print(d.groupby("call").ripley_ratio.describe()[["mean","min","max"]].round(3))
