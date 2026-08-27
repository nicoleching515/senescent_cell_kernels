#!/usr/bin/env python3
"""Phase 9 test A4 — sender clustering, Ripley-K ratio at 50 um against the N1 null.

`code/_ripley.py` verbatim, with the H1 cache substituted via `h1_sec` and the caller list
extended to the same three families A3 evaluates.  M1 reference: 1.11 (tierA_p95) /
1.26 (cdkn1a_pos) / 1.56 (senepy_p95) -- see results/phase3/ripley.csv.

Usage: python3 code/h1_a4_ripley.py [SPLN07 ...]
Writes results/phase9_h1/a4_ripley.csv
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_sec                      # noqa: F401  side effect: repoint the cache
import sasp_phase3 as P, run_phase3_nulls as RN
import h1_common as H
from scipy.spatial import cKDTree

secs = sys.argv[1:] or list(H.ALL_SECTIONS)
rows = []
for s in secs:
    sec = P.Sec(s); co = sec.coords.astype(float)
    elig = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    rng = np.random.default_rng(P.MASTER_SEED)
    for call in ("tierA_p90", "tierA_p95", "tierA_p99", "cdkn1a_pos",
                 "senepy_p90", "senepy_p95", "senepy_p99"):
        try:
            snd = sec.sender_mask(call)
        except Exception as e:
            print("  skip", s, call, e); continue
        if snd.sum() < 100:
            print("  skip", s, call, "only %d senders" % snd.sum()); continue
        obs = np.asarray(cKDTree(co[snd]).query_ball_point(
            co[snd], r=50, return_length=True, workers=-1)).mean() - 1
        nulls = []
        for _ in range(10):
            m = RN.permute_within_type(rng, snd, sec.celltype, elig)
            nulls.append(np.asarray(cKDTree(co[m]).query_ball_point(
                co[m], r=50, return_length=True, workers=-1)).mean() - 1)
        rows.append(dict(section=s, call=call, n_senders=int(snd.sum()),
                         prevalence=float(snd.mean()), obs=obs,
                         null=float(np.mean(nulls)), null_sd=float(np.std(nulls)),
                         ripley_ratio=obs / float(np.mean(nulls))))
        print(rows[-1], flush=True)
d = pd.DataFrame(rows)
d.to_csv(H.RESULTS + "/a4_ripley.csv", index=False)
print(d.groupby("call").ripley_ratio.describe()[["mean", "min", "max"]].round(3))
print("M1 reference: tierA_p95 1.11 / cdkn1a_pos 1.26 / senepy_p95 1.56")
print("wrote", H.RESULTS + "/a4_ripley.csv")
