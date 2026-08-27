#!/usr/bin/env python3
"""Phase 10 — H1 caller agreement BEFORE conditioning (global top-5 %), for Figure 6a.

Phase 9 produced only the depth- and type-MATCHED agreement on H1
(`results/phase9_h1/caller_agreement_depth_and_type_matched.csv`), which is the
pre-registered primary rule (PREREG §3.7 b).  Figure 6a asks for both arms "before and
after conditioning on cell type and depth decile", and the mouse arm's *before* rows exist
(`results/phase3/caller_pairwise_agreement_11sections.csv`, `threshold == "global_top5"`).
This file produces the H1 counterpart with the identical arithmetic:

    flag   = score > the 95th percentile of that score over all finite-scored cells
    jaccard = |A and B| / |A or B|
    chance  = p_A * p_B * n / |A or B|
    ratio   = jaccard / chance

verified to reproduce the mouse file's own rows to 3 dp from its stored counts.

Usage: python3 code/h1_caller_unconditioned.py
Writes results/phase10_h1/caller_pairwise_agreement_global_h1.csv
"""
import os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H
from caller_disagree_all import SCORES

RES10 = "/workspace/results/phase10_h1"
PROC = H.PROC + "/"
EXCL = set(H.EXCLUDE_TYPES)


def load(section):
    sen = pd.read_csv(PROC + "senders_h1_%s.csv" % section).set_index("cell_id")
    ds = PROC + "deepscence_h1_%s.csv" % section
    if os.path.exists(ds):
        sen = sen.join(pd.read_csv(ds).set_index("cell_id"))
    return sen[~sen.cell_type.isin(EXCL)].copy()


if __name__ == "__main__":
    rows = []
    for s in H.ALL_SECTIONS:
        df = load(s)
        if "deepscence_score" in df:
            df["abs_deepscence_score"] = df.deepscence_score.abs()
        scores = [c for c in list(SCORES) + ["abs_deepscence_score"]
                  if c in df and df[c].notna().any()]
        for a, b in itertools.combinations(scores, 2):
            va, vb = df[a].to_numpy(float), df[b].to_numpy(float)
            ok = np.isfinite(va) & np.isfinite(vb)
            if ok.sum() < 200:
                continue
            x, y = va[ok], vb[ok]
            fa = x > np.nanpercentile(x, 95)
            fb = y > np.nanpercentile(y, 95)
            uni = int((fa | fb).sum()); inter = int((fa & fb).sum())
            if uni == 0:
                continue
            jac = inter / uni
            ch = fa.mean() * fb.mean() * ok.sum() / uni
            rows.append(dict(section=s, threshold="global_top5", A=a, B=b,
                             n=int(ok.sum()), n_A=int(fa.sum()), n_B=int(fb.sum()),
                             n_both=inter, jaccard=round(float(jac), 5),
                             chance_jaccard=round(float(ch), 5),
                             ratio=round(float(jac / ch), 3) if ch else np.nan))
    d = pd.DataFrame(rows)
    os.makedirs(RES10, exist_ok=True)
    d.to_csv(RES10 + "/caller_pairwise_agreement_global_h1.csv", index=False)
    print(d.groupby(["A", "B"]).ratio.describe()[["count", "min", "50%", "max"]].to_string())
    print("->", RES10 + "/caller_pairwise_agreement_global_h1.csv")
