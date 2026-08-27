#!/usr/bin/env python3
"""Phase 9 test A5 — matched-decoy contrast, |SMD| <= 0.1 after matching.  GO/NO-GO.

Master Plan §8 Test 5.  The matcher is `phase3_core.match_decoys_section` verbatim -- greedy
1-1 nearest-neighbour propensity matching without replacement, within section and within cell
type, caliper 0.25 SD -- run on the H1 cache through `h1_sec`.  Matching set is the published
N2 set built by `phase3_core.build_blocks`: log density at 50 um, log transcript counts, the
arm's anatomical covariate (here the A6 red/white-pulp axis) and the 20-NN cell-type
composition vector.

PASS condition: max |SMD| after matching <= 0.1, per (section x sender call).
M1 reference: max |SMD| 0.0916 -> 0.0352, gate passes in 100 % of matches
(PREREG_PHASE8.md §3.8).

Usage: python3 code/h1_a5_matchbalance.py [SPLN07 ...]
Writes results/phase9_h1/a5_match_balance.csv and a5_smd_by_covariate.csv
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_sec                      # noqa: F401
import sasp_phase3 as P
import phase3_core as C
import h1_common as H

CALLS = ("tierA_p90", "tierA_p95", "tierA_p99", "cdkn1a_pos", "senepy_p95")
secs = sys.argv[1:] or list(H.ALL_SECTIONS)
rows, cov = [], []
for i, s in enumerate(secs):
    sec = P.Sec(s)
    blocks = C.build_blocks(sec, np.zeros(sec.n, bool))
    for j, call in enumerate(CALLS):
        snd = sec.sender_mask(call)
        if snd.sum() < 100:
            print("  skip", s, call, snd.sum()); continue
        seed = P.MASTER_SEED + 1000 * i + j
        m = C.match_decoys_section(sec, snd, blocks["Zmatch"], seed)
        rows.append(dict(section=s, call=call, seed=seed,
                         n_senders=int(snd.sum()), prevalence=float(snd.mean()),
                         n_matched=int(len(m["sender_matched"])),
                         match_rate=round(m["match_rate"], 6),
                         max_smd_before=round(m["max_smd_before"], 4),
                         max_smd_after=round(m["max_smd_after"], 4),
                         passes_A5=bool(m["max_smd_after"] <= 0.10)))
        for k, nm in enumerate(blocks["zmatch_cols"]):
            cov.append(dict(section=s, call=call, covariate=nm,
                            smd_before=round(float(m["smd_before"][k]), 4),
                            smd_after=round(float(m["smd_after"][k]), 4)))
        print(rows[-1], flush=True)
d = pd.DataFrame(rows); d.to_csv(H.RESULTS + "/a5_match_balance.csv", index=False)
pd.DataFrame(cov).to_csv(H.RESULTS + "/a5_smd_by_covariate.csv", index=False)
pd.set_option("display.width", 220)
print("\n=== A5 ===")
print(d.to_string(index=False))
print("\npass rate: %d/%d (%.1f%%)   max |SMD| after matching over all rows: %.4f"
      % (d.passes_A5.sum(), len(d), 100 * d.passes_A5.mean(), d.max_smd_after.max()))
print("M1 reference: 0.0916 -> 0.0352, 100 % of matches pass")
