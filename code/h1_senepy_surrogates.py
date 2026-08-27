#!/usr/bin/env python3
"""Phase 9 — re-derive the SenePy cross-tissue surrogate map for the H1 spleen label set,
on the hub set the MOUSE ARM ACTUALLY SCORED WITH.

WHY THIS EXISTS (a discrepancy found in Phase 9, reported, not patched over):

  `code/senepy_coverage_human.py` — the producer of the frozen P4 numbers in
  results/phase7_jobA/senepy_spleen_coverage.{csv,json} — reads the hub pickle
  `6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle` directly.  That is senepy's **v1** human hub
  set (65 hubs).

  `code/phase2_downstream.py` — the mouse arm's actual SenePy scorer — calls
  `senepy.load_hubs(species='Mouse')`, whose default is `sig_version='v2'`
  (`R1_TMS_HUBS_DICTIONARY_FILTERED.pickle`).

  So the frozen human coverage assessment and the mouse arm's estimator are not the same
  hub release.  The human v2 set has 64 hubs, 13 of v1's keys are absent from it and 12 keys
  are new, and hub sizes differ by up to 4x (intestine/macrophage/0: 3,147 genes in v1,
  713 in v2).  Two of the frozen surrogate assignments — 'Fibroblastic reticular cells' and
  'Fibroblasts', both mapped to ('skin','fibroblast',1) — DO NOT EXIST in v2.

  H1 is scored on **v2**, because matching the mouse arm's estimator is what makes the two
  arms comparable, and P4's own point is that SenePy is not the same estimator across arms.
  Both tables are written and both are reported.

The surrogate selection rule is `senepy_coverage_human.SURROGATE` verbatim (word-boundary
match on the SenePy `cell` string, best hub by on-panel gene count) with the same
MIN_ON_PANEL = 10.  No new pattern is introduced.

Writes results/phase9_h1/senepy_surrogates_v1_v2.csv
"""
import os, pickle, re, sys
import numpy as np, pandas as pd, senepy
sys.path.insert(0, "/workspace/code")
import h1_common as H
from senepy_coverage_human import SURROGATE          # the frozen pattern table
from markers_human_spleen import MARKERS

MIN_ON_PANEL = 10
D = "/usr/local/lib/python3.11/dist-packages/senepy/data/"

_, panel_names, _, _ = H.load_matrix("SPLN21", "gene")
PANEL = set(panel_names)

v1 = pickle.load(open(D + "6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle", "rb"))
v2 = senepy.load_hubs(species="Human").hubs        # the default the mouse arm used
print("v1 hubs %d   v2 hubs %d   keys only in v1 %d   only in v2 %d"
      % (len(v1), len(v2), len(set(v1) - set(v2)), len(set(v2) - set(v1))))


def best(hubs, pats):
    b = None
    for (t, c, h), genes in hubs.items():
        if not any(re.search(r"\b" + re.escape(p) + r"\b", c) for p in pats):
            continue
        onp = [g for g, _ in genes if g in PANEL]
        if b is None or len(onp) > b[3]:
            b = (t, c, h, len(onp), len(genes))
    return b


labels = list(MARKERS) + ["Plasma cells"]           # + the P6 exception label
rows = []
for lab in labels:
    pats = SURROGATE.get(lab, [])
    r = dict(cell_type=lab)
    for tag, hubs in (("v1", v1), ("v2", v2)):
        b = best(hubs, pats) if pats else None
        if b is None:
            r.update({f"{tag}_tissue": "", f"{tag}_cell": "", f"{tag}_hub": "",
                      f"{tag}_size": 0, f"{tag}_on_panel": 0, f"{tag}_usable": "no_hub"})
        else:
            t, c, h, o, n = b
            r.update({f"{tag}_tissue": t, f"{tag}_cell": c, f"{tag}_hub": h,
                      f"{tag}_size": n, f"{tag}_on_panel": o,
                      f"{tag}_usable": "yes" if o >= MIN_ON_PANEL else "below_min_on_panel"})
    r["same_hub"] = (r["v1_tissue"], r["v1_cell"], r["v1_hub"]) == \
                    (r["v2_tissue"], r["v2_cell"], r["v2_hub"])
    rows.append(r)

d = pd.DataFrame(rows)
os.makedirs(H.RESULTS, exist_ok=True)
d.to_csv(H.RESULTS + "/senepy_surrogates_v1_v2.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(d[["cell_type", "v1_cell", "v1_on_panel", "v1_usable",
         "v2_tissue", "v2_cell", "v2_hub", "v2_on_panel", "v2_usable", "same_hub"]].to_string())
for tag in ("v1", "v2"):
    u = (d[f"{tag}_usable"] == "yes").sum()
    print("%s: %d/%d labels usable at >=%d on-panel genes, %d with no hub at all"
          % (tag, u, len(d), MIN_ON_PANEL, (d[f"{tag}_usable"] == "no_hub").sum()))
print("wrote", H.RESULTS + "/senepy_surrogates_v1_v2.csv")
