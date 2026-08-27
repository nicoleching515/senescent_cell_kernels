#!/usr/bin/env python3
"""Phase 9 Job B step 1 (validation) — cross-check our marker-based cell types against the
DEPOSITORS' four-level annotation.

§14 assumes the bio collaborator "eyeballs fifty called cells against the image".  The images
were deliberately not downloaded (§12.3) and this deposit ships something stronger: four
nested annotation levels over 89-99 % of the cells (PREREG_PHASE8.md P21).  This script
characterises that filter (done in `h1_a1_geometry.py`) and then measures agreement.

Reported:
  * confusion of our fine and merged labels against Level_1..Level_4
  * per-label purity: for each of our labels, the modal depositor label and its share
  * per-depositor-label recall: for each depositor Level_3 label, the modal label we give it
  * adjusted Rand index and normalised mutual information, our labels vs each depositor level
  * which of our 23 labels are never realised, and which depositor labels we never recover
  * agreement restricted to cells BOTH label sets call non-low-quality

Usage: python3 code/h1_jobB_crosscheck.py [SPLN07 ...]
Writes results/phase9_h1/jobB_crosscheck_*.csv
"""
import sys, os
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
sys.path.insert(0, "/workspace/code")
import h1_common as H

PROC = H.PROC + "/"
LEVELS = ["Level_1_Annotations", "Level_2_Annotations",
          "Level_3_Annotations", "Level_4_Annotations"]
OURS = ["cell_type", "cell_type_merged"]


def run(section):
    ct = pd.read_csv(PROC + "celltypes_h1_%s.csv" % section).set_index("cell_id")
    ann = H.annotations(section).set_index("cell_id")
    j = ct.join(ann, how="inner")
    n_ours, n_dep, n_both = len(ct), len(ann), len(j)
    rows, purity, recall, conf = [], [], [], []
    for ocol in OURS:
        o = j[ocol].astype(str)
        for lcol in LEVELS:
            d = j[lcol].astype(str)
            # restrict to cells neither side calls low quality
            ok = (~o.isin(["Low_quality", "Unknown"])) & (d != "Low quality")
            rows.append(dict(section=section, our_label=ocol, depositor_level=lcol,
                             n_compared=int(len(j)), n_clean=int(ok.sum()),
                             ari_all=round(adjusted_rand_score(d, o), 4),
                             nmi_all=round(normalized_mutual_info_score(d, o), 4),
                             ari_clean=round(adjusted_rand_score(d[ok], o[ok]), 4),
                             nmi_clean=round(normalized_mutual_info_score(d[ok], o[ok]), 4),
                             n_our_labels=int(o.nunique()), n_dep_labels=int(d.nunique())))
        d3 = j["Level_3_Annotations"].astype(str)
        for lab, g in j.groupby(ocol):
            v = g["Level_3_Annotations"].astype(str).value_counts()
            purity.append(dict(section=section, our_label_set=ocol, our_label=lab,
                               n=int(len(g)), modal_depositor_L3=v.index[0],
                               modal_share=round(float(v.iloc[0] / len(g)), 4),
                               second=v.index[1] if len(v) > 1 else "",
                               second_share=round(float(v.iloc[1] / len(g)), 4)
                               if len(v) > 1 else 0.0))
        for lab, g in j.groupby("Level_3_Annotations"):
            v = g[ocol].astype(str).value_counts()
            recall.append(dict(section=section, our_label_set=ocol, depositor_L3=lab,
                               n=int(len(g)), modal_our_label=v.index[0],
                               modal_share=round(float(v.iloc[0] / len(g)), 4),
                               frac_unknown=round(float((g[ocol] == "Unknown").mean()), 4),
                               frac_lowquality=round(float((g[ocol] == "Low_quality").mean()), 4)))
        if ocol == "cell_type":
            c = pd.crosstab(d3, o)
            c.insert(0, "section", section)
            conf.append(c.reset_index())
    return (pd.DataFrame(rows), pd.DataFrame(purity), pd.DataFrame(recall),
            pd.concat(conf) if conf else pd.DataFrame(),
            dict(section=section, n_our=n_ours, n_depositor=n_dep, n_joined=n_both))


if __name__ == "__main__":
    secs = sys.argv[1:] or [s for s in H.ALL_SECTIONS
                            if os.path.exists(PROC + "celltypes_h1_%s.csv" % s)]
    R, P_, C, X, M = [], [], [], [], []
    for s in secs:
        print("...", s, flush=True)
        r, p, c, x, m = run(s)
        R.append(r); P_.append(p); C.append(c); X.append(x); M.append(m)
    pd.concat(R).to_csv(H.RESULTS + "/jobB_crosscheck_scores.csv", index=False)
    pd.concat(P_).to_csv(H.RESULTS + "/jobB_crosscheck_purity.csv", index=False)
    pd.concat(C).to_csv(H.RESULTS + "/jobB_crosscheck_recall.csv", index=False)
    pd.concat(X).to_csv(H.RESULTS + "/jobB_crosscheck_confusion_L3.csv", index=False)
    pd.DataFrame(M).to_csv(H.RESULTS + "/jobB_crosscheck_coverage.csv", index=False)
    d = pd.concat(R)
    pd.set_option("display.width", 250)
    print(d[d.our_label == "cell_type"].to_string(index=False))
    print("wrote", H.RESULTS + "/jobB_crosscheck_*.csv")
