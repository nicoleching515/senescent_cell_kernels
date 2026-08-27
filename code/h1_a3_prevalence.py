#!/usr/bin/env python3
"""Phase 9 test A3 — SENDER PREVALENCE.  The hard gate on the human arm.

Master Plan §8 Test 3 / Phase 7 §13 A3: a (section x cell type x caller x threshold) stratum
passes iff sender prevalence is in **1-20 %** with **>= 200 senders** and **>= 5,000
non-senders**.  §14 step 3 requires it PER CELL TYPE, not pooled.  PREREG_PHASE8.md §3.11
fixes that on H1 A3 is a *reported* quantity per section and per caller, not an exclusion
rule -- all 7 sections are analysed either way -- and that for SenePy it is evaluated on the
subset of cells that receive a score at all, **with both denominators stated**.

Thresholds evaluated (PREREG §3.7):
  tierA_p90 / tierA_p95 (PRIMARY) / tierA_p99   strict > within-cell-type percentile
  senepy_p90 / senepy_p95 / senepy_p99          same rule on the SenePy surrogate score
  cdkn1a_pos                                    CDKN1A counts > 0 -- the caller's own cutoff
  deepscence_own                                DeepScence's own binary call, if available
  tierApm_p95__<module>                         the D1 per-module sensitivity

NOTHING IS TUNED.  No threshold outside the frozen list is searched, and the outcome is
reported whichever way it falls.

Usage: python3 code/h1_a3_prevalence.py [SPLN07 ...]      (default: all sections present)
Writes results/phase9_h1/a3_prevalence_by_type.csv, a3_summary_by_caller.csv,
       a3_pooled_by_section.csv
"""
import sys, os, glob
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H

PROC = H.PROC + "/"
MIN_SENDERS, MIN_NONSENDERS = 200, 5000
LO, HI = 1.0, 20.0
EXCL = set(H.EXCLUDE_TYPES) | set(H.EXCLUDE_FROM_SENDERS)


def rows_for(section):
    sen = pd.read_csv(PROC + "senders_h1_%s.csv" % section)
    ds = PROC + "deepscence_h1_%s.csv" % section
    if os.path.exists(ds):
        sen = sen.merge(pd.read_csv(ds), on="cell_id", how="left")
    calls = {}
    for q in (90, 95, 99):
        calls["tierA_p%d" % q] = sen["sender_flag_p%d" % q].to_numpy(bool)
        if "senepy_flag_p%d" % q in sen:
            calls["senepy_p%d" % q] = sen["senepy_flag_p%d" % q].to_numpy(bool)
    calls["cdkn1a_pos"] = (sen.cdkn1a_counts > 0).to_numpy()
    for c in sen.columns:
        if c.startswith("tierApm_flag_p95__"):
            calls["tierApm_p95__" + c.split("__", 1)[1]] = sen[c].to_numpy(bool)
    if "deepscence_call" in sen:
        calls["deepscence_own"] = sen.deepscence_call.fillna(0).to_numpy(bool)
    if "deepscence_score" in sen:
        s = sen.deepscence_score.to_numpy(float)
        for q in (90, 95, 99):
            f = np.zeros(len(sen), bool)
            for c in pd.unique(sen.cell_type):
                if c in EXCL:
                    continue
                m = (sen.cell_type == c).to_numpy() & np.isfinite(s)
                if m.sum() < 100:
                    continue
                f[m] = s[m] > np.nanpercentile(s[m], q)
            calls["deepscence_p%d" % q] = f

    scored = {}                    # per-caller denominator restriction
    for k in calls:
        if k.startswith("senepy"):
            scored[k] = np.isfinite(sen.senepy_score.to_numpy(float))
        elif k.startswith("deepscence") and "deepscence_score" in sen:
            scored[k] = np.isfinite(sen.deepscence_score.to_numpy(float))
        else:
            scored[k] = np.ones(len(sen), bool)

    out, pooled = [], []
    for labcol, labname in (("cell_type", "fine"), ("cell_type_merged", "merged")):
        lab = sen[labcol].astype(str).to_numpy()
        for call, mask in calls.items():
            sc = scored[call]
            elig_all = (~np.isin(lab, list(EXCL)))
            for t in sorted(set(lab[elig_all])):
                m = (lab == t) & elig_all
                # denominator A: all cells of this type
                nA = int(m.sum()); sA = int((m & mask).sum())
                # denominator B: only cells this caller can score
                mB = m & sc
                nB = int(mB.sum()); sB = int((mB & mask).sum())
                pA = 100 * sA / nA if nA else np.nan
                pB = 100 * sB / nB if nB else np.nan
                out.append(dict(
                    section=section, label_set=labname, cell_type=t, call=call,
                    n_all=nA, n_senders_all=sA, n_nonsenders_all=nA - sA,
                    prevalence_pct_all=round(pA, 3) if nA else None,
                    n_scored=nB, n_senders_scored=sB, n_nonsenders_scored=nB - sB,
                    prevalence_pct_scored=round(pB, 3) if nB else None,
                    pass_band=bool(nB and LO <= pB <= HI),
                    pass_senders=bool(sB >= MIN_SENDERS),
                    pass_nonsenders=bool(nB - sB >= MIN_NONSENDERS),
                    passes_A3=bool(nB and LO <= pB <= HI and sB >= MIN_SENDERS
                                   and nB - sB >= MIN_NONSENDERS)))
            if labname == "fine":
                m = elig_all & sc
                s_ = int((m & mask).sum()); n_ = int(m.sum())
                pooled.append(dict(section=section, call=call, n_scored=n_, n_senders=s_,
                                   prevalence_pct=round(100 * s_ / n_, 3) if n_ else None,
                                   n_all=int(elig_all.sum()),
                                   frac_cells_scored=round(n_ / max(int(elig_all.sum()), 1), 4)))
    return out, pooled


if __name__ == "__main__":
    secs = sys.argv[1:] or [s for s in H.ALL_SECTIONS
                            if os.path.exists(PROC + "senders_h1_%s.csv" % s)]
    allr, allp = [], []
    for s in secs:
        r, p = rows_for(s); allr += r; allp += p
        print("A3 rows for", s, len(r), flush=True)
    d = pd.DataFrame(allr); d.to_csv(H.RESULTS + "/a3_prevalence_by_type.csv", index=False)
    pd.DataFrame(allp).to_csv(H.RESULTS + "/a3_pooled_by_section.csv", index=False)

    fine = d[d.label_set == "fine"]
    g = (fine.groupby("call")
              .agg(n_strata=("passes_A3", "size"),
                   n_pass=("passes_A3", "sum"),
                   n_pass_band=("pass_band", "sum"),
                   n_pass_senders=("pass_senders", "sum"),
                   n_pass_nonsenders=("pass_nonsenders", "sum"),
                   median_prevalence=("prevalence_pct_scored", "median"),
                   min_prevalence=("prevalence_pct_scored", "min"),
                   max_prevalence=("prevalence_pct_scored", "max"))
              .reset_index())
    g["pct_pass"] = (100 * g.n_pass / g.n_strata).round(1)
    # a cell type "qualifies" only if it passes in every one of the 7 sections
    qual = (fine.groupby(["call", "cell_type"]).passes_A3
                .agg(["sum", "size"]).reset_index()
                .rename(columns={"sum": "n_sections_pass", "size": "n_sections"}))
    qual["all_sections"] = qual.n_sections_pass == qual.n_sections
    g = g.merge(qual.groupby("call").all_sections.sum().rename("n_types_pass_all_sections"),
                on="call", how="left")
    g.to_csv(H.RESULTS + "/a3_summary_by_caller.csv", index=False)
    qual.to_csv(H.RESULTS + "/a3_types_by_caller.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("\n=== A3 SUMMARY (fine labels, scored denominator) ===")
    print(g.to_string(index=False))
    print("\nwritten to", H.RESULTS)
