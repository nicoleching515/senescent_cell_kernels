#!/usr/bin/env python3
"""Phase 10 ISO -- the intersected-vs-full comparison tables (test A8 §8.3, both arms).

Two products, both built by REUSING the frozen arithmetic rather than restating it:

1. `--stage headlines` -> results/phase10_h1_iso/iso_vs_full_headlines.csv
   One row per (arm, call, panel in {full, intersected}).  Every number comes from
   `h1_headlines.headlines(call, r3=..., r5=..., sections=...)`, which is `m1_headlines.py`'s
   formulas line for line, pointed at a results directory and a section list:

     arm=h1 panel=full          results/phase10_h1        sections = h1_common.ALL_SECTIONS
     arm=h1 panel=intersected   results/phase10_h1_iso    same 7 sections
     arm=m1 panel=full          results/phase3            sections = sasp_phase3.IN_BAND
     arm=m1 panel=intersected   results/phase10_m1_iso    same 6 sections

   `headlines()` filters `stratum == "all"` itself, so the M1 zonation strata are excluded
   on both M1 rows, matching `m1_headlines.py`.  Calls: H1 {tierAmg_p95, tierA_p95},
   M1 {tierA_p95}.

2. `--stage shift --arm {h1,m1}` -> results/phase10_<arm>_iso/iso_sender_shift.csv
   Per section and call: the Spearman correlation between the FULL-panel and
   INTERSECTED-panel Tier A score, and the Jaccard of the top-5 % sender SETS.  This is the
   measurement `h1_a8_crossarm.score_restricted` makes for H1 in Phase 9, extended to M1 and
   to the merged-label call.  The sender sets are taken from `sasp_phase3.Sec.sender_mask`
   on each cache, i.e. AFTER the frozen Low_quality/Unknown/Proliferating exclusion, so they
   are the sets the fits actually used -- not a re-derived percentile.

Usage
  python3 code/iso_compare.py --stage shift --arm h1
  python3 code/iso_compare.py --stage shift --arm m1
  python3 code/iso_compare.py --stage headlines
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")

H1_ISO = "/workspace/results/phase10_h1_iso"
M1_ISO = "/workspace/results/phase10_m1_iso"
H1_FULL = "/workspace/results/phase10_h1"
M1_FULL = "/workspace/results/phase3"

COLS = ["arm", "call", "panel", "n_fits", "n_reportable", "frac_reportable",
        "naive_amp_med", "ctrl_amp_med", "power80_bound",
        "SF_N2", "SF_N5", "SF_N6", "SF_N5+N6", "SF_N2+N5+N6",
        "SF_N2+N5+N6_iqr_lo", "SF_N2+N5+N6_iqr_hi",
        "lam_railed_frac", "lam_naive_median",
        "sender_prevalence_min", "sender_prevalence_max"]


def _row(arm, call, panel, r3, sections):
    import h1_headlines as HL
    o = HL.headlines(call, r3=r3, r5=r3, sections=sections)
    lo, hi = o["SF_N2+N5+N6_iqr"]
    d = {k: o.get(k) for k in COLS if k in o}
    d.update(arm=arm, call=call, panel=panel,
             **{"SF_N2+N5+N6_iqr_lo": lo, "SF_N2+N5+N6_iqr_hi": hi})
    return d


def stage_headlines():
    import h1_common as H
    import sasp_phase3 as P
    rows = []
    for call in ("tierAmg_p95", "tierA_p95"):
        rows.append(_row("h1", call, "full", H1_FULL, list(H.ALL_SECTIONS)))
        rows.append(_row("h1", call, "intersected", H1_ISO, list(H.ALL_SECTIONS)))
    rows.append(_row("m1", "tierA_p95", "full", M1_FULL, list(P.IN_BAND)))
    rows.append(_row("m1", "tierA_p95", "intersected", M1_ISO, list(P.IN_BAND)))
    df = pd.DataFrame(rows)[COLS]
    os.makedirs(H1_ISO, exist_ok=True)
    df.to_csv(H1_ISO + "/iso_vs_full_headlines.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False))
    return df


def _sec_pair(name, full_cache, iso_cache):
    """(full Sec, iso Sec) for one section, by rebinding sasp_phase3.CACHE3."""
    import sasp_phase3 as P
    P.CACHE3 = full_cache
    a = P.Sec(name)
    P.CACHE3 = iso_cache
    b = P.Sec(name)
    return a, b


def stage_shift(arm):
    import sasp_phase3 as P
    if arm == "h1":
        import h1_phase10_iso as BIND           # binds ISO cache + the tierAmg shim
        import h1_common as H
        sections, calls = list(H.ALL_SECTIONS), ["tierAmg_p95", "tierA_p95"]
        full_cache = H.PROC + "/cache3_h1"
        iso_cache = BIND.CACHE_ISO
        out = H1_ISO + "/iso_sender_shift.csv"
    else:
        import m1_phase10_iso as BIND
        sections, calls = list(P.IN_BAND), ["tierA_p95"]
        full_cache = "/workspace/data/processed/cache3"
        iso_cache = BIND.CACHE_ISO
        out = M1_ISO + "/iso_sender_shift.csv"

    rows = []
    for s in sections:
        a, b = _sec_pair(s, full_cache, iso_cache)
        assert np.array_equal(a.cell_id, b.cell_id), s + ": cell order differs"
        sa = pd.Series(a.tierA_score.astype(float))
        sb = pd.Series(b.tierA_score.astype(float))
        rho = float(sa.corr(sb, method="spearman"))
        rho_p = float(sa.corr(sb, method="pearson"))
        for call in calls:
            ma = a.sender_mask(call)
            mb = b.sender_mask(call)
            inter = int((ma & mb).sum()); union = int((ma | mb).sum())
            rows.append(dict(arm=arm, section=s, call=call, n_cells=int(a.n),
                             n_senders_full=int(ma.sum()), n_senders_iso=int(mb.sum()),
                             prevalence_full=round(100 * float(ma.mean()), 4),
                             prevalence_iso=round(100 * float(mb.mean()), 4),
                             jaccard_top5pct=round(inter / max(union, 1), 4),
                             n_intersection=inter, n_union=union,
                             spearman_tierA=round(rho, 4),
                             pearson_tierA=round(rho_p, 4)))
        a.z.close(); b.z.close()
    P.CACHE3 = iso_cache
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False))
    print("wrote", out)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["headlines", "shift"])
    ap.add_argument("--arm", choices=["h1", "m1"])
    a = ap.parse_args()
    if a.stage == "headlines":
        stage_headlines()
    else:
        assert a.arm, "--arm is required for --stage shift"
        stage_shift(a.arm)
