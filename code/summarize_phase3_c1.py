#!/usr/bin/env python3
"""Phase 7 / C1 — old vs new surviving fractions for N3 and N4.

Nothing here overwrites `sf_summary.csv`.  It writes
`results/phase3/sf_summary_c1.csv`, which carries the ORIGINAL
bounding-box N3/N4 rows (as published, from `perm_nulls.csv`) beside every
corrected in-tissue variant (from `perm_nulls_c1.csv`), on the same
reportable population.

Reportable population = exactly the one `summarize_phase3.py` uses: the six
in-band sections, the primary sender call, receiver-type x module fits whose
naive amplitude is positive and whose spatial block bootstrap CI excludes zero
(160 of 315).  The tile-scope variants are reported on the subset of those same
160 (section, receiver type, module) cells that still clears the 2,000-receiver
floor once the fit is restricted to solid tissue tiles.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN
import phase3_null_geom as GEO

RES = P.RESULTS
KEY = ["section", "celltype", "module"]

LABEL = {
    "N3_orig": "N3 torus shift, whole-section bounding box (ORIGINAL)",
    "N3_tile": "N3-tile torus shift inside solid-tissue tiles",
    "N3_occ": "N3-occ whole-section shift, <=5% of senders out of tissue",
    "N3_occ15": "N3-occ15 whole-section shift, <=15% out of tissue (suppl.)",
    "N3_swap": "N3-swap senders relocated to random real cell positions",
    "N3_snap": "N3-snap whole-section shift, snapped to nearest real cell (suppl.)",
    "N4_orig": "N4 rotation, whole-section bounding box (ORIGINAL)",
    "N4_tile": "N4-tile rotation inside solid-tissue tiles",
    "N4_occ": "N4-occ whole-section rotation, <=5% of senders out of tissue",
    "N4_occ15": "N4-occ15 whole-section rotation, <=15% out of tissue (suppl.)",
    "N4_swap": "N4-swap rotation, snapped to nearest real cell positions",
}


def reportable(call=RN.PRIMARY_CALL):
    """The population summarize_phase3.py reports: in-band, this sender call,
    positive naive amplitude whose spatial block bootstrap CI excludes zero."""
    mf = pd.read_csv(f"{RES}/main_fits.csv")
    return mf[(mf.call == call) & (mf.stratum == "all")
              & mf.section.isin(P.IN_BAND)
              & (mf.beta_naive > 0) & (mf.beta_base_lo > 0)][KEY].copy()


def _stats(v, p=None):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    q25, med, q75 = np.quantile(v, [.25, .5, .75])
    out = dict(n=int(v.size), q25=q25, median=med, q75=q75,
               mean=float(v.mean()),
               frac_le_0=float((v <= 0).mean()),
               frac_gt_05=float((v > .5).mean()))
    if p is not None:
        p = np.asarray(p, float)
        p = p[np.isfinite(p)]
        out["reject_rate_p05"] = float((p < 0.05).mean()) if p.size else np.nan
    return out


def n7_table():
    """The N7 sender-definition axis under the corrected nulls.

    N7 varies the sender DEFINITION, so it is the axis that tells us whether
    the C1 correction behaves the same way for every sender call -- which is
    what the Phase 7 freeze needs, because the Phase 7 primary sender
    definition is not the Phase 3 one.
    """
    # Phase 8 / 8.7: the axis now also carries the second pre-registered Tier A
    # variant (`tierApm_p95`, the per-module sender sets), whose files take the
    # `_pm` suffix because its sender mask is module-specific.
    parts, oparts = [], []
    for suf in ("_n7", "_pm"):
        f = f"{RES}/perm_nulls_c1{suf}.csv"
        if os.path.exists(f):
            parts.append(pd.read_csv(f))
        g = f"{RES}/perm_nulls{suf}.csv"
        if os.path.exists(g):
            oparts.append(pd.read_csv(g))
    if not parts:
        return pd.DataFrame()
    c1 = pd.concat(parts, ignore_index=True)
    o7 = (pd.concat(oparts, ignore_index=True) if oparts
          else pd.DataFrame(columns=["call"]))
    rows = []
    for call in sorted(c1.call.unique()):
        rep = reportable(call)
        oo = (rep.merge(o7[o7.call == call], on=KEY, how="inner")
              if len(o7) else pd.DataFrame())
        for nm, col in (("N3_orig", "N3"), ("N4_orig", "N4")):
            if oo.empty or f"{col}_sf" not in oo:
                continue
            st = _stats(oo[f"{col}_sf"], oo[f"{col}_p"])
            if st:
                rows.append(dict(call=call, variant=nm, label=LABEL[nm],
                                 scope="whole section",
                                 source="perm_nulls_n7.csv (ORIGINAL)", **st))
        for scope in ("full", "tile"):
            sub = c1[(c1.call == call) & (c1.scope == scope)]
            if sub.empty:
                continue
            m = rep.merge(sub, on=KEY, how="inner")
            names = GEO.TILE_NULLS if scope == "tile" else GEO.FULL_NULLS
            for nm in names:
                if f"{nm}_sf" not in m:
                    continue
                st = _stats(m[f"{nm}_sf"], m.get(f"{nm}_p"))
                if st:
                    rows.append(dict(
                        call=call, variant=nm, label=LABEL[nm],
                        scope="solid tiles" if scope == "tile" else "whole section",
                        source="perm_nulls_c1_n7.csv (C1 re-run)", **st))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/sf_summary_c1_n7.csv", index=False)
    return df


def main():
    rep = reportable()
    orig = pd.read_csv(f"{RES}/perm_nulls.csv")
    orig = orig[orig.call == RN.PRIMARY_CALL]
    c1 = pd.read_csv(f"{RES}/perm_nulls_c1.csv")
    dz = pd.read_csv(f"{RES}/null_destructiveness.csv")
    keep = dz.groupby("null", sort=False)[
        ["frac_retaining_a_neighbour", "real_median_nbrs", "null_median_nbrs",
         "median_displacement_um"]].median()

    rows = []
    m0 = rep.merge(orig, on=KEY, how="left")
    for nm, col in (("N3_orig", "N3"), ("N4_orig", "N4")):
        s = _stats(m0[f"{col}_sf"], m0[f"{col}_p"])
        fs = _stats(m0[f"{col}_full_sf"]) if f"{col}_full_sf" in m0 else None
        rows.append(dict(variant=nm, label=LABEL[nm], scope="whole section",
                         source="perm_nulls.csv (ORIGINAL, published)",
                         n_perm=int(m0.n_perm.median()), **s,
                         full_sf_median=fs["median"] if fs else np.nan,
                         full_sf_q25=fs["q25"] if fs else np.nan,
                         full_sf_q75=fs["q75"] if fs else np.nan,
                         full_sf_n=fs["n"] if fs else np.nan))

    for scope in ("full", "tile"):
        sub = c1[c1.scope == scope]
        if sub.empty:
            continue
        m = rep.merge(sub, on=KEY, how="inner")
        names = GEO.TILE_NULLS if scope == "tile" else GEO.FULL_NULLS
        for nm in names:
            if f"{nm}_sf" not in m:
                continue
            s = _stats(m[f"{nm}_sf"], m.get(f"{nm}_p"))
            if s is None:
                continue
            fs = _stats(m[f"{nm}_full_sf"]) if f"{nm}_full_sf" in m else None
            rows.append(dict(
                variant=nm, label=LABEL[nm],
                scope="solid tiles" if scope == "tile" else "whole section",
                source="perm_nulls_c1.csv (C1 re-run)",
                n_perm=int(m.n_perm.median()), **s,
                full_sf_median=fs["median"] if fs else np.nan,
                full_sf_q25=fs["q25"] if fs else np.nan,
                full_sf_q75=fs["q75"] if fs else np.nan,
                full_sf_n=fs["n"] if fs else np.nan))

    df = pd.DataFrame(rows)
    for c in ("frac_retaining_a_neighbour", "real_median_nbrs",
              "null_median_nbrs", "median_displacement_um"):
        df[c] = df.variant.map(keep[c])
    df.insert(0, "subset", "PRIMARY: in-band sections, tierA_p95")
    df.to_csv(f"{RES}/sf_summary_c1.csv", index=False)

    w = ["=" * 118,
         "PHASE 7 / C1 — N3 AND N4 SURVIVING FRACTIONS, ORIGINAL vs IN-TISSUE",
         "=" * 118,
         f"{'variant':62s} {'n':>4s} {'medSF':>7s} {'IQR':>17s} "
         f"{'<=0':>5s} {'rej':>5s} {'fullSF':>7s} {'keep':>5s} {'disp':>7s}"]
    for _, r in df.iterrows():
        w.append(f"{r.label:62s} {r.n:4d} {r['median']:7.3f} "
                 f"[{r.q25:7.3f},{r.q75:7.3f}] {r.frac_le_0:5.2f} "
                 f"{r.get('reject_rate_p05', np.nan):5.2f} "
                 f"{r.get('full_sf_median', np.nan):7.3f} "
                 f"{r.frac_retaining_a_neighbour:5.3f} "
                 f"{r.median_displacement_um:6.0f}um")
    w.append("")
    w.append("fullSF = the same surviving fraction with beta taken under the "
             "FULL N5+N6+zonation design at fixed lambda")
    w.append("keep = fraction of shifted senders retaining a real neighbour "
             "within 100 um (median over sections)")
    w.append("disp = median distance a sender is actually moved by the null "
             "(median over sections)")
    n7 = n7_table()
    if not n7.empty:
        w.append("")
        w.append("=" * 118)
        w.append("N7 SENDER-DEFINITION AXIS UNDER THE CORRECTED NULLS "
                 "(200 permutations, in-band sections)")
        w.append("=" * 118)
        w.append(f"{'call':12s} {'variant':10s} {'n':>4s} {'medSF':>7s} "
                 f"{'IQR':>17s} {'<=0':>5s} {'rej':>5s}")
        for _, r in n7.iterrows():
            w.append(f"{r.call:12s} {r.variant:10s} {r.n:4d} "
                     f"{r['median']:7.3f} [{r.q25:7.3f},{r.q75:7.3f}] "
                     f"{r.frac_le_0:5.2f} "
                     f"{r.get('reject_rate_p05', np.nan):5.2f}")
    txt = "\n".join(w)
    print(txt)
    with open(f"{RES}/summary_phase3_c1.txt", "w") as fh:
        fh.write(txt + "\n")
    print(f"\nwrote {RES}/sf_summary_c1.csv, sf_summary_c1_n7.csv "
          f"and summary_phase3_c1.txt")
    return df


if __name__ == "__main__":
    main()
