#!/usr/bin/env python3
"""Phase 8 — the variance-corrected null beside every C1 variant.

Writes NEW files only:

    results/phase3/sf_summary_var.csv        the comparison table
    results/phase3/summary_phase3_var.txt    the same, rendered
    results/phase3/var_pvalues.csv           the Mrkvicka et al. Monte Carlo
                                             test outcomes (RS_count, RS_ker)

Reportable population is EXACTLY the one `summarize_phase3.py` and
`summarize_phase3_c1.py` use -- imported from `summarize_phase3_c1`, not
re-implemented -- so every row of the table is over the same fits.
That population is derived from `main_fits.csv`, which the M1 re-run
regenerates, so the `n` column is printed rather than hard-coded.
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
import summarize_phase3_c1 as SC1

RES = P.RESULTS
KEY = ["section", "celltype", "module"]

LABEL = dict(SC1.LABEL)
LABEL["N3_var"] = ("N3-var VARIANCE-CORRECTED Euclidean shift "
                   "(Mrkvicka et al. 2021)")
LABEL["N4_var"] = ("N4-var VARIANCE-CORRECTED rotation "
                   "(this work, same principle)")

ORDER = ["N3_orig", "N3_orig_rerun", "N3_tile", "N3_occ", "N3_occ15",
         "N3_swap", "N3_snap", "N3_var", "N4_orig", "N4_orig_rerun",
         "N4_tile", "N4_occ", "N4_occ15", "N4_swap", "N4_var"]
# the C1 runner re-runs the ORIGINAL bounding-box nulls in the same job, so
# perm_nulls_c1.csv carries its own N3_orig/N4_orig.  They are kept under a
# distinct name so the published row and the re-run row do not collide.
RERUN = {"N3_orig": "N3_orig_rerun", "N4_orig": "N4_orig_rerun"}
for _k, _v in RERUN.items():
    LABEL[_v] = LABEL[_k].replace("(ORIGINAL)", "(ORIGINAL, re-run in the C1 job)")


def _load_destructiveness():
    a = pd.read_csv(f"{RES}/null_destructiveness.csv")
    b = pd.read_csv(f"{RES}/null_destructiveness_var.csv")
    cols = ["frac_retaining_a_neighbour", "real_median_nbrs",
            "null_median_nbrs", "median_displacement_um"]
    d = pd.concat([a[["null"] + cols], b[["null"] + cols]], ignore_index=True)
    return d.groupby("null", sort=False)[cols].median()


def main(var_file="perm_nulls_var.csv", full_file="perm_nulls_var_full200.csv"):
    rep = SC1.reportable()
    orig = pd.read_csv(f"{RES}/perm_nulls.csv")
    orig = orig[orig.call == RN.PRIMARY_CALL]
    c1 = pd.read_csv(f"{RES}/perm_nulls_c1.csv")
    va = pd.read_csv(f"{RES}/{var_file}")
    vf = (pd.read_csv(f"{RES}/{full_file}")
          if os.path.exists(f"{RES}/{full_file}") else None)
    keep = _load_destructiveness()

    rows = []
    # --- the ORIGINAL published bounding-box torus shift / rotation --------
    m0 = rep.merge(orig, on=KEY, how="left")
    for nm, col in (("N3_orig", "N3"), ("N4_orig", "N4")):
        s = SC1._stats(m0[f"{col}_sf"], m0[f"{col}_p"])
        fs = SC1._stats(m0[f"{col}_full_sf"])
        rows.append(dict(variant=nm, label=LABEL[nm], scope="whole section",
                         source="perm_nulls.csv (ORIGINAL, published)",
                         n_perm=int(m0.n_perm.median()), **s,
                         full_sf_median=fs["median"], full_sf_n=fs["n"]))

    # --- the C1 in-tissue variants ----------------------------------------
    for scope in ("full", "tile"):
        sub = c1[c1.scope == scope]
        if sub.empty:
            continue
        m = rep.merge(sub, on=KEY, how="inner")
        names = GEO.TILE_NULLS if scope == "tile" else GEO.FULL_NULLS
        for nm in names:
            if f"{nm}_sf" not in m:
                continue
            s = SC1._stats(m[f"{nm}_sf"], m.get(f"{nm}_p"))
            if s is None:
                continue
            fs = SC1._stats(m[f"{nm}_full_sf"]) if f"{nm}_full_sf" in m else None
            out_nm = RERUN.get(nm, nm)
            rows.append(dict(
                variant=out_nm, label=LABEL[out_nm],
                scope="solid tiles" if scope == "tile" else "whole section",
                source="perm_nulls_c1.csv (C1 re-run)",
                n_perm=int(m.n_perm.median()), **s,
                full_sf_median=fs["median"] if fs else np.nan,
                full_sf_n=fs["n"] if fs else np.nan))

    # --- the variance-corrected nulls -------------------------------------
    mv = rep.merge(va, on=KEY, how="inner")
    mvf = rep.merge(vf, on=KEY, how="inner") if vf is not None else None
    for nm in ("N3_var", "N4_var"):
        s = SC1._stats(mv[f"{nm}_sf"], mv[f"{nm}_cov_p_rscount"])
        if s is None:
            continue
        # window-matched SF as a RATIO OF MEANS, which is the convention the
        # tile scope already uses (beta_obs re-profiled on the restricted cell
        # set).  The per-draw mean-of-ratios column `{nm}_sf_wm` in
        # perm_nulls_var.csv is unstable when a single draw's beta_obs|W_i is
        # near zero and is kept only for audit.
        wm_ratio = ((mv[f"{nm}_beta_obs_wm"] - mv[f"{nm}_null_mean"])
                    / mv[f"{nm}_beta_obs_wm"])
        wm = SC1._stats(wm_ratio)
        fs = (SC1._stats(mvf[f"{nm}_full_sf"])
              if mvf is not None and f"{nm}_full_sf" in mvf else None)
        rows.append(dict(
            variant=nm, label=LABEL[nm],
            scope="whole section (tissue window W)",
            source=f"{var_file} (variance correction)",
            n_perm=int(mv.n_perm.median()), **s,
            full_sf_median=fs["median"] if fs else np.nan,
            full_sf_n=fs["n"] if fs else np.nan,
            sf_wm_median=wm["median"], sf_wm_q25=wm["q25"],
            sf_wm_q75=wm["q75"],
            frac_cells_retained=float(mv[f"{nm}_frac_cells_retained"].median()),
            full_sf_n_perm=(int(mvf.n_perm.median())
                            if mvf is not None else np.nan)))

    df = pd.DataFrame(rows)
    inv = {v: k for k, v in RERUN.items()}
    for c in ("frac_retaining_a_neighbour", "real_median_nbrs",
              "null_median_nbrs", "median_displacement_um"):
        df[c] = df.variant.map(lambda v: keep[c].get(inv.get(v, v), np.nan))
    df["variant"] = pd.Categorical(df.variant, ORDER, ordered=True)
    df = df.sort_values("variant").reset_index(drop=True)
    df.insert(0, "subset", "PRIMARY: in-band sections, tierA_p95, "
                           "A_SENDER_FINAL_strict (33 genes)")
    df.to_csv(f"{RES}/sf_summary_var.csv", index=False)

    # --- the Monte Carlo test outcomes, per fit ---------------------------
    pcols = [c for c in mv.columns
             if any(c.endswith(x) for x in
                    ("_p_rscount", "_p_rscount_1s", "_p_raw", "_p_naive"))
             or "_p_rsker" in c]
    mv[KEY + ["n", "lam", "beta_obs"] + pcols].to_csv(
        f"{RES}/var_pvalues.csv", index=False)

    w = ["=" * 124,
         "PHASE 8 — SURVIVING FRACTION UNDER EVERY N3/N4 VARIANT, "
         "INCLUDING THE MRKVICKA ET AL. (2021) VARIANCE CORRECTION",
         "=" * 124,
         f"{'variant':64s} {'n':>4s} {'medSF':>7s} {'IQR':>17s} "
         f"{'<=0':>5s} {'rej':>5s} {'fullSF':>7s} {'keep':>5s} {'disp':>8s}"]
    for _, r in df.iterrows():
        w.append(f"{r.label:64s} {r.n:4d} {r['median']:7.3f} "
                 f"[{r.q25:7.3f},{r.q75:7.3f}] {r.frac_le_0:5.2f} "
                 f"{r.get('reject_rate_p05', np.nan):5.2f} "
                 f"{r.get('full_sf_median', np.nan):7.3f} "
                 f"{r.frac_retaining_a_neighbour:5.3f} "
                 f"{r.median_displacement_um:6.0f}um")
    w += ["",
          "medSF = median surviving fraction (beta_obs - mean null beta)/beta_obs "
          "over the reportable fits (n column)",
          "rej   = rejection rate at p<0.05; for N3-var/N4-var this is the "
          "RS_count-standardized Monte Carlo test of Mrkvicka et al. (2021),",
          "        for every other row it is the repo's uncorrected "
          "permutation p",
          "keep  = fraction of shifted senders retaining a real cell within "
          "100 um (median over sections); for the var rows, of the RETAINED",
          "        shifted senders (those landing inside the tissue window W)",
          "disp  = median distance a sender is actually moved"]
    for _, r in df[df.variant.astype(str).str.endswith("_var")].iterrows():
        w.append(f"{r.variant}: window-matched SF median "
                 f"{r.sf_wm_median:.3f} [{r.sf_wm_q25:.3f},{r.sf_wm_q75:.3f}]; "
                 f"median fraction of receiver cells retained in W_i "
                 f"{r.frac_cells_retained:.3f}")
    txt = "\n".join(w)
    print(txt)
    with open(f"{RES}/summary_phase3_var.txt", "w") as fh:
        fh.write(txt + "\n")
    print(f"\nwrote {RES}/sf_summary_var.csv, summary_phase3_var.txt, "
          f"var_pvalues.csv")
    return df


if __name__ == "__main__":
    main()
