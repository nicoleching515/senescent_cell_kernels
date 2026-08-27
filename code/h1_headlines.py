#!/usr/bin/env python3
"""Phase 10 -- the H1 analogue of `code/m1_headlines.py`.

Same formulas, same column names, same reportable-fit filter, so the two arms' headline
vectors are produced by the same arithmetic and can be put side by side in the §17 table.
Every formula below is copied from `m1_headlines.py` line for line; the only differences
are the section list, the results directory, and `--call`, because H1's primary sender call
is `tierAmg_p95` (PREREG D-B) and M1's is `tierA_p95`.

  python3 code/h1_headlines.py --call tierAmg_p95
  python3 code/h1_headlines.py --call tierA_p95        # the frozen-literal sensitivity

Writes results/phase10_h1/h1_headlines_<call>.json and prints it.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H

R3 = "/workspace/results/phase10_h1"
R5 = "/workspace/results/phase10_h1"
SECTIONS = list(H.ALL_SECTIONS)


def headlines(call, r3=R3, r5=R5, sections=None):
    sections = sections or SECTIONS
    out = {"arm": "h1", "call": call, "n_sections": len(sections)}
    mf = pd.read_csv(f"{r3}/main_fits.csv")
    mf = mf[mf.section.isin(sections) & (mf.call == call) & (mf.stratum == "all")]
    rep = mf[(mf.beta_naive > 0) & (mf.sf_base.notna()) & (mf.beta_base_lo > 0)]
    out["n_fits"] = int(len(mf))
    out["n_reportable"] = int(len(rep))
    out["frac_reportable"] = float(len(rep) / max(len(mf), 1))
    out["sender_prevalence_min"] = float(mf.prevalence.min())
    out["sender_prevalence_max"] = float(mf.prevalence.max())
    out["naive_amp_med"] = float((rep.beta_naive / rep.sd_y).median())
    out["naive_abs_amp_med_allfits"] = float((mf.beta_naive / mf.sd_y).abs().median())
    ctrl = rep.beta_n2n5n6 / rep.sd_y
    out["ctrl_amp_med"] = float(ctrl.median())
    out["ctrl_amp_iqr"] = [float(ctrl.quantile(.25)), float(ctrl.quantile(.75))]
    se = ((rep.beta_n2n5n6_hi - rep.beta_n2n5n6_lo) / (2 * 1.959964) / rep.sd_y)
    out["ctrl_amp_se_med"] = float(se.median())
    out["power80_bound"] = float(2.802 * se.median())          # 1.96 + 0.842
    for k, lab in [("sf_n2", "N2"), ("sf_n5", "N5"), ("sf_n6", "N6"),
                   ("sf_zon", "zon"), ("sf_n6n5", "N5+N6"),
                   ("sf_n2n5n6", "N2+N5+N6")]:
        out[f"SF_{lab}"] = float(rep[k].median())
        out[f"SF_{lab}_iqr"] = [float(rep[k].quantile(.25)), float(rep[k].quantile(.75))]
    out["lam_railed_frac"] = float(mf.lam_railed.mean())
    out["lam_railed_lo"] = int((mf.lam_naive <= mf.lam_grid_lo + 1e-9).sum())
    out["lam_railed_hi"] = int((mf.lam_naive >= mf.lam_grid_hi - 1e-9).sum())
    out["lam_naive_median"] = float(mf.lam_naive.median())
    out["ctrl_pos_and_sig"] = int(((rep.beta_n2n5n6 > 0) & (rep.beta_n2n5n6_lo > 0)).sum())
    key = rep[["section", "celltype", "module"]].drop_duplicates()

    def _perm(fn, cols, prefix=""):
        p = os.path.join(r3, fn)
        if not os.path.exists(p):
            return
        d = pd.read_csv(p)
        d = d[d.section.isin(sections) & (d.call == call)]
        if "scope" in d.columns and "full" in set(d.scope):
            d = d[d.scope == "full"]
        d = d.merge(key, on=["section", "celltype", "module"])
        for c in cols:
            if c in d.columns:
                v = d[c].dropna()
                if len(v):
                    out[prefix + c.replace("_sf", "")] = float(v.median())
                    out[prefix + c.replace("_sf", "") + "_iqr"] = [
                        float(v.quantile(.25)), float(v.quantile(.75))]
                    out[prefix + c.replace("_sf", "") + "_n"] = int(len(v))

    _perm("perm_nulls.csv", ["N1_sf", "N3_sf", "N4_sf", "N1_full_sf"])
    _perm("perm_nulls_c1.csv",
          [f"{nm}_sf" for nm in ("N3_orig", "N3_tile", "N3_occ", "N3_occ15", "N3_swap",
                                 "N3_snap", "N4_orig", "N4_tile", "N4_occ", "N4_occ15",
                                 "N4_swap")], prefix="c1_")
    # the tile scope is a DIFFERENT (smaller) fit population and gets its own rows
    p = os.path.join(r3, "perm_nulls_c1.csv")
    if os.path.exists(p):
        d = pd.read_csv(p)
        d = d[d.section.isin(sections) & (d.call == call) & (d.scope == "tile")]
        d = d.merge(key, on=["section", "celltype", "module"])
        for nm in ("N3_tile", "N4_tile"):
            c = f"{nm}_sf"
            if c in d.columns:
                v = d[c].dropna()
                if len(v):
                    out[f"c1_{nm}"] = float(v.median())
                    out[f"c1_{nm}_iqr"] = [float(v.quantile(.25)), float(v.quantile(.75))]
                    out[f"c1_{nm}_n"] = int(len(v))
    _perm("perm_nulls_var.csv", ["N3_var_sf", "N4_var_sf"], prefix="var_")

    p = f"{r3}/poisson_fits.csv"
    if os.path.exists(p):
        d = pd.read_csv(p)
        for sub, tag in [("ALL sections", ""), ("tierA percentile", "_tierA")]:
            r = d[d.subset.str.contains(sub, na=False)]
            if len(r):
                out["poisson_slope" + tag] = float(r.slope.iloc[0])
                out["poisson_r2" + tag] = float(r.r2.iloc[0])
    p = f"{r5}/kernel_families.csv"
    if os.path.exists(p):
        d = pd.read_csv(p)
        d = d[d.call == call] if "call" in d.columns else d
        for dg in sorted(set(d.design.dropna())):
            s = d[d.design == dg]
            wins = (s.family == s.best_family).groupby(s.family).mean()
            out[f"kernel_win_{dg}"] = {k: round(float(v), 4) for k, v in wins.items()}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--call", default="tierAmg_p95")
    a = ap.parse_args()
    o = headlines(a.call)
    with open(f"{R3}/h1_headlines_{a.call}.json", "w") as fh:
        json.dump(o, fh, indent=1, default=float)
    print(json.dumps(o, indent=1, default=float))
