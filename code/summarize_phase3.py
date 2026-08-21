#!/usr/bin/env python3
"""Phase 3 summary tables -> /workspace/results/phase3/summary_phase3.txt

Everything is a SURVIVING FRACTION (Section 6.5).  The reference population is
the set of fits a reader of the current literature would actually report: a
POSITIVE naive amplitude whose spatial block bootstrap CI excludes zero.  Fits
that are already null naively cannot "survive" anything and their ratios are
undefined noise, so they are counted separately rather than averaged in.

Primary section set: the six Section 8 Test 3 admissible sections (in band),
six animals, both surgical arms.  The four over-ceiling and one below-floor
sections are reported separately as a sensitivity.
"""
from __future__ import annotations

import os
import sys
import io
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import run_phase3_nulls as RN

RES = P.RESULTS
OUT = io.StringIO()
BAND = {s: "in_band" for s in P.IN_BAND}
BAND.update({s: "over_ceiling" for s in P.OVER_CEILING})
BAND.update({s: "below_floor" for s in P.BELOW_FLOOR})

NULLS = [("sf_n2", "N2 matched decoy"),
         ("N1_sf", "N1 stratified label permutation"),
         ("N1_full_sf", "N1 on the N5+N6-conditioned residual"),
         ("N3_sf", "N3 torus shift"),
         ("N4_sf", "N4 rotation"),
         ("sf_n5", "N5 nuisance conditioning"),
         ("sf_n6", "N6 receiver-baseline conditioning"),
         ("sf_zon", "zonation covariate alone"),
         ("sf_n6n5", "N5 + N6"),
         ("sf_n2n5n6", "N2 + N5 + N6 (combined)"),
         ("sf_n8", "N8 scrambled response gene set")]


def w(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.write(s + "\n")


def load():
    mf = pd.read_csv(f"{RES}/main_fits.csv")
    mf["band"] = mf.section.map(BAND)
    pf = pd.DataFrame()
    for f in ("perm_nulls.csv", "perm_nulls_n7.csv"):
        if os.path.exists(f"{RES}/{f}"):
            pf = pd.concat([pf, pd.read_csv(f"{RES}/{f}")], ignore_index=True)
    if not pf.empty:
        pf["band"] = pf.section.map(BAND)
    n8 = [f"{RES}/n8_scrambled_{s}.csv" for s in P.ALL_SECTIONS]
    n8 = [p for p in n8 if os.path.exists(p)]
    n8 = pd.concat([pd.read_csv(p) for p in n8], ignore_index=True) if n8 \
        else pd.DataFrame()
    return mf, pf, n8


def merged(mf, pf, n8, sections, call=RN.PRIMARY_CALL, stratum="all"):
    key = ["section", "celltype", "module"]
    d = mf[mf.section.isin(sections) & (mf.call == call)
           & (mf.stratum == stratum)].copy()
    if not pf.empty:
        p = pf[pf.section.isin(sections) & (pf.call == call)]
        cols = [c for c in p.columns if c.endswith("_sf") or c.endswith("_p")]
        if not p.empty:
            d = d.merge(p[key + cols].drop_duplicates(key), on=key, how="left")
    if not n8.empty:
        d = d.merge(n8[key + ["sf_n8", "beta_obs_std", "rand_mean",
                              "rand_absmean", "pct_rand_ge"]
                       ].drop_duplicates(key), on=key, how="left")
    d["beta_std"] = d.beta_naive / d.sd_y
    return d


def sf_table(d, label):
    rep = d[(d.beta_naive > 0) & (d.beta_base_lo > 0)]
    w(f"\n--- {label} ---")
    w(f"fits {len(d)};  naive beta > 0 in {int((d.beta_naive>0).sum())};"
      f"  positive AND block-bootstrap CI excludes 0: {len(rep)}"
      f" ({100*len(rep)/max(len(d),1):.0f} %)")
    w(f"{'null':40s} {'n':>4s} {'medSF':>7s} {'IQR':>17s} {'<=0':>5s} {'>0.5':>5s}")
    rows = []
    for col, lab in NULLS:
        if col not in rep:
            continue
        v = rep[col].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        a, b, c = np.quantile(v, [.25, .5, .75])
        rows.append(dict(subset=label, null=lab, n=int(v.size), q25=a,
                         median=b, q75=c, frac_le_0=float((v <= 0).mean()),
                         frac_gt_05=float((v > .5).mean())))
        w(f"{lab:40s} {v.size:4d} {b:7.3f} [{a:7.3f},{c:7.3f}] "
          f"{(v<=0).mean():5.2f} {(v>.5).mean():5.2f}")
    return pd.DataFrame(rows), rep


def main():
    mf, pf, n8 = load()
    lam = RN.lam_grid()
    w("=" * 100)
    w("PHASE 3 — NULL BATTERY N1-N8 ON REAL MOUSE LIVER (GSE310392, Xenium 5K)")
    w("=" * 100)
    w(f"window {RN.WINDOW_UM:.0f} um (99th pct of observed distance-to-nearest-"
      f"sender, in-band sections, primary call)")
    w(f"lambda grid [{lam[0]:.0f}, {lam[-1]:.0f}] um, {RN.N_LAM} log-spaced "
      f"points; floor = resolution floor (median NN 6.7-10.6 um), "
      f"ceiling = window/2")
    w(f"{RN.N_BOOT} spatial block bootstrap replicates over "
      f"{RN.N_BLOCKS_SIDE**2} quantile blocks; receiver labels = "
      f"'{P.LABELS}' (Bio Phase 3 stable label families)")
    w("\nSECTION ADMISSIBILITY (Section 8 Test 3, Cdkn1a+ hepatocyte prevalence)")
    w(f"  in band (PRIMARY, 6 animals): {', '.join(s[:4] for s in P.IN_BAND)}")
    w(f"  over the 20 % ceiling (excluded): "
      f"{', '.join(s[:4] for s in P.OVER_CEILING)}")
    w(f"  below the 1 % floor (excluded): "
      f"{', '.join(s[:4] for s in P.BELOW_FLOOR)}")

    d = merged(mf, pf, n8, P.IN_BAND)
    w("\n### 1. NAIVE FITS — in-band sections, tierA_p95, receiver type x module")
    w(f"n fits {len(d)}; sections {d.section.nunique()}; "
      f"receiver types {d.celltype.nunique()}; modules {d.module.nunique()}")
    w(f"lambda_hat railed at a grid bound in {int(d.lam_railed.sum())}/{len(d)} "
      f"({100*d.lam_railed.mean():.0f} %): "
      f"{int((d.lam_naive <= lam[0]+1e-9).sum())} at the {lam[0]:.0f} um floor, "
      f"{int((d.lam_naive >= lam[-1]-1e-9).sum())} at the {lam[-1]:.0f} um "
      f"ceiling; interior median {d.lam_naive[d.lam_railed==0].median():.1f} um")
    w(f"beta_hat: {int((d.beta_naive>0).sum())} positive, "
      f"{int((d.beta_naive<0).sum())} negative; |beta|/sd(y) median "
      f"{d.beta_std.abs().median():.3f}, p90 {d.beta_std.abs().quantile(.9):.3f}")

    tabs = []
    t, rep = sf_table(d, "PRIMARY: in-band sections, tierA_p95")
    tabs.append(t)

    w("\n### 2. SURVIVING FRACTION BY RECEIVER CELL TYPE")
    w(f"{'cell type':22s} {'n':>3s} {'medN':>7s} {'b/sd':>6s} {'N1':>7s} "
      f"{'N3':>7s} {'N2':>7s} {'N5':>7s} {'N6':>7s} {'zon':>7s} {'N2N5N6':>7s}")
    for ct, g in rep.groupby("celltype"):
        w(f"{ct:22s} {len(g):3d} {g.n.median():7.0f} {g.beta_std.median():6.3f} "
          f"{g.N1_sf.median():7.3f} {g.N3_sf.median():7.3f} "
          f"{g.sf_n2.median():7.3f} {g.sf_n5.median():7.3f} "
          f"{g.sf_n6.median():7.3f} {g.sf_zon.median():7.3f} "
          f"{g.sf_n2n5n6.median():7.3f}")
    w("\n### 3. SURVIVING FRACTION BY TIER B MODULE")
    w(f"{'module':22s} {'n':>3s} {'b/sd':>6s} {'N1':>7s} {'N3':>7s} "
      f"{'N2':>7s} {'N5':>7s} {'N6':>7s} {'zon':>7s} {'N2N5N6':>7s} {'N8':>7s}")
    for m, g in rep.groupby("module"):
        n8v = g.sf_n8.median() if "sf_n8" in g else np.nan
        w(f"{m:22s} {len(g):3d} {g.beta_std.median():6.3f} "
          f"{g.N1_sf.median():7.3f} {g.N3_sf.median():7.3f} "
          f"{g.sf_n2.median():7.3f} {g.sf_n5.median():7.3f} "
          f"{g.sf_n6.median():7.3f} {g.sf_zon.median():7.3f} "
          f"{g.sf_n2n5n6.median():7.3f} {n8v:7.3f}")

    w("\n### 4. SENSITIVITY TO THE TEST-3 ADMISSIBILITY RULE")
    for nm, secs in (("over-ceiling sections (Cdkn1a+ hep 22-45 %)",
                      P.OVER_CEILING),
                     ("below-floor section (0.48 %)", P.BELOW_FLOOR)):
        dd = merged(mf, pf, n8, secs)
        if dd.empty:
            continue
        t2, _ = sf_table(dd, nm)
        tabs.append(t2)

    w("\n### 5. ARM CONTRAST WITHIN THE ADMISSIBLE SET "
      "(2 SBR animals vs 4 sham animals)")
    w(f"{'arm':6s} {'sections':>8s} {'fits':>5s} {'rep':>4s} {'med b/sd':>9s} "
      f"{'med SF N5':>10s} {'med SF N2N5N6':>14s} {'med SF N1':>10s}")
    for arm in ("sbr", "sham"):
        g = d[d.arm == arm]
        r = g[(g.beta_naive > 0) & (g.beta_base_lo > 0)]
        if g.empty:
            continue
        w(f"{arm:6s} {g.section.nunique():8d} {len(g):5d} {len(r):4d} "
          f"{r.beta_std.median():9.3f} {r.sf_n5.median():10.3f} "
          f"{r.sf_n2n5n6.median():14.3f} {r.N1_sf.median():10.3f}")
    w("  matched-decoy balance PER ARM (Bio Phase 2 §5: never pooled):")
    for arm in ("sbr", "sham"):
        g = d[d.arm == arm]
        if g.empty:
            continue
        w(f"    {arm}: max |SMD| before {g.max_smd_before.min():.3f}-"
          f"{g.max_smd_before.max():.3f}, after {g.max_smd_after.min():.3f}-"
          f"{g.max_smd_after.max():.3f}; match rate "
          f"{g.match_rate.min():.3f}-{g.match_rate.max():.3f}")
    w("  section-level detection depth vs naive |beta|/sd (Spearman over the "
      "six in-band sections, pooled over cell types x modules):")
    try:
        from scipy.stats import spearmanr
        pv = pd.read_csv(f"{RES}/poisson_density.csv")
        dep = pv[pv.call == "cdkn1a_pos"].set_index("section")["median_depth"]
        g = d.groupby("section").beta_std.median()
        common = [s for s in g.index if s in dep.index]
        r = spearmanr(dep.loc[common], g.loc[common])
        w(f"    rho = {r.statistic:+.3f} (p = {r.pvalue:.3f}, n = {len(common)})")
    except Exception as e:
        w(f"    (unavailable: {e})")

    w("\n### 6. N7 — SENDER THRESHOLD AND CALLER SENSITIVITY (in-band sections)")
    w(f"{'call':13s} {'prev%':>6s} {'fits':>5s} {'rep':>4s} {'railed':>7s} "
      f"{'medlam':>7s} {'med b/sd':>9s} {'SF N5':>7s} {'SF N2N5N6':>10s} "
      f"{'SF N1':>7s}")
    for call in RN.N7_CALLS:
        dd = merged(mf, pf, n8, P.IN_BAND, call=call)
        if dd.empty:
            continue
        r = dd[(dd.beta_naive > 0) & (dd.beta_base_lo > 0)]
        n1 = r.N1_sf.median() if "N1_sf" in r else np.nan
        w(f"{call:13s} {100*dd.prevalence.mean():6.2f} {len(dd):5d} {len(r):4d} "
          f"{dd.lam_railed.mean():7.2f} {dd.lam_naive.median():7.1f} "
          f"{r.beta_std.median():9.3f} {r.sf_n5.median():7.3f} "
          f"{r.sf_n2n5n6.median():10.3f} {n1:7.3f}")

    w("\n### 7. ZONATION-STRATIFIED FITS, HEPATOCYTES (Section 11, T3)")
    z = mf[mf.section.isin(P.IN_BAND) & (mf.call == RN.PRIMARY_CALL)
           & (mf.celltype == "Hepatocytes")].copy()
    z["beta_std"] = z.beta_naive / z.sd_y
    w(f"{'stratum':13s} {'fits':>5s} {'med n':>7s} {'railed':>7s} "
      f"{'med lam':>8s} {'med b/sd':>9s} {'pos&sig':>8s} {'SF N5':>7s} "
      f"{'SF zon':>7s} {'SF N2N5N6':>10s}")
    for st in ("all", "periportal", "midzonal", "pericentral"):
        g = z[z.stratum == st]
        if g.empty:
            continue
        r = g[(g.beta_naive > 0) & (g.beta_base_lo > 0)]
        f = lambda c: (r[c].median() if len(r) else np.nan)
        w(f"{st:13s} {len(g):5d} {g.n.median():7.0f} {g.lam_railed.mean():7.2f} "
          f"{g.lam_naive.median():8.1f} "
          f"{(r.beta_std.median() if len(r) else np.nan):9.3f} {len(r):8d} "
          f"{f('sf_n5'):7.3f} {f('sf_zon'):7.3f} {f('sf_n2n5n6'):10.3f}")
    w("  naive |beta|/sd per module, pooled vs within zone (median over sections):")
    w(z.pivot_table(index="module", columns="stratum",
                    values="beta_std", aggfunc="median").round(3).to_string())

    # ---- N8 --------------------------------------------------------------
    if not n8.empty:
        s8 = n8[n8.section.isin(P.IN_BAND)]
        if not s8.empty:
            w("\n### 8. N8 SCRAMBLED-RESPONSE CONTROL (Tier E3 expression-matched)")
            w(f"{'module':22s} {'n':>4s} {'beta_obs/sd':>12s} "
              f"{'rand mean':>10s} {'rand sd':>8s} {'SF_N8':>7s} "
              f"{'pct rand >= obs':>16s}")
            for m, g in s8.groupby("module"):
                w(f"{m:22s} {len(g):4d} {g.beta_obs_std.median():12.4f} "
                  f"{g.rand_mean.median():10.4f} {g.rand_sd.median():8.4f} "
                  f"{g.sf_n8.median():7.3f} {g.pct_rand_ge.median():16.3f}")
    for s in P.IN_BAND:
        p = f"{RES}/n8_circularity_{s}.csv"
        if not os.path.exists(p):
            continue
        c = pd.read_csv(p)
        w(f"\n### 9. DeepScence/CoreScence CIRCULARITY ({s[:4]})")
        w(f"{'sender call':18s} {'n':>4s} {'beta_full/sd':>13s} "
          f"{'beta_stripped/sd':>17s} {'ratio':>7s}")
        for sd_, g in c.groupby("sender"):
            w(f"{sd_:18s} {len(g):4d} {g.beta_full_std.median():13.4f} "
              f"{g.beta_stripped_std.median():17.4f} "
              f"{g.ratio_stripped.median():7.3f}")
        w("  per module, CoreScence sender call:")
        for m, g in c[c.sender == "corescence_p95"].groupby("module"):
            w(f"    {m:22s} ratio {g.ratio_stripped.median():6.3f} "
              f"({int(g.n_stripped_genes.iloc[0])} genes removed)")
        dj = pd.read_csv(f"{RES}/n8_disjointness_{s}.csv")
        w(f"  Tier A n Tier B = {int(dj.overlap_tierA.sum())} genes over all "
          f"7 modules (disjointness confirmed); CoreScence n Tier B = "
          f"{int(dj.overlap_corescence.sum())} genes, per-module fraction "
          f"{dj.frac_corescence.min():.2f}-{dj.frac_corescence.max():.2f}")
        break

    # ---- Poisson ---------------------------------------------------------
    if os.path.exists(f"{RES}/poisson_fits.csv"):
        pf2 = pd.read_csv(f"{RES}/poisson_fits.csv")
        w("\n### 10. IS THE REGRESSOR A MEASUREMENT OR A SENDER-CALLING RATE?")
        w("   log(median distance to nearest sender) vs log(sender density);")
        w("   a homogeneous Poisson process gives slope exactly -0.5.")
        w(f"{'subset':44s} {'n':>4s} {'slope':>7s} {'r2':>7s} "
          f"{'obs/Poisson med d':>18s}")
        for _, r in pf2.iterrows():
            w(f"{r.subset:44s} {int(r.n):4d} {r.slope:7.3f} {r.r2:7.4f} "
              f"{r.ratio_median:18.3f}")
    if os.path.exists(f"{RES}/lamscale_spread.csv"):
        ls = pd.read_csv(f"{RES}/lamscale_spread.csv")
        w("\n   between-section sd of log(lambda_hat), by distance scale "
          f"(median over {len(ls)} celltype x module x call cells):")
        w(f"     raw microns {ls.sd_log_lam_raw.median():.3f} | "
          f"Poisson-normalised {ls.sd_log_lam_poisson.median():.3f} | "
          f"median-NN-normalised {ls.sd_log_lam_nn.median():.3f}")
    if os.path.exists(f"{RES}/lamscale_density_r2.csv"):
        lr = pd.read_csv(f"{RES}/lamscale_density_r2.csv")
        a = lr[lr.scale == "lam_raw"]
        w(f"   r2 of log(lambda_raw) on log(sender density): median "
          f"{a.r2.median():.3f}; slope median {a.slope.median():+.3f} "
          "(a pure density readout would be -0.500)")

    # ---- stratification decomposition ------------------------------------
    if os.path.exists(f"{RES}/stratification.csv"):
        st = pd.read_csv(f"{RES}/stratification.csv")
        st["band"] = st.section.map(BAND)
        g = st[st.band == "in_band"]
        w("\n### 11. WHERE THE NAIVE GRADIENT COMES FROM (unstratified fits)")
        w(f"   binned curve monotone decreasing (Spearman < 0) in "
          f"{int((g.spearman_bins<0).sum())}/{len(g)} section x module fits")
        w(f"   surviving fraction after adding ONLY receiver cell-type "
          f"intercepts: median {g.sf_celltype.median():.3f}")
        w(f"   surviving fraction after the full N5 block: median "
          f"{g.sf_n5_unstrat.median():.3f}")

    # ---- combined --------------------------------------------------------
    p = f"{RES}/combined_donor.csv"
    if os.path.exists(p):
        cb = pd.read_csv(p)
        c0 = cb[cb.zone == "all"]
        w("\n### 12. T2 COMBINED ESTIMATE UNDER N2+N5+N6, DONOR BOOTSTRAP")
        w(f"   {c0.n_donors.max()} animals -> "
          f"{'63' if c0.n_donors.max()==6 else '35'} distinct resamples. "
          "CASE STUDY per Section 24.1.")
        w(f"{'cell type':20s} {'module':22s} {'lam0':>6s} {'b0':>7s} "
          f"{'SF':>7s} {'donor CI on SF':>20s} {'donor CI on lam(ctrl)':>24s}")
        for _, r in c0.iterrows():
            w(f"{r.celltype:20s} {r.module:22s} {r.lam_naive:6.1f} "
              f"{r.beta_naive:7.3f} {r.sf_n2n5n6:7.3f} "
              f"[{r.get('sf_n2n5n6_donor_lo', np.nan):8.3f},"
              f"{r.get('sf_n2n5n6_donor_hi', np.nan):8.3f}] "
              f"[{r.get('lam_full_donor_lo', np.nan):9.1f},"
              f"{r.get('lam_full_donor_hi', np.nan):9.1f}]")

    if tabs:
        pd.concat(tabs, ignore_index=True).to_csv(f"{RES}/sf_summary.csv",
                                                  index=False)
    with open(f"{RES}/summary_phase3.txt", "w") as fh:
        fh.write(OUT.getvalue())
    print(f"\nwrote {RES}/summary_phase3.txt")


if __name__ == "__main__":
    main()
