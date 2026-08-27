#!/usr/bin/env python3
"""Moran's I vs the A7 distance-to-sender kernel: the §29 objection-9 comparison.

Reads   results/moran/moran_fields.csv, moran_per_feature.csv.gz
        results/phase3/a7_control_probe_fits.csv, main_fits.csv
Writes  results/moran/moran_vs_a7.csv      one row per section x field
        results/moran/moran_pooled.csv     one row per field, pooled over sections
        results/moran/moran_sensitivity.csv  field x weights variant
        results/moran/moran_per_feature_summary.csv  Voyager reproduction
        results/moran/moran_verdict.txt

The question is not "is Moran's I significant" -- with n ~ 2.4e5 cells per
section SE(I) ~ 2e-3 and everything is significant.  The question is whether a
near-zero Moran's I on the controls coexists with a non-zero A7 kernel
amplitude on the SAME features, because that is the whole of the project's
"different question" defence.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import t as tdist, spearmanr

MOR = "/workspace/results/moran/"
RES = "/workspace/results/phase3/"
CALLS = ["tierA_p95", "cdkn1a_pos"]
DESIGNS = ["base", "n6", "n5", "n6n5", "n2"]
CONTROLS = ["neg_control_probe", "neg_control_codeword", "genomic_control",
            "all_controls", "neg_probe_rate"]


def add(df):
    for k in DESIGNS:
        df["e_" + k] = df["beta_" + k] / df.sd_y
    return df


def clustered(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return np.nan, np.nan, np.nan, np.nan
    mu, se = v.mean(), v.std(ddof=1) / np.sqrt(v.size)
    tc = tdist.ppf(.975, v.size - 1)
    return mu, mu - tc * se, mu + tc * se, 2 * tdist.sf(abs(mu / se), v.size - 1)


def main():
    M = pd.read_csv(MOR + "moran_fields.csv")
    A = add(pd.read_csv(RES + "a7_control_probe_fits.csv"))
    B = pd.read_csv(RES + "main_fits.csv")
    B = add(B[(B.stratum == "all") & B.call.isin(CALLS)].copy())
    B = B.rename(columns={"module": "response"})
    K = pd.concat([A[["section", "response"] + ["e_" + k for k in DESIGNS]],
                   B[["section", "response"] + ["e_" + k for k in DESIGNS]]])
    # A7 amplitude for one section = mean over sender calls and receiver types
    ksec = K.groupby(["section", "response"]).agg(
        n_fits=("e_base", "size"),
        **{"a7_" + k: ("e_" + k, "mean") for k in DESIGNS}).reset_index()
    ksec = ksec.rename(columns={"response": "field"})

    # ---- per section x field -------------------------------------------
    M2 = M.assign(key=M.weights + "_" + M.centring)
    wide = M2.pivot_table(index=["section", "arm", "band", "n_cells", "field",
                                 "kind", "frac_nonzero"],
                          columns="key", values="moran_I").reset_index()
    wide.columns.name = None
    # per-section A7 CI-exclusion rate on the naive design (the "is there a
    # gradient here at all" read, section by section)
    A["sig_base"] = (A.beta_base_lo * A.beta_base_hi) > 0
    sig = (A.groupby(["section", "response"])
           .sig_base.mean().reset_index()
           .rename(columns={"response": "field", "sig_base": "a7_frac_CI_excl_0"}))
    J = wide.merge(ksec, on=["section", "field"], how="left") \
            .merge(sig, on=["section", "field"], how="left")
    J.to_csv(MOR + "moran_vs_a7.csv", index=False)

    # ---- pooled over sections ------------------------------------------
    prim = M[(M.weights == "knn6")]
    rows = []
    for (fld, kind), g in prim.groupby(["field", "kind"]):
        r = dict(field=fld, kind=kind, n_sections=g.section.nunique())
        for cs in ("raw", "ctcentred"):
            h = g[g.centring == cs]
            mu, lo, hi, p = clustered(h.moran_I)
            r["I_%s_mean" % cs] = round(mu, 5)
            r["I_%s_lo" % cs] = round(lo, 5)
            r["I_%s_hi" % cs] = round(hi, 5)
            r["I_%s_p" % cs] = "%.3g" % p
            r["I_%s_min" % cs] = round(h.moran_I.min(), 5)
            r["I_%s_max" % cs] = round(h.moran_I.max(), 5)
            if cs == "raw":
                r["max_p_rand"] = "%.3g" % h.p_rand.max()
                r["max_p_sim"] = ("%.3g" % h.p_sim.max()
                                  if h.p_sim.notna().any() else "NA")
        kk = ksec[ksec.field == fld]
        for k in DESIGNS:
            mu, lo, hi, p = clustered(kk["a7_" + k])
            r["a7_%s_mean" % k] = round(mu, 5) if np.isfinite(mu) else np.nan
            r["a7_%s_lo" % k] = round(lo, 5) if np.isfinite(lo) else np.nan
            r["a7_%s_hi" % k] = round(hi, 5) if np.isfinite(hi) else np.nan
            r["a7_%s_p" % k] = ("%.3g" % p) if np.isfinite(p) else "NA"
        rows.append(r)
    Pl = pd.DataFrame(rows).sort_values(["kind", "field"])
    Pl.to_csv(MOR + "moran_pooled.csv", index=False)

    # ---- weights sensitivity -------------------------------------------
    sens = (M[M.centring == "raw"]
            .groupby(["field", "kind", "weights"])
            .agg(mean_neighbours=("mean_neighbours", "mean"),
                 islands=("islands", "max"),
                 I_mean=("moran_I", "mean"), I_min=("moran_I", "min"),
                 I_max=("moran_I", "max"))
            .reset_index().sort_values(["kind", "field", "mean_neighbours"]))
    sens.to_csv(MOR + "moran_sensitivity.csv", index=False)

    # ---- Voyager reproduction, per feature ------------------------------
    F = pd.read_csv(MOR + "moran_per_feature.csv.gz")
    fs = (F.groupby(["feature_type"])
          .agg(n_feat=("feature", "nunique"), n_rows=("moran_I", "size"),
               I_median=("moran_I", "median"), I_mean=("moran_I", "mean"),
               I_q05=("moran_I", lambda s: s.quantile(.05)),
               I_q95=("moran_I", lambda s: s.quantile(.95)),
               I_max=("moran_I", "max"),
               frac_I_gt_005=("moran_I", lambda s: float((s > 0.05).mean())),
               counts_median=("total_counts", "median"),
               nonzero_median=("frac_cells_nonzero", "median"))
          .reset_index().sort_values("I_median"))
    fs.to_csv(MOR + "moran_per_feature_summary.csv", index=False)

    # sparsity-matched genes: for each control feature, genes with total counts
    # in the same decile band -- does a gene as sparse as a control probe also
    # look "flat"?
    G = F[F.feature_type == "Gene Expression"]
    L = []
    L.append("SPARSITY-MATCHED GENE COMPARISON (primary weights knn6, all "
             "sections pooled)")
    L.append("  a control feature's near-zero I is only evidence of 'no spatial")
    L.append("  trend' if a GENE of the same sparsity would have shown one.")
    L.append("%-28s %8s %10s %10s   %8s %10s %10s" %
             ("feature class", "n", "med cts", "med I", "n genes",
              "med cts", "med I"))
    for ftp in ["Negative Control Probe", "Negative Control Codeword",
                "Genomic Control", "Unassigned Codeword"]:
        h = F[F.feature_type == ftp]
        lo, hi = h.total_counts.quantile([.1, .9])
        g = G[(G.total_counts >= lo) & (G.total_counts <= hi)]
        L.append("%-28s %8d %10.1f %+10.5f   %8d %10.1f %+10.5f" %
                 (ftp, len(h), h.total_counts.median(), h.moran_I.median(),
                  len(g), g.total_counts.median(), g.moran_I.median()))
    L.append("")

    # ---- verdict --------------------------------------------------------
    pool = Pl.set_index("field")
    L.insert(0, "")
    hdr = []
    hdr.append("MORAN'S I vs THE A7 DISTANCE KERNEL -- M1 (mouse), %d sections, "
               "n=%s cells" % (M.section.nunique(),
                               "/".join(str(int(x)) for x in
                                        sorted(M.n_cells.unique())[:2]) + "..."))
    hdr.append("primary weights: knn6 row-standardised.  Moran's I pooled as a "
               "section-clustered mean [95%% CI].")
    hdr.append("A7 amplitude = beta/sd_y, same responses, same cells, "
               "section-clustered mean.")
    hdr.append("")
    hdr.append("%-22s %-9s %24s %24s %24s" %
               ("field", "kind", "Moran I (raw)", "Moran I (ct-centred)",
                "A7 naive beta/sd"))
    for fld in CONTROLS + list(pool[pool.kind == "module"].index) + \
            list(pool[pool.kind == "technical"].index):
        if fld not in pool.index:
            continue
        r = pool.loc[fld]
        hdr.append("%-22s %-9s %+8.4f [%+.4f,%+.4f] %+8.4f [%+.4f,%+.4f] "
                   "%+8.4f [%+.4f,%+.4f]" %
                   (fld, r.kind, r.I_raw_mean, r.I_raw_lo, r.I_raw_hi,
                    r.I_ctcentred_mean, r.I_ctcentred_lo, r.I_ctcentred_hi,
                    r.a7_base_mean, r.a7_base_lo, r.a7_base_hi))
    hdr.append("")

    # agreement between the two statistics across fields
    sub = pool[pool.kind.isin(["control", "module"])]
    rho1, p1 = spearmanr(sub.I_raw_mean.abs(), sub.a7_base_mean.abs())
    rho2, p2 = spearmanr(sub.I_ctcentred_mean.abs(), sub.a7_base_mean.abs())
    hdr.append("Rank agreement across the %d control+module fields, |Moran I| vs "
               "|A7 naive amplitude|:" % len(sub))
    hdr.append("  raw        Spearman rho = %+.3f (p = %.3g)" % (rho1, p1))
    hdr.append("  ct-centred Spearman rho = %+.3f (p = %.3g)" % (rho2, p2))
    hdr.append("")
    # per-section agreement within the control family only
    cj = J[J.field.isin(CONTROLS)].dropna(subset=["a7_base"])
    if len(cj) > 4:
        r3, p3 = spearmanr(cj["knn6_raw"].abs(), cj["a7_base"].abs())
        hdr.append("Within the 5 control responses x %d sections (%d pairs): "
                   "Spearman rho = %+.3f (p = %.3g)"
                   % (cj.section.nunique(), len(cj), r3, p3))
    hdr.append("")

    # ---- discordance: near-zero Moran I with a non-null A7 fit ----------
    gene_med = F[F.feature_type == "Gene Expression"].moran_I.median()
    hdr.append("DISCORDANCE TEST -- can a field be flat by Moran's I and carry a")
    hdr.append("distance-to-sender gradient?  'near zero' = |I| below the median")
    hdr.append("Gene Expression feature (%.5f) under knn6." % gene_med)
    hdr.append("%-22s %-30s %8s %10s %10s" %
               ("field", "section", "Moran I", "A7 naive", "CI excl 0"))
    nd = 0
    for _, r in J[J.field.isin(CONTROLS)].sort_values(
            ["field", "section"]).iterrows():
        near = abs(r["knn6_raw"]) < gene_med
        strong = (r["a7_frac_CI_excl_0"] or 0) >= 0.25
        if near and strong:
            nd += 1
            hdr.append("%-22s %-30s %+8.5f %+10.4f %10.2f" %
                       (r.field, r.section, r["knn6_raw"], r["a7_base"],
                        r["a7_frac_CI_excl_0"]))
    hdr.append("  %d discordant section x response cells out of %d."
               % (nd, len(J[J.field.isin(CONTROLS)])))
    hdr.append("")

    txt = "\n".join(hdr + L)
    open(MOR + "moran_verdict.txt", "w").write(txt + "\n")
    pd.set_option("display.width", 250)
    print(txt)
    print()
    print(sens[sens.field.isin(CONTROLS + ["emt_ecm", "transcript_counts"])]
          .to_string(index=False))
    print()
    print(fs.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
