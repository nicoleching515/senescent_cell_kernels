"""Markdown tables for CS_PHASE4.md, generated from results/phase4/*.csv so the
report cannot drift from the numbers."""
from __future__ import annotations
import numpy as np, pandas as pd

P4 = "/workspace/results/phase4"
M_ORDER = ["COMMOT", "CellChat v2*", "SpaTalk*", "NCEM linear*"]
NULLS = ["N3_lig", "N4_lig", "N3_type", "N0_perm"]
NLAB = {"N3_lig": "N3 torus (ligand⁺)", "N4_lig": "N4 rotation (ligand⁺)",
        "N3_type": "N3t per-type torus", "N0_perm": "N0 full permutation"}
PAIRS = ["Ccl2->Ccr2", "Tnf->Tnfrsf1a/1b", "Tgfb1->Tgfbr1/2", "Il1a->Il1r1"]


def _md(df, floatfmt="%.3f"):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            df[c] = df[c].map(lambda v: "—" if not np.isfinite(v) else floatfmt % v)
    head = "| " + " | ".join([df.index.name or ""] + list(df.columns)) + " |"
    sep = "|" + "|".join(["---"] * (len(df.columns) + 1)) + "|"
    body = ["| " + " | ".join([str(i)] + [str(v) for v in row]) + " |"
            for i, row in zip(df.index, df.values)]
    return "\n".join([head, sep] + body)


def build():
    H = pd.read_csv(f"{P4}/headline.csv")
    V = pd.read_csv(f"{P4}/verdict.csv")
    CR = pd.read_csv(f"{P4}/score_rank_correlation.csv")
    A = H[H.pair == "ALL"]
    out = {}

    t = A.pivot(index="method", columns="null", values="sig_survival").reindex(M_ORDER)[NULLS]
    t.columns = [NLAB[c] for c in t.columns]; t.index.name = "method"
    out["sig_survival"] = _md(t)

    t = A.pivot(index="method", columns="null", values="score_sf_median").reindex(M_ORDER)[NULLS]
    t.columns = [NLAB[c] for c in t.columns]; t.index.name = "method"
    out["score_sf"] = _md(t)

    t = CR.pivot(index="method", columns="null", values="spearman").reindex(M_ORDER)[NULLS]
    t.columns = [NLAB[c] for c in t.columns]; t.index.name = "method"
    out["rank_corr"] = _md(t)

    t = A.pivot(index="method", columns="null", values="sig_survival_all").reindex(M_ORDER)[NULLS]
    t.insert(0, "real", A.groupby("method").real_sig_rate.first().reindex(M_ORDER))
    t.columns = ["real coordinates"] + [NLAB[c] for c in NULLS]; t.index.name = "method"
    out["sig_rate"] = _md(t)

    v = V.copy()
    v["null"] = v["null"].map(NLAB)
    v = v.set_index("method").loc[[m for m in M_ORDER if m in set(v.method)]]
    v = v.reset_index()[["method", "null", "score_sf", "sig_survival",
                         "null_sig_rate", "n_real_sig", "n_rep", "verdict"]]
    v.index = v.pop("method"); v.index.name = "method"
    out["verdict"] = _md(v)

    P = H[H.pair != "ALL"]
    for key, col in (("pair_sig_survival", "sig_survival"),
                     ("pair_score_sf", "score_sf_median")):
        rows = []
        for m in M_ORDER:
            for pr in PAIRS:
                g = P[(P.method == m) & (P.pair == pr)]
                r = {"pair": pr}
                for n in NULLS:
                    gg = g[g.null == n]
                    r[NLAB[n]] = float(gg[col].iloc[0]) if len(gg) else np.nan
                r["n real-sig"] = int(g.n_real_sig.max()) if len(g) else 0
                r["method"] = m
                rows.append(r)
        d = pd.DataFrame(rows)
        d.index = d.pop("method") + " · " + d.pop("pair")
        d.index.name = "method · LR pair"
        out[key] = _md(d)

    # positive controls
    try:
        PC = pd.read_csv(f"{P4}/positive_controls.csv")
        blocks = []
        for ctl, g in PC.groupby("control", sort=False):
            d = g.pivot_table(index="method", columns="truth",
                              values=["score_AB", "p_AB"], aggfunc="first")
            cols = list(dict.fromkeys(g.truth))
            rows2 = []
            for m in M_ORDER:
                if m not in set(g.method):
                    continue
                r = {"method": m}
                for c in cols:
                    h = g[(g.method == m) & (g.truth == c)]
                    if not len(h):
                        continue
                    sc, pv = h.score_AB.iloc[0], h.p_AB.iloc[0]
                    r[c] = (("—" if not np.isfinite(pv) else
                             ("score %.6g, p %.3g" % (sc, pv) if np.isfinite(sc)
                              else "p %.3g" % pv)))
                rows2.append(r)
            d2 = pd.DataFrame(rows2).set_index("method")
            d2.index.name = "method"
            blocks.append("**%s**\n\n%s" % (ctl, _md(d2)))
        out["positive_controls"] = "\n\n".join(blocks)
    except Exception as ex:
        out["positive_controls"] = f"_(pending: {ex})_"

    # the illustrative example: COMMOT's top Ccl2->Ccr2 calls, real vs N0
    try:
        I = pd.read_csv(f"{P4}/interactions.csv.gz")
        d = I[(I.method == "COMMOT") & (I.pair == "Ccl2->Ccr2") &
              (I.null == "N0_perm")]
        g = d.groupby(["sender", "receiver"]).agg(
            real=("real_sig", "mean"), null=("null_sig_frac", "mean"),
            n=("tile", "nunique"))
        g = g[g.n >= 4].sort_values(["real", "null"], ascending=False).head(8)
        e = pd.DataFrame({
            "sender → receiver": [f"{a} → {b}" for a, b in g.index],
            "tiles": g.n.values,
            "significant on real coordinates": g.real.values,
            "significant on N0-permuted coordinates": g.null.values})
        e.index = e.pop("sender → receiver"); e.index.name = "sender → receiver"
        out["commot_example"] = _md(e, "%.2f")
    except Exception as ex:
        out["commot_example"] = f"_(pending: {ex})_"

    try:
        MEC = pd.read_csv(f"{P4}/commot_mechanism.csv")
        MEC["mass_ratio"] = MEC.total_mass_perm / MEC.total_mass_real
        gb = MEC.groupby("pair")
        d = pd.DataFrame(index=PAIRS)
        d["communicating cell pairs, real"] = gb.n_edges_real.median().reindex(PAIRS).map("{:,.0f}".format)
        d["Jaccard of those pairs, real vs N0"] = gb.edge_jaccard.median().reindex(PAIRS).map("{:.4f}".format)
        d["transported mass, N0 ÷ real"] = gb.mass_ratio.median().reindex(PAIRS).map("{:.6f}".format)
        d["cluster-level Spearman, real vs N0"] = gb.cluster_spearman.median().reindex(PAIRS).map("{:.3f}".format)
        d["significant: real / N0 / shared"] = [
            "%d / %d / %d" % (gb.nsig_real.sum()[p], gb.nsig_perm.sum()[p],
                              gb.nsig_shared.sum()[p]) for p in PAIRS]
        d.index.name = "LR pair"
        out["commot_mechanism"] = _md(d)
    except Exception as e:                                # not run yet
        out["commot_mechanism"] = f"_(pending: {e})_"

    try:
        S = pd.read_csv(f"{P4}/ncem_radius_sweep.csv")
        best = S.loc[S.groupby(["tile", "cond"]).r2.idxmax()]
        order = ["real", "N3_lig", "N4_lig", "N3_type", "N0_perm"]
        d = pd.DataFrame({
            "median selected radius (µm)":
                best.groupby("cond").radius.median().reindex(order).map("{:.0f}".format),
            "range over %d tiles (µm)" % best.tile.nunique():
                best.groupby("cond").radius.agg(
                    lambda v: "%g–%g" % (v.min(), v.max())).reindex(order),
            "median R² at the selected radius":
                best.groupby("cond").r2.median().reindex(order).map("{:.4f}".format),
            "median R² over the whole sweep":
                S.groupby("cond").r2.median().reindex(order).map("{:.4f}".format),
        })
        d.index = ["real coordinates"] + [NLAB[c] for c in NULLS]
        d.index.name = "coordinates"
        out["ncem_lengthscale"] = _md(d, "%.4f")
    except Exception as e:
        out["ncem_lengthscale"] = f"_(pending: {e})_"

    return out


if __name__ == "__main__":
    for k, v in build().items():
        print(f"\n### {k}\n\n{v}")
