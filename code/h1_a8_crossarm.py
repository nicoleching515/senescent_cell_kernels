#!/usr/bin/env python3
"""Phase 9 test A8 — cross-arm comparability on the ortholog-intersected panel.

§13 A8 / PREREG §9.4: every cross-arm number is reported twice — on the ortholog-intersected
panel and on each arm's full panel — with the ortholog map pinned by version and archived.

MAP: genesets/mouse_human_orthologs_MGI.csv (pinned, MGI 1:1).  The intersected panel is the
set of HUMAN symbols reachable from the mouse panel and present on the human panel;
PREREG_PHASE8_genesets.md §  fixes it at **2,425 human symbols** (4,845 of 5,097 mouse panel
genes have a map row; 2,435 of those map onto the human panel; they land on 2,425 distinct
human symbols).  This script re-derives all four counts and fails loudly if any moves.

Part 1  panel and gene-set arithmetic, both arms, full vs intersected.
Part 2  the H1 sender call and the seven Tier B module scores RECOMPUTED with the panel
        restricted to the 2,425 intersection, and the resulting A3 prevalence compared to
        the full-panel call.  A caller whose prevalence or sender set moves when half the
        panel is removed is not comparable across arms, and that has to be measured.

Usage: python3 code/h1_a8_crossarm.py [--score SPLN07 ...]
Writes results/phase9_h1/a8_panel_arithmetic.csv, a8_ortho_sender_shift.csv
       and data/processed_h1/senders_h1_ortho_<sec>.csv
"""
import sys, os, csv, gzip, glob, argparse, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, "/workspace/code")
import h1_common as H

W = "/workspace"
GS_M, GS_H = W + "/genesets", W + "/genesets/human"
MODS = H.MODULES


def panels():
    mp = {r["gene_name"] for r in csv.DictReader(
        open(W + "/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv"))}
    p100 = [r["Gene"] for r in csv.DictReader(
        gzip.open(W + "/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz", "rt"))]
    geno = {g for g in p100 if ("_WT" in g or "_ALT" in g or "_del_" in g or "_splice_" in g)}
    M = (mp | set(p100)) - geno
    _, hn, _, _ = H.load_matrix("SPLN21", "gene")
    return M, set(hn)


def build():
    MPANEL, HPANEL = panels()
    ortho = {r["mouse_symbol"]: r["human_symbol"] for r in
             csv.DictReader(open(W + "/genesets/mouse_human_orthologs_MGI.csv"))}
    mapped = {g: ortho[g] for g in MPANEL if g in ortho}
    onto = {g: h for g, h in mapped.items() if h in HPANEL}
    INTER = set(onto.values())
    print("mouse panel                       %d   (expected 5097)" % len(MPANEL))
    print("human panel                       %d   (expected 5093)" % len(HPANEL))
    print("mouse genes with an MGI map row   %d   (expected 4845)" % len(mapped))
    print("  ... whose ortholog is on-panel  %d   (expected 2435)" % len(onto))
    print("distinct human symbols            %d   (expected 2425)" % len(INTER))
    assert (len(MPANEL), len(HPANEL), len(mapped), len(onto), len(INTER)) == \
           (5097, 5093, 4845, 2435, 2425), "the pinned panel arithmetic has moved"
    return MPANEL, HPANEL, INTER, onto


def rd(p):
    return {l.strip() for l in open(p) if l.strip()}


def arithmetic(MPANEL, HPANEL, INTER, onto):
    rows = []
    sets = [("A_SENDER_FINAL_strict", "A_SENDER_FINAL_strict")] + \
           [("B_" + m, "B_" + m) for m in MODS] + \
           [("A_sender_for_" + m, "A_sender_for_" + m) for m in MODS] + \
           [("E_negative_control_probes", "E_negative_control_probes")]
    m_inter = set(onto)                       # mouse symbols whose ortholog is on both panels
    for hn, mn in sets:
        hp = os.path.join(GS_H, hn + ".txt"); mp = os.path.join(GS_M, mn + ".txt")
        h = rd(hp) if os.path.exists(hp) else set()
        m = rd(mp) if os.path.exists(mp) else set()
        rows.append(dict(gene_set=hn,
                         human_full=len(h & HPANEL), human_intersected=len(h & INTER),
                         human_intersected_pct=round(100 * len(h & INTER) /
                                                     max(len(h & HPANEL), 1), 1),
                         mouse_full=len(m & MPANEL), mouse_intersected=len(m & m_inter),
                         mouse_intersected_pct=round(100 * len(m & m_inter) /
                                                     max(len(m & MPANEL), 1), 1)))
    return pd.DataFrame(rows)


def score_restricted(section, INTER):
    """Recompute Tier A and the seven modules with the panel cut to the intersection."""
    import scanpy as sc, anndata as ad
    sc.settings.n_jobs = 32; sc.settings.verbosity = 0
    ct = pd.read_csv(H.PROC + "/celltypes_h1_%s.csv" % section).set_index("cell_id")
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = pd.Index(bc).isin(ct.index)
    B = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]), var=pd.DataFrame(index=names))
    B = B[:, [g in INTER for g in B.var_names]].copy()
    B.layers["counts"] = B.X.copy()
    sc.pp.normalize_total(B); sc.pp.log1p(B)
    ct = ct.reindex(B.obs_names)
    lab = ct["cell_type"].astype(str)
    on = [g for g in H.gl("A_SENDER_FINAL_strict") if g in B.var_names]
    sc.tl.score_genes(B, on, score_name="_t", ctrl_size=200)
    tierA = B.obs["_t"].to_numpy()
    out = pd.DataFrame({"cell_id": B.obs_names, "cell_type": lab.values,
                        "tierA_score_ortho": np.round(tierA, 5)})
    for q in (90, 95, 99):
        f = np.zeros(B.n_obs, int)
        for c in pd.unique(lab):
            if c in ("Low_quality", "Unknown"):
                continue
            m = (lab == c).to_numpy()
            if m.sum() < 20:
                continue
            f[m] = (tierA[m] > np.percentile(tierA[m], q)).astype(int)
        out["sender_flag_ortho_p%d" % q] = f
    BM = {os.path.basename(p)[2:-4]: [l.strip() for l in open(p) if l.strip()]
          for p in sorted(glob.glob(GS_H + "/B_*.txt"))}
    for nm, gs in BM.items():
        o = [g for g in gs if g in B.var_names]
        sc.tl.score_genes(B, o, score_name="_m", ctrl_size=max(200, len(o) * 5))
        out["mod_ortho__" + nm] = np.round(B.obs["_m"].to_numpy(), 5)
    out.to_csv(H.PROC + "/senders_h1_ortho_%s.csv" % section, index=False)
    full = pd.read_csv(H.PROC + "/senders_h1_%s.csv" % section).set_index("cell_id") \
             .reindex(out.cell_id)
    rows = []
    for q in (90, 95, 99):
        a = full["sender_flag_p%d" % q].to_numpy(bool)
        b = out["sender_flag_ortho_p%d" % q].to_numpy(bool)
        rows.append(dict(section=section, call="tierA_p%d" % q,
                         n_cells=len(out), n_full=int(a.sum()), n_ortho=int(b.sum()),
                         prevalence_full=round(100 * a.mean(), 3),
                         prevalence_ortho=round(100 * b.mean(), 3),
                         jaccard=round(float((a & b).sum() / max((a | b).sum(), 1)), 4),
                         spearman_score=round(float(pd.Series(full.tierA_score.to_numpy())
                                                    .corr(pd.Series(tierA), method="spearman")), 4),
                         n_tierA_on_full=len([g for g in H.gl("A_SENDER_FINAL_strict")
                                              if g in names]),
                         n_tierA_on_ortho=len(on)))
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--score", nargs="*", default=None)
    a = ap.parse_args()
    MPANEL, HPANEL, INTER, onto = build()
    d = arithmetic(MPANEL, HPANEL, INTER, onto)
    d.to_csv(H.RESULTS + "/a8_panel_arithmetic.csv", index=False)
    pd.set_option("display.width", 220)
    print(); print(d.to_string(index=False))
    secs = a.score if a.score is not None else \
        [s for s in H.ALL_SECTIONS if os.path.exists(H.PROC + "/senders_h1_%s.csv" % s)]
    rows = []
    for s in secs:
        rows += score_restricted(s, INTER)
    if rows:
        pd.DataFrame(rows).to_csv(H.RESULTS + "/a8_ortho_sender_shift.csv", index=False)
    print("wrote", H.RESULTS + "/a8_*.csv")
