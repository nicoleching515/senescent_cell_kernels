#!/usr/bin/env python3
"""Phase 9 Job B step 2 — the three matrix-based sender callers and the seven Tier B module
scores for one H1 section.  (DeepScence is the fourth caller and runs separately, in
`code/h1_deepscence.py`, because it needs the whole memory budget to itself.)

This is `code/phase2_downstream.py` D-C and D-D transplanted to H1 with no threshold change:
  * Tier A     : scanpy.tl.score_genes(A_SENDER_FINAL_strict, ctrl_size=200), then
                 strict > NNth percentile WITHIN cell type, per section, types with >= 20 cells,
                 for q in {90, 95, 99}                                (PREREG §3.7 a)
  * tierApm    : the same, scored on each A_sender_for_<module> set   (PREREG §3.7 a, D1)
  * CDKN1A+    : CDKN1A raw counts > 0, no percentile, no stratification
  * SenePy     : senepy.score_hub per cell type on the declared cross-tissue surrogate hub,
                 >= 10 hub genes on panel (the mouse arm's MIN_ON_PANEL)
  * Tier B     : score_genes(module, ctrl_size=max(200, 5*len(on)))

Declared differences from the mouse call, none of them a threshold change:
  D-a  Liver zonation (the D-B block of phase2_downstream) is not computed here; the spleen
       analogue is test A6 and lives in `code/h1_a6_compartments.py`.
  D-b  `EXCLUDE_FROM_SENDERS = ('Proliferating',)` is a no-op on H1: the frozen spleen label
       set has no 'Proliferating' label (it was dropped by the >=4 on-panel marker gate at
       build time).  Recorded so the exclusion is not silently assumed to have bitten.
  D-c  SenePy runs on the **v2** human hubs, `senepy.load_hubs(species='Human')`, which is the
       same call form and the same default the mouse arm used.  See
       `code/h1_senepy_surrogates.py` for why that is not the v1 set the frozen coverage
       table was built from.

Usage: python3 code/h1_callers.py SPLN07 [SPLN14 ...]
Writes data/processed_h1/senders_h1_<sec>.csv and modules_h1_<sec>.csv.
"""
import sys, os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 32; sc.settings.verbosity = 0
sys.path.insert(0, "/workspace/code")
import h1_common as H

PROC = H.PROC + "/"
GS = H.GS_HUMAN + "/"
EXCLUDE_FROM_STRATA = {"Low_quality", "Unknown"}      # phase2_downstream.py, verbatim
A_STRICT = H.gl("A_SENDER_FINAL_strict")
A_PM = {os.path.basename(p)[len("A_sender_for_"):-4]:
        [l.strip() for l in open(p) if l.strip()]
        for p in sorted(glob.glob(GS + "A_sender_for_*.txt"))}
BMODS = {os.path.basename(p)[2:-4]: [l.strip() for l in open(p) if l.strip()]
         for p in sorted(glob.glob(GS + "B_*.txt"))}


def run(section):
    print("\n" + "=" * 90); print("SECTION", section); print("=" * 90, flush=True)
    ctdf = pd.read_csv(PROC + "celltypes_h1_%s.csv" % section).set_index("cell_id")
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = pd.Index(bc).isin(ctdf.index)
    B = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]),
                   var=pd.DataFrame(index=names))
    B.layers["counts"] = B.X.copy()
    sc.pp.normalize_total(B); sc.pp.log1p(B)
    ctdf = ctdf.reindex(B.obs_names)
    ct = ctdf["cell_type"].astype(str)
    ctm = ctdf["cell_type_merged"].astype(str)
    n = B.n_obs
    print("cells %d ; fine types %d ; merged types %d" % (n, ct.nunique(), ctm.nunique()),
          flush=True)

    # ---- Tier A (strict-33, PRIMARY) and the per-module sensitivity sets ----
    on = [g for g in A_STRICT if g in B.var_names]
    print("Tier A strict on panel: %d/%d" % (len(on), len(A_STRICT)))
    sc.tl.score_genes(B, on, score_name="_tierA", ctrl_size=200)
    tierA = B.obs["_tierA"].to_numpy()
    pm_scores = {}
    for m, gs in A_PM.items():
        o = [g for g in gs if g in B.var_names]
        sc.tl.score_genes(B, o, score_name="_pm", ctrl_size=200)
        pm_scores[m] = B.obs["_pm"].to_numpy()

    cdk = np.asarray(B.layers["counts"][:, B.var_names.get_loc("CDKN1A")].todense()).ravel()

    # ---- SenePy, cross-tissue surrogates on the v2 hub set ----
    import senepy
    HUBS = senepy.load_hubs(species="Human").hubs
    sur = pd.read_csv(H.RESULTS + "/senepy_surrogates_v1_v2.csv").set_index("cell_type")
    sp = np.full(n, np.nan); used = {}
    for cell_type in sorted(ct.unique()):
        if cell_type in EXCLUDE_FROM_STRATA or cell_type not in sur.index:
            continue
        r = sur.loc[cell_type]
        if r.v2_usable != "yes":
            continue
        m = (ct == cell_type).to_numpy()
        if m.sum() < 50:
            print("  senepy %-30s only %d cells, skipped" % (cell_type, m.sum())); continue
        hub = HUBS[(r.v2_tissue, r.v2_cell, int(r.v2_hub))]
        onp = [(g, i) for g, i in hub if g in B.var_names]
        if len(onp) < 10:
            print("  senepy %-30s %d hub genes on panel, skipped" % (cell_type, len(onp)))
            continue
        used[cell_type] = ("%s/%s/%d" % (r.v2_tissue, r.v2_cell, int(r.v2_hub)),
                           len(hub), len(onp), int(m.sum()))
        sp[m] = senepy.score_hub(B[m].copy(), onp, verbose=False)
    print("SenePy hubs used (v2, cross-tissue surrogates -- no spleen hub exists):")
    for k, v in used.items():
        print("   %-30s %-42s hub=%5d on-panel=%4d cells=%7d" % (k, v[0], v[1], v[2], v[3]))
    noscore = sorted(set(ct.unique()) - set(used) - EXCLUDE_FROM_STRATA)
    print("   NO SENEPY SCORE for %d labels: %s" % (len(noscore), ", ".join(noscore)))

    sen = pd.DataFrame({"cell_id": B.obs_names, "cell_type": ct.values,
                        "cell_type_merged": ctm.values,
                        "cdkn1a_counts": cdk, "cdkn1a_pos": (cdk > 0).astype(int),
                        "tierA_score": np.round(tierA, 5),
                        "senepy_score": np.round(sp, 5)})
    for m, v in pm_scores.items():
        sen["tierApm_score__" + m] = np.round(v, 5)

    def pct_flags(score, labels, q, min_cells=20):
        f = np.zeros(len(score), int)
        for c in pd.unique(labels):
            if c in EXCLUDE_FROM_STRATA:
                continue
            m = (labels == c).to_numpy()
            if m.sum() < min_cells:
                continue
            f[m] = (score[m] > np.nanpercentile(score[m], q)).astype(int)
        return f

    for q in (90, 95, 99):
        sen["sender_flag_p%d" % q] = pct_flags(tierA, ct, q)
        sen["sender_flag_merged_p%d" % q] = pct_flags(tierA, ctm, q)
    for m, v in pm_scores.items():
        sen["tierApm_flag_p95__" + m] = pct_flags(v, ct, 95)
    # SenePy call: sasp_phase3.Sec.sender_mask rule -- within cell type, >= 100 finite-score
    # sender-eligible cells, strict > percentile.
    for q in (90, 95, 99):
        f = np.zeros(n, int)
        for c in pd.unique(ct):
            if c in EXCLUDE_FROM_STRATA:
                continue
            m = ((ct == c).to_numpy()) & np.isfinite(sp)
            if m.sum() < 100:
                continue
            f[m] = (sp[m] > np.percentile(sp[m], q)).astype(int)
        sen["senepy_flag_p%d" % q] = f
    sen.to_csv(PROC + "senders_h1_%s.csv" % section, index=False)

    # ---- Tier B module scores ----
    mod = pd.DataFrame({"cell_id": B.obs_names})
    for name, genes in BMODS.items():
        o = [g for g in genes if g in B.var_names]
        sc.tl.score_genes(B, o, score_name="_m", ctrl_size=max(200, len(o) * 5))
        mod[name] = np.round(B.obs["_m"].to_numpy(), 5)
    mod.to_csv(PROC + "modules_h1_%s.csv" % section, index=False)
    json.dump({"senepy_hubs_used": {k: list(v) for k, v in used.items()},
               "senepy_no_score_labels": noscore,
               "tierA_on_panel": len(on), "tierA_total": len(A_STRICT),
               "n_cells": int(n)},
              open(PROC + "caller_meta_h1_%s.json" % section, "w"), indent=1)
    print("wrote senders_h1_%s.csv, modules_h1_%s.csv" % (section, section), flush=True)
    del B


if __name__ == "__main__":
    for s in sys.argv[1:]:
        run(s)
