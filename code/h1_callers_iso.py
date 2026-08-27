#!/usr/bin/env python3
"""Phase 10 ISO -- `code/h1_callers.py` re-run with the H1 panel cut to the
ortholog-intersected 2,425 human symbols (PREREG_PHASE8.md §9 item 4 / test A8).

WHY THE SUBSET HAPPENS FIRST.  `scanpy.tl.score_genes` draws its control genes from the
expression-binned background of whatever is in `adata.var`, and `normalize_total` divides by
the per-cell total over whatever is in `adata.var`.  Restricting the universe is therefore
the entire point of this script: the AnnData var is cut to the intersection BEFORE
`normalize_total`/`log1p`/`score_genes`, exactly as Phase 9's
`h1_a8_crossarm.score_restricted` does it.  Subsetting after scoring would change nothing.

IDENTICAL TO `code/h1_callers.py`, unchanged:
  * scoring convention: `sc.pp.normalize_total` -> `sc.pp.log1p` -> `sc.tl.score_genes`
  * Tier A `A_SENDER_FINAL_strict`, `ctrl_size=200`
  * Tier B module m, `ctrl_size=max(200, 5*len(on_panel_genes))`
  * sender flags: within-cell-type strict `>` percentile at q in {90,95,99}, excluding
    `Low_quality`/`Unknown`, skipping types with < 20 cells (`pct_flags`, copied verbatim
    from `code/h1_cache_extend.py`, which copied it verbatim from `code/h1_callers.py`)
  * both label families -- fine `cell_type` and merged `cell_type_merged`
  * CDKN1A raw counts from the counts layer.  CDKN1A IS on the intersected panel (checked,
    not assumed -- `iso_panel.human_symbols()` contains it; the script asserts this and
    reports it), so `cdkn1a_counts` is the same raw integer vector as the full-panel run.

DECLARED OMISSIONS (both out of scope for Phase 10 ISO, neither is a threshold change):
  * SenePy is NOT run.  Its hub gene sets are not the panel and the intersected panel would
    drop most hubs below the frozen >= 10-on-panel floor; the SenePy call is not among the
    calls the ISO fits use.  `senepy_score` is left absent from these CSVs and the ISO cache
    keeps the FULL-PANEL `senepy_score` key from the original .npz, unused.
  * the seven per-module Tier A sets `A_sender_for_<module>` are NOT scored, so
    `tierApm_pNN` is not an ISO call.  The ISO cache keeps the full-panel `tierApm__*` and
    `flag_pm_*` keys, unused.  Do not run `tierApm_*` against the ISO cache.

Usage: python3 code/h1_callers_iso.py SPLN07 [SPLN14 ...]
Writes data/processed_h1/senders_h1_iso_<sec>.csv, modules_h1_iso_<sec>.csv,
       caller_meta_h1_iso_<sec>.json  (all NEW paths).
"""
import sys, os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 4; sc.settings.verbosity = 0
sys.path.insert(0, "/workspace/code")
import h1_common as H
import iso_panel

PROC = H.PROC + "/"
GS = H.GS_HUMAN + "/"
EXCLUDE_FROM_STRATA = {"Low_quality", "Unknown"}      # h1_callers.py, verbatim
A_STRICT = H.gl("A_SENDER_FINAL_strict")
BMODS = {os.path.basename(p)[2:-4]: [l.strip() for l in open(p) if l.strip()]
         for p in sorted(glob.glob(GS + "B_*.txt"))}
ISO = iso_panel.human_symbols()


def pct_flags(score, labels, q, min_cells=20):
    """`code/h1_cache_extend.py::pct_flags`, verbatim (itself verbatim from h1_callers)."""
    f = np.zeros(len(score), int)
    for c in pd.unique(labels):
        if c in EXCLUDE_FROM_STRATA:
            continue
        m = (labels == c).to_numpy()
        if m.sum() < min_cells:
            continue
        f[m] = (score[m] > np.nanpercentile(score[m], q)).astype(int)
    return f


def run(section):
    print("\n" + "=" * 90); print("SECTION", section, "(ORTHOLOG-INTERSECTED PANEL)")
    print("=" * 90, flush=True)
    ctdf = pd.read_csv(PROC + "celltypes_h1_%s.csv" % section).set_index("cell_id")
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = pd.Index(bc).isin(ctdf.index)
    B = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]),
                   var=pd.DataFrame(index=names))
    n_full_panel = B.n_vars
    # ---- THE RESTRICTION.  Before normalisation, before scoring. ----
    B = B[:, [g in ISO for g in B.var_names]].copy()
    n_iso_panel = B.n_vars
    assert n_iso_panel == len(ISO & set(names)), section + ": iso panel subset lost genes"
    cdkn1a_on = "CDKN1A" in set(B.var_names)
    print("panel %d -> %d  (intersected; CDKN1A on intersected panel: %s)"
          % (n_full_panel, n_iso_panel, cdkn1a_on), flush=True)
    assert cdkn1a_on, "CDKN1A is not on the intersected panel -- cdkn1a_counts cannot be carried"
    B.layers["counts"] = B.X.copy()
    sc.pp.normalize_total(B); sc.pp.log1p(B)
    ctdf = ctdf.reindex(B.obs_names)
    ct = ctdf["cell_type"].astype(str)
    ctm = ctdf["cell_type_merged"].astype(str)
    n = B.n_obs
    print("cells %d ; fine types %d ; merged types %d"
          % (n, ct.nunique(), ctm.nunique()), flush=True)

    surv = {}
    # ---- Tier A ----
    on = [g for g in A_STRICT if g in B.var_names]
    on_full = [g for g in A_STRICT if g in names]
    surv["A_SENDER_FINAL_strict"] = dict(total=len(A_STRICT), on_full_panel=len(on_full),
                                         on_iso_panel=len(on))
    print("Tier A strict: %d total, %d on full panel, %d on intersected panel"
          % (len(A_STRICT), len(on_full), len(on)))
    sc.tl.score_genes(B, on, score_name="_tierA", ctrl_size=200)
    tierA = B.obs["_tierA"].to_numpy()

    cdk = np.asarray(B.layers["counts"][:, B.var_names.get_loc("CDKN1A")].todense()).ravel()

    sen = pd.DataFrame({"cell_id": B.obs_names, "cell_type": ct.values,
                        "cell_type_merged": ctm.values,
                        "cdkn1a_counts": cdk, "cdkn1a_pos": (cdk > 0).astype(int),
                        "tierA_score": np.round(tierA, 5)})
    for q in (90, 95, 99):
        sen["sender_flag_p%d" % q] = pct_flags(tierA, ct, q)
        sen["sender_flag_merged_p%d" % q] = pct_flags(tierA, ctm, q)
    sen.to_csv(PROC + "senders_h1_iso_%s.csv" % section, index=False)

    # ---- Tier B module scores ----
    mod = pd.DataFrame({"cell_id": B.obs_names})
    for name, genes in BMODS.items():
        o = [g for g in genes if g in B.var_names]
        o_full = [g for g in genes if g in names]
        surv["B_" + name] = dict(total=len(genes), on_full_panel=len(o_full),
                                 on_iso_panel=len(o))
        print("  B_%-24s %3d total, %3d full, %3d iso, ctrl_size=%d"
              % (name, len(genes), len(o_full), len(o), max(200, len(o) * 5)))
        sc.tl.score_genes(B, o, score_name="_m", ctrl_size=max(200, len(o) * 5))
        mod[name] = np.round(B.obs["_m"].to_numpy(), 5)
    mod.to_csv(PROC + "modules_h1_iso_%s.csv" % section, index=False)

    json.dump({"section": section, "panel_full": int(n_full_panel),
               "panel_intersected": int(n_iso_panel), "n_cells": int(n),
               "cdkn1a_on_intersected_panel": bool(cdkn1a_on),
               "gene_survival": surv,
               "senepy": "NOT RUN (declared)",
               "tierApm_A_sender_for_sets": "NOT SCORED (declared)",
               "sender_p95_fine": int(sen.sender_flag_p95.sum()),
               "sender_p95_merged": int(sen.sender_flag_merged_p95.sum())},
              open(PROC + "caller_meta_h1_iso_%s.json" % section, "w"), indent=1)
    print("wrote senders_h1_iso_%s.csv, modules_h1_iso_%s.csv" % (section, section),
          flush=True)
    del B


if __name__ == "__main__":
    for s in sys.argv[1:]:
        run(s)
