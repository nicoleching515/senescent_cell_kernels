#!/usr/bin/env python3
"""Phase 10 ISO -- the M1 mouse arm rescored on the ortholog-intersected panel.

This is the mouse counterpart of `code/h1_callers_iso.py` and the D-C/D-D blocks of
`code/phase2_downstream.py` with the panel cut to the 2,435 MOUSE symbols whose MGI
ortholog is on BOTH panels (`iso_panel.mouse_symbols()`; the same `set(onto)` convention
`h1_a8_crossarm.arithmetic` uses for its `mouse_intersected` column).

The var is subset BEFORE `normalize_total` -- `score_genes` draws its control genes from
the panel, so restricting the universe is the whole point.

IDENTICAL TO `code/phase2_downstream.py`, unchanged:
  * `sc.pp.normalize_total` -> `sc.pp.log1p` -> `sc.tl.score_genes`
  * Tier A `genesets/A_SENDER_FINAL_strict.txt`, `ctrl_size=200`
  * Tier B module m, `ctrl_size=max(200, 5*len(on_panel_genes))`
  * sender flags: within-cell-type strict `>` percentile at q in {90,95,99},
    `EXCLUDE_FROM_STRATA = {'Low_quality','Unknown'}`, strata with < 20 cells skipped
  * `Cdkn1a` raw counts from the counts layer.  Cdkn1a IS on the intersected panel
    (checked, not assumed; asserted below).

DECLARED CHOICES
  * Expression comes from `sasp_real.load_expression(sample)` (the loader
    `code/run_phase3_n8.py` and `sasp_phase3.prep` use), NOT `sasp_io.load`.  Cells are
    then restricted, IN CACHE ORDER, to the cell set of the existing
    `data/processed/cache3/<SAMPLE>.npz`, which is exactly the QC-passing labelled set
    `sasp_io.load` + the celltypes CSV produced.  The two loaders agree on that cell set.
  * CELL TYPES ARE NOT RE-ANNOTATED.  `celltype` (fine) and `celltype_merged` are read
    from the existing cache verbatim.
  * FINE-FAMILY FLAGS ONLY.  M1's frozen call thresholds within the FINE `cell_type`
    (`phase2_downstream.py`), and the mouse cache has no `flag_merged_*` key.  No merged
    flags are computed and none are written.
  * SenePy and the seven `A_sender_for_<module>` per-module Tier A sets are NOT scored,
    exactly as on the H1 ISO arm.  The ISO cache keeps the full-panel `senepy_score`,
    `tierApm__*` and `flag_pm_*` keys, unused; do not run those calls against it.
  * The loader's Gene Expression feature list on these samples carries 5,106 names, while
    the pinned mouse PANEL count (5,097) is the metadata-derived panel of
    `h1_a8_crossarm.panels()`.  All 2,435 intersected mouse symbols are present in the
    loader's var on every section, so the restriction is exact either way; the script
    asserts that and reports both numbers.

Usage: python3 code/m1_callers_iso.py <SAMPLE> [...]
Writes data/processed_m1_iso/senders_iso_<SAMPLE>.csv, modules_iso_<SAMPLE>.csv,
       caller_meta_iso_<SAMPLE>.json  (all NEW paths; data/processed/ is read-only here).
"""
import sys, os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 4; sc.settings.verbosity = 0
sys.path.insert(0, "/workspace/code")
import sasp_real as R
import sasp_phase3 as P
import iso_panel

OUT = "/workspace/data/processed_m1_iso"
CACHE_FULL = "/workspace/data/processed/cache3"
GS = "/workspace/genesets/"
EXCLUDE_FROM_STRATA = {"Low_quality", "Unknown"}      # phase2_downstream.py, verbatim


def gl(n):
    return [l.strip() for l in open(GS + n + ".txt") if l.strip()]


A_STRICT = gl("A_SENDER_FINAL_strict")
BMODS = {os.path.basename(p)[2:-4]: [l.strip() for l in open(p) if l.strip()]
         for p in sorted(glob.glob(GS + "B_*.txt"))}
assert sorted(BMODS) == sorted(P.MODULES), "Tier B module list moved"
ISO = iso_panel.mouse_symbols()


def pct_flags(score, labels, q, min_cells=20):
    """`phase2_downstream.py` D-C flag rule == `h1_callers.pct_flags`, verbatim."""
    f = np.zeros(len(score), int)
    for c in pd.unique(labels):
        if c in EXCLUDE_FROM_STRATA:
            continue
        m = (labels == c).to_numpy()
        if m.sum() < min_cells:
            continue
        f[m] = (score[m] > np.nanpercentile(score[m], q)).astype(int)
    return f


def run(sample):
    os.makedirs(OUT, exist_ok=True)
    print("\n" + "=" * 90); print("SAMPLE", sample, "(ORTHOLOG-INTERSECTED PANEL)")
    print("=" * 90, flush=True)
    z = np.load(os.path.join(CACHE_FULL, sample + ".npz"), allow_pickle=False)
    cid = z["cell_id"].astype(str)
    ct = pd.Series(z["celltype"].astype(str))
    ctm = pd.Series(z["celltype_merged"].astype(str))
    z.close()

    M, names, bc = R.load_expression(sample)
    n_loader_panel = len(names)
    pos = pd.Index(bc.astype(str)).get_indexer(pd.Index(cid))
    assert (pos >= 0).all(), sample + ": cache cells missing from the expression matrix"
    B = ad.AnnData(M[pos], obs=pd.DataFrame(index=cid), var=pd.DataFrame(index=names))
    del M
    # ---- THE RESTRICTION.  Before normalisation, before scoring. ----
    present = ISO & set(names)
    assert len(present) == len(ISO), \
        "%s: only %d/%d intersected mouse symbols on the loader var" % (
            sample, len(present), len(ISO))
    B = B[:, [g in ISO for g in B.var_names]].copy()
    assert B.n_vars == len(ISO)
    cdk_on = "Cdkn1a" in set(B.var_names)
    print("loader panel %d (pinned metadata panel 5097) -> %d intersected; "
          "Cdkn1a on intersected panel: %s" % (n_loader_panel, B.n_vars, cdk_on), flush=True)
    assert cdk_on, "Cdkn1a is not on the intersected panel"
    B.layers["counts"] = B.X.copy()
    sc.pp.normalize_total(B); sc.pp.log1p(B)
    n = B.n_obs
    print("cells %d ; fine types %d ; merged types %d" % (n, ct.nunique(), ctm.nunique()),
          flush=True)

    surv = {}
    on = [g for g in A_STRICT if g in B.var_names]
    on_full = [g for g in A_STRICT if g in names]
    surv["A_SENDER_FINAL_strict"] = dict(total=len(A_STRICT), on_full_panel=len(on_full),
                                         on_iso_panel=len(on))
    print("Tier A strict: %d total, %d on full panel, %d on intersected panel"
          % (len(A_STRICT), len(on_full), len(on)))
    sc.tl.score_genes(B, on, score_name="_tierA", ctrl_size=200)
    tierA = B.obs["_tierA"].to_numpy()
    cdk = np.asarray(B.layers["counts"][:, B.var_names.get_loc("Cdkn1a")].todense()).ravel()

    sen = pd.DataFrame({"cell_id": cid, "cell_type": ct.values,
                        "cell_type_merged": ctm.values,
                        "cdkn1a_counts": cdk, "cdkn1a_pos": (cdk > 0).astype(int),
                        "tierA_score": np.round(tierA, 5)})
    for q in (90, 95, 99):
        sen["sender_flag_p%d" % q] = pct_flags(tierA, ct, q)
    sen.to_csv(OUT + "/senders_iso_%s.csv" % sample, index=False)

    mod = pd.DataFrame({"cell_id": cid})
    for name, genes in BMODS.items():
        o = [g for g in genes if g in B.var_names]
        o_full = [g for g in genes if g in names]
        surv["B_" + name] = dict(total=len(genes), on_full_panel=len(o_full),
                                 on_iso_panel=len(o))
        print("  B_%-24s %3d total, %3d full, %3d iso, ctrl_size=%d"
              % (name, len(genes), len(o_full), len(o), max(200, len(o) * 5)))
        sc.tl.score_genes(B, o, score_name="_m", ctrl_size=max(200, len(o) * 5))
        mod[name] = np.round(B.obs["_m"].to_numpy(), 5)
    mod.to_csv(OUT + "/modules_iso_%s.csv" % sample, index=False)

    json.dump({"sample": sample, "panel_loader": int(n_loader_panel),
               "panel_pinned_metadata": 5097, "panel_intersected": int(B.n_vars),
               "n_cells": int(n), "cdkn1a_on_intersected_panel": bool(cdk_on),
               "gene_survival": surv, "senepy": "NOT RUN (declared)",
               "tierApm_A_sender_for_sets": "NOT SCORED (declared)",
               "merged_flags": "NOT COMPUTED (M1's frozen call is the fine family)",
               "sender_p95_fine": int(sen.sender_flag_p95.sum())},
              open(OUT + "/caller_meta_iso_%s.json" % sample, "w"), indent=1)
    print("wrote senders_iso_%s.csv, modules_iso_%s.csv" % (sample, sample), flush=True)
    del B


if __name__ == "__main__":
    for s in (sys.argv[1:] or list(P.IN_BAND)):
        run(s)
