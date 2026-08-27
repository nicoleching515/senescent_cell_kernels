#!/usr/bin/env python3
"""Phase 9 Job B step 2 — DeepScence on H1, NATIVE, no ortholog remapping.

This is §8's free experiment (PREREG_PHASE8.md §8): CoreScence is a human gene set, so on M1
it runs ortholog-remapped (4,845 of 5,097 panel genes) and on H1 it runs exactly as published.

Settings are the frozen ones (PREREG §3.9), identical to `code/run_deepscence_all.py`:
  denoise=False (PRIMARY, a chosen value), random_state=0, published CDKN1A anchor,
  >= 20 counts/cell.  The one deliberate difference is the panel: NATIVE human, no remap.

MEMORY.  DeepScence coerces X to dense and keeps `original = adata.copy()` alongside, so a
section costs several dense (n_cells x n_genes) float32 arrays at once.  M1's 83k-cell
sections needed ~16 GB of cgroup headroom; H1's are 220k-396k cells.  Sections are therefore
run ONE AT A TIME, largest first, and this script refuses to start a section whose projected
requirement exceeds the free headroom unless --subsample is given.  A section that does not
fit is recorded as not fitting; that is a reproducible result, not a failure to be hidden.

Usage: python3 code/h1_deepscence.py SPLN24 [--subsample N] [--headroom-factor F]
Writes data/processed_h1/deepscence_h1_<sec>.csv (+ _subsampled marker in the meta json).
"""
import sys, os, time, json, argparse
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, "/workspace/code")
import h1_common as H

import torch
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "32")))
sys.path.insert(0, "/workspace/code/_shims")
import DeepScence.api as api

CG = "/sys/fs/cgroup/"


def headroom_gb():
    cur = int(open(CG + "memory.current").read())
    mx = int(open(CG + "memory.max").read())
    return (mx - cur) / 2 ** 30, mx / 2 ** 30


def run(section, subsample=None, factor=5.0, force=False):
    out = H.PROC + "/deepscence_h1_%s.csv" % section
    if os.path.exists(out) and not force:
        print("%s: already present, skipping" % section); return
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = np.asarray(X.sum(1)).ravel() >= 20
    A = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]), var=pd.DataFrame(index=names))
    A.var_names_make_unique()
    free, mx = headroom_gb()
    need = factor * A.n_obs * A.n_vars * 4 / 2 ** 30
    print("%s: %d cells x %d NATIVE human genes | projected %.1f GB (%.1f x dense) | "
          "cgroup free %.1f of %.1f GB" % (section, A.n_obs, A.n_vars, need, factor, free, mx),
          flush=True)
    subsampled = False
    if subsample and A.n_obs > subsample:
        rng = np.random.default_rng(H.MASTER_SEED)
        sel = np.sort(rng.choice(A.n_obs, subsample, replace=False))
        A = A[sel].copy(); subsampled = True
        print("  SUBSAMPLED to %d cells (seed %d)" % (A.n_obs, H.MASTER_SEED), flush=True)
    elif need > free:
        print("  *** DOES NOT FIT: %.1f GB projected against %.1f GB free. "
              "Re-run with --subsample N. ***" % (need, free), flush=True)
        json.dump(dict(section=section, status="does_not_fit", n_cells=int(A.n_obs),
                       n_genes=int(A.n_vars), projected_gb=round(need, 1),
                       free_gb=round(free, 1)),
                  open(H.PROC + "/deepscence_meta_h1_%s.json" % section, "w"), indent=1)
        return
    ids = A.obs_names.to_numpy()
    t = time.time()
    res = api.DeepScence(A, denoise=False, verbose=False, random_state=0)
    mins = (time.time() - t) / 60
    pd.DataFrame({"cell_id": ids,
                  "deepscence_score": np.round(res.obs["ds"].to_numpy(), 5)}
                 ).to_csv(out, index=False)
    peak = int(open(CG + "memory.peak").read()) / 2 ** 30 if os.path.exists(CG + "memory.peak") else None
    json.dump(dict(section=section, status="ok", n_cells=int(A.n_obs),
                   n_genes=int(A.n_vars), subsampled=subsampled,
                   subsample_n=int(subsample) if subsampled else None,
                   seed=H.MASTER_SEED if subsampled else None,
                   minutes=round(mins, 1), denoise=False, random_state=0,
                   anchor="published CDKN1A", panel="native human, no ortholog remap",
                   log={k: (v.tolist() if hasattr(v, "tolist") else v)
                        for k, v in res.uns.get("log", {}).items()},
                   cgroup_peak_gb=round(peak, 1) if peak else None),
              open(H.PROC + "/deepscence_meta_h1_%s.json" % section, "w"), indent=1)
    print("wrote %s (%d cells, %.1f min)" % (os.path.basename(out), A.n_obs, mins), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sections", nargs="+")
    ap.add_argument("--subsample", type=int, default=None)
    ap.add_argument("--headroom-factor", type=float, default=5.0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    for s in a.sections:
        run(s, a.subsample, a.headroom_factor, a.force)
