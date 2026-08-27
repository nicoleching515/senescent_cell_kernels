#!/usr/bin/env python3
"""Phase 9 — DeepScence at its PUBLISHED default, `denoise=True`, on H1.

The H1 analogue of `code/run_deepscence_dca.py`.  It differs from that file in exactly two
places and both are the point of §8's free experiment:
  * the loader is `h1_common.load_matrix` -- the NATIVE human panel, no ortholog remap;
  * outputs go to data/processed_h1/.
Everything else -- the isolated Python 3.8 / TensorFlow 2.4.4 DCA venv reached through
`code/_shims_dca_bridge`, `random_state`, the >= 20-counts filter, the published CDKN1A
anchor -- is unchanged.  The main pinned 3.11 stack never imports TensorFlow.

This tests PREREG §8 predictions **P-vi** (denoise=True RAISES the score's depth loading, as
it does on 3 of 3 M1 sections) and **P-vii** (the denoise=True seed instability recurs:
>= 1 of 3 seeds gives a top-5 % sender set with Jaccard < 0.30 against the others).

Env: DCA_VENV_PYTHON, DCA_BRIDGE_SCRATCH, DCA_THREADS  (as for the mouse runner)
Usage: python3 code/h1_deepscence_dca.py [--subsample N] [--seed S] SPLN21 [...]
"""
import sys, os, time, json, argparse
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/workspace/code/_shims_dca_bridge")
import dca.api as _dca_api                                          # noqa: E402
assert "dca_bridge" in open(_dca_api.__file__).read(), \
    "the raising stub, not the bridge, is on sys.path -- refusing to run"

sys.path.insert(0, "/workspace/code")
import h1_common as H                                               # noqa: E402
import torch                                                        # noqa: E402
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "16")))
import DeepScence.api as api                                        # noqa: E402


def load_section(section):
    X, names, bc, _ = H.load_matrix(section, "gene")
    keep = np.asarray(X.sum(1)).ravel() >= 20
    A = ad.AnnData(X[keep], obs=pd.DataFrame(index=bc[keep]),
                   var=pd.DataFrame(index=names))
    A.var_names_make_unique()
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subsample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--denoise-false", action="store_true",
                    help="same code path and the SAME fixed subsample, with denoise=False. "
                         "This is how the seed-to-seed FLOOR is measured on identical cells; "
                         "without it a denoise=True vs denoise=False comparison is confounded "
                         "by which cells were drawn.")
    ap.add_argument("sections", nargs="+")
    a = ap.parse_args()
    for section in a.sections:
        stem = "nodn" if a.denoise_false else "dca"
        suffix = (stem if not a.subsample else "%s_sub%d" % (stem, a.subsample)) \
                 + ("" if a.seed == 0 else "_seed%d" % a.seed)
        out = H.PROC + "/deepscence_h1_%s_%s.csv" % (suffix, section)
        if os.path.exists(out):
            print("%s [%s]: already present, skipping" % (section, suffix), flush=True)
            continue
        A = load_section(section)
        if a.subsample and a.subsample < A.n_obs:
            # fixed subsampling seed, independent of --seed, so every seed sees the SAME
            # cells -- identical to the mouse runner
            idx = np.sort(np.random.default_rng(12345).choice(A.n_obs, a.subsample,
                                                              replace=False))
            A = A[idx].copy()
        ids = A.obs_names.to_numpy()
        counts = np.asarray(A.X.sum(1)).ravel().astype(np.int64)
        print("%s [%s]: %d cells x %d NATIVE human genes"
              % (section, suffix, A.n_obs, A.n_vars), flush=True)
        t = time.time()
        res = api.DeepScence(A, denoise=not a.denoise_false, verbose=False,
                             random_state=a.seed)
        mins = (time.time() - t) / 60
        pd.DataFrame({"cell_id": ids,
                      "deepscence_score": np.round(res.obs["ds"].to_numpy(), 5),
                      "counts": counts}).to_csv(out, index=False)
        log = {k: (v.tolist() if hasattr(v, "tolist") else v)
               for k, v in dict(res.uns["log"]).items()}
        json.dump(dict(section=section, config=suffix, n_cells=int(len(ids)),
                       n_genes=int(res.n_vars), denoise=not a.denoise_false, random_state=a.seed,
                       anchor="published_CDKN1A", panel="native human, no ortholog remap",
                       min_counts_per_cell=20, direction_log=log,
                       dca_bridge=res.uns.get("dca_bridge") or dict(_dca_api.LAST_RUN),
                       minutes=round(mins, 2)),
                  open(H.PROC + "/deepscence_meta_h1_%s_%s.json" % (suffix, section), "w"),
                  indent=1, default=str)
        print("wrote %s (%d cells, %.1f min)" % (os.path.basename(out), len(ids), mins),
              flush=True)


if __name__ == "__main__":
    main()
