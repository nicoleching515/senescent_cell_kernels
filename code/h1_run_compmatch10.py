#!/usr/bin/env python3
"""Phase 10 / roadmap 10.2 — the composition-matched rerun on the H1 arm.

`code/run_phase8_compmatch.py` is the producer and is NOT reimplemented here.  This file
exists for one reason: its `_job` is dispatched to loky workers, which are fresh
interpreters that import `run_phase8_compmatch` WITHOUT the H1 arm binding and would read
the mouse cache.  `_w` rebinds the arm inside the worker and then calls the frozen `_job`.

PI decision D15 is unchanged: the covariate-adjusted counterpart (`type_adj`, `comp_adj`,
`typecomp_adj`) is PRIMARY for any claim about how much of the gradient is composition, and
the matched-decoy protocol (`comp`, with `full` as the published-N2 comparison) is reported
ALONGSIDE every time -- its inertness is itself the finding (PREREG §10.8, P23/P25).
The five seeds are the frozen 20260901-05.

Usage:
  SASP_H1_UNFROZEN=1 python3 -u code/h1_run_compmatch10.py --calls tierAmg_p95,tierA_p95 \
      --n-jobs 8 --out-tag _h1
"""
import argparse, os, sys, time
import pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_phase10                              # noqa: F401  ARM BINDING -- must be first
from joblib import Parallel, delayed
import run_phase8_compmatch as CM

RES = h1_phase10.RESULTS10
CM.RES = RES


def _w(arm, cfg, sample, call, module, variants, seeds):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                          # noqa: F401
    import run_phase8_compmatch as _CM
    return _CM._job(arm, cfg, sample, call, module, variants, seeds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--sections", default=None)
    ap.add_argument("--calls", default=None)
    ap.add_argument("--modules", default=None)
    ap.add_argument("--variants", default=",".join(CM.VARIANTS))
    ap.add_argument("--seeds", default=",".join(str(s) for s in CM.COMPMATCH_SEEDS))
    ap.add_argument("--out-tag", default="_h1")
    a = ap.parse_args()

    cfg = CM.arm_config("h1")                 # the gate is exercised, not bypassed
    sections = (a.sections.split(",") if a.sections else cfg["sections"])
    sections = [s for s in sections
                if os.path.exists(os.path.join(cfg["cache"], "%s.npz" % s))]
    calls = (a.calls.split(",") if a.calls
             else [cfg["primary_call"], cfg["permodule_call"]])
    modules = a.modules.split(",") if a.modules else cfg["modules"]
    variants = a.variants.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    print("arm h1 | sections %s | calls %s | variants %s | seeds %s"
          % (sections, calls, variants, seeds), flush=True)

    jobs = CM.build_jobs("h1", cfg, sections, calls, modules, variants, seeds)
    t0 = time.time()
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(_w)(*j) for j in jobs)
    fits = pd.DataFrame([r for rs in out for r in rs]).rename(columns=CM.RENAME)
    fp = "%s/compmatch_fits%s.csv" % (RES, a.out_tag)
    sp = "%s/compmatch_reruns%s.csv" % (RES, a.out_tag)
    for p in (fp, sp):
        if os.path.exists(p):
            raise SystemExit("%s exists; refusing to overwrite. Use --out-tag" % p)
    fits.to_csv(fp, index=False)
    summ = CM.summarise(fits)
    summ.to_csv(sp, index=False)
    print("\nwrote %s %s\nwrote %s %s\nelapsed %.0fs"
          % (fp, fits.shape, sp, summ.shape, time.time() - t0))
    show = ["row_type", "call", "variant", "scope_kind", "seed", "n_reportable",
            "median_sf_matched", "median_comp_share", "median_match_rate",
            "median_max_smd_after"]
    with pd.option_context("display.width", 220):
        print(summ[summ.row_type == "summary"][show].to_string(index=False))


if __name__ == "__main__":
    main()
