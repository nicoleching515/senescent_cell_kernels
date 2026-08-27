#!/usr/bin/env python3
"""Phase 10 ISO -- run the primary fits on BOTH ortholog-intersected caches.

Every estimator is imported from `run_phase3_nulls` and called, never copied.  The only
choices made here are the arm binding, the section/call lists, the stratification mode, and
the seed basis:

  H1: 7 sections x {tierAmg_p95, tierA_p95}, strata_mode="none"   (as `h1_run_phase10.py`)
  M1: 6 in-band sections x {tierA_p95},      strata_mode="zonation" (M1's frozen mode)

SEEDS -- a declared judgement call.  The brief asks for
`RN._expand(sections, calls, P.MASTER_SEED, 1000, 1)` "identical to the full-panel run so
the N2 match draw matches".  `_expand`'s seed is `MASTER_SEED + 1000*i + 1*j` where i is the
section's index in the sections list and j the call's index in the CALLS LIST PASSED IN, so
the seed a call gets depends on how many other calls were run beside it.  The full-panel
runs used longer call lists (verified from their own logs):

  results/phase10_h1/main_fits.csv : sections = h1_common.ALL_SECTIONS (7),
      calls = run_phase3_nulls.ALL9_CALLS + h1_phase10.MERGED_CALLS (12)  -> tierA_p95 j=1,
      tierAmg_p95 j=10                                        (logs/phase10/main.log)
  results/phase3/main_fits.csv     : sections = sasp_phase3.ALL_SECTIONS (11),
      calls = run_phase3_nulls.ALL9_CALLS (9)                 -> tierA_p95 j=1
                                                               (logs/m1_main.log)

Passing only the two calls of interest would therefore give tierAmg_p95 j=0 and change its
seed, and would change M1's tierA_p95 seed by 1 -- i.e. the literal reading of the formula
CONTRADICTS its own parenthetical.  DEFAULT (`--seed-basis fullrun`) resolves it in favour
of the parenthetical: `_expand` is called with the full-panel run's own section and call
lists and the result is filtered to the calls we need, so every seed is byte-identical to
the full-panel fit it is compared against and the N2 decoy draw is the matched one.
`--seed-basis literal` gives the other reading; the seed actually used is recorded per fit
in `seeds_used.csv` under the results directory, so which was used is never in doubt.

Usage
  python3 -u code/run_phase10_iso.py --arm h1 --n-jobs 4
  python3 -u code/run_phase10_iso.py --arm m1 --n-jobs 4
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
from joblib import Parallel, delayed

H1_CALLS = ["tierAmg_p95", "tierA_p95"]
M1_CALLS = ["tierA_p95"]


# --- workers: loky gives a FRESH interpreter, so rebind the arm inside it ----
def _w_h1(sample, call, seed, module):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10_iso                                       # noqa: F401
    import run_phase3_nulls as _RN
    return _RN._section_job(sample, call, seed, "none", module)


def _w_m1(sample, call, seed, module):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import m1_phase10_iso                                       # noqa: F401
    import run_phase3_nulls as _RN
    return _RN._section_job(sample, call, seed, "zonation", module)


def _jobs(arm, basis):
    import run_phase3_nulls as RN
    import sasp_phase3 as P
    if arm == "h1":
        import h1_common as H
        secs, calls = list(H.ALL_SECTIONS), list(H1_CALLS)
        full_secs = list(H.ALL_SECTIONS)
        full_calls = list(RN.ALL9_CALLS) + ["tierAmg_p90", "tierAmg_p95", "tierAmg_p99"]
    else:
        secs, calls = list(P.IN_BAND), list(M1_CALLS)
        full_secs = list(P.ALL_SECTIONS)
        full_calls = list(RN.ALL9_CALLS)
    if basis == "literal":
        jobs = RN._expand(secs, calls, P.MASTER_SEED, 1000, 1)
    else:
        allj = RN._expand(full_secs, full_calls, P.MASTER_SEED, 1000, 1)
        want = set(secs), set(calls)
        jobs = [j for j in allj if j[0] in want[0] and j[1] in want[1]]
    assert len(jobs) == len(secs) * len(calls), \
        "expected %d jobs, got %d" % (len(secs) * len(calls), len(jobs))
    return secs, calls, jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=["h1", "m1"])
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--seed-basis", default="fullrun", choices=["fullrun", "literal"])
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    assert a.n_jobs <= 4, "resource budget: at most 4 workers"

    if a.arm == "h1":
        import h1_phase10_iso as BIND                            # noqa: F401
        worker, res = _w_h1, BIND.RESULTS10_ISO
    else:
        import m1_phase10_iso as BIND                            # noqa: F401
        worker, res = _w_m1, BIND.RESULTS10_ISO
    import sasp_phase3 as P
    secs, calls, jobs = _jobs(a.arm, a.seed_basis)
    print("arm=%s cache=%s -> %s" % (a.arm, P.CACHE3, res))
    print("sections=%d calls=%s seed_basis=%s jobs=%d"
          % (len(secs), calls, a.seed_basis, len(jobs)), flush=True)
    pd.DataFrame([dict(section=s, call=c, seed=sd, module=m or "")
                  for s, c, sd, m in jobs]).to_csv(
        res + "/seeds_used%s.csv" % a.tag, index=False)

    t0 = time.time()
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(worker)(s, c, sd, m) for s, c, sd, m in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(res + "/main_fits%s.csv" % a.tag, index=False)
    rep = df[(df.beta_naive > 0) & (df.sf_base.notna()) & (df.beta_base_lo > 0)]
    print(df.shape, "reportable(all strata)", len(rep))
    print("done in %.1f min -> %s/main_fits%s.csv"
          % ((time.time() - t0) / 60, res, a.tag), flush=True)


if __name__ == "__main__":
    main()
