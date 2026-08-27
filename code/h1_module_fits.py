#!/usr/bin/env python3
"""Phase 9 — the biological-module reference fits on H1.

A7 asks whether the control response is FLAT.  "Flat" is only falsifiable against two
things: the amplitude the same estimator returns for the seven Tier B modules on the same
cells, and the smallest amplitude a single fit could have resolved.  `summarize_a7.py`
supplies both from `results/phase3/main_fits.csv` on the mouse arm; H1 has no such file until
Phase 10, so this script produces the equivalent.

It does not reimplement anything: it calls `run_phase3_nulls._section_job` — the frozen
stage-2 driver — with the H1 cache in place.  Every parameter is the frozen one.
`strata_mode` is 'none', because the mouse 'zonation' stratification loops over the three
hepatocyte zonation labels and has no H1 counterpart at this stage.

Usage: python3 code/h1_module_fits.py [--calls tierA_p95,cdkn1a_pos] [--n-jobs 7]
Writes results/phase9_h1/h1_module_fits.csv
"""
import sys, argparse
import pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_sec                     # noqa: F401
from joblib import Parallel, delayed
import sasp_phase3 as P
import run_phase3_nulls as RN
import h1_common as H

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default="tierA_p95,cdkn1a_pos")
    ap.add_argument("--sections", default="all")
    ap.add_argument("--n-jobs", type=int, default=7)
    a = ap.parse_args()
    secs = list(H.ALL_SECTIONS) if a.sections == "all" else a.sections.split(",")
    calls = a.calls.split(",")
    jobs = [(s, c, P.MASTER_SEED + 1000 * i + j)
            for i, s in enumerate(secs) for j, c in enumerate(calls)]
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(RN._section_job)(s, c, sd, "none") for s, c, sd in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    df.to_csv(H.RESULTS + "/h1_module_fits.csv", index=False)
    print(df.shape, "->", H.RESULTS + "/h1_module_fits.csv")
    if "beta_naive" in df:
        rep = df[(df.beta_naive > 0) & (df.get("beta_base_lo", 0) > 0)]
        print("fits %d, reportable %d (%.1f%%)" % (len(df), len(rep),
                                                   100 * len(rep) / max(len(df), 1)))
