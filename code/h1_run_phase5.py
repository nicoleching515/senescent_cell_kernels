#!/usr/bin/env python3
"""Phase 10 — the Phase-5 producers on the H1 arm: kernel families, superposition vs
nearest, proximal vs downstream.

`run_phase5_kernels.py`, `run_phase5_super.py` and `phase5_common.py` are imported and
their own stage functions are called; nothing is reimplemented.  Three things must be
handled and are handled here rather than by editing those files:

 1. **They do not respect `sasp_phase3.RESULTS`.**  All three take their output directory
    from `phase5_common.RES5 = "/workspace/results/phase5"`, which is the MOUSE directory,
    and `run_phase5_kernels.py` has no `--tag`.  Running them unwrapped on H1 would
    silently overwrite the committed mouse Phase-5 CSVs.  This wrapper rebinds
    `phase5_common.RES5` **and** the already-captured `RES` in each runner module to
    `results/phase10_h1/` before anything is written.
 2. **The `--call` default is `run_phase3_nulls.PRIMARY_CALL = "tierA_p95"`**, the
    frozen-literal fine-label call.  H1's primary is `tierAmg_p95` (PREREG D-B), so the
    call is always passed explicitly and both are run.
 3. **loky workers are fresh interpreters** and would read the mouse cache.  The worker
    functions are wrapped so the arm is rebound inside the child.

Declared, not worked around: `sec.meta["condition"]` and `["week"]` are `"?"` / `NaN` on H1
because `sasp_real.parse_sample` cannot parse `SPLN07`, so the `arm` / `week` columns of
these outputs are meaningless on this arm and must not be read.  `run_phase5_super`'s N3
null is the whole-section BOUNDING-BOX torus shift: Phase 5 never adopted the Phase-7/C1
in-tissue geometry, on either arm, and that is a limitation of the T1 null rows on both.

Usage
  python3 -u code/h1_run_phase5.py --which kernels --stage section  --call tierAmg_p95
  python3 -u code/h1_run_phase5.py --which kernels --stage proxdown --call tierAmg_p95
  python3 -u code/h1_run_phase5.py --which super   --stage section  --call tierAmg_p95
"""
import argparse, sys, time
import pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_phase10                                   # noqa: F401  ARM BINDING first
import phase5_common as C5

RES = h1_phase10.RESULTS10
C5.RES5 = RES                                       # (1) redirect the output directory
import run_phase5_kernels as K                      # noqa: E402
import run_phase5_super as S                        # noqa: E402
K.RES = RES
S.RES = RES
import sasp_phase3 as P                             # noqa: E402
import h1_common as H                               # noqa: E402


def _bind():
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                # noqa: F401
    import phase5_common as _C5
    _C5.RES5 = "/workspace/results/phase10_h1"


def _w_kern(s, call, seed):
    _bind()
    import run_phase5_kernels as _K
    return _K._section_job(s, call, seed)


def _w_super(s, call, seed):
    _bind()
    import run_phase5_super as _S
    return _S._section_job(s, call, seed)


def _w_super_null(s, call, seed, ndraw):
    _bind()
    import run_phase5_super as _S
    return _S._null_job(s, call, seed, ndraw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", required=True, choices=["kernels", "super"])
    ap.add_argument("--stage", required=True,
                    choices=["section", "heldout", "proxdown", "nulls"])
    ap.add_argument("--call", default="tierAmg_p95")
    ap.add_argument("--n-jobs", type=int, default=7)
    ap.add_argument("--n-draw", type=int, default=5)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    from joblib import Parallel, delayed
    secs = list(H.ALL_SECTIONS)
    tag = a.tag or ("" if a.call == "tierAmg_p95" else "_" + a.call)
    t0 = time.time()
    print("which=%s stage=%s call=%s -> %s (tag %r)" % (a.which, a.stage, a.call, RES, tag),
          flush=True)

    if a.which == "kernels":
        if a.stage == "section":
            out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
                delayed(_w_kern)(s, a.call, P.MASTER_SEED + 300 * i)
                for i, s in enumerate(secs))
            df = pd.DataFrame([r for rs in out for r in rs])
            df["call"] = a.call
            df.to_csv("%s/kernel_families%s.csv" % (RES, tag), index=False)
        elif a.stage == "heldout":
            df = K.stage_heldout(secs, a.call, a.n_jobs)
            df["call"] = a.call
            df.to_csv("%s/kernel_heldout%s.csv" % (RES, tag), index=False)
        elif a.stage == "proxdown":
            df = K.stage_proxdown(secs, a.call, a.n_jobs)
            df["call"] = a.call
            df.to_csv("%s/proximal_vs_downstream%s.csv" % (RES, tag), index=False)
        else:
            raise SystemExit(a.stage)
    else:
        if a.stage == "section":
            out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
                delayed(_w_super)(s, a.call, P.MASTER_SEED + 500 * i)
                for i, s in enumerate(secs))
            df = pd.DataFrame([r for rs in out for r in rs])
            df["call"] = a.call
            df.to_csv("%s/super_section%s.csv" % (RES, tag), index=False)
        elif a.stage == "nulls":
            out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
                delayed(_w_super_null)(s, a.call, P.MASTER_SEED + 700 * i, a.n_draw)
                for i, s in enumerate(secs))
            df = pd.DataFrame([r for rs in out for r in rs])
            df["call"] = a.call
            df.to_csv("%s/super_nulls%s.csv" % (RES, tag), index=False)
        elif a.stage == "heldout":
            df = S.stage_heldout(secs, a.call, a.n_jobs)
            df["call"] = a.call
            df.to_csv("%s/super_heldout%s.csv" % (RES, tag), index=False)
        else:
            raise SystemExit(a.stage)
    print("done in %.1f min" % ((time.time() - t0) / 60), flush=True)


if __name__ == "__main__":
    main()
