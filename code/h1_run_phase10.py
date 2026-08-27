#!/usr/bin/env python3
"""Phase 10 — run H1 through the FROZEN Phase-3 pipeline.

Every estimator function is imported from `run_phase3_nulls` / `run_phase3_var` and called,
never copied.  The only H1-specific choices are declared here and nowhere else:

  * the arm binding is `code/h1_phase10.py` (cache, section list, results dir, and the one
    additional `tierAmg_pNN` call name);
  * `strata_mode = "none"` -- the mouse `zonation` stratification loops over the three
    hepatocyte zonation labels and has no H1 counterpart (same choice Phase 9 made in
    `code/h1_module_fits.py`);
  * `CURVE_TYPE = "B cells"` replaces the mouse `curves_for = [("Hepatocytes", j)]`.  B cells
    are realised in all 7 sections at 13.5-21.3 % (`CS_PHASE9_H1_AUDIT.md` §9.1); Hepatocytes
    do not exist on this arm, so the frozen literal would emit no curves at all.

SENDER CALLS.  PI decision, post-freeze and declared as such: on H1 the PRIMARY call is
`tierAmg_p95` -- the identical Tier A percentile rule evaluated at the MERGED label family,
the family the estimator actually stratifies receivers on -- and the frozen-literal
fine-label `tierA_p95` is retained as the sensitivity.  Both are computed everywhere and
reported side by side.  Reason: deviation H5, `CS_PHASE9_H1_AUDIT.md` §9.4 -- the frozen call
thresholds within FINE types while `sasp_phase3.LABELS = "merged"`, so 0-13.6 % of cells per
section are eligible receivers that can never be called senders, which zeroed the T/NK sender
set in SPLN43 (0 senders in 24,815 cells).  No threshold is tuned; only the label family the
identical rule is evaluated at.

Usage
  python3 -u code/h1_run_phase10.py --stage window
  python3 -u code/h1_run_phase10.py --stage main    --calls all12 --n-jobs 10
  python3 -u code/h1_run_phase10.py --stage perm    --calls primary2 --n-perm 1000 --n-jobs 7
  python3 -u code/h1_run_phase10.py --stage perm_c1 --calls primary2 --n-perm 1000 --n-jobs 7
  python3 -u code/h1_run_phase10.py --stage var     --calls primary2 --n-perm 1000 --n-jobs 7
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_phase10                              # noqa: F401  ARM BINDING -- must be first
from joblib import Parallel, delayed
import sasp_phase3 as P
import run_phase3_nulls as RN
import h1_common as H

RES = h1_phase10.RESULTS10
CURVE_TYPE = "B cells"
PRIMARY_H1 = "tierAmg_p95"          # PI decision (see docstring)
FROZEN_LITERAL = "tierA_p95"        # the frozen-literal sensitivity
MERGED_CALLS = list(h1_phase10.MERGED_CALLS)
ALL12 = list(RN.ALL9_CALLS) + MERGED_CALLS


def resolve(spec):
    if spec == "all12":
        return list(ALL12)
    if spec == "primary2":
        return [PRIMARY_H1, FROZEN_LITERAL]
    if spec == "all9":
        return list(RN.ALL9_CALLS)
    if spec == "merged":
        return list(MERGED_CALLS)
    return spec.split(",")


# --- workers: loky gives a FRESH interpreter, so rebind the arm inside it ----
def _w_main(sample, call, seed, module):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                          # noqa: F401
    import run_phase3_nulls as _RN
    return _RN._section_job(sample, call, seed, "none", module)


def _w_perm(sample, call, seed, n_perm, module):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                          # noqa: F401
    import sasp_phase3 as _P
    import run_phase3_nulls as _RN
    cf = [(CURVE_TYPE, j) for j in range(len(_P.MODULES))]
    return _RN._perm_job(sample, call, seed, n_perm, True, cf, module)


def _w_perm_c1(sample, call, seed, n_perm, module):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                          # noqa: F401
    import sasp_phase3 as _P
    import run_phase3_nulls as _RN
    cf = [(CURVE_TYPE, j) for j in range(len(_P.MODULES))]
    return _RN._perm_c1_job(sample, call, seed, n_perm, True, cf, module)


def _w_var(sample, call, seed, n_perm):
    import sys as _s; _s.path.insert(0, "/workspace/code")
    import h1_phase10                                          # noqa: F401
    import run_phase3_var as _RV
    return _RV._var_job(sample, call, seed, n_perm, True)


def stage_window(sections, calls):
    rows = []
    for s in sections:
        sec = P.Sec(s)
        for call in calls:
            for module in (P.MODULES if RN.is_permodule(call) else [None]):
                try:
                    snd = sec.sender_mask(call, module=module)
                except Exception as e:
                    print("  %s %s %s: %s" % (s, call, module or "", e)); continue
                if snd.sum() < 50:
                    continue
                d = P.dist_to_senders(sec.coords.astype(float), snd)
                recv = (~np.isin(sec.celltype, P.EXCLUDE_TYPES)) & (~snd)
                dr = d[recv]
                rows.append(dict(section=s, call=call, sender_module=module or "",
                                 n=int(recv.sum()), n_senders=int(snd.sum()),
                                 prevalence=float(snd.mean()),
                                 median_nn_um=float(sec.median_nn_um),
                                 **{("d_p%d" % q): float(np.percentile(dr, q))
                                    for q in (50, 75, 90, 95, 99)},
                                 d_max=float(dr.max()),
                                 frac_gt_80=float((dr > 80).mean()),
                                 frac_gt_100=float((dr > 100).mean()),
                                 frac_gt_150=float((dr > 150).mean())))
    df = pd.DataFrame(rows)
    df.to_csv(RES + "/window.csv", index=False)
    print(df.to_string(index=False))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["window", "main", "perm", "perm_c1", "var", "diag"])
    ap.add_argument("--calls", default="primary2")
    ap.add_argument("--sections", default="all")
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    secs = list(H.ALL_SECTIONS) if a.sections == "all" else a.sections.split(",")
    calls = resolve(a.calls)
    print("stage=%s sections=%d calls=%s -> %s" % (a.stage, len(secs), calls, RES),
          flush=True)
    t0 = time.time()

    if a.stage == "window":
        stage_window(secs, ALL12); return

    if a.stage == "main":
        jobs = RN._expand(secs, calls, P.MASTER_SEED, 1000, 1)
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_w_main)(s, c, sd, m) for s, c, sd, m in jobs)
        df = pd.DataFrame([r for rs in out for r in rs])
        df.to_csv(RES + "/main_fits%s.csv" % a.tag, index=False)
        rep = df[(df.beta_naive > 0) & (df.get("beta_base_lo", 0) > 0)]
        print(df.shape, "reportable", len(rep), flush=True)

    elif a.stage in ("perm", "perm_c1"):
        w = _w_perm if a.stage == "perm" else _w_perm_c1
        jobs = RN._expand(secs, calls, P.MASTER_SEED, 5000, 17)
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(w)(s, c, sd, a.n_perm, m) for s, c, sd, m in jobs)
        suff = "" if a.stage == "perm" else "_c1"
        pd.DataFrame([r for rs, _ in out for r in rs]).to_csv(
            RES + "/perm_nulls%s%s.csv" % (suff, a.tag), index=False)
        pd.DataFrame([r for _, cs in out for r in cs]).to_csv(
            RES + "/perm_curves%s%s.csv" % (suff, a.tag), index=False)

    elif a.stage == "var":
        jobs = RN._expand(secs, calls, P.MASTER_SEED, 5000, 17)
        out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
            delayed(_w_var)(s, c, sd, a.n_perm) for s, c, sd, _m in jobs)
        pd.DataFrame([r for rs, _ in out for r in rs]).to_csv(
            RES + "/perm_nulls_var%s.csv" % a.tag, index=False)
        pd.DataFrame([r for _, ds in out for r in ds]).to_csv(
            RES + "/perm_draws_var%s.csv" % a.tag, index=False)

    elif a.stage == "diag":
        import phase3_null_diag as ND
        import run_phase3_var as RV
        rows = []
        for s in secs:
            rows += ND.section_rows(s, call=calls[0])
        pd.DataFrame(rows).to_csv(RES + "/null_destructiveness%s.csv" % a.tag, index=False)
        rows = []
        for s in secs:
            rows += RV._diag_job(s, calls[0], 20, P.MASTER_SEED)
        pd.DataFrame(rows).to_csv(RES + "/null_destructiveness_var%s.csv" % a.tag,
                                  index=False)

    print("stage %s done in %.1f min" % (a.stage, (time.time() - t0) / 60), flush=True)


if __name__ == "__main__":
    main()
