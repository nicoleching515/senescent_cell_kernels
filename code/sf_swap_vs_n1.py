#!/usr/bin/env python3
"""Phase 8 / C1 -- is N3-swap the same operation as N1, for every sender call?

Reproduces `results/phase3/sf_summary_c1_swap_vs_n1.csv`, the table behind the
"N3-swap == N1" verdict in `reports/CS_PHASE8_C1_CLOSEOUT.md` (section 3) and
`reports/AUDIT_PHASE8_FACTCHECK.md:425`.  One row per sender call:

    N1            median surviving fraction of the within-celltype label
                  permutation (`perm_nulls{,_n7}.csv`)
    N3_swap       median SF of the C1 null that relocates senders to random
                  REAL cell positions (`perm_nulls_c1{,_n7}.csv`, scope=full)
    rho           Spearman rho between the two, computed PER FIT -- the medians
                  can agree while the fits disagree, so the correlation is the
                  claim, not the medians
    med_absdiff   median |N1_sf - N3_swap_sf| over the same fits
    N3_swap_full  the same SF with beta under the full N5+N6+zonation design;
                  only the primary call was run that way, so it is blank for
                  the five N7 calls

Usage (results dir is positional, defaults to the post-C6 tree):

    /workspace/envs/sasp311/bin/python /workspace/code/sf_swap_vs_n1.py \
        /workspace/results/phase3_pre_c6 > out.csv

RECONSTRUCTION.  Written on 2026-08-27 for the reproducibility repair recorded
in `reports/AUDIT_REPRODUCIBILITY.md` ("`swap_vs_n1` occurs nowhere in
`code/`").  The COMMITTED `results/phase3/sf_summary_c1_swap_vs_n1.csv` is the
PRE-C6 vintage -- it was never regenerated after the C6 gene-set change -- so
the byte-exact verification was run against `results/phase3_pre_c6`, with the
command above, and matched all six rows exactly.  Point it at
`/workspace/results/phase3` to get the post-C6 numbers, which differ.

The reportable population is imported from `summarize_phase3_c1.reportable`
rather than re-implemented, so this table is over exactly the fits the rest of
Phase 3 reports.  That function reads `main_fits.csv` out of its own module
global `RES`, which is pinned to the post-C6 tree, so the results directory is
injected by rebinding `SC1.RES` before calling it.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "/workspace/code")
import summarize_phase3_c1 as SC1

KEY = ["section", "celltype", "module"]

# Tier A percentiles first: the pre-registered axis, and where the Phase 7
# primary sender definition lives.  The two ad-hoc caller definitions follow.
CALLS = ["tierA_p90", "tierA_p95", "tierA_p99",
         "cdkn1a_pos", "senepy_p95", "senepy_p99"]

# tierA_p95 is the primary call; its nulls live in the primary files, the other
# five in the N7 sender-definition axis files.
PRIMARY = "tierA_p95"


def prevalence(res):
    """min/max sender prevalence (%) per call over the six in-band sections."""
    d = pd.read_csv(f"{res}/poisson_density.csv")
    g = d[d.band == "in_band"].groupby("call").sender_prevalence
    return (g.min() * 100).to_dict(), (g.max() * 100).to_dict()


def row(res, call, lo, hi):
    suf = "" if call == PRIMARY else "_n7"
    n1 = pd.read_csv(f"{res}/perm_nulls{suf}.csv")
    c1 = pd.read_csv(f"{res}/perm_nulls_c1{suf}.csv")
    n1 = n1[n1.call == call]
    # scope=="tile" rows fit the same null on solid-tissue tiles only; N1 has
    # no tile counterpart, so the comparison is whole-section throughout.
    c1 = c1[(c1.call == call) & (c1.scope == "full")]
    cols = [c for c in ("N3_swap_sf", "N3_swap_full_sf") if c in c1]
    m = (SC1.reportable(call)
         .merge(n1[KEY + ["N1_sf"]], on=KEY, how="inner")
         .merge(c1[KEY + cols], on=KEY, how="inner"))
    return dict(
        call=call,
        n=len(m),
        prevalence_pct="%.1f-%.1f" % (lo[call], hi[call]),
        N1=round(m.N1_sf.median(), 3),
        N3_swap=round(m.N3_swap_sf.median(), 3),
        rho=round(spearmanr(m.N1_sf, m.N3_swap_sf).statistic, 3),
        med_absdiff=round((m.N1_sf - m.N3_swap_sf).abs().median(), 4),
        N3_swap_full=(round(m.N3_swap_full_sf.median(), 3)
                      if "N3_swap_full_sf" in m else np.nan))


def main():
    res = sys.argv[1] if len(sys.argv) > 1 else "/workspace/results/phase3"
    SC1.RES = res                       # see the module docstring
    lo, hi = prevalence(res)
    df = pd.DataFrame([row(res, c, lo, hi) for c in CALLS])
    dest = sys.argv[2] if len(sys.argv) > 2 else sys.stdout
    df.to_csv(dest, index=False)
    return df


if __name__ == "__main__":
    main()
