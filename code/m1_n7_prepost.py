#!/usr/bin/env python3
"""Phase 8 / M1 -- the N7 sender-definition axis, pre-C6 beside post-C6.

Reproduces `results/phase3/m1_n7_prepost.txt`, the table quoted at
`reports/CORRECTIONS.md:725` (section 10, "The N7 sender-definition axis under
the corrected nulls, at 1,000 permutations").  One block per results tree:

    === pre-C6 (N7 at 200 permutations) ===        <- results/phase3_pre_c6
    === post-C6, this re-run (all at 1,000 ...     <- results/phase3

Regenerate the committed artefact with:

    /workspace/envs/sasp311/bin/python /workspace/code/m1_n7_prepost.py \
        --out /workspace/results/phase3/m1_n7_prepost.txt

(the table also goes to stdout, so a plain `> file` redirect is equivalent --
that is how the original was made, which is why no writer was ever committed).

RECONSTRUCTION.  This script did not exist when the artefact was written; it
was reconstructed on 2026-08-27 for the reproducibility repair recorded in
`reports/AUDIT_REPRODUCIBILITY.md` ("header and column block appear in no
tracked file").  It was verified byte-exact (`cmp`) against the committed
`results/phase3/m1_n7_prepost.txt` using the two results trees named above,
`results/phase3_pre_c6` for the first block and `results/phase3` for the
second.  It re-reads the frozen summary tables; it re-runs no permutation.

Why it reads the *summary* CSVs rather than `perm_nulls_c1{,_n7}.csv` directly:
`summarize_phase3_c1.py` already reduces each null to the median surviving
fraction over that call's own reportable population, and the committed table is
exactly those medians to three decimals.  Recomputing them here would duplicate
the reportable-population filter and risk drifting from it.
"""
from __future__ import annotations
import sys
import pandas as pd

PRE = "/workspace/results/phase3_pre_c6"
POST = "/workspace/results/phase3"

# Fixed presentation order.  It is NOT the alphabetical order the summary CSVs
# carry: the three Tier A percentile calls are the pre-registered axis and sit
# together, after the two ad-hoc caller definitions.
CALLS = ["cdkn1a_pos", "senepy_p95", "senepy_p99",
         "tierA_p90", "tierA_p95", "tierA_p99"]

# tierA_p95 is the PRIMARY call, so its C1 re-run lives in the primary files
# (`sf_summary_c1.csv` / `perm_nulls_c1.csv`), not in the N7 axis files.  It is
# also the only call that was at 1,000 permutations on both sides of C6, which
# is what makes it the one clean pre/post comparison.
PRIMARY = "tierA_p95"

VAR = ["N3_orig", "N3_tile", "N3_occ", "N3_occ15", "N3_swap", "N3_snap",
       "N4_orig", "N4_tile", "N4_occ", "N4_occ15", "N4_swap"]


def prevalence(res):
    """min/max sender prevalence (%) per call over the six in-band sections."""
    pd_ = pd.read_csv(f"{res}/poisson_density.csv")
    pd_ = pd_[pd_.band == "in_band"]
    g = pd_.groupby("call").sender_prevalence
    return (g.min() * 100).to_dict(), (g.max() * 100).to_dict()


def call_block(res, call):
    """(n_perm, fits, {variant: median SF}) for one sender call in one tree."""
    if call == PRIMARY:
        sf = pd.read_csv(f"{res}/sf_summary_c1.csv")
        perm = pd.read_csv(f"{res}/perm_nulls_c1.csv", usecols=["n_perm"])
    else:
        sf = pd.read_csv(f"{res}/sf_summary_c1_n7.csv")
        sf = sf[sf.call == call]
        perm = pd.read_csv(f"{res}/perm_nulls_c1_n7.csv", usecols=["n_perm"])
    # Both summary files carry the ORIGINAL bounding-box N3/N4 rows beside the
    # C1 re-run rows under the same variant names.  Keep only the C1 re-run:
    # this table is about the corrected in-tissue geometry throughout.
    sf = sf[sf.source.str.contains("c1")]
    med = dict(zip(sf.variant, sf["median"]))
    # `fits` is the reportable population size, which the tile-scope rows do
    # not share (they drop fits that fail the solid-tile receiver floor), so
    # take it from a whole-section row.
    n = int(sf.loc[sf.variant == "N3_occ", "n"].iloc[0])
    return int(perm.n_perm.median()), n, med


def block(res, title):
    lo, hi = prevalence(res)
    out = [title,
           "%10s %6s %5s %8s %8s" % ("call", "n_perm", "fits",
                                     "prev_lo", "prev_hi")
           + "".join("%*s" % (len(v) + 2, v) for v in VAR)]
    for call in CALLS:
        nperm, fits, med = call_block(res, call)
        out.append("%10s %6s %5d %8.2f %8.2f"
                   % (call, "[%d]" % nperm, fits, lo[call], hi[call])
                   + "".join("%*.3f" % (len(v) + 2, med[v]) for v in VAR))
    return out


def main():
    lines = block(PRE, "=== pre-C6 (N7 at 200 permutations) ===")
    lines.append("")
    lines += block(POST,
                   "=== post-C6, this re-run (all at 1,000 permutations) ===")
    txt = "\n".join(lines)
    print(txt)
    if "--out" in sys.argv:
        dest = sys.argv[sys.argv.index("--out") + 1]
        with open(dest, "w") as fh:
            fh.write(txt + "\n")


if __name__ == "__main__":
    main()
