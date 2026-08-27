#!/usr/bin/env python3
"""Phase 10 — carry the ALREADY-COMPUTED per-module and merged-label sender flags
into the H1 section cache.

Nothing is recomputed from expression and no threshold is touched.  Phase 9's
`code/h1_callers.py` already scored the seven `genesets/human/A_sender_for_<module>.txt`
sets and already wrote the merged-label Tier A flags into
`data/processed_h1/senders_h1_<sec>.csv`; `code/h1_prep_cache.py` simply did not copy
them into the .npz, so `sasp_phase3.Sec.sender_mask("tierApm_pNN")` raises and the
pre-registered D1 sensitivity axis (PREREG §3.7 a / decision D1) cannot run on H1.

What this adds to each existing cache file, leaving every existing key byte-identical:

  tierApm__<m>            the per-module Tier A score      (senders csv, verbatim)
  flag_pm_<m>_p95         the committed `tierApm_flag_p95__<m>` column, VERBATIM.
  flag_pm_<m>_p{90,99}    the same rule (`h1_callers.pct_flags`, within FINE cell type,
                          >= 20 cells, strict >) re-applied at 90 and 99 to the score as
                          STORED, which the senders csv rounds to 5 dp.  DECLARED: at p95,
                          where the committed unrounded flag exists to compare against, the
                          rounded recomputation differs on 0-12 cells of 200k-390k
                          (<= 0.006 %), all exact ties at the percentile boundary.  The
                          committed column is used at p95 so the primary axis is unaffected;
                          the p90/p99 flags carry that <= 0.006 % tie noise.  Reported, not
                          patched.  The per-section maximum is written to
                          results/phase10_h1/cache_extend_rounding.csv.
  flag_merged_p{90,95,99} the frozen Tier A percentile rule applied at the MERGED
                          label family (senders csv, verbatim).  This is deviation H5's
                          declared sensitivity (`CS_PHASE9_H1_AUDIT.md` §9.4); it is a
                          SENSITIVITY, never the primary, and no threshold is tuned.

Usage: python3 code/h1_cache_extend.py [--check] [SECTION ...]
"""
import sys, os, argparse
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H

H.RESULTS10 = "/workspace/results/phase10_h1"
os.makedirs(H.RESULTS10, exist_ok=True)

CACHE = H.PROC + "/cache3_h1"
EXCLUDE_FROM_STRATA = {"Low_quality", "Unknown"}      # h1_callers.py, verbatim


def pct_flags(score, labels, q, min_cells=20):
    """`code/h1_callers.py::pct_flags`, verbatim."""
    f = np.zeros(len(score), int)
    for c in pd.unique(labels):
        if c in EXCLUDE_FROM_STRATA:
            continue
        m = (labels == c).to_numpy()
        if m.sum() < min_cells:
            continue
        f[m] = (score[m] > np.nanpercentile(score[m], q)).astype(int)
    return f


def extend(section, check_only=False):
    p = os.path.join(CACHE, section + ".npz")
    z = np.load(p, allow_pickle=False)
    old = {k: z[k] for k in z.files}
    cid = old["cell_id"].astype(str)
    sen = pd.read_csv(H.PROC + "/senders_h1_%s.csv" % section).set_index("cell_id")
    sen.index = sen.index.astype(str)
    sen = sen.reindex(pd.Index(cid))
    assert sen["tierA_score"].notna().all(), section + ": senders csv misses cache cells"
    ct = sen["cell_type"].astype(str)

    new, rnd = {}, []
    for m in H.MODULES:
        s = sen["tierApm_score__" + m].to_numpy(np.float32)
        new["tierApm__" + m] = s
        ref = sen["tierApm_flag_p95__" + m].to_numpy().astype(bool)
        new["flag_pm_%s_p95" % m] = ref                    # committed column, verbatim
        d95 = int((pct_flags(s.astype(float), ct, 95).astype(bool) != ref).sum())
        rnd.append(dict(section=section, module=m, n_cells=len(cid),
                        n_p95_flag_diff_rounded_vs_committed=d95,
                        frac=d95 / len(cid)))
        for q in (90, 99):
            new["flag_pm_%s_p%d" % (m, q)] = pct_flags(s.astype(float), ct, q).astype(bool)
    for q in (90, 95, 99):
        new["flag_merged_p%d" % q] = (sen["sender_flag_merged_p%d" % q]
                                      .fillna(0).to_numpy().astype(bool))
    # sanity: the frozen fine-label Tier A rule reproduces the committed cache flag
    # to within the same rounding tie noise (the cache flag is the authoritative one and
    # is NOT rewritten).
    for q in (90, 95, 99):
        f = pct_flags(sen["tierA_score"].to_numpy(float), ct, q).astype(bool)
        d = int((f != old["flag_p%d" % q]).sum())
        rnd.append(dict(section=section, module="tierA_p%d" % q, n_cells=len(cid),
                        n_p95_flag_diff_rounded_vs_committed=d, frac=d / len(cid)))
        assert d / len(cid) < 1e-4, "%s: tierA_p%d rule drifted (%d cells)" % (section, q, d)

    msg = ("%s n=%d  +%d keys  senders p95 fine=%d merged=%d"
           % (section, len(cid), len(new), int(old["flag_p95"].sum()),
              int(new["flag_merged_p95"].sum())))
    pd.DataFrame(rnd).to_csv(H.RESULTS10 + "/cache_extend_rounding_%s.csv" % section,
                             index=False)
    if check_only:
        return "[check] " + msg + "  max_tie_diff=%d" % max(r["n_p95_flag_diff_rounded_vs_committed"] for r in rnd)
    d = dict(old); d.update(new)
    tmp = p + ".tmp.npz"
    np.savez_compressed(tmp, **d)
    zz = np.load(tmp, allow_pickle=False)
    for k in old:                       # every pre-existing key byte-identical
        assert np.array_equal(zz[k], old[k], equal_nan=True) if zz[k].dtype.kind == "f" \
            else np.array_equal(zz[k], old[k]), "%s: key %s changed" % (section, k)
    zz.close(); z.close()
    os.replace(tmp, p)
    return "[done] " + msg + "  max_tie_diff=%d" % max(r["n_p95_flag_diff_rounded_vs_committed"] for r in rnd)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("sections", nargs="*", default=None)
    a = ap.parse_args()
    for s in (a.sections or list(H.ALL_SECTIONS)):
        print(extend(s, a.check), flush=True)
