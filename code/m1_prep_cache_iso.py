#!/usr/bin/env python3
"""Phase 10 ISO -- build the ortholog-intersected M1 section cache.

Each `data/processed/cache3/<SAMPLE>.npz` is COPIED to
`data/processed_m1_iso/cache3_m1_iso/<SAMPLE>.npz` with ONLY the panel-dependent keys
replaced by the intersected-panel values `code/m1_callers_iso.py` wrote:

    tierA_score
    flag_p90 / flag_p95 / flag_p99        (FINE family -- M1's frozen call)
    mod__<module>   (all 7)
    cdkn1a_counts

There is no `flag_merged_*` key in the mouse cache and none is created: M1's frozen sender
call thresholds within the FINE `cell_type`.  EVERY OTHER KEY is carried through unchanged
and asserted byte-identical after the write, including `senepy_score`, `tierApm__*` and
`flag_pm_*`, which stay at their FULL-PANEL values and are NOT recomputed.  Do not run
`senepy_*` or `tierApm_*` against this cache.

`data/processed/cache3/` is opened read-only and never written.

Usage: python3 code/m1_prep_cache_iso.py [SAMPLE ...]   (default: sasp_phase3.IN_BAND)
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P

SRC = "/workspace/data/processed/cache3"
DST = "/workspace/data/processed_m1_iso/cache3_m1_iso"
CSV = "/workspace/data/processed_m1_iso"

REPLACED = (["tierA_score", "cdkn1a_counts"]
            + ["flag_p%d" % q for q in (90, 95, 99)]
            + ["mod__" + m for m in P.MODULES])
STALE_FULL_PANEL = ["senepy_score"] + ["tierApm__" + m for m in P.MODULES] + \
    ["flag_pm_%s_p%d" % (m, q) for m in P.MODULES for q in (90, 95, 99)]


def build(sample):
    os.makedirs(DST, exist_ok=True)
    src = os.path.join(SRC, sample + ".npz")
    dst = os.path.join(DST, sample + ".npz")
    z = np.load(src, allow_pickle=False)
    old = {k: z[k] for k in z.files}
    z.close()
    assert "flag_merged_p95" not in old, sample + ": mouse cache unexpectedly has merged flags"
    cid = old["cell_id"].astype(str)

    sen = pd.read_csv(CSV + "/senders_iso_%s.csv" % sample).set_index("cell_id")
    sen.index = sen.index.astype(str)
    sen = sen.reindex(pd.Index(cid))
    mods = pd.read_csv(CSV + "/modules_iso_%s.csv" % sample).set_index("cell_id")
    mods.index = mods.index.astype(str)
    mods = mods.reindex(pd.Index(cid))
    assert sen["tierA_score"].notna().all(), sample + ": iso senders csv misses cache cells"
    assert mods.notna().all().all(), sample + ": iso modules csv misses cache cells"
    assert (sen["cell_type"].astype(str).to_numpy() == old["celltype"]).all(), \
        sample + ": fine cell_type disagrees with the cache"
    assert (sen["cell_type_merged"].astype(str).to_numpy() == old["celltype_merged"]).all(), \
        sample + ": merged cell_type disagrees with the cache"

    new = dict(old)
    new["tierA_score"] = sen["tierA_score"].to_numpy(np.float32)
    new["cdkn1a_counts"] = sen["cdkn1a_counts"].to_numpy(np.float32)
    for q in (90, 95, 99):
        new["flag_p%d" % q] = sen["sender_flag_p%d" % q].fillna(0).to_numpy().astype(bool)
    for m in P.MODULES:
        new["mod__" + m] = mods[m].to_numpy(np.float32)

    assert set(new) == set(old), sample + ": the ISO cache changed the key set"
    tmp = dst + ".tmp.npz"
    np.savez_compressed(tmp, **new)
    zz = np.load(tmp, allow_pickle=False)
    n_ident = 0
    for k in old:
        if k in REPLACED:
            assert np.array_equal(zz[k], new[k], equal_nan=True) \
                if zz[k].dtype.kind == "f" else np.array_equal(zz[k], new[k]), \
                "%s: replaced key %s did not round-trip" % (sample, k)
            continue
        same = (np.array_equal(zz[k], old[k], equal_nan=True)
                if zz[k].dtype.kind == "f" else np.array_equal(zz[k], old[k]))
        assert same, "%s: carried-through key %r changed" % (sample, k)
        assert zz[k].dtype == old[k].dtype, "%s: dtype of %r changed" % (sample, k)
        n_ident += 1
    zz.close()
    os.replace(tmp, dst)

    dA = float(np.abs(new["tierA_score"].astype(float)
                      - old["tierA_score"].astype(float)).max())
    return ("[done] %s n=%d  replaced=%d  carried-identical=%d  p95 fine %d->%d  "
            "cdkn1a_delta=%d  max|dtierA|=%.4g"
            % (sample, len(cid), len(REPLACED), n_ident,
               int(old["flag_p95"].sum()), int(new["flag_p95"].sum()),
               int((old["cdkn1a_counts"] != new["cdkn1a_counts"]).sum()), dA))


if __name__ == "__main__":
    print("REPLACED keys      :", REPLACED)
    print("STALE (full-panel) :", len(STALE_FULL_PANEL), "keys, unused by the ISO fits")
    for s in (sys.argv[1:] or list(P.IN_BAND)):
        print(build(s), flush=True)
