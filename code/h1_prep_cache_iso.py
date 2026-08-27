#!/usr/bin/env python3
"""Phase 10 ISO -- build the ortholog-intersected H1 section cache.

Each `data/processed_h1/cache3_h1/<SEC>.npz` is COPIED to
`data/processed_h1/cache3_h1_iso/<SEC>.npz` with ONLY the panel-dependent keys replaced by
the intersected-panel values that `code/h1_callers_iso.py` wrote:

    tierA_score
    flag_p90 / flag_p95 / flag_p99
    flag_merged_p90 / flag_merged_p95 / flag_merged_p99
    mod__<module>       (all 7 of sasp_phase3.MODULES)
    cdkn1a_counts

EVERY OTHER KEY IS CARRIED THROUGH UNCHANGED AND ASSERTED BYTE-IDENTICAL after the write:
coords, cell_id, celltype, celltype_merged, densities, knn_idx, nn1_um, median_nn_um,
median_depth, ct_conf, transcript_counts, genes_detected, cell_area, nucleus_area,
n_nucleus_area_imputed, seg_levels, seg_code, zonation_score, compartment,
dist_to_boundary_um, dist_to_portal_triad_um -- and also `senepy_score`, `tierApm__*` and
`flag_pm_*`, which are the FULL-PANEL values and are NOT recomputed (see the declared
omissions in `code/h1_callers_iso.py`).  DO NOT run `senepy_*` or `tierApm_*` calls against
this cache: those keys are full-panel and would silently mix panels.

`data/processed_h1/cache3_h1/` is opened read-only and never written.

Usage: python3 code/h1_prep_cache_iso.py [SECTION ...]
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_common as H

SRC = H.PROC + "/cache3_h1"
DST = H.PROC + "/cache3_h1_iso"

REPLACED = (["tierA_score", "cdkn1a_counts"]
            + ["flag_p%d" % q for q in (90, 95, 99)]
            + ["flag_merged_p%d" % q for q in (90, 95, 99)]
            + ["mod__" + m for m in H.MODULES])

# keys that stay at their FULL-PANEL values, declared so the mixing is auditable
STALE_FULL_PANEL = ["senepy_score"] + ["tierApm__" + m for m in H.MODULES] + \
    ["flag_pm_%s_p%d" % (m, q) for m in H.MODULES for q in (90, 95, 99)]


def build(section):
    os.makedirs(DST, exist_ok=True)
    src = os.path.join(SRC, section + ".npz")
    dst = os.path.join(DST, section + ".npz")
    z = np.load(src, allow_pickle=False)
    old = {k: z[k] for k in z.files}
    z.close()
    cid = old["cell_id"].astype(str)

    sen = pd.read_csv(H.PROC + "/senders_h1_iso_%s.csv" % section).set_index("cell_id")
    sen.index = sen.index.astype(str)
    sen = sen.reindex(pd.Index(cid))
    mods = pd.read_csv(H.PROC + "/modules_h1_iso_%s.csv" % section).set_index("cell_id")
    mods.index = mods.index.astype(str)
    mods = mods.reindex(pd.Index(cid))
    assert sen["tierA_score"].notna().all(), section + ": iso senders csv misses cache cells"
    assert mods.notna().all().all(), section + ": iso modules csv misses cache cells"
    # the cell set and the fine labels must be identical to the full-panel cache
    assert (sen["cell_type"].astype(str).to_numpy() == old["celltype"]).all(), \
        section + ": fine cell_type disagrees with the cache"
    assert (sen["cell_type_merged"].astype(str).to_numpy() == old["celltype_merged"]).all(), \
        section + ": merged cell_type disagrees with the cache"

    new = dict(old)
    new["tierA_score"] = sen["tierA_score"].to_numpy(np.float32)
    new["cdkn1a_counts"] = sen["cdkn1a_counts"].to_numpy(np.float32)
    for q in (90, 95, 99):
        new["flag_p%d" % q] = sen["sender_flag_p%d" % q].fillna(0).to_numpy().astype(bool)
        new["flag_merged_p%d" % q] = (sen["sender_flag_merged_p%d" % q]
                                      .fillna(0).to_numpy().astype(bool))
    for m in H.MODULES:
        new["mod__" + m] = mods[m].to_numpy(np.float32)

    assert set(new) == set(old), section + ": the ISO cache changed the key set"
    tmp = dst + ".tmp.npz"
    np.savez_compressed(tmp, **new)
    zz = np.load(tmp, allow_pickle=False)
    n_ident = 0
    for k in old:                       # every non-replaced key byte-identical
        if k in REPLACED:
            assert np.array_equal(zz[k], new[k], equal_nan=True) \
                if zz[k].dtype.kind == "f" else np.array_equal(zz[k], new[k]), \
                "%s: replaced key %s did not round-trip" % (section, k)
            continue
        same = (np.array_equal(zz[k], old[k], equal_nan=True)
                if zz[k].dtype.kind == "f" else np.array_equal(zz[k], old[k]))
        assert same, "%s: carried-through key %r changed" % (section, k)
        assert zz[k].dtype == old[k].dtype, "%s: dtype of %r changed" % (section, k)
        n_ident += 1
    zz.close()
    os.replace(tmp, dst)

    dA = float(np.abs(new["tierA_score"].astype(float)
                      - old["tierA_score"].astype(float)).max())
    return ("[done] %s n=%d  replaced=%d  carried-identical=%d  "
            "p95 fine %d->%d  merged %d->%d  cdkn1a_delta=%d  max|dtierA|=%.4g"
            % (section, len(cid), len(REPLACED), n_ident,
               int(old["flag_p95"].sum()), int(new["flag_p95"].sum()),
               int(old["flag_merged_p95"].sum()), int(new["flag_merged_p95"].sum()),
               int((old["cdkn1a_counts"] != new["cdkn1a_counts"]).sum()), dA))


if __name__ == "__main__":
    print("REPLACED keys      :", REPLACED)
    print("STALE (full-panel) :", len(STALE_FULL_PANEL), "keys, unused by the ISO fits")
    for s in (sys.argv[1:] or list(H.ALL_SECTIONS)):
        print(build(s), flush=True)
