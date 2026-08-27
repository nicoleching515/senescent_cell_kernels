"""Phase 10 ISO -- bind the FROZEN Phase-3 estimator to the H1 arm on the
ORTHOLOG-INTERSECTED panel, and send its output to `results/phase10_h1_iso/`.

This is `code/h1_phase10.py` with two lines changed: the cache and the results directory.
Nothing else differs -- same estimator, same section list, same label family, same
`tierAmg_pNN` call shim (PREREG_PHASE8.md §13 / `CS_PHASE9_H1_AUDIT.md` §9.4).  The
`h1_phase10.py` binding is NOT imported, because importing it would install a `sender_mask`
shim closed over the full-panel cache path; the shim is reproduced here verbatim instead.

Import this BEFORE `run_phase3_nulls`, in the parent AND inside every joblib worker (loky
workers are fresh interpreters):

    import h1_phase10_iso           # noqa: F401  (side effect)
    import run_phase3_nulls as RN

WARNING.  `data/processed_h1/cache3_h1_iso/` carries FULL-PANEL `senepy_score`,
`tierApm__*` and `flag_pm_*` keys (see `code/h1_prep_cache_iso.py`).  Only `tierA_pNN` and
`tierAmg_pNN` are valid ISO calls; `senepy_*` and `tierApm_*` would silently mix panels.
"""
import os, sys
sys.path.insert(0, "/workspace/code")
import numpy as np
import h1_sec                      # noqa: F401  (Phase 9 arm binding: sections, labels)
import sasp_phase3 as P
import h1_common as H

CACHE_ISO = H.PROC + "/cache3_h1_iso"
RESULTS10_ISO = "/workspace/results/phase10_h1_iso"
os.makedirs(RESULTS10_ISO, exist_ok=True)
P.CACHE3 = CACHE_ISO
P.RESULTS = RESULTS10_ISO

VALID_ISO_CALLS = ("tierA_p90", "tierA_p95", "tierA_p99",
                   "tierAmg_p90", "tierAmg_p95", "tierAmg_p99")

# --- the `tierAmg_pNN` call shim, copied from `code/h1_phase10.py` -----------
MERGED_CALLS = ("tierAmg_p90", "tierAmg_p95", "tierAmg_p99")
MERGED_ALIAS = {"tierA_merged_p90": "tierAmg_p90",
                "tierA_merged_p95": "tierAmg_p95",
                "tierA_merged_p99": "tierAmg_p99"}

if not getattr(P.Sec, "_h1_merged_call_shim", False):
    _orig_sender_mask = P.Sec.sender_mask

    def sender_mask(self, call, module=None):
        call = MERGED_ALIAS.get(call, call)
        if call.startswith("tierAmg_p"):
            key = "flag_merged_p" + call[9:]
            if key not in self.z.files:
                raise ValueError("%s: cache has no %r" % (self.name, key))
            ok = ~np.isin(self.celltype,
                          P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
            return self.z[key].copy() & ok
        return _orig_sender_mask(self, call, module=module)

    P.Sec.sender_mask = sender_mask
    P.Sec._h1_merged_call_shim = True

assert P.CACHE3 == CACHE_ISO, "the ISO cache binding was overwritten"
