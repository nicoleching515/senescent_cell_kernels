"""Phase 10 — bind the FROZEN Phase-3 estimator to the H1 arm and send its output to
`results/phase10_h1/`.

Import this BEFORE `run_phase3_nulls` / `run_phase3_var` / `run_phase8_compmatch`, in the
parent AND inside every joblib worker (loky workers are fresh interpreters).

    import h1_phase10                # noqa: F401  (side effect)
    import run_phase3_nulls as RN

It does three things and nothing else:

 1. `import h1_sec` — Phase 9's binding: `sasp_phase3.CACHE3` -> the H1 cache,
    IN_BAND/ALL_SECTIONS -> the 7 spleen sections, CANON_TYPES -> the spleen label family.
 2. repoints `sasp_phase3.RESULTS` at `results/phase10_h1/` so no Phase-9 or mouse output
    can be overwritten.  (`h1_sec` points it at `results/phase9_h1/`.)
 3. adds ONE additional sender-call name, `tierAmg_pNN`, to `sasp_phase3.Sec.sender_mask`
    by wrapping it — the frozen file is not edited.  `tierAmg_pNN` is the identical Tier A
    percentile rule evaluated at the MERGED label family instead of the fine one; it is
    deviation H5's declared sensitivity (`CS_PHASE9_H1_AUDIT.md` §9.4, `PREREG_PHASE8.md`
    §13), it is NEVER the primary, and no threshold is changed.  The flags it reads were
    written by Phase 9's `code/h1_callers.py` and carried into the cache by
    `code/h1_cache_extend.py`.
"""
import os, sys
sys.path.insert(0, "/workspace/code")
import numpy as np
import h1_sec                      # noqa: F401  (Phase 9 arm binding)
import sasp_phase3 as P

RESULTS10 = "/workspace/results/phase10_h1"
os.makedirs(RESULTS10, exist_ok=True)
P.RESULTS = RESULTS10

MERGED_CALLS = ("tierAmg_p90", "tierAmg_p95", "tierAmg_p99")

if not getattr(P.Sec, "_h1_merged_call_shim", False):
    _orig_sender_mask = P.Sec.sender_mask

    def sender_mask(self, call, module=None):
        if call.startswith("tierAmg_p"):
            key = "flag_merged_p" + call[9:]
            if key not in self.z.files:
                raise ValueError("%s: cache has no %r; run code/h1_cache_extend.py"
                                 % (self.name, key))
            ok = ~np.isin(self.celltype,
                          P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
            return self.z[key].copy() & ok
        return _orig_sender_mask(self, call, module=module)

    P.Sec.sender_mask = sender_mask
    P.Sec._h1_merged_call_shim = True
