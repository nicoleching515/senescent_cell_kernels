"""Phase 10 ISO -- bind the FROZEN Phase-3 estimator to the M1 arm on the
ORTHOLOG-INTERSECTED panel, and send its output to `results/phase10_m1_iso/`.

The mouse analogue of `code/h1_phase10_iso.py`, built on the same pattern as
`code/h1_sec.py` + `code/h1_phase10.py`: rebind `sasp_phase3.CACHE3` and
`sasp_phase3.RESULTS`, leave every estimator untouched.  M1 keeps its FROZEN literals --
the fine/merged label sets, `IN_BAND`, `EXCLUDE_TYPES`, `LABELS="merged"` -- so nothing
here changes except which .npz the `Sec` accessor opens and where output lands.  No call
shim is added: M1's call is the frozen `tierA_pNN`.

Import this BEFORE `run_phase3_nulls`, in the parent AND inside every joblib worker.

WARNING.  `data/processed_m1_iso/cache3_m1_iso/` carries FULL-PANEL `senepy_score`,
`tierApm__*` and `flag_pm_*` keys (see `code/m1_prep_cache_iso.py`).  Only `tierA_pNN` is
a valid ISO call on this arm.
"""
import os, sys
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P

CACHE_ISO = "/workspace/data/processed_m1_iso/cache3_m1_iso"
RESULTS10_ISO = "/workspace/results/phase10_m1_iso"
os.makedirs(RESULTS10_ISO, exist_ok=True)
P.CACHE3 = CACHE_ISO
P.RESULTS = RESULTS10_ISO

SECTIONS = list(P.IN_BAND)          # the six in-band sections, frozen literal
VALID_ISO_CALLS = ("tierA_p90", "tierA_p95", "tierA_p99")

assert P.CACHE3 == CACHE_ISO, "the ISO cache binding was overwritten"
