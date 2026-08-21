"""Restore the NumPy <2.0 scalar-type aliases that NumPy 2 removed.

Why this exists
---------------
`commot` 0.0.3 -- the published package, which Phase 4 runs as released software
rather than reimplementing -- has `rho = np.Inf` as a *default argument* in
`commot/_optimal_transport/_usot.py`, so it raises at **import** time under the
numpy version this project pins (2.4.6):

    AttributeError: `np.Inf` was removed in the NumPy 2.0 release.

These are pure aliases (`np.Inf is np.inf` on 1.x), so restoring them changes no
numerical behaviour -- it only makes the released package importable.  Verified
by re-running a stored Phase 4 COMMOT job and diffing against the checkpoint
written under numpy 1.26.3; see CS_PHASE4 section 8.

Import this module *before* importing commot.
"""
import numpy as np

_ALIASES = {
    "Inf": "inf", "Infinity": "inf", "infty": "inf", "NINF": None,
    "NaN": "nan", "NAN": "nan",
    "float_": "float64", "complex_": "complex128",
    "unicode_": "str_", "string_": "bytes_",
    "bool8": "bool_", "int0": "intp", "uint0": "uintp",
}

for _name, _target in _ALIASES.items():
    if hasattr(np, _name):
        continue
    if _target is None:
        if _name == "NINF":
            np.NINF = -np.inf
        continue
    try:
        setattr(np, _name, getattr(np, _target))
    except AttributeError:
        pass
