"""Phase 9 — point the FROZEN Phase-3 estimator at the H1 cache.

Importing this module rebinds three module-level constants in `sasp_phase3` so that
`sasp_phase3.Sec`, `phase3_core.*` and `run_phase3_nulls.SectionFit` read the H1 section
cache instead of the mouse one.  No estimator code is copied, edited or reimplemented, and
nothing under data/processed/ or results/phase3/ is read or written by anything that imports
this.  Import it BEFORE constructing any Sec.

    import h1_sec          # noqa: F401  (side effect)
    import sasp_phase3 as P
    sec = P.Sec("SPLN07")
"""
import sys
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
import h1_common as H

P.CACHE3 = H.PROC + "/cache3_h1"
P.RESULTS = H.RESULTS
# All 7 H1 sections are analysed (PREREG §3.11: H1 has no analogue of the Test-3
# admissibility rule and none is invented).  Setting IN_BAND makes Sec.meta['band']
# read 'in_band' for every H1 section rather than the mouse default 'below_floor'.
P.IN_BAND = list(H.ALL_SECTIONS)
P.OVER_CEILING = []
P.BELOW_FLOOR = []
P.ALL_SECTIONS = list(H.ALL_SECTIONS)

# Merged-label families of the frozen spleen marker set (markers_human_spleen.MERGE),
# plus the standalone fine labels.  Used where the mouse arm uses CANON_TYPES_MERGED to
# harmonise the k-NN composition columns across pooled sections.
CANON_TYPES_MERGED_H1 = ("B cells", "CD4 T cells", "CD8 T cells", "Endothelial",
                         "Erythroid cells", "Low_quality", "Megakaryocytes",
                         "Mesothelial cells", "Mono/Mac/DC", "NK cells", "Neutrophils",
                         "Plasma cells", "Stromal")
P.CANON_TYPES_MERGED = CANON_TYPES_MERGED_H1
P.CANON_TYPES = CANON_TYPES_MERGED_H1
