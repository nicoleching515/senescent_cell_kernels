#!/usr/bin/env python3
"""Phase 10 ISO — the ortholog-intersected panel, as ONE importable definition.

PREREG_PHASE8.md §9 item 4 / test A8: every cross-arm quantity is reported twice, once on
each arm's full panel and once on the ortholog-intersected panel.  Phase 9's
`code/h1_a8_crossarm.py` already derives that panel and already asserts its arithmetic; this
module is a thin re-export of THAT derivation so the Phase-10 rescoring cannot drift onto a
second convention.  `h1_a8_crossarm.build()` is imported and called -- the mapping is not
re-derived here.

    human_symbols()  -> frozenset, 2425 HUMAN symbols  (`set(onto.values())`)
    mouse_symbols()  -> frozenset, 2435 MOUSE symbols  (`set(onto)`, the same convention
                        `h1_a8_crossarm.arithmetic` uses for its `m_inter` column)

Pinned counts, asserted on every derivation AND on every cache read:
    mouse panel 5097 | human panel 5093 | mouse genes with an MGI map row 4845 |
    of which on-panel human ortholog 2435 | distinct human symbols 2425

The derivation opens one Xenium .h5 (to read the human panel's feature names), which costs
~20 s and a few hundred MB.  Phase-10 runs it inside joblib workers, so the result is cached
to data/processed_h1/iso_panel_symbols.json (a NEW file; nothing frozen is touched).  The
cache stores the five counts and they are re-asserted on load, so a stale cache cannot pass
silently.  `--rebuild` forces the full derivation.

Usage: python3 code/iso_panel.py [--rebuild]
"""
from __future__ import annotations
import json, os, sys

sys.path.insert(0, "/workspace/code")

CACHE = "/workspace/data/processed_h1/iso_panel_symbols.json"
PINNED = dict(mouse_panel=5097, human_panel=5093, mouse_mapped=4845,
              mouse_onto_panel=2435, human_intersected=2425)

_MEM = None


def _assert_counts(c, where):
    bad = {k: (c.get(k), v) for k, v in PINNED.items() if c.get(k) != v}
    assert not bad, "the pinned panel arithmetic has moved (%s): %s" % (where, bad)
    return True


def build(verbose=True):
    """Re-derive from `h1_a8_crossarm.build()` and refresh the cache."""
    import h1_a8_crossarm as A8
    MPANEL, HPANEL, INTER, onto = A8.build()          # prints + asserts the five counts
    c = dict(mouse_panel=len(MPANEL), human_panel=len(HPANEL), mouse_mapped=None,
             mouse_onto_panel=len(onto), human_intersected=len(INTER))
    # `mapped` is internal to A8.build(); recompute it by the identical rule so the cache
    # carries all five pinned numbers.  (A8.build() has already asserted it.)
    import csv
    ortho = {r["mouse_symbol"]: r["human_symbol"] for r in csv.DictReader(
        open("/workspace/genesets/mouse_human_orthologs_MGI.csv"))}
    c["mouse_mapped"] = len({g for g in MPANEL if g in ortho})
    _assert_counts(c, "derivation")
    d = dict(counts=c, human=sorted(INTER), mouse=sorted(onto),
             map_file="genesets/mouse_human_orthologs_MGI.csv",
             source="code/h1_a8_crossarm.py::build()")
    tmp = CACHE + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(d, fh)
    os.replace(tmp, CACHE)
    if verbose:
        for k, v in PINNED.items():
            print("  ASSERT %-20s %5d == %5d  OK" % (k, c[k], v))
        print("  wrote", CACHE)
    return d


def _load():
    global _MEM
    if _MEM is not None:
        return _MEM
    if os.path.exists(CACHE):
        d = json.load(open(CACHE))
        _assert_counts(d["counts"], "cache " + CACHE)
        assert len(d["human"]) == PINNED["human_intersected"]
        assert len(d["mouse"]) == PINNED["mouse_onto_panel"]
    else:
        d = build(verbose=False)
    _MEM = d
    return d


def human_symbols() -> frozenset:
    return frozenset(_load()["human"])


def mouse_symbols() -> frozenset:
    return frozenset(_load()["mouse"])


if __name__ == "__main__":
    if "--rebuild" in sys.argv or not os.path.exists(CACHE):
        build()
    else:
        d = _load()
        print("cache", CACHE)
        for k, v in PINNED.items():
            print("  ASSERT %-20s %5d == %5d  OK" % (k, d["counts"][k], v))
    print("human_symbols() = %d ; mouse_symbols() = %d"
          % (len(human_symbols()), len(mouse_symbols())))
    print("CDKN1A on intersected human panel : %s" % ("CDKN1A" in human_symbols()))
    print("Cdkn1a on intersected mouse panel : %s" % ("Cdkn1a" in mouse_symbols()))
