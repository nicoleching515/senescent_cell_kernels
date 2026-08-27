#!/usr/bin/env python3
"""
A3 fallback screen -- run the FROZEN Phase 7 section 11 / audit test A2 disjointness gate
against a fallback candidate's real, measured panel.

This does NOT reimplement the gate. It reads code/gate_disjointness_human.py verbatim and
executes it with exactly two constants rebound:

    PANEL  <- the candidate's gene panel CSV (produced by code/screen_candidate_panels.py
              from the candidate's own cell_feature_matrix.h5)
    OUT    <- results/a3_fallback/gate_<SERIES>/

Every assertion, every Tier A / Tier B definition, every module, the per-module sender
sensitivity and the CoreScence circularity check are the frozen ones from
genesets/human/. Nothing under results/phase7_jobA/ is touched.

Usage: python3 code/gate_disjointness_candidate.py <SERIES> <panel_csv>
Exit status is the frozen gate's own: 0 = FROZEN CONFIGURATION PASS, 1 = FAIL.
"""
import os, re, sys

SRC = '/workspace/code/gate_disjointness_human.py'

def main(series, panel_csv):
    out = '/workspace/results/a3_fallback/gate_%s' % series
    os.makedirs(out, exist_ok=True)
    src = open(SRC).read()
    src, n1 = re.subn(r"^PANEL = \{r\['gene_name'\].*$",
                      "PANEL = {r['gene_name'] for r in csv.DictReader(open(%r))}" % panel_csv,
                      src, count=1, flags=re.M)
    src, n2 = re.subn(r"^OUT\s*=\s*.*$", "OUT   = %r" % out, src, count=1, flags=re.M)
    assert n1 == 1 and n2 == 1, 'substitution failed (%d, %d) -- the frozen gate changed shape' % (n1, n2)
    print('#' * 100)
    print('# A3 FALLBACK SCREEN: frozen A2 gate, candidate panel')
    print('#   candidate series : %s' % series)
    print('#   panel            : %s' % panel_csv)
    print('#   gate source      : %s (verbatim, PANEL and OUT rebound)' % SRC)
    print('#   outputs          : %s' % out)
    print('#   NOTE: the banner below is printed by the frozen script and says "H1 (GSE326743)".')
    print('#         The panel actually loaded is the candidate panel named above.')
    print('#' * 100)
    g = {'__name__': '__main__', '__file__': SRC}
    try:
        exec(compile(src, SRC + ' [panel=%s]' % series, 'exec'), g)
    except SystemExit as e:
        sys.exit(e.code)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
