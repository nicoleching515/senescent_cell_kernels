#!/usr/bin/env python3
"""
Phase 7 section 14 step 2 -- does SenePy actually cover SPLEEN? (Job A follow-on, task 5)

SenePy 1.0.1 as installed, inspected on disk. NO H1 DATA IS READ: the only inputs are the
package's own bundled hub pickles and the GSE326743 panel membership file.

The mouse arm scored SenePy with TISSUE-MATCHED mouse Liver hubs (phase2_downstream.py HUBMAP)
and required >= 10 hub genes on panel. This script asks the same two questions for the human arm:
  1. Does SenePy ship any SPLEEN hub?
  2. For each spleen cell type in markers_human_spleen.py, what is the best available surrogate
     hub, from which tissue, and how many of its genes are on the human panel?

Run: python3 /workspace/code/senepy_coverage_human.py
"""
import csv, json, os, pickle, re, sys
sys.path.insert(0, '/workspace/code')
import human_symbols as HS
from markers_human_spleen import MARKERS

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _pkgdata   # resolves DeepScence/senepy package data; see AUDIT_REPRODUCIBILITY D2
D = _pkgdata.senepy_data(verbose=False)
OUT = '/workspace/results/phase7_jobA'
MIN_ON_PANEL = 10                      # the mouse arm's own threshold (phase2_downstream.py)

hubs = pickle.load(open(D + '6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle', 'rb'))
meta = pickle.load(open(D + '6_HUMAN_HUBS_METADATA_FILTERED.pickle', 'rb'))

cov = {}
for (t, c, h), genes in hubs.items():
    onp = [g for g, _ in genes if HS.resolve(g)]
    key = (t, c)
    if key not in cov or len(onp) > cov[key][1]:
        cov[key] = (h, len(onp), len(genes))

tissues = sorted({t for t, _ in cov})
print('=' * 100)
print('SenePy 1.0.1 human hub inventory (bundled with the installed package)')
print('=' * 100)
print('hubs: %d   tissues: %d   (tissue, cell) pairs: %d' % (len(hubs), len(tissues), len(cov)))
print('tissues: %s' % ', '.join(tissues))
print('\n*** SPLEEN HUBS: %s ***' % ('present' if 'spleen' in tissues else 'NONE. SenePy ships no spleen signature.'))

# every hub SenePy has for an immune / stromal cell, with panel support
print('\nHubs for immune or stromal cell types, best hub per (tissue, cell):')
print('%-14s %-28s %5s %7s %9s' % ('tissue', 'cell', 'hub', 'size', 'on-panel'))
for (t, c), (h, o, n) in sorted(cov.items(), key=lambda kv: (-kv[1][1])):
    print('%-14s %-28s %5d %7d %9d' % (t, c, h, n, o))

# ---- surrogate mapping for the spleen label set -----------------------------------------
# Word-boundary matches on the SenePy `cell` string only; no judgement beyond the match.
# (Plain substring matching is wrong here: 'b cell' matches 'club cell'.)
SURROGATE = {
 'Red pulp macrophages':        ['macrophage'],
 'Monocytes':                   ['monocyte'],
 'cDC1':                        [], 'cDC2': [], 'pDC': [],
 'Follicular B cells':          ['b cell'],
 'Marginal zone B cells':       ['b cell'],
 'Germinal centre B cells':     ['b cell'],
 'Plasma cells':                ['plasma cell'],
 'CD4 T cells':                 ['t cell', 'lymphoid cell'],
 'CD8 T cells':                 ['t cell', 'lymphoid cell'],
 'NK cells':                    ['nk cell'],
 'Follicular dendritic cells':  [],
 'Fibroblastic reticular cells':['fibroblast'],
 'Sinusoidal endothelium':      ['endothelial'],
 'Endothelial cells':           ['endothelial'],
 'Lymphatic endothelium':       [],
 'Smooth muscle / capsule':     ['smooth muscle'],
 'Pericytes':                   ['pericyte'],
 'Fibroblasts':                 ['fibroblast'],
 'Erythroid cells':             [],
 'Megakaryocytes':              [],
 'Neutrophils':                 ['neutrophil'],
 'Mesothelial cells':           [],
 'Proliferating cells':         ['mitotic cell'],
}

rows, n_ok, n_none = [], 0, 0
print('\n### Spleen label set vs SenePy human hubs ###')
print('%-30s %-38s %9s %s' % ('spleen cell type', 'best surrogate hub (tissue / cell)', 'on-panel', 'usable?'))
for lab in list(MARKERS) + [l for l in SURROGATE if l not in MARKERS]:
    pats = SURROGATE.get(lab, [])
    best = None
    for (t, c), (h, o, n) in cov.items():
        if any(re.search(r'\b' + re.escape(p) + r'\b', c) for p in pats):
            if best is None or o > best[3]:
                best = (t, c, h, o, n)
    in_labelset = lab in MARKERS
    if best is None:
        n_none += in_labelset
        print('%-30s %-38s %9s %s' % (lab, 'NONE -- SenePy has no hub for this type', '-',
                                      'NO' + ('' if in_labelset else '  (label dropped anyway)')))
        rows.append(dict(cell_type=lab, in_label_set=in_labelset, hub_tissue='', hub_cell='',
                         hub_num='', hub_size='', on_panel=0, usable='no_hub',
                         tissue_matched='no'))
    else:
        t, c, h, o, n = best
        ok = o >= MIN_ON_PANEL
        n_ok += in_labelset and ok
        print('%-30s %-38s %9d %s' % (lab, '%s / %s (hub %d)' % (t, c, h), o,
                                      ('yes' if ok else 'NO, <%d on-panel' % MIN_ON_PANEL)
                                      + ('' if in_labelset else '  (label dropped anyway)')))
        rows.append(dict(cell_type=lab, in_label_set=in_labelset, hub_tissue=t, hub_cell=c,
                         hub_num=h, hub_size=n, on_panel=o,
                         usable='yes' if ok else 'below_min_on_panel', tissue_matched='no'))

print('\n### Verdict ###')
print('  SenePy ships NO spleen hub. Every hub above is a CROSS-TISSUE SURROGATE.')
print('  The mouse arm used TISSUE-MATCHED Liver hubs; the human arm cannot.')
print('  Of the %d assignable spleen labels: %d have a usable surrogate hub (>= %d on-panel genes),'
      % (len(MARKERS), n_ok, MIN_ON_PANEL))
print('  %d have no SenePy hub of any tissue.' % n_none)

with open(OUT + '/senepy_spleen_coverage.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['cell_type', 'in_label_set', 'hub_tissue', 'hub_cell',
                                       'hub_num', 'hub_size', 'on_panel', 'usable', 'tissue_matched'])
    w.writeheader()
    for r in rows:
        w.writerow(r)
json.dump(dict(senepy_tissues=tissues, spleen_hub=False, n_hubs=len(hubs),
               n_labels=len(MARKERS), n_usable_surrogate=n_ok, n_no_hub=n_none,
               min_on_panel=MIN_ON_PANEL),
          open(OUT + '/senepy_spleen_coverage.json', 'w'), indent=1)
print('\nWritten to %s/senepy_spleen_coverage.{csv,json}' % OUT)
