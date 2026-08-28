#!/usr/bin/env python3
"""
A3 fallback screen -- does SenePy ship TISSUE-MATCHED hubs for the fallback candidates?

Deviation P4 of reports/PREREG_PHASE8.md records that SenePy 1.0.1 ships 65 human hubs across
10 tissues and NO spleen hub, so on H1 SenePy is not the same estimator it is on M1 (which used
tissue-matched mouse Liver hubs). This script asks the same question for GSE336890 (kidney) and
GSE335963/GSE335962 (bone marrow), reusing code/senepy_coverage_human.py's method: inspect the
package's own bundled hub pickles, count hub genes surviving code/human_symbols.resolve()
against the candidate's measured panel, and apply the mouse arm's own MIN_ON_PANEL = 10.

NO candidate expression data is read: inputs are the senepy pickles and the panel CSVs written
by code/screen_candidate_panels.py.

Run: python3 /workspace/code/senepy_coverage_candidates.py
"""
import csv, json, os, pickle, sys
sys.path.insert(0, '/workspace/code')
import human_symbols as HS

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _pkgdata   # resolves DeepScence/senepy package data; see AUDIT_REPRODUCIBILITY D2
D   = _pkgdata.senepy_data(verbose=False)
OUT = '/workspace/results/a3_fallback'
MIN_ON_PANEL = 10
PANELS = {
    'H1_GSE326743_spleen': ('spleen', '/workspace/genesets/h1_candidate/GSE326743_gene_panel_5093.csv'),
    'GSE336890_kidney':    ('kidney', OUT + '/GSE336890_gene_panel_5101.csv'),
    'GSE335963_bonemarrow':('bone marrow', OUT + '/GSE335963_gene_panel_5101.csv'),
}

def rebind(panel_csv):
    rows = list(csv.DictReader(open(panel_csv)))
    HS.PANEL = {r['gene_name'] for r in rows}
    HS.ENSG2SYM = {r['gene_id']: r['gene_name'] for r in rows}
    return len(HS.PANEL)

hubs = pickle.load(open(D + '6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle', 'rb'))
tissues = sorted({t for t, _, _ in hubs})
print('SenePy human hub inventory: %d hubs, %d tissues' % (len(hubs), len(tissues)))
print('tissues: %s' % ', '.join(tissues))

summary = {}
for tag, (tissue, panel_csv) in PANELS.items():
    n = rebind(panel_csv)
    present = tissue in tissues
    print('\n' + '=' * 100)
    print('%s   panel = %s (%d genes)' % (tag, os.path.basename(panel_csv), n))
    print('  tissue-matched SenePy hub for %r: %s' % (tissue, 'PRESENT' if present else 'ABSENT'))
    rows = []
    for (t, c, h), genes in sorted(hubs.items()):
        if t != tissue:
            continue
        onp = [g for g, _ in genes if HS.resolve(g)]
        rows.append(dict(tissue=t, cell=c, hub=h, hub_size=len(genes), on_panel=len(onp),
                         usable=('yes' if len(onp) >= MIN_ON_PANEL else 'no')))
    if rows:
        print('  %-34s %5s %8s %9s %8s' % ('cell type (SenePy label)', 'hub', 'size', 'on-panel', 'usable'))
        for r in sorted(rows, key=lambda r: -r['on_panel']):
            print('  %-34s %5s %8d %9d %8s' % (r['cell'], r['hub'], r['hub_size'], r['on_panel'], r['usable']))
        print('  %d hubs, %d distinct cell types, %d usable at MIN_ON_PANEL=%d'
              % (len(rows), len({r['cell'] for r in rows}),
                 sum(r['usable'] == 'yes' for r in rows), MIN_ON_PANEL))
    else:
        print('  no hubs at all for this tissue -- every cell type needs a cross-tissue surrogate')
    summary[tag] = dict(tissue=tissue, panel=panel_csv, n_panel=n, tissue_hub_present=present,
                        n_hubs=len(rows), n_cell_types=len({r['cell'] for r in rows}),
                        n_usable=sum(r['usable'] == 'yes' for r in rows), hubs=rows)
    if rows:
        with open(OUT + '/senepy_coverage_%s.csv' % tag, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

summary['_inventory'] = dict(n_hubs=len(hubs), tissues=tissues)
json.dump(summary, open(OUT + '/senepy_coverage_candidates.json', 'w'), indent=1)
print('\nwrote %s/senepy_coverage_candidates.json' % OUT)
