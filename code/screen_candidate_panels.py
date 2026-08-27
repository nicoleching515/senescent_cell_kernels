#!/usr/bin/env python3
"""
A3 fallback screen -- Phase 7 section 12.1 step 2, applied to the two runner-up candidates
of reports/PHASE7_H1_SCREEN.md (GSE336890 kidney, GSE335963/GSE335962 bone marrow).

"Verify the panel ON THE DATA, not from the series title." Reads each downloaded
cell_feature_matrix.h5, counts features BOTH by the h5 `feature_type` column AND by the
feature-id prefix, checks the panel is byte-identical across samples (symmetric difference
== 0), and writes one gene_name/gene_id panel CSV per candidate in the SAME format as
genesets/h1_candidate/GSE326743_gene_panel_5093.csv so the frozen disjointness gate can be
run against it unmodified.

WHY BOTH COUNTS: two GSE336890 samples (Region01, Region15) label EVERY feature
'Gene Expression', controls included. The id prefixes (ENSG / NegControlProbe /
NegControlCodeword / Intergenic / UnassignedCodeword / DeprecatedCodeword) are intact in
those files, so the id prefix is the authoritative classifier and `feature_type` is not.
The gene panel is therefore defined as the ENSG-prefixed features.

Nothing is written outside data/raw_h1_candidates/ and results/a3_fallback/. No expression
value is read: only /matrix/features and /matrix/shape.

Run: python3 /workspace/code/screen_candidate_panels.py
"""
import csv, glob, json, os
import h5py

ROOT = '/workspace/data/raw_h1_candidates'
OUT  = '/workspace/results/a3_fallback'
REF  = '/workspace/genesets/h1_candidate/GSE326743_gene_panel_5093.csv'
# stock Prime 5K control complement, Phase 7 section 12.4 / genesets/h1_candidate/PROVENANCE.md
STOCK = {'NegControlCodeword': 609, 'NegControlProbe': 40, 'Intergenic': 21}
CLASS = ['ENSG', 'NegControlProbe', 'NegControlCodeword', 'Intergenic',
         'UnassignedCodeword', 'DeprecatedCodeword']

def features(path):
    with h5py.File(path, 'r') as f:
        g = f['matrix/features']
        rows = list(zip([x.decode() for x in g['name'][:]],
                        [x.decode() for x in g['id'][:]],
                        [x.decode() for x in g['feature_type'][:]]))
        shape = [int(v) for v in f['matrix/shape'][:]]
    return rows, shape

def klass(fid):
    for c in CLASS:
        if fid.startswith(c):
            return c
    return 'OTHER:' + fid[:12]

def main():
    os.makedirs(OUT, exist_ok=True)
    ref = {r['gene_name'] for r in csv.DictReader(open(REF))}
    summary, panels_by_series = {}, {}
    for series in sorted(d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d))):
        h5s = sorted(glob.glob(os.path.join(ROOT, series, '*cell_feature_matrix.h5')))
        if not h5s:
            continue
        print('=' * 100)
        print('%s -- %d cell_feature_matrix.h5' % (series, len(h5s)))
        per, panels, first_rows = {}, {}, None
        for p in h5s:
            rows, shape = features(p)
            s = os.path.basename(p).replace('_cell_feature_matrix.h5', '')
            if first_rows is None:
                first_rows = rows
            by_id, by_type = {}, {}
            for n, i, t in rows:
                by_id[klass(i)] = by_id.get(klass(i), 0) + 1
                by_type[t] = by_type.get(t, 0) + 1
            genes = sorted({n for n, i, t in rows if i.startswith('ENSG')})
            dup = by_id['ENSG'] - len(genes)
            panels[s] = set(genes)
            ft_ok = by_type.get('Gene Expression', 0) == by_id['ENSG']
            print('  %-22s features=%-6d cells=%-8d  gene panel (ENSG ids) = %-5d  dup symbols = %d'
                  % (s, shape[0], shape[1], len(genes), dup))
            print('      by feature id prefix : ' +
                  '  '.join('%s=%d' % (c, by_id.get(c, 0)) for c in CLASS))
            print('      by h5 feature_type   : ' +
                  '  '.join('%s=%d' % (k, v) for k, v in sorted(by_type.items())) +
                  ('' if ft_ok else '   <== feature_type column is UNUSABLE here '
                                    '(controls typed as Gene Expression)'))
            for c, want in STOCK.items():
                got = by_id.get(c, 0)
                print('      stock Prime 5K %-20s expect %3d  got %3d  %s'
                      % (c, want, got, 'MATCH' if got == want else 'MISMATCH'))
            per[s] = dict(file=os.path.basename(p), bytes=os.path.getsize(p),
                          n_features=shape[0], n_cells=shape[1], n_genes=len(genes),
                          n_duplicate_symbols=dup, by_id_prefix=by_id, by_feature_type=by_type,
                          feature_type_column_usable=ft_ok,
                          stock_control_complement_matches=all(
                              by_id.get(c, 0) == w for c, w in STOCK.items()))
        keys = sorted(panels)
        base = panels[keys[0]]
        sd = {k: len(base ^ panels[k]) for k in keys[1:]}
        identical = all(v == 0 for v in sd.values())
        print('  panel identical across %d samples (symmetric difference == 0): %s   %s'
              % (len(keys), identical, sd))
        print('  vs the frozen H1 panel (GSE326743, %d genes): %d shared, %d H1-only, %d candidate-only'
              % (len(ref), len(base & ref), len(ref - base), len(base - ref)))
        rows_csv = sorted({(n, i) for n, i, t in first_rows if i.startswith('ENSG')})
        csvp = os.path.join(OUT, '%s_gene_panel_%d.csv' % (series, len(rows_csv)))
        with open(csvp, 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(['gene_name', 'gene_id']); w.writerows(rows_csv)
        featp = os.path.join(OUT, '%s_panel_features.csv' % series)
        with open(featp, 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(['gene_name', 'gene_id', 'feature_type', 'id_class'])
            w.writerows(sorted((n, i, t, klass(i)) for n, i, t in first_rows))
        print('  wrote %s (%d rows)' % (csvp, len(rows_csv)))
        print('  wrote %s' % featp)
        panels_by_series[series] = base
        summary[series] = dict(per_sample=per, panel_identical_across_samples=identical,
                               symmetric_differences=sd, n_genes=len(base), panel_csv=csvp,
                               source_h5=os.path.basename(h5s[0]),
                               shared_with_H1=len(base & ref), h1_only=len(ref - base),
                               candidate_only=len(base - ref))
    ks = sorted(panels_by_series)
    if len(ks) == 2:
        a, b = panels_by_series[ks[0]], panels_by_series[ks[1]]
        core = a & b & ref
        print('=' * 100)
        print('cross-candidate: |%s ^ %s| = %d' % (ks[0], ks[1], len(a ^ b)))
        print('three-way core (H1 n candidate1 n candidate2) = %d genes' % len(core))
        summary['_cross'] = dict(symmetric_difference=len(a ^ b), three_way_core=len(core))
    json.dump(summary, open(os.path.join(OUT, 'panel_screen.json'), 'w'), indent=1)
    print('\nwrote %s/panel_screen.json' % OUT)

if __name__ == '__main__':
    main()
