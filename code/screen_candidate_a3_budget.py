#!/usr/bin/env python3
"""
A3 fallback screen -- the per-(donor x cell type) cell budget that upper-bounds audit test A3.

A3 (Master Plan Test 3): "at least one tissue-threshold combination in 1-20 %, with >= 200 senders
and >= 5,000 non-senders per donor", reported PER CELL TYPE. The prevalence itself cannot be
measured before the arm is annotated and scored, but the number of (donor x cell type) strata that
could EVER satisfy it is fixed by cell counts alone:

    non-sender floor : a (donor, type) stratum needs >= 5,000 cells of that type;
    sender floor     : at prevalence p it needs >= 200/p cells of that type
                       (1,000 at p = 20 %; 2,000 at 10 %; 4,000 at 5 %; 20,000 at 1 %).
    Both floors at once, at prevalence p: n_cells >= max(5000 / (1 - p), 200 / p).

This runs on GSE335963 (bone marrow) only, because it is the only candidate whose depositors
ship per-cell type labels. GSE336890 (kidney) ships no cell-type annotation at all -- its
cells_stats files carry CellID, Transcripts and Area and nothing else -- so for kidney only the
per-donor total is knowable before Job B, and that is reported as such.

Input: the depositors' *_Metadata.csv.gz (Annotation column). No expression value is read.
Run: python3 /workspace/code/screen_candidate_a3_budget.py
"""
import csv, glob, gzip, json, os
from collections import defaultdict

ROOT = '/workspace/data/raw_h1_candidates/GSE335963'
OUT  = '/workspace/results/a3_fallback'
P    = [0.20, 0.10, 0.05, 0.02, 0.01]

def need(p):
    return max(5000.0 / (1 - p), 200.0 / p)

def main():
    os.makedirs(OUT, exist_ok=True)
    per = defaultdict(lambda: defaultdict(int))     # donor -> type -> n
    sect = defaultdict(lambda: defaultdict(int))    # section -> type -> n
    cover = {}
    for f in sorted(glob.glob(ROOT + '/*_Metadata.csv.gz')):
        rows = list(csv.DictReader(gzip.open(f, 'rt')))
        s = os.path.basename(f).split('_', 1)[1].replace('_Metadata.csv.gz', '')
        donor = s.split('-')[0]
        cover[s] = len(rows)
        for r in rows:
            per[donor][r['Annotation']] += 1
            sect[s][r['Annotation']] += 1
    types = sorted({t for d in per.values() for t in d})
    print('GSE335963 / GSE335962 bone marrow -- depositor Annotation labels: %d types' % len(types))
    print('annotated cells per section: %s  (total %d)'
          % (cover, sum(cover.values())))
    print('\n### cells per (donor x cell type), 4 donors ###')
    donors = sorted(per)
    print('%-22s %s' % ('cell type', ' '.join('%9s' % d for d in donors)))
    for t in sorted(types, key=lambda t: -sum(per[d][t] for d in donors)):
        print('%-22s %s' % (t, ' '.join('%9d' % per[d][t] for d in donors)))
    print('%-22s %s' % ('TOTAL', ' '.join('%9d' % sum(per[d].values()) for d in donors)))

    print('\n### how many (donor x type) strata clear the A3 floors ###')
    print('%-8s %-10s %-14s %s' % ('p', 'need n>=', 'strata (of %d)' % (len(donors) * len(types)),
                                   'types clearing it in ALL 4 donors'))
    res = {}
    for p in P:
        n = need(p)
        ok = [(d, t) for d in donors for t in types if per[d][t] >= n]
        allfour = [t for t in types if all(per[d][t] >= n for d in donors)]
        print('%-8s %-10.0f %-14d %s' % ('%.0f%%' % (100 * p), n, len(ok),
                                         ', '.join(allfour) or '(none)'))
        res['p=%.2f' % p] = dict(cells_needed=round(n, 1), n_strata_ok=len(ok),
                                 strata_ok=['%s/%s' % x for x in ok],
                                 types_ok_in_all_donors=allfour)
    print('\n### the same, per SECTION (the pre-registration reports A3 per section) ###')
    secs = sorted(sect)
    for p in P:
        n = need(p)
        ok = [(s, t) for s in secs for t in types if sect[s][t] >= n]
        print('  p=%-5s need n>=%-7.0f strata clearing: %d of %d'
              % ('%.0f%%' % (100 * p), n, len(ok), len(secs) * len(types)))
        res['p=%.2f' % p]['n_section_strata_ok'] = len(ok)
        res['p=%.2f' % p]['n_section_strata_total'] = len(secs) * len(types)
    with open(OUT + '/GSE335963_donor_x_celltype.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['donor', 'cell_type', 'n_cells'])
        for d in donors:
            for t in types:
                w.writerow([d, t, per[d][t]])
    json.dump(dict(annotated_cells_per_section=cover, n_types=len(types), types=types,
                   per_donor_x_type={d: dict(per[d]) for d in donors}, a3_floors=res),
              open(OUT + '/GSE335963_a3_budget.json', 'w'), indent=1)
    print('\nwrote %s/GSE335963_donor_x_celltype.csv and GSE335963_a3_budget.json' % OUT)

if __name__ == '__main__':
    main()
