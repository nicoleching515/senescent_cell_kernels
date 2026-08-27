#!/usr/bin/env python3
"""
A3 fallback screen -- donor count and per-donor cell budget for the two runner-up candidates.

Master Plan Test 3 / audit test A3 requires "at least one tissue-threshold combination in 1-20%,
with >= 200 senders and >= 5,000 non-senders PER DONOR", reported per cell type. A3 cannot be
measured before the arm is annotated, but the CELL BUDGET that upper-bounds it can be, and it is
structural (cell counts only, no expression):

  * a donor with fewer than 5,000 cells of a type can never supply 5,000 non-senders of it;
  * at prevalence p, >= 200 senders of a type needs >= 200/p cells of that type in that donor
    (2,000 cells at p = 10 %, 4,000 at 5 %, 20,000 at 1 %).

Sources, all files on disk:
  GSE335963 (bone marrow): /matrix/shape of each cell_feature_matrix.h5; donor identity from the
      GSE335962 sample titles (CH02, CH15, NC03, NC03-2, NC05, NC05-2 -- the "-2" sections are
      REPEAT SECTIONS of the same donor).
  GSE336890 (kidney): each Xenium Region slide carries several patients, so the per-cell donor
      assignment is the depositors' GSE336890_cells_stats_dir.tar.gz; one CSV per
      (patient x kidney region) block. Block sums are cross-checked against the /matrix/shape
      of the Region h5 files that were downloaded.

Run: python3 /workspace/code/screen_candidate_cellbudget.py
"""
import csv, glob, json, os, re
from collections import defaultdict
import h5py

ROOT = '/workspace/data/raw_h1_candidates'
OUT  = '/workspace/results/a3_fallback'
THRESH = [0.20, 0.10, 0.05, 0.02, 0.01]

def h5_cells(p):
    with h5py.File(p, 'r') as f:
        return int(f['matrix/shape'][1])

def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []

    # ---------------- GSE335963 / GSE335962, bone marrow ----------------
    print('=' * 100)
    print('GSE335963 (Xenium subseries GSE335962) -- bone marrow biopsy')
    bm = []
    for p in sorted(glob.glob(ROOT + '/GSE335963/*cell_feature_matrix.h5')):
        s = os.path.basename(p).replace('_cell_feature_matrix.h5', '')
        section = s.split('_', 1)[1]
        donor = section.split('-')[0]
        n = h5_cells(p)
        bm.append(dict(series='GSE335963', section=section, donor=donor, n_cells=n))
        print('  %-10s donor=%-6s cells=%d' % (section, donor, n))
    per_donor = defaultdict(int)
    for r in bm:
        per_donor[r['donor']] += r['n_cells']
    print('  sections=%d  DONORS=%d  total cells=%d'
          % (len(bm), len(per_donor), sum(r['n_cells'] for r in bm)))
    for d in sorted(per_donor):
        print('    donor %-6s cells=%7d' % (d, per_donor[d]))
    for d, n in per_donor.items():
        rows.append(dict(series='GSE335963', tissue='bone marrow', donor=d, n_blocks=
                         sum(1 for r in bm if r['donor'] == d), n_cells=n))

    # ---------------- GSE336890, kidney ----------------
    print('=' * 100)
    print('GSE336890 -- kidney biopsy, per-cell donor assignment from cells_stats')
    blocks = []
    for p in sorted(glob.glob(ROOT + '/cells_stats/cells_stats_dir/Region_*/*_cells_stats.csv')):
        region = p.split('/')[-2]
        base = os.path.basename(p).replace('_cells_stats.csv', '')
        m = re.match(r'^(AIN|ATI|REF)_(C|M)_(\d+)_(\S+)$', base)
        dis, reg, idx, sid = m.groups()
        n = sum(1 for _ in open(p)) - 1
        blocks.append(dict(region=region, block=base, disease=dis, kidney_region=reg,
                           patient_idx=idx, specimen=sid, n_cells=n))
    print('  blocks=%d  total cells=%d' % (len(blocks), sum(b['n_cells'] for b in blocks)))
    # cross-check against the downloaded Region h5 files
    print('  cross-check, block sum vs /matrix/shape of the Region h5:')
    for p in sorted(glob.glob(ROOT + '/GSE336890/*cell_feature_matrix.h5')):
        rn = re.search(r'_Region(\d+)_', os.path.basename(p)).group(1)
        key = 'Region_%d' % int(rn)
        bs = sum(b['n_cells'] for b in blocks if b['region'] == key)
        h = h5_cells(p)
        print('    %-10s blocks=%7d  h5=%7d  diff=%+d %s'
              % (key, bs, h, bs - h, '' if bs == h else '  <-- cells with no donor assignment'))
    kd = defaultdict(lambda: dict(n_cells=0, n_blocks=0, regions=set(), disease=None))
    for b in blocks:
        k = '%s/%s' % (b['disease'], b['specimen'])
        kd[k]['n_cells'] += b['n_cells']; kd[k]['n_blocks'] += 1
        kd[k]['regions'].add(b['kidney_region']); kd[k]['disease'] = b['disease']
    print('  DONORS (distinct patient specimens) = %d' % len(kd))
    for k in sorted(kd, key=lambda k: -kd[k]['n_cells']):
        v = kd[k]
        print('    %-14s disease=%-4s blocks=%d regions=%-3s cells=%7d'
              % (k, v['disease'], v['n_blocks'], ','.join(sorted(v['regions'])), v['n_cells']))
        rows.append(dict(series='GSE336890', tissue='kidney', donor=k, n_blocks=v['n_blocks'],
                         n_cells=v['n_cells']))
    with open(OUT + '/GSE336890_donor_blocks.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0])); w.writeheader(); w.writerows(blocks)

    # ---------------- the A3 cell-budget arithmetic ----------------
    print('=' * 100)
    print('A3 cell budget. A donor needs >= 5,000 cells OF A TYPE for the non-sender floor, and')
    print('>= 200/p cells of that type for the >= 200-sender floor at prevalence p.')
    print('%-12s %-8s %8s %8s %8s   %s' % ('series', 'donors', 'min', 'median', 'max',
                                           'donors with >= 5,000 cells IN TOTAL'))
    budget = {}
    for series in ['GSE336890', 'GSE335963']:
        ns = sorted(r['n_cells'] for r in rows if r['series'] == series)
        med = ns[len(ns) // 2] if len(ns) % 2 else (ns[len(ns) // 2 - 1] + ns[len(ns) // 2]) / 2
        ge5k = sum(1 for n in ns if n >= 5000)
        print('%-12s %-8d %8d %8.0f %8d   %d/%d' % (series, len(ns), ns[0], med, ns[-1], ge5k, len(ns)))
        need = {('%.0f%%' % (100 * p)): int(-(-200 // p)) for p in THRESH}
        budget[series] = dict(n_donors=len(ns), per_donor_cells=ns, min=ns[0], median=med,
                              max=ns[-1], donors_ge_5000_total_cells=ge5k,
                              cells_needed_for_200_senders_at_prevalence=need)
        print('    cells of a type needed for >= 200 senders: ' +
              '  '.join('p=%s -> %d' % (k, v) for k, v in need.items()))
    # ---- what fraction of a donor must ONE cell type be, for that donor to clear both floors?
    # f_required(p) = max(5000/(1-p), 200/p) / n_cells_in_donor.  Measured cell counts only;
    # nothing is assumed about the composition.
    print()
    print('Required share of a donor held by a SINGLE cell type for that (donor x type) stratum')
    print('to clear both A3 floors -- max(5000/(1-p), 200/p) / donor cells:')
    print('%-12s %-16s %s' % ('series', 'donor cells', '  '.join('p=%.0f%%' % (100 * p) for p in THRESH)))
    frac = {}
    for series in ['GSE336890', 'GSE335963']:
        ns = sorted((r['n_cells'] for r in rows if r['series'] == series))
        _med = ns[len(ns) // 2] if len(ns) % 2 else (ns[len(ns) // 2 - 1] + ns[len(ns) // 2]) / 2
        for lab, n in (('min donor', ns[0]), ('median donor', _med), ('max donor', ns[-1])):
            f = [max(5000.0 / (1 - p), 200.0 / p) / n for p in THRESH]
            print('%-12s %-16s %s' % (series, '%s %.0f' % (lab, n),
                                      '  '.join('%6.1f%%' % (100 * x) for x in f)))
            frac['%s/%s' % (series, lab)] = dict(n_cells=n,
                                                 required_share={'p=%.2f' % p: round(x, 4)
                                                                 for p, x in zip(THRESH, f)})
    print('  (a share > 100 %% means NO cell type in that donor can satisfy A3 at that prevalence)')
    json.dump(frac, open(OUT + '/required_type_share.json', 'w'), indent=1)

    with open(OUT + '/cell_budget.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['series', 'tissue', 'donor', 'n_blocks', 'n_cells'])
        w.writeheader(); w.writerows(rows)
    json.dump(budget, open(OUT + '/cell_budget.json', 'w'), indent=1)
    print('\nwrote %s/cell_budget.csv and cell_budget.json' % OUT)

if __name__ == '__main__':
    main()
