#!/usr/bin/env python3
"""
Phase 7 test A6 -- the arm-specific anatomical covariate for the HUMAN SPLEEN arm (Job A
follow-on, tasks 2 and 4). Also evaluates the Tier E2 candidates against it.

A6 as written in the Phase 7 doc specifies a LUNG airway-to-alveolar axis. H1 is spleen
(reports/PHASE7_H1_SCREEN.md), so that instruction is void. The spleen analogue of liver zonation
is the RED PULP / WHITE PULP axis.

THIS SCRIPT BUILDS THE GENE-SIDE OF THE SPECIFICATION ONLY. It cannot build or validate the
covariate: the mouse zonation axis was derived from the GSM9295284 EXPRESSION MATRIX, and H1
expression is behind the section 15 freeze. Nothing here reads H1 data beyond panel membership.
The distance measures are specified in reports/BIO_PHASE7_JobA_FOLLOWON.md section 2 and are
implementable from cells.parquet coordinates plus the Job B cell-type calls, with no further
biological judgement.

Compartment marker sets come from markers_human_spleen.py, i.e. from the pinned CellMarker 2.0
table, and already carry the Tier A + Tier C over-adjustment guard.

Run: python3 /workspace/code/spec_a6_compartments_human.py
"""
import csv, glob, json, os, sys
from collections import OrderedDict
sys.path.insert(0, '/workspace/code')
import human_symbols as HS
from markers_human_spleen import MARKERS

HUM, OUT = '/workspace/genesets/human', '/workspace/results/phase7_jobA'

def gl(n):
    return {l.strip() for l in open(os.path.join(HUM, n)) if l.strip()}

TIER_A = set().union(*[gl(os.path.basename(f)) for f in glob.glob(HUM + '/A_*.txt')])
TIER_B = set().union(*[gl('B_%s.txt' % m) for m in
    ('tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
     'emt_ecm', 'oxidative_stress', 'secondary_senescence')])
TIER_C = gl('C_ligands.txt') | gl('C_receptors.txt')

# Compartment = union of the marker sets of the cell types that DEFINE that compartment.
COMPARTMENTS = OrderedDict([
 ('D_spleen_white_pulp_tzone',   ['CD4 T cells', 'CD8 T cells', 'Fibroblastic reticular cells']),
 ('D_spleen_white_pulp_follicle',['Follicular B cells', 'Germinal centre B cells']),
 ('D_spleen_marginal_zone',      ['Marginal zone B cells']),
 ('D_spleen_red_pulp',           ['Red pulp macrophages', 'Sinusoidal endothelium',
                                  'Erythroid cells']),
 ('D_spleen_capsule_trabecula',  ['Smooth muscle / capsule', 'Pericytes']),
])

print('=' * 100)
print('A6 SPLEEN COMPARTMENT MARKER SETS (gene side of the specification)')
print('=' * 100)
comp = OrderedDict()
for name, labs in COMPARTMENTS.items():
    missing = [l for l in labs if l not in MARKERS]
    g = sorted(set().union(*[set(MARKERS[l]) for l in labs if l in MARKERS]))
    comp[name] = g
    print('%-32s %2d genes from %-58s %s'
          % (name, len(g), ', '.join(l for l in labs if l in MARKERS),
             ('MISSING LABEL: ' + ','.join(missing)) if missing else ''))
    print('     %s' % ' '.join(g))

print('\n### Compartment x compartment overlap (must be small or the axis is not separable) ###')
hdr = 'compartment'.ljust(32) + ''.join(('C%d' % i).ljust(6) for i in range(1, len(comp) + 1))
print(hdr)
for i, (n1, g1) in enumerate(comp.items(), 1):
    print(('C%d %s' % (i, n1)).ljust(32) + ''.join(str(len(set(g1) & set(g2))).ljust(6) for g2 in comp.values()))

print('\n### Over-adjustment audit: compartment genes that are also Tier A / B / C ###')
for n, g in comp.items():
    print('  %-32s A:%d  B:%d  C:%d   %s'
          % (n, len(set(g) & TIER_A), len(set(g) & TIER_B), len(set(g) & TIER_C),
             ' '.join(sorted(set(g) & (TIER_A | TIER_B | TIER_C))) or '-'))
print('  (Tier A and Tier C are 0 by construction -- the guard is applied in the marker build.')
print('   Tier B members are FLAGGED, not removed: see build_markers_human_spleen.py.)')

for n, g in comp.items():
    with open(os.path.join(HUM, n + '.txt'), 'w') as fh:
        fh.write('\n'.join(g) + '\n')

# ---------------------------------------------------------------- Tier E2 candidates
print('\n' + '=' * 100)
print('TIER E2 CANDIDATES -- "a cell-type identity program unrelated to inflammation"')
print('=' * 100)
print('The mouse arm\'s E2 (hepatocyte identity) FAILED: its score correlated with the zonation')
print('axis at r = -0.67, worse than 99%% of random size-matched sets (genesets/README.md sec 6).')
print('The decisive test for any candidate here is the same correlation against the A6 axis, and')
print('that test NEEDS EXPRESSION -- it is post-freeze. What can be decided now is whether a')
print('candidate is confounded with the A6 axis BY CONSTRUCTION, i.e. whether its own genes are')
print('the compartment definition.\n')
CANDIDATES = OrderedDict([
 ('E2a_erythroid',            ['Erythroid cells']),
 ('E2b_smooth_muscle_capsule',['Smooth muscle / capsule']),
 ('E2c_pericyte',             ['Pericytes']),
 ('E2d_structural_union',     ['Smooth muscle / capsule', 'Pericytes']),
 ('E2e_sinusoidal_endothelium',['Sinusoidal endothelium']),
])
rows = []
print('%-30s %5s %4s %4s %4s  %s' % ('candidate', 'n', 'A', 'B', 'C', 'in which A6 compartment'))
for name, labs in CANDIDATES.items():
    g = sorted(set().union(*[set(MARKERS[l]) for l in labs if l in MARKERS]))
    inc = [n for n, cg in comp.items() if set(g) & set(cg)]
    rows.append(dict(candidate=name, n=len(g), n_tierA=len(set(g) & TIER_A),
                     n_tierB=len(set(g) & TIER_B), n_tierC=len(set(g) & TIER_C),
                     a6_compartments=';'.join(inc), genes=' '.join(g)))
    print('%-30s %5d %4d %4d %4d  %s' % (name, len(g), len(set(g) & TIER_A), len(set(g) & TIER_B),
                                         len(set(g) & TIER_C), ', '.join(inc) or 'NONE'))
with open(OUT + '/tierE2_candidates.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['candidate', 'n', 'n_tierA', 'n_tierB', 'n_tierC',
                                       'a6_compartments', 'genes'])
    w.writeheader()
    for r in rows:
        w.writerow(r)
json.dump(dict(compartments={k: v for k, v in comp.items()}, e2_candidates=rows),
          open(OUT + '/a6_compartments_and_E2.json', 'w'), indent=1)
print('\nWritten: %d D_spleen_*.txt in %s ; %s/tierE2_candidates.csv' % (len(comp), HUM, OUT))
