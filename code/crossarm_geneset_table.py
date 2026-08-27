#!/usr/bin/env python3
"""
Phase 8 freeze -- the two arms side by side, before and after C6. Feeds Phase 7 section 17 and
reports/PREREG_PHASE8_genesets.md.

Every number is recomputed here from files on disk: the two panels, the two pinned MSigDB
archives, the frozen human sets (genesets/human/), the frozen mouse sets (genesets/*.txt, the
PRE-C6 configuration the existing Phase 2-5 results were computed under) and the mouse C6 sets
(genesets/mouse_c6/).

It also quantifies where the two arms CANNOT be made identical, because a documented asymmetry is
acceptable and an undocumented one is not.

Run: python3 /workspace/code/crossarm_geneset_table.py
"""
import csv, gzip, json, os, sys
from collections import OrderedDict
sys.path.insert(0, '/workspace/code')
import human_symbols as HS

W = '/workspace'
RES = W + '/results/phase7_jobA'
MODS = ['tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
        'emt_ecm', 'oxidative_stress', 'secondary_senescence']

# ---------------- panels
mp = {r['gene_name'] for r in csv.DictReader(open(W + '/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv'))}
p100 = [r['Gene'] for r in csv.DictReader(gzip.open(W + '/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz', 'rt'))]
GENO = {g for g in p100 if ('_WT' in g or '_ALT' in g or '_del_' in g or '_splice_' in g)}
MPANEL = (mp | set(p100)) - GENO
HPANEL = HS.PANEL
ORTHO = {r['mouse_symbol']: r['human_symbol'] for r in
         csv.DictReader(open(W + '/genesets/mouse_human_orthologs_MGI.csv'))}

def rd(path):
    return {l.strip() for l in open(path) if l.strip()}

def m(n, c6=False):
    return rd('%s/genesets/%s%s.txt' % (W, 'mouse_c6/' if c6 else '', n))


# Phase 8 / 8.7: `genesets/mouse_c6/` was PROMOTED into `genesets/`, so
# `genesets/*.txt` is now the C6 state and can no longer serve as the "pre-C6"
# column.  The pre-C6 mouse sets are read from the tag that captures them.
import subprocess as _sp


def m_pre(n):
    r = _sp.run(['git', '-C', W, 'show', 'pre-c6-genesets:genesets/%s.txt' % n],
                capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit('cannot read pre-C6 %s from tag pre-c6-genesets: %s'
                         % (n, r.stderr.strip()))
    return {l.strip() for l in r.stdout.splitlines() if l.strip()}

def h(n, var=False):
    return rd('%s/genesets/human/%s%s.txt' % (W, 'variants/' if var else '', n))

def gap_split(mouse_set, human_set):
    """Split the arm-exclusive members of a pair of sets into MGI-MAP GAPS and REAL
    asymmetries, and return the map-gap-corrected shared count.

    A gap is a gene the two arms BOTH carry where the pinned 1:1 MGI report simply has no
    row: the mouse symbol is absent from the map and its upper-case form is a member of the
    human set (equivalently, the Title-case form of a human-only symbol is a member of the
    mouse set).  The project records this failure mode in `genesets/human/_symbol_resolutions.csv`
    (e.g. `Cdkn2b,CDKN2B,NOT in MGI map; upper-case form is a known HGNC symbol`) and already
    corrected for it on Tier C.  It is applied here to EVERY set, so the arms are described
    the same way everywhere -- see reports/AUDIT_PHASE8_FACTCHECK.md R2.
    """
    m2h = {ORTHO[g] for g in mouse_set if g in ORTHO}
    unmapped = sorted(g for g in mouse_set if g not in ORTHO)
    ho = sorted(human_set - m2h)
    ho_gap = [g for g in ho if g.capitalize() in mouse_set]
    ho_real = [g for g in ho if g.capitalize() not in mouse_set]
    mo_mapped = sorted(m2h - human_set)
    mo_gap = [g for g in unmapped if g.upper() in human_set]
    mo_real = [g for g in unmapped if g.upper() not in human_set]
    shared = len(m2h & human_set)
    shared_corr = shared + len(ho_gap)
    d = dict(n_mouse=len(mouse_set), n_human=len(human_set), n_mouse_mapped=len(m2h),
             n_overlap=shared, n_overlap_map_gap_corrected=shared_corr,
             human_only=' '.join(ho), human_only_map_gap=' '.join(ho_gap),
             human_only_real=' '.join(ho_real),
             mouse_only=' '.join(mo_mapped),
             mouse_unmapped=' '.join(unmapped),
             mouse_only_map_gap=' '.join(mo_gap), mouse_only_real=' '.join(mo_real),
             mouse_only_complete=' '.join(sorted(mo_mapped) + sorted(mo_real)))
    # arithmetic the previous version of this table failed: every non-shared member of each
    # arm must be accounted for by name.
    assert len(ho_real) == len(human_set) - shared_corr, (ho_real, shared_corr)
    assert len(mo_mapped) + len(mo_real) == len(mouse_set) - shared_corr, d
    assert len(mo_gap) == len(ho_gap), d
    return d



print('=' * 108)
print('CROSS-ARM GENE SET TABLE -- M1 (GSE310392, mouse liver) vs H1 (GSE326743, human spleen)')
print('=' * 108)
print('Panels: mouse %d genes (5,106 Gene Expression features minus %d genotyping probes) ; '
      'human %d genes' % (len(MPANEL), len(GENO), len(HPANEL)))

# ---------------- ortholog-intersected panel (Phase 7 test A8)
m2h_on = {g: ORTHO[g] for g in MPANEL if g in ORTHO}
XPANEL_H = {v for v in m2h_on.values() if v in HPANEL}
XPANEL_M = {k for k, v in m2h_on.items() if v in HPANEL}
# The pinned MGI report is one row per mouse gene but is NOT injective: 18,782 rows carry
# 17,609 distinct human symbols and 431 human symbols receive more than one mouse gene. So
# "mouse panel genes that map onto the human panel" (2,435) and "distinct human symbols they
# map to" (2,425) are DIFFERENT counts; the intersected panel is the human-symbol count.
# See reports/AUDIT_PHASE8_FACTCHECK.md M2.
from collections import Counter as _C
_n_h_distinct = len(set(ORTHO.values()))
_n_h_multi = sum(1 for _v in _C(ORTHO.values()).values() if _v > 1)
print('Ortholog map (MGI, pinned): %d rows, one per mouse gene, onto %d distinct human symbols '
      '(%d human symbols receive more than one mouse gene -- the map is many-to-one, not 1:1).'
      % (len(ORTHO), _n_h_distinct, _n_h_multi))
print('  %d of %d mouse panel genes have a row in the map; %d of those map onto the human panel, '
      'and they land on %d DISTINCT human symbols.'
      % (len(m2h_on), len(MPANEL), len(XPANEL_M), len(XPANEL_H)))
print('ORTHOLOG-INTERSECTED PANEL (test A8): %d genes (human symbols).' % len(XPANEL_H))
print('  Mouse-panel genes with NO human ortholog in the map: %d' % (len(MPANEL) - len(m2h_on)))
print('  Human-panel genes not reachable from the mouse panel: %d' % (len(HPANEL) - len(XPANEL_H)))

rows = []
def row(tier, name, mouse_pre, mouse_c6, human_pre, human_c6, note=''):
    r = OrderedDict(tier=tier, set=name,
                    mouse_preC6='' if mouse_pre is None else len(mouse_pre),
                    mouse_C6='' if mouse_c6 is None else len(mouse_c6),
                    human_preC6='' if human_pre is None else len(human_pre),
                    human_C6='' if human_c6 is None else len(human_c6))
    for lab, s in (('mouse_C6', mouse_c6), ('human_C6', human_c6)):
        pass
    r['orthologue_intersected_mouse_C6'] = ('' if mouse_c6 is None
                                            else len(set(mouse_c6) & XPANEL_M))
    r['orthologue_intersected_human_C6'] = ('' if human_c6 is None
                                            else len(set(human_c6) & XPANEL_H))
    r['note'] = note
    rows.append(r)
    return r

# Tier A. Human pre-C6 sizes come from the recorded gate run, recomputed here where files exist.
h_pre_strict = None       # superseded; recompute: A0_h minus union(B1..B6, B7v1)
Hsum = json.load(open(W + '/genesets/human/_test2_summary.json'))
A0_h = set(Hsum['A_ported_pre']) & HPANEL
B_h_c6 = {k: h('B_' + k) for k in MODS}
B_h_pre = dict(B_h_c6)
B_h_pre['secondary_senescence'] = h('B_secondary_senescence_v1_curated_ported', var=True)
h_pre_strict = A0_h - set().union(*B_h_pre.values())
h_c6_strict = h('A_SENDER_FINAL_strict')

Msum = json.load(open(RES + '/mouse_c6_summary.json'))
A0_m = set(Msum['A0'])
m_pre_strict = m_pre('A_SENDER_FINAL_strict')       # from tag pre-c6-genesets, n=25
m_c6_strict = m('A_SENDER_FINAL_strict', c6=True)

row('A', 'A0 candidate pool (pre-disjointness)', A0_m, A0_m, A0_h, A0_h,
    'unchanged by C6 -- C6 only changes B7')
row('A', 'A_SENDER_FINAL_strict (PRIMARY)', m_pre_strict, m_c6_strict, h_pre_strict, h_c6_strict,
    'C6 raises both arms to 33')
for k in MODS:
    row('A', 'A_sender_for_%s (sensitivity)' % k,
        None, m('A_sender_for_' + k, c6=True), None, h('A_sender_for_' + k), '')
for k in MODS:
    note = 'RE-SOURCED by C6' if k == 'secondary_senescence' else 'unchanged by C6'
    row('B', 'B_' + k, m_pre('B_' + k), m('B_' + k, c6=True),
        B_h_pre[k], B_h_c6[k], note)
row('C', 'C_ligands', m('C_ligands'), m('C_ligands'), h('C_ligands'), h('C_ligands'),
    'CXCL8/CXCR1 human-only; CXCL1/MMP3/TIMP1 mouse-only')
row('C', 'C_receptors', m('C_receptors'), m('C_receptors'), h('C_receptors'), h('C_receptors'), '')
row('E', 'E_housekeeping_planA', m('E_housekeeping_planA'), m('E_housekeeping_planA'),
    h('E_housekeeping_planA'), h('E_housekeeping_planA'), 'fails in BOTH arms')
row('E', 'E_housekeeping_expanded', m('E_housekeeping_expanded'), m('E_housekeeping_expanded'),
    h('E_housekeeping_expanded'), h('E_housekeeping_expanded'), 'human is thinner (8 vs 13)')

print('\n### Tier A and Tier B, mouse vs human, before and after C6 ###')
hdr = ('%-46s %10s %8s %10s %8s %10s %10s' %
       ('set', 'mouse_pre', 'mouse_C6', 'human_pre', 'human_C6', 'orth_M', 'orth_H'))
print(hdr); print('-' * len(hdr))
for r in rows:
    print('%-46s %10s %8s %10s %8s %10s %10s'
          % (r['set'], r['mouse_preC6'], r['mouse_C6'], r['human_preC6'], r['human_C6'],
             r['orthologue_intersected_mouse_C6'], r['orthologue_intersected_human_C6']))

with open(RES + '/crossarm_geneset_table.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- asymmetry audit
print('\n' + '=' * 108)
print('ASYMMETRY AUDIT -- where the arms cannot be made identical')
print('=' * 108)

def msig(species, name):
    d = json.load(open('%s/genesets/msigdb_%s/%s.json'
                       % (W, 'mouse_2026.1.Mm' if species == 'm' else 'human_2026.1.Hs', name)))
    k = list(d)[0]
    return set(d[k]['geneSymbols']), d[k]['systematicName']

asym = []
for name in ('SAUL_SEN_MAYO', 'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP'):
    gm, idm = msig('m', name)
    gh, idh = msig('h', name)
    gm2h = {ORTHO[g] for g in gm if g in ORTHO}
    inter = gm2h & gh
    print('\n%s' % name)
    print('  mouse %s n=%3d   human %s n=%3d' % (idm, len(gm), idh, len(gh)))
    print('  mouse mapped to human symbols via the pinned MGI map: %d of %d (%d have no 1:1 ortholog)'
          % (len(gm2h), len(gm), len(gm) - len(gm2h)))
    print('  overlap of the mapped mouse set with the human set: %d' % len(inter))
    print('  human-only members: %d   mouse-only (after mapping): %d'
          % (len(gh - gm2h), len(gm2h - gh)))
    print('  -> THE TWO SPECIES VERSIONS ARE NOT TRANSLATIONS OF EACH OTHER.')
    asym.append(dict(set=name, mouse_id=idm, human_id=idh, n_mouse=len(gm), n_human=len(gh),
                     n_mouse_mapped=len(gm2h), n_overlap=len(inter),
                     n_human_only=len(gh - gm2h), n_mouse_only=len(gm2h - gh)))

b7m = m('B_secondary_senescence', c6=True)
b7h = h('B_secondary_senescence')
g7 = gap_split(b7m, b7h)
print('\nFROZEN B7, arm to arm')
print('  mouse %d   human %d' % (g7['n_mouse'], g7['n_human']))
print('  mouse B7 mapped to human symbols: %d ; overlap with human B7: %d by the pinned map, '
      '%d after the map-gap correction (gaps: %s) ; on the ortholog-intersected panel: '
      'mouse %d / human %d'
      % (g7['n_mouse_mapped'], g7['n_overlap'], g7['n_overlap_map_gap_corrected'],
         g7['mouse_only_map_gap'] or '-', len(b7m & XPANEL_M), len(b7h & XPANEL_H)))
asym.append(dict(set='B7_FROZEN',
                 n_mouse_on_xpanel=len(b7m & XPANEL_M), n_human_on_xpanel=len(b7h & XPANEL_H),
                 **g7))

am, ah = m('A_SENDER_FINAL_strict', c6=True), h('A_SENDER_FINAL_strict')
ga = gap_split(am, ah)
print('\nFROZEN A_SENDER_FINAL_strict, arm to arm')
print('  mouse %d   human %d   (identical size, NOT identical membership)'
      % (ga['n_mouse'], ga['n_human']))
print('  mouse mapped to human symbols: %d ; overlap with the human set: %d by the pinned map'
      % (ga['n_mouse_mapped'], ga['n_overlap']))
print('  MGI-MAP GAPS (on BOTH arms, the pinned map lacks the row): %s'
      % (ga['mouse_only_map_gap'] or '-'))
print('  -> SHARED, map-gap corrected: %d of %d' % (ga['n_overlap_map_gap_corrected'],
                                                    ga['n_mouse']))
print('  mouse-only, complete (%d = %d - %d): %s'
      % (len(ga['mouse_only_complete'].split()), ga['n_mouse'],
         ga['n_overlap_map_gap_corrected'], ga['mouse_only_complete']))
print('    of which mapped: %s ; unmapped and genuinely mouse-only: %s'
      % (ga['mouse_only'], ga['mouse_only_real'] or '-'))
print('  human-only, real (%d): %s' % (len(ga['human_only_real'].split()),
                                       ga['human_only_real']))
print('  human-only, raw (uncorrected, DO NOT QUOTE): %s' % ga['human_only'])
asym.append(dict(set='A_SENDER_FINAL_strict_FROZEN', **ga))

print('\nTier C, both directions')
lm, lh = m('C_ligands'), h('C_ligands')
gc = gap_split(lm, lh)
lm2h = {ORTHO[g] for g in lm if g in ORTHO}
print('  mouse ligands %d -> mapped %d ; human ligands %d' % (len(lm), len(lm2h), len(lh)))
ho = sorted(lh - lm2h)
# Separate a REAL cross-species gap from a gap in the pinned MGI map: if the Title-case mouse form
# of the symbol is already in the mouse ligand list, the arms DO both carry the gene and only the
# map is missing the row.
map_gap = [g for g in ho if g.capitalize() in lm]
real = [g for g in ho if g.capitalize() not in lm]
print('  human-only ligands (raw): %s' % ' '.join(ho))
print('    of which MGI-MAP GAPS (gene is on BOTH panels, the pinned map lacks the row): %s'
      % (' '.join(map_gap) or '-'))
print('    of which REAL arm asymmetries: %s' % (' '.join(real) or '-'))
print('  mouse-only ligands (after mapping): %s' % ' '.join(sorted(lm2h - lh)))
print('  CXCL8 in the MGI ortholog map at all? %s'
      % ('yes' if 'CXCL8' in set(ORTHO.values()) else 'NO -- no mouse ortholog exists'))
assert map_gap == gc['human_only_map_gap'].split() and real == gc['human_only_real'].split()
print('  mouse-only ligands, complete (%d = %d - %d): %s'
      % (len(gc['mouse_only_complete'].split()), gc['n_mouse'],
         gc['n_overlap_map_gap_corrected'], gc['mouse_only_complete']))
asym.append(dict(set='C_ligands', **gc))

json.dump(dict(mouse_panel=len(MPANEL), human_panel=len(HPANEL),
               ortholog_intersected_panel=len(XPANEL_H),
               ortholog_intersected_panel_mouse_side=len(XPANEL_M),
               ortho_map_rows=len(ORTHO), ortho_map_distinct_human=_n_h_distinct,
               ortho_map_human_with_multiple_mouse=_n_h_multi,
               mouse_panel_no_ortholog=len(MPANEL) - len(m2h_on),
               human_panel_unreachable=len(HPANEL) - len(XPANEL_H),
               table=rows, asymmetry=asym),
          open(RES + '/crossarm_geneset_table.json', 'w'), indent=1)
print('\nWritten to %s/crossarm_geneset_table.{csv,json}' % RES)
