#!/usr/bin/env python3
"""
Phase 8 freeze -- MOUSE counterpart of the adopted C6 decision, so the two arms stay comparable.

The PI adopted the re-sourced B7 (SenMayo u Reactome SASP, minus the Tier A caller) for H1.
Applying it to the human arm alone would leave Phase 7 section 17 comparing a 116-gene human B7
against a 38-gene mouse B7 and a 33-gene human Tier A against a 25-gene mouse one. This script
rebuilds the MOUSE side by the SAME method and re-runs the section 11 gate on the mouse panel.

SOURCES -- pinned, nothing re-downloaded:
  genesets/msigdb_mouse_2026.1.Mm/   MSigDB release 2026.1.Mm, fetched 2026-08-20.
                                     It already contains MSigDB's OWN MOUSE versions of both
                                     source sets (SAUL_SEN_MAYO = MM16098, REACTOME SASP =
                                     MM14900), so NO ortholog mapping is needed for B7 -- the same
                                     discipline the mouse arm used for B1-B6.
  XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv + GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz
                                     5,006 + 100 = 5,106 Gene Expression features, of which 9 are
                                     genotyping probes -> 5,097 genes. Both counts asserted below.

WHY A SEPARATE OUTPUT DIRECTORY. The frozen mouse sets in genesets/*.txt are the ones every
existing Phase 2-5 result was computed under. Overwriting them would silently invalidate published
numbers. The C6 mouse sets are written to genesets/mouse_c6/ instead. Promoting them means
RE-FITTING the mouse arm, which is a PI decision, not a file copy.

NOTE, reported not fixed: code/build_genesets.py cannot be re-run as committed -- its SCRATCH
constant points at a per-session /tmp path that no longer exists, so MSigDB would load empty and
it would clobber genesets/*.txt with empty modules. This script reads the ARCHIVED pin instead,
and reproduces the mouse arm's published Tier A/B sizes exactly as a self-check.

Run: python3 /workspace/code/build_genesets_mouse_c6.py
"""
import csv, gzip, json, os
from collections import OrderedDict

W = '/workspace'
MS_DIR = W + '/genesets/msigdb_mouse_2026.1.Mm'
OUT = W + '/genesets/mouse_c6'
RES = W + '/results/phase7_jobA'
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ panel
panel_meta = {r['gene_name']: r for r in
              csv.DictReader(open(W + '/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv'))}
p100 = [r['Gene'] for r in csv.DictReader(
        gzip.open(W + '/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz', 'rt'))]
GENOTYPING = {g for g in p100 if ('_WT' in g or '_ALT' in g or '_del_' in g or '_splice_' in g)}
FEATURES = set(panel_meta) | set(p100)
PANEL = FEATURES - GENOTYPING
assert len(FEATURES) == 5106, len(FEATURES)
assert len(GENOTYPING) == 9, len(GENOTYPING)
assert len(PANEL) == 5097, len(PANEL)

# ARCHIVE INTEGRITY, reported not silently skipped. Three files in the pinned mouse archive are
# HTML error pages, not JSON: MSigDB has no MOUSE version of FRIDMAN_SENESCENCE_UP /
# FRIDMAN_SENESCENCE_DN / WP_NRF2_PATHWAY, and the 2026-08-20 fetch saved the 'not found' page.
# code/build_genesets.py hid this behind a bare `except: pass`. None of the three feeds any tier
# (B6 uses HALLMARK ROS + REACTOME detox + a CURATED NRF2 list, not WP_NRF2_PATHWAY), so no
# published mouse number is affected -- but the failure must be visible, not swallowed.
MS, MS_BAD = {}, []
for f in sorted(os.listdir(MS_DIR)):
    try:
        d = json.load(open(os.path.join(MS_DIR, f)))
    except json.JSONDecodeError:
        MS_BAD.append(f)
        continue
    k = list(d)[0]
    MS[k] = (set(d[k]['geneSymbols']), d[k]['systematicName'])

def on_panel(gs):
    return sorted(set(gs) & PANEL)

# ------------------------------------------------------------------ tier definitions
# Copied VERBATIM from code/build_genesets.py so the mouse arm is reproduced, not re-invented.
TIER_A_SUB = OrderedDict([
 ('A1_core_arrest', """Cdkn1a Cdkn2a Cdkn2b Trp53 Trp53i3 Gadd45a Gadd45b Gadd45g Cdkn1b Cdkn2c
    Cdkn2d Rb1 Rbl1 Rbl2 Ccnd1 Glb1 Cdkn1c""".split()),
 ('A2_proliferation_down', """Mki67 Top2a Pcna Ccnb1 Ccnb2 Ccna2 Cdk1 Birc5 Tyms Rrm2 Mcm2 Mcm3
    Mcm4 Mcm5 Mcm6 Mcm7 Aurkb Plk1 Foxm1 Cdc20 Bub1""".split()),
 ('A3_nuclear_chromatin', """Lmnb1 Lmnb2 Hmgb1 Hmgb2 Hmga1 Hist1h1c H1f0 Hp1bp3 Chaf1a Chaf1b
    Cbx5 Lbr Ehmt2""".split()),
 ('A4_dna_damage_response', """Atm Atr Chek1 Chek2 H2ax H2afx Trp53bp1 Mdm2 Rad51 Brca1 Brca2
    Parp1 Ercc1 Nbn Mre11a Rad50 Xrcc5 Xrcc6 Rpa1 Ddb2 Xpc Rif1 Mdc1 Ung Apex1""".split()),
 ('A5_senescence_curated_nonsecreted', """Tnfrsf10b Bcl2 Bcl2l1 Mdm4 Sirt1 Nfe2l2 Ezh2 Suv39h1
    Terf2 Tert Cdc25a E2f1 Foxo3 Trp63 Trp73 Nfatc1 Igfbp5 Tp53inp1 Zmat3 Bax Perp Aen Sesn1
    Sesn2 Eda2r Phlda3 Ccng1 Btg2 Plk3 Trim22 Nupr1""".split()),
])
TIER_B_SRC = OrderedDict([
 ('B1_tnfa_nfkb_proximal',  ['HALLMARK_TNFA_SIGNALING_VIA_NFKB']),
 ('B2_il6_jak_stat3',       ['HALLMARK_IL6_JAK_STAT3_SIGNALING']),
 ('B3_interferon_response', ['HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE']),
 ('B4_downstream_arrest',   ['HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MYC_TARGETS_V1']),
 ('B5_emt_ecm',             ['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']),
 ('B6_oxidative_stress',    ['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                             'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES']),
])
B6_EXTRA_NRF2 = ('Nqo1 Hmox1 Gclm Gclc Txnrd1 Srxn1 Slc7a11 Gsta3 Gstm1 Gstp1 Ftl1 Fth1 Mgst1 Cbr3 '
                 'Akr1b3 Keap1 Nfe2l2 Maff Mafg Sesn2 Txn1 Prdx1 Gpx4 Cat Sod1 Sod2 Gsr').split()
B7_V1_CUR = """Cdkn1a Cdkn2a Cdkn2b Trp53 Gadd45a Mdm2 Ccnd1 Glb1 Lmnb1 Hmgb1 Hmgb2 Serpine1 Il1a
 Il1b Il6 Cxcl1 Cxcl2 Ccl2 Ccl20 Mmp3 Mmp12 Timp1 Igfbp3 Igfbp7 Gdf15 Tnfrsf10b Bcl2l1 Trp53i3
 Sirt1 Ezh2 Suv39h1 Chek2 Atm Trp53bp1 Ddb2 Xpc Tgfb1 Thbs1 Nfkbia Junb Fos Egr1 Sod2 Cdkn1b
 Trp53inp1 Zmat3 Bax Eda2r Phlda3 Ccng1""".split()
# The section 10 sixteen, in mouse symbols (the gotcha map of section 10, read backwards).
A_S10_16_MOUSE = ('Cdkn1a Cdkn2a Cdkn2b Trp53 Trp53i3 Gadd45a Lmnb1 Hmgb1 Hmgb2 Atm Atr Chek1 '
                  'Chek2 H2ax Trp53bp1 Mdm2').split()
B7_SRC = ['SAUL_SEN_MAYO', 'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP']

print('=' * 100)
print('MOUSE C6 REBUILD -- GSE310392 panel: %d Gene Expression features, %d genotyping probes, '
      '%d genes' % (len(FEATURES), len(GENOTYPING), len(PANEL)))
print('MSigDB mouse pin: release 2026.1.Mm, fetched 2026-08-20 (%d archived sets)' % len(MS))
if MS_BAD:
    print('ARCHIVE INTEGRITY: %d archived files are HTML error pages, not JSON (no mouse version'
          ' exists upstream): %s' % (len(MS_BAD), ', '.join(MS_BAD)))
    print('   None feeds any tier; no published mouse number is affected. Reported, not hidden.')
print('=' * 100)

A_sub = OrderedDict((k, on_panel(v)) for k, v in TIER_A_SUB.items())
A0 = sorted(set().union(*A_sub.values()))
A_S10 = on_panel(A_S10_16_MOUSE)
print('\nTier A candidate pool A0 (on-panel, pre-disjointness): %d' % len(A0))
print('Section 10 sixteen in mouse symbols, on-panel: %d/16  (off-panel: %s)'
      % (len(A_S10), ' '.join(sorted(set(A_S10_16_MOUSE) - PANEL)) or '-'))

B, B_prov = OrderedDict(), OrderedDict()
for mod, srcs in TIER_B_SRC.items():
    g = set()
    for s in srcs:
        assert s in MS, s
        g |= MS[s][0]
    if mod == 'B6_oxidative_stress':
        g |= set(B6_EXTRA_NRF2)
    B[mod] = on_panel(g)
    B_prov[mod] = list(srcs) + (['CUR_NRF2_targets'] if mod == 'B6_oxidative_stress' else [])

B7v1 = on_panel(B7_V1_CUR)
_b7 = set()
for s in B7_SRC:
    assert s in MS, s
    _b7 |= MS[s][0]
B7_SOURCED_FULL = on_panel(_b7)
B['B7_secondary_senescence'] = sorted(set(B7_SOURCED_FULL) - set(A0))
B_prov['B7_secondary_senescence'] = B7_SRC + ['minus Tier A caller pool A0 (C6, adopted)']

print('\n### Self-check: does this reproduce the mouse arm as published (genesets/README.md)? ###')
PUB = {'B1_tnfa_nfkb_proximal': 126, 'B2_il6_jak_stat3': 68, 'B3_interferon_response': 100,
       'B4_downstream_arrest': 190, 'B5_emt_ecm': 125, 'B6_oxidative_stress': 31}
ok = True
for k, want in PUB.items():
    got = len(B[k])
    ok &= got == want
    print('  %-30s published %3d  recomputed %3d  %s' % (k, want, got, 'OK' if got == want else 'MISMATCH'))
print('  %-30s published %3d  recomputed %3d  %s' % ('B7 v1 (curated)', 38, len(B7v1),
                                                     'OK' if len(B7v1) == 38 else 'MISMATCH'))
print('  %-30s published %3d  recomputed %3d  %s' % ('A0 pre-disjointness', 74, len(A0),
                                                     'OK' if len(A0) == 74 else 'MISMATCH'))
assert ok and len(B7v1) == 38 and len(A0) == 74, 'mouse arm NOT reproduced -- stop and investigate'

print('\n### Mouse Tier B under C6 ###')
for k, v in B.items():
    print('  %-30s %3d on-panel   [%s]' % (k, len(v), '; '.join(B_prov[k])))
print('  B7 sourced union before caller subtraction: %d on-panel (SenMayo %d + Reactome SASP %d)'
      % (len(B7_SOURCED_FULL), len(on_panel(MS['SAUL_SEN_MAYO'][0])),
         len(on_panel(MS[B7_SRC[1]][0]))))
print('  B7 v1 curated (superseded, archived): %d' % len(B7v1))

# ------------------------------------------------------------------ section 11 gate, mouse
B_UNION = set().union(*B.values())
A_STRICT = sorted(set(A0) - B_UNION)
A_per_mod = OrderedDict((k, sorted(set(A0) - set(v))) for k, v in B.items())
MODS = list(B)

print('\n' + '=' * 100)
print('SECTION 11 GATE -- MOUSE ARM, C6 CONFIGURATION')
print('=' * 100)
gate_rows = []
for aname, A in (('A_S10_16_mouse', set(A_S10)), ('A_ported(A0)', set(A0))):
    c1 = len(A) >= 15
    print('\n### Tier A definition: %s  (|A & panel| = %d) ###' % (aname, len(A)))
    print('  [%s] len(A & panel) >= 15 : %d' % ('PASS' if c1 else 'FAIL', len(A)))
    c2 = True
    for k in MODS:
        okk = len(B[k]) >= 30
        c2 &= okk
        print('  [%s] len(%s) >= 30 : %d' % ('PASS' if okk else 'FAIL', k, len(B[k])))
    removed = OrderedDict((k, sorted(A & set(B[k]))) for k in MODS if A & set(B[k]))
    print('  [%s] A n B_k == 0 for all k : %d memberships'
          % ('PASS' if not removed else 'FAIL', sum(len(v) for v in removed.values())))
    for k, v in removed.items():
        print('     removed by %-30s (%2d): %s' % (k, len(v), ' '.join(v)))
    Afin = A - B_UNION
    c3 = len(Afin) >= 15
    print('  [%s] re-check len(A) >= 15 after removal : %d' % ('PASS' if c3 else 'FAIL', len(Afin)))
    print('     surviving: %s' % (' '.join(sorted(Afin)) or '(EMPTY)'))
    v = 'PASS' if (c1 and c2 and c3) else 'FAIL'
    print('  ==> MOUSE GATE VERDICT for %s: %s' % (aname, v))
    gate_rows.append(dict(arm='mouse', tierA=aname, n_A=len(A), A_ge15=c1, B_ge30=c2,
                          n_removed=len(A) - len(Afin), n_A_final=len(Afin), verdict=v))

print('\n### PRE-REGISTERED SENSITIVITY: per-module mouse sender sets ###')
pm_pass = True
for k in MODS:
    s_ = set(A_per_mod[k])
    okk = len(s_) >= 15 and not (s_ & set(B[k]))
    pm_pass &= okk
    print('  %-30s |A_mod|=%3d  A_mod n B_k=%d  %s' % (k, len(s_), len(s_ & set(B[k])),
                                                       'PASS' if okk else 'FAIL'))
print('  [%s] per-module gate' % ('PASS' if pm_pass else 'FAIL'))

# ------------------------------------------------------------------ panel-definition audit
# The 5,006-row pan-tissue metadata CSV is NOT "the panel": it is ONE of the two panel files.
# The measured panel is the h5 Gene Expression feature list, which genesets/README.md verified
# equals metadata(5,006) u custom100(100) exactly, 0 discrepancies either way. Dropping the 9
# genotyping probes (not genes) leaves 5,097. Every downstream script that touches the matrix
# (annotate_pipeline.py, build_random_null_sets.py) uses exactly that ENSMUSG-filtered set.
print('\n' + '=' * 100)
print('PANEL DEFINITION AUDIT (mouse)')
print('=' * 100)
csv_only = set(panel_meta)
print('  5K pan-tissue metadata CSV       : %d rows' % len(csv_only))
print('  GSE310392 custom add-on          : %d (%d genotyping probes, %d genes)'
      % (len(p100), len(GENOTYPING), len(p100) - len(GENOTYPING)))
print('  overlap between the two files    : %d' % len(csv_only & set(p100)))
print('  AUTHORITATIVE PANEL (h5 Gene Expression minus genotyping probes): %d' % len(PANEL))
print('  Module sizes differ between the two definitions where the add-on contributes:')
for k, v in B.items():
    a, b = len(v), len(set(v) & csv_only)
    if a != b:
        print('    %-30s %3d on the authoritative panel vs %3d on the CSV alone  (add-on genes: %s)'
              % (k, a, b, ' '.join(sorted(set(v) - csv_only))))
print('\n### MODULE MARGIN OVER THE >=30 FLOOR ###')
margins = []
for k, v in B.items():
    m_auth, m_csv = len(v) - 30, len(set(v) & csv_only) - 30
    margins.append(dict(arm='mouse', module=k, n_authoritative=len(v), margin_authoritative=m_auth,
                        n_csv_only=len(set(v) & csv_only), margin_csv_only=m_csv))
    flag = '  <-- NO MARGIN' if m_auth <= 1 else ''
    print('  %-30s n=%3d margin=%+3d   (CSV-only n=%3d margin=%+3d)%s'
          % (k, len(v), m_auth, len(set(v) & csv_only), m_csv, flag))

# ------------------------------------------------------------------ intersection matrices
A_ROWS = OrderedDict()
A_ROWS['A0_pre_disjointness'] = set(A0)
A_ROWS['A_S10_16_mouse (variant, FAILED)'] = set(A_S10)
for k in MODS:
    A_ROWS['A_sender_for_%s (sensitivity)' % k.split('_', 1)[1]] = set(A_per_mod[k])
A_ROWS['A_SENDER_FINAL_strict (PRIMARY)'] = set(A_STRICT)
BLAB = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
rows_ab = []
for an, av in A_ROWS.items():
    r = OrderedDict(tier_A_set=an, n_A_on_panel=len(av))
    tot = 0
    for k, lab in zip(MODS, BLAB):
        n = len(av & set(B[k]))
        tot += n
        r['%s_%s' % (lab, k.split('_', 1)[1])] = n
    r['total_memberships'] = tot
    rows_ab.append(r)
with open(RES + '/intersection_matrix_mouse.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_ab[0]))
    w.writeheader()
    for r in rows_ab:
        w.writerow(r)
rows_bb = []
for k, lab in zip(MODS, BLAB):
    r = OrderedDict(module='%s_%s' % (lab, k.split('_', 1)[1]), n_on_panel=len(B[k]))
    for k2, lab2 in zip(MODS, BLAB):
        r[lab2] = len(set(B[k]) & set(B[k2]))
    rows_bb.append(r)
with open(RES + '/intersection_matrix_mouse_BxB.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_bb[0]))
    w.writeheader()
    for r in rows_bb:
        w.writerow(r)
# same two matrices under the SUPERSEDED curated B7, so the figure can show what C6 changed
B_PRE = dict(B)
B_PRE['B7_secondary_senescence'] = B7v1
rows_bb_pre = []
for k, lab in zip(MODS, BLAB):
    r = OrderedDict(module='%s_%s' % (lab, k.split('_', 1)[1]), n_on_panel=len(B_PRE[k]))
    for k2, lab2 in zip(MODS, BLAB):
        r[lab2] = len(set(B_PRE[k]) & set(B_PRE[k2]))
    rows_bb_pre.append(r)
with open(RES + '/intersection_matrix_mouse_BxB_preC6.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_bb_pre[0]))
    w.writeheader()
    for r in rows_bb_pre:
        w.writerow(r)
print('\nwrote intersection_matrix_mouse.csv, _BxB.csv, _BxB_preC6.csv')

# ------------------------------------------------------------------ write
def write(n, g):
    with open(os.path.join(OUT, n + '.txt'), 'w') as fh:
        fh.write('\n'.join(g) + '\n')

write('A_SENDER_FINAL_strict', A_STRICT)
for k, v in A_per_mod.items():
    write('A_sender_for_%s' % k.split('_', 1)[1], v)
for k, v in B.items():
    write('B_%s' % k.split('_', 1)[1], v)
os.makedirs(OUT + '/variants', exist_ok=True)
for n, g in (('B_secondary_senescence_v1_curated', B7v1), ('A_S10_16_mouse', A_S10),
             ('B_secondary_senescence_sourced_before_caller_subtraction', B7_SOURCED_FULL)):
    with open(os.path.join(OUT, 'variants', n + '.txt'), 'w') as fh:
        fh.write('\n'.join(g) + '\n')

json.dump(dict(panel_features=len(FEATURES), panel_genes=len(PANEL),
               genotyping_probes=sorted(GENOTYPING),
               msigdb='2026.1.Mm (pinned 2026-08-20)', msigdb_bad_archive_files=MS_BAD,
               A0=A0, A_S10_16_on_panel=A_S10, A_SENDER_FINAL_strict=A_STRICT,
               A_per_module_sizes={k: len(v) for k, v in A_per_mod.items()},
               B_sizes={k: len(v) for k, v in B.items()}, B_prov=B_prov,
               B7_v1_curated_n=len(B7v1), B7_sourced_full_n=len(B7_SOURCED_FULL),
               gate=gate_rows, per_module_gate_pass=pm_pass, margins=margins,
               panel_csv_only=len(set(panel_meta))),
          open(RES + '/mouse_c6_summary.json', 'w'), indent=1)
print('\n|A_SENDER_FINAL_strict| mouse: %d  (was 25 under the curated B7)' % len(A_STRICT))
print('Written to %s (+ variants/) and %s/mouse_c6_summary.json' % (OUT, RES))
