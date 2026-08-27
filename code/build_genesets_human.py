#!/usr/bin/env python3
"""
SASP Spatial Response Kernel -- Phase 7 Part III / Job A (HUMAN arm, H1).
Human port of build_genesets.py. Same tier structure, same file naming, same output format,
so nothing downstream needs to change when it is pointed at genesets/human/ instead of genesets/.

Target panel: GSE326743 human Xenium = Xenium Prime 5K Human + 100 addon
              5,093 'Gene Expression' features, verified on GSM9638040's cell_feature_matrix.h5
              and cross-checked identical against GSM9638044 / GSM9638046.
              See genesets/h1_candidate/PROVENANCE.md.

Gene set sources (tags mirror genesets/README.md):
  MSIGDB_2026.1.Hs : MSigDB HUMAN collections (H / C2 / C5), fetched from gsea-msigdb.org
                     JSON API 2026-08-27, release 2026.1.Hs. Raw JSON archived in
                     genesets/msigdb_human_2026.1.Hs/. NOT recalled from memory.
  CUR              : the mouse arm's curated lists (genesets/README.md), translated to human
                     symbols through the pinned MGI ortholog map + HGNC. No new curation.
  PHASE7_S10       : the literal gene list printed in Phase 7 doc section 10.
  PANEL            : derived from the panel feature-type column (Tier E controls) or, for the
                     'Secreted' constraint, transferred from the MOUSE panel metadata through
                     the ortholog map (the human Prime 5K metadata csv is not on disk).
  ORTHO            : genesets/mouse_human_orthologs_MGI.csv (already pinned and used by
                     run_deepscence.py, caller_disagree2.py, run_phase3_n8.py).
  HGNC             : genesets/hgnc_pin/, for alias / previous-symbol -> panel-symbol resolution.

Run:  python3 /workspace/code/build_genesets_human.py
No network access required to re-run.
"""
import csv, gzip, json, os, glob, sys
from collections import OrderedDict

GS   = '/workspace/genesets'
OUT  = GS + '/human'
MSDIR = GS + '/msigdb_human_2026.1.Hs'
MSIGDB_RELEASE, MSIGDB_FETCHED = '2026.1.Hs', '2026-08-27'
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ panel
PANEL_ROWS = list(csv.DictReader(open(GS + '/h1_candidate/GSE326743_gene_panel_5093.csv')))
PANEL = {r['gene_name'] for r in PANEL_ROWS}
ENSG2SYM = {r['gene_id']: r['gene_name'] for r in PANEL_ROWS}
assert len(PANEL) == 5093, len(PANEL)

FEATURES = list(csv.DictReader(open(GS + '/h1_candidate/GSE326743_panel_features.csv')))
FEAT_BY_TYPE = OrderedDict()
for r in FEATURES:
    FEAT_BY_TYPE.setdefault(r['feature_type'], []).append(r['gene_name'])

# ------------------------------------------------------------------ symbol authorities
ORTHO = {r['mouse_symbol']: r['human_symbol']
         for r in csv.DictReader(open(GS + '/mouse_human_orthologs_MGI.csv'))}

# Symbol resolution lives in human_symbols.py so this build and build_markers_human_spleen.py
# resolve identically. See that module for why a naive symbol/alias/prev union is wrong.
sys.path.insert(0, '/workspace/code')
import human_symbols as HS
_APPROVED, _PREV, _ALIAS = HS.APPROVED, HS.PREV, HS.ALIAS
RESOLUTIONS = HS.RESOLUTIONS
resolve = HS.resolve

def m2h(mouse_symbols):
    """Mouse symbol list -> human symbols, MGI map first, upper-case fallback recorded."""
    out, unmapped = [], []
    for g in mouse_symbols:
        h = ORTHO.get(g)
        if h is None:
            cand = g.upper()
            if cand in PANEL or cand in _APPROVED or cand in _PREV or cand in _ALIAS:
                RESOLUTIONS.append((g, cand, 'NOT in MGI map; upper-case form is a known HGNC symbol'))
                h = cand
            else:
                unmapped.append(g)
                continue
        out.append(h)
    return out, unmapped

UNMAPPED = {}

def on_panel(gs):
    return sorted({r for r in (resolve(g) for g in gs) if r})

# ------------------------------------------------------------------ msigdb
MS = {}
for f in glob.glob(MSDIR + '/HALLMARK_*.json') + glob.glob(MSDIR + '/REACTOME_*.json') + \
        glob.glob(MSDIR + '/GOBP_*.json') + glob.glob(MSDIR + '/SAUL_*.json') + \
        glob.glob(MSDIR + '/FRIDMAN_*.json') + glob.glob(MSDIR + '/WP_*.json'):
    d = json.load(open(f))
    k = list(d)[0]
    MS[k] = (set(d[k]['geneSymbols']), d[k]['systematicName'], d[k].get('collection', ''))
H_ALL = json.load(open(MSDIR + '/_h.all.v2026.1.Hs.json'))   # full H collection, 50 sets

# ================================================================== TIER A
# A_PHASE7_16: the literal Phase 7 section 10 list. This is the doc's Tier A.
A_PHASE7_16 = """CDKN1A CDKN2A CDKN2B TP53 TP53I3 GADD45A LMNB1 HMGB1 HMGB2 ATM ATR CHEK1
   CHEK2 H2AX TP53BP1 MDM2""".split()
assert len(A_PHASE7_16) == 16

# The mouse arm's five Tier A subsets, ported. Mouse source: build_genesets.py TIER_A_SUB.
TIER_A_SUB_MOUSE = OrderedDict([
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
TIER_A_SUB = OrderedDict()
for k, v in TIER_A_SUB_MOUSE.items():
    h, um = m2h(v)
    TIER_A_SUB[k] = h
    if um:
        UNMAPPED[k] = um
# The section 10 list is a subset of the ported mouse subsets by construction; assert it.
_ported_union = set().union(*TIER_A_SUB.values())
A16_NOT_IN_PORT = sorted(set(A_PHASE7_16) - _ported_union)

# ================================================================== TIER B
# Source sets mirror the mouse arm exactly (cross-arm comparability, Phase 7 test A8).
# Phase 7 s10 names 7 HALLMARK sets; E2F_TARGETS and G2M_CHECKPOINT both feed B4, which is
# why 7 named HALLMARK sets + secondary_senescence = 7 modules, not 8.
TIER_B_SRC = OrderedDict([
 ('B1_tnfa_nfkb_proximal',  ['HALLMARK_TNFA_SIGNALING_VIA_NFKB']),
 ('B2_il6_jak_stat3',       ['HALLMARK_IL6_JAK_STAT3_SIGNALING']),
 ('B3_interferon_response', ['HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE']),
 ('B4_downstream_arrest',   ['HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MYC_TARGETS_V1']),
 ('B5_emt_ecm',             ['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']),
 ('B6_oxidative_stress',    ['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY']),
])
# Mouse B6 needed a rescue (union with REACTOME detox + curated NRF2) because HALLMARK ROS alone
# gave 17 on-panel. Whether the human panel needs the same rescue is decided BY THE NUMBER below,
# not assumed; the rescue components are kept ready and applied only if ROS alone misses >= 30.
B6_RESCUE_SRC = ['REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES']
B6_EXTRA_NRF2_MOUSE = ('Nqo1 Hmox1 Gclm Gclc Txnrd1 Srxn1 Slc7a11 Gsta3 Gstm1 Gstp1 Ftl1 Fth1 Mgst1 '
                       'Cbr3 Akr1b3 Keap1 Nfe2l2 Maff Mafg Sesn2 Txn1 Prdx1 Gpx4 Cat Sod1 Sod2 Gsr').split()
# B7 secondary senescence: the mouse arm's curated list, ported. Phase 7 s23/C6 flags this module
# as problematic (shares genes with the Tier A caller, sole citation is a meeting abstract).
# Port faithfully and REPORT the overlap; fixing it is C6 and out of scope for Job A.
B7_CUR_MOUSE = """Cdkn1a Cdkn2a Cdkn2b Trp53 Gadd45a Mdm2 Ccnd1 Glb1 Lmnb1 Hmgb1 Hmgb2 Serpine1 Il1a
 Il1b Il6 Cxcl1 Cxcl2 Ccl2 Ccl20 Mmp3 Mmp12 Timp1 Igfbp3 Igfbp7 Gdf15 Tnfrsf10b Bcl2l1 Trp53i3
 Sirt1 Ezh2 Suv39h1 Chek2 Atm Trp53bp1 Ddb2 Xpc Tgfb1 Thbs1 Nfkbia Junb Fos Egr1 Sod2 Cdkn1b
 Trp53inp1 Zmat3 Bax Eda2r Phlda3 Ccng1""".split()

# ================================================================== TIER C
# Mouse pairs (build_genesets.py TIER_C_PAIRS) ported, PLUS CXCL8 which has no mouse ortholog.
TIER_C_PAIRS = [
 ('IL6',    ['IL6R', 'IL6ST'],            'canonical SASP; drives B2. HUMAN receptor is IL6R (mouse Il6ra).'),
 ('CXCL8',  ['CXCR1', 'CXCR2'],           'HUMAN-ONLY. s10 addition: IL-8, NO MOUSE ORTHOLOG. CXCR1 also has no clean mouse ortholog.'),
 ('CXCL1',  ['CXCR2'],                    'kept for cross-arm parity; in mouse this was a CXCL8 analogue, in human it is a gene in its own right'),
 ('CXCL2',  ['CXCR2'],                    'kept for cross-arm parity; mouse analogue of CXCL8'),
 ('CXCL5',  ['CXCR2'],                    'kept for cross-arm parity; ELR+ CXC chemokine'),
 ('CCL2',   ['CCR2', 'ACKR3'],            'HIGHEST PRIORITY (Sec 9). ACKR3=CXCR7.'),
 ('CXCL12', ['CXCR4', 'ACKR3', 'DPP4'],   'DPP4 cleaves CXCL12 -> candidate range-limiting mechanism.'),
 ('TGFB1',  ['TGFBR1', 'TGFBR2'],         'contact-adjacent; short lambda expected'),
 ('IL1A',   ['IL1R1', 'IL1RAP'],          'IL-1alpha LARGELY MEMBRANE-BOUND -> SHORTEST lambda predicted'),
 ('IL1B',   ['IL1R1', 'IL1RAP'],          'secreted; longer lambda than IL1A expected'),
 ('TNF',    ['TNFRSF1A', 'TNFRSF1B'],     'drives B1'),
 ('IGFBP3', [],                           'SASP Atlas core secreted factor'),
 ('IGFBP7', [],                           'SASP Atlas core secreted factor'),
 ('GDF15',  [],                           'SASP Atlas core secreted factor'),
 ('MMP1',   [],                           'human single gene; mouse orthologs were Mmp1a/Mmp1b'),
 ('MMP3',   [],                           'SASP Atlas core secreted factor'),
 ('TIMP1',  [],                           'SASP Atlas core secreted factor'),
 ('THBS1',  ['CD47', 'SDC1'],             'paper-specific: THBS-mediated senescence network (Karpova 2026)'),
]

# ================================================================== TIER D
# The mouse D_zonation_* sets were DATA-DERIVED from the GSM9295284 expression matrix.
# H1 is SPLEEN and is still behind the section 15 freeze -- no H1 expression value has been read.
# So the arm-specific anatomical covariate CANNOT be built here; only the tissue-generic
# vascular/stromal landmark set can be ported. See the report.
D_LANDMARK_MOUSE = """Pecam1 Cdh5 Vwf Stab2 Lyve1 Kdr Eng Aqp1 Efnb2 Ephb4 Gja5 Krt19 Krt7 Krt8
   Krt18 Epcam Sox9 Spp1 Pkhd1 Onecut1 Hnf1b Cftr Acta2 Myh11 Tagln Des Dcn Col1a1 Col1a2 Lrat
   Pdgfrb Rgs5 Notch3""".split()
D_LIVER_ONLY = 'Krt19 Krt7 Pkhd1 Onecut1 Hnf1b Cftr Sox9 Lrat Stab2'.split()   # bile duct / liver-specific
D_COVARIATES = [
 'cell_type', 'density_25um', 'density_50um', 'density_100um', 'knn_composition',
 'total_counts', 'n_genes_detected', 'cell_area_um2', 'nucleus_area_um2', 'section',
 'segmentation_method', 'dist_to_tissue_boundary_um',
 'ANATOMICAL_ARM_SPECIFIC__spleen_red_pulp_white_pulp_axis__DEFERRED_see_report',
]

# ================================================================== TIER E
TIER_E_MOUSE = OrderedDict([
 ('E1_housekeeping_planA', 'Actb Gapdh Rpl13a Rps18 Tbp Ppia'.split()),
 ('E1b_housekeeping_expanded', """Actb Gapdh Rpl13a Rps18 Tbp Ppia Rplp0 Hprt Pgk1 Sdha Ywhaz
    Rpl11 Tubb5 B2m Eef1a1 Ubc Psmb2 Vcp Rab7 Chmp2a Emc7 Gpi1 Reep5 Snrpd3 Vps29 Cnot4
    Rer1 Tmem59 Actr2 Arpc2 Cfl1 Pfn1 Rack1 Hnrnpk Srp14 Sar1a Atp5f1b Cox4i1 Ndufb8""".split()),
])
# Phase 7 s10 gives the E1 housekeeping six directly in human symbols; use those verbatim.
E1_PHASE7 = 'ACTB GAPDH RPL13A RPS18 TBP PPIA'.split()

# ------------------------------------------------------------------ build
print('=' * 100)
print('PANEL: %d Gene Expression features (GSE326743, Xenium Prime 5K Human + 100 addon)' % len(PANEL))
for t, v in FEAT_BY_TYPE.items():
    print('   %-28s %5d' % (t, len(v)))
print('=' * 100)

# --- 'no secreted factors' constraint, transferred from the MOUSE panel metadata (deviation)
mouse_meta = {r['gene_name']: r for r in
              csv.DictReader(open('/workspace/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv'))}
SECRETED_HUMAN = set()
for mg, m in mouse_meta.items():
    if 'Secreted' in (m.get('location') or ''):
        h = ORTHO.get(mg)
        if h and h in PANEL:
            SECRETED_HUMAN.add(h)
print("Panel genes inheriting location='Secreted' from the MOUSE panel metadata via the "
      "ortholog map: %d" % len(SECRETED_HUMAN))

print('\n### MSigDB HUMAN sets used (release %s, fetched %s) ###' % (MSIGDB_RELEASE, MSIGDB_FETCHED))
for k in sorted(MS):
    print('  %-58s %-7s %-16s n=%3d  on-panel=%3d'
          % (k, MS[k][1], MS[k][2], len(MS[k][0]), len(on_panel(MS[k][0]))))
print('  full H collection archived: _h.all.v%s.json (%d sets)' % (MSIGDB_RELEASE, len(H_ALL)))

# --- Tier A on-panel
print('\n### PHASE 7 s10 TIER A (the literal 16) ###')
a16_panel, a16_off = [], []
for g in A_PHASE7_16:
    r = resolve(g)
    (a16_panel if r else a16_off).append(r or g)
print('  on-panel %d/16 : %s' % (len(a16_panel), ' '.join(sorted(a16_panel))))
print('  OFF-PANEL      : %s' % (' '.join(a16_off) or '-'))
print('  (s10 symbols not present in the ported mouse Tier A: %s)' % (' '.join(A16_NOT_IN_PORT) or 'none'))

print('\n### PORTED MOUSE TIER A subsets, on-panel (BEFORE disjointness enforcement) ###')
A_sub = OrderedDict((k, on_panel(v)) for k, v in TIER_A_SUB.items())
A_all = []
for k, v in A_sub.items():
    miss = sorted({g for g in TIER_A_SUB[k] if not resolve(g)})
    print('  %-38s %2d/%2d on-panel   missing: %s' % (k, len(v), len(set(TIER_A_SUB[k])), ','.join(miss) or '-'))
    if k in UNMAPPED:
        print('      %-34s no MGI ortholog: %s' % ('', ','.join(UNMAPPED[k])))
    A_all += v
A0 = sorted(set(A_all))
print('  TIER A TOTAL (union, on-panel, pre-disjointness): %d' % len(A0))
print('  Sanity: Tier A members inheriting Secreted (must be ~0): %s'
      % (' '.join(sorted(set(A0) & SECRETED_HUMAN)) or 'NONE'))

# --- Tier B
B, B_prov = OrderedDict(), OrderedDict()
for mod, srcs in TIER_B_SRC.items():
    g = set()
    for s in srcs:
        if s in MS:
            g |= MS[s][0]
        else:
            print('!! MISSING MSIGDB SET', s)
    B[mod] = on_panel(g)
    B_prov[mod] = list(srcs)
B6_ALONE = len(B['B6_oxidative_stress'])
if B6_ALONE < 30:
    g = set(MS['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'][0])
    for s in B6_RESCUE_SRC:
        g |= MS[s][0]
    nrf2, _um = m2h(B6_EXTRA_NRF2_MOUSE)
    g |= set(nrf2)
    B['B6_oxidative_stress'] = on_panel(g)
    B_prov['B6_oxidative_stress'] += B6_RESCUE_SRC + ['CUR_NRF2_targets (mouse arm rescue, ported)']
# ---- B7, FROZEN under the PI decision of 2026-08-27 (Phase 7 section 23 / C6) ----------------
# C6 as written -- "rebuild B7 without the genes it shares with the Tier A caller" -- DOES NOT
# WORK on the curated list: it leaves 24 genes against the section 10 caller and 12 against the
# ported caller, both below the section 11 floor of 30 (results/phase7_jobA/b7_c6_rebuild.csv).
# ADOPTED INSTEAD: re-source the module from MSigDB sets that carry peer-reviewed PMIDs, THEN
# subtract the caller. This also replaces the meeting-abstract citation that audit finding B8
# flagged (neretti2024dissecting), and it is closer to section 9's own definition of B7
# ("Tier A minus the calling genes") than the curated list ever was.
#   SAUL_SEN_MAYO                                            PMID 35974106 (saul2022senmayo)
#   REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP   Reactome R-HSA-2559582
# The subtracted caller is A0, the ported Tier A BEFORE disjointness enforcement, i.e. the whole
# candidate sender pool. That makes the fixed point self-consistent: B7 excludes all of A0, so
# A_SENDER_FINAL_strict (a subset of A0) is disjoint from B7 by construction.
B7_SRC = ['SAUL_SEN_MAYO', 'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP']
_b7src = set()
for _s in B7_SRC:
    _b7src |= MS[_s][0]
B['B7_secondary_senescence'] = sorted(set(on_panel(_b7src)) - set(A0))
B_prov['B7_secondary_senescence'] = B7_SRC + ['minus Tier A caller pool A0 (C6, adopted)']

# The superseded v1 stays on disk: section 23 / C6 requires BOTH versions be reported.
B7_h, B7_um = m2h(B7_CUR_MOUSE)
B7_V1_CURATED = on_panel(B7_h)
if B7_um:
    UNMAPPED['B7_secondary_senescence_v1'] = B7_um

B_UNION = set().union(*B.values())
B_UNION_NO_B7 = set().union(*[v for k, v in B.items() if k != 'B7_secondary_senescence'])
B_UNION_NO_B4_B7 = set().union(*[v for k, v in B.items()
                                 if k not in ('B4_downstream_arrest', 'B7_secondary_senescence')])

print('\n### TIER B modules on-panel ###')
for k, v in B.items():
    print('  %-30s %3d on-panel   [%s]' % (k, len(v), '; '.join(B_prov[k])))
print('  B6 note: HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY ALONE gives %d on-panel -> mouse-style '
      'rescue %s' % (B6_ALONE, 'APPLIED' if B6_ALONE < 30 else 'NOT NEEDED'))

# --- disjointness enforcement, exactly as the mouse arm: remove from SENDER, never response
collide = OrderedDict()
for k, v in B.items():
    ov = sorted(set(A0) & set(v))
    if ov:
        collide[k] = ov
print('\n### A x B COLLISIONS on the ported Tier A (genes removed from TIER A) ###')
for k, v in collide.items():
    print('  %-30s %3d : %s' % (k, len(v), ' '.join(v)))
A_FINAL_STRICT   = sorted(set(A0) - B_UNION)
A_FINAL_NO_B7    = sorted(set(A0) - B_UNION_NO_B7)
A_FINAL_NO_B4_B7 = sorted(set(A0) - B_UNION_NO_B4_B7)
print('\n  |A| pre  = %d' % len(A0))
print('  |A| after removing collisions with B1-B6 ONLY (B7 excluded) = %d' % len(A_FINAL_NO_B7))
print('  |A| after removing collisions with B1-B3,B5,B6 (no B4,B7)   = %d' % len(A_FINAL_NO_B4_B7))
print('  |A| after removing collisions with B1-B7 (STRICT)           = %d' % len(A_FINAL_STRICT))

A_per_module = OrderedDict((k, sorted(set(A0) - set(v))) for k, v in B.items())

# --- Tier C
print('\n### TIER C ligand-receptor, on-panel ###')
tierC_rows = []
for lig, recs, note in TIER_C_PAIRS:
    lp = resolve(lig)
    rp = [x for x in (resolve(r) for r in recs) if x]
    rm = [r for r in recs if not resolve(r)]
    tierC_rows.append((lig, lp, rp, rm, note))
    print('%-10s %-8s %-38s %s' % (lig, 'YES' if lp else 'no',
          '%d/%d [%s]%s' % (len(rp), len(recs), ','.join(rp), (' MISSING:' + ','.join(rm)) if rm else ''), note))
print('  CROSS-ARM ASYMMETRY: CXCL8 (+CXCR1) exist here and have NO mouse ortholog. The mouse arm '
      'substituted Cxcl1/Cxcl2/Cxcl5 -> Cxcr2. Any CXCL8 result is HUMAN-ONLY and must not be '
      'reported as replicating a mouse quantity.')

# --- Tier D
D = OrderedDict()
d_land, d_um = m2h([g for g in D_LANDMARK_MOUSE if g not in D_LIVER_ONLY])
_AB = set(A0) | B_UNION
D['D_landmark_vessel_stroma'] = [g for g in on_panel(d_land) if g not in _AB]
D_DROPPED = sorted(set(on_panel(d_land)) & _AB)
print('\n### TIER D, on-panel ###')
print('  D_landmark_vessel_stroma  %2d  (dropped by the over-adjustment guard: %s)'
      % (len(D['D_landmark_vessel_stroma']), ' '.join(D_DROPPED) or '-'))
print('  liver-only genes not ported: %s' % ' '.join(D_LIVER_ONLY))
print('  ANATOMICAL COVARIATE: DEFERRED. The mouse D_zonation_* sets were derived from the mouse '
      'EXPRESSION MATRIX; H1 is spleen and is behind the s15 freeze. Build the red-pulp/white-pulp '
      'axis post-download with the same procedure. Phase 7 test A6 as written specifies a LUNG '
      'axis and does not apply to this arm.')

# --- Tier E
E = OrderedDict()
E['E1_housekeeping_planA'] = on_panel(E1_PHASE7)
e1b, _ = m2h(TIER_E_MOUSE['E1b_housekeeping_expanded'])
E['E1b_housekeeping_expanded'] = on_panel(e1b)
E['E4_negative_control_probes']    = sorted(FEAT_BY_TYPE.get('Negative Control Probe', []))
E['E4_negative_control_codewords'] = sorted(FEAT_BY_TYPE.get('Negative Control Codeword', []))
E['E4_genomic_controls']           = sorted(FEAT_BY_TYPE.get('Genomic Control', []))
print('\n### TIER E controls ###')
for k in ('E1_housekeeping_planA', 'E1b_housekeeping_expanded'):
    src = E1_PHASE7 if k.endswith('planA') else e1b
    print('  %-32s %2d/%2d : %s' % (k, len(E[k]), len(set(src)), ' '.join(E[k])))
    print('     MISSING: %s' % ' '.join(sorted({g for g in src if not resolve(g)})))
for k in ('E4_negative_control_probes', 'E4_negative_control_codewords', 'E4_genomic_controls'):
    print('  %-32s %d features (PANEL feature_type)' % (k, len(E[k])))
print('  E2 cell-type identity: NOT BUILT. See report -- in spleen every abundant cell type is an '
      'immune cell, so the s10 "identity program unrelated to inflammation" has no valid instance '
      'on this panel. HGNC gene group "Hemoglobin subunits" puts only HBG1,HBZ on the panel.')
print('  E3 random size/expression-matched sets: NOT BUILT. build_random_null_sets.py needs the '
      'H1 expression matrix for the expression-matching bins; that is post-freeze.')

# --- SenMayo reference
sm_arrest = []
if 'SAUL_SEN_MAYO' in MS:
    sm = MS['SAUL_SEN_MAYO'][0]
    sm_p = set(on_panel(sm))
    sm_arrest = sorted(sm_p - SECRETED_HUMAN - B_UNION_NO_B7)
    print('\n### SenMayo audit -- SAUL_SEN_MAYO, MSigDB human %s ###' % MS['SAUL_SEN_MAYO'][1])
    print('  set size %d ; on-panel %d ; of those inheriting Secreted: %d'
          % (len(sm), len(sm_p), len(sm_p & SECRETED_HUMAN)))
    print('  SenMayo_arrest = on-panel MINUS secreted MINUS TierB(B1-B6) = %d genes: %s'
          % (len(sm_arrest), ' '.join(sm_arrest)))

# ------------------------------------------------------------------ write files
VAR = OUT + '/variants'
os.makedirs(VAR, exist_ok=True)

def write(name, genes):
    with open(os.path.join(OUT, name + '.txt'), 'w') as fh:
        fh.write('\n'.join(genes) + '\n')

def write_variant(name, genes):
    with open(os.path.join(VAR, name + '.txt'), 'w') as fh:
        fh.write('\n'.join(genes) + '\n')

A_sub_final = OrderedDict((k, sorted(set(v) - B_UNION)) for k, v in A_sub.items())
for k, v in A_sub_final.items():
    write('A_%s' % k.split('_', 1)[1], v)
# FROZEN PRIMARY sender definition (PI decision 2026-08-27).
write('A_SENDER_FINAL_strict', A_FINAL_STRICT)
# Reported but NOT used: the section 10 list that failed the gate, and two partial-removal
# sensitivity sets. They live in variants/ so the frozen directory is unambiguous.
write_variant('A_PHASE7_S10_16', sorted(a16_panel))
write_variant('A_SENDER_FINAL_noB7', A_FINAL_NO_B7)
write_variant('A_SENDER_FINAL_noB4_noB7', A_FINAL_NO_B4_B7)
for k, v in A_per_module.items():
    write('A_sender_for_%s' % k.split('_', 1)[1], v)
for k, v in B.items():
    write('B_%s' % k.split('_', 1)[1], v)
for k, v in D.items():
    write(k, v)
write('D_nuisance_covariates', D_COVARIATES)
for k, v in E.items():
    write('E_%s' % k.split('_', 1)[1], v)
write('C_ligands', sorted({p for _, p, _, _, _ in tierC_rows if p}))
write('C_receptors', sorted({r for _, _, rp, _, _ in tierC_rows for r in rp}))
if sm_arrest:
    write('A6_SenMayo_arrest_reference', sm_arrest)
write_variant('B_secondary_senescence_v1_curated_ported', B7_V1_CURATED)

with open(OUT + '/_symbol_resolutions.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['asked_symbol', 'resolved_to', 'route'])
    for r in sorted(set(RESOLUTIONS)):
        w.writerow(r)

json.dump({'panel_n': len(PANEL),
           'panel_feature_types': {t: len(v) for t, v in FEAT_BY_TYPE.items()},
           'msigdb_release': MSIGDB_RELEASE, 'msigdb_fetched': MSIGDB_FETCHED,
           'A_phase7_16_on_panel': sorted(a16_panel), 'A_phase7_16_off_panel': a16_off,
           'A_ported_pre': A0,
           'A_final_strict': A_FINAL_STRICT, 'A_final_noB7': A_FINAL_NO_B7,
           'A_final_noB4_noB7': A_FINAL_NO_B4_B7,
           'A_per_module_sizes': {k: len(v) for k, v in A_per_module.items()},
           'B_sizes': {k: len(v) for k, v in B.items()},
           'B_prov': B_prov, 'B6_hallmark_alone_on_panel': B6_ALONE,
           'collisions': collide, 'unmapped_mouse_symbols': UNMAPPED,
           'msigdb': {k: (MS[k][1], MS[k][2], len(MS[k][0]), len(on_panel(MS[k][0]))) for k in MS}},
          open(OUT + '/_test2_summary.json', 'w'), indent=1)
print('\nWritten to %s' % OUT)
