#!/usr/bin/env python3
"""
SASP Spatial Response Kernel -- Deliverable 2 (Gene set package) + Section 8 Test 2.
Target panel: GSE310392 mouse Xenium = Xenium Prime Mouse 5K Pan Tissue & Pathways (5006 genes)
              + GSE310392 custom 100-gene panel = 5106 'Gene Expression' features
              (verified against cell_feature_matrix.h5 of GSM9295284).
Gene set sources:
  MSIGDB_2026.1.Mm : MSigDB MOUSE collections, fetched from gsea-msigdb.org JSON API 2026-08-20.
                     Mouse sets are Alliance Genome Consortium ortholog mappings of the human
                     Hallmark sets (exactSource field states this). NOT hand-mapped.
  CUR              : curated by the biology collaborator from primary senescence literature.
  PANEL            : derived from the panel metadata file's own annotation columns
                     (location / protein_name / cellchat_pathway).
"""
import csv, gzip, json, os, glob
from collections import OrderedDict

# HAZARD FIX 2026-08-27 (Phase 8, D8). This pointed at a per-session /tmp
# scratchpad from a PRIOR session that no longer exists. Re-running the script
# then globbed ZERO MSigDB JSONs and, because the loop below swallowed every
# error with a bare `except`, silently overwrote genesets/*.txt with EMPTY
# Tier B modules and exited 0. The container has been wiped twice, which is
# exactly when someone re-runs a build script.
#
# Now: prefer the archived in-repo pin, and FAIL LOUDLY if no sets are found.
SCRATCH = os.environ.get(
    'MSIGDB_MOUSE_DIR',
    '/workspace/genesets/msigdb_mouse_2026.1.Mm')
OUT = '/workspace/genesets'
os.makedirs(OUT, exist_ok=True)

# HAZARD FIX 2026-08-27 (Phase 8, D8), second order. This script regenerates the
# PRE-C6 mouse gene sets. On 2026-08-27 the C6 re-sourced sets were promoted into
# genesets/ (B7 38->108, A_SENDER_FINAL_strict 25->33). A successful re-run would
# silently revert that promotion and desynchronise the mouse arm from the human
# one, breaking the Section 17 two-arm comparison. Refuse unless explicitly forced.
if os.path.exists(os.path.join(OUT, 'README_C6_mouse_provenance.md')) \
        and os.environ.get('ALLOW_OVERWRITE_C6') != '1':
    raise SystemExit(
        'FATAL: genesets/ currently holds the PROMOTED C6 mouse sets.\n'
        'This script writes the PRE-C6 versions and would silently revert them.\n'
        'Recover the pre-C6 state with `git checkout pre-c6-genesets -- genesets/` '
        'instead, or set ALLOW_OVERWRITE_C6=1 if reverting is genuinely intended.')

# ------------------------------------------------------------------ panel
panel_meta = {}
with open('/workspace/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv') as fh:
    for row in csv.DictReader(fh):
        panel_meta[row['gene_name']] = row
p100 = []
with gzip.open('/workspace/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz','rt') as fh:
    for row in csv.DictReader(fh):
        p100.append(row['Gene'])
GENOTYPING = [g for g in p100 if not g.startswith('ENSMUSG') and ('_WT' in g or '_ALT' in g or '_del_' in g or '_splice_' in g)]
PANEL = set(panel_meta) | set(p100)
assert len(PANEL) == 5106, len(PANEL)

def is_secreted(g):
    """Panel-metadata-driven enforcement of Sec 9's 'no secreted factors' sender constraint."""
    m = panel_meta.get(g)
    if not m: return False
    return 'Secreted' in (m.get('location') or '')

SECRETED_ON_PANEL = {g for g in PANEL if is_secreted(g)}

# ------------------------------------------------------------------ msigdb
MS = {}
_MS_ERRORS = []
_MS_FILES = sorted(glob.glob(SCRATCH + '/*.json'))
if not _MS_FILES:
    raise SystemExit(
        'FATAL: no MSigDB JSON found in %s\n'
        'Set MSIGDB_MOUSE_DIR or restore the archived pin. Refusing to run: '
        'continuing would overwrite genesets/*.txt with EMPTY Tier B modules.' % SCRATCH)
for f in _MS_FILES:
    try:
        d = json.load(open(f))
        k = list(d)[0]
        MS[k] = (set(d[k]['geneSymbols']), d[k]['systematicName'])
    except Exception as _e:
        _MS_ERRORS.append((f, repr(_e)))
        continue
    if False:
        pass

# ------------------------------------------------------------------ TIER A (sender: arrest + damage ONLY)
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

# ------------------------------------------------------------------ TIER B (response)
# B1..B6 come STRAIGHT from MSigDB mouse Hallmark. B7 is curated (see note in README).
TIER_B_SRC = OrderedDict([
 ('B1_tnfa_nfkb_proximal',  ['HALLMARK_TNFA_SIGNALING_VIA_NFKB']),
 ('B2_il6_jak_stat3',       ['HALLMARK_IL6_JAK_STAT3_SIGNALING']),
 ('B3_interferon_response', ['HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE']),
 ('B4_downstream_arrest',   ['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1']),
 ('B5_emt_ecm',             ['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']),
 # B6 RESCUE: HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY alone gives only 17 on-panel genes and
# FAILS the Sec 8 Test 2 >=30 requirement. Documented substitution: union with
# REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES + curated NRF2 target list -> 31 on-panel.
 ('B6_oxidative_stress',    ['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                             'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES']),
])
B6_EXTRA_NRF2 = ('Nqo1 Hmox1 Gclm Gclc Txnrd1 Srxn1 Slc7a11 Gsta3 Gstm1 Gstp1 Ftl1 Fth1 Mgst1 Cbr3 '
                 'Akr1b3 Keap1 Nfe2l2 Maff Mafg Sesn2 Txn1 Prdx1 Gpx4 Cat Sod1 Sod2 Gsr').split()
# B7 secondary senescence: curated readout of the neighbour BECOMING senescent.
B7_CUR = """Cdkn1a Cdkn2a Cdkn2b Trp53 Gadd45a Mdm2 Ccnd1 Glb1 Lmnb1 Hmgb1 Hmgb2 Serpine1 Il1a
 Il1b Il6 Cxcl1 Cxcl2 Ccl2 Ccl20 Mmp3 Mmp12 Timp1 Igfbp3 Igfbp7 Gdf15 Tnfrsf10b Bcl2l1 Trp53i3
 Sirt1 Ezh2 Suv39h1 Chek2 Atm Trp53bp1 Ddb2 Xpc Tgfb1 Thbs1 Nfkbia Junb Fos Egr1 Sod2 Cdkn1b
 Trp53inp1 Zmat3 Bax Eda2r Phlda3 Ccng1""".split()

# ------------------------------------------------------------------ TIER C
TIER_C_PAIRS = [
 ('Il6',    ['Il6ra','Il6st'],            'canonical SASP; drives B2. MOUSE receptor is Il6ra (human IL6R).'),
 ('Cxcl1',  ['Cxcr2'],                    'MOUSE FUNCTIONAL ANALOGUE OF HUMAN CXCL8/IL-8 (KC). Mouse has NO CXCL8 ortholog.'),
 ('Cxcl2',  ['Cxcr2'],                    'MOUSE FUNCTIONAL ANALOGUE OF HUMAN CXCL8/IL-8 (MIP-2). Mouse has NO CXCL8 ortholog.'),
 ('Cxcl5',  ['Cxcr2'],                    'mouse ELR+ CXC chemokine, CXCL8-analogous'),
 ('Ccl2',   ['Ccr2','Ackr3'],             'HIGHEST PRIORITY (Sec 9). Ackr3=CXCR7.'),
 ('Cxcl12', ['Cxcr4','Ackr3','Dpp4'],     'Dpp4 cleaves Cxcl12 -> candidate range-limiting mechanism. Paper reports Cxcl12-Cxcr4 conserved in mouse.'),
 ('Tgfb1',  ['Tgfbr1','Tgfbr2'],          'contact-adjacent; short lambda expected'),
 ('Il1a',   ['Il1r1','Il1rap'],           'IL-1alpha LARGELY MEMBRANE-BOUND -> SHORTEST lambda predicted (Sec 9 internal control anchor)'),
 ('Il1b',   ['Il1r1','Il1rap'],           'secreted; longer lambda than Il1a expected'),
 ('Tnf',    ['Tnfrsf1a','Tnfrsf1b'],      'drives B1'),
 ('Igfbp3', [],                           'SASP Atlas core secreted factor'),
 ('Igfbp7', [],                           'SASP Atlas core secreted factor'),
 ('Gdf15',  [],                           'SASP Atlas core secreted factor'),
 ('Mmp1a',  [],                           'mouse orthologs of human MMP1 are Mmp1a/Mmp1b'),
 ('Mmp1b',  [],                           'mouse orthologs of human MMP1 are Mmp1a/Mmp1b'),
 ('Mmp3',   [],                           'SASP Atlas core secreted factor'),
 ('Timp1',  [],                           'SASP Atlas core secreted factor'),
 ('Thbs1',  ['Cd47','Sdc1'],              'paper-specific: THBS-mediated senescence network (Karpova 2026)'),
]

# ------------------------------------------------------------------ TIER E
TIER_E = OrderedDict([
 ('E1_housekeeping_planA', 'Actb Gapdh Rpl13a Rps18 Tbp Ppia'.split()),
 ('E1b_housekeeping_expanded', """Actb Gapdh Rpl13a Rps18 Tbp Ppia Rplp0 Hprt Pgk1 Sdha Ywhaz
    Rpl11 Tubb5 B2m Eef1a1 Ubc Psmb2 Vcp Rab7 Chmp2a Emc7 Gpi1 Reep5 Snrpd3 Vps29 Cnot4
    Rer1 Tmem59 Actr2 Arpc2 Cfl1 Pfn1 Rack1 Hnrnpk Srp14 Sar1a Atp5f1b Cox4i1 Ndufb8""".split()),
 ('E2_hepatocyte_identity', """Alb Ttr Apoa1 Apoa2 Apob Apoc1 Apoc3 Apoe Trf Serpina1a Serpina1b
    Serpina1c Fga Fgb Fgg Ahsg Hp Hpx Tf Ambp Apoh Orm1 Itih3 Itih4 Vtn Kng1 Gc F2 Cps1
    Ass1 Tat Otc""".split()),
])

# ------------------------------------------------------------------ ZONATION (Sec 11)
# Rebuilt DATA-DRIVEN from GSM9295284 (7250_liver_sham) because the plan's human marker list is
# largely off-panel (Alb, Ass1, Sds, Oat all ABSENT from the Xenium Prime Mouse 5K panel).
# Method: seed axis = mean logNorm(Glul,Cyp2e1,Cyp1a2,Cyp2a5,Cyp27a1)
#                   - mean logNorm(Hal,Cps1,Arg1,Otc,Hamp,Igfbp2),
# computed in the top-60% hepatocyte-marker-scoring cells (n=142,143); Pearson r of every panel
# gene against that axis; retain |r| >= 0.20 and detection rate >= 5%.
# Genes belonging to ANY Tier A or Tier B set are EXCLUDED from the zonation covariate, so that
# conditioning on zonation does not partially condition on the response (over-adjustment).
# Full per-gene table: /workspace/results/zonation_gene_correlations_7250_sham.csv
#
# REPRODUCIBILITY GAP, recorded 2026-08-28 (AUDIT_REPRODUCIBILITY B7).  That per-gene table
# is a gene-set-DEFINING INPUT -- the two lines below turn it into D_zonation_pericentral and
# D_zonation_periportal -- and it has NO COMMITTED PRODUCER.  The method above pins the seed
# axis, the 60% cut (0.6 x 236,905 = 142,143 rows of celltypes_7250..., confirmed), the panel
# (all 5,106 genes) and the two thresholds, but NOT which hepatocyte-marker gene list scored
# the cells for that cut, nor how (sc.tl.score_genes and with what ctrl_size), nor whether
# `detection_rate` is over all cells or over the 142,143.  TIER_E['E2_hepatocyte_identity']
# below is the plausible candidate but is nowhere stated to be it.  A rewrite would land close
# but not exact, and a near-miss silently changes which genes end up in the zonation covariate.
# So the table is treated as a pinned input, like the MSigDB JSON: it is tracked, it is read,
# and it is not regenerated.  If it ever must be rebuilt, diff the resulting GENE LISTS, not
# the r values.
ZONATION_CSV = '/workspace/results/zonation_gene_correlations_7250_sham.csv'
ZON_R_MIN, ZON_DET_MIN = 0.20, 0.05
_zc = list(csv.DictReader(open(ZONATION_CSV)))
_pc = sorted(r['gene'] for r in _zc if float(r['corr_zonation']) >= ZON_R_MIN and float(r['detection_rate']) >= ZON_DET_MIN)
_pp = sorted(r['gene'] for r in _zc if float(r['corr_zonation']) <= -ZON_R_MIN and float(r['detection_rate']) >= ZON_DET_MIN)
ZONATION = OrderedDict([
 ('D_zonation_pericentral', _pc),
 ('D_zonation_periportal',  _pp),
 ('D_landmark_vessel_bileduct', """Pecam1 Cdh5 Vwf Stab2 Lyve1 Kdr Eng Aqp1 Efnb2 Ephb4 Gja5
    Krt19 Krt7 Krt8 Krt18 Epcam Sox9 Spp1 Pkhd1 Onecut1 Hnf1b Cftr Acta2 Myh11 Tagln Des Dcn
    Col1a1 Col1a2 Lrat Pdgfrb Rgs5 Notch3""".split()),
])

# ------------------------------------------------------------------ build + Test 2
def on_panel(gs): return sorted(set(gs) & PANEL)

A_sub = OrderedDict((k, on_panel(v)) for k, v in TIER_A_SUB.items())

B = OrderedDict()
B_prov = OrderedDict()
for mod, srcs in TIER_B_SRC.items():
    g = set()
    for s in srcs:
        if s in MS: g |= MS[s][0]
        else: print('!! MISSING MSIGDB SET', s)
    if mod == 'B6_oxidative_stress': g |= set(B6_EXTRA_NRF2)
    B[mod] = on_panel(g)
    B_prov[mod] = srcs + (['CUR_NRF2_targets'] if mod=='B6_oxidative_stress' else [])
B['B7_secondary_senescence'] = on_panel(B7_CUR)
B_prov['B7_secondary_senescence'] = ['CUR_secondary_senescence (Sec 9 B7)']

B_UNION = set().union(*B.values())
B_UNION_NO_B7 = set().union(*[v for k, v in B.items() if k != 'B7_secondary_senescence'])

print('='*100)
print('PANEL: %d Gene Expression features (%d from Mouse5K pan metadata, %d custom-100, %d genotyping probes)'
      % (len(PANEL), len(panel_meta), len(p100), len(GENOTYPING)))
print('Panel genes annotated location="Secreted": %d' % len(SECRETED_ON_PANEL))
print('='*100)
print('\n### MSigDB mouse sets used (release 2026.1.Mm, fetched 2026-08-20) ###')
for k in sorted(MS):
    n_panel = len(MS[k][0] & PANEL)
    print('  %-58s %s  n=%3d  on-panel=%3d' % (k, MS[k][1], len(MS[k][0]), n_panel))

print('\n### TIER A candidates on-panel (BEFORE disjointness enforcement) ###')
A_all = []
for k, v in A_sub.items():
    miss = sorted(set(TIER_A_SUB[k]) - PANEL)
    print('  %-38s %2d/%2d on-panel   missing: %s' % (k, len(v), len(set(TIER_A_SUB[k])), ','.join(miss) if miss else '-'))
    A_all += v
A0 = sorted(set(A_all))
print('  TIER A TOTAL (union, on-panel, pre-disjointness): %d' % len(A0))
sec_in_A = sorted(set(A0) & SECRETED_ON_PANEL)
print('  Sanity: Tier A members annotated Secreted by the panel (must be ~0): %s' % (sec_in_A or 'NONE'))

print('\n### TIER B modules on-panel ###')
for k, v in B.items():
    print('  %-30s %3d on-panel   [%s]' % (k, len(v), '; '.join(B_prov[k])))

# ---- disjointness enforcement (Sec 8 Test 2): remove from SENDER, never response
collide = {}
for k, v in B.items():
    ov = sorted(set(A0) & set(v))
    if ov: collide[k] = ov
print('\n### A x B COLLISIONS (genes removed from TIER A) ###')
for k, v in collide.items():
    print('  %-30s %3d : %s' % (k, len(v), ' '.join(v)))
A_FINAL_STRICT = sorted(set(A0) - B_UNION)
A_FINAL_NO_B7 = sorted(set(A0) - B_UNION_NO_B7)
print('\n  |A| pre  = %d' % len(A0))
print('  |A| after removing collisions with B1-B6 ONLY (B7 excluded) = %d' % len(A_FINAL_NO_B7))
print('  |A| after removing collisions with B1-B7 (STRICT)           = %d' % len(A_FINAL_STRICT))

print('\n### FULL TIER A x TIER B INTERSECTION MATRIX (post-removal, strict A) ###')
hdr = 'Tier A subset'.ljust(38) + ''.join(('%-14s' % k.split('_')[0]) for k in B)
print(hdr); print('-'*len(hdr))
A_sub_final = OrderedDict((k, sorted(set(v) - B_UNION)) for k, v in A_sub.items())
for k, v in A_sub_final.items():
    row = ('%s (n=%d)' % (k, len(v))).ljust(38)
    row += ''.join(('%-14d' % len(set(v) & set(bv))) for bv in B.values())
    print(row)
row = ('A_FINAL_STRICT (n=%d)' % len(A_FINAL_STRICT)).ljust(38)
row += ''.join(('%-14d' % len(set(A_FINAL_STRICT) & set(bv))) for bv in B.values())
print(row)

print('\n### SECTION 8 TEST 2 VERDICT ###')
c1 = len(A_FINAL_STRICT) >= 15
print('  [%s] len(A) >= 15                : %d' % ('PASS' if c1 else 'FAIL', len(A_FINAL_STRICT)))
allpass_b = True
for k, v in B.items():
    ok = len(v) >= 30
    allpass_b &= ok
    print('  [%s] len(B[%s]) >= 30 : %d' % ('PASS' if ok else 'FAIL', k, len(v)))
inter = sum(len(set(A_FINAL_STRICT) & set(v)) for v in B.values())
print('  [%s] A & B == 0                  : %d' % ('PASS' if inter==0 else 'FAIL', inter))

# ---- SenMayo audit (Sec 9 A6)
if 'SAUL_SEN_MAYO' in MS:
    sm = MS['SAUL_SEN_MAYO'][0]
    sm_p = sm & PANEL
    sm_sec = sm_p & SECRETED_ON_PANEL
    sm_arrest = sorted(sm_p - SECRETED_ON_PANEL - B_UNION_NO_B7)
    print('\n### SenMayo audit (Sec 9 A6) -- SAUL_SEN_MAYO, MSigDB mouse %s ###' % MS['SAUL_SEN_MAYO'][1])
    print('  mouse set size %d ; on-panel %d ; of those annotated Secreted by panel: %d (%.0f%%)'
          % (len(sm), len(sm_p), len(sm_sec), 100*len(sm_sec)/max(len(sm_p),1)))
    print('  SenMayo_arrest = on-panel MINUS secreted MINUS TierB(B1-B6) = %d genes: %s'
          % (len(sm_arrest), ' '.join(sm_arrest)))

# ---- Tier C
print('\n### TIER C ligand-receptor, on-panel ###')
print('%-10s %-8s %-34s %-8s %s' % ('LIGAND','on-panel','RECEPTORS(on-panel/total)','','NOTE'))
tierC_rows = []
for lig, recs, note in TIER_C_PAIRS:
    lp = lig in PANEL
    rp = [r for r in recs if r in PANEL]
    rm = [r for r in recs if r not in PANEL]
    tierC_rows.append((lig, lp, rp, rm, note))
    print('%-10s %-8s %-34s %s' % (lig, 'YES' if lp else 'no',
          '%d/%d [%s]%s' % (len(rp), len(recs), ','.join(rp), (' MISSING:'+','.join(rm)) if rm else ''), note))
print('  NOTE: Cxcl8/Il8 -- MOUSE HAS NO ORTHOLOG. Cxcr1 also absent from the mouse genome as a')
print('        clean ortholog; mouse ELR+ CXC signalling runs Cxcl1/Cxcl2/Cxcl5 -> Cxcr2.')

# ---- panel cellchat pathways relevant to Tier C
cc = {}
for g, m in panel_meta.items():
    v = m.get('cellchat_pathway') or ''
    for t in v.split(';'):
        t = t.strip()
        if t: cc.setdefault(t, []).append(g)
print('\n### PANEL-annotated cellchat pathways relevant to SASP (source=PANEL) ###')
for p in ['CXCL','CCL','IL6','IL1','TNF','TGFb','THBS','SPP1','IGF','GDF','ANNEXIN','COMPLEMENT']:
    if p in cc: print('  %-12s %2d : %s' % (p, len(cc[p]), ' '.join(sorted(cc[p]))))

# ---- Tier E
print('\n### TIER E negative controls, on-panel ###')
for k, v in TIER_E.items():
    op = on_panel(v)
    print('  %-28s %2d/%2d : %s' % (k, len(op), len(set(v)), ' '.join(op)))
    print('     MISSING: %s' % ' '.join(sorted(set(v)-PANEL)))

# ---- zonation
# over-adjustment guard: zonation covariate must not share genes with Tier A or Tier B
_AB = set(A0) | B_UNION
ZON_DROPPED = {k: sorted(set(v) & _AB) for k, v in ZONATION.items()}
ZONATION = OrderedDict((k, [g for g in v if g not in _AB]) for k, v in ZONATION.items())
print('\n### TIER D / ZONATION (Sec 11), on-panel ###')
print('  (genes dropped because they also appear in Tier A or Tier B -- over-adjustment guard)')
for k, v in ZON_DROPPED.items():
    print('    %-30s dropped %2d: %s' % (k, len(v), ' '.join(v) or '-'))
for k, v in ZONATION.items():
    op = on_panel(v)
    print('  %-30s %2d/%2d' % (k, len(op), len(set(v))))
    print('     ON-PANEL: %s' % ' '.join(op))
    print('     MISSING : %s' % ' '.join(sorted(set(v)-PANEL)))

# ---- write files
def write(name, genes):
    with open(os.path.join(OUT, name + '.txt'), 'w') as fh:
        fh.write('\n'.join(genes) + '\n')

for k, v in A_sub_final.items(): write('A_%s' % k.split('_',1)[1], v)
write('A_SENDER_FINAL_strict', A_FINAL_STRICT)
write('A_SENDER_FINAL_noB7', A_FINAL_NO_B7)
for k, v in B.items(): write('B_%s' % k.split('_',1)[1], v)
for k, v in TIER_E.items(): write('E_%s' % k.split('_',1)[1], on_panel(v))
for k, v in ZONATION.items(): write('D_%s' % k.split('_',1)[1], on_panel(v))
write('C_ligands', [l for l,p,_,_,_ in tierC_rows if p])
write('C_receptors', sorted({r for _,_,rp,_,_ in tierC_rows for r in rp}))
if 'SAUL_SEN_MAYO' in MS: write('A6_SenMayo_arrest_reference', sm_arrest)

json.dump({'panel_n':len(PANEL),
           'A_final_strict':A_FINAL_STRICT,'A_final_noB7':A_FINAL_NO_B7,
           'B_sizes':{k:len(v) for k,v in B.items()},
           'collisions':collide,
           'msigdb':{k:(MS[k][1],len(MS[k][0]),len(MS[k][0]&PANEL)) for k in MS}},
          open(OUT+'/_test2_summary.json','w'), indent=1)
print('\nWritten to %s' % OUT)
