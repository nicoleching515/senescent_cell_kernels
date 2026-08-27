#!/usr/bin/env python3
"""
Phase 7 section 11 -- the disjointness gate, HUMAN arm (H1 = GSE326743).

Runs the three assertions of section 11 verbatim:

    A   = tier_A & panel
    B_k = tier_B_module_k & panel        # for each of the 7
    assert len(A) >= 15
    assert all(len(B_k) >= 30 for k in modules)
    assert all(len(A & B_k) == 0 for k in modules)

and, per section 11, when A & B_k is non-empty the overlapping genes are removed FROM TIER A
(never from the response module) and len(A) >= 15 is re-checked. Every removed gene is printed
with the module that removed it.

The gate is run against BOTH candidate Tier A definitions, because they give different verdicts:
  * A_PHASE7_S10_16 -- the literal 16-symbol list printed in Phase 7 section 10.
  * A_ported        -- the mouse arm's five Tier A subsets ported to human symbols
                       (build_genesets_human.py), i.e. the like-for-like port of the mouse Tier A.

Also emits, per section 11 ("Print the full intersection matrix in the Methods"):
  results/phase7_jobA/intersection_matrix_human.csv   Tier A rows x Tier B columns, counts + genes
  results/phase7_jobA/intersection_matrix_BxB.csv     Tier B x Tier B, module cross-talk
  results/phase7_jobA/onpanel_counts_by_tier.csv      on-panel counts for every tier file
  results/phase7_jobA/gate_result_human.json          machine-readable verdict

and the CoreScence circularity check. The mouse-arm reference it is compared against is NOT
typed in here: it is derived from files by `code/corescence_circularity.py` (an earlier version
of this script asserted `24/35 = 69%`, whose denominator is reproducible under no mapping
convention -- see reports/AUDIT_PHASE8_FACTCHECK.md M1). CoreScence is a HUMAN gene set, so on
H1 it is evaluated natively with no ortholog remapping; the mouse arm needs the MGI map.

Run:  python3 /workspace/code/gate_disjointness_human.py
Exit status 0 = every assertion passed for at least one Tier A definition; 1 = the primary
(section 10) definition failed. Nothing is written outside results/phase7_jobA/.
"""
import csv, json, os, sys
from collections import OrderedDict

sys.path.insert(0, '/workspace/code')
import corescence_circularity as CC

GS    = '/workspace/genesets'
HUM   = GS + '/human'
OUT   = '/workspace/results/phase7_jobA'
CORE  = '/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv'
os.makedirs(OUT, exist_ok=True)

PANEL = {r['gene_name'] for r in csv.DictReader(open(GS + '/h1_candidate/GSE326743_gene_panel_5093.csv'))}

def gl(name):
    return [l.strip() for l in open(os.path.join(HUM, name + '.txt')) if l.strip()]

def glv(name):
    """reported-but-not-used sets, archived under genesets/human/variants/"""
    return [l.strip() for l in open(os.path.join(HUM, 'variants', name + '.txt')) if l.strip()]

MODULES = ['tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
           'emt_ecm', 'oxidative_stress', 'secondary_senescence']
BLAB    = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
B = OrderedDict((m, set(gl('B_' + m)) & PANEL) for m in MODULES)

A_SUBSETS = OrderedDict((k, set(gl('A_' + k)) & PANEL) for k in
                        ['core_arrest', 'proliferation_down', 'nuclear_chromatin',
                         'dna_damage_response', 'senescence_curated_nonsecreted'])
SUMMARY = json.load(open(HUM + '/_test2_summary.json'))

A_DEFS = OrderedDict([
    ('A_PHASE7_S10_16', set(glv('A_PHASE7_S10_16')) & PANEL),
    ('A_ported',        set(SUMMARY['A_ported_pre']) & PANEL),
])

BAR = '=' * 100
print(BAR)
print('PHASE 7 SECTION 11 DISJOINTNESS GATE -- HUMAN ARM H1 (GSE326743, %d Gene Expression features)'
      % len(PANEL))
print('MSigDB release %s fetched %s' % (SUMMARY['msigdb_release'], SUMMARY['msigdb_fetched']))
print(BAR)

# ---------------------------------------------------------------- the gate
results = OrderedDict()
for aname, A_raw in A_DEFS.items():
    print('\n### Tier A definition: %s   (|A & panel| = %d) ###' % (aname, len(A_raw)))
    c1 = len(A_raw) >= 15
    print('  [%s] assert len(A & panel) >= 15            : %d' % ('PASS' if c1 else 'FAIL', len(A_raw)))

    c2, b_rows = True, []
    for m, lab in zip(MODULES, BLAB):
        ok = len(B[m]) >= 30
        c2 &= ok
        b_rows.append((lab, m, len(B[m]), ok))
        print('  [%s] assert len(B[%s] & panel) >= 30  %-22s : %d'
              % ('PASS' if ok else 'FAIL', lab, m, len(B[m])))

    removed = OrderedDict()
    for m, lab in zip(MODULES, BLAB):
        ov = sorted(A_raw & B[m])
        if ov:
            removed[m] = ov
    c3 = not removed
    print('  [%s] assert len(A & B_k) == 0 for all k     : %d overlapping gene memberships'
          % ('PASS' if c3 else 'FAIL', sum(len(v) for v in removed.values())))

    A_final = set(A_raw)
    if removed:
        print('\n  Section 11 remedy: remove the overlapping genes FROM TIER A, then re-check >= 15.')
        for m, ov in removed.items():
            print('    removed by %-22s (%2d): %s' % (m, len(ov), ' '.join(ov)))
        A_final = A_raw - set().union(*removed.values())
        allremoved = sorted(A_raw - A_final)
        print('    UNIQUE genes removed from Tier A (%d): %s' % (len(allremoved), ' '.join(allremoved)))
    c1b = len(A_final) >= 15
    print('  [%s] re-check len(A) >= 15 after removal    : %d' % ('PASS' if c1b else 'FAIL', len(A_final)))
    print('    surviving Tier A: %s' % (' '.join(sorted(A_final)) or '(EMPTY)'))
    verdict = 'PASS' if (c1 and c2 and c1b) else 'FAIL'
    print('  ==> GATE VERDICT for %s: %s' % (aname, verdict))

    results[aname] = dict(n_A_raw=len(A_raw), assert_A_ge15=c1, assert_B_ge30=c2,
                          assert_disjoint=c3, removed={m: v for m, v in removed.items()},
                          n_removed=len(A_raw - A_final), n_A_final=len(A_final),
                          A_final=sorted(A_final), recheck_A_ge15=c1b, verdict=verdict,
                          B_sizes={m: len(B[m]) for m in MODULES})

# ------------------------------------------------- per-module sender sets: PRE-REGISTERED
# PI decision 2026-08-27: A_SENDER_FINAL_strict is PRIMARY; the per-module sets are the
# pre-registered SENSITIVITY analysis, so they are gated here too, not treated as informal.
# Disjointness is required between the sender score and the ONE response module fitted.
print('\n### PRE-REGISTERED SENSITIVITY: per-module sender sets ###')
print('%-24s %7s %6s %9s   %s' % ('module', '|A_mod|', '>=15?', 'A_mod&B_k', 'canonical markers retained'))
CANON = ['CDKN1A', 'CDKN2A', 'CDKN2B', 'TP53', 'MKI67', 'LMNB1', 'TOP2A', 'PCNA', 'MDM2',
         'ATM', 'GLB1', 'HMGB2']
per_mod, pm_pass = OrderedDict(), True
for m in MODULES:
    s_ = set(gl('A_sender_for_' + m)) & PANEL
    per_mod[m] = s_
    ov = len(s_ & B[m])
    ok = len(s_) >= 15 and ov == 0
    pm_pass &= ok
    print('%-24s %7d %6s %9d   %s' % (m, len(s_), 'PASS' if len(s_) >= 15 else 'FAIL', ov,
                                      ' '.join(g for g in CANON if g in s_) or '-'))
print('  [%s] per-module gate: all 7 have |A_mod| >= 15 AND A_mod n B_k == 0'
      % ('PASS' if pm_pass else 'FAIL'))

# ---------------------------------------------------------------- intersection matrices
rows = []
A_ROWS = OrderedDict()
for k, v in A_SUBSETS.items():
    A_ROWS['A_' + k] = v
A_ROWS['A_PHASE7_S10_16 (variant, FAILED)'] = A_DEFS['A_PHASE7_S10_16']
A_ROWS['A_ported_pre_disjointness'] = A_DEFS['A_ported']
A_ROWS['A_SENDER_FINAL_noB4_noB7 (variant)'] = set(glv('A_SENDER_FINAL_noB4_noB7')) & PANEL
A_ROWS['A_SENDER_FINAL_noB7 (variant)'] = set(glv('A_SENDER_FINAL_noB7')) & PANEL
for _m in MODULES:
    A_ROWS['A_sender_for_%s (sensitivity)' % _m] = set(gl('A_sender_for_' + _m)) & PANEL
A_ROWS['A_SENDER_FINAL_strict (PRIMARY)'] = set(gl('A_SENDER_FINAL_strict')) & PANEL

print('\n### FULL TIER A x TIER B INTERSECTION MATRIX (on-panel gene counts) ###')
hdr = 'Tier A set'.ljust(34) + 'n'.ljust(6) + ''.join(l.ljust(6) for l in BLAB) + 'total'
print(hdr)
print('-' * len(hdr))
for aname, aset in A_ROWS.items():
    line = aname.ljust(34) + str(len(aset)).ljust(6)
    tot = 0
    r = {'tier_A_set': aname, 'n_A_on_panel': len(aset)}
    for m, lab in zip(MODULES, BLAB):
        n = len(aset & B[m])
        tot += n
        line += str(n).ljust(6)
        r['%s_%s' % (lab, m)] = n
        r['%s_genes' % lab] = ';'.join(sorted(aset & B[m]))
    r['total_memberships'] = tot
    print(line + str(tot))
    rows.append(r)
cols = (['tier_A_set', 'n_A_on_panel'] +
        ['%s_%s' % (l, m) for l, m in zip(BLAB, MODULES)] +
        ['total_memberships'] + ['%s_genes' % l for l in BLAB])
with open(OUT + '/intersection_matrix_human.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)

print('\n### TIER B x TIER B INTERSECTION MATRIX (module cross-talk, on-panel) ###')
hdr = 'module'.ljust(24) + 'n'.ljust(6) + ''.join(l.ljust(6) for l in BLAB)
print(hdr)
print('-' * len(hdr))
bb = []
for m, lab in zip(MODULES, BLAB):
    line = ('%s %s' % (lab, m)).ljust(24) + str(len(B[m])).ljust(6)
    r = {'module': '%s_%s' % (lab, m), 'n_on_panel': len(B[m])}
    for m2, lab2 in zip(MODULES, BLAB):
        n = len(B[m] & B[m2])
        line += str(n).ljust(6)
        r[lab2] = n
    print(line)
    bb.append(r)
with open(OUT + '/intersection_matrix_BxB.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['module', 'n_on_panel'] + BLAB)
    w.writeheader()
    for r in bb:
        w.writerow(r)

# ---------------------------------------------------------------- on-panel counts per tier
# ---------------------------------------------------------------- module margins + pre-C6 BxB
print('\n### MODULE MARGIN OVER THE >=30 FLOOR (human) ###')
marg = []
for m, lab in zip(MODULES, BLAB):
    mg = len(B[m]) - 30
    marg.append(dict(arm='human', module='%s_%s' % (lab, m), n=len(B[m]), margin=mg))
    print('  %-30s n=%3d margin=%+3d%s' % (m, len(B[m]), mg, '  <-- NO MARGIN' if mg <= 1 else ''))
with open(OUT + '/module_margins_human.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['arm', 'module', 'n', 'margin'])
    w.writeheader()
    for r in marg:
        w.writerow(r)

B_PRE = dict(B)
B_PRE['secondary_senescence'] = set(glv('B_secondary_senescence_v1_curated_ported')) & PANEL
bbp = []
for m, lab in zip(MODULES, BLAB):
    r = {'module': '%s_%s' % (lab, m), 'n_on_panel': len(B_PRE[m])}
    for m2, lab2 in zip(MODULES, BLAB):
        r[lab2] = len(B_PRE[m] & B_PRE[m2])
    bbp.append(r)
with open(OUT + '/intersection_matrix_BxB_preC6.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['module', 'n_on_panel'] + BLAB)
    w.writeheader()
    for r in bbp:
        w.writerow(r)

print('\n### ON-PANEL COUNTS PER TIER FILE ###')
tier_rows = []
for f in sorted(os.listdir(HUM)):
    if not f.endswith('.txt'):
        continue
    name = f[:-4]
    g = gl(name)
    if name == 'D_nuisance_covariates':
        tier_rows.append(dict(file=name, tier='D', n_entries=len(g), n_on_panel='n/a',
                              note='covariate names, not genes'))
        continue
    if name.startswith('E_negative') or name.startswith('E_genomic'):
        tier_rows.append(dict(file=name, tier='E', n_entries=len(g), n_on_panel='n/a',
                              note='control features, not Gene Expression'))
        continue
    tier_rows.append(dict(file=name, tier=name[0], n_entries=len(g),
                          n_on_panel=len(set(g) & PANEL), note=''))
for r in tier_rows:
    print('  %-36s %-4s entries=%-5s on-panel=%-5s %s'
          % (r['file'], r['tier'], r['n_entries'], r['n_on_panel'], r['note']))
with open(OUT + '/onpanel_counts_by_tier.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['file', 'tier', 'n_entries', 'n_on_panel', 'note'])
    w.writeheader()
    for r in tier_rows:
        w.writerow(r)

# ---------------------------------------------------------------- CoreScence circularity
print('\n### CoreScence circularity ###')
try:
    _mouse_cc = CC.derive()
    _mouse_pre = CC.reference_string(_mouse_cc, 'pre_C6')
    _mouse_c6 = CC.reference_string(_mouse_cc, 'C6_promoted')
    print('  mouse arm, DERIVED (code/corescence_circularity.py, %s mapping): pre-C6 %s ; '
          'C6-promoted %s' % (_mouse_cc['convention_to_cite'], _mouse_pre, _mouse_c6))
except Exception as _e:                       # DeepScence uninstalled, or the mouse panel absent
    _mouse_pre = _mouse_c6 = None
    print('  !! mouse-arm reference NOT derivable here (%s: %s) -- it is therefore NOT reported. '
          'Never substitute a remembered number.' % (type(_e).__name__, _e))
print('  CITE THE FROZEN NUMBER. Re-sourcing B7 raised it -- see reports/BIO_PHASE8_FREEZE.md sec 3.')
core_res = None
if os.path.exists(CORE):
    cr = list(csv.DictReader(open(CORE)))
    occ = [r for r in cr if r['occurrence'] and float(r['occurrence']) >= 5]
    core_genes = sorted({r['gene_symbol'] for r in occ})
    core_panel = sorted(set(core_genes) & PANEL)
    # B7 was re-sourced on 2026-08-27 (35 -> 116 genes), which moves this number. Both are
    # computed here, from files, so the report can state which one to cite.
    B_OLD = dict(B)
    B_OLD['secondary_senescence'] = set(glv('B_secondary_senescence_v1_curated_ported')) & PANEL
    out = {}
    for tag, Bx in (('FROZEN (re-sourced B7, n=%d)' % len(B['secondary_senescence']), B),
                    ('SUPERSEDED (curated v1 B7, n=%d)' % len(B_OLD['secondary_senescence']), B_OLD)):
        inB = sorted({g for g in core_panel if any(g in Bx[m] for m in MODULES)})
        out[tag] = inB
        print('  %-42s %d/%d = %.0f%%' % (tag, len(inB), len(core_panel),
                                          100 * len(inB) / len(core_panel)))
    print('  CoreScence v2, occurrence >= 5 : %d genes (HUMAN symbols, NO remapping needed on H1)'
          % len(core_genes))
    print('  on the GSE326743 panel         : %d' % len(core_panel))
    frozen_tag = [t for t in out if t.startswith('FROZEN')][0]
    inB = out[frozen_tag]
    print('  circular genes under the FROZEN Tier B: %s' % ' '.join(inB))
    print('  per module (FROZEN):')
    for m, lab in zip(MODULES, BLAB):
        ov = sorted(set(core_panel) & B[m])
        print('    %-24s %2d / %d = %.2f' % (m, len(ov), len(core_panel), len(ov) / len(core_panel)))
    core_res = dict(source=CORE, n_occ_ge5=len(core_genes), n_on_panel=len(core_panel),
                    frozen_n_in_any_B=len(inB),
                    frozen_frac=round(len(inB) / len(core_panel), 4),
                    superseded_n_in_any_B=len(out[[t for t in out if t.startswith('SUPER')][0]]),
                    superseded_frac=round(len(out[[t for t in out if t.startswith('SUPER')][0]])
                                          / len(core_panel), 4),
                    genes_in_any_B_frozen=inB,
                    per_module_frozen={m: sorted(set(core_panel) & B[m]) for m in MODULES},
                    mouse_arm_reference=_mouse_pre,
                    mouse_arm_reference_C6=_mouse_c6,
                    mouse_arm_source='code/corescence_circularity.py -> '
                                     'results/phase7_jobA/corescence_circularity_mouse.json',
                    cite='frozen')
else:
    print('  !! %s NOT FOUND -- DeepScence is not installed, so this number cannot be computed '
          'from a file on disk and is NOT reported.' % CORE)

# ---------------------------------------------------------------- B7 overlap with the Tier A caller
print('\n### Phase 7 section 23 / C6 -- B7 secondary_senescence overlap with the Tier A caller ###')
for aname, A_raw in A_DEFS.items():
    ov = sorted(A_raw & B['secondary_senescence'])
    print('  %-18s |B7|=%d  shares %d genes with %s: %s'
          % (aname, len(B['secondary_senescence']), len(ov), aname, ' '.join(ov) or '-'))
print('  Mouse arm reference (section 23 C6): B7 shared 14 of 38 genes with the Tier A caller.')
print('  C6 ADOPTED 2026-08-27: B7 is now re-sourced from SAUL_SEN_MAYO + REACTOME_SASP minus the')
print('  Tier A caller pool, so the overlap is 0 by construction. The superseded curated v1 is')
print('  archived at genesets/human/variants/B_secondary_senescence_v1_curated_ported.txt and is')
print('  reported alongside, as section 23 / C6 requires.')

json.dump(dict(panel_n=len(PANEL), msigdb_release=SUMMARY['msigdb_release'],
               per_module_gate_pass=pm_pass,
               msigdb_fetched=SUMMARY['msigdb_fetched'],
               gate=results, per_module_sender_sizes={m: len(v) for m, v in per_mod.items()},
               B_sizes={m: len(B[m]) for m in MODULES}, corescence=core_res),
          open(OUT + '/gate_result_human.json', 'w'), indent=1)
print('\nWritten to %s' % OUT)

# Exit status reflects the FROZEN configuration: primary Tier A + the pre-registered
# per-module sensitivity. The section 10 sixteen is reported, not used, so it no longer
# determines the exit code.
FROZEN_OK = results['A_ported']['verdict'] == 'PASS' and pm_pass
print('\n### FROZEN CONFIGURATION VERDICT: %s ###' % ('PASS' if FROZEN_OK else 'FAIL'))
sys.exit(0 if FROZEN_OK else 1)
