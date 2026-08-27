#!/usr/bin/env python3
"""
Phase 7 section 23 / C6 -- rebuild Tier B module B7 `secondary_senescence` (HUMAN arm).

C6 asks for three things: rebuild B7 without the genes it shares with the Tier A caller, report
both versions, and replace the meeting-abstract citation. Job A promoted this from optional: B7
causes 11 of the 23 A-in-B collisions against the section 10 Tier A and shares 23 of its 35 genes
with the ported Tier A.

The Tier A choice (strict-21 vs per-module) is the PI's, so EVERY variant is built against BOTH
candidate callers and both are reported. The original file is left untouched.

Variants written to genesets/human/:
  B_secondary_senescence                     v1, unchanged (mouse arm's curated list, ported)
  B_secondary_senescence_C6_minusA_s10       v1 minus A_PHASE7_S10_16
  B_secondary_senescence_C6_minusA_ported    v1 minus the ported Tier A
  B_secondary_senescence_C6_sourced_s10      MSigDB-sourced rebuild minus A_PHASE7_S10_16
  B_secondary_senescence_C6_sourced_ported   MSigDB-sourced rebuild minus the ported Tier A

CITATION. v1 rests on Neretti 2024 (`neretti2024dissecting`), a ~250-word GSA meeting abstract --
audit finding B8, "DO NOT CITE AS A PAPER". The sourced rebuild replaces it with sets that carry
peer-reviewed PMIDs recorded in the archived MSigDB JSON:
  SAUL_SEN_MAYO                                        PMID 35974106  (saul2022senmayo, in references.bib)
  REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP (Reactome pathway R-HSA-2559582)
For the primary-versus-secondary DISTINCTION itself the audit's own recommended replacements are
already in references.bib: `martin2023modelling` and Acosta et al. 2013. This script does not
invent a citation; it points at entries that exist in the repo.

Run: python3 /workspace/code/rebuild_b7_secondary_senescence.py
"""
import csv, json, os, sys
sys.path.insert(0, '/workspace/code')
import human_symbols as HS

GS, HUM = '/workspace/genesets', '/workspace/genesets/human'
MS = GS + '/msigdb_human_2026.1.Hs'
OUT = '/workspace/results/phase7_jobA'
PANEL = HS.PANEL

def gl(n):
    return {l.strip() for l in open(os.path.join(HUM, n + '.txt')) if l.strip()}

def msig(n):
    d = json.load(open('%s/%s.json' % (MS, n)))
    k = list(d)[0]
    return set(d[k]['geneSymbols']), d[k].get('pmid') or '(none in MSigDB record)'

MODULES = ['tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
           'emt_ecm', 'oxidative_stress', 'secondary_senescence']
B = {m: gl('B_' + m) for m in MODULES}
# v1 is the SUPERSEDED curated module, archived by build_genesets_human.py. As of the PI decision
# of 2026-08-27 the FROZEN B_secondary_senescence is the re-sourced version; this script exists to
# satisfy section 23 / C6's requirement that BOTH be reported.
V1 = {l.strip() for l in open(HUM + '/variants/B_secondary_senescence_v1_curated_ported.txt') if l.strip()}
SUMMARY = json.load(open(HUM + '/_test2_summary.json'))
CALLERS = {'A_PHASE7_S10_16': {l.strip() for l in open(HUM + '/variants/A_PHASE7_S10_16.txt') if l.strip()},
           'A_ported': set(SUMMARY['A_ported_pre']) & PANEL}

B7v1 = V1
print('=' * 100)
print('C6 -- B7 secondary_senescence rebuild, human arm')
print('=' * 100)
print('v1 (unchanged, ported mouse curated list): %d on-panel genes' % len(B7v1))

senmayo, pm_sm = msig('SAUL_SEN_MAYO')
sasp, pm_sasp = msig('REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP')
print('  SAUL_SEN_MAYO            PMID %-10s %d genes, %d on-panel' % (pm_sm, len(senmayo), len(HS.on_panel(senmayo))))
print('  REACTOME_SASP            PMID %-10s %d genes, %d on-panel' % (pm_sasp, len(sasp), len(HS.on_panel(sasp))))
B7src = set(HS.on_panel(senmayo | sasp))
print('  sourced union on-panel   %d genes' % len(B7src))

rows = []
variants = {}
for cname, A in CALLERS.items():
    for vname, base in (('minusA', B7v1), ('sourced', B7src)):
        keep = sorted(base - A)
        tag = 'B_secondary_senescence_C6_%s_%s' % (vname, 's10' if cname.startswith('A_PHASE7') else 'ported')
        variants[tag] = (keep, cname, base)
        rows.append(dict(variant=tag, base='v1_curated' if vname == 'minusA' else 'SenMayo+ReactomeSASP',
                         caller=cname, n_base=len(base), n_shared_with_caller=len(base & A),
                         n_after=len(keep), passes_ge30='yes' if len(keep) >= 30 else 'NO'))

print('\n### B7 variants ###')
print('%-44s %-22s %-16s %5s %7s %6s %s' % ('variant', 'base', 'caller', 'base', 'shared', 'after', '>=30?'))
for r in rows:
    print('%-44s %-22s %-16s %5d %7d %6d %s'
          % (r['variant'], r['base'], r['caller'], r['n_base'], r['n_shared_with_caller'],
             r['n_after'], 'PASS' if r['passes_ge30'] == 'yes' else '**FAIL**'))

for tag, (keep, _, _) in variants.items():
    with open(os.path.join(HUM, 'variants', tag + '.txt'), 'w') as fh:
        fh.write('\n'.join(keep) + '\n')

# ---- what changes in the section 11 gate if B7 is swapped -------------------------------
print('\n### Section 11 gate under each B7 variant ###')
print('%-44s %-16s %6s %8s %10s %8s' % ('B7 used', 'Tier A', '|B7|', '|B7|>=30', 'A-in-B7', 'A_strict'))
gate = []
for tag in ['B_secondary_senescence (FROZEN, re-sourced)',
            'B_secondary_senescence_v1_curated_ported'] + list(variants):
    b7 = (B['secondary_senescence'] if tag.startswith('B_secondary_senescence (') else
          V1 if tag.endswith('v1_curated_ported') else
          {l.strip() for l in open(HUM + '/variants/' + tag + '.txt') if l.strip()})
    for cname, A in CALLERS.items():
        Bx = dict(B)
        Bx['secondary_senescence'] = b7
        union = set().union(*Bx.values())
        a_fin = A - union
        ov = len(A & b7)
        gate.append(dict(b7=tag, caller=cname, n_b7=len(b7), b7_ge30=len(b7) >= 30,
                         a_in_b7=ov, n_A_final=len(a_fin),
                         verdict='PASS' if (len(A) >= 15 and len(b7) >= 30 and len(a_fin) >= 15) else 'FAIL'))
        print('%-44s %-16s %6d %8s %10d %8d  %s'
              % (tag, cname, len(b7), 'PASS' if len(b7) >= 30 else 'FAIL', ov, len(a_fin), gate[-1]['verdict']))

with open(OUT + '/b7_c6_rebuild.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['b7', 'caller', 'n_b7', 'b7_ge30', 'a_in_b7', 'n_A_final', 'verdict'])
    w.writeheader()
    for r in gate:
        w.writerow(r)
json.dump(dict(variants=rows, gate=gate,
               citations=dict(v1='neretti2024dissecting (GSA meeting ABSTRACT -- audit B8, do not cite)',
                              senmayo_pmid=pm_sm, reactome_sasp_pmid=pm_sasp,
                              distinction_replacements=['martin2023modelling', 'acosta2013 (references.bib)'])),
          open(OUT + '/b7_c6_rebuild.json', 'w'), indent=1)
print('\nWritten to %s/b7_c6_rebuild.{csv,json} and %d comparison files in %s/variants' % (OUT, len(variants), HUM))
