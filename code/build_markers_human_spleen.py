#!/usr/bin/env python3
"""
Phase 7 Job A follow-on, task 1 -- HUMAN SPLEEN cell-type markers for the GSE326743 Xenium panel.

Replaces Phase 7 section 14 step 1's AT1/AT2/club/ciliated/alveolar-macrophage list, which was
written for a lung arm. H1 is spleen (reports/PHASE7_H1_SCREEN.md).

MARKER SOURCE: CellMarker 2.0 human, pinned at genesets/cellmarker_pin/ (60,877 rows, md5
recorded). NOT curated from memory. Every gene that survives is traceable to the PMIDs in that
table; the per-gene evidence table is written to genesets/human/markers_spleen_evidence.csv.

Selection, in order:
  1. CANDIDATES  rows whose cell_name matches the label's selectors. Spleen-tissue rows are used
                 alone where the label has enough of them; labels that do not clear MIN_SPLEEN_ROWS
                 fall back to all tissues and are FLAGGED, because a marker validated in another
                 tissue is weaker evidence here.
  2. EVIDENCE    a gene needs >= MIN_PMID distinct PMIDs behind it for that label.
  3. PANEL       the gene must resolve onto the 5,093-gene panel (human_symbols.resolve, so panel
                 legacy symbols such as H2AFX are matched, not missed).
  4. SPECIFICITY a gene claimed by more than MAX_LABELS of the candidate labels is dropped as
                 promiscuous. This is the mechanical version of annotate_pipeline's
                 DROP_NONSPECIFIC rule, applied to genes rather than to whole labels.
  5. GATE        a label needs >= MIN_MARKERS surviving genes (annotate_pipeline's own threshold)
                 or it is removed from the label set and reported as unresolvable on this panel.

WHAT THIS SCRIPT CANNOT DO: the mouse arm's DROP_NONSPECIFIC entries were justified by MEASURED
expression in GSM9295284. H1 is behind the section 15 freeze, so no H1 expression value may be
read. Every judgement here is therefore panel-membership and literature-evidence only, and the
label set must be re-gated against measured expression once the freeze is committed.

Writes:  code/markers_human_spleen.py            (static MARKERS/MERGE/DROP dicts, same format as
                                                  markers_mouse_liver.py -- consumed unchanged)
         genesets/human/markers_spleen_evidence.csv
Run:     python3 /workspace/code/build_markers_human_spleen.py
"""
import csv, gzip, os, sys, collections
sys.path.insert(0, '/workspace/code')
import human_symbols as HS

CM   = '/workspace/genesets/cellmarker_pin/CellMarker2.0_Human.2026-08-27.csv.gz'
OUTP = '/workspace/code/markers_human_spleen.py'
OUTE = '/workspace/genesets/human/markers_spleen_evidence.csv'

MIN_PMID_SPLEEN, MIN_PMID_FALLBACK = 1, 2
MAX_LABELS, MIN_MARKERS, MIN_SPLEEN_ROWS, MAX_PER_LABEL = 3, 4, 8, 15

# label -> list of cell_name selectors (lower-cased substring match on cell_name)
SELECTORS = collections.OrderedDict([
 ('Red pulp macrophages',        ['macrophage']),
 ('Monocytes',                   ['monocyte']),
 ('cDC1',                        ['conventional dendritic cell 1', 'cdc1', 'clec9a+ dendritic']),
 ('cDC2',                        ['conventional dendritic cell 2', 'cdc2', 'clec10a+ dendritic']),
 ('pDC',                         ['plasmacytoid dendritic']),
 ('Follicular B cells',          ['follicular b cell', 'naive b cell', 'mantle b cell']),
 ('Marginal zone B cells',       ['marginal zone b cell', 'marginal zone(mz)-like b cell']),
 ('Germinal centre B cells',     ['germinal center b cell', 'germinal centre b cell']),
 ('Plasma cells',                ['plasma cell', 'plasmablast']),
 ('CD4 T cells',                 ['cd4+ t cell', 'naive cd4 t cell', 'central memory cd4+ t cell',
                                  'cd4+ t helper cell', 't follicular helper']),
 ('CD8 T cells',                 ['cd8+ t cell', 'naive cd8 t cell', 'central memory cd8+ t cell',
                                  'effector memory cd8+ t cell', 'cytotoxic cd8+ t cell']),
 ('NK cells',                    ['natural killer cell', 'nk cell']),
 ('Follicular dendritic cells',  ['follicular dendritic cell']),
 ('Fibroblastic reticular cells',['fibroblastic reticular cell', 'reticular fibroblast']),
 ('Sinusoidal endothelium',      ['sinusoidal endothelial cell']),
 ('Endothelial cells',           ['endothelial cell']),
 ('Lymphatic endothelium',       ['lymphatic endothelial cell']),
 ('Smooth muscle / capsule',     ['smooth muscle cell', 'vascular smooth muscle cell(vsmc)',
                                  'capsular fibroblast']),
 ('Pericytes',                   ['pericyte']),
 ('Fibroblasts',                 ['fibroblast']),
 ('Erythroid cells',             ['erythroid cell', 'erythroid lineage cell', 'erythroblast',
                                  'red blood cell (erythrocyte)']),
 ('Megakaryocytes',              ['megakaryocyte']),
 ('Neutrophils',                 ['neutrophil', 'granulocyte']),
 ('Mesothelial cells',           ['mesothelial cell']),
 ('Proliferating cells',         ['dividing cell', 'dividing t cell']),
])
# selectors that must NOT match (substring veto), to keep tumour/progenitor rows out
VETO = ['cancer', 'lymphoma', 'malignant', 'tumor', 'tumour', 'progenitor', 'stem cell',
        'myeloid-derived suppressor', 'megakaryocyte erythroid', 'megakaryocyte-erythroid']

rows = list(csv.DictReader(gzip.open(CM, 'rt')))
print('CellMarker 2.0 human, pinned: %d rows' % len(rows))
spleen_rows = [r for r in rows if r['tissue_class'] == 'Spleen']
print('tissue_class == Spleen: %d rows, %d cell_name values'
      % (len(spleen_rows), len({r['cell_name'] for r in spleen_rows})))


def matches(r, sels):
    n = (r['cell_name'] or '').lower()
    if any(v in n for v in VETO):
        return False
    return any(s in n for s in sels)


# ---- candidate pools -------------------------------------------------------------------
# Stage A: spleen rows only. Stage B, ONLY if stage A cannot support the label: all tissues, with
# a stricter PMID bar, because a marker validated in another organ is weaker evidence here. Every
# label that reaches stage B is flagged in the evidence table and in the report.
def pool_to_genes(pool, min_pmid):
    g2p = collections.defaultdict(set)
    for r in pool:
        s_ = (r['Symbol'] or '').strip()
        if s_ and s_.lower() != 'nan':
            g2p[s_].add(r['PMID'])
    return {g: p for g, p in g2p.items() if len(p) >= min_pmid}


def to_panel(g2p):
    d = {}
    for g, p in g2p.items():
        r = HS.resolve(g)
        if r:
            d.setdefault(r, set()).update(p)
    return d


cand, scope, tier = collections.OrderedDict(), {}, {}
for lab, sels in SELECTORS.items():
    sp  = [r for r in spleen_rows if matches(r, sels)]
    allp = [r for r in rows if matches(r, sels)]
    ladder = [
        (to_panel(pool_to_genes(sp, MIN_PMID_SPLEEN)) if len(sp) >= MIN_SPLEEN_ROWS else {},
         'spleen-only (%d rows)' % len(sp)),
        (to_panel(pool_to_genes(allp, MIN_PMID_FALLBACK)),
         'ALL-TISSUE FALLBACK >=%d PMID (%d spleen rows)' % (MIN_PMID_FALLBACK, len(sp))),
        (to_panel(pool_to_genes(allp, 1)),
         'ALL-TISSUE FALLBACK >=1 PMID -- WEAKEST EVIDENCE (%d spleen rows)' % len(sp)),
    ]
    chosen = None
    for ti, (d, sc) in enumerate(ladder):
        if len(d) >= MIN_MARKERS:
            chosen = (d, sc, ti)
            break
    if chosen is None:                       # nothing clears the gate; keep the richest pool so
        ti = max(range(3), key=lambda i: len(ladder[i][0]))   # the drop is reported with content
        chosen = (ladder[ti][0], ladder[ti][1], ti)
    cand[lab], scope[lab], tier[lab] = chosen

# ---- specificity, arbitration, size cap ------------------------------------------------
panel_hit = cand
nlab = collections.Counter()
for lab, d in panel_hit.items():
    for g in d:
        nlab[g] += 1

# ARBITRATION. CellMarker rows are per-publication marker PANELS, not exclusive markers, so a gene
# lands under 'CD4+ T cell' merely because some paper's panel also stained CD8A. Fine labels inside
# one compartment must be separated by genes that actually differ between them, so within each
# compartment a gene is kept ONLY by the label with the most PMIDs behind it (ties keep both).
# 'Proliferating cells' is a STATE, not a lineage: it keeps a gene only if it is the global argmax,
# which is what removes CD3D/IL7R from it.
ARBITRATE = collections.OrderedDict([
 ('B compartment',    ['Follicular B cells', 'Marginal zone B cells', 'Germinal centre B cells',
                       'Plasma cells']),
 ('T/NK compartment', ['CD4 T cells', 'CD8 T cells', 'NK cells']),
 ('Myeloid compartment', ['Red pulp macrophages', 'Monocytes', 'cDC1', 'cDC2', 'pDC',
                          'Neutrophils']),
 ('Endothelial compartment', ['Sinusoidal endothelium', 'Endothelial cells',
                              'Lymphatic endothelium']),
 ('Stromal compartment', ['Fibroblastic reticular cells', 'Fibroblasts', 'Pericytes',
                          'Smooth muscle / capsule', 'Mesothelial cells']),
])
lost_to_arbitration = collections.defaultdict(list)
for grp, labs in ARBITRATE.items():
    labs = [l for l in labs if l in panel_hit]
    genes = set().union(*[set(panel_hit[l]) for l in labs]) if labs else set()
    for g in genes:
        sup = {l: len(panel_hit[l][g]) for l in labs if g in panel_hit[l]}
        if len(sup) < 2:
            continue
        # TIER PRECEDENCE FIRST. A label whose markers come from an all-tissue fallback must never
        # take a gene from a label with real spleen evidence: the bigger pool wins for reasons that
        # are about the database, not the biology. (Without this, Neutrophils -- 5 spleen rows, 49
        # all-tissue candidates -- strips CD14/FCAR/S100A12 from Monocytes and deletes that label.)
        best_tier = min(tier[l] for l in sup)
        sup = {l: v for l, v in sup.items() if tier[l] == best_tier}
        if len(sup) < 2:
            for l in labs:
                if g in panel_hit[l] and tier[l] != best_tier:
                    del panel_hit[l][g]
                    lost_to_arbitration[l].append('%s->%s(stronger evidence tier)'
                                                  % (g, [x for x in sup][0]))
            continue
        best = max(sup.values())
        for l, v in sup.items():
            if v < best:
                del panel_hit[l][g]
                lost_to_arbitration[l].append('%s->%s' % (g, max(sup, key=sup.get)))
glob_best = {}
for g in nlab:
    sup = {l: len(panel_hit[l][g]) for l in panel_hit if g in panel_hit[l]}
    if sup:
        glob_best[g] = max(sup.values())
for g in list(panel_hit.get('Proliferating cells', {})):
    sup = {l: len(panel_hit[l][g]) for l in panel_hit if g in panel_hit[l]}
    if len(sup) > 1 and sup['Proliferating cells'] < max(sup.values()):
        del panel_hit['Proliferating cells'][g]
        lost_to_arbitration['Proliferating cells'].append('%s->%s' % (g, max(sup, key=sup.get)))

# ---- over-adjustment guard --------------------------------------------------------------
# The cell-type call is a Tier D NUISANCE COVARIATE. Building it out of the genes whose spatial
# behaviour is the OUTCOME conditions the response on itself -- the same over-adjustment the mouse
# arm guarded against when it stripped Tier A/Tier B genes out of the zonation covariate
# (genesets/README.md section 6). Applied here at Tier A + Tier C.
# NOT applied at Tier B: Tier B is 623 on-panel genes and removing them destroys six labels,
# including both stromal labels the A6 red-pulp/white-pulp covariate depends on. Tier B overlaps
# are therefore FLAGGED in the evidence table (column tierB_member) and reported, not removed.
import glob as _glob
_G = '/workspace/genesets/human/'
def _gl(n):
    return {l.strip() for l in open(_G + n) if l.strip()}
TIER_A = set().union(*[_gl(os.path.basename(f)) for f in _glob.glob(_G + 'A_*.txt')])
# the SEVEN canonical modules only -- not the C6 B7 variants, which are written later and
# would otherwise make this build order-dependent
TIER_B = set().union(*[_gl('B_%s.txt' % m) for m in
    ('tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
     'emt_ecm', 'oxidative_stress', 'secondary_senescence')])
TIER_C = _gl('C_ligands.txt') | _gl('C_receptors.txt')
GUARD = TIER_A | TIER_C
guard_removed = collections.defaultdict(list)
for lab in panel_hit:
    for g in list(panel_hit[lab]):
        if g in GUARD:
            del panel_hit[lab][g]
            guard_removed[lab].append(g)
print('\nOver-adjustment guard (Tier A + Tier C genes removed from the cell-type covariate):')
for l, v in sorted(guard_removed.items()):
    print('  %-30s %s' % (l, ' '.join(sorted(v))))

MARKERS, evidence, dropped_thin = collections.OrderedDict(), [], collections.OrderedDict()
for lab, d in panel_hit.items():
    surv = [g for g in d if nlab[g] <= MAX_LABELS]
    # cap to the best-evidenced MAX_PER_LABEL genes so sets stay comparable in size (the mouse
    # sets are 4-16 genes) and one over-populated CellMarker entry cannot dominate scoring
    surv = sorted(surv, key=lambda g: (-len(d[g]), g))[:MAX_PER_LABEL]
    keep = sorted(surv)
    for g in sorted(d):
        evidence.append(dict(cell_type=lab, gene=g, n_pmid=len(d[g]),
                             n_labels_claiming=nlab[g],
                             kept='yes' if g in keep else
                                  ('no_promiscuous' if nlab[g] > MAX_LABELS else 'no_below_cap'),
                             tierB_member='yes' if g in TIER_B else 'no',
                             scope=scope[lab], pmids=';'.join(sorted(d[g]))))
    if len(keep) >= MIN_MARKERS:
        MARKERS[lab] = keep
    else:
        dropped_thin[lab] = keep
    print('  %-30s %-42s on-panel=%-3d kept=%-3d %s'
          % (lab, scope[lab], len(d), len(keep),
             '' if len(keep) >= MIN_MARKERS else '<-- DROPPED, < %d markers' % MIN_MARKERS))

print('\nArbitration (gene kept by the best-evidenced label of its compartment):')
for l, v in sorted(lost_to_arbitration.items()):
    print('  %-30s lost %2d: %s' % (l, len(v), ' '.join(sorted(v))))
print('\nLABEL SET: %d assignable, %d dropped for thin on-panel support' % (len(MARKERS), len(dropped_thin)))
for k, v in dropped_thin.items():
    print('  DROPPED %-30s only %d on-panel markers: %s' % (k, len(v), ','.join(v)))
print('\nPromiscuous genes removed (claimed by > %d labels):' % MAX_LABELS)
for g, n in sorted(nlab.items(), key=lambda x: -x[1]):
    if n > MAX_LABELS:
        print('   %-10s %d labels' % (g, n))

with open(OUTE, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['cell_type', 'gene', 'n_pmid', 'n_labels_claiming',
                                       'kept', 'tierB_member', 'scope', 'pmids'])
    w.writeheader()
    for r in evidence:
        w.writerow(r)
print('\nwrote %s (%d rows)' % (OUTE, len(evidence)))

# ---- emit the static module ------------------------------------------------------------
HDR = '''"""Human spleen cell-type markers, Xenium Prime 5K Human panel (GSE326743, 5,093 genes).

GENERATED FILE -- do not hand-edit. Rebuild with:
    python3 /workspace/code/build_markers_human_spleen.py

TAXONOMY SOURCE: Phase 7 Job A follow-on task 1 (the compartments the coordinator specified:
  red pulp, white pulp T zone / follicle / marginal zone / germinal centre, FDC, sinusoidal vs
  other endothelium, FRC, plasma, T/NK, B, mono/DC, erythroid, megakaryocyte, capsule/trabecular
  smooth muscle). It replaces section 14's lung list, which does not apply to this arm.
MARKER SOURCE: CellMarker 2.0 human, pinned at genesets/cellmarker_pin/ (md5 recorded).
  *** No marker gene here was written from memory. *** Each gene carries its PMIDs in
  genesets/human/markers_spleen_evidence.csv. Filters: >= %d PMID, on-panel, claimed by <= %d
  labels, capped at the 15 best-evidenced genes, and >= %d surviving genes per label.

*** NOT VALIDATED ON DATA. *** The mouse equivalent (markers_mouse_liver.py / annotate_pipeline.py)
had its non-specific labels removed on MEASURED expression. H1 is behind the section 15 freeze, so
no H1 expression value has been read and no such check was possible. Re-gate this label set
against measured expression before trusting the fine labels.
"""
''' % (MIN_PMID_SPLEEN, MAX_LABELS, MIN_MARKERS)

MERGE = collections.OrderedDict([
 ('B cells',      [l for l in ('Follicular B cells', 'Marginal zone B cells',
                               'Germinal centre B cells') if l in MARKERS]),
 ('T/NK cells',   [l for l in ('CD4 T cells', 'CD8 T cells', 'NK cells') if l in MARKERS]),
 ('Mono/Mac/DC',  [l for l in ('Red pulp macrophages', 'Monocytes', 'cDC1', 'cDC2', 'pDC')
                   if l in MARKERS]),
 ('Endothelial',  [l for l in ('Sinusoidal endothelium', 'Endothelial cells',
                               'Lymphatic endothelium') if l in MARKERS]),
 ('Stromal',      [l for l in ('Fibroblastic reticular cells', 'Fibroblasts', 'Pericytes',
                               'Smooth muscle / capsule', 'Follicular dendritic cells')
                   if l in MARKERS]),
])
MERGE = collections.OrderedDict((k, v) for k, v in MERGE.items() if len(v) >= 2)

with open(OUTP, 'w') as fh:
    fh.write(HDR)
    fh.write('MARKERS = {\n')
    for k, v in MARKERS.items():
        fh.write("'%s':'%s'.split(),\n" % (k, ' '.join(v)))
    fh.write('}\n')
    fh.write('# Labels removed by the >= %d on-panel marker gate, with what survived:\n' % MIN_MARKERS)
    for k, v in dropped_thin.items():
        fh.write('#   %-30s %s\n' % (k, ','.join(v) or '(none)'))
    fh.write('\n# Fine labels that share a compartment. annotate_pipeline.py recomputes the\n'
             '# assignment over each group\'s union of markers and writes it as cell_type_merged.\n')
    fh.write('MERGE = {\n')
    for k, v in MERGE.items():
        fh.write(" '%s':%r,\n" % (k, v))
    fh.write('}\n')
    fh.write("\n# No entries: DROP_NONSPECIFIC is an EXPRESSION-based judgement and H1 expression\n"
             "# is behind the section 15 freeze. Populate it after the freeze, from measured data.\n"
             "DROP_NONSPECIFIC = {}\n")
print('wrote %s (%d labels)' % (OUTP, len(MARKERS)))
