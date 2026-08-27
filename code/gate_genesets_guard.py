#!/usr/bin/env python3
"""Gene-set / panel change guard — Phase 8, hazard B6.

WHY THIS EXISTS. `B_oxidative_stress` clears the section-11 floor of >=30 by
exactly ONE gene on the mouse arm (31), and that gene is `Junb`, which comes
from GSE310392's 100-gene custom add-on rather than the stock Prime Mouse 5K
panel. Under the stock-panel-only definition it sits AT the floor, margin 0.
Any later edit that trims mouse B6 by one gene -- a QC filter, an MSigDB bump,
dropping the curated NRF2 component of the documented B6 rescue, or a change to
which panel file is authoritative -- fails the gate. The gate must therefore run
after ANY gene-set or panel change, not only at freeze time.

WHAT IT DOES. Two things, in this order:

  1. Reports which watched file changed since the recorded manifest, by SHA-256
     (`genesets/.geneset_manifest.json`). This is the *trigger* evidence: it
     names what moved.
  2. Runs the section-11 gate on BOTH arms and exits non-zero if either fails.
     This is the *test*.

     * HUMAN arm: delegates to `code/gate_disjointness_human.py` unchanged, as a
       subprocess. That script is the gate of record for H1 and already exits
       non-zero on failure; it is wired in here, not reimplemented.
     * MOUSE arm: the existing mouse gate lives inside
       `code/build_genesets_mouse_c6.py`, which gates the sets it is BUILDING in
       memory and writes files as a side effect. That is the wrong instrument
       for a post-change check: it cannot be run safely mid-freeze, and it does
       not test the promoted `genesets/*.txt` that the pipeline actually reads.
       So the mouse half here re-applies the same section-11 assertions,
       read-only, to the files on disk. The panel construction and its three
       count assertions (5,106 features / 9 genotyping probes / 5,097 genes) are
       kept identical to `build_genesets_mouse_c6.py` so that a drift in either
       one is a loud failure rather than a silent divergence.

Nothing outside `genesets/.geneset_manifest.json` is written, and that only on
--snapshot. Reads no expression data on either arm; the human half reads panel
membership only (`genesets/h1_candidate/`), which is the sanctioned section-12.1
screen, and never touches `data/raw_h1/`.

  python3 code/gate_genesets_guard.py             # verify + gate, exit 1 on gate failure
  python3 code/gate_genesets_guard.py --snapshot  # record the current state as expected
  python3 code/gate_genesets_guard.py --quiet     # only print on drift or failure
"""
import csv, glob, gzip, hashlib, json, os, subprocess, sys

W = '/workspace'
MANIFEST = os.path.join(W, 'genesets', '.geneset_manifest.json')
HUMAN_GATE = os.path.join(W, 'code', 'gate_disjointness_human.py')

MODULES = ['tnfa_nfkb_proximal', 'il6_jak_stat3', 'interferon_response', 'downstream_arrest',
           'emt_ecm', 'oxidative_stress', 'secondary_senescence']
BLAB = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
MIN_A, MIN_B = 15, 30

# Everything whose change can move the gate: both arms' gene sets, both arms'
# panel definitions, and the frozen manifest that records what is canonical.
WATCH_GLOBS = [
    'genesets/*.txt',
    'genesets/human/*.txt', 'genesets/human/variants/*.txt',
    'genesets/mouse_c6/*.txt', 'genesets/mouse_c6/variants/*.txt',
    'genesets/human/FROZEN_MANIFEST.csv',
    'genesets/h1_candidate/GSE326743_gene_panel_5093.csv',
    'XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv',
    'GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz',
]


def watched():
    out = set()
    for g in WATCH_GLOBS:
        out.update(os.path.relpath(p, W) for p in glob.glob(os.path.join(W, g)))
    return sorted(out)


def digest(rel):
    return hashlib.sha256(open(os.path.join(W, rel), 'rb').read()).hexdigest()


def snapshot():
    m = {f: digest(f) for f in watched()}
    json.dump(m, open(MANIFEST, 'w'), indent=1, sort_keys=True)
    print('snapshot: %d watched gene-set/panel files recorded -> %s' % (len(m), MANIFEST))


def drift():
    """(changed, added, removed) relative to the recorded manifest."""
    if not os.path.exists(MANIFEST):
        return None
    want = json.load(open(MANIFEST))
    now = {f: digest(f) for f in watched()}
    changed = sorted(f for f in want if f in now and now[f] != want[f])
    added = sorted(set(now) - set(want))
    removed = sorted(set(want) - set(now))
    return changed, added, removed


# ------------------------------------------------------------------ mouse panel
# Construction and assertions copied verbatim from code/build_genesets_mouse_c6.py
# so the two cannot silently disagree about what the mouse panel is.
def mouse_panel():
    panel_meta = {r['gene_name'] for r in
                  csv.DictReader(open(W + '/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv'))}
    p100 = [r['Gene'] for r in csv.DictReader(
            gzip.open(W + '/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz', 'rt'))]
    genotyping = {g for g in p100
                  if ('_WT' in g or '_ALT' in g or '_del_' in g or '_splice_' in g)}
    features = panel_meta | set(p100)
    panel = features - genotyping
    assert len(features) == 5106, 'mouse h5 feature count moved: %d' % len(features)
    assert len(genotyping) == 9, 'genotyping probe count moved: %d' % len(genotyping)
    assert len(panel) == 5097, 'authoritative mouse panel moved: %d' % len(panel)
    return panel, panel_meta


def gl(path):
    return {l.strip() for l in open(path) if l.strip()}


def mouse_gate(verbose=True):
    """Section 11, read-only, on the PROMOTED genesets/*.txt the pipeline reads."""
    panel, csv_only = mouse_panel()
    gs = os.path.join(W, 'genesets')
    B = {m: gl(os.path.join(gs, 'B_%s.txt' % m)) & panel for m in MODULES}
    A = gl(os.path.join(gs, 'A_SENDER_FINAL_strict.txt')) & panel
    Amod = {m: gl(os.path.join(gs, 'A_sender_for_%s.txt' % m)) & panel for m in MODULES}

    rows, ok = [], True
    if verbose:
        print('\n### SECTION 11 GATE -- MOUSE ARM, promoted genesets/*.txt, %d-gene panel ###'
              % len(panel))
    c = len(A) >= MIN_A
    ok &= c
    rows.append(('len(A_SENDER_FINAL_strict) >= %d' % MIN_A, len(A), c))
    for m, lab in zip(MODULES, BLAB):
        c = len(B[m]) >= MIN_B
        ok &= c
        margin = len(B[m]) - MIN_B
        note = '  <-- NO MARGIN (add-on gene: %s)' % ' '.join(sorted(B[m] - csv_only)) \
               if margin <= 1 else ''
        rows.append(('len(%s_%s) >= %d%s' % (lab, m, MIN_B, note), len(B[m]), c))
    for m, lab in zip(MODULES, BLAB):
        c = len(A & B[m]) == 0
        ok &= c
        rows.append(('A_strict n %s_%s == 0' % (lab, m), len(A & B[m]), c))
    for m, lab in zip(MODULES, BLAB):
        c = len(Amod[m]) >= MIN_A and not (Amod[m] & B[m])
        ok &= c
        rows.append(('sensitivity A_sender_for_%s: n>=%d and disjoint from %s'
                     % (m, MIN_A, lab), len(Amod[m]), c))
    if verbose:
        for label, n, c in rows:
            print('  [%s] %-72s : %d' % ('PASS' if c else 'FAIL', label, n))
        print('  ==> MOUSE GATE: %s' % ('PASS' if ok else 'FAIL'))
    return ok


def human_gate(verbose=True):
    r = subprocess.run([sys.executable, HUMAN_GATE], capture_output=True, text=True)
    ok = r.returncode == 0
    if verbose:
        print('\n### SECTION 11 GATE -- HUMAN ARM (delegated to gate_disjointness_human.py) ###')
        for line in r.stdout.strip().splitlines():
            if 'PASS' in line or 'FAIL' in line or 'VERDICT' in line or 'FROZEN' in line:
                print('  ' + line.strip())
        print('  ==> HUMAN GATE: %s (exit %d)' % ('PASS' if ok else 'FAIL', r.returncode))
    if not ok and not verbose:
        sys.stderr.write(r.stdout[-3000:] + r.stderr[-2000:])
    return ok


def main():
    args = set(sys.argv[1:])
    if '--snapshot' in args:
        snapshot()
        return 0
    quiet = '--quiet' in args

    d = drift()
    moved = False
    if d is None:
        print('WARNING: no manifest at %s; run --snapshot to record a baseline.' % MANIFEST)
    else:
        changed, added, removed = d
        moved = bool(changed or added or removed)
        if moved:
            print('GENE-SET / PANEL CHANGE DETECTED since the recorded manifest:')
            for f in changed:
                print('  CHANGED %s' % f)
            for f in added:
                print('  ADDED   %s' % f)
            for f in removed:
                print('  REMOVED %s' % f)
            print('Running the section-11 gate on both arms.')
        elif not quiet:
            print('no drift: %d watched gene-set/panel files match the manifest' % len(watched()))

    verbose = not quiet or moved
    m_ok = mouse_gate(verbose)
    h_ok = human_gate(verbose)
    if m_ok and h_ok:
        if verbose:
            print('\nGATE PASS (mouse + human). B6 margin above is the one to watch.')
        return 0
    sys.stderr.write('\nGATE FAIL: mouse=%s human=%s -- '
                     'a gene set or panel change has broken section 11. '
                     'Do not proceed to a fit or a freeze.\n'
                     % ('PASS' if m_ok else 'FAIL', 'PASS' if h_ok else 'FAIL'))
    return 1


if __name__ == '__main__':
    sys.exit(main())
