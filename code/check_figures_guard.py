#!/usr/bin/env python3
"""Figure integrity guard — Phase 8, D8/D10.

Two actors collided in figures/ on 2026-08-27: a committed baseline was
restored, then an agent regenerated figure2b/2c over it. Nothing warned. The
PI's figure policy is that COMMITTED figures are held at their committed state
and regenerated exactly once, from the frozen configuration, at task 8.7.

This makes that policy enforceable rather than advisory.

  python3 code/check_figures_guard.py           # verify, exit 1 on drift
  python3 code/check_figures_guard.py --snapshot  # record current as expected

PNG/CSV are compared by content hash. PDFs are compared with matplotlib's
embedded CreationDate/ModDate stripped, because those differ on every
regeneration even when the figure is identical.
"""
import hashlib, json, os, re, subprocess, sys

ROOT = '/workspace'
FIGS = os.path.join(ROOT, 'figures')
MANIFEST = os.path.join(ROOT, 'figures', '.committed_manifest.json')
DATESTAMP = re.compile(rb'/(CreationDate|ModDate)\s*\([^)]*\)')


def digest(path):
    data = open(path, 'rb').read()
    if path.endswith('.pdf'):
        data = DATESTAMP.sub(b'', data)
    return hashlib.sha256(data).hexdigest()


def tracked_figures():
    """Every figure artefact in figures/, tracked or not.

    GAP CLOSED 2026-08-27 (task 8.7 findings). This used `git ls-files`, so the
    18 untracked Phase 8 artefacts (figure2e, figure_gs1-4, figure_phase8_*) were
    outside its scope. During the M1 re-run another agent rewrote
    figure_gs2_crossarm_symmetry and figure_gs3_corescence_circularity and
    NOTHING WARNED, because untracked files were not being watched. The guard's
    whole purpose is detecting exactly that.

    revised_candidates/ is excluded: it holds deliberately superseded copies.
    """
    out = []
    for dirpath, dirnames, filenames in os.walk(FIGS):
        dirnames[:] = [d for d in dirnames if d != 'revised_candidates']
        for fn in filenames:
            if fn.rsplit('.', 1)[-1] in ('png', 'pdf', 'csv'):
                out.append(os.path.relpath(os.path.join(dirpath, fn), ROOT))
    return sorted(out)


def snapshot():
    m = {f: digest(os.path.join(ROOT, f)) for f in tracked_figures()
         if os.path.exists(os.path.join(ROOT, f))}
    json.dump(m, open(MANIFEST, 'w'), indent=1, sort_keys=True)
    print('snapshot: %d committed figures recorded -> %s' % (len(m), MANIFEST))


def verify():
    if not os.path.exists(MANIFEST):
        sys.exit('no manifest; run --snapshot first')
    expected = json.load(open(MANIFEST))
    drift, missing = [], []
    for f, want in expected.items():
        p = os.path.join(ROOT, f)
        if not os.path.exists(p):
            missing.append(f); continue
        if digest(p) != want:
            drift.append(f)
    if not drift and not missing:
        print('OK: all %d committed figures match (PDF date stamps ignored)' % len(expected))
        return 0
    for f in missing:
        print('MISSING: %s' % f)
    for f in drift:
        print('CHANGED: %s' % f)
    print('\nCommitted figures are held until task 8.7 regenerates them from the '
          'frozen configuration.\nRestore with:  git checkout -- figures/')
    return 1


if __name__ == '__main__':
    if '--snapshot' in sys.argv:
        snapshot()
    else:
        sys.exit(verify())
