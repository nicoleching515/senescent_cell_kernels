#!/usr/bin/env python3
"""Enforce the mechanisable half of PREREG_PHASE8.md section 10 -- "What may never be
reported, whatever the data says".

Until now all twelve prohibitions were honour-system prose: `reports/AUDIT_PREREG_VS_CODE.md`
PART 4 established that no literal, symbol or co-occurrence rule from section 10 appeared
anywhere in `code/`, `.githooks/` or `.claude/`, and that the one always-on hook
(`code/hook_geneset_gate.sh`) explicitly exits 0 on `*/results/*` and `*/reports/*` -- every
path where a violation would actually be written.

This checker is the mechanism.  It scans prose artefacts (`reports/*.md`, `README.md`,
`results/**/*.md`) and fires on:

  10.1   "negative-control probes" used for the POOLED `all_controls` response
  10.2   an affirmative caller-independence claim
  10.6   CXCL8/CXCR1 described as replicating a mouse result
  10.7   the circularity range 1.51-2.85x quoted
  10.8   the composition-matched SF 0.9837 quoted without 0.1461 beside it
  10.9   `mor` quoted without the `lib` result beside it
  10.10  a denoise=True number with no seed-stability companion
  10.11  the D2 `raw` rho_signed_dz_vs_depth values (-0.47, -0.16)
  10.12  any of the four struck sentences

10.3, 10.4 and 10.5 are semantic (an age claim, an MZ-specific claim, a species/tissue
attribution) and are NOT mechanised; `--list` prints them as REVIEW-ONLY so the list does
not imply a coverage it does not have.

Scoping.  A rule fires on a PARAGRAPH (blank-line-delimited block), not on a whole file, so
the co-occurrence rules (10.8, 10.9, 10.10) mean what the prereg means: the companion number
must appear BESIDE the forbidden one.

Two waivers, both deliberate and both narrow:
  * a paragraph that is *about* the prohibition -- it contains "prohibit", "must not",
    "may never", "never be quoted/reported", "struck", "forbidden", a section-10 citation, or
    an explicit waiver marker -- is meta-discussion, not a report of the number;
  * `PREREG_PHASE8.md` itself, which is where the rules and their forbidden literals live.
A per-line escape hatch exists for anything else: put `prohibition-waiver: 10.x <reason>`
in the paragraph.

Usage:
  check_prohibitions.py                 # default scan set, exit 1 on any violation
  check_prohibitions.py FILE [FILE...]  # scan specific files
  check_prohibitions.py --staged        # scan the staged versions of staged text files
  check_prohibitions.py --list          # print the rule table and stop
  check_prohibitions.py --self-test     # prove every mechanised rule fires, and that the
                                        # allowed companion form does NOT
"""
import os, re, subprocess, sys

ROOT = os.environ.get('SASP_ROOT', '/workspace')

# ---------------------------------------------------------------- rule helpers
def para_has(p, *pats):
    return all(re.search(x, p, re.I) for x in pats)


def _n(p):
    """Normalise a paragraph for literal matching: unicode minus/dashes -> ASCII,
    collapse whitespace, strip markdown emphasis.  A struck sentence must not slip
    through because it was re-typed with an en dash or wrapped across two lines."""
    p = p.replace('−', '-').replace('–', '-').replace('—', '-')
    p = p.replace('‘', "'").replace('’', "'")
    p = p.replace('“', '"').replace('”', '"')
    p = re.sub(r'[*`]', '', p)   # NB: not '_' -- rho_signed_dz_vs_depth must survive
    return re.sub(r'\s+', ' ', p)


STRUCK = [
    '0.93-1.22x of chance',
    'i.e. they are statistically independent',
    'Four of six pairs sit at 0.93-1.22x',
    'the one pair that looked concordant in sham is anti-concordant in SBR',
    "DeepScence's correlation with sequencing depth reverses sign between two sections of "
    'the same study',
]


def r_10_1(p, n):
    # the pooled all_controls response must be named "pooled negative-control features"
    if not re.search(r'all_controls', p, re.I):
        return None
    if re.search(r'negative[- ]control probes', p, re.I) and not re.search(r'pooled', p, re.I):
        return ('names the pooled `all_controls` response "negative-control probes"; '
                'section 10.1 requires "pooled negative-control features" '
                '(the 40 probes on their own are flat)')
    return None


def r_10_2(p, n):
    m = re.search(r'(?:are|is) statistically independent'
                  r'|callers? (?:are|were) independent'
                  r'|caller[- ]independen(?:ce|t)', n, re.I)
    if m and not re.search(r'\b(not|no|never|falsifi|reject|rule[sd]? out)\b',
                           n[max(0, m.start() - 120):m.end() + 120], re.I):
        return 'affirmative caller-independence claim (section 10.2 / P1)'
    return None


def r_10_6(p, n):
    if re.search(r'CXCL8|CXCR1', n) and re.search(r'replicat', n, re.I):
        return ('CXCL8/CXCR1 in a replication sentence; section 10.6 -- no mouse ortholog '
                'exists, so it can never replicate a mouse result')
    return None


def r_10_7(p, n):
    if re.search(r'1\.51\s*(?:-|to)\s*2\.85', n):
        return 'quotes the circularity range 1.51-2.85x (section 10.7)'
    if '1.51' in n and '2.85' in n:
        return 'quotes 1.51 and 2.85 together -- the withdrawn circularity range (section 10.7)'
    return None


def r_10_8(p, n):
    if '0.9837' in n and '0.1461' not in n:
        return ('quotes the composition-matched SF 0.9837 without the covariate-adjusted '
                '0.1461 beside it (section 10.8) -- alone it states the opposite of the data')
    return None


def r_10_9(p, n):
    if re.search(r'(?<![A-Za-z])mor(?![A-Za-z])', n) and re.search(r'normalis|normaliz|depth|caller', n, re.I):
        if not re.search(r'(?<![A-Za-z])lib(?![A-Za-z])', n):
            return ('quotes the `mor` normalisation result without the `lib` result beside it '
                    '(section 10.9)')
    return None


def r_10_10(p, n):
    if re.search(r'denoise\s*=\s*True', n):
        if not re.search(r'seed|stabilit|jaccard', n, re.I):
            return ('a denoise=True number with no seed-stability companion in the same '
                    'paragraph (section 10.10)')
    return None


def r_10_11(p, n):
    if re.search(r'rho_signed_dz_vs_depth', n) and re.search(r'(?<![\d.])-0\.47(?![\d])|(?<![\d.])-0\.16(?![\d])', n):
        if re.search(r'(?<![a-z])raw(?![a-z])', n, re.I):
            return ('quotes the D2 `raw` rho_signed_dz_vs_depth values (section 10.11) -- '
                    'the direction of numerical noise on a 0.0002-0.001 z-unit shift')
    return None


def r_10_12(p, n):
    for s in STRUCK:
        if s.lower() in n.lower():
            return 'reproduces a struck sentence (section 10.12): "%s"' % s[:60]
    return None


RULES = [
    ('10.1',  'name the pooled A7 response correctly',              r_10_1),
    ('10.2',  'no caller-independence claim',                       r_10_2),
    ('10.6',  'CXCL8/CXCR1 is not a mouse replication',             r_10_6),
    ('10.7',  '"1.51-2.85x" must not be quoted',                    r_10_7),
    ('10.8',  '0.9837 never without 0.1461',                        r_10_8),
    ('10.9',  '`mor` never without `lib`',                          r_10_9),
    ('10.10', 'no single-seed denoise=True number',                 r_10_10),
    ('10.11', 'D2 raw rho_signed_dz_vs_depth must not be quoted',   r_10_11),
    ('10.12', 'four struck sentences',                              r_10_12),
]

REVIEW_ONLY = [
    ('10.3', 'no age-stratified / young-vs-old claim on H1'),
    ('10.4', 'no marginal-zone-specific confirmatory claim'),
    ('10.5', 'no cross-arm species/tissue attribution'),
]

# A paragraph that is ABOUT the prohibition -- the rule itself, a withdrawal, a superseded
# passage kept as a record -- is not a report of the number.  The check also looks at the
# PRECEDING paragraph, because the "this is kept as a record" sentence is routinely the
# paragraph above the block quote it introduces.
META = re.compile(r'prohibit|must not|may never|(?:must|should|may|will) not be (?:quoted|reported)'
                  r'|never be (?:quoted|reported)|struck|forbidden|prohibition-waiver'
                  r'|section 10\.|§ ?10\.|10\.\d+\)'
                  r'|supersed|withdraw|withdrew|retract|no longer|does not survive'
                  r'|as the record|is a record|restate[ds]?|replaced by|kept below unaltered'
                  r'|correction c-|\[corrected|\[extended', re.I)

# PREREG_PHASE8.md is where the rules -- and therefore the forbidden literals -- live.
EXEMPT_FILES = {'reports/PREREG_PHASE8.md'}

SCAN_GLOBS = ['reports', 'results']
TEXT_EXT = {'.md', '.txt'}


def paragraphs(text):
    off, out = 0, []
    for block in re.split(r'\n\s*\n', text):
        line = text.count('\n', 0, off) + 1
        out.append((line, block))
        off += len(block) + 2
    return out


def check_text(rel, text):
    bad = []
    if rel.replace(os.sep, '/') in EXEMPT_FILES:
        return bad
    paras = paragraphs(text)
    prev_meta = False
    for i, (line, p) in enumerate(paras):
        if not p.strip():
            continue
        n = _n(p)
        meta = bool(META.search(n)) or prev_meta
        prev_meta = bool(META.search(n))
        for sec, _title, fn in RULES:
            msg = fn(p, n)
            if msg and not meta:
                bad.append((rel, line, sec, msg, _key(p)))
    return bad


def _key(p):
    """Identity of a violating paragraph, stable under re-wrapping and re-numbering, so the
    ratchet can tell a NEW violation from one that was already in HEAD."""
    return _n(p)[:400].lower()


def check_file(path, rel=None):
    rel = rel or os.path.relpath(path, ROOT)
    try:
        text = open(path, encoding='utf-8', errors='replace').read()
    except OSError as e:
        print('cannot read %s: %s' % (path, e), file=sys.stderr)
        return []
    return check_text(rel, text)


def default_files():
    out = []
    for top in SCAN_GLOBS:
        for dirpath, _d, names in os.walk(os.path.join(ROOT, top)):
            for f in sorted(names):
                if os.path.splitext(f)[1] in TEXT_EXT:
                    out.append(os.path.join(dirpath, f))
    rm = os.path.join(ROOT, 'README.md')
    if os.path.exists(rm):
        out.append(rm)
    return sorted(out)


def staged_files(new_only=True):
    """Check the staged content of staged prose files.

    With new_only (the pre-commit default) this is a RATCHET: a violation already present in
    HEAD's version of the same file does not block the commit, but any violation the commit
    ADDS does.  The corpus predates the checker and carries a backlog (see --backlog); the
    ratchet stops it growing without demanding a corpus-wide rewrite first, and without
    blocking unrelated work in files that already carry a backlog entry."""
    names = subprocess.run(['git', '-C', ROOT, 'diff', '--cached', '--name-only',
                            '--diff-filter=ACMR'],
                           capture_output=True, text=True).stdout.split()
    bad, carried = [], 0
    for rel in names:
        if os.path.splitext(rel)[1] not in TEXT_EXT:
            continue
        if not rel.startswith(('reports/', 'results/')) and rel != 'README.md':
            continue
        blob = subprocess.run(['git', '-C', ROOT, 'show', ':' + rel],
                              capture_output=True, text=True)
        if blob.returncode != 0:
            continue
        now = check_text(rel, blob.stdout)
        if not now:
            continue
        old = subprocess.run(['git', '-C', ROOT, 'show', 'HEAD:' + rel],
                             capture_output=True, text=True)
        known = {(v[2], v[4]) for v in check_text(rel, old.stdout)} if old.returncode == 0 else set()
        for v in now:
            if new_only and (v[2], v[4]) in known:
                carried += 1
            else:
                bad.append(v)
    if carried:
        print('section 10: %d pre-existing violation(s) carried unchanged in the staged files '
              '(not blocking; see check_prohibitions.py --backlog)' % carried)
    return bad


SELF_TEST = [
    ('10.1',  'The all_controls response is flat: negative-control probes at -0.0744.'),
    ('10.2',  'Across six pairs the callers are independent at chance overlap.'),
    ('10.6',  'CXCL8/CXCR1 replicates the mouse contact result in human spleen.'),
    ('10.7',  'The circularity figure is 1.51-2.85x across the two sections.'),
    ('10.8',  'Composition matching removes only 1.6% of the effect (SF 0.9837).'),
    ('10.9',  'Under mor the caller does not move, so normalisation cannot move it.'),
    ('10.10', 'With denoise=True the top-5% amplitude rises to 0.041 z-units.'),
    ('10.11', 'For the raw control rows rho_signed_dz_vs_depth is -0.47 and -0.16.'),
    ('10.12', 'Four of six pairs sit at 0.93-1.22x of chance.'),
]

# the same numbers reported the way section 10 REQUIRES them -- these must NOT fire
SELF_TEST_CLEAN = [
    'Composition matching leaves SF 0.9837 while the covariate-adjusted SF is 0.1461, '
    'i.e. type_adj 65.9% and typecomp_adj 85.4%.',
    'Under mor the caller barely moves; under lib the depth loading falls 74% and 100% of '
    'sender calls change.',
    'With denoise=True the amplitude is 0.041 at random_state=0; across three seeds the '
    'top-5% Jaccard is 0.000 for one of them.',
    'The pooled negative-control features (all_controls) sit at -0.0744.',
]


def self_test():
    ok = True
    for sec, text in SELF_TEST:
        hits = check_text('selftest/%s.md' % sec, text)
        got = {h[2] for h in hits}
        print('%-6s %s  %s' % (sec, 'FIRES ' if sec in got else 'MISSED', text[:64]))
        ok &= sec in got
    for text in SELF_TEST_CLEAN:
        hits = check_text('selftest/clean.md', text)
        print('clean  %s  %s' % ('OK    ' if not hits else 'FALSE-POSITIVE %s' % [h[2] for h in hits],
                                 text[:64]))
        ok &= not hits
    # the meta waiver must work, or the audit reports themselves become unshippable
    meta = check_text('selftest/meta.md',
                      'Section 10.7 prohibits quoting the range 1.51-2.85x.')
    print('waiver %s  meta paragraph about the rule' % ('OK    ' if not meta else 'BROKEN'))
    ok &= not meta
    print('SELF-TEST:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


def main(argv):
    if '--list' in argv:
        print('mechanised:')
        for sec, title, _f in RULES:
            print('  %-6s %s' % (sec, title))
        print('REVIEW-ONLY (semantic, not mechanised -- a human must check these):')
        for sec, title in REVIEW_ONLY:
            print('  %-6s %s' % (sec, title))
        return 0
    if '--self-test' in argv:
        return self_test()
    if '--backlog' in argv:
        bad = [v for f in default_files() for v in check_file(f)]
        print('section 10 backlog: %d pre-existing violation(s) in the committed corpus'
              % len(bad))
        for rel, line, sec, msg, _k in bad:
            print('  %s:%d  [%s] %s' % (rel, line, sec, msg))
        return 0
    if '--staged' in argv:
        bad = staged_files(new_only='--all-violations' not in argv)
        scope = 'staged text files'
    else:
        files = [a for a in argv if not a.startswith('-')]
        if files:
            bad = [v for f in files for v in check_file(os.path.abspath(f))]
            scope = '%d file(s)' % len(files)
        else:
            files = default_files()
            bad = [v for f in files for v in check_file(f)]
            scope = '%d prose artefacts' % len(files)
    if bad:
        print('SECTION 10 PROHIBITION VIOLATIONS (%d) in %s:' % (len(bad), scope),
              file=sys.stderr)
        for rel, line, sec, msg, _k in bad:
            print('  %s:%d  [%s] %s' % (rel, line, sec, msg), file=sys.stderr)
        print('\nIf the text is legitimately ABOUT the prohibition, say so in the paragraph '
              '(the words "prohibited"/"must not"/"struck", a section-10 citation, or\n'
              '"prohibition-waiver: 10.x <reason>").  Otherwise fix the text.', file=sys.stderr)
        return 1
    print('section 10: no prohibition violations in %s (%d mechanised rules; %d review-only)'
          % (scope, len(RULES), len(REVIEW_ONLY)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
