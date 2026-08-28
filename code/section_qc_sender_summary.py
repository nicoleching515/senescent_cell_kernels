#!/usr/bin/env python3
"""T-QC -- per-section QC + Cdkn1a sender summary for the 11 M1 (mouse GSE310392) sections.

WHAT THIS FILE IS
-----------------
This is the producer of /workspace/results/section_qc_sender_summary.csv, the table cited by
reports/BIO_PHASE3.md:436.  One row per M1 section, ordered (arm, timepoint_wk, section) --
the same ordering code/composition_all.py uses, which is why SBR precedes sham (capital 'S'
sorts first).

REGENERATE WITH
---------------
    /workspace/envs/sasp311/bin/python /workspace/code/section_qc_sender_summary.py \
        --out /workspace/results/section_qc_sender_summary.csv

Verify against the committed table without writing anything:

    /workspace/envs/sasp311/bin/python /workspace/code/section_qc_sender_summary.py --check

PROVENANCE -- READ THIS BEFORE TRUSTING THE TABLE
-------------------------------------------------
This script is a RECONSTRUCTION, written 2026-08-28 as part of the reproducibility repair.
The original producer of results/section_qc_sender_summary.csv was never committed.  Every
column below was reverse-engineered from the committed CSV plus the upstream code that made
its inputs (code/composition_all.py, code/phase2_downstream.py) and checked value-by-value
against the committed file.

Verified value-exact on all 11 sections (2026-08-28):
    section, arm, timepoint_wk
    cdkn1a_pos_pct_hepatocytes
    celltypes_in_1_20_band
    portal_triad_foci                (also cross-checked against logs/p2_rerun_m1.log and
                                      logs/downstream_p3_{1,2,3,4}.log, which print `nf`)
    corr_zonation_vs_dist_portal
    median_transcripts_per_cell
    median_cell_area_um2
    portal_triad_valid               (reproduces, but see the threshold caveat below --
                                      reproducing is not the same as being recovered)

Differ, and are EXPECTED to differ -- do not "fix" them:
    n_analysable, n_cdkn1a_pos, cdkn1a_pos_pct_all
    for sections 7239, 7435, 7248, 7001.  Those four were re-annotated after the committed
    CSV was written (commit "Restore ..."/reannotate.py lineage), so data/processed/
    celltypes_<tag>.csv today carries a different Low_quality/Unknown split than it did then.
    7435 is the loudest case: 139,768 analysable cells committed vs 151,291 today.  The
    recomputation is CORRECT for today's annotation; the committed number is correct for the
    annotation of the day.  The other seven sections still match exactly, which is what tells
    us the recipe itself is right rather than the annotation drift hiding a bug.

    cdkn1a_pos_pct_hepatocytes is unaffected even on those four sections, because the
    re-annotation moved cells between Low_quality/Unknown and the non-hepatocyte types, not
    into or out of Hepatocytes.

INVENTED VALUE -- portal_triad_valid
------------------------------------
No rule for portal_triad_valid is recorded anywhere in the repo, the logs, or the reports.
All that can be inferred from the committed CSV is that it is YES at corr = +0.198 (7361) and
NO at corr = +0.121 (7450), so the cut sits somewhere in (0.121, 0.198].  The constant
PORTAL_TRIAD_VALID_MIN_CORR below is a RECONSTRUCTION -- a midpoint I chose, NOT a recovered
value.  It reproduces all 11 committed YES/NO calls, but so would any other number in that
interval, and a 12th section could easily land between the real threshold and this one.  The
script prints this fact on every run so it cannot be quietly forgotten.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.cluster import DBSCAN

PROC = '/workspace/data/processed/'
RAW  = '/workspace/data/raw/'
COMMITTED = '/workspace/results/section_qc_sender_summary.csv'

# GEO sample table, copied verbatim from code/composition_all.py so the two tables cannot
# drift apart in which sections they consider "the M1 set".
META = {'7361_liver_sbr_Male_2-U1'  : ('7361', 'SBR',  2),
        '7352_liver_sham_Male_2-U1' : ('7352', 'sham', 2),
        '7448_liver_sbr_Male_10-U1' : ('7448', 'SBR',  10),
        '7450_liver_sbr_Male_10-U1' : ('7450', 'SBR',  10),
        '7435_liver_sham_Male_10-U1': ('7435', 'sham', 10),
        '7259_liver_sbr_Male_26-U1' : ('7259', 'SBR',  26),
        '7260_liver_sbr_Male_26-U1' : ('7260', 'SBR',  26),
        '7248_liver_sham_Male_26-U1': ('7248', 'sham', 26),
        '7250_liver_sham_Male_26-U1': ('7250', 'sham', 26),
        '7239_liver_sbr_Male_52-U1' : ('7239', 'SBR',  52),
        '7001_liver_sham_Male_52-U1': ('7001', 'sham', 52)}

# Same exclusion set as code/composition_all.py:39.  These two labels are annotation-quality
# artefacts, not cell types, so they cannot be receivers and must not sit in the denominator.
EXCL = ['Low_quality', 'Unknown']

# --- RECONSTRUCTION, NOT A RECOVERED VALUE.  See the module docstring. -------------------
# portal_triad_valid asks: did the bile-duct DBSCAN foci actually land on portal triads?  The
# check is that zonation increases away from them (pericentral = far from the portal triad),
# i.e. corr_zonation_vs_dist_portal must be convincingly positive.  The committed CSV pins the
# real cut only to the half-open interval (0.121, 0.198]; 0.15 is the midpoint I picked.
PORTAL_TRIAD_VALID_MIN_CORR = 0.15
PORTAL_TRIAD_VALID_IS_RECONSTRUCTION = True

# Committed column order.  Reproduced exactly so a diff against the committed file is a
# content diff and not a column-shuffle diff.
COLUMNS = ['section', 'arm', 'timepoint_wk', 'n_analysable', 'cdkn1a_pos_pct_all',
           'cdkn1a_pos_pct_hepatocytes', 'n_cdkn1a_pos', 'celltypes_in_1_20_band',
           'portal_triad_foci', 'corr_zonation_vs_dist_portal', 'portal_triad_valid',
           'median_transcripts_per_cell', 'median_cell_area_um2']


def section_row(tag):
    """Compute one section's row.  Reads only the parquet columns it needs -- cells.parquet
    carries 14 columns and we want 5, and this loop runs over 11 sections on a shared box."""
    sec, arm, wk = META[tag]

    # --- analysable set and Cdkn1a positivity -------------------------------------------
    ct = pd.read_csv(PROC + 'celltypes_%s.csv' % tag,
                     usecols=['cell_id', 'cell_type', 'cell_type_merged'])
    an = ct[~ct.cell_type_merged.isin(EXCL)]
    sd = pd.read_csv(PROC + 'senders_%s.csv' % tag,
                     usecols=['cell_id', 'cdkn1a_pos']).set_index('cell_id')
    # reindex, not merge: senders_*.csv is written from the same cell set, but if a
    # re-annotation ever drops cells we want NaN (and a visible count mismatch) rather than a
    # silent inner join that shrinks the denominator.
    pos_all = sd.cdkn1a_pos.reindex(an.cell_id).to_numpy()
    n_an = len(an)
    n_pos = int(np.nansum(pos_all))

    hep = (an.cell_type_merged == 'Hepatocytes').to_numpy()
    pos_hep = sd.cdkn1a_pos.reindex(an.cell_id[hep]).to_numpy()

    # --- prevalence band -----------------------------------------------------------------
    # "k/n": how many of the section's cell types put Cdkn1a prevalence inside the 1-20%
    # band that the sender calls assume.  n is the number of rows in test3_prevalence_*, i.e.
    # the number of cell types that survived that test's own size floor -- which is why the
    # denominator moves between sections (9..12).
    t3 = pd.read_csv(PROC + 'test3_prevalence_%s.csv' % tag)
    band = '%d/%d' % (int((t3.cdkn1a_in_1_20_band == 'YES').sum()), len(t3))

    # --- raw QC medians ------------------------------------------------------------------
    # Deliberately over ALL rows of cells.parquet, not the analysable subset: these are
    # instrument/segmentation QC numbers and must not be conditioned on the annotation, or
    # they would move every time the annotation is redone.
    cells = pq.read_table(RAW + tag + '/cells.parquet',
                          columns=['cell_id', 'x_centroid', 'y_centroid',
                                   'transcript_counts', 'cell_area']).to_pandas()
    med_tx = float(np.median(cells.transcript_counts.to_numpy()))
    med_area = float(np.median(cells.cell_area.to_numpy()))

    # --- zonation sanity check -----------------------------------------------------------
    # Hepatocytes only: zonation_score is standardised on hepatocytes upstream
    # (code/phase2_downstream.py), so the correlation is only meaningful within them.
    a = pd.read_csv(PROC + 'anatomy_%s.csv' % tag,
                    usecols=['cell_type', 'zonation_score', 'dist_to_portal_triad_um'])
    h = a[a.cell_type == 'Hepatocytes']
    corr = float(np.corrcoef(h.zonation_score.to_numpy(),
                             h.dist_to_portal_triad_um.to_numpy())[0, 1])

    # --- portal triad foci ---------------------------------------------------------------
    # Identical to the `nf` computed at code/phase2_downstream.py:57-61.  anatomy_*.csv keeps
    # the distance but not the cluster count, so it has to be recomputed here from the same
    # inputs: DBSCAN over the xy of the Biliary/ductular cells, using the FINE cell_type
    # column (not cell_type_merged), on the cells present in celltypes_*.csv.
    xy = cells.set_index('cell_id').reindex(ct.cell_id)[['x_centroid', 'y_centroid']].to_numpy()
    chol = xy[(ct.cell_type == 'Biliary/ductular').to_numpy()]
    nf = 0
    if len(chol) >= 30:                       # the >=30 floor is phase2_downstream's, kept
        lab = DBSCAN(eps=30, min_samples=10).fit(chol).labels_
        nf = int(len(np.unique(lab[lab >= 0])))   # -1 is DBSCAN noise, not a focus

    return {'section': sec, 'arm': arm, 'timepoint_wk': wk,
            'n_analysable': n_an,
            'cdkn1a_pos_pct_all': round(100.0 * n_pos / n_an, 2),
            'cdkn1a_pos_pct_hepatocytes': round(100.0 * float(np.nansum(pos_hep)) / int(hep.sum()), 2),
            'n_cdkn1a_pos': n_pos,
            'celltypes_in_1_20_band': band,
            'portal_triad_foci': nf,
            'corr_zonation_vs_dist_portal': round(corr, 3),
            'portal_triad_valid': 'YES' if corr >= PORTAL_TRIAD_VALID_MIN_CORR else 'NO',
            'median_transcripts_per_cell': med_tx,
            'median_cell_area_um2': round(med_area, 1)}


def build():
    rows = [section_row(tag) for tag in sorted(META)]
    D = pd.DataFrame(rows)[COLUMNS]
    # Same key as code/composition_all.py: arm then timepoint then section.  'SBR' < 'sham'
    # because capitals sort first -- that is why the committed file leads with SBR.
    return D.sort_values(['arm', 'timepoint_wk', 'section']).reset_index(drop=True)


def banner(stream):
    stream.write('# RECONSTRUCTION (2026-08-28). portal_triad_valid uses an INVENTED threshold:\n'
                 '#   PORTAL_TRIAD_VALID_MIN_CORR = %.3f, chosen from the interval (0.121, 0.198]\n'
                 '#   that the committed CSV pins it to. No rule for it is recorded anywhere.\n'
                 '# n_analysable / n_cdkn1a_pos / cdkn1a_pos_pct_all differ from the committed\n'
                 '#   file on 7239, 7435, 7248, 7001: those sections were re-annotated after it\n'
                 '#   was written. That difference is expected and must not be patched away.\n'
                 % PORTAL_TRIAD_VALID_MIN_CORR)


def check(D):
    """Cell-by-cell comparison against the committed CSV.  Prints MATCH/DIFF per column per
    section with the magnitude of every difference.  Exit status is 0 either way -- the known
    re-annotation differences are not a failure, and a non-zero exit would train people to
    ignore it."""
    C = pd.read_csv(COMMITTED, dtype={'section': str})
    D = D.copy(); D['section'] = D.section.astype(str)
    C = C.set_index('section'); Dx = D.set_index('section')

    if list(C.index) != list(Dx.index):
        print('SECTION ORDER/SET DIFFERS: committed %s vs rebuilt %s' % (list(C.index), list(Dx.index)))
    order = [s for s in C.index if s in Dx.index]

    print('\nper-column / per-section comparison vs %s' % COMMITTED)
    print('-' * 96)
    hdr = '%-30s ' % 'column' + ' '.join('%7s' % s for s in order)
    print(hdr); print('-' * 96)
    summary = []
    for col in COLUMNS[1:]:
        cells_out, ndiff = [], 0
        for s in order:
            a, b = C.loc[s, col], Dx.loc[s, col]
            if isinstance(a, str) or isinstance(b, str):
                same = str(a) == str(b)
                mag = '' if same else '%s->%s' % (a, b)
            else:
                same = bool(np.isclose(float(a), float(b), rtol=0, atol=1e-9))
                mag = '' if same else '%+g' % (float(b) - float(a))
            if same:
                cells_out.append('%7s' % '.')
            else:
                ndiff += 1
                cells_out.append('%7s' % 'DIFF')
            summary.append((col, s, same, mag))
        print('%-30s ' % col + ' '.join(cells_out) + ('   ALL MATCH' if ndiff == 0
                                                      else '   %d DIFF' % ndiff))
    print('-' * 96)
    print("('.' = exact match)\n")
    diffs = [t for t in summary if not t[2]]
    if diffs:
        print('differences, with magnitude (committed -> rebuilt):')
        for col, s, _, mag in diffs:
            print('  %-28s %-6s  %s' % (col, s, mag))
    else:
        print('no differences.')
    print('\n%d of %d cells match (%d columns x %d sections).'
          % (sum(1 for t in summary if t[2]), len(summary), len(COLUMNS) - 1, len(order)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', default=None,
                    help='write the CSV here; default is stdout. Nothing is written unless '
                         'this is given.')
    ap.add_argument('--check', action='store_true',
                    help='regenerate and compare cell-by-cell against %s; writes nothing.'
                         % COMMITTED)
    args = ap.parse_args()

    D = build()
    banner(sys.stderr)
    if args.check:
        check(D)
        return
    if args.out:
        D.to_csv(args.out, index=False)
        sys.stderr.write('wrote %s (%d rows)\n' % (args.out, len(D)))
    else:
        banner(sys.stdout)
        D.to_csv(sys.stdout, index=False)


if __name__ == '__main__':
    main()
