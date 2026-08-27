#!/usr/bin/env python3
"""Phase 8 / 8.7 -- the closing audit.  Every claim the report makes about the
state of the tree at the end of the re-run is produced here, from the files on
disk, in one pass."""
import glob, hashlib, json, os, subprocess, sys
import numpy as np, pandas as pd
R3 = '/workspace/results/phase3'
R5 = '/workspace/results/phase5'


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


print('=' * 78)
print('1. PERMUTATION COUNTS  (Sec 24.3 audit trail)')
print('=' * 78)
for f in sorted(glob.glob(f'{R3}/perm_nulls*.csv')):
    if 'sbrscope' in f:
        continue
    head = pd.read_csv(f, nrows=0).columns
    d = pd.read_csv(f, usecols=[c for c in ('n_perm', 'call', 'sender_set')
                                if c in head])
    ts = pd.Timestamp(os.path.getmtime(f), unit='s').strftime('%Y-%m-%d %H:%M')
    print('  %-26s n_perm=%-8s calls=%-58s %s'
          % (os.path.basename(f), sorted(d.n_perm.unique()),
             ','.join(sorted(d.call.unique()))[:58], ts))

print()
print('=' * 78)
print('2. PINNED FILES')
print('=' * 78)
PIN = {'results/phase3/perm_nulls.csv': '3b77aa1bba0712c205c5d9356654fb71',
       'results/phase3/sf_summary.csv': '69e3a1d3f60060deddcceba9896a7d31',
       'results/phase3/summary_phase3.txt': 'ecf86b9ca5460f31290e2f4c9e822ea2',
       'data/processed/deepscence_sham.csv': '8c4c52f5c1c7649d8c17d07010cc780c',
       'data/processed/deepscence_sbr.csv': 'b557e3dfb8eff517d040757c73f0a660'}
for p, want in PIN.items():
    got = md5('/workspace/' + p)
    print('  %-40s %s  %s' % (p, got,
                              'UNCHANGED' if got == want else 'SUPERSEDED (was %s)' % want))

print()
print('=' * 78)
print('3. HEADLINE VECTOR, pre-C6 vs post-C6  (same extractor, both trees)')
print('=' * 78)


def headline(r3, r5):
    out = subprocess.run([sys.executable, '/workspace/code/m1_headlines.py', r3, r5],
                         capture_output=True, text=True)
    return json.loads(out.stdout)


a = headline('/workspace/results/phase3_pre_c6', '/workspace/results/phase5_pre_c6')
b = headline(R3, R5)
KEYS = ['n_fits', 'n_reportable', 'naive_amp_med', 'ctrl_amp_med',
        'ctrl_amp_se_med', 'power80_bound', 'ctrl_pos_and_sig',
        'SF_N2', 'SF_N5', 'SF_N6', 'SF_zon', 'SF_N5+N6', 'SF_N2+N5+N6',
        'SF_N8', 'N1', 'N1_full', 'N3', 'N4', 'lam_railed_frac',
        'poisson_slope', 'poisson_r2', 'ripley_ratio_med',
        'c1_N3_tile', 'c1_N3_occ', 'c1_N3_occ15', 'c1_N3_swap', 'c1_N3_snap',
        'c1_N3_orig', 'c1_N4_tile', 'c1_N4_occ', 'c1_N4_occ15', 'c1_N4_swap',
        'c1_N4_orig', 'c1_N3_swap_fullsf']
print('  %-22s %14s %14s %10s' % ('quantity', 'pre-C6', 'post-C6', 'delta'))
print('  ' + '-' * 64)
for k in KEYS:
    x, y = a.get(k), b.get(k)
    if x is None or y is None:
        print('  %-22s %14s %14s' % (k, x, y)); continue
    print('  %-22s %14.4f %14.4f %+10.4f' % (k, x, y, y - x))
for k in ('ctrl_amp_iqr', 'SF_N2+N5+N6_iqr'):
    print('  %-22s %s -> %s' % (k, np.round(a.get(k, []), 4),
                                np.round(b.get(k, []), 4)))

print()
print('=' * 78)
print('4. THE SECOND PRE-REGISTERED TIER A VARIANT')
print('=' * 78)
m = pd.read_csv(f'{R3}/main_fits.csv')
IB = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
      "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
      "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
m = m[m.section.isin(IB) & (m.stratum == 'all')]
rep = m[(m.beta_naive > 0) & (m.beta_base_lo > 0)]
print('  %-14s %5s %9s %9s %9s' % ('call', 'rep', 'b/sd', 'SF_N5', 'SF_N2N5N6'))
for c in ['tierA_p90', 'tierApm_p90', 'tierA_p95', 'tierApm_p95',
          'tierA_p99', 'tierApm_p99']:
    r = rep[rep.call == c]
    if len(r):
        print('  %-14s %5d %9.3f %9.3f %9.3f'
              % (c, len(r), (r.beta_naive / r.sd_y).median(),
                 r.sf_n5.median(), r.sf_n2n5n6.median()))

print()
print('=' * 78)
print('5. FIGURE STATE')
print('=' * 78)
g = subprocess.run([sys.executable, '/workspace/code/check_figures_guard.py'],
                   capture_output=True, text=True)
print(' ', g.stdout.strip().replace('\n', '\n  '))
print('  guard exit code:', g.returncode)
for p in sorted(glob.glob('/workspace/figures/*.png')):
    print('  %s  %s  %s' % (md5(p), os.path.basename(p),
                            pd.Timestamp(os.path.getmtime(p), unit='s').strftime('%H:%M')))
print()
print('  data CSVs beside each figure:')
for p in sorted(glob.glob('/workspace/figures/*_data.csv')
                + glob.glob(f'{R3}/figure2*_data.csv')):
    print('    %-58s %s' % (p.replace('/workspace/', ''),
                            pd.Timestamp(os.path.getmtime(p), unit='s').strftime('%H:%M')))
