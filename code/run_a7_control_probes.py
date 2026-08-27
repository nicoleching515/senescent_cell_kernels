#!/usr/bin/env python3
"""Phase 8, task 8.5 -- audit test A7 on the MOUSE arm (M1).

Addendum §13, test A7: fit the kernel against NEGATIVE-CONTROL-PROBE counts
versus distance-to-nearest-sender.  It must be FLAT.  If it is not, there is a
spatial technical gradient in the assay itself and every biological kernel in the
paper is contaminated by it.  A7 is a go/no-go for H1, and the mouse half needs no
download: `data/raw/<section>/cells.parquet` already carries the per-cell control
tallies that `cell_feature_matrix.h5` breaks out as 40 `Negative Control Probe`,
609 `Negative Control Codeword` and 21 `Genomic Control` features.

METHOD -- the whole point is that NOTHING about the estimator changes.  This
script reuses `run_phase3_nulls.SectionFit` and `run_phase3_nulls.fit_cell`
verbatim; the only substitution is the response matrix Y (and hence the N6
neighbour baseline, which is recomputed from it by the same
`phase3_core.neighbour_baseline`).  Same 100 um window, same 40-point lambda grid,
same MIN_RECEIVERS=2000 floor, same nested designs, same 400-replicate spatial
block bootstrap, same sender calls.

RESPONSES (each z-scored over the section's analysis cells, so beta IS the
amplitude in response-SD units and is directly comparable to a module's
beta/sd_y):
  neg_control_probe     control_probe_counts       (40 Negative Control Probes)
  neg_control_codeword  control_codeword_counts    (609 Negative Control Codewords)
  genomic_control       genomic_control_counts     (21 Genomic Controls)
  all_controls          the sum of the three
  neg_probe_rate        control_probe_counts / (transcript_counts + 1)

Usage: run_a7_control_probes.py [--sections all|inband] [--calls a,b] [--n-jobs N]
"""
import os, sys, time, argparse, numpy as np, pandas as pd
sys.path.insert(0, '/workspace/code')
from joblib import Parallel, delayed
import sasp_phase3 as P
import phase3_core as C
import run_phase3_nulls as R

RES = '/workspace/results/phase3'
RAW = '/workspace/data/raw'

RESPONSES = ['neg_control_probe', 'neg_control_codeword', 'genomic_control',
             'all_controls', 'neg_probe_rate']
COLS = {'neg_control_probe': 'control_probe_counts',
        'neg_control_codeword': 'control_codeword_counts',
        'genomic_control': 'genomic_control_counts'}


def control_matrix(sec):
    """Per-cell control tallies for the cache's analysis cells, in cache order."""
    c = pd.read_parquet(os.path.join(RAW, sec.name, 'cells.parquet')).set_index('cell_id')
    c = c.reindex(pd.Index(sec.cell_id.astype(str)))
    assert c[COLS['neg_control_probe']].notna().all(), 'cache cell_id not in cells.parquet'
    raw = {k: c[v].to_numpy(float) for k, v in COLS.items()}
    raw['all_controls'] = sum(raw[k] for k in COLS)
    raw['neg_probe_rate'] = (raw['neg_control_probe']
                             / (c['transcript_counts'].to_numpy(float) + 1.0))
    Y = np.column_stack([raw[k] for k in RESPONSES])
    prov = [dict(section=sec.name, response=k, n_cells=len(c),
                 mean_per_cell=float(np.mean(raw[k])),
                 sd_per_cell=float(np.std(raw[k])),
                 frac_cells_nonzero=float(np.mean(raw[k] > 0)),
                 max_per_cell=float(np.max(raw[k])))
            for k in RESPONSES]
    # z-score so beta is an amplitude in response-SD units
    sd = Y.std(0); sd[sd < 1e-12] = 1.0
    return (Y - Y.mean(0)) / sd, pd.DataFrame(prov)


BINS = np.arange(0.0, R.WINDOW_UM + 1e-9, 5.0)


def _curves(sf, celltype='Hepatocytes'):
    """Binned control response vs distance-to-nearest-sender, raw and after
    residualising on the SAME full N6+N5 design the fits condition on.  This is
    what panel h of Figure 2 draws; it is the visual form of "must be flat"."""
    idx = sf.receivers(celltype)
    if idx.sum() < R.MIN_RECEIVERS:
        return pd.DataFrame()
    d = sf.d_obs[idx]
    k = np.clip(np.digitize(d, BINS) - 1, 0, len(BINS) - 2)
    rows = []
    for j, nm in enumerate(RESPONSES):
        y = sf.Y[idx, j].astype(float)
        X1, _, _, _ = R._designs(sf, idx, j)
        resid = P.FixedLambdaFitter(X1, y[:, None]).Y[:, 0]
        for b in range(len(BINS) - 1):
            m = k == b
            if m.sum() < 50:
                continue
            rows.append(dict(section=sf.sec.name, arm=sf.sec.meta['condition'],
                             call=sf.call, response=nm, celltype=celltype,
                             bin_lo=BINS[b], bin_hi=BINS[b + 1],
                             bin_mid=0.5 * (BINS[b] + BINS[b + 1]),
                             n=int(m.sum()),
                             mean_raw=float(y[m].mean()),
                             sem_raw=float(y[m].std(ddof=1) / np.sqrt(m.sum())),
                             mean_resid=float(resid[m].mean()),
                             sem_resid=float(resid[m].std(ddof=1) / np.sqrt(m.sum()))))
    return pd.DataFrame(rows)


def _job(sample, call, seed):
    t0 = time.time()
    sf = R.SectionFit(sample, call, seed)
    Y, prov = control_matrix(sf.sec)
    sf.Y = Y
    sf.NB = C.neighbour_baseline(sf.sec, Y, sf.sender)
    old = P.MODULES
    P.MODULES = RESPONSES                     # fit_cell reads P.MODULES[j] for the label
    try:
        rows = []
        types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
        for t in types:
            if sf.receivers(t).sum() < R.MIN_RECEIVERS:
                continue
            for j in range(len(RESPONSES)):
                rows.append(R.fit_cell(sf, t, j, seed + 7 * j, tag='A7'))
        cur = _curves(sf)
    finally:
        P.MODULES = old
    print('[a7] %s %s %d rows %.0fs' % (sample, call, len(rows), time.time() - t0), flush=True)
    return rows, prov, cur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sections', default='all')
    ap.add_argument('--calls', default='tierA_p95,cdkn1a_pos')
    ap.add_argument('--n-jobs', type=int, default=11)
    a = ap.parse_args()
    secs = P.ALL_SECTIONS if a.sections == 'all' else P.IN_BAND
    calls = a.calls.split(',')
    jobs = [(s, c, P.MASTER_SEED + 1000 * i + j)
            for i, s in enumerate(secs) for j, c in enumerate(calls)]
    out = Parallel(n_jobs=a.n_jobs, prefer='processes', verbose=5)(
        delayed(_job)(s, c, sd) for s, c, sd in jobs)
    df = pd.DataFrame([r for rs, _, _ in out for r in rs])
    df = df.rename(columns={'module': 'response'})
    df.to_csv(RES + '/a7_control_probe_fits.csv', index=False)
    pv = pd.concat([p for _, p, _ in out]).drop_duplicates(['section', 'response'])
    pv.to_csv(RES + '/a7_control_probe_provenance.csv', index=False)
    cu = pd.concat([c for _, _, c in out if len(c)])
    cu.to_csv(RES + '/a7_control_probe_curves.csv', index=False)
    print(cu.shape, '->', RES + '/a7_control_probe_curves.csv')
    print(df.shape, '->', RES + '/a7_control_probe_fits.csv')


if __name__ == '__main__':
    main()
