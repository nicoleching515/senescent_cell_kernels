#!/usr/bin/env python3
"""Phase 8, task 8.5 / C7-D2 -- seed-stability panel.

Four runs on ONE fixed 20,000-cell subsample of 7239 (subsampling seed 12345, independent
of the model seed, so all four see the same cells):

  raw seed0 vs raw seed1   DeepScence's own run-to-run spread at denoise=False
  dca seed0/1/2 pairwise   the spread with DCA in the loop (DeepScence hands its
                           random_state straight to dca(), so one seed moves both).
                           A third seed was added because the seed0-vs-seed1 top-5%
                           Jaccard came back at exactly 0 and a claim that strong should
                           not rest on one pair.
  raw seed0 vs dca seed0   the denoise effect, on the identical cells, at the same seed

Same metrics as analyse_d2_denoise.py so the numbers are directly comparable: Pearson,
Spearman, and Jaccard of the global top-5% call set.  Writes
results/phase8_d2/d2_stability.csv.
"""
import os, itertools
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

PROC = '/workspace/data/processed/'
OUT = '/workspace/results/phase8_d2/'
SEC = '7239_liver_sbr_Male_52-U1'
N = 20000
RUNS = {'raw_seed0': 'raw_sub%d' % N, 'raw_seed1': 'raw_sub%d_seed1' % N,
        'dca_seed0': 'dca_sub%d' % N, 'dca_seed1': 'dca_sub%d_seed1' % N,
        'dca_seed2': 'dca_sub%d_seed2' % N}


def main():
    v = {}
    for k, tag in RUNS.items():
        p = PROC + 'deepscence_%s_%s.csv' % (tag, SEC)
        if os.path.exists(p):
            v[k] = pd.read_csv(p).set_index('cell_id')['deepscence_score']
    rows = []
    for a, b in itertools.combinations(sorted(v), 2):
        j = pd.concat([v[a], v[b]], axis=1, join='inner').dropna()
        x = j.iloc[:, 0].to_numpy(float); y = j.iloc[:, 1].to_numpy(float)
        r = float(pearsonr(x, y).statistic)
        yy = y * (-1 if r < 0 else 1)
        fa = x > np.percentile(x, 95); fb = yy > np.percentile(yy, 95)
        rows.append(dict(pair='%s vs %s' % (a, b), n_cells=int(len(j)),
                         pearson_r=round(r, 5),
                         spearman_rho=round(float(spearmanr(x, yy).statistic), 5),
                         anchor_sign_flipped=bool(r < 0),
                         top5_jaccard=round(float((fa & fb).sum() / max((fa | fb).sum(), 1)), 4),
                         top5_n_changed=int((fa ^ fb).sum())))
    d = pd.DataFrame(rows)
    d.to_csv(OUT + 'd2_stability.csv', index=False)
    print(d.to_string(index=False))
    print('\nwrote ' + OUT + 'd2_stability.csv')


if __name__ == '__main__':
    main()
