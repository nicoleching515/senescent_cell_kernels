#!/usr/bin/env python3
"""Phase 8, task 8.5 / C7-D2 -- how far does the DeepScence score move when the depth
normalisation changes, and is the movement depth-dependent?

Compares each alternative configuration produced by run_deepscence_denoise_probe.py /
run_deepscence_dca.py against the COMMITTED `denoise=False` scores, cell by cell, on the
same cell set:

  raw  vs committed  -- determinism control.  Same code, same seed, same data: whatever
                        this shows is the run-to-run noise floor of a 300-epoch CPU
                        autoencoder, and no other comparison may be read below it.
  mor  vs committed  -- median-of-ratios depth normalisation before scoring (§6).
  lib  vs committed  -- exact library-size equalisation before scoring.
  ds10 vs committed  -- binomial downsampling of every cell to the 10th-percentile depth.
  dca  vs committed  -- the published default, denoise=True, real DCA 0.3.4.

DeepScence's score is a bottleneck pre-activation: its scale is arbitrary and its sign is
set per run by the CDKN1A anchor.  Raw Pearson r is therefore reported (a large negative r
means the anchor flipped, which is itself a result), but all difference metrics are
computed on z-scored, sign-aligned vectors so that "the score moved" means "the score
moved relative to its own spread", not "the run happened to pick the other sign".

Sender status is re-thresholded exactly as caller_disagree.py does it -- global top-5%,
and within cell type, and within (cell type x within-type depth decile) -- on the cell set
that pipeline uses (cell types joined, Low_quality/Unknown dropped).  Global top-1% is
added because Test 3's prevalence floor lives down there.

Writes results/phase8_d2/d2_{agreement,depth,deciles}.csv.
"""
import os, sys, glob, json
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

PROC = '/workspace/data/processed/'
RAW = '/workspace/data/raw/'
OUT = '/workspace/results/phase8_d2/'
# the two sections whose committed DeepScence scores are stored under the preserved
# original names (caller_disagree.py DS_ALIAS)
ALIAS = {'7250_liver_sham_Male_26-U1': 'sham', '7259_liver_sbr_Male_26-U1': 'sbr'}
EXCL = {'Low_quality', 'Unknown'}
MIN_STRATUM = 50
CONFIGS = ['raw', 'mor', 'lib', 'ds10', 'dca']


def baseline(section):
    return pd.read_csv(PROC + 'deepscence_%s.csv' % ALIAS.get(section, section),
                       usecols=['cell_id', 'deepscence_score']).set_index('cell_id')


def context(section):
    """cell type + transcript_counts, on the caller pipeline's cell set."""
    ct = pd.read_csv(PROC + 'celltypes_%s.csv' % ALIAS.get(section, section),
                     usecols=['cell_id', 'cell_type']).set_index('cell_id')
    cells = pd.read_parquet(RAW + section + '/cells.parquet',
                            columns=['cell_id', 'transcript_counts']).set_index('cell_id')
    df = ct.join(cells)
    return df[~df.cell_type.isin(EXCL)]


def top_flags(v, ct, tc, pct, mode):
    """mode: 'global' | 'within_type' | 'matched' (cell type x within-type depth decile)."""
    ok = np.isfinite(v)
    f = np.zeros(len(v), bool)
    if mode == 'global':
        f[ok] = v[ok] > np.nanpercentile(v[ok], 100 - pct)
        return f
    for c in pd.unique(ct):
        m0 = (ct == c) & ok
        if m0.sum() < MIN_STRATUM:
            continue
        if mode == 'within_type':
            f[m0] = v[m0] > np.nanpercentile(v[m0], 100 - pct)
            continue
        idx = np.where(m0)[0]
        try:
            d = pd.qcut(tc[idx], 10, labels=False, duplicates='drop')
        except Exception:
            d = np.zeros(len(idx), int)
        for k in np.unique(d):
            ii = idx[d == k]
            if len(ii) < MIN_STRATUM:
                continue
            f[ii] = v[ii] > np.nanpercentile(v[ii], 100 - pct)
    return f


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std()


def main():
    agree, depth, deciles, strength, comp = [], [], [], [], []
    for cfg in CONFIGS:
        for path in sorted(glob.glob(PROC + 'deepscence_%s_7*.csv' % cfg)):
            section = os.path.basename(path)[len('deepscence_%s_' % cfg):-4]
            if not os.path.isdir(RAW + section):
                continue                      # subsample probes etc.
            alt = pd.read_csv(path).set_index('cell_id')
            # how much depth variation did this configuration actually remove?  Without
            # this, "the score barely moved" cannot be distinguished from "the
            # normalisation barely did anything".
            if 'size_factor' in alt.columns and 'mapped_counts' in alt.columns:
                lc = np.log(alt.mapped_counts.to_numpy(float))
                lr = lc - np.log(alt.size_factor.to_numpy(float))
                strength.append(dict(
                    section=section, config=cfg,
                    sd_log_depth_before=round(float(lc.std()), 4),
                    sd_log_depth_after=round(float(lr.std()), 4),
                    frac_log_depth_variance_removed=round(
                        float(1 - lr.var() / lc.var()), 4),
                    depth_p90_over_p10_before=round(float(np.exp(
                        np.percentile(lc, 90) - np.percentile(lc, 10))), 3),
                    depth_p90_over_p10_after=round(float(np.exp(
                        np.percentile(lr, 90) - np.percentile(lr, 10))), 3)))
            base = baseline(section)
            ctx = context(section)
            j = base.join(alt, rsuffix='_alt', how='inner').join(ctx, how='inner').dropna(
                subset=['deepscence_score', 'deepscence_score_alt', 'cell_type',
                        'transcript_counts'])
            a = j.deepscence_score.to_numpy(float)
            b = j.deepscence_score_alt.to_numpy(float)
            tc = j.transcript_counts.to_numpy(float)
            ct = j.cell_type.to_numpy()
            r = float(pearsonr(a, b).statistic)
            rho = float(spearmanr(a, b).statistic)
            flip = bool(r < 0)
            # Sign alignment is only meaningful when the two vectors are actually
            # associated.  DeepScence's sign is set per run by the CDKN1A anchor, so a
            # large negative r means "the anchor flipped" and negating is right; an r near
            # zero means the runs found near-orthogonal axes, and which sign we negate is
            # arbitrary.  Flag it rather than let a coin-flip propagate into the Jaccards.
            ambiguous = bool(abs(r) < 0.1)
            za, zb = z(a), z(b * (-1 if flip else 1))
            d = zb - za
            row = dict(section=section, config=cfg, n_cells=int(len(j)),
                       # 8 dp, not 6: the determinism control lands at 0.99999978 and
                       # rounding it to 1.0 would erase the point of the control
                       pearson_r=round(r, 8), spearman_rho=round(rho, 8),
                       anchor_sign_flipped=flip,
                       sign_alignment_ambiguous=ambiguous,
                       abs_pearson_r=round(abs(r), 8),
                       mean_abs_dz=round(float(np.abs(d).mean()), 4),
                       p95_abs_dz=round(float(np.percentile(np.abs(d), 95)), 4))
            for mode in ('global', 'within_type', 'matched'):
                for pct in ((5, 1) if mode == 'global' else (5,)):
                    fa = top_flags(a, ct, tc, pct, mode)
                    fb = top_flags(b * (-1 if flip else 1), ct, tc, pct, mode)
                    k = '%s_top%d' % (mode, pct)
                    inter = int((fa & fb).sum()); union = int((fa | fb).sum())
                    row['%s_n_called' % k] = int(fa.sum())
                    row['%s_jaccard' % k] = round(inter / max(union, 1), 4)
                    row['%s_n_changed' % k] = int((fa ^ fb).sum())
                    row['%s_pct_of_called_changed' % k] = round(
                        100 * (fa.sum() - inter) / max(fa.sum(), 1), 2)
            agree.append(row)

            # ---- depth dependence ----
            rd_a = float(spearmanr(a, tc).statistic)
            rd_b = float(spearmanr(b * (-1 if flip else 1), tc).statistic)
            fa5 = top_flags(a, ct, tc, 5, 'global')
            fb5 = top_flags(b * (-1 if flip else 1), ct, tc, 5, 'global')
            depth.append(dict(
                section=section, config=cfg, n_cells=int(len(j)),
                median_transcript_counts=float(np.median(tc)),
                rho_depth_committed=round(rd_a, 4), rho_depth_alt=round(rd_b, 4),
                delta_rho_depth=round(rd_b - rd_a, 4),
                rho_signed_dz_vs_depth=round(float(spearmanr(d, tc).statistic), 4),
                rho_abs_dz_vs_depth=round(float(spearmanr(np.abs(d), tc).statistic), 4),
                pct_cells_flipping_top5=round(100 * float((fa5 ^ fb5).mean()), 4)))
            # cell-type composition of the top-5% calls, committed vs alternative.  The
            # project's standing claim about DeepScence is that its call set is
            # depth-dominated rather than type-driven; if the denoising step changes which
            # types get called, that claim is configuration-dependent and must say so.
            bg = pd.Series(ct).value_counts(normalize=True)
            for c in bg.index:
                mc = ct == c
                comp.append(dict(section=section, config=cfg, cell_type=c,
                                 bg_pct=round(100 * float(bg[c]), 3),
                                 call_pct_committed=round(100 * float(fa5[mc].sum()) / max(fa5.sum(), 1), 3),
                                 call_pct_alt=round(100 * float(fb5[mc].sum()) / max(fb5.sum(), 1), 3),
                                 enrich_committed=round(float(fa5[mc].sum()) / max(fa5.sum(), 1) / float(bg[c]), 3),
                                 enrich_alt=round(float(fb5[mc].sum()) / max(fb5.sum(), 1) / float(bg[c]), 3)))
            q = pd.qcut(tc, 10, labels=False, duplicates='drop')
            for k in np.unique(q):
                m = q == k
                deciles.append(dict(section=section, config=cfg, depth_decile=int(k) + 1,
                                    n=int(m.sum()),
                                    median_counts=float(np.median(tc[m])),
                                    mean_dz=round(float(d[m].mean()), 4),
                                    mean_abs_dz=round(float(np.abs(d[m]).mean()), 4),
                                    pct_called_committed=round(100 * float(fa5[m].mean()), 3),
                                    pct_called_alt=round(100 * float(fb5[m].mean()), 3),
                                    pct_flipped=round(100 * float((fa5[m] ^ fb5[m]).mean()), 3)))
            print('%-4s %-30s n=%7d r=%+.5f  J(top5)=%.4f  drho_depth=%+.4f'
                  % (cfg, section, len(j), r, row['global_top5_jaccard'],
                     rd_b - rd_a), flush=True)
    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame(agree).to_csv(OUT + 'd2_agreement.csv', index=False)
    pd.DataFrame(depth).to_csv(OUT + 'd2_depth.csv', index=False)
    pd.DataFrame(deciles).to_csv(OUT + 'd2_deciles.csv', index=False)
    pd.DataFrame(strength).to_csv(OUT + 'd2_normalisation_strength.csv', index=False)
    pd.DataFrame(comp).to_csv(OUT + 'd2_celltype_composition.csv', index=False)
    print('wrote %sd2_{agreement,depth,deciles,normalisation_strength,celltype_composition}.csv' % OUT)


if __name__ == '__main__':
    main()
