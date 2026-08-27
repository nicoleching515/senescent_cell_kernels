#!/usr/bin/env python3
"""Phase 8, task 8.5 / C7-D2 -- what does `denoise=False` cost?

Phase 7 addendum §6 path 3: DCA does not install (see reports/CS_PHASE8_D2_DENOISE.md
§3 for the attempt), so the cost of running DeepScence WITHOUT its denoising step is
measured directly: score the SAME cells with and without an explicit depth
normalisation applied before scoring, and report how far the score moves.

Why this file rather than an edit to run_deepscence_all.py:
  the committed DeepScence scores (data/processed/deepscence_*.csv) are inputs to a
  completed analysis (the 11-section coverage gate, task 8.3/8.4).  They must survive
  byte-identical.  This script writes NEW files under `deepscence_<config>_<section>.csv`
  and never touches the existing ones.

Three configurations, all `denoise=False`, `random_state=0`, published CDKN1A anchor,
MGI 1:1 ortholog-remapped panel, >=20 raw counts/cell:

  raw  -- byte-for-byte the committed configuration.  Run again on a few sections as a
          DETERMINISM CONTROL: without it, a score shift under `mor` cannot be
          separated from ordinary run-to-run drift of a 300-epoch CPU autoencoder.
  mor  -- median-of-ratios (DESeq2 "poscounts") size factors estimated on the
          ortholog-mapped panel and DIVIDED OUT of the counts before DeepScence sees
          them.  §6 names this estimator.  This is the "equivalent depth normalisation
          of our own".
  lib  -- per-cell library-size factors (mapped counts / median mapped counts) divided out
          before scoring.  After it every cell has exactly the median mapped depth, so
          depth is removed COMPLETELY and DeepScence's own size-factor offset collapses to
          1.  `mor` is the estimator §6 names; `lib` is the upper bound on what any
          rescaling normalisation can do, and the two together bracket the answer.
  ds10 -- binomial downsampling of every cell to a common depth (the 10th percentile of
          mapped counts).  The extreme depth control: after it, library size carries no
          information at all.  Bounds how far any depth normalisation could move the
          score.

WHAT `denoise=False` DOES AND DOES NOT SWITCH OFF.  Read DeepScence/io.py:normalize --
the published pipeline ALREADY calls sc.pp.normalize_total on the full panel, then
log1p, then scale, and it ALREADY passes per-cell size factors into the ZINB decoder as
an offset (network.py: `mu = mu * sf`).  So `denoise=False` does not remove depth
normalisation; it removes DCA's ZINB *imputation*.  `mor` therefore replaces
library-size normalisation with a composition-robust one rather than adding
normalisation where there was none, and `ds10` removes the depth signal outright.

Cell and gene sets are identical across configurations by construction: the >=20-count
cell filter and the ortholog map are applied to RAW counts before any normalisation, and
scaling a count matrix by a per-cell constant leaves the support (and hence
sc.pp.filter_genes(min_cells=1) inside DeepScence) unchanged.  Scores are therefore
comparable cell by cell.

Usage: run_deepscence_denoise_probe.py --config {raw,mor,ds10} <section_dir_name> [...]
"""
import sys, os, time, csv, json, argparse
import numpy as np, pandas as pd, anndata as ad, h5py
from scipy.sparse import csc_matrix, csr_matrix

import torch
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "16")))

sys.path.insert(0, '/workspace/code/_shims')
import DeepScence.api as api

RAW = '/workspace/data/raw/'
PROC = '/workspace/data/processed/'
META = '/workspace/results/phase8_d2/'
ORTH = '/workspace/genesets/mouse_human_orthologs_MGI.csv'

# genes used to ESTIMATE the median-of-ratios size factor.  DESeq2's estimator wants a
# gene's geometric mean to be meaningful; on a sparse Xenium panel a gene detected in a
# handful of cells gives a geometric mean driven by noise.  Genes detected in fewer than
# this fraction of cells are excluded FROM THE ESTIMATOR ONLY -- every gene is still
# normalised and still scored.
MOR_MIN_DETECT_FRAC = 0.05
MIN_COUNTS = 20          # identical to run_deepscence.py / run_deepscence_all.py
DS10_PCT = 10            # ds10 target depth = this percentile of mapped counts


def load_section(section):
    """Identical to run_deepscence_all.py's loader, up to and including the ortholog map."""
    f = h5py.File(RAW + section + '/cell_feature_matrix.h5', 'r')
    ft = np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
    ids = np.array([x.decode() for x in f['matrix/features/id'][:]])
    nm = np.array([x.decode() for x in f['matrix/features/name'][:]])
    bc = np.array([x.decode() for x in f['matrix/barcodes'][:]])
    keep = (ft == 'Gene Expression') & np.char.startswith(ids, 'ENSMUSG')
    M = csc_matrix((f['matrix/data'][:].astype(np.float32),
                    f['matrix/indices'][:].astype(np.int32),
                    f['matrix/indptr'][:].astype(np.int64)),
                   shape=tuple(f['matrix/shape'][:]))
    f.close()
    A = ad.AnnData(csr_matrix(M[np.where(keep)[0], :].T),
                   obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=nm[keep]))
    del M
    A = A[np.asarray(A.X.sum(1)).ravel() >= MIN_COUNTS].copy()
    orth = {r['mouse_symbol']: r['human_symbol'] for r in csv.DictReader(open(ORTH))}
    hs = np.array([orth.get(g, '') for g in A.var_names])
    A = A[:, hs != ''].copy()
    A.var_names = hs[hs != '']
    A.var_names_make_unique()
    return A


def mor_size_factors(X):
    """DESeq2 median-of-ratios, 'poscounts' variant, on a CSR count matrix.

    gm_g   = exp(mean over cells where x_ig > 0 of log x_ig)     (geometric mean, positives only)
    sf_i   = exp(median over genes g in G, x_ig > 0, of (log x_ig - log gm_g))
    then rescaled so median_i(sf_i) == 1, which keeps the matrix on its original count
    scale (only the per-cell tilt is removed, not the overall magnitude).

    G = genes detected in >= MOR_MIN_DETECT_FRAC of cells.
    """
    X = X.tocsr()
    n_cells = X.shape[0]
    logd = np.log(X.data.astype(np.float64))
    cols = X.indices
    n_pos = np.bincount(cols, minlength=X.shape[1]).astype(np.float64)
    sum_log = np.bincount(cols, weights=logd, minlength=X.shape[1])
    detect_frac = n_pos / n_cells
    use_gene = detect_frac >= MOR_MIN_DETECT_FRAC
    log_gm = np.full(X.shape[1], np.nan)
    log_gm[use_gene] = sum_log[use_gene] / n_pos[use_gene]

    ratio = logd - log_gm[cols]            # nan where the gene is not in G
    sf = np.empty(n_cells, np.float64)
    indptr = X.indptr
    for i in range(n_cells):
        r = ratio[indptr[i]:indptr[i + 1]]
        r = r[np.isfinite(r)]
        sf[i] = np.median(r) if r.size else np.nan
    sf = np.exp(sf)
    # a cell with no gene in G (possible only for very shallow cells) falls back to its
    # library-size factor, so it is still normalised rather than dropped.
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    fallback = ~np.isfinite(sf)
    if fallback.any():
        sf[fallback] = lib[fallback] / np.median(lib)
        sf[fallback & (sf <= 0)] = 1.0
    sf = sf / np.median(sf)
    return sf, int(use_gene.sum()), int(fallback.sum())


def downsample_to(X, target, rng):
    """Binomial thinning of every cell to `target` expected mapped counts.

    Cells already at or below `target` are left alone (thinning cannot add counts); the
    fraction so left is reported, because it is the honest limit of this control.
    """
    X = X.tocsr().copy()
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    p = np.clip(target / np.maximum(lib, 1e-9), 0.0, 1.0)
    indptr = X.indptr
    prow = np.repeat(p, np.diff(indptr))
    X.data = rng.binomial(X.data.astype(np.int64), prow).astype(np.float32)
    X.eliminate_zeros()
    return X, float((lib <= target).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True, choices=['raw', 'mor', 'lib', 'ds10'])
    ap.add_argument('--subsample', type=int, default=0,
                    help='score only N cells, chosen with a fixed seed independent of '
                         '--seed so every seed sees the SAME cells; for the seed-stability '
                         'panel, which has to fit alongside another agent in a 57.7 GiB cgroup')
    ap.add_argument('--seed', type=int, default=0,
                    help='DeepScence random_state.  0 = the committed setting.')
    ap.add_argument('sections', nargs='+')
    a = ap.parse_args()

    tag = a.config + ('' if not a.subsample else '_sub%d' % a.subsample) \
                   + ('' if a.seed == 0 else '_seed%d' % a.seed)
    for section in a.sections:
        out = PROC + 'deepscence_%s_%s.csv' % (tag, section)
        if os.path.exists(out):
            print('%s [%s]: already present, skipping' % (section, tag), flush=True)
            continue
        t0 = time.time()
        A = load_section(section)
        if a.subsample and a.subsample < A.n_obs:
            idx = np.sort(np.random.default_rng(12345).choice(A.n_obs, a.subsample,
                                                              replace=False))
            A = A[idx].copy()
        mapped_counts = np.asarray(A.X.sum(1)).ravel().astype(np.float64)
        meta = dict(section=section, config=tag, n_cells=int(A.n_obs),
                    n_genes=int(A.n_vars), denoise=False, random_state=a.seed,
                    anchor='published_CDKN1A', panel='ortholog_mapped_MGI_1to1',
                    min_counts_per_cell=MIN_COUNTS)
        sf = np.ones(A.n_obs)
        if a.config == 'mor':
            sf, n_gene_used, n_fallback = mor_size_factors(A.X)
            meta.update(mor_min_detect_frac=MOR_MIN_DETECT_FRAC,
                        mor_genes_in_estimator=n_gene_used,
                        mor_cells_fallback_to_libsize=n_fallback,
                        mor_sf_p01=float(np.percentile(sf, 1)),
                        mor_sf_p99=float(np.percentile(sf, 99)),
                        mor_sf_spearman_vs_mapped_counts=None)
            from scipy.stats import spearmanr
            meta['mor_sf_spearman_vs_mapped_counts'] = round(
                float(spearmanr(sf, mapped_counts).statistic), 4)
            Xn = A.X.tocsr().copy()
            Xn.data = (Xn.data / np.repeat(sf, np.diff(Xn.indptr))).astype(np.float32)
            A.X = Xn
            del Xn
        elif a.config == 'lib':
            sf = mapped_counts / np.median(mapped_counts)
            meta.update(lib_sf_p01=float(np.percentile(sf, 1)),
                        lib_sf_p99=float(np.percentile(sf, 99)))
            Xn = A.X.tocsr().copy()
            Xn.data = (Xn.data / np.repeat(sf, np.diff(Xn.indptr))).astype(np.float32)
            A.X = Xn
            del Xn
        elif a.config == 'ds10':
            target = float(np.percentile(mapped_counts, DS10_PCT))
            rng = np.random.default_rng(a.seed)
            Xd, frac_at_or_below = downsample_to(A.X, target, rng)
            A.X = Xd
            del Xd
            meta.update(ds10_target_counts=target,
                        ds10_frac_cells_left_untouched=round(frac_at_or_below, 4),
                        ds10_mean_counts_after=float(np.asarray(A.X.sum(1)).ravel().mean()))
        print('%s [%s]: %d cells x %d ortholog-mapped genes' %
              (section, a.config, A.n_obs, A.n_vars), flush=True)
        ids_kept = A.obs_names.to_numpy()
        t = time.time()
        res = api.DeepScence(A, denoise=False, verbose=False, random_state=a.seed)
        mins = (time.time() - t) / 60
        log = dict(res.uns['log'])
        meta['direction_log'] = {k: (v.tolist() if hasattr(v, 'tolist') else v)
                                 for k, v in log.items()}
        meta['minutes'] = round(mins, 2)
        meta['total_minutes'] = round((time.time() - t0) / 60, 2)
        pd.DataFrame({'cell_id': ids_kept,
                      'deepscence_score': np.round(res.obs['ds'].to_numpy(), 5),
                      'mapped_counts': mapped_counts.astype(np.int64),
                      'size_factor': np.round(sf, 6)}).to_csv(out, index=False)
        with open(META + 'runmeta_%s_%s.json' % (tag, section), 'w') as fh:
            json.dump(meta, fh, indent=1, default=str)
        print('wrote %s (%d cells, %.1f min)' % (os.path.basename(out), len(ids_kept), mins),
              flush=True)


if __name__ == '__main__':
    main()
