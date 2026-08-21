# SASP Spatial Response Kernel — Phase 1 (synthetic identifiability study)

Master Plan Section 22, Step 1. See `/workspace/reports/CS_PHASE1.md` for results.

## Modules
| file | what |
|---|---|
| `sasp_sim.py` | synthetic tissue generator: hard-core point process, GRF nuisance fields, controlled sender clustering, planted exponential kernel, sender-density-correlated confounder |
| `sasp_estimators.py` | naive / binned / nuisance-conditioned (N5) / matched-decoy (N2) estimators, propensity matching + SMD, `BlockProfiler` block bootstrap |
| `sasp_sweep.py` | seed-pinned, checkpointed sweep runner (joblib) |
| `make_curves.py` | representative binned curves for Figure 1c |
| `make_figure1.py` | Figure 1 + `figure1_data.csv` |
| `summarize.py` | the result tables in `/workspace/results/summary_tables.txt` |
| `sasp_palette.py` | validated chart palette |
| `diag.py`, `smoke.py` | small end-to-end diagnostics used to develop and sanity-check the pipeline |

## Reproduce everything
```bash
cd /workspace/code
python3 -u sasp_sweep.py --n-jobs 46      # ~9 min on 46 cores; resumes from checkpoints
python3 -u make_curves.py
python3 -u summarize.py
python3 -u make_figure1.py
```
Delete `/workspace/results/sweep/*.csv` first to force a full recompute.
Everything is seeded from `MASTER_SEED = 20260820` in `sasp_sweep.py`; repeated runs
are bit-identical.

## Gotchas
* Pin BLAS to one thread (`OMP_NUM_THREADS=1`); the fits use many small matmuls and
  multithreaded BLAS costs ~10x from thread contention. `sasp_sweep.py` sets this.
* `pandas.read_csv` parses the literal string `"null"` as `NaN`, which silently drops
  the `beta_true = 0` sweep. Always `df["sweep"] = df["sweep"].fillna("null")`.

## Phase 2 additions

| file | what |
|---|---|
| `sasp_real.py` | Xenium loader; Gene Expression features only; swappable sender/celltype/module/anatomy inputs |
| `prepare_samples.py` | one-time per-section cache (`/workspace/data/processed/cache/*.npz`); `Cached` accessor |
| `sasp_kernels.py` | all five Section 6.2 families, profiled fits, AIC, iid + spatial block bootstrap CIs, comparable `d_half` |
| `run_figure2a.py` | real-data fit grid (section x module x family), checkpointed |
| `make_figure2.py` | Figure 2a |
| `sasp_nulls.py` | nulls N1 (within-type label permutation), N3 (torus shift), N4 (rotation) |
| `sasp_phase1b.py` | kernel-family misspecification, superposition vs nearest, null size/power |

### Reproduce Phase 2
```bash
cd /workspace/code
python3 -u prepare_samples.py --n-jobs 4          # re-run with --force after Bio annotations land
python3 -u run_figure2a.py --sender tierA_p95 --n-jobs 24
python3 -u make_figure2.py
python3 -u sasp_phase1b.py --n-jobs 44
```

### Gotchas (Phase 2)
* Do NOT oversubscribe: with BLAS pinned to 1 thread, keep total workers <= cores.
  Running two 40-worker jobs at once dropped each worker to 17 % CPU.
* `pgrep -f <pattern>` inside a wait loop matches the wait loop's own command line.
  Two launches were lost to this before it was spotted.
* Spline knots are data-quantile based and MUST be carried from the fit to any later
  evaluation; recomputing them on a plotting grid silently yields a different basis.
* When a permutation null re-profiles lambda it compares different models; hold lambda
  at the observed value (see `sasp_nulls.run_nulls`, `*_fixedlam` outputs).


## Phase 3 additions (null battery N1-N8 on real tissue)

| file | what |
|---|---|
| `sasp_phase3.py` | Phase 3 per-section cache (`data/processed/cache3/*.npz`), `Sec` accessor, sender calls, Section 8 Test 3 section-admissibility lists (`IN_BAND` / `OVER_CEILING` / `BELOW_FLOOR`), merged receiver labels |
| `phase3_core.py` | matched decoys on real sections (same greedy propensity algorithm as `sasp_estimators.match_decoys`), Tier D covariate blocks, N6 neighbour baseline |
| `run_phase3_nulls.py` | stages `window` / `main` (N2,N5,N6,zonation + block bootstrap) / `perm` (N1,N3,N4) / `curves` |
| `run_phase3_n8.py` | N8: Tier A x Tier B disjointness, Tier E3 scrambled-response control, CoreScence circularity |
| `run_phase3_strat.py` | unstratified vs cell-type-conditioned vs N5 decomposition of the naive gradient |
| `run_phase3_attribution.py` | which N5 sub-block removes the amplitude |
| `run_phase3_combined.py` | T2 pooled estimate with a DONOR-level bootstrap (animal = block of `BlockProfiler`) |
| `run_phase3_poisson.py`, `run_phase3_lamscale.py` | is the regressor a measurement or a sender-calling rate |
| `_correlogram.py`, `_ripley.py` | where the real sections sit on the Figure 1 regime map |
| `summarize_phase3.py`, `make_figure2bc.py` | summary tables, Figures 2b/2c/2d |

Two methods were ADDED to `sasp_estimators.BlockProfiler` (additive, existing
behaviour unchanged): `beta_at(m, p, t)` and `beta2_at(m, p, t)` read beta at a
FIXED lambda grid index, which is what every null needs.

### Gotchas (Phase 3)
* `np.isin` on a numpy-unicode array against a pandas object Index falls back to an
  O(n*m) path: 15 min per section. Use `pd.Index(a).isin(b)`.
* `cmd 2>&1 | tail` in a background launcher buffers everything until exit.
* `pkill -f <pattern>` matches the launching shell itself (recorded in Phase 2, and
  it happened again).
* Pooling sections needs a CANONICAL cell-type list for the kNN-composition and
  segmentation covariate blocks, or the designs have different widths per section.

---

## Phase 4 — existing spatial CCC tools under the same nulls (Figure 4)

Report: `reports/CS_PHASE4.md`. Results: `results/phase4/`.

| file | what |
|---|---|
| `phase4_data.py` | LR-gene cache over the Phase 3 cell set (`/tmp/p4/cache4`), the four BIO_PHASE3-certified LR pairs, and the four coordinate nulls (`torus_shift` / `rotate_about_centroid` copied verbatim from `run_phase3_nulls.py`) |
| `phase4_tiles.py` | 1.2 mm contiguous tiles at native density, 3 per section — needed because commot 0.0.3 builds a dense n×n distance matrix |
| `phase4_methods.py` | COMMOT is called as installed; CellChat v2, SpaTalk and NCEM-linear are **reimplementations of the published statistic**, labelled as such everywhere |
| `phase4_run.py` | `--method commot\|cellchat\|spatalk\|ncem`, one checkpoint CSV per (method, tile, condition[, LR pair]) |
| `phase4_diag.py` | how destructive each null actually is (ligand⁺ cells retaining a neighbour after the shift) |
| `phase4_commot_mechanism.py` | cell-level vs cluster-level effect of coordinate permutation on COMMOT — the mechanism result |
| `phase4_decomp.py` | analytic spatial content of an edge-averaged LR score |
| `phase4_positive_control.py` | synthetic controls C1–C4: what each statistic can and cannot see, on data with a known answer |
| `phase4_summarize.py`, `phase4_tables.py`, `phase4_report.py` | aggregation, markdown tables, and rendering of `code/CS_PHASE4_template.md` into the report so no number is typed by hand |
| `make_figure4.py`, `make_figure4_supp.py`, `make_figure4_supp2.py` | Figure 4 and its two supplements |

### Gotchas (Phase 4)
* `nproc` says 48; `/sys/fs/cgroup/cpu.max` says **10.2 cores**. Size worker
  pools from the cgroup, not from `nproc`, or every worker runs at 20 % speed.
* `commot` 0.0.3 needs `np.Inf` → `np.inf` under NumPy 2, and materialises a
  dense n×n distance matrix — 105 GB for a 115 k-cell section.
* `ncem` is `python_requires <=3.10` and imports TensorFlow at package level;
  installing TF would pin `numpy<2` and break every other phase.
* `np.add.at` is unbuffered and dominated the NCEM design build; `np.bincount`
  on a flattened index is ~50× faster.
* `cKDTree.query(..., workers=-1)` is *slower* than `workers=1` for small
  queries inside an already-parallel loop (thread fan-out cost).
* `pkill -f <pattern>` matched the launching shell again. Use
  `ps -eo pid,args | grep ... | awk '{print $1}' | xargs kill`.
* Overwriting a log with `>` while orphaned workers still hold the old fd gives
  you a file padded with NULs at their old offsets.
