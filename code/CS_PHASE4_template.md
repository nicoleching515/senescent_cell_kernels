# CS Phase 4 — Existing spatial CCC tools under the same coordinate nulls (Figure 4)

Master Plan §23 (baselines), §25 Figure 4, §29 (positioning).
Deliverables: `figures/figure4.{png,pdf}` + `figures/figure4_data.csv`;
supplements `figures/figure4_supp_ncem_lengthscale.{png,pdf}` and
`figures/figure4_supp_commot_mechanism.{png,pdf}`;
`results/phase4/interactions.csv.gz` and
`results/phase4/{headline,headline_offdiagonal,headline_by_split,
verdict,score_rank_correlation,tiles,null_destructiveness,commot_mechanism,
cellchat_summary_statistic,ncem_radius_sweep,spatial_content_of_edge_scores,
positive_controls}.csv`;
`results/phase4/parts/` (per-job checkpoints, ~0.5 GB, git-ignored and
regenerable); code in `code/phase4_*.py` and `code/make_figure4*.py`; the report
template is `code/CS_PHASE4_template.md` (§6).

---

## 0. Headline

<!--HEADLINE-->

---

## 1. What ran as published software, and what did not

This is the first thing a reviewer will check, so it is the first thing stated.

| Method | Language | Status here | Why |
|---|---|---|---|
| **COMMOT 0.0.3** | Python | **PUBLISHED SOFTWARE** — `pip install commot`, run through `ct.tl.spatial_communication` + `ct.tl.cluster_communication` | installed cleanly; one source edit was required and is recorded in §1.1 |
| **CellChat v2** | R | **REIMPLEMENTATION of the published statistic** | no R interpreter in this image; see §1.3 for exactly how far the R route was scoped. The CellChat R package was never executed |
| **SpaTalk** | R | **REIMPLEMENTATION of the published statistic** | same |
| **NCEM (linear variant)** | Python | **REIMPLEMENTATION of the published linear model** | `ncem` declares `python_requires <=3.10` (we are on 3.11) and imports TensorFlow at package level; installing TF would have pinned `numpy<2` and broken the environment every other phase was produced under |

Everything labelled with an asterisk in Figure 4, in `results/phase4/*.csv`
(`method` column) and in this report is a reimplementation. **No result in this
phase should be attributed to the CellChat, SpaTalk or NCEM software.** What is
being tested is the *statistic* each of those tools computes, written out from
the published description, run under the same nulls as COMMOT. That is
scientifically informative — the failure modes below are properties of the
estimator, not of anybody's code — but it is not a software benchmark and is not
presented as one.

### 1.1 The one change made to COMMOT, recorded in full

`commot/_optimal_transport/_usot.py` used `np.Inf`, removed in NumPy 2.0, so the
package did not import. One substitution, `np.Inf` → `np.inf`, in one default
argument. Nothing else in the package was touched; the optimal-transport
solver, the cluster summarisation and the permutation test are the released
code.

### 1.2 A scalability result that is worth reporting on its own

`commot.tools._spatial_communication` line 390 materialises a **dense n×n
distance matrix** (`scipy.spatial.distance_matrix(spatial, spatial)`) whenever
`adata.obsp['spatial_distance']` is absent. For our smallest admissible section
(114,721 analysable cells) that is **105 GB**, and the process is OOM-killed in
about ten seconds on a 251 GB machine. Every one of the six sections is larger
than that, up to 207,785 cells (346 GB).

**COMMOT 0.0.3 as released cannot be run on a whole Xenium Prime section.**
That is not a criticism of the method — it predates 100k-cell-per-section
imaging — but it is a fact a benchmark should state, and it is the reason for
the tiling in §2.2.

### 1.3 The R route, scoped rather than hand-waved

There is no R interpreter in the image (`which R` and `which Rscript` are both
empty, `/usr/lib/R` does not exist). It is not unavailable in principle:
`apt-cache policy r-base-core` offers **4.1.2-1ubuntu2**, and both packages
declare a low floor — CellChat 2.2.0 `Depends: R (>= 3.6.0)`, SpaTalk 1.0
`Depends: R (>= 4.0.0)`. What is expensive is the dependency stack, which for
CellChat v2 includes Seurat 5.x (itself `R (>= 4.0.0)`, SeuratObject, Matrix),
NMF, `ComplexHeatmap` from a Bioconductor release matched to R 4.1, `presto` and
`circlize`, all built from source. Against the container's measured **10.2-core
CPU quota** (§2.7) — already saturated for the whole of this phase by the runs
in §3 — that build was judged out of budget and was not started.

That is a resource decision, not a claim that it cannot be done, and §7 lists it
as the highest-value follow-up: the route is open, and validating the two
reimplementations against the real packages on two or three tiles would convert
two asterisks into a software benchmark.

### 1.4 What we did not attempt at all

SCILD, RGAST, HARMONIC and NICHES are named in the Master Plan's related-work
section but not in §23's baseline list, and are not run here. `ncem` was
attempted four ways (PyPI, git, `--no-deps`, `--ignore-requires-python`); the
last installs but cannot import without TensorFlow, and TensorFlow pins
`numpy<2` against the `numpy 2.4.6` every other phase of this project was
produced under. That attempt is logged in `/tmp/inst/ncem*.log` and took under
25 minutes, inside the 45-minute installation budget.


---

## 2. Design

### 2.1 What is held fixed from earlier phases

* **Sections.** The six Section-8 Test-3 admissible sections, six animals, both
  arms: `7259`, `7260` (SBR, 26 wk), `7001` (sham, 52 wk), `7248` (sham, 26 wk),
  `7352` (sham, 2 wk), `7435` (sham, 10 wk). `7239`, `7448`, `7361`, `7450` are
  over the 20 % sender-prevalence ceiling and `7250` under the 1 % floor.
* **Cells.** The Phase 3 analysis cell set (`data/processed/cache3`), which is
  QC-passing cells carrying a Bio cell-type label, with `Low_quality` and
  `Unknown` dropped. Labels are `cell_type_merged` (Bio Phase 3 §1.1), the only
  label family comparable across sections.
* **Ligand–receptor pairs.** Exactly the four BIO_PHASE3 §3.4 certifies as
  defensible on this panel: `Ccl2→Ccr2`, `Tnf→Tnfrsf1a/Tnfrsf1b`,
  `Tgfb1→Tgfbr1/Tgfbr2`, `Il1a→Il1r1`. `Ackr3` is detected in ≤1 % of cells so
  CCL2/CXCL12→ACKR3 is not testable; `Il6` (≈3 senescent cells per section),
  `Cxcl2` (≈8) and `Mmp3` (≈3) carry no claims; mouse has no CXCL8 orthologue.
  Where a pair has two alternative receptors they are combined per cell as an
  **OR** (elementwise max of the log-normalised values), not as a complex.
* **Window.** 100 µm, the CS_PHASE3 §2 value, used as COMMOT's `dis_thr` and as
  the top of the NCEM radius sweep.
* **Geometry.** `scipy.spatial.cKDTree` throughout. No (n,n) matrix is built by
  any code in `code/phase4_*.py`. COMMOT builds one internally; see §1.2.
* **Seeds.** `MASTER_SEED = 20260820`, one derived stream per job; every job
  checkpoints to `results/phase4/parts/<method>_<tile>_<cond>[_p<pair>].csv`.

### 2.2 Tiles, and why

Because COMMOT cannot take a whole section (§1.2), every method is run on the
same **spatially contiguous 1.2 mm × 1.2 mm tiles at native cell density** —
three per section, chosen as the three densest non-overlapping tiles, 18 tiles
in total, 5,352–11,613 cells each (`results/phase4/tiles.csv`). Contiguous
tiles rather than random subsamples: a random subsample of 10 k cells from a
120 k-cell section thins local density by ~12×, which changes what any
distance-thresholded method sees. A tile does not.

Running *all four methods on identical cells* also means the differences in
Figure 4 are differences between methods, not between data scopes.

Cell types with fewer than 30 cells in a tile are dropped from that tile, since
no group-level statistic exists for them; that leaves 7–10 receiver/sender types
per tile.

### 2.3 The unit of analysis

One **interaction** = (tile, LR pair, sender cell type, receiver cell type).
Every method emits a score matrix and a p-value matrix over sender × receiver
types for one LR pair on one set of coordinates, which is what makes
"significant on real, still significant on shifted" comparable across methods.
Significance is `p < 0.05` for COMMOT, CellChat* and SpaTalk*, and BH-FDR
`q < 0.05` for NCEM* (NCEM reports FDR, so it gets FDR).

Sender = receiver (autocrine / homotypic) entries are legitimate outputs of all
four methods and are kept; `results/phase4/headline_offdiagonal.csv` repeats
every number with the diagonal removed. Every survival fraction moves by less
than 0.03 except NCEM* under N3t, which falls from 0.080 to 0.013 — i.e. the one
method that passes a null passes it *harder* off the diagonal. No conclusion
depends on the diagonal.

### 2.4 The coordinate nulls

Four, forming a ladder of severity. `torus_shift` and `rotate_about_centroid`
are **copied verbatim from `code/run_phase3_nulls.py`** so that Figure 4 sits on
exactly the same footing as Figure 2c.

| null | operation | what it preserves | what it destroys |
|---|---|---|---|
| **N3_lig** | torus-shift the **ligand⁺ cells** of that LR pair by one uniform random vector, wrapped on the tile | clustering of the ligand field; everything about the receiver field | alignment between the ligand field and everything else |
| **N4_lig** | rotate the ligand⁺ cells about the tile centroid, wrapped | same | same |
| **N3_type** | shift **each cell type by its own independent random vector** | each cell type's internal spatial clustering, exactly | all between-cell-type alignment, for every directed pair at once |
| **N0_perm** | permute all cell coordinates among cells | the point pattern and the marginal expression distribution | all association between position and cell identity — this is CellWHISPER's randomisation |

N3_lig/N4_lig are *our project's own N3/N4*, transposed to a ligand-receptor
setting: the "sender set" is the set of cells expressing that ligand. N3_type
exists because a method whose statistic is defined at the cell-type level (all
four of these are) is structurally almost blind to a ligand-field displacement —
that turns out to be the single most important finding in this phase, and it
would have been invisible with N3_lig alone. N0_perm is the destroy-everything
reference: any interaction still called significant under N0_perm is a Type I
error attributable to the method's own null model, not to the shuffle being
weak.

**The tiling makes the shift clean.** A torus shift on a whole section throws
~20 % of shifted cells into the void outside the tissue, which would weaken the
shifted condition for the wrong reason. On the tiles, which are solid tissue,
**100 % of shifted ligand⁺ cells retain a neighbour within 100 µm and the median
neighbour count is unchanged** (real 150.8, shifted 149.4;
`results/phase4/null_destructiveness.csv`). The shift changes alignment and
essentially nothing else.

`Tnf` is detected in 0.09–0.22 % of cells in the six sham tiles of `7248` and
`7435`, i.e. 6–15 ligand⁺ cells per tile. Below a floor of 20 ligand⁺ cells the
ligand-field nulls are not run and those cells are excluded from the N3_lig/
N4_lig rows (6 of 18 tiles for `Tnf` only). N3_type and N0_perm are unaffected.

### 2.5 The two numbers reported for every (method, null)

* **Significance survival** — of the interactions a method calls significant on
  **real** coordinates, what fraction it still calls significant on shuffled
  coordinates, averaged over replicates. This is the CellWHISPER quantity.
* **Score surviving fraction** — median over replicates of
  (null score ÷ real score). This is *our project's* metric (CS_PHASE3 §0), so
  the two lines of evidence are directly comparable.

The pair is what separates the two failure modes:

> **score SF ≈ 1 and significance survival ≈ 1** → the null is too weak *for
> that statistic*: the shuffled tissue really does still produce the number the
> method is measuring, so continuing to reject is correct behaviour on an
> uninformative null.
>
> **score SF ≈ 0 but significance survival ≈ 1** → the null is fine and the
> method's own significance test is miscalibrated: it keeps calling interactions
> after the quantity it is testing has collapsed.

### 2.6 The statistics, as implemented

**COMMOT (published).** `ct.tl.spatial_communication(dis_thr=100, pathway_sum=True,
heteromeric=False)` with a `df_ligrec` containing only that pair's rows, followed
by `ct.tl.cluster_communication(clustering='celltype', n_permutations=100)`
(the package default). The score is COMMOT's cluster–cluster communication
matrix — the mean of the cell-level optimal-transport communication matrix over
the sender×receiver block — and the p-value is COMMOT's own permutation of cell
labels holding the transport plan fixed. Defaults everywhere else.

**CellChat v2\* (reimplementation).** For cell groups i, j:
`P_ij = Hill(L_i · R_j) × s(d_ij)`, with `L_i`, `R_j` the **triMean**
(`(Q25 + 2·Q50 + Q75)/4`) of ligand / receptor log-expression in the group,
`Hill(x) = x/(0.5 + x)`, and the v2 spatial constraint
`s(d) = min(1, d₀/d)` truncated to zero beyond a 250 µm interaction range, where
`d_ij` is the median distance from a cell of type i to the nearest cell of type
j and `d₀` is the tile's median nearest-neighbour distance. p-values by
permutation of cell group labels, 100 permutations (CellChat's `nboot`
default), recomputing both the expression terms and the geometry each time.

> **CellChat's default group summary is undefined on this panel, and that is a
> result in itself.** `computeCommunProb` defaults to `type = "triMean"`, the
> Tukey trimean `(Q25 + 2·Q50 + Q75)/4`. Our four ligands are detected in
> 0.4–11 % of cells, so `Q25 = Q50 = Q75 = 0` and **the trimean of the ligand is
> exactly zero in every cell type in 18 of 18 tiles for `Ccl2`, `Il1a` and
> `Tnf`, and in 3 of 18 for `Tgfb1`** (`results/phase4/cellchat_summary_statistic.csv`).
> With `L_i = 0` the communication probability is identically zero and CellChat
> at its defaults reports **no communication at all** for three of the four
> senescence-relevant pairs on this data. The 10 %-truncated mean, CellChat's
> other documented option, is zero for `Tnf` too. The null comparison therefore
> uses `type = "mean"`, the remaining documented option, which is what a user
> analysing sparse imaging data would have to fall back to. Any CellChat number
> in Figure 4 is a `type = "mean"` number.

**SpaTalk\* (reimplementation).** A k-nearest-neighbour spatial graph (k = 10,
SpaTalk's default), edges kept in both orientations. For sender type A and
receiver type B, `LRscore = √L̄·√R̄ / (√L̄·√R̄ + 1)` with `L̄` the mean ligand
expression over the A-endpoints of A→B edges and `R̄` the mean receptor
expression over the B-endpoints. p-values by permutation of cell labels;
**200 permutations rather than SpaTalk's default 1,000**, for runtime — p is
resolved to 0.005, ample at α = 0.05.

**NCEM linear\* (reimplementation).** For each receiver type t, ordinary least
squares of the pair's receptor score on the counts of neighbours of each type s
within radius r: `y_i = μ_t + Σ_s β_{ts}·x_{is} + ε`. The coupling s→t is
called at BH-FDR `q < 0.05` on the Wald test for `β_{ts}`. The interaction
radius is chosen the way NCEM chooses it — the r maximising variance explained,
swept over {10, 15, 20, 30, 40, 50, 75, 100} µm — and the sweep is repeated
under every null, because a method whose *length scale* is unchanged by
coordinate shuffling has not measured a length scale. For a linear model this
is exact up to the optimiser (closed-form OLS instead of Keras SGD).

### 2.7 Replicates, and the honest accounting

| method | replicates per tile per null | per section | permutations inside the method |
|---|---|---|---|
| COMMOT | 10 | **30** | 100 (package default) |
| CellChat v2* | 25 | **75** | 100 (CellChat default) |
| SpaTalk* | 100 | **300** | 200 (reduced from 1,000) |
| NCEM linear* | 100 | **300** | — (analytic Wald + BH) |

The plan asks for ≥100 shift/rotation replicates per section. SpaTalk* and
NCEM* clear it by 3×. **CellChat\* gets 75 and COMMOT gets 30.** The reason is
arithmetic, and it is stated rather than hidden: one COMMOT run — one tile, one
LR pair, one coordinate set — costs ~8.7 s, the design needs 288 such jobs ×
10 replicates, and **the container's CPU quota is 10.2 cores, not the 48 that
`nproc` reports** (`/sys/fs/cgroup/cpu.max = 1020000 100000`, measured, not
assumed). At 30 replicates per section COMMOT alone is ~7 CPU-hours.

Nothing in the conclusion is replicate-limited. The headline quantity is an
average of a per-interaction survival fraction over ~10³ interactions; even at
COMMOT's 10 replicates per interaction the binomial contribution to the
average's standard error is under 0.01, an order of magnitude smaller than the
0.7–0.8 that separates the methods that fail from the one that does not.
Replicate count would matter if a survival fraction sat near a decision
boundary. None does.



<!--RESULTS-->

---

## 5. Limitations, stated before the conclusions are used

1. **Three of the four methods are reimplementations.** They test the
   *statistic*, not the software. A discrepancy between these results and the
   published packages would most likely mean the reimplementation is wrong, and
   the paper must say so. Only COMMOT is a software benchmark.
2. **COMMOT ran on 1.2 mm tiles, not whole sections**, because the released code
   cannot take a whole section (§1.2). Cluster-level statistics on 5–12 k cells
   are noisier than on 120 k, and a tile sees less anatomical heterogeneity than
   a section. All methods were run on the identical tiles so the *comparison* is
   unaffected, but the absolute significance rates are tile-scale numbers.
3. **COMMOT gets 30 shuffle replicates per section**, below the plan's 100
   (§2.7). CellChat* gets 75.
4. **The nulls move cells, not expression.** A cell carries its transcriptome
   and its type label wherever it is moved. Nulls that instead permuted
   expression would test a different hypothesis, and our own N1 (Phase 3) is
   that null.
5. **Four ligand–receptor pairs, one tissue, one species, one disease model.**
   `Tnf` is at 0.09–0.22 % detection in two sham sections, so its ligand-field
   nulls rest on 6 of 18 tiles being excluded and the rest carrying 20–131
   ligand⁺ cells. The CellWHISPER comparison is a comparison on
   *senescence-relevant pairs in mouse liver*, not a general replication of their
   benchmark.
6. **`p < 0.05` uncorrected for COMMOT, CellChat* and SpaTalk***, because that
   is what each method reports and what a user would act on. NCEM* gets BH-FDR
   because NCEM reports FDR. Applying BH to all four would lower every
   significance rate but not the survival fractions, which are ratios.
7. **Significance survival is not a false-positive rate.** It is the fraction of
   real-data calls reproduced on shuffled data. It equals an FPR only if one
   accepts the shuffled tissue as a true null, which §2.4 and CS_PHASE3 §5 both
   argue is exactly the assumption in question.

---

## 6. Reproduce

```bash
cd /workspace/code
python3 phase4_data.py                                  # LR-gene cache from the h5
python3 phase4_tiles.py                                 # results/phase4/tiles.csv
python3 phase4_diag.py                                  # null destructiveness
python3 phase4_run.py --method ncem     --reps 100 --jobs 2
python3 phase4_run.py --method spatalk  --reps 100 --jobs 2
python3 phase4_run.py --method cellchat --reps 25  --jobs 2
python3 phase4_run.py --method commot   --reps 10  --jobs 5
python3 phase4_summarize.py                             # interactions/headline/verdict
python3 phase4_commot_mechanism.py                      # commot_mechanism.csv
python3 phase4_decomp.py                                # spatial content of edge scores
python3 phase4_positive_control.py                      # synthetic controls C1-C4
python3 make_figure4.py                                 # figures/figure4.{png,pdf}
python3 make_figure4_supp.py                            # NCEM length-scale supplement
python3 make_figure4_supp2.py                           # COMMOT mechanism supplement
python3 phase4_report.py                                # fills reports/CS_PHASE4.md
```

`phase4_report.py` renders `code/CS_PHASE4_template.md` into
`reports/CS_PHASE4.md`, substituting every number from `results/phase4/`, so the
prose cannot drift from the data. Edit the template, not the report.

Every job is skipped if its checkpoint exists, so the run is restartable.
`MASTER_SEED = 20260820`; each job derives its own `default_rng` stream from it
and from the job index, so results are reproducible job by job and independent
of scheduling order.

**Environment note.** `commot==0.0.3` was installed from PyPI on 2026-08-20 and
patched as in §1.1; `ncem` is installed but unusable (§1). The COMMOT patch must
be reapplied on a fresh environment, or `numpy<2` pinned instead.

---

## 7. What I did not get to

* **CellChat v2 and SpaTalk as software.** `r-base-core` 4.1.2 is apt-installable
  here and both packages' declared R floors are met (§1.3); the cost is the
  Seurat/NMF/Bioconductor build. Running the real packages on two or three tiles
  and checking that the reimplementations reproduce their calls would convert
  two asterisks into a validation, and is the single highest-value remaining item
  in this phase.
* **NCEM as software**, which needs a Python ≤3.10 environment with TensorFlow —
  cheap in a separate venv, and worth doing for the same reason.
* **SCILD and HARMONIC**, the two most recent methods, are not in §23's baseline
  list and were not attempted.
* **Whole-section COMMOT**, which would need `cot_sparse` rewritten against a
  cKDTree sparse distance matrix. That is a fork of the method, not a run of it,
  and does not belong in a benchmark of the released tool.
* **A confounder-aware null for these tools.** CellWHISPER's own null brought
  their false-positive rate under 5 %. We have not run it. Our N1 (cell-type-
  stratified label permutation) is the nearest thing this project has, and
  §2.4's N3_type is a partial substitute; a proper comparison against
  CellWHISPER's null is future work.
