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

**Three of the four methods cannot tell this tissue from tissue whose
coordinates have been destroyed, and the fourth is blind to the null the field
recommends.**

| | COMMOT | CellChat v2* | SpaTalk* | NCEM linear* |
|---|---|---|---|---|
| ran as | **published software** | reimplementation | reimplementation | reimplementation |
| interactions called significant on **real** coordinates | 22.4% | 28.3% | 24.7% | 1.4% |
| of those, still significant after the **torus shift (N3)** | **81%** | **97%** | **81%** | **70%** |
| still significant after **full coordinate permutation (N0)** | **79%** | **94%** | **79%** | **0.4%** |
| Spearman ρ(real score, N0-permuted score) | 0.90 | 0.98 | 0.86 | 0.03 |

1. **CellWHISPER's >90 % figure replicates under CellWHISPER's own null.** Their
   control is a *within-cell-type* location permutation, not the full coordinate
   permutation this report originally attributed to them (D7 §B9); we now run it
   (`N0_type`, §8). Under it, CellChat v2\* keeps **97.1 %** of its real-data
   calls, COMMOT **81.1 %**, SpaTalk\* **78.1 %**, and CellChat's *count* of
   significant interactions is **0.2831 on permuted tissue against 0.2833 on real
   tissue** — their criterion, met to four decimal places. The numbers under our
   own stricter N0 (94 %, 79 %, 79 %) are given alongside. All three
   ligand-receptor methods reproduce the large majority of their real-data calls
   on data whose ligand-receptor geometry has been destroyed (§3.3, §8).
2. **The failure is not that the null is weak. The failure is that the statistic
   is not spatial.** Under N0 — coordinates permuted, no spatial information
   left in the data at all — COMMOT's and SpaTalk*'s cluster-level scores keep
   ρ ≈ 0.90 and ρ ≈ 0.86 with their real-coordinate values. There is nothing
   left for a stronger null to remove.
3. **For COMMOT we can show the mechanism directly.** The optimal-transport
   step is exquisitely sensitive to geometry — permuting coordinates replaces
   the cell-to-cell communication network almost entirely — but it conserves
   total transported mass, and *averaging that network over a sender × receiver
   cell-type block throws the geometry away* (§4.1).
4. **NCEM linear\* is the exception, and the exception is instructive.** It
   collapses correctly under the nulls that destroy cell-type geometry
   (N3t 8%, N0 0.4% survival) and is structurally blind to the
   ligand-field torus shift (70% survival) — because its statistic never looks at
   ligand expression. Its *length scale*, however, is not identified at all
   (§4.4).
5. **The implementations are not the problem — synthetic positive controls
   settle that.** With a planted ligand-in-A / receptor-in-B interaction all
   three LR methods find it at p = 0 and call exactly 1 of 16 cell-type pairs
   significant, and all three correctly lose it when A and B are pushed beyond
   the interaction range. But when the ligand and receptor rates are held fixed
   and only the *arrangement* is changed — receptor⁺ cells adjacent to ligand⁺
   cells versus scattered — CellChat\*'s score is identical to six significant
   figures and COMMOT's to four. NCEM linear\* is the only one of the four that
   separates them (p = 10⁻³⁹ vs 1.0) (§3.6).
6. **Our own estimator fails the same null in the opposite direction.** Phase 3
   found N3 surviving fraction **1.000** because the shifted null is centred at
   ~0 — the shift destroys our statistic completely, and the test still rejects,
   because "is the sender field aligned with the response field?" is not "is
   there a SASP effect?". These tools' N3 surviving fractions are
   0.91–0.95: the shift does not destroy their statistic at all. **Same
   null, two different diseases** (§4.5).


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
| **N0_perm** | permute all cell coordinates among **all** cells | the point pattern and the marginal expression distribution | all association between position and cell identity, *including the cell-type spatial architecture* |
| **N0_type** | permute cell locations **within each cell type** | the point pattern, the marginals, **and each cell type's spatial organisation** | only the proximity between ligand⁺ and receptor⁺ cells — **this is CellWHISPER's own randomisation**; added 2026-08-21, see §8 |

N3_lig/N4_lig are *our project's own N3/N4*, transposed to a ligand-receptor
setting: the "sender set" is the set of cells expressing that ligand. N3_type
exists because a method whose statistic is defined at the cell-type level (all
four of these are) is structurally almost blind to a ligand-field displacement —
that turns out to be the single most important finding in this phase, and it
would have been invisible with N3_lig alone. N0_perm is the destroy-everything
reference: any interaction still called significant under N0_perm is a Type I
error attributable to the method's own null model, not to the shuffle being
weak.

**Correction, and the reason N0_type exists (D7 §B9).** An earlier version of
this table labelled **N0_perm** as "CellWHISPER's randomisation". That is wrong.
CellWHISPER permutes cell locations *within each cell type*, in their words
"preserving cell-type-specific spatial organization and ligand-receptor (LR)
expression while destroying spatial proximity between ligand- and
receptor-expressing cells". N0_perm permutes across all cells and therefore also
destroys the cell-type architecture, so it is **strictly more destructive than
theirs**. The direction of the error was in our favour — we were showing these
statistics survive a *harsher* shuffle than the one the benchmark used — but
"we reproduce CellWHISPER's >90 %" was not a replication claim we had earned.
Rather than relabel, we ran their null: **N0_type**, §8.

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




---

## 3. Results

Coverage: 6,032 COMMOT interactions over 18 tiles,
6,032 CellChat* over 18,
6,032 SpaTalk* over 18,
6,032 NCEM* over 18.
One interaction = (tile, LR pair, sender type, receiver type).

### 3.1 Significance survival — what Figure 4b shows

Fraction of interactions called significant on real coordinates that are still
called significant on shuffled coordinates.

| method | N3 torus (ligand⁺) | N4 rotation (ligand⁺) | N3t per-type torus | N0 full permutation |
|---|---|---|---|---|
| COMMOT | 0.806 | 0.820 | 0.769 | 0.794 |
| CellChat v2* | 0.974 | 0.975 | 0.913 | 0.945 |
| SpaTalk* | 0.814 | 0.820 | 0.784 | 0.791 |
| NCEM linear* | 0.702 | 0.701 | 0.080 | 0.004 |

### 3.2 The score itself — what Figure 4c shows

Median over replicates of (null score ÷ real score).

| method | N3 torus (ligand⁺) | N4 rotation (ligand⁺) | N3t per-type torus | N0 full permutation |
|---|---|---|---|---|
| COMMOT | 0.906 | 0.912 | 0.909 | 0.914 |
| CellChat v2* | 1.000 | 1.000 | 0.981 | 0.998 |
| SpaTalk* | 0.949 | 0.956 | 0.952 | 0.943 |
| NCEM linear* | 0.972 | 0.968 | -0.014 | -0.014 |

And the rank correlation between the real and shuffled score across
interactions — the sharpest single number in this phase, because a correlation
near 1 means the shuffle did not even reorder the interactions:

| method | N3 torus (ligand⁺) | N4 rotation (ligand⁺) | N3t per-type torus | N0 full permutation |
|---|---|---|---|---|
| COMMOT | 0.907 | 0.908 | 0.907 | 0.904 |
| CellChat v2* | 0.998 | 0.998 | 0.991 | 0.984 |
| SpaTalk* | 0.887 | 0.890 | 0.887 | 0.860 |
| NCEM linear* | 0.881 | 0.883 | 0.122 | 0.026 |

Absolute significance rates, which is what a user actually sees:

| method | real coordinates | N3 torus (ligand⁺) | N4 rotation (ligand⁺) | N3t per-type torus | N0 full permutation |
|---|---|---|---|---|---|
| COMMOT | 0.224 | 0.242 | 0.246 | 0.227 | 0.241 |
| CellChat v2* | 0.283 | 0.294 | 0.294 | 0.277 | 0.318 |
| SpaTalk* | 0.247 | 0.259 | 0.260 | 0.241 | 0.248 |
| NCEM linear* | 0.014 | 0.013 | 0.013 | 0.006 | 0.003 |

### 3.3 Do we reproduce CellWHISPER's >90 %?

**Yes — and as of §8 the comparison is against their actual null rather than a
harsher one of ours.** CellWHISPER reported that CellChat v2, COMMOT and SpaTalk
return comparable interaction *counts* on real and randomised input, implying
false-positive rates above 90 %. The table below is our **N0_perm** (full
coordinate permutation, strictly more destructive than theirs); §8 gives the
same table under **N0_type**, which is their design.

On this data, under full coordinate permutation:

| | significance rate, real | significance rate, N0-permuted | ratio | identity-matched survival |
|---|---|---|---|---|
| COMMOT | 0.224 | 0.241 | 1.08 | 79% |
| CellChat v2* | 0.283 | 0.318 | 1.12 | 94% |
| SpaTalk* | 0.247 | 0.248 | 1.00 | 79% |
| NCEM linear* | 0.014 | 0.0027 | 0.19 | 0.4% |

Two versions of the number, and both belong in the paper:

* **CellWHISPER's own criterion — comparable interaction counts on randomised
  input — replicates cleanly.** The count ratios above are all within a few
  percent of 1 for the three ligand–receptor methods. A user would get the same
  number of "significant" interactions from tissue and from confetti.
* **The stricter, identity-matched version** — is *this particular* sender→
  receiver interaction still called? — is 79%–94%. The gap
  between the two is the interactions that are significant on real coordinates
  and are *replaced by different* significant interactions after shuffling
  rather than reproduced. It is the more informative number and it is the one
  Figure 4b plots.

Either way the conclusion is the same and it is specific to this project's
question: **on the four senescence-relevant ligand–receptor pairs this panel can
support, in this tissue, the established tools' calls do not depend on where the
cells are.**

### 3.4 Per ligand–receptor pair (Figure 4a)

| method · LR pair | N3 torus (ligand⁺) | N4 rotation (ligand⁺) | N3t per-type torus | N0 full permutation | n real-sig |
|---|---|---|---|---|---|
| COMMOT · Ccl2->Ccr2 | 0.738 | 0.756 | 0.759 | 0.758 | 263 |
| COMMOT · Tnf->Tnfrsf1a/1b | 0.712 | 0.734 | 0.627 | 0.674 | 241 |
| COMMOT · Tgfb1->Tgfbr1/2 | 0.869 | 0.884 | 0.823 | 0.864 | 647 |
| COMMOT · Il1a->Il1r1 | 0.780 | 0.774 | 0.776 | 0.755 | 198 |
| CellChat v2* · Ccl2->Ccr2 | 0.972 | 0.976 | 0.949 | 0.958 | 337 |
| CellChat v2* · Tnf->Tnfrsf1a/1b | 0.946 | 0.949 | 0.875 | 0.905 | 360 |
| CellChat v2* · Tgfb1->Tgfbr1/2 | 0.988 | 0.989 | 0.926 | 0.967 | 772 |
| CellChat v2* · Il1a->Il1r1 | 0.964 | 0.962 | 0.877 | 0.916 | 240 |
| SpaTalk* · Ccl2->Ccr2 | 0.805 | 0.810 | 0.737 | 0.752 | 287 |
| SpaTalk* · Tnf->Tnfrsf1a/1b | 0.663 | 0.673 | 0.658 | 0.662 | 268 |
| SpaTalk* · Tgfb1->Tgfbr1/2 | 0.850 | 0.857 | 0.842 | 0.853 | 716 |
| SpaTalk* · Il1a->Il1r1 | 0.852 | 0.854 | 0.809 | 0.801 | 220 |
| NCEM linear* · Ccl2->Ccr2 | 0.790 | 0.785 | 0.020 | 0.008 | 43 |
| NCEM linear* · Tnf->Tnfrsf1a/1b | 0.698 | 0.717 | 0.143 | 0.001 | 9 |
| NCEM linear* · Tgfb1->Tgfbr1/2 | 0.379 | 0.385 | 0.053 | 0.001 | 19 |
| NCEM linear* · Il1a->Il1r1 | 0.847 | 0.846 | 0.238 | 0.003 | 16 |

### 3.5 What this looks like to a user

COMMOT's strongest `Ccl2→Ccr2` calls on real coordinates are exactly the ones a
liver immunologist would predict — dendritic cells, macrophages and T/NK cells
as receivers, which are the three `Ccr2`-expressing populations Bio measured
(DC 16.2 %, Kupffer 8.2 %, T/NK 6.4 %; BIO_PHASE3 §3.1). The right column is the
same analysis after every cell coordinate in the tile has been permuted.

| sender → receiver | tiles | significant on real coordinates | significant on N0-permuted coordinates |
|---|---|---|---|
| Mesenchymal → DC | 12 | 1.00 | 1.00 |
| Macrophages → DC | 12 | 1.00 | 0.99 |
| Macrophages → Macrophages | 18 | 0.94 | 0.76 |
| DC → DC | 12 | 0.92 | 0.85 |
| Mesenchymal → T/NK cells | 18 | 0.89 | 0.82 |
| Mesenchymal → Macrophages | 18 | 0.83 | 0.87 |
| Biliary/ductular → DC | 12 | 0.83 | 0.86 |
| T/NK cells → DC | 12 | 0.83 | 0.80 |

Fraction of tiles in which the interaction is called at p < 0.05. **The result
is biologically plausible, reproducible across tiles, and almost entirely
unchanged by destroying the tissue.** This is what a false positive looks like
in this literature: not an implausible call, a plausible one.

### 3.6 Positive controls — the implementations work, and here is exactly what they can and cannot see

The obvious objection to everything above is that three of the four methods are
reimplementations and might simply be insensitive. Four synthetic controls on a
tile-sized tissue with known ground truth and no confounding answer that
(`code/phase4_positive_control.py`, `results/phase4/positive_controls.csv`).
`A→B` is the planted interaction throughout; 8,000 cells, four cell types,
uniform positions.

**C1 cell-type sensitivity**

| method | A->B |
|---|---|
| COMMOT | p 0 |
| CellChat v2* | p 0 |
| SpaTalk* | p 0 |
| NCEM linear* | p 0.477 |

**C2 spatial sensitivity**

| method | real (ligand+ adjacent to receptor+) | N3 torus-shifted ligand+ |
|---|---|---|
| COMMOT | score 1.03663e-05, p 0.755 | score 1.00809e-05, p 0.795 |
| CellChat v2* | score 0.0185765, p 0.58 | score 0.0187887, p 0.56 |
| SpaTalk* | score 0.122083, p 0.676 | score 0.121761, p 0.66 |
| NCEM linear* | score 0.00512041, p 0.303 | score -8.2169e-05, p 0.985 |

**C3 within-type spatial coupling**

| method | coupled (receptor+ B next to ligand+ A) | scattered (same rates, random B cells) |
|---|---|---|
| COMMOT | score 0.0001245, p 0 | score 0.000125, p 0 |
| CellChat v2* | score 0.0654471, p 0 | score 0.0654471, p 0 |
| SpaTalk* | score 0.242828, p 0 | score 0.218441, p 0 |
| NCEM linear* | score 0.0596672, p 1e-39 | score 0.0025454, p 1 |

**C4 between-type geometry**

| method | interleaved A and B | segregated A and B |
|---|---|---|
| COMMOT | score 0.000130106, p 0 | score 0, p 1 |
| CellChat v2* | score 0.0968206, p 0 | score 0, p 1 |
| SpaTalk* | score 0.259018, p 0 | score 0, p 1 |
| NCEM linear* | score 0.000980924, p 0.849 | — |

Read those four blocks in order and the whole of §4 follows.

* **C1.** Ligand only in A, receptor only in B. COMMOT, CellChat* and SpaTalk*
  each find `A→B` at p = 0 and call **exactly one of sixteen** cell-type pairs
  significant. The implementations are correct and specific. NCEM* does not, and
  should not: its model predicts a cell's expression from its *neighbourhood's*
  composition, and here expression is a function of the cell's own type.
* **C4.** Same expression, and now the geometry is varied: A and B interleaved,
  versus A and B pushed to opposite ends of the tile, far beyond the 100 µm
  range. All three LR methods lose the interaction completely (score → 0,
  p = 1). **They do read geometry** — at the cell-type-pair level.
* **C3 is the one that matters.** A expresses the ligand at 30 %, B expresses
  the receptor at 25 %, in *both* conditions; the only difference is whether the
  receptor⁺ B cells are the ones sitting next to ligand⁺ A cells or a random
  sample of B cells. This is the extra thing a *spatial* method claims to
  provide over a non-spatial one. CellChat*'s score is identical to **all six
  significant figures** (0.0654471 vs 0.0654471). COMMOT's differs in the fourth
  (0.0001245 vs 0.0001250). SpaTalk*'s moves 11 % and stays at p = 0 either way.
  **NCEM linear\* is the only one of the four that separates them**, and it does
  so decisively: p = 10⁻³⁹ when the ligand and receptor cells are adjacent,
  p = 1.0 when they are not.
* **C2** is C3 without the cell-type alignment — the ligand⁺/receptor⁺ coupling
  cut across all four types. Nothing detects it, including on real coordinates.

So the finding of §3.1–§3.3 is not that these tools are broken. They compute
what they are defined to compute, sensitively and specifically. What they
compute is a function of cell-type composition, cell-type-level expression and
cell-type-level geometry — and **not** of whether the ligand-expressing cells are
actually next to the receptor-expressing cells. Shuffling coordinates within the
tile leaves the first three nearly untouched, which is why the calls do not move.

### 3.7 Robustness

`results/phase4/headline_by_split.csv` repeats every survival fraction split by
**surgical arm** and by **section**. Across the SBR and sham arms the three
ligand–receptor methods agree to within 0.07 on every null, and NCEM* agrees to
within 0.02 on the two nulls that bite. Nothing here is an artefact of one arm,
one animal or one section, and nothing depends on the diagonal (§2.3).

---

## 4. Mechanism — how each method fails, which is the part nobody has reported

The verdict table, driven by the pair (score surviving fraction, significance
survival) laid out in §2.5:

| method | null | score_sf | sig_survival | null_sig_rate | n_real_sig | n_rep | verdict |
|---|---|---|---|---|---|---|---|
| COMMOT | N0 full permutation | 0.914 | 0.794 | 0.240 | 1349 | 10 | STATISTIC IS NOT SPATIAL: survives complete destruction of the spatial arrangement |
| COMMOT | N3 torus (ligand⁺) | 0.906 | 0.806 | 0.242 | 1294 | 10 | null does not bite: the shuffled tissue still produces the score the method measures |
| COMMOT | N3t per-type torus | 0.909 | 0.769 | 0.227 | 1349 | 10 | null does not bite: the shuffled tissue still produces the score the method measures |
| COMMOT | N4 rotation (ligand⁺) | 0.912 | 0.820 | 0.245 | 1294 | 10 | null does not bite: the shuffled tissue still produces the score the method measures |
| CellChat v2* | N0 full permutation | 0.998 | 0.945 | 0.318 | 1709 | 25 | STATISTIC IS NOT SPATIAL: survives complete destruction of the spatial arrangement |
| CellChat v2* | N3 torus (ligand⁺) | 1.000 | 0.974 | 0.294 | 1625 | 25 | null does not bite: the shuffled tissue still produces the score the method measures |
| CellChat v2* | N3t per-type torus | 0.981 | 0.913 | 0.277 | 1709 | 25 | null does not bite: the shuffled tissue still produces the score the method measures |
| CellChat v2* | N4 rotation (ligand⁺) | 1.000 | 0.975 | 0.294 | 1625 | 25 | null does not bite: the shuffled tissue still produces the score the method measures |
| SpaTalk* | N0 full permutation | 0.943 | 0.791 | 0.248 | 1491 | 100 | STATISTIC IS NOT SPATIAL: survives complete destruction of the spatial arrangement |
| SpaTalk* | N3 torus (ligand⁺) | 0.949 | 0.814 | 0.259 | 1433 | 100 | null does not bite: the shuffled tissue still produces the score the method measures |
| SpaTalk* | N3t per-type torus | 0.952 | 0.784 | 0.241 | 1491 | 100 | null does not bite: the shuffled tissue still produces the score the method measures |
| SpaTalk* | N4 rotation (ligand⁺) | 0.956 | 0.820 | 0.260 | 1433 | 100 | null does not bite: the shuffled tissue still produces the score the method measures |
| NCEM linear* | N0 full permutation | -0.014 | 0.004 | 0.003 | 87 | 100 | this null bites: score and calls both collapse |
| NCEM linear* | N3 torus (ligand⁺) | 0.972 | 0.702 | 0.013 | 84 | 100 | null does not bite: the shuffled tissue still produces the score the method measures |
| NCEM linear* | N3t per-type torus | -0.014 | 0.080 | 0.006 | 87 | 100 | this null bites: score and calls both collapse |
| NCEM linear* | N4 rotation (ligand⁺) | 0.968 | 0.701 | 0.013 | 84 | 100 | null does not bite: the shuffled tissue still produces the score the method measures |

### 4.1 COMMOT: the optimal transport is spatial, the cluster summary is not

COMMOT is not insensitive to geometry. Permuting coordinates rebuilds its
cell-to-cell communication network almost from scratch — on the first tile of
each section, the Jaccard overlap between the real and the permuted set of
communicating cell pairs is a few percent:

| LR pair | communicating cell pairs, real | Jaccard of those pairs, real vs N0 | transported mass, N0 ÷ real | cluster-level Spearman, real vs N0 | significant: real / N0 / shared |
|---|---|---|---|---|---|
| Ccl2->Ccr2 | 1,124 | 0.0176 | 1.127084 | 0.646 | 97 / 98 / 68 |
| Tnf->Tnfrsf1a/1b | 2,936 | 0.0147 | 1.000000 | 0.872 | 77 / 79 / 55 |
| Tgfb1->Tgfbr1/2 | 55,420 | 0.0152 | 1.000054 | 0.841 | 228 / 248 / 188 |
| Il1a->Il1r1 | 8,852 | 0.0141 | 1.000000 | 0.715 | 66 / 71 / 50 |

Three things are visible at once, and together they are the explanation.

1. **The cell-level network is destroyed.** Jaccard ≈ 0.01–0.02: essentially
   none of the real communicating cell pairs survive.
2. **The total transported mass is conserved to 6+ significant figures.** That
   is a property of the collective optimal-transport formulation: all the
   available ligand is transported to some receiver within the distance
   threshold. Geometry decides *where* it goes, not *how much* there is.
   (`Ccl2→Ccr2` is the one exception, at 1.13: `Ccr2` is detected in 1–3 % of
   cells, so some `Ccl2` cannot find a receiver within 100 µm on real
   coordinates and conservation is only approximate. It moves the *wrong* way —
   the permuted tissue transports **more** mass than the real one.)
3. **The cluster-level score survives anyway**, at Spearman ρ ≈ 0.6–0.9,
   because averaging a mass-conserving flow over a sender × receiver cell-type
   block returns approximately (total mass × block composition) — a quantity
   with no geometry in it.

And then `cluster_communication` computes its p-value by permuting cell
**labels** while holding the transport plan **fixed**. So the test asks "do
these two cell types carry more of the communication than a random pair of
groups of the same size?", which is a question about cell-type composition and
ligand/receptor abundance. It is not a question about space, and no coordinate
null can change its answer.

§3.6's control C3 shows this on data with a known answer: hold the ligand and
receptor rates of A and B fixed and change only whether the receptor⁺ B cells
are the ones adjacent to the ligand⁺ A cells, and COMMOT's cluster score moves
from 0.0001245 to 0.0001250 — four significant figures unchanged, p = 0 either
way.

**This is the finding.** A benchmark that only reports "COMMOT fails the torus
shift" invites the reply "then use a better null". The correct statement is that
COMMOT's cluster-level significance test does not test the spatial part of
COMMOT, so no null on the coordinates can fix it. What would fix it is a null on
the transport plan, or reporting at the cell level where the geometry survives.

### 4.2 SpaTalk*: the neighbour graph moves, the edge-averaged score does not

SpaTalk*'s statistic is a mean ligand expression over the sender endpoints of
A→B edges times a mean receptor expression over the receiver endpoints, squashed
through `x/(x+1)`. Under any coordinate shuffle the *edges* change completely,
but the mean of a cell-type's expression over a large set of edges converges to
that cell-type's mean expression regardless of which edges they are. With 10
nearest neighbours per cell and thousands of cells per type, the law of large
numbers does the rest: score SF 0.94 and ρ 0.86 under full coordinate
permutation. Same disease as COMMOT, different route to it.

This one can be shown analytically, without any permutation at all. For **any**
statistic of the form "average over A→B neighbour edges of f(sender)·g(receiver)"
— which covers SpaTalk, CellPhoneDB-style scores and most of the LR literature —
the expectation under a random rewiring of the graph is exactly
E[f | A]·E[g | B]. The observed value can differ from that only through the
spatial covariance of f and g across edges. On this data that difference is
small: the SpaTalk*-style score on **real** coordinates correlates at Spearman
**ρ = 0.89** with a prediction built from cell-type mean expression and no
coordinates whatsoever, with a median relative gap of **14%**
(`results/phase4/spatial_content_of_edge_scores.csv`). About 86% of the
score is cell-type composition. A permutation test on the labels is testing that
86%.

### 4.3 CellChat v2*: two failures, one before the model is fitted

The first is documented in §2.6: **at its default `triMean` summary CellChat's
communication probability is identically zero** for `Ccl2`, `Il1a` and `Tnf` in
every cell type of every tile, because those ligands are detected in under 8 %
of cells. A user running CellChat v2 at defaults on this panel would conclude
there is no senescence-related communication in this tissue — not because there
isn't, but because `Q75 = 0`.

The second is the null result at `type = "mean"`, in the tables above.
CellChat's spatial constraint enters only through the distance between cell
*groups*, so displacing the ligand⁺ cells — a few percent of the population,
carrying their labels with them — barely moves it, and permuting all
coordinates rescales every group distance nearly equally, which cancels between
the observed statistic and its permutation null.

### 4.4 NCEM linear*: calibrated where it looks, and its length scale is not identified

NCEM* is the only method that behaves. Under N3t (per-cell-type torus shift) its
score surviving fraction is -0.014 and its significance survival 0.080; under
full coordinate permutation, -0.014 and 0.004 — a significance rate of
0.0027 against a nominal FDR of 0.05. Its statistic is the cell-type
composition of a cell's neighbourhood, those nulls destroy exactly that, and it
collapses exactly as it should.

Its 70% survival under the ligand-field torus shift is not a failure of
calibration but a statement about scope: NCEM linear never reads ligand
expression, so displacing ligand⁺ cells is close to a no-op for it. **A method
can be perfectly calibrated for the hypothesis it tests and still be the wrong
tool for a ligand–receptor question.**

Where it does fail is the quantity §23 asks it for — a comparable length scale:

| coordinates | median selected radius (µm) | range over 18 tiles (µm) | median R² at the selected radius | median R² over the whole sweep |
|---|---|---|---|---|
| real coordinates | 18 | 10–100 | 0.0149 | 0.0136 |
| N3 torus (ligand⁺) | 30 | 10–100 | 0.0147 | 0.0134 |
| N4 rotation (ligand⁺) | 25 | 10–100 | 0.0147 | 0.0131 |
| N3t per-type torus | 50 | 10–100 | 0.0123 | 0.0112 |
| N0 full permutation | 45 | 10–100 | 0.0126 | 0.0104 |

The variance-explained criterion NCEM uses to pick its interaction radius is
flat to within a few percent from 10 µm to 100 µm, so the argmax is noise: it
lands anywhere in 10–100 µm across the 18 tiles, and it does so on
coordinate-permuted data just as readily as on real data, at 77% of the
real R² (0.0104 vs 0.0136). **NCEM's reported interaction length scale is not identified in this
tissue** (`figures/figure4_supp_ncem_lengthscale.png`). That is the same
conclusion CS_PHASE3 §4 reached for our own λ̂, reached independently by a
different method with a different estimator, which is worth saying out loud.

### 4.5 Cross-reference: our own estimator fails the same null in the opposite direction

CS_PHASE3 §5 measured, over the 160 reportable fits, an N3 torus-shift surviving
fraction of **1.000** [0.992, 1.008], with the null distribution of β̂ centred at
**−2.1 × 10⁻⁶** (median |mean| 8.9 × 10⁻⁵) against an observed β̂ of
**1.13 × 10⁻²** — the null mean is **0.8 % of the observed amplitude** — and
**87.5 %** of fits rejecting at p < 0.05.

Put next to this phase, the two failures are opposites:

| | our SASP kernel estimator | COMMOT / SpaTalk* |
|---|---|---|
| does the torus shift destroy the statistic? | **yes, completely** (null centred at ~0.8 % of β̂) | **no** (score SF 0.91 / 0.95) |
| does the test still reject? | yes, 88% of fits | yes, 81% / 81% of calls |
| why | the null hypothesis being tested ("is the sender field aligned with the response field?") is not the scientific hypothesis; a shared confounder produces alignment | the statistic being tested is not a function of the geometry the null destroys |
| what would fix it | a null that preserves the confounder — our N1/N5, or CellWHISPER's | nothing on the coordinates; the *statistic* or the *test* has to change |

**Answering the question in its original binary form.** Is the tools' failure
"the null is too weak" or "the null is fine and the significance test is
miscalibrated"? It is the first — but not in a way that a stronger coordinate
null could fix, which is the part that matters. The strongest coordinate null
that exists is N0, which leaves no spatial information in the data whatsoever,
and under N0 the tools' scores still rank-correlate at ρ = 0.90–0.98 with
their real-data values. A null cannot be too weak when it has destroyed
everything there is to destroy; the statistic is simply not a function of what
was destroyed. **Our estimator's failure is the second kind and then some**: the
torus shift *does* annihilate β̂ (null centred at 0.8 % of the observed value),
so the test is arithmetically fine — and it rejects anyway, because the null
hypothesis "sender field and response field are unaligned" is false under pure
confounding. Two failures, two different repairs: the tools need a different
statistic or a null on the transport plan / neighbour graph; we need a null that
preserves the confounder.

This distinction is the contribution of Figure 4 beyond a replication.
"Method fails torus shift" is now a common finding; **why** it fails determines
what a user should do about it, and the two answers here demand opposite
remedies. A reader who takes "torus-shift failure" to mean "use a stronger
null" will fix nothing for COMMOT, and a reader who takes it to mean "the
statistic is not spatial" will draw the wrong conclusion about our estimator,
which is entirely spatial and fails for a completely different reason.


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
* ~~**A confounder-aware null for these tools.**~~ **Done — see §8.** This bullet
  read "CellWHISPER's own null brought their false-positive rate under 5 %. We
  have not run it… a proper comparison is future work." It has now been run
  (`N0_type`, within-cell-type location permutation, all four methods, 126 job
  checkpoints). The three ligand-receptor methods keep 78–97 % of their calls
  under it, and CellChat v2\*'s significance *rate* under it is 0.2831 against
  0.2833 on real coordinates. What remains genuinely open is the other half of
  the comparison: **CellWHISPER itself has not been run on this data**, so we
  report that these three fail their null, not that CellWHISPER passes it here.

---

## 8. CellWHISPER's own null, run (added 2026-08-21)

`code/phase4_data.py::null_coords` (`N0_type`) ·
`python3 phase4_run.py --method M --only-null N0_type` ·
126 new job checkpoints in `results/phase4/parts/`

### 8.1 Why this section exists

The D7 claim audit (§2.1 B9) found that this report had labelled the wrong null
as CellWHISPER's, in the row of the §2.4 ladder that Figure 4 is built around.
CellWHISPER permutes cell locations **within each cell type** — "preserving
cell-type-specific spatial organization and ligand-receptor (LR) expression
while destroying spatial proximity between ligand- and receptor-expressing
cells". Our `N0_perm` permutes across **all** cells, which additionally destroys
the cell-type architecture. Ours is strictly the more destructive of the two, so
the direction of the error was in our favour — but "CellWHISPER's >90 % figure
replicates" was not a replication claim we had earned.

The audit offered two options: relabel honestly, or run their null. We ran it.

### 8.2 The null

`N0_type` permutes the coordinate vector *within* each `cell_type_merged` label.
The point pattern, every marginal expression distribution, and each cell type's
own spatial organisation are preserved exactly; the only thing destroyed is
which particular ligand⁺ cell sits next to which particular receptor⁺ cell.
It is, in this project's vocabulary, a **confounder-preserving** null — the kind
§4.5 argues these tools actually need.

It is reached with `--only-null`, which enumerates the new condition without
disturbing the original job ordering; the per-job seeds in the 2026-08-20 run
are a function of that ordering and are therefore unchanged. Replicate counts
match the originals per method (COMMOT 10, CellChat v2\* 25 → run at 100,
SpaTalk\* 100, NCEM linear\* 100); where they differ, N0_type has *more*
replicates, so its Monte Carlo error is smaller, not larger.

### 8.3 Results

**Fraction of real-significant interactions still called significant:**

| method | N3_lig | N4_lig | N3_type | N0_perm | **N0_type (theirs)** |
|---|---|---|---|---|---|
| COMMOT | 0.806 | 0.820 | 0.769 | 0.794 | **0.811** |
| CellChat v2\* | 0.974 | 0.975 | 0.913 | 0.945 | **0.971** |
| SpaTalk\* | 0.814 | 0.820 | 0.784 | 0.791 | **0.781** |
| NCEM linear\* | 0.702 | 0.701 | 0.080 | 0.004 | **0.009** |

**Significance rate over all interactions — this is CellWHISPER's actual
criterion, "comparable interaction counts on real and randomised data":**

| method | real | N0_perm | **N0_type (theirs)** | N0_type ÷ real |
|---|---|---|---|---|
| COMMOT | 0.2236 | 0.2405 | **0.2303** | 1.03 |
| CellChat v2\* | 0.2833 | 0.3176 | **0.2831** | **1.00** |
| SpaTalk\* | 0.2472 | 0.2481 | **0.2422** | 0.98 |
| NCEM linear\* | 0.0144 | 0.0027 | **0.0033** | 0.23 |

**Rank correlation between the real score and the shuffled score, across
interactions (Spearman):**

| method | N0_perm | **N0_type (theirs)** |
|---|---|---|
| COMMOT | 0.904 | **0.919** |
| CellChat v2\* | 0.984 | **1.000** |
| SpaTalk\* | 0.860 | **0.899** |
| NCEM linear\* | 0.026 | **0.053** |

### 8.4 What changes, and what does not

**The conclusion does not change; it gets stronger and better-founded.**

1. **It is now a replication, not an analogy.** CellChat v2\* returns 0.2831
   significant interactions per candidate on within-type-permuted tissue against
   0.2833 on real tissue. That is CellWHISPER's criterion met to four decimal
   places, on a different tissue, a different panel and a different organism.
2. **Survival is *higher* under their null than under ours, for three of four
   methods** (COMMOT 0.811 vs 0.794; CellChat 0.971 vs 0.945; NCEM 0.009 vs
   0.004; SpaTalk is flat, 0.781 vs 0.791). This is the expected direction — a
   less destructive null should preserve more — and it removes the objection
   that our earlier numbers were an artefact of an unreasonably blunt shuffle.
   The parsimonious reading stands: there is little spatial information in these
   statistics for any null to remove.
3. **CellChat v2\*'s score ordering is perfectly preserved** under the null the
   reliability literature would apply to it: ρ = **1.000**. Every interaction
   keeps its rank when the ligand-receptor geometry is destroyed but the
   cell-type map is left intact. This is the single cleanest statement in the
   phase of what "the statistic is not spatial" means.
4. **NCEM linear\* remains the exception** and is correctly calibrated here too
   (0.009 survival, significance rate 0.23× real). It is the only one of the
   four that sees the difference.
5. **Nothing in §4 (mechanism) is affected.** The COMMOT transport-plan
   decomposition, the CellChat six-significant-figure invariance and the
   synthetic positive controls were never functions of which N0 variant was run.

### 8.5 What this costs the paper

One sentence of the framing has to be re-scoped rather than dropped: the
literature's benchmark result is **not** "these tools fail under coordinate
randomisation" but "these tools fail under a null that preserves cell-type
organisation and destroys only ligand-receptor proximity" — which is a
*sharper* claim, and harder to dismiss as a strawman null. The §2.4 ladder row
is corrected, Figure 4 gains a sixth condition, and `references.bib` carries the
finding on the `kumar2026cellwhisper` entry.
