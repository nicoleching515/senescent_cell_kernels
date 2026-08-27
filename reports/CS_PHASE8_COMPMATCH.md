# CS Phase 8 / PI decision D15 — the composition-matched rerun protocol at 5 seeds

**Status: implemented, run on the mouse arm, built but deliberately NOT run on the human arm.**

Producer: **`code/run_phase8_compmatch.py`** (new).
Outputs: **`results/phase3/compmatch_reruns.csv`** (per-seed rows + summary) and
**`results/phase3/compmatch_fits.csv`** (fit-level audit trail), plus the two
per-stage files `compmatch_{fits,reruns}_tierA.csv` and
`compmatch_{fits,reruns}_tierApm.csv` they were merged from.
Driver: `code/_compmatch_chain.sh`. Logs: `logs/compmatch_{tierA,tierApm,merge}.log`.

Nothing in `results/phase3/` was overwritten; nothing in `figures/` was written;
`code/run_phase3_nulls.py`, `code/summarize_phase3*.py`, `code/make_figure*.py`,
`genesets/` and `data/processed/deepscence_*` were not modified — the runner
**imports** `run_phase3_nulls` and calls its estimator. `data/raw_h1/` was never
opened.

---

## 0. Lead

**1. The protocol, as specified, is inert — and that is the finding, not a bug.**
Composition-matched decoys reproduce **1.6 %** of the naive pooled amplitude
(SF **0.9837**, 95 % block-bootstrap [0.973, 0.994]) and **3.6 %** of it within
receiver cell type (SF **0.9647**). Matching balances well —
max \|SMD\| **0.092 → 0.035**, median match rate **0.99987** (min 0.99918), the Master Plan §8 Test 5
gate (\|SMD\| ≤ 0.1) passes in **100 %** of matches — and removes almost nothing.
This is `CS_PHASE3.md` §5's warning about N2 reproducing exactly, on a matching
set chosen to be *only* composition.

**2. The five seeds buy essentially nothing, and the pre-registration should say
so out loud.** Across the five frozen seeds the pooled median SF moves over
**0.98370 – 0.98397** (sd **0.00012**); within cell type, **0.9629 – 0.9648**
(sd **0.00083**). At a median match rate of 0.99987 (minimum 0.99918 over all 6,237 fits)
the greedy matcher has almost no freedom left to randomise, so the seed is close to a formality. **Five seeds is
four more than this protocol needs.** Freezing it as a five-seed protocol is
defensible only as an audit device.

**3. Report against interest: composition explains MORE of the naive gradient
than the §17 "66–76 %" range says.** Conditioning the pooled fit on the
receiver's own cell type removes **65.9 %** (SF 0.3414 [0.236, 0.402]) — which
reproduces the published 66 % almost exactly. But adding the 20-NN cell-type
composition vector on top takes it to **85.4 %** (SF 0.1461 [0.052, 0.246]),
**above the top of the published range**, and to **98.5 %** for
`downstream_arrest`. On the per-module Tier A sensitivity sets the same quantity
is **79.6 %**, and `downstream_arrest` goes to **114.6 %** — i.e. the pooled
gradient there is *entirely* composition and then some. The honest statement is
**66 % (own cell type) to 85 % (own cell type + neighbourhood composition)**,
not 66–76 %.

**4. So: does the composition confound explain the kernel? Yes, most of it — and
the matched rerun is not the instrument that shows it.** The comparison is
like-for-like: the `comp` matching set stratifies **exactly** on receiver cell
type and matches on the 20-NN composition vector, so `typecomp_adj` — cell-type
intercepts plus that same composition vector, entered as covariates — is its
exact regression counterpart. **The same variables, on the same fits, remove
1.6 % as a matched decoy set and 85.4 % as covariates: a factor of fifty.** The
regression is the one that is right about the confounding, because the
confounding acts through the *receiver's* dependence on composition, and
balancing sender against decoy does not touch that (`CS_PHASE3.md` §5, verbatim:
"matching balances the covariates between senders and decoys; it does not remove
the dependence of the response on those covariates at the receiver, which is
where the confounding acts"). Anyone reading only the matched number would
conclude composition is irrelevant.

**5. A §17 row that had no producer now has one.** `CS_PHASE8_M1_RERUN.md` §7
records that "Composition surrogate share 66–76 %" is emitted by no script in the
repository. `results/phase3/compmatch_reruns.csv` now emits it, four ways, with
its scope, its estimator and its reportable population stated. §5 proposes the
replacement wording.

**6. Validation.** On the 142 (section × receiver type × module) cells reportable
in both, `beta_naive` and `lam_naive` from this runner are **bit-identical** to
`results/phase3/main_fits.csv` (max \|Δβ\| = 0.000e+00; λ̂ equal in every cell),
and the published-N2 matching set rerun at the five seeds gives median SF
**0.9463 – 0.9499** against `main_fits.csv`'s **0.9490** on the same cells. The
harness is calling the same estimator on the same data.

---

## 0b. The numbers

Mouse arm M1, six Section 8 Test 3 admissible sections, post-C6 gene sets,
window 100 µm, λ fixed at λ̂_naive, 400 block-bootstrap replicates over 100
quantile blocks. Reportable population = positive naive amplitude whose
block-bootstrap CI excludes zero. All rows from
**`results/phase3/compmatch_reruns.csv`** (`row_type = summary`); per-fit detail
in **`results/phase3/compmatch_fits.csv`** (6,237 rows).

### PRIMARY sender set — `A_SENDER_FINAL_strict` (33 genes), call `tierA_p95`

| Control | scope | n fits | **SF** | across-seed [min, max] (sd) | **share removed** (IQR) | median 95 % CI on SF |
|---|---|---|---|---|---|---|
| **`comp` — composition-matched decoys (THE PROTOCOL)** | **pooled** | 165 | **0.9837** | [0.98370, 0.98397] (0.00012) | **1.6 %** [0.9, 3.3] | [0.973, 0.994] |
| `comp` | by cell type | 748 | 0.9647 | [0.9629, 0.9648] (0.00083) | 3.6 % [1.2, 8.5] | [0.895, 0.997] |
| `full` — published N2 matching set, same 5 seeds | pooled | 165 | 0.9855 | [0.98539, 0.98577] (0.00016) | 1.4 % [0.9, 3.7] | [0.979, 0.992] |
| `full` | by cell type | 748 | 0.9488 | [0.9471, 0.9505] (0.00153) | 5.1 % [2.2, 8.0] | [0.873, 0.991] |
| `comp_adj` — same variables, as covariates | pooled | 33 | 0.4989 | seed-free | **50.1 %** [33.0, 68.3] | [0.421, 0.606] |
| `comp_adj` | by cell type | 150 | 0.5023 | seed-free | 49.8 % [19.6, 71.6] | [0.105, 0.771] |
| `type_adj` — receiver cell-type intercepts | pooled | 33 | **0.3414** | seed-free | **65.9 %** [56.7, 80.2] | [0.236, 0.402] |
| `typecomp_adj` — cell type + 20-NN composition | pooled | 33 | **0.1461** | seed-free | **85.4 %** [72.6, 92.8] | [0.052, 0.246] |

### Per-module sensitivity sender sets — `A_sender_for_<module>.txt`, call `tierApm_p95`

| Control | scope | n fits | **SF** | across-seed [min, max] (sd) | **share removed** (IQR) |
|---|---|---|---|---|---|
| **`comp` (THE PROTOCOL)** | **pooled** | 172 | **0.9861** | [0.98589, 0.98625] (0.00013) | **1.4 %** [0.8, 2.8] |
| `comp` | by cell type | 790 | 0.9660 | [0.9656, 0.9682] (0.00106) | 3.4 % [0.8, 6.8] |
| `comp_adj` | pooled | 34 | 0.5140 | seed-free | 48.6 % [29.9, 65.4] |
| `comp_adj` | by cell type | 157 | 0.5693 | seed-free | 43.1 % [16.9, 66.3] |
| `type_adj` | pooled | 34 | 0.3350 | seed-free | **66.5 %** [58.5, 78.1] |
| `typecomp_adj` | pooled | 34 | 0.2042 | seed-free | **79.6 %** [71.8, 90.2] |

**The two pre-registered Tier A variants agree.** Every paired row differs by
less than 0.07 in SF, and the protocol's own number (`comp`, pooled) differs by
0.0023. The conclusion does not depend on which Tier A set is used.

### Per-seed rows, the protocol, PRIMARY call

| seed | scope | fits | reportable | SF | share | match rate | max \|SMD\| after |
|---|---|---|---|---|---|---|---|
| 20260901 | pooled | 42 | 33 | 0.983695 | 1.63 % | 0.999871 | 0.0357 |
| 20260902 | pooled | 42 | 33 | 0.983887 | 1.61 % | 0.999871 | 0.0340 |
| 20260903 | pooled | 42 | 33 | 0.983965 | 1.60 % | 0.999871 | 0.0357 |
| 20260904 | pooled | 42 | 33 | 0.983714 | 1.63 % | 0.999871 | 0.0356 |
| 20260905 | pooled | 42 | 33 | 0.983735 | 1.63 % | 0.999871 | 0.0353 |
| 20260901 | by cell type | 315 | 150 | 0.964668 | 3.53 % | 0.999887 | 0.0340 |
| 20260902 | by cell type | 315 | 149 | 0.962919 | 3.71 % | 0.999887 | 0.0315 |
| 20260903 | by cell type | 315 | 150 | 0.964666 | 3.53 % | 0.999887 | 0.0321 |
| 20260904 | by cell type | 315 | 149 | 0.964825 | 3.52 % | 0.999887 | 0.0327 |
| 20260905 | by cell type | 315 | 150 | 0.963679 | 3.63 % | 0.999887 | 0.0333 |

The only thing the seed moves visibly is **which fits are reportable** (149 vs
150 of 315), because `beta_base_lo` comes from the seeded bootstrap. The
estimate itself is stable in the third decimal place.

### Where the naive amplitude is, for reference

Pooled (unstratified) naive amplitude, PRIMARY call, reportable fits:
**0.643 response-sd** [IQR 0.483, 0.774], λ̂ railed in **76 %** of fits.
Within receiver cell type the same runner gives **0.326–0.333 response-sd** per
seed, reproducing `main_fits.csv`'s 0.329 and §17's 0.326. The pooled amplitude
is roughly twice the within-type amplitude — which is the composition effect,
stated as an amplitude instead of as a fraction.

### Heterogeneity behind the medians

Pooled composition share by response module, PRIMARY call
(`oxidative_stress` never enters: its pooled naive amplitude is negative in all
six sections, median −0.071):

| module | `type_adj` | `typecomp_adj` | `comp_adj` | **`comp` (matched)** |
|---|---|---|---|---|
| downstream_arrest | 0.853 | **0.985** | 0.536 | 0.013 |
| emt_ecm | 0.728 | 0.912 | 0.720 | 0.009 |
| il6_jak_stat3 | 0.619 | 0.710 | 0.334 | 0.017 |
| interferon_response | 0.668 | 0.836 | 0.334 | 0.023 |
| secondary_senescence | 0.714 | 0.911 | 0.633 | 0.016 |
| tnfa_nfkb_proximal | 0.650 | 0.743 | 0.470 | 0.018 |

By section (`type_adj` / `typecomp_adj`, pooled, PRIMARY): 7001 0.522 / 0.726,
7248 0.533 / 0.632, 7259 0.567 / 0.761, 7260 0.668 / 0.931, 7352 0.783 / 0.903,
7435 0.771 / 0.883. **The composition share is never below 0.52 in any section
and never below 0.59 in any module.**

---

## 1. The specification I found, and where

This is the first thing the task asked for, so it comes first, in full, with the
line references.

### 1.1 What the planning documents actually say

| Where | Text |
|---|---|
| `Phase7_Minimal_Human_Replication (1).md` **§15**, line 421 | "**Fixed and not retunable:** … the composition-matched rerun protocol at 5 seeds." |
| same, **§16 run order step 9**, line 449 | "**Composition-matched reruns**, 5 seeds, both arms." |
| same, **§20 day 9**, line 532 | "Composition-matched reruns, 5 seeds, both arms. **§8 comparison**" |
| same, **§17**, line 475 | table row "\| Composition surrogate share \| 66–76% \| \|" |
| `reports/PHASE8_ROADMAP_STATUS.md` line 126 | roadmap item 10.2, "Composition-matched reruns, 5 seeds, both arms" |

**That is the entire specification.** Five mentions, no method. No matching
variables, no matching rule, no seed values, no estimand, no output format, and
no statement of which of "composition matching" and "composition adjustment" is
meant. `code/composition_all.py` is unrelated (a descriptive per-section
composition table; no matching, no seeds, no reruns), exactly as
`reports/PREREG_PHASE8.md` §3.8 records.

### 1.2 What the surrounding documents constrain

These are not the protocol, but they fix everything the protocol must be
consistent with, and I reconstructed it from them rather than inventing one.

| Constraint | Source |
|---|---|
| The confound: k-NN cell-type composition is a Tier D nuisance covariate in the model of §6.1 | `SASP_Kernel_Master_Plan.md` line 224 (§6.1), line 487 (Tier D) |
| The only sender-vs-non-sender matching design in the plan is **N2**: "a non-senescent cell matched on cell type, local density (k-NN within 50 μm), and k-NN composition" | Master Plan §22 Step 3 N2, line 1011; §8 Test 5, line 366 |
| Matching rule: "nearest-neighbor matching on a propensity score, or coarsened exact matching" | Master Plan line 366 |
| Balance criterion: \|SMD\| ≤ 0.1 | Master Plan §8 Test 5; roadmap item 9.3 |
| Estimand: a **surviving fraction**, λ held fixed at λ̂_naive, never a permutation p-value | Master Plan §6.5; `CS_PHASE3.md` §0, §5 |
| **400 bootstrap replicates over 100 quantile blocks** — itself a §15 frozen parameter | Phase7 §15, line 421 |
| Reportable population for any SF median: positive naive amplitude whose block-bootstrap CI excludes zero | `code/summarize_phase3.py` docstring + `sf_table` |
| What "composition surrogate share 66–76 %" means: **66 %** = 1 − 0.344, the SF when receiver cell-type intercepts alone are added (`CS_PHASE3.md` §3.1); **76 %** = 0.212 / 0.260, the composition-only surrogate curve's share of the unstratified contact amplitude (`CS_PHASE5.md` §4). Two different measurements, on two different scales, reported as one range | `CS_PHASE3.md` §3.1; `CS_PHASE5.md` §4; `README.md` L170–179 |
| Matching in this dataset is known to be **nearly inert**: N2 balances to \|SMD\| ≤ 0.033 and still absorbs only ~6 %, while regressing on the same covariates absorbs 92 % | `CS_PHASE3.md` §5 N2, lines 285–306 |
| And N2 is inert **by calibration**, not by accident: on synthetic tissue it returns 0.934 with a planted real effect and 0.775 with none | `CS_PHASE3.md` §6, line 435 |

### 1.3 Where the documents are silent, and what I chose

Every one of these is a decision, not a reading. They are written into the
producer's docstring as `D15.1`–`D15.9` so the pre-registration can transcribe
them, and `reports/PREREG_PHASE8.md` §3.8's four `TBD` fields are filled in §6
below.

| # | Ambiguity | Choice | Why |
|---|---|---|---|
| **D15.1** | "Composition-matched" — matched *how*, and matching *what to what*? | The **N2 design** — senders vs non-sender decoys — with the matching set reduced to composition. | N2 is the only sender/non-sender matching design in the plan, it is already implemented, and the task's own framing ("match senders and non-senders on cell-type composition, re-fit") is N2 with a reduced matching set. |
| **D15.2** | Which variables? | Exact stratification on receiver **cell type**, plus 1-1 nearest-neighbour propensity matching without replacement, caliper **0.25 SD**, on the **20-NN cell-type composition vector** (`knn_frac_*`, one column per cell type in the section). Nothing else. | The point is to isolate composition. Leaving density, depth and anatomy in makes it the published N2 again. The published N2 matching set is run alongside as variant `full` so the isolation is measured, not asserted. |
| **D15.3** | Which five seeds? | **20260901, 20260902, 20260903, 20260904, 20260905.** | Derived by date from `sasp_phase3.MASTER_SEED = 20260820`, and chosen outside the range `run_phase3_nulls._expand` can reach (`MASTER_SEED + 1000·i + j`, and `+ PM_SEED_OFFSET`) so no rerun accidentally reproduces an existing N2 match. |
| **D15.4** | What does the seed control — matching, bootstrap, or both? | Both, but they are separable: the seed is passed to `match_decoys_section` (it fixes the order in which senders claim decoys in the greedy matcher) and also seeds the block bootstrap. **Point estimates are computed before the bootstrap**, so the between-seed spread reported here is matching variability and nothing else. | Otherwise "at 5 seeds" would measure bootstrap noise, which is already reported as a CI. |
| **D15.5** | Pooled or per receiver cell type? | **Both.** `scope = ALL` (pooled, unstratified) is primary — it is the fit the field reports and the one §3.1/§5 §4 measured composition against, so it is the row that fills the §17 cell. Per-cell-type rows are secondary. | A within-cell-type fit has already removed the between-type composition effect, so it cannot be the fit the §17 row is about. |
| **D15.6** | Which sender call? | `tierA_p95` on `A_SENDER_FINAL_strict` (33 genes, PRIMARY) **and** `tierApm_p95` on the seven `A_sender_for_*.txt` per-module sets. | The task requires both pre-registered Tier A variants. |
| **D15.7** | Which sections? | The six Section 8 Test 3 admissible ("in band") sections, `sasp_phase3.IN_BAND`. | Every other Phase 3 primary uses them. |
| **D15.8** | Which fits enter a median? | `beta_naive > 0` **and** `beta_base_lo > 0`, identical to `summarize_phase3.sf_table`. | So the medians here are directly comparable to the published SF table. |
| **D15.9** | Matching alone cannot answer the scientific question. | Three seed-free **covariate-adjustment** variants are computed at the same scopes, same λ̂, same bootstrap: `comp_adj` (20-NN composition as covariates), `type_adj` (receiver's own cell-type intercepts — the μ_{c_i} of §6.1, the block behind the "66 %"), `typecomp_adj` (both). | `CS_PHASE3.md` §5 is explicit that in this dataset matching balances covariates without removing the response's dependence on them. Reporting only the matched result would say "composition matching removes nothing" and a reader would hear "composition is not the confound". That would be the opposite of the truth. |

**Recorded silence.** The documents never say what the protocol's *output* is, so
the §17 row it is supposed to fill ("Composition surrogate share") had — before
this run — **no producer at all**; `reports/CS_PHASE8_M1_RERUN.md` §7 independently
reached the same conclusion ("no script emits a quantity in that range under that
name"). §5 below proposes what should go in that cell and how it must be labelled.

---

## 2. What was reused, and what is new

**Reused, unmodified, by import or subclass:**

| Component | File | Used for |
|---|---|---|
| `match_decoys_section` / `greedy_ps_match` | `code/phase3_core.py:26-107` | the matcher itself. Only its `Zmatch` argument changes. |
| `build_blocks` | `code/phase3_core.py:120-160` | the N5 covariate block, the `Zmatch` matching set and its column names — the composition columns are taken from it by name, never rebuilt |
| `SectionFit` | `code/run_phase3_nulls.py` | per-section setup: covariates, N6 neighbour baseline, quantile block ids, λ grid, receiver masks |
| `fit_cell` | `code/run_phase3_nulls.py` | **called unmodified.** Naive fit, shared-λ two-kernel decoy contrast, nested N5/N6/anatomy designs, and the 400-replicate spatial block bootstrap |
| `BlockProfiler` | `code/sasp_estimators.py:388-517` | the sufficient-statistic profiler and `beta_at` / `beta2_at` fixed-λ readouts |
| `dist_to_senders`, `block_ids`, `IN_BAND`, `MODULES`, `EXCLUDE_TYPES` | `code/sasp_phase3.py` | geometry and inventory |
| `run_phase3_attribution.py`'s `comp` sub-block definition | `code/run_phase3_attribution.py:38` | the `comp_adj` covariate set, extended to the pooled scope |
| `summarize_phase3.sf_table`'s reportable-population rule | `code/summarize_phase3.py` | reimplemented in two lines rather than imported, because importing that module reads `main_fits.csv`, which the concurrent M1 re-run was rewriting |

**New, and it is deliberately small:**

* `CompMatchFit(run_phase3_nulls.SectionFit)` — adds exactly one method,
  `rematch(variant, seed)`, which re-runs the same matcher with a different
  matching set and seed and refreshes the decoy distance field `fit_cell` reads.
* `adjust_rows(...)` — the D15.9 covariate-adjustment counterparts. One
  `BlockProfiler` per adjustment yields both the naive (`p = 1`) and the adjusted
  (`p = full`) fit, because `_acc(m, p)` slices `X[:, :p]`.
* `summarise(...)` — per-seed and across-seed aggregation.

No estimator, matcher, kernel, bootstrap or window rule is reimplemented.

**Column renaming.** `fit_cell` writes the decoy contrast under `n2`. Under
`rematch` those columns are no longer the published N2, so they are renamed
`*_matched` on the way out (`beta_n2 → beta_matched`, `sf_n2 → sf_matched`, and
so on). Nothing in the output can be confused with `main_fits.csv`'s `sf_n2`.

---

## 3. What was run

Post-C6 inputs throughout: `data/processed/{modules,senders}_*.csv` were rebuilt
by the M1 re-run at 05:48 UTC and `data/processed/cache3/*.npz` at 05:54 UTC;
this run started after both. System `python3`; the mid-build venv at
`/workspace/envs/sasp311` was not used. **No package was installed.**

```bash
# stage A — PRIMARY Tier A (A_SENDER_FINAL_strict, 33 genes)
python3 -u code/run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierA_p95 \
    --variants comp,full,comp_adj,type_adj,typecomp_adj --out-tag _tierA
# stage B — the seven per-module Tier A sensitivity sets
python3 -u code/run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierApm_p95 \
    --variants comp,comp_adj,type_adj,typecomp_adj --out-tag _tierApm
# merge into the deliverable
python3 -u code/run_phase8_compmatch.py --merge \
    results/phase3/compmatch_fits_tierA.csv,results/phase3/compmatch_fits_tierApm.csv
```

Frozen parameters carried through unchanged, none re-derived here: window
100 µm; λ grid 40 log-spaced points over [7 µm, 50 µm]; λ held fixed at
λ̂_naive inside every control; 400 bootstrap replicates over 100 quantile
blocks; `MIN_RECEIVERS = 2000`; caliper 0.25 SD.

---

## 4. The human arm

**H1 was not run and no H1 file was read.** The runner is arm-generic — an arm is
a section list, a cache directory, a label set and the two Tier A calls
(`ARMS` in the producer) — and `--arm h1` exits with:

```
REFUSING to run arm 'h1'.
  Phase7 s15: H1 is behind the pre-registration freeze and Phase 8's tag is not cut.
  The protocol is implemented and arm-generic; running it on H1 is roadmap item 10.2, after Phase 9.
  To run it after the freeze: populate ARMS['h1'] (cache + sections) and set SASP_H1_UNFROZEN=1.
```

Two conditions must both hold before it will run: `ARMS['h1']['sections']` must
be populated (it is empty, and the guard refuses an empty list separately), and
`SASP_H1_UNFROZEN=1` must be set. Running it on H1 in Phase 10 is then the same
command with `--arm h1`.

**What H1 will need that M1 has.** The protocol reads nothing mouse-specific of
its own; everything comes through `phase3_core.build_blocks`, which needs the H1
`Sec` cache to carry the arm-specific anatomical covariate in the `zonation_score`
slot (Phase7 §15's declared deviation: lung/spleen has no zonation) and a
`seg_code`. Neither affects the `comp` matching set, which is composition only —
so the **primary** variant of this protocol is the one least exposed to the
cross-arm covariate deviation, which is a point in its favour for a two-arm table.


---

## 5. What the §17 row should say, and what §15 should freeze

### 5.1 The "Composition surrogate share" row

The published range **66–76 %** is two different measurements on two different
scales reported as one interval: **66 %** is `1 − SF` on β̂ for receiver
cell-type intercepts (`CS_PHASE3.md` §3.1), **76 %** is a ratio of two *binned
curve amplitudes* — the composition-only surrogate against the unstratified
curve (`CS_PHASE5.md` §4). Neither is emitted as a number by any script: the
66 % had no file behind it at all (this run now produces it), and the 76 % is
derivable from `figures/figure2a_stratified_curves.csv` (`kind =
compositiononly` vs `kind = unstrat`, `celltype = all`) but is not computed
there. The second is not the estimand this protocol makes, so the row should be
split rather than patched:

| §17 row | Proposed M1 value | Source |
|---|---|---|
| **Composition share of the pooled naive amplitude — receiver cell type** | **65.9 %** (SF 0.3414 [0.236, 0.402]), n = 33 | `compmatch_reruns.csv`, `variant=type_adj, scope_kind=pooled, call=tierA_p95` |
| **Composition share — cell type + 20-NN neighbourhood composition** | **85.4 %** (SF 0.1461 [0.052, 0.246]), n = 33 | same file, `variant=typecomp_adj` |
| **Composition-matched rerun, 5 seeds (surviving fraction)** | **0.9837**, across-seed range 0.98370–0.98397 | same file, `variant=comp, scope_kind=pooled` |
| *(unchanged, different estimand)* composition-only surrogate **curve** amplitude, 76 % | not reproduced here — a curve-amplitude ratio, not an SF | `figures/figure2a_stratified_curves.csv` (`kind=compositiononly`), `figures/figure2a_amplitudes.csv` (`ALL RECEIVERS` = 0.2600), `CS_PHASE5.md` §4 |

**The "66–76 %" interval should not be carried into the pre-registration as
written.** Its lower end is reproduced exactly; its upper end is a different
estimator; and the quantity it was meant to bracket is **85 %**, above it.

### 5.2 Transcription block for `PREREG_PHASE8.md` §3.8

The four `TBD pending reports/CS_PHASE8_COMPMATCH.md` fields, filled:

| Field | Value |
|---|---|
| Producer script | **`code/run_phase8_compmatch.py`** (driver `code/_compmatch_chain.sh`) |
| Matching variables | Exact stratum: **receiver cell type**. Propensity score: the **20-NN cell-type composition vector** (`knn_frac_<type>`, one column per cell type present in the section, from `phase3_core.build_blocks`), plus cell-type dummies in the propensity model. **No** density, depth, anatomy or segmentation term — that is the `full` comparison variant, not the protocol |
| Matching rule | 1-1 nearest-neighbour on the propensity logit, **without replacement**, **within cell type and within section**, **caliper 0.25 SD** of the within-stratum score sd; `phase3_core.greedy_ps_match`. Estimator: shared-λ two-kernel sender-vs-decoy contrast (`BlockProfiler.beta2_at`) at λ = λ̂_naive; SF = β_sender / β_naive; 400 block-bootstrap replicates over 100 quantile blocks |
| Number of seeds | **5**, per §15 |
| The five literal seed values | **20260901, 20260902, 20260903, 20260904, 20260905** |
| Arms | Both. M1 run 2026-08-27; H1 built and gated behind `SASP_H1_UNFROZEN`, to be run as roadmap item 10.2 |
| Sections | `sasp_phase3.IN_BAND` (six Section 8 Test 3 admissible sections) |
| Sender calls | `tierA_p95` (PRIMARY) and `tierApm_p95` (seven per-module sets) |
| Scope | Pooled (primary, fills the §17 row) and per receiver cell type (secondary) |
| Reportable population | `beta_naive > 0 & beta_base_lo > 0`, as `summarize_phase3.sf_table` |
| Output | `results/phase3/compmatch_reruns.csv` (per-seed + summary), `results/phase3/compmatch_fits.csv` (fit level) |

### 5.3 A recommendation the PI should see before tagging

Freezing "the composition-matched rerun protocol at 5 seeds" is now possible and
this report supplies every value needed to do it. But the run shows the frozen
parameter is **nearly content-free as a scientific control**: at a median match rate of
0.99987 the seed moves the estimate by 1.2 × 10⁻⁴ (pooled), and the estimate itself
is 0.98, so the protocol certifies as "not composition" a gradient that is 66–85 %
composition by the regression on the same variables. Two options, both honest:

1. **Freeze it as written** (it is now implementable and reproducible), and state
   in the same breath that it is a **matching-balance audit**, not evidence about
   the confound — with the `type_adj` / `typecomp_adj` numbers reported beside it
   every time it appears. This is what §0b does.
2. **Freeze it and add the covariate-adjustment counterpart to the frozen list**,
   so the §15 item reads "the composition-matched rerun protocol at 5 seeds *and
   its covariate-adjusted counterpart*". This is the version I would argue for,
   because §17's row is about how much of the gradient is composition, and only
   the second design answers that.

Either way, `CS_PHASE3.md` §6's calibration is the reason: on synthetic tissue
N2-style matching returns **0.934 with a planted real effect and 0.775 with no
effect at all**. A design whose two answers are that close cannot discriminate,
and a composition-only matching set has strictly less to work with than N2 did.

---

## 6. Verification performed

| Check | Result |
|---|---|
| Same estimator, same data as the published pipeline | `beta_naive` and `lam_naive` **bit-identical** to `results/phase3/main_fits.csv` on all 142 (section × receiver type × module) cells reportable in both. `max abs diff = 0.000e+00`; λ̂ equal in every cell |
| Published N2 reproduced | `full` variant, like-for-like on those 142 cells: seed medians **0.9463, 0.9469, 0.9496, 0.9483, 0.9499** against `main_fits.csv`'s **0.9490**. The published value lies inside the five-seed range |
| Published naive amplitude reproduced | within-cell-type median **0.326–0.333** response-sd per seed vs `main_fits.csv` **0.3288** and §17's **0.326** |
| Published composition sub-block reproduced | `comp_adj`, by cell type, **0.502** against `results/phase3/attribution.csv`'s `sf_comp` on the same reportable rule, **0.4996** (and `CS_PHASE3.md` §3.2's **0.474** on the looser `beta_naive > 0` rule). `attribution.csv` is **pre-C6** (mtime 08-20) so exact equality was not expected |
| Published cell-type-intercept share reproduced | `type_adj` pooled **0.3414** against `CS_PHASE3.md` §3.1's **0.344** |
| Matching balance gate (Master Plan §8 Test 5) | \|SMD\| ≤ 0.1 in **100 %** of matches, both variants, both calls. Worst case over all 6,237 fits: **0.0759** (`tierApm_p95`, `comp`); median match rate 0.99987, minimum 0.99918 |
| Nothing in `results/phase3/` overwritten | Only six new files, all named `compmatch_*`; the runner **refuses** to write over an existing path and requires `--out-tag` |
| Figures untouched | `python3 code/check_figures_guard.py` → `OK: all 27 committed figures match (PDF date stamps ignored)`, before and after |
| Forbidden files unmodified | This task added exactly nine paths: `code/run_phase8_compmatch.py`, `code/_compmatch_chain.sh`, `reports/CS_PHASE8_COMPMATCH.md`, and the six `results/phase3/compmatch_*.csv`. `code/run_phase3_nulls.py` (mtime 05:52), `code/summarize_phase3.py` (06:30), `code/phase3_core.py` (08-20), `code/sasp_estimators.py` (08-20) and `genesets/A_SENDER_FINAL_strict.txt` (05:41) all predate this run (07:01–08:46) and carry no edit of mine |
| H1 | Not run. `data/raw_h1/` never opened. `--arm h1` exits on the guard |
| Environment | System `python3`. **No package installed.** `/workspace/envs/sasp311` not used |

## 7. What I did NOT do

* **No H1 run** — by instruction, and the guard is in the code so it cannot happen
  by accident.
* **No figure.** A three-panel figure (SF by control × scope; the five-seed
  spread; composition share by module) would carry §0b better than the tables,
  but `figures/` is frozen and `figures/revised_candidates/` is another agent's
  area this shift.
* **No `oxidative_stress` result at pooled scope.** Its pooled naive amplitude is
  negative in all six sections, so it never enters a reportable population. That
  is a property of the module, not of this protocol.
* **No zonation-stratified rows.** `run_phase3_nulls`'s `strata_mode="zonation"`
  path was not exercised; the protocol is defined on the `stratum = all` fits.
* **No donor bootstrap.** Section-level spatial block bootstrap only, as
  everywhere else in Phase 3.
