# CS Phase 2 — Naive estimation on real tissue (Figure 2a), plus the Phase 1 work still owed

**Status: complete.** Section 22 Step 2 delivered on 5 liver sections; the three
items owed from Phase 1 (kernel-family misspecification, the superposition variant,
and nulls N1/N3/N4) are also delivered.

| Deliverable | Location |
|---|---|
| Real-tissue loader, swappable annotations | `/workspace/code/sasp_real.py` |
| One-time per-section cache | `/workspace/code/prepare_samples.py` |
| All five Section 6.2 kernel families + AIC + CIs | `/workspace/code/sasp_kernels.py` |
| Real-data fit runner (checkpointed) | `/workspace/code/run_figure2a.py` |
| Nulls N1 / N3 / N4 | `/workspace/code/sasp_nulls.py` |
| Misspecification + superposition + nulls driver | `/workspace/code/sasp_phase1b.py` |
| **Figure 2a** | `/workspace/figures/figure2a.png`, `.pdf` |
| Figure 2a underlying data | `/workspace/figures/figure2a_curves.csv`, `figure2a_fits.csv` |
| Real fits (210 rows) and curves (1,050 rows) | `/workspace/results/real_fits.csv`, `real_curves.csv` |
| Phase 1b results | `/workspace/results/phase1b/{misspec,nulls,nulls_reprofiled_only}.csv` |

---

## Headline

The naive curve on real liver tissue is **exactly as clean and monotone as the
published gradients** — and essentially none of the quantities the field reports
from it are identified. λ̂ hits a grid bound in 66 % of fits, the choice of kernel
family moves d̂½ by a median factor of 5, the standard iid CI is 4–5× too narrow
(280× for the step family), and adjusting for transcript count, cell area and local
density alone removes 86 % of the effect. Separately, on synthetic data with a known
answer, **the torus-shift null the field recommends certifies pure confounding as a
real effect 98 % of the time.**

## Methods (Phase 2)

### Data
Xenium Prime Mouse 5K liver sections, one section per mouse. Coordinates are
already in **microns** (`x_centroid`, `y_centroid` in `cells.parquet`); the unit is
carried in every column name, per Section 8 Test 1.

| section | condition | week | cells (QC-passed) | median NN (µm) |
|---|---|---|---|---|
| 7239 | sbr | 52 | 83,392 | 8.50 |
| 7250 | sham | 26 | 236,906 | 9.67 |
| 7259 | sbr | 26 | 127,386 | 6.74 |
| 7352 | sham | 2 | 139,378 | 9.80 |
| 7361 | sbr | 2 | 193,984 | 8.77 |

Features: of 13,590 rows in `cell_feature_matrix.h5`, only the **5,106 Gene
Expression** features are used. Genomic Control (21), Negative Control Codeword
(609), Negative Control Probe (40), Unassigned Codeword (7,806) and Deprecated
Codeword (8) are dropped, as is any Gene Expression feature whose name matches a
genotyping-probe pattern (none matched on this panel). Cells with fewer than 20
transcripts are dropped.

**Median NN distance ranges from 6.74 to 9.80 µm across sections — a 45 % spread in
packing.** This is the resolution floor of Section 8 Test 1: no length scale below
~7–10 µm is reportable, and any cross-section contrast is confounded by packing
before any biology. Every length scale is therefore also reported in units of the
local median NN distance.

### Swappable annotation inputs
Sender calls, cell-type labels, anatomy and module scores are **inputs, not
hardcoded**. `sasp_real.load_sample` reads `celltypes_/anatomy_/senders_/modules_
{sample}.csv` from `/workspace/data/processed/` when present and otherwise makes a
clearly labelled provisional call, recorded in the `*_source` columns of every
results row so no figure can silently imply an annotation that did not happen.

At the time of this run the Bio agent's files were not yet on disk, so:
* **modules** — provisional: mean per-gene z-score (CP10K, log1p) over each Tier B
  set from `/workspace/genesets/B_*.txt`;
* **senders** — provisional: `tierA_p95`, the 95th percentile of the Tier A score
  built from `A_SENDER_FINAL_strict.txt` (25 genes). This gives exactly 5 %
  prevalence, inside the plan's 2–10 % sweet spot. `Cdkn1a > 0`, the simplest call
  and the one the source paper uses, was also computed and gives **0.68 %**
  prevalence on section 7250 — *below* the 1 % floor of Section 8 Test 3, so it is
  recorded but not used as the primary call;
* **cell types** — not available, so all fits are on a single unstratified
  stratum labelled `ALL`. Figure 2a is therefore **not** yet per receiver cell
  type; the code path for stratification exists and activates automatically.

### Pipeline
`prepare_samples.py` parses the h5 once per section, scores modules, calls senders,
builds cKDTree geometry (density at 25/50/100 µm, 20-NN indices, median NN
distance, distance to nearest sender for every sender call) and caches everything
to a ~5–15 MB `.npz`. Parsing costs 30–75 s per section; the kernel fits cost ~1 s,
so caching is what makes the fitting stage cheap and makes re-running after the Bio
annotations land a single cheap re-prepare. **cKDTree throughout; no (n,n) matrix
is ever formed** — a 238k-cell section would be 2.8 × 10^10 pairs.

### Estimation
Binned mean response in 10 µm bins to 300 µm (Section 22 Step 2), using the in-bin
mean distance as the regressor. All five Section 6.2 families are fitted by
profiling the nonlinear parameters on a grid and solving the linear parameters
exactly at each grid point: exponential, Gaussian, power law (2-D grid over λ and
p), step, and a cubic B-spline reference. Families are made comparable by reporting
**d½**, the distance at which the fitted kernel falls halfway from its value at
d = 0 to its far-field level.

Every fit reports **both** the iid asymptotic Gauss–Newton CI (what a standard NLS
package gives, and what the field quotes) and a **spatial block bootstrap** CI over
a 10 × 10 grid of quantile-defined spatial blocks, 300 replicates. Phase 1 measured
an SE understatement factor up to 7.9× for the iid CI under spatial confounding, so
reporting the iid CI alone here would be indefensible.

Three additional fits per section × module:
* distance rescaled to **local median NN units** (packing-invariant view);
* a **covariate-adjusted preview of null N5**, adding log transcript count, log
  cell area, density at 25/50/100 µm and the 1-NN distance to the design;
* the spline under the same adjustment.

All seeds pinned from `MASTER_SEED = 20260820` plus hashes of section, module and
stratum.

## Results

### 1. The naive gradient is exactly what the field reports — clean and monotone

Across 35 section × module combinations, the binned mean response falls monotonically
with distance to nearest sender: **Spearman ρ < 0 in 32/35**, and ρ ≤ −0.92 in 24/35.
The effect is large relative to sampling noise — median |Δ| between the first three
and last six bins is **0.056 module-score units against a median bin SEM of 0.0020,
a ratio of 26**. Every module shows it, in every section except 7239 (sbr, 52 wk, the
smallest at 83k cells), where three of seven modules are flat or reversed.

**This is Figure 2a and it looks precisely like the published gradients.** If we
stopped here we would report a clean SASP distance gradient in mouse liver.

### 2. λ is not identified, in two-thirds of fits

Fitting the exponential kernel on a log grid spanning 3–400 µm:

* **23 of 35 fits (66 %) land on a grid bound.** 8 rail to the 3 µm lower bound —
  *below the 6.7–9.8 µm median nearest-neighbour distance*, i.e. below the resolution
  floor Section 8 Test 1 defines. 15 rail to the 400 µm upper bound.
* Among the 12 interior fits, d̂½ has median 29 µm but ranges 3–123 µm.

The reason is concrete and measurable: **at 5 % sender prevalence, 99 % of cells lie
within 72–90 µm of a sender** (median 23–31 µm; fewer than 0.04 % of cells are beyond
150 µm). The plan's prescription of "10 µm bins out to 300 µm" is unreachable at this
prevalence — there is essentially no data past ~130 µm. Fitting an exponential with
λ up to 400 µm is extrapolating roughly five-fold beyond the observed range, and the
profile likelihood is correspondingly flat: for the best-behaved fit the deviance
changes by less than 12 units over λ ∈ [78, 200] µm.

**Sender prevalence sets the observable distance range, and the observable range caps
the resolvable λ.** Any long-range SASP length constant reported at a few-percent
sender prevalence is extrapolation, not measurement.

### 3. The kernel family you choose changes the answer by 5×

For the *same* data, median d̂½ by family: **step 26 µm, gaussian 44 µm, exponential
71 µm, powerlaw 75 µm, spline 136 µm.** Per section × module, the ratio of the largest
to the smallest d̂½ across the four parametric families has **median 5.0× (IQR
3.9–6.7)**. No published spatial-CCC result I am aware of reports this sensitivity.

**The nonparametric spline wins AIC in 34 of 35 fits.** Taken at face value that says
no parametric family is adequate. Phase 1b calibrates that claim and shows it means
something different — see §6.

### 4. The standard iid CI is 4–5× too narrow, and far worse for some families

Restricting to interior optima (a bound-railed λ makes the iid SE meaningless):
median SE understatement factor (block-bootstrap sd ÷ iid asymptotic SE) is
**exponential 5.0× (max 14.7), powerlaw 4.8×, gaussian 3.9×, step 280× (max 8149×)**.

This **confirms the Phase 1 synthetic prediction on real tissue** (Phase 1 reached
7.9×). The iid CI is not merely optimistic — for several sections it is nonsensical:
e.g. 7259/tnfa_nfkb_proximal gives an iid 95 % CI on λ of **[−220, 1020] µm**, a
negative length scale.

The honest block-bootstrap CIs are wide: median width 10 µm around a median d̂½ of
29 µm for interior fits, and **8 of 35 fits have a 95 % CI spanning more than a
10-fold range in d̂½.**

### 5. Most of the gradient is explained by transcript count, cell area and local density

The module scores correlate with log total transcript count at **r = 0.22–0.67**
(downstream_arrest 0.67, oxidative_stress 0.60, tnfa_nfkb_proximal 0.40), and the
Tier A sender score itself correlates with counts at r = 0.29. Distance to nearest
sender correlates with local density at r = −0.23.

Adding log counts, log cell area, density at 25/50/100 µm and the 1-NN distance to
the design — a preview of null N5 — **removes 86 % of β̂ (median |β̂| 0.266 → 0.056,
ratio 0.14)** and collapses the median d̂½ from **71 µm to 7 µm**, i.e. below the
resolution floor, with 66 % of adjusted fits below 10 µm. Per module the surviving
fraction ranges from 0.60 (emt_ecm) to 0.12 (il6_jak_stat3).

This is the Phase 1 density-confound mechanism appearing in real tissue, and it is
*before* any matched-decoy or torus-shift null is applied. The full null battery is
Phase 3, but the direction is already clear.

### 6. Length scales are not conserved across sections, even after packing normalisation

In units of the local median NN distance, median d̂½ per section is **3.2, 4.6, 8.8,
13.6 and 41.2 NN units** — a 13-fold spread across five mice. Packing normalisation
does not rescue comparability; if anything the sections disagree more in NN units than
in microns. With one section per mouse and five mice, this is a case study under
Section 24.1, and no cross-condition (sbr vs sham) contrast is attempted here: the
sections differ in median NN distance by 45 % and in cell count by 3×, which Phase 1
showed is sufficient on its own to manufacture a difference in λ̂.

---

## Part B — The Phase 1 work still owed

All three items are on synthetic tissue with a known planted answer, 15 replicates per
cell (12 for the nulls), in two regimes: **clean** (random senders, no latent
confounder, ℓ/λ = 0.25) and **confounded** (κ = 3, conf_strength = 1, ℓ/λ = 2).

### 7. Kernel-family misspecification — previously untested

Each of the four parametric families was planted in turn and all five fitted.

**Can AIC recover the true family?** Only sometimes, and the failure mode is
systematic:

| planted | picked exponential | gaussian | powerlaw | step | spline |
|---|---|---|---|---|---|
| exponential | **9** | 4 | 2 | 0 | 15 |
| gaussian | 0 | **13** | 0 | 0 | 17 |
| powerlaw | 13 | 1 | **2** | 0 | 14 |
| step | 0 | 0 | 0 | **19** | 11 |

(30 fits per row.) Power law and exponential are **badly confusable** — when a power
law is planted, AIC picks the exponential 13/30 times and the power law only 2/30.
Step is well discriminated. Gaussian is never mistaken for another parametric family.

**How wrong is the length scale when the wrong family is fitted?** Ratio of fitted to
true d½, median over replicates, excluding the spline:

* power law planted, exponential fitted: **2.05×** (clean) — you would report a
  length scale twice the truth;
* step planted, exponential fitted: **0.53×**;
* gaussian planted, powerlaw fitted: **0.49×** (confounded);
* the full range across all planted × fitted combinations is **0.36× to 2.67×**.

So family misspecification alone is worth a factor of ~2–3 in the reported length
scale — the same order as the 5× family spread measured on the real data (§3).

### 8. The spline "winning" is a marker of confounding, not of a non-parametric kernel

This is the finding that reinterprets the real-data result. The rate at which the
nonparametric spline wins AIC **when a parametric family is in fact the truth**:

| planted | clean regime | confounded regime |
|---|---|---|
| exponential | 0.07 | **0.93** |
| gaussian | 0.20 | **0.93** |
| powerlaw | 0.07 | **0.87** |
| step | 0.33 | 0.40 |

In the clean regime the spline rarely wins (0.07–0.33); in the confounded regime it
almost always does (0.87–0.93). **The spline wins because spatial confounding adds
structure no parametric decay kernel can absorb, not because the underlying kernel is
non-parametric.**

On the real data the spline wins **34 of 35** fits. Read against this calibration,
that is not evidence that the SASP kernel is exotic — **it is evidence that the real
sections sit in the confounded regime of the Figure 1 regime map.** This is an
independent line of evidence pointing the same way as §5.

### 9. Superposition vs nearest-sender IS identifiable (Section 6.3) — a positive result

Planting each and fitting both, model comparison is **perfect: 15/15 correct in every
regime**, with no errors in either direction.

| planted | regime | median AIC(superposition) − AIC(nearest) |
|---|---|---|
| nearest-sender | clean | +93 |
| nearest-sender | confounded | +177 |
| superposition | clean | **−664** |
| superposition | confounded | **−3316** |

**Unlike λ, the dose-versus-threshold question is answerable from a static snapshot,
and confounding does not degrade it** — it makes the discrimination easier, because
the superposition regressor is more distinguishable from the nuisance structure.
Given how little else in this project is identifiable, this deserves to be promoted:
it is a question the plan says "neither has been reported" for, and it is one we can
actually answer.

### 10. Nulls N1 / N3 / N4 — the torus shift does not do what the field thinks

Measured properly: **size** (rejection rate when β_true = 0, nominal 0.05) and
**power** (rejection rate when β_true = 1).

First, an implementation trap. Re-profiling λ inside the null — the obvious
implementation — compares two *different models*: when a null destroys the spatial
structure, λ̂ rails to the top of the grid, the kernel column becomes nearly constant,
and β̂ stops being an amplitude. That yields nonsense (N1 power = 0.00 in the
confounded regime, and a *negative* surviving fraction of −0.72). **Holding λ fixed at
the observed λ̂ fixes the power** — it becomes 1.00 everywhere. Both variants are
retained in `nulls_reprofiled_only.csv` and `nulls.csv`; the numbers below are the
correct fixed-λ version.

**Size (β_true = 0, nominal 0.05):**

| regime | N1 label permutation | N3 torus shift | N4 rotation |
|---|---|---|---|
| clean | 0.42 | 0.67 | 0.67 |
| confounded | 0.50 | **0.92** | **1.00** |

**Power (β_true = 1): 1.00 for all three nulls in both regimes.**

Every one of these nulls has badly inflated Type I error. But the p-value is the
wrong thing to look at, because the permutation null distributions are narrow, so a
negligible effect still clears the 95th percentile. The quantity that matters is the
one Section 6.5 already prescribes — the **surviving fraction**, (β̂_obs − null mean)
/ β̂_obs:

| regime | truth | N1 | N3 | N4 |
|---|---|---|---|---|
| clean | β = 1 (should be ≈1) | 0.86 | 0.96 | 0.91 |
| confounded | β = 1 (should be ≈1) | 0.69 | 0.99 | 0.96 |
| clean | **β = 0 (should be ≈0)** | 0.31 | **0.74** | **0.58** |
| confounded | **β = 0 (should be ≈0)** | **0.006** | **0.98** | **0.89** |

**This is the most important result in Part B.** With *no planted kernel whatsoever*,
in the confounded regime, the torus-shift null reports that **98 % of the observed
effect "survives"** — it certifies pure confounding as a genuine SASP effect. Rotation
does the same at 89 %.

The reason is structural, not a bug. N3 and N4 preserve sender clustering and receiver
autocorrelation and destroy only their **alignment**. But a shared confounder produces
alignment without any signalling. **The torus shift therefore tests "are the senders
aligned with the response field?", which is not the same question as "is there a SASP
effect?"** Passing the torus-shift null is not evidence of causation when senders and
response share a spatial cause.

**N1, cell-type-stratified label permutation, is the only null of the three that is
directionally trustworthy under confounding**: surviving fraction 0.006 when there is
no effect and 0.69 when there is. It is conservative — it discards about a third of a
real effect — but it does not manufacture one. Its behaviour is *worse* in the clean
regime (0.31 when β = 0), which is a further instance of the Phase 1 conclusion that
every control is regime-dependent.

---

## Recommendations carried into Phase 3

1. **Report the surviving fraction, never the permutation p-value.** All three nulls
   reject at 42–100 % under a true null; their p-values carry no information. The
   Section 6.5 prescription is empirically vindicated.
2. **Demote N3/N4 from "the strong null" to "a necessary but weak check."** They
   cannot separate confounding from causation, and on our synthetic data they endorse
   confounding at 89–98 %. If the paper reproduces the CellWHISPER finding that
   existing tools fail the torus shift, it should simultaneously report that *passing*
   the torus shift proves much less than the field assumes.
3. **Hold λ fixed at λ̂_obs inside every permutation null.** Re-profiling silently
   changes the model and destroys power.
4. **Report the family-spread on d̂½ as part of the estimate.** On the real data it is
   a factor of 5; misspecification alone contributes a factor of 2–3.
5. **Report the observable distance range next to every λ̂.** At 5 % prevalence, 99 %
   of cells are within 90 µm of a sender, so any λ̂ above ~40 µm is extrapolation.
6. **Promote the superposition-vs-nearest comparison.** It is the one estimand in this
   project that is cleanly identifiable, and it answers a question nobody has reported.

## What I did NOT get to

* **Receiver cell-type stratification.** The Bio agent's `celltypes_*.csv` had not
  landed, so every real-data fit is on a single unstratified stratum. Figure 2a is
  therefore per section × module, **not** per receiver cell type as Section 25
  specifies. The code path exists and activates automatically on re-run; this is the
  first thing to redo when the annotations arrive.
* **Bio's sender calls, module scores and anatomy.** All real-data numbers use
  provisional calls (`tierA_p95`, z-mean Tier B module scores) and must be re-run.
  Everything is a single `prepare_samples.py --force` plus a re-run of the fitter.
* **Held-out log-likelihood on left-out sections** (Section 24.6). AIC is reported;
  the leave-one-section-out cross-validation is implemented in `sasp_kernels.
  gaussian_heldout_ll` but not wired into the runner, because with one section per
  mouse and provisional per-section z-scored modules the pooled fit is not yet
  well-posed. Needs Bio's harmonised module scores first.
* **The remaining 6 of 11 liver sections**, which were still downloading.
* **Nulls N2, N5, N6, N7, N8 on real data.** §5 is only a *preview* of N5. The full
  battery is Phase 3.
* **Anatomical confounding (Section 8 Test 6 / Section 11 zonation).** The loader
  reads `anatomy_*.csv` when present; nothing has been fitted against it.
* **Misspecification under superposition truth** was run but is not reported above,
  and the spline's d½ summary is unreliable in the clean regime (it inflates 7–19×
  when the true kernel decays inside the first knot span); the spline should be
  compared by AIC, not by d½.
