# CS Phase 5 — Superposition, winner's curse, kernel families, and the corrected Figure 2a

**Status: complete, and T1 is now closed at three sender definitions (§10).
Verdict: the one estimand that was still plausibly positive is not identifiable
on this tissue either, and it now has a bound that holds at every caller. The
Phase 3 headline is not a selection artefact — the selection bias runs the other
way. NEW in §10: whether the superposition regressor earns its place at all is
sender-definition dependent — it does not under `tierA_p95` and does under
`cdkn1a_pos`/`senepy_p95` — so two of the T1 sentences below are re-scoped to
the caller they were measured at.**

Master Plan Sections 6.2, 6.3, 6.4, 24.6, 25. Everything below is fitted on the
six Test-3-admissible sections (six animals, both arms), receiver labels
`cell_type_merged`, `Low_quality`/`Unknown` excluded, sender call `tierA_p95`
— except **§10 (T1b)**, which re-runs T1 at `cdkn1a_pos` and `senepy_p95`.

| Deliverable | Location |
|---|---|
| Superposition regressor, held-out LL machinery | `/workspace/code/phase5_common.py` |
| T1 superposition vs nearest (3 stages) | `/workspace/code/run_phase5_super.py` |
| **T1b the same, at every sender call** | `run_phase5_super.py --call C --tag _C`, collated by `/workspace/code/summarize_super_callers.py` |
| T2 winner's curse (cross-fit + synthetic) | `/workspace/code/run_phase5_wc.py` |
| T4 kernel families + T5 proximal/downstream | `/workspace/code/run_phase5_kernels.py` |
| Why "the spline wins 34/35" does not reproduce | `/workspace/code/_spline_window_check.py` |
| Figures | `/workspace/code/make_phase5_figs.py` |
| Every number below, printed from the CSVs | `/workspace/results/phase5/summary_phase5.txt` |
| **Figure 2a (regenerated), Figure 3 (new)** | `/workspace/figures/figure2a.{png,pdf}`, `figure3.{png,pdf}` |
| iid vs block bootstrap SEs for the new estimands | `/workspace/code/_se_ratio_phase5.py` |
| Result tables (10 CSVs) | `/workspace/results/phase5/` |

No fitter was rewritten. `run_phase3_nulls.SectionFit` still owns the receiver
set, the 100 µm window, the [7, 50] µm λ grid, the N5/N6 covariate blocks and the
N2 matched decoys; `sasp_kernels.BasisBlockProfiler` still does every solve.

---

## 0. Headline

1. **Superposition vs nearest-sender is not identifiable under control.** Under
   N5 + N6 the median ΔAIC is **−0.060 per 1,000 cells** against a synthetic
   signature of **−57 to −251** when superposition is planted and **+15 to +17**
   when nearest-sender is planted — the same comparison that got **15/15** right
   in every synthetic regime. **§10 re-runs this at `cdkn1a_pos` and
   `senepy_p95`, which agree with Tier A at or below chance, and the bound holds
   at all three (−0.06, −0.15, −0.32 per 1,000 cells).** Held-out log-likelihood
   on left-out sections, the Section 24.6 criterion, is a coin flip *at
   `tierA_p95`*: superposition wins **50.8 %** of 252 folds — but **0.679 at
   `cdkn1a_pos` and 0.603 at `senepy_p95`** (§10.3). **The nearest-sender
   regressor never beats the no-kernel model under any caller**; the
   superposition regressor does not at `tierA_p95` (35 % of fits) and **does at
   the other two (61–64 %)**, with its λ̂ railed at the top of the grid where it
   helps — a regional density covariate, not a contact-scale kernel.
2. **At `tierA_p95` the verdict is not stable under the nulls.**
   Superposition "wins" 72.1 % of fits on real data, **66.6 % on N1-permuted
   senders and 55.5 % on N3 torus-shifted senders**. Most of the preference is
   reproduced by nulls that contain no signalling. What little the model
   selection is measuring is tissue geometry — local neighbour counting in a
   tissue of non-uniform density — not dose-versus-threshold SASP.
   **This sentence is caller-specific and §10 corrects it**: the null win
   fractions are near-identical at all three callers (N1 0.64–0.67, N3
   0.53–0.56), but the *observed* ΔAIC/1k is 2.7×/5.5× the N1/N3 null at
   `tierA_p95` and **17×/65× at `cdkn1a_pos`**. The win fraction alone hides
   that; the per-1k magnitude does not.
3. **The winner's curse is real, is worth +0.056 of surviving fraction, and has
   the opposite sign to the one the objection assumes.** Cross-fitting shows the
   in-sample SF(N5) is *inflated* by selection, not deflated. Correcting for it
   makes the Phase 3 negative **stronger**, not weaker. No correction to the
   headline is required; this becomes a robustness paragraph.
4. **Figure 2a is regenerated and the composition artefact is now visible on the
   figure.** The unstratified contact amplitude is +0.260 response-sd; replacing
   every cell's response by its cell type's section mean reproduces **+0.212 sd,
   i.e. 76 % of it**, with no signalling of any kind.
   > **⚠ THE "76 %" IS UNSOURCED AND THE ROW IS SUPERSEDED — 2026-08-27 (record
   > reconciliation).** 0.212 / 0.260 = **0.815**, not 0.76. The only two "ALL RECEIVERS"
   > amplitudes on disk are **0.260** (`results/phase5/summary_phase5.txt` T3, identical in
   > `results/phase5_pre_c6/`) and **0.266001** (`figures/figure2a_amplitudes.csv`); neither
   > yields 76 %, and a denominator of ≈ 0.279 would, which is in no file. **The 76 % and the
   > §17 "composition surrogate share 66–76 %" row are both withdrawn.** They are replaced by
   > the composition-matched decomposition on the same fits: **65.9 %** removed by the
   > receiver's own cell-type intercepts and **85.4 %** by cell type + the 20-NN composition
   > vector, against **1.6 %** for the same variables used as matched decoys
   > (`results/phase3/compmatch_reruns.csv`, `row_type == 'summary'`). Wherever the 1.6 %
   > appears the other two must appear beside it (`PREREG_PHASE8.md` §10.8). See
   > `reports/WRITING_PACK.md` §0.2 and `reports/CS_PHASE8_COMPMATCH.md`. Within receiver cell type
   the pooled amplitude falls to +0.091 sd, and it ranges from **+0.24
   (hepatocytes) to −0.06 (biliary/ductular)**.
5. **Under control no kernel family earns its place**, the family that comes
   closest is the **step/threshold** function (wins 90.8 % of fits by AIC and
   still beats the no-kernel model in only 55.9 %), family choice still moves
   d̂½ by a median **4.4×**, and the Phase 2 "spline wins 34/35" marker turns out
   to be **the receiver-cell-type composition artefact again**, not a window
   artefact: unstratified the spline wins 95 %, stratified it wins 21 %.
6. **λ_proximal vs λ_downstream is not estimable**, as expected. The donor
   bootstrap CI on λ_B1/λ_B4 **reaches both ends of the widest ratio the grid
   allows, [0.14, 7.14], in 3 of 6 receiver cell types**, and covers 1 in all
   six. Reported as a bound, not a point estimate.

---

## 1. A bug worth recording before any of the numbers

The first version of T1 drew the spatial block bootstrap **inside** each call to
the basis profiler, so the nearest-sender model and the superposition model were
evaluated on *different* resamples of blocks. The paired win fraction then
collapsed to exactly 0.500 by construction and the per-fit CI on ΔAIC came out at
±57 per 1,000 cells — wider than the planted-superposition signal itself. The
fix is one line of plumbing (`profile_basis` now takes the multiplicity matrix in
rather than drawing it), and it moves the paired win fraction from 0.500 to 0.700
and the CI from ±57 to [−0.86, +0.24].

This is worth a Methods sentence: **any paired model comparison under a block
bootstrap has to share the resample**, and the failure mode is silent — it
produces exactly the "we cannot tell them apart" answer that a negative result
would also produce.

---

## 2. T1 — Superposition vs nearest-sender (Section 6.3)

### 2.1 What was fitted

For each of the 315 (section × receiver cell type × module) cells:

```
nearest        r_i = mu_c + beta * exp(-d_i / lam)                + gamma'z_i
superposition  r_i = mu_c + beta * sum_{j in S} exp(-d_ij / lam)  + gamma'z_i
```

Both profiled over the identical [7, 50] µm grid, on the identical receivers
(the Phase 3 set: non-sender, labelled, within 100 µm of a sender), under the
identical design, with the identical number of free parameters. ΔAIC therefore
reduces exactly to `n · log(RSS_sup / RSS_near)` and cannot be a parameter-count
artefact. The superposition sum is built with a chunked `cKDTree`
sparse-distance-matrix truncated at 6 λ_max = 300 µm; no (n, n) matrix is formed
and the whole basis for a 238k-cell section costs ~9 s.

Three designs are reported: **naive** (intercept only, within cell type),
**ctrl** (N5 nuisance block + N6 receiver baseline), and **ctrl + N2** (each
model additionally carrying its own matched-decoy regressor at a shared λ).

### 2.2 The calibration, restated in a comparable unit

Phase 2's ΔAIC values (−664 to −3316) are not directly comparable to real fits
because the synthetic sections have ~11k cells and the real ones 3k–90k per fit.
Per 1,000 cells (`results/phase1b/misspec.csv`):

| planted | regime | median ΔAIC per 1,000 cells | verdict correct |
|---|---|---|---|
| nearest-sender | clean | **+15.2** | 15/15 |
| nearest-sender | confounded | **+17.0** | 15/15 |
| superposition | clean | **−57.2** | 15/15 |
| superposition | confounded | **−250.9** | 15/15 |

### 2.3 The real result

| design | median ΔAIC / 1k cells | IQR | superposition wins | paired bootstrap win fraction (median) | fits where the bootstrap CI excludes 0 |
|---|---|---|---|---|---|
| naive | **−0.448** | [−2.11, −0.03] | 0.784 | 0.840 | 0.216 |
| **ctrl (N5+N6)** | **−0.060** | [−0.238, +0.005] | 0.721 | 0.700 | **0.067** |
| ctrl + N2 | −0.098 | [−0.364, +0.003] | 0.733 | — | — |

**The bound.** Under full control the per-fit 95 % block-bootstrap interval on
ΔAIC per 1,000 cells is **[−0.86, +0.24]** in the median fit and reaches
**|ΔAIC/1k| ≤ 2.81 at the 90th percentile**. We can therefore exclude a
planted-superposition signature (−57 to −251) by a factor of **20–90 even at the
p90 fit**, and a planted-nearest signature (+15 to +17) by a factor of ~6.
Superposition is decisively preferred by the bootstrap (win fraction ≥ 0.95) in
**7.6 %** of fits; nearest-sender is decisively preferred in **0.0 %**.

In amplitude units the same statement is: the controlled superposition amplitude
is **−0.005 response-sd** (IQR −0.016 to +0.014, p90 +0.026) against a naive
+0.029, and the controlled nearest-sender amplitude is −0.012 sd against a naive
+0.193 sd. Both are an order of magnitude below the 0.203 sd that Phase 3
established this design can detect at 80 % power.

### 2.4 Neither model earns its place at all

This is the more damaging framing, and it is what makes the model selection moot:

| design | nearest vs covariates-only | superposition vs covariates-only |
|---|---|---|
| naive | ΔAIC **−10.1**, improves 78 % of fits | ΔAIC **−19.9**, improves 85 % |
| **ctrl** | ΔAIC **+2.75**, improves **17 %** | ΔAIC **+1.59**, improves **35 %** |
| N1 permuted senders | +3.05, 9 % | +2.47, 18 % |
| N3 shifted senders | +2.30, 28 % | +1.94, 30 % |

Under control, adding *either* kernel term makes AIC **worse** in the majority of
fits. Asking which of two kernels fits better, when neither improves on having no
kernel, is a question about the shape of noise.

### 2.5 Held-out log-likelihood on left-out sections (Section 24.6)

Pooled per receiver cell type across the six sections with responses z-scored
within section and both bases z-scored within section (mandatory — the
superposition sum scales with the section's sender density, so an unstandardised
amplitude is not transferable and the held-out comparison would be a density
comparison), leave-one-section-out, 252 folds:

| design | ΔLL per cell (sup − near) | superposition wins | nearest beats no-kernel | superposition beats no-kernel |
|---|---|---|---|---|
| naive | +2.3 × 10⁻⁵ | 0.544 | 0.754 | 0.698 |
| **ctrl** | **+8.1 × 10⁻⁷** | **0.508** | **0.532** | **0.548** |

Under control this is a coin flip on both questions. **This is the Section 24.6
criterion the plan actually specifies, and it is the cleanest single line in the
T1 result.**

### 2.6 Is the verdict stable under the nulls? No — and that is the answer

Five draws each of N3 (torus shift of the sender point set) and N1
(cell-type-stratified sender-label permutation), controlled design, 1,575
null fits each:

| data | median ΔAIC / 1k | superposition wins | per-draw range |
|---|---|---|---|
| **observed** | −0.060 | **0.721** | — |
| N1 permuted senders | −0.022 | **0.666** | 0.632 – 0.708 |
| N3 torus-shifted senders | −0.011 | **0.555** | 0.537 – 0.581 |

The observed superposition preference (0.721) is **barely outside the N1 range
and only 0.17 above the N3 value**, on data where the senders have been replaced
by a random subset of the same cell types. Two readings, both worth stating:

* **N1 > N3 is exactly what a geometric explanation predicts.** N1 keeps the
  sender point set inside the real tissue (senders are real cells, so the sum
  Σ exp(−d/λ) still tracks real local density); N3 moves the point set off the
  tissue's own density field. If the preference were about signalling, N1 —
  which destroys sender identity completely — should have killed it.
* **The superposition regressor is, under this control, largely a smoothed local
  density.** The N5 block already contains density at 25/50/100 µm and the 1-NN
  distance, which is why the preference falls from −0.448 to −0.060 per 1k cells
  when the control is applied. The correlation between the two regressors at
  their respective λ̂ is **0.656** [IQR 0.53, 0.86] — they are not even two very
  different descriptions.

**Conclusion for Section 6.3** (read with §10, which re-runs everything in this
subsection at `cdkn1a_pos` and `senepy_p95`; the bound holds at all three, the
null-instability reading does not). The dose-versus-threshold question is answerable
from a static snapshot *in principle* — Phase 1b demonstrated 15/15 — and it is
**not answerable on this tissue at this effect size**. The correct sentence for
the paper is not "superposition wins" but: *the one estimand our synthetic study
identified as cleanly identifiable is, on real tissue under control, separated
from its alternative by 0.06 AIC units per thousand cells, three orders of
magnitude below the planted signature, with a torus-shift null reproducing most
of the residual preference.* Reporting the 72 % win rate without the null column
would have been a manufactured positive.

### 2.7 iid vs spatial block intervals for the new estimands

Every interval above is a spatial block bootstrap. `_se_ratio_phase5.py` measures
what an iid asymptotic interval would have claimed instead, for the amplitude
β̂ at the profiled λ̂, over the same 315 fits (`results/phase5/se_ratio.csv`):

| design | regressor | SE ratio block/iid, median [IQR] | p90 | CI excludes 0: iid | block |
|---|---|---|---|---|---|
| naive | nearest | 1.20 [1.10, 1.38] | 1.57 | 0.781 | 0.737 |
| naive | superposition | **1.40** [1.20, 1.83] | **2.78** | 0.851 | 0.740 |
| ctrl | nearest | 1.04 [0.97, 1.12] | 1.20 | 0.175 | 0.143 |
| ctrl | superposition | 1.09 [1.00, 1.21] | 1.40 | 0.368 | 0.308 |

Two qualifications, because these numbers are smaller than the 4–8× quoted
elsewhere in the project and the difference is not a contradiction:

* **The 4–8× figure is for λ̂, not for β̂.** Phase 1's `se_ratio_naive` is
  `boot_sd(λ̂) / se_iid(λ̂)`; the table above is the amplitude at a held λ̂, which
  is a linear parameter and much better behaved.
* **This block bootstrap is within-section** (10 × 10 quantile blocks), so it
  cannot see between-animal variance. That is what the donor bootstrap in
  Figure 3a and in T5 is for, and it is where the intervals actually blow up —
  39 of 42 donor CIs on λ̂ span the entire admissible grid.

The interesting part is the **ctrl** rows: once the N5 block is in the model the
understatement factor falls to ~1.05, i.e. **most of the spatial dependence in
the residual was the nuisance structure itself**. That is a constructive
observation to carry to the field — conditioning on depth, size, density and
neighbourhood composition does not only remove 92 % of the amplitude, it also
removes most of the reason the naive standard error was wrong.

---

## 3. T2 — How much of SF(N5) = 0.084 is winner's curse?

### 3.1 The mechanism, and why the sign is not obvious

Phase 3 selects the 160 of 315 fits with `beta_naive > 0` and a block-bootstrap
CI excluding zero, then reports SF = β_controlled / β_naive **on the same data**.
The objection is that selection inflates β̂_naive, which sits in the
*denominator*, so SF is deflated by regression to the mean with no control doing
any work.

The objection is incomplete, because β̂_N5 and β̂_naive are nested fits on the
same cells and share most of their sampling noise. Selecting on a positive noise
draw in β̂_naive selects the same positive noise draw in β̂_N5, which pushes the
*ratio* toward 1. Whether selection deflates or inflates SF is therefore an
empirical question, and it has to be measured.

### 3.2 Cross-fitting on the real fits

Each fit's 100 spatial blocks are split at random into A and B; the Phase 3
selection rule (including its 200-replicate block-bootstrap CI) is applied
**within A**; SF is then estimated on **B** with B's own λ̂ and B's own β̂_naive,
which are statistically independent of the selection. The like-for-like
in-sample comparator is SF on **A**, at the same sample size. 20 random splits,
seeds pinned. A **placebo** row repeats the same in/out comparison over all 315
fits with no selection at all, which isolates the selection effect from any
in-sample/out-of-sample gap the estimator has on its own.

| selection fraction | SF(N5), in-sample | SF(N5), held-out | gap |
|---|---|---|---|
| 0.50 — selected on A (rate 0.421) | **0.1375** (sd 0.016) | **0.0794** (sd 0.025) | **+0.058** |
| 0.50 — **placebo**, all 315 fits | 0.1552 | 0.1534 | **+0.002** |
| 0.75 — selected on A (rate 0.481) | 0.0890 | 0.1624 | −0.073 |
| 0.75 — **placebo**, all 315 fits | 0.1152 | 0.2448 | −0.130 |

**Placebo-corrected winner's curse = +0.0563 at both split fractions** (0.058 −
0.002 and −0.073 − (−0.130)). The agreement of two splits with different
selection and estimation sample sizes to three decimal places is reassuring.
The same exercise on SF(N2+N5+N6) gives +0.055 and +0.057.

Without selection there is **no** in-sample/out-of-sample gap (+0.002), so the
estimator itself is unbiased in this respect; with selection there is a gap of
+0.056. That gap is the winner's curse, and **its sign is upward**: the reported
in-sample SF is too generous to the "the effect survives" reading.

The 0.75 rows also expose a fact worth its own sentence, because it is a trap for
anyone who reads surviving fractions across studies with different sample sizes:

| sample | median SF(N5), no selection |
|---|---|
| 100 blocks (full) | 0.096 |
| 50 blocks | 0.157 |
| 25 blocks | 0.247 |

**SF is pulled toward 1 as n falls**, because the shared noise in numerator and
denominator grows relative to the signal. A surviving fraction is not comparable
across sample sizes, and a small study will report a *higher* surviving fraction
than a large one for the same underlying truth.

### 3.3 Synthetic replication of the selection rule

Applying the identical rule to the Phase 1 sweep, where the truth is known
(`results/sweep_all.csv`):

| truth | selection rate | SF(N5), all runs | SF(N5), selected |
|---|---|---|---|
| planted effect (β = 1, 1,150 runs) | 1.000 | 0.8252 | **0.8252** |
| no effect (β = 0, 120 runs) | 0.867 | 0.4123 | **0.4454** |

With a planted effect every run passes the screen, so selection cannot move
anything. With no effect, selection moves SF **up** by +0.033 — the same
direction as the cross-fit, and about half the size.

### 3.4 What this means for the Phase 3 headline

**No correction is required.** Three statements, in the order a reviewer will
want them:

1. The measured selection contamination is **+0.056** of surviving fraction, at
   half-sample noise, and it points the wrong way for the objection: the honest
   cross-fitted SF is *lower* than the reported one.
2. The gap that would have to be explained away is **0.826 − 0.084 = 0.742** —
   the distance from the observed SF to the planted-effect calibration anchor.
   A +0.056 artefact of the wrong sign does not touch it.
3. Phase 3's own argument survives intact and is now quantified rather than
   asserted: N3 returns 1.000 and N2 0.943 on the identical selected set, so
   selection cannot be what drives N5 to 0.084.

This is a robustness paragraph, not an erratum. It should be written as one,
because "you selected on significance" is the first thing a reviewer will say.

---

## 4. T3 — Figure 2a, regenerated and stratified

`figures/figure2a.{png,pdf}`, data in `figures/figure2a_stratified_curves.csv`,
amplitudes in `figures/figure2a_amplitudes.csv`.

Ten panels: the unstratified curve the field plots, plus one panel per receiver
cell type, all seven Tier B modules on each, six admissible sections pooled,
non-uniform distance bins (the distance distribution is concentrated below 30 µm,
so uniform 5 µm bins spend three quarters of the axis on an empty tail), √-scaled
x-axis, and the < 10 µm region shaded because that is where the whole apparent
effect lives.

**The unstratified panel carries a dotted line that is new and is the point of
the figure**: each cell's response replaced by the mean response of its cell type
in that section, then binned identically. That curve contains no signalling, no
distance dependence, and no kernel — only composition.

Contact amplitude (first retained bin minus the 40–100 µm plateau, within-stratum
sd units):

| panel | amplitude |
|---|---|
| **unstratified** | **+0.260** |
| **composition-only surrogate** | **+0.212** — ~~76 %~~ **81.5 %** of the unstratified curve (0.212 / 0.260; the 76 % is unsourced — see the note in §4), per-module ratios 0.42–1.83 |
| pooled, within cell type | **+0.091** — a 65 % reduction, matching the 66 % Phase 3 measured on β̂ |
| Hepatocytes | +0.236 |
| Endothelial | +0.211 |
| Macrophages | +0.206 |
| DC | +0.186 |
| Proliferating | +0.185 |
| Mesenchymal | +0.081 |
| vSMCs | +0.056 |
| T/NK cells | **−0.039** |
| Biliary/ductular | **−0.059** |

Two things a reader should take from the panels. First, the sign is not even
consistent: two of nine receiver cell types show the response *rising* with
distance. Second, the annotation on each panel — **only 2.4 % of hepatocyte
receivers are within 10 µm of a sender**, 3.9 % of endothelial, 6.8 % of
macrophages — is the Phase 3 "the contact spike rests on a few dozen cells"
finding made visual across every cell type at once.

The old Figure 2a (Phase 2, unstratified, provisional module scores, 300 µm axis)
must not be used. It is superseded by this file at the same path.

---

## 5. T4 — Kernel families under full control (Section 6.2)

Five families — exponential, Gaussian, power law (p ∈ {0.5, 1, 2, 4}),
step/threshold, natural cubic spline (6 knots, dmax = window) — every one a
precomputed basis through the same `BasisBlockProfiler`, same cells, same design,
same λ grid, 315 fits stratified by receiver cell type.

### 5.1 AIC

| family | AIC win, naive | AIC win, **ctrl** | median ΔAIC vs no kernel, ctrl | beats no kernel, ctrl | median d̂½ (ctrl) | λ̂ railed (ctrl) |
|---|---|---|---|---|---|---|
| exponential | 0.143 | 0.006 | +2.75 | 0.165 | 10.0 µm | 0.689 |
| gaussian | 0.092 | 0.029 | +2.20 | 0.203 | 18.6 µm | 0.502 |
| power law | 0.063 | 0.000 | +4.35 | 0.098 | 4.6 µm | 0.711 |
| **step** | 0.489 | **0.908** | **−0.57** | **0.559** | 15.8 µm | 0.114 |
| spline | 0.213 | 0.057 | +6.64 | 0.137 | — | — |

**Under control the step/threshold function wins 91 % of fits, and it is the only
family whose median AIC beats having no kernel at all — by 0.57 units, i.e. it
beats no kernel in 55.9 % of fits, a coin flip.** The other four families make
AIC *worse* than the no-kernel model in 80–90 % of fits. The honest reading is
not "the response is thresholded"; it is that under control the only kernel shape
that is not actively harmful is the one that says the response is *not graded*,
and even that is not distinguishable from nothing.

**Family choice still moves the reported length scale.** Median d̂½ max/min across
families within a fit: **4.36× under control**, 3.20× naive (Phase 2 measured a
median 5× on the Phase 2 configuration). Conditioning does not reduce the spread
— if anything it widens it. Anyone reporting a SASP length constant without
reporting which family they fitted is reporting a number with a 4× discretion.

### 5.2 Held-out log-likelihood on left-out sections (Section 24.6)

252 leave-one-section-out folds, per family, per design:

| family | ctrl: beats no kernel | ctrl: held-out win share | naive: beats no kernel | naive: win share |
|---|---|---|---|---|
| exponential | 0.532 | 0.143 | 0.754 | 0.083 |
| gaussian | 0.536 | 0.163 | 0.762 | 0.111 |
| power law | 0.532 | 0.139 | 0.770 | 0.317 |
| step | 0.452 | 0.238 | 0.738 | 0.194 |
| spline | 0.429 | 0.317 | 0.702 | 0.294 |

Under control, **no family beats the covariates-only model in more than 53.6 % of
folds**, and step and spline are actively *worse* than no kernel (0.45, 0.43).
The held-out winner is essentially uniform over the five families (0.14 – 0.32):
model selection has no signal to work with. Naively every family beats no kernel
in 70–77 % of folds and the winner is again split between power law and spline —
i.e. **AIC and held-out log-likelihood disagree about which family wins**, which
is itself a reason not to report a family.

### 5.3 "The spline wins 34/35" is a composition artefact, not a shape finding

CS Phase 2 §8 reported the spline winning 34 of 35 real fits and read it, against
the synthetic calibration (spline wins 0.87–0.93 confounded, 0.07–0.33 clean), as
evidence that the sections sit in the confounded regime. Phase 5 stratifies by
receiver cell type and the spline wins only 21 %. That discrepancy needed
explaining, so I ran the 2 × 2 (`_spline_window_check.py`, naive design, same six
sections):

| configuration | spline wins | step wins |
|---|---|---|
| Phase 2 window 300 µm, λ ∈ [3, 400], **unstratified** (42 fits) | **0.881** | 0.024 |
| Phase 5 window 100 µm, λ ∈ [7, 50], **unstratified** (42 fits) | **0.952** | 0.024 |
| Phase 2 window 300 µm, λ ∈ [3, 400], **stratified** (322 fits) | 0.193 | 0.497 |
| Phase 5 window 100 µm, λ ∈ [7, 50], **stratified** (315 fits) | 0.213 | 0.489 |

**It is entirely a stratification effect and not at all a window effect.** The
spline was winning because it was the only family flexible enough to trace the
*receiver cell-type composition* curve — the same artefact Figure 2a is about.
Once the μ_c intercepts that Section 6.1 already specifies are in the model, the
spline's advantage vanishes.

This is a correction to CS Phase 2 §8 and it should be carried into the paper as
one: the Phase 2 confounding marker was reading a confounder, but not the one it
was attributed to. It does not change Phase 2's conclusion (the sections *are* in
the confounded regime — Phase 3 §4.4 establishes that with ℓ/λ ≈ 3.5–7 and a
Ripley-K ratio of 1.11), it changes the evidence offered for it.

---

## 6. T5 — λ_proximal vs λ_downstream (Section 6.4)

**Sender definition first, because the brief is right that it is the whole
question.** `tierA_score` is `scanpy.tl.score_genes` over
`genesets/A_SENDER_FINAL_strict.txt` — the **union-strict 25-gene** Tier A set —
computed once per cell and thresholded within cell type at p95
(`code/phase2_downstream.py`). The primary call `tierA_p95` used everywhere in
Phases 3 and 5 is therefore already the union-strict call, and B1 and B4 are read
out against **one identical sender set**. No per-module sender set enters any fit
in this report. (The `A_sender_for_<module>.txt` sets exist in `genesets/` but
are not wired into the Phase 3/5 pipeline at all.)

Pooled per receiver cell type over the six sections, responses z-scored within
section, N5 + N6 applied, donor bootstrap with the **animal** as the block
(6 donors, 2,000 replicates):

| receiver | λ̂_prox | λ̂_down | ratio | donor CI on the ratio | P(λ_prox > λ_down) | SF(ctrl) B1 | SF(ctrl) B4 |
|---|---|---|---|---|---|---|---|
| Biliary/ductular | 50.0 | 50.0 | 1.00 | [0.14, 3.45] | 0.08 | 0.010 | 0.397 |
| Endothelial | 7.6 | 50.0 | 0.15 | [0.14, 1.90] | 0.09 | −0.101 | 0.159 |
| Hepatocytes | 50.0 | 7.0 | 7.14 | **[0.14, 7.14]** | 0.49 | 0.225 | 0.043 |
| Macrophages | 50.0 | 7.0 | 7.14 | **[0.14, 7.14]** | 0.57 | 0.001 | −0.512 |
| Mesenchymal | 28.4 | 7.0 | 4.06 | **[0.14, 7.14]** | 0.46 | −0.025 | 0.332 |
| T/NK cells | 9.8 | 7.0 | 1.40 | [0.23, 5.07] | 0.53 | −0.104 | 0.644 |

**The honest expected outcome is the one that occurred.** The widest ratio the
[7, 50] µm grid can express is [0.14, 7.14]; **3 of 6 donor CIs reach both ends
of it**, all six include 1, and **nine of the twelve length constants sit on a grid rail** (three of six for B1, six of six for B4).
The bound to report is therefore:

> Over six receiver cell types in mouse liver, with one union-strict sender
> definition for both programmes, we cannot distinguish λ_proximal from
> λ_downstream anywhere in the range **0.14× to 7.1×** — the entire range the
> admissible λ grid permits. Section 6.4's containment question is not answerable
> from this design.

The point estimates are not merely uncertain, they are *incoherent*: the two
cell types with the most receivers (hepatocytes, macrophages) put λ_prox at the
grid ceiling and λ_down at the grid floor, while endothelium puts them the other
way round. That pattern is what a flat profile looks like when it is forced to
pick a maximiser, and it is why the ratio must not be quoted.

---

## 7. What changes in the paper

1. **Figure 2a is replaced** with the stratified version, and the
   composition-only dotted line goes in the caption as the primary claim: *76 %
   of the published-looking contact amplitude is reproduced by receiver cell-type
   composition with no signalling of any kind.*
   **⚠ Do not use "76 %" in that caption — see the note in §4.** The sourced ratio is
   **81.5 %** (0.212 / 0.260), and the claim itself is superseded by the
   composition-matched pair: **65.9 %** (receiver cell type) to **85.4 %** (cell type +
   20-NN composition), against **1.6 %** for the same variables as matched decoys.
2. **Figure 3 exists** (`figures/figure3.{png,pdf}`) with the four Section 25
   panels: (a) 39 of 42 donor CIs on λ̂ under control span all of [7, 50] µm;
   (b) no kernel family earns its place under control; (c) Section 6.4 not
   estimable, with the grid-limit bound drawn; (d) superposition vs nearest
   against the synthetic calibration and the two nulls, on a symlog axis that
   makes the three-orders-of-magnitude gap visible.
3. **The superposition result is reported as a bounded negative at three sender
   definitions** (see §10.4 for the exact wording, which supersedes this item's
   single-caller phrasing), **not omitted.**
   It is the strongest form of the paper's argument: *we identified, on synthetic
   data, the one thing that is cleanly identifiable from a snapshot (15/15 in
   every regime); we then measured it on real tissue under control and found it
   3 orders of magnitude below the identifiable signature, with a torus shift
   reproducing most of the residual.* A project that only reported the estimands
   that failed could be accused of choosing them; this one failed on the estimand
   it pre-registered as the winner.
4. **The winner's-curse paragraph goes in the Methods or a Limitations box**,
   with the +0.056 number, its sign, and the placebo. It pre-empts the most
   obvious objection to the headline and it costs three sentences.
5. **CS Phase 2 §8 needs a one-line correction** — the spline's dominance was
   receiver cell-type composition, not window or regime.
6. **Two new methodological warnings, both cheap and both transferable**: a
   paired model comparison under a block bootstrap must share the resample; and a
   surviving fraction is not comparable across sample sizes (0.096 → 0.157 →
   0.247 as n falls by 4×, with no change in truth).

---

## 8. Reproduce

```bash
cd /workspace/code
python3 -u run_phase5_super.py   --stage section --sections inband --n-jobs 6
python3 -u run_phase5_super.py   --stage nulls   --sections inband --n-jobs 6 --n-draw 5
python3 -u run_phase5_super.py   --stage heldout --sections inband --n-jobs 8
# T1b, the same three stages at the other two sender calls (~6 min each)
for C in cdkn1a_pos senepy_p95; do
  python3 -u run_phase5_super.py --stage section --sections inband --n-jobs 6 --call $C --tag _$C
  python3 -u run_phase5_super.py --stage nulls   --sections inband --n-jobs 6 --n-draw 5 --call $C --tag _$C
  python3 -u run_phase5_super.py --stage heldout --sections inband --n-jobs 8 --call $C --tag _$C
done
python3 -u summarize_super_callers.py
python3 -u run_phase5_wc.py      --stage synth
python3 -u run_phase5_wc.py      --stage crossfit --sections inband --n-jobs 6 --frac 0.50
python3 -u run_phase5_wc.py      --stage crossfit --sections inband --n-jobs 6 --frac 0.75 --tag _f75
python3 -u run_phase5_kernels.py --stage section  --sections inband --n-jobs 6
python3 -u run_phase5_kernels.py --stage heldout  --sections inband --n-jobs 8
python3 -u run_phase5_kernels.py --stage proxdown --sections inband
python3 -u _spline_window_check.py
python3 -u _se_ratio_phase5.py
python3 -u summarize_phase5.py
python3 -u make_phase5_figs.py --which 2a,3
```

All seeded from `MASTER_SEED = 20260820`; `OMP_NUM_THREADS=1` throughout;
`cKDTree` only, no (n, n) matrix anywhere. Total wall time ~12 min on 48 cores.
Workspace footprint added: ~5 MB of CSV plus two figures.

### Engineering notes

* `BasisBlockProfiler` blocks do not have to be spatial. Making the **section**
  the block turns leave-one-section-out into an indicator multiplicity vector,
  so the entire Section 24.6 held-out analysis costs one profiler build per
  (cell type × module × basis) instead of a refit per fold.
* The superposition basis for all 40 λ values costs one sparse distance matrix
  and 40 `np.bincount` calls; `np.add.at` for the same job is ~50× slower.
* Both candidate bases must be z-scored **within section** before pooling, or the
  pooled amplitude is a sender-density comparison.
* `joblib` with the loky backend re-imports the module in each worker, so a
  module-level constant set in `__main__` (e.g. a split fraction from argv) does
  not reach the workers. Pass it as a function argument.

---

## 9. What I did NOT get to

* **BH-FDR across programs × cell types (Section 24.5).** Still deliberately not
  reported, for the Phase 3 reason: the p-values carry no information.
* **The Section 23 method baselines** — COMMOT, SpaTalk, CellChat v2, NCEM — on
  real vs torus-shifted coordinates. Figure 4 still does not exist and is still
  the highest-value remaining experiment: our own N3 returns SF = 1.000, so if
  those tools also report ~100 % signal on shifted data the panel writes itself.
* **Tier C ligand-specific kernels** (`Ccl2`→`Ccr2`, `Tnf`→`Tnfrsf1a/1b`,
  `Tgfb1`→`Tgfbr1/2`, `Il1a`→`Il1r1`), which Bio's Deliverable 5 says are the
  only within-ligand fits still defensible.
* ~~**Superposition under the other five sender calls (N7).**~~ **Done — §10.**
  `cdkn1a_pos` and `senepy_p95` are now run in full (section, nulls, held-out),
  and the verdict is partly caller-dependent. Still open: `tierA_p90`,
  `tierA_p99`, `senepy_p99` (the prevalence axis rather than the definition
  axis), and **DeepScence**, which cannot be run at all on this design —
  `deepscence_score` is absent from `senders_{section}.csv` and the arm-level
  `deepscence_{arm}.csv` files cover only 1 of the 6 admissible sections.
* **A superposition kernel with a shape other than exponential.** Section 6.3
  only asks nearest vs sum, and the sum was taken over exponentials to keep the
  parameter count identical, but a Gaussian or step superposition would complete
  the 5 × 2 grid.
* **Per-(section × cell type) Test 3 admissibility**, still applied at section
  level on `Cdkn1a`+ hepatocytes.
* **The 0.75-split extrapolation of the winner's curse to full-sample noise.**
  The +0.056 is measured at half-sample noise and is an upper bound on the
  full-sample contamination; a proper extrapolation would need a third split
  fraction and a model for how the bias scales with SE.

---

## 10. T1b — superposition across sender definitions

**This closes the Section 9 loose end flagged in §9 ("T1 used `tierA_p95` only").
The verdict does NOT transfer unchanged. Two of the four T1 claims generalise,
two are specific to `tierA_p95`.** Same six admissible sections, same receivers
(`Low_quality`/`Unknown` excluded), same [7, 50] µm grid, same N5+N6 control,
same PAIRED block bootstrap (one resample of blocks shared by both bases), same
seeds. Code: `run_phase5_super.py --call C --tag _C` (a `--tag` argument is the
only change; the fitter is untouched). Tables:
`results/phase5/super_{section,nulls,heldout}_{cdkn1a_pos,senepy_p95}.csv`,
collated in **`results/phase5/super_by_caller.csv`** and
`results/phase5/summary_super_callers.txt`.

This matters because Bio established that the sender callers agree *at or below
chance* (Jaccard ratios 0.60×/0.88×/1.66×, all |Spearman| < 0.03, surviving
matching on cell type and depth decile). They are near-independent definitions
of "sender", so a result under one of them is not a result under the others.

**DeepScence could not be included.** `deepscence_score` is not a column of
`senders_{section}.csv`, and the two arm-level `deepscence_{arm}.csv` files
cover exactly one of the six admissible sections (7259; the sham file's cell ids
match ≤ 13 cells in any admissible sham section). Running it would mean
recomputing DeepScence on five sections, which is not this task.

### 10.1 The table

Controlled design (N5 + N6), per (section × receiver cell type × module) fit.
Synthetic calibration is unchanged and is the same unit: **+15.2 / +17.0** per
1,000 cells when nearest-sender is planted, **−57.2 / −250.9** when
superposition is planted, 15/15 correct verdicts.

| | `tierA_p95` | `cdkn1a_pos` | `senepy_p95` |
|---|---|---|---|
| fits | 315 | 322 | 315 |
| **median ΔAIC/1k [IQR]** | **−0.060** [−0.238, +0.005] | **−0.321** [−0.862, −0.033] | **−0.154** [−0.665, +0.006] |
| superposition wins | 0.721 | 0.789 | 0.733 |
| paired bootstrap win fraction (median) | 0.700 | 0.860 | 0.765 |
| paired 95 % CI on ΔAIC/1k (median fit) | [−0.86, +0.24] | [−1.62, +0.19] | [−1.20, +0.22] |
| CI excludes 0 | 0.067 | 0.258 | 0.171 |
| decisive for superposition (≥ 0.95) / for nearest (≤ 0.05) | 0.076 / 0.000 | 0.335 / 0.006 | 0.225 / 0.000 |
| superposition wins, **N1** permuted senders | 0.666 | 0.638 | 0.668 |
| superposition wins, **N3** torus-shifted senders | 0.555 | 0.532 | 0.525 |
| **nearest vs covariates-only** (ΔAIC, improves) | +2.75 (0.17) | +1.57 (0.34) | +1.37 (0.38) |
| **superposition vs covariates-only** (ΔAIC, improves) | +1.59 (0.35) | **−3.91 (0.64)** | **−1.77 (0.61)** |
| held-out LL, superposition win (252 folds) | 0.508 | **0.679** | **0.603** |
| held-out LL, nearest beats covariates-only | 0.532 | 0.627 | 0.611 |
| held-out LL, superposition beats covariates-only | 0.548 | 0.663 | 0.611 |
| controlled amplitude, response-sd (sup / near) | −0.005 / −0.012 | +0.009 / +0.057 | −0.011 / −0.021 |
| median λ̂ near / sup (µm) | 14.2 / 33.4 | 15.3 / 50.0 | 12.8 / 43.0 |

Null columns are 5 draws × 6 sections each (1,575–1,610 null fits per caller);
per-draw ranges are in `super_by_caller.csv`.

### 10.2 What generalises

1. **The direction is the same at every caller**: superposition is preferred to
   nearest-sender, in 72–79 % of fits, at every sender definition and in every
   design (naive, ctrl, ctrl+N2).
2. **The magnitude is 2–3 orders of magnitude below the planted signature at
   every caller.** −0.060, −0.154 and −0.321 per 1,000 cells against −57 to −251
   when superposition is *actually* the generating mechanism. The worst case
   (`cdkn1a_pos`, p90 of the per-fit bootstrap bound, |ΔAIC/1k| ≤ 5.27) still
   excludes the planted-superposition signature by a factor of **11–48** and the
   planted-nearest signature (+15 to +17) by a factor of ~3. **The Section 6.3
   bound — dose-versus-threshold is not identifiable on this tissue at this
   effect size — holds at all three sender definitions.** This is the claim the
   paper makes, and it survives.
3. **Nearest-sender never earns its place.** At no caller does the
   nearest-sender kernel beat the covariates-only model: +2.75 / +1.57 / +1.37
   AIC, improving 17 % / 34 % / 38 % of fits. The Section 6.1 kernel is the one
   the field actually fits, and it is dominated by its own covariate block under
   every definition of "sender" we have.

### 10.3 What does NOT generalise — and this is a finding

**Phase 5 §0.1's "neither regressor beats the no-kernel model at all" and
§2.5's "held-out log-likelihood is a coin flip" are `tierA_p95` statements, not
general ones.**

At `cdkn1a_pos` and `senepy_p95` the **superposition** regressor does beat the
covariates-only model: median ΔAIC **−3.91** (improves 64 % of fits) and
**−1.77** (61 %), against +1.59 (35 %) at `tierA_p95`. The same reversal appears
on held-out sections, which is the Section 24.6 criterion and cannot be a
parameter-count artefact: superposition wins **0.679** and **0.603** of the 252
leave-one-section-out folds, against 0.508 at `tierA_p95`, and it beats
covariates-only in 66 % / 61 % of folds against 55 %.

**This is not null-reproducible, and I checked before writing it that way.** The
null win fractions are essentially caller-invariant (N1 0.64–0.67, N3 0.53–0.56
at all three callers), and on the continuous scale the null ΔAIC/1k medians are
−0.018 to −0.024 (N1) and −0.003 to −0.005 (N3) everywhere. So at `cdkn1a_pos`
the observed −0.321 is **17× the N1 null and 65× the N3 null**, and at
`senepy_p95` the observed −0.154 is 7× and 50×; at `tierA_p95` the observed
−0.060 is only **2.7× and 5.5×**. The same holds for "beats no kernel":
superposition improves 64 % / 61 % of real fits against 21–27 % of null fits.
**At `tierA_p95` the preference was mostly reproduced by the nulls; at the other
two callers it is not.** The Phase 5 §2.6 sentence "most of the preference is
reproduced by nulls that contain no signalling" is therefore true at
`tierA_p95` and false at `cdkn1a_pos`/`senepy_p95`, and the report should not
have generalised it. It now does not.

Two things this is **not**:

* **It is not gene-set circularity.** `Cdkn1a` is a member of 4 of the 7 Tier B
  response modules, so the obvious explanation for `cdkn1a_pos` is that the
  sender flag is one of the response genes. It does not survive the split: in
  the four modules containing `Cdkn1a` the controlled ΔAIC/1k is −0.254 and
  superposition beats no-kernel in 0.61 of fits; in the three that do **not**
  contain it (`emt_ecm`, `il6_jak_stat3`, `oxidative_stress`) it is **−0.370**
  and **0.67**. The effect is if anything *stronger* where the circularity is
  absent, and `senepy_p95` — which has no `Cdkn1a` in its sender call at all —
  shows the same split in the same direction (−0.094 → −0.267). (`Cdkn1a` is not
  in `A_SENDER_FINAL_strict.txt` either, so `tierA_p95`, the caller with the
  *weakest* effect, is also one without the overlap.) Numbers in
  `results/phase5/super_by_caller_diag.csv`.
* **It is not a contact-scale dose response.** In the fits where superposition
  beats covariates-only, λ̂_sup sits on the **50 µm grid ceiling in 57–60 % of
  fits at all three callers (0.600 / 0.602 / 0.565), median λ̂_sup = 50.0 µm** —
  the broadest smoothing
  the admissible grid allows. What the superposition regressor is adding over
  the N5 block is a *regional* sender-density field, not a signalling kernel at
  the length scale Section 6 is about. That is consistent with §2.6's reading
  (the regressor is largely a smoothed local density) but sharpens it: the
  smoothing that helps is at the top of the grid, i.e. the model wants a
  covariate, not a kernel.

The mechanical reason the two callers differ from `tierA_p95` is visible in the
sender geometry: `cdkn1a_pos` and `senepy_p95` are sparser and more dispersed
(median distance-to-nearest-sender 42.7 and 39.4 µm against 30.2 µm at
`tierA_p95`, `results/phase3/window.csv`), so a 50 µm-smoothed sum over them
carries regional structure that the N5 density block (25/50/100 µm densities of
*all* cells, plus 1-NN distance) does not already contain. At `tierA_p95` the
senders are dense enough that the N5 block absorbs it.

### 10.4 What to write in the paper

The T1 paragraph must be re-scoped, and it becomes a **better** paragraph, not a
worse one, because the caller dependence is one of the paper's own themes
(Section 9 A7, *your sender call is a choice, not a fact*):

> The one estimand our synthetic study identified as cleanly identifiable
> (15/15 in every regime) separates the two kernels, on real tissue under full
> control, by 0.06–0.32 AIC units per thousand cells at three near-independent
> sender definitions — two to three orders of magnitude below the −57 to −251
> planted signature, with the per-fit block-bootstrap bound excluding that
> signature by ≥ 11× even at the worst caller's 90th percentile fit. Whether the
> *superposition regressor earns its place at all* is itself sender-definition
> dependent: it does not beat a covariates-only model under the Tier A call
> (35 % of fits, held-out win 0.508) and does under Cdkn1a⁺ and SenePy
> (61–64 %, held-out win 0.60–0.68), while the nearest-sender kernel — the one
> the field fits — never beats it under any call. Where the superposition term
> does help, its length constant rails at the top of the admissible grid, i.e.
> it is acting as a regional density covariate rather than a contact-scale dose.

Concretely:

1. **Keep** the headline bound (§0.1) — it is now established at three callers
   and is stronger for it.
2. **Re-scope** the "neither kernel beats no kernel" sentence to the
   nearest-sender kernel (true everywhere) plus a caller-dependence clause for
   superposition.
3. **Re-scope** the "not stable under nulls" sentence (§0.2) to `tierA_p95`, and
   report the null-excess ratios (2.7×/5.5× at `tierA_p95` vs 17×/65× at
   `cdkn1a_pos`) rather than the win fractions alone — the win fraction is the
   less informative of the two statistics and it is what made the three callers
   look alike.
4. **Add one row to the caller-dependence table in Section 9 A7**: model
   selection between kernel families is sender-definition dependent, with the
   held-out win fraction moving 0.508 → 0.679 across callers that agree at
   chance.
5. Figure 3d gains two more (observed, N1, N3) triples — the panel is already
   built to take them.

> **Applied 2026-08-21.** All five are done. (1)–(4) are in `README.md`: the
> superposition paragraph is re-scoped to the nearest-sender kernel with a
> caller-dependence clause, the null-excess ratios (2.7×/5.5× at `tierA_p95` vs
> 17×/65× at `cdkn1a_pos`) are reported alongside the win fractions, and the
> caller-dependence table gained the model-selection row (held-out win fraction
> 0.508 → 0.679). (5) `code/make_phase5_figs.py` now builds panel (d) from all
> three callers — twelve boxes, `figures/figure3.{png,pdf}` regenerated.
