# How Far Does Senescence Actually Reach?

**Confounding and identifiability in estimating SASP spatial response kernels.**

Target: NeurIPS 2026 workshop — ml4spatialbio (primary), ICBINB-BIO (secondary).
Deadline: 2026-08-29 AoE. Full spec in [`SASP_Kernel_Master_Plan.md`](SASP_Kernel_Master_Plan.md).

---

## The question

A senescent cell secretes SASP factors that alter its neighbours, and that
influence must fall off with distance. This has been shown descriptively
(Ma et al., *Cell* 2024;187(24):7025–7044.e34, PMID 39500323 — commonly
miscited as "Zhao et al."; Zhao L is 14th of 47 authors) — and shown on
**Stereo-seq bins, not segmented cells**: 1,535,191 spots at ~1,450
genes/spot, in young versus old male mice across nine tissues (CNGB STOMICS
`STDS0000247`). The closest prior distance-gradient result in this field was
therefore measured at a resolution that blurs the very distance it reports.
It has never been estimated as a parameter with error bars,
never tested against confounder-aware nulls, and never subjected to an
identifiability analysis — despite contemporaneous benchmarks reporting that
the class of methods used to make such claims returns comparable interaction
counts on real data and on data whose ligand–receptor proximity has been
destroyed by permuting cell locations *within cell type*.

This is an **evaluation and identifiability** project, not a new-method project.

---

## Status

| Phase | What | State |
|---|---|---|
| 0 | Data acquisition, storage, environment | complete — [report](reports/PHASE0_DATA_AND_ENV.md) |
| 1 (Bio) | Day-1 audit, gene sets, Test 1–3 | complete — [report](reports/BIO_PHASE1.md) |
| 1 (CS) | Synthetic identifiability study, Figure 1 | complete — [report](reports/CS_PHASE1.md) |
| 2 (Bio) | Cell typing, anatomy, sender calls, modules | complete — [report](reports/BIO_PHASE2.md) |
| 2 (CS) | Naive estimation on real tissue, Figure 2a | complete — [report](reports/CS_PHASE2.md) |
| 3 (Bio) | All 11 sections annotated, ligand–receptor audit, caller analysis | complete — [report](reports/BIO_PHASE3.md) |
| 3 (CS) | Null battery N1–N8 on real data — **the go/no-go** | complete — [report](reports/CS_PHASE3.md) |
| 4 (CS) | Figure 4: existing CCC tools under the same nulls | complete — [report](reports/CS_PHASE4.md) |
| 5 (CS) | Superposition, winner's curse, kernel families, Figures 2a + 3 | complete — [report](reports/CS_PHASE5.md) |
| 6 (Audit follow-up) | D7 corrections applied; CellWHISPER's own null run; environment re-verified | complete — [CS_PHASE4 §8](reports/CS_PHASE4.md), [CS_PHASE1 §8](reports/CS_PHASE1.md) |

Master-plan gates: **Day 3 met** (synthetic recovery + audit).
**Day 6 met three days early** (Figures 1 and 2a exist). The null battery, budgeted
for Days 6–7, is also complete.

All **11 liver sections** (11 animals, 2/10/26/52 weeks, SBR vs sham) are processed
and annotated.

---

## Headline findings so far

All numbers are reproducible from `results/`; see the linked reports for method.

**The kernel does not survive the null battery.** Across **153 of 315**
reportable (section × receiver type × module) fits — not every receiver type is
present in every section — the surviving fraction
under combined control (N2+N5+N6) is **0.088**, **IQR across fits
[−0.017, 0.234]**, with 30% at or
below zero. 91% of the naive amplitude is nuisance. Controlled amplitude is
0.029 response-sd against a naive 0.329.
**Bound: controlled amplitude ≤ 0.18 response-sd** (80% power at 0.183 sd).

> *Updated 2026-08-27 (record reconciliation).* These were the **pre-C6** figures —
> 160 of 315, SF 0.082 [−0.099, 0.249], 35 % ≤ 0, controlled 0.027, naive 0.326,
> bound 0.203. The frozen post-C6 vector is `results/phase3/m1_final_audit.txt` §3.
> **The bracket is an inter-quartile range across fits, not a confidence interval**
> — `results/phase3/sf_summary.csv` carries `q25 / median / q75` and no CI column.

**That is a conclusion, not a failure to detect — because the synthetic study
calibrates it.** With a *planted* effect, N5 leaves 0.826 [0.734, 0.921]; with
*no* effect it leaves 0.412. The observed value under **N5 alone, 0.084** — the
like-for-like comparison, since the synthetic calibration is N5 — is below 100%
of the 600 planted-effect runs and sits at the 23rd percentile of the no-effect
runs. (The **0.082** above is the combined N2+N5+N6 figure; `CS_PHASE3.md` §0
carries both.) The
tissue behaves like the synthetic no-effect case.

**The estimator is correct.** With a planted kernel and no nuisance,
λ̂ = 30.53 [29.14, 31.91] against λ_true = 30, CI coverage 0.95.

**No single control is safe across regimes, and the regime is not identifiable
from the data.** Matched-decoy (N2) does not restore calibration; nuisance
conditioning (N5) swaps positive bias for negative and is *worse than nothing*
at zero confounder strength. The 4-estimator spread brackets λ_true in only
25% of grid cells — it is a lower bound on uncertainty, not a coverage interval.

**The torus shift certifies confounding as real.** With β_true = 0 in the
confounded regime, N3 returns a **98%** surviving fraction (N4 rotation 89%),
versus 99.5% when β_true = 1 — it cannot tell signal from confounding. N1
(cell-type-stratified label permutation) is the only null that discriminates
(0.006 vs 0.69). All three reject at 42–100% under a true null, so **p-values
are uninformative here; report surviving fractions.** N3 is the null the
current reliability literature recommends.

**This reproduces on real tissue: N3 surviving fraction 1.000 [0.992, 1.008]**,
with 88% rejecting at p < 0.05. The mechanism is visible in the data — the
torus-shifted null is centred at ~0 (`N3_null_mean` 8×10⁻⁶ against
`beta_obs` 3.4×10⁻³), so it reproduces none of the confounding. Shifting senders
to arbitrary coordinates destroys the sender–neighbourhood association that
generates the spurious signal, making the null far too weak. N1, which keeps
senders at real cell locations within cell type, recovers ~71%.

**The established tools cannot tell this tissue from confetti.** COMMOT (run as
released software), plus faithful, clearly-labelled reimplementations of the
CellChat v2\*, SpaTalk\* and NCEM-linear\* statistics, on the six admissible sections
and the four ligand–receptor pairs this panel supports. Of the interactions each
calls significant on real coordinates, **79 % (COMMOT), 94 % (CellChat\*) and 79 %
(SpaTalk\*) are still called significant after every cell coordinate is
permuted** — and their scores rank-correlate at ρ = 0.86–0.98 with their
real-data values on that same destroyed tissue. The failure is not a weak null:
there is nothing left in N0 to destroy. The parsimonious reading is that the
statistic carries little or no spatial information.
For COMMOT we show the mechanism directly — coordinate permutation replaces its
entire cell-to-cell communication network (Jaccard 0.015) while conserving
transported ligand mass to six figures, so the cell-type-level summary it
actually tests barely moves (ρ = 0.78 — the per-pair cluster-level Spearman on
the mechanism tile; the ρ = 0.90 quoted later is the different, all-interaction
figure over 6,032 interactions), and its permutation test holds the transport
plan fixed. Synthetic positive controls rule out the obvious
objection: the same implementations find a planted interaction at p = 0 with
1-in-16 specificity and lose it correctly when the two cell types are separated —
but hold the ligand and receptor rates fixed and change only whether the
receptor⁺ cells are the ones next to the ligand⁺ cells, and CellChat\*'s score is
identical to six significant figures. NCEM-linear\* is the only one of the four
that sees the difference (p = 10⁻³⁹ vs 1.0) and the only one correctly calibrated
under the nulls, but its reported interaction length scale is not identified.
Details in [CS_PHASE4](reports/CS_PHASE4.md); Figure 4.

**Naive CIs are indefensible on this data.** iid standard errors understate by
up to 7.9× in synthetic and 4–5× on real tissue.

**The naive real-data curve looks exactly like published gradients** — monotone
in 32/35 fits, effect 26× the bin SEM — while λ̂ pins to a grid bound in **63%**
of fits (200 of 315 on the final admissible set; the 66% figure is the earlier
naive Phase 2 population), kernel family moves d̂½ by a median 5× naive and
**4.4× under control**, and adjusting for counts, area and
density *alone* removes 86% of β̂.

**The one estimand that was cleanly identifiable in synthetic tissue is not
identifiable here either.** Superposition vs nearest-sender (Section 6.3) was
15/15 correct in every synthetic regime. Under full control on real tissue the
median ΔAIC is **−0.06 per 1,000 cells** against a planted-superposition
signature of −57 to −251, held-out log-likelihood on left-out sections is a coin
flip at `tierA_p95` (50.8 % of 252 folds; 0.679 and 0.603 at the two other
callers), and the nulls reproduce most of the residual preference (superposition
"wins" 72 % of fits observed, 67 % on permuted senders, 56 % on shifted senders
at `tierA_p95`). **The nearest-sender regressor — the one the field fits — never
beats the no-kernel model under any sender call.** Whether the *superposition*
regressor earns its place is itself sender-definition dependent, which is a
result rather than an inconsistency; see below.

Report the magnitude, not the win fraction. The observed ΔAIC/1k is only
**2.7× / 5.5×** the N1 / N3 null at `tierA_p95` but **17× / 65×** at
`cdkn1a_pos`, while the null win fractions are near-identical at all three
callers (N1 0.64–0.67, N3 0.53–0.56). The win fraction hides the caller
dependence; the per-1,000-cell magnitude shows it.

**The negative is not a selection artefact.** Cross-fitting the Phase 3
significance screen (select on half the spatial blocks, estimate on the other
half) puts the winner's curse at **+0.056** of surviving fraction — and in the
*upward* direction, so correcting for it makes the negative stronger. A
surviving fraction is also not comparable across sample sizes: the same fits give
SF(N5) 0.096 / 0.157 / 0.247 at 100 / 50 / 25 spatial blocks.

**Sender signal is far weaker than the confound it must be separated from.**
k-NN spatial autocorrelation: zonation 0.379–0.752, `Cdkn1a` 0.0085–0.0194,
Tier A sender score 0.0040–0.0124, permutation null ~0.0006. Real, but one to
two orders of magnitude below the confound.

**Zonation is not the liver confound.** The master plan rated this the
highest-probability failure mode for liver. The zonation covariate alone removes
~0% in hepatocytes (SF 1.043), and the kernel does not vanish within zones
because zonation does not appear to be driving it. The confound is dominated by
technical and geometric terms:
transcript depth and cell size alone leave 0.288, local density 0.219, kNN
composition 0.474, anatomy 0.810, segmentation method 0.998. Separately, **66%
of the unstratified gradient is receiver cell-type composition** — a
composition-only surrogate, replacing each cell's response with its cell type's
section mean, reproduces **76%** of the contact amplitude (median; per-module
ratios span **0.42–1.83**) with no signalling whatsoever. Stratified within cell type the amplitude falls from +0.260 to
+0.091 response-sd, and two of nine receiver types go *negative*.

**Matched decoys balance beautifully and bound nothing.** |SMD| ≤ 0.033 after
matching, surviving fraction 0.943. The synthetic study explains why it carries
no information: N2 returns 0.934 under a real effect and 0.775 under none — it
barely moves either way.

**The regressor is largely a sender-calling-rate readout.** Across 77
section × sender-definition combinations, log(median distance to nearest sender)
on log(sender density) gives slope **−0.524, r² = 0.984** against a Poisson
prediction of −0.50; Ripley-K at 50 µm is 1.11 at `tierA_p95` (1.26 at
`cdkn1a_pos`, 1.56 at `senepy_p95`), so senders are nearly a random thinning at
every caller. Density normalisation does *not* rescue λ̂ (between-section sd of
log λ̂: 0.726 raw, 0.705 Poisson-normalised). Sender prevalence in turn tracks
section median transcripts/cell at ρ = +0.94.

**No kernel earns its place at a biological length scale.** Under full control
the nearest-sender regressor never beats a covariates-only model at any sender
definition (ΔAIC +2.75 / +1.57 / +1.37). The superposition regressor does beat
it under two of the three near-independent sender calls (`cdkn1a_pos` −3.91,
64% of fits, held-out win 0.679; `senepy_p95` −1.77, 61%, 0.603) but not under
`tierA_p95` (+1.59, 35%, 0.508) — a sender-definition dependence that is itself
a result (Section 9 A7). Where it does win, **λ̂ rails at the 50 µm grid ceiling
in 57–60% of fits**: it is acting as a regional sender-density covariate, not a
signalling kernel. It is not `Cdkn1a` circularity — the effect is *stronger* in
the modules that exclude `Cdkn1a`, and `senepy_p95` reproduces the split.
Superposition vs
nearest — the one estimand the synthetic study showed is *cleanly* identifiable
(15/15 correct, |ΔAIC| 15–251 per 1,000 cells when a truth is planted) — comes
out at **−0.060 per 1,000 cells** on real tissue, two orders of magnitude below
the smallest planted signal. Held-out log-likelihood makes it a coin flip
(superposition wins 50.8% of 252 folds). And the preference is not stable under
nulls: superposition "wins" 72.1% observed, **66.6% on permuted senders and
55.5% on torus-shifted senders**, so it is tracking geometry — the superposition
regressor is largely a smoothed local density (r = 0.66 with the nearest-sender
regressor).

**The negative is not a winner's-curse artefact — selection biases it the other
way.** Cross-fitting with a matched no-selection placebo puts the contamination
at **+0.0563**, i.e. the in-sample surviving fraction is *inflated*: β̂_N5 and
β̂_naive share sampling noise, so selecting a positive noise draw lifts numerator
and denominator alike and pulls the ratio toward 1. Two independent split
fractions give raw gaps of +0.058 and −0.073 yet the same corrected +0.0563,
because the placebo absorbs SF's sample-size dependence. Correcting would make
the negative *stronger*. (Caution for anyone reusing SF: it is not comparable
across sample sizes — 0.096 / 0.157 / 0.247 at 100 / 50 / 25 blocks.)

**Established CCC tools do not fail our way — they fail worse.** Run on the
same tiles with coordinates fully permuted, so that *no spatial information
remains*, the three ligand–receptor methods' scores still rank-correlate with
their real-data values at **ρ = 0.90 (COMMOT), 0.98 (CellChat v2\*), 0.86
(SpaTalk\*)**, and 79% / 94% / 79% of real-data-significant interactions stay
significant. Their significance *rates* on randomised input match the real ones
(COMMOT 0.241 vs 0.224; CellChat\* 0.318 vs 0.283; SpaTalk\* 0.248 vs 0.247),
consistent with CellWHISPER's criterion on senescence-relevant pairs.

**And it replicates under CellWHISPER's *own* null, which we now run.** Their
control permutes locations *within each cell type* — preserving every cell
type's spatial organisation and destroying only ligand–receptor proximity —
whereas our N0 permutes across all cells and is strictly more destructive. Under
their design (`N0_type`), survival is **higher**, as a less destructive null
should be: **COMMOT 0.811, CellChat v2\* 0.971, SpaTalk\* 0.781** (NCEM linear\*
0.009). CellChat\*'s significance *rate* is **0.2831 on permuted tissue against
0.2833 on real tissue** — their criterion met to four decimal places — and its
score ordering is preserved **exactly**, ρ = **1.000**. This removes the
objection that the earlier numbers came from an unreasonably blunt shuffle:
there is little spatial information in these statistics for any null to remove.
Details in [CS_PHASE4](reports/CS_PHASE4.md) §8; Figure 4, orange bars.
We report that these three methods fail CellWHISPER's null — not that
CellWHISPER passes it here, which would need running their tool.

**The statistic carries little spatial information — and for COMMOT we show the
mechanism in the published code.** Coordinate permutation replaces **98.5%** of the cell-to-cell
communication network (edge Jaccard **0.014–0.018**, mean 0.0154) while the collective optimal
transport **conserves transported ligand mass to seven significant figures**
(e.g. 863.815754 → 863.815745). Averaging that flow over a sender × receiver
cell-type block therefore returns ≈ mass × composition, which is geometry-free;
`cluster_communication` then permutes cell labels with the transport plan held
fixed, so it never tests space at all. Synthetic positive controls close the
loop: holding ligand/receptor rates fixed and changing *only* whether receptor⁺
cells sit beside ligand⁺ cells leaves CellChat\*'s score identical to six
significant figures and COMMOT's to four. NCEM-linear\* is the only one of the
four that separates them (p = 10⁻³⁹ vs 1.0) and the only one calibrated under
these nulls — though its own reported length scale is not identified.

**Two opposite failures, two opposite repairs.** Our surviving fraction of 1.000
arises because the torus shift *annihilates* β̂ (null centred at 0.8% of the
observed value) and the test rejects anyway: the null hypothesis "sender field
unaligned with response field" is false under pure confounding, so we need a
**confounder-preserving null**. These tools' score surviving fractions of
0.91–1.00 arise because the shift does not touch their statistic at all: they
need a **different statistic**, or a null on the transport plan / neighbour
graph. Reading "fails the torus shift" as "use a stronger null" fixes nothing
for COMMOT.

\* reimplementation of the published statistic, not the published software —
see [CS_PHASE4.md](reports/CS_PHASE4.md). COMMOT 0.0.3 is the published package.

**Segmentation is not the problem here.** Section 8 Test 1's assignment rate,
measured on 129,104,526 transcripts from admissible section 7259: **88.27%
assigned to a cell, 11.72% unassigned** (11.72% at Q≥20; 81.68% of transcripts
pass Q≥20; 33.82% of assigned transcripts are nuclear). The plan's >30%
bleed-through threshold is comfortably cleared, so the confounding documented
here is *not* attributable to transcript bleed-through.

**Senescence calling is a choice, not a fact — but the callers are not
independent, and saying they were was wrong.** Four sender definitions —
DeepScence (run with `denoise=False`, a forced deviation from the published
configuration, because its DCA dependency needs an obsolete TensorFlow stack,
now scored on **all eleven sections**), SenePy, `Cdkn1a`⁺ and a curated
arrest-and-damage score — overlap **above chance** after conditioning on cell
type *and* transcript-depth decile: pooled **1.21× chance (z = 21.9,
p = 1.8 × 10⁻¹⁰⁶)** over the four non-circular pairs, eleven sections, frozen
strict-33 Tier A (`results/phase3/caller_coverage_gate_headline.csv`). The
effect is **small and certain**: 21 % above chance, not 2× above it.

> **Superseded, and by how much.** This paragraph read "**0.93–1.22× of chance
> for four of six pairs** … i.e. statistically independent", from a
> **two-section** DeepScence base and the pre-C6 25-gene Tier A. Both premises
> moved. At full coverage on the same gene sets the pooled value went
> 1.040 → 1.129; on the frozen 33-gene Tier A it is 1.212. It is above chance on
> the *published two-section base* too (**1.131×, p = 6.5 × 10⁻⁹** — the 4-pair basis,
> the same basis as the 1.212 above; *corrected 2026-08-27, this read "1.13×,
> p = 4 × 10⁻⁸", which pairs a 4-pair ratio with a 3-pair p-value. On the 3-pair basis
> the consistent pair is 1.128, p = 4.4 × 10⁻⁸ → 1.212, p = 1.8 × 10⁻⁹⁴*), so coverage is not
> the only reason. The defensible restatement is **"weakly but genuinely
> dependent, in a direction each pair's depth loading explains"**
> (`reports/CS_PHASE8_CALLERS.md` §3, `reports/CORRECTIONS.md` §2).
>
> Two per-pair numbers in the old text are also superseded. **DeepScence vs
> `Cdkn1a`⁺ was quoted at 1.51–2.85×**; over eleven sections it is 0.963–2.849
> with a median of **1.071** and a pooled **1.255**, so the quoted range
> overstates the measured circularity by about two-fold. The pair is still
> circular by construction — `CDKN1A` is DeepScence's own sign anchor — and is
> still excluded from every pooled number. **SenePy vs DeepScence was read as
> "concordant in sham (2.15×), anti-concordant in SBR (0.38×)"**; at full
> coverage it is 0.33–0.55× in **ten of eleven** sections. There is no arm
> effect — there is one anomalous section, 7250. DeepScence's CoreScence set is **79%** circular
with our response modules (26 of the 33 CoreScence genes reachable on the mouse
panel are in ≥1 Tier B module; corrected 2026-08-27 from "69%", whose denominator
of 35 was a typed-in literal — `reports/BIO_PHASE2.md` §4.2 correction box). Because of the `denoise=False` deviation — ~~DCA
denoising is precisely the step that would normalise depth~~ — these characterise
**DeepScence as we could run it on this panel, not DeepScence as published**.

> **Corrected 2026-08-27.** The struck rationale is **refuted by measurement**
> (`PREREG_PHASE8.md` P29). DCA denoising does not normalise depth — it **raises**
> the depth loading, on three of three sections and both arms: Spearman ρ of score
> against transcript counts 0.3891 → **0.6404** (7239, ×1.65), 0.3176 → **0.5314**
> (7259, ×1.67), 0.4096 → **0.5419** (7352, ×1.32) — i.e. **×1.32–1.67**, not
> "roughly doubles" either (`results/phase8_d2/d2_depth.csv`, columns
> `rho_depth_committed` / `rho_depth_alt`, `config == "dca"`). **Whatever DCA
> contributes here, it is not depth normalisation.** The `denoise=False` choice
> itself is unaffected and is now a *chosen* value rather than a forced deviation
> (P22): DCA installed and ran. The sentence that survives is the one after the
> dash — these characterise DeepScence as we could run it on this panel.
What the choice propagates into:

| Quantity | Moves across callers that agree **weakly but genuinely above chance** ***[header corrected 2026-08-27 — it read "callers that agree at chance", which restates the superseded claim the correction box three paragraphs above retires. Frozen: pooled **1.212× chance, z = 21.92, p = 1.84e-106** over 11 sections, `results/phase3/caller_coverage_gate_headline.csv`]*** |
|---|---|
| Controlled surviving fraction | varies 3–5×; no caller exceeds 0.29 |
| Does the superposition regressor beat covariates-only? | **no** at `tierA_p95` (35 % of fits) — **yes** at `cdkn1a_pos` / `senepy_p95` (61–64 %) |
| Held-out win fraction, model selection between kernel families | **0.508 → 0.679** |
| Observed ΔAIC/1k as a multiple of the N1 / N3 null | 2.7× / 5.5× → 17× / 65× |

Model selection between kernel families is therefore sender-definition
dependent, not just parameter estimation. The bound itself is not: it holds at
all three callers.

---

## What this paper contributes

Framed for ml4spatialbio, whose call asks for *benchmarks and evaluation
standards specific to spatial tasks* and *uncertainty-aware spatial models
biologists can trust*. The contribution is the evaluation framework; the
negative senescence result is the demonstration that it works, not the headline.

1. **An identifiability regime map with planted ground truth.** Where a spatial
   response kernel is recoverable from a static snapshot and where it is not, as
   a function of sender clustering, baseline autocorrelation length, sender
   prevalence and confounder strength — with bias *and* CI coverage reported.
2. **A calibration protocol for null batteries.** Interpret a real-data
   surviving fraction by reference to synthetic runs with a *known* planted
   effect and with none. This converts "the effect shrank under control" into a
   quantified statement about which hypothesis the data resemble. We are not
   aware of prior work in spatial CCC doing this, and it is reusable by anyone.
3. **Evidence that a widely recommended null is a calibration failure.** The
   torus shift returns a surviving fraction of ~1.0 both in synthetic data with
   no planted effect and on real tissue, with the mechanism identified: the
   shifted null is centred at zero because relocating senders destroys the
   sender–neighbourhood association that produces the spurious signal.
4. **A structural identifiability result about the estimand itself.**
   Distance-to-nearest-sender is a sender-calling-rate readout to r² = 0.98,
   and sender rate tracks per-cell detection depth (median transcripts/cell —
   not "sequencing depth"; this is an imaging assay) at ρ = 0.94 within the SBR
   arm, n = 6. Pooled across all 11 sections ρ = 0.16, so this is suggestive,
   not established. This applies to every
   method in this literature that regresses response on distance-to-nearest-X.
5. **A reusable null battery and a negative result with a quantified bound** on
   1,834,806 cells across 11 sections from 11 animals (1,036,459 in the six
   Test-3-admissible sections).

---

## Dataset

**GEO GSE310392** — the *mouse* arm of Karpova et al., *Cell Genomics*
2026;6(2):101133 (PMID 41576948). The human arm is on the SenNet portal and the
mCRC arm on HTAN; only the mouse data went to GEO.

- Xenium Prime Mouse 5K, **5,106 Gene Expression features** (verified from the
  matrices, not the platform record)
- Mouse IFALD model: 75% small-bowel resection (`sbr`) vs `sham`, male C57BL/6J
- Timepoints 2 / 10 / 26 / 52 weeks; **12 sections from 11 animals**
- 1,834,806 cells total; 1,036,459 in the six Test-3-admissible sections
- median NN cell distance 6.7–9.7 µm (the resolution floor)

Note the plan's Section 7 originally described this accession as a 43-donor
human aging atlas. That was wrong — right paper, wrong accession — and
Section 7 needs rewriting accordingly.

---

## Repository layout

```
code/           simulator, estimators, null battery, gene-set builder, fetch scripts
code/_shims/    compatibility shims (numpy-2 aliases for commot; DCA for DeepScence)
genesets/       Tiers A–E, one mouse symbol per line, + archived MSigDB JSON
results/        sweep outputs, fits, null distributions, summary tables
figures/        Figures 1, 2a-2d, 3, 4 (+ underlying CSVs)
reports/        phase reports — the project narrative and every number's provenance
references.bib  all 30 references in the D7 audit's corrected form, with the
                audit's findings carried inline as % AUDIT: flags
data/           NOT tracked; see below
```

Not tracked, and why: `data/raw` and `data/interim` (re-downloadable from GEO);
`data/processed/cache3` (166 MB of per-cell `.npz`, regenerable by
`code/sasp_phase3.py`); the Phase 4 job checkpoints under `results/phase4/parts`
(~650 MB, regenerable by `code/phase4_*.py` — the aggregated tables beside them
carry every number in the report); and `logs/`. Everything a number in a report
depends on **is** tracked.

---

## Reproducing the data

Do **not** download `GSE310392_RAW.tar`, and do not download the whole
per-sample tarball. Each 12–31 GB archive contains a `morphology.ome.tif`
(4–11 GB), a `transcripts.parquet` (1.9 GB) and sometimes an 11.7 GB processed
`.rds` — none of which are needed. The three files that matter total 76–130 MB.

```bash
bash code/fetch_xenium_bundle.sh GSM9295284 \
  GSM9295284_7250_liver_sham_Male_26-U1.tar.gz /workspace/data/raw
```

It streams the archive, extracts only `cell_feature_matrix.h5`, `cells.parquet`
and `cell_boundaries.parquet`, and tears the pipeline down as soon as they stop
growing. Bandwidth is traded for disk deliberately: egress is free on RunPod and
the network volume is quota-capped.

Gene sets rebuild offline from the archived MSigDB JSON — no network needed:

```bash
python3 code/build_genesets.py
```

## Environment

There are **two** interpreters in this project, and the difference matters.

**1. The main environment** — a pip venv at `/workspace/envs/sasp311`, Python 3.11,
described exactly by `requirements.txt` (144 packages; the 21 `==` pins are the
load-bearing ones):

```bash
python3.11 -m venv /workspace/envs/sasp311
/workspace/envs/sasp311/bin/pip install -r requirements.txt
```

Every driver in `code/*.sh` sources `code/_env.sh`, which resolves that
interpreter (`SASP_PYTHON`), puts it first on `PATH` so a bare `python3` inside a
driver reaches it too, and **fails loudly if it is missing** rather than falling
back. Run Python scripts directly the same way:

```bash
/workspace/envs/sasp311/bin/python code/summarize_phase3.py
# or:  . code/_env.sh && python3 code/summarize_phase3.py
```

Do **not** use the system `python3`. On the development container `/usr/bin/python3`
loads `/usr/local/lib/python3.11/dist-packages`, a *second* complete scientific stack
of ~249 packages that no manifest here describes and that has been wiped twice. Every
pin currently matches in both stacks, so nothing has diverged numerically yet — but
only the venv is described by `requirements.txt`, and only the venv is rebuildable.
Override with `SASP_PYTHON=/path/to/python` if your environment lives elsewhere.

Eight scripts still read package data by absolute path from the overlay
(`/usr/local/lib/python3.11/dist-packages/{DeepScence,senepy}/data/...`); see
`reports/REPRODUCIBILITY_REPAIR.md`. `pip install -r requirements.txt` into a venv does
**not** populate that path, so those eight need `DeepScence`/`senepy` importable from
the venv and the constant edited, or the packages present at that path.

**2. The DCA environment** — a separate CPython 3.8 venv, needed only for the D2
`denoise=True` arm (TensorFlow < 2.5 ships no wheel for 3.11):

```bash
DCA_ENV_ROOT=/some/scratch bash code/setup_dca_env.sh
export DCA_ENV_ROOT=/some/scratch      # the D2 runners honour this
```

Python 3.11 for everything else. Developed on 48 cores / 251 GB RAM (the container
cgroup, not `free`, is the real ceiling). The workload is CPU- and RAM-bound; no GPU
is used at any point.

---

## Reproducibility notes

- All seeds are pinned; sweeps checkpoint incrementally to `results/`.
- Gene sets are pinned to **MSigDB 2026.1.Mm**, fetched 2026-08-20, with the
  raw JSON archived in `genesets/msigdb_mouse_2026.1.Mm/` so results do not
  drift with MSigDB releases.
- Cell-type annotation is a function over a sample directory with a `--relabel`
  mode that reuses saved Leiden clustering.
- Neighbour work uses `scipy.spatial.cKDTree` throughout. An (n,n) distance
  matrix is never constructed — a 238k-cell section is 2.8 × 10¹⁰ pairs.

## Claim audit

An independent audit ([D7](reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md)) checked all
30 references in the master plan and every claim in this repository.
**No reference was fabricated — all 30 exist.** 12 verified cleanly, 5 were
incomplete but not wrong, and 13 carried a genuine error: 10 bibliographic and
3 where the cited paper does not say what was attributed to it.

The structural finding is worth stating plainly: **the phase reports hedge
carefully and volunteer limitations against their own interest; almost every
overreach in this project was introduced by this README compressing them.**
Ten such errors have been corrected here, including a wrong first author on the
closest prior work, an unsupported cell count, a withdrawn arm restriction, and
an FPR claim that one of our own reports explicitly disclaims.

**Closed since the audit** (2026-08-21):

- ~~The four containment mechanisms attributed to Martin et al.~~ **Corrected.**
  They are the master plan's gloss, not that paper's — `degrad` and `refractor`
  occur zero times in the full text. D6 §A.3 is re-attributed, a new §A.3(e)
  treats Martin et al.'s three *actual* resolutions (SASP-poor secondary cells as
  a firebreak, time-delayed induction, secondary cells secreting less), and the
  master plan §3/§6.4/§12 carry the correction inline. No result moved.
- ~~DeepScence claims need the `denoise=False` caveat~~ **Restated above**, with
  the two-of-eleven-sections scope and the superseding depth-matched ratios.
- ~~The CellWHISPER null is mislabelled~~ **Fixed by running theirs.** Their
  control permutes locations *within cell type*; our `N0_perm` permuted across
  all cells and is strictly more destructive. We now run the real thing —
  `N0_type` — for all four methods. See [CS_PHASE4](reports/CS_PHASE4.md) §8.
- ~~Environment drift unverified~~ **Verified.** See the environment note below.
- ~~Ma et al. described only as a distance-gradient result~~ **Corrected (B2/B3).**
  It is **Stereo-seq bins, not segmented cells**, and **two age groups**, not a
  time course. Stated in the opening, and in master plan §4.1 / §7 Rank 3 — which
  also had the data pointed at GSA rather than CNGB STOMICS. This *helps* the gap
  argument, so it is now said out loud rather than left implicit.
- ~~Numeric inconsistencies across reports (M1–M10)~~ **Reconciled in this file.**
  160 → *160 of 315*; 0.082 vs 0.084 labelled (combined vs N5-alone); grid-rail
  66% → **63%** (final admissible set); kernel-family 5× naive / **4.4×** under
  control; edge Jaccard **0.014–0.018**, 98.7% → **98.5%**; Ripley-K given at all
  three callers; composition surrogate given with its **0.42–1.83** spread; the
  two different COMMOT ρ's (0.78 mechanism-tile vs 0.90 all-interaction)
  distinguished; and the `\*` reimplementation marker is now on **every**
  CellChat\*/SpaTalk\*/NCEM\* mention, per CS_PHASE4 §5.1.
- ~~Figure 4 plotted call survival against CellWHISPER's FPR line unlabelled~~
  **Fixed.** The panel now states that our bars are call survival, and that
  survival equals an FPR only if the shuffle is a true null — the assumption
  under test. CS_PHASE4 §5 item 7 disclaims the conflation; the figure now
  carries the disclaimer too.
- ~~DeepScence and senepy lost with the container~~ **Reinstalled** at the pinned
  versions; every import in `code/` now resolves except `commot`, by design.

**Still outstanding — these need a human, not another agent:**

- ~~Two claims sit behind publisher paywalls~~ **The CellWHISPER one is now
  resolved (2026-08-27, citation audit).** The v1 full text was rendered from
  bioRxiv and the sentence reads verbatim: *"CellChat v2, COMMOT, and SpaTalk
  predict similar numbers of interactions on real (blue) and randomized (orange)
  data, indicating poor specificity and false positive rates (FPR) >90%. In
  contrast, CellWHISPER produced markedly fewer interactions on randomized data
  compared to real data, suggesting FPR < 5%."* Note the form: it is an
  **interaction-count ratio from which they infer an FPR**, not a measured type-I
  error — always write "implying". The same reading also confirmed that
  **CellWHISPER contains no torus shift** (the strings *torus*/*toroidal* do not
  occur) and uses no negative control probes; see `reports/CITATION_AUDIT.md`.
  **Still outstanding:** the **Ma et al. distance-gradient sentence**
  (cell.com/sciencedirect 403). The abstract-level claim — SSSs "serve as
  epicenters for heightened inflammation that compromises surrounding cells in a
  distance-dependent manner" — *is* verified from PubMed; the finer per-pathway
  monotone-with-distance claim is not, and is not presented here as verified.
- The **2026 Author Correction on Acosta et al. 2013** (doi
  10.1038/s41556-026-01959-z, *Nat Cell Biol* 28(6):1343) exists — re-confirmed
  via Europe PMC — but its **content is still unread**; nature.com redirects
  anonymous fetches to an auth endpoint and Europe PMC holds metadata only.
  Someone with Nature access must confirm it does not touch the transwell result.
- **M7 / M11 need a decision, not a lookup.** Cell counts differ by ≤1,100
  between reports (7259: 128,030 vs 127,386; 7250: 237,982 / 236,906 / 236,905),
  almost certainly pre- vs post-QC — state the QC rule once and use one number in
  Methods. And two framing sentences do rhetorical work with no citation
  (`BIO_PHASE2` §2.3 "that is textbook IFALD"; `BIO_PHASE3` §2.2 "consistent with
  the aged-liver ductular reaction literature") — the second carries the n = 1
  52-week sham animal and needs a reference if it reaches the Discussion.
- **Reference 20 is a meeting abstract**, not a paper (Neretti N, *Innov Aging*
  2024;8(Suppl 1):351 — ~250 words, single author, no methods, not peer
  reviewed), and it is the sole support for the primary/secondary senescence
  distinction underlying module **B7** `secondary_senescence`. B7 also shares 14
  of its 38 genes with the Tier A caller, so it is partly circular and
  `CS_PHASE3` §7 already recommends it never be the primary caller. Either
  replace the citation with a peer-reviewed source or drop the module's
  interpretive weight. Flagged in `references.bib` as do-not-cite.

All 30 references now live in [`references.bib`](references.bib) in the audit's
corrected form, with 25 inline `% AUDIT:` flags carrying each finding into the
write-up.

## Known limitations

- Mouse surgical IFALD model, not human aging. The framing changed with the
  dataset and the write-up must reflect that.
- 11 animals, of which only **6 sections are admissible** under Section 8
  Test 3 — **2 SBR and 4 sham**. An early reading that the sham arm fails Test 3
  came from one shallow section (7250, 0.48%); across all 11, four of five sham
  sections clear the floor while four of six SBR sections exceed the 20% ceiling.
- Erythroid, mast and lymphatic-endothelial cells are **not callable** on this
  panel and were dropped rather than assigned by noise.
- `Krt7` and `Cftr` do not detect, so the biliary compartment is reported as
  `Biliary/ductular` rather than as cholangiocytes.
- Enforcing Tier A ∩ Tier B = ∅ removes `Cdkn1a`, `Cdkn2a`, `Trp53`, `Lmnb1`
  and `Mki67` from the sender set. Headline results use per-module sender sets;
  cross-module comparisons use the union-strict 25-gene set.

### Environment drift — flagged honestly, now settled empirically

`numpy` was upgraded 1.26.3 → 2.4.6 partway through the project when an agent
pulled a transitive dependency, so Phase 1 ran under a different numpy than
later phases. This was recorded as unverified. On **2026-08-21** the working
container was reset and the whole scientific stack had to be reinstalled from
`requirements.txt`, which made the check free to run. It was run, at both ends:

**Phase 1 — the synthetic sweep.** Re-run in full under numpy 2.4.6 and diffed
against the stored results: **all 141,160 numeric values across all 43 per-config
files are bit-identical**, with identical shapes and column sets
(`results/repro_2026-08-21/`, `code/_repro_check.py`). Seeded reproduction holds
across the version boundary.

Two things the check turned up that are worth stating rather than hiding:

- **Column order** differs between runs in 35 of 43 files. The frames are
  assembled from per-run dicts, so order is not stable. Values are unaffected;
  anything consuming these CSVs must select by name, not position.
- The aggregate `results/sweep_all.csv` agrees with the re-run to ≤ 3×10⁻¹⁶
  relative (1 ulp) in 0.82 % of values rather than exactly. This is **not** a
  numpy effect: the stored aggregate differs from *its own* per-config
  checkpoints by exactly the same amount, an artefact of the original sweep
  having crashed and been resumed, so part of it was assembled in memory and
  part re-read from CSV. No reported figure moves at 1 ulp; the per-config files
  are the primary artefact and they are exact.

**Phase 4 — COMMOT, run as released software.** `commot` 0.0.3 cannot even be
*imported* under the pinned numpy 2.4.6 (`rho = np.Inf` as a default argument;
`np.Inf` was removed in NumPy 2.0), which means the original Phase 4 COMMOT
results were necessarily produced under numpy 1.x. `code/_shims/np2_compat.py`
restores the removed aliases — pure aliases, so nothing numerical changes — and
a stored checkpoint was then recomputed under 2.4.6: **81/81 rows, every p-value
bit-identical, max |diff| on the score column 3.4×10⁻²¹.**

`requirements.txt` now pins the verified state and documents both checks.
