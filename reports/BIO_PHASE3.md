# BIO PHASE 3 — All Sections, Composition Validation, Tier C, Caller Forensics

**Biology collaborator · 2026-08-20 · SASP Spatial Response Kernel**
Tasks T1 (annotate remaining sections), T2 (composition across all sections),
T3 (Deliverable 5 — ligand–receptor plausibility, §12), T4 (characterise caller disagreement, §10/§9-A7).

---

## 0. Headline

1. **All 11 liver sections of GSE310392 are annotated and carry the full five-file per-cell
   deliverable.** The animal-level bootstrap of §24.1 is unblocked: n = 11 mice, one section each,
   6 SBR / 5 sham, 1–2 mice per arm per timepoint.
2. **My Phase 2 headline result did not survive the full cohort, and I am withdrawing half of it.**
   The arm contrast is real, large and replicated. The *monotone SBR time course* and the *flat
   sham baseline* were both artefacts of which four sections happened to be on disk. §2.
3. **Deliverable 5 verdict: no fitted λ ordering across Tier C ligands can be interpreted as
   diffusion range on this panel.** The distance regressor is a detection-rate readout —
   log(median distance to nearest ligand⁺ cell) vs log(ligand⁺ cell density) gives
   **r² = 0.987, slope −0.54 against a Poisson prediction of −0.50**, across 6 SBR sections × 14
   ligands. The Section 9 internal control cannot separate diffusion from probe detection. §3.
4. **The ~99% caller disagreement is now explained mechanistically, and it is publishable.**
   Within cell type, the four sender definitions load on *opposite ends of the sequencing-depth
   distribution*: SenePy and `Cdkn1a`⁺ call the top depth quintile (2.3–2.5× enriched), Tier A and
   DeepScence call the bottom (1.3–1.9×). But conditioning on cell type **and** depth decile does
   **not** make them agree (Jaccard ratios 0.93–1.22× of chance). Depth explains *where* each
   caller looks; it does not explain the disagreement. §4.
5. **New confound the CS lead must handle:** within the SBR arm, *section-level* `Cdkn1a`⁺
   prevalence tracks *section-level* median transcripts per cell at Spearman **ρ = +0.94, p = 0.005**
   (n = 6). Animal-level "senescent burden" is not separable from animal-level detection depth
   without an explicit covariate. §5.

---

## 1. T1 — all 11 sections annotated

| Section | Mouse | Arm | Week | Cells (QC-pass) | Analysable | Unknown | median conf. |
|---|---|---|---|---|---|---|---|
| 7361 | 7361 | SBR | 2 | 193,983 | 161,579 | 7.9% | 0.604 |
| 7448 | 7448 | SBR | 10 | 187,536 | 162,303 | **8.4%** | 0.434 |
| 7450 | 7450 | SBR | 10 | 93,197 | 81,428 | 0.0% | — |
| 7259 | 7259 | SBR | 26 | 127,386 | 114,721 | 0.9% | 0.838 |
| 7260 | 7260 | SBR | 26 | 202,016 | 180,669 | 3.6% | 0.585 |
| 7239 | 7239 | SBR | 52 | 83,392 | 76,229 | 0.0% | 0.856 |
| 7352 | 7352 | sham | 2 | 139,378 | 122,380 | 1.8% | 0.940 |
| 7435 | 7435 | sham | 10 | 172,218 | 151,291 | 0.0% | — |
| 7248 | 7248 | sham | 26 | 224,921 | 207,785 | 0.0% | 0.849 |
| 7250 | 7250 | sham | 26 | 236,905 | 218,273 | 0.0% | 0.966 |
| 7001 | 7001 | sham | 52 | 165,961 | 159,279 | 0.0% | 0.590 |

The 12th GEO sample (`7239_tumor_sbr_Male_52-U1`) is a **tumour** section from the same mouse as
7239 and was deliberately not downloaded — different tissue, out of scope, and it would break the
one-section-per-animal structure.

Each section has `celltypes_`, `anatomy_`, `senders_`, `modules_`, `test3_prevalence_` CSVs in
`/workspace/data/processed/`, keyed on `cell_id`, identical schema to Phase 2. Backward-compatible
`sbr`/`sham` aliases for 7259/7250 are retained as symlinks.

### 1.1 One pipeline change: a merged label set (`cell_type_merged`)

7001 exposed a defect that would have silently corrupted T2. Its clusters 10 and 11 (19,420 cells,
11.7% of the section) score `Hepatic stellate cells` 2.65 / `Portal fibroblasts` 2.64 — a margin of
0.02 against a `MIN_MARGIN` of 0.20 — so both fell to `Unknown` and 7001 came out with **0.0%
stellate cells** while every other section has 7–17%. That is an artefact, not biology.

Three label families are not stably separable on this panel and the winner flips between sections:

| Family | Members | Evidence of instability |
|---|---|---|
| Mesenchymal | stellate / portal fibroblast / pericyte | all `Col1a1`/`Col1a2`/`Pdgfrb`⁺; 7001 margin 0.02 |
| Macrophages | Kupffer / inflammatory macs | 7448 splits 0.2% / 10.6%; 7450 splits 14.1% / 0.0% |
| Endothelial | LSEC / central-venous LSEC / portal endothelial | margins 0.6–1.4, portal endothelial called in only 1 of 11 sections |

`annotate_pipeline.py` now writes **both**: `cell_type` (fine, unchanged) and `cell_type_merged`,
assigned as the **max over each group's member scores** — "these are one compartment, take whichever
member evidences it". A union-of-markers mean was tried first and failed: it dilutes 2.65 to ≈2.0
and loses to a spurious `Mesothelial` runner-up at 1.92. Re-running is cheap (`--relabel`, ~90 s per
section, Leiden reused). **Use `cell_type_merged` for anything compared across sections; use
`cell_type` where the fine distinction is well separated within one section.**

After the merge, `Unknown` is 0.0% in 9 of 11 sections. The two remaining are real ambiguity, not a
bug, and I have left them as `Unknown` rather than force a call:

* **7448 cluster 4 (14,262 cells, 7.6%)** scores `Biliary/ductular` 0.93 vs `Hepatocytes` 0.82.
  This is precisely the ductular-metaplasia intermediate the Phase 2 `Biliary/ductular` naming
  decision was about — hepatocytes part-way through the transition. Assigning it to either side
  would be a fabrication.
* **7361 (7.9%)** is a mixture of the same and a DC/B/T-ambiguous cluster.

### 1.2 7001 checked by hand, because it is the single 52-week sham

7001 reports 14.6% `Biliary/ductular` where the 2–26 wk shams have 2.9–4.3%, which drives the whole
"sham is not flat" conclusion in §2. I verified it is not an annotation failure. Per-cluster mean
log-normalised expression, 7001 clusters 8+9 (23,191 cells):

`Epcam` 1.53/1.84, `Sox9` 0.81/0.71, `Pkhd1` 1.12/1.08, `Spp1` 3.29/2.93 — and
`Msln` 0.00–0.02, `Upk3b` 0.00, `Wt1` 0.01–0.02, so **not** mesothelium (mesothelial was only the
runner-up label). Section-wide, 7001 vs the 26 wk sham 7250: `Epcam` 0.268 vs 0.024 (11×),
`Col1a1` 0.209 vs 0.012 (17×), `Spp1` 0.738 vs 0.086 (9×). The biliary/mesenchymal expansion in the
aged control liver is in the counts, not in the classifier. **It rests on one animal and one
section and must be labelled as such.**

---

## 2. T2 — composition across all 11 sections

Full table: `/workspace/results/composition_by_arm_timepoint.csv`
Figure: `/workspace/figures/fig_phase3_composition.png`
Per-cell-type % is of **analysable** cells (excluding `Low_quality` and `Unknown`), because those
two labels swing 0–13% between sections for reasons unrelated to biology. Both denominators are in
the CSV.

| Section | Arm | Wk | Hepatocytes | Biliary/ductular | Mesenchymal | Macrophages | Endothelial | T/NK |
|---|---|---|---|---|---|---|---|---|
| 7361 | SBR | 2 | 50.39 | 2.71 | 7.27 | 16.81 | 14.67 | 2.31 |
| 7448 | SBR | 10 | 24.77 | 20.73 | 13.46 | 10.80 | 13.25 | 3.85 |
| 7450 | SBR | 10 | 38.24 | 10.01 | 10.98 | 14.13 | 16.03 | 3.41 |
| 7259 | SBR | 26 | 26.24 | 28.60 | 13.18 | 11.52 | 11.30 | 5.50 |
| 7260 | SBR | 26 | 31.03 | 16.12 | 16.91 | 10.58 | 13.08 | 3.13 |
| 7239 | SBR | 52 | 35.73 | 21.95 | 11.38 | 10.05 | 10.89 | 2.97 |
| 7352 | sham | 2 | 60.85 | 3.36 | 6.95 | 7.36 | 13.11 | 2.51 |
| 7435 | sham | 10 | 58.77 | 2.94 | 7.70 | 8.47 | 16.60 | 1.79 |
| 7248 | sham | 26 | 51.83 | 4.05 | 13.53 | 6.92 | 16.77 | 2.50 |
| 7250 | sham | 26 | 55.69 | 3.19 | 7.75 | 8.79 | 18.59 | 1.81 |
| 7001 | sham | 52 | 38.89 | 14.56 | 12.19 | 11.95 | 12.76 | 3.26 |

### 2.1 What holds

**The arm contrast is real, large and replicated at every timepoint from 10 weeks on.**
`Biliary/ductular` is the discriminating compartment:

| Week | SBR | sham | separation |
|---|---|---|---|
| 2 | 2.71 | 3.36 | none |
| 10 | 20.73, 10.01 | 2.94 | 3.4–7.1× |
| 26 | 28.60, 16.12 | 4.05, 3.19 | 4.0–9.0× |
| 52 | 21.95 | 14.56 | 1.5× |

Every SBR section from 10 weeks on (5/5) is above every 2–26 wk sham section (4/4), with no overlap.
Hepatocyte fraction separates the same way from 10 weeks on (SBR 24.8–38.2% vs sham 51.8–58.8%); at
2 weeks the SBR section is still at 50.4% and overlaps the sham range, as it biologically should.
**The ductular reaction is the strongest, most reproducible feature of this dataset**, and it is
recovered by a pipeline that knows nothing about arm or timepoint.

Also stable, and useful as a technical control: `Endothelial` sits at 10.9–18.6% in all 11
sections regardless of arm or timepoint (CV 17%). The annotation is not simply tracking overall
tissue quality.

### 2.2 What breaks — and this corrects my Phase 2 report

**The monotone SBR time course does not survive.** Phase 2 reported hepatocytes falling
42.0 → 33.4 → 23.6% across 2/10/26 wk. With both replicates at 10 and 26 wk:

* Hepatocytes: Spearman ρ vs timepoint = **−0.35, p = 0.49**. Not monotone.
* `Biliary/ductular`: ρ = +0.77, p = 0.076 — a large jump from week 2 to week 10, then flat/noisy.
* Only `Macrophages` (ρ = −0.85, p = 0.031) and `Endothelial` (ρ = −0.88, p = 0.020) are nominally
  monotone, and with n = 6 sections, 10 compartments tested and no multiplicity correction those
  p-values should not be believed.

**Between-replicate spread at the same timepoint is as large as the between-timepoint spread.**
At 10 wk SBR, hepatocytes are 24.8% (7448) vs 38.2% (7450); `Biliary/ductular` 20.7% vs 10.0%. At
26 wk SBR, 26.2% vs 31.0% and 28.6% vs 16.1%. The Phase 2 "monotone decline" was one section per
timepoint drawing a line through single points, and the 2 wk section (7361) — which is close to
sham, as it biologically should be two weeks after surgery — was doing most of the work.

**Sham is not flat either.** Hepatocytes 60.9 → 58.8 → 51.8 / 55.7 → 38.9 (ρ = −0.98, p = 0.005),
and `Biliary/ductular` 3.4 → 2.9 → 4.0 / 3.2 → 14.6. **The 52-week sham looks partly like an
injured liver.** Weeks 2–26 of sham *are* flat (hepatocytes 51.8–60.9, biliary 2.9–4.3); the change
is entirely the 52 wk animal. Read biologically that is an ageing effect (these mice are ~60 weeks
old at sacrifice) and it is consistent with the aged-liver ductular reaction literature; read
statistically it is **n = 1 animal**, and the SBR and sham arms *converge* at 52 weeks
(hepatocytes 35.7 vs 38.9; biliary 22.0 vs 14.6).

### 2.3 The honest statement for the paper

> Across 11 sections from 11 mice, cell-type composition separates the two surgical arms cleanly
> from 10 weeks onward — the ductular-reaction compartment is 3–9× larger in SBR than in
> timepoint-matched sham with no overlap — but it does not resolve a within-arm time course. With
> one to two animals per arm per timepoint, between-animal variation at a fixed timepoint is as
> large as the variation across timepoints, and the longest-surviving control animal shows an
> injury-like phenotype of its own.

**Consequence for §24.1.** Treat *arm* as the contrast with animal as the unit (n = 6 vs 5) and
report timepoint as a covariate or a stratification, not as an ordered dose. Do not bootstrap a
per-timepoint λ: at 2 wk and 52 wk that is n = 1 per arm.

---

## 3. T3 — Deliverable 5: Tier C ligand–receptor plausibility

Outputs: `/workspace/results/phase3/tierC_expression_by_celltype.csv`,
`tierC_ligand_identifiability.csv`, `tierC_meta.json`;
figure `/workspace/figures/fig_phase3_tierC_identifiability.png`.
Computed on all 6 SBR sections plus the 26 wk sham; per-cell-type numbers below are SBR 7259.

### 3.1 Receptor side — good coverage, one important exception

Detection rate (% of cells positive), SBR 7259, selected receivers:

| Receptor | Hepatocytes | Kupffer | LSEC | Stellate | T/NK | DC | Biliary/duct. |
|---|---|---|---|---|---|---|---|
| `Tnfrsf1a` | **54.9** | 40.1 | 45.9 | 39.3 | 27.7 | 34.1 | 25.2 |
| `Il6st` (gp130) | **33.6** | 18.3 | **63.6** | 36.0 | 16.3 | 21.2 | 21.5 |
| `Il1r1` | **26.7** | 6.3 | 16.3 | 20.3 | 7.2 | 5.9 | 8.6 |
| `Tgfbr2` | 24.9 | 27.8 | 26.0 | 27.3 | 28.8 | 22.8 | 11.6 |
| `Il6ra` | 25.3 | 19.8 | 4.6 | 8.3 | 6.8 | 17.5 | 5.8 |
| `Tnfrsf1b` | 26.4 | **38.5** | 7.7 | 26.7 | 25.4 | 25.4 | 8.3 |
| `Cxcr4` | 7.7 | 29.0 | 3.8 | 7.1 | 23.7 | 30.7 | 3.9 |
| `Dpp4` | 13.5 | 3.3 | **19.7** | 4.1 | 3.8 | 15.5 | 3.5 |
| `Ccr2` | 0.7 | 8.2 | 0.5 | 1.4 | 6.4 | **16.2** | 0.8 |
| **`Ackr3` (CXCR7)** | **0.2** | **0.3** | **1.0** | **0.7** | **0.3** | **0.3** | **0.2** |

**Correction to my Phase 1 report.** I wrote that the CCL2/CXCR7/DPP4 axis was "fully covered"
because `Ccr2`, `Ackr3`, `Cxcr4` and `Dpp4` are all on the panel. Being on the panel is not the
same as being detected. **`Ackr3` is detected in ≤1.0% of cells in every cell type of every
section.** CCL2→CCR2 is testable (DC 16.2%, Kupffer 8.2%, T/NK 6.4%); **CCL2→ACKR3 and
CXCL12→ACKR3 are not.** The Ackr3 arm of the 2026 brain paracrine-senescence result cannot be
replicated here, and Deliverable 5 should say so rather than list Ackr3 as covered.
`Cxcl12`→`Cxcr4` and the `Dpp4` range-limiting hypothesis remain testable.

### 3.2 Ligand side — three canonical SASP ligands are effectively absent

% of analysable cells positive, SBR 7259 / sham 7250:

| Ligand | SBR | sham | ligand⁺ **and** `Cdkn1a`⁺ cells in the whole SBR section |
|---|---|---|---|
| `Cxcl12` | 49.9 | 63.9 | 3,161 |
| `Igfbp3` | 16.7 | 11.0 | 757 |
| `Tgfb1` | 6.8 | 7.1 | 386 |
| `Cxcl5` | 6.6 | 0.17 | 283 |
| `Cxcl1` | 5.8 | 0.83 | 474 |
| `Thbs1` | 5.7 | 1.74 | 263 |
| `Ccl2` | 4.6 | 0.98 | 274 |
| `Il1a` | 2.7 | 4.6 | 220 |
| `Timp1` | 2.1 | 0.17 | 110 |
| `Tnf` | 1.1 | 0.21 | 52 |
| `Gdf15` | 0.44 | 0.43 | 114 |
| **`Il6`** | **0.087** | 0.096 | **3** |
| **`Cxcl2`** | **0.083** | 0.023 | **8** |
| **`Mmp3`** | **0.071** | 0.038 | **3** |

`Il6` — the canonical SASP ligand, and the one that motivates Tier B module B2 — is detected in
**3 senescent cells in a 127,386-cell section**. `Cxcl2` and `Mmp3` likewise. No kernel can be fit
from a ligand's own expression for these three; the B2 IL6/JAK/STAT3 readout can still be scored in
receivers, but **it cannot be tied back to an Il6-expressing sender in this data**.

At the other extreme `Cxcl12` is expressed by 50–64% of all cells (constitutive in liver
endothelium). Distance to the nearest `Cxcl12`⁺ cell has a median of 5–8 µm — below the
resolution floor (median NN cell distance 6.7 µm). It has no dynamic range.

Ligand enrichment in `Cdkn1a`⁺ senders over `Cdkn1a`⁻ cells **of the same cell type** is real but
modest: `Il1a` 1.5–2.9×, `Ccl2` 1.3–2.7×, `Cxcl1` 1.0–1.8×, `Tgfb1` 1.0–1.7×, `Igfbp3` ~1.0–1.3×.
Senders are enriched for SASP ligands; they are nowhere near being defined by them.

### 3.3 The internal control fails before any model is fitted

Section 9 asks for an ordering check: membrane-bound `Il1a` should fit the shortest λ, small
diffusible chemokines the longest. That check is only meaningful if the *regressor* — distance to
the nearest ligand-expressing sender — carries ligand-specific spatial information. It does not.

Across **6 SBR sections × 14 on-panel Tier C ligands** (84 points):

> **log₁₀(median distance to nearest ligand⁺ cell) = −0.54 · log₁₀(ligand⁺ cell density) + 2.65,
> r² = 0.987.** Poisson theory for randomly placed points predicts slope exactly −0.50.

The distance regressor is, to 98.7% of its variance, a deterministic function of how often the gene
is detected, and the fitted slope says the ligand⁺ cells are placed essentially at random at the
scale that sets that distance. Detection rate on Xenium is set by expression level and probe
efficiency — not by diffusion range.

The consequence is visible directly in the ranks. Median distance to the nearest ligand-expressing
sender, and the rank each ligand would imply (1 = shortest ⇒ shortest fitted λ):

| Section | `Il1a` | `Ccl2` | `Cxcl1` | `Cxcl5` | `Cxcl2` | `Il1a` rank |
|---|---|---|---|---|---|---|
| 7361 (2 wk) | 70.7 | 59.7 | 103.4 | 269.0 | 190.2 | 2 |
| 7448 (10 wk) | 90.6 | 112.7 | 110.8 | 182.8 | 772.5 | **1** |
| 7450 (10 wk) | 87.9 | 103.9 | 156.8 | 277.3 | 645.4 | **1** |
| 7259 (26 wk) | 137.8 | 123.9 | 94.9 | 113.3 | 697.3 | **4** |
| 7260 (26 wk) | 119.7 | 105.0 | 99.5 | 195.5 | 465.9 | 3 |
| 7239 (52 wk) | 73.3 | 90.5 | 70.7 | 102.2 | 391.5 | 2 |

`Il1a`'s rank moves between 1st and 4th across sections of the same arm. `Cxcl1` — predicted
*longest* — is 1st in two sections. And `Cxcl1` vs `Cxcl2`, two genes encoding functionally
interchangeable analogues of the same absent human ligand (CXCL8), differ by **4–9×** in every
section, purely because `Cxcl1` is detected in 5.8% of cells and `Cxcl2` in 0.083%.

### 3.4 Verdict (Deliverable 5, one paragraph)

**No.** On this panel and in this tissue, a fitted λ ordering across Tier C ligands would not be
biologically interpretable, and the Section 9 internal control cannot do the job it was designed
for. Three independent reasons, in order of severity. First, the regressor is a detection-rate
readout: distance to the nearest ligand-expressing sender is determined to r² = 0.987 by the
ligand's detection frequency, with a Poisson slope, so any cross-ligand λ ordering recapitulates
probe sensitivity rather than diffusion physics — and empirically the ordering is unstable, with
membrane-bound `Il1a` ranking anywhere from 1st to 4th across sections of the same arm, and two
functional analogues of the same ligand (`Cxcl1`, `Cxcl2`) separated by 4–9×. Second, the control's
sharpest form is unavailable: `Il1b` is off-panel, so membrane-bound versus secreted signalling
through the shared receptor `Il1r1` cannot be contrasted, and mouse has no CXCL8 orthologue, so the
"longest-range" pole is represented only by three analogues that disagree with each other. Third,
the ligands that matter most are not measurable: `Il6` is detected in 3 senescent cells per section,
`Cxcl2` in 8, `Mmp3` in 3, while `Cxcl12` is expressed by half of all cells and its
nearest-neighbour distance (5–8 µm) sits below the segmentation resolution floor of 6.7 µm; and
`Ackr3`, which Section 9 flags as highest priority, is detected in ≤1% of cells everywhere.
**What can still be defended** is a *within-ligand* statement — for a single ligand–receptor pair
with adequate detection (`Ccl2`→`Ccr2` into DC/Kupffer/T-NK, `Tnf`→`Tnfrsf1a/1b` into hepatocytes,
`Tgfb1`→`Tgfbr1/2`, `Il1a`→`Il1r1` into hepatocytes at 26.7%) — whether the response falls off with
distance, with the density of that ligand's own cells conditioned on, and reported as a single λ
with its confidence interval rather than as a rank against other ligands. The cross-ligand ordering
test should be reported as **attempted and failed for an identifiability reason**, which is squarely
the kind of result this paper exists to report.

---

## 4. T4 — what the callers are actually picking up

Outputs: `/workspace/results/phase3/caller_technical_loading.csv`, `caller_technical_loading2.csv`,
`caller_celltype_composition.csv`, `caller_strata.csv`, `caller_within_type_depth_bias.csv`,
`caller_pairwise_agreement.csv`, `caller_depth_stratified_agreement.csv`,
`caller_agreement_depth_and_type_matched.csv`; figure `fig_phase3_caller_depth.png`.

Note on thresholding: the stored `sender_flag_p*` columns are *within-cell-type* percentiles and are
therefore flat across cell types **by construction**. Every caller is re-thresholded here at a
common rule so the comparison is like-for-like; both global and within-type versions are reported.

### 4.1 The stable answer: each caller selects a different end of the depth distribution

Enrichment of the top-5% calls by transcript-count quintile, **computed within cell type** so cell
type cannot explain it (1.0 = no bias):

| Caller | Q1 low | Q2 | Q3 | Q4 | Q5 high | Q5/Q1 |
|---|---|---|---|---|---|---|
| **Tier A score** (sham / SBR) | 1.58 / 1.79 | 1.22 / 1.28 | 1.02 / 0.92 | 0.76 / 0.69 | 0.43 / 0.31 | **0.27 / 0.17** |
| **DeepScence** | 1.78 / 1.33 | 1.20 / 1.11 | 0.96 / 1.03 | 0.67 / 0.87 | 0.39 / 0.67 | **0.22 / 0.50** |
| **SenePy** | 0.12 / 0.09 | 0.36 / 0.39 | 0.72 / 0.78 | 1.27 / 1.20 | 2.53 / 2.54 | **21× / 27×** |
| **`Cdkn1a` > 0** | 0.38 / 0.28 | 0.77 / 0.55 | 1.02 / 0.75 | 1.23 / 1.17 | 1.61 / 2.26 | **4.2× / 8.2×** |

This is consistent in both arms and it is the cleanest single explanation of the disagreement:
**two callers select the bottom of the depth distribution and two select the top of it, inside every
cell type.** They are close to disjoint by construction.

`Cdkn1a` > 0 having a 4–8× depth gradient is not a bug in anyone's method — it is the fact that on
a counting assay, "gene X detected" is partly "this cell had more transcripts". Every binary marker
call on Xenium carries it.

### 4.2 Cell type: DeepScence and SenePy are each dominated by one artefact

Enrichment of global top-5% calls by cell type:

| | Hepatocytes | Kupffer | T/NK | B | Proliferating |
|---|---|---|---|---|---|
| DeepScence, sham | **0.03** | 4.62 | 2.35 | 2.28 | 6.22 |
| DeepScence, SBR | **2.46** | 0.08 | 0.14 | 0.16 | 2.20 |
| SenePy, sham | 0.26 | **9.50** | 0.00 | 0.20 | 0.00 |
| SenePy, SBR | 0.16 | **8.30** | 0.00 | 0.00 | 0.00 |
| Tier A, sham | 0.66 | 1.00 | 3.78 | 3.15 | 6.23 |
| `Cdkn1a`⁺, sham | 0.66 | 1.02 | 1.56 | 1.04 | **23.4** |

* **SenePy's Kupffer enrichment is a hub-size artefact, and it is measurable.**
  Per-cell-type mean SenePy score vs the number of that hub's genes on our panel:
  Kupffer 7.03 (226 genes), T/NK 1.20 (83), Hepatocytes 1.00 (54), LSEC 0.78 (62), B 1.35 (68) —
  **r = 0.992**. The hub score scales with hub size, so it is not comparable across cell types and a
  global threshold puts 83–95% of its calls into one type. **SenePy must be thresholded within cell
  type.** Separately, SenePy has no hub for `Biliary/ductular`, mesenchymal, vSMC, mesothelial,
  portal-endothelial or proliferating cells here — it cannot score 16–45% of a section at all.
* **DeepScence's cell-type profile is its depth profile, and its sign is not stable between
  sections.** In sham it calls small, low-count cells (ρ with transcript counts −0.35, with cell
  area −0.43) and therefore essentially no hepatocytes (0.03×) in a hepatocyte-dominated tissue. In
  SBR it does the opposite (ρ +0.32, +0.32) and calls hepatocytes at 2.46×.

### 4.3 Where the DeepScence instability comes from — and a caveat I own

To separate "the gene set loads on depth" from "the model loads on depth", I scored the same cells
with a **naive CoreScence score** (occurrence ≥ 5, up-genes minus down-genes, 17 up / 14 down on the
mouse panel), computed here:

| Score | ρ with transcript counts, sham | ρ, SBR |
|---|---|---|
| naive CoreScence (its own gene set) | **+0.070** | **+0.326** |
| DeepScence `ds` | **−0.350** | **+0.318** |

In SBR the two agree (+0.33 vs +0.32, Spearman between them +0.387). In sham DeepScence **inverts
relative to its own gene set**. DeepScence fixes the sign of its bottleneck node by correlating it
with `CDKN1A`; in 7250, `Cdkn1a` is detected in 0.48% of hepatocytes, which is a very weak anchor.
That is a plausible and, I think, likely mechanism — but I am labelling it an **inference**, not a
measurement, because I did not instrument the model internals.

**The caveat I have to state loudly:** we run DeepScence with `denoise=False`, a forced deviation
(its DCA dependency needs an obsolete TensorFlow stack), and DCA denoising is precisely the step
that would normalise depth. So §4.2–4.3 characterise *DeepScence as we could run it on this panel*,
not DeepScence as published. It should be reported in exactly those words.

### 4.4 But depth does not explain the disagreement

If depth were the whole story, matching on it should produce agreement. It does not. Top-5%
recomputed **within each cell type × within-type depth decile** (54–120 strata), Jaccard as a ratio
to chance:

| Pair | sham | SBR |
|---|---|---|
| Tier A vs SenePy | 0.93 | 0.94 |
| Tier A vs DeepScence | 1.04 | 1.22 |
| Tier A vs `Cdkn1a`⁺ | 0.96 | 1.06 |
| SenePy vs `Cdkn1a`⁺ | 1.37 | 0.98 |
| SenePy vs DeepScence | 2.15 | **0.38** |
| DeepScence vs `Cdkn1a`⁺ | 2.85 | 1.51 |

Four of six pairs sit at **0.93–1.22× of chance** — statistical independence — after conditioning on
both cell type and depth. Partial Spearman given log counts moves nothing (e.g. Tier A vs SenePy
−0.024 → −0.037). The one pair that looked concordant in sham (SenePy vs DeepScence, 2.15×,
ρ = +0.12) is **anti**-concordant in SBR (0.38×, ρ = −0.25). The only pair that is consistently
above chance is DeepScence vs `Cdkn1a`⁺, which is expected because `CDKN1A` is DeepScence's own
direction anchor.

### 4.5 The finding, stated for the paper (§9 A7, §10)

> Four published or panel-standard ways of calling a senescent cell — DeepScence, SenePy, a
> disjoint arrest-and-damage score, and `CDKN1A` positivity — were applied to the same 127k–237k
> cells. Their top-5% calls overlap at 0.93–1.22× of chance after conditioning on cell type and
> sequencing depth, i.e. they are statistically independent. What each one *does* select is
> identifiable and technical: within cell type, SenePy and `CDKN1A`⁺ are enriched 2.3–2.5× in the
> highest transcript-count quintile and depleted 8–11× in the lowest, while the arrest score and
> DeepScence run the other way; SenePy's cross-cell-type score scales with the number of its hub
> genes on the panel (r = 0.992) and so is not comparable across cell types; and DeepScence's
> correlation with sequencing depth reverses sign between two sections of the same study. A
> senescence call on targeted spatial data is therefore not a noisy measurement of one latent
> state — it is a choice of which end of the detection-depth distribution to name senescent, and
> the length constant that follows inherits that choice.

This is a headline-grade result in its own right and it is exactly what §9 A7 asked for, quantified.

---

## 5. New confound: section-level burden tracks section-level depth

`/workspace/results/section_qc_sender_summary.csv`

| Section | Arm | Wk | median transcripts/cell | `Cdkn1a`⁺ hepatocytes % |
|---|---|---|---|---|
| 7259 | SBR | 26 | 446 | 9.64 |
| 7260 | SBR | 26 | 596 | 10.46 |
| 7450 | SBR | 10 | 743 | 22.71 |
| 7448 | SBR | 10 | 754 | 25.61 |
| 7361 | SBR | 2 | 786 | 25.34 |
| 7239 | SBR | 52 | 824 | **45.04** |
| 7250 | sham | 26 | 511 | **0.48** |
| 7352 | sham | 2 | 646 | 7.20 |
| 7435 | sham | 10 | 766 | 2.31 |
| 7001 | sham | 52 | 840 | 8.94 |
| 7248 | sham | 26 | 968 | 4.92 |

Within the SBR arm, Spearman ρ (median transcripts/cell, `Cdkn1a`⁺ hepatocyte %) = **+0.943,
p = 0.005** (n = 6). Within sham, ρ = +0.50, p = 0.39 (n = 5). Pooled ρ = +0.16 (arm dominates).

I am **not** claiming this is purely technical — deeper sections may also be less injured, with
larger, more transcript-rich hepatocytes, so depth and phenotype are entangled. But the practical
consequence is unambiguous: **animal-level senescent burden cannot be compared across sections
without median transcripts per cell as a covariate.** Add it to Tier D at the section level, not
only the cell level.

Two further corrections to Phase 2 that follow from the full cohort:

* **§8 Test 3 does not fail in sham.** Phase 2 said sham hepatocytes are at 0.48%, below the 1%
  floor, and that estimation must therefore happen in SBR. With all five sham sections the
  hepatocyte values are 7.20 / 2.31 / 4.92 / **0.48** / 8.94 — **4 of 5 pass**, and 7250 (the one
  section Phase 2 had) is the outlier and also the shallowest sham section. Sham can be an
  estimation arm for 4 of 5 animals. The SBR-arm preference still stands on other grounds (larger
  dynamic range, more cell types in band), but the stated reason was wrong.
* **7239 (SBR, 52 wk) is at 45.0% `Cdkn1a`⁺ hepatocytes — above §8's 20% ceiling.** In that section
  distance-to-nearest-sender is near zero everywhere and λ is unidentifiable by construction.
  Exclude 7239 from the hepatocyte kernel fit, or threshold it higher and say so.

---

## 6. Anatomy: the SBR portal-triad failure replicates in all 5 SBR sections

Phase 2 flagged `dist_to_portal_triad_um` as invalid in SBR from one section. Now measured in all 11:

| Section | Arm | Wk | biliary foci | corr(zonation, dist-to-portal-triad) in hepatocytes |
|---|---|---|---|---|
| 7352 | sham | 2 | 115 | **+0.265** |
| 7435 | sham | 10 | 140 | **+0.247** |
| 7248 | sham | 26 | 226 | **+0.343** |
| 7250 | sham | 26 | 213 | **+0.255** |
| 7001 | sham | 52 | 439 | **+0.318** |
| 7361 | SBR | 2 | 118 | +0.198 |
| 7450 | SBR | 10 | 152 | +0.121 |
| 7448 | SBR | 10 | 715 | −0.078 |
| 7259 | SBR | 26 | 725 | +0.000 |
| 7260 | SBR | 26 | 444 | +0.005 |
| 7239 | SBR | 52 | 330 | −0.085 |

**All 5 sham sections give +0.25 to +0.34 with 115–439 foci. All SBR sections from 10 weeks on give
−0.09 to +0.12.** The only SBR section that behaves is 7361, the 2-week section that is also the one
without a ductular reaction. The landmark collapses exactly when and where the ductular reaction
appears. Phase 2's n = 1 warning is now a fully replicated result, and it is a good figure panel:
one architectural covariate that is valid in controls and destroyed by the disease being studied.
**Use `zonation_score`, per §11, for anything in the SBR arm.**

---

## 7. Deliverables produced in Phase 3

| Path | Contents |
|---|---|
| `data/processed/{celltypes,anatomy,senders,modules,test3_prevalence}_*.csv` | **all 11 sections**; `celltypes_*` now carries `cell_type_merged` + `cell_type_merged_confidence` |
| `results/composition_by_arm_timepoint.csv` | T2 table, merged + fine labels, both denominators |
| `results/section_qc_sender_summary.csv` | per-section depth, burden, portal-triad validity |
| `results/phase3/tierC_expression_by_celltype.csv` | Tier C ligand/receptor detection + sender enrichment, per cell type, 7 sections |
| `results/phase3/tierC_ligand_identifiability.csv` | ligand⁺ density, distance distributions, per section |
| `results/phase3/caller_*.csv` (8 files) | T4 forensics |
| `figures/fig_phase3_composition.png` | T2 |
| `figures/fig_phase3_caller_depth.png` | T4 |
| `figures/fig_phase3_tierC_identifiability.png` | T3 |
| `code/{tierC_lr,caller_disagree,caller_disagree2,composition_all,make_phase3_figs}.py` | new; `annotate_pipeline.py`, `sasp_io.py` patched |

## 8. What the CS lead must know

1. **All 11 sections are ready.** Re-run `prepare_samples.py --force`; nine new sections are now
   available to the fit grid and the animal-level bootstrap.
2. **Use `cell_type_merged` for anything cross-section** (§1.1). The fine labels flip between
   sections for three families and will manufacture spurious between-animal variance.
3. **Withdraw the Phase 2 monotone-progression claim** (§2.2). Arm contrast: yes, strong. Time
   course: not resolvable at n = 1–2 per timepoint. Do not write it into the abstract.
4. **Deliverable 5 verdict is negative and that is the result** (§3.4). Do not fit a cross-ligand λ
   ordering as a validation. Fit within-ligand kernels for the four pairs with adequate detection,
   and report the identifiability failure as a finding.
5. **Add section-level median transcripts per cell to Tier D** (§5). Within SBR it correlates with
   burden at ρ = +0.94.
6. **Exclude 7239 from hepatocyte kernel fits** (45% `Cdkn1a`⁺, above §8's ceiling) — or state the
   exclusion.
7. **Correction: sham does not fail Test 3** (§5); 7250 was an unrepresentative section.
8. **Threshold SenePy within cell type only** (§4.2); its raw score scales with hub size at r = 0.992
   and it cannot score 16–45% of a section at all.
9. **`dist_to_portal_triad_um` is invalid in every SBR section from 10 weeks on** — now replicated
   5/5 (§6). Sham is valid 5/5.
10. **The caller forensics (§4) are, in my view, the strongest single result Phase 3 produced** and
    they serve the paper's thesis directly. They belong in a main figure, not the supplement.

### Open / not done
* Transcript assignment rate (§8 Test 1) still unmeasured — needs one `transcripts.parquet`.
* DeepScence scores exist for the two 26 wk sections only; the other nine would take ~50 min each.
  Given §4.3 I would spend that compute elsewhere unless the CS lead wants the second caller on
  every animal.
* B7 split-half `A_call`/`A_readout` design still not built (carried over from Phase 1).
* `naive_corescence` is computed inside `caller_disagree2.py` but not written per cell; say the word
  and it becomes a column in `senders_*.csv`.
