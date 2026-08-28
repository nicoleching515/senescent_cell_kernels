# How Far Does Senescence Actually Reach?

Senescent cells secrete inflammatory factors that are believed to push neighbouring
cells toward senescence — a local contagion. The field assumes this and has never
measured its range. This repository asks whether the range is **estimable at all**
from the data that exists: a single spatial snapshot per tissue, with every cell's
position and expression, no time course and no perturbation.

The answer is that it is not, at achievable power — and the interesting part is
what we found while establishing that.

---

## Methods

**The estimand.** A distance-decay kernel: response amplitude in a receiver cell as
a function of its distance to the nearest senescent ("sender") cell, with a length
constant λ. Five kernel families, selected by strict argmin AIC.

**The problem.** Senders cluster; cell types sit near their own kind; sequencing
depth varies across a section. All three produce a distance-dependent response
signal with no signalling in it. So the analysis is a confounding problem, and it
sits inside the **spatial confounding** literature (Hodges & Reich 2010; Dupont
*et al.*, *Spatial+*, 2022) — one of our nuisance adjustments is Spatial+'s
residualise-the-covariate move, arrived at independently.

**The null battery, N1–N8.** Each removes one candidate explanation:

| | |
|---|---|
| N1 | within-cell-type label permutation |
| N2 | matched-decoy contrast |
| N3 / N4 | coordinate shift / rotation of the sender set |
| N5 / N6 | technical-covariate and neighbourhood-baseline adjustment |
| N7 | sender-definition axis |
| N8 | gene-set disjointness |

**The estimator's calibration is measured, not assumed.** The same kernel is fitted
to the assay's own negative-control features — probes and codewords with no biology
by construction. This is a negative-control *outcome* in the sense of Lipsitch
*et al.* (2010), not a generic QC statistic: the estimand's own estimator is pointed
at a target whose true value is zero.

**Sender definitions.** Four independent callers — a curated arrest-and-damage score,
DeepScence, SenePy, and `CDKN1A`⁺ positivity — with Tier A (sender-defining) held
strictly disjoint from the seven Tier B response modules, enforced by a gate that
exits non-zero.

**Pre-registration.** Every parameter was frozen and committed before the human data
was analysed (see *Why we froze*, below). Deviations are declared, never silent:
`reports/PREREG_PHASE8.md` carries eleven dated corrections and a deviation table.

---

## Data

Two arms, both 10x Xenium Prime 5K in situ spatial transcriptomics.

| | Mouse (M1) | Human (H1) |
|---|---|---|
| Accession | GEO `GSE310392` | GEO `GSE326743` |
| Tissue | Liver, sham vs SBR fibrosis | Spleen, normal |
| Sections | 11 | 7 |
| Cells | 1,834,806 (1,036,459 admissible) | 2,207,593 |
| Panel | 5,097 genes (5,006 stock + 91 custom) | 5,093 genes |
| Controls | 609 codewords / 40 probes / 21 genomic | identical |

H1 was selected by a panel-first screen of all **132 series** on GEO platform
`GPL33762`. Title-only screening finds 5 candidates; full-metadata screening finds
19; only one is normal tissue with ≥3 donors. Every panel was verified **on the
data**, not from the series description. The ortholog-intersected panel used for
cross-arm comparison is **2,425 genes**.

Raw bundles are not tracked. `code/fetch_h1_geo.sh` and `code/fetch_xenium_bundle.sh`
re-fetch them; `reports/PHASE7_H1_SCREEN.md` records the screen and its audit trail.

---

## Significant findings

### 1. The kernel is not estimable at achievable power — and this replicates

Controlled amplitude **0.0288** response-SD on the mouse arm against its own
80 %-power detectable bound of **0.1833**; **0.0307** against **0.1094** on the human
arm. **Zero of seven response modules cross the bound on either arm**, on either
sender call, on either gene panel.

The fitted length constant is **λ̂ = 14.73 µm, IQR [7.0, 50.0], with 60 % of fits
pinned to a grid bound** — the distribution, not the point estimate, is the result.

The null replicates across species and tissue. **The confound structure does not**,
which means a confound characterisation published on one dataset cannot be carried
to another.

### 2. The field's standard control does not work — shown three independent ways

Matching senders to comparable non-senders is the usual defence against confounding.
It fails here, and the three demonstrations are independent:

- **Composition-matched decoys** removed **1.6 %** of the naive amplitude where the
  *same variables entered as covariates* removed **85.4 %** — a factor of 52.
- **On the assay's own negative controls**, matched decoys left **86 %** of a purely
  technical distance gradient intact; covariate adjustment removed it.
- **On synthetic data with a planted kernel**, matched decoys gave *worse* confidence-
  interval coverage than no correction at all: **0.35**, against **0.51** naive and
  **0.85** nuisance-conditioned, and worse bias than naive in 8 of 20 regimes.

The third is decisive, because with a known planted answer the objection *"your
covariates removed real signal"* does not apply. It had been sitting in the project's
own headline synthetic figure since the first phase, unread.

### 3. A measured false-positive rate

**9–16 % against a nominal 5 %**, obtained from assay-internal null features — no
simulation, no held-out data, no randomisation.

### 4. The standard spatial null is undefined on real tissue

Toroidal shift nulls require a rectangular window (**Lotwick & Silverman, 1982**).
Tissue sections are not rectangles: **23 %** of shifted senders land outside the
tissue. Tiling the shift — the intuitive fix — makes it *worse*: measured type-I
error up to **2.35× nominal** on an irregular window, exactly as Mrkvička *et al.*
(2021) predicts for a union-of-rectangles window. The published variance correction
holds nominal and is now the primary null.

The finding is not that this is unknown — it is forty years old — but that it is
being violated in current spatial-omics practice, and here is the price.

### 5. Senescence callers are not interchangeable

Four callers on the same cells agree at **1.21× chance** after conditioning on cell
type and sequencing depth — weakly dependent, not independent. Each selects a
different end of the detection-depth distribution. One caller's polarity is anchored
on a single gene; another is not seed-reproducible on human data (top-5 % Jaccard
**0.157–0.698** between seeds, and one of five seeds returned an inverted score).

Sender prevalence tracks **sequencing depth, not age**: ρ = **0.82** in human spleen
(against ρ = −0.04 for donor age), reproducing ρ = 0.94 in mouse liver.

---

## Thought process

**Start from identifiability, not effect size.** The first question was not "how big
is the kernel" but "under what conditions could this kernel be recovered at all."
Figure 1 is a synthetic regime map with a planted ground truth, built before any
real data was fitted. It is also where finding 2 was hiding.

**Report against interest, by default.** Every correction applied to this project
moved *against* its own result, and none changed the conclusion. That is the strongest
robustness statement available, and it is deliberate: the winner's-curse correction
went the wrong way, zonation failed to be the confound the plan predicted, and
`denoise=False` was declared rather than buried.

**Freeze before looking.** See below.

**Numbers come from files, not from documents.** This was learned the hard way. A
length constant quoted in six places turned out to be derived from the very claim it
supported; a circularity denominator existed in no file; a bibliography was written
from recall and had 41 wrong author names. Every headline number is now re-derivable
by a committed producer, and `code/_repro_artefacts.sh --check` verifies 15 of them
against their committed copies.

**Make the failure modes loud.** Several scripts silently destroyed their own outputs
and exited 0 — an empty section list wrote a 1-byte `main_fits.csv`; a superseded
figure producer replaced a committed figure while its data CSVs came back identical.
Those paths now refuse. The prohibitions in the pre-registration are mechanically
checked rather than trusted.

**Prefer an honest gap to a plausible reconstruction.** Where a producer could not be
faithfully rebuilt, none was written and the gap is recorded at the read site.

---

## Why we froze

The mouse arm was analysed first and its results were known. Everything decided after
that point — thresholds, kernel families, null definitions, gene-set membership,
sender callers — could have been chosen, consciously or not, to make the human arm
agree with it. A reader cannot distinguish a prediction from a rationalisation after
the fact, and neither can the person who made it.

So every parameter was fixed, committed, and tagged (`phase8-frozen`) **before the
human expression data was read**. Only panel membership was inspected beforehand,
which the protocol explicitly sanctions as a screening step.

It earned its keep immediately. On the human arm the sender-prevalence gate came
within one cell-count floor of failing, the fine/merged label interaction left cells
uncallable, and the anatomical covariate validated weakly. Each was a live temptation
to adjust a threshold. Because the parameters were frozen, each became a **declared
deviation with both versions reported** instead of a silent fix — and the pre-
registration's own replication criteria were found to be **non-exhaustive**, a defect
that is now on the record rather than resolved in whichever direction suited the
result.

The freeze also failed usefully in one respect: the first tag shipped code that could
not reproduce its own results, because working-tree edits had never been committed.
That is recorded in the re-cut tag's message rather than erased.

---

## Repository layout

```
code/        producers, guards, and verifiers
genesets/    frozen Tier A-E sets, both arms, with provenance pins
data/raw*/   Xenium bundles (untracked; see Reproducing the data)
results/     all analysis outputs
figures/     58 guarded artefacts + their *_data.csv
reports/     phase reports, audits, corrections, pre-registration
```

## Reproducing

```bash
. code/_env.sh                        # names the interpreter; fails loudly if absent
python3 code/check_figures_guard.py   # 58 artefacts, content-hashed
python3 code/gate_genesets_guard.py   # Tier A/B disjointness, both arms
python3 code/phase10_verify_report.py # 73 headline checks
python3 code/h1_verify_report.py      # 28 headline checks
bash    code/_repro_artefacts.sh --check   # 15 artefacts vs committed copies
python3 code/check_prohibitions.py <paths> # section-10 prohibitions
```

## Known limitations

- The `cache3` → `main_fits.csv` step is verified on sampled jobs, not end to end.
- Two cited inputs have no producer and are documented rather than reconstructed:
  `zonation_gene_correlations_7250_sham.csv` (a gene-set-defining input) and
  `cellchat_summary_statistic.csv`.
- The human arm's median nearest-neighbour distance sits **below** the frozen λ-grid
  floor in every section. Declared, not patched.
- SenePy ships no spleen signature; cross-tissue surrogates are used and named.
- The two arms use different DeepScence estimators — a multi-seed consensus on H1,
  the frozen single-seed score on M1 — because the frozen configuration is not
  seed-reproducible on human data. The asymmetry is reported.
- Species and tissue both differ between arms, so **no cross-arm difference can be
  attributed to either**.
- 21 pre-existing prohibition violations remain in committed prose; the checker
  ratchets on new text only.
