# How Far Does Senescence Actually Reach?

## Confounding and Identifiability in Estimating SASP Spatial Response Kernels

**Master project document: research plan, data audit protocol, biology scope, and compute setup**

**Target:** NeurIPS 2026 workshop paper
**Primary venue:** ml4spatialbio (ML for Spatially Resolved High-dimensional Biology), Paris
**Secondary venue:** ICBINB-BIO (Failure Modes of AI in Biology), Sydney
**Team:** 1 CS lead (Nicole Ching) + 1 biology collaborator
**Deadline:** August 29, 2026, 11:59 PM AoE
**Last updated:** August 17, 2026

---

## Contents

**Part 0 — Orientation**
- 0.1 Do These Six Things Today
- 0.2 Deadline Reality Check
- 0.3 The One Thing That Kills This Project

**Part I — The Science**
- 1. Project Summary
- 2. Novelty Verdict: What Is Taken and What Is Left
- 3. Biology Primer
- 4. State of the Field
- 5. The Gap We Fill
- 6. Formal Problem Statement

**Part II — Data and Biology**
- 7. Where to Look for Data
- 8. The Day 1 Audit: Six Executable Tests
- 9. Gene Set Architecture
- 10. Sender Calling
- 11. The Liver Zonation Problem
- 12. Biology Collaborator Scope of Work
- 13. Handoff Contract

**Part III — Compute (RunPod)**
- 14. Why RunPod, and How to Use It Correctly
- 15. Account and Pod Setup
- 16. Environment and Docker Image
- 17. Storage and Data Transfer
- 18. Making the Compute Small
- 19. Cost Control
- 20. Failure Modes
- 21. AWS Fallback

**Part IV — Experiments**
- 22. Experimental Design
- 23. Baselines and Null Battery
- 24. Statistical Analysis Plan
- 25. Figures

**Part V — Execution**
- 26. Twelve-Day Timeline
- 27. Budget
- 28. Risks and Mitigations
- 29. Positioning and Reviewer Objections
- 30. Paper Outline
- 31. Key References
- 32. Reading List

---

# Part 0 — Orientation

## 0.1 Do These Six Things Today

Before reading the rest of this document. Roughly four hours.

1. **Open GEO and confirm GSE310392 is downloadable**, and check whether `cell_feature_matrix.h5` and `cells.parquet` are deposited per sample. (Section 7)
2. **Create a RunPod account and add credits.** Ten minutes. (Section 15)
3. **Email Caltech research computing** asking about interactive nodes with 64+ GB RAM. Free beats paid, and you want the reply in your inbox by Day 3.
4. **Confirm your biology collaborator** and send them Section 12. Get a yes or no on availability for eight of the next twelve days, in writing.
5. **Read Zhao et al. (2024), *Cell*.** It is the closest prior work and it will change how you frame the paper. (Section 32)
6. **Create the OpenReview submission** as a placeholder on the ml4spatialbio site so you know the template and page limit before you write anything.

## 0.2 Deadline Reality Check

NeurIPS 2026 announced accepted workshops on August 10, 2026. The suggested contribution deadline is **August 29, 2026 AoE**, confirmed on the ml4spatialbio, ICBINB-BIO, and Interpretability for Discovery sites. Author notification is September 29 and NeurIPS states it cannot be extended.

**Twelve days, including writing.**

Three consequences that shape every decision below:

- **Scope is fixed and small.** Four figures. One primary tissue plus one replication cohort. Do not add a third.
- **Nothing on the critical path may wait on anyone's email.** Every dataset and tool named here is public today.
- **There are two hard gates.** Day 3: synthetic recovery working and data audited. Day 6: Figure 1 and Figure 2a exist. Miss either and stop, rather than submitting something thin.

An honest alternative: if the Day 1 audit comes back badly, do not force August 29. Build the full version, add a second tissue and a proper benchmark release, and target ICLR 2027 workshops or the spring ML-for-biology circuit. A rushed workshop paper with a weak causal claim is worse for you than no paper, because it becomes the first thing reviewers find when they search your name.

## 0.3 The One Thing That Kills This Project

**Circular measurement.** If the gene set used to call a cell senescent overlaps with the gene set used to measure its neighbor's response, you get a beautiful, smooth, highly significant decay curve that measures nothing but spatial autocorrelation of one score against itself.

This is the default outcome, not an edge case, because the field's most popular senescence gene set is mostly SASP genes. **SenMayo is 125 genes and 83 of them are SASP factors**; the remainder are 20 transmembrane and 22 intracellular proteins. Score senders with SenMayo, measure any SASP-adjacent response in neighbors, and two-thirds of your sender definition is the thing you claim to be measuring downstream.

The governing rule for the whole project:

> **Senders are defined by arrest and damage. Receivers are measured by secretion-response. The two gene lists must have zero intersection, and the intersection matrix must be printed in the Methods to prove it.**

Section 9 is built entirely around enforcing this.

---

# Part I — The Science

## 1. Project Summary

**The question.** A senescent cell secretes SASP factors that alter neighboring cells. That influence must fall off with distance. Everyone assumes this; one paper has measured it descriptively; nobody has estimated it as a parameter with error bars, and nobody has checked whether the estimate survives the confounds known to inflate every other spatial cell-cell communication method.

**Goal.** Estimate the senescence spatial response kernel: for a given SASP-driven transcriptional program, how does a receiver cell's response scale with distance to the nearest senescent sender? Recover a length constant λ with confidence intervals, per program, per receiver cell type, per tissue. Then determine how much of that estimate is real and how much is an artifact of senescent cells sitting in structurally unusual neighborhoods.

**Why this framing rather than "build a better spread model."** A better spread model is not the bottleneck. At least five published methods already model ligand diffusion and decay from spatial data. The bottleneck is that a January 2026 benchmark found CellChat v2, COMMOT, and SpaTalk report similar interaction counts on real and coordinate-randomized data, implying false positive rates above 90%. A June 2026 method paper showed that even weak spatial autocorrelation inflates type I error across every competing ligand-receptor method tested. If those tools are that badly calibrated, every published statement about how far senescence reaches is standing on sand, and the first useful contribution is to find out by how much.

**The ML contribution.** An identifiability analysis. Senescent cells are spatially clustered, not randomly placed, so distance-to-nearest-senescent-cell is correlated with local density, cell-type composition, tissue architecture, and the receiver's own baseline state. Under what sender-clustering regimes is the true kernel recoverable from a static spatial snapshot, and where does it become unidentifiable? Answered with a synthetic-ground-truth study, then applied as calibrated caution to real tissue.

**The biology contribution.** A confounder-controlled estimate of how far senescence reaches, and whether the answer is conserved across tissues or specific to senotype. Plus a data-grounded constraint on the open containment question: Martin et al. (*Aging Cell* 2023) showed that current mechanistic understanding cannot explain why senescence spread stays local rather than running away through the tissue. A measured kernel narrows the candidate mechanisms.

**Deliverable.** A 4–9 page workshop paper with 4 main figures, code, and a reusable null battery.

## 2. Novelty Verdict: What Is Taken and What Is Left

Stated plainly: **partially novel, and less novel than it first looks.** The naive version has been done. The version in this document has not.

### 2.1 Taken

| Piece | Who | Status |
|---|---|---|
| Distance-dependent SASP gradients around senescent cells in tissue | Zhao et al., *Cell* 2024. Defined "senescence-sensitive spots," ranked genes by proximity, showed SASP scores rise and TNF signaling falls with distance, across multiple aged mouse organs | **Done, descriptively.** This is in *Cell*. You cannot claim the observation |
| Mapping senescent cells and their neighbors in spatial transcriptomics | scDOT, *Genome Biology* 2024 (Bar-Joseph lab). Optimal transport plus deconvolution; spatial organization of senescent cells in lung, candidate interaction genes | Done |
| Explicit ligand diffusion, competitive binding, and concentration decay in a fitted model | SCILD, *Communications Biology* 2026. Single-cell resolution, interpretable optimization, in silico perturbation | Done, generically |
| Learning niche-composition effects on expression at characteristic length scales | NCEM, *Nature Biotechnology*. GNN; recovered length scales matching known communication mechanisms | Done, generically |
| Distance-weighted attention for multi-scale communication | RGAST, bioRxiv 2024 | Done |
| Mathematical model of senescence spread and the containment paradox | Martin et al., *Aging Cell* 2023 | Done; the paradox is left open, and the paper declines to say whether spread is finite |
| Classifying primary vs secondary (bystander) senescence with ML | PMC11689308; paracrine spread between human brain cell types, bioRxiv Feb 2026 | Done |
| Boolean network model of SASP onset after DNA damage | *PLOS Comput Biol* 2017 | Done |
| Senescent cell detection in spatial data | DeepScence, *Cell Genomics* Dec 2025 | Done, and you will use it |

### 2.2 Left

1. **No confounder-controlled estimate.** The *Cell* 2024 gradient analysis has no negative controls, does not test against matched non-senescent decoys, does not control local density or cell-type composition, and reports no length constant with uncertainty. Given the >90% FPR finding, the prior probability that a naive gradient survives a proper null is not high.
2. **No identifiability analysis.** Every method above estimates decay from a static snapshot. A July 2026 arXiv paper asks in general when static spatial patterns identify diffusion dynamics and finds real limits. Nobody has asked this for the senescence case, where sender clustering is severe and non-random by construction.
3. **Signal decay and response decay are conflated.** Ligand concentration falling off is physics. Cells becoming less responsive is biology. Every current method estimates their product and reports one number.
4. **No cross-tissue comparison of λ.** SenNet's central finding is that senescence is not one state. If senotypes differ, their reach probably differs too. Untested.
5. **No benchmark or evaluation standard.** ml4spatialbio explicitly lists benchmarks, datasets, and evaluation standards for spatial tasks as a topic of interest.

### 2.3 The framing decision that follows

Pursue this as an **evaluation and identifiability paper**, not a new-method paper.

Framed as "we built a better SASP spread model," a reviewer who knows SCILD, NCEM, and COMMOT rejects it in ten minutes. Framed as "here is what the field's central spatial claim about senescence looks like under proper nulls, and here is when it is even recoverable in principle," you have something nobody has, at a venue that explicitly wants it.

## 3. Biology Primer

For the CS lead, and for the biology collaborator on the modeling terms.

**SASP (senescence-associated secretory phenotype).** The cytokines, chemokines, proteases, and growth factors secreted by senescent cells. IL-6, IL-8/CXCL8, CCL2, TGF-β, GDF15, and the IGFBPs are the recurring names. The SASP Atlas (Basisty et al., *PLoS Biology* 2020) is the reference proteomic catalog and is the right citation for what is actually secreted, as opposed to which transcripts go up.

**Paracrine / bystander / secondary senescence.** A senescent cell can push a healthy neighbor into senescence via secreted factors. Demonstrated with transwell experiments (Acosta et al. 2013), so physical contact is not required. Secondary senescent cells are transcriptomically distinct from primary ones, particularly in their own SASP output.

**Senotype.** SenNet's term for the fact that senescence is not one state. A senescent hepatocyte, astrocyte, and B cell share almost no markers. SenCat, profiling 14 primary human cell types across 30+ senescence paradigms, found no universal marker.

**Why decay must exist.** If every senescent cell reliably senesced its neighbors, senescence would propagate without limit. It does not. Martin et al. built a minimal mathematical model showing that current mechanistic understanding of SASP diffusion and binding does not explain how spread stays local — while stating explicitly that the model *"can not determine … whether the spread of senescence is controlled (finite) or uncontrolled."* Their own proposed resolution comes from the secondary senescent cells: that they are poor SASP producers and act as a firebreak, and that induction is time-delayed rather than instantaneous. The field's standing candidates for what limits the *ligand field* — rapid degradation, a response threshold, immune clearance of senders, receptor-level refractoriness in receivers — are **not** from that paper.

> **Correction 2026-08-21 (D7 §B4b).** An earlier version of this section attributed those four mechanisms to Martin et al. Verified against the open-access full text (`PMC10410058`): `degrad` and `refractor` occur zero times, and immune clearance appears only as the mechanism the paper brackets out. The four remain a useful hypothesis space; they are the plan's gloss, not a citation.

**Length scales to expect.** Cytokine signaling in tissue is generally thought to act over tens of microns, not hundreds. For calibration: a 2026 human endometrium study measured nearest-neighbor distances from senescent cells to macrophages and monocytes at roughly 45 ± 20 μm and 45 ± 25 μm, to NK cells at 53 ± 23 μm, to T-helper cells at 102 ± 42 μm, and to B cells at 211 ± 66 μm. A single cell diameter is roughly 10–20 μm. **If your fitted λ comes out at 500 μm, suspect the model, not the biology.**

**Spatial transcriptomics platforms, and why the choice matters.**
- *Single-cell resolution* (Xenium, MERFISH, CosMx, seqFISH+): individual cells with coordinates, targeted panel of hundreds to thousands of genes. **You need this.** Distance is only meaningful if you know where individual cells are.
- *Spot-based* (Visium, Visium HD, Slide-seq): whole or near-whole transcriptome, but each spot mixes multiple cells. Distance is blurred by spot size and requires deconvolution, which adds its own spatially structured confound.

**Spatial autocorrelation.** Neighboring cells resemble each other for reasons unrelated to signaling: shared microenvironment, shared clonal origin, technical diffusion of transcripts, segmentation bleed-through. This is the largest source of false positives in the field and the thing your nulls must beat.

## 4. State of the Field

### 4.1 Senescence in space

The 2023 *Nature Aging* SenNet review, "Spatial mapping of cellular senescence," is the field's own statement of the problem and is your framing citation. The June 2026 *Cell* package includes a spatial multi-omics atlas of human lymph nodes across ages 18–86 (Farzad et al.), which found senescent-like B cells shifting from interfollicular zones into germinal centers with age. The SenNet portal listed 1,753 public human and mouse datasets across 15 organs and 6 assay types as of January 2026.

**Ma et al. (*Cell* 2024) is the closest scientific precedent and you must engage with it directly.** [Corrected 2026-08-21, D7 §B1/B2/B3: cited throughout earlier drafts as "Zhao et al."; the first author is **Ma S** and Zhao L is 14th of 47. The platform is **Stereo-seq, spot/bin level, not single-cell**, and the design is **two age groups (young, old)**, not a time course.] They profiled **young versus old male mice across nine tissues** on Stereo-seq bins, defined senescence-sensitive spots, and showed that SASP score, TNF signaling, ATP biosynthesis, and cell-cycle genes all vary monotonically with distance from those spots, consistently across organs. Their conclusion, that senescent foci act as epicenters compromising surrounding cells in a distance-dependent manner, is exactly the phenomenon you propose to quantify. Your contribution is quantification with controls, not discovery.

### 4.2 Spatial CCC methods

NCEM (*Nature Biotechnology*), COMMOT (*Nature Methods*), SpaTalk, CellChat v2, GCNG, DT-CCC, GraphST, NICHES, RGAST, and SCILD (*Communications Biology* 2026, which fits explicit ligand diffusion, competitive binding, and concentration decay in one optimization). HARMONIC (bioRxiv Jan 2026) adds H&E histology to condition on tissue context and reduce false positives.

### 4.3 The reliability problem, which is why this project exists

- **CellWHISPER** (bioRxiv Jan 2026) benchmarked CellChat v2, COMMOT, and SpaTalk on coordinate-randomized data and found they report similar interaction counts on real and randomized input, implying **false positive rates above 90%**. Their own confounder-aware null brought this under 5%.
- **CONCISE** (bioRxiv June 2026) showed that introducing even weak spatial autocorrelation (a = 0.1) into one gene inflated type I error for every competing spatial ligand-receptor method tested.
- A March 2026 *Briefings in Bioinformatics* benchmark of cell-type-specific spatially variable gene methods found rotation invariance unresolved and cell-type confounding a central open challenge.

The field has, in the last eight months, started admitting that spatial CCC inference is poorly calibrated. Nobody has yet asked what that means for the specific claims the senescence field has built on top of these tools.

## 5. The Gap We Fill

> The claim that senescent cells influence their neighbors in a distance-dependent manner is central to senescence biology and has been shown descriptively. It has never been estimated as a parameter with uncertainty, never tested against confounder-aware nulls, and never subjected to an identifiability analysis, despite contemporaneous evidence that the class of methods used to make such claims has false-positive rates above 90% under randomization.

Three contributions:

1. **A senescence response kernel estimator** with an explicit, interpretable length constant, fit per SASP-driven program, per receiver cell type, per tissue, with donor-level confidence intervals.
2. **A null battery** designed for the sender-clustering confound specifically: matched-decoy senders, cell-type-stratified label permutation, torus shift, rotation, density and composition conditioning, receiver-baseline conditioning, threshold sensitivity, and gene-set disjointness. Reported as a *surviving fraction* of the naive estimate.
3. **An identifiability study on synthetic tissue** with a planted ground-truth kernel, characterizing the regimes of sender density, sender clustering, and noise where λ is recoverable and where it is not. Fully under your control and independent of data access, which is why it is scheduled first.

## 6. Formal Problem Statement

Let a tissue section contain cells $i = 1 \dots N$ at positions $x_i \in \mathbb{R}^2$, each with cell type $c_i$ and expression $y_i$. Let $S$ be the set of senescent senders. For each non-sender cell $i$:

$$d_i = \min_{j \in S} \|x_i - x_j\|$$

and let $r_i$ be the score of a SASP-driven target program in cell $i$.

### 6.1 Model

$$r_i = \mu_{c_i} + \beta \cdot K_\lambda(d_i) + \gamma^\top z_i + \varepsilon_i$$

where $K_\lambda$ is a decay kernel, $z_i$ collects nuisance covariates (local cell density, k-NN cell-type composition, total counts, donor, section, distance to tissue boundary, anatomical landmarks), and $\mu_{c_i}$ is a receiver-cell-type intercept.

### 6.2 Kernel families

Fit all, compare by AIC and held-out log-likelihood on left-out sections:

- **Exponential:** $K_\lambda(d) = e^{-d/\lambda}$. Natural form for degradation-limited diffusion.
- **Gaussian:** $K_\lambda(d) = e^{-d^2/2\lambda^2}$.
- **Power law:** $K_\lambda(d) = (1 + d/\lambda)^{-p}$. Heavier tail; distinguishes local from long-range.
- **Step / threshold:** $K(d) = \mathbb{1}[d < \lambda]$. Tests whether the response is graded at all.
- **Nonparametric spline in $d$.** Reference for whether any parametric family is adequate.

### 6.3 Aggregate-sender variant

Nearest-sender distance discards the fact that a cell with five senescent neighbors sees more signal than a cell with one. Also fit superposition:

$$r_i = \mu_{c_i} + \beta \sum_{j \in S} K_\lambda(\|x_i - x_j\|) + \gamma^\top z_i + \varepsilon_i$$

If superposition wins, that is evidence for a dose-like SASP effect. If nearest-neighbor wins, that is evidence for a threshold effect. Either result is interesting and neither has been reported.

### 6.4 Separating signal decay from response decay

Fit the kernel twice: once with $r_i$ = a **receptor-proximal** program (immediate-early NF-κB targets, which should track ligand concentration closely) and once with $r_i$ = a **downstream** program (cell-cycle arrest, secondary SASP).

If $\lambda_{\text{proximal}} > \lambda_{\text{downstream}}$, the response is thresholded rather than graded, which speaks directly to the containment question (not, as earlier drafts had it, to a threshold mechanism proposed by Martin et al. — see the correction in §3). This comparison is cheap and nobody has published it.

### 6.5 The quantity to report

Not just $\hat\lambda$, but for each null in the battery:

$$\text{Surviving fraction} = \frac{\hat\beta_{\text{controlled}}}{\hat\beta_{\text{naive}}}$$

with bootstrap confidence intervals. Two numbers, always, side by side.

---

# Part II — Data and Biology

## 7. Where to Look for Data

Work down this list. Stop when something passes Section 8.

### Rank 1 — GSE310392, human liver senescence atlas

**"Cellular senescence in human liver under normal aging and cancer," *Cell Genomics*, January 2026.** 43 normal human livers spanning ages and fibrosis stages, plus 24 colorectal cancer liver metastases. Xenium spatial transcriptomics with the 5K panel, single-cell multiome, and CODEX protein imaging.

Why this is the top candidate:

- **Single-cell resolution.** Xenium, not Visium spots.
- **5,000-gene panel.** Enough to build a sender score and a fully disjoint response score. Most targeted panels cannot do this.
- **43 donors.** You can bootstrap over donors, which almost no spatial study allows.
- **Senescence already annotated.** CDKN1A+ senescent hepatocytes, fibroblasts, cholangiocytes, endothelial cells; two senescent CAF populations in tumor stroma (CDKN2A/2B and CDKN1A/SERPINE1).
- **They left your question open.** They compared clustered versus isolated senescent hepatocytes, examined the surrounding microenvironment, found trends toward more hepatic stellate cells and inflammatory macrophages near clustered senescent cells, and **explicitly reported this as falling short of significance**. They wrote "possible paracrine dissemination of senescence" as a suggestion, not a finding.

That last point is the whole opportunity. They had the right data, saw the trend, and could not call it. A properly powered, confounder-controlled kernel estimate is exactly the missing analysis.

**Caution: read Section 11 before touching this dataset.** They report CDKN1A+ hepatocytes localize periportally, within 100–150 μm of the portal triad in young individuals. That is a severe anatomical confound and it must be modeled explicitly.

### Rank 2 — SenNet portal (`data.sennetconsortium.org`)

1,753 public human and mouse datasets across 15 organs and 6 assay types as of January 2026. Filter for spatial assays at single-cell resolution. The consortium exists to map senescence, so annotation quality should exceed generic atlases. API docs at `docs.sennetconsortium.org/apis`.

Specific targets: **Farzad et al. (2026) human lymph node spatial multi-omics, ages 18–86**, from the June 2026 *Cell* package. Lymph node has strong architectural structure, which is a confound but a well-defined one you can control for. Also check lung parenchyma, prefrontal cortex, and liver TMC releases.

### Rank 3 — Ma et al. mouse aging atlas

*Cell* 2024. Nine tissues, **two age groups (young, old)** — not "multiple life stages". This is the study whose gradient you would be re-examining.

> **Corrected 2026-08-21 (D7 §B2/B3).** This rank was written on two wrong premises. The data are **Stereo-seq at spot/bin level** — 1,535,191 spots at ~1,450 genes/spot — **not single-cell resolution**, so the plan under "If deposited at single-cell resolution, reproducing their distance-ranked gradient … is your strongest Figure 2" **was never executable as written**. And the data are in **CNGB STOMICS DB (`STDS0000247`)**, neither GEO nor GSA. This is in our favour and should be said out loud: the field's closest distance-gradient result was measured at a resolution that blurs distance, which §3 of this plan warns is unsuitable for kernel estimation.

### Rank 4 — 10x Genomics public Xenium datasets

Free, immediate, no application. Xenium Prime 5K Human Pan Tissue and Pathways Panel datasets exist for breast, ovarian, and cervical cancer among others. Tumor tissue has high senescent burden. Downside: no aging axis, no senescence annotation, and cancer architecture is chaotic in ways that complicate the nuisance model.

**Use as a replication cohort, not a primary.** Showing your kernel estimate reproduces in an independently generated dataset with a different panel is worth a lot at review.

### Rank 5 — CELLxGENE Discover, HuBMAP

Broad spatial holdings, slow to search. Do not spend Day 1 here.

### What not to use

- **Visium (non-HD) as primary.** 55 μm spots containing multiple cells. Your λ is plausibly 20–100 μm, so spot size is comparable to the quantity you are estimating, and deconvolution error adds a spatially structured confound. Visium HD (2 μm bins) is acceptable with binning to cell-like units; regular Visium is not.
- **Any single-section dataset.** No donor replication means no valid uncertainty, and your central claim is about uncertainty.

## 8. The Day 1 Audit: Six Executable Tests

Run in order. Each is a small script. Budget four to six hours. Record every number; they go in Methods and the supplement.

### Test 1 — Resolution and segmentation quality

- Confirm one row per cell with (x, y) coordinates. **Confirm the units are microns**, not pixels. Put the unit in every column name.
- Report median transcripts per cell and median genes detected per cell.
- Report **median nearest-neighbor cell distance**. This is your resolution floor. You cannot claim to resolve λ below it. If median NN distance is 15 μm, do not report λ = 8 μm.
- Check what fraction of transcripts are assigned to a cell. Unassigned rates above roughly 30% mean bleed-through, which manufactures spatial autocorrelation.

**Pass:** single-cell resolution confirmed, NN distance recorded, assignment rate documented.

### Test 2 — Panel adequacy and disjointness

The test people skip and the one that most often fails.

```python
A = sender_genes & panel_genes
B = response_genes & panel_genes

assert len(A) >= 15
assert len(B) >= 30        # per response module you plan to fit
assert len(A & B) == 0     # non-negotiable
```

If `A & B` is non-empty, remove overlapping genes from the **sender** set, not the response set, then re-check `len(A) >= 15`. If you cannot reach 15 disjoint sender genes, the panel is inadequate. Next dataset.

**Pass:** ≥15 sender genes, ≥30 genes per response module, empty intersection, all on-panel.

### Test 3 — Sender prevalence

Score senescence, threshold, count.

- **Below ~1% of cells:** too few senders. Distance-to-nearest-sender is dominated by a handful of foci and the donor bootstrap is uninformative.
- **Above ~20%:** distance-to-nearest-sender is near zero everywhere. No dynamic range in the independent variable; λ is unidentifiable by construction.
- **Sweet spot: 2–10%.**

Report prevalence at the 90th, 95th, and 99th percentile thresholds and at DeepScence's own cutoff, **per cell type**, not just overall.

**Pass:** at least one tissue-threshold combination in 1–20%, with ≥200 senders and ≥5,000 non-senders per donor.

### Test 4 — Sender clustering

Compute Ripley's K or the pair correlation function for sender positions against (a) a Poisson null and (b) a null that permutes sender labels within cell type.

- **Near-random senders:** best case, λ most identifiable.
- **Strongly clustered senders:** expected, and it is the confound the paper is about. Record the statistic; it determines which synthetic regime in Figure 1 you are actually in.

If senders form a handful of large blobs rather than many small foci, note that your effective sample size is the number of blobs, not the number of cells.

**Pass:** none, this is a measurement.

### Test 5 — The decoy contrast test (the go/no-go)

The most important hour of Day 1.

1. For each sender, select a matched decoy: a non-senescent cell of the same type, matched on local density (k-NN count within 50 μm) and k-NN cell-type composition. Use nearest-neighbor matching on a propensity score, or coarsened exact matching.
2. Compute the distribution of distance-to-nearest-sender and distance-to-nearest-decoy across all remaining cells.
3. Compare.

**If the distributions are nearly identical, you have a contrast and the project works.** Any difference in response can then be attributed to sender identity rather than neighborhood structure.

**If they differ substantially,** senders occupy structurally distinct neighborhoods that decoy matching cannot equalize. Report the imbalance, add the imbalanced covariates to the nuisance model, re-check. If it still fails, this tissue cannot support the estimate, which is itself a finding worth a paragraph.

**Pass:** standardized mean difference < 0.1 on all matching covariates after matching.

### Test 6 — Anatomical confounding

- Compute distance from each cell to relevant landmarks: tissue boundary, large vessels, portal triads (liver), germinal centers (lymph node), airways (lung).
- Regress sender status on those distances. Report pseudo-R².
- **If sender status is largely predicted by distance to one structure, then distance-to-sender is a proxy for distance-to-that-structure**, and any kernel you fit is measuring anatomy.

Not a reason to abandon the dataset. A reason to include the anatomical distance as a covariate and report the kernel with and without it. But you must know before you start.

**Pass:** none, this is a measurement, but it determines your nuisance model.

### Audit summary table

| Test | Metric | Value | Pass? |
|---|---|---|---|
| 1 | Median NN cell distance (μm) | | |
| 1 | Transcript assignment rate | | |
| 2 | Sender genes on panel | | ≥15 |
| 2 | Response genes on panel (per module) | | ≥30 |
| 2 | Sender ∩ response | | must be 0 |
| 3 | Sender prevalence at 95th pct | | 1–20% |
| 3 | Senders per donor | | ≥200 |
| 4 | Ripley's K deviation from Poisson | | record |
| 5 | Max SMD after decoy matching | | <0.1 |
| 6 | Pseudo-R² of anatomy → sender status | | record |

## 9. Gene Set Architecture

Five tiers. Tiers A and B must not intersect. Print the intersection matrix in the supplement.

### Tier A — Sender definition (arrest and damage only)

Design constraint: **no secreted factors.** Every gene here is about the cell's own cycle arrest or damage state, not about what it is releasing.

**A1. Core arrest markers (up):**
`CDKN1A` (p21), `CDKN2A` (p16), `CDKN2B` (p15), `TP53`, `TP53I3`, `GADD45A`
Note: `SERPINE1` (PAI-1) is borderline because it is secreted. Exclude it if you use ECM response modules.

**A2. Proliferation markers (down):**
`MKI67`, `TOP2A`, `PCNA`, `CCNB1`, `CCNA2`, `CDK1`, `BIRC5`, `TYMS`, `RRM2`, `MCM2`–`MCM7`

**A3. Nuclear envelope and chromatin (down):**
`LMNB1`, `HMGB1`, `HMGB2`, `HIST1H1` family

**A4. DNA damage response:**
`ATM`, `ATR`, `CHEK1`, `CHEK2`, `H2AFX`, `TP53BP1`, `MDM2`

**A5. Curated sets to intersect with the above:**
- **SenCat ML-refined transcriptomic signature** (Anerillas et al., *Mol Cell* 2026). Newest, derived across 14 human primary cell types and 30+ paradigms. Preferred.
- **SenePy** cell-type-specific weighted signatures, explicitly built to recapitulate in vivo senescence better than in vitro-derived signatures. Use the signature matching your receiver cell type. (Verify the exact signature counts per species against the paper.)
- **CellAge** and **GenAge** from the Human Aging Genomic Resources.
- **SeneQuest**, filtered as DeepScence did: retain a gene if reported by ≥15 publications, or by ≥4 publications with ≥70% agreement on direction.
- **CoreScence**, the curated set underlying DeepScence, built by integrating nine published sets.

**A6. Handle carefully or exclude:**
- **SenMayo (125 genes).** 83 are SASP factors. **Do not use as the sender score.** For comparison only, define `SenMayo_arrest = SenMayo \ SASP_genes` and report both the reduced set size and the result.
- **SenSig (Cherry et al.).** Over a thousand genes. Too large; it will overlap with everything.
- **Reactome R-HSA-2559582 (SASP).** This is a *response* set, not a sender set. Note also that the SenMayo authors reported GSEA on R-HSA-2559582 failed to detect an age-related increase in their cohorts, so it is a weak discriminator.

**A7. The heterogeneity problem, for your Limitations section.** A 2026 *FEBS Open Bio* comparison of nine senescence resources found substantial variability in set size and limited intersection, with a core of only 22 genes shared. SenCat found no universal marker across 14 cell types. **Your sender call is a choice, not a fact.** Test 3 and the threshold sensitivity null (N7) are how you handle it.

### Tier B — Response definition (what the neighbor does)

Fit the kernel separately per module. The proximal-versus-downstream comparison is one of your cheapest novel results.

**B1. Receptor-proximal / immediate-early** (should track ligand concentration most closely):
`HALLMARK_TNFA_SIGNALING_VIA_NFKB`
Key members to check on-panel: `NFKBIA`, `NFKB1`, `NFKB2`, `RELB`, `TNFAIP3`, `JUNB`, `FOS`, `FOSB`, `EGR1`, `ZFP36`, `SOD2`, `ICAM1`, `CXCL2`, `PTGS2`, `BIRC3`

**B2. IL-6 / JAK / STAT3 axis:**
`HALLMARK_IL6_JAK_STAT3_SIGNALING`
Key members: `SOCS3`, `STAT3`, `IL6R`, `IL6ST`, `OSMR`, `JAK1`, `JAK2`, `PIM1`, `MYD88`

**B3. Interferon response:**
`HALLMARK_INTERFERON_GAMMA_RESPONSE`, `HALLMARK_INTERFERON_ALPHA_RESPONSE`
Included because cross-model work on aging in foundation models found TNF/NF-κB and type-II IFN-γ to be the two submodules on which two independent models agreed. If those are the robust aging inflammation axes, they are the right response readouts here.

**B4. Downstream arrest** (should show a shorter λ if the response is thresholded):
`HALLMARK_E2F_TARGETS` (expect down near senders), `HALLMARK_G2M_CHECKPOINT` (down), `HALLMARK_MYC_TARGETS_V1`

**B5. Fibrotic / ECM remodeling:**
`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`; `COL1A1`, `COL1A2`, `COL3A1`, `FN1`, `TIMP1`, `ACTA2`

**B6. Oxidative stress:**
`HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY`; NRF2 targets `NQO1`, `HMOX1`, `GCLM`, `TXNRD1`

**B7. Secondary senescence** (the strongest biological claim, if supportable):
Does a neighbor near a sender become senescent itself? Score neighbors with Tier A **excluding** the genes used to call the original senders, or better, use a secondary-senescence-specific signature. Primary and secondary senescent cells are transcriptomically distinct, particularly in SASP output. If your kernel fits for this readout, that is the direct measurement of senescence spread and it belongs in the abstract.

**Verify all MSigDB set sizes and membership against the current release.** Do not cite gene counts from memory, including from this document.

### Tier C — Ligand-receptor pairs (mechanistic interpretation)

| Ligand (sender) | Receptor (receiver) | Notes |
|---|---|---|
| `IL6` | `IL6R`, `IL6ST` | Canonical SASP; drives B2 |
| `CXCL8` (IL-8) | `CXCR1`, `CXCR2` | Canonical SASP |
| `CCL2` | `CCR2`, `ACKR3` (CXCR7) | **Highest priority.** Key SASP factor in paracrine spread; CXCR7 selectively expressed in brain cell types capable of secondary senescence |
| `CXCL12` | `CXCR4`, `ACKR3` | `DPP4` cleaves CXCL12 and is itself senescence-associated; a candidate range-limiting mechanism |
| `TGFB1` | `TGFBR1`, `TGFBR2` | Contact-adjacent, short range expected |
| `IL1A`, `IL1B` | `IL1R1`, `IL1RAP` | IL-1α is largely membrane-bound, predicting very short λ |
| `TNF` | `TNFRSF1A`, `TNFRSF1B` | Drives B1 |
| `IGFBP3`, `IGFBP7`, `GDF15`, `MMP1`, `MMP3`, `TIMP1` | various | SASP Atlas core secreted factors |

**A built-in internal control.** Ligands differ in expected range. Membrane-bound IL-1α should give the shortest λ; freely diffusible small chemokines like CXCL8 the longest. **If your fitted λ values do not respect that ordering, your estimate is picking up something other than diffusion**, which is exactly the kind of finding this paper exists to report.

### Tier D — Nuisance covariates

These go in $z_i$. Underspecifying this is how the naive estimate gets inflated.

- Receiver cell type (random or fixed intercept)
- Local cell density: k-NN count within 25, 50, and 100 μm
- k-NN cell-type composition vector (within 50 μm)
- Total transcript count and genes detected per cell (technical)
- Donor, section, batch
- Distance to tissue boundary
- **Tissue-specific anatomical landmarks** (Section 11)

### Tier E — Negative control gene sets

- **Housekeeping:** `ACTB`, `GAPDH`, `RPL13A`, `RPS18`, `TBP`, `PPIA`. A kernel fit against these should be flat. If it is not, you have a technical gradient (transcript density, tissue thickness, hybridization efficiency) and everything else is suspect.
- **Cell-type identity program unrelated to inflammation:** e.g. `ALB`, `TTR`, `APOA1` in hepatocytes. Should be flat after conditioning on cell type.
- **Random matched sets:** 500 random gene sets matched to each response module on size and mean expression. Gives an empirical null for effect size.

## 10. Sender Calling

Four options, in the order to try them.

**1. DeepScence** (*Cell Genomics*, Dec 2025; Qu et al., Duke; 5(12):101035). A deep neural network built on the CoreScence curated gene set, explicitly designed for **single-cell and spatial** data, validated on Visium and Stereo-seq and on simulated Xenium panels with varying numbers of senescence genes on-panel. The only tool built for your data type. **Start here.**

**2. SenePy.** Cell-type-specific weighted signatures built from in vivo data and shown to recapitulate in vivo senescence better than in vitro-derived signatures. Use as the independent second opinion. If DeepScence and SenePy disagree on more than roughly 30% of calls, report the disagreement and run the full analysis under both.

**3. SenCat ML signatures.** Newest (June 2026), SenNet-derived, transcriptome and proteome. Third opinion.

**4. SenCID.** Comparison only, with the limitation stated openly: DeepScence's authors note SenCID was trained on in vitro bulk data and that it is questionable whether it transfers to in vivo settings, and that its models may not apply directly to Xenium-style data with different distributions and gene counts.

**Reporting requirement:** report the headline kernel under at least two independent sender-calling methods. If λ moves substantially between them, that is a limitation stated in the abstract, not a footnote.

## 11. The Liver Zonation Problem

Worth its own section because it is the clearest illustration of what the paper is about, and because it applies directly to the Rank 1 dataset.

The liver is zonated. Hepatocytes near the portal triad and near the central vein have systematically different programs. Periportal: `ASS1`, `SDS`, `HAL`, `CPS1`, `ALB`. Pericentral: `GLUL`, `CYP2E1`, `CYP1A2`, `OAT`. This gradient is one of the strongest spatial signals in the tissue and has nothing to do with senescence.

The *Cell Genomics* liver paper reports that **CDKN1A+ hepatocytes preferentially localize near the portal triad, within 100–150 μm, in young individuals.**

Put those together:

1. Senders cluster periportally.
2. Periportal hepatocytes have a distinct program for zonation reasons.
3. Therefore distance-to-nearest-sender correlates strongly with distance-to-portal-triad.
4. Therefore **any response program that varies with zonation produces a clean, monotone, highly significant "decay curve" around senescent cells that is entirely zonation.**

Several Tier B modules plausibly vary with zonation. Oxidative stress (B6) certainly does; pericentral hepatocytes are the site of CYP450-driven oxidative metabolism.

**What to do:**

- Include distance-to-portal-triad and distance-to-central-vein as explicit covariates. Compute from vessel annotation, or derive a continuous zonation score from the marker sets above.
- Run a **zonation-matched decoy**: match on zonation score in addition to cell type and density.
- Report the kernel with and without the zonation covariate. **That difference is a figure panel**, and it is the clearest demonstration in the paper of why the null battery is necessary.
- Stratify: fit separately within periportal, midzonal, and pericentral hepatocytes. If the kernel exists only in the pooled analysis and vanishes within zones, it was zonation.

Every tissue has a version of this. Lymph node has follicle versus interfollicular structure, and the Farzad atlas specifically reports senescent-like B cells shifting between those compartments with age. Lung has airway-to-alveolar gradients. Identify the dominant architectural axis on Day 1 (Test 6) and build the covariate before fitting anything.

## 12. Biology Collaborator Scope of Work

Roughly 40–50 hours over twelve days. Send this section to them today and get a written yes or no.

### Deliverable 1 — Dataset selection memo (Day 1, ~6 h)

Work Section 7's ranked list. For each candidate record: platform, resolution, donor count, section count, panel size, tissue, disease state, accession. Download and verify the top candidate. Produce the Section 8 audit table jointly with the CS lead.
**Output:** one page naming the primary dataset, the replication dataset, and why.

### Deliverable 2 — Gene set package (Days 1–2, ~8 h)

Assemble Tiers A–E as plain text, one gene per line, human and mouse symbols as needed. Intersect each with the panel and report on-panel counts. **Compute and report the full Tier A × Tier B intersection matrix**; resolve overlap by removing from Tier A. Pull current MSigDB Hallmark sets rather than trusting any list, including this one. Map Ensembl IDs to symbols consistently (Xenium outputs use symbols; gene sets are mixed).
**Output:** `genesets/` directory plus README with source, version, and date for every set, and the intersection matrix.

### Deliverable 3 — Anatomical annotation (Days 2–3, ~8 h)

The deliverable that most directly determines whether the paper is credible. Identify the dominant architectural axis. Define the landmark (portal triads and central veins for liver, follicle boundaries for lymph node, airway epithelium for lung). Either annotate from the H&E or DAPI image, or derive a continuous positional score from marker genes.
**Output:** per-cell table with anatomical position covariates, plus a note on what the axis is and why it matters.

### Deliverable 4 — Sender call validation (Days 2–4, ~6 h)

Run DeepScence and SenePy independently; report agreement. **Eyeball the calls against tissue images.** Are senescent hepatocytes where the literature says? Are the senescent fibroblasts actually fibroblasts? A biologist looking at fifty called cells catches errors no metric will. Flag implausible cell types.
**Output:** agreement statistics, plausibility assessment, list of cell types to exclude with reasons.

### Deliverable 5 — Ligand-receptor plausibility check (Days 6–8, ~6 h)

Once a fitted λ exists: which Tier C ligands are expressed in senders and receptors in receivers, on this panel? Is λ consistent with the ligand's expected range? Does the ordering across ligands respect known biology (Section 9, Tier C)?
**Output:** table of candidate pairs with expression evidence, and a one-paragraph verdict on whether the length scale is biologically sensible.

### Deliverable 6 — Containment interpretation and Discussion draft (Days 8–10, ~8 h)

Take the fitted kernel back to the containment question. Which limiting mechanisms are consistent — rapid degradation, response threshold, immune clearance, receiver refractoriness (the field's candidate list, **not** Martin et al.'s; see §3 correction) — and separately, what do Martin et al.'s own proposals (SASP-poor secondary cells as a firebreak, time-delayed induction) require that this design can or cannot test? If λ_proximal > λ_downstream, argue what that means for graded versus thresholded response. Draft the Discussion and the biological half of Limitations.
**Output:** Discussion draft.

### Deliverable 7 — Claim audit (Days 10–12, ~6 h)

Read every biological sentence against a source. Verify every citation resolves and says what the manuscript says it says. Flag any correlational result described in causal language.
**Output:** marked-up manuscript.

### What they do not do

Write model code, own figure generation, or decide the statistical plan. They own figure *content* review.

## 13. Handoff Contract

Agree on these by end of Day 1 and write them into `data/README.md`. Format disagreements cost more days than anything else on a short project.

| Artifact | Owner | Format | Due |
|---|---|---|---|
| Cell table | Bio → CS | Parquet/CSV: `cell_id, x_um, y_um, cell_type, donor_id, section_id, total_counts, n_genes` | Day 2 |
| Expression matrix | Bio → CS | `.h5ad`, raw counts in `.X`, cell IDs matching the cell table exactly | Day 2 |
| Gene sets | Bio → CS | `genesets/{tier}_{name}.txt`, one symbol per line, plus README with sources and intersection matrix | Day 2 |
| Anatomical covariates | Bio → CS | CSV: `cell_id, dist_to_landmark_um, zonation_score, compartment_label` | Day 3 |
| Sender labels | Bio → CS | CSV: `cell_id, deepscence_score, senepy_score, sender_flag_p90, sender_flag_p95, sender_flag_p99` | Day 3 |
| Interpretation inputs | CS → Bio | CSV: `cell_id, distance_bin, module_score, fitted_residual` | Day 6 |
| Fitted kernel table | CS → Bio | CSV: `tissue, cell_type, module, kernel_family, lambda_hat, ci_low, ci_high, null_surviving_fraction` | Day 8 |

**Coordinate units: microns, always.** Xenium outputs are in microns; some pipelines convert to pixels. Confirm on Day 1 and put the unit in every column name.

---

# Part III — Compute (RunPod)

## 14. Why RunPod, and How to Use It Correctly

**Verdict: use RunPod, but rent a low-end GPU pod and ignore the GPU. Do not use RunPod's CPU pods.**

### Why not the CPU pods

RunPod's CPU instances are compute-biased, not memory-biased. Their naming tells you the ratios: `cpu3g-8-32` is 8 vCPU with 32 GB, and `cpu5c-4-8` is 4 vCPU with 8 GB. That is 4 GB per vCPU on general purpose and 2 GB on compute-optimized. Container disk auto-sizes at vCPU × 10 GB (or × 15 GB for cpu5c) with a hard cap.

Your workload wants roughly 8 GB per vCPU. RunPod's CPU tiers force you to buy cores you will not use to reach the RAM you need, and may not get you to 128 GB at all.

### Why a GPU pod instead

RunPod ships generous **system** RAM with GPU pods; their own comparison writeups cite 167–283 GB of system RAM on shared datacenter cards. Renting a low-end GPU pod purely as a big-memory Linux box is often cheaper per GB of RAM than their CPU tiers and competitive with AWS.

As of mid-2026, Secure Cloud rates around **RTX A5000 ~$0.27/hr** and **A40 ~$0.44/hr** were being quoted, against an AWS `r7i.2xlarge` (8 vCPU, 64 GB) at ~$0.53/hr. **Shop the RunPod listings sorted by system RAM, not by GPU.** You never call CUDA. Verify current rates on RunPod's pricing page; these move.

### Why RunPod over AWS at all

You already know it from the EpiBERT project. Time to first working shell is roughly ten minutes on RunPod versus two to three hours on AWS once you account for IAM, security groups, key pairs, and VPC. On a twelve-day clock with a Day 3 gate, that is the strongest argument in RunPod's favor and it has nothing to do with hardware. Egress is also free, where AWS charges it, and Jupyter is built in rather than requiring an SSH tunnel.

### Three structural differences you must plan around

**1. Storage billing is inverted.** Container disk runs roughly $0.10/GB-month while running and roughly $0.20/GB-month while stopped. **You pay more when idle**, which is backwards from EBS and directly punishes the stop-overnight discipline. Solution: keep the container disk small, put all data on a **network volume** (from roughly $0.05/GB-month), and **terminate** rather than stop when done for the day. The network volume survives termination.

**2. Community Cloud can preempt without notice.** For a twelve-day deadline with multi-hour permutation runs, either pay for Secure Cloud or checkpoint every 50 permutations (Section 18.4). Do not run the null battery unattended on Community Cloud.

**3. Pods are containers, not VMs.** The mental model is "rebuild from an image," not "come back to my machine tomorrow." Bake the conda environment into a custom Docker image on Day 1 or you will reinstall scanpy and squidpy every morning.

## 15. Account and Pod Setup

### 15.1 Account

1. Sign up at `runpod.io`. Add credits. **You need at least one hour's worth of balance to deploy a pod**, so add more than you think.
2. Set up billing notifications in the console.
3. Generate an SSH key if you do not have one:

```bash
ssh-keygen -t ed25519 -C "sasp-kernel"
cat ~/.ssh/id_ed25519.pub
```

Paste the public key into RunPod → Settings → SSH Public Keys. This enables real SSH rather than the browser terminal.

### 15.2 Network volume (do this before the pod)

RunPod → Storage → Network Volume → Create.

- Name: `sasp-data`
- Size: **300 GB**
- Datacenter: pick one and remember it. **Pods can only mount a network volume in the same datacenter**, which constrains which GPUs are available to you. If the datacenter you pick has no cheap high-RAM pods, delete and recreate elsewhere.

At roughly $0.05/GB-month, 300 GB is about $15/month, about $6 for twelve days. It persists across pod terminations, which is the entire point.

### 15.3 Deploy the pod

RunPod → Pods → Deploy.

| Setting | Value |
|---|---|
| Cloud type | **Secure Cloud** (not Community; preemption is not worth the savings here) |
| Datacenter | Same as your network volume |
| GPU | Cheapest card with the most system RAM. Sort the list by RAM |
| Template | `runpod/pytorch` or any Ubuntu base image |
| Container disk | **20 GB** (small; billed at the higher rate) |
| Network volume | `sasp-data`, mounted at `/workspace` |
| Ports | Leave defaults; add 8888 only if you want Jupyter |

Deploy, then connect:

```bash
ssh root@<pod-ip> -p <pod-ssh-port> -i ~/.ssh/id_ed25519
```

RunPod gives you the exact command in the pod's Connect panel. The IP and port change on every deploy.

### 15.4 Sizing guidance

**RAM binds, not cores.** A Xenium 5K section with 500,000 cells has roughly 300–800 nonzero entries per cell in the count matrix, so 2–4 GB per section as a sparse matrix. Ten sections in memory is 20–40 GB before you do anything, and scanpy operations that densify a slice multiply that. **Budget 3× your raw data size.**

**Cores matter for the null battery.** Eight nulls × 1,000 permutations across tissues, cell types, and response modules is embarrassingly parallel and scales linearly with cores up to memory limits.

Strategy: a modest pod (≥64 GB RAM) for development, a larger one for the day you run the full null battery. Terminate between sessions; the network volume keeps your data.

## 16. Environment and Docker Image

### 16.1 First-time setup on the pod

```bash
apt update && apt install -y build-essential git tmux htop unzip wget

# Miniforge (conda-forge default, avoids Anaconda licensing questions)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /workspace/miniforge3
/workspace/miniforge3/bin/conda init bash
source ~/.bashrc

mamba create -p /workspace/envs/sasp python=3.11 -y
mamba activate /workspace/envs/sasp

mamba install -y -c conda-forge \
    numpy scipy pandas pyarrow \
    scanpy anndata squidpy spatialdata spatialdata-io \
    scikit-learn statsmodels \
    matplotlib seaborn \
    joblib tqdm \
    pointpats geopandas shapely \
    jupyterlab

pip install formulaic
```

**Install into `/workspace`, not the container disk.** `/workspace` is the network volume and survives pod termination. Anything in `/root` or `/` does not.

What each package buys you: `pointpats` gives Ripley's K for Test 4; `squidpy` gives spatial neighbor graphs and plotting you would otherwise write; `spatialdata-io` reads Xenium output bundles directly.

Verify:

```bash
python -c "import scanpy, squidpy, spatialdata; print(scanpy.__version__, squidpy.__version__)"
free -h && nproc && df -h /workspace
```

### 16.2 Build a custom image (Day 1, thirty minutes, saves hours)

Because pods are ephemeral, bake the environment into an image once.

`Dockerfile`:

```dockerfile
FROM mambaorg/micromamba:1.5-jammy

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git tmux htop openssh-server wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /workspace
```

`env.yaml`:

```yaml
name: base
channels: [conda-forge]
dependencies:
  - python=3.11
  - numpy
  - scipy
  - pandas
  - pyarrow
  - scanpy
  - anndata
  - squidpy
  - spatialdata
  - spatialdata-io
  - scikit-learn
  - statsmodels
  - matplotlib
  - seaborn
  - joblib
  - tqdm
  - pointpats
  - geopandas
  - shapely
  - jupyterlab
  - pip
  - pip:
      - formulaic
```

Build and push from your laptop:

```bash
docker build -t <dockerhub-user>/sasp:v1 .
docker push <dockerhub-user>/sasp:v1
```

Then specify `<dockerhub-user>/sasp:v1` as the pod template. New pods are ready in under a minute.

## 17. Storage and Data Transfer

### 17.1 Layout on the network volume

```
/workspace/
  data/
    raw/          # untouched downloads
    interim/      # cleaned AnnData, one .h5ad per section
    processed/    # scored, annotated, ready to model
  results/        # kernel fits, null distributions, tables
  code/           # git repo
  figures/
  envs/           # conda environment
```

**Put `/workspace/code` under git and push to GitHub from the pod.** If the volume is lost, data is re-downloadable; code you wrote at 2 AM is not.

### 17.2 Downloading from GEO

```bash
mkdir -p /workspace/data/raw && cd /workspace/data/raw

wget -r -np -nH --cut-dirs=5 -R "index.html*" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310392/suppl/"
```

Check the series page in a browser first to see what is deposited and how large. If the deposit is per-sample rather than series-level, iterate over GSM accessions.

**Critical download tip.** Xenium output bundles include `transcripts.parquet`, one row per detected transcript, which can be tens of gigabytes per section. **You do not need it.** You need `cell_feature_matrix.h5` (or the MEX directory) and `cells.parquet`, typically a few hundred megabytes per section together. Downloading only these turns a multi-terabyte problem into a manageable one. Only pull `transcripts.parquet` if you decide to re-segment, which you should not do in twelve days.

### 17.3 Getting results out

RunPod egress is free, which is a real advantage over AWS.

```bash
# From your laptop
scp -P <port> -r root@<pod-ip>:/workspace/results ./results
scp -P <port> -r root@<pod-ip>:/workspace/figures ./figures
```

Or push results to GitHub from the pod if they are small enough (kernel fit tables and null distributions usually are).

## 18. Making the Compute Small

Most of the estimated CPU hours disappear if the code is written correctly. Three things do almost all the work.

### 18.1 KD-tree for neighbors, always

**A 500,000-cell section has 1.25 × 10¹¹ pairwise distances, which is a terabyte in float64.** `scipy.spatial.cKDTree` turns this into seconds and a few hundred megabytes. Getting this wrong is the single most likely way to repeatedly crash a pod.

```python
from scipy.spatial import cKDTree
import numpy as np

# coords: (n_cells, 2) array in MICRONS
tree_senders = cKDTree(coords[sender_mask])

# Distance to nearest sender, for every cell
d_nearest, _ = tree_senders.query(coords, k=1)

# Local density: neighbors within 50 um
tree_all = cKDTree(coords)
density_50 = np.array([len(x) for x in tree_all.query_ball_point(coords, r=50.0)]) - 1
```

Never build an `(n, n)` distance matrix. Not once, not for a subset "just to check."

### 18.2 Vectorize the permutations

The naive null battery is 8 nulls × 1,000 permutations × tissues × cell types × modules, each with a model refit. That is 8,000+ refits and it is avoidable.

For **label permutation** nulls (N1, N7), the coordinates never change. Precompute the geometry once, and a permutation becomes re-indexing rather than re-querying the tree. For binned-mean statistics, a thousand permutations is one sparse matrix multiply against an `(n_cells × 1000)` permuted-label matrix. Seconds, not hours.

Only **torus shift** (N3) and **rotation** (N4) genuinely change the geometry and require rebuilding the tree. Those are the expensive ones. Concentrate your parallelism there.

### 18.3 Parallelize the expensive nulls

```python
from joblib import Parallel, delayed

def one_shift(seed):
    rng = np.random.default_rng(seed)
    shifted = torus_shift(sender_coords, rng, bounds)
    return fit_statistic(coords, shifted, response, covariates)

results = Parallel(n_jobs=-2, prefer="processes", batch_size=8)(
    delayed(one_shift)(s) for s in range(1000)
)
```

`n_jobs=-2` leaves one core free so the machine stays responsive. Watch memory in `htop`: joblib process workers each get a copy of anything they close over. If memory spikes, pass large arrays through a memory-mapped file rather than closing over them.

### 18.4 Checkpoint incrementally

Non-negotiable if you use Community Cloud, and good practice regardless.

```python
import os
os.makedirs("/workspace/results/nulls", exist_ok=True)

for block in range(0, 1000, 50):
    out = Parallel(n_jobs=-2)(delayed(one_shift)(s) for s in range(block, block + 50))
    np.save(f"/workspace/results/nulls/shift_{block:04d}.npy", np.array(out))
```

A preemption or crash then costs one block, not the whole run.

### 18.5 Always work inside tmux

If your SSH connection drops, anything not in tmux dies with it. Campus wifi will drop.

```bash
tmux new -s sasp        # start
# ctrl-b then d         # detach
tmux attach -t sasp     # reattach
```

For long unattended runs:

```bash
nohup python -u scripts/run_nulls.py > logs/nulls_$(date +%F_%H%M).log 2>&1 &
```

`-u` forces unbuffered output so the log updates live.

### 18.6 Profile before scaling up

Run the full pipeline on **one section** and time each step. If a step takes 40 minutes on one section it takes hours on ten, and you should fix it rather than rent a bigger pod. Renting your way out of an `O(n²)` loop does not work.

## 19. Cost Control

### Rules

1. **Terminate pods when you finish for the day.** On RunPod, terminating is correct because stopped container disk costs *more* than running container disk. The network volume keeps your data.
2. **Keep the container disk at 20 GB.** All real data goes on the network volume.
3. **Check the billing page every other day.** Two minutes; catches surprises early.
4. **Set a phone reminder** at 9 PM: "is the pod terminated." Automation is better; a reminder you actually see beats a script you will not write.

### Rough total

| Item | Estimate |
|---|---|
| Low-end GPU pod, development, ~35 h | ~$15 |
| Larger pod, pipeline + null battery, ~25 h | ~$25 |
| 300 GB network volume, 12 days | ~$6 |
| Container disk, 20 GB | ~$1 |
| Egress | $0 (free on RunPod) |
| **Subtotal** | **~$47** |
| **With 50% debugging buffer** | **~$70** |

All prices approximate as of August 17, 2026. Verify on RunPod's pricing page.

## 20. Failure Modes

| Symptom | Cause | Fix |
|---|---|---|
| Pod freezes, SSH hangs, job dies silently | Out of memory, OOM killer | `dmesg \| tail`. Chunk the data or size up |
| Environment gone after redeploy | Installed to container disk, not `/workspace` | Install to `/workspace/envs`, or use the custom image |
| Pod terminated mid-run without warning | Community Cloud preemption | Use Secure Cloud; checkpoint every 50 permutations |
| Cannot mount network volume | Pod is in a different datacenter | Deploy in the volume's datacenter, or recreate the volume |
| Job dies when laptop sleeps | Not in tmux | Use tmux |
| Storage bill higher than expected | Large container disk, or stopped-pod disk rate | Shrink container disk; terminate rather than stop |
| Memory explodes during permutations | joblib workers copying large arrays | Memory-mapped arrays, or the threading backend where work releases the GIL |
| Kernel fit is instant, neighbor computation takes an hour | Pairwise distance matrix | `cKDTree` |
| `df` shows a small disk | Writing to container disk instead of `/workspace` | Check your paths |

## 21. AWS Fallback

If RunPod does not work out, AWS EC2 is the alternative. Sections 16 through 18 transfer unchanged; only the provisioning differs.

- **Instance:** `r7i.2xlarge` (8 vCPU, 64 GB, ~$0.53/hr) for development, `r7i.4xlarge` (16 vCPU, 128 GB, ~$1.06/hr) for the pipeline, `r7i.8xlarge` for the null battery. Avoid `t3`/`t2` entirely; they throttle when CPU credits run out, which is exactly during a long permutation run.
- **Storage:** 300 GB gp3, roughly $0.08/GB-month. Do not accept the default 8 GB root volume; change it at launch.
- **Before launching anything:** set a Budgets alert at $150 with a forecasted-spend trigger, and a CloudWatch billing alarm at $50 (must be created in `us-east-1`).
- **Region:** `us-west-2` (Oregon), closest to Pasadena.
- **Security group:** SSH port 22 from **My IP** only, never `0.0.0.0/0`.
- **Do not allocate an Elastic IP.** Free while attached, billed while unattached, which is exactly when you have stopped the instance overnight.
- On AWS, **stop** rather than terminate for the day (stopped instances cost nothing for compute), the opposite of the RunPod rule.
- **Egress costs money on AWS.** Pull down figures and tables, not raw data.
- Snapshot an AMI after the environment is built so relaunching takes two minutes instead of thirty.

Two other fallbacks worth pursuing in parallel:

- **Caltech research computing.** Free, and this workload is small by HPC standards. Email on Day 0.
- **Your own machine.** With 32 GB or more of RAM you can likely do the whole project on one tissue locally. Run 200 permutations during development and scale to 1,000 only for final numbers.

---

# Part IV — Experiments

## 22. Experimental Design

Ownership tags: **[CS]**, **[Bio]**.

### Step 0 — Day 1 feasibility audit [both]

Section 8. Stop if Tests 2, 3, or 5 fail.

### Step 1 — Synthetic identifiability study [CS] — do this first

Scheduled before any real-data work because it is the ML contribution, it is fully under your control, and it produces a publishable figure even if every dataset falls through.

1. **Generate synthetic tissue.** Place N cells by a spatial point process. Designate a sender subset with controlled clustering: Poisson for random placement, through Thomas or Matérn cluster processes for increasingly aggregated placement, with a clustering parameter you sweep.
2. **Plant a known kernel.** $r_i = \mu + \beta e^{-d_i/\lambda_{\text{true}}} + \gamma^\top z_i + \varepsilon_i$ with known $\lambda_{\text{true}}$.
3. **Add realistic nuisance.** Spatial autocorrelation in the baseline (Gaussian random field), cell-type-dependent baselines, count noise, and a **sender-density-correlated confounder**. That last one is the key: make the neighborhoods around senders genuinely different for reasons unrelated to signaling.
4. **Sweep and measure recovery.** Vary sender prevalence, sender clustering, confounder strength, baseline autocorrelation length relative to $\lambda_{\text{true}}$, and N. Record bias and CI coverage of $\hat\lambda$.
5. **Report the regime map.** Expected headline: when baseline autocorrelation length is comparable to or larger than $\lambda_{\text{true}}$, the kernel is not identifiable without a matched-decoy control, and naive estimation is confidently wrong.

**This step is the paper's spine.** It tells the reader how much to believe any real-data number, including yours.

### Step 2 — Naive estimation on real tissue [CS]

Reproduce the standard analysis, deliberately without controls, so you have the thing you are going to test.

- Bin cells by distance to nearest sender (10 μm bins out to 300 μm).
- Plot mean response per bin, per receiver cell type, per tissue.
- Fit each kernel family; report $\hat\lambda_{\text{naive}}$ and model comparison.

Expected: a clean-looking monotone decay comparable to the published gradients. This is Figure 2a and it is what you are about to attack.

### Step 3 — The null battery [CS]

Each null applied independently, then in combination. Report $\hat\beta$ and $\hat\lambda$ under each.

**N1 — Cell-type-stratified label permutation.** Reassign sender labels at random among cells of the same cell type. Preserves composition and architecture; destroys sender-specific signal. 1,000 repeats.

**N2 — Matched-decoy senders (the critical control).** For each true sender, select a non-senescent cell matched on cell type, local density (k-NN within 50 μm), and k-NN composition. Recompute the kernel using decoys as senders. **Any decay appearing around decoys is not a SASP effect.** Report $\hat\beta_{\text{true}} - \hat\beta_{\text{decoy}}$ as the corrected effect. **This is the single most important number in the paper.**

**N3 — Torus shift.** Translate the sender coordinate set by a random vector with wraparound. Preserves sender clustering structure and receiver autocorrelation; destroys their alignment. 1,000 repeats. This is the null CellWHISPER showed leading methods fail.

**N4 — Rotation.** Rotate sender coordinates about the tissue centroid by random angles. Complements N3 for anisotropic tissue.

**N5 — Nuisance conditioning.** Include Tier D covariates. Report $\hat\beta$ before and after.

**N6 — Receiver baseline conditioning.** Regress out a spatially smoothed estimate of the receiver's expected response given its neighborhood *excluding* sender identity. Targets spatial autocorrelation directly.

**N7 — Sender threshold sensitivity.** Re-run at the 90th, 95th, and 99th percentile thresholds, and under both DeepScence and SenePy calls.

**N8 — Gene-set disjointness.** Confirm and report zero Tier A ∩ Tier B intersection. Additionally re-run with a scrambled response gene set matched on size and expression level.

### Step 4 — Kernel comparison and cross-tissue analysis [CS + Bio]

- Fit all kernel families under full control (N2 + N5 + N6). Report $\hat\lambda$ with donor-bootstrap CIs.
- Compare $\lambda_{\text{proximal}}$ vs $\lambda_{\text{downstream}}$.
- Compare λ across tissues and receiver cell types. Test conservation with a mixed-effects model treating tissue as a random effect: is there a common λ, or is λ tissue-specific? Senotype heterogeneity predicts the latter.
- Nearest-sender vs superposition model comparison.

### Step 5 — Biological interpretation [Bio]

- Which Tier C ligand-receptor pairs are expressed in sender and receiver populations at the fitted length scale?
- Does λ fall in the plausible range for cytokine signaling in tissue (tens of microns), or imply something else (long-range vesicle transport, systemic factors, or a missed confound)?
- Fold the estimate back into Martin et al.'s framework: which containment mechanisms remain consistent, which are excluded?

## 23. Baselines and Null Battery

Reviewers at ml4spatialbio will ask what you compare against. Reviewers at ICBINB will ask whether your nulls are strong enough.

**Method baselines. Run them, report them, do not claim to beat them:**

- **COMMOT**, **SpaTalk**, and **CellChat v2** on the same data for senescence-relevant ligand-receptor pairs. Then run each on **torus-shifted coordinates**. If they report comparable signal on shifted data, you have independently reproduced the CellWHISPER finding on senescence data specifically, which is a figure.
- **NCEM** linear variant, which gives a directly comparable niche-effect estimate with a length scale.
- **Naive binned mean**, the Zhao et al.-style analysis, as the "what the field currently does" reference.
- **Nonparametric spline**, to test whether any parametric kernel is justified.

**Do not attempt to beat SCILD or COMMOT on their own terms.** You are not proposing a better inference engine; you are asking what any of these estimates mean. State this explicitly in the Introduction so reviewers do not evaluate you on the wrong axis.

## 24. Statistical Analysis Plan

1. **Donor and section are the unit of replication, not the cell.** Bootstrap over donors. With fewer than 3 donors, label the result a case study in the figure title.
2. **Mixed-effects fits** with random intercepts for donor and section throughout.
3. **All nulls at 1,000+ permutations**, with empirical p-values and the full null distribution shown, not just a p-value.
4. **Report the surviving fraction with CIs for every null**, side by side with the naive estimate. Two numbers, always.
5. **Multiple testing:** BH-FDR across programs × cell types × tissues.
6. **Model selection:** AIC and held-out log-likelihood on left-out sections, not in-sample fit.
7. **Coverage check on synthetic data:** report whether your CIs achieve nominal coverage under the confounded regimes. Most methods in this space have never done this; it is cheap and strong.
8. **Prespecify the analysis before looking at real data.** Write the analysis script against synthetic data in Step 1, then run it once on real data. Estimating a curve from a snapshot has an enormous garden of forking paths.

## 25. Figures

Four main figures.

**Figure 1 — Identifiability regime map (synthetic).**
(a) Bias in $\hat\lambda$ as a function of sender clustering and baseline autocorrelation length. (b) CI coverage across the same grid. (c) Recovery with and without the matched-decoy control.
*Expected: naive estimation is badly biased and overconfident once baseline autocorrelation length approaches the true kernel length; the decoy control restores approximate calibration.*

**Figure 2 — The naive decay curve and what happens to it.**
(a) Binned response vs distance-to-nearest-sender, per tissue, looking exactly like the published gradients. (b) Same curve overlaid with the matched-decoy curve and the torus-shift null band. (c) Surviving fraction of $\hat\beta$ under each null.
*The money figure. Whatever the surviving fraction turns out to be, it is a result.*

**Figure 3 — Controlled kernel estimates.**
(a) $\hat\lambda$ with donor-bootstrap CIs, per program × receiver cell type × tissue. (b) Kernel family comparison. (c) Proximal vs downstream length constants. (d) Nearest-sender vs superposition comparison.

**Figure 4 — Existing tools under the same nulls.**
COMMOT, SpaTalk, and CellChat v2 signal on real vs torus-shifted coordinates for senescence-relevant pairs.
*If they fail, that is a direct senescence-specific replication of the CellWHISPER result and it justifies the whole paper. If they pass, that is also worth reporting and strengthens confidence in the published literature.*

**Supplementary:** full sweeps, all null distributions, coverage tables, gene sets and intersection matrix, sensitivity analyses, zonation-stratified fits.

---

# Part V — Execution

## 26. Twelve-Day Timeline

Today is Day 0 (August 17). Submission is Day 12 (August 29).

| Day | CS Lead | Bio Collaborator |
|---|---|---|
| **0** | Section 0.1 action list. RunPod account, network volume, Docker image build | Confirm availability in writing. Begin dataset scouting |
| **1** | Deploy pod. Set up spatial stack. Write the point-process simulator | Dataset audit (Deliverable 1). Panel and disjointness checks (Test 2). Fill audit table |
| **2** | Synthetic generator complete: clustering sweep, autocorrelated baseline, confounder injection. Run Tests 1, 4, 5, 6 | Gene set package (Deliverable 2). Score senders; check prevalence (Test 3) |
| **3** | **GATE:** synthetic recovery working, naive estimator implemented, audit passed. If the audit failed, proceed synthetic-only and reframe (Section 28) | Anatomical annotation (Deliverable 3). Deliver cleaned dataset with sender labels. Start Related Work draft |
| **4** | Run the full synthetic sweep. **Figure 1** | Sender call validation (Deliverable 4). Finish Related Work. Draft Introduction |
| **5** | Naive estimation on real data, all tissues and cell types. **Figure 2a** | Sanity-check sender calls against images. Flag implausible cell types |
| **6** | **GATE:** Figures 1 and 2a exist. Implement null battery N1–N8. Matched-decoy matching is the fiddly one; budget the day | Begin ligand-receptor plausibility (Deliverable 5) |
| **7** | Run all nulls, 1,000 permutations each. **Figures 2b, 2c** | Continue interpretation. Draft the biology results narrative |
| **8** | Controlled kernel fits, family comparison, proximal vs downstream, cross-tissue mixed model, zonation-stratified fits. **Figure 3** | Containment interpretation (Deliverable 6) |
| **9** | Run COMMOT / SpaTalk / CellChat on real and shifted coordinates. **Figure 4** | Draft Discussion and Limitations |
| **10** | All statistics, coverage checks, final figures. **Choose which abstract version to use** | Full editing pass on biological claims |
| **11** | Write Methods and Results. Code cleanup, seed pinning, README | Write Abstract. Assemble supplementary |
| **12** | Full read-through. Fix the one obviously wrong thing. **Submit by 11:59 PM AoE** | Claim audit (Deliverable 7). Verify every citation resolves |

**Critical path:** Day 1 audit → Day 3 gate → Day 6 null battery → everything.

**Non-negotiables:**
- Start writing on Day 3, not Day 10.
- Write **two abstracts on Day 3**, one for "the effect largely survives" and one for "it largely does not." Choose on Day 10. This is the single best insurance against a bad result costing you the deadline.
- **Figure 1 is synthetic and depends on nothing external.** If real data collapses on Day 3, the paper becomes "Identifiability limits for spatial response kernel estimation, with senescence as the motivating case," which is still a legitimate ml4spatialbio submission. Do not skip Step 1 to save time; it is what makes the deadline survivable.

## 27. Budget

| Item | Estimate |
|---|---|
| RunPod compute (Section 19) | ~$70 |
| Network volume, 12 days | included above |
| Egress | $0 |
| **Total** | **~$70** |

CPU-bound, not GPU-bound. This is genuinely a cheap project. If Caltech HPC comes through, it is free.

## 28. Risks and Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| **No spatial dataset with both single-cell resolution and a usable disjoint response panel** | **High.** Most likely failure mode; targeted panels are small and often lack downstream programs | Day 1 audit, Test 2. GSE310392's 5K panel is specifically why it ranks first. Fallback: Visium HD with deconvolution and a declared limitation, or synthetic-only reframe |
| **The naive decay does not survive the nulls at all** | Medium-High | **This is a result, not a failure.** It is the ICBINB paper. Two abstracts on Day 3 |
| **The decay survives everything and looks clean** | Low-Medium | Also a result, and a more comfortable one: the first properly controlled estimate of senescence reach. Report λ with CIs and cross-tissue comparison |
| **Sender-score circularity** | High if careless | Enforce and report disjointness (Test 2, N8). Run the scrambled-response control. A reviewer will check this |
| **Zonation or another anatomical axis drives the whole effect** | **High for liver specifically** | Section 11. Include the covariate, run zonation-matched decoys, stratify within zones, and make the with/without comparison a figure panel |
| **Senders confined to one anatomical compartment** | Medium | Detected by Test 6 on Day 1. If true, add the covariate; if the covariate absorbs everything, drop that tissue |
| **Too few donors for valid uncertainty** | Medium | GSE310392's 43 donors is why it ranks first. Otherwise report section-level bootstrap and label it a case study |
| **Twelve days is not enough** | **High** | Day 3 and Day 6 gates are real. If Figure 1 is not producible by Day 4, stop and target ICLR 2027 workshops or a spring venue with the full version |
| **Reviewers read it as "yet another spatial CCC method"** | Medium | Frame as evaluation and identifiability from the first sentence of the abstract. State explicitly you are not proposing to beat SCILD or COMMOT |
| **Pod preemption during the null battery** | Low on Secure Cloud | Secure Cloud plus checkpointing every 50 permutations (Section 18.4) |

## 29. Positioning and Reviewer Objections

### Venues

**Primary: ml4spatialbio (Paris, December).** Their call lists modeling cell-cell communication and tissue dynamics, building interpretable and uncertainty-aware spatial models biologists can trust, and designing benchmarks and evaluation standards specific to spatial tasks. Their overview says outright that spatial methods remain immature and there is little agreement on how to evaluate spatial models. This project is written to that sentence. Deadline August 29 AoE, notification September 29, submissions via OpenReview.

**Secondary: ICBINB-BIO (Sydney).** Topics include out-of-distribution generalization and domain shift, causal mechanisms versus spurious correlation, and interpretability and trustworthy biological inference. If the nulls kill the effect, submit here instead; the venue exists for exactly that result.

### Anticipated objections and responses

1. **"Zhao et al. already showed the distance gradient in *Cell*."** Yes, descriptively, without negative controls or uncertainty. We reproduce their qualitative finding and then quantify how much survives matched-decoy and torus-shift nulls. That is a different claim.
2. **"SCILD, COMMOT, and NCEM already model distance decay."** They do, and we run them as baselines. Our question is whether the estimates any of them produce are identifiable from a static snapshot given sender clustering. We evaluate the estimand, not compete on the estimator.
3. **"Your senescence calls are inferred, not experimental."** True and unavoidable in tissue. We report sensitivity across thresholds and across two independent signature sources and show the conclusions do not depend on the choice. If they do, we say so.
4. **"Only one or two tissues."** Twelve-day workshop paper, and the synthetic study covers the generality question that more tissues would only partly address.
5. **"Why not a time course or a perturbation?"** Because that data does not exist publicly at spatial resolution. We name it as the decisive future experiment: induce senescence at a known focus, profile spatially at multiple timepoints, and fit the kernel with the confound removed by design.
6. **"Isn't the decoy control too conservative?"** Possibly, and we quantify it: on synthetic data with a known planted kernel, we report how much true signal the decoy control removes. That number is in Figure 1c.

## 30. Paper Outline

Target 8 pages plus appendix; confirm against the workshop template on Day 0.

1. **Abstract.** Two versions written Day 3, chosen Day 10.
2. **Introduction.** Senescence spreads locally; how far is a central open question with therapeutic stakes. Existing evidence is descriptive. Contemporaneous benchmarks show spatial CCC inference is poorly calibrated. We ask whether the senescence kernel is estimable at all, and what it is once controlled.
3. **Related Work.** Spatial senescence mapping (SenNet review, Zhao et al., scDOT, Farzad et al., DeepScence). Spatial CCC methods (NCEM, COMMOT, SCILD, RGAST, SpaTalk, HARMONIC). Reliability and benchmarking (CellWHISPER, CONCISE, ctSVG benchmark). Mechanistic models of senescence spread (Martin et al., Boolean SASP models).
4. **Methods.** Kernel model and families; sender and response definitions with disjointness; null battery; synthetic generator; statistics.
5. **Results.**
   - 5.1 Identifiability regimes on synthetic tissue (Fig. 1)
   - 5.2 The naive gradient and its behavior under nulls (Fig. 2)
   - 5.3 Controlled kernel estimates across tissues and programs (Fig. 3)
   - 5.4 Existing tools under the same nulls (Fig. 4)
6. **Discussion.** What the surviving estimate implies for senescence containment. What it implies for how spatial senescence claims should be reported. Limitations: inferred senescence calls, static snapshots, targeted panels, one or two tissues, no perturbation. The decisive future experiment.
7. **Appendix.** Full sweeps, null distributions, coverage tables, gene sets, sensitivity analyses.

## 31. Key References

Verify every bibliographic detail before submission. DOIs given where confirmed as of August 17, 2026; several years and page numbers should be double-checked against the publisher record.

**Senescence in space**
1. Zhao et al. (2024). Spatial transcriptomic landscape unveils immunoglobin-associated senescence as a hallmark of aging. *Cell*. S0092-8674(24)01201-7 — **closest prior work**
2. Spatial mapping of cellular senescence: emerging challenges and opportunities (2023). *Nature Aging*. doi:10.1038/s43587-023-00446-6
3. Cellular senescence in human liver under normal aging and cancer (2026). *Cell Genomics*. S2666-979X(25)00389-1; GEO GSE310392 — **your primary dataset**
4. scDOT: optimal transport for mapping senescent cells in spatial transcriptomics (2024). *Genome Biology*. doi:10.1186/s13059-024-03426-0
5. Farzad et al. (2026). A spatial multi-omics atlas of immunosenescence reveals germinal-center B cell alteration in human lymph nodes. *Cell Press Blue* 1(4):100053
6. Suryadevara, Farzad et al. (2026). Charting human cellular senescence in aging and disease. *Cell*. S0092-8674(26)00587-8
7. Spatial distribution of senescent cells and their proximity to immune subsets in the human endometrium (2026) — source of the 45–211 μm nearest-neighbor calibration figures

**Senescence detection**
8. Qu et al. (2025). Single-cell and spatial detection of senescent cells using DeepScence. *Cell Genomics* 5(12):101035. doi:10.1016/j.xgen.2025.101035 — **your sender-calling tool**
9. Anerillas et al. (2026). SenCat. *Molecular Cell*. doi:10.1016/j.molcel.2026.05.017
10. SenePy (2025). Unveiling the cell-type-specific landscape of cellular senescence through single-cell transcriptomics. *Nature Communications*. PMID 39987255
11. Saul et al. (2022). SenMayo: a new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nature Communications*. doi:10.1038/s41467-022-32552-1
12. SenCID (2024). *Cell Metabolism*. S1550-4131(24)00088-3
13. Hughes et al. (2025). SenPred. *Genome Medicine*. doi:10.1186/s13073-024-01418-0
14. Ntintas et al. (2026). Overview of molecular signatures of senescence and associated resources: pros and cons. *FEBS Open Bio*. doi:10.1002/2211-5463.70134

**SASP and paracrine spread**
15. Martin L, Schumacher L, Chandra T (2023). Modelling the dynamics of senescence spread. *Aging Cell* 22(8):e13892. doi:10.1111/acel.13892; PMC10410058 — **the containment paradox** (the paradox is theirs; the four containment mechanisms often attached to it here are not — see §3 correction)
16. Basisty et al. (2020). A proteomic atlas of senescence-associated secretomes. *PLoS Biology* 18(1):e3000599
17. Acosta et al. (2013). Paracrine senescence transmission
18. Characterizing the SASP-dependent paracrine spreading of senescence between human brain cell types (2026). bioRxiv, February 10, 2026 — CCL2, CXCR7, DPP4
19. A model of the onset of the SASP after DNA-damage-induced senescence (2017). *PLOS Computational Biology*. doi:10.1371/journal.pcbi.1005741
20. Dissecting the heterogeneity of senescence: primary and secondary senescent states. PMC11689308

**Spatial CCC methods**
21. Fischer et al. NCEM: Modeling intercellular communication in tissues using spatial graphs of cells. *Nature Biotechnology*. doi:10.1038/s41587-022-01467-z
22. Cang et al. COMMOT: Screening cell-cell communication in spatial transcriptomics via collective optimal transport. *Nature Methods*
23. SCILD (2026). Advancing spatial cellular communication inference with ligand diffusion and transport model. *Communications Biology*. doi:10.1038/s42003-025-09413-w
24. SpaTalk. Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data. *Nature Communications*
25. RGAST: A relational graph attention network for multi-scale cell-cell communication inference. bioRxiv 2024
26. HARMONIC: Histology-aware graph for modeling intercellular communication in spatial transcriptomics. bioRxiv, January 2026

**Reliability, calibration, identifiability — the framing citations**
27. CellWHISPER: Inference of direct cell-cell communication from spatial transcriptomics. bioRxiv, January 2026 — **FPR >90% for CellChat v2, COMMOT, SpaTalk under randomization**
28. CONCISE: Spatial co-expression and cell-cell communication inference from spatially resolved transcriptomics. bioRxiv, June 2026
29. Benchmarking cell-type-specific spatially variable gene detection methods (2026). *Briefings in Bioinformatics* 27(2):bbag190
30. Identifiability limits of physics-informed inference for spatial stochastic dynamics from static snapshots (2026). arXiv 2607.01749

## 32. Reading List

Six papers, in this order, roughly one working day. Do not skip this to start coding sooner; three of them will change your design.

1. **Zhao et al. (2024), *Cell*.** The closest prior result. Understand exactly what they claimed and how they measured it, because you are re-examining it.
2. **"Cellular senescence in human liver under normal aging and cancer" (2026), *Cell Genomics*.** Your Rank 1 dataset. Read the spatial sections carefully, especially the clustered-versus-isolated senescent hepatocyte analysis and the microenvironment comparison that fell short of significance.
3. **CellWHISPER (bioRxiv, Jan 2026).** The >90% FPR benchmark. Your framing citation and the source of the torus-shift null.
4. **CONCISE (bioRxiv, June 2026).** Spatial autocorrelation inflating type I error. Second framing citation.
5. **DeepScence (*Cell Genomics*, Dec 2025).** Your sender-calling tool, the CoreScence gene set, and a useful survey of the nine published senescence gene sets and their disagreement.
6. **Martin et al. (2023), *Aging Cell*.** The containment paradox your measurement speaks to.

Optional: the 2023 *Nature Aging* SenNet review for the field's own framing, and the 2026 *FEBS Open Bio* signature comparison for the gene-set heterogeneity numbers you will cite in Limitations.
