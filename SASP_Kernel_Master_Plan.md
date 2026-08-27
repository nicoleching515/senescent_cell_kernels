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
5. **Read Ma et al. (2024), *Cell*.** It is the closest prior work and it will change how you frame the paper. (Section 32) [Corrected 2026-08-27: read "Zhao et al."; the first author is **Ma S** and Zhao L is 14th of 47 — see §4.1.]
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
| Distance-dependent SASP gradients around senescent cells in tissue | **Ma et al.**, *Cell* 2024 [corrected 2026-08-27 from "Zhao et al."; see §4.1]. Defined "senescence-sensitive spots," ranked genes by proximity, showed SASP scores rise and TNF signaling falls with distance, across multiple aged mouse organs | **Done, descriptively.** This is in *Cell*. You cannot claim the observation |
| Mapping senescent cells and their neighbors in spatial transcriptomics | scDOT, *Genome Biology* 2024 (Bar-Joseph lab). Optimal transport plus deconvolution; spatial organization of senescent cells in lung, candidate interaction genes | Done |
| Explicit ligand diffusion, competitive binding, and concentration decay in a fitted model | SCILD, *Communications Biology* 2026. Single-cell resolution, interpretable optimization, in silico perturbation | Done, generically |
| Learning niche-composition effects on expression at characteristic length scales | NCEM, *Nature Biotechnology*. GNN; recovered length scales matching known communication mechanisms | Done, generically |
| Distance-weighted attention for multi-scale communication | RGAST, bioRxiv 2024 | Done |
| Mathematical model of senescence spread and the containment paradox | Martin et al., *Aging Cell* 2023 | Done; the paradox is left open, and the paper declines to say whether spread is finite |
| Classifying primary vs secondary (bystander) senescence with ML | PMC11689308; paracrine spread between human brain cell types, bioRxiv Feb 2026 | Done |
| Boolean network model of SASP onset after DNA damage | *PLOS Comput Biol* 2017 | Done |
| Senescent cell detection in spatial data | DeepScence, *Cell Genomics* Dec 2025 | Done, and you will use it |

### 2.2 Left

1. **No confounder-controlled estimate.** The *Cell* 2024 gradient analysis has no negative controls, does not test against matched non-senescent decoys, does not control local density or cell-type composition, and reports no length constant with uncertainty. Given that contemporaneous benchmarks report comparable interaction counts on permuted coordinates (implying FPR >90%), the prior probability that a naive gradient survives a proper null is not high.
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

The 2023 *Nature Aging* SenNet review, "Spatial mapping of cellular senescence," is the field's own statement of the problem and is your framing citation. The June 2026 SenNet package across Cell Press titles includes a spatial multi-omics atlas of human lymph nodes across ages 18–86 (Farzad et al., ***Cell Press Blue*** 1(4):100053 — **not *Cell*** itself; the *Cell* item in that package is the Suryadevara et al. **Review**), which found senescent-like B cells shifting from interfollicular zones into germinal centers with age. The SenNet portal hosted **2,041** public human and mouse datasets across 15 organs and 6 general assay types **as of April 2026** (Börner et al., *SenNet Portal: Build, Optimization and Usage*, bioRxiv 2026.02.06.704469).

> **Correction 2026-08-27 (citation audit CIT-3/CIT-4).** This paragraph previously said "1,753 … as of January 2026" and attributed Farzad et al. to *Cell*. The 1,753/January figure **could not be verified**: the abstracts of both v1 (2026-02-10) and v2 (2026-05-04) of the portal preprint state 2,041 as of April 2026, and no retrievable version states 1,753. The ages-18–86 and germinal-centre claims **are** verified against the published Farzad et al. article.

**Ma et al. (*Cell* 2024) is the closest scientific precedent and you must engage with it directly.** [Corrected 2026-08-21, D7 §B1/B2/B3: cited throughout earlier drafts as "Zhao et al."; the first author is **Ma S** and Zhao L is 14th of 47. The platform is **Stereo-seq, spot/bin level, not single-cell**, and the design is **two age groups (young, old)**, not a time course.] They profiled **young versus old male mice across nine tissues** on Stereo-seq bins, defined senescence-sensitive spots, and showed that SASP score, TNF signaling, ATP biosynthesis, and cell-cycle genes all vary monotonically with distance from those spots, consistently across organs. Their conclusion, that senescent foci act as epicenters compromising surrounding cells in a distance-dependent manner, is exactly the phenomenon you propose to quantify. Your contribution is quantification with controls, not discovery.

### 4.2 Spatial CCC methods

NCEM (*Nature Biotechnology*), COMMOT (*Nature Methods*), SpaTalk, CellChat v2, GCNG, DT-CCC, GraphST, NICHES, RGAST, and SCILD (*Communications Biology* 2026, which fits explicit ligand diffusion, competitive binding, and concentration decay in one optimization). HARMONIC (bioRxiv Jan 2026) adds H&E histology to condition on tissue context and reduce false positives.

### 4.3 The reliability problem, which is why this project exists

- **CellWHISPER** (bioRxiv Jan 2026) benchmarked CellChat v2, COMMOT, and SpaTalk against a null that **permutes cell locations within each cell type** — preserving cell-type expression and cell-type-specific spatial organization, destroying only ligand⁺/receptor⁺ proximity — and found they report similar interaction counts on real and permuted input, **implying** false positive rates above 90%. Their own confounder-aware null brought this under 5%. [Corrected 2026-08-27, citation audit CIT-1/CIT-5: this was written as "coordinate-randomized data", which understates their design; and the >90% figure is an **interaction-count ratio from which the authors infer an FPR**, not a measured type-I error against a nominal level. Their null is this project's **N1**, not N3 — see §22 Step 3.]
- **CONCISE** (bioRxiv June 2026) showed that introducing even weak spatial autocorrelation (a = 0.1) into one gene inflated type I error for every competitor it tested. [Precision, 2026-08-27 citation audit, per D7 §S16: those competitors are **MERINGUE, SpatialDM, Copulacci and LIANA+** — spatial co-expression / bivariate methods. **Do not use CONCISE to support a claim about the CellChat/COMMOT/SpaTalk family**; that is CellWHISPER's claim, not theirs.]
- A March 2026 *Briefings in Bioinformatics* benchmark of cell-type-specific spatially variable gene methods found rotation invariance unresolved. [Corrected 2026-08-27, per D7 §S17: the clause "and cell-type confounding a central open challenge" **is not the authors' framing** — cell-type confounding is a documented weakness of **one** method (Celina) in that benchmark, not a stated open challenge for the field. Rotation invariance is exact.]

The field has, in the last eight months, started admitting that spatial CCC inference is poorly calibrated. Nobody has yet asked what that means for the specific claims the senescence field has built on top of these tools.

## 5. The Gap We Fill

> The claim that senescent cells influence their neighbors in a distance-dependent manner is central to senescence biology and has been shown descriptively. It has never been estimated as a parameter with uncertainty, never tested against confounder-aware nulls, and never subjected to an identifiability analysis, despite contemporaneous evidence that the class of methods used to make such claims returns comparable interaction counts on permuted coordinates, implying false-positive rates above 90%.

*(Wording tightened 2026-08-27, citation audit CIT-5: "has false-positive rates above 90% under randomization" asserted a measured type-I error. CellWHISPER's number is inferred from a count ratio. Keep the "implying" construction everywhere.)*

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

**N2 — Matched-decoy senders.** For each true sender, select a non-senescent cell matched on cell type, local density (k-NN within 50 μm), and k-NN composition. Recompute the kernel using decoys as senders. **Any decay appearing around decoys is not a SASP effect** — but the converse does not follow, and that is the finding. Report $\hat\beta_{\text{true}} - \hat\beta_{\text{decoy}}$ **beside**, never instead of, the covariate-adjusted estimate.

> **Correction 2026-08-27 (plan update D12/D13, PI-decision series).** This paragraph read *"the critical control"* and closed *"**This is the single most important number in the paper.**"* **The decoy contrast is still central to the paper — as a negative result about the method, not as the primary estimate.** Three independent measurements now say that a matched decoy set removes almost nothing of what the *same variables entered as covariates* remove.
>
> - **Composition-matched decoys (D15 protocol, 5 seeds).** The protocol reproduces **98.4 %** of the naive pooled amplitude — SF **0.9837** [0.973, 0.994], i.e. **1.6 % removed** — while cell-type intercepts plus the 20-NN composition vector *as covariates on the same fits* remove **85.4 %** (SF **0.1461** [0.052, 0.246]). A factor of fifty from the same variables. Matching is not the problem: max \|SMD\| goes 0.092 → 0.035, median match rate 0.99987, and the §8 Test 5 gate passes in 100 % of matches. `results/phase3/compmatch_reruns.csv`; `reports/CS_PHASE8_COMPMATCH.md` §0.
> - **The assay's own negative-control features (A7).** Naive gradient **−0.0744** SD [−0.1306, −0.0182], p = 0.0145 on the pooled control set; **N5 removes it** (+0.0038, p = 0.72; with N6, +0.0053, p = 0.60); **N2 does not** (−0.0642, p = 0.0124 — 86 % undiminished). `results/phase3/a7_summary.csv`.
> - **Our own synthetic study, which already said so.** In Figure 1a/1b the matched-decoy design tracks the naive fit rather than correcting it: \|relative bias\| in $\hat\lambda$ is smaller than naive in only 12 of 20 grid cells and *larger* in 8, worst-case 2.27 against naive 2.03 and nuisance-conditioned 0.33; over the eight cells where baseline autocorrelation length is ≥ 2λ, mean CI coverage is **0.51 naive, 0.35 matched-decoy, 0.85 nuisance-conditioned** against a nominal 0.95. `figures/figure1_data.csv`.
>
> **The mechanism, in one sentence, and it is what makes this transferable:** matching balances the covariates between senders and decoys; it does not remove the dependence of the *response* on those covariates at the receiver, which is where the confounding acts (`reports/CS_PHASE3.md` §5). On A7 the confounder is named — per-cell detection efficiency, which N5 models directly and a propensity match on neighbourhood covariates cannot see.
>
> **Consequences, both pre-registered.** The primary controlled estimate is the **N2 + N5 + N6** fit, and **no naive or N2-only kernel may be reported on this platform**; wherever the matched SF 0.9837 appears, `type_adj` (65.9 %) and `typecomp_adj` (85.4 %) appear beside it (`reports/PREREG_PHASE8.md` §10 rule 8).

**N3 — Torus shift.** Translate the sender coordinate set by a random vector with wraparound. Preserves sender clustering structure and receiver autocorrelation; destroys their alignment. 1,000 repeats. The random-shift/toroidal null originates with **Lotwick & Silverman (1982)**; **Mrkvička et al. (2021)**, *Spatial Statistics* 42:100430, show that the torus correction "makes a crack in the autocorrelation structure" and is therefore **liberal**, that it "assumes a rectangular observation window", and that extending it to "windows which are finite unions of (aligned) rectangles" — i.e. our tiling — "would increase the amount of cracks in the autocorrelation structure [and] increase the liberality of the test." Their **variance correction**, which computes the statistic on $W \cap (W+v)$ and standardizes because different shifts retain different amounts of data, "can be applied in case of general (compact) observation windows" and is the classical remedy for our situation. Run it as an N3 variant, or state why tiling is preferred over it.

> **Correction 2026-08-27 (citation audit CIT-1).** This paragraph previously read "This is the null CellWHISPER showed leading methods fail." **It is not.** The CellWHISPER v1 full text was retrieved and read: their null "shuffles locations among cells of same type", the strings *torus*, *toroidal* and *wraparound* do not occur anywhere in the paper, and it uses no negative control probes. Their null is **N1**, not N3. N3's correct attribution is Lotwick & Silverman (1982) for the construction and Mrkvička et al. (2021) for its failure mode on non-rectangular windows.

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

- **COMMOT**, **SpaTalk**, and **CellChat v2** on the same data for senescence-relevant ligand-receptor pairs. Then run each on **torus-shifted coordinates** *and* on **CellWHISPER's own null — locations permuted within cell type**. If they report comparable signal on either, that is a figure. [Corrected 2026-08-27, citation audit CIT-1: this read "you have independently reproduced the CellWHISPER finding" of the torus-shifted run. **The torus-shifted run reproduces nothing of CellWHISPER's** — they never ran a torus shift. Only the within-cell-type permutation (`N0_type`, run in CS_PHASE4 §8) is a replication of their design.]
- **NCEM** linear variant, which gives a directly comparable niche-effect estimate with a length scale.
- **Naive binned mean**, the Ma et al.-style analysis, as the "what the field currently does" reference. [Corrected 2026-08-27: "Zhao et al." — the first author is Ma S; see §4.1.]
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
> **Correction 2026-08-27 (plan update D12/D13).** The second clause is **contradicted by this figure's own committed data.** In `figures/figure1_data.csv` the matched-decoy design's \|relative bias\| in $\hat\lambda$ is smaller than the naive fit's in only **12 of 20** grid cells and *larger* in 8 (worst case 2.27 vs naive 2.03), and over the eight cells where baseline autocorrelation length is ≥ 2λ mean CI coverage is **0.51 naive, 0.35 matched-decoy, 0.85 nuisance-conditioned** against a nominal 0.95. The clause that survives is: *nuisance conditioning restores approximate calibration; the decoy control does not.* This is the synthetic-ground-truth arm of the N2-vs-N5 result (§22 Step 3 N2, §29 objection 6, §30 5.6). Figure 1 itself was not regenerated at task 8.7 and does not need to be — the panels already show this; the caption and the surrounding text do not.

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

#### ⚠ Venue reassessment — evidence for the PI, no decision taken

*(Added 2026-08-27, plan update D12/D13. **The assignment above is unchanged and remains the PI's call.** This subsection exists because the nulls did in fact kill the effect, which is the condition the "Secondary" paragraph itself names, and because an independent novelty review argues the primary/secondary assignment is now backwards. Source: `reports/NOVELTY_ASSESSMENT.md` §5, whose venue facts were retrieved from the two CfP pages on 2026-08-27, not recalled.)*

**The facts, as retrieved.**

| | ml4spatialbio 2026 (Paris, Dec) | ICBINB-BIO 2026 (Sydney, Dec 11–12) |
|---|---|---|
| Length | **4 pages** max, excl. references + optional appendix | **8 pages** (full track); 4-page "tiny" track also offered |
| Archival | **Non-archival; concurrent submission explicitly permitted** | not verified |
| Deadline | CfP page now shows a **window, Aug 29 – Sept 4 AoE**, not the hard Aug 29 recorded above | **Aug 29, 2026, 11:59 AoE** |
| Topic match | cell-cell communication; interpretable and uncertainty-aware spatial models; benchmarks and evaluation standards for spatial tasks. **Negative results not mentioned** | causal mechanisms vs spurious correlation; failure under weak or confounded supervision; uncertainty, calibration and decision-aware reliability; explicitly wants "candid failure analysis, negative results … and benchmarks that support the claims" |

**The argument for inverting.** §29's assignment was made before the outcome was known and was right for "the effect largely survives". The outcome is naive **0.329** → controlled **0.029**, nothing above a **0.183** detectable bound. ICBINB-BIO then matches topic-by-topic: a confounded-supervision failure, a spurious-correlation-versus-causal-mechanism story, a calibration measurement (9–16 % against 5 % nominal), and a candid self-correction of the project's own motivating claim. It also gives 8 pages, which is the only length at which the null battery, the identifiability study, the N2-vs-N5 result and Figures 1 and 4 all fit — and §2.5 of the novelty review is blunt that with **no published length constant for senescence spatial influence to contradict**, the paper's force rests entirely on Figures 1 and 4, so squeezing them is the one thing that cannot be done.

**The argument for keeping ml4spatialbio in play regardless.** It is non-archival and permits concurrent submission, so a 4-page distillation there is compatible with an 8-page submission elsewhere *on ml4spatialbio's side*. The right 4-page paper for that venue is not this one compressed: it is **"negative-control probes as a calibration instrument for distance-dependent spatial models"** — the A7 construction, the N2-vs-N5 result, the measured FPR and the torus/tiling methods note — which lands on "benchmarks, datasets and evaluation standards specific to spatial tasks" and is positive-valence (here is a tool) even though the underlying result is negative.

**Two things to check before acting, both unresolved here.** (i) **ICBINB-BIO's dual-submission policy was not verified** — ml4spatialbio's permission is one-sided evidence. (ii) The Sept 4 ml4spatialbio window has to hold, and if it does it gives six days after the ICBINB deadline to cut the short version from the long one, which is strictly easier than the reverse ordering.

### Anticipated objections and responses

1. **"Ma et al. already showed the distance gradient in *Cell*."** Yes, descriptively, without negative controls or uncertainty, and on Stereo-seq bins rather than segmented cells. We reproduce their qualitative finding and then quantify how much survives the full battery — nuisance and receiver-baseline conditioning (N5, N6) and the **variance-corrected** random shift (N3-var), with the matched-decoy contrast reported alongside as a diagnostic rather than as the correction (§22 Step 3 N2, and objection 6 below). *(Wording updated 2026-08-27: this read "how much survives matched-decoy and torus-shift nulls", both of which are now the wrong primaries.)* The answer is **0.088** [−0.017, 0.234] of the naive amplitude, i.e. a controlled **0.029** [−0.007, 0.084] response-sd against a detectable bound of **0.183** at 80 % power over 0–100 µm (`results/phase3/m1_final_audit.txt`). That is a different claim, and the resolution argument is the one to lead with: their measurement is at Stereo-seq bin level, which blurs the very quantity being estimated. [Corrected 2026-08-27: this objection was written as "Zhao et al."; the first author is **Ma S** and Zhao L is 14th of 47.]
2. **"SCILD, COMMOT, and NCEM already model distance decay."** They do, and we run them as baselines. Our question is whether the estimates any of them produce are identifiable from a static snapshot given sender clustering. We evaluate the estimand, not compete on the estimator.
3. **"Your senescence calls are inferred, not experimental."** True and unavoidable in tissue. We report sensitivity across thresholds (N7) and across signature sources — and the closing clause this response used to carry, *"if they do, we say so"*, is now the operative half of it. **They partly do.** *(Rewritten 2026-08-27, plan update D12/D13, once tasks 8.4/8.5/8.7 reported. The struck version claimed "two independent signature sources" and "the conclusions do not depend on the choice"; the first word is false and the second needs splitting in two.)* Four things go in the paper rather than waiting to be found:
   - **The callers are not independent, and the published sentence was wrong on the published data.** Depth- and cell-type-matched, the top-5 % calls of DeepScence, SenePy, a disjoint arrest-and-damage score and `Cdkn1a`⁺ overlap at **0.751–2.198× of chance** across eleven mouse liver sections, pooling to **1.212×** (Mantel–Haenszel, z = 21.9, p = 1.8 × 10⁻¹⁰⁶). On the **published two-section base**, with coverage held fixed, replacing the hollow 25-gene Tier A with the frozen 33-gene one already moves the three-pair pooled ratio from **1.030 (p = 0.20)** to **1.128 (p = 4.4 × 10⁻⁸)**. Independence dies from the sender-set fix; coverage then makes it certain. `results/phase3/caller_coverage_gate_headline.csv`, `caller_agreement_matched_significance_{verify2sec,2sec_c6,11sections}.csv`; `reports/SUBMISSION_PATCH_2026-08-29.md` §3.1.
   - **The dependence is weak, heterogeneous, and mechanistically different at each end** — 0.737× (SenePy vs DeepScence, z = −15.1) to 1.471× (arrest vs `Cdkn1a`⁺, above chance 11/11, whose cause is the p53 axis and not shared membership) — which is *not* what one latent state looks like. Do **not** ship the older "weak and technical rather than latent, and each pair's direction is predicted by its depth loading" wording: the depth-direction rule is refuted at the level of independent units (pair-level exact permutation p = 0.30; within-pair Spearman negative in 5 of 5). `SUBMISSION_PATCH_2026-08-29.md` §2.
   - **DeepScence's polarity is a `CDKN1A` anchor, and the anchor is weak or reversed in four of eleven sections** once transcript depth is partialled out (ρ spans −0.024 to +0.182; fold-split sign stability falls to 0.60 in one section). Two caller-free anchors — an eight-gene proliferation set and `Lmnb1` — agree with each other in 11 of 11 and with the published anchor in 10 of 11; the exception is 7250, which was half of the published two-section base. Frame this as a documented consequence of DeepScence's own sign-anchoring step and a cross-section comparability caveat for its users, **not** as a defect we discovered. `reports/CS_PHASE8_CALLERS.md` §5.3, `results/phase3/deepscence_anchor_decisions.csv`.
   - **"Two independent signature sources" is not available on the human arm at all.** SenePy ships no spleen signature: of 22 spleen cell types, 0 are matched, 15 take a surrogate and 7 are unscoreable (`results/phase7_jobA/senepy_spleen_coverage.csv`). It is a declared deviation with no fix.

   **What survives, and it is the sentence to lead with:** the *conclusion* does not depend on the caller even though the *agreement claim* did. Across seven sender definitions spanning 0.50–8.96 % prevalence the in-tissue correction is unmoved (N3-tile **0.960–1.001**; `results/phase3/sf_summary_c1_n7.csv`, `m1_n7_prepost.txt`), and the controlled amplitude is **0.029** [−0.007, 0.084] against a detectable bound of **0.183** at 80 % power. Say "the bound is invariant to the sender definition; our caller-independence claim was not, and we withdraw it" — never "the sources agree."
4. **"Only one or two tissues."** Workshop paper, and the synthetic study covers the generality question that more tissues would only partly address. *(Updated 2026-08-27.)* Concretely: the mouse arm M1 is **one tissue, eleven liver sections, six admissible under §8 Test 3**, and the human arm H1 is **one tissue, seven donor spleens** (GSE326743, Xenium Prime 5K; acquired and structurally verified, but it runs in Phase 10 behind the freeze tag and no H1 estimate exists yet) which is a replication of the *geometry*, not an ageing result and not a second liver — per PI decision D4, age enters as a continuous covariate only and no age-stratified claim is made. Two species and two tissues are confounded with each other by design and we say so rather than implying a cross-tissue λ comparison we cannot make (`reports/PHASE7_H1_SCREEN.md`; addendum §18 outcome B).
5. **"Why not a time course or a perturbation?"** Because that data does not exist publicly at spatial resolution. We name it as the decisive future experiment: induce senescence at a known focus, profile spatially at multiple timepoints, and fit the kernel with the confound removed by design.
6. **"Isn't the decoy control too conservative?"** **No — on this assay it is insufficient, and that inversion is one of the paper's results.** *(Rewritten 2026-08-27, plan update D12/D13. The struck version answered "possibly, and we quantify it … that number is in Figure 1c", which assumed the answer was "yes". Figure 1c is still the right panel for how much *planted* signal a decoy contrast removes; it is the wrong panel for the question a reviewer is actually asking, which is how much *confounding* it removes.)* Two independent instruments, plus our own synthetic grid, say a matched decoy set fails to remove what covariate adjustment removes:
   - **A7, the negative-control-outcome kernel.** The raw assay carries a distance gradient on the pooled control features — **−0.0744 SD [−0.1306, −0.0182], p = 0.0145**. **N5 removes it** (+0.0038, p = 0.72). **N2 matched decoys do not** (−0.0642, p = 0.0124 — 86 % undiminished). `results/phase3/a7_summary.csv`. Two qualifications travel with this every time: the number is on the **pooled** control set (40 probes + 609 codewords + 21 genomic controls), never call it a negative-control-*probe* number, since the 40 named probes alone are flat (−0.0225, p = 0.129) and they are the pre-registered primary A7 response; and `neg_probe_rate` is flat naively (+0.0113, p = 0.232), so the gradient is per-cell **detection efficiency** projected onto distance-to-sender, which is precisely why a match on neighbourhood covariates cannot see it.
   - **The composition-matched protocol.** Matched decoys remove **1.6 %** of the naive amplitude (SF 0.9837 [0.973, 0.994]) where **the same variables as covariates on the same fits** remove **85.4 %** (SF 0.1461 [0.052, 0.246]) — a factor of fifty, with matching quality that passes the §8 Test 5 gate in 100 % of matches. `results/phase3/compmatch_reruns.csv`.
   - **Figure 1a/1b, which already said it.** Matched-decoy \|relative bias\| in $\hat\lambda$ beats naive in only 12 of 20 synthetic grid cells and is worse in 8; mean CI coverage where baseline autocorrelation length ≥ 2λ is 0.51 naive / 0.35 matched-decoy / **0.85 nuisance-conditioned**. `figures/figure1_data.csv`.

   **The consequence is pre-registered and must be in the text, not only in the repo: no naive or N2-only kernel may be reported on this platform** (`reports/PREREG_PHASE8.md` §10 rule 8), and wherever the matched SF appears, `type_adj` (65.9 %) and `typecomp_adj` (85.4 %) appear beside it. This overturns the project's own prior — see the correction in §22 Step 3 N2 — and an external novelty review rated it the strongest of the project's six findings, so it should be argued, not buried.
7. **"Torus shifts on irregular windows are a solved problem — Lotwick & Silverman 1982, Mrkvička et al. 2021 variance correction. Why did you not use the standard fix?"** *(Added 2026-08-27, citation audit Task 2.)* We cite both, we state that toroidal shifts are undefined on non-rectangular windows, and we report that Mrkvička et al. specifically predict that a union-of-aligned-rectangles extension — our tiling — increases the liberality of the test. **We ran it** *(updated 2026-08-27; this response previously read "we either run the variance-corrected shift as an N3 variant or justify tiling against it explicitly")*, implemented from the retrieved full text, and it is now the **primary** N3: **N3-var 0.996 [IQR 0.975, 1.007]** against N3-tile 0.971 and the published bounding-box 0.999, with N4-var 0.985 against N4-tile 0.924 and a published 0.947. The window-matched N3-var is 0.995, so the agreement is not an artefact of computing statistics on different amounts of data (`reports/CS_PHASE8_TORUS_VAR.md` §1–§2, `results/phase3/sf_summary_var.csv`). **And we then measured what the violation costs**, which is the part no spatial-omics paper has: on an irregular synthetic window under the null of independence, the **tiled** torus rejects at **0.080–0.118 against a nominal 0.05 — up to 2.4×** — while the variance correction (RS_count) holds **0.033–0.060** across every window and correlation scale (`results/phase3/var_sim_calibration.csv`). Reported against interest, that means **our own C1 correction replaced a liberal test with a more liberal one**, which is exactly why N3-var and not N3-tile is presented as primary. Two caveats travel with it: the calibration study is synthetic (Gaussian fields, 100 sampling points) and establishes the *direction* of the tiling effect, not a type-I error number for the Phase 3 fits — the instrument for that is A7; and on the real data the tiled null looks slightly *more* conservative (SF 0.971 vs 0.999), because at a 1,200 µm tile side and a pooled λ̂ of 15.7 µm the seams are ~76 λ̂ apart. **Presenting the N3 degeneracy as a discovery rather than as the first quantification on real tissue of a known classical pathology is the single easiest way to lose a statistician reviewer.**
8. **"You have no spatial-statistics citations. This is the spatial-confounding problem and it has a literature."** *(Added 2026-08-27.)* It is, and we now cite it: Hodges & Reich (2010) for the problem, Dupont et al. (2022) for the residualise-the-covariate remedy that N6 implements, and Khan & Calder (2022) / Zimmerman & Ver Hoef (2022) for why we do not simply orthogonalise. Three sentences in Methods.
9. **"Voyager already computes Moran's I on Xenium negative control probes and concludes there is no technical spatial trend. What is new?"** *(Added 2026-08-27.)* We cite Voyager and Ren et al. (2025) approvingly, report our **own** Moran's I on the controls alongside the kernel amplitude, and state the distinction: a near-zero *global* autocorrelation does not preclude a projection onto one *specific* covariate. Our construction is a **negative control outcome for the estimand itself** (Lipsitch et al. 2010) — refitting the model under test — not a generic spatial statistic. **The sentence "nobody in this literature reports it" must be struck wherever it appears.** **⚠ RESOLVED 2026-08-27 — Moran's I has been run; see reports/CS_PHASE8_MORAN.md §5 for drop-in text. Original marker: 2026-08-27: the Moran's I half of this response is a promise, not a result.** The string does not appear anywhere in `code/` or `results/` — nothing in the repo computes it. Until it does, this objection is answered by assertion. It is one cheap computation on features already in the M1 bundles, and the novelty review rates it the **highest-priority open gap** (`reports/NOVELTY_ASSESSMENT.md` §4 O1). What *is* measured and can be quoted today: the kernel amplitude on the pooled controls (−0.0744, p = 0.0145), the flat 40-probe primary response (−0.0225, p = 0.129), and the estimator's measured type-I error of **9–16 %** against a 5 % nominal (9–15 % on the four count-based responses; the 16 % upper end is `neg_probe_rate`, whose denominator is itself an N5 column and is not a clean null) — `results/phase3/a7_summary.csv`.

## 30. Paper Outline

Target 8 pages plus appendix; confirm against the template on Day 0. *(2026-08-27: 8 pages is the ICBINB-BIO full track. ml4spatialbio is 4 pages, at which only a subset of §5 fits — see the venue reassessment subsection in §29 before laying this out.)*

1. **Abstract.** Two versions written Day 3, chosen Day 10.
2. **Introduction.** Senescence spreads locally; how far is a central open question with therapeutic stakes. Existing evidence is descriptive. Contemporaneous benchmarks show spatial CCC inference is poorly calibrated. We ask whether the senescence kernel is estimable at all, and what it is once controlled. *(Added 2026-08-27: the introduction must also set up the two claims the paper is now actually strongest on — that a matched-decoy contrast is **not** the conservative option on imaging ST, and that random-shift nulls are being used on windows where the classical literature says they are undefined — and it must state, in the introduction rather than in Limitations, that no published length constant with uncertainty exists for this quantity, so the contribution is a **bound and an identifiability argument**, not the overturning of a number.)*
3. **Related Work.** Spatial senescence mapping (SenNet review, Ma et al., scDOT, Farzad et al., DeepScence). Spatial CCC methods (NCEM, COMMOT, SCILD, RGAST, SpaTalk, HARMONIC). Reliability and benchmarking (CellWHISPER, CONCISE, ctSVG benchmark). Mechanistic models of senescence spread (Martin et al., Boolean SASP models). **Spatial statistics** — random-shift nulls and their window requirements (Lotwick & Silverman 1982; Mrkvička et al. 2021), spatial confounding and the restricted-spatial-regression critique (Hodges & Reich 2010; Dupont et al. 2022; Khan & Calder 2022; Zimmerman & Ver Hoef 2022), and negative-control outcomes (Lipsitch et al. 2010), plus the existing negative-control-probe spatial diagnostics we are *not* the first to run (Voyager; Ren et al. 2025). [Added 2026-08-27, citation audit Task 2 — this paragraph previously had no spatial-statistics citations at all, and "Zhao et al." should have read "Ma et al."] **Three further additions, 2026-08-27 (plan update D12/D13), each closing a named objection:** senescence-caller disagreement is thoroughly established and must be cited as prior work rather than framed as a finding — add **ICE** (*Genome Biology* 2026, doi:10.1186/s13059-026-03997-0), **markeR** (*NAR Genomics and Bioinformatics* 2026, 8(2):lqag057) and bioRxiv 2026.01.02.697374, whose *"apparent concordance in prior studies may reflect circular validation"* is the closest published statement to our own §5.9; **stAge** (bioRxiv 2025.11.23.689860, spatial biological-age hotspots) is the nearest contemporaneous "aging has spatial gradients" claim and is absent from `references.bib`; and DeepScence's `CDKN1A` sign-anchoring step must be cited as **documented method behaviour**, so that §5.9's anchor instability reads as a user-facing comparability caveat rather than as a bug we found.
4. **Methods.** Kernel model and families; sender and response definitions with disjointness; null battery; synthetic generator; statistics.
5. **Results.** *(5.1–5.4 unchanged in scope; 5.5–5.9 added 2026-08-27, plan update D12/D13, because the outline predated half the work. Figures 5 and 6 are defined in `Phase7_Minimal_Human_Replication (1).md` §19.)*
   - 5.1 Identifiability regimes on synthetic tissue (Fig. 1). **Carry panels a and b explicitly, not only c:** matched-decoy relative bias in $\hat\lambda$ beats naive in only 12 of 20 grid cells and mean CI coverage where baseline autocorrelation length ≥ 2λ is 0.51 naive / 0.35 matched-decoy / 0.85 nuisance-conditioned. This is where 5.6 is first visible, on ground truth.
   - 5.2 The naive gradient and its behavior under nulls (Fig. 2a–d). Naive amplitude **0.329**; surviving fraction under N2+N5+N6 **0.088** [−0.017, 0.234].
   - 5.3 Controlled kernel estimates across tissues and programs (Fig. 3). Controlled amplitude **0.029** [−0.007, 0.084] against a detectable bound of **0.183** at 80 % power **over 0–100 µm** — state the distance range in the same sentence as the bound, and justify the 100 µm window (λ̂ ≈ 15.7 µm pooled, so ≈ 6λ; the endometrium nearest-neighbour calibration of 45–211 µm) or "we found nothing" is answerable with "you looked in the wrong place".
   - 5.4 Existing tools under the same nulls (Fig. 4). Keep CellWHISPER's ">90 %" as an **inferred** count ratio, and keep the null attribution straight: their null is our N1.
   - **5.5 The estimand's own estimator, refit on the assay's negative controls (Fig. 2h).** The naive kernel returns **−0.0744 SD** [−0.1306, −0.0182], p = 0.0145 on the pooled control features; the 40 named negative-control **probes** — the pre-registered primary response — are **flat** (−0.0225, p = 0.129); `neg_probe_rate` is flat naively (+0.0113, p = 0.232), so the gradient is per-cell detection efficiency projected onto distance-to-sender. Report a Moran's I of our own controls beside the amplitude (**not yet computed — see objection 9**) and the measured type-I error, **9–16 % against 5 % nominal** (9–15 % on the four count-based responses), as one Methods sentence and one table row, with its three caveats: which response, powered only pooled (a single fit resolves ±0.134 SD), and it is the *estimator's* rate, not the reportable-fit filter's (3.0–13.3 %). `results/phase3/a7_summary.csv`.
   - **5.6 A matched-decoy contrast does not remove what covariate adjustment removes.** The N2-vs-N5 inversion, from three directions: A7 (N5 +0.0038, p = 0.72; N2 −0.0642, p = 0.0124), the composition-matched protocol (1.6 % vs **85.4 %** from the same variables on the same fits — a factor of fifty), and Figure 1a/1b. Named mechanism in one sentence; consequence pre-registered (no naive or N2-only kernel). **An external novelty review rates this the project's strongest contribution — give it its own subsection and a table, not a footnote.**
   - **5.7 Random-shift nulls on a non-convex window.** Frame as **import and quantification, never as discovery**: the rectangular-window requirement is Lotwick & Silverman (1982) and the variance correction for general compact windows is Mrkvička et al. (2021). Ours are (i) the cost on real tissue — **35.5 %** of N3's shifted senders land outside the tissue under a bounding-box wrap and 19.9 % of N4's (`1 − frac_in_occupancy`; do **not** merge this with the 22.8 % / 8.0 % `1 − frac_retaining_a_neighbour` column), and at a ≤5 % out-of-tissue tolerance only 1–66 of 38,080–108,375 candidate offsets are admissible with one section admitting only the identity; (ii) a **direct calibration measurement**, which no spatial-omics paper has — the tiled torus rejects at **0.080–0.118 against a nominal 0.05, up to 2.4×**, while the variance correction holds **0.033–0.060**, so our own C1 correction replaced a liberal test with a more liberal one and **N3-var (0.996), not N3-tile (0.971), is the primary corrected null**; and (iii) the **FFT enumeration** — "which translations keep ≥ x % of a point set inside a mask?" is a circular cross-correlation, so one `rfft2` gives the exact admissible set over all offsets at once. **Lead the methods paragraph with the FFT trick; it converts a rediscovery into a tool.**
   - **5.8 Two-arm replication: mouse liver and human spleen (Figs. 5 and 6).** Fig. 5: surviving fraction by null, both arms side by side; the two-arm comparison table as a forest plot; the geometric predictions (Poisson identity r² 0.984, grid-railing rate) against their expected values in both arms. Fig. 6: DeepScence native (human) versus ortholog-remapped (mouse), caller agreement before and after conditioning on cell type and depth decile, and whether the anchor instability transfers. **Status: the mouse arm is complete and the human arm is acquired and structurally verified but runs in Phase 10, behind the `phase8-frozen` tag. If it has not run by submission, say so and report the mouse arm alone — do not imply a replication that has not happened.**
   - **5.9 What we withdrew, and why that is a result.** The caller-agreement claim: "statistically independent" is falsified (**1.212×**, p = 1.8 × 10⁻¹⁰⁶ at eleven sections; already **1.128×**, p = 4.4 × 10⁻⁸ on the published two sections once the sender set is fixed), and the cause is a contaminated sender/response split, not a power problem. Do **not** claim novelty for "senescence callers disagree" — DeepScence, SenCID, SenePy, ICE, markeR and Ntintas et al. all report it; claim only the coverage-and-definition sensitivity, and only with the decomposition attached.
6. **Corrections, pre-registration, and reproducibility.** *(New section, 2026-08-27. The project's single strongest structural defence against "garden of forking paths" is currently invisible in the outline.)* The pre-registration (`reports/PREREG_PHASE8.md`) and the `phase8-frozen` tag; the corrections ledger (`reports/CORRECTIONS.md`) giving, for every headline number, the pre value, the post value and the attributed cause; and the fact that **every correction moved against interest and none changed the conclusion**. Cite the pre-registration in the paper, not only in the repo.
7. **Discussion.** What the bound implies for senescence containment — it is a *constraint* on Martin et al.'s model, which is the most useful thing a negative result can be. What it implies for how spatial senescence claims should be reported: run the negative-control-outcome kernel, do not trust a matched-decoy contrast alone, and do not use a torus shift on a tissue section. **State plainly that there is no published length constant with uncertainty for senescence spatial influence, so this bound contradicts no prior number — a strength and a risk, and the reason Figures 1 and 4 are load-bearing rather than supporting.** Limitations: inferred and non-independent senescence calls, static snapshots, targeted panels, one mouse tissue plus one human tissue confounded with species, `denoise=False` as the frozen primary with `denoise=True` as the published-default sensitivity, and no perturbation. The decisive future experiment.
8. **Appendix.** Full sweeps, null distributions, coverage tables, gene sets and the intersection matrix (Tier A ∩ ∪B = 0 — but note the seven Tier B modules are **not** mutually disjoint, 18 of 21 pairs share genes, and no text may claim otherwise), sensitivity analyses, the eleven-variant shift-null family with its destructiveness diagnostics, the five-seed composition-matched protocol, the per-section caller tables, and the deviation tables.

## 31. Key References

Verify every bibliographic detail before submission. DOIs given where confirmed as of August 17, 2026; several years and page numbers should be double-checked against the publisher record.

**Senescence in space**
1. **Ma S**, Ji Z, Zhang B, … Zhang W, Gu Y, Liu G-H (2024). Spatial transcriptomic landscape unveils immunoglobin-associated senescence as a hallmark of aging. *Cell* 187(24):7025–7044.e34. doi:10.1016/j.cell.2024.10.019; PMID 39500323 — **closest prior work**. [Corrected 2026-08-27: listed here as "Zhao et al." Zhao L is 14th of 47 authors. Verified from PubMed and Crossref: nine tissues, male mice, Stereo-seq. Abstract states SSSs "serve as epicenters for heightened inflammation that compromises surrounding cells in a distance-dependent manner"; the finer per-pathway monotone-with-distance claim is still **not rendered from the publisher** — see references.bib AUDIT (3.1).]
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
17. Acosta JC, Banito A, Wuestefeld T, *et al.* (2013). **A complex secretory program orchestrated by the inflammasome controls paracrine senescence.** *Nature Cell Biology* 15(8):978–990. doi:10.1038/ncb2784 — [corrected: "Paracrine senescence transmission" is not this paper's title. **See also the 2026 Author Correction**, *Nat Cell Biol* 28(6):1343, doi:10.1038/s41556-026-01959-z, whose content is still unread.]
18. **Russo T, Riessland M** (2026). Characterizing the SASP-dependent paracrine spreading of senescence between human brain cell types. ***Aging Cell*** **25(8):e70673**. doi:10.1111/acel.70673; PMID 42601837 — CCL2, CXCL12/DPP4, MIF, ACKR3/CXCR7. [Corrected 2026-08-27: cite the peer-reviewed version, not the February 2026 bioRxiv preprint listed here. First author is **Taylor** Russo, verified from PubMed.]
19. A model of the onset of the SASP after DNA-damage-induced senescence (2017). *PLOS Computational Biology*. doi:10.1371/journal.pcbi.1005741
20. Dissecting the heterogeneity of senescence: primary and secondary senescent states. PMC11689308

**Spatial CCC methods**
21. Fischer et al. NCEM: Modeling intercellular communication in tissues using spatial graphs of cells. *Nature Biotechnology*. doi:10.1038/s41587-022-01467-z
22. Cang et al. COMMOT: Screening cell-cell communication in spatial transcriptomics via collective optimal transport. *Nature Methods*
23. SCILD (2026). Advancing spatial cellular communication inference with ligand diffusion and transport model. *Communications Biology*. doi:10.1038/s42003-025-09413-w
24. SpaTalk. Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data. *Nature Communications*
25. Gong Y, Yuan X, Yu Z (2026). **Empowering multifaceted analysis of spatial transcriptomics data with RGAST.** *Briefings in Bioinformatics* 27(3):bbag298. doi:10.1093/bib/bbag298 — [corrected 2026-08-27: the title/year pairing listed here never existed. The published framing is **representation learning benchmarked on spatial domain identification**, with CCC as one downstream task; **do not cite it as a CCC method paper.**]
26. HARMONIC: Histology-aware graph for modeling intercellular communication in spatial transcriptomics. bioRxiv, January 2026

**Reliability, calibration, identifiability — the framing citations**
27. Kumar A, Rivera F, Aggarwal B, Zhang N, Coskun AF, Sinha S. CellWHISPER: Inference of direct cell-cell communication from spatial transcriptomics. bioRxiv 2026.01.07.697982; v2 retitled *CellWHISPER disentangles direct cell–cell communication from structural proximity*. doi:10.64898/2026.01.07.697982 — CellChat v2, COMMOT and SpaTalk "predict similar numbers of interactions on real … and randomized … data, indicating poor specificity and false positive rates (FPR) >90%". **The randomization is a within-cell-type location permutation (= our N1). There is no torus shift in this paper.** The >90% is inferred from an interaction-count ratio, not measured against a nominal level.
28. CONCISE: Spatial co-expression and cell-cell communication inference from spatially resolved transcriptomics. bioRxiv, June 2026
29. Benchmarking cell-type-specific spatially variable gene detection methods (2026). *Briefings in Bioinformatics* 27(2):bbag190
30. **Gu R, Zhang RZ, Miles CE** (2026). Identifiability limits of physics-informed inference for spatial stochastic dynamics from static snapshots. arXiv:2607.01749 — [author forenames corrected 2026-08-27 from the arXiv record: **Rujie** Gu, **Ray Zirui** Zhang.]

**Spatial statistics — random-shift nulls and spatial confounding** *(added 2026-08-27, citation audit Task 2. Before this pass §31 contained no spatial-statistics methods paper at all. Full verification notes and the "what each supports" line live in `references.bib`.)*

31. **Lotwick HW, Silverman BW (1982).** Methods for analysing spatial processes of several types of points. *JRSS-B* 44(3):406–413. doi:10.1111/j.2517-6161.1982.tb01221.x — **the origin of the toroidal/random-shift null. This is what N3 cites, not CellWHISPER.**
32. **Mrkvička T, Dvořák J, González JA, Mateu J (2021).** Revisiting the random shift approach for testing in spatial statistics. *Spatial Statistics* 42:100430. doi:10.1016/j.spasta.2020.100430; arXiv:1911.00240 — torus correction is **liberal**; it "assumes a rectangular observation window"; extending it to unions of aligned rectangles (our tiling) "would increase the amount of cracks … [and] increase the liberality"; the **variance correction** is the remedy for general compact windows. **Directly answers reviewer objection O3, and directly challenges our tiling choice.**
33. **Hodges JS, Reich BJ (2010).** Adding spatially-correlated errors can mess up the fixed effect you love. *The American Statistician* 64(4):325–334. doi:10.1198/tast.2010.10052 — names the **spatial confounding** problem this project is an instance of.
34. **Dupont E, Wood SN, Augustin NH (2022).** Spatial+: a novel approach to spatial confounding. *Biometrics* 78(4):1279–1290. doi:10.1111/biom.13656 — residualising the covariate against the spatial field, which is what N6 does.
35. **Khan K, Calder CA (2022).** Restricted spatial regression methods: implications for inference. *JASA* 117(537):482–494. doi:10.1080/01621459.2020.1788949 — and **Zimmerman DL, Ver Hoef JM (2022)**, *The American Statistician* 76(2):159–167, doi:10.1080/00031305.2021.1946149 — why we do **not** simply orthogonalise against the spatial field. One Limitations sentence.
36. **Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010).** Negative controls: a tool for detecting confounding and bias in observational studies. *Epidemiology* 21(3):383–388. doi:10.1097/EDE.0b013e3181d61eeb — **the correct framing for A7**: a negative control *outcome* for the estimand itself.
37. **`spatstat.random::rshift.ppp` reference manual** — "The window must be a rectangle. Toroidal shifts are undefined if the window is non-rectangular." The one-line software statement of the N3-occ problem.
38. **Moses L, Einarsson PH, … Pachter L (2023).** Voyager: exploratory single-cell genomics data analysis with geospatial statistics. bioRxiv 2023.07.20.549945 (**still a preprint**) — and **Ren P *et al.* (2025)**, *Nat Commun* 16, doi:10.1038/s41467-025-64292-3. **Both already compute Moran's I on negative control probes.** Cite them approvingly; the claim that nobody reports a negative-control-probe spatial diagnostic is false.

## 32. Reading List

Six papers, in this order, roughly one working day. Do not skip this to start coding sooner; three of them will change your design.

1. **Ma et al. (2024), *Cell* 187(24):7025–7044.e34.** The closest prior result. Understand exactly what they claimed and how they measured it, because you are re-examining it. [Corrected 2026-08-27 from "Zhao et al."]
2. **"Cellular senescence in human liver under normal aging and cancer" (2026), *Cell Genomics*.** Your Rank 1 dataset. Read the spatial sections carefully, especially the clustered-versus-isolated senescent hepatocyte analysis and the microenvironment comparison that fell short of significance.
3. **CellWHISPER (bioRxiv, Jan 2026).** The >90% FPR benchmark and your framing citation. **Its null is a within-cell-type location permutation — this project's N1. It is not the source of the torus-shift null and contains no torus shift at all** [corrected 2026-08-27, citation audit CIT-1]. For N3, read **Lotwick & Silverman (1982)**, JRSS-B 44(3):406–413, and **Mrkvička et al. (2021)**, *Spatial Statistics* 42:100430 — the second is about precisely our failure mode and should be item 3b on this list.
4. **CONCISE (bioRxiv, June 2026).** Spatial autocorrelation inflating type I error. Second framing citation.
5. **DeepScence (*Cell Genomics*, Dec 2025).** Your sender-calling tool, the CoreScence gene set, and a useful survey of the nine published senescence gene sets and their disagreement.
6. **Martin et al. (2023), *Aging Cell*.** The containment paradox your measurement speaks to.

Optional: the 2023 *Nature Aging* SenNet review for the field's own framing, and the 2026 *FEBS Open Bio* signature comparison for the gene-set heterogeneity numbers you will cite in Limitations.
