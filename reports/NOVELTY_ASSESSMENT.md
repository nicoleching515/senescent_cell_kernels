# Novelty, significance and relevance assessment

> ## ⚠ CORRECTED 2026-08-27 — one recommendation in this report is FALSIFIED
>
> This report advises the project to argue that **Moran's I and the A7
> distance-to-sender kernel "ask different questions" and disagree** (§2.1 point 3
> and §4 O1). **Moran's I has since been run** (`reports/CS_PHASE8_MORAN.md`,
> `results/moran/`) and **the disagreement claim is falsified**: across the 12
> control and module fields, |Moran's I| and |A7 naive amplitude| rank together at
> **Spearman rho = +0.895 raw, +0.944 cell-type-centred**. Writing "the two tests
> disagree" would be visible to any reviewer in a single plot.
>
> **The replacement is stronger, because it is measured rather than asserted.**
> Using the real lambda-hat: the entire A7 gradient contributes **Delta-I = 2.2e-4,
> i.e. 0.83% of the observed control Moran's I**, and the smallest kernel amplitude
> Moran's I could resolve is **0.362 SD** — five times the A7 control gradient
> (0.074 SD) and **larger than the project's own naive biological amplitude
> (0.291 SD)**. I verified that bound from `results/moran/moran_kernel_power.csv`.
> **Moran's I could not have detected the paper's headline effect either.**
>
> So "different question" survives — but justified by **power, not orthogonality**.
> Use that framing. Everything else in this report stands.
>
> **A second, genuinely new point came out of the same run:** Voyager's
> controls-vs-genes contrast is an **abundance** contrast. Genes matched to the
> controls on total counts give identical Moran's I (-0.00018 vs -0.00012); a
> negative control probe carries ~21 counts per section against a median gene's
> 5,885. The per-feature statistic has no power at control abundance.


**Prepared:** 2026-08-27. **Scope:** the five findings named in the assessment brief, judged
against the retrieved literature, not recall.
**Author:** research-analyst pass. Read-only; no repo file other than this one was modified.

---

## 0. How to read this

Every novelty judgement below is tagged one of three ways:

| Tag | Meaning |
|---|---|
| **PRIOR FOUND** | I retrieved a specific paper/resource that does the thing, or does enough of it to threaten the claim |
| **SEARCHED, NONE FOUND** | I ran named searches and got nothing; the search strings are listed so the gap is auditable |
| **COULD NOT VERIFY** | The source exists but I could not retrieve it (paywall / rate limit); stated as unresolved, not guessed |

**What I could not retrieve** (bioRxiv returned HTTP 429 repeatedly, Springer/Nature/Cell
returned 303/403): full text of CONCISE (bioRxiv 2026.06.22.733860), the *Genome Biology* 2026
CCI benchmark (10.1186/s13059-026-04063-5), SpatialArtifacts (bioRxiv 2026.05.15.725260),
SPLISOSM full text, and the 10x "Metrics Matter" blog. For each I have a secondary description
from search-result summaries and I flag it. CellWHISPER, the Voyager Xenium vignette, the
*Nature Communications* 2025 platform benchmark, the eLife off-target paper, and the *Cell
Genomics* liver paper I did retrieve in full or near-full.

---

## 1. LEAD WITH THESE — prior work that undercuts a current claim

Ordered by how much damage each does if a reviewer finds it and you have not.

### U1. The negative-control-probe spatial diagnostic already exists, in the field's most-used ESDA tutorial

**This directly falsifies the sentence in `reports/SUBMISSION_PATCH_2026-08-29.md` §6:
"nobody in this literature reports it."**

The Voyager Xenium vignette (Pachter lab) computes **Moran's I on the negative control probes
and codewords** and reads it exactly as this project proposes to:

> "generally the negative controls are tightly clustered around 0, while the real genes have
> positive Moran's I, **which means there is generally no technical artifact spatial trend**."

It also already reports the *confound mechanism* this project rediscovers:

> "Generally there are more negative control spots in the region with higher cell and transcript
> density … There don't visually seem to be regions with more negative control spots not
> accounted for by cell and transcript density."

and the 2022 Xenium preview artifact:

> "the first Xenium preview data from 2022 has a region in the top left corner with more negative
> control probes detected that seems like an artifact."

Source: <https://pachterlab.github.io/voyager/articles/vig5_xenium.html>; method paper Moses L,
Einarsson PH, Jackson K, Luebbert L, Booeshaghi AS, Antonsson S, Bray N, Melsted P, Pachter L,
*Voyager: exploratory single-cell genomics data analysis with geospatial statistics*, bioRxiv
2023.07.20.549945 (PMC10461913) — **still a preprint** as of this search, which matters for how
you cite it.

**Second instance, peer-reviewed:** Ren P, Zhang R, Wang Y, *et al.* (2025). *Systematic
benchmarking of high-throughput subcellular spatial transcriptomics platforms across human
tumors.* **Nature Communications 16**, doi:10.1038/s41467-025-64292-3:

> "Spatial autocorrelation analysis using Moran's I revealed stronger aggregation of negative
> control signals in CosMx 6K …, indicating higher background interference."

**Consequence.** "Nobody reports the negative-control-probe spatial diagnostic" is false and must
be struck. What survives is narrower and still real — see §2.1.

### U2. The torus-shift-on-irregular-windows problem is 40+ years old in spatial statistics, and the fix you did not use is 5 years old

Finding 2 is a **rediscovery**, and the classical literature has a *better* answer than tiling.

- **Lotwick HW, Silverman BW (1982),** *Methods for analysing spatial processes of several types
  of points*, JRSS-B — the origin of the toroidal shift test. The window must be a rectangle;
  toroidal shifts are **undefined** for non-rectangular windows. This is standard enough to be in
  the `spatstat` documentation for `rshift.ppp` (`edge="torus"` requires a rectangular window).
- **Mrkvička T, Dvořák J, González JA, Mateu J (2021),** *Revisiting the random shift approach
  for testing in spatial statistics*, **Spatial Statistics 42**:100430; arXiv:1911.00240. This
  paper is about precisely your failure mode. Its findings, in its own terms: torus correction
  glues distant parts of the data together, introduces artificial cracks in the correlation
  structure and makes the test **liberal**; the proposed **variance correction** instead shifts in
  Euclidean geometry and *discards the part shifted outside W*, "**can also be used for irregular
  windows**", and — the sentence that matters most for you — "different amounts of data are
  dropped for different shift vectors and the variance of the test statistic varies greatly, so
  the variance needs to be standardized before performing the test."
- **Dvořák J, Mrkvička T, *et al.* (2022)** show the variance-corrected shift test holds the
  nominal level even under strong autocorrelation.

Your N3-occ is a hard **retention constraint** (≤5% out of tissue) where the literature uses a
**soft drop-and-standardize**. The degeneracy you found (admissible set collapses to near-identity;
one section admits only the identity) is the predictable consequence of imposing a hard constraint
where the standard method imposes a soft one. A reviewer from spatial statistics will see this
immediately.

**Consequence.** Do not present §3 of `CS_PHASE7_C1.md` as a discovery. Present it as *the first
demonstration that a known classical pathology is binding on real tissue sections at the
magnitudes spatial-omics papers actually use* — and either implement variance correction or state
why tiling is preferred over it. See §2.2.

### U3. `CellWHISPER` is cited in the plan as "the source of the torus-shift null". It is not.

`SASP_Kernel_Master_Plan.md` §32 item 3 says CellWHISPER is "the source of the torus-shift null."
I retrieved the preprint. Its null is a **within-cell-type location permutation**:

> "To assess statistical significance, CellWHISPER constructs a null model by permuting cell
> identities within each cell type," which "preserves cell-type–specific gene expression
> distributions and spatial organization while disrupting spatial proximity between signaling
> cells."

**No torus shift or toroidal translation appears anywhere in the paper, and it does not use
negative control probes or blank barcodes.** Citation: Kumar A, Rivera F, Aggarwal B, Zhang N,
Coskun A, Sinha S (2026), *CellWHISPER: Inference of Direct Cell-Cell Communication from Spatial
Transcriptomics*, bioRxiv 2026.01.07.697982 (also SSRN 6815051, v2 retitled *CellWHISPER
disentangles direct cell-cell communication from structural proximity*).

Their null is **your N1**, not your N3. Fix the citation before submission — this is exactly the
kind of error a reviewer who is a CellWHISPER author will catch, and it is load-bearing because
N3 is the null your Figure 4 is built on.

Second nuance on the same citation: the ">90% FPR" is **not a measured false-positive rate**. It
is the observation that competing tools "reported similar numbers of CCC mechanisms in the
randomized and original datasets, suggesting a high false positive rate." Quoting it as a measured
FPR overstates it. Say "reported comparable interaction counts on randomized coordinates,
implying FPR >90%" — which is what the plan's §1 actually does; keep that wording and do not let
it drift.

### U4. Senescence-caller disagreement is thoroughly established, including the circularity point

Finding 4's *headline* (callers disagree) is not novel by a wide margin:

- Qu *et al.* (2025), **DeepScence**, *Cell Genomics* 5(12):101035 — itself built on a systematic
  comparison of **nine published senescence gene sets** and their disagreement; CoreScence is the
  39 genes reported by ≥5 of them.
- **SenCID**, *Cell Metabolism* 2024 — "none of the feature genes selected by [six] senescence
  identification models were shared by all six models, even for … CDKN2A (p16) and CDKN1A (p21)."
- **SenePy**, *Nature Communications* 2025 (doi:10.1038/s41467-025-57047-7) — the cell-type-specific
  landscape argument, i.e. that a single signature cannot transfer.
- **ICE**, *Genome Biology* 2026 (doi:10.1186/s13059-026-03997-0) — reports partial overlap
  (88% of ICE cells also called by GSEA/SenSig, each with non-overlapping remainder).
- **markeR**, *NAR Genomics and Bioinformatics* 2026, 8(2):lqag057 — a toolkit whose stated purpose
  includes surfacing method disagreement.
- **Ntintas *et al.* (2026)**, *FEBS Open Bio*, doi:10.1002/2211-5463.70134 — the signature
  pros-and-cons overview already in `references.bib`.
- *Integrative transcriptomic identification of cellular senescence beyond marker limitations*,
  bioRxiv 2026.01.02.697374 — goes further than anything in the repo: "curated gene sets show
  **opposing enrichment patterns** in experimentally defined senescent cells, suggesting apparent
  concordance in prior studies may reflect **circular validation**." This is the closest published
  statement to your own §0.3 circularity thesis and it should be cited; not citing it looks like a
  gap in the related work.

### U5. The DeepScence polarity flip is a documented step of the method, not a discovery

DeepScence's own procedure: the bottleneck neuron with the highest mean absolute correlation to
CoreScence genes is selected, then **"DeepScence calculates the Pearson correlation between the
output senescence score and the expression of the CDKN1A gene … If the correlation is negative,
the sign of the output score is flipped."** So the sign is set by a per-run data-dependent decision,
and an inversion between two runs is the *documented behaviour of the anchoring rule under a weak
CDKN1A signal*, not a surprise. Frame accordingly (§2.4).

### U6. "Zhao et al." is the wrong first author, and the repo already knows this

The brief still says "Zhao et al. (Cell, 2024)". The paper is **Ma S, Ji Z, Zhang B, … Liu G-H
(2024), Cell 187(24):7025–7044, doi:10.1016/j.cell.2024.10.019** — Zhao L is 14th of 47. The repo
corrected this on 2026-08-21 (D7 §B1/B2/B3) and `references.bib` keys it `ma2024spatial`. Verified
independently here. Make sure no draft text still says Zhao.

### U7. "Spatial confounding" is a named, mature statistics literature the project does not cite

The project's core methodological framing — a covariate of interest that is smoothly spatially
varying and collinear with the spatial random field, so that its coefficient is not identified —
is the textbook **spatial confounding** problem: Hodges JS & Reich BJ (2010) and restricted
spatial regression; Dupont E, Wood SN, Augustin NH, *Spatial+: a novel approach to spatial
confounding*, **Biometrics** 78(4):1279 (2022); and the subsequent critical literature (Hanks
*et al.* 2015; Khan & Calder 2022; Zimmerman & Ver Hoef 2022) showing RSR uncertainty
quantification is severely anticonservative. `references.bib` has **none of these** (30 entries;
none is a spatial-statistics methods paper). At an ML venue you may get away with it; at any venue
with a statistician on the PC you will not. Two or three of these citations in the Methods buy a
lot of credibility and cost nothing.

---

## 2. Per-finding assessment

### 2.1 Finding 1 — the negative-control-probe kernel test

**Verdict: PARTIALLY NOVEL. The headline claim as currently written is false; a narrower claim is
genuinely new and is the strongest thing in the project.**

**What is not novel** (see U1 in full):
- Negative control probes/codewords as a technical-noise QC feature: universal. 10x documents
  per-cell NCP metrics; Xenium Prime ships 40 negative control probes, 609 codewords, 21 genomic
  controls; SpaceTrooper (bioRxiv 2025.12.24.696336) computes NCP proportions per cell across
  Xenium/CosMx/MERFISH; the OSTA Bioconductor book has a QC chapter on them.
- **Spatial** analysis of negative controls: **done**. Voyager computes Moran's I on them and reads
  a near-zero value as "no technical artifact spatial trend"; Ren *et al.* 2025 (*Nat Commun*) uses
  Moran's I on negative-control signal to compare platform backgrounds.
- The *mechanism* you find (control counts track cell/transcript density): stated qualitatively in
  the Voyager vignette.
- Negative controls as an FDR device: 10x's own "Metrics Matter" position, and the general
  statistical machinery — Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010), *Negative controls: a
  tool for detecting confounding and bias in observational studies*, **Epidemiology** 21(3):383;
  and control-gene calibration in genomics (RUV; Gagnon-Bartsch & Speed; Wang, Zhao, Hastie & Owen,
  *Confounder adjustment in multiple hypothesis testing*, arXiv:1508.04178 / the RUV unification
  literature, which also warns that control-based calibration is **anti-conservative when there
  are few control genes** — you have 40).

**What I searched for and did not find** (searches: `"negative control probes" fitted against
distance spatial transcriptomics`; `negative control probes as null outcome spatial neighborhood
analysis distance to cell type`; `negative control probes used as null hypothesis test spatially
variable genes false positive rate`; `torus shift null spatial omics negative control`):

1. **No paper fits the estimand's own estimator to the negative control features.** Everyone
   computes a *generic* spatial statistic (Moran's I, per-cell rate) on negative controls. Nobody
   runs the actual model under test — here, a distance-to-nearest-sender kernel regression — with
   negative controls as the response. That is a **negative control outcome for the specific
   estimand**, which is Lipsitch's construction, and I found no instance of it in spatial omics.
2. **No paper reports which null removes the control gradient.** Your N2-vs-N5 result
   (`N2` leaves −0.061 SD, p=0.020, ~80% undiminished; `+N6+N5` gives +0.007 [−0.011, +0.025],
   p=0.41) is the most transferable thing in the whole repo and I found nothing like it anywhere.
3. Moran's I ≈ 0 (Voyager's test) and a distance-kernel amplitude of −0.070 SD are **not the same
   test**. A signal with no global autocorrelation can still project onto a specific covariate.
   That is your defence against U1 and you should state it explicitly with a Moran's I of your own
   controls alongside, so the reader can see the two tests disagree.

**Framing that survives review.** Kill "nobody reports it". Write instead:

> Negative control probes are routinely used to quantify background rate, and their global spatial
> autocorrelation has been used as a platform QC metric [Voyager; Ren 2025]. We instead use them
> as a **negative-control outcome for the estimand itself** [Lipsitch 2010]: we refit the identical
> distance-to-sender kernel with control counts as the response. The naive estimator returns
> −0.070 SD (p=0.023), a quarter of the naive biological amplitude in the same fits, so a naive
> Xenium distance kernel reports this in part. The full nuisance design removes it
> (+0.007, p=0.41); **a matched-decoy contrast does not** (−0.061, p=0.020).

**The one thing you must say in the same breath, or a reviewer will say it for you.**
`neg_probe_rate` (probe counts ÷ transcript counts) is **flat naively** (+0.014, p=0.079). So the
gradient is a **per-cell detection-efficiency/size effect projected onto distance-to-sender**, not
a spatial gradient in probe binding. Reporting "the raw assay is not flat" without that sentence
invites the reading that you have rediscovered "bigger cells have more counts". Reporting *with*
it is stronger, because it names the confounder and explains exactly why matched-decoy matching
on neighbourhood covariates cannot catch it.

**Significance and audience.** High and broad. This is a two-line diagnostic that any Xenium/CosMx/
MERFISH paper fitting distance-dependent anything can run for free, and the N2-insufficiency result
contradicts the standard prior that propensity-matched decoys are the conservative option. It also
directly contradicts your own `SASP_Kernel_Master_Plan.md` §23, which calls the matched-decoy
number "the single most important number in the paper" — that sentence has to change too.

### 2.2 Finding 2 — torus-shift nulls are degenerate on non-convex tissue

**Verdict: NOT NOVEL IN STATISTICS. NOVEL IN SPATIAL OMICS AS A DEMONSTRATION. Frame as import,
not discovery. The FFT admissible-set computation is a small genuine methods contribution.**

**PRIOR FOUND** (full detail in U2): Lotwick & Silverman (1982) established the toroidal shift and
its rectangular-window requirement; `spatstat`'s `rshift` documents that torus shifts are undefined
on non-rectangular windows; **Mrkvička *et al.* (2021), Spatial Statistics 42:100430** diagnoses
the liberality of torus correction, proposes the **variance correction** explicitly *because* it
"can also be used for irregular windows", and notes the exact trade-off you hit — different shifts
drop different amounts of data, so the test statistic's variance must be standardized. Dvořák
*et al.* (2022) show it attains nominal level under strong autocorrelation.

**SEARCHED, NONE FOUND**: no spatial-omics paper warning about torus/random-shift nulls on
non-convex tissue. Searches: `torus shift null spatial omics tissue shape irregular cells shifted
outside tissue`; `spatial transcriptomics null model permuting cell coordinates tissue boundary
irregular shape rotation null`. The nearest thing in spatial omics is the *rotation-invariance*
literature (PMC12318016; STANCE, *Nat Commun* 2025 doi:10.1038/s41467-025-57117-w), which is about
orientation-dependence of SVG tests, not about shift nulls leaving the tissue. The Voyager vignette
notes in passing that "for some types of spatial neighborhood graphs, such as the k nearest
neighbor graph, these cells [outside tissue] will also affect spatial analysis" and filters them —
awareness of the mechanism, not of the null-model degeneracy.

**Two things here are genuinely yours and should carry the section:**
1. The **quantification on real tissue**: 23% of shifted cells in the void under a bounding-box
   wrap; at ≤5% out-of-tissue tolerance only 1–63 of 38,080–108,375 candidate offsets are
   admissible (0.001–0.17%), median displacement **27 µm** against a median λ̂ of **12.8 µm**, and
   one of six sections admits only the identity (SF = −0.000 by construction). Nobody has published
   these numbers for spatial omics.
2. The **FFT trick**: "which translations keep ≥x% of a point set inside a mask?" is a circular
   cross-correlation, so one `rfft2` gives the exact admissible set over all offsets at once.
   Rejection sampling would never have revealed that the admissible set is one offset wide. This is
   a small, clean, reusable contribution and it is the reason the degeneracy was detectable at all.
   **Lead the methods paragraph with it** — it converts a rediscovery into a tool.

**The framing that will get you past a statistician reviewer:**

> Toroidal shift nulls require a rectangular window (Lotwick & Silverman 1982) and are known to be
> liberal because of the seam they introduce; the accepted remedy for irregular windows is the
> variance correction of Mrkvička et al. (2021), which drops the shifted-out portion and
> standardizes the variance. Tissue sections are non-convex, and spatial-omics practice has adopted
> shift nulls without this correction. We quantify the cost: [numbers], and give an exact
> FFT-based enumeration of the admissible offset set. We adopt tiling because [reason];
> we report the variance-corrected shift as a sensitivity analysis.

**Risk if you do not do this.** As written, `CS_PHASE7_C1.md` §6.3 says "this is a finding in its
own right, and it is transferable." A reviewer with a spatial-statistics background reads that as
not knowing the field. The same paragraph with the two 1982/2021 citations attached reads as
competence. The evidence cost is two citations.

**Recommendation with teeth:** implement Mrkvička variance correction as one more N3 variant
before submission if there is time. It is the null the classical literature says to use, it has no
retention/displacement trade-off, and if it returns SF ≈ 0.95–1.00 like N3-tile and N3-occ15 do,
your Contribution-3 claim becomes essentially unassailable. If it does not, you need to know now.

### 2.3 Finding 3 — a measured false-positive rate of 9–16% against 5% nominal

**Verdict: NOVEL IN CONSTRUCTION, MODEST IN WEIGHT. Do not report it standalone; report it as the
quantitative output of Finding 1.**

**Answering the brief's actual question — how often does the spatial CCC literature report a
measured FPR or calibration check?** Honest answer: **rarely, but no longer never, and the
"never" window closed about eight months ago.**

- **CellWHISPER** (bioRxiv Jan 2026) reports one, in the "interaction counts on randomized vs. real
  coordinates" sense — >90% for CellChat v2 / COMMOT / SpaTalk, <5% for itself. Retrieved and
  verified. Note it is a count-ratio argument, not a nominal-vs-empirical type-I-error rate.
- **CONCISE** (bioRxiv Jun 2026) reports type I error inflation for all competitors in simulation
  with spatially autocorrelated genes, including SpatialDM. **COULD NOT VERIFY** the full text
  (429s); I have the abstract-level description from two independent search summaries and the
  SpatialDM comparison quote.
- **SpatialDM** (*Nat Commun* 2023, doi:10.1038/s41467-023-39608-w) derives an analytic null and is
  therefore calibration-aware by construction, though its null assumes expression is independent of
  location — the assumption CONCISE attacks.
- **The March/April 2026 benchmarks** — *Benchmarking tools for deciphering cellular crosstalk in
  spatially-resolved transcriptomics*, **Genome Biology 2026**, doi:10.1186/s13059-026-04063-5
  (nine CCI methods, simulations + nine real ST datasets across Visium/Stereo-seq/Xenium), and the
  *Briefings in Bioinformatics* ctSVG benchmark 27(2):bbag190. **COULD NOT VERIFY** the Genome
  Biology full text (Springer 303). Search summaries describe F1/accuracy/pathway-relevance metrics
  and "realistic simulation settings"; I could not confirm they report a calibrated type-I error.
- **NCEM** (*Nat Biotechnol* 2022, doi:10.1038/s41587-022-01467-z) and **COMMOT** (*Nat Methods*
  2023) — I found no calibration or FPR analysis in either; searched
  `NCEM length scale validation criticism`, `COMMOT false positive calibration`.

**So the honest sentence is:** "Calibration is now being reported for the first time by
2026-vintage methods papers, each for its own estimator; it is essentially absent from the
established tools and from applied papers using them."

**What is actually new in yours.** Every FPR above is obtained by **randomizing real data**. Yours
is obtained from **assay-internal features with a known-zero effect, with no randomization at
all** — 9.1 / 10.3 / 10.9 / 12.7 / 16.4% CI-exclusion across five control families under the full
N6+N5 design. That is a different and in some ways stronger instrument: it measures the estimator's
type I error *in situ*, against the real spatial structure, including whatever the randomization
would have destroyed. I found no prior instance of this in spatial omics.

**Weaknesses a reviewer will press, and they are real:**
- It is **self-calibration** — an FPR for your own estimator, on your own data, which is the same
  move CellWHISPER makes for itself. It does not generalize to COMMOT/CellChat unless you also run
  those on the control probes, which you could and probably should (it is nearly free and it makes
  Figure 4 far stronger).
- 40 negative control probes is few, and control-gene calibration is **known to be anti-conservative
  with few controls** (the RUV literature). You use 5 families and 825 fits, which mitigates it,
  but say so.
- Your own §4.2 concedes a single A7 fit resolves only ±0.137 SD, 1.8× the conditioned biological
  amplitude. The test is only powered **pooled across sections**. That caveat has to travel with
  the number every time it is quoted, or the 9–16% figure is not defensible.
- The two readings you offer ("bootstrap FPR at 2–3× nominal" vs "residual confounding N5 misses")
  are not distinguished by the data. Say which you believe and why, or a reviewer will read the
  ambiguity as the number being uninterpretable.

**Framing.** One sentence in Methods and one row in a table, attached to Finding 1. Not an abstract
claim. "Our estimator's type I error, measured on assay-internal null features under the reported
design, is 9–16% against a 5% nominal — so we treat nominal CIs as optimistic by 2–3× throughout,
and the detectable bound in §X is stated accordingly."

### 2.4 Finding 4 — caller disagreement is coverage-sensitive

**Verdict: THE HEADLINE IS OLD (U4). THE COVERAGE-SENSITIVITY OBSERVATION IS, AS A STATISTICAL
PHENOMENON, ALSO OLD. THE SPECIFIC MEASUREMENT IS NEW BUT WEAK, AND IT IS CURRENTLY CONFOUNDED.
Report it as a self-correction and a methods caveat, not as a contribution.**

**Prior work on the headline:** U4, seven citations. Do not claim novelty for "senescence callers
disagree."

**Prior work on the coverage sensitivity:** I found no paper reporting section-count sensitivity of
senescence-caller agreement (searched `senescence caller agreement across sections coverage`;
`cell type annotation tool concordance depends on number of samples`). But the *statistical*
content — that a chance-corrected agreement ratio estimated on two units is unstable, prevalence-
dependent, and can invert — is classical: **Feinstein AR & Cicchetti DV (1990), *High agreement but
low kappa: I. The problems of two paradoxes*, J Clin Epidemiol 43(6):543–549**, and Cicchetti &
Feinstein's companion; kappa's marginal/prevalence dependence is standard textbook material.
A 1.693 → 0.737 inversion between an n=2 and an n=11 estimate is what small-sample instability
looks like. It is not a new phenomenon; it is your own earlier analysis being underpowered.

**The confound you must resolve before publishing the number.** `reports/SUBMISSION_PATCH_2026-08-29.md`
carries its own warning banner: the pooled ratio moved from 1.129× to **1.212×**, the band from
0.700–1.711 to **0.751–2.198**, and Tier A vs SenePy from "0.914×, below chance 11/11" to
"**0.972×, p=0.104 — NOT significantly below chance**", **because the Tier A gene set changed from
25 to 33 genes**, not because coverage changed. So two variables moved at once. As of that file the
M1 re-run was still in flight. **A reviewer will ask whether the effect is coverage or gene-set
definition, and right now the repo cannot answer.** The 6-section in-band subset (1.115,
p=2.6×10⁻¹⁹) partially addresses "is it the excluded sections", not "is it the gene set".

**Required evidence:** recompute the 2-section and 11-section agreement **with the frozen strict-33
Tier A on both**, so coverage is the only thing that varies. That is the entire claim. If it has
already landed from the M1 re-run, use those numbers; if not, the claim is not yet supportable.

**On the DeepScence polarity flip.** See U5 — the CDKN1A sign-anchoring is a documented step of the
method, so an inversion is the anchoring rule behaving as specified under a weak or inverted CDKN1A
correlation. Two consequences:
- Do **not** frame it as a defect in DeepScence. Frame it as: *DeepScence scores are anchored per
  run by the sign of their correlation with CDKN1A; in one of our two surgical arms that
  correlation inverted, so scores are not comparable across sections without re-anchoring. Users
  applying DeepScence across sections/batches should check the anchor.* That is a legitimate,
  citable practical caveat and I found no prior report of it (searched `DeepScence cross-sample
  score comparability sign flip`).
- Report the per-section CDKN1A correlation itself. Without it the claim is unfalsifiable and looks
  like a pipeline bug.

**Framing overall.** This finding's real value is that the project **caught and corrected its own
motivating claim before submission**. At ICBINB that is a virtue and should be said out loud; the
`SUBMISSION_PATCH` §"upgrade, not damage control" argument is correct and well made. At
ml4spatialbio it is a paragraph in Limitations at most — it is not spatial, not ML, and there are
only 4 pages.

### 2.5 Finding 5 — the central negative result

**Verdict: THE NEGATIVE CLAIM IS NOVEL AND WELL-POSITIONED. NOBODY HAS ESTIMATED A SENESCENCE
LENGTH CONSTANT WITH UNCERTAINTY, SO THERE IS NOTHING TO CONTRADICT — WHICH IS ALSO THE RISK.**

**Who has claimed to measure senescence spatial spread, and how you land against them:**

| Prior work | What they actually claim | How your result lands |
|---|---|---|
| **Ma S *et al.* (2024), Cell 187(24):7025–7044** (the "Zhao et al." of the brief) | Senescence-sensitive spots; SASP score peaks at SSS centre and "gradually decreased outward to peripheral cells in a distance-dependent manner"; TNF signalling decreases with distance; SSSs are "epicenters … compromising surrounding cells in a distance-dependent manner". **Stereo-seq spots/bins, 1,535,191 spots, ~1,450 genes/spot, two age groups.** No kernel, no length constant, no negative control, no permutation null that I could find. | **Direct engagement, and the resolution argument is your best card.** Their measurement is at bin resolution, which blurs the very quantity being estimated. You should say so plainly — the repo's own D7 correction already frames this well: "the field's closest distance-gradient result was measured at a resolution that blurs distance." You do not contradict their observation; you show it does not survive controls at single-cell resolution. |
| **Karpova/*Cell Genomics* 2026 liver (GSE310392)** | Verified in full text: senescent hepatocytes "were more often adjacent to one another than to normal hepatocytes"; clustered vs isolated differ in CXCL2/IGFBP2/SDC4 (Wilcoxon FDR<0.05); on the microenvironment — **"Although falling short of significance, we observed trends of increased content of HSCs and inflammatory macrophages"**; "possible paracrine dissemination of senescence". CellChat v2.1.2 + Squidpy co-occurrence; **no permutation/null model for spatial proximity that the paper describes.** | The plan's characterization is accurate. Your negative result is *consistent* with theirs — they saw a trend and declined to call it. Frame as: the analysis they said was missing, run properly, returns a bound. This is friendly, not adversarial, and it is the strongest possible position against the group whose data you use. |
| **scDOT** (*Genome Biology* 2024, doi:10.1186/s13059-024-03426-0) | Maps senescent cells and identifies neighbouring cells / candidate interaction genes. No decay parameter. | No conflict. |
| **COMMOT** (*Nat Methods* 2023), **SCILD** (*Commun Biol* 2026, doi:10.1038/s42003-025-09413-w — published 7 Jan 2026, integrates ligand diffusion, competitive binding and concentration decay), **NCEM** (*Nat Biotechnol* 2022; "length scales that are characteristic for known communication mechanisms") | All estimate distance decay generically. None reports calibration or a type-I-error check that I could find. | Your claim is about **identifiability of the estimand**, not about beating them. Plan §29 objection 2 already has the right answer; keep it verbatim. |
| **Martin, Schumacher & Chandra (2023), Aging Cell 22(8):e13892** | The containment paradox; explicitly declines to say whether spread is finite. | Your bound is a *constraint* on their model, which is the most useful thing a negative result can be. This is the discussion paragraph that makes the paper feel like biology rather than an audit. |
| **stAge** (bioRxiv 2025.11.23.689860, *Multi-tissue spatial transcriptomics reveals biological age hotspots*) | **Not in `references.bib` and it should be.** Localized transcriptomic age; "hotspots of accelerated aging and coldspots"; "robust spatial gradients of biological age" across mouse and human tissues. | The nearest contemporaneous competitor claim of the form "aging/senescence has spatial gradients". Add it to Related Work and say what your bound does and does not say about it. |

**SEARCHED, NONE FOUND:** any published length constant, decay rate, or λ with uncertainty for
senescence spatial influence (searches: `senescence length constant exponential decay distance
kernel fit spatial transcriptomics lambda micron`; `senescence spatial spread distance micrometers
measured quantify length scale`). The only quantitative distances in the literature are
**nearest-neighbour distances** (the 2026 human endometrium study, 45–211 µm by immune subtype,
already in §3 of the plan) and the **in vitro ~1 mm** spread from a seeded circle of senescent
cells (Nelson *et al.* 2012, *Aging Cell* 11(2):345, senescence-induced senescence; discussed in
Martin 2023). **Your §5 contribution-1 claim is therefore intact.**

**The structural risk of a negative result with no prior number to contradict.** You are not
overturning a published λ — there isn't one. So the paper's force rests entirely on (a) the
identifiability analysis showing the estimand is *not recoverable* rather than merely *not
recovered here*, and (b) the demonstration that existing tools would have returned something.
Figure 1 (synthetic identifiability map) and Figure 4 (COMMOT/SpaTalk/CellChat under the same
nulls) are therefore not supporting material — **they are what makes this a result rather than a
null experiment**. Prioritize them over anything else if time runs short.

**The 1 mm problem.** Nelson 2012's in vitro spread is ~1 mm; your fitting window is 100 µm and
your median λ̂ is 12.8 µm. A reviewer will ask whether you looked far enough. You need one sentence
stating the distance range over which the bound applies — "no effect above 0.203 response-sd at 80%
power **over 0–100 µm**" — and one sentence on why 100 µm (cytokine signalling is expected over
tens of µm; the endometrium nearest-neighbour distances are 45–211 µm; a 100 µm window with
λ̂ ≈ 13 µm is ~8λ). Without that, "we found nothing" is answerable with "you looked in the wrong
place."

---

## 3. Ranked by strength of contribution

Ranking is novelty × defensibility × transferability, not by how much of the paper each occupies.

| # | Finding | Why here |
|---|---|---|
| **1** | **The N2-vs-N5 result inside Finding 1** — a matched-decoy contrast does **not** remove a technical gradient that a covariate block does | Genuinely unreported, mechanistically explained (per-cell detection efficiency, which propensity matching on neighbourhood covariates cannot see), directly transferable to every imaging-ST distance analysis, and it **overturns the project's own prior** that N2 is the conservative gold standard. Highest surprise-per-word in the repo. Costs one table. |
| **2** | **Finding 5, the central negative result / bound** | The headline, and the reason the paper exists. Novel because no λ with uncertainty exists to contradict — which is also why it needs Figures 1 and 4 to carry it. Highly relevant to senescence biology (Martin's containment question), to SenNet, and to anyone building on Ma 2024. Defensibility depends entirely on the identifiability argument landing. |
| **3** | **Finding 2, torus degeneracy** | Practically valuable and vividly quantified, and the FFT enumeration is a real tool. Demoted from #1 because the statistical content is 1982-and-2021 vintage and the classical fix (variance correction) is better than the one the repo adopted. Reframed as an import + demonstration + tool it is strong; presented as a discovery it is a liability. |
| **4** | **Finding 1's diagnostic per se** (fit the kernel to negative controls) | Novel as a construction, but its two nearest neighbours (Voyager's Moran's I on negative controls; Lipsitch's negative-control outcome) are close enough that it cannot stand alone. It earns its place as the vehicle for #1 and #5. |
| **5** | **Finding 3, the 9–16% FPR** | Novel instrument, modest weight: self-calibration, few control probes, only powered when pooled, and two competing interpretations left undistinguished. One Methods sentence. Would rise a rank if you also ran COMMOT/CellChat/SpaTalk on the control probes — that turns a self-check into a benchmark. |
| **6** | **Finding 4, caller coverage sensitivity** | Headline is old (7 citations), the statistical phenomenon is older (Feinstein & Cicchetti 1990), and the specific measurement is currently confounded with the Tier A 25→33 gene-set change. Real value is as an honest self-correction. Not a contribution; a Limitations paragraph and, at ICBINB, a good anecdote. |

---

## 4. Anticipated reviewer objections, and the evidence each needs

| # | Objection | Evidence needed | Status in repo |
|---|---|---|---|
| O1 | *"Voyager already computes Moran's I on Xenium negative controls and concludes there is no technical spatial trend. What is new?"* | Report **your own Moran's I on the controls alongside the kernel amplitude**, and state that a near-zero global autocorrelation does not preclude a projection onto a specific covariate. Cite Voyager and Ren 2025 approvingly. | **NOT DONE.** Highest-priority gap. |
| O2 | *"Your −0.070 SD is just cell size / sequencing depth. You have rediscovered normalization."* | `neg_probe_rate` flat naively (+0.014, p=0.079) → the gradient is detection efficiency, not probe binding; and that is precisely why N2 misses it. | **AVAILABLE**, `CS_PHASE8_CALLERS.md` §4.1. Must be promoted into the main text, not left in an appendix. |
| O3 | *"Torus shifts on irregular windows are a solved problem — Lotwick & Silverman 1982, Mrkvička et al. 2021 variance correction. Why did you not use the standard fix?"* | Both citations, plus either an implemented variance-corrected N3 variant or an explicit justification for tiling. | **NOT DONE.** Second-highest priority; the variance-corrected run is cheap. |
| O4 | *"You cite CellWHISPER as the source of the torus-shift null. It permutes within cell type."* | Fix the citation everywhere (plan §32, Related Work, Figure 4 caption). | **NOT DONE.** Verified error. Cheap to fix, expensive to be caught on. |
| O5 | *"λ̂ = 12.8 µm is below your segmentation resolution and your window is 100 µm, but senescence spreads ~1 mm in vitro (Nelson 2012). You looked in the wrong place."* | State the bound's distance range explicitly; justify 100 µm against the endometrium 45–211 µm nearest-neighbour calibration and λ̂ ≈ 13 µm (≈8λ of window); ideally show the fit is stable to a wider window. | **PARTIALLY** — the calibration numbers are in plan §3; the range statement is not attached to the bound. |
| O6 | *"Your caller-agreement change is the gene set, not the coverage — your own patch file says Tier A went 25→33 and the numbers moved."* | 2-section vs 11-section agreement recomputed **on the frozen strict-33 Tier A on both**. | **NOT DONE / pending M1.** Until this exists, do not put the coverage claim in an abstract. |
| O7 | *"The DeepScence sign flip is the method's documented CDKN1A anchoring rule, not a finding."* | Per-section CDKN1A correlation; reframe as a cross-section comparability caveat for users. | **NOT DONE.** |
| O8 | *"This is a negative result about your estimator, not about biology."* | Figure 1 (identifiability regimes, synthetic ground truth) + Figure 4 (COMMOT/SpaTalk/CellChat under the same nulls) + CI-coverage on synthetic data (plan §24.7). | Fig 1 **DONE**; Fig 4 **DONE** but rests on the mis-cited null; coverage check status unclear. |
| O9 | *"Garden of forking paths — you ran a null battery and reported the one that gave zero."* | `PREREG_PHASE8.md` + `phase8-frozen` tag + the P2/P3 pre-registered statements + `CORRECTIONS.md`. | **STRONG.** This is the project's best structural defence; cite the pre-registration in the paper itself, not just the repo. |
| O10 | *"n = 2 donors over 55; one mouse tissue."* | D4 (age as continuous covariate, no stratified claim) is the right call; say the two-arm replication is a geometry replication, not an ageing result. | **DONE** as a decision; needs to be in the text. |
| O11 | *"Senescence-caller disagreement is well known — see DeepScence, SenCID, SenePy, ICE, markeR, Ntintas."* | Cite them. Claim only the coverage-sensitivity increment, and only once O6 is answered. | **NOT DONE** — `references.bib` has DeepScence/SenCID/SenePy/Ntintas but the framing still reads as if disagreement were a finding. ICE, markeR and the Jan-2026 circular-validation preprint are missing. |
| O12 | *"You have no spatial-statistics citations. Spatial confounding is a whole literature."* | Hodges & Reich 2010; Dupont, Wood & Augustin, *Spatial+*, Biometrics 78(4):1279 (2022); the RSR-critique line. | **NOT DONE.** Three citations, large credibility return. |

---

## 5. Venue fit, given a negative headline

**Verified facts.** ml4spatialbio 2026 (<https://imsb-uke.github.io/ml4spatialbio-2026/cfp.html>):
**4 pages maximum** excluding references and optional appendix, NeurIPS 2026 style,
**non-archival** (concurrent submission elsewhere explicitly permitted), double-blind, deadline
listed as **August 29 – September 4, 2026 AoE** — i.e. the CfP page now shows a window extending
to Sept 4, not a hard Aug 29. Topics verbatim include *"Modeling cell-cell communication and tissue
dynamics over time"*, *"Building interpretable and uncertainty-aware spatial models that biologists
can trust"*, and *"Designing benchmarks, datasets, and evaluation standards specific to spatial
tasks"*. **Negative results are not mentioned.**

ICBINB-BIO (<https://icbinb-bio.github.io/>): NeurIPS 2026, Sydney, Dec 11–12; deadline
**Aug 29, 2026 11:59 AoE**; **full (8 pages) and tiny (4 pages)** tracks. Topic list includes
*"Causal mechanisms versus spurious correlation"*, *"Failure under weak or confounded supervision"*,
*"Uncertainty, calibration, and decision-aware reliability"*, and *"Deployment-relevant evaluation
beyond benchmarks"*. Its stated wants: *"candid failure analysis, negative results, unexpectedly
strong simple baselines, and benchmarks that support the claims"*.

**Assessment. The primary/secondary assignment in the master plan is backwards for the result you
actually got.**

- The plan (§29) was written before the outcome was known and assigned ml4spatialbio primary. That
  was right for "the effect largely survives." It is not right for "naive 0.326 → controlled 0.027,
  nothing above a 0.203-sd bound."
- **ICBINB-BIO is a line-by-line match** for what you have: a confounded-supervision failure, a
  spurious-correlation-vs-causal-mechanism story, a calibration measurement (9–16% vs 5%), and a
  candid self-correction (Finding 4). It gives you **8 pages**, which is the only way all five
  findings plus the null battery plus two figures fit. A negative headline is the point of the
  venue, not a liability to be managed.
- **ml4spatialbio is still worth submitting to, as a different, shorter paper.** At 4 pages you
  cannot tell this story; you can tell **one** of them well. The right 4-page paper for that venue
  is *"Negative-control probes as a calibration instrument for distance-dependent spatial models"* —
  Finding 1 + Finding 3 + the tiling/torus methods note — which lands squarely on "designing
  benchmarks, datasets and evaluation standards specific to spatial tasks" and "interpretable and
  uncertainty-aware spatial models." That framing is positive-valence (here is a tool) even though
  the underlying result is negative, which is the correct register for a venue that does not
  advertise for negative results.
- ml4spatialbio is **non-archival and explicitly permits concurrent submission**, so submitting the
  4-page distillation there and the 8-page paper to ICBINB-BIO is allowed on ml4spatialbio's side.
  **I could not verify ICBINB-BIO's dual-submission policy** — check it before doing this.
- The Sept 4 ml4spatialbio window, if it holds, gives you six days after the ICBINB deadline to
  cut the 4-page version from the 8-page one. That ordering is strictly easier than the reverse.

**Recommendation:** ICBINB-BIO primary (8-page track); ml4spatialbio secondary as a 4-page
methods-instrument distillation. This inverts plan §29 and should be a deliberate, recorded PI
decision, not a drift.

---

## 6. Bibliography of everything retrieved for this assessment

Grouped by which finding it bears on. All URLs were fetched or returned by search on 2026-08-27.

**Negative controls, QC and calibration (Findings 1, 3)**
- Moses L, Einarsson PH, Jackson K, Luebbert L, Booeshaghi AS, Antonsson S, Bray N, Melsted P, Pachter L. *Voyager: exploratory single-cell genomics data analysis with geospatial statistics.* bioRxiv 2023.07.20.549945; PMC10461913. Xenium vignette: <https://pachterlab.github.io/voyager/articles/vig5_xenium.html> — **the key undercutting source.**
- Ren P, Zhang R, Wang Y, *et al.* *Systematic benchmarking of high-throughput subcellular spatial transcriptomics platforms across human tumors.* Nat Commun 16 (2025). doi:10.1038/s41467-025-64292-3; PMC12534522.
- Hallinan C, Ji HJ, Tsou E, Salzberg SL, Fan J. *Evidence of off-target probe binding affecting 10x Genomics Xenium gene panels compromise accuracy of spatial transcriptomic profiling.* eLife 14:RP107070 (2026). Confirms NCPs are **not** analysed spatially there: "relying solely on such probes for error detection may be insufficient."
- *SpaceTrooper: a quality control framework for imaging-based spatial omics data.* bioRxiv 2025.12.24.696336. Per-cell NCP proportion QC across Xenium/CosMx/MERFISH.
- *Orchestrating Spatial Transcriptomics Analysis with Bioconductor* (OSTA), QC chapter: <https://lmweber.org/OSTA/pages/img-quality-control.html>; bioRxiv 2025.11.20.688607.
- 10x Genomics, *Calculating Negative Control Metrics* and *Metrics Matter: false discovery rate and target specificity* (blog; **not retrieved**, 429).
- Lipsitch M, Tchetgen Tchetgen E, Cohen T. *Negative controls: a tool for detecting confounding and bias in observational studies.* Epidemiology 21(3):383–388 (2010).
- Gerard D & Stephens M / Wang J, Zhao Q, Hastie T, Owen AB — RUV / control-gene calibration line; PMC10751021, arXiv:1705.08393. Notes control-based calibration is anti-conservative with few controls.
- ResolVI (bioRxiv 2025.01.20.634005) — models unspecific background and the "diffusion phenomenon" in imaging ST; relevant as the tool that would *absorb* your confound.
- SpatialArtifacts (bioRxiv 2026.05.15.725260) — MAD + morphology artifact detection; does **not** appear to use NCPs (**COULD NOT VERIFY**, 429).

**Random-shift nulls (Finding 2)**
- Lotwick HW, Silverman BW. *Methods for analysing spatial processes of several types of points.* JRSS-B 44(3):406–413 (1982).
- Mrkvička T, Dvořák J, González JA, Mateu J. *Revisiting the random shift approach for testing in spatial statistics.* Spatial Statistics 42:100430 (2021); arXiv:1911.00240.
- Dvořák J, Mrkvička T, *et al.* Nonparametric testing of the dependence structure among points–marks–covariates (arXiv:2005.01019) and covariate significance under nuisance covariates (arXiv:2210.05424; JCGS 2024).
- `spatstat.random::rshift.ppp` documentation — torus edge correction requires a rectangular window.
- Rotation-invariance in SVG detection: PMC12318016; STANCE, Nat Commun (2025) doi:10.1038/s41467-025-57117-w.

**Spatial CCC methods, benchmarks, calibration (Findings 3, 5)**
- Kumar A, Rivera F, Aggarwal B, Zhang N, Coskun A, Sinha S. *CellWHISPER.* bioRxiv 2026.01.07.697982 (v2: *CellWHISPER disentangles direct cell-cell communication from structural proximity*); SSRN 6815051. **Null = within-cell-type permutation; no torus shift; no negative control probes.**
- *CONCISE: spatial co-expression and cell-cell communication inference.* bioRxiv 2026.06.22.733860 (**COULD NOT VERIFY** full text).
- *Benchmarking tools for deciphering cellular crosstalk in spatially-resolved transcriptomics.* Genome Biology (2026), doi:10.1186/s13059-026-04063-5. Nine CCI methods (**COULD NOT VERIFY** full text).
- Li Z *et al.* *SpatialDM.* Nat Commun 14 (2023), doi:10.1038/s41467-023-39608-w.
- Fischer DS, Schaar AC, Theis FJ. *NCEM.* Nat Biotechnol (2022), doi:10.1038/s41587-022-01467-z.
- Cang Z *et al.* *COMMOT.* Nat Methods (2023), doi:10.1038/s41592-022-01728-4.
- *SCILD: advancing spatial cellular communication inference with ligand diffusion and transport model.* Commun Biol (2026), doi:10.1038/s42003-025-09413-w, published 7 Jan 2026.
- Hodges JS & Reich BJ (2010); Dupont E, Wood SN, Augustin NH. *Spatial+.* Biometrics 78(4):1279 (2022) — the spatial-confounding literature the project does not cite.
- *Identifiability limits of physics-informed inference for spatial stochastic dynamics from static snapshots.* arXiv:2607.01749 (2026) — already in the plan; verified real, and its "distributed sources are non-identifiable, a point source restores identifiability" result is directly quotable for your §2.2 argument.

**Senescence (Findings 4, 5)**
- Ma S, Ji Z, Zhang B, *et al.* (47 authors), Liu G-H. *Spatial transcriptomic landscape unveils immunoglobin-associated senescence as a hallmark of aging.* Cell 187(24):7025–7044 (2024). doi:10.1016/j.cell.2024.10.019. Data: CNGB STOMICS STDS0000247.
- *Cellular senescence in human liver under normal aging and cancer.* Cell Genomics (2026), S2666-979X(25)00389-1; PMC12903365; GEO GSE310392. **Verified quotes** on clustering, "falling short of significance", "possible paracrine dissemination", CellChat v2.1.2 + Squidpy.
- Qu *et al.* *DeepScence.* Cell Genomics 5(12):101035 (2025). CoreScence = 39 genes from ≥5 of nine gene sets; **CDKN1A sign-anchoring documented**.
- *SenePy.* Nat Commun (2025), doi:10.1038/s41467-025-57047-7; PMID 39987255.
- *SenCID.* Cell Metabolism (2024), S1550-4131(24)00088-3. No feature gene shared by all six SID models.
- *ICE.* Genome Biology (2026), doi:10.1186/s13059-026-03997-0. **Not in references.bib.**
- *markeR.* NAR Genomics and Bioinformatics 8(2):lqag057 (2026). **Not in references.bib.**
- *Integrative transcriptomic identification of cellular senescence beyond marker limitations.* bioRxiv 2026.01.02.697374 — "opposing enrichment patterns … apparent concordance may reflect circular validation." **Not in references.bib; closest published statement to the project's own §0.3 thesis.**
- Nelson G *et al.* *A senescent cell bystander effect: senescence-induced senescence.* Aging Cell 11(2):345–349 (2012); PMC3488292. Source of the ~1 mm in vitro spread.
- Martin L, Schumacher L, Chandra T. *Modelling the dynamics of senescence spread.* Aging Cell 22(8):e13892 (2023); PMC10410058.
- *Multi-tissue spatial transcriptomics reveals biological age hotspots in mouse and human aging* (stAge). bioRxiv 2025.11.23.689860. **Not in references.bib; should be.**
- Feinstein AR, Cicchetti DV. *High agreement but low kappa: I. The problems of two paradoxes.* J Clin Epidemiol 43(6):543–549 (1990).

**Venues**
- ml4spatialbio 2026 CfP: <https://imsb-uke.github.io/ml4spatialbio-2026/cfp.html>
- ICBINB-BIO 2026: <https://icbinb-bio.github.io/>
