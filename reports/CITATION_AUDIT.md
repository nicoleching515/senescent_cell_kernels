# Citation audit

**Prepared:** 2026-08-27. **Scope:** `references.bib` (all entries), `SASP_Kernel_Master_Plan.md`
§4 / §29 / §31 / §32, and every report carrying the misattribution named below.
**Method:** every verdict here rests on a source retrieved in-session on 2026-08-27 —
Crossref, the NCBI PubMed E-utilities, the bioRxiv API, PMC, arXiv, or the publisher page.
Nothing was verified from recall. Where retrieval failed, the row says **UNVERIFIABLE** and
says why; it does not say CONFIRMED.

**Files changed:** `references.bib`, `SASP_Kernel_Master_Plan.md`, `README.md`,
`Phase7_Minimal_Human_Replication (1).md`, `reports/CS_PHASE1.md`, `reports/CS_PHASE2.md`,
`reports/CS_PHASE3.md`, `reports/SUBMISSION_PATCH_2026-08-29.md`, and this file.
**Not touched:** `results/`, `figures/`, `data/`, `genesets/`, `code/`. No commit, no tag.

---

## 0. Lead with these

### 0.1 The CellWHISPER torus misattribution — CONFIRMED as an error, now corrected

**CellWHISPER is not the source of the torus-shift null, and contains no torus shift.**

I retrieved the v1 full text independently (`biorxiv.org/content/10.64898/2026.01.07.697982v1.full`,
rendered 2026-08-27 after the earlier HTTP 429 block lifted). Verbatim:

> "To assess significance, CellWHISPER shuffles locations among cells of same type (preserving
> cell-type expression and spatial distributions), recomputes N_ijkl and uses the resulting null
> to compute a z-score."

and

> "permuting cell locations within each cell type, preserving cell-type-specific spatial
> organization and LR expression while destroying spatial proximity between ligand- and
> receptor-expressing cells."

The strings **torus**, **toroidal** and **wraparound** do not occur anywhere in the paper, and it
uses no negative control probes or blank barcodes. **Their null is this project's N1**, the
cell-type-stratified label permutation — not N3.

**The correct attribution for the torus-shift null, stated plainly:**

> The random-shift / toroidal null originates with **Lotwick & Silverman (1982)**, *Methods for
> analysing spatial processes of several types of points*, JRSS-B 44(3):406–413,
> doi:10.1111/j.2517-6161.1982.tb01221.x. Its rectangular-window requirement, its liberality, and
> the variance correction that replaces it on irregular windows are **Mrkvička, Dvořák, González &
> Mateu (2021)**, *Revisiting the random shift approach for testing in spatial statistics*,
> *Spatial Statistics* 42:100430, doi:10.1016/j.spasta.2020.100430.

This matters because **N3 is what Figure 4 rests on**. The figure itself is unaffected — the
statistics were never CellWHISPER's — but the sentence attributing the null to them was load-bearing
prose in six places and would have been caught by any reviewer who had read the preprint.

**Where it was, and what it now says:**

| File | Was | Now |
|---|---|---|
| `SASP_Kernel_Master_Plan.md` §22 Step 3 (N3) | "This is the null CellWHISPER showed leading methods fail." | Attributed to Lotwick & Silverman (1982) / Mrkvička et al. (2021), with the §2.1.4 quotation on tiling, plus a dated correction block |
| `SASP_Kernel_Master_Plan.md` §32 item 3 | "Your framing citation and the source of the torus-shift null." | "Its null is a within-cell-type location permutation — this project's N1. It is not the source of the torus-shift null and contains no torus shift at all" |
| `SASP_Kernel_Master_Plan.md` §31 ref 27 | "FPR >90% … under randomization" | Full author list, DOI, the verbatim FPR sentence, and the N1-not-N3 statement |
| `SASP_Kernel_Master_Plan.md` §23 baselines | torus-shifted run "independently reproduced the CellWHISPER finding" | Only the within-cell-type permutation (`N0_type`) replicates their design; the torus run reproduces nothing of theirs |
| `SASP_Kernel_Master_Plan.md` §4.3 | "benchmarked … on coordinate-randomized data" | "against a null that permutes cell locations within each cell type" |
| `reports/CS_PHASE1.md` §"not implemented" | "it is the null CellWHISPER showed leading methods fail" | corrected inline, dated |
| `reports/CS_PHASE2.md` finding 2 | "reproduces the CellWHISPER finding that existing tools fail the torus shift" | corrected inline, dated |
| `reports/CS_PHASE3.md` §on N3/N4 | same phrasing | corrected inline, dated |
| `references.bib` `kumar2026cellwhisper` | B9 note only | new `AUDIT (CIT-1)` block with the verbatim null sentence and the N1-not-N3 statement |

`reports/CS_PHASE4.md` and `README.md` were **already correct** on this point — the D7 §B9 pass had
caught the "coordinate randomisation" half of the error and the repo went on to actually run
`N0_type`. What D7 did **not** catch is that the plan simultaneously credited them with the *torus*
null. Two independent misreadings of the same paper, one fixed and one not.

### 0.2 The ">90% FPR" characterisation — no drift found in the reports, one drift in the plan

CellWHISPER's own sentence, now rendered and verified verbatim:

> "CellChat v2, COMMOT, and SpaTalk predict similar numbers of interactions on real (blue) and
> randomized (orange) data, indicating poor specificity and false positive rates (FPR) >90%. In
> contrast, CellWHISPER produced markedly fewer interactions on randomized data compared to real
> data, suggesting FPR < 5%."

So the ">90% FPR" phrasing is **theirs**, but the argument behind it is an **interaction-count
ratio between real and permuted input from which they infer an FPR** — not a measured type-I error
against a nominal level. The repo's characterisation:

- **Plan §1, §4.3, §31 ref 27** — use "implying". **Correct, kept.**
- **Plan §5 (the gap statement)** — read "the class of methods … **has** false-positive rates above
  90% under randomization". **This had drifted into asserting a measured rate.** Rewritten to
  "returns comparable interaction counts on permuted coordinates, implying false-positive rates
  above 90%", with a dated note.
- **Plan §2.2** — "Given the >90% FPR finding". Softened to the count-ratio formulation.
- **`CS_PHASE4.md` §5 item 7 and Figure 4's caption** already disclaim the conflation explicitly
  ("survival equals an FPR only if the shuffle is a true null — the assumption under test").
  **No change needed. That disclaimer is the best-written sentence on this topic in the repo and
  should travel to the paper verbatim.**

`README.md`'s "still outstanding" list said the >90% sentence had never been rendered from the
publisher. **That is now resolved** and the README says so.

### 0.3 The audit's own finding: 19 of 32 bib entries carry invented author forenames

This is the answer to Task 3's question *"is the error symptomatic?"*, and it is worse than the two
known cases suggested.

The D7 pass checked journal, volume, issue, pages and DOI against publisher records — and those are
overwhelmingly right. **It did not check forenames.** Comparing every entry against *both* the
Crossref deposit and the NCBI PubMed record (which agree in every case) gives **41 wrong given names
across 19 of the 32 entries.** Selected examples, all corrected:

| Key | Was | Is |
|---|---|---|
| `martin2023modelling` | Martin, **Luke** | Martin, **Lucy** |
| `kumar2026cellwhisper` | **Abhishek** Kumar; **Fernando** Rivera Moctezuma; Bhav**ya** Aggarwal; **Nan** Zhang | **Anurendra** Kumar; **Felix** Rivera; Bhav**ay** Aggarwal; **Nicholas** Zhang |
| `yao2026benchmarking` | **Hong** Yao; **Shuang** Mu; **Fang** He; **Zhixiang** Fang | **Hui** Yao; **Shuai** Mu; **Fei** He; **Zhaoyuan** Fang |
| `ntintas2026overview` | **Odysseas** Ntintas; **Sofia** Vagena; **Panagiotis** Pantelis; **George** Theocharous; **Russell** Petty | **Orestis**; **Sylvia**; **Pavlos**; **Giorgos**; **Russel** |
| `zhao2026concise` | **Jinyu** Zhao; **Xiaoyu** Shan; **Guanjue** Wang; Rui **B.** Chang | **Jia** Zhao; **Xinning** Shan; **Gefei** Wang; **Rui** Chang |
| `qu2025deepscence` | **Yi** Qu; **Bing** Ji; **Rui** Dong | **Yilong** Qu; **Beijie** Ji; **Runze** Dong |
| `russo2026characterizing` | **Tommaso** Russo | **Taylor** Russo |
| `meyer2017model` | **Pascal** Meyer; **Andreas** Burkovski | **Patrick** Meyer; **Andre** Burkovski |
| `gu2026identifiability` | **Ruiyang** Gu; **Ruo-Zhen** Zhang | **Rujie** Gu; **Ray Zirui** Zhang |
| `wang2026harmonic` | **Xiao** Wang; **Chenyu** Tao; **Yun** Jiang | **Xiaofei** Wang; **Chenyang** Tao; **Yixuan** Jiang |
| `ma2024spatial` | Zhang, **Beibei**; Gu, **Yingjie** | Zhang, **Bin**; Gu, **Ying** |
| `yu2026scild` | **Jiaxin** Yu; **Jun** Zhao | **Jiating** Yu; **Jinyue** Zhao |
| `gong2026rgast` | **Yuxuan** Gong; **Zhaoyang** Yu | **Yuqiao** Gong; **Zhangsheng** Yu |
| `shao2022spatalk` | Xu, **Xiaoyan** | Xu, **Xiao** (the bib had conflated Xiaoyan **Lu**, author 4, with Xiao **Xu**, author 11) |
| plus | `karpova` (Cheng-Wei→**Chien-Wei** Peng), `farzad` (Yang→**Yao** Lu), `suryadevara` (Mei→**Mingyu** Yang), `parvanov` (Magdalena→**Margarita** Ruseva), `sanborn` (Shanshan→**Shang** Gao) | |

**Diagnosis.** The forenames are plausible, correctly-transliterated, wrong. That is the signature
of a bibliography written from recall and then spot-checked on the fields a DOI resolver returns.
The Zhao-vs-Ma error and the CellWHISPER error are the same failure surfacing on the two fields
that happened to be load-bearing. **The remedy is procedural: copy author names from the Crossref
or PubMed record; never type them.** A note to that effect is now in the `references.bib` header.

Clean on authors: `gurkar`, `nguyen`, `anerillas`, `saul`, `tao`, `hughes`, `basisty`,
`acosta2013complex`, `acosta2026correction`, `fischer`, `cang`, `neretti`, `borner`.

---

## 1. What was added: the spatial-statistics block

`references.bib` previously contained **zero** spatial-statistics methods papers — verified by
searching the file for Lotwick, Silverman, Mrkvička, Hodges, Reich, Dupont, Baddeley and Ripley
before this pass. Eleven entries were added in a new block at the end of the file. Each carries a
`% SUPPORTS:` line naming the specific claim it is for, and a `% VERIFICATION:` line naming what was
retrieved and what was not.

| Key | Reference | What claim in this project it supports |
|---|---|---|
| `lotwick1982methods` | Lotwick HW, Silverman BW (1982), *JRSS-B* 44(3):406–413, doi:10.1111/j.2517-6161.1982.tb01221.x | **The origin of N3.** This is what the torus-shift null cites instead of CellWHISPER. Without it, §22 Step 3 has no attribution at all. |
| `mrkvicka2021revisiting` | Mrkvička T, Dvořák J, González JA, Mateu J (2021), *Spatial Statistics* 42:100430 | (a) Torus correction is **liberal** — it "makes a crack in the autocorrelation structure". (b) It "assumes a rectangular observation window". (c) **The direct objection to the project's tiling choice**: extending it "to windows which are finite unions of (aligned) rectangles … would increase the amount of cracks in the autocorrelation structure. Subsequently, it would increase the liberality of the test." (d) The **variance correction** — statistic computed on W ∩ (W+v), then standardized because different shifts retain different amounts of data — which "can be applied in case of general (compact) observation windows". Answers reviewer objection O3, and reframes CS_PHASE7_C1 §3 from discovery to first-quantification-on-real-tissue. |
| `spatstat2026rshift` | `spatstat.random::rshift.ppp` reference manual | The one-line software statement of the N3-occ problem: *"The window must be a rectangle. Toroidal shifts are undefined if the window is non-rectangular."* Cheapest possible citation for the constraint. |
| `baddeley2015spatial` | Baddeley A, Rubak E, Turner R (2015), *Spatial Point Patterns*, CRC, doi:10.1201/b19708 | The standard textbook behind `spatstat`. Included for completeness; **cite the `rshift` manual, not the book**, for the rectangular-window sentence — the book's own text was not consulted. |
| `hodges2010adding` | Hodges JS, Reich BJ (2010), *The American Statistician* 64(4):325–334 | **Names the problem this project is an instance of.** A spatially structured exposure (distance-to-nearest-sender) collinear with a spatially structured random field is textbook spatial confounding. Cite in Methods where N5/N6 are introduced. |
| `dupont2022spatialplus` | Dupont E, Wood SN, Augustin NH (2022), *Biometrics* 78(4):1279–1290 | The modern remedy, and the sentence that makes N6 legible to a statistician: Spatial+ "reduces the sensitivity of the estimates to smoothing by replacing the covariates by their residuals after spatial dependence has been regressed away." That residualise-the-covariate move **is** N6. |
| `khan2022restricted` | Khan K, Calder CA (2022), *JASA* 117(537):482–494 | Why the project does **not** simply orthogonalise the kernel against the spatial field: "RSR methods will typically perform worse than nonspatial methods", and the problem "cannot be fixed with a selection of 'better' spatial basis vectors". One Limitations sentence. |
| `zimmerman2022deconfounding` | Zimmerman DL, Ver Hoef JM (2022), *The American Statistician* 76(2):159–167 | The frequentist half of the same Limitations sentence. Cite with Khan & Calder, not alone. |
| `lipsitch2010negative` | Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010), *Epidemiology* 21(3):383–388 | **The correct framing for A7.** A7 refits the identical distance-to-sender kernel with control-probe counts as the response — a **negative control outcome for the estimand itself**. This paper formalises the construction ("We distinguish 2 types of negative controls (exposure controls and outcome controls) … and identify the conditions for the use of such negative controls to detect confounding"). It is also what separates A7 from Voyager/Ren, who compute a *generic* statistic on controls rather than refitting the model under test. |
| `moses2023voyager` | Moses L, … Pachter L (2023), bioRxiv 2023.07.20.549945 — **still a preprint** | The prior art that **falsifies** "nobody in this literature reports it". The Xenium vignette: "generally the negative controls are tightly clustered around 0, while the real genes have positive Moran's I, which means there is generally no technical artifact spatial trend." Must be cited approvingly, with our own Moran's I reported alongside the kernel amplitude. |
| `ren2025systematic` | Ren P *et al.* (2025), *Nat Commun* 16, doi:10.1038/s41467-025-64292-3 | The peer-reviewed instance of the same: "Spatial autocorrelation analysis using Moran's I revealed stronger aggregation of negative control signals in CosMx 6K." |

Corresponding entries 31–38 were added to master plan §31, spatial statistics was added to the
Related Work paragraph in §30, and three new reviewer objections (7, 8, 9) were added to §29 —
the torus/variance-correction objection, the spatial-confounding objection, and the Voyager
objection, each with the evidence it needs.

**Considered and not added:** Hanks *et al.* (2015), *Environmetrics* 26(4):243–254,
doi:10.1002/env.2331 — the geostatistical (continuous-space) RSR critique, which is arguably closer
to this project's setting than the areal papers. Bibliographic record verified; **its abstract is
not available from Crossref or Semantic Scholar and I could not render the full text**, so it is
named here rather than added. Add it if a human can read it.

---

## 2. Per-citation verdicts

Verdict applies to **what the repo cites the work for**, not merely to whether the work exists.

- **CONFIRMED** — the cited work says what it is cited for, verified from a retrieved source.
- **MISATTRIBUTED** — it does not.
- **UNVERIFIABLE** — the source could not be retrieved; stated as unverified, not assumed correct.

Bibliographic corrections applied in this pass are noted in the last column; they do not change
the verdict unless the attributed *claim* was wrong.

### 2.1 Load-bearing for a figure or a contribution claim

| Key | Cited for | Verdict | Note |
|---|---|---|---|
| `kumar2026cellwhisper` | **the torus-shift null (N3)** — Figure 4's null | **MISATTRIBUTED** | Their null is a within-cell-type permutation = **N1**. No torus anywhere in the paper. **Corrected everywhere; see §0.1.** |
| `kumar2026cellwhisper` | ">90% FPR for CellChat v2 / COMMOT / SpaTalk" | **CONFIRMED, with a caveat that must travel** | Full text rendered 2026-08-27; sentence verbatim. But it is an **interaction-count ratio**, not a measured type-I error. Plan §5 had drifted; fixed. |
| `kumar2026cellwhisper` | author names | — | 4 of 6 forenames were wrong. Corrected. |
| `ma2024spatial` | "closest prior work", the distance gradient | **CONFIRMED at abstract level; the finer claim UNVERIFIABLE** | PubMed abstract verified verbatim: "nine tissues in male mice", and SSSs "serve as epicenters for heightened inflammation that compromises surrounding cells in a distance-dependent manner." The per-pathway monotone-with-distance claim (SASP score, TNF, ATP biosynthesis, cell cycle) is **still not rendered** — cell.com/sciencedirect return 403, not in PMC or Europe PMC. **A human must open the PDF.** |
| `ma2024spatial` | cited project-wide as "Zhao et al." | **MISATTRIBUTED (already known, incompletely fixed)** | D7 fixed the bib and §4.1. **Four more instances survived** in §23, §29 objection 1, §30 Related Work and §31 ref 1. All now corrected. |
| `gu2026identifiability` | identifiability of distributed sources | **CONFIRMED (existence + title + authors)** | arXiv:2607.01749, 2026-07-02, verified from the arXiv abs page metadata. Forenames were wrong (Rujie Gu; Ray Zirui Zhang); corrected. The arXiv API returned rate-limited/empty, so the *content* claim rests on the plan's earlier reading, not on a rendering in this pass. |
| `nguyen2024scdot` | maps senescent cells, no decay parameter | **CONFIRMED** | Record and abstract verified. No conflict with our claim. |
| `martin2023modelling` | the containment paradox | **CONFIRMED** | The D7 §B4b correction (four mechanisms not in the paper) is intact everywhere I checked; no drift. **First author is Lucy Martin, not Luke.** Corrected. |

### 2.2 Assertions about what a named tool does

| Key | Cited for | Verdict | Note |
|---|---|---|---|
| `qu2025deepscence` | sender caller; CoreScence = 39 genes from ≥5 of nine sets; CDKN1A sign anchoring | **CONFIRMED** | Verbatim from the paper: "DeepScence calculates the Pearson correlation between the output senescence score and the expression of the CDKN1A gene … If the correlation is negative, the sign of the output score is flipped." So the observed polarity flip is **documented method behaviour**, not a defect — frame it as a cross-section comparability caveat. |
| `cang2023commot` | optimal-transport CCC, distance decay | **CONFIRMED** | Record verified. The `np.Inf` / numpy-2 import note is a repo fact, unaffected. |
| `fischer2023ncem` | niche-effect estimate "with a length scale" (plan §23) | **PARTIALLY CONFIRMED** | The abstract confirms a GNN-based niche-composition model and nothing about length scales. The "characteristic length scales" phrasing is in the body, which I did **not** render. Do not quote a length-scale claim from NCEM without opening it. |
| `shao2022spatalk` | knowledge-graph CCC baseline | **CONFIRMED** | Abstract verified. One forename corrected. |
| `yu2026scild` | "fits explicit ligand diffusion, competitive binding, and concentration decay in one optimization" | **CONFIRMED verbatim** | "SCILD integrates ligand diffusion, competitive ligand-receptor binding, and concentration decay into a unified optimization model." Two forenames corrected. |
| `wang2026harmonic` | "adds H&E histology to condition on tissue context and reduce false positives" | **CONFIRMED** | bioRxiv abstract: integrates H&E; "Significant refinement of false-positive/negative predictions was observed compared to ST-only CCC tools." **Title was wrong** — the deposited title has no "HARMONIC:" prefix; corrected. Three forenames corrected. |
| `gong2026rgast` | listed under "Spatial CCC methods" | **MISATTRIBUTED (D7 §S13, previously flagged, now applied)** | The published paper is representation learning benchmarked on spatial **domain identification**, with CCC as one downstream task. §31 ref 25 still carried a title/year pairing that never existed; replaced with the published record. |
| `zhao2026concise` | type-I inflation at a = 0.1 across "every competing spatial ligand-receptor method" | **CONFIRMED for their actual competitors; MISATTRIBUTED if read as the LR-tool family** | Their competitors are MERINGUE, SpatialDM, Copulacci, LIANA+ — co-expression/bivariate methods, not CellChat/COMMOT/SpaTalk. §4.3 now says so. Four forenames corrected. |
| `yao2026benchmarking` | "rotation invariance unresolved **and cell-type confounding a central open challenge**" | **MISATTRIBUTED in the second clause (D7 §S17, previously flagged, now applied)** | Rotation invariance is exact. Cell-type confounding is a documented weakness of **one** method (Celina), not an authorial open-challenge claim. §4.3 corrected. Four forenames corrected. |
| `tao2024sencid` | caller disagreement (no shared feature gene across six models) | **CONFIRMED** | Record verified; **DOI was missing** and has been added (10.1016/j.cmet.2024.03.009). |
| `sanborn2025senepy` | cell-type-specific senescence landscape | **CONFIRMED** | Record verified. One forename corrected. |
| `hughes2025senpred` | senescence classifier | **CONFIRMED** | The D7 "title incomplete" flag is **resolved** — full title recovered and inserted. |
| `anerillas2026sencat` | "14 primary human cell types across 30+ senescence paradigms, no universal marker" | **CONFIRMED verbatim** | "we profiled the transcriptomes and proteomes in 14 different primary human cell types undergoing over 30 senescence paradigms … senescent cells from all primary cell types did not share a single unique marker." |
| `saul2022senmayo` | SenMayo, 125 genes = 83+20+22 | **CONFIRMED** | Unchanged from D7. |
| `ntintas2026overview` | signature pros-and-cons overview | **CONFIRMED** | Abstract verified. **Five of seven forenames were wrong**; corrected. |

### 2.3 Assertions of priority or novelty

| Claim | Where | Verdict | Action taken |
|---|---|---|---|
| "Nobody in this literature reports [the negative-control-probe spatial diagnostic]" | `SUBMISSION_PATCH_2026-08-29.md` §6; `Phase7_Minimal_Human_Replication (1).md` A7 row | **MISATTRIBUTED / FALSE** | Voyager's Xenium vignette and Ren *et al.* 2025 both compute Moran's I on negative controls; both quotations retrieved and verified. **Struck in both files** with the narrower surviving claim substituted (refitting *the estimand's own estimator* to the controls — a Lipsitch negative-control **outcome**). New §29 objection 9 carries the defence. |
| "Nobody has published a length constant for senescence spatial influence" | plan §5 contribution 1 | **NOT CONTRADICTED** | I found no counterexample. Reported as *searched, none found*, not as verified-absent. |
| The N3 degeneracy on non-convex tissue is a finding in its own right | `CS_PHASE7_C1.md` §6.3 | **PRIOR EXISTS** | Not corrected here (`CS_PHASE7_C1.md` was outside the files I touched), but Mrkvička et al. §2.1.4 predicts exactly this, and the master plan §22 and §29 now carry the citation and the objection. **Reframe §6.3 before submission: first quantification on real tissue, not discovery.** |
| CellWHISPER is the source of the torus-shift null | plan §22, §23, §32; three phase reports | **MISATTRIBUTED** | See §0.1. |

### 2.4 The rest

| Key | Verdict | Note |
|---|---|---|
| `gurkar2023spatial` | **CONFIRMED** | SenNet review; record clean, authors clean. |
| `karpova2026cellular` | **CONFIRMED** | Primary dataset source. The D7 accession warning (GSE310392 is the **mouse** arm) stands. One forename corrected. |
| `farzad2026spatial` | **CONFIRMED** | Ages 18–86 and the interfollicular→germinal-centre shift both verified. **But the plan called it "the June 2026 *Cell* package"** — it is *Cell Press Blue* 1(4):100053; the *Cell* item in that package is Suryadevara et al. §4.1 corrected. One forename corrected. |
| `suryadevara2026charting` | **CONFIRMED, and the D7 flag is RESOLVED** | PubMed types it `Journal Article` + **`Review`**, and the abstract is four sentences with no data. **It is a Review. Do not cite it as primary evidence.** One forename corrected. |
| `parvanov2025spatial` | **CONFIRMED** | The five nearest-neighbour distances stand. **Title was incomplete** — the published title ends "…During the Implantation Window"; corrected. One forename corrected. |
| `basisty2020proteomic` | **CONFIRMED** | SASP Atlas. Clean. |
| `acosta2013complex` | **CONFIRMED for the transwell claim** | The D7 note stands. **§31 ref 17 gave the title as "Paracrine senescence transmission", which is not this paper's title** — corrected to the published title with the DOI. |
| `acosta2026correction` | **UNVERIFIABLE (content)** | Existence re-confirmed via Crossref and PubMed (25 authors, *Nat Cell Biol* 28(6):1343). **Content still unread** — nature.com redirects anonymous fetches to `idp.nature.com`. A human with Nature access must confirm it does not touch the transwell result. |
| `russo2026characterizing` | **CONFIRMED** | Peer-reviewed *Aging Cell* 25(8):e70673. **§31 ref 18 still cited the February 2026 preprint**; corrected. First author is **Taylor** Russo, not Tommaso. |
| `meyer2017model` | **CONFIRMED** | Record clean. Two forenames corrected. |
| `neretti2024dissecting` | **MISATTRIBUTED (D7 §B8, still unresolved)** | Verified again: *Innovation in Aging* 8(Supplement_1):351–351, **single author, conference abstract**. It remains the only support for the primary-vs-secondary senescence distinction underwriting Tier B module B7. **Replace or drop the claim** — this is a live liability, not a formatting nit. |
| `borner2026sennet` | **MISATTRIBUTED (title) + UNVERIFIABLE (the number)** | (a) The title was recorded as "The SenNet Data Portal"; the deposited title is **"SenNet Portal: Build, Optimization and Usage"**. Corrected. (b) **The "1,753 datasets … as of January 2026" figure could not be verified.** The bioRxiv API abstract for **both** v1 (2026-02-10) and v2 (2026-05-04) reads "As of April 2026, the portal hosts **2,041** publicly available human and mouse datasets across 15 organs using 6 general assay types." No retrievable version states 1,753 or January 2026. §4.1 now cites 2,041/April 2026 with a dated correction. |
| `yao2026benchmarking`, `zhao2026concise`, `gong2026rgast` | see §2.2 | |

---

## 3. What I could not verify — stated as unverified

1. **Ma et al. (2024), the per-pathway distance-gradient sentence.** Abstract-level claim verified;
   the figure-panel claim (SASP score / TNF / ATP biosynthesis / cell-cycle genes all monotone with
   distance from SSSs) is **not rendered**. cell.com and sciencedirect return 403; not in PMC or
   Europe PMC. *A human must open this PDF through the Caltech proxy before camera-ready.*
2. **Acosta et al. 2026 Author Correction, content.** Exists; unread. nature.com auth redirect.
3. **Hodges & Reich (2010), primary text.** Bibliographic record verified from Crossref; the
   abstract and full text were **not** rendered (Taylor & Francis 403; Semantic Scholar holds no
   abstract for the DOI). Its content is confirmed only indirectly, from Dupont et al. (2022)'s
   citation of it. Flagged as such inside `references.bib`.
4. **Zimmerman & Ver Hoef (2022), primary text.** Record verified; content from publisher-page and
   Semantic Scholar summaries only. Cite alongside Khan & Calder (whose abstract *was* retrieved),
   not alone.
5. **Lotwick & Silverman (1982), primary text.** Record verified from Crossref and the Wiley
   landing page. The *toroidal shift* attribution is verified **via Mrkvička et al. (2021) §1**,
   which quotes p.410 of the original directly. The JRSS-B full text itself is behind JSTOR/Wiley
   and was not rendered. The published abstract mentions only second-moment estimation and
   empty-space methods — **do not attribute anything beyond the toroidal shift from the abstract**.
6. **Hanks et al. (2015).** Record verified; no abstract available anywhere I could reach. Named in
   §1 above rather than added to the bib.
7. **The SenNet "1,753 / January 2026" figure.** See §2.4. Not in any retrievable version.
8. **NCEM's "characteristic length scales" phrasing.** In the body, not the abstract; not rendered.
9. **`gu2026identifiability` content.** Title, authors and date verified from the arXiv abs page;
   the abstract text was not retrieved this pass (arXiv API rate-limited).

One phrasing caution that is not an error but will become one if copied: **`NOVELTY_ASSESSMENT.md`
§U2 quotes Mrkvička et al. as saying the variance correction "can also be used for irregular
windows".** That is a paraphrase. The paper's actual wording is *"can be applied in case of general
(compact) observation windows"*. Quote the paper, not the summary. Noted in `references.bib`.

---

## 4. Recommended next actions, in order of exposure

1. **Reframe `CS_PHASE7_C1.md` §6.3** from discovery to first-quantification, with the Lotwick /
   Mrkvička citations attached. That file was outside this pass's edit scope. It is now the last
   place the torus finding is presented as novel statistics.
2. **Replace or drop `neretti2024dissecting`.** A GSA conference abstract cannot be the sole support
   for a Tier B module. This has been flagged since D7 and has not moved.
3. **Run the Mrkvička variance-corrected shift as an N3 variant**, or write the sentence explaining
   why tiling is preferred — noting that Mrkvička et al. §2.1.4 specifically predicts that a
   union-of-aligned-rectangles construction *increases* liberality. If it returns SF ≈ 0.95–1.00
   like N3-tile does, Contribution 3 becomes very hard to attack.
4. **Report our own Moran's I on the negative control probes** next to the kernel amplitude, and
   cite Voyager and Ren approvingly. This is the single highest-value unbuilt piece of evidence in
   the repo, and the sentence it defends is currently false.
5. **Library access, one sitting:** Ma et al. figure panel, and the Acosta 2026 Author Correction.
6. **Procedural:** copy author names from Crossref/PubMed. Never type them.
