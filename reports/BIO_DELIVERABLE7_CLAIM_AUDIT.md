# BIO DELIVERABLE 7 — CLAIM AUDIT

**Biology collaborator · 2026-08-21 · SASP Spatial Response Kernel**
Master plan §12 Deliverable 7: *"Read every biological sentence against a source.
Verify every citation resolves and says what the manuscript says it says. Flag any
correlational result described in causal language."*

**This is an audit, not a fix.** Nothing in `/workspace/reports/` or `README.md` was
edited. Every item below is a recommendation for the CS lead to accept or reject.

**Epistemic convention used throughout, and please hold me to it:**

| Tag | Means |
|---|---|
| **VERIFIED** | I fetched a publisher/index record and it matches. |
| **ERROR** | I fetched a record and it *contradicts* what we wrote. This is a positive finding of wrongness. |
| **UNVERIFIED** | I could not find a record. **This is not evidence the citation is wrong.** Search queries are listed in §3 so someone can repeat them. |
| **PARTIAL** | Bibliographic details check out; the *attributed claim* does not, or could not be reached. |

---

## 1. Reference-by-reference verification

**Scoreboard for §31's 30 references.**

| Outcome | Count | Refs |
|---|---|---|
| **VERIFIED** — exists; bibliography correct as given; attribution correct | **12** | 2, 4, 5, 6, 8, 11, 13, 14, 21, 23, 26, 30 |
| **INCOMPLETE, not wrong** — truncated/paraphrased title or missing journal-level fields; must be completed before submission | **5** | 12, 16, 22, 24, 28 |
| **ERROR — bibliographic** | **10** | 1, 3\*, 7, 9, 10, 17, 18, 19, 20, 25 |
| **ERROR — attributed claim** (bibliography fine, the paper does not say what we say it says) | **3** | 15, 27, 29 |
| **UNRESOLVED — could not find at all** | **0** | — |
| **Total** | **30** | |

**13 of 30 carry a genuine error; 12 verified cleanly; 5 are merely incomplete; none is
missing.**

\* Ref 3's *bibliography* is fully verified; the ERROR is the accession attribution that
this project already discovered and that motivated this audit.

**Every reference in §31 exists.** Nothing in the list is fabricated. That is worth
stating plainly, because the GSE310392 episode reasonably created a suspicion of
invention, and invention is not what went wrong. What went wrong is **author
attribution, publication year, document type, and — three times — what the paper
actually says**.

### 1.1 Reference-by-reference

Legend: **V** = verified · **E** = error I confirmed · **P** = partial (bibliography
fine, attribution off) · **I** = incomplete but not wrong.

| # | Citation as given in §31 | ? | What is wrong | Corrected form |
|---|---|---|---|---|
| 1 | Zhao et al. (2024). *Spatial transcriptomic landscape unveils immunoglobin-associated senescence…* Cell. S0092-8674(24)01201-7 | **E** | **Wrong first author.** First author is **Ma S**; "Zhao L" is 14th of ~50. Platform is **Stereo-seq, spot-based**, not single-cell. Data are at **CNGB/STOmicsDB**, not GEO/GSA. Two age groups, not "life stages". Title spelling "immunoglob**in**" is correct as published — do not fix it. | Ma S, Ji Z, Zhang B, … Qu J, Zhang W, Gu Y, Liu GH. *Cell.* 2024;187(24):7025–7044.e34. doi:10.1016/j.cell.2024.10.019. PMID 39500323. Data: STOmicsDB STDS0000247. |
| 2 | *Spatial mapping of cellular senescence…* (2023). Nature Aging. doi:10.1038/s43587-023-00446-6 | **V** | Nothing. Is a SenNet review, as described. Fields missing only. | Gurkar AU, Gerencser AA, Mora AL, … Passos JF. *Nat Aging.* 2023;3(7):776–790. PMID 37400722. |
| 3 | *Cellular senescence in human liver under normal aging and cancer* (2026). Cell Genomics. S2666-979X(25)00389-1; **GEO GSE310392** | **E** | Bibliography **verified**; 43 livers, 24 CRC mets, Xenium + multiome + CODEX all confirmed. The **accession attribution is wrong**: GSE310392 is the **mouse** Xenium arm. Human data are on the SenNet portal (WashU provider group); mCRC on HTAN. All three GEO accessions in the paper (GSE293958, GSE310392, GSE311064) are mouse. | Karpova A, Li X, Peng CW, … Fields RC, Ding L. *Cell Genom.* 2026;6(2):101133. doi:10.1016/j.xgen.2025.101133. PMID 41576948. PMC12903365. |
| 4 | scDOT (2024). Genome Biology. doi:10.1186/s13059-024-03426-0 | **V** | Nothing. | Nguyen ND, Rosas L, Khaliullin T, … Lugo-Martinez J, Bar-Joseph Z. *Genome Biol.* 2024;25(1):288. PMID 39516853. |
| 5 | Farzad et al. (2026). *A spatial multi-omics atlas of immunosenescence…* **Cell Press Blue** 1(4):100053 | **V** | Nothing — **"Cell Press Blue" is a real journal** (Cell Press/Elsevier, ISSN 3051-3839). Volume/issue/article all correct. Ages 18–86 and the extrafollicular→germinal-centre shift both confirmed. | Farzad N, Enninful A, Lu Y, et al. *Cell Press Blue.* 2026;1(4):100053. doi:10.1016/j.cpblue.2026.100053. |
| 6 | Suryadevara, Farzad et al. (2026). *Charting human cellular senescence in aging and disease.* Cell. S0092-8674(26)00587-8 | **V** | Nothing. But it is a **5-page item** (189(12):3501–3505) — check whether it is a research article or a Perspective before citing it as primary evidence. | Suryadevara V, Farzad N, Yang M, … Robbins P, Fan R. *Cell.* 2026;189(12):3501–3505. doi:10.1016/j.cell.2026.05.028. PMID 42276030. |
| 7 | *Spatial distribution of senescent cells … human endometrium* (2026) | **E** | **Year is 2025, not 2026.** Journal omitted entirely: *Diagnostics (Basel)* (MDPI). Title truncated. **All five distance figures verified verbatim.** Caveats: p16 IHC on serial sections, not transcriptomic; mid-luteal implantation-window endometrium from 68 IVF patients. | Parvanov D, Ganeva R, Ruseva M, et al. *Diagnostics (Basel).* 2025;15(21):2679. doi:10.3390/diagnostics15212679. PMID 41225972. |
| 8 | Qu et al. (2025). DeepScence. Cell Genomics 5(12):101035 | **V** | Nothing. Defines **CoreScence** (39 genes) and surveys **nine** senescence gene sets with documented disagreement, as attributed. | Qu Y, Ji B, Dong R, … Nixon AB, Ji Z. *Cell Genom.* 2025;5(12):101035. PMID 41061702. |
| 9 | Anerillas et al. (2026). *SenCat.* Molecular Cell. doi:10.1016/j.molcel.2026.05.017 | **E** | Title is not "SenCat" — that is the resource name. **14 cell types / >30 paradigms / no shared marker all verified.** | Anerillas C, Altés G, Gresova K, … Basisty N, Gorospe M. *SenCat: Cataloging human cell senescence through multi-omic profiling of multiple senescent primary cell types.* Mol Cell. 2026;86(13):2605–2616.e8. PMID 42276073. |
| 10 | SenePy (2025). Nature Communications. PMID 39987255 | **E** | Title truncated — ends *"…using SenePy."* Otherwise correct. | Sanborn MA, Wang X, Gao S, Dai Y, Rehman J. *Nat Commun.* 2025;16(1):1884. doi:10.1038/s41467-025-57047-7. |
| 11 | Saul et al. (2022). SenMayo. Nature Communications. doi:10.1038/s41467-022-32552-1 | **V** | Nothing. **"125 genes, 83 SASP, 20 transmembrane, 22 intracellular" is exactly the paper's own breakdown** (83+20+22 = 125). Independently confirmed. | Saul D, Kosinsky RL, Atkinson EJ, … Farr JN, Khosla S. *Nat Commun.* 2022;13(1):4827. PMID 35974106. |
| 12 | SenCID (2024). Cell Metabolism. S1550-4131(24)00088-3 | **I** | Journal and year correct; **no title or authors were given at all**. | Tao W, Yu Z, Han JJ. *Single-cell senescence identification reveals senescence heterogeneity, trajectory, and modulators.* Cell Metab. 2024;36(5):1126–1143.e5. PMID 38604170. |
| 13 | Hughes et al. (2025). SenPred. Genome Medicine. doi:10.1186/s13073-**024**-01418-0 | **V** | **The year/DOI mismatch is NOT an error.** BMC/Springer DOI suffixes encode submission year. 2025 is correct. Title incomplete. | Hughes BK, Davis A, Milligan D, … Bishop CL. *Genome Med.* 2025;17(1):2. PMID 39810225. |
| 14 | Ntintas et al. (2026). FEBS Open Bio. doi:10.1002/2211-5463.70134 | **V** | Nothing — surname spelling, year, and the **nine resources / 22-gene core** all correct. | Ntintas OA, Vagena S, Pantelis P, Theocharous G, Petty R, Evangelou K, Gorgoulis VG. *FEBS Open Bio.* 2026;16(5):821–836. PMID 41045047. |
| 15 | Martin et al. (2023). *Modelling the dynamics of senescence spread.* Aging Cell. doi:10.1111/acel.13892 | **P** | **Bibliography fully correct.** The **attributed containment mechanisms are not the paper's** — see §2.1 B4b. `degrad` and `refractor` appear **zero times** in the full text. | Martin L, Schumacher L, Chandra T. *Aging Cell.* 2023;22(8):e13892. PMID 37288475. PMC10410058. |
| 16 | Basisty et al. (2020). PLoS Biology 18(1):e3000599 | **I** | Title truncated — ends *"…for aging biomarker development."* All other fields correct. | Basisty N, Kale A, Jeon OH, … Campisi J, Schilling B. *PLoS Biol.* 2020;18(1):e3000599. PMID 31945054. |
| 17 | Acosta et al. (2013). *Paracrine senescence transmission* | **E** | **Not the paper's title**; no journal, volume, pages or DOI given. Transwell / contact-not-required claim **verified verbatim**. ⚠️ A **2026 Author Correction exists and I could not read it** — see §2.2 S0d. | Acosta JC, Banito A, Wuestefeld T, et al. *A complex secretory program orchestrated by the inflammasome controls paracrine senescence.* Nat Cell Biol. 2013;15(8):978–990. doi:10.1038/ncb2784. PMID 23770676. |
| 18 | *…SASP-dependent paracrine spreading … human brain cell types* (2026). bioRxiv, **February 10, 2026** — CCL2, CXCR7, DPP4 | **E** | Posting date is **2026-02-12** (the "02.10" is the DOI submission stamp). **Now peer-reviewed** — cite the journal version. Molecule list incomplete: the headline set is four (**MIF omitted**), and DPP4's role is specifically CXCL12 cleavage. | Russo T, Riessland M. *Aging Cell.* 2026;25(8):e70673. doi:10.1111/acel.70673. PMID 42601837. |
| 19 | *A model of the onset of the SASP…* (2017). PLOS Comp Biol. doi:10.1371/journal.pcbi.1005741 | **E** | Title paraphrased — the real one spells out *"senescence associated secretory phenotype"*. No authors or article number. | Meyer P, Maity P, Burkovski A, … Kestler HA, Scharffetter-Kochanek K. *PLoS Comput Biol.* 2017;13(12):e1005741. PMID 29206223. |
| 20 | *Dissecting the heterogeneity of senescence: primary and secondary senescent states.* PMC11689308 | **E** | **Not a paper.** PMCID resolves to a **~250-word GSA meeting abstract**, single author, no methods, no data, not peer reviewed. Title also truncated. See §2.1 B8. | Neretti N. *Innov Aging.* 2024;8(Suppl 1):351. doi:10.1093/geroni/igae098.1145. **Replace with a peer-reviewed source.** |
| 21 | Fischer et al. NCEM. Nature Biotechnology. doi:10.1038/s41587-022-01467-z | **V** | Nothing; fields missing only. Note the year is **2023**, not 2022 as the DOI suffix suggests. | Fischer DS, Schaar AC, Theis FJ. *Nat Biotechnol.* 2023;41(3):332–336. PMID 36302986. |
| 22 | Cang et al. COMMOT. Nature Methods | **I** | Nothing wrong; year, volume, pages, DOI all absent. | Cang Z, Zhao Y, Almet AA, … Atwood SX, Nie Q. *Nat Methods.* 2023;20(2):218–228. doi:10.1038/s41592-022-01728-4. PMID 36690742. |
| 23 | SCILD (2026). Communications Biology. doi:10.1038/s42003-**025**-09413-w | **V** | **The year/DOI mismatch is NOT an error** — Springer manuscript-ID artefact. 2026 is correct. Attributed diffusion/competitive-binding/decay claim matches the abstract nearly verbatim. | Yu J, Zhao J, Ren T, Sun D, Wu LY. *Commun Biol.* 2026;9(1):133. PMID 41501529. |
| 24 | SpaTalk. Nature Communications | **I** | Title omits trailing *"with SpaTalk"*; no year, authors, volume or article number. | Shao X, Li C, Yang H, … Xu X, Fan X. *Nat Commun.* 2022;13(1):4429. doi:10.1038/s41467-022-32111-8. PMID 35908020. |
| 25 | RGAST: *A relational graph attention network for multi-scale cell-cell communication inference.* bioRxiv **2024** | **E** | **The title/year pairing never existed.** v1 (2024-08-10) has a different title; the CCC title is **v2 (2025-05-26)**; the published version is a **third** title. And the published framing is representation learning benchmarked on **spatial domain identification**, with CCC one downstream task — see §2.2 S13. | Gong Y, Yuan X, Yu Z. *Empowering multifaceted analysis of spatial transcriptomics data with RGAST.* Brief Bioinform. 2026;27(3):bbag298. doi:10.1093/bib/bbag298. |
| 26 | HARMONIC. bioRxiv, January 2026 | **V** | Nothing. Attributed H&E / false-positive-reduction claim confirmed from the abstract. Not yet journal-published. | Wang X, Tao C, Jiang Y, et al. bioRxiv 2026.01.22.701166, posted 2026-01-23. doi:**10.64898**/2026.01.22.701166. |
| 27 | **CellWHISPER.** bioRxiv, January 2026 — *FPR >90% for CellChat v2, COMMOT, SpaTalk under randomization* | **I / P** | **The preprint is real** and the >90% / <5% figures are substantiated. Two problems: (a) the **bioRxiv DOI prefix for 2026 preprints is `10.64898`, not `10.1101`** — if the .bib says 10.1101 it resolves to nothing; (b) **their null is not "coordinate randomization"** — it permutes locations *within each cell type*. See §2.1 B9, which is the most consequential single finding in this audit. | Kumar A, Rivera Moctezuma F, Aggarwal B, Zhang N, Coskun AF, Sinha S. bioRxiv 2026.01.07.697982, posted 2026-01-08. doi:10.64898/2026.01.07.697982. (v2, 2026-04-03, retitled *"CellWHISPER disentangles direct cell–cell communication from structural proximity"*.) |
| 28 | **CONCISE.** bioRxiv, June 2026 | **I** | Title omits trailing *"with CONCISE"*; DOI prefix `10.64898`. **The a = 0.1 claim is verified verbatim from the PMC full text.** But the competitors are **MERINGUE, SpatialDM, Copulacci, LIANA+** — spatial co-expression / bivariate methods, **not** CellChat/COMMOT/SpaTalk. See §2.2 S16. | Zhao J, Shan X, Wang G, Chu T, Lin C, Chang RB, Zhao H. *…with CONCISE.* bioRxiv 2026.06.22.733860, posted 2026-06-28. doi:10.64898/2026.06.22.733860. PMID 42395397. PMC13320749. |
| 29 | *Benchmarking cell-type-specific spatially variable gene detection methods* (2026). Briefings in Bioinformatics 27(2):**bbag190** | **E** | **The "bbag190" suspicion is refuted — the ID is correct.** BiB article IDs run bbad=2023, bbae=2024, bbaf=2025, **bbag=2026**; vol 27 / issue 2 / March 2026 all check out. The error is in the **attributed claim**: rotation invariance is exact, but cell-type confounding is a documented weakness of **one method (Celina)**, not an authorial "central open challenge". See §2.2 S17. | Yao H, Mu S, He F, Fang Z. *Brief Bioinform.* 2026;27(2):bbag190. doi:10.1093/bib/bbag190. |
| 30 | *Identifiability limits of physics-informed inference…* (2026). **arXiv 2607.01749** | **V** | **The "anomalous" arXiv ID resolves exactly.** I fetched `arxiv.org/abs/2607.01749` directly. 2607 = July 2026, in the past; five post-dot digits is correct modern format. | Gu R, Zhang RZ, Miles CE. arXiv:2607.01749 [q-bio.QM; physics.bio-ph; stat.ML], submitted 2026-07-02. |

### 1.2 Non-reference factual claims audited alongside §31

| Claim | Where | ? | Finding |
|---|---|---|---|
| SenNet portal lists **1,753** public human+mouse datasets, **15** organs, **6** assay types as of January 2026 | §4.1, §7 Rank 2 | **P** | The **numbers are verbatim correct**, but they are **mis-sourced**: they come from the SenNet Portal preprint (Börner K, Blood PD, Silverstein JC, et al., bioRxiv doi:10.64898/2026.02.06.704469, PMC12918883), **not from the portal**, which is a JavaScript app that will not corroborate them. And **v2 of that preprint already updates the figure to 2,041 as of April 2026.** Cite the preprint, keep the "as of January 2026" qualifier. |
| `docs.sennetconsortium.org/apis` exists | §7 Rank 2 | **V** | Exists; documents Entity, Ingest, Search, UUID and Ontology/UBKG APIs. |
| ml4spatialbio deadline 2026-08-29 **AoE**; call asks for *"benchmarks, datasets, and evaluation standards specific to spatial tasks"* and *"interpretable and uncertainty-aware spatial models that biologists can trust"* | `README.md` | **V** | I fetched the CFP. Deadline is **August 29, 2026, 11:59 PM AoE** — the README's "AoE" is right (one search snippet said UTC; the CFP page says AoE, and the CFP page wins). Both quoted topic lines are verbatim from the call. 4-page limit, non-archival, double-blind. |
| RunPod RTX A5000 ~$0.27/hr, A40 ~$0.44/hr Secure Cloud | §14 | **V** | A5000 Secure $0.27 (Community $0.16); A40 Secure $0.44 (Community $0.35). |
| Container disk ~$0.10/GB-mo running, **~$0.20/GB-mo stopped** | §14, §19, §20 | **E** | **Two products conflated.** Container disk is $0.10/GB-mo **while running and is not billed when stopped — it is erased.** The $0.10-running / $0.20-stopped pair is the **volume disk**. The §14/§19 "you pay more when idle, which is backwards from EBS" argument therefore applies to volume disk, not container disk. |
| Network volume from ~$0.05/GB-mo | §14 | **P** | $0.05 is the **>1 TB** tier; under 1 TB it is **$0.07/GB-mo**. Our 300 GB volume bills at $0.07, so §27's "~$6 for 300 GB / 12 days" is understated by ~40%. |
| `cpu3g-8-32` = 8 vCPU/32 GB; `cpu5c-4-8` = 4 vCPU/8 GB; disk auto-sizes vCPU×10 GB (×15 for cpu5c) | §14 | **V** | Exact. Format is `cpu[gen][type]-[vCPU]-[RAM]`. |
| RunPod GPU pods ship **167–283 GB** system RAM on shared datacenter cards | §14 | **V** | Range is exact: 283 GB (B200), 276 GB (H200), 188 GB (RTX PRO 6000), 167 GB (RTX 6000 Ada). |
| AWS r7i.2xlarge ~$0.53/hr; r7i.4xlarge ~$1.06/hr; 300 GB gp3 ~$0.08/GB-mo | §21 | **V** | $0.5292/hr, $1.0584/hr (identical us-east-1 / us-west-2); gp3 $0.08/GB-mo with 3,000 IOPS + 125 MB/s baseline included. |
| Egress is free on RunPod | §14, §19, `README.md` | **V** | *"…billed by the second for compute and storage, with no fees for data ingress or egress."* |

---

## 2. Claims to fix, by severity

Severity definitions: **BLOCKING** = would be a factual error in the submitted
paper, or a claim whose collapse takes a headline with it. **SHOULD-FIX** = an
overreach, a stale number, or an unsourced assertion that a reviewer will
challenge. **MINOR** = internal inconsistency or imprecision.

### 2.1 BLOCKING

**B1 — "Zhao et al. (2024)" is the wrong first author for the closest-prior-work
citation.**
*Where:* `README.md` L14 and L280-ish framing; master plan §4.1, §7 Rank 3, §31 ref 1,
§32 reading-list item 1.
*Finding:* **VERIFIED ERROR.** Europe PMC record for DOI `10.1016/j.cell.2024.10.019`
gives the author string beginning `Ma S, Ji Z, Zhang B, Geng L, Cai Y, …`. There *is*
a `Zhao L`, at **position 14 of 48**. The corresponding authors are Qu J, Zhang W,
Gu Y and Liu GH.
*Why it is blocking:* this is cited in the first paragraph of the README as the paper
we are re-examining, and it is §31's flagged **"closest prior work."** Getting the
first author of your closest prior work wrong is the single most visible citation
error a reviewer in this subfield can find.
*Corrected form:* `Ma S, Ji Z, Zhang B, et al. Spatial transcriptomic landscape
unveils immunoglobin-associated senescence as a hallmark of aging. Cell.
2024;187(24):7025–7044.e34. doi:10.1016/j.cell.2024.10.019. PMID 39500323.`
Note the title really is spelled **"immunoglobin"** (not "immunoglobulin") in the
published title, while the abstract body uses "immunoglobulin". Do not "correct" it.

**B2 — Ma et al. is Stereo-seq spot data, not single-cell resolution, and the
master plan's own platform guidance says that matters.**
*Where:* master plan §7 Rank 3 ("If deposited at single-cell resolution…"), §4.1.
*Finding:* **VERIFIED.** CNGB STOMICS `STDS0000247` records the platform as
**Stereo-seq on DNBSEQ-T1**, at spot/bin level, nine tissues, and the paper reports
1,535,191 spots at ~1,450 genes/spot. Master plan §3 explicitly warns that spot-based
data blurs distance and is unsuitable for kernel estimation.
*Consequence, and it is in our favour:* the closest prior distance-gradient result in
the field was measured on **binned spots, not segmented cells**. That strengthens our
framing and should be said explicitly rather than left as an implication. It also
means §7 Rank 3 ("reproducing their distance-ranked gradient … is your strongest
Figure 2") was never executable as written.
*Also:* §7 Rank 3 says "Chinese consortium data is sometimes in GSA rather than GEO."
It is in **CNGB STOMICS DB (`STDS0000247`)**, neither GEO nor GSA. Fix the pointer.

**B3 — "multiple life stages" is wrong; Ma et al. profiled two age groups.**
*Where:* master plan §7 Rank 3 "Multiple organs, multiple life stages"; §4.1 "across
life stages".
*Finding:* **VERIFIED ERROR.** The STOMICS record lists exactly two stages,
**"Young, Old"**, male mice. "Across life stages" implies a time course we do not
have evidence for. Reword to "young versus old male mice across nine tissues."

**B4 — the README asserts ">90% false-positive rates" as a published fact, while
CS_PHASE4 §5.7 explicitly says our version of that number is not a false-positive
rate.**
*Where:* `README.md` L14–18 ("despite contemporaneous benchmarks reporting >90%
false-positive rates for the class of methods used to make such claims");
`CS_PHASE4.md` headline item 1 ("CellWHISPER's >90 % figure replicates for CellChat
v2\* (94%)"); versus `CS_PHASE4.md` §5 limitation 7: *"Significance survival is not a
false-positive rate. It is the fraction of real-data calls reproduced on shuffled
data. It equals an FPR only if one accepts the shuffled tissue as a true null, which
§2.4 and CS_PHASE3 §5 both argue is exactly the assumption in question."*
*Finding:* this is an **internal contradiction inside one report**, and the README
propagates the weaker-provenance version. The whole paper argues that treating a
coordinate-shuffled tissue as a true null is the mistake — so we cannot simultaneously
report survival-on-shuffled-data as an FPR.
*Suggested rewording (README):* "…despite a contemporaneous benchmark reporting that
leading spatial CCC tools return comparable interaction counts on real and
coordinate-randomised input."
*Suggested rewording (CS_PHASE4 headline 1):* "CellWHISPER's qualitative finding
reproduces: 94% (CellChat v2\*), 79% (COMMOT), 79% (SpaTalk\*) of real-data-significant
calls are still significant on fully permuted coordinates. We report this as
*significance survival*, not as a false-positive rate — see §5.7."
*Dependency:* whether CellWHISPER states a ">90%" number at all for competing methods
is a §1 item; see §3 for what I could and could not confirm.

**B4b — the containment paradox is attributed to Martin et al. with four mechanisms
that are not in the paper.**
*Where:* master plan §3 "Why decay must exist": *"Martin et al. showed that plausible
parameter values for ligand diffusion and binding do not by themselves produce local
containment, so something else limits spread: **rapid degradation, a response
threshold, immune clearance of senders, or receptor-level refractoriness in
receivers.**"* Also §31 ref 15 (flagged "**the containment paradox**"), §32 reading-list
item 6, and §6.4, which frames the λ_proximal / λ_downstream comparison as speaking
"directly to Martin's containment question."
*Finding:* **VERIFIED ERROR, from the full text** (Europe PMC `PMC10410058`, the
open-access XML of Martin L, Schumacher L, Chandra T. *Aging Cell* 2023;22(8):e13892).
1. The **paradox itself is real and correctly attributed.** Abstract: *"In the absence
   of the immune system, senescence could theoretically spread infinitely from one cell
   to another, but this contradicts experimental evidence."* Significance statement:
   *"…current understanding fails to explain how senescence can spread in a controlled
   and local way."*
2. **The four resolving mechanisms are not the paper's.** In the full text the string
   `degrad` occurs **zero times** and `refractor` occurs **zero times**. Immune
   clearance appears only as the mechanism the paper deliberately brackets out:
   *"Currently, there are no hypothesised mechanisms to contain the spread of senescence
   in the absence of the immune system."* A ligand-binding threshold N_D exists as a
   model parameter but is never offered as the containment mechanism.
3. **The paper's actual three answers are different ones:** juxtacrine secondary
   senescent cells acting as a firebreak because they do not produce SASP (*"analogous
   to removing trees to prevent the spread of a forest fire"*); time-delayed dynamic
   induction (*"dynamic, time-dependent paracrine signalling prevents the uncontrolled
   spread of senescence"*); and secondary senescent cells secreting **fewer** SASP
   molecules than primary ones.
4. **The strength of the claim is also overstated.** The paper says explicitly:
   *"Using our minimal model we can find parameter regimes that lead to senescence
   spread from a single cell. However, we can not determine from this model whether the
   spread of senescence is controlled (finite) or uncontrolled."* So "plausible
   parameter values … do **not** produce local containment" asserts more than the
   minimal model establishes.
*Why blocking for biology specifically:* this is the paper the Discussion is supposed
to speak to, and §6.4's whole rationale (*"If λ_proximal > λ_downstream, the response
is thresholded rather than graded, which speaks directly to Martin's containment
question"*) is built on the **threshold** mechanism, which is the one attribution that
is closest to being defensible but is still not the paper's proposed resolution.
Fortunately `CS_PHASE5.md` §6.4 concluded the comparison is not estimable and reported
it as a bound, so **no result depends on this** — but the Introduction/Discussion text
does, and it must be rewritten.
*Suggested rewording:* "Martin et al. showed with a mathematical model that current
mechanistic understanding of SASP diffusion and binding does not explain how senescence
spread stays local, and proposed that containment arises from the properties of
secondary senescent cells themselves — that they are poor SASP producers, and that
induction is time-delayed rather than instantaneous."

**B5 — README's "Known limitations" repeats a claim that BIO_PHASE3 and CS_PHASE3
both formally withdrew.**
*Where:* `README.md` L355-ish: *"the sham arm fails Test 3 (hepatocytes 0.48%, below
the 1% floor), so estimation is SBR-only."*
*Finding:* **VERIFIED ERROR against our own results.** `BIO_PHASE3.md` §5 corrects
this: sham hepatocyte `Cdkn1a`⁺ is 7.20 / 2.31 / 4.92 / **0.48** / 8.94 across the five
sham sections — **4 of 5 pass**, and 7250 (the only sham section Phase 2 had) is the
outlier. `CS_PHASE3.md` §1(a) goes further: the admissible set is **six sections, two
SBR and four sham**, because four of six SBR sections fail the *20% ceiling* and 7250
fails the floor. So the headline analysis is majority-sham, and the README says the
opposite.
*Corrected form:* "Admissibility follows §8 Test 3 prevalence, not surgical arm: four
of six SBR sections exceed the 20% ceiling and one sham section falls below the 1%
floor, leaving six admissible sections from six animals (2 SBR, 4 sham), with arm as a
contrast."
*Why blocking:* a reviewer who reads the README and then the Methods will find the two
disagree about which arm the paper is about.

**B6 — "sender rate tracks sequencing depth at ρ = 0.94" is listed as a
contribution, and it is an n = 6 within-arm correlation that BIO_PHASE3 explicitly
declined to call technical.**
*Where:* `README.md` contribution 4 ("Distance-to-nearest-sender is a
sender-calling-rate readout to r² = 0.98, and sender rate tracks sequencing depth at
ρ = 0.94. This applies to every method in this literature that regresses response on
distance-to-nearest-X.").
*Finding:* three separate problems.
1. `BIO_PHASE3.md` §5 states ρ = +0.943, p = 0.005, **n = 6, within the SBR arm only**;
   within sham ρ = +0.50, p = 0.39, and **pooled ρ = +0.16**. The README quotes the
   most favourable of three numbers with no n and no arm.
2. BIO_PHASE3 §5 adds an explicit hedge the README drops: *"I am **not** claiming this
   is purely technical — deeper sections may also be less injured, with larger, more
   transcript-rich hepatocytes, so depth and phenotype are entangled."*
3. **"Sequencing depth" is the wrong term for Xenium.** This is imaging-based in-situ
   hybridisation; there is no sequencing. The measured quantity is *median transcripts
   detected per cell*.
*Suggested rewording:* "Within the SBR arm, section-level `Cdkn1a`⁺ prevalence tracks
section-level median transcripts detected per cell (Spearman ρ = +0.94, p = 0.005,
n = 6 sections); the association is weaker in sham (ρ = +0.50, n = 5) and absent when
arms are pooled (ρ = +0.16). Depth and injury severity are entangled and we do not
separate them; the practical consequence is that animal-level burden must not be
compared across sections without a depth covariate."
Also drop or heavily qualify *"This applies to every method in this literature that
regresses response on distance-to-nearest-X."* The Poisson `d ≈ 0.4697 ρ^(−1/2)`
identity does generalise — that part is mathematics, not an empirical claim — but the
depth–burden coupling is one dataset, one panel, one disease model.
*Suggested split:* keep the r² = 0.98 sender-rate result as the general claim (it is a
geometric identity confirmed across 77 section × caller combinations); demote the
ρ = 0.94 depth coupling to a dataset-specific caution.

**B7 — the DeepScence caller results are reported without the caveat their own
report insists on, and rest on 2 of 11 sections.**
*Where:* `README.md` L296-ish: *"Senescence calling is a choice, not a fact.
DeepScence, SenePy and the curated Tier A score agree at or below chance…"*
*Finding:* **VERIFIED ERROR of provenance.** `BIO_PHASE3.md` §4.3 says, in terms:
*"we run DeepScence with `denoise=False`, a forced deviation (its DCA dependency needs
an obsolete TensorFlow stack), and DCA denoising is precisely the step that would
normalise depth. So §4.2–4.3 characterise DeepScence as we could run it on this panel,
not DeepScence as published. It should be reported in exactly those words."* The
README reports it in none of those words. `BIO_PHASE3.md` "Open / not done" also
records that **DeepScence scores exist for the two 26 wk sections only**.
> **⚠ SUPERSEDED 2026-08-27 (Phase 8).** The suggested rewording below quotes
> **0.93–1.22×** and **1.51–2.85×**. Both are properties of the **two-section,
> pre-C6 Tier A** base and no longer hold. At 11-section coverage on the frozen
> strict-33 Tier A the pooled agreement is **1.212×**, and the circular
> DeepScence vs `Cdkn1a`⁺ pair is **weaker**, not stronger (median 1.071 over
> eleven sections, vs the 1.51–2.85 quoted here). The README has already been
> corrected; this block is left **unedited as a record of what was proposed at
> the time**. Do not apply it. Current numbers:
> `results/phase3/caller_coverage_gate{,_headline}.csv`, which carry six
> explicitly-labelled bases. See `reports/CS_PHASE8_CALLERS.md` and
> `reports/SUBMISSION_PATCH_2026-08-29.md`.

*Additionally:* the README quotes the **Phase 2 sham-only** Jaccard ratios
(0.60× / 0.88× / 1.66×) which were superseded by the Phase 3 cell-type- **and**
depth-matched ratios (0.93–1.22× of chance for four of six pairs, but 1.51–2.85× for
DeepScence vs `Cdkn1a`⁺, which is above chance and expected because `CDKN1A` is
DeepScence's own sign anchor). Quoting the stale numbers understates the nuance and
overstates the uniformity.
*Suggested rewording:* "Four sender definitions — DeepScence (run with `denoise=False`,
a forced deviation from the published configuration, on two sections), SenePy,
`Cdkn1a`⁺ and a curated arrest-and-damage score — overlap at 0.93–1.22× of chance for
four of six pairs after conditioning on cell type and transcript-depth decile. The one
consistently above-chance pair, DeepScence vs `Cdkn1a`⁺, is expected: `CDKN1A` is
DeepScence's sign anchor."

**B9 — we labelled the wrong null as CellWHISPER's, and it is the null Figure 4 is
built around. This is the most consequential finding in the audit.**
*Where:* `CS_PHASE4.md` §2.4, the null-ladder table, final row: *"**N0_perm** | permute
all cell coordinates among cells | … | **this is CellWHISPER's randomisation**."* Also
§2.5 (*"This is the CellWHISPER quantity"*), §3.3 (*"Do we reproduce CellWHISPER's
>90 %?"*), the Phase 4 headline, and `README.md` (*"reproducing CellWHISPER's criterion
directly on senescence-relevant pairs"*).
*Finding:* **ERROR of attribution.** CellWHISPER's randomised control **permutes cell
locations within each cell type**, in their words *"preserving cell-type-specific
spatial organization and ligand-receptor (LR) expression while destroying spatial
proximity between ligand- and receptor-expressing cells."* Our **N0_perm permutes all
coordinates among all cells**, which additionally destroys the cell-type spatial
architecture. **The two nulls are not the same, and ours is strictly the more
destructive of the two.**
*What this does and does not damage:*
- It does **not** damage the conclusion. Because our null destroys strictly more, our
  survival rates are if anything *conservative* relative to theirs — we are showing the
  statistics survive an even harsher shuffle. The direction of the error is in our
  favour.
- It **does** mean the sentence "CellWHISPER's >90% figure replicates" is not accurate
  as a replication claim. **We did not run their null.** Neither N0_perm (too
  destructive) nor N3_type (a per-cell-type rigid *shift*, which preserves each type's
  internal geometry exactly, rather than a within-type *permutation*, which destroys it)
  is their design.
- The **numerically closest** thing we have is N3_type: COMMOT 0.769, CellChat v2\*
  0.913, SpaTalk\* 0.784, NCEM linear\* 0.080. Note **CellChat still exceeds 90% there**,
  so the headline comparison to their >90% survives under the closer null too.
*Two options for the CS lead, and the first is cheap:*
1. **Run the within-cell-type coordinate permutation.** It is a small change to
   `phase4_run.py`'s null ladder — permute coordinates among cells sharing a
   `cell_type_merged` label — and it converts an attribution error into an actual
   replication. Given the machinery already exists, this looks like well under a day.
   It also fits the paper's own argument: a within-type permutation is precisely a
   *confounder-preserving* null of the kind CS_PHASE4 §4.5 says these tools need.
2. **If not run, relabel honestly.** Change §2.4's row to *"N0_perm | permute all cell
   coordinates among cells | … | all association between position and cell identity.
   This is strictly more destructive than CellWHISPER's control, which permutes
   locations within cell type."* And change the headline from "CellWHISPER's >90 %
   figure replicates" to *"we reproduce CellWHISPER's qualitative finding under a
   strictly harsher shuffle."*
*Related wording fix in the plan:* §4.3 and §31 ref 27 both say CellWHISPER benchmarked
the tools on *"coordinate-randomized data"* / *"under randomization."* Both should say
**within-cell-type location permutation**. As agent-verified: describing their control
as crude coordinate randomisation understates its design and invites the reviewer
objection that the benchmark was too blunt — which is the opposite of the point we want
to make with it.

**B8 — §31 ref 20 is a conference abstract, not a paper, and it is the only citation
for the primary-versus-secondary senescence distinction.**
*Where:* §31 ref 20, "Dissecting the heterogeneity of senescence: primary and secondary
senescent states. PMC11689308"; the claim it supports is master plan §3 *"Secondary
senescent cells are transcriptomically distinct from primary ones, particularly in
their own SASP output."*
*Finding:* **VERIFIED ERROR of document type.** PMC11689308 resolves to
**Neretti N.** *Innovation in Aging* 2024;8(Suppl 1):351, doi:10.1093/geroni/igae098.1145
— a ~250-word **GSA Annual Scientific Meeting abstract** from "SESSION 3020
(BIOLOGICAL SCIENCES INVITED SYMPOSIUM)", single author, no methods, no data, not peer
reviewed. The title is also truncated; the real one continues *"…Profiling Primary and
Secondary Senescent States Using Single-Cell and Spatial Transcriptomics."*
*Why blocking:* the primary/secondary distinction is not a decorative claim for us — it
underwrites Tier B module B7 (`secondary_senescence`), which is one of seven readouts,
and it is the biological premise for §6.4. Citing a meeting abstract for it will not
survive review.
*Fix:* replace with a peer-reviewed source. Martin et al. 2023 (ref 15) itself
distinguishes primary from secondary and states the reduced-SASP property, so it can
carry part of this; Acosta et al. 2013 (ref 17) carries the paracrine-induction half.

### 2.2 SHOULD-FIX

**S0a — the human-endometrium calibration is a 2025 paper in a journal we do not name,
and its senescence call is not transcriptomic.**
*Where:* master plan §3 "Length scales to expect"; §31 ref 7.
*Finding:* **the five numbers are VERIFIED exactly**, quoted verbatim from the paper:
*"Macrophages (CD68+) and monocytes (CD14+) were positioned closest to senescent cells,
with mean distances of 45 ± 20 μm and 45 ± 25 μm, respectively. NK cells (CD56+) and
total T cells (CD3+) were located at intermediate distances (53 ± 23 μm and 62 ± 29 μm).
In contrast, T-helper cells (CD4+) were positioned significantly farther away
(102 ± 42 μm; p < 0.05), while B cells (CD79α+) displayed the greatest separation from
senescent cells (211 ± 66 μm; p < 0.01)."* But three bibliographic facts are wrong or
missing: the year is **2025, not 2026**; the journal is **Diagnostics (Basel)**, an MDPI
title, which §31 omits entirely; and the title is truncated — it ends *"…During the
Implantation Window."*
*Two methodological caveats that matter more than the bibliography,* because this is our
only external length-scale anchor: (i) senescence was called by **p16
immunohistochemistry on adjacent serial sections** digitally aligned in HALO, not by a
transcriptomic senescence score, so it is not the same measurement as ours; (ii) the
cohort is **mid-luteal endometrium from 68 IVF patients during the implantation
window** — a physiologically senescence-rich, actively decidualising tissue. These
distances are not a general tissue constant and should not be presented as one.
*Corrected form:* Parvanov D, Ganeva R, Ruseva M, et al. Spatial Distribution of
Senescent Cells and Their Proximity to Immune Subsets in the Human Endometrium During
the Implantation Window. *Diagnostics (Basel).* 2025;15(21):2679.
doi:10.3390/diagnostics15212679. PMID 41225972.

**S0b — the Karpova framing omits the half of their comparison that *was*
significant.**
*Where:* master plan §7 Rank 1: *"They compared clustered versus isolated senescent
hepatocytes, examined the surrounding microenvironment, found trends toward more hepatic
stellate cells and inflammatory macrophages near clustered senescent cells, and
explicitly reported this as falling short of significance… They had the right data, saw
the trend, and could not call it."*
*Finding:* **PARTIAL — half right.** The **niche-composition** comparison did fall short:
*"Although falling short of significance, we observed trends of increased content of
HSCs and inflammatory macrophages and decreased content of non-inflammatory
macrophages."* But a separate **gene-expression** comparison of clustered versus
isolated senescent hepatocytes **did** reach significance for **CXCL2, IGFBP2 and SDC4**
(Wilcoxon, FDR < 0.05). Writing "they could not call it" without that qualification
misrepresents their result, and a reviewer who knows the paper will notice.
*Suggested rewording:* "Their niche-composition comparison between clustered and
isolated senescent hepatocytes fell short of significance, though a parallel
gene-expression comparison did reach it for three genes. Neither is a distance-resolved
estimate with uncertainty, which is the gap we address."

**S0c — the periportal localisation we treat as a standing anatomical confound is
age-conditional in the source.**
*Where:* master plan §11 and §7 Rank 1 Caution: *"They report CDKN1A+ hepatocytes
localize periportally, within 100–150 μm of the portal triad in young individuals."*
*Finding:* **VERIFIED, and §11's later use drops the qualifier.** The paper: *"In young
individuals, CDKN1A+ hepatocytes were preferentially located near the portal triad
(<100–150 μm)… However, when senescence exceeded ∼15% of total hepatocytes in older
donors, this periportal enrichment diminished."* §7 keeps "in young individuals"; §11's
argument chain ("1. Senders cluster periportally.") drops it and treats periportal
clustering as a general property.
*Why this matters to us specifically:* our sections run at 2.3–45% `Cdkn1a`⁺
hepatocytes (`BIO_PHASE3.md` §5) — i.e. mostly **above** the ~15% threshold at which the
source says the enrichment disappears — and the finding was human, not mouse
(`BIO_PHASE1.md` §2.1 already flags the species transfer). This is a second, independent
reason §11's zonation prediction failed on our data, and it is a better one than "we
measured it and it wasn't there." Worth saying.

**S0d — Acosta et al. 2013 carries a 2026 Author Correction whose content I could
not read.**
*Where:* master plan §3 (transwell / contact-not-required claim) and §31 ref 17.
*Finding:* the **attributed claim is VERIFIED** from the full text
(`PMC3732483`): *"To test whether soluble factors mediate paracrine senescence, we used
transwell inserts that ensure physical separation of the cells."* The citation itself is
badly incomplete — "Paracrine senescence transmission" is not the paper's title and no
journal, volume or DOI is given. Correct form: Acosta JC, Banito A, Wuestefeld T, et al.
A complex secretory program orchestrated by the inflammasome controls paracrine
senescence. *Nat Cell Biol.* 2013;15(8):978–990. doi:10.1038/ncb2784. PMID 23770676.
**UNVERIFIED and flagged for a human:** *Nature Cell Biology* published an **Author
Correction on 2026-04-28** (doi:10.1038/s41556-026-01959-z, PMID 42050148). I confirmed
it exists via Crossref and Europe PMC but **could not read what it corrects** —
nature.com redirects anonymous fetches to `idp.nature.com`, and the Europe PMC record
carries no abstract and no PMCID. Because this is the foundational citation for
"physical contact is not required", **someone with library access should read the
correction before submission.** I am not asserting anything is wrong with it.

**S1 — "The statistic is not spatial" is stronger than CS_PHASE4's own section
heading, which says the opposite about half of COMMOT.**
*Where:* `README.md` L96 and L205 (twice, both in bold).
*Finding:* `CS_PHASE4.md` §4.1 is titled **"COMMOT: the optimal transport is spatial,
the cluster summary is not,"** and the body says *"COMMOT is not insensitive to
geometry. Permuting coordinates rebuilds its cell-to-cell communication network almost
from scratch."* The README compresses that into a universal negative.
*Suggested rewording:* "The cluster-level summary these methods test is close to
geometry-free, even where the underlying cell-to-cell computation is not."

**S2 — mass conservation is stated without its one large exception.**
*Where:* `README.md` L207-ish: *"conserves transported ligand mass to seven significant
figures (e.g. 863.815754 → 863.815745)."*
*Finding:* `CS_PHASE4.md` §4.1's table gives N0÷real transported mass of 1.000000,
1.000054, 1.000000 — and **1.127084 for `Ccl2→Ccr2`**, a 13% change, with the report
giving the reason (Ccr2 detected in 1–3% of cells, so some Ccl2 cannot find a receiver
within 100 µm on real coordinates). The README states the conservation as unqualified.
Note also the README says "**six** figures" at L94 and "**seven** significant figures"
at L207 for the same quantity.
*Suggested rewording:* "conserves transported ligand mass to six or more significant
figures for three of the four pairs (the exception, `Ccl2→Ccr2` at 1.13, is itself
diagnostic: `Ccr2` is too rare for all `Ccl2` to find a receiver on real coordinates)."

**S3 — "Zonation is not the liver confound" cites only the most favourable
stratum.**
*Where:* `README.md` L141–143: *"The zonation covariate alone removes ~0% in
hepatocytes (SF 1.043), and the kernel does not vanish within zones because zonation
was never driving it."*
*Finding:* `CS_PHASE3.md` §0 table gives **zonation covariate alone, median SF 0.843
[0.510, 0.992] over all fits**, and §3.2 gives 1.043 in hepatocytes *and*
**0.244 in `Biliary/ductular`**, where the report notes the zonation score is acting as
a proxy for ductular-reaction geography. So zonation removes ~16% overall and ~76% in
one receiver type. The headline is defensible *for hepatocytes*, which is the cell type
§11's argument is actually about — but the README does not say so.
*Suggested rewording:* "Zonation is not the dominant confound where §11 predicted it
would be: entered as a covariate it removes essentially nothing in hepatocytes
(SF 1.043) and 16% overall (SF 0.843), though it removes 76% in `Biliary/ductular`,
where the score is a proxy for ductular-reaction geography rather than for zonation."
*Causal-language note:* "**because zonation was never driving it**" is a causal
statement resting on a covariate-adjustment result. See §4, flag C2.

**S4 — the master plan asserts a zonation loading for B6 that our own data
contradicts.**
*Where:* master plan §11: *"Several Tier B modules plausibly vary with zonation.
Oxidative stress (B6) certainly does; pericentral hepatocytes are the site of
CYP450-driven oxidative metabolism."*
*Finding:* **contradicted by our own measurement.** `BIO_PHASE1.md` §2.5 measured
B6 oxidative_stress at **r = −0.025** with the zonation axis within hepatocytes
(n = 142,143), the *smallest* absolute loading of any of the seven Tier B modules. The
mechanistic reasoning about CYP450 is sound; the transcriptional module we built from
it does not carry the loading. Either drop the "certainly does", or state that the
prediction was tested and not borne out — the latter is more useful, since §11 is the
plan's showcase for why the null battery exists and this is a second instance of the
plan's liver prediction failing.

**S5 — PHASE0_DATA_AND_ENV.md still says the CCL2/CXCR7/DPP4 axis is "fully
covered"; BIO_PHASE3 refuted it and the refutation is not stamped on the earlier
file.**
*Where:* `PHASE0_DATA_AND_ENV.md` §8.1: *"Tier C receptor coverage includes `Ccr2`,
`Ackr3`, `Cxcr4` and `Dpp4`, so the CCL2 / CXCR7 / DPP4 axis is fully covered."*
*Finding:* `BIO_PHASE3.md` §3.1 issues an explicit correction: **`Ackr3` (CXCR7) is
detected in ≤1.0% of cells in every cell type of every section**, so CCL2→ACKR3 and
CXCL12→ACKR3 are not testable here. On-panel ≠ detected.
*Fix:* add a correction stamp to PHASE0 §8.1 pointing at BIO_PHASE3 §3.1, the same way
BIO_PHASE2 and BIO_PHASE3 stamp their own corrections. The phase reports are the
provenance record; leaving a refuted claim unmarked in the earliest one is the
mechanism by which the GSE310392 error propagated for weeks.

**S6 — "~1.3M cells" does not reconcile with any cell count in any report.**
*Where:* `README.md` L272 (*"a negative result with a quantified bound on 1.3M cells
across 11 sections from 11 animals"*) and L286 (*"~1.3M cells"*).
*Finding:* summing `BIO_PHASE3.md` §1: QC-pass cells across the 11 liver sections is
**1,826,893**; analysable (excluding `Low_quality`/`Unknown`) is **1,635,937**. The six
*admissible* sections that the bound is actually computed on total **936,125**
analysable cells. **None of these is 1.3M**, and I could not locate a derivation for
1.3M anywhere in `reports/`, `results/` filenames, or the README itself.
*Fix:* replace with whichever number is meant, and say which denominator it is. My
recommendation, because it is the honest one and the one the bound belongs to:
"…on 936k analysable cells in six admissible sections from six animals, drawn from
1.83M cells across 11 sections from 11 animals."
*Flagged as SHOULD-FIX rather than MINOR because it appears twice, once inside a
contributions list.*

**S7 — README claims a segmentation result from one section as a property of the
dataset.**
*Where:* `README.md` L~330: *"Segmentation is not the problem here … 88.27% assigned
to a cell … so the confounding documented here is not attributable to transcript
bleed-through."*
*Finding:* the measurement is on **one section (7259)**, as the README itself says. The
conclusion ("the confounding documented here is not attributable to bleed-through") is
drawn for all six admissible sections. Also note assignment rate varies with cell
density, and BIO_PHASE1/CS_PHASE2 record a 45% spread in packing across sections, so
one section is a weak basis for generalising.
*Suggested rewording:* "…measured on one admissible section (7259). Transcript
bleed-through is therefore not a plausible driver of the confounding *in that section*;
we did not measure the other five."

**S8 — an unsourced claim about Xenium panel design, stated as a hard rule.**
*Where:* `BIO_PHASE1.md` §3.3 and `PHASE0_DATA_AND_ENV.md` §8.2: *"Xenium omits
high-expressors to avoid optical crowding."*
*Finding:* **PARTIAL.** 10x Genomics' Xenium Panel Designer documentation *recommends
against* including genes that are moderately-to-highly expressed across many cell
types, on exactly the optical-crowding grounds we give, and names `ACTB`, `B2M`,
`EEF1A1`, `GAPDH`, `NEAT1` as examples. But that is a design recommendation, not a
guarantee that any given panel omits them. Our *observation* — that `Actb`, `Gapdh`,
`Rpl13a`, `Rps18`, `Ppia` and `Alb` are absent from this specific panel — is our own
measurement and stands.
*Suggested rewording:* "absent from this panel, consistent with 10x's own panel-design
guidance to avoid broadly high-expressed genes because of optical crowding
(Xenium Panel Designer documentation)."

**S9 — an unsourced biological claim doing real work in the sender-calling
recommendation.**
*Where:* `BIO_PHASE2.md` §4.4: *"`Cdkn1a` is induced in cycling cells, so it is not
senescence-specific. Recommend excluding `Proliferating` from sender calls."*
*Finding:* the claim is **well supported in the literature** (p21/CDKN1A has
CDK-independent roles, is expressed in proliferating cells, functions at the G2/M
transition, and the magnitude/timing of a p21 pulse — not its presence — distinguishes
proliferative from senescent fates), but we assert it with no citation while using it
to justify an exclusion that changes the sender set. Add a citation.

**S10 — "functionally interchangeable" overstates the Cxcl1/Cxcl2 relationship.**
*Where:* `BIO_PHASE3.md` §3.3: *"`Cxcl1` vs `Cxcl2`, two genes encoding functionally
interchangeable analogues of the same absent human ligand (CXCL8)."*
*Finding:* the underlying species point is **correct and well handled** (mouse has no
CXCL8 orthologue; `Cxcl1`/`Cxcl2` are the conventional functional analogues, both
signalling through Cxcr2 — and CS_PHASE4 §158 and BIO_PHASE3 §3.4 both state it
correctly). But "functionally interchangeable" is stronger than the analogy supports:
they are non-redundant paralogues with distinct induction kinetics.
*Suggested rewording:* "`Cxcl1` and `Cxcl2`, the two mouse chemokines conventionally
used as functional analogues of human CXCL8, which mouse lacks."
*This is the strongest argument in §3.3 and it does not need the overstatement.*

**S11 — a human gene symbol used for a mouse result in the paragraph explicitly
written "for the paper".**
*Where:* `BIO_PHASE3.md` §4.5, the block-quoted paragraph intended for the manuscript:
*"…and `CDKN1A` positivity…"* and *"…SenePy and `CDKN1A`⁺ are enriched 2.3–2.5×…"*.
*Finding:* the data are mouse. The symbol is `Cdkn1a`. The same paragraph also says
"sequencing depth" for an imaging assay (see B6.3).
*Fix:* `Cdkn1a` throughout; "transcript detection depth" or "transcripts detected per
cell" instead of "sequencing depth". This one matters more than usual because the
paragraph is pre-drafted manuscript text.

**S12 — human-derived length-scale calibration applied to mouse liver without
comment (master plan only).**
*Where:* master plan §3 "Length scales to expect", the 45 ± 20 / 45 ± 25 / 53 ± 23 /
102 ± 42 / 211 ± 66 µm figures from the 2026 human endometrium study, used as the
sanity check for a fitted λ ("If your fitted λ comes out at 500 µm, suspect the model,
not the biology").
*Finding:* **Good news — no report or README result depends on this.** I grepped all
nine phase reports and the README; the endometrium calibration is never invoked. The
exposure is confined to the plan.
*If it is carried into the manuscript*, it needs two qualifications: it is **human**
endometrium, and mouse hepatocytes are smaller than human hepatocytes, so a
cell-diameter-anchored length scale does not transfer unchanged. Our own resolution
floor (median NN 6.7–9.8 µm) is the mouse-specific anchor and should be used instead.
*Also:* master plan §3 "SASP" lists **IL-8/CXCL8** among "the recurring names" — a
human-only ligand in a section that now introduces a mouse study. Add the mouse
caveat there, since BIO_PHASE3 §3.4 already established it downstream.

**S13 — RGAST is not primarily a cell-cell-communication method, and the title/year
pairing in §31 does not exist.**
*Where:* §31 ref 25; master plan §4.2 lists RGAST among spatial CCC methods.
*Finding:* **VERIFIED ERROR.** Three titles exist for one work. bioRxiv **v1**
(2024-08-10) is *"RGAST: Relational Graph Attention Network for Spatial Transcriptome
Analysis."* The CCC-inference subtitle we quote belongs to **v2, posted 2025-05-26**.
The peer-reviewed version is a **third** title: Gong Y, Yuan X, Yu Z. *Empowering
multifaceted analysis of spatial transcriptomics data with RGAST.* Brief Bioinform.
2026;27(3):bbag298. doi:10.1093/bib/bbag298. So "bioRxiv 2024" plus the CCC title is a
pairing that never existed. The published framing is a relational graph-attention
**autoencoder for representation learning**, benchmarked headline-first on **spatial
domain identification**; CCC is one downstream application among clustering, SVG
detection, trajectory inference and 3D reconstruction. Listing it flatly as a CCC method
overstates its focus.

**S14 — the 2026 brain paracrine-senescence preprint is superseded and its molecule
list is incomplete.**
*Where:* §31 ref 18; the CCL2/CXCR7/DPP4 axis that `PHASE0_DATA_AND_ENV.md` §8.1 and
`BIO_PHASE3.md` §3.1 both build on.
*Finding:* the bioRxiv posting date is **2026-02-12**, not February 10 (the "2026.02.10"
in the DOI is the submission stamp, not the posting date), and it has since been
**peer-reviewed**: Russo T, Riessland M. *Characterizing the SASP-Dependent Paracrine
Spreading of Senescence Between Human Brain Cell Types.* Aging Cell. 2026;25(8):e70673.
doi:10.1111/acel.70673. PMID 42601837. Cite the journal version. All three attributed
molecules verify from the full text, but the paper's headline set is **four**: **MIF is
omitted from our citation**, and DPP4's role is specifically **cleavage/inactivation of
CXCL12**, which a bare gene list loses.
*Relevance to us:* `BIO_PHASE3.md` §3.1 correctly reports that we cannot test the ACKR3
arm (`Ackr3` ≤ 1.0% detection everywhere). That statement should now name the journal
version, and should note the axis is a **human brain** result being asked about in
**mouse liver**.

**S16 — CONCISE's benchmark is against co-expression methods, not against the tools
we test, and the plan's sentence blurs that.**
*Where:* master plan §4.3: *"CONCISE (bioRxiv June 2026) showed that introducing even
weak spatial autocorrelation (a = 0.1) into one gene inflated type I error for **every
competing spatial ligand-receptor method** tested."*
*Finding:* **the a = 0.1 claim is VERIFIED verbatim** from the PMC full text
(PMC13320749), Results, "CONCISE achieves superior false positive rate control in real
data permutation studies," Scenario 2: *"Even when one gene retained its original
counts, and only weak spatial autocorrelation (a_spatial = 0.1) was introduced to the
other, all competing methods except CONCISE were substantially confounded."* Every
element checks — weak, a = 0.1, one gene, all competitors.
**But the competitors are MERINGUE, SpatialDM, Copulacci and LIANA+** — spatial
co-expression and bivariate-statistic methods. Only two of the four are
ligand-receptor tools, and **none of them is CellChat, COMMOT or SpaTalk**. Calling them
"every competing spatial ligand-receptor method" invites the reader to think CONCISE
independently indicted the same three tools CellWHISPER did. It did not.
*Suggested rewording:* "CONCISE showed that introducing even weak spatial
autocorrelation (a_spatial = 0.1) into one of two genes inflated type I error for every
competing method it tested (MERINGUE, SpatialDM, Copulacci, LIANA+)."
*Note the citation gain:* naming the four makes the claim checkable and costs half a
line. Do it.

**S17 — "cell-type confounding a central open challenge" is our framing, not the
benchmark's.**
*Where:* master plan §4.3, third bullet.
*Finding:* **rotation invariance is exact** — abstract point (iv): *"rotation invariance
still warrants further investigation"*, and the body reports that no method achieved
rotational invariance, with CTSV worst and STANCE failing in practice despite a
theoretical guarantee. **Cell-type confounding is overstated.** The paper reports it as
a specific weakness of the *best-performing* method: abstract point (v), *"Celina
appears to have a relative superiority in most metrics, though it tends to generate
spurious signals affected by nontargeted cell types"*, with the body giving the
mechanism (*"residual spatial signals may be incorrectly attributed to the target cell
type"* when colocalised cell types are not captured by fixed effects). That is a real,
citable confounding result — it is simply not an authorial declaration that cell-type
confounding is the field's central open problem.
*Suggested rewording:* "A March 2026 *Briefings in Bioinformatics* benchmark of six
cell-type-specific spatially variable gene methods across 46 real datasets found that
no method achieved rotation invariance, and that the best-performing method (Celina)
generates spurious signals attributable to non-target cell types."
*Minor:* the article page shows a publication date of 27 April 2026 against a March 2026
issue cover date. This is routine OUP behaviour, but if you want to be bulletproof,
cite "27(2):bbag190" without a month.

**S18 — the SenNet dataset count is correct but attributed to the wrong source, and
it has already moved.**
*Where:* master plan §4.1 and §7 Rank 2: *"The SenNet portal listed 1,753 public human
and mouse datasets across 15 organs and 6 assay types as of January 2026."*
*Finding:* **the numbers are verbatim correct** — but they come from the **SenNet Portal
preprint** (Börner K, Blood PD, Silverstein JC, et al., *SenNet Portal: Build,
Optimization and Usage*, bioRxiv doi:10.64898/2026.02.06.704469, PMC12918883), which
states: *"As of January 2026, the portal hosts 1,753 publicly available human and mouse
datasets across 15 organs using 6 general assay types."* They do **not** come from the
portal itself, which is a JavaScript single-page app that renders nothing to a fetcher
and will not corroborate the figure. Worse, **v2 of that same preprint already updates
the count to 2,041 as of April 2026**, so the number is drifting under us.
*Fix:* cite the preprint, not `data.sennetconsortium.org`, and keep the "as of January
2026" qualifier, which makes the statement true and durable.

**S15 — do not cross-attribute the "39 genes" and "22 genes" senescence-gene-set
cores.** DeepScence (ref 8) surveys **nine** gene sets and finds **39 genes** reported
by ≥5 of them out of 2,966 total. Ntintas et al. (ref 14) also compares **nine
resources** — a *different* nine (SenMayo, Fridman, Casella; GO BP / Reactome / KEGG
senescence; GenAge, CellAge, SeneQuest) — and finds a core of **22 genes** in ≥5 of
them. Both the "nine resources" and the small-core claims in §31 refs 8 and 14 are
**VERIFIED**, but they are two different measurements and the numbers must not be
swapped. Our own `BIO_PHASE2.md` §4.2 correctly uses DeepScence's 39 (35 on panel).

### 2.3 MINOR

**M0 — truncated or paraphrased titles in §31.** All resolve; all need completing
before submission. Ref 9: the title is not "SenCat" — that is the resource name; the
paper is *"SenCat: Cataloging human cell senescence through multi-omic profiling of
multiple senescent primary cell types."* Ref 10 (SenePy): title ends *"…using SenePy."*
Ref 16 (Basisty): title ends *"…for aging biomarker development."* Ref 24 (SpaTalk):
title ends *"…with SpaTalk."* Ref 19: the real title spells out *"senescence associated
secretory phenotype"*, not "SASP", and the correct form is Meyer P, Maity P, Burkovski A,
et al. *PLoS Comput Biol.* 2017;13(12):e1005741. Ref 12 (SenCID) had **no title or
authors at all**; it is Tao W, Yu Z, Han JJ. *Single-cell senescence identification
reveals senescence heterogeneity, trajectory, and modulators.* Cell Metab.
2024;36(5):1126–1143.e5. doi:10.1016/j.cmet.2024.03.009. PMID 38604170.

**M0b — two apparent year/DOI conflicts in §31 are not errors; do not "fix" them.**
Ref 13 (SenPred): *Genome Medicine* 2025;17(1):2 with DOI `10.1186/s13073-**024**-01418-0`
is correct — BMC/Springer DOI suffixes encode the submission year, not the issue year.
Ref 23 (SCILD): *Communications Biology* 2026;9(1):133 with DOI
`10.1038/s42003-**025**-09413-w` is likewise correct; Crossref and Europe PMC both give
2026-01-07. Both years as written in §31 are right.

**M1 — arithmetic in the README's fit denominator.** `README.md` L~90: "160 reportable
fits (6 admissible sections × 9 receiver types × 7 modules)" implies 378. `CS_PHASE3.md`
§0 says **160 of 315**. Not every receiver type is present in every section. Write
"160 of 315 (section × receiver type × module) fits."

**M2 — 0.082 vs 0.084 in adjacent sentences.** README L~88–95 gives the combined
N2+N5+N6 SF as **0.082** and two sentences later calls "the observed **0.084**" the
number compared to synthetic. Both are correct (`CS_PHASE3.md` §0: N5 alone = 0.084,
N2+N5+N6 = 0.082) but the reader cannot tell. Label them.

**M3 — kernel-family effect on d̂½: 5× vs 4.4×.** README L~116 says "kernel family
moves d̂½ by a median 5×" (`CS_PHASE2.md`, naive fits); `CS_PHASE5.md` §0.5 revises it
to **4.4×** under control. Both are real, of different populations. Say which.

**M4 — COMMOT ρ appears as 0.78 and 0.90 in the same README.** L~92 ("the cell-type-level
summary it actually tests barely moves (ρ = 0.78)") and L~200 ("ρ = 0.90 (COMMOT)").
`CS_PHASE4.md` §4.1 gives per-pair cluster-level Spearman of 0.646/0.872/0.841/0.715
(mean ≈ 0.77) on the mechanism tile; the headline table gives ρ = 0.90 over all 6,032
interactions. Different quantities, same symbol, no distinguishing label.

**M5 — edge Jaccard: 0.015 / 0.012–0.018 / "98.7%".** `CS_PHASE4.md` §4.1 gives
0.0176 / 0.0147 / 0.0152 / 0.0141, i.e. **0.014–0.018**, mean 0.0154 ⇒ **98.5%**
replaced. The README's "0.012–0.018" and "98.7%" are both slightly outside the table.

**M6 — λ̂ grid-rail rate: 66% vs 63%.** README L~115 says 66% (that is `CS_PHASE2.md`,
naive, before annotations landed); `CS_PHASE3.md` §4.1 says **200 of 315 = 63%** on the
final admissible set. Use the Phase 3 number for anything describing the final fits.

**M7 — cell counts differ by ≤1,100 between reports.** 7259: 128,030 (`BIO_PHASE1.md`)
vs 127,386 (`BIO_PHASE2/3`, `CS_PHASE2`). 7250: 237,982 (`BIO_PHASE1.md`) /
236,906 (`CS_PHASE2.md`) / 236,905 (`BIO_PHASE2/3`). Almost certainly pre- vs
post-QC and an off-by-one in a filter; state the QC rule once and use one number.

**M8 — Ripley-K generalisation.** README: *"Ripley-K at 50 µm is 1.11, so senders are
nearly a random thinning."* `CS_PHASE3.md` §4.2 gives 1.109 for `tierA_p95`, **1.263**
for `cdkn1a_pos` and **1.556** for `senepy_p95`. The conclusion holds for all three
(all sit at the bottom of the Phase 1 clustering axis) but the single number is the
most favourable one.

**M9 — composition surrogate range omitted.** README: "reproduces 76% of the contact
amplitude." `CS_PHASE5.md` §4 adds "per-module ratios **0.42–1.83**". The median is a
fair summary but the spread is a factor of four and belongs in the sentence.

**M10 — provenance asterisks are attached to some claims and not others.** The README
carries the `*` footnote correctly at L~200 but the earlier bullet (L~85–100) states
"CellChat 0.318 vs 0.283" and "94% (CellChat)" with the asterisk only on some mentions.
Since `CS_PHASE4.md` §5.1 says "**No result in this phase should be attributed to the
CellChat, SpaTalk or NCEM software**", the marker should be on every occurrence in the
public-facing file, not most.

**M11 — unsourced biological framing sentences.** Two assertions do real rhetorical
work with no citation: `BIO_PHASE2.md` §2.3 *"That is textbook IFALD"* and
`BIO_PHASE3.md` §2.2 *"consistent with the aged-liver ductular reaction literature."*
Both are plausible and I did not find anything contradicting either, but if the second
survives into the Discussion it needs a reference, because it is the sentence carrying
the n = 1 52-week sham animal.

---

## 3. Citations that do not resolve

**None. All 30 references in §31 exist and were located.** Every one was matched to a
publisher, index or preprint-server record. There are no phantom citations in the
reference list.

That said, "located" is not the same as "read end to end," and three items rest on
weaker evidence than the rest. **Failure to reach a full text is not evidence that a
claim is wrong.** Here is exactly where the evidence thins, and what was tried.

### 3.1 Ref 1 (Ma et al., *Cell* 2024) — bibliography verified, distance-gradient claim substantiated but not rendered

The **bibliographic record is fully verified** (Europe PMC core record, DOI
`10.1016/j.cell.2024.10.019`), and that record is what establishes the first-author
error. What could **not** be confirmed from a page rendered in this session is the
*content* attribution: that SASP score, TNF signalling, ATP biosynthesis and cell-cycle
genes all vary monotonically with distance from senescence-sensitive spots.

*Tried:* WebFetch on `cell.com/cell/fulltext/S0092-8674(24)01201-7` (**HTTP 403**, three
attempts, different prompts); WebFetch on
`sciencedirect.com/science/article/pii/S0092867424012017` (**403**); Europe PMC core
record (confirms `isOpenAccess: No`, `inEPMC: No`, `inPMC: No` — there is no free full
text to reach); `db.cngb.org/stomics/datasets/STDS0000247` (dataset metadata only);
the MGI-tech commentary page (does not cover the distance analyses); Google Scholar on
`"senescence-sensitive spots" "ATP biosynthetic"` (snippet only); plus five WebSearch
queries targeting the exact phrases.

*What came back:* multiple independent search retrievals surfaced what reads as verbatim
cell.com full text — *"proximity to the SSS correlated with higher SASP scores"*;
*"TNF-signaling gradually decreases with increasing distance from the SSS"*;
*"expression levels of genes associated with ATP biosynthetic process and cell cycle
were restored with increasing distance from the SSS"*; and that SSSs *"serve as
epicenters for heightened inflammation that compromises surrounding cells in a
distance-dependent manner"*, observed consistently across organs.

**Verdict: substantiated, not personally verified.** All four attributed gradients and
the "epicentre" framing appear. Since this is our closest prior work and the sentence
our Introduction is built against, **someone should pull the PDF through Caltech's
library and confirm the figure panel** before camera-ready.

### 3.2 Ref 27 (CellWHISPER) — exists beyond doubt; the >90% sentence not rendered

**Existence is not in question.** Confirmed via three independent non-HTML routes:
`api.biorxiv.org/details/biorxiv/10.64898/2026.01.07.697982` (both versions, full
metadata and abstract), Europe PMC (`PPR1216598`, bioRxiv, 2026-01-08), and Semantic
Scholar (`CorpusId 284585047`).

*What could not be reached:* the full text. **bioRxiv returned HTTP 429 to every
attempt** — roughly fifteen tries across five URL paths (`v1.full`, `v1.full.pdf`,
`v2.full`, the `/content/biorxiv/early/2026/01/08/…full.pdf` path, and the JATS
`.source.xml`), two user agents, both WebFetch and curl, with backoff spacing of 30–45 s
over about ten minutes. This is a shared-IP rate limit, not a bad URL. **CellWHISPER is
in neither PMC nor Europe PMC full text** (`inEPMC: N`, `hasPDF: N`,
`isOpenAccess: N`), so there is no mirror.

*What came back:* two independent search retrievals, run with different query phrasings,
both returned the search engine's extraction of the bioRxiv full text and both
independently reported: *"CellChat v2, COMMOT, and SpaTalk predicted similar numbers of
interactions on real and randomized data, indicating poor specificity and false positive
rates (FPR) >90%"* and *"CellWHISPER produced markedly fewer interactions on randomized
data compared to real data, suggesting FPR < 5%."* The three tool names and both
percentages match what we claim. The abstract, which **was** fetched directly, is fully
consistent (*"prone to high error rates due to spatial structure in distribution of cell
types and gene expression"*; *"achieves strict error control"*).

**Verdict: substantiated, one notch below the standard applied to CONCISE.** Two
independent retrievals of the primary source agreeing on all three numbers is good
evidence. It is not the same as having loaded the page. **Given that this is the paper
the entire framing rests on — it is the first citation in our README's opening
paragraph — one human should open that PDF once from an unblocked network.**

*Also flagged:* a **name collision**. `CellWhisperer` (epigen/CellWhisperer, Bock lab,
*Nature Biotechnology* 2025, doi:10.1038/s41587-025-02857-9) is an unrelated multimodal
chat tool for scRNA-seq. Different paper, different capitalisation. Check the .bib does
not cross-contaminate.

### 3.3 Ref 17 (Acosta et al. 2013) — the paper verified, its 2026 correction not read

The paper and the transwell claim are **verified** from `PMC3732483`. What could not be
read is the **Author Correction** *Nature Cell Biology* published on **2026-04-28**
(doi:10.1038/s41556-026-01959-z, PMID 42050148). Its existence is confirmed via Crossref
and Europe PMC; its **content is not**. Tried: nature.com (redirects anonymous fetches to
`idp.nature.com`); the Europe PMC core record for the correction DOI (returns the title
*"Author Correction: A complex secretory program orchestrated by the inflammasome
controls paracrine senescence"* and the target article, but **no abstract, no PMCID, no
full text**).

**I am not asserting anything is wrong with Acosta et al.** I am flagging that a
correction exists on a foundational citation and that nobody on this project has read
it. It needs one human with library access.

### 3.4 Live portal state — could not be read

`data.sennetconsortium.org` is a JavaScript single-page app that renders "Loading,
please wait…" to a non-executing fetcher, and its search API endpoint 404'd at the
guessed path. **I cannot tell you what the portal shows today.** The 1,753 figure was
verified against the SenNet Portal preprint instead (see §2.2 S18), which is where it
actually comes from.

---

## 4. Causal-language flags

Deliverable 7's explicit remit. My rule for judging: **is the causal verb licensed by
an intervention (including a simulation where the mechanism is planted by
construction), or only by an adjustment/conditioning result on observational data?**
Planted-truth simulation earns causal language. Covariate adjustment does not.

| # | Claim, as written | Where | Evidence | Earned? | Proposed hedge |
|---|---|---|---|---|---|
| **C1** | "local cell density **mechanically confounds** distance-to-nearest-sender" | `CS_PHASE1.md` §1 | Synthetic tissue with the density effect **planted by construction**; λ̂ biased +45% with random senders and no latent field | **YES — fully earned.** In the simulator this is a mechanism, not an inference. Leave it. | none |
| **C2** | "This is the Phase 1 mechanism … **appearing in real tissue**" | `CS_PHASE3.md` §3.2 | Conditioning result: the density sub-block alone leaves SF 0.219 | **NO — this is adjustment, not identification.** The synthetic mechanism is *consistent with* the real-data pattern; it is not shown to be the operative one. | "the real-data pattern is what the Phase 1 density mechanism would produce: conditioning on local density alone removes 78% of the amplitude" |
| **C3** | "zonation **was never driving it**" | `README.md` L143 | Zonation covariate alone leaves SF 1.043 in hepatocytes (0.843 overall, 0.244 in biliary) | **NO — a null adjustment result cannot establish that a variable was never causal.** It establishes that the covariate as constructed removes nothing. Our zonation score is itself a derived marker score with known limits (`BIO_PHASE1.md` §4.1 flags `Pigr` bile-duct contamination). | "conditioning on our zonation score removes essentially none of the hepatocyte amplitude, so zonation as we can measure it is not the operative confound here" |
| **C4** | "The confound **is** technical and geometric" | `README.md` L143; `CS_PHASE3.md` §3.2 heading | Sub-block attribution: tech 0.288, dens 0.219, comp 0.474, anat 0.810, seg 0.998 | **PARTIAL.** The ordering is a real, replicated result and the language is defensible *as an attribution*, but "is" claims identification. Note also the cumulative sequence (0.288 → 0.303 → 0.085 → 0.044 → 0.056) is order-dependent, which the report says but the README does not. | "the amplitude is removed principally by technical and geometric covariates and only marginally by anatomical ones" |
| **C5** | "**The statistic is not spatial**" | `README.md` L96, L205 | ρ(real, N0) = 0.86–0.98 for three methods; COMMOT code read directly | **PARTIAL, and self-contradicted** — see S1. Earned for the *cluster-level summary*, and for COMMOT it is earned mechanistically from the released code. Not earned as a blanket statement, and two of the three methods are reimplementations. | see S1 |
| **C6** | "The established tools **cannot tell this tissue from confetti**" | `README.md` L85 | Significance survival 79/94/79% under N0 | **NO, as phrased.** Three of four are reimplementations; COMMOT ran on 1.2 mm tiles, not sections; and `CS_PHASE4.md` §5.7 says survival is not an FPR. The rhetorical compression drops all three qualifications at once. | "On these tiles and these four ligand–receptor pairs, three of the four statistics reproduce most of their real-data calls on data with no spatial information left in it" |
| **C7** | "The torus shift **certifies** confounding as real" | `README.md` L~72 | Synthetic: β_true = 0 in the confounded regime, N3 returns 98% SF | **YES — earned.** Planted ground truth; "certifies" describes what the null does, not a claim about tissue. Leave it. | none |
| **C8** | "The regressor **is** largely a sender-calling-rate readout" | `README.md` L~228; `CS_PHASE3.md` §4.2 | slope −0.524, r² = 0.984 across 77 combinations against the Poisson identity 0.4697ρ^(−1/2) | **YES — earned.** This is a geometric identity confirmed empirically, not a causal inference from covariates. Leave it. | none |
| **C9** | "This **applies to every method** in this literature that regresses response on distance-to-nearest-X" | `README.md` contribution 4 | one dataset, one panel | **PARTIAL.** The Poisson identity generalises; the empirical r² and the depth coupling do not. See B6. | "The geometric part of this — that median distance-to-nearest-X is pinned to X's calling rate by the Poisson identity — applies to any method regressing response on distance-to-nearest-X; the strength of the coupling in a given dataset is empirical" |
| **C10** | "It is **not** `Cdkn1a` circularity — the effect is *stronger* in the modules that exclude `Cdkn1a`" | `README.md` L~250 | Module-level contrast plus `senepy_p95` reproducing the split | **PARTIAL.** A negative-control result of the right shape, but circularity has other routes than `Cdkn1a` alone — `BIO_PHASE2.md` §4.2 documents 24 of 35 CoreScence genes sitting inside Tier B modules. The claim rules out one route. | "`Cdkn1a` circularity does not explain it: the effect is stronger in modules that exclude `Cdkn1a`, and `senepy_p95` reproduces the split. Other circularity routes are not excluded." |
| **C11** | "The negative **is not** a selection artefact" | `README.md` L128, L186 | Cross-fit with a matched no-selection placebo; contamination +0.0563 | **YES — earned, and unusually well.** The placebo design is an intervention on the selection step, and the sign is reported against interest. Leave it. | none |
| **C12** | "This is a finding the plan did not anticipate and **it applies to every real dataset**" | `CS_PHASE1.md` §1 | one simulator | **NO.** Universal generalisation from a single generative model whose density–response coupling was planted. | "…and it will apply to any dataset in which local cell density is associated with the response, which we expect to be common" |
| **C13** | "The ductular reaction … is recovered by a pipeline that **knows nothing about** arm or timepoint" | `BIO_PHASE3.md` §2.1 | Blind-annotation composition table | **YES — earned.** This is the correct form of the argument (the classifier had no access to the labels) and it is stated precisely. Leave it. | none |
| **C14** | "SenePy's Kupffer enrichment **is** a hub-size artefact" | `BIO_PHASE3.md` §4.2 | mean score vs on-panel hub size, r = 0.992 across five cell types | **PARTIAL.** r = 0.992 on **n = 5 points** is a strong pattern but a weak sample, and hub size may covary with cell-type transcriptional complexity. The operational recommendation (threshold within cell type) is right regardless. | "SenePy's cross-cell-type score scales almost perfectly with the number of that hub's genes on our panel (r = 0.992, n = 5 cell types), so it is not comparable across cell types" |
| **C15** | "DeepScence's CoreScence set **is** 69% circular with our response modules" | `README.md` L~300; `BIO_PHASE2.md` §4.2 | ~~24 of 35~~ **26 of 33** on-panel CoreScence genes are members of ≥1 Tier B module | **NO on the number as audited, YES on the claim.** ⚠ **Superseded 2026-08-27** (`AUDIT_PHASE8_FACTCHECK.md` M1): this row verified `24/35` against the report, not against files. The denominator 35 is reproducible under no mapping convention; 31 of the 39 CoreScence genes are reachable on the mouse panel through the pinned MGI map and 33 with the documented Title-case fallback. The circularity claim survives and strengthens: **79%**, not 69%. "Circular" is the (correct) interpretation but should be introduced as such once. | "…shares 26 of the 33 CoreScence genes reachable on our panel (79%) with at least one Tier B response module — i.e. scoring senders with DeepScence and reading out a Tier B response is partly circular by construction" |

**One general observation, and it is the most useful thing in this section.** The phase
reports are consistently well hedged — `CS_PHASE4.md` §5, `CS_PHASE5.md` §6.4,
`BIO_PHASE3.md` §4.3 and §5 all volunteer limitations against their own interest, and
`CS_PHASE5.md` §1 records a bug that changed a headline. **Almost every overreach in
this audit is introduced by the README**, which compresses hedged report findings into
unhedged headline sentences. The fix is therefore concentrated: audit the README against
the reports, rather than re-litigating the reports.

---

## 5. Priority order for the CS lead

Ordered by what breaks if it is not fixed, not by effort.

| Rank | Item | Type | Why first |
|---|---|---|---|
| 1 | **B1** — "Zhao et al." → **Ma et al.** | Citation | Wrong first author on the closest prior work, cited in the README's opening sentence. Costs one line to fix; costs the paper's credibility if it ships. |
| 2 | **B9** — N0_perm mislabelled as CellWHISPER's null | Method attribution | Figure 4's headline is a replication claim about a null we did not run. Either run the within-type permutation (cheap, and it strengthens the figure) or relabel. |
| 3 | **B4** — README asserts ">90% FPR" that CS_PHASE4 §5.7 says is not an FPR | Internal contradiction | The paper's central argument is that shuffled tissue is not a valid null. We cannot use survival-on-shuffled-data as an FPR while arguing that. |
| 4 | **B5** — README says estimation is SBR-only; it is 2 SBR + 4 sham | Internal contradiction | README and Methods disagree about which arm the paper is about. |
| 5 | **B4b** — Martin et al. containment mechanisms are not the paper's | Citation content | The Discussion is written against this paper. No *result* depends on it (CS_PHASE5 §6.4 reported a bound), but the prose does. |
| 6 | **B6** — ρ = 0.94 depth coupling promoted to a contribution | Overreach | n = 6, one arm, hedge removed, and "sequencing depth" is wrong for an imaging assay. |
| 7 | **B7** — DeepScence claims without the `denoise=False` caveat, on 2 of 11 sections | Provenance | BIO_PHASE3 §4.3 asked for specific words; the README uses none of them. |
| 8 | **B8** — ref 20 is a GSA meeting abstract | Citation | Sole source for the primary/secondary distinction underwriting module B7. |
| 9 | **B2, B3** — Ma et al. is Stereo-seq spot data, two age groups | Citation content | Both *help* our framing. Fixing them is free and makes the gap argument stronger. |
| 10 | **S0a–S0d, S13–S18** | Citations | Year, journal, title, and attribution corrections across nine references. Mechanical. |
| 11 | **S1–S12, C1–C15** | Hedging | Concentrated almost entirely in the README (see §4's closing note). |
| 12 | **M0–M11**, container-disk pricing | Consistency | Cheap, and the numeric inconsistencies are the kind a careful reviewer greps for. |

### What this audit did **not** find

Worth recording, because negative audit findings are also findings.

- **No fabricated references.** All 30 exist.
- **No fabricated numbers in the results.** Every headline figure I traced from the
  README resolved to a table in a phase report. The discrepancies in §2.3 are rounding,
  stale-version and different-population issues, not invention — with the single
  exception of "~1.3M cells" (S6), whose provenance I could not find at all.
- **No surviving human-gene-symbol errors in the results pipeline.** The species switch
  was handled carefully: `BIO_PHASE3.md` §3.3–3.4 and `CS_PHASE4.md` both state
  correctly that mouse has no CXCL8 orthologue and reason from `Cxcl1`/`Cxcl2`; the
  human symbols that remain (`BIO_PHASE2.md` §4.2) are correctly marked as DeepScence's
  own human core set. The two exceptions are cosmetic and listed at S11.
- **No case where the human-derived endometrium length-scale calibration leaked into a
  mouse result.** It appears only in the master plan; no report invokes it.
- **The gene-set provenance is clean.** `BIO_PHASE1.md` §2.2's claim that gene sets are
  real MSigDB 2026.1.Mm collections, fetched live and archived, with nothing hand-mapped
  or recalled from memory, is exactly the discipline §9 asked for and is the reason
  there is nothing to audit there.
- **The phase reports hedge well.** `CS_PHASE4.md` §5, `CS_PHASE5.md` §0.6 and §6.4,
  `BIO_PHASE3.md` §2.2, §4.3 and §5 all volunteer limitations against their own
  interest, and `CS_PHASE5.md` §1 records a bug that changed a headline. **The overreach
  is in the README, not in the reports.** Audit the README against the reports and most
  of §2.2, §2.3 and §4 resolves at once.

---

*Deliverable 7 complete. Nothing in this audit was applied to any other file; every
change is the CS lead's call.*
