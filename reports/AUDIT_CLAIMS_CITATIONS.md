# Audit — claims and citations, 2026-08-27 (post-correction re-verification)

**Scope.** (1) re-verify all 43 `references.bib` entries against retrieved Crossref
and PubMed records, independently; (2) verify the new spatial-statistics block
including each `% SUPPORTS:` line; (3) sweep every report and planning document for
the three falsified claims and the forbidden wordings; (4) test each
superseded-marker's scope against the text beneath it.

**Method.** Crossref REST (`api.crossref.org/works/<doi>`) for all 43; NCBI E-utilities
`efetch` for all 25 entries carrying a PMID; DataCite + the arXiv abs page for the one
arXiv DOI; Europe PMC for all 10 PMCIDs; and, for the `SUPPORTS` lines, the primary
sources themselves — the Mrkvička arXiv PDF, PMC full texts for Dupont and Ren, the
Voyager Xenium vignette, the `spatstat.random` `rshift.ppp` help page, the CellWHISPER
v1 full text, and the OpenAlex abstract record for Hodges & Reich. Nothing below rests
on recall.

**Read-only run.** No file other than this one was modified.

---

# 1. LEAD — live forbidden claims still in documents

Ordered by what a reviewer would hit first.

### L1. `references.bib:783` — the bibliography itself instructs the writer to make claim (b)

The `% SUPPORTS:` annotation on `moses2023voyager` ends:

> Report our own Moran's I on the controls next to the kernel amplitude so the reader
> can see **the two tests disagree**.

That is forbidden-claims item **21** verbatim, and it is the *same sentence* that
`WRITING_PACK.md` §3.1(b) records as struck from `NOVELTY_ASSESSMENT.md` §2.1 point 3.
It is unmarked, it is phrased as an instruction, and it sits in the block that was
**added by the citation audit on 2026-08-27** — i.e. the falsified sentence was
re-introduced into a new file after being struck from the old one. Measured, the two
statistics rank together at ρ = +0.895 / +0.944 (`results/moran/moran_verdict.txt`).
The rest of that entry's SUPPORTS text is correct and verified; only this clause is dead.

### L2. `reports/BIO_PHASE3.md` — carries all four struck sentences, with **no correction marker anywhere in the file**

`grep` for `SUPERSEDED|superseded|Corrected 2026|⚠|WITHDRAWN|pre-C6` over the whole
file returns **zero** matches. Live:

- §0 headline item 4 (L26): *"conditioning on cell type **and** depth decile does **not**
  make them agree (Jaccard ratios 0.93–1.22× of chance)"* — items 13/14.
- §4.4 (L408): *"Four of six pairs sit at **0.93–1.22× of chance** — statistical
  independence"* — items 13/14.
- §4.4 (L410–411): *"The one pair that looked concordant in sham … is **anti**-concordant
  in SBR"* — item 15.
- §4.4 table (L406): the `1.51 / 2.85` cell that is the origin of the forbidden
  circularity range — item 8.
- **§4.5 (L417–424), headed "The finding, stated for the paper (§9 A7, §10)"** — a
  blockquote drafted as manuscript text containing, in one paragraph, *"overlap at
  0.93–1.22× of chance … i.e. they are **statistically independent**"* **and**
  *"DeepScence's correlation with sequencing depth **reverses sign between two sections**
  of the same study"* — items 13 and 16 together.

This is the single largest concentration of forbidden text in the repo, and it is the
only such document with no banner at all. `PREREG_PHASE8.md` §10.12 strikes exactly
these four sentences; nothing in `BIO_PHASE3.md` says so.

### L3. `reports/BIO_PHASE2.md` §4.3 — the independence claim in bold, unmarked

L239: *"**The three definitions are statistically independent — they are not noisy
measurements of one latent senescent state.**"* Plus §0 headline item 2:
*"the sender callers do not agree with each other at better than chance"*, and the §4.3
heading *"⚠️ The sender callers agree at or below chance"*. The ⚠️ marks are the
authors' own emphasis on the finding, not correction markers — the file's only real
correction box is in §4.2 and covers the 69 % → 79 % circularity denominator. Forbidden
item **3** / **13**.

### L4. `reports/CS_PHASE8_M1_RERUN.md:307` — the A7 response-naming error, in a row marked **FILL IN**

> \| **Negative-control-probe kernel (must be flat)** \| *new* \| **NOT FLAT**: −0.074 SD,
> p = 0.015 naive; N5 removes it (+0.004, p = 0.72), **N2 does not** (−0.064, p = 0.012).
> Measured FPR **9–16 %** against a 5 % nominal \| **FILL IN — and it is a negative result
> about the platform** \|

−0.074 is the pooled `all_controls` response. Calling it a *probe* kernel and asserting
**NOT FLAT** is forbidden item **2** / `PREREG_PHASE8.md` §10.1 / audit R1: the 40 named
probes are the pre-registered primary null and are **flat** (−0.0225, p = 0.129), so on its
own primary response M1's A7 **passes**. The row is explicitly staged for the manuscript
("FILL IN"), and the document carries no banner covering it.

### L5. `reports/NOVELTY_ASSESSMENT.md` §U2 (L332) — the Mrkvička misquotation, in quotation marks

> Mrkvička *et al.* (2021) … proposes the **variance correction** explicitly *because* it
> "can also be used for irregular windows"

Forbidden item **32**. The `mrkvicka2021revisiting` entry's own CAUTION line names this
exact location. I extracted the arXiv:1911.00240 PDF: the paper's §2.1.4 sentence is
*"The approaches using minus correction and variance correction can be applied in case of
general (compact) observation windows."* The quoted wording does not occur in the paper.

### L6. `reports/NOVELTY_ASSESSMENT.md` §U1 (L300–307) — a drop-in manuscript paragraph carrying claim (c) and pre-C6 digits

> Negative control probes are routinely used to quantify background rate … We instead use
> **them** as a negative-control outcome for the estimand itself [Lipsitch 2010]: we refit
> the identical distance-to-sender kernel with control counts as the response. The naive
> estimator returns **−0.070 SD (p=0.023)** …

The antecedent of "them" is *negative control probes*; the number is the pooled
`all_controls` response, and it is the **pre-C6 05:19** vintage (frozen: −0.0744,
p = 0.0145). Forbidden item **2**, plus a stale vintage. `CORRECTIONS.md` §17.4 already
flagged this ("**⚠ The banner ends 'Everything else in this report stands.' It does not.**")
and it has not been fixed.

### L7. `Phase7_Minimal_Human_Replication (1).md:187` — "1.51–2.85×" live below a banner that does not reach it

§4 D-d: *"The one pair consistently above chance is DeepScence vs `Cdkn1a`⁺ at
**1.51–2.85×**"* — forbidden item **8**. See §4 below: the file's SUPERSEDED banner names
§4 only in respect of item **D-b**.

### L8. `SASP_Kernel_Master_Plan.md:284` (§7 Rank 2) — the unverifiable SenNet figure, still live

> 1,753 public human and mouse datasets across 15 organs and 6 assay types **as of
> January 2026**.

`BIO_DELIVERABLE7_CLAIM_AUDIT.md` names **two** locations for this figure — "§4.1 **and
§7 Rank 2**" — and only §4 was corrected (L182–184). I re-checked the bioRxiv API today:
the abstracts of **both** v1 (2026-02-10) and v2 (2026-05-04) read *"As of April 2026, the
portal hosts 2,041 publicly available human and mouse datasets across 15 organs using 6
general assay types."* No retrievable version states 1,753 or January 2026.

### L9. Smaller live remnants

| where | what | should read |
|---|---|---|
| `README.md:331` | table header *"Moves across callers that **agree at chance**"* — three paragraphs below the correction box that restates them as *above* chance | "callers that agree weakly but genuinely above chance" |
| `Phase7_Minimal_Human_Replication (1).md:21` | the SUPERSEDED **banner itself** asserts DCA *"**roughly doubles** the depth loading on 3 of 3 sections"* | ×1.32–1.67, i.e. **32–67 %** (`PREREG_PHASE8.md` P29; `CORRECTIONS.md` §1739) |
| `reports/CS_PHASE8_C1_CLOSEOUT.md:195, 237, 294` | *"the median λ̂ (12.8 µm)"* ×3 | **14.7 µm** (`PREREG_PHASE8.md` §0.0 C-7, which claims the withdrawal is "recorded project-wide") |
| `reports/CS_PHASE8_C1_CLOSEOUT.md:234–236` | *"1 to 63 of 38,080 to 108,375"*; *"median displacement 27 µm (N3)"* | **1–66**; **28 µm** (C-6, and `AUDIT_PHASE8_FACTCHECK.md`) |
| `reports/WRITING_PACK.md:1280` (§4.1), `:893`, `:896`, `:1513`; `reports/PLAN_UPDATE_D12_D13.md:216, 277` | tiled-torus inflation quoted as **2.4×** | the repo standardised on **2.35×** (0.1175 / 0.05 exactly) at `WRITING_PACK.md:1481` and `CORRECTIONS.md:131` — the writing pack contradicts itself in four places |
| `reports/BIO_DELIVERABLE6_DISCUSSION.md:72` | *"agree at or below chance (Jaccard **0.93–1.22×** of chance …)"* | see §4 — the file's only banner is scoped to the pre-C6 SF vector |

### Checked and clean

The following are **not** defects — every instance is either quoted-and-retired, correctly
scoped, or a different quantity:

- *"the two statistics disagree"* at `WRITING_PACK.md:657` and `CS_PHASE8_MORAN.md:176` —
  both are the narrow, correct, reported-against-interest statement that `genomic_control`
  is A7's third-largest amplitude but Moran's smallest I. Explicitly scoped to one
  response's rank. Keep.
- *"at chance (1.025, z = 1.4, n.s.)"* (`SUBMISSION_PATCH:180`, `WRITING_PACK:1042`,
  `CS_PHASE8_CALLERS:598`, `COMPLETED_TASKS:510`, `CORRECTIONS:1136`) — the |score|-ranking
  sensitivity result, a measured finding, not the independence claim. Keep.
- **λ̂ = 15.7 µm is fully swept.** Every surviving instance is a withdrawal record.
- **Brackets as CIs**: `PREREG_PHASE8.md` §5 (PROVISIONAL + C-5) and §6 R1 (inline C-4
  correction: *"'Paired-bootstrap interval' is a misnomer … read R1 as: the IQR"*) are both
  properly labelled, as are the master plan §30 5.2/5.3, `SUBMISSION_PATCH` §8,
  `PLAN_UPDATE`, `RECORD_RECONCILIATION` §157 and `WRITING_PACK` §0.7/§185. I found no
  unlabelled SF or amplitude bracket presented as a confidence interval.
- `CS_PHASE8_CALLERS.md`'s probe-naming error is **fixed** — `CORRECTIONS.md` §17.5's
  "`:271, :283, :290, :348, :435` … still live and still needs fixing" is now **stale**;
  the only surviving mention (`:477`) is correct.
- `NOVELTY_ASSESSMENT.md` §2.1 point 3 **is** struck in the body (`~~…~~` + a replacement
  paragraph). `WRITING_PACK.md` §3.1(b)'s parenthetical "still unmarked in that document's
  body" is stale.
- `PREREG_PHASE8.md` §10.1 para 1's pre-C6 digits carry an inline `[Corrected 2026-08-27 …
  Frozen values: …]` bracket. Covered.

---

# 2. Citation defects — wrong, not merely unverified

### C1. `ren2025systematic` — **WRONG last author**

```
author  = {Ren, Pengfei and Zhang, Rui and Wang, Yunfeng and others and Zhang, Zongxu},
```

The `… and others and X` convention asserts **X** as the final author. Crossref
(doi:10.1038/s41467-025-64292-3) and the PMC12534522 `citation_author` metadata both list
19 authors ending **Zhe Zhang; Zexian Zeng**. **Zongxu Zhang is 8th of 19.** The bib names
the wrong senior author on the only peer-reviewed prior-art citation for the
negative-control-probe diagnostic — the citation `WRITING_PACK` item 20 and
`NOVELTY_ASSESSMENT` §U1 both lean on. Also incomplete: no `number = {1}`, no
`pages`/article number (**9232**), no `pmid` (**41107232**, now confirmed).

Correct form:

```
author = {Ren, Pengfei and Zhang, Rui and Wang, Yunfeng and others and Zhang, Zhe and Zeng, Zexian},
number = {1},
pages  = {9232},
pmid   = {41107232},
```

### C2. `hughes2025senpred` — the **audit annotation** is wrong (the entry is right)

> `% AUDIT: title is INCOMPLETE here -- the audit confirmed the record but did not supply
> the full title. Complete it before camera-ready.`

The title in the entry matches Crossref and PubMed (PMID 39810225) **word for word**,
including the trailing "…for the detection of an in vivo senescent cell burden". Nothing
is missing. The note should be deleted; leaving it standing invites someone to "fix" a
correct title.

### C3. `hodges2010adding` — the "NOT VERIFIED" caveat is now over-cautious

The entry says its content is confirmed only indirectly via Dupont. I retrieved the full
abstract today from the OpenAlex record for doi:10.1198/tast.2010.10052 (Taylor & Francis
still returns 403). It states verbatim: *"We show how to avoid this spatial confounding by
restricting the spatial random effect to the orthogonal complement (residual space) of the
fixed effects, **which we call restricted spatial regression**."* The `SUPPORTS` claim —
"this is the paper that named it and proposed restricted spatial regression" — is
**directly confirmed**. The Dupont cross-reference is also verbatim-correct
(*"spatial confounding is widely acknowledged as an issue that affects spatial models in
general (see, eg, Hodges and Reich, 2010; Paciorek, 2010)"*, PMC10084199 §1).

### No other citation is wrong

42 of 43 entries reconcile on authors, title, venue, volume/issue/pages and DOI against a
retrieved record. Full table in §3.

---

# 3. Per-entry verdicts (43 entries)

`CR` = Crossref record retrieved and reconciled; `PM` = PubMed `efetch` record retrieved
and reconciled independently; `PMC` = PMCID resolved via Europe PMC.

| key | verdict | evidence |
|---|---|---|
| `ma2024spatial` | **CONFIRMED** | CR+PM. Cell 187(24):7025–7044.e34. 47 authors; the bib's tail `Qu, Jing / Zhang, Weiqi / Gu, Ying / Liu, Guang-Hui` is exactly Crossref's last four. First author **Ma S** ✓ |
| `gurkar2023spatial` | **CONFIRMED** | CR+PM. Nat Aging 3(7):776–790. PM PublicationType **Review** ✓ (bib note "SenNet review") |
| `karpova2026cellular` | **CONFIRMED** | CR+PM+PMC12903365. Cell Genomics 6(2):101133 |
| `nguyen2024scdot` | **CONFIRMED** | CR+PM. Genome Biol 25(1), art. 288 |
| `farzad2026spatial` | **CONFIRMED** | CR only (no PubMed record). Cell Press Blue 1(4):100053, issued 2026-07; 33 authors, head `Farzad / Enninful / Lu` ✓ |
| `suryadevara2026charting` | **CONFIRMED** | CR+PM. Cell 189(12):3501–3505. PM types = `Journal Article`, **`Review`** — the CIT-4 resolution is correct |
| `parvanov2025spatial` | **CONFIRMED** | CR+PM. Diagnostics 15(21):2679, issued **2025**-10-23 ✓ (the B-year fix holds) |
| `qu2025deepscence` | **CONFIRMED** | CR+PM. Cell Genomics 5(12):101035 |
| `anerillas2026sencat` | **CONFIRMED** | CR+PM. Mol Cell 86(13):2605–2616.e8; 29 authors, tail `Basisty / Gorospe` ✓ |
| `sanborn2025senepy` | **CONFIRMED** | CR+PM. Nat Commun 16(1), art. 1884; all 5 forenames exact |
| `saul2022senmayo` | **CONFIRMED** | CR+PM. Nat Commun 13(1), art. 4827. PM gives "Robyn **Laura** Kosinsky"; bib's "Robyn L." is a correct abbreviation |
| `tao2024sencid` | **CONFIRMED** | CR+PM. Cell Metab 36(5):1126–1143.e5. "Jing-Dong Jackie Han" = CR "Jing-Dong J. Han" ✓ |
| `hughes2025senpred` | **CONFIRMED** (see **C2**: its own AUDIT note is wrong) | CR+PM. Genome Med 17(1), art. 2; title complete and exact |
| `ntintas2026overview` | **CONFIRMED** | CR+PM. FEBS Open Bio 16(5):821–836; all 7 forenames exact |
| `martin2023modelling` | **CONFIRMED** | CR+PM+PMC10410058. Aging Cell 22(8):e13892. **Lucy** Martin ✓ |
| `basisty2020proteomic` | **CONFIRMED** | CR+PM. PLoS Biol 18(1):e3000599 |
| `acosta2013complex` | **CONFIRMED** | CR+PM+PMC3732483. Nat Cell Biol 15(8):978–990 |
| `acosta2026correction` | **CONFIRMED** | CR+PM. Nat Cell Biol 28(6):1343; PM type **Published Erratum** ✓ (content still unread — nature.com auth wall; the bib's flag stands) |
| `russo2026characterizing` | **CONFIRMED** | CR+PM. Aging Cell 25(8):e70673; 2 authors ✓ |
| `meyer2017model` | **CONFIRMED** | CR+PM. PLoS Comput Biol 13(12):e1005741 |
| `neretti2024dissecting` | **CONFIRMED as an abstract** | CR: Innovation in Aging 8(Supplement_1):351–351, **single author**; PMC11689308 resolves with no DOI/PMID. The DO-NOT-CITE-AS-A-PAPER flag is accurate |
| `fischer2023ncem` | **CONFIRMED** | CR+PM. Nat Biotechnol 41(3):332–336; CR issued 2022-10-27 (online) / issue year 2023 — the bib's note is right |
| `cang2023commot` | **CONFIRMED** | CR+PM. Nat Methods 20(2):218–228 |
| `yu2026scild` | **CONFIRMED** | CR+PM. Commun Biol 9(1), art. 133; all 5 forenames exact |
| `shao2022spatalk` | **CONFIRMED** | CR+PM. Nat Commun 13(1), art. 4429 |
| `gong2026rgast` | **CONFIRMED** | CR. Brief Bioinform 27(3), art. **bbag298**, issued 2026-05 — the `bbag`=2026 reasoning checks out |
| `wang2026harmonic` | **CONFIRMED** | CR posted-content, openRxiv, DOI prefix **10.64898** ✓, posted 2026-01-23 ✓; 11 authors, tail `Mou / Li` ✓ |
| `kumar2026cellwhisper` | **CONFIRMED — including every annotation** | CR + bioRxiv API + **v1 full text rendered in-session**. Authors `Anurendra Kumar, Felix Rivera, Bhavay Aggarwal, Nicholas Zhang, Ahmet (F.) Coskun, Saurabh Sinha` confirmed against both the Crossref deposit and the page's own `citation_author` metadata. v1 title is the ALL-CAPS "CELLWHISPER: INFERENCE OF DIRECT CELL-CELL COMMUNICATION…"; v2 (2026-04-03) is the retitled form the entry uses ✓. See §5 for the quote-level verification |
| `zhao2026concise` | **CONFIRMED** | CR + PM (PMID 42395397, type `Preprint`) + PMC13320749; 7 forenames exact |
| `yao2026benchmarking` | **CONFIRMED** | CR. Brief Bioinform 27(2), art. bbag190. `Hui Yao / Shuai Mu / Fei He / Zhaoyuan Fang` exact — the CIT-2 forename fix holds |
| `gu2026identifiability` | **CONFIRMED** | arXiv abs page + **DataCite** (`10.48550/arxiv.2607.01749`). Crossref 404 is expected and is not a defect — arXiv DOIs are registered with DataCite. Authors `Gu, Rujie / Zhang, Ray Zirui / Miles, Christopher E.`; submitted 2026-07-02; subjects q-bio.QM, physics.bio-ph, stat.ML — all exactly as the note states |
| `borner2026sennet` | **CONFIRMED**; CIT-3 caveat re-verified | CR posted-content + PMC12918883 (PMID **41727059**, absent from the entry). Title = deposited title ✓. Both v1 and v2 abstracts say **2,041 / April 2026**; neither says 1,753 or January 2026 |
| `lotwick1982methods` | **CONFIRMED (record)** | CR: JRSS-B 44(3):406–413, 1982, doi 10.1111/j.2517-6161.1982.tb01221.x. Crossref's abstract mentions only second-moment estimation and empty-space methods — exactly as the entry's caveat says. **Full text still not rendered** (paywall); the p.410 quotation remains Mrkvička's |
| `mrkvicka2021revisiting` | **CONFIRMED — all four SUPPORTS claims verbatim** | CR + arXiv:1911.00240 PDF extracted in-session. See §5 |
| `hodges2010adding` | **CONFIRMED — and now directly, not second-hand** (see **C3**) | CR + OpenAlex abstract + Dupont's citing sentence |
| `dupont2022spatialplus` | **CONFIRMED — both quotes verbatim** | CR + PMC10084199 abstract rendered. See §5 |
| `khan2022restricted` | **CONFIRMED — both quotes verbatim** | CR (JASA 117(537):482–494) + Semantic Scholar abstract |
| `zimmerman2022deconfounding` | **CONFIRMED — now directly** | CR (Amer Statist 76(2):159–167) + Semantic Scholar abstract, which states the inferiority result the SUPPORTS line claims. The entry's "search-result summaries" caveat can be upgraded |
| `baddeley2015spatial` | **CONFIRMED (record)** | CR book record, Chapman and Hall/CRC 2015, doi 10.1201/b19708, 3 authors. Book text not consulted — the entry says so and redirects to `spatstat2026rshift`, which is the right call |
| `spatstat2026rshift` | **CONFIRMED — quote verbatim** | help page re-rendered today; see §5 |
| `lipsitch2010negative` | **CONFIRMED — quote verbatim** | CR + PM (PMID 20335814) + PMC3053408. The erratum the note flags **exists**: Crossref confirms *Erratum: Negative Controls…*, Epidemiology 21(4):589, doi 10.1097/EDE.0b013e3181e4bfd7 (content still unread) |
| `moses2023voyager` | **CONFIRMED — both quotes verbatim; still a preprint** (but see **L1**) | CR posted-content only; bioRxiv API `"published": "NA"` for both v1 and v2; PMC10461913 → PMID 37645732. v2 (2023-08-20) author list `Moses, L.; Einarsson, P. H.; Jackson, K. C.; Luebbert, L.; Booeshaghi, A. S.; Antonsson, S. E.; Bray, N.; Melsted, P.; Pachter, L.` matches the entry's expanded forenames. 10.1101 prefix is correct for a 2023 preprint |
| `ren2025systematic` | **WRONG — see C1** | CR + PMC12534522. Record otherwise correct; both PMC quotes verbatim |

**Tally: 42 CONFIRMED, 1 WRONG, 0 unverifiable as records.** Two entries carry residual
content-level items that are correctly flagged in the file and remain open by design
(`ma2024spatial` figure panel, `acosta2026correction` full text) plus one that remains
genuinely unreadable (`lotwick1982methods` full text).

**On the forename question specifically:** I re-checked the given names of all 43 entries
against Crossref, and 25 of them against PubMed independently. **Every forename in the
file is now correct**, including all four in `kumar2026cellwhisper` and `Martin, Lucy`.
The earlier pass's arithmetic was wrong; its output was not.

---

# 4. Superseded-marker scope

Verified each marker against the text beneath it.

### Too narrow — forbidden or superseded content reads as live

| file | marker | what it misses |
|---|---|---|
| `Phase7_Minimal_Human_Replication (1).md` L3–27 | "SUPERSEDED IN PART", enumerating §10, §13 A6, §17, §23/C6, **§4 (D-b)**, §15 | **§4 D-d** (L187) carries "1.51–2.85×" — item 8. The banner names §4 only for D-b, so D-d reads live. Separately, the banner's own §4 (D-b) bullet asserts the superseded *"roughly doubles"* |
| `reports/NOVELTY_ASSESSMENT.md` L3–29 | "one recommendation in this report is FALSIFIED … **Everything else in this report stands.**" | It does not. §U1's drop-in paragraph (**L6** above) and §U2's Mrkvička misquote (**L5**) are both wrong. `CORRECTIONS.md` §17.4 says exactly this and the banner was never amended |
| `reports/CS_PHASE8_C1_CLOSEOUT.md` L38–49 | scoped to "**the surviving fractions** this report tabulates" | leaves live: "1 to 63 of 38,080 to 108,375" (frozen **1–66**), "median displacement 27 µm (N3)" (frozen **28 µm**), and "median λ̂ (12.8 µm)" ×3 (withdrawn for **14.7 µm** by `PREREG_PHASE8.md` C-7, which claims the withdrawal is "recorded project-wide") |
| `SASP_Kernel_Master_Plan.md` L184 | "This **paragraph** previously said 1,753 … as of January 2026" | §7 Rank 2 (L284) is the *other* location `BIO_DELIVERABLE7_CLAIM_AUDIT.md` names, and it still says 1,753 / January 2026 (**L8**) |
| `reports/BIO_DELIVERABLE6_DISCUSSION.md` L22 | "⚠ PRE-C6 DIGITS", scoped to the SF / amplitude / power vector | §A.1 item 2 (L72) *"agree at or below chance (Jaccard 0.93–1.22× of chance …)"* is a **claim**, not a digit, and is not covered |
| `README.md` L302–323 | "Superseded, and by how much. **This paragraph** read …" | the table header at L331, three paragraphs later, still says "callers that **agree at chance**" |

### Adequate

- `reports/CS_PHASE8_CALLERS.md` L3–33 — the strongest marker in the repo: names the file
  pre-C6 throughout, gives a pre-C6 → frozen substitution table including the A7 rows, kills
  the §3 drafted paragraph by name, and pre-empts two specific downstream misuses. Nothing
  beneath it reads live.
- `reports/PREREG_PHASE8.md` §0.0 C-1/C-4/C-5/C-6/C-7 — each names its target section and
  gives the frozen replacement inline; §10.1 and §6 R1 both carry the correction *in the
  paragraph*, which is the pattern that works.
- `reports/NOVELTY_ASSESSMENT.md` §4 status box (L556–577) — supersedes O1/O3/O4/O11/O12 and
  O5's λ̂ explicitly, so the "NOT DONE" cells beneath it and O1's now-dead
  orthogonality advice are covered. (Only §U1 and §U2, outside its scope, are not.)
- `reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` — the scope extension recorded at
  `RECORD_RECONCILIATION.md:339` is present and does cover the *Additionally* paragraph.

### Zero-marker documents carrying forbidden claims

`reports/BIO_PHASE3.md` (**L2**) and `reports/BIO_PHASE2.md` (**L3**). Neither has a
correction banner of any kind.

---

# 5. The spatial-statistics block — SUPPORTS lines, verified against the primary source

Every claim below was checked against a source retrieved today, not against recall.

**`mrkvicka2021revisiting` — all four SUPPORTS claims verbatim, and the section numbers
are right.** Extracted arXiv:1911.00240 to text and located each quotation:

1. §1, attribution + p.410: *"In point process literature the random shift approach was
   suggested already in Lotwick & Silverman (1982) who claim on p.410 that '… it seems
   intuitively that any statistics calculated after wrapping onto the torus will show less
   discrepancy from independence, and therefore that spurious significance should not be
   introduced.'"* — **exact**.
2. §2.1.1, liberality: *"Such a toroidal shift makes a crack in the autocorrelation
   structure of Ψ which causes the liberality of the procedure…"* — **exact**.
3. §2.1.4 "Shape of the observation window": the whole rectangular-window / union-of-aligned-
   rectangles / increased-liberality passage — **exact, word for word**, including
   *"The approaches using minus correction and variance correction can be applied in case of
   general (compact) observation windows."*
4. §2.1.3, variance correction: *"the values T₀,…,T_N are now computed from different amounts
   of data (windows W, W₁,…,W_N are different) and hence are not directly comparable"* —
   **exact**.

The section headings in the PDF are literally `1 Introduction`, `2.1.1 Torus correction`,
`2.1.3 Variance correction`, `2.1.4 Shape of the observation window`. The entry's CAUTION
about the `NOVELTY_ASSESSMENT` paraphrase is correct — and that paraphrase is **still live**
(**L5**).

**`lotwick1982methods`** — record confirmed; the entry's honest "NOT VERIFIED" on the full
text is accurate and should stay. The published abstract covers only second-moment
estimation and empty-space methods, so the entry's warning not to over-attribute from the
abstract is right.

**`spatstat2026rshift`** — I re-rendered
`search.r-project.org/CRAN/refmans/spatstat.random/html/rshift.ppp.html`. The sentence
*"The window must be a rectangle. Toroidal shifts are undefined if the window is
non-rectangular."* appears **verbatim**, in the `edge="torus"` description. The page has
**no References section** and the strings `Lotwick` and `Silverman` do not occur on it —
exactly as the entry's note says. Preferring this over `baddeley2015spatial` for the
constraint is the correct call.

**`dupont2022spatialplus`** — both quotes verbatim from the PMC10084199 abstract
(*"collinearity between covariates and spatial effects can lead to significant bias in
effect estimates"*; *"reduces the sensitivity of the estimates to smoothing by replacing the
covariates by their residuals after spatial dependence has been regressed away"*). The
SUPPORTS line's analogy to N6 is a fair reading, not a stretch.

**`hodges2010adding`** — see **C3**. Directly confirmed; the caveat is now over-cautious.

**`khan2022restricted`** and **`zimmerman2022deconfounding`** — both abstracts retrieved;
both SUPPORTS lines are supported, including Zimmerman's "inferior both to the original
spatial model and to the nonspatial model, and inferior predictive inference", which the
abstract states in those terms. The instruction to cite them together is sound: Khan &
Calder is the Bayesian/RSR-class result, Zimmerman & Ver Hoef the frequentist one.

**`baddeley2015spatial`** — record confirmed; the entry correctly declines to cite the book
for the constraint.

**`lipsitch2010negative`** — the quoted sentence (*"We distinguish 2 types of negative
controls (exposure controls and outcome controls) …"*) is verbatim in the PubMed abstract.
The SUPPORTS framing — A7 as a negative-control **outcome** for the estimand, as opposed to
Voyager/Ren's generic spatial statistic on controls — is the correct use of this paper. The
flagged erratum exists and is confirmed unread.

**`moses2023voyager`** — both Xenium-vignette quotations are **verbatim** in the live
vignette (*"As expected, generally the negative controls are tightly clustered around 0,
while the real genes have positive Moran's I, which means there is generally no technical
artifact spatial trend."* and the cell-density passage). Preprint status re-confirmed.
**The entry's closing instruction is a forbidden claim — see L1.**

**`ren2025systematic`** — both PMC quotations are **verbatim** in PMC12534522 §2 (*"Spatial
autocorrelation analysis using Moran's I revealed stronger aggregation of negative control
signals in CosMx 6K"*; *"CosMx 6K detected a higher total number of transcripts but exhibited
reduced spatial variation and elevated negative control signals compared to Xenium 5K"*).
The support is real; **the author line is wrong — see C1.**

**`kumar2026cellwhisper`** — I rendered the v1 full text
(`biorxiv.org/content/10.64898/2026.01.07.697982v1.full`) and confirmed all four audit
annotations:

- Null, verbatim: *"To assess significance, CellWHISPER shuffles locations among cells of
  same type (preserving cell-type expression and spatial distributions), recomputes N_ijkl
  and uses the resulting null to compute a z-score."* — this project's **N1**, not N3.
- The strings **`torus`**, **`toroidal`**, **`wraparound`/`wrap-around`**, **`blank barcode`**
  and **`negative control`** return **zero hits** in the full text. The mis-attribution fix
  is correct and the correct attribution is `lotwick1982methods` + `mrkvicka2021revisiting`.
- FPR, verbatim: *"CellChat v2, COMMOT, and SpaTalk predict similar numbers of interactions
  on real (blue) and randomized (orange) data, indicating poor specificity and false positive
  rates (FPR) >90%. In contrast, CellWHISPER produced markedly fewer interactions on
  randomized data compared to real data, suggesting FPR < 5%."* The entry's insistence on
  "implying", not "measured", is right — and the sentence also carries a scope caveat the
  entry does not: the benchmark ran on *"a downsampled hippocampal subset (~5,000 cells) due
  to computational constraints of competing methods"*. Worth adding if the >90 % figure is
  quoted.

**Verdict on Task 2: none of the eleven spatial-statistics entries is decorative.** Each
`SUPPORTS` line states a claim the source actually makes. The two defects in the block are
`ren2025systematic`'s author line (**C1**) and the forbidden sentence inside
`moses2023voyager`'s annotation (**L1**).

---

# 6. What I could not check

Stated separately from "searched and found nothing" and from "retrieved and it disagrees".

- **`lotwick1982methods` full text** — JSTOR/Wiley paywall; no open copy. The p.410
  quotation is verified only as Mrkvička quotes it. Unchanged from the entry's own note.
- **`acosta2026correction` content** — nature.com redirects anonymous fetches; Europe PMC
  holds metadata only. Existence, type (Published Erratum) and pagination confirmed; whether
  it touches the transwell result is **still unknown**.
- **`ma2024spatial` distance-gradient panel** — cell.com and ScienceDirect both 403; not in
  PMC or Europe PMC. The entry's "ONE HUMAN MUST OPEN THIS PDF" flag stands.
- **`lipsitch2010negative` erratum content** — record confirmed, text not retrieved.
- **`baddeley2015spatial` book text** — not consulted; not needed, since the constraint is
  cited to `spatstat2026rshift`.
- **Whether `results/` files back the digits quoted in the documents** — out of scope here;
  I verified *wording* against the prohibitions lists, not numbers against CSVs. Where I
  name a frozen value it is taken from `PREREG_PHASE8.md` §0.0 or `WRITING_PACK.md` §3,
  which cite the files.
- **`data/raw_h1/`, `results/phase9_h1/`, `data/processed_h1/`, `code/h1_*`** — not read, per
  instruction; the human arm is another agent's.

---

# 7. Recommended order of fixes

1. **`references.bib:783`** — delete the "so the reader can see the two tests disagree"
   clause from the `moses2023voyager` SUPPORTS block and replace it with the power
   framing. The bibliography is the one file that gets copied into the manuscript wholesale.
2. **`reports/BIO_PHASE3.md`** — add a banner in the `CS_PHASE8_CALLERS.md` style, and
   strike §4.5's blockquote by name. It is the only zero-marker document with a
   drafted-for-the-paper paragraph made of struck sentences.
3. **`reports/BIO_PHASE2.md` §4.3** — banner.
4. **`reports/CS_PHASE8_M1_RERUN.md:307`** — rename the row to "pooled negative-control
   features", state the 40 probes as flat, and drop "NOT FLAT" as the unqualified verdict.
5. **`references.bib`** — fix `ren2025systematic`'s last author; delete the stale
   `hughes2025senpred` title note; optionally upgrade the `hodges2010adding` and
   `zimmerman2022deconfounding` verification lines.
6. **`reports/NOVELTY_ASSESSMENT.md`** — amend the head banner (it currently says
   "Everything else in this report stands"), fix §U2's misquote, restate §U1's paragraph.
7. Widen the four narrow markers listed in §4; correct `SASP_Kernel_Master_Plan.md:284`;
   standardise the remaining `2.4×` to `2.35×` in `WRITING_PACK` and `PLAN_UPDATE`.
