# BIO PHASE 1 — Dataset Identity, Gene Sets, Zonation

**Biology collaborator · 2026-08-20 · SASP Spatial Response Kernel**
Scope: Master plan §7 (data), §8 (Day-1 audit), §9 (gene sets), §11 (zonation), §12 (Deliverables 1–2).

---

## 0. Headline

1. **The plan's Section 7 Rank 1 entry is wrong about the accession, right about the paper.**
   GSE310392 is the **mouse** arm of Karpova et al. It is not, and never was, the 43-donor human
   atlas. **Verified**, not inferred — from the paper's own data-availability statement.
2. **Section 8 Test 2 PASSES on all three criteria** against the real mouse 5K panel — but only
   after a design change to how disjointness is enforced. The literal rule passes numerically while
   destroying the sender score.
3. **The single biggest threat is no longer panel adequacy. It is replication units: n = 11 mice,
   one section each, versus the 43 donors the plan budgeted for.**

---

## 1. Section 8 audit summary table

| Test | Metric | Value | Pass? |
|---|---|---|---|
| 1 | Median NN cell distance (µm) | **6.74** (sbr) / **9.66** (sham) | ✅ resolution floor recorded |
| 1 | Coordinate units | **microns**, `x_centroid`/`y_centroid` | ✅ confirmed |
| 1 | Median transcripts/cell | 446 (sbr) / 511 (sham) | ✅ |
| 1 | Median genes/cell | **288** (sham) | ✅ |
| 1 | Cells / features | 128,030 / 237,982 · **5,106 genes** | ✅ |
| 1 | Transcript **assignment rate** | **NOT MEASURED** — needs `transcripts.parquet` (1.85 GB/sample), not downloaded | ⚠️ **unmeasured, not passed** |
| 2 | Sender genes on panel | **25** union-strict; **37–74** per-module | ✅ ≥15 |
| 2 | Response genes on panel (per module) | 126 / 68 / 100 / 190 / 125 / **31** / 38 | ✅ ≥30 (B6 only after documented substitution) |
| 2 | Sender ∩ response | **0** | ✅ must be 0 |
| 3 | Sender prevalence (`Cdkn1a` > 0) | **4.22%** sbr / **0.68%** sham | ⚠️ sbr in band; **sham below the 1% floor** |
| 3 | Sender prevalence (Tier A score, fixed threshold) | p95: 8.46% sbr / 5.00% sham | ✅ in 1–20% |
| 3 | Senders per donor | ≫200 both sections | ✅ |
| 4 | Ripley's K vs Poisson | not run (CS lead) | — |
| 5 | Max SMD after decoy matching | not run (CS lead) | — |
| 6 | Pseudo-R² anatomy → sender | not run; **zonation axis built and validated** (§5) | — |

Test 1 numbers are the CS lead's; reproduced here for the Methods table. Tests 4–6 are CS-lead
tasks and are not blocked by anything in this report.

---

## 2. WHAT PASSED

### 2.1 Dataset identity is now fully resolved (Deliverable 1)

**GSE310392 = "Cellular senescence in human liver under normal aging and cancer [Xenium]",
organism *Mus musculus*, 12 samples, GPL33896, Washington University (Ding lab),
submitted 2025-11-18, updated 2026-04-22, PMID 41576948.**

The GEO **title** says "human liver" because it inherits the paper title; the **organism field,
every sample, and every gene ID are mouse**. Both facts are on the same page. That is the whole
source of the confusion in the plan.

Paper: Karpova A, …, Ding L. *Cell Genomics* 2026; 6(2):101133. doi:10.1016/j.xgen.2025.101133.
Its data-availability statement splits the deposit:

| Data | Where | Access |
|---|---|---|
| **Human** (43 livers; 25 Xenium sections) | **SenNet portal**, provider group = Washington University | not GEO |
| mCRC (24 metastases) | HTAN DCC, HTAN WUSTL Atlas | not GEO |
| **Mouse Xenium** | **GEO GSE310392** ← what we have | public |
| Mouse snRNA-seq | GEO GSE311064, GSE293958 | public |

**The model.** 75% proximal small-bowel resection (SBR) in male C57BL/6J mice at 8 weeks, versus
sham surgery; a model of **IFALD (intestinal-failure-associated liver disease)**, confirmed by the
Zenodo deposit's own directory name `IFALD_mouse`. It is **not** a CCl4 fibrosis model — do not
describe it as one.

**Full sample inventory (all 12, fetched from GEO):**

| GSM | Title | Tissue | Arm | Timepoint |
|---|---|---|---|---|
| GSM9295276 | 7361_liver_sbr_Male_2-U1 | liver | sbr | 2 wk |
| GSM9295277 | 7352_liver_sham_Male_2-U1 | liver | sham | 2 wk |
| GSM9295278 | 7448_liver_sbr_Male_10-U1 | liver | sbr | 10 wk |
| GSM9295279 | 7450_liver_sbr_Male_10-U1 | liver | sbr | 10 wk |
| GSM9295280 | 7435_liver_sham_Male_10-U1 | liver | sham | 10 wk |
| **GSM9295281** | **7259_liver_sbr_Male_26-U1** | liver | sbr | 26 wk | ← on disk |
| GSM9295282 | 7260_liver_sbr_Male_26-U1 | liver | sbr | 26 wk |
| GSM9295283 | 7248_liver_sham_Male_26-U1 | liver | sham | 26 wk |
| **GSM9295284** | **7250_liver_sham_Male_26-U1** | liver | sham | 26 wk | ← on disk |
| GSM9295285 | 7239_liver_sbr_Male_52-U1 | liver | sbr | 52 wk |
| GSM9295286 | 7239_tumor_sbr_Male_52-U1 | **tumor** | sbr | 52 wk |
| GSM9295287 | 7001_liver_sham_Male_52-U1 | liver | sham | 52 wk |

**11 distinct mice, 12 sections** (mouse 7239 contributes both a liver and a tumor section).
6 sbr liver, 5 sham liver, 1 tumor. Per timepoint: 2 wk (1v1), 10 wk (2v1), 26 wk (2v2), 52 wk (1v1).

**Can it serve the Rank 1 role? Partly — with three explicit downgrades.**

| Plan assumed | Reality | Consequence |
|---|---|---|
| 43 human donors → donor bootstrap (§24.1) | **11 mice, 1 section each** | §24.1 donor-level CIs rest on n=11, or **n=4 within the 26 wk timepoint**. This is the binding constraint on the whole paper. |
| Human **aging** axis | **IFALD injury-duration axis**, 2/10/26/52 wk post-surgery | Substitutable but *not the same claim*. It is disease progression, not chronological aging. Must be reworded everywhere. |
| Senescence **already annotated** (CDKN1A+ hepatocytes etc.) | **No annotations for mouse Xenium.** The Zenodo deposit contains `human_Xenium` merge scripts and mouse snRNA scripts but **zero mouse Xenium code** (verified by exhaustive grep). | We must call senescence ourselves (DeepScence/SenePy, §10). Adds Days 2–4 work; does not block. |
| Human 5K panel | **Mouse 5K panel + 100 custom** = 5,106 genes | Fine — see §2.2. Requires mouse orthologs throughout. |

The paper's periportal-CDKN1A finding (§11) was reported in **human**; whether it holds in this
mouse model is an open question we can test, not an assumption we can import.

### 2.2 Panel adequacy — Test 2 passes on all three criteria

Panel verified feature-by-feature: the union of the two panel files reproduces the h5
`Gene Expression` feature list **exactly**, 5,106 genes, 0 discrepancies either way. Both panels
apply to **both** our samples (the custom 100 adds 91 genes + 9 genotyping probes; **zero overlap**
with the 5K panel).

| Criterion | Required | Actual | Verdict |
|---|---|---|---|
| `len(A)` | ≥15 | **25** union-strict / 37–74 per-module | **PASS** |
| `len(B)` per module | ≥30 | 126, 68, 100, 190, 125, **31**, 38 | **PASS** |
| `len(A∩B)` | 0 | **0** | **PASS** |

Full intersection matrix, per-gene collision lists, and provenance: `/workspace/genesets/README.md`.

Gene sets are **real MSigDB mouse collections, release 2026.1.Mm, fetched live 2026-08-20** from the
MSigDB JSON API, with MSigDB's own Alliance-Genome ortholog mapping. **Nothing is hand-mapped and
nothing is recalled from memory** — Section 9's verification requirement is discharged. Raw JSON
archived in `genesets/msigdb_mouse_2026.1.Mm/`.

### 2.3 Tier C internal control is intact

`Ccl2`→`Ccr2`+`Ackr3`, `Cxcl12`→`Cxcr4`+`Ackr3`+`Dpp4`, `Il1a`→`Il1r1`+`Il1rap`,
`Tnf`→`Tnfrsf1a/1b`, `Il6`→`Il6ra`+`Il6st`, `Tgfb1`→`Tgfbr1/2` — **all on-panel**. The Section 9
λ-ordering control (membrane-bound Il1a shortest, diffusible chemokines longest) is runnable.

### 2.4 A validated zonation covariate now exists

Section 11 calls zonation the clearest liver confound. The plan's marker list is human and mostly
off-panel (`Alb`, `Ass1`, `Sds`, `Oat` all absent). Rebuilt from the real sham matrix and validated:
pericentral score r = **+0.897**, periportal r = **−0.737** against the axis, versus a null of 300
random size-60 sets (mean |r| = 0.114, p90 = 0.261). Method, caveats and the 42 retained markers:
`genesets/README.md` §6; per-gene table `/workspace/results/zonation_gene_correlations_7250_sham.csv`.

### 2.5 Good news on the zonation confound — measured, not assumed

Section 11 predicts response modules will covary with zonation and manufacture fake decay curves.
Measured directly, as correlation of each **module score** with the zonation axis within
hepatocytes (n=142,143, sham):

| Set | r with zonation | Verdict |
|---|---|---|
| A_SENDER_FINAL_strict | **+0.058** | clean |
| B1 tnfa_nfkb | −0.045 | clean |
| B2 il6_jak_stat3 | −0.041 | clean |
| B3 interferon | −0.217 | inside null p90 (0.261) |
| B4 downstream_arrest | −0.080 | clean |
| B5 emt_ecm | −0.139 | clean |
| B6 oxidative_stress | −0.025 | clean |
| B7 secondary_senescence | −0.051 | clean |
| D_zonation_pericentral | **+0.897** | (positive control) |
| **E_hepatocyte_identity** | **−0.666** | ⚠️ see §3.3 |

**No Tier B module carries an unusual intrinsic zonation loading.** Section 11's fear is not
realised at the module level. Note the limit: this measures *intrinsic* loading only. Zonation can
still bite through the sender-clustering path (if senders sit periportally), which is Test 6 and
still must be run.

### 2.6 The sbr-vs-sham dose-response positive control survives

The CS lead flagged that burden and packing density covary, weakening this control. Measured:

| | sbr (7259) | sham (7250) | ratio |
|---|---|---|---|
| `Cdkn1a`-positive cells | **4.22%** | **0.68%** | **6.2×** |
| Tier A score, fixed threshold p95 | 8.46% | 5.00% | 1.69× |
| median NN distance | 6.74 µm | 9.66 µm | 1.43× (density) |

**The `Cdkn1a` burden contrast (6.2×) is four times larger than the packing-density contrast
(1.43×).** The positive control survives on the single-gene marker. It is much weaker on the
composite Tier A score (1.69×, same order as density) — so **use `Cdkn1a`+ prevalence, not the Tier
A composite, for the burden contrast**, and condition on local density regardless. n=1 vs n=1 here;
confirm on GSM9295282/GSM9295283 before relying on it.

---

## 3. WHAT FAILED

### 3.1 The literal Section 8 disjointness rule produces a hollow sender score — **design flaw in the plan**

Enforcing "A disjoint from ∪B" removes **49 of 74** Tier A candidates. The surviving 25 pass the
≥15 bar but contain **no `Cdkn1a`, no `Cdkn2a`, no `Trp53`, no `Lmnb1`, no `Mki67`** and still
contain `Ccnb1`/`Foxm1`, whose expected direction is opposite. It is a numerically passing,
biologically meaningless sender definition.

Root cause, **verified against MSigDB 2026.1.Mm**: `Cdkn1a` is a member of
`HALLMARK_TNFA_SIGNALING_VIA_NFKB`, `HALLMARK_INTERFERON_GAMMA_RESPONSE` and the E2F/G2M/MYC union.
And Section 9's own B4 (`E2F_TARGETS`/`G2M_CHECKPOINT`) **is** the arrest program that Tier A2
measures — the plan asks for two disjoint sets that are definitionally the same biology.

**Resolution (adopted).** Disjointness is only required between the sender score and the **one**
response module being fitted. Per-module sender sets `A_sender_for_<module>.txt` are all ≥37 genes,
all disjoint from their own readout, and 5 of 7 retain the canonical markers. Use these as primary;
`A_SENDER_FINAL_strict` as the conservative sensitivity analysis. **This needs the CS lead's sign-off
because it is a deviation from §8 as written, and the intersection matrix in Methods must show both.**

This is also a *result*, not just an inconvenience: it is a concrete, quantified instance of the
Section 0.3 circularity that the paper exists to expose, and it belongs in the paper.

### 3.2 B6 oxidative stress fails as specified

`HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` gives **17** on-panel genes — **FAIL** against ≥30.
Fixed by a documented union with `REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES` + curated NRF2
targets → **31**. It is the weakest module and should be reported last, with its substitution stated.

### 3.3 Tier E negative controls fail as specified — **both of them**

- **Housekeeping: 1/6 on-panel.** Only `Tbp`. `Actb`, `Gapdh`, `Rpl13a`, `Rps18`, `Ppia` are all
  absent — Xenium omits high-expressors to avoid optical crowding. Expanded to 13 on-panel
  stably-expressed genes; usable but thin.
- **Cell-type identity control is worse than absent, it is actively misleading.** Section 9 says
  `ALB`/`TTR`/`APOA1` "should be flat after conditioning on cell type". Measured: its score
  correlates with zonation at **r = −0.666**, worse than 99% of random size-matched panel sets.
  In liver, the hepatocyte identity program *is* a zonation readout. **A flat-kernel test against
  this set would be a false alarm generator.** Do not use it as a negative control; it is a
  *positive* control for zonation. (`Alb` itself is absent from the panel.)
- **Substitute (recommended deviation):** per-cell background columns already in `cells.parquet`
  (`control_probe_counts`, `genomic_control_counts`, `control_codeword_counts`,
  `unassigned_codeword_counts`). They measure decoding background directly. **Constraint:** ~0.17
  counts/cell, so they work **binned** (≈50 µm) as a spatial-gradient control, not per cell. Use
  the 13-gene housekeeping set for the per-cell flat-kernel test. Report both, and state that
  Section 9's Tier E could not be run as written.

### 3.4 Transcript assignment rate — unmeasured

§8 Test 1 asks what fraction of transcripts are assigned to cells (>30% unassigned ⇒ bleed-through,
which manufactures spatial autocorrelation). It needs `transcripts.parquet` (1.85 GB/sample), not
downloaded. The control/unassigned-**codeword** fractions (0.1–0.2%) are a different quantity —
probe background, not assignment rate. **My judgement: pull it for one sample.** Bleed-through is
the exact failure mode this paper claims to police; leaving the number blank while asserting a
controlled λ is the kind of gap a reviewer will find. One sample is enough and fits the quota.

---

## 4. WHAT IS UNCERTAIN

1. **Zonation axis from one section.** Derived from sham only, using a marker-score gate rather
   than real cell-type calls. Recompute per section; check marker stability. `Pigr` (also a
   cholangiocyte marker) loads strongly periportal and may signal bile-duct contamination of the gate.
2. **The 6.2× `Cdkn1a` burden ratio is n=1 vs n=1.** Confirm on GSM9295282/83 before building a
   figure on it.
3. **`Cdkn1a` detection in sham is 0.68%** — below Test 3's 1% floor. Senescence calling in the
   *control* arm may be underpowered, which matters because sham is the reference for every contrast.
4. **B7 is not independent of Tier A by construction** (Section 9 defines it that way). Needs the
   split-half `A_call`/`A_readout` design with a fixed seed. Not yet built.
5. **Mouse-vs-human transfer.** All of §9/§11 was written for human. Ortholog mapping is handled,
   but the *biology* claims (CDKN1A+ periportal localisation) were human findings.
6. **`Il1b` off-panel** weakens the sharpest form of the Tier C internal control (membrane-bound
   Il1a vs secreted Il1b through the same receptor).
7. **SenMayo audit.** Mouse SAUL_SEN_MAYO is 117 genes, 90 on-panel, of which **48 (53%)** are
   annotated `Secreted` by the panel. Directionally consistent with the plan's "83 of 125 are SASP"
   but **not equal to it** — the plan's figure is a human-set claim and our 53% is a
   panel-annotation proxy, not the same measurement. Do not cite one as confirming the other.

---

## 5. WHAT THE CS LEAD MUST KNOW

1. **Rewrite Section 7 Rank 1.** It is a mouse IFALD/SBR surgical model, n=11 mice, injury-duration
   axis, no senescence annotations — not a 43-donor human aging atlas. Every downstream sentence
   about "43 donors", "aging", and "senescence already annotated" needs changing. The human data is
   on the **SenNet portal** (WashU provider group), not GEO, if we want it.
2. **Sign off on the per-module sender-set deviation (§3.1).** Blocks the kernel fit. Both the
   union-strict and per-module intersection matrices must appear in Methods.
3. **§24.1 donor bootstrap: n = 11 mice, or n = 4 at the 26 wk timepoint.** Not 43. This is the
   binding constraint on the paper's central claim (uncertainty). Decide now whether to pool
   timepoints (confounds progression) or accept n=4.
4. **Tier E must be rebuilt (§3.3)** — both Section 9 controls fail, and the cell-type-identity one
   fails in a direction that would generate false positives.
5. **Pull `transcripts.parquet` for one sample** (§3.4).
6. Add `segmentation_method` (3 levels, spatially patterned) to Tier D nuisance covariates.
7. Density confound is real: matched-decoy selection **within sample only**; consider reporting
   distances in units of local median NN distance as a sensitivity analysis.
8. Drop the 9 genotyping probes before scoring.
9. Resolution floor is **6.7–9.7 µm**; do not report λ below it. Expected λ is tens of µm, so there
   is headroom.

### Top 3 threats to the project

1. **Replication units.** n=11 mice / n=4 at 26 wk, against a plan built on 43 donors. The paper's
   whole selling point is uncertainty done properly, and the uncertainty is now dominated by a very
   small number of animals. Consider adding the SenNet human data or a 10x public Xenium cohort as
   the replication cohort (§7 Rank 4).
2. **Circularity is structural, not incidental.** `Cdkn1a` sits inside the canonical inflammatory
   Hallmark sets. Any naive sender/response pairing in this field is measurably circular. We can
   handle it — and should *report* it — but it means the headline λ depends on a gene-set design
   choice that must be pre-registered and defended.
3. **The claim has changed under us.** "How far does senescence reach in aging human liver" has
   become "…in a mouse model of surgical intestinal-failure liver disease." That is a defensible
   paper, but it is a different one, and the framing, title and abstract must change before writing
   starts — not after.

---

## 6. Deliverables produced

| Path | Contents |
|---|---|
| `/workspace/genesets/*.txt` | Tiers A–E + zonation, one mouse symbol per line |
| `/workspace/genesets/README.md` | Provenance, versions, dates, full intersection matrix, deviations |
| `/workspace/genesets/msigdb_mouse_2026.1.Mm/` | Archived MSigDB JSON (offline reproducibility) |
| `/workspace/code/build_genesets.py` | Rebuilds everything; no network needed |
| `/workspace/results/zonation_gene_correlations_7250_sham.csv` | Per-gene zonation correlation, all 5,106 |
| `/workspace/reports/BIO_PHASE1.md` | This report |

Still owed (Deliverables 3–4): per-cell anatomical covariate table, sender-call validation.
