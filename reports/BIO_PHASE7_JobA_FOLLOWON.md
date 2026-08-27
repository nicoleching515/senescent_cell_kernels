# BIO PHASE 7 / Job A follow-on — Spleen Markers, the A6 Covariate, C6, Tier E2, SenePy, and the Tier A Fork

**Biology collaborator · 2026-08-27 · SASP Spatial Response Kernel**
Companion to `reports/BIO_PHASE7_JobA.md`. Six assigned tasks, all against H1 = **GSE326743,
7 normal human spleens**, Xenium Prime 5K Human + 100 addon, 5,093 `Gene Expression` features.

**Freeze discipline.** Nothing in this block read any expression value, cell record or annotation
from `/workspace/data/raw_h1/`. The only H1-derived input is panel membership
(`genesets/h1_candidate/`), the sanctioned §12.1 screen. Two of the six tasks are therefore
delivered as *specifications with the validation step named and deferred*, and said so plainly
rather than being filled in. Nothing under `code/run_deepscence*`, `data/processed/deepscence_*`,
`results/phase3/` or `code/make_figure4.py` was written to.

---

> **PARTLY SUPERSEDED — corrections, 2026-08-27.** The PI adopted the re-sourced B7 recommended in
> §3 below and froze Tier A with strict as primary. Current values are in
> `reports/PREREG_PHASE8_genesets.md` and `reports/BIO_PHASE8_FREEZE.md`. What moved:
> `B_secondary_senescence` is now **116** genes (human) and **108** (mouse, `genesets/mouse_c6/`);
> `A_SENDER_FINAL_strict` is **33** on both arms; the CoreScence circularity figure moved
> **76% → 88%** (cite 88%). §6's Tier A memo is **resolved** — strict is primary, per-module is the
> pre-registered sensitivity. §1, §2, §4 and §5 stand as written.

## 0. Headline

1. **A spleen marker set exists and is not curated from memory.** 22 assignable cell types, every
   gene traced to PMIDs in a pinned CellMarker 2.0 download. §1.
2. **Three cell types the panel cannot resolve**, named: **follicular dendritic cells** (3 on-panel
   markers), **plasma cells** (3, and it fails by exactly one gene), **proliferating cells** (3,
   once the over-adjustment guard removes `MKI67`/`CDK1`/`HMGB2`). §1.3.
3. **C6's literal prescription breaks B7.** Subtracting the Tier A caller from the curated B7
   leaves **24 genes** (§10 caller) or **12** (ported caller) — both **fail the ≥30 floor**.
   Re-sourcing B7 from SenMayo ∪ Reactome SASP and *then* subtracting gives **119 / 116**, passes,
   is disjoint by construction, and replaces the meeting-abstract citation with PMID 35974106.
   Under the ported Tier A it also raises the strict sender set from **21 to 33**. §3.
4. **SenePy ships no spleen signature at all.** 65 human hubs across 10 tissues; spleen is not one
   of them. Every hub usable on H1 is a **cross-tissue surrogate**, where the mouse arm used
   tissue-matched Liver hubs. 15 of 22 labels get a surrogate; **7 get nothing**. §5.
5. **Tier E2 should be dropped for this arm, on the record.** Every candidate identity program is,
   by construction, a member of one A6 compartment — the exact failure that killed the mouse arm's
   E2 (r = −0.67 with the zonation axis). §4.
6. **The A6 covariate is specified and ready to pre-register**, and it is the cleanest structural
   match either arm has had. It cannot be built or validated pre-freeze. §2.
7. **The §14 "eyeball fifty cells against the image" check is not available** — morphology images
   were deliberately not downloaded. A replacement is proposed. §1.4.

Two new pinned sources: **CellMarker 2.0 human** (`genesets/cellmarker_pin/`, md5
`c7a1b764b66cb3a3c16cfac428160f72`) and, from the previous block, HGNC. **`openpyxl 3.1.5` was
installed with pip** to read the CellMarker `.xlsx`; it is not in `requirements.txt`. Reported
because environment pinning is a tracked concern. The build itself reads the pinned `.csv.gz`, so
nothing downstream depends on `openpyxl`.

---

## 1. Task 1 — spleen cell-type markers

```bash
python3 /workspace/code/build_markers_human_spleen.py   # -> code/markers_human_spleen.py
```

`code/markers_human_spleen.py` is a **generated** file in the exact format of
`code/markers_mouse_liver.py` — a `MARKERS` dict of `'Label':'G1 G2 ...'.split()`, plus `MERGE`
and `DROP_NONSPECIFIC` in the shape `annotate_pipeline.py` already consumes. Nothing in the
pipeline needs changing.

### 1.1 How the genes were chosen

Source: **CellMarker 2.0 human**, 60,877 rows, `tissue_class == 'Spleen'` = 447 rows over 62
`cell_name` values. Five mechanical filters, in order, each printed in
`results/phase7_jobA/build_markers_human_spleen.log`:

| Filter | Rule |
|---|---|
| Evidence ladder | spleen rows (≥1 PMID) → all-tissue (≥2 PMID) → all-tissue (≥1 PMID). First tier that reaches 4 on-panel genes wins; the tier used is recorded per label. |
| Panel | must resolve onto the 5,093 genes via `code/human_symbols.py` (so panel legacy symbols are matched, not missed). |
| Specificity | a gene claimed by more than 3 labels is dropped as promiscuous. Removed: `CD14 CCR7 SELL ACTA2 CD34 PECAM1`. |
| Arbitration | CellMarker rows are per-publication marker *panels*, not exclusive markers. Within a compartment a gene is kept only by the best-evidenced label — **and a label on a weaker evidence tier can never take a gene from a label with real spleen evidence**. Without the tier rule, Neutrophils (5 spleen rows, 49 all-tissue candidates) strips `CD14`/`FCAR`/`S100A12` from Monocytes and deletes that label. |
| Over-adjustment guard | Tier A and Tier C genes are removed from every marker set. |

Arbitration did the work it was meant to: `CD8A`→CD8 T, `CD4`→CD4 T, `KLRD1`→NK, `PROX1`→
Endothelial, `RGS5`/`PDGFRB`/`MCAM`→Pericytes, `ACTA2`/`MYH11`→Smooth muscle. Full ledger in the
log.

**The over-adjustment guard is a deliberate deviation from the mouse arm.** The cell-type call is
a Tier D nuisance covariate; building it from the genes whose spatial behaviour is the outcome
conditions the response on itself. The mouse arm applied this guard to the zonation covariate
(`genesets/README.md` §6) but **not** to its cell-type markers — `markers_mouse_liver.py` calls
`'Proliferating'` with `Mki67 Top2a Pcna Ccnb1 Ccna2 Birc5 Aurkb Cdk1`, all of which are Tier A or
Tier B. Applying it here is a correction, and it is why `Proliferating cells` drops out (§1.3).
It is applied at Tier A ∪ Tier C only: extending it to Tier B removes 58 marker slots and destroys
six labels including **both stromal labels the A6 covariate depends on**. Tier B overlaps are
therefore flagged in `genesets/human/markers_spleen_evidence.csv` (`tierB_member`), not removed.

### 1.2 On-panel coverage per cell type

22 assignable labels. `scope` records which evidence tier the label reached.

| Cell type | markers | evidence scope |
|---|---|---|
| Red pulp macrophages | 15 | spleen-only (37 rows) |
| Monocytes | 5 | spleen-only (19) |
| cDC1 | 4 | spleen-only (10) |
| cDC2 | 5 | spleen-only (17) |
| pDC | 8 | spleen-only (11) |
| Follicular B cells | 9 | spleen-only (22) |
| Marginal zone B cells | 8 | **all-tissue, ≥1 PMID — weakest** (6 spleen rows) |
| Germinal centre B cells | 4 | **all-tissue, ≥2 PMID** (0 spleen rows) |
| CD4 T cells | 14 | spleen-only (37) |
| CD8 T cells | 15 | spleen-only (42) |
| NK cells | 15 | spleen-only (37) |
| Fibroblastic reticular cells | 4 | **all-tissue, ≥1 PMID — weakest** (0 spleen rows) |
| Sinusoidal endothelium | 9 | **all-tissue, ≥1 PMID — weakest** (0 spleen rows) |
| Endothelial cells | 4 | spleen-only (11) |
| Lymphatic endothelium | 4 | all-tissue, ≥2 PMID (0) |
| Smooth muscle / capsule | 6 | all-tissue, ≥2 PMID (0) |
| Pericytes | 12 | all-tissue, ≥2 PMID (0) |
| Fibroblasts | 5 | spleen-only (13) |
| Erythroid cells | 5 | all-tissue, ≥2 PMID (10 spleen rows) |
| Megakaryocytes | 15 | **all-tissue, ≥1 PMID — weakest** (10 spleen rows) |
| Neutrophils | 15 | all-tissue, ≥2 PMID (5 spleen rows) |
| Mesothelial cells | 15 | **all-tissue, ≥1 PMID — weakest** (0 spleen rows) |

**Read the scope column before trusting a fine label.** Nine of 22 rest on markers validated in
other organs, and four of those on the weakest tier. `Marginal zone B cells` — the compartment
that most distinguishes spleen from every other lymphoid organ — has only **6** spleen rows in
CellMarker 2.0 and had to fall back. Its surviving set (`CD180 CR2 EBF1 FCRL4 ITGAE NOTCH2 PAX5
TNFRSF13C`) is plausible but is not spleen-validated evidence.

`MERGE` groups (what `cell_type_merged` will be computed over): B cells, T/NK cells, Mono/Mac/DC,
Endothelial, Stromal.

### 1.3 What the panel cannot resolve — named

| Cell type | on-panel markers surviving | why it matters |
|---|---|---|
| **Follicular dendritic cells** | 3 — `CR1 CR2 FCER2` | FDCs *are* the follicle's stromal scaffold. Their loss means the B-follicle compartment in §2 is defined by B cells alone, with no independent stromal confirmation. CellMarker has only 3 FDC rows in spleen and 9 in total. |
| **Plasma cells** | 3 — `JCHAIN MZB1 XBP1` | **Fails by exactly one gene.** These three are highly specific. `annotate_pipeline.py`'s `MIN_MARKERS = 4` is what drops it. A PI decision: lower the threshold for this label, or let plasma cells fold into the B compartment. Flagging rather than quietly re-tuning the threshold. |
| **Proliferating cells** | 3 — `CD3D IL7R TUBB` | Only after the guard removes `MKI67 CDK1 HMGB2`. Keeping it would mean labelling cells "proliferating" with the same genes that define Tier A and Tier B — the mouse arm did exactly this. Its residual markers are T-cell genes, so it was not a clean label anyway. |

Two more labels ship with visible CellMarker contamination that arbitration could not fix because
they have no compartment partner on the same evidence tier: **Megakaryocytes** retains
`CD79A CPA3 CCL5 CXCR4`, and **Mesothelial cells** retains `CENPF MYH11`. Both are weakest-tier
labels with zero spleen rows. Treat them as provisional; they are the first two to re-gate on
measured expression.

### 1.4 The §14 "eyeball fifty cells against the image" check

**Not available as written.** §12.2's skip list excluded `morphology.ome.tif.gz`, so there is no
image to eyeball against — confirmed in `reports/PHASE7_H1_SCREEN.md` §4 (21 files skipped).

Proposed replacement, in priority order:

1. **Check against the depositors' own labels.** `annotations.csv.gz` ships four nested annotation
   levels. This is a *stronger* check than eyeballing fifty cells: it is an external, whole-section
   label set produced by people who had the images. Compute the confusion matrix of
   `cell_type_merged` against `Level_1`/`Level_2`, report adjusted Rand and per-label recall.
   **This is post-freeze work** and its caveat is already recorded: their labels cover fewer cells
   than the matrix (SPLN07: 239,167 of 249,420), so their QC filter must be characterised before
   the two label sets are compared, or the disagreement will be read as annotation error when it
   is really a filtering difference.
2. **Spatial-coherence sanity check, image-free.** Splenic architecture is macroscopic: B follicles
   are compact 200–500 µm islands, the T zone wraps the central arteriole, red pulp is the
   continuous background. A correct annotation must show follicular B cells forming compact
   spatial clusters, not a salt-and-pepper scatter. Quantify with Ripley's K per label using the
   existing `code/_ripley.py` — machinery the project already has, and a check that catches the
   failure mode eyeballing is meant to catch.
3. If the PI wants the literal check, `morphology.ome.tif.gz` can still be downloaded for **one**
   section. It is a §12.2 deviation and should be declared as one.

---

## 2. Task 2 — the A6 anatomical covariate for spleen (specification)

**Status: specified, not built.** The mouse zonation covariate was derived from the **expression
matrix** (`genesets/README.md` §6: seed axis, per-gene correlation, |r| ≥ 0.20 retention). H1
expression is behind the §15 freeze, so the equivalent derivation cannot run and **no part of this
covariate has been validated**. What follows is precise enough to pre-register and to implement
without further biological judgement.

### 2.1 Compartment marker sets — built and on disk

```bash
python3 /workspace/code/spec_a6_compartments_human.py
```

Written to `genesets/human/`, all carrying the Tier A + Tier C guard:

| File | n | Defined by |
|---|---|---|
| `D_spleen_white_pulp_tzone` | 24 | CD4 T, CD8 T, fibroblastic reticular cells |
| `D_spleen_white_pulp_follicle` | 13 | follicular B, germinal centre B |
| `D_spleen_marginal_zone` | 8 | marginal zone B |
| `D_spleen_red_pulp` | 28 | red pulp macrophages, sinusoidal endothelium, erythroid |
| `D_spleen_capsule_trabecula` | 18 | smooth muscle / capsule, pericytes |

**The five compartments are pairwise disjoint — every off-diagonal entry of the compartment ×
compartment matrix is 0** (`spec_a6_compartments_human.log`). Tier A and Tier C overlap is 0 by
construction; Tier B overlap is 6 / 3 / 1 / 3 / 6 genes respectively, flagged not removed.

### 2.2 The covariate itself

Liver zonation is a **continuous** axis and was carried as a continuous covariate. The spleen
axis is **not** continuous end-to-end: white pulp is a set of discrete islands embedded in a
continuous red pulp background. So it takes three parts, not one, all computed from
`cells.parquet` `x_centroid`/`y_centroid` (microns, verified) plus the Job B cell-type calls:

**(a) `compartment` — a categorical label, per cell.** Five levels: `white_pulp_follicle`,
`white_pulp_tzone`, `marginal_zone`, `red_pulp`, `capsule_trabecula`. Assignment is *spatial*,
not per-cell transcriptional, because a single T cell in the red pulp is not the T zone:

1. Score each cell on all five compartment marker sets (`sc.tl.score_genes`, matched control sets
   — the same call `annotate_pipeline.py` already makes).
2. Smooth each of the five scores over a **50 µm** radius neighbourhood (the mouse arm's own
   density radius set; 25/50/100 µm are already Tier D covariates).
3. Assign the cell to the argmax **smoothed** compartment score, with the same margin rule the
   pipeline uses elsewhere (`MIN_MARGIN = 0.20`) and `compartment = 'ambiguous'` otherwise.
4. **Fixed follow-up, no judgement:** connected components of `white_pulp_follicle` ∪
   `white_pulp_tzone` cells at 50 µm linkage form **white pulp nodules**. Discard components with
   fewer than 100 cells as noise. Each surviving component gets an id and a centroid.

**(b) `dist_to_white_pulp_um` — the continuous axis, and the direct analogue of the zonation
score.** For every cell, the Euclidean distance in microns to the nearest white-pulp-nodule
**boundary** (not centroid — nodules vary in size by an order of magnitude and a centroid distance
would confound size with position). Signed: negative inside a nodule, positive in red pulp. This
is the single covariate to enter the kernel fit as zonation did.

**(c) Two supporting distances, same construction as the mouse arm's landmark set.**
`dist_to_marginal_zone_um` — distance to the nearest `marginal_zone` cell, which is where blood
enters and the compartment most likely to carry a real senescence gradient.
`dist_to_capsule_um` — distance to the tissue boundary, computed from the alpha-shape of all
QC-passing cells. Note that §10's Tier D already lists "distance to tissue boundary", so
`dist_to_capsule_um` **is** that covariate and must not be entered twice.

### 2.3 Validation, deferred, and its exact form

Run the mouse arm's §6 validation unchanged, substituting the axis:

1. Pearson r of every panel gene against `dist_to_white_pulp_um`, restricted to cells with ≥20
   counts; retain |r| ≥ 0.20 and detection ≥ 5%; **exclude any gene in Tier A or Tier B** (the
   over-adjustment guard).
2. Compare the recovered gene set against the seed compartment sets. The seeds are inside their
   own scores, so the **non-seed discoveries are the evidence**, exactly as the mouse report says.
3. Null: 300 random size-matched on-panel sets, report mean |r| and p90.
4. Per-section stability: recompute on all 7 sections and report the spread. The mouse arm derived
   zonation from **one** section and flagged that as caveat (i); with 7 spleens there is no excuse
   for repeating it.

**If step 2 fails** — i.e. the axis recovers nothing beyond its own seeds — A6 has no covariate and
that must be reported as a failed go/no-go, not patched.

### 2.4 Consequences for the frozen document

§13 test A6 currently reads "distance to the nearest conducting airway… alveolar-vs-airway-vs-
vascular compartment label… fibroblastic focus". **Every word of that is void for this arm** and
must be replaced with §2.2 before `PREREG_PHASE7.md` is committed. §17's cross-arm table needs the
same edit: the M1 column says "zonation", the H1 column must say "red pulp / white pulp axis", and
the two are *analogous, not comparable* — a λ measured against distance-to-white-pulp is not the
same physical quantity as a λ measured against a pericentral–periportal gradient.

---

## 3. Task 3 — C6: rebuilding B7 `secondary_senescence`

```bash
python3 /workspace/code/rebuild_b7_secondary_senescence.py
```

The original `B_secondary_senescence.txt` (35 genes) is **untouched**. Four variants were written
alongside it, and both candidate Tier A callers were used throughout, because the caller choice is
the PI's.

### 3.1 C6's literal prescription fails

| B7 variant | base | caller | base n | shared | after | ≥30? |
|---|---|---|---|---|---|---|
| `..._C6_minusA_s10` | v1 curated | `A_PHASE7_S10_16` | 35 | 11 | **24** | **FAIL** |
| `..._C6_minusA_ported` | v1 curated | `A_ported` | 35 | 23 | **12** | **FAIL** |
| `..._C6_sourced_s10` | SenMayo ∪ Reactome SASP | `A_PHASE7_S10_16` | 123 | 4 | **119** | PASS |
| `..._C6_sourced_ported` | SenMayo ∪ Reactome SASP | `A_ported` | 123 | 7 | **116** | PASS |

**"Rebuild without the shared genes" does not work on the curated list.** It was only 35 genes
on-panel; removing what it shares with the caller drops it below the §11 floor under either
caller — decisively so under the ported one (12). C6 as written trades one gate failure for
another.

### 3.2 The rebuild that does work

Re-source the module first, then subtract:

`B7_sourced = on-panel( SAUL_SEN_MAYO ∪ REACTOME_SASP ) − Tier A caller`
= 99 ∪ 27 → **123 on-panel** → **119** (§10 caller) / **116** (ported caller).

This is also closer to §9's *own* definition of B7 — "Tier A minus the calling genes" — than the
curated list ever was.

### 3.3 What it changes in the §11 gate

| B7 used | Tier A | \|B7\| | ≥30 | A∩B7 | \|A_strict\| | verdict |
|---|---|---|---|---|---|---|
| original v1 | §10 sixteen | 35 | PASS | 11 | 1 | FAIL |
| original v1 | ported | 35 | PASS | 23 | **21** | PASS |
| `C6_minusA_ported` | ported | 12 | **FAIL** | 0 | 33 | FAIL |
| **`C6_sourced_ported`** | **ported** | **116** | **PASS** | **0** | **33** | **PASS** |
| `C6_sourced_s10` | §10 sixteen | 119 | PASS | 0 | 5 | FAIL |

The sourced rebuild **raises the strict sender set from 21 to 33 genes** while making B7 exactly
disjoint. It does not rescue the §10 sixteen — nothing can, because `B4_downstream_arrest` alone
removes 8 of its 14 on-panel genes.

### 3.4 The citation

v1 rests on `neretti2024dissecting` — a ~250-word GSA meeting abstract, single author, no methods,
which audit finding B8 marks **"DO NOT CITE AS A PAPER"**. Replacements, all already in
`references.bib`, none invented here:

- **Module content:** `saul2022senmayo` — **PMID 35974106**, taken from the archived MSigDB record
  for `SAUL_SEN_MAYO`, not from recall. Plus Reactome pathway R-HSA-2559582 (the MSigDB record for
  `REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP` carries no PMID; cite the pathway and
  the MSigDB release).
- **The primary-versus-secondary distinction itself:** `martin2023modelling` and Acosta et al.
  2013, which are the audit's own recommended replacements at `BIO_DELIVERABLE7_CLAIM_AUDIT.md`
  §B8.

`references.bib` was **not** edited — removing an entry the audit report cross-references is a
call for whoever owns the bibliography.

---

## 4. Task 4 — Tier E2, resolved

**Recommendation, on the record: drop E2 for the human arm.** Reasoning, with the numbers:

The mouse arm's E2 (hepatocyte identity) did not merely underperform, it **failed**: its score
correlated with the zonation axis at **r = −0.67**, worse than 99% of random size-matched sets
(`genesets/README.md` §6). The cause was structural — in liver, the dominant cell-type identity
program *is* the anatomical axis. Spleen reproduces that structure exactly:

| Candidate | n | Tier A | Tier B | Tier C | A6 compartment it belongs to |
|---|---|---|---|---|---|
| `E2a_erythroid` | 5 | 0 | 1 | 0 | **red pulp** |
| `E2b_smooth_muscle_capsule` | 6 | 0 | 1 | 0 | capsule/trabecula |
| `E2c_pericyte` | 12 | 0 | 5 | 0 | capsule/trabecula |
| `E2d_structural_union` | 18 | 0 | 6 | 0 | capsule/trabecula |
| `E2e_sinusoidal_endothelium` | 9 | 0 | 0 | 0 | **red pulp** |

**Every candidate is a member of an A6 compartment by construction.** `E2a` and `E2e` are
disqualified outright: erythroid and sinusoidal identity *are* the red pulp definition, so they
would repeat the mouse failure exactly. `E2b`/`E2c`/`E2d` sit on the capsule/trabecular
compartment, which is off the red–white axis — but that compartment is captured by
`dist_to_capsule_um`, which §10 already lists as a Tier D nuisance covariate, so a "flat control"
built from it is confounded with a covariate the model conditions on.

There is also no size to work with: the largest defensible candidate is 18 genes, and the mouse
E2 had 32 with only 9 on-panel.

**What carries the control load instead**, and why it is enough:

1. **The 40 negative control probes** (`E_negative_control_probes.txt`), plus 609 negative control
   codewords and 21 genomic controls. These have **no biology by construction** — strictly better
   than any identity program — and audit test **A7** already requires them to be fitted against
   distance-to-sender and shown flat. As in the mouse arm they are sparse and must be used binned.
2. **`E_housekeeping_expanded`** (8 on-panel) for the per-cell flat-kernel test. Thin — the mouse
   arm had 13 — and it should be reported as thin.
3. **E3, the 500 size- and expression-matched random sets**, which is the real empirical null and
   is unaffected by this decision. It remains unbuildable pre-freeze (it bins genes by mean
   expression on the actual matrix).

**If the PI wants an E2 anyway**, `E2d_structural_union` (18 genes) is the least-bad, and it must
be put through the *same* test the mouse E2 failed — correlation against the A6 axis versus 300
random size-matched sets — and dropped if it behaves like the hepatocyte set. That test is
post-freeze.

---

## 5. Task 5 — does SenePy cover spleen? No.

```bash
python3 /workspace/code/senepy_coverage_human.py    # reads the installed package only
```

SenePy 1.0.1's bundled human hub file carries **65 hubs across 10 tissues**:
`blood, bone marrow, heart, hippocampus, intestine, kidney, liver, lung, skin, tongue`.

**Spleen is not among them.** There is no spleen signature to use.

The mouse arm scored SenePy with **tissue-matched** mouse Liver hubs
(`phase2_downstream.py` `HUBMAP`: `Liver/hepatocyte`, `Liver/endothelial cell of hepatic
sinusoid`, `Liver/Kupffer cell`). **The human arm cannot do this.** Every hub available for H1 is
a cross-tissue surrogate, and that asymmetry is a pre-registration fact, not a day-7 discovery.

Applying the mouse arm's own threshold (≥10 hub genes on panel):

| Spleen cell type | best surrogate hub | on-panel genes | usable |
|---|---|---|---|
| Red pulp macrophages | intestine / macrophage | 940 | yes |
| Monocytes | lung / monocyte | 167 | yes |
| Follicular B / Marginal zone B / Germinal centre B | blood / memory b cell | 267 | yes — **but one hub for all three** |
| CD4 T cells / CD8 T cells | lung / t cell | 412 | yes — **one hub for both** |
| NK cells | blood / nk cell | 107 | yes |
| Fibroblastic reticular cells / Fibroblasts | skin / fibroblast | 482 | yes — one hub for both |
| Sinusoidal endothelium / Endothelial cells | skin / endothelial cell | 164 | yes — one hub for both |
| Smooth muscle / capsule | heart / smooth muscle cell | 27 | yes |
| Pericytes | heart / pericyte | 16 | yes |
| Neutrophils | blood / neutrophil | 381 | yes |
| **cDC1, cDC2, pDC** | — | — | **no hub in any tissue** |
| **Lymphatic endothelium** | — | — | **no hub** |
| **Erythroid cells** | — | — | **no hub** |
| **Megakaryocytes** | — | — | **no hub** |
| **Mesothelial cells** | — | — | **no hub** |

**15 of 22 labels get a usable surrogate; 7 get nothing.** And the 15 collapse onto far fewer
distinct hubs than labels: the entire B compartment is scored by one blood memory-B hub, both T
subsets by one lung T-cell hub. SenePy's selling point is cell-type *specificity*; on this arm it
delivers, at best, compartment-level specificity borrowed from other organs.

**Consequences to pre-register now:**

- SenePy on H1 is **not** the same estimator as SenePy on M1. Any cross-arm comparison of the
  SenePy-called sender population is comparing a tissue-matched caller against a cross-tissue
  surrogate. `senepy_p95` appears throughout `run_phase3_*.py` and `summarize_super_callers.py`;
  every H1 number from it needs that caveat attached.
- The seven no-hub labels get **no SenePy score at all** (`NaN`), exactly as the mouse code
  already handles (`phase2_downstream.py` skips a cell type below the threshold). Sender
  prevalence under SenePy will therefore be computed over a *subset* of the section, and test A3's
  "≥200 senders, ≥5,000 non-senders" must be evaluated on that subset, not on all cells.
- Consider demoting SenePy from "caller 2" to a sensitivity analysis on this arm, and promoting
  the `CDKN1A`⁺ call — which needs no tissue-matched resource — into the primary trio. That is a
  §15 decision and is flagged, not taken.

---

## 6. Task 6 — the Tier A fork (memo, not a decision)

**Not resolved here.** Both options below are live; the §10 sixteen is not, because it already
failed A2/§11 (14 on-panel, 1 gene surviving disjointness). Numbers use the **sourced C6 B7**
where noted, since §3 changes both options in the same direction.

### Option 1 — `A_SENDER_FINAL_strict`, one sender score, disjoint from all seven modules

| | |
|---|---|
| **Size** | 21 genes; **33** with the C6 sourced B7 |
| **Buys** | One number per section. A reviewer checking §11 sees a single clean matrix of zeros. Nothing about the sender score changes when the readout changes, so the seven module fits are directly comparable to each other. Simplest possible Methods paragraph. |
| **Costs** | The 21-gene version is biologically hollow — no `CDKN1A`, `CDKN2A`, `TP53`, `LMNB1`, `MKI67`, and it still contains `CCNB1`/`FOXM1`, whose direction is opposite to the rest. The mouse arm reached the same 25-gene set and its own README calls it "numerically passing but biologically hollow". The C6 rebuild improves this materially (33 genes) but does **not** restore the canonical markers: `B4_downstream_arrest` is what removes them, and B4 is untouched by C6. |
| **Breaks it if** | A reviewer asks "why is your senescence score missing every senescence marker?" There is no good answer. |

### Option 2 — per-module sender sets, disjoint from the one module being fitted

| | |
|---|---|
| **Size** | 77 / 81 / 80 / 36 / 79 / 79 / 58 genes for B1…B7 (mouse: 70 / 74 / 73 / 37 / 73 / 71 / 55) |
| **Buys** | Every set clears ≥15 and retains the canonical markers (all 12 for B2, B5, B6). Statistically it is the correct requirement — disjointness is only needed between the sender score and the readout actually paired with it. |
| **Costs** | Seven sender scores, so "prevalence of senescent cells" is not one number. Effect sizes across modules are no longer strictly comparable because the regressor differs. Requires a clear Methods paragraph explaining why this is not p-hacking; the §11 assertion as literally written is not what is being run, and that must be declared as a deviation. |
| **Breaks it if** | A reviewer reads §11 literally and treats a per-module sender as moving the goalposts. Pre-registering the choice defuses this; deciding it after seeing H1 does not. |

### Cross-arm symmetry — the useful part

**Both options can be run symmetrically on both arms, and the files already exist.**

| | M1 mouse | H1 human |
|---|---|---|
| strict | `genesets/A_SENDER_FINAL_strict` (25) | `genesets/human/A_SENDER_FINAL_strict` (21; 33 with C6) |
| per-module | `genesets/A_sender_for_<module>` (37–74) | `genesets/human/A_sender_for_<module>` (36–81) |

Neither option forces an asymmetry, and the per-module sizes track closely between arms
(36–81 vs 37–74), so a per-module design does not make the human arm structurally different.
Re-running the mouse arm under either choice needs **no new mouse gene sets** — only a re-fit.

Three things that do bear on the decision:

1. **The C6 rebuild helps Option 1 more than Option 2** (strict 21→33; per-module sets are already
   large). If C6 is adopted, Option 1's main weakness narrows but does not close.
2. **Test A8 is unaffected either way.** Both options are gene sets on each arm's own panel, and
   the ortholog-intersected comparison is computed the same way.
3. **The mouse arm's own README already recommends Option 2** ("use these") and reports Option 1
   as the conservative sensitivity analysis. Choosing Option 1 for H1 would put the two arms'
   primary analyses out of step unless M1 is re-fitted to match — which is cheap, but must be
   decided before the freeze, not after.

**What would make one option unusable:** if the PI wants a single reported "senescent cell
prevalence" per donor for the age-axis analysis (17–59), Option 2 does not give one without an
extra arbitrary choice of which module's sender set defines prevalence. That is the strongest
argument for Option 1, and it is worth weighing against the hollow-marker objection.

---

## 7. Files

```
code/human_symbols.py                       shared panel-symbol resolver (extracted; outputs verified identical)
code/build_markers_human_spleen.py          builds the spleen marker set
code/markers_human_spleen.py                GENERATED -- 22 labels, markers_mouse_liver.py format
code/spec_a6_compartments_human.py          A6 compartment gene sets + Tier E2 candidate audit
code/rebuild_b7_secondary_senescence.py     C6
code/senepy_coverage_human.py               SenePy spleen coverage
genesets/cellmarker_pin/                    CellMarker 2.0 human, md5-pinned + PROVENANCE.md
genesets/human/markers_spleen_evidence.csv  265 rows, per-gene PMIDs and filter decisions
genesets/human/D_spleen_*.txt               5 compartment marker sets
genesets/human/B_secondary_senescence_C6_*.txt   4 C6 variants (original untouched)
results/phase7_jobA/*.log                   every run above
results/phase7_jobA/b7_c6_rebuild.{csv,json}
results/phase7_jobA/senepy_spleen_coverage.{csv,json}
results/phase7_jobA/tierE2_candidates.csv
results/phase7_jobA/a6_compartments_and_E2.json
```

## 8. What I could not do

- **Build or validate the A6 covariate** (§2) — needs the expression matrix, behind the freeze.
- **Build Tier E3**, the 500 expression-matched random sets — same reason.
- **Populate `DROP_NONSPECIFIC`** in the spleen marker module — it is an expression-based judgement
  and the mouse equivalent was made on measured data. It is empty, and labelled as empty.
- **Run the §14 image check** — the images were deliberately not downloaded (§1.4).
- **Compare against the depositors' `annotations.csv.gz`** — post-freeze, and it needs their QC
  filter characterised first.
- **Resolve the Tier A fork** (§6) or the plasma-cell `MIN_MARKERS` question (§1.3) — both are PI
  decisions and are laid out, not taken.
