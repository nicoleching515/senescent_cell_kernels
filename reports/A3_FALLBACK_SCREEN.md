# A3 fallback screen — is there a replacement human arm if H1 fails A3?

**Date:** 2026-08-27
**Scope:** contingency screen of the two runner-up candidates in `reports/PHASE7_H1_SCREEN.md` §2,
run before audit test **A3** is measured on H1, so that §18 outcome C can be decided with the
fallback known rather than guessed.
**Method:** the original screen's own, §12.1 — *verify the panel on the data, not the title* —
plus the frozen §11 / A2 disjointness gate run against each candidate's real, measured panel.
**Nothing under `data/raw_h1/` was read.** No expression value was read from either candidate:
only `/matrix/features`, `/matrix/shape`, the depositors' cell-level *label* columns, and
`cells.parquet` column names.

---

## 1. Verdict

**There is no adequate like-for-like fallback.** Neither candidate can replace H1 while preserving
the pre-registered human-arm claims. The ranking below is a ranking of *degraded substitutes*.

| Rank | Series | A2 (go/no-go) | Donors | Verdict |
|---|---|---|---|---|
| **1** | **GSE336890** — human kidney biopsy | **PASS** | **20** | **Viable but materially degraded.** Prepare this one. Costs: no chronological age axis at all, 15 of 20 donors are diseased, no depositor cell types, per-donor cell budget 12–20× smaller than H1, and the A6 anatomical covariate is confounded with donor in 18 of 20 donors. |
| **2** | GSE335963 — human bone marrow biopsy | **PASS** | **4** | **Not viable.** A2 passes, but the arm fails the A3 *cell budget* structurally, before any measurement: **no cell type in any donor can supply ≥ 200 senders and ≥ 5,000 non-senders at a prevalence of 2 %**, and **no cell type clears the floors in all four donors at any prevalence**. 4 donors is one above the bootstrap minimum. |

**The load-bearing point for the PI.** A2 — the gate §18 outcome C actually names — is *not* what
kills either candidate. Both pass it cleanly, with the same 33-gene strict Tier A on-panel and all
seven Tier B modules above the ≥ 30 floor. What kills them is the very test that would have
triggered the fallback:

- The three panels are near-identical: the **three-way core is 5,006 genes** (H1 5,093; both
  candidates 5,101; H1 ∩ kidney = 5,008, H1 ∩ bone marrow = 5,023). The sender callers are the
  frozen ones. **An A3 failure on H1 caused by the caller producing out-of-band prevalence on a
  Prime 5K panel is therefore likely to reproduce on both candidates** — the panel and the caller
  are held nearly constant, only the tissue changes.
- And each candidate *adds* a failure mode H1 does not have. On H1 a cell type needs only
  **2.4–2.5 %** of a donor's cells to clear both A3 floors at p = 5–10 %. On the kidney median
  donor it needs **29–31 %**; on the bone-marrow median donor **7.7–8.1 %**, and measurement
  against the depositors' own labels shows only 9 of 88 (donor × type) strata reach it.

So: if A3 fails on H1, redirecting to either candidate is *more* likely to fail A3 again than to
succeed, and the honest reading is that **an A3 failure ends the human arm rather than moving it**.
GSE336890 is worth preparing only as a reduced, disease-defined arm that the PI has explicitly
accepted as such — not as a substitute that preserves the pre-registration.

---

## 2. Two corrections to the runner-up table in `PHASE7_H1_SCREEN.md` §2

Both rows of that table state a **sample** count where a **donor** count is needed, and both are
wrong in the direction that matters.

| Series | Table says | Measured here | Command |
|---|---|---|---|
| GSE336890 | "9 human kidney biopsies" | **9 Xenium *Region slides*, carrying 22 tissue blocks from 20 distinct patient specimens** — the series' own overall design says "each Xenium run (Region) imaged one slide that could contain multiple tissue sections from different patients" | `python3 code/screen_candidate_cellbudget.py` |
| GSE335963 | "39 human bone marrow samples … h5 + parquet" | **A SuperSeries.** 39 GSM total, but only **6** are GPL33762 Xenium (`SuperSeries of: GSE335962`); the other 33 are 10x scRNA-seq/TCR on GPL34284 and ship no h5/parquet. The 6 Xenium sections come from **4 donors** — `NC03`/`NC03-2` and `NC05`/`NC05-2` are repeat sections of the same donor | `grep '^!Series_relation' data/raw_h1_candidates/GSE335963.soft.txt` and the sample table below |

Consequence: the table's own next row, "GSE335962 | 6 | Bone marrow biopsy (CHIP) | Subset of the
above", is the *whole* of GSE335963's Xenium content. Choosing GSE335963 over GSE335962 gains
nothing. The 132-series platform count was re-verified today and is unchanged
(`<Count>132</Count>`, `results/a3_fallback/gpl33762_count.xml`).

---

## 3. Panel verified on the data (§12.1 step 2)

Downloaded with `bash code/fetch_h1_candidates.sh` (4 + 6 `cell_feature_matrix.h5`, ~141 MB total
including metadata), verified with `python3 code/screen_candidate_panels.py`
(log: `results/a3_fallback/panel_screen.log`, machine-readable: `panel_screen.json`).

| Feature class | H1 GSE326743 | GSE336890 kidney | GSE335963 bone marrow |
|---|---|---|---|
| Gene Expression (ENSG-prefixed) | 5,093 | **5,101** | **5,101** |
| Negative Control Codeword | 609 | **609** | **609** |
| Negative Control Probe | 40 | **40** | **40** |
| Genomic Control (`Intergenic`) | 21 | **21** | **21** |
| Unassigned Codeword | 695 | 582 | 630 |
| Deprecated Codeword | 3,291 | 3,291 | 3,291 |
| duplicate gene symbols | 0 | **0** | **0** |
| panel identical across samples | yes (3 checked) | **yes, 4 of 4, symmetric difference 0** | **yes, 4 of 4, symmetric difference 0** |

Both candidates carry the **stock Prime 5K control complement 609 / 40 / 21 exactly**, so the 40
negative control probes that audit test **A7** requires are present in both. Both are one stock
panel across samples, not per-sample custom designs. H1's per-feature-class figures in the table
above are quoted from `genesets/h1_candidate/PROVENANCE.md`; they were **not** re-derived here,
because `data/raw_h1/` is out of scope for this screen.

**A data defect in GSE336890 that a title-level or `feature_type`-level screen would have gotten
wrong.** In two of the four kidney samples (`Region01`, `Region15`) the h5 `feature_type` column
labels **every one of the 9,644 features `Gene Expression`**, controls included — a naive count
returns a 9,644-gene "panel". The feature *ids* are intact in those files
(`ENSG…` / `NegControlProbe…` / `NegControlCodeword…` / `Intergenic…` / `UnassignedCodeword…` /
`DeprecatedCodeword…`), and restricting to the ENSG prefix recovers exactly 5,101 genes in all
four samples with symmetric difference 0. `code/screen_candidate_panels.py` therefore classifies
by id prefix and reports both counts; the panel CSVs it writes are the ENSG-restricted ones.

**Panel relationships (measured).** Kidney ^ bone marrow = 164 (82 genes each way); H1 vs each
candidate = 5,008 and 5,023 shared. **Three-way core = 5,006 genes.** All three are the same stock
predesigned base with a different ~90–95-gene add-on, which is the structural reason §1's
"the panel is nearly held constant" argument holds.

---

## 4. The A2 disjointness gate — go/no-go — **PASS on both**

Run with `code/gate_disjointness_candidate.py`, which does **not** reimplement the gate: it reads
`code/gate_disjointness_human.py` verbatim and executes it with exactly two constants rebound
(`PANEL` → the candidate panel CSV, `OUT` → `results/a3_fallback/gate_<SERIES>/`). Every assertion,
Tier A/B definition, module and the per-module sender sensitivity are the frozen ones in
`genesets/human/`. `results/phase7_jobA/` was not written to.

```
python3 code/gate_disjointness_candidate.py GSE336890 results/a3_fallback/GSE336890_gene_panel_5101.csv
python3 code/gate_disjointness_candidate.py GSE335963 results/a3_fallback/GSE335963_gene_panel_5101.csv
```
Logs: `results/a3_fallback/gate_GSE336890.log`, `gate_GSE335963.log`; verdicts in each
`gate_<SERIES>/gate_result_human.json`. **Both exit 0: `### FROZEN CONFIGURATION VERDICT: PASS ###`.**

| Quantity | H1 (frozen) | GSE336890 | GSE335963 |
|---|---|---|---|
| panel used | 5,093 | 5,101 | 5,101 |
| **strict Tier A `A_SENDER_FINAL_strict` on-panel** | **33 / 33** | **33 / 33** | **33 / 33** |
| Tier A `A_ported`: raw → removed → final | 81 → 48 → **33** | 81 → 48 → **33** | 81 → 48 → **33** |
| per-module sender gate (all 7, ≥ 15 and disjoint) | PASS | **PASS** | **PASS** |
| `A_PHASE7_S10_16` (reported, not used) | FAIL (5) | FAIL (5) | FAIL (5) |

Tier B module sizes on panel (floor is ≥ 30):

| Module | H1 | GSE336890 | GSE335963 |
|---|---|---|---|
| B1 tnfa_nfkb_proximal | 120 | 113 | 116 |
| B2 il6_jak_stat3 | 71 | 69 | 71 |
| B3 interferon_response | 126 | 121 | 122 |
| B4 downstream_arrest | 231 | 231 | 231 |
| B5 emt_ecm | 113 | 101 | 104 |
| B6 oxidative_stress | 36 | **36** | **36** |
| B7 secondary_senescence | 116 | 108 | 113 |

Every module clears the floor on both candidates; the tightest, `oxidative_stress`, has the same
+6 margin it has on H1. Per-module sender-set sizes are identical to H1's on both candidates
(77 / 81 / 80 / 36 / 79 / 79 / 81), because the 33-gene strict Tier A and all seven per-module
sender sets are fully on-panel in all three.

CoreScence circularity, computed by the same frozen code path: H1 29/33 = 87.9 %; kidney
27/34 = 79.4 %; bone marrow 28/35 = 80.0 %. No candidate makes the circularity problem worse.

**A2 is not the discriminator.** It passes on both. This is worth stating plainly because §18
outcome C is triggered by A2 or A5, and a reader could otherwise assume the runners-up were held
back by the gate. They are not.

---

## 5. Donors and the A3 cell budget

`python3 code/screen_candidate_cellbudget.py` → `results/a3_fallback/cell_budget.{csv,json,log}`,
`GSE336890_donor_blocks.csv`, `required_type_share.json`.

### 5.1 GSE336890 — kidney, 20 donors

Donor identity is **not** the sample count. The depositors ship
`GSE336890_cells_stats_dir.tar.gz`, one CSV per (patient × kidney region) block; there are 22
blocks over 9 Region slides carrying **20 distinct patient specimens**, **358,005 cells** total.
The block sums reconcile against the h5 files exactly for two of the four Regions downloaded
(Region_1 121,864 = 121,864; Region_15 98,714 = 98,714) and to within 4 and 148 cells for the
other two (Region_2, Region_14) — i.e. the per-cell donor assignment is essentially complete.

- **8 AIN, 7 ATI, 5 reference** donors. **The only non-diseased donors are the 5 references**
  (L102, L111, L112, L115, L117).
- Per-donor cells: **min 2,799, median 17,972, max 40,589**. 19 of 20 donors have ≥ 5,000 cells
  *in total*.
- Cortex/medulla: 15 cortex blocks (249,362 cells), 7 medulla (108,643). **Only 2 of 20 donors
  contribute both compartments.**
- **No cell-type annotation is deposited**: the cells_stats CSVs carry `CellID, Transcripts, Area`
  and nothing else.
- `cells.parquet` carries every Tier D nuisance column the plan asks for — verified on
  `GSM9844157_Region02_cells.parquet.gz`: `transcript_counts, control_probe_counts,
  genomic_control_counts, control_codeword_counts, unassigned_codeword_counts,
  deprecated_codeword_counts, total_counts, cell_area, nucleus_area, nucleus_count,
  segmentation_method`; centroids in microns (x 178–9,528; y 8–2,189).

### 5.2 GSE335963 — bone marrow, 4 donors

| GSM | Section | Donor | CHIP | Sex | Cells (h5 `/matrix/shape`) |
|---|---|---|---|---|---|
| GSM9824312 | CH02 | CH02 | CHIP | M | 61,583 |
| GSM9824313 | CH15 | CH15 | CHIP | M | 40,264 |
| GSM9824314 | NC03 | NC03 | Non-CHIP | M | 23,992 |
| GSM9824315 | NC03-2 | NC03 | Non-CHIP | M | 51,464 |
| GSM9824316 | NC05 | NC05 | Non-CHIP | F | 35,084 |
| GSM9824317 | NC05-2 | NC05 | Non-CHIP | F | 69,901 |

**6 sections, 4 donors, 282,288 cells.** Per-donor: CH02 61,583; CH15 40,264; NC03 75,456;
NC05 104,985. 2 CHIP vs 2 non-CHIP, 3 male / 1 female. **No age is deposited** for any donor.

Unlike kidney, the depositors ship cell-level labels (`*_Metadata.csv.gz`, `Annotation` column):
**22 cell types**, covering **183,127 of 282,288 cells (64.9 %)** — a QC filter of the same kind
H1's `annotations.csv.gz` applies, and one that would have to be characterised before use.

### 5.3 The A3 arithmetic

A3 (Master Plan Test 3) passes only if some threshold gives prevalence in 1–20 % **with ≥ 200
senders and ≥ 5,000 non-senders per donor**, reported **per cell type**. That fixes a floor on
cell counts alone, independent of any measurement: a (donor × type) stratum needs

    n ≥ max( 5000 / (1 − p) , 200 / p )
      = 6,250 at p = 20 % ; 5,556 at 10 % ; 5,263 at 5 % ; 10,000 at 2 % ; 20,000 at 1 %.

**Required share of a donor held by a single cell type**, from measured per-donor cell counts
(`results/a3_fallback/required_type_share.json`; H1 row computed from the per-sample counts in
`PHASE7_H1_SCREEN.md` §3, whose stated total 2,207,593 was re-checked and is correct):

| Arm / donor | p = 20 % | p = 10 % | p = 5 % | p = 2 % | p = 1 % |
|---|---|---|---|---|---|
| **H1** min donor (220,435 cells) | 2.8 % | 2.5 % | 2.4 % | 4.5 % | 9.1 % |
| **H1** max donor (396,173) | 1.6 % | 1.4 % | 1.3 % | 2.5 % | 5.1 % |
| GSE336890 min donor (2,799) | 223 % | 199 % | 188 % | 357 % | 715 % |
| GSE336890 **median donor (17,972)** | **34.8 %** | **30.9 %** | **29.3 %** | **55.6 %** | 111 % |
| GSE336890 max donor (40,589) | 15.4 % | 13.7 % | 13.0 % | 24.6 % | 49.3 % |
| GSE335963 min donor (40,264) | 15.5 % | 13.8 % | 13.1 % | 24.8 % | 49.7 % |
| GSE335963 **median donor (68,520)** | **9.1 %** | **8.1 %** | **7.7 %** | 14.6 % | 29.2 % |
| GSE335963 max donor (104,985) | 6.0 % | 5.3 % | 5.0 % | 9.5 % | 19.1 % |

A share above 100 % means **no** cell type in that donor can satisfy A3 at that prevalence: the
smallest kidney donor (REF/L112, 2,799 cells) is out at every prevalence, and the median kidney
donor is out at p = 1 %.

**Bone marrow, measured rather than bounded** (`python3 code/screen_candidate_a3_budget.py` →
`GSE335963_a3_budget.json`, `GSE335963_donor_x_celltype.csv`). Using the depositors' own 22 labels,
of **88 (donor × type) strata**:

| p | cells needed | strata clearing (of 88) | types clearing in **all 4** donors |
|---|---|---|---|
| 20 % | 6,250 | 6 | none |
| 10 % | 5,556 | 8 | none |
| 5 % | 5,263 | 9 | none |
| **2 %** | **10,000** | **0** | none |
| 1 % | 20,000 | 0 | none |

The largest single stratum is Erythroid (Cycling) in CH02 at 9,155 cells — short of the 10,000 the
2 % floor needs. Scaling every count up by that donor's annotation coverage (an optimistic bound
that assumes Job B labels every unlabelled cell with the same composition) improves this to
17 of 88 at p = 5 % and **6 of 88 at p = 2 %**, still with **no cell type clearing in all four
donors at any prevalence**. The Master Plan's own sweet spot is 2–10 %. **Bone marrow fails the A3
budget structurally.**

For kidney the equivalent measurement is impossible before Job B, because no cell types are
deposited. What can be said from measured counts: at the sweet-spot p = 5 %, a kidney cell type
must be **≈ 29 % of the median donor** to clear both floors — attainable in principle by proximal
tubule in a cortex biopsy and by essentially nothing else, in the larger donors only. Kidney's
plausible A3 outcome is therefore *one or two dominant epithelial types in perhaps half the
donors*, against H1 where a 2.4 % share suffices and many types qualify.

---

## 6. The biological axis, honestly

Both candidates are disease-defined. Neither has a chronological ageing axis: **age is not
deposited for any donor of either series** — GSE336890's characteristics are
`tissue / tissue preservation method / disease / kidney region`, GSE335962's are
`tissue / treatment / gender`. §12.1's stated reason for choosing H1 over exactly these two —
"a real chronological human ageing axis (17–59), which §12.1 argues is what the senescence field
actually cares about and which the mouse arm cannot give" — is simply gone under either.

### 6.1 What H1's §15 declared deviations would have to become

| H1 deviation, as declared | Under GSE336890 (kidney) | Under GSE335963 (bone marrow) |
|---|---|---|
| "Tissue differs between arms, so no cross-arm difference can be attributed to species or tissue" | Stands, and **must be strengthened**: liver-vs-kidney, mouse-vs-human, **and normal-vs-diseased**. Three confounded axes, not two. | Same, with CHIP-vs-non-CHIP as the third axis. |
| "No disease/treatment contrast at all; all seven are normal" (H1) | **Inverts.** 15 of 20 donors are diseased (8 AIN, 7 ATI); only 5 references. A new deviation is required: a surviving distance-to-sender kernel in AIN/ATI tissue cannot be separated from acute inflammatory injury. The depth-vs-prevalence reframing §12.1 made for H1 has to be redone. | **Inverts.** 2 CHIP vs 2 non-CHIP. With n = 2 per group the disease contrast is not analysable, so the arm carries the confound without being able to test it. |
| "A real chronological ageing axis, 17–59" | **Withdrawn — no age at all.** | **Withdrawn — no age at all.** |
| "The old end of the age axis is thin (n = 2 over 55)" | Superseded; there is no axis to be thin. | Superseded. |
| "Sex is unbalanced (4M/3F), race unknown for two donors" | **Worse: sex is not deposited at all** for any of the 20 donors. | 3M / 1F over 4 donors. |
| A6: "spleen has a clean anatomical axis, red pulp vs white pulp, supplying the arm-specific covariate test A6 demands" | **Partial analogue, seriously compromised** — see §6.2. | **Weak analogue** — see §6.2. |
| P4: "SenePy ships no spleen signature; it is not the same estimator across arms" | **Improved but not resolved** — see §6.3. | **Barely improved** — see §6.3. |
| §12.3 deviation: "the depositors ship `annotations.csv.gz`, four nested levels — an external label set to check Job B against" | **Lost.** No cell-type annotation is deposited. | **Retained**, 22 types, but on only 64.9 % of cells. |

### 6.2 Is there an A6 anatomical covariate?

**Kidney: yes in name, compromised in practice.** Cortex vs medulla is a textbook anatomical axis
and the depositors label it per block, so every cell already carries it — no marker-based
compartment score has to be built, which is *easier* than H1's red/white pulp job (task 9.5).
But: **only 2 of 20 donors contribute both cortex and medulla** (ATI/7997 and REF/L102). For the
other 18 the compartment is a property of the donor, so the A6 covariate is **confounded with
donor and cannot be fitted within section**, which is what the kernel machinery needs — the role
liver zonation plays in M1 and red/white pulp plays in H1 is a *within-section* gradient. Getting a
within-section covariate in kidney means building a new marker-based compartment set
(glomerulus / proximal tubule / distal-collecting / interstitium), i.e. redoing the §14 A6 spec
work from scratch for a new tissue.

**Bone marrow: a within-section covariate exists but is tiny.** The `QuPath_Annotation` column
gives, over all 6 sections: `Other` 76,978, `Tissue` 58,945, `Ignore*` 39,531, `Fiber-enriched`
5,721, `Bone` 1,952 (`results/a3_fallback/GSE335963_qupath_regions.json`). The endosteal
(`Bone`) compartment — the biologically meaningful niche axis — is **1.07 % of annotated cells**,
and `Fiber-enriched` is concentrated in one section (4,157 of 5,721 in CH15). 21.6 % of annotated
cells are `Ignore*`. This will not support an A6 covariate at the standard H1 sets.

### 6.3 Does SenePy have tissue-matched signatures?

Deviation P4 records SenePy 1.0.1 as 65 human hubs over 10 tissues with **no spleen hub**, so on
H1: 0 tissue-matched, 15 cross-tissue surrogates, 7 labels with no hub at all. Re-derived here for
the candidates with `python3 code/senepy_coverage_candidates.py`
(→ `senepy_coverage_candidates.json`, `senepy_coverage_GSE*.csv`), using
`code/senepy_coverage_human.py`'s method and the mouse arm's own `MIN_ON_PANEL = 10`. The 10
tissues are: blood, **bone marrow**, heart, hippocampus, intestine, **kidney**, liver, lung, skin,
tongue.

| | H1 spleen | GSE336890 kidney | GSE335963 bone marrow |
|---|---|---|---|
| tissue-matched hub present | **no** | **yes** | **yes** |
| hubs for that tissue | 0 | 5 | **1** |
| distinct cell types covered | 0 | **3** (macrophage, epithelial cell, endothelial cell) | **1** (t cell) |
| usable at ≥ 10 genes on panel | 0 | 4 of 5 (69, 23, 19, 11 genes on panel) | 1 (269 genes on panel) |

So the P4 asymmetry is **narrowed by kidney and essentially untouched by bone marrow**. Kidney gains
tissue-matched hubs for exactly three generic cell classes — every tubular subtype, podocyte,
mesangial, fibroblast, T/B/NK label would still take a cross-tissue surrogate, and all the
epithelial subtypes would collapse onto one generic kidney "epithelial cell" hub, reproducing the
surrogate-collapse problem P4 already documents for spleen ("one blood memory-B hub scores all
three B labels"). Bone marrow's single tissue-matched hub is a T-cell hub — nothing for erythroid,
neutrophil, GMP, megakaryocyte, MSC or osteo lineages, which are 70 % of the tissue.

### 6.4 What would have to be rebuilt

Either candidate invalidates the tissue-specific half of the frozen human gene sets. Of the 43
files in `genesets/human/FROZEN_MANIFEST.csv`, **5 are spleen-specific**
(`D_spleen_red_pulp`, `D_spleen_white_pulp_follicle`, `D_spleen_white_pulp_tzone`,
`D_spleen_marginal_zone`, `D_spleen_capsule_trabecula`), plus `markers_spleen_evidence.csv` and
`code/markers_human_spleen.py`'s 22-type label set and the A6 compartment spec
(`results/phase7_jobA/a6_compartments_and_E2.json`). Tiers A, B, C and E — everything the A2 gate
tests — port unchanged, which is why A2 passes. The tissue tier does not port at all. That is a
fresh Job A for a new tissue, a new marker-evidence build, a new A6 spec, and a **re-freeze**:
`genesets/` is under a SHA-256 manifest with a pre-commit hook and a PostToolUse hook, so the
`phase8-frozen` pre-registration would have to be superseded, not amended.

---

## 7. Recommendation

1. **Run A3 on H1 as planned.** Nothing in this screen argues for pre-emptive redirection; both
   candidates are worse on A3 than H1 is, by construction.
2. **If A3 fails on H1, do not treat redirection as the default remedy.** Establish *why* it failed
   first. If prevalence was out of the 1–20 % band, that is a caller-on-Prime-5K property and the
   candidates share 98 % of the panel and 100 % of the caller definitions — expect it to recur.
   Only if the failure is specifically spleen-compositional (e.g. no non-lymphoid type reaches the
   sender floor) does a tissue change plausibly help, and then only GSE336890.
3. **If the PI directs a fallback anyway, it is GSE336890**, prepared explicitly as a reduced
   disease-defined arm: no age axis, 5 reference donors, A6 to be rebuilt as a within-section
   marker-based compartment set, Job B annotation from scratch with no depositor labels to check
   against, and A3 expected to qualify only one or two dominant epithelial types.
4. **GSE335963 should be struck from the runner-up list**, and the table corrected: it is a
   SuperSeries whose Xenium content is 6 sections from 4 donors, and it fails the A3 cell budget
   before any measurement.
5. **The §18 outcome C paragraph is stronger for this screen having been run**: the two runners-up
   **pass** the gate that outcome C names and fail instead on donor-level power, axis and cell
   budget. "The alternatives clear the disjointness gate and still cannot support the test" is a
   sharper data-availability finding than "no replacement was found".
6. **One caveat on the wider claim.** `PHASE7_H1_SCREEN.md` §3 states H1 is the "only normal-tissue
   Prime 5K human deposit with ≥ 3 donors found in the 132". This screen did **not** re-verify that
   across all 19 Prime-5K series — and it did find that the donor counts in that table are sample
   counts, wrong in both directions for the two rows checked. The nearest challenger in the table,
   GSE319763 (3 human, "Lung, normal + AML-infiltrated"), would need the same donor-level check
   before the "only one" wording goes into a paper.

---

## 8. Provenance and hygiene

- **New code** (screening only, nothing existing modified): `code/fetch_h1_candidates.sh`,
  `code/screen_candidate_panels.py`, `code/gate_disjointness_candidate.py`,
  `code/screen_candidate_cellbudget.py`, `code/screen_candidate_a3_budget.py`,
  `code/senepy_coverage_candidates.py`.
- **New outputs:** `results/a3_fallback/` (panel CSVs, gate outputs under `gate_GSE*/`, cell budgets,
  SenePy coverage, logs). `results/phase3/`, `figures/`, `genesets/` and existing `code/` untouched;
  `results/phase7_jobA/` not written to.
- **Downloads:** `data/raw_h1_candidates/`, 141 MB, added to `.gitignore` under the same policy as
  `data/raw_h1/`. Re-fetch with `bash code/fetch_h1_candidates.sh`. **`data/raw_h1/` was not read.**
- **No package was installed.** DeepScence 1.0.0 and senepy 1.0.1 were already present and were
  used as installed.
- **Figures guard:** `python3 code/check_figures_guard.py` → `OK: all 52 committed figures match`.
- **Memory:** peak observed `/sys/fs/cgroup/memory.current` 4.91 GB against a `memory.max` of
  57.74 GB.
- **No commit, no push, no tag.** `reports/SUBMISSION_PATCH_2026-08-29.md` shows as modified in
  `git status`; that change was not made by this screen and has been left alone.
