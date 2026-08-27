# Phase 7 — Human Replication, Minimal Path

> ## ⚠ SUPERSEDED IN PART — read before using any number or dataset name here
>
> **H1 is NOT "human aging lung".** The GEO screen (132 GPL33762 series) selected
> **GSE326743 — 7 normal human SPLEENS**, ages 17–59, Prime 5K + 100 addon,
> panel verified on the data at **5,093 genes**. `GSE335761`, this document's
> primary candidate, correctly FAILED the screen at 386 targets. Every "aging
> lung" reference below — including §12.1, §13 test A6 and the §17 table — names
> the wrong tissue. See `reports/PHASE7_H1_SCREEN.md`.
>
> **Also superseded by Phase 8:**
> - **§10's 16-gene Tier A FAILS the §11 disjointness gate** (14 on-panel, 13
>   inside Hallmark modules, `ATR` alone surviving). Replaced by a strict-33 set,
>   with per-module sets as pre-registered sensitivity.
> - **§13 test A6** is written for lung. The spleen analogue is a red-pulp /
>   white-pulp axis.
> - **§17's table is entirely pre-C6** and its numbers no longer hold.
> - **§23/C6's stated method does not work** (yields 12/24 genes, below the floor).
> - **§4 (D-b)'s premise is REFUTED**: DCA does not normalise depth — measured, it
>   roughly doubles the depth loading on 3 of 3 sections.
> - **§15's composition-matched protocol had no implementation** and, once built,
>   is **inert** (1.6% vs 85.4% from the same variables as covariates).
>
> This document is retained as the **plan of record for Phase 7**. Current state:
> `reports/PHASE8_ROADMAP_STATUS.md`; frozen parameters: `reports/PREREG_PHASE8.md`.


## Fix the null, fix the sender caller, annotate a human Xenium Prime 5K dataset, run the frozen pipeline

**Addendum to:** `SASP_Kernel_Master_Plan.md`
**Supersedes:** the SenNet Phase 7 draft and the three-arm Phase 6R/7 plan
**Status:** post-submission. Start after the NeurIPS workshop deadline (Aug 29, 2026)
**Scope:** ~10 working days, 2 arms, 5 steps
**Last updated:** August 26, 2026

---

## Contents

**Part 0 — Orientation**
- 0.1 What This Is
- 0.2 The Two Kinds of Annotation
- 0.3 Why Only Two Corrections
- 0.4 The Rule

**Part I — Step 1: Fix the Torus Shift (C1)**
- 1. The Problem
- 2. The Fix
- 3. What Each Outcome Means

**Part II — Step 2: Fix the Sender Caller (C7)**
- 4. The Four Problems
- 5. D1 — Complete Coverage
- 6. D2 — Resolve Denoise
- 7. D3 — Fix the Sign Anchor
- 8. The Free Experiment

**Part III — Step 3: Annotate the Panel (Job A)**
- 9. Gene Sets in Human Symbols
- 10. The Disjointness Gate

**Part IV — Step 4: Acquire and Annotate the Human Arm**
- 11. What to Download
- 12. Download Discipline
- 13. Acquisition Audit
- 14. Annotate the Cells (Job B)

**Part V — Step 5: Run and Compare**
- 15. Freeze and Pre-Register
- 16. Run Order
- 17. Two-Arm Comparison Table
- 18. Outcomes
- 19. Figures

**Part VI — Logistics**
- 20. Timeline
- 21. Compute
- 22. Risks
- 23. The Optional Menu

---

# Part 0 — Orientation

## 0.1 What This Is

Five steps, in order:

1. **Fix the torus shift** so the null your third contribution rests on is not mis-specified.
2. **Fix DeepScence** so the sender caller works properly before you point it at new data.
3. **Annotate the human panel** — build Tiers A–E in human symbols and prove Tier A ∩ Tier B = ∅.
4. **Download and annotate a human Xenium Prime 5K dataset** — call senescent cells with DeepScence, since the public data ships with no senescence labels.
5. **Run the frozen pipeline and report both arms side by side.**

Roughly ten working days. Two arms: your existing mouse liver (M1) and one human dataset (H1).

## 0.2 The Two Kinds of Annotation

The word "annotate" has been doing double duty and it caused confusion. Two different jobs:

**Job A — annotate the panel.** Take the 5,001 human genes and assign them to Tier A (sender definition: arrest and damage only), the seven Tier B modules (response), Tier C (ligand–receptor), Tier D (nuisance covariates), Tier E (controls). Gene-level, one-time per species, a day of work. Part III.

**Job B — annotate the cells.** For every cell in the dataset, produce a senescence score, so you have senders to measure distance from. Public Xenium deposits carry no senescence labels, so you generate them. Part IV §14.

**Job B is what DeepScence does.** Sender calling and cell annotation are the same step. That is why Part II is not an add-on — it is the tool that produces the annotation on the new data, and it currently has problems you would inherit.

## 0.3 Why Only Two Corrections

The full audit found seven issues. Two are load-bearing:

- **C1 (torus shift)** because contribution 3 currently has a competing explanation documented in your own Phase 4 appendix.
- **C7 (DeepScence)** because it is the sender caller for the new arm.

The other five (C2 through C6) are real but cosmetic or presentational, and they are in §23 as a menu. Do them if there is time; do not let them delay the replication.

## 0.4 The Rule

> **Freeze the pipeline and commit the pre-registration before you download the human data.**

You know what mouse said. Every decision made afterward is contaminated by that knowledge. The mouse paper's credibility rests substantially on reporting against interest — the winner's-curse correction going the wrong direction, zonation failing to be the confound the plan predicted, the `denoise=False` deviation declared rather than buried. Preserve that by deciding before you can see which way it cuts.

---

# Part I — Step 1: Fix the Torus Shift (C1)

## 1. The Problem

`code/run_phase3_nulls.py` computes the torus shift over the whole-section bounding box:

```python
lo = sf.coords.min(0); hi = sf.coords.max(0)
def torus_shift(rng, pts, lo, hi):
    span = hi - lo
    return lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span
```

A liver section is not a rectangle. Wrapping over its bounding box places shifted senders where there are no cells.

**Your own `CS_PHASE4.md` §2.4 says exactly this:**

> "A torus shift on a whole section throws ~20 % of shifted cells into the void outside the tissue, which would weaken the shifted condition for the wrong reason."

Phase 4 fixed it by tiling onto solid tissue and verified the fix (100% of shifted ligand⁺ cells retained a neighbour within 100 µm; median neighbour count real 150.8 vs shifted 149.4). **Phase 3 was never re-run.**

`CS_PHASE3.md` argues N3 SF = 1.000 is structural — the shift destroys alignment but not shared causes. That may be true. But the artifactual explanation is sitting in your repository, fixed elsewhere, and a reviewer who reads both reports will find it.

## 2. The Fix

Implement three variants and report all three:

| Variant | Implementation |
|---|---|
| **N3-tile** | Torus shift within Phase 4-style solid-tissue tiles |
| **N3-occ** | Whole-section shift, rejected and resampled if >5% of senders land outside a 25 µm occupancy grid |
| **N3-swap** | Relocate each sender to a randomly chosen **real cell position**, preserving the count. Every shifted sender is in tissue by construction |

**Report the destructiveness diagnostic for each**, exactly as Phase 4 does: fraction of shifted senders retaining a neighbour within 100 µm, and median neighbour count real versus shifted. Write to `results/phase3/null_destructiveness.csv`.

**Apply the same fix to N4 (rotation).** `rotate_about_centroid` uses the identical `% span` and has the identical problem.

**Also, ten-minute fix while you are in there:** `make_figure4.py` line 42 hardcodes the Phase 3 values in an `OURS` dict. They currently match, but the figure will silently go wrong after this re-run. Make it read from `results/phase3/sf_summary.csv`.

## 3. What Each Outcome Means

- **SF stays ≈ 1.0 under all three variants** → the structural claim is vindicated and much stronger than it is now. Likely, and worth the day.
- **SF drops materially under N3-occ or N3-swap** → the original result was partly an artifact. Report the corrected value, scale the calibration-failure framing back to what it supports, and note that the void effect is itself a warning to anyone running torus shifts on non-convex tissue. Still a contribution, just smaller.

**Cost:** one day of implementation, under a day of compute.

---

# Part II — Step 2: Fix the Sender Caller (C7)

DeepScence is the only senescence caller built for single-cell and spatial data, so it is the right primary sender definition. Right now it is the weakest link in the project.

## 4. The Four Problems

**D-a — Coverage. Scores exist for two of eleven sections.** `BIO_PHASE3.md` §5: the 26-week pair only, the other nine ~50 minutes each. **Every DeepScence-based number in the paper rests on one sham section (7250) and one SBR section (7259)** — and 7250 is the section Test 3 excludes for falling below the 1% prevalence floor. The caller-agreement result that motivates the paper is a two-section result.

**D-b — `denoise=False`.** DCA needs an obsolete TensorFlow stack. Scores are computed on undenoised counts, and DCA denoising is precisely the step that would normalize depth — the confound under investigation.

**D-c — Ortholog remapping.** DeepScence ships a **human** CoreScence v2 gene set and hardcodes `CDKN1A` for direction fixing. `run_deepscence.py` maps mouse symbols to human 1:1 orthologs via the MGI `HOM_MouseHumanSequence` report, the same resource the DeepScence authors use in their own Zenodo notebook. **4,845 of 5,097 panel genes map; 252 are dropped.**

**D-d — Sign anchoring.** DeepScence fixes its bottleneck sign by correlating with `CDKN1A`. Two consequences:

1. **Circularity with your `Cdkn1a`⁺ caller is built in.** The one pair consistently above chance is DeepScence vs `Cdkn1a`⁺ at 1.51–2.85×, and `BIO_PHASE3.md` §4.4 already identifies why. That pair is not evidence of agreement.
2. **The polarity flips between arms.** Correlation between the `ds` score and its own gene set: **−0.350 in sham, +0.318 in SBR.** In sham it inverts relative to its own gene set.

## 5. D1 — Complete Coverage

**Non-negotiable, and it is only compute.** Run DeepScence on all 11 M1 sections and every section of H1. At ~50 min/section, roughly 9 hours for M1, parallel across sections. One overnight run.

Then regenerate at full coverage: `caller_pairwise_agreement.csv`, `caller_agreement_depth_and_type_matched.csv`, `caller_celltype_composition.csv`, `caller_within_type_depth_bias.csv`, `caller_technical_loading.csv`.

**Report the two-section values alongside the eleven-section values** so readers can see what the old base gave. **If agreement rises above chance at full coverage, your motivating claim weakens and must be restated.** That is a real possibility and it is why this comes before the human arm, not after.

## 6. D2 — Resolve Denoise

In order, stop when one works:

1. **Isolated environment.** A separate conda env pinned to DCA's TensorFlow version, running only the denoising step, handing back an `.h5ad`. Cleanest; does not contaminate the main pinned stack.
2. **Container.** DCA in its own Docker image, same handoff.
3. **Accept and quantify.** If neither works, keep `denoise=False` and **measure what it costs**: score the same cells with and without an equivalent depth normalization of your own (scran, or median-of-ratios before scoring) and report how much the score moves. That converts an unmeasured caveat into a measured one.

**Report both configurations for at least one section** whichever path lands.

## 7. D3 — Fix the Sign Anchor

1. **Re-anchor outside the circular set.** Fix the bottleneck sign using an anchor in neither Tier A nor any Tier B module — `Lmnb1` (down in senescence), proliferation markers (down), or the consensus of the other three callers. **Report both the published `CDKN1A` anchor and the re-anchored version**; name the primary in the pre-registration.
2. **Add a sign-invariant summary.** Rank by |score| as well as by score, so nothing silently depends on a polarity that flips.
3. **Promote the flip to Results.** "A published senescence caller's polarity inverts between the two surgical arms of the same study, because its sign is anchored on a single gene whose behaviour differs between them" is concrete, verifiable, and useful. It belongs in the paper, not in §4.3 of a phase report.

## 8. The Free Experiment

**CoreScence is human.** On mouse it runs remapped; **on H1 it runs natively, exactly as published, with no remapping.**

That separates two explanations currently confounded:

| Observation | Mouse arm only | Both arms |
|---|---|---|
| Sign instability across conditions | Artifact of remapping or `denoise=False` | Property of the published tool |
| Depth-dominated cell-type profile | Same | Property of the published tool |
| Near-chance agreement with other callers | Same | Property of the published tool |

**Pre-register the prediction before running H1.** Either answer publishes: you have found a real limitation of a widely adopted tool, or you have shown your mouse adaptation caused it, which strengthens everything else by removing a confound.

You get this for free because you are running DeepScence on both arms anyway.

## 9. Reporting Standard

Every DeepScence number carries four attributes: **coverage** (how many sections), **denoise state**, **anchor** (published `CDKN1A` or re-anchored), **panel** (native or ortholog-mapped, with the mapping rate). One row in the Methods table, not scattered caveats.

**All of Part II must be resolved before the freeze**, because the sender-caller definition is a frozen parameter.

---

# Part III — Step 3: Annotate the Panel (Job A)

## 10. Gene Sets in Human Symbols

Rebuild Tiers A–E from the mouse versions. Owner: biology collaborator, Deliverable 2. One day.

**Tier A — senders (arrest and damage only, no secreted factors):**
`CDKN1A`, `CDKN2A`, `CDKN2B`, `TP53`, `TP53I3`, `GADD45A`, `LMNB1`, `HMGB1`, `HMGB2`, `ATM`, `ATR`, `CHEK1`, `CHEK2`, `H2AX`, `TP53BP1`, `MDM2`

**Tier B — response, seven modules.** Pull the current **human** MSigDB Hallmark collection (the H collection, not MH), pinned by version and date, raw JSON archived in the repo the way `msigdb_mouse_2026.1.Mm` is: `HALLMARK_TNFA_SIGNALING_VIA_NFKB`, `HALLMARK_IL6_JAK_STAT3_SIGNALING`, `HALLMARK_INTERFERON_GAMMA_RESPONSE`, `HALLMARK_E2F_TARGETS`, `HALLMARK_G2M_CHECKPOINT`, `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`, `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY`, plus `secondary_senescence`.

**Tier C — ligand–receptor.** Same pairs, plus one addition: **`CXCL8` (IL-8) exists in human and has no mouse ortholog.** In mouse you used `Cxcl1`/`Cxcl2` as functional analogues. In human you get the real ligand. Declare the asymmetry.

**Tier D — nuisance.** Cell type, density at 25/50/100 µm, k-NN composition, total counts, genes detected, cell and nucleus area, section, segmentation method, distance to tissue boundary, **plus the arm-specific anatomical covariate** (§13, test A6).

**Tier E — controls.** Housekeeping (`ACTB`, `GAPDH`, `RPL13A`, `RPS18`, `TBP`, `PPIA`), a cell-type identity program unrelated to inflammation, 500 random size- and expression-matched sets, **and the negative control probes** (§12).

**Symbol gotchas going the other way:** `Trp53` → `TP53`, `Il6ra` → `IL6R`, `H2ax` → `H2AX`, `Ackr3` → `ACKR3`.

## 11. The Disjointness Gate

```python
A   = tier_A & panel
B_k = tier_B_module_k & panel     # for each of the 7

assert len(A) >= 15
assert all(len(B_k) >= 30 for k in modules)
assert all(len(A & B_k) == 0 for k in modules)
```

If any `A ∩ B_k` is non-empty, remove the overlapping genes **from Tier A**, then re-check `len(A) >= 15`.

**Print the full intersection matrix in the Methods.** A reviewer will check this, and CoreScence was 69% circular with your response modules last time.

**A 5K panel clears this comfortably.** Anything under ~1,000 genes almost certainly does not, which is why the plan is built on Xenium Prime 5K rather than on finding human liver.

---

# Part IV — Step 4: Acquire and Annotate the Human Arm

## 12. What to Download

### 12.1 The dataset

**H1 = human aging lung parenchyma. Primary candidate: GEO `GSE335761`, "Human Cellular Clock of the Aging Lung Parenchyma spatial transcriptomics data."**

**Why lung and not liver.** Liver would be the cleaner contrast — matched tissue means species is the only variable, and every liver-specific asset you built (cell-type markers, the zonation covariate, periportal stratification) carries over. **It is not available at adequate quality.** The audit trail:

| Human liver Xenium option | Why it fails |
|---|---|
| Human liver, Xenium Human Multi-Tissue and Cancer Panel | **377-gene panel** (474 with the cancer add-on), a cell-typing panel not a pathway panel; **median 80 transcripts/cell** against your mouse 448, a 5.6× depth deficit in a project whose central finding is a depth confound at ρ = 0.94; **2 donors**; healthy vs cancer, no aging or fibrosis axis; older chemistry generation |
| `GSE286254` pediatric liver, age-related signatures | Age axis runs the wrong direction; low senescent burden will likely fail Test A3 |
| `GSE332850` HCC, CAR-T trial | Tumour plus therapy exposure, confounded |
| HTAN mCRC liver arm | Untried — check the open-access tier (§16 parallel track) |
| SenNet WashU human liver arm | Not retrievable on this timeline |

**Run the 377-gene check anyway.** Thirty minutes, and it earns a Limitations sentence rather than asserting one:

> A matched-tissue human replication was attempted. The only public human liver Xenium data uses a 377-gene cell-typing panel at 80 median transcripts per cell, which cannot support the disjointness requirement between sender-defining and response-measuring gene sets.

That is a documented data-availability finding, and it tells the field exactly where the gap is.

**Why aging lung is a good arm, not a consolation prize.**

- **It is the human aging axis the mouse arm could not give you.** Your mouse contrast is time-post-surgery with sham controls ageing in parallel. This is chronological human aging, which is what the senescence field actually cares about.
- **Lung is a high-senescent-burden tissue.** Better odds of landing in the Test A3 1–20% prevalence band than healthy liver.
- **Two variables change instead of one, and that is the stronger generality claim.** "The confound reproduces across species *and* tissues" answers the one-dataset objection better than "it reproduces in the same tissue." You give up attribution — you cannot say whether a difference is species or tissue — and you must state that plainly.

**But the panel comes first, and tissue comes second.**

Two human Xenium datasets have now been screened and both fail on panel size:

| Checked | Panel | Verdict |
|---|---|---|
| Human liver, Multi-Tissue and Cancer Panel | 377 genes (474 with add-on), 80 median transcripts/cell, 2 donors | **Fail** |
| `GSE335761` aging lung, `GSM9820340` | **386 targets** (Xenium Human Lung panel + 97 custom, Design ID 77VFDX, v1 chemistry) | **Fail A2 almost certainly** — Tier A may scrape through on the custom add-on, but 7 Tier B modules at 30 genes each will not come out of a 386-gene cell-typing panel |

That is a pattern, not bad luck: **most public human Xenium predates Prime 5K.** So the search is panel-first.

**Search protocol, GEO only.**

1. Query `GPL33762` series titles and summaries for **"Xenium 5K"**, **"Prime 5K"**, **"5K panel"**. These exist — `GSE322974` ("Spatial Transcriptomic Analysis of a Human Bone Marrow Clot Biopsy Using Xenium 5K") is proof.
2. **Confirm on the data, not the title.** Pull one `cell_feature_matrix.h5` per candidate — under 10 MB — and count:

```python
import h5py, numpy as np
f  = h5py.File("<sample>_cell_feature_matrix.h5")
ft = np.array([x.decode() for x in f['matrix/features/feature_type'][:]])
nm = np.array([x.decode() for x in f['matrix/features/name'][:]])
P  = set(nm[ft == 'Gene Expression'])
print(len(P))    # ~5000 -> Prime 5K.  300-500 -> v1, discard
```

3. **Donor count** from the series matrix. ≥3 or the donor bootstrap is gone.
4. **Disjointness** (§11) on the surviving panels.

A dozen candidates screen in an afternoon at ~10 MB each.

**Drop the tissue constraint during the search.** You have already ruled out liver and hit a wall on lung. Prime 5K human deposits are rare enough that you cannot also require a specific organ. Take whatever passes panel and donor screening; aging lung remains the preference, not the requirement.

**One route worth twenty minutes: SenNet data is on GEO.** `GSM9820340`'s title is `SNT228_DHGQ_454_...` — a SenNet ID. SenNet TMCs deposit to GEO directly, so search sample titles for `SNT`. That reaches consortium data through an archive you already know how to use, without the portal.

**The fallback if nothing passes: symmetric `CORE`.**

If no human Prime 5K deposit with ≥3 donors exists, take the best available human panel and define

```
CORE = ortholog_map( panel_H1 ∩ panel_M1 )
```

then **re-run the mouse arm on `CORE` as well**. You would drop from seven Tier B modules to perhaps two or three and lower the per-module threshold from 30 to ~15. That costs real power, but it is symmetric, declarable, and makes the two arms strictly comparable. **Report the mouse arm twice — full panel and `CORE`** — so the cost of the reduction is visible rather than hidden.

**Tier A ∩ Tier B = ∅ never relaxes.** Everything else is negotiable; that is not.

Decide the `CORE` parameters and write them into the pre-registration **before** looking at any human result.

**Platform handle:** `GPL33762` (Xenium In Situ Analyzer: *Homo sapiens*), 1,044 samples across 117 series as of August 2026. Pull the full accession list:

```bash
wget "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL33762&targ=self&view=brief&form=text" \
  -O gpl33762_samples.txt
```

The mouse equivalent, if you later add a second mouse arm, is `GPL33896`.

### 12.2 The exact files

Download the **Xenium Output Bundle**, not the Explorer subset if you can choose — the Explorer subset ships `.zarr.zip` rather than `.h5`/`.parquet` and you would need a reader.

**Take, per sample:**

| File | Why |
|---|---|
| `cell_feature_matrix.h5` (or the MEX directory) | Counts |
| `cells.parquet` (or `cells.csv.gz`) | Coordinates **in microns**, cell IDs, area, transcript counts |
| `cell_boundaries.parquet` | Segmentation polygons; needed for the occupancy grid in C1 |
| `gene_panel.json` | Panel membership, for Job A |
| `experiment.xenium` | Run metadata, chemistry, segmentation method |

**Total: on the order of 100–500 MB per sample.**

**Skip:**

| File | Size | Why not |
|---|---|---|
| `transcripts.parquet` | ~2 GB+ | One row per transcript. You need none of it unless re-segmenting |
| `morphology*.ome.tif`, `morphology_focus/` | 4–11 GB | Imaging, not needed |
| `analysis/`, `analysis.zarr.zip` | large | vendor onboard clustering; you do your own |
| Any `.rds` | up to ~12 GB | R objects you will not load |

`code/fetch_xenium_bundle.sh` already implements streaming extraction of exactly these files. Point it at the GEO per-GSM URLs rather than rewriting.

### 12.3 What you do NOT download

**Senescence labels.** They do not exist for this data. You generate them in §14. That is the whole point of Part II.

**Cell type labels.** Also not provided in usable form. Your `annotate_celltypes.py` / `annotate_pipeline.py` pipeline handles this; run it with human markers.

### 12.4 Also grab

- **Human MSigDB Hallmark collection**, current release, raw JSON, archived in `genesets/`.
- **The negative control probe list** from `gene_panel.json` — pre-designed Prime 5K panels ship 609 negative control codewords, 40 negative control targets, and 21 genomic control probes. You need these for §13 test A7 and they are already in the bundle.

## 13. Acquisition Audit

Master plan Section 8, unchanged, plus two new tests.

| Test | What | Pass condition |
|---|---|---|
| **A1** | Resolution, segmentation, assignment rate | Single-cell confirmed. Record median NN distance (M1: 6.7–9.7 µm) and assignment rate (M1: 88.27%) |
| **A2** | **Panel adequacy and Tier A ∩ Tier B = ∅** | §11. Go/no-go |
| **A3** | Sender prevalence | 1–20%, ≥200 senders, ≥5,000 non-senders |
| **A4** | Ripley's K, sender clustering | Record; compare to M1 (1.11 / 1.26 / 1.56 at the three callers) |
| **A5** | **Matched-decoy contrast test** | \|SMD\| ≤ 0.1 after matching. Go/no-go |
| **A6** | Anatomical confounding | **Not liver zonation.** Lung: the **airway-to-alveolar axis** — distance to the nearest conducting airway and to the nearest large vessel, plus an alveolar-vs-airway-vs-vascular compartment label. In fibrotic lung add distance to the nearest fibroblastic focus or honeycomb region. Build the covariate before fitting anything |
| **A7** | **NEW: negative-control-probe kernel** | Fit the kernel against negative control probe counts vs distance-to-sender. **Must be flat.** If not, there is a spatial technical gradient and every biological kernel is contaminated. Cleaner than housekeeping genes because these have no biology by construction. **[Corrected 2026-08-27, citation audit: "Nobody in this literature reports it" is false — Voyager's Xenium vignette and Ren et al., *Nat Commun* 16 (2025) both compute Moran's I on negative control probes. What is unreported is refitting *the estimand's own estimator* to the controls, i.e. a Lipsitch (2010) negative-control **outcome**.]** |
| **A8** | Cross-arm comparability | Report every cross-arm number on the ortholog-intersected panel **and** on each arm's full panel. Pin the ortholog map by version and archive it |

## 14. Annotate the Cells (Job B)

The step you meant by "the panel isn't marked for senescence."

1. **Cell types first.** Run `annotate_pipeline.py` with human lung markers. Expected compartments: **AT1** (`AGER`, `PDPN`), **AT2** (`SFTPC`, `SFTPB`, `NAPSA`), **club/secretory** (`SCGB1A1`, `SCGB3A2`), **ciliated** (`FOXJ1`, `TPPP3`), **basal** (`KRT5`, `TP63`), **capillary and general endothelial** (`CLDN5`, `PECAM1`, `CA4`, `EDNRB`), **lymphatic endothelial** (`PROX1`, `LYVE1`), **fibroblasts** (`COL1A1`, `PDGFRA`, `PDGFRB`), **smooth muscle** (`ACTA2`, `DES`), **alveolar macrophages** (`MARCO`, `FABP4`), **monocyte-derived macrophages** (`SPP1`, `TREM2`), **T/NK**, **B**, **mast** (`TPSAB1`).

   **AT2 cells are the senescence compartment to watch.** AT2 senescence is the best-established senescent population in aging and fibrotic human lung, so report it as a named receiver type rather than folding it into "epithelial." Your bio collaborator eyeballs fifty called cells against the image before you trust any of it.

2. **Sender scores, four callers, exactly as in mouse:**
   - **DeepScence** — the primary. **Native run, no ortholog remapping** (§8). Resolved denoise state and anchor from Part II.
   - **SenePy** — human cell-type-specific signatures.
   - **`CDKN1A`⁺** — the deliberately crude baseline.
   - **Tier A curated arrest score** — the p90/p95/p99 percentile calls.

3. **Threshold and check prevalence** (test A3), per cell type, not pooled.

4. **Caller agreement**, conditioned on cell type and transcript-depth decile, exactly as `BIO_PHASE3.md` §4.4 does for mouse. **This is half of §8's experiment.**

---

# Part V — Step 5: Run and Compare

## 15. Freeze and Pre-Register

**After Parts I–III are done and M1 is re-run. Before downloading H1.**

```bash
git tag -a phase7-frozen -m "Frozen for human replication, post C1 and C7"
git push --tags
```

**Fixed and not retunable:** λ grid bounds relative to resolution floor and window; window = 99th percentile of distance-to-nearest-sender; kernel families and the selection rule; null battery N1–N8 **with the corrected N3/N4**; the surviving fraction definition; 400 bootstrap replicates over 100 quantile blocks; Tier A–E membership; the four sender callers and their thresholds; the composition-matched rerun protocol at 5 seeds.

**`PREREG_PHASE7.md` contents:** the frozen tag hash; the H1 dataset identifier; the fixed-parameter list; the deviation table (§below); a named primary outcome; a replication criterion written before seeing data; the §18 outcome table; **and the §8 prediction about DeepScence.**

**Declared deviations:**

| Deviation | Reason |
|---|---|
| Human symbols throughout | Species |
| `CXCL8` added to Tier C | No mouse ortholog existed |
| Arm-specific anatomical covariate | Lung has no zonation; the airway-to-alveolar axis replaces it. **Report both arms with and without the anatomical term** so the without-version is strictly comparable |
| **Reduced `CORE` gene set, if invoked** | Best available human panel may cap at ~400 genes. Both arms re-run on `CORE`; mouse reported twice |
| **Tissue differs between arms** | Liver is unavailable at adequate panel depth (§12.1). Species and tissue both vary, so **no cross-arm difference can be attributed to either.** State this in the abstract, not the limitations |
| DeepScence native in H1, remapped in M1 | **This is §8's experiment**, not a flaw |
| DeepScence coverage 2/11 → 11/11 in M1 | C7/D1. Report the two-section values alongside |
| Donor count below 3 in H1 | **Section-level bootstrap only, labelled as such.** No donor bootstrap |
| Corrected N3/N4 | C1. Applies to both arms |

## 16. Run Order

1. **C1** — three torus variants, destructiveness diagnostics, N4 too. Re-run Phase 3 N3/N4.
2. **C7** — D1 coverage overnight, D2 denoise, D3 re-anchor. Regenerate all five caller tables.
3. **Re-run M1** end to end on the corrected pipeline. Regenerate Figures 2 and 4. Record what moved in `reports/CORRECTIONS.md`.
4. **Job A** — human gene sets, disjointness matrix.
5. **Freeze and pre-register.**
6. **Download H1.** Audit A1–A8. Stop on A2 or A5 failure.
7. **Job B** — cell types, four sender callers, prevalence, caller agreement.
8. **Run H1** through the frozen pipeline: naive, N1–N8, controlled fits, kernel families, superposition vs nearest, proximal vs downstream.
9. **Composition-matched reruns**, 5 seeds, both arms.
10. **§8 comparison**, §17 table, figures, write-up.

## 17. Two-Arm Comparison Table

The paper's new centrepiece.

| Quantity | M1 mouse liver (SBR fibrosis) | H1 human aging lung |
|---|---|---|
| Platform / panel | Prime 5K Mouse + 100 custom | Prime 5K Human |
| Sections / donors | 11 / 11, 6 admissible | |
| Cells | 1,834,806 total; 1,036,459 admissible | |
| Median NN distance (µm) | 6.7–9.7 | |
| Transcript assignment rate | 88.27% | |
| Sender prevalence, `tierA_p95` | | |
| **DeepScence coverage** | 2/11 → **11/11** | |
| **DeepScence panel** | ortholog-remapped, 4,845/5,097 | **native** |
| **DeepScence sign vs own gene set** | −0.350 sham / +0.318 SBR | |
| Caller agreement, depth- and type-matched | 0.93–1.22× chance (2-section base) → restate at 11 | |
| Depth ↔ prevalence ρ | 0.94 (SBR, n=6); 0.16 pooled | |
| Naive amplitude (response-sd) | 0.326 | |
| Controlled amplitude (N2+N5+N6) | 0.027 | |
| SF, N2+N5+N6 | 0.082 [−0.099, 0.249] | |
| SF, N5 alone | 0.084 | |
| **SF, N3 corrected (tile / occ / swap)** | *C1 pending* | |
| SF, N2 matched decoy | 0.943 | |
| Composition surrogate share | 66–76% | |
| λ̂ railed at a grid bound | 63% | |
| Poisson identity r² | 0.984 | |
| **Negative-control-probe kernel (must be flat)** | *new* | *new* |
| Detectable bound (response-sd, 80% power) | 0.203 | |

**Note which quantities are geometric predictions** — the Poisson identity and the grid-railing rate — and **should** replicate in any tissue. If one of them does not, a pipeline is broken and finding that out is worth the phase on its own.

## 18. Outcomes

Decide before you look.

**A — Both arms agree, no kernel above bound.** The strongest version. Two species, two tissues, two labs: the distance-dependent SASP response is not separable from technical and compositional confounding at achievable power. With the corrected null and the certification protocol, this is main-track or methods-journal material.

**B — H1 shows a surviving kernel.** More interesting, harder to write. Report both arms, do not drop mouse, hypothesize why: **species, tissue** (these two are confounded by design and you cannot separate them), prevalence, panel, sender-caller behaviour, architecture, donor variance. Name the follow-up that would distinguish them: a human liver arm at adequate panel depth, which does not currently exist publicly.

**C — H1 fails A2 or A5.** A data-availability finding, one honest paragraph, and move to the next candidate in the screen.

**D — The confound structure differs but the null result does not.** The confound is context-specific, so nobody can trust a published characterization from a different dataset. Strengthens the case for running the battery every time.

**E — C1 changes M1's N3 result.** Handle it in the correction ledger and rewrite contribution 3 to what the corrected numbers support. This is why C1 comes first.

**F — DeepScence's instability appears in H1 too.** §8. A real limitation of a published tool, verifiable by anyone.

## 19. Figures

**Figure 1** — unchanged. Synthetic identifiability regime map.

**Figure 2, revised** — naive curve and null battery, now with the **corrected** N3/N4 across three variants, the destructiveness diagnostic, and the **negative-control-probe kernel**.

**Figure 3, revised** — controlled kernel estimates, both arms.

**Figure 4, revised** — CCC tools under the nulls, reading constants from `results/` rather than hardcoded.

**Figure 5, new — Two-arm replication.** (a) Surviving fraction by null, mouse and human side by side. (b) The §17 table as a forest plot. (c) Geometric predictions (Poisson identity, grid-railing) against their expected values in both arms.

**Figure 6, new — DeepScence native vs remapped.** (a) Caller agreement at full coverage, both arms, before and after conditioning on cell type and depth decile. (b) Sign-vs-own-gene-set correlation per arm and condition. (c) Whether the instability transfers. *This carries §8 and it is the most transferable single result in the replication.*

---

# Part VI — Logistics

## 20. Timeline

Ten working days.

| Day | CS Lead | Bio Collaborator |
|---|---|---|
| **1** | **C1**: implement N3-tile, N3-occ, N3-swap, and the N4 equivalents. Destructiveness diagnostics. De-hardcode `make_figure4.py` | **Job A**: human Tiers A–E, pull current human MSigDB Hallmark, archive raw JSON |
| **2** | **C1** re-run of Phase 3 N3/N4. **C7/D2**: attempt the isolated DCA environment. **C7/D3**: implement re-anchoring | **Job A**: intersection matrix, on-panel counts per tier. Confirm §11 gate passes on the human panel |
| **2 (overnight)** | **C7/D1: DeepScence on all 11 M1 sections**, both anchors, resolved denoise state (~9 h, parallel) | — |
| **3** | Regenerate all five caller-agreement tables at 11-section coverage. Compare to the 2-section values. **Gate: did the headline move?** | Interpret the coverage change |
| **4** | **Re-run M1** end to end. Regenerate Figures 2 and 4. Write `reports/CORRECTIONS.md` | Review what the corrections changed; flag anything implausible |
| **5** | **Freeze and pre-register.** Screen `GSE335761` and the four lung alternatives on panel size, donor count, disjointness. Also: 30-min 377-gene liver check for the Limitations sentence; email Ding/Karpova; check HTAN open tier. Download the winner. Audit A1–A8 | Screen candidates on §11 from `gene_panel.json` before download |
| **6** | Handle any format surprises. Verify A2, A5, A7 | **Job B**: cell types with human markers. Eyeball fifty calls against the image |
| **7** | **Job B**: four sender callers, prevalence per cell type (**report AT2 separately**), caller agreement conditioned on type and depth decile | A6: build the airway-to-alveolar covariate — distance to nearest conducting airway and large vessel, compartment label |
| **8** | Run H1 through the frozen pipeline: naive, N1–N8, controlled fits, kernel families | Ligand–receptor plausibility, now with real `CXCL8` |
| **9** | Composition-matched reruns, 5 seeds, both arms. **§8 comparison** | Interpret against the pre-registered prediction |
| **10** | Figures 5 and 6; revise 2, 3, 4. §17 table. Write the corrections and replication sections | Claim audit on every new statement; verify citations |

**Gates.** Day 3: DeepScence at full M1 coverage and the correction ledger records whether the caller-agreement headline moved. Day 5: freeze committed, or you are no longer running a replication. Day 6: A2 and A5 passed, or switch datasets.

## 21. Compute

CPU-bound except DeepScence, which wants a GPU but will run on CPU. The DeepScence full-coverage run is the only heavy step: ~9 hours for M1, parallel across sections.

RunPod section of the master plan applies unchanged: install to `/workspace` on a network volume, bake a Docker image, terminate rather than stop, checkpoint every 50 permutations, work in tmux. 300 GB volume is sufficient for two arms given the download discipline in §12.

## 22. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| **C1 materially changes N3** | Medium | Why it goes first. Rewrite contribution 3 to the corrected number; still a contribution |
| **Full DeepScence coverage changes the caller-agreement headline** | **Medium** | Two sections is a thin base. Report old and new side by side. If agreement rises above chance, restate the motivating claim honestly |
| DCA will not install in any environment | Medium-High | §6 path 3: measure what `denoise=False` costs rather than leaving it unquantified |
| No Prime 5K human deposit with ≥3 donors exists | **Medium-High** | The symmetric `CORE` fallback (§12.1): intersect the best available human panel with the mouse panel, re-run **both** arms on the reduced set, report mouse twice. Declared before the freeze |
| `GSE335761` fails panel or donor screening | **Medium** | GPL33762 includes v1 280–480 gene panels. Four ranked lung alternatives in §12.1; screen five in an afternoon |
| Tissue differs between arms, so differences are unattributable | **High, by construction** | Declared in the abstract. The generality claim is stronger; the attribution claim is gone. Do not pretend otherwise |
| Lung fibrosis architecture destabilizes the nuisance model | Medium | Prefer the aging-lung dataset over the fibrosis ones for the primary. Report with and without the anatomical term |
| A deposit ships `.zarr` rather than `.h5`/`.parquet` | Low | GEO Xenium deposits are normally per-file `.h5` and `.parquet.gz`; if one is not, budget an hour for a reader |
| A2 fails on the human panel | Low | 5K clears it comfortably. If it fails, the Tier B module definitions need trimming, not the panel |
| **You retune after seeing H1** | **Medium-High, the real risk** | The freeze exists for this. Deviations go in the table with a reason; report both versions |

## 23. The Optional Menu

Do these only if time allows. None blocks the replication.

| ID | Issue | Cost |
|---|---|---|
| **C2** | Report the unselected all-315 fit population alongside the selected 160, and make the **amplitude difference** in response-sd the primary quantity rather than the ratio (SF is unstable at small β and not comparable across sample sizes: 0.096 / 0.157 / 0.247 at 100 / 50 / 25 blocks) | Hours, no re-run |
| **C3** | Three different value pairs describe `beta_obs` vs `N3_null` across the README and `make_figure4.py`. Pick one, define it once in Methods | 1 hour |
| **C4** | Three of four CCC tools are reimplementations. **Free fix:** restructure Figure 4 to lead with the COMMOT mechanism result, the only one in published code. **Optional:** run real CellChat v2 in R on one pair, one tile | 0–2 days |
| **C5** | Five open D7 items: the Ma et al. and CellWHISPER sentences need publisher access via the Caltech proxy; the Acosta 2013 Author Correction content is unread; M7/M11 cell counts differ by ≤1,100 (state the QC rule once); two uncited framing sentences | 1 day, mostly library access |
| **C6** | Module B7 `secondary_senescence` shares 14 of 38 genes with the Tier A caller and its sole citation is a meeting abstract. Rebuild without the shared genes, report both, and replace the citation | Hours |
| **M2** | A second **mouse** arm via GEO platform `GPL33896`. Would make it three arms and answer "one dataset" more completely | +5 days |
