# Phase 9 — H1 acquisition audit (A1–A8) and cell annotation (Job B)

**Arm:** H1 = GEO **GSE326743**, 7 normal human spleens, ages 17 / 31 / 32 / 32 / 37 / 57 / 59,
4 M / 3 F, Xenium Prime 5K Human + 100-gene add-on, **5,093-gene panel**, **2,207,593 cells**.
**Date:** 2026-08-27. **Freeze:** `phase8-frozen` = `926439629a07269a32c93f998da0f6e1cd20933c`;
`reports/PREREG_PHASE8.md` committed. This is the first work permitted to read H1 expression data.

**Every number below was produced by a script in this report and is in a file named beside it.**
New code is `code/h1_*.py`; new outputs are `results/phase9_h1/` and `data/processed_h1/`.
Nothing under `results/phase3/`, `figures/`, `genesets/` or `data/processed/` was written.

---

## 0. Verdicts, in the order the audit ran

| Test | Question | Verdict |
|---|---|---|
| **A1** | resolution, segmentation, assignment rate | **PASS**, with one declared discrepancy: H1's median NN distance is **5.45–6.29 µm, below the frozen λ-grid floor of 7.0 µm** on **every** section. Reported, not patched (PREREG §3.1) |
| **A2** | Tier A ∩ Tier B = ∅ on the real panel | **PASS** — re-run of `code/gate_disjointness_human.py` is **byte-identical** to the frozen log |
| **A3** | **sender prevalence — the hard gate** | **PASS.** At the primary call `tierA_p95`, **60 of 76 (78.9 %)** of (section × fine cell type) strata pass, and **4 of 9 merged cell-type families clear all three floors in all 7 sections** — **5 of 9** once the sender call is taken at the same label family the estimator strat​ifies on (§9.4, a declared sensitivity with no threshold tuned). **No failure at the primary call is a prevalence-band failure; every one is a cell-count floor** |
| **A4** | Ripley's K, sender clustering | **PASS (recorded, no threshold).** `tierA_p95` **1.012** (M1 1.136), `cdkn1a_pos` **1.161** (M1 1.263), `senepy_p95` **1.647** (M1 1.556) |
| **A5** | matched-decoy balance, \|SMD\| ≤ 0.1 | **PASS, 35/35** (7 sections × 5 calls). Max \|SMD\| after matching **0.0933** (M1 0.0352) |
| **A6** | red/white-pulp anatomical covariate | **BUILT; VALIDATION IS WEAK AND PARTLY FAILS.** Real spatial structure (Moran's I 0.22–0.49) but the independent follicle-distance check has the **wrong sign in 2 of 7 sections** and never exceeds +0.33. Reported as a limitation |
| **A7** | negative-control kernel, must be flat | **PASS on the pre-registered primary response.** The 40 negative control probes are flat, naive **−0.0118 [−0.0340, +0.0103] p = 0.24**, conditioned **−0.0028 [−0.0263, +0.0206] p = 0.78**. The **pooled** control set is **not** flat naive and **N2 leaves 98 % of it** — an independent replication of M1's P2/P24 |
| **A8** | cross-arm comparability | **PASS.** The pinned panel arithmetic reproduces exactly (5,097 / 5,093 / 4,845 / 2,435 / **2,425**); every panel-dependent quantity is reported on both panels |

**§18 outcome:** neither A2 nor A5 fails, so **outcome C is not triggered**. The human arm proceeds
to Phase 10.

**Reported against interest, the four things that went worst:**
1. **A6's covariate does not validate cleanly** (§7). It is the spleen analogue of liver zonation
   and it is weaker than liver zonation was.
2. **Job B recovers 10–15 of 23 labels per section and the label set is section-dependent**
   (§9). Against the depositors' own 25-label annotation the agreement is **ARI 0.35–0.45**
   (0.55–0.66 on cells neither side calls low-quality). Our fine labels lose the T-cell
   compartment entirely in 4 of 7 sections.
3. **A frozen-pipeline interaction that costs whole compartments their sender call** (§9.4): the
   Tier A call thresholds within **fine** cell types while the estimator strat​ifies receivers on
   **merged** labels, so 0–13.6 % of cells per section are eligible but can never be called
   senders. On H1 that zeroed the T/NK sender set in SPLN43.
4. **The DeepScence caller we froze as PRIMARY is not reproducible across seeds on this arm**
   (§10.5). At full section size, `denoise=False` at `random_state = 1` reproduces the committed
   `random_state = 0` score at **Pearson r = 0.372** and its top-5 % sender set at
   **Jaccard = 0.211**, against an M1 floor of 0.9955 / 0.761. The qualitative §8 findings survive
   the seed; the magnitudes do not.

---

## 1. A3 — SENDER PREVALENCE. The hard gate. **PASS.**

**Rule (Master Plan §8 Test 3 / §13 A3, unchanged):** a stratum passes iff prevalence is in
**1–20 %** with **≥ 200 senders** and **≥ 5,000 non-senders**, **per cell type**, evaluated at
p90 / p95 / p99 and at each caller's own cutoff. Per PREREG §3.11, on H1 A3 is a **reported**
quantity, not an exclusion rule: all 7 sections are analysed either way. For SenePy it is
evaluated on the subset of cells that receive a score, **with both denominators stated**.

```
python3 code/h1_a3_prevalence.py          # -> results/phase9_h1/a3_prevalence_by_type.csv
                                          #    a3_summary_by_caller.csv, a3_types_by_caller.csv
                                          #    a3_pooled_by_section.csv
```

### 1.1 The per-cell-type table — primary call `tierA_p95`, merged label family

The merged label family is the one the estimator strat​ifies on (`sasp_phase3.LABELS = "merged"`),
so this is the table that governs. Cells are PASS / fail; blank = the label is not realised in
that section.

| merged cell type | SPLN07 | SPLN14 | SPLN21 | SPLN24 | SPLN30 | SPLN43 | SPLN44 | sections present | passes |
|---|---|---|---|---|---|---|---|---|---|
| **B cells** | PASS | PASS | PASS | PASS | PASS | PASS | PASS | 7 | **7** |
| **Endothelial** | PASS | PASS | PASS | PASS | PASS | PASS | PASS | 7 | **7** |
| **Mono/Mac/DC** | PASS | PASS | PASS | PASS | PASS | PASS | PASS | 7 | **7** |
| **Stromal** | PASS | PASS | PASS | PASS | PASS | PASS | PASS | 7 | **7** |
| **T/NK cells** (frozen fine-label call) | PASS | **fail** | PASS | PASS | PASS | **fail** | **fail** | 7 | 4 |
| **T/NK cells** (`tierA_merged_p95`, §9.4) | PASS | PASS | PASS | PASS | PASS | PASS | PASS | 7 | **7** |
| Plasma cells | fail | PASS | PASS | PASS | PASS | PASS | PASS | 7 | 6 |
| Neutrophils | PASS | — | PASS | PASS | fail | fail | PASS | 6 | 4 |
| Erythroid cells | PASS | — | fail | — | — | — | PASS | 3 | 2 |
| Megakaryocytes | fail | — | — | — | — | — | fail | 2 | 0 |

**Both T/NK rows are shown deliberately.** The frozen call is computed within *fine* cell types
and the estimator strat​ifies receivers on *merged* labels, so in three sections whole T-cell
clusters are eligible receivers that can never be called senders — 23 senders in 39,158 cells
(SPLN14), **0 in 24,815** (SPLN43), 295 in 33,482 (SPLN44). Applying the **identical percentile
rule** at the merged label family fixes it. **No threshold was tuned**; §9.4 has the mechanism and
the counts. The PI must choose which row is primary before Phase 10; both are in
`a3_prevalence_by_type.csv`. Every other row in this table is identical under either call.

**Why each failure fails.** Under the coherent merged-family call every failure is a
**cell-count** failure, never a prevalence-band failure:

| section | type | n scored | senders | non-senders | prevalence % | band | ≥200 senders | ≥5,000 non-senders |
|---|---|---|---|---|---|---|---|---|
| SPLN07 | Megakaryocytes | 684 | 35 | 649 | 5.12 | ✓ | ✗ | ✗ |
| SPLN07 | Plasma cells | 4,735 | 237 | 4,498 | 5.01 | ✓ | ✓ | ✗ (short by 502) |
| SPLN21 | Erythroid cells | 584 | 30 | 554 | 5.14 | ✓ | ✗ | ✗ |
| SPLN30 | Neutrophils | 1,287 | 65 | 1,222 | 5.05 | ✓ | ✗ | ✗ |
| SPLN43 | Neutrophils | 2,688 | 135 | 2,553 | 5.02 | ✓ | ✗ | ✗ |
| SPLN44 | Megakaryocytes | 761 | 38 | 723 | 4.99 | ✓ | ✗ | ✗ |

The three T/NK failures under the **frozen fine-label** call are the exception and they *are*
band failures — 0.059 %, 0.000 % and 0.881 %, all **below** the 1 % floor — which is exactly the
signature of the §9.4 mechanism rather than of a thin population: 39,158 / 24,815 / 33,482 cells
is not a cell-budget problem.

### 1.2 All callers, all thresholds

`results/phase9_h1/a3_summary_by_caller.csv`, fine labels, scored denominator, 76 strata.

| call | strata passing | % | pass band | pass ≥200 senders | pass ≥5,000 non-senders | median prevalence % | range % |
|---|---|---|---|---|---|---|---|
| **`tierA_p95` (PRIMARY)** | **60 / 76** | **78.9** | 76 | 63 | 60 | 5.00 | 4.99 – 5.14 |
| `tierA_p90` | 59 / 76 | 77.6 | 76 | 69 | 59 | 10.00 | 9.98 – 10.10 |
| `tierA_p99` | 29 / 76 | 38.2 | 74 | 30 | 60 | 1.00 | 0.99 – 1.09 |
| `tierApm_p95__<module>` (D1 sensitivity, all 7) | 60 / 76 | 78.9 | 76 | 63 | 60 | 5.00 | 4.99 – 5.14 |
| `senepy_p95` | 55 / 76 | 72.4 | 65 | 58 | 55 | 5.00 | 4.97 – 5.05 |
| `senepy_p90` | 55 / 76 | 72.4 | 65 | 62 | 55 | 10.00 | 9.51 – 10.07 |
| `senepy_p99` | 24 / 76 | 31.6 | 59 | 30 | 55 | 1.00 | 0.99 – 1.09 |
| **`cdkn1a_pos` (own cutoff)** | **42 / 76** | **55.3** | **64** | 50 | 60 | **2.69** | **0.16 – 13.28** |

**Read this table honestly.** For the three percentile callers the 1–20 % band is satisfied **by
construction** — the caller *defines* prevalence — so for them A3 is purely a **cell-budget** test.
The only caller for which the band is an empirical result is **`CDKN1A`⁺**, and it passes the band
in 64 of 76 strata; **the 12 that fail are all *below* the 1 % floor (0.155 – 0.903 %) and none is
above the 20 % ceiling** — the opposite of the M1 failure mode, where 4 of 11 sections were
excluded for `Cdkn1a`⁺ prevalence *above* 20 %. At p99 both percentile callers fail mostly on the
≥ 200-sender floor, which is arithmetic, not biology.

### 1.3 The `CDKN1A`⁺ prevalence is a depth readout, not an age readout

Pooled `CDKN1A`⁺ prevalence per donor, against that donor's median transcript count and age
(`a3_pooled_by_section.csv` × `a1_sections.csv`):

| section | age | median transcripts | `CDKN1A`⁺ prevalence % |
|---|---|---|---|
| SPLN30 | 57 | 42 | 1.40 |
| SPLN43 | 31 | 50 | 2.02 |
| SPLN14 | 37 | 56 | 1.53 |
| SPLN21 | 32 | 62 | 2.29 |
| SPLN07 | 17 | 75 | 1.95 |
| SPLN44 | 59 | 127 | 7.13 |
| SPLN24 | 32 | 273 | 7.41 |

**Spearman(median depth, prevalence) = +0.821, p = 0.023. Spearman(age, prevalence) = −0.036,
p = 0.94.** The mouse arm found ρ = 0.94 between depth and `Cdkn1a`⁺ prevalence across the SBR
sections; this reproduces it in a **normal human ageing cohort**, where the crude caller's
prevalence varies 5.3× across donors and tracks sequencing depth, not chronological age.
(PREREG §10.3 forbids any age-stratified claim; this is a **negative** statement about the caller,
not an age result.)

### 1.4 SenePy's two denominators (PREREG §3.7, P4)

SenePy has **no spleen hub**; every score is a cross-tissue surrogate, and 7 of 23 spleen labels
get no hub in any tissue. Fraction of eligible cells that receive **any** SenePy score, per section
(`a3_pooled_by_section.csv`): **92.5 % – 100 %** (SPLN14 100 %, SPLN30 100 %, SPLN21 98.2 %,
SPLN24 94.6 %, SPLN44 93.4 %, SPLN43 92.5 %). Both denominators are carried in every row of
`a3_prevalence_by_type.csv` (`*_all` and `*_scored`).

**A discrepancy found and reported, not patched** (`code/h1_senepy_surrogates.py`,
`results/phase9_h1/senepy_surrogates_v1_v2.csv`): the frozen coverage table
`results/phase7_jobA/senepy_spleen_coverage.csv` was built by reading senepy's **v1** human hub
pickle directly, while the mouse arm scored with `senepy.load_hubs(species=...)`, whose default is
**v2**. The two releases differ: 64 vs 65 hubs, 13 keys only in v1, 12 only in v2, hub sizes
differing up to 4×, and **two of the frozen surrogate assignments — `Fibroblastic reticular cells`
and `Fibroblasts`, both mapped to `('skin','fibroblast',1)` — do not exist in v2 at all.**
H1 is scored on **v2**, because matching the mouse arm's estimator is what makes the arms
comparable. The headline of P4 is unchanged on v2 (**16 of 23 labels usable, 7 with no hub of any
tissue: cDC1, cDC2, pDC, lymphatic endothelium, erythroid, megakaryocytes, mesothelial**), but
4 of the 16 get a *different* surrogate hub than the frozen table records.

### 1.5 What A3 does **not** say

- It does not say the prevalence band is informative for percentile callers (§1.2).
- It does not license the fine label set. **Six of the 23 labels are never realised in any of the
  7 sections** — `Lymphatic endothelium`, `Marginal zone B cells`, `Mesothelial cells`,
  `Pericytes`, `Sinusoidal endothelium`, `pDC` — and **three more are realised but never pass A3
  in any section** at `tierA_p95`: `CD8 T cells`, `Megakaryocytes`, `cDC1`. **Nine of 23 labels
  therefore contribute nothing.** A3 qualifies **compartments**, not fine types.
- The A3 fallback screen's warning (`reports/A3_FALLBACK_SCREEN.md` §1) — that an A3 failure would
  likely reproduce on the runner-up panels — **did not need to be acted on.** The screen's own
  arithmetic predicted this: at p = 5 % an H1 cell type needs only 2.4 % of a donor's cells, and
  the qualifying compartments hold 10–21 % each.

---

## 2. A1 — resolution, segmentation, assignment rate

```
python3 code/h1_a1_geometry.py       # -> results/phase9_h1/a1_sections.csv, a1_segmentation.csv,
                                     #    a1_depositor_annotation.csv
python3 code/h1_a1_assignment.py     # -> results/phase9_h1/a1_assignment_rate.csv
```

### 2.1 Resolution — and the one thing that is materially different from M1

| section | age | cells | QC-pass (≥20 counts & ≥5 genes) | median NN, all cells | median NN, QC-passed | median transcripts | median genes |
|---|---|---|---|---|---|---|---|
| SPLN07 | 17 | 249,420 | 227,360 (91.2 %) | 5.89 µm | 6.11 µm | 75 | 65 |
| SPLN14 | 37 | 329,371 | 283,628 (86.1 %) | 5.78 | 6.12 | 56 | 47 |
| SPLN21 | 32 | 220,435 | 196,142 (89.0 %) | 5.82 | 6.10 | 62 | 53 |
| SPLN24 | 32 | 396,173 | 393,202 (99.3 %) | 5.45 | 5.45 | **273** | 212 |
| SPLN30 | 57 | 366,199 | 291,577 (79.6 %) | 5.83 | 6.29 | **42** | 37 |
| SPLN43 | 31 | 331,582 | 270,472 (81.6 %) | 5.58 | 5.98 | 50 | 42 |
| SPLN44 | 59 | 314,413 | 299,897 (95.4 %) | 5.77 | 5.89 | 127 | 102 |
| **all 7** | | **2,207,593** | **1,962,278 (88.9 %)** | **5.45 – 5.89** | **5.45 – 6.29** | **42 – 273** | |

Single-cell resolution is confirmed: median NN 5.45–6.29 µm, 5th–95th percentile 3.85–8.35 µm,
median cell area 41–45 µm², median nucleus count 1, 1.5–4.0 % of cells multinucleate.

**DECLARED DEVIATION, not patched — the λ-grid floor.** PREREG §3.1 freezes
`LAM_LO_FLOOR = 7.0 µm` as "the resolution floor: at or below the median NN distance of every
section", true on M1 (6.74–10.61 µm) for most sections. **On H1 it is above every section's median
NN distance.** §3.1 anticipates exactly this: *"If H1's median NN distance is materially different
from M1's, that is reported as a limitation, not patched."* It is not patched. **Consequence:** the
grid cannot resolve a length constant shorter than 7 µm on H1, i.e. shorter than ~1.2 median cell
spacings, and any H1 λ̂ railing at the low bound must be read as "at or below the grid floor", not
as an estimate. H1 spleen is a denser tissue than mouse liver and this is a real, arm-level
difference.

**Per-donor depth heterogeneity is the other structural fact.** Median transcripts per cell span
**6.5×** across the 7 donors (42 → 273). SPLN24 is an outlier in every direction — deepest, largest,
99.25 % QC-pass, most Leiden clusters, most cell types recovered. Any cross-donor comparison on
this arm must carry depth as a covariate, which the frozen N5 block does (`log_counts`,
`log_genes`).

### 2.2 Segmentation

Three methods in every section (`a1_segmentation.csv`), with a stable mix:
**interior stain (18S) 51.4–69.7 %**, **boundary stain (ATP1A1+CD45+E-Cadherin) 29.0–47.4 %**,
**nucleus expansion 5.0 µm 1.1–2.4 %**. The frozen N5 block carries two segmentation dummies
(`phase3_core.build_blocks`) and the H1 cache records 3 levels per section, so the covariate is
populated identically to M1.

`nucleus_area` is `NaN` for 0.47–1.29 % of cells (no segmented nucleus). `build_blocks` takes
`log1p(nucleus_area)`, so the cache median-fills those and records the count
(`n_nucleus_area_imputed`, 1,121–1,973 per section). Declared.

### 2.3 Assignment rate — a declared acquisition deviation

§12.3 says do **not** download `transcripts.parquet`, and the H1 acquisition did not. A1 asks for
the assignment rate, which no other deposited file carries — there is no `metrics_summary.csv` in
this deposit (checked against the GEO filelist). **Three sections were therefore fetched
(925 MB, `data/raw_h1_transcripts/`) purely for this test.** The method is
`code/assignment_rate.py` verbatim.

| section | transcripts | assigned | assigned % | unassigned % | Q≥20 assigned % | Q≥20 assigned in nucleus % | verdict |
|---|---|---|---|---|---|---|---|
| SPLN07 | 27,431,028 | 25,348,827 | **92.41** | 7.59 | 92.44 | 64.21 | PASS |
| SPLN21 | 21,277,179 | 19,837,302 | **93.23** | 6.77 | 93.18 | 58.03 | PASS |
| SPLN30 | 21,088,882 | 19,858,420 | **94.17** | 5.83 | 94.21 | 76.49 | PASS |

**M1 reference: 88.27 % assigned on section 7259.** H1 is **4–6 points better**, well inside the
plan's 30 %-unassigned bleed-through threshold. The decompressed parquet files are deleted after
use; the `.gz` are kept for reproducibility.

### 2.4 The depositors' QC filter, characterised (PREREG P21)

`annotations.csv.gz` covers **89.3–99.6 %** of the matrix cells. What it drops is a low-content
population, and it is nearly the same population our own QC drops:

| section | annotated | % | median transcripts, annotated | median transcripts, dropped | % of dropped cells that pass **our** QC | annotated cells failing **our** QC |
|---|---|---|---|---|---|---|
| SPLN07 | 239,167 | 95.9 | 79 | **7** | 0.07 % | 11,814 |
| SPLN14 | 302,488 | 91.8 | 60 | **9** | 3.55 % | 19,814 |
| SPLN21 | 209,089 | 94.9 | 65 | **7** | 0.11 % | 12,959 |
| SPLN24 | 394,493 | 99.6 | 274 | **5** | 0.00 % | 1,291 |
| SPLN30 | 326,986 | 89.3 | 46 | **7** | 0.16 % | 35,471 |
| SPLN43 | 298,476 | 90.0 | 57 | **8** | 0.45 % | 28,153 |
| SPLN44 | 306,712 | 97.6 | 131 | **8** | 0.96 % | 6,889 |

**Their filter is a superset-compatible low-count filter: 96.4–100 % of the cells they drop also
fail our ≥20-counts/≥5-genes rule.** The asymmetry runs the other way — they *keep* 1,291–35,471
cells per section that our QC drops, and they label 0.15–0.43 % of the cells they keep
`Low quality` at Level 1. The two label sets are therefore compared on the intersection, and §9
reports agreement both on all joined cells and on cells neither side calls low-quality.

---

## 3. A2 — disjointness gate. **PASS, byte-identical to the frozen log.**

```
python3 code/gate_disjointness_human.py > logs/phase9/a2_gate_rerun.log
diff results/phase7_jobA/gate_disjointness_human.log logs/phase9/a2_gate_rerun.log   # identical
```

Exit 0, `### FROZEN CONFIGURATION VERDICT: PASS ###`, and the re-run log is **identical** to the
committed frozen log after blank-line normalisation. Tier A `A_SENDER_FINAL_strict` = 33/33
on-panel; Tier B = 120 / 71 / 126 / 231 / 113 / 36 / 116, all ≥ 30; A ∩ B_k = 0 for all seven; all
seven per-module sender sets ≥ 15 and disjoint from their own module.
`results/phase7_jobA/` was not modified (`git status` clean for that path).

CoreScence circularity on the human panel, from the same run: **29 of 33 on-panel CoreScence genes
sit in ≥ 1 Tier B module = 88 %**, `B_secondary_senescence` alone accounting for 18. This confirms
**prediction P-iv natively** — 88 % is measured on the native human panel with no ortholog
remapping, so it cannot be an artefact of the mouse mapping.

---

## 4. A4 — Ripley's K, sender clustering

```
python3 code/h1_a4_ripley.py      # -> results/phase9_h1/a4_ripley.csv
```

`code/_ripley.py` verbatim (K at r = 50 µm against the N1 within-cell-type label-permutation null,
10 permutations), with the H1 cache substituted through `code/h1_sec.py`.

| call | H1 mean | H1 range (7 sections) | **M1 mean** | M1 range (11 sections) |
|---|---|---|---|---|
| `tierA_p90` | 1.010 | 0.994 – 1.030 | — | — |
| **`tierA_p95`** | **1.012** | 0.988 – 1.036 | **1.136** | 1.009 – 1.363 |
| `tierA_p99` | 1.028 | 0.943 – 1.084 | — | — |
| **`cdkn1a_pos`** | **1.161** | 1.118 – 1.218 | **1.263** | 1.032 – 1.684 |
| `senepy_p90` | 1.417 | 1.202 – 1.929 | — | — |
| **`senepy_p95`** | **1.647** | 1.295 – 2.415 | **1.556** | 1.272 – 1.936 |
| `senepy_p99` | 2.423 | 1.559 – 4.085 | — | — |

**The ordering of the three callers is identical across the two arms** — SenePy ≫ `CDKN1A`⁺ >
Tier A — which is a property of the callers, not of liver or spleen. The Tier A caller produces
**essentially no excess clustering on H1** (1.012, and 3 of 7 sections below 1.0), against 1.136 on
M1. Since the null is a within-cell-type permutation, this says the H1 Tier A sender set is no more
spatially aggregated than a random draw from the same cells of the same types — a *weaker* starting
point for a distance-to-sender analysis than M1 had, and it belongs in the Phase 10 interpretation.

---

## 5. A5 — matched-decoy contrast. **GO/NO-GO: PASS, 35/35.**

```
python3 code/h1_a5_matchbalance.py    # -> results/phase9_h1/a5_match_balance.csv,
                                      #    a5_smd_by_covariate.csv
```

`phase3_core.match_decoys_section` verbatim — greedy 1-1 nearest-neighbour propensity matching
without replacement, within section and within cell type, caliper 0.25 SD — on the published N2
matching set built by `phase3_core.build_blocks` (log density at 50 µm, log transcript counts, the
arm's anatomical covariate = the A6 pulp axis, and the 20-NN cell-type composition vector).

| call | matches | pass | max \|SMD\| **before** | max \|SMD\| **after** | min match rate |
|---|---|---|---|---|---|
| `tierA_p90` | 7 | 7 | 0.322 | **0.0177** | 0.9999 |
| **`tierA_p95`** | 7 | **7** | 0.483 | **0.0284** | 0.9998 |
| `tierA_p99` | 7 | 7 | 0.824 | **0.0620** | 1.0000 |
| `cdkn1a_pos` | 7 | 7 | 0.721 | **0.0348** | 0.9995 |
| `senepy_p95` | 7 | 7 | **2.011** | **0.0933** | **0.8264** |
| **all** | **35** | **35 (100 %)** | 2.011 | **0.0933** | 0.826 |

**M1 reference: 0.0916 → 0.0352, 100 % of matches pass.** H1 passes with the same margin at the
primary call. Two rows are worth flagging rather than burying:

- **`senepy_p95` is the worst-balanced call on both arms** and is the only one whose match rate
  falls below 1.0 (0.826–0.891 in SPLN43 and SPLN14). Its pre-matching imbalance is **2.011 SD**,
  four times the primary call's. It passes, but it is matching a much more extreme sender set, and
  its post-match 0.0933 is the closest any H1 row comes to the 0.10 threshold.
- `tierA_p99` reaches 0.062 because it has ~2,000 senders to match rather than ~10,000.

---

## 6. A7 — the negative-control kernel. **PASS on the pre-registered primary response.**

```
python3 code/h1_module_fits.py --n-jobs 4     # the biological reference: run_phase3_nulls._section_job
python3 code/h1_a7_controls.py --n-jobs 4     # -> a7_control_probe_{fits,provenance,curves}.csv
python3 code/h1_a7_summarize.py               # -> a7_summary.csv, a7_verdict.txt
```

The estimator is `run_phase3_nulls.SectionFit` / `fit_cell` **verbatim** — 100 µm window, 40-point
λ grid, `MIN_RECEIVERS = 2000`, nested base → +N6 → +N6+N5 designs, 400-replicate spatial block
bootstrap over 100 quantile blocks — with only the response matrix Y substituted. 7 sections ×
2 calls (`tierA_p95`, `cdkn1a_pos`) = **98 control fits per response** and **686 module fits**.
CIs are **section-clustered** over the 7 sections, which PREREG §6 R5 fixes as the only adequately
powered form on this arm.

**Per feature family, separately, as the task requires.** Amplitudes are `β̂ / sd(y)` in
response-SD units; the CI is the section-clustered 95 % interval on the mean.

| response | design | clustered mean [95 % CI] | p | median \|β\|/sd | CI excl. 0 |
|---|---|---|---|---|---|
| **40 negative control probes — PRE-REGISTERED PRIMARY** | naive | **−0.0118 [−0.0340, +0.0103]** | **0.239** | 0.0412 | 0.163 |
| | +N6 | −0.0118 [−0.0341, +0.0105] | 0.242 | 0.0415 | 0.153 |
| | **+N6+N5** | **−0.0028 [−0.0263, +0.0206]** | **0.779** | 0.0346 | 0.133 |
| | N2 matched decoy | −0.0120 [−0.0340, +0.0099] | 0.229 | 0.0417 | 0.153 |
| 609 negative control codewords | naive | +0.0013 [−0.0153, +0.0179] | 0.857 | 0.0362 | 0.114 |
| | +N6+N5 | +0.0013 [−0.0146, +0.0172] | 0.851 | 0.0397 | 0.068 |
| **21 genomic controls** | naive | **−0.0296 [−0.0532, −0.0060]** | **0.0221** | 0.0589 | 0.235 |
| | **N2 matched decoy** | **−0.0289 [−0.0516, −0.0062]** | **0.0207** | 0.0583 | 0.194 |
| | **+N6+N5** | **−0.0021 [−0.0168, +0.0126]** | **0.735** | 0.0386 | 0.102 |
| **pooled control features (`all_controls`)** | naive | **−0.0337 [−0.0530, −0.0145]** | **0.0052** | 0.0518 | 0.255 |
| | **N2 matched decoy** | **−0.0331 [−0.0515, −0.0148]** | **0.0045** | 0.0510 | 0.235 |
| | **+N6+N5** | **−0.0051 [−0.0156, +0.0054]** | **0.280** | 0.0358 | 0.112 |
| `neg_probe_rate` (**not a clean null** — its denominator is an N5 column) | naive | −0.0095 [−0.0271, +0.0082] | 0.238 | 0.0414 | 0.133 |
| **BIOLOGICAL Tier B modules (reference, 686 fits)** | naive | **+0.0392 [+0.0052, +0.0732]** | 0.0304 | 0.0977 | 0.491 |
| | **+N6+N5** | **+0.0147 [+0.0015, +0.0279]** | 0.0344 | 0.0459 | 0.198 |

**Verdict.** The pre-registered primary response — the 40 negative control probes, which
`PREREG_PHASE8_genesets.md` §11 designates the primary technical null — **is flat on H1** under
both the naive and the full N6+N5 design. **A7 passes.**

**Three further results from the same run, all reported rather than buried.**

1. **The pooled control set is NOT flat naive on H1 either, and the N2 matched-decoy design leaves
   98 % of it.** −0.0337 naive → −0.0331 under N2 (a 1.8 % reduction) → −0.0051 under N5
   (85 % removed, p = 0.28). On M1 the same comparison was −0.070 → −0.061 (13 % removed) →
   +0.007. **This is an independent replication of P2 and P24 in a second species, a second
   tissue and a second laboratory.** PREREG §10.1's prohibition — no naive and no N2-only kernel
   may be reported as a distance effect — is now supported on both arms.
2. **The family carrying the non-flatness is different between the arms, and that must be stated.**
   On M1 the pooled response was carried by the 609 codewords (~73 % of control counts,
   0.0428 counts/cell). On H1 the codewords are **62× sparser (0.00069 counts/cell, 0.068 % of
   cells non-zero)** and are flat; the pooled non-flatness is carried by the **21 genomic
   controls**, which are 90 % of H1's control counts (0.0243 of 0.0269 per cell). Ten of 98
   codeword fits have `sd_y` exactly 0 — no codeword count at all among that stratum's receivers —
   and are dropped and counted (`n_dropped_degenerate`) rather than allowed to explode a ratio.
   The **H1 codeword row is therefore a low-power null, not a demonstration of flatness.**
3. **Power.** The smallest amplitude a single H1 A7 fit could resolve is a median CI half-width of
   **0.0775 SD naive / 0.0828 SD conditioned**; the biological modules sit at 0.0977 naive /
   0.0459 conditioned. So one fit cannot separate the conditioned biological amplitude from zero —
   exactly why R5 requires the pooled, section-clustered form. Pooled, the conditioned control
   residual is −0.0051 [−0.0156, **+0.0054**] and the conditioned biological mean is
   **+0.0147** [+0.0015, +0.0279]: **the control interval's upper limit excludes the biological
   mean**, so the residual technical gradient does not account for it.

---

## 7. A6 — the red/white-pulp covariate. Built; **its validation is weak and partly fails.**

```
python3 code/h1_a6_compartments.py SPLN07 ...   # -> data/processed_h1/anatomy_h1_<sec>.csv,
                                                #    results/phase9_h1/a6_summary.csv,
                                                #    a6_validation_by_label_<sec>.csv
```

§13's A6 text is written for lung and is void for this arm. The gene side was frozen before the tag
as the five `genesets/human/D_spleen_*.txt` sets (`code/spec_a6_compartments_human.py`); this is
the half that needed expression.

**Construction, mirroring `phase2_downstream.py`'s D-B block one for one.**
mouse `zon = score(pericentral) − score(periportal)`, z-scored **on hepatocytes**, tertiles.
H1 `pulp = score(D_spleen_red_pulp) − mean(score(follicle), score(tzone))`, z-scored **on all
analysis cells**, tertiles → white_pulp / intermediate / red_pulp. `score_genes` `ctrl_size = 200`,
identical. **Declared difference:** the mouse axis is standardised on the one parenchymal type that
carries it; spleen has no such type — red and white pulp are made of *different* cells — so
standardising on any one type would define the axis out of existence. All five compartment sets are
fully on panel (28 / 13 / 24 / 8 / 18 genes). Also built: `dist_to_boundary_um` (occupancy grid →
closing/fill/opening → EDT, `GRID_UM = 25.0`, verbatim) and `dist_to_follicle_um`
(DBSCAN eps = 30, min_samples = 10 on follicular + germinal-centre B cells → centroids), which
occupies the mouse `dist_to_portal_triad_um` slot.

| section | follicle foci | **V1** corr(axis, dist-to-follicle) | V1s, 20-NN smoothed axis | **V4** Moran's I (20-NN) | **V5** corr(axis, log counts) | median dist-to-boundary |
|---|---|---|---|---|---|---|
| SPLN07 | 113 | +0.026 | +0.048 | 0.290 | +0.071 | 539 µm |
| SPLN14 | 178 | **−0.013** | **−0.022** | 0.248 | −0.024 | 682 |
| SPLN21 | 42 | **+0.172** | **+0.328** | 0.247 | +0.123 | 550 |
| SPLN24 | 177 | +0.002 | +0.004 | **0.492** | −0.105 | 654 |
| SPLN30 | 237 | **−0.052** | **−0.103** | 0.221 | +0.096 | 586 |
| SPLN43 | 102 | +0.142 | +0.257 | 0.278 | −0.085 | 600 |
| SPLN44 | 111 | +0.078 | +0.145 | 0.273 | +0.145 | 675 |

**Honest reading.**
- **The axis has real spatial structure.** Moran's I on the 20-NN graph is **0.22–0.49** in every
  section — this is not a per-cell noise score.
- **The independent validation is weak and inconsistent.** V1 — the exact analogue of the mouse
  check `corr(zonation, dist_to_portal_triad) > 0` — is **negative in 2 of 7 sections** and never
  exceeds +0.17 raw (+0.33 smoothed). The mouse arm's own range on the same check was **−0.085 to
  +0.343** over 11 sections, so H1 is *comparable* to M1 rather than worse — but neither is a
  strong validation, and it should not be presented as one.
- **The failures are where the annotation is worst.** The two negative sections, SPLN30 and
  SPLN14, are the two shallowest (42 and 56 median transcripts) and are the two where DBSCAN found
  the most, smallest follicle foci (237 and 178) — i.e. the follicular-B call fragmented and the
  "distance to nearest follicle" measure degraded. This is a Job B failure propagating into A6, not
  an independent second failure.
- **Depth loading, against interest.** V5 is −0.11 to +0.15, i.e. small but not zero and
  **sign-inconsistent across sections**, which is what a weak score looks like.

**Consequence carried into Phase 10:** the pulp axis occupies the `zonation` slot of the frozen N5
block and the N2 matching set, exactly as liver zonation does. It is used because the design
requires an arm-specific anatomical covariate, and its weakness is a stated limitation on any claim
that H1's anatomy has been conditioned away. `a6_validation_by_label_<sec>.csv` carries the mean
axis value per **depositor Level-3 label** as an external check.

---

## 8. A8 — cross-arm comparability

```
python3 code/h1_a8_crossarm.py    # -> results/phase9_h1/a8_panel_arithmetic.csv,
                                  #    a8_ortho_sender_shift.csv
```

### 8.1 The pinned arithmetic reproduces exactly

The script asserts all five counts and fails loudly if any moves:

```
mouse panel 5097 | human panel 5093 | mouse genes with an MGI map row 4845
  ... whose ortholog is on-panel 2435 | distinct human symbols 2425
```

Map: `genesets/mouse_human_orthologs_MGI.csv` (pinned, MGI 1:1). **Ortholog-intersected panel =
2,425 human symbols**, exactly as `PREREG_PHASE8_genesets.md` fixes it.

### 8.2 Gene sets on the full panel vs the intersection, both arms

| gene set | human, full panel | human, intersected | % | mouse, full panel | mouse, intersected | % |
|---|---|---|---|---|---|---|
| **`A_SENDER_FINAL_strict`** | 33 | **26** | 78.8 | 33 | **27** | 81.8 |
| `B_tnfa_nfkb_proximal` | 120 | 90 | 75.0 | 126 | 88 | 69.8 |
| `B_il6_jak_stat3` | 71 | 59 | 83.1 | 68 | 58 | 85.3 |
| `B_interferon_response` | 126 | 82 | 65.1 | 100 | 80 | 80.0 |
| `B_downstream_arrest` | 231 | 138 | 59.7 | 190 | 137 | 72.1 |
| `B_emt_ecm` | 113 | 87 | 77.0 | 125 | 86 | 68.8 |
| **`B_oxidative_stress`** | 36 | **19** | **52.8** | 31 | **18** | **58.1** |
| `B_secondary_senescence` | 116 | 88 | 75.9 | 108 | 85 | 78.7 |

The 26-vs-27 asymmetry on Tier A is the documented MGI map gap (one mouse symbol maps to a human
symbol that is itself on the intersection, but two mouse symbols collide on one human symbol).
**`B_oxidative_stress` loses nearly half its membership on the intersected panel on both arms** —
it was already the module with the smallest margin over the ≥30 floor, and on the intersection it
is 19 and 18. Any cross-arm oxidative-stress claim on the intersected panel rests on 19 genes.

### 8.3 What the panel restriction does to the H1 sender call

Tier A rescored with the panel cut to the 2,425 intersection (26 of 33 Tier A genes survive), same
percentile rule (`a8_ortho_sender_shift.csv`, 7 sections):

| call | Spearman(score full, score intersected) | sender-set Jaccard | prevalence |
|---|---|---|---|
| `tierA_p90` | 0.848 – 0.934 (mean 0.877) | 0.583 – 0.633 (mean 0.614) | unchanged by construction |
| **`tierA_p95`** | 0.848 – 0.934 (mean **0.877**) | **0.522 – 0.575 (mean 0.554)** | unchanged |
| `tierA_p99` | 0.848 – 0.934 | 0.423 – 0.481 (mean 0.455) | unchanged |

**Removing 7 of 33 Tier A genes changes 45 % of the sender set at p95** while the underlying score
correlates at ρ = 0.88. This is the number to quote whenever a cross-arm sender-based comparison is
made on the intersected panel: the two panels do not select the same cells, and the pre-registered
requirement to report every cross-arm number twice is not a formality.

---

## 9. Job B — cell types, and the cross-check against the depositors

```
bash code/h1_run_annotate_queue.sh SPLN24 SPLN30 ...   # code/h1_annotate.py per section
python3 code/h1_jobB_crosscheck.py                     # -> jobB_crosscheck_*.csv
```

`code/h1_annotate.py` is `code/annotate_pipeline.py` transplanted with **no threshold changed** —
`MIN_MARKERS`, `MIN_DET`, `MIN_Z`, `MIN_MARGIN`, `LOWQ_FRAC`, `RES`, `MIN_COUNTS`, `MIN_GENES`,
`LOWQ_RESCUE_*` are **imported from that file** so they cannot drift. Three declared differences:
the human panel loader; the label set; and the liver hepatocyte-floor assertion replaced by a
**reported** haematopoietic fraction plus this cross-check (the images were deliberately not
downloaded, and the depositors' four-level annotation is a stronger check than §14's eyeball).

**Label set: 23 = the frozen 22 + Plasma cells.** `code/markers_human_spleen.py` is a generated file
under `genesets/.geneset_manifest.json` and was **not edited**; the pre-registered plasma-cell
exception (PREREG **P6** / PI decision D2b — `MIN_MARKERS = 4` waived for the three highly specific
on-panel markers `JCHAIN`, `MZB1`, `XBP1`) is applied in `h1_common.marker_set()`. All 23 labels
clear their marker floor on the H1 panel; none is dropped.

### 9.1 What the pipeline actually recovers

Leiden clusters at the frozen `RES = 1.0`: **13, 14, 14, 14, 18, 19, 20** (M1: 19–27). Fine labels
realised per section: **10–15 of 23**. Mean composition over the 7 sections:

| merged label | mean % | per-section range |
|---|---|---|
| Endothelial | 17.5 | 15.7 – 20.9 |
| B cells | 17.4 | 13.5 – 21.3 |
| Mono/Mac/DC | 14.3 | 11.3 – 18.2 |
| T/NK cells | 13.1 | 9.2 – 18.6 |
| Stromal | 11.9 | 10.6 – 14.6 |
| **Unknown** | **9.7** | **0.2 – 19.8** |
| **Low_quality** | **7.4** | **0.0 – 12.5** |
| Neutrophils | 3.9 | 0.0 – 10.7 |
| Plasma cells | 3.8 | 2.1 – 5.6 |
| Erythroid cells | 1.0 | 0.0 – 4.5 |
| Megakaryocytes | 0.1 | 0.0 – 0.3 |

The five large compartments are stable across donors to within ±3 points. **`Unknown` and
`Low_quality` are not**: `Low_quality` is 0 % in SPLN14 and SPLN30 (no cluster fell below the
frozen 50 %-of-median-counts rule) and 12.5 % in SPLN21, and `Unknown` runs 0.2–19.8 %. The fine
labels are worse: `CD4 T cells` is 0 % in 4 of 7 sections, `NK cells` 0 % in 3, `Neutrophils`
0 % in 1 and 10.7 % in another. **This is the same instability the mouse arm documented
(`BIO_PHASE3.md` §1.1) and the reason `sasp_phase3.LABELS = "merged"` is the default; on H1 it is
larger.**

### 9.2 Against the depositors' four levels

`jobB_crosscheck_scores.csv`, 7 sections. "clean" = cells neither label set calls low-quality or
Unknown.

| our labels | vs depositor level | ARI (all) | ARI (clean) | NMI (clean) |
|---|---|---|---|---|
| `cell_type` (fine) | Level_3 (25 labels) | 0.352 – 0.451, median **0.412** | 0.552 – 0.663, median **0.589** | median 0.560 |
| `cell_type_merged` | Level_3 | 0.369 – 0.454, median 0.388 | 0.438 – 0.579, median **0.522** | median 0.533 |

**Where it agrees** (SPLN21, `jobB_crosscheck_recall.csv`, merged labels): depositor B cells →
our B cells 91.5 %; sinusoidal endothelium → our Endothelial 77.8 %; neutrophils → Neutrophils
77.2 %; monocytes → Mono/Mac/DC 84.6 %; red pulp macrophages → Mono/Mac/DC 71.8 %; vascular smooth
muscle → Stromal 93.5 %; arterial endothelium → Stromal 85.9 %; T cells → T/NK 64.1 %;
gdT/NK/T → T/NK 74.0 %.

**Where it fails, plainly.**
- **The fine label set loses the T-cell compartment.** At the fine level, depositor `T cells`
  (11,032 cells in SPLN21) are **63.7 % `Unknown`** and `gdT cells, NK cells, and T cells` are
  59.3 % `Unknown`, because `CD4 T cells` and `CD8 T cells` tie inside the frozen `MIN_MARGIN =
  0.20`. The **merged** label set recovers them (`Unknown` falls to 0.05 %), which is why the
  merged family is the operative one.
- **Our `Smooth muscle / capsule` is largely arterial endothelium** in the depositors' view:
  83.9 % of their `Arterial endothelial cells (KDR+ CD34+)` land in it, and only 40.6 % of our
  label is what we named it after.
- **Our `Endothelial cells` is really sinusoidal endothelium** (78.2 % of the label), and the
  distinct `Sinusoidal endothelium` label is never realised — the two are not separable at
  `RES = 1.0` on this panel.
- **Depositor `Smooth muscle` is 24.8 % `Low_quality` to us**, and `Erythroid` only 33.2 % ours.
- Of the depositors' 25 Level-3 labels, the ones we never recover as a modal assignment include
  the four reticular-cell populations, white-pulp macrophages, tingible-body macrophages,
  follicular dendritic cells, promyelocytes and metamyelocytes.

**Interpretation for the paper.** Job B recovers the five major spleen compartments reliably and
does not resolve the fine architecture the depositors annotated with images and a purpose-built
taxonomy. Every H1 result stratified by cell type is a result about **compartments**, and the
marginal-zone-specific claims PREREG P7 already restricted to exploratory status cannot be
supported at all: `Marginal zone B cells` is never realised as a fine label in any of the 7
sections.

### 9.3 The four sender callers

| caller | status on H1 | note |
|---|---|---|
| **Tier A curated, `A_SENDER_FINAL_strict` (33 genes)** | 33/33 on panel, all 7 sections | frozen strict-33 **primary**; the seven `A_sender_for_<module>` sets run as the pre-registered D1 sensitivity and give an A3 result identical to the primary to 3 decimal places |
| **`CDKN1A`⁺** | all 7 sections | prevalence 0.16–13.3 % per stratum; depth-driven (§1.3) |
| **SenePy** | all 7 sections, **cross-tissue surrogates only** | **no spleen hub exists.** 16 of 23 labels get a surrogate; **7 get nothing**: cDC1, cDC2, pDC, lymphatic endothelium, erythroid, megakaryocytes, mesothelial. 92.5–100 % of eligible cells scored. Run on **v2** hubs — see §1.4 |
| **DeepScence** | see §10 | **native, no ortholog remapping** — §8's free experiment |

### 9.4 A pipeline interaction that costs whole compartments their sender call

**Found in Phase 9, reported, and not silently fixed.** The frozen Tier A call thresholds within
**fine** cell types (`phase2_downstream.py`, reproduced verbatim in `code/h1_callers.py`), while the
frozen receiver stratification and the sender-eligibility mask use **merged** labels
(`sasp_phase3.LABELS = "merged"`, `Sec.sender_mask`). A cluster that is `Unknown` at the fine level
but resolved at the merged level is therefore an **eligible cell that can never be called a
sender**. Cells in that state, per section:

| section | fine-`Unknown` but merged-assigned | % of QC-passed cells |
|---|---|---|
| SPLN07 | 0 | 0.00 |
| SPLN14 | 38,701 | **13.64** |
| SPLN21 | 20,155 | **10.28** |
| SPLN24 | 7,386 | 1.88 |
| SPLN30 | 0 | 0.00 |
| SPLN43 | 24,815 | **9.17** |
| SPLN44 | 27,588 | **9.20** |

Effect at `tierA_p95` on the merged `T/NK cells` stratum: **23 senders in 39,158 cells (SPLN14),
0 in 24,815 (SPLN43), 295 in 33,482 (SPLN44)** — i.e. T/NK fails A3 in 3 of 7 sections for a purely
mechanical reason. Applying the **identical percentile rule at the merged label family**
(`tierA_merged_p95`, a **declared sensitivity, not a new threshold** — no parameter is tuned) makes
T/NK pass in **7 of 7**. Both rows are in `a3_prevalence_by_type.csv` and both are shown in
§1.1. The same interaction exists on M1 (`BIO_PHASE3.md` §1.1
documents 7001 losing its stellate compartment to `Unknown`) but is smaller there. **Recommendation
for the PI, flagged not taken:** decide before Phase 10 which label family the sender call is
defined on, and record it as a deviation either way.

---

## 10. DeepScence, native — §8's free experiment

```
bash code/h1_run_deepscence_supervisor.sh   # code/h1_deepscence.py, 2 sections at a time
python3 code/h1_deepscence_anchor.py    # -> results/phase9_h1/deepscence_anchor_h1.csv
python3 code/h1_caller_agreement.py     # -> results/phase9_h1/caller_*.csv
```

**Coverage: 7 / 7 sections, 1,962,278 cells, no subsampling.** Settings are the frozen ones
(PREREG §3.9) — `denoise=False` (the chosen primary), `random_state=0`, published `CDKN1A`
anchor, ≥ 20 counts/cell — with the one deliberate difference that makes this the experiment:
**the panel is native human, all 5,093 genes, no ortholog remapping.** On M1 the same caller runs
on 4,845 ortholog-remapped genes. Wall time 16.5–41 min per section.

The four attributes every DeepScence number carries (§9 reporting standard):
**coverage 7/7 · denoise=False · published CDKN1A anchor · native human panel, no remap.**

> **Read §10.1–10.3 with §10.5.** Every DeepScence number below is at `random_state = 0`. On a
> full H1 section a different seed reproduces the score at Pearson r = 0.372 and the top-5 %
> sender set at Jaccard = 0.211 — far below the M1 floor of 0.9955 / 0.761. The *directions*
> below survive the seed; the *magnitudes* do not.

### 10.1 The five §8 predictions that Phase 9 can test

| # | prediction | falsifier | **H1 result** | verdict |
|---|---|---|---|---|
| **P-i** | score correlates **positively** with transcript counts in ≥ 5 of 7 sections, magnitude comparable to M1's +0.29…+0.56 | null or negative in ≥ 3 of 7 | **positive in 7 of 7: +0.182, +0.183, +0.229, +0.296, +0.302, +0.314, +0.354** | **CONFIRMED** — but at the *bottom* of M1's range (M1 median 0.385, max 0.558) |
| **P-ii** | the published `CDKN1A` anchor is weak, unstable or inverted — depth-partialled fold-split sign stability < 0.90 — in ≥ 1 of 7 | all 7 sections ≥ 0.90 | **stability = 1.000 in all 7 sections**; depth-partialled ρ = **+0.191 … +0.254**, positive everywhere; `CDKN1A` ranks 1st–7th of the 33 on-panel CoreScence genes | **FALSIFIED** |
| **P-iii** | sign-invariant (\|score\|) depth- and type-matched agreement with Tier A is above chance, pooled, ratio > 1.10 | pooled ≤ 1.05, or below chance | pooled **1.102**, z = 5.67, above chance in **5 of 7** sections | **CONFIRMED, but only just** — 1.102 against a 1.10 threshold |
| **P-iv** | CoreScence circularity with the frozen Tier B modules is ≈ 88 % natively | materially below 88 % | **29 of 33 on-panel CoreScence genes are in ≥ 1 Tier B module = 88.0 %**, measured on the native human panel | **CONFIRMED** |
| **P-v** | within cell type, DeepScence is bottom-of-the-depth-distribution selecting (Q5/Q1 < 1) in ≥ 5 of 7 | Q5/Q1 > 1 in ≥ 3 of 7 | **Q5/Q1 < 1 in 6 of 7**: 0.244, 0.404, 0.416, 0.457, 0.496, 0.778, and 1.170 (SPLN43) | **CONFIRMED** |

**Verdict on §8: four of five confirmed, one falsified, and the falsified one is the one that
matters most for the tool's defence.** The depth loading, the bottom-of-depth selection and the
circularity are **properties of the published tool on its own species and its own panel**, not of
our mouse adaptation. The **sign anchor, however, is not broken on H1**: the published `CDKN1A`
anchor decides the polarity the same way in 20 of 20 random folds in every section, and the
depth-partialled correlation is positive everywhere. **The M1 polarity problem (P12, P13) is
therefore at least partly a property of the ortholog-remapped mouse adaptation.** That is stated
plainly because it removes a confound from *our* result, not from the tool's.

**P-ii in full, because a falsified prediction deserves its numbers.**
`results/phase9_h1/deepscence_anchor_h1.csv`, produced by `code/h1_deepscence_anchor.py`, which
imports `partial_spearman` from the mouse producer `code/deepscence_reanchor.py` so the definition
cannot drift. All correlations are **depth-partialled** — the rank of log transcript counts is
linearly removed from both variables — because on a 5K panel every anchor candidate is partly a
detection-rate readout. Stability = the fraction of 20 random equal folds whose sign of the
partialled correlation matches the whole-section sign.

| section | ρ(score, counts) | `CDKN1A` ρ partial | **`CDKN1A` stability** | 8-gene proliferation ρ partial | prolif. stability | `LMNB1` ρ partial | `CDKN1A` rank among the 33 on-panel CoreScence genes |
|---|---|---|---|---|---|---|---|
| SPLN07 | +0.319 | +0.230 | **1.00** | +0.025 | 1.00 | +0.203 | 2 |
| SPLN14 | +0.204 | +0.205 | **1.00** | +0.039 | 1.00 | +0.213 | 4 |
| SPLN21 | +0.312 | +0.236 | **1.00** | +0.010 | **0.75** | +0.233 | 1 |
| SPLN24 | +0.229 | +0.191 | **1.00** | +0.031 | 1.00 | +0.204 | 5 |
| SPLN30 | +0.193 | +0.222 | **1.00** | +0.031 | 1.00 | +0.222 | 7 |
| SPLN43 | +0.317 | +0.220 | **1.00** | +0.045 | 1.00 | +0.227 | 6 |
| SPLN44 | +0.372 | +0.254 | **1.00** | +0.030 | 1.00 | +0.245 | 2 |

The published `CDKN1A` anchor decides the polarity the same way in **20 of 20 folds in every
section**, its depth-partialled correlation with the score is positive everywhere, and `CDKN1A`
sits in the top 7 of the 33 on-panel CoreScence genes by correlation with the score in all seven.
By contrast the **D3 primary alternative anchor** — the 8-gene proliferation set, chosen on M1
precisely because the published one misbehaved — carries almost no signal here (ρ +0.010 … +0.045)
and is itself unstable in one section (0.75). **On the human arm the published anchor is the
better one.** `LMNB1` tracks `CDKN1A` closely but remains a secondary only: PREREG **P12** excludes
it because it is a member of `B_downstream_arrest` and `B_secondary_senescence`.

### 10.2 The caller's depth loading, both arms, same estimator

Spearman ρ of each score against per-cell transcript counts, `caller_technical_loading.csv`
(H1) and `results/phase3/caller_technical_loading_11sections.csv` (M1):

| caller | **H1** min – median – max (7) | **M1** min – median – max (11) |
|---|---|---|
| `tierA_score` | −0.174 – **−0.023** – +0.098 | −0.061 – **+0.056** – +0.125 |
| `cdkn1a_counts` | +0.056 – **+0.084** – +0.165 | +0.015 – **+0.192** – +0.481 |
| `senepy_score` | +0.040 – **+0.234** – +0.332 | +0.097 – **+0.168** – +0.513 |
| `deepscence_score` | +0.182 – **+0.296** – +0.354 | −0.350 – **+0.385** – +0.558 |

Within-cell-type depth-quintile enrichment, Q5/Q1 (`caller_within_type_depth_bias.csv`):

| caller | **H1** range (7) | **M1** median (11) |
|---|---|---|
| `tierA_score` | 0.200 – 0.397 | 0.230 |
| `deepscence_score` | 0.244 – 1.170 | 0.607 |
| `cdkn1a_counts` | 2.978 – 8.337 | 7.161 |
| `senepy_score` | **28.5 – 224.7** | **24.1** |

**The four callers occupy the same four positions on the depth axis in both arms.** Tier A and
DeepScence select from the *bottom* of the within-type depth distribution; `CDKN1A`⁺ and SenePy
from the *top*, SenePy overwhelmingly so — its top-5 % call is 28–225× enriched in the deepest
within-type quintile relative to the shallowest, on H1 even more extremely than on M1. **No
SenePy-called H1 quantity can be reported without this number beside it**, alongside the two
caveats P4 and P8 already impose.

### 10.3 Caller agreement, conditioned on cell type and transcript-depth decile

Primary rule (PREREG §3.7 b): top-5 % recomputed inside each (cell type × within-type transcript
count decile), strata ≥ 50 cells, 10 deciles. Pooled as Σ`n_both` / Σ`exp_both_stratified`, z from
the pooled sd, exactly as PREREG P1 defines it. `caller_agreement_pooled.csv`,
`caller_agreement_matched_significance.csv`.

| pair | **H1 pooled** | z | above chance in | per-section ratios | **M1 (11 sections, post-C6)** |
|---|---|---|---|---|---|
| Tier A × DeepScence | **1.602** | 38.9 | 6 / 7 | 0.97 1.17 1.31 1.31 1.53 1.79 2.82 | 1.288 |
| Tier A × `CDKN1A`⁺ | **1.081** | 3.05 | 5 / 7 | 0.94 0.99 1.09 1.11 1.13 1.22 1.33 | 1.471 |
| **Tier A × SenePy** | **0.874** | **−7.96** | **0 / 7** | 0.79 0.80 0.86 0.87 0.90 0.92 0.98 | 0.972 (n.s.) |
| Tier A × \|DeepScence\| (sign-invariant) | 1.102 | 5.67 | 5 / 7 | 0.60 0.62 1.01 1.12 1.15 1.40 2.26 | 1.285 |
| **DeepScence × `CDKN1A`⁺ (CIRCULAR — excluded from every pooled claim)** | **6.436** | 204.8 | **7 / 7** | 3.46 5.99 6.42 7.38 8.13 9.01 10.12 | median 1.071, pooled 1.255 |
| SenePy × `CDKN1A`⁺ | 1.505 | 18.7 | 7 / 7 | — | — |
| SenePy × DeepScence | 1.093 | 5.89 | 6 / 7 | — | — |

**Three things follow, and one of them settles an open item in the pre-registration.**

1. **The caller-independence claim is dead on this arm too.** No pair is at chance. PREREG §10.2
   and P1 already struck the sentence "their top-5 % calls overlap at 0.93–1.22× of chance … i.e.
   they are statistically independent"; H1 confirms it independently, with Tier A × DeepScence at
   **1.602** — higher than M1's 1.288.
2. **PREREG open item 7 is answered.** That item asks whether the "one pair sits *below* chance in
   every section" sentence survives, because under the post-C6 M1 configuration Tier A × SenePy
   weakened to 0.972, n.s. **On H1 it is 0.874, z = −7.96, and it is below chance in 7 of 7
   sections.** The sentence survives on the human arm and is load-bearing for the same reason it
   always was: a pair that is reliably *anti*-concordant is sharing a technical variable with
   opposite loadings (Tier A selects shallow cells, SenePy deep ones — §10.2), not a latent
   biological state.
3. **The DeepScence–`CDKN1A`⁺ circularity is five times worse natively than it is remapped.**
   Pooled **6.436** on H1 against **1.255** on M1, above chance in 7 of 7 sections and never below
   3.46. DeepScence anchors its sign on `CDKN1A` and its gene set contains `CDKN1A`; on the native
   human panel that circularity is far stronger than the ortholog-remapped mouse run suggested.
   **This pair is excluded from every pooled number** (`BIO_PHASE3.md` §4.4) and is reported alone.
   PREREG §10.7 forbids quoting "1.51–2.85×"; the H1 figure to quote is the pooled 6.436 with its
   7-section range 3.46–10.12.

### 10.4 P-vi and P-vii — `denoise=True` on H1

The isolated DCA environment was rebuilt for this arm (`DCA_ENV_ROOT=/tmp/dca_env bash code/setup_dca_env.sh`, ~4 min; DCA 0.3.4 under TensorFlow 2.4.4 / Keras 2.4.3 in a
CPython 3.8.19 venv). **The main pinned 3.11 stack was not modified and has no TensorFlow.**
TensorFlow 2.4 again could not use this box's GPU — it wants CUDA 11 and the box has an RTX PRO
4500 Blackwell with CUDA 12 libraries — so it ran on CPU. Producer: `code/h1_deepscence_dca.py`
(the H1 analogue of `code/run_deepscence_dca.py`, differing only in the native-panel loader and
the output path). Design copied from the mouse D2 run: **one fixed 20,000-cell subsample**
(subsampling seed 12345, independent of `--seed`, so every run sees the same cells), three
DeepScence seeds, nothing else changed, plus the `denoise=False` companion on the identical cells
as the seed-to-seed floor. Analysis: `code/h1_d2_analyse.py` →
`results/phase9_h1/d2_depth.csv`, `d2_stability.csv`.

**P-vi — the direction does NOT reproduce on H1, on the one section tested.**

| scope | n | ρ(score, counts) `denoise=False` | ρ `denoise=True` | Δρ | ratio | sender-set Jaccard, False vs True |
|---|---|---|---|---|---|---|
| **SPLN21, full section** | 196,142 | **+0.3122** | **+0.1017** | **−0.2104** | **0.33×** | **0.016** |
| SPLN21, 20,000-cell panel | 20,000 | +0.2697 | −0.0972 | −0.3669 | −0.36× | 0.030 |
| **M1 reference** (3 full sections) | 75k–115k | 0.3891 / 0.3176 / 0.4096 | 0.6404 / 0.5314 / 0.5419 | **+0.13 … +0.25** | **1.32–1.67×** | 0.118–0.280 |

On M1, `denoise=True` **raised** the caller's depth loading by 1.32–1.67× on three of three
sections, which is what refuted §4 (D-b)'s premise (PREREG **P29**). **On H1 it lowers it, by a
factor of three, on the one section run.** Stated precisely, because the prediction has a numeric
falsifier: **P-vi's falsifier is Δρ ≤ 0 in ≥ 5 of 7 sections and Phase 9 ran 1 of 7.** So P-vi is
**not formally falsified** — it is contradicted on the only section tested, and the remaining six
are a Phase-10 item. What is already established either way is that **denoising is not
depth-neutral on either arm**: it moves the caller's depth loading by a factor of 1.3–1.7 up on
mouse and 3× down on human, and it changes **98.4 %** of the H1 sender set (Jaccard 0.016).
PREREG **P29**'s conclusion — "whatever DCA contributes here, it is not depth normalisation" —
survives; the *sign* of what it contributes does not transfer across arms and the paper must say
so rather than generalising the mouse direction.

**P-vii — the `denoise=True` seed instability recurs, and so does something worse.**

| config | seeds | Pearson r | top-5 % Jaccard | cells changing status (of 20,000) |
|---|---|---|---|---|
| `denoise=True` | 0 vs 1 | 0.3125 | **0.0444** | 1,830 |
| `denoise=True` | 0 vs 2 | 0.2848 | **0.0510** | 1,806 |
| `denoise=True` | 1 vs 2 | 0.1710 | **0.0482** | 1,816 |
| `denoise=False` | 0 vs 1 | 0.3829 | **0.2006** | 1,327 |
| `denoise=False` | 0 vs 2 | 0.2565 | **0.1788** | 1,392 |
| `denoise=False` | 1 vs 2 | 0.7128 | **0.3683** | 921 |
| **M1 reference** | `False` 0 vs 1 | **0.9955** | **0.761** | 272 |
| **M1 reference** | `True` 0 vs 1 | 0.5703 | 0.000 | 2,000 |

**P-vii is CONFIRMED** — all three `denoise=True` seed pairs are below the 0.30 Jaccard the
prediction names, against a falsifier of "all three agree at ≥ 0.60".

**But the honest headline is bigger than P-vii and runs against our own configuration.** On the
identical 20,000 cells, the **frozen primary `denoise=False` configuration is itself seed-unstable
on the native human panel**: Pearson r 0.257–0.713 and top-5 % Jaccard 0.179–0.368, against M1's
floor of r = 0.9955 and Jaccard = 0.761 at the same panel size and the same design. **The
seed-to-seed floor that every M1 D2 effect was read against does not hold on H1.**

Two things this does and does not license. It **does** mean any DeepScence-derived H1 number is
conditional on `random_state = 0` and that the arm's DeepScence results carry a reproducibility
caveat the mouse arm's do not. It does **not** yet mean the seven full-section runs in §10.1–10.3
are unstable: those are 196k–393k cells, 10–20× the panel, and instability generally falls with n.
A direct full-section seed check is therefore run and reported in §10.5 rather than inferred.

### 10.5 The full-section seed check — and the caveat it puts on §10.1–10.3

Because §10.4's 20,000-cell panel could have been a small-sample artefact, the frozen primary
configuration was re-run on a **full** section at a different seed: `denoise=False`,
`random_state=1`, SPLN21, all 196,142 cells, nothing else changed
(`code/h1_deepscence_dca.py --denoise-false --seed 1 SPLN21`,
`data/processed_h1/deepscence_h1_nodn_seed1_SPLN21.csv`).

| comparison | n | Pearson r | Spearman | top-5 % Jaccard | cells changing status |
|---|---|---|---|---|---|
| **H1 SPLN21, `denoise=False`, seed 0 vs seed 1, FULL section** | 196,142 | **0.372** | 0.416 | **0.211** | 12,779 |
| H1 SPLN21, same, 20,000-cell panel | 20,000 | 0.383 | 0.328 | 0.201 | 1,327 |
| **M1 reference, `denoise=False`, seed 0 vs 1, 20,000-cell panel** | 20,000 | **0.9955** | — | **0.761** | 272 |
| M1 determinism control, *same* seed, full section | 75k / 115k | **0.9999991 / 0.9999999** | — | — | 24 / 2 |

**It is not a small-sample artefact.** At full section size the frozen primary DeepScence
configuration reproduces at **r = 0.372** and **top-5 % Jaccard = 0.211** across seeds on H1,
against **r = 0.9955 / Jaccard = 0.761** on M1 at the same panel size. **Stated against interest:
the seven H1 DeepScence scores in §10.1–10.3 are conditional on `random_state = 0`, and a
different seed would move the identity of roughly 80 % of the top-5 % sender set.**

**What survives the seed and what does not.** Re-deriving the §10 statistics on the seed-1 score
for SPLN21:

| quantity | seed 0 | seed 1 | conclusion |
|---|---|---|---|
| ρ(score, transcript counts) — **P-i** | +0.312 | +0.231 | **direction and rough magnitude survive** |
| within-type Q5/Q1 — **P-v** | 0.496 | 0.292 | **direction survives** (both ≪ 1) |
| Tier A × DeepScence matched ratio | 1.174 (z 3.5) | **2.967 (z 40.0)** | direction survives; **magnitude moves 2.5×** |

So the **qualitative** §8 findings — positive depth loading, bottom-of-depth-distribution
selection, above-chance dependence with Tier A — are seed-robust. **The magnitudes in §10.2 and
§10.3 are not**, and every H1 DeepScence number in this report must be read as "at
`random_state = 0`", with the seed-to-seed spread above as its error bar. This is a stronger
version of PREREG **P26** (which restricted the caveat to `denoise=True`): **on the native human
panel, the caller we froze as PRIMARY is itself not reliably reproducible across seeds.**

---

## 11. Deviations and discrepancies found in Phase 9

Numbered `H*` so they cannot be confused with the gene-set `D*` rows or the Phase-8 `P*` rows.

| # | Deviation / discrepancy | Evidence, and what was done |
|---|---|---|
| **H1** | **H1's median NN distance (5.45–6.29 µm) is below the frozen λ-grid floor of 7.0 µm in every section.** | `a1_sections.csv`. PREREG §3.1 pre-registers exactly this contingency and forbids patching it. **Not patched.** Consequence: no H1 length constant below 7 µm is resolvable, and low-bound railing must be read as "at or below the grid floor". §2.1 |
| **H2** | **`transcripts.parquet` was downloaded for 3 sections, against §12.3.** | A1 asks for the assignment rate and no other deposited file carries it (no `metrics_summary.csv` in this deposit). 925 MB to `data/raw_h1_transcripts/`; the decompressed parquet is deleted after use. Nothing else reads these files. §2.3 |
| **H3** | **The plasma-cell `MIN_MARKERS` exception is applied in code, not in the frozen marker file.** | PREREG **P6** / PI decision D2b admits `Plasma cells` on 3 markers. `code/markers_human_spleen.py` is under `genesets/.geneset_manifest.json` and was **not edited**; the exception lives in `h1_common.marker_set()`. Label set = 23. §9 |
| **H4** | **The frozen SenePy coverage table was built on senepy hub release v1; the mouse arm scored on v2.** | `code/senepy_coverage_human.py` reads `6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle` directly; `phase2_downstream.py` calls `senepy.load_hubs(...)`, default `sig_version='v2'`. 13 v1 keys absent from v2, 12 new, hub sizes up to 4× different, and **two frozen surrogate assignments do not exist in v2**. H1 is scored on **v2** to match the mouse estimator. Both tables written to `senepy_surrogates_v1_v2.csv`. P4's headline is unchanged. §1.4 |
| **H5** | **The Tier A sender call is defined on FINE cell types while the estimator strat​ifies on MERGED labels.** | 0–13.6 % of cells per section are eligible-but-uncallable; T/NK got 0 senders in 24,815 cells in SPLN43. Both the frozen fine-label call and an identical-rule merged-label call are reported; **no threshold was tuned**. A PI decision is needed before Phase 10. §9.4 |
| **H6** | **`nucleus_area` is NaN for 0.47–1.29 % of cells and is median-filled in the H1 cache.** | `build_blocks` takes `log1p(nucleus_area)`; NaN would poison the N5 block. Count recorded per section as `n_nucleus_area_imputed` (1,121–1,973). §2.2 |
| **H7** | **Three DeepScence run-metadata JSONs (SPLN21, SPLN24, SPLN30) are truncated.** | A numpy-float serialisation bug in `code/h1_deepscence.py` raised inside `json.dump` **after** the score CSV was written. The **scores are complete and valid**; only DeepScence's internal `log` (node choice, `reverse` flag, per-node CoreScence metrics) was lost for those three. Fixed (`default=str`); the four later sections have valid metadata. Recoverable only by re-running (~30 min per section) and **not re-run**. §10 |
| **H8** | **Three concurrent DeepScence processes OOM-kill each other at this 57.7 GB cgroup ceiling; two is the maximum.** | Recorded because it is a reproducibility fact, not a mishap: at 3 × ~13 GB anonymous plus page cache the kernel killed the newest processes. A supervisor holding 2 concurrent runs completed all 7 sections without subsampling. The first headroom check budgeted against `memory.current`, which counts reclaimable page cache and wrongly refused sections that fit; it now budgets against `anon` in `memory.stat`. §12 |
| **H10** | **The frozen PRIMARY DeepScence configuration is not reproducible across seeds on H1**, at full section size. | §10.5. `denoise=False`, `random_state` 0 vs 1, all 196,142 cells of SPLN21, nothing else changed: Pearson r **0.372**, top-5 % Jaccard **0.211**, 12,779 cells changing status — against an M1 floor of r 0.9955 / Jaccard 0.761 at the same panel size and a same-seed determinism control of r = 0.9999991. Every H1 DeepScence number in this report is "at `random_state = 0`". Directions survive the seed; magnitudes do not. This extends PREREG **P26**, which restricted the caveat to `denoise=True` |
| **H11** | **DeepScence's own binary cutoff (`binarize=True`) was NOT run**, so A3's "each caller's own cutoff" is satisfied for DeepScence by the p90/p95/p99 percentile calls only. | `binarize=True` costs 50 extra permuted forward passes per section and is not in the frozen settings list (PREREG §3.9), so the M1 arm did not run it either and there would be no like-for-like comparison. Declared rather than silently substituted |
| **H9** | **The section-to-section heterogeneity of this deposit is large and is not a nuisance to be averaged over.** | Median transcripts per cell span **6.5×** (42–273), QC-pass 79.6–99.3 %, Leiden clusters 13–20, realised fine labels 10–15, `Unknown` 0.2–19.8 %, `Low_quality` 0.0–12.5 %. Depth is in the frozen N5 block, but any cross-donor H1 statement must say which sections carry it. §2.1, §9.1 |

---

## 12. Environment, resources and provenance

**Interpreter:** the pinned Python 3.11 stack (`requirements.txt`); `/workspace/envs/sasp311`
present and unchanged. **No package was installed into it.** The only environment built was the
isolated DCA venv (`bash code/setup_dca_env.sh` with `DCA_ENV_ROOT=/tmp/dca_env`, outside the
repo) — DCA 0.3.4 / TensorFlow 2.4.4 / Keras 2.4.3 under CPython 3.8.19, exactly the recipe
PREREG §12 records. `py-spy` was pip-installed once for a stack dump and is not used by any
producer.

**Memory.** The ceiling is a **57.7 GiB cgroup** (`memory.max` = 61,999,996,928 B), not what
`free` reports. Observed peaks, from `memory.stat`'s `anon` (the binding quantity; `memory.current`
includes reclaimable page cache):

| stage | peak anonymous | note |
|---|---|---|
| `h1_annotate.py`, 196k-cell section | ~23 GB | scale + arpack PCA |
| `h1_annotate.py`, **393k-cell section (SPLN24)** | **~41 GB** | fits; the largest single job |
| `h1_deepscence.py`, one 393k-cell section | ~19 GB | full section, **no subsampling needed** |
| two concurrent `h1_deepscence.py` | ~28 GB | the sustainable configuration |
| **three concurrent `h1_deepscence.py`** | **> 57 GB → OOM** | see H8 |

**DeepScence fitted every H1 section at full size.** The pre-registration's worry that H1 sections
(220k–396k cells) might not fit did not materialise: the largest, SPLN24 at 393,202 cells ×
5,093 genes, ran in 41 min at ~19 GB anonymous. **No section was subsampled**, and the
`--subsample` path in `code/h1_deepscence.py` was never taken for a production score.

**Wall time**, one pass: annotation 10–17 min/section (SPLN24 16 min); callers + A6 + cache
4–6 min/section; module fits 10 min for 14 jobs at 4-way; A7 2.5 min for 14 jobs at 4-way;
DeepScence 16.5–41 min/section; DCA `denoise=True` 2.4–2.9 min to denoise 20,000 cells on CPU.

**Figures guard:** `python3 code/check_figures_guard.py` → `OK: all 52 committed figures match`.

**Files written by Phase 9** — all new paths:
- `code/h1_common.py`, `h1_sec.py`, `h1_a1_geometry.py`, `h1_a1_assignment.py`, `h1_annotate.py`,
  `h1_callers.py`, `h1_a6_compartments.py`, `h1_prep_cache.py`, `h1_a3_prevalence.py`,
  `h1_a4_ripley.py`, `h1_a5_matchbalance.py`, `h1_a7_controls.py`, `h1_a7_summarize.py`,
  `h1_a8_crossarm.py`, `h1_module_fits.py`, `h1_jobB_crosscheck.py`, `h1_caller_agreement.py`,
  `h1_senepy_surrogates.py`, `h1_deepscence.py`, `h1_deepscence_anchor.py`,
  `h1_deepscence_dca.py`, `h1_d2_analyse.py`, `h1_run_annotate_queue.sh`, `h1_run_stage2.sh`,
  `h1_run_deepscence_supervisor.sh`, `h1_run_dca_panel.sh`
- `results/phase9_h1/` — 30+ CSVs, all named in this report
- `data/processed_h1/` — per-section `celltypes_h1_*`, `senders_h1_*`, `modules_h1_*`,
  `anatomy_h1_*`, `deepscence_h1_*`, `cache3_h1/*.npz`
- `logs/phase9/`

**Not written:** `results/phase3/`, `figures/`, `genesets/`, `data/processed/`,
`results/phase7_jobA/` (the A2 re-run writes there and was verified to leave it byte-identical).

**Reuse rather than reimplementation.** A5, A7 and the module reference fits call
`phase3_core.match_decoys_section`, `run_phase3_nulls.SectionFit`, `fit_cell` and `_section_job`
**verbatim**; `code/h1_sec.py` only rebinds `sasp_phase3.CACHE3` and the section lists so those
functions read the H1 cache. `code/h1_annotate.py` **imports** every threshold from
`code/annotate_pipeline.py` rather than restating it. `code/h1_a4_ripley.py` is `code/_ripley.py`
with the cache substituted. `code/h1_caller_agreement.py` imports `_within_type_flags`,
`_matched_flags` and `_stratified_null` from `code/caller_disagree_all.py`.

---

## 13. What Phase 10 inherits

1. **The gates are clear.** A2 and A5 pass, so §18 outcome C is not triggered. A3 passes at the
   primary call on five compartments in all seven sections. A7's pre-registered primary null is
   flat.
2. **Two PI decisions are needed before the Phase-10 fits, and neither is mine to take.**
   (a) **H5** — which label family the Tier A sender call is defined on. Both are computed and
   both are in `a3_prevalence_by_type.csv`. (b) Whether SenePy stays in the primary trio, given
   that it now carries **four** caveats on this arm: no spleen hub (P4), the 100 µm window
   truncating its receivers (P8), a within-type depth enrichment of **28–225×** (§10.2), and the
   v1/v2 hub-release discrepancy (H4).
3. **`results/phase9_h1/h1_module_fits.csv` already exists** — 686 fits over 7 sections × 2 sender
   calls × the seven Tier B modules, produced by the frozen `_section_job` as the A7 reference.
   **It is NOT the Phase-10 primary outcome and R1 is deliberately not evaluated here**, because
   R1 requires the full N7 sender axis, the N1/N3/N4 perturbation nulls and the paired-bootstrap
   interval on the median SF, none of which Phase 9 ran. For orientation only, and marked as such:
   227 of 686 fits are reportable (33.1 %, against M1's 160/315 = 50.8 %), λ̂ rails at a grid bound
   in **73.9 %** of fits (M1 63 %), and the median SF under N2+N5+N6 over reportable fits is 0.373
   pooled over both calls. **These are provisional orientation numbers, not the replication test.**
4. **The A6 covariate is the weakest input Phase 10 will use** (§7) and its weakness must appear
   wherever an H1 anatomy-conditioned quantity does.
5. **Composition-matched rerun (roadmap 10.2).** `code/run_phase8_compmatch.py --arm h1` is gated
   on `ARMS['h1']['sections']` being populated and `SASP_H1_UNFROZEN=1`. The H1 cache now exists in
   the format it needs (`data/processed_h1/cache3_h1/`), reached through `code/h1_sec.py`.
6. **P-vi and P-vii are answered on one section; the remaining six are Phase 10** if the PI wants
   them. The environment recipe and both producers are in place
   (`code/setup_dca_env.sh`, `code/h1_deepscence_dca.py`, `code/h1_d2_analyse.py`).
   **P-vi's own falsifier needs 5 of 7 sections and Phase 9 ran 1**, so the prediction is
   contradicted-but-not-falsified and must be described that way until the other six run.
7. **A third PI decision, added by §10.5: whether DeepScence stays in the H1 caller set at all.**
   It is not seed-reproducible on this arm at full section size (H10). The options are to keep it
   with the seed caveat attached to every number, to report it only as a multi-seed ensemble, or
   to demote it. Phase 9 takes none of them.
8. **Three DeepScence run-metadata files are truncated (H7).** If DeepScence's internal node
   choice and `reverse` flag are wanted for SPLN21 / SPLN24 / SPLN30, those three sections must be
   re-run (~30 min each). The scores themselves are complete.
