# Pre-registration — Gene Sets (Tiers A–E), both arms

**Component file for `PREREG_PHASE8.md`, which the PI assembles and commits. This file is the
gene-set section only.**

**Biology collaborator · frozen 2026-08-27**
**M1** = GEO GSE310392, mouse liver, Xenium Prime Mouse 5K + 100 custom, **5,097 genes**
(5,106 `Gene Expression` features − 9 genotyping probes).
**H1** = GEO GSE326743, 7 normal human spleens (ages 17/31/32/32/37/57/59, 4M/3F), Xenium Prime
5K Human + 100 addon, **5,093 genes**.

Every count below is recomputed from files on disk by the scripts named in §7 and is reproduced in
`results/phase7_jobA/*.log`. No gene list in this document was written from memory.

**Freeze status of the data.** No H1 expression value, cell record or annotation has been read.
Panel membership only (`genesets/h1_candidate/`), which is the sanctioned §12.1 screen.

---

## 1. What is frozen, and where

> **Rule, enforced on disk:** every `.txt` directly in `genesets/human/` is FROZEN and is what the
> pipeline reads. Everything in `genesets/human/variants/` is **reported but not used**. The same
> rule applies to `genesets/mouse_c6/`. Machine-readable index with SHA-256 per file:
> `genesets/human/FROZEN_MANIFEST.csv` (35 frozen files, 8 variants).

| Arm | Frozen directory | Status |
|---|---|---|
| H1 human | `genesets/human/` | **Frozen. This is the configuration to run.** |
| M1 mouse | `genesets/` (promoted from `genesets/mouse_c6/`) | **Frozen and PROMOTED**, PI decision D5, 2026-08-27 05:41. All 15 mouse Tier A/B files in `genesets/` are byte-identical to `genesets/mouse_c6/`; three changed: `A_SENDER_FINAL_strict` 25→33, `A_sender_for_secondary_senescence` 55→74, `B_secondary_senescence` 38→108. Promotion means the mouse arm is **being re-fitted** (roadmap item 8.7). |
| M1 mouse, pre-C6 | tag **`pre-c6-genesets`** (`e789372`) | The provenance of all published Phase 2–5 mouse numbers. `genesets/*.txt` is no longer the pre-C6 state — read it from the tag (`git show pre-c6-genesets:genesets/…`), which is what `code/crossarm_geneset_table.py` and `code/corescence_circularity.py` do. |

---

## 2. Tier A — sender

### 2.1 Declaration

**PRIMARY: `A_SENDER_FINAL_strict`, n = 33 on both arms.** One sender definition, disjoint from
all seven response modules.

**PRE-REGISTERED SENSITIVITY: the seven `A_sender_for_<module>` sets**, each disjoint from the one
module it is paired with. These are not an informal alternative; they are gated in the same run
and both gates must pass.

**Rationale, recorded before any H1 expression was read.** H1's design is a continuous 17–59 age
axis over 7 donors. One sender definition yields one prevalence per donor and one regression
against age; seven per-module sender definitions yield seven prevalences per donor and no coherent
age model. The C6 decision (§4) raised the strict set from 21 to 33 genes, clearing the ≥15 floor
by more than 2×, which blunts the statistical-power argument that previously favoured per-module.
Both are reported.

**What the primary costs, stated against interest.** The 33 surviving genes still contain no
`CDKN1A`, `CDKN2A`, `TP53`, `LMNB1` or `MKI67`. `B4_downstream_arrest` removes them, and C6 does
not touch B4. The strict sender score is a DNA-damage / p53-effector / replicative-senescence score,
not a score built on the field's canonical arrest markers. The per-module sensitivity retains all
12 canonical markers for B2, B5, B6 and B7, which is exactly why it is pre-registered alongside.

### 2.2 Frozen membership

**H1 human — `genesets/human/A_SENDER_FINAL_strict.txt` (33):**

```
ATM ATR BAX BCL2 BCL2L1 CCNB1 CDKN1C CDKN2B DDB2 EHMT2 ERCC1 FOXM1 FOXO3 GADD45G GLB1 MDC1 MDM2
MDM4 NFATC1 PARP1 PHLDA3 RAD51 RB1 RBL2 SIRT1 TERF2 TERT TNFRSF10B TP53BP1 TP63 TP73 XPC XRCC5
```

**M1 mouse — `genesets/mouse_c6/A_SENDER_FINAL_strict.txt` (33):**

```
Atm Atr Bax Bcl2 Bcl2l1 Ccnb1 Cdkn1c Cdkn2b Ehmt2 Ercc1 Foxm1 Foxo3 Glb1 H2afx Hmga1 Lmnb2 Mdm2
Mdm4 Nfatc1 Parp1 Plk3 Rad51 Rb1 Rbl2 Rif1 Sesn1 Sirt1 Terf2 Tert Trp53bp1 Trp63 Trp73 Xrcc5
```

**Identical size, not identical membership** — see the asymmetry table in §6.

### 2.3 Per-module sensitivity sizes

| Module | M1 mouse | H1 human |
|---|---|---|
| B1 tnfa_nfkb_proximal | 70 | 77 |
| B2 il6_jak_stat3 | 74 | 81 |
| B3 interferon_response | 73 | 80 |
| B4 downstream_arrest | 37 | 36 |
| B5 emt_ecm | 73 | 79 |
| B6 oxidative_stress | 71 | 79 |
| B7 secondary_senescence | 74 | 81 |

All fourteen clear ≥15 and every one is disjoint from its own readout. B7's per-module set now
equals the full candidate pool A0 on both arms, because C6 makes B7 disjoint from A0 by
construction.

### 2.4 Component subsets (human, post-disjointness)

`A_core_arrest` 6, `A_proliferation_down` 2, `A_nuclear_chromatin` 1, `A_dna_damage_response` 11,
`A_senescence_curated_nonsecreted` 13. Reference set: `A6_SenMayo_arrest_reference` 33.

---

## 3. Tier B — response, seven modules

| Module | Source sets | M1 mouse | H1 human |
|---|---|---|---|
| B1 tnfa_nfkb_proximal | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 126 | 120 |
| B2 il6_jak_stat3 | `HALLMARK_IL6_JAK_STAT3_SIGNALING` | 68 | 71 |
| B3 interferon_response | `HALLMARK_INTERFERON_GAMMA_RESPONSE` ∪ `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 100 | 126 |
| B4 downstream_arrest | `HALLMARK_E2F_TARGETS` ∪ `HALLMARK_G2M_CHECKPOINT` ∪ `HALLMARK_MYC_TARGETS_V1` | 190 | 231 |
| B5 emt_ecm | `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` | 125 | 113 |
| B6 oxidative_stress | `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` ∪ `REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES` ∪ curated NRF2 targets | 31 | 36 |
| **B7 secondary_senescence** | **`SAUL_SEN_MAYO` ∪ `REACTOME_SASP`, minus the Tier A caller pool A0** | **108** | **116** |

§10 names seven HALLMARK sets; `E2F_TARGETS` and `G2M_CHECKPOINT` both feed B4, so seven named
sets plus `secondary_senescence` make **seven** modules, not eight.

**B6 needed the same documented rescue in both arms.** `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY`
alone gives 17 on-panel in mouse and 18 in human, both below the ≥30 floor.

### 3.1 ⚠ B6 has essentially no margin — anything that trims it fails the gate

Margin over the ≥30 floor (`results/phase7_jobA/module_margins_human.csv`, and the mouse audit in
`build_genesets_mouse_c6.log`):

| Module | M1 mouse | H1 human |
|---|---|---|
| B1 tnfa_nfkb_proximal | +96 | +90 |
| B2 il6_jak_stat3 | +38 | +41 |
| B3 interferon_response | +70 | +96 |
| B4 downstream_arrest | +160 | +201 |
| B5 emt_ecm | +95 | +83 |
| **B6 oxidative_stress** | **+1** (n=31) | **+6** (n=36) |
| B7 secondary_senescence | +78 | +86 |

**`B6_oxidative_stress` is the only module anywhere near the floor, and on the mouse arm it clears
it by a single gene.** Consequences that must be pre-registered:

- **Any later change that removes even one gene from mouse B6 fails the §11 gate.** That includes
  a future over-adjustment guard, a QC filter that drops low-detection genes, an MSigDB release
  bump, or dropping the curated NRF2 component of the documented B6 rescue.
- The gate script must therefore be **re-run after any change to a gene set or a panel**, not only
  at freeze time. It exits non-zero on failure precisely so this cannot pass silently.
- B6 is already the module carrying a documented substitution (§3). Its result should be read as
  the weakest of the seven and reported as such, independently of whether it passes.

### 3.2 Which panel definition is authoritative

The mouse arm has **two** panel files, and they answer different questions. Recomputing B6 against
each gives **31** or **30** — both pass, but the difference matters and is documented here so it
cannot be rediscovered as an inconsistency.

| Definition | n | What it is |
|---|---|---|
| `XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv` | 5,006 | The **stock** Prime Mouse 5K pan-tissue panel only. Does **not** include GSE310392's custom add-on. |
| + `GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz` | 5,106 | The union. `genesets/README.md` verified this equals the `cell_feature_matrix.h5` `Gene Expression` feature list **exactly, 0 discrepancies either way**. The two files are disjoint (overlap = 0). |
| **minus the 9 genotyping probes** | **5,097** | **AUTHORITATIVE.** Genotyping probes (`Jak2_WT_1940`, `Kras_ALT_235:GA`, …) are not genes. |

**The 5,097-gene set is authoritative for the gate**, on three grounds: it is what the instrument
measured; it is what `annotate_pipeline.py` and `build_random_null_sets.py` already select
(`feature_type == 'Gene Expression'` **and** id starts with `ENSMUSG`); and the 5,006-row CSV is
not a panel definition at all — it is one of the two files whose union is the panel. The counts are
asserted in `code/build_genesets_mouse_c6.py`, which stops if they drift.

The custom add-on contributes to **every** module, not just B6:

| Module | authoritative (5,097) | CSV alone (5,006) | add-on genes |
|---|---|---|---|
| B1 | 126 | 119 | Areg Fos Il6 Jun Junb Klf6 Nfkbia |
| B2 | 68 | 65 | Il6 Jun Reg1 |
| B3 | 100 | 97 | Cfh Il6 Nfkbia |
| B4 | 190 | 189 | Hmgb2 |
| B5 | 125 | 115 | Acta2 Areg Il6 Jun Plod2 Pvr Spp1 Tgfbi Thy1 Vim |
| **B6** | **31** | **30** | **Junb** |
| B7 | 108 | 104 | Areg Il6 Jun Spp1 |

So the coordinator's 30 is correct *for the stock-panel-only definition*, and the one gene at issue
is **`Junb`**. Under the CSV-only definition mouse B6 sits **exactly at the floor with margin 0**.
The human arm has no such ambiguity: GSE326743 ships a single panel and the 5,093 genes were read
straight from the `.h5` feature table.

---

## 4. C6 — the adopted B7, and the method that was replaced

**§23/C6 prescribes: "rebuild B7 without the shared genes." That method fails.** Applied to the
curated module it leaves **24** genes against the §10 caller and **12** against the ported caller,
both below the §11 floor of 30 (`results/phase7_jobA/b7_c6_rebuild.csv`).

**Adopted instead:** re-source the module from MSigDB sets carrying peer-reviewed identifiers,
*then* subtract the caller pool. On the human arm 123 on-panel → **116**; on the mouse arm 113 →
**108**. Both pass, both are disjoint from the caller by construction, and the construction is
closer to §9's own definition of B7 ("Tier A minus the calling genes") than the curated list was.

The subtracted caller is **A0, the Tier A candidate pool before disjointness enforcement**. That
makes the fixed point self-consistent: `A_SENDER_FINAL_strict` ⊂ A0, so it is disjoint from B7
whichever Tier A option is fitted.

**Citation.** The superseded module rested on `neretti2024dissecting`, a ~250-word GSA meeting
abstract that audit finding B8 marks "DO NOT CITE AS A PAPER". Replacements, all already in
`references.bib`, none invented: `saul2022senmayo` (**PMID 35974106**, taken from the archived
MSigDB record, not from recall) for the module content, Reactome pathway R-HSA-2559582 for the
SASP half, and `martin2023modelling` + Acosta et al. 2013 for the primary-versus-secondary
distinction itself — the audit's own recommended replacements.

**Both versions are on disk, as §23/C6 requires:**
`genesets/human/variants/B_secondary_senescence_v1_curated_ported.txt` (35) and
`genesets/mouse_c6/variants/B_secondary_senescence_v1_curated.txt` (38).

**What C6 costs, stated up front.** Re-sourcing B7 makes it *more* entangled with the other
response modules: B7∩B1 rises from 10 to 26 genes, B7∩B5 from 6 to 20. And it raises the
DeepScence circularity figure — see §5. The trade was accepted because circularity between the
**sender** and a readout is a validity threat, while correlation **among readouts** is a
multiple-comparison and interpretation issue that the intersection matrix already discloses.

---

## 5. CoreScence circularity — the number to cite is 88%

DeepScence is the sender caller for Job B. Its CoreScence gene set (occurrence ≥5, 39 genes, 33
on the human panel and 33 reachable on the mouse panel) overlaps the response modules:

| Configuration | in ≥1 Tier B module | fraction |
|---|---|---|
| M1 mouse, pre-C6 (published reference) | 26 / 33 | 79% |
| **M1 mouse, FROZEN re-sourced B7** | **29 / 33** | **88%** |
| H1 human, **superseded** curated B7 | 25 / 33 | 76% |
| **H1 human, FROZEN re-sourced B7** | **29 / 33** | **88%** |

> **Correction, 2026-08-27.** The mouse row read `24 / 35 = 69%` in every earlier draft of this
> document, in `code/make_figure_genesets.py` and in `code/gate_disjointness_human.py`. It was a
> typed-in literal produced by no script, in violation of §7b, and the denominator 35 is
> reproducible under **no** mapping convention. The number is now derived from files by
> `code/corescence_circularity.py` → `results/phase7_jobA/corescence_circularity_mouse.json`.
> Of the 39 CoreScence genes, **31** are reachable on the 5,097-gene mouse panel through the
> pinned 1:1 MGI map and **33** with the project's own documented Title-case fallback for
> symbols the map has no row for (`CDKN2B`, `CXCL1` — the same map-gap class as §6 below);
> 6 (`CXCL8 HMGB1 IGFBP5 IGFBP7 MIF TNFRSF10C`) are off the mouse panel entirely. The
> numerator 24 was right for the strict convention (24/31 = 77%); the pair `24/35` mixed the
> two. **33 is the denominator the mouse arm actually ran under**: it is what
> `code/run_phase3_n8.py::corescence_mouse` implements and the value of `corescence_on_panel`
> in the committed `results/phase3/n8_disjointness_*.csv`, whose per-module overlaps
> (10/9/5/5/0/14/8, union 26) reproduce exactly. Source: `reports/AUDIT_PHASE8_FACTCHECK.md` M1.

**Cite 88%.** It is the number for the configuration that will actually be run. Report 79% and 76%
alongside so the movement is visible: the rise is a direct, disclosed consequence of the C6
decision, because SenMayo and CoreScence are both senescence signatures and B7 grew from 35 to 116
genes (human) and 38 to 108 (mouse). B7 alone now accounts for 18 of the 33 (55%) on **both** arms.

**The cost of C6 is smaller than previously stated, and it is now measured within an arm.**
The old "69% → 76% → 88%" narrative compared a mouse baseline with two human configurations and
put the movement at 19 points. Measured properly it is **+9 points on the mouse arm (79 → 88)**
and **+12 on the human arm (76 → 88)**, and the two arms land on the same 29/33. C6 still raises
circularity, the strip-and-refit sensitivity is still required, and the pre-C6 mouse arm was
already substantially more circular than the paper has so far reported (79%, not 69%) — which
makes the *published* mouse result more circular, not less, and is a correction against interest.

Per module under the frozen Tier B, H1 human: B7 0.55, B4 0.30, B5 0.27, B1 0.21, B3 0.18,
B2 0.12, B6 0.00. M1 mouse: B7 0.55, B4 0.30, B5 0.27, B1 0.24, B2 0.15, B3 0.15, B6 0.00.

**Required consequence:** the strip-and-refit sensitivity (fit each module with the CoreScence-shared
genes removed, report both amplitudes) is **part of the frozen run order**, not an afterthought.
CoreScence is a human gene set, so on H1 it runs natively with no ortholog remapping — the 88%
cannot be attributed to mapping loss.

---

## 6. Cross-arm comparability (Phase 7 test A8)

Full table: `results/phase7_jobA/crossarm_geneset_table.csv`.

**Ortholog-intersected panel: 2,425 genes.** Of 5,097 mouse panel genes, 4,845 have a row in the
pinned MGI map; **2,435** of those map onto the 5,093-gene human panel, and they land on **2,425
distinct human symbols** — that human-symbol count is the intersected panel. 252 mouse panel
genes have no ortholog in the map; 2,668 human panel genes are unreachable from the mouse panel.
**The map is not 1:1**: its 18,782 rows (one per mouse gene) carry 17,609 distinct human symbols,
and 431 human symbols receive more than one mouse gene. Earlier drafts said "4,845 have a 1:1
human ortholog … and 2,425 of those land on the human panel", which conflates the mouse-side
count with the human-side one; both are correct numbers of different things
(`reports/AUDIT_PHASE8_FACTCHECK.md` M2). **Every cross-arm number must be reported both on this 2,425-gene intersection and on each
arm's full panel**, as A8 requires. On-panel-and-on-intersection counts are the last two columns
of the CSV.

### Where the arms cannot be made identical — declared

| Asymmetry | Quantified | Status |
|---|---|---|
| `SAUL_SEN_MAYO` is not a translation | mouse MM16098 n=117, human M45803 n=124; 111 of 117 mouse members map; 110 overlap; **14 human-only, 1 mouse-only** | Declared. Both arms use MSigDB's own species version, as the mouse arm did for B1–B6. |
| `REACTOME_SASP` is far larger in human | mouse MM14900 n=**40**, human M27187 n=**111**; all 38 mappable mouse members are inside the human set; **73 human-only** | Declared, and the main reason human B7 (116) exceeds mouse B7 (108). |
| Frozen B7 membership | mouse 108, human 116; 104 mouse members map, **85 overlap by the pinned map, 88 after the map-gap correction** (gaps: `Ccl3`, `Cxcl2`, `Cxcl3`). Mouse-only, complete (20 = 108 − 88): the 19 mapped ones plus `Cxcl1`. Human-only, real (28 = 116 − 88). On the intersected panel, mouse 85 / human 88 — a **different** quantity that coincides numerically with the corrected overlap. | Declared. |
| Frozen Tier A membership | both 33; 31 mouse members map; **26 overlap by the pinned map, 27 after the map-gap correction** (gap: `Cdkn2b`/`CDKN2B`, which is a member of **both** frozen Tier A sets). Mouse-only, complete (6 = 33 − 27): `HMGA1 LMNB2 PLK3 RIF1 SESN1 H2afx`. Human-only, real (6 = 33 − 27): `DDB2 GADD45G MDC1 PHLDA3 TNFRSF10B XPC` | Declared. Equal size is a coincidence of the two panels, not identity. **Corrected 2026-08-27:** earlier drafts listed `CDKN2B` as human-only and gave a `mouse_only` list of 5 where 33 − 26 = 7 members are unaccounted for. `CDKN2B` is in both arms and is a map gap, not an asymmetry — the same failure mode this table flags two rows down for `CXCL2`/`CXCL5` (`reports/AUDIT_PHASE8_FACTCHECK.md` R2). `H2afx` **is** genuinely mouse-only: human `H2AFX` sits in B4 `downstream_arrest` and is removed from human Tier A by disjointness. |
| `CXCL8` / `CXCR1` | **no mouse ortholog exists** — `CXCL8` is absent from the MGI map entirely | Declared. Human-only quantity; must never be reported as replicating a mouse result. M1 used `Cxcl1`/`Cxcl2`/`Cxcl5`→`Cxcr2` as analogues. |
| `IL1B` | on the human panel, **off** the mouse panel | Declared. The sharp §9 internal control — membrane-bound `IL1A` vs secreted `IL1B` through the same `IL1R1` — is runnable on H1 and was **not** runnable on M1. |
| `MMP3`, `TIMP1` | on the mouse panel, **off** the human panel | Declared. Mouse-only quantities. |
| `CXCL2`, `CXCL5` (Tier C); `CDKN2B` (Tier A); `Ccl3`, `Cxcl2`, `Cxcl3` (B7) | on **both** panels and in both arms' sets; the pinned MGI map simply lacks the rows | **Map gap, not a biological asymmetry.** Do not report these as arm-exclusive. `code/crossarm_geneset_table.py::gap_split` now applies this split to **every** set and asserts that both arms' non-shared members are accounted for by name. |
| Housekeeping control | mouse 13/39 on-panel, human **8/39** | Declared. Human control is thinner; report as such. |
| Tissue | mouse liver, human spleen | Declared already in §15: no cross-arm difference can be attributed to species or tissue. Now more emphatic — the plan assumed lung. |

---

## 7. Pins, and how to rebuild

| Pin | Detail |
|---|---|
| MSigDB human | release **2026.1.Hs**, fetched **2026-08-27**. `genesets/msigdb_human_2026.1.Hs/` — 26 per-set JSON + the complete H collection bundle (50 sets). All 27 files valid JSON. |
| MSigDB mouse | release **2026.1.Mm**, fetched **2026-08-20**. `genesets/msigdb_mouse_2026.1.Mm/` — pre-existing pin, **not re-downloaded**. **3 of its 26 files are HTML error pages, not JSON** (`FRIDMAN_SENESCENCE_UP/DN`, `WP_NRF2_PATHWAY`): MSigDB has no mouse version of those sets. None feeds any tier, so no published mouse number is affected. `code/build_genesets.py` hid this behind a bare `except: pass`; the new scripts report it. |
| HGNC | downloaded **2026-08-27**, md5 `2d741e796d5538cc48d3696452237781`. `genesets/hgnc_pin/`. Symbol/alias/previous-symbol → Ensembl → panel symbol. Dated monthly archive URLs returned 404, so the pin is by download date + md5. |
| CellMarker 2.0 human | downloaded **2026-08-27**, md5 `c7a1b764b66cb3a3c16cfac428160f72`, 60,877 rows. `genesets/cellmarker_pin/`. Sole source of the spleen cell-type markers. |
| MGI orthologs | `genesets/mouse_human_orthologs_MGI.csv`, pre-existing pin already used by `run_deepscence.py`. **18,782 rows, one per mouse gene, onto 17,609 distinct human symbols — it is many-to-one, not 1:1**, and it has documented gaps (§6). |
| DeepScence CoreScence | `DeepScence/data/coreGS_v2.csv` from the installed package, occurrence ≥5. |
| Mouse panel | `XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv` ∪ `GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz`, minus 9 genotyping probes = **5,097 genes**. See §3.2 for why this, not the 5,006-row CSV, is authoritative. |
| Human panel | `genesets/h1_candidate/GSE326743_gene_panel_5093.csv`, read from the `.h5` feature table = **5,093 genes**. Single panel, no ambiguity. |

**Rebuild order** (all offline except the two 2026-08-27 downloads, which are archived):

```
python3 code/build_genesets_human.py
python3 code/rebuild_b7_secondary_senescence.py
python3 code/build_markers_human_spleen.py
python3 code/spec_a6_compartments_human.py
python3 code/senepy_coverage_human.py
python3 code/build_genesets_mouse_c6.py
python3 code/crossarm_geneset_table.py
python3 code/corescence_circularity.py         # mouse CoreScence anchor, derived not typed
python3 code/gate_disjointness_human.py        # exit 0 = frozen configuration passes
python3 code/freeze_manifest_human.py
python3 code/make_figure_genesets.py           # Methods figures, from the CSVs above
```

`code/build_genesets_mouse_c6.py` **asserts** that it reproduces the mouse arm's published sizes
(B1 126, B2 68, B3 100, B4 190, B5 125, B6 31, B7 v1 38, A0 74) before applying C6, and stops if
it does not.

---

## 7b. Methods figures

Built by `code/make_figure_genesets.py` from the CSVs above — nothing is recomputed or typed in.
(This was violated once, by gs3's mouse bar; see the correction box in §5. The producing script is
`code/corescence_circularity.py` and the claim now holds.)
Each writes `.png`, `.pdf` and a `_data.csv` beside it in `figures/`. Styling uses
`code/sasp_palette.py` (`apply_style`, sequential `SEQ`), matching every other `make_figure*.py`.

| Figure | Shows | Source CSVs |
|---|---|---|
| `figures/figure_gs1_intersection_matrix` | **§11 disjointness, both arms.** (a) candidate Tier A definitions × Tier B before removal, including the §10 list that failed; (b) the frozen sets × Tier B with the **gate cells outlined in green** — the PRIMARY row against all seven modules, each per-module sensitivity row against its own module — every outlined cell 0; (c) Tier B × Tier B cross-talk. Exact zeros are painted flat neutral so 0 is unambiguously distinguishable from small-but-nonzero. | `intersection_matrix_human.csv`, `intersection_matrix_BxB.csv`, `intersection_matrix_mouse.csv`, `intersection_matrix_mouse_BxB.csv`, and the two `_preC6` variants |
| `figures/figure_gs2_crossarm_symmetry` | **§17 cross-arm symmetry**, every tier, pre- and post-C6, plus the ortholog-intersected counts against the 2,425-gene panel. The dark inner segment on the frozen bars is how many genes the *other* arm also has, map-gap corrected: **27 of 33** for Tier A (26 by the pinned map), **88 of 108/116** for B7 (85 by the pinned map). The `_data.csv` now carries those plotted values and the whole footnote block, which previously lived only in the JSON. | `crossarm_geneset_table.csv/.json` |
| `figures/figure_gs3_corescence_circularity` | **mouse 79% → 88%, human 76% → 88%**, with the per-module breakdown showing B7 alone contributing **18/33** on both arms. | `corescence_circularity_mouse.json`, `gate_result_human.json` |
| `figures/figure_gs4_senepy_coverage` | **SenePy spleen coverage**: 22 labels vs the best available hub, and the hub collapse (one blood memory-B hub for all three B labels, one lung T-cell hub for both T subsets). | `senepy_spleen_coverage.csv/.json` |

`figure_gs1` **asserts in code** that every gate cell is 0 on both arms; the figure cannot be
produced from a configuration that fails the gate.

---

## 8. §11 gate result — frozen configuration

`results/phase7_jobA/gate_disjointness_human.log`, exit status 0.

| Assertion | H1 human | M1 mouse (C6) |
|---|---|---|
| `len(A_strict & panel) >= 15` | **33** PASS | **33** PASS |
| `len(B_k & panel) >= 30`, k = 1…7 | 120 / 71 / 126 / 231 / 113 / 36 / 116 — PASS ×7 | 126 / 68 / 100 / 190 / 125 / 31 / 108 — PASS ×7 |
| `len(A_strict & B_k) == 0`, k = 1…7 | **0** ×7 PASS | **0** ×7 PASS |
| Per-module: `len(A_mod) >= 15` and `A_mod ∩ B_k = ∅` | PASS ×7 | PASS ×7 |
| **Verdict** | **PASS** | **PASS** |

Full intersection matrix (Tier A × Tier B, counts and per-cell gene lists) and the Tier B × Tier B
cross-talk matrix: `results/phase7_jobA/intersection_matrix_human.csv` and
`intersection_matrix_BxB.csv`.

**Tier B modules are not disjoint from each other and are not required to be**: B7∩B1 = 26,
B1∩B3 = 23, B7∩B5 = 20, B2∩B3 = 22. Seven modules are not seven independent tests. This must be
one sentence in Methods.

---

## 9. Tier C — ligand–receptor

`C_ligands` 13 human / 14 mouse; `C_receptors` 15 / 15. Human-only: **`CXCL8`** with `CXCR1` and
`CXCR2` (§10 addition; no mouse ortholog for either), and **`IL1B`**. Mouse-only: `CXCL1`, `MMP3`,
`TIMP1`, `SDC1`. See §6.

---

## 10. Tier D — nuisance

`D_nuisance_covariates` (13 covariate **names**, not genes): cell type, density at 25/50/100 µm,
k-NN composition, total counts, genes detected, cell and nucleus area, section, segmentation
method, distance to tissue boundary, plus the arm-specific anatomical covariate.

**Arm-specific anatomical covariate.** M1: liver zonation (already derived, `genesets/D_zonation_*`).
H1: the **red pulp / white pulp axis**. §13's A6 as written specifies a lung airway-to-alveolar
axis and is **void for this arm**. The gene side is frozen —
`D_spleen_white_pulp_tzone` (24), `D_spleen_white_pulp_follicle` (13), `D_spleen_marginal_zone` (8),
`D_spleen_red_pulp` (28), `D_spleen_capsule_trabecula` (18), pairwise disjoint, Tier A/C overlap 0
by construction. The covariate itself — `compartment`, `dist_to_white_pulp_um`,
`dist_to_marginal_zone_um`, `dist_to_capsule_um` — is **specified but not built**, because the
mouse equivalent was derived from the expression matrix. Full specification:
`reports/BIO_PHASE7_JobA_FOLLOWON.md` §2. Also frozen: `D_landmark_vessel_stroma` (13).

**Cell-type marker set for H1** (`code/markers_human_spleen.py`, 22 labels, every gene traced to
PMIDs in the CellMarker pin). Three labels the panel cannot resolve are declared, not hidden:
**follicular dendritic cells** (3 markers), **plasma cells** (3), **proliferating cells** (3 after
the over-adjustment guard). See §12.

---

## 11. Tier E — controls

| Set | H1 | Status |
|---|---|---|
| `E_negative_control_probes` | **40** | Primary technical null. Audit test **A7** fits these against distance-to-sender and requires flat. |
| `E_negative_control_codewords` | 609 | A7. |
| `E_genomic_controls` | 21 | A7. |
| `E_housekeeping_expanded` | 8 (mouse 13) | Per-cell flat-kernel control. **Thin — report as thin.** |
| `E_housekeeping_planA` | 1 (`TBP`) | §10's six; five are off-panel. **Unusable as written, in both arms.** |
| **E2 cell-type identity** | **DROPPED** | See §12 and the deviation table. |
| E3 500 matched random sets | **not built** | Requires the H1 expression matrix to bin genes by mean expression. Post-freeze; `code/build_random_null_sets.py` re-points unchanged. |

---

## 12. Deviations from the Phase 7 doc — complete table

| # | Deviation | Reason |
|---|---|---|
| **D1** | **§10's Tier A (the literal 16 symbols) is NOT the sender definition.** Replaced by `A_SENDER_FINAL_strict`. | It fails A2/§11 twice: only **14** of 16 are on the human panel (`TP53I3`, `HMGB1` absent, verified by Ensembl ID), which is below the ≥15 floor before disjointness is tested; then 9 of the 14 collide with a module and **5 survive**. It fails identically on the mouse panel (13/16 on-panel, 5 survive). Archived at `genesets/human/variants/A_PHASE7_S10_16.txt` and reported as the pre-registered definition that failed. |
| **D2** | **`H2AX` is written as `H2AFX`.** | The panel ships the legacy symbol for `ENSG00000188486`. Resolved through HGNC and logged in `genesets/human/_symbol_resolutions.csv`. **Not dropped.** |
| **D3** | `TP53I3` and `HMGB1` absent from Tier A. | Genuinely not on the panel — checked by literal symbol, by HGNC approved-symbol → Ensembl (`ENSG00000115129`, `ENSG00000189403`) → panel `gene_id`, and by direct grep. No substitute invented. |
| **D4** | **§23/C6's stated method replaced.** "Rebuild B7 without the shared genes" → re-source from SenMayo ∪ Reactome SASP, then subtract. | C6 as written leaves 24 / 12 genes, below the §11 floor of 30. §4. Applied to **both** arms. |
| **D5** | **Tier E2 dropped for the human arm.** | §10 asks for "a cell-type identity program unrelated to inflammation". In spleen every abundant cell type is an immune cell, and **all five candidates are members of an A6 compartment by construction** — erythroid and sinusoidal identity *are* the red pulp definition. That is the exact failure that killed the mouse E2 (r = −0.67 with the zonation axis, worse than 99% of random size-matched sets). The 40 negative control probes are a strictly better technical null and A7 already requires them. `results/phase7_jobA/tierE2_candidates.csv`. |
| **D6** | **A6 becomes the red pulp / white pulp axis.** | §13's A6 specifies a lung airway-to-alveolar axis plus fibroblastic-focus distance. **H1 is spleen.** Every word of the lung version is void for this arm. §10 above. |
| **D7** | **SenePy is not the same estimator across arms.** | §13. |
| **D8** | **Tier C asymmetry, both directions.** `CXCL8`+`CXCR1` and `IL1B` human-only; `CXCL1`, `MMP3`, `TIMP1`, `SDC1` mouse-only. | §10 adds `CXCL8`; no mouse ortholog exists. The others are panel differences. `CXCL2`/`CXCL5` are on both panels and are a **map gap, not an asymmetry**. §6. |
| **D9** | B3 keeps `INTERFERON_ALPHA_RESPONSE` and B4 keeps `MYC_TARGETS_V1`, neither named in §10. | Mouse-identical source sets, for cross-arm comparability (A8). Dropping them would change human B3 to 113 and B4 to 178 and make the arms non-comparable. |
| **D10** | B6 rescued in both arms, **and it has no margin**. | `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` alone gives 17 (mouse) / 18 (human) on-panel, below ≥30. After the documented rescue: 31 (mouse, **margin +1**) and 36 (human, +6). Every other module clears the floor by ≥38. Anything that later trims B6 fails the gate — see §3.1. |
| **D17** | The mouse panel is defined as the **5,097-gene** set (both panel files, minus 9 genotyping probes), not the 5,006-row stock-panel CSV. | The two files are disjoint and their union was verified to equal the `.h5` feature list exactly. The CSV alone omits GSE310392's custom add-on, which contributes to all seven modules and puts B6 at exactly 30 — the floor, margin 0. §3.2. |
| **D11** | "No secreted factors" enforced by transferring the mouse panel's `location` annotation through the ortholog map. | The human Prime 5K panel-metadata CSV is not on disk. 263 human panel genes inherit `Secreted`; **0** Tier A members do. |
| **D12** | Cell-type markers carry a Tier A + Tier C over-adjustment guard; the mouse arm's did not. | The cell-type call is a Tier D nuisance covariate; building it from the outcome genes conditions the response on itself. `markers_mouse_liver.py` calls `'Proliferating'` with `Mki67 Top2a Pcna Ccnb1 Ccna2 Birc5 Aurkb Cdk1`, all Tier A/B. This is a correction. Guard not extended to Tier B: that removes 58 marker slots and destroys six labels including both stromal labels A6 depends on, so Tier B overlaps are flagged (`markers_spleen_evidence.csv`), not removed. |
| **D13** | E3 random matched sets not built; A6 covariate not built; `DROP_NONSPECIFIC` empty; §14's image check not run. | All require data behind the §15 freeze, or files deliberately not downloaded (`morphology.ome.tif.gz`). Each is named with what unblocks it. |
| **D14** | Tier E gains three control-feature files (40 / 609 / 21). | §12.4 and A7 need them; they have no mouse counterpart file. |
| **D15** | ~~Mouse C6 sets are in `genesets/mouse_c6/`, not promoted into `genesets/`.~~ **SUPERSEDED 2026-08-27 05:41 by PI decision D5: the C6 sets ARE promoted.** `genesets/*.txt` is now the C6 state (15/15 byte-identical to `genesets/mouse_c6/`); the pre-C6 state is the tag `pre-c6-genesets` (`e789372`). | All published Phase 2–5 mouse numbers were computed under the pre-C6 sets, so promotion required re-fitting the mouse arm — that re-fit is roadmap item 8.7 and is running. Anything that needs the pre-C6 sets must read them from the tag, not from `genesets/`. |
| **D16** | `kneed==0.8.6` and `openpyxl==3.1.5` added to `requirements.txt`. | `kneed` is a **hard** import of DeepScence (`DeepScence/io.py: from kneed import KneeLocator`) and was missing from the pins, so the environment could not rebuild itself. §14. |

---

## 13. SenePy — caller 2 is not the same estimator across arms

**Declared deviation, not a footnote.** SenePy 1.0.1's bundled human hub file carries **65 hubs
across 10 tissues**: blood, bone marrow, heart, hippocampus, intestine, kidney, liver, lung, skin,
tongue. **Spleen is absent. Liver is present.**

M1 scored SenePy with **tissue-matched** mouse Liver hubs (`phase2_downstream.py` `HUBMAP`:
`Liver/hepatocyte`, `Liver/endothelial cell of hepatic sinusoid`, `Liver/Kupffer cell`). H1 cannot.
Applying the mouse arm's own threshold (≥10 hub genes on panel), of the 22 spleen labels:

- **15 get a cross-tissue surrogate**, and they collapse onto far fewer hubs than labels: the
  entire B compartment (follicular, marginal zone, germinal centre) is scored by **one** blood
  memory-B hub; both T subsets by **one** lung T-cell hub; sinusoidal and general endothelium by
  **one** skin endothelial hub; FRC and fibroblasts by **one** skin fibroblast hub.
- **7 get no SenePy score at all**: **cDC1, cDC2, pDC, lymphatic endothelium, erythroid cells,
  megakaryocytes, mesothelial cells.** SenePy has no hub for these in any tissue.

**What it costs, concretely.** SenePy sender prevalence on H1 is computed over a **subset** of the
section — every cell assigned to one of those 7 labels is `NaN`, exactly as the mouse code already
handles. **Audit test A3's "1–20% prevalence, ≥200 senders, ≥5,000 non-senders" must therefore be
evaluated on that subset for the SenePy caller, and on all cells for the other callers**, and the
two denominators must be stated separately. Any cross-arm comparison of a SenePy-called quantity
compares a tissue-matched caller against a cross-tissue surrogate, and every H1 number derived from
`senepy_p90/p95/p99` (`run_phase3_*.py`, `summarize_super_callers.py`, `sasp_phase3.py`) carries
that caveat.

**Open recommendation for the PI:** demote SenePy to a sensitivity analysis on H1 and promote the
`CDKN1A`⁺ call — which needs no tissue-matched resource — into the primary trio. Flagged, not taken.

---

## 14. Environment

`kneed==0.8.6` and `openpyxl==3.1.5` were added to `requirements.txt` on 2026-08-27.

**The reproducibility claim was broken until now.** `kneed` is a hard import of DeepScence, so
`import DeepScence` raised `ModuleNotFoundError` on a clean rebuild from the pins — i.e. the pinned
environment could not reconstruct the sender caller that produces `deepscence_score`. It had to be
installed by hand after a container reset. `openpyxl` is needed only to *regenerate* the CellMarker
pin, not to re-run any human-arm build, but it was likewise unpinned. Both are now in the file,
each in a commented block explaining what it is for.

---

## 15. Two open items for the PI

**(a) Plasma cells — recommend lowering `MIN_MARKERS` to 3 for this label only.**
Its surviving on-panel markers are `JCHAIN`, `MZB1`, `XBP1`. It fails the label gate by exactly one
gene, and only after the over-adjustment guard removes `MKI67` (which was CellMarker contamination
from "dividing plasma cell" rows, not a plasma-cell marker). Those three genes are among the most
specific markers in immunology and are not shared with any other label in the set. Plasma cells are
also a genuine senescence-relevant population in spleen and a plausible receiver. Against: a
3-gene score is noisier and `MIN_MARKERS = 4` was set for a reason. **Recommendation: admit plasma
cells with a 3-marker exception, recorded as such, rather than folding them silently into the B
compartment — a plasma cell is not a B cell for receiver purposes.** This is the PI's to overturn.

**(b) Marginal zone B cells — recommend pre-registering the label but NOT a marginal-zone-specific
claim.** CellMarker 2.0 has only **6** spleen rows for this type, below the 8-row threshold, so its
markers (`CD180 CR2 EBF1 FCRL4 ITGAE NOTCH2 PAX5 TNFRSF13C`) come from the **weakest** evidence
tier — an all-tissue fallback at ≥1 PMID. It is the compartment that most defines spleen, and
`D_spleen_marginal_zone` (8 genes) is one of the five A6 compartments, which makes it structurally
load-bearing. Keeping it costs nothing if it is honestly scoped; dropping it would leave the A6
axis without its middle term. **Recommendation: keep the label and the compartment, and
pre-register that any claim specific to the marginal zone is exploratory and conditional on the
label surviving the post-freeze re-gate against measured expression and against the depositors'
own annotations.** Do not pre-register a confirmatory marginal-zone hypothesis.
