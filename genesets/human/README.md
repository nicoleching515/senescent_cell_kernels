# Gene Set Package — HUMAN arm (H1) · Phase 7 Part III / Job A

**Owner:** biology collaborator · **Date built:** 2026-08-27 · **Species: HUMAN**
**Target panel:** GSE326743 Xenium = *Xenium Prime 5K Human* + 100-gene addon = **5,093
`Gene Expression` features**, verified on `GSM9638040_cell_feature_matrix.h5` and cross-checked
identical against GSM9638044 and GSM9638046 (symmetric difference 0). Provenance and the
pre-freeze caveat: `../h1_candidate/PROVENANCE.md`.

This is the human port of `../README.md` (mouse). File names, formats and tier semantics are
identical so downstream code needs no change — point it at `genesets/human/` instead of
`genesets/`.

Rebuild: `python3 /workspace/code/build_genesets_human.py`
Gate:    `python3 /workspace/code/gate_disjointness_human.py`
Neither needs network access to re-run.

---

## 1. Provenance of every set

| Tag | Meaning |
|---|---|
| `MSIGDB_2026.1.Hs` | MSigDB **human** collections, release **2026.1.Hs**, fetched from `gsea-msigdb.org/gsea/msigdb/human/geneset/<NAME>.json` on **2026-08-27**. Raw JSON archived in `../msigdb_human_2026.1.Hs/`, together with the complete **H** collection bundle `_h.all.v2026.1.Hs.json` (50 sets) from `data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/`. |
| `CUR` | The mouse arm's curated lists (`../README.md`), translated to human symbols. No new curation was performed for the human arm. |
| `PHASE7_S10` | The literal gene list printed in Phase 7 §10. |
| `PANEL` | Derived from the panel's own `feature_type` column (Tier E controls). |
| `ORTHO` | `../mouse_human_orthologs_MGI.csv` — the map already pinned and used by `run_deepscence.py`, `caller_disagree2.py`, `run_phase3_n8.py`. |
| `HGNC` | `../hgnc_pin/` (downloaded 2026-08-27, md5 recorded) — alias / previous-symbol → Ensembl → panel-symbol resolution only. |

### Symbol resolution

Two hops, both from files, never from memory:

1. **mouse → human**: the MGI ortholog map. The Phase 7 §10 gotchas verify against it —
   `Trp53→TP53`, `Il6ra→IL6R`, `H2ax→H2AX`, `Ackr3→ACKR3`, plus `Trp53bp1→TP53BP1`,
   `Mre11a→MRE11`, `Mmp1a`/`Mmp1b→MMP1`, `Gpi1→GPI`, `Ftl1→FTL`, `Rab7→RAB7A`, `Tubb5→TUBB`.
2. **human symbol → panel symbol**: HGNC. The panel ships some legacy symbols, so a literal
   string match under-reports membership. Rule: *a symbol that is itself an approved HGNC symbol
   only ever resolves to its own Ensembl ID*; previous symbols and unambiguous aliases are used
   only for names HGNC does not recognise as approved. Without that rule the naive
   symbol/alias/prev union silently produced wrong genes (`CCNL1`→`MCM2`, `MIF`→`AMH`,
   `OPA1`→`MED12`) and injected four spurious Tier A × Tier B collisions.

Every resolution actually applied is logged in `_symbol_resolutions.csv` (21 rows).
**`H2AX` → `H2AFX`** (both `ENSG00000188486`) is the one that matters — see §2.

MSigDB human sets used, with systematic IDs and on-panel counts:

| Set | ID | collection | size | on-panel |
|---|---|---|---|---|
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | M5890 | H | 200 | 120 |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | M5897 | H | 87 | 71 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | M5913 | H | 200 | 113 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | M5911 | H | 97 | 45 |
| HALLMARK_E2F_TARGETS | M5925 | H | 200 | 112 |
| HALLMARK_G2M_CHECKPOINT | M5901 | H | 200 | 117 |
| HALLMARK_MYC_TARGETS_V1 | M5926 | H | 200 | 80 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | M5930 | H | 200 | 113 |
| HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | M5938 | H | 49 | 18 |
| REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES | M27244 | C2:CP:REACTOME | 37 | 15 |
| SAUL_SEN_MAYO (SenMayo) | M45803 | C2:CGP | 124 | 99 |
| REACTOME_SASP | M27187 | C2:CP:REACTOME | 111 | 27 |
| GOBP_CELLULAR_SENESCENCE | M11558 | C5:GO:BP | 109 | 61 |

(Full table, including the reference sets not used in any tier, in `_test2_summary.json`.)

---

## 2. Tier A — sender

Two definitions are shipped because they give **different Phase 7 §11 verdicts**:

| File | n on-panel | §11 gate |
|---|---|---|
| `A_PHASE7_S10_16` — the literal §10 sixteen | **14** | **FAIL** (14 < 15 before disjointness; 1 gene survives it) |
| `A_ported_pre_disjointness` (in `_test2_summary.json`) — mouse Tier A ported | 81 | PASS → `A_SENDER_FINAL_strict` n=21 |

**`H2AX`.** Not present under that symbol. It **is** on the panel as the legacy symbol
**`H2AFX`** (`ENSG00000188486`, resolved through HGNC, logged in `_symbol_resolutions.csv`).
It was **not** dropped.

**`TP53I3` and `HMGB1` are genuinely absent.** Checked three ways: literal symbol, HGNC approved
symbol → Ensembl (`ENSG00000115129`, `ENSG00000189403`) → panel `gene_id`, and a direct grep of
the panel file. Neither Ensembl ID appears among the 5,093. `TP53I3` also has no MGI mouse
ortholog record, so it was already a hand-added symbol in the mouse arm.

Files (post-removal subsets, mirroring the mouse layout):
`A_core_arrest` (4), `A_proliferation_down` (2), `A_nuclear_chromatin` (1),
`A_dna_damage_response` (6), `A_senescence_curated_nonsecreted` (8),
`A_SENDER_FINAL_strict` (21), `A_SENDER_FINAL_noB7` (33), `A_SENDER_FINAL_noB4_noB7` (75),
`A_sender_for_<module>` (36–81) — **use these**, as in the mouse arm —
`A_PHASE7_S10_16` (14), `A6_SenMayo_arrest_reference` (33).

The mouse arm's "no secreted factors" constraint was enforced by the mouse panel metadata's
`location` column. The human Prime 5K metadata CSV is not on disk, so the annotation was
**transferred through the ortholog map** (263 panel genes inherit `Secreted`). 0 Tier A members
carry it. This is a documented deviation.

---

## 3. Tier B — response

`B_tnfa_nfkb_proximal` (120), `B_il6_jak_stat3` (71), `B_interferon_response` (126),
`B_downstream_arrest` (231), `B_emt_ecm` (113), `B_oxidative_stress` (36),
`B_secondary_senescence` (35). All seven clear the ≥30 floor.

Source sets are **identical to the mouse arm's** so the two arms are comparable (Phase 7 test A8).
Phase 7 §10 names seven HALLMARK sets; `E2F_TARGETS` and `G2M_CHECKPOINT` both feed
`B4_downstream_arrest`, which is why seven named HALLMARK sets plus `secondary_senescence` make
**seven** modules, not eight.

**B6 needed the same rescue as mouse.** `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` alone gives
**18** on-panel and fails ≥30. The mouse arm's documented substitution (∪
`REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES` ∪ the curated NRF2 target list, ported)
gives **36**. The rescue is applied automatically and only when the Hallmark set alone misses the
floor — see `B6_hallmark_alone_on_panel` in `_test2_summary.json`.

**B7 is the mouse arm's curated list, ported, unfixed.** Phase 7 §23/C6 flags it. Overlap with
the Tier A caller on this panel: **23 of 35** against the ported Tier A, **11 of 35** against the
§10 sixteen (mouse reference: 14 of 38). Reported, not fixed — C6 is out of scope for Job A.

---

## 4. Tier C — ligand–receptor (`C_ligands` 13, `C_receptors` 15)

Same pairs as mouse, **plus `CXCL8`**.

| Ligand | on-panel | Receptors (on-panel) | Note |
|---|---|---|---|
| **`CXCL8`** | ✅ | `CXCR1` ✅ `CXCR2` ✅ | **§10 addition. No mouse ortholog** — and `CXCR1` has no clean mouse ortholog either. **HUMAN-ONLY quantity.** |
| `IL6` | ✅ | `IL6R` ✅ `IL6ST` ✅ | human receptor `IL6R`, mouse `Il6ra` |
| `CXCL1` | ❌ **absent** | `CXCR2` ✅ | was on the mouse panel |
| `CXCL2` | ✅ | `CXCR2` ✅ | |
| `CXCL5` | ✅ | `CXCR2` ✅ | |
| `CCL2` | ✅ | `CCR2` ✅ `ACKR3` ✅ | highest priority; fully covered |
| `CXCL12` | ✅ | `CXCR4` ✅ `ACKR3` ✅ `DPP4` ✅ | range-limiting mechanism testable |
| `TGFB1` | ✅ | `TGFBR1` ✅ `TGFBR2` ✅ | |
| `IL1A` | ✅ | `IL1R1` ✅ `IL1RAP` ✅ | membrane-bound → shortest λ |
| `IL1B` | ✅ | `IL1R1` ✅ `IL1RAP` ✅ | **on-panel here; it was OFF-panel in mouse** |
| `TNF` | ✅ | `TNFRSF1A` ✅ `TNFRSF1B` ✅ | |
| `IGFBP3` ✅ · `GDF15` ✅ · `THBS1` ✅ (`CD47` ✅, `SDC1` ❌) | | | |
| `IGFBP7` ❌ · `MMP1` ❌ · `MMP3` ❌ · `TIMP1` ❌ | | | off-panel |

**Cross-arm asymmetry, declared.** `CXCL8`/`CXCR1` exist only in the human arm; the mouse arm
substituted `Cxcl1`/`Cxcl2`/`Cxcl5` → `Cxcr2`. Any CXCL8 result is human-only and must not be
reported as replicating a mouse quantity. Conversely `CXCL1`, `MMP3` and `TIMP1` were on the
mouse panel and are **not** on the human one, so those three mouse ligands have no human
counterpart here.

**The §9 Tier C internal control is *stronger* here than in mouse:** `IL1B` is on-panel, so the
sharp version of the contrast — membrane-bound `IL1A` vs secreted `IL1B` through the *same*
receptor `IL1R1` — is runnable on H1 and was not runnable on M1.

---

## 5. Tier D — nuisance

`D_nuisance_covariates` is the §10 covariate list as **column names, not genes** (cell type,
density at 25/50/100 µm, k-NN composition, total counts, genes detected, cell and nucleus area,
section, segmentation method, distance to tissue boundary, plus the arm-specific anatomical
covariate).

`D_landmark_vessel_stroma` (14) is the mouse `D_landmark_vessel_bileduct` set minus the nine
liver-specific bile-duct/stellate genes, with the same over-adjustment guard applied (`ACTA2`,
`PDGFRB`, `SPP1` dropped for also being Tier A or Tier B members).

**The arm-specific anatomical covariate is DEFERRED, and this is a real gap.** The mouse
`D_zonation_*` sets were derived **from the expression matrix** (mouse README §6). H1 is spleen
and is still behind the §15 freeze — no H1 expression value has been read. Build the
red-pulp/white-pulp axis post-download with the identical procedure. **Phase 7 test A6 as written
specifies a lung airway-to-alveolar axis and does not apply to this arm.**

---

## 6. Tier E — controls

| File | n | Status |
|---|---|---|
| `E_housekeeping_planA` | **1/6** (`TBP`) | The §10 six are `ACTB GAPDH RPL13A RPS18 TBP PPIA`; five are absent, exactly as in mouse. **Unusable as written.** |
| `E_housekeeping_expanded` | **8/39** (`HPRT1 PGK1 RAB7A SDHA SNRPD3 TBP TUBB YWHAZ`) | Thinner than mouse (13). Report as such. |
| `E_negative_control_probes` | **40** | §12.4 / test A7. The real negative-control target probes. |
| `E_negative_control_codewords` | **609** | §12.4 / test A7. |
| `E_genomic_controls` | **21** | §12.4 / test A7. |
| E2 cell-type identity | **NOT BUILT** | In spleen every abundant cell type is an immune cell, so §10's "identity program unrelated to inflammation" has no valid instance. HGNC gene group *Hemoglobin subunits* puts only `HBG1`, `HBZ` on the panel. Choose it after Job B cell typing. |
| E3 random matched sets | **NOT BUILT** | `build_random_null_sets.py` needs the H1 expression matrix for the expression-matching bins. Post-freeze. |

Xenium panels deliberately omit high-expressors, which is why the housekeeping control fails in
both arms. Use the negative control probes (per-cell counts in `cells.parquet`) as the technical
null, binned, exactly as recommended for the mouse arm.

---

## 7. Section 11 gate

See `/workspace/reports/BIO_PHASE7_JobA.md` and
`/workspace/results/phase7_jobA/` (`gate_result_human.json`, `intersection_matrix_human.csv`,
`intersection_matrix_BxB.csv`, `onpanel_counts_by_tier.csv`, plus both run logs).
