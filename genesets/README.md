# Gene Set Package — SASP Spatial Response Kernel (Deliverable 2)

**Owner:** biology collaborator · **Date built:** 2026-08-20 · **Species: MOUSE**
**Target panel:** GSE310392 Xenium = *Xenium Prime Mouse 5K Pan Tissue & Pathways* (5,006 genes)
∪ GSE310392 custom 100-gene panel = **5,106 `Gene Expression` features**, verified feature-by-feature
against `cell_feature_matrix.h5` of GSM9295284 (union of the two panel files reproduces the h5
feature list exactly; 0 discrepancies either way). 9 of the custom 100 are genotyping probes
(`Jak2_WT_1940`, `Kras_ALT_235:GA`, `Pkd1_del_ex2-4_ALT_inv`, `Brca1_del_ex5-6_ALT_inv`, …) and are
**not** genes — exclude them from all scoring.

> **Species note.** These sets are MOUSE. The master plan Section 9 lists HUMAN symbols. Mouse
> symbols are Title-case, and several are not simple case conversions: `TP53`→`Trp53`,
> `TP53BP1`→`Trp53bp1`, `TP53I3`→`Trp53i3`, `IL6R`→`Il6ra`, `MMP1`→`Mmp1a`/`Mmp1b`.
> **Mouse has no `CXCL8`/IL-8 ortholog** — see Tier C.

---

## 1. Provenance of every set

| Tag | Meaning |
|---|---|
| `MSIGDB_2026.1.Mm` | Fetched live from the MSigDB JSON API (`gsea-msigdb.org/gsea/msigdb/mouse/geneset/<NAME>.json`) on **2026-08-20**. Release **2026.1.Mm**. These are MSigDB's *own* mouse collections — each set's `exactSource` field reads "Alliance Genome Consortium orthologs based on <human set>". **We did not hand-map orthologs.** Raw JSON archived in `msigdb_mouse_2026.1.Mm/`. |
| `CUR` | Curated by the biology collaborator from primary senescence literature, 2026-08-20. |
| `PANEL` | Derived from the panel metadata file's own annotation columns (`location`, `protein_name`, `cellchat_pathway`). |
| `DATA` | Derived from the real GSM9295284 expression matrix (zonation only — see §6). |

MSigDB sets used, with systematic IDs and on-panel counts:

| Set | ID | size | on-panel |
|---|---|---|---|
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | MM3860 | 196 | 126 |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | MM3866 | 85 | 68 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | MM3878 | 188 | 94 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | MM3877 | 94 | 35 |
| HALLMARK_E2F_TARGETS | MM3886 | 199 | 91 |
| HALLMARK_G2M_CHECKPOINT | MM3868 | 195 | 100 |
| HALLMARK_MYC_TARGETS_V1 | MM3887 | 197 | 62 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | MM3889 | 194 | 125 |
| HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | MM3895 | 48 | 17 |
| REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES | MM14944 | 34 | 11 |
| SAUL_SEN_MAYO (SenMayo) | MM16098 | 117 | 90 |
| REACTOME_SASP | MM14900 | 40 | 23 |
| GOBP_CELLULAR_SENESCENCE | MM10157 | 93 | 57 |

Section 9's instruction *"verify all MSigDB set sizes and membership against the current release"*
is **satisfied**: nothing here is recalled from memory. Note the plan's own hint that Hallmark sets
are "200 genes" does not hold for the mouse collections (ortholog mapping loses members).

---

## 2. Section 8 Test 2 — result

| Criterion | Required | Actual | Verdict |
|---|---|---|---|
| `len(A)` | ≥ 15 | **25** (union-strict) | **PASS** |
| `len(B[m])` per module | ≥ 30 | 126 / 68 / 100 / 190 / 125 / 31 / 38 | **PASS** (all 7) |
| `len(A ∩ B)` | == 0 | **0** | **PASS** |

### Tier A × Tier B intersection matrix (after enforcing disjointness by removal from Tier A)

| Tier A subset | B1 | B2 | B3 | B4 | B5 | B6 | B7 |
|---|---|---|---|---|---|---|---|
| A1_core_arrest (n=3) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A2_proliferation_down (n=2) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A3_nuclear_chromatin (n=3) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A4_dna_damage_response (n=7) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A5_senescence_curated (n=10) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **A_SENDER_FINAL_strict (n=25)** | **0** | **0** | **0** | **0** | **0** | **0** | **0** |

Pre-removal Tier A had 74 on-panel candidates. **49 were removed as collisions.** Which module
killed what:

| Module | n removed | genes |
|---|---|---|
| B4 downstream_arrest | 37 | Apex1 Aurkb Birc5 Brca1 Brca2 Bub1 Cbx5 Ccna2 Ccnd1 Cdc20 Cdc25a Cdk1 **Cdkn1a Cdkn1b Cdkn2a Cdkn2c** Chek1 Chek2 E2f1 Ezh2 Hmgb2 Lbr **Lmnb1** Mcm2 **Mki67** Mre11a Nbn Pcna Plk1 Rad50 Rbl1 Rpa1 Suv39h1 Top2a **Trp53** Ung Xrcc6 |
| B7 secondary_senescence | 19 | Atm Bax Bcl2l1 Ccnd1 Cdkn1a Cdkn1b Cdkn2a Cdkn2b Chek2 Ezh2 Gadd45a Glb1 Hmgb2 Lmnb1 Mdm2 Sirt1 Suv39h1 Trp53 Trp53bp1 |
| B1 tnfa_nfkb | 4 | Ccnd1 **Cdkn1a** Gadd45a Nfe2l2 |
| B6 oxidative_stress | 3 | Cdkn2d Nfe2l2 Sesn2 |
| B3 interferon_response | 1 | **Cdkn1a** |
| B5 emt_ecm | 1 | Gadd45a |

### ⚠️ Read this before using `A_SENDER_FINAL_strict`

The literal Section 8 rule ("A disjoint from the union of all B") technically passes, but the
surviving 25 genes are

```
Atr Bcl2 Ccnb1 Cdkn1c Ehmt2 Ercc1 Foxm1 Foxo3 H2afx Hmga1 Lmnb2 Mdm4 Nfatc1 Parp1
Plk3 Rad51 Rb1 Rbl2 Rif1 Sesn1 Terf2 Tert Trp63 Trp73 Xrcc5
```

— containing **no `Cdkn1a`, no `Cdkn2a`, no `Trp53`, no `Lmnb1`, no `Mki67`**, i.e. none of the
markers that actually define senescence, and it still contains `Ccnb1`/`Foxm1` whose direction is
opposite to the rest. **This is a numerically passing but biologically hollow sender score.**
Verified fact (MSigDB 2026.1.Mm, not memory): `Cdkn1a` is a member of
`HALLMARK_TNFA_SIGNALING_VIA_NFKB`, `HALLMARK_INTERFERON_GAMMA_RESPONSE`, and the E2F/G2M/MYC
union. The field's canonical senescence marker is inside the field's canonical inflammatory
response sets.

**Recommended design instead — per-module sender sets** (`A_sender_for_<module>.txt`).
Disjointness is only *statistically* required between the sender score and the **one** response
readout being fitted, not against the union of readouts you never pair it with:

| Response module | \|A_mod\| | ≥15? | canonical markers retained |
|---|---|---|---|
| B1 tnfa_nfkb_proximal | 70 | PASS | Cdkn2a Cdkn2b Trp53 Mki67 Lmnb1 Top2a Pcna Mdm2 Atm Glb1 Hmgb2 |
| B2 il6_jak_stat3 | 74 | PASS | + Cdkn1a (all 12) |
| B3 interferon_response | 73 | PASS | Cdkn2a Cdkn2b Trp53 Mki67 Lmnb1 Top2a Pcna Mdm2 Atm Glb1 Hmgb2 |
| B4 downstream_arrest | 37 | PASS | Cdkn2b Mdm2 Atm Glb1 |
| B5 emt_ecm | 73 | PASS | all 12 incl. Cdkn1a |
| B6 oxidative_stress | 71 | PASS | all 12 incl. Cdkn1a |
| B7 secondary_senescence | 55 | PASS | Mki67 Top2a Pcna |

Every per-module set clears ≥15 and every one is disjoint from its own readout.
Report `A_SENDER_FINAL_strict` as the conservative sensitivity analysis, not the primary.

---

## 3. Files

### Tier A — sender (arrest + damage only, NO secreted factors)
`A_core_arrest`, `A_proliferation_down`, `A_nuclear_chromatin`, `A_dna_damage_response`,
`A_senescence_curated_nonsecreted` (post-removal subsets, `CUR`)
`A_SENDER_FINAL_strict` (n=25, union-strict), `A_SENDER_FINAL_noB7` (n=33),
`A_SENDER_FINAL_noB4_noB7` (n=68), `A_sender_for_<module>` (n=37–74) — **use these**.
`A6_SenMayo_arrest_reference` (n=21) — the Section 9 A6 comparison set.

The Section 9 "no secreted factors" constraint was enforced *mechanically*, not by eye: a gene was
rejected from Tier A if the panel metadata `location` column contains `Secreted`. 0 Tier A members
carry that annotation. `Serpine1` excluded per Section 9.

### Tier B — response
`B_tnfa_nfkb_proximal` (126), `B_il6_jak_stat3` (68), `B_interferon_response` (100),
`B_downstream_arrest` (190), `B_emt_ecm` (125), `B_oxidative_stress` (31),
`B_secondary_senescence` (38).

**Documented deviation, B6.** `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` alone gives only **17**
on-panel genes and **FAILS** the ≥30 bar. Substitution: ∪ `REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES`
∪ a curated NRF2-target list → **31**, with only 6 genes shared with B1–B5. Rejected alternatives:
`GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS` (149 on-panel but 33 shared with B1–B5 — too diluted).

**Caveat, B7.** Section 9 defines B7 as "Tier A minus the calling genes", so it is *not* independent
of Tier A by construction. Recommended: split Tier A into disjoint halves `A_call` / `A_readout` by
random assignment with a fixed seed, and report the split-half stability. The B7 file here is a
curated secondary-senescence list, not that split.

### Tier C — ligand–receptor (`C_ligands`, `C_receptors`)

| Ligand | on-panel | Receptors (on-panel) | Note |
|---|---|---|---|
| `Il6` | ✅ | `Il6ra` ✅ `Il6st` ✅ | mouse receptor is **Il6ra**, not IL6R |
| **`Cxcl8`** | **— (no mouse ortholog)** | — | **Mouse has NO CXCL8/IL-8 gene.** Functional analogues below. |
| `Cxcl1` (KC) | ✅ | `Cxcr2` ✅ | CXCL8 analogue |
| `Cxcl2` (MIP-2) | ✅ | `Cxcr2` ✅ | CXCL8 analogue |
| `Cxcl5` | ✅ | `Cxcr2` ✅ | CXCL8 analogue |
| `Ccl2` | ✅ | `Ccr2` ✅ `Ackr3`(CXCR7) ✅ | Section 9 highest priority — **fully covered** |
| `Cxcl12` | ✅ | `Cxcr4` ✅ `Ackr3` ✅ `Dpp4` ✅ | range-limiting mechanism testable; paper reports Cxcl12–Cxcr4 conserved in mouse |
| `Tgfb1` | ✅ | `Tgfbr1` ✅ `Tgfbr2` ✅ | short λ expected |
| `Il1a` | ✅ | `Il1r1` ✅ `Il1rap` ✅ | membrane-bound → **shortest λ**; anchor of the Section 9 internal control |
| `Il1b` | ❌ **absent** | `Il1r1` ✅ `Il1rap` ✅ | see below |
| `Tnf` | ✅ | `Tnfrsf1a` ✅ `Tnfrsf1b` ✅ | drives B1 |
| `Igfbp3` ✅ · `Gdf15` ✅ · `Mmp3` ✅ · `Timp1` ✅ · `Thbs1` ✅ (`Cd47`,`Sdc1` ✅) | | | SASP Atlas / paper-specific |
| `Igfbp7` ❌ · `Mmp1a`/`Mmp1b` ❌ | | | off-panel |

**The Section 9 Tier C internal control survives.** The ordering it needs
(membrane-bound `Il1a` shortest → diffusible chemokines longest) is testable: `Il1a`+`Il1r1` are
on-panel, and so are `Ccl2`/`Cxcl12`/`Cxcl1`/`Cxcl2` with all their receptors.
**Weakened by:** `Il1b` is off-panel, so the sharpest version of the contrast — membrane-bound
Il1a vs. secreted Il1b through the *same* receptor `Il1r1` — cannot be run. Substitute the
Il1a-vs-Cxcl1/Ccl2 contrast and state the limitation.

`PANEL`-annotated CellChat pathway members are also available and were cross-checked
(CXCL: 8 receptors; CCL: 12; IL1: 5; TNF: 2; TGFb: 6; THBS: 7).

### Tier E — negative controls
`E_housekeeping_planA` — **1/6 on-panel** (`Tbp` only; `Actb` `Gapdh` `Rpl13a` `Rps18` `Ppia` all
absent). Xenium panels deliberately omit high-expressors to avoid optical crowding. **This control
as written in Section 9 is unusable.**
`E_housekeeping_expanded` — **13/39** (`Actr2 Cnot4 Gpi1 Hnrnpk Hprt Pgk1 Rab7 Rack1 Rpl11 Sdha Tbp Vcp Ywhaz`).
Usable as the per-cell flat-kernel control. n=13 is thin; report it as such.
`E_hepatocyte_identity` — **9/32** (`Alb`, `Ttr`, `Apoa1`… `Alb` is ABSENT from a liver panel).
**⚠️ DO NOT USE as a flat control** — see §6: its score correlates with the zonation axis at
**r = −0.67**, worse than 99% of random size-matched panel sets. In liver, "cell-type identity
program" *is* a zonation readout. This is a genuine defect in Section 9's Tier E design for liver.

**Substitute technical control (recommended, deviation from Section 9).** `cells.parquet` carries
per-cell `control_probe_counts`, `genomic_control_counts`, `control_codeword_counts`,
`unassigned_codeword_counts`, `deprecated_codeword_counts`. These measure optical/decoding
background directly and are a *better* technical null than housekeeping genes.
**Constraint:** they are extremely sparse — summed over 237,982 sham cells, non-deprecated control
counts total 41,614, i.e. ~0.17 counts/cell, so a **per-cell** control score is ~0 and cannot be
regressed per cell. Use them **binned** (e.g. 50 µm hexbins) as a spatial-gradient control, and use
`E_housekeeping_expanded` for the per-cell flat-kernel test. Report both.

### Tier D — zonation and landmarks
`D_zonation_pericentral` (20), `D_zonation_periportal` (22), `D_landmark_vessel_bileduct` (22).

---

## 4. Handoff notes for the CS lead

- **Coordinate units are microns.** Name every derived column with the unit (`x_um`, `y_um`,
  `dist_to_portal_um`). Section 13.
- Drop the 9 genotyping probes before any scoring.
- `segmentation_method` in `cells.parquet` takes 3 values and is spatially patterned — add it to
  Tier D nuisance covariates.
- Use the per-module sender sets as primary; `A_SENDER_FINAL_strict` as sensitivity.

## 5. Reproduce

```bash
python3 /workspace/code/build_genesets.py
```
Needs `/workspace/results/zonation_gene_correlations_7250_sham.csv` (§6) and the archived MSigDB
JSON in `msigdb_mouse_2026.1.Mm/`. No network access required to re-run.

## 6. Zonation derivation (`DATA`)

The plan's Section 11 marker list is **human** and largely **off-panel**: periportal `Ass1` ❌
`Sds` ❌ `Alb` ❌ (`Hal` ✅ `Cps1` ✅); pericentral `Oat` ❌ (`Glul` ✅ `Cyp2e1` ✅ `Cyp1a2` ✅).
Rebuilt from data:

1. Seed axis on GSM9295284 (sham), log1p CP-median normalised, cells with ≥20 counts (99.5%):
   `mean(Glul,Cyp2e1,Cyp1a2,Cyp2a5,Cyp27a1) − mean(Hal,Cps1,Arg1,Otc,Hamp,Igfbp2)`.
2. Restricted to the top 60% of cells by a 19-gene hepatocyte marker score (n = 142,143).
3. Pearson r of every panel gene against the axis. Retain |r| ≥ 0.20, detection ≥ 5%.
4. **Over-adjustment guard:** any gene also in Tier A or Tier B is dropped, so conditioning on
   zonation does not partially condition on the response. (Dropped: `Cfh`, `Igfbp2` periportal;
   `Acta2 Col1a1 Col1a2 Pdgfrb Spp1 Tagln` from the landmark set.)

Full per-gene table: `/workspace/results/zonation_gene_correlations_7250_sham.csv`.

**Validation:** the resulting pericentral score correlates with the seed axis at r = **+0.897** and
the periportal score at r = **−0.737**, versus a null of 300 random size-60 on-panel sets
(mean|r| = 0.114, p90 = 0.261). The axis is real and strong, and recovers textbook zonation without
being told it: `Cyp1a2` +0.77, `Glul` +0.76, `Cyp2e1` +0.74 pericentral; `Arg1` −0.76, `Cps1` −0.57,
`Hal` −0.49 periportal. Newly recovered on-panel markers not in the plan's list: pericentral
`Cyp2c29 Cyp3a11 Pon1 Lect2 Cyb5a Sult2a8 Glud1 Por Prodh Nrn1 C6 Cyp7a1 Cyp7b1 Tbx3 Hpd`;
periportal `Pigr Apoa4 Acly Etnppl Cp Pklr Sfxn1 Hamp Ak4 Gm2a Hc Pzp Vtn Orm1`.

**Caveats.** (i) Derived from ONE section (sham); recompute per section and check stability before
using across samples. (ii) The seed genes are inside their own scores — the +0.897/−0.737 figures
are partly circular; the non-seed discoveries are the real evidence. (iii) "Hepatocyte" here is a
marker-score threshold, not a clustered annotation — replace with real cell-type calls when
Deliverable 4 lands. (iv) `Pigr` is also a cholangiocyte marker; its strong periportal loading may
partly reflect bile-duct contamination of the hepatocyte gate.
