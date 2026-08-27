# BIO PHASE 7 / Job A — Human Panel Annotation, Tiers A–E, and the §11 Disjointness Gate

**Biology collaborator · 2026-08-27 · SASP Spatial Response Kernel**
Phase 7 Part III (§10, §11). Human arm H1 = GEO **GSE326743**, 7 human spleens (FFPE),
Xenium Prime 5K Human + 100 addon, **5,093 `Gene Expression` features**.
Deliverables: `genesets/human/` (33 gene/covariate files + README), `code/build_genesets_human.py`,
`code/gate_disjointness_human.py`, `results/phase7_jobA/`.

> **SUPERSEDED IN PART — corrections, 2026-08-27.** The PI adopted the re-sourced B7 (C6) and
> froze Tier A. Three numbers below moved as a result and the current values are in
> `reports/PREREG_PHASE8_genesets.md` and `reports/BIO_PHASE8_FREEZE.md`:
> `B_secondary_senescence` **35 → 116**; `A_SENDER_FINAL_strict` **21 → 33**; the §10 Tier A now
> loses **9** genes rather than 13 and leaves **5** survivors (`ATM ATR CDKN2B MDM2 TP53BP1`), not
> `ATR` alone — it still fails the gate. The CoreScence circularity figure moved **76% → 88%**;
> cite 88%. Everything else in this report stands.
>
> **Follow-on block: `reports/BIO_PHASE7_JobA_FOLLOWON.md`** — spleen cell-type markers, the A6
> red-pulp/white-pulp covariate specification, the C6 B7 rebuild, the Tier E2 decision, SenePy's
> spleen coverage, and the Tier A decision memo. Two results there change numbers in this report:
> the C6-sourced B7 raises `A_SENDER_FINAL_strict` from 21 to 33 genes, and SenePy ships no spleen
> signature at all.

---

## 0. Headline

1. **The §11 gate FAILS on Tier A exactly as §10 defines it, and it fails twice over.** The 16
   symbols printed in §10 put **14** genes on the GSE326743 panel — below the `len(A) >= 15`
   floor **before disjointness is even tested**. Then 13 of those 14 collide with a Tier B module,
   and the §11 remedy ("remove the overlapping genes from Tier A") leaves **`ATR`, alone**. §2.
2. **The gate PASSES on the ported mouse Tier A**, which is the like-for-like port Job A was
   actually asked to build: 81 on-panel → 60 removed as collisions → **21** survive, disjoint from
   all seven modules. All seven Tier B modules clear ≥30 (120 / 71 / 126 / 231 / 113 / 36 / 35).
   §2, §3.
3. **CoreScence is circular on both arms at a similar rate: 25 of 33 on-panel genes = 76%**
   are members of ≥1 Tier B module here, against **26 of 33 = 79%** on the pre-C6 mouse arm. And
   this is the *native* number — CoreScence is a human gene set, so on H1 there is no ortholog
   remapping to blame it on. §4. **Corrected 2026-08-27:** this bullet read "mouse arm: 24/35 =
   69%", and therefore claimed human was *more* circular than mouse. The mouse denominator 35 was
   a typed-in literal reproducible under no mapping convention; derived from files it is 33
   (`code/corescence_circularity.py`), the human arm is if anything *slightly less* circular
   pre-C6, and the direction of that comparison does not survive. Under the frozen C6 sets the two
   arms are identical at 29/33 = 88%.
4. **`H2AX` was not dropped.** It is on the panel under the legacy symbol **`H2AFX`**
   (`ENSG00000188486`). `TP53I3` and `HMGB1` are genuinely absent — verified by Ensembl ID, not by
   symbol match. §1.2.
5. **Two things will break the human arm as currently planned.** (i) Phase 7 test **A6 specifies a
   lung airway-to-alveolar axis**; H1 is **spleen**, so the arm-specific anatomical covariate has
   to be a red-pulp/white-pulp axis, and it cannot be built before the freeze because the mouse
   equivalent was derived from the expression matrix. (ii) Tier E's "cell-type identity program
   unrelated to inflammation" has **no valid instance in spleen** — every abundant cell type there
   is an immune cell. §5, §6.

MSigDB pinned at release **2026.1.Hs**, fetched **2026-08-27**, raw JSON archived in
`genesets/msigdb_human_2026.1.Hs/` (26 per-set files in the mouse archive's format, plus the
complete **H** collection bundle `_h.all.v2026.1.Hs.json`, 50 sets).

---

## 1. What was built

`genesets/human/` mirrors `genesets/` file-for-file and format-for-format (one symbol per line,
same `A_`/`B_`/`C_`/`D_`/`E_` prefixes, same `_test2_summary.json`), so downstream code needs no
change beyond the directory. Full per-tier documentation: `genesets/human/README.md`.

```bash
python3 /workspace/code/build_genesets_human.py      # builds genesets/human/
python3 /workspace/code/gate_disjointness_human.py   # runs §11, writes results/phase7_jobA/
```

Both re-run offline. Logs of the runs that produced every number below:
`results/phase7_jobA/build_genesets_human.log`, `results/phase7_jobA/gate_disjointness_human.log`.

### 1.1 Sources — nothing here is recalled from memory

| Tag | Source |
|---|---|
| `MSIGDB_2026.1.Hs` | `gsea-msigdb.org/gsea/msigdb/human/geneset/<NAME>.json`, fetched 2026-08-27; plus `data.broadinstitute.org/.../release/2026.1.Hs/h.all.v2026.1.Hs.json`. Archived raw. |
| `PHASE7_S10` | the literal lists printed in Phase 7 §10 |
| `CUR` | the mouse arm's curated lists (`genesets/README.md`), translated — no new curation |
| `PANEL` | the panel's own `feature_type` column (Tier E controls) |
| `ORTHO` | `genesets/mouse_human_orthologs_MGI.csv`, already pinned and used by `run_deepscence.py` |
| `HGNC` | `genesets/hgnc_pin/` (downloaded 2026-08-27, md5 `2d741e796d5538cc48d3696452237781`) |

### 1.2 Symbol resolution, and the `H2AX` question

Resolution is two hops, both from files. Mouse→human uses the MGI map; the §10 gotchas all verify
against it (`Trp53→TP53`, `Il6ra→IL6R`, `H2ax→H2AX`, `Ackr3→ACKR3`), as do `Trp53bp1→TP53BP1`,
`Mre11a→MRE11`, `Mmp1a`/`Mmp1b→MMP1`, `Gpi1→GPI`, `Ftl1→FTL`, `Rab7→RAB7A`, `Tubb5→TUBB`.
Human symbol→panel symbol uses HGNC, because the panel ships some legacy symbols and a literal
string match under-reports membership.

**`H2AX` → `H2AFX`.** The panel carries `H2AFX,ENSG00000188486`. HGNC's approved symbol for
`ENSG00000188486` is `H2AX`, with `H2AFX` as its previous symbol. The build resolves `H2AX`
through the Ensembl ID and writes `H2AFX` — the symbol the data will actually use. It is in
`A_dna_damage_response`'s input, in `A_PHASE7_S10_16`, and it is one of the genes the gate later
removes for colliding with B4. Logged in `genesets/human/_symbol_resolutions.csv`.

**`TP53I3` and `HMGB1` are absent, checked three ways** — literal symbol, HGNC approved symbol →
Ensembl (`ENSG00000115129`, `ENSG00000189403`) → panel `gene_id`, and a direct grep of the panel
CSV. Neither Ensembl ID is among the 5,093. `TP53I3` additionally has no MGI mouse-ortholog record
at all, so it was already a hand-added symbol on the mouse side.

**One defect found and fixed during the build.** A naive "symbol ∪ alias ∪ previous symbol → Ensembl"
index silently produced *wrong genes*: `CCNL1`→`MCM2`, `MIF`→`AMH`, `MBP`→`MBL2`, `OPA1`→`MED12`,
`SFN`→`REXO2` and seven more, because those names appear in some other gene's alias list. That
injected four spurious Tier A × Tier B collisions (`MCM2` into B1 among them) and inflated three
module sizes. The rule that fixes it: **a symbol that is itself an approved HGNC symbol only ever
resolves to its own Ensembl ID**; aliases are consulted only for names HGNC does not recognise as
approved, and only when unambiguous. Anyone porting tiers to a third species must apply the same
rule — the failure is silent and it moves the collision matrix.

---

## 2. §11 — the disjointness gate

Assertions run verbatim as printed in §11, against **both** candidate Tier A definitions, because
they disagree. Full output: `results/phase7_jobA/gate_disjointness_human.log`;
machine-readable: `gate_result_human.json`.

### 2.1 Tier A = the literal §10 sixteen — **FAIL**

| Assertion | Required | Actual | Verdict |
|---|---|---|---|
| `len(A & panel)` | ≥ 15 | **14** | **FAIL** |
| `len(B_k & panel)` ×7 | ≥ 30 | 120 / 71 / 126 / 231 / 113 / 36 / 35 | PASS (all 7) |
| `len(A & B_k) == 0` ×7 | 0 | **23 gene memberships** | **FAIL** |
| re-check `len(A) >= 15` after removal | ≥ 15 | **1** | **FAIL** |

Off-panel before anything else: `TP53I3`, `HMGB1`.

Genes removed from Tier A by the §11 remedy, and the module that removed each:

| Module | n | genes removed from Tier A |
|---|---|---|
| B1 tnfa_nfkb_proximal | 2 | `CDKN1A` `GADD45A` |
| B3 interferon_response | 1 | `CDKN1A` |
| B4 downstream_arrest | 8 | `CDKN1A` `CDKN2A` `CHEK1` `CHEK2` `H2AFX` `HMGB2` `LMNB1` `TP53` |
| B5 emt_ecm | 1 | `GADD45A` |
| B7 secondary_senescence | 11 | `ATM` `CDKN1A` `CDKN2A` `CDKN2B` `CHEK2` `GADD45A` `HMGB2` `LMNB1` `MDM2` `TP53` `TP53BP1` |

**Unique genes removed (13):** `ATM CDKN1A CDKN2A CDKN2B CHEK1 CHEK2 GADD45A H2AFX HMGB2 LMNB1
MDM2 TP53 TP53BP1`. **Surviving Tier A: `ATR`.** A one-gene sender score is not a sender score.

This is not a panel problem. §11 predicts "a 5K panel clears this comfortably", and the *panel*
does — the modules are large and healthy. What fails is that §10's Tier A **is** the canonical
senescence marker list, and the canonical senescence markers are inside the canonical response
sets. Same finding as the mouse arm (`genesets/README.md` §2), reproduced independently in human,
and sharper because §10's Tier A has no slack: 16 genes cannot lose 13 and clear a floor of 15.

### 2.2 Tier A = the ported mouse Tier A — **PASS**

| Assertion | Required | Actual | Verdict |
|---|---|---|---|
| `len(A & panel)` | ≥ 15 | **81** | PASS |
| `len(B_k & panel)` ×7 | ≥ 30 | 120 / 71 / 126 / 231 / 113 / 36 / 35 | PASS (all 7) |
| `len(A & B_k) == 0` ×7 | 0 | 77 memberships → **60 unique genes removed** | enforced by removal |
| re-check `len(A) >= 15` after removal | ≥ 15 | **21** | PASS |

Removals by module: B1 4, B3 1, **B4 45**, B5 2, B6 2, **B7 23** (full lists in the log and in
`intersection_matrix_human.csv`). Surviving strict Tier A (n=21):

```
ATR BCL2 CCNB1 CDKN1C EHMT2 ERCC1 FOXM1 FOXO3 GADD45G MDC1 MDM4 NFATC1
PARP1 RAD51 RB1 RBL2 TERF2 TERT TP63 TP73 XRCC5
```

Same caveat as mouse: numerically passing, biologically hollow — no `CDKN1A`, `CDKN2A`, `TP53`,
`LMNB1`, `MKI67`; and `CCNB1`/`FOXM1` point the wrong way. **Use the per-module sender sets as
primary**, as the mouse arm concluded. Disjointness is only *statistically* required between the
sender score and the one response readout being fitted:

| Response module | \|A_mod\| | ≥15? | canonical markers retained |
|---|---|---|---|
| B1 tnfa_nfkb_proximal | 77 | PASS | CDKN2A CDKN2B TP53 MKI67 LMNB1 TOP2A PCNA MDM2 ATM GLB1 HMGB2 |
| B2 il6_jak_stat3 | 81 | PASS | + CDKN1A (all 12) |
| B3 interferon_response | 80 | PASS | CDKN2A CDKN2B TP53 MKI67 LMNB1 TOP2A PCNA MDM2 ATM GLB1 HMGB2 |
| B4 downstream_arrest | 36 | PASS | CDKN2B MDM2 ATM GLB1 |
| B5 emt_ecm | 79 | PASS | all 12 incl. CDKN1A |
| B6 oxidative_stress | 79 | PASS | all 12 incl. CDKN1A |
| B7 secondary_senescence | 58 | PASS | MKI67 TOP2A PCNA |

Every per-module set clears ≥15 and every one is disjoint from its own readout. The pattern is
almost identical to mouse (mouse: 70/74/73/37/73/71/55).

### 2.3 Recommendation

Report `A_PHASE7_S10_16` as the **pre-registered definition that failed**, with the removal table
above — that is exactly the reporting-against-interest the mouse paper's credibility rests on.
Run the analysis on the per-module sender sets, with `A_SENDER_FINAL_strict` (n=21) as the
conservative sensitivity analysis. Do **not** quietly widen §10's Tier A until it passes and then
present it as the plan.

---

## 3. The full intersection matrix (§11: "print it in the Methods")

`results/phase7_jobA/intersection_matrix_human.csv` (counts **and** the gene lists per cell);
`intersection_matrix_BxB.csv`; `onpanel_counts_by_tier.csv`.

### 3.1 Tier A × Tier B, on-panel gene counts

| Tier A set | n | B1 | B2 | B3 | B4 | B5 | B6 | B7 | total |
|---|---|---|---|---|---|---|---|---|---|
| A_core_arrest | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A_proliferation_down | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A_nuclear_chromatin | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A_dna_damage_response | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A_senescence_curated_nonsecreted | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **A_PHASE7_S10_16** | **14** | **2** | **0** | **1** | **8** | **1** | **0** | **11** | **23** |
| A_ported_pre_disjointness | 81 | 4 | 0 | 1 | 45 | 2 | 2 | 23 | 77 |
| A_SENDER_FINAL_noB4_noB7 | 75 | 0 | 0 | 0 | 42 | 0 | 0 | 20 | 62 |
| A_SENDER_FINAL_noB7 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 12 | 12 |
| **A_SENDER_FINAL_strict** | **21** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** |

B1 tnfa_nfkb_proximal · B2 il6_jak_stat3 · B3 interferon_response · B4 downstream_arrest ·
B5 emt_ecm · B6 oxidative_stress · B7 secondary_senescence.

### 3.2 Tier B × Tier B (module cross-talk)

| module | n | B1 | B2 | B3 | B4 | B5 | B6 | B7 |
|---|---|---|---|---|---|---|---|---|
| B1 tnfa_nfkb_proximal | 120 | 120 | 15 | 23 | 5 | 15 | 2 | 10 |
| B2 il6_jak_stat3 | 71 | 15 | 71 | 22 | 1 | 6 | 1 | 3 |
| B3 interferon_response | 126 | 23 | 22 | 126 | 3 | 5 | 1 | 3 |
| B4 downstream_arrest | 231 | 5 | 1 | 3 | 231 | 3 | 0 | 11 |
| B5 emt_ecm | 113 | 15 | 6 | 5 | 3 | 113 | 0 | 6 |
| B6 oxidative_stress | 36 | 2 | 1 | 1 | 0 | 0 | 36 | 0 |
| B7 secondary_senescence | 35 | 10 | 3 | 3 | 11 | 6 | 0 | 35 |

§11 does not require the modules to be disjoint from *each other*, and they are not: B1∩B3 = 23,
B1∩B2 = 15, B1∩B5 = 15. Same structure as mouse. Worth one sentence in Methods so a reviewer does
not read seven modules as seven independent tests.

### 3.3 On-panel counts per tier

| Tier | files | on-panel |
|---|---|---|
| A | 12 gene files | 1–81 (see table above); `A6_SenMayo_arrest_reference` 33 |
| B | 7 modules | 120 / 71 / 126 / 231 / 113 / 36 / 35 |
| C | `C_ligands` 13, `C_receptors` 15 | |
| D | `D_landmark_vessel_stroma` 14; `D_nuisance_covariates` 13 covariate names (not genes) | |
| E | `E_housekeeping_planA` 1, `E_housekeeping_expanded` 8; controls: 40 probes / 609 codewords / 21 genomic | |

---

## 4. CoreScence circularity — comparable on the two arms

§11 flags that "CoreScence was 69% circular with your response modules last time". That mouse
number, 24 of 35 on-panel CoreScence genes (`BIO_PHASE2.md` §4.2, audited as claim C15), **does
not reproduce** — see the correction box below; the mouse figure is **26 of 33 = 79%**.
Recomputed natively on H1 — **CoreScence is a human gene set, so on this arm it runs with no
ortholog remapping at all** (Phase 7 §7 D-c):

| Arm | CoreScence occ≥5 | on-panel | in ≥1 Tier B module | fraction |
|---|---|---|---|---|
| M1 mouse, pre-C6 (reference, after mouse remapping) | 39 | **33** | **26** | **79%** |
| **H1 human**, superseded B7 | 39 | **33** | **25** | **76%** |
| M1 mouse, frozen C6 B7 | 39 | 33 | 29 | **88%** |
| H1 human, frozen C6 B7 | 39 | 33 | 29 | **88%** |

> **Correction, 2026-08-27** (`reports/AUDIT_PHASE8_FACTCHECK.md` M1). The mouse row read
> `35 / 24 / 69%`, sourced from `BIO_PHASE2.md` §4.2. 35 is reproducible under no mapping
> convention: 31 of the 39 CoreScence genes are reachable on the 5,097-gene mouse panel through
> the pinned MGI map, 33 with the project's documented Title-case fallback (`CDKN2B`, `CXCL1`),
> and 6 are off the mouse panel entirely. 24/31 = 77% strict; **26/33 = 79%** under the
> convention `run_phase3_n8.py::corescence_mouse` implements and the committed
> `results/phase3/n8_disjointness_*.csv` records (`corescence_on_panel = 33`, per-module
> 10/9/5/5/0/14/8). Derived by `code/corescence_circularity.py`.

Circular genes: `BRCA1 BUB1B CCL2 CCNA2 CDK1 CDKN1A CDKN2A CDKN2B CXCL8 FAS FGF2 GDF15 HELLS
HMGB2 ICAM1 IGFBP2 IGFBP3 IL1A IL6 LMNB1 MDM2 SERPINE1 STAT1 TGFB1 VEGFA`.

Per module, H1 human under the superseded B7 (fraction of the 33 on-panel CoreScence genes):
B7 0.39, B4 0.30, B5 0.27, B1 0.21, B3 0.18, B2 0.12, B6 0.00. Under the frozen C6 B7 the human
per-module fractions are B7 0.55, B4 0.30, B5 0.27, B1 0.21, B3 0.18, B2 0.12, B6 0.00 and the
mouse ones B7 0.55, B4 0.30, B5 0.27, B1 0.24, B2 0.15, B3 0.15, B6 0.00
(`results/phase7_jobA/corescence_circularity_mouse.json`,
`figures/figure_gs3_corescence_circularity_data.csv`).

**So the answer to "is the same true here" is: yes, at essentially the same rate.** The mouse
figure could not be blamed on the mouse→human remapping, and the human number confirms it: 76%
natively against a corrected mouse 79% pre-C6, and 88% on both arms under the frozen C6 sets.
(The earlier claim that the human arm was *higher* rested on the irreproducible 69%.)

Since DeepScence is the sender caller for Job B, the phase-3 circularity
sensitivity (fit the module with the shared genes stripped, report both amplitudes) must be
repeated on H1, not inherited.

**B7 / §23-C6, reported not fixed.** `B_secondary_senescence` (35 on-panel) shares **23 genes**
with the ported Tier A caller and **11** with the §10 sixteen (mouse reference: 14 of 38). It is
the single largest contributor to the §10 Tier A gate failure after B4. Rebuilding it is C6 and
out of scope for Job A.

---

## 5. Every deviation from the mouse tiers, with its reason

| # | Deviation | Reason |
|---|---|---|
| D1 | Tier A shipped in **two** definitions (`A_PHASE7_S10_16` n=14, ported n=81) | They give opposite §11 verdicts. Hiding either would be reporting toward interest. §2. |
| D2 | `H2AX` written as **`H2AFX`** | The panel uses the legacy symbol for `ENSG00000188486`. Resolved through HGNC, logged; not dropped. §1.2. |
| D3 | `TP53I3`, `HMGB1` absent from Tier A | Genuinely not on the panel (checked by Ensembl ID). No substitute invented. |
| D4 | "No secreted factors" enforced by **transferring** the mouse panel's `location` annotation through the ortholog map | The human Prime 5K panel-metadata CSV is not on disk. 263 panel genes inherit `Secreted`; 0 Tier A members do. |
| D5 | B3 keeps `INTERFERON_ALPHA_RESPONSE`, B4 keeps `MYC_TARGETS_V1` — neither is named in §10 | Mouse-identical source sets, for cross-arm comparability (test A8). §10 names 7 HALLMARK sets; E2F+G2M both feed B4, so 7 names + secondary_senescence = **7** modules, not 8. Dropping the two extras would change B3 to 113 and B4 to 178 and make the arms non-comparable. |
| D6 | B6 rescued the same way as mouse | `HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY` alone = **18** on-panel, fails ≥30. ∪ `REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES` ∪ ported NRF2 targets = **36**. Applied automatically and only on failure. |
| D7 | **Tier C gains `CXCL8` + `CXCR1`** | §10. No mouse ortholog for either. **Human-only quantity — must not be reported as replicating a mouse result.** |
| D8 | Tier C loses `CXCL1`, `MMP3`, `TIMP1`, and `SDC1` | On the mouse panel, off the human one. Mouse-only quantities. |
| D9 | Tier C gains `IL1B` | Off-panel in mouse, on-panel here. The §9 internal control (membrane-bound `IL1A` vs secreted `IL1B` through the same receptor `IL1R1`) is **runnable on H1 and was not runnable on M1**. |
| D10 | `D_landmark_vessel_bileduct` → `D_landmark_vessel_stroma` (14) | Nine bile-duct/stellate genes are liver-specific and were not ported. Same over-adjustment guard applied (`ACTA2`, `PDGFRB`, `SPP1` dropped). |
| D11 | **No anatomical covariate gene set** | The mouse `D_zonation_*` sets were derived from the mouse **expression matrix**. H1 is spleen and behind the §15 freeze. Build the red-pulp/white-pulp axis post-download with the identical procedure. |
| D12 | **Tier E2 (cell-type identity) not built** | §10 asks for an identity program "unrelated to inflammation". In spleen every abundant cell type is an immune cell. HGNC gene group *Hemoglobin subunits* puts only `HBG1`, `HBZ` on the panel — 2 genes, unusable. Choose after Job B. |
| D13 | **Tier E3 (500 random matched sets) not built** | `build_random_null_sets.py` bins genes by mean expression on the real matrix. Post-freeze. |
| D14 | Tier E4 controls added as three files (40 / 609 / 21) | §12.4 and test A7 need them; they have no mouse counterpart file. |

Nothing was installed with pip. The build uses only `csv`, `gzip`, `json`, `os`, `glob`,
`collections` from the standard library — no numpy, pandas, scanpy or h5py — so the pinned
environment was neither touched nor exercised by this job.

---

## 6. What looks like it will break the human arm

1. **Test A6 is written for the wrong tissue.** §13 A6 specifies a lung airway-to-alveolar axis
   plus fibroblastic-focus distance. H1 is **spleen**. The covariate has to be red pulp vs white
   pulp (plus distance to the nearest large vessel / capsule), and §13's pass condition and the
   §17 cross-arm table both need rewording before the freeze. This is a **pre-registration
   document edit, and it must happen before §15**, not after H1 is downloaded.
2. **The §10 Tier A cannot be the pre-registered sender definition** — it fails A2/§11 outright.
   Decide between (a) the per-module sender design and (b) a widened Tier A, **and write the
   decision into the freeze**, because deciding after seeing H1 is exactly the retuning risk §22
   names as medium-high.
3. **Tier E has no usable per-cell flat control on this panel.** `E_housekeeping_planA` is 1/6
   (`TBP` only) and `E_housekeeping_expanded` is 8/39, thinner than mouse's 13. The negative
   control probes (40) are the better control and A7 already requires them, but they are sparse
   and will need binning, as in mouse. Plan on the control probes, not on housekeeping genes.
4. **DeepScence's circularity carries into Job B at 76%.** The sender caller for H1 shares three
   quarters of its on-panel gene set with the response modules it will be regressed against. The
   strip-and-refit sensitivity has to be part of the frozen run order, not an afterthought.
5. **Three Tier C ligands are arm-exclusive in each direction** (`CXCL8`/`CXCR1` human-only;
   `CXCL1`/`MMP3`/`TIMP1` mouse-only). Test A8 asks for every cross-arm number on the
   ortholog-intersected panel **and** each arm's full panel — for Tier C that intersection is
   materially smaller than either arm, and the λ-ordering comparison must be run on it explicitly.
6. **`GSE326743` is a normal-spleen ageing series (ages 17–59), not an injury model.** Nothing in
   Job A depends on it, but §13 A3 requires 1–20% sender prevalence with ≥200 senders and ≥5,000
   non-senders. That is the next go/no-go and it has not been tested — no H1 expression value has
   been read, per the §15 freeze.

---

## 7. Files

```
genesets/human/                       33 .txt files, mouse layout, + README.md,
                                      _test2_summary.json, _symbol_resolutions.csv
genesets/msigdb_human_2026.1.Hs/      26 per-set JSON + _h.all.v2026.1.Hs.json (50 sets)
genesets/hgnc_pin/                    HGNC subset + PROVENANCE.md (md5-pinned)
code/build_genesets_human.py          builds the tiers
code/gate_disjointness_human.py       runs §11, exits 1 (the §10 Tier A fails)
results/phase7_jobA/gate_result_human.json
results/phase7_jobA/intersection_matrix_human.csv
results/phase7_jobA/intersection_matrix_BxB.csv
results/phase7_jobA/onpanel_counts_by_tier.csv
results/phase7_jobA/build_genesets_human.log
results/phase7_jobA/gate_disjointness_human.log
```

Nothing under `code/run_deepscence*`, `data/processed/deepscence_*` or `results/phase3/` was
written to. `results/phase3/summary_phase3.txt` and the DeepScence package's `coreGS_v2.csv` were
read only.
