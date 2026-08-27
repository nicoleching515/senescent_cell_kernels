# MOUSE arm under the adopted C6 decision — NOT YET PROMOTED

Built 2026-08-27 by `code/build_genesets_mouse_c6.py`. Panel: GSE310392 Xenium Prime Mouse 5K
+ 100 custom = 5,106 `Gene Expression` features, of which 9 are genotyping probes → **5,097 genes**
(asserted in the script, not assumed).

**Why this directory exists.** The PI adopted the re-sourced B7 for the human arm. Applying it to
one arm only would make Phase 7 §17 compare a 116-gene human B7 against a 38-gene mouse B7. These
are the mouse counterparts, built by the *same* method from the *same* pinned sources.

**Why it is not in `genesets/`.** Every Phase 2–5 mouse result on disk was computed under the
**pre-C6** sets in `genesets/*.txt`. Overwriting them would silently invalidate published numbers.
Promoting this directory means **re-fitting the mouse arm**, which is a PI decision, not a file
copy.

| | pre-C6 (`genesets/`) | C6 (`genesets/mouse_c6/`) |
|---|---|---|
| `B_secondary_senescence` | 38 (curated, Neretti meeting abstract) | **108** (SAUL_SEN_MAYO ∪ REACTOME_SASP, minus the Tier A caller pool) |
| `A_SENDER_FINAL_strict` | 25 | **33** |
| B1–B6 | 126 / 68 / 100 / 190 / 125 / 31 | unchanged |
| per-module senders | 70 / 74 / 73 / 37 / 73 / 71 / 55 | 70 / 74 / 73 / 37 / 73 / 71 / **74** |

Source sets are MSigDB's **own mouse** versions from the existing pin
(`../msigdb_mouse_2026.1.Mm/`, release 2026.1.Mm, fetched 2026-08-20): `SAUL_SEN_MAYO` = MM16098
(117 genes, 90 on-panel) and `REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP` = MM14900
(40 genes, 23 on-panel). **No ortholog mapping was used for B7** — the same discipline the mouse
arm used for B1–B6. Nothing was re-downloaded.

**Self-check.** The script reproduces the mouse arm's published Tier A/B sizes exactly
(B1 126, B2 68, B3 100, B4 190, B5 125, B6 31, B7 v1 38, A0 74) before applying C6, and asserts
it. If that assertion ever fails, the rebuild is wrong and the script stops.

**§11 gate on the mouse arm, C6 configuration:** `A_ported(A0)` **PASS** (74 → 33 after removal);
per-module gate **PASS** (all seven ≥15 and disjoint from their own readout); the §10 sixteen in
mouse symbols **FAILS** exactly as in human (13/16 on panel, 5 survive).

`variants/` holds the superseded curated B7 (38), the sourced union before caller subtraction
(113), and the §10 sixteen in mouse symbols — all reported, none used.
