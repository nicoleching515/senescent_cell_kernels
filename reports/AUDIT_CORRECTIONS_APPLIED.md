# Corrections applied from `AUDIT_PHASE8_FACTCHECK.md` — reporting layer, gene-set arm

**Agent:** corrections agent (gene-set / figure_gs* / pre-registration scope).
**Date:** 2026-08-27, 06:45–07:10 UTC. **Repo:** `/workspace`.
**Method:** every number below was re-derived with system `python3` from files on disk before it
was written into a report. Nothing was installed. No `git commit`, no `git push`, no tags.
`data/raw_h1/` was not read or touched.

**Moving target, recorded.** `results/phase3/a7_control_probe_{fits,curves,provenance}.csv` were
rewritten at **06:52**, mid-pass, by the M1 re-run. `a7_summary.csv` and `a7_verdict.txt` are still
at **05:19** and every A7 number quoted below was read from that 05:19 summary and re-checked at
07:05. **They will move when the re-run lands** — the per-family *ordering* (probes flat, codewords
and genomic controls not) is the finding; the decimals are provisional, as `PREREG_PHASE8.md` P2
already says. Nothing else I read changed during the pass.

**Figure guard, before and after:** `python3 code/check_figures_guard.py` →
`OK: all 27 committed figures match (PDF date stamps ignored)`, exit 0, both times.
`figure_gs2` and `figure_gs3` are untracked, so regenerating them does not disturb the
committed set — see the M7 note below on what that guard does and does not establish.

---

## 0. One-line summary

The audit is right that the science reproduces and the reporting layer does not. Two errors were
real and both are now derived from files rather than asserted: `CDKN2B` was never human-only
(**shared Tier A is 27, not 26**) and the mouse CoreScence anchor `24/35 = 69%` is not reproducible
under any convention (**it is 26/33 = 79%**). Both corrections run *against* the authors' interest
in the sense that matters — the published mouse arm was more circular than reported — while making
the *incremental* cost of C6 about half what the roadmap claimed. Nine more audit items in my
scope are applied in §3; the four the brief reserves for the M1 re-run agent — plus five neighbours
I found in the same territory — are handed over untouched, with the exact edit each needs, in §4;
and five points where I disagree with the audit are argued from files in §5.

---

## 1. R2 — `CDKN2B` is in **both** frozen Tier A sets

**Refutation confirmed, and it propagates further than the audit found.**

```
grep -c Cdkn2b genesets/A_SENDER_FINAL_strict.txt        -> 1   (mouse Tier A, promoted C6)
grep -c CDKN2B genesets/human/A_SENDER_FINAL_strict.txt  -> 1   (human Tier A)
grep -c Cdkn2b genesets/mouse_human_orthologs_MGI.csv    -> 0   (the pinned map has no row)
genesets/human/_symbol_resolutions.csv:3
  Cdkn2b,CDKN2B,NOT in MGI map; upper-case form is a known HGNC symbol
```

### What changed

| Quantity | Before | After | Settled by |
|---|---|---|---|
| Tier A shared members | **26** of 33 | **26 by the pinned map, 27 after the map-gap correction** | the two `.txt` files above + `_symbol_resolutions.csv` |
| Tier A `human_only` | `CDKN2B DDB2 GADD45G MDC1 PHLDA3 TNFRSF10B XPC` (7) | `human_only_map_gap = CDKN2B`; `human_only_real = DDB2 GADD45G MDC1 PHLDA3 TNFRSF10B XPC` (6 = 33 − 27) | same |
| Tier A `mouse_only` | `HMGA1 LMNB2 PLK3 RIF1 SESN1` (5, and 2 members unaccounted for) | `mouse_only_complete = HMGA1 LMNB2 PLK3 RIF1 SESN1 H2afx` (**6** = 33 − 27); `mouse_only_map_gap = Cdkn2b`, `mouse_only_real = H2afx` | same |
| **B7 shared members** | **85** | **85 by the pinned map, 88 after the map-gap correction** (gaps `Ccl3 Cxcl2 Cxcl3`) | `genesets/B_secondary_senescence.txt`, `genesets/human/B_secondary_senescence.txt` |
| B7 `mouse_only` | 19 mapped (9 members unaccounted for) | `mouse_only_complete` = 19 mapped + `Cxcl1` = **20** = 108 − 88 | same |
| Tier C `mouse_only` | `MMP3 TIMP1` (2, one unaccounted for) | `mouse_only_complete = MMP3 TIMP1 Cxcl1` (**3** = 14 − 11) | `genesets/C_ligands.txt`, `genesets/human/C_ligands.txt` |

`H2afx` is a **real** mouse-only member and is correctly retained as such: human `H2AFX` is in
`genesets/human/B_downstream_arrest.txt` and is removed from human Tier A by disjointness
(verified directly).

### How it is now enforced

`code/crossarm_geneset_table.py` gained `gap_split(mouse_set, human_set)`, which applies the
map-gap split the file already applied to Tier C to **every** pair of arm-matched sets, and
**asserts** the arithmetic that failed before:

```python
assert len(ho_real)  == len(human_set) - shared_corr
assert len(mo_mapped) + len(mo_real) == len(mouse_set) - shared_corr
assert len(mo_gap) == len(ho_gap)
```

so a `mouse_only` list that does not account for every non-shared member by name can no longer be
written. Re-run: `python3 code/crossarm_geneset_table.py` → exit 0. The `table` block and every
scalar in `crossarm_geneset_table.{csv,json}` are **byte-identical** to before (verified by diff);
only the `asymmetry` block changed.

### Figure

`figures/figure_gs2_crossarm_symmetry.{png,pdf}` and `_data.csv` regenerated via
`code/make_figure_genesets.py::gs2` (existing producer, `sasp_palette.apply_style`, both formats,
CSV beside). The dark inner segments now plot the map-gap-corrected counts — Tier A **27** (was
26), B7 **88** (was 85) — and the caption gives both numbers with the gap named.

### Reports updated

`PREREG_PHASE8_genesets.md` §6 (both asymmetry rows + the map-gap row now lists all six gap
symbols), §7b figure row; `BIO_PHASE8_FREEZE.md` summary item 4, §5, figure list;
`PHASE8_ROADMAP_STATUS.md` figure table.

---

## 2. M1 — the mouse CoreScence anchor `24/35 = 69%`

**Reclassified from UNVERIFIABLE to REFUTED** (see §4.2 for why), corrected, and made
reproducible.

### Re-derivation

New file `code/corescence_circularity.py` derives it from disk and writes
`results/phase7_jobA/corescence_circularity_mouse.json`:

```
CoreScence v2 occurrence >= 5 : 39 human symbols
mouse panel (D17 authoritative): 5097
on the ortholog-mapped mouse panel: strict 31, with the documented Title-case fallback 33
                                    (adds CDKN2B CXCL1)
not on the mouse panel under either convention: CXCL8 HMGB1 IGFBP5 IGFBP7 MIF TNFRSF10C

configuration  mapping   in >=1 Tier B module
pre_C6         strict    24/31 = 77.4%
pre_C6         fallback  26/33 = 78.8%
C6_promoted    strict    28/31 = 90.3%
C6_promoted    fallback  29/33 = 87.9%
```

**The convention to cite is the fallback one (33).** It is what `code/run_phase3_n8.py::corescence_mouse`
implements — i.e. what the mouse arm's Phase 3 circularity results were actually computed under —
and the script now **asserts** agreement with the project's own committed evidence:

- `git show HEAD:results/phase3/n8_disjointness_7259_….csv` → `corescence_on_panel = 33`, per-module
  overlaps `downstream_arrest 10, emt_ecm 9, il6_jak_stat3 5, interferon_response 5,
  oxidative_stress 0, secondary_senescence 14, tnfa_nfkb_proximal 8` — sum 51, union **26**.
  My re-derivation reproduces all seven cells exactly. (`logs/phase3_summary.log` §9 prints the
  same 51 and the same per-module "genes removed".)
- `logs/ds_smoke.log`: `CoreScence(occurrence>=5): 39 genes; on our ortholog-mapped panel: 31`.
- `logs/caller2.log`: `CoreScence occ>=5 : 23 up / 16 down human ; 17 up / 14 down on mouse panel`
  (= 31).

No file anywhere in the repository contains 35.

### What changed

| Where | Before | After |
|---|---|---|
| `code/make_figure_genesets.py` L245 | `bars = [("M1 mouse\npublished", 24, 35, …)` — a literal | read from `corescence_circularity_mouse.json` |
| `code/gate_disjointness_human.py` L291 | `mouse_arm_reference='24/35 = 69%'` — a literal | `mouse_arm_reference='26/33 = 79%'`, plus `mouse_arm_reference_C6='29/33 = 88%'` and `mouse_arm_source` |
| `reports/BIO_PHASE2.md` §4.2 (the origin) | 39 genes, **35** on panel; 24 of 35 = **69%** | 39 genes, **33** on panel; **26 of 33 = 79%**, with a correction box and the two added genes (`Cdkn2b`, `Cxcl1`) marked in the per-module table |
| `reports/BIO_PHASE7_JobA.md` §4 | "worse in human than in mouse"; mouse `39 / 35 / 24 / 69%` | "comparable on the two arms"; four-row table, mouse pre-C6 **26/33 = 79%** |
| `reports/BIO_PHASE8_FREEZE.md` §4, summary item 3 | 69% → 76% → 88% | mouse **79% → 88%**, human 76% → 88% |
| `reports/PREREG_PHASE8_genesets.md` §5, §7b | 24/35 = 69% | 26/33 = 79%, four-row table, correction box, and the mouse C6 row |
| `reports/PREREG_PHASE8.md` §13 | "69 % in the mouse arm" | "**79 %** in the pre-C6 mouse arm … under the promoted C6 sets the mouse arm is **also 29/33 = 88 %**" |
| `README.md` L293 | "CoreScence set is 69% circular" | "**79%** circular", with the correction pointer |
| `reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` C15 | "YES on the number" (verified `24/35` against the report) | marked superseded; verdict restated |
| `reports/PHASE8_ROADMAP_STATUS.md` | figure row "69 -> 76 -> 88%", decision D6 | corrected, with the magnitude change called out |

### Does the story change? Yes, and it should be said plainly

**It does, in two directions.**

1. **The published mouse arm is more circular than the paper has been saying** — 79%, not 69%.
   DeepScence's own gene set shares 26 of the 33 genes it can see on our panel with at least one
   Tier B response module. That is a correction against interest and it strengthens, not weakens,
   the reason DeepScence is a comparison caller rather than the primary one.
2. **The advertised "cost of C6" was roughly double the truth, and it was measured across arms.**
   The old "69% → 76% → 88%" put a mouse baseline next to two human configurations and read as a
   19-point rise. Measured within an arm it is **+9 points on mouse (79 → 88)** and **+12 on human
   (76 → 88)**, and both arms land on exactly 29/33 = 88%, with B7 alone contributing 18/33 on
   both. A smaller cost is still a real cost: C6 raises circularity on both arms, 88% is still the
   number to cite, and PI decision D6 (strip-and-refit in the frozen run order) is unaffected —
   only the magnitude quoted for it moves.

`figures/figure_gs3_corescence_circularity` now has **four** bars (mouse pre-C6 / mouse C6 / human
superseded / human frozen) so the cost is legible within each arm, and its `_data.csv` carries a
`source` column naming the producing file for every row plus both arms' per-module breakdowns.
This is a deliberate choice over the audit's alternative ("drop the 69% bar"): the derivation was
cheap and a derived bar is more defensible than a missing one.

---

## 3. Other audit items applied in this pass

### M5 — `figure_gs2`'s headline numbers were not in its `_data.csv` — **FIXED**
The dark segments and the whole footnote block came from `crossarm_geneset_table.json`; the CSV
held only `orthologue_intersected_*`, a different quantity. `figure_gs2_crossarm_symmetry_data.csv`
now carries a `dark_segment_shared_with_other_arm` column and long-format
`tier ∈ {PLOTTED, FOOTNOTE}` rows with `metric`/`value` for every caption number. The coincidence
M5 warns about is now **worse** (the corrected B7 shared count, 88, equals human B7 ∩ intersected
panel, 88) so the caption states explicitly that panel (c) is a different quantity.

### M2 — "1:1 human ortholog" and "2,425 of those land on the human panel" — **FIXED**
Re-derived: the pinned map has **18,782 rows onto 17,609 distinct human symbols, with 431 human
symbols receiving more than one mouse gene** — it is many-to-one, not 1:1. And **2,435** mouse
panel genes map onto the human panel while they land on **2,425** distinct human symbols. Both
counts are right; the sentence conflated them. `PREREG_PHASE8_genesets.md` §6 and the MGI pin row
of §7 corrected; `code/crossarm_geneset_table.py` now prints both counts and writes
`ortholog_intersected_panel_mouse_side`, `ortho_map_rows`, `ortho_map_distinct_human`,
`ortho_map_human_with_multiple_mouse` into the JSON.

### R7 — the pre-registration said the mouse C6 sets are unpromoted — **FIXED**
Verified: all **15/15** files in `genesets/mouse_c6/` are byte-identical to `genesets/`, and
`git diff --stat pre-c6-genesets -- genesets/` shows exactly the three expected changes
(25→33, 55→74, 38→108). `PREREG_PHASE8_genesets.md` §1 table and deviation **D15** now record the
promotion (PI decision D5, 05:41) and point anything needing the pre-C6 sets at the tag
`pre-c6-genesets`, which is what `crossarm_geneset_table.py` and `corescence_circularity.py` do.

### R1 repeated elsewhere — **FIXED in the two files I own**
Verified independently against `results/phase3/a7_summary.csv` (stable, 05:19): `base` design
`all_controls` −0.0697 p = 0.023; **`neg_control_probe` −0.0177 [−0.0453, +0.0099], p = 0.183**;
codewords −0.0549 p = 0.039; genomic −0.0337 p = 0.0039.
- `PHASE8_ROADMAP_STATUS.md` row 8.5b now says "pooled control features", gives the probe-only
  number, states that the pre-registered primary A7 response therefore **passes** on M1, and
  records that A7 ran on pre-C6 senders and must be re-run after 8.7.
- `PREREG_PHASE8.md` §10.1 item 1 and deviation **P2** gained the per-family breakdown and the
  same naming rule.

**Fork taken:** the audit offered "rename everywhere and amend the pre-registration" or "report
probe-only as the A7 result and pooled as a powered supplement". I took the second. Amending a
pre-registered primary null after seeing that it passes is exactly the move a reviewer will
object to; keeping `E_negative_control_probes` as the primary A7 response and reporting the
pooled families as a higher-powered supplement costs nothing and is defensible.

### M3 — "five control responses" are four overlapping views of one quantity — **partially applied**
Applied only where it touches the pre-registration and the roadmap headline:
`PREREG_PHASE8.md` **P3** and the roadmap's "A7 gave a second unplanned result" now quote
**9–13 %** with 16 % named as the `neg_probe_rate` outlier (the one response whose denominator is
an N5 column). Verified from `a7_summary.csv`: `frac_CI_excludes_zero` under `n6n5` =
0.091 / 0.103 / 0.109 / 0.127 for the clean nulls and 0.164 for `neg_probe_rate`. The instances in
`CS_PHASE8_CALLERS.md` are deferred (§4 below).

### M7 — the guard covers 27 of 45 figure artefacts — **DISCLOSED, not "fixed"**
The roadmap's figure-policy section now states the guard's scope explicitly: it enumerates via
`git ls-files figures/`, so the 18 untracked Phase 8 outputs (including all of `figure_gs1`–`gs4`)
are outside it, and `figures/.committed_manifest.json` is gitignored and re-baselineable. Closing
it properly means committing the Phase 8 figures before the freeze tag — a PI action, folded into
the M8 recommendation below.

### M10 — H1 cell counts "bracket" the mouse range — **FIXED**
Re-derived from the 11 mouse `.h5` barcode axes: **83,621–237,982** (total 1,834,806). The human
range 220k–396k sits above it and overlaps only at the top. `reports/PHASE7_H1_SCREEN.md` now says
"sit **above** the mouse range … overlap it only at the top", and notes that this is a point in
H1's favour but still requires explicit depth matching in the two-arm comparison.

### R8 — "`build_genesets.py` cannot be re-run as committed … Not fixed" — **FIXED in the roadmap**
Verified: `HEAD`'s copy has the dead `SCRATCH = /tmp/claude-0/-workspace/e6da7884…` path; the
working-tree copy (uncommitted, `M code/build_genesets.py`) defaults to
`/workspace/genesets/msigdb_mouse_2026.1.Mm` and raises `SystemExit` on a zero-JSON glob. The
"Known risks, live" row now reads **"CLOSED IN THE WORKING TREE, still open at `HEAD`"** — the
risk is real only until that edit is committed, which is not something I can do here.

### M8 — the evidence is untracked — **recommended, not done**
Confirmed in spirit: `results/phase7_jobA/`, `genesets/human/`, `genesets/mouse_c6/` including
`FROZEN_MANIFEST.csv`, and every file I corrected under `results/phase7_jobA/` are untracked. I did
**not** `git add` anything: staging files while five agents are writing the tree would create a
worse hazard than it removes, and the brief forbids committing. **Recommendation for the PI, before
`git tag phase8-frozen`:** `git add` `results/phase7_jobA/`, `results/phase3/a7_*`,
`results/phase3/*_c1*.csv`, `genesets/human/`, `genesets/mouse_c6/`, `figures/figure_gs*`,
`figures/figure2e*`, `figures/figure_phase8_*`, and `code/build_genesets.py`. A freeze over
untracked files is not a freeze.

---

## 4. Deferred to the M1 re-run agent

These fall inside `results/phase3/`, `figures/figure2*`, `CS_PHASE8_M1_RERUN.md`,
`CS_PHASE8_CALLERS.md`, `CS_PHASE7_C1.md`, `CS_PHASE8_C1_CLOSEOUT.md` or `CORRECTIONS.md`, which
that agent owns. **I did not touch any of them.** Each is stated with the exact edit it needs so
they can be applied in one pass.

### R3 — "23 % of shifted senders in the void" is the wrong column
**Where:** `CS_PHASE7_C1.md` §0 and §6.2; `CS_PHASE8_C1_CLOSEOUT.md` §4.1 footnote ‡ and the §4.2
Figure 2(e) caption; `PHASE8_ROADMAP_STATUS.md` "Known risks"/completed row ("23% void, not ~20%").
**Correction:** 23 % is `1 − frac_retaining_a_neighbour` — the fraction of shifted senders with no
real cell inside the 100 µm fitting window. The out-of-tissue measure is a different column of the
same file, `frac_in_occupancy`: **N3 35.5 %, N4 19.9 %** (against 22.8 % / 8.3 % for retention).
Either quote 35.5 % / 19.9 % and call it out-of-tissue, or keep 23 % / 8 % and restate it as
"lost every real neighbour within 100 µm". The bug is worse than the paper is about to claim, so
the fix runs in the safe direction. **Do not do both in one sentence.**
`python3 -c "import pandas as pd;d=pd.read_csv('results/phase3/null_destructiveness.csv');print(d.groupby('null')[['frac_retaining_a_neighbour','frac_in_occupancy']].median())"`
*(Note: I left the roadmap's "23% void" row alone precisely because it belongs to this item — one
owner, one pass.)*

### R4 — the figure2b/2c hash narrative in `CS_PHASE8_CALLERS.md` §4.4
**Correction:** the quoted hashes `df2f0afd…` / `983b47b2…` are
`figures/revised_candidates/figure2b_REVISED.png` and `figure2c_REVISED.png`, not `figures/`.
The committed and current files are `5ecd9ad1029851c1955fc938abf9444c` (2b) and
`a232a529a566d7c0680e04840dd07a9b` (2c). Three sub-claims need fixing: only **2d** regenerated
byte-identically to the repository (`figures/revised_candidates/README.md` says 2b/2c "Content
differs? YES"); the 05:22:41 rewrite is the documented PI-directed restore, not an unexplained
event; and the agent's regeneration was itself reverted at 05:29, so what is on disk is the
committed baseline, not what the report says it left.

### R5 — "under +N6+N5 every control family's clustered mean is indistinguishable from zero"
**Where:** `CS_PHASE8_CALLERS.md` §4.1 statement 2. **Confirmed refuted by me** against
`a7_summary.csv`: `neg_probe_rate` under `n6n5` is **+0.0108 [+0.0021, +0.0195], p = 0.0204** —
the CI excludes zero. **Correction:** "every count-based control family (4 of 5)". The report
already discloses this two paragraphs later, so only the word "every" is wrong.

### R6 — "the reportable-fit filter admits two to three times more fits than its nominal rate implies"
**Where:** `CS_PHASE8_CALLERS.md` §4.3 warning. **Correction:** 9–16 % is the *two-sided* 95 %-CI
exclusion rate under the full `n6n5` design; the reportable filter is
`beta_naive > 0 AND beta_base_lo > 0` (`code/summarize_phase3_c1.py::reportable`), one-sided and on
the **naive** design. Measured, it runs **3.0–6.7 %** on the full design (essentially nominal) and
1.2–13.3 % on the naive one. Delete the sentence or replace it with the measured 3.0–6.7 %. The
"2–3× nominal" bound on the *estimator* is fine; it does not transfer to the *filter*.

### R9 — the C1 reports no longer match `null_destructiveness.csv`
**Correction:** the diagnostic was regenerated at 05:58:41 (`logs/m1_nulldiag.log`) under the
promoted C6 sender definition, i.e. after both C1 reports were written. Section-geometry columns
are unchanged; every sender-dependent number moved (N3-orig displacement 2,974 → 2,910 µm;
N3-tile 489 → 479; N4-orig retention 0.917 → 0.920; N3-occ15 0.966/304 → 0.969/317;
N4-swap 3,038 → 2,980; 7259 admissible N3 offsets 63 → 66 of 38,080; 7352 9 → 10 of 74,178;
7259 admissible N4 angles 12 → 13 of 720; real median neighbours 140.5 → 140.0). Conclusions
unchanged; `figure2e_data.csv` carries the older set and must be regenerated with the citing text.

### Also in that territory, from my pass
- **M3, caller half.** The "0.34–0.51" and "9–16 %" ranges in `CS_PHASE8_CALLERS.md` should be
  labelled as ranges over correlated statistics, with the clean-null FPR subset given as
  9.1–12.7 %. I applied this to `PREREG_PHASE8.md` P3 and the roadmap only.
- **M4.** `CS_PHASE8_CALLERS.md` §4.1: "+0.019 SD … to +0.29 SD … a quarter of the naive
  biological gradient (+0.291 SD)" places a bin range and an amplitude in adjacent sentences with
  near-identical numerals. The span is ~0.27 SD, essentially all of the biological curve; the
  "quarter" is 0.070/0.291. Separate the two sentences.
- **M6.** `results/phase3/figure2c_data.csv` has 19 rows on disk against 11 lines at `HEAD`; the
  committed `figure2c.png` plots the 10 published nulls. Regenerate the pair together.
- **M9.** A7 (05:19) and the audited C1 (04:29) are keyed to the pre-C6 sender definition;
  `tierA_p95` moved when Tier A went 25 → 33 at 05:41. `CS_PHASE8_C1_CLOSEOUT.md` §5.1 says so for
  C1. See §5.4 below for the part of M9 I disagree with.
- **`README.md` L285–292** still carries the superseded 2-section caller-agreement numbers
  ("0.93–1.22× of chance for four of six pairs", "1.51–2.85×", "2.15× / 0.38×"). I corrected only
  the CoreScence sentence in that paragraph; the caller numbers belong to the 8.4 restatement.

---

## 5. Where I disagree with the audit

### 5.1 "complete `mouse_only` to 7 entries" (audit §6 fix 3) — **wrong, it is 6**
The audit's own R2 says the corrected shared count is **27** of 33. Then 33 − 27 = **6**
mouse-only members, not 7. The 7 comes from the *uncorrected* 26 and cannot coexist with the
correction that produced it — `Cdkn2b` cannot be both "shared, because it is a map gap" and
"mouse-only, to make the arithmetic close". Applied: `mouse_only_complete` has 6 entries
(`HMGA1 LMNB2 PLK3 RIF1 SESN1 H2afx`), asserted in code against `n_mouse − n_overlap_map_gap_corrected`.

### 5.2 M1 is REFUTED, not UNVERIFIABLE
The audit classes `24/35 = 69%` as unverifiable "because a mapping convention I have not
reconstructed might yield 35". No such convention exists, and the project's own committed files say
so: `corescence_on_panel = 33` in `git show HEAD:results/phase3/n8_disjointness_*.csv` (all nine
sections), `31` in `logs/ds_smoke.log`, `17 up / 14 down = 31` in `logs/caller2.log`. Every route
from the 39-gene CoreScence set to the 5,097-gene mouse panel yields 31 (strict MGI) or 33
(with the documented Title-case fallback); 6 genes are off the mouse panel entirely under both.
Treating it as merely unverified would have left the literal in place. It is refuted, and it is
fixed.

### 5.3 The map-gap correction cannot stop at Tier A
The audit applies R2 to Tier A only. Applying it there and nowhere else leaves the same
inconsistency one row down: frozen **B7** has three more map gaps (`Ccl3`, `Cxcl2`, `Cxcl3`) worth
**85 → 88** shared, and Tier C's `mouse_only` was incomplete by one (`Cxcl1`) in exactly the way
Tier A's was incomplete by two. `gap_split` is therefore applied to every arm-matched pair and
asserts completeness. Note that `Cxcl2`/`Cxcl5` were *already* declared map gaps by the project for
Tier C, so this is the project's own convention applied consistently, not a new one.

### 5.4 M9 is half stale
"`CS_PHASE8_C1_CLOSEOUT.md` §5.1 says so for C1; **nothing says so for A7**" is no longer true —
`PREREG_PHASE8.md` deviation **P2** already states "**Magnitudes are PROVISIONAL** — A7 was run at
05:19 on the pre-C6 sender calls; the three qualitative findings are what is frozen. A7 must be
re-run after 8.7." That file was rewritten inside the audit window (06:23–06:34, per the audit's
own moving-target log), so the auditor most likely read an earlier copy. What *was* genuinely
missing is the roadmap's row 8.5b, which marked 8.5b DONE with no re-run listed; that is now fixed.
The substance of M9 stands.

### 5.5 M5's "coincides numerically on one cell" is now two cells, by construction
After the correction, Tier A's shared count (27) equals `orthologue_intersected_mouse_C6` (27) and
B7's corrected shared count (88) equals `orthologue_intersected_human_C6` (88). These remain
different quantities. Rather than renaming columns mid-freeze I made the `figure_gs2` caption say
so in as many words and put both the pinned-map and corrected counts in the `_data.csv`. A future
pass should rename the panel-(c) columns to something that cannot be mistaken for an overlap.

---

## 6. Complete list of files changed

**Code (4 files, 1 new):**
`code/corescence_circularity.py` (**new**, derives the mouse CoreScence anchor; runnable, writes
`results/phase7_jobA/corescence_circularity_mouse.json`), `code/crossarm_geneset_table.py`
(`gap_split` + M2 wording + new JSON fields), `code/gate_disjointness_human.py` (mouse reference
derived, not asserted; guarded so a missing DeepScence install reports "NOT derivable" instead of
substituting a remembered number), `code/make_figure_genesets.py` (gs2 dark segments + data CSV,
gs3 four bars + provenance column, docstrings).

**Figures (regenerated via their existing producer, `.png` + `.pdf` + `_data.csv`):**
`figures/figure_gs2_crossarm_symmetry`, `figures/figure_gs3_corescence_circularity`.
`figure_gs1` and `figure_gs4` were **not** regenerated — nothing in them changed.

**Results (untracked, gene-set arm only):**
`results/phase7_jobA/crossarm_geneset_table.{csv,json,log}` (CSV table byte-identical; only the
asymmetry block and new scalars changed), `results/phase7_jobA/gate_result_human.json` +
`gate_disjointness_human.log` (only the three `corescence.mouse_arm_*` fields changed — the five
other outputs of that script are byte-identical, verified by `cmp`),
`results/phase7_jobA/corescence_circularity_mouse.json` (**new**).

**Reports:** `PREREG_PHASE8_genesets.md`, `PREREG_PHASE8.md`, `PHASE8_ROADMAP_STATUS.md`,
`BIO_PHASE8_FREEZE.md`, `BIO_PHASE7_JobA.md`, `BIO_PHASE2.md`,
`BIO_DELIVERABLE7_CLAIM_AUDIT.md`, `PHASE7_H1_SCREEN.md`, `README.md`, and this file.

**Not touched:** `results/phase3/`, `figures/figure2*`, `CS_PHASE8_M1_RERUN.md`,
`CS_PHASE8_CALLERS.md`, `CS_PHASE7_C1.md`, `CS_PHASE8_C1_CLOSEOUT.md`, `CORRECTIONS.md`,
`data/raw_h1/`, `SUBMISSION_PATCH_2026-08-29.md`, `genesets/*` (no gene list was edited),
`Phase7_Minimal_Human_Replication (1).md` (the PI's input spec — its "CoreScence was 69% circular
last time" line is quoted and corrected in `BIO_PHASE7_JobA.md` §4 rather than rewritten).

**Reproduce everything in this file:**
```
python3 code/corescence_circularity.py
python3 code/crossarm_geneset_table.py
python3 code/gate_disjointness_human.py        # exit 0
python3 -c "import importlib.util as u;s=u.spec_from_file_location('m','/workspace/code/make_figure_genesets.py');m=u.module_from_spec(s);s.loader.exec_module(m);m.gs2();m.gs3()"
python3 code/check_figures_guard.py            # OK: all 27 committed figures match
```
