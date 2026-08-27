# Recording four PI decisions in the frozen pre-registration

**Date:** 2026-08-27. **Author:** documentation agent, under the PI's four decisions.
**File edited:** `reports/PREREG_PHASE8.md` — **and nothing else.** No file under `results/`,
`figures/`, `code/`, `genesets/`, `data/` or any other report was written or modified; `data/raw_h1/`
was not read; nothing was written under `results/phase9_h1/` or `results/phase10*`.
**Freeze:** `phase8-frozen`, `d04691e2692a7be8d1ff676d2fb74ad9d1df049d` — see item 9 of §3 below,
which is not the hash §1 of the pre-registration records.
**Gate:** `python3 code/gate_genesets_guard.py` → **exit 0**, mouse PASS + human PASS, "no drift:
96 watched gene-set/panel files match the manifest".

---

## 1. What was recorded, and where

A new **§0.1 PI DECISION BLOCK — dated 2026-08-27**, immediately after the §0.0 correction block,
in the same form §0.0 uses: a four-row summary table, then one subsection per decision, then a
list of open items, then a closing statement of what the decisions do and do not change. The
`D-A … D-D` namespace is declared new and disambiguated from the four `D*`/`P*`/`H*` series
already in play.

**No silent rewrite, and no deletion.** Every affected site keeps its original wording and carries
a dated `***[…]***` marker beside it. Verified mechanically after the edit: of the six existing
lines the diff touches, all six still contain their original text verbatim (`git diff -U0` old
lines re-found as substrings of the new file). The two markers that would have split a table
(§3.7 and §3.9) were moved below their table, with a one-line pointer left in the row.

| Decision | Recorded in | Sites marked inline |
|---|---|---|
| **D-A** DeepScence on H1 → five-seed consensus + dispersion; M1 keeps its single-seed score | §0.1 D-A | §3.9 (`random_state` row + marker under the table), §7 outcome row **F**, §9 item 4, §10.10 |
| **D-B** Fine and merged both computed; **merged primary on H1**, fine = frozen-literal sensitivity | §0.1 D-B | §3.7 (`tierA_pNN` row + marker under the table), §5, §6 (before R2, covering R1/R2), §13 head-note |
| **D-C** **P-ii falsified**, reported at the front of §8, with what it does *not* touch | §0.1 D-C | §8 preamble (above the prediction table), §8 "Either answer publishes", §13 **P12**, §13 **P13** |
| **D-D** **P7 retired as UNANSWERABLE**, original text preserved | §0.1 D-D | §10.4, §13 **P7** |
| — | §14 | closing status note: two of Phase 9's three requested PI decisions are now taken; the SenePy one (item 4) is **not** |

---

## 2. Every number, and the file behind it

Nothing below is quoted from a document. Each row was re-read from the named file by the command
recorded beside it in §0.1.

| Number as recorded | File | Command / key |
|---|---|---|
| H1 seed check: *r* = **0.3719**, Spearman 0.4157, top-5 % Jaccard **0.2107**, **12,779** cells changing, n = **196,142** | `results/phase9_h1/d2_stability.csv` | row `SPLN21, "denoise=False, FULL section", 0, 1` |
| H1 20,000-cell panel companion: 0.3829 / 0.2006 / 1,327 | " | row `SPLN21, denoise=False, 0, 1` |
| M1 seed-to-seed floor: *r* = **0.99553**, Jaccard **0.7606**, 272 changing | `results/phase8_d2/d2_stability.csv` | row `raw_seed0 vs raw_seed1` |
| M1 determinism control: *r* = **0.99999913** (7239, 75,384 cells, 24 calls moving) and **0.99999995** (7259, 114,721, 2) | `results/phase8_d2/d2_agreement.csv` | rows `config == "raw"` |
| H1 depth loading survives the seed: ρ 0.3122 (seed 0) → 0.2308 (seed 1) | `results/phase9_h1/d2_depth.csv` | row `SPLN21, "full section, seed check"`, cols `rho_denoise_False`, `rho_seed1` |
| The five seeds **20260901–20260905** | `code/run_phase8_compmatch.py` | line 165, `COMPMATCH_SEEDS`; already frozen at §3.8 |
| Eligible-but-uncallable per section: **0.00 / 13.64 / 10.28 / 1.88 / 0.00 / 9.17 / 9.20 %** (counts 0 / 38,701 / 20,155 / 7,386 / 0 / 24,815 / 27,588) | `data/processed_h1/celltypes_h1_*.csv` × `data/processed_h1/annotation_meta_h1_*.json` | fine label in {Unknown, unknown, Low_quality} AND merged label not; denominator `n_cells_qc` — full one-liner in §0.1 D-B |
| T/NK at `tierA_p95` vs `tierA_merged_p95`: SPLN14 23 → 1,958; **SPLN43 0 in 24,815** → 1,241; SPLN44 295 → 1,675; prevalence 0.059 / 0.000 / 0.881 % → ≈ 5.00 % | `results/phase9_h1/a3_prevalence_by_type.csv` | `label_set == "merged"`, `cell_type == "T/NK cells"` |
| `CDKN1A` fold-split sign stability **1.00 in all 7 sections** | `results/phase9_h1/deepscence_anchor_h1.csv` | col `stab_cdkn1a` |
| **20** folds | `code/h1_deepscence_anchor.py` | line 103, `def fold_stability(v, k=20, seed=0)` |
| depth-partialled ρ(`CDKN1A`, score) **+0.1911 … +0.2540**, positive in 7/7; `CDKN1A` rank **1st–7th** of 33 | `results/phase9_h1/deepscence_anchor_h1.csv` | cols `rho_partial_cdkn1a`, `rank_cdkn1a_in_core` |
| D3 8-gene proliferation anchor for contrast: ρ **+0.0097 … +0.0451**, stability **0.75** in SPLN21 | " | cols `rho_partial_prolif`, `stab_prolif` |
| `LMNB1` +0.2029 … +0.2449; consensus anchor +0.0642 … +0.1695 | " | cols `rho_partial_lmnb1`, `rho_partial_consensus` |
| Native CoreScence circularity **29 / 33 = 0.8788 (88 %)** | `results/phase7_jobA/gate_result_human.json` | `corescence.frozen_n_in_any_B`, `.n_on_panel`, `.frozen_frac`; independently re-printed by the gate guard run at the end of this task |
| DeepScence × `CDKN1A`⁺ pooled **6.436**, z **204.83**, above chance **7/7**, `circular = True` | `results/phase9_h1/caller_agreement_pooled.csv` | row `deepscence_score, cdkn1a_counts` |
| its per-section range **3.459 – 10.115** | `results/phase9_h1/caller_agreement_matched_significance.csv` | same pair, col `ratio_stratified` |
| P-iii pooled **1.102**, z 5.67 (against its own 1.10 threshold) | `results/phase9_h1/caller_agreement_pooled.csv` | row `tierA_score, abs_deepscence_score` |
| Tier A × SenePy **0.874**, z **−7.96**, below chance 7/7 | " | row `tierA_score, senepy_score` |
| Tier A × DeepScence SPLN21 seed-0 ratio **1.174**, z 3.54 | `results/phase9_h1/caller_agreement_matched_significance.csv` | `SPLN21, tierA_score, deepscence_score` |
| `Marginal zone B cells` **never realised in any of the 7 sections** | `data/processed_h1/celltypes_h1_*.csv` | union of `cell_type` over the seven files → 19 labels, MZ absent; label *is* in `code/markers_human_spleen.py:27` |
| P7's "**6** CellMarker spleen rows" — verified, not taken on trust | `genesets/human/markers_spleen_evidence.csv` | 8 marker rows for the label, each stamped `ALL-TISSUE FALLBACK >=1 PMID -- WEAKEST EVIDENCE (6 spleen rows)`; same tier in `results/phase7_jobA/build_markers_human_spleen.log` line 18 |
| SenePy within-type depth enrichment **Q5/Q1 = 28.5 – 228.1×** | `results/phase9_h1/caller_within_type_depth_bias.csv` | recomputed; see §3 item 8 — the audit quotes 28.5–224.7 from unrounded values |
| A6 axis composition (used to flag P7's stated rationale) | `code/h1_a6_compartments.py:85-86` | `pulp = score(D_spleen_red_pulp) − mean(follicle, tzone)`; `D_spleen_marginal_zone` is scored but is **not** a term in the axis |

---

## 3. What these decisions contradict, and what is still open

All nine are recorded **inside** §0.1 as well, so they cannot be lost by anyone reading only the
pre-registration. Ranked by how much damage an unnoticed one would do.

1. **D-B contradicts the literal text of §5, §6 R1 and §6 R2, which name `tierA_p95`.** This is the
   most consequential of the four. R1 is the pre-registered **primary replication criterion** and
   it names the call. Under D-B the H1 primary call is `tierA_merged_p95`. The threshold, the
   direction and M1's own outcome are untouched — only which H1 fit population R1 is scored on
   moves. **The PI has not said whether R1 is scored on the merged call alone, or on both with the
   fine call as the declared sensitivity.** Marked at §5 and §6; **unresolved.**
2. **D-B leaves the N7 sender axis unspecified.** `N7_CALLS` (`run_phase3_nulls.py:68-69`) is a
   frozen six-call list, nine with `TIERA_PM_CALLS`, and `tierA_merged_p95` is in neither.
   Replace, or add? **Unresolved**, and it is a frozen-list question, not an implementation detail.
3. **D-A and D-C interact and the interaction is unresolved.** Every H1 DeepScence quantity in
   `CS_PHASE9_H1_AUDIT.md` §10 — including D-C's own 6.436 circularity, the P-ii stabilities, P-i's
   depth loadings, P-iii's 1.102 — is computed at `random_state = 0`. The moment the H1 primary is
   the five-seed consensus, those become numbers on a **non-primary** estimator. Directions are
   seed-robust and the P-ii falsification is a *direction* (sign stability 1.00 in 7/7, and a
   consensus can only be more stable than its members, not less), so **D-C's verdict is safe**;
   its magnitudes are not.
4. **P-iii is confirmed by 0.002 and is seed-fragile.** §8 registers "pooled ratio > 1.10"; H1
   returns 1.102. On an estimator whose full-section seed-to-seed Jaccard is 0.211, that is not a
   margin. **P-iii must be re-scored on the consensus before it is reported either way** — and it
   could go to either verdict.
5. **D-A does not define the consensus.** The pooling rule (per-cell mean of z-scored per-seed
   scores? rank-mean? median?) and the dispersion statistic reported beside it are **not**
   specified by the decision, and I did not invent them. They must be written down **before** the
   consensus is computed, or the estimator is chosen after seeing the data.
6. **D-A and D-B each create a cross-arm estimator asymmetry.** M1: single-seed DeepScence,
   fine-label sender call. H1: five-seed consensus, merged-label sender call. §9 item 4 requires
   every cross-arm number twice **on panel grounds**; these are two further axes on which the arms
   are not the same measurement. Whether M1 is re-called at the merged family for comparability is
   **unresolved**. Marked at §9 item 4.
7. **P7's stated reason for keeping the MZ label was already wrong.** P7 argues that dropping the
   label "would leave the A6 axis without its middle term". The H1 axis is
   `score(D_spleen_red_pulp) − mean(score(follicle), score(tzone))`
   (`code/h1_a6_compartments.py:85-86`): `D_spleen_marginal_zone` is scored as one of five
   compartments but **is not a term in the axis**. Retiring P7 costs the axis nothing. **Flagged,
   not corrected** — it is P7's own wording and D-D preserves it deliberately.
8. **Two numbers in the evidence chain do not reproduce from a file.**
   (a) `CS_PHASE9_H1_AUDIT.md` §10.5 quotes the Tier A × DeepScence matched ratio moving
   **1.174 → 2.967** between seeds on SPLN21. The seed-0 value is in
   `caller_agreement_matched_significance.csv`; **the seed-1 value is in no file under
   `results/`**, so §0.1 does not rely on it and says so. It should be regenerated into a file or
   dropped from the audit. (b) §10.2's SenePy Q5/Q1 range **28.5–224.7** does not reproduce from
   the committed `caller_within_type_depth_bias.csv`, which is rounded to 3 dp and gives
   **28.5–228.1** (SPLN24 4.105/0.018). §14's new note quotes the file-derived range and states the
   discrepancy. Neither changes any conclusion; both are the exact failure mode the project has
   twice caught in itself.
9. **§1's frozen tag hash no longer matches the tag** — recorded as an observation, **no decision
   taken, §1 left unedited.** §1 records `926439629a07269a32c93f998da0f6e1cd20933c` at
   2026-08-27 15:32 UTC. On disk `phase8-frozen` is an annotated tag dated 2026-08-27 19:55:58 UTC,
   subject "Phase 8 frozen (**RE-CUT**): pipeline and pre-registration for the human replication",
   resolving to **`d04691e2692a7be8d1ff676d2fb74ad9d1df049d`**; 9264396 is an ancestor of it.
   `CS_PHASE9_H1_AUDIT.md` also cites 9264396 as the freeze. The re-cut is outside these four
   decisions and needs the PI's own record of why the tag moved and which commit the freeze is.
10. **P-vi is untouched and must keep its awkward description.** Contradicted on the one section
   run, but its own falsifier needs 5 of 7 and Phase 9 ran 1 — "contradicted, not falsified".
   None of D-A–D-D changes that.
11. **§14 item 4 (SenePy) is still open.** Phase 9 asked for three PI decisions; D-A and D-B take
   two. SenePy now carries four caveats on H1 (no spleen hub, the 100 µm truncation, Q5/Q1
   28.5–228.1×, and the v1/v2 hub-release discrepancy). **No decision is recorded for it here**,
   and §14's new note says so explicitly rather than letting the silence read as a decision.

---

## 4. What was deliberately **not** done

- **No number was changed anywhere in the frozen text**, and no C-1…C-11 correction was touched.
- **§1's tag hash was not corrected** (item 9) — outside the four decisions.
- **`genesets/`, `code/`, `results/`, `figures/`, `data/` were not modified.** The frozen marker
  file `code/markers_human_spleen.py` still carries `Marginal zone B cells`, as P7's retirement is
  a claim decision, not a pipeline change; the label set stays at 23.
- **The consensus estimator was not implemented or specified** beyond what the PI decided (item 5).
- **P-iii, P-i, P-iv and P-v verdicts were not restated** as if they were computed on the new
  primary estimator (item 3).
