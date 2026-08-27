# Plan update D12/D13 — §22 N2, §25 Fig. 1, §29 and §30 brought onto the settled numbers

**Date:** 2026-08-27. **Scope:** `SASP_Kernel_Master_Plan.md` only, plus this file.
Nothing in `results/`, `figures/`, `code/`, `genesets/`, `data/` or any other report
was written. `data/raw_h1/` was not opened. No commit, no tag.

**PI decision this discharges:** D12/D13 (PI-decision series) — *"hold paper-level
edits until the runs land; §29 objections 3 and 6 and the §30 outline stay stale by
choice. Re-raise once 8.7 and 8.5 report."* Tasks 8.5, 8.5b and 8.7 have all
reported, so the hold is lifted and the edits are made.

**Provenance rule applied throughout.** Every number written into the plan was read
from a file in `results/` or `figures/`, not transcribed from a report's prose. The
caller-agreement decomposition (1.030 → 1.128 → 1.212) was **independently
re-derived here** by Mantel–Haenszel pooling straight from
`caller_agreement_matched_significance_{verify2sec,2sec_c6,11sections}.csv`; it
reproduces `caller_coverage_gate_headline.csv` and
`SUBMISSION_PATCH_2026-08-29.md` §3.1 exactly (3-pair: 1.0299 p=0.203 → 1.1283
p=4.40e-8 → 1.2122 p=1.8e-94; 4-pair: 1.0398 → 1.1312 p=6.5e-9 → 1.2120 p=1.8e-106).

---

## 0. LEAD — still factually wrong, and I could not fix it

Ordered by how much damage each does if the PI ships without seeing it.

**0.1 §29 objection 9 promises a Moran's I that does not exist.** The response says
"we … report our **own** Moran's I on the controls alongside the kernel amplitude".
`grep -ril moran` over `code/` and `results/` returns **nothing**. The objection is
currently answered by assertion. This is the novelty review's **highest-priority
open gap** (`NOVELTY_ASSESSMENT.md` §4 O1) and it is the designated defence against
the one prior-art hit that directly falsifies a project claim (Voyager's Xenium
vignette). I marked the gap inside objection 9 with a ⚠ STATUS line rather than
quietly leaving the promise standing, but **fixing it requires a run, not an edit.**
It is cheap — the control features are already in the M1 bundles.

**0.2 §30 5.8 now describes a two-arm replication that has not happened.** Figures 5
and 6 and the human arm are Phase 10 work, and Phase 9/10 cannot begin until the
`phase8-frozen` tag and pre-registration are committed — which is item 8.9, owned by
the **PI**, and is the only thing blocking it (`PHASE8_ROADMAP_STATUS.md`). H1
(GSE326743, 7 human spleens) is acquired and structurally verified; no H1 estimate
exists. I wrote the subsection with an explicit status caveat — *"if it has not run
by submission, say so and report the mouse arm alone"* — but an outline that lists
work not yet done is a live hazard at a 4- or 8-page limit.

**0.3 `Phase7_Minimal_Human_Replication (1).md` §17 still calls H1 "human aging
lung".** H1 is **spleen**. The same table also carries pre-C6 values throughout
(0.326 / 0.027 / 0.203 / SF-N2 0.943 / "composition surrogate share 66–76 %" /
"0.93–1.22× chance"), every one of which is superseded (0.329 / 0.029 / 0.183 /
0.9516 / 66–85 % / 0.751–2.198×). **That file is outside my edit scope**, so it is
reported, not fixed. It is the file §30 5.8 points at for the Figure 5/6
definitions, so a reader following the pointer lands on stale numbers and the wrong
organ.

**0.4 §31 was not extended, deliberately.** §30's Related Work now names four
sources the paper needs — ICE (*Genome Biology* 2026), markeR (*NAR GAB* 2026),
bioRxiv 2026.01.02.697374 (the "circular validation" preprint), and stAge (bioRxiv
2025.11.23.689860). I verified that **none of the four is in `references.bib`**
(43 entries) or in §31. I did **not** add §31 entries for them: I have no repo file
carrying their bibliographic records, and this repo has already been burned once by
a bibliography written from recall (19 of 32 entries, 41 wrong forenames,
`CITATION_AUDIT.md`). They need to be retrieved and entered into `references.bib`
with the same `% SUPPORTS:` discipline as entries 31–38, then mirrored into §31.

**0.5 The §29 "Venues" paragraph still records "Deadline August 29 AoE" for
ml4spatialbio.** The CfP now shows a **window, Aug 29 – Sept 4 AoE**
(`NOVELTY_ASSESSMENT.md` §5). I recorded that in the new reassessment subsection and
**left the Venues paragraph untouched**, because the brief reserves the venue
decision — including its stated deadline — to the PI. If the Sept 4 window holds it
changes the sequencing option materially, so it should not stay buried in a
subsection.

**0.6 Two smaller staleness items I did not touch.** §1 still says *"a 4–9 page
workshop paper with **4 main figures**"* — it is six (Figs. 5 and 6 are defined in
the addendum §19), and the page target now depends on a venue decision that is open.
§25's Figure 2 caption still describes panel (b) as carrying "the torus-shift null
band"; the primary N3 is now **N3-var**, not the bounding-box torus that caption was
written against, and Figure 2 has grown panels e–h. Both are one-line fixes I left
for whoever owns §25/§1 so as not to widen the edit surface further.

---

## 1. What changed, section by section

### 1.1 §22 Step 3, null N2 — Task 2 ("the single most important number")

**Struck:** *"N2 — Matched-decoy senders (**the critical control**) … Report
β̂_true − β̂_decoy as the corrected effect. **This is the single most important number
in the paper.**"*

**Written:** the decoy contrast is *still central, as a negative result about the
method, not as the primary estimate*; report it **beside**, never instead of, the
covariate-adjusted estimate. Backed by a three-source correction block:

| evidence | matched decoys remove | same variables as covariates remove | file |
|---|---|---|---|
| composition-matched protocol, 5 seeds | **1.6 %** (SF 0.9837 [0.973, 0.994]) | **85.4 %** (SF 0.1461 [0.052, 0.246]) | `results/phase3/compmatch_reruns.csv` |
| A7 negative-control features | ~0 % (N2 −0.0642, p = 0.0124 vs naive −0.0744, p = 0.0145) | all of it (N5 +0.0038, p = 0.72) | `results/phase3/a7_summary.csv` |
| Figure 1 synthetic grid | see §2.3 below | — | `figures/figure1_data.csv` |

Plus the named mechanism (`CS_PHASE3.md` §5 verbatim: matching balances covariates
between senders and decoys but not the response's dependence on them at the
**receiver**), and the two pre-registered consequences: **no naive or N2-only kernel
may be reported**, and `type_adj` (65.9 %) / `typecomp_adj` (85.4 %) travel wherever
0.9837 does (`PREREG_PHASE8.md` §10 rule 8).

*Note on section numbering: the brief calls this "§23's sentence". The sentence is
in **§22 Step 3** (the null battery); §23 is Baselines. The sentence edited is the
one quoted in the brief and in `NOVELTY_ASSESSMENT.md` §2.1.*

### 1.2 §29 objection 3 — "your senescence calls are inferred"

**Struck:** *"We report sensitivity across thresholds and across **two independent
signature sources** and show the conclusions do not depend on the choice. If they
do, we say so."* The first claim is false and the promise is now due.

**Written:** four things stated rather than left to a reviewer —
(i) the callers are **not independent**: 0.751–2.198× of chance over eleven
sections, pooled **1.212×** (z = 21.9, p = 1.8e-106), and on the *published two
sections* the sender-set fix alone moves 1.030 (p = 0.20) → **1.128 (p = 4.4e-8)**;
(ii) the dependence spans 0.737× to 1.471× with a different identifiable cause at
each end, and the older "direction predicted by depth loading" plank is **refuted**
(pair-level permutation p = 0.30; within-pair Spearman negative in 5 of 5) and must
not ship; (iii) DeepScence's polarity is a `CDKN1A` anchor that is weak or reversed
in **four of eleven** sections once depth is partialled out (ρ −0.024…+0.182), to be
framed as documented method behaviour, not as a defect we found; (iv) SenePy ships
**no spleen signature** (0 matched / 15 surrogate / 7 none of 22 cell types), so the
second source does not exist on the human arm.

Then the sentence that survives: the **bound is invariant to the sender
definition** — N3-tile **0.960–1.001** across seven sender definitions spanning
0.50–8.96 % prevalence — *"the bound is invariant to the caller; our
caller-independence claim was not, and we withdraw it."*

Sources: `caller_coverage_gate_headline.csv`,
`caller_agreement_matched_significance_*.csv`, `deepscence_anchor_decisions.csv`,
`CS_PHASE8_CALLERS.md` §5.3, `senepy_spleen_coverage.csv`, `sf_summary_c1_n7.csv`,
`m1_n7_prepost.txt`.

### 1.3 §29 objection 6 — "isn't the decoy control too conservative?"

**The answer inverts.** Struck: *"Possibly, and we quantify it … That number is in
Figure 1c."* Written: **no — on this assay it is insufficient**, with the A7 and
composition-matched evidence, the Figure 1a/1b synthetic arm, and the pre-registered
consequence. Two qualifications are written in because they are the two ways this
result gets misquoted: the −0.0744 is on the **pooled** control set and never on the
40 named probes (which are flat, −0.0225, p = 0.129, and are the pre-registered
*primary* A7 response); and `neg_probe_rate` is flat naively (+0.0113, p = 0.232),
so the gradient is **detection efficiency**, which is exactly why an N2 match cannot
see it.

### 1.4 §29 objections 1, 4, 7, 9 — the other four, checked as instructed

- **1 (Ma et al.).** The attribution was already right in §29 — the repo corrected
  Zhao → **Ma S** (Zhao L is 14th of 47) and §29 carried the correction note. What
  was stale was the *answer*: "we quantify how much survives **matched-decoy and
  torus-shift nulls**" names the two nulls that are no longer the primaries.
  Rewritten to N5/N6 + **N3-var**, with the matched-decoy contrast demoted to
  diagnostic, and the actual answer attached (SF **0.088**, IQR across fits
  [−0.017, 0.234]; controlled **0.029**, IQR across fits [−0.007, 0.084], vs a
  **0.183** bound **over 0–100 µm**). *(Bracket type labelled 2026-08-27: both are
  IQRs across the 153 reportable fits, not confidence intervals —
  `results/phase3/sf_summary.csv` carries `q25/median/q75` and no CI column.)*
  Verdict on "which is right": **Ma et al.** — `references.bib` keys it
  `ma2024spatial`, PMID 39500323, and the novelty review verified it independently.
- **4 (only one or two tissues).** Made concrete and honest: one mouse tissue,
  eleven liver sections, six §8-Test-3-admissible; one human tissue, seven donor
  spleens, **acquired but not yet run**; species and tissue confounded by design;
  D4 (age continuous, no stratified claim).
- **7 (torus / Mrkvička).** Was *"we either run the variance-corrected shift … or
  justify tiling"*. **We ran it.** N3-var **0.996** [0.975, 1.007] vs N3-tile 0.971
  and published bbox 0.999; N4-var 0.985; window-matched N3-var 0.995. Plus the
  calibration measurement: tiled torus **0.080–0.118 against nominal 0.05 (up to
  2.35×)** while RS_count holds **0.033–0.060** — i.e. reported against interest,
  **C1 replaced a liberal test with a more liberal one**, which is why N3-var is
  primary. Both honest caveats carried (synthetic study; on real data the tile looks
  slightly *more* conservative because seams sit **~81 λ̂** apart — at the sourced
  pooled λ̂ of **14.7 µm**, not the withdrawn 15.7 µm; corrected 2026-08-27).
  `sf_summary_var.csv`, `var_sim_calibration.csv`, `CS_PHASE8_TORUS_VAR.md` §1–§4.
- **9 (Voyager).** Content unchanged — it was already correct — with a ⚠ STATUS flag
  added: see §0.1. The quotable A7 numbers were attached so the objection is not
  answered on an unrun computation alone.
- **2, 5, 8** were checked and left **verbatim**. Objection 2 is the one the novelty
  review explicitly says to keep; 5 and 8 rest on nothing that moved.

### 1.5 §29 Venues — Task 4, flagged subsection added, decision untouched

New `#### ⚠ Venue reassessment — evidence for the PI, no decision taken`, inserted
after the Secondary paragraph and before the objections. It carries the retrieved
facts as a table (4 pages non-archival with concurrent submission permitted and a
window to ~Sept 4, vs 8 pages and an Aug 29 AoE deadline with a topic list that
matches line by line), the argument for inverting, the argument for keeping
ml4spatialbio as a **different, shorter paper** (the negative-control-probe
calibration instrument), and the two unresolved checks — **ICBINB-BIO's
dual-submission policy was not verified**, and the Sept 4 window has to hold. The
Primary/Secondary paragraphs above it are **byte-unchanged**.

### 1.6 §30 — Task 3, outline extended

- Header: page target now points at the venue subsection (8 pp = ICBINB full track;
  4 pp fits only a subset of §5).
- **Introduction**: must set up the two claims the paper is actually strongest on
  (N2 is not conservative; shift nulls on windows where they are undefined) and must
  state **in the introduction** that no published length constant with uncertainty
  exists — the contribution is a bound plus an identifiability argument, not the
  overturning of a number.
- **Related Work**: four missing citations named (ICE, markeR, bioRxiv
  2026.01.02.697374, stAge) and the instruction not to claim novelty for "callers
  disagree"; DeepScence's `CDKN1A` anchoring to be cited as documented behaviour.
- **Results**: 5.1–5.4 kept, with 5.1 told to carry panels a/b (not only c) and 5.3
  told to state the 0–100 µm range in the same sentence as the bound. **New 5.5**
  negative-control-outcome kernel + measured FPR with its three caveats; **new 5.6**
  the N2-vs-N5 result, given its own subsection because an external review ranks it
  #1 of six; **new 5.7** shift nulls, framed as import-and-quantify with the three
  genuinely-ours items (35.5 %/19.9 % out-of-tissue — *not* to be merged with the
  22.8 %/8.0 % retention column; the 2.4×-nominal calibration; the FFT enumeration,
  which leads the methods paragraph); **new 5.8** two-arm replication and Figs. 5–6
  with the status caveat; **new 5.9** what was withdrawn and why that is a result.
- **New section 6, "Corrections, pre-registration, and reproducibility"** — the
  pre-registration, the `phase8-frozen` tag and `CORRECTIONS.md`, and the fact that
  every correction moved against interest and none changed the conclusion. Discussion
  and Appendix renumbered 7 and 8.
- **Discussion**: the bound as a *constraint* on Martin et al.; the no-prior-number
  point stated as both strength and risk, and hence why Figures 1 and 4 are
  load-bearing; limitations updated (non-independent callers, `denoise=False` as
  frozen primary).
- **Appendix**: adds the eleven-variant shift-null family, the five-seed composition
  protocol, the per-section caller tables, the deviation tables — and the warning
  that the **seven Tier B modules are not mutually disjoint** (18 of 21 pairs share
  genes; only A ∩ ∪B = 0 is gated), so no text may claim otherwise.

---

## 2. Found and not in the brief

### 2.1 The plan still said "Zhao et al." in two places — fixed

§0.1 item 5 ("Read Zhao et al.") and the §2.1 "Taken" table. Everywhere else in the
file the 2026-08-21 D7 correction to **Ma S** had been applied, including §29
objection 1, so these were the last two. Fixed with the plan's own bracketed
correction convention. This is the one place the brief's Task 1 wording ("make §29
consistent") turned out to be satisfied already, with the inconsistency living
*outside* §29.

### 2.2 §25's Figure 1 caption is contradicted by Figure 1's own committed data — fixed

The caption read *"the decoy control restores approximate calibration."* It does not.
From `figures/figure1_data.csv` (the committed data behind the committed figure):

| | naive | matched decoy (N2) | nuisance-conditioned (N5) |
|---|---|---|---|
| \|relative bias\| in λ̂ smaller than naive | — | **12 of 20 cells** (larger in 8) | — |
| worst-case \|relative bias\| | 2.03 | **2.27** | **0.33** |
| mean CI coverage where ℓ ≥ 2λ (8 cells) | 0.51 | **0.35** | **0.85** |

**This is a third, independent line of evidence for the N2-vs-N5 inversion, and it
was sitting in the project's own headline synthetic figure the whole time.** It is
stronger than the two measured lines in one specific way: it is on **planted ground
truth**, so it is not open to the "your covariates removed real signal" reading. It
is now cited in §22 N2, §29 objection 6 and §30 5.6, and the caption carries a
correction note. `make_figure1.py` confirms `decoyS` is the matched-decoy (N2)
design and `nuis` the covariate design, so the column semantics are not assumed.
**Figure 1 itself needs no regeneration** — panels a and b already plot this; only
the prose was wrong.

### 2.3 Three claims in the brief that I verified rather than transcribed

- "1.128× at p = 4.4e-8 on the published two-section base **once the sender set is
  fixed**" — reproduced exactly by independent Mantel–Haenszel pooling (§0 above).
  Note the brief's 1.212× and this 1.128× are on the **three-pair** basis; the
  gate file's headline rows are the **four-pair** basis (1.212 / 1.131). Both are in
  the plan text with the basis named, because mixing them is exactly the bug
  `caller_coverage_gate.csv` was fixed for (ledger row 34).
- "weak or reversed in four of eleven sections" — confirmed against
  `CS_PHASE8_CALLERS.md` §5.3 item 4, which names the four (7248, 7435, 7352, 7001)
  and gives the depth-partialled ρ range −0.024…+0.182.
- "up to 2.4× nominal" — confirmed from `var_sim_calibration.csv` directly:
  irregular window, s = 0.30, 4×4 tiling **0.1175** against nominal 0.05.

### 2.4 Numbers in the brief that are one revision behind, and what I wrote instead

The brief quotes A7 as **−0.0744, p = 0.0145** and **N2 −0.0642, p = 0.0124** — both
current, both used. It also quotes "pooled 1.212× post-C6", current. No conflict
found. Where the brief and a report disagreed with a CSV, the CSV won; that happened
once, on the FPR upper end (`CS_PHASE8_CALLERS.md` quotes 9–16 % from the pre-C6
run; `a7_summary.csv`'s post-C6 `n6n5` column gives 0.091 / 0.103 / 0.109 / 0.145 /
0.164, so 9–16 % is still right and 9–15 % is right for the four count-based
responses — which is how it is written).

---

## 3. What was *not* changed, and why

- **The venue decision.** Reserved to the PI. Evidence added, assignment untouched.
- **§29 objections 2, 5, 8.** Nothing they rest on moved.
- **§25 Figures 2/3/4 captions, §1's "4 main figures", §24, §26–§28.** Outside the
  brief; the two stale items found are listed in §0.6.
- **Any file other than `SASP_Kernel_Master_Plan.md` and this report.** In
  particular `figures/`, `results/`, `code/`, `genesets/` and the other reports are
  untouched; `figures/.committed_manifest.json` was not re-snapshotted and
  `check_figures_guard.py --snapshot` was not run.
