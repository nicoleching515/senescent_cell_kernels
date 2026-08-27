# Phase 8+ — roadmap and live status

**Created:** 2026-08-27. **Supersedes** the "Phase 7" numbering for all *remaining* work.
**Last updated:** 2026-08-27 11:20 UTC — ALL QUEUED WORK COMPLETE; awaiting PI

> ## ON RETURN, READ THESE FOUR, IN THIS ORDER
>
> 0. **`reports/COMMIT_PLAN_FOR_PI.md`** — 590 changed paths, commit order, and the
>    destructive-`git checkout` hazard to clear first.
>
> 1. **`reports/COMPLETED_TASKS.md`** — 77 verified entries with the evidence file
>    for each, plus every correction made to my own earlier statements.
> 2. **`reports/CORRECTIONS.md`** — what moved and why (final pass in progress).
> 3. **`reports/SUBMISSION_PATCH_2026-08-29.md`** — **the deadline item.** The
>    manuscript is not in this repo; this must be applied by hand.
>
> **⚠ A `git checkout` IS CURRENTLY DESTRUCTIVE.** `git show HEAD:code/build_genesets.py`
> contains **none** of today's guards — restoring `code/` from HEAD reinstates the
> version that silently overwrites `genesets/*.txt` with EMPTY Tier B modules.
> And `results/phase3_pre_c6/` — the **sole copy** of the baseline the corrections
> ledger compares against — is **98 files, 0 tracked**.
>
> **The one thing only you can do:** **482** untracked + 107 modified files hold the
> entire Phase 8 evidence base — `genesets/human/`, `genesets/mouse_c6/`, most of
> `results/`. **Nothing is committed, so `phase8-frozen` would reference nothing
> immutable.** `PREREG_PHASE8.md` is complete apart from tag hashes and is waiting
> on exactly this.
>
> **The headline did not move.** Controlled amplitude **0.029** against a
> detectable bound of **0.183**; SF **0.088**. §18 outcome **A** stands. Every
> correction today moved against interest and none changed the conclusion.



> ## ⚠ Deviation-ID namespace collision — READ BEFORE CITING ANY "D" NUMBER
>
> **`D1`–`D14` currently mean two different things in two different series.**
>
> - **This file's `D` series = PI DECISIONS** (D5 = promote `mouse_c6`, D15 = implement
>   the composition-matched protocol, D16 = freeze the 100 µm window, …).
> - **`PREREG_PHASE8_genesets.md` §12's `D` series = GENE-SET DEVIATIONS**, an
>   unrelated numbering that also runs D1–D17.
> - `PREREG_PHASE8.md` uses a third series, **`P1`–`P22`**, for its own deviations.
>   That one does not collide.
>
> **Renumbering is queued, not done.** Five agents are still appending rows to these
> documents, so renaming now guarantees rework. It will be done once the runs land,
> before the freeze tag — the gene-set series will move to a `G` prefix, leaving `D`
> to mean "PI decision" everywhere.
>
> **Until then: always name the series, never cite a bare `D<n>`.**

## Renumbering

The Phase 7 addendum (`Phase7_Minimal_Human_Replication (1).md`) defined five steps.
Steps 1-3 are complete and keep their Phase 7 names (`CS_PHASE7_C1.md`,
`BIO_PHASE7_JobA.md`, `BIO_PHASE7_JobA_FOLLOWON.md`) because they are already
written and cross-referenced. **Everything remaining is renumbered:**

| New | Was | Content |
|---|---|---|
| **Phase 8** | P7 Step 5 (freeze) + new | M1 comparability re-run, freeze, pre-registration. Tag `phase8-frozen`. Pre-reg file `PREREG_PHASE8.md` |
| **Phase 9** | P7 Step 4 | H1 acquisition audit A1-A8, Job B (cell types + 4 sender callers) |
| **Phase 10** | P7 Step 5 (run) | H1 through the frozen pipeline, two-arm comparison, Figures 5 and 6 |

**Nothing in Phase 9 or 10 may begin until Phase 8's tag and pre-registration are committed.**

## Status board

Legend: DONE / RUN (in progress) / WAIT (blocked) / TODO

### Completed

| Item | Status | Note |
|---|---|---|
| Environment rebuild after container reset | DONE | Whole scientific stack was lost a 2nd time; rebuilt from pins. `kneed` + `openpyxl` were missing from `requirements.txt` |
| C1 - N3/N4 in-tissue nulls, 3 variants each | DONE | Contribution 3 survives. N3-tile 0.974 vs published 1.000 |
| C1 - destructiveness diagnostics | DONE | `results/phase3/null_destructiveness.csv`. Bbox bug worse than recorded: **35.5% out of tissue** for N3 and 19.9% for N4 (`1 - frac_in_occupancy`), not ~20%. The often-quoted 23%/8% is a *different* column, `1 - frac_retaining_a_neighbour` = shifted senders left with no real cell inside the 100um window (22.8% / 8.0%). Audit item R3 - do not merge the two in one sentence |
| C1 - `make_figure4.py` de-hardcoded | DONE | Byte-identical figure verified |
| Job A - human Tier A-E + gate | DONE | §10 Tier A FAILS the gate; replaced |
| Job A follow-on - spleen markers, A6 spec, C6, E2, SenePy | DONE | 22 spleen cell types; SenePy has no spleen signature |
| GEO screen, 132 series on GPL33762 | DONE | `PHASE8`-relevant record in `reports/PHASE7_H1_SCREEN.md` |
| H1 acquisition (7 samples, 525 MB) | DONE | Panel verified 5,093 genes on data. **Structural checks only** - freeze intact |
| **C7/D1 - DeepScence 11/11 section coverage** | **DONE** | All 9 new sections scored, 1.47M cells. 3 OOM-killed at 5-way concurrency, recovered sequentially. Original 2 preserved unchanged |
| Bio - frozen sets, both arms, gate PASSES | DONE | Strict Tier A 33 on BOTH arms; all 7 modules >=30; disjoint. Independently reverified |
| Mouse C6 counterparts (`genesets/mouse_c6/`) | DONE | Mouse B7 38->108, strict Tier A 25->33. Built from the EXISTING pinned mouse MSigDB, no re-download |
| `PREREG_PHASE8_genesets.md` + FROZEN_MANIFEST | DONE | 35 frozen / 8 variants, SHA-256 per file, 16-row deviation table |

### Running

| Item | Status | Note |
|---|---|---|

| CS - close 4 C1 gaps (N7, full_sf, Fig 2c, §17 footnotes) | RUN | N7 matters most: it varies the sender definition, which just changed |


### Figures (§19)

Palette is fixed: `code/sasp_palette.py` + `apply_style()`. Every figure writes a
`*_data.csv` alongside so every plotted number is auditable. `.png` AND `.pdf`
(compare PNGs for reproducibility - matplotlib date-stamps PDFs).

| Figure | State | Note |
|---|---|---|
| Fig 1 - synthetic identifiability map | DONE, unchanged | §19 says unchanged |
| Fig 2 revised - null battery + corrected N3/N4 + destructiveness | DONE | 2a/2d byte-identical; 2b/2c revised; 2e-h NEW. lambda-hat + 100um lines drawn; N3-swap flagged as = N1; 2h an explicit PENDING box |
| Fig 3 revised - controlled kernels, both arms | WAIT | Needs H1 |
| Fig 4 revised - CCC tools, de-hardcoded | DONE | Byte-identical PNG verified |
| **NEW: intersection / disjointness matrix** | DONE | §11 requires it in Methods. A x B all-zero panel + B x B where C6's cost shows |
| **NEW: cross-arm gene-set symmetry** | DONE, **corrected 2026-08-27** | Pre/post C6, ortholog-intersected. **27**-of-33 Tier A overlap must be visible (26 by the pinned MGI map + `CDKN2B`, which both arms carry and the map has no row for). B7: **88** of 108/116, not 85. `_data.csv` now carries the plotted headline numbers (was JSON-only) |
| **NEW: CoreScence circularity, mouse 79 -> 88%, human 76 -> 88%** | DONE, **corrected 2026-08-27** | Real cost of C6, disclosed. Was labelled "69 -> 76 -> 88%"; the mouse anchor was a typed-in literal and is now derived by `code/corescence_circularity.py` |
| **NEW: SenePy spleen coverage** | DONE | 22 types: 0 matched / 15 surrogate / 7 none. Evidence for the §15 deviation |
| Fig 2h - negative-control-probe kernel | RUN (mouse half) | M1 bundles already ship 40 neg-control probes, so A7 mouse side needs NO download. H1 half stays behind the freeze |
| **Fig 5 NEW - two-arm replication** | WAIT | Phase 10 |
| **Fig 6 NEW - DeepScence native vs remapped** | WAIT | 6a partly unblocked by D1; rest needs H1 |

### Phase 8 - remaining before the freeze

| # | Item | Owner | Blocked by |
|---|---|---|---|
| 8.1 | ~~Mouse re-sourced B7 + rebuilt mouse strict Tier A~~ | Bio | **DONE** |
| 8.2 | ~~Two-arm gene-set symmetry table~~ | Bio | **DONE** - symmetry NOT achievable, quantified |
| 8.3 | ~~Regenerate 5 caller tables at 11-section coverage~~ | CS | **DONE**. Producer of 2 of the 6 tables was NEVER COMMITTED; reconstructed as `code/caller_disagree_all.py`, `--verify` reproduces all 6 committed 2-section tables exactly |
| 8.4 | **Gate: did the caller-agreement headline move?** | CS + Bio | **DONE - YES, IT ROSE ABOVE CHANCE.** Pooled 1.03x (p=0.20) -> 1.118x (p=1.4e-30). "Statistically independent" must go. See below |
| 8.5 | C7/D2 - resolve `denoise=False` | CS | **DONE - DCA INSTALLED.** §6 path 1 landed: DCA 0.3.4 / TF 2.4.4 in an isolated py3.8 venv, `denoise=True` run on 2 full sections. Cost of `denoise=False` is large and the OPPOSITE sign to §4 (D-b): denoising RAISES the depth loading (rho 0.32->0.53, 0.39->0.64), makes the top-5% calls ~100% hepatocyte, and is seed-unstable (1 of 3 seeds gave a disjoint sender set). Recommend freezing `denoise=False` as primary with `denoise=True` as the published-default sensitivity. See `reports/CS_PHASE8_D2_DENOISE.md` |
| 8.6 | ~~C7/D3 re-anchor~~ | CS | **DONE**. `Lmnb1` (which §7 proposes) is itself in 2 Tier B modules - used a disjoint 8-gene proliferation set. Re-anchor = 1 sign bit/section, D1 untouched |
| 8.5b | ~~A7 neg-control-probe kernel, MOUSE half~~ | CS | **DONE - GO, conditionally.** Raw assay is NOT flat on the **POOLED** control features (`all_controls` = 40 probes + 609 codewords + 21 genomic controls): -0.070 SD, p=0.023. **The 40 negative-control probes ALONE are flat (-0.018 [-0.045, +0.010], p=0.183)** - and they are the pre-registered primary A7 null (`PREREG_PHASE8_genesets.md` sec 11, Phase 9 item 9.4), so A7 passes on its own primary response; the gradient sits in the codewords (-0.055, p=0.039) and genomic controls (-0.034, p=0.0039). N5 removes it; N2 does NOT. Never report a naive or N2-only kernel. **Never call -0.070 a negative-control-PROBE number** (`AUDIT_PHASE8_FACTCHECK.md` R1; `results/phase3/a7_summary.csv`). **A7 was run 05:19 on PRE-C6 sender calls and must be re-run after 8.7** (`PREREG_PHASE8.md` P2) |
| 8.7 | M1 end-to-end re-run | CS | **RUNNING.** D5 unblocked it: `pre-c6-genesets` tagged, C6 sets promoted, gate re-verified PASS |
| 8.8 | `reports/CORRECTIONS.md` - what moved and why | CS | folded into 8.7 |
| 8.9 | `PREREG_PHASE8.md` + `git tag phase8-frozen` | **PI** | **PREREG IS DRAFTED AND COMPLETE** — 882 lines, only tag hashes outstanding. Blocked solely on the PI creating the tag and committing the (currently untracked) evidence base |

### Phase 9 - H1, after the freeze

| # | Item | Owner |
|---|---|---|
| 9.1 | Audit A1 (resolution, segmentation, assignment rate), A4 (Ripley's K), A8 (cross-arm on ortholog-intersected panel) | CS |
| 9.2 | **A2 gate** - disjointness on the real panel (pre-verified: passes) | CS |
| 9.3 | **A5 gate** - matched-decoy contrast, \|SMD\| <= 0.1 | CS |
| 9.4 | **A7 - negative-control-probe kernel, must be flat** | CS |
| 9.5 | A6 - build + validate the red/white pulp covariate from the spec | Bio + CS |
| 9.6 | Job B step 1 - cell types with the 22-type spleen marker set | Bio |
| 9.7 | Job B - cross-check against the depositors' 4-level annotations | Bio |
| 9.8 | Job B step 2 - 4 sender callers (DeepScence native, SenePy w/ surrogates, CDKN1A+, Tier A) | CS |
| 9.9 | A3 - prevalence per cell type, 1-20% band | CS |
| 9.10 | Caller agreement conditioned on type + depth decile | CS |

**Stop conditions:** A2 or A5 failure -> §18 outcome C, a data-availability finding.

### Phase 10 - run and compare

| # | Item |
|---|---|
| 10.1 | H1 through the frozen pipeline: naive, N1-N8, controlled fits, kernel families, superposition vs nearest, proximal vs downstream |
| 10.2 | Composition-matched reruns, 5 seeds, both arms |
| 10.3 | §8 DeepScence native-vs-remapped comparison against the pre-registered prediction |
| 10.4 | §17 two-arm table |
| 10.5 | Figures 5 and 6 new; 2, 3, 4 revised |
| 10.6 | Claim audit, citation verification |

## THE 8.4 GATE RESULT — the motivating claim must be restated

Independently reverified against `results/phase3/caller_coverage_gate{,_headline}.csv`.

| Basis | 3-pair band | Pooled | z | p |
|---|---|---|---|---|
| 2-section (published) | **0.932-1.221** (= the published 0.93-1.22x) | 1.03x | 1.27 | 0.20 |
| 11-section (D1) | **0.700-1.711** | **1.118x** | 11.5 | 1.4e-30 |
| 6-section (in-band) | 0.775-1.374 | 1.115x | 8.99 | 2.6e-19 |

- **Tier A vs DeepScence is above chance in 11 of 11 sections** (pooled 1.248). That single pair breaks the independence claim.
- Tier A vs SenePy stays **below** chance in 11 of 11 (0.914).
- **SenePy vs DeepScence FLIPS**: 1.693 at 2 sections -> 0.737 at 11. The old base was not merely underpowered, it was unrepresentative.
- The circular DeepScence-`Cdkn1a`+ pair (excluded from every pooled number) is *weaker* than published: median 1.071 vs the quoted 1.51-2.85x.

**Effect is small, confidence is decisive.** "Statistically independent" is dead;
"weakly but genuinely dependent, in a direction each pair's depth loading
explains" is the defensible restatement (drafted verbatim in `CS_PHASE8_CALLERS.md` §3).

## A7 gave a second unplanned result

The estimator's **false-positive rate is 9-13% against a 5% nominal** - the first
direct measurement, obtained free from the control features. (The often-quoted
"9-16%" spans five responses that are **four overlapping views of one quantity**:
`all_controls` is the sum of the probe, codeword and genomic responses and
`neg_probe_rate` is a ratio of two of them. The 16% upper end is `neg_probe_rate`
alone - the one response whose denominator is an N5 column and which is not a
clean null. Clean-null subset: 0.091 / 0.103 / 0.109 / 0.127.
`AUDIT_PHASE8_FACTCHECK.md` M3, `results/phase3/a7_summary.csv`.)

## Figure policy (set by the PI 2026-08-27)

**Committed figures are NOT to be changed.** `figures/` is held at its committed
state; regenerated versions live in `figures/revised_candidates/` with `_REVISED`
suffixes and a README carrying the regeneration ledger. Everything is regenerated
once, from the frozen configuration, at 8.7.

**Guard scope, stated so it is not over-read:** `code/check_figures_guard.py` enumerates via
`git ls-files figures/`, so it covers the **27 committed** figure artefacts of 45 on disk. The 18
**untracked** outputs - `figure2e.{png,pdf}`, all of `figure_gs1`-`gs4`, `figure_phase8_callers`,
`figure_phase8_d3` and their `_data.csv`s, i.e. precisely the new Phase 8 figures - are outside it,
and `figures/.committed_manifest.json` is gitignored and re-baselineable with `--snapshot`.
"`OK: all 27 committed figures match`" means the committed set is intact; it is **not**
"all figures verified" (`AUDIT_PHASE8_FACTCHECK.md` M7). Committing the Phase 8 figures before the
freeze tag is what would close this.

**Note:** an agent regenerated `figure2b/2c` after the baseline was restored, and
restated them as its own outputs. Restored again; PNGs verified identical to the
archived revisions, so nothing was lost. **More than one actor writes `figures/`
and nothing warns** - this needs a guard before 8.7.

## PI DECISIONS TAKEN — 2026-08-27

| # | Decision | Consequence |
|---|---|---|
| **Direction** | **STAY on the original research question.** "How far does senescence signalling reach?" is answered **in the negative, as a bound**: naive amplitude 0.326 -> controlled 0.027; SF 0.082 [-0.099, 0.249]; nothing above the 0.203 response-sd detectable bound at 80% power | The paper keeps its question. The audit battery (A7, corrected nulls, FPR 9-16%) is supporting evidence for *why* the bound is the honest answer, not a replacement thesis |
| **D2** | `denoise=False`: **both paths in parallel** - measure the cost now, attempt DCA opportunistically | Freeze is never blocked on DCA. A measured caveat + a documented install failure satisfies §6 |
| **D5** | **Tag then promote.** `git tag pre-c6-genesets` created; C6 mouse sets promoted into `genesets/` | Only 3 files changed: B7 38->108, strict Tier A 25->33, A_sender_for_secondary_senescence 55->74. Gate re-verified PASS on the authoritative 5,097-gene panel. 8.7 unblocked |
| **D4** | Age is a **continuous covariate only**. No young-vs-old contrast, no age-stratified prevalence claim | Honest about n=2 above 55. H1's value is a human replication of the geometry, not an ageing result |
| **D1** | Tier A: **strict-33 PRIMARY**, per-module sets as pre-registered sensitivity | Both carried through 8.7 (+~15 min, not a second pass) |
| **D6** | CoreScence 88% circular: strip-and-refit is primary | Settled by evidence once independence also fell |
| **D7** | Mouse panel = **5,097** (stock 5,006 + 91 add-on genes; 9 genotyping probes excluded) | Resolved. B6 margin is exactly 1 gene (`Junb`) |

## Still open for the PI

**All decisions are now taken.** D8/D10/D2b/D3b implemented; D11/D12/D13/D14/D9 decided 06:15.

| # | Decision taken | Action |
|---|---|---|
| **D11** | **Patch the claim, submit Aug 29 as planned** | `reports/SUBMISSION_PATCH_2026-08-29.md` — search strings, replacement text, verified numbers. **The manuscript is not in this repo; the PI must apply it by hand.** Framing note: the restatement is STRONGER than what it replaces — "independent" was an accept-the-null argument on 2 sections; 1.12x at p=1.4e-30 with one pair significantly below chance is a positive finding on 11 |
| **D12/D13** | **Hold paper-level edits until the runs land** | §29 objections 3 and 6, and the §30 outline, stay stale by choice so nothing is written against moving numbers. **Re-raise once 8.7 and 8.5 report** |
| **D14** | **Re-run N7 at 1,000 permutations** | §24.3 compliance. Message sent to the running M1 agent to catch it before that stage; ~55 min vs 10.7 at 200 |
| **D9** | **Isolated `/workspace` venv is sufficient** | Docker unavailable in-pod. `/workspace/envs/sasp311` building, fully isolated, on the network volume per master plan §16.1 |

## Superseded decisions list

| # | Decision | Default if unanswered |
|---|---|---|
| D1 | Which Tier A is primary - **proposed strict-33**, per-module as sensitivity | Proceed as proposed |
| D2 | Plasma cells: 3 markers vs `MIN_MARKERS=4` | Excluded (fails by one gene) |
| D3 | Marginal zone B cells pre-registered as a named receiver type on 6 CellMarker spleen rows? | Included, flagged as weak evidence |
| D4 | Is a 2-donor-over-55 age axis enough to pre-register an age-stratified claim? | Report age as continuous, no stratified claim |
| D5 | **Promote `genesets/mouse_c6/` into `genesets/`?** Doing so invalidates every Phase 2-5 mouse result, which were computed pre-C6 | Not promoted; M1 re-fit required either way (8.7) |
| D6 | CoreScence circularity rose under C6 (as posed: "69% -> 88%"; **corrected 2026-08-27**: mouse 79% -> 88%, human 76% -> 88% - the 69% was a typed-in literal, `AUDIT_PHASE8_FACTCHECK.md` M1). Accept as a disclosed cost, or strip-and-refit? | Strip-and-refit added to the frozen run order. Decision unaffected; the *magnitude* of the cost is roughly half what was stated |
| ~~D7~~ | **RESOLVED.** Mouse panel = stock 5,006 + 100-gene add-on = 5,106 h5 features; minus 9 genotyping probes = **5,097 authoritative**. Verified: the two CSVs are disjoint and their union IS the h5 feature list. B6 = 31, margin of exactly **1** (`Junb`, from the add-on) | Asserted in code; numbered pre-reg item D10 |

## Known risks, live

| Risk | State |
|---|---|
| Container reset wipes the stack a 3rd time | **Materialised twice.** Pins now being completed; consider a baked image |
| Freeze not committed while H1 data sits on disk | **Open now.** Structural reads only so far |
| C6 adopted asymmetrically breaks §17 | **Being fixed** (8.1) |
| SenePy is not the same estimator across arms | Declared deviation; no fix exists |
| Full DeepScence coverage moves the caller-agreement headline | **Unknown until 8.4.** Now UNBLOCKED - D1 is complete. This is the one that could force a restatement of the paper's motivating claim |
| **`code/make_figure2.py` is a superseded dual producer of figure 2a** | Running it silently replaced the committed figure while its data CSVs came back identical, so nothing warned. Live producer is `make_phase5_figs.py --which 2a`. **Do not run `make_figure2.py` during the M1 re-run** |
| `code/build_genesets.py` cannot be re-run as committed | **CLOSED IN THE WORKING TREE, still open at `HEAD`.** `HEAD`'s copy has `SCRATCH` pointing at a dead per-session /tmp path; re-running *that* copy would glob 0 JSONs and overwrite `genesets/*.txt` with EMPTY Tier B modules. The uncommitted working-tree copy (`M code/build_genesets.py`) defaults `SCRATCH` to `/workspace/genesets/msigdb_mouse_2026.1.Mm` and raises `SystemExit` with an explicit message if it globs zero JSONs. **The risk is therefore live only until that edit is committed** (`AUDIT_PHASE8_FACTCHECK.md` R8) |
| **B6 `oxidative_stress` has a margin of exactly 1 gene** on the mouse arm (31 vs floor 30) | The gate must be re-run after ANY gene-set or panel change, not only at freeze. It exits non-zero so it cannot pass silently |
| 3 files in the pinned mouse MSigDB archive are HTML error pages, not JSON | Confirmed by the bio agent, hidden since 2026-08-20 by a bare `except: pass`. No tier uses them, so no published number is affected |
