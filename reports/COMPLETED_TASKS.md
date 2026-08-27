# Completed tasks ledger — Phase 8

Appended as work lands, so the PI's return summary is **evidence-backed rather
than reconstructed**. Every row names the artefact that proves it.

**Session opened 2026-08-27 ~03:05 UTC. PI away 06:30–~14:30 UTC** with authority
delegated for decisions, judged on preserving the research mission and on
cost-effectiveness.

---

## Completed and independently verified

| # | Task | Evidence | Verified how |
|---|---|---|---|
| 1 | **Environment rebuilt after 2nd container wipe** | `logs/env_rebuild_2026-08-27.log` | Whole scientific stack was gone (numpy back to base 1.26.3). Rebuilt from pins |
| 2 | **`kneed` + `openpyxl` added to `requirements.txt`** | `requirements.txt:46-57` | `kneed` is a HARD DeepScence import (`io.py:17`) — the pins could not rebuild a working DeepScence. Verified by import |
| 3 | **GEO scope: 132 GPL33762 series screened** | `reports/PHASE7_H1_SCREEN.md` | E-utilities; 19 with real Prime-5K evidence. Screen validated against a known answer (GSE335761 correctly fails at 386 targets) |
| 4 | **H1 selected: GSE326743**, 7 human spleens | `reports/PHASE7_H1_SCREEN.md` | Only normal-tissue Prime-5K deposit with >=3 donors in all 132 |
| 5 | **H1 panel verified ON DATA: 5,093 genes** | `genesets/h1_candidate/` | Byte-identical across 3 samples (symmetric diff 0); 609/40/21 controls match the §12.4 expectation |
| 6 | **H1 acquired, 28 files / 525 MB** | `code/fetch_h1_geo.sh`, `data/raw_h1/` | 21 transcripts/morphology files skipped per §12.2. cell_id sets match h5<->parquet in all 7; coords in microns |
| 7 | **DeepScence D1: 2/11 -> 11/11 coverage** | `data/processed/deepscence_*.csv` | 1.47M cells scored. 3 sections OOM-killed at 5-way concurrency, recovered sequentially. Original 2 preserved unchanged |
| 8 | **C1: in-tissue N3/N4 nulls, 3 variants each** | `reports/CS_PHASE7_C1.md` | Bbox bug real and worse than recorded (23% void, not ~20%). **Contribution 3 survives**: N3-tile 0.974 vs published 1.000 |
| 9 | **C1: N3-occ shown degenerate; N3-swap == N1** | `results/phase3/null_destructiveness.csv` | occ admits only near-identity offsets (27um vs median lambda-hat **14.7um** — this row originally said 12.8um, a FOURTH value then in circulation, corrected 2026-08-27). swap: rho=0.948 vs published N1 |
| 10 | **N7: correction invariant to sender definition** | `results/phase3/sf_summary_c1_n7.csv` | N3-tile 0.974-1.006 across 6 callers, prevalence 0.5-9% |
| 11 | **`make_figure4.py` de-hardcoded** | `code/make_figure4.py` | All 7 constants re-derive; regenerated PNG byte-identical |
| 12 | **Job A: human Tiers A-E built** | `genesets/human/`, `reports/BIO_PHASE7_JobA.md` | §10's 16-gene Tier A **FAILS the gate** (14 on-panel, 13 inside Hallmark modules, ATR alone survives). Independently recomputed |
| 13 | **Spleen marker set, 22 cell types** | `code/markers_human_spleen.py` | Every gene traced to PMIDs in a pinned CellMarker 2.0 download |
| 14 | **A6 covariate specified (red/white pulp)** | `genesets/human/D_spleen_*.txt` | Replaces §13's lung airway-alveolar axis. Build deferred — needs the matrix, behind the freeze |
| 15 | **SenePy has NO spleen signature** | `results/phase7_jobA/senepy_spleen_coverage.csv` | Verified directly: 65 human hubs, 10 tissues, liver present, spleen absent. 7 of 22 cell types unscoreable |
| 16 | **C6: B7 re-sourced, both arms** | `genesets/`, `genesets/mouse_c6/` | §23's stated method FAILS (12/24 genes, below floor). Re-sourced passes: mouse 108 / human 116 |
| 17 | **Cross-arm symmetry quantified as NOT achievable** | `results/phase7_jobA/crossarm_geneset_table.csv` | Root cause: `REACTOME_SASP` is 40 mouse vs 111 human genes. The two 33-gene Tier A sets share only 26 |
| 18 | **Mouse panel resolved: 5,097** | verified from `.h5` + both CSVs | stock 5,006 + 100 add-on = 5,106 features, disjoint, union == h5 exactly; minus 9 genotyping probes |
| 19 | **Gate PASSES both arms** | `results/phase7_jobA/gate_result_human.json` | Independently recomputed: A=33 both arms, all 7 modules >=30, all disjoint |
| 20 | **`pre-c6-genesets` tag + C6 promotion** | `git tag pre-c6-genesets` | Only 3 files changed. Baseline proven recoverable byte-identically |
| 21 | **8.4 GATE ANSWERED: independence FALSIFIED** | `results/phase3/caller_coverage_gate*.csv` | Like-for-like 1.030 -> **1.118, p=1.44e-30**. The paper's motivating claim must be restated |
| 22 | **A7 mouse: raw assay is NOT flat** | `reports/CS_PHASE8_CALLERS.md` §4 | -0.070 SD (p=0.023). N5 removes it; **N2 does NOT**. Never report a naive or N2-only kernel |
| 23 | **FPR measured: 9-16% vs 5% nominal** | `results/phase3/a7_summary.csv` | Free from the negative-control probes. Magnitudes PROVISIONAL (pre-C6) |
| 24 | **D3 re-anchor; `Lmnb1` rejected** | `reports/CS_PHASE8_CALLERS.md` | §7's proposed anchor is itself in 2 Tier B modules. 7250 is the sole disagreeing section — 4 published anomalies converge there |
| 25 | **4 silent-failure hazards fixed** | `code/build_genesets.py`, `make_figure2.py`, `check_figures_guard.py` | Each falsification-tested, not asserted |
| 26 | **Figure policy enforced** | `figures/revised_candidates/`, `.committed_manifest.json` | 27 committed figures held; guard proven to catch drift and clear on restore |
| 27 | **8 figures produced** | `figures/figure_gs1-4`, `figure2e` | All from existing CSVs, repo palette, `*_data.csv` beside each |
| 28 | **Container-wipe root cause found** | master plan §16.1 | `/` is an ephemeral overlay; `/workspace` is the network volume. The plan says install to `/workspace` **in bold** — it was not followed |
| 29 | **`.gitignore` gaps closed** | `.gitignore` | `envs/` (6.6GB), `__pycache__` (719 dirs), `data/raw_h1/` (525MB) were all committable |
| 30 | **Gate wired to fire on any change** | `code/gate_genesets_guard.py`, `.githooks/pre-commit` | Falsification-tested by monkeypatch |
| 31 | **`PREREG_PHASE8.md` drafted, 770 lines** | `reports/PREREG_PHASE8.md` | Every fixed parameter read from code with file:line. Tag hash left explicitly TBD |

| 32 | **8.7 M1 end-to-end re-run — COMPLETE** | `reports/CS_PHASE8_M1_RERUN.md`, `CORRECTIONS.md` | 05:45-09:08 UTC (~3h25m), all stages. `a7_summary.csv` and `sf_summary_c1.csv` now 09:06 (post the 05:48 sender rebuild). `n_perm=1000` in all six perturbation-null files. Verified directly |
| 33 | **N7 at 1,000 perms — DONE** | `results/phase3/perm_nulls_c1_n7.csv` (07:47) | `n_perm` now reads **1000** across all 5 non-primary calls. §24.3 satisfied. Verified directly |
| 34 | **`caller_coverage_gate.csv` mixed-basis bug FIXED** | now carries 4 explicitly-labelled bases | Was pooling a pre-C6 2-section base against post-C6 11-section rows |
| 35 | **"22 of 33 above chance" corrected to 20 of 33** | `CS_PHASE8_M1_RERUN.md` §0.3 | Earlier report was wrong |

## The three results that matter most

**⚠ These three are from a partially-superseded state and are PROVISIONAL until
the re-run finishes and the numbers are re-verified.**

**1. The headline did NOT move — the negative result is robust.**
Controlled amplitude **0.029** response-sd against a design detectable bound of
**0.183** at 80% power (pre-C6: 0.027 vs 0.203). SF under N2+N5+N6 **0.088**
[-0.017, 0.234] (pre-C6 0.082). **§18 outcome A unchanged. C1 verdict unchanged.**
Every correction moved against interest and none changed the conclusion — which
is the strongest position a negative result can be in.

**2. The independence claim was an ARTEFACT OF A CIRCULAR SENDER DEFINITION,
not a power problem.** On the *published two-section base*, once the circular
Tier A is replaced, agreement is already **1.13x at p=4e-8**. The claim does not
survive even without the coverage fix. This is a sharper and more damaging
finding than "more sections revealed dependence" — the original analysis was
wrong on its own data, because its sender set was 13/14 circular with the
response modules. There is **no longer any caller pair reliably below chance**.

**3. Circularity was manufacturing signal.** The re-sourced B7
`secondary_senescence` — the module whose name promises a paracrine effect — lost
its gradient once it was made *actually* paracrine: 36 -> **22** reportable fits,
naive amplitude 0.342 -> **0.246**. The pre-C6 B7 shared 27 arrest genes with the
sender-adjacent literature. **The gradient it showed was partly a gradient in
arrest.** B7 is the only module whose response score changed and it accounts for
the entire net loss of reportable fits (160 -> 153).


| 36 | **Independent novelty assessment, literature-verified** | `reports/NOVELTY_ASSESSMENT.md` (592 lines) | Overturned three of the project's own novelty claims — see below |
| 38 | **Independent fact-check audit** | `reports/AUDIT_PHASE8_FACTCHECK.md` | **Science reproduces; reporting layer does not.** Gate, panel arithmetic, C1 battery, A7 stats all re-derive exactly (66/66 SF cells, 30x4 A7 stats, 6 figure->data chains, 0 mismatches). **9 claims refuted, 10 true-but-misleading** |
| 39 | **A7 response-naming error corrected** | `SUBMISSION_PATCH_2026-08-29.md` §4a | Manuscript-bound. See below |
| 37 | **Caller finding decomposed into its two causes** | `caller_coverage_gate_headline.csv`, 6 labelled bases | Circular Tier A fix: 1.04 -> 1.131 (p=6.5e-9) at 2 sections. Coverage: 1.131 -> 1.212. **Circularity, not power, is the main cause** |

## Novelty assessment — three claims overturned before submission

| Claim | Verdict |
|---|---|
| "Nobody reports the negative-control-probe spatial diagnostic" | **FALSE.** Voyager's Xenium vignette computes Moran's I on negative control probes; Ren et al., *Nat Commun* 16 (2025) does it peer-reviewed. **What IS novel:** fitting *the estimand's own estimator* to the controls (a Lipsitch 2010 negative-control OUTCOME, not generic QC), and the **N2-vs-N5 result** — matched decoys leave the technical gradient ~80% intact while covariates remove it. Rated the project's **strongest** contribution, and it contradicts the master plan §23's claim that the decoy contrast is "the single most important number in the paper" |
| "Torus degeneracy on non-convex tissue is new" | **Novel in spatial omics, 1982-vintage in statistics.** Lotwick & Silverman (1982); Mrkvicka et al. (2021) proposes the variance correction *specifically for irregular windows*. Implementation dispatched — using the field-standard fix removes the obvious reviewer objection |
| "Nobody reports a measured FPR" | **No longer safe** — CellWHISPER (Jan 2026) and CONCISE (Jun 2026) report calibration. Ours is distinct in being measured from assay-internal nulls, but it is 40 probes, only powered pooled: one Methods sentence, not an abstract claim |
| **Verified citation error** | **CellWHISPER is NOT the source of the torus-shift null.** Its null is a within-cell-type location permutation — this project's **N1**. N3 is what Figure 4 rests on. Citation audit dispatched |
| **Bibliography gap** | `references.bib` has 30 entries and **zero spatial-statistics methods papers**, though the project is structurally an instance of the *spatial confounding* literature |
| **Central negative result** | **INTACT.** Searched and found **no published length constant for senescence spatial influence** — the claim has no prior to contradict |
| **Venue** | The plan's primary/secondary is **backwards** for a negative result. ICBINB-BIO matches topic-by-topic and gives 8 pages; ml4spatialbio is 4 pages, non-archival, permits concurrent submission, window to ~Sept 4 |

| 40 | **Citation audit — bibliography was written from recall** | `reports/CITATION_AUDIT.md`, `references.bib` | **19 of 32 entries carried invented author forenames — 41 wrong given names**, verified against Crossref AND PubMed independently (they agree). 49 author lines corrected |
| 41 | **CellWHISPER misattribution corrected in 9 places** | plan §4.3/§22/§23/§31/§32, `CS_PHASE1/2/3.md`, `references.bib` | CellWHISPER's null is a within-cell-type location permutation (this project's **N1**). "torus"/"toroidal"/"wraparound" appear **nowhere** in the paper. Correct attribution: **Lotwick & Silverman (1982)** |
| 42 | **11 spatial-statistics references added** | `references.bib` | The gap was structural: the project is an instance of the *spatial confounding* literature and cited none of it. Each entry carries a `% SUPPORTS:` line so nothing is decorative |
| 43 | **">90% FPR" drift corrected** | plan §5, §2.2 | CellWHISPER's figure is an interaction-count ratio from which an FPR is *inferred*. §1/§4.3/§31 correctly said "implying"; §5 had drifted to asserting it |

## The citation audit's two most consequential findings

**1. `Martin, Luke` was actually `Martin, Lucy`.** A real researcher was
misgendered in the bibliography. This was not an isolated typo — it is one of 41
wrong given names across 19 of 32 entries. The pattern (plausible, correctly
transliterated, wrong) is diagnostic of a bibliography written from recall rather
than from retrieved records. An earlier audit (D7) had checked journal, volume,
pages and DOI — but not forenames.

**2. Mrkvicka et al. (2021) §2.1.4 predicts that THIS PROJECT'S TILING MAKES THE
PROBLEM WORSE.** A union-of-aligned-rectangles window — exactly what
`phase4_tiles.py` produces — is predicted to *increase* the liberality of the
torus correction, not fix it. **N3-tile (SF 0.974) is the variant C1 currently
presents as the corrected primary.** Relayed immediately to the running
variance-correction agent with instructions to verify it from the paper and to
lead with the result either way.

**Reframing this forces:** the torus finding is NOT "we discovered torus shifts
break on non-convex tissue" — known since 1982 — but "a 40-year-old documented
limitation is being violated in current spatial-omics practice, and we quantify
what it costs." That is still publishable; it is a different claim.

| 44 | **Audit corrections applied** | `reports/AUDIT_CORRECTIONS_APPLIED.md` | Guard clean before and after (27/27). Agent also **disagreed with the audit on 3 points, with evidence** — the audit's own fix list was internally inconsistent |
| 45 | **`CDKN2B` map-gap fix, extended beyond Tier A** | `code/crossarm_geneset_table.py`, `figure_gs2` | Shared Tier A **27 not 26**; `mouse_only` = 6. Same logic applied uniformly: **B7 goes 85 -> 88 shared**, Tier C was short by `Cxcl1`. New `gap_split()` **asserts** both arms' non-shared members are accounted for by name, so the arithmetic that failed cannot be written again |
| 46 | **CoreScence circularity re-derived from files** | `code/corescence_circularity.py`, `figure_gs3` | See correction below. Gate and figure now **compute** it instead of asserting it |
| 49 | **DCA INSTALLED AND RAN — §6 path 1/2 partially landed** | `data/processed/deepscence_dca_7239_liver_sbr_Male_52-U1.csv` | Full-section DCA **succeeded on 7239**; subsample-20k arm also succeeded; 7259 OOM-killed against the 57.7 GB cgroup. §22 rated install failure Medium-High — **it did not fail**. The `denoise=False` deviation can be *removed* for at least one section, not merely quantified |
| 50 | **M1 re-run on WAVE 2** | `logs/m1_chain.log` | WAVE1 done; now running `tierApm_p95` (the pre-registered per-module sensitivity arm) at 1,000 perms |
| 47 | **Variance-corrected torus null implemented from source** | `reports/CS_PHASE8_TORUS_VAR.md`, `results/phase3/*_var*.csv` | Mrkvicka §2.1.3 + §2.2 (RS_count, RS_ker) implemented from the retrieved full text, quotes in the module docstring. **Nothing improvised** |
| 48 | **Type-I error MEASURED for every torus variant** | `results/phase3/var_sim_calibration.csv` | The paper's own §5 design plus an irregular-window panel. **This is the result that upgrades contribution 3** |
| 51 | **Composition-matched protocol (D15) implemented and run** | `reports/CS_PHASE8_COMPMATCH.md`, `results/phase3/compmatch_*.csv` | Built from near-zero spec — the phrase appears 5 times in the planning docs and nowhere in code. 9 ambiguities resolved and recorded as D15.1-D15.9 |
| 52 | **`PREREG_PHASE8.md` COMPLETE except tag hashes** | 882 lines, 15 sections, 25 P-rows + D15.1-D15.9 | Only 7 `TBD`s remain and **all are tag commit hashes** that cannot exist until the PI creates the tag. The freeze document is ready |
| 53 | **§10 rule 8 added: the matched-decoy number may never be reported alone** | `PREREG_PHASE8.md:704` | Wherever 0.9837 appears, `type_adj` (65.9%) and `typecomp_adj` (85.4%) must appear beside it. Quoting the first without the second states the opposite of what the data says |
| 54 | **The inertness was PREDICTABLE from the project's own Phase 3 report** | `reports/CS_PHASE3.md:435` | Verified: N2 matched decoy returns **0.934 [0.869, 0.988] with a planted effect** and **0.775 [0.564, 0.924] with none** — a gap of 0.159 with **overlapping CIs**. A design that cannot discriminate, documented since Phase 3 and **unconnected to the frozen parameter until now** |

| 55 | **Figure regeneration pass executed — the one the policy allows** | `figures/revised_candidates/README.md` (regeneration record) | Guard before: OK 27/27. After: exit 1 listing exactly figures 2a-2d, 3, 4 and their data CSVs. Headlines did NOT move: controlled amplitude **0.029** vs a detectable bound of **0.183**; SF **0.088**. §18 outcome A stands |
| 56 | **Figure guard extended to untracked artefacts** | `code/check_figures_guard.py` | Was `git ls-files`, so the 18 untracked Phase 8 figures were unwatched — and two WERE rewritten mid-run with no warning. Now walks all of `figures/`: **46 artefacts, up from 27**. Proven to catch untracked drift, then re-snapshotted to the post-8.7 frozen state |
| 57 | **A7 re-run under C6 senders** | `results/phase3/a7_summary.csv` (09:06) | Raw assay still **not flat** (-0.074 SD, p=0.015); N5 removes it; **N2 does not**; FPR unchanged at 9-16%. One caveat *disappeared* (`neg_probe_rate` under N5, p=0.016 -> 0.183) and its absence was recorded rather than left silent |
| 58 | **8.7 closed out: `fig_phase3_*` on the frozen basis** | `figures/fig_phase3_caller_depth.png` (new hash `8a68d699…`) | **My instruction was based on an overstated premise.** Only ONE of the three reads a caller table. The agent checked *producers*, not filenames |
| 59 | **Deferred fact-check items applied** | `reports/AUDIT_CORRECTIONS_APPLIED.md` §4 | R3, R4, R5, R6, R9, M4, README. M6/M9 already resolved by 8.7 (verified) |
| 60 | **N3-var / N4-var adopted as the primary null** | `reports/CS_PHASE8_TORUS_VAR.md`, §17 rows | **N3-var 0.996** [0.975, 1.007], **N4-var 0.985** [0.958, 1.003] vs N3-tile 0.971 and a published 0.999 |
| 61 | **Figure guard re-snapshotted post-8.7** | `figures/.committed_manifest.json` | **52 artefacts** now protected (27 -> 46 -> 52 as untracked figures and new `.pdf`/`_data.csv` outputs were brought in scope) |
| 62 | **Superseded claim-audit suggestion marked, not rewritten** | `reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` | Carried "0.93-1.22x / 1.51-2.85x" as a **live suggested rewording** for a README that has since been corrected. Marked SUPERSEDED with current values; the historical text left intact, because a claim-audit's value is its record of what was proposed |
| 63 | **Aug 29 submission patch re-derived against settled numbers** | `reports/SUBMISSION_PATCH_2026-08-29.md` (630 lines) | Every §1-§3 number re-derived by Mantel-Haenszel pooling **directly from the per-section CSVs**, independently of the gate files — all match |
| 64 | **8.5 / C7-D2 RESOLVED — DCA installed and ran** | `reports/CS_PHASE8_D2_DENOISE.md` | §6 **path 1 landed, not the fallback.** DCA 0.3.4 / TF 2.4.4 in an isolated CPython 3.8.19 venv, out of process. Three full M1 sections + a 20k three-seed panel. Main env verified **still has no TensorFlow**; all 11 committed `deepscence_*.csv` byte-identical |
| 65 | **Repo cleanup (safe items)** | — | **2,425 `__pycache__` dirs removed, 18G -> 16G.** Ignore rules verified by creating a test path, not by testing a directory just deleted |
| 66 | **D2 folded into the pre-registration — framing INVERTED, not softened** | `reports/PREREG_PHASE8.md` (993 lines, 29 P-rows) | §3.9 now opens *"RESOLVED. DCA installed, ran, and lost on the merits."* The old sentence — *"DeepScence as we can run it on this panel, not DeepScence as published"* — is **withdrawn**: we ran it as published and chose against it |
| 67 | **Three new prohibitions added** | `PREREG_PHASE8.md` §10.9-§10.11 | No `denoise=True` number from a single seed without its seed-stability companion; `mor` may never be cited as evidence normalisation cannot move this caller; `rho_signed_dz_vs_depth` on `raw` rows must not be quoted |
| 68 | **§8's confound halved — two new falsifiable H1 predictions** | `PREREG_PHASE8.md` §8 P-vi, P-vii | DCA now runs on **both** arms, so `denoise=False` can no longer explain a mouse-only artefact. Registered: that denoising raises depth loading natively too, and that the seed instability recurs. Each states its own falsifier |
| 69 | **D12/D13: master plan §22, §29, §30 corrected** | `reports/PLAN_UPDATE_D12_D13.md` | §23's "single most important number" struck; objections 3 and 6 rewritten; 1, 4, 7, 9 updated; 2, 5, 8 left verbatim. Venue evidence added as a **flagged subsection**, Primary/Secondary paragraphs byte-unchanged (the PI's call) |
| 70 | **Phase 7 plan marked superseded in part** | `Phase7_Minimal_Human_Replication (1).md` | It still called H1 **"human aging lung"** in 5 places — H1 is 7 human **spleens**. Banner added listing all seven superseded claims; the document is retained as the plan of record |
| 71 | **Moran's I dispatched** | pending | §29 objection 9 **promised an analysis that does not exist** (`grep -ril moran` returns nothing) — and it is the objection answering the one piece of prior art that directly falsifies a project novelty claim |
| 72 | **Moran's I run — §29 objection 9 now answered by measurement** | `reports/CS_PHASE8_MORAN.md`, `results/moran/` | 11 sections, all 13,590 features individually, 7 weights variants, 999 permutations. Producer validated against `esda.Moran` to **0.00e+00**. No package installed |
| 73 | **A falsified recommendation corrected before it was acted on** | `reports/NOVELTY_ASSESSMENT.md` (banner) | That report told the project to argue Moran's I and A7 **disagree**. They do not |
| 74 | **New methodological point against the prior art** | `reports/CS_PHASE8_MORAN.md` | **Voyager's controls-vs-genes contrast is an ABUNDANCE contrast** |
| 75 | **`CORRECTIONS.md` finalised: 659 -> 1,966 lines** | `reports/CORRECTIONS.md` | Front matter A-D + six new sections. Leads with what did NOT move, re-derived by running `code/m1_headlines.py`: controlled amplitude **0.0288** vs an 80%-power bound of **0.1833**, SF **0.0885**, naive 0.3288, 153 fits. The bound **tightened** from 0.203 |
| 76 | **`kumar2026cellwhisper` forenames actually fixed** | `references.bib` | The citation audit's own note said they were corrected. **They were not.** Abhishek->Anurendra, Fernando->Felix, Bhavya->Bhavay, Nan->Nicholas — in the project's most load-bearing citation, while `SUBMISSION_PATCH` told the PI the bibliography was clean |
| 77 | **Falsified Moran framing removed from the DEADLINE document** | `reports/SUBMISSION_PATCH_2026-8-29.md` | §9 still instructed *"so the reader can see the two tests disagree"* — in the one document applied to the manuscript by hand, two days out. Replaced with the measured power argument |
| 78 | **Freeze document corrected before tagging** | `reports/PREREG_PHASE8.md` | **5 corrections.** It carried my "roughly doubles" overstatement in 2 places and 3 **pre-C6** A7 values. Frozen values re-read by me from `a7_summary.csv` (09:06): codewords **-0.0604** (was -0.0549), genomic **-0.0307** (was -0.0337), probes **-0.0225** p=0.129 (was -0.0177). Residual stale digits: **0** || 83 | **`WRITING_PACK.md` built — 1,389 lines** | `reports/WRITING_PACK.md` | 53 marked re-derivations, a 33-item forbidden-claims checklist, 22 tabulated document disagreements with the authoritative file named, and a per-subsection source map with reproduce commands |

## ⚠ THE PACK FOUND A CIRCULAR NUMBER, AND TWO MORE ERRORS OF MINE

**1. `lambda-hat = 15.7 um "pooled"` was UNSOURCED and CIRCULAR — NOW RESOLVED
(see the lambda-hat entry further down: authoritative value **14.7321 um**,
emitted by `summary_phase3.txt` line 90, with mandatory IQR [7.0, 50.0] and 60%
railed). 15.7 traced to a pre-C6 *interior* median over 441 pseudo-replicated
rows. The original diagnosis is kept below as the record of how it was found.**
It appears in six places including §30 5.3 and §29 objection 7. **No file emits
it.** Its only apparent origin is a back-derivation from the torus report's own
*"2,215 um = 141x the pooled lambda-hat"* — i.e. derived from the claim it is used
to support. Verified alternatives, all real: **16.07** um (median over the 153
reportable fits), **14.73** (all 315), **17.1** (interior median — what the
summariser actually prints). Two dependent claims inherit the problem: "141x
lambda-hat" and "seams ~76 lambda-hat apart". **This must be resolved before it
reaches a paper.**

**2. The brackets I have been quoting are IQRs, NOT confidence intervals.**
Verified: `sf_summary.csv` carries `q25 / median / q75` and **no CI column**. So
`0.088 [-0.017, 0.234]` and `0.029 [-0.007, 0.084]` are **interquartile ranges
across fits**. I wrote them as CIs repeatedly in this session, and §30 presents
them unlabelled. An IQR and a CI mean entirely different things.

**3. My Spearman rho was aggregation-dependent and I presented it as definitive.**
I reported **+0.923** and said I had verified it. It does reproduce — but only
under my aggregation (median per field, knn6). The same data gives **+0.8951**
under the pack's method and **+0.7104** per-row with no aggregation. **All three
are defensible; none is "the" value.** The real finding is that the statistic
swings 0.71-0.94 with aggregation choice, so **the method must be stated
alongside it.** The falsification it supports (Moran and A7 do NOT disagree)
survives at every value.

**Also unsourced:** the "76 %" composition-surrogate share (0.212/0.260 = **0.815**;
no denominator on disk yields 76 %) — superseded anyway by the 65.9 %/85.4 % pair.

## Document disagreements the pack tabulated (22 total)

- **`PHASE8_ROADMAP_STATUS.md` contradicts itself** — my banner is post-C6 and
  correct; the PI-decisions table, the 8.4 gate table and the 8.5b A7 block are
  all pre-C6.
- **`PREREG_PHASE8.md` §10.1 mixes vintages inside one item.**
- **§30 5.9 mixes bases in one sentence** — 1.212 (4-pair) paired with 1.128
  (3-pair). Consistent pairs are 1.128 -> 1.212 or 1.131 -> 1.212.
- **Audit R3 was never applied to `CS_PHASE8_TORUS_VAR.md`**, which still says
  "23 % in the void".
- **Audit R5 is now moot** — `neg_probe_rate` under N6+N5 no longer excludes zero
  (p 0.199), so "every control family is flat under +N6+N5" is now *true*.
- **"Naive biological amplitude" has FOUR values in circulation** (0.277 / 0.291 /
  0.312 / 0.314) — two estimators x two vintages.
- Two files still carry live forbidden numbers: `BIO_DELIVERABLE7_CLAIM_AUDIT.md`
  L285 and `CS_PHASE7_C1.md` §6.3.

## ✅ ALL OF THE ABOVE WORKED, 2026-08-27 — `reports/RECORD_RECONCILIATION.md`

The four findings and the 22 disagreements have been resolved against files. The
record above is left as written; the outcomes are:

1. **lambda-hat = 14.7 um.** The pooled median of `lam_naive` over the 315 primary
   fits, printed by `code/summarize_phase3.py:221` into
   `results/phase3/summary_phase3.txt` sec 6, `tierA_p95` row, column `medlam`.
   **15.7 is withdrawn everywhere.** Reasons for choosing the pooled median over
   16.07 / 14.99 / 17.1: it is the only one an emitted file carries; the
   pre-registration deliberately does **not** make lambda-hat an estimand, so the
   tiebreak falls to the frozen code; it shares the 315-fit denominator with the
   60 % railing rate it must always travel with; and it is a valid median of a
   censored sample (the median order statistic lies in the interior) whereas the
   interior median discards 60 % of fits non-randomly and is biased upward.
   **Dependents: 2,215 um = 150x (was 141x); seams ~81 lambda-hat apart (was ~76);
   the 100 um window spans ~6.8 lambda-hat (was "6 lambda"). All three strengthen.**
   The nearest reproducible match to 15.7 anywhere is 15.716 — the pre-C6 *interior*
   median over the zonation-stratified rows, which triple-counts hepatocytes.
2. **The brackets are IQRs and stay IQRs.** The pre-registered bootstrap emits
   *per-fit* CIs only (median span [-0.415, +0.381]); it emits no interval on the
   median across fits, so a genuine CI cannot be computed without a new run.
   Relabelled in eight documents; the pre-registration's "paired-bootstrap
   interquartile range" and criterion R1's "paired-bootstrap interval" are corrected
   by dated note.
3. **Spearman rho: state the aggregation, and state that the falsification is not
   fragile.** +0.8951 (clustered mean per field, knn6 raw) is the frozen value;
   +0.9441 cell-type-centred; +0.9231 median-per-field; +0.7104 per-row. All four
   positive and significant.
4. **Naive biological amplitude: +0.2767**, the section-clustered signed mean under
   `design = base`. 0.3120 is the median |beta|/sd on the same fits and must be named
   when used; 0.2914 and 0.314 are pre-C6.
5. **The 76 % composition surrogate is withdrawn** in `CS_PHASE5.md` sec 4 and
   replaced by the 65.9 % / 85.4 % composition-matched pair.
6. **The two files with live forbidden numbers are marked.**
   `BIO_DELIVERABLE7_CLAIM_AUDIT.md`'s existing marker was scoped to "the suggested
   rewording below" and did **not** cover the *Additionally* paragraph that repeats
   both forbidden figures; its scope is now extended. `CS_PHASE7_C1.md` sec 6 item 3
   now carries a banner recording that the discovery framing is superseded and that
   tiling is no longer the recommended remedy.


| 80 | **`phase8-frozen` TAGGED** | `926439629a07269a32c93f998da0f6e1cd20933c` | Verified before tagging: figure guard 52/52, gene-set gate exit 0 both arms, all three pinned hashes. Pre-registration hashes filled; working tree clean |
| 81 | **Hazard fix committed (`11fa773`) + temp hook retired (`9264396`)** | proven by running `git checkout -- code/` against HEAD | The PI's own commit (`1351ce8`) had **left the guard file out** — my `skip-worktree` protection stopped it being staged. Guards now in history; hazard permanently closed |
| 82 | **A3 fallback screen** | `reports/A3_FALLBACK_SCREEN.md`, `results/a3_fallback/` | **No adequate like-for-like fallback exists** |

## THE A3 SCREEN'S STRATEGIC FINDING: a failure likely ENDS the human arm

**A2 is not the discriminator — both candidates PASS the frozen gate** (strict Tier
A 33/33, all seven modules above the floor, `oxidative_stress` at 36 on both).
That was the obvious screen and it does not separate them.

**What does:** I verified the three panels share a **5,006-gene three-way core**
(H1 5,093; kidney and BM 5,101 each; H1∩kidney **5,008**, H1∩BM **5,023**). The
sender callers are frozen. **So an A3 failure on H1 caused by out-of-band
prevalence on a Prime 5K panel is likely to REPRODUCE on both candidates** — they
are ~98% the same instrument.

Each also adds a cell-budget failure H1 does not have. To clear both A3 floors at
p=5%, a cell type needs **2.4%** of an H1 donor, **29%** of the median kidney
donor, **7.7%** of the median BM donor.

**Bone marrow fails structurally, and I verified it** from
`GSE335963_a3_budget.json`: at p=0.02, **`n_strata_ok: 0`** of 88 donor x type
strata, and **`types_ok_in_all_donors: []` at every prevalence level tested**.

**Consequence for planning: an A3 failure most likely ends the human arm rather
than redirecting it.** That is worth knowing *before* A3 runs, not after.

## ⚠ MY OWN SCREEN WAS WRONG IN BOTH DIRECTIONS

`PHASE7_H1_SCREEN.md` §2 gave **sample** counts where **donor** counts were needed:

| I reported | Actually |
|---|---|
| GSE336890 "9" | 9 Region **slides** carrying **20 patient specimens** (8 AIN / 7 ATI / only **5** reference) — *more* than I said |
| GSE335963 "39" | A **SuperSeries**. Only **6 of 39** GSM are Xenium at all, from **4 donors** (two repeat sections) — far *fewer* than I said |

Neither error changed the H1 selection, but both would have misled anyone using
that table as a fallback list. Corrected in place with a banner.

**Also gone:** §12.1's stated reason for preferring H1 over these two was its age
axis. **Neither candidate deposits age at all**, so that comparison was never
available to make.

**And a data defect worth remembering:** two kidney samples type **all 9,644
features** as `Gene Expression`. Only the ENSG id prefix recovers the real
5,101-gene panel — a title-trusting screen would have mis-sized that panel by 89%.


| 79 | **Commit plan written for the PI** | `reports/COMMIT_PLAN_FOR_PI.md` | 590 changed paths categorised, a 6-step commit order, pre-tag verification commands, and §0 leading with the destructive-checkout hazard. **Nothing executed** — no commits, pushes or tags, per standing instruction |


## ⚠ MORE OF MY OWN OVERSTATEMENTS, CAUGHT BY THE FINAL PASS

The corrections agent ran seven parallel verifiers and several of **my** claims did
not survive contact with the files:

| I said | Actually |
|---|---|
| DCA "roughly **doubles**" the depth loading | **x1.32-1.67** (7239 x1.65, 7259 x1.67, **7352 x1.32**). Verified. My overstatement is now also in `PREREG_PHASE8.md` P29 |
| `NOVELTY_ASSESSMENT` §2.1 **and** §4 O1 were falsified | **Only §2.1 point 3.** O1 asks for a statement the Moran run *endorses*. **My banner overstates its own correction** |
| A7 per-response: -0.0177 / -0.0549 / -0.0337 | Those are the **pre-C6 05:19** file. Frozen 09:06: **-0.0225 / -0.0604 / -0.0307**. Stale in **eight** documents |
| Naive biological amplitude 0.291 SD | Frozen: **0.277** (and 0.036 -> 0.031) |
| "**49** author lines corrected" | **18.** I counted diff lines including context |
| B7 "absorbed 27 arrest genes" | 27 genes, **19** from the Tier A candidate pool |
| Tiled torus "**2.4x** nominal" | **2.35x**; RS_count over all 8 cells is **0.033-0.060**, not 0.040-0.060 |
| CoreScence 69% "verified from git history" | **Impossible** — both scripts are untracked. And 35 *does* exist in the repo: it is the superseded human B7 size |
| H1 called "aging lung" in **5** places | **15 lines across 7 sections** |
| **194** untracked files | **482.** `git status` shows 202 because it collapses directories |

**Conclusions survive throughout; digits did not.** That is the distinction the
ledger now draws explicitly.

## ⚠ A `git checkout` IS CURRENTLY DESTRUCTIVE

Verified: `git show HEAD:code/build_genesets.py` contains **none** of the guards
added today. Restoring `code/` from HEAD would reinstate the version that
**silently overwrites `genesets/*.txt` with EMPTY Tier B modules**. And
`results/phase3_pre_c6/` — the **sole copy** of the baseline this entire ledger
compares against — is **98 files, 0 tracked**.

**Committing is not housekeeping. Until it happens, the repo can destroy its own
evidence base with one ordinary git command.**


## The Moran's I run overturned the project's planned differentiation — and replaced it with something better

**Falsified as written — and I have now VERIFIED this directly** from
`results/moran/moran_vs_a7.csv` (12 fields), not merely accepted it as reported:
**rho = +0.923 raw, +0.944 cell-type-centred**. The centred figure matches the
agent's exactly; my raw figure differs slightly from its +0.895, almost certainly
an aggregation choice (medians across sections at knn6 weights) — immaterial to
the conclusion, and recorded rather than smoothed. `tnfa_nfkb_proximal` ranks top
on both statistics and `neg_probe_rate` bottom on both. `NOVELTY_ASSESSMENT.md` §2.1 point 3 and §4 O1 both told the
project to write *"the two tests disagree"* — **that would have been visible to
any reviewer in a single plot.**

**The replacement is stronger because it is measured.** I verified the bound from
`results/moran/moran_kernel_power.csv`:

- the entire A7 gradient contributes **Delta-I = 2.2e-4 — 0.83% of the observed
  control Moran's I**
- the smallest kernel amplitude Moran's I can resolve is **0.362 SD**
- against the A7 control gradient of **0.074 SD** (5x smaller) and the project's
  own naive biological amplitude of **0.291 SD**

**0.362 > 0.291: Moran's I could not have detected the paper's headline effect
either.** So "different question" survives — justified by **power, not
orthogonality**. That is a quantified argument replacing an assertion that was
about to be wrong.

**A second, genuinely new point against the prior art.** Voyager reads near-zero
Moran's I on controls as "no technical artifact spatial trend". But **genes
matched to the controls on total counts give identical values** (-0.00018 vs
-0.00012). A negative control probe carries ~21 counts per section; the median
gene carries 5,885. **Voyager's controls-vs-genes contrast is an abundance
contrast**, and the per-feature statistic has no power at control abundance. Holds
under CP10K+log1p too.

**A7's split reproduces from an independent statistic** with no kernel, lambda,
sender call or nuisance design: `all_controls` +0.0455, codewords +0.0421,
**probes +0.0058 (7x smaller)** — the same ordering. **Against interest:** Moran's
I does *not* call the probes flat (pooled CI excludes zero), and it ranks
`genomic_control` lowest where A7 ranks it third.


## THE MOST IMPORTANT FINDING OF THE AUTONOMOUS WINDOW

**A third, independent line of evidence for the N2-vs-N5 inversion was sitting in
the project's own headline synthetic figure the entire time.**

From `figures/figure1_data.csv` — committed, unchanged, **planted ground truth**.
I recomputed all of it:

| Where ell/lambda >= 2 (n=8 cells) | mean CI coverage |
|---|---|
| naive | 0.51 |
| **matched-decoy** | **0.35 — WORSE THAN DOING NOTHING** |
| nuisance-conditioned | **0.85** |

|relative bias| on lambda: matched-decoy beats naive in only **12 of 20** grid
cells and is **worse in 8**. Worst case: decoy **2.27**, naive **2.02**,
nuisance **0.33**.

**Why this is the strongest of the three lines:**

| Evidence | Matched decoys | Covariates | Ground truth? |
|---|---|---|---|
| A7 (technical gradient) | leaves ~80% intact | N5 removes it | no |
| Composition-matched | removes 1.6% | removes 85.4% | no |
| **Figure 1 (synthetic)** | **coverage 0.35, worse than naive** | **0.85** | **YES — planted** |

The first two are open to "your covariates removed real signal." **The third is
not**: the kernel is planted, so the truth is known. And it shows matched-decoy
is not merely *insufficient* — it is **actively worse than doing nothing** on
coverage, and worse than naive on bias in 8 of 20 regimes.

**The figure needs no regeneration.** Only its §25 caption was wrong, and that is
now corrected. The result has been in the repo since Phase 1.


## The nuance the agent kept, which I would have been tempted to drop

**The denoised score is simultaneously *better* by DeepScence's own internal
criterion (0.126 -> 0.470) and *more* depth-confounded.** The tool's own quality
metric and the confound axis point in **opposite directions**.

That is a far more interesting finding than "denoising is bad", and it is the
version a reviewer should see: a published tool, run exactly as published,
improves on its own stated objective while getting worse on the axis the science
depends on. The agent flagged it as "the uncomfortable half you'd want a reviewer
to see" rather than reporting only the half that favours the decision.

**Seed instability verified by me** from `d2_stability.csv`: `top5_jaccard`
**exactly 0.0000** on two seed pairs at r = 0.5703 / 0.5732, with
`anchor_sign_flipped = False` throughout — so it is genuine instability in *which
cells are called*, not a polarity flip. **And no diagnostic fires.** The
determinism control on the same file shows r = 0.999999, so the base is
re-derivable.


## §4 (D-b)'s PREMISE IS REFUTED BY MEASUREMENT

§4 (D-b) states DCA *"is precisely the step that would normalize depth — the
confound under investigation."* **Measured, it roughly doubles the depth loading.**
Verified by me from `results/phase8_d2/d2_depth.csv`:

| section | denoise=False | denoise=True | delta |
|---|---|---|---|
| 7239 | 0.389 | **0.640** | +0.251 |
| 7259 | 0.318 | **0.531** | +0.214 |
| 7352 | 0.410 | **0.542** | +0.132 |

Three of three sections, both arms. The top-5% call set becomes **100% hepatocyte**
on all three; sender-set agreement between configurations is Jaccard **0.12-0.28**.

**DECISION TAKEN (delegated authority): freeze `denoise=False` as PRIMARY — now a
CHOSEN value — with `denoise=True` as the published-default sensitivity.** This
converts a **caveat** ("we deviated because DCA would not install") into a
**finding** ("we ran the published default and it is more depth-loaded and
seed-unstable, so we chose against it"). Strictly the stronger position.

## Three further results from D2

1. **The published default is not reliably reproducible.** One of three seeds
   returned a top-5% sender set **perfectly disjoint** from the other two
   (r = 0.57, Jaccard **0.000**), against Jaccard 0.76 for `denoise=False` across
   seeds. A second concrete instability in the same tool, beside the D3 polarity
   flip.
2. **§6's own named estimator was too weak to test its own question.**
   Median-of-ratios *poscounts* removed only **11-24%** of log-depth variance and
   returned a clean null on four sections — a null about a weak normalisation. A
   `lib` configuration removing 100% cuts depth loading by **74% and 93%**.
3. **The frozen base is RE-DERIVABLE, not merely intact.** Same seed reproduces
   the committed scores at **r = 0.9999991 / 0.9999999** — including
   `deepscence_sbr.csv`, written 2026-08-20 **before two container rebuilds**.
   Given that 194 untracked files currently hold the evidence base, a demonstrated
   re-derivation carries real weight.

**Resource fact worth keeping:** DeepScence holds **five** dense
`n_cells x 4,845` float32 arrays at once, so **~16 GB of cgroup headroom per
83k-cell section** is the working requirement. Six D2 jobs were OOM-killed before
that was established.


## ⚠ A PHRASE I USED ALL SESSION WAS WRONG: "circular Tier A"

I repeatedly told the PI the independence claim "rested on a **circular sender
definition**". That misplaces the defect.

**The defect was in the RESPONSE module.** Pre-C6 B7 `secondary_senescence` had
absorbed **27 arrest genes**. Because the §11 gate removes A∩B overlaps *from Tier
A*, a B7 stuffed with arrest genes stripped the sender set — leaving Tier A a
**hollow 25-gene remnant by subtraction**. Verified: B7 **38 -> 108** and Tier A
**25 -> 33** when B7 was re-sourced.

This matters for what the paper says went wrong: not "our sender set was circular"
but "**a response module had absorbed the sender's genes, and the disjointness
gate resolved that by gutting the sender**." Different diagnosis, different lesson.

**Also disclosed, against interest:** C6 moved one circularity measure the *other*
way — CoreScence circularity **79% -> 88%**, and Tier A ∩ CoreScence **2/25 ->
4/33** (adding `Cdkn2b`, `Mdm2`), which touches the arrest-vs-DeepScence headline
pair.

## The restatement's second plank is REFUTED — and the agent said so

I told the PI the drafted restatement in `CS_PHASE8_CALLERS.md` §3 could be used
"verbatim". **Part of it cannot.**

- **Plank 1 (one pair below chance) — replaced, from evidence already on disk.**
  **SenePy vs DeepScence: pooled 0.737, z = -15.08, below chance in 10 of 11
  sections.** I verified it is **byte-identical pre- and post-C6** because it does
  not involve Tier A — which is exactly what makes it robust. Strengthens to 0.495
  under the caller-free re-anchor. Caveat travels with it: under |score| ranking
  the same pair is at chance (1.025, n.s.), so the anti-concordance is about
  **polarity**.
- **Plank 2 ("the direction of each pair is predicted by its depth loading") —
  REFUTED by three tests.** Pair-level exact permutation over the five independent
  pairs **p = 0.30**; within-pair across 11 sections the depth-loading product is
  *negatively* rank-correlated with agreement in **5 of 5** pairs (rho -0.16 to
  -0.70); pooled continuous version null (rho = +0.096, p = 0.49). Broken by the
  largest value in the matrix — arrest vs `Cdkn1a`+ at **1.471**, above chance in
  11/11, whose cause is **biological** (four of the eight re-entering genes sit on
  the p53 axis that induces `Cdkn1a`).
- **The weaker claim the data actually supports:** the dependence is weak,
  heterogeneous, and mechanistically *different at each end* (0.74x to 1.47x) —
  still incompatible with one latent state. The two-depth-camp fact survives at
  11/11 for all four callers.

**The agent was asked not to manufacture a replacement plank if one could not be
found honestly. It found one for the first and refused for the second.**


## Two more corrections — one to me, one an agent made to itself

**My "all three `fig_phase3_*` are on a stale basis" was wrong about two of three.**
The agent checked the *producers* rather than the filenames:

| figure | gene-set dependent? | outcome |
|---|---|---|
| `fig_phase3_composition` | no — cell types only | byte-identical: a **genuine** pass |
| `fig_phase3_caller_depth` | **yes** | re-pointed to the frozen table, **CHANGED** |
| `fig_phase3_tierC_identifiability` | **no, verified** | byte-identical: genuine |

**The Tier C case nearly fooled it:** `tierC_lr.py` reads `sender_flag_p95` at
line 48 — a column that *did* change — but **never uses it**. Dead code. Only
re-running the producer and diffing showed that the read was inert. Reasoning
from "it reads a changed input" would have given the wrong answer.

**The agent also audited the auditor and won.** Fact-check item **R6**'s numbers
(3.0-6.7% / 1.2-13.3%) were correct for the **pre-C6 A7 file of 05:19**, which
the re-run superseded at 06:52. On the frozen 825 control fits the filter admits
**4.8% on the full design — essentially nominal, identical across all five
families**. The substance held (the "2-3x nominal" bound is on the *estimator*,
not the *filter*, and that sentence is withdrawn), but the audit's figures were
stale and it said so.

**And it corrected its own misattribution unprompted**, with a clean statement of
the reasoning failure: timestamps told it *that* something wrote the files, never
*who* — it inferred from the one job it could see and did not check.


## ⚠ ANOTHER OF MY OWN ERRORS, CAUGHT BY AN AGENT

**My figure regeneration ledger exempted figures 2a and 2d as "not null-dependent,
verified byte-identical". The exemption was invalid.** Both statements are true
*of the null correction* and false of a **gene-set** change. Figure 2a is binned
response vs distance-to-nearest-sender; 2d is median distance vs sender density.
**Both are functions of the sender set, and the sender set changed 25 -> 33.**
Both were regenerated and both moved.

**The subtler trap the agent caught:** `figure2a` **caches its own input**
(`figure2a_stratified_curves.csv`, rebuilt only if absent). Left in place, 2a
would have regenerated from **pre-C6 senders** and returned byte-identical —
**reading as a passing reproducibility check when it was a stale-input result.**
The cache is now deleted before the figure stage.

**Pinned-file hashes updated after the legitimate 8.7 supersession:**
`sf_summary.csv` -> `a5ccc9b0e81f4c335e8039e975ec1975`;
`summary_phase3.txt` -> `dc92ddc6605eef52f6359aeab4e16fd7`;
`perm_nulls.csv` -> `d906394958dbe1b99981756290c511fa` (unchanged since earlier).


## THE SECOND-STRONGEST RESULT: matched decoys fail where covariates work — now shown TWICE, independently

**The §15 protocol as written is inert.** Composition-matched decoys remove
**1.6%** of the pooled naive amplitude (3.5% within cell type). Verified by me
from `compmatch_reruns.csv`: five seeds (20260901-05) give median SF
**0.9837-0.9839**, seed sd **1.2e-4**. The seeds buy nothing — at a 0.99987 match
rate the greedy matcher has no freedom left.

**The like-for-like covariate counterpart — same variables, same fits — removes
85.4%.** Cell-type intercepts alone give **65.9%**, reproducing the published 66%
almost exactly. Never below 0.52 in any section.

**A factor of fifty from the same variables.** §15's frozen parameter, as
written, would certify as "not composition" a gradient that is 66-85% composition.

**Why this matters beyond the fix — it CONVERGES WITH A7:**

| Analysis | Matched decoys | Covariates |
|---|---|---|
| **A7** (technical gradient) | leaves ~80% intact | N5 removes it |
| **Composition** (this) | removes 1.6% | removes 85.4% |

**Two independent analyses, same conclusion: matched-decoy designs systematically
fail to remove what covariate adjustment removes.** An external novelty review had
already rated the N2-vs-N5 result the project's strongest contribution; this
doubles its evidentiary base and makes it a general claim rather than one
observation.

**DECISION TAKEN (delegated authority): freeze BOTH variants.** Covariate-adjusted
as primary, matched-decoy reported alongside — its inertness is itself a finding.
Additive, not a substitution; same report-both discipline used for the two Tier A
variants and for B7. Freezing the matched-decoy version alone would license a
claim the data contradicts, which is the opposite of what a pre-registration is
for. Routed to the pre-registration agent with the nine D15.x ambiguities to
carry into the deviation table.


## The strongest result of the autonomous window

**Contribution 3 stands, and is now defensible on calibration rather than assertion.**

The variance correction **agrees** with N3-tile on the data: N3-var **0.996**,
N3-tile 0.971, window-matched 0.995; beta-hat and lambda-hat match `perm_nulls_c1.csv`
to 1e-16. So the corrected answer does not move.

**But the calibration simulation found that C1 replaced a liberal test with a
MORE liberal one.** Type-I error at nominal 5%, irregular window (independently
verified by me from `var_sim_calibration.csv`):

| variant | type-I error | verdict |
|---|---|---|
| **tiled torus (4x4)** | 0.048 - **0.118** | up to **2.4x nominal** |
| whole-window torus | 0.033 - 0.073 | liberal |
| **RS_count (variance-corrected)** | **0.040 - 0.055** | **calibrated** |
| drop without standardising | 0.003 - 0.020 | over-conservative |

**The standardisation, not the dropping, is what works** — which is exactly what
Mrkvicka §2.1.4 predicts and what the project's tiling missed.

**Recommendation adopted: N3-var / N4-var become primary.** Not because the
number moved (it did not) but because **N3-tile is the one variant with a
published theoretical prediction against it**, now confirmed empirically.

**Reported against interest, twice:** on real data the direction is the
*opposite* (N3-tile 0.971 < bbox 0.999; rejection 0.801 < 0.824) and the report
explains why (tile seams **81** lambda-hat apart at the sourced pooled lambda-hat of
14.7 um; the "76" this entry recorded rested on the withdrawn 15.7 um) rather than
reconciling it away; and
RS_count's own 1/n_i variance assumption was validated on this data (slope
-0.451 / -0.492 against a predicted -0.5).

**Caveat recorded:** the reportable population is now **153, not 160** — the M1
re-run regenerated `main_fits.csv` — so C1 variants were recomputed on the
current population and are **not cell-by-cell comparable** with `CS_PHASE7_C1.md`.
The covariate-adjusted `*_full_sf` for var is at 200 perms, not 1,000.


## ⚠ A NUMBER I RELAYED REPEATEDLY WAS WRONG — CORRECTED

I reported "CoreScence circularity rose **69% -> 88%**, a real cost of C6" several
times, and used it to support PI decision D6. **The 69% was fabricated.**

- `24/35 = 69%` was a **typed-in literal** in two scripts. The project's own
  committed `results/phase3/n8_disjointness_*.csv` says
  **`corescence_on_panel = 33`**; `logs/ds_smoke.log` says **31** under the strict
  ortholog map. **No convention yields 35.** I verified both independently.
- Re-derived from files: **26/33 = 79%** pre-C6, **29/33 = 88%** under C6.
- **Two errors compounded.** First, a fabricated baseline. Second, the
  "69 -> 76 -> 88" story compared a *mouse* baseline against two *human*
  configurations — apples to oranges.

**Corrected, measured within each arm:** the cost of C6 is **+9 points (mouse,
79->88)** and **+12 (human)**. Roughly **half** what I attributed to it. Both arms
land on 29/33.

**Direction of the correction cuts both ways, and both halves matter:** the
*published mouse arm is MORE circular than reported* (79%, not 69%) — against
interest — while the *cost of C6* is smaller than claimed. **PI decision D6
(strip-and-refit primary) is unaffected**: 88% is still high enough to justify it.
`figure_gs3` now shows four bars so the within-arm cost is legible.


## ⚠ TWO THINGS NEEDING THE PI ON RETURN

**1. The entire Phase 8 evidence base is UNTRACKED.** `genesets/human/` and
`genesets/mouse_c6/` are untracked; only 266 of 1,238 files under `results/` are
tracked. Every file backing C1, A7 and the gene-set freeze exists **only in the
working tree**. A container wipe — which has happened **twice today** — would not
destroy them (they are on the network volume) but nothing is version-pinned, so
the pre-registration cannot honestly reference immutable artefacts. **This must
be committed before the `phase8-frozen` tag means anything.** I have not
committed anything, per standing instruction.

**3. OUTSTANDING (minor): `reports/CS_PHASE7_C1.md` carries no attribution for
the torus null at all** — it neither cites Lotwick & Silverman nor wrongly cites
CellWHISPER. Not a false claim, but a gap; the incoming variance-correction
report should supersede its framing.

**2. A manuscript-bound error was caught with two days to the deadline.** The
submission patch stated a "-0.070 SD gradient **in negative-control probes**".
That is the pooled `all_controls` response. The 40 negative-control probes —
which `PREREG_PHASE8_genesets.md` §11 designates the **primary** A7 null — are
**flat** (-0.0177, p=0.183). On the pre-registered response, M1's A7 **passes
naively**. Corrected in `SUBMISSION_PATCH_2026-08-29.md` §4a with the full
per-response breakdown. The "assay is not flat" conclusion survives on codewords
(-0.0549, p=0.039) and genomic controls (-0.0337, p=0.0039); the *name* on it did
not. The N2-vs-N5 result is unaffected and remains the strongest finding.


## Errors caught before they propagated

| What | How it was caught |
|---|---|
| **Submission patch went stale within 15 min** | Caller numbers changed at 06:09 under the new Tier A. Patch marked DO-NOT-APPLY before the PI could hand-apply it to a manuscript due Aug 29 |
| **`caller_coverage_gate.csv` mixes bases** | Pre-C6 2-section rows pooled against post-C6 11-section rows. Routed to the M1 agent |
| **Composition-matched protocol has NO implementation** | Listed in §15 as a *frozen parameter*. Now being implemented (D15) |
| **Window is a 100um literal, not the 99th percentile §15 claims** | Frozen as-is with provenance recorded (D16) |
| **B6 one-gene-trim claim was false** | 31 vs floor 30 — a one-gene trim passes. Takes two |
| **My own venv fix was wrong** | `--system-site-packages` meant nothing was installed into it; every import still resolved to the container. Rebuilt isolated |
| **Statusline shipped depending on `jq`, which is absent** | Directory rendered empty, script exited 0, nothing would have warned |
| **Deviation-ID namespace collision (D1-D14 mean two different things)** | Found across 13 files / 142 references. Renumbering DEFERRED (documents still moving); disambiguation note added to the status board so nothing is miscited meanwhile |
| **I quoted a memory ceiling 4.3x too high all session** | The cgroup limit is **57.7 GB**, not the 251 GB `free -g` reports (that is HOST memory; the container is capped). Every "190 GB available / no OOM risk" statement I made was wrong. **This is what actually caused the DeepScence OOM kills** — I ran 5 sections concurrently against a budget a quarter the size I believed. The D2 agent's own logging (`cgroup free 25 GB`) had the right number before I did. Recurring job now reads `/sys/fs/cgroup/memory.current` against `memory.max` |
| **My own monitoring was silently broken all session** | `find` here is **`bfs`, not GNU findutils**, and it REJECTS `-newermt "-60 minutes"` — erroring out and yielding nothing. Every "0 files written" I reported was a **false negative**; the real count was 138 files in the last hour alone. Caused one false stall alarm. Fixed by switching to a reference-file method, and the recurring job now carries the warning so it cannot recur |
| **An agent invented a figure, then caught itself** | A pre-registration draft carried an invented "inner ~55-85%" bracket; the agent replaced it with measured values from `results/phase3/window.csv` (7.4-21.5% beyond 80um) before reporting |

---
*Rows are appended as work lands. In-flight items are NOT listed here.*

## A DEFECT INSIDE THE FROZEN PRE-REGISTRATION

**There is no bootstrap in the reported bracket at all.** `summarize_phase3.py:99`
is a plain `np.quantile([.25,.5,.75])` over per-fit point estimates. So
`PREREG_PHASE8.md` §5's *"paired-bootstrap interquartile range"* and replication
criterion **R1**'s *"paired-bootstrap interval"* are **misnomers in the frozen
document**.

The pre-registered 400x100-block bootstrap emits **per-fit** CIs only (median span
[-0.415, +0.381]); there is **no interval on the median across fits**. Handled as
a dated correction that leaves R1's threshold and M1's outcome untouched — the
right treatment for a frozen artefact. **A genuine CI would require a new run**,
and the agent stopped and said so rather than manufacturing one.

**Consequence: every `0.088 [-0.017, 0.234]`-style bracket in this project is an
IQR across fits, not a confidence interval.** Relabelled in eight documents.

