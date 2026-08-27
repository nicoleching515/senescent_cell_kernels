# CS Phase 8 / task 8.7 — the M1 end-to-end re-run

**Status: complete. Verdict: the contribution stands, and the two things that
moved both moved against us and are reported first.**

Four changes were carried together for the first time: the corrected N3/N4
in-tissue nulls (C1), the promoted C6 mouse gene sets (D5), both pre-registered
Tier A variants on the sender axis, and DeepScence at 11/11 section coverage
(C7/D1). `reports/CORRECTIONS.md` (task 8.8) is the ledger of what moved and
why; this report is the run record — what was executed, on what, with what
verification.

---

## 0. Lead with the bad news

**1. Caller agreement rose again, and the last pair that was below chance no
longer is.** Task 8.4 reported the depth- and type-matched agreement rising from
1.030× chance to 1.118× at full coverage. Under the promoted strict-33 Tier A it
is **1.212× (z = 21.9, p = 1.8 × 10⁻¹⁰⁶)** on the four-pair basis, and — the
part that hurts — **on the *published two-section base* it is already 1.13× at
p = 4 × 10⁻⁸**. The independence claim does not survive even without the
coverage fix. `CS_PHASE8_CALLERS.md`'s line *"Tier A vs SenePy is the one pair
that does not move … two callers are genuinely close to disjoint. That is worth
keeping"* must be withdrawn: post-C6 that pair is **0.972, z = −1.63, p = 0.10**,
above chance in 4 of 11 sections. There is no longer any pair reliably below
chance.

**2. The module whose name promised a paracrine effect lost its gradient when it
was made actually paracrine.** The re-sourced B7 `secondary_senescence` goes from
36 to **22** reportable fits and its naive amplitude from 0.342 to **0.246**. B7
is the only module whose response score changed, and it accounts for the entire
net loss of reportable fits (160 → 153). The pre-C6 B7 shared 27 arrest genes
with the sender-adjacent literature; the gradient it showed was partly a gradient
in arrest.

**3. Two published/earlier numbers are wrong and are corrected here.**
`CS_PHASE8_CALLERS.md` §2.1's "22 of 33 sections above chance" is **20 of 33**.
`caller_coverage_gate.csv` was pooling a pre-C6 two-section base against a
post-C6 eleven-section base; it now carries four explicitly-labelled bases.

**4. Two of the three pinned files were legitimately superseded and say so.**
`perm_nulls.csv` `3b77aa1b…` → `d9063949…`, and `sf_summary.csv` /
`summary_phase3.txt` with it. The reason and the two recovery paths are in
`CORRECTIONS.md` §6. The two DeepScence files are byte-identical.

**And the headline that did not move.** Re-verified after every stage landed,
from the completed tree (`results/phase3/m1_final_audit.txt`):

| | pre-C6 | **post-C6, final** |
|---|---|---|
| reportable fits | 160 | **153** |
| naive amplitude (response-sd) | 0.326 | **0.329** |
| **controlled amplitude, N2+N5+N6** | **0.027** [−0.028, 0.090] | **0.029** [−0.007, 0.084] |
| **SF, N2+N5+N6** | **0.082** [−0.099, 0.249] | **0.088** [−0.017, 0.234] |
| **detectable bound, 80 % power** | **0.203** | **0.183** |
| controlled fits positive with CI excluding 0 | 15 / 160 | 13 / 153 |

The controlled amplitude remains **far below** the bound, so §18 outcome **A**
stands, and the C1 verdict — the in-tissue correction does not move the
answer — stands for **both** pre-registered Tier A definitions, and now also
under the variance-corrected null that replaces N3-tile as primary
(N3-var 0.996). None of these three moved between the partial state reported
mid-run and the completed tree.

**One thing did change after this report was first written, and it is a
correction to correction C1 itself.** N3-tile is no longer the primary corrected
null: a calibration simulation shows tiling is *more* liberal than the
whole-section torus it replaced (up to 2.4× nominal on an irregular window),
so the field-standard variance correction N3-var / N4-var is presented instead.
The measured answer is unchanged; the defensibility is not. See
`CORRECTIONS.md` §8.1.

---

## 1. What was run, and what was deliberately not

Per `CS_PHASE8_C1_CLOSEOUT.md` §5.2's stage table. `python3` is the system
interpreter throughout; the persistent env at `/workspace/envs/sasp311` was not
used, per the coordinator's instruction not to switch mid-run.

| stage | script / invocation | re-run? | why |
|---|---|---|---|
| Xenium ingest | `prepare_samples.py` | **no** | raw data unchanged |
| cell-type annotation | `annotate_pipeline.py` | **no** | markers independent of Tier A/B |
| DeepScence | `run_deepscence_all.py` | **no** | D1 finished at 11/11; it is an input |
| Phase 4 CCC battery | `phase4_run.py` and friends | **no** | depends on Tier C pairs and cell types only |
| anatomy / senders / module scores | `phase2_downstream.py`, 11 sections | **yes** | gene sets changed |
| Phase 3 cache | `sasp_phase3.prep(force=True)`, 11 sections | **yes** | new sender flags and module scores |
| window | `--stage window --sections all` | yes | new sender sets |
| main fits | `--stage main --sections all --calls all9 --n-jobs 24` | **yes** | 297 jobs (11 × 9 calls, per-module fanned out over 7 modules) |
| published bbox nulls | `--stage perm --sections inband --calls tierA_p95 --n-perm 1000` | **yes** | needed so `sf_summary.csv` is single-basis (`CORRECTIONS.md` §6) |
| corrected nulls, primary | `--stage perm_c1 --sections inband --calls tierA_p95 --n-perm 1000` | **yes** | C1 |
| corrected nulls, N7 | `--stage perm_c1 --calls <5 non-primary> --n-perm 1000` | **yes** | **1,000 perms, not 200** — §24.3, per the PI's decision this pass |
| bbox nulls, N7 | `--stage perm --calls <5 non-primary> --n-perm 1000` | yes | " |
| second Tier A variant | `--stage perm_c1 / perm --calls tierApm_p95 --n-perm 1000 --tag _pm` | **yes** | pre-registered sensitivity |
| curves | `--stage curves` | yes | |
| N8, stratification, attribution, combined, Poisson, λ-scale, Ripley, correlogram | the eight small scripts | yes | |
| destructiveness diagnostic | `phase3_null_diag.py` | yes | |
| A7 negative-control-probe kernel | `run_a7_control_probes.py`, `summarize_a7.py` | **yes** | it uses the sender call, which changed |
| caller tables + gate + D3 | `caller_disagree_all.py --all`, `summarize_caller_coverage.py`, `caller_disagree_d3.py` | **yes** | they read `tierA_score`, which changed |
| Phase 5 kernels and superposition | `run_phase5_*`, `_se_ratio_phase5.py`, `_spline_window_check.py` | **yes** | |
| summaries | `summarize_phase3{,_c1}.py`, `summarize_phase5.py`, `summarize_super_callers.py` | yes | |
| figures | one pass — see §5 | yes | |
| `code/make_figure2.py` | — | **NOT RUN** | superseded producer of figure 2a; it now refuses to run |

**Permutation counts are stated, not assumed.** Every perturbation-null file
written by this re-run carries `n_perm = 1000`: `perm_nulls.csv`,
`perm_nulls_c1.csv`, `perm_nulls_n7.csv`, `perm_nulls_c1_n7.csv`,
`perm_nulls_pm.csv`, `perm_nulls_c1_pm.csv`. The pre-C6 N7 files were at 200;
they are not carried forward.

---

## 2. Code changes made for this run

All additive; nothing that produces an existing number was altered except where
noted.

| file | change |
|---|---|
| `phase2_downstream.py` | scores the seven `A_sender_for_<module>.txt` sets and writes `tierA_<module>_score` + `sender_flag_<module>_p{90,95,99}`. Existing columns and their order are untouched, and the block sits **before** the DeepScence merge so row alignment cannot drift. Asserts the per-module set names match the Tier B module list |
| `sasp_phase3.py` | `prep()` caches the per-module scores and flags; `Sec.sender_mask(call, module=None)` gains the `tierApm_p*` family and raises rather than guesses if the cache predates Phase 8; `Sec.has_permodule` |
| `run_phase3_nulls.py` | `TIERA_PM_CALLS`, `ALL9_CALLS`, `is_permodule()`, `_js()`, `_expand()`; `SectionFit(..., module=)`; `_section_job` / `_perm_job` / `_perm_c1_job` take a module and restrict to it; `--calls all9|tierA_pm`; `--tag` for output suffixes; `fit_cell` and both perm rows record `sender_set`. **Seeds for the six existing calls are unchanged** — `_expand` reproduces `MASTER_SEED + step_i·i + step_j·j` for non-per-module calls at the same section and call ordering, so old and new are directly comparable |
| `summarize_phase3.py` | §6's N7 table iterates both Tier A variants |
| `summarize_caller_coverage.py` | four explicitly-labelled bases; every row carries `tierA_definition` and `tierA_n_genes` (see `CORRECTIONS.md` §5.1) |
| `crossarm_geneset_table.py` | reads the pre-C6 mouse sets from `git show pre-c6-genesets:` instead of `genesets/*.txt`, which the D5 promotion turned into the post-C6 state (§4 below) |
| `make_figure2bc.py` | writes `figure2b_data.csv` and `figure2d_data.csv`; those two panels had no data CSV |
| `make_phase5_figs.py` | writes `figure3_data.csv` |

New helper scripts, both read-only over results:
`code/m1_headlines.py` (the §17 headline vector from one results tree) and
`code/m1_compare_modules.py` (per-module / per-cell-type / per-call pre-vs-post).

**No packages were installed.** The environment rebuilt from `requirements.txt`
was sufficient for every stage: numpy 2.4.6, pandas 2.3.3, scipy 1.17.1,
scikit-learn 1.9.0, anndata 0.12.19, scanpy 1.11.5, joblib 1.5.3, h5py 3.16.0,
matplotlib 3.11.1, senepy, on Python 3.11.10.

---

## 3. The control that makes this an attribution rather than a guess

`phase2_downstream.py` re-run on all 11 sections, compared cell-by-cell against
the pre-C6 outputs:

* `cdkn1a_counts`, `cdkn1a_pos`, `senepy_score`, `zonation_score`,
  `compartment_label`, `dist_to_boundary_um`, `dist_to_portal_triad_um` and
  **six of the seven Tier B module scores** come back **bit-identical**
  (max \|Δ\| = 0) across a container rebuild;
* only `secondary_senescence` (max \|Δ\| 0.43) and `tierA_score` (0.244) moved.

So the environment is not a confound and exactly two inputs changed. Every
number in `CORRECTIONS.md` is attributed on that basis.

The disjointness gate was **re-verified independently** against the panel derived
from the data (`cell_feature_matrix.h5`: 5,106 `Gene Expression` features minus
9 genotyping probes = 5,097): Tier A strict 33/33 on panel, all seven Tier B
≥ 30 on panel (`oxidative_stress` still on a margin of exactly 1: 31), A ∩ ∪B = 0,
every per-module set ≥ 15 and disjoint from its own module. **GATE PASS.**
Reported against interest: the Tier B modules are **not** mutually disjoint —
18 of 21 module pairs share genes. No gate requires it, but no text may claim it.

---

## 4. Two silent breakages the promotion caused, found and fixed

**4.1 `crossarm_geneset_table.py` would have mislabelled its own pre/post
columns.** It read the "mouse pre-C6" column from `genesets/*.txt` and the
"mouse C6" column from `genesets/mouse_c6/*.txt`. After the D5 promotion those
are the *same files*, so `mouse_pre` would have silently reported the C6 values
(Tier A 33 instead of 25; B7 108 instead of 38) and the cross-arm symmetry figure
would have shown no C6 effect at all. It now reads the pre-C6 sets from
`git show pre-c6-genesets:genesets/<name>.txt` and fails loudly if the tag is
unreadable. Re-run: `mouse_pre` 25 / 38, `mouse_C6` 33 / 108, as it should be.
**All four gene-set figures (`figure_gs1`–`gs4`) then reproduce byte-identically**,
which confirms no published gene-set number moved and that the earlier run of
that script happened before the promotion.

**4.2 `figure2a` caches its own input.** `make_phase5_figs.py --which 2a` rebuilds
`figures/figure2a_stratified_curves.csv` only if the file is absent. The sender
set changed, so leaving the Aug-20 cache in place would have regenerated
figure 2a from **pre-C6 senders** and reproduced it byte-identically — looking
like a passing reproducibility check. The cache is deleted before the figure
stage runs. This is the same failure mode as `make_figure2.py`: an output that
comes back identical because the *input* was stale, with nothing warning.

**Related, and stated because the regeneration ledger says otherwise.**
`figures/revised_candidates/README.md` exempts figures 2a and 2d as "not
null-dependent; verified byte-identical". That verification was against the
**null correction**, not against the **gene sets**. Both panels are
sender-dependent — 2a is the binned response versus distance to nearest sender,
2d is median distance to nearest sender versus sender density over
section × sender-definition combinations — so both are expected to move under a
new Tier A. §5 records what actually happened to them.

---

## 5. The figures — one pass, from the frozen configuration

`python3 code/check_figures_guard.py` **before**: `OK: all 27 committed figures
match`, exit 0. **After**: exit 1, listing exactly the intended changes and
nothing else.

| figure | ledger said | what happened | new PNG md5 |
|---|---|---|---|
| Figure 1 | no | **unchanged**, as §19 says | `1eaa67c094eb2916f3f041dbf4c94566` |
| **Figure 2a** | *exempt* | **CHANGED** — the exemption was verified against the null correction, not the gene sets, and 2a is sender-dependent | `332fbca7b7e8eb69b08abb5beab184dd` |
| Figure 2b | yes | CHANGED — corrected N3/N4 bands, new gene sets | `fbe42009f8d73d96af6921c4e404369e` |
| Figure 2c | yes | CHANGED — corrected variants, new gene sets | `9519bcef8f2f749790f35786129f0666` |
| **Figure 2d** | *exempt* | **CHANGED** — same reason as 2a; it plots median distance to nearest sender against sender density | `a63bb8ff710699dc1dfb36f456e82c1a` |
| Figure 2e (panels e–h) | already new | regenerated, A7 panel h filled from the C6 re-run | `b1e791607f0d7c82fa7f6084cacdb21e` |
| Figure 3 | yes | CHANGED — mouse-only; the two-arm version is Phase 10 | `668225a45adf56a9486ff013c0e70109` |
| Figure 4 | yes | CHANGED — de-hardcoded source plus new gene sets | `d44fac63411d6c30a42c40894a287f17` |
| `figure4_supp_*` | — | unchanged (Phase 4 not re-run) | — |
| `fig_phase3_*` (3) | not listed | **unchanged** — see the caveat below | — |
| `figure_phase8_callers` | — | regenerated for the four-basis gate | `a77a3aaf8c10037d241f914a44e8fe01` |
| `figure_phase8_d3` | — | regenerated | `205191aab7d42791f1035aebb19f5a52` |
| `figure_gs1`, `gs4` | — | **byte-identical** | `9d5b5b4a…`, `44e4ca2d…` |

**Two figures in the ledger's exempt list did change, and the ledger's reasoning
is why.** It exempts 2a and 2d as "not null-dependent; verified byte-identical".
Both statements are true *about the null correction* and false about a gene-set
change: 2a is the binned response versus distance to nearest sender and 2d is
median distance to nearest sender versus sender density, so both move whenever
the sender set moves. The ledger has been amended.

**Every regenerated figure now has a `*_data.csv`.** Figures 2b, 2d and 3 had
none; producers were extended to write `results/phase3/figure2b_data.csv`,
`results/phase3/figure2d_data.csv` and `figures/figure3_data.csv`. PNGs were
compared, not PDFs — matplotlib date-stamps PDFs, and the guard strips
`CreationDate`/`ModDate` before hashing them.

**Reported against interest — three caveats on the figure pass.**

1. **`fig_phase3_caller_depth`, `fig_phase3_composition` and
   `fig_phase3_tierC_identifiability` are unchanged because they read the
   *committed two-section* caller tables**, which are pre-C6. They regenerate
   byte-identically, which looks like a passing reproducibility check and is
   really a stale-input result. They are not in the regeneration ledger, so they
   were left on their committed basis rather than silently re-pointed. The
   frozen-set equivalents exist (`caller_*_2sec_c6.csv`) and
   `figure_phase8_callers` carries the eleven-section version. **The PI should
   decide whether these three move to the frozen basis.**
2. **Another actor wrote `figures/` during this run.** `figure_gs2_crossarm_symmetry`
   and `figure_gs3_corescence_circularity` were byte-identical when I verified
   them at 06:06 and had different content by 06:55–06:57. **I first attributed
   this to the concurrent C7/D2 DeepScence job; that was wrong.** It was the
   fact-check-corrections agent, deliberately regenerating both to apply the
   `CDKN2B` map-gap fix and the re-derived CoreScence circularity — a legitimate
   correction, not drift. I inferred the culprit from the one other job I could
   see running and did not check; the timestamps told me *that* something wrote
   the files, never *who*. The detection was sound and the attribution was not,
   and the two should not have travelled in one sentence.
   The underlying gap was real: those two files are **uncommitted**, so
   `check_figures_guard.py` as it stood did not cover them and nothing warned.
   **Now fixed** — the guard walks all of `figures/` rather than `git ls-files`
   and covers **46 artefacts instead of 27**, with the manifest re-snapshotted to
   the post-8.7 state as the protected baseline.
3. **`make_figure_phase8_callers.py` crashed on the first attempt** and was
   fixed, not worked around: it selected gate rows by the two old basis labels,
   which the four-basis fix renamed. Panel (b) now draws the coverage effect at a
   **fixed** sender definition and keeps the published pre-C6 two-section value
   as a separate marker, so the figure can no longer imply that a pre-C6 point
   and a post-C6 point are the same comparison.

---

## 6. Section 17 — the two-arm table's mouse column, confirmed or corrected

Every row checked against a file this re-run produced. `code/m1_headlines.py`
generates the numeric rows from one results tree and reproduces **every** pre-C6
published value exactly when pointed at `results/phase3_pre_c6/`, which is what
licenses the "was" column.

| Quantity | §17 as written | **re-run value** | verdict |
|---|---|---|---|
| Platform / panel | Prime 5K Mouse + 100 custom | 5,106 `Gene Expression` features; **5,097** after removing 9 genotyping probes | **HOLDS**, add the 5,097 |
| Sections / donors | 11 / 11, 6 admissible | 11 / 11, 6 admissible | HOLDS |
| Cells | 1,834,806 total; 1,036,459 admissible | **1,826,893** total; **1,031,880** in-band, in the Phase 3 analysis set | **CORRECT IT** — the §17 figures are 7,913 / 4,579 higher than the set the fits actually use (cells that pass Xenium QC *and* carry a cell-type label). State which set is meant |
| Median NN distance (µm) | 6.7–9.7 | **6.74–10.61** over all 11 sections, and the same over the 6 in-band | **CORRECT IT.** 6.7–9.7 is the range over the four SBR sections of the original SBR-only scoping, not over the sections the paper reports |
| Transcript assignment rate | 88.27 % | 88.27 % (7259; unchanged, ingest not re-run) | HOLDS |
| **Sender prevalence, `tierA_p95`** | *(blank)* | **4.04–4.48 %** per in-band section, mean **4.29 %** | **FILL IN** |
| DeepScence coverage | 2/11 → 11/11 | 2/11 (364,291 cells) → **11/11 (1,826,894 cells)** | HOLDS, add the cell counts |
| DeepScence panel | ortholog-remapped, 4,845/5,097 | 4,845 / 5,097 = **95.06 %** | HOLDS |
| DeepScence sign vs own gene set | −0.350 sham / +0.318 SBR | unchanged — DeepScence was not re-run | HOLDS, but §4.2/§4.5's "reverses between arms" is a **7250-only** effect at 11 sections (`CS_PHASE8_CALLERS.md` §2.3) |
| **Caller agreement, depth- and type-matched** | 0.93–1.22× chance (2-section base) → restate at 11 | **0.751–2.198**, pooled **1.212×**, z 21.9, p 1.8e-106 (11 sec, frozen Tier A). Coverage effect at fixed Tier A: 1.040 → 1.129 pre-C6; 1.131 → 1.212 post-C6 | **RESTATE** — see `CORRECTIONS.md` §2 and §5.1 |
| Naive amplitude (response-sd) | 0.326 | **0.329** | HOLDS to 2 dp |
| Controlled amplitude (N2+N5+N6) | 0.027 | **0.029** | HOLDS |
| SF, N2+N5+N6 | 0.082 [−0.099, 0.249] | **0.088** [−0.017, 0.234] | HOLDS (the bracket is the IQR, not a CI — label it) |
| SF, N5 alone | 0.084 | **0.115** | **UPDATE** — same conclusion, 3 pp higher |
| SF, N2 matched decoy | 0.943 | **0.952** | UPDATE |
| λ̂ railed at a grid bound | 63 % | **60 %** | UPDATE |
| **Detectable bound (response-sd, 80 % power)** | 0.203 | **0.183** | **UPDATE** — the bound tightened |

The remaining rows, now that every stage has landed:

| Quantity | §17 as written | **re-run value** | verdict |
|---|---|---|---|
| **SF, N3 corrected — PRIMARY N3-var** | *C1 pending* | **0.996** [IQR 0.975, 1.007]; full design 0.997 | **FILL IN** |
| SF, N3 supporting variants (tile / occ15 / swap) | — | 0.971 / 0.940 / 0.695†; N3-occ = 0.302 at the literal 5 % tolerance | supporting |
| **SF, N4 corrected — PRIMARY N4-var** | *C1 pending* | **0.985** [IQR 0.958, 1.003]; full design 0.999 | **FILL IN** |
| SF, N4 supporting variants (tile / occ15 / swap) | — | 0.924 / 0.883 / 0.946; N4-occ = 0.183 | supporting |
| Null destructiveness, N3 published → corrected | 0.772 → 1.000 | **0.772 → 1.000** | **HOLDS EXACTLY** |
| Poisson identity r² | 0.984 | **0.984** (slope −0.525 against the geometric −0.500) | HOLDS |
| λ̂ interior median | — | 15.6 → **17.1 µm** | note if quoted |
| **Negative-control-probe kernel (must be flat)** | *new* | **NOT FLAT**: −0.074 SD, p = 0.015 naive; N5 removes it (+0.004, p = 0.72), **N2 does not** (−0.064, p = 0.012). Measured FPR **9–16 %** against a 5 % nominal | **FILL IN — and it is a negative result about the platform** |
| Composition surrogate share | 66–76 % | **now sourced** — **65.9 %** from receiver cell-type intercepts alone (reproducing the published 66 % almost exactly), **85.4 %** with the 20-NN composition vector added | **SPLIT THE ROW** |

**The primary corrected null is N3-var / N4-var, not N3-tile.** The
variance-corrected random shift (Mrkvička et al. 2021: shift on the plane, drop
what leaves the tissue window, standardise by retained sample size) replaces
N3-tile as the presented variant. The answer does not move — 0.996 against
N3-tile's 0.971 and a published 0.999 — but a direct type-I-error simulation
shows **tiling is more liberal than the whole-window torus it replaced, in 7 of
8 window × correlation-scale cells, up to 2.4× nominal (0.118 against 0.05)**,
while the variance-corrected estimator holds 0.040–0.060 everywhere. C1 replaced
a liberal test with a more liberal one, and no surviving fraction could have
revealed that. Full reasoning and both caveats in `CORRECTIONS.md` §8.1 and
`reports/CS_PHASE8_TORUS_VAR.md`.

† N3-swap is a label permutation, not a torus shift: it reproduces N1 (0.695
against N1's 0.707) and conditioning on the N5+N6+zonation block moves it to
**1.003**. The footnote drafted in `CS_PHASE8_C1_CLOSEOUT.md` §4.1 stands as
written with these three numbers substituted.

**The §17 "composition surrogate share" row — resolved, by someone else.**
I raised it as untraceable: no script emitted a quantity in that range under that
name, and the closest candidate (`attribution.csv`'s `sf_comp`, median **0.603**
both pre- and post-C6) did not match. `reports/CS_PHASE8_COMPMATCH.md` has since
sourced it, and my "untraceable" verdict was too strong — the producer existed,
it was the *pairing of two different measurements into one range* that had no
producer:

* **66 %** = 1 − 0.344, the surviving fraction when **receiver cell-type
  intercepts alone** are added (`CS_PHASE3.md` §3.1). Recomputed there as
  **65.9 %** (SF 0.3414 [0.236, 0.402]) — the published figure reproduces.
* **76 %** = 0.212 / 0.260, the **composition-only surrogate curve's** share of
  the unstratified contact amplitude (`CS_PHASE5.md` §4) — a different
  measurement on a different scale.
* Adding the 20-NN composition vector on top of cell-type intercepts reaches
  **85.4 %** (SF 0.1461 [0.052, 0.246]), **above** the published range entirely.

So "66–76 %" is not a range of one quantity; it is two numbers from two analyses
printed with an en-dash between them. **The pre-registration should split the row
rather than patch the range** — 66 % (own cell type) and 85 % (own cell type +
neighbourhood composition) are the defensible pair, with the 76 % figure kept
separately and labelled as the Phase 5 surrogate-curve share.

---

## 7. The second pre-registered Tier A variant — how it was implemented

Decision D1 makes `A_SENDER_FINAL_strict` (33 genes) **primary** and the seven
`A_sender_for_<module>.txt` sets the **pre-registered sensitivity analysis**.
That second variant is not another threshold on the same score: its sender set
**depends on which response module is being fitted**, because disjointness is
required only against that module's own readout.

The pipeline had no way to express that. It now does:

* `phase2_downstream.py` scores all seven per-module sets and writes
  `tierA_<module>_score` and `sender_flag_<module>_p{90,95,99}`;
* `Sec.sender_mask(call, module)` resolves `tierApm_pNN`;
* every Phase 3 stage that takes a sender call now takes an optional module, and
  the per-module calls **fan out over (section × module)** instead of
  (section × call) — 7× more jobs, each doing 1/7 of the fits — because a single
  `SectionFit` can no longer serve all seven modules.

**Reported against interest — the recorded cost estimate was wrong.**
`CS_PHASE8_C1_CLOSEOUT.md` §5.2 costed "both Tier A variants" at about 15 extra
minutes, on the assumption that they add three ordinary sender calls (main
66 → 99 jobs, N7 30 → 48). That accounting misses the module fan-out. In
practice the `main` stage went to **297** jobs, not 99, and the per-module
perturbation nulls are 42 jobs for one percentile rather than 6. The main-fit
half is still cheap; the perturbation-null half is not, and it is the reason
this pass runs the sensitivity variant's nulls at the **primary percentile p95
only** rather than at all three. `tierApm_p90` and `tierApm_p99` are carried
through the conditioning nulls (`main_fits.csv`) but not the perturbation nulls,
and that limit is stated here rather than left to be inferred from a missing
column.

**Note on the bio collaborator's recommendation.** `genesets/README.md` argues
the opposite of decision D1 — it calls the strict set "a numerically passing but
biologically hollow sender score" and recommends the per-module sets as primary.
D1 went the other way. The re-run shows the choice does not change the answer
(`CORRECTIONS.md` §3.3), which is the useful thing to be able to say either way,
but the disagreement should be visible in the pre-registration rather than
resolved silently.

---

## 8. Verification performed

* **Environment reproducibility.** Six of seven Tier B module scores and every
  anatomical/technical covariate reproduce **bit-identically** through a
  container rebuild (§3). This is the strongest statement available that the
  differences reported here are the gene sets and nothing else.
* **Caller code reproducibility.** `caller_disagree_all.py --verify`, pointed at
  the **pre-C6** `senders_*.csv` held aside before the re-run, reproduces all six
  committed two-section tables cell-by-cell: `VERIFY: PASS`. The code is
  unchanged; only `tierA_score` is.
* **Gene-set gate.** Re-verified independently against the panel derived from
  the h5, not from a recorded number. PASS on every criterion (§3).
* **Gene-set figures.** `figure_gs1`–`gs4` reproduce **byte-identically** after
  `crossarm_geneset_table.py` was fixed, confirming that no published gene-set
  number moved.
* **Destructiveness is geometric, as §17 claims.** `phase3_null_diag.py` re-run
  on the new sender sets: N3 bounding-box retention **0.772 → 0.772**, N4
  **0.917 → 0.920**, every in-tissue variant still 0.954–1.000, median real
  neighbour count 140.5 → 140.0. The §17 row "0.772 → 1.000" **holds exactly**.
* **Seed continuity.** `_expand()` reproduces the pre-Phase-8 seed formula for
  the six existing sender calls at the same section and call ordering, so the
  pre/post comparison is not confounded by a different random stream.
* **Figure guard.** `python3 code/check_figures_guard.py` run **before** the
  figure stage: `OK: all 27 committed figures match`. Run again after — §5.
* **Pinned files.** Two of three deliberately superseded, with both recovery
  paths verified; the two DeepScence files byte-identical (`CORRECTIONS.md` §6).
* **`data/raw_h1/` was not read or written.** `make_figure2.py` was not run.
  No commits, no pushes, no tags.

---

## 9. What this leaves for the freeze (8.9)

1. **The motivating claim must be restated once more.** `CS_PHASE8_CALLERS.md`
   §11's drafted restatement was written against the pre-C6 numbers. Under the
   frozen sets the effect is larger (1.212× rather than 1.118×) and, crucially,
   **no caller pair is reliably below chance any more**. The sentence about Tier A
   and SenePy being "genuinely close to disjoint" has to go with the sentence
   about statistical independence.
2. **`CS_PHASE8_CALLERS.md` §2.1 needs one number corrected**: 22 of 33 → **20
   of 33** (§5.2 of `CORRECTIONS.md`).
3. **`figures/revised_candidates/README.md`'s exemption for figures 2a and 2d is
   not valid for a gene-set change**, only for a null change. Both were
   regenerated here and the ledger has been updated with what actually happened.
4. **Decision D5's cost is now measured, not estimated.** Promoting C6 bought a
   biologically defensible B7 and a non-hollow Tier A; it cost a weaker apparent
   gradient in B7, seven reportable fits, and a further rise in caller agreement
   that makes the independence framing harder. Every one of those is in
   `CORRECTIONS.md`.
5. **`code/build_genesets.py` still cannot be re-run**, and the promotion has now
   made a second script depend on the tag rather than on `genesets/` for its
   pre-C6 column (§4.1). The freeze should record that `pre-c6-genesets` is
   load-bearing for provenance, not merely archival.
6. **The §17 "composition surrogate share 66–76 %" row has no producer.** It must
   be sourced or dropped before `PREREG_PHASE8.md`.

---

## 10. Reproduce

```bash
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# 1. gene sets -> sender calls, module scores, anatomy  (11 sections, ~5 min)
python3 -u phase2_downstream.py sham sbr \
    7001_liver_sham_Male_52-U1 7239_liver_sbr_Male_52-U1 7248_liver_sham_Male_26-U1 \
    7260_liver_sbr_Male_26-U1 7352_liver_sham_Male_2-U1 7361_liver_sbr_Male_2-U1 \
    7435_liver_sham_Male_10-U1 7448_liver_sbr_Male_10-U1 7450_liver_sbr_Male_10-U1

# 2. Phase 3 cache, forced                                        (~2 min)
python3 -u _run_prep.py

# 3. window, main fits (both Tier A variants), destructiveness
python3 -u run_phase3_nulls.py --stage window --sections all
python3 -u run_phase3_nulls.py --stage main --sections all --calls all9 --n-jobs 24
python3 -u phase3_null_diag.py

# 4. nulls, 1,000 permutations everywhere
python3 -u run_phase3_nulls.py --stage perm    --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
bash _m1_rerun_stage2.sh          # N7 at 1000, then tierApm_p95 at 1000, then curves

# 5. the eight small Phase 3 scripts (N8 fanned out over sections)
bash _m1_rerun_stage3b.sh

# 6. Phase 5
bash _m1_rerun_stage4.sh

# 7. callers, A7, summaries, and the ONE figure pass
python3 -u caller_disagree_all.py --all
bash _m1_rerun_stage5.sh
python3 code/check_figures_guard.py

# NOT make_figure2.py -- superseded producer of figure 2a
```

Pre-C6 baselines for every comparison in `CORRECTIONS.md`:
`results/phase3_pre_c6/`, `results/phase5_pre_c6/`, and `git tag pre-c6-genesets`.
---

## 11. Run supervision — an engineering note worth recording

Two mistakes were made supervising this run, and both are the kind that quietly
corrupt a pipeline rather than crash it.

**`pgrep -f` matches the shell that contains the pattern.** The first stage
driver waited on `pgrep -f "run_phase3_nulls.py --stage main"`. Every monitoring
shell whose own command line contained that string matched it, so the driver
waited on itself and the downstream stages never started. Detected because
`main_fits.csv` was written while the driver's log stayed empty. Fixed by
waiting on **log markers and file mtimes**, never on process-name patterns.

**A stage that looks serial can be the critical path.** `run_phase3_n8.py` has
no internal parallelism — 200 random gene sets per section, nine sections in a
`for` loop — and with `OMP_NUM_THREADS=1` set for the permutation stages it was
running at roughly one section per fifteen minutes, which would have made N8
alone longer than the entire rest of the re-run. It was stopped and fanned out
over sections (`_m1_rerun_stage3b.sh`), two BLAS threads each. Nothing was lost:
N8 writes one file per section and re-running a section is idempotent.

Neither affected a number. Both are recorded because the project has now been
bitten three times by *silent* pipeline faults — `make_figure2.py` writing the
same path as the live producer, `figure2a` caching its own input, and
`crossarm_geneset_table.py` reading the promoted sets as its own "pre" column.
The common shape is: **the output looks right because the input was stale.**
`code/check_figures_guard.py` now covers the figure half of that class; the data
half is still covered only by discipline.
---

## 12. Permutation counts — the §24.3 audit trail

Every perturbation-null file this re-run produced, with the count actually
recorded in the file rather than the count intended:

| file | `n_perm` | sender calls |
|---|---|---|
| `perm_nulls.csv` | **1000** | `tierA_p95` |
| `perm_nulls_c1.csv` | **1000** | `tierA_p95` |
| `perm_nulls_n7.csv` | **1000** | `tierA_p90`, `tierA_p99`, `cdkn1a_pos`, `senepy_p95`, `senepy_p99` |
| `perm_nulls_c1_n7.csv` | **1000** | the same five |
| `perm_nulls_pm.csv` | **1000** | `tierApm_p95` (per-module sender sets) |
| `perm_nulls_c1_pm.csv` | **1000** | `tierApm_p95` |

The pre-C6 N7 files were at 200. **That makes the N7 pre/post comparison
confounded** — five of the six calls changed permutation count *and* gene sets
at once — and `CORRECTIONS.md` §10 says so rather than attributing the movement
to the gene sets. The primary call `tierA_p95` was at 1,000 in both trees and is
the only clean row on that axis.

---

## 13. Closing state of the tree

`results/phase3/m1_final_audit.txt` regenerates everything below from the files
on disk.

**The three pinned files.** Two superseded deliberately, one preserved; the two
DeepScence files preserved.

```
d906394958dbe1b99981756290c511fa  results/phase3/perm_nulls.csv       SUPERSEDED (was 3b77aa1bba0712c205c5d9356654fb71)
a5ccc9b0e81f4c335e8039e975ec1975  results/phase3/sf_summary.csv       SUPERSEDED (was 69e3a1d3f60060deddcceba9896a7d31)
dc92ddc6605eef52f6359aeab4e16fd7  results/phase3/summary_phase3.txt   SUPERSEDED (was ecf86b9ca5460f31290e2f4c9e822ea2)
8c4c52f5c1c7649d8c17d07010cc780c  data/processed/deepscence_sham.csv  UNCHANGED
b557e3dfb8eff517d040757c73f0a660  data/processed/deepscence_sbr.csv   UNCHANGED
```

The pre-C6 content of all three is recoverable byte-identically from
`git show pre-c6-genesets:` and from `results/phase3_pre_c6/`.

**Wall clock.** 05:45 → 09:08 UTC, about **3 h 25 m** against the ~2 h 30 m
estimate. The overrun is entirely the decision to run every perturbation null at
1,000 permutations rather than 200 for the N7 axis and 1,000 rather than nothing
for the second Tier A variant: those two stages alone were 05:55 → 09:06.

**Not touched.** `data/raw_h1/` was neither read nor written. `make_figure2.py`
was not run. No commits, no pushes, no tags. No packages installed.

---

## 14. Follow-up pass, 2026-08-27 09:10–09:25 UTC

Three tasks closing out 8.7, plus one correction to this report.

### 14.1 The three `fig_phase3_*` figures — and my open item was overstated

I raised all three as sitting on the stale pre-C6 caller tables. **Only one of
them reads a caller table at all**, which I established by reading and re-running
their producers rather than by inspecting the figure scripts' filenames:

| figure | reads | gene-set dependent? | outcome |
|---|---|---|---|
| `fig_phase3_composition` | `results/composition_by_arm_timepoint.csv` — cell types only | **no** | **byte-identical**, and this is a genuine reproducibility pass |
| `fig_phase3_caller_depth` | `caller_within_type_depth_bias.csv` (committed, pre-C6) | **yes** | **re-pointed to `caller_within_type_depth_bias_2sec_c6.csv`; CHANGED** |
| `fig_phase3_tierC_identifiability` | `tierC_ligand_identifiability.csv` | **no, verified** | **byte-identical**, genuine |

The Tier C case is the one worth stating. `tierC_lr.py` *does* read
`senders_*.csv` and *does* pull `sender_flag_p95` at line 48 — which changed with
Tier A 25 → 33, so it looked gene-set dependent. **But `p95` is never used
anywhere else in the file**: it is dead. The senders it actually uses are
`cdkn1a_pos`, which is gene-set independent. Re-running `tierC_lr.py` on the
frozen senders confirms it: `tierC_ligand_identifiability.csv` and
`tierC_expression_by_celltype.csv` come back **`DataFrame.equals` identical**,
sender-dependent columns included. So Fig C's byte-identical regeneration is a
real pass, not a stale-input artefact. **My open item was right about one figure
of three and wrong about the other two**, and the difference was only visible by
re-running the producers.

Only `tierA_score` rows differ between the committed and frozen depth-bias tables
(max \|Δ\| **0.069**); `deepscence_score`, `senepy_score` and `cdkn1a_counts` are
identical. Tier A's Q5/Q1 bottom-selection is unchanged in direction.

All three now also emit a `*_data.csv` and a `.pdf`, which none of them did.

**Guard after this pass: exit 1, `CHANGED: figures/fig_phase3_caller_depth.png`
and nothing else.** New hash `8a68d69907a755e8773c089df065776a` (was
`51f9a93c4ae3004cd2b60d031f67764d`). `fig_phase3_composition.png`
`cb22425e055a2149c0d4e73d67636dba` and `fig_phase3_tierC_identifiability.png`
`6e2f06632dc9f7690b8b8f0a63a890ae` are unchanged. `--snapshot` was **not** run.

### 14.2 Deferred fact-check items applied

| item | where | applied |
|---|---|---|
| **R3** | `CS_PHASE7_C1.md` §1 + §6.2, `CS_PHASE8_C1_CLOSEOUT.md` §4.1 ‡ and §4.2 caption, `PHASE8_ROADMAP_STATUS.md` | **yes** — verified: out-of-tissue is `1 − frac_in_occupancy` = **35.5 % (N3) / 19.9 % (N4)**; the 23 %/8 % figures are `1 − frac_retaining_a_neighbour` = **22.8 % / 8.0 %**, "lost every real neighbour within 100 µm". Each site now states one quantity with its column name |
| **R4** | `CS_PHASE8_CALLERS.md` §4.4 | **yes** — confirmed by `md5sum`: `df2f0afd…`/`983b47b2…` are the `revised_candidates/*_REVISED.png` files; `git show HEAD:figures/figure2{b,c}.png` gives `5ecd9ad1…`/`a232a529…`. Only **2d** reproduced the repository. The "anomaly at 05:22:41" was the documented restore |
| **R5** | `CS_PHASE8_CALLERS.md` §4.1 stmt 2 | **yes** — "every" → "every count-based (4 of 5)". Also recorded that the exception has **disappeared** in the frozen re-run (`neg_probe_rate` n6n5 p 0.020 → 0.199) |
| **R6** | `CS_PHASE8_CALLERS.md` §4.3 | **yes, sentence withdrawn** — see below |
| **R9** | both C1 reports | **yes** — correction box in each |
| **M4** | `CS_PHASE8_CALLERS.md` §4.1 | **yes** — bin span (~0.27 SD) and amplitude ratio (0.070/0.291) separated |
| **M6** | `figure2c_data.csv` / `figure2c.png` | **already resolved by 8.7** — regenerated together at 09:07, 19 rows against 11 at `HEAD` |
| **M9** | A7, C1 keyed to pre-C6 senders | **already resolved by 8.7** — A7 fits 06:52, `a7_summary.csv` and `sf_summary_c1.csv` 09:06 |
| README L285–292 | caller-agreement paragraph | **yes** — restated at 1.212×, with the superseded 1.51–2.85× and 2.15×/0.38× figures corrected |

**One site left for its owner.** `reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` L285
still carries "0.93–1.22× of chance for four of six pairs, but 1.51–2.85× for
DeepScence vs `Cdkn1a`⁺" as a **live** recommendation — it is that report's
*suggested rewording* for the README, so applying the README fix without fixing
it leaves the two documents contradicting each other. It is bio-owned and not in
the audit's §4 list, so I did not edit it. **It needs the same restatement.**
Every other surviving instance of those strings in the repo is a quotation
inside a correction box, which is intended.

**R6 — the audit was right in substance and its numbers were stale.** It gave
3.0–6.7 % (full design) and 1.2–13.3 % (naive). Those are correct for the A7 file
it had, the **pre-C6 run of 05:19**, which this re-run superseded at 06:52 when
the sender call changed. Measured on the 825 control fits of the frozen run, the
reportable filter admits **4.8 % on the full design — identical across all five
control families, i.e. essentially nominal — and 3.0–13.3 % on the naive
design**. The conclusion is unchanged and stronger: the "2–3× nominal" bound
applies to the *estimator*, not to the *filter*, and the filter sentence is
withdrawn.

### 14.3 N3-var / N4-var adopted as the primary corrected null

See `CORRECTIONS.md` §8.1 and `reports/CS_PHASE8_TORUS_VAR.md`. **N3-var 0.996,
N4-var 0.985**, against N3-tile 0.971 / N4-tile 0.924 and a published
0.999 / 0.947. The answer does not move; the calibration does. A type-I-error
simulation on Mrkvička's own §5 design, extended with an irregular window,
measures **tiling at 0.040–0.118 against a 5 % nominal (up to 2.4×)**, the
whole-window torus at 0.033–0.073, and the variance-corrected estimator at
**0.040–0.060**. Correction C1 replaced a liberal test with a more liberal one —
a fact no surviving fraction could expose, because SF is not a rejection rate.
N3-tile, occ, occ15, swap and snap stay in the battery and in Figure 2c as
supporting variants.

### 14.4 A correction to §5 of this report

I attributed the 06:55–06:57 rewrite of `figure_gs2` and `figure_gs3` to the
concurrent DeepScence DCA job. **That was wrong** — it was the
fact-check-corrections agent, deliberately applying the `CDKN2B` map-gap fix and
the re-derived CoreScence circularity. I inferred the culprit from the one other
job I could see running and did not check; the timestamps told me *that*
something wrote the files, never *who*. The detection was sound, the attribution
was not, and the two should not have travelled in one sentence. The gap the
detection exposed was real and is now closed: the guard walks all of `figures/`
and covers 46 artefacts against 27.
