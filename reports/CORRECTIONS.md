# CORRECTIONS — what moved today, why, and what it costs the paper

**Task 8.8. The ledger the addendum §16 step 3 asks for: for every headline
number, the pre value, the post value, and the attributed cause.** Where a number
moved for more than one reason, both are named rather than the tidier one being
picked.

**Scope.** §0–§12 are the ledger of the 8.7 M1 end-to-end re-run (companion to
`reports/CS_PHASE8_M1_RERUN.md`, which is the run record). **§13–§18 were added
at 10:50 UTC and cover the six findings that landed after that re-run** — the
composition-matched result, Figure 1's planted ground truth, the corrected caller
diagnosis, the D2/DCA resolution, Moran's I, and the corrections to the record
itself. Every number in this file was read from a named file, and the ones that
carry the argument were re-derived independently for this pass, because two
audits have already caught invented figures and a bibliography written from
recall.

---

## A. THE HEADLINE DID NOT MOVE — read this first

**Every correction made today moved against interest. None of them changed the
conclusion.** That is the frame for everything else in this file, and it is the
strongest position a negative result can be in.

Re-derived for this ledger by running `code/m1_headlines.py` directly against
`results/phase3/`:

| quantity | pre-C6 (published) | **post-C6 (frozen)** | verdict |
|---|---|---|---|
| **controlled amplitude, N2+N5+N6** | 0.0272 | **0.0288** response-sd | unchanged in substance and sign |
| **detectable bound, 80 % power** | 0.2032 | **0.1833** response-sd | **tightened** |
| ratio of the two | 7.5× | **6.4×** | the bound is still far above the estimate |
| **SF under N2+N5+N6** | 0.0822 [−0.099, 0.249] | **0.0885** [−0.017, 0.234] | unchanged |
| naive amplitude, median \|β̂\|/sd(y) | 0.3262 | **0.3288** | unchanged |
| controlled fits positive with CI excluding 0 | 15 / 160 | **13 / 153** | unchanged |

**Controlled amplitude 0.029 against a detectable bound of 0.183 at 80 % power;
SF under N2+N5+N6 = 0.088. The addendum's §18 outcome A stands.** *No
distance-dependent SASP kernel is identifiable at achievable power* — unchanged
in substance, in sign, and now on a **tighter** bound than the published one.

**It survived, specifically:** a full gene-set re-sourcing (§0, §3), a container
rebuild (§1), a corrected in-tissue null battery (§8), a *second* correction that
replaced the corrected null again (§8.1), a doubling of caller coverage (§2), a
second pre-registered sender definition (§3.3), and every item in §13–§17 below.
The C1 verdict is likewise unchanged.

*(Note on cross-references: **§ numbers without a document name are sections of
this file.** References to the addendum `Phase7_Minimal_Human_Replication (1).md`
and to `SASP_Kernel_Master_Plan.md` are named as such — the two numbering schemes
collide at §17 and §18, which is why.)*

**And it is now defended better than it was.** The geometric predictions the
addendum's §17 says *should* replicate both do: Poisson identity slope
**−0.5249** against a geometric −0.5, r² **0.9843**. Six of seven Tier B module scores and every
anatomical and technical covariate return **bit-identical** across a container
rebuild (§1). The attribution in this file is measured, not assumed.

---

## B. Every number that moved, in one table

Each row points at the section that carries the derivation and the file.

| # | quantity | pre | post | attributed cause | § |
|---|---|---|---|---|---|
| 1 | mouse B7 `secondary_senescence` | **38** genes | **108** | C6 re-sourcing from an arrest-flavoured set to a paracrine secretome set | §0.1 |
| 2 | strict Tier A (PRIMARY sender set) | **25** | **33** | **consequence of row 1**, not an edit to Tier A: the gate builds Tier A by *subtracting* Tier B, so shrinking B7's arrest content released 8 genes | §0.1, §14.2 |
| 3 | `A_sender_for_secondary_senescence` | 55 | 74 | same mechanism | §0.1 |
| 4 | reportable fits (in-band, `tierA_p95`) | **160** | **153** | B7 alone loses 14; the other six modules gain 7 | §3, §3.1 |
| 5 | B7 reportable fits | **36** | **22** | B7's response score changed — the only module whose did | §3.1 |
| 6 | B7 naive amplitude | **0.342** | **0.246** | making B7 *actually* paracrine **reduced** its distance gradient | §3.1 |
| 7 | caller agreement, 3-pair pooled, 2-section **published base** | **1.030** (p = 0.20) | **1.128** (z 5.47, **p = 4.4e-8**) | the frozen Tier A alone — sections held fixed | §5.1, §14.1 |
| 8 | caller agreement, 3-pair pooled, 11 sections | 1.030 | **1.212** (p = 1.8e-94) | gene sets **and** coverage — the two are separated in §5.1, not merged | §2.1, §5.1 |
| 9 | Tier A vs SenePy | 0.914, below chance 11/11 | **0.972, n.s., above chance 4/11** | new Tier A. **No Tier A pair is reliably below chance any more** | §2.2 |
| 10 | Tier A vs `Cdkn1a`⁺ | 1.171 | **1.471**, 11/11 | 4 of the 8 re-entering genes sit on the p53 axis. **Biological correlation, not circularity** — `Cdkn1a` is still not in Tier A | §2.2, §14.2 |
| 11 | "22 of 33 sections above chance" | 22 | **20** | arithmetic error in `CS_PHASE8_CALLERS.md` §2.1, **still uncorrected there** | §5.2, §18.3 |
| 12 | CoreScence circularity, mouse | "**69 %**" | **79 % → 88 %** | the 69 % was a **typed-in literal** (`24/35`) whose denominator no mapping convention reproduces. Re-derived 26/33 pre-C6, 29/33 post | §17.7 |
| 13 | CoreScence circularity, human | "76 %" | **76 % → 88 %** | measured within-arm; both arms land on 29/33 | §17.7 |
| 14 | cost of C6 to circularity | "+19 points" | **+9 mouse, +12 human** | the old story compared a *mouse* baseline against two *human* configurations | §17.7 |
| 15 | Tier A ∩ CoreScence | 2/25 | **4/33** (`Cdkn2b`, `Mdm2`) | same 8 re-entrants as row 2 | §17.7 |
| 16 | SF, N5 alone | 0.084 | **0.115** | new Tier A sender set | §3 |
| 17 | SF, N2 matched decoy | 0.943 | **0.952** | new sender set → new matched decoys | §3 |
| 18 | λ̂ railed at a grid bound | 63 % | **60 %** | new sender set | §3 |
| 19 | block-bootstrap SE of the controlled amplitude | 0.0725 | **0.0654** | fewer, tighter reportable fits | §3 |
| 20 | primary corrected null | N3-tile / N4-tile | **N3-var / N4-var** | a **calibration** correction, not a numerical one: tiling is *more* liberal than the whole-window torus it replaced | §8.1 |
| 21 | median NN distance | 6.7–9.7 µm | **6.74–10.61 µm** | the published range was over 4 SBR sections, not the 11 the paper reports. **Never right for the stated scope** | §7 |
| 22 | cells | 1,834,806 / 1,036,459 | **1,826,893 / 1,031,880** | the analysis set is QC-passing **and** cell-type-labelled | §7 |
| 23 | composition surrogate share | "66–76 %" | **66 % → 85 %**, split into two rows | not one quantity: 66 % is 1−SF on β̂, 76 % is a curve-amplitude ratio. The **pairing** had no producer | §7, §18.4 |
| 24 | A7 per-response amplitudes | pre-C6 vintage in 8 documents | **post-C6 (09:06)** | A7 was re-run under the frozen sender calls; the digits in circulation are the 05:19 file | §17.6 |
| 25 | A7 biological-module amplitude | 0.291 / 0.036 | **0.277 / 0.031** | same cause as row 24 | §16.2, §17.6 |
| 26 | DeepScence depth loading, `denoise=True` | *unmeasured* | **+32 % to +67 %** over `denoise=False`, 3 of 3 sections | §4 (D-b) predicted the opposite sign | §15 |
| 27 | H1 tissue | "human aging lung" | **7 human spleens** (GSE326743) | wrong organ on **15 lines across 7 sections** of the plan of record | §17.3 |
| 28 | bibliography | 32 entries, 19 with invented forenames | **43 entries**, 18 author lines corrected | written from recall; D7 had checked DOIs but not forenames | §17.1 |
| 29 | torus-null attribution | CellWHISPER | **Lotwick & Silverman (1982)** | CellWHISPER's null is a within-cell-type location permutation — this project's **N1**. Corrected in 9 places | §17.2 |

---

## C. The four findings that changed what the paper can claim

Full derivations in the sections named. In brief:

1. **Matched decoys fail where covariates work — now shown three times
   independently** (§13). A7: N5 removes the platform's technical gradient, N2
   does not. Composition matching: **1.6 % against 85.4 %** from the *same
   variables* — a factor of fifty. And **Figure 1's own committed synthetic data,
   on planted ground truth**: matched-decoy λ CI coverage **0.346** against naive
   **0.513** and nuisance-conditioned **0.854** — *worse than doing nothing* — with
   |relative bias| worse than naive in **8 of 20** grid cells. **The third is the
   strongest, because the kernel is planted**, so "your covariates removed real
   signal" does not apply. It has been in the project's headline figure since
   Phase 1.
2. **Caller independence is dead, and the cause is not what I first reported**
   (§14). It fails on the **published two-section base** once the sender set is
   frozen: **1.128×, p = 4.4 × 10⁻⁸**. The defect was in the **response** module —
   pre-C6 B7 had absorbed 27 genes, 19 of them from the Tier A candidate pool, and
   the disjointness gate resolved that by **gutting Tier A to a 25-gene remnant
   by subtraction**. Not "a circular sender definition".
3. **§4 (D-b)'s premise is refuted** (§15). DCA installed and ran — §6 path 1 in
   substance, not the fallback — and **raises** the depth loading by 32–67 % on 3
   of 3 sections. `denoise=False` is now a **chosen** value, not a deviation. The
   published default is also seed-unstable: top-5 % Jaccard **exactly 0.000**
   between two seed pairs, **with no diagnostic firing**.
4. **The torus finding is not novel statistics** (§8.1, §17.2) — Lotwick &
   Silverman (1982). What is new is the **measurement**: tiled torus runs at up to
   **2.35× nominal** type-I error while the variance correction holds 0.033–0.060.
   **C1 replaced a liberal test with a more liberal one**, and no surviving
   fraction could have revealed it.

**Add to those, from §16:** Moran's I was run and **falsified an argument the
project was about to make** — that Moran's I and A7 "disagree". They rank
together at ρ = +0.895. The replacement is stronger because it is measured:
Moran's I could not have detected the paper's own headline effect either.

---

## D. What is still wrong or unverified

§18, in full, and it is not short. The two lines that matter most:

* **Nothing is committed.** The last commit is 2026-08-21; **482 untracked and
  107 modified files** hold the entire Phase 8 evidence base, including
  `genesets/human/`, `genesets/mouse_c6/`, `results/phase3_pre_c6/` (the sole copy
  of the baseline this whole ledger compares against) and 23 of 34 files in
  `reports/` — **this one included**. `phase8-frozen` would reference nothing
  immutable. Worse: **a `git checkout` is currently destructive**, because
  `code/build_genesets.py` at `HEAD` would overwrite `genesets/*.txt` with empty
  Tier B modules and the guard exists only in the working tree.
* **`SUBMISSION_PATCH_2026-08-29.md` §9 still instructs the falsified Moran's I
  framing**, in the one document the PI applies to the manuscript by hand, two
  days before the deadline.

---

## 0. What changed upstream, and the causal chain between the changes

Four things were carried together for the first time. They are **not
independent**, and the ledger below depends on knowing which caused which.

| # | Change | Where |
|---|---|---|
| 1 | Corrected N3/N4 in-tissue nulls (C1) | `run_phase3_nulls.py --stage perm_c1` |
| 2 | Promoted C6 mouse gene sets (decision D5) | `genesets/` — 3 of 15 files changed |
| 3 | Both pre-registered Tier A variants on the sender axis | new `tierApm_p*` calls |
| 4 | DeepScence at 11/11 section coverage (C7/D1) | `data/processed/deepscence_*.csv` |

### 0.1 The gene-set change is one change, not three

`git diff pre-c6-genesets -- genesets/` touches exactly three files, and the
second and third are **consequences of the first**:

| file | pre-C6 | post-C6 | added | removed |
|---|---|---|---|---|
| `B_secondary_senescence.txt` (B7) | 38 | **108** | 97 | 27 |
| `A_SENDER_FINAL_strict.txt` (Tier A, PRIMARY) | 25 | **33** | 8 | 0 |
| `A_sender_for_secondary_senescence.txt` | 55 | **74** | 19 | 0 |

B7 was re-sourced from an arrest-flavoured set to a genuine paracrine
secondary-senescence/SASP set. It **lost** 27 canonical arrest genes —
`Atm Bax Bcl2l1 Ccnd1 Cdkn1a Cdkn1b Cdkn2a Cdkn2b Chek2 Egr1 Ezh2 Fos Gadd45a
Glb1 Hmgb2 Junb Lmnb1 Mdm2 Nfkbia Sirt1 Sod2 Suv39h1 Tgfb1 Thbs1 Timp1 Trp53
Trp53bp1` — and **gained** 97 secretome genes (`Ccl*`, `Cxcl*`, `Mmp*`, `Il*`,
`Vegf*`, `Tnf`, `Igf1`, `Hgf`, …).

Tier A is defined by removal of Tier B collisions. Eight of the 27 genes B7 gave
up (`Atm Bax Bcl2l1 Cdkn2b Glb1 Mdm2 Sirt1 Trp53bp1`) collide with **no** Tier B
module any more, so they re-enter Tier A. **The strict Tier A grew 25 → 33
because B7 shrank its overlap, not because anyone added genes to Tier A.** The
same mechanism gives `A_sender_for_secondary_senescence` 55 → 74.

Consequence for attribution: any number that moved has at most two proximate
causes — **a different B7 response score**, or **a different Tier A sender
score** — and both trace to the single B7 re-sourcing.

### 0.2 The other twelve gene sets are byte-identical

Verified against `git show pre-c6-genesets:` for all 15 mouse files;
12 unchanged, 3 as tabulated. `genesets/mouse_c6/*.txt` is byte-identical to the
promoted `genesets/*.txt` for all 15 files, so the promotion is complete and
faithful.

### 0.3 The disjointness gate, re-verified independently

Re-run here against the authoritative mouse panel derived from the data, not
from a recorded number: `cell_feature_matrix.h5` carries 5,106 `Gene Expression`
features, minus 9 `ENSMUSG`-less genotyping probes = **5,097 panel genes**
(decision D7).

| criterion | required | measured | verdict |
|---|---|---|---|
| `len(Tier A strict)` | ≥ 15 | **33**, all 33 on panel | PASS |
| `len(B[m])` per module | ≥ 30 | 190 / 125 / 68 / 100 / **31** / 108 / 126 | PASS (all 7 on panel) |
| `len(A ∩ ∪B)` | 0 | **0** | PASS |
| per-module `A_sender_for_m` | ≥ 15 and `∩ B[m] = 0` | 37–74 genes, all intersections 0 | PASS |

`B_oxidative_stress` still passes on a margin of exactly **1** gene (31 vs the
floor of 30) — the live risk the roadmap records. It is unaffected by C6.

**Reported against interest:** the Tier B modules are **not** mutually disjoint —
18 of the 21 unordered module pairs share at least one gene. Section 8 Test 2
never required that, and no gate fails, but any statement of the form "the seven
response modules are disjoint" is false and must not be made.

---

## 1. The re-run reproduces everything the gene sets did not touch

This is the control that makes the rest of the ledger an attribution rather than
a guess. `phase2_downstream.py` was re-run on all 11 sections in the rebuilt
environment. Compared cell-by-cell against the pre-C6 outputs
(`data/processed/` as captured before the re-run):

| quantity | max |Δ| over 93,197 cells of 7450 | attributed to |
|---|---|---|
| `cdkn1a_counts`, `cdkn1a_pos` | **0** | — |
| `senepy_score` | **0** | — |
| `zonation_score`, `compartment_label`, `dist_to_boundary_um`, `dist_to_portal_triad_um` | **0** | — |
| module scores: `downstream_arrest`, `emt_ecm`, `il6_jak_stat3`, `interferon_response`, `oxidative_stress`, `tnfa_nfkb_proximal` | **0** | — |
| module score `secondary_senescence` | **0.43** | B7 re-sourced, 38 → 108 |
| `tierA_score` | **0.244** | strict Tier A 25 → 33 |

Six of seven module scores and every anatomical and technical covariate come
back **bit-identical** across a container rebuild. The environment is not a
confound, and exactly two inputs changed.

---

## 2. The caller-agreement gate moved again, and it moved further against us

**This is the most important line in the ledger.** Task 8.4's gate result was
computed on the pre-C6 Tier A. The promoted Tier A changes three of the six
caller pairs, and it pushes the headline further from independence.

Attribution is exact, because `caller_disagree_all.py` was first re-verified
against the **pre-C6** `senders_*.csv` (held aside before the re-run): all six
committed two-section tables reproduce cell-by-cell, `VERIFY: PASS`. The code
did not change; only `tierA_score` did.

### 2.1 Headline (`results/phase3/caller_coverage_gate_headline.csv`)

| basis | band | median | pooled ratio | z | p | above chance |
|---|---|---|---|---|---|---|
| 2-section (published) | 0.932 – 1.369 | 1.010 | 1.040 | 1.76 | 0.078 | 4/8 |
| 11-section, **pre-C6 Tier A** (task 8.4) | 0.700 – 1.711 | 1.156 | **1.129** | 13.35 | 1.1e-40 | 29/44 |
| 11-section, **post-C6 Tier A** (this re-run) | **0.751 – 2.198** | **1.190** | **1.212** | **21.92** | **1.8e-106** | **35/44** |
| 6-section in-band, post-C6 | 0.811 – 1.565 | 1.168 | 1.167 | 13.09 | 3.6e-39 | 18/24 |

### 2.2 Per pair — the three non-Tier-A pairs are bit-identical

| pair | pre-C6 pooled (11 sec) | post-C6 pooled | above chance | cause |
|---|---|---|---|---|
| Tier A vs `Cdkn1a`⁺ | 1.171 (z 7.1) | **1.471 (z 19.5)** | 9/11 → **11/11** | Tier A gained 8 p53-axis genes |
| Tier A vs DeepScence | 1.248 (z 16.6) | **1.288 (z 19.2)** | 11/11 → **11/11** | " |
| Tier A vs SenePy | 0.914 (z −4.9) | **0.972 (z −1.6, p = 0.10)** | 0/11 → **4/11** | " |
| SenePy vs `Cdkn1a`⁺ | 1.211 | 1.211 | unchanged | — |
| SenePy vs DeepScence | 0.737 | 0.737 | unchanged | — |
| *(circular)* DeepScence vs `Cdkn1a`⁺ | 1.255 | 1.255 | unchanged | — |

**Two statements from `CS_PHASE8_CALLERS.md` must now be withdrawn.**

1. **"Tier A vs SenePy is the one pair that does not move: 0.935 → 0.914, below
   chance on all eleven sections. Two callers are genuinely close to disjoint.
   That is worth keeping."** It does not survive. Post-C6 it is **0.972,
   z = −1.63, p = 0.10, above chance in 4 of 11** — statistically
   indistinguishable from independence, not below it. There is no longer any
   pair that is reliably *below* chance.
2. The per-pair number **"Tier A vs `Cdkn1a`⁺ ≈ 1.17×"** is now **1.47×**, the
   largest non-circular agreement in the matrix and larger than the circular
   DeepScence–`Cdkn1a`⁺ pair (1.255).

**The same three pairs on the two-section base, so the gene-set effect is
visible without the coverage change** (`caller_coverage_gate.csv`, rows tagged
`2-section`): Tier A vs `Cdkn1a`⁺ **1.017 → 1.300**, Tier A vs DeepScence
**1.103 → 1.179**, Tier A vs SenePy **0.935 → 1.007**. The three non-Tier-A
pairs are identical on that base too.

**Why Tier A vs `Cdkn1a`⁺ jumped, stated plainly.** The eight genes that
re-entered Tier A are `Atm Bax Bcl2l1 Cdkn2b Glb1 Mdm2 Sirt1 Trp53bp1` —
four of them (`Atm`, `Mdm2`, `Trp53bp1`, `Bax`) sit on the p53 axis that induces
`Cdkn1a`. The pre-C6 Tier A was, in the bio collaborator's own words in
`genesets/README.md`, "a numerically passing but biologically hollow sender
score" containing no `Cdkn1a`, no `Cdkn2a`, no `Trp53`, no `Lmnb1`, no `Mki67`.
Making it less hollow necessarily made it correlate more with a `Cdkn1a`-based
call. **`Cdkn1a` itself is still not in Tier A — disjointness holds — so this is
biological correlation, not circularity.** But it is a real cost of C6 and the
independence framing pays it.

### 2.3 What did NOT move in the caller tables

* `caller_technical_loading_11sections.csv`: SenePy, DeepScence and `Cdkn1a`⁺
  rows bit-identical. Tier A's own depth loading **improved**: Spearman ρ with
  transcript counts 0.008–0.163 → **−0.061–0.125**; with cell area
  −0.023–0.164 → −0.082–0.123. The promoted Tier A is *less* technically loaded.
* `caller_within_type_depth_bias_11sections.csv`: §4.1's qualitative claim
  survives; Tier A's Q5/Q1 depth enrichment goes 0.174–0.291 → **0.146–0.317**
  (still bottom-selecting in all eleven sections). The quoted range must be
  updated.
* The 7250-only anomalies (DeepScence depth-correlation sign flip, hepatocyte
  under-calling, SenePy-vs-DeepScence 2.150) are untouched: they are DeepScence
  properties and DeepScence did not change.

---

## 3. Phase 3 — the headline null battery

Primary sender call `tierA_p95`, six Test-3-admissible sections, receiver
labels `merged`, λ held at λ̂, 400 spatial block bootstrap replicates over 100
quantile blocks. Reportable population unchanged in definition: positive naive
amplitude whose block-bootstrap CI excludes zero.

Both columns are produced by the same extractor, `code/m1_headlines.py`, run
once on `results/phase3_pre_c6/` and once on `results/phase3/`. It reproduces
every published pre-C6 number exactly, which is what licenses the comparison.

| quantity | pre-C6 (published) | post-C6 (this re-run) | attributed cause |
|---|---|---|---|
| fits | 315 | 315 | — |
| **reportable fits** | **160** | **153** | B7's new score is a weaker naive gradient (§3.1) |
| naive amplitude, median \|β̂\|/sd(y) | **0.326** | **0.329** | net of two offsetting moves (§3.1) |
| **controlled amplitude (N2+N5+N6)** | **0.027** [−0.028, 0.090] | **0.029** [−0.007, 0.084] | as above |
| block-bootstrap SE of the controlled amplitude | 0.0725 | 0.0654 | fewer, tighter reportable fits |
| **detectable bound, 80 % power** | **0.203** | **0.183** | " |
| controlled fits positive with CI excluding 0 | 15 / 160 | 13 / 153 | " |
| SF, N2 matched decoy | 0.943 | **0.952** | new Tier A sender set → new decoys |
| SF, N5 alone | 0.084 | **0.115** | new Tier A sender set |
| SF, N6 alone | 0.486 | 0.471 | " |
| SF, zonation alone | 0.843 | 0.843 | — |
| SF, N5+N6 | 0.082 | 0.084 | " |
| **SF, N2+N5+N6** | **0.082** [−0.099, 0.249] | **0.088** [−0.017, 0.234] | " |
| SF, N1 stratified label permutation | 0.716 | **0.707** | " |
| SF, N1 on the N5+N6-conditioned residual | 0.987 | 0.989 | " |
| SF, N3 torus shift (`perm_nulls.csv`) | 1.000 | **0.999** | " |
| SF, N4 rotation (`perm_nulls.csv`) | 0.964 | **0.947** | " |
| SF, N8 scrambled response gene set | 0.925 | **0.916** | new Tier A sender set; the Tier E3 decoys are unchanged |
| λ̂ railed at a grid bound | 63 % | **60 %** | " |

**The contribution stands, and stands slightly more comfortably.** The
controlled amplitude 0.029 is still far below the design's 80 %-power detectable
bound of 0.183 response-sd, so §17's negative result — *no distance-dependent
SASP kernel is identifiable at achievable power* — is unchanged in substance and
in sign. The bound itself tightened from 0.203 to 0.183.

### 3.1 Where the movement actually is: one module

Per module, over the reportable fits (`results/phase3/m1_prepost_main_fits.txt`):

| module | reportable pre → post | naive \|β̂\|/sd pre → post | SF N5+N6 pre → post |
|---|---|---|---|
| downstream_arrest | 7 → 8 | 0.246 → 0.202 | 0.753 → 0.666 |
| emt_ecm | 33 → 35 | 0.538 → 0.452 | 0.018 → 0.018 |
| il6_jak_stat3 | 30 → 32 | 0.318 → 0.349 | −0.015 → 0.084 |
| interferon_response | 21 → 21 | 0.294 → 0.318 | 0.161 → 0.087 |
| oxidative_stress | 0 → 0 | — | — |
| **secondary_senescence (B7)** | **36 → 22** | **0.342 → 0.246** | 0.090 → 0.138 |
| tnfa_nfkb_proximal | 33 → 35 | 0.305 → 0.312 | 0.101 → 0.064 |

**B7 alone loses 14 reportable fits; the other six modules gain 7 between them,
for a net loss of 7.** B7 is the only module whose *response score* changed. The re-sourced B7 — a paracrine
secretome set rather than an arrest set — has a **weaker naive distance
gradient**, so fewer of its fits clear the "positive amplitude, CI excludes
zero" bar. Every other module moved only through the changed **sender** set, and
those moves are small and unsigned.

This matters for the paper's framing. B7 was the module whose name most directly
promised a *paracrine* effect, and making it actually paracrine **reduced** the
apparent distance dependence. That is the honest direction: the pre-C6 B7 shared
27 arrest genes with the sender-adjacent literature, and the gradient it showed
was partly a gradient in arrest, not in secondary senescence.

### 3.2 The reportable population, receiver type by receiver type

No receiver type gains or loses its qualitative behaviour. The two visible
shifts are **vSMCs 12 → 6 reportable fits** and **Endothelial naive amplitude
0.367 → 0.255**; both are sender-set effects (their module scores are
bit-identical) and both leave the controlled SF within its pre-C6 IQR.

### 3.3 The N7 sender-definition axis, and the second Tier A variant

Both pre-registered Tier A variants are now on the axis, so it carries **nine**
sender calls, not six. `tierA_p*` is scored on `A_SENDER_FINAL_strict` (33
genes, **PRIMARY**, roadmap decision D1); `tierApm_p*` is scored on the
per-module `A_sender_for_<module>.txt` sets (**pre-registered sensitivity**), so
its sender mask is module-specific and every stage carries the module down to
`Sec.sender_mask`.

| call | reportable pre → post | naive \|β̂\|/sd pre → post | SF N5 pre → post | SF N2+N5+N6 pre → post |
|---|---|---|---|---|
| tierA_p90 | 198 → 181 | 0.395 → 0.392 | 0.069 → 0.122 | 0.088 → 0.094 |
| **tierA_p95 (PRIMARY)** | 160 → 153 | 0.326 → 0.329 | 0.084 → 0.115 | **0.082 → 0.088** |
| tierA_p99 | 91 → 93 | 0.320 → 0.320 | 0.246 → 0.207 | 0.247 → 0.188 |
| cdkn1a_pos | 203 → 192 | 0.506 → 0.462 | 0.385 → 0.388 | 0.203 → 0.205 |
| senepy_p95 | 155 → 146 | 0.476 → 0.475 | 0.302 → 0.296 | 0.219 → 0.208 |
| senepy_p99 | 107 → 103 | 0.516 → 0.516 | 0.444 → 0.446 | 0.283 → 0.284 |
| **tierApm_p90** *(new)* | — → 185 | — → 0.381 | — → 0.108 | — → 0.086 |
| **tierApm_p95** *(new)* | — → 160 | — → 0.297 | — → 0.136 | — → **0.096** |
| **tierApm_p99** *(new)* | — → 92 | — → 0.320 | — → 0.243 | — → 0.193 |

**The sensitivity variant agrees with the primary at every percentile.** The
controlled surviving fraction is 0.086 / 0.096 / 0.193 against the primary's
0.094 / 0.088 / 0.188, and the naive amplitudes are within 0.03 sd. The choice
between the two pre-registered Tier A definitions does **not** change the
conclusion, which is what the pre-registration needs to be able to say.

For the record, the two definitions call **genuinely different cells**: on 7259
the per-module p95 masks overlap the strict p95 mask by only 46–72 %
(`downstream_arrest` 0.724, `emt_ecm`/`il6_jak_stat3`/`secondary_senescence`
0.459). The agreement in outcome is not an artefact of the sets being nearly the
same.

---

## 4. DeepScence / CoreScence — the §9 reporting standard, as one Methods row

Addendum §9 requires four attributes on every DeepScence number: **coverage**,
**denoise state**, **anchor**, **panel** (native or ortholog-mapped, with the
mapping rate). One row, not scattered caveats.

| attribute | published (2-section) | this re-run (11-section) |
|---|---|---|
| **coverage** | 2 of 11 M1 sections — 7250 (`deepscence_sham.csv`) and 7259 (`deepscence_sbr.csv`), **364,291 cells** | **11 of 11** M1 sections, **1,826,894 cells** |
| **denoise state** | `denoise=False` (DCA/TensorFlow unavailable in the pinned stack) | `denoise=False`, identical setting; C7/D2 quantification is a separate task (8.5) |
| **anchor** | published `CDKN1A` | published `CDKN1A` (**primary**); the D3 re-anchored variant is a declared sensitivity, `deepscence_d3_<section>.csv` |
| **panel** | ortholog-mapped, **4,845 of 5,097** mouse panel genes = **95.06 %** MGI 1:1 mapping rate | identical: 4,845 / 5,097, 95.06 % |
| other fixed settings | `random_state=0`, ≥20 counts/cell | identical |

**The two original files are byte-identical after the re-run**, as required:

```
8c4c52f5c1c7649d8c17d07010cc780c  data/processed/deepscence_sham.csv   (7250)
b557e3dfb8eff517d040757c73f0a660  data/processed/deepscence_sbr.csv    (7259)
```

The nine added sections are the coverage change and nothing else: settings are
identical to `run_deepscence.py`, and the 7250/7259 rows inside the
eleven-section caller tables reproduce the committed two-section rows exactly.

**Two-section values are reported alongside eleven-section values throughout**
(§2.1 above, and `caller_coverage_gate.csv` carries both bases in every row), as
§5 of the addendum requires.

---

## 5. Two errors in earlier Phase 8 reports, found while doing this

### 5.1 `caller_coverage_gate.csv` was pooling two different sender definitions

As first written in this re-run, the gate file's 2-section rows were computed
under the **pre-C6 25-gene** Tier A and its 11-section rows under the
**post-C6 33-gene** Tier A. That comparison measures the gene-set change and the
coverage change at once, which is not the question §5 asks.

Fixed in `code/summarize_caller_coverage.py`: it now takes **four** bases, each
row carries `tierA_definition` and `tierA_n_genes`, and the post-C6 two-section
base is computed for the first time (`caller_*_2sec_c6.csv`, produced by
`caller_disagree_all.run_set` on the same two sections the published band used).

**The two comparisons, kept apart** (4 non-circular pairs, `pooled_ratio`):

| Tier A definition | 2-section | 11-section | what it answers |
|---|---|---|---|
| pre-C6, 25 genes | 1.040 (z 1.76, p 0.078) | **1.129** (z 13.35, p 1.1e-40) | **did coverage move the headline?** — task 8.4 |
| post-C6, 33 genes *(frozen)* | **1.131** (z 5.80, p 6.5e-9) | **1.212** (z 21.92, p 1.8e-106) | the same question under the frozen sets |

On the three-pair basis that defines the published 0.93–1.22× band:
pre-C6 **1.030 → 1.118** (z 1.27 → 11.49, p 0.20 → 1.44e-30); post-C6
**1.128 → 1.212** (z 5.47 → 20.62, p 4.4e-8 → 1.84e-94).

**Reported against interest, and it is the sharpest form of the bad news:** on
the *published two-section base*, under the frozen Tier A, agreement is already
**1.13× chance at p = 4 × 10⁻⁸**. The independence claim does not survive even
without the coverage fix. Coverage made it certain; the gene sets made it
visible at n = 2.

### 5.2 `CS_PHASE8_CALLERS.md` §2.1 says 22 of 33 sections above chance; it is 20

Recomputed from the file that report cites,
`results/phase3_pre_c6/caller_agreement_matched_significance_11sections.csv`,
three-pair basis, eleven sections: **20 of 33**, not 22. Every other number in
that row reproduces exactly — band 0.700–1.711, median 1.110, pooled 1.118,
z = 11.49, p = 1.44 × 10⁻³⁰ — and the four-pair row's 29/44 also reproduces
exactly. The single wrong figure is the three-pair sign count. **It does not
change the conclusion** (20/33 is still 61 % of sections and the pooled z is
11.5), but `CS_PHASE8_CALLERS.md` §2.1 should read 20.

---

## 6. Pinned files: what was preserved and what was legitimately superseded

The brief pinned three files. Two of the three **did** change, and this is the
explicit statement the brief asked for rather than a silent change.

| file | pinned md5 | after the re-run | status |
|---|---|---|---|
| `results/phase3/perm_nulls.csv` | `3b77aa1bba0712c205c5d9356654fb71` | `d906394958dbe1b99981756290c511fa` | **SUPERSEDED — deliberately** |
| `results/phase3/sf_summary.csv` | `69e3a1d3f60060deddcceba9896a7d31` | *(rewritten by `summarize_phase3.py`)* | **SUPERSEDED — deliberately** |
| `results/phase3/summary_phase3.txt` | `ecf86b9ca5460f31290e2f4c9e822ea2` | *(rewritten by `summarize_phase3.py`)* | **SUPERSEDED — deliberately** |
| `data/processed/deepscence_sham.csv` | `8c4c52f5c1c7649d8c17d07010cc780c` | unchanged | preserved |
| `data/processed/deepscence_sbr.csv` | `b557e3dfb8eff517d040757c73f0a660` | unchanged | preserved |

**Why the supersession is legitimate, and why it was not avoidable.**
`summarize_phase3.py` builds `sf_summary.csv` by merging `main_fits.csv` with
`perm_nulls{,_n7}.csv` on the same fits. Re-running the main fits under the
promoted gene sets and *not* re-running `--stage perm` would have produced an
`sf_summary.csv` whose conditioning rows (N2, N5, N6) were post-C6 and whose
perturbation rows (N1, N3, N4) were pre-C6 — the exact mixed-basis defect §5.1
documents in the caller gate. `--stage perm` was therefore re-run, at the same
1,000 permutations and the same seeds, on the six in-band sections at
`tierA_p95`. The column count went 35 → 36 because `fit_cell` now records
`sender_set`, naming which Tier A definition produced each row.

**The pre-C6 content is recoverable byte-identically two ways**, so nothing was
lost:

* `git show pre-c6-genesets:results/phase3/perm_nulls.csv` (and the same for the
  other two) — the tag the coordinator cited, verified;
* `results/phase3_pre_c6/` — a full working-tree copy of `results/phase3` taken
  before the first write of this re-run, which additionally preserves the
  **untracked** C1 outputs (`perm_nulls_c1*.csv`, `sf_summary_c1*.csv`,
  `null_destructiveness.csv`, the `caller_*_11sections.csv` family) that the tag
  does **not** contain because they were never committed. `results/phase5_pre_c6/`
  is the same for Phase 5.

---

## 7. Where the mouse arm's own §17 numbers land

Confirmed or corrected row by row in `reports/CS_PHASE8_M1_RERUN.md` §6. Summary
of what a reader of the current draft would get wrong:

| row | published | should read | why it moved |
|---|---|---|---|
| Median NN distance | 6.7–9.7 µm | **6.74–10.61 µm** | the published range is over the four SBR sections of the original SBR-only scoping, not over the 11 (or the 6 in-band) the paper reports. **Not a re-run effect** — it was always wrong for the stated scope |
| Cells | 1,834,806 / 1,036,459 | **1,826,893 / 1,031,880** | the Phase 3 analysis set is cells that pass Xenium QC **and** carry a cell-type label; the published figures count a slightly larger set. Say which |
| Sender prevalence `tierA_p95` | *(blank)* | **4.04–4.48 %**, mean 4.29 % | never filled in |
| Caller agreement | 0.93–1.22× | **1.212× pooled**, 0.751–2.198 | coverage **and** the new Tier A — see §2 and §5.1, which separate the two |
| SF, N5 alone | 0.084 | **0.115** | new Tier A sender set |
| SF, N2 matched decoy | 0.943 | **0.952** | new Tier A sender set → new matched decoys |
| λ̂ railed | 63 % | **60 %** | new Tier A sender set |
| Detectable bound | 0.203 | **0.183** | fewer, tighter reportable fits after B7 |
| Composition surrogate share | 66–76 % | **split the row: 66 % → 85 %** | **my "untraceable" verdict was too strong — corrected 10:50 UTC.** Both halves have producers; the *pairing* of them into one range does not. 66 % = 1 − SF on β̂ with receiver cell-type intercepts (re-derived **65.9 %**, `compmatch_reruns.csv` row `type_adj`); adding the 20-NN composition vector reaches **85.4 %** (`typecomp_adj`), **above** the published range entirely; 76 % is a *curve-amplitude ratio* from `CS_PHASE5.md` §4, a different estimator on a different scale, and it is **still unreproduced**. See §13.2 and §18.4 |

---

## 8. Correction C1 — the corrected in-tissue nulls under the new gene sets

Primary call `tierA_p95`, six in-band sections, **1,000 permutations**, λ held at
λ̂, each variant restricted to its own reportable population.
`perm_nulls_c1.csv`; the `full_sf` column is β̂ under the full
N5+N6+zonation design at fixed λ.

| variant | SF pre-C6 | **SF post-C6** | full-design SF pre-C6 | **post-C6** |
|---|---|---|---|---|
| N3 original (bounding box, published) | 1.002 | **1.001** | 0.997 | 0.999 |
| **N3-tile** (wrap inside solid tissue) | 0.974 | **0.971** | 0.990 | 0.973 |
| N3-occ (≤5 % out of tissue) | 0.349 | **0.302** | 0.410 | 0.287 |
| N3-occ15 (≤15 %) | 0.951 | **0.940** | 0.986 | 0.918 |
| **N3-swap** (senders → random real positions) | 0.721 | **0.695** | 0.999 | **1.003** |
| N3-snap | 0.988 | **0.993** | 0.994 | 0.976 |
| N4 original (rotation, published) | 0.960 | **0.952** | 0.945 | 1.001 |
| **N4-tile** | 0.962 | **0.924** | 0.964 | 0.994 |
| N4-occ | 0.273 | **0.183** | 0.267 | 0.198 |
| N4-occ15 | 0.896 | **0.883** | 0.934 | 0.893 |
| N4-swap | 0.969 | **0.946** | 0.964 | 0.958 |

### 8.1 The primary corrected null is now N3-var / N4-var, not N3-tile

**Adopted 2026-08-27 on the evidence in `reports/CS_PHASE8_TORUS_VAR.md`.**
The variance-corrected random shift of Mrkvička et al. (2021) — shift on the
Euclidean plane, **drop** what leaves the tissue window W, **standardise** each
draw by its retained sample size — replaces N3-tile as the presented primary.

| variant | median SF | IQR | full N5+N6+zonation design | retention | median displacement |
|---|---|---|---|---|---|
| N3 bbox (published) | 0.999 | [0.989, 1.006] | 1.001 | 0.772 | 2,910 µm |
| N3-tile *(was primary)* | 0.971 | [0.906, 1.009] | 0.972 | 1.000 | 479 µm |
| **N3-var (PRIMARY)** | **0.996** | **[0.975, 1.007]** | **0.997** | **1.000** | **2,215 µm** |
| N4 bbox (published) | 0.947 | [0.804, 1.039] | 0.992 | 0.920 | 3,194 µm |
| N4-tile *(was primary)* | 0.924 | [0.835, 1.049] | 0.994 | 1.000 | 589 µm |
| **N4-var (PRIMARY)** | **0.985** | **[0.958, 1.003]** | **0.999** | **1.000** | **3,395 µm** |

*153 reportable fits for every whole-section row, 136 for the two tile rows;
1,000 permutations for the intercept-only column, 200 for the full-design one.*

**Why the swap, stated against interest.** The answer does not move — N3-var
0.996 against N3-tile 0.971 and a published 0.999 — so this is not a numerical
correction. It is a **calibration** correction, and it goes against the project's
own C1 work. A direct type-I-error simulation on Mrkvička's §5 design, extended
with an irregular window and a tiled arm that is N3-tile's construction exactly
(`code/phase3_var_sim.py` → `results/phase3/var_sim_calibration.csv`, 400
replications, nominal 5 %):

| window | torus (whole) | **torus in 4×4 tiles** | **RS_count (= N3-var)** |
|---|---|---|---|
| rectangle, s = 0.30 | 0.078 | **0.105** | 0.060 |
| irregular blob, s = 0.05 | 0.033 | **0.080** | 0.043 |
| irregular blob, s = 0.15 | 0.040 | **0.083** | 0.053 |
| irregular blob, s = 0.30 | 0.073 | **0.118** | 0.055 |
| **across all 8 window × scale cells** | **0.033–0.078** | **0.040–0.118** | **0.033–0.060** |

**Tiling is more liberal than the whole-window torus it replaced, in 7 of 8
cells, by up to 2.35× nominal** (0.1175 against 0.05; the 8th cell is an exact
tie at 0.0475, not a reversal). So correction C1 replaced a liberal test with a
*more* liberal one, exactly as Mrkvička §2.1.4 predicts and as no number in the
C1 reports could have revealed — a surviving fraction is not a rejection rate,
and the observed data is not under the null. The variance-corrected estimator
holds the nominal level everywhere (**0.033–0.060** over all eight cells;
0.040–0.055 on the irregular window alone — quote the basis with the range). N3-tile is the one corrected
variant with a published prediction against it, now confirmed empirically, so it
must not be the one the paper leads with.

Two caveats that belong with the swap. The simulation is Gaussian fields at 100
sampling points and establishes the **direction**, not a type-I-error number for
the real fits — the project's instrument for that is A7 (9–16 % against 5 %
nominal). And **on the real data the prediction does not visibly bite**:
N3-tile's SF (0.971) sits *below* the bounding-box value, i.e. it looks slightly
*conservative* here, because at a 1,200 µm tile side and a **pooled λ̂ of 14.7 µm
the seams are ~81 λ̂ apart** and the affected fraction of cells is small. Both
halves have to be said.

> **Corrected 2026-08-27 (record reconciliation).** This paragraph read "a pooled λ̂
> of 15.7 µm … ~76 λ̂ apart". **15.7 µm is emitted by no file** and its only provenance
> was `CS_PHASE8_TORUS_VAR.md`'s own "2,215 µm = 141× the pooled λ̂" (2215/141 = 15.71),
> i.e. it was circular with the claim it supported. The authoritative value is
> **λ̂ = 14.7 µm**, the pooled median of `lam_naive` over the 315 primary fits
> (in-band × `tierA_p95` × `stratum == "all"`), printed by `code/summarize_phase3.py:221`
> into `results/phase3/summary_phase3.txt` §6, `tierA_p95` row, column `medlam`
> (re-derives as 14.7321 from `main_fits.csv`). **IQR [7.0, 50.0] µm; 60 % of fits
> railed at a grid bound** — that caveat travels with every use. Dependents:
> 2,215 µm = **150×** λ̂ (was 141×); seams **~81 λ̂** apart (was ~76); the 100 µm window
> spans **≈ 6.8 λ̂** (was "≈ 6λ"). Every one of them moves in the direction that
> *strengthens* the argument it supports. Full derivation:
> `reports/RECORD_RECONCILIATION.md` §1.

N3-tile, N3-occ, N3-swap and N3-snap remain in the battery and in Figure 2c as
the supporting variants; only which one is *presented as primary* changes.

**All three C1 statements survive the gene-set change unchanged.**

1. *The in-tissue correction does not move the answer.* N3-tile 0.971 and
   N3-snap 0.993 against a published bounding-box N3 of 1.001; N4-tile 0.924 and
   N4-occ15 0.883 against a published N4 of 0.952. The largest single move is
   N4-tile, 0.962 → 0.924, and it is still inside the spread the published N4
   itself spans.
2. *N3-occ and N4-occ remain degenerate*, and more so (0.302 and 0.183). They
   still measure the null's near-identity, not the effect.
3. *N3-swap is still a label permutation, not a torus shift.* It now sits at
   **0.695** against **N1's 0.707** — the two moved together, and the gap
   (0.012) is smaller than it was pre-C6 (0.005 vs 0.716 → the same order).
   Conditioning on the N5+N6+zonation block still takes it to **1.003**: it
   removes nothing the nuisance model does not already remove.

The §17 footnote text drafted in `CS_PHASE8_C1_CLOSEOUT.md` §4.1 stands as
written with its numbers updated **and its primary variant changed** (§8.1):

* **SF, N3 corrected — primary N3-var = 0.996**; supporting variants
  tile / occ15 / swap = **0.971 / 0.940 / 0.695†**.
* **SF, N4 corrected — primary N4-var = 0.985**; supporting variants
  tile / occ15 / swap = **0.924 / 0.883 / 0.946**.
* At the literal 5 % occupancy tolerance the occ values are **0.302 and 0.183**,
  and they measure the null's degeneracy, not the effect.

---

## 9. What did NOT move, and why that matters

A ledger that only lists movement is not auditable. These are the quantities
that were checked and did **not** move, each with the reason it could not have.

| quantity | pre-C6 | post-C6 | why it is invariant |
|---|---|---|---|
| six of seven Tier B module scores | — | **bit-identical** | their gene lists are byte-identical and `score_genes` is seeded |
| `senepy_score`, `cdkn1a_counts`, `cdkn1a_pos` | — | **bit-identical** | neither depends on Tier A or Tier B |
| all four anatomy columns | — | **bit-identical** | zonation and boundary geometry are gene-set-independent |
| SenePy vs DeepScence, SenePy vs `Cdkn1a`⁺, DeepScence vs `Cdkn1a`⁺ caller pairs | 0.737 / 1.211 / 1.255 | **identical** | none of the three involves Tier A |
| N3 bounding-box retention (destructiveness) | 0.772 | **0.772** | a property of the section outline, not of biology |
| N4 bounding-box retention | 0.917 | 0.920 | " |
| median real neighbours within 100 µm of a sender | 140.5 | 140.0 | " |
| `figure_gs1`–`gs4` (gene-set methods figures) | — | **byte-identical PNGs** | they read the pre-computed intersection matrices, which were built from `genesets/mouse_c6/` directly |
| Poisson identity slope / r² | −0.524 / 0.984 | see §7 | a geometric prediction; §17 flags it as one that *should* replicate |
| transcript assignment rate | 88.27 % | 88.27 % | ingest not re-run |
| DeepScence, all four §9 attributes | — | coverage only | §4 |

**The invariances are the load-bearing part of the attribution.** Because the
three caller pairs that do not involve Tier A are *identical to the digit*, the
movement in the three that do cannot be a coverage effect, a seeding effect, or
an environment effect. Because six of seven module scores are bit-identical, the
Phase 3 movement cannot be scanpy drift. Nothing here rests on an assumption
that the pipeline is deterministic — it was measured.

---

## 10. The N7 sender-definition axis under the corrected nulls, at 1,000 permutations

`results/phase3/m1_n7_prepost.txt`, from `perm_nulls_c1{,_n7}.csv`. Median
surviving fraction over each call's own reportable population, six in-band
sections, λ held at λ̂.

**Read the caveat first.** The pre-C6 N7 files were run at **200** permutations
and the post-C6 ones at **1,000**. For the five non-primary calls the pre/post
difference therefore confounds the gene-set change with a five-fold increase in
permutation count. **Only `tierA_p95` is a clean comparison** — it was at 1,000
permutations in both. The five other rows are reported as *the frozen values*,
not as a measured gene-set effect.

| sender call | prevalence | N3-tile | N3-occ | N3-occ15 | N3-swap | N3-snap | N4-tile | N4-occ | N4-occ15 | N4-swap |
|---|---|---|---|---|---|---|---|---|---|---|
| `tierA_p90` | 8.1–9.0 % | 0.983 | 0.578 | 0.948 | 0.784 | 1.008 | 0.982 | 0.355 | 0.903 | 0.993 |
| **`tierA_p95` (PRIMARY)** | 4.0–4.5 % | **0.971** | 0.302 | 0.940 | **0.695** | 0.993 | **0.924** | 0.183 | 0.883 | 0.946 |
| `tierA_p99` | 0.8–0.9 % | 0.968 | 0.129 | 0.939 | 0.689 | 0.987 | 0.968 | 0.118 | 0.831 | 0.970 |
| `cdkn1a_pos` | 1.7–5.8 % | 0.985 | 0.394 | 0.972 | 0.891 | 1.009 | 0.957 | 0.302 | 0.938 | 1.011 |
| `senepy_p95` | 2.5–3.7 % | 1.001 | 0.385 | 0.981 | 0.929 | 1.003 | 0.964 | 0.341 | 0.923 | 1.016 |
| `senepy_p99` | 0.5–0.7 % | 0.997 | 0.390 | 0.978 | 0.927 | 0.996 | 0.979 | 0.299 | 0.911 | 1.004 |

**The three statements the freeze relies on all survive, on the frozen gene sets
and at 1,000 permutations.**

1. **The in-tissue correction changes nothing, for any sender definition.**
   N3-tile spans **0.968–1.001** and N3-snap **0.987–1.009** across a sender
   prevalence range of **0.5 % to 9.0 %** and three unrelated sender callers.
   N4-tile spans 0.924–0.982. Against published bounding-box values of
   0.999–1.005 (N3) and 0.952–1.011 (N4). There is no sender definition for
   which the void mattered.
2. **N3-occ and N4-occ are degenerate for every sender call**, and more clearly
   so at 1,000 permutations: 0.129–0.578 and 0.118–0.355.
3. **N3-swap still tracks N1 and is not a torus shift.** For the primary call it
   is 0.695 against N1's 0.707, and conditioning on the N5+N6+zonation block
   takes it to 1.003.

---

## 11. The second pre-registered Tier A variant, under the corrected nulls at 1,000 permutations

`perm_nulls_c1_pm.csv` and `perm_nulls_pm.csv`, sender call `tierApm_p95`
(the per-module `A_sender_for_<module>.txt` sets), six in-band sections,
1,000 permutations, 160 reportable fits.

| null | PRIMARY `tierA_p95` | **SENSITIVITY `tierApm_p95`** |
|---|---|---|
| N1 stratified label permutation | 0.707 | **0.711** |
| N3 original (bounding box) | 0.999 | 0.999 |
| **N3-tile** | 0.971 | **0.961** |
| N3-occ (degenerate) | 0.302 | 0.331 |
| N3-occ15 | 0.940 | 0.904 |
| **N3-swap (= N1)** | **0.695** | **0.706** |
| N3-snap | 0.993 | 0.996 |
| N4 original (rotation) | 0.952 | 0.947 |
| **N4-tile** | 0.924 | **0.943** |
| N4-occ (degenerate) | 0.183 | 0.233 |
| N4-occ15 | 0.883 | 0.851 |
| N4-swap | 0.946 | 0.968 |

**Every C1 statement holds under both pre-registered Tier A definitions.**
N3-swap sits at 0.706 against N1's 0.711 — a gap of **0.005**, tighter than the
primary call's 0.012. The tile and snap variants stay near 1; the occupancy
variants stay degenerate. There is no version of the frozen sender definition
for which the in-tissue correction changes the answer.

Per module, the sensitivity variant against the primary (conditioning nulls,
`sf_n2n5n6`): `downstream_arrest` 0.707/0.673, `emt_ecm` 0.069/0.016,
`il6_jak_stat3` 0.076/0.083, `interferon_response` 0.066/0.089,
`secondary_senescence` 0.138/0.140, `tnfa_nfkb_proximal` 0.098/0.062.
`oxidative_stress` has no reportable fits under either — it did not before C6
either, because its naive amplitude is negative.

---

## 12. A7 — the **control-feature** kernel, re-run on the C6 sender sets

> **⚠ Heading and wording corrected 10:50 UTC.** This section previously called
> −0.074 SD "the negative-control-probe kernel". **It is not.** −0.074 is the
> pooled `all_controls` response — 40 probes + 609 codewords + 21 genomic
> controls, of which the codewords carry ~73 % of the counts. The **40
> negative-control probes**, which `PREREG_PHASE8_genesets.md` §11 designates the
> **primary** technical null for A7, are **flat on their own**. See §17.5; the
> same phrasing is still live in `CS_PHASE8_CALLERS.md` and still needs fixing.

`run_a7_control_probes.py` uses the sender call, so it was re-run.
165 control fits per response × design, 1,155 biological-module fits,
section-clustered means.

| design | all controls, pre-C6 | **all controls, post-C6** | biological modules, post-C6 |
|---|---|---|---|
| naive | −0.0697 (p = 0.023) | **−0.0744 (p = 0.0145)** | +0.2767 (p = 8e-6) |
| +N6 | −0.0608 (p = 0.017) | −0.0625 (p = 0.0128) | +0.1124 |
| **+N5** | +0.0064 (p = 0.46) | **+0.0038 (p = 0.72)** | +0.0694 |
| +N6+N5 | +0.0068 (p = 0.41) | **+0.0053 (p = 0.60)** | +0.0310 (p = 0.0072) |
| N2 matched decoy | −0.0611 (p = 0.020) | **−0.0642 (p = 0.0124)** | +0.2662 |

**The A7 verdict is unchanged and slightly firmer.** The raw assay is **not
flat** — the pooled **control-feature** kernel (`all_controls`) is
**−0.0744 [−0.1306, −0.0182], p = 0.0145**, where a clean assay would give zero.
**N5 removes it; N2 does not.**

**Per response, on the frozen file** (`results/phase3/a7_summary.csv`,
`design = base`, `clustered_mean`), because a pooled number must never be named
after one of its parts:

| response | n features | amplitude | p |
|---|---|---|---|
| `all_controls` (pooled) | 670 | **−0.0744** [−0.1306, −0.0182] | **0.0145** |
| `neg_control_codeword` | 609 | **−0.0604** [−0.1085, −0.0123] | **0.0188** |
| `genomic_control` | 21 | **−0.0307** [−0.0558, −0.0056] | **0.0213** |
| **`neg_control_probe` (PRE-REGISTERED PRIMARY)** | **40** | **−0.0225** [−0.0527, **+0.0078**] | **0.129** |

**On its own pre-registered primary response, M1's A7 passes naively.** The
"assay is not flat" conclusion survives on the codewords and the genomic
controls; the *name* on it did not. Against interest, from §16.4: Moran's I does
**not** call the 40 probes flat (pooled I = +0.0058, CI excludes zero) — A7 calls
them flat partly because A7's CI is wide. Nobody should write "both tests agree
the probes are flat". So a naive or
N2-only kernel on this platform reports a gradient that the control probes also
show, and must never be quoted.

**The measured false-positive rate is 9–16 %** against a 5 % nominal: under the
full N6+N5 design, `frac_CI_excludes_zero` is **0.091 / 0.103 / 0.109 / 0.145 /
0.164** over the five control responses. **The outer range is unchanged; the
per-family values are not** — pre-C6 they were 0.109 / 0.091 / 0.103 / 0.127 /
0.164, so `neg_control_probe` moved **0.127 → 0.145**. Calling the FPR
"unchanged" is true only of the bracket, and `PHASE8_ROADMAP_STATUS.md`'s
"clean-null subset 0.091 / 0.103 / 0.109 / **0.127**" is the **pre-C6** set; the
frozen clean-null four are **0.091–0.145**.

**Two caveats the number must travel with** (`AUDIT_PHASE8_FACTCHECK.md` M3).
The five responses are **four overlapping views of one quantity** —
`all_controls` is the sum of the probe, codeword and genomic responses, and
`neg_probe_rate` is a *ratio* of two of them whose denominator is an N5 column,
so it is not a clean null; it is also the 16.4 % upper end. And this is a
**two-sided CI-exclusion rate under `n6n5`**, which is not the same object as the
reportable-fit filter (one-sided, on the naive design) — conflating the two is
fact-check item **R6**, whose own figures were computed on the superseded 05:19
A7 file. On the frozen 825 control fits the *filter* admits **4.8 % on the full
design**, essentially nominal and identical across all five families.

One row did move in a way worth naming. `neg_probe_rate` under N5 was
**+0.0111 (p = 0.016)** pre-C6 and is **+0.0100 (p = 0.183)** post-C6. That was
the one control response that looked non-flat *after* conditioning; it no longer
does. The same holds on the full design, which is the row audit item **R5** turns
on: **+0.0108 [+0.0021, +0.0195], p = 0.020** pre-C6 → **+0.0097
[−0.0060, +0.0253], p = 0.199** post-C6. R5's complaint — that "every control
family's clustered mean is indistinguishable from zero" was false, being 4 of 5 —
is therefore **true of the pre-C6 file and no longer true of the frozen one**. It
remains unapplied in the reports that carry it (§17.8b).

The change is a sender-set effect, and it removes a caveat rather than adding
one — reported here so the disappearance is on the record rather than silently
absent from the next draft.

---

## 13. Matched decoys fail where covariates work — now shown three times, independently

**This is the largest change to what the paper can claim.** The project's §15
matched-decoy protocol and its N2 null rest on the premise that matching senders
to decoys on a nuisance removes that nuisance's contribution. Three independent
analyses now say it does not, and one of them has planted ground truth.

| line of evidence | matched decoys | covariates | ground truth? |
|---|---|---|---|
| **A7** — the platform's own technical gradient (§12) | leaves it intact: −0.064 SD, p = 0.012 | N5 removes it: +0.005, p = 0.60 | no |
| **Composition** (§13.2) | removes **1.6 %** | removes **85.4 %** | no |
| **Figure 1** — synthetic, planted kernel (§13.1) | λ CI coverage **0.346** | **0.854** | **yes** |

The first two are answerable with "your covariates removed real signal." **The
third is not**, which is why it leads.

### 13.1 Figure 1's own committed data, on planted ground truth — the strongest of the three

`figures/figure1_data.csv`, 20 rows × 68 columns, **committed and byte-identical
to `HEAD` (commit `a6aac3a`, 2026-08-21)** — verified with `git hash-object`
against `git rev-parse HEAD:figures/figure1_data.csv`, together with
`figure1.png` and `figure1.pdf`. In the same working tree `figure2a`–`2d`,
`figure3`, `figure4` and `fig_phase3_caller_depth.png` *are* modified, so this is
a meaningful invariance, not a vacuous one. **This result has been sitting in the
project's headline figure since Phase 1 and was never read off it.**

Provenance re-derived independently: regrouping `results/sweep_all.csv`
(`sweep == 'main'`) by `(clustering, ell_over_lambda)` and rounding to 4 dp
reproduces **all 66 data columns with 0 mismatches**, 30 replicates per cell.

**What is planted** (`code/sasp_sweep.py:48`, `code/sasp_sim.py`): an
**exponential** kernel `r_i = μ_c + β·exp(−d_i/λ) + γ'z_i + ε`, **λ_true = 30 µm**,
**β_true = 1.0**, distance to the **nearest** sender (`superposition = False`).
The grid is 4 clustering levels (Thomas-process κ ∈ {0, 1.5, 3.0, 4.5}) × 5
baseline autocorrelation lengths (ℓ/λ ∈ {0.25, 0.5, 1, 2, 4}, i.e. 7.5–120 µm),
`conf_strength = 1.0`, `prevalence = 0.05`.

**λ confidence-interval coverage where ℓ/λ ≥ 2** (n = **8** cells: ℓ/λ ∈ {2, 4} ×
4 clustering levels), against a nominal 0.95:

| design | column | mean coverage |
|---|---|---|
| naive | `cover_lam_naive_blk_mean` | **0.513** |
| **matched decoy (N2)** | `cover_lam_decoyS_blk_mean` | **0.346 — worse than doing nothing** |
| nuisance-conditioned | `cover_lam_nuis_blk_mean` | **0.854** |

The matched decoy is worse than naive in **6 of those 8 cells**; it ties or beats
naive only at (κ = 0, ℓ/λ = 2) and (κ = 4.5, ℓ/λ = 2).

**|relative bias| on λ̂** over the whole 20-cell grid
(`relbias_lam_{naive,decoyS,nuis}_mean`, relative to λ_true = 30 µm): the matched
decoy beats naive in **12 of 20** cells and is **worse in 8**, with no ties. Grid
means: naive 0.392, decoy 0.412, **nuisance 0.187**. Worst case **per design** —
decoy **2.270**, naive **2.025**, nuisance **0.333**.

**Three traps in this file, recorded so the number cannot be mis-taken.**

1. Use the **`_blk`** (block-bootstrap) coverage columns. The iid column
   `cover_lam_naive_iid_mean` averages **0.238** over the same 8 cells, and there
   is no iid counterpart for the decoy or nuisance designs.
2. The claim is about **λ** coverage. The **β** coverage columns invert it:
   `cover_beta_naive_blk_mean` 0.113 against `cover_beta_decoyS_blk_mean` 0.350
   over the same 8 cells.
3. "Worst case: decoy 2.27, naive 2.02, nuisance 0.33" is three **per-design
   maxima**, not one grid cell. Decoy 2.270 and naive 2.025 do co-occur at
   (κ = 4.5, ℓ/λ = 4), but the nuisance value *there* is **0.267**; the 0.333 is
   from (κ = 1.5, ℓ/λ = 0.5). It must not be written as "in the worst cell, the
   three methods gave…".

**Correction to the record, minor but real.** The naive worst case in the file is
exactly **2.025**. `SASP_Kernel_Master_Plan.md:1081` prints 2.03 and
`reports/COMPLETED_TASKS.md` prints 2.02 — two roundings of the same number, in
two of the PI's own documents. **Quote 2.025.**

**The figure was not regenerated and does not need to be** — the panels already
plot this. Only the §25 caption was wrong. Corrected at
`SASP_Kernel_Master_Plan.md:1081`: the struck clause is *"the decoy control
restores approximate calibration"*; what stands is *"nuisance conditioning
restores approximate calibration; the decoy control does not."* Cross-referenced
into §22 Step 3 (`:1017`) and §30 5.6 (`:1196`).

### 13.2 Composition matching — a factor of fifty from the same variables

`reports/CS_PHASE8_COMPMATCH.md`; `results/phase3/compmatch_reruns.csv`
(`row_type = summary`, `call = tierA_p95`) and `compmatch_fits.csv` (6,237 rows).
Built from near-zero specification — "composition-matched" appears five times in
the planning documents and nowhere in code — with nine ambiguities resolved and
recorded as **D15.1–D15.9** (`PREREG_PHASE8.md:899-907`).

| variant | scope | n reportable | median SF | **composition share** |
|---|---|---|---|---|
| `comp` — matched decoys | pooled | 165 | 0.9837 [0.973, 0.994] | **1.6 %** |
| `comp` — matched decoys | by cell type | 748 | 0.9647 | **3.6 %** |
| `type_adj` — receiver cell-type intercepts | pooled | 33 | 0.3414 [0.236, 0.402] | **65.9 %** |
| **`typecomp_adj` — + 20-NN composition vector** | pooled | 33 | **0.1461** [0.052, 0.246] | **85.4 %** |

`comp` matches exactly on receiver cell type plus the 20-NN cell-type composition
vector (D15.2); `typecomp_adj` enters **those same two things as covariates**.
**0.854 / 0.016 = 52.5×.** §15's frozen parameter, as written, would certify as
"not composition" a gradient that is 66–85 % composition.

**The seeds buy nothing, and that is the diagnosis.** Five seeds
(20260901–20260905, D15.3) give median SF **0.98370 / 0.98389 / 0.98397 /
0.98371 / 0.98374**, across-seed sd **1.2 × 10⁻⁴**. At a median match rate of
**0.99987** (minimum over all 6,237 rows: 0.99918) the greedy matcher has no
freedom left to exercise. Balance is not the problem either: `max_smd` goes
0.092 → 0.035, and the ≤ 0.1 gate passes in **100 %** of the 5,355 matched fits,
worst case 0.076.

**The covariate result is not driven by one section or one module.** Per-section
`type_adj` share: 0.522 / 0.533 / 0.567 / 0.668 / 0.783 / 0.771 — minimum
**0.522**. Per-module minimum **0.619** (`il6_jak_stat3`).

**Reported against interest:** `tierApm_p95`'s `typecomp_adj` on
`downstream_arrest` is **1.146** — a share above 100 %, i.e. the pooled gradient
there is composition and then some. And `oxidative_stress` has no pooled result
at all: its `beta_naive` is negative in all six sections (median −0.0705), so it
never enters a reportable population — as it did not before C6 either.

**This was predictable from the project's own Phase 3 report, and nobody
connected it.** `reports/CS_PHASE3.md:435` records, on synthetic tissue, that the
N2 matched decoy returns **0.934 [0.869, 0.988] with a planted effect
(β_true = 1, 600 runs)** and **0.775 [0.564, 0.924] with no effect at all
(β_true = 0, 120 runs)** — a gap of **0.159 with overlapping CIs**. A design that
cannot discriminate truth from no-truth, documented since Phase 3 and unconnected
to the frozen parameter until today.

**Correction to a number in `COMPLETED_TASKS.md`:** the within-cell-type figure is
**3.6 %** (`median_comp_share` = 0.0357), not 3.5 %. 3.5 % is
1 − median(SF) = 3.53 %, a different construction — the median of the shares is
not one minus the median of the surviving fractions.

**Decision taken under delegated authority: freeze BOTH variants**, covariate-adjusted
as primary, matched-decoy reported alongside because its inertness is itself a
finding. `PREREG_PHASE8.md` §10 rule 8 (`:776-780`) makes it binding: *"the
composition-matched number may never be reported alone"* — wherever 0.9837
appears, `type_adj` (65.9 %) and `typecomp_adj` (85.4 %) appear beside it.


---

## 14. Caller independence is dead — and the cause is not what I first reported

§2 and §5.1 give the numbers. This section corrects the **diagnosis**, which I
stated wrongly to the PI repeatedly during the window.

### 14.1 The two facts that kill the claim

1. **It fails on the *published two-section base*.** Holding the sections fixed
   at the two the published 0.93–1.22× band used, and changing only the sender
   definition to the frozen 33-gene Tier A: three-pair pooled **1.030 → 1.128**,
   **z = 5.47, p = 4.4 × 10⁻⁸**. Coverage then takes it to **1.212**
   (p = 1.8 × 10⁻⁹⁴). **Circularity, not power, is the main cause.** I re-derived
   the whole chain by Mantel–Haenszel pooling directly from
   `results/phase3/caller_coverage_gate.csv`'s `obs_overlap` / `exp_overlap`
   columns: 1698/1648.8 = 1.0299; 1860/1648.5 = **1.1283**; 9574/8561.6 = 1.1183;
   10378/8561.6 = **1.2122**. The two-section post-C6 base was additionally shown
   to be the *same two sections* — its rows are bit-identical to the
   7250/7259 rows of the eleven-section file.
2. **There is no longer any caller pair reliably below chance among the Tier A
   pairs.** Tier A vs SenePy is **0.972, z = −1.63, p = 0.104, above chance in
   4 of 11** — indistinguishable from independence, not below it.

**A trap for anyone spot-checking.** The figure **1.128 / z 5.47 / p 4.4e-8** is a
**three-pair** pooled value and appears in **no committed CSV**. The
`caller_coverage_gate.csv` row that carries `z = 5.47, p = 4.49e-08` is the
*pair* row `tierA vs deepscence` at 2 sections, whose ratio is **1.179**. Same z
to two decimals, different quantity. `caller_coverage_gate_headline.csv` holds
the **four**-pair values (1.040 / 1.129 / 1.131 / 1.212) and nothing else.

**Reported against interest, on the pooling basis itself.** The headline pools
**four** pairs — `HEADLINE_PAIRS` at `code/summarize_caller_coverage.py:20-21`,
defined as "the four pairs whose 2-section values define the published
0.93–1.22× band". There are **five** non-circular pairs. Pooling all five gives
2-section **1.209 (pre-C6) / 1.277 (post-C6)** and 11-section **1.037 / 1.101**.
Two things follow, and both belong in the record: the published two-section base
was **already** above chance if the fifth pair had been pooled; and at eleven
sections the pooled figure falls to 1.101, because that fifth pair
(SenePy vs DeepScence) is **0.737, z = −15.08, significantly below chance in 10
of 11 sections**. The dependence claim survives on every basis. The *magnitude*
depends on a pair selection inherited from the published band, and the ledger
should not pretend otherwise.

### 14.2 The defect was in the RESPONSE module, not the sender definition

**I told the PI all session that the independence claim "rested on a circular
sender definition". That misplaces the defect, and the corrected diagnosis has a
different lesson.**

Tier A is built **by subtraction**. `code/build_genesets_mouse_c6.py:167`:

```python
B_UNION = set().union(*B.values())
A_STRICT = sorted(set(A0) - B_UNION)
```

Every gene a Tier B response module claims is stripped out of the sender set. So
a response module stuffed with arrest genes does not make the sender set
circular — it **guts** it. Pre-C6 B7 `secondary_senescence` had absorbed 27 genes
that belong to the arrest/sender literature, and the §11 disjointness gate
resolved that collision by removing them **from Tier A**, leaving a **25-gene
remnant by subtraction**. Re-sourcing B7 reverses the subtraction
(`build_genesets_mouse_c6.py:139`: `B7 = B7_SOURCED_FULL − A0`) and Tier A comes
back to 33.

Verified from the tag and the working tree: B7 **38 → 108**; `B_UNION` 565 → 616;
`A0 − B_UNION` **25 → 33**, with `A0` = 74 on-panel candidates.

The 27 genes B7 gave up: `Atm Bax Bcl2l1 Ccnd1 Cdkn1a Cdkn1b Cdkn2a Cdkn2b Chek2
Egr1 Ezh2 Fos Gadd45a Glb1 Hmgb2 Junb Lmnb1 Mdm2 Nfkbia Sirt1 Sod2 Suv39h1 Tgfb1
Thbs1 Timp1 Trp53 Trp53bp1`. The 8 that re-entered Tier A, with **no gene
removed**: `Atm Bax Bcl2l1 Cdkn2b Glb1 Mdm2 Sirt1 Trp53bp1`.

**Correction to my own phrasing, again.** "B7 had absorbed **27 arrest genes**"
overstates it: only **19 of the 27 are in `A0`**, the Tier A candidate pool. The
other eight — `Egr1 Fos Junb Nfkbia Sod2 Tgfb1 Thbs1 Timp1` — are SASP and
immediate-early genes that could never have entered Tier A. **Write: "B7 had
absorbed 27 genes, 19 of them from the Tier A candidate pool, of which 8 were
being held out of Tier A by B7 alone."**

**`Cdkn1a` is still not in Tier A — and the reason is stronger than stated.**
`Cdkn1a` *is* in `A0`, and it *was* one of the 27 B7 gave up. It stays out of
Tier A because it remains inside **three other** Tier B modules —
`B_tnfa_nfkb_proximal`, `B_interferon_response`, `B_downstream_arrest`. Tier A's
disjointness from `Cdkn1a` therefore survives C6 for three independent reasons,
not because B7 happens to still hold it.

**The lesson for the paper, restated:** not *"our sender set was circular"* but
***"a response module had absorbed the sender's genes, and the disjointness gate
resolved that by gutting the sender."*** Different diagnosis, different lesson,
and it is the one that generalises to anyone building disjoint marker tiers by
subtraction.

### 14.3 The second plank of the drafted restatement is refuted

I told the PI the restatement drafted in `CS_PHASE8_CALLERS.md` §3 could be used
"verbatim". Part of it cannot.

* **Plank 1 — "one pair is below chance" — replaced, from evidence already on
  disk.** The surviving pair is **SenePy vs DeepScence: pooled 0.737, z = −15.08,
  significantly below chance in 10 of 11 sections**. It is **byte-identical pre-
  and post-C6** because it does not involve Tier A, which is exactly what makes
  it robust. It strengthens to **0.495 (z = −28.9, 0/11 above chance)** under the
  caller-free re-anchor (`prolif` and `lmnb1` give identical decisions). **Caveat
  that must travel with it:** under |score| ranking the same pair is at chance
  (**1.025**, n.s.), so the anti-concordance is about **polarity**, not overlap.
* **Plank 2 — "the direction of each pair is predicted by its depth loading" —
  REFUTED by three tests.** (i) Exact permutation over the five non-circular
  pairs (2 same-camp vs 3 opposite, C(5,2) = 10 assignments): statistic +0.1895,
  **one-sided p = 0.300**, and the observed assignment ranks only 3rd of 10.
  (ii) Within-pair across 11 sections, the Q5/Q1 depth-loading product is
  *negatively* rank-correlated with agreement in **5 of 5** pairs: −0.164, −0.555,
  **−0.700**, −0.214, −0.527. (iii) Pooled continuous version over all 55
  pair-sections: **ρ = +0.096, p = 0.488**. It is broken by the largest value in
  the matrix — Tier A vs `Cdkn1a`⁺ at **1.471**, above chance 11/11, whose cause
  is **biological**: four of the eight genes that re-entered Tier A (`Atm`,
  `Mdm2`, `Trp53bp1`, `Bax`) sit on the p53 axis that induces `Cdkn1a`.
* **The weaker claim the data does support:** the dependence is weak,
  heterogeneous, and mechanistically *different at each end* (0.74× to 1.47×) —
  still incompatible with one latent state. The two-depth-camp fact survives at
  11/11 for all four callers: SenePy Q5/Q1 **10.58–41.74** and `Cdkn1a`⁺
  **4.19–42.36** (top-selecting) against Tier A **0.146–0.317** and DeepScence
  **0.218–0.795** (bottom-selecting).

**Provenance warning on those three tests.** **No script in the repo computes
them.** They exist as prose in `SUBMISSION_PATCH_2026-08-29.md:195-208` and
`COMPLETED_TASKS.md` (the caller-restatement section). They *are* reproducible — they were re-derived for
this ledger from `results/phase3/caller_within_type_depth_bias_11sections.csv`
(Q5/Q1 `enrichment`, pivoted section × caller) crossed with
`caller_agreement_matched_significance_11sections.csv` — but **do not cite a
producing script, because there is none.** If the pre-registration is to reference
these, a producer has to be written.


---

## 15. §4 (D-b)'s premise is refuted by measurement — `denoise=False` is now a *chosen* value

`reports/CS_PHASE8_D2_DENOISE.md`; `results/phase8_d2/`. This closes task 8.5 /
C7-D2, and it closes it the opposite way round from how the plan expected.

**§4 (D-b) states that DCA "is precisely the step that would normalize depth —
the confound under investigation."** It is not. It **raises** the depth loading,
on three of three sections.

`results/phase8_d2/d2_depth.csv`, rows `config = dca`. The metric
(`code/analyse_d2_denoise.py:51-55,159-167`) is the **Spearman ρ between the
DeepScence score and the per-cell `transcript_counts` column of
`data/raw/<section>/cells.parquet`** — panel-wide Xenium transcripts, not the
4,845 mapped genes used for normalisation — on the caller cell set with
`Low_quality`/`Unknown` dropped and the alternative score sign-aligned.

| section | `denoise=False` | `denoise=True` | Δ | ratio |
|---|---|---|---|---|
| 7239 | 0.3891 | **0.6404** | +0.2512 | ×1.65 |
| 7259 | 0.3176 | **0.5314** | +0.2138 | ×1.67 |
| 7352 | 0.4096 | **0.5419** | +0.1323 | **×1.32** |

**Correction to my own framing.** I have been saying denoising "roughly doubles"
the depth loading, and `PREREG_PHASE8.md:886` (P29) says so too. It does not:
the ratios are **×1.32 to ×1.67**. The D2 report's own wording — "raises it by
+0.13, +0.21, +0.25" — is the accurate one. **Write "raises the depth loading by
32–67 %, on three of three sections."**

### 15.1 It ran as published — path 1, not the fallback

`dca_venv_pip_freeze.txt`: **DCA 0.3.4, TensorFlow 2.4.4, Keras 2.4.3**,
numpy 1.19.5, anndata 0.7.8, scanpy 1.8.2 — 83 packages in an isolated
**CPython 3.8.19** interpreter (`dca_venv_python.txt`), driven out of process.
Per-run evidence in `runmeta_dca_7259/7352_*.json` and the three
`dca_sub20000*` metas, each carrying a `dca_bridge` block; 7239's bridge block is
null but `dca_worker_meta_7239_*.json` and `logs/deepscence_d2/dca_7239_*.log`
record the same stack at 83,392 cells / 8.58 min.

**The main environment still has no TensorFlow** — verified live:
`importlib.util.find_spec('tensorflow')` returns `None`, and `import keras`
(3.15.1, installed 03:09 UTC, before this work) raises
`ModuleNotFoundError: No module named 'tensorflow'`. The standing check is
recorded at `requirements.txt:41-46`.

All three DCA sections are **full sections**, not subsamples — `n_cells` 83,392
(7239) / 127,386 (7259) / 139,378 (7352). 7259 was OOM-killed at full section
twice against the **57.74 GiB** cgroup (`/sys/fs/cgroup/memory.max` =
61,999,996,928 B) before the surviving 111.8-minute run.

**Stated precisely, because "path 1 landed" is a slight overstatement.** The plan
(`Phase7_Minimal_Human_Replication (1).md:202`) specifies a separate **conda**
environment handing back an **`.h5ad`**. Conda/mamba/micromamba are absent from
this pod, so what was built is a python-build-standalone **venv** handing back
**`.npz`/`.npy`** (anndata on-disk format incompatibility). The report documents
both deviations in §3/§3.1. **It is path 1 in substance, with two documented
mechanism deviations** — not the fallback, and not literally what §4 wrote.

### 15.2 The published default is also seed-unstable, and nothing warns

`results/phase8_d2/d2_stability.csv`, the 20k-cell three-seed panel on 7239:

| pair | Pearson r | top-5 % Jaccard | cells changed | anchor sign flipped |
|---|---|---|---|---|
| dca seed0 vs seed1 | 0.5703 | **0.0000** | 2,000 | **False** |
| dca seed1 vs seed2 | 0.5732 | **0.0000** | 2,000 | **False** |
| dca seed0 vs seed2 | 0.9824 | 0.6653 | — | False |
| **raw seed0 vs seed1** | 0.9955 | **0.7606** | — | False |

Two seed pairs return top-5 % sender sets that are **perfectly disjoint**, against
0.761 for `denoise=False`. `anchor_sign_flipped = False` throughout, so this is
genuine instability in *which cells are called*, not a polarity flip.

**"No diagnostic fires" is literal.** `DeepScence/io.py:216-236`
(`fix_score_direction`) does `np.argmax(corr_metrics)` over the two bottleneck
nodes with **no tie test, no margin threshold, no floor on the metric and no
warning path**; `api.py:26` globally suppresses `UserWarning`. Seed 1's log
(`logs/deepscence_d2/stab_dca_seed1.log`) is clean — same "capturing scores in 2
neurons", normal training, normal write. The near-tie is visible only *post hoc*
in `runmeta_dca_sub20000_seed1_*.json`: `corr_metrics = [0.3732, 0.3354]`, against
seed 0 `[0.4636, 0.1791]` and seed 2 `[0.4634, 0.1517]`. Seed 1 also carries
`reverse: false` where seeds 0 and 2 carry `true`.

### 15.3 The uncomfortable half — and a caveat that weakens it

The denoised score is **better by DeepScence's own internal criterion** while
being *more* depth-confounded. The criterion is `direction_log.corr_metrics` for
the chosen node — the mean |Pearson| between the node's score and the
log-normalised, scaled expression of the CoreScence `occurrence ≥ 5` genes that
survive the low-variable filter (`io.py:148-183,209-212`; 31 genes).

7259: **0.125 → 0.469**. 7239: **0.158 → 0.497**. 7352: 0.127 → 0.673.

**Two corrections to how this has been reported.** The widely-quoted
"**0.126** → 0.470" baseline is not in any file: `runmeta_raw_7259_*.json` gives
**0.12548**. And the D2 report's §0 prose attributes 7239's baseline to the `raw`
control at 0.156, but raw 7239 is **0.15776**; 0.15644 is the **`mor`** run. Use
**0.125 → 0.469** and **0.158 → 0.497**.

**And the comparison is not apples-to-apples.** `DeepScence/api.py:68-71` runs
`dca()` *before* `read_dataset()`, so under `denoise=True` `corr_metrics` is
computed against the **DCA-imputed** matrix while under `denoise=False` it is
computed against real counts. Imputation mechanically inflates gene–score
correlations. **The "better by its own criterion" half is therefore partly an
artefact of the metric being measured on the smoothed matrix**, and the D2 report
does not say so. The finding survives — the tool improves on its stated objective
while getting worse on the axis the science depends on — but it is weaker than it
has been presented.

### 15.4 The rest of D2, and one rounding to fix

* The top-5 % call set becomes **99.97–100.0 % hepatocyte** (`d2_celltype_composition.csv`,
  `call_pct_alt`: 99.973 / 99.965 / 100.0, from a committed 71.5 / 64.5 / 97.2).
  "100 % on all three" is a rounding — 7239 keeps DC at 0.027 %, 7259 keeps
  biliary and proliferating at 0.017 % each.
* Sender-set agreement between configurations, `global_top5_jaccard`:
  **0.118 / 0.141 / 0.280**.
* **§6's own named estimator was too weak to test its own question.** Median-of-ratios
  (`poscounts`) removes **11.4–23.6 %** of log-depth variance
  (`d2_normalisation_strength.csv`); a `lib` configuration removes **100 %** and
  cuts the depth loading by **74.2 %** (0.318 → 0.082) and **92.9 %**
  (0.410 → 0.029).
* **The frozen base is re-derivable, not merely intact.** `sha256sum -c
  results/phase8_d2/committed_deepscence_sha256.txt` returns **OK on all 11**
  committed `deepscence_*.csv`, and same-seed re-runs reproduce them at
  r = **0.99999913** (7239) and **0.99999995** (7259) — including
  `deepscence_sbr.csv`, written 2026-08-20 17:57, **before two container
  rebuilds**. With 202 untracked files holding the evidence base (§18), a
  demonstrated re-derivation carries real weight.
* **Resource fact, corrected.** DeepScence holds five full-size dense copies at
  once (`api.py:60-71`, `io.py:44`, `io.py:64-65`). At 83,392 × 4,845 float32
  that is 1.62 GB each, **~8.1 GB**, not 16. The **~16 GB per 83k-cell section**
  figure is **empirical** — `code/run_d2_stream.sh:9-13` records an OOM at 11 GB
  free — and must not be presented as arithmetic from the five arrays.

**Decision taken under delegated authority: freeze `denoise=False` as PRIMARY, now
a *chosen* value, with `denoise=True` as the published-default sensitivity.** This
converts a caveat ("we deviated because DCA would not install") into a finding
("we ran the published default and it is more depth-loaded and seed-unstable, so
we chose against it"). `PREREG_PHASE8.md` §3.9 now opens *"RESOLVED. DCA
installed, ran, and lost on the merits"*, and the old sentence — *"DeepScence as
we can run it on this panel, not DeepScence as published"* — is **withdrawn**.
Three prohibitions follow (§10.9–§10.11): no `denoise=True` number from a single
seed without its seed-stability companion; `mor` may never be cited as evidence
that normalisation cannot move this caller; `rho_signed_dz_vs_depth` on `raw`
rows must not be quoted.

**Also: §8's confound is halved.** DCA now runs on **both** arms, so
`denoise=False` can no longer explain a mouse-only artefact. Two falsifiable H1
predictions were registered on the back of it (`PREREG_PHASE8.md` §8 P-vi,
P-vii), each stating its own falsifier.


---

## 16. Moran's I — an argument the project was about to make, falsified before it was made

`reports/CS_PHASE8_MORAN.md`; `results/moran/`. Run on all 11 M1 sections, all
**13,590** `cell_feature_matrix.h5` features individually, 16 aggregate fields,
**7 weights variants**, 2 centrings, **999 permutations**; `data/raw_h1/` untouched
and the figure written under `results/`, not `figures/`, so the guard is
undisturbed (re-run at the end: OK, 52/52). The producer was validated against
`esda.Moran` to **d = 0.00 × 10⁰** on I and `z_rand`; I re-ran
`code/run_moran_controls.py --validate` independently and it passes.

**Why it was run at all:** master plan §29 objection 9 promised *"we report our
own Moran's I on the controls alongside the kernel amplitude"* and
`grep -ril moran` over `code/` and `results/` returned **nothing**. The objection
was answered by assertion, and it is the designated defence against the one piece
of prior art that directly falsifies a project novelty claim.

### 16.1 The falsified recommendation

`reports/NOVELTY_ASSESSMENT.md` §2.1 point 3 ends: *"…with a Moran's I of your own
controls alongside, **so the reader can see the two tests disagree.**"*

**They do not disagree.** Across the twelve control + module fields, |Moran's I|
and |A7 naive amplitude| rank together at **Spearman ρ = +0.895 (p = 8.4e-5)**
raw and **+0.944 (p = 3.9e-6)** cell-type-centred
(`results/moran/moran_verdict.txt`; I recomputed both from `moran_pooled.csv`
and got 0.8951 / 0.9441). Written as advised, **a reviewer with `results/moran/`
in hand would have found the opposite in a single plot.** The report now carries
a correction banner at the top saying so.

**Correction to how this has been relayed, including by me.** The banner and my
own summaries say `NOVELTY_ASSESSMENT.md` §4 **O1** was falsified too. **It was
not.** O1's "evidence needed" cell asks the project to *"report your own Moran's I
on the controls alongside the kernel amplitude, and state that a near-zero global
autocorrelation does not preclude a projection onto a specific covariate"* — a
statement the Moran run **endorses** (§4.2 conclusion 3). What actually happened
to O1 is narrower and should be written that way: its *premise* is gone (the
aggregate control field is **not** near-zero, I = +0.0455) and its status moves
from NOT DONE to done. Only §2.1 point 3 was falsified.

### 16.2 The replacement, which is stronger because it is measured

`results/moran/moran_kernel_power.csv`, 22 section × sender-call cells
(11 sections × `tierA_p95`, `cdkn1a_pos`). Medians recomputed by me:

| quantity | median | range |
|---|---|---|
| λ̂ (naive, `all_controls`) | 35.6 µm | 7.0 – 50.0 |
| Moran's I of the kernel covariate itself | +0.797 | +0.412 – +0.852 |
| **ΔI contributed by the entire A7 gradient** | **2.20 × 10⁻⁴** | 1.4e-6 – 2.6e-3 |
| **ΔI as a fraction of the observed control I** | **0.83 %** | 0.002 % – 6.1 % |
| SE(I) on the same cells | 1.50 × 10⁻³ | 1.2 – 2.3e-3 |
| **β̂ that Moran's I could just detect (ΔI = 2 SE)** | **0.362 SD** | 0.308 – 1.070 |

**0.362 SD is larger than anything this paper measures.** The A7 control gradient
is 0.074 SD (4.9× smaller). The naive biological amplitude on the frozen tree is
**0.277 SD** (`results/phase3/a7_summary.csv`, `BIOLOGICAL MODULES (reference)` /
`base`) and the Phase 3 headline naive amplitude is **0.329**
(`code/m1_headlines.py`) — **Moran's I could not have detected the paper's own
headline effect either.** So "different question" survives, justified by **power,
not orthogonality**. That is a quantified argument replacing an assertion that was
about to be wrong.

**Number correction.** The Moran report, the novelty banner,
`CS_PHASE8_CALLERS.md` §4.1 and `PREREG_PHASE8.md:751` all quote the biological
amplitude as **0.291 SD naive / 0.036 SD conditioned**. Those are
`results/phase3_pre_c6/a7_summary.csv` — the **superseded** pre-C6 snapshot this
ledger's §6 documents. The frozen values are **0.277** and **0.031**. The banner
additionally attributes 0.291 to `moran_kernel_power.csv`, where it does not
appear. The conclusion is unaffected (0.362 > 0.277 as well as > 0.291); the
digits are not.

### 16.3 A second, genuinely new point against the prior art

Voyager's Xenium vignette reads a near-zero per-feature Moran's I on negative
controls as "no technical artifact spatial trend". On this data that reproduces
exactly — median per-feature I of **−0.00012** for the 40 negative control probes,
−0.00004 for the 609 codewords, −0.00029 for the 21 genomic controls, against
**+0.00999** for the 5,106 real genes.

**But gene-expression features matched to the controls on total counts give
indistinguishable values** — **−0.00018** for the probe-matched set, −0.00012 for
codeword-matched, −0.00025 for genomic-matched; if anything slightly *more*
negative than the controls themselves. A negative control probe carries a median
of **21** counts per section; the median gene carries **5,884.5** — a **280×** gap.
**Voyager's controls-vs-genes contrast is an abundance contrast**, and the
per-feature statistic has no power at control abundance. It holds under
CP10K + log1p as well (probes −0.000108, codewords −0.000042, genomic −0.000282,
genes +0.007176 over 3 sections).

*Caveat on that last sentence:* `moran_per_feature_lognorm.csv.gz` carries no
`total_counts` column, so **the count-matched contrast itself was not recomputed
under CP10K + log1p** in any shipped artefact. The per-class medians are what
reproduce. The matched contrast does hold when the counts are joined back in from
`moran_per_feature.csv.gz`, but that derivation is not in `results/` and must be
regenerated before it is cited there.

**Aggregate the same features the way A7 does and the controls are not flat.** The
per-cell sum over all 670 control features has Moran's I = **+0.0455
[+0.0302, +0.0609]** (section-clustered, k = 6 NN), above **86.8 %** of the 56,166
gene × section values — verified by recomputation from `moran_per_feature.csv.gz`.

**A7's three-way split reproduces from an independent statistic** with no kernel,
no λ, no sender call and no nuisance design: `all_controls` **+0.0455**,
codewords **+0.0421**, probes **+0.0058** — the codeword/probe ratio is **7.2×**,
against A7's own 2.7×. Same ordering, different instrument.

### 16.4 Reported against interest

* **Moran's I does not call the probes flat.** `moran_pooled.csv`
  `neg_control_probe` `I_raw_lo` = **+0.0033** — the pooled CI excludes zero.
* **It ranks `genomic_control` lowest** (I = +0.0042, below even `neg_probe_rate`
  at +0.0047) where A7 ranks it third of five by |amplitude|.
* **Within the control family the two statistics are uncorrelated** — ρ = +0.155,
  p = 0.26 over 55 section × response pairs — which is what the power calculation
  predicts, and is the honest form of "these are different questions".
* **The 14 discordant cells must not be quoted as the disagreement the plan
  wanted.** A7's own CI-exclusion rate on a response with no biology is inflated,
  so a 25–50 % exclusion rate in one section is partly the estimator's type-I
  error. *Mis-citation to fix:* the report compares those cells against
  "9–16 % against a 5 % nominal", but the discordance table's
  `a7_frac_CI_excl_0` is the **naive** design (`code/summarize_moran.py:73`),
  whose control CI-exclusion rate in `results/phase3/a7_verdict.txt` is
  **19–43 %**. The correct comparator makes the report's own argument stronger.
* **"No ordering changes across the 7 weights variants" is false as written.**
  14 of 120 field pairs flip; `genomic_control` vs `neg_probe_rate` swap between
  knn10 and knn20 and stay swapped for all three distance bands, and at band100
  `all_controls` (0.0353) *exceeds* the Tier B module `downstream_arrest`
  (0.0289), so §0's "between a fifth and a half of the Tier B modules" is
  knn6-specific. The orderings the argument rests on — controls ≪ modules,
  probes ≪ codewords, `all_controls` above the median gene — hold at every
  variant. **Write "no ordering that carries the argument changes."**


---

## 17. Corrections to the record itself

A corrections ledger has to include the ones that are embarrassing. These are
corrections to the project's own reporting layer — the bibliography, the plan
documents, the audit trail and my own statements — rather than to a measurement.
The pattern across all of them is the same: **the science reproduces; the
reporting layer did not.**

### 17.1 The bibliography was written from recall

`reports/CITATION_AUDIT.md`; `references.bib`.

**19 of 32 entries carried invented author forenames — 41 wrong given names.**
The pattern (plausible, correctly transliterated, wrong) is diagnostic of a
bibliography written from memory rather than from retrieved records. An earlier
audit (gene-set deviation series D7) had checked journal, volume, pages and DOI
and passed it — **it never checked forenames**.

The one to put in front of the PI: **`Martin, Luke` was actually `Martin, Lucy`**
(`references.bib:251`, now `Martin, Lucy and Schumacher, Linus and Chandra,
Tamir`; the pre-audit form is in `git show a6aac3a:references.bib`). A real
researcher was misgendered in a bibliography that was about to be submitted.

**Three counts in circulation are wrong and should not be carried forward.**

| as circulated | corrected | why |
|---|---|---|
| "30 entries" (`references.bib:23`, `COMPLETED_TASKS.md` novelty table, `NOVELTY_ASSESSMENT.md:205`) | **32 entries before, 43 now** | 30 was D7's count of *references in master-plan §31*, not bib entries. `grep -c '^@'` gives 43 now and 32 at `a6aac3a` |
| "49 author lines corrected" (`COMPLETED_TASKS.md` row 40) | **18 author lines, 37 given names** | there are only 43 author lines in the file; a keyed diff of old vs current gives 18 changed |
| "verified against Crossref AND PubMed independently (they agree)" | **"against the Crossref deposit, and against PubMed where a record exists"** | 6 of the 19 corrected entries are preprints/manuals with no PubMed record — 12 of the 41 names. `gu2026identifiability` rests on the arXiv abstract page alone |

**⚠ And one correction was never actually applied.**
`references.bib:445` — **`kumar2026cellwhisper` still carries its four wrong
forenames** (`Abhishek / Fernando / Bhavya / Nan`), byte-identical to the
pre-audit version, **while lines 482–485 of the same entry assert they were
"corrected above"**. So the true tally is 37 of 41 given names fixed in 18 of 19
entries; four remain wrong, **in the single most load-bearing citation in the
project**. `SUBMISSION_PATCH_2026-08-29.md:508-514` tells the author the
bibliography is now clean and names this entry among the fixed. **It is not.
Fix the entry or withdraw that line before the patch is hand-applied.**

Also: 11 spatial-statistics references were added, closing a structural gap — the
project is an instance of the *spatial confounding* literature and cited none of
it. All 11 carry a `% SUPPORTS:` line so nothing is decorative; 10 of 11 also
carry `% VERIFICATION:` (`spatstat2026rshift` puts its provenance in `note`
instead). Four sources the plan §30 and `SUBMISSION_PATCH` §7.4 both say the
paper needs — ICE (*Genome Biology* 2026), markeR (*NAR GAB* 2026), bioRxiv
2026.01.02.697374, stAge (bioRxiv 2025.11.23.689860) — are **still absent**, and
were deliberately not added from recall.

### 17.2 CellWHISPER was miscited as the source of the torus null

**CellWHISPER's null is a within-cell-type location permutation — this project's
N1.** The strings *torus*, *toroidal* and *wraparound* appear nowhere in that
paper. N3 is what Figure 4 rests on, so the citation was attached to the wrong
null entirely.

**Corrected in 9 places**, each carrying a dated `citation audit CIT-1` note:
master plan §4.3 (`:194`), §22 Step 3 (`:1025`), §23 (`:1056`), §31 ref 27
(`:1264`), §32 item 3 (`:1286`); `reports/CS_PHASE1.md:275`,
`CS_PHASE2.md:348-349`, `CS_PHASE3.md:329`; and the `AUDIT (CIT-1)` block at
`references.bib:451-467`. A repo-wide sweep finds no residual file making the
attribution. *(`CITATION_AUDIT.md` §0.1's prose says "six places" while its own
table lists nine — six carried the torus-attribution sentence specifically, three
carried the adjacent "coordinate-randomized" error. **Use 9.**)*

**Correct attribution: Lotwick & Silverman (1982)**, *JRSS-B* 44(3):406–413,
doi:10.1111/j.2517-6161.1982.tb01221.x (`references.bib:559`).
⚠ **Verified indirectly only** — the 1982 full text was never rendered
(JSTOR/Wiley paywall); the attribution rests on Mrkvička et al. (2021) §1, which
quotes p.410. The bib entry says so. Do not present it as read from primary.

A related drift was corrected in the same pass: CellWHISPER's ">90 % FPR" figure
is an **interaction-count ratio from which an FPR is inferred**. Plan §1, §4.3 and
§31 correctly said "implying"; **§5 had drifted to asserting it**, and is now
`"implying false-positive rates above 90 %"` (`:202-204`, CIT-5), with §2.2
(`:144`) softened to match.

### 17.3 H1 is seven human spleens, not "human aging lung"

`Phase7_Minimal_Human_Replication (1).md` — the plan of record for the human arm —
described H1 as **human aging lung**. H1 is **GSE326743, 7 normal human spleens**
(ages 17–59, 4 M / 3 F, Xenium Prime 5K + 100 add-on, panel verified **on data**
at 5,093 genes, 2,207,593 cells; `reports/PHASE7_H1_SCREEN.md`). No file in the
repo disagrees on the accession or the donor count.

**Correction to the count I have been giving.** "5 places" is wrong under every
reading. The literal string "human aging lung" occurs on **2** body lines
(`:282`, `:482`); "aging lung" on 4; and **15 distinct body lines across 7
sections** (§12.1, §13/A6, §14, §15, §17, §20, §22) describe H1 as lung, airway
or alveolar. **Write: "H1 was described as lung on 15 lines across 7 sections."**

**And the fix was a banner, not an edit.** All 15 lines are still present
verbatim; a superseded banner was prepended at `:3-26`. That was the right call —
the document is retained as the plan of record — but it must be stated that way,
because a reader following §30 5.8's pointer into §17 still lands on the lung
text *and* on pre-C6 numbers (0.326 / 0.027 / 0.203 / SF-N2 0.943 /
"66–76 %" / "0.93–1.22×"), every one of which this ledger supersedes. Two further
accuracy notes: the banner enumerates **6 bullets plus a lead claim**, not seven
items; and its own scope line names only §12.1, §13 and §17, leaving §14's lung
marker panel (AT1/AT2/`SFTPC`/alveolar macrophages, `:420-422`), §15, §20 and §22
unflagged.

**Downstream consequences already absorbed:** the A6 covariate moved from the
lung airway-to-alveolar axis to spleen **red pulp / white pulp**
(`genesets/human/D_spleen_*.txt`, five sets: t-zone 24, follicle 13, marginal
zone 8, red pulp 28, capsule/trabecula 18; gene-set deviation D6 —
*"H1 is spleen. Every word of the lung version is void for this arm"*). And
**SenePy has no spleen signature at all**: 65 human hubs across 10 tissues, liver
present, spleen absent; of 22 spleen cell types, **0 matched / 15 cross-tissue
surrogate / 7 unscoreable**.

### 17.4 `NOVELTY_ASSESSMENT.md` advised an argument that Moran's I then falsified

Covered in §16.1. Recorded here because it is a correction to *advice the project
was about to act on*, caught by running the measurement instead of taking the
advice. The report now carries a correction banner (`:3-10`).

**⚠ The banner ends "Everything else in this report stands." It does not.**
Two things in the same report are still wrong: §2.1's framing box (`:250`, `:261`)
and objection O2 (`:507`) still attach the **−0.070 SD** gradient to
**"negative control probes"** — the exact response-naming error of §17.5 — and
§U2 (`:287`) misquotes Mrkvička as *"can also be used for irregular windows"*
where the paper says *"can be applied in case of general (compact) observation
windows"* (flagged at `CITATION_AUDIT.md:288-292`, still live).

For the record, the three novelty claims the report overturned, all of which
stand: the negative-control-probe spatial diagnostic is **not** new (Voyager's
Xenium vignette; Ren et al., *Nat Commun* 16, 2025) — what is new is fitting *the
estimand's own estimator* to the controls, a Lipsitch (2010) negative-control
**outcome**, plus the N2-vs-N5 result; the torus degeneracy is **1982-vintage
statistics**; and "nobody reports a measured FPR" is **no longer safe**
(CellWHISPER Jan 2026, CONCISE Jun 2026). The central negative result is
**intact** — no published length constant for senescence spatial influence exists
to contradict. On venue, the plan's primary/secondary is **backwards** for a
negative result: ICBINB-BIO matches topic-by-topic and gives 8 pages;
ml4spatialbio is 4 pages, non-archival, permits concurrent submission. *(The
report could not verify ICBINB-BIO's dual-submission policy.)*

### 17.5 The A7 response-naming error — and this file was one of the offenders

`SUBMISSION_PATCH_2026-08-29.md` stated a **"−0.070 SD gradient in
negative-control probes"**. That is the pooled `all_controls` response — 40
probes + 609 codewords + 21 genomic controls, of which the codewords carry ~73 %
of the counts. **The 40 negative-control probes are the response
`PREREG_PHASE8_genesets.md` §11 designates the PRIMARY technical null for test
A7, and on their own they are flat.** So on its own pre-registered primary
response, **M1's A7 passes naively**. The "assay is not flat" verdict survives on
the codewords and the genomic controls; the *name* on it did not.

**⚠ `reports/CORRECTIONS.md` §12 — this file — carried the same error**, in its
heading ("A7 — the negative-control-probe kernel") and in its text ("the
negative-control-probe kernel is −0.074 SD at p = 0.015"). `SUBMISSION_PATCH` §4a
calls it out by name. **Both are corrected in this pass.** So is the same
phrasing's presence in `CS_PHASE8_CALLERS.md` (`:271, :283, :290, :348, :435`),
which is **still live and still needs fixing.**

### 17.6 The A7 per-response numbers in circulation are the pre-C6 vintage

**Found while writing this ledger, and it is the widest-spread staleness in the
repo.** A7 was re-run at 09:06 under the frozen C6 sender calls
(`results/phase3/a7_summary.csv`). Five documents still quote the **05:19 pre-C6**
file (`results/phase3_pre_c6/a7_summary.csv`), in most cases while citing the
post-C6 path.

| response, `design = base`, `clustered_mean` | pre-C6 (as quoted) | **post-C6 (frozen)** |
|---|---|---|
| `all_controls` | −0.0697 [−0.1276, −0.0118], p = 0.023 | **−0.0744 [−0.1306, −0.0182], p = 0.0145** |
| `neg_control_probe` (the 40) | −0.0177 [−0.0453, +0.0099], p = 0.183 | **−0.0225 [−0.0527, +0.0078], p = 0.129** |
| `neg_control_codeword` (609) | −0.0549 [−0.1064, −0.0034], p = 0.0389 | **−0.0604 [−0.1085, −0.0123], p = 0.0188** |
| `genomic_control` (21) | −0.0337 [−0.0538, −0.0136], p = 0.00389 | **−0.0307 [−0.0558, −0.0056], p = 0.0213** |
| `BIOLOGICAL MODULES`, naive / conditioned | 0.2914 / 0.0356 | **0.2767 / 0.0310** |

Quoted stale at `PREREG_PHASE8.md:756-758` and `:751`,
`AUDIT_CORRECTIONS_APPLIED.md:202-203`, `AUDIT_PHASE8_FACTCHECK.md:52-54`,
`COMPLETED_TASKS.md` ("two things needing the PI" item 2),
`CS_PHASE8_CALLERS.md` §4.1, `CS_PHASE8_MORAN.md`
§4.2 and §5, and the `NOVELTY_ASSESSMENT.md` banner. `PHASE8_ROADMAP_STATUS.md`
and `PREREG_PHASE8.md` P2 at least label the magnitudes **PROVISIONAL**;
`SUBMISSION_PATCH` §4a has already been re-derived on the post-C6 file.

**Every conclusion survives the update, in direction and in sign.** The probes are
flat under both vintages (p = 0.183 → 0.129); the codewords and genomic controls
are significant under both; `all_controls` is not flat under either. **Only the
digits are stale — but a ledger that prints pre-C6 digits beside a post-C6 file
path is a fresh, checkable error, so pick one vintage and cite the file that
holds it.** This ledger uses the post-C6 column throughout.

### 17.7 The CoreScence circularity anchor was a typed-in literal

I reported "CoreScence circularity rose **69 % → 88 %**, a real cost of C6"
several times, and used it to support PI decision D6. **The 69 % was fabricated.**

* `24/35 = 69 %` was a **hand-typed literal in two scripts** —
  `code/make_figure_genesets.py:245` (`bars = [("M1 mouse\npublished", 24, 35, …)]`)
  and `code/gate_disjointness_human.py:291`
  (`mouse_arm_reference='24/35 = 69%'`, a *string*). The numerator 24 is right;
  the denominator is not. The project's own committed
  `results/phase3/n8_disjointness_*.csv` (all 9, git-tracked) carry
  **`corescence_on_panel = 33`**, and `logs/ds_smoke.log:2` gives **31** under the
  strict ortholog map.
* Re-derived from files by `code/corescence_circularity.py::derive()`, matching
  the committed `results/phase7_jobA/corescence_circularity_mouse.json`
  byte for byte:

| configuration | mapping | on panel | in ≥1 Tier B | % |
|---|---|---|---|---|
| pre-C6 | strict MGI | 31 | 24 | 77.4 |
| **pre-C6** | fallback | 33 | **26** | **79 %** |
| C6-promoted | strict MGI | 31 | 28 | 90.3 |
| **C6-promoted** | fallback | 33 | **29** | **88 %** |

* **Two errors compounded.** A fabricated baseline; and a "69 → 76 → 88" story
  that compared a *mouse* baseline against two *human* configurations. Measured
  **within each arm**, the cost of C6 is **+9 points on mouse (79 → 88)** and
  **+12 on human (76 → 88)** — roughly **half** what I attributed to it. Both arms
  land on **29/33 = 0.8788**, identical to four decimal places.
* **The direction cuts both ways and both halves matter.** The *published mouse
  arm is more circular than reported* — 79 %, not 69 % — which is against
  interest; while the *cost of C6* is smaller than claimed. It also flips
  `BIO_PHASE7_JobA.md` §4's conclusion from "worse in human than mouse" to
  "comparable on the two arms". **PI decision D6 (strip-and-refit primary) is
  unaffected**: 88 % is high enough to justify it either way.
  `figure_gs3` now shows four bars so the within-arm cost is legible, and the gate
  and figure **compute** the number instead of asserting it.
* **Tier A ∩ CoreScence went 2/25 → 4/33** — `{Foxm1, Parp1}` gaining `Cdkn2b`
  and `Mdm2`. Both are also among the 8 Tier A re-entrants of §14.2, so the two
  stories join here: making the sender set less hollow also made it more
  CoreScence-overlapping, and that touches the arrest-vs-DeepScence caller pair.

**Two provenance claims I have made about this must be softened.**

1. **"Verified from git history" — it cannot be.** Both scripts carrying the
   literal are **untracked** and appear in no commit, so there is no recoverable
   "before" state. The evidence is contemporaneous and adequate — the fact-check
   audit quoted both lines with line numbers (`AUDIT_PHASE8_FACTCHECK.md:243-256`)
   and both files now carry self-reporting comments — but write *"per the M1 audit
   item and the scripts' own annotations"*, not *"from git history"*.
2. **"A denominator that exists in no file" — too strong.** 35 does exist in the
   project, as the size of the *superseded human* B7:
   `figures/figure_gs3_corescence_circularity_data.csv` row 3
   (`H1 human superseded B7 (n=35)`) and `genesets/human/README.md:119`. That is
   the most likely provenance of the mix-up. The defensible claim, which does
   hold, is narrower: **35 is reproducible under no CoreScence-on-panel mapping
   convention** — every committed denominator is 31 (strict) or 33 (fallback).

### 17.8 Two column conflations that are still live in the reports

**(a) "23 % of shifted senders in the void" is the wrong column** (fact-check item
R3, **not yet applied**). 22.8 % and 8.0 % are
`1 − frac_retaining_a_neighbour` — shifted senders that lost **every** real
neighbour inside the 100 µm window. The **out-of-tissue** fraction is
`1 − frac_in_occupancy`. Medians over the six in-band sections of
`results/phase3/null_destructiveness.csv`, recomputed here:

| null | out of tissue | lost every neighbour |
|---|---|---|
| N3 (bounding-box torus) | **35.5 %** | 22.8 % |
| N4 (rotation) | **19.9 %** | 8.0 % |

Both quantities are real and both matter; they must never be merged into one
sentence. Still live in `CS_PHASE7_C1.md` §0/§6.2, `CS_PHASE8_C1_CLOSEOUT.md`
§4.1 and its §4.2 caption, `PHASE8_ROADMAP_STATUS.md`, **and in
`CS_PHASE8_TORUS_VAR.md` §1 and §10** — which conflates them despite
`CS_PHASE7_C1.md:95-103` carrying an explicit box warning not to. The §10
submission-framing sentence should read **35.5 %**.

**(b) "Every control family's clustered mean is indistinguishable from zero" is
false** (R5, **not yet applied**). It is 4 of 5: `neg_probe_rate` under `n6n5`
was +0.0108 [+0.0021, +0.0195], p = 0.0204 pre-C6. *(Post-C6 it is
+0.0097 [−0.0060, +0.0253], p = 0.199 — the caveat has since disappeared, which
§12 records.)*

### 17.9 Corrections I made to my own statements during the window

Collected so the PI can see what was wrong while decisions were being taken on it.

| what I said | what is true |
|---|---|
| "the independence claim rested on a **circular sender definition**" | the defect was in the **response** module; the gate gutted the sender by subtraction (§14.2) |
| "B7 had absorbed **27 arrest genes**" | 27 genes, **19** of them from the Tier A candidate pool (§14.2) |
| "CoreScence circularity rose **69 % → 88 %**" | **79 % → 88 %** on mouse, **76 % → 88 %** on human (§17.7) |
| "denoising **roughly doubles** the depth loading" | raises it by **32–67 %**, on 3 of 3 sections (§15) |
| "all three `fig_phase3_*` figures are on a stale basis" | **one** of three. Only `fig_phase3_caller_depth` reads a caller table; the agent checked *producers*, not filenames |
| "figures 2a and 2d are exempt — not null-dependent, verified byte-identical" | true of the **null** correction, false of a **gene-set** change. Both are functions of the sender set, which went 25 → 33. Both regenerated, both moved |
| "the composition surrogate share 66–76 % is **untraceable**" | too strong. Both numbers have producers; what has no producer is the **pairing** of them into one range (§7, §18) |
| "H1 was called lung in **5 places**" | **15 lines across 7 sections** (§17.3) |
| "the memory ceiling is ~251 GB" | the cgroup limit is **57.7 GiB**; `free -g` reports HOST memory. This is what caused the DeepScence OOM kills — five sections were run concurrently against a quarter of the believed budget |
| "0 files written in the last hour" (repeatedly) | a **false negative** every time. `find` here is **`bfs`**, which rejects `-newermt "-60 minutes"` and errored out silently; the real count was 138 files in one hour. Caused one false stall alarm |
| "the fact-check audit's R6 numbers stand" | they were correct for the **pre-C6** A7 file of 05:19 and stale after 06:52 (§18) |

**Two traps caught by agents that are worth keeping as method notes.**
`figure2a` **caches its own input** (`figure2a_stratified_curves.csv`, rebuilt only
if absent) — left in place it would have regenerated from **pre-C6 senders** and
come back byte-identical, **reading as a passing reproducibility check when it was
a stale-input result**. And `code/tierC_lr.py:48` *reads* `sender_flag_p95` — a
column that did change — but never uses it; only re-running the producer and
diffing showed the read was inert. **Reasoning from "it reads a changed input"
gives the wrong answer in both directions.**


---

## 18. What is still wrong, unverified, or unfinished

Measured at **2026-08-27 10:50 UTC**, not transcribed. Ordered by what it costs
if it is missed.

### 18.1 P0 — the evidence base is not committed, and that is the only thing the PI can fix

**The last commit is 2026-08-21.** Every artefact produced today is uncommitted.

| directory | files on disk | tracked | untracked (not ignored) | modified |
|---|---|---|---|---|
| `results/` | 1,295 | **266** | 255 | 67 |
| `genesets/` | 172 | 68 | **104** | 3 |
| `figures/` | 63 | 27 | 25 | 16 |
| `data/processed/` | 137 | 40 | 0 | 0 |
| `code/` | 154 | 82 | **72** | 11 |
| `reports/` | 34 | 11 | **23** | 5 |
| **whole repo** | 35,445 | 503 | **482** | **107** |

**⚠ `git status` will not show you this.** It collapses wholly-untracked
directories into one entry, so `git status --porcelain | grep -c '^??'` returns
**202**, not 482. The file-level count is
`git ls-files --others --exclude-standard | wc -l` = **482**. Use the second one
when deciding what to `git add`.

`genesets/human/` (49 files) and `genesets/mouse_c6/` (19 files) are **entirely
untracked**. The figure "194 untracked" that has been circulating is **stale by
2.5×** — it is now 482 repo-wide; its companion "107 modified" is still exactly
right, which is why the pair reads as current.

**Three things not previously recorded:**

1. **23 of 34 files in `reports/` are untracked**, including `PREREG_PHASE8.md`
   and **this file**. The documentary record is as unpinned as the evidence. So
   is half the code: **72 of 154 files under `code/` are untracked**, including
   both scripts that carried the fabricated CoreScence literal (§17.7) — which is
   why that correction has no recoverable "before" state.
2. **`results/phase3_pre_c6/` — 98 files, 0 tracked** — is the *sole copy* of the
   pre-C6 baseline that licenses every pre/post comparison in this ledger. The
   `pre-c6-genesets` tag does not contain it, because those files were never
   committed. It exists only on the network volume.
3. **A `git checkout` is currently destructive.** `code/build_genesets.py` **at
   `HEAD`** has `SCRATCH` pointing at a dead per-session `/tmp` path from a
   different session; running HEAD's copy globs zero JSONs and overwrites
   `genesets/*.txt` with **empty Tier B modules**. The `SystemExit`-on-zero-glob
   guard exists **only in the uncommitted working tree**. So the untracked state
   is not merely a provenance gap — it is the only thing between a routine
   checkout and destruction of the frozen gene sets.

A container wipe would not destroy any of this (`/workspace` is the network
volume; `/` is the ephemeral overlay, which is what was lost **twice today**), but
nothing is version-pinned, so **`phase8-frozen` would reference nothing
immutable** and the pre-registration cannot honestly cite immutable artefacts.
This is fact-check item **M8**, the only audit item no one has applied, and it
blocks 8.9, all of Phase 9 and 10, and the A6 build. I have committed nothing,
per standing instruction.

`PREREG_PHASE8.md` has **7 occurrences of `TBD`, of which 5 are live** (lines 28,
58, 66, 67, 68) — **4 tag commit hashes and one `Tag date (UTC)`**, not five
hashes. `git tag -l` shows only `pre-c6-genesets`; `phase8-frozen` does not exist.

### 18.2 P1 — manuscript-bound, and the deadline is in two days

* **⚠ `SUBMISSION_PATCH_2026-08-29.md` §9 still instructs the falsified framing.**
  It tells the PI to write *"so the reader can see the two tests disagree"* —
  the sentence `CS_PHASE8_MORAN.md` §4.1 measured to be false (ρ = +0.895 /
  +0.944) and says *must not be written*. **This is the most dangerous unpatched
  line in the repo, because it sits in the one document the PI applies by hand.**
  The same falsified text is unedited in `NOVELTY_ASSESSMENT.md` §2.1 point 3 and
  §4 O1 — corrected by a top banner only, body untouched, and O1 still reads
  "NOT DONE" for a run that has completed.
* **`references.bib:445` — `kumar2026cellwhisper` still carries four invented
  forenames** while its own audit note says they were corrected (§17.1), and
  `SUBMISSION_PATCH` §508-514 tells the author the bibliography is clean.
* **The patch cannot be verified from here.** The manuscript is not in this repo;
  15 search-and-replace strings plus five sections must be applied by hand and
  nothing can check the result. It has already gone stale once **within 15
  minutes** and been corrected three times.
* **The venue decision is open and the two options have different deadlines.**
  ICBINB-BIO Aug 29 11:59 AoE (8 pp) against ml4spatialbio Aug 29 – **Sept 4**
  AoE (4 pp, non-archival, concurrent submission permitted). Plan §29 still
  records a hard "Deadline August 29 AoE" for ml4spatialbio. ICBINB-BIO's
  dual-submission policy **was never verified**.
* **Four reviewer objections are genuinely NOT DONE and are cheap:** O4
  (CellWHISPER mis-cited as the torus-null source — Figure 4 rests on it), O7
  (the DeepScence sign flip is the documented `CDKN1A` anchoring rule), O11
  (caller disagreement framed as a finding), O12 (spatial-statistics citations).
  O1, O3 and O6 are *labelled* NOT DONE but their runs have landed — those labels
  are stale, not the work.
* **Two live citation liabilities.** `neretti2024dissecting` is a **single-author
  conference abstract** and is the sole support for the primary/secondary
  distinction that underwrites Tier B module B7 — flagged since D7 and unmoved.
  `acosta2026correction`'s content is **unread** (Nature paywall) and needs a
  human with library access. ICE, markeR, stAge and the Jan-2026
  circular-validation preprint are still absent from `references.bib`.

### 18.3 Numbers that are still stale in live documents

| where | says | should say |
|---|---|---|
| `CS_PHASE8_CALLERS.md:116,127` §2.1 | 22 of 33 above chance | **20 of 33** (§5.2 of this file; recomputed twice) |
| `CS_PHASE8_CALLERS.md` §4.1 | contradicts itself 15 lines apart: R5 box says the `neg_probe_rate` exception is gone (p = 0.199); the closing paragraph still asserts it is non-zero at p = 0.020 | pick one — post-C6 it is **p = 0.199** |
| `CS_PHASE8_CALLERS.md` §5.5 | D3 Tier A column 1.248 / 1.237 / 0.935 | **1.288 / 1.256 / 0.907** — the committed `caller_agreement_matched_d3_11sections.csv` was regenerated on post-C6 Tier A |
| `CS_PHASE8_TORUS_VAR.md` §9 | certifies `sf_summary.csv` = `69e3a1d3…`, `summary_phase3.txt` = `ecf86b9c…` as byte-identical to their pins | those are the **committed** hashes; the working tree is `a5ccc9b0…` and `dc92ddc6…` (both 09:06). True when written, stale now |
| `CS_PHASE8_TORUS_VAR.md` §1, §10 | "23 % of shifted senders in the void" | **35.5 %** — §17.8(a). The §10 submission-framing sentence is the one that matters |
| `CS_PHASE8_TORUS_VAR.md` §4 | tiled torus "0.080–0.118" on the irregular window; "up to 2.4× nominal"; raw rejection "0.802" | **0.048–0.118** over that window's four cells (0.080–0.118 only for s ≥ 0.05); **2.35×**; **0.801** |
| `CS_PHASE8_TORUS_VAR.md` §6 | four excluded variance cases "come back positive (+1.05 to +2.01)" | **three** are positive; the fourth is −0.717. And the "22,387–27,567 cells" gloss describes three of the four |
| `PREREG_PHASE8.md:22` | "8.7 is still in flight" | 8.7 completed at **09:08**; the file was written at 09:40 |
| `PREREG_PHASE8.md` | **20 PROVISIONAL marks** citing `a7_summary.csv` @05:19, `sf_summary_c1.csv` @04:41, `summary_phase3.txt` @08-20 | all four files are now **09:06**. Its own action item 6 (`:992`) — "replace every PROVISIONAL mark" — is undone. P3 states the 9–16 % FPR with **no** PROVISIONAL label while its sibling P2 carries one |
| `PREREG_PHASE8.md:886` (P29) | denoising "roughly doubles" the depth loading | raises it by **32–67 %** (§15) |
| `PHASE8_ROADMAP_STATUS.md` | FPR "9–13 %, clean-null subset … 0.127"; caller pooled 1.118; guard "27 of 45"; Tier A vs SenePy "below chance 11/11 (0.914)" | **9–16 %** post-C6 with the clean-null four at 9.1–14.5 %; **1.212** on the frozen basis; guard **52/52**; Tier A vs SenePy **0.972, n.s.** |
| `CS_PHASE7_C1.md` | 160-fit tables; §6.3 presents the N3 degeneracy as a **discovery**; **no attribution for the torus null at all** | the population is now **153**; Mrkvička §2.1.4 predicts the degeneracy, so it is the **first quantification on real tissue**, not a discovery. The file was outside every audit's edit scope and is now the last place the torus finding reads as novel statistics |
| `README.md:285-292`, `:293` | superseded 2-section caller numbers; the 69 % circularity bar | §5.1 and §17.7 |
| `SUBMISSION_PATCH` §7.4, `NOVELTY_ASSESSMENT` U7 | `references.bib` contains no spatial-statistics papers | it does — 11 were added (§17.1) |

**Counts that four different files state four different ways.** Figure-guard
artefacts: 27 / 45 / 46 / **52** (truth: **52**; I ran the guard's own digest
function — 52 manifest entries, 52 in-scope files, 0 drift, 0 missing, 0
unguarded). `references.bib`: 30 / 32 / **43** (truth: 43). `PREREG_PHASE8.md`
length: 770 / 882 / **993**. DCA completed sections: 1 / 2 / **3**.
`COMPLETED_TASKS.md` self-reports **77 entries against 74 numbered rows**.

### 18.4 Items where the reports contradict each other and I have settled it

* **The composition surrogate share.** **This file's §7 is one of the two
  offenders.** §7 says "66–76 % — **untraceable**; do not carry it into the
  pre-registration", while `PREREG_PHASE8.md` P25 has already carried it on the
  grounds that the prohibition lifted. **Two live documents, opposite
  instructions.** Settled: it has a producer —
  `code/run_phase8_compmatch.py` → `results/phase3/compmatch_reruns.csv`, rows
  `type_adj` (**65.9 %**) and `typecomp_adj` (**85.4 %**). What has **no**
  producer is the *pairing* of two different estimators into one range: 66 % is
  1 − SF on β̂; **76 % is a curve-amplitude ratio from `CS_PHASE5.md` §4, a
  different estimator on a different scale, and it is still unreproduced.** The
  row must be **split**, not patched, and the honest bracket is **66 → 85 %**.
  §7 is corrected in this pass.
* **Caller agreement 1.118 vs 1.129 vs 1.212.** The frozen row — 11 sections,
  post-C6 33-gene Tier A, four pairs — is **1.212, z = 21.92, p = 1.84e-106,
  band 0.751–2.198, 35/44 above chance**. 1.129 is the 11-section **pre-C6** row.
  **1.118 appears in no row of `caller_coverage_gate_headline.csv`** — it is a
  three-pair computation, and `PHASE8_ROADMAP_STATUS.md` quotes it as though it
  were the headline.
* **Tier A cross-arm share 26 vs 27.** Recomputed: **27** — 26 via the pinned MGI
  map plus `CDKN2B`, which both frozen sets carry and the map has no row for.
  Always quote it with the convention. `COMPLETED_TASKS.md` row 17 says 26; row 45 of
  the same file already fixes it to 27.
* **A7 `neg_control_probe`.** Two rows of the same table were crossed: the
  p = 0.183 in circulation is the post-C6 **`neg_probe_rate`** p-value, not the
  probe's. See §17.6.

### 18.5 Analysis that is genuinely incomplete

* **A6 spleen covariate — DEFERRED behind the freeze.** The spec is complete
  (`genesets/human/D_spleen_*.txt`); the build needs the expression matrix.
  Roadmap 9.5, and nothing in Phase 9 may begin until the tag exists. Also behind
  the freeze: Tier E3, `DROP_NONSPECIFIC`, and the comparison against the
  depositors' `annotations.csv.gz`.
* **Torus — two items NOT DONE:** the N7 sender axis under N3-var, and
  shift-radius sensitivity. And the covariate-adjusted `*_full_sf` for
  N3-var/N4-var is at **200 permutations, not 1,000**, where every other row now
  meets the §24.3 standard.
* **D2 — the seed panel is one subsample, one section, three seeds.** Enough to
  show the instability happens, not enough to say how often. Two `mor` sections
  (7250, 7260) were OOM-killed and never retried; `mor 7352` completed at
  **09:42**, *after* the `d2_*.csv` tables were written at 09:34, so those tables
  state `mor` coverage as 4 sections when 5 now exist on disk. Regenerate before
  citing `mor` coverage.
* **`B_oxidative_stress` clears the ≥ 30 floor by exactly one gene (31)** — the
  margin gene is `Junb`, from the 100-gene add-on, so the margin is **+1 on the
  authoritative 5,097-gene panel and +0 on the stock-CSV panel definition**. The
  gate must be re-run after **any** gene-set or panel change, not only at freeze;
  it exits non-zero so it cannot pass silently. It has **no reportable fits**
  under either Tier A definition, because its naive amplitude is negative.
* **Human Tier A — status corrected.** The current human Tier A is **33 genes and
  the gate PASSES on both arms**, frozen and primary. What fails is §10's
  *sixteen*-gene proposal: 14 of 16 on-panel human, 13 of 16 mouse, and **five
  survivors — `ATM ATR CDKN2B MDM2 TP53BP1` — not `ATR` alone.** The "ATR alone"
  figure in `COMPLETED_TASKS.md` row 12 is the pre-C6 v1-B7 number, and the failure is
  now attributable to **B4**, not B7. "13 inside Hallmark modules" appears
  nowhere and should not be repeated.
* **Cross-arm symmetry is NOT achievable, and that is quantified, not a gap.**
  `REACTOME_SASP` is **40 mouse against 111 human** genes (verified from the
  pinned MSigDB JSONs, MM14900 / M27187); the two 33-gene Tier A sets share
  **27**.
* **`BIO_DELIVERABLE7_CLAIM_AUDIT.md` is unapplied by construction** — "this is
  an audit, not a fix": **10 BLOCKING, 18 SHOULD-FIX, 13 MINOR, 9 causal-language
  flags**, of which only B7 and C15 carry supersession banners. The highest are
  **B1** (wrong first author on the closest prior work — `Zhao L` is author 14 of
  48), **B9** (CellWHISPER's null is our N1, not N3 — "we did not run their
  null"), **B4** (the README asserts a ">90 % FPR" that `CS_PHASE4` §5.7 says is
  not an FPR).
* **The 7250 anomalies are resolved but not propagated.** All four — the
  DeepScence depth-correlation sign flip, hepatocyte under-calling (0.027 against
  1.38–2.46), SenePy-vs-DeepScence at 2.150 against 0.33–0.55 in ten of eleven
  sections, and the circularity pair overstated ~2× (1.51–2.85× quoted against a
  median 1.071) — have **one cause: an inverted DeepScence anchor sign in section
  7250**. That section was **half the published two-section base.** The published
  sentences this falsifies are listed under "what must be struck" and **have not
  been struck.**
* **`S6: "~1.3M cells"` is untraceable** — the fact-check audit's one unexplained
  provenance failure. The candidates are 1,826,893 / 1,635,937 / 936,125; none is
  1.3M.

### 18.6 Bookkeeping that will embarrass under review

* **The deviation-ID namespace collision is real and unfixed.** Three series:
  `PHASE8_ROADMAP_STATUS.md`'s `D` = **PI decisions** (D5, D15, D16);
  `PREREG_PHASE8_genesets.md` §12's `D` = **gene-set deviations**, D1–D17;
  `PREREG_PHASE8.md`'s `P1`–`P22`, which does not collide. Renumbering is
  **queued, not done** — the gene-set series moves to a `G` prefix before the
  freeze tag. **Until then, always name the series; never cite a bare `D<n>`.**
  *Caveat on the stated scale:* "13 files / 142 references" is asserted only at
  `COMPLETED_TASKS.md`'s "errors caught before they propagated" table and
  **could not be reproduced**; explicit discussion
  appears in 5 files. The collision also occurs **intra-file** — in
  `CS_PHASE8_M1_RERUN.md`, `D1` means the DeepScence coverage task at line 86 and
  a PI decision at line 354.
* **Audit tallies, for the record.** `AUDIT_PHASE8_FACTCHECK.md`: **9 REFUTED
  (R1–R9), 10 TRUE-BUT-MISLEADING (M1–M10), 8 CONFIRMED, 5 DEFERRED, 3
  UNVERIFIABLE** — 35 verdict-bearing items over **34 distinct claims**, because
  M1 is double-filed as both misleading and unverifiable. Two internal
  mismatches: it says "six figure/data chains re-derived" where its own table has
  nine rows, five caveat-free. `AUDIT_CORRECTIONS_APPLIED.md`: **8 fully applied,
  1 partial (M3), 1 disclosed-not-fixed (M7), 1 recommended-not-done (M8), 10
  deferred, 5 disagreed-with** — and its own §0 summary miscounts twice.
  **The corrections agent disagreeing with the audit five times, with evidence,
  is a feature, not a defect**: the audit's own fix list was internally
  inconsistent (e.g. it asked for Tier A `mouse_only` = 7 where 33 − 27 = 6, and
  classed the 69 % literal UNVERIFIABLE where it is REFUTED).
* **Fact-check R6, the item an agent corrected the auditor on.** Its figures were
  computed on the **pre-C6** A7 file of 05:19 and superseded by the re-run at
  06:52. On the frozen 825 control fits the reportable filter admits **4.8 % on
  the full design — essentially nominal, and identical across all five
  families**. The substance holds (the "2–3× nominal" bound is on the *estimator*,
  not the *filter*, and that sentence is withdrawn) but neither audit file
  carries a supersession banner.
* **`COMPLETED_TASKS.md` is not a picture of open work by construction** — its
  own closing line says "in-flight items are NOT listed here" — and several rows
  (12, 17, 23) are stale **in place**, superseded by later rows in the same file
  without being amended. It is also still being appended to as this is written,
  so **it must be cited by row or section, never by line number**; this ledger
  does that throughout for exactly that reason.
