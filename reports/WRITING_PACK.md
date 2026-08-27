# WRITING PACK — SASP-kernel paper

**Created 2026-08-27.** **There is no manuscript.** Nothing has been drafted. This file is
the reference a person writes the paper *from*: every number the paper needs, organised by the
§30 outline, each with its value, its source file, and the command that produces it; the
claims that may never be written; the framing decisions already settled; and what is not yet
measured.

It **supersedes the framing** of `reports/SUBMISSION_PATCH_2026-08-29.md`, which was written as
a patch against a draft that does not exist. Its *content* has been mined into this file and
re-verified against `results/` where possible.

**Method note.** Every number below is marked:

- **[V]** — I re-derived it from `results/` in this session and it reproduced. The derivation
  command is given.
- **[F]** — read directly from a named result file, not independently re-derived.
- **[UNSOURCED]** — appears in project documents but **cannot be traced to any file**, and I
  could not re-derive it. **Do not put an UNSOURCED number in the paper.**

---

# 0. LEAD: what is still UNSOURCED

This project has twice shipped numbers that turned out to be recalled rather than measured.
These are the ones that would embarrass the paper now.

## 0.1 λ̂ = 15.7 µm "pooled" — **UNSOURCED — WITHDRAWN 2026-08-27; the value is 14.7 µm**

**Where it appears:** `SASP_Kernel_Master_Plan.md` §30 5.3 ("λ̂ ≈ 15.7 µm pooled, so ≈ 6λ") and
§29 objection 7; `reports/CS_PHASE8_TORUS_VAR.md` L105, L291, L323; `reports/CORRECTIONS.md`
L637; `reports/SUBMISSION_PATCH_2026-08-29.md` L481, L502.

**Status.** No file emits it. I could not reproduce 15.7 from `results/phase3/main_fits.csv`
under any subset I tried. Its only apparent provenance is a back-derivation from the torus
report's own sentence "2,215 µm … 141× the pooled λ̂" (2215/141 = 15.71), i.e. it is circular
with the claim it supports.

> ## ✅ RESOLVED 2026-08-27 — **λ̂ = 14.7 µm**, and 15.7 is withdrawn
>
> **The authoritative definition is the pooled median of `lam_naive` over the 315 primary
> fits** (in-band × `tierA_p95` × `stratum == "all"`) — **14.7 µm** (14.7321). It is printed
> by `code/summarize_phase3.py:221` into `results/phase3/summary_phase3.txt` §6, the
> `tierA_p95` row, column `medlam`. **IQR [7.0, 50.0] µm; 60 % of fits railed** — that caveat
> is not optional.
>
> **Why this and not one of the other three** (full reasoning in
> `reports/RECORD_RECONCILIATION.md` §1):
> 1. **It is emitted.** 16.07 and 14.99 are emitted by no file. Replacing an invented number
>    with an unemitted one is the same failure one step later.
> 2. **The pre-registration freezes the summariser, not a λ̂ estimand.** `PREREG_PHASE8.md`
>    §5 explicitly declines to make a length constant the estimand — *"λ̂ rails at a grid
>    bound in a majority of M1 fits, so a fitted length constant is not the estimand"* — and
>    its benchmark table quotes only the **railing rate** for λ̂. There is no pre-registered
>    λ̂ to match, so the tiebreak falls to what the frozen code prints.
> 3. **It shares a denominator with its mandatory caveat.** The 60 % railing rate is over the
>    315. Quoting a λ̂ over 153 beside a railing rate over 315 mixes bases.
> 4. **It is a valid median of a censored sample; the interior median is not.** 103 fits are
>    censored at the 7 µm floor and 86 at the 50 µm ceiling, but the median order statistic
>    falls in the **interior** (the 55th of the 126 unrailed values), so 14.73 is a legitimate
>    median. The interior median discards 60 % of the sample non-randomly, at both ends
>    asymmetrically, and is biased upward — and it is unstable to the population choice
>    (17.06 over 315 vs 14.99 over 153) in a way the pooled median is not.
> 5. **"Pooled" becomes true.** The phrase in circulation is "the pooled λ̂"; 14.73 is
>    literally the pooled median. 15.7 was never pooled anything.
>
> **Where 15.7 probably came from.** The closest reproducible match anywhere in either tree
> is **15.716 µm** — the *interior* median over the **pre-C6** `tierA_p95` in-band fits
> **including the zonation-stratified rows** (441 rows, so hepatocytes are counted four times
> over). Pre-C6, interior, and pseudo-replicated: not a defensible estimand under any
> reading, so 15.7 is withdrawn either way.
> ```bash
> python3 -c "
> import pandas as pd
> d=pd.read_csv('results/phase3_pre_c6/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
> p=d[(d.call=='tierA_p95')&(d.sec.isin(['7259','7260','7001','7248','7352','7435']))]
> print(p[p.lam_railed==0].lam_naive.median(), len(p))"   # -> 15.716  441
> ```
>
> **Dependents, re-derived** (all move in the direction that *strengthens* the claim):
> N3-var's 2,215 µm displacement is **150×** λ̂ (was 141×); at a 1,200 µm tile side the seams
> are **~81 λ̂** apart (was ~76); the 100 µm window spans **≈ 6.8 λ̂** (was "≈ 6λ"), i.e. 14×
> the 7 µm floor but only **2×** the 50 µm ceiling — say that too. N4-var's 3,395 µm is 230×
> λ̂; N3-occ's 28 µm is 1.9× λ̂.

**Verified alternatives, pick one and name the estimand [V]:**

| definition | value | source |
|---|---|---|
| median λ̂_naive over the **153 primary reportable fits** | **16.07 µm** | `results/phase3/main_fits.csv` |
| median λ̂_naive over all **315** primary fits | **14.73 µm** | same |
| **interior** median (railed fits excluded) over 315 — *the value the summariser prints* | **17.1 µm** | `results/phase3/summary_phase3.txt` §1 |
| interior median over the 153 reportable | 14.99 µm | `results/phase3/main_fits.csv` |
| median λ̂ in the N7 table row `tierA_p95` | 14.7 µm | `results/phase3/summary_phase3.txt` §6 |

IQR is **[7.0, 50.0] µm** for every one of these — i.e. **60 % of fits are railed at a grid
bound** (`m1_final_audit.txt`, `lam_railed_frac = 0.6000` [F]). That fact must travel with any
λ̂ the paper quotes; "pooled median" and "interior median" are materially different quantities.

```bash
python3 - <<'PY'
import pandas as pd, numpy as np
d=pd.read_csv('results/phase3/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
inband=['7259','7260','7001','7248','7352','7435']
p=d[(d.call=='tierA_p95')&(d.sec.isin(inband))&(d.stratum=='all')]
r=p[(p.beta_naive>0)&(p.beta_base_lo>0)]
print(len(p), len(r), p.lam_naive.median(), r.lam_naive.median(),
      p[p.lam_railed==0].lam_naive.median(), np.percentile(r.lam_naive,[25,75]))
PY
```

**Consequences — now recomputed.** The two dependent claims become **"N3-var displaces 150×
the pooled λ̂"** and **"at a 1,200 µm tile side the seams are ~81 λ̂ apart"**
(`CS_PHASE8_TORUS_VAR.md` §2, §4, both corrected 2026-08-27). Neither reverses; both get
slightly stronger.

## 0.2 The "76 %" composition-surrogate share — **UNSOURCED denominator**

`reports/CS_PHASE5.md` §4 L71–72, L384: *"the unstratified contact amplitude is +0.260
response-sd; the composition-only surrogate reproduces +0.212 sd, i.e. 76 % of it."*
0.212 / 0.260 = **0.815**, not 0.76. The only two "ALL RECEIVERS" amplitudes on disk are
**0.260** (`results/phase5/summary_phase5.txt` T3, identical in `results/phase5_pre_c6/`) and
**0.266001** (`figures/figure2a_amplitudes.csv`) [V]. Neither yields 76 %. A denominator of
≈ 0.279 would, and it is in no file.

**This row is superseded anyway.** §17's "composition surrogate share 66–76 %" is replaced by
the composition-matched decomposition (§ 5.6 below): **65.9 % (receiver's own cell type) to
85.4 % (cell type + 20-NN composition)**, both from `results/phase3/compmatch_reruns.csv` [V].
**Do not carry 66–76 % forward.**

## 0.3 CoreScence mouse circularity "24/35 = 69 %" — **UNSOURCED denominator** (already known)

`AUDIT_PHASE8_FACTCHECK.md` M1: the 69 % is a typed-in literal in
`code/make_figure_genesets.py:245` and a **string** in `code/gate_disjointness_human.py:291`.
Numerator 24 reproduces; denominator 35 does not. Corrected, now derived by
`code/corescence_circularity.py`: **mouse 79 % → 88 %, human 76 % → 88 %**. The real cost of
C6 is about half what the original narrative said. Human halves reproduce exactly
(25/33 = 75.8 %, 29/33 = 87.9 %); the mouse bar in `figures/figure_gs3_*` came from the
literal. **Never write "69 % → 76 % → 88 %".**

## 0.4 Spearman ρ for Moran's I vs A7 — **aggregation-dependent; use +0.895 and say which**

> ## ✅ RESOLVED 2026-08-27 — four aggregations, all reproduce, the falsification survives all
>
> +0.923 **is** reproducible — under a different aggregation. All four re-derived:
>
> | aggregation | ρ | p | emitted by |
> |---|---|---|---|
> | section-clustered **mean** per field, knn6 **raw**, 12 control+module fields | **+0.8951** | 8.37e-05 | `moran_verdict.txt`; `code/summarize_moran.py:183` — **frozen; quote this** |
> | section-clustered **mean** per field, knn6 **cell-type-centred**, 12 fields | **+0.9441** | 3.93e-06 | same, `:184` — **frozen** |
> | **median** per field, knn6 raw, 12 fields | +0.9231 | 1.86e-05 | no file |
> | **per-row**, no aggregation, 132 pairs | +0.7104 | 1.43e-21 | no file |
>
> **Every value is positive and significant, so "Moran's I and the A7 kernel do not disagree"
> holds at every aggregation. Say that explicitly** — otherwise a reader who recomputes it a
> different way will think the finding is fragile. Only the digit is aggregation-dependent.
> **The aggregation must be stated in the same clause as the number, everywhere.**

`SUBMISSION_PATCH_2026-08-29.md` §9 L606 says "+0.923 raw". I recomputed from **both**
`results/moran/moran_pooled.csv` and `results/moran/moran_vs_a7.csv` field means and got
**ρ = +0.8951, p = 8.37 × 10⁻⁵** [V] — matching `results/moran/moran_verdict.txt`.
The centred value **+0.9441, p = 3.93 × 10⁻⁶** is agreed by every source [V].
**Authoritative file: `results/moran/moran_verdict.txt`.**

## 0.5 The detectable bound's distance range is an editorial join, not a file value

The bound **0.1833** is in `results/phase3/m1_final_audit.txt` as `power80_bound` [F,V]. No
file attaches a distance range to it. "**over 0–100 µm**" (§30 5.3) is correct as an inference
— every fit is capped at `WINDOW_UM = 100.0` (`code/run_phase3_nulls.py:59`, enforced at
`:169`) — but it is the writer joining two facts. Say it that way; do not cite a file for the
range.

## 0.6 The 3-pair caller-agreement numbers are emitted by no file

`1.030 / 1.118 / 1.128 / 1.212` (3-pair basis) appear in `PHASE8_ROADMAP_STATUS.md` and
`SASP_Kernel_Master_Plan.md` §30 5.9 but in **no CSV**: both copies of
`caller_coverage_gate_headline.csv` carry `n_pairs = 4` only. **I re-derived all four exactly
[V]** by Mantel–Haenszel pooling of the significance CSVs (command in §5.9 below). Either use
the 4-pair values, which *are* file-emitted, or cite the significance CSVs plus the pooling
formula. Never mix bases inside one sentence (see §5.9).

## 0.7 Bracket type: SF and amplitude brackets are **IQRs, not CIs**

`results/phase3/sf_summary.csv` has columns `q25, median, q75` — **there is no CI column** [V].
So `0.088 [−0.017, 0.234]` and `0.029 [−0.007, 0.084]` are **inter-quartile ranges across
fits**, not confidence intervals. `CS_PHASE8_M1_RERUN.md` §6 says so explicitly ("HOLDS — the
bracket is the IQR, not a CI — label it"). §30 5.2/5.3 currently present them unlabelled.
**Label every such bracket "IQR across fits".** Genuine CIs exist only for the
composition-matched SFs and the A7 clustered means.

> ## ✅ RESOLVED 2026-08-27 — they are IQRs; the pre-registered bootstrap cannot make them CIs
>
> **The bracket has no bootstrap in it at all.** `code/summarize_phase3.py:99` is
> `np.quantile(v, [.25, .5, .75])` over the per-fit SF **point estimates**, written to
> `sf_summary.csv` as `q25 / median / q75`. `m1_final_audit.txt` names the two headline
> brackets `ctrl_amp_iqr` and `SF_N2+N5+N6_iqr`. **The §3.6 bootstrap — 400 replicates over
> 100 quantile blocks — emits per-fit CIs only** (`sf_n2n5n6_lo/hi` in `main_fits.csv`, whose
> median span is [−0.415, +0.381], two orders of magnitude wider than the IQR). **It emits no
> interval on the median across fits, so a genuine CI cannot be computed without a new run**
> — and `results/` is frozen. **Therefore: relabel everywhere, do not fabricate a CI.**
>
> `PREREG_PHASE8.md` §5 calls it a *"paired-bootstrap interquartile range"* and §6 R1 writes
> the replication criterion on a *"paired-bootstrap interval"*. Both are misnomers and are
> corrected by dated note in that file's new §0.0 (item C-4). **R1 reads: the IQR across the
> reportable fits includes 0 and its upper quartile is below 0.50.** M1's own outcome is
> unchanged.
>
> Relabelled in: `SASP_Kernel_Master_Plan.md` §29/§30, `README.md`,
> `SUBMISSION_PATCH_2026-08-29.md` §8, `PLAN_UPDATE_D12_D13.md`,
> `PHASE8_ROADMAP_STATUS.md`, `BIO_DELIVERABLE6_DISCUSSION.md`, `PREREG_PHASE8.md` §0.0.

---

# 1. THE AUTHORITATIVE TREE, AND THE TWO VINTAGES

Almost every disagreement in this repo is one thing: **pre-C6 vs post-C6**.

| | pre-C6 (superseded) | **post-C6 (FROZEN, authoritative)** |
|---|---|---|
| Tier A strict | 25 genes | **33 genes** |
| results tree | `results/phase3_pre_c6/`, `results/phase5_pre_c6/` | **`results/phase3/`, `results/phase5/`, `results/moran/`, `results/phase8_d2/`** |
| A7 files | 05:19–05:20 | **09:06** |
| SF / summaries | — | **09:06** |
| `main_fits.csv` | — | **06:13** |
| git tag of the old gene sets | `pre-c6-genesets` @ `e789372`, 05:41:13 | — |

**Rule for the writer: if a number came from a document rather than from `results/phase3/`,
`results/phase5/`, `results/moran/` or `results/phase8_d2/`, check its vintage before using
it.** The following documents are wholly or partly pre-C6 and their digits must not be copied:
`reports/CS_PHASE7_C1.md` (all SFs), `reports/CS_PHASE8_CALLERS.md` (every Tier-A-involving
number), `reports/NOVELTY_ASSESSMENT.md` (all A7 and torus numbers),
`reports/AUDIT_PHASE8_FACTCHECK.md` R1/R5/C5/C6 digits, and
`reports/PHASE8_ROADMAP_STATUS.md`'s 8.4 gate table, 8.5b A7 block and "N3-tile 0.974" row.

`PHASE8_ROADMAP_STATUS.md`'s **top-of-file banner** (0.029 / 0.183 / 0.088) is post-C6 and
correct; its **PI-decisions table** (0.326 / 0.027 / 0.082 / 0.203) is pre-C6 and stale. Both
are in the same file.

---

# 2. NUMBERS BY §30 OUTLINE SECTION

## 2.0 Cohort and panel descriptors (Methods §4, and every Results section)

### M1 — mouse liver, GSE310392, Xenium Prime 5K

| quantity | value | source | mark |
|---|---|---|---|
| sections / donors | **11 / 11** (one section per mouse) | `results/composition_by_arm_timepoint.csv` | V |
| arms | **6 SBR / 5 sham** | same | V |
| timepoints | **2, 10, 26, 52 weeks**; all Male | same | V |
| total cells | **1,826,893** | same (`n_cells.sum()`) | V |
| analysable (QC-pass **and** cell-type-labelled) | **1,635,937** | same (`n_analysable.sum()`) | V |
| **in-band (6 admissible sections), cells** | **1,031,880** | same, restricted to the 6 | V |
| per-section cells | 83,392 – 236,905 | same | V |
| §17 as written — **wrong, correct it** | 1,834,806 / 1,036,459 | — | — |
| admissible sections (Test 3, `Cdkn1a`⁺ hepatocyte prevalence 1–20 %) | **7259, 7260, 7001, 7248, 7352, 7435** | `results/phase3/summary_phase3.txt` header | F |
| excluded, over the 20 % ceiling | 7239, 7448, 7361, 7450 | same | F |
| excluded, below the 1 % floor | 7250 (0.48 %) | same | F |
| median nearest-neighbour distance | **6.74 – 10.61 µm** over all 11 | `results/phase3/window.csv` (`median_nn_um`) | F |
| §17 as written — **wrong** (SBR-only scoping range) | 6.7 – 9.7 µm | — | — |
| transcript assignment rate (7259) | 88.27 % | `code/assignment_rate.py` | F |
| **panel** | h5 Gene Expression **5,106** = stock 5,006 ∪ add-on 100, disjoint; minus 9 non-ENSMUSG genotyping probes ⇒ **5,097 authoritative** | `cell_feature_matrix.h5`; audit C2 | F |
| control features | 40 negative-control probes / 609 negative-control codewords / 21 genomic controls | audit C2, `a7_verdict.txt` | F |
| **sender prevalence, `tierA_p95`, in-band** | **4.04 – 4.48 %**, mean **4.29 %** | `results/phase3/window.csv`, `summary_phase3.txt` §6 | F |
| prevalence across the 6 sender calls | 0.63 – 8.58 % | `summary_phase3.txt` §6 | F |

```bash
python3 -c "
import pandas as pd; c=pd.read_csv('results/composition_by_arm_timepoint.csv')
print(len(c), c.mouse.nunique(), c.arm.value_counts().to_dict(), c.n_cells.sum(), c.n_analysable.sum())
inb=[7259,7260,7001,7248,7352,7435]; print(c[c.section.isin(inb)].n_cells.sum())"
```

### Gene sets (Methods; Appendix)

| quantity | value | source | mark |
|---|---|---|---|
| Tier A strict, on panel | **33 / 33** both arms (pre-C6: 25) | `genesets/A_SENDER_FINAL_strict.txt` (33 lines) | V |
| Tier B modules, on panel (mouse) | `tnfa_nfkb_proximal` 126, `il6_jak_stat3` 68, `interferon_response` 100, `downstream_arrest` 190, `emt_ecm` 125, **`oxidative_stress` 31**, `secondary_senescence` 108 | audit C1; line counts of `genesets/B_*.txt` agree | V |
| all seven Tier B ≥ 30 | yes; `B_oxidative_stress` = **31**, margin of **exactly 1 gene** (`Junb`, from the add-on) | audit C1/C2 | F |
| **Tier A ∩ ⋃Tier B** | **0**, all 7 modules, both arms — gate PASS | `code/gate_genesets_guard.py` | F |
| **but** the seven Tier B modules are **not mutually disjoint** | **18 of 21 module pairs share genes** | `figures/figure_gs1_intersection_matrix_data.csv` | F |
| per-module sender sets (`tierApm`) | mouse 70/74/73/37/73/71/74; all ≥ 15 and disjoint from their own module | audit C1 | F |
| B7 `secondary_senescence` under C6 | mouse 38 → **108** | `git diff pre-c6-genesets -- genesets/` | F |
| cross-arm Tier A overlap | **27 of 33** (26 by the pinned MGI map + `CDKN2B`, which both arms carry and the map has no row for) | audit R2 | F |
| cross-arm B7 | **88** of 108/116 on the ortholog-intersected panel | `figures/figure_gs2_crossarm_symmetry_data.csv` | F |
| CoreScence circularity, post-C6 | mouse **88 %**, human **88 %** (pre-C6: mouse 79 %, human 76 %) | `code/corescence_circularity.py` | F |
| SenePy spleen coverage | 22 spleen cell types: **0 matched / 15 surrogate / 7 none** | `figures/figure_gs4_senepy_coverage_data.csv` | F |

### H1 — human spleen, GSE326743 (**acquired and structurally verified only; not analysed**)

| quantity | value | source | mark |
|---|---|---|---|
| donors | **7** (SPLN07/14/21/24/30/43/44) | `reports/PHASE7_H1_SCREEN.md` §3 | F |
| ages / sex | 17, 31, 32, 32, 37, 57, 59; **4 M / 3 F**; two donors over 55 | same | F |
| total cells | **2,207,593**; per sample 220,435 – 396,173 | same | F |
| panel, **verified on the data in all 7 samples** | **5,093** Gene Expression (all ENSG, 0 duplicate symbols) | `genesets/h1_candidate/GSE326743_gene_panel_5093.csv`; audit C3 | F |
| control complement | 609 codewords / 40 probes / 21 genomic — identical to M1 | audit C3 | F |
| screen scope | 132 GPL33762 series screened; 19 with real Prime-5K evidence; GSE326743 the **only normal-tissue Prime-5K deposit with ≥ 3 donors** | `reports/PHASE7_H1_SCREEN.md` | F |
| **do not write** | "the human per-sample range brackets the mouse range" — it sits **above** it and overlaps only at the top (audit M10) | — | — |

---

## §5.1 — Identifiability regimes on synthetic tissue (Figure 1)

Grid: **20 cells** = 4 clustering levels (κ = 0, 1.5, 3.0, 4.5) × 5 autocorrelation ratios
(ℓ/λ_true = 0.25, 0.5, 1, 2, 4), 30 seeds each.
**Source: `figures/figure1_data.csv` (20 rows × 68 cols).** Producer `code/make_figure1.py`;
sweep `code/sasp_sweep.py`; rolled up in `results/summary_tables.txt`.

| quantity | value | mark |
|---|---|---|
| **mean λ̂ CI coverage where ℓ ≥ 2λ (8 of 20 cells)** — naive | **0.5125** | V |
| — matched-decoy | **0.346** | V |
| — nuisance-conditioned | **0.854** | V |
| mean coverage over all 20 cells — naive / decoy / nuisance | 0.670 / 0.630 / 0.702 | V |
| naive **iid** CI coverage, all 20 | 0.410 | V |
| matched-decoy \|relative bias in λ̂\| beats naive in | **12 of 20** grid cells | V |
| within ℓ ≥ 2λ, matched-decoy coverage **strictly worse** than naive in | **5 of 8** (6 of 8 if the κ=0, ℓ/λ=4 tie at 0.000 counts) — **CORRECTIONS.md §13.1 says "6 of 8"; strictly it is 5** | V |
| block sd / iid SE inflation | 1.46 – 7.37 across the grid | F |
| clean-regime sanity check | "pure" regime coverage 0.95 naive / 0.95 block; "easy" 0.80 / 0.60 | F |

```bash
python3 - <<'PY'
import pandas as pd
d=pd.read_csv('figures/figure1_data.csv'); s=d[d.ell_over_lambda>=2]
print(len(s), s.cover_lam_naive_blk_mean.mean(), s.cover_lam_decoyS_blk_mean.mean(), s.cover_lam_nuis_blk_mean.mean())
print("decoy relbias better in", (d.relbias_lam_decoyS_mean.abs()<d.relbias_lam_naive_mean.abs()).sum(), "of", len(d))
PY
```

**Why this matters:** §5.1 is where §5.6 is first visible, **on planted ground truth**. Carry
panels a and b, not only c.

---

## §5.2 — The naive gradient and its behaviour under nulls (Figure 2a–d)

**Primary analysis population, frozen:** in-band 6 sections, sender call `tierA_p95`,
receiver type × Tier B module, stratum `all` ⇒ **315 fits**; **reportable** =
`beta_naive > 0 AND beta_base_lo > 0` (`code/summarize_phase3_c1.py::reportable`) ⇒
**153 fits** [V].

| quantity | value | source | mark |
|---|---|---|---|
| **naive amplitude** (median β̂/sd(y) over reportable fits) | **0.3288 → 0.329** | `results/phase3/m1_final_audit.txt`, `main_fits.csv` | V |
| naive median \|β\|/sd over **all 315** | 0.238 | `summary_phase3.txt` §1 | V |
| naive β̂ > 0 in | 216 of 315 | same | V |
| **surviving fraction under N2+N5+N6** | **0.0885 → 0.088**, IQR **[−0.017, 0.234]** | `results/phase3/sf_summary.csv` | V |
| λ̂ railed at a grid bound | **60 %** (103 at the 7 µm floor, 86 at the 50 µm ceiling) | `summary_phase3.txt` §1 | F |
| pre-C6 comparators (**do not quote**) | naive 0.3262, SF 0.0822 [−0.099, 0.249] | `m1_final_audit.txt` | V |

### The full null battery, primary population (n = 153)

**Source: `results/phase3/sf_summary.csv` / `results/phase3/summary_phase3.txt` §1.**
Brackets are **IQRs**. `frac ≤ 0` and `frac > 0.5` are the two extra columns.

| null | median SF | IQR | ≤0 | >0.5 |
|---|---|---|---|---|
| N1 stratified label permutation | **0.707** | [0.416, 0.865] | 0.04 | 0.69 |
| N1 on the N5+N6-conditioned residual | 0.989 | [0.897, 1.092] | 0.04 | 0.95 |
| **N2 matched decoy** | **0.952** | [0.919, 0.978] | 0.00 | 1.00 |
| N3 torus shift (published bounding box) | **0.999** | [0.989, 1.006] | 0.00 | 1.00 |
| N4 rotation (published bounding box) | **0.947** | [0.804, 1.039] | 0.02 | 0.92 |
| **N5 nuisance conditioning** | **0.115** | [−0.034, 0.258] | 0.29 | 0.10 |
| N6 receiver-baseline conditioning | 0.471 | [0.183, 0.742] | 0.10 | 0.45 |
| zonation covariate alone | 0.843 | [0.555, 0.974] | 0.06 | 0.78 |
| N5 + N6 | 0.084 | [−0.013, 0.233] | 0.30 | 0.08 |
| **N2 + N5 + N6 (headline)** | **0.088** | [−0.017, 0.234] | 0.30 | 0.08 |
| N8 scrambled response gene set | 0.916 | [0.817, 1.070] | 0.01 | 0.97 |

N3/N4 corrected variants are in §5.7. N7 (sender axis) is in §5.3.

**Producers.**
```bash
python3 -u code/run_phase3_nulls.py --stage window  --sections all
python3 -u code/run_phase3_nulls.py --stage main    --sections all --calls all9 --n-jobs 24
python3 -u code/run_phase3_nulls.py --stage perm    --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
python3 -u code/run_phase3_nulls.py --stage perm_c1 --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
python3 -u code/summarize_phase3.py        # -> sf_summary.csv, summary_phase3.txt
python3 -u code/summarize_phase3_c1.py     # -> sf_summary_c1.csv, summary_phase3_c1.txt
python3 -u code/m1_final_audit.py          # -> m1_final_audit.txt  (the headline vector)
```
**Every permutation null in the frozen tree is at n_perm = 1000** — verified in
`m1_final_audit.txt` §1 for `perm_nulls{,_c1,_n7,_c1_n7,_pm,_c1_pm,_var}.csv`; the sole
exception is `perm_nulls_var_full200.csv` (the covariate-adjusted N3-var/N4-var family, 200)
[F]. Pre-C6 N7 files were at 200, which makes the **N7 pre/post comparison confounded**
(gene sets and permutation count moved together); only `tierA_p95` is clean on that axis.

### Where the naive gradient comes from (`summary_phase3.txt` §11)
- Binned curve monotone decreasing in **32 of 42** section × module fits.
- SF after adding **only receiver cell-type intercepts**: median **0.348**.
- SF after the full N5 block: median **0.042**.

### Ripley / Poisson geometry (`summary_phase3.txt` §10, `results/phase3/poisson_fits.csv`)
- log(median distance to nearest sender) on log(sender density): slope **−0.525**, **r² 0.9843**
  (homogeneous Poisson predicts exactly −0.5); obs/Poisson median distance ratio **1.042** [F].
- Ripley ratio median **1.1237** [F] (`m1_final_audit.txt`).
- λ̂ is **not** a density readout: r² of log(λ̂_raw) on log(sender density) median **0.016**,
  slope median **−0.091** against −0.500 for a pure density readout [F].

---

## §5.3 — Controlled kernel estimates (Figure 3)

| quantity | value | source | mark |
|---|---|---|---|
| **controlled amplitude** (median β̂_N2N5N6 / sd(y), reportable fits) | **0.0288 → 0.029** | `results/phase3/main_fits.csv`, `m1_final_audit.txt` | V |
| its bracket **(IQR, not a CI)** | **[−0.0066, +0.0845]** | `m1_final_audit.txt` (`ctrl_amp_iqr`) | V |
| median SE of the controlled amplitude | **0.0654** | `m1_final_audit.txt` (`ctrl_amp_se_med`) | F |
| **detectable bound at 80 % power** | **0.1833 → 0.183** = 2.802 × 0.0654 (z_{0.975}+z_{0.80}) | `m1_final_audit.txt` (`power80_bound`) | V |
| ratio, bound ÷ estimate | **6.4×** (pre-C6: 7.5×) | `CORRECTIONS.md` §A | V |
| controlled fits positive **and** CI excludes 0 | **13 of 153** (pre-C6: 15 of 160) | `m1_final_audit.txt` | F |
| distance range | fits are hard-capped at **100 µm** (`WINDOW_UM`, `run_phase3_nulls.py:59`, enforced `:169`) | code | F |
| pre-C6 comparators (**do not quote**) | 0.027, bound 0.203 | `m1_final_audit.txt` | V |

```bash
python3 - <<'PY'
import pandas as pd, numpy as np
d=pd.read_csv('results/phase3/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
p=d[(d.call=='tierA_p95')&(d.sec.isin(['7259','7260','7001','7248','7352','7435']))&(d.stratum=='all')]
r=p[(p.beta_naive>0)&(p.beta_base_lo>0)]
print("naive", np.median(r.beta_naive/r.sd_y), "ctrl", np.median(r.beta_n2n5n6/r.sd_y))
print("SF", np.median(r.sf_n2n5n6), np.percentile(r.sf_n2n5n6,[25,75]))
PY
```

### Justifying the 100 µm window (needed in the same sentence as the bound)
- **Frozen as a literal, not a runtime rule** (PI decision D16, `PREREG_PHASE8.md` §3.2). The
  "99th percentile of distance-to-nearest-sender" is the *provenance* of the number 100, not
  something the code computes.
- Measured 99th percentiles, `results/phase3/window.csv` [F]: `tierA_p95` in-band
  **76.0–112.1 µm** (median over all 11: 96.1); `cdkn1a_pos` 76.3–160.9; `senepy_p95`
  **118.3–186.8 µm** — all six in-band SenePy values **exceed the cap**.
- Share of receivers the cap discards, in-band: `tierA_p95` **0.8–7.0 %** beyond 80 µm and
  0.00–0.11 % beyond 150 µm; `senepy_p95` **7.4–21.5 %** and 0.2–2.3 % [F]. A SenePy-called
  kernel is fitted on a materially truncated distance distribution (deviation P8).
- External calibration: the endometrium nearest-neighbour figures **45–211 µm**
  (`references.bib`, §31 ref 7). **Retrieved, not re-derived.**
- λ̂ for the window-justification claim: **RESOLVED — λ̂ = 14.7 µm pooled, so the window is
  ≈ 6.8 λ̂** (`results/phase3/summary_phase3.txt` §6, `tierA_p95`, `medlam`). Write it with the
  railing caveat in the same sentence: IQR [7.0, 50.0] µm, 60 % railed, so the window is 14×
  the 7 µm floor and **2×** the 50 µm ceiling. See §0.1.

### λ grid, frozen (`PREREG_PHASE8.md` §3.1, `code/run_phase3_nulls.py:59-93`)
`WINDOW_UM = 100.0`; floor **7.0 µm** (resolution floor, a literal that does **not** adapt per
section); ceiling **window/2 = 50.0 µm**; **40 log-spaced points**. "Railed" = exact index
equality at either end (`:239`). Three other non-equivalent railing definitions exist in the
codebase — any railing rate in the paper is the Phase-3 one.

### Kernel families (Figure 3; `results/phase5/kernel_families.csv`, `summary_phase5.txt` T4)
Five families, frozen: exponential, gaussian, powerlaw (p ∈ {0.5,1,2,4}), step, spline
(6 knots). Selection rule: **strict argmin AIC**, no ΔAIC band, no within-section CV.

| design | AIC win rate | beats the no-kernel (covariates-only) model |
|---|---|---|
| naive | step 0.508, spline 0.248, exponential 0.121, gaussian 0.089, powerlaw 0.035 | step 0.898, gaussian 0.778, exponential 0.743 |
| **ctrl (the reportable one)** | **step 0.898**, spline 0.054, gaussian 0.032, exponential 0.013, powerlaw 0.003 | **step 0.511, gaussian 0.241, exponential 0.197, spline 0.152, powerlaw 0.130** |

d̂½ spread across families within a fit: **ctrl median 6.00×**, naive 3.20× [F].
Held-out (leave-one-section-out, 252 folds), ctrl: **spline and step have negative median
ΔLL vs covariates-only**; no family beats it in more than 54 % of folds [F].
**Under control, no kernel family earns its place.**

### Superposition vs nearest sender (`results/phase5/summary_phase5.txt` T1)
- Synthetic calibration: correct verdict **1.00** in all four planted/confound conditions.
- Real data, ctrl design, n = 315: ΔAIC/1k median **−0.109**, superposition wins **0.762**;
  paired block-bootstrap win fraction median 0.730, **decisive (≥0.95) in only 0.130**,
  decisive for nearest in **0.000**.
- Amplitude bound, controlled: superposition median **+0.0076** IQR [−0.0163, +0.0185];
  nearest **−0.0029**. Naive: **+0.0277** and **+0.1472**.
- Does either kernel beat covariates-only under ctrl? nearest ΔAIC **+2.67** (improves 0.20 of
  fits); superposition **+0.74** (0.45). **No.** [F]

### Winner's curse (`summary_phase5.txt` T2)
Cross-fitting at selection fraction 0.50: placebo-corrected winner's curse **+0.0513** of
surviving fraction; the in-sample SF(N5) is **inflated** by selection. Correcting for it makes
the negative **stronger**. SF(N5) also depends on sample size: 100 blocks **0.144**, 50 blocks
0.200, 25 blocks 0.278 [F].

---

## §5.4 — Existing tools under the same nulls (Figure 4)

**Source: `results/phase4/headline.csv`** (100 rows; the `pair == 'ALL'` rows below).
Producers `code/phase4_run.py` → `phase4_summarize.py` → `code/make_figure4.py`.
`make_figure4.py` was **de-hardcoded**: `load_ours()` re-derives all seven previously literal
values and the regenerated PNG is byte-identical (`figure4.png` md5
`000f34051112aff4fed293fe7a5b25c2`) [F].

Significant-interaction **survival rate** under each null (fraction of real-significant
interactions still significant on nulled data), ALL pairs:

| method | N0_perm (=N1) | N0_type | N3_lig | N3_type | N4_lig |
|---|---|---|---|---|---|
| COMMOT | 0.794 | 0.811 | 0.806 | 0.769 | 0.820 |
| **CellChat v2\*** | **0.945** | **0.971** | **0.974** | 0.913 | **0.975** |
| NCEM linear\* | 0.004 | 0.009 | 0.702 | 0.080 | 0.701 |
| SpaTalk\* | 0.791 | 0.781 | 0.814 | 0.784 | 0.820 |

Real-significant rate on real data: COMMOT 0.224, CellChat 0.283, NCEM 0.0144, SpaTalk 0.247.
Score surviving fractions (median) are in the same file. [F]

**Citation discipline for this section (mandatory):**
- CellWHISPER's **">90 % FPR" is an inferred interaction-count ratio, not a measured type-I
  error against a nominal level.** Use "implying", never "measured".
- **CellWHISPER's null is a within-cell-type location permutation = this project's N1.**
  It contains **no torus shift**; the strings *torus*, *toroidal*, *wraparound* do not occur in
  the paper (`CITATION_AUDIT.md` §0.1). The torus null is **Lotwick & Silverman (1982)**.
  This was corrected in nine places; **`reports/CS_PHASE7_C1.md` §6.3 is the last surviving
  place where the torus finding is still presented as novel statistics.**

---

## §5.5 — The estimand's own estimator, refit on the assay's negative controls (Fig. 2h)

**Source: `results/phase3/a7_summary.csv` and `a7_verdict.txt`, both 09:06 (FROZEN).**
825 control fits (165 × 5 responses) + 1,155 biological-module fits, same estimator
(100 µm window, 40-point λ grid, MIN_RECEIVERS = 2000, 400-replicate spatial block bootstrap),
11 sections × 2 sender calls. `clustered_mean` is a section-clustered mean with a t-CI.

### Per control family, `design = base` (naive) — **the frozen 09:06 values** [V]

| response | n features | amplitude (β/sd_y) | 95 % CI | p | flat? |
|---|---|---|---|---|---|
| `all_controls` (**pooled**) | 670 | **−0.0744** | [−0.1306, −0.0182] | **0.0145** | **no** |
| `neg_control_codeword` | 609 | **−0.0604** | [−0.1085, −0.0123] | 0.0188 | no |
| `genomic_control` | 21 | **−0.0307** | [−0.0558, −0.0056] | 0.0213 | no |
| **`neg_control_probe` (the 40)** | 40 | **−0.0225** | [−0.0527, **+0.0078**] | **0.129** | **YES — flat** |
| `neg_probe_rate` (ratio) | — | **+0.0113** | [−0.0085, +0.0310] | 0.232 | yes |
| BIOLOGICAL MODULES (reference) | — | **+0.2767** | [+0.2032, +0.3502] | 8 × 10⁻⁶ | — |

### `all_controls` across designs [V]

| design | amplitude | 95 % CI | p |
|---|---|---|---|
| base (naive) | −0.0744 | [−0.1306, −0.0182] | 0.0145 |
| +N6 neighbour baseline | −0.0625 | [−0.1085, −0.0164] | 0.0128 |
| **+N5 technical covariates** | **+0.0038** | [−0.0186, +0.0261] | **0.715** |
| +N6+N5 (full nuisance) | +0.0053 | [−0.0162, +0.0268] | 0.595 |
| **N2 matched-decoy contrast** | **−0.0642** | [−0.1113, −0.0172] | **0.0124** |

Biological modules for comparison: base **+0.2767**, N6 +0.1124, N5 +0.0694,
**N6+N5 +0.0310**, N2 +0.2662.

**Pre-C6 values that are still in circulation and must NOT be used** (`results/phase3_pre_c6/a7_summary.csv`, 05:19) [V]:
all_controls base −0.0697 [−0.1276, −0.0118] p = 0.023; N2 −0.0611 p = 0.0204; N6+N5 +0.0068
p = 0.411; probes −0.0177 [−0.0453, +0.0099] p = 0.183; codewords −0.0549 p = 0.039; genomic
−0.0337 p = 0.0039; biological base +0.2914, N6+N5 +0.0356.
**`PREREG_PHASE8.md` §10.1 paragraph 1 and `PHASE8_ROADMAP_STATUS.md` 8.5b both quote this
pre-C6 set.**

### Sparsity — why power is the binding constraint (`a7_verdict.txt`) [F]
`neg_control_probe` **0.0067 counts/cell, 0.65 %** of cells non-zero;
`neg_control_codeword` 0.0428 / 3.90 %; `genomic_control` 0.0094 / 0.90 %;
`all_controls` 0.0588 / 5.26 %. Codewords carry ~73 % of the pooled counts.

### The measured false-positive rate

**Definition:** `frac_CI_excludes_zero` — the fraction of the 165 per-fit two-sided 95 % CIs
that exclude zero, on a response with no biology, under the **full N6+N5 design**. [V]

| response | FPR under N6+N5 | under N5 |
|---|---|---|
| `neg_control_codeword` | **0.091** | 0.121 |
| `all_controls` | **0.103** | 0.121 |
| `genomic_control` | **0.109** | 0.109 |
| `neg_control_probe` | **0.145** | 0.139 |
| `neg_probe_rate` | **0.164** | 0.158 |

**Quote as: 9–16 % against a 5 % nominal; 9–15 % on the four count-based responses.**
This matches `SASP_Kernel_Master_Plan.md` §30 5.5. The pre-C6 list was
0.091 / 0.103 / 0.109 / **0.127** / 0.164 — `PHASE8_ROADMAP_STATUS.md`'s "9–13 %" and its
"clean-null subset 0.091 / 0.103 / 0.109 / 0.127" are the **pre-C6** figures.

**Three caveats that must travel with it** (§30 5.5, and audit M3):
1. **Which response.** The five responses are **four overlapping views of one quantity** —
   `all_controls` is the sum of probe + codeword + genomic, and `neg_probe_rate` is a ratio of
   two of them. This is a range over correlated statistics, **not five replications**.
2. **Powered only pooled.** A single fit resolves **±0.1346 SD naive / ±0.1336 SD
   conditioned** (median CI half-width) [F]; pooled CI half-width ±0.018 SD. The conditioned
   biological median \|β\|/sd is **0.0795**, so one fit resolves **1.7×** that.
3. **It is the estimator's rate, not the reportable-fit filter's.** The reportable filter
   (`beta_naive > 0 AND beta_base_lo > 0`, one-sided, naive design) admits **3.0–13.3 %**
   across the five control responses; on the full N6+N5 design it admits **4.8 %, identical
   across all five — essentially nominal.** [V]

```bash
python3 - <<'PY'
import pandas as pd
f=pd.read_csv('results/phase3/a7_control_probe_fits.csv')
for m,g in f.groupby('response'):
    print(m, round(((g.beta_naive>0)&(g.beta_base_lo>0)).mean(),4),
             round(((g.beta_n6n5>0)&(g.beta_n6n5_lo>0)).mean(),4))
PY
```

### Moran's I of our own controls — **now computed** (§30 5.5 says "not yet computed"; it is)

**Source: `results/moran/`.** 11 M1 sections, primary weights **k = 6 NN row-standardised**,
999 conditional permutations, section-clustered mean with t-CI.
Producers: `code/run_moran_controls.py --sections all --n-jobs 3 --perms 999`,
`code/run_moran_lognorm.py`, `code/moran_kernel_power.py`, `code/summarize_moran.py`,
`code/make_moran_figure.py`. Validation against `esda`/`libpysal`: agreement **0.00e+00** on
I and on `z_rand`. [F]

Aggregated per control family (`results/moran/moran_pooled.csv`) [V]:

| field | Moran's I raw | 95 % CI | clustered p | max per-section p_rand |
|---|---|---|---|---|
| `all_controls` | **+0.0455** | [+0.0302, +0.0609] | 6.0e-05 | **9.0e-13** |
| `neg_control_codeword` | +0.0421 | [+0.0281, +0.0561] | 5.4e-05 | 5.2e-12 |
| `neg_control_probe` | **+0.0058** | [+0.0033, +0.0084] | 4.2e-04 | **0.537 (n.s. per section)** |
| `genomic_control` | +0.0042 | [+0.0022, +0.0062] | 9.1e-04 | **0.802 (n.s.)** |
| `neg_probe_rate` | +0.0047 | [−0.0010, +0.0105] | 0.095 | **0.965 (n.s.)** |
| Tier B modules | +0.085 to +0.244 | — | — | — |
| `density_50um` (technical) | +0.9483 | [+0.9362, +0.9605] | — | — |

**The sparsity control — the genuinely new point.** Voyager's controls-vs-genes contrast is an
**abundance** contrast. Per-feature medians (`moran_verdict.txt`) [F]:

| control class | n | med counts | med I | **sparsity-matched genes** | med counts | med I |
|---|---|---|---|---|---|---|
| Negative Control Probe | 440 | 21 | −0.00012 | 4,495 | 32 | **−0.00018** |
| Negative Control Codeword | 6,699 | 7 | −0.00004 | 1,179 | 18 | **−0.00012** |
| Genomic Control | 231 | 63 | −0.00029 | 8,735 | 62 | **−0.00025** |

A median Gene Expression feature carries **5,885 counts per section against a probe's 21 —
280×.** The per-feature statistic has no power at control abundance.

### **The Moran power bound** (`results/moran/moran_kernel_power.csv`, 22 rows = 11 sections × 2 calls)

Model: ΔI = β_z² · Var(k) · I(k), k_i = exp(−d_i/λ̂). [V]

| quantity | median | range |
|---|---|---|
| ΔI contributed by the **whole** A7 gradient | **2.20 × 10⁻⁴** | 1.4e-06 – 2.6e-03 |
| ΔI as a fraction of the observed control I | **0.83 %** | 0.002 % – 6.1 % |
| SE(I) on the same cell set | 1.50 × 10⁻³ | 0.0012 – 0.0023 |
| **smallest β Moran's I could resolve (ΔI = 2·SE)** | **0.362 SD** | **0.308 – 1.070** |
| — restricted to `tierA_p95` | 0.338 SD | 0.308 – 0.418 |

**0.362 SD exceeds the project's own naive biological amplitude (0.277) and is 5× the A7
control gradient (0.074).** Moran's I could not have detected the paper's headline effect.

**Rank agreement (`moran_verdict.txt`, re-derived) [V]:** over the 12 control + module fields,
\|Moran I\| and \|A7 naive amplitude\| rank together — **ρ = +0.895 (p = 8.37e-05) raw,
+0.944 (p = 3.93e-06) cell-type-centred**. Within the 5 control responses × 11 sections
(55 pairs) they are **uncorrelated, ρ = +0.155, p = 0.259** — the honest form of "different
questions". Discordance test: **14 of 55** section × response cells are near-zero by Moran and
carry an A7 gradient.

**Two disagreements, report them against interest:** (a) Moran's I does **not** call the 40
probes flat — its pooled CI excludes zero (+0.0058 [+0.0033, +0.0084]); A7 calls them flat
only because A7's CI is wide. (b) `genomic_control` is A7's third-largest amplitude but
Moran's **smallest** I; the two statistics disagree on that response's rank.

### ⚠ The "naive biological amplitude" has four values in circulation
**Decide once, name the estimator.** [V]

> ## ✅ RESOLVED 2026-08-27 — **+0.2767**, the section-clustered signed mean, `design = base`
>
> **Authoritative: `naive biological amplitude = +0.2767` response-SD**, defined as the
> **section-clustered signed mean of β̂/sd(y) over the 1,155 biological-module fits under
> `design = base`** — `results/phase3/a7_summary.csv`, row `BIOLOGICAL MODULES (reference)`,
> column `clustered_mean`, frozen 09:06. Its conditioned counterpart is **+0.0310** (`n6n5`).
>
> **Why the clustered mean and not the median |β|/sd:** the clustered mean is A7's own
> primary statistic — it is the only one of the two that carries an interval and a p-value
> (`clustered_lo/hi/p`, a section-clustered t-CI), it is the statistic every A7 verdict in the
> paper is decided on (flat / not flat, per family and per design), and it is signed, so it
> is comparable with the control gradients it is contrasted against. The median |β|/sd is
> unsigned and interval-free.
>
> **The companion, when you want it: 0.3120** = median |β|/sd over the *same* 1,155 fits
> (`median_abs_amplitude`), conditioned counterpart **0.0795**. Legitimate, different
> estimand — **name it explicitly whenever you use it.**
>
> **0.2914 and 0.314 are the pre-C6 vintages of those two estimators. Do not use them.**
> The power argument is unaffected: 0.362 SD exceeds all four.
>
> ```bash
> python3 -c "
> import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv')
> print(d[d.response.str.startswith('BIO')][['design','clustered_mean','median_abs_amplitude']])"
> ```

| value | what it is | file |
|---|---|---|
| **0.2767** | section-clustered signed mean, `design=base`, **FROZEN** | `results/phase3/a7_summary.csv` |
| 0.2914 | the same, **pre-C6** | `results/phase3_pre_c6/a7_summary.csv` |
| 0.3120 | median **\|β\|/sd** over the 1,155 module fits, frozen | `results/phase3/a7_verdict.txt` |
| 0.314 | the same, pre-C6 | `PREREG_PHASE8.md` P2 |

Conditioned counterparts: **0.0310** (clustered mean, frozen), 0.0356 (pre-C6),
**0.0795** (median \|β\|/sd, frozen), 0.077 (pre-C6). The power argument is unaffected — 0.362
exceeds all four — but the number must be consistent within the paper.

---

## §5.6 — A matched-decoy contrast does not remove what covariate adjustment removes

**The strongest contribution. Three independent lines** (`CORRECTIONS.md` §13 / §C.1):

| line | matched decoys | covariates | ground truth? |
|---|---|---|---|
| **Figure 1 — synthetic, planted kernel** | λ̂ CI coverage **0.346** | **0.854** | **YES** |
| **A7 — the platform's own technical gradient** | N2 leaves it: **−0.0642, p = 0.0124** | N5 removes it: **+0.0038, p = 0.715** | no |
| **Composition-matched protocol** | removes **1.6 %** | removes **85.4 %** | no |

*(naive coverage in line 1 is 0.5125, so the matched decoy is worse than doing nothing.)*
**Lead with Figure 1.** The other two are answerable with "your covariates removed real
signal"; the planted one is not. `NOVELTY_ASSESSMENT.md` knows only line 2 and ranks it #1 of
six anyway; `SUBMISSION_PATCH` §0.4 knows two. **Three is the current count.**

### The composition-matched contrast, in full

**Source: `results/phase3/compmatch_reruns.csv`, `row_type == 'summary'`.** Fit-level audit
trail: `compmatch_fits.csv` (6,237 rows) = merge of `_tierA` + `_tierApm`.
Producer `code/run_phase8_compmatch.py`; driver `code/_compmatch_chain.sh`.

**All four rows below are the same 42 fits / 33 reportable, pooled scope, same variables.** [V]

| control | SF | 95 % CI | **composition share removed** |
|---|---|---|---|
| `comp` — 20-NN composition as **matched decoys** (the protocol) | **0.9837** | [0.973, 0.994] | **1.6 %** |
| `comp_adj` — the same variables as **covariates** | 0.4989 | [0.421, 0.606] | 50.1 % |
| `type_adj` — receiver's own cell-type intercepts | **0.3414** | [0.236, 0.402] | **65.9 %** |
| `typecomp_adj` — cell type **+** 20-NN composition | **0.1461** | [0.052, 0.246] | **85.4 %** |

**85.4 / 1.6 ≈ 53 — "a factor of fifty".**
Within receiver cell type (748 reportable fits) the matched protocol removes **3.6 %**
(SF 0.9647); the published N2 matching set (`full`) removes 1.4 % pooled / 5.1 % by cell type.

Per-module sensitivity sender sets (`tierApm_p95`, `A_sender_for_<module>.txt`, n = 34):
`comp` **0.9861** (1.4 %), `comp_adj` 0.5140, `type_adj` **0.3350** (66.5 %),
`typecomp_adj` **0.2042** (79.6 %) [V]. **Every paired row differs by < 0.07 in SF.**

**Five seeds, frozen: 20260901–20260905.** Across-seed spread of the protocol's own number:
pooled **0.98370–0.98397, sd 1.2 × 10⁻⁴**; by cell type 0.9629–0.9648, sd 8.3 × 10⁻⁴ [V].
The covariate-adjusted variants are seed-free. The only thing the seed moves visibly is which
fits are reportable (149 vs 150 of 315).

**Matching balance:** max \|SMD\| **0.092 → 0.035**; median match rate **0.99987**; the §8
Test-5 gate (\|SMD\| ≤ 0.1) passes in **100 %** of matches; worst case over all 6,237 fits
0.0759. Caliper 0.25 SD; exact stratification on receiver cell type; 1-1 NN propensity match
without replacement.

```bash
python3 -c "
import pandas as pd
d=pd.read_csv('results/phase3/compmatch_reruns.csv')
s=d[(d.row_type=='summary')&(d.scope_kind=='pooled')&(d.call=='tierA_p95')]
print(s[['variant','n_reportable','median_sf_matched','median_sf_matched_lo','median_sf_matched_hi','median_comp_share']].to_string(index=False))"
```

**Producers.**
```bash
python3 -u code/run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierA_p95 \
    --variants comp,full,comp_adj,type_adj,typecomp_adj --out-tag _tierA
python3 -u code/run_phase8_compmatch.py --arm m1 --n-jobs 6 --calls tierApm_p95 \
    --variants comp,comp_adj,type_adj,typecomp_adj --out-tag _tierApm
python3 -u code/run_phase8_compmatch.py --merge \
    results/phase3/compmatch_fits_tierA.csv,results/phase3/compmatch_fits_tierApm.csv
```
`--arm h1` **refuses to run** until Phase 8's tag is cut and `SASP_H1_UNFROZEN=1` is set.

### The named mechanism, in one sentence
Per-cell **detection efficiency** is a property of the cell, not of its neighbourhood.
Propensity matching on neighbourhood covariates cannot see it; entering the same variables as
covariates can. That is why N2 leaves the A7 gradient 86 % undiminished while N5 removes it
completely, and why the same 20-NN composition vector removes 1.6 % as a decoy set and 85.4 %
as a regressor.

---

## §5.7 — Random-shift nulls on a non-convex window

**Frame as import and quantification, never as discovery.** See §7.1.

### The full variant table
**Source: `results/phase3/sf_summary_var.csv` / `results/phase3/summary_phase3_var.txt`.**
n = **153** reportable fits on every whole-section row, **136** on the two tile rows (17 fall
below the 2,000-receiver floor once the fit is restricted to tiles). 1,000 permutations
throughout. **Brackets are IQRs.** [F, table read directly; the two headline rows spot-checked V]

| variant | median SF | IQR | keeps a real neighbour ≤100 µm | median displacement | full N5+N6+zon design |
|---|---|---|---|---|---|
| N3 torus, whole-section bounding box (**published**) | **0.999** | [0.989, 1.006] | 0.772 | 2,910 µm | 1.001 |
| N3, same, re-run in the C1 job | 1.001 | [0.991, 1.008] | 0.772 | 2,910 µm | 0.999 |
| N3-tile, torus inside solid-tissue tiles | **0.971** | [0.906, 1.009] | 1.000 | 479 µm | 0.972 |
| N3-occ, ≤5 % of senders out of tissue | **0.302** | [0.000, 0.734] | 1.000 | **28 µm** | 0.287 |
| N3-occ15, ≤15 % out (suppl.) | 0.940 | [0.761, 1.011] | 0.969 | 317 µm | 0.917 |
| N3-swap, senders to random real cell positions | **0.695** | [0.392, 0.872] | 1.000 | 3,241 µm | **1.003** |
| N3-snap, shift then snap to nearest cell (suppl.) | 0.993 | [0.950, 1.017] | 1.000 | 2,977 µm | 0.976 |
| **N3-var — variance-corrected Euclidean shift (PRIMARY)** | **0.996** | **[0.975, 1.007]** | **1.000** | **2,215 µm** | **0.997** |
| N4 rotation, bounding box (**published**) | **0.947** | [0.804, 1.039] | 0.920 | 3,194 µm | 0.992 |
| N4, same, re-run in the C1 job | 0.952 | [0.796, 1.048] | 0.920 | 3,194 µm | 1.001 |
| N4-tile | **0.924** | [0.835, 1.049] | 1.000 | 589 µm | 0.994 |
| N4-occ | **0.183** | [0.000, 0.559] | 1.000 | **25 µm** | 0.198 |
| N4-occ15 (suppl.) | 0.883 | [0.690, 0.989] | 0.969 | 320 µm | 0.893 |
| N4-swap | 0.946 | [0.786, 1.035] | 1.000 | 2,980 µm | 0.958 |
| **N4-var (PRIMARY)** | **0.985** | **[0.958, 1.003]** | **1.000** | **3,395 µm** | **0.999** |

Window-matched N3-var **0.995 [0.975, 1.008]**; N4-var 0.985. The `full_sf` values for the two
`var` rows come from `perm_nulls_var_full200.csv` at **200** permutations, not 1,000.
**Do not read the `rej` column across rows** — for the `var` rows it is the RS_count Monte
Carlo test of Mrkvička et al., for every other row the repo's uncorrected permutation p.

**N3-var, not N3-tile, is the primary corrected N3** (`CS_PHASE8_TORUS_VAR.md` §2,
`CORRECTIONS.md` §8.1). N3-tile is now a supporting variant.

**Superseded values still live in documents:** `CS_PHASE7_C1.md` §0 is keyed to **160**
reportable fits (pre-C6): N3-tile 0.974 (→ **0.971**), published N3 1.000 (→ **0.999**),
N3-occ 0.349 (→ **0.302**), N3-swap 0.721 (→ **0.695**), N4-tile 0.962 (→ **0.924**), N4-occ
0.273 (→ **0.183**), N1 0.716 (→ **0.707**).

### N3-swap ≡ N1 — flag it
N3-swap 0.695 against N1 0.707; per-fit Spearman **ρ = 0.948**, median \|difference\| 0.0087;
under the full covariate design it moves to **1.003**. It is a label permutation, not a torus
shift, and must be flagged as such in Figure 2c. [F]

### The destructiveness diagnostic — **two columns, never merged**

**Source: `results/phase3/null_destructiveness.csv`, medians over the 6 in-band sections.** [V]

| null | **out of tissue** = `1 − frac_in_occupancy` | **lost every real neighbour ≤100 µm** = `1 − frac_retaining_a_neighbour` |
|---|---|---|
| N3 bounding box | **35.5 %** | **22.8 %** |
| N4 bounding box | **19.9 %** | **8.0 %** |

**These are different quantities in the same file. Never state both in one sentence, and never
call 22.8 % / 8.0 % "out of tissue".** (Audit R3; the fix was applied to `CS_PHASE7_C1.md`,
`CS_PHASE8_C1_CLOSEOUT.md` and `PHASE8_ROADMAP_STATUS.md` but **not** to
`CS_PHASE8_TORUS_VAR.md`, whose §1 table and §10 framing blockquote still say "23 % in the
void" / "8 % in the void".)

Neighbour thinning: median real cells within 100 µm **140.0 → 119.7** under published N3
(−14.5 %), → 129.7 under N4, → 139.2 under N3-tile (−2.6 %), → **133.5 under N3-var (−4.6 %)** [F].
N3-var retains **59.2 %** of senders and **60.3 %** of receiver cells in W — the 41 % cost the
variance standardisation exists to absorb.

```bash
python3 -c "
import pandas as pd
d=pd.read_csv('results/phase3/null_destructiveness.csv'); d['s']=d.section.str.split('_').str[0]
d=d[d.s.isin(['7259','7260','7001','7248','7352','7435'])]
print(d.groupby('null')[['frac_in_occupancy','frac_retaining_a_neighbour','median_displacement_um']].median())"
```

### The admissible-offset enumeration — **frozen values are 1–66, 28 µm**

`results/phase3/null_destructiveness.csv`, N3_occ, in-band [V]:
per section **66, 6, 1, 9, 10, 12** admissible offsets of **38,080 – 108,375** candidates
(0.001 %–0.17 %); median displacement **28.3 µm**. **Section 7001 admits only the identity**
(1 of 108,375, displacement 0 µm) and returns SF = −0.000 by construction. N4_occ: **1–13 of
720** angles.

**The pre-C6 file gives 1–63 and 27.4 µm** (`results/phase3_pre_c6/null_destructiveness.csv`) [V].
`CS_PHASE8_TORUS_VAR.md` §10 and `NOVELTY_ASSESSMENT.md` §2.2 quote **1–63 / 27 µm** — stale.
`SUBMISSION_PATCH` §6's **1–66 / 28 µm** is correct.

### The FFT trick — **lead the methods paragraph with it**
"Which translations keep ≥ x % of a point set inside a mask?" is a **circular
cross-correlation**, so a single `rfft2` gives the exact admissible offset set over all offsets
at once. Implementation: `code/phase3_null_geom.py`; described in `CS_PHASE7_C1.md` §7.
`CS_PHASE8_TORUS_VAR.md` mentions it but carries no method description — **write it from
`phase3_null_geom.py`.**

### The direct calibration measurement — **the part no spatial-omics paper has**

**Source: `results/phase3/var_sim_calibration.csv`**, producer `code/phase3_var_sim.py`.
Design: 100 sampling points, two independent stationary Gaussian fields with isotropic
exponential correlation exp(−r/s), 199 shifts on a disk of radius ½, **400 replications**
(Monte Carlo SE ≈ 0.011). Fields synthesised on a 3×-larger domain and cropped, so the
realisation is **not** periodic. Irregular window = three overlapping discs minus a bite,
**74 %** of the bounding box, inside the 0.658–0.858 range the six real sections show. [V]

Rejection rate at nominal 0.05, under the null of independence:

| window | s | whole-window torus | **torus in 4×4 tiles** | 8×8 tiles | **RS_count (variance correction)** | RS_ker | drop, no standardisation |
|---|---|---|---|---|---|---|---|
| rectangle | 0.02 | 0.035 | 0.040 | 0.035 | **0.033** | 0.038 | 0.005 |
| rectangle | 0.05 | 0.055 | 0.065 | 0.080 | **0.055** | 0.050 | 0.013 |
| rectangle | 0.15 | 0.048 | 0.063 | 0.063 | **0.035** | 0.013 | 0.013 |
| rectangle | 0.30 | **0.078** | **0.105** | 0.055 | **0.060** | 0.035 | 0.018 |
| irregular | 0.02 | 0.048 | 0.048 | 0.040 | **0.040** | 0.038 | 0.003 |
| irregular | 0.05 | 0.033 | **0.080** | 0.073 | **0.043** | 0.033 | 0.003 |
| irregular | 0.15 | 0.040 | **0.083** | 0.050 | **0.053** | 0.033 | 0.008 |
| irregular | 0.30 | 0.073 | **0.118** | 0.085 | **0.055** | 0.040 | 0.020 |

- **"Tiled torus rejects at 0.080–0.118, up to 2.4× nominal"** is the **irregular-window 4×4
  rows at s ≥ 0.05**. The s = 0.02 irregular row is **0.048**. **Write it as
  "0.080–0.118 for s ≥ 0.05"** or "0.048–0.118 across all four correlation scales". 0.118/0.05
  = 2.36; `CORRECTIONS.md` §C.4 rounds to 2.35×, `SUBMISSION_PATCH` to 2.4×. [V]
- **RS_count holds 0.033–0.060 across every window and every correlation scale** — never
  outside ±0.010 of nominal. **`CS_PHASE8_M1_RERUN.md` §6 and §14.3 say "0.040–0.060"; the
  true minimum is 0.033.** [V]
- Whole-window torus: **0.033–0.078 overall**; 0.033–0.073 on the irregular window;
  **0.073–0.078 at strong autocorrelation**. [V]
- Tiling is more liberal than the whole-window torus in **7 of the 8 cells** (1 tie). [V]

**Two caveats that must travel with it:** (i) the study is **synthetic** — Gaussian fields,
100 sampling points — and establishes the *direction* of the tiling effect, not a type-I error
number for the Phase 3 fits; the instrument for that is A7 (§5.5). (ii) On the real data the
tiled null looks slightly **more conservative** (SF 0.971 vs 0.999; raw rejection 0.801 vs
0.824), because at a 1,200 µm tile side the seams are far apart relative to λ̂.

**Reported against interest: our own C1 correction replaced a liberal test with a more liberal
one.** That is exactly why N3-var and not N3-tile is primary.

### Variance-correction validation (`results/phase3/var_variance_check.csv`) [F]
RS_count predicts slope −0.5 in log sd(T \| n-bin) on log n. Measured: **N3-var median −0.451**
(range −0.585 … −0.122), **N4-var −0.492** (−0.644 … −0.119), over the 12 of 16 cases with
≥2× dynamic range in n_i. Estimator cross-check: β̂_obs and λ̂ in `perm_nulls_var.csv` agree
with `perm_nulls_c1.csv` to max \|Δβ\| = 1.04 × 10⁻¹⁶ and max \|Δλ\| = 0 over all 315 rows.

**Producers.**
```bash
python3 -u code/run_phase3_var.py --stage diag --n-rep 20 --n-jobs 6
python3 -u code/run_phase3_var.py --stage perm --n-perm 1000 --n-jobs 6 --no-full
python3 -u code/run_phase3_var.py --stage perm --n-perm 200  --n-jobs 6 --tag _full200
python3 -u code/phase3_var_validate.py
python3 -c "import sys;sys.path.insert(0,'code');import phase3_var_sim as S;S.main(n_rep=400,n_shift=199,n_jobs=8)"
python3 -u code/summarize_phase3_var.py
```

---

## §5.8 — Two-arm replication (Figures 5 and 6) — **NOT RUN. See §8.**

Only the mouse arm exists. Every number the section needs from the human arm is unmeasured.
What **is** available now:

- The geometric predictions the section would compare across arms: Poisson identity
  **r² = 0.984, slope −0.525**; grid-railing rate **60 %**; Ripley ratio **1.124** — mouse only [F].
- H1 descriptors (§2.0) — acquisition and structural verification only.
- The pre-registered replication criterion and §18 outcome table: `PREREG_PHASE8.md` §6, §7.
- The §8 DeepScence native-vs-remapped prediction, stated before H1 runs: `PREREG_PHASE8.md` §8.

**If the human arm has not run by submission, say so and report the mouse arm alone. Do not
imply a replication that has not happened.**

---

## §5.9 — What we withdrew, and why that is a result

### The six labelled bases — **the authoritative table**

**Source: `results/phase3/caller_coverage_gate_headline.csv` (6 rows, `n_pairs = 4`).**
Producer: `code/caller_disagree_all.py --all` → `code/summarize_caller_coverage.py`.
Pooling: Mantel–Haenszel over sections, stratum-exact, after conditioning on cell type and
within-cell-type sequencing depth. **All six rows verified [V].**

| basis (literal string in the CSV) | Tier A | n_sec | band | median | **pooled** | z | p | above chance | sign-test p |
|---|---|---|---|---|---|---|---|---|---|
| `2-section, pre-C6 Tier A (PUBLISHED)` | 25 | 2 | 0.932–1.369 | 1.010 | **1.040** | 1.76 | 0.078 | 4/8 | 1.000 |
| `11-section, pre-C6 Tier A (task 8.4)` | 25 | 11 | 0.700–1.711 | 1.156 | **1.129** | 13.35 | 1.12e-40 | 29/44 | 0.0488 |
| `2-section, post-C6 Tier A (FROZEN)` | 33 | 2 | 0.979–1.442 | 1.099 | **1.131** | 5.80 | 6.5e-09 | 6/8 | 0.289 |
| **`11-section, post-C6 Tier A (FROZEN)` — USE THIS** | 33 | 11 | **0.751–2.198** | 1.190 | **1.212** | 21.92 | **1.84e-106** | 35/44 | 1.06e-04 |
| `6-section, in-band only, pre-C6 Tier A (task 8.4)` | 25 | 6 | 0.775–1.374 | 1.131 | 1.115 | 8.99 | 2.56e-19 | 16/24 | 0.152 |
| `6-section, in-band only, post-C6 Tier A (FROZEN)` | 33 | 6 | 0.811–1.565 | 1.168 | 1.167 | 13.09 | 3.63e-39 | 18/24 | 0.0227 |

The **four pooled pairs** are: Tier A × SenePy, Tier A × DeepScence, Tier A × `Cdkn1a`⁺,
SenePy × `Cdkn1a`⁺. **Excluded from every pooled number:** DeepScence × `Cdkn1a`⁺ (circular)
and SenePy × DeepScence. [V]

### The 3-pair basis (the one that literally defines the published "0.93–1.22×" band)
Emitted by **no file**; re-derived here by MH pooling of the three Tier-A pairs [V]:

| basis | pooled | z | p | band |
|---|---|---|---|---|
| 2-section pre-C6 (published) | **1.0299** | 1.27 | 0.203 | 0.932–1.221 |
| 11-section pre-C6 | **1.1182** | 11.49 | 1.44e-30 | 0.700–1.711 |
| **2-section post-C6 (frozen)** | **1.1283** | 5.47 | **4.40e-08** | 0.982–1.442 |
| **11-section post-C6 (frozen)** | **1.2122** | 20.62 | **1.84e-94** | 0.751–2.198 |

```bash
python3 - <<'PY'
import pandas as pd, numpy as np
from math import erfc, sqrt
P3=['tierA_score|senepy_score','tierA_score|deepscence_score','tierA_score|cdkn1a_counts']
P4=P3+['senepy_score|cdkn1a_counts']
def pool(d,pairs):
    d=d[(d.A+'|'+d.B).isin(pairs)]
    obs,exp=d.n_both.sum(),d.exp_both_stratified.sum()
    z=(obs-exp)/np.sqrt((d.sd_both_stratified**2).sum())
    return obs/exp, z, erfc(abs(z)/sqrt(2)), d.ratio_stratified.min(), d.ratio_stratified.max()
for f in ['results/phase3/caller_agreement_matched_significance_11sections.csv',
          'results/phase3/caller_agreement_matched_significance_2sec_c6.csv',
          'results/phase3_pre_c6/caller_agreement_matched_significance_11sections.csv',
          'results/phase3_pre_c6/caller_agreement_matched_significance_verify2sec.csv']:
    d=pd.read_csv(f); print(f.split('/')[-1], "3p", pool(d,P3), "4p", pool(d,P4))
PY
```

### ⚠ Do not mix bases inside one sentence
`SASP_Kernel_Master_Plan.md` §30 5.9 pairs **"1.212×, p = 1.8 × 10⁻¹⁰⁶"** (a **4-pair** value)
with **"1.128×, p = 4.4 × 10⁻⁸"** (a **3-pair** value). The consistent pairs are:
**3-pair: 1.128 (p = 4.4e-8) → 1.212 (p = 1.8e-94)**; **4-pair: 1.131 (p = 6.5e-9) → 1.212
(p = 1.8e-106)**. Pick one pair-count per sentence.
`PHASE8_ROADMAP_STATUS.md`'s 8.4 table labels all three of its rows "3-pair band" but its
6-section row is a **4-pair** number lifted from the headline CSV.

### The decomposition to tell — **the sender-set fix kills independence at n = 2**

| step | what changes | 4-pair pooled | p |
|---|---|---|---|
| published state | 2 sections, pre-C6 25-gene Tier A | 1.040 | 0.078 |
| ① fix the contaminated sender/response split (coverage held at 2) | 2 sections, frozen 33-gene | **1.131** | **6.5 × 10⁻⁹** |
| ② add coverage (Tier A held frozen) | 11 sections | **1.212** | **1.8 × 10⁻¹⁰⁶** |

**Say it in that order: the published sentence was wrong on the published data.** "More
sections revealed dependence" is the wrong story.

### Per pair, **frozen Tier A, 11 sections** (`results/phase3/caller_coverage_gate.csv`) [V]

| pair | pooled | z | p | min–max | above chance | sig. above | sig. below |
|---|---|---|---|---|---|---|---|
| Tier A vs `Cdkn1a`⁺ | **1.471** | 19.45 | 2.7e-84 | 1.085–2.198 | **11/11** | 9 | 0 |
| Tier A vs DeepScence | **1.288** | 19.23 | 2.1e-82 | 1.096–1.660 | **11/11** | 10 | 0 |
| SenePy vs `Cdkn1a`⁺ | 1.211 | 7.43 | 1.1e-13 | 0.842–1.391 | 9/11 | 7 | 0 |
| **Tier A vs SenePy** | **0.972** | −1.63 | **0.104** | 0.751–1.179 | **4/11** | 1 | 3 |
| **SenePy vs DeepScence** | **0.737** | −15.08 | 2.3e-51 | 0.332–2.150 | 1/11 | 1 | **10** |
| *(circular, never pooled)* DeepScence vs `Cdkn1a`⁺ | 1.255 | 10.53 | 6.2e-26 | **0.963–2.849**, **median 1.071** | 7/11 | 5 | 0 |

On the published two-section base, frozen sets: Tier A × `Cdkn1a`⁺ 1.017 → **1.300**;
Tier A × DeepScence 1.103 → **1.179**; Tier A × SenePy 0.935 → **1.007**. [V]

**Pre-C6 per-pair values still quoted in `PHASE8_ROADMAP_STATUS.md` — stale:**
Tier A × DeepScence 1.248; Tier A × SenePy **0.914, below chance in 11 of 11**. The frozen
value is **0.972, p = 0.104, above chance in 4 of 11** — the "one pair genuinely below chance"
plank is **dead**.

### SenePy vs DeepScence — the replacement plank
Byte-identical pre- and post-C6 (neither caller involves Tier A). Published reading was
"2.15× in sham, 0.38× in SBR = concordant in one arm, anti-concordant in the other". At full
coverage: **0.33–0.55× in ten of eleven sections**, single exception 2.150 in **7250**.
**There is no arm effect. There is one anomalous section.** Re-anchored on the caller-free
8-gene proliferation set it strengthens to **0.495 (z = −28.9)**. **Mandatory caveat:** ranking
by \|score\| puts the same pair **at chance, 1.025, z = 1.4, n.s.** — the anti-concordance is
about polarity, not about which cells are extreme. [F]

### Depth camps — the identifiable cause at each end
`results/phase3/caller_within_type_depth_bias_11sections.csv`, within-cell-type Q5/Q1
transcript-depth enrichment, all 11 sections [F]:

| caller | Q5/Q1 | direction | sections |
|---|---|---|---|
| SenePy | **10.58 – 41.74** | top-selecting | 11/11 |
| `Cdkn1a`⁺ | **4.19 – 42.36** | top-selecting | 11/11 |
| Tier A arrest score | **0.146 – 0.317** | bottom-selecting | 11/11 |
| DeepScence | **0.218 – 0.795** | bottom-selecting | 11/11 |

**But the rule "the direction of each pair is predicted by its depth loading" is REFUTED**:
pair-level exact permutation **p = 0.30**; within-pair Spearman ρ **−0.16 to −0.70 in all five
pairs** (Tier A × `Cdkn1a`⁺ ρ = −0.700, p = 0.016 — the wrong direction); pooled continuous
ρ = +0.096, p = 0.49. The Mann–Whitney same-camp-vs-opposite result (1.258 vs 0.990,
p = 0.008) is **pseudo-replicated — do not quote it.** [F]

### DeepScence anchor instability (`results/phase3/deepscence_anchor_decisions.csv`) [F]
Depth-partialled ρ with the published `CDKN1A` anchor is **negative in 7248 (−0.012) and 7435
(−0.024)** and effectively zero in 7352 (+0.0021); fold-split sign stability drops to **0.60**
in 7352 and 0.85 in 7248, against 0.95–1.00 for the proliferation anchor. Detection rates:
`Cdkn1a` in **0.68–21.6 %** of cells; the 8-gene proliferation set in 3.5–12.8 %; `Lmnb1` in
2.2–8.9 %. Re-anchoring changes **1 sign bit per section** and leaves D1 untouched.
**Cite DeepScence's `CDKN1A` sign-anchoring as documented method behaviour**, so this reads as
a comparability caveat rather than a bug we found.

### `denoise=False` — the frozen primary, and its measured cost (`results/phase8_d2/`)
Producers `code/setup_dca_env.sh`, `code/run_deepscence_dca.py`,
`code/run_deepscence_denoise_probe.py`, `code/analyse_d2_denoise.py`,
`code/analyse_d2_stability.py`, `code/report_d2_tables.py`. DCA 0.3.4 / TF 2.4.4 in an
isolated py3.8 venv; main environment untouched. [F]

| | 7259 | 7352 | 7239 |
|---|---|---|---|
| depth loading ρ(score, counts), `denoise=False` → `True` | **0.318 → 0.531** | **0.410 → 0.542** | **0.389 → 0.640** |
| hepatocyte % of top-5 % calls | 64.5 → **100.0** | 97.2 → **100.0** | 71.5 → **100.0** |
| Pearson r between the two configurations | 0.614 | 0.671 | 0.718 |
| global top-5 % Jaccard | 0.141 | 0.118 | 0.280 |
| % of the committed sender set not called | 75.3 % | 78.8 % | 56.2 % |

**Seed instability at the published default** (`d2_stability.csv`, one fixed 20,000-cell
subsample of 7239): seeds 0 and 2 agree at Jaccard **0.665**; **seed 1 shares no cells with
either — Jaccard 0.000, 2,000 of 2,000 cells changing status, twice.** `denoise=False` agrees
seed-to-seed at **Jaccard 0.76–0.99**. DeepScence hands `random_state` straight to `dca()`, so
one seed moves both stages.

**Normalisation (`d2_normalisation_strength.csv`):** §6's named estimator `mor` removes only
**11.4–23.6 %** of log-depth variance on this panel (4 sections); `lib` removes **100 %**
(p90/p10 depth ratio 11.0 → 1.000), moves the depth loading by **−74 % (7259)** and **−93 %
(7352)**, and changes **46 %** and **100 %** of sender calls. Determinism floor: `raw` vs
committed at the same seed, Pearson r = 0.9999991–0.9999999, 2–24 cells of 75k–115k changed. [F]

---

## §6 — Corrections, pre-registration, reproducibility

| item | value | source |
|---|---|---|
| pre-registration | `reports/PREREG_PHASE8.md`, 993 lines; complete apart from tag hashes | — |
| frozen tag | `phase8-frozen` (and `pre-c6-genesets` @ `e789372`, 05:41:13) | — |
| corrections ledger | `reports/CORRECTIONS.md`, §B = 29 rows of pre / post / attributed cause | — |
| **the frame sentence** | "**Every correction made today moved against interest. None of them changed the conclusion.**" | `CORRECTIONS.md` §A L22–24 |
| gene-set freeze manifest | `genesets/human/FROZEN_MANIFEST.csv`, **43 rows, 35 FROZEN / 8 variants**, SHA-256 per file, 0 mismatches | audit C1 |
| deviation tables | `PREREG_PHASE8.md` §13 (P1–P28) and `PREREG_PHASE8_genesets.md` §12 | — |
| figure guard | **`git ls-files figures/` = 52; 52 artefacts on disk** — the Phase 8 figures have been committed, so the guard now covers everything | verified here [V] |
| **stale guard claim** | "all 27 committed figures match" / "27 of 45" — **do not quote as current** | audit M7 |

**Deviation-ID namespace collision — binding on the writer.** `D1`–`D17` means **three
different things**: PI decisions (`PHASE8_ROADMAP_STATUS.md`), gene-set deviations
(`PREREG_PHASE8_genesets.md` §12), and — separately, without collision — `P1`–`P28`
(`PREREG_PHASE8.md`). **Always name the series; never cite a bare `D<n>`.**

**Environment (Methods).** numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.0,
anndata 0.12.19, scanpy 1.11.5, joblib 1.5.3, h5py 3.16.0, matplotlib 3.11.1, libpysal 4.14.1,
esda 2.8.2, on CPython 3.11.10; DCA 0.3.4 / TF 2.4.4 / numpy 1.19.5 in an isolated CPython
3.8.19 venv. Memory ceiling is a **57.7 GiB cgroup**. [F]

---

## §7 (Discussion) and §8 (Appendix) — numbers already listed above

Discussion needs: the bound and its distance range (§5.3); the three N2-vs-N5 lines (§5.6); the
torus attribution (§5.7); the caller decomposition (§5.9); the `denoise` limitation (§5.9);
**and the statement that no published length constant with uncertainty exists for this
quantity, so the bound contradicts no prior number.**

Appendix needs: the full null battery (§5.2), the eleven-variant shift-null family (§5.7), the
five-seed composition protocol (§5.6), the per-section caller tables
(`caller_coverage_gate.csv`), the intersection matrix
(`figures/figure_gs1_intersection_matrix_data.csv`, **196 cells, every Tier A × Tier B cell
zero, but 18 of 21 Tier B pairs non-disjoint**), and the deviation tables.

---

# 3. THE FORBIDDEN-CLAIMS CHECKLIST

Run this against any draft. Items 1–12 are `PREREG_PHASE8.md` §10 verbatim in substance;
13–20 are the falsifications and the corrections audits. Every item is a consequence of a
measurement already made, and they bind **both** arms.

| # | ✗ Never write | ✓ Write instead | source |
|---|---|---|---|
| **1** | any **naive** distance kernel, or any **N2-only** kernel, reported as a distance effect | the N5 technical covariate block is required; report only the conditioned fit | PREREG §10.1 |
| **2** | "**−0.074 SD in negative-control probes**" (or −0.070, or any pooled number called a *probe* number) | "**−0.0744 SD in pooled negative-control features**", and separately: "the 40 named negative-control probes — the pre-registered primary response — are **flat**, −0.0225, p = 0.129" | PREREG §10.1; audit R1 |
| **3** | any **caller-independence** claim | see item 13 | PREREG §10.2 |
| **4** | any **age-stratified or young-vs-old** claim on H1 | age is a continuous covariate only; **two donors are over 55** | PREREG §10.3, PI decision D4 |
| **5** | any **marginal-zone-specific confirmatory** claim | exploratory only | PREREG §10.4 |
| **6** | any cross-arm difference **attributed to species or tissue** | they are confounded by design — mouse liver against human spleen. **This belongs in the abstract, not in Limitations** | PREREG §10.5 |
| **7** | `CXCL8`/`CXCR1` as **replicating a mouse result** | no mouse ortholog exists; `MMP3`/`TIMP1` are mouse-only; `CXCL2`/`CXCL5` are on **both** panels and are a **map gap**, not a biological asymmetry | PREREG §10.6 |
| **8** | the circularity figure "**1.51–2.85×**" | over 11 sections the DeepScence × `Cdkn1a`⁺ pair is **0.963–2.849, median 1.071, pooled 1.255**; both published values were the two largest of the eleven | PREREG §10.7 |
| **9** | the composition-matched number **alone** | wherever the matched **1.6 %** (SF 0.9837) appears, `type_adj` **65.9 %** and `typecomp_adj` **85.4 %** appear beside it. Quoting the first without the second **states the opposite of what the data says** | PREREG §10.8, P23/P25 |
| **10** | `mor` as evidence that normalisation cannot move this caller | §6's named estimator removes only **11.4–23.6 %** of log-depth variance on this panel; wherever it appears, the `lib` result (100 % removed; depth loading −74 % and −93 %; 46 % and 100 % of sender calls changed) appears beside it | PREREG §10.9, P27 |
| **11** | any `denoise=True` number from **a single seed** without its seed-stability companion | one of three seeds gave a top-5 % sender set **disjoint** from the other two (Jaccard 0.000) | PREREG §10.10, P26 |
| **12** | `rho_signed_dz_vs_depth` for the D2 `raw` control rows (**−0.47, −0.16**) | it is the direction of numerical noise on a shift of 0.0002–0.001 z-units | PREREG §10.11, P28 |
| **13** | "0.93–1.22× of chance … i.e. they are **statistically independent**" | see the restatement below | PREREG §10.12 |
| **14** | "**Four of six pairs** sit at 0.93–1.22×" | a property of the two-section base only | PREREG §10.12 |
| **15** | "the one pair that looked **concordant in sham** is anti-concordant in SBR" | there is **no arm effect**; there is one anomalous section (7250) | PREREG §10.12 |
| **16** | "DeepScence's correlation with sequencing depth **reverses sign between two sections** of the same study" | it is a **7250-only** effect at 11 sections | PREREG §10.12 |
| **17** | "one of the four pairs sits **below chance at 0.91×** in all eleven sections" | **DEAD**: frozen value is **0.972, z = −1.63, p = 0.104, above chance in 4 of 11** | `SUBMISSION_PATCH` §2.3 |
| **18** | "**the direction of each pair is predicted by its depth loading**" | **REFUTED**: permutation p = 0.30; within-pair ρ negative in all five pairs | `SUBMISSION_PATCH` §2.3 |
| **19** | "we **discovered** that torus shifts break on non-convex tissue" / "this is a finding in its own right" | see item 21 | `SUBMISSION_PATCH` §0.5, §6 |
| **20** | "**nobody in this literature reports** the negative-control-probe spatial diagnostic" | **FALSE** — Voyager and Ren et al. (2025) both do. See §4.2 | `CITATION_AUDIT.md`; `NOVELTY_ASSESSMENT.md` §U1 |
| **21** | "**the two tests disagree**" (Moran's I vs the A7 kernel) | **FALSIFIED** — they rank together at ρ = +0.895 / +0.944. See §4.3 | `NOVELTY_ASSESSMENT.md` banner |
| **22** | "the reportable-fit filter admits **2–3× more fits** than its nominal rate implies" | **NOT SUPPORTED** — the filter measures **3.0–13.3 %** on the naive design and 4.8 % on the full design; the 9–16 % is the estimator's two-sided CI-exclusion rate | audit R6 |
| **23** | "**circular sender set**" | "**contaminated sender/response split**" or "**hollow sender set**" | `SUBMISSION_PATCH` §0.3 |
| **24** | "CellWHISPER's torus-shift null" / attributing N3 to CellWHISPER | CellWHISPER's null is a **within-cell-type location permutation = our N1**; the torus shift is **Lotwick & Silverman (1982)** | `CITATION_AUDIT.md` §0.1 |
| **25** | CellWHISPER's ">90 % FPR" as a **measured** rate | it is an **inferred interaction-count ratio**; use "implying" | `CITATION_AUDIT.md` §0.2 |
| **26** | "**23 % / 8 % out of tissue**" | out of tissue = `1 − frac_in_occupancy` = **35.5 % / 19.9 %**; 22.8 % / 8.0 % = `1 − frac_retaining_a_neighbour` = "lost every real cell within the 100 µm window". **Never merge the two in one sentence** | audit R3 |
| **27** | "CoreScence circularity **69 % → 76 % → 88 %**" | **mouse 79 % → 88 %, human 76 % → 88 %**; the 69 % was a typed-in literal | audit M1 |
| **28** | any claim that the seven **Tier B modules are mutually disjoint** | Tier A ∩ ⋃B = 0, **but 18 of 21 Tier B pairs share genes** | audit C1 |
| **29** | "the human per-sample cell range **brackets** the mouse range" | it sits **above** it and overlaps only at the top | audit M10 |
| **30** | "all figures verified" from the guard's output | the guard enumerates `git ls-files figures/`; state the count it actually covers | audit M7 |
| **31** | a bracket on an SF or amplitude presented as a **confidence interval** | it is an **IQR across fits**; `sf_summary.csv` has no CI column | `CS_PHASE8_M1_RERUN.md` §6 |
| **32** | Mrkvička "can also be used for **irregular windows**" | the paper's wording is "**can be applied in case of general (compact) observation windows**" | `CITATION_AUDIT.md` §3 |
| **33** | a bare `D<n>` deviation ID | name the series (PI decision / gene-set deviation / `PREREG_PHASE8.md` P-series) | `PHASE8_ROADMAP_STATUS.md` |

## 3.1 The three falsified claims, with their correct wording

### (a) Caller independence — **FALSIFIED**

**Struck** (four sentences, PREREG §10.12): *"0.93–1.22× of chance … i.e. they are
statistically independent"*; *"Four of six pairs sit at 0.93–1.22×"*; *"the one pair that
looked concordant in sham is anti-concordant in SBR"*; *"DeepScence's correlation with
sequencing depth reverses sign between two sections of the same study."*
Also struck, later: the "0.91× below chance in all eleven" plank and the depth-loading
direction rule (items 17–18).

**Correct wording (one sentence, from `SUBMISSION_PATCH` §2.1 — the only version whose planks
all survive):**

> After conditioning on cell type and within-cell-type sequencing depth, four senescence
> callers overlap at **0.75–2.20× of chance across eleven sections, pooling to 1.21×
> (Mantel–Haenszel, z = 21.9)**; they are not independent, but the dependence spans
> **anti-concordance at 0.74× to concordance at 1.47×** with a different identifiable cause at
> each end, which is not what one latent state looks like.

**And the framing that makes it a positive finding, not an accept-the-null:** the sender-set
fix kills independence **on the published two sections** (1.030 → 1.128, p = 4.4 × 10⁻⁸);
coverage then makes it certain. The restatement is **stronger** than what it replaces.

**Do not ship `CS_PHASE8_CALLERS.md` §3's drafted paragraph** — it is pre-C6 and two of its
clauses are dead.

### (b) "The two tests disagree" (Moran's I vs the A7 kernel) — **FALSIFIED**

**Struck** (`NOVELTY_ASSESSMENT.md` §2.1 point 3, still unmarked in that document's body):
*"…state it explicitly with a Moran's I of your own controls alongside, so the reader can see
the two tests disagree."*

**Correct wording:** the two statistics **agree** — over the 12 control and module fields,
\|Moran's I\| and \|A7 naive amplitude\| rank together at **ρ = +0.895 raw, +0.944
cell-type-centred**. "Different question" survives, but justified by **power, not
orthogonality**:

> The entire A7 gradient contributes **ΔI = 2.2 × 10⁻⁴, i.e. 0.83 % of the observed control
> Moran's I**, and the smallest kernel amplitude Moran's I could resolve on these sections is
> **0.362 SD** — five times the A7 control gradient (0.074 SD) and larger than the project's
> own naive biological amplitude (0.277 SD). **Moran's I could not have detected the paper's
> headline effect either.** Within the control family, section by section, the two statistics
> are uncorrelated (ρ = +0.155, p = 0.26), which is what the power calculation predicts.

**Scope note:** `CORRECTIONS.md` §16.1 narrows the falsification — **only §2.1 point 3 was
falsified; §4 O1 was not.** O1's premise (that the aggregate control field is near zero) is
simply gone, since I = +0.0455.

### (c) The A7 response-naming error — **WRONG ATTRIBUTION**

**Struck:** *"a −0.070 SD gradient (p = 0.023) **in negative-control probes**"* — and every
variant of it, including `CORRECTIONS.md` §12's *"the negative-control-probe kernel is −0.074
SD"*, which is the same error with the new number.

**Correct wording:**

> The assay carries a distance gradient in the **codewords and genomic controls**, and in the
> pooled control set (**−0.0744 SD [−0.1306, −0.0182], p = 0.0145**) — but **not** in the 40
> named negative-control probes (**−0.0225 [−0.0527, +0.0078], p = 0.129**), which are the
> sparsest of the three families (0.0067 counts/cell, 0.65 % of cells non-zero) and
> correspondingly the least powered. Report it as "**pooled negative-control features**", never
> as "negative-control probes"; and say which response every number comes from. **Do not let
> "controls" stand unqualified anywhere.**

**Why it bites:** `PREREG_PHASE8_genesets.md` §11 designates `E_negative_control_probes` the
**primary technical null**, and Phase 9 item 9.4 repeats it. **On the pre-registered primary
response, M1's A7 passes naively.** Writing "the assay is not flat in negative-control probes"
states the opposite of the pre-registered test's own result.

---

# 4. FRAMING DECISIONS ALREADY SETTLED

Do not relitigate these mid-draft.

## 4.1 The torus finding is **import-and-quantify**, not discovery

**Verdict (`NOVELTY_ASSESSMENT.md` §2.2):** *"NOT NOVEL IN STATISTICS. NOVEL IN SPATIAL OMICS
AS A DEMONSTRATION. Frame as import, not discovery. The FFT admissible-set computation is a
small genuine methods contribution."*

- The rectangular-window requirement is **Lotwick & Silverman (1982)**, JRSS-B 44(3):406–413.
  It is standard enough to be in the `spatstat` manual: *"The window must be a rectangle.
  Toroidal shifts are undefined if the window is non-rectangular."*
- The remedy for general compact windows is the **variance correction of Mrkvička et al.
  (2021)**, *Spatial Statistics* 42:100430.
- Mrkvička **specifically predicts against our tiling**: extending the approach to windows that
  are finite unions of aligned rectangles *"would increase the amount of cracks in the
  autocorrelation structure. Subsequently, it would increase the liberality of the test of
  independence using random shifts with torus correction."* N3-tile is that construction.
  Presenting as the corrected primary the one variant with a published prediction against it is
  an avoidable own goal — which is why **N3-var is primary**.

**The claim, stated correctly:** *a 40-year-old documented limitation is being violated in
current spatial-omics practice, and we quantify what it costs.* **Three things in that sentence
are ours:** (i) the cost measured on real tissue (35.5 % / 19.9 % out of tissue; 1–66 of
38,080–108,375 offsets admissible; one section admitting only the identity); (ii) a **direct
calibration measurement**, which no spatial-omics paper has (tiled torus 0.080–0.118 vs
nominal 0.05, up to 2.4×; RS_count holds 0.033–0.060); (iii) the **FFT enumeration** — **lead
the methods paragraph with it; it converts a rediscovery into a tool.**

*"Presenting the N3 degeneracy as a discovery rather than as the first quantification on real
tissue of a known classical pathology is the single easiest way to lose a statistician
reviewer."*

## 4.2 The negative-control-probe **diagnostic** is not novel; refitting the estimand's own estimator is

**Prior art, cite it approvingly.** The Voyager Xenium vignette (Moses et al., bioRxiv
2023.07.20.549945, **still a preprint**): *"generally the negative controls are tightly
clustered around 0, while the real genes have positive Moran's I, which means there is
generally no technical artifact spatial trend."* And **Ren P et al. (2025)**, *Nat Commun* 16,
doi:10.1038/s41467-025-64292-3: *"Spatial autocorrelation analysis using Moran's I revealed
stronger aggregation of negative control signals in CosMx 6K …"*

**What survives, exactly (`NOVELTY_ASSESSMENT.md` §2.1 item 1):**

> **No paper fits the estimand's own estimator to the negative control features.** Everyone
> computes a *generic* spatial statistic (Moran's I, per-cell rate) on negative controls.
> Nobody runs the actual model under test — here, a distance-to-nearest-sender kernel
> regression — with negative controls as the response. That is a **negative control outcome for
> the specific estimand**, which is **Lipsitch's construction** (Lipsitch, Tchetgen Tchetgen &
> Cohen 2010, *Epidemiology* 21(3):383–388), and I found no instance of it in spatial omics.

**Second, genuinely new point from the Moran run:** Voyager's controls-vs-genes contrast is an
**abundance** contrast. Genes matched to the controls on total counts give **identical**
Moran's I (−0.00018 vs −0.00012); a negative control probe carries ~21 counts per section
against a median gene's 5,885. **The per-feature statistic has no power at control abundance.**

## 4.3 "Different question" is justified by **power, not orthogonality**

See §3.1(b). This replaces an asserted orthogonality with a measured power bound, and it is
stronger for it. `NOVELTY_ASSESSMENT.md` carries the correction banner; its **body still
contains the falsified sentence unmarked** at §2.1 point 3. A drafter reading that document
linearly will hit the dead sentence.

## 4.4 N2-vs-N5 is the strongest contribution — **three independent lines**

An external novelty review ranks it **#1 of six** contributions: *"Genuinely unreported,
mechanistically explained …, directly transferable to every imaging-ST distance analysis, and
it overturns the project's own prior that N2 is the conservative gold standard. Highest
surprise-per-word in the repo. Costs one table."* **Give it its own subsection and a table,
not a footnote.**

**The three lines (`CORRECTIONS.md` §13 / §C.1):**

1. **Figure 1 — synthetic, planted kernel, ground truth.** λ̂ CI coverage where ℓ ≥ 2λ:
   matched-decoy **0.346** vs nuisance-conditioned **0.854**, with naive at 0.5125 — so the
   matched decoy is **worse than doing nothing**. **This is the strongest of the three, because
   the kernel is planted. It has been in the project's headline figure since Phase 1 and was
   never read off it. Lead with it.**
2. **A7 — the platform's own technical gradient.** N2 leaves it: **−0.0642, p = 0.0124**
   (86 % undiminished). N5 removes it: **+0.0038, p = 0.715**.
3. **The composition-matched protocol.** The same variables, on the same fits, remove
   **1.6 %** as a matched decoy set and **85.4 %** as covariates — **a factor of fifty**.

*"The first two are answerable with 'your covariates removed real signal.' The third is not,
which is why it leads."*

**Note the counts in circulation:** `NOVELTY_ASSESSMENT.md` knows one line;
`SUBMISSION_PATCH` §0.4 says "twice, independently". **Three is current.**

## 4.5 Venue — evidence gathered, **decision is the PI's and is not taken**

`NOVELTY_ASSESSMENT.md` §5 argues the master plan's primary/secondary assignment is **backwards
for the result actually obtained**, and recommends **ICBINB-BIO primary (8-page full track),
ml4spatialbio secondary (4-page methods-instrument distillation)**. Verified facts:
ml4spatialbio is **4 pages, non-archival, concurrent submission explicitly permitted**, window
**Aug 29 – Sept 4 AoE**; ICBINB-BIO is **8 pages (full) or 4 (tiny)**, deadline **Aug 29 11:59
AoE**, and explicitly wants candid failure analysis and negative results. ICBINB-BIO's
**dual-submission policy was NOT verified.** `SASP_Kernel_Master_Plan.md` §29 states the
assignment is unchanged and remains the PI's call. **This must be a recorded decision, not a
drift.**

## 4.6 Other settled positions

- **The research question is unchanged.** "How far does senescence signalling reach?" is
  answered **in the negative, as a bound**. The audit battery is *supporting evidence for why
  the bound is the honest answer*, not a replacement thesis. (PI decision "Direction".)
- **§18 outcome A stands:** *no distance-dependent SASP kernel is identifiable at achievable
  power* — now on a **tighter** bound (0.183) than the published one (0.203).
- **Tier A strict-33 is PRIMARY**; the seven per-module sets are the pre-registered
  sensitivity. (PI decision D1.)
- **`denoise=False` is the frozen primary**, with `denoise=True` as the published-default
  sensitivity, carrying its seed-stability companion. (Task 8.5.)
- **Age is a continuous covariate only.** (PI decision D4.)
- **CoreScence 88 % circular: strip-and-refit is primary.** (PI decision D6.)
- **The introduction must state, in the introduction and not in Limitations, that no published
  length constant with uncertainty exists for this quantity** — so the contribution is a bound
  and an identifiability argument, not the overturning of a number. This is a strength and a
  risk, and it is why Figures 1 and 4 are load-bearing rather than supporting.
- **Do not claim novelty for "senescence callers disagree."** DeepScence, SenCID, SenePy, ICE,
  markeR and Ntintas et al. all report it. Claim only the **coverage-and-definition
  sensitivity**, and only with the decomposition attached. Cite **ICE** (*Genome Biology* 2026,
  doi:10.1186/s13059-026-03997-0), **markeR** (*NAR Genomics and Bioinformatics* 2026,
  8(2):lqag057) and bioRxiv 2026.01.02.697374 (*"apparent concordance in prior studies may
  reflect circular validation"* — the closest published statement to our §5.9).

---

# 5. WHAT IS NOT YET MEASURED

**Nothing in this section may be written as if it existed.**

| item | status | blocked by |
|---|---|---|
| **Everything H1 / the human arm** | **NOT RUN.** Data acquired (28 files, 525 MB) and **structurally verified only** — panel, cell_id sets, coordinate units, annotation levels. **No expression value has been read from `data/raw_h1/`.** | Phase 8 tag + pre-registration commit |
| Phase 9: A1 (resolution, segmentation, assignment rate), A4 (Ripley's K), A8 (cross-arm on the ortholog-intersected panel) | NOT RUN | the freeze |
| Phase 9: **A2 gate** (disjointness on the real panel — pre-verified as passing) | NOT RUN | the freeze |
| Phase 9: **A5 gate** (matched-decoy contrast, \|SMD\| ≤ 0.1) | NOT RUN | the freeze |
| Phase 9: **A7 on H1** (negative-control-probe kernel, must be flat) | NOT RUN — the **mouse** half is complete | the freeze |
| Phase 9: A6 red/white pulp covariate | spec written, not built | the freeze |
| Phase 9: Job B — 22-type spleen cell typing; cross-check against depositors' 4-level annotations; 4 sender callers | NOT RUN | the freeze |
| Phase 9: A3 prevalence per cell type (1–20 % band); caller agreement conditioned on type + depth decile, human | NOT RUN | the freeze |
| Phase 10: H1 through the frozen pipeline (naive, N1–N8, controlled fits, kernel families, superposition vs nearest, proximal vs downstream) | NOT RUN | Phase 9 |
| Phase 10: composition-matched reruns on H1, 5 seeds, both arms | NOT RUN — the protocol is implemented and arm-generic; `--arm h1` **refuses to run** | Phase 9 |
| Phase 10: the §8 DeepScence native-vs-remapped comparison against the pre-registered prediction | NOT RUN | Phase 9 |
| **The §17 two-arm table** | **ONE COLUMN ONLY.** Every mouse cell is filled; every human cell is empty | Phase 10 |
| **Figure 3 revised** (controlled kernels, both arms) | WAIT | H1 |
| **Figure 5** (two-arm replication: SF by null both arms; forest plot; geometric predictions) | **DOES NOT EXIST** | Phase 10 |
| **Figure 6** (DeepScence native vs ortholog-remapped; caller agreement before/after conditioning; whether the anchor instability transfers) | **DOES NOT EXIST**; 6a partly unblocked by the D1 coverage work | Phase 10 |
| Figure 2h, human half | mouse half complete (`figure2e`); human half behind the freeze | the freeze |
| A7 on the M1 **N7 sender axis** under N3-var | not run — N3-var covers `tierA_p95` on the 6 in-band sections only | — |
| Shift-disk radius sensitivity; `RS_var` (the third variance estimator); the point-process branch (cross-K, cross-nearest-neighbour, RS_K/RS_G) | not run | — |
| Running COMMOT / CellChat / SpaTalk **on the control probes** | not done — would turn the FPR self-check into a benchmark and raise §5.5's rank | — |

**Stop conditions (pre-registered):** A2 or A5 failure ⇒ **§18 outcome C**, a data-availability
finding.

**Citations that are verified as records but whose content is UNREAD** (`CITATION_AUDIT.md` §3)
— do not paraphrase their content: Ma et al. (2024) per-pathway distance-gradient panel;
Acosta et al. 2026 Author Correction; Hodges & Reich (2010) primary text; Zimmerman & Ver Hoef
(2022) primary text; Lotwick & Silverman (1982) primary text (**the toroidal-shift attribution
is verified only via Mrkvička §1 — attribute nothing beyond it**); Hanks et al. (2015); the
SenNet "1,753 / January 2026" figure (**unverifiable; the portal says 2,041 as of April 2026**);
NCEM's "characteristic length scales" phrasing; `gu2026identifiability` content.
Also: `neretti2024dissecting` is a **single-author conference abstract** — replace or drop the
claim; `suryadevara2026charting` is a **Review** — do not cite as primary evidence.

---

# 6. SOURCE MAP — which report and which result files back each subsection

| §30 subsection | primary report(s) | result files | producing code |
|---|---|---|---|
| **Methods / cohort** | `PREREG_PHASE8.md` §2–§3, `PHASE7_H1_SCREEN.md`, `PHASE0_DATA_AND_ENV.md` | `results/composition_by_arm_timepoint.csv`, `results/phase3/window.csv`, `genesets/h1_candidate/` | `code/prepare_samples.py`, `code/annotate_pipeline.py`, `code/screen_candidate_panels.py` |
| **Methods / gene sets** | `PREREG_PHASE8_genesets.md`, `BIO_PHASE8_FREEZE.md`, `BIO_PHASE7_JobA*.md` | `genesets/*.txt`, `genesets/human/FROZEN_MANIFEST.csv`, `figures/figure_gs1..gs4_*_data.csv` | `code/build_genesets*.py`, `code/gate_genesets_guard.py`, `code/gate_disjointness_human.py`, `code/make_figure_genesets.py`, `code/corescence_circularity.py`, `code/crossarm_geneset_table.py` |
| **5.1 identifiability (Fig. 1)** | `CS_PHASE1.md`, `CS_PHASE2.md` | `figures/figure1_data.csv`, `results/summary_tables.txt`, `results/sweep_all.csv`, `results/phase1b/` | `code/sasp_sim.py`, `code/sasp_sweep.py`, `code/sasp_phase1b.py`, `code/make_figure1.py` |
| **5.2 naive gradient + nulls (Fig. 2a–d)** | `CS_PHASE3.md`, `CS_PHASE8_M1_RERUN.md` | `results/phase3/main_fits.csv`, `sf_summary.csv`, `summary_phase3.txt`, `perm_nulls*.csv`, `m1_final_audit.txt`, `figures/figure2a_*.csv`, `results/phase3/figure2b_data.csv`, `figure2c_data.csv`, `figure2d_data.csv` | `code/run_phase3_nulls.py`, `code/summarize_phase3.py`, `code/m1_headlines.py`, `code/m1_final_audit.py`, `code/make_phase5_figs.py --which 2a`, `code/make_figure2bc.py` — **never `code/make_figure2.py`** |
| **5.3 controlled kernels (Fig. 3)** | `CS_PHASE3.md`, `CS_PHASE5.md`, `CS_PHASE8_M1_RERUN.md` | `results/phase3/m1_final_audit.txt`, `main_fits.csv`, `results/phase5/kernel_families.csv`, `kernel_heldout.csv`, `summary_phase5.txt`, `results/phase3/window.csv`, `lamscale*.csv` | `code/run_phase5_kernels.py`, `code/summarize_phase5.py`, `code/run_phase3_lamscale.py`, `code/make_phase5_figs.py` |
| **5.4 existing tools (Fig. 4)** | `CS_PHASE4.md`, `CITATION_AUDIT.md` §0.1–0.2 | `results/phase4/headline.csv`, `headline_by_split.csv`, `headline_offdiagonal.csv`, `ncem_radius_sweep.csv`, `commot_mechanism.csv`, `figures/figure4_data.csv` | `code/phase4_run.py`, `code/phase4_summarize.py`, `code/make_figure4.py` (`load_ours()`), `code/make_figure4_supp*.py` |
| **5.5 estimator on the controls (Fig. 2h) + Moran** | `CS_PHASE8_CALLERS.md` §4 *(pre-C6 digits)*, **`CS_PHASE8_MORAN.md`**, `AUDIT_PHASE8_FACTCHECK.md` R1/M3/R6 | **`results/phase3/a7_summary.csv`**, `a7_verdict.txt`, `a7_control_probe_{fits,curves,provenance}.csv`, **`results/moran/moran_{pooled,vs_a7,sensitivity,per_feature_summary,kernel_power}.csv`**, `moran_verdict.txt`, `results/phase3/figure2e_data.csv` | `code/run_a7_control_probes.py`, `code/summarize_a7.py`, `code/run_moran_controls.py`, `code/run_moran_lognorm.py`, `code/moran_kernel_power.py`, `code/summarize_moran.py` |
| **5.6 N2 vs N5** | **`CS_PHASE8_COMPMATCH.md`**, `CORRECTIONS.md` §13, `NOVELTY_ASSESSMENT.md` §3 | **`results/phase3/compmatch_reruns.csv`**, `compmatch_fits.csv` (+`_tierA`,`_tierApm`), `results/phase3/a7_summary.csv`, `figures/figure1_data.csv`, `results/phase3/attribution.csv` | `code/run_phase8_compmatch.py`, `code/_compmatch_chain.sh`, `code/run_phase3_attribution.py` |
| **5.7 random-shift nulls** | **`CS_PHASE8_TORUS_VAR.md`**, `CS_PHASE7_C1.md` *(pre-C6 digits)*, `CS_PHASE8_C1_CLOSEOUT.md` | **`results/phase3/sf_summary_var.csv`**, `summary_phase3_var.txt`, `null_destructiveness{,_var}.csv`, `var_sim_calibration.csv`, `var_pvalues.csv`, `var_variance_check.csv`, `perm_nulls_var{,_full200}.csv`, `perm_draws_var.csv`, `sf_summary_c1.csv` | `code/phase3_null_geom.py`, `code/phase3_null_diag.py`, `code/phase3_null_var.py`, `code/run_phase3_var.py`, `code/summarize_phase3_var.py`, `code/phase3_var_validate.py`, `code/phase3_var_sim.py` |
| **5.8 two-arm (Figs. 5, 6)** | `PREREG_PHASE8.md` §6–§8, `PHASE7_H1_SCREEN.md`, `Phase7_Minimal_Human_Replication (1).md` §19 | **none — not run.** Mouse-side geometry: `results/phase3/poisson_fits.csv`, `ripley.csv` | Phase 10 |
| **5.9 what we withdrew** | **`CS_PHASE8_CALLERS.md`** *(pre-C6 digits)*, `SUBMISSION_PATCH_2026-08-29.md` §2–§3, `CS_PHASE8_D2_DENOISE.md`, `CORRECTIONS.md` §14 | **`results/phase3/caller_coverage_gate{,_headline}.csv`**, `caller_agreement_matched_significance_{verify2sec,2sec_c6,11sections}.csv`, `caller_within_type_depth_bias_11sections.csv`, `caller_technical_loading_11sections.csv`, `caller_celltype_composition_11sections.csv`, `caller_agreement_matched_d3_11sections.csv`, `deepscence_anchor_decisions.csv`, `figures/figure_phase8_callers_data.csv`, `figure_phase8_d3_data.csv`, **`results/phase8_d2/d2_*.csv`** | `code/caller_disagree_all.py --all` / `--verify`, `code/summarize_caller_coverage.py`, `code/deepscence_reanchor.py`, `code/caller_disagree_d3.py`, `code/analyse_d2_denoise.py`, `code/analyse_d2_stability.py`, `code/report_d2_tables.py` |
| **6 corrections & pre-registration** | **`CORRECTIONS.md`**, `PREREG_PHASE8.md`, `PREREG_PHASE8_genesets.md`, `AUDIT_PHASE8_FACTCHECK.md`, `AUDIT_CORRECTIONS_APPLIED.md`, `COMPLETED_TASKS.md` | `results/phase3/m1_{final_audit,prepost_main_fits,n7_prepost}.txt`, `results/phase3_pre_c6/` (the baseline the ledger compares against), `figures/.committed_manifest.json` | `code/m1_final_audit.py`, `code/m1_compare_modules.py`, `code/check_figures_guard.py` |
| **7 discussion** | `BIO_DELIVERABLE6_DISCUSSION.md`, `BIO_DELIVERABLE7_CLAIM_AUDIT.md`, `NOVELTY_ASSESSMENT.md` §2–§3, `SASP_Kernel_Master_Plan.md` §29 | — | — |
| **8 appendix** | `CS_PHASE3.md`, `CS_PHASE5.md`, `PREREG_PHASE8_genesets.md` §12 | `results/phase3/n8_*.csv`, `stratification.csv`, `correlogram.csv`, `tierC_*.csv`, `figures/figure_gs1_intersection_matrix_data.csv`, all `sf_summary*` | as above |

## 6.1 Reproduce the whole M1 arm

```bash
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -u phase2_downstream.py sham sbr <11 section ids>      # gene sets -> sender calls, module scores
python3 -u _run_prep.py                                        # Phase 3 cache, forced
python3 -u run_phase3_nulls.py --stage window --sections all
python3 -u run_phase3_nulls.py --stage main   --sections all --calls all9 --n-jobs 24
python3 -u phase3_null_diag.py
python3 -u run_phase3_nulls.py --stage perm    --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
bash _m1_rerun_stage2.sh     # N7 at 1000, tierApm_p95 at 1000, curves
bash _m1_rerun_stage3b.sh    # the eight small Phase 3 scripts (N8 fanned out over sections)
bash _m1_rerun_stage4.sh     # Phase 5
python3 -u caller_disagree_all.py --all
bash _m1_rerun_stage5.sh     # A7, summaries, the ONE figure pass
python3 code/check_figures_guard.py
# NOT make_figure2.py -- superseded producer of figure 2a; it now refuses to run
```
Baselines for every pre/post comparison: `results/phase3_pre_c6/`, `results/phase5_pre_c6/`,
and `git tag pre-c6-genesets`. `MASTER_SEED = 20260820`.
**`results/phase3_pre_c6/` is the sole copy of that baseline — 98 files, 0 tracked.**

---

# 7. RESIDUAL DISAGREEMENTS BETWEEN DOCUMENTS, AND THE AUTHORITATIVE FILE FOR EACH

> ## ✅ ALL 22 WORKED, 2026-08-27 — see `reports/RECORD_RECONCILIATION.md`
>
> The table below is the *diagnosis*. The fixes have now been applied to the documents, and
> the per-file record of what changed and why is in `reports/RECORD_RECONCILIATION.md`.
> Status by row:
>
> | rows | status |
> |---|---|
> | **1, 4, 5, 8, 9, 12** — pre-C6 digits in `PHASE8_ROADMAP_STATUS.md`, `PREREG_PHASE8.md`, `CS_PHASE7_C1.md`, `NOVELTY_ASSESSMENT.md` | **FIXED.** Frozen values in place; pre-C6 digits struck through or marked, never deleted |
> | **2, 3, 13** — caller bases, the dead Tier A × SenePy plank, "22 of 33" | **FIXED.** `CS_PHASE8_CALLERS.md` now carries a whole-file pre-C6 banner and the arithmetic error (20 of 33, and 26 of 33 frozen) is corrected |
> | **6** — naive biological amplitude | **RESOLVED: +0.2767**, section-clustered signed mean. See §5.5 |
> | **7** — Spearman ρ | **RESOLVED: +0.8951 frozen, aggregation stated everywhere.** +0.923 reproduces as the median-per-field aggregation. See §0.4 |
> | **10, 21** — RS_count 0.033–0.060; tiled-torus inflation | **FIXED.** 0.040–0.060 corrected in `CS_PHASE8_M1_RERUN.md` §6/§14.3; the inflation is standardised repo-wide on **2.35×**, the exact 0.1175 / 0.05 — "2.4×" and "2.36×" both rounded it up |
> | **11, 16** — "23 % / 8 % in the void"; N4 8.3 % | **FIXED.** Audit R3's column fix applied to `CS_PHASE8_TORUS_VAR.md` (§1 table and §10 blockquote) and to `NOVELTY_ASSESSMENT.md` §2.2; audit R3's own 8.3 % marked superseded to 8.0 % |
> | **14** — figure guard | **CLOSED.** `git ls-files figures/` = 52, on disk 52, `check_figures_guard.py` exits 0 with "all 52 committed figures match". `PHASE8_ROADMAP_STATUS.md`'s "27 of 45" paragraph is replaced |
> | **15** — matched decoy worse in "6 of 8" | **FIXED to 5 of 8 strictly** (1 exact tie at 0.000), `CORRECTIONS.md` §13.1 |
> | **17** — audit R5 | **MOOT, and said so.** `neg_probe_rate` under N6+N5 is +0.0097 [−0.0060, +0.0253], p = 0.199. "Every control family is flat under +N6+N5" is now **true**; the hedge is withdrawn in `CS_PHASE8_CALLERS.md` §4.1 and `AUDIT_PHASE8_FACTCHECK.md` R5 |
> | **18, 19** — stale "NOT DONE" statuses; the `SUBMISSION_PATCH` §9 warning | **CLOSED** in `NOVELTY_ASSESSMENT.md` §4 and `CORRECTIONS.md` §1/§18.2 |
> | **22** — the composition-surrogate row | **WITHDRAWN.** The 76 % is 0.212/0.260 = 0.815 with an unsourced denominator; `CS_PHASE5.md` §4 now says so, and the row is replaced by the 65.9 % / 85.4 % composition-matched pair |
> | **20** — bibliography counts | **HALF RESOLVED.** The count is **43 entries** — `grep -c '^@' references.bib` — so "32 entries" is a stale vintage and `CORRECTIONS.md` §18.3's "truth: 43" is right. The other half ("41 wrong given names across 19 of 32 entries") is an **audit finding against the 32-entry vintage** and cannot be restated for the 43-entry file without re-checking every author line against Crossref/PubMed. **Not done here** — see `RECORD_RECONCILIATION.md` §6 |


| # | disagreement | authoritative source | resolution |
|---|---|---|---|
| 1 | headline vector: 0.326/0.027/0.082/0.203 vs **0.329/0.029/0.088/0.183** | `results/phase3/m1_final_audit.txt` | the second is post-C6; the first is pre-C6 (still in `PHASE8_ROADMAP_STATUS.md`'s PI-decisions table) [V] |
| 2 | caller agreement: 1.118 (p 1.4e-30) vs 1.128 (4.4e-8) vs **1.212** (1.8e-106 / 1.8e-94) | `results/phase3/caller_coverage_gate_headline.csv` | six different bases; see §5.9. 1.118 is 11-section **pre-C6, 3-pair**; 1.212 is 11-section **post-C6** [V] |
| 3 | Tier A × SenePy: 0.914, below chance 11/11 vs **0.972, p = 0.104, 4/11** | `results/phase3/caller_coverage_gate.csv` | pre-C6 vs frozen. The plank is dead [V] |
| 4 | A7 pooled: −0.070 (p 0.023) vs **−0.0744 (p 0.0145)**; probes −0.018 (p 0.183) vs **−0.0225 (p 0.129)** | `results/phase3/a7_summary.csv` (09:06) | pre-C6 vs frozen. `PREREG_PHASE8.md` §10.1 para 1 is pre-C6; para 2 is frozen [V] |
| 5 | FPR "9–13 %" vs **"9–16 % (9–15 % count-based)"** | `results/phase3/a7_summary.csv` | 9–13 % rests on the pre-C6 probe value 0.127; frozen is 0.145 [V] |
| 6 | naive biological amplitude: 0.277 / 0.291 / 0.312 / 0.314 | `results/phase3/a7_summary.csv` (**0.2767**, clustered mean) and `a7_verdict.txt` (**0.3120**, median \|β\|/sd) | **two estimators, two vintages.** Name the estimator; the power argument holds for all four [V] |
| 7 | Spearman ρ raw: **+0.895** vs +0.923 | `results/moran/moran_verdict.txt` | **+0.8951**, reproduced from two independent files. +0.923 is not reproducible [V] |
| 8 | N3-tile 0.974 vs **0.971**; published N3 1.000 vs **0.999** | `results/phase3/sf_summary_c1.csv`, `sf_summary_var.csv` | 160-fit (pre-C6) vs **153-fit (frozen)** population [V] |
| 9 | admissible offsets **1–66 / 28 µm** vs 1–63 / 27 µm | `results/phase3/null_destructiveness.csv` | frozen vs pre-C6 [V] |
| 10 | RS_count calibration **0.033–0.060** vs 0.040–0.060 | `results/phase3/var_sim_calibration.csv` | the minimum is 0.033 (rectangle, s = 0.02). `CS_PHASE8_M1_RERUN.md` §6/§14.3 is wrong [V] |
| 11 | "23 %/8 % in the void" (`CS_PHASE8_TORUS_VAR.md` §1, §10; `NOVELTY_ASSESSMENT.md` §2.2) vs **35.5 %/19.9 %** | `results/phase3/null_destructiveness.csv` | audit R3 was applied to three documents but **not** to `CS_PHASE8_TORUS_VAR.md` [V] |
| 12 | "3-pair band" label on a 4-pair row | `results/phase3/caller_coverage_gate_headline.csv` (`n_pairs = 4`) | `PHASE8_ROADMAP_STATUS.md`'s 8.4 table mislabels its 6-section row [V] |
| 13 | "22 of 33 sections above chance" | re-derivation from the significance CSV | it is **20 of 33** (0/11 + 11/11 + 9/11). `CS_PHASE8_CALLERS.md` §2.1 arithmetic error |
| 14 | figure guard "27 of 45" vs "all 52" | `git ls-files figures/` = **52**, on-disk = 52 | the Phase 8 figures were committed; audit M7's scope criticism is closed [V] |
| 15 | matched decoy worse than naive in "6 of those 8" cells | `figures/figure1_data.csv` | **5 of 8 strictly**; 6 of 8 only if the κ=0, ℓ/λ=4 tie at 0.000 counts [V] |
| 16 | N4 "lost all neighbours" 8.3 % vs **8.0 %** | `results/phase3/null_destructiveness.csv` | 0.0805 ⇒ **8.0 %**. Audit R3's 8.3 % predates the 05:58 regeneration [V] |
| 17 | `neg_probe_rate` under N6+N5: "CI excludes zero, p = 0.020" | `results/phase3/a7_summary.csv` | **frozen: +0.0097 [−0.0060, +0.0253], p = 0.199 — CI now includes zero.** Audit R5 is moot; "every control family is indistinguishable from zero under +N6+N5" is now **true** [V] |
| 18 | `NOVELTY_ASSESSMENT.md` §4 O4/O11/O12 marked "NOT DONE" | `CITATION_AUDIT.md` §0.1, §1 | done: nine locations fixed, eleven spatial-statistics entries added. **Do not carry the "NOT DONE" statuses forward** |
| 19 | `CORRECTIONS.md` §D says `SUBMISSION_PATCH` §9 "still instructs the falsified Moran framing" | `SUBMISSION_PATCH_2026-08-29.md` L606 | the dated correction block **is** there; the warning has been overtaken |
| 20 | bibliography "32 entries, 19 with invented forenames" vs "43 entries, 18 author lines corrected" | `references.bib` | count it before citing either. **41 wrong given names across 19 of 32 entries** is the audit's finding; the remedy is procedural — copy author names from Crossref/PubMed, never type them |
| 21 | tiled-torus inflation 2.4× vs 2.35× | `results/phase3/var_sim_calibration.csv` | 0.118/0.05 = **2.36**. Either rounding is defensible; be consistent [V] |
| 22 | `CS_PHASE8_M1_RERUN.md` §6 says the composition-surrogate row **is** sourced; §9 item 6 says it is not | — | §9 is stale relative to §6. Both are superseded by §0.2 above: **the 76 % is unsourced and the row is replaced** |

**Two files still carry live stale numbers that nobody has fixed:**
`reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` L285 carries *"0.93–1.22× of chance for four of six
pairs, but 1.51–2.85× for DeepScence vs `Cdkn1a`⁺"* as a **live recommendation** — both halves
are forbidden (checklist items 8, 13, 14). And `reports/CS_PHASE7_C1.md` §6.3 is the last place
the torus finding is presented as novel statistics (checklist item 19).

---

*Written 2026-08-27. Every number marked [V] was re-derived from `results/` in the session that
produced this file; the derivation commands are inline. Nothing in `results/`, `figures/`,
`code/`, `genesets/`, `data/` or any other report was modified.*
