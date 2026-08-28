# Phase 10 — the two-arm comparison

**Arms.** **M1** = GSE310392, mouse liver (SBR fibrosis), Xenium Prime Mouse 5K + 100 custom,
5,097-gene panel, 11 sections / 11 mice, **6 admissible** under §8 Test 3.
**H1** = GSE326743, human spleen, Xenium Prime Human 5K + 100 add-on, 5,093-gene panel,
7 sections / 7 donors, ages 17–59, **all 7 analysed**.
**Freeze.** `phase8-frozen` = **`d04691e2692a7be8d1ff676d2fb74ad9d1df049d`**, re-cut
2026-08-27 19:55:58 UTC. (The pre-registration's §1 previously recorded the first cut,
`9264396…`, which is an ancestor of it; `PREREG_PHASE8.md` §0 now carries the current hash.)
**Date.** 2026-08-27.

> **Every number in this report was read from a file by the command printed beside it.**
> Where a producer has not finished, the row says **NOT RUN** and why, rather than being
> omitted. The machine-readable version of the §17 table is
> `results/phase10_h1/two_arm_table.csv`, one row per (quantity, arm, panel) with the
> absolute source path and the exact filter in every row:
> `python3 code/two_arm_table.py`.

---

## 0. Which §18 outcome the two arms jointly support

**The two arms jointly support outcome A together with outcome D, with F holding
orthogonally. Outcomes B, C and E are not triggered.** The pre-registered replication
criteria do not, however, resolve cleanly, and that is stated first because it is the part
that runs against the tidy version of this result.

| §18 | Statement | Verdict | Basis |
|---|---|---|---|
| **A** | Both arms agree, no kernel above bound | **SUPPORTED** | M1's median controlled amplitude under N2+N5+N6 is **0.0288** response-SD against its own 80 %-power detectable bound of **0.1833**; H1's is **0.0307** against its own bound of **0.1094**. **0 of 7 Tier B modules exceed the bound on either arm.** |
| **B** | H1 shows a surviving kernel (criterion R4) | **NOT TRIGGERED** | R4 requires *both* that the SF interval exclude 0 *and* that the controlled amplitude exceed H1's own bound in ≥ 2 of the 7 Tier B modules. The first holds; **the second fails 0 of 7**, at both sender calls. |
| **C** | H1 fails A2 or A5 | **NOT TRIGGERED** | Phase 9: A2 byte-identical to the frozen log; A5 passes 35/35, max \|SMD\| after matching 0.0933. |
| **D** | The confound structure differs but the null result does not | **SUPPORTED, and this is the informative half** | Seven structural differences, §4. The largest: the residual after the full nuisance design is **0.158 on H1 against 0.088 on M1**, the reportable-fit fraction is **16.9 % against 48.6 %**, and λ̂ rails at the **ceiling** on H1 (44.9 % of fits) where it rails at the **floor** on M1 (32.7 %). |
| **E** | C1 changes M1's N3 result | **ALREADY RESOLVED, unchanged** | It does not. N3-var 0.9960, N4-var 0.9851, N3-tile 0.9706 against the published bounding-box 0.9989 / 0.9472. |
| **F** | DeepScence's instability appears in H1 too | **SUPPORTED, and it is worse natively** | Full-section seed 0 vs 1 on H1: Pearson *r* **0.372**, top-5 % Jaccard **0.211**, against an M1 floor of **0.9955 / 0.7606** and an M1 same-seed determinism control of *r* = 0.9999991. |

### 0.1 The pre-registered criteria do not resolve cleanly, and neither reading is chosen here

**R2 passes decisively.** H1's median controlled amplitude under N2+N5+N6 is **0.0307**
response-SD (primary call) / **0.0345** (frozen-literal call) against **H1's own** 80 %-power
detectable bound of **0.1094** / **0.1083** — a factor of 3.1–3.6 below it — and **no Tier B
module's median controlled amplitude exceeds that bound**, on either call.

**R1 fails, narrowly, on one of its two clauses.** R1 requires the inter-quartile range of the
SF under N2+N5+N6 across reportable fits to **include 0** *and* to have an **upper quartile
below 0.50**. On H1 the IQR is **[+0.0074, +0.3837]** at the primary call and
**[+0.0040, +0.4759]** at the frozen-literal call. The upper-quartile clause passes on both.
**The "includes 0" clause fails on both — by 0.007 and 0.004 of a surviving fraction**, i.e.
by less than one percent of the naive amplitude.

**R4, the pre-registered non-replication criterion, is not met either.** So H1 satisfies
neither the replication criterion as literally written nor the non-replication criterion.
**R1 and R4 are not exhaustive, and this data falls in the gap between them. That is a
pre-registration defect, not a result, and it is reported as one rather than resolved by
picking whichever reading is convenient.** The PI must decide which of the two clauses of R1
is the criterion. Both readings are in `results/phase10_h1/two_arm_table.csv`; nothing below
depends on the choice, because R2 and R4 agree with each other and with M1.

**Which fit population R1 was scored on is also an open pre-registration question**
(`PREREG_PHASE8.md` §0.1 open item 1). §5, §6 R1 and §6 R2 name `tierA_p95`; PI decision D-B
makes the merged call primary on H1. **Both are computed and both are reported everywhere in
this document**, and they agree: SF 0.1582 vs 0.1630, controlled amplitude 0.0307 vs 0.0345,
0 of 7 modules over the bound either way.

### 0.2 What the answer is, in one paragraph

On two species, two tissues, two laboratories and two independently annotated panels, the
distance-to-nearest-sender amplitude is **large before conditioning and small after it**, and
after the full nuisance design it sits **well below what either arm is powered to detect**.
M1 removes 91 % of the naive amplitude with N2+N5+N6; H1 removes 84 %. The arms do **not**
agree on how much residual is left (0.088 vs 0.158) or on which nuisance carries it, and
they do not agree on the geometry of the fitted length constant at all. **The null result
replicates; the confound structure does not.** That combination is outcome D on top of
outcome A, and it is the argument for running the whole battery on every dataset rather than
importing a characterisation from another one.

---

## 1. What was run, and on what

### 1.1 The H1 arm through the frozen pipeline

The estimator is not reimplemented anywhere in Phase 10. `code/h1_phase10.py` rebinds three
module-level constants — the section cache, the section list and the results directory — and
`code/h1_run_phase10.py` then calls `run_phase3_nulls`'s own job functions
(`_section_job`, `_perm_job`, `_perm_c1_job`) and `run_phase3_var._var_job`. This is the same
code path M1's re-run used (`reports/CS_PHASE8_M1_RERUN.md`), including the reconstructed
`perm_c1` stage, `--tag` and the `all9` call expansion.

```
python3 -u code/h1_run_phase10.py --stage window
python3 -u code/h1_run_phase10.py --stage main    --calls all12 --n-jobs 8
python3 -u code/h1_run_phase10.py --stage perm    --calls primary2 --n-perm 1000 --n-jobs 7
python3 -u code/h1_run_phase10.py --stage perm_c1 --calls primary2 --n-perm 1000 --n-jobs 14
python3 -u code/h1_run_phase10.py --stage var     --calls primary2 --n-perm 1000 --n-jobs 7
python3 -u code/h1_run_phase10.py --stage diag    --calls tierAmg_p95
python3 -u -c "import h1_phase10, run_phase3_poisson as RP; RP.CALLS += ['tierAmg_p90','tierAmg_p95','tierAmg_p99']; RP.main()"
```

**Three H1-specific choices, all declared here and nowhere else.**

1. `strata_mode = "none"`. The mouse `zonation` stratification loops over the three
   *hepatocyte* zonation labels; there is no hepatocyte on this arm. Same choice Phase 9 made
   in `code/h1_module_fits.py`.
2. `CURVE_TYPE = "B cells"` replaces the mouse `curves_for = [("Hepatocytes", j)]`. B cells
   are realised in all 7 sections at 13.5–21.3 % of cells; the frozen literal would have
   emitted no binned curves at all.
3. **N3-var / N4-var are the primary corrected null**, as the calibration study concluded:
   the tiled torus reaches **2.35× nominal type-I error** on an irregular window
   (0.1175 against 0.05) while the variance correction holds nominal (0.033–0.060)
   — `reports/CS_PHASE8_TORUS_VAR.md` §4, `results/phase3/var_sim_calibration.csv`. The tile
   variants are still computed and reported beside them.

**The sender calls.** All nine frozen calls plus the three merged-label calls were run at
stage `main` (12 calls × 7 sections = 210 jobs, 4,109 fits). `tierAmg_pNN` **is** PREREG
decision D-B's `tierA_merged_pNN` — identical rule, identical flags — and both spellings are
accepted by the shim; the short form is used in the output files only because
`sasp_phase3.Sec.sender_mask` dispatches on the `tierA_p` prefix and a name beginning
`tierA_` is one edit away from being read as the fine-label call.

**Declared, and raised rather than fixed** (`PREREG_PHASE8.md` §0.1 open item 2):
`tierA_merged_p95` is in neither `N7_CALLS` nor its nine-call extension `ALL9_CALLS`. Whether
it **replaces** or is **added to** the frozen N7 axis is a frozen-list question. **The frozen
lists were not edited.** `code/h1_run_phase10.ALL12` is a local extension in Phase-10 code,
and every stage that iterates a frozen call list (`run_phase5_*`, `run_phase3_poisson`) was
passed the call explicitly instead.

**Cache extension, declared.** The H1 cache carried no per-module Tier A columns and no
merged-label flags, so the pre-registered D1 sensitivity axis could not run.
`code/h1_cache_extend.py` copies flags **Phase 9 had already computed** into the `.npz` —
nothing is recomputed from expression, no threshold is touched, and every pre-existing key is
asserted byte-identical after the rewrite. One rounding discrepancy is declared: the senders
CSV stores the per-module score rounded to 5 dp, so a recomputation of the p95 flag from the
stored score differs from the committed unrounded flag on **at most 20 cells of 291,577
(0.0069 %)**, all exact ties at the percentile boundary. The committed column is used at p95;
p90 and p99 carry that ≤ 0.007 % tie noise. Per-section counts:
`results/phase10_h1/cache_extend_rounding_*.csv`.

### 1.2 The fitting window is still valid on H1, and the λ-grid floor is still violated

`results/phase10_h1/window.csv`. At the primary call the frozen **100 µm** window retains
**99.914–99.997 %** of receivers on every H1 section (`frac_gt_100` = 2.6 × 10⁻⁵ – 8.6 × 10⁻⁴;
99th percentile of distance-to-nearest-sender 43.8–61.7 µm), against 99.1 % on M1 — so the
window is, if anything, more comfortable here.
`python3 -c "import pandas as pd; w=pd.read_csv('results/phase10_h1/window.csv'); print(w[w.call=='tierAmg_p95'][['section','n','d_p99','frac_gt_100']].to_string(index=False))"`

**The λ-grid floor is not.** `LAM_LO_FLOOR = 7.0 µm` is defined as "at or below the median
nearest-neighbour distance of every section". On H1 the median NN distance is
**5.45–5.89 µm over all cells and 5.45–6.29 µm over QC-passed cells — below 7.0 µm in every
section.** PREREG §3.1 anticipates exactly this and **forbids patching it**. It is not
patched. **Consequence, carried into every H1 λ̂ below: the grid cannot resolve a length
constant shorter than 7 µm, i.e. shorter than ~1.2 median cell spacings, and any H1 λ̂ sitting
at the low bound must be read as "at or below the grid floor", never as an estimate.**
On H1 that affects **86 of 343 fits (25.1 %)** at the primary call.


---

## 2. The §17 two-arm comparison table

The paper's new centrepiece. Every cell is generated from
`results/phase10_h1/two_arm_table.csv` by `code/two_arm_table_md.py`, so no number here
can drift from the file it came from; that CSV carries the absolute source path and the
exact filter for every value. **Both H1 sender calls are shown side by side throughout**
(PI decision D-B: `tierAmg_p95` = `tierA_merged_p95` is primary, the frozen-literal
`tierA_p95` is the sensitivity); "—" in the last column means the two calls give the same
value. Quantities that are properties of the arm rather than of the call appear once.

| Quantity | M1 mouse liver (SBR fibrosis) | H1 human spleen — PRIMARY `tierAmg_p95` | H1 — frozen-literal `tierA_p95` |
|---|---|---|---|
| Platform / panel | Xenium Prime Mouse 5K + 100 custom; 5,097 panel genes | Xenium Prime Human 5K + 100 addon; 5,093 panel genes | — |
| Sections / donors | 11 / 11 mice; 6 admissible (§8 Test 3) | 7 / 7 donors, ages 17-59, 4M/3F; all 7 analysed | — |
| Cells, total | 1,826,893 | 2,207,593 | — |
| Cells, admissible sections | 1,031,880 | 1,962,278 | — |
| Median NN distance (um) | 6.74 - 10.61 | 5.45 - 5.89 (all cells); 5.45 - 6.29 (QC-passed) | — |
| Transcript assignment rate | 88.27 % (section 7259 only) | 92.41 - 94.17 % (3 sections) | — |
| Sender prevalence, primary call (%) | 4.041 - 4.480 | 3.837 - 4.553 | 3.379 - 4.554 |
| Fits / reportable | 315 / 153 (48.6 %) | 343 / 58 (16.9 %) | 343 / 61 (17.8 %) |
| Naive amplitude, median beta/sd(y) over reportable fits | 0.3288 | 0.1843 | 0.1807 |
| Controlled amplitude (N2+N5+N6), response-SD | 0.0288  IQR [-0.0066, 0.0845] | 0.0307  IQR [0.0028, 0.0666] | 0.0345  IQR [0.0006, 0.0707] |
| 80 %-power detectable bound, response-SD | 0.1833 | 0.1094 | 0.1083 |
| Controlled amplitude vs own bound | 0.0288 vs 0.1833 -> BELOW | 0.0307 vs 0.1094 -> BELOW | 0.0345 vs 0.1083 -> BELOW |
| SF, N2 matched decoy | 0.9516  IQR [0.9191, 0.9776] | 0.9786  IQR [0.9463, 0.9940] | 0.9739  IQR [0.9472, 0.9885] |
| SF, N5 alone | 0.1150  IQR [-0.0344, 0.2583] | 0.1264  IQR [-0.0245, 0.4016] | 0.1745  IQR [-0.0253, 0.4767] |
| SF, N6 alone | 0.4708  IQR [0.1829, 0.7416] | 0.7395  IQR [0.4594, 0.9387] | 0.8007  IQR [0.5037, 0.9616] |
| SF, anatomical covariate alone | 0.8429  IQR [0.5553, 0.9740] | 0.9307  IQR [0.7509, 0.9933] | 0.9526  IQR [0.8021, 1.0022] |
| SF, N5+N6 | 0.0844  IQR [-0.0126, 0.2333] | 0.1608  IQR [0.0047, 0.3869] | 0.1628  IQR [0.0070, 0.4766] |
| SF, N2+N5+N6 (PRIMARY OUTCOME) | 0.0885  IQR [-0.0166, 0.2338] | 0.1582  IQR [0.0074, 0.3837] | 0.1630  IQR [0.0040, 0.4759] |
| lambda-hat railed at a grid bound | 60.0 % (103 at the 7 um floor, 86 at the 50 um ceiling); median lam-hat 14.73 um | 70.0 % (86 at the 7 um floor, 154 at the 50 um ceiling); median lam-hat 38.16 um | 74.1 % (88 at the 7 um floor, 166 at the 50 um ceiling); median lam-hat 40.49 um |
| SF, N1 | 0.7074  IQR [0.4163, 0.8648]  n=153 | 0.3319  IQR [0.0474, 0.6583]  n=58 | 0.2623  IQR [0.0164, 0.6172]  n=61 |
| SF, N3 | 0.9989  IQR [0.9889, 1.0058]  n=153 | 0.9929  IQR [0.9770, 1.0020]  n=58 | 0.9975  IQR [0.9791, 1.0064]  n=61 |
| SF, N4 | 0.9472  IQR [0.8037, 1.0389]  n=153 | 0.9561  IQR [0.8119, 1.1317]  n=58 | 1.0181  IQR [0.8607, 1.2696]  n=61 |
| SF, C1 N3_tile (tile scope) | 0.9706  IQR [0.9061, 1.0089]  n=136 | 0.9507  IQR [0.8521, 1.0750]  n=58 | 0.9654  IQR [0.8490, 1.0769]  n=61 |
| SF, C1 N4_tile (tile scope) | 0.9243  IQR [0.8346, 1.0489]  n=136 | 0.9063  IQR [0.7836, 1.1417]  n=58 | 0.9101  IQR [0.7690, 1.1372]  n=61 |
| SF, N3_var | 0.9960  IQR [0.9754, 1.0074]  n=153 | 1.0000  IQR [0.9597, 1.0480]  n=58 | 1.0055  IQR [0.9530, 1.0501]  n=61 |
| SF, N4_var | 0.9851  IQR [0.9581, 1.0033]  n=153 | 0.9740  IQR [0.9455, 1.0118]  n=58 | 0.9788  IQR [0.9364, 1.0190]  n=61 |
| SF, C1 N3_orig | 1.0011  IQR [0.9915, 1.0079]  n=153 | 0.9921  IQR [0.9744, 1.0007]  n=58 | 0.9980  IQR [0.9890, 1.0065]  n=61 |
| SF, C1 N3_occ | 0.3024  IQR [0.0000, 0.7338]  n=153 | 0.5515  IQR [0.4041, 0.7390]  n=58 | 0.6080  IQR [0.4037, 0.7972]  n=61 |
| SF, C1 N3_swap | 0.6948  IQR [0.3920, 0.8721]  n=153 | 0.3374  IQR [0.0515, 0.6604]  n=58 | 0.2897  IQR [0.0039, 0.6609]  n=61 |
| SF, C1 N3_snap | 0.9934  IQR [0.9499, 1.0166]  n=153 | 1.0017  IQR [0.9857, 1.0185]  n=58 | 1.0050  IQR [0.9766, 1.0270]  n=61 |
| SF, C1 N4_orig | 0.9520  IQR [0.7965, 1.0480]  n=153 | 0.9611  IQR [0.8118, 1.1357]  n=58 | 1.0152  IQR [0.8607, 1.2708]  n=61 |
| SF, C1 N4_occ | 0.1833  IQR [0.0000, 0.5586]  n=153 | 0.4549  IQR [0.2784, 0.6329]  n=58 | 0.4815  IQR [0.3230, 0.6723]  n=61 |
| SF, C1 N4_swap | 0.9457  IQR [0.7861, 1.0345]  n=153 | 0.9867  IQR [0.8506, 1.1561]  n=58 | 1.0228  IQR [0.8807, 1.2231]  n=61 |
| Poisson identity (R3a): slope, r2 | -0.5249, r2 = 0.9843 | -0.5124, r2 = 0.9783 | — |
| R3c null destructiveness: senders keeping a real neighbour <=100 um | N3_occ 1.000; N3_occ15 0.969; N3_orig 0.772; N3_snap 1.000; N3_swap 1.000; N3_tile 1.000; N4_occ 1.000; N4_occ15 0.969; N4_orig 0.920; N4_swap 1.000; N4_tile 1.000 | N3_occ 0.998; N3_occ15 0.947; N3_orig 0.801; N3_snap 1.000; N3_swap 1.000; N3_tile 1.000; N4_occ 0.998; N4_occ15 0.957; N4_orig 0.923; N4_swap 1.000; N4_tile 1.000 | — |
| A7 PRIMARY: 40 negative control probes (naive) | -0.0225 [-0.0527, +0.0078] p = 0.129 | -0.0118 [-0.0340, +0.0103] p = 0.239 | — |
| A7 PRIMARY: 40 negative control probes (N6+N5) | +0.0015 [-0.0160, +0.0191] p = 0.849 | -0.0028 [-0.0263, +0.0206] p = 0.779 | — |
| A7 pooled negative-control features (naive) | -0.0744 [-0.1306, -0.0182] p = 0.0145 | -0.0337 [-0.0530, -0.0145] p = 0.00516 | — |
| A7 pooled negative-control features (N6+N5) | +0.0053 [-0.0162, +0.0268] p = 0.595 | -0.0051 [-0.0156, +0.0054] p = 0.28 | — |
| Composition-matched SF [comp] | 0.9837  95% CI [0.9730, 0.9942]  share removed 1.6 %  n_rep 165 | 0.9988  95% CI [0.9668, 1.0136]  share removed 0.2 %  n_rep 59 | 0.9872  95% CI [0.9465, 1.0168]  share removed 1.3 %  n_rep 61 |
| Composition-matched SF [full] | 0.9855  95% CI [0.9792, 0.9923]  share removed 1.4 %  n_rep 165 | 0.9445  95% CI [0.9255, 0.9636]  share removed 5.5 %  n_rep 59 | 0.9563  95% CI [0.9314, 0.9749]  share removed 3.9 %  n_rep 61 |
| Composition-matched SF [comp_adj] | 0.4989  95% CI [0.4215, 0.6063]  share removed 50.1 %  n_rep 33 | 0.1023  95% CI [-0.0469, 0.2202]  share removed 89.8 %  n_rep 12 | 0.1182  95% CI [-0.0356, 0.2200]  share removed 88.2 %  n_rep 12 |
| Composition-matched SF [type_adj] | 0.3414  95% CI [0.2356, 0.4016]  share removed 65.9 %  n_rep 33 | 0.2958  95% CI [0.1731, 0.4256]  share removed 70.4 %  n_rep 12 | 0.3349  95% CI [0.1649, 0.4372]  share removed 66.5 %  n_rep 12 |
| Composition-matched SF [typecomp_adj] (PRIMARY of the pair) | 0.1461  95% CI [0.0520, 0.2461]  share removed 85.4 %  n_rep 33 | 0.0698  95% CI [-0.1056, 0.1889]  share removed 93.0 %  n_rep 12 | 0.0821  95% CI [-0.0700, 0.1882]  share removed 91.8 %  n_rep 12 |
| Ripley's K at 50 um vs a within-type permutation null [tierA_p95] | 1.085 (range 1.009 - 1.160) over the 6 in-band sections; 1.136 (1.009 - 1.363) over all 11 | 1.012 (range 0.988 - 1.036) | — |
| Ripley's K at 50 um vs a within-type permutation null [cdkn1a_pos] | 1.319 (range 1.158 - 1.684) over the 6 in-band sections; 1.263 (1.032 - 1.684) over all 11 | 1.161 (range 1.118 - 1.218) | — |
| Ripley's K at 50 um vs a within-type permutation null [senepy_p95] | 1.593 (range 1.305 - 1.936) over the 6 in-band sections; 1.556 (1.272 - 1.936) over all 11 | 1.647 (range 1.295 - 2.415) | — |
| A5 matched-decoy balance, max \|SMD\| after matching | 0.0352 (100 % pass) | 0.0933 (35/35 matches pass \|SMD\| <= 0.1) | — |
| Caller agreement, depth- and type-matched, pooled [senepy_score vs cdkn1a_counts] | 1.211 (z 7.43) | 1.505 (z 18.70) | — |
| Caller agreement, depth- and type-matched, pooled [senepy_score vs deepscence_score] | 0.737 (z -15.08) | 1.093 (z 5.89) | — |
| Caller agreement, depth- and type-matched, pooled [tierA_score vs cdkn1a_counts] | 1.471 (z 19.45) | 1.081 (z 3.05) | — |
| Caller agreement, depth- and type-matched, pooled [tierA_score vs deepscence_score] | 1.288 (z 19.23) | 1.602 (z 38.86) | — |
| Caller agreement, depth- and type-matched, pooled [tierA_score vs senepy_score] | 0.972 (z -1.63) | 0.874 (z -7.96) | — |
| Caller agreement, depth- and type-matched, pooled [deepscence_score vs cdkn1a_counts] | 1.255 (z 10.53)  CIRCULAR | 6.436 (z 204.83)  CIRCULAR | — |
| Caller agreement, depth- and type-matched, pooled [tierA_score vs abs_deepscence_score] | — | 1.102 (z 5.67) | — |
| Caller depth loading, Spearman(score, transcript counts) [tierA_score] | 0.0561  (min -0.0609, max 0.1252, n=11) | -0.0227  (min -0.1738, max 0.0979, n=7) | — |
| Caller depth loading, Spearman(score, transcript counts) [cdkn1a_counts] | 0.1923  (min 0.0152, max 0.4814, n=11) | 0.0837  (min 0.0563, max 0.1653, n=7) | — |
| Caller depth loading, Spearman(score, transcript counts) [senepy_score] | 0.1676  (min 0.0968, max 0.5134, n=11) | 0.2339  (min 0.0405, max 0.3323, n=7) | — |
| Caller depth loading, Spearman(score, transcript counts) [deepscence_score] | 0.3846  (min -0.3500, max 0.5577, n=11) | 0.2963  (min 0.1822, max 0.3540, n=7) | — |
| Kernel family AIC win rate, controlled design | step 0.898; spline 0.054; gaussian 0.032; exponential 0.013; powerlaw 0.003 | step 0.956; spline 0.032; gaussian 0.012; exponential 0.000; powerlaw 0.000 | — |
| Superposition beats nearest-sender (controlled design) | 0.762 of fits; paired block-bootstrap win fraction median 0.730 | 0.746 of fits; paired block-bootstrap win fraction median 0.730 | — |
| Proximal vs downstream lambda ratio: identifiability | 6 receiver types; CI reaches a grid bound in 6, BOTH bounds in 2 | 6 receiver types; CI reaches a grid bound in 6, BOTH bounds in 3 | — |

*Every cell above is a row of `results/phase10_h1/two_arm_table.csv`, which carries the absolute source path and the exact filter beside each value. Regenerate with `python3 code/two_arm_table.py && python3 code/two_arm_table_md.py`. Brackets on surviving fractions and amplitudes are **inter-quartile ranges across the arm's reportable fits**, not confidence intervals (PREREG correction C-4); the only genuine CIs here are on the composition-matched SFs and the A7 section-clustered means.*

**The three rows to read first.** *Controlled amplitude vs own bound* is the §18 outcome-A
test and it is BELOW on both arms and both calls. *SF, N2+N5+N6* is the pre-registered
primary outcome. *SF, N1* is where the arms differ most, and §4 explains why.

**Two rows of the original §17 draft are deliberately absent.**
*(i)* "DeepScence sign vs own gene set: −0.350 sham / +0.318 SBR" is the **two-section**
framing PREREG §10.12 strikes. Measured over all 11 M1 sections the score's correlation
with transcript counts is **−0.3500 … +0.5577, median +0.3846, positive in 10 of 11**, and
the single negative section is **7250**, which the §8 Test-3 prevalence floor excludes
anyway. *(ii)* "Composition surrogate share 66–76 %" has no producer: the 66 % traces to a
superseded `summary_phase3.txt` line and the 76 % is a ratio whose numerator appears in no
results file. It is replaced by §7's `type_adj` / `typecomp_adj` rows, which are measured.
---

## 3. §8 — the pre-registered DeepScence prediction, against what H1 returned

**§8 exists to separate exactly two explanations of DeepScence's behaviour on M1: is it a
property of *our mouse adaptation* (ortholog remapping, `denoise=False`), or a property of
*the published tool*? On the axis that matters most for the tool's defence, it returned the
first answer, and that is the experiment working.**

### 3.1 P-ii is FALSIFIED, and it leads

**Prediction (registered before any H1 expression value was read):** the published `CDKN1A`
anchor is *weak, unstable or inverted* — depth-partialled fold-split sign stability < 0.90 —
in **≥ 1 of 7** sections. **Falsifier, also written in advance: all 7 sections at ≥ 0.90.**

**Measured: stability = 1.000 in all seven sections.** The anchor decides the score's polarity
the same way in **20 of 20 random folds in every section**, its depth-partialled correlation
with the score is **positive in all seven (+0.1911 … +0.2540)**, and `CDKN1A` ranks **1st to
7th of the 33 on-panel CoreScence genes** in every section.

```
python3 -c "import pandas as pd; print(pd.read_csv('results/phase9_h1/deepscence_anchor_h1.csv')[['section','rho_partial_cdkn1a','stab_cdkn1a','rho_partial_prolif','stab_prolif','rank_cdkn1a_in_core']].to_string(index=False))"
```

**The consequence, stated plainly because it runs against our own result: M1's polarity flip
(§13 P12, P13) is at least partly an artefact of our own ortholog remapping, not a defect in
the published tool.** On the native human panel the published anchor is the *better* anchor,
and the D3 alternative anchor we chose on M1 — the 8-gene proliferation set — is the one that
carries almost no signal here (depth-partialled ρ **+0.0097 … +0.0451**) and is itself
**unstable in one section (0.75, SPLN21)**. On M1 the published anchor's depth-partialled
correlation is **negative in 2 of 11 sections** (7248 −0.0118, 7435 −0.0241), ≈ 0 in a third
(7352 +0.0021), and its fold-split stability falls to **0.60** in 7352 and 0.85 in 7248
(`results/phase3/deepscence_anchor_decisions.csv`, columns `prho_ds_Cdkn1a`, `stab_Cdkn1a`).

**What this falsification does NOT reach, said in the same breath so the narrowing is visible
rather than implied.** All three of the following are measured on the **native** H1 run and
are **unaffected by ortholog remapping**:

1. **The 88 % CoreScence circularity.** 29 of 33 on-panel CoreScence genes sit in ≥ 1 frozen
   Tier B module = **0.8788**, on the native human panel, with `B_secondary_senescence` alone
   accounting for 18 of them.
   `python3 -c "import json; print(json.load(open('results/phase7_jobA/gate_result_human.json'))['corescence'])"`
2. **The seed instability.** *r* = 0.3719, top-5 % Jaccard = 0.2107 at full section size on
   the native panel (`results/phase9_h1/d2_stability.csv`).
3. **DeepScence × `CDKN1A`⁺ at 6.436 pooled, natively** — z = 204.8, above chance in 7 of 7
   sections, per-section range 3.459–10.115, against **1.255** pooled on the ortholog-remapped
   mouse arm. *The circularity is five times worse natively than the remapped run suggested.*
   `python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_agreement_pooled.csv'); print(d[(d.A=='deepscence_score')&(d.B=='cdkn1a_counts')].to_string(index=False))"`
   PREREG §10.7's prohibition stands: "1.51–2.85×" may never be quoted; the H1 figure is the
   pooled 6.436 with its range, and **this pair is excluded from every pooled cross-caller
   number**.

**And it does not vacate P12.** `LMNB1` remains excluded as the D3 primary anchor because it
is a member of `B_downstream_arrest` and `B_secondary_senescence` — a gene-set-membership fact
about our own modules, not a claim about anchor performance. On H1 `LMNB1` tracks `CDKN1A`
closely (`rho_partial_lmnb1` +0.2029 … +0.2449) and stays a secondary for that reason alone.

### 3.2 The other six predictions

| # | Prediction | Falsifier | H1 result | Verdict |
|---|---|---|---|---|
| **P-i** | score correlates **positively** with transcript counts in ≥ 5 of 7 sections, magnitude comparable to M1's | null or negative in ≥ 3 of 7 | **positive in 7 of 7**, ρ **+0.1822 … +0.3540, median +0.2963**; M1 over 11 sections **−0.3500 … +0.5577, median +0.3846** | **CONFIRMED**, at the bottom of M1's range |
| **P-ii** | anchor weak/unstable/inverted in ≥ 1 of 7 | all 7 at ≥ 0.90 | **1.000 in all 7** | **FALSIFIED** (§3.1) |
| **P-iii** | sign-invariant (\|score\|) depth- and type-matched agreement with Tier A above chance, pooled ratio **> 1.10** | pooled ≤ 1.05, or below chance | pooled **1.102**, z = 5.67, above chance in 5 of 7 | **CONFIRMED BY 0.002. Treat as MARGINAL** — see below |
| **P-iv** | CoreScence circularity ≈ 88 % natively | materially below 88 % | **29/33 = 0.8788**, native panel, no remapping | **CONFIRMED** |
| **P-v** | within cell type, bottom-of-depth selecting (Q5/Q1 < 1) in ≥ 5 of 7 | Q5/Q1 > 1 in ≥ 3 of 7 | **< 1 in 6 of 7**: 0.244, 0.404, 0.416, 0.457, 0.496, 0.778, and 1.170 (SPLN43) | **CONFIRMED** |
| **P-vi** | `denoise=True` **raises** the depth loading, as on 3 of 3 M1 sections | Δρ ≤ 0 in ≥ 5 of 7 | see §3.4 | **see §3.4** |
| **P-vii** | ≥ 1 of 3 `denoise=True` seeds gives a top-5 % set with Jaccard < 0.30 against the others | all three ≥ 0.60 | all three pairs **0.0444 / 0.0510 / 0.0482** | **CONFIRMED** |

**P-iii is confirmed by 0.002 and must never be quoted without that.** The registered
threshold is > 1.10 and the measurement is 1.102, on an estimator whose full-section
seed-to-seed top-5 % Jaccard is 0.211. **A margin of 0.002 on a call set that ~80 % turns over
between seeds is not a margin.** It is flagged as marginal wherever it appears, and it must be
re-evaluated on the five-seed consensus before it is reported as a verdict either way
(`PREREG_PHASE8.md` §0.1 open item 5).

### 3.3 Every H1 DeepScence magnitude above is on a NON-PRIMARY estimator

PI decision **D-A** makes the H1 DeepScence score a **five-seed consensus**. Every number in
§3.1–3.2 — including P-ii's stabilities, P-i's loadings, P-iii's 1.102 and D-C's 6.436 — was
computed at `random_state = 0`, i.e. on the estimator that D-A supersedes. **The directions
are seed-robust; the magnitudes are not.** Re-deriving the §10 statistics on a seed-1 score
for SPLN21 gave: ρ(score, counts) +0.312 → +0.231 (direction survives), within-type Q5/Q1
0.496 → 0.292 (direction survives, both ≪ 1). **Treat every magnitude in this section as
provisional pending the consensus.** §3.5 reports how far the consensus got.

### 3.4 P-vi is now FALSIFIED too, and it points the same way as P-ii

Phase 9 ran `denoise=True` on **one** of seven sections, so P-vi was *contradicted but not
falsified*: its falsifier is **Δρ ≤ 0 in ≥ 5 of 7 sections**. Phase 10 ran the remaining six.

**Declared approximation.** A full-section `denoise=True` run is a DCA autoencoder pass over
200–400k cells on CPU — TensorFlow 2.4 cannot use this box's CUDA-12 GPU — and was not
affordable beside the rest of Phase 10. The six new sections were therefore run on the **same
fixed 20,000-cell panel design the mouse D2 study used** (subsampling seed 12345, independent
of `--seed`, so every configuration sees identical cells), with the `denoise=False` companion
on the identical cells. **P-vi is evaluated at 7 of 7 on the 20,000-cell panel and at 1 of 7
at full section size.** On SPLN21, where both exist, they agree in sign and roughly in
magnitude (Δρ **−0.2104** full, **−0.3669** panel), which is why the panel is used — but it is
an approximation and is reported as one.

```
bash code/h1_run_dca_panel_all.sh        # 6 sections x {denoise=True, denoise=False}
python3 code/h1_d2_analyse.py            # -> results/phase9_h1/d2_depth.csv (rows APPENDED)
```

| section | ρ(score, counts) `denoise=False` | ρ `denoise=True` | **Δρ** | sender-set Jaccard, False vs True |
|---|---|---|---|---|
| SPLN07 | +0.2806 | **−0.1838** | **−0.4644** | 0.0627 |
| SPLN14 | +0.1695 | +0.1715 | **+0.0020** | 0.0283 |
| SPLN21 | +0.2697 | −0.0972 | **−0.3669** | 0.0304 |
| SPLN24 | +0.2237 | **−0.1564** | **−0.3801** | 0.0188 |
| SPLN30 | +0.1796 | +0.1349 | **−0.0447** | 0.0401 |
| SPLN43 | +0.0807 | **−0.2638** | **−0.3445** | 0.0695 |
| SPLN44 | +0.3709 | −0.0609 | **−0.4317** | 0.0152 |
| **M1 reference, 3 full sections** | 0.3891 / 0.3176 / 0.4096 | 0.6404 / 0.5314 / 0.5419 | **+0.13 … +0.25 (×1.32–1.67)** | 0.118 – 0.280 |

**Its seed-stability companion, carried with the number as PREREG §10.10 / P26 requires:**
every `denoise=True` value in the table above is a **single seed (`random_state = 0`)**, and on
the same 20,000-cell panel the three `denoise=True` seeds agree at top-5 % **Jaccard 0.0444 /
0.0510 / 0.0482** (P-vii) — i.e. the `denoise=True` call set is essentially unreproducible
across seeds, and only the **direction** of Δρ should be read from these rows, never the
magnitude. That is also why P-vi's falsifier was written on a sign and a count of sections
rather than on an effect size.

**Δρ ≤ 0 in 6 of 7 sections. The falsifier is ≥ 5 of 7. P-vi is FALSIFIED.** The one exception,
SPLN14, is +0.0020 — a ratio of 1.012, i.e. no change rather than the predicted rise. In four
of the seven sections `denoise=True` does not merely lower the depth loading, it **inverts its
sign**.

**P-vi's own pre-registered text states what its falsification means**: *"which would make the
M1 result an artefact of ortholog remapping"*. So **two of the seven §8 predictions are
falsified, and both falsifications point the same way** — part of what M1 attributed to the
published tool is attributable to our own mouse adaptation. That is the experiment doing its
job, and it is reported at the front rather than in a table.

**What survives regardless, and is now measured on 7 of 7 sections: denoising is not depth
normalisation on either arm.** It moves the caller's depth loading by ×1.32–1.67 **up** on
mouse and inverts or flattens it on human, and it changes **93.1 – 98.5 %** of the H1 sender
set in **every** section (Jaccard 0.0152 – 0.0695). PREREG **P29**'s conclusion — "whatever DCA
contributes here, it is not depth normalisation" — stands on both arms; **the *sign* of what it
contributes does not transfer, and the paper must say so rather than generalising the mouse
direction.**

**Declared file modification.** `results/phase9_h1/d2_depth.csv` was **appended to** by this
run: 6 rows added, every pre-existing row byte-identical (`git diff` shows 6 insertions and 0
deletions), and `d2_stability.csv` unchanged. It is the only Phase-9 output Phase 10 touched.

### 3.5 The five-seed DeepScence consensus (PI decision D-A) — specification and status

**The specification was written down and committed BEFORE any consensus value was read**, which
is what D-A requires of a producer-level choice (`code/h1_deepscence_consensus.py`, committed at
`7edf07c`, before the first seeded run finished). In order:

0. **SIGN ALIGNMENT FIRST.** DeepScence fixes its bottleneck sign by correlating with `CDKN1A`,
   and that polarity is known to flip — it inverted between the two surgical arms on M1
   (−0.350 sham / +0.318 SBR). **Z-scoring does not fix a sign flip**, and a per-cell median
   over a mixture of two polarities is meaningless. Each seed's polarity is measured against
   the published anchor — the depth-partialled Spearman of its score with `CDKN1A` counts,
   using `deepscence_reanchor.partial_spearman` **imported** so the definition cannot drift —
   and a seed whose anchor correlation is negative while the majority are positive is
   multiplied by −1. **The number of seeds flipped is reported, per section, whatever it is.**
   The raw pairwise Pearson matrix is written beside it as an independent check. Phase 9 found
   the anchor stable in 20/20 folds in all seven sections, so **zero flips is the expectation —
   it is verified, not assumed.**
1. z-score each aligned seed score **within its section**.
2. **consensus = the per-cell MEDIAN** across seeds. Median, not mean, for robustness to one
   divergent seed — which is not hypothetical: on M1 one of three `denoise=True` seeds gave a
   top-5 % set **perfectly disjoint** from the other two (Jaccard 0.000).
3. threshold as normal: the frozen within-cell-type percentile rule at p90 / p95 / p99.

**Dispersion is reported twice, because the two diverge.** On M1 the same comparison gave
*r* = 0.99553 and Jaccard = 0.7606: **a score that looks stable can sit on a call-set that is
not, and the call-set is what downstream analysis consumes.**
- **score**: the per-cell **interquartile range** across seeds, summarised per section;
- **call-set**: the **mean pairwise top-5 % Jaccard** across seeds, **with its range**.

**Seeds: the five already-frozen composition-matched values 20260901–05** (PREREG §3.8). No new
seed value is introduced. Configuration otherwise frozen: `denoise=False`, published `CDKN1A`
anchor, ≥ 20 counts/cell, native human panel, **full sections, no subsampling**.


*Status at the time of writing — see §12 for the final count.* The five-seed panel is
**35 full-section runs** (7 sections × 5 seeds) and only **two fit concurrently** in the
57.7 GB cgroup. Measured rate under the Phase-10 load: **25–65 min per run**, i.e. roughly
**1.5 h per completed seed across all 7 sections**. The seed count was **not** silently
reduced; the producer pools over whatever exists, names the estimator `consensus_k<n>`, writes
`n_seeds` into every output row, and **refuses to pool a section with fewer than 3 seeds**.
The final per-section coverage is in `results/phase10_h1/deepscence_consensus_coverage.csv`.

**Until the consensus exists, no H1 DeepScence magnitude in this report is on the primary
estimator**, and every one of them is marked as such (§3.3). The directions — P-i's positive
depth loading, P-ii's stable anchor, P-v's bottom-of-depth selection, P-vi's inverted denoise
effect — are seed-robust and are not affected.

---

## 4. Where the two arms differ — the evidence for outcome D

**The null result replicates. The confound structure does not.** Seven differences, each read
from a file. None of them may be attributed to species or to tissue: mouse liver against human
spleen confounds the two by design (PREREG §10.5).

| # | Quantity | M1 mouse liver | H1 human spleen | Source |
|---|---|---|---|---|
| 1 | **Reportable fits** (naive β > 0 **and** block-bootstrap CI excludes 0) | **153 / 315 = 48.6 %** | **58 / 343 = 16.9 %** (merged call) · 61 / 343 = 17.8 % (frozen-literal) | `main_fits.csv`, both arms |
| 2 | **SF under N2+N5+N6** (the primary outcome) | **0.0885** IQR [−0.0166, +0.2338] | **0.1582** IQR [+0.0074, +0.3837] | `sf_summary.csv` / `h1_headlines_*.json` |
| 3 | **SF under N1**, the within-type sender-label permutation | **0.7074** — a random relabelling of senders within cell type reproduces **29 %** of the naive amplitude | **0.3319** — the same relabelling reproduces **67 %** of it | `perm_nulls.csv`, both arms |
| 4 | **Ripley's K at 50 µm**, `tierA_p95` senders vs a within-type permutation null | **1.085** (1.009–1.160) over the 6 in-band sections, **1.136** (1.009–1.363) over all 11 | **1.012** (0.988–1.036); **3 of 7 sections below 1.0** | `results/phase3/ripley.csv`, `results/phase9_h1/a4_ripley.csv` |
| 5 | **λ̂ geometry** | 60.0 % railed: **32.7 % at the 7 µm floor**, 27.3 % at the 50 µm ceiling; median λ̂ **14.73 µm** | 70.0 % railed: 25.1 % at the floor, **44.9 % at the ceiling**; median λ̂ **38.16 µm** | `main_fits.csv`, both arms |
| 6 | **Local cell density**: median real cells within the 100 µm window of a sender | **140** | **476** — 3.4× denser | `null_destructiveness.csv`, `real_median_nbrs`, both arms |
| 7 | **Sequencing depth**: median transcripts per cell, per section | **446 – 968** (in-band 6) | **42 – 273** — 3 to 10× shallower, and spanning **6.5×** within the arm | `results/section_qc_sender_summary.csv`; `results/phase9_h1/a1_sections.csv` |

**Rows 3 and 4 are one mechanism, and it is the most consequential difference.**
*(Read the surviving fraction as PREREG §3.5 defines it for a perturbation null:
SF = (β̂_obs − mean β̂_null)/β̂_obs, so **SF ≈ 1 means the null does NOT reproduce the observed
amplitude** and SF ≈ 0 means it reproduces it entirely.)* On H1 the
Tier A sender set is **not spatially aggregated beyond a random draw from the same cells of
the same types** (Ripley 1.012, and below 1.0 in 3 of 7 sections). A within-cell-type
label permutation therefore produces a distance-to-nearest-sender field that is nearly the
observed one, and **two-thirds of H1's naive amplitude is reproduced by randomly relabelling
which cells are senders — against 29 % on M1.** **The H1 arm is a weaker starting point for a
distance-to-sender analysis than M1 was, and the ~0.16 residual it carries after the full
nuisance design must be read against that**, not as a stronger signal. It is the same amplitude structure with less of it attributable to sender geometry.

**Row 5 is a real inversion.** M1's fitted length constants pile up at the *short* end of the
grid; H1's pile up at the *long* end. This is not a resolution artefact in the usual direction:
H1's median NN distance is **below** the grid floor (§1.2), so if anything H1 has *more* room
at the short end than M1. It is the opposite behaviour.

**Row 7 is why the arms cannot be compared on any depth-sensitive quantity without the N5
block.** `CDKN1A`⁺ prevalence on H1 tracks **sequencing depth, not age**:
Spearman(median depth, prevalence) = **+0.857, p = 0.0137** across the 7 donors;
Spearman(age, prevalence) = **−0.036, p = 0.939** — reproducing M1's ρ = **+0.943**
(SBR, n = 6, `results/section_qc_sender_summary.csv`) in a normal human ageing cohort where
the crude caller's prevalence varies 5.6× across donors (1.11 % – 6.21 %).

```
python3 -c "
import pandas as pd; from scipy.stats import spearmanr
d=pd.read_csv('results/phase10_h1/poisson_density.csv'); a=pd.read_csv('results/phase9_h1/a1_sections.csv')
c=d[d.call=='cdkn1a_pos'][['section','median_depth','sender_prevalence']].merge(a[['section','age']],on='section')
print(spearmanr(c.median_depth,c.sender_prevalence)); print(spearmanr(c.age,c.sender_prevalence))"
```

**PREREG §10.3 forbids any age-stratified claim on H1; this is a negative statement about the
caller, not an age result.**

### 4.1 A7 — the technical null is not flat in the same way on the two arms

The pre-registered primary technical null — the **40 negative control probes** — **is flat on
both arms**, naive and conditioned, so **A7 passes on both**:

| response | design | M1 | H1 |
|---|---|---|---|
| **40 negative control probes (PRE-REGISTERED PRIMARY)** | naive | −0.0225 [−0.0527, +0.0078] p = 0.129 | **−0.0118 [−0.0340, +0.0103] p = 0.239** |
| | +N6+N5 | +0.0015 [−0.0160, +0.0191] p = 0.849 | **−0.0028 [−0.0263, +0.0206] p = 0.779** |
| **pooled negative-control features** | naive | **−0.0744 [−0.1306, −0.0182] p = 0.0145** | **−0.0337 [−0.0530, −0.0145] p = 0.0052** |
| | N2 matched decoy | −0.0642 (86 % undiminished) | −0.0331 (98 % undiminished) | 
| | +N6+N5 | +0.0053 [−0.0162, +0.0268] p = 0.595 | −0.0051 [−0.0156, +0.0054] p = 0.280 |

```
python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); print(d[d.design.isin(['base','n2','n6n5'])].to_string())"
python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/a7_summary.csv'); print(d[d.design.isin(['base','n2','n6n5'])].to_string())"
```

**Two things replicate and one differs.**
- **Replicates:** the raw assay is **not** flat on its pooled control features on either arm,
  and a **matched-decoy contrast is not a substitute for the technical covariate block** —
  N2 leaves **86 %** of the gradient on M1 and **98 %** on H1, while N5 removes it entirely on
  both. This is an independent replication of P2 and P24 in a second species, a second tissue
  and a second laboratory, and it is what makes PREREG §10.1's prohibition — no naive and no
  N2-only kernel may be reported as a distance effect — binding on both arms.
- **Differs:** *which family carries the non-flatness*. On M1 it is the **609 negative-control
  codewords** (~73 % of control counts). On H1 the codewords are **62× sparser**
  (0.00069 counts/cell, 0.068 % of cells non-zero) and are flat, and the non-flatness is
  carried by the **21 genomic controls**, which are 90 % of H1's control counts. **The H1
  codeword row is a low-power null, not a demonstration of flatness**, and 10 of 98 codeword
  fits had `sd_y` exactly 0 and were dropped and counted rather than allowed to explode a
  ratio.

### 4.2 The anatomical covariate is the weakest input on H1, and it does not matter to the answer

A6's red/white-pulp axis occupies the frozen N5 block's `zonation` slot. **Its independent
validation does not hold: the follicle-distance check has the wrong sign in 2 of 7 sections**
(SPLN14 −0.013, SPLN30 −0.052) and never exceeds +0.172 raw. The axis does have real spatial
structure (Moran's I 0.221–0.492 on the 20-NN graph in every section), and M1's own equivalent
check ranged −0.085 to +0.343 over 11 sections — so H1 is *comparable* to M1 rather than
worse, and **neither should be presented as a validation.**

§15 requires the result both ways. `code/h1_a6_sensitivity.py` adds the one design the frozen
nesting cannot express — the primary design with the anatomical term removed and everything
else kept — using the same `BlockProfiler`, the same cells and the same λ index:

| design | H1 primary call (`tierAmg_p95`) | H1 frozen-literal (`tierA_p95`) |
|---|---|---|
| N6 only (no anatomy, no other N5) | 0.7395 [0.4594, 0.9387] | 0.8007 [0.5037, 0.9616] |
| N6 + the anatomical term only | 0.6947 [0.3820, 0.9117] | 0.7438 [0.4135, 0.9345] |
| N2+N5+N6 **without** the anatomical term | **0.1563** [0.0065, 0.4091] | **0.1576** [0.0166, 0.4685] |
| N2+N5+N6 **with** it (frozen, PRIMARY) | **0.1582** [0.0074, 0.3837] | **0.1630** [0.0040, 0.4759] |
| the anatomical term alone | 0.9307 [0.7509, 0.9933] | 0.9526 [0.8021, 1.0022] |

```
python3 code/h1_a6_sensitivity.py --calls tierAmg_p95,tierA_p95 --n-jobs 4
```

**The weakest input in the design moves the primary outcome by 0.002–0.005 of a surviving
fraction, and the controlled amplitude by 0.0005–0.002 response-SD.** The A6 limitation is
real and is declared, but it does not reach the conclusion. Stated in that direction because
it is the direction the evidence goes, not because it is convenient: had it gone the other
way the H1 primary outcome would have had to be withdrawn.

---

## 5. The ortholog-intersected panel — every cross-arm number, twice

PREREG §9 item 4 requires every cross-arm number on **both** each arm's full panel and the
**2,425-gene ortholog-intersected panel**. A8 found the restriction is not cosmetic:
**it changes 45 % of the H1 sender set at p95 while the underlying score correlates at
ρ = 0.88.** Phase 10 therefore did the restriction as a **refit on both arms**, not as a
diagnostic.

**Method.** Tier A and the seven Tier B modules were rescored with the AnnData `var` restricted
to the intersected symbols **before** normalisation, so `score_genes` draws its control genes
from the restricted universe — which is the point of the restriction. New caches
(`cache3_h1_iso`, `cache3_m1_iso`) are copies of the frozen caches with only those keys
replaced; **52 of 67 keys (H1) and 51 of 63 (M1) are asserted byte-identical** after the
rewrite. `run_phase3_nulls._section_job` is then run unchanged, with the seeds taken from the
full-panel run's own section/call lists so **every N2 match draw is identical to the
full-panel fit it is compared against** (recorded per fit in `seeds_used.csv`).

```
python3 code/iso_panel.py                       # asserts 5097 / 5093 / 4845 / 2435 / 2425
python3 code/h1_callers_iso.py ; python3 code/h1_prep_cache_iso.py
python3 code/m1_callers_iso.py ; python3 code/m1_prep_cache_iso.py
python3 code/run_phase10_iso.py --arm h1 --n-jobs 4 ; python3 code/run_phase10_iso.py --arm m1 --n-jobs 4
python3 code/iso_compare.py                     # -> results/phase10_h1_iso/iso_vs_full_headlines.csv
```

### 5.1 What the restriction does to the sender call, on BOTH arms

`results/phase10_h1_iso/iso_sender_shift.csv`, `results/phase10_m1_iso/iso_sender_shift.csv`.
A8 measured this on H1 only; here it is extended to M1.

| arm | Spearman(full score, intersected score) | top-5 % sender-set Jaccard | sender set that **turns over** |
|---|---|---|---|
| **M1** (6 in-band, `tierA_p95`) | 0.8756 – 0.8882 | **0.4586 – 0.4903** | **51 – 54 %** |
| **H1** (7 sections, `tierAmg_p95`) | 0.8480 – 0.9340 | **0.5192 – 0.5722** | **43 – 48 %** |

**Half the sender set is a different set of cells on the intersected panel, on both arms,
while the score correlates at ρ ≈ 0.88.** Any cross-arm sender-based comparison must carry
this number.

### 5.2 The primary outcome on both panels, both arms

`results/phase10_h1_iso/iso_vs_full_headlines.csv`. The three `full` rows reproduce
`h1_headlines_*.json` and `m1_headlines.py` field for field.

| arm | call | panel | reportable | naive amp | **controlled amp** | **own 80 % bound** | **SF N2+N5+N6** | IQR | λ̂ railed | median λ̂ |
|---|---|---|---|---|---|---|---|---|---|---|
| M1 | `tierA_p95` | full | 153 / 315 | 0.3288 | **0.0288** | 0.1833 | **0.0885** | [−0.0166, +0.2338] | 60.0 % | 14.73 µm |
| M1 | `tierA_p95` | **intersected** | 133 / 315 | 0.3010 | **0.0285** | 0.1987 | **0.1054** | [−0.0470, +0.3048] | 71.1 % | 15.80 µm |
| H1 | `tierAmg_p95` | full | 58 / 343 | 0.1843 | **0.0307** | 0.1094 | **0.1582** | [+0.0074, +0.3837] | 70.0 % | 38.16 µm |
| H1 | `tierAmg_p95` | **intersected** | 56 / 343 | 0.1864 | **0.0346** | 0.1136 | **0.1785** | [+0.0063, +0.3623] | 72.0 % | 33.37 µm |
| H1 | `tierA_p95` | full | 61 / 343 | 0.1807 | **0.0345** | 0.1083 | **0.1630** | [+0.0040, +0.4759] | 74.1 % | 40.49 µm |
| H1 | `tierA_p95` | **intersected** | 58 / 343 | 0.1844 | **0.0337** | 0.1170 | **0.2236** | [**−0.0045**, +0.4372] | 74.9 % | 38.11 µm |

**The conclusion is unchanged on both panels on both arms: the controlled amplitude is 3.1 –
7.0 × below each arm's own detectable bound, and 0 of 7 Tier B modules crosses it anywhere.**

**One thing the intersected panel does change, and it is against the tidy version of §0.1.**
R1's "IQR includes 0" clause is **panel- and call-dependent on H1**: it fails on the full
panel at both calls, fails on the intersected panel at the merged call, and **passes on the
intersected panel at the frozen-literal call** (IQR [−0.0045, +0.4372]). M1's IQR includes 0
on both panels. **A criterion whose verdict flips with a gene-panel restriction that leaves
the score at ρ = 0.88 is a criterion sitting on noise**, which is the same conclusion §0.1
reaches from the size of the margin. It is reported, not resolved.

**Declared limitation of the ISO caches.** `senepy_score`, `tierApm__*` and `flag_pm_*` are
carried through **unchanged from the full-panel cache**, because SenePy and the seven
per-module Tier A sets were deliberately not rescored. Those keys are full-panel values
sitting in an intersected-panel cache: **the `senepy_*` and `tierApm_*` calls must not be run
against `cache3_h1_iso` or `cache3_m1_iso`.** Nothing mechanically blocks it; it is declared
in both `*_prep_cache_iso.py` docstrings and in a `STALE_FULL_PANEL` list.

**Also declared:** `sasp_real.load_expression` returns 5,106 Gene Expression names on the mouse
arm, not the 5,097 of the metadata-derived panel used by `h1_a8_crossarm.panels()`. All 2,435
intersected mouse symbols are present in the loader's `var` on every section (asserted), so
the restriction is exact either way, but the two panel definitions are not the same object.

---

## 6. Kernel families, superposition vs nearest, proximal vs downstream

`code/h1_run_phase5.py` runs `run_phase5_kernels` and `run_phase5_super` on H1.
**A hazard was found and handled rather than tripped over:** none of the three Phase-5
producers respects `sasp_phase3.RESULTS` — all take their output directory from
`phase5_common.RES5 = "/workspace/results/phase5"`, the mouse directory, and
`run_phase5_kernels.py` has no `--tag`. **Running them unwrapped on H1 would have silently
overwritten the committed mouse Phase-5 CSVs.** The wrapper rebinds the output directory and
the already-captured `RES` in each runner module, passes `--call` explicitly (their default is
`run_phase3_nulls.PRIMARY_CALL`, the fine-label call, not H1's primary), and rebinds the arm
inside the loky workers.

```
python3 -u code/h1_run_phase5.py --which kernels --stage section  --call tierAmg_p95 --n-jobs 4
python3 -u code/h1_run_phase5.py --which kernels --stage proxdown --call tierAmg_p95 --n-jobs 4
python3 -u code/h1_run_phase5.py --which super   --stage section  --call tierAmg_p95 --n-jobs 4
```

### 6.1 Kernel families — the shape replicates

Five families on the same cells under the same design and the same profiler
(`exponential`, `gaussian`, `powerlaw`, `step`, `spline`), compared by AIC.
Fraction of fits each family wins, and fraction that beat a covariates-only model:

| design | family | M1 wins | H1 wins (`tierAmg_p95`) | H1 wins (`tierA_p95`) | M1 beats cov-only | H1 beats cov-only |
|---|---|---|---|---|---|---|
| **ctrl** | **step** | **0.898** | **0.956** | **0.948** | 0.511 | **0.536** |
| ctrl | spline | 0.054 | 0.032 | 0.029 | 0.152 | 0.105 |
| ctrl | gaussian | 0.032 | 0.012 | 0.017 | 0.241 | 0.155 |
| ctrl | exponential | 0.013 | 0.000 | 0.006 | 0.197 | 0.122 |
| ctrl | powerlaw | 0.003 | 0.000 | 0.000 | 0.130 | 0.073 |
| naive | step | 0.508 | 0.676 | 0.685 | 0.898 | 0.749 |

**Under the controlled design a step function beats the exponential kernel the whole model is
written in, on ~90–96 % of fits, on both arms.** That is a geometric statement, not a
biological one, and it replicates across species, tissue and laboratory. It also says the same
thing the railing rate says: what the data supports is a *near/far contrast*, not a length
constant.

**Held out across sections, no family beats a covariates-only model on either arm.**
Leave-one-section-out log-likelihood per cell against covariates-only, controlled design
(`kernel_heldout.csv`, 294 H1 folds and 252 M1 folds per family):

| family | M1 median `dll_vs_cov` | M1 folds > 0 | H1 median | H1 folds > 0 |
|---|---|---|---|---|
| gaussian | +0.0913 | **0.544** | −0.0185 | 0.486 |
| exponential | +0.0583 | 0.520 | +0.0106 | **0.514** |
| powerlaw | +0.0462 | 0.540 | +0.0006 | 0.500 |
| step | −0.0604 | 0.480 | −0.1703 | 0.429 |
| spline | −0.3864 | 0.437 | −0.2978 | 0.422 |

**The best any family manages is 54.4 % of folds on M1 and 51.4 % on H1 — a coin flip.**
The family that wins *within* a section (step) is among the worst *across* sections on both
arms, which is the signature of in-sample selection, not of a transferable kernel.

### 6.2 Superposition vs nearest — replicates to three decimal places

Nearest = `exp(−d_i/λ)` on the distance to the single nearest sender; superposition =
`Σ_j exp(−‖x_i − x_j‖/λ)` over all senders. Same λ grid, same design, same cells, same
parameter count, so `ΔAIC = n·log(RSS_sup/RSS_near)` exactly.

| design | quantity | M1 | H1 (`tierAmg_p95`) | H1 (`tierA_p95`) |
|---|---|---|---|---|
| naive | superposition wins (fraction of fits) | 0.844 | 0.834 | 0.843 |
| naive | ΔAIC per 1,000 cells, median | −0.765 | −0.481 | −0.326 |
| **ctrl** | **superposition wins** | **0.762** | **0.746** | **0.732** |
| **ctrl** | **paired block-bootstrap win fraction, median** | **0.730** | **0.730** | **0.730** |
| ctrl+N2 | superposition wins | 0.752 | 0.776 | 0.778 |

**Superposition beats nearest-sender on three-quarters of controlled fits on both arms, and
the paired block-bootstrap win fraction is 0.730 on both to three decimals.** The direction
replicates; the strength of preference (ΔAIC per 1,000 cells) is 1.6–2.3× weaker on H1.
**Neither basis beats a covariates-only model on either arm**, which is why this is a
statement about model form and not a rescue of the kernel.

### 6.3 Proximal vs downstream — unidentified on both arms

`results/phase10_h1/proximal_vs_downstream.csv`, `results/phase5/proximal_vs_downstream.csv`.
The λ of `tnfa_nfkb_proximal` against the λ of `downstream_arrest`, fitted against the **same**
sender call so a difference cannot be a sender-definition difference. The grid ratio bounds are
[0.140, 7.143].

| arm | receiver types fitted | CI reaching a grid bound | CI reaching **both** bounds |
|---|---|---|---|
| M1 | 6 | 6 / 6 | **2 / 6** |
| H1 (`tierAmg_p95`) | 6 | 6 / 6 | **3 / 6** |

On H1 the point ratios are 7.143 (B cells, Endothelial, Plasma cells — all railed), 1.398
(Stromal), 1.230 (T/NK) and 1.000 (Mono/Mac/DC); every interval spans at least one grid bound
and three span both. **The ordering of proximal and downstream length constants is not
identified on either arm**, and the fitted λ̂s themselves rail at 7 µm or 50 µm in almost every
cell. No proximal-vs-downstream claim is supportable from either arm.

---

## 7. The composition-matched rerun, five frozen seeds, both arms

PI decision **D15**: the **covariate-adjusted counterpart is PRIMARY** for any claim about how
much of the gradient is composition, and the **matched-decoy protocol is reported alongside,
every time**, because §15 specifies it and **its inertness is itself the finding** (P23/P25).
PREREG §10.8 forbids quoting the matched number alone.

Seeds: the five frozen values **20260901–20260905**. The M1 arm was run in Phase 8; the H1 arm
is roadmap item 10.2 and is run here, through the same frozen producer — `_job` is
`run_phase8_compmatch._job`, unchanged; `code/h1_run_compmatch10.py` exists only to rebind the
arm inside the loky workers. `ARMS['h1']` stays `frozen=True`, so `SASP_H1_UNFROZEN=1` is
still required and the gate is exercised rather than bypassed.

```
SASP_H1_UNFROZEN=1 python3 -u code/h1_run_compmatch10.py \
    --calls tierAmg_p95,tierA_p95,tierApm_p95 --n-jobs 6 --out-tag _h1
```

### 7.1 Pooled scope, the fit the §17 composition row is about

| variant | what it is | **M1** SF (share removed) | **H1 `tierAmg_p95`** (PRIMARY call) | **H1 `tierA_p95`** | **H1 `tierApm_p95`** (D1) |
|---|---|---|---|---|---|
| **`comp`** | the §15 matched-decoy protocol | 0.9837 [0.9730, 0.9942] — **1.6 %** | **0.9988 [0.9668, 1.0136] — 0.19 %** | 0.9872 [0.9465, 1.0168] — 1.3 % | 0.9926 [0.9489, 1.0162] — 0.7 % |
| `full` | the published N2 matching set | 0.9855 [0.9792, 0.9923] — 1.4 % | 0.9445 [0.9255, 0.9636] — 5.5 % | 0.9563 [0.9314, 0.9749] — 3.9 % | 0.9477 [0.9249, 0.9658] — 5.2 % |
| `comp_adj` | the same 20-NN composition vector as **covariates** | 0.4989 [0.4215, 0.6063] — 50.1 % | 0.1023 [−0.0469, 0.2202] — **89.8 %** | 0.1182 [−0.0356, 0.2200] — 88.2 % | 0.1138 [−0.0741, 0.2104] — 88.6 % |
| `type_adj` | receiver cell-type intercepts | 0.3414 [0.2356, 0.4016] — 65.9 % | 0.2958 [0.1731, 0.4256] — 70.4 % | 0.3349 [0.1649, 0.4372] — 66.5 % | 0.3120 [0.0685, 0.4279] — 68.8 % |
| **`typecomp_adj`** | **both — PRIMARY** | **0.1461 [0.0520, 0.2461] — 85.4 %** | **0.0698 [−0.1056, 0.1889] — 93.0 %** | **0.0821 [−0.0700, 0.1882] — 91.8 %** | **0.0647 [−0.1538, 0.1495] — 93.5 %** |

n reportable: M1 165 (`comp`/`full`, 5 seeds) and 33 (adjusted, seed-free); H1 59–65 and 12–13.
Match rate 1.000 and the balance gate |SMD| ≤ 0.1 passes in **100 %** of matches on H1 (max
|SMD| after matching 0.0369 at the primary call). Between-seed spread of the pooled matched SF
on H1: **sd 4.3 × 10⁻⁴, range 0.99790 – 0.99900** — the same order as M1's 1.2 × 10⁻⁴.

### 7.2 What replicates

**The inertness of the matched-decoy protocol replicates, and is more extreme on H1.**
On M1 the same variables remove **1.6 %** as a matching set and **85.4 %** as covariates. On H1
they remove **0.19 %** as a matching set and **93.0 %** as covariates — a factor of **490**
between two analyses of the same variables on the same fits. **This is the cross-arm
replication of P23/P25, and it is the strongest single argument in the report that a
matched-decoy design is not a substitute for conditioning.** It is also the second independent
demonstration of that on H1: A7 found the matched-decoy contrast leaves **98 %** of the
technical gradient in the negative-control response (§4.1) while N5 removes it entirely.

**Composition is the largest single nuisance on both arms, and it is larger on H1.**
Receiver cell type alone removes 65.9 % (M1) and 70.4 % (H1); the 20-NN composition vector
alone removes 50.1 % (M1) and 89.8 % (H1); both together remove 85.4 % and 93.0 %.

**The by-cell-type scope tells the same story**, and is reported because the pooled scope is a
single fit: `comp` removes 3.6 % on H1 (289 reportable fits) against 3.5 % on M1 (748), while
`comp_adj` removes 52.6 % on H1 (57 fits) against 50.2 % on M1 (150).

---

## 8. The constraints, reported rather than worked around

Each of these was given as something to respect, not solve. None was tuned around, and each
is stated where it bites.

| # | Constraint | How it is handled here |
|---|---|---|
| **1** | **The frozen PRIMARY DeepScence configuration is not seed-reproducible on H1** — `random_state=1` reproduces the committed `random_state=0` at *r* = **0.372**, top-5 % Jaccard **0.211**, against M1's floor of 0.9955 / 0.761 | **No DeepScence magnitude is quoted on H1 without this caveat.** §3.3 states it once for the whole section, Figure 6 states it on its own face, and PI decision **D-A** replaces the estimator with a five-seed consensus (§3.5). Directions survive the seed; magnitudes do not |
| **2** | **A live, unresolved frozen-pipeline interaction**: the Tier A call thresholds within *fine* cell types while the estimator stratifies receivers on *merged* labels, so **0–13.6 % of cells per section are eligible but uncallable**, producing **0 senders in 24,815 cells** for T/NK in SPLN43 | **Both calls are computed and reported side by side everywhere** — every table, every figure, the §17 table and the ISO refit. Nothing is tuned: `tierAmg_pNN` is the *identical percentile rule at the other label family*. The uncallable fraction is quantified per section below. The two calls agree on every conclusion |
| **3** | **A6's covariate does not validate cleanly** — the follicle-distance check has the wrong sign in 2 of 7 sections | **Reported with and without the anatomical term**, as §15 requires (§4.2). It moves the primary outcome by 0.002–0.005 SF |
| **4** | **Marginal-zone B cells are never realised** | **P7's exploratory MZ claims cannot be supported at all, and no MZ result is reported in either direction.** PI decision **D-D** retires P7 as **UNANSWERABLE** — never testable on this panel, which is a different statement from tested and null, and is recorded as the different statement it is. The label is one of six never realised in any of the 7 sections |
| **5** | **H1's median NN distance (5.45–6.29 µm) is below the frozen 7.0 µm λ-grid floor in every section** | **Declared, not patched** (§1.2). PREREG §3.1 pre-registers this contingency and forbids patching it. Consequence carried into every λ̂: 86 of 343 fits (25.1 %) rail at a floor the tissue is finer than |
| **6** | **SenePy has no spleen hub**; caller 2 is not the same estimator across arms | Already a declared deviation. Every H1 SenePy score is a **cross-tissue surrogate**; 7 of 23 labels get no hub in any tissue; 92.5–100 % of eligible cells are scored, with both denominators carried. Its within-type depth enrichment on H1 is **Q5/Q1 = 28.5–224.7**, the most extreme of the four callers on either arm, and no SenePy-called H1 quantity is reported without it |

### 8.1 The eligible-but-uncallable fraction, per section (constraint 2)

| section | fine-`Unknown` but merged-assigned | QC-passed cells | % eligible but uncallable |
|---|---|---|---|
| SPLN07 | 0 | 227,360 | 0.00 |
| SPLN14 | 38,701 | 283,628 | **13.64** |
| SPLN21 | 20,155 | 196,142 | **10.28** |
| SPLN24 | 7,386 | 393,202 | 1.88 |
| SPLN30 | 0 | 291,577 | 0.00 |
| SPLN43 | 24,815 | 270,472 | **9.17** |
| SPLN44 | 27,588 | 299,897 | **9.20** |

```
python3 -c "
import pandas as pd, json
BAD={'Low_quality','Unknown','unknown'}
for s in ['SPLN07','SPLN14','SPLN21','SPLN24','SPLN30','SPLN43','SPLN44']:
    d=pd.read_csv(f'data/processed_h1/celltypes_h1_{s}.csv')
    n=json.load(open(f'data/processed_h1/annotation_meta_h1_{s}.json'))['n_cells_qc']
    u=int(((d.cell_type.isin(BAD))&(~d.cell_type_merged.isin(BAD))).sum())
    print(s,u,n,round(100*u/n,2))"
```

**The cost in senders, at the primary threshold, on the merged `T/NK cells` stratum:**
SPLN14 **23 senders in 39,158 cells (0.059 %)** under the frozen fine-label call against 1,958
(5.000 %) under the merged one; SPLN43 **0 in 24,815 (0.000 %)** against 1,241 (5.001 %);
SPLN44 295 in 33,482 (0.881 %) against 1,675 (5.003 %). All three failures are **band**
failures below the 1 % floor — the signature of the mechanism, not of a thin population.

**What the two calls do to the answer: essentially nothing.** SF under N2+N5+N6 0.1582
(merged) vs 0.1630 (fine); controlled amplitude 0.0307 vs 0.0345 against bounds of 0.1094 and
0.1083; **0 of 7 modules over the bound either way**; 58 vs 61 reportable fits of 343.

---

## 9. Deviations and discrepancies added by Phase 10

Numbered `T*` so they cannot be confused with the gene-set `D*` rows, the Phase-8 `P*` rows,
the Phase-9 `H*` rows or the PI decision series `D-A … D-D`.

| # | Deviation / discrepancy | Evidence, and what was done |
|---|---|---|
| **T1** | **The H1 section cache was extended with per-module and merged-label sender flags.** | Without them the pre-registered D1 sensitivity axis and PI decision D-B's primary call could not run at all. `code/h1_cache_extend.py` copies flags **Phase 9 had already computed** (`data/processed_h1/senders_h1_*.csv`) into the `.npz`; nothing is recomputed from expression and no threshold is touched. Every pre-existing key is asserted byte-identical after the rewrite |
| **T2** | **A 5-dp rounding discrepancy in the per-module flags.** | The senders CSV stores the per-module score rounded to 5 dp, so recomputing the p95 flag from the stored score differs from the committed unrounded flag on **at most 20 cells of 291,577 (0.0069 %)**, all exact ties at the percentile boundary. The committed column is used at p95; p90 and p99 carry that ≤ 0.007 % tie noise. Per-section counts in `results/phase10_h1/cache_extend_rounding_*.csv`. **Reported, not patched** |
| **T3** | **`tierA_merged_p95` is in neither `N7_CALLS` nor `ALL9_CALLS`.** | Whether it replaces or is added to the frozen N7 axis is a frozen-list question (`PREREG_PHASE8.md` §0.1 open item 2). **The frozen lists were not edited.** `h1_run_phase10.ALL12` is a local extension in Phase-10 code, and every stage that iterates a frozen list was passed `--call` explicitly instead. **Raised, not resolved** |
| **T4** | **`phase3_null_diag.py` and `run_phase3_var._diag_job` derived a per-section seed as `int(sample[:4])`, which raises on `SPLN07`.** | Replaced by `phase3_null_diag.section_offset()`, which returns `int(sample[:4])` for the mouse names — **bit-identical to the previous behaviour, so every committed `results/phase3/null_destructiveness.csv` value is reproducible unchanged** — and a stable crc32 offset for non-numeric names |
| **T5** | **The three Phase-5 producers hardcode `results/phase5` and ignore `sasp_phase3.RESULTS`.** | Running them unwrapped on H1 would have **silently overwritten the committed mouse Phase-5 CSVs**, and `run_phase5_kernels.py` has no `--tag` to disambiguate. `code/h1_run_phase5.py` rebinds `phase5_common.RES5` and each runner's captured `RES`. **This is a latent hazard in the mouse producers and is recorded as one** |
| **T6** | **`sasp_real.parse_sample` cannot parse H1 section names**, so `Sec.meta["condition"]` is `"?"` and `["week"]` is `NaN`. | The `arm` and `week` columns of every H1 fit table are therefore **meaningless on this arm and must not be read**. No H1 analysis uses them. Declared rather than back-filled with an invented value |
| **T7** | **The ISO caches carry stale full-panel `senepy_score`, `tierApm__*` and `flag_pm_*`.** | SenePy and the seven per-module Tier A sets were deliberately not rescored on the intersected panel. Those keys are full-panel values sitting in an intersected-panel cache. Declared in both `*_prep_cache_iso.py` docstrings and in a `STALE_FULL_PANEL` list; **the `senepy_*` and `tierApm_*` calls must not be run against `cache3_h1_iso` / `cache3_m1_iso`**. Nothing mechanically blocks it |
| **T8** | **`_expand`'s seed rule is index-within-call-list, so the ISO refit's seeds depended on which call list was passed.** | Resolved in favour of reproducing the full-panel run's own seeds, so **every N2 match draw in the ISO refit is identical to the full-panel fit it is compared against**; the seed actually used is recorded per fit in `seeds_used.csv`, and `--seed-basis literal` gives the other reading. **A judgement call, recorded as one** |
| **T9** | **`code/run_phase5_wc.py:161` seeds with `hash(ct)`, which is `PYTHONHASHSEED`-randomised.** | Its per-cell seeds are therefore **not reproducible across invocations on either arm**, and the operator precedence is `seed + 7*j + ((101*hash(ct)) % 9973)`, probably not what was intended. **Pre-existing; found while assessing the Phase-5 producers for H1; the winner's-curse stage was not run on H1** (§11) |
| **T10** | **The five-seed DeepScence consensus (PI decision D-A) is compute-bound.** | See §3.5. Measured rate under the Phase-10 load: **~65 min per full-section run at 2 concurrent** (the hard maximum — 3 OOM-kill each other at this 57.7 GB cgroup), i.e. **~19 h for the 35 runs**, against ~8 h uncontended. **The seed count was NOT silently reduced.** The consensus is computed over whatever completed and is named `consensus_k<n>` with the count in every output row |
| **T11** | **Phase 9 deviation H8 recurred, and the first version of the Phase-10 DeepScence supervisor did not notice.** | Its concurrency guard counted processes but did **not** check memory headroom, so when a run was OOM-killed inside the 90 s settling window the count fell and the loop launched the next one into the same wall. **Five runs of seed 20260901 — SPLN44, SPLN14, SPLN30, SPLN43, SPLN24 — were launched at 90-second intervals between 20:57 and 21:03 UTC and all five were killed**, leaving no output and no error in the log (an OOM kill is silent to the shell). Found at 22:54 by noticing that the seed-major ordering had been violated in `logs/phase10/deepscence_seeds.log`. **Fixed and rerun**, and then fixed twice more when the fix itself failed — all three are recorded because they are the kind of failure that looks like a result: (i) require **≥ 22 GB of *anonymous* headroom** (`memory.stat`, not `memory.current`, not `free`) before each launch — 22 rather than 17 because SPLN24 peaks near 19 GB; (ii) re-derive the outstanding list every pass so a killed run is **retried**, not skipped; (iii) **`pgrep -fc PAT \|\| echo 0` emits `"0\n0"` when there is no match** — `pgrep -c` prints `0` *and* exits 1, so the `\|\| echo 0` fires as well — and `[ "0\n0" -ge 2 ]` is a shell **error**, not false, so the concurrency guard stopped guarding; piping through `wc -l` always yields one integer; (iv) **pgrep matches `/proc/pid/cmdline` with NUL separators mapped to spaces and that string can carry a trailing separator, so a `$` end-anchor never matches** — which launched `seed 20260902 SPLN24` **twice**, two minutes apart, and the pair blew the ceiling and one was OOM-killed. **The cost is ~2 h of wall clock and it is the direct reason §3.5 reports fewer completed seeds than five** |
| **T12** | **A committed mouse artefact was truncated by a concurrent Phase-10 process and had to be restored from git.** | `results/phase3/window.csv` lost its **135 out-of-band rows** at 22:09:14 UTC — the five sections excluded by §8 Test 3 (7239, 7250, 7361, 7448, 7450) — leaving only the 162 in-band rows, i.e. someone re-ran a window stage against the mouse results directory with a restricted section list. **Found by `git status`, restored with `git checkout -- results/phase3/window.csv`** (297 rows, 11 sections, md5 `4be813a9962af1da01358eb002b342f3`), and `results/phase3/`, `results/phase5/`, `genesets/`, `data/processed/` and `results/phase7_jobA/` are clean again. **I could not identify which process wrote it**, and that is the reportable part: `code/check_figures_guard.py` makes the figure policy enforceable, and **there is no equivalent guard over `results/`** — the directory this project's numbers actually come from. Recommended, not done here: a content-hash manifest over `results/phase3/` and `results/phase5/` on the same pattern |
---

## 10. Figures

Per PREREG §11: every figure writes a `*_data.csv` beside it carrying **every plotted number**,
both `.png` and `.pdf` are emitted, and the palette is the project's pre-validated categorical
theme reused through `code/sasp_palette.apply_style` — not re-chosen. Neither producer computes
a statistic; every value is read from a results file.

**Figure 5 — Two-arm replication** (`code/make_figure5.py` → `figures/figure5.{png,pdf}`,
`figures/figure5_data.csv`).
(a) Surviving fraction by null, M1 and H1 side by side, thirteen nulls, point = median over
that arm's reportable fits and bar = the IQR across fits — **labelled on the axis as an IQR
and not a confidence interval**, because that confusion is exactly what PREREG correction C-4
had to fix. (b) The §17 table as a forest plot: each Tier B module's controlled amplitude
against **its own arm's** 80 %-power detectable bound, drawn as a dotted line per arm — this
is the R2/R4 evaluation, and nothing crosses the bound on either arm. (c) The two geometric
predictions: the Poisson identity with the −1/2 line, and the λ̂ railing rate split into floor
and ceiling, which is where the arms visibly disagree.

**Figure 6 — DeepScence native vs remapped** (`code/make_figure6.py` → `figures/figure6.*`,
`figures/figure6_data.csv`). This carries §8.
(a) Caller agreement at full coverage, both arms, **before and after** conditioning on cell
type × depth decile — four series, log axis, chance at 1. The circular pair is drawn and
labelled CIRCULAR. (b) The sign anchor per section per arm, which is where P-ii was falsified.
(c) Score agreement against call-set agreement across seeds — the panel exists because the two
are **not** the same quantity, and the figure shows M1 at (0.9955, 0.7606) and H1 at
(0.372, 0.211).

**Figure guard.** `python3 code/check_figures_guard.py` reported **OK: all 52 committed
figures match** before Phase 10 began and again after Figures 5 and 6 were written, because
the guard verifies the files in its manifest and new files are outside it. The manifest was
re-snapshotted **once, at the end of the phase**; the before/after counts are in §12.

---

## 11. What was NOT run, and why

Listed rather than omitted, because an unrun test that looks like a passed one is the failure
mode this project has been correcting all week.

| Item | Status | Reason |
|---|---|---|
| **N8 on H1** (gene-set disjointness / scrambled response / CoreScence circularity) | **NOT RUNNABLE on H1, by three independent blockers** | (i) `run_phase3_n8.py:138` reads mouse expression through `sasp_real.load_expression`, whose path and file layout do not exist for H1; (ii) `:44` reads the **mouse** gene-set directory, so on a human panel the Tier A and Tier B sets would be mouse symbols with ~zero overlap; (iii) `:170` needs `genesets/E3_random_matched/`, which **has no human counterpart** — `code/build_genesets_human.py:369` prints *"E3 random size/expression-matched sets: NOT BUILT"*. Building them now would be creating a gene set after the pre-registration freeze, so **it was not done**. Two of N8's three products exist for H1 by other routes and are reported: **disjointness** = test A2 (byte-identical re-run of the frozen gate, PASS) and **CoreScence circularity** = 29/33 = 88.0 % measured natively (§3.1). **The scrambled-response product has no H1 counterpart and none is claimed** |
| **Winner's-curse cross-fit (`run_phase5_wc.py`) on H1** | **NOT RUN** | Its per-cell seeds are `PYTHONHASHSEED`-randomised (T9), so an H1 run would not be reproducible without pinning the interpreter's hash seed; and the M1 result it would be compared against was produced under an unpinned seed. Reporting a non-reproducible cross-arm comparison would be worse than reporting none |
| **`super --stage nulls` and `super --stage heldout` on H1** | **NOT RUN** | Budgeted at ~30 min/section and ~20 min respectively against a box already carrying the perturbation nulls, the composition-matched rerun and the DeepScence seed series. `kernel --stage heldout` **was** run, and it is the held-out test that bears on the conclusion |
| **P-vi at full section size on 6 of 7 H1 sections** | **NOT RUN at full size; run on the 20,000-cell panel** | A full-section `denoise=True` run is a DCA autoencoder pass over 200–400k cells on CPU (TensorFlow 2.4 cannot use this box's CUDA-12 GPU). §3.4 reports what the panel design gives and marks it as an approximation |
| **DeepScence `binarize=True`** | **NOT RUN** (unchanged from Phase 9, H11) | It costs 50 extra permuted forward passes per section, is not in the frozen settings list (PREREG §3.9), and the mouse arm did not run it either — so there would be no like-for-like comparison. **Declared rather than silently substituted**: A3's "each caller's own cutoff" is satisfied for DeepScence by the percentile calls only |
| **M1 re-called at the merged label family** | **NOT RUN** | `PREREG_PHASE8.md` §0.1 open item 3: whether M1 is re-called for cross-arm comparability, or the asymmetry is declared and carried, is not resolved. **The asymmetry is declared and carried.** The same fine/merged interaction exists on M1 and is smaller (`BIO_PHASE3.md` §1.1) |
| **M1's DeepScence score re-run at five seeds** | **DELIBERATELY NOT RUN** | PI decision D-A: M1 keeps its frozen single-seed score, because on M1 the single-seed score is not the source of instability (same-seed determinism control *r* = 0.99999913 / 0.99999995; seed-to-seed floor *r* = 0.99553). **The consequence is that the two arms do not use the same DeepScence estimator, and that asymmetry is itself a reportable finding** |
| **The three truncated DeepScence run-metadata JSONs** (SPLN21 / SPLN24 / SPLN30, Phase 9 H7) | **NOT RE-RUN** | The scores are complete and valid; only DeepScence's internal node choice and `reverse` flag were lost. Recovering them costs ~30 min per section and buys nothing this phase needs |

---

## 12. Environment, resources and provenance

**Interpreter.** The pinned Python 3.11 stack (`requirements.txt`). **No package was installed
into it.** The only other environment used is the isolated DCA venv built in Phase 9
(`/tmp/dca_env`, DCA 0.3.4 / TensorFlow 2.4.4 / Keras 2.4.3 under CPython 3.8.19), reached
through `code/_shims_dca_bridge` over a subprocess; **the main stack still has no TensorFlow
and never imports one**. TensorFlow 2.4 again could not use this box's GPU — it wants CUDA 11
and the box has an RTX PRO 4500 Blackwell with CUDA 12 — so DCA ran on CPU.

**Memory.** The ceiling is a **57.7 GiB cgroup** (`memory.max` = 61,999,996,928 B), and the
binding quantity is `anon` in `/sys/fs/cgroup/memory.stat`, **not** `memory.current` (which
counts reclaimable page cache) and **not** `free`. Observed:

| stage | concurrency | peak anonymous |
|---|---|---|
| `--stage main`, 12 calls × 7 sections | 8 workers | ~3.5 GB total |
| `--stage perm`, 1,000 permutations | 7 workers | ~8 GB total |
| `--stage perm_c1`, 1,000 permutations, 11 variants | 14 workers | ~14 GB total |
| `run_phase8_compmatch._job` | 6 workers | ~5 GB total |
| **`h1_deepscence_dca.py`, one full H1 section** | 1 | **~12–13 GB** |
| **three concurrent full-section DeepScence runs** | 3 | **> ceiling → OOM (deviation H8, recurred as T11)** |

**Wall time, one pass.** `main` 9.1 min (210 jobs, 8-way) · `perm` 66.1 min (14 jobs, 7-way,
1,000 permutations) · `perm_c1` 100.9 min (14 jobs, 14-way) · `var` see §2 · composition-matched
~100 min (63 jobs, 6-way) · Phase 5 kernels `section` 3.2–4.5 min, `proxdown` 1.3 min, `super`
`section` ~5 min, kernels `heldout` ~12 min · the ISO refit on both arms ~35 min ·
**one full-section DeepScence run 45–65 min under Phase-10 load** (vs 16.5–41 min uncontended
in Phase 9).

**`find` is `bfs` on this box and rejects relative timestamps**; every file listing here uses
explicit paths or `git status`.

**Provenance discipline.** Every number in this report was read from a file by a command
printed beside it, and `results/phase10_h1/two_arm_table.csv` carries the absolute source path
and the exact filter for every cell of the §17 table. §2's markdown is **generated** from that
CSV by `code/two_arm_table_md.py`, so it cannot drift from it.

**`python3 code/phase10_verify_report.py` re-reads every headline figure in this report from
its file and exits non-zero if any has drifted — 68 checks, 0 failed.** It is the Phase-10
analogue of Phase 9's `code/h1_verify_report.py` and exists for the same reason: this project
has repeatedly found its own written numbers not to match the files behind them
(`reports/AUDIT_PHASE8_FACTCHECK.md`). It covers both arms and both H1 sender calls — the fit
counts, the naive and controlled amplitudes, each arm's own 80 %-power bound, the count of
Tier B modules above it, SF under N1 / N3 / N4 / N3-var / N4-var / N3-tile / N4-tile / N2+N5+N6,
the composition-matched SFs, A7's primary response, the Poisson identity, the intersected-panel
SFs, the Phase-5 kernel and superposition results, and the §8 predictions P-ii, P-iii, P-iv and
P-vi.

**What Phase 10 wrote.** New paths only, with one declared exception:
- `code/h1_phase10.py`, `h1_run_phase10.py`, `h1_cache_extend.py`, `h1_headlines.py`,
  `h1_a6_sensitivity.py`, `h1_caller_unconditioned.py`, `h1_run_phase5.py`,
  `h1_run_compmatch10.py`, `h1_deepscence_consensus.py`, `h1_run_deepscence_seeds.sh`,
  `h1_run_dca_panel_all.sh`, `_h1_phase5_chain.sh`, `two_arm_table.py`,
  `two_arm_table_md.py`, `make_figure5.py`, `make_figure6.py`, `iso_panel.py`,
  `h1_callers_iso.py`, `h1_prep_cache_iso.py`, `m1_callers_iso.py`, `m1_prep_cache_iso.py`,
  `h1_phase10_iso.py`, `m1_phase10_iso.py`, `run_phase10_iso.py`, `iso_compare.py`
- `results/phase10_h1/`, `results/phase10_h1_iso/`, `results/phase10_m1_iso/`,
  `data/processed_h1/cache3_h1_iso/`, `data/processed_m1_iso/`, `logs/phase10/`
- `figures/figure5.*`, `figures/figure6.*` — **new files, outside the guard manifest until it
  was re-snapshotted once at the end**

**Modified, and declared:**
- `data/processed_h1/cache3_h1/*.npz` — 31 keys added per section by
  `code/h1_cache_extend.py`; **every pre-existing key asserted byte-identical** (T1, T2)
- `results/phase9_h1/d2_depth.csv` — **6 rows appended, 0 deletions, every pre-existing row
  byte-identical** (§3.4)
- `code/phase3_null_diag.py`, `code/run_phase3_var.py` — the per-section seed rule made
  arm-generic, **bit-identical for mouse sample names** (T4)
- `code/run_phase8_compmatch.py` — `ARMS['h1']` populated as that file's own docstring
  prescribes; it stays `frozen=True` and still requires `SASP_H1_UNFROZEN=1`
- `.gitignore` — `data/processed_m1_iso/` added, matching `data/processed_h1/`

**Not written by Phase 10:** `results/phase3/`, `results/phase5/`, `genesets/`,
`data/processed/`, `results/phase7_jobA/`. `git status` is clean for all five at the end of
the phase.

**Two concurrency facts, recorded because this repo had two agents in it at once.**
*(i)* `results/phase3/window.csv` **was** truncated by some process at 22:09:14 UTC and was
restored from git — deviation **T12**, and the reason that entry recommends a content-hash
guard over `results/` on the pattern of `code/check_figures_guard.py`.
*(ii)* A concurrent reproducibility agent committed `results/phase7_jobA/gate_result_human.json`
during this phase (commit `47c1bcd`, recording the venv path the gate resolves). **The
`corescence` block this report cites is unchanged by it** — `n_on_panel` 33,
`frozen_n_in_any_B` 29, `frozen_frac` 0.8788, with `B_secondary_senescence` accounting for 18
of the 29 — re-read after that commit. That agent also added `code/check_prohibitions.py` and
a pre-commit hook; **this report passes it with 0 violations** (`python3
code/check_prohibitions.py --backlog`), after one real fix: §3.4's `denoise=True` table needed
its seed-stability companion in the same paragraph (PREREG §10.10 / P26).

