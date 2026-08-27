# AUDIT — INDEPENDENT RE-DERIVATION OF EVERY HEADLINE NUMBER IN `reports/WRITING_PACK.md`

**Date 2026-08-27. Auditor: independent numbers pass.** Read-only except this file.
Nothing in `results/`, `figures/`, `code/`, `genesets/`, `data/` or any other report was
modified. `data/raw_h1/`, `results/phase9_h1/`, `data/processed_h1/` and `code/h1_*` were not
read.

**Method.** Every number below was recomputed from the files under `results/`, `figures/` and
`genesets/` without reading the value out of the report that quotes it. Where the pack supplies
a derivation command I wrote my own rather than running theirs, except where noted. The pack
carries **53 `[V]` marks** — `grep -o '\[V\]' reports/WRITING_PACK.md | wc -l` → 53, so the
"53 marked derivations" self-count is exact.

**Bottom line.** The **headline vector, λ̂ and its dependents, the caller decomposition, the
A7 per-family table, the composition ladder, the torus variant table, the calibration study,
the Moran result and the FPR all reproduce to the digit.** The paper's conclusion is not at
risk from anything in this audit. What is wrong is a **layer of provenance and scope
statements wrapped around correct numbers** — eight refutations, all of them in that layer,
plus one live contradiction inside the pre-registration.

---

# 0. LEAD — REFUTED and UNREPRODUCIBLE

## R1. `denoise=False` seed-to-seed agreement is **not** "Jaccard 0.76–0.99" — REFUTED

**Pack §5.9:** *"`denoise=False` agrees seed-to-seed at **Jaccard 0.76–0.99**."*

`results/phase8_d2/d2_stability.csv` contains exactly **one** `denoise=False` seed pair:

```
raw_seed0 vs raw_seed1, n=20000, pearson_r 0.99553, top5_jaccard 0.7606, top5_n_changed 272
```

There is no second raw-vs-raw row and no Jaccard of 0.99 anywhere in the file. **The "0.99"
is the Pearson r; the sentence merges a correlation and a set overlap into one range.**

`PREREG_PHASE8.md` P26 states it correctly — *"r 0.9955 / Jaccard 0.761 for `denoise=False`
across seeds"* — so the pack **degraded a correct statement**. The `denoise=True` half of the
same sentence is right: seeds 0 and 2 at Jaccard **0.6653**, seed 1 at **0.0000** against both,
2,000 of 2,000 cells changing, twice.

**Correct wording:** *"`denoise=False` reproduces across the two seeds at Pearson r = 0.996 and
top-5 % Jaccard 0.76."*

```bash
python3 -c "import pandas as pd; print(pd.read_csv('results/phase8_d2/d2_stability.csv').to_string(index=False))"
```

## R2. "RS_count … never outside ±0.010 of nominal" — REFUTED by the same sentence's own range

**Pack §5.7:** *"RS_count holds **0.033–0.060** across every window and every correlation scale
— never outside **±0.010** of nominal."*

Nominal is 0.05, so ±0.010 is [0.040, 0.060]. `results/phase3/var_sim_calibration.csv`,
`rs_count_reject_05`, all eight cells: **0.0325, 0.0550, 0.0350, 0.0600, 0.0400, 0.0425,
0.0525, 0.0550**. Two cells (rectangle s = 0.02 → **0.0325**; rectangle s = 0.15 → **0.0350**)
are below 0.040. The minimum is **0.0175 below nominal**, not 0.010.

The **range** (0.033–0.060) is correct and I confirm it; the **±0.010 gloss contradicts it.**
The Monte-Carlo SE at 400 replicates is √(.05·.95/400) = **0.0109**, so 0.0325 is 1.6 SE low —
statistically unremarkable, but the sentence as written is false.

**Correct wording:** *"RS_count holds 0.033–0.060 across every window and every correlation
scale, within 1.6 Monte-Carlo SE of nominal everywhere."*

## R3. The `figure4.png` md5 is stale and was transcribed from a pre-C6 report — REFUTED

**Pack §5.4, marked `[F]`:** *"the regenerated PNG is byte-identical (`figure4.png` md5
`000f34051112aff4fed293fe7a5b25c2`)"*.

```
md5sum figures/figure4.png  ->  d44fac63411d6c30a42c40894a287f17
```

`d44fac63…` is also what `results/phase3/m1_final_audit.txt` §5 records. The quoted hash
`000f3405…` occurs in exactly two places in the repo: `reports/CS_PHASE7_C1.md` L308 (a
**pre-C6** document) and the pack line that copied it. **It was read from a report, not from a
file — the exact failure mode the `[F]` mark is supposed to exclude.**

The underlying claim is nevertheless **sound**: `sha256sum figures/figure4.png` =
`76718f1ec5c1625aa977d6df8d831fd23cfe33f993a8b206ce09a17860044829`, which matches
`figures/.committed_manifest.json`, and `python3 code/check_figures_guard.py` exits **0** with
*"OK: all 52 committed figures match"*. **Drop the md5 or replace it with the on-disk one.**

## R4. "the guard now covers everything" — REFUTED; 10 uncovered figure artefacts on disk

**Pack §6:** *"**`git ls-files figures/` = 52; 52 artefacts on disk** — the Phase 8 figures have
been committed, so the guard now covers everything."*

- `git ls-files figures/` = **52** ✔
- top-level files in `figures/` = **53** — the 52 plus `.committed_manifest.json`, which is the
  guard's own manifest, not an artefact. So "52 artefacts on disk" is defensible.
- **but `figures/` also holds an untracked subdirectory `figures/revised_candidates/` with 10
  further figure artefacts** — `figure2a_REVISED.pdf`, `figure2b_REVISED.{pdf,png}`,
  `figure2c_REVISED.{pdf,png}`, `figure2d_REVISED.pdf`, `figure4_REVISED.{pdf,png}`,
  `figure4_data_REVISED.csv`, `README.md`, timestamped 05:22–09:15. `git ls-files
  figures/revised_candidates | wc -l` = **0**.

The guard enumerates `git ls-files figures/`, so it does not see them. **"Covers everything" is
exactly the overreach the pack's own checklist item 30 forbids** ("state the count it actually
covers"). Say: *"the guard covers the 52 committed artefacts and passes; a further 10 candidate
figures in `figures/revised_candidates/` are untracked and outside its scope."*

## R5. The intersection-matrix appendix line — "196 cells, every Tier A × Tier B cell zero" — REFUTED as written

**Pack §5.9 appendix / §2.0:** *"the intersection matrix
(`figures/figure_gs1_intersection_matrix_data.csv`, **196 cells**, every Tier A × Tier B cell
zero, but 18 of 21 Tier B pairs non-disjoint)"*.

The file has **224 rows**, in three panels × two arms:

| panel | rows | non-zero |
|---|---|---|
| `BxB` | 84 (42 mouse, 42 human) | mouse 36 of 42 ordered = **18 of 21 unordered pairs** ✔ |
| `A0xB` (pre-disjointness / failed variants) | 28 | 18 |
| `frozenAxB` | 112 | **55** |

Within `frozenAxB`, the eight row-labels are `A_SENDER_FINAL_strict (PRIMARY)` plus the seven
`A_sender_for_<module> (sensitivity)` sets. **`A_SENDER_FINAL_strict` × all seven modules is
0 in both arms** — max 0, sum 0, 14 cells. The gate PASS is real. The 55 non-zero cells are
**the per-module sensitivity sender sets**, which are only required to be disjoint from *their
own* module. Nothing is broken, but "every Tier A × Tier B cell zero" describes 14 of 224
cells and reads as describing all of them. **196 is reproducible from no subset of this file.**

**Correct wording:** *"…224 cells across three panels; the 14 `A_SENDER_FINAL_strict` × Tier B
cells are all zero in both arms, and 18 of 21 Tier B pairs are non-disjoint."*

## R6. "pooled CI half-width ±0.018 SD" — subset-dependent, and wrong for the response the sentence is about

**Pack §5.5, caveat 2:** *"A single fit resolves ±0.1346 SD naive / ±0.1336 SD conditioned
(median CI half-width); **pooled CI half-width ±0.018 SD**."*

The per-fit half-widths are right (`a7_summary.csv`, `median_CI_halfwidth`, `all_controls`:
base **0.1346**, n6n5 **0.1336** — note `a7_verdict.txt` prints **0.1345** for the naive one,
a third-decimal disagreement between the two frozen files).

The pooled half-width, `(clustered_hi − clustered_lo)/2`, is **not 0.018 in general and not
0.018 for `all_controls`**:

| | base | n6 | n5 | n6n5 | n2 |
|---|---|---|---|---|---|
| `all_controls` | **0.0562** | 0.0461 | 0.0224 | **0.0215** | 0.0471 |
| `neg_control_codeword` | 0.0481 | 0.0398 | 0.0190 | 0.0181 | 0.0410 |
| `neg_control_probe` | 0.0303 | 0.0294 | **0.0176** | **0.0176** | 0.0270 |
| `genomic_control` | 0.0251 | 0.0251 | 0.0233 | 0.0233 | 0.0239 |
| `neg_probe_rate` | 0.0198 | 0.0198 | 0.0156 | 0.0157 | 0.0195 |

Range over the 25 control rows: **0.0156 – 0.0562**. The 0.018 is the *conditioned probe /
codeword* row. For the number the caveat sits beside — `all_controls`, `base`, −0.0744 — the
pooled half-width is **0.0562**, three times the quoted figure. **Name the response and design,
or give the range.**

## R7. Four `run_phase3_nulls.py` line citations are wrong

Every other line citation I checked is correct — `summarize_phase3.py:99` really is
`np.quantile(v,[.25,.5,.75])`; `summarize_phase3.py:221` really is the `medlam` print;
`summarize_moran.py:183/184` really are the two `spearmanr` calls;
`summarize_phase3_c1.py::reportable` really is `beta_naive > 0 & beta_base_lo > 0`. But all four
citations into `code/run_phase3_nulls.py` (541 lines, unmodified vs HEAD) miss:

| pack says | actual |
|---|---|
| `WINDOW_UM`, `:59` | **`:52`** (`:59` is blank) |
| window "enforced at `:169`" | **`:143`** (`self.d_obs <= WINDOW_UM`); `:169` is `def fit_cell` |
| railing "exact index equality at either end (`:239`)" | **`:204`** (`lam_railed=int(t0==0 or t0==sf.lam.size-1)`) |
| λ grid "`:59-93`" | **`:52–72`** (`LAM_LO_FLOOR = 7.0` at `:53`, `lam_grid()` at `:65–72`) |

**The substance is fully confirmed:** `WINDOW_UM = 100.0`, `LAM_LO_FLOOR = 7.0`, ceiling
`dmax/2 = 50.0`, `N_LAM = 40` log-spaced, `MIN_RECEIVERS = 2000`, `N_BOOT = 400`, railing =
index equality at either end. Only the pointers are wrong — and a wrong pointer into frozen
code is how this project produced three of its previous errors.

## R8. The pack's map of where the pre-C6 A7 set is still in circulation is wrong in **both** named places, and misses four live ones

**Pack §5.5:** *"`PREREG_PHASE8.md` §10.1 paragraph 1 and `PHASE8_ROADMAP_STATUS.md` 8.5b both
quote this pre-C6 set."*

- `PHASE8_ROADMAP_STATUS.md` **8.5b is fully corrected** — frozen digits in bold, pre-C6 values
  struck through, with a dated reconciliation note. It is **not** in circulation.
- `PREREG_PHASE8.md` **§10.1 ¶1 carries a dated `[Corrected 2026-08-27 — see §0.0 item C-1]`
  block** naming every frozen replacement. Marked, not live.

**The four that *are* live and unmarked:**

| location | stale content | frozen |
|---|---|---|
| `PREREG_PHASE8.md` **P2** (:902) | −0.070 [−0.128,−0.012] p 0.023; N2 −0.061 p 0.020; N5 +0.007 p 0.41; biological +0.314 / +0.077 — **and** *"A7 was run at 05:19 on the pre-C6 sender calls … **A7 must be re-run after 8.7**"* | −0.0744 [−0.1306,−0.0182] p 0.0145; −0.0642 p 0.0124; +0.0038 p 0.715; 0.3120 / 0.0795. **A7 *was* re-run — the frozen file is 09:06.** The instruction is stale, not just the digits |
| `PREREG_PHASE8.md` **P3** (:903) | *"0.091 / 0.103 / 0.109 / **0.127** / 0.164 … the clean-null subset is **9.1–12.7 %** … **Quote the range as 9–13 %** with 16 % as the `neg_probe_rate` outlier"* | probe is **0.145**; count-based subset **9.1–14.5 %**. **P3's own header says "9–16 %" while its body instructs "9–13 %"** |
| `PREREG_PHASE8.md` **P24** (:924) | *"N2 matched decoys leave the technical gradient **~80 % intact** (−0.061 of −0.070 SD)"* | −0.0642 / −0.0744 = **86.3 %**. The file's own C-1 row fixes §10.1 to 86 % and leaves P24 at 80 % |
| `NOVELTY_ASSESSMENT.md` :305, and **§4 O2** (:582) | −0.070 SD in the body unmarked; O2: *"`neg_probe_rate` flat naively (**+0.014, p=0.079**)"* | **+0.0113, p = 0.232.** The claim survives (n.s. either way) but p moved by a factor of three |

**And `COMPLETED_TASKS.md`:161 asserts of `PREREG_PHASE8.md`: "Residual stale digits: 0."**
That is refuted by P2, P3 and P24 in the same file.

**P3 is the one that bites.** It is a *pre-registered deviation* instructing the writer to quote
**9–13 %**, which contradicts `WRITING_PACK` §5.5, §7 row 5, and the frozen `a7_summary.csv`.
A drafter reading the pre-registration as authority will publish the pre-C6 range.

## R9. "IQR is [7.0, 50.0] µm for every one of these" — false for two of the five λ̂ definitions

**Pack §0.1**, immediately under the five-row alternatives table.

| definition | n | median | **IQR** |
|---|---|---|---|
| all 315 primary fits | 315 | 14.732 | **[7.00, 50.00]** ✔ |
| 153 reportable | 153 | 16.069 | **[7.00, 50.00]** ✔ |
| interior of 315 (the summariser's printed value) | 126 | 17.063 | **[10.67, 28.68]** ✘ |
| interior of 153 | 67 | 14.994 | **[10.63, 25.43]** ✘ |
| N7 table row `tierA_p95` | = the 315 | 14.7 | **[7.00, 50.00]** ✔ |

Interior medians exclude the railed fits **by construction**, so their quartiles cannot be at
the rails. The claim is true for the three definitions that matter and false for the two the
pack is arguing against. **The authoritative choice is unaffected** — 14.7321 over 315 does
carry [7.0, 50.0] and 60 % railed — but the blanket "every one of these" is wrong.

## R10 (internal contradiction). 2.35× vs 2.36× — the pack disagrees with itself

- **§5.7 bullet:** *"0.118/0.05 = 2.36; `CORRECTIONS.md` §C.4 rounds to 2.35×, `SUBMISSION_PATCH`
  to 2.4×."* — implies **2.36 is exact and 2.35 is a rounding**.
- **§7 resolution banner:** *"the inflation is standardised repo-wide on **2.35×, the exact
  0.1175 / 0.05** — '2.4×' and '2.36×' both rounded it up."*

The file value is `torus_tile4x4_reject_05 = 0.1175` (irregular, s = 0.30). **0.1175 / 0.05 =
2.35 exactly.** §7 is right; §5.7 has it backwards. Delete the §5.7 clause.

---

# 1. TRUE BUT MISLEADING

| # | statement | why it misleads | say instead |
|---|---|---|---|
| M1 | §5.6: *"**All four rows below are the same 42 fits / 33 reportable**, pooled scope"* | `compmatch_reruns.csv` shows `comp` and `full` pooled at **`n_fits` 210, `n_reportable` 165** — because they are 42 × 5 seeds and 33 × 5. Anyone checking the file against the sentence concludes the populations differ, and the "factor of fifty" is the load-bearing claim | *"the same 42 fits / 33 reportable **per seed**; the matched rows are pooled over the five frozen seeds"* |
| M2 | §5.5 caveat 3: the reportable filter *"admits **4.8 %**, identical across all five"* | the value is **0.0485** — 4.9 % at the precision the adjacent "3.0–13.3 %" is given at. Truncation, not rounding | 4.9 %, or 0.0485 |
| M3 | §5.2 / §5.8: Poisson identity *"slope −0.525, r² 0.9843"* | **subset-dependent and the subset is not named.** `summary_phase3.txt` §10 gives four: ALL sender definitions (−0.525, 0.9843, n=77); `cdkn1a_pos` only (**−0.486, 0.9982**); Tier A percentile calls (−0.517, 0.9958); **in-band only (−0.529, 0.9833)**. §5.8 offers this as *the cross-arm geometric prediction* — it must not be a moving target | *"−0.525, r² 0.984, over all 77 section × sender-definition cells"* |
| M4 | §5.4: *"Real-significant rate on real data: COMMOT 0.224 …"* | also null-dependent: the ligand-permuted nulls run on **5,629** interactions, the others on **6,032**, so COMMOT is 0.2236 vs **0.2299**, CellChat 0.2833 vs 0.2887, NCEM 0.01442 vs 0.01492, SpaTalk 0.2472 vs 0.2546 | give one basis, or the two-value range |
| M5 | §5.5: *"**0.362 SD** … Moran's I could not have detected the paper's headline effect"* | 0.362 is the **median over 22 section × call rows spanning 0.308–1.070**; restricted to `tierA_p95` it is 0.308–0.418. The paper's **naive** amplitude is 0.329, which sits *inside* that range — in the best sections Moran could have resolved it. The claim is safe against the *controlled* amplitude (0.029) by two orders of magnitude | *"at the median section Moran's I could not have detected it (0.362 SD vs 0.277); in the best-powered section the bound falls to 0.308 SD, still above the conditioned amplitude by 10×"* |
| M6 | `CS_PHASE8_MORAN.md`:246 *"A7 amplitude β_z used −0.070, range −0.279 – +0.164"* | reads as the **pre-C6 pooled A7 value**. It is not: it is the **median of the 22 per-section `a7_beta_z`** in `moran_kernel_power.csv` (−0.0704, min −0.2787, max +0.1636). A coincidence, but a reviewer who has read the corrections ledger will flag it | label it *"median per-section A7 β_z"* |
| M7 | §5.3: *"no family beats it in more than **54 %** of folds"* | gaussian is **0.544** | "in more than 55 %", or quote 0.544 |
| M8 | §5.5 caveat 1 and forbidden item 2 are right that `all_controls` = probe + codeword + genomic; but §5.5's per-family table gives `n features` = 670 for the pooled row | 609 + 21 + 40 = 670 ✔ — confirmed, no issue. Recorded only because it is the arithmetic the whole R1 correction rests on | — |

---

# 2. CONFIRMED — priority 1, the headline vector

```bash
python3 - <<'PY'
import pandas as pd, numpy as np
d=pd.read_csv('results/phase3/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
p=d[(d.call=='tierA_p95')&(d.sec.isin(['7259','7260','7001','7248','7352','7435']))&(d.stratum=='all')]
r=p[(p.beta_naive>0)&(p.beta_base_lo>0)]
se=(r.beta_n2n5n6_hi-r.beta_n2n5n6_lo)/(2*1.959963985)/r.sd_y
print(len(p), len(r), np.median(r.beta_naive/r.sd_y), np.median(r.beta_n2n5n6/r.sd_y),
      np.median(se), (1.959963985+0.8416212336)*np.median(se),
      np.median(r.sf_n2n5n6), np.percentile(r.sf_n2n5n6,[25,75]),
      int(((r.beta_n2n5n6>0)&(r.beta_n2n5n6_lo>0)).sum()))
PY
```

| quantity | pack | **re-derived** | verdict |
|---|---|---|---|
| primary fits | 315 | **315** | CONFIRMED |
| reportable (`beta_naive>0 & beta_base_lo>0`) | 153 | **153** | CONFIRMED |
| **naive amplitude** | 0.3288 → 0.329 | **0.328795** | CONFIRMED |
| **controlled amplitude** | 0.0288 → 0.029 | **0.028827** | CONFIRMED |
| ctrl amplitude IQR | [−0.0066, +0.0845] | **[−0.006613, +0.084465]** | CONFIRMED |
| median SE of controlled amplitude | 0.0654 | **0.065427** | CONFIRMED |
| **80 %-power detectable bound** | 0.1833 = 2.802 × 0.0654 | **0.183299**; z₀.₉₇₅+z₀.₈₀ = **2.801585** | CONFIRMED |
| bound ÷ estimate | 6.4× | **6.359×** | CONFIRMED |
| **SF under N2+N5+N6** | 0.0885 → 0.088 | **0.088468** | CONFIRMED |
| its IQR | [−0.017, 0.234] | **[−0.016643, 0.233783]** | CONFIRMED |
| controlled fits positive **and** CI excludes 0 | 13 of 153 | **13** | CONFIRMED |

Independent of `m1_final_audit.txt`, which agrees on every one of these. The full null battery
(§5.2) reproduces from `sf_summary.csv` row for row: N1 0.707 [0.416,0.865]; N1-on-residual
0.989; N2 0.952 [0.919,0.978]; N3 0.999 [0.989,1.006]; N4 0.947 [0.804,1.039]; N5 0.115
[−0.034,0.258]; N6 0.471 [0.183,0.742]; zonation 0.843; N5+N6 0.084 [−0.013,0.233]; N8 0.916.
`frac_le_0` / `frac_gt_05` columns match the ≤0 / >0.5 columns. **`sf_summary.csv` has columns
`subset,null,n,q25,median,q75,frac_le_0,frac_gt_05` — there is no CI column, and
`summarize_phase3.py:99` is a bare `np.quantile(v,[.25,.5,.75])` over per-fit point estimates.
§0.7 (the brackets are IQRs, and no bootstrap can retroactively make them CIs) is CONFIRMED.**

Where the naive gradient comes from (`summary_phase3.txt` §11): monotone in **32/42** ✔;
cell-type intercepts alone **0.348** ✔; full N5 block **0.042** ✔. Ripley ratio **1.1237** ✔.
λ̂-is-not-density: r² median **0.016**, slope median **−0.091** ✔.

---

# 3. CONFIRMED — priority 2, λ̂ = 14.7321 µm and its dependents

| definition | pack | **re-derived** |
|---|---|---|
| pooled median over the **315** | **14.73** | **14.7321271** ✔ |
| median over the 153 reportable | 16.07 | **16.068978** ✔ |
| interior median over 315 (what the summariser prints in §1) | 17.1 | **17.062698** ✔ |
| interior median over 153 | 14.99 | **14.994266** ✔ |
| N7 table row `tierA_p95`, `medlam` | 14.7 | **14.7** ✔ (`summarize_phase3.py:221` = `dd.lam_naive.median()`, group = the 315) |
| railing rate | 60 % | **0.6000 over 315** (0.5621 over 153) ✔ |
| railed counts | 103 floor / 86 ceiling | **103 / 86**, 126 unrailed ✔ |

**The censoring argument in §0.1 point 4 is exact.** The median of 315 is the 158th order
statistic; 103 fits sit at the 7 µm floor, so it is the **55th of the 126 unrailed values** —
strictly interior. 14.7321 is a legitimate median of a censored sample.

**Dependents, all CONFIRMED:**

| claim | arithmetic |
|---|---|
| N3-var 2,215 µm = **150 λ̂** | 2215 / 14.7321 = **150.35** ✔ |
| 1,200 µm tile side = **~81 λ̂** | 1200 / 14.7321 = **81.46** ✔ |
| 100 µm window = **≈ 6.8 λ̂** | 100 / 14.7321 = **6.788** ✔ |
| = **14×** the 7 µm floor, **2×** the 50 µm ceiling | 14.29, 2.00 ✔ |
| N4-var 3,395 µm = **230 λ̂** | 230.45 ✔ |
| N3-occ 28 µm = **1.9 λ̂** | 1.901 ✔ |

**The withdrawal of 15.7 is CONFIRMED, and so is its diagnosed provenance.** 2215/141 =
**15.709** (the circular back-derivation). And the pack's alternative account reproduces
exactly: over `results/phase3_pre_c6/main_fits.csv`, `tierA_p95`, six in-band sections,
**including the zonation-stratified rows (441 rows: 315 `all` + 3 × 42 hepatocyte strata)**,
the interior median is **15.7164**. Pre-C6, interior, and pseudo-replicated — not a defensible
estimand. The pre-C6 `stratum=='all'` values are 13.303 pooled / 15.566 interior, neither of
which is 15.7 either. **REFUTED as an estimand; its arithmetic origin CONFIRMED.**

The mandatory caveat — IQR **[7.0, 50.0]**, **60 % railed** — is CONFIRMED for the chosen
definition. See **R9** for the two rows where the blanket IQR claim fails.

---

# 4. CONFIRMED — priority 3, the caller-agreement decomposition

I re-implemented the Mantel–Haenszel pooling from the significance CSVs without reading
`summarize_caller_coverage.py`'s `pooled()` first, then checked mine against it:
`obs = Σ n_both`, `exp = Σ exp_both_stratified`, `z = (obs−exp)/√(Σ sd_both_stratified²)`,
two-sided normal p, `n_above_chance = #(ratio_stratified > 1)`, sign test = binomial.

## All six labelled bases, **4-pair** — every cell reproduces `caller_coverage_gate_headline.csv`

| basis | n_sec | n_val | pooled | z | p | band | median | above | sign p |
|---|---|---|---|---|---|---|---|---|---|
| 2-section, pre-C6 (PUBLISHED) | 2 | 8 | **1.0398** | 1.76 | 0.078 | 0.932–1.369 | 1.010 | 4/8 | 1.000 |
| 11-section, pre-C6 | 11 | 44 | **1.1292** | 13.35 | 1.12e-40 | 0.700–1.711 | 1.155 | 29/44 | 0.0488 |
| 2-section, post-C6 (FROZEN) | 2 | 8 | **1.1312** | 5.80 | 6.50e-09 | 0.979–1.442 | 1.099 | 6/8 | 0.289 |
| **11-section, post-C6 (FROZEN)** | 11 | 44 | **1.2120** | 21.92 | **1.84e-106** | 0.751–2.198 | 1.190 | 35/44 | 1.06e-04 |
| 6-section in-band, pre-C6 | 6 | 24 | **1.1146** | 8.99 | 2.56e-19 | 0.775–1.374 | 1.131 | 16/24 | 0.152 |
| 6-section in-band, post-C6 | 6 | 24 | **1.1669** | 13.09 | 3.63e-39 | 0.811–1.565 | 1.168 | 18/24 | 0.0227 |

## The **3-pair** basis — emitted by no file; all four re-derived exactly

| basis | pack | **re-derived** | z | p | band |
|---|---|---|---|---|---|
| 2-section pre-C6 (published) | 1.0299 | **1.0299** | 1.27 | 0.203 | **0.932–1.221** |
| 11-section pre-C6 | 1.1182 | **1.1182** | 11.49 | 1.44e-30 | 0.700–1.711 |
| 2-section post-C6 (frozen) | 1.1283 | **1.1283** | 5.47 | 4.40e-08 | 0.982–1.442 |
| 11-section post-C6 (frozen) | 1.2122 | **1.2122** | 20.62 | 1.84e-94 | 0.751–2.198 |

**The 3-pair 2-section pre-C6 band is [0.932, 1.221] — i.e. the published "0.93–1.22×" band is
literally the three Tier-A pairs.** §0.6's identification is CONFIRMED.

**One thing the pack does not say, and should:** at 11 sections post-C6 the 3-pair and 4-pair
**ratios coincide to three decimals** (1.2122 vs 1.2120 → both "1.212"); only the **p-values**
diverge (1.84e-94 vs 1.84e-106). So the "never mix bases in one sentence" rule is about the
exponent, not the ratio — worth stating, because a reader who checks only the ratio will
conclude the warning is pedantic.

**Also confirmed:** "four of six pairs sit at 0.93–1.22×" **is** true on the 2-section base at
the per-pair pooled level (0.935, 1.017, 1.103, 1.168 in; 1.693, 2.05 out) and **is** false at
11 sections (only 1.211 and 0.972 in — 2 of 6). Checklist item 14's scoping is exact.

## Per pair, frozen Tier A, 11 sections (`caller_coverage_gate.csv`) — all six CONFIRMED

Tier A × `Cdkn1a`⁺ **1.471** (z 19.45, p 2.74e-84, 1.085–2.198, 11/11, 9 sig above, 0 below);
Tier A × DeepScence **1.288** (19.23, 2.13e-82, 1.096–1.660, 11/11, 10, 0);
SenePy × `Cdkn1a`⁺ **1.211** (7.43, 1.11e-13, 0.842–1.391, 9/11, 7, 0);
Tier A × SenePy **0.972** (−1.63, **0.104**, 0.751–1.179, **4/11**, 1, 3);
SenePy × DeepScence **0.737** (−15.08, 2.33e-51, 0.332–2.150, 1/11, 1, **10**);
DeepScence × `Cdkn1a`⁺ (circular) **1.255** (10.53, 6.24e-26, **0.963–2.849, median 1.071**).

Two-section frozen: Tier A × `Cdkn1a`⁺ 1.017 → **1.300**; × DeepScence 1.103 → **1.179**;
× SenePy 0.935 → **1.007**. ✔

**The dead plank is confirmed dead:** pre-C6 Tier A × SenePy pooled **0.914** with
`n_sections_above_chance = 0`, i.e. below chance in 11 of 11; frozen **0.972, p = 0.104, above
chance in 4 of 11**.

**"22 of 33" → 20 of 33** is CONFIRMED: pre-C6 3-pair 11-section gives **20/33** (0/11 Tier A ×
SenePy + 11/11 × DeepScence + 9/11 × `Cdkn1a`⁺); the frozen 3-pair figure is **26/33**.

**Depth camps** (`caller_within_type_depth_bias_11sections.csv`, Q5/Q1 of `enrichment` within
cell type), all four CONFIRMED and all 11/11 in direction:
SenePy **10.583–41.739** top; `Cdkn1a`⁺ **4.190–42.360** top; Tier A **0.146–0.317** bottom;
DeepScence **0.218–0.795** bottom.

---

# 5. CONFIRMED — priority 4, A7 per control family (frozen 09:06)

`results/phase3/a7_summary.csv` is 09:06; `results/phase3_pre_c6/a7_summary.csv` is 05:19.
825 control fits (**165 × 5 responses**) + **1,155** module fits, 11 sections × 2 sender calls,
9 receiver types — confirmed from `a7_control_probe_fits.csv` (825 rows, 5 responses, 2 calls,
11 sections, 9 cell types).

**I re-derived the clustered mean independently** — mean over 11 per-section means of β/sd_y
with a t₁₀ CI — and got `all_controls` base **−0.0744 [−0.1306, −0.0182], p = 0.01451** and
`neg_control_probe` base **−0.0225 [−0.0527, +0.0078], p = 0.1287**, matching the file to four
decimals. The "section-clustered mean with a t-CI" definition is therefore confirmed, not just
the digits.

| response | n feat | β/sd_y | 95 % CI | p | flat? |
|---|---|---|---|---|---|
| `all_controls` | 670 | **−0.0744** | [−0.1306, −0.0182] | **0.0145** | no |
| `neg_control_codeword` | 609 | **−0.0604** | [−0.1085, −0.0123] | 0.0188 | no |
| `genomic_control` | 21 | **−0.0307** | [−0.0558, −0.0056] | 0.0213 | no |
| `neg_control_probe` | 40 | **−0.0225** | [−0.0527, **+0.0078**] | **0.129** | **YES** |
| `neg_probe_rate` | — | **+0.0113** | [−0.0085, +0.0310] | 0.232 | yes |
| BIOLOGICAL MODULES | 1,155 fits | **+0.2767** | [+0.2032, +0.3502] | 8.0e-06 | — |

`all_controls` across designs: base −0.0744 (0.0145); n6 −0.0625 [−0.1085,−0.0164] (0.0128);
**n5 +0.0038 [−0.0186,+0.0261] (0.715)**; n6n5 +0.0053 [−0.0162,+0.0268] (0.595);
**n2 −0.0642 [−0.1113,−0.0172] (0.0124)**. Modules: base +0.2767, n6 +0.1124, n5 +0.0694,
**n6n5 +0.0310**, n2 +0.2662. **All CONFIRMED.**

N2 leaves **86.3 %** of the gradient (−0.0642 / −0.0744) — the pack's "86 %" is right and
`PREREG_PHASE8.md` P24's "~80 %" is the pre-C6 ratio (see R8).

**Naive biological amplitude — the four-value resolution CONFIRMED.** Frozen clustered signed
mean **0.2767**; frozen median |β|/sd **0.3120** (`a7_verdict.txt`); pre-C6 **0.2914** and
**0.3143** (`PREREG` rounds to 0.314). Conditioned counterparts **0.0310** / **0.0795** frozen,
**0.0356** / **0.0770** pre-C6. The Moran power bound 0.362 exceeds all four ✔.

**Sparsity CONFIRMED** (`a7_verdict.txt`): `all_controls` 0.0588 c/cell, 5.26 % non-zero;
codeword 0.0428 / 3.90 %; genomic 0.0094 / 0.90 %; probe **0.0067 / 0.65 %**. Codewords carry
0.0428/0.0588 = **72.8 %** of the pooled counts ("~73 %" ✔).

**Pre-C6 values the pack lists as forbidden — all CONFIRMED as the 05:19 file's contents:**
all_controls base −0.0697 [−0.1276,−0.0118] p 0.023; n2 −0.0611 p 0.0204; n6n5 +0.0068 p 0.411;
probe −0.0177 [−0.0453,+0.0099] p 0.183; codeword −0.0549 p 0.0389; genomic −0.0337 p 0.00389;
modules base +0.2914, n6n5 +0.0356. **But see R8** — the pack's map of where they are still in
circulation is wrong.

**Audit R5's moot-ness CONFIRMED:** frozen `neg_probe_rate` under n6n5 is **+0.0097
[−0.0060, +0.0253], p = 0.199** — CI includes zero, so "every control family is flat under
+N6+N5" is now simply true. Pre-C6 it was +0.0108 [+0.0021, +0.0195], p = 0.0204.

---

# 6. CONFIRMED — priority 5, the composition ladder

`results/phase3/compmatch_reruns.csv`, `row_type == 'summary'`, `call == 'tierA_p95'`,
`scope_kind == 'pooled'`:

| variant | SF | 95 % CI | share removed |
|---|---|---|---|
| `comp` (20-NN composition as **matched decoys**) | **0.983735** | [0.973037, 0.994166] | **1.63 %** |
| `comp_adj` (same variables as **covariates**) | **0.498892** | [0.421492, 0.606284] | **50.11 %** |
| `type_adj` (receiver's own cell type) | **0.341450** | [0.235641, 0.401643] | **65.86 %** |
| `typecomp_adj` (type + 20-NN composition) | **0.146139** | [0.052005, 0.246070] | **85.39 %** |

**85.39 / 1.63 = 52.5 — "a factor of fifty"** ✔. Within cell type the matched protocol removes
**3.57 %** (SF 0.964666, 748 reportable) ✔; the published N2 set (`full`) removes **1.45 %**
pooled and **5.06 %** by cell type ✔.

`tierApm_p95` sensitivity (n = 34 reportable): comp **0.986066** (1.39 %), comp_adj **0.514046**,
type_adj **0.335009** (66.50 %), typecomp_adj **0.204167** (79.58 %) ✔. Largest paired
difference is typecomp_adj, **0.0581 < 0.07** ✔.

Seeds **20260901–20260905** ✔. Across-seed spread: pooled **0.983695–0.983965, sd 1.20e-04**;
by cell type **0.962919–0.964825, sd 8.26e-04** ✔. Adjusted variants are seed-free (`seed = -1`,
`n_seeds = 1`) ✔.

Balance: median max |SMD| **0.091612 → 0.035222**; median match rate **0.999871**; SMD gate pass
**1.000**; worst `max_smd_after` over all **6,237** rows of `compmatch_fits.csv` = **0.0759** ✔.
Caliper 0.25 SD, 400 bootstrap replicates, 100 blocks, 100 µm window ✔.

See **M1** for the one thing to reword.

---

# 7. CONFIRMED — priority 6, the torus variant table and calibration

**All 15 rows × 5 columns of `results/phase3/sf_summary_var.csv` reproduce.** Spot values:
N3 published **0.9989** [0.9889,1.0058], keeps-a-neighbour **0.7716**, displacement **2,909.9**,
full **1.0015**; N3-var **0.9960** [0.9754,1.0074], **1.0000**, **2,214.9**, full **0.9970**;
N3-occ **0.3024** [~0, 0.7338], displacement **28.30**, full **0.2874**; N4-occ **0.1833**
[~0, 0.5586], **24.69**, full **0.1979**; N4-var **0.9851** [0.9581,1.0033], **3,395.2**, full
**0.9986**. Window-matched N3-var **0.9948** [0.9755,1.0077]; N4-var **0.9849** ✔.

n = **153** on whole-section rows and **136** on the two tile rows — **17** fits drop below the
2,000-receiver floor ✔. `full_sf_n_perm = 200` on the two `var` rows only ✔.
Real-data tiled-vs-whole rejection **0.8015 vs 0.8235** ✔.

**Destructiveness, two columns never merged** (`null_destructiveness.csv`, medians over the six
in-band sections):

| null | out of tissue = 1 − `frac_in_occupancy` | lost every real neighbour ≤ 100 µm |
|---|---|---|
| N3 bbox | **35.53 %** | **22.84 %** |
| N4 bbox | **19.93 %** | **8.05 %** |

Both CONFIRMED, and the 8.0 % (not 8.3 %) is confirmed. Neighbour thinning 140.0 → 119.7 under
N3 (−14.5 %), → 129.7 N4, → 139.2 of 143.0 under N3-tile (−2.6 %), → 133.5 under N3-var
(−4.6 %) ✔.

**Admissible-offset enumeration CONFIRMED:** frozen N3_occ per in-band section **66, 6, 1, 9, 10,
12** of **38,080–108,375** candidates (0.00092 %–0.173 %), median displacement **28.30 µm**;
**7001 admits only the identity** (1 of 108,375, displacement 0). N4_occ **1–13 of 720**.
Pre-C6: **1–63**, median **27.42 µm**, N4_occ 1–12 ✔.

**Calibration study — all 48 cells of `var_sim_calibration.csv` CONFIRMED**, including every
rounding in the pack's table. Derived claims:

- tiled torus, irregular window, s ≥ 0.05: **0.0800, 0.0825, 0.1175** → "0.080–0.118" ✔;
  s = 0.02 irregular is **0.0475** → "0.048" ✔.
- whole-window torus overall **0.0325–0.0775** ("0.033–0.078") ✔; irregular **0.0325–0.0725**
  ("0.033–0.073") ✔; at s = 0.30 **0.0725–0.0775** ("0.073–0.078") ✔.
- **tiling more liberal than whole-window in 7 of 8 cells, 1 exact tie** (irregular, s = 0.02,
  both 0.0475) ✔.
- irregular window is **74 %** of the bounding box, and the six real sections show
  `tissue_frac_bbox` **0.6576–0.8579** ("0.658–0.858") ✔.
- 100 sampling points, 199 shifts, 400 replications ✔; MC SE **0.0109** ✔.
- **RS_count range 0.033–0.060 ✔ — but see R2 for the "±0.010" clause.**
- **See R10 for 2.35 vs 2.36×.**

**Variance-correction validation CONFIRMED** (`var_variance_check.csv`): restricting to the
**12 of 16** cases with `n_i_max/n_i_min ≥ 2`, N3-var median slope **−0.4506** (range −0.585 …
−0.122) and N4-var **−0.4923** (−0.644 … −0.119), against the RS_count prediction of −0.5.

---

# 8. CONFIRMED — priority 7, the Moran result and its power bound

**Rank agreement, all four aggregations reproduced from `moran_vs_a7.csv` (132 field × section
rows, 12 fields, 11 sections):**

| aggregation | ρ | p |
|---|---|---|
| section-clustered **mean** per field, knn6 raw, 12 fields | **+0.8951** | **8.37e-05** |
| section-clustered mean, knn6 cell-type-centred | **+0.9441** | **3.93e-06** |
| **median** per field, knn6 raw | **+0.9231** | **1.86e-05** |
| **per-row**, no aggregation, 132 pairs | **+0.7104** | **1.43e-21** |
| 5 controls × 11 sections, 55 pairs | **+0.1548** | **0.259** |

The first two also reproduce independently from `moran_pooled.csv`. **§0.4's resolution is
CONFIRMED in full, including the diagnosis that +0.923 is the median-per-field aggregation.**
Every value positive; the falsification of "the two tests disagree" survives all four.

**Pooled Moran per family CONFIRMED** (`moran_pooled.csv`): `all_controls` **+0.04552
[+0.03018, +0.06087]**, clustered p **6.0e-05**, max per-section `p_rand` **9.01e-13**;
codeword **+0.04213** [+0.02813,+0.05614], 5.35e-05, 5.20e-12; probe **+0.00583**
[+0.00331,+0.00835], 4.23e-04, **0.537**; genomic **+0.00421** [+0.00219,+0.00623], 9.13e-04,
**0.802**; `neg_probe_rate` **+0.00475** [−0.00100,+0.01049], **0.0954**, **0.965**; Tier B
modules **+0.0852 to +0.2437**; `density_50um` **+0.9483** [+0.9362,+0.9605]. ✔

**Both "report against interest" disagreements CONFIRMED.** (a) Moran's pooled CI on the 40
probes **excludes** zero while A7's does not. (b) `genomic_control` is A7's **third-largest**
|amplitude| of the five (0.0744 > 0.0604 > **0.0307** > 0.0225 > 0.0113) and Moran's
**smallest** I of the five (0.00421 < 0.00475 < 0.00583).

**Sparsity-matched comparison CONFIRMED** (`moran_per_feature_summary.csv` +
`moran_verdict.txt`): probe 440 features / 21 counts / I −0.00012 vs 4,495 matched genes / 32
counts / **−0.00018**; codeword 6,699 / 7 / −0.00004 vs 1,179 / 18 / **−0.00012**; genomic 231 /
63 / −0.00029 vs 8,735 / 62 / **−0.00025**. Median Gene Expression feature carries **5,884.5**
counts per section against a probe's **21.0** → **280×** ✔ (`counts_median`, `Gene Expression`
row). Discordance **14 of 55** ✔.

**Power bound CONFIRMED, and the model identity verified rather than assumed.**
`moran_kernel_power.csv`, 22 rows = 11 sections × 2 calls:

| quantity | median | range |
|---|---|---|
| ΔI from the whole A7 gradient | **2.1994e-04** | 1.367e-06 – 2.558e-03 |
| ΔI as a fraction of observed control I | **0.826 %** | 0.0021 % – 6.08 % |
| SE(I) | **1.4984e-03** | 1.226e-03 – 2.254e-03 |
| **β Moran could resolve (ΔI = 2·SE)** | **0.3622 SD** | **0.3082 – 1.0697** |
| — `tierA_p95` only | **0.3381** | 0.3082 – 0.4175 |

I verified `dI = β_z²·Var(k)·I(k)` to max abs error **9.4e-17** and
`β_min = √(2·SE / (Var(k)·I(k)))` to **1.3e-14** across all 22 rows — the model in the pack is
the model in the file. 0.362 / 0.0744 = **4.87** ("5×" is a round-up). See **M5** and **M6**.

---

# 9. CONFIRMED — priority 8, the FPR, and which rate it is

**Definition confirmed from the fit-level file, not from the summary.** `frac_CI_excludes_zero`
in `a7_summary.csv` reproduces exactly as the fraction of the 165 per-fit **two-sided** 95 %
bootstrap CIs (`beta_<design>_lo > 0 or beta_<design>_hi < 0`) in
`a7_control_probe_fits.csv`:

| response | base | n6 | n5 | **n6n5** | n2 |
|---|---|---|---|---|---|
| `neg_control_codeword` | 0.430 | 0.327 | 0.121 | **0.091** | 0.400 |
| `all_controls` | 0.430 | 0.327 | 0.121 | **0.103** | 0.394 |
| `genomic_control` | 0.194 | 0.194 | 0.109 | **0.109** | 0.194 |
| `neg_control_probe` | 0.315 | 0.309 | 0.139 | **0.145** | 0.297 |
| `neg_probe_rate` | 0.212 | 0.212 | 0.158 | **0.164** | 0.212 |

**"9–16 % against a 5 % nominal; 9–15 % on the four count-based responses" — CONFIRMED**
(0.091–0.164 and 0.091–0.145). The pre-C6 multiset is **{0.091, 0.103, 0.109, 0.127, 0.164}**,
so "9–13 %" is the pre-C6 range and rests on the pre-C6 probe value 0.127 ✔. *(One nuance: the
pre-C6 responses are attached to different values — pre-C6 `all_controls` is 0.091 and
`neg_control_codeword` is 0.109, the reverse of frozen. Quoting the pre-C6 list as an ordered
per-response list is wrong even before the vintage issue.)*

**The estimator-vs-filter distinction — CONFIRMED, and this is the right resolution.**
The reportable-fit filter is one-sided on the naive design (`beta_naive > 0 AND
beta_base_lo > 0`, `summarize_phase3_c1.py:46-52`):

```
all_controls          naive-filter 0.1333   n6n5-filter 0.0485
genomic_control       naive-filter 0.0303   n6n5-filter 0.0485
neg_control_codeword  naive-filter 0.1333   n6n5-filter 0.0485
neg_control_probe     naive-filter 0.0848   n6n5-filter 0.0485
neg_probe_rate        naive-filter 0.0667   n6n5-filter 0.0485
```

**3.03 %–13.33 % on the naive design; 4.85 %, identical across all five, on the full design.**
So checklist item 22 is correct: the 9–16 % is the **estimator's two-sided CI-exclusion rate**
under the full N6+N5 design, and the **reportable filter is essentially nominal**. The two are
different statistics on different designs with different sidedness, and the pack keeps them
apart correctly. See **M2** for the 4.8 / 4.9 truncation.

---

# 10. CONFIRMED — everything else I checked

- **Cohort** (`composition_by_arm_timepoint.csv`): 11 sections / 11 mice; **6 SBR / 5 sham**;
  total **1,826,893**; analysable **1,635,937**; in-band six **1,031,880**; per-section
  **83,392–236,905**. All ✔. The §17 figures the pack calls wrong (1,834,806 / 1,036,459) are
  indeed in no file.
- **Window** (`window.csv`, 297 rows): median NN **6.74–10.61 µm** over 11 ✔; `tierA_p95`
  in-band `d_p99` **76.0–112.1**, median over 11 **96.1** ✔; `cdkn1a_pos` **76.3–160.9** ✔;
  `senepy_p95` **118.3–186.8** — **all six in-band values exceed the 100 µm cap** ✔.
  Discarded share: `tierA_p95` **0.75–7.04 %** beyond 80 µm, **0.00–0.11 %** beyond 150;
  `senepy_p95` **7.45–21.46 %** and **0.19–2.29 %** ✔. Sender prevalence `tierA_p95` in-band
  **4.04–4.48 %, mean 4.29 %** ✔.
- **Gene sets**: `A_SENDER_FINAL_strict.txt` = **33** lines ✔; Tier B line counts **126 / 68 /
  100 / 190 / 125 / 31 / 108** ✔, so `B_oxidative_stress` clears 30 by **one gene** ✔;
  18 of 21 Tier B pairs non-disjoint ✔ (see **R5** for the framing).
- **CoreScence circularity** (`figure_gs3_..._data.csv`): mouse **26/33 = 78.8 % → 29/33 =
  87.9 %**; human **25/33 = 75.8 % → 29/33 = 87.9 %**. The pack's "mouse 79 % → 88 %, human
  76 % → 88 %" ✔, and "69 % → 76 % → 88 %" is in no file ✔.
- **Figure 1 / §5.1** (`figures/figure1_data.csv`, 20 × 68): 8 of 20 cells at ℓ ≥ 2λ ✔;
  coverage naive **0.5125**, matched-decoy **0.3458**, nuisance **0.8542** ✔; over all 20
  **0.6700 / 0.6300 / 0.7017** ✔; naive iid **0.4100** ✔; decoy relbias better in **12 of 20**
  ✔; decoy coverage **strictly worse than naive in 5 of 8**, with **one exact tie** — so "6 of
  8" only if the tie counts ✔.
- **Phase 4 / §5.4** (`results/phase4/headline.csv`, `pair == 'ALL'`): every survival rate in
  the pack's table reproduces to three decimals for all four methods × five nulls. See **M4**.
- **Kernel families / §5.3** (`summary_phase5.txt` T4): naive AIC wins step **0.508**, spline
  0.248, exp 0.121, gauss 0.089, powerlaw 0.035 ✔; beats-no-kernel naive step 0.898, gauss
  0.778, exp 0.743 ✔. ctrl wins step **0.898**, spline 0.054, gauss 0.032, exp 0.013, powerlaw
  0.003 ✔; beats-no-kernel ctrl step **0.511**, gauss 0.241, exp 0.197, spline 0.152, powerlaw
  0.130 ✔. d̂½ spread ctrl **6.00×**, naive **3.20×** ✔. Held-out ctrl (252 folds): spline
  **−0.386** and step **−0.060** median ΔLL, both negative ✔ (see **M7**).
- **Superposition / §5.3** (T1): ctrl ΔAIC/1k median **−0.1089**, superposition wins **0.762**,
  paired win fraction median **0.730**, decisive ≥0.95 in **0.130**, decisive for nearest in
  **0.000** ✔; amplitude bound ctrl sup **+0.0076** IQR [−0.0163,+0.0185], nearest **−0.0029**;
  naive **+0.0277** / **+0.1472** ✔; neither beats covariates-only under ctrl (nearest
  **+2.67**, improves 0.20; sup **+0.74**, improves 0.45) ✔; synthetic correct verdict **1.00**
  in all four conditions ✔.
- **Winner's curse** (T2): placebo-corrected **+0.0513** at selection fraction 0.50 ✔ — and the
  pack does name the fraction, which matters because at 0.75 it is **+0.0739**. SF(N5) by
  sample size **0.1440 / 0.1998 / 0.2778** ✔.
- **D2 / `denoise`** (`d2_depth.csv`, `d2_agreement.csv`, `d2_normalisation_strength.csv`):
  depth loading 7259 **0.3176 → 0.5314**, 7352 **0.4096 → 0.5419**, 7239 **0.3891 → 0.6404** ✔;
  Pearson r **0.614 / 0.671 / 0.718** ✔; global top-5 % Jaccard **0.1412 / 0.1184 / 0.2804** ✔;
  committed set not called **75.26 / 78.82 / 56.21 %** ✔. `mor` removes **11.41–23.57 %** of
  log-depth variance over four sections ✔; `lib` removes **100 %** (p90/p10 11.0 → 1.000),
  moves the depth loading **−74.2 %** (7259) and **−92.9 %** (7352), changes **46.4 %** and
  **100 %** of sender calls ✔. Determinism floor r **0.99999913 / 0.99999995**, **24 of 75,384**
  and **2 of 114,721** cells changed ✔. **See R1 for the one number that is wrong.**
- **§0.2, the 76 % composition surrogate**: 0.212 / 0.260 = **0.8154**, not 0.76; the only two
  "ALL RECEIVERS" amplitudes on disk are **0.260** (`summary_phase5.txt` T3) and **0.266001**
  (`figures/figure2a_amplitudes.csv`); an implied denominator of **0.2789** appears in no file.
  **UNSOURCED — CONFIRMED as unsourced.**
- `references.bib`: `grep -c '^@'` = **43** ✔.
- Environment: cgroup ceiling `memory.max` = **61,999,996,928 B = 57.7 GiB** ✔.

---

# 11. WHAT THIS AUDIT DID NOT COVER

- The human arm in every form. `data/raw_h1/`, `results/phase9_h1/`, `data/processed_h1/` and
  `code/h1_*` were not read, per instruction. §2.0's H1 descriptors (7 donors, 2,207,593 cells,
  the 5,093-feature panel) are therefore **UNVERIFIED HERE**, not refuted — they are marked
  `[F]` in the pack and sourced to `PHASE7_H1_SCREEN.md`, itself a report.
- Every `[F]` value sourced to a **report** rather than a `results/` file. Those are one
  transcription away from `results/` by construction and are exactly where R3 was found.
- The bibliography's author lines. The **count** is confirmed at 43; the "41 wrong given names
  across 19 of 32 entries" finding is against a 32-entry vintage and cannot be restated for the
  43-entry file without a Crossref/PubMed pass. The pack says this and is right to.
- `results/phase3/n8_*.csv`, `stratification.csv`, `correlogram.csv`, `tierC_*`,
  `attribution.csv`, `lamscale*.csv`, the `sbrscope` variants, and the phase-4 supplementary
  sweeps.
- Whether the frozen tag's contents match the working tree. `git status` shows
  `code/h1_*` untracked (another agent's active work) and `figures/revised_candidates/`
  untracked; nothing under `results/phase3/`, `results/phase5/`, `results/moran/` or
  `results/phase8_d2/` is modified relative to HEAD.

---

# 12. THE FIX LIST, IN ORDER OF COST TO THE PAPER

| # | fix | where |
|---|---|---|
| 1 | **`PREREG_PHASE8.md` P3 instructs "9–13 %" off the pre-C6 probe value.** It is a pre-registered deviation and will be read as authority. Correct or mark it | `PREREG_PHASE8.md` P3 |
| 2 | `PREREG_PHASE8.md` **P2** carries pre-C6 digits *and* the false instruction "A7 must be re-run after 8.7"; **P24** carries "~80 %" against the frozen 86.3 % | `PREREG_PHASE8.md` P2, P24 |
| 3 | `NOVELTY_ASSESSMENT.md` §4 **O2** quotes `neg_probe_rate` at "+0.014, p = 0.079" (pre-C6); frozen is **+0.0113, p = 0.232** | `NOVELTY_ASSESSMENT.md` :582, :305 |
| 4 | Delete or replace the `figure4.png` md5; it is a pre-C6 hash from `CS_PHASE7_C1.md` | `WRITING_PACK.md` §5.4 |
| 5 | Fix "`denoise=False` Jaccard 0.76–0.99" → "r 0.996, Jaccard 0.76" | `WRITING_PACK.md` §5.9 |
| 6 | Drop "never outside ±0.010 of nominal" from the RS_count sentence | `WRITING_PACK.md` §5.7 |
| 7 | Restate the figure guard as covering 52, and name `figures/revised_candidates/` as outside it | `WRITING_PACK.md` §6 |
| 8 | "pooled CI half-width ±0.018 SD" → name the response/design, or give 0.016–0.056 | `WRITING_PACK.md` §5.5 |
| 9 | "the same 42 fits / 33 reportable" → add "per seed" | `WRITING_PACK.md` §5.6 |
| 10 | Intersection matrix: 224 cells, and scope "every Tier A × Tier B cell zero" to the PRIMARY strict set | `WRITING_PACK.md` appendix |
| 11 | Fix the four `run_phase3_nulls.py` line numbers (52 / 143 / 204 / 52–72) | `WRITING_PACK.md` §5.3, §0.5 |
| 12 | Remove the §5.7 clause implying 2.36 is exact; 2.35 is | `WRITING_PACK.md` §5.7 |
| 13 | Scope the λ̂ IQR sentence to the non-interior definitions | `WRITING_PACK.md` §0.1 |
| 14 | Name the subset on the Poisson slope/r²; state the Moran bound as a median with its range | `WRITING_PACK.md` §5.2, §5.5 |
| 15 | `COMPLETED_TASKS.md`:161's "Residual stale digits: 0" for `PREREG_PHASE8.md` is false | `COMPLETED_TASKS.md` |

---

*Every number in this file was recomputed from `results/`, `figures/` or `genesets/` in this
session. Nothing outside `reports/AUDIT_NUMBERS_FINAL.md` was written.*
