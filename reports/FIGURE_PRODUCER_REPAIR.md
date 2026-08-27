# FIGURE PRODUCER REPAIR

**Date:** 2026-08-27. **Scope:** the figure *producers* in `code/`, against
`reports/AUDIT_FIGURES.md`. **Guard:** `python3 code/check_figures_guard.py`
was `OK: all 52 committed figures match` before this work and is
`OK: all 52 committed figures match` after `--snapshot`; **23 of the 52
artefacts changed and were re-snapshotted, 29 are byte-identical.**

---

## 0. Lead — what changed and what was rebuilt from nothing

**Figures whose artefacts changed (23 files, 11 figures):** `figure2a`,
`figure2b`, `figure2c`, `figure2e`, `figure3`, `figure4`,
`figure4_supp_ncem_lengthscale`, `fig_phase3_caller_depth`,
`fig_phase3_tierC_identifiability`, plus the `figure3_data.csv`,
`figure4_data.csv`, `figure2a_amplitudes.csv` and two `fig_phase3_*_data.csv`
backing files.

**Figures whose artefacts did NOT change (29 files):** `figure1`, `figure2d`
(byte-identical after a full producer rebuild), `fig_phase3_composition`
(byte-identical, PNG and CSV), `figure4_supp_commot_mechanism`,
`figure_gs1`–`gs4`, `figure_phase8_callers`, `figure_phase8_d3`, and the three
orphaned `figures/figure2a_*.csv`.

**Producers rebuilt from nothing (no producer existed anywhere):**

| producer | what it now makes | faithfulness check |
|---|---|---|
| **`code/make_figure2e.py`** — NEW FILE | `figure2e.png/.pdf` + `results/phase3/figure2e_data.csv` | every row of the committed CSV reproduces from upstream **before** anything was redrawn: 11/11 variants on all geometry/SF columns to 1e-12; panel-h curve 20/20 bins, max \|Δ\| 1e-16, `n` exact |
| **`fig3()` → `figures/figure3_data.csv`** in `make_phase5_figs.py` | the six + two extra source frames the panels draw | see §4 |
| **`fig2b()` → `results/phase3/figure2b_data.csv`** in `make_figure2bc.py` | 1200 rows, 5 series | **1200/1200 rows reproduce, max \|Δ\| = 0.0 on every column** |
| **`fig2d()` → `results/phase3/figure2d_data.csv`** in `make_figure2bc.py` | 81 rows | **81/81 exact; `figure2d.png` regenerates BYTE-IDENTICAL** |
| **the three `fig_phase3_*` `.pdf` and `_data.csv`** in `make_phase3_figs.py` | added in `1351ce8` without touching the producer | `fig_phase3_composition.png` and its `_data.csv` regenerate **byte-identical**; tierC CSV identical except a 2.8e-14 float-associativity residual |

**Producers that disagreed with their own artefact and were corrected**
(re-running them regressed the figure; the guard cannot see this):
`make_figure4.py` (§1), `make_figure2bc.py` (§3), `make_phase3_figs.py` (§6),
`make_figure4_supp.py` (§5).

**What could not be done, stated plainly:** the frozen primary null cannot be
drawn as a **band** in `figure2b`. `results/phase3/perm_curves*.csv` contains no
variance-corrected draws at all, so no per-bin interval for N3-var/N4-var exists
anywhere in the repository. I did not substitute another null; the figure now
says this on its face and points to 2c for the surviving fraction.

---

## 1. PRIORITY 1 — `make_figure4.py` was regressing the figure to pre-C6

**Confirmed exactly as the audit describes, by running it.** With the committed
producer, `figure4_data.csv` came back with two rows changed and 136 unchanged:

```
c, our SASP kernel estimator, N3_lig   0.998871 (n=153)  ->  1.0001 (n=160)
c, our SASP kernel estimator, N4_lig   0.947193 (n=153)  ->  0.9641 (n=160)
```

`grep -rn load_ours code/` was empty and `git log -S load_ours -- code/make_figure4.py`
was empty: the de-hardcoded version documented in `CS_PHASE7_C1.md:18,295` and
`WRITING_PACK.md:488,1431` was never committed. It is re-implemented as
`load_ours()`, which rebuilds the reportable population from
`results/phase3/main_fits.csv` ∩ `perm_nulls.csv` using
`summarize_phase3.py:85`'s own definition verbatim — in-band sections,
`tierA_p95`, `stratum=all`, `beta_naive > 0` **and** `beta_base_lo > 0` — giving
153 fits, and recomputes all seven literals:

| literal | was (pre-C6, 160 fits) | now (153 fits) |
|---|---|---|
| `N3_sf` | 1.0001 | **0.998870982287264** |
| `N4_sf` | 0.9641 | **0.9471932849110392** |
| `N1_sf` | 0.7160 | **0.7073616928464024** |
| `N3_reject` | 0.875 | **0.8235294117647058** |
| `N4_reject` | 0.900 | **0.8496732026143791** |
| `beta_obs` | 0.01126 | **0.0101481078482851** |
| `N3_null_abs` | 8.91e-5 | **8.959122385730753e-05** |
| `n_interactions` | 160 | **153** |

**Verification demanded by the brief:** with `load_ours()` in place and *no other
change*, regenerating produced `figure4.png` and `figure4_data.csv`
**byte-identical to the committed artefacts** (guard: `OK: all 52`). So the
committed artefact was produced by exactly this data-driven logic, the artefact
was right, the producer was wrong, and the reconstruction is not a plausible
substitute — it is the thing. That state is commit `9f5ee5a`, which changes no
artefact at all.

`figures/revised_candidates/figure4_data_REVISED.csv` was **not** used; it
carries the older 1.0001/160 values, and its README's "max 2.34e-5" claim is
wrong by ~700× (the real move is 0.9641 → 0.9472).

`figure4.png` did change later, in commit `10462c9`, for the reasons in §2.

---

## 2. PRIORITY 3 — the frozen primary null now appears in every figure that reports nulls

`sf_summary_var.csv` PRIMARY: **N3-var 0.995966 [0.9754, 1.0074]**,
**N4-var 0.985068 [0.9581, 1.0033]**, n = 153 each. These were the only 2 of the
15 variants in that file that no figure plotted.

| figure | before | after |
|---|---|---|
| `figure2c` | 19 rows, no var | **21 rows**: 10 published, 9 in-tissue C1, and a third group below its own divider — N3-var / N4-var, labelled ★ PRIMARY, in a distinct colour |
| `figure2e` (e, f, g) | 11 variants, no var | **13**: N3-var after N3-snap, N4-var last, green, bold row labels |
| `figure4` panel (c) | 2 bars labelled "N3 torus shift" / "N4 rotation" with **no variant qualifier** | **4 bars**: `N3 bbox`, `N4 bbox`, `★N3-var`, `★N4-var`, each with its own tick label; `figure4_data.csv` gains `cond = N3_var_lig / N4_var_lig` |
| `figure2b` | 3 bands, no var | **still 3 bands** — a var band does not exist (see §0). The figure now states that and gives the SF instead |
| `figure3` panel (d) | "N3 torus-shift" per caller | **unchanged.** `results/phase5/super_nulls*.csv` carries `null ∈ {N1, N3}` only; there is no variance-corrected superposition-vs-nearest run, so this panel cannot gain the primary null without new computation. Flagged, not faked. |

If `sf_summary_var.csv` cannot supply the pair, `make_figure4.py` and
`make_figure2e.py` **print a warning / exit** rather than silently drawing the
bounding-box pair alone.

---

## 3. PRIORITY 2 — `figure2e` had no producer, and drew a withdrawn number

`code/make_figure2e.py` is new. It rebuilds `results/phase3/figure2e_data.csv`
from `sf_summary_var.csv`, `sf_summary.csv`, `null_destructiveness{,_var}.csv`,
`main_fits.csv`, `a7_summary.csv` and `a7_control_probe_{fits,curves}.csv`, then
draws the figure from that CSV. Nothing is typed in.

**(f) One line, three numbers → one number with its distribution.** The figure
drew **"median λ̂ 15.7 µm"** (the withdrawn interior-only median, 15.682110,
recorded as unsourced and circular in `COMPLETED_TASKS.md:165-169` and
`CS_PHASE8_TORUS_VAR.md:7-9`), while its written caption said **12.8 µm**. Both
are gone. The line is now **λ̂ = 14.7321271090776 µm** — median `lam_naive` over
the 315 primary in-band `tierA_p95` `stratum=all` fits, which
`summary_phase3.txt` prints as `medlam 14.7` — and it is drawn **with
"IQR [7.0, 50.0], 60 % railed"** in the panel, because 189 of those 315 fits sit
on a grid bound and a point estimate alone misrepresents that.

**(h) One scope clause over two populations → two scopes, both read from files.**
The green box said *"…six in-band sections, hepatocytes"* and then printed β/sd
numbers that come from `a7_control_probe_fits.csv`, which spans **11 sections and
9 cell types** (165 control fits, 1,155 biological). Both scopes are now computed
and printed separately, and the panel additionally carries the per-family naive
amplitudes its own "40 probes + 609 codewords + 21 genomic controls" implies
(probe **−0.0225**, codeword **−0.0604**, genomic **−0.0307**), which were in
`a7_summary.csv` but in no figure.

**Correction to the brief.** The brief says panel (h) "carries nine superseded
numbers." It does not — the *figure* already carried the current values
(−0.0744 / [−0.1306, −0.0182] / p=0.0145 / −0.0642 / p=0.0124 / +0.0053
[−0.0162, +0.0268] / p=0.595 / +0.2767 → +0.0310). The nine superseded numbers
are in the **written caption** at `CS_PHASE8_CALLERS.md:474-485`, which is what
audit §2.7 says. I fixed the caption (§7) rather than "fixing" a figure that was
already right.

---

## 4. PRIORITY 4 — `figure2a` mislabelled three modules

`make_phase5_figs.py:49-55` was a three-way rotation of the numbering in
`build_genesets.py:105-118` (B5 = `emt_ecm`, B6 = `oxidative_stress`,
B7 = `secondary_senescence`), which Master Plan §9:455-461 also gives.
Corrected. `ORDER` now follows B1..B7 as well, so the legend reads in order;
that reassigns per-module line colours in `figure2a` and marker colours / row
order within each cell type in `figure3(a)`. **No plotted value moves.**

**Did anything else inherit it?** No. I checked every other site:
`make_figure_genesets.py:36-38`, `gate_genesets_guard.py:52`,
`build_genesets_mouse_c6.py:146,247`, `gate_disjointness_human.py:61`,
`build_genesets_human.py:317` and `BIO_PHASE7_JobA.md:217` all carry the correct
mapping. Figure 2a's legend was the only place in the repository using the wrong
one — which is exactly why it went unnoticed.

**`figure3_data.csv`** had no producer and its `source` column covered only the
`call = tierA_p95` files, leaving 6 of panel (d)'s 12 real boxplots unbacked.
`fig3()` now writes it from the frames the panels actually draw:

| panel | source | committed | now |
|---|---|---|---|
| 3a | `combined_donor.csv` | 63 | **42** (`zone == "all"`; the 21 zonation rows are not drawn) |
| 3b | `kernel_families.csv` | 3150 | 3150 |
| 3c | `proximal_vs_downstream.csv` | 6 | 6 |
| 3d | `misspec.csv` | 1440 | **480** (`fit_family == "exponential"`, which is what panel d pivots on) |
| 3d | `super_section.csv` / `super_nulls.csv` | 945 / 3150 | 945 / 3150 |
| 3d | `super_section_cdkn1a_pos.csv` / `_senepy_p95.csv` | **absent** | **966 / 945** |
| 3d | `super_nulls_cdkn1a_pos.csv` / `_senepy_p95.csv` | **absent** | **3220 / 3150** |

`figure2a_amplitudes.csv` changes **only in the 16th significant digit**
(`0.26600088809815814` → `...803`). That is **not** caused by the relabel: I ran
the unmodified `HEAD` producer against the same cached curves and it also emits
`...803`, so the committed CSV was written by a different numpy/pandas vintage.
All ten amplitudes agree to 15 digits.

`make_phase5_figs.py` also now **prints** when it reuses the cached
`figure2a_stratified_curves.csv`, and takes `--rebuild`. Previously the reuse was
silent and correctness depended on `_m1_rerun_stage5.sh:22` deleting the file
first.

---

## 5. PRIORITY 5

### 5.1 `make_figure2bc.py` silent degradation — fixed first, as instructed

Lines 40, 105 and 108 read `perm_curves.csv` / `perm_nulls.csv` behind
`if os.path.exists(...) else pd.DataFrame()`, and the panels then dropped the
null band and every null row **with no warning**: a missing input produced a
complete-looking figure with the null absent. Every input now goes through
`_need()`, which exits naming the file and what it supplies. The same guard now
covers the `n8_scrambled_*.csv` glob, a null name that is present in the file but
carries no rows, and a null series with no finite values.

### 5.2 `figure2b` / `figure2c`: producer now emits what the artefacts show

The committed producer drew **1** band in 2b and **10** rows in 2c; the artefacts
show **3** and **19**. Restored from `perm_curves_c1.csv` / `perm_nulls_c1.csv`.
Verified before redrawing: `figure2b_data.csv` **1200/1200 rows, max \|Δ\| = 0.0
on every column, 5/5 series**; `figure2c_data.csv` **all 19 committed rows exact**
(then + 2 primary rows). The titles the artefacts carried and the producer did
not are now in the producer.

### 5.3 `figure4_supp_ncem_lengthscale`: stale against its own inputs

`make_figure4_supp.py:41-45` typed in the arrow endpoints `0.0104`/`0.0136` and
the label `"77%"`. Those matched a 2026-08-20 vintage of
`results/phase4/parts/`; the parts were regenerated 2026-08-21 and **575 of 720
r² values moved**, while the figure and its own rebuilt CSV stayed at the
previous evening — and because the producer overwrites
`results/phase4/ncem_radius_sweep.csv` unconditionally, that CSV was the only
surviving record of what had been drawn. The three numbers are now read off the
same frame the panel plots:

| | committed CSV | live `parts/` |
|---|---|---|
| median r² `real` | 0.01364 | **0.01364** (unchanged — only the nulls moved) |
| median r² `N0_perm` | 0.01044 | **0.01015** |
| annotation | 76.6 % → "77%" | **74.4 %** |

Regenerating also brings panel b onto the current parts: median selected radius
N3 30→20, N4 25→20, N3t 50→45, N0 45→50; `real` unchanged at 17.5.
`results/phase4/ncem_radius_sweep.csv` was rewritten (the producer's normal
behaviour). **`CS_PHASE4.md:710-722` still reproduces the stale table and now
disagrees with both the figure and the data** — not repaired here.

### 5.4 `fig_phase3_caller_depth`: pre-C6 table, and 2 of 11 sections

`make_phase3_figs.py:36` read `caller_within_type_depth_bias.csv`, which is
byte-identical to `results/phase3_pre_c6/`. I merged the committed
`fig_phase3_caller_depth_data.csv` against both candidates:

```
vs caller_within_type_depth_bias.csv          40 rows, max |Δ| enrichment 0.069
vs caller_within_type_depth_bias_2sec_c6.csv  40 rows, max |Δ| enrichment 0.0
```

So the artefact is the **C6** table and the producer read the pre-C6 one;
re-running it regressed e.g. SBR Tier A Q5 from 0.262 back to 0.311. The
producer now reads `_2sec_c6.csv`.

The base was **2 of 11** sections while `caller_within_type_depth_bias_11sections.csv`
(220 rows, 11 sections) sat unused. A **third panel** now plots it — one thin
line per section plus the across-section median per caller. The two committed
panels are unchanged. The backing CSV carries both bases with a `basis` column.

### 5.5 Tier C title

`fig_phase3_tierC_identifiability`'s title claimed "6 SBR sections x 14 Tier C
ligands" — 84 points — for a panel that draws **81**; three (section, ligand)
cells have a zero median distance and cannot go on a log axis. The counts are now
computed at draw time and the drop is stated.

---

## 6. Artefact hash table (sha256, first 12; PDFs with date stamps stripped)

| artefact | before | after |
|---|---|---|
| `fig_phase3_caller_depth.png` | `4de3ef2df3f2` | `d6148ad9092f` |
| `fig_phase3_caller_depth.pdf` | `a8eb0c3a42bb` | `186649e352e7` |
| `fig_phase3_caller_depth_data.csv` | `d68782b38030` | `d8668a8bd736` |
| `fig_phase3_tierC_identifiability.png` | `53a2d596c356` | `82a70088f770` |
| `fig_phase3_tierC_identifiability.pdf` | `9cce859dead7` | `1d53e58aed45` |
| `fig_phase3_tierC_identifiability_data.csv` | `428cdd019f96` | `9b9925ca51cb` |
| `figure2a.png` | `eb0d07191ddd` | `1cb1ddb9aa6a` |
| `figure2a.pdf` | `2efcc417b8eb` | `5664613c2a36` |
| `figure2a_amplitudes.csv` | `bb8cae2fe870` | `f75db4a7be6b` |
| `figure2b.png` | `a7e7bebd6f4e` | `11089d3e9f26` |
| `figure2b.pdf` | `fa0fba9a9d0d` | `87a67b0bbacd` |
| `figure2c.png` | `b36f4edd585c` | `0a81ac0d0f20` |
| `figure2c.pdf` | `045eb259edea` | `acaf2f030540` |
| `figure2e.png` | `d9c907ff4d3c` | `103f899be0a1` |
| `figure2e.pdf` | `69e6433b61a0` | `b4ccab2a560d` |
| `figure3.png` | `aa7ebb78c580` | `e0a4f81093c4` |
| `figure3.pdf` | `2c511faa1a0b` | `ac844f7b4734` |
| `figure3_data.csv` | `4cc2ed68b37b` | `97f013bd6bb4` |
| `figure4.png` | `76718f1ec5c1` | `f929f92b0563` |
| `figure4.pdf` | `67565e10e9e7` | `8c041f4b3733` |
| `figure4_data.csv` | `565fa8194b49` | `7e96ee247748` |
| `figure4_supp_ncem_lengthscale.png` | `895d6a03094a` | `a02ef255a680` |
| `figure4_supp_ncem_lengthscale.pdf` | `5dd9e24b2cc5` | `070f51d35ae0` |

Reasons, one line each: **2a** legend B5/B6/B7 + colour order; **2b** the
"no var band" note (bands and curves pixel-unchanged); **2c** two primary rows +
divider; **2e** rebuilt from a new producer (λ̂, primary null, panel-h scope);
**3** 3a colour/row order from `ORDER`; **3_data** panel-d callers now backed;
**4** four bars in panel c + the CIT-5 "implying" hedge in panel b; **4_supp**
current parts + computed annotation; **caller_depth** C6 table + 11-section
panel; **tierC** computed title; **amplitudes/tierC_data** float residuals only
(1e-16 / 2.8e-14).

Outside `figures/`, these results files were rewritten by their own producers:
`results/phase3/figure2b_data.csv`, `figure2c_data.csv` (+2 rows),
`figure2d_data.csv` (81/81 exact), `figure2e_data.csv` (λ̂ row corrected, var and
per-family rows added), and `results/phase4/ncem_radius_sweep.csv`.

**Guard:** before — `OK: all 52 committed figures match`. After the
regenerations and before `--snapshot` — `23 CHANGED, 29 match`, exit 1, exactly
the intended set. After `--snapshot` — `snapshot: 52 committed figures recorded`,
then `OK: all 52 committed figures match`, exit 0. Re-running every touched
producer a second time reproduces all 52 artefacts identically, so the pipeline
is idempotent.

---

## 7. Captions brought onto the figures they describe

- `CS_PHASE8_C1_CLOSEOUT.md` §4.2: (c) 160 → **153** reportable fits; (c) names
  the new primary group; (b) records why no var band exists; (f) λ̂ **14.7321 µm
  with IQR [7.0, 50.0] and 60 % railed**; (h) rewritten from *"Pending: …
  requires the H1 arm"* — which described a panel that has been filled with the
  M1 result — with the curve's and the amplitudes' scopes stated separately. The
  closing claim that pending status is marked "in the docstring of `fig2e()`" is
  retracted; no such function existed.
- `CS_PHASE8_CALLERS.md` §4.5 (the panel-h caption): all **nine** superseded
  numbers replaced with the `a7_summary.csv` values the figure has been drawing,
  the one-clause-two-populations scope error fixed, and the per-family breakdown
  added. Superseded text retained in a `<details>` block.
- `SASP_Kernel_Master_Plan.md` §25: the Figure 2 caption described **one** null
  band for a figure showing **three** and never mentioned the primary null
  (`PLAN_UPDATE_D12_D13.md:76-77`, still open); the Figure 4 caption still called
  a torus-shift failure *"a direct senescence-specific replication of the
  CellWHISPER result"* — the CIT-1 error `NOVELTY_ASSESSMENT.md:584` (O4) records
  as "NOT DONE. Verified error" — and listed three methods for a figure that
  plots four. Both corrected in place with a dated correction note.
- Inside `figure4` itself, panel b's *"CellWHISPER's reported >90% FPR"* now uses
  the **"implying"** construction that `Master_Plan.md:204` (CIT-5) mandates
  everywhere; this string was the last site in the repo without it.

---

## 8. Known-remaining, not repaired here

Honest list, so none of it is mistaken for done.

1. **`figure3` panel (d) has no primary null and cannot get one** without a new
   run: `results/phase5/super_nulls*.csv` carries `null ∈ {N1, N3}` only.
2. **`figure2b` cannot carry a primary-null band** — no `_var` draws in
   `perm_curves*.csv` (§0). Stated on the figure.
3. **`figure1` panel (b)** still plots three coverage heatmaps and not
   `cover_lam_nuis_blk`, while §25's correction says "the panels already show
   this" (audit §2.3). Adding a fourth heatmap changes a figure the audit
   certifies as consistent; left for an explicit decision.
4. **`figure4_supp_commot_mechanism`** still annotates *"Ccr2 is detected in
   ~1–3 % of cells"* with no file behind it, and `make_figure4_supp2.py:46`'s
   `ylim(0.9, 1.25)` still drops section 7248's `Ccl2→Ccr2 mass_ratio = 1.4599`
   off-canvas with no clip marker — the single largest violation of the panel's
   own title. Outside this task's five priorities; both are real defects.
5. **`make_figure_genesets.py:238,334`** still write *"through the pinned **1:1**
   MGI map"* into two figure captions; `AUDIT_CORRECTIONS_APPLIED.md` M2
   re-derived it as many-to-one (431 human symbols receive >1 mouse gene).
6. **`figure_phase8_callers` / `figure_phase8_d3`** still have colliding panel
   titles rendering text illegible, and the former's docstring promises 95 %
   intervals panel (b) never draws.
7. **Report-side staleness the audit lists** and I did not touch:
   `CS_PHASE4.md:710-722` (now disagrees with the regenerated ncem figure),
   `CS_PHASE4.md:730`, `BIO_PHASE3.md:26,192,333,405-411,419`,
   `CS_PHASE5.md:87-90`, `WRITING_PACK.md:287,493,893,896,1133-1134,1280,1513`,
   `PLAN_UPDATE_D12_D13.md:216,277`, `COMPLETED_TASKS.md:630`,
   `CS_PHASE3.md:527`, `AUDIT_PHASE8_FACTCHECK.md:471`.
8. **The stale-input caches** at `run_figure2a.py:227-228`, `sasp_phase3.py:99-101`
   and `caller_disagree_all.py:5-9` are unchanged. Only
   `make_phase5_figs.py:412-413` was made loud (§4).
9. **Guard scope.** `check_figures_guard.py` still watches only `figures/`, so
   `results/phase3/figure2{b,c,d,e}_data.csv` and
   `results/phase4/ncem_radius_sweep.csv` — which this work rewrote — remain
   outside it, and it still compares artefacts to themselves rather than to a
   producer. §5.3 is precisely the drift it cannot see. Extending it was not in
   scope.
10. **`figures/figure2a_curves.csv` and `figure2a_fits.csv`** are still orphaned
    provisional artefacts under the guard, written by `make_figure2.py`, which is
    not figure2a's producer.

---

## 9. Commits (none pushed; no tag created or moved)

```
9f5ee5a  Figure 4: read our own estimator's numbers from results/, not literals
39e7791  figure2e: write the producer that never existed, and fix what it drew
4f56a0d  figure2a: fix the B5/B6/B7 mislabel; give figure3_data.csv a producer
7a0ba99  make_figure2bc: stop degrading silently; restore the bands and rows; add the primary null
10462c9  Figure 4 panel c: draw the frozen primary null, and name every variant
5812690  figure4_supp_ncem: compute the panel-a annotation instead of typing it
a8b5641  fig_phase3_*: read the C6 table, plot all 11 sections, emit the PDFs and CSVs
1ee3a14  Bring figure 2e's two written captions onto the figure they describe
8678808  Master Plan §25: correct the Figure 2 and Figure 4 captions
```

`data/`, `genesets/`, `results/phase9_h1/`, `data/processed_h1/` and `code/h1_*`
were not touched. `figures/.committed_manifest.json` is gitignored, so the new
snapshot is local state, not a commit.
