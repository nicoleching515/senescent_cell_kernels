# AUDIT — figures against the data they claim to show

**Scope:** all 19 figures in `figures/` at `phase8-frozen`, traced to producer → `*_data.csv` →
`results/`, plus every caption that describes them.
**Date:** 2026-08-27. **Read-only.** No figure was regenerated, modified or deleted; no
`--snapshot` was run. `python3 code/check_figures_guard.py` → `OK: all 52 committed figures match`,
exit 0. Working tree `figures/` is clean against HEAD.

---

## 0. Verdict table

| Figure | Verdict |
|---|---|
| `figure2e` | **UNTRACEABLE + INCONSISTENT** — no producer anywhere; plots a superseded λ̂ and a wrong scope clause |
| `figure4` | **INCONSISTENT** — committed producer cannot reproduce the committed artefact; plots a non-primary null |
| `figure4_supp_ncem_lengthscale` | **INCONSISTENT** — stale against its own current inputs |
| `fig_phase3_caller_depth` | **INCONSISTENT** — producer reads a pre-C6 table the artefact does not use; 2-of-11-section base |
| `figure2c` | **INCONSISTENT** — committed producer emits 10 of the 19 plotted rows; omits the primary null |
| `figure2b` | **INCONSISTENT** — committed producer emits 1 of the 3 plotted bands; omits the primary null |
| `figure2a` | **INCONSISTENT** — legend mislabels B5/B6/B7; two orphaned provisional `figure2a_*.csv` in `figures/` |
| `figure3` | **INCONSISTENT** — `figure3_data.csv` has no producer and omits 2 of the 3 callers panel (d) plots |
| `figure1` | **CONSISTENT** (one caption over-claim: panel b does not show the nuisance-conditioned coverage the correction cites) |
| `figure2d` | CONSISTENT |
| `fig_phase3_composition` | CONSISTENT (data); `.pdf`/`_data.csv` have no producer |
| `fig_phase3_tierC_identifiability` | CONSISTENT (data); `.pdf`/`_data.csv` have no producer; two loose caption claims |
| `figure4_supp_commot_mechanism` | CONSISTENT (one untraceable annotation; one silently clipped point) |
| `figure_phase8_callers` | CONSISTENT (data exact); report §2.2 contradicts the figure; title collision |
| `figure_phase8_d3` | CONSISTENT (data exact); title collision |
| `figure_gs1_intersection_matrix` | CONSISTENT |
| `figure_gs2_crossarm_symmetry` | CONSISTENT (CDKN2B defect confirmed fixed) |
| `figure_gs3_corescence_circularity` | CONSISTENT (typed-in CoreScence literal confirmed removed) |
| `figure_gs4_senepy_coverage` | CONSISTENT |

---

## 1. Untraceable and inconsistent — lead with these

### 1.1 `figure2e` has no producer anywhere in the repo

`grep -rl figure2e /workspace --include=*.py` returns exactly one file: `code/check_figures_guard.py`.
There is no `fig2e()`, no `make_figure2e.py`, nothing. `reports/CS_PHASE8_C1_CLOSEOUT.md:305` cites
"the docstring of `fig2e()`" — that function does not exist. The figure cannot be regenerated,
audited by re-run, or corrected.

Its backing data (`results/phase3/figure2e_data.csv`) does exist and its 44 rows do trace upstream.

### 1.2 `figure2e` panel (f) draws a superseded λ̂ **in the figure itself**

`figure2e_data.csv` row `reference_median_lambda_hat_um = 15.682110`, and the PNG draws a black
vertical line labelled **"median λ̂ 15.7 µm"**.

The authoritative frozen value is **14.7321 µm** — `results/phase3/main_fits.csv`, median `lam_naive`
over the 315 in-band `tierA_p95` `stratum=all` fits; `summary_phase3.txt` §6 prints `medlam 14.7`.
I recomputed it: **14.7321**.

15.6821 is an *interior-only* median (a superseded definition that drops the 60 % railed fits) and is
the "15.7" that `reports/COMPLETED_TASKS.md:165-169` and `CS_PHASE8_TORUS_VAR.md:7-9` record as
having been **withdrawn as unsourced and circular**. The withdrawal reached eight documents; it did
not reach the figure. The line is load-bearing: panel (f)'s whole argument is that N3-occ/N4-occ move
senders "barely further than λ̂".

**And the figure's own written caption gives a third value.** `reports/CS_PHASE8_C1_CLOSEOUT.md`
§4.1-§4.2 — the panel-by-panel caption for this figure — says the reference line is
*"the median λ̂ (**12.8 µm**)"* at lines **195, 237 and 294**. So one drawn line carries three
numbers: **12.8** in its caption, **15.7** in the figure, **14.7321** in the frozen data. `12.8` is
the "FOURTH value then in circulation" that `COMPLETED_TASKS.md:24` corrected on 2026-08-27.

The same §4.2 caption is also stale on panel (h): it reads *"**(h)** *Pending:* the
negative-control-probe kernel (Section 13, test A7)"*, while the committed figure's panel h is filled
("M1 DONE, H1 pending") with the A7 result. Panel (e)'s caption numbers do check out
(N3-orig 0.772 → 22.8 %, N4-orig 0.920 → 8.0 %).

### 1.3 `figure2e` panel (h) states the wrong scope for the numbers it prints

The green box reads: *"A7: 40 probes + 609 codewords + 21 genomic controls, **six in-band sections,
hepatocytes**. β/sd, section-clustered mean [95% CI]: naive −0.0744 … Tier B modules, same fits:
+0.2767 → +0.0310."*

- The **curve** (red/green lines) is six in-band sections, hepatocytes. Verified: pooling
  `a7_control_probe_curves.csv` over `P.IN_BAND` reproduces `sf_median`/`sf_published` exactly
  (bin 2.5 µm: 0.034049 / 0.168863).
- The **β/sd numbers** are not. They come from `results/phase3/a7_summary.csv`, which reports
  `n_sections = 11` and `n_fits = 165` (all controls) / `1155` (biological). `a7_control_probe_fits.csv`
  spans **11 sections and 9 cell types**, not six sections and hepatocytes.

One scope clause is applied to two different populations. `AUDIT_PHASE8_FACTCHECK.md:494-496`
lists these panels as "Not audited."

Brackets are correctly typed here: `summarize_a7.py:44` computes a genuine section-clustered 95 %
t-interval, so panel (h)'s "[95% CI]" is right, and panel (g) correctly says "line = IQR". The
IQR-vs-CI misnomer recorded in `COMPLETED_TASKS.md` applies to `summarize_phase3.py:99`
(a plain `np.quantile([.25,.5,.75])`), and no figure mislabels those. Note the asymmetry a reader
must not collapse: **Figure 3's λ brackets are genuine donor-bootstrap CIs** (panels a and c), while
the SF and amplitude brackets quoted beside them everywhere else are IQRs across fits. No λ bracket
anywhere in the repo is labelled a confidence interval; eight documents say so explicitly.

Adjacent: the brief's A7 per-family frozen values **−0.0604 / −0.0307 / −0.0225**
(`neg_control_codeword` / `genomic_control` / `neg_control_probe`, `design=base`) are all present in
`a7_summary.csv` and all verified — but **none of them is in `figure2e_data.csv` or in the figure**,
which carries only the pooled `all_controls` −0.0744. The per-family breakdown the caption's
"40 probes + 609 codewords + 21 genomic controls" implies is not plotted.

### 1.4 `figure4`: the committed producer cannot reproduce the committed figure

`code/make_figure4.py:41-44` carries a typed-in dict:

```
# CS_PHASE3 §0 / results/phase3/perm_nulls.csv, the 160 reportable fits
OURS = dict(N3_sf=1.0001, N4_sf=0.9641, N1_sf=0.7160,
            N3_reject=0.875, N4_reject=0.900,
            beta_obs=0.01126, N3_null_abs=8.91e-5)
```

`N3_sf`/`N4_sf` are plotted; `n_interactions=160` is written at line 103. The **committed**
`figures/figure4_data.csv` instead ends:

```
c,our SASP kernel estimator (CS_PHASE3),ALL,N3_lig,...,0.998870982287264,153
c,our SASP kernel estimator (CS_PHASE3),ALL,N4_lig,...,0.9471932849110392,153
```

Those are `results/phase3/sf_summary.csv` PRIMARY on the current **153**-fit basis, independently
recomputed from `main_fits.csv` ∩ `perm_nulls.csv`. The PNG's bars are annotated **1.00 / 0.95**;
the producer's 0.9641 would print 0.96. 136 of 138 rows reproduce; these 2 do not.

All seven literals are stale (`N1_sf` 0.7160→0.7074, `N3_reject` 0.875→0.8235, `N4_reject`
0.900→0.8497, `beta_obs` 0.01126→0.010148, `N3_null_abs` 8.91e-5→8.96e-5).
`reports/CS_PHASE8_TORUS_VAR.md:76`: *"The reportable population is now 153, not 160."*
`grep -rn load_ours code/` → zero hits: the de-hardcoded producer documented in
`CS_PHASE7_C1.md:18,295` and `WRITING_PACK.md:488,1431` **is not in the repo**.

**Re-running `code/make_figure4.py` today silently regresses the figure to the pre-C6 160-fit basis.**

`figures/revised_candidates/figure4_data_REVISED.csv` carries the **older** 1.0001/160 values despite
its name; `revised_candidates/README.md`'s claim that the two differ by "max 2.34e-5" is wrong by
~700× (the real move is 0.9641→0.9472, visible at the 2-dp bar label).

### 1.5 `figure4_supp_ncem_lengthscale` is stale against its own inputs

`code/make_figure4_supp.py:26-29` globs `results/phase4/parts/ncem_sweep_*.csv` and overwrites
`results/phase4/ncem_radius_sweep.csv`. The parts are dated **2026-08-21 03:10**; the CSV and the
PNG/PDF are **2026-08-20 21:40** — the inputs were regenerated after the figure. 575 of 720 r²
values differ. Median selected radius, plotted vs current parts: N3 **30 → 20**, N4 **25 → 20**,
N3t **50 → 45**, N0 **45 → 50**. The panel-a arrow literals `0.0104 → 0.0136`, "**77%**"
(lines 41-45) recompute to 0.01015/0.01364 = **74.4 %** from the live parts.

Every plotted number traces to a file — but only to the producer's **own output**, which is now the
sole surviving record of what was drawn. `CS_PHASE4.md:710-722` reproduces the stale table, so
report and figure agree with each other and not with the data.

### 1.6 `figure2b` / `figure2c`: the committed producers emit a fraction of what is plotted

`code/make_figure2bc.py` draws **one** null band in 2b (`per.null == "N3"`) and **ten** rows in 2c.
The committed artefacts show **three** bands and **nineteen** rows, under titles the producer does
not contain ("…bounding-box and in-tissue (correction C1)"; "…with the corrected in-tissue N3/N4").

`results/phase3/figure2b_data.csv` carries five series (`observed`, `matched_decoy_N2`, `null_N3`,
`null_N3_snap`, `null_N3_occ15`); `figure2c_data.csv` carries 19 rows. The producer writes
`figure2c_data.csv` only — nothing in `code/` writes `figure2b_data.csv` or `figure2d_data.csv`,
contradicting `CS_PHASE8_M1_RERUN.md:128` ("`make_figure2bc.py` writes `figure2b_data.csv` and
`figure2d_data.csv`").

All plotted values themselves verify **exactly**: observed/decoy/N3 vs `curves.csv` and
`perm_curves.csv` max |Δ| = 0.0; snap/occ15 vs `perm_curves_c1.csv` max |Δ| = 0.0; 2c's first
10 rows reproduce from the producer's own logic to 9 dp (n = 153 of 315 confirmed); the 9 variant
rows equal `sf_summary_var.csv`.

Silent-degradation guards in the same file (`:40`, `:105`, `:108`) mean a missing `perm_curves.csv`
or `perm_nulls.csv` produces a figure with **no null band and no null rows, and no warning**.

### 1.7 `fig_phase3_caller_depth`: producer reads the pre-C6 table, artefact uses the C6 one

10 of 40 values in `fig_phase3_caller_depth_data.csv` (all `tierA_score`, max |Δ| 0.069) are absent
from `results/phase3/caller_within_type_depth_bias.csv`, the file
`code/make_phase3_figs.py:36` reads. They match `caller_within_type_depth_bias_2sec_c6.csv` exactly
(40/40). `caller_within_type_depth_bias.csv` is byte-identical to
`results/phase3_pre_c6/caller_within_type_depth_bias.csv`. Pixel readout of the PNG (SBR Tier A Q5 =
0.263) confirms the **C6** value 0.262, not the pre-C6 0.311. Re-running the committed producer
regresses the figure.

The base is **2 of 11 sections** (7250 sham, 7259 SBR — `caller_disagree.py:22`). The 11-section
table `caller_within_type_depth_bias_11sections.csv` exists and is unused.

### 1.8 `figure3_data.csv` has no producer and does not cover panel (d)

Nothing in `code/` writes `figures/figure3_data.csv` (`make_phase5_figs.py` writes PNG/PDF only),
contradicting `CS_PHASE8_M1_RERUN.md:129,230`. Its `source` column is only
`combined_donor.csv` / `kernel_families.csv` / `proximal_vs_downstream.csv` / `super_section.csv` /
`super_nulls.csv` / `misspec.csv` — all `call = tierA_p95`. Panel (d) plots **three** sender callers
(tierA_p95, Cdkn1a+, SenePy p95), so 6 of its 12 real boxplots have no rows in the backing CSV;
`super_{section,nulls}_{cdkn1a_pos,senepy_p95}.csv` are read by the producer and absent from the CSV.
Panel 3a also carries 63 rows when 42 are plotted (the 21 zonation rows are not in the figure).

Everything present verifies exactly (3a vs `combined_donor.csv` max |Δ| = 0.0; 3c vs
`proximal_vs_downstream.csv` identical; 3b row count matches `kernel_families.csv`). All three
panel titles recompute exactly: **34 of 42** CIs span [7,50]; step wins **89.8 %** / beats no kernel
**51.1 %**; **2 of 6** ratio CIs reach the grid bound.

### 1.9 The `fig_phase3_*` `.pdf` and `_data.csv` have no producer

`code/make_phase3_figs.py` writes three PNGs and nothing else (lines 33, 61, 91). The three PDFs and
three `_data.csv` were added in commit `1351ce8` without touching the producer, and carry columns
the producer never computes (`tierA_basis`, `is_section9_control`, `poisson_prediction_um`).
`CS_PHASE8_M1_RERUN.md:616` asserts a change the committed script does not contain.

---

## 2. Wrong captions on correct figures

### 2.1 Figure 2a's legend mislabels three of the seven Tier B modules

`code/build_genesets.py:105-118` is authoritative: **B5 = emt_ecm, B6 = oxidative_stress,
B7 = secondary_senescence**. Master Plan §9 (lines 455-461) and every report agree.

`code/make_phase5_figs.py:49-55` labels them:

```
"secondary_senescence": "2nd senescence (B5)",   # is B7
"emt_ecm":              "EMT/ECM (B6)",          # is B5
"oxidative_stress":     "Oxidative stress (B7)"  # is B6
```

These strings are the figure legend and the per-panel line labels (lines 167, 205), so **Figure 2a's
legend is the only place in the repo using this numbering**. The consequence is direct: the module a
reader identifies as "B7" in Figure 2a is oxidative stress (31 genes); the re-sourced B7
(`secondary_senescence`, 108 genes — the frozen configuration's headline module rebuild) appears as
"B5".

### 2.2 §25's Figure 4 caption still miscites CellWHISPER

`SASP_Kernel_Master_Plan.md:1092`, verbatim:

> *If they fail, that is a direct senescence-specific replication of the CellWHISPER result and it
> justifies the whole paper.*

Figure 4's bars are the **torus shift**. The CIT-1 correction is applied at §22 Step 3 (`:1025`),
§23 (`:1056`), §31 item 27 (`:1264`), §32 item 3 (`:1286`) and in `CITATION_AUDIT.md` — all state
CellWHISPER's null is a **within-cell-type permutation (our N1/N0_type)**, that the strings *torus*
and *toroidal* do not occur in their paper, and that *"the torus-shifted run reproduces nothing of
CellWHISPER's."* §25's caption was missed by that sweep.
`reports/NOVELTY_ASSESSMENT.md:584` (O4): **"NOT DONE. Verified error."**

The same caption lists three methods; the figure plots **four** (NCEM linear* is absent from the
caption). Inside the figure the fix *is* applied correctly (`make_figure4.py:143-149`), minus the
mandated "implying" hedge for the >90 % FPR (`Master_Plan.md:204`).

### 2.3 §25's Figure 1 correction over-claims what the panels show

The 2026-08-27 correction (`:1081`) says *"the panels already show this."* Its numbers are exact —
I reproduced all of them from `figures/figure1_data.csv`: **12 of 20** cells better, **8** worse,
worst 2.2699 vs naive 2.0250; over the 8 cells with ℓ/λ ≥ 2, coverage **0.5125 naive / 0.3458
matched-decoy / 0.8542 nuisance-conditioned**.

But the surviving clause is *"nuisance conditioning restores approximate calibration"*, and
`make_figure1.py:120-122` plots only three coverage heatmaps: naive-iid, naive-block,
matched-decoy-block. **`cover_lam_nuis_blk` is in `figure1_data.csv` and is not in panel (b).** The
0.85 is in the data; it is not in the figure. (Panel (a) *does* show all three for bias.)

### 2.4 §25's Figure 2 caption is two configurations out of date

`:1084` still reads *"(b) Same curve overlaid with the matched-decoy curve and the torus-shift null
band."* The committed figure2b shows **three** bands (bounding-box, snapped, occupancy-screened),
and the frozen primary null is N3-var, which appears in neither. Flagged at
`PLAN_UPDATE_D12_D13.md:76-77` and still open.

### 2.5 §25 captions 4 of the 19 committed figures

`grep -c "figure_gs\|figure_phase8\|fig_phase3\|figure2d\|figure2e\|figure4_supp" SASP_Kernel_Master_Plan.md`
→ **0**. §25 says "Four main figures" and describes Figures 1–4 only. figure2d, figure2e,
figure_gs1–gs4, figure_phase8_callers, figure_phase8_d3, fig_phase3_×3 and figure4_supp_×2 —
12 of 19 — are captioned only inside reports.

### 2.6 `figure4_supp_commot_mechanism` carries an untraceable annotation and clips its worst point

`make_figure4_supp2.py:68-70` annotates *"Ccr2 is detected in ~1–3 % of cells"*. No file in
`results/` carries a Ccr2 detection fraction. (The companion "within 100 µm" is traceable to
`phase4_commot_mechanism.py:30`, `dis_thr=100.`)

`:46` sets `ax.set_ylim(0.9, 1.25)`; section 7248's `Ccl2→Ccr2 mass_ratio = 1.4599` is drawn
off-canvas with no clip marker. The PNG shows 5 Ccl2 points, not 6 — and the dropped point is the
single largest violation of the panel's own title, *"but the transported mass is conserved."*

### 2.7 `figure2e` panel (h): every number in its written caption is superseded

`reports/CS_PHASE8_CALLERS.md:474-485` is the caption for this panel. Against
`results/phase3/a7_summary.csv` (which the figure plots):

| caption | figure / `a7_summary.csv` |
|---|---|
| amplitude **−0.070** SD, 95 % CI **[−0.128, −0.012]**, **p = 0.023** | **−0.0744**, **[−0.1306, −0.0182]**, **p = 0.0145** |
| matched decoy **−0.061**, **p = 0.020** | **−0.0642**, **p = 0.0124** |
| N5 **+0.007 [−0.011, +0.025]**, **p = 0.41** | **+0.0038 [−0.0186, +0.0261]**, **p = 0.715** (N6+N5 **+0.0053 [−0.0162, +0.0268]**, p = 0.595 — the value the figure prints) |
| naive Tier B **+0.291** SD | **+0.2767** |
| conditioned Tier B **+0.036** SD | **+0.0310** |

Nine numbers, nine superseded. `CORRECTIONS.md:1628` records the +0.291/+0.036 → +0.2767/+0.0310
move. This caption is also where the figure's wrong "six Test-3-admissible sections, hepatocytes"
scope clause (§1.3) originates.

### 2.8 `figure2c`'s written caption says 160 fits; the figure says 153

`CS_PHASE8_C1_CLOSEOUT.md:284` — *"**(c)** Surviving fraction of β̂ under each null over the **160**
reportable fits"* — and `CS_PHASE3.md:527` — *"160 fits per null"*. The committed figure's own title
reads **"n = 153 of 315 fits"**, `figure2c_data.csv` carries `n = 153` on every row, and I
reproduced 153 from `main_fits.csv` ∩ `perm_nulls.csv`. `CS_PHASE8_TORUS_VAR.md:76`: *"The reportable
population is now 153, not 160."* Same 160 is frozen into `make_figure4.py`'s comment and its
`n_interactions=160` (§1.4).

### 2.9 `figure_gs2` and `figure_gs3` render "the pinned 1:1 MGI map" into the figure

`make_figure_genesets.py:238` and `:334` both write *"through the pinned **1:1** MGI map"* into the
caption block. `AUDIT_CORRECTIONS_APPLIED.md:184-188` (M2) re-derived the map as **many-to-one**:
18,782 rows onto 17,609 distinct human symbols, **431 human symbols receiving more than one mouse
gene** — and marks M2 **FIXED**. The fix reached `PREREG_PHASE8_genesets.md` §6; the two figure
strings were not changed. (Everything else in these two figures verifies exactly.)

### 2.10 `figure_phase8_callers`' docstring promises intervals the code never draws

`make_figure_phase8_callers.py:6` — *"(b) the same pooled, two-section base vs eleven-section base,
**with 95% intervals**"*. Lines 84-99 plot three markers and a straight connector; no interval is
drawn anywhere in panel (b). The same file writes the `_data.csv`, so the docstring is the
provenance note a reader gets.

### 2.11 Report captions contradicted by the figures they describe

- **`CS_PHASE8_CALLERS.md:203-205`** — *"Tier A vs SenePy is the one pair that does not move:
  0.935 → 0.914, below chance on all eleven sections"* — is falsified by
  `figure_phase8_callers` panel (a), which plots the **frozen** basis: 0.972, z = −1.63 (not
  significant), 4 of 11 sections above 1.0, range 0.751–2.198 (report says 0.70–1.71).
  §2.2's whole "11-section pooled" column is the pre-C6 basis; the figure is post-C6.
  `caller_coverage_gate.csv` carries both.
- **`BIO_PHASE3.md:333`** tabulates the pre-C6 Tier A depth row (`0.27 / 0.17`) directly above the
  figure that now plots `0.269 / 0.146`. `CORRECTIONS.md:325-327` flagged the move; the update never
  reached this file.
- **`BIO_PHASE3.md:26, 405-411, 419`** quote the two-section agreement ratios (`2.15`, `0.38`,
  `2.85`, `1.51`, `0.93–1.22×`) that `SUBMISSION_PATCH_2026-08-29.md:88-90` withdrew in favour of
  the 11-section pooled **1.212×**. That patch's own instruction — *"any caption quoting an
  agreement ratio is not [fine] … Check the three captions"* — has not been applied here.
- **`WRITING_PACK.md:1133-1134`** cites `figure_gs1_intersection_matrix_data.csv` as *"196 cells,
  every Tier A × Tier B cell zero"*. The CSV has **224** rows, and only the gate cells are zero
  (e.g. `A_sender_for_tnfa_nfkb_proximal × B4 = 43` human, 35 mouse). This revives the exact error
  `BIO_PHASE8_FREEZE.md:218-222` records as caught and corrected in the figure itself. Same "196" at
  `AUDIT_PHASE8_FACTCHECK.md:471`.
- **`WRITING_PACK.md:287`** reports cross-arm B7 as "88 of 108/116 on the ortholog-intersected
  panel" — the conflation `figure_gs2`'s own caption explicitly warns against (88 is the
  map-gap-corrected shared count; the panel quantity is 85 mouse / 88 human, coinciding only for
  human).
- **`WRITING_PACK.md:493`** heads its Figure 4 column `N0_perm (=N1)`, reinstating the equivalence
  `CS_PHASE4.md:234-243` corrected under D7 §B9 (CellWHISPER's null is `N0_type`).
- **`CS_PHASE4.md:730`** writes *"1.000 [0.992, 1.008]"* with **no bracket label**, against
  `RECORD_RECONCILIATION.md` rule 2; the sentence is also stale on all four numbers
  (160 → 153 fits; 1.000 [0.992, 1.008] → 0.999 [0.989, 1.006]; 87.5 % → 82.4 %; 1.13e-2 → 1.01e-2).
- **`WRITING_PACK.md:893, 896, 1280, 1513`** still carry the tiled-torus inflation as **2.4×**;
  `CORRECTIONS.md:1878` standardises repo-wide on **2.35×** (0.1175/0.05 exactly), and 896/1513
  derive "2.36" from a pre-rounded 0.118 — the exact error `CORRECTIONS.md:131` names. No figure
  quotes this number.
- **`fig_phase3_tierC_identifiability`**'s own title says *"6 SBR sections × 14 Tier C ligands"*
  (84 points) for **81** points — `Cxcl12` is missing from three sections. `BIO_PHASE3.md:192` claims
  the Tier C outputs cover "6 SBR sections plus the 26 wk sham";
  `tierC_ligand_identifiability.csv` contains no sham section.
- **`CS_PHASE5.md:87-90`** describes Figure 3(b) as *"step wins **90.8 %** of fits by AIC and still
  beats the no-kernel model in only **55.9 %**"*, with family choice moving d̂½ by **4.4×**. The
  committed figure's runtime-computed title says **90 % / 51 %**, and I recomputed **89.8 % / 51.1 %**
  from `results/phase5/kernel_families.csv`. `WRITING_PACK.md:459-461` has the current 0.898/0.511.
  Because that title is computed at draw time, whichever `kernel_families.csv` is on disk silently
  wins — the report is the only place the old vintage is visible.
- **`Master_Plan.md:1081`** prints the naive worst-case relative bias as **2.03**;
  `CORRECTIONS.md:961-964` says the file value is exactly **2.025** and instructs *"Quote 2.025."*
  `COMPLETED_TASKS.md` prints 2.02. Three roundings of one number, one of them inside a caption.
  (The figure itself prints `+202%`, which is right.)
- **λ̂'s mandatory IQR companion has two upper bounds in circulation:** `[7.0, 50.0] µm, 60 % railed`
  at 15 sites (frozen), against `12.8 µm, IQR 7.0–48.2` at `CS_PHASE7_C1.md:232` (pre-C6 — the same
  line that feeds `figure2e`'s 12.8 caption, §1.2).
- **The tiled-torus 2.4× survives uncorrected at four sites** — `PLAN_UPDATE_D12_D13.md:216`, `:277`;
  `WRITING_PACK.md:1280`; `COMPLETED_TASKS.md:630` — against the repo-wide **2.35×**
  (`CORRECTIONS.md:1878`). No figure quotes it.
- **`figure2a`'s amplitude has two emitted values:** `figure2a_amplitudes.csv` **0.266001** and
  `summary_phase5.txt` T3 **0.260**; the withdrawn "76 %" composition share needs a denominator of
  ≈0.279 that no file emits (`CS_PHASE5.md:73-84`, `WRITING_PACK.md:117-130`; sourced ratio **81.5 %**).

---

## 3. Frozen configuration

Frozen = strict-33 Tier A · re-sourced B7 · 11-section DeepScence coverage · N3-var/N4-var primary.

**Strict-33 Tier A and re-sourced B7 landed.** `sf_summary_var.csv` labels its PRIMARY subset
*"in-band sections, tierA_p95, A_SENDER_FINAL_strict (33 genes)"*. Mouse/human `.txt` gene counts
(A 33/33; B1–B7 126/68/100/190/125/31/108 and 120/71/126/231/113/36/116) match `figure_gs1`/`gs2`
exactly. `genesets/human/B_secondary_senescence.txt` is set-identical to the
`_C6_sourced_ported` variant from `rebuild_b7_secondary_senescence.py`
(`results/phase7_jobA/b7_c6_rebuild.csv`: n=116, a_in_b7=0, PASS). `data/processed/cache3/*.npz`
(05:54) postdates the genesets (05:41) and the sender tables (05:48–05:53), so the phase-3 chain is
current.

**N3-var / N4-var — the primary null — is plotted by no figure.**
`results/phase3/sf_summary_var.csv` carries `N3_var 0.995966 [0.9754, 1.0074]` and
`N4_var 0.985068 [0.9581, 1.0033]`, and `perm_nulls_var.csv` / `perm_draws_var*.csv` back them.

- `figure2b` bands: N3 bbox, N3-snap, N3-occ15. No var.
- `figure2c` rows: N3/N4 orig, tile, occ, occ15, swap, snap. **N3-var and N4-var are the only two of
  the 15 variants in `sf_summary_var.csv` that were left out.**
- `figure2e` bars: same 11 variants, no var.
- `figure3` panel (d): "N3 torus-shift" per caller. No var.
- `figure4` panel (c): the bbox N3/N4 (0.9989 / 0.9472) vs the primary N3-var/N4-var
  (0.9960 / 0.9851), labelled "N3 torus shift" / "N4 rotation" with no variant qualifier, so a
  reader cannot tell which was used.
- `perm_curves*.csv` contains no var draws at all, so a **band** for the primary null does not exist.

`NOVELTY_ASSESSMENT.md:584` (O3) still records the variance-corrected run as "NOT DONE"; it has in
fact been run — it simply never reached a figure.

**11-section coverage landed unevenly.** `fig_phase3_composition` (11), `figure2d` (11 × 7 = 77
combinations), `figure2e` panel h's β/sd (11) are on 11 sections. `figure2a`, `figure2b`, `figure2c`,
`figure2e` panels e–g and `figure3` are on the six Test-3-admissible sections (correct and stated).
`fig_phase3_caller_depth` is on **2**, and the 11-section table exists unused.
`figure_gs4` is SenePy-only and carries no section count — "11-section DeepScence coverage" does not
apply to it.

**The three `fig_phase3_*` re-points:** `composition` and `tierC` never read a caller table and were
already on-config. `caller_depth` landed in the artefact but not in the producer (§1.7), and its
2-section base means the coverage leg never landed at all.
`SUBMISSION_PATCH_2026-08-29.md:573-579` says all three read caller tables; only one ever did, and
`CS_PHASE8_M1_RERUN.md:588-598` states the correct version — the two report families contradict each
other in the record.

---

## 4. Stale-input caching, beyond `figure2a`

| site | pattern |
|---|---|
| `make_phase5_figs.py:412-413` | the known one — `figure2a_stratified_curves.csv` read if present, rebuilt only if absent. Still live in the committed code; `_m1_rerun_stage5.sh:22` `rm -f`s it as a workaround, so correctness depends on running the wrapper |
| `run_figure2a.py:227-228` | per-`(section, module)` checkpoint skip on `results/real/fit_*.csv`. Those 70 files are frozen at 2026-08-20 16:33 and would never rebuild on a changed sender definition |
| `sasp_phase3.py:99-101` | `cache3/{sample}.npz` rebuilt only with `--force`. Feeds `figure2a`, `figure2b`, `figure2c`, `figure2e`, `figure3`. Currently fresh (05:54 > genesets 05:41), but a geneset edit without `--force` is invisible |
| `caller_disagree_all.py:5-9` | the same effect by naming convention: refreshed results land under a `_11sections` suffix while the figure keeps reading the frozen filename. `SUBMISSION_PATCH_2026-08-29.md:576-578` names this as *"they regenerated byte-identically, which looked like a passing reproducibility check and was a stale-input result"* |
| `make_figure4_supp.py:26-29` | not a cache — it rebuilds unconditionally — but it is the only record of its own plotted state, which is how §1.5 happened |
| `make_figure2bc.py:40,105,108` | silent-degradation guards: a missing null file yields a complete-looking figure with the null silently absent |
| **`figures/figure2a_curves.csv`, `figures/figure2a_fits.csv`** | **orphaned superseded artefacts under the guard.** Written by `make_figure2.py`, which is *not* figure2a's producer (`WRITING_PACK.md:1429`: "**never `code/make_figure2.py`**"). Their own provenance columns read `sender_source = "provisional:Cdkn1a>0 (Bio senders_*.csv pending)"`, `module_source = "provisional:z-mean of Tier B sets"`, `celltype_source = "provisional:unstratified"`, over **5 sections that are not the six in-band ones** (7239, 7250, 7361 in; 7260, 7001, 7248, 7435 out). Anyone tracing "figure2a" to a `figure2a_*.csv` in `figures/` lands here |

---

## 5. Guard scope

`check_figures_guard.py` verifies **52 artefacts, all inside `figures/`** — verified OK, and the
manifest covers every `.png/.pdf/.csv` on disk with nothing unwatched. But the backing data for
figures 2b, 2c, 2d, 2e, `figure_phase8_callers` and `figure_phase8_d3` lives in
`results/phase3/figure*_data.csv`, and `figure4_supp_ncem`'s lives in
`results/phase4/ncem_radius_sweep.csv` — **all outside the guard**. §1.5 is precisely a drift the
guard cannot see. The guard also compares artefacts to themselves, never to a producer, which is why
§1.4, §1.6 and §1.7 (producer ≠ artefact) all pass it.

---

## 6. What verified clean

- **`figure1`** — `figure1_data.csv` reproduces from `results/sweep_all.csv` bit-exactly
  (`max abs diff 0.0`, 600 main rows, 30 seeds/cell). All 120 heatmap cells in panels (a) and (b)
  match the CSV. Panels c1/c3/d draw on `figure1c_curves.csv` and the `clean`/`prev`/`conf` sweeps,
  which are in `sweep_all.csv` but not in `figure1_data.csv`.
- **`figure2a`** — all 10 contact amplitudes and all 10 "% of receivers < 10 µm" annotations
  reproduce from `figure2a_stratified_curves.csv` to 4 dp. Note the collision worth watching:
  the ALL-RECEIVERS panel prints "+0.27 sd" from **0.2660** (a curve amplitude), which rounds to the
  same as the unrelated authoritative naive biological amplitude **+0.2767** (section-clustered
  signed mean, `a7_summary.csv`). `CS_PHASE5.md:76` and `WRITING_PACK.md:123` already separate them;
  the figure does not say which it is.
- **`figure2d`** — 77 section × sender-definition combinations, slope −0.525, r² 0.984, exact
  against `poisson_fits.csv` row 0.
- **`figure_gs1`–`gs4`** — every A×B and B×B cell recomputes from the frozen `.txt` sets (0
  mismatches); all typed literals in `make_figure_genesets.py` check against a file; no caching.
  **CDKN2B**: present in both `genesets/A_SENDER_FINAL_strict.txt:8` (`Cdkn2b`) and
  `genesets/human/A_SENDER_FINAL_strict.txt:8` (`CDKN2B`); the figure now reads *"Human-only: DDB2
  GADD45G MDC1 PHLDA3 TNFRSF10B XPC"* plus *"plus CDKN2B, which both arms carry and the map has no
  row for"* — **defect fixed**. `figure_gs3`'s typed CoreScence literal is gone; its numbers
  recompute from `DeepScence/data/coreGS_v2.csv`.
- **`figure_phase8_callers` / `figure_phase8_d3`** — every plotted value exact against
  `results/phase3/*`; the `BAND = (0.93, 1.22)` literal is the min/max of
  `caller_agreement_depth_and_type_matched.csv`. Both PNGs have colliding panel titles
  (`make_figure_phase8_callers.py:109-111` vs `:123`; `make_figure_phase8_d3.py:50-53` vs `:95-96`),
  rendering text illegible in the committed images.
- **`figure4` panels a/b/c** — 136 of 138 rows exact against `headline.csv` / `interactions.csv.gz`.
