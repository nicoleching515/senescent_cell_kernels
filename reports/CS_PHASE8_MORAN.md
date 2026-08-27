# CS Phase 8 — Moran's I on the negative controls: the measured answer to §29 objection 9

**Run:** 2026-08-27. **Arm:** M1 (mouse) only, all 11 sections. `data/raw_h1/` was not
touched (§15 freeze). **Nothing was installed** — `libpysal` 4.14.1 and `esda` 2.8.2 were
already in the environment. **Nothing under `results/phase3/`, `figures/`, `genesets/`,
`data/processed/` or `SASP_Kernel_Master_Plan.md` was modified**; all new output is under
`results/moran/`. `python3 code/check_figures_guard.py` re-run at the end: **OK, all 52
committed figures match.**

---

## 0. The answer, in four sentences

1. **Voyager's result reproduces exactly on this data.** Per feature, the 40 negative
   control probes have a median Moran's I of **−0.00012**, the 609 codewords **−0.00004**,
   the 21 genomic controls **−0.00029**, against **+0.00999** for the 5,106 real genes.
   Read Voyager's way, this data says "no technical artifact spatial trend."
2. **That reading does not survive a sparsity control, and this is a new finding against
   the prior art.** Gene-expression features matched to the controls on total counts have a
   median Moran's I of **−0.00018 / −0.00012 / −0.00025** — indistinguishable from the
   controls, if anything slightly *more* negative. The genes that carry positive I are the
   abundant ones (median 5,885 counts per section vs **21** for a negative control probe, a
   280× gap). **The per-feature Moran's I that Voyager actually computes has no power at
   negative-control abundance**, so a near-zero value there is not evidence of anything.
3. **Aggregate the same features the way A7 does and the controls are not flat at all.**
   The per-cell sum over all 670 control features has Moran's I = **+0.0455
   [+0.0302, +0.0609]** (section-clustered, k = 6 NN) — above **86.8 %** of the
   56,166 gene × section Moran's I values and between a fifth and a half of the Tier B biological modules (+0.085 to +0.244).
   Voyager's conclusion is an artefact of feature abundance and of the level of aggregation,
   not a property of the assay.
4. **But — and this is the part that cuts against the project — the two statistics do not
   disagree the way `NOVELTY_ASSESSMENT.md` §2.1 point 3 assumes.** Across the twelve
   control-plus-module fields, |Moran's I| and |A7 naive amplitude| rank together at
   **Spearman ρ = +0.895 (p = 8.4 × 10⁻⁵)** raw, **+0.944** cell-type-centred — in both cases
   the **section-clustered mean per field** under knn6 weights (§4.1 gives the other
   aggregations, all of which agree in sign and significance). The sentence
   *"so the reader can see the two tests disagree"* is **not supported by this data and must
   not be written.**

**Does the differentiation survive?** Yes, but only in a restated, narrower and
*quantitative* form — and the restatement is stronger than the original claim, because it is
a measured power statement rather than a conceptual assertion. See §4.

---

## 1. What was run

| Piece | File |
|---|---|
| Producer, 11 sections × 7 weights variants × 2 centrings × 16 fields, + per-feature over all 13,590 h5 features | `code/run_moran_controls.py` → `results/moran/moran_fields.csv`, `moran_per_feature.csv.gz` |
| CP10K + log1p sensitivity (Voyager normalises `logcounts`), 3 sections | `code/run_moran_lognorm.py` → `moran_per_feature_lognorm.csv.gz` |
| The power calculation that settles the objection | `code/moran_kernel_power.py` → `moran_kernel_power.csv` |
| Summariser and the A7 join | `code/summarize_moran.py` → `moran_pooled.csv`, `moran_vs_a7.csv`, `moran_sensitivity.csv`, `moran_per_feature_summary.csv`, `moran_verdict.txt` |
| Figure (under `results/`, **not** `figures/`, so the guard is untouched) | `code/make_moran_figure.py` → `moran_summary.{png,pdf}` |

**Cells.** Exactly the analysis cell set the Phase-3 / A7 estimator uses — `sasp_phase3.Sec`,
i.e. QC-passing cells carrying a Bio cell-type label. 83,392–236,905 cells per section.

**Features.** Read from `cell_feature_matrix.h5` directly: 5,106 Gene Expression, 40 Negative
Control Probe, 609 Negative Control Codeword, 21 Genomic Control, 7,806 Unassigned Codeword,
8 Deprecated Codeword. The per-class sums were asserted equal, cell for cell, to
`cells.parquet`'s `control_probe_counts` / `control_codeword_counts` /
`genomic_control_counts` in every section (`assert` in `run_section`), so the aggregate
responses are byte-identical to A7's.

**Weights.** Primary **k = 6 nearest neighbours, row-standardised** (median 1-NN distance on
these sections is 9.5–10 µm, so k = 6 spans roughly the first cell shell, ~20 µm).
Sensitivity: k ∈ {4, 10, 20} and binary distance bands of **30, 50 and 100 µm** — the 100 µm
band is matched to the A7 fitting window, which answers the "you used a different scale"
objection directly. Islands (cells with no neighbour in the band) at most 143 / 61 / 29 cells
per section for the three bands; they are reported, not dropped.

**Inference.** Analytic *z* under the randomisation assumption (Cliff & Ord), the same
quantity `esda` reports as `z_rand`, plus **999 conditional permutations** for the 16
aggregate fields under the primary weights. **Validation:** `run_moran_controls.py --validate`
recomputes I, S0, S1, S2 and `z_rand` with `libpysal.weights.KNN` + `esda.moran.Moran` on a
4,000-cell subsample and agrees to **0.00e+00** on I and on `z_rand`. Output pasted in §7.

**Read the p-values with care.** At n ≈ 10⁵ cells, SE(I) ≈ 1.5 × 10⁻³, so almost anything
non-zero is "significant". Effect size is what carries the argument here, not p. The
exceptions are informative: `neg_control_probe`, `genomic_control` and `neg_probe_rate` each
have at least one section with p_rand ≥ 0.54 / 0.80 / 0.97 (permutation p agrees: 0.58 /
0.81 / 0.96), i.e. genuinely indistinguishable from zero *within a section*.

---

## 2. Per feature — the Voyager test, reproduced and then broken

`moran_per_feature_summary.csv`, figure panel (a). All 11 sections pooled; one row per
feature per section, primary weights, raw counts.

| feature class | n (feat × section) | median I | 5th–95th pct | max I | frac I > 0.05 | median total counts | median frac cells non-zero |
|---|---|---|---|---|---|---|---|
| Gene Expression | 56,166 | **+0.00999** | −0.0004, +0.0872 | +0.705 | **11.6 %** | 5,885 | 3.29 % |
| Genomic Control | 231 | −0.00029 | −0.0007, +0.0058 | +0.0137 | 0 % | 63 | 0.038 % |
| Negative Control Probe | 440 | −0.00012 | −0.0004, −0.00003 | +0.0201 | 0 % | 21 | 0.013 % |
| Negative Control Codeword | 6,699 | −0.00004 | −0.0002, −0.000006 | +0.0833 | 0.06 % | 7 | 0.004 % |
| Unassigned Codeword | 85,866 | −0.00004 | −0.0002, −0.000006 | +0.167 | 0.04 % | 6 | 0.004 % |
| Deprecated Codeword | 88 | −0.00002 | −0.0003, +0.0942 | +0.161 | 12.5 % | 26 | 0.016 % |

(251 codeword and 3,920 unassigned-codeword feature × section cells are constant-zero and
return NaN; they are excluded from the medians rather than counted as I = 0.)

**This is Voyager's figure.** "The negative controls are tightly clustered around 0, while
the real genes have positive Moran's I."

**And here is why it does not license Voyager's conclusion.** For each control class,
compare against genes whose total counts fall inside the class's own 10th–90th percentile
band (`moran_verdict.txt`):

| control class | n | median counts | median I | matched genes | median counts | median I |
|---|---|---|---|---|---|---|
| Negative Control Probe | 440 | 21 | −0.00012 | 4,495 | 32 | **−0.00018** |
| Negative Control Codeword | 6,699 | 7 | −0.00004 | 1,179 | 18 | **−0.00012** |
| Genomic Control | 231 | 63 | −0.00029 | 8,735 | 62 | **−0.00025** |
| Unassigned Codeword | 85,866 | 6 | −0.00004 | 763 | 15 | **−0.00011** |

A gene as sparse as a negative control probe looks exactly as "flat" as the probe does. The
contrast Voyager reads as *controls vs genes* is, on this data, **abundance vs abundance**.
This is a reportable methodological finding about a widely used diagnostic, it is friendly
rather than adversarial to the Pachter lab (the vignette is a tutorial, not a claim), and it
belongs in the paper regardless of what happens to the rest of the argument.

**Normalisation makes no difference.** Voyager computes Moran's I on the SFE `logcounts`
assay. Recomputing every feature on CP10K + log1p over three sections (below-floor, in-band,
over-ceiling) gives median I of −0.00011 (probes), −0.00004 (codewords), −0.00028 (genomic),
+0.00718 (genes) — the same picture. `moran_per_feature_lognorm.csv.gz`.

---

## 3. Aggregated — the A7 responses, side by side with the A7 kernel

`moran_pooled.csv`, figure panel (b). Moran's I pooled as a **section-clustered mean with a
t-CI on the 11 section values**, matching how `summarize_a7.py` pools amplitudes. A7 column
is `beta_base / sd_y`, the naive design, same responses, same cells, from
`results/phase3/a7_control_probe_fits.csv` and `main_fits.csv`.

| field | Moran's I (raw) | Moran's I (cell-type-centred) | A7 naive β/sd_y |
|---|---|---|---|
| `neg_control_probe` | **+0.0058 [+0.0033, +0.0083]** | +0.0054 [+0.0027, +0.0080] | **−0.0225 [−0.0527, +0.0078]** |
| `neg_control_codeword` | **+0.0421 [+0.0281, +0.0561]** | +0.0391 [+0.0246, +0.0536] | **−0.0604 [−0.1085, −0.0123]** |
| `genomic_control` | +0.0042 [+0.0022, +0.0062] | +0.0040 [+0.0019, +0.0062] | −0.0307 [−0.0558, −0.0056] |
| `all_controls` | **+0.0455 [+0.0302, +0.0609]** | +0.0424 [+0.0265, +0.0583] | **−0.0744 [−0.1306, −0.0182]** |
| `neg_probe_rate` | **+0.0047 [−0.0010, +0.0105]** | +0.0043 [−0.0013, +0.0098] | **+0.0113 [−0.0085, +0.0310]** |
| `downstream_arrest` | +0.0852 [+0.0645, +0.1058] | +0.0587 | +0.0989 [+0.0279, +0.1700] |
| `secondary_senescence` | +0.0872 [+0.0648, +0.1096] | +0.0721 | +0.2270 [+0.1263, +0.3277] |
| `il6_jak_stat3` | +0.1205 [+0.0769, +0.1640] | +0.1001 | +0.5159 [+0.3915, +0.6403] |
| `oxidative_stress` | +0.1502 [+0.1085, +0.1918] | +0.0810 | −0.2160 [−0.2948, −0.1372] |
| `interferon_response` | +0.1719 [+0.1343, +0.2095] | +0.1382 | +0.3944 [+0.2693, +0.5196] |
| `tnfa_nfkb_proximal` | +0.1903 [+0.1190, +0.2616] | +0.1695 | +0.5728 [+0.4231, +0.7226] |
| `emt_ecm` | +0.2437 [+0.1910, +0.2964] | +0.1431 | +0.3437 [+0.2558, +0.4315] |
| `transcript_counts` | +0.1375 [+0.1100, +0.1651] | +0.1320 | — |
| `cell_area` | +0.0855 [+0.0640, +0.1071] | +0.0696 | — |
| `unassigned_codeword` | +0.2098 [+0.1555, +0.2642] | +0.2039 | — |
| `density_50um` | +0.9483 [+0.9362, +0.9605] | +0.8664 | — |

### 3.1 The A7 split reproduces, as an ordering, not as a dichotomy

The brief asked whether Moran's I reproduces A7's split — probes flat, codewords and the
pooled set not. **It does, and more sharply.** The probes carry **7×** less spatial
autocorrelation than the codewords (+0.0058 vs +0.0421) where A7's amplitudes differ by only
**2.7×** (0.0225 vs 0.0604). The pooled `all_controls` field tracks the codewords, because
609 codewords dominate 40 probes in the sum. `neg_probe_rate` — probes divided by transcript
counts — is the one response whose Moran's I CI **includes zero**, exactly as it is the one
response A7 calls flat naively. That is four independent corroborations of A7's response-level
structure from a statistic with no kernel, no λ, no sender call and no nuisance design.

**Report against interest:** the split reproduces as an *ordering of magnitudes*, not as a
flat/not-flat verdict. Moran's I does **not** call the 40 probes flat — pooled over 11
sections its CI excludes zero (+0.0058 [+0.0033, +0.0083]). A7 calls them flat only because
A7's CI is wide (±0.030 pooled, ±0.126 per fit). Anyone writing "both tests agree the probes
are flat" would be overstating it. The correct sentence is *both tests put the probes an
order of magnitude below the codewords.*

Second discrepancy, also against interest: `genomic_control` is A7's third-largest amplitude
(−0.0307, CI excludes zero) but Moran's **smallest** I (+0.0042). The two statistics disagree
on that one response's rank.

### 3.2 Weights sensitivity — the conclusion does not depend on the graph

`moran_sensitivity.csv`, mean over 11 sections:

| field | knn4 | **knn6** | knn10 | knn20 | band30 (16.1 nb) | band50 (43.4 nb) | band100 (166 nb) |
|---|---|---|---|---|---|---|---|
| `all_controls` | 0.0452 | **0.0455** | 0.0447 | 0.0423 | 0.0463 | 0.0421 | 0.0353 |
| `neg_control_codeword` | 0.0415 | **0.0421** | 0.0415 | 0.0395 | 0.0428 | 0.0393 | 0.0332 |
| `neg_control_probe` | 0.0064 | **0.0058** | 0.0057 | 0.0053 | 0.0059 | 0.0056 | 0.0046 |
| `genomic_control` | 0.0045 | **0.0042** | 0.0037 | 0.0033 | 0.0037 | 0.0033 | 0.0026 |
| `neg_probe_rate` | 0.0050 | **0.0047** | 0.0037 | 0.0029 | 0.0035 | 0.0024 | 0.0021 |
| `emt_ecm` (positive control) | 0.2656 | **0.2437** | 0.2184 | 0.1881 | 0.1928 | 0.1542 | 0.1055 |
| `transcript_counts` | 0.1177 | **0.1375** | 0.1464 | 0.1363 | 0.1418 | 0.1267 | 0.1002 |

Everything decays smoothly with neighbourhood size, as it must; **no ordering changes**, and
the 100 µm band — the A7 window itself — still gives `all_controls` = 0.0353, above the
median gene. The result is a property of the data, not of the weights choice.

---

## 4. The comparison that answers objection 9 — and what it does to the claim

### 4.1 The two statistics agree at the field level. Say so.

Across the twelve control + module fields, |Moran's I| against |A7 naive amplitude|:

* raw: **Spearman ρ = +0.895, p = 8.4 × 10⁻⁵**
* cell-type-centred: **ρ = +0.944, p = 3.9 × 10⁻⁶**

**Both are the section-clustered *mean* per field, knn6 row-standardised weights, over the
12 control + module fields** — `code/summarize_moran.py:183-184` → `moran_verdict.txt`.
**The aggregation must be stated wherever this ρ appears, because ρ moves with it**
(added 2026-08-27, record reconciliation):

| aggregation | ρ | p | emitted by |
|---|---|---|---|
| clustered **mean** per field, knn6 raw, 12 fields | **+0.8951** | 8.37e-05 | `moran_verdict.txt` — **frozen; the one to quote** |
| clustered **mean** per field, knn6 cell-type-centred, 12 fields | **+0.9441** | 3.93e-06 | same — **frozen** |
| **median** per field, knn6 raw, 12 fields | +0.9231 | 1.86e-05 | no file; re-derivable from `moran_vs_a7.csv` |
| **per-row**, no aggregation, 12 fields × 11 sections = 132 pairs | +0.7104 | 1.43e-21 | no file; same |

**All four are positive and significant, so the falsification is not fragile: at every
defensible aggregation, Moran's I and the A7 kernel agree rather than disagree.** Only the
digit is aggregation-dependent. Do not let a reader infer that the finding turns on the
choice. (The separate within-control-family value, ρ = +0.155, p = 0.259 over the 5 control
responses × 11 sections, is a **different subset**, not a fifth aggregation of this one.)

Panel (b) of the figure shows it: modules top-right, controls bottom-left, monotone between.
**The plan's implied defence — "the two tests disagree, look" — is false on this data.**
`NOVELTY_ASSESSMENT.md` §2.1 point 3 tells the project to write exactly that sentence. It
must not be written. If it were, a reviewer with `results/moran/` in hand would find the
opposite in one plot.

### 4.2 …but they agree for a reason that does not rescue Moran's I

The agreement is *not* Moran's I detecting the kernel. That can be computed, and it was
(`moran_kernel_power.py`, `moran_kernel_power.csv`, 22 section × sender-call cells).

Under the A7 model y = Xγ + β·k, k_i = exp(−d_i/λ) with d the distance to the nearest sender,
if the rest of y is spatially white then the gradient's contribution to Moran's I is
β_z² · Var(k) · I(k). Measured on the real receiver sets, with the real λ̂ per fit:

| quantity | median over 22 section × call | range |
|---|---|---|
| λ̂ (naive, `all_controls`) | 35.6 µm | 7.0 – 50.0 |
| Moran's I of the kernel covariate itself, I(k) | **+0.797** | +0.412 – +0.852 |
| Var(k) | 0.0309 | 0.0085 – 0.0408 |
| A7 amplitude β_z used | −0.070 | −0.279 – +0.164 |
| **ΔI contributed by the whole A7 gradient** | **2.20 × 10⁻⁴** | 1.4 × 10⁻⁶ – 2.6 × 10⁻³ |
| **ΔI as a fraction of the observed control I** | **0.83 %** | 0.002 % – 6.1 % |
| SE(I) on the same cell set | 1.50 × 10⁻³ | 1.2 – 2.3 × 10⁻³ |
| **β_z that Moran's I could just detect (ΔI = 2 SE)** | **0.362 SD** | 0.308 – 1.070 |

Three consequences, and they are the substance of the answer to objection 9:

1. **The gradient A7 measures accounts for under 1 % of the Moran's I on the same field.**
   The +0.0455 is not the confound A7 found; it is bulk short-range clustering of detection
   efficiency (`density_50um` I = +0.948, `transcript_counts` I = +0.138). The two statistics
   co-vary because one confound drives both through different channels, not because Moran's I
   can see the projection.
2. **Moran's I is blind to a distance-to-sender kernel at the amplitudes that matter here.**
   The smallest amplitude it could resolve is **0.362 SD**. The A7 control gradient is
   **0.074 SD** (5× too small). The paper's own **naive biological amplitude is 0.277 SD**
   and its **conditioned amplitude is 0.031 SD** — both the section-clustered signed mean of
   β/sd(y) over the 1,155 biological-module fits (`results/phase3/a7_summary.csv`, row
   `BIOLOGICAL MODULES (reference)`, `design = base` / `n6n5`), and both **frozen post-C6**.
   *(Corrected 2026-08-27, record reconciliation: this read 0.291 / 0.036, the pre-C6 vintage
   of the same estimator. The companion estimator on the same fits, median |β|/sd, gives
   0.312 naive / 0.0795 conditioned — name whichever you use. The power argument is unaffected:
   0.362 exceeds all of them.)* Moran's I could not detect the project's
   *headline effect* either, let alone the confound. A statistic that cannot see the estimand
   is not a test of the estimand.
3. **Therefore the "different question" claim is true, but its correct justification is
   power, not orthogonality.** It is not that a near-zero I *precludes nothing in principle*;
   it is that on this data a global autocorrelation statistic would need the effect to be an
   order of magnitude larger before it moved at all.

### 4.3 The discordance count, reported honestly

Among the 55 section × control-response cells, **14** have |Moran's I| below the median gene
(0.00999) while ≥ 25 % of that cell's A7 fits have naive CIs excluding zero — including 8 of
the 11 sections for `neg_control_probe`. That looks like the disagreement the plan wanted.
**It should not be quoted as such.** A7's own measured CI-exclusion rate on a response with
no biology is **9–16 % against a 5 % nominal** (`a7_verdict.txt`), so a 25–50 % exclusion rate
in a single section is partly the estimator's own inflated type-I error, not proof of a
gradient Moran's I missed. The clean statement is the power calculation in §4.2, which does
not depend on any single fit being right. Within the control family, section by section,
Moran's I and the A7 amplitude are **uncorrelated** (ρ = +0.155, p = 0.26, 55 pairs) — which
is what §4.2 predicts and is the honest form of "these are different questions."

---

## 5. Drop-in replacement for §29 objection 9

Offered as text, not applied — `SASP_Kernel_Master_Plan.md` is outside this task's edit scope.
It removes the ⚠ STATUS marker and the unsupported "the two tests disagree" line.

> **9. "Voyager already computes Moran's I on Xenium negative control probes and concludes
> there is no technical spatial trend. What is new?"** We cite Voyager and Ren et al. (2025)
> approvingly, and we ran their test. Per feature it reproduces: median Moran's I is −0.00012
> for the 40 negative control probes and −0.00004 for the 609 codewords, against +0.00999 for
> the 5,106 genes (11 sections, k = 6 NN row-standardised; `results/moran/`). **It does not
> survive a sparsity control:** genes matched to the controls on total counts have a median I
> of −0.00018, i.e. the *controls vs genes* contrast is an *abundance* contrast — a negative
> control probe carries 21 counts per section against 5,885 for the median gene. Aggregated
> the way our estimand's response is defined — the per-cell sum over all 670 control features —
> the same data gives I = **+0.0455 [+0.0302, +0.0609]**, above 87 % of the per-gene values, so this
> assay does have a control-side spatial trend. Our construction is not the Moran's I but a
> **negative control outcome for the estimand itself** (Lipsitch et al. 2010): the identical
> distance-to-sender kernel, refit with control counts as the response, giving −0.074 SD
> naively, +0.007 SD under the full nuisance design and −0.061 SD under a matched-decoy
> contrast. The two statistics are not redundant, and the reason is measurable: the entire
> fitted gradient contributes **0.83 %** of the observed Moran's I, and the smallest kernel
> amplitude Moran's I could resolve on these sections is **0.36 SD**, five times our control
> gradient and larger than our own naive biological amplitude of 0.29 SD. **The sentence
> "nobody in this literature reports it" must be struck wherever it appears.**

Two further edits this run forces:

* `SASP_Kernel_Master_Plan.md` §30 5.5 says "*Report a Moran's I of our own controls beside
  the amplitude (**not yet computed — see objection 9**)*". The parenthesis can be replaced
  with `+0.0455 [+0.0302, +0.0609] pooled, +0.0058 on the 40 probes; results/moran/`.
* `reports/NOVELTY_ASSESSMENT.md` §2.1 point 3 and §4 O1 both instruct the project to say
  "the two tests disagree". **That instruction is now falsified** (ρ = +0.90 across fields)
  and should be replaced with the power statement.

---

## 6. Limits of this run

* **Mouse only.** H1 is untouched by design. The sparsity argument will transfer (the human
  panel has the same control-feature design), but the numbers will not.
* **One graph family.** k-NN and distance bands, row-standardised. Not tested: Delaunay,
  Gabriel, kernel-decay weights. Voyager's Xenium vignette uses its own graph; the seven
  variants here span 4–166 mean neighbours without changing an ordering, so this is a
  robustness gap, not a live risk.
* **Global Moran's I only.** No local Moran / LISA, no Geary's C, no correlogram beyond what
  `results/phase3/correlogram.csv` already holds (module ρ at 10 µm = 0.06–0.45, ℓ = 53–114 µm;
  consistent with the module I values here). Ren et al. (2025) compare platforms; that
  comparison is out of scope with one platform.
* **The power calculation assumes the non-kernel part of y is spatially white.** It is not —
  that is the whole point of the +0.0455 — but the assumption enters only in *attributing* the
  observed I, and it makes the "ΔI is 0.83 % of I" statement **conservative** in the direction
  that matters: any spatial structure in the residual raises the denominator, so the true
  fraction is at most this. The β_min figure does not depend on the assumption at all.
* **`neg_probe_rate` has a ratio denominator that is itself an N5 column.** Its near-zero
  Moran's I inherits that caveat exactly as its A7 amplitude does.

---

## 7. Verification log

```
$ python3 code/run_moran_controls.py --validate
libpysal S0=4000.0000 mine S0=4000.0000  S1 1213.2222/1213.2222  S2 16390.0000/16390.0000
  control_probe_counts   esda I=-0.001837 mine -0.001837 (d=0.00e+00) | z_rand esda -0.1949 mine -0.1949 (d=0.00e+00) | p_rand esda 0.8455 mine 0.8455 | esda p_sim 0.264
  transcript_counts      esda I=+0.117326 mine +0.117326 (d=0.00e+00) | z_rand esda +13.5416 mine +13.5416 (d=0.00e+00) | p_rand esda 8.878e-42 mine 8.878e-42 | esda p_sim 0.001
  emt_ecm                esda I=+0.030744 mine +0.030744 (d=0.00e+00) | z_rand esda +3.5646 mine +3.5646 (d=0.00e+00) | p_rand esda 0.0003645 mine 0.0003645 | esda p_sim 0.001
VALIDATION PASS

$ python3 code/run_moran_controls.py --sections all --n-jobs 3 --perms 999
[Parallel(n_jobs=3)]: Done  11 out of  11 | elapsed: 12.3min finished
(2464, 20) -> /workspace/results/moran/moran_fields.csv
(149490, 10) -> /workspace/results/moran/moran_per_feature.csv.gz

$ python3 code/check_figures_guard.py
OK: all 52 committed figures match (PDF date stamps ignored)
```

Per-class h5 sums were asserted equal to the `cells.parquet` tallies in all 11 sections
(hard `assert` in `run_moran_controls.run_section`; the run would have aborted otherwise).
Peak cgroup memory during the run: 13.9 GB of the 57.7 GB ceiling, at `--n-jobs 3`.
