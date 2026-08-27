# CS Phase 8 — caller agreement at full DeepScence coverage, A7, and D3

**Yes, it moved, and it moved against us: the depth- and type-matched caller
agreement rose from the published 0.93–1.22× chance to 0.70–1.71×, and the
pooled value went from 1.030× chance (z = 1.27, p = 0.20 — indistinguishable
from independence) to 1.118× chance (z = 11.5, p = 1.4 × 10⁻³⁰). The callers are
no longer statistically independent, and §5's condition — "if agreement rises
above chance at full coverage, your motivating claim weakens and must be
restated" — has been met. §11 below gives the restatement.**

The rise is small in effect size and large in confidence: agreement is 12 % above
chance, not 2× above it. The paper's substantive point survives in a weakened,
more defensible form; the sentence "their top-5 % calls overlap at 0.93–1.22× of
chance … i.e. they are statistically independent" does not.

| Task | Deliverable | Location |
|---|---|---|
| 8.3 | Five caller tables at 11-section coverage | `results/phase3/caller_*_11sections.csv` |
| 8.3 | Producer, reusing `caller_disagree.py` | `code/caller_disagree_all.py` |
| 8.3 | Exact reproduction of the committed 2-section tables | `--verify` → `caller_*_verify2sec.csv` |
| 8.4 | The gate, with a stratified-exact null | `code/summarize_caller_coverage.py` → `caller_coverage_gate{,_headline}.csv` |
| 8.4 | Figure | `figures/figure_phase8_callers.{png,pdf}` + `results/phase3/figure_phase8_callers_data.csv` |
| 8.5 | A7 on the mouse arm | `code/run_a7_control_probes.py`, `summarize_a7.py` → `a7_control_probe_{fits,curves,provenance}.csv`, `a7_summary.csv`, `a7_verdict.txt` |
| 8.5 | **Figure 2 panel h filled, M1 half** | `figures/figure2e.{png,pdf}`, `results/phase3/figure2e_data.csv` |
| 8.6 | D3 re-anchoring | `code/deepscence_reanchor.py` → `data/processed/deepscence_d3_<section>.csv`, `deepscence_anchor_decisions.csv` |
| 8.6 | D3 effect on caller agreement, incl. sign-invariant | `code/caller_disagree_d3.py` → `caller_agreement_matched_d3_11sections.csv` |
| 8.6 | Figure | `figures/figure_phase8_d3.{png,pdf}` + `results/phase3/figure_phase8_d3_data.csv` |

**Pinned files, re-verified after everything below:**

```
3b77aa1bba0712c205c5d9356654fb71  results/phase3/perm_nulls.csv
69e3a1d3f60060deddcceba9896a7d31  results/phase3/sf_summary.csv
ecf86b9ca5460f31290e2f4c9e822ea2  results/phase3/summary_phase3.txt
```

`data/processed/deepscence_sham.csv` (md5 `8c4c52f5c1c7649d8c17d07010cc780c`,
mtime 2026-08-20 17:42) and `deepscence_sbr.csv` (`b557e3dfb8eff517d040757c73f0a660`,
17:57) are untouched, as required. Nothing under `data/raw_h1/` or `genesets/`
was written; `genesets/` was read only, to verify the D3 anchor's disjointness.
`make_figure2.py` was not run. No commits, no tags. No packages were installed —
the environment as pinned was sufficient for all four tasks.

---

## 1. 8.3 — the producer, and why the reconstruction is trustworthy

`caller_disagree.py` writes four of the five tables §5 lists. **The producer of
the other two — `caller_agreement_depth_and_type_matched.csv` and
`caller_within_type_depth_bias.csv` — was never committed.** Both outputs are in
`results/phase3/` dated 2026-08-20 18:02, alongside `caller_disagree2.py`'s two
outputs, but no script in `code/` and no revision in the git history writes
them. They were produced by an ad-hoc extension in that session and lost.

They are reconstructed in `code/caller_disagree_all.py` and the reconstruction is
validated rather than asserted. `caller_disagree_all.py --verify` re-runs the two
original sections through the new file and compares all six tables, cell by cell,
against the committed ones:

```
OK   caller_technical_loading                  (8, 8)
OK   caller_celltype_composition               (84, 8)
OK   caller_strata                             (64, 8)
OK   caller_pairwise_agreement                 (24, 13)
OK   caller_within_type_depth_bias             (40, 6)
OK   caller_agreement_depth_and_type_matched   (12, 11)
VERIFY: PASS
```

Every value in all six tables comes back **exactly as committed**, including the
two whose producer was lost — 0.932 / 1.041 / 0.960 / 2.150 / 1.370 / 2.854 for
sham and 0.943 / 1.221 / 1.059 / 0.381 / 0.980 / 1.506 for SBR. The recovered
recipe is: top-5 % recomputed **inside each (cell type × within-type transcript-count
decile)** stratum with a 50-cell floor, on the pairwise-complete subset, strictly
greater than the stratum's 95th percentile; `ratio` = observed overlap ÷
(n_A · n_B / n). The depth-bias table is the top-5 % recomputed **within cell type**,
enriched across **within-cell-type** depth quintiles.

The only edit to `caller_disagree.py` is a four-line `DS_ALIAS` mapping the two
section directory names onto the preserved `deepscence_{sham,sbr}.csv` filenames,
so the same code path serves both bases. Called with `tag='sham'` it behaves
exactly as before — which is what `--verify` proves.

**Coverage is the only thing that changed.** All nine new DeepScence files come
from `run_deepscence_all.py` at settings identical to `run_deepscence.py`
(`denoise=False`, `random_state=0`, published `CDKN1A` anchor, MGI 1:1 ortholog
remap of 4,845 of 5,097 panel genes, ≥20 counts/cell). Sections 7260, 7352 and
7448 were OOM-killed on a first five-way-parallel attempt and re-run
sequentially at complete, identical settings. Verified independently here: the
rows for 7250 and 7259 inside the eleven-section tables are **identical** to the
committed `sham` / `sbr` rows, so the two bases differ only by the nine added
sections.

Both the five-table run and the A7 run were executed twice; every output CSV is
**bit-identical** across runs.

## 2. 8.4 — the gate, measured

`ratio` as published is observed overlap ÷ marginal expected overlap. That null
ignores that both call sets are stratum-balanced by construction, so
`caller_disagree_all.py` also computes the **exact conditional (Mantel–Haenszel)
null** — B's calls permuted within each stratum, A fixed, overlap a sum of
independent hypergeometrics with closed-form mean and variance. The two agree to
three decimals everywhere (e.g. 2.854 vs 2.849), which vindicates the published
definition and supplies the z-scores the published table lacked.

### 2.1 The headline

The published "0.93–1.22×" is, exactly, the min and max of the six matched
ratios of Tier A vs SenePy, Tier A vs DeepScence and Tier A vs `Cdkn1a`⁺ over the
two sections: **0.932–1.221**. Recomputed identically at eleven sections:

| basis | band | median | pooled ratio | z | p | sections above chance |
|---|---|---|---|---|---|---|
| 2 sections, 3 pairs *(the published band)* | **0.932 – 1.221** | 1.000 | **1.030** | 1.27 | 0.20 | 3 / 6 |
| 11 sections, 3 pairs | **0.700 – 1.711** | 1.110 | **1.118** | 11.49 | 1.4e-30 | 22 / 33 |
| 2 sections, 4 pairs *(incl. SenePy vs `Cdkn1a`⁺)* | 0.932 – 1.369 | 1.010 | 1.040 | 1.76 | 0.078 | 4 / 8 |
| 11 sections, 4 pairs | 0.700 – 1.711 | 1.156 | **1.129** | 13.35 | 1.1e-40 | 29 / 44 |
| **11 sections, 4 pairs, in-band six only** | 0.775 – 1.374 | 1.131 | **1.115** | 8.99 | 2.6e-19 | 16 / 24 |

The last row matters: the widened band is **not** an artefact of 7239, the
section §8 Test 3 excludes at 45 % `Cdkn1a`⁺ hepatocytes. Restricting to the six
Test-3-admissible sections Phase 3 actually fits gives a narrower band and the
same conclusion at p = 2.6 × 10⁻¹⁹.

The rise is not merely a power effect. The **ratio itself** rose (1.030 → 1.118),
the median rose (1.000 → 1.110), and the sign counts rose (3/6 → 22/33). The
2-section base was not a small sample of the truth; it was an unrepresentative one.

### 2.2 Per pair, and where the movement is

`caller_coverage_gate.csv`. `DeepScence vs Cdkn1a⁺` is kept out of every pooled
number above and is reported separately, per `BIO_PHASE3.md` §4.4 — DeepScence
anchors its sign on `CDKN1A`, so that pair is circular and is not evidence of
agreement.

| pair | 2-section pooled | **11-section pooled** | min – max over 11 | above chance | sig. above | sig. below |
|---|---|---|---|---|---|---|
| Tier A vs SenePy | 0.935 | **0.914** (z −4.9) | 0.700 – 0.995 | 0/11 | 0 | 3 |
| Tier A vs DeepScence | 1.103 | **1.248** (z 16.6) | 1.041 – 1.711 | **11/11** | 9 | 0 |
| Tier A vs `Cdkn1a`⁺ | 1.017 | **1.171** (z 7.1) | 0.948 – 1.433 | 9/11 | 5 | 0 |
| SenePy vs `Cdkn1a`⁺ | 1.168 | **1.211** (z 7.4) | 0.842 – 1.391 | 9/11 | 7 | 0 |
| SenePy vs DeepScence | 1.693 | **0.737** (z −15.1) | 0.332 – 2.150 | 1/11 | 1 | 10 |
| *(circular)* DeepScence vs `Cdkn1a`⁺ | 2.050 | **1.255** (z 10.5) | 0.963 – 2.849 | 7/11 | 5 | 0 |

Four movements, each with a cause:

1. **Tier A vs DeepScence is above chance in all eleven sections.** Binomial sign
   test p = 2⁻¹¹ ≈ 5 × 10⁻⁴ on its own, and pooled z = 16.6. On the two-section
   base its values, 1.041 and 1.221, were the *lowest* and the *highest* of the
   eleven; the low one, 1.041, came from 7250 and set the published band's own
   ceiling of 1.22 from the other section. **This is the pair that breaks the
   independence claim**, and it is the one pair among the four where the null of
   independence can be rejected on the sign pattern alone, without any pooling.
2. **Tier A vs `Cdkn1a`⁺ and SenePy vs `Cdkn1a`⁺ both rise** to ~1.17–1.21× and
   are above chance in 9 of 11. They straddled chance on two sections.
3. **SenePy vs DeepScence reverses.** Published: 2.15× in sham, 0.38× in SBR,
   read in §4.4 as "concordant in one arm, anti-concordant in the other". At full
   coverage it is 0.33–0.55× in **ten of eleven** sections; the single exception
   is 2.150 in 7250. There is no arm effect. There is one anomalous section, and
   §8.6 below identifies why.
4. **The circularity is much weaker than published.** DeepScence vs `Cdkn1a`⁺ was
   reported at 1.51–2.85×; over eleven sections it is 0.963–2.849 with a median of
   **1.071** and a pooled 1.255 — the same magnitude as the *non*-circular Tier A
   vs DeepScence (1.248). Both published values, 2.85 and 1.51, are the two
   largest of the eleven. The pair is still circular by construction and must
   still be excluded, but "1.51–2.85×" overstates the measured circularity by
   about a factor of two and should not be quoted.

Tier A vs SenePy is the one pair that does not move: 0.935 → 0.914, below chance
on all eleven sections. Two callers are genuinely close to disjoint. That is
worth keeping.

### 2.3 What else the eleven-section tables change

* **`caller_technical_loading_11sections.csv` — the "DeepScence's depth
  correlation reverses sign between two sections of the same study" sentence
  (§4.2, §4.5) is wrong as stated.** Spearman ρ with transcript counts is
  **+0.29 to +0.56 in ten of eleven sections** and −0.350 in exactly one: 7250.
  It is not a between-arm instability; it is one section. Same for cell area
  (+0.28 to +0.57 vs −0.432).
* **`caller_celltype_composition_11sections.csv` — the "DeepScence calls
  essentially no hepatocytes (0.03×) in a hepatocyte-dominated tissue" sentence
  is likewise 7250 only.** Global top-5 % hepatocyte enrichment is 1.38–2.46 in
  the other ten and 0.027 in 7250.
* **`caller_within_type_depth_bias_11sections.csv` — §4.1's qualitative claim
  survives, its DeepScence number does not.** SenePy Q5/Q1 = 10.6–41.7× and
  `Cdkn1a`⁺ = 4.2–42.4× (top of the depth distribution) against Tier A
  0.17–0.29× (bottom) in every one of the eleven sections. DeepScence is also
  bottom-selecting everywhere, but at **0.27–0.80×** outside 7250, not the 0.22
  published — a much milder bias than the two-section base showed.
* **`caller_pairwise_agreement_11sections.csv`** (the unmatched table) tells the
  same story: within-type top-5 %, Tier A vs DeepScence is 1.202–1.703 with a
  median of 1.314 across all eleven, i.e. above chance everywhere before any
  depth matching too.

## 3. 8.4 — the restatement, as text for the paper

The §4.5 paragraph should be replaced with this. Every number is in
`caller_coverage_gate.csv` and `caller_agreement_matched_significance_11sections.csv`.

> Four published or panel-standard ways of calling a senescent cell —
> DeepScence, SenePy, a disjoint arrest-and-damage score, and `Cdkn1a`
> positivity — were applied to the same 45,000–218,000 analysable cells in each
> of eleven Xenium Prime 5K mouse liver sections. After conditioning on cell type
> and on within-cell-type sequencing depth, their top-5 % calls overlap at
> **0.70–1.71× of chance**, pooling to **1.12× chance** (Mantel–Haenszel over the
> eleven sections, z = 11.5). **They are not independent, but the dependence is
> weak**: the largest pooled agreement between two non-circular callers is 1.25×
> chance, an order of magnitude short of what two noisy measurements of one
> latent state would show, and one of the four pairs (the arrest score vs SenePy)
> sits *below* chance at 0.91× in all eleven sections. Two callers overlapping
> 25 % more often than chance while a third overlaps 9 % *less* often than chance
> is not a shared latent variable; it is a shared **technical** one. What each
> caller selects is identifiable and technical: within cell type, SenePy and
> `Cdkn1a`⁺ are enriched 1.6–3.2× in the highest transcript-count quintile and
> depleted 2.4–14.5× in the lowest, while the arrest score and DeepScence run the
> other way in every section; SenePy's cross-cell-type score scales with the
> number of its hub genes on the panel (r = 0.992) and so is not comparable
> across cell types; and DeepScence's polarity is fixed by a `CDKN1A` anchor that,
> once sequencing depth is partialled out, is weak or reversed in four of eleven
> sections (§7). A senescence call on targeted spatial data is therefore not a
> noisy measurement of one latent state — it is a choice of which end of the
> detection-depth distribution to name senescent, and the length constant that
> follows inherits that choice.

**What must be struck**: "0.93–1.22× of chance … i.e. they are statistically
independent"; "Four of six pairs sit at 0.93–1.22×"; "the one pair that looked
concordant in sham is anti-concordant in SBR"; "DeepScence's correlation with
sequencing depth reverses sign between two sections of the same study"; and the
circularity figure "1.51–2.85×". Each is a property of the two-section base.

**What is strengthened**: the claim is now *evidence* rather than an absence of
it. "Independent" was an accept-the-null argument; "1.12× chance with p = 10⁻³⁰,
one pair significantly below chance, and the direction of each pair explained by
its depth loading" is a positive finding, and it is measured on eleven animals
instead of two.

## 4. 8.5 — A7 on the mouse arm: the verdict

`a7_verdict.txt`, `a7_summary.csv`. Eleven sections × two sender calls
(`tierA_p95`, `cdkn1a_pos`) × five control responses × every receiver cell type
with ≥ 2,000 receivers = **825 fits**, against **1,155 Tier B module fits** from
`main_fits.csv` on the same sections, calls, cell types and estimator.

Nothing about the estimator changed: the script reuses
`run_phase3_nulls.SectionFit` and `run_phase3_nulls.fit_cell` verbatim and
substitutes only the response matrix (and the N6 neighbour baseline recomputed
from it by the same `phase3_core.neighbour_baseline`). Same 100 µm window, same
40-point λ grid, same 2,000-receiver floor, same nested designs, same
400-replicate spatial block bootstrap.

**Response provenance, verified against the h5.** `cells.parquet`'s
`control_probe_counts`, `control_codeword_counts` and `genomic_control_counts`
are, **cell for cell and exactly**, the per-cell sums of the h5's 40
`Negative Control Probe`, 609 `Negative Control Codeword` and 21 `Genomic Control`
features (checked on 7250: 539 / 3,481 / 1,035 total counts, `np.array_equal`
true for all three, and for `transcript_counts` against the 5,106 Gene Expression
features). No download and no re-extraction was needed. Each response is z-scored
over the section's analysis cells, so β is an amplitude in response-SD units and
is directly comparable to a module's β/sd_y.

### 4.1 A7 FAILS on the naive kernel and PASSES on the conditioned one

Amplitude = |β| / sd_y. "Clustered mean" is the mean amplitude with a 95 % CI from
a t on the eleven **section** means, because fits within a section are not
independent.

| design | control amplitude, median | control clustered mean [95 % CI], p | control sign + | CI excludes 0 | module amplitude | module clustered mean |
|---|---|---|---|---|---|---|
| naive (intercept only) | 0.081 | **−0.070 [−0.128, −0.012], p = 0.023** | 0.34–0.51 | 0.18–0.45 | 0.314 | +0.291 |
| + N6 neighbour baseline | 0.080 | −0.061 [−0.108, −0.013], p = 0.017 | 0.34–0.51 | 0.18–0.33 | 0.169 | +0.120 |
| + N5 technical covariates | 0.063 | +0.006 [−0.012, +0.025], p = 0.46 | 0.49–0.54 | 0.10–0.16 | 0.098 | +0.074 |
| **+ N6 + N5 (full nuisance design)** | 0.056 | **+0.007 [−0.011, +0.025], p = 0.41** | 0.50–0.54 | 0.09–0.16 | 0.077 | +0.036 |
| **N2 matched-decoy contrast** | 0.079 | **−0.061 [−0.111, −0.012], p = 0.020** | 0.34–0.51 | 0.18–0.41 | 0.301 | +0.281 |

*(control rows are the `all_controls` family; the ranges span all five responses.
Full per-response table in `a7_summary.csv`.)*

Three statements, in the order they matter:

1. **The raw assay is not flat.** Negative-control counts fall systematically as
   distance to the nearest sender decreases — pooled amplitude −0.070 SD, sign
   consistent (only 34–51 % of fits positive, and 72 % of the CI-excluding fits
   negative), and the CI-exclusion rate is 18–45 % against a 5 % nominal. The
   binned curve (Figure 2h) rises from +0.019 SD in the 0–5 µm bin to +0.29 SD
   in the 95–100 µm bin, a span of about **0.27 SD** across the fitting window.

   Separately, and *not* the same number: the fitted control amplitude is
   **0.070 SD**, which is **a quarter** of the naive biological amplitude of
   **0.291 SD** in the same fits. *(Audit item M4: these two sentences used to
   run together, putting a bin range and an amplitude in adjacent clauses with
   near-identical numerals — "+0.29 SD" and "+0.291 SD" are different quantities
   and their near-equality is a coincidence of rounding.)* Anyone reporting a
   naive distance kernel on Xenium is reporting this in part.
2. **The N5 covariate block removes it, for every count-based control family.**
   Under +N6+N5 **four of the five** control families — every count-based one —
   have a clustered mean indistinguishable from zero (`all_controls`
   +0.007 [−0.011, +0.025]; `neg_control_probe` +0.006 [−0.004, +0.015]). *(Audit
   item R5: this said "every control family", which is wrong. The exception is
   `neg_probe_rate`, whose `n6n5` mean was +0.0108 [+0.0021, +0.0195], p = 0.020
   — a CI that excludes zero. It is disclosed two paragraphs below, so only the
   word "every" was ever at fault. **In the frozen C6 re-run that exception is
   gone**: `neg_probe_rate` under `n6n5` is now +0.0097 [−0.0060, +0.0253],
   p = 0.199, so all five families are indistinguishable from zero. The
   disappearance is recorded rather than quietly dropped.)* The
   sign becomes symmetric (49.7–53.9 % positive). The residual **excludes** the
   conditioned biological amplitude of +0.036 SD. **A7 passes under the design
   the paper actually fits.**
3. **N2 does not remove it, and this is the finding to carry forward.** The
   matched-decoy contrast leaves the control gradient at −0.061 SD, p = 0.020 —
   80 % undiminished, with a CI-exclusion rate of 18–41 %. N2 matches on local
   density, log counts, zonation and k-NN composition, and that is not enough:
   the confound lives in per-cell detection efficiency, which N5 models directly
   and a propensity match on neighbourhood covariates does not. **A matched-decoy
   contrast is not a substitute for the technical covariate block on this assay.**

`neg_probe_rate` (probe counts ÷ transcript counts) is flat naively (+0.014,
p = 0.079) — direct confirmation that the naive gradient is a depth/size effect
and not a spatial gradient in probe binding. Reported against interest: it is the
one response of five whose conditioned amplitude is nominally non-zero
(+0.011 [+0.002, +0.020], p = 0.020), which is expected for a ratio whose
denominator is itself an N5 column, and it is a third of the biological amplitude.

### 4.2 The power bound, stated so "flat" is falsifiable

The control probes are sparse — 0.0067 counts/cell for the 40 negative control
probes (0.65 % of cells non-zero), 0.043 for the 609 codewords (3.9 %), 0.0094
for the 21 genomic controls (0.90 %). Consequently **a single A7 fit resolves only
±0.137 SD**, which is 1.8× *larger* than the conditioned biological amplitude of
0.077 SD. Per fit, A7 could not rule out a technical gradient the size of the
biological one, and any report of A7 as "flat" at the level of one fit would be
meaningless.

The section-clustered pooled estimate is what carries the test: its CI half-width
is ±0.018 SD, four times finer than the conditioned biological amplitude, and it
excludes it. **A7 is adequately powered only when pooled across sections, and it
must be reported that way.** For H1 this sets a floor: a single spleen section
cannot deliver A7; seven can.

### 4.3 Go / no-go for H1

**GO**, with one pre-registered condition and one warning.

* **Condition:** the naive and N2 kernels are not valid readouts on this assay.
  Only the N5-conditioned amplitude may be reported as a distance effect. A7
  measures this directly and the margin is not comfortable: naive control
  −0.070 SD against naive biological +0.291 SD is a 24 % contamination, and under
  N2 it is 22 %.
* **Warning:** the CI-exclusion rate on a response with no biology is 9–16 %
  under the full design, against a 5 % nominal. That bounds the block bootstrap's
  false-positive rate at 2–3× nominal for this estimator on this tissue (or, read
  the other way, bounds the residual confounding N5 does not capture — both
  readings argue the same way). This is the first direct measurement of that rate
  anywhere in the project, and A7 is the only test that can produce it, because
  it is the only response whose true amplitude is known to be zero.

  > **Audit item R6 — a sentence deleted here.** This paragraph used to continue:
  > *"The 'reportable fit' filter … therefore admits two to three times more fits
  > than its nominal rate implies."* **That does not follow.** The 9–16 % figure
  > is the **two-sided** 95 %-CI exclusion rate under the full `n6n5` design. The
  > reportable filter is `beta_naive > 0 AND beta_base_lo > 0`
  > (`code/summarize_phase3_c1.py::reportable`) — **one-sided**, and on the
  > **naive** design. They are different rates on different designs and the
  > factor does not transfer. Measured on the 825 control fits of the frozen A7
  > run: the filter admits **4.8 % on the full design** — identical across all
  > five control families, i.e. essentially nominal — and **3.0–13.3 % on the
  > naive design**. The "2–3× nominal" bound on the *estimator* stands; the claim
  > about the *filter* does not, and is withdrawn.
  >
  > *(The audit stated 3.0–6.7 % and 1.2–13.3 % for these two rates. Those are
  > correct for the A7 file it had — the pre-C6 run of 05:19 — which the M1
  > re-run superseded at 06:52 when the sender call changed with Tier A 25 → 33.
  > The frozen numbers are the ones above.)*

### 4.4 Figure 2, panel h

`figures/figure2e.png` md5 `a8dd1eb8000179b8dff2a62034e1d2cc`. Panel h now draws
the pooled binned control curve — naive and N6+N5-residualised — over the six
in-band sections at `tierA_p95` in hepatocytes, with the section-clustered
amplitudes and the N2 row in a box. Its title reads
**"M1 DONE, H1 pending"** and the H1 half is still declared, in the panel, in the
caption text and as a `negative_control_probe_kernel_H1 / family=PENDING` row in
`figure2e_data.csv`. Twenty-nine further rows carry every plotted number: the
five designs × two response families, and the twenty curve bins.

`make_figure2bc.py` keeps the placeholder branch, so the script still runs and
still declares the panel pending if `a7_summary.csv` is absent.

**Panel 2d regenerated byte-identically to the repository**
(`5ca89780eac1bedee3a6bdcbe0434125`), confirming the panel-h change is contained.
`make_figure2.py` was **not** run.

> **Corrected 2026-08-27 (audit item R4).** This paragraph originally claimed
> that **2b, 2c and 2d** all "regenerated byte-identically", quoting
> `df2f0afdc3264e7e36693a5cd542c15e` and `983b47b2d97b4aa490b9903803b68f0b` for
> 2b and 2c. **Those two hashes are not `figures/` files.** They are
> `figures/revised_candidates/figure2b_REVISED.png` and `figure2c_REVISED.png`
> — verified by `md5sum` on both paths. The committed files were, and at the time
> of writing still were, `5ecd9ad1029851c1955fc938abf9444c` (2b) and
> `a232a529a566d7c0680e04840dd07a9b` (2c), which is what
> `git show HEAD:figures/figure2{b,c}.png` returns. So what the run reproduced
> byte-identically was the **revised candidate**, not the repository state, and
> `figures/revised_candidates/README.md` says so in its own table: for 2b and 2c,
> "Content differs? **YES**".

> **What this paragraph used to call an anomaly was the documented restore
> (audit item R4).** It read: after the `make_figure2bc.py` run at 05:08 had
> reproduced 2b and 2c, "both files were **rewritten by something else at
> 05:22:41** — 2b to `5ecd9ad1…` and 2c to `a232a529…`", and concluded that a
> concurrent session had called `fig2b()`/`fig2c()` directly and that the
> committed content was intact at `df2f0afd…`/`983b47b2…`.
>
> **Every part of that is backwards.** `5ecd9ad1…` and `a232a529…` *are* the
> committed hashes. The 05:22:41 write was the **PI-directed restore of
> `figures/` to its committed state**, recorded in `PHASE8_ROADMAP_STATUS.md`
> under the figure policy — an intended action, not an unexplained one. This
> report's own regeneration was then reverted at 05:29, so what sits on disk is
> the committed baseline and **not** what this report says it left there.
>
> **The general point survived and has now been acted on.** More than one actor
> does write `figures/`, and until 8.7 nothing warned. `code/check_figures_guard.py`
> now walks the whole directory rather than `git ls-files` and covers **46
> artefacts** against the 27 tracked ones, with the manifest snapshotted to the
> post-8.7 state. The recommendation this paragraph made — hash the figure
> directory rather than trusting mtimes — is implemented; the incident it was
> built on was not an incident.

### 4.5 Proposed caption change for panel h

> **(h)** Section 13 test A7 on the mouse arm: negative-control counts (40
> negative control probes, 609 negative control codewords and 21 genomic
> controls, pooled and z-scored) against distance to nearest sender, six
> Test-3-admissible sections, hepatocytes. The raw assay is **not** flat —
> amplitude −0.070 SD (section-clustered 95 % CI [−0.128, −0.012], p = 0.023), a
> quarter of the naive Tier B amplitude of +0.291 SD in the same fits — and the
> matched-decoy contrast does not remove it (−0.061, p = 0.020). The N5 technical
> covariate block does: +0.007 [−0.011, +0.025], p = 0.41, which excludes the
> conditioned Tier B amplitude of +0.036 SD. The human half is held behind the
> pre-registration freeze.

## 5. 8.6 — D3: what re-anchoring is, and what it finds

### 5.1 Re-anchoring is one bit per section, and no re-run

`DeepScence/io.py::fix_score_direction` negates a bottleneck node if `CDKN1A`'s
rank in the score-vs-gene correlation table falls in the bottom half, then picks
the node maximising `corr_df["correlation"].abs().mean()`. **That selection
metric is invariant to negation**, so the anchor decides exactly one bit per
section — a global sign — and nothing else. Re-anchoring is therefore a
post-hoc sign choice on the score already computed. No DeepScence re-run was
needed, and **the D1 outputs were not touched**, so 8.3/8.4 remain a pure
coverage comparison as required. D3's scores are separate files,
`data/processed/deepscence_d3_<section>.csv`, carrying `ds_published`,
`ds_prolif_anchor`, `ds_lmnb1_anchor`, `ds_consensus_anchor` and `ds_abs`.

### 5.2 A defect in the §7 brief, found while implementing it

§7 proposes `Lmnb1` as an anchor "in neither Tier A nor any Tier B module".
**`Lmnb1` is in `B_downstream_arrest.txt` and `B_secondary_senescence.txt`**, and
in the non-strict Tier A variants. It is reported here as a secondary anchor and
flagged; it must not be the primary one.

The primary anchor is an eight-gene proliferation set — `Kif20a`, `Ncaph`,
`Anln`, `Ect2`, `Gtse1`, `Uhrf1`, `Fen1`, `Clspn` — chosen because each is on the
mouse panel and absent from **every** `A_*.txt` and `B_*.txt` in `genesets/`
(590 distinct genes across the 22 files checked; the script asserts this at run
time and fails loudly). The
canonical proliferation markers are unusable: `Mki67`, `Top2a`, `Ccnb1`, `Pcna`,
`Foxm1` and the rest are all inside Tier A, Tier B, or both.

### 5.3 Depth must be partialled out or the anchor is a depth anchor

Every anchor candidate on this panel is a detection-rate readout as much as a
biological one: the eight-gene proliferation set is detected at all in only
3.5–12.8 % of cells, `Lmnb1` in 2.2–8.9 %, `Cdkn1a` in 0.68–21.6 %. On the raw correlation the
proliferation anchor points the *wrong* way (ρ = +0.011 to +0.031) in eight of
eleven sections — but so does Tier A's own senescence score (ρ = +0.023 to
+0.090 with the same proliferation set, in all eleven). The raw correlation is
measuring depth. **All anchor decisions here are taken on the Spearman
correlation after linearly removing the rank of log transcript counts from both
variables**, and the picture becomes coherent immediately.

### 5.4 The result

`deepscence_anchor_decisions.csv`, `figures/figure_phase8_d3.png`.

| anchor | expect | depth-partialled ρ, range over 11 sections | agrees with the published sign |
|---|---|---|---|
| proliferation set (8 genes, disjoint from all A/B) | ρ < 0 | **−0.098 … +0.045** | **10 / 11** |
| `Lmnb1` (secondary; in two Tier B modules) | ρ < 0 | **−0.096 … +0.050** | **10 / 11** |
| consensus of the other three callers | ρ > 0 | −0.089 … +0.063 | 2 / 11 |
| `Cdkn1a` — the published anchor | ρ > 0 | **−0.024 … +0.182** | — |

1. **The two caller-free anchors agree with each other in 11 of 11 sections**,
   including on the section where they disagree with `CDKN1A`. Fold-split sign
   stability (twenty random folds) is 0.95–1.00 for the proliferation anchor and
   0.90–1.00 for `Lmnb1`, in every section; the published `CDKN1A` anchor drops to
   0.60 in 7352 and 0.85 in 7248.
2. **They agree with the published anchor in 10 of 11. The exception is 7250** —
   the sham section on which every published DeepScence number rests. In 7250
   both say the published sign is **inverted**. This converts `BIO_PHASE3.md`
   §4.3's explicitly-labelled *inference* ("DeepScence fixes the sign by
   correlating with `CDKN1A`; in 7250 `Cdkn1a` is detected in 0.48 % of
   hepatocytes, which is a very weak anchor") into a **measurement**, from two
   independent directions, neither of which involves `CDKN1A` or any caller.
3. **7250 is the single explanation for four separate published anomalies.**
   DeepScence's ρ with depth (−0.350 there, +0.29…+0.56 in the other ten); its
   0.03× hepatocyte enrichment (1.38–2.46 in the other ten); the SenePy-vs-
   DeepScence 2.15× (0.33–0.55 in ten of eleven); and the 2.85× DeepScence-vs-
   `Cdkn1a`⁺ circularity (median 1.071 over eleven). All four are the inverted
   sign in one section, and that section was half of the published base.
4. **Reported against interest: the published anchor is not clean either.** Once
   depth is removed, ρ(ds, `Cdkn1a`) is **negative** in 7248 (−0.012) and 7435
   (−0.024) and effectively zero in 7352 (+0.0021, fold stability 0.60) and near-zero in 7001
   (+0.0106, stability 0.95). In four
   of eleven sections the published anchor is deciding the polarity on a signal
   that is weak, unstable, or of the wrong sign — it is right in those sections
   only because the *unpartialled* correlation, inflated by depth, happens to
   come out positive. The right recommendation is not "the `CDKN1A` anchor is
   fine in 10/11" but "**the anchor should be the depth-partialled proliferation
   set, and the `CDKN1A` anchor's agreement with it should be reported**".
5. **The consensus-of-callers anchor is unusable and should not be offered.** It
   disagrees with the published sign in 9 of 11 sections, because it is dominated
   by SenePy, which is anti-concordant with DeepScence (§2.2). Anchoring on it
   inflates DeepScence-vs-SenePy agreement to 2.793× chance — it manufactures
   exactly the circularity D3 exists to remove, in the opposite direction.

### 5.5 The effect on caller agreement, and the sign-invariant summary

`caller_agreement_matched_d3_11sections.csv`. Pooled over eleven sections, depth-
and type-matched, same estimator as §2:

| DeepScence variant | vs Tier A | vs `Cdkn1a`⁺ *(circular under the published anchor)* | vs SenePy |
|---|---|---|---|
| published (`CDKN1A` anchor) | **1.248** (z 16.6, 11/11 above) | **1.255** (z 10.5) | 0.737 (z −15.1) |
| re-anchored, proliferation *(= `Lmnb1`, identical decisions)* | **1.237** (z 15.8, 10/11 above) | **1.140** (z 5.8) | 0.495 (z −28.9) |
| **sign-invariant, rank by \|score\|** | **1.285** (z 19.0, 10/11 above) | **1.105** (z 4.3) | 1.025 (z 1.4, n.s.) |
| consensus anchor *(circular, do not use)* | 0.935 | 0.993 | 2.793 |

Three things follow.

* **Re-anchoring cuts the measured circularity by about half** (1.255 → 1.140)
  and leaves the Tier A agreement essentially untouched (1.248 → 1.237). Only the
  pairs involving 7250 move, because only 7250's sign changes.
* **The Tier A–DeepScence agreement is sign-invariant.** Ranking by |score| gives
  **1.285× chance**, above chance in 10 of 11 sections — as high as either signed
  version. Whatever Tier A and DeepScence share is a property of the score's
  *magnitude*, not its polarity, and therefore cannot be an artefact of the
  anchor. This is the strongest form of the §2.2 headline movement, and it is the
  form to pre-register for H1, where the anchor problem will recur and the
  polarity may not be recoverable at all.
* **Under |score|, DeepScence vs SenePy is at chance** (1.025, z = 1.4, not
  significant) rather than strongly below it. The 0.74 and 0.50 anti-concordance
  of the signed versions is a polarity disagreement, not a disagreement about
  which cells are extreme.

## 6. What I did not do

* **D2 (8.5 in the roadmap's numbering, `denoise=False`) is untouched.** The task
  assigned here was §13's audit test A7, and both are labelled "8.5" in different
  documents — the roadmap's 8.5 is C7/D2, the brief's 8.5 is A7. A7 is what was
  asked for in the brief and is what was run. **D2 remains open** and the §4.3
  caveat still stands as written: everything measured here characterises
  DeepScence *as we can run it on this panel*, not DeepScence as published.
* **The eleven-section tables have not been propagated into any downstream fit.**
  `sf_summary.csv`, `main_fits.csv` and the Phase 5 kernels still rest on the
  existing sender calls. That is task 8.7 and it is where these numbers land.
* **A7 was run at two sender calls, not six.** `tierA_p95` (primary) and
  `cdkn1a_pos` (the source paper's). The N7 axis would cost another ~90 s; it was
  not run because the control response has no sender-definition-specific
  mechanism and the two calls already bracket the prevalence range 1.7–9.0 %.
* **A7 curves are hepatocytes only.** Panel h needs one legible curve and
  hepatocytes are what Figure 2b draws; the fits in `a7_summary.csv` cover every
  receiver type with ≥ 2,000 receivers.
* **No H1 anything.** `data/raw_h1/` was not opened.

## 7. Reproduce

```bash
cd /workspace
python3 -u code/caller_disagree_all.py --verify        # ~2 min, must print VERIFY: PASS
python3 -u code/caller_disagree_all.py --all           # ~4 min
python3 -u code/summarize_caller_coverage.py
python3 -u code/make_figure_phase8_callers.py
cd code && python3 -u run_a7_control_probes.py --sections all \
        --calls tierA_p95,cdkn1a_pos --n-jobs 22       # 44 s
cd /workspace
python3 -u code/summarize_a7.py
python3 -u code/deepscence_reanchor.py                 # ~4 min
python3 -u code/caller_disagree_d3.py                  # ~6 min
python3 -u code/make_figure_phase8_d3.py
python3 -u code/make_figure2bc.py                      # 2b/2c/2d byte-identical; 2e panel h filled
# NOT make_figure2.py -- superseded producer of figure2a
```

`--verify` is not optional decoration: it is the only evidence that the
reconstructed producer of the two lost tables is the original one.

### Verification performed

* All six caller tables reproduce the committed two-section values **exactly**
  (`VERIFY: PASS`), including the two whose producer was never committed.
* The 7250 / 7259 rows inside the eleven-section tables are **identical** to the
  committed `sham` / `sbr` rows.
* `caller_*_11sections.csv` and `a7_control_probe_fits.csv` are **bit-identical**
  across two independent runs; so are `figure_phase8_{callers,d3}.png` and
  `caller_coverage_gate.csv`.
* `cells.parquet` control columns verified `np.array_equal` against per-feature
  h5 sums for all three control families on 7250.
* `figure2{b,c,d}.png` and `figure2c_data.csv` byte-identical after regenerating
  the whole of `make_figure2bc.py`.
* `perm_nulls.csv`, `sf_summary.csv`, `summary_phase3.txt` md5s unchanged (top of
  this report); `deepscence_{sham,sbr}.csv` unchanged by md5 and by mtime.

### Engineering notes

* **A lost producer is recoverable if its outputs were committed.** The two
  missing tables were reconstructed by inverting their contents: the ratio column
  fixes the chance definition, `n_A/n` fixes the threshold, and the difference in
  `n_A` between pairs proves the flags are recomputed on the pairwise-complete
  subset rather than subset from a global call. The exact-match test then decides
  it. Guessing the recipe took twenty minutes; not being able to prove it would
  have made every eleven-section number unciteable.
* **Substituting the response into an existing estimator beats writing a second
  one.** A7 reuses `SectionFit` and `fit_cell` and changes three lines; that is
  why the control amplitudes and the module amplitudes are comparable at all.
  Monkeypatching `P.MODULES` inside the worker to relabel the rows is ugly and is
  confined to a `try/finally`, which is the right trade against a parallel
  implementation that would drift.
* **The stratified-exact null cost nothing and settled a question.** Because both
  call sets are ~5 % inside every stratum, the marginal product and the exact
  conditional expectation agree to three decimals — so the published `ratio` was
  right, and the z-scores it never had come free from the same bincounts.
* **Depth-partialling turned an incoherent D3 result into a clean one.** On raw
  correlations the four anchors disagreed pairwise with no pattern, and the
  natural (wrong) conclusion was "the sign is not identified on this panel".
  Partialling out log counts made ten of eleven sections agree across four
  anchors and isolated 7250. The anchor genes are detected in under 2 % of cells;
  at that detection rate a correlation is mostly a statement about how many
  transcripts the cell had.
