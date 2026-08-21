# CS Phase 1 — Synthetic Identifiability Study (Master Plan Section 22, Step 1)

**Status: complete. The Day 3 gate (Section 26) is met.**
1,270 simulated tissue sections across six sweeps; 30 seeds per grid cell in the main
grid; all seeds pinned; Figure 1 produced.

| Deliverable | Location |
|---|---|
| Generator | `/workspace/code/sasp_sim.py` |
| Estimators (naive, N2 matched-decoy, N5 nuisance) | `/workspace/code/sasp_estimators.py` |
| Sweep runner (checkpointed, seed-pinned) | `/workspace/code/sasp_sweep.py` |
| Figure 1 + supporting scripts | `/workspace/code/make_figure1.py`, `make_curves.py`, `summarize.py`, `sasp_palette.py` |
| Per-block checkpoints (41 blocks) | `/workspace/results/sweep/` |
| Tidy results, 1,270 runs x 115 columns | `/workspace/results/sweep_all.csv` |
| Summary tables | `/workspace/results/summary_tables.txt` |
| **Figure 1** | `/workspace/figures/figure1.png`, `.pdf` |
| Figure 1 underlying data | `/workspace/figures/figure1_data.csv` |

Simulation constants: λ_true = 30 µm, β_true = 1.0, 2000 x 2000 µm window,
**n ≈ 11,900 cells, median nearest-neighbour distance 11.7 µm** (inside the 10–20 µm
cell-diameter range of Section 3, so λ_true sits comfortably above the resolution
floor of Section 8 Test 1). Sender prevalence 5 % unless swept.

---

## 1. Does recovery work in the easy regime? Yes — but only when the regime is
## *actually* easy, and "random senders, no latent confounder" is not enough.

Two clean configurations, 40 seeds each:

| configuration | λ̂_naive (mean ± sd) | β̂_naive | coverage (block CI) |
|---|---|---|---|
| **PURE** — no nuisance at all: no GRF baseline, no density effect on response, no count effect, senders a random subset | **30.53 ± 4.46** vs λ_true = 30 | **1.008** vs β_true = 1.0 | **λ 0.95, β 1.00** |
| **EASY** — random senders (κ = 0), no latent confounder, ℓ/λ = 0.25, but the ordinary nuisance terms present | **43.55 ± 9.47** (+45 %) | 1.184 (+18 %) | λ 0.60, β 0.30 |

**The estimator is correct.** When the planted kernel is the only structure in the
data, λ̂ = 30.5 ± 4.5 against a truth of 30, β̂ = 1.008 against 1.0, and the block
bootstrap CI covers at 0.95 and 1.00. Nothing downstream is an artefact of a broken
fit.

**But "easy" in the sense the plan means it is already biased by +45 %.** The reason
is a confound that needs no latent field at all: **local cell density mechanically
confounds distance-to-nearest-sender.** In a denser neighbourhood the nearest sender
is closer, and density also shifts the response. Random senders and a short baseline
correlation length do not protect you from it. Nuisance conditioning (N5) removes
essentially all of it in this regime — λ̂_nuis = 29.76 [28.30, 31.23], β̂ = 1.00
[0.99, 1.02]. **The matched-decoy control does not** (λ̂ = 44.85, β̂ = 1.18). Note this is regime-specific: at
conf_strength = 0 with clustered senders, N5 is the *worst* of the four (Section 2).

This is a finding the plan did not anticipate and it applies to every real dataset.

---

## 2. The regime map (Figure 1a, 1b)

Rows: sender clustering κ (realised label-permutation Ripley ratio at 50 µm in
parentheses). Columns: baseline autocorrelation length ℓ / λ_true.

**Relative bias in λ̂ (naive):**

| κ (Ripley) | ℓ/λ=0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| 0.0 (1.1–3.0) | −1 % | −7 % | **+22 %** | **+58 %** | **+114 %** |
| 1.5 (4.2–7.1) | −11 % | −26 % | −13 % | **+25 %** | **+141 %** |
| 3.0 (8.8–10.8) | −1 % | −26 % | −26 % | −6 % | **+48 %** |
| 4.5 (10.9–12.6) | +10 % | −10 % | −22 % | −13 % | **+202 %** |

**95 % CI coverage of λ, naive estimator with the iid asymptotic CI** (nominal 0.95):

| κ | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| 0.0 | 0.93 | 0.67 | 0.53 | **0.00** | **0.00** |
| 1.5 | 0.60 | 0.13 | 0.40 | 0.37 | 0.10 |
| 3.0 | 0.83 | 0.33 | 0.27 | 0.50 | 0.33 |
| 4.5 | 0.73 | 0.47 | 0.40 | 0.33 | 0.27 |

Mean coverage over the whole grid: **naive iid 0.41, naive block bootstrap 0.67,
matched-decoy block bootstrap 0.63, nuisance-conditioned block bootstrap 0.71.**

### The headline the plan predicted: **confirmed for the first clause, refuted for the second.**

* **"Naive estimation is badly biased once ℓ approaches or exceeds λ_true" — CONFIRMED,
  strongly.** Averaged over κ, λ̂ goes 29.8 → 24.8 → 27.0 → 34.7 → 67.9 µm as ℓ/λ goes
  0.25 → 4. At ℓ/λ = 4 the naive estimate is roughly **double** the truth.
* **"…and overconfident" — CONFIRMED.** Coverage of the standard iid CI falls
  monotonically with ℓ: 0.78 → 0.40 → 0.40 → 0.30 → **0.18**. The direct measure is the
  **SE understatement factor** (block-bootstrap sd ÷ iid asymptotic SE), which rises
  from 1.5 to **7.9**: a standard nonlinear-least-squares CI is up to eight times too
  narrow. Bias and over-confidence are separate failures — switching to a spatial
  block bootstrap recovers coverage from 0.41 to 0.67 on average, but cannot fix the
  bias, and at ℓ/λ ≥ 2 with κ = 0 even the block bootstrap covers 0.00.
* **"…and the decoy control restores approximate calibration" — NOT CONFIRMED.**
  The matched-decoy control is nearly inert on λ̂. Its bias grid is almost identical
  to the naive one (+99 % vs +114 % at κ=0, ℓ/λ=4; +227 % vs +202 % at κ=4.5), and its
  mean coverage (0.63) is no better than the naive block bootstrap (0.67) and *worse*
  at ℓ/λ = 4 (0.18 vs 0.42).

### The actual headline: **no single control is safe across regimes**

An earlier draft of this report concluded that "N5, not N2, is what rescues λ̂."
**That was wrong, and the correct statement is both stronger and more useful.**
N5 does not remove the bias; it *replaces a positive bias with a negative one*, and
which one you get depends on a regime parameter that is not observable in real tissue.

Mean relative bias in λ̂, pooled over κ:

| estimator | ℓ/λ=0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| naive | −1 % | −17 % | −10 % | +16 % | **+126 %** |
| N2 matched-decoy | 0 % | −14 % | −2 % | +34 % | **+141 %** |
| N5 nuisance-conditioned | −12 % | **−29 %** | **−26 %** | −16 % | +10 % |
| N2 + N5 | −12 % | **−29 %** | **−27 %** | −16 % | +11 % |

Three things follow, and none of them is "use N5":

1. **The best estimator alternates.** Ranking by *absolute* pooled bias, the winner at
   ℓ/λ = 0.25, 0.5, 1, 2, 4 is N5, N2, N2, N2+N5, N5 respectively. No control dominates
   at more than two of the five settings.
2. **N5 can be actively worse than doing nothing.** In the confounder-strength sweep at
   **conf_strength = 0** — no latent confounder at all — N5 returns λ̂ = **47.6 µm**
   against a truth of 30, while the naive estimator returns 36.7 and N2 returns 40.9.
   Conditioning on density and kNN composition over-adjusts when those covariates are
   partly *downstream* of proximity to a sender (a cell near a sender has that sender in
   its own 50 µm neighbourhood, so density and composition are mediators as well as
   confounders). Real Tier-D covariates have exactly this defect.
3. **Even the spread across all four estimators is not a valid uncertainty interval.**
   Across the 20 main-grid cells the four estimates span a median of 8.0 µm
   (range 1.8–60.1 µm), and that range **brackets λ_true in only 25 % of cells.**

**Recommended practice, and the claim this phase supports:** report naive, N2, N5 and
N2+N5 **side by side and treat the spread as the reportable uncertainty**, while stating
explicitly that the spread is a *lower bound* on uncertainty — it is not calibrated and
covered the truth only a quarter of the time here. A single "controlled" number, from
any of these controls, is not defensible from a static snapshot. For an identifiability
paper this is a better result than "use N5": the finding is that **every available
control is regime-dependent, and the regime is not identifiable from the data.**

I am reporting this as it came out. The simulation was not re-tuned to agree with the
plan's expectation, and the erroneous "use N5" reading was corrected against the same
tables rather than defended.

### Why the decoy control under-performs — a diagnosable mechanism
The synthetic data lets us audit the latent confounder `u` directly (no estimator ever
sees it). Averaged over the main grid, propensity matching on the observable covariates
closes only **29 %** of the sender-to-population gap in `u`. Matched decoys sit in
neighbourhoods that *look* like sender neighbourhoods on density, kNN composition and
counts, but are still substantially less extreme on the thing that actually drives the
response. The decoy kernel therefore absorbs only part of the confound.

**Crucially, this happens while matching passes the plan's own quality bar.**
Max SMD after matching is 0.05–0.07 across the entire grid, and **95.2 % of runs pass
Section 8 Test 5 (max SMD < 0.1)**. *Passing Test 5 is not sufficient for the decoy
control to work.* That is a caveat the plan should carry into the real-data analysis:
SMD balance on observables says nothing about balance on unobservables.

### Where λ becomes unidentifiable
* **ℓ/λ_true ≥ 2** — the primary failure axis. λ̂ is pulled *toward* ℓ (biased low when
  ℓ < λ_true, high when ℓ > λ_true), because the confounder's own spatial scale sets the
  scale of the spurious bump around senders. At ℓ/λ = 4 no estimator except N5 is usable
  and no CI construction covers.
* **Sender prevalence below ~1 %** — at 0.5 % prevalence (59 senders) the between-seed
  sd of λ̂ is **16.8 µm** on a truth of 30, versus 4.5 µm at 2–5 %. λ is not so much
  biased as simply unresolvable. This reproduces Section 8 Test 3's lower bound.
* **Sender prevalence above ~20 %** — median distance-to-nearest-sender falls to 24 µm
  at 20 % and 17 µm at 30 %, i.e. to about the median NN cell distance, and λ̂_naive
  climbs to 32.6 and 37.3. Dynamic range in the independent variable disappears, again
  as Section 8 Test 3 predicts. **The plan's 2–10 % sweet spot is confirmed**, though
  note λ̂_naive is mildly biased *low* (−15 %) throughout it under this confounder.
* **Strong sender clustering does *not* by itself destroy λ̂** — it inflates variance and
  degrades coverage, and at κ = 4.5 it interacts badly with long ℓ (+202 %), but κ alone
  at short ℓ is comparatively benign. Clustering matters mainly through the
  *effective sample size*: at κ = 4.5 the senders occupy a handful of foci.
* **N is not the binding constraint.** Across 2,950 → 23,500 cells λ̂ and coverage are
  flat (0.90/0.95/0.75/0.95). More cells do not buy identifiability; they buy precision
  around a biased point.

---

## 3. Coverage numbers, collected (Section 24.7)

Mean over the 20-cell main grid, 30 seeds each:

| CI construction | coverage of λ | note |
|---|---|---|
| naive, iid asymptotic (the field standard) | **0.41** | falls to 0.18 at ℓ/λ = 4 |
| naive, spatial block bootstrap | 0.67 | fixes the SE, not the bias |
| matched-decoy (N2), block bootstrap | 0.63 | no better than naive |
| nuisance-conditioned (N5), block bootstrap | 0.71 | 0.90 at ℓ/λ = 4 |

Coverage of β is far worse than coverage of λ. Under the **β_true = 0 null** — no
planted kernel whatsoever — the naive estimator returns β̂ = 0.63 to 1.86 depending on
regime, i.e. **an effect as large as a real one, in 90–100 % of replicates**, and its
block-bootstrap CI covers the truth (zero) only 3–27 % of the time. The matched-decoy
control reduces but does not eliminate this (β̂ = 0.62–1.07). In the *clean* null
(κ = 0, no latent confounder) nuisance conditioning brings β̂ down to 0.11 while the
decoy control leaves it at 0.62. **β is the fragile parameter; λ is comparatively
robust.** Any real-data claim should lead with λ and treat β as provisional.

---

## 4. How much true signal does the decoy control cost? (Section 29, objection 6)

Reported as β̂_decoy-corrected ÷ β̂_naive, with a real planted kernel present:

| regime | β̂ naive | β̂ decoy-corrected | fraction of naive effect retained |
|---|---|---|---|
| clean, no nuisance at all | 1.01 [0.99, 1.02] | 1.01 [0.99, 1.02] | **0.999** |
| clean, no latent confounder | 1.18 [1.15, 1.22] | 1.18 [1.15, 1.22] | **1.001** |
| ℓ/λ = 0.25 (confounder on) | 1.13 [1.09, 1.17] | 1.13 [1.09, 1.17] | 0.996 |
| ℓ/λ = 1 | 1.86 [1.79, 1.92] | 1.70 [1.63, 1.77] | 0.911 |
| ℓ/λ = 4 | 2.07 [1.94, 2.20] | 1.75 [1.64, 1.87] | 0.844 |

**Answer: essentially nothing — the decoy control removes 0.0–0.1 % of a genuine
planted effect when there is no confounding, and 9–16 % when there is.**
The decoy control is **not** over-conservative in the two-kernel regression form.
Its problem is the opposite: it is **under-powered**, removing far less confounding
than is present (β̂ stays at 1.75 when the truth is 1.0).

### A methodological caveat on "the single most important number in the paper"
The plan prescribes reporting **β_true − β_decoy** from two separate fits. That scalar
is **not safe as written**. When the two fits are free to choose their own length
scales, the decoy fit runs off to λ_decoy ≈ 240–300 µm in long-range-confounded
regimes, and subtracting a long-range amplitude from a short-range one is
apples-to-oranges: the number comes out at **−1.44 [−1.76, −1.12]** at ℓ/λ = 4, i.e.
a strongly *negative* "corrected effect" that is an artefact of scale mismatch, not
over-conservatism. Constraining λ_decoy = λ_sender (implemented as
`BlockProfiler.fit2_shared`) makes the contrast interpretable and is what the numbers
above use. **Recommendation: report the shared-scale two-kernel β_s, and report
β_true − β_decoy only alongside both fitted λ's so the reader can see whether the
subtraction is meaningful.**

---

## 5. Engineering notes (Section 18 compliance)

* `scipy.spatial.cKDTree` for all neighbour work; **no (n,n) distance matrix anywhere**.
* Decoy matching runs in 1-D propensity space by sort + `searchsorted`, O(n log n).
* The block bootstrap exploits the fact that a replicate is an integer multiplicity per
  block and every quantity the fits need is a cross-product: per-block Gram matrices are
  computed once, and a replicate is a weighted sum of 64 small matrices plus a
  vectorised 2-D grid solve. **A 400-replicate block bootstrap of a two-kernel profiled
  nonlinear fit costs ~0.2 s instead of ~30 s.** This is what made 1,270 runs x 400
  bootstrap replicates x 4 estimators affordable (~9 min wall on 46 cores).
* BLAS pinned to one thread per worker; parallelism across runs (`joblib`, n_jobs=46).
  Multithreaded BLAS on these small matrices was a 10x *slowdown* from thread contention.
* Every (sweep, config) block checkpointed to `/workspace/results/sweep/` and skipped on
  restart. This was exercised for real: a `ZeroDivisionError` on the β_true = 0 sweep
  killed the run after 35 blocks, and the resume re-used all 35 checkpoints in 0.3 s.
* All seeds pinned from `MASTER_SEED = 20260820` plus (sweep id, config id, replicate).
  The entire figure is reproducible from `sasp_sweep.py`.
* Total disk written: ~5 MB. Nothing large was written to `/workspace`.

Two implementation choices worth flagging in Methods:
1. The binned estimator uses the **in-bin mean distance**, not the bin midpoint.
   Midpoints inject a ~20 % upward bias in λ̂ (Jensen, plus the near-empty 0–10 µm bin
   under a hard-core point process) that would otherwise be misread as confounding bias.
2. Cell-type intercepts µ_c are always in the design, since Section 6.1 makes them part
   of the model rather than a control.

---

## 6. What I did NOT get to

* **Other kernel families** (Section 6.2): Gaussian, power law, step/threshold,
  nonparametric spline. Only the exponential kernel is planted and fitted. Model
  *misspecification* is therefore untested — an obvious and cheap extension, and one a
  reviewer will ask about.
* **The superposition / aggregate-sender variant** (Section 6.3). Only
  nearest-sender is implemented. Under clustered senders these differ most, so this is
  the highest-value missing piece.
* **Proximal vs downstream programs** (Section 6.4).
* **Nulls N1, N3, N4, N6, N7, N8.** Only N2 and N5 are implemented. N3 (torus shift) is
  particularly worth adding to the synthetic study, since the generator is already on a
  torus and it is the null CellWHISPER showed leading methods fail.
* **Coarsened exact matching** is implemented (`match_method="cem"`) but was not swept;
  all reported results use propensity-score matching.
* **Anisotropic / anatomically structured tissue.** The confounder is an isotropic GRF;
  real tissue has vessels, zonation and boundaries (Sections 11 and 8 Test 6).
* **Multiple donors / sections and the donor bootstrap** (Section 24.1). Each run is one
  section; the block bootstrap stands in for the donor bootstrap.
* No real data was touched.

## 7. Recommended changes to the plan, based on these results

1. **Drop the promise that any single control restores calibration.** N2 is nearly
   inert on λ̂; N5 trades positive bias for negative bias and is *worse than naive* when
   there is no latent confounder; the best of the four alternates with ℓ/λ. Report all
   four side by side and present the spread as a lower bound on uncertainty. This should
   be the paper's stated position, not a caveat buried in the supplement.
2. **Add "SMD < 0.1 is necessary, not sufficient" to the Day 1 audit interpretation.**
   95 % of runs pass Test 5 and the decoy control still fails to remove the confounding.
3. **Report λ̂ as the headline and β̂ as provisional.** The β̂ = 0.63–1.86 result under a
   zero-effect null is the strongest single number produced in this phase.
4. **Add the baseline-autocorrelation length to the Day 1 audit.** ℓ/λ̂ is the axis that
   determines whether any real-data estimate is trustworthy, and it is measurable on
   real data (a variogram or Moran's-I length scale of the response residual). Without
   it there is no way to place a real dataset on Figure 1.

---

## Methods (Phase 1)

### Synthetic tissue generator (`/workspace/code/sasp_sim.py`)

**Cell placement.** Matérn type-II hard-core process on a 2000 x 2000 µm periodic
window: a dense Poisson proposal (0.020 points/µm²) is thinned by deleting any point
with a neighbour inside 9 µm carrying a smaller mark. A raw Poisson process produces
nearest-neighbour distances far below one cell diameter, which is not tissue. The
realised section has **n ≈ 11,900 cells and a median nearest-neighbour distance of
11.7 µm**, inside the 10–20 µm range Section 3 gives for a cell diameter and therefore
a resolution floor consistent with Section 8 Test 1. Neighbour work uses
`cKDTree.query_pairs` / `.query` / `.query_ball_point`; no (n,n) distance matrix is
constructed anywhere in the codebase.

**Random fields.** Latent fields are generated by spectral synthesis on a 1024²
periodic FFT grid using the Whittle–Matérn ν = ½ spectrum, so the covariance is
exp(-r/ℓ) and **ℓ is an e-folding length directly comparable to λ_true** of the
exponential kernel. Values are bilinearly interpolated to cell positions.

**The confounder.** A single latent field `u` simultaneously drives
(i) sender propensity (ψ·u), (ii) local cell density (cells are retained with
probability ∝ exp(0.55·u), so density genuinely varies with the niche),
(iii) cell-type composition (type logits load on u), and (iv) the baseline response
(η·u). `u` is never exposed to any estimator; density, kNN composition and counts are
its only observable proxies. `conf_strength` scales ψ and η together. This is
Section 22 Step 1.3(d) — "make the neighbourhoods around senders genuinely different
for reasons unrelated to signalling."

**Sender clustering.** `clustering = κ` mixes a Thomas cluster field (45 parents,
σ = 70 µm) into the log sender propensity; κ = 0 gives senders that are a random
subset. Exact prevalence is achieved by Gumbel top-k weighted sampling without
replacement. Realised clustering is reported as a label-permutation Ripley ratio at
50 µm — the same statistic Section 8 Test 4(b) prescribes for real tissue, so a real
dataset can be located on the Figure 1 axis.

**Response.** r_i = µ_{c_i} + β·exp(-d_i/λ_true) + η·u_i + 0.55·b_i
+ 0.25·z_dens,i + 0.30·z_counts,i + ε_i, with d_i the distance to the **nearest**
sender (Section 6.1), `b` an independent GRF with the same correlation length ℓ,
lognormal total counts, and heteroscedastic noise sd_i = 0.80·sqrt(T_med/T_i).
λ_true = 30 µm and β_true = 1.0 throughout.

### Estimators (`/workspace/code/sasp_estimators.py`)

* **naive** — cell-level profiled nonlinear least squares of
  r = X θ + β exp(-d_sender/λ), with X = intercept + receiver cell-type dummies (µ_c is
  part of the model in Section 6.1, not a control). λ is profiled on a 96-point log
  grid with the linear block solved exactly, then parabolically refined.
* **naive_bin** — the binned version the field plots: 10 µm bins to 300 µm
  (Section 22 Step 2), weighted by bin count. The **in-bin mean distance** is used as
  the regressor rather than the bin midpoint; midpoints inject a ~20 % upward bias in
  λ̂ from Jensen's inequality plus the near-empty 0–10 µm bin under a hard-core
  process, and that artefact must not contaminate a figure about confounding.
* **nuis (N5)** — the same fit with Tier-D style covariates: local density within
  50 µm, kNN cell-type composition, log total counts, distance to tissue boundary.
* **decoy (N2)** — for each sender, a non-sender of the **same cell type** matched on
  density-50, kNN composition and log counts by propensity-score nearest-neighbour
  matching, greedy 1-1 without replacement within cell type, caliper 0.25 sd. Matching
  runs in 1-D propensity space by sort + `searchsorted`, so it is O(n log n).
  Standardized mean differences are reported before and after (Section 8 Test 5).
  The corrected fit is a two-kernel regression
  r = X θ + β_s K(d_sender; λ_s) + β_d K(d_decoy; λ_d), so β_s is the sender effect
  over and above whatever appears around a matched non-sender. Two variants:
  **free λ_d**, and **shared λ_d = λ_s** (`decoyS`). The shared version is the
  apples-to-apples form of the plan's β_true − β_decoy; with both scales free, the
  decoy fit runs off to λ_d ≈ 300 µm in long-range-confounded regimes and subtracting
  a long-range amplitude from a short-range one is not a meaningful contrast.
* The plan's literal scalar **β_true − β_decoy** is also reported, from a separate
  binned fit using decoys as senders.

### Uncertainty
Two CIs per estimate: the **iid asymptotic Gauss–Newton CI** (what a standard NLS
package reports, and what the field would quote) and a **spatial block bootstrap**
(8 x 8 = 64 blocks of 250 µm, 400 replicates). Reporting both lets the coverage
failure be decomposed into bias versus understated standard error; their ratio is
the "SE understatement factor".

### Compute
Section 18 rules are hard requirements and are met. The block bootstrap exploits the
fact that a replicate is fully described by an integer multiplicity per block and that
every quantity the fits need is a cross-product: per-block Gram matrices are computed
once and each replicate is a weighted sum of 64 small matrices plus a vectorised 2-D
grid solve. **A 400-replicate block bootstrap of a two-kernel profiled nonlinear fit
costs ~0.2 s instead of ~30 s.** BLAS is pinned to one thread per worker and
parallelism is across runs (`joblib`, 46 workers). Every (sweep, config) block is
checkpointed to `/workspace/results/sweep/` on completion and skipped on restart, so a
crash costs one block. All seeds are pinned: every tissue and analysis is seeded from
(MASTER_SEED = 20260820, sweep id, config id, replicate), so the whole figure is
reproducible from `sasp_sweep.py` alone.

---

## 8. Reproducibility check across the numpy boundary (added 2026-08-21)

The README flagged, correctly, that `numpy` moved 1.26.3 → 2.4.6 partway through
the project and that Phase 1 therefore ran under a different numpy than the later
phases, with the seeded-reproduction question left open. On 2026-08-21 the
working container was reset and the stack was reinstalled from
`requirements.txt`, which made the check nearly free. The claim in the last
paragraph of Methods — *"the whole figure is reproducible from `sasp_sweep.py`
alone"* — has now been tested rather than asserted.

**Procedure.** `code/_sweep_repro.py` is `sasp_sweep.py` with `RESULTS` and the
aggregate path redirected to `results/repro_2026-08-21/`; nothing else differs.
All six sweeps (`main`, `prev`, `conf`, `clean`, `null`, `nsize`) were re-run
from scratch under numpy 2.4.6 / scipy 1.17.1 / pandas 2.3.3, then compared with
`code/_repro_check.py`.

**Result.**

| | |
|---|---|
| files compared | 43 / 43, none missing, none extra |
| shapes identical | 43 / 43 |
| column *sets* identical | 43 / 43 |
| column *order* identical | **8 / 43** |
| numeric values compared | 141,160 |
| **bit-identical** | **141,160 (100.0000 %)** |
| max relative difference | 0.0 |

Seeded reproduction holds exactly across the version boundary. `Figure 1` and
every number in §2–§4 above regenerate from `sasp_sweep.py` under either numpy.

**Two caveats found in the process, both worth keeping:**

1. **Column order is not stable between runs** (35 of 43 files). The result
   frames are assembled from per-run dicts, so ordering depends on which keys the
   estimator emitted first. Values are unaffected. Anything downstream must
   select columns **by name**; `results/sweep_all.csv` should never be indexed
   positionally.
2. **`results/sweep_all.csv` is 1 ulp off its own checkpoints**, in 0.82 % of
   values (max absolute 5.7×10⁻¹⁴, max relative 2.8×10⁻¹⁶). This is *not* a numpy
   effect. Rebuilding the aggregate from the stored per-config checkpoints
   reproduces the same 0.82 % / 5.7×10⁻¹⁴ discrepancy against the stored
   aggregate, i.e. the difference is internal to the original run: it crashed on
   the `null` sweep (`ZeroDivisionError`, `logs/sweep_2026-08-20_1537.log`) and
   was resumed, so some blocks were concatenated from memory and others re-read
   from CSV. The per-config files are the primary artefact and they are exact; no
   reported quantity moves at the fourteenth decimal place. Left as-is rather
   than regenerated, so that the stored aggregate remains the one the phase
   reports were written against.

Full per-file table: `results/repro_2026-08-21/sweep_file_diff.csv`.
