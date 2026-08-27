# CS Phase 3 — The full null battery on real tissue, and the GO/NO-GO

**Status: complete. Verdict: NO-GO on the length constant, GO on the paper.**
Master Plan Sections 22 Step 3, 23, 24, 25; Section 28 row 2 ("the naive decay does
not survive the nulls at all — **this is a result, not a failure**") is the row we
landed on, and it is now quantified.

| Deliverable | Location |
|---|---|
| Phase 3 data prep, estimator core, section admissibility | `/workspace/code/sasp_phase3.py`, `phase3_core.py` |
| Null battery runner (N1–N6, curves, permutations) | `/workspace/code/run_phase3_nulls.py` |
| N8 (disjointness, scrambled response, circularity) | `/workspace/code/run_phase3_n8.py` |
| T2 combined estimate + donor bootstrap | `/workspace/code/run_phase3_combined.py` |
| Density/Poisson identifiability test | `/workspace/code/run_phase3_poisson.py`, `run_phase3_lamscale.py` |
| Confounder attribution, composition decomposition | `/workspace/code/run_phase3_attribution.py`, `run_phase3_strat.py` |
| Summary tables | `/workspace/results/phase3/summary_phase3.txt` |
| **Figure 2b, 2c, 2d** | `/workspace/figures/figure2{b,c,d}.{png,pdf}` |
| All result tables (27 CSVs) | `/workspace/results/phase3/` |

---

## 0. Headline — the surviving fractions

Six sections, six mice, both surgical arms, nine receiver cell types × seven Tier B
modules, primary sender call `tierA_p95`. Reference population = the **160 of 315
fits that a reader of the current literature would actually report**: a positive
naive amplitude whose spatial block-bootstrap CI excludes zero.

Surviving fraction = β̂ under the null ÷ β̂ naive, both at the same λ̂ (Section 6.5;
CS Phase 2 §10 — never the permutation p-value).

| null | median SF | IQR | fraction ≤ 0 |
|---|---|---|---|
| **N3 torus shift** | **1.000** | [0.992, 1.008] | 0.00 |
| **N4 rotation** | **0.964** | [0.867, 1.042] | 0.00 |
| **N2 matched decoy** (shared-λ two-kernel) | **0.943** | [0.905, 0.980] | 0.00 |
| **N8 scrambled response gene set** | **0.935** | [0.819, 1.054] | 0.01 |
| zonation covariate alone (Section 11) | 0.843 | [0.510, 0.992] | 0.04 |
| **N1 cell-type-stratified label permutation** | **0.716** | [0.403, 0.882] | 0.09 |
| **N6 receiver-baseline conditioning** | **0.486** | [0.222, 0.756] | 0.12 |
| **N5 nuisance conditioning** | **0.084** | [−0.130, 0.261] | **0.38** |
| N5 + N6 | 0.082 | [−0.099, 0.252] | 0.35 |
| **N2 + N5 + N6 (the combined estimate)** | **0.082** | [−0.099, 0.249] | **0.35** |
| N1 applied *to the N5+N6-conditioned residual* | 0.987 | [0.921, 1.058] | 0.03 |

**Read the table top to bottom and the whole paper is in it.** Every null the field
treats as the strong one — torus shift, rotation, matched decoy — removes between
0 % and 6 % of the effect. Nuisance conditioning removes **92 %**, and in 38 % of
fits it drives the amplitude to zero or below. The three nulls that agree the
effect is real are the three that cannot tell confounding from causation.

**The kernel does not survive in either arm.** After N2 + N5 + N6 the median
amplitude is **0.027 response-sd at contact** (IQR −0.028 to 0.090, p90 0.176)
against a naive **0.326 sd**, and only **15 of 160** controlled fits have a
positive amplitude whose block bootstrap excludes zero.

**The bound.** The block-bootstrap SE of the controlled amplitude is 0.073 sd, so
this design has 80 % power at **0.203 response-sd**. We can therefore exclude a
controlled SASP amplitude larger than **~0.20 response standard deviations at cell
contact**, for any of 7 Tier B programs in any of 9 receiver cell types in mouse
liver at 5 % sender prevalence. We cannot exclude something smaller, and we cannot
report a length constant at all (§4).

---

## 1. What changed from the Phase 3 brief, and why

Two constraints in my brief were wrong and I re-scoped before finalising.

**(a) The arm restriction was wrong in both directions.** With all 11 sections
annotated (Bio Phase 3), Section 8 Test 3 admissibility does not follow the
surgical arm. `Cdkn1a`+ hepatocyte prevalence:

* **over the 20 % ceiling — excluded:** 7239 (45.0 %), 7448 (25.6 %), 7361 (25.3 %),
  7450 (22.7 %). Above the ceiling, distance-to-nearest-sender is near zero
  everywhere and λ is unidentifiable *by construction*.
* **in band — the primary set:** 7260 (10.5 %), 7259 (9.6 %), 7001 (8.9 %),
  7352 (7.2 %), 7248 (4.9 %), 7435 (2.3 %). Two SBR animals, four sham.
* **below the 1 % floor — excluded:** 7250 (0.48 %).

So four of six SBR sections fail the ceiling and four of five sham sections pass.
The primary analysis is the **six in-band sections from six animals across both
arms**, with arm as a contrast. That exclusion rule is a Methods sentence, and it
is worth stating that it *changes which arm the paper is about*.

The conclusion is not sensitive to it. Re-running the battery on the excluded
sections gives SF(N2+N5+N6) = **0.098** on the four over-ceiling sections and
**0.124** on the below-floor section, against 0.082 in band.

**(b) Receiver labels.** Bio Phase 3 §1.1 showed three label families are not
stably separable on this panel and the winner flips between sections (7001 lost its
entire mesenchymal compartment to `Unknown` on a 0.02 margin). All Phase 3 fits use
`cell_type_merged`. `Low_quality` and `Unknown` are excluded from receivers, and
`Proliferating` from sender calls, throughout.

**(c) Section-level detection depth.** Within SBR, section-level `Cdkn1a`+ burden
tracks section-level median transcripts per cell at ρ = **+0.943** (p = 0.005,
n = 6); I reproduce that exactly. Across all 11 sections it is +0.345 (p = 0.30) and
within sham +0.500 (p = 0.39), so it is an SBR-arm phenomenon in this cohort.
With one section per animal, section-level depth is perfectly collinear with the
section indicator and cannot be identified alongside section fixed effects; what I
did instead is (i) put **cell-level** transcript counts, genes detected, cell area
and nucleus area in the N5 block — where they turn out to be the single largest
confounder (§3) — and (ii) check the section-level association directly: over the
six in-band sections, median naive |β̂|/sd vs section median depth gives
ρ = +0.143 (p = 0.79). The depth confound acts within sections, not between them.

---

## 2. The window, fixed first

CS Phase 2 §2 showed the plan's "10 µm bins out to 300 µm" is unreachable. Measured
again on the admissible set at the primary call, pooled over 848,596 receivers:
median 32.6 µm, p90 64.7, p95 75.6, **p99 98.3**, and 99.8 % of receivers within
120 µm.

**Window = 100 µm** (99.1 % of receivers retained; 96.5 % at `cdkn1a_pos`).
**λ grid = [7, 50] µm, 40 log-spaced points.** The floor is the resolution floor of
Section 8 Test 1 — the median nearest-neighbour distance is 6.7–10.6 µm across
sections, so nothing below ~7 µm is a reportable length scale. The ceiling is
window/2, because an exponential with λ > d_max/2 is not distinguishable from a
linear trend over the observed range. Everything below is reported with the fraction
of fits that rail against those bounds, which is the honest way to say "not
identified".

---

## 3. Where the naive gradient actually comes from

### 3.1 More than half of it is receiver cell-type composition

Fitting **unstratified**, exactly as Phase 2 did and as the field does, the binned
curve is monotone decreasing in **33 of 42** section × module fits — the clean
published-looking gradient, reproduced on Bio's real module scores.

* Adding **only receiver cell-type intercepts** (the μ_{c_i} that Section 6.1
  already specifies): surviving fraction **0.344**.
* Adding the full N5 block on top: **0.019**.

So **66 % of the naive amplitude is the fact that the mix of cell types changes with
distance from a sender**, before any signalling. Phase 2's Figure 2a was measuring
composition. This is the single most important correction Bio's cell-type labels
bought us.

### 3.2 Of what remains, the confound is technical and geometric — not anatomical

Within receiver cell type, each N5 sub-block alone (median over 223 fits, in-band):

| N5 sub-block | contents | SF alone |
|---|---|---|
| **tech** | log counts, log genes detected, log cell area, log nucleus area | **0.288** |
| **dens** | density at 25/50/100 µm, 1-NN distance | **0.219** |
| comp | 20-NN receiver cell-type composition | 0.474 |
| anat | zonation score (+quadratic), log distance to tissue boundary | 0.810 |
| seg | `segmentation_method`, 3 levels | 0.998 |

Cumulative, added in that order: 0.288 → 0.303 → **0.085** → 0.044 → 0.056.

Two things follow, and one of them contradicts the plan.

1. **Transcript depth and local density each remove ~75 % of the effect on their
   own.** This is the Phase 1 mechanism — local cell density mechanically confounds
   distance-to-nearest-sender — appearing in real tissue, now with a technical
   partner: Bio Phase 3 §4 shows the Tier A caller is enriched 1.6× in the *lowest*
   depth quintile within cell type, so senders are systematically shallow cells and
   their neighbourhoods are systematically different for reasons that are not
   biology.
2. **Zonation is not the dominant confound.** Section 11 and Section 28 both
   predicted it would be, for liver specifically. It is not: the zonation covariate
   alone leaves SF = 0.843 overall and **1.043 in hepatocytes** — for the cell type
   the zonation argument is actually about, conditioning on zonation removes
   *nothing*. (It matters elsewhere: `Biliary/ductular` SF_zon = 0.244, where the
   score is a proxy for the ductular-reaction geography rather than for zonation.)
   We should report this as a prediction of the plan that the data refuted.

---

## 4. λ is not identified, and the regressor is a sender-calling rate

### 4.1 λ̂ rails

Of 315 primary fits, **200 (63 %) rail at a grid bound**: 114 at the 7 µm resolution
floor, 86 at the 50 µm ceiling. The interior median is 15.6 µm. Under the full
control the donor-bootstrap CI on λ̂ spans **the entire grid, [7, 50] µm, in 39 of 42
pooled cell-type × module fits** (§7). There is no length constant to report.

### 4.2 The structural reason: distance-to-nearest-sender ≈ ρ^(−1/2)

Bio Phase 3 §3.3 found that for Tier C ligands, log(median distance to the nearest
ligand+ cell) regressed on log(ligand+ density) gives r² = 0.987 with a slope of
−0.54 against a Poisson prediction of −0.50. **It generalises to the main estimand.**

For a homogeneous Poisson process of intensity ρ, the distance from an arbitrary
point to the nearest event has median exactly √(ln2/πρ) = 0.4697 ρ^(−1/2).
Measured across **77 (section × sender-definition) combinations** — 11 sections ×
7 sender calls:

| subset | n | slope | r² | observed ÷ Poisson median d |
|---|---|---|---|---|
| all sections × all sender definitions | 77 | **−0.524** | **0.984** | 1.034 |
| `cdkn1a_pos` only, across sections | 11 | −0.486 | **0.998** | 1.034 |
| Tier A percentile calls | 33 | −0.515 | **0.998** | 0.985 |
| in-band sections only | 42 | −0.528 | 0.983 | 1.033 |

**98 % of the variance in the model's independent variable, across every section and
every sender definition we have, is one number: how many cells were called
senders.** The observed/Poisson ratio of 1.03 (range 0.92–1.33) says the called
senders are placed essentially at random at the scale that sets this distance, which
the direct clustering measurement confirms: the Ripley-K ratio at 50 µm against the
N1 null is **1.109** for `tierA_p95` (1.04–1.19 across sections), 1.263 for
`cdkn1a_pos`, 1.556 for `senepy_p95` — the very bottom of the Phase 1 clustering
axis, where κ = 0 corresponded to 1.1–3.0. This is **Figure 2d**.

The practical consequence is visible in N7: moving the Tier A threshold from p90 to
p99 changes sender prevalence 8.6 % → 0.86 %, and the median λ̂ moves 7.0 → 19.7 µm.
Nothing biological changed; the calling rate did.

### 4.3 Density normalisation does not rescue it

Following the brief's step 3, I refitted every (section × call × receiver type ×
module) cell on three distance scales with the same relative window and grid: raw
microns, distance in units of the section's Poisson expectation
d_pois = 0.4697 ρ^(−1/2), and distance in units of the section's median
nearest-neighbour distance. Between-section spread of log λ̂ (sd, median over 42
celltype × module × call cells):

* raw microns **0.726** | Poisson-normalised **0.705** | median-NN-normalised **0.696**

Normalisation changes nothing, because λ̂ is not tracking density either: the r² of
log λ̂_raw on log sender density is **0.023**, slope −0.097 (a pure density readout
would give −0.500). The honest statement is therefore *sharper* than "λ is a density
readout": **the regressor is a density readout to r² = 0.98, and λ̂ is not even that —
it is a factor-of-two-in-either-direction wander over the admissible grid.** Phase 2
found the same for packing normalisation. This is the same conclusion Bio reached
for the cross-ligand ordering, arrived at independently for the main estimand, and
it belongs in the abstract.

### 4.4 Where the real sections sit on the Figure 1 regime map

Two measurements place them, both new here:

* **Sender clustering**: Ripley-K ratio at 50 µm = 1.11 (primary call) — the κ ≈ 0
  row of the Phase 1 grid.
* **Baseline autocorrelation length ℓ**: the spatial correlogram of each Tier B
  module after removing receiver cell-type means, fitted with an exponential,
  gives ℓ = **53–109 µm** across modules in the in-band sections (zonation score
  ℓ = 92 µm with ρ(50 µm) = 0.49; the Tier A sender score has ρ(10 µm) = 0.017,
  i.e. essentially none).

With λ̂ at an interior median of 15.6 µm, that is **ℓ/λ ≈ 3.5–7**, at or beyond the
right-hand edge of the Phase 1 grid, which stopped at ℓ/λ = 4 where naive bias was
+114 % at κ = 0 and iid CI coverage was **0.00**. The real sections sit in the
regime the synthetic study identified as non-identifiable, and we can now say so
with a measurement rather than an assumption.

---

## 5. Null by null

### N1 — cell-type-stratified label permutation (1,000 permutations, λ held fixed)

**SF = 0.716 [0.403, 0.882].** Weighting this most, as instructed, requires saying
what it does and does not license here.

First, it is well posed on this data in a way it was not guaranteed to be: because
senders are almost unclustered (§4.2), the permuted distance distribution closely
matches the observed one (7259: median 24.0 → 22.8 µm, p99 73.0 → 65.8 µm), so the
N1 contrast is not comparing different distance ranges.

Second — and this is the part that matters — **N1 is not protective in this
dataset.** Phase 2 established that under a *shared latent field* confounder, N1 is
the only null that discriminates (0.006 at β = 0 vs 0.69 at β = 1). But the
confounder here is not a shared latent field: it is a property of the sender cells
themselves (they are shallow, they sit in denser neighbourhoods). Permuting sender
labels destroys that property along with any signal, so a large SF is exactly what
confounding of this kind predicts. The observed 0.716 sits inside the Phase 2
"real effect" range (0.69–0.86) and tells us nothing, because the mechanism Phase 2
calibrated it against is not the mechanism operating here.

The informative version is **N1 applied to the residual after N5 + N6**:
**SF = 0.987 [0.921, 1.058]**. The small amount that survives conditioning *is*
sender-identity-specific — permuting sender labels removes essentially all of it.
That is the one positive signal in this report, and §6 bounds how much to believe it.

### N2 — matched decoys, within section and within cell type

Implemented as the **shared-λ two-kernel form** (`BlockProfiler.fit2_shared`), not
the plan's literal β_true − β_decoy, which Phase 1 showed returns −1.44 from scale
mismatch. Decoys matched on cell type, local density at 50 µm, log transcript
counts, 20-NN cell-type composition and — per Section 11 — the zonation score.

**SF = 0.943 [0.905, 0.980].** Matching quality, reported **per arm** (Bio Phase 2
§5), is excellent: match rate 1.000 in every section; max |SMD| falls from
0.412–0.459 to **0.013–0.033** in SBR and from 0.274–0.387 to **0.019–0.032** in
sham. Every matched-decoy analysis in this literature would report those numbers as
success.

**Passing SMD < 0.1 bounds nothing.** We match to |SMD| ≤ 0.033 on five covariate
families and the decoy kernel still absorbs only 6 % of the effect, while
*regressing on the same covariates* absorbs 92 %. Matching balances the covariates
between senders and decoys; it does not remove the dependence of the response on
those covariates at the receiver, which is where the confounding acts. This is the
real-tissue version of Phase 1's finding that N2 is nearly inert on λ̂, and it is a
concrete, quantified warning for a control the plan calls "the single most important
number in the paper".

### N3 / N4 — torus shift and rotation: **calibration failures, not checks**

**N3 SF = 1.000 [0.992, 1.008]. N4 SF = 0.964 [0.867, 1.042].** Zero of 160 fits
have SF ≤ 0 under either.

This is a headline, not a footnote. On synthetic tissue with **no planted kernel at
all**, in the confounded regime, Phase 2 measured the torus shift reporting a **98 %**
surviving fraction and rotation **89 %** — they certify pure confounding as a real
effect. On real tissue they now report **100 %** and **96 %** for an effect that
nuisance conditioning shows to be 92 % nuisance. Their Type I error was 0.92 and
1.00 respectively on synthetic data at β = 0; here **88 % and 90 %** of fits reject at
p < 0.05 (N1: 77 %).

The reason is structural: N3 and N4 preserve sender clustering and receiver
autocorrelation and destroy only their *alignment*. A shared cause produces
alignment without signalling, and a sender-intrinsic confounder (depth, size) is not
touched by moving the senders at all. **The torus shift tests "are the senders
aligned with the response field?", which is not "is there a SASP effect?"** If the
paper shows that existing tools fail the torus shift,
it must simultaneously report that *passing* it proves very little.

> **Correction 2026-08-27 (citation audit CIT-1).** This sentence read "reproduces the
> CellWHISPER finding that existing tools fail the torus shift." CellWHISPER contains no
> torus shift; its null permutes cell locations within each cell type, which is this
> project's N1. The torus-shift failure reported here is this project's own result. The
> null itself is Lotwick & Silverman (1982); Mrkvička et al. (2021), *Spatial Statistics*
> 42:100430, show it is liberal and that it assumes a rectangular window.

Figure 2b makes this visible: the observed binned curve, the matched-decoy curve and
the torus-shift 95 % band lie on top of each other in every panel.

### N5 — nuisance conditioning

Covariates: receiver cell type (by stratification), log transcript counts, log genes
detected, log cell area, log nucleus area, density at 25/50/100 µm, 1-NN distance,
20-NN cell-type composition, `segmentation_method` (3 levels), zonation score and its
square, log distance to tissue boundary. `dist_to_portal_triad_um` is **not** used —
Bio's SBR failure now replicates 5/5.

**SF = 0.084 [−0.130, 0.261]**, with 38 % of fits at or below zero. Attribution in
§3.2.

### N6 — receiver-baseline conditioning

For each receiver, the mean module score over its 20 nearest neighbours **excluding
senders** and excluding itself, entered as a covariate. **SF = 0.486 [0.222, 0.756].**
This control is deliberately aggressive — a genuine SASP effect on a neighbourhood
would partly be absorbed by it — so it is reported as a bracket, not as the estimate.
It is also nearly redundant with N5 (N5 alone 0.084, N5+N6 0.082).

### N7 — sender threshold and caller sensitivity

Six sender definitions across the three callers Bio validated, all on the in-band set:

| call | prevalence | λ̂ railed | median λ̂ | naive β/sd | SF N5 | SF N2+N5+N6 | SF N1 |
|---|---|---|---|---|---|---|---|
| `tierA_p90` | 8.58 % | 0.71 | 7.0 µm | 0.395 | 0.069 | 0.088 | 0.777 |
| `tierA_p95` | 4.29 % | 0.63 | 13.3 µm | 0.326 | 0.084 | 0.082 | 0.716 |
| `tierA_p99` | 0.86 % | 0.57 | 19.7 µm | 0.320 | 0.246 | 0.247 | 0.728 |
| `cdkn1a_pos` | 4.29 % | 0.73 | 9.6 µm | 0.506 | 0.385 | 0.203 | 0.911 |
| `senepy_p95` | 3.18 % | 0.77 | 18.3 µm | 0.476 | 0.302 | 0.219 | 0.895 |
| `senepy_p99` | 0.63 % | 0.64 | 18.2 µm | 0.516 | 0.444 | 0.283 | 0.886 |

The surviving fraction is **caller-dependent by a factor of 3–5** (0.082 to 0.283)
and λ̂ moves by a factor of 2.8 with the threshold. Given Bio's finding that the
three callers agree at or below chance even after matching on cell type *and* depth
decile, "λ under caller X" is a property of caller X. No caller yields a controlled
SF above 0.29.

### N8 — disjointness, scrambled response, and the DeepScence circularity

**(a) Disjointness confirmed.** Tier A ∩ Tier B = **0 genes** across all seven
modules, on-panel, in every section.

**(b) Scrambled response — the most damaging single result after N5.** Using Bio's
Tier E3 sets (500 expression-matched random gene sets per module, matched on size and
on mean log-normalised expression in 20 quantile bins), I scored 200 per module with
the same estimator scanpy's `score_genes` uses, on the same cells, at the same λ̂,
and fitted the same naive kernel. Per module, over the in-band sections:

| module | observed β/sd | random-set mean | random-set sd | **% of random sets with β ≥ observed** |
|---|---|---|---|---|
| emt_ecm | 0.419 | −0.143 | 0.157 | **0.0 %** |
| tnfa_nfkb_proximal | 0.313 | 0.052 | 0.163 | **5.0 %** |
| secondary_senescence | 0.190 | 0.005 | 0.147 | **7.0 %** |
| il6_jak_stat3 | 0.193 | 0.026 | 0.151 | **12.5 %** |
| interferon_response | 0.218 | 0.078 | 0.151 | **20.0 %** |
| downstream_arrest | 0.067 | 0.001 | 0.148 | **37.5 %** |
| oxidative_stress | −0.220 | −0.110 | 0.186 | **77.0 %** |

**A random, expression-matched gene set of the same size produces a distance
gradient as strong as the real Tier B module between 5 % and 77 % of the time.**
Only `emt_ecm` is clearly outside the random distribution. Note that the
*mean-based* surviving fraction (median 0.935) is the wrong statistic here: the
random-set distribution is roughly centred on zero, so subtracting its mean tells you
nothing — its **spread** (sd 0.15–0.19 β/sd units) is what matters, and it is the
same size as the observed effects. The table above is the number to report.

**(c) DeepScence/CoreScence circularity, quantified.** Building a CoreScence sender
call (v2, occurrence ≥ 5, human→mouse via the MGI ortholog report Bio used;
33 genes on panel) and fitting the same Tier B readouts:

* CoreScence senders give a **1.7× larger** naive amplitude than Tier A senders
  (β/sd 0.239 vs 0.143) on identical receivers — the circular caller looks stronger.
* CoreScence ∩ Tier B = **51 gene memberships**, 0–37 % of a module.
* Stripping the shared genes out of the response module changes the CoreScence-sender
  amplitude by a median ratio of **0.993** overall — but by **0.698 for
  `secondary_senescence`**, the module that shares 14 of its 38 genes. So the
  circularity is concentrated, not diffuse: **31 % of a DeepScence-sender →
  secondary-senescence fit is literally the same genes on both sides**, while the
  other six modules are largely unaffected by gene sharing (their inflation must come
  from the caller's own depth loading, which Bio measured separately).

The operational recommendation stands and is now quantified: DeepScence is a
comparison method here, never the primary caller, and `secondary_senescence` should
not be read out against a DeepScence sender call at all.

---

## 6. Is the residual real? The calibration that answers it

The residual after N2+N5+N6 is small but sender-specific (N1 removes 99 % of it). Two
calibrations bound how much to believe it, using the Phase 1 sweep where the truth is
known. Surviving fraction β̂_controlled / β̂_naive, measured on synthetic tissue:

| control | with a **planted real effect** (β_true = 1, 600 runs) | with **no effect at all** (β_true = 0, 120 runs) |
|---|---|---|
| N5 nuisance conditioning | **0.826** [0.734, 0.921] | **0.412** [0.151, 0.720] |
| N2 matched decoy (shared λ) | 0.934 [0.869, 0.988] | 0.775 [0.564, 0.924] |
| N2 + N5 | 0.832 [0.733, 0.929] | 0.357 [0.014, 0.666] |

**N5 does not destroy real effects.** At zero confounder strength it returns
β̂/β_true = 0.999; averaged over the whole confounded grid it still returns 1.333× the
true β. If the liver effect were genuine, N5 should have left ~83 % of the naive
amplitude. It left **8.4 %** — a value **below 100 % of the 600 synthetic runs with a
planted effect** (none goes that low) and at the **23rd percentile of the 120 runs
with no planted effect at all**.

Symmetrically, N2's synthetic behaviour (0.934 real, 0.775 null) explains why its
real-data value of 0.943 carries no information: N2 barely moves under either truth.

**Conclusion.** The naive amplitude in this tissue behaves like the synthetic
no-effect case, not like the synthetic real-effect case. The residual that survives
conditioning is 0.027 response-sd, below the 0.203 sd this design can detect at 80 %
power, and within the spread that expression-matched random gene sets generate.
I would not claim it.

---

## 7. T2 — the combined estimate, and why it is a case study

Pooled across the six in-band animals with section fixed effects, responses z-scored
within section, N2 + N5 + N6 applied, donor bootstrap implemented by making the
**animal** the block of the block profiler (2,000 replicates over 6 donors → 63
distinct resamples). Six receiver types × seven modules; `combined_donor.csv`.

* SF(N2+N5+N6) for hepatocytes ranges −0.090 (`interferon_response`) to +0.240
  (`tnfa_nfkb_proximal`); the donor CI includes zero for six of seven modules.
  `tnfa_nfkb_proximal` is the one exception, SF = 0.240 [0.083, 0.365].
* **The donor bootstrap CI on λ̂ under the full control spans the entire admissible
  grid, [7, 50] µm, in 39 of 42 cell-type × module cells.** There is no length
  constant, and it is not a matter of needing more animals: the profile is flat.
* Across all six receiver types, **18 of 42** controlled SFs are negative.

Per Section 24.1 this is **labelled a case study**: six animals give 63 distinct
bootstrap resamples and a lumpy interval. I have not manufactured a tighter one.

**Arm contrast within the admissible set** (2 SBR animals vs 4 sham):

| arm | animals | naive β/sd | SF N5 | SF N2+N5+N6 | SF N1 |
|---|---|---|---|---|---|
| SBR | 2 | 0.300 | 0.036 | 0.053 | 0.551 |
| sham | 4 | 0.392 | 0.113 | 0.104 | 0.759 |

The arm does not rescue the effect; if anything the admissible SBR sections are
*worse*. With two SBR animals this is a description, not a test.

---

## 8. T3 — zonation-stratified fits (Section 11)

Hepatocytes only, in-band sections, split by Bio's hepatocyte-tertile
`compartment_label`:

| stratum | fits | median n | λ̂ railed | median λ̂ | naive β/sd | SF N5 | SF zonation | SF N2+N5+N6 |
|---|---|---|---|---|---|---|---|---|
| pooled | 42 | 63,848 | 0.74 | 8.1 µm | 0.273 | 0.124 | **1.043** | 0.107 |
| periportal | 42 | 21,610 | 0.86 | 20.6 µm | 0.333 | 0.202 | 0.995 | 0.173 |
| midzonal | 42 | 21,300 | 0.79 | 7.0 µm | 0.337 | 0.227 | 0.987 | 0.197 |
| pericentral | 42 | 20,937 | 0.76 | 16.2 µm | 0.309 | 0.022 | 0.988 | −0.034 |

**The answer to Section 11's question is: the kernel does not exist pooled, and it
does not vanish within zones either — because zonation was never the confounder.**
The naive amplitude is if anything *larger* within zones (0.309–0.337 sd) than pooled
(0.273 sd), the zonation covariate alone removes ~0 % in hepatocytes (SF 0.99–1.04),
and the controlled amplitude stays in the same 0.02–0.20 range in every zone. The
Section 11 hypothesis — "any response program that varies with zonation produces a
clean decay curve that is entirely zonation" — is **testable, tested, and false for
this tissue**. That is a figure panel worth keeping precisely because it is a
negative: the anatomical confound everyone would have adjusted for is not the one
that matters, and the ones that do matter (sequencing depth, cell size, local
density, neighbourhood composition) are the ones nobody adjusts for.

Per-module amplitudes pooled vs within-zone are in `summary_phase3.txt` §7;
`oxidative_stress` is the one module with a consistently **negative** amplitude
(β/sd −0.18 to −0.26 in every zone), i.e. oxidative-stress score *rises* with
distance from a sender — the direction Section 11 predicts for a pericentral
program, and the only place the zonation story shows up at all.

---

## 9. Figures delivered

* **Figure 2b** (`figures/figure2b.png`) — binned response vs distance for
  hepatocytes in all six admissible sections, two modules, with the matched-decoy
  curve and the torus-shift 95 % band overlaid. The three curves are
  indistinguishable outside the first 5 µm bin — which contains a **median of 44
  receivers, 0.48 % of the data** (5–63 hepatocytes per section). The entire visible
  "contact spike" rests on a few dozen cells.
* **Figure 2c** (`figures/figure2c.png`) — the surviving-fraction table above as a
  dot-and-IQR plot, 160 fits per null, ordered from the nulls that remove nothing to
  the nulls that remove everything.
* **Figure 2d** (`figures/figure2d.png`, new) — median distance to nearest sender vs
  sender density, 77 section × sender-definition points, against the
  homogeneous-Poisson slope of −1/2. This is the identifiability result and I think
  it is the strongest panel in the paper.

Figure 2a from Phase 2 must be **regenerated**: it was fitted unstratified on
provisional module scores, and §3.1 shows two thirds of what it displays is cell-type
composition. The code path exists; see "not got to".

---

## 10. What this means for the paper

1. **Lead with the surviving fractions and Figure 2d.** The result is: a clean,
   monotone, highly significant distance gradient exists in mouse liver; it is 66 %
   receiver cell-type composition and 92 % nuisance; it passes the torus shift,
   rotation and matched-decoy controls essentially intact; and the model's
   independent variable is a sender-calling rate to r² = 0.98.
2. **Report N3/N4 as measured calibration failures**, with the synthetic evidence
   (98 %/89 % surviving fraction at β_true = 0) attached, not as passed checks.
3. **Report that matching to |SMD| ≤ 0.03 bounds nothing**, with the number.
4. **The negative is bounded**: controlled amplitude ≤ 0.20 response-sd at contact at
   80 % power; no λ identified anywhere on [7, 50] µm.
5. **Section 11's zonation prediction was wrong for this tissue**, and the confound
   that matters is technical. That is a more useful warning to the field than the one
   the plan anticipated, because depth and segmentation area are not covariates that
   spatial-CCC tools currently offer.
6. Section 28 pre-authorised this outcome. It is the ICBINB / ml4spatialbio paper.

---

## 11. Reproduce

```bash
cd /workspace/code
python3 -u -c "import sasp_phase3 as P; from joblib import Parallel, delayed; \
  print(Parallel(n_jobs=6)(delayed(P.prep)(s, True) for s in P.ALL_SECTIONS))"
python3 -u run_phase3_nulls.py --stage window  --sections all
python3 -u run_phase3_nulls.py --stage main    --sections all    --calls all --n-jobs 22
python3 -u run_phase3_nulls.py --stage perm    --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6
python3 -u run_phase3_nulls.py --stage perm    --sections inband \
        --calls tierA_p90,tierA_p99,cdkn1a_pos,senepy_p95,senepy_p99 --n-perm 200 --n-jobs 12
python3 -u run_phase3_nulls.py --stage curves  --sections inband --calls tierA_p95 --n-jobs 6
python3 -u run_phase3_n8.py          --sections inband
python3 -u run_phase3_strat.py
python3 -u run_phase3_attribution.py
python3 -u run_phase3_combined.py    --sections inband
python3 -u run_phase3_poisson.py
python3 -u run_phase3_lamscale.py
python3 -u summarize_phase3.py
python3 -u make_figure2bc.py
```

Everything is seeded from `MASTER_SEED = 20260820`; `OMP_NUM_THREADS=1` throughout
(Phase 1 gotcha: multithreaded BLAS costs ~10× on these many small solves).
cKDTree everywhere; no (n, n) matrix is ever formed. Total workspace footprint added:
~200 MB of CSV/npz, well inside the 20 GB quota.

### Engineering notes worth carrying forward

* `np.isin` on numpy-unicode against a pandas object Index silently falls back to an
  O(n·m) path and cost 15 minutes per section before it was spotted. Use
  `pd.Index(a).isin(b)`.
* `cmd 2>&1 | tail` in a background launcher swallows all output until the process
  exits — two debugging rounds were lost to it.
* `pkill -f <pattern>` matches the launching shell's own command line (Phase 2
  recorded this; it happened again and killed the parent).
* `BlockProfiler._acc(m, p)` slices `X[:, :p]`, so one profiler yields every *nested*
  design for free. Two extra methods (`beta_at`, `beta2_at`) were added to read β at a
  **fixed** grid index, which is what every null needs.
* Pooling sections requires the k-NN-composition and segmentation covariate blocks to
  be built from a **canonical** label list, or the design matrices have different
  widths per section.

---

## 12. What I did NOT get to

* **Figure 2a is not regenerated.** It is still the Phase 2 version: unstratified,
  provisional module scores. Given §3.1 it now overstates the gradient by ~3×. This is
  the first thing to redo and it is a 10-minute run.
* **Kernel families other than the exponential.** Every Phase 3 fit is exponential.
  Phase 2 measured a 5× spread in d̂½ across families on this data; that spread should
  be re-measured under the full control, and the AIC/spline comparison redone with
  Bio's labels. `sasp_kernels.py` supports all five families and is wired for it.
* **Superposition vs nearest-sender on real tissue.** Phase 1b showed this is the one
  cleanly identifiable estimand (15/15 in every regime) and Phase 2 recommended
  promoting it. It has not been run on the real sections. **This is the highest-value
  remaining analysis** and I would do it before anything else — it is the only place
  a positive result is still plausible.
* **λ_proximal vs λ_downstream (Section 6.4).** Not attempted; with 63 % of λ̂ railing
  the comparison is not currently well posed.
* **Held-out log-likelihood on left-out sections (Section 24.6).** Still not wired in.
  AIC only.
* **BH-FDR across programs × cell types (Section 24.5).** Deliberately not reported.
  All three permutation nulls reject at 80–93 % of fits here and at 42–100 % under a
  synthetic true null; a BH correction on p-values that carry no information would be
  worse than silence. Surviving fractions with bootstrap CIs are reported instead.
* **The Section 23 method baselines** — COMMOT, SpaTalk, CellChat v2, NCEM — on real
  and torus-shifted coordinates. Figure 4 does not exist. Given that our own N3 result
  is 1.000, this is now a high-value experiment: if those tools also report ~100 %
  signal on shifted data, Figure 4 writes itself.
* **Per-cell-type Test 3 admissibility.** The 1–20 % rule was applied at the section
  level using `Cdkn1a`+ hepatocytes, as instructed. Non-hepatocyte receivers have
  their own prevalence bands (`test3_prevalence_*.csv`) and a per-(section × cell type)
  rule would change which cells enter a few fits.
* **`deepscence_score` for the SBR sections.** Only `deepscence_sham.csv` exists, so
  the circularity analysis used a CoreScence score I rebuilt from the package's own
  `coreGS_v2.csv`. That is arguably cleaner (it isolates the gene set from the
  autoencoder), but it is not DeepScence's own output.
* **Tier C ligand-specific kernels.** Bio's Deliverable 5 verdict is that the
  cross-ligand λ ordering is uninterpretable; the within-ligand fits they say are
  still defensible (`Ccl2`→`Ccr2`, `Tnf`→`Tnfrsf1a/1b`, `Tgfb1`→`Tgfbr1/2`,
  `Il1a`→`Il1r1`) have not been run.
