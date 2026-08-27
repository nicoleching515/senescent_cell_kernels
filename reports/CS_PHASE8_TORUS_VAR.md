# CS Phase 8 — the variance-corrected random shift null (N3-var / N4-var)

**Status: complete.**

> ## ⚠ CORRECTED 2026-08-27 (record reconciliation) — two fixes applied to this file
>
> **1. λ̂ = 15.7 µm was never sourced and is withdrawn.** Every "141× the pooled λ̂ of
> 15.7 µm" and "~76 λ̂ apart" in this file was computed against a value that no file
> emits; its only provenance was this file's own "2,215 µm = 141× λ̂" (2215/141 = 15.71),
> i.e. it was circular with the claim it supported. The authoritative value is
> **λ̂ = 14.7 µm**, the pooled median of `lam_naive` over the 315 primary fits
> (in-band × `tierA_p95` × `stratum == "all"`), printed by `code/summarize_phase3.py:221`
> into `results/phase3/summary_phase3.txt` §6, `tierA_p95` row, column `medlam`;
> re-derives as 14.7321 µm from `results/phase3/main_fits.csv`. **IQR [7.0, 50.0] µm,
> 60 % of fits railed at a grid bound** — that caveat travels with every use.
> Downstream: 2,215 µm = **150×** λ̂ (was 141×); the 1,200 µm tile seams are **~81 λ̂**
> apart (was ~76); N4-var's 3,395 µm = 230× λ̂. **No conclusion changes — every
> dependent claim moves in the direction that strengthens it.** See
> `reports/RECORD_RECONCILIATION.md` §1.
>
> **2. Audit R3's column fix, never applied here, now is.** The §1 table and the §10
> framing blockquote said "23 % in the void" / "8 % in the void". Out of tissue is
> `1 − frac_in_occupancy` = **35.5 % (N3) / 19.9 % (N4)**; 22.8 % / 8.0 % is
> `1 − frac_retaining_a_neighbour` — "left with no real cell inside the 100 µm window".
> Two different columns of `results/phase3/null_destructiveness.csv`; never merge them
> in one sentence. The §10 blockquote's "1–63" admissible offsets was also the pre-C6
> file's value and is now **1–66**.

**Verdict, in one line: the field-standard fix agrees with the project's own numbers on
the real data — N3-var 0.996 against N3-tile 0.971 and a published 1.000 — but a direct
synthetic calibration study says N3-tile is the wrong variant to present as primary,
because tiling is measurably *more* liberal than the whole-section torus it replaced.
Present N3-var as primary.**

Nothing was overwritten. Every number below comes from a file this task produced, listed
in §9 with its md5 and the command that produced it. `results/phase3/sf_summary.csv` and
`summary_phase3.txt` are byte-identical to their pinned values (§9). No package was
installed. `figures/` was not written; `python3 code/check_figures_guard.py` passes (§9).

---

## 1. The comparison table

Six in-band sections, primary sender call `tierA_p95` on the frozen strict-33 Tier A
(`genesets/A_SENDER_FINAL_strict.txt`, 33 lines, md5 `e5daf2c73d704e5bdf0afe216998fa12`),
**1,000 permutations** for every row (§24.3 compliant), λ held at λ̂.
SF = (β̂_obs − mean β̂_null)/β̂_obs.

Source: `results/phase3/sf_summary_var.csv`, produced by `code/summarize_phase3_var.py`.

| variant | what the move is | in tissue? | median displacement | senders keeping a neighbour ≤100 µm | **median SF** | IQR | SF, full N5+N6+zonation design |
|---|---|---|---|---|---|---|---|
| **N3 bbox (ORIGINAL, published)** | wrap on the bounding box | no — **35.5 % out of tissue** | 2,910 µm | 0.772 | **0.999** | [0.989, 1.006] | 1.001 |
| N3 bbox (re-run in the C1 job) | " | " | 2,910 µm | 0.772 | 1.001 | [0.991, 1.008] | 0.999 |
| **N3-tile** | wrap inside solid-tissue tiles | yes | 479 µm | 1.000 | **0.971** | [0.906, 1.009] | 0.972 |
| **N3-occ** | bbox wrap, ≤5 % senders out of tissue | yes | **28 µm** | 1.000 | **0.302** | [0.000, 0.734] | 0.287 |
| N3-occ15 | ≤15 % out of tissue *(suppl.)* | mostly | 317 µm | 0.969 | 0.940 | [0.761, 1.011] | 0.917 |
| **N3-swap** | senders → random real cell positions | yes | 3,241 µm | 1.000 | **0.695** | [0.392, 0.872] | 1.003 |
| N3-snap | bbox wrap, then snap to nearest cell *(suppl.)* | yes | 2,977 µm | 1.000 | 0.993 | [0.950, 1.017] | 0.976 |
| **N3-var — VARIANCE CORRECTION** | **Euclidean shift, drop what leaves W, standardise** | **yes, by construction** | **2,215 µm** | **1.000** | **0.996** | **[0.975, 1.007]** | **0.997** |
| **N4 bbox (ORIGINAL, published)** | rotate on the bounding box | no — **19.9 % out of tissue** | 3,194 µm | 0.920 | **0.947** | [0.804, 1.039] | 0.992 |
| N4 bbox (re-run in the C1 job) | " | " | 3,194 µm | 0.920 | 0.952 | [0.796, 1.048] | 1.001 |
| **N4-tile** | rotate inside solid-tissue tiles | yes | 589 µm | 1.000 | **0.924** | [0.835, 1.049] | 0.994 |
| **N4-occ** | ≤5 % out of tissue | yes | **25 µm** | 1.000 | **0.183** | [0.000, 0.559] | 0.198 |
| N4-occ15 | ≤15 % out of tissue *(suppl.)* | mostly | 320 µm | 0.969 | 0.883 | [0.690, 0.989] | 0.893 |
| N4-swap | rotate, then snap to real cell positions | yes | 2,980 µm | 1.000 | 0.946 | [0.786, 1.035] | 0.958 |
| **N4-var — VARIANCE CORRECTION** | **rotate, drop what leaves W, standardise** | **yes, by construction** | **3,395 µm** | **1.000** | **0.985** | **[0.958, 1.003]** | **0.999** |

*n = 153 reportable fits for every whole-section row, 136 for the two tile rows (17 fall
below the 2,000-receiver floor once the fit is restricted to tiles). The destructiveness
columns are medians over the six sections.*

**Two things about this table differ from `CS_PHASE7_C1.md` §0 and neither is a
correction to it.**

1. **The reportable population is now 153, not 160.** It is derived from
   `results/phase3/main_fits.csv`, which the M1 end-to-end re-run regenerated at 06:13
   (md5 `35b4de9b9f5b0386779545ba3b4f42b4` when this table was built). Every row above —
   including all eleven C1 variants — was recomputed on that same current population, so
   the table is internally consistent; it is *not* comparable cell-by-cell with the
   160-fit table in `CS_PHASE7_C1.md`. The C1 conclusions are unchanged in substance:
   N3-tile 0.974 → **0.971**, N3-orig 1.000 → **0.999**, N3-occ 0.349 → **0.302**,
   N3-swap 0.721 → **0.695**, N4-tile 0.962 → **0.924**.
2. **The `rej` column is not comparable across rows** and is therefore not in the table
   above. For N3-var/N4-var the p-value is the RS_count-standardised Monte Carlo test the
   published method prescribes; for every other row it is the repo's uncorrected
   permutation p. The numbers are in `sf_summary_var.csv` (`reject_rate_p05`): 0.908 /
   0.889 for N3-var / N4-var, 0.824 / 0.850 for the bbox originals, 0.801 / 0.846 for the
   tile variants. Do not read the difference as a power comparison.

---

## 2. Recommendation: **present N3-var as primary**

**The variance correction does *not* disagree with N3-tile's 0.974.** N3-var's 0.996 sits
between N3-tile (0.971) and the published bounding-box value (0.999), and the three
differ by less than a sixth of N3-tile's own IQR width. Reported against interest, since
the brief anticipated the opposite: on this data the choice of N3 variant does not move
Contribution 3 at all. The window-matched version of N3-var — the observed statistic
re-fitted on the *same* retained window W ∩ (W+v), which is the convention N3-tile
already uses — is **0.995 [0.975, 1.008]**, so the agreement is not an artefact of
comparing statistics computed on different amounts of data.

The recommendation is therefore about **defensibility, not about the number**:

1. **N3-var is the null the classical literature prescribes for this window.** Lotwick &
   Silverman (1982) established the toroidal shift and its rectangular-window
   requirement; Mrkvička et al. (2021) §2.1.4 states in as many words that "the
   approaches using minus correction and variance correction can be applied in case of
   general (compact) observation windows", where the torus correction cannot. A liver
   section is a general compact window. This closes objection **O3** in
   `reports/NOVELTY_ASSESSMENT.md` §4 outright rather than answering it with a
   justification.
2. **N3-tile is the one variant the same paper predicts is *worse* than what it replaced,
   and §4 below measures that it is.** Mrkvička et al. §2.1.4, verbatim: *"While it is
   possible to extend the approach to windows which are finite unions of (aligned)
   rectangles, such a procedure would increase the amount of cracks in the
   autocorrelation structure. Subsequently, it would increase the liberality of the test
   of independence using random shifts with torus correction."* N3-tile is that
   construction — and in fact a stronger form of it, because it wraps inside each tile
   separately, so the seam count scales with the number of tiles (5–23 per section)
   rather than with the section boundary. Presenting as the corrected primary the one
   variant with a published prediction against it is an avoidable, and findable, own
   goal.
3. **N3-var carries no retention/displacement trade-off.** It is the only whole-section
   variant that is simultaneously (a) in tissue by construction, (b) clustering- and
   autocorrelation-preserving (it is a rigid Euclidean translation, with no wrap and
   therefore no seam anywhere), and (c) actually displaces the senders — 2,215 µm median,
   22× the 100 µm fitting window and **150× the pooled λ̂ of 14.7 µm** (IQR 7.0–50.0 µm;
   60 % of fits railed). N3-occ fails (c);
   N3-swap fails (b); N3-tile satisfies all three but introduces 4 seams per tile.

**Concretely, for the paper.** Report N3-var as the primary N3 and N4-var as the primary
N4. Keep N3-tile, N3-occ, N3-occ15, N3-swap and N3-snap as the sensitivity family they
already are, and keep the bounding-box row as the published-for-reference row. The
sentence that replaces `CS_PHASE7_C1.md` §6 item 1:

> Toroidal shifts require a rectangular window (Lotwick & Silverman 1982) and are known
> to be liberal because of the seam they introduce; the accepted remedy for irregular
> windows is the variance correction of Mrkvička et al. (2021), which shifts in Euclidean
> geometry, discards the part shifted outside W, and standardises the statistic because
> different shifts retain different amounts of data. We adopt it. The surviving fraction
> is 0.996 [IQR 0.975–1.007] under the variance-corrected shift and 0.985 [0.958–1.003]
> under the variance-corrected rotation, against 0.999 and 0.947 for the whole-section
> bounding-box versions that spatial-omics practice actually uses. Five further in-tissue
> variants span 0.302–0.993 and are reported as a sensitivity family.

And `CS_PHASE7_C1.md` §6 item 5's row `SF, N3 corrected (tile / occ / swap)` should
become `SF, N3 corrected (var / tile / occ15 / swap)` = `0.996 / 0.971 / 0.940 / 0.695`.

---

## 3. What was implemented, and from which source

**Source retrieved:** Mrkvička T, Dvořák J, González JA, Mateu J (2021), *Revisiting the
random shift approach for testing in spatial statistics*, **Spatial Statistics 42**:100430;
preprint arXiv:1911.00240, full text retrieved from `ar5iv.labs.arxiv.org/html/1911.00240`
on 2026-08-27. §2.1.3 (the construction), §2.1.4 (window shape), §2.2 (the random-field
case, RS_count and RS_ker), Theorem 1 (the 1/n justification) and §5 (the shift
distribution) were read in full and are quoted verbatim in the module docstring of
`code/phase3_null_var.py`. **Nothing here was improvised from the brief's summary.**

### 3.1 The construction, exactly as the paper specifies it

§2.1.3: *"For the shift vectors v₁,…,v_N denote W_i = W ∩ (W + v_i), i = 1,…,N, the
smaller window where both the information about Φ and (Ψ + v_i) is available. The first
step consists of producing the simulated values T_i = T(Φ|_{W_i}, (Ψ+v_i)|_{W_i}; W_i)…
Thus in the second step, the observed value T₀ and the simulated values T₁,…,T_N are
standardized to have zero mean and equal variance, i.e. we subtract the overall mean
T̄ = 1/(N+1) Σ_{i=0}^{N} T_i and divide by √var(T_i): **S_i = (T_i − T̄)/√var(T_i),
i = 0,…,N**."*

§2.2: *"a natural choice of the test statistic is the sample covariance cov(Φ(X), Ψ(X))…
Hence var(T_i) ≈ C/n_i where n_i is the number of sampling locations in W_i… thus setting
**var(T_i) ≈ 1/n_i** stabilizes the variance of S_i. We denote such variance correction
approach **RS_count**."*

§2.2, eqs (1)–(2), the alternative: *"var̂(T_i) = Σ_{k=0}^{N} (T_k − T̄)² · w_ik"* with
Nadaraya–Watson weights *"w_ik = K(‖v_i − v_k‖/h)/Σ_j K(‖v_i − v_j‖/h)"*, K Epanechnikov
— **RS_ker**. Implemented too, as a check on RS_count (§6).

§5/§6, the shift distribution: *"shift vectors with uniform distribution on a disk
centered in the origin and having radius 1/2"* for W = [0,1]², and *"random shift vectors
distributed uniformly on a disk with radius 250 m"* for a 1000×500 m² window — i.e. a
disk of radius **half the shorter side of the window**. Adopted verbatim: radius
1,738–3,918 µm across the six sections (`perm_nulls_var.csv:shift_radius_um`).

### 3.2 How the Phase 3 estimand maps onto the paper's random-field case

The Phase 3 statistic is the amplitude β of exp(−d/λ) with d the distance to the nearest
sender. Write Ψ(x) = exp(−d(x, S)/λ) for the sender set S. Then

    d(x, S + v) = d(x − v, S)   ⟹   (Ψ + v)(x) = Ψ(x − v),

so **translating the senders is exactly translating the field Ψ**, and this is the
paper's random-field case with Φ = the module score field at the cell centroids, Ψ as
above, and T = cov(Φ(X), Ψ(X)) — the paper's own §2.2 choice of statistic. RS_count's
1/n_i therefore applies without modification, with n_i the number of *receiver cells* in
W_i. This mapping is why the correction transfers at all, and it is stated in the module
docstring so a reviewer can check it.

**W is the tissue mask**, built with exactly the construction `phase3_null_geom.Geom`
already uses (25 µm occupancy grid, 5×5 binary closing, hole fill), so N3-var sees the
same window as N3-occ. Tissue fraction of the bounding box 0.658–0.858.

**Retention.** A receiver cell x is kept iff x ∈ W and x − v ∈ W. Every real cell is in
W by construction, so this reduces to x − v ∈ W. Ψ(x − v) is then read from the *full*
observed sender set, which is the correct reading of "(Ψ+v)|_{W_i}": at x ∈ W_i the field
value is Ψ(x − v) and x − v lies in W where Ψ is observed. Senders translated out of W
contribute to no retained cell — **that is the drop**, and it is soft: the shift vector is
never rejected, only the translated-out data is discarded. This is precisely the contrast
with N3-occ, which imposes a *hard* ≥95 % retention constraint on the *offset* and
therefore collapses to near-identity offsets.

### 3.3 N4-var is our extension, and is labelled as such

Mrkvička et al. define the construction for **translations only**. N4-var applies the
same drop-and-standardise principle to a rotation: for R_θ about the section centroid,
W_i = W ∩ R_θ(W), a cell x is retained iff x ∈ W and R_θ⁻¹x ∈ W, and Ψ_θ(x) = Ψ(R_θ⁻¹x).
RS_count is unchanged because it depends on the moved window only through n_i. This is
stated as an extension in the module docstring and should be stated as one in the paper.

### 3.4 What the runner reports, and why more than one SF

`perm_nulls_var.csv` carries, per fit and per null:

* `_sf` — (β̂_obs − mean β̂_null)/β̂_obs with β̂_obs on the **full** window, i.e. the same
  definition every other variant uses. This is the headline column.
* `_sf_wm` and `_beta_obs_wm` — the **window-matched** version: β̂_obs re-fitted on each
  draw's own W_i. The summariser reports it as a ratio of means
  (0.995 for N3-var, 0.985 for N4-var); the per-draw mean-of-ratios column in the CSV is
  unstable when a single draw's β̂_obs|W_i is near zero and is kept only for audit.
* `_cov_p_rscount`, `_beta_p_rscount` — the paper's Monte Carlo test,
  p = (1 + #{i ≥ 1 : |S_i| ≥ |S₀|})/(1 + N), on the sample-covariance and on the β scale.
* `_cov_p_rsker{10,20,30}` — RS_ker at bandwidths 0.1/0.2/0.3 of the largest drawn shift.
* `_cov_p_raw` — the same drawn statistics with **no** standardisation (drop-only).
* `_frac_cells_retained`, `_n_valid` — how much data each draw kept, and how many draws
  produced a statistic at all.

**Cross-check that the estimator is unchanged.** β̂_obs and λ̂ in `perm_nulls_var.csv`
agree with `perm_nulls_c1.csv` (scope `full`) to **max |Δβ| = 1.04 × 10⁻¹⁶ and
max |Δλ| = 0.0 over all 315 rows**. The fit population, the λ grid, the receiver
definition and the covariate design are imported from `run_phase3_nulls.py`, not
re-implemented. All 1,000 draws produced a valid statistic for every one of the 153
reportable fits (`_n_valid` min = 1000); the 50-cell floor never bound.

---

## 4. **The coordinator's question, answered with a measurement: tiling is more liberal**

The audit's reading of §2.1.4 is correct — the verbatim text is quoted in §2 item 2 above,
and I read it in the retrieved full text rather than taking it second-hand. But a
prediction about liberality is a statement about **type I error**, and type I error cannot
be read off the real-data table: SF is not a rejection rate, and the observed data is not
under the null. So I measured it directly.

`code/phase3_var_sim.py` runs a reduced version of the paper's own §5 design — W = [0,1]²,
100 sampling points, two **independent** stationary Gaussian fields with isotropic
exponential correlation exp(−r/s), 199 shifts uniform on a disk of radius ½, 400
replications — and adds two things the paper does not have: a second window that is
**irregular** (three overlapping discs minus a bite, 74 % of the bounding box, inside the
0.658–0.858 range the six real sections show), and a **tiled** torus arm that is N3-tile's
construction exactly (one shared offset, wrapped inside each cell of a k×k partition,
restricted to the sub-squares entirely inside the window). Fields are synthesised on a
3×-larger domain and cropped, so the realisation is **not** periodic — otherwise the torus
correction would be valid by construction and the whole comparison would be vacuous.

Rejection rate at a nominal 5 %, under the null of independence
(`results/phase3/var_sim_calibration.csv`; Monte Carlo SE ≈ 0.011 at 400 replications):

| window | corr. scale s | torus (whole window) | **torus in 4×4 tiles** | torus in 8×8 tiles | **RS_count** | RS_ker | drop, no standardisation |
|---|---|---|---|---|---|---|---|
| rectangle (paper's design) | 0.02 | 0.035 | 0.040 | 0.035 | **0.033** | 0.038 | 0.005 |
| rectangle | 0.05 | 0.055 | 0.065 | 0.080 | **0.055** | 0.050 | 0.013 |
| rectangle | 0.15 | 0.048 | 0.063 | 0.063 | **0.035** | 0.013 | 0.013 |
| rectangle | 0.30 | **0.078** | **0.105** | 0.055 | **0.060** | 0.035 | 0.018 |
| irregular blob | 0.02 | 0.048 | 0.048 | 0.040 | **0.040** | 0.038 | 0.003 |
| irregular blob | 0.05 | 0.033 | **0.080** | 0.073 | **0.043** | 0.033 | 0.003 |
| irregular blob | 0.15 | 0.040 | **0.083** | 0.050 | **0.053** | 0.033 | 0.008 |
| irregular blob | 0.30 | **0.073** | **0.118** | 0.085 | **0.055** | 0.040 | 0.020 |

**Four readings, and the first is the one that matters.**

1. **Mrkvička §2.1.4's prediction holds, measurably.** Tiling is more liberal than the
   whole-window torus in **7 of the 8 cells**, and on the irregular window the 4×4 tiling
   is at **0.048–0.118** across all four correlation scales — **0.080–0.118 for s ≥ 0.05**,
   with 0.048 at s = 0.02 — against a nominal 0.05, while the whole-window torus is at
   0.033–0.073. At the strongest autocorrelation the tiled test rejects a true null
   **2.35×** as often as it should (0.1175 / 0.05 exactly). *(Corrected 2026-08-27, record
   reconciliation: "0.080–0.118 on the irregular window" omitted the s = 0.02 cell at 0.048
   — quote the basis with the range; and "2.4×" was a rounding up of 2.35.)* **The project's own correction C1 replaced a liberal
   test with a more liberal one.** That is the finding this section exists to record, and
   it is against the project's interest.
2. **RS_count holds the nominal level throughout** — 0.033–0.060 across every window and
   every correlation scale, ~~never outside ±0.010 of nominal~~ **within 1.6 Monte-Carlo SE of
   nominal everywhere**. This reproduces the paper's
   headline claim and is the check that our `rs_count` is the paper's estimator and not
   something adjacent to it.
   ***[Corrected 2026-08-27 (remediation pass) — the "±0.010" gloss **contradicted the range in
   its own sentence** and is struck; this is the site the same wording in `WRITING_PACK.md` §5.7
   was taken from. Nominal is 0.05, so ±0.010 is [0.040, 0.060], but
   `results/phase3/var_sim_calibration.csv`, column `rs_count_reject_05`, all eight cells read
   **0.0325, 0.0550, 0.0350, 0.0600, 0.0400, 0.0425, 0.0525, 0.0550** — two are below 0.040
   (rectangle s = 0.02 at **0.0325**; rectangle s = 0.15 at **0.0350**) and the largest shortfall
   is **0.0175**. **The range 0.033–0.060 is correct and is confirmed; only the gloss was false.**
   The Monte-Carlo SE at 400 replicates is √(0.05·0.95/400) = **0.0109**, so 0.0325 is **1.6 SE**
   low — statistically unremarkable, which is why the replacement states the SE bound. The item's
   conclusion (that our `rs_count` is the paper's estimator) is unaffected.
   `AUDIT_NUMBERS_FINAL.md` R2.]***
3. **The torus correction is liberal at strong autocorrelation** (0.073–0.078 at s = 0.30),
   reproducing the paper's own diagnosis.
4. **The standardisation is doing the work, not the dropping.** Drop-and-*don't*-
   standardise rejects at 0.003–0.020 — badly conservative, because the small-|W_i| draws
   dominate the tails of the null distribution. This is exactly the failure §2.1.3 says
   the standardisation exists to fix, and it means a naive "just discard what leaves the
   tissue" null is not a substitute for the published method.

**The honest caveat.** This is a synthetic study on Gaussian fields with a bounded
correlation range, at 100 sampling points; the real fits have 2,000–31,000 receivers and a
non-Gaussian response. It establishes the *direction* of the tiling effect and validates
the implementation. It does not license a numerical claim about the type I error of the
Phase 3 fits — the project's instrument for that is the A7 negative-control-probe family
(9–16 % against 5 % nominal), which is a different measurement of a different thing.

**And the empirical result on the real data does not follow the prediction.** N3-tile's
SF (0.971) is *below* the bounding-box value (0.999), and its raw rejection rate (**0.801**,
per §2 of this report; the 0.802 written here was a typo and the two sections disagreed)
is *below* the bounding-box one (0.824) — i.e. on this data the tiled null looks slightly
more conservative, not more liberal. Two reasons that is not a contradiction. First, SF
and a real-data rejection rate are not type I error. Second, liberality from seams scales
with the fraction of cells lying within one correlation length of a seam: at a 1,200 µm
tile side and a pooled λ̂ of 14.7 µm the seams are **~81 λ̂** apart, so the damage per seam is
real but the affected fraction of cells is small. The prediction is directionally
confirmed in the regime where it can be measured and is quantitatively small in ours. Say
both; a reviewer who knows §2.1.4 needs to see that we knew it too.

---

## 5. The destructiveness diagnostic, same columns as every other variant

`code/run_phase3_var.py --stage diag` → `results/phase3/null_destructiveness_var.csv`
(20 draws per section per null; the same construction as `code/phase3_null_diag.py`).

| section | null | shift-disk radius | senders retained in W | receiver cells retained in W_i | **retaining a neighbour ≤100 µm** | real → null median neighbours | **median displacement** |
|---|---|---|---|---|---|---|---|
| 7259 | N3-var | 1,738 µm | 0.600 | 0.614 | **1.000** | 259.0 → 247.9 | 1,222 µm |
| 7260 | N3-var | 3,612 µm | 0.549 | 0.538 | 0.995 | 169.0 → 157.4 | 2,599 µm |
| 7001 | N3-var | 3,600 µm | 0.639 | 0.626 | **1.000** | 105.0 → 97.3 | 2,149 µm |
| 7248 | N3-var | 3,918 µm | 0.579 | 0.592 | 0.996 | 138.0 → 128.5 | 2,605 µm |
| 7352 | N3-var | 2,922 µm | 0.624 | 0.618 | **1.000** | 142.0 → 136.8 | 1,814 µm |
| 7435 | N3-var | 3,229 µm | 0.585 | 0.591 | **1.000** | 130.0 → 130.3 | 2,281 µm |
| **median, N3-var** | | | **0.592** | **0.603** | **1.000** | **140.0 → 133.5** | **2,215 µm** |
| **median, N4-var** | | | **0.841** | **0.836** | **1.000** | **140.0 → 135.7** | **3,395 µm** |

Read against the C1 table in `CS_PHASE7_C1.md` §0:

* **Retention is 1.000**, as for N3-tile / N3-swap / N3-snap, and against **0.772** for the
  published bounding-box shift. The void problem is solved.
* **Neighbour thinning is 140.0 → 133.5 (−4.6 %)**, against 140.0 → 119.7 (−14.5 %) for the
  bounding-box shift and 143.0 → 139.2 (−2.6 %) for N3-tile. The residual thinning is
  because the retained shifted senders sit preferentially in W ∩ (W+v), whose boundary is
  denser in low-density edge tissue than the section interior.
* **Displacement is 2,215 µm**, against 28 µm for N3-occ — 22× the 100 µm fitting window
  and **150× the pooled λ̂ of 14.7 µm** (IQR 7.0–50.0 µm). N3-var does not have N3-occ's
  degeneracy, and the reason is structural, not a matter of tuning: the constraint is on
  the *data retained*, not on the *offset admitted*.
* **The cost is 41 % of the receiver cells** (median 0.603 retained under N3-var, 0.836
  under N4-var, which retains more because a rotation about the centroid overlaps itself
  more than a translation does). That cost is what the variance standardisation exists to
  absorb, and §6 shows it is absorbed.

For completeness, the column the other variants report over *all* shifted senders — before
the drop — is `frac_retaining_a_neighbour_all_shifted`: 0.583–0.665 for N3-var (median 0.640). That is the
number a torus-free Euclidean shift *would* have if you kept the translated-out senders,
and it is why keeping them is not an option.

---

## 6. Does RS_count's var(T_i) ≈ 1/n_i hold on our data?

The 1/n_i factor is justified by the paper's Theorem 1 under assumptions our data does not
literally satisfy (two independent *stationary Gaussian-ish* fields with compactly
supported autocovariance). `code/phase3_var_validate.py` tests it directly: 2 sections,
4 reportable fits each, 400 draws per null, storing per-draw (T_i, n_i), then regressing
log sd(T | n-bin) on log n. RS_count predicts a slope of **−0.5**.

`results/phase3/var_variance_check.csv`, restricted to the 12 of 16 cases where the draws
give at least a 2× dynamic range in n_i (below that the slope is not identified):

| null | median slope | range | median Spearman ρ(n_i, |T_i − T̄|) |
|---|---|---|---|
| N3-var | **−0.451** | −0.585 … −0.122 | −0.122 |
| N4-var | **−0.492** | −0.644 … −0.119 | −0.166 |

**The assumption holds well enough to use.** Reported against interest: for four of the
sixteen cases — all N4-var on section 7260, where the rotation retains 22,387–27,567 cells,
a range of only 1.23× — the slope is unidentified and comes back positive (+1.05 to +2.01).
Those are excluded above and they are excluded because the *design* cannot resolve the
slope there, not because they disagree.

**RS_ker agrees with RS_count on the real data**, which is the paper's own recommended
consistency check. Rejection rates over the 153 reportable fits
(`results/phase3/var_pvalues.csv`):

| null | RS_count (2-sided) | RS_count (1-sided) | RS_count on β | RS_ker h=0.1 | h=0.2 | h=0.3 | drop, no standardisation | repo's naive p |
|---|---|---|---|---|---|---|---|---|
| N3-var | 0.908 | 0.961 | 0.935 | 0.712 | 0.824 | 0.863 | 0.732 | 0.817 |
| N4-var | 0.889 | 0.941 | 0.922 | 0.725 | 0.771 | 0.765 | 0.850 | 0.889 |

RS_ker is uniformly lower than RS_count and rises toward it as the bandwidth widens, which
is the expected direction: at a narrow bandwidth the Nadaraya–Watson variance is estimated
from few neighbouring shifts and is noisy, which inflates √var̂ for the extreme draws. The
paper reaches the same conclusion — §7: *"The best results were achieved with the
asymptotic order of variance 1/n of the sample covariance, i.e. the RS_count method. Thus
this method can be recommended… both for its simplicity and its performance."* We report
RS_count as primary for the same reason.

---

## 7. Per-section results

`results/phase3/perm_nulls_var.csv`, median over each section's reportable fits.

| section | fits | λ̂ median | shift radius | tissue frac. | **N3-var SF** | **N4-var SF** | N3-var cells retained | RS_count reject rate |
|---|---|---|---|---|---|---|---|---|
| 7001 sham 52 wk | 40 | 22.4 µm | 3,600 µm | 0.735 | 0.996 | 0.983 | 0.591 | 0.950 |
| 7248 sham 26 wk | 24 | 14.6 µm | 3,918 µm | 0.858 | 0.994 | 0.984 | 0.563 | 0.875 |
| 7259 sbr 26 wk | 25 | 12.8 µm | 1,738 µm | 0.673 | 1.004 | 0.983 | 0.628 | 0.920 |
| 7260 sbr 26 wk | 31 | 50.0 µm | 3,612 µm | 0.705 | 0.996 | 0.956 | 0.552 | 0.903 |
| 7352 sham 2 wk | 11 | 10.0 µm | 2,922 µm | 0.658 | 0.965 | 0.982 | 0.590 | 0.909 |
| 7435 sham 10 wk | 22 | 7.2 µm | 3,229 µm | 0.695 | 0.991 | 1.009 | 0.638 | 0.864 |

**No section is degenerate.** Compare `CS_PHASE7_C1.md` §3, where section 7001's N3-occ
admitted only the identity offset and returned SF = −0.000 by construction: 7001 under
N3-var returns 0.996 on 40 fits. The per-section spread for N3-var is 0.965–1.004, tighter
than every other whole-section variant in the table.

---

## 8. What I did **not** do

* **Only the primary sender call.** `perm_nulls_var.csv` covers `tierA_p95` on the six
  in-band sections. The N7 sender-definition axis (five further calls) has **not** been run
  against N3-var/N4-var. The C1 closeout showed the correction is invariant to the sender
  definition across a 0.5–9.0 % prevalence range; that is evidence, not proof, that
  N3-var would behave the same.
* **No sensitivity to the shift-disk radius.** The paper's half-the-shorter-side convention
  was adopted without varying it. A smaller radius retains more data and displaces less; at
  the limit it approaches N3-occ's failure mode. The current radius displaces **150× λ̂**, so
  the null is comfortably decoupled, but the sensitivity is not measured.
* **The covariate-adjusted family is at 200 permutations, not 1,000.** The `*_full_sf`
  column for N3-var/N4-var (0.997 / 0.999) comes from `perm_nulls_var_full200.csv`, which
  refits the full N5+N6+zonation design on each draw's own retained window — a QR per
  (fit, draw) that costs ~5× the intercept-only loop. The headline `_sf` family is at the
  full 1,000. Every other row's `full_sf` is at 1,000.
* **`RS_var`** — the paper's third variance estimator, which plugs fitted variograms into
  the exact variance formula — is not implemented. The paper itself calls it "too complex
  and time consuming in order to be generally applied" (§7) and recommends RS_count.
* **The point-process branch** (§2.3: cross-K, cross-nearest-neighbour distance,
  RS_K/RS_G) is not implemented and does not apply: our statistic is a random-field sample
  covariance, not a point-process summary function. The paper is explicit that the variance
  correction works *well* in the random-field case and *badly* in the point-process case
  (§7), so the mapping in §3.2 is load-bearing and should be stated in the Methods.
* **No figure was made or changed.** `figures/` was not written. Figure 2's null battery
  panels would need a new `perm_curves_var.csv`; the runner does not accumulate binned
  curves. That is a figure job for 8.7, and every number it needs is in `sf_summary_var.csv`.
* **No commit, no push, no tag.**

---

## 9. Files, hashes, and how to reproduce

**New code (nothing existing was modified):**

| File | What |
|---|---|
| `code/phase3_null_var.py` | `VarGeom` (tissue window W, shift/rotation draws, pull-backs), `rs_count`, `rs_ker`, `mc_pvalue`, `sample_cov`. Carries the verbatim Mrkvička quotes. |
| `code/run_phase3_var.py` | `--stage perm` (the null) and `--stage diag` (destructiveness). Imports `SectionFit`, `_designs`, `_expand`, `WINDOW_UM`, `MIN_RECEIVERS`, `PRIMARY_CALL` from `run_phase3_nulls.py`. |
| `code/summarize_phase3_var.py` | The comparison table. Imports `reportable` and `_stats` from `summarize_phase3_c1.py`. |
| `code/phase3_var_validate.py` | The var(T_i) ≈ 1/n_i check (§6). |
| `code/phase3_var_sim.py` | The synthetic calibration study (§4). |

**New results (all new filenames; nothing overwritten):**

```
32397d5d6b626bb59924fdbc04669dac  results/phase3/perm_nulls_var.csv          (1,000 perms, headline)
afc30d27028e444bd7f309d36451bb99  results/phase3/perm_nulls_var_full200.csv  (200 perms, covariate-adjusted)
e788f0d824264bb1ff8a2e499809c772  results/phase3/perm_draws_var.csv          (per-draw geometry)
0d2e556082bdd9a7d8bf767cfd27fbd8  results/phase3/perm_draws_var_full200.csv
0dd0d3f75cf9f7b3d5d6f3cedb37d321  results/phase3/null_destructiveness_var.csv
3d53f7586c53acdfd80b55f2364daad5  results/phase3/sf_summary_var.csv
cc3c19b9e1149e75ae82af35a9fbcc48  results/phase3/summary_phase3_var.txt
f3e4fdc8623259e8d42de91d1c88c6b8  results/phase3/var_pvalues.csv
983301c7863888fbfca5a41baba77efe  results/phase3/var_sim_calibration.csv
19c2a600d284a91db524ae03fa5f5be0  results/phase3/var_variance_check.csv
```

**Inputs, hashed at the moment the table was built** (the M1 re-run writes these; if they
move, `summarize_phase3_var.py` must be re-run):

```
35b4de9b9f5b0386779545ba3b4f42b4  results/phase3/main_fits.csv
d906394958dbe1b99981756290c511fa  results/phase3/perm_nulls.csv
0318737e85a9b98e1b3bdd1461880ef5  results/phase3/perm_nulls_c1.csv
423d9ae3c345b57d4312c5f85da97558  results/phase3/null_destructiveness.csv
e5daf2c73d704e5bdf0afe216998fa12  genesets/A_SENDER_FINAL_strict.txt   (33 genes)
```

**Pinned files, re-verified after everything above** (unchanged from
`CS_PHASE8_C1_CLOSEOUT.md`):

```
69e3a1d3f60060deddcceba9896a7d31  results/phase3/sf_summary.csv
ecf86b9ca5460f31290e2f4c9e822ea2  results/phase3/summary_phase3.txt
```

`python3 code/check_figures_guard.py` → `OK: all 27 committed figures match (PDF date
stamps ignored)`, exit 0, run after all of the above.

**Reproduce:**

```bash
cd /workspace
python3 -u code/run_phase3_var.py --stage diag --n-rep 20 --n-jobs 6            # ~40 s
python3 -u code/run_phase3_var.py --stage perm --n-perm 1000 --n-jobs 6 \
        --no-full                                                              # 604-1003 s/section, 6 in parallel (16.8 min)
python3 -u code/run_phase3_var.py --stage perm --n-perm 200 --n-jobs 6 \
        --tag _full200                                                         # 28.5 min, covariate-adjusted
python3 -u code/phase3_var_validate.py                                         # ~3 min
python3 -c "import sys;sys.path.insert(0,'code');import phase3_var_sim as S;S.main(n_rep=400,n_shift=199,n_jobs=8)"
python3 -u code/summarize_phase3_var.py
python3 code/check_figures_guard.py
```

`MASTER_SEED = 20260820`; `--stage perm` derives its per-section seeds through
`run_phase3_nulls._expand(sections, calls, MASTER_SEED, 5000, 17)`, the same call
`--stage perm_c1` uses, so N3-var and the C1 variants sit on the same seed grid. System
`python3` (numpy 2.4.6, scipy 1.17.1, pandas 2.3.3); the `envs/sasp311` venv was not used.
No package was installed.

---

## 10. Citations to add, and the framing to use

`references.bib` has no spatial-statistics methods paper (`NOVELTY_ASSESSMENT.md` §U7).
Three entries make this section defensible, and they are the *only* thing standing between
the current text and objection O3:

* **Lotwick HW, Silverman BW (1982).** *Methods for analysing spatial processes of several
  types of points.* JRSS-B 44(3):406–413. doi:10.1111/j.2517-6161.1982.tb01221.x.
  The origin of the toroidal shift test and of the rectangular-window requirement.
* **Mrkvička T, Dvořák J, González JA, Mateu J (2021).** *Revisiting the random shift
  approach for testing in spatial statistics.* Spatial Statistics 42:100430;
  arXiv:1911.00240. The liberality diagnosis, §2.1.4's union-of-rectangles warning, and the
  variance correction implemented here.
* **`spatstat.random::rshift.ppp`** — `edge="torus"` requires a rectangular window. Worth a
  parenthetical because it is what a reviewer will check first.

**Framing, per the coordinator's instruction and `NOVELTY_ASSESSMENT.md` §2.2.** The torus
finding is **not** "we discovered that torus shifts break on non-convex tissue" — that has
been known since 1982. It is:

> A 40-year-old, documented limitation of the toroidal shift is being violated in current
> spatial-omics practice. We quantify what the violation costs on real tissue (**35.5 %
> of shifted senders land out of tissue**, and 22.8 % are left with no real cell inside the
> 100 µm window — two different columns, never merged; a 14.5 % thinning of the shifted
> sender's receiver neighbourhood; at a ≤5 % out-of-tissue tolerance only **1–66** of
> 38,080–108,375 candidate offsets are admissible and one of six sections admits only the
> identity), we give an
> exact FFT enumeration of the admissible offset set, and we adopt the remedy the classical
> literature prescribes for irregular windows — the variance correction of Mrkvička et al.
> (2021). Under it the surviving fraction is 0.996, so the conclusion is unchanged.

**And fix the null's attribution.** `SASP_Kernel_Master_Plan.md` §32 item 3 credits
CellWHISPER with the torus-shift null. CellWHISPER's null is a within-cell-type location
permutation — i.e. this project's **N1**, not N3 (`NOVELTY_ASSESSMENT.md` §U3, verified
against the preprint there). The torus shift is Lotwick & Silverman (1982). That citation
is load-bearing because Figure 4 is built on N3. It is outside this task's scope to edit
the plan, but it must be fixed before submission.

---

## 11. Summary for the roadmap status board

| Item | Status | Note |
|---|---|---|
| N3-var / N4-var implemented per Mrkvička et al. (2021) §2.1.3+§2.2 (RS_count), full text retrieved | DONE | `code/phase3_null_var.py`; N4-var flagged as our extension |
| Destructiveness diagnostic, same columns as C1 | DONE | Retention **1.000**, 140.0 → 133.5 neighbours, displacement **2,215 µm** |
| Surviving fraction, 1,000 perms, strict-33 Tier A | DONE | **N3-var 0.996 [0.975, 1.007]; N4-var 0.985 [0.958, 1.003]** on 153 fits |
| Comparison table, 7 variants + var | DONE | `results/phase3/sf_summary_var.csv` |
| Primary-variant recommendation | **N3-var / N4-var** | Not because the number moved — it did not — but because §2.1.4 predicts and §4 measures that N3-tile is *more* liberal than the null it replaced |
| Does the variance correction disagree with N3-tile's 0.974? | **NO** | 0.996 vs 0.971 on the current 153-fit population; window-matched 0.995. Contribution 3 stands as written |
| RS_count's 1/n_i assumption, tested on our data | DONE | slope −0.451 (N3-var) / −0.492 (N4-var) against a predicted −0.5 |
| Implementation validated against the paper's own simulation | DONE | RS_count 0.033–0.060 at nominal 0.05; torus 0.073–0.078 at strong autocorrelation; **tiled torus 0.080–0.118 on an irregular window** |
| N7 axis under N3-var | NOT DONE | see §8 |
| Shift-radius sensitivity | NOT DONE | see §8 |
