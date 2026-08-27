# CS Phase 7 / C1 — The torus shift, confined to tissue

**Status: complete. Verdict: the Phase 3 N3 result survives. Contribution 3 stands
as written.**

Phase 7 Part I, Step 1. Addresses the correction `Phase7_Minimal_Human_Replication.md`
§1–§3 opens with: `run_phase3_nulls.py` wrapped the N3 torus shift and the N4 rotation
on the whole-section **bounding box**, which on a non-convex liver section throws
shifted senders into empty space. `CS_PHASE4.md` §2.4 identified this and Phase 4
fixed it on solid-tissue tiles; **Phase 3 was never re-run.** It has now been re-run.

| Deliverable | Location |
|---|---|
| In-tissue null geometry (all variants) | `/workspace/code/phase3_null_geom.py` |
| Destructiveness diagnostic | `/workspace/code/phase3_null_diag.py` → `results/phase3/null_destructiveness.csv` |
| Corrected N3/N4 runner | `/workspace/code/run_phase3_nulls.py --stage perm_c1` → `results/phase3/perm_nulls_c1.csv` |
| Old-vs-new surviving fractions | `/workspace/code/summarize_phase3_c1.py` → `results/phase3/sf_summary_c1.csv`, `summary_phase3_c1.txt` |
| Figure 4 de-hardcoded | `/workspace/code/make_figure4.py` (`load_ours()`) |

**Nothing was overwritten.** `results/phase3/perm_nulls.csv`, `sf_summary.csv` and
`summary_phase3.txt` are byte-for-byte as published. Every corrected number lives in a
new `*_c1` file, and the published bounding-box null was **re-run in the same job**, at
the same 1,000 permutations, so old and new sit on the same page.

---

> ### Correction box — audit item R9, added 2026-08-27
>
> **The destructiveness numbers in this report are keyed to the PRE-C6 sender
> definition and no longer match `results/phase3/null_destructiveness.csv`.**
> That file was regenerated at 05:58:41 (`logs/m1_nulldiag.log`) during the M1
> end-to-end re-run (task 8.7), under the promoted C6 gene sets, where the
> strict Tier A went 25 → 33 genes and `tierA_p95` therefore selects different
> cells.
>
> **Section-geometry columns are unchanged** — occupied fraction of bbox, tissue
> fraction, solid-tile counts and coverage are properties of the section, not of
> the sender set. **Every sender-dependent number moved, all of them slightly:**
> N3-orig median displacement 2,974 → **2,910 µm**; N3-tile 489 → **479**;
> N4-orig retention 0.917 → **0.920**; N3-occ15 0.966 / 304 µm → **0.969 / 317 µm**;
> N4-swap 3,038 → **2,980 µm**; admissible N3 offsets for 7259 63 → **66** of
> 38,080 and for 7352 9 → **10** of 74,178; admissible N4 angles for 7259
> 12 → **13** of 720; median real neighbours within 100 µm 140.5 → **140.0**.
>
> **No conclusion in this report changes**, and the corrected surviving
> fractions on the frozen sets are in `reports/CORRECTIONS.md` §8 and §10.
> `results/phase3/figure2e_data.csv` and `figures/figure2e.png` were regenerated
> together in the 8.7 figure pass and carry the current values.

---

## 0. Headline

Six in-band sections, primary sender call `tierA_p95`, the same 160 reportable fits
(positive naive amplitude whose spatial block bootstrap CI excludes zero), 1,000
permutations, λ held at λ̂. Surviving fraction = (β̂_obs − mean β̂_null) / β̂_obs.

| variant | what it does | senders keeping a neighbour ≤100 µm | median displacement | **median SF** | IQR |
|---|---|---|---|---|---|
| **N3 (ORIGINAL, published)** | wrap on the bounding box | **0.772** | 2,974 µm | **1.000** | [0.992, 1.008] |
| N3 (ORIGINAL, re-run here) | " | 0.772 | 2,974 µm | 1.002 | [0.994, 1.010] |
| **N3-tile** | wrap inside solid-tissue tiles | **1.000** | 489 µm | **0.974** | [0.913, 1.008] |
| **N3-occ** | ≤5 % of senders out of tissue | **1.000** | **27 µm** | **0.349** | [0.073, 0.764] |
| N3-occ15 | ≤15 % out of tissue *(suppl.)* | 0.966 | 304 µm | 0.951 | [0.834, 1.027] |
| **N3-swap** | senders → random real cell positions | **1.000** | 3,252 µm | **0.721** | [0.406, 0.879] |
| N3-snap | bbox wrap, then snap to nearest real cell *(suppl.)* | 1.000 | 2,981 µm | 0.988 | [0.943, 1.026] |
| **N4 (ORIGINAL, published)** | rotate on the bounding box | **0.917** | 3,232 µm | **0.964** | [0.867, 1.042] |
| N4 (ORIGINAL, re-run here) | " | 0.917 | 3,232 µm | 0.960 | [0.870, 1.038] |
| **N4-tile** | rotate inside solid-tissue tiles | **1.000** | 596 µm | **0.962** | [0.866, 1.028] |
| **N4-occ** | ≤5 % out of tissue | **1.000** | **25 µm** | **0.273** | [0.000, 0.592] |
| N4-occ15 | ≤15 % out of tissue *(suppl.)* | 0.969 | 308 µm | 0.896 | [0.755, 0.998] |
| **N4-swap** | rotate, then snap to real cell positions | **1.000** | 3,038 µm | **0.969** | [0.862, 1.030] |

*(`sf_summary_c1.csv`; destructiveness columns are medians over the six sections from
`null_destructiveness.csv`. Tile rows are 144 of the 160 fits — 16 fall below the
2,000-receiver floor once the fit is restricted to tiles.)*

**Which §3 outcome occurred: the first.** Every variant that is *actually a torus shift
or a rotation*, confined to tissue, and that *actually moves the senders further than
the 100 µm fitting window*, returns a surviving fraction indistinguishable from the
published one:

> **N3-tile 0.974, N3-snap 0.988, N3-occ15 0.951; N4-tile 0.962, N4-swap 0.969,
> N4-occ15 0.896 — against a published N3 1.000 and N4 0.964.**

The two variants where the surviving fraction collapses do so for reasons that have
nothing to do with the void, and §3 below shows exactly what those reasons are. The
structural claim in `CS_PHASE3.md` — that N3 and N4 preserve sender clustering and
receiver autocorrelation and destroy only their *alignment*, so a shared cause survives
them intact — **is vindicated and is now stronger than it was**, because the obvious
competing explanation has been measured and eliminated.

**The bug was real, and it was worse than the appendix said.** `CS_PHASE4.md` §2.4
estimated "~20 % of shifted cells into the void". On the Phase 3 sections it is
**35.5 %** for N3 and **19.9 %** for N4 (`null_destructiveness.csv`,
`1 − frac_in_occupancy`, median over the six sections).

> **Corrected 2026-08-27 (audit item R3).** This paragraph originally read
> "**23 %** for N3 … and **8 %** for N4", quoting `1 − frac_retaining_a_neighbour`.
> That column is **not** the out-of-tissue fraction: it is the fraction of shifted
> senders left with **no real cell inside the 100 µm fitting window**, which is
> 22.8 % for N3 and 8.0 % for N4. The out-of-tissue fraction is
> `1 − frac_in_occupancy`, **35.5 % and 19.9 %**. Both quantities are real and
> both are in the same file; only the label was wrong. The bug is therefore
> **worse** than this report first claimed, so the correction runs in the safe
> direction. Never state both numbers as one quantity in one sentence.
It simply did not matter here.

---

## 1. What was wrong, measured

`run_phase3_nulls.py` used

```python
lo = sf.coords.min(0); hi = sf.coords.max(0)
span = hi - lo
return lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span
```

A liver section fills only part of its bounding box. Measured on the 25 µm occupancy
grid (`null_destructiveness.csv`):

| section | cells | occupied fraction of bbox | tissue fraction of bbox | N3: senders keeping a neighbour | N4: same |
|---|---|---|---|---|---|
| 7259 sbr 26 wk | 127,386 | 0.656 | 0.673 | 0.757 | 0.835 |
| 7260 sbr 26 wk | 202,016 | 0.656 | 0.705 | 0.790 | 0.917 |
| 7001 sham 52 wk | 165,961 | 0.644 | 0.735 | 0.787 | 0.918 |
| 7248 sham 26 wk | 224,921 | 0.784 | 0.858 | 0.916 | 0.944 |
| 7352 sham 2 wk | 139,378 | 0.620 | 0.658 | 0.726 | 0.863 |
| 7435 sham 10 wk | 172,218 | 0.662 | 0.695 | 0.740 | 0.936 |

"Tissue fraction" is the 25 µm occupancy grid after a 5×5 binary closing and hole fill,
so an empty 25 µm cell that is a sinusoid or a vessel lumen counts as tissue and only
the space outside the section counts as void.

Under the published N3, the **median number of real cells within 100 µm of a sender
falls from 140.5 to 119.7** — a 15 % thinning of the shifted sender's receiver
neighbourhood that has nothing to do with biology. Under N4 it falls from 140.5 to
129.3. Phase 4's tiles, for comparison, gave 100 % retention and 150.8 → 149.4
(re-derived from `results/phase4/null_destructiveness.csv`).

---

## 2. The variants

All in `code/phase3_null_geom.py`. Each returns shifted sender coordinates; nothing else
in the estimator changes, and λ is held at λ̂ exactly as before.

**N3-tile / N4-tile — Phase 4 solid-tissue tiles.** Tiles are chosen by
`phase4_tiles.tiles_for` itself — the Phase 4 routine, reused verbatim through a
three-attribute shim, side 1,200 µm, sliding grid at 600 µm, greedy densest
non-overlapping — and kept only when ≥98 % of the tile is tissue. That yields **5–23
solid tiles per section covering 46–69 % of cells** (`null_destructiveness.csv`). One
shared random offset (N3) or angle (N4) per permutation, wrapped inside each tile about
that tile's own origin/centre. Senders *and* receivers are both restricted to the
tiles, and β̂_obs is re-profiled on that restricted cell set, so the ratio is
apples-to-apples — this is Phase 4's design, transposed.

**N3-occ / N4-occ — whole-section move, occupancy-rejected.** A move is accepted only
if ≥95 % of the shifted senders land in an occupied 25 µm cell. Rejection sampling is
not needed: on the 25 µm lattice the in-occupancy fraction of *every* candidate offset
is one circular cross-correlation, so a single FFT gives the exact admissible set for
all 38,080–108,375 candidate offsets at once. Rotations are screened over 720 angles
directly. `N3-occ15`/`N4-occ15` relax the tolerance to 15 % and are reported because
the 5 % criterion turns out to be degenerate (§3).

**N3-swap — senders relocated to random real cell positions.** Drawn without
replacement from the sender-eligible cells (`Low_quality`, `Unknown`, `Proliferating`
excluded, as everywhere else). In tissue by construction, sender count preserved.
**N3-swap does not preserve sender clustering**, so it is not a torus shift; see §3.
The literal swap is orientation-free and therefore has no distinct rotation analogue,
so **N4-swap is defined as rotate-then-snap**: rotate about the centroid, then move each
rotated sender to the nearest real cell position. `N3-snap` is the same construction
applied to the translation and is reported as a supplement, because it is the only
variant that is both a full-section move *and* in tissue by construction *and*
clustering-preserving.

---

## 3. The two variants where SF collapses, and why it is not the void

### N3-occ / N4-occ: the specified 5 % criterion is degenerate on this geometry

Requiring 95 % of shifted senders to stay in tissue admits **1 to 63 of 38,080–108,375
candidate offsets** (0.001 %–0.17 %), and **the admissible offsets are all
near-identity**:

| section | admissible offsets (N3-occ) | of | median displacement | admissible angles (N4-occ) | of 720 | median displacement |
|---|---|---|---|---|---|---|
| 7259 | 63 | 38,080 | 71 µm | 12 | 720 | 52 µm |
| 7260 | 6 | 90,168 | 21 µm | 3 | 720 | 10 µm |
| **7001** | **1** | 108,375 | **0 µm** | **1** | 720 | **0 µm** |
| 7248 | 9 | 99,224 | 30 µm | 3 | 720 | 18 µm |
| 7352 | 9 | 74,178 | 25 µm | 5 | 720 | 32 µm |
| 7435 | 12 | 95,571 | 31 µm | 7 | 720 | 35 µm |

For section 7001 the only offset that keeps 95 % of senders in tissue is **the identity**,
and its measured surviving fraction is exactly −0.000. That is not an artifact being
removed; that is a null that does not move anything, so β̂_null = β̂_obs and SF → 0 by
construction. Median displacement over all six sections is **27 µm (N3-occ) and 25 µm
(N4-occ)** — *inside* the 100 µm fitting window, and comparable to λ̂ itself (median λ̂
12.8 µm, IQR 7.0–48.2).

**The finding here is about the null, not about the effect: on a non-convex tissue
section a whole-section torus shift cannot simultaneously stay in tissue and displace
anything. The two requirements are incompatible, and the 5 % criterion resolves the
incompatibility by silently choosing not to move.** Relaxing to 15 % (N3-occ15,
displacement 304 µm, 96.6 % retention) recovers **SF = 0.951**, and this is the honest
reading of the occupancy family.

### N3-swap: it is not a torus shift, it is N1

N3-swap scatters senders over uniformly chosen real cell positions, which destroys
sender clustering as well as alignment. So it should behave like a label permutation,
and it does — to a degree that leaves no room for interpretation:

* N1 (cell-type-stratified label permutation), published: **median SF 0.716**
* N3-swap, this run: **median SF 0.721**
* per-fit Spearman ρ(N1_sf, N3_swap_sf) = **0.948**; median |difference| = **0.0087**

The 0.28 that N3-swap removes and N3 does not is **the sender clustering and the
composition structure, not the void**. N1 already reports it and `CS_PHASE3.md` §0
already lists it. N3-swap is a useful confirmation that N1 and the shift nulls are
measuring different things; it is not a corrected N3.

---

## 4. Old vs new, per section

Median SF over the reportable fits in each section (`perm_nulls_c1.csv`):

| variant | 7001 | 7248 | 7259 | 7260 | 7352 | 7435 |
|---|---|---|---|---|---|---|
| N3 original (re-run) | 1.003 | 0.999 | 1.004 | 0.995 | 1.001 | 1.003 |
| **N3-tile** | 0.975 | 0.980 | 1.026 | 0.914 | 0.989 | 0.960 |
| N3-occ *(degenerate)* | −0.000 | 0.395 | 0.902 | 0.276 | 0.688 | 0.742 |
| N3-occ15 | 0.953 | 1.043 | 0.989 | 0.789 | 0.947 | 0.978 |
| N3-swap *(= N1)* | 0.822 | 0.678 | 0.744 | 0.429 | 0.696 | 0.906 |
| N3-snap | 1.020 | 0.903 | 1.015 | 0.973 | 0.970 | 0.999 |
| N4 original (re-run) | 0.937 | 1.024 | 0.949 | 0.822 | 1.074 | 1.000 |
| **N4-tile** | 0.882 | 0.976 | 1.065 | 0.928 | 1.006 | 0.936 |
| N4-occ *(degenerate)* | 0.000 | 0.179 | 0.726 | 0.206 | 0.434 | 0.650 |
| N4-occ15 | 0.793 | 1.012 | 0.960 | 0.756 | 0.907 | 0.997 |
| N4-swap | 0.990 | 0.955 | 1.030 | 0.814 | 0.984 | 1.011 |

The re-run of the original null reproduces the published value (median 1.002 against
1.000; 0.960 against 0.964). It is not bitwise identical and was not expected to be:
the corrected runner draws from the same seed but consumes the random stream in a
different order, so this is a statistical, not an exact, reproduction. Rejection rates
at p < 0.05 are also unchanged: 0.88 for N3 original, 0.83 for N3-tile, against the
published 0.875.

---

## 5. Figure 4 de-hardcoded, and the proof it is faithful

`make_figure4.py` line 42 carried

```python
OURS = dict(N3_sf=1.0001, N4_sf=0.9641, N1_sf=0.7160,
            N3_reject=0.875, N4_reject=0.900,
            beta_obs=0.01126, N3_null_abs=8.91e-5)
```

It is now `load_ours()`, which reads the three surviving fractions from
`results/phase3/sf_summary.csv` (PRIMARY subset) and re-derives the rejection rates and
the amplitudes from `perm_nulls.csv` joined to `main_fits.csv` over the same reportable
population `summarize_phase3.py` defines. All seven constants come back **exactly**:

```
N3_sf 1.0001120  N4_sf 0.9641234  N1_sf 0.7160375  n_fits 160
N3_reject 0.875  N4_reject 0.900  beta_obs 0.0112607  N3_null_abs 8.9053e-05
```

**Verification, done before any re-run touched the inputs:** regenerating the figure
with the derived constants rounded to the 4 decimals at which they had been hardcoded
gives `figure4.png` and `figure4_data.csv` that are **byte-identical** to the committed
ones (md5 `000f34051112aff4fed293fe7a5b25c2` and `bc6fb1d6c2d05dc36a278d5c6b44e1b9`).
Without that rounding the only difference in `figure4_data.csv` is the two
`our SASP kernel estimator` rows, which now carry full precision — max absolute change
**2.34 × 10⁻⁵**, no other cell in the 138-row table altered. (`figure4.pdf` differs on
every regeneration regardless; matplotlib stamps a creation date into it.)

The committed figure is unchanged in content because `sf_summary.csv` is unchanged: the
corrected nulls went to `sf_summary_c1.csv`. When Phase 7 §16 step 3 re-runs M1 and
regenerates `sf_summary.csv`, Figure 4 will now follow it.

---

## 6. What this means for the paper

1. **Contribution 3 survives as written.** The N3 surviving fraction is 1.000 as
   published and 0.974 when the shift is confined to solid tissue; the rotation is 0.964
   published and 0.962 on tiles. The claim that the torus shift and the rotation are
   *calibration failures rather than checks* does not depend on the bounding-box defect.
2. **Say so with the diagnostic attached.** The paper should now state that the
   whole-section shift did put 35.5 % of shifted senders outside the tissue
   (and left 22.8 % with no real neighbour inside the fitting window), that this was
   corrected three ways, and that the result did not move. That is a stronger sentence
   than the one currently in `CS_PHASE3.md`, and it pre-empts the reviewer who reads
   `CS_PHASE4.md` §2.4.
3. **The occupancy result is a finding in its own right, and it is transferable.** On a
   non-convex section, "torus shift, but keep the senders in tissue" is not
   simultaneously satisfiable: at a 5 % out-of-tissue tolerance the admissible offsets
   are near-identity, and for one of six sections the only admissible offset is the
   identity. Anyone running a torus shift on non-convex tissue is choosing between a
   null that leaves the tissue and a null that does not move — **tiling is the way out,
   and this is the concrete argument for it.** It belongs in the Methods, and it is
   §3's "warning to anyone running torus shifts on non-convex tissue" arriving from the
   *other* direction than expected.
4. **Do not present N3-swap as a corrected torus shift.** It reproduces N1 to
   ρ = 0.948 and to a median absolute difference of 0.0087. Reporting it as an N3
   variant would double-count N1 and would misattribute a composition/clustering effect
   to the void correction.
5. **Phase 7 §17's row `SF, N3 corrected (tile / occ / swap)` should read
   `0.974 / 0.951* / 0.721†`** with the footnotes `* at a 15 % out-of-tissue tolerance;
   the specified 5 % tolerance admits only near-identity offsets (SF 0.349) and is not a
   null` and `† not clustering-preserving; equals N1 (0.716)`.
6. **Outcome E of §18 did not occur.** C1 does not change M1's N3 result and no entry is
   needed in the correction ledger beyond the record that the check was run.

---

## 7. Reproduce

```bash
cd /workspace/code
python3 -u phase3_null_diag.py                                   # ~2 min
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband \
        --calls tierA_p95 --n-perm 1000 --n-jobs 6               # 569-671 s/section, 6 in parallel
python3 -u summarize_phase3_c1.py
python3 -u make_figure4.py
```

`MASTER_SEED = 20260820` as everywhere else; `--stage perm_c1` derives its per-section
seeds the same way `--stage perm` does. cKDTree throughout, no (n, n) matrix, ~0.5 MB of
new CSV.

### Notes worth carrying forward

* **The FFT trick.** "Which translations keep ≥ x % of a point set inside a mask?" is a
  circular cross-correlation of the point-count grid with the mask. One `rfft2` answers
  it for every offset at once, exactly, in milliseconds — where naive rejection sampling
  would need ~10³ draws per accepted offset and would never have revealed that the
  admissible set is 1 offset wide.
* **`tiles_for` reuses cleanly** on a Phase 3 section through a shim exposing only
  `coords`, `lo`, `hi`; there was no reason to reimplement Phase 4's tile chooser and no
  reason to import `phase4_data`.
* **A null that guarantees "in tissue" by construction is not automatically a good
  null.** Three of the eleven variants here are 100 % in tissue and two of those are
  degenerate or are a different null altogether. Both extra diagnostic columns —
  `median_displacement_um` and `n_admissible_moves` — exist because the retention
  fraction alone would have hidden it, and the honest table needs all three.
* `N3-snap` keeps 100 % of senders in tissue but drops the median neighbour count from
  140.5 to 120.0, because snapping piles the void-bound senders onto the cells at the
  tissue edge. Retention and neighbour count are not the same diagnostic.

---

## 8. What I did NOT do

* **Only the primary sender call.** `perm_nulls_c1.csv` covers `tierA_p95` on the six
  in-band sections. The N7 sender-definition axis (`perm_nulls_n7.csv`, five further
  calls at 200 permutations) has not been re-run against the corrected nulls.
* **No covariate-adjusted null.** The original `--stage perm` also emits
  `*_full_sf` (β under the full N5+N6 design at fixed λ). The C1 stage computes the
  intercept-only version only, which is the one `sf_summary.csv` and `CS_PHASE3.md` §0
  report. Adding it costs another factor of ~2 in the permutation loop.
* **Figure 2c is not regenerated.** It draws the N3/N4 bands from `perm_curves.csv`,
  which the C1 stage does not write. Phase 7 §19 wants Figure 2 revised with the
  corrected nulls and the destructiveness diagnostic; that is a figure job, not a
  compute job, and the numbers it needs are all in `sf_summary_c1.csv` and
  `null_destructiveness.csv`.
* **`sasp_nulls.py` (synthetic tissue) is untouched.** Its `torus_shift`/`rotate` wrap on
  a square simulation window that *is* fully occupied, so the defect does not arise
  there and the Phase 1/2 size-and-power numbers are unaffected.
