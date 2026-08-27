# CS Phase 8 — C1 closeout: the null battery is now freezable

**Status: complete. Verdict: the C1 conclusion is invariant to the sender
definition and to covariate adjustment. N3/N4 can be frozen as corrected.**

Companion to `reports/CS_PHASE7_C1.md`, which is unchanged and still carries the
primary C1 result. This report closes the four gaps that report's §8 listed, and
answers the two scheduling questions the coordinator asked. New artefacts take
Phase 8 naming; the `*_c1` result-file family keeps its name, because these are
further outputs of the same runner and renaming them would break the
cross-references already verified.

| Deliverable | Location |
|---|---|
| N7 sender-definition axis, corrected nulls | `results/phase3/perm_nulls_c1_n7.csv`, `sf_summary_c1_n7.csv` |
| Covariate-adjusted `*_full_sf` family | `results/phase3/perm_nulls_c1.csv` (23 new columns), `sf_summary_c1.csv` |
| N3-swap ≡ N1, all six sender calls | `results/phase3/sf_summary_c1_swap_vs_n1.csv` |
| Binned null curves, corrected nulls | `results/phase3/perm_curves_c1.csv`, `perm_curves_c1_n7.csv` |
| Figure 2b, 2c revised; **2e–h new** | `figures/figure2{b,c,e}.{png,pdf}`, `results/phase3/figure2{c,e}_data.csv` |
| Runner and figure changes | `code/run_phase3_nulls.py`, `code/summarize_phase3_c1.py`, `code/make_figure2bc.py` |

**Pinned files, re-verified after everything below:**

```
3b77aa1bba0712c205c5d9356654fb71  results/phase3/perm_nulls.csv
69e3a1d3f60060deddcceba9896a7d31  results/phase3/sf_summary.csv
ecf86b9ca5460f31290e2f4c9e822ea2  results/phase3/summary_phase3.txt
```

**And the C1 numbers already reported did not move.** `--stage perm_c1` was
extended with the covariate-adjusted fit and the curve accumulation without
disturbing the order in which it consumes the random stream, so the re-run
reproduces all 22 previously reported columns of `perm_nulls_c1.csv`
**bit-identically** across all 595 rows (max absolute difference 0.0, `beta_obs`
identical). Every number in `CS_PHASE7_C1.md` still stands as written.

---

> ### Correction box — audit item R9, added 2026-08-27
>
> `results/phase3/null_destructiveness.csv` was regenerated at 05:58:41 during
> the M1 re-run (task 8.7) under the promoted C6 gene sets, i.e. **after this
> report was written**. Section-geometry columns are unchanged; every
> sender-dependent number moved slightly (N3-orig displacement 2,974 → 2,910 µm,
> N3-tile 489 → 479, N4-orig retention 0.917 → 0.920, real median neighbours
> 140.5 → 140.0). **The §4.1 row "0.772 → 1.000" holds exactly** — N3 bounding-box
> retention is 0.772 before and after. The surviving fractions this report
> tabulates are superseded by `reports/CORRECTIONS.md` §8 (primary), §10 (the N7
> axis at 1,000 permutations) and §11 (the second Tier A variant); the
> conclusions are unchanged.

---

## 1. N7 — the correction is invariant to the sender definition

This was the gap that mattered for the freeze, because Phase 7's primary sender
definition is not Phase 3's. Six sender calls, six in-band sections, 200
permutations, all eleven null variants (`sf_summary_c1_n7.csv`). Median
surviving fraction over each call's own reportable population:

| sender call | reportable fits | sender prevalence | N3 orig | **N3-tile** | N3-occ | N3-occ15 | N3-swap | N3-snap |
|---|---|---|---|---|---|---|---|---|
| `tierA_p90` | 198 | 8.1–9.0 % | 0.999 | **0.989** | 0.673 | 0.961 | 0.811 | 1.013 |
| `tierA_p95` *(primary)* | 160 | 4.0–4.5 % | 1.002 | **0.974** | 0.349 | 0.951 | 0.721 | 0.988 |
| `tierA_p99` | 91 | 0.8–0.9 % | 0.992 | **0.980** | 0.252 | 0.920 | 0.755 | 0.984 |
| `cdkn1a_pos` | 203 | 1.7–5.8 % | 0.999 | **0.990** | 0.419 | 0.981 | 0.893 | 1.009 |
| `senepy_p95` | 155 | 2.5–3.7 % | 0.996 | **1.006** | 0.406 | 0.993 | 0.926 | 0.999 |
| `senepy_p99` | 107 | 0.5–0.7 % | 1.006 | **1.005** | 0.445 | 0.987 | 0.933 | 0.988 |

| sender call | N4 orig | **N4-tile** | N4-occ | N4-occ15 | N4-swap |
|---|---|---|---|---|---|
| `tierA_p90` | 0.968 | **0.990** | 0.489 | 0.920 | 0.996 |
| `tierA_p95` | 0.960 | **0.962** | 0.273 | 0.896 | 0.969 |
| `tierA_p99` | 0.958 | **0.963** | 0.121 | 0.891 | 0.964 |
| `cdkn1a_pos` | 1.002 | **0.962** | 0.350 | 0.955 | 1.011 |
| `senepy_p95` | 1.009 | **0.976** | 0.409 | 0.951 | 1.016 |
| `senepy_p99` | 0.999 | **0.995** | 0.315 | 0.890 | 1.014 |

**Three statements the freeze can rely on.**

1. **The in-tissue correction changes nothing, for any sender definition.**
   N3-tile spans 0.974–1.006 and N4-tile 0.962–0.995 across a sender prevalence
   range of **0.5 % to 9.0 %** and across three unrelated sender callers
   (curated Tier A percentiles, a single-gene call, SenePy hub scores). The
   published bounding-box values span 0.992–1.006 and 0.958–1.009 over the same
   grid. There is no sender definition for which the void mattered.
2. **N3-occ is degenerate for every sender call** — 0.252–0.673 for N3-occ and
   0.121–0.489 for N4-occ — always with a median displacement below or barely
   above λ̂. It is a property of the tissue
   geometry, not of the sender set.
3. **N3-swap tracks N1 for every sender call, and most tightly for the Tier A
   percentile calls** — which is where the Phase 7 primary sender definition
   lives (`sf_summary_c1_swap_vs_n1.csv`):

| sender call | N1 | N3-swap | Spearman ρ per fit | median abs. difference |
|---|---|---|---|---|
| `tierA_p90` | 0.777 | 0.811 | 0.922 | 0.0107 |
| `tierA_p95` | 0.716 | 0.721 | **0.948** | **0.0087** |
| `tierA_p99` | 0.728 | 0.755 | **0.983** | 0.0252 |
| `cdkn1a_pos` | 0.911 | 0.893 | 0.830 | 0.0397 |
| `senepy_p95` | 0.895 | 0.926 | 0.434 | 0.0711 |
| `senepy_p99` | 0.886 | 0.933 | 0.511 | 0.0548 |

   Reported against interest: the identity is **not** uniform. It is near-exact
   for the three Tier A calls and only directional for SenePy (ρ 0.43–0.51),
   because SenePy's senders are more concentrated within particular cell types,
   so permuting labels *within* type (N1) and relocating senders *across* types
   (swap) stop being the same operation. The correct statement is therefore:
   **N3-swap is a label permutation, stratified by cell type only to the extent
   that the sender call already is.** It is still not a torus shift.

---

## 2. The covariate-adjusted `*_full_sf` family, and what it settles

`--stage perm_c1` now also reports β̂ under the full N5+N6+zonation design at
fixed λ — the `*_full_sf` family `--stage perm` has always emitted. All eleven
variants, primary call, 160 reportable fits (`sf_summary_c1.csv`):

| variant | SF (intercept only) | **SF (full N5+N6+zonation design)** |
|---|---|---|
| N3 original, published | 1.000 | 0.991 |
| N3 original, re-run | 1.002 | 0.997 |
| **N3-tile** | 0.974 | **0.990** |
| N3-occ | 0.349 | 0.410 |
| N3-occ15 | 0.951 | 0.986 |
| **N3-swap** | **0.721** | **0.999** |
| N3-snap | 0.988 | 0.994 |
| N4 original, published | 0.964 | 0.940 |
| N4 original, re-run | 0.960 | 0.945 |
| **N4-tile** | 0.962 | 0.963 |
| N4-occ | 0.273 | 0.267 |
| N4-occ15 | 0.896 | 0.934 |
| N4-swap | 0.969 | 0.964 |

**The new fact is the N3-swap row.** Conditioning on the nuisance block moves
N3-swap from 0.721 to **0.999** — the entire 0.28 that N3-swap removes is
absorbed by the covariates that N5 and N6 already model. It removes nothing the
nuisance model does not already remove. That is the third independent line of
evidence, after ρ = 0.948 against N1 and the sender-call table above, that
N3-swap's drop is composition and clustering, not spatial alignment. It also
mirrors `CS_PHASE3.md` §0's row "N1 applied to the N5+N6-conditioned residual =
0.987": N1 behaves the same way for the same reason.

N3-occ and N4-occ stay low under the full design (0.410, 0.267), as they must —
covariate adjustment cannot rescue a null that does not move.

---

## 3. Figure 2, fully revised and audited panel by panel

| panel | script | depends on the nulls? | outcome |
|---|---|---|---|
| **2a** | `make_phase5_figs.py --which 2a` | no | **regenerated, byte-identical** — `69453a1ec3fe58a0411b3de5cc235009` |
| **2b** | `make_figure2bc.py` | yes | **revised** — in-tissue null bands added |
| **2c** | `make_figure2bc.py` | yes | **revised** — nine corrected variants below a divider |
| **2d** | `make_figure2bc.py` | no | **regenerated, byte-identical** — `5ca89780eac1bedee3a6bdcbe0434125` |
| **2e–h** | `make_figure2bc.py` | yes | **new** — the destructiveness audit and the pending A7 slot |

Current PNG hashes: 2b `df2f0afdc3264e7e36693a5cd542c15e`, 2c
`983b47b2d97b4aa490b9903803b68f0b`, 2e `98534d630e2ac16aee94ff727fe464b7`,
figure 4 `1ca9459f4706da45c2d15256705c56a4`. PDFs are excluded from every
reproducibility comparison here because matplotlib stamps a creation date into
them; PNGs are the check.

**A trap found while doing this, worth recording.** Figure 2a has *two*
producers: `make_figure2.py` (the superseded Phase 2 version) and
`make_phase5_figs.py --which 2a` (the live, receiver-type-stratified one). Both
write `figures/figure2a.png`. Running `make_figure2.py` silently replaced the
committed figure with the stale version — its two data CSVs came back
byte-identical, so nothing warned. The committed file was restored from git and
then reproduced byte-identically from `make_phase5_figs.py`, which is the
producer of record. **Do not run `make_figure2.py` during the M1 re-run.**

**2b — the naive gradient with in-tissue bands.** The N3 bounding-box 95 % band
is kept as a filled grey region and two in-tissue bands are drawn over it as
dotted envelopes: N3-snap (whole-section shift, 100 % in tissue,
clustering-preserving) and N3-occ15. All three lie on top of each other and on
the observed and matched-decoy curves in every panel — which is the same
statement `CS_PHASE3.md` made, now made with a null that stays in the tissue.
The tile variants are deliberately *not* drawn here: their fit is restricted to
solid tiles, so their binned curve sits on a different receiver set and would
not be comparable on these axes. They appear in 2c and 2g, where the ratio is
scope-matched.

**2c — the surviving-fraction dot plot.** The ten published rows are unchanged
and unmoved; the nine corrected variants sit below a labelled divider. `N3-occ`
and `N4-occ` carry `(degenerate)` in the row label and `N3-swap` carries
`(= N1)`, so the two misreadable numbers are annotated on the axis itself, not
only in the caption.

**2e–h — the destructiveness audit.** Panel **e** is the retention fraction with
the two published bounding-box rows in red, bold, and separated by a rule from
the in-tissue rows: 0.772 against 1.000. Panel **f** is the median displacement
on a log axis with **two reference lines drawn — the median λ̂ (12.8 µm) and the
100 µm fitting window** — which is what makes N3-occ's 27 µm and N4-occ's 25 µm
legible as disqualifying rather than merely small. Panel **g** is old-vs-new
surviving fraction, bar = median, black line = IQR, **◆ = the published value**,
with a dash-dot line at N1 = 0.716 and an arrow labelling N3-swap as ≡ N1.
Panel **h** is a **declared placeholder** for the negative-control-probe kernel:
the axes are drawn, the title reads `PENDING`, and the box states the test, the
pass condition, and why it is not filled.

Every plotted number is in `results/phase3/figure2e_data.csv`, including the two
reference-line values and an explicit `negative_control_probe_kernel /
family=PENDING` row so the empty slot is in the data table too. Styling is
`sasp_palette.apply_style` throughout; no new palette, fonts or idiom.

**One thing panel h can stop waiting for.** The M1 bundles already ship the
probes — `data/raw/<sample>/cell_feature_matrix.h5` carries **40 Negative
Control Probe, 609 Negative Control Codeword and 21 Genomic Control** features
alongside the 5,106 Gene Expression ones. Test A7 is behind the freeze only for
H1. The mouse half is runnable now, needs no download, and would fill half of
panel h before the pre-registration rather than after. I have not run it — it is
a new analysis, not one of the four gaps — but it is a few hours at most and I
would take it if it is wanted.

---

## 4. Proposed §17 rows and Figure 2 caption, as text

### 4.1 §17, the two-arm comparison table

Replace the single placeholder row with these three, verbatim:

| Quantity | M1 mouse liver (SBR fibrosis) | H1 human |
|---|---|---|
| **SF, N3 corrected (tile / occ / swap)** | **0.974 / 0.951\* / 0.721†** | |
| **SF, N4 corrected (tile / occ / swap)** | **0.962 / 0.896\* / 0.969** | |
| **Null destructiveness, N3 published → corrected** | **0.772 → 1.000** shifted senders keeping a real neighbour within 100 µm‡ | |

> **\*** N3-occ and N4-occ as specified — reject any move that puts more than 5 %
> of the shifted senders outside the 25 µm occupancy grid — are **not
> implementable on a liver section**. The criterion admits 1 to 63 of 38,080 to
> 108,375 candidate translations and 1 to 12 of 720 candidate rotations, all
> near-identity: median displacement 27 µm (N3) and 25 µm (N4), inside the
> 100 µm fitting window and close to the median λ̂ of 12.8 µm. For section 7001
> the only admissible translation, and the only admissible rotation, is the
> identity, whose surviving fraction is 0.000 by construction. The values quoted
> here are the same nulls at a 15 % tolerance (median displacement 304 µm and
> 308 µm; 96.6 % and 96.9 % of senders retaining a real neighbour within
> 100 µm). At the literal 5 % tolerance the values are **0.349 and 0.273**, and
> they measure the null's degeneracy, not the effect. This holds for all six
> sender definitions (0.121–0.673).
>
> **†** N3-swap relocates each sender to a uniformly chosen real cell position,
> which destroys sender clustering as well as sender–response alignment. **It is
> therefore not a torus shift and 0.721 is not a corrected N3.** It reproduces
> N1, the cell-type-stratified label permutation: median SF 0.721 against N1's
> 0.716, per-fit Spearman ρ = 0.948, median absolute difference 0.0087; and
> conditioning on the N5+N6+zonation nuisance block moves it to **0.999**, i.e.
> it removes nothing the nuisance model does not already remove. The
> clustering-preserving in-tissue variants are N3-tile (0.974) and N3-snap
> (0.988). The rotation family's swap is defined as rotate-then-snap, since a
> relocation to random real positions is orientation-free.
>
> **‡** Measured over the six admissible sections, 20 draws each, against every
> real cell position (`results/phase3/null_destructiveness.csv`). The published
> bounding-box torus shift left **35.5 %** of shifted senders outside the tissue
> (`1 − frac_in_occupancy`) and left **22.8 %** with no real cell at all inside
> the 100 µm fitting window (`1 − frac_retaining_a_neighbour`), thinning the
> median sender neighbourhood from 140.5 real cells to 119.7; for the rotation
> the two figures are **19.9 %** and **8.0 %**, 140.5 → 129.3. *(Audit item R3:
> this footnote previously gave 23 % / 8 % as the out-of-tissue fraction. Those
> are the lost-every-neighbour figures. The two must not be conflated.)* All in-tissue variants retain
> **96.6–100 %**, with median neighbour counts of 120.0–142.4 — the low end being
> the two snap variants, which pile void-bound senders onto tissue-edge cells and
> so keep every sender in tissue at the cost of thinning its neighbourhood.
> N3-tile is 139.7 against a real 144.0. **This is a geometric quantity** in the §17 sense —
> it depends on the section outline and the sender point pattern, not on
> biology — so it should be reported for H1 and should behave the same way for
> any non-convex section.

### 4.2 Figure 2 caption, revised

> **Figure 2 | The naive distance gradient and the null battery, with the
> coordinate nulls confined to tissue.**
> **(a)** Naive binned response versus distance to nearest sender, by receiver
> cell type, six Test-3-admissible sections. **(b)** The naive gradient for two
> Tier B modules in hepatocytes, with the matched-decoy (N2) curve and the
> torus-shift null band: the published bounding-box band (grey fill) and two
> in-tissue bands (dotted) that keep every shifted sender in the section. All
> three coincide with the observed curve. **(c)** Surviving fraction of β̂ under
> each null over the 160 reportable fits (positive naive amplitude, spatial block
> bootstrap CI excluding zero); the ten published nulls above the divider, the
> nine in-tissue N3/N4 variants below it. **(d)** Median distance to nearest
> sender versus sender density against the homogeneous-Poisson slope of −1/2.
> **(e)** Fraction of shifted senders retaining a real cell within the 100 µm
> fitting window: the published bounding-box nulls (red) leave **22.8 %** of
> N3's senders and **8.0 %** of N4's with no real cell inside that window — and,
> separately, **35.5 %** and **19.9 %** of them outside the tissue altogether;
> every in-tissue variant retains
> 96.6–100 %. **(f)** Median distance each null actually moves a sender, with the
> median λ̂ (12.8 µm) and the 100 µm window marked: the occupancy-screened nulls
> at the specified 5 % tolerance move senders *less far than λ̂*, so they are in
> tissue by being near-identity. **(g)** Surviving fraction, published value (◆)
> against the corrected variants; N3-swap coincides with N1 (dash-dot line) and
> is a label permutation, not a torus shift. **(h)** *Pending:* the
> negative-control-probe kernel (Section 13, test A7) — counts of negative
> control probes against distance to nearest sender, which must be flat. Requires
> the H1 arm and is held behind the pre-registration freeze.

The negative-control-probe slot is the only part of §19's Figure 2 brief not
delivered, and it is marked pending in the figure, in the caption, in
`figure2e_data.csv` and in the docstring of `fig2e()`.

---

## 5. Answers to the two scheduling questions

### 5.1 Does the incoming M1 comparability re-run change how these four gaps were closed?

**No, and N7 in particular was worth doing now rather than after.** The
justification is in §1: what the N7 run establishes is that the *geometric*
conclusion — the void mattered, the correction does not move the answer, the 5 %
occupancy criterion is degenerate, the swap is a label permutation — holds
across sender prevalences from 0.5 % to 9.0 % and across three unrelated sender
callers. Those are properties of the section outline and the sender point
pattern. A new mouse strict Tier A at p90/p95/p99 will produce a sender set
whose prevalence lands inside the band already covered, so the *conclusion*
transfers.

What will need redoing is the **numeric SF values**, because both the sender
flags and the module scores change. That is the M1 re-run, and it was always
going to include the N3/N4 columns. Concretely: re-running N7 now costs 10.7
minutes and buys the evidence the freeze needs; not re-running it would have
left the battery frozen on an untested axis. Nothing here has to be thrown away.

One caveat worth pre-registering: the N3-swap ≡ N1 identity is tight for Tier A
percentile calls (ρ 0.92–0.98) and loose for SenePy (ρ 0.43–0.51). If the Phase 7
per-module sender sets are much more cell-type-concentrated than Tier A
percentiles, expect the swap and N1 to separate further. That does not affect
the freeze — the swap is not the corrected N3 either way — but the footnote in
§4.1 should say "reproduces N1 for the Tier A calls" rather than "reproduces N1",
and it does.

### 5.2 What does an M1 end-to-end re-run cost?

Wall-clock on this box at the modest parallelism used throughout (n_jobs ≤ 6 for
the permutation stages, ≤ 22 for the fits), with the DeepScence job finished.
"Source" says where each number comes from; nothing here is a guess dressed as a
measurement.

| stage | script | re-run? | wall | source |
|---|---|---|---|---|
| Xenium ingest | `prepare_samples.py` | **no** | — | raw data unchanged |
| cell-type annotation | `annotate_pipeline.py` | **no** | — | cell-type markers are independent of Tier A/B |
| DeepScence | `run_deepscence_all.py` | **no** | — | D1 just finished at 11/11; it is an input |
| Phase 4 CCC battery | `phase4_run.py` and friends | **no** | — | depends on Tier C pairs and cell types only |
| anatomy / sender calls / module scores | `phase2_downstream.py` | **yes** | **~14 min** | log mtimes, 11 sections in 4 batches, 17:49→18:03 |
| Phase 3 cache | `sasp_phase3.prep` | **yes** | **~1 min** | measured today: 8.3 s (127 k cells), 16.4 s (225 k cells) per section, 6-way |
| window | `--stage window` | yes | ~2 min | one tree query per section × call; not separately timed |
| main fits | `--stage main` | **yes** | **16.9 min** at 66 jobs; **~25 min** at 99 | `logs/phase3_main3.log`, n_jobs=22 |
| published bbox nulls | `--stage perm` | optional | ~16 min | `logs/phase3_perm2.log`, 831–967 s/section, 6-way |
| **corrected nulls, primary** | `--stage perm_c1` | **yes** | **20.2 min** | measured today, 1,000 perms, `do_full` + curves |
| **corrected nulls, N7** | `--stage perm_c1` | **yes** | **10.7 min** at 30 jobs; **~17 min** at 48 | measured today, 200 perms |
| curves | `--stage curves` | yes | ~2 min | `logs/phase3_curves.log` |
| N8, stratification, attribution, combined, Poisson, λ-scale, Ripley, correlogram | eight scripts | yes | **~35 min** total | log mtimes 18:27→19:03 |
| destructiveness diagnostic | `phase3_null_diag.py` | yes | ~2 min | measured today |
| summaries and Figures 2a–2e, 4 | four scripts | yes | ~5 min | measured today |
| Phase 5 kernels and superposition | `run_phase5_*` | **yes** | ~12 min | log mtimes 19:56→20:13 |

**Total: about 2 h 00 m** as scoped today (summing the "yes" rows), and **about
2 h 30 m** if both Tier A variants expand the sender-call axis from 6 calls to 9
(main 66→99 jobs, N7 30→48 jobs) *and* the published bounding-box nulls are
re-run alongside the corrected ones for the old-vs-new table.

**Schedule it as half a day, not a day**, with three caveats:

1. **The critical path starts at the bio agent.** `phase2_downstream.py` cannot
   start until the re-sourced mouse B7 and the rebuilt mouse strict Tier A are in
   `genesets/`. Everything else is strictly downstream of it. The compute is not
   the constraint; the gene sets are.
2. **Both Tier A variants are cheap.** Adding the second sender definition costs
   about 15 minutes across `main` and `perm_c1`, not a second pass.
3. **Do not run `make_figure2.py`** (§3). And Phase 4 does not need to re-run,
   but `make_figure4.py` should, since it now reads its constants from
   `results/phase3/sf_summary.csv` and will follow the new values automatically.

---

## 6. Correction to `CS_PHASE7_C1.md` §8

That report's §8 said "Figure 2c is not regenerated. It draws the N3/N4 bands
from `perm_curves.csv`". **That is wrong on both halves.** Figure 2c is the
surviving-fraction dot plot and reads `perm_nulls.csv`; it is **Figure 2b** that
consumes `perm_curves.csv` for the null band. Both are now regenerated, and the
`perm_curves_c1.csv` pass that 2b needed is written by `--stage perm_c1`.
`CS_PHASE7_C1.md` is left byte-unchanged as instructed; the correction lives
here.

The other three §8 items are closed by §1 and §2 above. The fourth —
`sasp_nulls.py` on synthetic tissue — stands: its window is a fully occupied
square, the defect does not arise, and the Phase 1/2 size-and-power numbers are
unaffected.

---

## 7. Reproduce

```bash
cd /workspace/code
python3 -u phase3_null_diag.py
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband \
        --calls tierA_p95 --n-perm 1000 --n-jobs 6                    # 20.2 min
python3 -u run_phase3_nulls.py --stage perm_c1 --sections inband \
        --calls tierA_p90,tierA_p99,cdkn1a_pos,senepy_p95,senepy_p99 \
        --n-perm 200 --n-jobs 6                                       # 10.7 min
python3 -u summarize_phase3_c1.py
python3 -u make_figure2bc.py          # 2b, 2c, 2d, 2e-h
python3 -u make_phase5_figs.py --which 2a
python3 -u make_figure4.py
# NOT make_figure2.py -- it is the superseded producer of figure2a (section 3)
```

### Verification performed

* `perm_nulls_c1.csv` re-run with the two new features: **all 22 previously
  reported columns bit-identical over 595 rows.**
* `figures/figure2a.png` → `69453a1ec3fe58a0411b3de5cc235009`, and
  `figures/figure2d.png` → `5ca89780eac1bedee3a6bdcbe0434125`: **byte-identical
  to the committed versions**, confirming neither depends on the nulls.
* `results/phase3/{perm_nulls,sf_summary,summary_phase3}` md5s unchanged
  (top of this report).
* No commits, no tags. Nothing under `code/run_deepscence*`,
  `data/processed/deepscence_*` or `genesets/` was read or written.

### Engineering notes

* **Extending a permutation runner without moving its numbers** means leaving the
  order of `rng` consumption alone. The covariate-adjusted fit and the curve
  accumulation were added *inside* the existing per-null loop, after the draw,
  so every variant still sees the same sequence of offsets and angles. The
  bit-identity check above is what makes that claim rather than assumes it.
* The covariate-adjusted fit costs about **+25 %** wall-clock, not the 10× a
  naive flop count suggests: `FixedLambdaFitter` orthonormalises the design once,
  so each null costs two (n × p) matvecs in BLAS, against an `exp` and two dot
  products for the intercept-only version.
* Figure 2a having two producers cost twenty minutes and a `git checkout`. When
  two scripts write the same output path, the one that ran last wins and nothing
  says so.
