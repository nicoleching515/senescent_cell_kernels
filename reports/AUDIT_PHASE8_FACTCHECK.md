# AUDIT — Phase 8 independent fact-check

**Auditor:** independent verification agent. **Date:** 2026-08-27, 06:15–06:40 UTC.
**Method:** every claim below was re-derived from files on disk with `python3`, not read
from a report. Where a report and a file disagree, the file wins and the command that
settles it is given. Nothing outside this file was written.

**Bottom line.** The gene-set freeze, the panel arithmetic, the C1 null battery, the A7
statistics and the H1 acquisition all reproduce exactly from the data. Six figure/data
chains were re-derived value-by-value and match. **Nine claims are refuted and ten more
are true-but-misleading.** One refuted claim is in `SUBMISSION_PATCH_2026-08-29.md`, i.e.
it is scheduled to reach a manuscript in two days.

---

## 0. Moving-target log

Three concurrent agents rewrote files during this audit. Recorded per the brief.

| File | md5 at audit start (06:15) | md5 at audit end (06:40) | Consequence |
|---|---|---|---|
| `results/phase3/perm_nulls_c1.csv` | `93194e2e0ce2` (04:29:55) | **`0318737e85a9` (06:25:53)** | **Changed mid-audit.** All C1 §4 verifications below used the 04:29 version. |
| `results/phase3/perm_curves_c1.csv` | `8ad03e8fe86a` | `e864a8ab2e84` (06:25:53) | Changed mid-audit. |
| `results/phase3/caller_coverage_gate{,_headline}.csv` | 06:09:34 | 06:20:15 | Deferred per brief. |
| `results/phase3/n8_*` (12 files) | 2026-08-20 | 06:20–06:27 | In flight. |
| `reports/{PHASE8_ROADMAP_STATUS,PREREG_PHASE8,CORRECTIONS,CS_PHASE8_M1_RERUN}.md` | — | all rewritten 06:23–06:34 | Roadmap gained a "D-namespace collision" banner mid-audit. |
| New files appeared | — | `reports/COMPLETED_TASKS.md`, `reports/NOVELTY_ASSESSMENT.md` | Not audited. |
| `results/phase3/null_destructiveness.csv` | `423d9ae3c345` (05:58:41) | unchanged | **Already regenerated at 05:58**, i.e. *after* both C1 reports were written. See R9. |
| `results/phase3/{main_fits,perm_nulls,perm_curves}.csv` | modified vs `HEAD` | still modified | M1 re-run. Pre-re-run versions recovered via `git show HEAD:…` where needed. |

Files verified as *stable* through the audit: `sf_summary_c1.csv`, `sf_summary_c1_n7.csv`,
`sf_summary_c1_swap_vs_n1.csv`, `sf_summary.csv`, `a7_*.csv`, `a7_verdict.txt`, all
`genesets/`, all `figures/*_data.csv`, `figure2c_data.csv`, `figure2e_data.csv`.

---

# 1. REFUTED

### R1 — "a −0.070 SD distance gradient (p = 0.023) in negative-control probes"
**Where:** `reports/SUBMISSION_PATCH_2026-08-29.md` L112–113 (manuscript-bound);
`reports/CS_PHASE8_CALLERS.md` §4.1 statement 1 and §4.5 caption; `PHASE8_ROADMAP_STATUS.md`
row 8.5b ("Raw assay is NOT flat (-0.070 SD, p=0.023)").

**Status: REFUTED as stated.** −0.070 is the `all_controls` response — the pooled sum of
the 40 Negative Control Probes **plus the 609 Negative Control Codewords plus the 21
Genomic Controls**, and the codewords carry 73 % of those counts. The 40 negative-control
probes **on their own are flat**:

| response (design `base`) | clustered mean | 95 % CI | p |
|---|---|---|---|
| `all_controls` (quoted as "−0.070") | −0.0697 | [−0.1276, −0.0118] | 0.023 |
| **`neg_control_probe` (the 40 probes)** | **−0.0177** | **[−0.0453, +0.0099]** | **0.183** |
| `neg_control_codeword` (609) | −0.0549 | [−0.1064, −0.0034] | 0.039 |
| `genomic_control` (21) | −0.0337 | [−0.0538, −0.0136] | 0.0039 |

`python3 -c "import pandas as pd;d=pd.read_csv('results/phase3/a7_summary.csv');print(d[d.design=='base'][['response','clustered_mean','clustered_lo','clustered_hi','clustered_p']])"`

The binned curve tells the same story: `all_controls` runs +0.019 → +0.293 SD across the
window, `neg_control_probe` runs −0.053 → +0.084.

This matters because `PREREG_PHASE8_genesets.md` §11 designates
`E_negative_control_probes` (n = 40) the **"Primary technical null. Audit test A7 fits
these against distance-to-sender and requires flat"**, and Phase 9 item 9.4 repeats it.
**On the pre-registered primary response, M1's A7 passes naively.** The report's own §4.5
caption is honest ("40 negative control probes, 609 negative control codewords and 21
genomic controls, pooled"); the headline, the roadmap and the submission patch are not.

**Fix:** either say "pooled negative-control features" everywhere and amend the
pre-registration, or report the probe-only number (−0.018, p = 0.18) as the A7 result and
the pooled one as a powered supplement. The substantive conclusion ("the raw assay is not
flat") survives on the codewords and genomic controls; the *name* on it does not.

### R2 — "Human-only: CDKN2B DDB2 GADD45G MDC1 PHLDA3 TNFRSF10B XPC"
**Where:** `PREREG_PHASE8_genesets.md` §6 asymmetry table; `results/phase7_jobA/crossarm_geneset_table.json`
→ `asymmetry[A_SENDER_FINAL_strict_FROZEN].human_only`; rendered verbatim into the caption
of `figures/figure_gs2_crossarm_symmetry.png`.

**Status: REFUTED.** `Cdkn2b` is a member of the **mouse** frozen Tier A set. It is listed
as human-only only because the pinned MGI map has no `Cdkn2b` row — a map gap the project
records in its own `genesets/human/_symbol_resolutions.csv`
(`Cdkn2b,CDKN2B,NOT in MGI map; upper-case form is a known HGNC symbol`).

```
grep -c Cdkn2b genesets/A_SENDER_FINAL_strict.txt          # 1  (mouse Tier A)
grep -c CDKN2B genesets/human/A_SENDER_FINAL_strict.txt    # 1  (human Tier A)
grep -c Cdkn2b genesets/mouse_human_orthologs_MGI.csv      # 0  (map gap)
```

Two consequences:
1. The **true shared count is 27 of 33, not 26**, once the same map-gap correction the
   document applies to Tier C (`human_only_map_gap: "CXCL2 CXCL5"`, `human_only_real:
   "CXCL8 IL1B"`) is applied to Tier A. The document warns about exactly this failure mode
   two rows above and then commits it on its headline set.
2. The `mouse_only` list is **arithmetically incomplete**: 33 − 26 = 7 mouse genes are
   unmatched but only 5 are listed. The two omitted are `Cdkn2b` and `H2afx`, the two that
   fail to map. (`H2afx` *is* genuinely mouse-only: human `H2AFX` sits in
   `genesets/human/B_downstream_arrest.txt` and is removed from human Tier A by
   disjointness. That is a real asymmetry and is correctly absent from the human set.)

### R3 — "the whole-section shift put 23 % of shifted senders in the void"
**Where:** `CS_PHASE7_C1.md` §0 ("it is **23 %** for N3 … and **8 %** for N4") and §6.2;
`CS_PHASE8_C1_CLOSEOUT.md` §4.1 footnote ‡ and the §4.2 Figure 2(e) caption
("the published bounding-box nulls (red) leave 23 % of N3's senders and 8 % of N4's
outside the tissue"); `PHASE8_ROADMAP_STATUS.md` ("23% void, not ~20%").

**Status: REFUTED.** 23 % is `1 − frac_retaining_a_neighbour`, i.e. the fraction of shifted
senders with **no real cell within the 100 µm fitting window**. The out-of-tissue fraction
is a *different column of the same file*, `frac_in_occupancy`:

| null | 1 − retention (quoted as "void") | **1 − frac_in_occupancy (actually out of tissue)** |
|---|---|---|
| N3 bounding-box | 22.8 % | **35.5 %** |
| N4 bounding-box | 8.3 % | **19.9 %** |

`python3 -c "import pandas as pd;d=pd.read_csv('results/phase3/null_destructiveness.csv');print(d.groupby('null')[['frac_retaining_a_neighbour','frac_in_occupancy']].median())"`

The report's own §1 table header is correct ("N3: senders keeping a neighbour"); the
translation into "in the void" in §0, §6.2 and the figure caption is not. The bug is
**worse** than the paper is about to claim (35 % vs 23 %), so the correction runs against
the authors' interest in the safe direction — but the sentence as drafted is factually
wrong about what was measured, and Figure 2(e)'s caption states it as a measurement.

### R4 — "figure2b.png `df2f0afd…`, figure2c.png `983b47b2…` … the committed content is intact"
**Where:** `CS_PHASE8_CALLERS.md` §4.4, including the block-quoted anomaly narrative.

**Status: REFUTED.** Those hashes belong to `figures/revised_candidates/figure2b_REVISED.png`
and `figure2c_REVISED.png`, not to `figures/`. The committed and current figures are:

```
$ md5sum figures/figure2b.png figures/figure2c.png
5ecd9ad1029851c1955fc938abf9444c  figures/figure2b.png
a232a529a566d7c0680e04840dd07a9b  figures/figure2c.png
$ git show HEAD:figures/figure2b.png | md5sum   # 5ecd9ad1029851c1955fc938abf9444c
$ md5sum figures/revised_candidates/figure2b_REVISED.png  # df2f0afdc3264e7e36693a5cd542c15e
```

Three sub-claims fail:
- "Panels 2b, 2c and 2d regenerated **byte-identically**" — only 2d matches the repository.
  `figures/revised_candidates/README.md` states plainly that 2b and 2c **do** differ in
  content ("Content differs? **YES**").
- The 05:22:41 rewrite is described as unexplained ("no process of mine ran in that window
  … most likely a concurrent session calling `fig2b()`/`fig2c()` directly"). It is the
  documented PI-directed restore to the committed baseline, written up at 05:23 in
  `figures/revised_candidates/README.md` and again in `PHASE8_ROADMAP_STATUS.md`.
- "so the committed content is intact" — the agent's regeneration was itself reverted at
  05:29; what is on disk is the committed baseline, not what the report says it left.

The repository state is correct and well-documented. **The report's account of it is not**,
and the two documents currently disagree about which figure2b/2c is canonical.

### R5 — "under +N6+N5 every control family's clustered mean is indistinguishable from zero"
**Where:** `CS_PHASE8_CALLERS.md` §4.1 statement 2.

> **⚠ SUPERSEDED 2026-08-27 (record reconciliation) — R5 IS NOW MOOT; THE CLAIM IS TRUE.**
> This item was audited against the **pre-C6 05:19** `a7_summary.csv`. In the frozen 09:06
> file, `neg_probe_rate` under `n6n5` is **+0.0097 [−0.0060, +0.0253], p = 0.199** — the CI
> **includes** zero. All five control families are now indistinguishable from zero under
> +N6+N5, so "every control family is flat under +N6+N5" is correct as written and needs no
> hedge. Do **not** carry forward the "every count-based control family (4 of 5)" wording, and
> do not carry forward the stale caveat that one response is nominally non-zero.
> `python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); print(d[(d.response=='neg_probe_rate')&(d.design=='n6n5')][['clustered_mean','clustered_lo','clustered_hi','clustered_p']])"`
> The R5 finding below is retained as the record of what was true of the pre-C6 file.

**Status at audit time (pre-C6 file): REFUTED, and self-contradicted two paragraphs later.** `neg_probe_rate` under
`n6n5` is **+0.0108 [+0.0021, +0.0195], p = 0.0204** — CI excludes zero. The report
discloses this immediately afterwards ("Reported against interest: it is the one response
of five whose conditioned amplitude is nominally non-zero"). The word "every" in statement 2
should be "every count-based control family (4 of 5)". Verified in `a7_summary.csv`.

### R6 — "the reportable-fit filter admits two to three times more fits than its nominal rate implies"
**Where:** `CS_PHASE8_CALLERS.md` §4.3 warning. **Also live in the frozen
`reports/PREREG_PHASE8.md` §3.6** — see the dated correction block at the head of that file.

> **⚠ DIGITS SUPERSEDED 2026-08-27 (record reconciliation).** R6's **verdict stands** — the
> claim is not supported — but the table below is computed from the **pre-C6** fits file. The
> frozen `results/phase3/a7_control_probe_fits.csv` gives, for the reportable filter as defined
> (naive, one-sided): `all_controls` **13.3 %**, `neg_control_codeword` **13.3 %**,
> `neg_control_probe` **8.5 %**, `neg_probe_rate` **6.7 %**, `genomic_control` **3.0 %** —
> i.e. **3.0–13.3 %**, not 1.2–13.3 %. On the full `n6n5` design it is **4.8 %, identical
> across all five responses — essentially nominal**, not 3.0–6.7 %. Quote the frozen figures.

**Status: REFUTED by the same file.** The 9–16 % figure is a **two-sided** 95 %-CI exclusion
rate under the full `n6n5` design. The reportable-fit filter is
`beta_naive > 0 AND beta_base_lo > 0` (`code/summarize_phase3_c1.py::reportable`) — one-sided,
on the **naive** design. Measured directly on the null responses:

| response | reportable filter as defined (naive, one-sided) | same filter on `n6n5` | two-sided `n6n5` CI-exclusion (the quoted 9–16 %) |
|---|---|---|---|
| `all_controls` | 13.3 % | 4.8 % | 9.1 % |
| `genomic_control` | 1.2 % | 4.2 % | 10.3 % |
| `neg_control_codeword` | 13.3 % | 6.7 % | 10.9 % |
| `neg_control_probe` | 7.9 % | 3.0 % | 12.7 % |
| `neg_probe_rate` | 6.7 % | 5.5 % | 16.4 % |

`python3 -c "import pandas as pd;f=pd.read_csv('results/phase3/a7_control_probe_fits.csv');print(f.groupby('response').apply(lambda g:((g.beta_naive>0)&(g.beta_base_lo>0)).mean()))"`

Neither column is a uniform 2–3×. On the full design the reportable filter runs at
**3.0–6.7 %, i.e. essentially nominal**; on the naive design it spans 1.2–13.3 %. The
"2–3× nominal" bound on the *estimator* (R-item aside, that part is fine) does not
transfer to the *filter*, and the sentence should be deleted or replaced with the
measured 3.0–6.7 %.

### R7 — "`genesets/mouse_c6/` — Frozen, NOT YET PROMOTED" / "`genesets/*.txt` — Untouched"
**Where:** `PREREG_PHASE8_genesets.md` §1 table and deviation D15.

**Status: REFUTED — stale since 05:41.** PI decision D5 promoted the C6 sets. All 15
mouse Tier A/B files in `genesets/` are now byte-identical to `genesets/mouse_c6/`, and
three of them are modified against `HEAD`:

```
git diff --stat pre-c6-genesets -- genesets/
  A_SENDER_FINAL_strict.txt              25 -> 33
  A_sender_for_secondary_senescence.txt  55 -> 74
  B_secondary_senescence.txt             38 -> 108
```

The promotion itself is exactly as D5 describes (3 files, those three size changes, tag
`pre-c6-genesets` at `e789372`, 05:41:13) — **CONFIRMED**. It is the pre-registration
component file that is now wrong, and it is the document the PI is about to commit.

### R8 — "`code/build_genesets.py` cannot be re-run as committed … Not fixed"
**Where:** `PHASE8_ROADMAP_STATUS.md`, "Known risks, live".

**Status: REFUTED — the risk is closed in the working tree.** `HEAD`'s copy indeed carries
a dead path (`SCRATCH = '/tmp/claude-0/-workspace/e6da7884…/msigdb_mouse'`), but the
working-tree copy (uncommitted, `M code/build_genesets.py`) defaults to
`/workspace/genesets/msigdb_mouse_2026.1.Mm` and raises `SystemExit` with an explicit
message if it globs zero JSONs. `git show HEAD:code/build_genesets.py | grep -n SCRATCH`
vs `grep -n SCRATCH code/build_genesets.py`.

### R9 — `CS_PHASE7_C1.md` §1 and §3 no longer match `null_destructiveness.csv`
**Status: REFUTED against the current file.** The diagnostic was regenerated at
**05:58:41** (`logs/m1_nulldiag.log`), *after* both C1 reports, under the **promoted C6
sender definition**. Section-geometry columns are unchanged (`n_cells`, `occ_frac_bbox`,
`tissue_frac_bbox` all match to 4 dp); every sender-dependent column moved:

| quantity | report | current CSV |
|---|---|---|
| N3-orig median displacement | 2,974 µm | 2,910 µm |
| N3-tile median displacement | 489 µm | 479 µm |
| N4-orig retention | 0.917 | 0.920 |
| N3-occ15 retention / displacement | 0.966 / 304 µm | 0.969 / 317 µm |
| N4-swap displacement | 3,038 µm | 2,980 µm |
| §3: 7259 admissible N3 offsets | 63 of 38,080 | **66** of 38,080 |
| §3: 7352 admissible N3 offsets | 9 of 74,178 | **10** of 74,178 |
| §3: 7259 admissible N4 angles | 12 of 720 | **13** of 720 |
| §1: real median neighbours | 140.5 | 140.0 |

The **conclusions** are untouched — the geometric argument (occupancy-screened offsets are
near-identity; 7001 admits only the identity, 1 of 108,375, displacement 0) holds in both
versions. But the numbers in the two C1 reports and in the proposed §17 footnote are now
provably not the numbers in the file they cite, and `figure2e_data.csv` carries the older
set. This is a mid-freeze citation hazard, not a scientific error.

---

# 2. TRUE BUT MISLEADING

These are arithmetically correct statements supporting an overreaching or unsupported
conclusion. This is the class an external reviewer will find.

**M1 — "24/35 = 69 %", the mouse CoreScence anchor, is a typed-in literal and is not
reproducible.** `PREREG_PHASE8_genesets.md` §7b promises the gene-set figures are built
"from the CSVs above — nothing is recomputed or typed in." `code/make_figure_genesets.py`
L245 reads `bars = [("M1 mouse\npublished", 24, 35, PAL.MUTED), …]`, and
`code/gate_disjointness_human.py` L291 writes `mouse_arm_reference='24/35 = 69%'` as a
string. No script derives it. My independent re-derivation from
`DeepScence/data/coreGS_v2.csv` (occ ≥ 5, n = 39), the MGI pin and the mouse panel gives
**24/31 = 77.4 %** (strict MGI mapping) or **26/33 = 78.8 %** (with the project's own
title-case fallback for `Cdkn2b`/`Cxcl1`) — the numerator agrees, the denominator does not,
and the percentage is ~9 points higher. The human halves reproduce exactly (25/33 = 75.8 %,
29/33 = 87.9 %). Under the *promoted* C6 mouse sets I get 28/31 = 90.3 % (or 29/33 = 87.9 %).
**The "69 % → 76 % → 88 %" narrative rests on an unsourced anchor**; if the mouse baseline
is really ~78 %, the "real cost of C6" is roughly half what Figure gs3's leftmost bar
implies. `figures/figure_gs3_corescence_circularity_data.csv` row 1 exists but is
populated from the same literal.

**M2 — "1:1 human ortholog" and "2,425 of those land on the human panel."**
`PREREG_PHASE8_genesets.md` §6. The map is **not** 1:1: 18,782 rows carry only 17,609
distinct human symbols, and 431 human symbols receive more than one mouse gene. And
**2,435** mouse panel genes have an ortholog on the human panel; 2,425 is the count of
distinct *human* symbols. Both counts are defensible; the sentence conflates them.
(5,097 / 4,845 / 252 / 2,668 / 2,425 all reproduce exactly — see C8.)

**M3 — A7's "five control responses" are four overlapping views of one quantity.**
`all_controls` is the sum of `neg_control_probe`, `neg_control_codeword` and
`genomic_control`; `neg_probe_rate` is a ratio of two of them. Ranges quoted as spanning
"all five responses" ("sign +: 0.34–0.51", "9–16 %") are ranges over correlated statistics,
not five replications. In particular the **16 % upper end of the false-positive rate comes
entirely from `neg_probe_rate`** — the one response the report itself says is not a clean
null, because its denominator is an N5 column. The clean-null subset gives **9.1–12.7 %**.

**M4 — the binned-curve sentence invites a category error.** §4.1: "The binned curve rises
from +0.019 SD in the 0–5 µm bin to +0.29 SD in the 95–100 µm bin. **This gradient is a
quarter of the naive biological gradient (+0.291 SD)**." The 0.019→0.293 span is ~0.27 SD;
the "quarter" refers to 0.070/0.291. Placing a bin-range and an amplitude in adjacent
sentences with near-identical numerals (0.29 vs 0.291) reads as though the control curve
spans a quarter of the biological one when it spans essentially all of it.

**M5 — `figure_gs2`'s headline numbers are not in its `_data.csv`.**
`PHASE8_ROADMAP_STATUS.md` figure policy: "Every figure writes a `*_data.csv` alongside so
every plotted number is auditable." The dark inner segments (26 of 33, 85 of 108/116) and
the entire footnote block are drawn from `crossarm_geneset_table.json`'s `asymmetry` array;
`figures/figure_gs2_crossarm_symmetry_data.csv` contains only
`orthologue_intersected_{mouse,human}_C6` (27 / 26 and 85 / 88), which are a *different*
quantity that coincides numerically on one cell. A reader auditing the CSV would conclude
the figure plots 27, not 26.

**M6 — `results/phase3/figure2c_data.csv` describes a figure that is not committed.**
The on-disk CSV has 19 rows (10 published nulls + 9 corrected variants); the committed
`figure2c.png` plots the 10. `git show HEAD:results/phase3/figure2c_data.csv | wc -l` → 11
lines. Under the stated policy that the CSV is the audit trail for the figure, the pair is
currently inconsistent.

**M7 — the figure guard covers 27 of 45 figure artefacts.** `code/check_figures_guard.py`
enumerates via `git ls-files figures/`, so the 18 **untracked** outputs — `figure2e.{png,pdf}`,
all of `figure_gs1`–`gs4`, `figure_phase8_callers`, `figure_phase8_d3` and their `_data.csv`s
— are outside it. Those are precisely the new Phase 8 figures. The manifest itself
(`figures/.committed_manifest.json`) is gitignored and re-baselineable with `--snapshot`,
so a drifted figure can be blessed rather than caught. **The guard runs clean:
`python3 code/check_figures_guard.py` → `OK: all 27 committed figures match`, exit 0**, and
I independently confirmed all 27 against `git show HEAD:` — but "all figures verified" is
not what it establishes.

**M8 — none of the evidence for C1, A7 or the gene-set freeze is under version control.**
`results/` has 1,238 files on disk and 266 tracked. Every file this audit was asked to check
is untracked: `null_destructiveness.csv`, `sf_summary_c1*.csv`, `perm_nulls_c1*.csv`,
`a7_*.csv/txt`, and the whole of `genesets/human/` and `genesets/mouse_c6/` including
`FROZEN_MANIFEST.csv`. That is how `null_destructiveness.csv` could be rewritten at 05:58
with nothing to diff against (R9), and it means the "freeze" currently protects files that
`git` does not track.

**M9 — A7 and C1 are computed on the *pre-C6* sender definition.** `data/processed/senders_*.csv`
and `modules_*.csv` were regenerated 05:47–05:53 and the Phase 3 caches at 05:54; A7 ran
05:04–05:19 and the audited `perm_nulls_c1.csv` at 04:29. Since `tierA_p95` is defined by
the Tier A score, and Tier A changed 25 → 33 genes at 05:41, **both results are keyed to a
sender set the frozen configuration no longer uses.** `CS_PHASE8_C1_CLOSEOUT.md` §5.1 says
so for C1; nothing says so for A7, and `PHASE8_ROADMAP_STATUS.md` marks 8.5b **DONE** with
no re-run listed while quoting −0.070 into the submission patch.

*Provisional, mid-run:* recomputing C1 against the 06:25 `perm_nulls_c1.csv` and the 06:13
`main_fits.csv` (153 reportable fits, down from 160) gives N3-orig 1.001, **N3-tile 0.971**,
N3-occ 0.302, N3-occ15 0.940, N3-swap 0.695, N3-snap 0.993, N4-tile 0.924, N4-occ 0.183.
**The C1 conclusion survives C6**, but every quoted decimal will move.

**M10 — H1 "Per-sample cell counts (220k–396k) bracket the mouse range (84k–238k)."**
Both ranges are exactly right (mouse raw: 83,621–237,982). They do not bracket; the human
range sits above the mouse one and overlaps only at the top.

---

# 3. CONFIRMED

Compact; each row was re-derived, not read.

### C1 — Gene sets and the disjointness gate — **all CONFIRMED**
Re-derived from `genesets/` and the panels, ignoring every log and report:

| Assertion | mouse (`genesets/`, promoted) | human (`genesets/human/`) |
|---|---|---|
| strict Tier A on panel | **33**, 0 off-panel | **33**, 0 off-panel |
| `len(B_k ∩ panel) ≥ 30`, k = 1…7 | 126/68/100/190/125/**31**/108 | 120/71/126/231/113/**36**/116 |
| `A_strict ∩ B_k == 0` | **0 × 7** | **0 × 7** |
| per-module `A_mod` sizes | 70/74/73/37/73/71/74 | 77/81/80/36/79/79/81 |
| per-module `A_mod ∩ B_k == 0` | **0 × 7** | **0 × 7** |
| **Gate** | **PASS** | **PASS** |

Matches `PREREG_PHASE8_genesets.md` §2.3, §3, §3.1 and §8 to the digit.
`genesets/*.txt` are byte-identical to `genesets/mouse_c6/*.txt` (15/15).
`genesets/human/FROZEN_MANIFEST.csv`: 43 rows, **35 FROZEN / 8 variants**, and I recomputed
all 43 SHA-256 prefixes and gene counts — **0 mismatches, 0 missing**.
MSigDB pins: human 27 JSON, **0 invalid**; mouse 26 JSON, **3 invalid HTML**
(`FRIDMAN_SENESCENCE_UP/DN`, `WP_NRF2_PATHWAY`) exactly as declared.
Variant sizes confirmed: `A_PHASE7_S10_16` 14, `A_S10_16_mouse` 13,
`B_secondary_senescence_v1_curated_ported` 35, `…_v1_curated` 38.
Tier E confirmed: 40 / 609 / 21 / 8 / 1. Tier C: 13 ligands / 15 receptors (human).

### C2 — Mouse panel arithmetic — **CONFIRMED exactly**
```
h5 Gene Expression features        5,106   (identical name set across all 11 sections)
stock CSV 5,006  ∪  add-on 100   =  5,106   overlap = 0, union == h5 list, 0 discrepancies either way
non-ENSMUSG "Gene Expression" ids       9   Brca1×2, Jak2×2, Kras×3, Pkd1×2
AUTHORITATIVE                       5,097   add-on contributes 91 real genes
```
Per-module add-on contribution table (§3.2) reproduces **row for row**: B1 +7, B2 +3, B3 +3,
B4 +1, B5 +10, **B6 +1**, B7 +4. **`B_oxidative_stress` is 31 on the 5,097 panel and 30 on
the stock CSV alone, and the single gene is `Junb`** — CONFIRMED
(`sorted(B6 - stock_csv) == ['Junb']`, and `Junb` is in the add-on file).
Also confirmed that `annotate_pipeline.py:101` and `build_random_null_sets.py:22` do select
`feature_type=='Gene Expression' & id.startswith('ENSMUSG')`, and that
`build_genesets_mouse_c6.py` asserts 5106 / 9 / 5097 and `B7v1==38 and A0==74`.

### C3 — Human panel — **CONFIRMED, and on all 7 samples, not 3**
Every `.h5` in `data/raw_h1/` gives Gene Expression **5,093**, Negative Control Codeword
**609**, Negative Control Probe **40**, Genomic Control **21**, plus 695 Unassigned and
3,291 Deprecated. The Gene Expression name set is identical across all seven (the
PROVENANCE note claims a 3-way check; the 7-way check also passes), and equals
`genesets/h1_candidate/GSE326743_gene_panel_5093.csv` exactly.

### C4 — Cross-arm counts — **CONFIRMED** (except R2's `human_only` list and M2's wording)
5,097 mouse panel → 4,845 with an MGI ortholog → **2,425 distinct human symbols on the
5,093 panel**; 252 mouse genes unmapped; 2,668 human genes unreachable. Tier A: 31 of 33
mouse members map, **26 overlap**. B7: mouse 108 / human 116, 104 map, **85 overlap**, and
on the intersected panel mouse 85 / human 88 — all matching
`crossarm_geneset_table.{csv,json}` and `figure_gs2…_data.csv`.
§8's Tier-B cross-talk headline numbers on the human arm (B7∩B1 = 26, B1∩B3 = 23,
B7∩B5 = 20, B2∩B3 = 22) and §4's "B7∩B1 rises 10 → 26, B7∩B5 6 → 20" all reproduce.

### C5 — C1 null corrections — **CONFIRMED against `sf_summary_c1*.csv`**
Every number in `CS_PHASE7_C1.md` §0 and §4 and in `CS_PHASE8_C1_CLOSEOUT.md` §1 and §2
reproduces. I rebuilt the reportable population from scratch
(`git show HEAD:results/phase3/main_fits.csv`, `call==tierA_p95 & stratum=='all' &
in-band & beta_naive>0 & beta_base_lo>0` → **160 fits**, tile scope **144**) and recomputed
all **66 cells** of the §4 per-section table from `perm_nulls_c1.csv` @ `93194e2e0ce2`:
**0 mismatches.** Rejection rates 0.875 (N3-orig) and 0.826 (N3-tile) confirmed.

- **N3-tile 0.9745 against a published N3 1.00011** — CONFIRMED.
- **N3-occ is degenerate** — CONFIRMED. Admissible offsets 1–66 of 38,080–108,375;
  median displacement 27 µm (N3) / 25 µm (N4), below the 100 µm window and comparable to
  median λ̂ = 12.83 µm.
- **Section 7001: the only admissible offset is the identity, SF = −0.000** — CONFIRMED
  three ways: `n_admissible_moves == 1`, `median_displacement_um == 0.0`, and
  `np.allclose(N3_occ_null_mean, beta_obs) == True`, giving a median `N3_occ_sf` of
  −0.000000. Same for N4 (1 of 720, 0 µm, 0.000).
- **N3-swap ≡ N1** — CONFIRMED. `sf_summary_c1_swap_vs_n1.csv`: N1 0.716, N3-swap 0.721,
  **ρ = 0.948**, median |Δ| = 0.0087 for `tierA_p95`; and the covariate-adjusted
  `full_sf` moves N3-swap 0.721 → **0.999**.
- **The SenePy caveat is real and correctly reported against interest** — CONFIRMED:
  `senepy_p95` ρ = **0.434**, `senepy_p99` ρ = **0.511**, against 0.922 / 0.948 / 0.983 for
  the three Tier A calls. The N7 table (11 variants × 6 calls) reproduces cell-for-cell
  from `sf_summary_c1_n7.csv`, including the reportable-fit counts 198/160/91/203/155/107.
- The three pinned md5s at the top of the closeout: `sf_summary.csv` (`69e3a1d3f600…`) and
  `summary_phase3.txt` (`ecf86b9ca546…`) still hold; **`perm_nulls.csv` no longer does**
  (the M1 re-run rewrote it at 06:13 — expected, DEFERRED).

### C6 — A7 statistics — **CONFIRMED, every cell**
I recomputed all **30 rows × 4 statistics** of `a7_summary.csv` from
`a7_control_probe_fits.csv` — section-clustered mean, one-sample t p-value on the 11
section means, `frac_CI_excludes_zero` and `frac_positive`: **0 mismatches**. Design
counts confirmed: 165 control fits per response × 5 = **825**, and **1,155** module fits.
The three headline numbers are exactly as quoted (subject to R1's naming):
−0.0697 p = 0.023 (base), +0.0068 p = 0.411 (n6n5), −0.0611 p = 0.0204 (n2).
CI half-widths 0.1363 / 0.1366 per fit and ±0.0177 pooled — confirmed. Sparsity
0.0067 / 0.0428 / 0.0094 counts per cell and 0.65 / 3.90 / 0.90 % non-zero — confirmed as
the unweighted mean over the 11 sections in `a7_control_probe_provenance.csv`.
Contamination ratios 0.070/0.291 = 24 % and 0.061/0.281 = 22 % — confirmed.

**The provenance claim is the strongest thing in the report, and it holds.** I re-derived
it from scratch on 7250: summing the h5's 40 / 609 / 21 / 5,106 features per cell and
comparing to `cells.parquet` gives **539 / 3,481 / 1,035 / 147,922,404** with
`np.array_equal == True` for all four, and `cell_id` order-identical.

### C7 — H1 acquisition integrity — **CONFIRMED, and I checked more than was claimed**
All 7 samples: `cell_id` sets in `.h5` and `cells.parquet` match exactly **and in the same
order**; centroids in microns (x 5–6,364, y 7–6,512 across the series, i.e. "~5–6,500 µm");
all six Tier D nuisance columns present; `annotations.csv.gz` carries `Level_1`–`Level_4`;
SPLN07 annotations cover 239,167 of 249,420 cells as stated. **28 files present.**
Per-sample cell counts match the §3 table exactly and total **2,207,593**. The mouse
comparator "11 sections, 1,834,806 cells" is also exact (summed from the 11 `.h5`
barcode axes). No expression value was read from `data/raw_h1/`.

### C8 — Figures — **six chains fully traced, 0 mismatches**

| Figure | `_data.csv` | Verification | Result |
|---|---|---|---|
| `figure1` | `figure1_data.csv` | regrouped `results/sweep_all.csv` (`sweep=='main'`) by `(clustering, ell_over_lambda)`, mean/std, round 4 | 20 rows × **66 columns**, 0 mismatches |
| `figure2a` | `figure2a_amplitudes.csv` | re-ran `_amp()` and the panel loop over `figure2a_stratified_curves.csv` | all **10** amplitudes exact |
| `figure2c` | `figure2c_data.csv` | rows 0–9 ≡ `sf_summary.csv` PRIMARY; rows 10–18 ≡ `sf_summary_c1.csv` | exact (but see M6) |
| `figure2e` | `figure2e_data.csv` | 11 null rows ≡ `sf_summary_c1.csv`; 10 A7 rows ≡ `a7_summary.csv`; **20 curve bins** recomputed n-weighted from `a7_control_probe_curves.csv` | 0 mismatches (destructiveness columns stale per R9) |
| `figure4` | `figure4_data.csv` | 96 panel-a and 20 panel-b values vs `results/phase4/headline.csv`; 2 panel-c values vs `sf_summary.csv` | **118/118 exact** |
| `figure_gs1` | `…_data.csv` | recomputed **196** intersection cells (A×B and B×B, both arms) directly from the `.txt` files | 0 mismatches; every gate cell 0 |
| `figure_gs2` | `…_data.csv` | intersected-panel counts recomputed from the MGI pin + both panels | exact (27/26, 85/88); see M5 |
| `figure_gs3` | `…_data.csv` | human rows ≡ `gate_result_human.json` (29/33 = 87.9 %, 25/33 = 75.8 %, per-module 18/10/9/7/6/4/0) | human exact; **mouse bar unsourced — M1** |
| `figure_gs4` | `…_data.csv` | joined to `senepy_spleen_coverage.csv` on all 5 shared columns | identical; **22 labels = 0 matched / 15 surrogate / 7 none** confirmed, and the hub collapse (one blood memory-B hub for all 3 B labels, one lung T-cell hub for both T subsets) |

**Guard:** `python3 code/check_figures_guard.py` → `OK: all 27 committed figures match
(PDF date stamps ignored)`, **exit 0**. I independently confirmed all 27 against
`git show HEAD:` as well. Scope caveat in M7.

**Figures with no `_data.csv` at all** (untraceable to a file by construction, though most
predate Phase 8): `fig_phase3_caller_depth`, `fig_phase3_composition`,
`fig_phase3_tierC_identifiability`, `figure2b`, `figure2d`, `figure3`,
`figure4_supp_commot_mechanism`, `figure4_supp_ncem_lengthscale`. `figure2b` is the one
that matters — it is a Phase 8 revision target and its null bands have no audit trail.

---

# 4. DEFERRED (in flight)

Per the brief, and confirmed to be moving:

- **All caller-agreement numbers** — the 8.4 gate table (0.932–1.221 / 1.03× / p = 0.20 vs
  0.700–1.711 / 1.118× / p = 1.4e-30), the 11-of-11 Tier A vs DeepScence claim, the
  SenePy-vs-DeepScence flip, `CS_PHASE8_CALLERS.md` §1–§3, `SUBMISSION_PATCH` §1–§3.
  `caller_coverage_gate{,_headline}.csv` moved at 06:09 **and again at 06:20**, and eight
  `*_2sec_c6.csv` tables appeared at 06:19. Not audited.
- `figures/figure_phase8_callers{,_data}` and `figure_phase8_d3{,_data}` — downstream of
  the above. Not audited.
- Anything downstream of `results/phase3/{main_fits,perm_nulls,perm_curves}.csv`,
  `results/phase5/*` and `data/processed/senders_*|modules_*` — all rewritten today.
  The reportable population has already moved 160 → 153.
- `reports/{CORRECTIONS,PREREG_PHASE8,COMPLETED_TASKS,NOVELTY_ASSESSMENT}.md` — rewritten
  during the audit window.
- `denoise=False` / DCA (8.5) — running; `logs/deepscence_d2_*` active at 06:25.

---

# 5. UNVERIFIABLE

- **M1's "24/35 = 69 %"** — no producing script exists; see M1 above. Classed unverifiable
  rather than refuted because a mapping convention I have not reconstructed might yield 35.
- **The GEO screen** (132 series on GPL33762, 19 with Prime-5K evidence, `GSE335761`'s
  386-target panel) — requires network access. The *outcome* of the screen is verified on
  the acquired data (C3, C7).
- **SenMayo / REACTOME_SASP source-set sizes** (mouse MM16098 n = 117 vs human M45803
  n = 124; MM14900 n = 40 vs M27187 n = 111) — these live in the pinned JSONs; I verified
  the pins' validity and file counts but did not re-derive the per-set memberships, since
  the downstream B7 sizes (108 / 116, 85 overlap) already reproduce exactly.

---

# 6. What I would fix before the Aug 29 submission

Ordered by exposure.

1. **R1.** `SUBMISSION_PATCH_2026-08-29.md` L112–113 sends "−0.070 SD … in negative-control
   probes" to the manuscript. Change to "in pooled negative-control features (40 probes,
   609 codewords, 21 genomic controls)", and add that the 40 probes alone are flat
   (−0.018, p = 0.18). A reviewer who requests the per-family breakdown — and the
   pre-registration promises one — finds this immediately.
2. **R3.** "23 % in the void" appears in a proposed Figure 2 caption. Either quote 35 %
   (the occupancy measure) or restate as "23 % lost every real neighbour within 100 µm".
3. **R2.** Remove `CDKN2B` from `human_only` in `crossarm_geneset_table.json` — the figure
   caption is generated from that field — and complete `mouse_only` to 7 entries. State the
   Tier A overlap as 26 by the pinned map / 27 after the documented map-gap correction.
4. **R6, R5.** Two sentences in `CS_PHASE8_CALLERS.md` §4.1 and §4.3 that the same file
   contradicts.
5. **M1.** Either derive the mouse CoreScence baseline or drop the 69 % bar from
   `figure_gs3` and say the mouse arm's value was not recomputed under C6.
6. **R7.** `PREREG_PHASE8_genesets.md` §1 and D15 still say the mouse sets are unpromoted.
   This is the file the PI is about to commit as the pre-registration.
7. **M9.** Decide whether A7 re-runs on the C6 senders before the freeze. `8.5b` is marked
   DONE and its number is already in the submission patch.
8. **M8.** `git add` the Phase 7/8 result files and `genesets/human/`, `genesets/mouse_c6/`
   before tagging `phase8-frozen`. A freeze over untracked files is not a freeze — R9
   happened inside this audit window.
