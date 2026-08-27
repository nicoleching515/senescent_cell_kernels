# Audit — does `PREREG_PHASE8.md` describe the code that actually runs?

**Scope:** every frozen parameter in `reports/PREREG_PHASE8.md` §3 and in the §15 list of
`Phase7_Minimal_Human_Replication (1).md:447`, checked against the code at tag `phase8-frozen`
(commit `9264396`). Read-only audit; no pipeline stage was re-run. Working tree is 12 commits
past the tag (`8f69ca6`); `git diff phase8-frozen HEAD` touches no file named below except
`reports/PREREG_PHASE8.md` itself, so every verdict applies to the frozen tag.

**Method note.** Line numbers below were taken with `cat -n` / `grep -n` on the files as they
stand at the tag. The pre-registration's own citations were *not* trusted, and a large fraction
of them do not resolve — see F-1 and Appendix A.

---

# PART 0 — The finding that conditions every other one

## F-1 (BLOCKING). The frozen tag ships **pre-8.7 code alongside post-8.7 results**. The scripts that produced the frozen numbers were never committed.

Six core Phase-3 files are frozen at their **2026-08-21** state (commit `a6aac3a`), while the
Phase-8 scripts that depend on them were added **2026-08-27** (commit `1351ce8`):

| file | last commit at the tag |
|---|---|
| `code/run_phase3_nulls.py` | `a6aac3a` 2026-08-21 |
| `code/sasp_phase3.py` | `a6aac3a` 2026-08-21 |
| `code/summarize_phase3.py` | `a6aac3a` 2026-08-21 |
| `code/phase2_downstream.py` | `a6aac3a` 2026-08-21 |
| `code/caller_disagree.py` | `a6aac3a` 2026-08-21 |
| `code/phase3_null_geom.py` | `1351ce8` 2026-08-27 |
| `code/summarize_phase3_c1.py` | `1351ce8` 2026-08-27 |
| `code/run_phase3_var.py` | `1351ce8` 2026-08-27 |
| `code/run_phase8_compmatch.py` | `1351ce8` 2026-08-27 |
| `code/caller_disagree_all.py` | `1351ce8` 2026-08-27 |

The M1 re-run (task 8.7) modified `run_phase3_nulls.py` on 08-27. `reports/CS_PHASE8_M1_RERUN.md:124`
itemises those modifications:

> `run_phase3_nulls.py` | `TIERA_PM_CALLS`, `ALL9_CALLS`, `is_permodule()`, `_js()`, `_expand()`;
> `SectionFit(..., module=)`; `_section_job` / `_perm_job` / `_perm_c1_job` take a module and
> restrict to it; `--calls all9|tierA_pm`; `--tag` for output suffixes; `fit_cell` and both perm
> rows record `sender_set`

**None of that is in the committed file.** Every uniform-`15:15` mtime in `code/` is a checkout
artefact (see the "Retire the temporary post-checkout guard" commit); the 08-27 edits to the six
core files were lost and the 08-21 versions were restored, then tagged.

### Six independent proofs, each on its own sufficient

**F-1a — `--stage perm_c1` does not exist.** `code/run_phase3_nulls.py:531-540` dispatches only
`window | main | perm | curves`, else `raise SystemExit(f"unknown stage {a.stage}")`.
`code/_m1_rerun_stage2.sh:11,23` and `code/_m1_driver_b.sh:5` invoke `--stage perm_c1`.
`results/phase3/perm_nulls_c1.csv`, `perm_nulls_c1_n7.csv`, `perm_nulls_c1_pm.csv` exist and are
tracked, `logs/m1_perm_c1.log` shows six `[perm_c1]` job lines — **and no file in the tree writes
them.** `perm_nulls_c1.csv` is *read* by `summarize_phase3_c1.py:133`, `summarize_phase3_var.py:62`,
`m1_headlines.py:60`; it is written by nothing.
→ **The corrected N3/N4 null battery — a headline contribution — has no producer at the tag.**

**F-1b — `--tag` does not exist.** `code/run_phase3_nulls.py:516-521` defines only
`--stage --n-jobs --sections --calls --n-perm`. `code/_m1_rerun_stage2.sh:24,28` pass `--tag _pm`;
argparse would abort with "unrecognized arguments".

**F-1c — `RN._expand` does not exist.** `code/run_phase3_var.py:237` calls
`RN._expand(sections, calls, P.MASTER_SEED, 5000, 17)`. `grep -n "_expand" code/run_phase3_nulls.py`
returns nothing → `AttributeError` on import-time-to-run. **`run_phase3_var.py` is the producer of
N3-var / N4-var, which correction C-5 declares the *primary* corrected pair.** It cannot run at the
tag.

**F-1d — `SectionFit` does not accept `module=`.** `code/run_phase3_nulls.py:118-119`:
`def __init__(self, sample: str, call: str, seed: int, types=None, labels: str = None)`.
`code/run_phase8_compmatch.py:249-251` calls `super().__init__(sample, call, seed, types=types,
labels=labels, module=module)` → `TypeError`. **The composition-matched protocol of §3.8 cannot be
instantiated at the tag.**

**F-1e — `tierApm_*` sender calls are unsupported.** `code/sasp_phase3.py:224-239` handles only
`tierA_p*`, `cdkn1a_pos`, `senepy_p*`; anything else hits `raise ValueError(call)` at `:239`.
`"tierApm_p95".startswith("tierA_p")` is **False**. `git log -S tierApm -- code/sasp_phase3.py`
returns nothing — the branch has never existed in that file. Yet `tierApm_p90/p95/p99` rows are
present in `results/phase3/main_fits.csv` and `window.csv`, and `_m1_rerun_stage2.sh:24,28`,
`_compmatch_chain.sh:12`, `run_phase8_compmatch.py:186,198` all depend on the call.
Correspondingly, `code/phase2_downstream.py` contains **no** `tierApm` flag construction and no
reference to `A_sender_for_*` (`grep` returns nothing), although §3.7 cites
`phase2_downstream.py:107-119` for exactly that.

**F-1f — the result files carry columns the frozen code never writes.**
`results/phase3/main_fits.csv` and `results/phase3/a7_control_probe_fits.csv` both have a
`sender_set` column; `fit_cell` (`run_phase3_nulls.py:169-272`) never emits it.
`results/phase3/window.csv` has **297 rows = 11 sections × 27** (6 non-per-module calls + 3
`tierApm` calls × 7 modules). The committed `stage_window` (`:79-107`) loops sections × `N7_CALLS`
with no module axis and can emit at most **66** rows — and would raise on the first `tierApm` call.

### Consequence

The pre-registration's central promise — "Every value in this section was read from the code or
from a result file on disk, with the file named" (§3 preamble) — does not hold for
`run_phase3_nulls.py`, `sasp_phase3.py`, `summarize_phase3.py` or `phase2_downstream.py`. A
reviewer who clones the tag and tries to reproduce §3.4's corrected nulls, §3.8's composition
protocol, or any `tierApm` sensitivity **cannot run any of them**. This is not a documentation
defect; it is a reproducibility failure of the freeze itself.

**Recommended remedy (out of this audit's scope to apply):** recover the 08-27 working copies of
the four core files, commit them, and re-cut or annotate the tag. If they are unrecoverable, the
freeze must say so explicitly, because the alternative is a tag that certifies code which cannot
produce the results shipped beside it.

---

# PART 1 — MISMATCHES

Ordered by consequence. F-1 above is the first and largest.

## F-2. §3.4 and correction C-5 contradict each other on which corrected null is primary, and the primary one is absent from the frozen battery table.

§3.4 (the frozen N1–N8 table) ends its N3/N4 discussion with:

> **`*_tile` … These are the clustering-preserving in-tissue variants and are the ones to cite as
> the corrected N3/N4.**

Correction C-5 (§0.0), scoped to the §5 *benchmark table*, says in passing:

> **N3-var 0.996 / N4-var 0.985 are now the primary corrected pair**, not the tile variants.

**§3.4 carries no correction marker** (§3.6, §5 and §6 all do). The N3-var/N4-var construction —
the Mrkvicka et al. (2021) variance correction, `code/phase3_null_var.py`, `RADIUS_FRAC = 0.5`
(`:108`), `MIN_N_RETAINED = 50` (`:110`) — **appears nowhere in §3.4, nowhere in §3.5's SF
definitions, and nowhere in the §15 fixed-parameter list.** `phase3_null_geom.ALL_NULLS`
(`:59`) contains no `*_var` member either.

The frozen battery table therefore designates as "corrected" a variant pair the project has since
demoted, and omits the pair it has promoted. Verified values (`results/phase3/sf_summary_var.csv`):
N3_var median **0.995966**, N4_var **0.985068** (both n = 153, n_perm = 1000); N3_tile **0.970552**,
N4_tile **0.924266** (n = 136). C-5's digits are correct; §3.4's designation is stale.

**Verdict: MISMATCHES** — internal contradiction, uncorrected at the site that freezes the battery.

## F-3. §3.6: "blocks are equal-count and none is empty" is false on M1's own data.

§3.6 freezes "**400 replicates over 100 quantile blocks**" and explains:

> a 10 × 10 grid of **marginal quantiles of the x and y cell-centroid coordinates**, computed over
> **all cells of the section**, so blocks are equal-count and none is empty.

`sasp_phase3.block_ids` (**`:259-265`**, not the cited `:290-296`, which is
`FixedLambdaFitter.beta`) takes the marginal deciles of x and of y independently. That makes the
**marginals** equal-count; the 2-D intersections are neither equal-count nor guaranteed non-empty.
Measured directly on the six in-band cached sections (read-only, no fit run):

| section | n cells | populated blocks | min / median / max block |
|---|---|---|---|
| 7259 | 127,386 | **99** / 100 | 0 / 1,267 / 2,678 |
| 7260 | 202,016 | **97** / 100 | 0 / 2,020 / 4,122 |
| 7001 | 165,961 | 100 / 100 | 2 / 1,658 / 2,775 |
| 7248 | 224,921 | 100 / 100 | 775 / 2,203 / 3,492 |
| 7352 | 139,378 | **98** / 100 | 0 / 1,354 / 2,840 |
| 7435 | 172,218 | **99** / 100 | 0 / 1,709 / 2,750 |

Four of six sections have empty blocks **before** any cell-type restriction, and block occupancy
spans 0–4,122. §3.6's "Declared:" note attributes the ≤ 100 effective count to *cell-type-restricted
fits only*; the shortfall is already present at the section level.

What is actually resampled: `run_phase3_nulls.py:190-191` re-indexes to populated blocks
(`ub, bid = np.unique(bid, return_inverse=True)`), and `:241-242` draws
`rng.multinomial(nb, np.full(nb, 1.0/nb))` over that **nb ≤ 100**. So the frozen parameter "100
quantile blocks" is a **nominal grid size, not the number of blocks resampled over**, and it is
recorded as a literal in the outputs (`run_phase8_compmatch.py:433,471`:
`n_blocks=RN.N_BLOCKS_SIDE ** 2`) rather than measured. `compmatch_reruns.csv` accordingly reports
`n_blocks = 100` for every row.

**Verdict:** `N_BOOT = 400` **MATCHES** (`run_phase3_nulls.py:57`, used at `:240-241`);
`N_BLOCKS_SIDE = 10` **MATCHES** (`:56`); "100 quantile blocks" **MISMATCHES** as a description of
what the bootstrap is computed over; "equal-count and none is empty" **MISMATCHES** — false as
stated.

## F-4. §3.11's md5 of the pinned baseline is both wrong and malformed, and the file is not the pinned one.

§3.11: *"Frozen from `results/phase3/summary_phase3.txt` (pinned, md5 `ecf86b9ca5460f31290e2f4c9e822ea2`)"*.

- The quoted digest is **31 hex characters**. No MD5 has 31 characters; it cannot be a valid digest.
- Actual: `md5sum results/phase3/summary_phase3.txt` → **`dc92ddc6605eef52f6359aeab4e16fd7`**.
- §0's provisional table records the file as *"written 08-20 19:03, pre-C6 … the **published
  baseline, deliberately unchanged**"*. Its mtime is **2026-08-27 09:06** and its content is the
  post-C6 8.7 re-run: line 19 reads `fits 315; naive beta > 0 in 216; positive AND
  block-bootstrap CI excludes 0: 153 (49 %)`.
- Consequently §3.5's *"the same 315/**160** appears in the pinned `results/phase3/summary_phase3.txt`"*
  is **false** — the file says 315/**153** — and §5's benchmark table's stated source ("from the
  pinned `summary_phase3.txt`, primary call, **160** reportable fits") no longer resolves to
  anything on disk.

The *substance* is fine: C-5 already carries the corrected 153, and I re-derived 315 fits /
153 reportable / medlam 14.7321 µm / 60 % railed (189/315) / SF IQR [−0.0166, 0.0885, 0.2338]
directly from `main_fits.csv`. What is wrong is the pinning mechanism: **the file the freeze pins by
hash was overwritten, and the hash it pins by is not a hash.**

**Verdict: MISMATCHES.** §3.11's admissibility lists themselves MATCH (below).

## F-5. §3.4 mislabels the null-destructiveness quantity, understating the out-of-tissue fraction by ~12 points.

§3.4: *"`*_orig` — the published bounding-box wrap. It leaves **23 %** of N3's shifted senders and
**8 %** of N4's **outside the tissue** (`results/phase3/null_destructiveness.csv`)."*

Measured from that file, in-band × `tierA_p95` (medians over 6 sections):

| variant | `frac_retaining_a_neighbour` | `frac_in_occupancy` |
|---|---|---|
| N3_orig | 0.7716 → **23 % non-retaining** | 0.6447 → **35.5 % out of tissue** |
| N4_orig | 0.9195 → **8 % non-retaining** | 0.8007 → **20.0 % out of tissue** |

The 23 % / 8 % are `1 − frac_retaining_a_neighbour`, i.e. *shifted senders with no real neighbour
inside the 100 µm window* — a different and weaker statement than "outside the tissue". The true
out-of-tissue fractions are **35.5 %** and **20.0 %**. The mislabel makes the published nulls look
less destructive than they are, in the direction that favours the paper's own correction argument.

**Verdict: MISMATCHES** (numbers correct, quantity misnamed).

## F-6. Criterion R3(c)'s quoted M1 range is not reproducible, and the criterion's variant set is unspecified.

§6 R3(c): *"every in-tissue N3/N4 variant retains **≥ 95 %** of shifted senders with a real neighbour
inside the 100 µm window (M1: **96.6–100 %**), and the published bounding-box variants do not
(M1: 77 % / 92 %)."*

`results/phase3/null_destructiveness.csv`, in-band × `tierA_p95`, `frac_retaining_a_neighbour`:

| variant | min | max |
|---|---|---|
| N3_tile | 0.99867 | 1.00000 |
| N4_tile | 0.99882 | 1.00000 |
| N3_occ | 0.99912 | 1.00000 |
| N4_occ | 0.99855 | 1.00000 |
| N3_swap / N3_snap / N4_swap | 1.00000 | 1.00000 |
| **N3_occ15** | **0.95712** | 0.99566 |
| **N4_occ15** | **0.95446** | 0.99297 |
| N3_orig | 0.72942 | 0.91599 |
| N4_orig | 0.82718 | 0.94332 |

- The **77 % / 92 %** for the bounding-box variants **MATCHES** (medians 0.7716 / 0.9195).
- The **96.6–100 %** does not reproduce under either reading. Including `*_occ15`, the minimum is
  **95.45 %**; excluding it (§3.4 calls `occ15` "supplementary"), the minimum is **99.86 %**.
- §3.4, §3.5 and §6 never define which variants are "in-tissue", so R3(c) is a pre-registered
  replication criterion whose population is not pinned. Under the inclusive reading M1 clears its
  own ≥ 95 % bar by **0.4 percentage points** — a margin the freeze does not disclose.

**Verdict: MISMATCHES** (quoted range unreproducible; criterion population unspecified).

## F-7. §3.7's `TIERA_PM_CALLS` and the "extension to nine" describe code that does not exist.

§3.4's N7 row: *"`run_phase3_nulls.py:68-69` (`N7_CALLS`), extended with `TIERA_PM_CALLS` at `:77-78`"*.
§3.7: *"**six** calls, extended to nine with `TIERA_PM_CALLS` (`:77-78`) under decision D1."*

`grep -rn TIERA_PM_CALLS code/` → **no match anywhere in the tree.** The only occurrence in git
history is in `reports/CS_PHASE8_M1_RERUN.md:124`, describing the uncommitted 8.7 edits (F-1).

`N7_CALLS` at **`:61-62`** (not `:68-69`) is the six-call list and is correct. The three
`tierApm_p90/p95/p99` calls that do appear in `main_fits.csv` reached the pipeline **through the
`--calls` CLI flag** (`code/_m1_rerun_stage2.sh:24,28`), against a `sender_mask` that at the tag
would reject them (F-1e). The "nine-call N7 axis" is therefore frozen in prose only: **nothing in
the code encodes it, and nothing at the tag can run it.**

**Verdict: MISMATCHES.**

## F-8. §3.8's seed-choice justification cites a function that does not exist, and the guarantee it asserts does not hold against the committed code.

§3.8: *"**20260901, 20260902, 20260903, 20260904, 20260905** — chosen outside the range
`run_phase3_nulls._expand` can reach, so no rerun accidentally reproduces an existing N2 match."*

- The five literal seeds **MATCH**: `code/run_phase8_compmatch.py:165`,
  `COMPMATCH_SEEDS = (20260901, 20260902, 20260903, 20260904, 20260905)`; confirmed end-to-end in
  `results/phase3/compmatch_reruns.csv` (`seeds` column = `20260901|…|20260905`, `n_seeds = 5`).
- `run_phase3_nulls._expand` **does not exist** (F-1c), so the justification is unverifiable as
  written.
- Against the seed formulas that *are* committed, the guarantee fails. `stage_perm`
  (`run_phase3_nulls.py:461-462`) uses `P.MASTER_SEED + 5000*i + 17*j` with
  `MASTER_SEED = 20260820` (`sasp_phase3.py:30`). For section index `i = 0` and call index `j = 5`
  — reachable with `--calls all`, since `N7_CALLS` has six entries — the seed is
  `20260820 + 85 = ` **`20260905`**, exactly the fifth compmatch seed, and it is passed straight
  into `SectionFit → match_decoys_section` as the N2 matching seed. The claim "no rerun
  accidentally reproduces an existing N2 match" is therefore not established for the committed
  code. (In the invocations actually used, `_m1_rerun_stage2.sh` passes five N7 calls, so `j ≤ 4`
  and the collision is not realised — but that is an accident of invocation, not the stated
  property.)

**Verdict: MISMATCHES** (seed values MATCH; the justification does not).

## F-9. §11's figure policy makes two claims that are simply untrue.

- *"`code/make_figure2.py` must never be run … **It now refuses to run.**"* — **there is no guard.**
  `grep -n "sys.exit\|raise\|SystemExit\|assert\|refuse" code/make_figure2.py` returns nothing. The
  file is 235 lines ending `if __name__ == "__main__": main()`, and `main()` unconditionally does
  `fig.savefig(f"{FIG}/figure2a.png")` / `.pdf` at `:227-228`. Running it silently overwrites the
  committed `figure2a` — precisely the collision §11 says the policy prevents.
  `git log -S'SystemExit' -- code/make_figure2.py` is empty: the refusal was never written.
  The false claim is repeated at `PREREG_PHASE8.md:848`, `WRITING_PACK.md:1459`,
  `CS_PHASE8_M1_RERUN.md:105`.
- *"a content-hash guard over the **27** committed figures"* — `figures/.committed_manifest.json`
  holds **52 entries** (19 png + 19 pdf + 14 csv) = **19 figures**. `check_figures_guard.py` runs
  clean (52/52, exit 0). `WRITING_PACK.md:1110` already flags the 27 as stale; it was never
  propagated back into the prereg. (PDF date stripping is real:
  `check_figures_guard.py:23,28-29`.)
- *"Every figure writes a `*_data.csv` beside it"* — false for **9 of 19** PNGs
  (`figure2a`, `figure2b/c/d/e`, `figure4_supp_commot_mechanism`, `figure4_supp_ncem_lengthscale`,
  `figure_phase8_callers`, `figure_phase8_d3`). `figure4_supp_commot_mechanism` has **no CSV
  anywhere** (`code/make_figure4_supp2.py` contains no `to_csv`). Six figures' data CSVs live in
  `results/phase3/`, outside the manifest — `grep -c 'results/phase3' figures/.committed_manifest.json`
  = **0**, so those plotted numbers are covered by no hash guard at all.
- **`check_figures_guard.py` is not wired to anything.** It appears in neither
  `.githooks/pre-commit` nor `.claude/settings.json`. It is manual-only, and `--snapshot`
  re-baselines it (the repo has to warn against that in prose,
  `SUBMISSION_PATCH_2026-08-29.md:582`). §11's "the guard makes the policy enforceable rather than
  advisory" is not true as installed.

**Verdict: MISMATCHES** (three claims).

## F-10. §3.9: DeepScence coverage cell count, and a half-asserted D3 anchor condition.

- *"Coverage: **11 / 11** M1 sections, **1.47 M** cells"* — 11/11 **MATCHES**; the cell count
  **MISMATCHES**. The eleven committed `deepscence_*.csv` in `data/processed/` total **1,826,894**
  cells. 1,462,603 is the sum of the **nine newly scored** sections only, omitting the two
  preserved base files (`deepscence_sham.csv` 236,905 + `deepscence_sbr.csv` 127,386 = 364,291).
  The repo's own `CS_PHASE8_M1_RERUN.md:284` and `CORRECTIONS.md:451` already say 1,826,894; §3.9
  inherited the stale figure from `PHASE8_ROADMAP_STATUS.md:84`.
- *"the script asserts this at run time"* (that every anchor gene is **on the mouse panel** and
  absent from every `A_*.txt`/`B_*.txt`) — **only half is asserted.**
  `code/deepscence_reanchor.py:71` asserts the Tier A/B disjointness
  (`assert not bad, 'proliferation anchor overlaps Tier A/B: %s' % bad`, called from `main()` at
  `:126`). **Panel membership is never asserted**; it is silently filtered at `:92`
  (`want = [g for g in genes if g in set(nm[ft == 'Gene Expression'])]`) and `:139`
  (`prolif = np.mean([z(L[g]) for g in PROLIF if g in L], 0)`). If a member fell off the panel the
  anchor would be computed from seven genes with no signal. Both are `assert`s, i.e. void under
  `python -O`.
  The gene list itself **MATCHES** verbatim: `deepscence_reanchor.py:54`,
  `PROLIF = ['Kif20a','Ncaph','Anln','Ect2','Gtse1','Uhrf1','Fen1','Clspn']`.
- *Citation:* the `4,845 of 5,097` rate is **true** (independently re-derived from the h5 feature
  table and the MGI map; 4,845 appears in `logs/ds_smoke.log:1` and every `logs/deepscence_d1/*.log`)
  but is computed by **neither** cited script — `run_deepscence_all.py:58` prints only the mapped
  count; 5,097 comes from `code/deepscence_smoke.py:19`. And **`logs/` is entirely untracked**
  (`git ls-files logs/` → 0 files), so `logs/ds_smoke.log` — the only cited evidence for the
  mapping rate — **does not exist at the tag.**
- *Worth flagging, not a §3.9 mismatch:* "MGI **1:1**" is looser than stated.
  `genesets/mouse_human_orthologs_MGI.csv` has 18,782 rows / 18,782 unique mouse symbols but only
  17,609 unique human symbols — 1:1 on the mouse side only. On the panel the 4,845 mapped probes
  collapse to **4,827** distinct human symbols; 16 human symbols receive 2–3 mouse probes, and
  `var_names_make_unique()` (`run_deepscence_all.py:56`) renames the extras to `SYMBOL-1`, which
  then match nothing in CoreScence. No committed script builds that CSV.

## F-11. §3.10 / §4 give the seven Tier B counts as an unlabelled list in non-canonical module order.

Prereg: *"Tier B: seven modules. M1 **126 / 68 / 100 / 190 / 125 / 31 / 108**; H1 **120 / 71 / 126 /
231 / 113 / 36 / 116**."* The canonical order used by `sasp_phase3.MODULES` (`:42-43`) and by every
result file is `downstream_arrest, emt_ecm, il6_jak_stat3, interferon_response, oxidative_stress,
secondary_senescence, tnfa_nfkb_proximal`. Measured from the files the pipeline reads
(`phase2_downstream.py:11`, `run_phase3_n8.py:44` → `genesets/*.txt`):

| module | mouse | human |
|---|---|---|
| downstream_arrest | **190** | **231** |
| emt_ecm | **125** | **113** |
| il6_jak_stat3 | **68** | **71** |
| interferon_response | **100** | **126** |
| oxidative_stress | **31** | **36** |
| secondary_senescence | **108** | **116** |
| tnfa_nfkb_proximal | **126** | **120** |

The **multisets match exactly** in both arms; the ordering does not. A reader mapping the frozen
list onto `MODULES` order mis-assigns six of the seven counts (e.g. reads `downstream_arrest = 126`
where it is 190). Low-severity but trivially fixable, and it sits in the section that freezes gene
sets.

## F-12. Line-citation drift: every citation into `run_phase3_nulls.py` is wrong; two point past end-of-file.

The file is **541 lines**. §3.2 cites `run_phase3_nulls.py:585` and §3.5 cites `:678` — **both
beyond EOF.** The remaining ~24 citations into that file are off by a consistent +7 (constants) to
+25/+34 (bodies), the signature of the uncommitted 8.7 version (F-1). Full table in Appendix A.

By contrast, citations into the files that *were* committed on 08-27 resolve almost exactly:
`run_phase5_kernels.py` (:51, :52, :53, :54, :79, :122, :124, :169-170, :176, :191-192, :213, :361,
:366) all land on the right lines; `phase3_null_geom.py` (:33, :51-54, :57-61, :193-197, :207-216,
:235-240, :246-247, :249-252, :258-273) all land; `phase3_core.py:146-154` and `:164-173` land
exactly; `summarize_phase3.py:85`, `:99` and `:221` land exactly. This split is itself the strongest
corroboration of F-1.

Smaller citation errors elsewhere (substance correct, line wrong):
`summarize_phase3_c1.py` "200 permutations" is at **`:206`**, not `:191-192`;
`reportable()` is `:50-55`, not `:50-52`; `caller_disagree.py` SCORES is `:36` not `:41`, global
top-5 % `:53-57` not `:60-62`, within-type `:59-63` not `:65-68`;
`phase2_downstream.py` tierA scoring is `:63` not `:73`, flags `:88-95` not `:98-105`;
`sasp_phase3.py` sender-mask branches are `:225-237`, not `:247-268`; `block_ids` is `:259-265`,
not `:290-296`.

---

# PART 2 — MATCHES

Everything below was checked against the code and, where the prereg quotes a number, re-derived
from the tracked result file.

## §3.1 λ grid — MATCHES (with the file's own declarations intact)

| claim | verdict | evidence |
|---|---|---|
| `WINDOW_UM = 100.0` | MATCHES | `run_phase3_nulls.py:52` |
| `LAM_LO_FLOOR = 7.0` | MATCHES | `:53` |
| `N_LAM = 40` | MATCHES | `:55` |
| grid = `exp(linspace(log(7), log(dmax/2), 40))` | MATCHES | `:72` |
| ceiling = window/2 = 50 | MATCHES | `:72`; `main_fits.csv` `lam_grid_lo ≡ 7.0`, `lam_grid_hi ≡ 50.0` on all 7,094 rows |
| 40 log-spaced points | MATCHES | `:72` |
| `med_nn` accepted and never used; floor is a literal, does not adapt per section | MATCHES | `:65` signature, `:72` body — `med_nn` unreferenced |
| "median NN 6.7–10.6 µm; 7 µm at or below it for most but not all" | MATCHES | `poisson_density.csv`: 6.74–10.61 over 11 sections; exactly **1** section below 7.0 |
| railing = `int(t0 == 0 or t0 == sf.lam.size - 1)` | MATCHES (line `:204`, cited `:239`) | `run_phase3_nulls.py:204` |
| three other non-equivalent railing definitions | MATCHES | `run_phase5_kernels.py:124` exact value equality; `:361` `lam[0]*1.001` (0.1 %); `sasp_estimators.py:155` `1.02/0.98` (2 %) |

## §3.2 Fitting window — MATCHES (the correction is accurate)

The known-wrong item is now correctly stated. **`window = 100 µm, fixed` is what the code does**,
and the 99th-percentile rule is correctly recorded as provenance only.

| claim | verdict | evidence |
|---|---|---|
| hard cap on receiver distance | MATCHES (line `:143`, cited `:169`) | `run_phase3_nulls.py:143` `& (self.d_obs <= WINDOW_UM)` |
| same literal at `run_phase3_n8.py:224` | MATCHES exactly | `idx = base & (~smask) & (dd <= RN.WINDOW_UM)` |
| `run_phase3_lamscale.py:52` | MATCHES exactly | `& (d <= RN.WINDOW_UM))` |
| `make_phase5_figs.py:94` | MATCHES exactly | `& np.isfinite(d) & (d <= RN.WINDOW_UM))` |
| `phase4_data.py:36` re-declared literal | MATCHES exactly | `WINDOW_UM = 100.0` |
| `phase3_null_diag.py:31` re-declared literal | MATCHES exactly | `WINDOW_UM = 100.0` |
| COMMOT threshold `phase4_run.py:61` | MATCHES exactly | `COMMOT_DIS_THR = D.WINDOW_UM` |
| spline `dmax` `run_phase5_kernels.py:79, 213` | MATCHES exactly | both `dmax=RN.WINDOW_UM` |
| `:585` (fourth in-file site) | **MISMATCH** — past EOF | file is 541 lines |
| `d_p99 = np.percentile(dr, 99)` in `stage_window` | MATCHES (line `:103`, cited `:128`) | — |
| `window.csv` = 297 rows, 11 sections | MATCHES the file | but unproducible by the frozen code — see F-1f |
| **every percentile in the §3.2 table** | **MATCHES** | `tierA_p95` in-band 76.0–112.1, med(11) 96.1; `cdkn1a_pos` 76.3–160.9, 95.7; `senepy_p95` 118.3–186.8, 151.0; `tierA_p90` 66.3; `tierA_p99` 230.4; `senepy_p99` 354.6 |
| "all six in-band `senepy_p95` p99s exceed the cap" | MATCHES | min 118.3 > 100 |
| truncation shares 7.4–21.5 % / 0.2–2.3 % (senepy) vs 0.8–7.0 % / 0.0–0.11 % (tierA) | MATCHES | `frac_gt_80` 0.074–0.215 & 0.008–0.070; `frac_gt_150` 0.0019–0.0229 & 0.0000–0.0011 |

## §3.3 Kernel families and selection rule — MATCHES in full

| claim | verdict | evidence |
|---|---|---|
| five families | MATCHES | `sasp_kernels.py:81` `FAMILIES = ("exponential","gaussian","powerlaw","step","spline")`; `run_phase5_kernels.py:52` `FAMS` identical |
| forms `exp(-d/λ)`, `exp(-0.5(d/λ)²)`, `(1+d/λ)^-p`, `1[d<λ]` | MATCHES | `sasp_kernels.py:48-78` |
| powerlaw `p ∈ {0.5, 1.0, 2.0, 4.0}` | MATCHES | `run_phase5_kernels.py:53` `P_POW` |
| spline cubic B-spline, `N_KNOTS = 6`, knots at data quantiles, `dmax = WINDOW_UM` | MATCHES | `:54`, `:78-79`; `sasp_kernels.py:64-78` |
| **strict argmin AIC** | **MATCHES** | `run_phase5_kernels.py:176` `best = min(FAMS, key=lambda f: res[f]["aic"])` |
| `aic = n·log(rss/n) + 2(k_par+1)` | MATCHES | `:122` |
| `d_aic_vs_best` recorded per row | MATCHES | `:191-192` |
| "win" rate scored as exact ties only | MATCHES | `summarize_phase5.py:215` `AIC_win=("d_aic_vs_best", lambda v: (v == 0).mean())` |
| **no BIC, no within-section CV, no ΔAIC band** | **MATCHES** | exhaustive grep over `run_phase5_kernels.py` / `summarize_phase5.py`: no ΔAIC tolerance anywhere |
| held-out LOSO Gaussian log-likelihood as second criterion | MATCHES | `run_phase5_kernels.py:203-264`; `phase5_common.py:146-186`; λ on training folds `phase5_common.py:168`; winner `pp.idxmax(axis=1)` `summarize_phase5.py:236` |
| both criteria under `naive` and `ctrl` | MATCHES | `run_phase5_kernels.py:169-170` |
| `N_BOOT = 200` in phase 5 | MATCHES | `:51` |

*Note (not a mismatch):* `sasp_kernels.py:45` defines a **different** 8-value `P_GRID`
(0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0). It is not used by the phase-5 family selection, which
takes `P_POW` from `run_phase5_kernels.py:53`. Worth a one-line declaration so nobody quotes the
wrong grid.

*Tie-breaking:* `min()` returns the first minimum in `FAMS` order, so exact AIC ties resolve to
`exponential`. Not stated in §3.3; harmless in practice but it is a selection rule.

## §3.4 Null battery — MATCHES except F-2 and F-5

| null | verdict | evidence |
|---|---|---|
| **N1** sender-label permutation stratified by cell type, among sender-eligible cells | MATCHES | `permute_within_type` `:324-332` (cited `:379-386`); invoked `:399` (cited `:464-465`); `eligible = ~isin(celltype, EXCLUDE_TYPES + EXCLUDE_FROM_SENDERS)` `:392` |
| **N2** greedy 1-1 propensity match, within section **and** within cell type, caliper **0.25 SD**, contrast = `fit2_shared` | MATCHES | `phase3_core.greedy_ps_match:26-27` (`caliper_sd=0.25`), strata = cell type `phase3_core.py:97-98`, per-section by construction; `run_phase3_nulls.py:216` `fit2_shared`, `:230-232` `beta2_at` |
| **N3 / N4** — *implemented* torus / rotation | MATCHES the *original* only | `torus_shift` `:334-337`, `rotate_about_centroid` `:339-343`. **`phase3_null_geom` is not imported by `run_phase3_nulls.py`** — the corrected variants are unreachable from it (F-1a) |
| corrected-variant tuples, frozen exactly as strings | MATCHES | `phase3_null_geom.py:57-58` — `TRANSLATION` 6 strings, `ROTATION` 5 strings, verbatim |
| `OCC_TOL = 0.05`, `OCC_TOL_RELAXED = 0.15`, `GRID_UM = 25.0`, `TILE_SOLID = 0.98` | MATCHES exactly | `:51-54` |
| `*_orig` at `:193-197`; `*_occ` at `:207-216`; `N3_swap` `:235-240`; `N3_snap` `:246-247`; `N4_swap` = rotate-then-snap `:249-252`; `*_tile` `:258-273` | MATCHES exactly | — |
| `N4_snap` in header prose but not in `ROTATION`; rotation family has five | MATCHES | `phase3_null_geom.py:33` prose; `:58` tuple has 5 |
| **N5** 11 named columns + kNN composition + 2 segmentation dummies | MATCHES exactly | `phase3_core.py:146-154`; column names verbatim |
| **N6** mean over k = 20 NN excluding senders and self | MATCHES | `phase3_core.py:164-174`; `sasp_phase3.py:139-141` `tree.query(coords, k=21)` then `nn_idx[:, 1:]` → 20, self excluded |
| **N7** six calls in `N7_CALLS` | MATCHES (line `:61-62`) | extension to nine → **F-7** |
| **N8** three products; `N_RAND = 200`; expression-matched sets; same cells, same λ | MATCHES | `run_phase3_n8.py:36` `N_RAND = 200`; disjointness `:151-162`; scrambled `:163-204` reading `genesets/E3_random_matched/<module>.tsv` head(200), same `ii`, same `lam`; circularity `:206-240` |
| corrected variants run at **1,000** permutations | MATCHES | `sf_summary_c1.csv`: 13 rows, `n_perm = 1000` on every row |
| declared "200 permutations" header inconsistency | MATCHES (line `:206`) | `summarize_phase3_c1.py:206` |

## §3.5 Surviving fraction — MATCHES

| claim | verdict | evidence |
|---|---|---|
| conditioning nulls: ratio at the **same** λ index `t0` | MATCHES | `run_phase3_nulls.py:235` `out[f"sf_{k}"] = v / b["base"]`, with `t0 = base["t"]` at `:202` and every numerator via `beta_at(..., t0)` `:222-232` |
| perturbation nulls: `(β_obs − mean(β_null)) / β_obs` | MATCHES (line `:431`, cited `:498`) | — |
| "identically at `:678`" | **MISMATCH** — past EOF | file is 541 lines |
| N8 scrambled `sf_n8 = (bstd − br.mean()) / bstd` | MATCHES exactly | `run_phase3_n8.py:203` |
| bootstrap CI on SF is the ratio of **paired** draws | MATCHES | `:263-269` — `sfv = boot[k] / denom` with `denom = boot["base"]`, index-aligned per replicate |
| reportable filter `beta_naive > 0 & beta_base_lo > 0` | MATCHES exactly | `summarize_phase3.py:85`; `summarize_phase3_c1.py:50-55`; `run_phase8_compmatch.py:534` |
| prior floor `MIN_RECEIVERS = 2000` | MATCHES | `:58`, enforced `:182` (cited `:217`), also `:284`, `:370`, `:487` |
| 315 fits at the primary call | MATCHES | re-derived from `main_fits.csv`: in-band × `tierA_p95` × `stratum=="all"` → **315** |
| 160 reportable | superseded → **153** | see F-4; C-5 already carries 153 |

## §3.6 Bootstrap — `400` MATCHES, `100` is nominal (F-3)

| claim | verdict | evidence |
|---|---|---|
| `N_BOOT = 400` | MATCHES | `:57`; loop `:241` `for r in range(N_BOOT)` |
| `N_BLOCKS_SIDE = 10`, `self.nb = 100` | MATCHES | `:56`, `:129-130` |
| resampling by `rng.multinomial` | MATCHES (line `:241-242`) | over **populated** blocks, not 100 |
| CI at 2.5 / 97.5 percentiles | MATCHES (line `:260-261`) | `np.quantile(v, .025/.975)` |
| blocks are marginal x/y quantiles over all cells of the section | MATCHES | `sasp_phase3.py:259-265` |
| "equal-count, none empty" | **MISMATCHES** | F-3 |
| other stages use different counts | MATCHES | `run_phase5_kernels.py:51` 200; `run_phase5_wc.py:57` 200 / `N_SPLIT = 20`; `run_phase5_kernels.py:366` `for _ in range(2000)` |
| **A7 FPR list `0.091 / 0.103 / 0.109 / 0.145 / 0.164`** (correction C-2) | **MATCHES** | `a7_summary.csv`, `design == "n6n5"`, `frac_CI_excludes_zero`: codeword 0.091, all_controls 0.103, genomic 0.109, **neg_control_probe 0.145**, neg_probe_rate 0.164 |
| **C-3's replacement: filter admits 3.0–13.3 % naive, 4.8 % on N6+N5** | **MATCHES** | `a7_control_probe_fits.csv`: naive `(beta_naive>0)&(beta_base_lo>0)` → 0.1333 / 0.0303 / 0.1333 / 0.0848 / 0.0667; n6n5 `(beta_n6n5>0)&(beta_n6n5_lo>0)` → **0.0485 on all five** |
| C-3's withdrawal of the forbidden sentence is in place | MATCHES | correction marker present at `PREREG_PHASE8.md:333-339` |

## §3.7 Sender callers — masks MATCH, `tierApm` does not (F-7)

| claim | verdict | evidence |
|---|---|---|
| three score families enter the fits, DeepScence is not among them | MATCHES | `sasp_phase3.py:225-239` — three branches only |
| `tierA_pNN`: strict `>` NNth pct **within cell type**, per section, types ≥ 20 cells | MATCHES | `sasp_phase3.py:225-226` reads `flag_p{NN}`; flags built `phase2_downstream.py:88-95`, `if m.sum()<20: continue`, `tierA[m] > np.percentile(tierA[m], q_)` |
| `cdkn1a_pos`: `cdkn1a_counts > 0`, no percentile, no stratification | MATCHES | `sasp_phase3.py:227-228` |
| `senepy_pNN`: strict `>` NNth pct within cell type, ≥ 100 finite-score sender-eligible cells | MATCHES | `sasp_phase3.py:229-237`, `if sel.sum() < 100: continue` |
| `tierApm_pNN` | **MISMATCHES** | no branch; `raise ValueError` at `:239` — F-1e, F-7 |
| `tierA_score = score_genes(A_SENDER_FINAL_strict, ctrl_size=200)` | MATCHES (line `:63`, cited `:73`) | `phase2_downstream.py:63` |
| all masks finally `& ok`, excluding `Low_quality/Unknown/unknown` and `Proliferating` | MATCHES | `sasp_phase3.py:33`, `:35`, `:224`, `:240` |
| **four scores** in the caller-agreement analysis | MATCHES | `caller_disagree_all.py:37` exactly; `caller_disagree.py:36` |
| three thresholding rules, all top-5 %, per section | MATCHES | global `caller_disagree.py:53-57`; within-type `:59-63` (strata ≥ 50 at `:62`); depth+type matched `caller_disagree_all.py:84-109` |
| `N_DEPTH_DECILES = 10`, strata ≥ 50 | MATCHES exactly | `caller_disagree_all.py:38-39` |
| depth+type-matched rule is primary | MATCHES as a stated convention | not machine-enforced |
| SenePy on H1: `spleen_hub: false`, 65 hubs, 22 labels, 15 surrogate, 7 no-hub | MATCHES | `results/phase7_jobA/senepy_spleen_coverage.json` |

## §3.8 Composition-matched protocol — every frozen value MATCHES the result file

Seeds, matching set, rule, variants, scope, gate and all fifteen M1 numbers were re-read from
`results/phase3/compmatch_reruns.csv` (`row_type = summary`, `call = tierA_p95`,
`scope_kind = pooled`) and reproduce exactly:

| field | frozen | measured |
|---|---|---|
| five literal seeds | 20260901–05 | `seeds = 20260901\|…\|20260905`, `n_seeds = 5` |
| matching set = 9 `knn_frac_*` columns | 9 named | `n_match_cols = 9`, `matched_on` identical, verbatim |
| caliper 0.25 SD | ✓ | `caliper_sd = 0.25`; `run_phase8_compmatch.py:172` `CALIPER_SD = 0.25` |
| 400 boot / 100 blocks / window 100 | ✓ | `n_boot = 400`, `n_blocks = 100`, `window_um = 100.0` |
| `comp` n = 165, SF **0.9837**, removed **1.6 %**, CI [0.973, 0.994], spread 0.98370–0.98397 (1.2e-4) | ✓ | 165, 0.983735, 0.016265, [0.973037, 0.994166], 0.983695–0.983965, sd 1.20e-4 |
| `full` 165, 0.9855, 1.4 %, [0.979, 0.992], 0.98539–0.98577 (1.6e-4) | ✓ | 165, 0.985495, 0.014493, [0.979163, 0.992347], 0.985385–0.985773, sd 1.64e-4 |
| `comp_adj` 33, 0.4989, 50.1 %, [0.421, 0.606] | ✓ | 33, 0.498892, 0.501108, [0.421492, 0.606284] |
| `type_adj` 33, **0.3415**, 65.9 %, [0.236, 0.402] | ✓ | 33, 0.341450, 0.658550, [0.235641, 0.401643] |
| `typecomp_adj` 33, **0.1461**, 85.4 %, [0.052, 0.246] | ✓ | 33, 0.146139, 0.853861, [0.052005, 0.246070] |
| median match rate 0.999871 | ✓ | 0.999871 |
| balance gate max\|SMD\| 0.0916 → 0.0352, passes 100 % | ✓ | 0.091612 → 0.035222, `frac_smd_gate_pass = 1.0` |
| `SMD_GATE = 0.10` | ✓ | `run_phase8_compmatch.py:173` |
| sections = `sasp_phase3.IN_BAND` (6) | ✓ | `:185`; `n_sections = 6` |
| five variants | ✓ | `:170` `VARIANTS = ("comp","full","comp_adj","type_adj","typecomp_adj")` |
| **H1 gated**: refuses unless sections populated **and** `SASP_H1_UNFROZEN=1` | **MATCHES** | `:191-200` `sections=[]`, `frozen=True`; `:206-215` `arm_config` raises `SystemExit` |
| reportable population identical to `summarize_phase3.sf_table` | MATCHES | `:534` `(d.beta_naive > 0) & (d.beta_base_lo > 0)` |
| pooled is primary; per-type secondary | MATCHES | `:359`, `:465` `scope_kind` |

**Two defects only:** the seed-choice justification (F-8) and the fact that the script cannot run
against the frozen `run_phase3_nulls.py` (F-1d).

## §3.9 DeepScence — MATCHES except F-10

| claim | verdict | evidence |
|---|---|---|
| `denoise=False`, `random_state=0` | MATCHES, cited lines exact | `run_deepscence_all.py:60` and `run_deepscence.py:40`, both literally `api.DeepScence(A, denoise=False, verbose=False, random_state=0)` |
| the two scripts byte-identical in settings | MATCHES | same min-count filter, same ENSMUSG/GeneExpression filter, same ortholog CSV, same 5-dp rounding. (Non-settings difference: `run_deepscence_all.py:24-25` pins `torch.set_num_threads`; `run_deepscence.py` does not — a float-reduction-order nuance, not a settings mismatch) |
| anchor = published `CDKN1A`, not overridden | MATCHES | `DeepScence/io.py:203` hardcodes the `CDKN1A` row; `DeepScence/api.py:101` calls `fix_score_direction` unconditionally; `DeepScence.DeepScence` exposes no anchor parameter and neither run script passes one |
| minimum counts ≥ 20 | MATCHES, cited lines exact | `run_deepscence_all.py:54`, `run_deepscence.py:35` |
| ortholog-remapped, 4,845 of 5,097 | numbers MATCH; citation defective | F-10 |
| coverage 11/11 | MATCHES | 11 committed `deepscence_*.csv` |
| coverage 1.47 M cells | **MISMATCHES** | 1,826,894 — F-10 |
| `denoise=True` runnable from committed code | MATCHES | `setup_dca_env.sh`, `dca_denoise_worker.py`, `run_deepscence_dca.py` (`denoise=True` at `:68`), `_shims_dca_bridge/dca/{__init__,api}.py`, `run_deepscence_denoise_probe.py` — all tracked; `results/phase8_d2/runmeta_dca_7239_*.json` records `denoise: true`, 83,392 × 4,845, 16.4 min |
| D3 anchor gene list | MATCHES verbatim | `deepscence_reanchor.py:54` |
| "the script asserts this at run time" | **PARTIAL → MISMATCHES** | F-10 |

## §3.10 / §4 Gene sets as loaded at runtime — MATCHES

| claim | verdict | evidence |
|---|---|---|
| Tier A PRIMARY `A_SENDER_FINAL_strict`, **n = 33 on both arms** | MATCHES | `genesets/A_SENDER_FINAL_strict.txt` 33; `genesets/human/A_SENDER_FINAL_strict.txt` 33 |
| seven `A_sender_for_<module>` sensitivity sets present | MATCHES | 7 files in `genesets/` |
| seven Tier B modules, both arms | MATCHES (multisets) | counts table in F-11; **order defect** F-11 |
| B6 mouse = **31**, margin one over the 30 floor, and the margin gene is `Junb` | MATCHES | `B_oxidative_stress.txt` = 31 lines; `grep -c Junb` = 1 |
| Tier C ligand/receptor, Tier D nuisance, Tier E controls present | MATCHES | `C_ligands.txt`, `C_receptors.txt`, `D_*`, `E_*`, `E3_random_matched/` |
| promoted C6 sets are what the pipeline reads | MATCHES | all 7 `genesets/B_*.txt` byte-identical (`cmp`) to `genesets/mouse_c6/B_*.txt` |
| `genesets/human/FROZEN_MANIFEST.csv` — **35 frozen, 8 variants** | MATCHES exactly | 43 rows: `Counter({'FROZEN': 35, 'variant (reported, not used)': 8})`; `genesets/human/variants/` holds 8 files |
| `genesets/.geneset_manifest.json` — **96 watched files** | MATCHES exactly | 96 entries |
| loaders read `genesets/*.txt` | MATCHES | `phase2_downstream.py:10-11`; `run_phase3_n8.py:43-44` |

*Note:* `E_negative_control_probes.txt` exists only under `genesets/human/`. §10.1 designates it the
pre-registered **primary technical null** via `PREREG_PHASE8_genesets.md §11`; on the mouse arm the
equivalent response is built from the panel's control features by `run_a7_control_probes.py`, not
from a gene-set file. Not a mismatch, but the asymmetry should be stated where the designation is
made.

## §3.11 Section admissibility — lists MATCH, pin does not (F-4)

| claim | verdict | evidence |
|---|---|---|
| in band, 6: 7259, 7260, 7001, 7248, 7352, 7435 | MATCHES exactly | `sasp_phase3.py:67-69` |
| over ceiling, 4: 7239, 7448, 7361, 7450 | MATCHES exactly | `:70-71` |
| below floor: 7250 | MATCHES exactly | `:72` |
| pinned by md5 | **MISMATCHES** | F-4 |

## Corrections in §0.0 — all seven re-derived and MATCHING

| item | verdict | re-derivation |
|---|---|---|
| **C-1** A7 frozen digits | MATCHES | `a7_summary.csv` `all_controls`: naive −0.0744 [−0.1306, −0.0182] p 0.0145; n2 −0.0642 [−0.1113, −0.0172] p 0.0124 (0.0642/0.0744 = **86 %** undiminished); n5 +0.0038 [−0.0186, +0.0261] p 0.715 |
| **C-2** FPR list | MATCHES | see §3.6 table |
| **C-3** forbidden claim withdrawn; 3.0–13.3 % / 4.8 % | MATCHES | see §3.6 table |
| **C-4** the bracket is an IQR, not a bootstrap interval | **MATCHES exactly** | `summarize_phase3.py:99` `a, b, c = np.quantile(v, [.25, .5, .75])`; `sf_summary.csv` columns are `subset, null, n, q25, median, q75, frac_le_0, frac_gt_05` — **no CI column**; the per-fit `sf_n2n5n6_lo/hi` in `main_fits.csv` are the only bootstrap output and carry no interval on the median across fits. Primary row: n = 153, q25 −0.016643, median 0.088468, q75 0.233783 |
| **C-5** frozen benchmark replacements | MATCHES | re-derived from `main_fits.csv`: amplitude 0.32879, SF 0.088468 IQR [−0.01664, 0.23378], SF N5 0.11502, railed 60.00 % (189/315), reportable 153; `m1_final_audit.txt:31-32` power80 0.1833, SE 0.0654; tile 0.9706/0.9243; N3-var 0.99597 / N4-var 0.98507 |
| **C-6** destructiveness digits | MATCHES | `null_destructiveness.csv`, in-band × `tierA_p95`: N3_occ `n_admissible_moves` 1–66 of `n_candidate_moves` 38,080–108,375; median displacement N3_occ **28.30 µm**, N4_occ **24.69 µm** |
| **C-7** λ̂ = 14.7 µm | MATCHES exactly | `median(lam_naive)` over the 315 primary fits = **14.7321**; printed by `summarize_phase3.py:221` (`dd.lam_naive.median()`) into §6 `medlam`; IQR/railing: 60 % railed confirmed |

## §6 replication criteria — R3(a) MATCHES

`results/phase3/poisson_fits.csv`: "ALL sections × ALL sender definitions" slope **−0.5249**,
r² **0.9843** ≥ 0.95 ✓. (In-band-only subset: slope −0.5292, r² 0.9833 — also clears.)
R3(c) → F-6. R1's corrected reading (IQR includes 0, q75 < 0.50) is satisfied on M1 by
[−0.0166, 0.2338] ✓.

---

# PART 3 — UNVERIFIABLE without a run

1. **Effective block count per fit.** The nominal 100 is contradicted at section level (F-3), but
   the per-fit populated-block count after the cell-type and window restrictions is not recorded in
   any output. Settling it requires instrumenting `fit_cell` and re-running — out of scope. **The
   right fix is to emit `nb` as a column**, not to re-run for this audit.
2. **§3.8's bit-identity claim** — *"`beta_naive` and `lam_naive` are bit-identical to
   `main_fits.csv` on all 142 cells reportable in both (max \|Δβ\| = 0.000e+00)"*. Both files are on
   disk and the merge is checkable, but the claim as stated is about a harness that cannot be
   executed at the tag (F-1d), so re-establishing it means running the recovered 8.7 code. Not
   attempted.
3. **The five per-seed pooled SFs** (0.983695 / 0.983887 / 0.983965 / 0.983714 / 0.983735) are
   bracketed exactly by the summary row's `sf_across_seed_min/max` (0.983695 / 0.983965), so they
   are consistent; reading them individually needs the `row_type = seed` rows, which I did not
   enumerate.
4. **Whether the 08-27 versions of the four core scripts are recoverable.** They are not in git
   history (`git log -S` finds the symbols only in a report). Whether a copy survives in a container
   layer, a `__pycache__`, or a log is a recovery question, not a code-reading one.

---

# PART 4 — Are the §10 prohibitions enforceable?

**Answer: no. All twelve are honour-system statements. There is no mechanism anywhere in the repo
that could catch a violation of any of them.**

Whole-repo sweep of `code/`, `.githooks/`, `.claude/`:

- `grep -rniE 'may never be reported|must not be quoted|prohibit|forbidden'` → **0 hits**
- `CXCL8`, `CXCR1`, `MMP3`, `TIMP1`, `type_adj`, `typecomp_adj`, `jaccard`, `seed_stability` → **0
  occurrences in `code/`**
- the literals `0.9837`, `0.3415`, `0.1461`, `1.51`, `2.85` → **0 occurrences in any script, hook or
  config**
- no `.github/`, no Makefile, no `tox.ini`/`pytest.ini`/`pyproject.toml`, no grep-based text linter
  over `reports/` or `figures/`

| § | Prohibition | Enforced? | Evidence / gap |
|---|---|---|---|
| 10.1 | No naive / N2-only kernel as a distance effect; name the response "pooled negative-control features" | **No** | nothing reads `a7_summary.csv` as a gate |
| 10.2 | No caller-independence claim | **No** | prose only |
| 10.3 | No age-stratified H1 claim | **No** | prose only |
| 10.4 | No MZ-specific confirmatory claim | **No** | prose only |
| 10.5 | No cross-arm species/tissue attribution | **No** | prose only |
| 10.6 | `CXCL8`/`CXCR1` not a mouse replication; `CXCL2`/`CXCL5` a map gap | **No** | symbols absent from `code/`; the disjointness gate checks set sizes and overlaps, never ortholog provenance |
| 10.7 | "1.51–2.85×" must not be quoted | **No** | fixed literal, ungrepped; `caller_coverage_gate.csv` regenerable and unguarded |
| **10.8** | Composition-matched **0.9837 never alone** | **No** | a two-number co-occurrence rule — the cheapest possible check — and it does not exist. Per the prereg's own wording, quoting the first alone "states the opposite of what the data says" |
| **10.9** | `mor` never as evidence normalisation can't move the caller | **No** | no pairing check with the `lib` result; the only `mor` hit in `code/` is the word "memory" in a comment |
| **10.10** | No single-seed `denoise=True` number without its seed-stability companion | **Partially, wrong layer** | `code/_shims/dca/api.py:2` raises `RuntimeError("DCA denoising is not installed…")`, which blocks *computing* denoise in the main env but does nothing to stop an already-computed single-seed number being quoted. `analyse_d2_stability.py` computes the companion but asserts nothing |
| **10.11** | D2 `raw` `rho_signed_dz_vs_depth` (−0.47, −0.16) must not be quoted | **No — actively worse** | `analyse_d2_denoise.py:168` emits the column and `report_d2_tables.py:52,68` tabulates it, **including the `raw` rows**. The pipeline hands the author the forbidden number, pre-formatted, in a report table |
| 10.12 | Four struck sentences | **No** | four exact strings; no search for them anywhere. They survive only in `PREREG_PHASE8.md:833-837` |

### Guards that do exist, and their wiring

| guard | wired? | detail |
|---|---|---|
| `code/gate_genesets_guard.py` | **Yes, properly** | `git config --get core.hooksPath` = `.githooks`; `.githooks/pre-commit` executable, invokes it at `:30`; also a Claude Code `PostToolUse` hook on `Write\|Edit\|NotebookEdit` via `.claude/settings.json:9` → `code/hook_geneset_gate.sh`. Runs clean, exit 0 |
| `code/gate_disjointness_human.py` | Yes, transitively | subprocessed from `gate_genesets_guard.py:157`; own `sys.exit(0 if FROZEN_OK else 1)` at `:340` |
| **`code/check_figures_guard.py`** | **NO — wired to nothing** | absent from `.githooks/pre-commit` and `.claude/settings.json`. Manual-only, and `--snapshot` silently re-baselines it. §11's claim that it "makes the policy enforceable rather than advisory" is false as installed |
| `code/m1_final_audit.py` | **not a gate** | zero `sys.exit`/`assert`/`raise`; prints `SUPERSEDED` for drifted md5s (`:38-41`) and always exits 0 |
| `code/_repro_check.py`, `code/phase3_var_validate.py` | **not gates** | zero exits/asserts; manual, argv-driven |

**Coverage hole in the one always-on hook:** `code/hook_geneset_gate.sh:34` explicitly `exit 0`s on
`*/results/*` and `*/reports/*`. **Every path where a §10 violation would actually be written is
deliberately excluded from the only hook that always fires.**

`reports/WRITING_PACK.md` (1,526 lines) contains many `python3 -c` verification recipes, but they
are snippets for a human to paste: no driver, no exit codes, no record of which passed, and its
`[F]`/`[V]` provenance flags are hand-typed. It also restates the false "it now refuses to run"
at `:1459`.

### Cheapest fixes, in order

1. **§10.11** — highest accident risk, because a committed script re-emits the forbidden values
   into a report table. Either drop the `raw` rows from `report_d2_tables.py` output or mark them.
2. **§10.8 / §10.9 / §10.7 / §10.12** — all four are literal-string or co-occurrence rules over
   `reports/`. One ~40-line linter wired into `.githooks/pre-commit` would cover all of them.
3. **Wire `check_figures_guard.py` into `pre-commit`** and remove `--snapshot` from the default
   path.
4. **Remove the `results/`+`reports/` exclusion** from `hook_geneset_gate.sh:34`, or add a second
   hook that covers those paths.
5. §10.1–10.6 are semantic and genuinely hard to mechanise. They should be **labelled as
   review-only** rather than sitting in a list that implies the guards cover them.

---

# Appendix A — citation resolution table, `code/run_phase3_nulls.py` (541 lines)

| prereg cite | actual | content |
|---|---|---|
| `:59-62` constants | **`:52-58`** | `WINDOW_UM … MIN_RECEIVERS` |
| `:59` `WINDOW_UM` | **`:52`** | `WINDOW_UM = 100.0` |
| `:63-64` `N_BLOCKS_SIDE`/`N_BOOT` | **`:56-57`** | |
| `:65` `MIN_RECEIVERS` | **`:58`** | |
| `:68-69` `N7_CALLS` | **`:61-62`** | |
| `:77-78` `TIERA_PM_CALLS` | **absent** | F-7 |
| `:86-93` `lam_grid` | **`:65-72`** | |
| `:100-136` `stage_window` | **`:79-107`** | |
| `:128` `d_p99` | **`:103`** | |
| `:155-156` blocks built | **`:129-130`** | |
| `:159-160` N2 match | **`:133-134`** | |
| `:169` window cap | **`:143`** | |
| `:217` `MIN_RECEIVERS` enforced | **`:182`** | |
| `:226-227` block re-index | **`:190-191`** | |
| `:236` `sf_{k}` | **`:235`** | |
| `:239` `lam_railed` | **`:204`** | |
| `:251-253`, `:265-267` N2 fits | **`:216`, `:230-232`, `:250-252`** | |
| `:276-277` `rng.multinomial` | **`:241-242`** | |
| `:295-296` CI percentiles | **`:260-261`** | |
| `:297-304` paired SF bootstrap | **`:263-269`** | |
| `:379-386` `permute_within_type` | **`:324-332`** | |
| `:389-399` torus | **`:334-337`** | |
| `:402-408` rotation | **`:339-343`** | |
| `:464-465` N1 invoked | **`:399`** | |
| `:498` perturbation SF | **`:431`** | |
| `:585` window cap (4th site) | **past EOF** | |
| `:678` perturbation SF (2nd site) | **past EOF** | |

---

*Audit performed 2026-08-27. Read-only; this file is the only artefact written. No pipeline stage
was executed; `data/raw_h1/`, `results/phase9_h1/`, `data/processed_h1/` and `code/h1_*` were not
touched.*
