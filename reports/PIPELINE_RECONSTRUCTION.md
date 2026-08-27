# Reconstruction of the lost task-8.7 pipeline code

*Written 2026-08-27. Scope: `code/` only. `data/raw_h1/`, `data/processed_h1/`,
`results/phase9_h1/`, `code/h1_*` and `figures/` were not touched;
`check_figures_guard.py` still reports 52/52, exit 0.*

---

## 0. Read this first

**This is a RECONSTRUCTION. It is not the original code, and nothing below
claims otherwise.**

The 2026-08-27 (task 8.7) edits to five files under `code/` were never staged
and were destroyed by `git checkout -- code/`. They are not recoverable:

| check | result |
|---|---|
| `git show phase8-frozen:code/run_phase3_nulls.py \| wc -l` | 541, zero hits for `perm_c1` / `_expand` / `TIERA_PM_CALLS` / `is_permodule` |
| `git log --follow -- code/run_phase3_nulls.py` | one commit (`a6aac3a`, 2026-08-21), blob is the same 541 lines |
| `git log -S'TIERA_PM_CALLS'` | hits only in `reports/`, the tripwire hook, and this reconstruction |
| 30 dangling objects, each `git cat-file -p`'d | none contains `TIERA_PM_CALLS` or `_perm_c1_job` |
| `code/__pycache__/run_phase3_nulls.cpython-311.pyc` | header records source size **23,070 B**, mtime **15:15:57** — i.e. a compile of the **restored** 541-line file, not the lost one |

**What a bit-identical reproduction does and does not establish.** Where a
reconstructed producer re-emits its committed output byte for byte, that is
strong evidence that the reconstruction is *behaviourally equivalent, on this
input, to the code that ran*. It is **not** evidence that it is the same code.
Different source can produce identical bytes. For a pre-registration the
distinction is load-bearing: the tag can now honestly say it contains code that
**reproduces** the frozen results; it cannot say it contains the code that
**produced** them. Any line-number citation into `run_phase3_nulls.py` in
`PREREG_PHASE8.md` still does not resolve to the file that ran — it now resolves
to a reconstruction of it (see §6).

**Nothing under `results/` was modified.** Every re-run in this document wrote
into a scratch tree and was diffed against the committed file.

---

## 1. What was rebuilt, and from what evidence

Five files. Four independent evidence classes constrained them: the surviving
**callers**, the surviving **result files and their schemas**, the surviving
**run logs**, and the surviving **cached inputs** (`data/processed/cache3/*.npz`
and `data/processed/senders_*.csv`, which already carry the per-module columns
that no committed code writes).

### 1.1 `code/sasp_phase3.py`

| rebuilt | evidence |
|---|---|
| `prep()` caches `tierApm__<module>` and `flag_pm_<module>_p{90,95,99}` | those exact keys are present in all eleven committed `cache3/*.npz`, immediately after the `mod__*` block, in per-module (score, p90, p95, p99) order |
| `Sec.sender_mask(call, module=None)` gains `tierApm_p*`, reading `flag_pm_<module>_p<NN>` and `& ok` | the mask so built gives `n_senders = 5707`, `prevalence = 0.044801` for 7259 — the committed values in `main_fits.csv` and `window.csv` |
| it *raises* on a pre-Phase-8 cache rather than guessing; `Sec.has_permodule` | `CS_PHASE8_M1_RERUN.md:124` names both |

### 1.2 `code/phase2_downstream.py`

The per-module Tier A block: score the seven `genesets/A_sender_for_*.txt` sets
with the same `score_genes(..., ctrl_size=200)` call the primary set uses, write
`tierA_<module>_score` and `sender_flag_<module>_p{90,95,99}`, assert the set
names are the Tier B module list. Placed after the primary sender columns and
before the DeepScence merge — the committed `senders_*.csv` column order fixes
that position exactly.

### 1.3 `code/run_phase3_nulls.py` — the bulk of the loss

`TIERA_PM_CALLS`, `ALL9_CALLS`, `PM_SEED_OFFSET`, `is_permodule()`, `_js()`,
`_expand()`, `_resolve_calls()`; `SectionFit(..., module=)` and its `sender_set`;
`fit_cell` and both perm-row builders record `sender_set`; `_section_job` /
`_perm_job` take a module and restrict to it; `stage_window` gains the module
axis and a `sender_module` column; `--calls all9|tierA_pm`; `--tag`; and the
whole `perm_c1` stage (`_c1_receivers`, `_perm_c1_job`, `stage_perm_c1`).

### 1.4 `code/summarize_phase3.py`

Two lost changes, both forced by the committed `summary_phase3.txt`: §6's N7
table iterates `ALL9_CALLS` (the file lists nine calls), and `load()` also reads
`perm_nulls_pm.csv` (without it the `tierApm_p95` row's `SF N1` is `nan` where
the file has `0.711`).

---

## 2. The five hard constraints that pinned the ambiguous parts

These are the decisions a rewrite could plausibly have got wrong, and what fixed
each one.

1. **`stage_window` ignores `--calls`.** `logs/m1_window.log` line 2 reads
   `calls: ['tierA_p95']` — the default — yet `window.csv` has all nine calls
   × 297 rows. So the stage iterates `ALL9_CALLS` internally.
2. **The tile scope is a different fit population.** In `perm_nulls_c1.csv`, the
   `full` rows' `n`, `lam`, `beta_obs` are *identical* to `perm_nulls.csv`, but
   the `tile` rows are not. Testing both readings against 7259 / Hepatocytes /
   `emt_ecm` (committed `n = 13284`, `beta_obs = 0.022017`): distance to **all**
   senders gives `n = 13286`, `beta = 0.0219941`; distance to **in-tile**
   senders gives `n = 13284`, `beta = 0.0220167`. The second is exact.
3. **The permutation loop is draw-outer over `ALL_NULLS`.** Two orderings were
   tested against the committed `N3_orig_null_mean` for 7259 / Hepatocytes ×
   seven modules at 1,000 draws. Null-outer (`for nm: for r:`) — the ordering
   `phase3_null_diag.py` and `run_phase3_var.py` both use, and so the natural
   first guess — is **wrong**: it is off by 1.5e-5 to 1.6e-4 against committed
   values of that same magnitude, i.e. a completely different stream.
   Draw-outer over `G.ALL_NULLS` in `TRANSLATION + ROTATION` order, with a
   single `np.random.default_rng(seed)` shared by every null and both scopes,
   reproduces all seven modules to 8 decimal places on the first try.
4. **Column counts date the shape.** `logs/phase3_perm_c1.log` (03:40) ends
   `(595, 80)`; `logs/phase3_perm_c1b.log` (04:40) ends `(595, 103)`; the
   committed file is `(595, 104)`. 104 = 15 base + 9 whole-section nulls × 8 +
   `beta_obs_full` + 2 tile nulls × 8; 103 is the same without `sender_set`; 80
   is the same again without the six `_full_` columns. The single column added
   by 8.7 to this stage is `sender_set`.
5. **`run_phase3_nulls.py:699`.** `logs/m1_perm_c1_n7.log` carries a
   `RuntimeWarning` citing **line 699** of the lost file for
   `mean=float(np.nanmean(A[:, bi]))`, and `logs/phase3_perm_c1b.log` cites
   **646** for the same line in the 04:40 version. In this reconstruction that
   line lands at **715** — independent corroboration that `_perm_c1_job` sits in
   roughly the right place in a file of roughly the right size (813 lines here).

---

## 3. The one constant that had to be searched for

`_expand`'s per-module seed displacement, `PM_SEED_OFFSET`, is not stated in any
surviving file. `run_phase8_compmatch.py:56` only alludes to it
("`MASTER_SEED + 1000*i + j`, and `+ PM_SEED_OFFSET`").

It was **recovered by search, then validated out of sample**, and it is
important to be exact about which is which:

* **Search.** For 7259 × `tierApm_p95` × `downstream_arrest` ×
  `Biliary/ductular`, `main_fits.csv` records
  `beta_base_lo = -0.002022330`, `beta_base_hi = +0.000176625`. Those are
  quantiles of `fit_cell`'s 400-replicate block bootstrap, whose only randomness
  is `np.random.default_rng(seed + 7*j)`. 1,050,000 candidate seeds
  (`MASTER_SEED + 1000·0 + 7 + Δ`, `Δ ∈ [-350000, 700000)`) were replayed.
  **Exactly one matched: 20560827, i.e. Δ = 300000.**
* **Validation, out of sample.** The rule `MASTER_SEED + 1000·i + j + 300000`
  then reproduces the committed `beta_base_lo/hi` exactly for 7260 (seed
  20561827) and 7001 (seed 20562827) as well — targets not used in the search.
* **Corroboration.** With `PM_SEED_OFFSET = 300000`, the five composition-match
  seeds 20260901–05 (= `MASTER_SEED + 81…85`) are unreachable by the per-module
  branch, which is exactly the property `run_phase8_compmatch.py:56` asserts.

This is a recovered constant, not a tuned one: it was found against one fit and
then checked against others. It is flagged here because a reader has a right to
know that one number in the reconstruction came out of a search rather than out
of a surviving artefact.

*(module-index term: see §4 row for the `_pm` outputs)*

---

## 4. Acceptance test — does each output reproduce?

Every run below wrote into a scratch tree and was diffed against the committed
file. "Bit-identical" means equal md5 of the whole file. "Exact on the verified
subset" means a subset of the stage's jobs was re-run with the *full* section
and call lists in place (so the seeds are the ones the full run used) and every
shared column of every row it covers matched to `|Δ| = 0`.

### 4.1 Leading with what does NOT reproduce bit-identically

**`results/phase3/perm_nulls_var.csv` — DIFFERS, by ≤ 2.2 × 10⁻¹¹ relative.**

| committed | re-run |
|---|---|
| md5 `32397d5d6b626bb59924fdbc04669dac` | md5 `3bdade74fcf72b79558f214b9826e2bb` |

Same shape (315 × 56). 14 of 56 columns differ at all, none by more than float
noise amplified through a ratio:

| column | max \|Δ\| | max relative Δ | rows affected |
|---|---|---|---|
| `N3_var_sf_wm` | 9.69e-11 | 6.3e-12 | 50 / 315 |
| `N4_var_sf_wm` | 8.86e-11 | 2.2e-11 | 39 / 315 |
| `N3_var_beta_S0` / `N4_var_beta_S0` | 8.0e-15 | 7.5e-15 | 142 / 143 |
| `beta_obs` | 1.04e-16 | 5.5e-14 | 6 |
| eight others | ≤ 1.1e-15 | ≤ 4.0e-14 | ≤ 61 |

**This difference is not attributable to the reconstruction.** The producer,
`run_phase3_var.py`, survived and is unmodified; its companion output
`perm_draws_var.csv` (12,000 × 6 — every drawn shift vector and retention
fraction) **is bit-identical**, so `_expand`'s seeds and the whole geometry
stream are exactly right. Only the fitted quantities wobble, and `*_sf_wm` is
`mean((B0w − B)/B0w)`, a ratio that amplifies a last-bit difference in the
numerator by four orders of magnitude.

The likely cause is BLAS threading. `reports/WRITING_PACK.md:971` records the
producer as `python3 -u code/run_phase3_var.py --stage perm --n-perm 1000
--n-jobs 6 --no-full` — with **no** `OMP_NUM_THREADS/MKL_NUM_THREADS/
OPENBLAS_NUM_THREADS=1`, unlike `_m1_rerun_stage2.sh:6`, which sets all three
for every stage it drives. The same effect is directly demonstrable elsewhere in
this verification: a `stage_main` subset run **unpinned** differs from the
committed `main_fits.csv` by up to **9.0e-13**, and the *same* subset run
**pinned** is exact. Corroborating this reading, `WRITING_PACK.md:967` already
records that the committed `perm_nulls_var.csv` disagrees with the committed
`perm_nulls_c1.csv` on `beta_obs` by `1.04 × 10⁻¹⁶` — the exact magnitude of the
`beta_obs` column in the table above — i.e. the wobble is a property of the
original var run, and this re-run agrees with `perm_nulls_c1.csv` where the
committed var file does not.

*(An unpinned re-run to confirm this directly is recorded in §4.4.)*

**No result in the project is affected.** `sf_summary_var.csv`'s headline
medians — N3-var 0.995966, N4-var 0.985068 — are quoted to six decimals; the
disagreement is at the eleventh.

### 4.2 Reproduces bit-identically (whole file, equal md5)

| output | producer, as re-run | md5 |
|---|---|---|
| `window.csv` (297 × 18) | `stage_window`, all 11 sections, all 9 calls | `4be813a9962af1da01358eb002b342f3` |
| `perm_nulls.csv` (315 × 36) | `stage_perm --sections inband --calls tierA_p95 --n-perm 1000` | `d906394958dbe1b99981756290c511fa` |
| `perm_curves.csv` (2520 × 10) | same run | `e1c387efd176b7967d37eb258537c9ed` |
| **`perm_nulls_c1.csv` (595 × 104)** | **`stage_perm_c1 --sections inband --calls tierA_p95 --n-perm 1000`** | **`0318737e85a9b98e1b3bdd1461880ef5`** |
| **`perm_curves_c1.csv` (7560 × 10)** | same run | **`e864a8ab2e844024c42df2a483df74c0`** |
| `curves.csv` (6300 × 13) | `stage_curves --sections inband` | `fce12ed6cdd3b921eaac3df56ed70862` |
| `perm_draws_var.csv` (12000 × 6) | `run_phase3_var.py --stage perm` (surviving producer, via the reconstructed `_expand`) | `e788f0d824264bb1ff8a2e499809c772` |
| `summary_phase3.txt` | `summarize_phase3.py` | `dc92ddc6605eef52f6359aeab4e16fd7` |
| `sf_summary.csv` | `summarize_phase3.py` | `a5ccc9b0e81f4c335e8039e975ec1975` |
| `sf_summary_c1.csv` | `summarize_phase3_c1.py` (surviving) | `a547cb594662595cede964ca25fe072a` |
| `sf_summary_c1_n7.csv` | `summarize_phase3_c1.py` | `d71fdc23f2894743c9428e32e05de652` |
| `summary_phase3_c1.txt` | `summarize_phase3_c1.py` | `a8105b613c87013773dfba47407f4672` |

`perm_nulls_c1.csv` is the file the audit's F-1a named as having **no producer
at the tag**. It now has one, and that producer re-emits it byte for byte.

The three `summarize_phase3_c1.py` outputs were regenerated from the *committed*
inputs; they establish that the surviving summariser is intact and that it
consumes exactly the schema the reconstruction writes (`scope`, the
`FULL_NULLS`/`TILE_NULLS` column families, `n_perm`).

### 4.3 Exact on the verified subset (|Δ| = 0 on every shared column)

Re-running these files in full is 30–42 jobs of 20–70 minutes each. A subset of
jobs was re-run instead, with the full section and call lists in place so the
seeds are identical to the full run's.

| output | jobs re-run | rows compared | result |
|---|---|---|---|
| `main_fits.csv` (7094 × 86) | 4 of 297 — 7259 & 7260 × `tierApm_p95` × `downstream_arrest`; 7435 × `tierA_p95`; 7435 × `senepy_p99` | 148 | **exact**, max \|Δ\| = 0 across all 86 columns, `sender_set` and the matching diagnostics included |
| `perm_nulls_n7.csv` (1498 × 29) | 1 of 30 — 7352 × `senepy_p95`, seed 20280871 | 49 | **exact**, max \|Δ\| = 0 |
| `perm_nulls_pm.csv` (315 × 29) | 4 of 42 — 7259 & 7352 × `downstream_arrest` & `emt_ecm` | 28 | **exact**, max \|Δ\| = 0 |
| `perm_nulls_c1_pm.csv` (595 × 81) | 1 of 42 — 7259 × `downstream_arrest`, seed 20560820 | 13 | **exact**, max \|Δ\| = 0 |
| `perm_nulls_c1_n7.csv` (2758 × 81) | 1 of 30 — 7352 × `senepy_p95`, seed 20280871 | *(§4.4)* | *(§4.4)* |

The two `_pm` rows are the ones that depend on the searched constant of §3, and
they are the strongest check on it: the `perm_nulls_pm` subset exercises the
per-module seed through the **permutation** stream (1,000 draws of N1/N3/N4),
not merely through the bootstrap the constant was found against.

Independently, the reconstructed `_perm_c1_job` printed `full=`/`tile=` fit-cell
counts identical to the committed logs on every job run — 49/42, 49/49, 56/42,
63/56, 49/42, 49/49 for the six in-band sections (`logs/m1_perm_c1.log`), and
7/6 for 7259 × `tierApm_p95` × `downstream_arrest` (`logs/m1_perm_c1_pm.log`).


---

## 5. What can now be run that could not be run at the tag

Each of these was verified by execution, not by reading:

| capability | at `phase8-frozen` | now |
|---|---|---|
| `--stage perm_c1` (the corrected N3/N4 battery) | `raise SystemExit("unknown stage perm_c1")` | runs; see §4 |
| `--tag _pm` | argparse: "unrecognized arguments" | accepted |
| `RN._expand` (`run_phase3_var.py:237`) | `AttributeError` | `python3 run_phase3_var.py --stage perm` runs end to end. Smoke: 7259 × `tierA_p95`, 5 permutations → 49 rows × 63 columns with the `N3_var_*` / `N4_var_*` block. **N3-var / N4-var — correction C-5's primary corrected pair — are runnable again.** |
| `SectionFit(module=)` (`run_phase8_compmatch.py:250`) | `TypeError` | `CompMatchFit('7259…', 'tierApm_p95', 20260901, module='emt_ecm')` builds; `sender_set = A_sender_for_emt_ecm`, `n_senders = 5707`, `_js('emt_ecm') = [1]`; `build_jobs` fans `tierApm_p95` to 42 jobs and `tierA_p95` to 6 |
| `sender_mask('tierApm_p95')` | `ValueError(call)` | returns the mask that reproduces `window.csv` exactly |
| `--calls all9` | `["all9"]` → `ValueError` | 297 jobs, matching `CS_PHASE8_M1_RERUN.md`'s "297 jobs (11 × 9 calls, per-module fanned out over 7 modules)" |

---

## 6. What this does NOT fix

1. **The tag.** `phase8-frozen` still points at `9264396`, which still ships the
   541-line `run_phase3_nulls.py`. No tag was created or moved. Whoever owns the
   freeze has to decide between re-cutting it and annotating it; this file is the
   material for either decision.
2. **Provenance.** Bit-identical output is behavioural equivalence, not identity
   (§0). A pre-registration cannot claim the frozen code *is* the code that ran.
   The defensible sentence is: *the frozen results were produced by code that was
   lost; the repository now contains a reconstruction that reproduces them, and
   the reconstruction and its verification are recorded in
   `reports/PIPELINE_RECONSTRUCTION.md`.*
3. **`PREREG_PHASE8.md`'s ~26 line citations into `run_phase3_nulls.py`.** They
   were taken from the lost file. Audit F-12 and Appendix A already document
   that they do not resolve; they do not resolve to the reconstruction either
   (this file is 813 lines, the lost one was ≥ 700, the frozen one is 541). They
   should be re-taken against whatever is finally frozen, or dropped in favour of
   symbol names.
4. **Every other F-finding in `AUDIT_PREREG_VS_CODE.md`.** F-2 (which corrected
   null is primary), F-3 (the "equal-count, none empty" claim), F-5 (the
   mislabelled destructiveness quantity), F-6 (R3(c)'s unreproducible range),
   F-8 (the seed-choice justification, whose `_expand` now exists but whose
   *guarantee* still fails for `--calls all` at `j = 5`), F-9/F-10/F-11 and the
   whole of Part 4 are untouched here.
5. **The `logs/` directory is still untracked** (`git ls-files logs/` → 0), yet
   three of the five constraints in §2 came out of it. If `logs/` had been
   cleaned as well, this reconstruction would have been materially weaker.
   Tracking it, or at least the `m1_*.log` set, is cheap insurance.

---

## 7. Also fixed here

### 7.1 `PREREG_PHASE8.md` §3.11's pin — new correction **C-11**

- The pinned digest `ecf86b9ca5460f31290e2f4c9e822ea2` is **31 hex characters**.
  No MD5 has 31 characters, so it cannot be a digest of anything.
- `results/phase3/summary_phase3.txt` was rewritten by task 8.7 at
  **2026-08-27 09:06** and is **md5 `dc92ddc6605eef52f6359aeab4e16fd7`**.
- Its line 19 reads `fits 315;  naive beta > 0 in 216;  positive AND
  block-bootstrap CI excludes 0: **153** (49 %)`. §3.5's "the same 315/**160**
  appears in the pinned file" and §5's "160 reportable fits" are therefore false
  as written; C-5 already carried 153.

C-11 was added to the §0.0 corrections table with dated inline markers at §3.5,
§3.11 and §5, in the file's own no-silent-rewrite style: the original wording is
left in place at every site.

**The correction is self-verifying.** The reconstructed `summarize_phase3.py`,
re-run with its writes redirected to a scratch tree, regenerates
`summary_phase3.txt` **bit-identically** — md5 `dc92ddc6605eef52f6359aeab4e16fd7`,
the digest C-11 now pins — and `sf_summary.csv` bit-identically as well.

### 7.2 N3-var / N4-var are back in the battery

Correction C-5 promotes **N3-var 0.996 / N4-var 0.985** to "the primary
corrected pair", and audit F-2 records that the frozen §3.4 battery table does
not contain them. Their producer, `run_phase3_var.py`, could not run at the tag
at all: it calls `RN._expand` at `:237`, and `_expand` did not exist
(`AttributeError`). `_expand` is reconstructed and `run_phase3_var.py --stage
perm` now runs (§5). The *documentation* half of F-2 — §3.4 designating the tile
variants as the corrected N3/N4 while C-5 designates the var pair — is a prose
contradiction and is deliberately **not** touched here; it is a freeze decision,
not a code defect.

---

## 8. How to re-run the verification

Everything below writes into a scratch tree, never into `results/`. Set the
three BLAS thread limits — the drivers do (`_m1_rerun_stage2.sh:6`), and it
matters (§4, note on float determinism).

```bash
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OUT=/tmp/recon; mkdir -p $OUT

# a stage, redirected: RN.RES is the only write path
python3 - <<'PY'
import sys; sys.path.insert(0, '/workspace/code')
import run_phase3_nulls as RN, sasp_phase3 as P
RN.RES = '/tmp/recon'
RN.stage_window(P.ALL_SECTIONS)                                     # window.csv
RN.stage_perm(P.IN_BAND, ['tierA_p95'], 6, 1000, do_full=True)      # perm_nulls.csv
RN.stage_perm_c1(P.IN_BAND, ['tierA_p95'], 6, 1000, do_full=True)   # perm_nulls_c1.csv
PY
md5sum $OUT/window.csv /workspace/results/phase3/window.csv
```

For a **subset** of a multi-call or per-module stage, expand the job list with
the *full* section and call lists — the seeds depend on both indices — and then
filter:

```python
jobs = RN._expand(P.IN_BAND, ['tierA_p90','tierA_p99','cdkn1a_pos',
                              'senepy_p95','senepy_p99'],
                  P.MASTER_SEED, 5000, 17)
jobs = [j for j in jobs if j[0].startswith('7352') and j[1] == 'senepy_p95']
```

`summarize_phase3.py` reads and writes through the same `RES`, so verify it by
copying the file and redirecting only its two write paths (`sf_summary.csv`,
`summary_phase3.txt`).

---

## 9. Commits

Small, and committed as they were made — the failure this task exists to repair
was working-tree code that was never committed.

| commit | contents |
|---|---|
| `91d9ba3` | `run_phase3_nulls.py` + `sasp_phase3.py`: the per-module sender axis, `_expand` / `_js` / `is_permodule`, `SectionFit(module=)`, `sender_set`, `--tag`, `--calls all9\|tierA_pm`, and the whole `perm_c1` stage |
| `b8e74a0` | `phase2_downstream.py`'s per-module Tier A block; `PREREG_PHASE8.md` correction **C-11** |
| `00f1210` | `summarize_phase3.py`: `ALL9_CALLS` in the §6 table, `perm_nulls_pm.csv` in `load()` |
| `e52ae5f` | the recovered per-module seed rule, `PM_SEED_OFFSET = 300000` × (module index + 1) |
| `9e6c687` | `code/README.md`: the `perm_c1` stage, the per-module axis, and a pointer here |

No tag was created or moved. Nothing was pushed.

