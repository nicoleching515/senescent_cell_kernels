# Reproducibility repair, 2026-08-27/28

What an independent party can and cannot do after cloning this repository and trying to
regenerate the headline result. Written against `reports/AUDIT_REPRODUCIBILITY.md`, which is
the specification for this work; its finding labels (B1-B10, D1-D8) are used throughout.

Every claim below was exercised. Where a producer is said to reproduce a committed file, the
comparison was run in this session and the command is given. Where something was only written
and not verified, it says so.

---

## 1. What an independent party STILL cannot regenerate

Ranked by how much it costs them.

### 1.1 Nothing in the pipeline has been re-run end to end, so the headline numbers are still only *re-derivable in principle*

This repair made the producers exist and proved that the cheap ones reproduce their committed
bytes. It did **not** re-run `--stage main`, `--stage perm` or `--stage perm_c1` — that is
hours of compute on a shared 57.7 GB cgroup with another agent live in it, and the task
forbade it. So the chain from `data/processed/cache3/*.npz` to `main_fits.csv`,
`perm_nulls.csv` and `perm_nulls_c1.csv` remains **unverified by execution**. Everything
downstream of those three files is now verified; those three are not.

`code/_m1_stage1.sh --run` is the command that would settle it. Its dry run passes preflight
today. Only a full re-run can close this.

### 1.2 `results/zonation_gene_correlations_7250_sham.csv` — a gene-set-defining input with no producer, and one input that is not written down

`code/build_genesets.py:166` reads it as `ZONATION_CSV` and turns it into
`D_zonation_pericentral` and `D_zonation_periportal`. The method is documented at
`build_genesets.py:157-165` and pins the seed axis, the 60 % cut (0.6 x 236,905 = 142,143
rows — confirmed arithmetic), the 5,106-gene panel and both thresholds. It does **not** pin
which hepatocyte-marker list scored the cells for that cut, nor whether it was
`sc.tl.score_genes` and with what `ctrl_size`, nor whether `detection_rate` is over all cells
or over the 142,143. `TIER_E['E2_hepatocyte_identity']` is the plausible candidate and is
nowhere stated to be it.

A rewrite would land close but not exact, and a near miss silently changes gene-set
membership. **No producer was written.** The gap is now recorded in a comment at the read
site, and the file is treated as a pinned input like the MSigDB JSON. If it is ever rebuilt,
diff the resulting gene lists, not the r values.

### 1.3 `results/phase4/cellchat_summary_statistic.csv` — producible, but only behind a pipeline stage

Every primitive exists (`phase4_methods._group_summary` already implements trimean/truncmean/
mean; `phase4_tiles.all_tiles()` yields the 18 tiles; the committed 144 rows = 18 tiles x 4
LR pairs x 2 roles is consistent). But `/tmp/p4/cache4` no longer exists, so a verification
run needs `phase4_data.build()` for six sections first — 10-30 minutes of pipeline. Out of
scope by instruction. **No producer was written**, because writing one that cannot be
verified against the committed file is exactly the failure mode this repair exists to end.

### 1.4 `results/section_qc_sender_summary.csv` reproduces 120 of 132 cells, and one of its columns rests on an invented threshold

It now has a producer (`code/section_qc_sender_summary.py`, new) and a `--check` mode that
prints a per-column, per-section grid. Seven of the eleven rows come back byte-identical.
The twelve cells that differ are `n_analysable`, `n_cdkn1a_pos` and `cdkn1a_pos_pct_all` on
7239, 7435, 7248 and 7001 — sections re-annotated **after** the committed CSV was written
(7001: 139,859 against 159,279 analysable cells today). That is drift in the inputs, not a
defect in the recipe, and `cdkn1a_pos_pct_hepatocytes` still matches exactly on all four,
which says the re-annotation moved cells between `Low_quality`/`Unknown` and non-hepatocyte
types. No attempt was made to force those four to match.

`portal_triad_valid`'s rule is recorded nowhere. The committed table pins the threshold only
to (0.121, 0.198]; the script uses 0.15 and states, on every run, that this is a midpoint and
not a recovered value. Any value in that interval reproduces all eleven calls, so that row
matching is **not** evidence the threshold is right — a twelfth section landing in the
interval would expose it.

### 1.5 `results/phase8_d2/dca_worker_meta_*.json` — unreconstructable as committed

`dca_denoise_worker.py` writes its meta into a tempdir that `_shims_dca_bridge/dca/api.py`
then discards; the two committed files were copied out by hand. Even with a copy step, the
content cannot come back: `minutes` is wall clock, and `out_min`/`out_max`/`out_mean` come
from a 300-epoch TensorFlow 2.4 autoencoder. One of the two also carries `n_genes = 4601`
against 4845 for the full run, i.e. a smoke-run gene filter whose driver is not in the repo.
Treat as historical provenance.

### 1.6 `results/phase3/m1_final_audit.txt` section 5, and the 14 `.log` files under `results/`

Section 5 ("FIGURE STATE") embeds the live figure-guard verdict and the md5 + mtime of every
figure. The committed copy was captured mid-8.7 **while the guard was failing**, so it can
never come back; today the same producer prints `OK: all 52 committed figures match`.
Sections 1-4, which carry every number the reports cite from that file, do reproduce
byte-for-byte (§2.1). The 14 `.log` files are stdout captures with no numbers in them; no
producer was written for them.

### 1.7 The environment is documented but was not rebuilt from the manifest

`requirements.txt` now names the interpreter and every driver resolves it. Nobody built a
fresh venv from `requirements.txt` in this session and re-ran anything under it — the existing
`/workspace/envs/sasp311` was used throughout. Four imported packages remain deliberately
unpinned (`torch`, `esda`, `libpysal`, `setuptools`), and `torch` in particular installs
`2.4.1+cu124` here, a build served only from `download.pytorch.org`, so `pip install -r
requirements.txt` resolves a different one. DeepScence scores are a function of the torch
version. **Unchanged by this repair, and still a live divergence risk (D3b).**

### 1.8 Seeds still depend on `--calls` argument order (B5)

`run_phase3_nulls.py` seeds each job `base + step_i*i + step_j*j` over the *positions* of the
section and the call. Running the same calls in a different order is a different permutation
null, with no warning. `run_phase3_nulls.py` and `sasp_phase3.py` were off limits in this
session (Phase 10 imports them live), so this is **documented, not fixed**: the keyword forms
`--calls all9 | all | tierA_pm` expand from constants in the script and are order-stable, and
`code/_m1_stage1.sh` uses them and says so. A free-form comma list is still a hazard.

### 1.9 ~20 pre-existing section-10 violations in the committed corpus

The new checker finds them; the pre-commit ratchet stops the count growing but does not fix
what is already there. They are prose fixes for the authors, not code fixes, and rewriting
scientific prose was not in scope. Full list in §4.3.

---

## 2. What now regenerates, and how that was checked

### 2.1 The B7 producers — 15 artefacts, one command

```
bash code/_repro_artefacts.sh --check     # regenerates into a temp dir, diffs, writes nothing
bash code/_repro_artefacts.sh --write     # regenerates in place
```

Result today: **match=15 differ=0 skipped=0**, in about a minute.

| Artefact | Producer | Status |
|---|---|---|
| `results/phase3/caller_*_2sec_c6.csv` (7 files) | `caller_disagree_all.py --set 2sec_c6` | **byte-identical**, all seven |
| `results/phase3/m1_final_audit.txt` | `m1_final_audit.py` | **sections 1-4 byte-identical** (incl. `power80_bound` 0.1833 and `ctrl_amp_med` 0.0288); section 5 is live state, see §1.5 |
| `results/phase3/m1_prepost_main_fits.txt` | `m1_compare_modules.py results/phase3_pre_c6 results/phase3` | **byte-identical** |
| `results/phase3/m1_n7_prepost.txt` | **new** `code/m1_n7_prepost.py` | **byte-identical** |
| `results/phase3/sf_summary_c1_swap_vs_n1.csv` | **new** `code/sf_swap_vs_n1.py results/phase3_pre_c6` | **byte-identical** |
| `results/phase8_d2/d2_tables.md` | `report_d2_tables.py --out ...` | **byte-identical** to the regenerated committed copy |
| `results/phase8_d2/committed_deepscence_sha256.txt` | `sha256sum ...` (in the driver) | **byte-identical** |
| `results/phase8_d2/dca_venv_pip_freeze.txt` | DCA venv `pip freeze` | **byte-identical** |
| `results/phase8_d2/dca_venv_python.txt` | DCA venv `python -VV` | **byte-identical** |
| `results/a3_fallback/gpl33762_count.xml` | `curl` eutils (in the driver, `--with-network`) | command committed; **byte-identity checked today by a separate agent against the live NCBI response, not re-checked by me** |

One more producer was written but is **partial by design**:
`code/section_qc_sender_summary.py` (see §1.4). `_repro_artefacts.sh --check` runs its
`--check` and prints the result as information, without counting it in the totals.

Three more B7 rows were already closed by earlier commits and are **verified here for the
first time**, by re-running the producers with `to_csv`/`savefig` redirected to a scratch
directory so nothing in the tree moved: `results/phase3/figure2b_data.csv`,
`figure2c_data.csv`, `figure2d_data.csv` (from `make_figure2bc.py`) and `figure2e_data.csv`
(from `make_figure2e.py`) all come back **byte-identical**.

### 2.2 The caller-coverage headline now has a runnable producer (B3, B7)

`caller_disagree_all.py:62` called `CD.DS_ALIAS`, which did not exist — the one-line edit its
own docstring says was made was never committed, so **every** invocation raised
`AttributeError`. `DS_ALIAS` is now in `caller_disagree.py` as `{v: k for k, v in SAMP.items()}`,
and `caller_disagree_all.py` has a real CLI (`--set {2sec_c6,2sec,11sections}`, `--suffix`,
`--out-dir`). Running `--set 2sec_c6` into a scratch directory reproduces all seven committed
`caller_*_2sec_c6.csv` byte-for-byte. Those seven are the frozen post-C6 base of
`summarize_caller_coverage.py:67`, hence of `caller_coverage_gate_headline.csv` and of
`README.md:299`.

The same script used to drop DeepScence silently when `data/processed/deepscence_*.csv` was
absent (gitignored, so absent on a fresh clone) and exit 0 with a different caller set. It now
refuses, and names the file it wanted (D6).

### 2.3 The head of the M1 pipeline exists (B6)

`code/_m1_stage1.sh`. Every other M1 driver *waits* for stage 1 — on `pgrep`, on line counts
in `logs/m1_{main,perm,perm_c1}.log`, or on a PID — and on a clean machine every one of those
waits falls through immediately. This script runs the three primary stages, writes exactly the
three log files the chain waits on, and hands over to `_m1_chain.sh`.

The invocations were recovered from two records that agree: the `set -x` trace in
`logs/m1_chain.log` plus the section/call/`n_jobs`/`n_perm` headers of `logs/m1_main.log`,
`logs/m1_perm.log`, `logs/m1_perm_c1.log`; and `WRITING_PACK.md:378` /
`CS_PHASE8_M1_RERUN.md:478`.

It is a **dry run by default** and preflights three things:

- the project interpreter, via `code/_env.sh`;
- a `cache3` npz for **all eleven** sections. `run_phase3_nulls.py`'s own guard (added
  earlier today) fires only when the filter empties the list, and `--stage main` writes
  `main_fits.csv` with no tag and no merge — so a *partial* cache still silently replaces the
  nine-call table with a smaller one;
- that `run_phase3_nulls.py`/`sasp_phase3.py` still carry `perm_c1`, `--tag`, `all9`,
  `is_permodule`, `_expand` and the `tierApm_p*` sender mask.

That last check was verified in both directions: it passes on the current tree, and against
the pre-repair tag commit `9264396` it reports **all five missing**.

**Note on B1/B2/B3.** B1 and B2 were repaired before the tag was re-cut: `phase8-frozen` now
points at `d04691e`, whose `run_phase3_nulls.py` has `perm_c1`, `--tag`, `all9`,
`is_permodule` and `_expand`, and whose `sasp_phase3.py` resolves `tierApm_p*`. The audit's
B1/B2 describe the **original** `9264396` tag commit. B3 was still open this morning and is
closed by this work.

---

## 3. Environment (B9, B10, D2, D3)

**One place names the interpreter.** `code/_env.sh` sets `SASP_PYTHON`
(default `/workspace/envs/sasp311/bin/python`), prepends its `bin/` to `PATH` so a bare
`python3` inside a driver also lands there, sets the single-thread BLAS variables, and
**exits 1 with rebuild instructions** if the interpreter is missing — it never falls back to
the system python. 24 tracked drivers source it (`code/h1_*` deliberately untouched: another
agent is live in them). Verified: sourcing puts `sys.executable` in the venv;
`SASP_PYTHON=/nonexistent` exits 1 with the message.

**The commot shim is now on all three call sites (B9).** `phase4_commot_mechanism.py` and
`phase4_positive_control.py` imported `commot` with no shim on any path, so neither could run
under the pinned `numpy==2.4.6` — and the former is the sole producer of
`results/phase4/commot_mechanism.csv`, the Figure 4 mechanism result. Verified by execution:
`import phase4_commot_mechanism` now succeeds, the shimmed import inside `_commot()` succeeds,
and a bare `import commot` still raises `AttributeError: np.Inf`.

**Package data no longer comes from the container overlay (D2).** Eight scripts hardcoded
`/usr/local/lib/python3.11/dist-packages/{DeepScence,senepy}/data/...`. `code/_pkgdata.py`
resolves it in order: `$SASP_DEEPSCENCE_DATA`/`$SASP_SENEPY_DATA`, then `importlib` in the
running interpreter, then the project venv, then the overlay last — and raises with all four
candidates listed rather than guessing. The two copies of `coreGS_v2.csv` are byte-identical
(`md5 b981b9e9e730217d339306a709ada201`), so no number moves; the only visible change is
`gate_result_human.json`'s `corescence.source`, which now records the venv path. Verified:
`gate_genesets_guard.py` still exits 0, mouse and human, under **both** interpreters. Note
the consequence: that one field now records whichever interpreter ran the gate, so running it
with `/usr/bin/python3` rewrites the line back to the overlay path. The committed value is the
venv one; every count in the file is interpreter-independent.

**The DCA runners honour `DCA_ENV_ROOT` (D3).** Four of them hardcoded a per-session `/tmp`
path that dies with the pod. The old path survives as the last fallback and a missing
interpreter is now fatal with build instructions.

**Documentation.** `README.md` describes both environments, the overlay hazard, how to run a
script directly, and the four unpinned imports; `requirements.txt` carries the same statement
next to the pins.

---

## 4. The section-10 prohibitions are now enforced (PART 4 of `AUDIT_PREREG_VS_CODE.md`)

### 4.1 The checker

`code/check_prohibitions.py` mechanises nine of the twelve — 10.1, 10.2, 10.6, 10.7, 10.8,
10.9, 10.10, 10.11, 10.12 — over `reports/`, `results/` and `README.md`, at **paragraph**
scope, so the co-occurrence rules mean what the prereg means: the companion number must appear
beside the forbidden one. 10.3, 10.4 and 10.5 are semantic; `--list` prints them as
REVIEW-ONLY rather than leaving them in a list that implies a coverage the code does not have.

Two waivers, both narrow: a paragraph that is *about* the prohibition (the rule itself, a
withdrawal, a superseded passage kept as a record) is meta-discussion, and `PREREG_PHASE8.md`
itself, where the forbidden literals live by definition. A per-paragraph escape hatch,
`prohibition-waiver: 10.x <reason>`, exists for anything else.

`--self-test` proves every rule fires on a synthetic violation **and** that the compliant form
of the same number does not, and that the meta waiver works. It passes.

### 4.2 Wiring, and proof that it fires

- `.githooks/pre-commit` runs `--staged`. It is a **ratchet**: a violation already in HEAD's
  copy of a file does not block, one the commit *adds* does.
- `code/hook_geneset_gate.sh` — the only always-on hook — used to `exit 0` on `*/results/*`
  and `*/reports/*`, i.e. every path where a violation would be written. Those paths now run
  the checker (exit 2 = blocking feedback); gene sets still run the section-11 gate.
- `code/recut_phase8_gate.sh` runs the self-test, reports the backlog, and requires
  `_repro_artefacts.sh --check` to pass before the tag may be re-cut.

Exercised, not assumed:

1. A deliberate two-violation file made the PostToolUse hook **exit 2** and the pre-commit
   hook **block the commit** (HEAD unmoved). Test file removed.
2. Appending to a file that already carries backlog entries does **not** block (3 carried,
   0 new).
3. The ratchet blocked one of my own commits — the README's "the D2 `denoise=True` arm
   (TensorFlow < 2.5 ...)" tripped 10.10 on two version numbers. The rule was narrowed to
   require a value in [0,1) or a percentage; the self-test still passes.
4. It caught a violation the **pipeline itself was generating**: see §4.4.

### 4.3 The backlog: 22 pre-existing violations at the time of writing, unfixed

Not blocking, listed so they are not lost. `python3 code/check_prohibitions.py --backlog`
reproduces the current list, which moves as authors fix them — it stood at 22 when this table
was built and at 21 an hour later, after the Phase 10 author put the seed-stability companion
back beside three `denoise=True` numbers in `CS_PHASE10_TWO_ARM.md` (commit `2e46e8c`) in
response to the new pre-commit ratchet. That is the mechanism working on live text.

| rule | count | files |
|---|---|---|
| 10.8 (the matched number quoted without its covariate-adjusted companion) | 8 | `CS_PHASE8_COMPMATCH.md` (3), `CORRECTIONS.md`, `COMPLETED_TASKS.md`, `WRITING_PACK.md`, `PLAN_UPDATE_D12_D13.md`, `AUDIT_PREREG_VS_CODE.md`, and `results/phase3_pre_c6/summary_phase3.txt` (a frozen snapshot — cannot be edited) |
| 10.9 (`mor` without the `lib` result beside it) | 4 | `CS_PHASE8_D2_DENOISE.md` |
| 10.10 (a single-seed denoised number with no stability companion) | 3 | `CORRECTIONS.md` (2), `COMPLETED_TASKS.md` |
| 10.2 (an affirmative caller-independence claim) | 3 | `CS_PHASE5.md`, `WRITING_PACK.md`, `SUBMISSION_PATCH_2026-08-29.md` |
| 10.11 (the withheld D2 `raw` depth correlation) | 1 | `CS_PHASE8_D2_DENOISE.md` — reads as a caution about the value, but quotes it |
| 10.7 (the withdrawn circularity range) | 1 | `BIO_PHASE3.md` |

Some of these are near-misses of the meta waiver rather than real breaches; each needs a human
read. The rule is now mechanical, so the argument is at least about specific paragraphs.

### 4.4 A committed script was handing the author a forbidden number

`report_d2_tables.py` printed `rho_signed_dz_vs_depth` for the D2 `raw` rows — the exact pair
section 10.11 forbids — pre-formatted, in a report table. Those cells are now `WITHHELD` with
the reason printed beside them, `--emit-prohibited` exists for internal diagnostics and
refuses to write a file, and `--out` gives `results/phase8_d2/d2_tables.md` a producer it never
had. The regenerated file passes the checker. Its section H now also states how many
`runmeta_*.json` it globbed: the row count is a function of disk state, and a rebuilt tree
renders one more row (`mor`/7352) than the committed copy carried.

---

## 5. Silent degradation made loud (D1, D4, D6)

| Was | Now |
|---|---|
| `figure2a` regenerated byte-identically from its own **tracked, guard-watched** cache, so `check_figures_guard.py`'s pass on it was vacuous (D1) | `make_phase5_figs.py` **refuses** (SystemExit, naming the offending input and both mtimes) when any `cache3` npz, gene set or `senders_*.csv` is newer than the cache. `--rebuild` recomputes, `--accept-stale` overrides. Verified: a synthetic old cache produces the refusal (55 newer inputs); the real cache is fresh (0), so the committed tree still runs |
| `corescence_circularity.load()` silently returned the tracked JSON, so the human gene-set gate compared against a committed cache, never a fresh derivation (D1) | it says so on stderr, and honours `fresh=True` / `SASP_CORESCENCE_FRESH=1` |
| `caller_disagree_all.py` dropped DeepScence silently when its gitignored input was absent (D6) | refuses, naming the missing file and the rebuild script; `--allow-missing-deepscence` is explicit |
| `figures/.committed_manifest.json` was gitignored, so on a clone the figure guard exits "no manifest; run `--snapshot` first" and `--snapshot` re-baselines whatever is on disk. The tag's own "52/52 match" claim was unverifiable from the tag (D4) | tracked. Verified: the guard passes at 52/52 against the committed manifest |
| `report_d2_tables.py` section H silently changed row count with disk state | prints the glob count |
| `run_figure2a.py` skipped any (section, module) whose `fit_*.csv` checkpoint exists — and 35 of those checkpoints are **tracked**, so a fresh clone recomputes nothing | prints how many committed checkpoints it is reusing, and `--force` recomputes. Exercised: the notice fires |

Left alone deliberately: `sasp_sweep.py`'s checkpoint skip (has `--force`, documented in
`code/README.md`), and the `prepare_samples.py` / `phase4_data.py` / `run_deepscence_*.py`
skip-if-exists paths (their outputs are gitignored, so a fresh clone is safe). `sasp_phase3.py:99`
and `run_phase3_nulls.py` were off limits this session.

---

## 6. What the tag would capture if re-cut today

Added since `d04691e`: the `DS_ALIAS` edit; a real CLI and a DeepScence guard on
`caller_disagree_all.py`; `code/m1_n7_prepost.py`; `code/sf_swap_vs_n1.py`;
`code/_repro_artefacts.sh`; `code/_env.sh`; `code/_pkgdata.py`; `code/_m1_stage1.sh`;
`code/check_prohibitions.py`; `code/section_qc_sender_summary.py`; the shim on two `commot`
call sites; the DCA runners' honouring
of `DCA_ENV_ROOT`; `figures/.committed_manifest.json`; the section-10 wiring in
`.githooks/pre-commit` and `hook_geneset_gate.sh`; the extended re-cut gate; and the
environment documentation.

Still **not** in it: `logs/` (227 files, gitignored — and the only record from which the
stage-1 invocations were recovered; `code/_m1_stage1.sh` now carries them in tracked form, but
the raw logs remain the primary evidence), `data/processed/cache3/`, the venvs, and a
zonation-table producer.

`git push` was not run — there are no credentials in this environment, and no tag was created
or moved. `code/recut_phase8_gate.sh` currently reports every check in its new section 5 as
`ok`.

---

## 7. Method, and what was not done

Static reading plus execution of the cheap producers. No pipeline stage was run. Nothing under
`results/` was written except `results/phase8_d2/d2_tables.md` (regenerated deliberately, one
withheld cell and one added provenance row) and `results/phase7_jobA/gate_result_human.json`
(one path line, from re-running the gate). Verification runs wrote into a scratch directory
via `--out-dir` or a redirected `to_csv`; `git status` was clean between commits.

Two facts worth recording for whoever audits this next:

- `results/phase3/sf_summary_c1_swap_vs_n1.csv` is a **pre-C6 artefact**. The committed
  post-C6 copy is byte-identical to the `phase3_pre_c6` one. `sf_swap_vs_n1.py` therefore
  regenerates it from `results/phase3_pre_c6`; run with no argument it produces the post-C6
  numbers, which differ. `CS_PHASE8_C1_CLOSEOUT.md:109` quotes the pre-C6 values and
  `CORRECTIONS.md:756` the post-C6 ones.
- `results/real_fits.csv` and `results/real_curves.csv` reproduce in **content but not in
  byte order**: `run_figure2a.py` concatenates in `os.listdir` order, so the row order is
  filesystem order. Verified: 280 and 1,050 rows, equal after sorting, and
  `sorted(os.listdir(...))` does not reproduce the committed order either. Recorded in a
  comment at the concat; anyone diffing those two files must sort first.
