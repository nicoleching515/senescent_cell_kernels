# Reproducibility audit of `phase8-frozen`

Static audit plus cheap import/path/CLI checks. **Nothing was re-run**; where only a
re-run could settle a question it is said so. Read-only except this file.

- Tag: `phase8-frozen` → commit `926439629a07269a32c93f998da0f6e1cd20933c` (annotated tag object `78e42de`).
- Working tree at audit time: `HEAD` = `8f69ca6`, 12 commits ahead of the tag, clean apart
  from the human arm (`code/h1_*`, `data/processed_h1/`, `results/phase9_h1/`), which is out of scope.
- The three checksums in the tag annotation still match on disk
  (`perm_nulls.csv d906394…`, `sf_summary.csv a5ccc9b…`, `summary_phase3.txt dc92ddc…`).
  **The frozen artefacts are intact. The problem is regenerating them.**

---

## Verdict

An independent party who clones this repository at `phase8-frozen` **receives every
cited result file** (555 tracked files under `results/`; every file cited by a report
exists on disk and is tracked). They **cannot regenerate a substantial part of it**, and
in several cases their first attempt will **destroy the frozen artefacts and exit 0**.

The root cause of the largest finding is a single one: **`code/run_phase3_nulls.py` in the
repository is an older revision than the one that produced the Phase 8 results.** It was
committed exactly once (`a6aac3a`) and never updated, while every script that depends on it
was. The working copy that grew the missing features was lost when the container was wiped.

---

## BREAKS reproducibility

### B1. The committed `run_phase3_nulls.py` cannot produce the results attributed to it

`code/run_phase3_nulls.py` is byte-identical at `phase8-frozen` and `HEAD`, has one commit
in its entire history, and that version lacks four things the frozen tree depends on.
Verified by execution against `/workspace/envs/sasp311/bin/python`:

```
$ python run_phase3_nulls.py --stage perm_c1 --sections … --calls tierA_p95
unknown stage perm_c1

$ python run_phase3_nulls.py --stage perm … --calls tierApm_p95 --tag _pm
error: unrecognized arguments: --tag _pm
```

| Missing | Referenced by (tracked) | Frozen artefacts with no working producer |
|---|---|---|
| `--stage perm_c1` | `_m1_rerun_stage2.sh`, `WRITING_PACK.md:380,1452` | `perm_nulls_c1.csv`, `perm_curves_c1.csv`, `perm_nulls_c1_n7.csv`, `perm_curves_c1_n7.csv`, `perm_nulls_c1_pm.csv`, `perm_curves_c1_pm.csv` |
| `--tag` | `_m1_rerun_stage2.sh:24,28` | `perm_nulls_pm.csv`, `perm_curves_pm.csv` |
| `is_permodule`, `_js` | `run_phase8_compmatch.py:446,483,540` | `compmatch_fits_*.csv` |
| `_expand` | `run_phase3_var.py:237` (its docstring line 19 states `_expand` "comes from" `run_phase3_nulls`) | `perm_nulls_var.csv`, `perm_draws_var.csv`, `perm_nulls_var_full200.csv`, `perm_draws_var_full200.csv` |

A static cross-module attribute check over all 158 tracked scripts found exactly five
unresolved references; four are the above, the fifth is B2.

The apparent coverage is a trap: `run_phase3_nulls.py:468` writes
`f"{RES}/perm_nulls{tag}.csv"`, which *looks* like it produces the whole `perm_nulls*`
family. It does not — the only assignment to `tag` is line 467,
`tag = "" if calls == [PRIMARY_CALL] else "_n7"`. The committed script can emit exactly
`perm_nulls{,_n7}.csv`, `perm_curves{,_n7}.csv`, `main_fits.csv`, `window.csv`, `curves.csv`.
Nothing else.

**C1 (the corrected in-tissue N3/N4 nulls) and the second pre-registered Tier A variant
are two of the four changes the tag annotation names as the point of the freeze. Neither
has a runnable producer in the tag.** `perm_nulls_c1.csv` is cited by twelve reports
including `PREREG_PHASE8.md`, and is a hard input of `summarize_phase3_c1.py:133`.

**Recovery search — negative.** The lost revisions are not in this repository by any route:
`find /` returns exactly one `run_phase3_nulls.py`; `git stash list` is empty; the reflog holds
only the 12 post-tag commits; the one dangling commit (`5dd8309`, an abandoned commit from
today's reconciliation) carries a `run_phase3_nulls.py` byte-identical to `HEAD`'s; and a grep
of all 30 dangling objects for `perm_c1`, `is_permodule`, `DS_ALIAS`, `def _expand` and
`tierApm_p95` finds nothing. Only a copy held outside this workspace, or a re-derivation from
the specification in the reports, can restore them. **This is the single highest-value action
available.** Until then, only a re-run could confirm the frozen numbers reproduce.

### B2. `sasp_phase3.py` cannot evaluate three of the nine sender calls in `main_fits.csv`

`code/sasp_phase3.py:217` `Sec.sender_mask()` handles `tierA_p*`, `cdkn1a_pos`, `senepy_p*`,
and `raise ValueError(call)` otherwise. `tierApm_p95` does not match `tierA_p` and falls to
the raise. But the frozen `results/phase3/main_fits.csv` — the file carrying the headline
controlled amplitude 0.0288 and naive 0.3288 — contains nine calls:

```
tierA_p90 812  cdkn1a_pos 812  tierApm_p90 812  tierApm_p95 805  tierA_p95 805
senepy_p95 798  tierApm_p99 773  tierA_p99 770   senepy_p99 707
```

Three (`tierApm_p90/p95/p99`) are unreachable from the committed code. The per-module gene
sets they need (`genesets/A_sender_for_*.txt`) are present and tracked; only the resolution
logic is missing. `run_phase8_compmatch.py` carries its own per-module handling and could
be the donor, but that script is itself broken by B1.

Consequently the producer command the project documents for this file —
`WRITING_PACK.md:378,1449`, `run_phase3_nulls.py --stage main --sections all --calls all9` —
also fails: `all9` is not `all`, so it is split to `["all9"]` and hits the same raise.

### B3. `caller_disagree_all.py` — the reconstructed caller-agreement producer — does not run

`code/caller_disagree_all.py:62` calls `CD.DS_ALIAS.get(tag, tag)`. `code/caller_disagree.py`
has no `DS_ALIAS` (its top-level names are `EXCL, OUT, PROC, RAW, SAMP, run, …`). The file's
own docstring, line 12, says: *"The only edit made to `caller_disagree.py` is `DS_ALIAS`."*
**That edit was never committed.** `_load()` is on the path of both `--verify` and `--all`,
so every invocation raises `AttributeError`.

This is a direct extension of the previously-flagged defect: the reconstruction that closed
the two-producerless-tables finding is itself not runnable against the committed tree.

Unlike B1 and B2, **B3 is unambiguously repairable in one line**: `caller_disagree.py:22` has
`SAMP={'sham':'7250_liver_sham_Male_26-U1','sbr':'7259_liver_sbr_Male_26-U1'}`, and
`deepscence_reanchor.py:62` already carries the exact inverse the caller needs,
`DS_ALIAS = {'7250_liver_sham_Male_26-U1': 'sham', '7259_liver_sbr_Male_26-U1': 'sbr'}`.
(Not applied — this audit is read-only.)

### B4. A first-attempt re-run silently destroys the frozen headline tables and exits 0

`run_phase3_nulls.py:527` filters the section list to those with a `cache3` npz:

```python
secs = [s for s in secs if os.path.exists(os.path.join(P.CACHE3, f"{s}.npz"))]
```

`data/processed/cache3/` is **gitignored** (correctly — 166 MB, regenerable). So on a fresh
clone `secs` is `[]`. There is no guard. `stage_main` then builds `jobs=[]`,
`df = pd.DataFrame([])`, and executes `df.to_csv(f"{RES}/main_fits.csv")` — verified to write
a 1-byte file containing a single newline — then prints `(0, 0)` and returns. Exit status 0.

The same path destroys `perm_nulls.csv`, `perm_curves.csv` (`stage_perm`) and `curves.csv`
(`stage_curves`). All four are tracked; `main_fits.csv` and `perm_nulls.csv` carry the headline
and two of the three checksums in the tag annotation.

This is exactly the `build_genesets.py` defect — dead input path, bare fallthrough, empty
overwrite, exit 0 — on the headline files rather than the gene sets. Note that
`run_phase3_n8.py:252` handles the identical situation correctly:

```python
missing = [s for s in secs if not os.path.exists(...)]
if missing: raise SystemExit(f"no Phase 3 cache for: {missing}")
```

`run_phase8_compmatch.py:697` has a partial guard (`raise SystemExit("no cached sections for
this arm")`). `run_phase3_nulls.py` — the one that matters most — has none.

Secondary, same file: `main_fits.csv` is written with **no tag and no merge**. Re-running
`--stage main` with any subset of calls or sections replaces the nine-call table with a
smaller one, silently.

### B5. Seeds are a function of argument order

`stage_main:302` seeds `P.MASTER_SEED + 1000*i + j` and `stage_perm:461` seeds
`P.MASTER_SEED + 5000*i + 17*j`, where `i` indexes the **position** in the section list and
`j` the **position** in the `--calls` list. Running a subset of calls, or the same calls in a
different order, produces different seeds and therefore different permutation nulls — with no
warning. Section order is safe (`P.IN_BAND` etc. are hardcoded in `sasp_phase3.py:67`), but
call order is free-form on the command line and is not recorded for `--stage main`.

Related: `stage_perm`'s `do_full=(calls == [PRIMARY_CALL])` silently switches off the
full-scope N1 null and changes the output filename to `_n7` based on the argument list.

### B6. No committed script launches the head of the M1 pipeline

Every tracked M1 driver **waits for** the primary stages rather than starting them:
`_m1_chain.sh`, `_m1_pipeline_driver.sh`, `_m1_driver_a.sh`, `_m1_driver_b.sh` block on
`pgrep`; `_m1_final.sh` blocks on `grep` of `logs/m1_*.log`; `_m1_driver_b2.sh` blocks on a
PID passed as `$1`. The exact invocation of `run_phase3_nulls.py --stage main` — its
`--calls`, `--n-perm`, `--n-jobs` — appears in no committed file. (`--stage perm`/`perm_c1`
args are partially recoverable from `_m1_driver_b.sh`'s pgrep pattern.)

On a clean machine `pgrep` matches nothing and `grep` on the absent, gitignored
`logs/m1_chain.log` errors — so **the wait loops fall through immediately** and stages 2–5
run against absent or stale inputs. `_m1_final.sh` additionally hardcodes
`START=1787812466`, an epoch stamp from the original run, as its file-freshness threshold.

The real invocations survive only in `logs/` — 225 files, gitignored, on a workspace that
has been wiped twice. `logs/m1_perm_c1.log` line 1 confirms the `perm_c1` stage existed and
ran at `perms=1000` over the six in-band sections.

### B7. Files with no committed producer at all

Confirmed by a full sweep of the 234 concrete `results/` paths cited by `reports/*.md` and
`README.md` (~205 have a tracked producer; 23 do not; **0 have an untracked producer** — the
only untracked code is `code/h1_*`).

| File | Carries | Status |
|---|---|---|
| `results/phase3/caller_*_2sec_c6.csv` (7 tables) | the **FROZEN** post-C6 base of `summarize_caller_coverage.py:67`, hence `caller_coverage_gate_headline.csv`, hence `README.md:299` and *"independence FALSIFIED: 1.030 → 1.118, p=1.44e-30"* | `caller_disagree_all.py` defines `run_set(specs, suffix)` but `__main__` only exposes `--verify` and `--all`; nothing passes `suffix='_2sec_c6'`. Produced by an ad-hoc call that was never committed — **and would not run today anyway (B3)** |
| `results/phase3/m1_final_audit.txt` | **`power80_bound` = 0.1833** and the controlled amplitude 0.0288 (`WRITING_PACK.md:167,409,412`), `lam_railed_frac`, the §3 frozen vector | `m1_final_audit.py` is tracked but prints to stdout only. The redirect is not in the repo, and no committed driver invokes the script |
| `results/phase3/m1_n7_prepost.txt` | the N7 pre/post comparison (`CORRECTIONS.md:725`) | header and column block appear in no tracked file |
| `results/phase3/m1_prepost_main_fits.txt` | per-module reportable fits (`CORRECTIONS.md:376`) | stdout of `m1_compare_modules.py`; no committed redirect |
| `results/phase3/figure2b_data.csv`, `figure2d_data.csv` | Figure 2b/2d plotted values | `make_figure2bc.py` writes **only** `figure2c_data.csv` (line 156). **`CS_PHASE8_M1_RERUN.md:128` states it writes 2b and 2d; it does not.** That extension was never committed |
| `results/phase3/figure2e_data.csv` | every plotted number of Figure 2h (`CS_PHASE8_CALLERS.md:50,429`) | `figure2e` appears once in tracked code, in a comment. No `fig2e()` exists, though `CS_PHASE8_C1_CLOSEOUT.md:305` cites its docstring |
| `results/phase3/sf_summary_c1_swap_vs_n1.csv` | the **N3-swap ≡ N1** C1 verdict (N1 0.716, N3-swap 0.721, ρ 0.948) | `swap_vs_n1` occurs nowhere in `code/` |
| `results/phase4/cellchat_summary_statistic.csv` | *"in 3 of 18 for Tgfb1"* (`CS_PHASE4.md:304`) | the only one of 12 top-level `results/phase4/*.csv` with no writer |
| `results/section_qc_sender_summary.csv` | per-section depth/burden/portal-triad validity (`BIO_PHASE3.md:436`) | `section_qc` occurs nowhere in `code/` |
| `results/zonation_gene_correlations_7250_sham.csv` | per-gene zonation correlation, all 5,106 (`BIO_PHASE1.md:137`) | **a gene-set-defining input** — `build_genesets.py:166` reads it as `ZONATION_CSV` — with no committed derivation |
| `results/phase8_d2/dca_venv_pip_freeze.txt` | the DCA environment manifest | `setup_dca_env.sh` builds the venv but never dumps `pip freeze`. **Low severity** — the file's 83 pins are complete and tracked (see D3), so only its regeneration is uncaptured, not its content |
| `results/phase8_d2/committed_deepscence_sha256.txt`, `dca_worker_meta_*.json`, `d2_tables.md` | D2 provenance | `dca_denoise_worker.py` writes its meta to a tempdir; nothing copies it out. `report_d2_tables.py` prints to stdout |
| `results/a3_fallback/gpl33762_count.xml` | `<Count>132</Count>` (`A3_FALLBACK_SCREEN.md:60`) | no eutils call in `code/` |
| 14 `.log` files under `results/phase7_jobA/` and `results/a3_fallback/` | provenance only, no numbers | stdout redirects, no committed driver. Low severity |

### B8. The committed M1 pipeline is not a regeneration path

`_m1_final.sh` → `_m1_rerun_stage5.sh` never invokes `m1_final_audit.py`, `m1_headlines.py`,
`m1_compare_modules.py`, `caller_disagree_all.py`, or `run_a7_control_probes.py`, and emits no
figure2b/2d/2e data. Run end to end on a clean clone it would leave every artefact in B7
untouched — including the file holding 0.0288 and 0.1833 — while `summarize_caller_coverage.py`
(which *is* in stage 5) would consume an input nothing in the repo can rebuild.

---

## DEGRADES reproducibility

### D1. `figure2a` caches its own input, and the cache is a tracked, guarded artefact

`make_phase5_figs.py:411-421`:

```python
f = f"{FIG}/figure2a_stratified_curves.csv"
if os.path.exists(f):
    df = pd.read_csv(f)
else:
    df = build_fig2a_data(P.IN_BAND); df.to_csv(f, index=False)
fig2a(df)
```

`figures/figure2a_stratified_curves.csv` is **tracked** and is one of the 52 artefacts
`check_figures_guard.py` watches. So `git clone` + `make_phase5_figs.py --which 2a`
regenerates Figure 2a byte-identically from the committed cache regardless of any upstream
change, and the guard then reports OK. **The guard's pass on figure2a is vacuous.** There is
no `--force`. The sole mitigation is one `rm -f` line in `_m1_rerun_stage5.sh:22`; anyone
invoking the figure script directly gets the stale figure silently.

Same pattern, lower stakes:
- `corescence_circularity.load()` (line 153) returns the tracked
  `results/phase7_jobA/corescence_circularity_mouse.json` if present rather than re-deriving.
  `gate_disjointness_human.py:42` imports it, so **the human gene-set gate compares against a
  committed cache, never a fresh derivation.**
- `sasp_sweep.py:141` skips any `(sweep, cfg)` whose checkpoint exists — and all 43
  `results/sweep/*.csv` are tracked, so a default `python3 sasp_sweep.py` on a fresh clone
  recomputes **nothing** and just re-concatenates. Mitigated: `--force` exists and
  `code/README.md:27` documents deleting the checkpoints. This is why the Phase 1
  re-verification needed a separate `results/repro_2026-08-21/`.
- `prepare_samples.py:27`, `phase4_data.py:42`, `run_deepscence_all.py:38`,
  `run_deepscence_dca.py:55`, `run_deepscence_denoise_probe.py:173` all skip when the output
  exists. Their outputs are gitignored, so a fresh clone is safe; a partially populated tree
  is not.

### D2. Eight tracked scripts read package data from the container overlay, not the project env

```
/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv
/usr/local/lib/python3.11/dist-packages/senepy/data/
```

hardcoded in `caller_disagree2.py:22`, `corescence_circularity.py:49`,
`deepscence_reanchor.py:49`, `deepscence_smoke.py:21`, `gate_disjointness_human.py:47`,
`run_phase3_n8.py:39`, `senepy_coverage_candidates.py:21`, `senepy_coverage_human.py:21`.

This is precisely what master plan §16.1 forbids in bold — the overlay is what got wiped
twice. The identical file exists inside the project env at
`/workspace/envs/sasp311/lib/python3.11/site-packages/DeepScence/data/coreGS_v2.csv`, and no
script points there. A clone plus `pip install -r requirements.txt` into a venv — the
procedure `README.md:446` documents — does not populate `/usr/local/lib/python3.11/dist-packages`
at all.

Affected outputs include `run_phase3_n8.py`'s N8 control (the scrambled-response SF 0.916 in
`WRITING_PACK.md`) and `gate_disjointness_human.py`, which is the human half of
`gate_genesets_guard.py` and therefore of the pre-commit hook and of the gate the freeze
claims passed.

### D3. The DCA venv is on the container overlay, in a session-scoped scratchpad

`code/setup_dca_env.sh` is genuinely good: `DCA_ENV_ROOT` is a required parameter, the
CPython 3.8.19 download is SHA-256 pinned, versions are exact, the PyYAML 5.4.1 deviation is
explained, and `_shims_dca_bridge/dca/api.py:28` fails loudly if `DCA_VENV_PYTHON` is unset.

But four tracked runner scripts — `run_d2_dca.sh:3-4`, `run_d2_dca_stream.sh:6`,
`run_d2_stability.sh:15`, `run_d2_stream.sh:31` — hardcode

```
/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt/v38/bin/python
```

a per-session path under `/tmp`. The venv is there now (CPython 3.8.19 confirmed) but it is
on the overlay and outside `/workspace`, so it does not survive a pod restart, and the path
is meaningless to anyone else. `setup_dca_env.sh` itself is correct — it *requires*
`DCA_ENV_ROOT` — so the fix is to make the runners honour it too.

**The isolation is documented well enough to reconstruct.** `results/phase8_d2/dca_venv_pip_freeze.txt`
is tracked, in the tag, and complete: 83 fully pinned packages including `DCA==0.3.4`,
`tensorflow==2.4.4`, `Keras==2.4.3`, `PyYAML==5.4.1`. Together with `setup_dca_env.sh`'s
SHA-256-pinned CPython 3.8.19 download and the `_shims_dca_bridge` subprocess protocol, an
independent party has what they need. The only gap is that no committed script *writes* that
manifest (B7) and the runners hardcode a dead path.

### D3b. `torch` is imported directly but not pinned

`run_deepscence_all.py:25`, `run_deepscence_dca.py:37` and `run_deepscence_denoise_probe.py:56`
all call `torch.set_num_threads(...)` at import time. `torch` appears nowhere in
`requirements.txt`. It is present as `2.4.1+cu124` in both the venv and the overlay, and would
arrive transitively via DeepScence, but unpinned — and DeepScence scores are a function of the
torch version. Same defect class as the `kneed` and `openpyxl` omissions that have already
bitten.

### D4. `figures/.committed_manifest.json` is gitignored

`.gitignore:93` excludes it; it is untracked. `check_figures_guard.py:62-63` exits
`'no manifest; run --snapshot first'` without it, and `--snapshot` on a fresh clone records
whatever is present as correct — a vacuous baseline. **The tag annotation's
"figure guard 52/52 committed artefacts match" cannot be re-verified from the tag.**
(All 52 figure artefacts themselves are tracked and the guard's walk covers them, so only the
expected-hash file is missing. Contrast `genesets/.geneset_manifest.json`, which *is* tracked,
and whose guard only warns when it is absent rather than refusing to run.)

### D5. `data/interim/` is empty and two scripts read from it

`annotate_celltypes.py:96` writes interim `.h5ad`; `reannotate.py:20` and `compact_h5ad.py:7`
read them. `*.h5ad` is gitignored and `data/interim/` is empty on disk, so `reannotate.py`
cannot run until `annotate_celltypes.py` has. This ordering is not documented anywhere.

### D6. `caller_disagree*.py` silently drop DeepScence when its input is absent

`caller_disagree.py:30`, `caller_disagree2.py:42`, `caller_disagree_all.py:63`,
`phase2_downstream.py:97` all wrap the DeepScence join in `if os.path.exists(ds):`.
`data/processed/deepscence_*.csv` is gitignored, so on a fresh clone the caller-agreement
tables regenerate **without DeepScence**, silently, with a different caller set, exit 0.
Similarly `analyse_d2_stability.py:36` computes pairwise correlations over only whichever
run files happen to exist, so the D2 stability table's row count depends on disk state.

### D7. Silent partial figures

`make_figure2bc.py`'s `__main__` (lines 191-196) runs each panel only `if os.path.exists(...)`,
and `fig2c` merges the null columns only `if not pf.empty` then skips absent columns with
`if col not in rep: continue` (line 130). Without `perm_nulls.csv` it emits a plausible-looking
Figure 2c with fewer null rows, plus a correspondingly short `figure2c_data.csv`, exit 0.
`m1_headlines.py` likewise omits ~8 headline keys rather than failing when inputs are absent.

### D8. Documented environment setup conflicts with the actual environment

- Master plan §16.1 prescribes conda/mamba at `/workspace/envs/sasp` with an **unpinned**
  conda-forge list including `squidpy`, `spatialdata`, `spatialdata-io`, `geopandas`,
  `shapely`, `jupyterlab`, plus `pip install formulaic`. The real environment is a **pip
  venv** at `/workspace/envs/sasp311`, and `requirements.txt` pins none of those and states
  `formulaic: not installed`. §16.2's Dockerfile/`env.yaml` repeats the conda recipe.
- `README.md:446` gives only `pip install -r requirements.txt` — no venv path, no note that
  the shell's `python3` is `/usr/bin/python3` (3.11.10) and **not** the project environment
  (`VIRTUAL_ENV` is unset), and no mention of the separate DCA environment. Every driver
  script calls bare `python3`.
- `README.md:440` documents `python3 code/build_genesets.py` as the gene-set rebuild. That
  script now (correctly) refuses to run whenever `genesets/README_C6_mouse_provenance.md`
  exists — which it does, and is tracked. The README instruction exits 1 on a fresh clone.
  The failure is loud and the message is good; the instruction is simply stale.

---

## COSMETIC

- `core.hooksPath = .githooks` is repo-local git config, not part of the tree. A clone gets
  `.githooks/pre-commit` but the hook does not fire until the documented
  `git config core.hooksPath .githooks` is run. The hook's own comment says so.
- The promotion of `genesets/mouse_c6/*.txt` into `genesets/*.txt` is unscripted, but both
  sets are tracked and verified byte-identical (15/15), and `build_genesets_mouse_c6.py`
  regenerates `mouse_c6/` from the in-repo MSigDB pin. Reproducible in effect.
- `README.md:184` still asserts a composition-only surrogate *"reproduces **76%** of the
  contact amplitude"*. Commit `fb23a7e` withdrew that number as unsourced, but only in
  `CS_PHASE5.md`. Claim consistency rather than reproducibility; flagging in passing.
- `results/phase3_pre_c6/` (98 files) and `results/phase5_pre_c6/` (21) have no writer by
  design — they are tracked frozen snapshots, recoverable via the `pre-c6-genesets` tag as
  `CS_PHASE8_M1_RERUN.md:501` documents. Regenerating them needs a full re-run under the old
  gene sets, not a script. Worth stating because `summarize_caller_coverage.py`'s **PUBLISHED**
  row depends on them.

---

## What the tag actually captures

`phase8-frozen` **does** contain every cited result file, all 52 figure artefacts, all 172
gene-set files, both MSigDB pins, the 40 tracked `data/processed/` summary tables, the
drivers, and `requirements.txt`. The three annotation checksums verify.

Essential things **not** in it:

1. **The `run_phase3_nulls.py` revision that produced the results** (B1) — and with it the
   `tierApm_*` sender logic (B2) and the `DS_ALIAS` edit to `caller_disagree.py` (B3). These
   are not merely undocumented; they are absent from the repository's entire history.
2. **The head of the M1 pipeline** (B6) — no committed script starts stage 1.
3. **`figures/.committed_manifest.json`** (D4) — the tag's own figure-guard claim is
   unverifiable from the tag.
4. **`logs/`** — 225 files, gitignored, and the only surviving record of the actual
   invocations, including proof that `perm_c1` ran at 1,000 permutations.
5. **The DCA venv** (D3) — reconstructable from `setup_dca_env.sh`, but the tracked runners
   point at a dead session path.
6. **Twelve commits of corrections made after the tag was cut** — including a withdrawn
   claim, a corrected p-value bound, and `RECORD_RECONCILIATION.md`. `origin/main` is at
   `926439…`, the tag commit, so **none of the post-freeze corrections have been pushed.**
   They exist only in this workspace, which has been wiped twice. Reports inside the tag
   still carry numbers that the reconciliation subsequently corrected.

## Things that are done well

Worth recording so they are not lost in a later cleanup: `build_genesets.py`'s two hazard
guards; `build_genesets_mouse_c6.py` reading the archived pin and reporting the three HTML
error pages instead of swallowing them; `run_phase3_n8.py:252`'s loud missing-cache check;
`run_phase8_compmatch.py`'s `SASP_H1_UNFROZEN` freeze guard; `_shims_dca_bridge`'s loud
failure on an unset interpreter; `setup_dca_env.sh`'s checksum-pinned interpreter download;
`gate_genesets_guard.py` warning rather than refusing when its manifest is absent;
`sasp_sweep.py`'s stable `cfg_id` seeding; and the tracked `.gitignore` rationale comments,
which are unusually good.

---

## Method, and one side effect

Static analysis plus cheap checks: `git ls-files` / `git ls-tree` inventories; a grep of every
`results/` path cited by `reports/*.md` and `README.md` against every write call in
`git ls-files code/`; an AST-based cross-module attribute check over all 158 tracked scripts;
an AST-based check of every CLI flag in every tracked `.sh` driver against the target script's
`argparse`; existence checks on all 139 hardcoded `/workspace/...` paths and all quoted
relative paths; and argument-parse-only invocations of `run_phase3_nulls.py` against
non-existent section names (which compute nothing) to confirm B1 by execution rather than by
reading. No pipeline stage was run. No package was installed.

**Side effect, disclosed:** confirming B1/B2/B3 required importing `run_phase3_nulls`,
`sasp_phase3` and `caller_disagree`. Python rewrote their `code/__pycache__/*.pyc` from the
committed sources. `__pycache__/` is gitignored and regenerable and no source file was touched,
but if a stale `.pyc` from the original Phase 8 run had been sitting there it would have been a
recovery route, and it is now overwritten. Worth knowing before the next such audit.
