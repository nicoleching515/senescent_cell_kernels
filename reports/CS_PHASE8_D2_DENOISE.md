# Phase 8, task 8.5 — C7/D2: resolving `denoise=False`

**Owner:** CS. **Date:** 2026-08-27. **Spec:** `Phase7_Minimal_Human_Replication (1).md`
§4 (D-b), §6 (D2), §9. **Status of 8.5:** resolved.

---

## The two answers, first

### 1. Did DCA install? **YES — §6 path 1 landed, and DeepScence has been run at its published default on real M1 data.**

DCA 0.3.4 runs under TensorFlow 2.4.4 / Keras 2.4.3 in an isolated Python 3.8.19 venv, out
of process, handing a denoised matrix back to the pinned Python 3.11 stack over plain
`.npz`/`.npy` files. **The main environment was not modified and has no TensorFlow in it.**
`denoise=True` completed on **three full M1 sections** — 7239 and 7259 on the SBR arm
(7259 is the SBR half of the preserved two-section base) and 7352 on the sham arm — plus a
20,000-cell subsample panel at three seeds. §22 rated this Medium-High to fail; it did not fail.
Recipe: `code/setup_dca_env.sh`. Rebuild takes about four minutes.

### 2. How much does `denoise=False` cost? **A lot — and in the opposite direction to the one §4 (D-b) assumes.**

Committed `denoise=False` vs published `denoise=True`, same cells, same seed, same anchor,
same panel:

| | 7259 (SBR 26 wk) | 7352 (**sham** 2 wk) | 7239 (SBR 52 wk) |
|---|---|---|---|
| cells compared | 114,721 | 122,380 | 75,384 |
| median transcript counts | 491 | 710 | 911 |
| Pearson *r* between the two score vectors | **0.614** | **0.671** | **0.718** |
| Spearman ρ | 0.621 | 0.673 | 0.693 |
| Global top-5% sender set, Jaccard | **0.141** | **0.118** | **0.280** |
| cells changing global top-5% status | 8,634 | 9,646 | 4,238 |
| % of the committed sender set not called | **75.3%** | **78.8%** | **56.2%** |
| Within-cell-type top-5%, Jaccard | 0.073 | 0.090 | 0.132 |
| Cell-type × depth-decile matched top-5%, Jaccard | 0.092 | 0.077 | 0.142 |
| Global top-1%, Jaccard | 0.041 | 0.058 | 0.157 |
| **ρ(score, transcript_counts)** | 0.318 → **0.531** (**+0.214**) | 0.410 → **0.542** (**+0.132**) | 0.389 → **0.640** (**+0.251**) |
| ρ(signed score shift, depth) | +0.337 | +0.256 | +0.313 |
| hepatocyte % of the top-5% calls | 64.5% → **100.0%** | 97.2% → **100.0%** | 71.5% → **100.0%** |
| biliary/ductular enrichment | 0.99 → **0.001** | 0.16 → **0.000** | 1.20 → **0.000** |
| DeepScence's own gene-set metric (`corr_metrics`, chosen node) | 0.126 → **0.470** | — → **0.673** | 0.158 → **0.497** |

**§4 (D-b)'s premise is wrong and should be corrected in the paper.** It says *"DCA
denoising is precisely the step that would normalize depth — the confound under
investigation."* It is not. Switching denoising **on** raises the score's correlation with
sequencing depth on **all three sections, both arms** (+0.13, +0.21, +0.25 in Spearman ρ,
across a 1.9× spread in median depth), pushes deep cells up and shallow cells down
(ρ(Δz, depth) +0.26 to +0.34), and collapses the global top-5% call set onto hepatocytes
**alone — 100.0% on every section**, with every other cell type, biliary/ductular included,
dropping to zero or near it. At the same time the denoised score is a **better** senescence score by
DeepScence's own internal criterion (mean |correlation| with the CoreScence genes, 0.156 →
0.497, and 0.126 → 0.470 on 7259; the `denoise=False` figures come from the `raw` control
runs, and 7352 has no `raw` run so only its denoised value, 0.673, is quoted). **The
published default is simultaneously more faithful to its gene set and more confounded with
depth**, which is the uncomfortable shape of this result: the step that makes the score
look more like a senescence score is the same step that makes it look more like a depth
score.

A small piece of corroboration from the other direction: removing depth (`lib`, 7259) moves
that same internal metric slightly *down*, 0.126 → 0.121. Part of DeepScence's agreement
with its own gene set on this panel is depth.

### 3. And a third answer nobody asked for: **the published default is not reliably reproducible.**

One fixed 20,000-cell subsample, three seeds, changing nothing else (DeepScence hands its
`random_state` straight to `dca()`, so one seed moves both stages):

| pair | Pearson *r* | top-5% Jaccard | cells changing status |
|---|---|---|---|
| `denoise=False`, seed 0 vs seed 1 | **0.9955** | **0.761** | 272 |
| `denoise=True`, seed 0 vs seed 2 | 0.9824 | 0.665 | 402 |
| `denoise=True`, seed 0 vs seed 1 | **0.570** | **0.000** | 2,000 |
| `denoise=True`, seed 1 vs seed 2 | **0.573** | **0.000** | 2,000 |

Seeds 0 and 2 land in the same place. **Seed 1 lands somewhere else entirely** — its top-5%
sender set is *completely disjoint* from both others, Jaccard exactly 0. Nothing about the
seed-1 run looks broken: its denoised matrix has the same global statistics as seed 0's
(mean 0.2628 vs 0.2631, max 513 vs 486), and it produced a normal-looking score
distribution. Its internal gene-set metric does hint at what happened: at seeds 0 and 2 one
bottleneck node dominates (0.464 and 0.463, against 0.179 and 0.152 for the other), while at
seed 1 the two nodes are nearly tied (0.373 vs 0.335). The denoised representation at seed 1
is not one where a single senescence axis stands out — but the code still picks a node and
returns a score, with no diagnostic that anything is different.

State it precisely, because it is a strong claim: **one of three seeds at the published
default produced a sender set sharing no cells with the other two**, while ~~three
`denoise=False` runs agreed at Jaccard 0.76–0.99~~ **the two `denoise=False` runs agreed at
Pearson r = 0.996 and top-5 % Jaccard 0.76**.

***[Corrected 2026-08-27 (remediation pass) — two errors in one clause, and this is the site
`WRITING_PACK.md` §5.9 inherited them from. **(1) "0.76–0.99" merges a correlation and a set
overlap into one range.** `results/phase8_d2/d2_stability.csv` contains exactly **one**
`denoise=False` seed pair — `raw_seed0 vs raw_seed1`, n = 20,000, `pearson_r` **0.99553**,
`top5_jaccard` **0.7606**, `top5_n_changed` 272 — and **no Jaccard of 0.99 appears anywhere in the
file**: the 0.99 is the Pearson r. **(2) "three runs" is wrong — there are two.** The only
`denoise=False` seeds on disk are `runmeta_raw_sub20000_7239_liver_sbr_Male_52-U1.json` and
`runmeta_raw_sub20000_seed1_7239_liver_sbr_Male_52-U1.json`; there is no `raw_seed2`
(`grep -rho "raw_seed[0-9]" results/phase8_d2/ | sort -u` → `raw_seed0`, `raw_seed1`).
`PREREG_PHASE8.md` P26 states it correctly ("r 0.9955 / Jaccard 0.761"). **The finding is
unaffected and if anything sharper**: the `denoise=True` outlier's Jaccard of 0.000 is being
compared against a two-seed floor of 0.76, not a three-seed range.
`AUDIT_NUMBERS_FINAL.md` R1.]***

The seed-0 numbers in §2 are therefore
the *favourable* case for `denoise=True` — seed 2 reproduces them almost exactly
(`dca_seed2` vs `raw_seed0`: *r* 0.724, J 0.324, against seed 0's 0.724 / 0.320) — and the
unfavourable case is that the caller does not converge to a single answer at all.

### 4. Removing depth ourselves moves the score too — in the **opposite** direction.

§6 path 3 was run at two strengths, and the pair is the cleanest evidence in this report.

| configuration | log-depth variance removed | Pearson *r* vs committed | top-5% Jaccard | **Δρ(score, depth)** |
|---|---|---|---|---|
| `raw` (control, 2 sections) | 0% | 0.9999991 / 0.9999999 | 0.9937 / 0.9997 | −0.0006 / −0.0000 |
| `mor`, median-of-ratios (4 sections) | **11–24%** | 0.9979–0.9997 | 0.876–0.960 | −0.009 … +0.005 |
| `lib`, exact library size (2 sections) | **100%** (p90/p10 depth ratio 11.0 → 1.000) | **0.843** / **−0.017** | **0.366** / **0.000** | **−0.236 / −0.381** |
| `dca`, published default (3 sections) | n/a — imputation, not rescaling | 0.61–0.72 | 0.118–0.280 | **+0.132 … +0.251** |

So the `mor` null is not evidence that normalisation cannot move this caller. It is
evidence that **the estimator §6 names is too weak on a sparse 5K panel to test the
question** — *poscounts* conditions on detection and throws away most of the depth signal,
leaving the p90/p10 depth ratio at ~9–10× instead of 11×. Take depth out completely and
the score moves a great deal on both sections: 46.4% (7259) and 100% (7352) of the committed
sender set stops being called, and the caller's depth loading falls from ρ = 0.318 to
**0.082** and from ρ = 0.410 to **0.029** — reductions of 74% and 93%. The call set also
*broadens* on both: on 7259 hepatocyte enrichment 2.46 → 1.63, biliary/ductular 0.99 →
**1.53**, and every minor type roughly doubles (stellate 0.09 → 0.24, Kupffer 0.08 → 0.16,
LSEC 0.25 → 0.51, T/NK 0.14 → 0.34, B 0.16 → 0.44); 7352 does the same, more sharply
(LSEC 0.07 → 0.45, stellate 0.07 → 0.48, Kupffer 0.03 → 0.30). §4 carries the caveat that
belongs on 7352's *r* = −0.017.

**Put the two together and the direction is unambiguous.** Removing depth de-confounds the
caller and spreads its calls across cell types. The published denoising step does the exact
reverse: it raises the depth loading on every section tested and concentrates the calls on
the deepest cell type.
Whatever DCA is contributing to DeepScence on this data, it is not depth normalisation.

### What this means for the freeze

`denoise` is not a free parameter that can be left at `False` with a footnote. It changes
which cells are senders — 75–79% of them on two of the three sections at the global
threshold, and 77–86% within cell type — and it changes the caller's technical loading,
which is the property the paper's argument turns on. **Recommendation: freeze `denoise=False` as primary and report
`denoise=True` as the published-default sensitivity**, for the reasons in §6 — chiefly that
the default is both more depth-loaded and seed-unstable. That is now a *chosen* value with
a stated reason, not a forced one, which is what §9 requires.

---

## 1. What was run, in §9's four attributes

Every configuration shares: **denoise state** as marked; **anchor** = published `CDKN1A`;
**panel** = ortholog-mapped, MGI `HOM_MouseHumanSequence` 1:1 pairs, 4,845 of 5,097 mouse
panel genes mapped onto 4,827 distinct human symbols; `random_state=0` unless stated;
≥20 raw counts/cell.

| tag | denoise | counts fed to DeepScence | **coverage** | file |
|---|---|---|---|---|
| *committed* | False | raw | **11/11 M1 sections** (task 8.3 / D1) | `data/processed/deepscence_<section>.csv`; `deepscence_{sham,sbr}.csv` for 7250/7259 |
| `dca` | **True**, real DCA 0.3.4 | raw; DCA denoises inside DeepScence | **3 sections** (7259, 7352, 7239 — both arms) + a 20k subsample at 3 seeds | `deepscence_dca_<section>.csv` |
| `mor` | False | median-of-ratios (DESeq2 *poscounts*) factors divided out | **4 sections** (7259, 7450, 7239, 7248) | `deepscence_mor_<section>.csv` |
| `lib` | False | library-size factors divided out — depth removed exactly | **2 sections** (7259 SBR, 7352 sham) | `deepscence_lib_<section>.csv` |
| `raw` | False | raw — **byte-identical configuration to the committed run** | **2 sections** (7239, 7259) + 20k subsample at 2 seeds | `deepscence_raw_<section>.csv` |

The four `mor` sections span median transcript counts **491 → 1138**, i.e. essentially the
whole M1 depth range (section medians run 446 → 968 panel-wide), and include **both**
sections of the preserved two-section base (7259 directly; 7250's `mor` run was OOM-killed
and 7248, the other 26-week sham, covers the deep end).

**Nothing committed was overwritten.** `deepscence_sham.csv`, `deepscence_sbr.csv` and the
nine `deepscence_<section>.csv` files are untouched. SHA-256 of all eleven was recorded at
the start of this task in `results/phase8_d2/committed_deepscence_sha256.txt` and
re-verified at the end — `sha256sum -c` returns OK on all eleven. Every new output carries
a configuration prefix.

**A panel fact for the Methods table.** DeepScence's CoreScence v2 set at its default
`occurrence >= 5` is 39 genes (23 up, 16 down). After the ortholog remap, **31 of those 39
are on the Xenium mouse panel** (`CDKN1A` among them; the run log prints "Using 31 genes in
the gene set for scoring"). Every DeepScence number in this project rests on those 31
genes.

---

## 2. The determinism control — the floor everything else is read against

`raw` is the same code, the same seed and the same data as the committed run:

| | Pearson *r* | top-5% Jaccard | cells changed | Δρ(score, depth) |
|---|---|---|---|---|
| `raw` vs committed, **same seed**, 7239 | **0.9999991** | **0.9937** | 24 of 75,384 | −0.0006 |
| `raw` vs committed, **same seed**, 7259 | **0.9999999** | **0.9997** | 2 of 114,721 | −0.0000 |
| `raw` seed 0 vs seed 1 (20k subsample) | **0.9955** | **0.7606** | 272 of 20,000 | — |

So the pipeline is effectively deterministic at fixed seed, and its seed-to-seed spread is
*r* ≈ 0.996 / Jaccard ≈ 0.76. Those are the two floors. Read against them:

* `mor` (*r* 0.998–0.9997, J 0.876–0.960) sits **at or above the seed-to-seed floor** — no
  detectable effect.
* `dca` (*r* 0.61–0.72, J 0.12–0.28) sits **far below both floors** — a real, large effect.
* `lib` (*r* 0.843 and −0.017, J 0.366 and 0.000, Δρ −0.236 and −0.381) is likewise far
  below both floors — also real, also large, and in the opposite direction.
* `dca` **seed 1 vs seeds 0 and 2** (*r* ≈ 0.57, J 0.000) is *itself* below the
  `dca`-vs-committed comparison, which is the §3 finding.

These controls settle a second question in passing, and it is one the freeze needs. **7259's
committed score is `deepscence_sbr.csv`, written 2026-08-20, before two container
rebuilds** — and it reproduces today at *r* = 0.9999999 with 2 of 114,721 calls moving. The
preserved two-section base is not merely intact on disk; it is re-derivable from raw data
in the rebuilt environment. The same holds for 7239 from the 2026-08-27 D1 run.

One caution on reading these rows: `rho_signed_dz_vs_depth` for the `raw` comparisons
(−0.47 and −0.16 in `d2_depth.csv`) is the correlation of a shift whose mean magnitude is
0.001–0.0002 z-units. It is the direction of numerical noise, not a depth effect, and
should not be quoted.

---

## 3. §6 path 1 — the isolated environment: what failed first, and what worked

DCA 0.3.4 requires `tensorflow>=2.0,<2.5` and `keras>=2.4,<2.6`.

**(a) venv on the project's Python 3.11 — fails. This is the failure `requirements.txt` records.**

```
$ python3.11 -m venv v311 && ./v311/bin/pip install "dca==0.3.4"
ERROR: Could not find a version that satisfies the requirement tensorflow<2.5,>=2.0 (from dca)
       (from versions: 2.12.0rc0, ..., 2.21.0)
ERROR: No matching distribution found for tensorflow<2.5,>=2.0
```

The oldest TensorFlow with a cp311 wheel is 2.12. No overlap with DCA's range.

**(b) venv on the system Python 3.10 — fails twice over.**

```
$ python3.10 -m venv v310
... The virtual environment was not created successfully because ensurepip is not available.
$ python3.10 -m ensurepip --version
/usr/bin/python3.10: No module named ensurepip
```

and even with pip bootstrapped it could not have worked. From the PyPI JSON API:
`tensorflow 2.4.4` publishes exactly nine files, all `cp36`/`cp37`/`cp38`; `2.3.4` is
cp36/37/38; `2.0.4` is cp36/37. **No TensorFlow below 2.5 has a wheel for any Python above
3.8.** Ubuntu 22.04 has no `python3.8` package (`apt-cache policy python3.8` →
`Candidate: (none)`).

**(c) a downloaded, self-contained CPython 3.8 — works.** `python-build-standalone`
`cpython-3.8.19+20240726-x86_64-unknown-linux-gnu-install_only`, SHA-256 checked against
the release's own `.sha256` (`e81ea4dd…b2765`), unpacked into scratch, venv created inside
it. `pip install dca==0.3.4` then resolves cleanly and pulls TensorFlow 2.4.4, Keras 2.4.3,
numpy 1.19.5, anndata 0.7.8, scanpy 1.8.2. Full manifest (83 packages):
`results/phase8_d2/dca_venv_pip_freeze.txt`.

**One forced deviation from DCA's own dependency list.** `dca.api` imports `dca.train` →
`dca.hyper` → `kopt` unconditionally, and `kopt` 0.1.0 calls `yaml.load()` without a
`Loader`, which PyYAML 6 removed:

```
File ".../kopt/config.py", line 60, in <module>
    _config = yaml.load(open(_config_path))
TypeError: load() missing 1 required positional argument: 'Loader'
```

PyYAML is pinned to 5.4.1 **inside the isolated venv only**. It affects `kopt`'s config
parsing; `kopt` is the hyperparameter-search wrapper and is never reached on the `dca()`
path.

**§6 path 2 (container) was not attempted, because it cannot be on this box:** `docker`,
`podman`, `apptainer` and `singularity` are all absent, as are `conda`, `mamba` and
`micromamba`. Path 1 made it moot.

### 3.1 How the two environments talk

`DeepScence/api.py` does `from dca.api import dca` at import time and, with
`denoise=True`, calls `dca(adata, random_state=random_state)` expecting it to mutate
`adata.X` in place. `code/_shims_dca_bridge/dca/api.py` honours that contract while the
real DCA executes in the other interpreter:

* counts out via `scipy.sparse.save_npz`, denoised matrix back via `numpy.save`;
* **deliberately not `.h5ad`** — the on-disk format changed between anndata 0.7.8 (the
  newest anndata that co-installs with TF 2.4) and the pinned anndata 0.12.19;
* the worker calls `dca(A, random_state=seed)` and nothing else, so every DCA
  hyperparameter keeps the value DeepScence itself would have used (`mode='denoise'`,
  `ae_type='nb-conddisp'`, `hidden_size=(64,32,64)`, `epochs=300`, `batch_size=32`,
  `early_stop=15`, RMSprop);
* the bridge asserts the returned shape matches, and `run_deepscence_dca.py` refuses to
  start unless the bridge — not the raising stub in `code/_shims/dca` — is the module that
  was imported.

`code/_shims/dca/api.py`, the stub that raises, is **unchanged**. Nothing that does not
explicitly ask for the bridge can pick it up.

### 3.2 The main environment was not modified

Verified after the fact: `tensorflow` is **not importable** in the pinned Python 3.11
environment, and the pinned stack is at its pins — numpy 2.4.6, scipy 1.17.1, pandas 2.3.3,
anndata 0.12.19, scanpy 1.11.5, torch 2.4.1+cu124, h5py 3.16.0, scikit-learn 1.9.0. Nothing
was added to `/usr/local/lib/python3.11/dist-packages` while this task ran (`find -newermt`
returns nothing). `keras` 3.15.1 *is* present in the main environment, but it was installed
at 03:09 UTC, before this task started, it is not mine, and it is non-functional there in
any case (it imports `tensorflow`, which is absent).

### 3.3 What it costs to run

| section | cells | DCA denoise | whole `denoise=True` call | same section, `denoise=False` |
|---|---|---|---|---|
| 7239 | 83,392 | 8.6 min | 16.4 min | 5.2 min |
| 7259 | 127,386 | 53.3 min | ~112 min | ~30 min |
| 20k subsample | 20,000 | 17.7–36.4 min | — | ~2 min |

The 7259 and subsample figures are inflated by CPU contention with two other agents, not by
the section size alone. TensorFlow 2.4 could **not** use this box's GPU: it wants CUDA 11
(`libcudart.so.11.0`, `libcublas.so.11`, `libcudnn.so.8`), none of which is present, so it
falls back to CPU with `Skipping registering GPU devices`. A CUDA-11 image would make this
much cheaper if the run is ever repeated at scale.

---

## 4. §6 path 3 — measuring the cost by normalisation instead

Two normalisations were applied to the counts *before* DeepScence saw them, and the same
cells re-scored:

* **`mor`** — DESeq2 median-of-ratios, *poscounts* variant (geometric mean over cells where
  a gene is non-zero; per-cell median of ratios over the genes it detects), estimated on
  genes detected in ≥5% of cells, rescaled to median 1. The estimator §6 names.
* **`lib`** — per-cell library size divided out exactly, so every cell ends at the median
  mapped depth and DeepScence's internal size-factor offset collapses to 1. The upper bound
  on what any rescaling normalisation can do.

**A correction to the framing, which matters for reading these.** From
`DeepScence/io.py:normalize`, the published pipeline **already** depth-normalises: it calls
`sc.pp.normalize_total` over the full panel, then `log1p`, then `scale`, and it already
passes per-cell size factors into the ZINB decoder as an offset (`network.py`:
`mu = mu * sf`). `denoise=False` does **not** remove depth normalisation; it removes DCA's
ZINB *imputation*. `mor` and `lib` replace library-size normalisation with a
composition-robust one and an exact one; they do not add normalisation where there was
none.

**`mor` result: no detectable effect.** Four sections spanning median depth 491–1138,
Pearson *r* 0.9979–0.9997, top-5% Jaccard 0.876–0.960, |Δρ(score, depth)| ≤ 0.009 — all
inside the seed-to-seed floor of §2.

**But `mor` only removed 11–24% of the log-depth variance** (`d2_normalisation_strength.csv`):
the 90th/10th-percentile depth ratio falls from ~11× to only ~8.7–10.0×. That is the
*poscounts* estimator behaving as it must on a sparse 5K panel — it conditions on
detection, and a shallow cell's non-zero entries are mostly 1s, so the median of ratios
barely registers its depth. **The `mor` null is a null about a weak normalisation, not
about normalisation.**

**`lib` result: a large effect, in the de-confounding direction.** With depth removed
exactly (verified: post-normalisation sd of log depth = 0.0000, p90/p10 ratio = 1.000):

| | 7259 (SBR 26 wk) | 7352 (sham 2 wk) |
|---|---|---|
| Pearson *r* / Spearman ρ vs committed | 0.843 / 0.822 | **−0.017 / 0.013** |
| global top-5% Jaccard | 0.366 (46.4% of the committed sender set not called) | **0.000 (100%)** |
| within-type top-5% / matched top-5% Jaccard | 0.379 / 0.379 | 0.049 / 0.046 |
| **ρ(score, transcript_counts)** | 0.318 → **0.082** (**−0.236**) | 0.410 → **0.029** (**−0.381**) |
| hepatocyte enrichment in the calls | 2.458 → 1.628 | 1.598 → 1.405 |
| biliary/ductular enrichment | 0.992 → **1.528** | 0.156 → **0.340** |
| stellate / Kupffer / LSEC / T-NK / B enrichment | 0.093 / 0.076 / 0.247 / 0.136 / 0.155 → 0.238 / 0.157 / 0.506 / 0.342 / 0.440 | 0.071 / 0.033 / 0.071 / 0.091 / 0.114 → 0.478 / 0.300 / 0.453 / 0.365 / 0.400 |

Both sections agree on the two things that matter: **removing depth cuts DeepScence's depth
loading by 74% and 93%**, and it **spreads the calls off hepatocytes and onto every minor
type** — the opposite of what the published denoising step does on the same sections.

**7352 needs a caveat, and it is an interesting one.** There, `lib` did not merely shift
the score, it produced one essentially *uncorrelated* with the committed score (*r* =
−0.017) whose top-5% is perfectly disjoint from it. Two independent 5% sets would be
expected to overlap by ~306 cells, not 0, so this is a structured alternative, not noise.
Its `corr_metrics` shows why: the two bottleneck nodes come back nearly tied (0.150 vs
0.127), and the chosen node's 0.150 is **as good as the committed configuration's ~0.15–0.16**
by DeepScence's own criterion. Take depth away on this section and there are two
equally-defensible senescence axes; DeepScence picks one and reports it with no indication
that the other exists. That is the same failure mode as the seed-1 DCA run in §3, reached
from the opposite direction.

Two consequences for reading the numbers. First, `pearson_r = −0.017` is **not** a polarity
flip in the D3 sense — the analysis negates whenever *r* < 0, which is right for a genuine
anchor flip and arbitrary when *r* ≈ 0, so `d2_agreement.csv` carries a
`sign_alignment_ambiguous` column and it is `True` on exactly this row and no other. Second,
the 7352 `lib` Jaccards quantify "a completely different answer", not "a shifted answer",
and should be quoted that way.

---

## 5. Coverage, and why it is not 11/11 for every configuration

Three things bounded the sweep, and they should be recorded because the next agent will hit
them too.

1. **The container cgroup is 57.7 GiB** — `/sys/fs/cgroup/memory.max` = 61,999,996,928
   bytes — not the 251 GB that
   `free` reports for the host. DeepScence holds five dense `n_cells × 4,845` float32
   arrays of a section at once — the input, `original`, the `read_dataset` copy, the
   `normalize` copy and a `raw_counts` layer. An 83k-cell section was OOM-killed with 11 GB
   free; **~16 GB free is the working requirement** for a section of that size, more for
   the 200k+ ones. Read `memory.current` against `memory.max`; `free` is misleading inside
   this container and was the reason for the first two failures.
2. **Two other agents share that 57.7 GiB.** From ~06:40 the M1 re-run held ~30 GB across
   ten `run_phase3_n8.py` workers. Four of my jobs were OOM-killed in that window
   (`mor 7250`, `mor 7260`, `raw 7239` on its first attempt, and `dca 7259` on its first
   attempt). Everything was then serialised behind a memory guard rather than contend with
   a job on the freeze's critical path, which is why the wall-clock figures in §3.3 are
   what they are. `dca 7259` succeeded on the second attempt; `raw 7239` succeeded on the
   second attempt; `mor 7250` and `mor 7260` were not retried.
3. **Compute was deliberately re-aimed once DCA proved to work.** §6 says "stop when one
   works", which makes path 1 the deliverable and path 3 the control. Six `mor` sections
   were dropped so that a **second full `denoise=True` section**, the **three-seed
   stability panel** and the **`lib` bracket** could run instead. That was the right trade:
   `mor` is a null, and a fifth null section adds less than a second section of an effect
   that is real.

**Still running at hand-off.** `dca 7352` and `lib 7352` both completed and are included
above. Only `mor 7352` — a fifth section of a configuration that is already a null on four —
was still queued in `code/run_d2_stream.sh`. When it lands, re-run

```bash
python code/analyse_d2_denoise.py && python code/report_d2_tables.py
bash   code/run_d2_stability.sh && python code/analyse_d2_stability.py   # the seed panel > results/phase8_d2/d2_tables.md
```

and the tables in §7 update. **Nothing in §0 depends on it**: the conclusions rest on three
full `dca` sections spanning both arms, two `lib` sections spanning both arms, four `mor`
sections, two full determinism controls and the three-seed panel. If those jobs are unwanted, kill `run_d2_stream.sh`; they write
only to `data/processed/deepscence_{dca,lib,mor}_7352_*.csv` and
`results/phase8_d2/`, never to `results/phase3/` or `figures/`.

---

## 6. Recommendation for the freeze

1. **Report `denoise` as a chosen, frozen parameter**, folded into §9's "denoise state"
   attribute — but as a value that was *selected*, with the reason below, rather than one
   that was forced by an install failure. The install failure no longer exists.
2. **Freeze `denoise=False` as primary.** Not because DCA is unavailable — it is available
   now — but because the published default (a) **increases** the caller's technical
   loading on all three sections measured and on both arms (ρ with transcript counts
   0.32→0.53, 0.41→0.54, 0.39→0.64),
   (b) collapses the sender set onto a single cell type, and (c) **does not reliably
   reproduce across seeds** (one of three seeds gave a top-5% set disjoint from the other
   two, against Jaccard 0.76 for `denoise=False`). For a paper whose subject is whether
   senescence callers are measuring depth, adopting the more depth-loaded and less
   reproducible configuration as primary would need a positive argument, and there is none.
   Note what this does *not* license: `denoise=False` is not the de-confounded option
   either. The `lib` result says an explicit depth normalisation would cut the caller's
   depth loading from ρ 0.32 to 0.08 and change 46% of its sender calls. The frozen caller
   is depth-loaded; it is simply less depth-loaded than the published default.
3. **Report `denoise=True` alongside it as the published-default sensitivity**, which is
   what §6's last line asks for. §0 is that sensitivity.
4. **Correct §4 (D-b) in the paper.** "DCA denoising is precisely the step that would
   normalize depth" is now measured to be false in this data: denoising raises the depth
   loading on three of three sections, on both arms. The caveat as currently written misdescribes the deviation.
5. **Promote the seed instability to Results.** "A published senescence caller's default
   denoising step returned, from one of three random seeds on the same cells, a sender set
   sharing no cells with the other two" is concrete, verifiable and useful — the same
   family of finding as the D3 polarity flip (task 8.6), and it belongs in the same place.
   It also deserves a note to the DeepScence authors, together with the observation in §4
   that the denoising step raises rather than lowers the score's depth loading.
6. **Consider `lib` as a declared sensitivity too.** Two sections, both arms: "the caller's
   depth loading falls by 74% and 93% if depth is removed before scoring, and its sender set
   changes by 46% and 100%" is the most direct evidence in the project that DeepScence's
   sender set is substantially a depth ranking. It is not a frozen candidate — on 7352 it
   lands on an axis uncorrelated with the committed one — but that instability is itself
   the finding: strip depth out and the caller has no single answer.
7. **Carry the DCA environment forward to H1.** §8's free experiment — CoreScence runs
   natively on human data — is now runnable at the published default on both arms, which
   removes `denoise=False` from the list of things that could explain a mouse-only
   artifact. Budget CPU time, or a CUDA-11 image, accordingly.

---

## 7. Tables

Generated by `code/report_d2_tables.py` from `results/phase8_d2/d2_*.csv`; no number below
is hand-transcribed. `sec` is the section number; `median_transcript_counts` is the
panel-wide median over the caller cell set (cell types joined, `Low_quality`/`Unknown`
dropped), which is why it differs slightly from the section-wide median.

### A. Score agreement, per section (committed `denoise=False` vs alternative)

| config | sec  | median_transcript_counts | n_cells | pearson_r  | spearman_rho | anchor_sign_flipped | sign_alignment_ambiguous | mean_abs_dz | p95_abs_dz |
|--------|------|--------------------------|---------|------------|--------------|---------------------|--------------------------|-------------|------------|
| dca    | 7259 | 491                      | 114721  | 0.6139     | 0.6211       | False               | False                    | 0.7009      | 1.7284     |
| dca    | 7352 | 710                      | 122380  | 0.6712     | 0.6727       | False               | False                    | 0.6605      | 1.5087     |
| dca    | 7239 | 911                      | 75384   | 0.7183     | 0.693        | False               | False                    | 0.5924      | 1.4894     |
| lib    | 7259 | 491                      | 114721  | 0.8428     | 0.8225       | False               | False                    | 0.4534      | 1.1207     |
| lib    | 7352 | 710                      | 122380  | -0.017     | 0.0132       | True                | True                     | 1.1772      | 2.5732     |
| mor    | 7259 | 491                      | 114721  | 0.99929147 | 0.99927706   | False               | False                    | 0.0278      | 0.0783     |
| mor    | 7450 | 869                      | 81428   | 0.99968460 | 0.99961343   | False               | False                    | 0.0182      | 0.0522     |
| mor    | 7239 | 911                      | 75384   | 0.99921725 | 0.99923602   | False               | False                    | 0.0289      | 0.084      |
| mor    | 7248 | 1138                     | 195443  | 0.99791975 | 0.99832700   | False               | False                    | 0.0461      | 0.1377     |
| raw    | 7259 | 491                      | 114721  | 0.99999995 | 0.99999861   | False               | False                    | 0.0002      | 0.0005     |
| raw    | 7239 | 911                      | 75384   | 0.99999913 | 0.99999292   | False               | False                    | 0.0011      | 0.0024     |

### B. Sender-status change at the operative thresholds

| config | sec  | global_top5_n_called | global_top5_jaccard | global_top5_n_changed | global_top5_pct_of_called_changed | global_top1_jaccard | within_type_top5_jaccard | matched_top5_jaccard |
|--------|------|----------------------|---------------------|-----------------------|-----------------------------------|---------------------|--------------------------|----------------------|
| dca    | 7259 | 5736                 | 0.1412              | 8634                  | 75.26                             | 0.0408              | 0.0733                   | 0.0916               |
| dca    | 7352 | 6119                 | 0.1184              | 9646                  | 78.82                             | 0.0584              | 0.0898                   | 0.0768               |
| dca    | 7239 | 3770                 | 0.2804              | 4238                  | 56.21                             | 0.1573              | 0.1316                   | 0.142                |
| lib    | 7259 | 5736                 | 0.3664              | 5320                  | 46.37                             | 0.3372              | 0.3792                   | 0.379                |
| lib    | 7352 | 6119                 | 0                   | 12236                 | 100                               | 0                   | 0.0486                   | 0.0458               |
| mor    | 7259 | 5736                 | 0.9162              | 502                   | 4.38                              | 0.8851              | 0.9306                   | 0.9299               |
| mor    | 7450 | 4072                 | 0.96                | 166                   | 2.04                              | 0.9451              | 0.9528                   | 0.9459               |
| mor    | 7239 | 3770                 | 0.92                | 314                   | 4.16                              | 0.8733              | 0.9152                   | 0.9154               |
| mor    | 7248 | 9773                 | 0.8756              | 1296                  | 6.63                              | 0.8127              | 0.881                    | 0.881                |
| raw    | 7259 | 5736                 | 0.99970000          | 2                     | 0.02                              | 0.99830000          | 0.99900000               | 0.99970000           |
| raw    | 7239 | 3770                 | 0.99370000          | 24                    | 0.32                              | 0.99210000          | 0.99810000               | 0.99840000           |

### C. Depth dependence of the shift

| config | sec  | median_transcript_counts | rho_depth_committed | rho_depth_alt | delta_rho_depth | rho_signed_dz_vs_depth | rho_abs_dz_vs_depth |
|--------|------|--------------------------|---------------------|---------------|-----------------|------------------------|---------------------|
| dca    | 7259 | 491                      | 0.3176              | 0.5314        | 0.2138          | 0.3369                 | -0.0004             |
| dca    | 7352 | 710                      | 0.4096              | 0.5419        | 0.1323          | 0.2561                 | -0.2248             |
| dca    | 7239 | 911                      | 0.3891              | 0.6404        | 0.2512          | 0.3132                 | -0.1799             |
| lib    | 7259 | 491                      | 0.3176              | 0.0819        | -0.2357         | -0.3967                | -0.1539             |
| lib    | 7352 | 710                      | 0.4096              | 0.0289        | -0.3807         | -0.2918                | -0.0747             |
| mor    | 7259 | 491                      | 0.3176              | 0.3129        | -0.0047         | -0.1349                | 0.0682              |
| mor    | 7450 | 869                      | 0.5577              | 0.5591        | 0.0014          | 0.1517                 | -0.0762             |
| mor    | 7239 | 911                      | 0.3891              | 0.3804        | -0.0087         | -0.2022                | 0.0054              |
| mor    | 7248 | 1138                     | 0.4571              | 0.4616        | 0.0046          | 0.0799                 | -0.1335             |
| raw    | 7259 | 491                      | 0.3176              | 0.3176        | -0              | -0.1583                | -0.08               |
| raw    | 7239 | 911                      | 0.3891              | 0.3885        | -0.0006         | -0.4673                | -0.0209             |

### D. How much depth variation each configuration actually removed

| config | sec  | sd_log_depth_before | sd_log_depth_after | frac_log_depth_variance_removed | depth_p90_over_p10_before | depth_p90_over_p10_after |
|--------|------|---------------------|--------------------|---------------------------------|---------------------------|--------------------------|
| lib    | 7259 | 0.9249              | 0                  | 1                               | 11                        | 1                        |
| lib    | 7352 | 0.8127              | 0                  | 1                               | 7.674                     | 1                        |
| mor    | 7239 | 0.9257              | 0.8093             | 0.2357                          | 11.207                    | 8.667                    |
| mor    | 7248 | 0.9415              | 0.8603             | 0.165                           | 10.608                    | 8.768                    |
| mor    | 7259 | 0.9249              | 0.8705             | 0.1141                          | 11                        | 9.991                    |
| mor    | 7450 | 0.963               | 0.8808             | 0.1635                          | 11.534                    | 9.664                    |
| raw    | 7239 | 0.9257              | 0.9257             | 0                               | 11.207                    | 11.207                   |
| raw    | 7259 | 0.9249              | 0.9249             | 0                               | 11                        | 11                       |

### E. Per-configuration summary across sections

| config | n_sections | cells  | min_r      | median_r   | max_r      | min_J5     | median_J5  | max_J5     | median_drho_depth | max_abs_drho_depth |
|--------|------------|--------|------------|------------|------------|------------|------------|------------|-------------------|--------------------|
| dca    | 3          | 312485 | 0.6139     | 0.6712     | 0.7183     | 0.1184     | 0.1412     | 0.2804     | 0.2138            | 0.2512             |
| lib    | 2          | 237101 | -0.017     | 0.4129     | 0.8428     | 0          | 0.1832     | 0.3664     | -0.3082           | 0.3807             |
| mor    | 4          | 466976 | 0.99790000 | 0.99930000 | 0.99970000 | 0.8756     | 0.9181     | 0.96       | -0.0016           | 0.0087             |
| raw    | 2          | 190105 | 1          | 1          | 1          | 0.99370000 | 0.99670000 | 0.99970000 | -0.0003           | 0.0006             |

### F. Score shift by within-section depth decile (mean z-score change)


**raw**

| depth_decile | 7239    | 7259 |
|--------------|---------|------|
| 1            | 0.0006  | 0    |
| 2            | 0.0006  | 0    |
| 3            | 0.0005  | 0    |
| 4            | 0.0005  | 0    |
| 5            | 0.0004  | 0    |
| 6            | 0.0002  | 0    |
| 7            | -0.0001 | 0    |
| 8            | -0.0006 | 0    |
| 9            | -0.001  | 0    |
| 10           | -0.0011 | 0    |

**mor**

| depth_decile | 7239    | 7248    | 7259    | 7450    |
|--------------|---------|---------|---------|---------|
| 1            | 0.0087  | -0.0171 | 0.0055  | -0.0032 |
| 2            | 0.0097  | -0.005  | 0.0047  | -0.0015 |
| 3            | 0.0094  | -0.0006 | 0.0041  | -0.0014 |
| 4            | 0.0071  | -0.0006 | 0.0033  | -0.0012 |
| 5            | 0.0039  | 0.001   | 0.0026  | -0.0013 |
| 6            | 0.0001  | 0.0034  | 0.0013  | -0.0017 |
| 7            | -0.0041 | 0.0052  | -0.0022 | -0.0009 |
| 8            | -0.009  | 0.0054  | -0.0066 | 0.0009  |
| 9            | -0.011  | 0.0052  | -0.01   | 0.0039  |
| 10           | -0.0149 | 0.0031  | -0.0028 | 0.0065  |

**lib**

| depth_decile | 7259    | 7352    |
|--------------|---------|---------|
| 1            | 0.3862  | 0.7226  |
| 2            | 0.208   | 0.4024  |
| 3            | 0.1358  | 0.259   |
| 4            | 0.0929  | 0.1372  |
| 5            | 0.0582  | 0.0463  |
| 6            | 0.0063  | -0.044  |
| 7            | -0.0614 | -0.1335 |
| 8            | -0.1435 | -0.2607 |
| 9            | -0.2961 | -0.4452 |
| 10           | -0.3905 | -0.6892 |

**dca**

| depth_decile | 7239    | 7259    | 7352    |
|--------------|---------|---------|---------|
| 1            | -0.3038 | -0.2784 | -0.464  |
| 2            | -0.2444 | -0.2556 | -0.2967 |
| 3            | -0.1995 | -0.2028 | -0.1673 |
| 4            | -0.1935 | -0.1985 | -0.039  |
| 5            | -0.1185 | -0.1786 | 0.0486  |
| 6            | -0.0289 | -0.1331 | 0.1351  |
| 7            | 0.1282  | -0.0525 | 0.1793  |
| 8            | 0.2912  | 0.0838  | 0.2256  |
| 9            | 0.3432  | 0.4057  | 0.2285  |
| 10           | 0.3278  | 0.8126  | 0.1531  |

### G. Top-5% call rate by depth decile, committed vs alternative


**raw** (mean over 2 section(s))

| depth_decile | pct_called_committed | pct_called_alt | pct_flipped |
|--------------|----------------------|----------------|-------------|
| 1            | 4.14                 | 4.146          | 0.006       |
| 2            | 4.033                | 4.06           | 0.026       |
| 3            | 3.406                | 3.413          | 0.02        |
| 4            | 3.228                | 3.241          | 0.014       |
| 5            | 2.622                | 2.634          | 0.011       |
| 6            | 2.593                | 2.593          | 0           |
| 7            | 2.966                | 2.972          | 0.006       |
| 8            | 5.3                  | 5.28           | 0.033       |
| 9            | 10.216               | 10.196         | 0.02        |
| 10           | 11.499               | 11.468         | 0.031       |

**mor** (mean over 4 section(s))

| depth_decile | pct_called_committed | pct_called_alt | pct_flipped |
|--------------|----------------------|----------------|-------------|
| 1            | 2.911                | 2.88           | 0.266       |
| 2            | 3.108                | 3.094          | 0.298       |
| 3            | 2.826                | 2.896          | 0.306       |
| 4            | 3.024                | 3.062          | 0.33        |
| 5            | 3.006                | 3.082          | 0.325       |
| 6            | 3.772                | 3.794          | 0.378       |
| 7            | 4.728                | 4.767          | 0.389       |
| 8            | 6.894                | 6.847          | 0.492       |
| 9            | 9.808                | 9.695          | 0.692       |
| 10           | 9.934                | 9.894          | 0.826       |

**lib** (mean over 2 section(s))

| depth_decile | pct_called_committed | pct_called_alt | pct_flipped |
|--------------|----------------------|----------------|-------------|
| 1            | 4.04                 | 6.628          | 5.938       |
| 2            | 4.034                | 5.837          | 6.79        |
| 3            | 3.924                | 5.014          | 6.775       |
| 4            | 4.374                | 5.538          | 7.805       |
| 5            | 4.202                | 5.535          | 8.112       |
| 6            | 4.198                | 5.464          | 8.231       |
| 7            | 4.219                | 4.558          | 7.518       |
| 8            | 5.092                | 4.171          | 7.178       |
| 9            | 7.717                | 4.382          | 7.733       |
| 10           | 8.2                  | 2.845          | 7.103       |

**dca** (mean over 3 section(s))

| depth_decile | pct_called_committed | pct_called_alt | pct_flipped |
|--------------|----------------------|----------------|-------------|
| 1            | 3.488                | 0.594          | 3.634       |
| 2            | 3.935                | 1.403          | 4.369       |
| 3            | 3.805                | 1.665          | 4.42        |
| 4            | 3.994                | 2.029          | 4.768       |
| 5            | 3.682                | 2.328          | 4.573       |
| 6            | 3.719                | 2.782          | 4.893       |
| 7            | 4.017                | 3.507          | 5.406       |
| 8            | 5.579                | 5.991          | 7.656       |
| 9            | 8.699                | 11.397         | 12.299      |
| 10           | 9.089                | 18.332         | 18.098      |

### G2. Cell-type enrichment of the top-5% calls, raw (mean over 2 section(s))

| cell_type              | bg_pct | enrich_committed | enrich_alt | delta  |
|------------------------|--------|------------------|------------|--------|
| Hepatocytes            | 31.187 | 2.218            | 2.214      | -0.004 |
| Biliary/ductular       | 25.398 | 1.094            | 1.1        | 0.007  |
| Hepatic stellate cells | 12.344 | 0.07             | 0.07       | 0      |
| LSECs                  | 10.716 | 0.138            | 0.138      | 0      |
| Kupffer cells          | 10.276 | 0.052            | 0.052      | 0      |
| T/NK cells             | 4.252  | 0.104            | 0.106      | 0.002  |
| DC                     | 2.556  | 0.032            | 0.032      | 0      |
| vSMCs                  | 2.217  | 0.024            | 0.024      | 0      |
| B-cells                | 0.94   | 0.1              | 0.1        | 0      |
| Central venous LSECs   | 0.877  | 0.121            | 0.121      | 0      |
| Proliferating          | 0.754  | 1.286            | 1.286      | 0      |
| Mesothelial cells      | 0.062  | 0                | 0          | 0      |

### G2. Cell-type enrichment of the top-5% calls, mor (mean over 4 section(s))

| cell_type              | bg_pct | enrich_committed | enrich_alt | delta  |
|------------------------|--------|------------------|------------|--------|
| Hepatocytes            | 38.929 | 2.12             | 2.094      | -0.026 |
| Biliary/ductular       | 16.278 | 0.934            | 0.934      | -0     |
| LSECs                  | 13.508 | 0.11             | 0.101      | -0.009 |
| Hepatic stellate cells | 10.935 | 0.069            | 0.067      | -0.002 |
| Kupffer cells          | 10.508 | 0.054            | 0.053      | -0.001 |
| T/NK cells             | 3.642  | 0.095            | 0.089      | -0.005 |
| vSMCs                  | 2.854  | 0.078            | 0.081      | 0.003  |
| DC                     | 1.846  | 0.021            | 0.022      | 0.001  |
| Proliferating          | 1.642  | 0.912            | 0.835      | -0.077 |
| Central venous LSECs   | 1.068  | 0.11             | 0.105      | -0.005 |
| B-cells                | 0.887  | 0.131            | 0.127      | -0.004 |
| Mesothelial cells      | 0.651  | 0.03             | 0.03       | 0      |

### G2. Cell-type enrichment of the top-5% calls, lib (mean over 2 section(s))

| cell_type              | bg_pct | enrich_committed | enrich_alt | delta  |
|------------------------|--------|------------------|------------|--------|
| Hepatocytes            | 43.546 | 2.028            | 1.516      | -0.512 |
| Biliary/ductular       | 15.978 | 0.574            | 0.934      | 0.36   |
| LSECs                  | 12.204 | 0.159            | 0.48       | 0.32   |
| Hepatic stellate cells | 10.064 | 0.082            | 0.358      | 0.276  |
| Kupffer cells          | 9.439  | 0.054            | 0.228      | 0.174  |
| T/NK cells             | 4.004  | 0.114            | 0.354      | 0.24   |
| vSMCs                  | 3.378  | 0.053            | 0.14       | 0.087  |
| DC                     | 1.958  | 0.044            | 0.126      | 0.083  |
| B-cells                | 0.622  | 0.134            | 0.42       | 0.286  |
| Proliferating          | 0.497  | 1.138            | 0.44       | -0.698 |

### G2. Cell-type enrichment of the top-5% calls, dca (mean over 3 section(s))

| cell_type              | bg_pct | enrich_committed | enrich_alt | delta  |
|------------------------|--------|------------------|------------|--------|
| Hepatocytes            | 41.075 | 2.012            | 2.74       | 0.728  |
| Biliary/ductular       | 18.052 | 0.781            | 0          | -0.781 |
| LSECs                  | 11.513 | 0.116            | 0          | -0.116 |
| Hepatic stellate cells | 10.546 | 0.071            | 0          | -0.071 |
| Kupffer cells          | 9.305  | 0.046            | 0          | -0.046 |
| T/NK cells             | 3.671  | 0.099            | 0          | -0.099 |
| vSMCs                  | 2.798  | 0.038            | 0          | -0.038 |
| DC                     | 2.205  | 0.036            | 0.003      | -0.032 |
| Central venous LSECs   | 0.877  | 0.121            | 0          | -0.121 |
| B-cells                | 0.818  | 0.104            | 0          | -0.104 |
| Proliferating          | 0.638  | 0.884            | 0.01       | -0.874 |
| Mesothelial cells      | 0.062  | 0                | 0          | 0      |

### H. Run provenance

| config             | section | n_cells | n_genes | denoise | anchor           | node | reverse | minutes |
|--------------------|---------|---------|---------|---------|------------------|------|---------|---------|
| dca                | 7239    | 83392   | 4845    | True    | published_CDKN1A | 0    | True    | 16.4    |
| dca                | 7259    | 127386  | 4845    | True    | published_CDKN1A | 0    | True    | 111.8   |
| dca                | 7352    | 139378  | 4845    | True    | published_CDKN1A | 0    | True    | 16.35   |
| dca_sub20000       | 7239    | 20000   | 4845    | True    | published_CDKN1A | 0    | True    | 21.18   |
| dca_sub20000_seed1 | 7239    | 20000   | 4845    | True    | published_CDKN1A | 0    | False   | 47.58   |
| dca_sub20000_seed2 | 7239    | 20000   | 4845    | True    | published_CDKN1A | 0    | True    | 31.35   |
| lib                | 7259    | 127386  | 4845    | False   | published_CDKN1A | 0    | True    | 61.08   |
| lib                | 7352    | 139378  | 4845    | False   | published_CDKN1A | 0    | False   | 8.6     |
| mor                | 7239    | 83392   | 4845    | False   | published_CDKN1A | 0    | True    | 5.16    |
| mor                | 7248    | 224922  | 4845    | False   | published_CDKN1A | 0    | True    | 30.75   |
| mor                | 7259    | 127386  | 4845    | False   | published_CDKN1A | 0    | True    | 44.13   |
| mor                | 7450    | 93197   | 4845    | False   | published_CDKN1A | 0    | True    | 17.85   |
| raw                | 7239    | 83392   | 4845    | False   | published_CDKN1A | 0    | True    | 37.55   |
| raw                | 7259    | 127386  | 4845    | False   | published_CDKN1A | 0    | True    | 40.75   |
| raw_sub20000       | 7239    | 20000   | 4845    | False   | published_CDKN1A | 0    | True    | 10.23   |
| raw_sub20000_seed1 | 7239    | 20000   | 4845    | False   | published_CDKN1A | 0    | True    | 10.27   |

### I. Seed-stability panel (one fixed 20,000-cell subsample of 7239)

| pair                   | n_cells | pearson_r | spearman_rho | anchor_sign_flipped | top5_jaccard | top5_n_changed |
|------------------------|---------|-----------|--------------|---------------------|--------------|----------------|
| dca_seed0 vs dca_seed1 | 20000   | 0.5703    | 0.6113       | False               | 0            | 2000           |
| dca_seed0 vs dca_seed2 | 20000   | 0.9824    | 0.984        | False               | 0.6653       | 402            |
| dca_seed0 vs raw_seed0 | 20000   | 0.7244    | 0.7062       | False               | 0.3201       | 1030           |
| dca_seed0 vs raw_seed1 | 20000   | 0.6887    | 0.673        | False               | 0.261        | 1172           |
| dca_seed1 vs dca_seed2 | 20000   | 0.5732    | 0.6121       | False               | 0            | 2000           |
| dca_seed1 vs raw_seed0 | 20000   | 0.5085    | 0.4805       | False               | 0.0309       | 1880           |
| dca_seed1 vs raw_seed1 | 20000   | 0.5066    | 0.4744       | False               | 0.0515       | 1804           |
| dca_seed2 vs raw_seed0 | 20000   | 0.7241    | 0.7003       | False               | 0.3236       | 1022           |
| dca_seed2 vs raw_seed1 | 20000   | 0.6883    | 0.6665       | False               | 0.2563       | 1184           |
| raw_seed0 vs raw_seed1 | 20000   | 0.9955    | 0.9954       | False               | 0.7606       | 272            |

(`results/phase8_d2/d2_stability.csv`, produced by `code/analyse_d2_stability.py`.)


---

## 8. Files and code

**New code** (all new files; nothing existing was edited except as noted):

| file | what it does |
|---|---|
| `code/setup_dca_env.sh` | builds the isolated CPython 3.8 + TF 2.4.4 venv that runs DCA. Never run against the main environment. |
| `code/dca_denoise_worker.py` | runs inside that venv; counts in as `.npz`, denoised matrix out as `.npy`. |
| `code/_shims_dca_bridge/dca/api.py` | a `dca.api.dca` that satisfies DeepScence's import and executes the real DCA out of process. **Separate directory from `code/_shims/dca`, which still raises and is unchanged.** |
| `code/run_deepscence_dca.py` | DeepScence at the published default, `denoise=True`, through the bridge. |
| `code/run_deepscence_denoise_probe.py` | DeepScence at `denoise=False` with `raw` / `mor` / `lib` / `ds10` pre-normalisation. |
| `code/analyse_d2_denoise.py` | the comparison: correlations, sender-status change at the operative thresholds, depth dependence, cell-type composition, normalisation strength. |
| `code/analyse_d2_stability.py` | the seed-stability panel (table I). |
| `code/run_d2_stability.sh` | launches that panel; small enough to run under contention. |
| `code/report_d2_tables.py` | renders the tables in §7 from the CSVs. |
| `code/run_d2_stream.sh` | the serialised runner (and the record of why it is serialised). |
| `code/run_d2_queue.sh`, `code/run_d2_queue_seq.sh`, `code/run_d2_probe_stream.sh`, `code/run_d2_dca_stream.sh`, `code/run_d2_dca.sh` | superseded launchers, kept because their headers record the memory findings that forced each retreat. |

**Outputs**

| file | content |
|---|---|
| `data/processed/deepscence_dca_<section>.csv` | scores at the published default `denoise=True` |
| `data/processed/deepscence_{mor,lib,raw}_<section>.csv` | `denoise=False` under each pre-normalisation; each carries `mapped_counts` and the per-cell `size_factor` used |
| `results/phase8_d2/d2_agreement.csv` | per section × configuration: Pearson/Spearman vs the committed scores, sign-flip flag, mean/95th-percentile \|Δz\|, and Jaccard / n-changed / %-of-calls-changed at global top-5%, global top-1%, within-type top-5%, and cell-type × depth-decile matched top-5% |
| `results/phase8_d2/d2_depth.csv` | ρ(score, transcript counts) committed vs alternative and its change; ρ(Δz, depth) signed and absolute |
| `results/phase8_d2/d2_deciles.csv` | per within-section depth decile: mean Δz, top-5% call rate committed vs alternative, % of cells flipping |
| `results/phase8_d2/d2_normalisation_strength.csv` | how much log-depth variance each normalisation actually removed |
| `results/phase8_d2/d2_celltype_composition.csv` | cell-type enrichment of the top-5% calls, committed vs alternative |
| `results/phase8_d2/d2_stability.csv` | the seed panel: every pairwise comparison among `raw` seed 0/1 and `dca` seed 0/1/2 on the fixed 20k subsample |
| `results/phase8_d2/d2_tables.md` | §7, regenerable |
| `results/phase8_d2/runmeta_<config>_<section>.json` | per run: cell/gene counts, denoise state, anchor, panel, the CDKN1A direction log (`node`, `reverse`), runtime, and for `mor` the size-factor diagnostics |
| `results/phase8_d2/dca_worker_meta_<section>.json` | the DCA side: TF/Keras/numpy/anndata/scanpy/Python versions, thread count, seed, runtime, and min/max/mean of the denoised matrix |
| `results/phase8_d2/dca_venv_pip_freeze.txt`, `dca_venv_python.txt` | the isolated environment, exactly |
| `results/phase8_d2/committed_deepscence_sha256.txt` | SHA-256 of all eleven committed DeepScence score files, so it is checkable that none was touched |

**Reproducing it**

```bash
DCA_ENV_ROOT=/some/scratch bash code/setup_dca_env.sh          # ~4 min
export DCA_VENV_PYTHON=/some/scratch/v38/bin/python
export DCA_BRIDGE_SCRATCH=/some/scratch/bridge
python code/run_deepscence_dca.py 7239_liver_sbr_Male_52-U1    # denoise=True
python code/run_deepscence_denoise_probe.py --config mor 7239_liver_sbr_Male_52-U1
python code/analyse_d2_denoise.py && python code/report_d2_tables.py
bash   code/run_d2_stability.sh && python code/analyse_d2_stability.py   # the seed panel
```

Run one section at a time, and check `/sys/fs/cgroup/memory.current` against
`/sys/fs/cgroup/memory.max` first — not `free`. DeepScence needs ~16 GB of cgroup headroom
for an 83k-cell section and more for the 200k+ ones.

## 9. What is not covered

* **Coverage is partial by configuration.** The committed `denoise=False` scores exist for
  11/11 M1 sections (task 8.3/D1); the configurations compared here do not, and §5 says
  why. Table A names exactly which sections each configuration covers. Nothing in §0 is
  claimed beyond the sections listed there.
* **`ds10`** (binomial downsampling of every cell to a common depth) is implemented in
  `run_deepscence_denoise_probe.py` and was **not run**. It was dropped when memory
  contention arrived: it destroys ~90% of the counts in deep cells, so a large score change
  under it would have been unattributable between "depth removed" and "information
  removed", and `lib` answers the depth question cleanly without that confound.
* **The seed panel is one subsample, one section, three seeds.** It is enough to show that
  a disjoint outcome occurs; it is not enough to estimate how often. The two full `dca`
  sections were each run once, at `random_state=0` — which seeds 0 and 2 agree is a
  representative solution, but seed 1 shows it is not the only one.
* **`lib` has two sections and they differ in magnitude a lot** (*r* 0.843 vs −0.017). Both
  agree on the direction and on the collapse in depth loading; neither is enough to say how
  often the near-orthogonal outcome on 7352 occurs.
* **`mor` has no section at the deep end above 1138 median counts**, and `mor 7250`/`mor
  7260` were lost to OOM and not retried.
* **H1 is untouched.** `data/raw_h1/` was not read. Mouse only, per the §15 freeze.
