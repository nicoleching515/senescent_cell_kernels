# Seven open flags in the frozen pre-registration — how each was resolved

**Date:** 2026-08-27. **Author:** documentation agent, resolving the seven open flags of
`reports/PREREG_DECISIONS_PHASE10.md` §3 against `reports/PREREG_PHASE8.md`.

**Files edited:** `reports/PREREG_PHASE8.md` (new **§0.2 FLAG-RESOLUTION BLOCK**, plus dated inline
markers at eleven sites) and this report. **Nothing else was written.** No file under `results/`,
`figures/`, `genesets/`, `data/`, `code/` or any other report was created or modified;
`data/raw_h1/` was not read; nothing was written under `results/phase9_h1/` or `results/phase10*`.

**Form:** identical to §0.0's `C-1 … C-11` and §0.1's `D-A … D-D` — **original wording preserved at
every site, dated `***[…]***` markers added beside it, no silent rewrite and no deletion.** Verified
mechanically after the edit: of the three existing lines the diff modifies in place (the §3.4 N7
row, §8's P-iii row, §13's P7 row), each old line is an exact prefix of its new line with only a
marker appended inside the cell; every other change is a pure insertion. Markdown table integrity
re-checked across the whole file: **0 blocks with mismatched column counts.**

**Gates, run at the end:** `python3 code/gate_genesets_guard.py` → **exit 0**, "no drift: 96 watched
gene-set/panel files match the manifest", MOUSE GATE PASS + HUMAN GATE PASS.
`python3 code/check_figures_guard.py` → **exit 0**, "OK: all **52** committed figures match".

**Commits (not pushed; no tag created or moved):** `4bdc1c9` (the §0.2 block), `d274ba3` (the inline
site markers), plus this report.

---

## Read this first — where I did not do what the ruling said

### F2: the ruling's substance is applied; the code edit it implies is **REJECTED**

**Ruling:** "`tierA_merged_p95` is in neither `N7_CALLS` nor its nine-call extension. Add it, do not
replace anything… Verify it actually runs after the edit. If adding it changes any existing output,
stop and report rather than proceeding."

**What I found:** the addition **already exists**, in the right place, in exactly the
append-don't-replace form the ruling asks for — and making it instead in `run_phase3_nulls.py`
**would change existing output**, which is the ruling's own stop condition.

- The H1 arm binding appends the merged calls **after** the nine:
  `code/h1_run_phase10.py:48` — `ALL12 = list(RN.ALL9_CALLS) + MERGED_CALLS`; the alias and the
  `sender_mask` wrapper are at `code/h1_phase10.py:42-62`, which deliberately does **not** edit the
  frozen `sasp_phase3.py`.
- **It runs.** All twelve calls, `tierAmg_p95` (= `tierA_merged_p95`) among them, are present in
  `results/phase10_h1/window.csv` and `results/phase10_h1/main_fits.csv`:
  `python3 -c "import pandas as pd; [print(f, sorted(pd.read_csv(f, low_memory=False)['call'].unique())) for f in ['results/phase10_h1/window.csv','results/phase10_h1/main_fits.csv']]"`
  → `['cdkn1a_pos','senepy_p95','senepy_p99','tierA_p90','tierA_p95','tierA_p99','tierAmg_p90','tierAmg_p95','tierAmg_p99','tierApm_p90','tierApm_p95','tierApm_p99']` in both.
- **Why editing the frozen list would have changed existing output.** `_expand`
  (`run_phase3_nulls.py:123-136`) seeds each job by the **index of the call in the list it is
  given**: `sd = int(base) + step_i * i + step_j * j` at `:131`. Since
  `ALL9_CALLS = N7_CALLS + TIERA_PM_CALLS` (`:77`), a seventh entry in `N7_CALLS` moves
  `tierApm_p90/p95/p99` from `j = 6,7,8` to `j = 7,8,9` and **changes the seed of every per-module
  null job**. Separately, the M1 cache holds no `flag_merged_pNN` array, so `--calls all` on the
  mouse arm would raise on the new name.
- Operational note: `run_phase3_nulls.py` and `code/h1_phase10.py` are imported by the Phase-10 job
  running during this task (`h1_run_phase10.py --stage perm_c1`) and by every fresh `loky` worker it
  spawns, so editing them mid-run is unsafe independently of the above.

**Recorded in** §0.2 F2 and at §3.4's N7 row. **No code was changed.**

### F3: the flag is broader than the evidence — three quantities are **not** marked provisional

The flag says "every H1 DeepScence magnitude in the Phase 9 audit is at `random_state = 0`". The
ruling ("mark every such magnitude PROVISIONAL") is applied, but "be precise" required excluding:

1. **P-iv's circularity, 29/33 = 0.8788.** It is a gene-list membership computation — DeepScence's
   shipped `coreGS_v2.csv` ∩ panel ∩ frozen Tier B — in which **no score is produced and no
   `random_state` enters**. `python3 -c "import json; print(json.load(open('results/phase7_jobA/gate_result_human.json'))['corescence'])"`
   → `n_on_panel 33`, `frozen_n_in_any_B 29`, `frozen_frac 0.8788`. **P-iv is unaffected by D-A.**
2. **Tier A × SenePy 0.874 (z −7.96)** and **SenePy Q5/Q1 28.5–228.1**: no DeepScence term.
3. **The seed-check statistics themselves** (*r* 0.3719, Jaccard 0.2107, 12,779 cells; ρ 0.3122 →
   0.2308): statements *about* the seed pair, not statements at one seed.

And one qualification recorded **with** the ruling rather than against it: **`stab_cdkn1a` = 1.00 is
itself computed on the seed-0 score**, so it is listed among the provisional magnitudes even though
D-C's verdict, which it supports, is retained as a direction. The argument for retaining it — that a
consensus cannot be less stable than its members, and that 1.00 in 7/7 is two folds in twenty above
the falsifier — is **asserted, not measured**, and §0.2 says so, so the consensus run can check it
rather than inherit it.

### F7: one number beyond the ruling, flagged rather than applied

The same audit table (`CS_PHASE9_H1_AUDIT.md` §10.5) carries a **second** seed-1 quantity with no
producer — the within-type **Q5/Q1 = 0.292** for P-v — by the identical check
(`code/h1_d2_analyse.py:92-116` emits only `rho_seed1`, into `results/phase9_h1/d2_depth.csv`). The
ruling names only 2.967. I recorded 0.292 as a flagged extension, **not** as part of the ruling, and
did not act on it further; the pre-registration uses it nowhere — 0.292 appears in `PREREG_PHASE8.md` only inside §0.2 F7's own flag.

---

## The seven, one by one

| # | Ruling | Accepted / rejected | Backing file |
|---|---|---|---|
| **F1** | R1/R2 scored on each arm's primary call | **ACCEPTED in full** | `results/phase9_h1/a3_prevalence_by_type.csv`; `results/phase3/m1_final_audit.txt` §3 (via §0.0 C-5) |
| **F2** | Add `tierA_merged_p95`, replace nothing | **Substance ACCEPTED; the implied edit to the frozen list REJECTED** — already added in the arm binding | `code/h1_run_phase10.py:48`, `code/h1_phase10.py:42-62`, `code/run_phase3_nulls.py:70-77,123-136`; `results/phase10_h1/{window,main_fits}.csv` |
| **F3** | Mark every H1 DeepScence magnitude PROVISIONAL | **ACCEPTED, with three exclusions named** (above) | `code/h1_deepscence.py:76`; `caller_technical_loading.csv`, `caller_within_type_depth_bias.csv`, `deepscence_anchor_h1.csv`, `caller_agreement_pooled.csv`, `caller_agreement_matched_significance.csv` (all `results/phase9_h1/`) |
| **F4** | P-iii MARGINAL, pending the consensus | **ACCEPTED in full** | `results/phase9_h1/caller_agreement_pooled.csv`, `results/phase9_h1/d2_stability.csv` |
| **F5** | All three asymmetries in one place, with costs | **ACCEPTED in full** | `results/phase9_h1/a8_ortho_sender_shift.csv`, `d2_stability.csv`, `a3_prevalence_by_type.csv`; `results/phase8_d2/d2_agreement.csv` |
| **F6** | Record the correction; retiring P7 costs A6 nothing | **ACCEPTED, verified in source** | `code/h1_a6_compartments.py:57-58, 85-86` |
| **F7** | Strike 2.967; correct Q5/Q1 to 228.1 | **ACCEPTED in full**, plus one flagged extension | `results/phase9_h1/caller_agreement_matched_significance.csv`, `d2_depth.csv`, `caller_within_type_depth_bias.csv`; `code/h1_d2_analyse.py:92-116` |

### F1 — R1 and R2 are scored on each arm's primary call

**Done:** §0.2 F1 states that the primary outcome, R1 and R2 are scored on the primary call of the
arm being scored — `tierA_p95` on M1, `tierA_merged_p95` on H1 — with the frozen-literal fine call
computed and reported on H1 as the declared sensitivity. Marked at **§5** and **§6** beneath the
existing D-B markers; §0.1 open item 1 is annotated as closed.

**Stated explicitly, as required, that this resolves an ambiguity rather than relaxing a
criterion**, with three checkable legs: (i) no threshold, bracket or direction moves (R1 is still
"IQR includes 0 **and** q75 < 0.50" after §0.0 C-4; R2 is still against **that arm's own** 80 %
bound); (ii) M1's outcome is untouched — 0.088, IQR [−0.017, 0.234]; 0.029 against 0.183 (§0.0 C-5);
(iii) **the moved call is not the easier call** — the merged call **adds** senders everywhere the
fine call fell short and removes none: T/NK 23 → 1,958 (SPLN14), **0 → 1,241** (SPLN43), 295 →
1,675 (SPLN44), 322 → 1,328 (SPLN21), and ≤ 1 cell of movement in SPLN07/24/30, so if anything it
makes an IQR that **excludes** 0 — non-replication — easier to observe. Command in §0.2 F1. **H1's
actual outcome is deliberately not stated.**

### F2 — see above. Recorded at §3.4's N7 row: on H1 the axis is nine frozen calls **plus** three merged calls, twelve in all; the M1 axis is untouched.

### F3 — the provisional set, enumerated

§0.2 F3 lists every affected number with its file and command: P-i's depth loadings
**+0.1822 … +0.3540**; P-v's within-type Q5/Q1 **0.244 … 1.169**; every `rho_partial_*` and rank in
`deepscence_anchor_h1.csv` (CDKN1A +0.1911 … +0.2540, LMNB1 +0.2029 … +0.2449, prolif
+0.0097 … +0.0451, consensus +0.0642 … +0.1695, ranks 1st–7th, `stab_prolif` 0.75, and
`stab_cdkn1a` itself); and every agreement ratio with a DeepScence term (**6.436** z 204.83, 2.069,
1.602, **1.102** z 5.67, 1.093, 0.998, 0.290; per-section **3.459–10.115** and 0.972–2.816). All
descend from `data/processed_h1/deepscence_h1_<section>.csv`, written by `code/h1_deepscence.py:76`
at `random_state=0`. Marked at **§8's preamble** and **§10.10**. **D-C's verdict stands.**

### F4 — P-iii is MARGINAL

Marked in **§8's P-iii row** and beneath the table. The evidence added beyond the flag: the
registered rule confirms at > 1.10 and falsifies at ≤ 1.05, so **(1.05, 1.10] is pre-registered as
neither** — H1's **1.102** sits 0.002 above the top of a band the pre-registration itself declared
indeterminate, on an estimator whose full-section seed-to-seed top-5 % Jaccard is **0.2107**. The
outcome is not restated in either direction; no threshold is moved and the falsifier is not
weakened.

### F5 — the three asymmetries, with costs

§0.2 F5 is a single table: **panel** (M1 remapped 4,845/5,097 vs H1 native 5,093 — cutting H1 to
the 2,425-gene intersection leaves 26 of 33 Tier A genes and moves the `tierA_p95` sender set at
Jaccard **0.5222–0.5747**, `a8_ortho_sender_shift.csv`); **DeepScence estimator** (D-A — any
cross-arm difference smaller than H1's between-seed spread is uninterpretable, and **there is no
mitigation on the M1 side**); **Tier A call level** (D-B — different label granularity, partially
mitigated because the fine call is also computed on H1, **not** mitigated on M1, where no merged
counterpart exists; §0.1 open item 3 stays open). Closing statement: a cross-arm DeepScence number
carries axes 1+2, a cross-arm sender-call number carries axes 1+3, **and no cross-arm number in the
project is affected by fewer than two of the three.** Marked at **§9 item 4**.

### F6 — P7's justification, corrected

Verified in the source before recording: `code/h1_a6_compartments.py:57-58` scores five compartments
including `D_spleen_marginal_zone`; `:85-86` builds the axis as
`pulp = score(D_spleen_red_pulp) − 0.5·(score(follicle) + score(tzone))`; and
`grep -n "marginal_zone" code/h1_a6_compartments.py` returns **line 58 only** — the label is scored
and enters the producer nowhere else. So P7's "would leave the A6 axis without its middle term" was
already inaccurate before H1 was read, and **retiring P7 costs the A6 axis nothing**, which
strengthens D-D. Promoted from §0.1's flag to a recorded correction inside **§13's P7 row**, with
P7's original wording preserved as D-D requires; no gene set, threshold or scored compartment
changes and `genesets/human/D_spleen_marginal_zone.txt` is not edited.

### F7 — one number struck, one confirmed

**(a) 2.967 struck as unsourced, no substitute.** `caller_agreement_matched_significance.csv` has no
seed column and gives SPLN21 Tier A × DeepScence = **1.174** (z 3.54) at seed 0 only;
`grep -rl "2\.967" results/phase9_h1/` → no match; the seed-check producer emits only `rho_seed1`.
**Consequence applied inside the pre-registration:** §0.1 D-A's own clause *"moves by a factor of
about 2.5"* is derived from the struck number and **is struck with it** (marker at the site). D-A's
reason survives on file-backed evidence alone (*r* 0.3719, Jaccard 0.2107, 12,779 cells; 0.3122 →
0.2308).

**(b) Q5/Q1 = 28.5–228.1**, per section 74.6 / 30.8 / 50.6 / **228.1** / **28.5** / 98.0 / 70.3
(`caller_within_type_depth_bias.csv`; SPLN24 = 4.105/0.018). **Propagation checked by command, not
assumed:** `grep -rn "224\.7" reports/ code/ figures/` returns exactly two hits —
`CS_PHASE9_H1_AUDIT.md:765` (the audit's own table, outside this task's edit scope and left for its
owner) and `PREREG_PHASE8.md:1403`, where it appears only as a citation beside the corrected value.
**The wrong value reached no claim, no figure and no code path.** Marked at **§14**; SenePy's status
(§14 item 4) remains open.

---

## What was deliberately not done

- **No code was modified.** F2 needed none, and the two files it would have touched are in the
  import path of a Phase-10 job that was running.
- **`CS_PHASE9_H1_AUDIT.md` was not edited**, so 2.967, 0.292 and 28.5–224.7 still stand in it. The
  pre-registration now records that they are struck / unsourced / superseded; correcting the audit
  itself belongs to its owner.
- **§1's tag hash is still uncorrected and §1 is still unedited** (flag closed elsewhere, not mine).
- **The consensus pooling rule is still unspecified**, and §0.2 repeats that it must be written down
  before the consensus is read.
- **No H1 outcome was stated or restated** — not R1's, not R2's, not P-iii's.
