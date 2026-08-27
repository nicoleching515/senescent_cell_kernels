# RECORD RECONCILIATION — the authoritative value for each contested quantity

**Written 2026-08-27**, against the specification in `reports/WRITING_PACK.md` (§0, §7 and the
33-item forbidden-claims checklist), `reports/PHASE8_ROADMAP_STATUS.md` and
`reports/CORRECTIONS.md`. Repo at tag `phase8-frozen`.

**Scope.** Reports and planning documents only. **Nothing under `results/`, `figures/`,
`code/`, `genesets/` or `data/` was modified**, no run was launched, `data/raw_h1/` was not
read, and no tag was created or moved. `python3 code/check_figures_guard.py` passes:
`OK: all 52 committed figures match`, exit 0.

**Discipline held throughout.** Every replacement number below comes from a file, with the
command. Historical audit records (`AUDIT_*`, `BIO_DELIVERABLE7_CLAIM_AUDIT.md`, the
`COMPLETED_TASKS.md` findings log) are **marked superseded, not rewritten** — their value is
that they capture what was believed at the time. The frozen pre-registration is corrected by a
**dated correction block plus inline markers**, never by a silent edit.

---

# 1. λ̂ — LEAD ITEM

## 1.1 The authoritative value

> **λ̂ = 14.7 µm** (14.7321), defined as the **median `lam_naive` over the 315 primary fits** —
> in-band six sections × sender call `tierA_p95` × `stratum == "all"`.
>
> **Source, emitted:** `results/phase3/summary_phase3.txt` §6, the `tierA_p95` row, column
> `medlam`, written by `code/summarize_phase3.py:221`.
>
> **Mandatory companion, always in the same sentence: IQR [7.0, 50.0] µm and 60 % of fits
> railed at a grid bound** (103 at the 7 µm floor, 86 at the 50 µm ceiling).

```bash
python3 - <<'PY'
import pandas as pd, numpy as np
d = pd.read_csv('results/phase3/main_fits.csv'); d['sec'] = d.section.str.split('_').str[0]
p = d[(d.call=='tierA_p95') &
      (d.sec.isin(['7259','7260','7001','7248','7352','7435'])) & (d.stratum=='all')]
print(len(p), p.lam_naive.median(), np.percentile(p.lam_naive,[25,75]), p.lam_railed.mean())
PY
# -> 315  14.7321271090776  [ 7. 50.]  0.6
```

## 1.2 Why 15.7 µm had to go, and where it came from

`λ̂ = 15.7 µm "pooled"` appeared in six documents (`SASP_Kernel_Master_Plan.md` §29 obj 7 and
§30 5.3; `CS_PHASE8_TORUS_VAR.md` ×4; `CORRECTIONS.md`; `SUBMISSION_PATCH_2026-08-29.md` ×2;
`PLAN_UPDATE_D12_D13.md`). **No file emits it**, and its only apparent provenance was the
torus report's own sentence *"2,215 µm … 141× the pooled λ̂"* — 2215 / 141 = 15.71 — i.e. it
was derived from the claim it was used to support.

**A likely true origin, found here.** The closest reproducible match anywhere in either results
tree is **15.716 µm**: the ***interior*** median (railed fits excluded) over the **pre-C6**
`tierA_p95` in-band fits **including the zonation-stratified rows** — 441 rows, so hepatocyte
fits are counted four times (`all`, `periportal`, `midzonal`, `pericentral`).

```bash
python3 -c "
import pandas as pd
d=pd.read_csv('results/phase3_pre_c6/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
p=d[(d.call=='tierA_p95')&(d.sec.isin(['7259','7260','7001','7248','7352','7435']))]
print(round(p[p.lam_railed==0].lam_naive.median(),3), len(p))"
# -> 15.716 441
```

Pre-C6, interior, and pseudo-replicated. **Not a defensible estimand under any reading, so 15.7
is withdrawn whichever way it arose.** No post-C6 definition I could construct returns it.

## 1.3 Why the pooled median, and not the other three

The four candidates all reproduce exactly:

| definition | value | emitted by a file? |
|---|---|---|
| **median over all 315 primary fits ("pooled")** | **14.7321** | **YES** — `summary_phase3.txt` §6, `medlam` |
| interior median (railed excluded) over 315 | 17.0627 | **YES** — `summary_phase3.txt` §1, "interior median 17.1 µm" |
| median over the 153 **reportable** fits | 16.0690 | no |
| interior median over the 153 reportable | 14.9943 | no |

**The reasoning, in order of weight.**

1. **It is emitted.** 16.07 and 14.99 are computed by nothing. This task exists because a
   number was invented; adopting an unemitted one would be the same failure one step later,
   and it would silently break the next time the pipeline is re-run.
2. **The pre-registration freezes the summariser, not a λ̂ estimand — so there is nothing to
   pre-register-match, and the tiebreak falls to the frozen code.** `PREREG_PHASE8.md` §3.1
   freezes the λ grid (7 µm floor, 50 µm ceiling, 40 log-spaced points) and the railing
   definition, and nothing else about λ̂. §5 then **explicitly declines** to make a length
   constant the estimand: *"λ̂ rails at a grid bound in a majority of M1 fits, so a fitted
   length constant is not the estimand; the amplitude and its survival under conditioning
   are."* Its M1 benchmark table quotes only the **railing rate** for λ̂ — no point value.
   **Consequence for the writer: any λ̂ in the paper is descriptive, never an estimand, and
   must be labelled with its population and its railing rate.**
3. **It shares a denominator with the caveat that must travel with it.** The 60 % railing rate
   is computed over the 315. Quoting a λ̂ over 153 beside a railing rate over 315 is exactly
   the base-mixing the writing pack forbids in §5.9.
4. **It is a valid median of a censored sample; the interior median is not.** 103 of 315 fits
   are censored at the floor and 86 at the ceiling, but with 32.7 % below and 27.3 % above,
   **the median order statistic falls in the interior** (the 55th of the 126 unrailed values).
   So 14.73 is a legitimate median despite the censoring. The *interior* median instead
   discards 60 % of the sample non-randomly and asymmetrically, which biases it upward — and
   it is unstable to the population choice (17.06 over 315 against 14.99 over 153) in a way
   the pooled median is not.
5. **The word "pooled" becomes true.** The phrase in circulation is "the pooled λ̂". 14.73 is
   literally the pooled median. 15.7 was never pooled anything.

**Against interest, stated plainly.** 14.7 is the *smallest* of the four candidates and the
furthest from the withdrawn 15.7. It makes the window-justification claim read "≈ 6.8 λ̂"
rather than the neater "≈ 6λ", and it is not the value nearest the number being replaced.

## 1.4 Every dependent claim, re-derived

| claim | was | **is** | ratio |
|---|---|---|---|
| N3-var median displacement 2,215 µm, in λ̂ units | 141× | **150×** | 2215 / 14.7321 = 150.35 |
| tile seams at a 1,200 µm tile side | ~76 λ̂ apart | **~81 λ̂ apart** | 1200 / 14.7321 = 81.45 |
| 100 µm fitting window ("≈ 6λ") | ≈ 6λ | **≈ 6.8 λ̂** | 100 / 14.7321 = 6.79 |
| N4-var median displacement 3,395 µm | — | **230× λ̂** | 3395 / 14.7321 = 230.45 |
| N3-occ median displacement 28 µm | — | **1.9× λ̂** | 28.3 / 14.7321 = 1.92 |

**Nothing load-bearing reverses. All three of the claims that matter move in the direction that
strengthens the argument they support**, so no conclusion is at risk:

- **The torus argument gets stronger.** N3-var displaces senders 150 λ̂, not 141 λ̂ — further
  decoupled from the correlation length, which is the whole point of the variant.
- **The tiling caveat gets stronger.** Seams ~81 λ̂ apart, not ~76, means an even smaller share
  of cells sits within a correlation length of a seam — which is precisely why the tiled null
  looks *more* conservative on real data than Mrkvička's prediction implies.
- **The window justification survives but must be written more carefully.** "≈ 6.8 λ̂" is still
  comfortably wide. **But the honest form adds the railing caveat in the same sentence: the
  100 µm window is 14× the 7 µm grid floor and only 2× the 50 µm grid ceiling, and 60 % of
  fits sit at one of those two bounds.** That sentence is now in `SASP_Kernel_Master_Plan.md`
  §30 5.3 and in `NOVELTY_ASSESSMENT.md` O5.

**One further λ̂ vintage found and corrected:** `PREREG_PHASE8.md` §13 deviation **P15** says
*"close to the median λ̂ of 12.8 µm"*. That is a **pre-C6 reportable-fit** median (pre-C6
`tierA_p95` reportable = 12.70; the pre-C6 `senepy_p95` reportable is 12.78). Handled as
correction **C-6/C-7** in that file's new §0.0; P15's verdict is unaffected.

---

# 2. The brackets are IQRs — and a genuine CI cannot be computed without a new run

## 2.1 The finding

`results/phase3/sf_summary.csv` carries `subset, null, n, q25, median, q75, frac_le_0,
frac_gt_05` — **no CI column**. The bracket is produced by `code/summarize_phase3.py:99`:

```python
a, b, c = np.quantile(v, [.25, .5, .75])
```

over the **per-fit SF point estimates**. There is **no bootstrap in it at all**.
`results/phase3/m1_final_audit.txt` names the two headline brackets `ctrl_amp_iqr` and
`SF_N2+N5+N6_iqr` — the file itself already calls them IQRs.

**So `0.088 [−0.017, 0.234]` and `0.029 [−0.007, 0.084]` are inter-quartile ranges across the
153 reportable fits.**

## 2.2 Can the pre-registered bootstrap give a real CI? No.

`PREREG_PHASE8.md` §3.6 freezes **400 replicates over 100 quantile blocks**. That machinery
emits **per-fit** CIs — `sf_n2n5n6_lo` / `sf_n2n5n6_hi` in `main_fits.csv` — and those are
enormous: median span **[−0.415, +0.381]**, two orders of magnitude wider than the IQR, because
they are intervals on a *single fit*, not on the median across fits.

```bash
python3 -c "
import pandas as pd, numpy as np
d=pd.read_csv('results/phase3/main_fits.csv'); d['sec']=d.section.str.split('_').str[0]
p=d[(d.call=='tierA_p95')&(d.sec.isin(['7259','7260','7001','7248','7352','7435']))&(d.stratum=='all')]
r=p[(p.beta_naive>0)&(p.beta_base_lo>0)]
print(np.median(r.sf_n2n5n6), np.percentile(r.sf_n2n5n6,[25,75]))
print('median per-fit bootstrap CI:', np.median(r.sf_n2n5n6_lo), np.median(r.sf_n2n5n6_hi))"
# -> 0.0885 [-0.0166 0.2338] ; median per-fit CI: -0.4153 0.3813
```

**No file emits an interval on the median across fits.** Producing one would require a new
block-bootstrap re-aggregated across fits — i.e. a re-run, which is out of scope and would
write to `results/`. **Per the rules, I stop and say so.** The resolution is therefore
**relabel, do not fabricate**.

Genuine CIs *do* exist in this project, for two quantities only, and both should be named as
CIs to keep the contrast visible: the **composition-matched SFs**
(`compmatch_reruns.csv`, `median_sf_matched_lo/hi`) and the **A7 section-clustered means**
(`a7_summary.csv`, `clustered_lo/hi/p`, a t-CI over sections).

## 2.3 A defect this exposed in the frozen pre-registration

`PREREG_PHASE8.md` §5 specifies the primary estimand as reported *"with its **paired-bootstrap
interquartile range**"*, and §6's replication criterion **R1** is written on a
*"**paired-bootstrap interval** that includes 0"*. **Both are misnomers** — the bracket is an
IQR from `np.quantile` with no bootstrap in it. Because R1 is a pre-registered replication
criterion, this is not cosmetic: it is corrected by dated note (item **C-4**), restated as
*"the IQR across the reportable fits includes 0 and its upper quartile is below 0.50"*, and
**neither the threshold nor M1's own outcome changes** (IQR [−0.017, 0.234] includes 0;
q75 = 0.234 < 0.50).

## 2.4 Relabelled in

`SASP_Kernel_Master_Plan.md` §29 (obj 1, obj 6) and §30 5.2/5.3 · `README.md` ·
`SUBMISSION_PATCH_2026-08-29.md` §8 · `PLAN_UPDATE_D12_D13.md` ·
`PHASE8_ROADMAP_STATUS.md` (PI-decisions table) · `BIO_DELIVERABLE6_DISCUSSION.md` ·
`PREREG_PHASE8.md` §0.0 (C-4). `CS_PHASE3.md` already labelled its column "IQR" and needed no
change.

---

# 3. Naive biological amplitude — one value, one estimator

> **Authoritative: +0.2767 response-SD**, the **section-clustered signed mean of β̂/sd(y) over
> the 1,155 biological-module fits under `design = base`**.
> **Source:** `results/phase3/a7_summary.csv`, row `BIOLOGICAL MODULES (reference)`,
> `design = base`, column `clustered_mean` — the frozen 09:06 file.
> **Conditioned counterpart: +0.0310** (`design = n6n5`).

```bash
python3 -c "
import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv')
print(d[d.response.str.startswith('BIO')][['design','clustered_mean','clustered_lo','clustered_hi','median_abs_amplitude']].round(4).to_string(index=False))"
```

**Why the clustered mean and not the median |β|/sd.** The clustered mean is A7's own primary
statistic: it is the only one of the two that carries an interval and a p-value
(`clustered_lo/hi/p`, a section-clustered t-CI), every A7 verdict in the paper — flat or not
flat, per control family and per design — is decided on it, and it is **signed**, so it is
directly comparable with the negative control gradients it is contrasted against. The median
|β|/sd is unsigned and interval-free.

| value | what it is | verdict |
|---|---|---|
| **0.2767** | section-clustered signed mean, `design = base`, **FROZEN** | **AUTHORITATIVE** |
| 0.3120 | median \|β\|/sd over the same 1,155 fits, frozen (`median_abs_amplitude`) | **legitimate companion — name the estimator explicitly whenever used** (conditioned: 0.0795) |
| 0.2914 | the clustered mean, **pre-C6** | **do not use** |
| 0.314 | the median \|β\|/sd, **pre-C6** | **do not use** |

**The power argument is unaffected: 0.362 SD exceeds all four.** Corrected in
`CS_PHASE8_MORAN.md` §4.2 (0.291/0.036 → 0.277/0.031), `NOVELTY_ASSESSMENT.md` banner (which
additionally attributed 0.291 to `moran_kernel_power.csv`, where it does not appear) and
`SUBMISSION_PATCH_2026-08-29.md` §9. `CS_PHASE8_CALLERS.md`'s 0.314 / 0.291 are covered by that
file's new whole-file pre-C6 banner.

---

# 4. Spearman ρ — aggregation-dependent, and the falsification is not

All four values in circulation reproduce **exactly**, from the same data:

| aggregation | ρ | p | emitted by |
|---|---|---|---|
| section-clustered **mean** per field, knn6 **raw**, 12 control+module fields | **+0.8951** | 8.37e-05 | `results/moran/moran_verdict.txt`; `code/summarize_moran.py:183` — **FROZEN; quote this** |
| section-clustered **mean** per field, knn6 **cell-type-centred**, 12 fields | **+0.9441** | 3.93e-06 | same, `:184` — **FROZEN** |
| **median** per field, knn6 raw, 12 fields | +0.9231 | 1.86e-05 | no file; re-derivable from `moran_vs_a7.csv` |
| **per-row**, no aggregation, 12 fields × 11 sections = 132 pairs | +0.7104 | 1.43e-21 | no file; same |

```bash
python3 -c "
import pandas as pd; from scipy.stats import spearmanr
P=pd.read_csv('results/moran/moran_pooled.csv'); s=P[P.kind.isin(['control','module'])]
print(spearmanr(s.I_raw_mean.abs(), s.a7_base_mean.abs()))
print(spearmanr(s.I_ctcentred_mean.abs(), s.a7_base_mean.abs()))"
```

**So +0.923 is not "not reproducible" — it is the median-per-field aggregation.** The writing
pack §0.4's diagnosis is corrected accordingly.

> **The falsification survives at every aggregation, and this must be said out loud.** ρ is
> **positive and significant under all four** (+0.71 to +0.94; the largest p of the four is
> 8.4 × 10⁻⁵).
> "Moran's I and the A7 kernel do not disagree" is therefore **not a fragile finding** — only
> the digit is aggregation-dependent. Any writer who states a bare ρ invites a reader to
> recompute it another way, get a different number, and conclude the finding turns on the
> choice. **State the aggregation in the same clause as the number, and state that the
> conclusion holds at all of them.**

**Not a fifth aggregation, and not a substitute:** within the 5 control responses × 11 sections
(55 pairs) ρ = **+0.155, p = 0.259**. That is a *different subset* of the same data — the
honest quantitative form of "the two statistics ask different questions", and the thing the
power argument predicts.

Applied in: `CS_PHASE8_MORAN.md` §0.4 and §4.1 (with the full table),
`NOVELTY_ASSESSMENT.md` banner, `SUBMISSION_PATCH_2026-08-29.md` §9 (which said "+0.923 raw"),
`WRITING_PACK.md` §0.4.

---

# 5. Every file changed, and why

## 5.1 Live working documents — corrected in place, with dated notes

| file | what changed |
|---|---|
| `SASP_Kernel_Master_Plan.md` | §29 obj 7: λ̂ 15.7 → **14.7**, seams ~76 → **~81 λ̂**. §30 5.2/5.3: brackets labelled **IQR across fits**; the window justification restated as **≈ 6.8 λ̂** with the railing caveat. §30 5.9: **no longer mixes bases** — now 4-pair throughout (1.131, p = 6.5e-9 → 1.212, p = 1.8e-106), with the consistent 3-pair alternative named. Tiled-torus inflation → **2.35×** with its basis (s ≥ 0.05). §29 obj 1/6 brackets labelled |
| `reports/CS_PHASE8_TORUS_VAR.md` | New head banner. **Audit R3's column fix finally applied**: §1 table "23 % / 8 % in the void" → **35.5 % / 19.9 % out of tissue**, and §10's framing blockquote likewise, with the 22.8 % retention column named separately. λ̂ dependents at L105/L291/L323/L408 → **150× / ~81 λ̂ / 150×**. §10's "1–63" admissible offsets → **1–66**. §4's irregular-window range now carries its basis (**0.048–0.118** overall, 0.080–0.118 for s ≥ 0.05) and **2.35×**; the file's internal 0.802 vs 0.801 rejection-rate contradiction resolved to **0.801**, its own §2 value |
| `reports/PHASE8_ROADMAP_STATUS.md` | **The self-contradiction is closed.** PI-decisions "Direction" row: pre-C6 0.326 / 0.027 / 0.082 / 0.203 → frozen **0.329 / 0.029 / 0.088 / 0.183**, brackets labelled IQR. The **8.4 gate table** marked pre-C6 with the frozen row given, and its "3-pair band" mislabel on the 6-section (4-pair) row corrected. **Tier A × SenePy 0.914 / below chance 11 of 11 struck as DEAD** → 0.972, p = 0.104, 4 of 11. The **8.5b A7 block** re-dated to the frozen 09:06 digits with the pre-C6 ones struck through, and **audit R5 recorded as moot**. FPR "9–13 %" → **9–16 %** with the pre-C6 0.127 → **0.145**, plus the forbidden "2–3× more fits" claim withdrawn. Figure-guard paragraph "27 of 45" → **52 / 52, closed**. N3-tile 0.974 → **0.971** |
| `reports/CS_PHASE8_MORAN.md` | ρ aggregation named at §0.4 and tabulated in full at §4.1; naive/conditioned biological amplitude 0.291 / 0.036 → **0.277 / 0.031** with the estimator named |
| `reports/NOVELTY_ASSESSMENT.md` | §2.1 point 3's falsified sentence **struck in the body** (it was banner-only before, so a linear reader hit a dead sentence); banner amplitude 0.291 → **0.277** with its mis-attribution to `moran_kernel_power.csv` noted; ρ aggregation table added; §2.2's "23 % in the void" → **35.5 % out of tissue**, "1–63 / 27 µm / λ̂ 12.8" → **1–66 / 28 µm / λ̂ 14.7**; O5's λ̂ corrected and given the railing caveat; §4's **O1 / O3 / O4 / O11 / O12 "NOT DONE" statuses superseded** |
| `reports/SUBMISSION_PATCH_2026-08-29.md` | λ̂ 15.7 → **14.7** (×2), seams **~81 λ̂**, ρ "+0.923 raw" → **+0.895** with the aggregation stated and the four-aggregation table added, amplitude 0.291 → **0.277**, §8 brackets labelled IQR, inflation → **2.35×** |
| `reports/PLAN_UPDATE_D12_D13.md` | brackets labelled IQR; seams ~76 → **~81 λ̂**; inflation → **2.35×** |
| `reports/CORRECTIONS.md` | §8.1's λ̂ paragraph corrected with a dated block giving the full derivation; §13.1's "worse in **6 of those 8** cells" → **5 of 8 strictly** (one exact 0.000 tie), with the command; **2.35× standardised** (§C.4 and §11); the two warnings that `SUBMISSION_PATCH` §9 still instructs the falsified Moran framing **closed**; §18.3's "still stale" list annotated with what is fixed and what remains |
| `reports/CS_PHASE8_M1_RERUN.md` | RS_count **0.040–0.060 → 0.033–0.060** in §6 and §14.3 (the minimum is the rectangle at s = 0.02); inflation → **2.35×** in three places; §9 item 6 marked **stale relative to its own §6**, and the 66–76 % row withdrawn outright |
| `reports/CS_PHASE5.md` | §4's **"76 % of the contact amplitude" withdrawn** in all three places it appears — the sourced ratio is 0.212 / 0.260 = **81.5 %**, and the claim is superseded by the composition-matched **65.9 % / 85.4 %** pair |
| `reports/BIO_DELIVERABLE6_DISCUSSION.md` | pre-C6 headline vector marked with its frozen replacement (bound **0.183**, tighter); the 76 % pointed at `CS_PHASE5.md` §4's withdrawal |
| `README.md` | headline block updated to the frozen vector (**153 of 315**, SF 0.088 IQR [−0.017, 0.234], 30 % ≤ 0, 91 % nuisance, controlled 0.029, bound **≤ 0.18**), bracket labelled IQR, pre-C6 values kept as a dated note; the caller paragraph's **1.13× / p = 4e-8 base-mixing** corrected to 1.131 / 6.5e-9 |
| `reports/WRITING_PACK.md` | the spec now carries each resolution beside its diagnosis: §0.1 λ̂ (with the 15.716 provenance finding), §0.4 ρ, §0.7 brackets, §5.5 amplitude, §5.3/§5.7 dependents, and a **status line for each of the 22 disagreements** |

## 5.2 The frozen pre-registration — dated correction block, no silent edits

`reports/PREREG_PHASE8.md` gains a new **§0.0 CORRECTION BLOCK** with seven items, and each
affected site keeps its original wording with an **inline dated marker** beside it:

- **C-1 — §10.1 item 1 mixes vintages inside one item.** Paragraph 1 quotes the pre-C6 05:19
  A7 file (−0.070 p = 0.023; N2 −0.061 p = 0.020; N5 +0.007 p = 0.41; conditioned biological
  +0.036) while paragraph 2 of the *same item* quotes the frozen 09:06 file. Frozen:
  **−0.0744 p = 0.0145 / −0.0642 p = 0.0124 (86 %, not 80 %, undiminished) / +0.0038 p = 0.715
  / +0.0310**.
- **C-2 — §3.6's FPR list** 0.091 / 0.103 / 0.109 / **0.127** / 0.164 is pre-C6 → **0.145**.
- **C-3 — §3.6 carries a live forbidden claim.** *"The reportable-fit filter therefore admits
  two to three times more fits than its nominal rate implies"* is checklist item 22 / audit R6.
  Withdrawn: the filter admits **3.0–13.3 %** naively and **4.8 %, identical across all five
  responses, on the full N6+N5 design**.
- **C-4 — "paired-bootstrap interquartile range" / "paired-bootstrap interval"** in §5 and R1.
  See §2.3 above.
- **C-5 — the §5 M1 benchmark table's frozen replacements**, including λ̂ railing 63 % → **60 %**
  and reportable 160 → **153**.
- **C-6 — deviation P15's pre-C6 digits**: 1–63 → **1–66**, 27/25 µm → **28/25 µm**,
  0.349/0.273 → **0.302/0.183**.
- **C-7 — the λ̂ withdrawal, recorded project-wide**, including P15's "median λ̂ of 12.8 µm".

**None of these changes a pre-registered decision, threshold, stop condition or the §18
outcome.** C-3 removes a forbidden claim; C-4 renames a bracket without moving a number; the
rest are pre-C6 → frozen substitutions the file had already flagged PROVISIONAL.

## 5.3 Historical records — marked superseded, never rewritten

| file | marker added |
|---|---|
| `reports/AUDIT_PHASE8_FACTCHECK.md` | **R5 marked MOOT**: `neg_probe_rate` under N6+N5 is now **+0.0097 [−0.0060, +0.0253], p = 0.199**, so "every control family is flat under +N6+N5" is **true** and the "4 of 5" hedge must not be carried forward. **R6's verdict stands, its digits superseded** (3.0–13.3 % naive; 4.8 % on the full design). Its own **8.3 %** N4 retention figure noted as predating the 05:58 regeneration → **8.0 %** |
| `reports/BIO_DELIVERABLE7_CLAIM_AUDIT.md` | **the existing marker did not suffice.** It was scoped to "the suggested rewording below", which left the *Additionally* paragraph at L294–300 reading as a live recommendation while repeating both forbidden figures (0.93–1.22× for four of six pairs — item 14; 1.51–2.85× — item 8). Scope **explicitly extended** to cover it, with the frozen replacements given. Text unedited |
| `reports/CS_PHASE7_C1.md` | **the existing marker did not suffice either.** It covered `null_destructiveness.csv` only. New banner adds (a) **every SF in the file is pre-C6 on the 160-fit population** — 0.974/1.000/0.349/0.721/0.962/0.273/0.716 → **0.971/0.999/0.302/0.695/0.924/0.183/0.707** — and (b) **§6 item 3, the last place the torus finding reads as novel statistics** (checklist item 19), whose recommendation that "tiling is the way out" is now against the settled position: Mrkvička predicts against tiling, our own calibration confirms it, and **N3-var is primary**. Text unedited |
| `reports/CS_PHASE8_CALLERS.md` | whole-file **pre-C6 banner** with the frozen substitution table, the instruction not to ship §3's drafted paragraph, and the note that §4.3's R6 sentence is unsupported. Two in-place fixes: the **arithmetic error 22 of 33 → 20 of 33** (0/11 + 11/11 + 9/11; **26 of 33** on the frozen base), which `CORRECTIONS.md` §5.2 recorded as "still uncorrected there"; and §4.1's **self-contradiction about `neg_probe_rate`** 15 lines apart, resolved to the frozen p = 0.199 |
| `reports/COMPLETED_TASKS.md` | the findings log left as written, with the six outcomes appended |

---

# 6. What I could not resolve, and what it would take

1. **A genuine confidence interval on the median SF or the median controlled amplitude.**
   Blocked by design, not by effort: the pre-registered 400 × 100-block bootstrap emits
   per-fit CIs only. **What it would take:** a new block-bootstrap that resamples blocks and
   re-computes the *median across fits* on each replicate, written to a new `sf_summary_ci.csv`
   — i.e. a run and a write to `results/`, both out of scope here, and both changes to a
   pre-registered procedure that would need a declared deviation. **Until then the brackets are
   IQRs and must be labelled as such.** The pre-registration's R1 has been restated accordingly
   without changing its threshold.

2. **The bibliography's author lines** (writing-pack disagreement 20). The *count* is settled:
   **43 entries**, `grep -c '^@' references.bib` — so "32 entries" is a stale vintage. The other
   half, "41 wrong given names across 19 of 32 entries", is an audit finding against the
   32-entry file and **cannot be restated for the 43-entry file without re-checking every author
   line against Crossref/PubMed.** That is a citation-verification task with network access, not
   a number-reconciliation task. The remedy the audit prescribes is procedural and stands: copy
   author names from Crossref/PubMed, never type them.

3. **Five rows of `CORRECTIONS.md` §18.3 that need re-derivation, not editing** — left open and
   flagged in place:
   - `CS_PHASE8_TORUS_VAR.md` §9's md5 certification of `sf_summary.csv` / `summary_phase3.txt`
     is against the **committed** hashes; the working tree is at the 09:06 files. True when
     written. Fixing it means re-certifying, which needs the pin policy decided.
   - `CS_PHASE8_TORUS_VAR.md` §6's "four excluded variance cases come back positive
     (+1.05 to +2.01)" — three are positive, the fourth is −0.717. Needs `var_variance_check.csv`
     re-read case by case.
   - `CS_PHASE8_CALLERS.md` §5.5's D3 Tier A column (1.248 / 1.237 / 0.935 →
     1.288 / 1.256 / 0.907). Covered by the new whole-file banner but not corrected in place.
   - `PREREG_PHASE8.md`'s **20 PROVISIONAL marks** and P29's "roughly doubles" (→ 32–67 %).
     These are frozen-document edits that belong in a second correction pass with the PI, not in
     a reconciliation sweep; §0.0 records the ones that carry numbers.
   - `README.md:293`'s 69 % circularity bar (→ mouse **79 % → 88 %**, human 76 % → 88 %).

4. **The N3-tile raw rejection rate 0.801 vs 0.802.** `CS_PHASE8_TORUS_VAR.md` disagreed with
   itself (§2 said 0.801, §4 said 0.802) and `CORRECTIONS.md` §18.3 and `COMPLETED_TASKS.md`
   both say 0.801, so §4 is resolved to **0.801** as an internal-consistency fix. **I could not
   independently re-derive either figure** — the obvious reconstruction from
   `perm_nulls_c1.csv` joined to the reportable fits does not reproduce it, so the definition
   used is not recoverable from the file names alone. **What it would take:** the producing
   snippet, or a line reference in `summarize_phase3_c1.py`. Flagged rather than asserted.

5. **λ̂ for the H1 arm.** Not applicable yet — the human arm has not run. But the definition
   settled here is arm-generic and should be applied identically, **including the railing rate**,
   which `PREREG_PHASE8.md` §6 R3(b) already pre-registers as a reported diagnostic with no
   threshold.

---

# 7. Standing rules this reconciliation adds

1. **λ̂ = 14.7 µm, pooled median over the 315 primary fits, and it never travels without
   "IQR [7.0, 50.0] µm, 60 % railed".** λ̂ is descriptive, never an estimand — the
   pre-registration says so.
2. **Every SF and amplitude bracket is an "IQR across fits" and is written as one.** The only
   genuine CIs in this project are the composition-matched SFs and the A7 clustered means.
3. **Naive biological amplitude = +0.2767, the section-clustered signed mean.** If you use
   0.3120 instead, write "median |β|/sd" beside it.
4. **Every Spearman ρ carries its aggregation in the same clause, and every statement of the
   Moran falsification says it holds at all four aggregations.**
5. **Tiled-torus inflation is 2.35×** — the exact 0.1175 / 0.05. "2.4×" rounded in our favour.
6. **Never merge `1 − frac_in_occupancy` (35.5 % / 19.9 %, out of tissue) with
   `1 − frac_retaining_a_neighbour` (22.8 % / 8.0 %, lost every real neighbour in the window)
   in one sentence.**
7. **Never mix the 3-pair and 4-pair caller bases inside one sentence.** The 4-pair values are
   the ones a file emits; the 3-pair values must be cited as an MH pooling of the significance
   CSVs.
