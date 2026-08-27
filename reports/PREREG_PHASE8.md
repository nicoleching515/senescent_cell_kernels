# Pre-registration — Phase 8 freeze (M1 mouse + H1 human)

**DRAFT for the PI to assemble and commit.** Written by the biology collaborator, 2026-08-27.
Nothing in this file has been committed and no tag exists yet.

**Companion file, not duplicated here:** `reports/PREREG_PHASE8_genesets.md` is the gene-set
section (Tiers A–E, both arms, the §11 gate result, the pins, and the gene-set deviation table
D1–D17). This file covers everything else §15 asks for and references that file rather than
restating it.

**Freeze status of the H1 data.** No H1 expression value, cell record or annotation has been read
by any bio agent. Panel membership only (`genesets/h1_candidate/`), which is the sanctioned §12.1
step-2 screen. `data/raw_h1/` is on disk ahead of the tag; §13 below records that as a declared
deviation.

---

## 0.0 CORRECTION BLOCK — dated 2026-08-27 (record reconciliation)

**This file is the frozen pre-registration. Nothing below is rewritten silently.** The
following items are corrected here, with the original wording left in place at each site and
an inline dated marker pointing back to this block. Full derivations and commands:
`reports/RECORD_RECONCILIATION.md`.

| # | Where | What is wrong | Correction, with source |
|---|---|---|---|
| **C-1** | **§10.1, paragraph 1** | **Mixes vintages inside one item.** Paragraph 1 quotes the **pre-C6 05:19** A7 file (−0.070 p = 0.023 naive; −0.061 p = 0.020 N2; +0.007 p = 0.41 N5; conditioned biological +0.036 SD) while paragraph 2 of the *same item* quotes the **frozen 09:06** file (−0.0604 / −0.0307 / −0.0225). A reader taking both paragraphs as one measurement is comparing two runs. | Frozen `results/phase3/a7_summary.csv`, `all_controls`: naive **−0.0744 [−0.1306, −0.0182], p = 0.0145**; N2 **−0.0642 [−0.1113, −0.0172], p = 0.0124** (86 %, not 80 %, undiminished); N5 **+0.0038 [−0.0186, +0.0261], p = 0.715**; conditioned biological amplitude **+0.0310** (clustered mean), not +0.036. Paragraph 2's per-family digits are already correct. |
| **C-2** | **§3.6** | Quotes the FPR list **0.091 / 0.103 / 0.109 / 0.127 / 0.164** — the 0.127 is the **pre-C6** `neg_control_probe` value. | Frozen: **0.091 / 0.103 / 0.109 / 0.145 / 0.164**. The range "9–16 % against a 5 % nominal" is unchanged; the four count-based responses span **9–15 %**. |
| **C-3** | **§3.6** | *"The reportable-fit filter therefore admits two to three times more fits than its nominal rate implies, in both arms."* **This is a forbidden claim** (audit R6, writing-pack checklist item 22) and it is live in the frozen text. | The 9–16 % is the **estimator's two-sided 95 %-CI exclusion rate**, not the filter's. The reportable filter (`beta_naive > 0 AND beta_base_lo > 0`, one-sided, naive design) admits **3.0–13.3 %** across the five control responses and **4.8 %, identical across all five — essentially nominal — on the full N6+N5 design** (`results/phase3/a7_control_probe_fits.csv`). Delete the sentence; keep the FPR statement. |
| **C-4** | **§5 (primary outcome) and §6 R1** | The primary estimand is specified as *"reported with its paired-bootstrap interquartile range"*, and R1's criterion is written on a *"paired-bootstrap interval"*. **There is no bootstrap in that bracket.** | The bracket is an **inter-quartile range across the reportable fits** of the per-fit SF point estimates: `code/summarize_phase3.py:99`, `np.quantile(v, [.25, .5, .75])` → `results/phase3/sf_summary.csv` columns `q25 / median / q75`. **There is no CI column in that file.** The pre-registered 400-replicate / 100-quantile-block bootstrap (§3.6) emits **per-fit** CIs (`sf_n2n5n6_lo/hi` in `main_fits.csv`, median span [−0.415, +0.381]) — it emits **no interval on the median across fits**. **R1 is therefore to be read as: the IQR across reportable fits includes 0 and its upper quartile is below 0.50.** Nothing about M1's outcome changes (IQR [−0.017, 0.234] includes 0; q75 = 0.234). Genuine CIs in this project exist only for the composition-matched SFs (`compmatch_reruns.csv`, `median_sf_matched_lo/hi`) and the A7 section-clustered means. |
| **C-5** | **§5, M1 benchmark table** | Marked **PROVISIONAL** and pre-C6, as intended — but it is being copied. | Frozen replacements, `results/phase3/m1_final_audit.txt` §3: naive amplitude **0.329** (was 0.326); SF N2+N5+N6 **0.088**, IQR **[−0.017, 0.234]** (was 0.082 [−0.099, 0.249]); SF N5 alone **0.115** (was 0.084); SF N1 **0.707** (was 0.716); SF N2 **0.952** (was 0.943); controlled amplitude **0.029** (was 0.027); 80 %-power bound **0.183**, SE 0.0654 (was 0.203, SE 0.073); λ̂ railed **60 % (189/315)** (was 63 %, 200/315); reportable fits **153** (was 160); corrected N3/N4 tile **0.971 / 0.924** (was 0.974 / 0.962) — and **N3-var 0.996 / N4-var 0.985 are now the primary corrected pair**, not the tile variants. |
| **C-6** | **§13, deviation P15** | *"close to the median λ̂ of 12.8 µm"*, and *"1–63 of 38,080–108,375"*, *"27 µm (N3) and 25 µm (N4)"*, *"0.349 / 0.273"*. All pre-C6. | Frozen `results/phase3/null_destructiveness.csv`: **1–66** admissible offsets of 38,080–108,375; median displacement **28 µm (N3-occ)** and **25 µm (N4-occ)**; SFs **0.302 / 0.183**. And the pooled λ̂ to compare against is **14.7 µm**, not 12.8 — see C-7. P15's verdict (the 5 % variants are degenerate, not a corrected null) is unaffected. |
| **C-7** | **project-wide, recorded here** | **λ̂ = 15.7 µm "pooled" is unsourced and is withdrawn.** It appears in six documents, is emitted by no file, and was back-derived from the claim it supported. This file never quotes 15.7, but its P15 quotes a different λ̂ vintage (12.8 µm). | **Authoritative: λ̂ = 14.7 µm**, the pooled median of `lam_naive` over the **315** primary fits (in-band × `tierA_p95` × `stratum == "all"`), printed by `code/summarize_phase3.py:221` into `results/phase3/summary_phase3.txt` §6, `tierA_p95` row, column `medlam`; re-derives as 14.7321 µm from `main_fits.csv`. **IQR [7.0, 50.0] µm, 60 % railed** — that caveat travels with it. This file's §5 position is unchanged and is the reason the value is only ever descriptive: **λ̂ is deliberately not the pre-registered estimand.** |
| **C-8** | **§13, deviation P2** | **Pre-C6 digits, and one stale *instruction*.** P2 quotes the **05:19** A7 file throughout (naive `all_controls` −0.070 [−0.128, −0.012] p 0.023; N2 −0.061 [−0.111, −0.012] p 0.020; N5 +0.007 [−0.011, +0.025] p 0.41; the 40 probes −0.018 p 0.18; biological +0.314 naive / +0.077 conditioned) and closes *"Magnitudes are PROVISIONAL — A7 was run at 05:19 on the pre-C6 sender calls … **A7 must be re-run after 8.7**."* **A7 *was* re-run.** `stat -c '%y' results/phase3/a7_summary.csv` → `2026-08-27 09:06:16`, against `2026-08-27 05:19:12` for `results/phase3_pre_c6/a7_summary.csv`. The instruction is discharged, not outstanding, and the magnitudes are no longer provisional. | Frozen `results/phase3/a7_summary.csv`, `all_controls`: naive **−0.0744 [−0.1306, −0.0182], p = 0.0145**; N2 **−0.0642 [−0.1113, −0.0172], p = 0.0124**; N5 **+0.0038 [−0.0186, +0.0261], p = 0.715**. `neg_control_probe` naive **−0.0225 [−0.0527, +0.0078], p = 0.129** — **still flat, so the pre-registered primary A7 response still passes** and P2's finding still rests on the codewords and genomic controls. Biological modules `median_abs_amplitude` **0.3120** naive / **0.0795** conditioned (clustered means +0.2767 / +0.0310). **All three qualitative findings in P2 are unchanged.** Command: `python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); print(d[d.design.isin(['base','n2','n5'])].to_string())"` |
| **C-9** | **§13, deviation P3** | **P3 instructs the writer to publish a pre-C6 range — this is the one correction in this block that changes an instruction rather than a digit.** Its body lists **0.091 / 0.103 / 0.109 / 0.127 / 0.164**; the **0.127** is the pre-C6 `neg_control_probe` value (the same stale digit C-2 corrects in §3.6, left uncorrected here), from which P3 derives *"the clean-null subset is **9.1–12.7 %**"* and then instructs *"**Quote the range as 9–13 %** with 16 % as the `neg_probe_rate` outlier."* **The instruction is wrong**, and it contradicts P3's own header, which already says 9–16 %. | Frozen `results/phase3/a7_summary.csv`, `design == "n6n5"`, `frac_CI_excludes_zero`: **0.091** (`neg_control_codeword`) / **0.103** (`all_controls`) / **0.109** (`genomic_control`) / **0.145** (`neg_control_probe`) / **0.164** (`neg_probe_rate`). The clean-null (count-based) subset is therefore **9.1–14.5 %**. **The instruction is corrected to read: quote the range as 9–15 % on the four count-based responses, with 16 % as the `neg_probe_rate` outlier.** The header ("9–16 % against a 5 % nominal") is correct as written and does not change. Note also that the pre-C6 list is wrong in more than one place: pre-C6 `all_controls` is 0.091 and `neg_control_codeword` is 0.109 — **the reverse of frozen** — so quoting that list as an ordered per-response list was already wrong before the vintage issue. Command: `python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); print(d[d.design=='n6n5'][['response','frac_CI_excludes_zero']].to_string(index=False))"` |
| **C-10** | **§13, deviation P24** | *"N2 matched decoys leave the technical gradient **~80 % intact** (−0.061 of −0.070 SD)"* — the pre-C6 ratio. C-1 corrects the same quantity to 86 % at §10.1 and leaves P24 at 80 %. | Frozen `results/phase3/a7_summary.csv`, `all_controls`: **−0.0642 / −0.0744 = 86.3 %**, with N5 still removing the gradient entirely (**+0.0038, p = 0.715**). P24's convergence argument is **strengthened, not weakened** — the matched-decoy design leaves more of the technical gradient than the pre-C6 file said, which is the direction P24 argues. |
| **C-11** | **§3.11 (the pin); §3.5 and §5 (the 160)** | **The pin is not a hash, and the file it pins was overwritten.** §3.11 pins `results/phase3/summary_phase3.txt` by md5 `ecf86b9ca5460f31290e2f4c9e822ea2` — that string is **31 hex characters**, which no MD5 can be, and it is not this file's digest under any vintage. The file itself was rewritten by task 8.7 at **2026-08-27 09:06**, so §3.5's *"the same 315/160 appears in the pinned `results/phase3/summary_phase3.txt`"* and §5's *"from the pinned `results/phase3/summary_phase3.txt`, primary call, 160 reportable fits"* no longer resolve to anything on disk. | The file on disk is **md5 `dc92ddc6605eef52f6359aeab4e16fd7`** (32 hex), mtime **2026-08-27 09:06**, and its line 19 reads `fits 315;  naive beta > 0 in 216;  positive AND block-bootstrap CI excludes 0: 153 (49 %)`. **The pinned digest is replaced by that one, and the reportable count everywhere in this file is 315 / 153, not 315 / 160** — consistent with C-5, which already carries 153. §3.11's three admissibility lists are unaffected and still re-verify exactly against `sasp_phase3.py:67-72`. Command: `md5sum results/phase3/summary_phase3.txt; sed -n 19p results/phase3/summary_phase3.txt` |

**None of these corrections changes a pre-registered decision, a threshold, a stop condition,
or the §18 outcome.** C-3 removes a forbidden claim; C-4 renames a bracket without moving a
number; **C-9 corrects a drafting *instruction* — P3 told the writer to publish 9–13 %, the
pre-C6 range — and C-8 discharges a stale instruction ("A7 must be re-run after 8.7": it was);**
the rest are pre-C6 → frozen digit substitutions the file already flagged as
PROVISIONAL. **No original wording has been deleted anywhere in this file; every site carries
its original text with a dated inline marker beside it.**

---

## 0.1 PI DECISION BLOCK — dated 2026-08-27 (four post-freeze decisions, D-A … D-D)

**This file is the frozen pre-registration and the freeze is real: `phase8-frozen` is committed and
pushed. Nothing below is rewritten silently.** The four decisions recorded here are taken by the PI
**after** the freeze and **after** H1 expression data was first read (Phase 9,
`reports/CS_PHASE9_H1_AUDIT.md`). They are recorded in exactly the form §0.0 uses: **the original
wording stays in place at every affected site**, and each site carries a dated inline marker
pointing back to this block. **Two of them replace a frozen estimator on one arm and one retires a
pre-registered question, so all four are declared deviations, not clarifications.** Every number
below was re-read from the file named beside it, by the command given; none is quoted from a
document. What these decisions contradict and do **not** resolve is listed at the end of this
block, and in `reports/PREREG_DECISIONS_PHASE10.md`.

*(Namespace note: **`D-A` … `D-D` is a new series.** It is not the **PI decision** D1–D16 series of
`PHASE8_ROADMAP_STATUS.md`, not the gene-set deviation series D1–D17 of
`PREREG_PHASE8_genesets.md` §12, not the Phase-8 `P*` rows of §13, and not the Phase-9 `H*` rows of
`CS_PHASE9_H1_AUDIT.md` §11.)*

| # | Decision | Kind | Sites marked |
|---|---|---|---|
| **D-A** | **DeepScence on H1 becomes a five-seed consensus score, reported with its between-seed dispersion.** M1 keeps its frozen single-seed score, so **the two arms use different DeepScence estimators** and that asymmetry is itself reportable | **New estimator, not in the frozen list.** Declared deviation | §3.9, §7 row F, §8, §9 item 4, §10.10 |
| **D-B** | **Sender labelling: both the fine-label and the merged-label Tier A call are computed; the MERGED call is PRIMARY on H1.** The frozen fine-label call is retained as the **frozen-literal sensitivity** | **Declared post-freeze choice** between two already-computed calls; no threshold tuned | §3.7, §5, §6 (R1/R2), §13 head-note |
| **D-C** | **P-ii is FALSIFIED and it leads the §8 reporting.** DeepScence's published `CDKN1A` anchor is stable in all seven H1 sections, so **M1's polarity flip is partly an artefact of our own ortholog remapping**, not a defect of the published tool | **A pre-registered prediction falsified and reported as such** | §8, §13 P12/P13 |
| **D-D** | **P7 (marginal-zone B cells) is retired as UNANSWERABLE** — never testable on this panel, as distinct from tested and null | **Question withdrawn**, original P7 text preserved | §10.4, §13 P7 |

---

### D-A — DeepScence on H1 moves to a five-seed consensus score

**The decision.** On the **H1 arm only**, the DeepScence score is no longer the single
`random_state = 0` run frozen in §3.9. The H1 primary is a **consensus over five seeds**, and the
**between-seed dispersion is reported beside every consensus number** — the spread is the finding,
not a nuisance to be averaged away. **This is a new estimator that is not in the frozen list**, and
it is recorded here as a deviation with its reason and its evidence.

**The seed set is already frozen.** The five values are the composition-matched seeds of §3.8 —
**20260901, 20260902, 20260903, 20260904, 20260905** — chosen there to be outside the range
`run_phase3_nulls._expand` can reach. **No new seed value is introduced by this decision.**
Command: `grep -n "^COMPMATCH_SEEDS" code/run_phase8_compmatch.py` → line 165,
`COMPMATCH_SEEDS = (20260901, 20260902, 20260903, 20260904, 20260905)`.

**The reason: the frozen single-seed configuration is not seed-reproducible on H1.** Full section,
`denoise=False`, SPLN21, all 196,142 cells, nothing changed but `random_state` 0 → 1:

| comparison | n | Pearson *r* | Spearman | top-5 % Jaccard | cells changing status | file |
|---|---|---|---|---|---|---|
| **H1 SPLN21, seed 0 vs 1, FULL section** | 196,142 | **0.3719** | 0.4157 | **0.2107** | 12,779 | `results/phase9_h1/d2_stability.csv` |
| H1 SPLN21, same, 20,000-cell panel | 20,000 | 0.3829 | 0.3282 | 0.2006 | 1,327 | " |
| **M1 seed-to-seed floor (M1's own `denoise=False` 0 vs 1)** | 20,000 | **0.99553** | 0.99543 | **0.7606** | 272 | `results/phase8_d2/d2_stability.csv` |
| M1 determinism control, **same** seed, full sections 7239 / 7259 | 75,384 / 114,721 | **0.99999913 / 0.99999995** | — | — | 24 / 2 | `results/phase8_d2/d2_agreement.csv` |

Commands:
`python3 -c "import pandas as pd; print(pd.read_csv('results/phase9_h1/d2_stability.csv').to_string(index=False))"`;
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase8_d2/d2_stability.csv'); print(d[d.pair=='raw_seed0 vs raw_seed1'].to_string(index=False))"`;
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase8_d2/d2_agreement.csv'); print(d[d.config=='raw'][['section','n_cells','pearson_r','global_top5_n_changed']].to_string(index=False))"`.

**r = 0.372 and Jaccard 0.211 against a floor of 0.9955 / 0.761.** It is not a small-sample
artefact: the 20,000-cell panel and the full 196,142-cell section give the same answer.
**Directions survive the seed; magnitudes do not** — on SPLN21 the depth loading moves only
+0.3122 → +0.2308 (`results/phase9_h1/d2_depth.csv`, `rho_seed1`), while the Tier A × DeepScence
matched agreement ratio moves by a factor of about 2.5 (see the open items below — that second
figure is quoted in `CS_PHASE9_H1_AUDIT.md` §10.5 and **has no producer file**, so it is *not*
relied on here).

***[Struck 2026-08-27 — see §0.2 item **F7(a)**. **The clause "moves by a factor of about 2.5" is STRUCK**, together with the 2.967 it is derived from: no file under `results/` emits a seed-1 Tier A × DeepScence matched ratio, and the producer of the full-section seed check (`code/h1_d2_analyse.py:92-116`) emits only `rho_seed1` into `results/phase9_h1/d2_depth.csv`. **No value is substituted** — a number that is needed needs a producer. **D-A's reason is unaffected**, because the evidence it rests on is file-backed without it: *r* = 0.3719, Spearman 0.4157, top-5 % Jaccard 0.2107 and 12,779 cells changing status (`results/phase9_h1/d2_stability.csv`), and the depth loading 0.3122 → 0.2308 (`results/phase9_h1/d2_depth.csv`).]***

**M1 keeps its frozen single-seed score, and the asymmetry is reportable.** The M1 harness
reproduces its committed scores at *r* = **0.99999913** and **0.99999995** at the same seed, and
its seed-to-seed floor is *r* = 0.99553 / Jaccard 0.7606 — so on M1 the single-seed score is not
the source of instability and there is no reason to change it. **The consequence is stated rather
than hidden: after this decision the two arms do not use the same DeepScence estimator** (M1:
single run, `random_state = 0`; H1: five-seed consensus). **That asymmetry is itself a reportable
finding** — it exists because the published tool is reproducible on the ortholog-remapped mouse
panel and is not reproducible on its own native human panel — and it must be stated wherever a
cross-arm DeepScence quantity is reported, alongside §9 item 4's requirement that every cross-arm
number appear on both panels.

**§3.9's reporting standard is extended from four attributes to five.** Every DeepScence number
now carries **coverage, denoise state, anchor, panel (native or ortholog-mapped, with the mapping
rate), and the seed configuration with its dispersion** — `random_state = 0`, single run, on M1;
five-seed consensus with its between-seed spread on H1.

**What this decision deliberately does not fix, and must be fixed before Phase 10 computes it:**
the **pooling rule** for the consensus (per-cell mean of z-scored per-seed scores, rank-mean, or
median) and the **dispersion statistic** reported beside it are not specified here and are not
invented here. Whatever is chosen is a producer-level choice that must be written down before the
consensus is read, not after.

---

### D-B — Fine and merged labelling: both computed, **MERGED named primary for H1**

**The decision.** The Tier A sender call is computed at **both** label families on H1. The
**merged-label call `tierA_merged_p95` is PRIMARY**, because it matches the label family the
estimator actually stratifies on (`sasp_phase3.LABELS = "merged"`). The **frozen fine-label
`tierA_p95` call is retained and reported as the frozen-literal sensitivity.** **This is a declared
post-freeze choice**, not a clarification of the frozen text: §3.7 freezes the percentile rule
*within cell type* as `phase2_downstream.py` implements it, which is the fine family.

**The reason — the frozen combination leaves cells eligible but uncallable.** The call thresholds
within **fine** cell types while the estimator's receiver stratification and sender-eligibility
mask use **merged** labels, so a cluster that is `Unknown` at the fine level but resolved at the
merged level is an eligible cell that can never be called a sender. Quantified, per section:

| section | fine-`Unknown` but merged-assigned | QC-passed cells | **% eligible but uncallable** |
|---|---|---|---|
| SPLN07 | 0 | 227,360 | **0.00** |
| SPLN14 | 38,701 | 283,628 | **13.64** |
| SPLN21 | 20,155 | 196,142 | **10.28** |
| SPLN24 | 7,386 | 393,202 | **1.88** |
| SPLN30 | 0 | 291,577 | **0.00** |
| SPLN43 | 24,815 | 270,472 | **9.17** |
| SPLN44 | 27,588 | 299,897 | **9.20** |

i.e. **0–13.6 % of cells per section.** Command (denominator is `n_cells_qc` from the annotation
metadata, not a count taken from a document):

```
python3 -c "
import pandas as pd, json
BAD={'Low_quality','Unknown','unknown'}
for s in ['SPLN07','SPLN14','SPLN21','SPLN24','SPLN30','SPLN43','SPLN44']:
    d=pd.read_csv(f'data/processed_h1/celltypes_h1_{s}.csv')
    n=json.load(open(f'data/processed_h1/annotation_meta_h1_{s}.json'))['n_cells_qc']
    u=int(((d.cell_type.isin(BAD))&(~d.cell_type_merged.isin(BAD))).sum())
    print(s,u,n,round(100*u/n,2))"
```

**The cost in senders, at the primary threshold, on the merged `T/NK cells` stratum**
(`results/phase9_h1/a3_prevalence_by_type.csv`, `label_set == "merged"`; command:
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/a3_prevalence_by_type.csv'); m=d[(d.label_set=='merged')&(d.cell_type=='T/NK cells')&(d.call.isin(['tierA_p95','tierA_merged_p95']))]; print(m[['section','call','n_all','n_senders_all','prevalence_pct_all','passes_A3']].to_string(index=False))"`):

| section | cells | senders, **frozen fine-label call** | prevalence % | senders, **merged-label call** | prevalence % |
|---|---|---|---|---|---|
| SPLN14 | 39,158 | 23 | 0.059 | 1,958 | 5.000 |
| **SPLN43** | **24,815** | **0** | **0.000** | 1,241 | 5.001 |
| SPLN44 | 33,482 | 295 | 0.881 | 1,675 | 5.003 |

**Zero senders in 24,815 cells** in SPLN43, and all three failures are *band* failures below the
1 % floor — the signature of the mechanism, not of a thin population. Under the merged call T/NK
passes A3 in 7 of 7 sections and every other row of the A3 table is unchanged; **no threshold is
tuned — it is the identical percentile rule applied at the other label family.**

**What is reported.** Both calls, every time, exactly as P11 already requires for the two Tier A
variants: the merged call as primary and the frozen-literal fine call beside it, with the
uncallable fraction above quoted wherever the difference between them matters.

---

### D-C — **P-ii is FALSIFIED. This is the headline of the §8 experiment, and it is reported as one.**

**A pre-registered prediction was made, it was falsified, and it is reported at the front rather
than in a table.** §8 predicted that DeepScence's published `CDKN1A` anchor would be *weak,
unstable or inverted* — depth-partialled fold-split sign stability < 0.90 — in **≥ 1 of 7** H1
sections, with the falsifier written in advance as "all 7 sections have stability ≥ 0.90".

**Measured: stability = 1.000 in all seven sections — the anchor decides the polarity the same way
in 20 of 20 random folds in every section.** `results/phase9_h1/deepscence_anchor_h1.csv`,
produced by `code/h1_deepscence_anchor.py`, whose `fold_stability(v, k=20, seed=0)` at line 103
defines the 20 folds and which imports `partial_spearman` from the mouse producer
`code/deepscence_reanchor.py` so the definition cannot drift.

| quantity | value across the 7 sections | file |
|---|---|---|
| `CDKN1A` fold-split sign stability | **1.00 in all 7** | `results/phase9_h1/deepscence_anchor_h1.csv`, `stab_cdkn1a` |
| depth-partialled ρ(`CDKN1A`, score) | **+0.1911 … +0.2540**, positive in all 7 | " , `rho_partial_cdkn1a` |
| `CDKN1A` rank among the 33 on-panel CoreScence genes | **1st – 7th** in all 7 | " , `rank_cdkn1a_in_core` |
| the D3 alternative 8-gene proliferation anchor, for contrast | ρ +0.0097 … +0.0451, and **unstable (0.75) in SPLN21** | " , `rho_partial_prolif`, `stab_prolif` |

Command:
`python3 -c "import pandas as pd; print(pd.read_csv('results/phase9_h1/deepscence_anchor_h1.csv')[['section','rho_partial_cdkn1a','stab_cdkn1a','rho_partial_prolif','stab_prolif','rank_cdkn1a_in_core']].to_string(index=False))"`.

**The consequence, stated plainly because it runs against our own result: M1's polarity flip
(§13 P12, P13) is at least partly an artefact of our own ortholog remapping, not a defect in the
published tool.** §8 exists to separate exactly two explanations — "a property of our mouse
adaptation" against "a property of the published tool" — and on this axis it returned the first.
On the native human panel the published anchor is the *better* anchor, and the D3 alternative
anchor chosen on M1 is the one that carries almost no signal and is itself unstable in a section.

**What this falsification does NOT touch.** The narrowing is explicit, because a falsified
prediction is not a general acquittal. All three of the following are measured on the **native** H1
run and are **unaffected by ortholog remapping**:

1. **The 88 % CoreScence circularity.** **29 of 33** on-panel CoreScence genes sit in ≥ 1 frozen
   Tier B module = **0.8788**, measured on the native human panel.
   `results/phase7_jobA/gate_result_human.json`, keys `corescence.frozen_n_in_any_B` = 29,
   `corescence.n_on_panel` = 33, `corescence.frozen_frac` = 0.8788. Command:
   `python3 -c "import json; print(json.load(open('results/phase7_jobA/gate_result_human.json'))['corescence'])"`.
   This is §8's P-iv, **CONFIRMED**.
2. **The seed instability (D-A above).** *r* = 0.3719 / Jaccard 0.2107 across seeds at full section
   size on the native panel, `results/phase9_h1/d2_stability.csv`.
3. **DeepScence × `CDKN1A`⁺ agreement at 6.436 pooled, natively.** The circular pair is
   **6.436**, z = 204.83, above chance in **7 of 7** sections, per-section range **3.459 – 10.115**
   — against a pooled **1.255** on the ortholog-remapped mouse arm.
   `results/phase9_h1/caller_agreement_pooled.csv` (row `deepscence_score, cdkn1a_counts`,
   `circular = True`) and `caller_agreement_matched_significance.csv` for the per-section range.
   Commands:
   `python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_agreement_pooled.csv'); print(d[(d.A=='deepscence_score')&(d.B=='cdkn1a_counts')].to_string(index=False))"`;
   `python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_agreement_matched_significance.csv'); m=d[(d.A=='deepscence_score')&(d.B=='cdkn1a_counts')]; print(m[['section','ratio_stratified','z']].to_string(index=False))"`.
   **§10.7's prohibition is unchanged** — "1.51–2.85×" may still never be quoted; the H1 figure is
   the pooled 6.436 with its 3.46–10.12 range, and this pair stays excluded from every pooled
   cross-caller claim (`BIO_PHASE3.md` §4.4).

**And it does not vacate P12.** P12 excludes `LMNB1` as the D3 primary anchor because `LMNB1` is a
member of `B_downstream_arrest` and `B_secondary_senescence` — a **gene-set membership** fact about
our own modules, not a claim about anchor performance. `LMNB1` tracks `CDKN1A` closely on H1
(`rho_partial_lmnb1` +0.2029 … +0.2449) and remains a **secondary** for that reason alone.

---

### D-D — P7 (marginal-zone B cells) is retired as **UNANSWERABLE**

**The decision.** P7's exploratory marginal-zone claims are **withdrawn as unanswerable on this
panel**, by dated note. **The PI's grounds are that the question was never testable on this
dataset — which is a different statement from "tested and null", and is recorded as the different
statement it is.** No MZ result is being reported as negative; no MZ result exists to report.

**The evidence.** The label is **never realised in any of the seven H1 sections.** It is present in
the frozen marker set (`code/markers_human_spleen.py:27`, 8 markers, and in the merged B
compartment at `:52`) but no cluster in any section is annotated to it. Command:

```
python3 -c "
import pandas as pd
labs=set()
for s in ['SPLN07','SPLN14','SPLN21','SPLN24','SPLN30','SPLN43','SPLN44']:
    labs |= set(pd.read_csv(f'data/processed_h1/celltypes_h1_{s}.csv').cell_type.unique())
print('Marginal zone B cells' in labs, sorted(labs))"
```

→ `False`. It is one of six labels never realised in any section
(`Lymphatic endothelium`, **`Marginal zone B cells`**, `Mesothelial cells`, `Pericytes`,
`Sinusoidal endothelium`, `pDC` — `CS_PHASE9_H1_AUDIT.md` §1.5).

**P7's original text is preserved in §13 in full, and deliberately so.** It records that the weak
evidence tier was flagged **in advance**: CellMarker 2.0 has only **6** spleen rows for this type,
so its eight markers come from the all-tissue fallback at ≥ 1 PMID — the weakest tier the marker
build defines. That is verifiable, and it is the part of P7 worth keeping:
`grep -c "^Marginal zone B cells" genesets/human/markers_spleen_evidence.csv` → **8** marker rows,
each carrying the provenance string `ALL-TISSUE FALLBACK >=1 PMID -- WEAKEST EVIDENCE (6 spleen rows)`
(command: `grep "^Marginal zone B cells" genesets/human/markers_spleen_evidence.csv`), and
`results/phase7_jobA/build_markers_human_spleen.log` line 18 records the same tier at build time.
**A pre-registration that flagged its own weakest input in advance, and then found that input
produced nothing, is a record worth keeping intact.**

**What is retired and what is not.** Retired: every exploratory marginal-zone-*specific* claim, and
P7's own conditional ("conditional on the label surviving the post-freeze re-gate" — it did not
survive). **Not retired, and unchanged:** §10.4's prohibition on any MZ-specific confirmatory
claim (it was already forbidden, and remains so a fortiori); the frozen marker file itself, which
is under `genesets/.geneset_manifest.json` and **is not edited**; and the frozen gene set
`genesets/human/D_spleen_marginal_zone.txt` (8 genes), which is scored as one of the five A6
compartments — see the open items below for why P7's stated reason for keeping the label does not
survive contact with how A6 was actually built.

---

### Open items these four decisions create and do **not** resolve

Recorded here rather than left for a reader to notice. A known inconsistency is better than one
nobody has seen.

1. **D-B contradicts the literal text of §5, §6 R1 and §6 R2**, which name the primary call as
   `tierA_p95`. Under D-B the H1 primary is `tierA_merged_p95`. **R1's threshold, its direction and
   M1's own outcome are untouched**; what moves is which H1 fit population R1 is evaluated on. The
   PI has not stated whether R1 is evaluated on the merged call with the fine call as sensitivity
   (the reading consistent with D-B) or on both with neither named first. **Not resolved here.**
2. **D-B leaves the N7 sender axis unspecified.** `N7_CALLS` (`run_phase3_nulls.py:68-69`,
   six calls, nine with `TIERA_PM_CALLS`) does not contain `tierA_merged_p95`. Whether the merged
   call **replaces** `tierA_p95` in the N7 axis on H1 or is **added** to it is a frozen-list
   question this decision does not answer. **Not resolved here.**
3. **D-B creates a second cross-arm estimator asymmetry**, alongside D-A's. M1's published sender
   call is the fine-label one and the same fine/merged interaction exists on M1 (smaller —
   `BIO_PHASE3.md` §1.1). Whether M1 is re-called at the merged family for cross-arm comparability,
   or the asymmetry is declared and carried, is **not resolved here.**
4. **D-A and D-C interact, and the interaction is not resolved.** Every H1 DeepScence number in
   `CS_PHASE9_H1_AUDIT.md` §10 — including D-C's own 6.436 circularity, the P-ii anchor
   stabilities, P-i's depth loadings and P-iii's pooled 1.102 — is computed at `random_state = 0`.
   Once the H1 primary is the five-seed consensus, those quantities are **stated on a
   non-primary estimator** until they are recomputed. The directions are seed-robust; the
   magnitudes are not.
5. **P-iii's verdict is seed-fragile and sits on its own threshold.** §8 registers "pooled ratio
   > 1.10" and H1 returns **1.102** (`caller_agreement_pooled.csv`, `tierA_score` ×
   `abs_deepscence_score`, z = 5.67). A margin of 0.002 on an estimator whose full-section
   seed-to-seed Jaccard is 0.211 is not a margin. **P-iii must be re-evaluated on the consensus
   score before it is reported either way.**
6. **P7's stated reason for keeping the MZ label does not match how A6 was built.** P7 argues that
   dropping the label "would leave the A6 axis without its middle term". The H1 axis is
   `pulp = score(D_spleen_red_pulp) − mean(score(follicle), score(tzone))`
   (`code/h1_a6_compartments.py:85-86`): **`D_spleen_marginal_zone` is scored as one of the five
   compartments but does not enter the axis at all.** So retiring P7 costs the A6 axis nothing,
   and P7's justification for keeping the label was already inaccurate before H1 was read.
   Flagged, **not corrected** — it is P7's own wording and D-D preserves it.
7. **P-vi remains contradicted-but-not-falsified** and is untouched by all four decisions: its
   falsifier needs 5 of 7 sections and Phase 9 ran 1 (`CS_PHASE9_H1_AUDIT.md` §10.4). It must keep
   being described that way.
8. **One figure in the evidence chain has no producer file.** `CS_PHASE9_H1_AUDIT.md` §10.5 quotes
   the Tier A × DeepScence matched ratio moving 1.174 → **2.967** between seeds on SPLN21. The
   seed-0 value 1.174 is in `results/phase9_h1/caller_agreement_matched_significance.csv`; **the
   seed-1 value is in no file in `results/`**, so it is not relied on anywhere in this block. It
   should be regenerated into a file or dropped from the audit.
9. **Recorded observation, no decision taken: §1's tag hash no longer matches the tag.** §1 records
   the `phase8-frozen` commit as `926439629a07269a32c93f998da0f6e1cd20933c`, 2026-08-27 15:32 UTC.
   The tag on disk is an annotated tag dated 2026-08-27 19:55:58 UTC, subject
   "Phase 8 frozen (**RE-CUT**): pipeline and pre-registration for the human replication", and it
   resolves to **`d04691e2692a7be8d1ff676d2fb74ad9d1df049d`**; 9264396 is an ancestor of it.
   Commands: `git rev-list -n1 phase8-frozen`;
   `git for-each-ref refs/tags/phase8-frozen --format='%(objecttype) %(taggerdate) %(subject)'`.
   **§1 is deliberately left unedited** — the re-cut is outside these four decisions and needs the
   PI's own record of why the tag moved.

***[Status of these nine open items, 2026-08-27 — see §0.2. **Seven of the nine are now closed** and the closures are recorded in the new §0.2 block, not here: item 1 → **F1** (R1/R2 are scored on each arm's primary call), item 2 → **F2** (the merged call is ADDED to the H1 sender axis, never substituted, and it already runs), item 3 → **F5** (registered with D-A's asymmetry and the panel axis, with the cost of each), item 4 → **F3** (every H1 DeepScence magnitude marked PROVISIONAL, with the affected numbers enumerated and the unaffected ones named), item 5 → **F4** (P-iii is **MARGINAL**, not confirmed, pending the consensus), item 6 → **F6** (recorded as a correction: retiring P7 costs the A6 axis nothing), item 8 → **F7(a)** (2.967 struck as unsourced). **Two remain exactly as written above:** item 7 (P-vi stays contradicted-but-not-falsified) and item 9 (**§1's tag hash — still an observation, §1 still unedited**). The consensus pooling rule that D-A leaves unspecified is also still unspecified, and must still be written down before the consensus is read.]***

**None of these four decisions changes a pre-registered threshold, a stop condition, the §18
outcome table, or any M1 number.** D-A and D-B change *which estimator* an H1 quantity is computed
with, on one arm, with both alternatives computed and reported; D-C reports a registered prediction
as falsified and states precisely what the falsification does not reach; D-D withdraws a question
that was never answerable on this panel and keeps its original text. **No original wording has been
deleted anywhere in this file; every affected site carries its original text with a dated inline
marker beside it.**

---

## 0.2 FLAG-RESOLUTION BLOCK — dated 2026-08-27 (F1 … F7, closing seven of the nine open flags)

**Same rules as §0.0 and §0.1, and for the same reason: this file is the frozen pre-registration
and nothing in it is rewritten silently.** Every item below leaves the original wording in place at
each site and adds a dated inline marker beside it. Every number was re-read from the file named
beside it, by the command given; none is quoted from a document. `F1 … F7` is a **new series**,
disjoint from `C-1 … C-11` (§0.0), from `D-A … D-D` (§0.1), and from every `D*`/`P*`/`H*` series
already in play; it numbers the **flags raised against §0.1 itself** in
`reports/PREREG_DECISIONS_PHASE10.md` §3, and it exists so that a reader of the pre-registration
alone can see how each was closed.

**Two of the nine flags are closed elsewhere and are not restated here:** the consensus pooling
rule (§3 item 5 there; §0.1 D-A already records that it is deliberately unspecified and must be
written down before the consensus is read) and §1's stale tag hash (§3 item 9 there; §0.1 open item
9 records it as an observation and **§1 is still deliberately unedited**).

| # | Flag, as raised in `PREREG_DECISIONS_PHASE10.md` §3 | §0.1 open item | Resolution | Sites marked |
|---|---|---|---|---|
| **F1** | R1/R2 name `tierA_p95`; D-B makes `tierA_merged_p95` the H1 primary | 1 | **Ruling applied.** R1 and R2 are scored on **each arm's own primary call**; the frozen-literal fine-level call is the declared sensitivity. **An ambiguity is resolved; no criterion is relaxed** | §6 (below the existing D-B marker), §5 |
| **F2** | `tierA_merged_p95` is in neither `N7_CALLS` nor `ALL9_CALLS` | 2 | **Ruling applied as an ADDITION — but in the H1 arm binding, not in the frozen mouse list, and it was already there.** Editing `run_phase3_nulls.py` is **rejected**: it would change existing per-module null seeds and break the M1 arm | §3.4 (N7 row + marker below the table) |
| **F3** | Every H1 DeepScence magnitude is at `random_state = 0` | 4 | **Ruling applied**, with the affected numbers enumerated exactly. D-C's **verdict** stands as a direction; its **magnitudes**, including `stab_cdkn1a`, are PROVISIONAL | §8 preamble, §10.10 |
| **F4** | P-iii is confirmed by 0.002 on a seed-unstable estimator | 5 | **Ruling applied.** P-iii is **MARGINAL, pending re-score on the consensus**; its outcome is not restated either way | §8 (P-iii row + marker below the table) |
| **F5** | D-A and D-B add two cross-arm asymmetries on top of the panel axis | 3, and §9 item 4 | **Ruling applied.** All three are registered here in one place, each with what it costs the cross-arm comparison | §9 item 4 |
| **F6** | P7's own justification (the A6 "middle term") was already wrong | 6 | **Ruling applied, verified in the source.** Retiring P7 costs the A6 axis nothing — this **strengthens** D-D. P7's original wording is still preserved | §13 P7 row |
| **F7** | Two audit numbers do not reproduce from any file | 8 | **Ruling applied.** (a) the seed-1 ratio **2.967 is struck as unsourced**, and D-A's own "about 2.5" sentence is struck with it; (b) SenePy Q5/Q1 is **28.5–228.1** from the committed file, and the wrong value did not propagate | §0.1 D-A (below the seed-instability paragraph), §14 |

---

### F1 — R1 and R2 are scored on each arm's PRIMARY call

**The resolution.** **R1 and R2 are scored on the primary sender call of the arm being scored** —
`tierA_p95` on M1, and on H1 `tierA_merged_p95` under §0.1 **D-B** — with the **frozen-literal
fine-level call `tierA_p95` computed and reported on H1 as the declared sensitivity**, every time,
exactly as D-B already requires. §0.1 open item 1 recorded that the PI had not stated which of the
two calls R1 is scored on. **This states it.**

**This resolves an ambiguity; it does not relax a criterion.** Three separate things are unchanged
and are listed separately so the claim can be checked rather than taken:

1. **No threshold moves.** R1 is still "the IQR across the arm's reportable fits includes 0 **and**
   its upper quartile is below 0.50" (as §0.0 **C-4** already corrected the bracket's name), and R2
   is still "the median controlled `|β| / sd(y)` is below **that arm's own** 80 %-power detectable
   bound". Neither number, and neither direction, is touched.
2. **M1's outcome is unchanged**, because M1's primary call is unchanged: SF N2+N5+N6 **0.088**,
   IQR **[−0.017, 0.234]** — includes 0, upper quartile 0.234 < 0.50 — and controlled amplitude
   **0.029** against the 80 %-power bound **0.183** (§0.0 **C-5**, from
   `results/phase3/m1_final_audit.txt` §3).
3. **The moved call is not the easier call.** The merged call **adds** senders where the frozen
   combination left cells eligible but uncallable; it removes none. On the merged `T/NK cells`
   stratum, senders go 23 → **1,958** (SPLN14), **0 → 1,241** (SPLN43) and 295 → **1,675**
   (SPLN44), and A3 goes from failing to passing in those three sections. Of the four sections that
   already passed, three move by at most one cell (SPLN07 1,407 → 1,406, SPLN24 3,658 → 3,657,
   SPLN30 1,900 → 1,900) and SPLN21 rises from 322 senders (1.211 %) to 1,328 (4.996 %) — i.e. the
   merged call adds senders in every section where the fine call was short of the percentile and
   takes none away anywhere. Command:
   `python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/a3_prevalence_by_type.csv'); m=d[(d.label_set=='merged')&(d.cell_type=='T/NK cells')&(d.call.isin(['tierA_p95','tierA_merged_p95']))]; print(m[['section','call','n_all','n_senders_all','prevalence_pct_all','passes_A3']].to_string(index=False))"`.
   A **larger** sender population and a **larger** reportable-fit population is, if anything, more
   able to produce an IQR that **excludes** 0 — i.e. more able to make H1 **fail** R1. **Scoring R1
   on the merged call cannot be characterised as choosing the call that passes.** Which way it
   actually falls is a Phase-10 outcome and is deliberately not stated here.

**What is still not resolved by F1:** whether **M1** is re-called at the merged family for
cross-arm comparability (§0.1 open item 3, and F5 below). F1 fixes how R1/R2 are *scored*; it does
not remove the asymmetry.

---

### F2 — `tierA_merged_p95` is ADDED to the H1 sender axis, and it already runs

**The ruling — add it, never replace an existing call — is applied. The code edit implied by the
flag is rejected, because the addition already exists in the correct place and making it in
`run_phase3_nulls.py` instead would change existing output.** This is recorded as a disagreement
with the letter of the flag and an agreement with its substance.

**Where the addition lives.** The frozen mouse list is untouched:
`run_phase3_nulls.py:70-71` is still the six-call `N7_CALLS`, `:76` the three-call
`TIERA_PM_CALLS`, `:77` `ALL9_CALLS = N7_CALLS + TIERA_PM_CALLS`. The H1 arm **appends** the three
merged calls **after** the nine:

- `code/h1_phase10.py:42-44` — `MERGED_ALIAS = {"tierA_merged_p90": "tierAmg_p90", "tierA_merged_p95": "tierAmg_p95", "tierA_merged_p99": "tierAmg_p99"}`; both spellings are accepted and they are the identical percentile rule at the merged label family, reading the `flag_merged_pNN` flags written by `code/h1_callers.py` and carried into the cache by `code/h1_cache_extend.py`.
- `code/h1_phase10.py:47-62` — wraps `sasp_phase3.Sec.sender_mask` rather than editing the frozen file.
- `code/h1_run_phase10.py:48` — `ALL12 = list(RN.ALL9_CALLS) + MERGED_CALLS`. **Appended, so the nine frozen calls keep index 0…8 and every existing seed.**

Commands: `sed -n '68,78p' code/run_phase3_nulls.py`; `sed -n '42,62p' code/h1_phase10.py`;
`sed -n '46,50p' code/h1_run_phase10.py`.

**Verified that it runs, from output rather than from intent.** All twelve calls, `tierAmg_p95`
among them, are present in the Phase-10 H1 outputs:
`python3 -c "import pandas as pd; [print(f, sorted(pd.read_csv(f, low_memory=False)['call'].unique())) for f in ['results/phase10_h1/window.csv','results/phase10_h1/main_fits.csv']]"`
→ `['cdkn1a_pos', 'senepy_p95', 'senepy_p99', 'tierA_p90', 'tierA_p95', 'tierA_p99', 'tierAmg_p90', 'tierAmg_p95', 'tierAmg_p99', 'tierApm_p90', 'tierApm_p95', 'tierApm_p99']` in both.

**Why adding it to the frozen list instead would have changed existing output — the test the
ruling itself sets.** `_expand` (`run_phase3_nulls.py:123-136`) seeds each job by the **index of
the call in the list it is given**: `sd = int(base) + step_i * i + step_j * j` at `:131`. Because
`ALL9_CALLS = N7_CALLS + TIERA_PM_CALLS`, a seventh entry in `N7_CALLS` moves
`tierApm_p90/p95/p99` from `j = 6, 7, 8` to `j = 7, 8, 9` and **changes the seed of every
per-module null job**. Separately, the mouse cache has no `flag_merged_pNN` array, so
`--calls all` on M1 would raise on the new name. Both are exactly the "if adding it changes any
existing output, stop" condition, and both are avoided by the append-in-the-arm-binding form that
is already in place. *(Operational note, recorded because it is a real constraint and not a
preference: `run_phase3_nulls.py` is imported by the Phase-10 job running at the time of writing,
and by every fresh `loky` worker it spawns, so editing it mid-run is not a safe operation either.)*

**What this does not do.** It does not decide whether the merged call **replaces** `tierA_p95` in
the reported N7 sensitivity axis — it does not: **nine frozen calls plus three merged calls, twelve
in total, all reported.** And it changes nothing on the mouse arm, whose N7 axis is still six calls
(nine with `TIERA_PM_CALLS`).

---

### F3 — every H1 DeepScence **magnitude** is PROVISIONAL at `random_state = 0`

**The resolution.** Every H1 DeepScence **magnitude** in `CS_PHASE9_H1_AUDIT.md` §10 is marked
**PROVISIONAL, pending the five-seed consensus of §0.1 D-A**. All of them descend from the single
seed-0 score files `data/processed_h1/deepscence_h1_<section>.csv`, produced by
`code/h1_deepscence.py:76` — `res = api.DeepScence(A, denoise=False, verbose=False, random_state=0)`
(command: `sed -n '74,78p' code/h1_deepscence.py`) — and consumed by
`code/h1_deepscence_anchor.py:60` and `code/h1_caller_agreement.py:34`.

**Precisely which numbers are PROVISIONAL** (each re-read from the file named):

| quantity | value | file | command key |
|---|---|---|---|
| P-i depth loadings, 7 sections | **+0.1822 … +0.3540** (median +0.2963) | `results/phase9_h1/caller_technical_loading.csv` | `caller == "deepscence_score"`, `spearman_vs_transcript_counts` |
| P-v within-type Q5/Q1, 7 sections | **0.244, 0.403, 0.416, 0.456, 0.496, 0.778, 1.169** | `results/phase9_h1/caller_within_type_depth_bias.csv` | `caller == "deepscence_score"`, `enrichment` Q5 ÷ Q1 |
| anchor ρ partial, `CDKN1A` / `LMNB1` / 8-gene prolif / consensus | **+0.1911 … +0.2540** / +0.2029 … +0.2449 / +0.0097 … +0.0451 / +0.0642 … +0.1695 | `results/phase9_h1/deepscence_anchor_h1.csv` | `rho_partial_*` |
| `CDKN1A` rank among the 33 on-panel CoreScence genes; prolif stability in SPLN21 | 1st – 7th; **0.75** | " | `rank_cdkn1a_in_core`, `stab_prolif` |
| **`CDKN1A` fold-split stability itself** | **1.00 ×7** | " | `stab_cdkn1a` — see the note below |
| every pooled agreement ratio with a DeepScence term | **6.436** (z 204.83), 2.069, 1.602, **1.102** (z 5.67), 1.093, 0.998, 0.290 | `results/phase9_h1/caller_agreement_pooled.csv` | rows containing `deepscence_score` or `abs_deepscence_score` |
| the per-section ranges behind them | DeepScence × `CDKN1A`⁺ **3.459 – 10.115**; Tier A × DeepScence 0.972 – 2.816 (SPLN21 **1.174**, z 3.54) | `results/phase9_h1/caller_agreement_matched_significance.csv` | `ratio_stratified` |

Commands:
`python3 -c "import pandas as pd; print(pd.read_csv('results/phase9_h1/caller_agreement_pooled.csv').to_string(index=False))"`;
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_within_type_depth_bias.csv'); s=d[d.caller=='deepscence_score'].pivot_table(index='section',columns='within_type_depth_quintile',values='enrichment'); print((s['Q5']/s['Q1']).round(3).to_string())"`;
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_technical_loading.csv'); print(d[d.caller=='deepscence_score'].to_string(index=False))"`.

**Precisely which numbers are NOT affected, because the flag as raised is broader than the
evidence.** Marking these PROVISIONAL would be wrong:

1. **P-iv's circularity, 29 / 33 = 0.8788.** It is a **gene-list membership** computation —
   DeepScence's shipped `coreGS_v2.csv` intersected with the panel and with the frozen Tier B
   modules — and **no score is computed and no `random_state` enters it**.
   `python3 -c "import json; print(json.load(open('results/phase7_jobA/gate_result_human.json'))['corescence'])"`
   → `n_on_panel = 33`, `frozen_n_in_any_B = 29`, `frozen_frac = 0.8788`. P-iv is unaffected by
   D-A and stays **CONFIRMED**.
2. **Tier A × SenePy 0.874, z = −7.96, below chance 7/7**, and the SenePy within-type Q5/Q1
   **28.5–228.1** (F7b): no DeepScence term.
3. **Every M1 DeepScence number.** M1 keeps its frozen single-seed estimator (D-A).
4. **The seed-check statistics themselves** — *r* 0.3719, Spearman 0.4157, top-5 % Jaccard 0.2107,
   12,779 cells changing status, and ρ 0.3122 → 0.2308 — which are statements *about* the seed
   pair, not statements at one seed.

**D-C's verdict stands, and one honest qualification is recorded with it.** P-ii is falsified as a
**direction**: the published `CDKN1A` anchor decides the polarity the same way in 20 of 20 folds in
**all seven** sections, against a falsifier written in advance as "all 7 sections have stability
≥ 0.90". **But `stab_cdkn1a` is itself computed on the seed-0 score**, so it is listed above with
the provisional magnitudes. Retaining the verdict rests on the argument stated in §0.1 — that a
five-seed consensus cannot be less stable than its members, and that 1.00 in 7/7 sits two folds in
twenty above the falsifier — and **that argument is asserted, not measured**. It is recorded as an
argument so that the consensus run can check it rather than inherit it.

---

### F4 — P-iii is **MARGINAL**, pending re-score on the consensus

**The resolution.** **P-iii's outcome is not restated in either direction until the five-seed
consensus of §0.1 D-A exists.** Until then it is recorded as **MARGINAL**, not as CONFIRMED.

**The evidence, and why "confirmed" does not survive it.**

- The registered rule (§8) confirms P-iii at **pooled ratio > 1.10** and falsifies it at **≤ 1.05,
  or below chance**. **The interval (1.05, 1.10] is pre-registered as neither** — the prediction
  has an explicit indeterminate band, written before the data.
- H1 returns **1.102**, z = 5.67, above chance in **5 of 7** sections
  (`results/phase9_h1/caller_agreement_pooled.csv`, row `tierA_score` × `abs_deepscence_score`;
  command in F3 above). **That is 0.002 above the top of the band the pre-registration already
  declared indeterminate.**
- The estimator carrying it reproduces across seeds at **top-5 % Jaccard 0.2107** and *r* =
  **0.3719** at full section size (`results/phase9_h1/d2_stability.csv`, row
  `SPLN21, "denoise=False, FULL section", 0, 1`), against an M1 floor of 0.7606 / 0.99553.
- **A margin of 0.002 on that estimator is not a margin**, and the nearest thing to a direct test
  points the same way: the one Tier A × DeepScence matched ratio that was re-derived at another
  seed is the one number this block strikes for having no producer (F7a) — so **there is not even
  an unsourced reassurance left** that the pooled ratio is seed-stable.

**What is unchanged.** The threshold is not moved, the prediction is not rewritten, and the
falsifier is not weakened. **P-iii is simply not scored yet.** It may return CONFIRMED, MARGINAL or
indeterminate on the consensus, and all three are publishable.

---

### F5 — the three cross-arm asymmetries, in one place, with what each costs

**Registered together so that a reader meets the full set at once rather than discovering them one
at a time.** §9 item 4 already required every cross-arm number to be reported twice on panel
grounds; **there are now three axes on which the two arms are not the same measurement.**

| # | Axis | M1 | H1 | Source | **What it costs the cross-arm comparison** |
|---|---|---|---|---|---|
| **1** | **Panel** (pre-freeze, §9 item 4) | 5,097-gene mouse panel, DeepScence ortholog-remapped onto 4,845 of them | native 5,093-gene human panel, no remapping | `genesets/h1_candidate/GSE326743_gene_panel_5093.csv`; `results/phase9_h1/a8_panel_arithmetic.csv`; `results/phase9_h1/a8_ortho_sender_shift.csv` | Any cross-arm difference is confounded with the panel. It is **measurable, and it is large**: cutting H1 to the 2,425-gene ortholog intersection leaves **26 of the 33** Tier A genes and moves the sender set itself — Jaccard between the full-panel and intersected-panel `tierA_p95` sender sets is **0.5222 – 0.5747** across the 7 sections (`a8_ortho_sender_shift.csv`, `n_tierA_on_full`, `n_tierA_on_ortho`, `jaccard`). **Mitigation, already frozen: every cross-arm number twice, on both panels (test A8).** |
| **2** | **DeepScence estimator** (§0.1 **D-A**) | single run, `random_state = 0` | five-seed consensus (20260901–05) with its between-seed dispersion | `code/h1_deepscence.py:76`; §3.8 `COMPMATCH_SEEDS` | **No cross-arm DeepScence magnitude is like-for-like**: the H1 side is an average of five draws and carries a dispersion the M1 side does not have. **Any cross-arm DeepScence difference smaller than H1's own between-seed spread is uninterpretable**, and that spread is not small — at full section size two seeds of the *same* configuration agree at *r* = 0.3719 and top-5 % Jaccard 0.2107 (`results/phase9_h1/d2_stability.csv`). **No mitigation exists on the M1 side**: M1 is not re-run as a consensus, because M1 reproduces at the determinism floor (0.99999913 / 0.99999995, `results/phase8_d2/d2_agreement.csv`) and has nothing to average. The asymmetry is therefore **carried, not removed** — and it is itself a finding (the tool is reproducible on the remapped mouse panel and not on its own native human panel). |
| **3** | **Tier A call level** (§0.1 **D-B**) | fine-label `tierA_p95` | merged-label `tierA_merged_p95` primary, fine call as the frozen-literal sensitivity | `results/phase9_h1/a3_prevalence_by_type.csv` | The two arms' sender populations are defined at **different label granularity**, so a cross-arm sender-prevalence or SF difference partly reflects the call level rather than the biology: on H1 the merged call restores **0–13.6 % of cells per section** that were eligible but uncallable, and on merged `T/NK cells` it moves SPLN43 from **0 senders in 24,815 cells** to 1,241 (§0.1 D-B, F1 above). **Partial mitigation: the fine call is computed on H1 too**, so a fine-vs-fine comparison is always available. **The merged counterpart is NOT computed on M1**, so a merged-vs-merged comparison is not — §0.1 open item 3 stays open. |

**Stated plainly, because the point of the register is that no reader has to assemble it:** a
cross-arm **DeepScence** number carries axes 1 **and** 2; a cross-arm **sender-call** number
carries axes 1 **and** 3; **no cross-arm number in this project is affected by fewer than two of
the three.** Every one of them must be reported with the axes it carries named.

---

### F6 — P7's own justification was already wrong, and that **strengthens** D-D

**The correction.** P7 (§13) argues that dropping `Marginal zone B cells` "would leave the A6 axis
without its middle term". **That is not how the A6 axis is built.** Verified in the source, not
taken from a document:

- `code/h1_a6_compartments.py:57-58` — `COMPARTMENTS = ["D_spleen_red_pulp", "D_spleen_white_pulp_follicle", "D_spleen_white_pulp_tzone", "D_spleen_marginal_zone", "D_spleen_capsule_trabecula"]`. Five compartments are **scored**, marginal zone among them.
- `code/h1_a6_compartments.py:85-86` — `white = 0.5 * (sco["D_spleen_white_pulp_follicle"] + sco["D_spleen_white_pulp_tzone"])` then `pulp = sco["D_spleen_red_pulp"] - white`. **The axis has three terms and `D_spleen_marginal_zone` is not one of them.**
- `grep -n "marginal_zone" code/h1_a6_compartments.py` → **line 58 only**: the label appears in the scored list and **nowhere else in the producer**.

Command: `sed -n '57,59p;85,87p' code/h1_a6_compartments.py; grep -n "marginal_zone" code/h1_a6_compartments.py`.

**The consequence, recorded in the direction it actually runs: retiring P7 costs the A6 axis
nothing.** §0.1 D-D's retirement of P7 is therefore **strengthened**, not weakened — the one
argument P7 offered for keeping the label was already inaccurate **before** any H1 expression value
was read, and the label then turned out never to be realised in any of the seven sections. §0.1
open item 6 flagged this and deliberately did not correct it; **F6 records it as a correction**,
while **P7's original wording stays in place**, exactly as D-D requires. What is corrected is a
**justification**, not a decision and not a number: no threshold, no gene set and no scored
compartment changes, and `genesets/human/D_spleen_marginal_zone.txt` is **not** edited.

---

### F7 — two audit numbers that do not reproduce from a file

**(a) The seed-1 matched ratio 2.967 is STRUCK as unsourced, and no value is substituted.**
`CS_PHASE9_H1_AUDIT.md` §10.5 quotes the Tier A × DeepScence matched ratio moving
1.174 → **2.967 (z 40.0)** between seeds on SPLN21. The seed-**0** value is in a file —
`results/phase9_h1/caller_agreement_matched_significance.csv`, row `SPLN21, tierA_score,
deepscence_score`, `ratio_stratified = 1.174`, `z = 3.54`. **The seed-1 value is in no file under
`results/`, and no producer emits it:** that file has no seed column at all, and the producer of
the full-section seed check, `code/h1_d2_analyse.py:92-116`, writes exactly one seed-1 quantity —
`rho_seed1` — into `results/phase9_h1/d2_depth.csv`. Commands:
`head -1 results/phase9_h1/caller_agreement_matched_significance.csv`;
`head -1 results/phase9_h1/d2_depth.csv`;
`grep -rl "2\.967" results/phase9_h1/` → no match.
**It is struck rather than recomputed: a number that is needed needs a producer, not a
reconstruction.**

**Consequence inside this file, applied rather than left implicit:** §0.1 **D-A**'s own sentence
*"while the Tier A × DeepScence matched agreement ratio moves by a factor of about 2.5"* is derived
from the struck number and **is struck with it** (marked at the site). D-A's conclusion does not
depend on it: the seed-sensitivity evidence that survives is entirely file-backed — *r* 0.3719,
Spearman 0.4157, top-5 % Jaccard 0.2107, 12,779 cells changing status
(`results/phase9_h1/d2_stability.csv`) and the depth loading 0.3122 → 0.2308
(`results/phase9_h1/d2_depth.csv`).

**One further number in the same audit table fails the same test, and is flagged rather than
used.** §10.5's seed-1 within-type Q5/Q1 for P-v, **0.292**, is emitted by no producer either — by
the same `h1_d2_analyse.py:92-116` check. **This file quotes it nowhere**, and it must not be
quoted until a producer writes it. *(Recorded as an extension of the ruling, not as part of it.)*

**(b) SenePy within-type depth enrichment is Q5/Q1 = 28.5 – 228.1×, the committed file's value.**
`CS_PHASE9_H1_AUDIT.md` §10.2 quotes **28.5–224.7**; the committed
`results/phase9_h1/caller_within_type_depth_bias.csv` is rounded to 3 dp and gives **228.1** for
SPLN24 (4.105 / 0.018). Per section: SPLN07 74.6, SPLN14 30.8, SPLN21 50.6, **SPLN24 228.1**,
SPLN30 **28.5**, SPLN43 98.0, SPLN44 70.3. Command:
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_within_type_depth_bias.csv'); s=d[d.caller=='senepy_score'].pivot_table(index='section',columns='within_type_depth_quintile',values='enrichment'); print((s['Q5']/s['Q1']).round(1).to_string())"`.
**The corrected value is the one already used in §14** (the SenePy open item), and it is confirmed
here rather than changed.

**Propagation check, run rather than assumed.** `grep -rn "224\.7" reports/ code/ figures/` returns
**two** hits and no others: `CS_PHASE9_H1_AUDIT.md:765`, the audit's own table — **not editable
from this task, and left for the audit's owner** — and `PREREG_PHASE8.md:1403`, where it appears
**only** as a citation of the audit's value alongside the corrected one. **The wrong value did not
propagate into any claim, any figure or any code path**, and no conclusion anywhere depends on the
difference between 224.7 and 228.1: SenePy's depth enrichment is a two-orders-of-magnitude caveat
either way.

---

**None of F1 … F7 changes a pre-registered threshold, a stop condition, the §18 outcome table, or
any M1 number.** F1 states which call an existing criterion is scored on; F2 records an addition to
the H1 sender axis that was already implemented as an addition; F3 and F4 **withhold** magnitudes
and one verdict pending the consensus, which is strictly more conservative than §0.1 was; F5
assembles three already-declared asymmetries into one register with their costs; F6 corrects a
justification while preserving the wording it corrects; F7 removes one unsourced number and
confirms one corrected number against its file. **No original wording has been deleted anywhere in
this file; every affected site carries its original text with a dated inline marker beside it.**

---

## 0. What is provisional in this draft

Two agents were running while this was first drafted — the **M1 end-to-end re-run (task 8.7)** and
the **D2 denoise-cost measurement (task 8.5)**. **D2 has since landed and is resolved**
(`reports/CS_PHASE8_D2_DENOISE.md`); 8.7 is still in flight. Every number below is marked with
where it came from, and anything 8.7 will move is marked **PROVISIONAL**. Nothing provisional may
be frozen as if it were settled.

| Item | State |
|---|---|
| Frozen tag hash | **`d04691e2692a7be8d1ff676d2fb74ad9d1df049d`** — tag `phase8-frozen`, **RE-CUT 2026-08-27 19:55:58 UTC**. The first cut, `926439629a07269a32c93f998da0f6e1cd20933c` (15:32 UTC), shipped code that could not produce its own results and was replaced; it is an ancestor of the current tag. §1 |
| Composition-matched rerun protocol | **RESOLVED and FILLED (PI decision D15): freeze both variants.** Implemented as `code/run_phase8_compmatch.py`, run on M1, gated on H1. All four previously-`TBD` fields transcribed from `reports/CS_PHASE8_COMPMATCH.md`. **The protocol is inert (1.6 %) where the same variables as covariates remove 85.4 %** — P23. §3.8 |
| Fitting window | **RESOLVED (PI decision D16): frozen as `100 µm, fixed`**, with the 99th-percentile rule recorded as provenance. §3.2 |
| Every pooled number in `results/phase3/caller_coverage_gate.csv` | **PROVISIONAL** — the file currently mixes a pre-C6 2-section base with a post-C6 11-section base. The M1 re-run agent is recomputing both bases under one sender definition and labelling the basis in the output. §13 P1 |
| D2 `denoise` | **RESOLVED (task 8.5): DCA installed and ran; `denoise=False` frozen as a *chosen* value, `denoise=True` as the published-default sensitivity.** No longer provisional and no longer a caveat. §3.9, P22, P26–P28 |

**What 8.7 has and has not yet touched**, by file mtime at the time of writing (2026-08-27
~06:14 UTC). `phase2_downstream.py` rebuilt the sender calls and module scores on the promoted C6
gene sets at **05:48–05:53**, so anything written before then is a **pre-C6** number:

| Artefact | Written | Configuration | Status here |
|---|---|---|---|
| `results/phase3/summary_phase3.txt` (pinned) | 08-20 19:03 | pre-C6 | The **published baseline**, deliberately unchanged. Source of the §5 benchmark values, all marked **PROVISIONAL** |
| `results/phase3/sf_summary_c1.csv` (corrected N3/N4) | 04:41 | **pre-C6** | **PROVISIONAL.** The C1 *conclusion* is invariant across six sender calls spanning 0.5–9.0 % prevalence, so it transfers; the numeric SF values do not |
| `results/phase3/a7_summary.csv` (A7) | 05:19 | **pre-C6 sender calls** | **PROVISIONAL** in magnitude. The three qualitative A7 findings — raw assay not flat, N5 removes it, N2 does not — are what §10.1 and P2 freeze |
| `results/phase3/caller_*_11sections.csv` | 06:07 | **post-C6** | New. See the P1 table in §13 — the headline moved a second time |
| `results/phase3/caller_*_verify2sec.csv` | 06:03 | **pre-C6** | Reproduces the committed published tables. **Do not pool it against the 06:07 files** — see P1 caution (i) |
| `results/phase3/main_fits.csv` | 06:13, still being written | **post-C6** | 8.7 in flight |

Nothing provisional may be frozen as if it were settled, and every provisional row must be re-read
before the tag (§14 item 6).

---

## 1. Frozen tag

```
git tag -a phase8-frozen -m "Frozen for human replication, post C1, C7 and the C6 gene-set decision"
```

**Tag commit hash: `d04691e2692a7be8d1ff676d2fb74ad9d1df049d`** (`phase8-frozen`, **re-cut** 2026-08-27 19:55:58 UTC; supersedes `926439629a07269a32c93f998da0f6e1cd20933c` of 15:32 UTC).

This is deliberately left empty. The tag does not exist at the time of writing and inventing a
hash would be worse than a gap. The PI fills it at the moment of tagging, together with:

| Field | Value |
|---|---|
| Tag name | `phase8-frozen` |
| Tag commit hash | **`d04691e2692a7be8d1ff676d2fb74ad9d1df049d`** (re-cut; supersedes `926439629a07269a32c93f998da0f6e1cd20933c`) |
| Tag date (UTC) | **2026-08-27 15:32 UTC** |
| Predecessor tag | `pre-c6-genesets` — captures the mouse gene sets **before** the C6 sets were promoted into `genesets/`. Already created (PI decision D5). Its hash: **`c002ddda0d1e6402a85bf96abfcd2d3a6165287f`** |

**Nothing in Phase 9 or Phase 10 may begin until this tag and this file are committed.**

---

## 2. Datasets

| Arm | Accession | Tissue | Design | Panel |
|---|---|---|---|---|
| **M1 mouse** | **GEO GSE310392** | liver | 11 Xenium Prime 5K sections, sham / SBR fibrosis | Xenium Prime Mouse 5K pan-tissue + 100-gene custom add-on = 5,106 `Gene Expression` features − 9 genotyping probes = **5,097 genes** |
| **H1 human** | **GEO GSE326743** | **normal spleen** | 7 donors, ages 17 / 31 / 32 / 32 / 37 / 57 / 59, 4 M / 3 F, all normal FFPE | Xenium Prime 5K Human standard + 100 add-on, **5,093 genes**, verified on the `.h5` feature table, byte-identical across SPLN07 / SPLN30 / SPLN44 |

H1 sample-level detail, the 132-series GEO screen that selected it, and the 19 Prime-5K candidates
it was chosen from: `reports/PHASE7_H1_SCREEN.md`.

**H1 sections and cells:** 7 sections, 2,207,593 cells (SPLN07 249,420 / SPLN14 329,371 /
SPLN21 220,435 / SPLN24 396,173 / SPLN30 366,199 / SPLN43 331,582 / SPLN44 314,413).
Source: `reports/PHASE7_H1_SCREEN.md` §3, from the `.h5` shapes at acquisition.

**Control features present on both panels** (Prime 5K stock complement, confirmed on the data of
both arms): 609 negative control codewords, 40 negative control probes, 21 genomic controls. Audit
test A7 needs the 40 negative control probes and they are there.

---

## 3. Fixed parameters — not retunable after the tag

Every value in this section was read from the code or from a result file on disk, with the file
named. Where the planning document's wording and the implementation differ, the **implementation**
is what is frozen and the difference is declared.

### 3.1 λ grid

`code/run_phase3_nulls.py:86-93` (`lam_grid`), constants at `:59-62`.

```python
WINDOW_UM    = 100.0
LAM_LO_FLOOR = 7.0     # the resolution floor: at or below the median NN distance of every section
N_LAM        = 40
def lam_grid(dmax=WINDOW_UM, med_nn=0.0):
    return np.exp(np.linspace(np.log(LAM_LO_FLOOR), np.log(dmax / 2.0), N_LAM))
```

- **Floor = 7.0 µm**, the resolution floor. Median nearest-neighbour distance is 6.7–10.6 µm across
  sections, so 7 µm is at or below it for most but **not** all of them.
- **Ceiling = window / 2 = 50.0 µm.** An exponential with λ > dmax/2 is not distinguishable from a
  linear trend over [0, dmax].
- **40 log-spaced points.** Confirmed in `results/phase3/main_fits.csv`: `lam_grid_lo` ≡ 7.0 and
  `lam_grid_hi` ≡ 50.0 on every row.
- **Declared:** `lam_grid`'s `med_nn` argument is accepted and never used. The floor is a literal
  and does **not** adapt per section. This is frozen as-is; it must not be quietly "fixed" for H1,
  because that would make the two arms non-comparable. If H1's median NN distance is materially
  different from M1's, that is reported as a limitation, not patched.
- **"Railed at a grid bound"** is frozen as the Phase-3 definition, `run_phase3_nulls.py:239`:
  `lam_railed = int(t0 == 0 or t0 == sf.lam.size - 1)` — exact index equality at either end.
  *Declared:* three other non-equivalent definitions exist in the codebase
  (`run_phase5_kernels.py:124` exact value equality; `:361` a 0.1 % tolerance for the proximal-vs-
  downstream control; `sasp_estimators.py:155` a 2 % tolerance, synthetic only). Any railing rate
  quoted in the paper is the Phase-3 one unless it says otherwise.

### 3.2 Fitting window

**Frozen value: 100 µm, applied as a hard cap on receiver distance in every fit.**
`code/run_phase3_nulls.py:59` (`WINDOW_UM = 100.0`), enforced at `:169`
(`& np.isfinite(self.d_obs) & (self.d_obs <= WINDOW_UM)`), and again at `:585`,
`run_phase3_n8.py:224`, `run_phase3_lamscale.py:52`, `make_phase5_figs.py:94`; re-declared as the
same literal in `phase4_data.py:36` and `phase3_null_diag.py:31`; it is also COMMOT's distance
threshold (`phase4_run.py:61`) and the spline `dmax` (`run_phase5_kernels.py:79, 213`).

**RESOLVED by the PI (decision D16): freeze the code's actual behaviour.** §15 says
"window = 99th percentile of distance-to-nearest-sender". **The code does not compute a percentile
at run time**, and it never has. The 99th-percentile rule is recorded here as the **provenance** of
the literal — it is how 100 was chosen, once, from the six Test-3-admissible sections at the
primary sender call — and **not** as a runtime rule. **What is pre-registered is `window = 100 µm,
fixed`.** Every published result already used 100 µm, so freezing the literal preserves continuity
and re-runs nothing.

*(Namespace note: "D16" here is the **PI decision** series of `PHASE8_ROADMAP_STATUS.md`, not the
gene-set deviation D16 of `PREREG_PHASE8_genesets.md` §12, which is the `kneed`/`openpyxl` pin.)*

The percentile measurement itself is `stage_window`
(`run_phase3_nulls.py:100-136`, `d_p99 = float(np.percentile(dr, 99))` at `:128`), written to
`results/phase3/window.csv` (297 rows = 11 sections × 27 section/call/module combinations).

Measured 99th percentiles from that file:

| sender call | in-band (6 sections) range | median over all 11 |
|---|---|---|
| `tierA_p95` (**primary**) | **76.0 – 112.1 µm** | 96.1 |
| `cdkn1a_pos` | 76.3 – 160.9 µm | 95.7 |
| `senepy_p95` | 118.3 – 186.8 µm | 151.0 |
| `tierA_p90` | — | 66.3 |
| `tierA_p99` | — | 230.4 |
| `senepy_p99` | — | 354.6 |

So 100 µm sits close to the in-band `tierA_p95` median (96.1) and **truncates** the window for the
sparser calls. **This is a real limitation of the SenePy caller under the frozen window and is kept
in the deviation table (P8) rather than buried here: under `senepy_p95` all six in-band 99th
percentiles — 118.3 to 186.8 µm — exceed the 100 µm cap.** Quantified from the same file over those six sections, the
share of receivers the cap discards is bracketed by `frac_gt_80` and `frac_gt_150`: **7.4–21.5 %
lie beyond 80 µm and 0.2–2.3 % beyond 150 µm under `senepy_p95`**, against **0.8–7.0 %** and
**0.0–0.11 %** under the primary `tierA_p95`. A SenePy-called kernel is therefore fitted on a
materially truncated distance distribution, which bounds what it can claim about reach. **Pre-registered consequence:** the window is a fixed 100 µm on both arms,
it is not re-derived from H1, and any SenePy-called H1 quantity carries **two** caveats — this
truncation, and the cross-tissue surrogate-hub problem of §3.7. `results/phase3/window.csv` is regenerated for H1 and reported, but does not change the cap.
*(The M1 copy of `window.csv` was being rewritten by the 8.7 re-run as this was written; the ranges
above are from the file as it stood and are **PROVISIONAL** to the extent the new sender calls move
them. The 100 µm cap itself is a literal and does not move.)*

### 3.3 Kernel families and the selection rule

**Families (5), frozen:** `code/sasp_kernels.py:81` and `code/run_phase5_kernels.py:52` —
`exponential`, `gaussian`, `powerlaw`, `step`, `spline`.

| family | form | source |
|---|---|---|
| exponential | `exp(-d/λ)` | `sasp_kernels.py:48-78` |
| gaussian | `exp(-0.5 (d/λ)²)` | " |
| powerlaw | `(1 + d/λ)^(-p)`, `p ∈ {0.5, 1.0, 2.0, 4.0}` (`run_phase5_kernels.py:53`) | " |
| step | `1[d < λ]` | " |
| spline | cubic B-spline, `N_KNOTS = 6` (`run_phase5_kernels.py:54`), knots at data quantiles, `dmax = WINDOW_UM` | " |

**Selection rule, frozen: strict argmin AIC. No BIC, no within-section cross-validation, and no
ΔAIC tolerance band.**
`run_phase5_kernels.py:176` — `best = min(FAMS, key=lambda f: res[f]["aic"])`, with
`aic = n·log(rss/n) + 2(k_par + 1)` at `:122`. `d_aic_vs_best` is recorded per row (`:191-192`);
`summarize_phase5.py:215` scores a family's "win" rate as `(d_aic_vs_best == 0).mean()`, i.e. exact
ties only. If a ΔAIC < 2 band was ever intended, it is not in the code and is **not** frozen.

**Second, independent and separately reported criterion:** leave-one-**section**-out held-out
Gaussian log-likelihood (`stage_heldout`, `run_phase5_kernels.py:203-264`; core in
`phase5_common.py:146-186`; λ chosen on the training folds at `phase5_common.py:168`; winner
`pp.idxmax(axis=1)` on `ll_per_cell`, `summarize_phase5.py:236`). Both criteria are run under
`design="naive"` and `design="ctrl"` (the full N5+N6 design), `run_phase5_kernels.py:169-170`.
**Per §3.6 below, only the `ctrl` result is reportable as a distance effect on this assay.**

### 3.4 Null battery N1–N8, with the corrected N3/N4

| Null | Definition as implemented | Source |
|---|---|---|
| **N1** | Sender-label permutation, stratified by cell type, among sender-eligible cells | `run_phase3_nulls.py:379-386` (`permute_within_type`), invoked `:464-465` |
| **N2** | Matched decoys: greedy 1-1 propensity matching within section and within cell type, caliper **0.25 SD**; contrast is the shared-λ two-kernel fit `fit2_shared` | `phase3_core.py:26-108`; `run_phase3_nulls.py:159-160, 251-253, 265-267` |
| **N3** | Torus shift of the sender coordinate set — six variants, below | `run_phase3_nulls.py:389-399`; `phase3_null_geom.py` |
| **N4** | Rotation of the sender set about its centroid — five variants, below | `run_phase3_nulls.py:402-408`; `phase3_null_geom.py` |
| **N5** | Nuisance-covariate conditioning: 11 named columns + k-NN cell-type composition + 2 segmentation dummies. Columns: `log_counts, log_genes, log_area, log_nucarea, log_dens25, log_dens50, log_dens100, nn1_um, zonation, zonation_sq, log_dist_boundary` | `phase3_core.py:146-154` |
| **N6** | Receiver neighbour baseline: mean response over the k = 20 nearest neighbours, excluding senders and self | `phase3_core.py:164-173` |
| **N7** | Sender-definition sensitivity axis: the same fits repeated over the sender calls | `run_phase3_nulls.py:68-69` (`N7_CALLS`), extended with `TIERA_PM_CALLS` at `:77-78` ***[on H1 extended again 2026-08-27 — see the F2 marker under this table and §0.2]*** |
| **N8** | Three products: gene-set disjointness; scrambled-response (real module vs **200** expression-matched random gene sets, same cells, same λ); CoreScence/DeepScence circularity | `run_phase3_n8.py:151-240`, `N_RAND = 200` at `:36` |

***[Post-freeze addition 2026-08-27 — see §0.2 item **F2**. **On H1 the N7 axis is the nine frozen calls PLUS the three merged-label calls — twelve in total. Nothing is replaced and nothing is dropped.** `tierA_merged_pNN` (spelled `tierAmg_pNN` in the Phase-10 output files; both spellings are accepted) is the identical Tier A percentile rule at the merged label family, and under §0.1 **D-B** `tierA_merged_p95` is the H1 PRIMARY. **The frozen mouse list above is unchanged and is deliberately not edited**: `N7_CALLS` is still the six calls at `run_phase3_nulls.py:70-71` and `ALL9_CALLS = N7_CALLS + TIERA_PM_CALLS` at `:77`. The addition is made in the H1 arm binding — `code/h1_phase10.py:42-62` (alias + `sender_mask` wrapper) and `code/h1_run_phase10.py:48`, `ALL12 = list(RN.ALL9_CALLS) + MERGED_CALLS`, **appended after the nine so every existing call keeps its index and therefore its seed** (`_expand` seeds by call index, `run_phase3_nulls.py:131`). Verified runnable from output, not from intent: all twelve calls, `tierAmg_p95` among them, appear in `results/phase10_h1/window.csv` and `results/phase10_h1/main_fits.csv`. **The M1 N7 axis is untouched.**]***

**Corrected N3/N4 variant identifiers, frozen exactly as these strings**
(`phase3_null_geom.py:57-61`):

```
TRANSLATION = ("N3_orig", "N3_tile", "N3_occ", "N3_occ15", "N3_swap", "N3_snap")
ROTATION    = ("N4_orig", "N4_tile", "N4_occ", "N4_occ15", "N4_swap")
```

- `*_orig` — the published bounding-box wrap (`:193-197`). It leaves 23 % of N3's shifted senders
  and 8 % of N4's outside the tissue (`results/phase3/null_destructiveness.csv`).
- `*_occ` — accept an offset/angle only if ≥ 95 % of shifted senders land in an occupied 25 µm grid
  cell (`OCC_TOL = 0.05`, `:52`; `:207-216`). **Degenerate on a liver section** and reported as
  such, not as a corrected null.
- `*_occ15` — the same at ≥ 85 % (`OCC_TOL_RELAXED = 0.15`, `:53`).
- `N3_swap` — relocate senders to random real cell positions without replacement (`:235-240`).
  **This is a label permutation, not a torus shift**, and reproduces N1; see §3.5 and the footnotes
  in `reports/CS_PHASE8_C1_CLOSEOUT.md` §4.1.
- `N3_snap` — `orig`, then snap each sender to the nearest real cell (`:246-247`).
- `*_tile` — wrap inside the Phase-4 solid-tissue tiles, `TILE_SOLID = 0.98`, `GRID_UM = 25.0`
  (`:51-54, 258-273`). **These are the clustering-preserving in-tissue variants and are the ones to
  cite as the corrected N3/N4.**
- **Declared:** `N4_snap` appears in the header prose of `phase3_null_geom.py:33` but is **not** in
  the `ROTATION` tuple and is not implemented; `N4_swap` is itself defined as rotate-then-snap
  (`:249-252`). The rotation family has five variants, not six.

Permutation count for the corrected variants: **1,000** (`results/phase3/sf_summary_c1.csv`,
`n_perm` column, 13 rows). *Declared inconsistency:* `summarize_phase3_c1.py:191-192` prints a
header saying "200 permutations" for the N7 sub-table; the data column says 1000 for the primary
call and the N7 axis was run at 200 (`reports/CS_PHASE8_C1_CLOSEOUT.md` §7). The **data column**
is authoritative per row.

### 3.5 Surviving fraction

**Two definitions, by null family. Both are frozen; which one applies is determined by the null,
not by the analyst.**

*Conditioning nulls* (N2, N5, N6, zonation) — ratio of β̂ at the **same** λ index `t0`, where `t0`
is the naive λ̂ (`run_phase3_nulls.py:236, 256-270`):

```python
out[f"sf_{k}"] = v / b["base"]
```

*Perturbation nulls* (N1, N3, N4, and the corrected variants) — `run_phase3_nulls.py:498`, and
identically at `:678`:

```python
o[f"{nm}_sf"] = float((rec["beta_obs"] - v.mean()) / rec["beta_obs"])
```

*N8 scrambled* — `run_phase3_n8.py:203`: `sf_n8 = (bstd - br.mean()) / bstd`.

Bootstrap CI on SF is the ratio of **paired** bootstrap draws (`:297-304`).

**The "reportable fit" filter — the denominator of every SF quoted:** a fit is reportable iff its
naive amplitude is positive **and** its spatial block-bootstrap CI excludes zero.
`summarize_phase3.py:85` — `rep = d[(d.beta_naive > 0) & (d.beta_base_lo > 0)]`; identically
`summarize_phase3_c1.py:50-52`. A prior floor also applies: **`MIN_RECEIVERS = 2000`**
(`run_phase3_nulls.py:65`, enforced at `:217`).

On M1 this gives 315 fits and **160 reportable** at the primary call
(`results/phase3/main_fits.csv`, in-band × `tierA_p95` × `stratum == "all"`; the same 315/160
appears in the pinned `results/phase3/summary_phase3.txt`). **PROVISIONAL** — 8.7 recomputes it.

***[Corrected 2026-08-27 — see §0.0 item C-11 (and C-5). 8.7 did recompute it: the count is
**315 fits / 153 reportable**, re-derived from `results/phase3/main_fits.csv` and printed at
line 19 of `results/phase3/summary_phase3.txt`. The "same 315/160 appears in the pinned file"
claim is **withdrawn**: that file was overwritten at 2026-08-27 09:06 and now reads 153.]***

**Pre-registered consequence, from A7 (§3.6): the reportable-fit filter's nominal 5 % is wrong on
this assay. Its measured false-positive rate is 9–16 %.**

### 3.6 Bootstrap

**400 replicates over 100 quantile blocks.** `run_phase3_nulls.py:63-64` —
`N_BLOCKS_SIDE = 10`, `N_BOOT = 400`; blocks built at `:155-156`, `self.nb = N_BLOCKS_SIDE ** 2`
= 100; resampling at `:276-277` (`rng.multinomial`), CI at the 2.5 / 97.5 percentiles (`:295-296`).

**What the blocks are quantiles of** — `sasp_phase3.py:290-296`: a 10 × 10 grid of **marginal
quantiles of the x and y cell-centroid coordinates**, computed over **all cells of the section**,
so blocks are equal-count and none is empty.

**Declared:** per fit, blocks are re-indexed to populated blocks only
(`run_phase3_nulls.py:226-227`), so the effective block count is ≤ 100 for cell-type-restricted
fits. Other stages use different counts and are not covered by this parameter:
`run_phase5_kernels.py:51` `N_BOOT = 200`; `run_phase5_wc.py:57` `N_BOOT = 200`, `N_SPLIT = 20`;
the proximal-vs-downstream control draws 2,000 donor blocks with **sections** as blocks
(`run_phase5_kernels.py:366`).

**A7's finding on the bootstrap, pre-registered as a limit on interpretation.**
`results/phase3/a7_summary.csv`, design `n6n5`, `frac_CI_excludes_zero` over the five control
responses: **0.091, 0.103, 0.109, 0.127, 0.164** — i.e. **9–16 % against a 5 % nominal**, on
responses whose true amplitude is known to be zero. The reportable-fit filter therefore admits two
to three times more fits than its nominal rate implies, in both arms. This is stated in Methods,
not in a footnote.

***[Corrected 2026-08-27 — see §0.0 items C-2 and C-3. The list is pre-C6: the frozen values are
0.091 / 0.103 / 0.109 / **0.145** / 0.164, and the range 9–16 % is unchanged. The sentence "the
reportable-fit filter therefore admits two to three times more fits than its nominal rate
implies" is **withdrawn** — it is a forbidden claim (audit R6). The filter admits 3.0–13.3 % on
the naive design and 4.8 %, essentially nominal, on the full N6+N5 design; 9–16 % is the
estimator's two-sided CI-exclusion rate, which is a different statistic.]***

### 3.7 Sender callers and their thresholds

**§15 says "the four sender callers and their thresholds". The phrase is ambiguous in the code and
the pre-registration resolves it explicitly, because two different things are called a caller.**

**(a) Sender *masks* that enter the distance fits — three score families, not four**
(`sasp_phase3.py:233-271`). DeepScence is **not** among them.

| call | rule | granularity | source |
|---|---|---|---|
| `tierA_pNN`, **`tierA_p95` is PRIMARY** | strict `>` NNth percentile of `tierA_score`, **within cell type**, per section; cell types with ≥ 20 cells | per section × cell type | `sasp_phase3.py:256-257`; flags built at `phase2_downstream.py:98-105` for `q ∈ {90, 95, 99}` ***[M1 primary; on H1 superseded 2026-08-27 — see the D-B marker under this table and §0.1]*** |
| `cdkn1a_pos` | `cdkn1a_counts > 0`. No percentile, no stratification | per cell | `sasp_phase3.py:258-259` |
| `senepy_pNN` | strict `>` NNth percentile of the SenePy hub score, **within cell type**, per section; ≥ 100 finite-score sender-eligible cells | per section × cell type | `sasp_phase3.py:260-268` |
| `tierApm_pNN` (**pre-registered sensitivity**, PI decision D1) | the same rule scored on `A_sender_for_<module>.txt`, so the mask depends on the response module | per section × cell type × module | `sasp_phase3.py:247-255`; flags at `phase2_downstream.py:107-119` |

***[Post-freeze decision 2026-08-27 — see §0.1 item **D-B**. "Within cell type" in the `tierA_pNN` row above is the **fine** label family, while the estimator stratifies receivers and builds the sender-eligibility mask on **merged** labels (`sasp_phase3.LABELS = "merged"`), leaving **0–13.6 % of cells per section eligible but uncallable** — including **0 senders in 24,815 cells** for merged `T/NK cells` in SPLN43. **On H1 the PRIMARY call is the identical percentile rule applied at the merged family, `tierA_merged_p95`**; the frozen fine-label `tierA_p95` call is retained and reported as the **frozen-literal sensitivity**. No threshold is tuned. This row is unchanged for M1.]***

`tierA_score` is `scanpy.tl.score_genes(..., ctrl_size=200)` on `A_SENDER_FINAL_strict`
(`phase2_downstream.py:73`). All masks are finally `& ok`, excluding `Low_quality`, `Unknown`,
`unknown` (`sasp_phase3.py:33`) and `Proliferating` (`:35`).

**N7 axis, frozen:** `N7_CALLS = ["tierA_p90", "tierA_p95", "tierA_p99", "cdkn1a_pos",
"senepy_p95", "senepy_p99"]` (`run_phase3_nulls.py:68-69`) — **six** calls, extended to nine with
`TIERA_PM_CALLS` (`:77-78`) under decision D1. The freeze does **not** claim four.

**(b) The four *scores* in the caller-agreement analysis** —
`caller_disagree.py:41`, `caller_disagree_all.py:37`:
`['tierA_score', 'senepy_score', 'deepscence_score', 'cdkn1a_counts']`. Three thresholding rules,
all top-5 %, all per section: global (`caller_disagree.py:60-62`); within cell type, strata ≥ 50
cells (`:65-68`); and **depth- and type-matched** — recomputed inside each
(cell type × within-type transcript-count **decile**), `N_DEPTH_DECILES = 10`, strata ≥ 50 cells
(`caller_disagree_all.py:83-109`). **The depth- and type-matched rule is the primary one for every
agreement number.**

**SenePy on H1 is not the same estimator as SenePy on M1.** SenePy 1.0.1 ships 65 hubs across 10
tissues and **spleen is not among them**; liver is (`results/phase7_jobA/senepy_spleen_coverage.json`:
`spleen_hub: false`, `n_hubs: 65`, `n_labels: 22`, `n_usable_surrogate: 15`, `n_no_hub: 7`,
`min_on_panel: 10`). Of the 22 spleen labels, **0 get a tissue-matched hub, 15 get a cross-tissue
surrogate, and 7 get no SenePy score at all** (cDC1, cDC2, pDC, lymphatic endothelium, erythroid
cells, megakaryocytes, mesothelial cells). See §13 deviation P4 and
`reports/PREREG_PHASE8_genesets.md` §13.

### 3.8 Composition-matched rerun protocol at 5 seeds — **FROZEN, BOTH VARIANTS**

**PI decision D15, taken under delegated authority: freeze the matched-decoy protocol *and* its
covariate-adjusted counterpart.** §15's item now reads "the composition-matched rerun protocol at
5 seeds **and its covariate-adjusted counterpart**". Implemented and run on the mouse arm;
built, arm-generic and **gated** on the human arm. Full specification, provenance and results:
**`reports/CS_PHASE8_COMPMATCH.md`**.

*(Namespace note: "D15" here is the **PI decision** series of `PHASE8_ROADMAP_STATUS.md`, not the
gene-set deviation D15 of `PREREG_PHASE8_genesets.md` §12 — the `genesets/mouse_c6/` promotion,
itself superseded by P10.)*

**Which is primary.** The **covariate-adjusted counterpart is PRIMARY** for any claim about how
much of the gradient is composition. The **matched-decoy protocol is reported alongside**, every
time, because §15 specifies it and **its inertness is itself a finding** (P23). This is the same
report-both discipline already applied to the two Tier A variants (P11) and to B7 (§4 of the
gene-set file).

| Field | Frozen value |
|---|---|
| **Producer script** | **`code/run_phase8_compmatch.py`** (driver `code/_compmatch_chain.sh`) |
| **Matching variables** | Exact stratum: **receiver cell type**. Propensity score: the **20-NN cell-type composition vector** (`knn_frac_<type>`, one column per cell type present in the section, taken by name from `phase3_core.build_blocks`), plus cell-type dummies in the propensity model. **No density, depth, anatomy or segmentation term** — adding them is the `full` comparison variant, not the protocol. On M1 the matching set is 9 columns: `knn_frac_{B-cells, Biliary/ductular, DC, Endothelial, Hepatocytes, Macrophages, Mesenchymal, Proliferating, T/NK cells}` (`compmatch_reruns.csv`, `matched_on`) |
| **Matching rule** | 1-1 nearest-neighbour on the propensity logit, **without replacement**, **within cell type and within section**, **caliper 0.25 SD** of the within-stratum score sd (`phase3_core.greedy_ps_match`). Estimator: shared-λ two-kernel sender-vs-decoy contrast (`BlockProfiler.beta2_at`) at **λ = λ̂_naive**; SF = β_sender / β_naive; **400 block-bootstrap replicates over 100 quantile blocks** (§3.6), window 100 µm (§3.2) |
| **Number of seeds** | **5** |
| **The five literal seed values** | **20260901, 20260902, 20260903, 20260904, 20260905** — chosen outside the range `run_phase3_nulls._expand` can reach, so no rerun accidentally reproduces an existing N2 match |
| **What the seed controls** | Both the greedy matcher's claim order and the block bootstrap, **but point estimates are computed before the bootstrap**, so the reported between-seed spread is matching variability only |
| **Covariate-adjusted counterpart (PRIMARY), 3 seed-free variants** | `type_adj` — receiver's own cell-type intercepts; `comp_adj` — the 20-NN composition vector as covariates; `typecomp_adj` — both. Same scopes, same λ̂, same bootstrap. `typecomp_adj` is the exact regression counterpart of the matching set |
| **Comparison variant** | `full` — the published N2 matching set (adds `log_dens50`, `log_counts`, `zonation`), run at the same five seeds, so the isolation of composition is measured rather than asserted |
| **Scope** | **Pooled (unstratified) is primary** — it is the fit the §17 composition row is about; per receiver cell type is secondary |
| **Sections** | `sasp_phase3.IN_BAND` — the six §8 Test 3 admissible sections |
| **Sender calls** | `tierA_p95` on `A_SENDER_FINAL_strict` (**PRIMARY**) and `tierApm_p95` on the seven per-module sets (the D1 sensitivity) |
| **Reportable population** | `beta_naive > 0 & beta_base_lo > 0`, identical to `summarize_phase3.sf_table` (§3.5) |
| **Balance gate** | Master Plan §8 Test 5, \|SMD\| ≤ 0.1. On M1: max \|SMD\| **0.0916 → 0.0352**, gate passes in **100 %** of matches |
| **Output** | `results/phase3/compmatch_reruns.csv` (per-seed + summary), `results/phase3/compmatch_fits.csv` (fit level) |
| **Arms** | **Both.** M1 run 2026-08-27. **H1 is gated**: `--arm h1` refuses unless `ARMS['h1']['sections']` is populated *and* `SASP_H1_UNFROZEN=1` is set. Running it on H1 is roadmap item 10.2, after Phase 9 |

**Frozen M1 values** (`results/phase3/compmatch_reruns.csv`, `row_type = summary`,
`call = tierA_p95`, `scope_kind = pooled`; independently re-read from the file for this
pre-registration):

| Variant | n reportable | SF | share removed | 95 % CI on SF | across-seed range (sd) |
|---|---|---|---|---|---|
| **`comp` — the protocol** | 165 | **0.9837** | **1.6 %** | [0.973, 0.994] | 0.98370 – 0.98397 (1.2 × 10⁻⁴) |
| `full` — published N2 set | 165 | 0.9855 | 1.4 % | [0.979, 0.992] | 0.98539 – 0.98577 (1.6 × 10⁻⁴) |
| `comp_adj` | 33 | 0.4989 | 50.1 % | [0.421, 0.606] | seed-free |
| `type_adj` | 33 | **0.3415** | **65.9 %** | [0.236, 0.402] | seed-free |
| **`typecomp_adj` — PRIMARY** | 33 | **0.1461** | **85.4 %** | [0.052, 0.246] | seed-free |

Median match rate **0.999871**; the five per-seed pooled SFs are 0.983695 / 0.983887 / 0.983965 /
0.983714 / 0.983735.

**The harness is the published estimator.** `beta_naive` and `lam_naive` are **bit-identical** to
`results/phase3/main_fits.csv` on all 142 cells reportable in both (max |Δβ| = 0.000e+00), and the
`full` variant reproduces the published N2 (five-seed medians 0.9463–0.9499 against 0.9490).

### 3.9 DeepScence run settings

Frozen exactly as currently invoked. `run_deepscence_all.py` (11 sections) and `run_deepscence.py`
(the preserved 2-section run) are byte-identical in settings.

| Setting | Value | Source |
|---|---|---|
| denoise | **`denoise=False` — PRIMARY, and a *chosen* value, not a forced one.** `denoise=True` (the published default) is the pre-registered sensitivity. See the D2 block below | `run_deepscence_all.py:60`, `run_deepscence.py:40`; `reports/CS_PHASE8_D2_DENOISE.md` |
| `random_state` | `0` | same lines ***[M1 only after 2026-08-27 — see the D-A marker under this table and §0.1]*** |
| Anchor | **published `CDKN1A`** (DeepScence's own `io.py::fix_score_direction`), not overridden | `run_deepscence_all.py:12-13` |
| Minimum counts per cell | **≥ 20** | `run_deepscence_all.py:54`, `run_deepscence.py:35` |
| Panel, M1 | ortholog-remapped, MGI 1:1, **4,845 of 5,097** panel genes map | `run_deepscence_all.py:32-34, 55-56`; `logs/ds_smoke.log` |
| Panel, H1 | **native — no remapping.** This is §8's experiment | §8, and `reports/PREREG_PHASE8_genesets.md` §5 |
| Coverage | **11 / 11 M1 sections** (was 2 / 11), 1.47 M cells | task C7/D1 |

***[Post-freeze decision 2026-08-27 — see §0.1 item **D-A**. The `random_state` row above is unchanged for **M1**, which reproduces at the determinism floor (*r* = 0.99999913 / 0.99999995, `results/phase8_d2/d2_agreement.csv`). On **H1** the single `random_state = 0` run is **no longer the primary estimator**: it reproduces at *r* = 0.3719 and top-5 % Jaccard 0.2107 across seeds at full section size (`results/phase9_h1/d2_stability.csv`), against an M1 seed-to-seed floor of 0.99553 / 0.7606. The H1 primary is a **five-seed consensus at the frozen seeds 20260901–05, reported with its between-seed dispersion**. The two arms therefore use **different DeepScence estimators**, and that asymmetry is reportable. **The four attributes every DeepScence number carries become five**: coverage, denoise state, anchor, panel — and the seed configuration with its dispersion.]***

**The four attributes every DeepScence number carries** (§9 reporting standard, one Methods row,
not scattered caveats): **coverage, denoise state, anchor, panel (native or ortholog-mapped, with
the mapping rate)**.

**D2 (`denoise`) — RESOLVED. DCA installed, ran, and lost on the merits.**

**§6 path 1 landed; the fallback was not needed.** DCA 0.3.4 runs under TensorFlow 2.4.4 / Keras
2.4.3 in an isolated CPython **3.8.19** venv, out of process, handing a denoised matrix back to the
pinned 3.11 stack over plain `.npz`/`.npy`. **The main environment was not modified and has no
TensorFlow in it**, and every committed `deepscence_*.csv` is byte-identical. `denoise=True`
completed on **three full M1 sections, both arms** — 7239 and 7259 (SBR) and 7352 (sham) — plus a
20,000-cell three-seed panel. Recipe: `code/setup_dca_env.sh`; rebuild ≈ 4 minutes.

**The sentence this pre-registration used to carry — "DeepScence as we can run it on this panel,
not DeepScence as published" — is withdrawn. We ran it as published and chose against it.**

**Reason 1 — §4 (D-b)'s premise is refuted by measurement.** The planning document states that DCA
denoising "is precisely the step that would normalize depth — the confound under investigation".
**It is not: denoising RAISES the caller's depth loading by x1.32-1.67**, on three of three sections and
on both arms (`results/phase8_d2/d2_depth.csv`, Spearman ρ of score against transcript counts):

| section | arm | `denoise=False` | `denoise=True` | Δρ |
|---|---|---|---|---|
| 7239 | SBR | 0.3891 | **0.6404** | **+0.2512** |
| 7259 | SBR | 0.3176 | **0.5314** | **+0.2138** |
| 7352 | **sham** | 0.4096 | **0.5419** | **+0.1323** |

It also collapses the global top-5 % sender set onto **100.0 % hepatocytes on all three sections**
(from 64.5 / 71.5 / 97.2 %), with biliary/ductular enrichment falling to ~0, and it changes which
cells are senders: sender-set Jaccard against the committed run is **0.118–0.280**, i.e. **56–79 %
of the committed sender set is not called**.

**Reason 2 — the published default is not reliably reproducible.** One fixed 20,000-cell
subsample, three seeds, nothing else changed (`results/phase8_d2/d2_stability.csv`):

| pair | Pearson *r* | top-5 % Jaccard | cells changing status |
|---|---|---|---|
| `denoise=False`, seed 0 vs 1 | **0.9955** | **0.761** | 272 |
| `denoise=True`, seed 0 vs 2 | 0.9824 | 0.665 | 402 |
| **`denoise=True`, seed 0 vs 1** | **0.5703** | **0.000** | 2,000 |
| **`denoise=True`, seed 1 vs 2** | **0.5732** | **0.000** | 2,000 |

**One of three seeds returned a top-5 % sender set sharing no cells at all with the other two**,
with no diagnostic that anything was different. §10.10 and P26.

**Stated against interest: `denoise=False` is not the de-confounded option either.** It is simply
*less* depth-loaded than the published default. An explicit library-size normalisation (`lib`,
below) cuts the caller's depth loading by 74 % and 93 % and changes 46 % and 100 % of its sender
calls. **The frozen caller is depth-loaded.**

**Also frozen: the two §6-path-3 normalisation configurations, with `mor`'s inadequacy stated.**
§6 names median-of-ratios (*poscounts*); measured on this sparse 5K panel it removes only
**11.4–23.6 %** of log-depth variance (p90/p10 depth ratio 11.0 → 8.7–10.0,
`d2_normalisation_strength.csv`) and returns a clean null on four sections — **a null about a weak
normalisation, not about the caller.** The `lib` configuration removes **100 %** (p90/p10 → 1.000)
and moves the caller a great deal: ρ 0.3176 → **0.0819** on 7259 and 0.4096 → **0.0289** on 7352
(−74 %, −93 %), with 46.4 % and 100 % of the committed sender set no longer called, and the call
set **broadening** across cell types rather than concentrating. Both configurations are in the
frozen list; **`mor` may not be reported as evidence that normalisation cannot move this caller.**

**Determinism control — a positive result, and the floor everything else is read against.** Re-run
at the same seed, the harness reproduces the committed scores at Pearson *r* = **0.9999991**
(7239, 24 of 75,384 calls moving) and **0.9999999** (7259, 2 of 114,721). 7259's committed score is
`deepscence_sbr.csv`, **written 2026-08-20, before two container rebuilds** — so **the preserved
two-section base is not merely intact on disk, it is re-derivable from raw data in the rebuilt
environment.** Every effect above is read against these two floors and against the seed-to-seed
floor (*r* ≈ 0.996, Jaccard ≈ 0.76). P28.

**Frozen consequence for H1.** The DCA environment carries forward, so §8's free experiment runs at
the published default on **both** arms — which removes `denoise=False` from the list of things that
could explain a mouse-only artefact. Budget CPU time or a CUDA-11 image: TensorFlow 2.4 could not
use this box's GPU (wants CUDA 11; absent) and fell back to CPU.

**D3 re-anchoring, frozen.** The primary alternative anchor is an **eight-gene proliferation set**
(`Kif20a Ncaph Anln Ect2 Gtse1 Uhrf1 Fen1 Clspn`), chosen because every member is on the mouse
panel and absent from **every** `A_*.txt` and `B_*.txt`; the script asserts this at run time.
**`Lmnb1`, which §7 of the planning document proposes, is itself in `B_downstream_arrest` and
`B_secondary_senescence` and must not be the primary anchor** — it is reported as a secondary. All
anchor decisions are taken on the Spearman correlation **after linearly removing the rank of log
transcript counts from both variables**; on the raw correlation every anchor on this panel is a
detection-rate readout. Source: `reports/CS_PHASE8_CALLERS.md` §5,
`results/phase3/deepscence_anchor_decisions.csv`.

**Sign-invariant summary, pre-registered as the form to carry to H1**: rank by `|score|` as well as
by score. On M1 the Tier A–DeepScence agreement is `1.285×` chance under `|score|`, above chance in
10 of 11 sections — as high as either signed version, so it cannot be an artefact of the anchor.
H1's anchor problem will recur and its polarity may not be recoverable at all.

### 3.10 Gene sets — Tier A–E membership

**Frozen in `reports/PREREG_PHASE8_genesets.md`, which is the gene-set section of this
pre-registration.** Summary of what that file fixes:

- **Tier A PRIMARY: `A_SENDER_FINAL_strict`, n = 33 on both arms**, disjoint from all seven
  response modules (PI decision D1).
- **Pre-registered sensitivity: the seven `A_sender_for_<module>` sets**, each disjoint from the
  one module it is paired with. Gated in the same run; both gates must pass.
- **Tier B: seven modules.** M1 126 / 68 / 100 / 190 / 125 / **31** / 108; H1 120 / 71 / 126 / 231 /
  113 / 36 / 116.
- **Tier C** ligand–receptor, **Tier D** nuisance (13 covariate names, arm-specific anatomical
  covariate: liver zonation for M1, red pulp / white pulp for H1), **Tier E** controls.
- **Frozen directories:** `genesets/` (mouse, now carrying the promoted C6 sets) and
  `genesets/human/`. `variants/` is reported but never used.
  `genesets/human/FROZEN_MANIFEST.csv` — 35 frozen files, 8 variants, SHA-256 per file.
- **§11 gate: PASS on both arms.** Re-verified for this draft; see §4.

### 3.11 Section admissibility (M1)

Frozen from `results/phase3/summary_phase3.txt` (pinned, md5 `ecf86b9ca5460f31290e2f4c9e822ea2`),
§8 Test 3 on `Cdkn1a`⁺ hepatocyte prevalence:

***[Corrected 2026-08-27 — see §0.0 item C-11. The quoted digest is 31 hex characters and so
cannot be an MD5; the pinned file was also overwritten at 2026-08-27 09:06. **The pin is
md5 `dc92ddc6605eef52f6359aeab4e16fd7`.** The three lists below are unaffected and re-verify
exactly against `sasp_phase3.py:67-72`.]***

- **In band, PRIMARY, 6 animals:** 7259, 7260, 7001, 7248, 7352, 7435.
- Over the 20 % ceiling, excluded: 7239, 7448, 7361, 7450.
- Below the 1 % floor, excluded: 7250.

**H1 has no analogue of this rule** and none is invented: all 7 spleens are normal, there is no
disease contrast, and the admissibility rule is `Cdkn1a`⁺-prevalence-based. **Pre-registered: all
7 H1 sections are analysed; the A3 prevalence band (1–20 %, ≥ 200 senders, ≥ 5,000 non-senders) is
reported per section and per caller, and is a *reported* quantity, not an exclusion rule, on H1.**
For the SenePy caller, A3 is evaluated on the subset of cells that receive a SenePy score at all,
and the two denominators are stated separately (§3.7).

---

## 4. The §11 disjointness gate is now enforced automatically

`B_oxidative_stress` clears the ≥ 30 floor by exactly one gene on the mouse arm (**31**), and that
gene is **`Junb`**, from GSE310392's 100-gene custom add-on. Under the stock-panel-only definition
(5,006 rows) mouse B6 is **30** — the floor, margin 0.

**New in this draft: `code/gate_genesets_guard.py`** runs the gate on **both arms** and exits
non-zero if either fails, and it is wired to fire automatically:

- **`.githooks/pre-commit`** (installed via `git config core.hooksPath .githooks`) blocks any commit
  that stages a file under `genesets/` or either mouse panel file until the gate passes.
- **`.claude/settings.json`** declares a `PostToolUse` hook on `Write|Edit|NotebookEdit` running
  `code/hook_geneset_gate.sh`, so an agent editing a gene set fires the gate immediately.
- A SHA-256 manifest of 96 watched gene-set and panel files
  (`genesets/.geneset_manifest.json`) names *which* file moved.

The human half **delegates to `code/gate_disjointness_human.py` unchanged**, as a subprocess — that
script is the gate of record and already exits non-zero. The mouse half re-applies the same §11
assertions read-only to the **promoted `genesets/*.txt`** that the pipeline actually reads; the
existing mouse gate lives inside `build_genesets_mouse_c6.py`, which gates the sets it is *building
in memory* and writes files as a side effect, so it is the wrong instrument for a post-change check.

**Current verdict: PASS on both arms.** Mouse: A = 33; B = 126 / 68 / 100 / 190 / 125 / **31** /
108; A ∩ B_k = 0 for all seven; all seven per-module sensitivity sets ≥ 15 and disjoint from their
own module. Human: delegated gate exits 0, frozen configuration PASS.

**Correction to `reports/PREREG_PHASE8_genesets.md` §3.1 and deviation D10.** That file states that
"any later change that removes even one gene from mouse B6 fails the §11 gate". **Measured, it does
not.** On the authoritative 5,097-gene panel B6 has n = 31 against a floor of 30, so removing one
gene leaves exactly 30 and the gate still **passes**; removing two fails. The one-gene statement is
true only under the CSV-only 5,006-gene definition, where B6 is already 30. Verified by
falsification against the guard: dropping 1 gene from mouse B6 → gate PASS; dropping 2 → gate FAIL;
injecting `Junb` into `A_SENDER_FINAL_strict` → gate FAIL. **The correct pre-registered statement
is: mouse B6 has a margin of exactly one gene, it is the only module anywhere near the floor, it is
the module carrying a documented source substitution, and it should be read as the weakest of the
seven independently of whether it passes.**

---

## 5. Named primary outcome

**PRIMARY OUTCOME — the reach bound.**

> The **surviving fraction of the naive distance-to-nearest-sender amplitude under the combined
> N2 + N5 + N6 nuisance design**, at the primary sender call `tierA_p95` on the strict-33 Tier A,
> taken as the **median over the arm's reportable fits** (positive naive amplitude, spatial
> block-bootstrap CI excluding zero, ≥ 2,000 receivers), reported with its **inter-quartile range
> across those fits** *[wording corrected 2026-08-27 — see §0.0 item C-4; this read
> "paired-bootstrap interquartile range", and there is no bootstrap in that bracket]*; **together with the controlled amplitude `|β| / sd(y)` in response-SD units,
> compared against that arm's own 80 %-power detectable bound.**

One number per arm, one sender definition, one design. The per-module Tier A sets are the
pre-registered sensitivity and are reported alongside (decision D1), not instead.

***[Post-freeze decision 2026-08-27 — see §0.1 item **D-B**. On **H1** the sender definition this outcome is computed at is **`tierA_merged_p95`**, the same percentile rule at the label family the estimator stratifies on; the frozen `tierA_p95` fine-label call is reported beside it as the frozen-literal sensitivity. **The estimand, the design and the M1 numbers are unchanged.** §0.1 open item 1 records that this contradicts the literal wording of §5 and of R1/R2 below, and that the PI has not stated which of the two calls R1 is scored on.]***

***[Resolved 2026-08-27 — see §0.2 item **F1**. **The primary outcome, and R1 and R2 with it, are scored on the PRIMARY call of the arm being scored** — `tierA_p95` on M1, `tierA_merged_p95` on H1 — with the frozen-literal fine-label call computed and reported on H1 as the declared sensitivity. **This resolves an ambiguity; it does not relax the outcome definition**: no threshold, bracket or direction moves, and M1's numbers are untouched.]***

**Why this and not a length constant.** λ̂ rails at a grid bound in a majority of M1 fits, so a
fitted length constant is not the estimand; the amplitude and its survival under conditioning are.
This is the operational form of the question the PI has kept — *how far does senescence signalling
reach* — answered as a bound.

**M1 benchmark values, PROVISIONAL** (pre-C6; from the pinned
`results/phase3/summary_phase3.txt`, primary call, 160 reportable fits; task 8.7 recomputes all of
them under the promoted C6 sets):

***[Corrected 2026-08-27 — see §0.0 items C-11 and C-5. The pin is md5
`dc92ddc6605eef52f6359aeab4e16fd7`, not the 31-character string in §3.11, and the primary call
has **153** reportable fits, not 160. The frozen replacements for every value in the table
below are in C-5.]***

| Quantity | M1 value | Source |
|---|---|---|
| Naive amplitude, median `\|β\|/sd(y)` | **0.326** | `summary_phase3.txt` §6, `tierA_p95` row |
| SF, N2 + N5 + N6 | **0.082 [−0.099, 0.249]** | `summary_phase3.txt` §1 |
| SF, N5 alone | 0.084 | " |
| SF, N1 | 0.716 | " |
| SF, N2 matched decoy | 0.943 | " |
| Controlled amplitude | 0.027 response-sd | `reports/CS_PHASE3.md` §, derived |
| 80 %-power detectable bound | **0.203 response-sd** (SE 0.073) | `reports/CS_PHASE3.md`, `reports/BIO_DELIVERABLE6_DISCUSSION.md` |
| λ̂ railed at a grid bound | 63 % (200/315) | `summary_phase3.txt` §1 |
| Corrected N3 / N4 (tile) | **0.974 / 0.962** | `results/phase3/sf_summary_c1.csv` |

---

## 6. Replication criterion — written before any H1 outcome is computed

H1 is scored against these four criteria, fixed now. **All thresholds are absolute, not
"whatever M1 turns out to be."**

**R1 — primary (the bound replicates).** H1 replicates M1 iff the H1 median SF under
N2 + N5 + N6, at `tierA_p95`, over H1's reportable fits, has a **paired-bootstrap interval that
includes 0** *and* an **upper limit below 0.50**. (M1: 0.082 [−0.099, 0.249] — includes 0, upper
limit 0.249.)

***[Corrected 2026-08-27 — see §0.0 item C-4. "Paired-bootstrap interval" is a misnomer: the
bracket is the **inter-quartile range across the reportable fits** (`sf_summary.csv` `q25/q75`,
`np.quantile`), and the pre-registered bootstrap emits no interval on the median across fits.
**Read R1 as: the IQR across H1's reportable fits includes 0 and its upper quartile is below
0.50.** The criterion, the threshold and M1's own outcome are unchanged; the frozen M1 values
are 0.088, IQR [−0.017, 0.234].]***

***[Post-freeze decision 2026-08-27 — see §0.1 item **D-B**. R1 and R2 name the call `tierA_p95`. On **H1** the primary sender call is **`tierA_merged_p95`** — the identical percentile rule at the label family the estimator stratifies on — with the frozen fine-label call reported beside it. **The thresholds (IQR includes 0; upper quartile < 0.50; amplitude below H1's own 80 %-power bound) and M1's outcome do not change**; only the H1 fit population does. §0.1 open item 1 records that the PI has not stated whether R1 is scored on the merged call alone or on both.]***

***[Resolved 2026-08-27 — see §0.2 item **F1**. **R1 and R2 are scored on the PRIMARY sender call of the arm being scored: `tierA_p95` on M1, `tierA_merged_p95` on H1** (§0.1 D-B), with the frozen-literal fine-label call computed and reported on H1 as the **declared sensitivity**, every time. **This resolves the ambiguity §0.1 open item 1 recorded; it does not relax R1 or R2.** No threshold moves — R1 is still "IQR includes 0 **and** upper quartile < 0.50" and R2 is still "median controlled `\|β\|/sd(y)` below **that arm's own** 80 %-power bound"; M1's outcome is unchanged (0.088, IQR [−0.017, 0.234]; 0.029 against 0.183 — §0.0 C-5); and the merged call **adds** senders and fits rather than pruning them (T/NK: 23 → 1,958 SPLN14, **0 → 1,241** SPLN43, 295 → 1,675 SPLN44; `results/phase9_h1/a3_prevalence_by_type.csv`), which if anything makes an IQR that **excludes** 0 — i.e. non-replication — easier to observe, not harder. **Which way H1 actually falls is a Phase-10 outcome and is not stated here.**]***

**R2 — amplitude (nothing above the bound).** H1's median controlled amplitude `|β|/sd(y)` under
N2 + N5 + N6 is **below H1's own 80 %-power detectable bound**, computed from H1's own standard
error by the same rule that produced M1's 0.203. R2 is evaluated with H1's bound, not M1's, because
the arms differ in n, prevalence and depth.

**R3 — geometric predictions (must replicate in any tissue).** These are properties of point
patterns, not of biology, and a failure means a pipeline is broken:
 (a) **Poisson identity** — median distance-to-nearest-sender against sender density fits the
 homogeneous-Poisson slope of −1/2 with **r² ≥ 0.95** (M1: 0.984);
 (b) **λ̂ railing rate** is reported for both arms; no threshold is set, because it is diagnostic;
 (c) **null destructiveness** — every in-tissue N3/N4 variant retains **≥ 95 %** of shifted senders
 with a real neighbour inside the 100 µm window (M1: 96.6–100 %), and the published bounding-box
 variants do not (M1: 77 % / 92 %).

**R4 — non-replication (outcome B).** H1 fails to replicate iff the H1 SF interval under
N2 + N5 + N6 **excludes 0** *and* the controlled amplitude exceeds H1's own detectable bound in
**≥ 2 of the 7** Tier B modules. On R4 the mouse arm is **not** dropped: both arms are reported and
the candidate explanations are named in advance (species, tissue — confounded by design and not
separable — prevalence, panel, sender-caller behaviour, architecture, donor variance), together
with the follow-up that would distinguish them (a human liver arm at adequate panel depth, which
does not currently exist publicly).

**R5 — hard prerequisites, checked before R1–R4 are computed.** A2 (disjointness on the real H1
panel — pre-verified: passes), A5 (matched-decoy contrast, |SMD| ≤ 0.1), and **A7 (negative-control
kernel flat under the N5-conditioned design, pooled over sections)**. A7 on H1 is **only adequately
powered pooled across sections** — one spleen section cannot deliver it; seven can (M1's per-fit
resolution is ±0.137 SD, 1.8× larger than the conditioned biological amplitude, while the
section-clustered pooled CI half-width is ±0.018 SD). Failure of A2 or A5 → **outcome C**.

**Bootstrap scope on H1.** Section-level bootstrap only, labelled as such. Donor count is 7 and the
donor bootstrap is not run.

---

## 7. §18 outcome table — decided before looking

| | Outcome | What it means, and what we do |
|---|---|---|
| **A** | **Both arms agree, no kernel above bound.** | The strongest version. Two species, two tissues, two labs: the distance-dependent SASP response is not separable from technical and compositional confounding at achievable power. Reported with the corrected null battery, A7, and the measured 9–16 % false-positive rate as the supporting evidence for *why* the bound is the honest answer. |
| **B** | **H1 shows a surviving kernel** (criterion R4). | Report both arms; do not drop mouse. Hypothesise: species, **tissue** (confounded by design, not separable), prevalence, panel, sender-caller behaviour, architecture, donor variance. Name the distinguishing follow-up: a human liver arm at adequate panel depth. |
| **C** | **H1 fails A2 or A5.** | A data-availability finding, one honest paragraph, and move to the next candidate in the §12.1 screen (runners-up: GSE336890 kidney, GSE335963 bone marrow). |
| **D** | **The confound structure differs but the null result does not.** | The confound is context-specific, so no published characterisation from a different dataset can be trusted. Strengthens the case for running the battery every time. |
| **E** | **C1 changes M1's N3 result.** | Handled in the correction ledger. **Already resolved: it does not.** N3-tile 0.974 and N4-tile 0.962 against the published 1.000 and 0.964, invariant across six sender calls spanning 0.5–9.0 % prevalence (`results/phase3/sf_summary_c1_n7.csv`). *Those values are **PROVISIONAL** — pre-C6 — but the invariance across a 0.5–9.0 % prevalence range and three unrelated sender callers is a property of the section outline and the sender point pattern, so the conclusion transfers even though the numbers will move.* Contribution 3 survives; what changes is that N3-occ is declared degenerate and N3-swap is declared a label permutation, not a corrected torus shift. |
| **F** | **DeepScence's instability appears in H1 too.** | §8's experiment. A real limitation of a widely used tool, verifiable by anyone. See §8 for the prediction and its falsifier. ***[2026-08-27, §0.1 **D-A** and **D-C**: outcome F is realised, and in a **different form** from the one anticipated. The instability that appears natively is **seed** instability in the frozen PRIMARY `denoise=False` configuration (*r* = 0.372, top-5 % Jaccard 0.211 at full section size), not the anchor instability P-ii predicted — **P-ii is falsified**. F is reported jointly with D-A's estimator change.]*** |

Outcomes are not mutually exclusive: A or B is the primary axis, D, E and F are orthogonal and any
combination may be reported.

---

## 8. The §8 DeepScence prediction — stated before H1 runs

**CoreScence is a human gene set.** On M1 it runs ortholog-remapped (4,845 of 5,097 panel genes);
**on H1 it runs natively, exactly as published, with no remapping.** That separates two explanations
currently confounded: is DeepScence's observed behaviour a property of *our mouse adaptation*
(remapping loss, `denoise=False`), or a property of *the published tool*?

**D2 has already removed one half of that confound.** The DCA environment now exists
(`code/setup_dca_env.sh`, §12), so `denoise=True` is runnable on **both** arms — `denoise=False`
can no longer be offered as the explanation for a mouse-only artefact, and the H1 comparison is
run at both denoise states. It also supplies a sixth prediction, below, from a direction M1 has
already measured.

***[Post-freeze result 2026-08-27 — see §0.1 item **D-C**. **P-ii below is FALSIFIED and it leads the reporting of this section.** The published `CDKN1A` anchor is stable **20/20 folds in all seven H1 sections** (`results/phase9_h1/deepscence_anchor_h1.csv`, `stab_cdkn1a` = 1.00 ×7), so **M1's polarity flip (P12, P13) is partly an artefact of our own ortholog remapping, not a defect in the published tool** — which is exactly the separation this section was designed to make. **What that does not touch**, all measured on the **native** run and so unaffected by remapping: the **88 %** CoreScence circularity (29/33, `results/phase7_jobA/gate_result_human.json`), the **seed instability** (§0.1 D-A), and DeepScence × `CDKN1A`⁺ at **6.436 pooled** natively (`results/phase9_h1/caller_agreement_pooled.csv`). Every verdict below is at `random_state = 0` and must be re-read against D-A.]***

***[Extended 2026-08-27 — see §0.2 item **F3**. **Every H1 DeepScence MAGNITUDE below and in `CS_PHASE9_H1_AUDIT.md` §10 is PROVISIONAL pending the five-seed consensus** — they all descend from the single seed-0 score files `data/processed_h1/deepscence_h1_<section>.csv` (`code/h1_deepscence.py:76`, `random_state=0`). That covers P-i's depth loadings (+0.1822 … +0.3540), P-v's within-type Q5/Q1 (0.244 … 1.169), every anchor ρ partial in `deepscence_anchor_h1.csv` (including `stab_cdkn1a` itself), and every pooled or matched agreement ratio with a DeepScence term (6.436, 2.069, 1.602, 1.102, 1.093, 0.998, 0.290, and the per-section ranges 3.459–10.115 and 0.972–2.816). **Three things are explicitly NOT provisional**, because the flag is broader than the evidence: **P-iv's 29/33 = 0.8788 circularity** is a gene-list membership computation in which no score is produced and no `random_state` enters (`results/phase7_jobA/gate_result_human.json`); Tier A × SenePy **0.874** carries no DeepScence term; and the seed-check statistics themselves (*r* 0.3719, Jaccard 0.2107) are statements *about* the seed pair. **D-C's verdict stands as a direction**; §0.2 F3 records that `stab_cdkn1a` = 1.00 is itself a seed-0 magnitude and that the argument for retaining the verdict is asserted, not measured.]***

**PREDICTION, registered before any H1 expression value is read: it is a property of the tool, and
it will recur natively in H1.** Concretely, and each falsifiable on its own:

| # | Prediction on H1 (7 sections) | Falsified if |
|---|---|---|
| **P-i** | DeepScence's score correlates **positively** with per-cell transcript counts (Spearman ρ) in **≥ 5 of 7** sections, magnitude comparable to M1's +0.29 … +0.56 | ρ is null or negative in ≥ 3 of 7 |
| **P-ii** | The published `CDKN1A` anchor is **weak, unstable or inverted** — depth-partialled fold-split sign stability < 0.90 — in **≥ 1 of 7** sections | all 7 sections have stability ≥ 0.90 |
| **P-iii** | The **sign-invariant** (`\|score\|`) depth- and type-matched agreement with the Tier A caller is **above chance, pooled** (ratio > 1.10) | pooled ratio ≤ 1.05, or below chance ***[MARGINAL, not scored yet — 2026-08-27, see the F4 marker under this table and §0.2]*** |
| **P-iv** | CoreScence circularity with the frozen Tier B modules is **≈ 88 %** natively | materially below 88 % on the native panel |
| **P-v** | Within cell type, DeepScence is **bottom-of-the-depth-distribution selecting** (Q5/Q1 < 1) in ≥ 5 of 7 sections, as it is in all 11 M1 sections | Q5/Q1 > 1 in ≥ 3 of 7 |
| **P-vi** *(new, from D2)* | Switching `denoise=True` on H1 **raises** the score's depth loading rather than lowering it, as it does on 3 of 3 M1 sections (Δρ +0.13 … +0.25) — i.e. §4 (D-b)'s premise fails natively too | Δρ ≤ 0 in ≥ 5 of 7 H1 sections, which would make the M1 result an artefact of ortholog remapping |
| **P-vii** *(new, from D2)* | The `denoise=True` seed instability recurs: **≥ 1** of 3 seeds on an H1 section gives a top-5 % sender set with Jaccard < 0.30 against the others | all three seeds agree at Jaccard ≥ 0.60, as `denoise=False` does on M1 |

**The evidence already leans one way, and this is recorded so the prediction is not read as
neutral.** CoreScence is **88 % circular with the frozen Tier B modules on the human panel**
(29 of 33 on-panel CoreScence genes are in ≥ 1 Tier B module —
`results/phase7_jobA/gate_result_human.json`), against 76 % under the superseded curated B7 and
**79 %** in the pre-C6 mouse arm (26 of the 33 CoreScence genes reachable on the mouse panel;
`results/phase7_jobA/corescence_circularity_mouse.json`). Under the promoted C6 sets the mouse arm
is **also 29/33 = 88 %**. **That 88 % is measured natively and cannot be attributed to ortholog
remapping.** B7 alone accounts for 18 of the 33 on both arms. *(Corrected 2026-08-27 from
"24/35 = 69 %", a typed-in literal reproducible under no mapping convention —
`reports/AUDIT_PHASE8_FACTCHECK.md` M1.)*

**Either answer publishes.** If P-i…P-v hold, we have found a real limitation of a widely adopted
tool, in the species it was written for. If they fail, our mouse adaptation caused the M1 behaviour,
which strengthens everything else by removing a confound — and we say so.

***[2026-08-27 — see §0.1 item **D-C**. **Both halves of this sentence came true, on different predictions**: P-i, P-iv and P-v hold natively (a real limitation of the tool) and **P-ii is falsified** (our mouse adaptation caused the M1 polarity behaviour). **We say so, at the front of §8 and not in a table.** P-iii is confirmed only at 1.102 against its own 1.10 threshold and is seed-fragile — §0.1 open item 5.]***

***[Corrected 2026-08-27 — see §0.2 item **F4**. **P-iii is recorded as MARGINAL, and its outcome is not stated in either direction until it has been re-scored on the five-seed consensus.** The sentence above calls it "confirmed only at 1.102"; **on this estimator that is not a confirmation.** The registered rule confirms at > 1.10 and falsifies at ≤ 1.05, so **(1.05, 1.10] is pre-registered as neither** — and H1's 1.102 (z = 5.67, above chance 5 of 7; `results/phase9_h1/caller_agreement_pooled.csv`, row `tierA_score` × `abs_deepscence_score`) sits **0.002 above the top of that indeterminate band**, on an estimator whose full-section seed-to-seed top-5 % Jaccard is **0.2107** (`results/phase9_h1/d2_stability.csv`). **No threshold is moved and the falsifier is not weakened; P-iii is simply not scored yet**, and may return CONFIRMED, MARGINAL or indeterminate on the consensus. P-i, P-ii, P-iv and P-v are not affected by this marker.]***

**Consequence already frozen regardless of the outcome:** the **strip-and-refit sensitivity** (fit
each Tier B module with the CoreScence-shared genes removed, report both amplitudes) is part of the
frozen run order, not an afterthought (PI decision D6).

---

## 9. Analysis and run order after the tag

1. **Phase 9** — H1 audit A1, A2 (gate), A3, A4, A5 (gate), A6 (build and validate the red/white
   pulp covariate), A7 (must be flat, pooled), A8 (cross-arm on the 2,425-gene ortholog-intersected
   panel). **Stop on A2 or A5 failure → outcome C.**
2. **Job B** — cell types with the 22-label spleen marker set; cross-check against the depositors'
   four-level `annotations.csv.gz` (which covers fewer cells than the matrix, so its QC filter is
   characterised first); then the sender callers.
3. **Phase 10** — H1 through the frozen pipeline: naive, N1–N8, controlled fits, kernel families,
   superposition vs nearest, proximal vs downstream. Then the §8 comparison, the §17 two-arm table,
   Figures 5 and 6, and the revised 2, 3, 4.
4. **Every cross-arm number is reported twice**, on the 2,425-gene ortholog-intersected panel and
   on each arm's full panel (test A8).

   ***[Post-freeze decisions 2026-08-27 — see §0.1 items **D-A** and **D-B**. Two **estimator** asymmetries now sit across the arms in addition to the panel axis: DeepScence is a single `random_state = 0` run on M1 and a five-seed consensus on H1 (D-A), and the Tier A sender call is fine-label on M1 and merged-label on H1 (D-B). **Both are declared, both have the alternative computed, and neither is a panel effect** — so a cross-arm DeepScence or sender-call number carries the estimator asymmetry as well as the two-panel requirement.]***

   ***[Registered together 2026-08-27 — see §0.2 item **F5**, which is the single place all **three** cross-arm asymmetries and their costs are recorded, so that no reader meets them one at a time. **(1) Panel** (this item): M1 runs DeepScence ortholog-remapped onto 4,845 of 5,097 genes, H1 natively on 5,093 — measurable and large, since cutting H1 to the 2,425-gene intersection leaves **26 of the 33** Tier A genes and moves the `tierA_p95` sender set at Jaccard **0.5222 – 0.5747** (`results/phase9_h1/a8_ortho_sender_shift.csv`); mitigated by the two-panel rule above. **(2) DeepScence estimator** (D-A): single seed on M1, five-seed consensus on H1 — **any cross-arm DeepScence difference smaller than H1's own between-seed spread is uninterpretable** (*r* 0.3719 / Jaccard 0.2107, `results/phase9_h1/d2_stability.csv`), and **there is no mitigation on the M1 side**, which has nothing to average (determinism floor 0.99999913 / 0.99999995). **(3) Tier A call level** (D-B): fine on M1, merged primary on H1 — the sender populations are defined at different label granularity (the merged call restores 0–13.6 % of cells per section, and SPLN43's merged T/NK goes from **0 senders in 24,815 cells** to 1,241); **partially** mitigated because the fine call is computed on H1 too, **not** mitigated on the M1 side, where no merged counterpart is computed (§0.1 open item 3, still open). **A cross-arm DeepScence number carries axes 1 and 2; a cross-arm sender-call number carries axes 1 and 3; no cross-arm number in this project is affected by fewer than two of the three**, and each must be reported with the axes it carries named.]***

---

## 10. What may never be reported, whatever the data says

These are consequences of measurements already made, and they bind both arms.

1. **No naive distance kernel, and no N2-only kernel, may be reported as a distance effect.**
   ***[Corrected 2026-08-27 — see §0.0 item C-1. The digits in this paragraph are the pre-C6
   05:19 A7 file; paragraph 2 below is the frozen 09:06 file. Frozen values: naive −0.0744
   [−0.1306, −0.0182] p = 0.0145; N2 −0.0642 [−0.1113, −0.0172] p = 0.0124, 86 % undiminished;
   N5 +0.0038 [−0.0186, +0.0261] p = 0.715; conditioned biological amplitude +0.0310.]***
   A7 measured the mouse arm's negative-control response against distance to nearest sender:
   naive clustered mean **−0.070 SD [−0.128, −0.012], p = 0.023**, and the matched-decoy contrast
   leaves **−0.061 [−0.111, −0.012], p = 0.020** — 80 % undiminished. The N5 technical covariate
   block removes it completely (**+0.007 [−0.011, +0.025], p = 0.41**), and that residual excludes
   the conditioned biological amplitude of +0.036 SD. Source: `results/phase3/a7_summary.csv`,
   `all_controls` rows. **The raw assay is not flat; a matched-decoy contrast is not a substitute
   for the technical covariate block on this assay.**
   **Name the response correctly.** `all_controls` is the **pooled** control feature set — 40
   negative-control probes + 609 negative-control codewords + 21 genomic controls, and the
   codewords carry ~73 % of the counts. Per family under `base`: codewords −0.0604 [−0.1085,
   −0.0123] p = 0.0188; genomic controls −0.0307 [−0.0558, −0.0056] p = 0.0213; and the **40
   negative-control probes on their own are flat, −0.0225 [−0.0527, +0.0078], p = 0.129.** Since
   `PREREG_PHASE8_genesets.md` §11 designates `E_negative_control_probes` the *primary technical
   null* and Phase 9 item 9.4 repeats it, **the pre-registered primary A7 response passes on M1**;
   the non-flatness is carried by the codewords and genomic controls. Report it as "pooled
   negative-control features", never as "negative-control probes"
   (`reports/AUDIT_PHASE8_FACTCHECK.md` R1).
2. **No caller-independence claim.** See §13 P1.
3. **No age-stratified or young-vs-old claim on H1.** Age is a continuous covariate only (PI
   decision D4). Two donors are over 55.
4. **No marginal-zone-specific confirmatory claim.** Exploratory only (PI decision D3b, §13 P7).
   ***[Post-freeze decision 2026-08-27 — see §0.1 item **D-D**. **Strengthened, not relaxed: P7 is retired as UNANSWERABLE.** The `Marginal zone B cells` label is **never realised in any of the 7 H1 sections**, so the exploratory MZ claims have no support to be exploratory about. **Never testable on this panel — which is not the same as tested and null**, and is recorded as the different statement it is. This prohibition stands a fortiori.]***
5. **No cross-arm difference attributed to species or tissue.** They are confounded by design —
   mouse liver against human spleen. This belongs in the abstract, not the limitations.
6. **`CXCL8`/`CXCR1` results must never be reported as replicating a mouse result** — no mouse
   ortholog exists. Conversely `MMP3` and `TIMP1` are mouse-only. `CXCL2`/`CXCL5` are on **both**
   panels and are a **map gap**, not a biological asymmetry.
7. **The circularity figure "1.51–2.85×" must not be quoted.** Measured over 11 sections the
   DeepScence-vs-`Cdkn1a`⁺ pair is 0.963–2.849 with a **median of 1.071** and a pooled 1.255; both
   published values were the two largest of the eleven (`results/phase3/caller_coverage_gate.csv`).
8. **The composition-matched number may never be reported alone.** The matched-decoy SF of 0.9837
   ("composition removes 1.6 %") and the covariate-adjusted SF of 0.1461 ("composition removes
   85.4 %") come from the **same variables on the same fits**. Quoting the first without the second
   states the opposite of what the data says. Wherever the matched number appears, `type_adj`
   (65.9 %) and `typecomp_adj` (85.4 %) appear beside it. P23, P25.
9. **`mor` may never be reported as evidence that normalisation cannot move this caller.** §6's
   named estimator removes only 11.4–23.6 % of log-depth variance on this panel; its null is a null
   about a weak normalisation. Wherever it appears, the `lib` result (100 % removed; depth loading
   −74 % and −93 %; 46 % and 100 % of sender calls changed) appears beside it. P27.
10. **No `denoise=True` number may be reported from a single seed without its seed-stability
   companion.** One of three seeds gave a top-5 % sender set disjoint from the other two
   (Jaccard 0.000). Any published-default sensitivity result carries that fact. P26.
   ***[Extended 2026-08-27 — see §0.1 item **D-A**. **On H1 this applies to `denoise=False` as well**, i.e. to the frozen PRIMARY configuration: at full section size the seed-0 and seed-1 scores agree at *r* = 0.3719 and top-5 % Jaccard 0.2107 (`results/phase9_h1/d2_stability.csv`). **No single-seed H1 DeepScence number may be reported at all** once the five-seed consensus exists; until it does, every such number is labelled `random_state = 0` and carries this spread. M1 is unaffected.]***
    ***[Extended 2026-08-27 — see §0.2 item **F3**. Until the consensus exists, every H1 DeepScence **magnitude** already published in `CS_PHASE9_H1_AUDIT.md` §10 is **PROVISIONAL** and must be labelled so wherever it is quoted — the depth loadings, the within-type Q5/Q1 values, every anchor ρ partial (and `stab_cdkn1a` itself), and every agreement ratio with a DeepScence term. **P-iv's 29/33 = 0.8788 circularity is exempt**: it is a gene-list membership fact with no score and no seed in it.]***
11. **`rho_signed_dz_vs_depth` for the D2 `raw` control rows (−0.47, −0.16) must not be quoted** —
   it is the direction of numerical noise on a shift of 0.0002–0.001 z-units. P28.
12. **Four sentences from the two-section base must be struck** (§13 P1): "0.93–1.22× of chance …
   i.e. they are statistically independent"; "Four of six pairs sit at 0.93–1.22×"; "the one pair
   that looked concordant in sham is anti-concordant in SBR"; "DeepScence's correlation with
   sequencing depth reverses sign between two sections of the same study".

---

## 11. Figure policy, frozen

Committed figures are held at their committed state and regenerated **exactly once**, from the
frozen configuration, at task 8.7. Regenerated candidates live in `figures/revised_candidates/`
with `_REVISED` suffixes and a README carrying the regeneration ledger.
`code/check_figures_guard.py` is a content-hash guard over the 27 committed figures (PDF date
stamps stripped before hashing). **Two actors have already collided in `figures/` and nothing
warned**; the guard makes the policy enforceable rather than advisory.
**`code/make_figure2.py` must never be run** — it is a superseded second producer of `figure2a`;
the live producer is `make_phase5_figs.py --which 2a`. It now refuses to run.

Every figure writes a `*_data.csv` beside it so every plotted number is auditable, and both `.png`
and `.pdf` are emitted; **PNGs are the reproducibility comparison** because matplotlib date-stamps
PDFs.

---

## 12. Environment

Pinned in `requirements.txt`. `kneed==0.8.6` and `openpyxl==3.1.5` were added on 2026-08-27:
`kneed` is a **hard** import of DeepScence (`DeepScence/io.py`), so the pinned environment could
not reconstruct the sender caller that produces `deepscence_score` — the reproducibility claim was
broken until then. A persistent interpreter lives at `/workspace/envs/sasp311/bin/python` on the
network volume, which survives container resets; the environment has been lost twice.

**Second, isolated environment for DCA.** `denoise=True` cannot run in the pinned stack (DCA 0.3.4
needs TensorFlow 2.4 / Keras 2.4, which will not sit alongside the 3.11 scientific stack). It runs
in a separate **CPython 3.8.19** venv built by **`code/setup_dca_env.sh`** (≈ 4 minutes), out of
process, exchanging matrices as plain `.npz`/`.npy`. **The main environment has no TensorFlow and
was not modified.** Pins archived at `results/phase8_d2/dca_venv_pip_freeze.txt` and
`dca_venv_python.txt`. TensorFlow 2.4 could not use this box's GPU — it wants CUDA 11
(`libcudart.so.11.0`, `libcublas.so.11`, `libcudnn.so.8`), none present — so it runs on CPU:
8.6 min to denoise an 83k-cell section, 16.4 min for the whole `denoise=True` call against 5.2 min
for `denoise=False`.

**Memory — a reproducibility-relevant resource fact, recorded for Methods.** The container's real
ceiling is a **57.7 GiB cgroup** (`/sys/fs/cgroup/memory.max` = 61,999,996,928 bytes), **not** the
~251 GB `free` reports for the host; `free` is misleading inside the container and
`memory.current` must be read against `memory.max`. **DeepScence holds five dense
`n_cells × 4,845` float32 arrays at once**, so an 83k-cell section was OOM-killed with 11 GB free
and **~16 GB of cgroup headroom is the working requirement** for a section of that size — more for
the 200k+ ones. Six D2 jobs were OOM-killed before this was established, and three D1 sections were
OOM-killed at five-way concurrency earlier in the day and recovered sequentially. **Anyone
reproducing this must budget headroom per section, not per box.**

**Declared:** three files in `genesets/msigdb_mouse_2026.1.Mm/` are HTML error pages, not JSON
(`FRIDMAN_SENESCENCE_UP`, `FRIDMAN_SENESCENCE_DN`, `WP_NRF2_PATHWAY`) — MSigDB has no mouse version
of those sets. **No tier uses any of the three**, so no published mouse number is affected. The
human archive is clean, 27 of 27 valid.

---

## 13. Deviation table — Phase 8 additions

`reports/PREREG_PHASE8_genesets.md` §12 carries the **gene-set** deviations, **D1–D17, seventeen
rows** as it now stands (the roadmap's "16-row" description predates the D17 panel-definition row).
Those are not restated here. The rows below are everything that has landed since, numbered `P*` so
they cannot be confused with the `D*` gene-set rows.

***[2026-08-27 — the four post-freeze PI decisions **D-A … D-D** of §0.1 are deviations too, and are deliberately **not** renumbered into this `P*` series: the `P*` rows are Phase-8 additions taken **before** any H1 expression value was read, and D-A–D-D were taken **after**. §0.1 is the register for them. The rows below that they touch — **P7** (retired by D-D), **P12** and **P13** (narrowed by D-C) — carry their own inline markers and their original text.]***

| # | Deviation | Reason / evidence |
|---|---|---|
| **P1** | **Caller independence is falsified. The motivating claim is restated, not defended.** | See the dedicated table immediately below — the number moved twice, for two different reasons, and the two must not be conflated. |
| **P2** | **The raw mouse assay is not flat, so no naive or N2-only kernel may be reported.** | A7, first direct measurement, on the **pooled** control features (`all_controls` = 40 probes + 609 codewords + 21 genomic controls; the 40 probes **alone** are flat at −0.018, p = 0.18, so the pre-registered primary A7 response passes and this finding rests on the codewords and genomic controls — §10.1, `AUDIT_PHASE8_FACTCHECK.md` R1): naive control amplitude **−0.070 SD [−0.128, −0.012], p = 0.023**, a quarter of the naive biological +0.314 SD in the same fits; **N2 leaves −0.061 [−0.111, −0.012], p = 0.020**; **N5 removes it: +0.007 [−0.011, +0.025], p = 0.41**, excluding the conditioned biological +0.077 SD. 825 control fits vs 1,155 module fits, same estimator. `results/phase3/a7_summary.csv`, `a7_verdict.txt`. §10.1. **Magnitudes are PROVISIONAL** — A7 was run at 05:19 on the pre-C6 sender calls; the three qualitative findings are what is frozen. A7 must be re-run after 8.7. ***[Corrected 2026-08-27 — see §0.0 item C-8. Every digit in this row is the pre-C6 05:19 A7 file. Frozen `results/phase3/a7_summary.csv` (09:06): naive `all_controls` **−0.0744 [−0.1306, −0.0182], p = 0.0145**; N2 **−0.0642 [−0.1113, −0.0172], p = 0.0124**; N5 **+0.0038 [−0.0186, +0.0261], p = 0.715**; the 40 probes alone **−0.0225 [−0.0527, +0.0078], p = 0.129** (still flat, so the primary A7 response still passes); biological **0.3120** naive / **0.0795** conditioned. **The closing instruction is discharged: A7 WAS re-run after 8.7 — `stat -c '%y' results/phase3/a7_summary.csv` → `2026-08-27 09:06:16` — so "Magnitudes are PROVISIONAL" and "A7 must be re-run after 8.7" no longer hold.** The three qualitative findings are unchanged.]*** |
| **P3** | **The estimator's false-positive rate is 9–16 % against a 5 % nominal.** | Same A7 run, full N6+N5 design, `frac_CI_excludes_zero` on responses with known-zero amplitude: 0.091 / 0.103 / 0.109 / 0.127 / 0.164 across the five control families (`results/phase3/a7_summary.csv`). **The five families are not five replications**: `all_controls` is the sum of the probe, codeword and genomic responses and `neg_probe_rate` is a ratio of two of them, so this is a range over correlated statistics. The **16 % upper end comes entirely from `neg_probe_rate`**, the one response whose denominator is an N5 column and which the caller report itself says is not a clean null; the clean-null subset is **9.1–12.7 %** (`AUDIT_PHASE8_FACTCHECK.md` M3). Quote the range as 9–13 % with 16 % as the `neg_probe_rate` outlier. Either a bootstrap FPR above nominal, or residual confounding N5 does not capture — both readings argue the same way. Stated in Methods. ***[Corrected 2026-08-27 — see §0.0 item C-9. The **0.127** is the pre-C6 `neg_control_probe` value; frozen is **0.145**, so the five families are **0.091 / 0.103 / 0.109 / 0.145 / 0.164** and the clean-null subset is **9.1–14.5 %**, not 9.1–12.7 %. **The instruction in this row is corrected to read: quote the range as 9–15 % on the four count-based responses, with 16 % as the `neg_probe_rate` outlier.** This row's header ("9–16 % against a 5 % nominal") is correct and unchanged; it was the body that carried the pre-C6 range. The pre-C6 list is also mis-ordered — pre-C6 `all_controls` is 0.091 and `neg_control_codeword` 0.109, the reverse of frozen.]*** |
| **P4** | **SenePy ships no spleen signature; it is not the same estimator across arms.** | `results/phase7_jobA/senepy_spleen_coverage.json`: 65 hubs, 10 tissues, `spleen_hub: false`. Of 22 spleen labels: **0 tissue-matched, 15 cross-tissue surrogate, 7 with no hub in any tissue.** M1 used tissue-matched mouse **Liver** hubs. Surrogates collapse: one blood memory-B hub scores all three B labels, one lung T-cell hub both T subsets. A3 must therefore be evaluated on the scored subset for SenePy and on all cells for the other callers, with both denominators stated. *Open recommendation, flagged not taken: demote SenePy to sensitivity on H1 and promote the `CDKN1A`⁺ call into the primary trio.* |
| **P5** | **Age is a continuous covariate only.** No young-vs-old contrast, no age-stratified prevalence claim. | PI decision D4. Ages 17/31/32/32/37/57/59 — five under 40, **two over 55**. A sparse continuum, not a two-group design. H1's value is a human replication of the geometry, not an ageing result. |
| **P6** | **Plasma cells are admitted by an explicit 3-marker exception to `MIN_MARKERS = 4`.** | PI decision D2b. Surviving on-panel markers `JCHAIN`, `MZB1`, `XBP1` — among the most specific markers in immunology, shared with no other label in the set; the label fails by exactly one gene, and only after the over-adjustment guard removes `MKI67` (which was CellMarker contamination from "dividing plasma cell" rows). **Recorded as an exception with its reason, so it does not read as the threshold having been 3 all along.** The alternative — folding plasma cells into the B compartment — would misattribute their signal: a plasma cell is not a B cell for receiver purposes. A 3-gene score is noisier and the label is reported as such. |
| **P7** | **Marginal zone B cells: label and compartment kept; every MZ-*specific* claim is exploratory.** | PI decision D3b. CellMarker 2.0 has only **6** spleen rows for this type, below the 8-row threshold, so its markers come from the weakest evidence tier (all-tissue fallback at ≥ 1 PMID). The compartment is structurally central to spleen and `D_spleen_marginal_zone` is one of the five A6 compartments, so dropping it would leave the A6 axis without its middle term. **No confirmatory marginal-zone hypothesis is pre-registered**; any MZ-specific claim is conditional on the label surviving the post-freeze re-gate against measured expression and against the depositors' own annotations. ***[Post-freeze decision 2026-08-27 — see §0.1 item **D-D**. **P7 IS RETIRED AS UNANSWERABLE.** The condition this row sets — "conditional on the label surviving the post-freeze re-gate" — was not met: `Marginal zone B cells` is **never realised in any of the 7 H1 sections** (verified over `data/processed_h1/celltypes_h1_*.csv`; it is one of six labels never realised, `CS_PHASE9_H1_AUDIT.md` §1.5). **The PI's grounds are that the question was never testable on this panel — not that it was tested and null.** **This row's original text is preserved deliberately**, because it records that the weakest evidence tier was flagged **in advance**: the 6 CellMarker spleen rows are verifiable in `genesets/human/markers_spleen_evidence.csv` (8 marker rows, each stamped `ALL-TISSUE FALLBACK >=1 PMID -- WEAKEST EVIDENCE (6 spleen rows)`) and in `results/phase7_jobA/build_markers_human_spleen.log` line 18. **Caveat on this row's own reasoning, flagged and not corrected**: the A6 axis is `score(D_spleen_red_pulp) − mean(score(follicle), score(tzone))` (`code/h1_a6_compartments.py:85-86`), so `D_spleen_marginal_zone` is scored as one of the five compartments but is **not** a term in the axis — dropping the label leaves the axis intact. §0.1 open item 6.]*** ***[Corrected 2026-08-27 — see §0.2 item **F6**. **The caveat above is promoted from a flag to a recorded correction, and it runs in D-D's favour: retiring P7 costs the A6 axis nothing.** Verified in the source rather than taken from a document — `code/h1_a6_compartments.py:57-58` scores five compartments including `D_spleen_marginal_zone`, `:85-86` builds the axis as `pulp = score(D_spleen_red_pulp) − 0.5·(score(follicle) + score(tzone))`, and `grep -n "marginal_zone" code/h1_a6_compartments.py` returns **line 58 only** — the label appears in the scored list and nowhere else in the producer. **So this row's stated reason for keeping the label ("dropping it would leave the A6 axis without its middle term") was already inaccurate before any H1 expression value was read**, and D-D's retirement is **strengthened**, not weakened. **This row's original wording is preserved deliberately**, as D-D requires: what is corrected is a justification, not a decision and not a number — no threshold, no gene set and no scored compartment changes, and `genesets/human/D_spleen_marginal_zone.txt` is not edited.]*** |
| **P8** | **RESOLVED (PI decision D16): the fitting window is pre-registered as `100 µm, fixed`. The 99th-percentile rule of §15 is recorded as provenance, not as a runtime computation.** | §3.2. The code has never computed a percentile at fit time (`run_phase3_nulls.py:59`, capped at `:169`); §15's wording describes how 100 was chosen once. Every published result already used 100 µm, so freezing the literal preserves continuity and re-runs nothing. Identical on both arms; `window.csv` is regenerated for H1 and reported but does not change the cap. **Kept in this table because it bounds the SenePy results:** under `senepy_p95` **all six in-band 99th percentiles (118.3–186.8 µm) exceed the cap**, and 7.4–21.5 % of receivers lie beyond 80 µm with 0.2–2.3 % beyond 150 µm, against 0.8–7.0 % and 0.0–0.11 % under the primary `tierA_p95` (`results/phase3/window.csv`). A SenePy-called kernel is fitted on a materially truncated distance distribution, and that limitation is reported wherever a SenePy-called reach quantity is. |
| **P9** | **RESOLVED (PI decision D15): the composition-matched rerun protocol at 5 seeds is frozen — together with its covariate-adjusted counterpart — and is implemented.** | §3.8. At the time of the first draft **no code implemented it**: the phrase appeared only in the planning document (five mentions, no method — no matching variables, no rule, no seeds, no estimand, no output format), and neither `composition_all.py` nor `run_phase5_super.py`'s 5-draw geometry null is the protocol. It is now `code/run_phase8_compmatch.py`, run on M1 and gated on H1. **§15's item is amended to read "the composition-matched rerun protocol at 5 seeds *and its covariate-adjusted counterpart*"** — an addition, not a substitution. All four previously-`TBD` fields are filled in §3.8 from `reports/CS_PHASE8_COMPMATCH.md`. Nine specification ambiguities had to be resolved to implement it at all; they are carried as **D15.1–D15.9** below. |
| **P10** | **The mouse C6 gene sets are promoted into `genesets/`. Supersedes gene-set deviation D15.** | PI decision D5, "tag then promote". `git tag pre-c6-genesets` captures the prior state. Three files changed: B7 38 → 108, strict Tier A 25 → 33, `A_sender_for_secondary_senescence` 55 → 74. Gate re-verified PASS on the authoritative 5,097-gene panel. **Every published Phase 2–5 mouse number was computed pre-C6 and is being re-fitted at task 8.7.** |
| **P11** | **Tier A: strict-33 primary, per-module sets as pre-registered sensitivity.** | PI decision D1. Both are gated in the same run and both gates must pass. Cost, stated against interest: the 33 survivors contain no `CDKN1A`, `CDKN2A`, `TP53`, `LMNB1` or `MKI67` — it is a DNA-damage / p53-effector / replicative-senescence score, not a score on the field's canonical arrest markers. The per-module sensitivity retains all 12 canonical markers for B2, B5, B6 and B7, which is exactly why it is pre-registered alongside. |
| **P12** | **`Lmnb1` is not usable as the D3 re-anchor**, contrary to §7 of the planning document. | It is a member of `B_downstream_arrest` and `B_secondary_senescence`. The primary anchor is an 8-gene proliferation set disjoint from every `A_*.txt` and `B_*.txt`, asserted at run time; `Lmnb1` is reported as a secondary. All anchor decisions are depth-partialled. §3.9. ***[Narrowed 2026-08-27 — see §0.1 item **D-C**. **This row stands, for the reason it gives.** `Lmnb1`/`LMNB1` is excluded because it is a **member** of `B_downstream_arrest` and `B_secondary_senescence` — a gene-set membership fact about our own modules, not a claim about anchor performance — and on H1 it tracks `CDKN1A` closely (`rho_partial_lmnb1` +0.2029 … +0.2449, `results/phase9_h1/deepscence_anchor_h1.csv`) yet remains a secondary on that ground alone. **What is narrowed is the premise that the published anchor misbehaves**: on H1 it is stable 20/20 folds in all 7 sections and the 8-gene proliferation alternative carries almost no signal (ρ +0.0097 … +0.0451) and is itself unstable in SPLN21 (0.75). **On the human arm the published anchor is the better one**, so the D3 re-anchoring is an M1 measure, not a general correction to the tool.]*** |
| **P13** | **The consensus-of-callers anchor is unusable and is not offered.** | It disagrees with the published sign in 9 of 11 M1 sections because SenePy dominates it, and anchoring on it inflates DeepScence-vs-SenePy agreement to 2.793× chance — it manufactures the exact circularity D3 exists to remove, in the opposite direction. ***[Narrowed 2026-08-27 — see §0.1 item **D-C**. Unchanged as a statement about **M1**. On H1 the consensus anchor is also the weakest of the three tested (`rho_partial_consensus` +0.0642 … +0.1695 against `CDKN1A`'s +0.1911 … +0.2540, `results/phase9_h1/deepscence_anchor_h1.csv`), so it is not offered on either arm — but the **reason** the published anchor needed replacing at all is now known to be **partly our ortholog remapping**: P-ii is falsified natively. The circularity this row worries about is not narrowed — DeepScence × `CDKN1A`⁺ is **6.436** pooled on the native H1 panel against 1.255 remapped on M1.]*** |
| **P14** | **`N4_snap` is documented but not implemented.** The rotation family has five variants. | `phase3_null_geom.py:33` (prose) vs `:58` (`ROTATION`, 5 entries). `N4_swap` is itself rotate-then-snap. |
| **P15** | **N3-occ / N4-occ at the specified 5 % tolerance are degenerate and are not a corrected null.** | On a liver section the criterion admits 1–63 of 38,080–108,375 candidate translations, all near-identity: median displacement 27 µm (N3) and 25 µm (N4), inside the 100 µm window and close to the median λ̂ of 12.8 µm; for section 7001 the only admissible translation is the identity. Values quoted (0.951 / 0.896) are the 15 %-tolerance variants; at the literal 5 % they are **0.349 / 0.273** and measure the null's degeneracy, not the effect. ***[Corrected 2026-08-27 — see §0.0 items C-6 and C-7: frozen values are 1–66 admissible offsets, 28 µm (N3-occ) and 25 µm (N4-occ), SFs 0.302 / 0.183, against a pooled λ̂ of 14.7 µm. The verdict is unaffected.]*** Holds for all six sender calls (0.121–0.673). |
| **P16** | **N3-swap reproduces N1 and is not a corrected torus shift.** | Median SF 0.721 vs N1's 0.716, per-fit Spearman ρ = 0.948, median absolute difference 0.0087; conditioning on the N5+N6+zonation block moves it to **0.999**, i.e. it removes nothing the nuisance model does not already remove. The identity is tight for the Tier A percentile calls (ρ 0.92–0.98) and only directional for SenePy (ρ 0.43–0.51), so the footnote reads "reproduces N1 **for the Tier A calls**". |
| **P17** | **The §11 gate now runs automatically after any gene-set or panel change**, not only at freeze. | §4. `code/gate_genesets_guard.py` + a git `pre-commit` hook + a `PostToolUse` hook. The mouse B6 margin is the reason. |
| **P18** | **Mouse B6's margin is one gene, but a one-gene trim does not fail the gate — a two-gene trim does.** Corrects `PREREG_PHASE8_genesets.md` §3.1 / D10. | §4, verified by falsification against the guard. The one-gene statement holds only under the CSV-only 5,006-gene panel definition. |
| **P19** | **`code/build_genesets.py` cannot be re-run as committed and now refuses to run.** | Its `SCRATCH` constant pointed at a dead per-session `/tmp` path, so it globbed zero MSigDB JSONs and would have **silently overwritten `genesets/*.txt` with EMPTY Tier B modules, exiting 0**. It now refuses to run when the archived MSigDB pin is missing, and refuses to revert the promoted C6 sets. |
| **P20** | **`data/raw_h1/` is on disk ahead of the tag.** | §0.4 and §15 require the freeze committed before the human data is downloaded. What has been read is panel membership, file structure, cell counts and coordinate ranges — the §12.1 step-2 screen plus a structural integrity check. **No expression value, senescence score, cell-type assignment or outcome-bearing quantity has been computed.** The pre-registration is still writable without contamination, but it must be committed before Job B or any A3/A5 test runs against these files. |
| **P21** | **The depositors ship usable cell-type annotations, contrary to §12.3.** | `annotations.csv.gz` carries four nested levels (`Level_1`–`Level_4`). This does not replace Job B — the plan's own pipeline still runs — but it gives an external label set to check the marker-based annotation against. The depositors' annotations cover **fewer cells than the matrix** (SPLN07: 239,167 of 249,420), so their QC filter must be characterised before the two label sets are compared. |
| **P22** | **RESOLVED, and the premise inverted: `denoise=False` is a CHOSEN value, not a forced deviation.** | §3.9. **DCA installed and ran** — §6 path 1, not the fallback: DCA 0.3.4 under TensorFlow 2.4.4 / Keras 2.4.3 in an isolated CPython 3.8.19 venv, out of process, on three full M1 sections across both arms plus a three-seed panel. **The install failure that justified the original caveat no longer exists.** What was "we deviated because DCA would not install" is now "**we ran the published default and chose against it, for measured reasons**". `denoise=False` is frozen as PRIMARY and `denoise=True` as the pre-registered published-default sensitivity, satisfying §9's requirement that the denoise state be a stated choice. **The sentence "DeepScence as we can run it on this panel, not DeepScence as published" is withdrawn from this pre-registration.** `reports/CS_PHASE8_D2_DENOISE.md`. |
| **P23** | **The composition-matched protocol as §15 specifies it is INERT — and freezing it alone would have licensed a false claim.** | §3.8. Composition-matched decoys remove **1.6 %** of the pooled naive amplitude (SF 0.9837 [0.973, 0.994]; 3.6 % within cell type). The **same variables, on the same fits, entered as covariates, remove 85.4 %** (`typecomp_adj`, SF 0.1461 [0.052, 0.246]) — **a factor of fifty**. Matching is not failing to balance: max \|SMD\| 0.0916 → 0.0352 and the §8 Test 5 gate passes in **100 %** of matches. It balances sender against decoy without touching the *receiver's* dependence on composition, which is where the confounding acts (`CS_PHASE3.md` §5, verbatim). **The five seeds buy nothing**: at a median match rate of 0.999871 the greedy matcher has no freedom left, and the pooled SF moves by 1.2 × 10⁻⁴ across all five. **Reported, not hidden — §15's parameter as written would certify as "not composition" a gradient that is 66–85 % composition.** The matched number may never appear without the covariate-adjusted number beside it (§10.9). |
| **P24** | **This converges with A7 from an independent direction, and the convergence is the result.** | A7: N2 matched decoys leave the technical gradient **~80 % intact** (−0.061 of −0.070 SD) while the N5 covariate block removes it entirely (+0.007, p = 0.41) — P2. Composition matching: matched decoys remove **1.6 %** where the same variables as covariates remove **85.4 %** — P23. **Two independent analyses, different responses, different confounds, same conclusion: matched-decoy designs systematically fail to remove what covariate adjustment removes on this assay.** The N2-vs-N5 result was already rated the project's strongest contribution by external novelty review; this doubles its evidentiary base and is pre-registered as a claim in its own right, to be tested again on H1. ***[Corrected 2026-08-27 — see §0.0 item C-10. "~80 % intact (−0.061 of −0.070 SD)" is the pre-C6 ratio; frozen `results/phase3/a7_summary.csv` gives **−0.0642 of −0.0744 = 86.3 %**, with N5 still removing the gradient entirely (+0.0038, p = 0.715). The convergence claim is unaffected and slightly strengthened.]*** |
| **P25** | **§17's "composition surrogate share 66–76 %" is confirmed at its lower bound and exceeded at its upper.** | **65.9 %** from receiver cell-type intercepts alone (`type_adj`) reproduces the published 66 % almost exactly and now, for the first time, **has a producer** — `CS_PHASE8_M1_RERUN.md` §7 had recorded that no script emitted the range. But the full composition vector reaches **85.4 %** (`typecomp_adj`), **above the top of the published range**, and 98.5 % for `downstream_arrest`; on the per-module Tier A sets it is 79.6 %, with `downstream_arrest` at 114.6 %. The published interval is also **two different estimators on two different scales** — 66 % is `1 − SF` on β̂, 76 % is a ratio of binned *curve* amplitudes (`CS_PHASE5.md` §4). **The row is split, not patched**, per `CS_PHASE8_COMPMATCH.md` §5.1. Honest statement: **66 % (own cell type) to 85 % (own cell type + neighbourhood composition)**. The composition share is never below 0.52 in any section and never below 0.59 in any module. |
| **P26** | **The published default is not reliably reproducible across seeds. Promote to Results, beside the D3 polarity flip.** | One fixed 20,000-cell subsample, three seeds, nothing else changed (`results/phase8_d2/d2_stability.csv`): **one of three seeds returned a top-5 % sender set perfectly disjoint from the other two** — Pearson *r* **0.5703** and **0.5732**, **Jaccard exactly 0.000**, 2,000 of 2,000 cells changing status — against *r* 0.9955 / Jaccard **0.761** for `denoise=False` across seeds. Nothing about the outlier run looks broken; its denoised matrix has the same global statistics and it returned a normal-looking score distribution, **with no diagnostic that anything was different**. Its internal gene-set metric shows why: at the two agreeing seeds one bottleneck node dominates (0.464, 0.463 vs ~0.15–0.18 for the other), while at the outlier the two are nearly tied (0.373 vs 0.335) — no single senescence axis stands out, and the code picks a node anyway. **This is a second concrete instability in the same published tool**, in the same family as the D3 polarity flip (P12, P13), and belongs in the same section of the paper. It also warrants a note to the DeepScence authors, together with P29. |
| **P27** | **§6's named normalisation estimator was too weak to test its own question. Both configurations are frozen, with `mor`'s inadequacy stated.** | §6 names median-of-ratios (*poscounts*). On this sparse 5K panel it removes only **11.4–23.6 %** of log-depth variance (p90/p10 depth ratio 11.0 → 8.7–10.0) and returns a clean null on four sections — **a null about a weak normalisation, not about the caller**. An added `lib` configuration removes **100 %** (p90/p10 → 1.000) and moves the caller a great deal: depth loading ρ 0.3176 → **0.0819** (7259) and 0.4096 → **0.0289** (7352), i.e. **−74 % and −93 %**, with **46.4 % and 100 %** of the committed sender set no longer called and the call set **broadening** across cell types (biliary/ductular 0.99 → 1.53; stellate, Kupffer, LSEC, T/NK, B all roughly doubling). `results/phase8_d2/d2_normalisation_strength.csv`, `d2_depth.csv`. **`mor` may never be reported as evidence that normalisation cannot move this caller.** `lib` is a declared sensitivity, not a frozen candidate — on 7352 it lands on an axis uncorrelated with the committed one (*r* = −0.017), and that instability is itself the finding: strip depth out and the caller has no single answer. |
| **P28** | **Determinism control — a positive result: the frozen two-section base is re-derivable, not merely intact on disk.** | Re-run at the same seed, the harness reproduces the committed scores at Pearson *r* = **0.9999991** (7239; 24 of 75,384 calls moving) and **0.9999999** (7259; 2 of 114,721). **7259's committed score is `deepscence_sbr.csv`, written 2026-08-20, before two container rebuilds** — so the preserved base reproduces from raw data in the rebuilt environment. This matters more than it would otherwise: this session found **194 untracked files** holding the project's evidence base, so a *demonstrated re-derivation* is worth more than a hash. These two figures, plus the seed-to-seed floor (*r* ≈ 0.996, Jaccard ≈ 0.76), are the floors every D2 effect above is read against. **Caution:** `rho_signed_dz_vs_depth` for the `raw` rows (−0.47, −0.16) is the direction of numerical noise on a shift of 0.0002–0.001 z-units and **must not be quoted**. |
| **P29** | **§4 (D-b)'s stated rationale is refuted by measurement and must be corrected in the paper.** | The planning document says DCA denoising "is precisely the step that would normalize depth — the confound under investigation". **Measured, it RAISES the depth loading by x1.32-1.67** (7239 x1.65, 7259 x1.67, 7352 x1.32 -- 'roughly doubles' was an overstatement by the coordinator and is corrected here): Spearman ρ of score against transcript counts 0.3891 → **0.6404** (7239), 0.3176 → **0.5314** (7259), 0.4096 → **0.5419** (7352) — **three of three sections, both arms**, across a 1.9× spread in median depth (`results/phase8_d2/d2_depth.csv`). It also collapses the global top-5 % call set onto **100.0 % hepatocytes on all three sections** and leaves sender-set agreement with the committed run at Jaccard **0.118–0.280** (56–79 % of the committed senders no longer called). The uncomfortable shape of the result, stated: the denoised score is simultaneously a **better** senescence score by DeepScence's own internal criterion (mean \|corr\| with CoreScence 0.126 → 0.470 on 7259) **and** a more depth-confounded one. **Whatever DCA contributes here, it is not depth normalisation.** |

### D15.1–D15.9 — the nine specification decisions behind the composition-matched protocol

**Why these are in the pre-registration.** §15 names the protocol in five places and specifies its
method in none: no matching variables, no matching rule, no seed values, no estimand, no output
format, and no statement of whether "composition matching" or "composition adjustment" was meant.
Implementing it therefore required nine decisions that are **reconstructions, not readings**. They
are recorded in the producer's docstring as `D15.1`–`D15.9` and carried here verbatim in substance,
so a reader can see exactly which parts of §15 were followed and which were rebuilt.

| # | Ambiguity in §15 | Decision taken | Basis |
|---|---|---|---|
| **D15.1** | "Composition-matched" — matched *how*, and matching *what to what*? | The **N2 design** — senders vs non-sender decoys — with the matching set reduced to composition alone. | N2 is the only sender/non-sender matching design in the plan (Master Plan §22 Step 3, §8 Test 5) and it is already implemented. |
| **D15.2** | Which variables? | Exact stratification on receiver **cell type**, plus 1-1 nearest-neighbour propensity matching on the **20-NN cell-type composition vector**, caliper 0.25 SD. **Nothing else.** | The point is to isolate composition; leaving density, depth and anatomy in makes it the published N2 again. The published N2 set is run alongside as variant `full`, so the isolation is **measured, not asserted**. |
| **D15.3** | Which five seeds? | **20260901–20260905.** | Derived by date from `sasp_phase3.MASTER_SEED = 20260820`, and deliberately outside the range `run_phase3_nulls._expand` can reach, so no rerun accidentally reproduces an existing N2 match. |
| **D15.4** | Does the seed control matching, the bootstrap, or both? | **Both, but separably**: point estimates are computed **before** the bootstrap, so the reported between-seed spread is matching variability and nothing else. | Otherwise "at 5 seeds" would measure bootstrap noise, which the CI already reports. |
| **D15.5** | Pooled, or per receiver cell type? | **Both; pooled is primary.** | A within-cell-type fit has already removed the between-type composition effect, so it cannot be the fit the §17 composition row is about. |
| **D15.6** | Which sender call? | `tierA_p95` on the strict-33 (**PRIMARY**) **and** `tierApm_p95` on the seven per-module sets. | Both pre-registered Tier A variants (P11). They agree: every paired row differs by < 0.07 in SF, and the protocol's own number by 0.0023. |
| **D15.7** | Which sections? | The six §8 Test 3 admissible sections, `sasp_phase3.IN_BAND`. | Every other Phase 3 primary uses them. |
| **D15.8** | Which fits enter a median? | `beta_naive > 0` **and** `beta_base_lo > 0`, identical to `summarize_phase3.sf_table`. | So these medians are directly comparable to the published SF table (§3.5). |
| **D15.9** | **Matching alone cannot answer the scientific question.** | Three seed-free **covariate-adjustment** variants added at the same scopes, same λ̂, same bootstrap: `comp_adj`, `type_adj`, `typecomp_adj`. | `CS_PHASE3.md` §5 is explicit that in this dataset matching balances covariates **without** removing the response's dependence on them. Reporting only the matched result would say "composition matching removes nothing" and a reader would hear "composition is not the confound" — the opposite of the truth. **This decision is what PI decision D15 ratified into the frozen list.** |

**Recorded silence.** The documents never state the protocol's *output*, so the §17 row it was
meant to fill had no producer at all before this run — reached independently by
`reports/CS_PHASE8_M1_RERUN.md` §7. See P25.

**Calibration, and why the inertness was foreseeable.** On synthetic tissue, N2-style matching
returns **0.934 with a planted real effect and 0.775 with none** (`CS_PHASE3.md` §6). A design
whose two answers are that close cannot discriminate, and a composition-only matching set has
strictly less to work with than N2 did.

### P1 in full — the caller-agreement headline moved twice

**Quantity:** depth- and type-matched top-5 % overlap of the three non-circular Tier A pairs
(Tier A vs SenePy, vs DeepScence, vs `Cdkn1a`⁺), pooled as Σ`n_both` / Σ`exp_both_stratified`,
z from the pooled `sd_both_stratified`. The circular DeepScence-vs-`Cdkn1a`⁺ pair is excluded from
every pooled number, per `BIO_PHASE3.md` §4.4. All three rows recomputed for this draft directly
from `results/phase3/caller_agreement_matched_significance_{11sections,verify2sec}.csv`.

| Configuration | Pooled | z | p | Band | Median | Above chance |
|---|---|---|---|---|---|---|
| **(a)** 2 sections, **pre-C6** Tier A — *the published base* | **1.030×** | 1.27 | 0.203 | 0.932 – 1.221 | 1.000 | 3 / 6 |
| **(b)** 11 sections, **pre-C6** Tier A — *the like-for-like coverage comparison, and the basis of the 8.4 gate decision* | **1.118×** | 11.49 | 1.44 × 10⁻³⁰ | 0.700 – 1.711 | 1.110 | 20 / 33 |
| **(c)** 11 sections, **post-C6** Tier A — *the frozen configuration, task 8.7, 06:07 UTC.* **PROVISIONAL** | **1.212×** | 20.62 | 1.84 × 10⁻⁹⁴ | 0.751 – 2.198 | 1.179 | 26 / 33 |

**(a) → (b) is the answer to the 8.4 gate question — "did full DeepScence coverage move the
headline?" — and it is the pair that must be cited for the restatement.** Both bases are pre-C6, so
it holds the gene sets fixed and varies only DeepScence coverage, which is exactly the comparison
§5 asked for. **The like-for-like answer is 1.030 → 1.118, p = 1.44 × 10⁻³⁰.**

**(c) answers a different question and must never be quoted as the gate result.** It mixes the
coverage change with a **sender-definition** change (the promoted C6 33-gene strict Tier A), and
the two cannot be separated inside it. Both belong in the pre-registration, and the distinction is
stated here explicitly rather than left for a reader to infer from mtimes. **The
motivating sentence "their top-5 % calls overlap at 0.93–1.22× of chance … i.e. they are
statistically independent" is dead** and the defensible restatement — "weakly but genuinely
dependent, in a direction each pair's depth loading explains" — is drafted verbatim in
`reports/CS_PHASE8_CALLERS.md` §3. Tier A vs DeepScence is above chance in **11 of 11** sections
under (b), which breaks the independence claim on its sign pattern alone, without pooling.

**(b) → (c) is the C6 gene-set change, not coverage, and every figure in it is PROVISIONAL.** The 8.7 re-run rebuilt
`tierA_score` on the promoted 33-gene strict Tier A, so every Tier A pair moved: vs DeepScence
1.248 → **1.288** (11/11 above chance, 10 significantly); vs `Cdkn1a`⁺ 1.171 → **1.471** (11/11);
vs SenePy 0.914 → **0.972**, which is **no longer significantly below chance** (z −1.63, p = 0.104).
`results/phase3/caller_coverage_gate.csv`, rewritten 06:09 UTC.

⚠ **Two cautions on `results/phase3/caller_coverage_gate.csv` as it now stands.**
 (i) **Confirmed by the coordinator.** Its "2-section (published base)" rows come from
 `*_verify2sec.csv`, last written **06:03** with the **pre-C6** sender scores, while its 11-section
 rows were written **06:07/06:09** from the **post-C6** ones. **The file therefore compares two
 different gene-set configurations in one table, and every pooled number in it — including the
 1.212 recorded above — is provisional until it is regenerated.** The M1 re-run agent has been
 directed to recompute both bases under the same sender definition and to label the basis
 explicitly in the output. Row (b) above is the correct like-for-like comparison and was recomputed
 here for that reason.
 (ii) Under (c) the claim "one pair sits *below* chance in all eleven sections" — which is load-
 bearing for the restatement, because it is what distinguishes a shared **technical** variable from
 a shared latent state — **weakens to 0.972, n.s.** If (c) is confirmed when 8.7 finishes, that
 sentence must be re-derived or dropped. **It must not be carried over from the pre-C6 text.**

⚠ `reports/CS_PHASE8_CALLERS.md` §2.1 states **22/33** sections above chance for row (b); the file
gives **20/33** (`ratio_stratified > 1`, equivalently `z > 0`). Use 20/33.

---

## 14. Open items the PI must close before tagging

**Closed since the first draft:** the fitting window (**PI decision D16** — frozen as the 100 µm
literal, §3.2 / P8); the composition-matched rerun protocol (**PI decision D15** — implemented,
and frozen **together with its covariate-adjusted counterpart**, §3.8 / P9 / P23–P25, with the nine
reconstruction decisions carried as D15.1–D15.9); and **D2 `denoise`** (task 8.5 — DCA installed
and ran, `denoise=False` frozen as a *chosen* value with `denoise=True` as the published-default
sensitivity, §3.9 / P22 / P26–P29). All three are recorded above; none is open.

**Two items D2 adds to the paper rather than to this table**, flagged so they are not lost:
the **seed instability** (P26) is promoted to Results beside the D3 polarity flip, and both it and
the refutation of §4 (D-b) (P29) warrant a note to the DeepScence authors.

| # | Item | Where |
|---|---|---|
| 1 | **Fill the tag hash and the `pre-c6-genesets` hash.** | §1 |
| 2 | ~~Transcribe the composition-matched protocol's four fields~~ — **DONE.** §3.8 is filled from `reports/CS_PHASE8_COMPMATCH.md` and re-verified against `results/phase3/compmatch_reruns.csv`. What remains is editorial: **§17's "composition surrogate share 66–76 %" row must be split**, not patched — see P25. | §3.8, P25 |
| 3 | **Confirm the primary outcome and the R1–R5 thresholds** as stated, or amend them now. They are numeric and absolute by design. | §5, §6 |
| 4 | **Confirm the SenePy recommendation** — demote to sensitivity on H1 and promote `CDKN1A`⁺ into the primary trio, or leave as is. Flagged, not taken. Note that SenePy now carries **two** independent caveats: no spleen hub (P4) and the 100 µm window truncating 7–22 % of its receivers (P8). | §3.7, P4, P8 |
| 5 | **Accept the three corrections to already-written text**: the mouse-B6 one-gene statement (P18), the 22/33 → 20/33 count (P1), and the mixed-configuration `caller_coverage_gate.csv` (P1 caution i — since confirmed by the coordinator and being regenerated). | §4, P1 |
| 6 | **Re-run the gate and re-read every remaining PROVISIONAL value once 8.7 finishes** — in particular A7 and the corrected N3/N4, both computed on the pre-C6 sender calls, and every pooled number in `caller_coverage_gate.csv` — and replace every **PROVISIONAL** mark. | §0, §5, P1 |
| 7 | **Decide whether the "one pair below chance in all eleven sections" sentence survives.** Under the post-C6 configuration Tier A vs SenePy is 0.972, n.s. That sentence is load-bearing for the restatement — it is what separates a shared *technical* variable from a shared latent state. | P1 caution (ii) |

***[Post-freeze status note, 2026-08-27 — see §0.1. `CS_PHASE9_H1_AUDIT.md` §13 item 2 and item 7 name **three** PI decisions
needed before the Phase-10 fits. **Two are now taken:** which label family the Tier A sender call is defined on (**D-B** —
merged is primary on H1, the frozen fine-label call is the sensitivity), and whether DeepScence stays in the H1 caller set
(**D-A** — it stays, as a five-seed consensus with its dispersion reported, not as the frozen single-seed score).
**The third is still open and is item 4 above:** whether SenePy stays in the primary trio. It now carries **four** caveats on
H1, not two — no spleen hub (P4), the 100 µm window truncating its receivers (P8), a within-type depth enrichment of
**Q5/Q1 = 28.5–228.1×** across the 7 sections, recomputed here from
`results/phase9_h1/caller_within_type_depth_bias.csv` (`CS_PHASE9_H1_AUDIT.md` §10.2 quotes 28.5–224.7 from the unrounded
values; the committed file is rounded to 3 dp and gives 228.1 for SPLN24) — command:
`python3 -c "import pandas as pd; d=pd.read_csv('results/phase9_h1/caller_within_type_depth_bias.csv'); s=d[d.caller=='senepy_score'].pivot_table(index='section',columns='within_type_depth_quintile',values='enrichment'); print((s['Q5']/s['Q1']).round(1))"`, and the v1/v2 hub-release discrepancy
(`CS_PHASE9_H1_AUDIT.md` H4, `results/phase9_h1/senepy_surrogates_v1_v2.csv`). **Item 7 is answered on the human arm**: Tier A
× SenePy is **0.874, z = −7.96, below chance in 7 of 7 H1 sections** (`results/phase9_h1/caller_agreement_pooled.csv`), so the
sentence survives on H1 even though the post-C6 M1 configuration weakened it to 0.972 n.s. **None of D-A–D-D closes item 4.**]***

***[Confirmed 2026-08-27 — see §0.2 item **F7(b)**. The SenePy Q5/Q1 range quoted in the note above, **28.5–228.1×**, is re-confirmed as the committed file's value and is the one to use; the audit's **28.5–224.7** is superseded. Per section, `results/phase9_h1/caller_within_type_depth_bias.csv`: SPLN07 74.6, SPLN14 30.8, SPLN21 50.6, **SPLN24 228.1** (4.105/0.018), SPLN30 **28.5**, SPLN43 98.0, SPLN44 70.3. **Propagation checked, not assumed:** `grep -rn "224\.7" reports/ code/ figures/` returns exactly two hits — `CS_PHASE9_H1_AUDIT.md:765` (the audit's own table, left for the audit's owner) and this note, where it appears only as a citation beside the corrected value. **The wrong value reached no claim, no figure and no code path**, and no conclusion turns on the difference. **Item 4 is still open** — F7(b) corrects one of its four caveats' digits, it does not decide SenePy's status.]***
