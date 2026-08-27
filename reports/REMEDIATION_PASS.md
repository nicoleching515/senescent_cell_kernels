# REMEDIATION PASS — 2026-08-27

Remediation of the defects found by `reports/AUDIT_NUMBERS_FINAL.md` and
`reports/AUDIT_CLAIMS_CITATIONS.md`. **Documents and `references.bib` only.** Nothing under
`results/`, `figures/`, `code/`, `genesets/` or `data/` was modified by this pass; `data/raw_h1/`,
`results/phase9_h1/`, `data/processed_h1/` and `code/h1_*` were not read or touched. No analysis
was re-run. No tag was created or moved; nothing was pushed.

**Method.** Every replacement number below was derived from a file in this session, with the
command recorded beside it. Nothing was quoted out of a report — that failure mode is what this
pass exists to remove, and §6 records the one place it nearly recurred.

---

# 0. STILL LIVE AFTER THIS PASS — read this first

## 0.1 The figures guard does **not** pass, and that is not this pass's doing

**A task-8.7 figure regeneration is in progress in this working tree, and the guard is failing
against a manifest that has not been re-snapshotted.** The guard passed at the start of this pass
("OK: all 52 committed figures match", exit 0) and again after the PRIORITY 3 commit. It began
failing after commit **`39e7791` "figure2e: write the producer that never existed, and fix what it
drew"** — another agent's work — and the count has grown as that work continues:

```
17:24  CHANGED: figure2e.{pdf,png}                                     (2 artefacts)
...    CHANGED: figure2a.{pdf,png}, figure2a_amplitudes.csv, figure2b.{pdf,png},
                figure2c.{pdf,png}, figure2e.{pdf,png}, figure3.{pdf,png},
                figure3_data.csv, figure4.{pdf,png}, figure4_data.csv,
                figure4_supp_ncem_lengthscale.{pdf,png}                (17 artefacts)
exit=1
```

**I did not touch `figures/` and did not run `--snapshot`.** Reconciling the manifest belongs to
whoever owns task 8.7 regeneration — it is their re-snapshot to make, and the count will keep
moving until they make it. **Until then the guard's output must not be quoted as a pass anywhere**,
and the `WRITING_PACK.md` figure-guard row now says exactly that.

Related, and the reason §6 exists: **`figures/figure4.png` was regenerated mid-pass** — md5
`d44fac63411d6c30a42c40894a287f17` at first observation, **`fdead29871b61481e297951dfea75b3d`**
minutes later.

## 0.2 Seven numbers-audit items outside my brief are still live in `WRITING_PACK.md`

These are the audit's **"true but misleading"** class (§1 of `AUDIT_NUMBERS_FINAL.md`) plus R6 and
R9. They were **not** in this pass's scope and I did **not** edit them. Each is one edit away; the
corrected value is derived below so nobody has to re-derive it.

| audit | site | what is wrong | derived correction, and the command |
|---|---|---|---|
| **R6** | `WRITING_PACK.md:608` | *"pooled CI half-width **±0.018 SD**"* sits beside the `all_controls` / `base` number, where the pooled half-width is **0.0562** — **three times** the quoted figure. 0.018 is the *conditioned probe/codeword* row | Over the 25 control rows the pooled half-width `(clustered_hi − clustered_lo)/2` spans **0.016–0.056**. Name the response and design, or give the range. `python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); c=d[d.response!='BIOLOGICAL MODULES (reference)']; h=(c.clustered_hi-c.clustered_lo)/2; print(h.min(), h.max())"` → `0.0155 0.0562` |
| **R9** | `WRITING_PACK.md:96` | *"IQR is [7.0, 50.0] µm for **every one of these**"* — false for the two **interior** definitions, whose quartiles cannot sit at the rails by construction (interior of 315: [10.67, 28.68]; interior of 153: [10.63, 25.43]) | True for the three definitions that matter, including the authoritative one (315 fits, median **14.7321**, IQR **[7.0, 50.0]**, **60 %** railed — re-derived this session). Scope the sentence to the non-interior definitions |
| **M1** | `WRITING_PACK.md:749` | *"All four rows below are the **same 42 fits / 33 reportable**"* — `comp` and `full` are pooled at **210 / 165** (= 42 × 5 seeds and 33 × 5); only the adjusted variants are 42 / 33 | Say *"the same 42 fits / 33 reportable **per seed**; the matched rows are pooled over the five frozen seeds"*. `python3 -c "import pandas as pd; r=pd.read_csv('results/phase3/compmatch_reruns.csv'); r=r[(r.row_type=='summary')&(r.call=='tierA_p95')&(r.scope_kind=='pooled')]; print(sorted(set(zip(r.variant,r.n_fits,r.n_reportable))))"` → `comp 210/165, full 210/165, comp_adj 42/33, type_adj 42/33, typecomp_adj 42/33` |
| **M2** | `WRITING_PACK.md:612` | the reportable filter *"admits **4.8 %**"* — the value is **0.0485**, i.e. **4.9 %** at the precision the adjacent "3.0–13.3 %" is given at. Truncation, not rounding | write **4.9 %**, or **0.0485** |
| **M3** | `WRITING_PACK.md:397, :978` | Poisson identity *"slope −0.525, r² 0.9843"* is **subset-dependent and the subset is not named**; `summary_phase3.txt` §10 gives four (all sender definitions −0.525/0.9843 n=77; `cdkn1a_pos` −0.486/0.9982; Tier A percentile −0.517/0.9958; in-band −0.529/0.9833). §5.8 offers it as *the cross-arm geometric prediction*, so it must not move | *"−0.525, r² 0.984, over all **77 section × sender-definition cells**"* |
| **M7** | `WRITING_PACK.md:476` | *"no family beats it in more than **54 %** of folds"* — gaussian is **0.544** | "more than **55 %**", or quote **0.544**. `grep -n gaussian results/phase5/summary_phase5.txt` → `0.091  0.544  0.115` |
| **M5** | `WRITING_PACK.md:667` | **largely already handled** — the pack does give 0.362 with its 0.308–1.070 range at :667. The remaining exposure is :670 and :1282, which state 0.362 bare | when 0.362 is stated bare, attach *"median over 22 section × call rows; 0.308–1.070"* |

## 0.3 Deliberately not fixed

- **`lotwick1982methods`'s p.410 quotation** is still verified only as Mrkvička quotes it (paywall).
  The entry says so; that is honest and stays.
- **`hodges2010adding` / `zimmerman2022deconfounding` verification-line upgrades** — the claims
  audit lists these as *optional*, and I did not re-retrieve those two records myself. I do not
  restate another auditor's retrieval as my own. Left as they are.
- **Anything requiring a re-run.** Nothing in this pass needed one, and none was performed.

---

# 1. PRIORITY 1 — the live error inside the frozen pre-registration

**Commit `96a0cee`** (`reports/PREREG_PHASE8.md`, `reports/COMPLETED_TASKS.md`).

Handled exactly as §0.0's C-1/C-2/C-3 were: **no original wording deleted anywhere**, three new
dated rows in the §0.0 correction block, and an inline dated marker at each site pointing back.

**The file the corrections come from.**

```bash
python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/a7_summary.csv'); print(d[d.design=='n6n5'][['response','frac_CI_excludes_zero']].to_string(index=False))"
#   neg_control_codeword 0.091 | all_controls 0.103 | genomic_control 0.109
#   neg_control_probe    0.145 | neg_probe_rate 0.164
stat -c '%y' results/phase3/a7_summary.csv          # 2026-08-27 09:06:16  (frozen)
stat -c '%y' results/phase3_pre_c6/a7_summary.csv   # 2026-08-27 05:19:12  (pre-C6)
```

| new row | site | what was wrong | correction |
|---|---|---|---|
| **C-9** | **P3** | The one that bites. Body lists `0.091 / 0.103 / 0.109 / **0.127** / 0.164`, derives *"the clean-null subset is **9.1–12.7 %**"*, and **instructs** *"Quote the range as **9–13 %**"*. 0.127 is the pre-C6 `neg_control_probe` value | Frozen probe is **0.145**, so the clean-null subset is **9.1–14.5 %**. **The instruction now reads: quote 9–15 % on the four count-based responses, with 16 % as the `neg_probe_rate` outlier.** P3's header ("9–16 %") was already correct and is unchanged. Also recorded: the pre-C6 list is **mis-ordered** — pre-C6 `all_controls` is 0.091 and `neg_control_codeword` 0.109, the reverse of frozen — so quoting it as an ordered per-response list was wrong even before the vintage issue |
| **C-8** | **P2** | Every digit is the 05:19 file, **and** it closes with a stale *instruction*: *"Magnitudes are PROVISIONAL — A7 was run at 05:19 … **A7 must be re-run after 8.7**"* | Frozen `all_controls`: naive **−0.0744 [−0.1306, −0.0182], p = 0.0145**; N2 **−0.0642 [−0.1113, −0.0172], p = 0.0124**; N5 **+0.0038 [−0.0186, +0.0261], p = 0.715**; the 40 probes **−0.0225 [−0.0527, +0.0078], p = 0.129** — **still flat, so the pre-registered primary A7 response still passes**; biological **0.3120** naive / **0.0795** conditioned. **A7 *was* re-run: the frozen file is 09:06.** The instruction is discharged and the magnitudes are no longer provisional. All three qualitative findings unchanged |
| **C-10** | **P24** | *"N2 matched decoys leave the technical gradient **~80 % intact** (−0.061 of −0.070 SD)"* — the pre-C6 ratio. C-1 fixed the same quantity to 86 % at §10.1 and left P24 at 80 % | **−0.0642 / −0.0744 = 86.3 %.** P24's convergence argument is **strengthened**, not weakened |

The block's closing paragraph was amended: it previously said the corrections were all "pre-C6 →
frozen digit substitutions", which C-8 and C-9 are not — they change *instructions*.

**Swept afterwards** for every known pre-C6 pattern in the file (`0.127`, `−0.070`, `0.0697`,
`~80 %`, `12.8`, `0.0549`, `0.0337`, `0.0177`, `0.2914`, `0.0356`, `9–13 %`). Every surviving hit
is inside a correction row or under a dated marker (§3.6 under C-2, §10.1 ¶1 under C-1, P15 under
C-6/C-7).

**`reports/COMPLETED_TASKS.md`:161 — retracted in place.** Row 78's *"Residual stale digits: **0**"*
was false. It is struck, with the reason recorded against its author: the zero was true only of the
three patterns searched (codeword / genomic / probe amplitudes), not of the file. *A residual-count
is only as good as its pattern list, and a zero must state the patterns it was searched over.* The
row had also lost its line break and was fused to row 83; that is repaired.

---

# 2. PRIORITY 2 — the forbidden claim re-introduced by the correction pass

**Commit `5fef337`** (`references.bib`).

`references.bib:783`, the `% SUPPORTS:` block on `moses2023voyager` **added by the citation audit
itself**, ended: *"Report our own Moran's I on the controls next to the kernel amplitude so the
reader can see the two tests disagree."* That is the exact sentence struck from
`NOVELTY_ASSESSMENT.md` §2.1 point 3 as falsified — re-introduced into a new file after being
struck from the old one.

**Replaced with the power framing, not the orthogonality framing.** All values re-derived:

```bash
# the two statistics AGREE -- |Moran I| vs |A7 naive amplitude|, 12 control+module fields
python3 -c "
import pandas as pd; from scipy.stats import spearmanr
d=pd.read_csv('results/moran/moran_pooled.csv'); s=d[d.kind.isin(['control','module'])]
print(len(s), spearmanr(s.I_raw_mean.abs(), s.a7_base_mean.abs()), spearmanr(s.I_ctcentred_mean.abs(), s.a7_base_mean.abs()))"
# 12  rho=+0.8951 p=8.37e-05   rho=+0.9441 p=3.93e-06   -- reproduces moran_verdict.txt exactly
```

Note for the record: the numbers audit describes these as "section-clustered mean per field, knn6
raw" but **omits the absolute-value step**. Without `.abs()` they come out at +0.566 / +0.657. The
statistic is `|Moran's I|` against `|A7 naive amplitude|` (`code/summarize_moran.py:183-184`).

```bash
python3 -c "
import pandas as pd, numpy as np; d=pd.read_csv('results/moran/moran_kernel_power.csv')
print(len(d), np.median(d.dI_as_frac_of_I_controls), d.dI_as_frac_of_I_controls.min(), d.dI_as_frac_of_I_controls.max())
print(np.median(d.beta_min_visible_to_moran), d.beta_min_visible_to_moran.min(), d.beta_min_visible_to_moran.max())"
# 22 rows; dI/I median 0.00826 (2.13e-05 .. 0.0608); beta_min median 0.3622 (0.3082 .. 1.0697)
```

The block now states: the whole A7 gradient moves the controls' Moran's I by a median **0.83 %** of
its observed value, and Moran's I cannot resolve a kernel amplitude below a median **0.362 SD**
(range 0.308–1.070) — larger than every A7 amplitude the project reports, and 12× the conditioned
amplitude. A near-zero global Moran's I is therefore **consistent** with the A7 gradient, not
evidence against it.

A `% FORBIDDEN` note now records the re-introduction by name so it cannot happen a third time.

**Swept:** every surviving *"two tests disagree"* in the repo is a prohibition, a strikethrough, or
a correction record. None reads as an instruction.

---

# 3. PRIORITY 3 — the two zero-marker reports

**Commit `9ee4172`** (`reports/BIO_PHASE3.md`, `reports/BIO_PHASE2.md`).

Findings were **not** rewritten. Each file now carries a head banner that names the offending text
**line by line**, plus an inline dated strike **at every site** — because six markers elsewhere were
found scoped too narrowly, and a banner that names a section without naming the sentence is how
that happens.

**Frozen replacements, all re-derived:**

```bash
python3 -c "import pandas as pd; print(pd.read_csv('results/phase3/caller_coverage_gate_headline.csv').to_string(index=False))"
python3 -c "import pandas as pd; g=pd.read_csv('results/phase3/caller_coverage_gate.csv'); print(g[['basis','pair','n_sections','ratio_min','ratio_median','ratio_max','pooled_ratio','pooled_z','n_sections_above_chance']].to_string(index=False))"
python3 -c "
import pandas as pd; d=pd.read_csv('results/phase3/caller_within_type_depth_bias_11sections.csv')
p=d.pivot_table(index=['section','caller'],columns='within_type_depth_quintile',values='enrichment')
print((p.Q5/p.Q1).groupby('caller').agg(['min','max','count']))"
```

- 11-section post-C6 (FROZEN): band **0.751–2.198**, pooled **1.212**, z **21.92**, p **1.84e-106**,
  above chance in **35 of 44**.
- 2-section post-C6 (FROZEN): pooled **1.131**, z **5.80**, p **6.5e-09** — vs pre-C6 **1.040**,
  z 1.76, p 0.078 **on the same two sections**.
- `deepscence vs cdkn1a`, 11 sections: **0.963–2.849, median 1.071, pooled 1.255**, above chance in
  **7 of 11**.
- `tierA vs senepy`, 11 sections, C6: **0.972, z −1.63, p 0.104, above chance in 4 of 11.**
- `senepy vs deepscence`, 11 sections: 0.332–2.150, median 0.506, pooled **0.737, z −15.08**.
- Depth camps, Q5/Q1, **11 of 11 in direction for every caller**: SenePy **10.58–41.74**,
  `Cdkn1a`⁺ **4.19–42.36** (top); Tier A **0.146–0.317**, DeepScence **0.218–0.795** (bottom).

## `reports/BIO_PHASE3.md` — five sites marked

| site | struck | why |
|---|---|---|
| §0 headline item 4 | *"conditioning … does **not** make them agree (0.93–1.22× of chance)"* | items 13/14 |
| §4.4 | *"Four of six pairs sit at 0.93–1.22× — statistical independence"* | items 13/14; a **two-section** property — at 11 sections only 2 of 6 fall in that band |
| §4.4 table | the `2.85 / 1.51` cell | item 8; **the two published values are the two largest of the eleven** |
| §4.4 | *"concordant in sham … **anti**-concordant in SBR"* | item 15 — **no arm effect**; one anomalous section (**7250**) |
| §4.4 | *"The only pair consistently above chance is DeepScence vs `Cdkn1a`⁺"* | **not in either audit — found by my own sweep.** At 11 sections `tierA vs cdkn1a` (1.471) and `tierA vs deepscence` (1.288) are above chance in **11 of 11**, while `deepscence vs cdkn1a` is only 7 of 11. The circularity *explanation* stands; the "only pair" claim does not |
| **§4.5 in full** | the blockquote headed *"The finding, stated for the paper"* | **STRUCK IN FULL.** Drafted as manuscript text; welds items 13, 14 **and** 16 into one paragraph. Replacement named: `SUBMISSION_PATCH_2026-08-29.md` §2.1. The text is kept unaltered beneath the strike **as a record, not a draft** |

**What still stands** is stated explicitly: §§1–3 and §5 are untouched; the §4 depth-camp
*mechanism* survives coverage 11/11 (only its four-section magnitudes, 2.3–2.5× / 8–11×, are
superseded); and — stated plainly — **absence of a correction is not a verification.**

## `reports/BIO_PHASE2.md` — banner + §4.3 correction box in the file's own §4.2 style

Struck: §0 headline item 2's *"do not agree … at better than chance"*, §4.3's *"**The three
definitions are statistically independent**"*, and *"two of three pairs overlap less than random"*.
The §4.3 heading is annotated rather than rewritten.

**The decisive number:** on the frozen C6 33-gene Tier A over **the same two sections**, **none** of
the three pairs is below chance — `tierA vs senepy` **1.007** (was 0.935), `tierA vs cdkn1a`
**1.300** (was 1.017), `senepy vs cdkn1a` **1.168** (unchanged; it does not involve Tier A). **The
sender-set defect, not coverage, is what killed independence.** §4.3's closing sentence — that λ is
a property of the sender definition chosen — is called out as **unaffected and still standing**.

---

# 4. PRIORITY 4 — the wrong citation on load-bearing prior art

**Commit `5fef337`** (`references.bib`). Both records retrieved by me in-session; nothing rests on
the other audit's retrieval.

## `ren2025systematic` — wrong senior author

```bash
curl -s 'https://api.crossref.org/works/10.1038/s41467-025-64292-3'          # 19 authors
curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=41107232&retmode=xml'
```

Both agree exactly: **19 authors**, #1 Ren Pengfei, **#8 Zhang, Zongxu**, #18 Zhang, Zhe,
**#19 Zeng, Zexian**. Volume **16**, issue **1**, article number **9232**, PMID **41107232**.

The entry read `... and others and Zhang, Zongxu`, which under the `and others and X` convention
**asserts Zongxu Zhang as the senior author**. He is eighth. Corrected to
`... and others and Zhang, Zhe and Zeng, Zexian`; `number`, `pages` and `pmid` added. A
`% CORRECTED` note records both retrievals and why it mattered — this is the only peer-reviewed
prior-art citation for the negative-control-probe diagnostic.

## `hughes2025senpred` — the `% AUDIT` note was the error

```bash
curl -s 'https://api.crossref.org/works/10.1186/s13073-024-01418-0'
curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=39810225&retmode=xml'
```

The title in the entry matches Crossref and PubMed **word for word**, including the trailing
*"…for the detection of an in vivo senescent cell burden"*. The note claiming the title was
incomplete is removed, replaced by a dated record of the two retrievals — leaving it standing
invited someone to "fix" a correct title.

`grep -c '^@' references.bib` = **43**, unchanged.

---

# 5. PRIORITY 5 — the remaining live items

## 5.1 From the numbers audit (`WRITING_PACK.md`) — commits `388edec`, `4aab386`, `3e4ffbf`

| audit | correction | derivation |
|---|---|---|
| **R3** | figure4.png md5 — **hash removed, not replaced.** See §6 | `md5sum figures/figure4.png` |
| **R5** | ~~"196 cells, every Tier A × Tier B cell zero"~~ → **224 cells across three panels; the 14 `A_SENDER_FINAL_strict` × Tier B cells are all zero in both arms**. The gate PASS is real but describes **14 of 224**; 55 of the 112 `frozenAxB` cells are non-zero (the per-module *sensitivity* sender sets, required disjoint only from their own module). **196 is reproducible from no subset of the file** | `python3 -c "import pandas as pd; d=pd.read_csv('figures/figure_gs1_intersection_matrix_data.csv'); f=d[d.panel=='frozenAxB']; print(len(d), d.panel.value_counts().to_dict(), (f.n_frozen!=0).sum(), f[f.row.str.contains('strict')].n_frozen.abs().max())"` → `224 {frozenAxB:112, BxB:84, A0xB:28} 55 0` |
| **R2** | ~~"never outside ±0.010 of nominal"~~ → **"within 1.6 Monte-Carlo SE of nominal everywhere"**. The gloss contradicted the range in its own sentence: two of eight cells are below 0.040 and the largest shortfall is **0.0175**. MC SE = √(0.05·0.95/400) = **0.0109** → 1.6 SE. **The range 0.033–0.060 is correct and confirmed** | `rs_count_reject_05` = 0.0325, 0.0550, 0.0350, 0.0600, 0.0400, 0.0425, 0.0525, 0.0550 |
| **R1** | ~~"`denoise=False` … Jaccard 0.76–0.99"~~ → **"Pearson r = 0.996 and top-5 % Jaccard 0.76"**. The 0.99 is the **Pearson r**; the sentence merged a correlation and a set overlap | `d2_stability.csv` has exactly one `denoise=False` pair: `raw_seed0 vs raw_seed1`, r **0.99553**, `top5_jaccard` **0.7606**. No Jaccard of 0.99 exists in the file |
| **R7** | four line citations → **`:52` / `:143` / `:204` / `:52–72`** (were :59 / :169 / :239 / :59–93). `:59` is **blank**, `:169` is `def fit_cell`, `:239` is `keys = list(b)`. **Substance fully confirmed and unchanged** | `grep -n 'WINDOW_UM\|LAM_LO_FLOOR\|N_LAM\|lam_grid\|d_obs <= \|lam_railed' code/run_phase3_nulls.py` |
| **R4** | ~~"the guard now covers everything"~~ → **it covers the 52 committed artefacts**; `figures/revised_candidates/` holds **10 untracked files (9 `_REVISED` artefacts + a README)** it cannot see | `git ls-files figures/` = 52; `ls figures/revised_candidates/ \| wc -l` = 10; `git ls-files figures/revised_candidates \| wc -l` = **0** |
| **R10** | ~~"0.118/0.05 = 2.36; `CORRECTIONS.md` rounds to 2.35"~~ — backwards. **0.1175 / 0.05 = 2.35 exactly**; 2.36 comes from rounding 0.1175 twice | `torus_tile4x4_reject_05` max = **0.1175** (irregular, s = 0.30) |

## 5.2 From the claims audit — commits `9308cf8`, `388edec`

**`2.4×` → `2.35×`** at the four `WRITING_PACK.md` sites (`:893`, `:896`, `:1280`, `:1513`) **and**
the two in `PLAN_UPDATE_D12_D13.md` (`:216`, `:277`) — the latter named by the audit but not by my
brief; `:277` quoted the correct file value `0.1175` and then rounded it wrong in the same
sentence. `:1513`'s *"either rounding is defensible"* is retracted: only one is.

**`Phase7_Minimal_Human_Replication (1).md`** — the banner named §4 **for D-b only**, leaving D-d
live. It now names **§4 (D-b) and §4 (D-d)**. `:187`'s **"1.51–2.85×"** struck (item 8) with the
11-section values, and the banner's *own* **"roughly doubles"** corrected to **×1.32–1.67**
(`d2_depth.csv`: 0.3891→0.6404 ×1.65; 0.3176→0.5314 ×1.67; 0.4096→0.5419 ×1.32).

**`README.md`** — the table header *"Moves across callers that **agree at chance**"*, three
paragraphs below the box that restates them as *above* chance, now reads **"agree weakly but
genuinely above chance"** with the frozen 1.212× / z 21.92 inline. **Additionally (not in either
audit, found while there):** the same paragraph carried the refuted rationale *"DCA denoising is
precisely the step that would normalise depth"*. Struck, with P29's measurement beside it.

**`reports/CS_PHASE8_M1_RERUN.md:307`** — the **FILL-IN** row staged for the manuscript called a
pooled number a *probe* number and asserted **NOT FLAT**. Renamed to **"Pooled negative-control-
feature kernel"**; it now leads with **the 40 named probes are FLAT (−0.0225 [−0.0527, +0.0078],
p = 0.129) and are the pre-registered primary null, so on its own primary response M1's A7
PASSES**, then gives the pooled −0.0744 attributed to the codewords and genomic controls. The
FILL-IN instruction now carries an explicit prohibition against filling it in the old way. *(The
row's digits were already frozen; only the naming and the verdict were wrong.)*

**`reports/NOVELTY_ASSESSMENT.md`** — banner's *"Everything else in this report stands"* **struck**
(it does not), and the banner now enumerates the three items that do not stand:
- **§U1's drop-in manuscript paragraph** struck and replaced. Two defects: the antecedent of *"we
  instead use **them**"* is *negative control probes* while the number is the **pooled** response
  (item 2), and the digits are the **pre-C6 05:19** vintage. The replacement gives the pooled
  **−0.0744 [−0.1306, −0.0182], p = 0.0145**, attributes it to codewords (−0.0604) and genomic
  controls (−0.0307), states the **40 probes flat**, and gives N5 **+0.0038, p = 0.715** / N2
  **−0.0642, p = 0.0124 (86.3 % left)**.
- **§U1 and §4 O2's `neg_probe_rate` "+0.014, p = 0.079"** → frozen **+0.0113 [−0.0085, +0.0310],
  p = 0.232**. The claim survives — not significant either way — but **p moved by a factor of
  three**.
- **§U2's Mrkvička misquotation** (item 32), which was inside quotation marks. **I retrieved
  arXiv:1911.00240 and extracted it myself** (`curl -sL https://arxiv.org/pdf/1911.00240`, then
  decompressed the content streams). §2.1.4 reads verbatim: *"The approaches using minus correction
  and variance correction **can be applied in case of general (compact) observation windows**."*
  The string *"irregular window"* **does not occur in the paper**. The correction also flags the
  weaker logical claim: the paper says the approaches *can be applied* to general windows — it does
  **not** say the variance correction was proposed *because* of them.

**`reports/CS_PHASE8_C1_CLOSEOUT.md`** — box widened from *"the surviving fractions this report
tabulates"* to name all four pre-C6 digits that read as live beneath it, each corrected in place:
*"median λ̂ (12.8 µm)"* ×3 → **14.7 µm**; *"1 to 63"* → **1 to 66**; *"27 µm (N3)"* → **28 µm**
(N4-occ's 25 µm was already correct).

```bash
python3 -c "import pandas as pd; d=pd.read_csv('results/phase3/null_destructiveness.csv'); n=d[d['null']=='N3_occ']; print(sorted(n.n_admissible_moves), n.median_displacement_um.median())"
# [1, 6, 9, 10, 12, 66]  28.302669...      (N4_occ: [1,3,3,5,7,13], median 24.686)
python3 -c "import pandas as pd,numpy as np; d=pd.read_csv('results/phase3/main_fits.csv'); d['s']=d.section.str.split('_').str[0]; p=d[(d.call=='tierA_p95')&(d.s.isin(['7259','7260','7001','7248','7352','7435']))&(d.stratum=='all')]; print(len(p), np.median(p.lam_naive), np.percentile(p.lam_naive,[25,75]), p.lam_railed.mean())"
# 315  14.7321271  [7. 50.]  0.6
```

**`SASP_Kernel_Master_Plan.md:284` (§7 Rank 2)** — the second of the two locations
`BIO_DELIVERABLE7_CLAIM_AUDIT.md` names for the SenNet figure; §4.1 was corrected and **this one was
missed, precisely because the §4.1 marker was scoped to "this paragraph"**. Corrected to **2,041
datasets as of April 2026**, re-verified by me:

```bash
curl -s 'https://api.biorxiv.org/details/biorxiv/10.64898/2026.02.06.704469'
# v1 (2026-02-10) AND v2 (2026-05-04) both: "As of April 2026, the portal hosts 2,041
# publicly available human and mouse datasets across 15 organs using 6 general assay types."
```

No retrievable version states 1,753 or January 2026. The Rank 2 entry also repeated the CIT-4
attribution of Farzad et al. to *Cell* (it is *Cell Press Blue* 1(4):100053); flagged in the same
box. **And in `references.bib`:** the older `% AUDIT` note instructing *"Keep the 'as of January
2026' qualifier"* is **marked superseded** — its premise, that v1 states 1,753, is refuted by the
retrieved v1 abstract, and it directly contradicted the CIT-3 note below it.

**`reports/BIO_DELIVERABLE6_DISCUSSION.md:72` (§A.1 item 2)** — *"agree at or below chance (Jaccard
0.93–1.22×…)"* is **a claim, not a digit**, so the file's *"⚠ PRE-C6 DIGITS"* banner — scoped to the
SF / amplitude / power vector — never covered it. Corrected with its own box. **The argument of
item 2 is explicitly preserved**: "the sender call is a choice" does not need independence, only
disagreement about which cells are senescent, which still holds at ~99 %.

## 5.3 Found by my own post-fix sweep — commit `4aab386`

Fixing `WRITING_PACK.md` alone would have left the same defects live in the reports it copied them
from. All three are now corrected at source:

- **`reports/CS_PHASE8_TORUS_VAR.md:295`** carried the same refuted **"never outside ±0.010 of
  nominal"**.
- **`reports/CS_PHASE8_D2_DENOISE.md:87`** carried **"three `denoise=False` runs agreed at Jaccard
  0.76–0.99"** — **two** errors: 0.99 is a Pearson r, and there are **two** `denoise=False` seeds,
  not three (`grep -rho "raw_seed[0-9]" results/phase8_d2/ | sort -u` → `raw_seed0`, `raw_seed1`;
  the only runmeta files are `..._sub20000_...` and `..._sub20000_seed1_...`). The finding is
  sharpened, not weakened: the `denoise=True` outlier's Jaccard 0.000 is measured against a
  two-seed floor of 0.76.
- **`reports/CS_PHASE7_C1.md:308`** — the origin of the stale md5. Its two banners did not cover it.
  Marked, with its byte-identity claim **left standing as the record of the check it documents**.

---

# 6. THE ONE PLACE THIS PASS NEARLY REPEATED THE DEFECT IT WAS FIXING

R3 exists because a hash was **read from a report instead of a file**. My first correction replaced
`000f3405…` with the on-disk `d44fac63…`. Minutes later, task 8.7 regenerated the figure:

```
md5sum figures/figure4.png   ->  d44fac63411d6c30a42c40894a287f17   (first observation)
md5sum figures/figure4.png   ->  fdead29871b61481e297951dfea75b3d   (minutes later)
```

Replacing one stale literal with another would have reproduced the defect on a slower fuse. **The
hash is therefore removed from `WRITING_PACK.md` §5.4 rather than replaced**, and both that site and
`CS_PHASE7_C1.md:308` now instruct the reader to establish byte-identity by running
`python3 code/check_figures_guard.py` against `figures/.committed_manifest.json` and quoting its
**dated output**, never a literal. The audit offered "drop the md5 **or** replace it"; under an
active regeneration only one of those is stable.

The same reasoning applies to the guard's pass/fail state: the `WRITING_PACK.md` figure-guard row
now says it **must be re-checked at drafting time, not quoted from the pack** — see §0.1.

---

# 7. SCOPE AND VERIFICATION

**Files changed by this pass** (documents and `references.bib` only):

```
Phase7_Minimal_Human_Replication (1).md   reports/CS_PHASE8_C1_CLOSEOUT.md
README.md                                 reports/CS_PHASE8_D2_DENOISE.md
SASP_Kernel_Master_Plan.md                reports/CS_PHASE8_M1_RERUN.md
references.bib                            reports/CS_PHASE8_TORUS_VAR.md
reports/BIO_DELIVERABLE6_DISCUSSION.md    reports/NOVELTY_ASSESSMENT.md
reports/BIO_PHASE2.md                     reports/PLAN_UPDATE_D12_D13.md
reports/BIO_PHASE3.md                     reports/PREREG_PHASE8.md
reports/COMPLETED_TASKS.md                reports/WRITING_PACK.md
reports/CS_PHASE7_C1.md
```

Verified with `git diff --name-only` per commit: **no file under `results/`, `figures/`, `code/`,
`genesets/` or `data/` was touched by any commit of this pass.** `data/raw_h1/`,
`results/phase9_h1/`, `data/processed_h1/` and `code/h1_*` were neither read nor written. No
analysis was re-run; no tag created or moved; nothing pushed. Markdown table rows were checked for
field-count integrity after every edit.

**Concurrency note.** Three other audits ran against this working tree throughout. `git diff` over
the pass window therefore shows files this pass did not touch (`code/h1_*`, `figures/figure2*`,
`figures/figure3*`, `results/phase9_h1/*`, `code/make_figure*`, `.githooks/`). Those belong to other
agents. The figure-guard failure in §0.1 and the figure4 regeneration in §6 are consequences of
that concurrency, not of this pass.
