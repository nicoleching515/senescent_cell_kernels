# BIO PHASE 2 — Cell Types, Anatomy, Sender Calls, Module Scores

**Biology collaborator · 2026-08-20 · SASP Spatial Response Kernel**
Deliverables D-A (cell types), D-B (anatomy), D-C (senders), D-D (module scores), plus Tier E3.

---

> ## ⚠ SUPERSEDED IN PART — 2026-08-27 (remediation pass). §4.3's INDEPENDENCE CLAIM IS FALSIFIED.
>
> **This is a historical phase report, kept as the record of what was found on 2026-08-20 on the two
> published sections with the pre-C6 25-gene Tier A. Its findings are not rewritten.** Until this
> banner was added, the file's only correction box was §4.2's (the 69 % → 79 % circularity
> denominator); the ⚠️ marks elsewhere are the authors' own emphasis on a finding, **not correction
> markers**, and §4.3 has been read as live.
>
> ### Superseded — do not quote
>
> | where | what is dead | frozen replacement, and its file |
> |---|---|---|
> | **§0 headline item 2**, the clause *"the sender callers do **not** agree with each other at better than chance"* | the conclusion | **They do.** `results/phase3/caller_coverage_gate_headline.csv`: on **the same two sections**, replacing the pre-C6 25-gene Tier A with the frozen C6 33-gene set moves the pooled agreement from **1.040× chance (z = 1.76, p = 0.078)** to **1.131× (z = 5.80, p = 6.5e-09)**; at 11 sections it is **1.212× (z = 21.92, p = 1.84e-106)**, band 0.751–2.198. |
> | **§4.3 heading** *"⚠️ The sender callers agree at or below chance"* | the heading itself | See above. Read the heading as *"the sender callers agree at or below chance **on the two published sections, under the pre-C6 sender set**"*. |
> | **§4.3 body**, *"**The three definitions are statistically independent — they are not noisy measurements of one latent senescent state.**"* | **STRUCK.** Forbidden-claims item **3** / **13**; `PREREG_PHASE8.md` §10.12 | The dependence is real and its *shape* is the finding. Replacement wording, `reports/SUBMISSION_PATCH_2026-08-29.md` §2.1: *"…four senescence callers overlap at **0.75–2.20× of chance across eleven sections, pooling to 1.21× (Mantel–Haenszel, z = 21.9)**; they are not independent, but the dependence spans **anti-concordance at 0.74× to concordance at 1.47×** with a different identifiable cause at each end, which is not what one latent state looks like."* |
> | **§4.3 body**, *"two of three pairs overlap *less* than random"* | the count | On the frozen C6 set over the **same two sections**, **none** of the three is below chance: `tierA vs senepy` **1.007** (was 0.935), `tierA vs cdkn1a` **1.300** (was 1.017), `senepy vs cdkn1a` **1.168** — `results/phase3/caller_coverage_gate.csv`, `n_sections == 2`. |
>
> ### What still stands
>
> - **The ~99 % disagreement itself is not disputed** — the callers select almost disjoint cell sets.
>   What is falsified is the inference from that to *independence*: a large disagreement and an
>   above-chance overlap are not in conflict, and the frozen numbers show both.
> - **§4.2's circularity correction (26 of 33 = 79 %) stands**, and so do §§1–3, §4.1, §4.4 and §5.
>   Neither audit contradicted them. They were also not independently re-derived in this pass —
>   absence of a correction here is not a verification.
> - **The cause of the move is a sender-set defect, not coverage.** The pre-C6 Tier A used here was
>   25 genes; the frozen C6 set is 33. The fix alone kills independence on these two sections.
>
> Derivations: `reports/AUDIT_CLAIMS_CITATIONS.md` L3; `reports/PREREG_PHASE8.md` §10.12;
> `reports/WRITING_PACK.md` §3 items 3, 13, 14. Commands: `reports/REMEDIATION_PASS.md`.


## 0. Headline

1. **Cell typing failed twice before it worked.** The full before/after is in §2; the correction is
   on the record, not silently overwritten. Final sham composition is 51.3% hepatocytes, 0%
   Unknown, median confidence 0.966 (first attempt: 22.9% hepatocytes, 24.2% Unknown, 0.241).
2. **DeepScence and SenePy both installed and ran.** Section 10 gets its two independent
   sender-calling methods after all. But see §4: ~~**the sender callers do not agree with each
   other at better than chance**~~ ***[SUPERSEDED 2026-08-27 — see the banner at the head of this
   file. On these same two sections the frozen C6 Tier A gives **1.131× chance, z = 5.80,
   p = 6.5e-09**; at 11 sections **1.212×, z = 21.92**
   (`results/phase3/caller_coverage_gate_headline.csv`). They agree above chance.]***, and
   DeepScence's own gene set is **79%** circular with our
   response modules (**corrected 2026-08-27** from 69%; see the correction box in §4.2).
3. **The two arms are not equally analysable.** Section 8 Test 3 **passes for every SBR cell type**
   (hepatocytes 9.6% Cdkn1a+) and **fails for sham** (hepatocytes 0.48%). The paper's estimates
   must come from the SBR arm.
4. **The strongest threat is new and quantitative:** sender scores carry almost no spatial
   structure (Moran's I 0.004–0.019) against a zonation field at 0.38–0.75. §5.

---

## 1. What was produced

Per sample, keyed on `cell_id`, in `/workspace/data/processed/`:

| File | Contents |
|---|---|
| `celltypes_{s}.csv` | `cell_id, leiden, cell_type, cell_type_confidence` |
| `anatomy_{s}.csv` | `cell_id, cell_type, zonation_score, compartment_label, dist_to_boundary_um, dist_to_portal_triad_um` |
| `senders_{s}.csv` | `cell_id, cell_type, cdkn1a_counts, cdkn1a_pos, tierA_score, senepy_score, sender_flag_p90/p95/p99` |
| `modules_{s}.csv` | `cell_id` + 7 Tier B module scores |
| `test3_prevalence_{s}.csv` | Section 8 Test 3 table, per cell type |
| `composition_all_sections.csv` | composition across all 6 annotated sections |
| `annotation_meta_{s}.json` | label set, dropped types + reasons, parameters, composition |

`genesets/E3_random_matched/` — Tier E third component: **500 expression-matched random sets per
response module** (7 modules), plus `_expression_bins.csv`. Matched on size and on mean
log-normalised expression in 20 quantile bins, drawn from the same bins, excluding the module's own
genes. Stored as gene lists rather than per-cell scores (3,500 × 237k floats would be 3.3 GB).

**Cell typing runs as a function over a sample directory**, per your note:
`python3 code/annotate_pipeline.py <sample_dir> ... [--relabel]`. `--relabel` reuses saved Leiden
labels and re-runs only the assignment step, which is what made the three correction rounds
affordable (~90 s instead of ~35 min per section). **6 of the sections on disk are annotated**;
7352 (sham 2wk) was still clustering at write time.

---

## 2. D-A — cell typing: two failures and the fix (ON THE RECORD)

### 2.1 Before / after, sham (GSM9295284)

| Cell type | v1 (broken) | v3 (final) |
|---|---|---|
| Hepatocytes | **22.9%** | **51.3%** |
| LSECs | 1.0% | 14.7% |
| Kupffer cells | 2.4% | 8.1% |
| Hepatic stellate cells | — | 7.1% |
| Low_quality | — | 7.9% |
| Biliary/ductular | 2.9% | 2.9% |
| vSMCs | 2.9% | 2.9% |
| T/NK cells | — | 1.7% |
| Portal endothelial | — | 1.4% |
| Central venous LSECs | — | 1.0% |
| B-cells | — | 0.4% |
| Mesothelial | 0.2% | 0.25% |
| Proliferating | 0.2% | 0.23% |
| **Erythroid cells** | **28.4%** | **not callable — dropped** |
| **Mast** | **6.0%** | **not callable — dropped** |
| **Lymphatic endothelial** | **8.7%** | **not callable — dropped** |
| **Unknown** | **24.2%** | **0%** |
| median `cell_type_confidence` | **0.241** | **0.966** |

### 2.2 What was actually wrong, in three rounds

**Round 1 — the z-score bug (your diagnosis was right).** v1 z-scored each cell type's *aggregate*
score across clusters. That normalises away magnitude, so a 2-gene set that is barely expressed
anywhere still produces z≈+1.3 in *some* cluster and wins the argmax. Erythroid (only `Slc4a1`,
`Klf1` on panel) and Mast (only `Kit`, `Ms4a2`, `Tpsb2`) took 34% of the section on no evidence.
Fixed by (a) z-scoring **each marker gene** across clusters then averaging, so no single gene can
carry a set; (b) a **marker-support gate** — fewer than 4 informative on-panel markers and the type
is removed from the label set entirely; (c) requiring detection ≥5% and a margin ≥0.2 over the
runner-up.

**Round 2 — non-specific markers.** Lymphatic endothelium passed the count gate (4 on-panel
markers) but is still not callable. Measured in GSM9295284: `Ccl21a` and `Mmrn1` mean expression
**0.00 in every cluster**; `Prox1` is expressed by hepatocytes (0.11–0.19) as highly as anywhere;
`Lyve1` and `Flt4` peak in the LSEC clusters. Your independent check put lymphatic at 13% — that
number is `Lyve1`/`Flt4` cross-reactivity with LSECs, not a lymphatic population. Dropped with the
measurement recorded. **This also resolves the LSEC/lymphatic inversion you flagged**: v3 gives
LSECs 14.7% and lymphatic 0%.

**Round 3 — cholangiocytes, and a real biological finding underneath.** v3 initially put **25.8% of
the SBR section into "Cholangiocytes"**. It is not an annotation error and it is not cholangiocytes.
Measured across every cluster of both sections: **`Krt7` ≤ 0.03 and `Cftr` ≤ 0.06 — the two
definitive bile-duct markers effectively do not detect on this panel** — and `Krt19` reaches only
0.22 in the best sham cluster. What the panel resolves is a `Sox9`+/`Epcam`+/`Pkhd1`+/`Spp1`+
compartment containing **both** true cholangiocytes **and** hepatocytes undergoing ductular
metaplasia. They are not separable here, so the label is now **`Biliary/ductular`**, named for the
compartment rather than overstating it. My first attempted fix (dropping the injury-inducible
`Pigr`/`Spp1`/`Mmp7`/`Cldn4`/`Cldn7`) did **not** work, because ductular-reaction cells genuinely
express `Sox9`/`Epcam`/`Pkhd1` too — which is the point.

**`Low_quality` is a label, not a guess.** Clusters whose median counts fall below 50% of the
section median are segmentation fragments (sham cl13: 141 counts vs section median 513, median
nucleus area 18 µm² vs 37 µm²). They are labelled `Low_quality` and must be **excluded from receiver
analyses**. One exception is coded: a low-count cluster with unambiguous identity keeps its label —
mesothelial cells are thin and flat and were being discarded at z=+4.3, margin=3.7.

### 2.3 Independent validation: the annotation reproduces the IFALD phenotype

Composition across the 6 annotated sections (%; `composition_all_sections.csv`):

| Section | n | Hepatocytes | Biliary/ductular | LSECs | Kupffer | Stellate | T/NK | Low_qual |
|---|---|---|---|---|---|---|---|---|
| 7361 sbr 2wk | 193,983 | 41.97 | 2.26 | 11.57 | 6.62 | 6.06 | 1.92 | 8.78 |
| 7450 sbr 10wk | 93,197 | 33.41 | 8.75 | 14.00 | 12.34 | 9.60 | 2.98 | 12.63 |
| 7259 sbr 26wk | 127,386 | 23.63 | **25.75** | 10.18 | 10.37 | 11.87 | 4.95 | 9.00 |
| 7239 sbr 52wk | 83,392 | 32.66 | 20.07 | 9.16 | 8.17 | 10.40 | 2.72 | 9.60 |
| 7435 sham 10wk | 172,218 | 51.63 | 2.58 | 7.22 | 7.30 | 6.76 | 1.57 | 12.15 |
| 7250 sham 26wk | 236,905 | 51.31 | 2.93 | 14.72 | 8.10 | 7.14 | 1.67 | 7.86 |

Nothing in the pipeline knows about surgery, timepoint or arm. Yet **sham is flat across timepoints**
(hepatocytes 51.6/51.3, biliary 2.58/2.93) while **SBR shows a monotone, dose-dependent injury
phenotype**: ductular reaction 2.3 → 8.8 → 25.8%, stellate expansion 6.1 → 11.9%, immune
infiltration T/NK 1.9 → 5.0%, hepatocyte loss 42.0 → 23.6%. That is textbook IFALD and it is the
strongest evidence the annotation is now right.

### 2.4 Types this panel cannot call — a publishable limitation

Erythroid (`Hba-a1`,`Hbb-bs`,`Alas2`,`Gypa` all off-panel), Mast (`Cpa3`,`Cma1`,`Mcpt4` off-panel),
Neutrophils (`S100a8`,`S100a9`,`Retnlg`,`Mpo`,`Elane` off-panel), lymphatic endothelium
(non-specific), and true cholangiocytes vs ductular metaplasia (`Krt7`/`Cftr` non-detecting).
Reported as absent rather than guessed at.

---

## 3. D-B — anatomical covariates

Zonation score from the Phase 1 data-derived marker sets, standardised on hepatocytes;
`compartment_label` by hepatocyte tertiles; `dist_to_boundary_um` from a 25 µm occupancy grid with
closing/fill/opening and a Euclidean distance transform; `dist_to_portal_triad_um` from bile-duct
(`Biliary/ductular`) foci found by DBSCAN (eps=30 µm, min=10) — cKDTree throughout, never an N×N
matrix (§18.1).

**Sham validates cleanly.** 213 portal-triad foci. Median distance-to-portal-triad rises monotonically
periportal → midzonal → pericentral: **180.9 → 249.2 → 302.5 µm**, and
`corr(zonation, dist_to_portal_triad) = +0.255` within hepatocytes — the biologically correct sign,
recovered from two independent constructions (a marker score and a spatial landmark). Biliary cells
sit 34.8 µm from their own foci, as they should.

**⚠️ The portal-triad landmark is INVALID in the SBR arm.** In 7259 the same procedure finds 725
foci at median distance 59.9 µm, and `corr(zonation, dist_to_portal_triad) = +0.000`. The ductular
reaction disperses biliary cells throughout the parenchyma, so "distance to nearest biliary focus"
stops meaning "distance to portal triad". **Do not use `dist_to_portal_triad_um` for SBR sections.**
Use the continuous zonation score, which Section 11 explicitly permits as the alternative.

Zonation architecture is itself degraded by injury: Moran's I of the zonation score is **0.752 in
sham vs 0.379 in SBR**. Zonation-matched decoys will therefore behave differently between arms, and
the matching quality must be reported per arm, not pooled.

---

## 4. D-C — sender calls, and a serious problem

### 4.1 Both Section 10 tools are available

DeepScence 1.0.0 and senepy 1.0.1 both installed from PyPI. Two forced, documented deviations for
DeepScence: it ships a **human** core set and hardcodes `CDKN1A`, so mouse symbols are renamed to
1:1 human orthologs via the MGI `HOM_MouseHumanSequence` report — **the same resource the paper's own
authors used** (their notebook reads it directly); 4,845/5,097 panel genes map. And `denoise=False`,
because DeepScence's DCA denoising step needs an obsolete TensorFlow stack that will not install.
DeepScence was still running at write time; `senders_*.csv` gains a `deepscence_score` column when
it lands. SenePy uses cell-type-matched **mouse liver** hubs: hepatocyte hub 1 (143 genes, 54
on-panel), LSEC hub 0 (205/62), Kupffer hub 0 (638/226), plus Lymphoid T and B hubs.

### 4.2 ⚠️ DeepScence's own gene set is circular with our response modules

CoreScence v2 at the occurrence ≥5 threshold DeepScence uses: 39 genes, **33** on our panel.
**26 of those 33 (79%) are members of at least one Tier B response module.**

> **Correction, 2026-08-27** — `reports/AUDIT_PHASE8_FACTCHECK.md` M1, applied in
> `reports/AUDIT_CORRECTIONS_APPLIED.md`. This paragraph originally read "39 genes, 35 on our
> panel. 24 of those 35 (69%)". **The denominator 35 is wrong and is reproducible under no
> mapping convention.** Of the 39 CoreScence genes, 31 are reachable on the 5,097-gene mouse
> panel through the pinned 1:1 MGI map (`logs/ds_smoke.log`: "on our ortholog-mapped panel: 31";
> `logs/caller2.log`: "17 up / 14 down on mouse panel"), 33 with the Title-case fallback the
> project applies for symbols the map has no row for (`CDKN2B`, `CXCL1`), and 6 —
> `CXCL8 HMGB1 IGFBP5 IGFBP7 MIF TNFRSF10C` — are not on the mouse panel at all. The numerator 24
> below is correct for the strict 31 (24/31 = 77%); the module table that follows lists exactly
> those 24 genes. Under the 33-gene convention — the one `code/run_phase3_n8.py::corescence_mouse`
> implements and the one behind `corescence_on_panel = 33` in the committed
> `results/phase3/n8_disjointness_*.csv` — it is **26 of 33 = 79%**, the two extra circular genes
> being `Cdkn2b` and `Cxcl1`. Derived by `code/corescence_circularity.py`. **This makes the
> circularity finding stronger, not weaker.**

(**Bold** = added by the 2026-08-27 correction: the two genes the pinned MGI map has no row for.)

| Tier B module | CoreScence genes it contains |
|---|---|
| secondary_senescence | Ccl2 Cdkn1a Cdkn2a **Cdkn2b** **Cxcl1** Gdf15 Hmgb2 Igfbp3 Il1a Il6 Lmnb1 Mdm2 Serpine1 Tgfb1 |
| downstream_arrest | Brca1 Bub1b Ccna2 Cdk1 Cdkn1a Cdkn2a Hells Hmgb2 Lmnb1 Tgfb1 |
| emt_ecm | Fas Fgf2 Igfbp2 Igfbp3 Il6 Jun Serpine1 Tgfb1 Vegfa |
| tnfa_nfkb_proximal | Cdkn1a **Cxcl1** Icam1 Il1a Il6 Jun Serpine1 Vegfa |
| il6_jak_stat3 | Fas Il6 Jun Stat1 Tgfb1 |
| interferon_response | Cdkn1a Fas Icam1 Il6 Stat1 |

The up-direction half is dominated by secreted SASP factors — `IL6`, `IGFBP3`, `SERPINE1`,
`IGFBP1`, `FGF2`, `CCL2`, `IGFBP2`, `GDF15`, `IGF1`, `TGFB1`, `CXCL8`, `CXCL1`, `IL1A`.
**Scoring senders with DeepScence and measuring a SASP-adjacent response in neighbours is exactly
the Section 0.3 failure mode, quantified.** Section 10 says "start here"; on this project
DeepScence must be a *comparison* method, and the circularity reported, not the primary caller.

Comparison of the three sender definitions against the Tier B union (on-panel genes):

| Sender definition | on-panel | overlaps Tier B | % |
|---|---|---|---|
| DeepScence CoreScence (occ≥5) | 33 | 26 | **79%** |
| SenePy liver hepatocyte hub 1 | 54 | 13 | **24%** |
| Tier A union-strict | 25 | 0 | **0%** (by construction) |

SenePy is much the cleaner of the two published tools.

### 4.3 ⚠️ The sender callers agree at or below chance ***[on the two published sections, under the pre-C6 25-gene Tier A — SUPERSEDED 2026-08-27, see the head banner and the box below]***

Sham, within-cell-type p95 flags, and Spearman on the continuous scores:

| Pair | observed Jaccard | chance Jaccard | ratio | Spearman ρ |
|---|---|---|---|---|
| TierA vs SenePy | 0.0134 | 0.0223 | **0.60×** | −0.024 |
| TierA vs Cdkn1a+ | 0.0053 | 0.0060 | **0.88×** | +0.003 |
| SenePy vs Cdkn1a+ | 0.0097 | 0.0058 | 1.66× | +0.012 |

All three |ρ| < 0.03. Section 10 tolerates ~30% disagreement and asks for the analysis under both
callers if it is exceeded; here disagreement is **~99%** and ~~two of three pairs overlap *less* than
random~~. ~~**The three definitions are statistically independent — they are not noisy measurements of
one latent senescent state.**~~ Whatever λ we report will be a property of the sender definition we
chose. That has to be in the abstract, not a footnote.

> **Correction, 2026-08-27 (remediation pass)** — see the banner at the head of this file.
> **The independence sentence is STRUCK** (forbidden-claims items 3 and 13; `PREREG_PHASE8.md`
> §10.12), and so is "two of three pairs overlap less than random".
>
> The table above is the **pre-C6 25-gene Tier A on the two published sections**, and both defects
> matter. Re-run on the frozen C6 33-gene Tier A over **the same two sections**
> (`results/phase3/caller_coverage_gate.csv`, `n_sections == 2`, `tierA_definition ==
> "A_SENDER_FINAL_strict C6"`), **none of the three pairs is below chance**: `tierA vs senepy`
> **1.007** (was 0.935), `tierA vs cdkn1a` **1.300** (was 1.017), `senepy vs cdkn1a` **1.168**
> (unchanged — it does not involve Tier A). Pooled over the four gated pairs the two sections move
> from **1.040× chance (z = 1.76, p = 0.078)** to **1.131× (z = 5.80, p = 6.5e-09)**, and at 11
> sections to **1.212× (z = 21.92, p = 1.84e-106)**, band 0.751–2.198
> (`results/phase3/caller_coverage_gate_headline.csv`).
>
> **What survives, and it is stronger than what it replaces:** the ~99 % disagreement is real, but
> it does not imply independence — the callers overlap **above** chance while selecting almost
> disjoint sets. The wording to use is `reports/SUBMISSION_PATCH_2026-08-29.md` §2.1. The closing
> sentence of this subsection — that λ is a property of the sender definition chosen — is
> **unaffected and still stands**; if anything the frozen numbers sharpen it, since changing the
> Tier A definition alone moved every Tier-A-involving pair.

### 4.4 Section 8 Test 3 — passes for SBR, fails for sham

`Cdkn1a`+ prevalence per cell type (the percentile flags are 5% *by construction* and carry no
prevalence information — they are within-type rankings):

| Cell type | sham | SBR |
|---|---|---|
| Hepatocytes | **0.48% LOW** | **9.64% YES** |
| LSECs | 1.23% YES | 5.50% YES |
| Kupffer cells | 0.75% LOW | 2.93% YES |
| Hepatic stellate cells | 0.34% LOW | 2.31% YES |
| Biliary/ductular | 1.11% YES | 1.69% YES |
| T/NK cells | 1.14% YES | 2.03% YES |
| B-cells | 0.76% LOW | 1.94% YES |
| Proliferating | 17.09% YES | 19.29% YES |

**Every SBR cell type is in the 1–20% band; hepatocytes at 9.6% sit in Section 8's 2–10% sweet spot.
In sham, the main receiver population is at 0.48%, below the 1% floor.** The headline estimates must
come from the SBR arm; sham is a comparison arm, not an estimation arm.

Note `Proliferating` at 17–19% — `Cdkn1a` is induced in cycling cells, so it is not senescence-specific.
Recommend excluding `Proliferating` from sender calls, or reporting it separately.

---

## 5. ⚠️ THE BIGGEST THREAT: senders carry almost no spatial signal

k-NN (k=20) Moran's-I-style spatial autocorrelation, computed with cKDTree:

| Variable | sham | SBR |
|---|---|---|
| `dist_to_portal_triad_um` | 0.995 | — |
| **zonation score** | **0.752** | **0.379** |
| emt_ecm module | 0.161 | — |
| secondary_senescence module | 0.091 | — |
| interferon_response module | 0.089 | — |
| oxidative_stress module | 0.078 | — |
| tnfa_nfkb_proximal module | 0.051 | — |
| **`Cdkn1a` counts** | **0.0085** | **0.0194** |
| **Tier A sender score** | **0.0040** | **0.0124** |
| permutation null | −0.0006 | +0.0006 |

The sender signal is real — 20–30× the permutation null in SBR — but it is **one to two orders of
magnitude weaker than the zonation field it has to be separated from**. Two consequences:

1. Senders are only very weakly clustered, so distance-to-nearest-sender will have limited dynamic
   range, and Test 4 (Ripley's K, yours) should be run expecting a small effect.
2. Combined with your synthetic result — +45% bias in λ even in the easy regime, and matched decoys
   closing only 16–49% of the confounder gap while passing SMD<0.1 — the honest prior is that a
   naive λ here will be **mostly** density and zonation. Which is the paper's thesis, but it means
   the headline may have to be a negative result with a bound, not a positive λ.

I am **not** saying the estimate is impossible. I am saying: the go/no-go on Figure 2 should be
whether the kernel survives the null battery in the SBR arm, and that should be checked before
anyone writes an abstract claiming a length constant.

---

## 6. What the CS lead must know

1. **Use the SBR arm for estimation** (Test 3 passes); sham is the contrast arm.
2. **Exclude `Low_quality` (8–13% per section) and `Unknown` from receiver analyses.** They are in
   the CSVs so the exclusion is auditable rather than silent.
3. **`dist_to_portal_triad_um` is invalid for SBR sections** (§3). Use `zonation_score`.
4. **Report λ under at least two sender definitions and expect them to disagree** (§4.3). Given
   §4.2, do not use DeepScence as the primary caller with B1/B2/B3/B5/B7 readouts.
5. **Zonation strength differs by arm** (Moran I 0.75 sham vs 0.38 SBR) — report decoy-matching
   balance per arm, and per your own finding, report it as a **lower bound** on residual confounding.
6. `Proliferating` cells are `Cdkn1a`-high for cell-cycle reasons; handle separately.
7. Tier E3 nulls are gene lists, not scores — score them at fit time.
8. Per-module sender sets for single-module fits, union-strict 25-gene set for the §6.4
   proximal-vs-downstream comparison, exactly as you specified.

### Open / not done
- 7352 (sham 2wk) still clustering; 4 sections not yet on disk.
- DeepScence scores still computing; `deepscence_score` column lands when they finish.
- Anatomy/senders/modules currently generated for the two 26wk sections only; the other sections
  have cell types and need one `phase2_downstream.py` pass each.
- Transcript assignment rate still owed by you.

### Correction to my Phase 1 report
Phase 1 said the "no secreted factors" Tier A constraint was enforced *mechanically* via the panel's
`location` column. That is weaker than it sounded and I should flag it: the column only covers the
5,006-gene 5K panel, not the 100 custom genes (`Il6`, `Igfbp7`, `Igfbp5` have no annotation at all),
and it is single-valued and idiosyncratic — `Il1a` is annotated `Nucleus`, `Vegfa` `Cell membrane`,
`Icam1` `Membrane`. Tier A is still clean **by inspection** (all 25 members are intracellular), but
the automated filter should not be described as having verified that.
