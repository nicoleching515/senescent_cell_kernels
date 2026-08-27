# Submission patch — caller independence, for the Aug 29 deadline

> # ⚠ REFRAMED 2026-08-27 — THERE IS NO MANUSCRIPT
>
> This file was written as **corrections to a draft**. **No draft exists.** Nothing
> has been written yet, so there is nothing to find-and-replace.
>
> **What it is actually good for:** the numbers, the falsified claims, and the
> wordings that must not be used. All of that is being consolidated into
> **`reports/WRITING_PACK.md`**, organised by the §30 outline and traced to source
> files. **Use the writing pack when drafting. Keep this file only as the record
> of what was found wrong and when.**
>
> Note that two of its own sections were themselves corrected after being written
> (§4a's A7 attribution, §9's Moran framing) — it went stale twice while the
> numbers settled. That is the reason for the consolidation.


**Status: APPLY. Rewritten 2026-08-27 against the settled post-C6 numbers.**
This supersedes the "DO NOT APPLY — NUMBERS SUPERSEDED" hold of 06:09. Task 8.7
(the M1 end-to-end re-run) is **complete**, all stages landed, and every number
in §1–§3 below was read from `results/`, not from any report's prose.

**Decision D11 (PI-decision series): patch the claim, submit Aug 29 as planned.**
The manuscript is not in this repo, so this cannot be applied automatically.
Everything needed to apply it by hand is below.

---

## 0. Read this first — five things a skim will get wrong

**0.1 The independence claim dies from the sender-set fix, not from adding
sections.** This is the single most important change from the previous draft of
this patch, which framed the finding as "more sections revealed dependence".
That framing is wrong and it is the *weaker* half of the truth. On the
**published two-section base**, holding coverage fixed at exactly the two
sections the paper reports, replacing the published 25-gene Tier A with the
frozen 33-gene one moves agreement from **1.030× chance (p = 0.20 — consistent
with independence)** to **1.128× chance (p = 4.4 × 10⁻⁸)**. Independence is
already dead at n = 2. Adding the other nine sections then takes it to 1.212×.
Say it in that order: the published sentence was wrong on the published data.

**0.2 The "technical not latent" argument does NOT survive, and must not be
shipped in its published form.** The old restatement turned on *"one of the four
pairs sits below chance at 0.91× in all eleven sections"* and on *"the direction
of each pair is predicted by its transcript-depth loading."* Post-C6 that pair
is **0.972, z = −1.63, p = 0.10**, and the depth-direction rule is **refuted at
the level of independent units** (§2.3). A weaker, defensible claim is derived in
§2 and it is the one to use. **No replacement plank was manufactured**; one
genuine, pre-existing plank (a fifth pair that *is* reliably below chance) was
promoted, with its caveat attached.

**0.3 "Circular Tier A" is the wrong phrase and must not reach the manuscript.**
The defect being fixed was in the **response** module, not the sender set: the
pre-C6 B7 `secondary_senescence` had absorbed 27 canonical arrest genes, and
Tier A — which is defined by *subtraction* of Tier B collisions — was left as a
25-gene remnant containing no `Cdkn1a`, `Cdkn2a`, `Trp53`, `Lmnb1` or `Mki67`,
i.e. "a numerically passing but biologically hollow sender score"
(`genesets/README.md`). Re-sourcing B7 to a genuine secretome set returned eight
of those genes to Tier A. Write "**contaminated sender/response split**" or
"**hollow sender set**", never "circular sender set" — and see §3.4, because C6
moved one circularity measure in the *other* direction.

**0.4 The N2-vs-N5 result is now supported twice, independently, and it is the
project's strongest contribution.** An external novelty review ranked it #1 of
six findings. The composition-matched protocol (`CS_PHASE8_COMPMATCH.md`) then
found the same effect from the opposite direction: matched decoys remove
**1.6 %** of the naive amplitude where **the same variables entered as
covariates** remove **85.4 %** — a factor of fifty, on the same fits. §5.

**0.5 Contribution 3 stands, but its framing changes.** "We discovered that torus
shifts break on non-convex tissue" is **Lotwick & Silverman (1982)** and must be
struck. The claim is: a 40-year-old documented limitation is being violated in
current spatial-omics practice, and we quantify what it costs. §6.

---

## 1. Find these strings in the draft

Every string below is a property of the superseded two-section base, the
superseded 25-gene Tier A, or a claim the files no longer support.

| # | Search for | Why it must go | Replace with |
|---|---|---|---|
| S1 | `0.93–1.22` / `0.93-1.22` / `0.93 to 1.22` | Min/max of six two-section ratios under the hollow Tier A. Under the frozen Tier A on **the same two sections** it is 0.982–1.442; at eleven sections **0.751–2.198** | `0.75–2.20×` (11 sections) |
| S2 | `statistically independent` | Falsified. p = 1.8 × 10⁻¹⁰⁶ at 11 sections; p = 6.5 × 10⁻⁹ on the published two | delete; use §2 |
| S3 | `Four of six pairs` | Two-section artefact, and the pooled band is over **four** headline pairs, not six | see §2 |
| S4 | `the one pair that looked concordant in sham is anti-concordant in SBR` | Two-section artefact. SenePy vs DeepScence is 0.33–0.55× in **ten of eleven** sections; the 2.15× is section 7250 alone. There is no arm effect | §2, sentence 4 |
| S5 | `DeepScence's correlation with sequencing depth reverses sign between two sections` | Two-section artefact. ρ with transcript counts is **+0.29 to +0.56 in ten of eleven** sections and −0.350 in 7250 only | "…is inverted in one of eleven sections (7250), where the `CDKN1A` anchor is weakest" |
| S6 | `1.51–2.85` / `2.85×` / `1.51×` | The circular DeepScence–`Cdkn1a`⁺ pair. Over eleven sections its median is **1.071** and its pooled value **1.255**. Both published values are the two largest of the eleven — "1.51–2.85×" overstates the measured circularity ~2× | "median 1.071× over eleven sections (still excluded as circular by construction)" |
| S7 | `at chance` / `no better than chance` / `indistinguishable from chance` | Any variant of the independence phrasing | §2 |
| S8 | `2.3–2.5×` *(highest transcript-count quintile)* and `8–11×` *(lowest)* | Two-section values. At eleven sections: enriched **1.6–3.2×** in Q5, depleted **2.4–14.5×** in Q1 | as stated |
| S9 | `127k–237k cells` / `127,000–237,000` | Two-section cell range | **45,000–218,000 pairwise-complete cells in each of eleven sections** |
| S10 | `the arrest score and DeepScence run the other way` | **Keep the claim, check the wording** — it survives and is now stronger (11/11 sections, §2.2). Make sure it says *within cell type*; the *global* depth correlation of DeepScence is **positive** | "…run the other way **within cell type**, in every one of the eleven sections" |
| S11 | `Zhao et al.` | Wrong first author | **Ma S, Ji Z, Zhang B, … Liu G-H (2024), Cell 187(24):7025–7044**, `ma2024spatial` |
| S12 | `CellWHISPER` *(anywhere near "torus", "toroidal", "shift null")* | CellWHISPER's null is a within-cell-type location permutation — this project's **N1**, not N3. §7 | Lotwick & Silverman (1982); Mrkvička et al. (2021) |
| S13 | `nobody in this literature reports it` *(negative-control probes)* | False. §8 | the narrower claim in §8 |
| S14 | `we discovered` / `a finding in its own right` *(torus / non-convex tissue)* | 1982 prior art. §6 | the import-and-quantify framing in §6 |
| S15 | `in negative-control probes` *(attached to −0.070 or −0.074)* | Wrong response family. §4a | "in pooled negative-control features" |

**Also check:** the abstract, the introduction, the Methods caller table, §29
objection 2, and **any figure caption quoting an agreement ratio**
(`figure_phase8_callers` carries the current values; `fig_phase3_caller_depth`,
`fig_phase3_composition` and `fig_phase3_tierC_identifiability` are still drawn
from the *committed two-section, pre-C6* caller tables — see §7.4).

---

## 2. Replacement text

### 2.1 The paragraph (replaces `BIO_PHASE3.md` §4.5 and its manuscript copy)

> Four published or panel-standard ways of calling a senescent cell —
> DeepScence, SenePy, a disjoint arrest-and-damage score, and `Cdkn1a`
> positivity — were applied to the same 45,000–218,000 pairwise-complete cells
> in each of eleven Xenium Prime 5K mouse liver sections. After conditioning on
> cell type and on within-cell-type sequencing depth, their top-5 % calls overlap
> at **0.75–2.20× of chance**, pooling to **1.21× chance** (Mantel–Haenszel over
> eleven sections, z = 21.9). **They are not independent.** The dependence is
> weak in effect size — the largest agreement between two non-circular callers is
> 1.47× chance, far short of what two noisy measurements of one latent state
> would show — but its more informative property is that **it is not one
> mechanism**. The four callers split cleanly into two camps by which end of the
> within-cell-type depth distribution they select: SenePy and `Cdkn1a`⁺ are
> enriched 1.6–3.2× in the highest transcript-count quintile and depleted
> 2.4–14.5× in the lowest, while the arrest score and DeepScence run the other
> way, in every one of the eleven sections. The two callers that select opposite
> ends of that distribution most strongly, SenePy and DeepScence, overlap at
> **0.74× of chance** (z = −15.1; below chance in ten of eleven sections, and
> significantly so in ten), and re-anchoring DeepScence on a caller-free
> proliferation set strengthens this to 0.50×. The largest agreement, by
> contrast, is between the arrest score and `Cdkn1a`⁺ (1.47×), and it has a named
> biological cause rather than a technical one: four of the arrest score's genes
> (`Atm`, `Mdm2`, `Trp53bp1`, `Bax`) sit on the p53 axis that induces `Cdkn1a`,
> which is itself excluded from the score by construction. What each caller
> selects is separately identifiable: SenePy's cross-cell-type score scales with
> the number of its hub genes on the panel and so is not comparable across cell
> types; and DeepScence's polarity is fixed by a `CDKN1A` anchor that, once
> sequencing depth is partialled out, is weak or reversed in four of eleven
> sections. A single latent senescent state would produce uniform positive
> dependence among all four callers. What is measured instead spans 0.74× to
> 1.47×, with a different, identifiable reason at each end. A senescence call on
> targeted spatial data is therefore not a noisy measurement of one latent
> state — it is a choice of which end of the detection-depth distribution to name
> senescent, and the length constant that follows inherits that choice.

### 2.2 One-sentence core, if space is tight

> After conditioning on cell type and within-cell-type sequencing depth, four
> senescence callers overlap at 0.75–2.20× of chance across eleven sections,
> pooling to 1.21× (Mantel–Haenszel, z = 21.9); they are not independent, but the
> dependence spans anti-concordance at 0.74× to concordance at 1.47× with a
> different identifiable cause at each end, which is not what one latent state
> looks like.

### 2.3 Why the old argument was re-derived and not renumbered — the audit trail

The published restatement rested on two planks. **Both were tested against the
post-C6 files and one failed outright.**

**Plank 1 — "one of the four pairs sits below chance at 0.91× in all eleven
sections." DEAD.** That pair (arrest score vs SenePy) is now **0.972, z = −1.63,
p = 0.104, above chance in 4 of 11 sections** — statistically indistinguishable
from independence, not below it. It moved because Tier A did.

**Plank 1′ — the replacement, which was already on disk and is not new
analysis.** A **fifth non-circular pair, SenePy vs DeepScence, is reliably below
chance**: pooled **0.737, z = −15.08, p = 2.3 × 10⁻⁵¹**, below chance in 10 of 11
sections and *significantly* below in 10. It was excluded from the published band
only because the published band was defined as the four pairs whose two-section
values produced "0.93–1.22", not for any question of validity. Two properties
make it usable where the dead plank was not:
* it is **byte-identical pre- and post-C6** (it does not involve Tier A), so it
  cannot be an artefact of the gene-set change; and
* re-anchoring DeepScence on the caller-free eight-gene proliferation set
  **strengthens** it to 0.495 (z = −28.9), so it is not an artefact of the
  published `CDKN1A` anchor either.

  **The caveat that must travel with it, every time:** ranking by |score| rather
  than by signed score puts the same pair **at chance** (1.025, z = 1.4, n.s.).
  The anti-concordance is about **polarity**, not about which cells are extreme.
  State this in the same breath or a reviewer will state it for you.

**Plank 2 — "the direction of each pair is predicted by its transcript-depth
loading." REFUTED as a general rule; retained only as the two-camp observation.**

The two-camp fact itself is solid and post-C6:

| caller | within-type Q5/Q1 depth enrichment, 11 sections | direction | sections |
|---|---|---|---|
| SenePy | 10.58 – 41.74 | top-selecting | **11/11** |
| `Cdkn1a`⁺ | 4.19 – 42.36 | top-selecting | **11/11** |
| arrest score (Tier A) | 0.146 – 0.317 | bottom-selecting | **11/11** |
| DeepScence | 0.218 – 0.795 | bottom-selecting | **11/11** |

*(`results/phase3/caller_within_type_depth_bias_11sections.csv`, Q5/Q1 of the
enrichment column, pivoted per section × caller.)*

But the rule "same camp ⇒ above chance, opposite camps ⇒ below chance" does not
hold, and the exception is the largest value in the matrix:

| pair | camps | pooled ratio | rule predicts |
|---|---|---|---|
| arrest vs `Cdkn1a`⁺ | **opposite** | **1.471** (11/11 above) | below chance — **WRONG** |
| arrest vs DeepScence | same (bottom) | 1.288 (11/11 above) | above — right |
| SenePy vs `Cdkn1a`⁺ | same (top) | 1.211 (9/11 above) | above — right |
| arrest vs SenePy | opposite | 0.972 | below — borderline |
| SenePy vs DeepScence | opposite | 0.737 | below — right |

Three tests, all against the rule:
1. **Pair-level exact permutation** over the five independent non-circular pairs
   (2 same-camp vs 3 opposite): **p = 0.30**. The observed assignment is not even
   the best of the ten possible ones.
2. **Within-pair, across the eleven sections**, the product of the two callers'
   log Q5/Q1 loadings is **negatively** rank-correlated with the agreement ratio
   in **all five** pairs (Spearman ρ = −0.16 to −0.70) — the opposite of the
   predicted direction. The one nominally significant case (arrest vs `Cdkn1a`⁺,
   ρ = −0.700, p = 0.016) points the wrong way.
3. The continuous version pooled over all 55 pair-sections is null
   (Spearman ρ = +0.096, p = 0.49).

A Mann–Whitney over the 55 pair-section values *does* separate same-camp
(median 1.258) from opposite-camp (median 0.990) at p = 0.008 — **do not quote
it.** It is pseudo-replicated: 55 values from five independent pairs.

**So the claim that must be dropped is "the dependence is weak and technical
rather than latent, and the direction of each pair is predicted by its depth
loading."** The claim the data supports is the one drafted in §2.1: the
dependence is weak, heterogeneous, and mechanistically *different* at each end —
which is still incompatible with one latent state, and is the honest version.

---

## 3. Numbers to quote, with provenance

All from `results/phase3/caller_coverage_gate.csv` (per pair) and
`caller_coverage_gate_headline.csv` (the pooled band). **Both files now carry
explicitly-labelled bases** — every row names its Tier A definition and gene
count, so a pre-C6 value can no longer be silently pooled against a post-C6 one.
The headline file carries **six** bases; the per-pair file carries **four** (it
omits the two six-section in-band rows). Independently re-derived here from
`caller_agreement_matched_significance_{verify2sec,2sec_c6,11sections}.csv` by
Mantel–Haenszel pooling: **every value below reproduces exactly.**

### 3.1 The decomposition — quote this, and quote it in this order

The two effects are different and the paper must not conflate them.

| step | what changes | pooled ratio (4 headline pairs) | z | p |
|---|---|---|---|---|
| published state | 2 sections, pre-C6 25-gene Tier A | **1.040** | 1.76 | **0.078** |
| **① fix the contaminated sender/response split** *(coverage held at 2 sections)* | 2 sections, **frozen 33-gene Tier A** | **1.131** | 5.80 | **6.5 × 10⁻⁹** |
| **② add coverage** *(Tier A held at the frozen set)* | **11 sections**, frozen Tier A | **1.212** | 21.92 | **1.8 × 10⁻¹⁰⁶** |

On the **three-pair basis that literally defines the published "0.93–1.22×"
band** (arrest vs SenePy, arrest vs DeepScence, arrest vs `Cdkn1a`⁺), the same
decomposition is: **1.030 (p = 0.20) → 1.128 (p = 4.4 × 10⁻⁸) → 1.212
(p = 1.8 × 10⁻⁹⁴).** This is the version that maps one-to-one onto the sentence
being struck, and it is the sharper one: **step ① alone kills independence, on
the published data.**

For completeness, the same two steps under the *pre-C6* Tier A give
1.040 → 1.129 (coverage effect at the old sender set), so coverage matters under
either definition — but it is not what breaks the claim.

### 3.2 Band and pooled values, all six bases

`caller_coverage_gate_headline.csv`, four headline pairs
(arrest × {SenePy, DeepScence, `Cdkn1a`⁺}, SenePy × `Cdkn1a`⁺):

| basis | Tier A | band | median | pooled | z | p | above chance |
|---|---|---|---|---|---|---|---|
| 2-section (PUBLISHED) | pre-C6, 25 | 0.932–1.369 | 1.010 | 1.040 | 1.76 | 0.078 | 4/8 |
| 11-section (task 8.4) | pre-C6, 25 | 0.700–1.711 | 1.156 | 1.129 | 13.35 | 1.1 × 10⁻⁴⁰ | 29/44 |
| **2-section (FROZEN)** | **C6, 33** | **0.979–1.442** | **1.099** | **1.131** | **5.80** | **6.5 × 10⁻⁹** | 6/8 |
| **11-section (FROZEN) — USE THIS** | **C6, 33** | **0.751–2.198** | **1.190** | **1.212** | **21.92** | **1.8 × 10⁻¹⁰⁶** | **35/44** |
| 6-section in-band | pre-C6, 25 | 0.775–1.374 | 1.131 | 1.115 | 8.99 | 2.6 × 10⁻¹⁹ | 16/24 |
| **6-section in-band (FROZEN)** | **C6, 33** | **0.811–1.565** | **1.168** | **1.167** | **13.09** | **3.6 × 10⁻³⁹** | 18/24 |

Sign test on the 11-section frozen row: **p = 1.06 × 10⁻⁴**. The last row is the
one that answers "is this driven by the sections Test 3 excludes?" — **no**:
restricting to the six Test-3-admissible sections Phase 3 actually fits gives a
narrower band and the same conclusion at p = 3.6 × 10⁻³⁹.

### 3.3 Per pair, frozen Tier A, eleven sections

| pair | pooled | z | p | min–max | above chance | sig. above | sig. below |
|---|---|---|---|---|---|---|---|
| arrest vs `Cdkn1a`⁺ | **1.471** | 19.45 | 2.7 × 10⁻⁸⁴ | 1.085–2.198 | **11/11** | 9 | 0 |
| arrest vs DeepScence | **1.288** | 19.23 | 2.1 × 10⁻⁸² | 1.096–1.660 | **11/11** | 10 | 0 |
| SenePy vs `Cdkn1a`⁺ | 1.211 | 7.43 | 1.1 × 10⁻¹³ | 0.842–1.391 | 9/11 | 7 | 0 |
| arrest vs SenePy | **0.972** | −1.63 | **0.104** | 0.751–1.179 | 4/11 | 1 | 3 |
| **SenePy vs DeepScence** | **0.737** | **−15.08** | **2.3 × 10⁻⁵¹** | 0.332–2.150 | 1/11 | 1 | **10** |
| *(circular, never pooled)* DeepScence vs `Cdkn1a`⁺ | 1.255 | 10.53 | 6.2 × 10⁻²⁶ | 0.963–2.849 | 7/11 | 5 | 0 |

**The three pairs that do not involve the arrest score are byte-identical pre-
and post-C6.** Only the three `tierA_*` pairs moved, and all three moved in the
same direction the frozen sets predict. On the published two-section base the
frozen sets give arrest vs `Cdkn1a`⁺ **1.017 → 1.300**, arrest vs DeepScence
**1.103 → 1.179**, arrest vs SenePy **0.935 → 1.007**.

### 3.4 Two disclosures that belong in the same paragraph as §3.3

**(a) Why arrest vs `Cdkn1a`⁺ jumped, stated before a reviewer asks.** The eight
genes that re-entered Tier A are `Atm Bax Bcl2l1 Cdkn2b Glb1 Mdm2 Sirt1
Trp53bp1`; four sit on the p53 axis that induces `Cdkn1a`. `Cdkn1a` itself is
**still not in Tier A** — disjointness is re-verified and passes — so this is
biological correlation, not shared membership. It is nonetheless a real cost of
the C6 decision and the independence framing pays it.

**(b) C6 moved one circularity measure the wrong way, and the affected pair is a
headline one.** CoreScence circularity (the share of DeepScence's own anchor gene
set falling inside Tier A ∪ Tier B) rose from **79 % to 88 %** on the mouse arm
(`code/corescence_circularity.py`; the previously-quoted "69 %" was a typed-in
literal with an irreproducible denominator). Direct Tier A ∩ CoreScence
membership went from **2 of 25** (`Foxm1`, `Parp1`) to **4 of 33** (adding
`Cdkn2b`, `Mdm2`). DeepScence's *sign anchor* is `CDKN1A`, which remains outside
Tier A, so the arrest-vs-DeepScence pair (1.288) is not circular by
construction — but two shared genes out of 39 CoreScence members is a real
channel and should be stated rather than found.

---

## 4. Two other claims to check while you are in there

### 4a. The A7 / decoy claim — CORRECTED, and re-run under the C6 senders

**An earlier version of this file said "a −0.070 SD gradient (p = 0.023) *in
negative-control probes*". That attribution is WRONG and must not ship.** It was
caught by independent audit (`AUDIT_PHASE8_FACTCHECK.md` R1) and is preserved
here as a correction, now updated to the post-C6 values. Note that
`CORRECTIONS.md` §12 still carries the old phrasing ("the negative-control-probe
kernel is −0.074 SD") — **that is the same error and must not be copied.**

−0.074 is the **pooled `all_controls`** response — 40 negative-control probes
**plus** 609 negative-control codewords **plus** 21 genomic controls, of which
the codewords are ~73 % of the counts. Per response, naive design, post-C6
(`results/phase3/a7_summary.csv`, `design == 'base'`):

| Response | Gradient | 95 % CI | p | Flat? |
|---|---|---|---|---|
| `all_controls` (pooled) | **−0.0744** | [−0.1306, −0.0182] | **0.0145** | no |
| `neg_control_codeword` (609) | −0.0604 | [−0.1085, −0.0123] | 0.0188 | no |
| `genomic_control` (21) | −0.0307 | [−0.0558, −0.0056] | 0.0213 | no |
| **`neg_control_probe` (the 40)** | **−0.0225** | **[−0.0527, +0.0078]** | **0.129** | **YES — flat** |
| `neg_probe_rate` (ratio) | +0.0113 | [−0.0085, +0.0310] | 0.232 | yes |

**Why this matters more than a naming slip:** `PREREG_PHASE8_genesets.md` §11
designates those **40 negative-control probes the *primary* A7 null**, and Phase 9
item 9.4 repeats it. **On the pre-registered response, M1's A7 passes naively.**
Writing "the assay is not flat in negative-control probes" states the opposite of
the pre-registered test's own result.

**What survives, stated accurately:** the assay carries a distance gradient in
the **codewords and genomic controls**, and in the pooled control set — but *not*
in the 40 named negative-control probes, which are the sparsest of the three
(0.0067 counts/cell, 0.65 % of cells non-zero) and correspondingly the least
powered. Say that, and say which response each number comes from. **Do not let
"controls" stand unqualified anywhere.**

**The N2-vs-N5 result is unaffected and is the finding to carry forward.** On
`all_controls`, post-C6: N5 removes the gradient (**+0.0038, p = 0.72**; with N6,
+0.0053, p = 0.60), N2 matched decoys do **not** (**−0.0642, p = 0.0124** —
86 % undiminished). Same pattern on codewords and genomic controls. It is
undefined on the 40 probes, because there is no gradient there to remove.
**Consequence stands: never report a naive or N2-only kernel on this platform.**

**One caveat disappeared under the C6 senders — recorded so its absence is not
silent.** `neg_probe_rate` under N5 was **+0.0111 [+0.0025, +0.0197], p = 0.016**
pre-C6 — the one control response that looked non-flat *after* conditioning. It
is now **+0.0100 [−0.0056, +0.0255], p = 0.183**. The published caveat
("the one response of five whose conditioned amplitude is nominally non-zero")
is **no longer true and should be deleted, not softened.** This is a sender-set
effect; it removes a caveat rather than adding one.

**Measured false-positive rate: unchanged at 9–16 %.** Under the full N6+N5
design the two-sided 95 %-CI exclusion rate is 0.091 / 0.103 / 0.109 / 0.145 /
0.164 across the five control responses. Two things must travel with it:
* the **16 % upper end is `neg_probe_rate` alone**, the one response whose
  denominator is itself an N5 column and which is not a clean null. On the four
  count-based responses the range is **9.1–14.5 %**. Quote "9–16 % (9–15 % on the
  four count-based responses)".
* it is a **two-sided CI-exclusion rate for the estimator**, not the rate of the
  pipeline's "reportable fit" filter. That filter is one-sided on the naive
  design and measures **3.0–13.3 %** on the same null responses. The published
  sentence "the reportable-fit filter admits two to three times more fits than
  its nominal rate implies" is **not supported** and must be struck
  (`AUDIT_PHASE8_FACTCHECK.md` R6).
* A7 is powered **only pooled across sections** — a single fit resolves only
  ±0.134 SD, which is 1.7× the conditioned biological amplitude of 0.079. That
  caveat has to travel with the number every time.

**One more sentence to check while you are in §4a.** "The binned curve rises from
+0.019 SD in the 0–5 µm bin to +0.29 SD in the 95–100 µm bin. This gradient is a
quarter of the naive biological gradient (+0.291 SD)." The 0.019→0.29 span is a
**bin range**; the "quarter" is the ratio 0.074/0.277 of two **amplitudes**.
Putting them in adjacent sentences with near-identical numerals reads as though
the control curve spans a quarter of the biological one when it spans essentially
all of it. Separate the two sentences (`AUDIT_PHASE8_FACTCHECK.md` M4).

### 4b. Signature-source sensitivity

If the draft promises that conclusions do not depend on the choice of senescence
signature, that promise is now broken and the paper must say so — its own §29
wording already commits to this ("If they do, we say so"). SenePy also ships no
spleen signature at all, which matters for the human arm (22 spleen cell types:
0 matched, 15 surrogate, 7 none).

---

## 5. The N2-vs-N5 result — now supported twice, independently

**Give this the prominence it has earned.** An external novelty review ranked it
the project's **strongest contribution** of six ("highest surprise-per-word in
the repo"), on the grounds that it is unreported, mechanistically explained,
directly transferable to every imaging-ST distance analysis, and **overturns the
project's own prior** that a matched-decoy contrast is the conservative option.
`SASP_Kernel_Master_Plan.md` §23 calls the matched-decoy number "the single most
important number in the paper" — **that sentence has to change too.**

It now rests on two independent measurements pointing the same way:

| instrument | matched decoys remove | the same variables as covariates remove | ratio |
|---|---|---|---|
| **A7**, negative-control features (`a7_summary.csv`) | ~0 % of the −0.074 SD control gradient (N2: −0.064, p = 0.012) | **all of it** (N5: +0.004, p = 0.72) | — |
| **Composition-matched protocol** (`compmatch_reruns.csv`) | **1.6 %** of the naive amplitude (SF 0.9837 [0.973, 0.994]) | **85.4 %** (SF 0.1461 [0.052, 0.246]) | **≈ 50×** |

The composition arm is a **like-for-like** comparison and that is what makes it
evidence rather than an anecdote: the `comp` matching set stratifies exactly on
receiver cell type and matches on the 20-NN composition vector, and
`typecomp_adj` enters *those same variables* as covariates on *the same fits*.
Matching balances well — max |SMD| 0.092 → 0.035, median match rate 0.99987, the
§8 Test 5 gate (|SMD| ≤ 0.1) passes in **100 %** of matches — and removes almost
nothing.

**The mechanism, in one sentence, which is what makes it transferable:** matching
balances the covariates between senders and decoys; it does not remove the
dependence of the *response* on those covariates **at the receiver**, which is
where the confounding acts. On A7 the specific confounder is named: per-cell
detection efficiency, which N5 models directly and a propensity match on
neighbourhood covariates cannot see.

**Two consequences for the manuscript, both of which must be written down:**
1. **Never report a naive or N2-only kernel on this assay.**
2. The §17 "composition surrogate share **66–76 %**" row — which had no producer
   and is flagged as untraceable — is superseded. The honest statement is
   **66 % (receiver's own cell type) to 85 % (own cell type + 20-NN
   neighbourhood composition)**, from `compmatch_reruns.csv`, which now emits it
   four ways with its scope, estimator and reportable population stated.

---

## 6. Contribution 3 — it stands, and the framing changes

**The result is unchanged and is now better defended.** The variance-corrected
random shift, which is the null the classical literature prescribes for a
non-rectangular window, agrees with the project's own numbers:
**N3-var = 0.996 [IQR 0.975, 1.007]** against N3-tile 0.971 and a published
bounding-box 0.999; N4-var 0.985 against N4-tile 0.924 and a published 0.947.
The window-matched version of N3-var is 0.995, so the agreement is not an
artefact of comparing statistics computed on different amounts of data.

**What must be struck:** any sentence of the form "we discovered that torus
shifts break on non-convex tissue", or "this is a finding in its own right".
Toroidal shifts require a rectangular window — **Lotwick & Silverman (1982)**,
JRSS-B 44(3):406–413 — and the variance correction that replaces them on
irregular windows is **Mrkvička, Dvořák, González & Mateu (2021)**, *Spatial
Statistics* 42:100430. `spatstat`'s `rshift.ppp` documents the requirement.

**What to say instead:** a 40-year-old documented limitation is being violated in
current spatial-omics practice, and we quantify what it costs. Three things in
that sentence are genuinely ours:
1. **The cost, measured on real tissue.** Under a bounding-box wrap, **35.5 %**
   of N3's shifted senders land outside the tissue and 19.9 % of N4's
   (`1 − frac_in_occupancy`). *Do not quote the 23 % / 8 % figures as
   "out of tissue"* — those are a **different column**,
   `1 − frac_retaining_a_neighbour`, i.e. shifted senders left with no real cell
   inside the 100 µm window (22.8 % / 8.0 %). The two must not travel in one
   sentence (`AUDIT_PHASE8_FACTCHECK.md` R3). At a ≤5 % out-of-tissue tolerance
   only 1–66 of 38,080–108,375 candidate offsets are admissible, median
   displacement 28 µm against a pooled λ̂ of 15.7 µm, and one section admits only
   the identity.
2. **A direct calibration measurement of what the violation costs**, which no
   spatial-omics paper has. On an irregular synthetic window under the null of
   independence, the **tiled** torus rejects at **0.080–0.118 against a nominal
   0.05 — up to 2.4× nominal** — while the variance correction (RS_count) holds
   0.033–0.060 across every window and correlation scale. The whole-window torus
   is at 0.033–0.073. Reported against interest: **this means the project's own
   C1 correction replaced a liberal test with a more liberal one**, and it is why
   **N3-var, not N3-tile, should be presented as the primary corrected null.**
3. **The FFT enumeration.** "Which translations keep ≥ x % of a point set inside
   a mask?" is a circular cross-correlation, so one `rfft2` gives the exact
   admissible set over all offsets at once. Rejection sampling would never have
   revealed that the admissible set is one offset wide. Lead the methods
   paragraph with this — it converts a rediscovery into a tool.

**Two honest caveats to keep attached.** The calibration study is synthetic
(Gaussian fields, 100 sampling points) and establishes the *direction* of the
tiling effect, not a type-I error number for the Phase 3 fits — the project's
instrument for that is A7. And on the real data the tiled null looks slightly
*more* conservative (SF 0.971 vs 0.999), because at a 1,200 µm tile side and
λ̂ ≈ 15.7 µm the seams are ~76 λ̂ apart, so the affected fraction of cells is
small. Say both; a reviewer who knows Mrkvička §2.1.4 needs to see that we did.

---

## 7. Citation corrections

**7.1 CellWHISPER is not the source of the torus-shift null.** Its null is a
within-cell-type location permutation — **this project's N1, not N3** — and the
strings "torus", "toroidal" and "wraparound" do not occur in the paper. This was
load-bearing prose in six places and **N3 is what Figure 4 rests on**. The figure
is unaffected (the statistics were never theirs); the attribution sentence is
not. Correct attribution: Lotwick & Silverman (1982) / Mrkvička et al. (2021).

**7.2 Keep the ">90 % FPR" wording as an inference, not a measurement.**
CellWHISPER's argument is an interaction-count ratio between real and permuted
input from which they *infer* an FPR. Write "reported comparable interaction
counts on randomized coordinates, implying FPR > 90 %". Do not let it drift into
"has a measured false-positive rate of >90 %".

**7.3 Nineteen of 32 bibliography entries carried invented author forenames** —
41 wrong given names, now corrected in `references.bib` against both the Crossref
deposit and the PubMed record. **If the manuscript's bibliography was typed
rather than exported from `references.bib`, re-export it.** Includes
`martin2023modelling` (Luke → **Lucy**), `qu2025deepscence` (Yi → **Yilong**),
`ma2024spatial`, `kumar2026cellwhisper` (4 of 6 wrong), `ntintas2026overview`
(5 of 7 wrong).

**7.4 Citations the reviewers will expect and `references.bib` does not have.**
Adding them is cheap and each closes a named objection: Lotwick & Silverman
(1982) and Mrkvička et al. (2021) for the shift null; Lipsitch, Tchetgen Tchetgen
& Cohen (2010) for the negative-control-outcome construction; Hodges & Reich
(2010) and Dupont, Wood & Augustin (*Spatial+*, Biometrics 78(4):1279, 2022) for
spatial confounding, of which `references.bib` currently has **none**; and, for
the caller section, ICE (*Genome Biology* 2026), markeR (*NAR GAB* 2026) and the
bioRxiv 2026.01.02.697374 preprint whose "apparent concordance may reflect
circular validation" is the closest published statement to this project's own
thesis. **Do not claim novelty for "senescence callers disagree"** — DeepScence,
SenCID, SenePy, ICE, markeR and Ntintas et al. all report it.

---

## 8. What does NOT need to change

- **The central negative result stands, on both pre-registered Tier A
  definitions.** Post-C6: naive amplitude **0.329**, controlled (N2+N5+N6)
  **0.029** [−0.007, 0.084], SF **0.088** [−0.017, 0.234], detectable bound at
  80 % power **0.183**, 13 of 153 controlled fits positive with a CI excluding
  zero. The controlled amplitude remains far below the bound, so §18 outcome
  **A** stands. *(These are small updates from the published 0.326 / 0.027 /
  0.082 / 0.203 / 15 of 160 — check whether the draft quotes the old ones.)*
- **C1's verdict stands.** The corrected in-tissue nulls return the published
  value across the whole family, and the variance correction confirms it
  (§6). Invariant across six sender callers and prevalences from 0.5 % to 9 %.
- **The disjointness gate passes** on the frozen configuration, both arms,
  re-verified independently against the panel derived from the data (5,097 mouse
  genes). Tier A strict 33/33 on panel; all seven Tier B ≥ 30; A ∩ ∪B = 0.
  **But:** the seven Tier B modules are **not** mutually disjoint — 18 of 21
  module pairs share genes. No gate requires it; no text may claim it.
- **Figures.** Figure 1 is unchanged. Figures 2a, 2b, 2c, 2d, 2e, 3 and 4 were
  regenerated once, from the frozen configuration, at task 8.7, and **have
  moved** — if the draft embeds pre-8.7 versions, re-export them. The regeneration
  ledger's exemption of 2a and 2d ("not null-dependent") was valid for the null
  correction and **not** for a gene-set change; both are sender-dependent and both
  changed.
- **Three figures were on a stale basis and were re-pointed at 09:18 by a
  concurrent agent — verify before use.** `fig_phase3_caller_depth`,
  `fig_phase3_composition` and `fig_phase3_tierC_identifiability` read the
  *committed two-section, pre-C6* caller tables; they regenerated byte-identically,
  which looked like a passing reproducibility check and was a stale-input result.
  They now read `caller_*_2sec_c6.csv`, i.e. the **frozen** sets — but still on the
  **two-section** base, not the eleven-section one `figure_phase8_callers` carries.
  Depth-enrichment panels are fine on that base; **any caption quoting an agreement
  ratio is not**, and would contradict §3 of this patch. Check the three captions.
- **Do not run `python3 code/check_figures_guard.py --snapshot`.** The guard now
  walks all of `figures/` (46 artefacts, up from 27 via `git ls-files`) and its
  manifest is current at the post-8.7 state. `--snapshot` would re-baseline it
  and bless whatever is on disk.

---

## 9. Optional, if space allows — two free results

Both are already measured, transferable, and cost nothing to state.

- **Negative-control-outcome kernel (A7).** Refit the estimand's own estimator —
  the distance-to-nearest-sender kernel — with control counts as the response.

  > **Correction (citation audit).** "Nobody in this literature reports it" is
  > **false as written and must be struck.** The Voyager Xenium vignette computes
  > Moran's I on negative control probes and codewords — "generally the negative
  > controls are tightly clustered around 0, while the real genes have positive
  > Moran's I, which means there is generally no technical artifact spatial
  > trend" — and Ren *et al.*, *Nat Commun* 16 (2025),
  > doi:10.1038/s41467-025-64292-3, does the same peer-reviewed: "Spatial
  > autocorrelation analysis using Moran's I revealed stronger aggregation of
  > negative control signals in CosMx 6K." Both quotations retrieved and verified.
  > **What survives** is the narrower claim: nobody refits *the estimand's own
  > estimator* with control counts as the response. That is a negative control
  > *outcome* in the sense of Lipsitch *et al.* (2010), and it is a different test
  > from a global Moran's I — a signal with no global autocorrelation can still
  > project onto a specific covariate. **Report your own Moran's I on the controls
  > alongside the kernel amplitude**.
    >
  > **CORRECTED 2026-08-27 — DO NOT CLAIM THE TWO TESTS DISAGREE.** The claim that the two tests disagree is **FALSIFIED**. Moran's I has now been run (`reports/CS_PHASE8_MORAN.md`): across 12 control and module fields the two statistics rank **together**, Spearman rho **+0.923 raw / +0.944 cell-type-centred**, re-derived from `results/moran/moran_vs_a7.csv`. A reviewer would see it in a single plot. **The correct defence is POWER, not orthogonality:** the entire A7 gradient contributes **0.83%** of the observed control Moran's I, and the smallest amplitude Moran's I can resolve is **0.362 SD** — larger than the project's own naive biological amplitude (**0.291 SD**). **Moran's I could not have detected the headline effect either.** This gap is now CLOSED.

  **And say this in the same breath, or a reviewer will say it for you:**
  `neg_probe_rate` (probe counts ÷ transcript counts) is **flat naively**
  (+0.011, p = 0.23), so the gradient is a per-cell **detection-efficiency**
  effect projected onto distance-to-sender, not a spatial gradient in probe
  binding — which is exactly why an N2 match on neighbourhood covariates cannot
  catch it. Reporting "the raw assay is not flat" without that sentence invites
  the reading that you have rediscovered "bigger cells have more counts".

- **A measured false-positive rate of 9–16 % against a 5 % nominal**, obtained
  free from the same control features, with no randomization at all. Every
  competing FPR in this literature is obtained by randomizing real data; this one
  is measured *in situ* against the real spatial structure. Report it as one
  Methods sentence and one table row attached to the A7 result, **not** as an
  abstract claim, and carry the three caveats in §4a (which response, powered only
  pooled, and it is the estimator's rate rather than the filter's).

---

## 10. What changed in this document, and why

Recorded so the next reader does not have to diff it. This file has now been
corrected three times; two of those were caught by independent audit.

| # | Was | Now | Why |
|---|---|---|---|
| 1 | §1–§3 computed under the pre-C6 25-gene Tier A | Recomputed under the frozen strict-33 Tier A, all six bases | Task 8.7 landed |
| 2 | "more sections revealed dependence" | "the sender-set fix kills independence on the published two sections; coverage then makes it certain" | §3.1 — the decomposition was not previously separated |
| 3 | Restatement rested on "one of four pairs below chance at 0.91×" | Plank dead (0.972, p = 0.10); replaced by SenePy vs DeepScence at 0.737 (z = −15.1), with its absolute-score caveat | §2.3 |
| 4 | "the direction of each pair is predicted by its depth loading" | **Refuted** at pair level (permutation p = 0.30; within-pair ρ negative in 5/5). Replaced by the heterogeneous-mechanism claim | §2.3 |
| 5 | "circular Tier A" | "contaminated sender/response split" / "hollow sender set", plus the disclosure that C6 *raised* CoreScence circularity 79 → 88 % | §0.3, §3.4 |
| 6 | "−0.070 SD in negative-control probes" | "−0.074 SD in pooled negative-control features"; the 40 probes alone are flat (−0.023, p = 0.129) | audit R1 — was manuscript-bound |
| 7 | "the reportable-fit filter admits 2–3× more fits than nominal" | Struck; the filter measures 3.0–13.3 %, the 9–16 % is the estimator's two-sided CI-exclusion rate | audit R6 |
| 8 | "nobody in this literature reports it" (A7) | Struck; narrowed to the negative-control-*outcome* construction | novelty review U1 |
| 9 | §5 "N3-tile 0.974 vs published 1.000" and "figures regenerated at 8.7" | N3-**var** 0.996 is now the recommended primary; the figures **have** been regenerated and moved | torus study + 8.7 completion — §5/§6 were stale on these two points despite being marked current |
| 10 | Nothing on N2-vs-N5's second confirmation | §5, with the 1.6 % vs 85.4 % comparison | `CS_PHASE8_COMPMATCH.md` |

**Provenance of §1–§3.** Every value was re-derived here from
`caller_agreement_matched_significance_{verify2sec,2sec_c6,11sections}.csv` by
Mantel–Haenszel pooling and from
`caller_within_type_depth_bias_11sections.csv`, independently of
`caller_coverage_gate{,_headline}.csv`, and matches those files exactly. The
three statistical tests in §2.3 are new computations reported here in full so
they can be re-run.
