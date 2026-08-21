# BIO Deliverable 6 — Containment Interpretation and Discussion Draft

**Biology collaborator · 2026-08-21 · SASP Spatial Response Kernel**
Master Plan §12 Deliverable 6, §3 (biology primer), §6.3/§6.4, §11, §29, §30 item 6.

---

## Preface: the deliverable as specified no longer exists, and the replacement is harder

Section 12 Deliverable 6 says: *"Take the fitted kernel back to Martin et al.'s framework.
Which limiting mechanisms are consistent: rapid degradation, response threshold, immune
clearance, receiver refractoriness? If λ_proximal > λ_downstream, argue what that means for
graded versus thresholded response."*

**There is no fitted kernel.** Under combined control (N2 + N5 + N6) the surviving fraction of
the naive amplitude is **0.082 [−0.099, 0.249]**, with 35 % of fits at or below zero
(CS Phase 3 §0). Ninety-two per cent of the naive amplitude is nuisance. The controlled
amplitude is **0.027 response-sd at contact** against a naive 0.326, and the block-bootstrap
SE of 0.073 sd puts 80 % power at **0.203 response-sd** — so the reportable quantity is a
**bound of ≤ 0.20 response-sd at cell contact**, not a length constant (CS Phase 3 §0). Under
full control the donor-bootstrap CI on λ̂ spans the entire admissible grid, [7, 50] µm, in
**39 of 42** cell-type × module fits (CS Phase 3 §7). And §6.4's proximal-versus-downstream
comparison — the one Deliverable 6 explicitly asks me to interpret — returns a λ_prox/λ_down
ratio whose donor CI **reaches both ends of the widest ratio the grid can express, [0.14, 7.14],
in 3 of 6 receiver cell types, and includes 1 in all six** (CS Phase 5 §6). That question is
not answerable from this design and no argument should be built on the point estimates, which
are internally incoherent: hepatocytes and macrophages put λ_prox at the grid ceiling and
λ_down at the floor, endothelium the other way round.

The negative is calibrated rather than merely absent. On synthetic tissue with a **planted**
effect, N5 nuisance conditioning leaves **0.826 [0.734, 0.921]** of the naive amplitude; with
**no** effect it leaves **0.412 [0.151, 0.720]**. The observed **0.084** sits **below 100 % of
the 600 planted-effect runs** and at the **23rd percentile of the 120 no-effect runs**
(CS Phase 3 §6). The tissue behaves like the synthetic no-effect case. Selection does not
rescue it: cross-fitting puts the winner's curse at **+0.056** of surviving fraction, in the
*upward* direction, so correcting for it makes the negative stronger (CS Phase 5 §3).

So the question this deliverable actually has to answer is: **what does a non-detection with
this bound license us to say about containment?** That is what follows. Section A is the core
interpretation; Section B is the Discussion draft; Section C the biological half of Limitations;
Section D the decisive future experiment; Section E the engagement with Zhao et al.

Every number below traces to a report in `/workspace/reports/` or to the master plan and is
attributed inline. Facts the project does not own are marked `[VERIFY]` and listed in §F.

---

# A. Containment interpretation

## A.1 What was measured, precisely

We measured, for each non-sender cell, the score of one of seven Tier B response programs as a
function of the Euclidean distance to the nearest transcriptionally called senescent cell,
within receiver cell type, with nuisance covariates conditioned out, in mouse liver.

Four qualifiers are load-bearing and every claim below inherits them.

1. **The readout is transcriptional.** Tier B modules are MSigDB Hallmark-derived gene sets
   (31–190 genes on-panel; Bio Phase 1 §1) scored per cell. We did not measure secreted protein,
   receptor occupancy, pathway phosphorylation, chromatin state, metabolism or motility.
2. **The sender call is transcriptional and is a choice.** Three near-independent callers —
   a curated Tier A arrest/damage score, SenePy, and `Cdkn1a` positivity — agree at or below
   chance (Jaccard 0.93–1.22× of chance after conditioning on cell type **and** depth decile,
   all |Spearman| < 0.03; Bio Phase 3 §4.4). Sender prevalence at the primary call is 4.29 %
   (CS Phase 3 §7 N7).
3. **The observable distance range is [≈7, ≈100] µm.** λ grid [7, 50] µm; window 100 µm, the
   99th percentile of observed distance-to-nearest-sender (CS Phase 3 §2).
4. **The bound is on amplitude, at 80 % power, in within-stratum response-sd units**, computed
   on exponential fits under N2+N5+N6.

## A.2 What length scales we could and could not have detected

This is the part of the interpretation the field usually leaves implicit, and it does most of
the work here.

| Range | Status | Why |
|---|---|---|
| **< ~7 µm** (sub-cellular-diameter; juxtacrine, contact) | **Blind** | The resolution floor is the median nearest-neighbour cell distance, **6.7 µm (SBR 7259) to 9.7 µm (sham 7250)** (Bio Phase 1 §1). The λ grid floor of 7 µm is set by it. The first distance bin (< 5 µm) holds a **median of 44 receivers per section, 0.48 % of the data, 5–63 hepatocytes** (CS Phase 3 §9). Membrane-bound IL-1α signalling, or any contact-dependent mechanism, lives entirely inside our blind spot. |
| **~7–50 µm** (one to three cell diameters; the "tens of microns" of cytokine signalling `[VERIFY]`; the 45 µm senescent-to-macrophage and 53 µm senescent-to-NK distances of the 2026 endometrium study, §3 `[VERIFY]`) | **This is where our power lies** | The full λ grid, 40 log-spaced points. Median distance-to-nearest-sender **32.6 µm**, p90 64.7, p95 75.6 (CS Phase 3 §2). Only **2.4 % of hepatocyte receivers, 3.9 % of endothelial, 6.8 % of macrophage** sit within 10 µm of a sender (CS Phase 5 §4). Bound: **≤ 0.20 response-sd**. |
| **~50–100 µm** (the 102 µm senescent-to-T-helper distance `[VERIFY]`) | **Data present, λ not identifiable** | The grid ceiling is window/2, because an exponential with λ > d_max/2 is not distinguishable from a linear trend over the observed range (CS Phase 3 §2). We can say the *amplitude* contrast across this range is bounded; we cannot fit a length constant in it. |
| **> ~100 µm** (the 211 µm senescent-to-B-cell distance `[VERIFY]`; the 100–150 µm periportal `CDKN1A`⁺ localisation reported for human liver, §11 `[VERIFY]`) | **Structurally unobservable at this sender prevalence** | p99 of distance-to-nearest-sender is **98.3 µm** and **99.8 % of receivers are within 120 µm** (CS Phase 3 §2). There is no far-field. |

The last row is the single most consequential fact in this deliverable and it deserves to be
stated as a biological observation rather than a data limitation. **At 4–5 % sender prevalence
there is no unexposed tissue.** Half of all receiver cells sit within 33 µm — two to three cell
diameters — of a cell we call senescent, and essentially all of them within 100 µm. Every
contrast we can form is a contrast *between two exposed states*, not between exposed and
unexposed. That is not an accident of this dataset: it is the generic geometry of any tissue at
a senescent burden in the Test 3 admissible band, and it follows from the same Poisson relation
that governs the regressor. Median distance to nearest sender = 0.4697 ρ_S^(−1/2), verified on
our own data at **r² = 0.984 across 77 section × sender-definition combinations, slope −0.524
against the theoretical −0.50, observed/Poisson ratio 1.03** (CS Phase 3 §4.2). Taking our
measured cell density from the median nearest-neighbour distance (ρ ≈ (0.4697/6.74)² ≈ 4,900
cells/mm² in 7259) and a 4.29 % sender prevalence, the relation predicts a median
distance-to-nearest-sender of **32.5 µm** against the 32.6 µm measured. The geometry is fully
determined by the calling rate.

**Consequence for interpretation: a flat profile over [7, 100] µm has two readings that this
design cannot separate — no response, or a response that is uniform because the whole
observable field is inside it.** Everything in §A.3 is written with that fork held open.

## A.3 The four candidate limiting mechanisms, one at a time

Martin et al. (*Aging Cell* 2023) built a minimal mathematical model of senescence spread and
showed that current mechanistic understanding of SASP diffusion and binding does not explain how
that spread stays local. The paper is explicit about the limits of its own claim: *"Using our
minimal model we can find parameter regimes that lead to senescence spread from a single cell.
However, we can not determine from this model whether the spread of senescence is controlled
(finite) or uncontrolled."*

> **Correction, from the Deliverable 7 audit (§B4b).** The first draft of this section, following
> master plan §3, attributed the four mechanisms below — rapid degradation, a response threshold,
> immune clearance of senders, receptor-level refractoriness — to Martin et al. **They are not
> that paper's.** In the open-access full text (Europe PMC `PMC10410058`) the string `degrad`
> occurs zero times and `refractor` occurs zero times, and immune clearance appears only as the
> mechanism the paper deliberately brackets out: *"Currently, there are no hypothesised mechanisms
> to contain the spread of senescence in the absence of the immune system."* A ligand-binding
> threshold exists as a model parameter but is never offered as the containment mechanism. The
> four are the master plan's gloss, and more broadly the field's standing candidate list. They
> remain the right hypothesis space for *our* measurement — each makes a distinct prediction about
> the shape and amplitude of the profile we can observe — so §A.3(a)–(d) stand as analysis, now
> under their correct provenance. Martin et al.'s own three proposed resolutions are different
> ones, and **§A.3(e)** takes them separately.

A fitted λ would have spoken to these by locating the containment scale. A null with a bound
speaks to them differently, and — this is the point — it discriminates *by what each mechanism
predicts about the shape and amplitude of the observable profile*, which is a different cut
through the hypothesis space.

### (a) Rapid ligand degradation — **disfavoured in the regime we can see; untouched below it**

Degradation-limited containment is the mechanism our design is best placed to test, because a
degradation-limited exponential field has precisely the functional form we fitted, over
precisely the range where our data are densest. If containment in this tissue were set by
ligand degradation with an effective length constant anywhere in **7–50 µm**, and if that field
produced a transcriptional consequence in the receiver larger than **0.20 response-sd at
contact**, we would have detected it. We did not — the controlled amplitude is 0.027 sd, the
naive gradient loses 92 % of itself to nuisance, and the residual behaves like the synthetic
no-effect case.

Two escapes remain open and both must be stated.

* **Sub-resolution degradation.** If λ < 7 µm the entire field lives inside the first distance
  bin, which contains 0.48 % of the data. Very rapid degradation is not disfavoured by our
  result; it is *invisible* to it. Given that the plan's own Tier C reasoning predicts the
  shortest λ for membrane-bound IL-1α, and that `Il1a` is the best-detected short-range ligand
  we have (2.7 % of cells, `Il1r1` on 26.7 % of hepatocytes; Bio Phase 3 §3.1–3.2), a
  contact-range mechanism is entirely live and entirely unaddressed.
* **Small amplitude.** Degradation-limited containment with a transcriptional footprint below
  0.20 sd is untouched. A response of 0.05 sd is biologically real and beyond this design.

So the honest statement is a **joint** exclusion: *not* "degradation does not limit spread", but
"degradation-limited containment with an effective range of 7–50 µm and a transcriptional
consequence above 0.20 response-sd in any of seven Hallmark-derived programs, in any of nine
receiver cell types, is excluded in this tissue".

### (b) Response threshold in the receiver — **partially constrained; the discriminators are dead**

A thresholded response predicts a step: uniform inside a radius, uniform outside, with the
boundary set by where ligand concentration crosses threshold. Our result constrains this only
where the boundary would fall inside the observable window. A step of amplitude > 0.20 sd with
its edge anywhere in 7–50 µm is excluded on the same amplitude grounds as (a). A threshold
whose radius is below 7 µm, or above ~100 µm — where every cell is already above threshold and
the profile is flat because the tissue is saturated — is not excluded at all, and the second of
those is exactly the ambiguity §A.2 flagged.

More importantly, **the two estimands designed to adjudicate graded versus thresholded response
both failed on identifiability grounds, and that failure is itself the finding.**

* **Superposition versus nearest-sender (§6.3).** This is the one estimand our synthetic study
  showed to be *cleanly* identifiable — 15/15 correct verdicts in every regime, |ΔAIC| of 15–251
  per 1,000 cells when a truth is planted (CS Phase 5 §0, §10.1). On real tissue under full
  control it comes out at **−0.060 per 1,000 cells** at `tierA_p95`, **−0.321** at `cdkn1a_pos`
  and **−0.154** at `senepy_p95` — two to three orders of magnitude below the planted signature,
  with the per-fit bootstrap bound excluding it by ≥ 11× even at the worst caller's 90th
  percentile fit (CS Phase 5 §10.2). Where the superposition regressor does earn its place
  (at two of three callers), **λ̂ rails at the 50 µm grid ceiling in 57–60 % of fits**: it is
  acting as a regional sender-density covariate, not a contact-scale dose (CS Phase 5 §10.3).
  Dose-versus-threshold is not answerable here.
* **λ_proximal versus λ_downstream (§6.4).** Not estimable, as above: ratio CI [0.14, 7.14] in
  3 of 6 receiver types, includes 1 in all six, 9 of 12 length constants on a grid rail.
* **Kernel family.** Under full control the **step/threshold** function is the family that comes
  closest — it wins 90.8 % of fits by AIC — but it still beats the *no-kernel* model in only
  55.9 % of fits, and family choice moves d̂½ by a median **4.4×** (CS Phase 5 §0.5). The
  Phase 2 "spline wins 34/35" marker turned out to be the receiver cell-type composition
  artefact (unstratified the spline wins 95 %, stratified 21 %). The mild preference for a step
  shape is therefore not evidence for a thresholded biological response; it is what AIC does
  when no kernel earns its place and the cheapest one wins by default.

**Verdict: untouched, and we should say we cannot adjudicate it rather than that we favour a
threshold.**

### (c) Immune clearance of senders — **untouched by construction; and its prediction is the flat profile we observe**

Clearance is a mechanism on the *lifetime of the sender*, not on the *range of the ligand*. A
static snapshot cannot observe lifetime. Our cohort has timepoints (2/10/26/52 weeks) but 1–2
animals per arm per timepoint, and the monotone time-course claim was withdrawn when the full
cohort arrived (Bio Phase 3 §0.2, §2.2). There is no temporal leverage in this design.

There is, however, one thing worth saying, and it cuts in an unexpected direction. If
containment were set by clearance rather than by range, then around each *surviving* sender the
ligand field would be at its full, diffusion-limited extent — which is the very thing Martin
et al. show is *not* local. A broad field relative to a 100 µm window predicts a **spatially
flat** response profile within the window, because there is no dose contrast left to measure.
That is what we observe. Our result is therefore **consistent with** clearance-limited
containment — but it is equally consistent with no effect at all, and the design contains no
unexposed reference that would separate them. This is the strongest argument in the whole
deliverable for why the decisive experiment must engineer sparse senders (§D).

We did not measure sender-to-immune-cell nearest-neighbour distances, so we cannot place this
tissue against the endometrium calibration figures (~45 µm to macrophages/monocytes, 53 µm to
NK) `[VERIFY]`. That is a cheap analysis and a reasonable reviewer request; it is not in the
current result set. What we do have is that the immune receiver types show no controlled
response gradient either: macrophages SF(ctrl) 0.001 for the receptor-proximal module, T/NK
−0.104, and unstratified contact amplitudes of +0.206 (macrophages) and **−0.039** (T/NK) —
the latter one of two receiver types where the response *rises* with distance (CS Phase 5 §4, §6).

### (d) Receiver refractoriness / receptor-level limitation — **the mechanism our null most directly favours, with a sharp caveat**

Of the four, receiver-side limitation is the only one whose prediction is *exactly* what we
observe: little or no transcriptional response even at contact, at any distance, regardless of
ligand concentration. A null result raises the relative standing of any hypothesis that predicts
null, and this one predicts null unconditionally rather than conditionally on geometry.

The panel lets us go one step further and ask *where* a receiver-side block would have to sit,
and the answer constrains it usefully. **It is not receptor absence for the canonical axes.**
In SBR section 7259, `Tnfrsf1a` is detected in **54.9 %** of hepatocytes, `Il6st` (gp130) in
**33.6 %** (and 63.6 % of LSEC), `Il1r1` in **26.7 %**, `Tgfbr2` in **24.9 %**, `Tnfrsf1b` in
**38.5 %** of Kupffer cells, `Ccr2` in **16.2 %** of DC (Bio Phase 3 §3.1). Hepatocytes are
receptor-competent for TNF, IL-6, IL-1 and TGF-β. If receiver refractoriness is the mechanism,
it operates **downstream of receptor expression** — desensitisation, negative feedback,
chromatin state, or a signalling block — not by the receptor not being there.

There is one clean exception and it is worth reporting. **`Ackr3` (CXCR7) is detected in ≤ 1.0 %
of cells in every cell type in every section** (Bio Phase 3 §3.1). ACKR3 is the receptor the
2026 brain paracrine-senescence work implicates in secondary-senescence competence (plan §31
ref 18) `[VERIFY: that ACKR3 is the receptor that work identifies as conferring competence]`.
In this tissue that axis is simply not available. A receiver-side explanation for the absence of
CCL2/CXCL12-driven secondary senescence in mouse liver is therefore directly supported by
measurement, and it corrects my own Phase 1 report, which listed the CCL2/CXCR7/DPP4 axis as
"fully covered" on the basis of panel membership rather than detection.

**Verdict: favoured, and refined — receptor presence is not limiting for the main SASP axes in
hepatocytes, so a receiver-side mechanism must act downstream of it, except for the ACKR3 arm
where the receptor genuinely is not there.**

### (e) Martin et al.'s own three resolutions — **two untestable here, one moved upstream**

Martin et al. resolve the paradox from the properties of the *secondary* senescent cells rather
than from the ligand field. Their three answers are: (i) juxtacrine secondary senescent cells act
as a firebreak because they do not themselves produce SASP (*"analogous to removing trees to
prevent the spread of a forest fire"*); (ii) induction is time-delayed rather than instantaneous,
so that *"dynamic, time-dependent paracrine signalling prevents the uncontrolled spread of
senescence"*; and (iii) secondary senescent cells secrete fewer SASP molecules than primary ones.

All three turn on a **distinction our sender call cannot make.** Tier A calls senescence from
arrest and damage genes and has no primary-versus-secondary axis. The one module that nominally
carries that axis, B7 `secondary_senescence`, shares 14 of its 38 genes with the caller and is
therefore partly circular (`CS_PHASE3.md` §7 recommends it never be used as the primary caller),
and the primary/secondary distinction itself rests on a single meeting abstract
(`BIO_DELIVERABLE7_CLAIM_AUDIT.md` §B8). We cannot label a cell as secondary, so we cannot test a
firebreak built out of secondary cells.

What our result does is **move the question upstream.** Mechanisms (i) and (iii) both presuppose a
working primary→neighbour induction step that then fails to propagate further: the first ring
lights up, the second does not. Our bound constrains that first ring — the one the firebreak is
supposed to stop — to below 0.20 response-sd at contact range, in the modules where induction
would be transcriptionally visible. If the first ring is already empty at that amplitude, then
containment *in the regime we can observe* does not require a firebreak; something limits the
field before secondary cells would get the chance to contain it. This does not contradict Martin
et al. Their model is not transcriptional, its containment scale may sit below our seven-micron
resolution floor, and the model itself declines to say whether spread is finite. It does mean the
firebreak is not the mechanism our data call for.

Mechanism (ii), time delay, is untouched by construction. A static snapshot cannot see induction
kinetics. Our four timepoints (2 / 10 / 26 / 52 weeks) are spaced orders of magnitude coarser than
the delays the model turns on; they bracket the arrest and secondary-senescence timescale but
cannot resolve a delayed induction step.

Mechanism (iii) is adjacent to §A.4 below and inherits its caveat in full: the ligand-transcript
evidence is thin for *all* our senders, and we can separate neither primary from secondary nor
genuinely low secretion from a detection limit.

## A.4 A possibility neither list enumerates: the senders may not be secreting

All four mechanisms in §A.3(a)–(d) presuppose that senders secrete, and Martin et al.'s (iii)
presupposes that primary senders do. This tissue raises the prior question, and it needs to be
handled with more care than any other statement in this document.

In SBR 7259 (127,386 cells), among cells that are both ligand⁺ and `Cdkn1a`⁺, **`Il6` is
detected in 3 cells, `Cxcl2` in 8, `Mmp3` in 3** (Bio Phase 3 §3.2). `Il6` — the canonical SASP
ligand, and the ligand that motivates Tier B module B2 — is detected in 0.087 % of all
analysable cells. Ligand enrichment in `Cdkn1a`⁺ senders over `Cdkn1a`⁻ cells of the same cell
type is real but modest: `Il1a` 1.5–2.9×, `Ccl2` 1.3–2.7×, `Cxcl1` 1.0–1.8×, `Tgfb1` 1.0–1.7×.
Senders are enriched for SASP ligands; they are nowhere near being defined by them.

**This is not evidence that the senders are not secreting, and it must not be written as if it
were.** Three reasons. First, non-detection on a targeted imaging panel is a statement about
probe sensitivity at that transcript's abundance, not about transcription — `Cxcl1` and `Cxcl2`
encode functionally interchangeable analogues of the same absent human ligand and differ by
4–9× in nearest-neighbour distance in every section purely because one is detected in 5.8 % of
cells and the other in 0.083 % (Bio Phase 3 §3.3). Second, cytokine mRNAs are characteristically
low-abundance and short-lived relative to the potency of the secreted protein `[VERIFY]`, so
mRNA counts are a poor proxy for secretory flux. Third, and most fundamentally, the SASP is
defined proteomically — Basisty et al.'s SASP Atlas is the reference for what is *secreted*, as
opposed to which transcripts go up (plan §3) — and the two are not in register.

What we can say is the weaker, defensible version: **on this platform, in this tissue, the
transcriptional evidence that the called senders are actively producing the canonical SASP
ligands is thin, and this is a hypothesis the design cannot separate from a detection limit.**
It belongs in Limitations, not in the interpretation of containment.

## A.5 The verdict, assembled

Every row below is conditional on the same clause: *given that any response would have to be
transcriptionally visible in a Hallmark-derived module at > 0.20 response-sd*. That conditional
is not a hedge; it is the actual scope of the measurement, and dropping it converts a defensible
statement into an indefensible one.

| Mechanism | Our bound | Reasoning |
|---|---|---|
| Rapid ligand degradation, effective range 7–50 µm | **Disfavoured** | This is exactly the form and range we fitted, with our densest data. Amplitude bounded at 0.20 sd; controlled amplitude 0.027 sd; the residual matches the synthetic no-effect calibration. |
| Rapid degradation, range < 7 µm (juxtacrine/contact) | **Untouched** | Below the resolution floor. The < 5 µm bin is 0.48 % of the data. |
| Response threshold with boundary in 7–50 µm | **Disfavoured** | A step > 0.20 sd anywhere in the window would have been seen; under control the step family beats the no-kernel model in only 55.9 % of fits. |
| Response threshold with boundary outside the window, or graded-versus-thresholded as a question | **Untouched — and not adjudicable here** | Superposition vs nearest is 2–3 orders of magnitude below its planted signature at three callers; λ_prox/λ_down CI spans the full [0.14, 7.14] grid range. |
| Immune clearance of senders | **Untouched; its signature is consistent with what we see** | A lifetime mechanism, invisible to a snapshot. Clearance-limited containment predicts a field broad relative to the window and therefore a flat profile — which is also what "no effect" predicts. No unexposed reference exists to separate them. |
| Receiver refractoriness downstream of receptor expression | **Favoured** | The only candidate predicting null unconditionally. Receptor presence is *not* limiting: `Tnfrsf1a` 54.9 %, `Il6st` 33.6 %, `Il1r1` 26.7 %, `Tgfbr2` 24.9 % on hepatocytes. |
| Receiver refractoriness by receptor absence, ACKR3 axis specifically | **Favoured, directly measured** | `Ackr3` ≤ 1.0 % in every cell type in every section. |
| Senders not secreting (in neither list; cf. Martin's mechanism (iii)) | **Raised, not established** | `Il6` in 3 senescent cells per section, `Cxcl2` in 8, `Mmp3` in 3 — inseparable from a detection limit on a targeted panel. |
| Anything acting on protein, metabolism, chromatin, motility, or on a timescale a snapshot averages over | **Entirely untouched** | Not measured. |
| *Martin (i)* — SASP-poor secondary cells as a firebreak | **Not testable here; and not required in the regime we can see** | No primary/secondary axis in the sender call (§A.3(e)). The first ring the firebreak is meant to stop is itself bounded at 0.20 sd. |
| *Martin (ii)* — time-delayed induction | **Untouched by construction** | Static snapshot; timepoints 2–52 weeks are far coarser than the model's delays. |
| *Martin (iii)* — secondary cells secrete less than primary | **Not testable here** | Requires the same missing label, and inherits the §A.4 detection-limit caveat. |

**The containment paradox is sharpened rather than resolved.** In a tissue where 4–5 % of cells
carry a senescence call and every remaining cell sits within ~100 µm of one — half of them
within 33 µm, two to three cell diameters — the inflammatory, interferon, arrest, ECM and
oxidative-stress programs of the neighbours are, after conditioning, essentially indifferent to
how close the nearest senescent cell is. Whatever limits spread in mouse liver holds the
transcriptional bystander effect below 0.20 response-sd at contact range. That is a
quantitative containment statement, it is the kind the containment paradox asks for, and it is
obtainable from a negative result in a way it is not obtainable from a fitted λ alone.

---

# B. Discussion draft

*Prose, for a co-author to edit. Target ~700 words in the paper; this draft runs longer so that
material can be cut rather than invented.*

---

**Discussion.**

We set out to estimate the length constant of the senescence-associated secretory phenotype's
spatial response kernel, and we report instead a bound and an identifiability failure. In mouse
liver, across eleven animals and six sections admissible under a prevalence-based inclusion
rule, the naive relationship between a neighbour's inflammatory or arrest program and its
distance to the nearest senescent cell is exactly what the field reports: monotone in 32 of 35
section-by-module fits, with a rank correlation below −0.92 in 24 of them, and an effect
twenty-six times the standard error of the distance bins. Under control that relationship
loses ninety-two per cent of its amplitude. Two thirds of what the unstratified curve displays
is receiver cell-type composition — replacing every cell's response with its cell type's section
mean, an operation containing no signalling of any kind, reproduces seventy-six per cent of the
contact amplitude — and of what remains, the largest single contributors are transcript depth,
cell and nucleus area, and local cell density, none of which is a covariate that current spatial
communication tools offer. The controlled amplitude is 0.027 response standard deviations at
cell contact, against a design that has eighty per cent power at 0.203, and no length constant
is identified anywhere on the admissible seven-to-fifty-micron grid: the donor-bootstrap
interval on λ̂ spans the entire grid in thirty-nine of forty-two fits.

That this is a bound rather than a failure to look is established by calibration rather than by
assertion. On synthetic tissue with a planted effect, the same nuisance conditioning leaves
eighty-three per cent of the naive amplitude; with no effect planted it leaves forty-one per
cent. The observed value of eight per cent lies below every one of six hundred planted-effect
runs and at the twenty-third percentile of the no-effect runs. Correcting for the winner's curse
introduced by our own significance screen moves the estimate further in the negative direction,
not towards rescue.

**What this implies for senescence containment.** Martin and colleagues showed with a minimal
mathematical model that current mechanistic understanding of SASP diffusion and binding does not
explain how the spread of senescence stays local, and proposed that containment arises from the
properties of the secondary senescent cells themselves — that they are poor SASP producers, and
that induction is time-delayed rather than instantaneous. The standing candidate mechanisms for
what limits the *ligand field* are different ones: rapid ligand degradation, a receiver response
threshold, immune clearance of senders, and receiver refractoriness. A
fitted length constant would have located the containment scale. A calibrated non-detection
discriminates differently, and in some respects more sharply, because each mechanism makes a
different prediction about the *shape and amplitude* of the profile we can observe. Our bound
disfavours containment set by ligand degradation over an effective range of seven to fifty
microns, if that field has any transcriptional consequence larger than a fifth of a response
standard deviation: that is the regime our estimator was built for, where our data are densest,
and it is empty. It disfavours a thresholded response whose boundary falls within the same
range. It leaves untouched every mechanism operating below our seven-micron resolution floor,
which is set by the median nearest-neighbour distance between cells and which excludes
contact-range and membrane-bound signalling from consideration entirely. It leaves untouched, by
construction, any mechanism that acts on sender lifetime rather than ligand range, since a
static snapshot cannot see lifetime. And it is most consistent with a receiver-side limitation —
the one candidate that predicts an absent response unconditionally rather than conditionally on
geometry. Here the panel adds a constraint worth having: receptor presence is not what limits
the response, because more than half of hepatocytes carry detectable `Tnfrsf1a` and a quarter to
a third carry `Il6st`, `Il1r1` and `Tgfbr2`, so any receiver-side block must act downstream of
receptor expression. The single clear exception is `Ackr3`, detected in one per cent or fewer of
cells of every type in every section, which removes the CCL2/CXCL12–ACKR3 axis implicated in
paracrine senescence elsewhere from consideration in this tissue on measured grounds.

Martin et al.'s own resolutions we cannot test, because our sender call has no primary-versus-
secondary axis and a static snapshot has no access to induction kinetics. We note only that a
firebreak made of SASP-poor secondary cells presupposes a first ring of induced neighbours for it
to contain, and it is that first ring which our bound constrains.

One structural caveat governs all of it, and it is a property of the biology rather than of our
sample size. At the senescent burden that makes a kernel estimable at all, there is no unexposed
tissue. Half of all receiver cells lie within thirty-three microns of a senescent cell and
ninety-nine point eight per cent lie within a hundred and twenty microns; the distance to the
nearest sender is determined to r² = 0.98 by how many cells were called senescent, with the
slope a homogeneous Poisson process predicts. Every contrast available in such a tissue is
between two exposed states, so a flat profile is consistent both with no response and with a
response that is uniform because the entire observable field is already inside it. Distinguishing
those requires a design in which some tissue is genuinely far from any sender, which in turn
requires either engineered focal senescence or a burden below the level at which enough senders
exist to fit. We regard resolving that tension as the central design problem for the next
generation of these experiments, and we specify it concretely below.

**Nothing here speaks to protein.** We measured transcriptional programs in receivers with a
targeted five-thousand-gene panel. A SASP field can exist, diffuse, bind and act without
producing a detectable shift in a Hallmark-derived module score, and the secretome and the
transcriptome are known not to be in register — the SASP is defined proteomically for precisely
that reason. Immediate-early transcriptional responses are also transient, so a snapshot taken
weeks after injury can read zero over a field that is continuously active. "We did not detect a
transcriptional response gradient at seven to a hundred microns, above 0.20 response standard
deviations" is the entire claim, and it should not be paraphrased as an absence of paracrine
senescence.

**What this implies for how spatial senescence claims should be reported.** Three of the
practices this literature treats as controls do not, on our data, control anything. A torus
shift returns a surviving fraction of 1.000 and rotation 0.964 on an effect that nuisance
conditioning shows to be ninety-two per cent nuisance, and on synthetic tissue with no planted
effect at all they return 0.98 and 0.89 — they certify pure confounding as real, because they
test whether the sender field is aligned with the response field, which is not the same question
as whether there is a SASP effect. Matched decoys balanced to a maximum standardised mean
difference of 0.033 across five covariate families absorb six per cent of the effect while
regression on those same covariates absorbs ninety-two per cent; passing a balance diagnostic
bounds nothing, because matching equalises covariates between senders and decoys without
removing the dependence of the response on those covariates at the receiver. And an
expression-matched random gene set of the same size produces a distance gradient as strong as
the real response module between five and seventy-seven per cent of the time, depending on the
module. We therefore recommend that any distance-to-senescent-cell analysis report, as a
minimum: the surviving fraction of the effect after conditioning on receiver cell type,
transcript depth, segmentation area and local density, with a calibration against synthetic data
in which the truth is known; the fraction of fits whose length constant rails against a grid
bound; the relationship between the distance regressor and the sender calling rate; and the
result under at least two independent sender definitions. On our data the surviving fraction
varies by a factor of three to five across sender callers that agree with each other at chance,
and no caller yields a controlled surviving fraction above 0.29.

**Limitations.** [Section C follows here in the paper.]

---

# C. Limitations — the biological half

*Prose, as it should appear in the paper.*

---

Our senescence calls are inferred from transcriptional signatures, not established
experimentally, and on this data that is more consequential than the usual caveat admits. Four
published or panel-standard ways of calling a senescent cell — a curated arrest-and-damage
score disjoint from the response modules, SenePy, `Cdkn1a` positivity, and the CoreScence set
underlying DeepScence — applied to the same 127,000 to 237,000 cells, produce top-five-per-cent
calls that overlap at 0.93 to 1.22 times chance after conditioning on cell type *and* on
sequencing-depth decile. They are, statistically, independent definitions. What each one selects
is identifiable and technical: within cell type, SenePy and `Cdkn1a`⁺ are enriched two-and-a-half
fold in the highest transcript-count quintile and depleted eight- to elevenfold in the lowest,
while the arrest score and DeepScence run the other way. A senescence call on targeted spatial
data is therefore not a noisy measurement of one latent state but a choice of which end of the
detection-depth distribution to name senescent, and any length constant inherits that choice.
We report every headline result under three near-independent callers for this reason. We also
note that DeepScence's own CoreScence gene set overlaps our response modules in fifty-one gene
memberships, and that stripping the shared genes reduces the CoreScence-sender amplitude for the
secondary-senescence readout by a factor of 0.698 — thirty-one per cent of that fit is literally
the same genes on both sides — so DeepScence is used here as a comparison method and never as
the primary caller.

The design is a static snapshot. Two of the four containment mechanisms we discuss — immune
clearance of senders and receiver adaptation — act on timescales a snapshot integrates over, and
we have no leverage on either. Our cohort spans four timepoints from two to fifty-two weeks, but
with one or two animals per arm per timepoint, and we withdraw as unsupported the monotone
time-course we reported from a partial cohort.

The panel is targeted. Five thousand one hundred and six gene-expression features are enough for
every response module to clear a thirty-gene floor, but they are not enough to see the SASP
itself: `Il6` is detected in three senescent cells in a 127,000-cell section, `Cxcl2` in eight
and `Mmp3` in three, while `Cxcl12` is expressed by half of all cells and its
nearest-neighbour distance of five to eight microns falls below the segmentation resolution
floor. `Ackr3`, the receptor with the strongest prior claim on secondary-senescence competence,
is detected in one per cent or fewer of cells everywhere. Consequently no kernel can be tied
back to an interpretable ligand: the distance to the nearest ligand-expressing cell is determined
to r² = 0.987 by that ligand's detection frequency, with a Poisson slope, so a cross-ligand
ordering of length constants would recapitulate probe sensitivity rather than diffusion range —
and empirically it does, with membrane-bound `Il1a` ranking anywhere from first to fourth across
sections of the same arm. The internal control specified in our own design was attempted and
failed for an identifiability reason. Three cell types the analysis should have contained cannot
be called on this panel and were reported as absent rather than guessed at: erythroid cells,
mast cells and lymphatic endothelium; neutrophils likewise. `Krt7` and `Cftr` do not detect, so
we report a `Biliary/ductular` compartment rather than distinguishing true cholangiocytes from
ductular metaplasia — a distinction that matters in this model, since the ductular reaction is
the phenotype.

The biology is one tissue in one species under one insult. This is mouse liver in a surgical
model of intestinal-failure-associated liver disease induced by seventy-five per cent
small-bowel resection, not human tissue and not natural aging, and the injury actively destroys
one of the anatomical covariates the design depends on: the portal-triad landmark, which
validates cleanly in all five sham sections with a zonation correlation of +0.25 to +0.34,
collapses to −0.09 to +0.12 in every resected section from ten weeks on, because the ductular
reaction disperses biliary cells through the parenchyma. Eleven animals contributed one section
each; four sections exceed the twenty-per-cent prevalence ceiling above which
distance-to-nearest-sender is near zero everywhere and the length constant is unidentifiable by
construction, and one falls below the one-per-cent floor, leaving six sections from six animals
in the primary analysis. Re-running the battery on the excluded sections gives surviving
fractions of 0.098 and 0.124 against 0.082 in band, so the exclusion does not drive the result,
but a six-donor bootstrap yields sixty-three distinct resamples and a lumpy interval, and we
label the pooled estimate a case study accordingly. One further confound is specific to this
cohort and generalises beyond it: within the resected arm, section-level senescent burden tracks
section-level median transcripts per cell at a Spearman correlation of +0.94. Animal-level
senescent burden and animal-level detection depth are not separable here, and we do not claim
that the highest-burden section, at forty-five per cent `Cdkn1a`⁺ hepatocytes, is the most
senescent rather than the most deeply sequenced.

Finally, one confound we can rule out. Segmentation bleed-through is the standard explanation
for spurious spatial autocorrelation, and it is absent here: on 129,104,526 transcripts from an
admissible section, 88.27 per cent are assigned to a cell and 11.72 per cent are unassigned,
against a thirty-per-cent threshold, with the rate unchanged by quality filtering. The
confounding we document is not transcripts leaking between cells.

---

# D. The decisive future experiment

*Master Plan §29 objection 5. The plan already names it — "induce senescence at a known focus,
profile spatially at multiple timepoints, and fit the kernel with the confound removed by
design." What follows is what that has to look like given what we now know.*

## D.1 The sender label must be non-transcriptional

This is the requirement everything else depends on, and our own results say why. Distance to the
nearest sender is a **sender-calling-rate readout to r² = 0.984** across 77 section ×
sender-definition combinations, with a slope of −0.524 against the homogeneous-Poisson −0.50
(CS Phase 3 §4.2). Sender prevalence in turn tracks section median transcripts per cell at
**ρ = +0.94** within the SBR arm (Bio Phase 3 §5). Chained together, the model's independent
variable is, to two decimal places of its variance, a function of how deeply the section was
sequenced. No amount of covariate adjustment fixes an independent variable that *is* the
confounder; and moving thresholds does not help — going from p90 to p99 on the same score
changes prevalence from 8.6 % to 0.86 % and median λ̂ from 7.0 to 19.7 µm with nothing biological
having changed (CS Phase 3 §7).

The fix is to define senders by construction rather than by score: an inducible focal
senescence-driving transgene, a lineage or reporter label read out independently of the
expression assay, or engraftment of a labelled senescent population at a known site. With sender
identity assigned by design, sender position no longer depends on detection depth, the r² = 0.98
pathology disappears, and — importantly — the confounder-preserving null problem changes
character. The reason the torus shift fails on our data is that the confounder is
*sender-intrinsic*: our Tier A callers select cells that are 1.6× enriched in the lowest depth
quintile within their own cell type, so the senders are systematically shallow cells sitting in
systematically different neighbourhoods, and relocating them does not touch that. An orthogonally
labelled sender has no such intrinsic property, and a torus shift becomes a legitimate null again.

## D.2 Sender prevalence: the number, and the tension it creates

The experiment needs receivers that are genuinely far from any sender. Using the Poisson relation
we verified on our own data (median distance to nearest sender = 0.4697 ρ_S^(−1/2), r² = 0.984,
observed/Poisson 1.03) and our own measured cell density (from median nearest-neighbour distances
of 6.74 µm in SBR 7259 and 9.66 µm in sham 7250, i.e. ρ ≈ 4,900 and ≈ 2,400 cells/mm²), the
required sender prevalence for a given median distance-to-nearest-sender is:

| target median distance to nearest sender | required sender prevalence (at ρ ≈ 4,900 cells/mm²) | (at ρ ≈ 2,400 cells/mm²) |
|---|---|---|
| 32.6 µm (**what we had**, 4.29 % prevalence) | 4.3 % | — |
| 100 µm | **0.45 %** | 0.93 % |
| 200 µm | **0.11 %** | 0.23 % |
| 300 µm | **0.05 %** | 0.10 % |

*(Derived arithmetic, not a measured quantity; inputs are the two measured medians and the
verified Poisson relation.)*

**This exposes a genuine tension in the design rules we ourselves used.** Section 8 Test 3
requires sender prevalence in the 1–20 % band for the kernel to be estimable — but the table
shows that a *randomly dispersed* sender population at ≥ 1 % prevalence cannot produce a median
distance-to-nearest-sender beyond about 100 µm. The admissible band and the far-field
requirement are almost mutually exclusive under dispersion.

**The resolution is the focus, and this is precisely why "induce at a known focus" is the right
design and not merely a convenient one.** Concentrating the sender mass decouples "enough
senders to have signal" from "enough distance range to have a lever". A single focus of a few
hundred senescent cells in a 200,000-cell section is well under 1 % prevalence yet furnishes a
distance axis running from contact out to millimetres, with a large, genuinely unexposed
receiver population at the far end. Our own senders are the opposite: a Ripley-K ratio at 50 µm
of **1.11** against the stratified-permutation null (1.04–1.19 across sections) puts them at the
very bottom of the clustering axis — they are close to a random thinning of the tissue
(CS Phase 3 §4.2), which is the worst case for this estimand.

Practical specification: **one or a few well-separated foci per section, each of at least a few
hundred labelled senescent cells, in tissue large enough that at least a quarter of receivers of
each cell type sit beyond five expected length constants from the nearest focus** — on the
field's own "tens of microns" calibration `[VERIFY]`, beyond roughly 250 µm.

## D.3 Detection depth for the response, and what it has to buy

Two separate requirements, and the second is the one usually missed.

*Sensitivity.* If the readout is to include the SASP ligands themselves or the immediate-early
NF-κB targets, a panel that detects `Il6` in three senescent cells per 127,000-cell section is
not adequate. The requirement is either single-cell-resolution whole-transcriptome capture, or a
targeted panel with materially better sensitivity for low-abundance cytokine transcripts than
the 446–968 median transcripts and 288 median genes per cell we worked with. Note that
spot-based whole-transcriptome platforms are not a solution: mixing multiple cells per spot
reintroduces the receiver-composition confound that accounts for two thirds of our naive
gradient, and deconvolution adds its own spatially structured confound.

*Statistical resolution, and why not to solve it with n.* Our block-bootstrap SE on the
controlled amplitude is 0.073 response-sd, from 848,596 pooled receivers and a median of 63,848
cells per hepatocyte fit. Since the SE scales as n^(−1/2), resolving 0.05 response-sd at the same
power would need roughly sixteen times the cells — which is the wrong lever to pull. The right
lever is the estimand: with a focus design the contrast becomes *exposed versus unexposed*
rather than *near versus slightly less near*, so the quantity being estimated is the full effect
rather than the difference between two saturated states. That change is worth more than an order
of magnitude of sample size.

*Depth as a covariate, still.* Even with sender identity assigned by design, per-cell transcript
depth, genes detected, cell area and nucleus area must remain in the model: *within* receiver
cell type, that block on its own removes 71 % of the naive effect (SF 0.288), and the
local-density block on its own removes 78 % (SF 0.219) (CS Phase 3 §3.2). Section-level median transcripts per cell should be a covariate at
the section level too.

## D.4 Timepoints, replication and controls

Multiple timepoints are required for two independent reasons: immune clearance and receiver
adaptation are lifetime and adaptation mechanisms invisible in a snapshot, and immediate-early
transcriptional responses are transient enough that a single late sample can read zero over an
active field. Sampling should bracket both regimes — hours to days for receptor-proximal
NF-κB targets, weeks for arrest and secondary senescence.

Replication has to be at the animal level and it has to be adequate: our six-donor bootstrap
yields 63 distinct resamples and an interval we were obliged to label a case study. Enough
animals per timepoint that the donor bootstrap is not the limiting resolution is the criterion,
and our experience says one or two is not it.

Two controls are non-negotiable. First, a **sham-induction focus** — the same procedure, same
label, no senescence programme — which supplies the unexposed absolute reference the gradient
design lacks, and turns "flat" into an interpretable statement rather than an ambiguous one.
Second, the **calibrated null battery** run as we ran it: surviving fractions rather than
p-values (all three of our permutation nulls reject at 42–100 % of fits under a synthetic true
null, so p-values carry no information here), reported against synthetic runs with a planted
effect and with none, so that a shrinkage can be read as evidence rather than as disappointment.

---

# E. Engagement with Zhao et al., *Cell* 2024

*Master Plan §29 objection 1 is the frame: "Yes, descriptively, without negative controls or
uncertainty. We reproduce their qualitative finding and then quantify how much survives. That is
a different claim."*

## E.1 What they reported, and what we reproduce

Zhao et al. profiled aging mouse tissue across life stages, defined senescence-sensitive spots,
and showed that SASP score, TNF signalling, ATP biosynthesis and cell-cycle genes vary
monotonically with distance from those spots, consistently across organs, concluding that
senescent foci act as epicentres compromising surrounding cells in a distance-dependent manner
(plan §4.1) `[VERIFY: platform, organ list, exact definition of a senescence-sensitive spot,
and whether any spatial permutation or coordinate-randomisation control was reported]`.

**We reproduce the qualitative phenomenon cleanly.** Across 35 section × module combinations the
binned mean response falls monotonically with distance to the nearest sender in **32 of 35**
fits, with Spearman ρ ≤ −0.92 in 24 of them; the median difference between the first three and
last six bins is 0.056 module-score units against a median bin SEM of 0.0020, a ratio of **26**
(CS Phase 2 §1). Every module shows it, in every section but one. Had we stopped there we would
have reported a clean SASP distance gradient in mouse liver, and the figure would have looked
like the published ones.

## E.2 What we then show, in our data

Under control, 92 % of that amplitude is nuisance, and the decomposition is specific.
**Sixty-six per cent of the unstratified gradient is receiver cell-type composition** — adding
receiver cell-type intercepts alone takes the surviving fraction to 0.344 (CS Phase 3 §3.1) — and
a composition-only surrogate, in which each cell's response is replaced by its cell type's
section mean and which therefore contains no signalling, distance dependence or kernel of any
kind, reproduces **76 %** of the contact amplitude (CS Phase 5 §4). Within receiver cell type the
pooled amplitude falls from +0.260 to +0.091 response-sd, and **two of nine receiver types go
negative** (T/NK −0.039, Biliary/ductular −0.059). Of what survives stratification, transcript
depth and segmentation area alone leave 0.288 and local density alone 0.219 (CS Phase 3 §3.2).
The contact spike itself rests on very little data: the first bin holds a median of 44 receivers,
0.48 % of the data, 5 to 63 hepatocytes per section (CS Phase 3 §9).

## E.3 What this does and does not say about their finding

**It is not a refutation, and the paper must not read as one.** Our result is a measurement in a
different system by a different method, and four differences are individually sufficient to
break any direct inference from ours to theirs. They studied aging across life stages in
multiple organs; we studied a surgical model of intestinal-failure-associated liver disease in
one organ. They used a different platform, and if it is spot-based `[VERIFY]` then the unit of
observation is not a cell and the confound structure is not the same one we decomposed. Their
senescence definition — senescence-sensitive spots — is not our per-cell transcriptional caller,
and our own work shows that near-independent callers on the same cells give surviving fractions
differing by a factor of three to five. And their readouts include quantities we did not score.
A real distance-dependent effect in aged multi-organ mouse tissue is entirely compatible with our
bound: ours is ≤ 0.20 response-sd at contact for seven Hallmark-derived programs in nine receiver
types in mouse liver at 4–5 % sender prevalence, and it is not a bound on their tissues.

**What it does say is about the inference, not the observation.** The specific analytical
pattern — bin a response score by distance to a senescence-defined locus, observe monotone
decline, conclude distance-dependent influence — is not on its own evidence of a signalling
gradient. We produced that pattern at 26× the bin standard error in a tissue where it does not
survive conditioning, and we showed that a surrogate containing no signalling whatsoever
reproduces three quarters of it. The step from the observation to the mechanistic conclusion
therefore needs controls, and on our data the controls the field currently treats as strong do
not supply them: a torus shift returns a surviving fraction of 1.000 and rotation 0.964 on an
effect that is 92 % nuisance, and on synthetic tissue with no planted effect at all they return
0.98 and 0.89 (CS Phase 3 §0, §5). Passing a coordinate-shift null would not have licensed the
claim either.

**And there is a constructive asymmetry worth stating plainly.** The composition confound was
visible to us *only* because we had single-cell resolution with per-cell type labels and could
build the composition-only surrogate. On a spot-based platform, if that is what was used
`[VERIFY]`, both the senescence-sensitive-spot call and the neighbouring spots' scores are
composition-weighted mixtures, and neither the original authors nor a reader could have checked
this from the published analysis without deconvolution — which introduces its own spatially
structured confound (plan §3). This is not a criticism of their execution; it is an argument that
the necessary control is cheap, specifiable, and was not available to be reported. Our
contribution is to specify it: report the surviving fraction after conditioning on receiver
cell-type composition, transcript depth, segmentation area and local density, calibrated against
synthetic data where the truth is known, alongside the raw gradient. Where that has been done —
here — the gradient does not survive; where it has not, the gradient is not yet evidence either
way.

The fair one-sentence version for the paper: *we reproduce the qualitative distance gradient that
motivated this literature, and we show that in our tissue it is 66 % receiver cell-type
composition and 92 % nuisance overall — which does not refute the published result in aged
multi-organ tissue on a different platform, but does mean the observation alone cannot carry the
mechanistic conclusion, and specifies the control that would.*

---

# F. Items marked `[VERIFY]`

For the Deliverable 7 claim audit. None of these is used to support a load-bearing claim without
its marker.

1. ~~**Martin et al. (*Aging Cell* 2023) enumeration.**~~ **RESOLVED — verified error; corrected
   in place.** D7 §B4b checked the open-access full text (`PMC10410058`): the four mechanisms are
   the master plan §3's gloss, not Martin et al.'s. `degrad` and `refractor` occur zero times.
   §A.3 has been re-headed and re-attributed, §A.3(e) now treats Martin et al.'s three actual
   resolutions, and the §A.5 table and §B draft are corrected. The mechanism *analysis* was
   unaffected — only its provenance.
2. **"Cytokine signalling in tissue is generally thought to act over tens of microns."** Plan §3
   states it without a citation. Used in §A.2 and §D.2.
3. **The 2026 human endometrium nearest-neighbour figures** — 45 ± 20 µm (macrophages),
   45 ± 25 µm (monocytes), 53 ± 23 µm (NK), 102 ± 42 µm (T-helper), 211 ± 66 µm (B cells) — and
   whether these are means ± sd. Plan §3 and §31 ref 7. Used as the calibration ladder in §A.2.
4. **The 100–150 µm periportal `CDKN1A`⁺ hepatocyte localisation** attributed to the *Cell
   Genomics* liver paper in plan §11 — and specifically that it refers to the **human** arm,
   whereas GSE310392 is the mouse arm (Phase 0 §5). Cited in §A.2 only as a length scale we
   could not test.
5. **ACKR3/CXCR7 as the receptor conferring secondary-senescence competence** in the 2026 brain
   paracrine-senescence preprint (plan §31 ref 18). Load-bearing for the one directly measured
   receiver-side claim in §A.3(d).
6. **Cytokine mRNAs are characteristically low-abundance and short-lived relative to the potency
   of the secreted protein.** My own background claim, not a project measurement. Used in §A.4 to
   argue *against* over-reading our own ligand non-detections, so the risk of error is
   conservative, but it should still be sourced or cut.
7. **Zhao et al. (*Cell* 2024): platform (spot-based vs single-cell resolution), organ list,
   exact definition of a "senescence-sensitive spot", and whether any coordinate-randomisation,
   torus-shift or permutation control was reported.** Section E depends on all four, and the
   composition argument in §E.3 depends specifically on the platform.
8. **Acosta et al. 2013 transwell result** as the canonical demonstration that paracrine
   senescence does not require contact (plan §3, §31 ref 17) — full bibliographic details are
   missing from the plan's reference list.
9. **Basisty et al. SASP Atlas (*PLoS Biology* 2020) as the reference for what is secreted, as
   opposed to which transcripts go up** (plan §3, §31 ref 16). Used to draw the
   transcript-vs-secretome distinction in §A.4 and §B.
10. **The derived sender-prevalence table in §D.2** is arithmetic, not measurement. It combines
    two measured medians (6.74 / 9.66 µm) with the Poisson relation verified at r² = 0.984, and
    should be checked by whoever owns the Methods, or dropped to a single illustrative row.
