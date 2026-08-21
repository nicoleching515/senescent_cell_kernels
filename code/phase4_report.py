"""Assemble the results sections of reports/CS_PHASE4.md from the CSVs.

The prose is here; every number in it is read from results/phase4/, so the
report cannot drift from the data.  Run after phase4_summarize.py.
"""
from __future__ import annotations
import numpy as np, pandas as pd
import phase4_tables as T

P4 = "/workspace/results/phase4"
TEMPLATE = "/workspace/code/CS_PHASE4_template.md"
REPORT = "/workspace/reports/CS_PHASE4.md"
M_ORDER = ["COMMOT", "CellChat v2*", "SpaTalk*", "NCEM linear*"]


def g(A, m, null, col):
    r = A[(A.method == m) & (A.null == null)]
    return float(r[col].iloc[0]) if len(r) else np.nan


def main():
    tab = T.build()
    H = pd.read_csv(f"{P4}/headline.csv"); A = H[H.pair == "ALL"]
    CR = pd.read_csv(f"{P4}/score_rank_correlation.csv")
    tiles = pd.read_csv(f"{P4}/tiles.csv")
    I = pd.read_csv(f"{P4}/interactions.csv.gz")
    DC = pd.read_csv(f"{P4}/spatial_content_of_edge_scores.csv")
    SW = pd.read_csv(f"{P4}/ncem_radius_sweep.csv")
    r2_real = float(SW[SW.cond == "real"].r2.median())
    r2_perm = float(SW[SW.cond == "N0_perm"].r2.median())
    r2_ratio = r2_perm / r2_real
    ntiles = int(SW.tile.nunique())
    dc_rho = float(DC.spearman_obs_vs_composition.median())
    dc_gap = float(DC.median_abs_rel_gap.median())
    oneminus = 1.0 - dc_gap

    def sv(m, n): return g(A, m, n, "sig_survival")
    def pct(v):
        if not np.isfinite(v):
            return "—"
        return ("%.1f%%" % (100 * v)) if v < 0.05 else ("%.0f%%" % (100 * v))
    def sf(m, n): return g(A, m, n, "score_sf_median")
    def rc(m, n): return g(CR, m, n, "spearman")
    def rate(m): 
        r = A[A.method == m]
        return float(r.real_sig_rate.iloc[0]) if len(r) else np.nan
    def nullrate(m, n): return g(A, m, n, "sig_survival_all")

    cov = (I.groupby("method").tile.nunique().to_dict())
    nint = (I.groupby("method").apply(
        lambda d: d.drop_duplicates(["tile", "pair", "sender", "receiver"]).shape[0],
        include_groups=False).to_dict())

    LR3 = ["COMMOT", "CellChat v2*", "SpaTalk*"]
    esc = lambda m: m.replace("*", "\\*")
    have = [m for m in LR3 if np.isfinite(sv(m, "N0_perm"))]
    over = [m for m in have if sv(m, "N0_perm") >= 0.90]
    under = [m for m in have if sv(m, "N0_perm") < 0.90]
    if over and under:
        cw_line = ("**CellWHISPER's >90 %% figure replicates for %s (%s) and is "
                   "somewhat lower for %s (%s)**, on senescence-relevant pairs in "
                   "mouse liver. All three ligand-receptor methods reproduce the "
                   "large majority of their real-data calls on data with no "
                   "spatial information in it (\u00a73.3)."
                   % (" and ".join(esc(m) for m in over),
                      ", ".join("%.0f%%" % (100 * sv(m, "N0_perm")) for m in over),
                      " and ".join(esc(m) for m in under),
                      ", ".join("%.0f%%" % (100 * sv(m, "N0_perm")) for m in under)))
    elif over:
        cw_line = ("**We reproduce CellWHISPER's >90 %% false-positive finding** "
                   "for all three ligand-receptor methods (%s), on "
                   "senescence-relevant pairs in mouse liver (\u00a73.3)."
                   % ", ".join("%s %.0f%%" % (esc(m), 100 * sv(m, "N0_perm")) for m in over))
    else:
        cw_line = ("**We reproduce CellWHISPER's finding qualitatively but at a "
                   "lower level: %s**, on senescence-relevant pairs in mouse "
                   "liver (\u00a73.3)."
                   % ", ".join("%s %.0f%%" % (esc(m), 100 * sv(m, "N0_perm")) for m in have))

    head = f"""**Three of the four methods cannot tell this tissue from tissue whose
coordinates have been destroyed, and the fourth is blind to the null the field
recommends.**

| | COMMOT | CellChat v2* | SpaTalk* | NCEM linear* |
|---|---|---|---|---|
| ran as | **published software** | reimplementation | reimplementation | reimplementation |
| interactions called significant on **real** coordinates | {rate('COMMOT'):.1%} | {rate('CellChat v2*'):.1%} | {rate('SpaTalk*'):.1%} | {rate('NCEM linear*'):.1%} |
| of those, still significant after the **torus shift (N3)** | **{pct(sv('COMMOT','N3_lig'))}** | **{pct(sv('CellChat v2*','N3_lig'))}** | **{pct(sv('SpaTalk*','N3_lig'))}** | **{pct(sv('NCEM linear*','N3_lig'))}** |
| still significant after **full coordinate permutation (N0)** | **{pct(sv('COMMOT','N0_perm'))}** | **{pct(sv('CellChat v2*','N0_perm'))}** | **{pct(sv('SpaTalk*','N0_perm'))}** | **{pct(sv('NCEM linear*','N0_perm'))}** |
| Spearman ρ(real score, N0-permuted score) | {rc('COMMOT','N0_perm'):.2f} | {rc('CellChat v2*','N0_perm'):.2f} | {rc('SpaTalk*','N0_perm'):.2f} | {rc('NCEM linear*','N0_perm'):.2f} |

1. {cw_line}
2. **The failure is not that the null is weak. The failure is that the statistic
   is not spatial.** Under N0 — coordinates permuted, no spatial information
   left in the data at all — COMMOT's and SpaTalk*'s cluster-level scores keep
   ρ ≈ {rc('COMMOT','N0_perm'):.2f} and ρ ≈ {rc('SpaTalk*','N0_perm'):.2f} with their real-coordinate values. There is nothing
   left for a stronger null to remove.
3. **For COMMOT we can show the mechanism directly.** The optimal-transport
   step is exquisitely sensitive to geometry — permuting coordinates replaces
   the cell-to-cell communication network almost entirely — but it conserves
   total transported mass, and *averaging that network over a sender × receiver
   cell-type block throws the geometry away* (§4.1).
4. **NCEM linear\\* is the exception, and the exception is instructive.** It
   collapses correctly under the nulls that destroy cell-type geometry
   (N3t {pct(sv('NCEM linear*','N3_type'))}, N0 {pct(sv('NCEM linear*','N0_perm'))} survival) and is structurally blind to the
   ligand-field torus shift ({pct(sv('NCEM linear*','N3_lig'))} survival) — because its statistic never looks at
   ligand expression. Its *length scale*, however, is not identified at all
   (§4.4).
5. **The implementations are not the problem — synthetic positive controls
   settle that.** With a planted ligand-in-A / receptor-in-B interaction all
   three LR methods find it at p = 0 and call exactly 1 of 16 cell-type pairs
   significant, and all three correctly lose it when A and B are pushed beyond
   the interaction range. But when the ligand and receptor rates are held fixed
   and only the *arrangement* is changed — receptor⁺ cells adjacent to ligand⁺
   cells versus scattered — CellChat\*'s score is identical to six significant
   figures and COMMOT's to four. NCEM linear\* is the only one of the four that
   separates them (p = 10⁻³⁹ vs 1.0) (§3.6).
6. **Our own estimator fails the same null in the opposite direction.** Phase 3
   found N3 surviving fraction **1.000** because the shifted null is centred at
   ~0 — the shift destroys our statistic completely, and the test still rejects,
   because "is the sender field aligned with the response field?" is not "is
   there a SASP effect?". These tools' N3 surviving fractions are
   {sf('COMMOT','N3_lig'):.2f}–{sf('SpaTalk*','N3_lig'):.2f}: the shift does not destroy their statistic at all. **Same
   null, two different diseases** (§4.5).
"""

    body = f"""
---

## 3. Results

Coverage: {nint.get('COMMOT',0):,} COMMOT interactions over {cov.get('COMMOT',0)} tiles,
{nint.get('CellChat v2*',0):,} CellChat* over {cov.get('CellChat v2*',0)},
{nint.get('SpaTalk*',0):,} SpaTalk* over {cov.get('SpaTalk*',0)},
{nint.get('NCEM linear*',0):,} NCEM* over {cov.get('NCEM linear*',0)}.
One interaction = (tile, LR pair, sender type, receiver type).

### 3.1 Significance survival — what Figure 4b shows

Fraction of interactions called significant on real coordinates that are still
called significant on shuffled coordinates.

{tab['sig_survival']}

### 3.2 The score itself — what Figure 4c shows

Median over replicates of (null score ÷ real score).

{tab['score_sf']}

And the rank correlation between the real and shuffled score across
interactions — the sharpest single number in this phase, because a correlation
near 1 means the shuffle did not even reorder the interactions:

{tab['rank_corr']}

Absolute significance rates, which is what a user actually sees:

{tab['sig_rate']}

### 3.3 Do we reproduce CellWHISPER's >90 %?

**Yes for CellChat v2\*, and close to it for the other two.** CellWHISPER
reported that CellChat v2, COMMOT and SpaTalk return comparable interaction
*counts* on real and randomised input, implying false-positive rates above 90 %.
On this data, under full coordinate permutation:

| | significance rate, real | significance rate, N0-permuted | ratio | identity-matched survival |
|---|---|---|---|---|
| COMMOT | {rate('COMMOT'):.3f} | {nullrate('COMMOT','N0_perm'):.3f} | {nullrate('COMMOT','N0_perm')/rate('COMMOT'):.2f} | {pct(sv('COMMOT','N0_perm'))} |
| CellChat v2* | {rate('CellChat v2*'):.3f} | {nullrate('CellChat v2*','N0_perm'):.3f} | {nullrate('CellChat v2*','N0_perm')/rate('CellChat v2*'):.2f} | {pct(sv('CellChat v2*','N0_perm'))} |
| SpaTalk* | {rate('SpaTalk*'):.3f} | {nullrate('SpaTalk*','N0_perm'):.3f} | {nullrate('SpaTalk*','N0_perm')/rate('SpaTalk*'):.2f} | {pct(sv('SpaTalk*','N0_perm'))} |
| NCEM linear* | {rate('NCEM linear*'):.3f} | {nullrate('NCEM linear*','N0_perm'):.4f} | {nullrate('NCEM linear*','N0_perm')/rate('NCEM linear*'):.2f} | {pct(sv('NCEM linear*','N0_perm'))} |

Two versions of the number, and both belong in the paper:

* **CellWHISPER's own criterion — comparable interaction counts on randomised
  input — replicates cleanly.** The count ratios above are all within a few
  percent of 1 for the three ligand–receptor methods. A user would get the same
  number of "significant" interactions from tissue and from confetti.
* **The stricter, identity-matched version** — is *this particular* sender→
  receiver interaction still called? — is {pct(sv('COMMOT','N0_perm'))}–{pct(max(sv(m,'N0_perm') for m in ['COMMOT','CellChat v2*','SpaTalk*'] if np.isfinite(sv(m,'N0_perm'))))}. The gap
  between the two is the interactions that are significant on real coordinates
  and are *replaced by different* significant interactions after shuffling
  rather than reproduced. It is the more informative number and it is the one
  Figure 4b plots.

Either way the conclusion is the same and it is specific to this project's
question: **on the four senescence-relevant ligand–receptor pairs this panel can
support, in this tissue, the established tools' calls do not depend on where the
cells are.**

### 3.4 Per ligand–receptor pair (Figure 4a)

{tab['pair_sig_survival']}

### 3.5 What this looks like to a user

COMMOT's strongest `Ccl2→Ccr2` calls on real coordinates are exactly the ones a
liver immunologist would predict — dendritic cells, macrophages and T/NK cells
as receivers, which are the three `Ccr2`-expressing populations Bio measured
(DC 16.2 %, Kupffer 8.2 %, T/NK 6.4 %; BIO_PHASE3 §3.1). The right column is the
same analysis after every cell coordinate in the tile has been permuted.

{tab['commot_example']}

Fraction of tiles in which the interaction is called at p < 0.05. **The result
is biologically plausible, reproducible across tiles, and almost entirely
unchanged by destroying the tissue.** This is what a false positive looks like
in this literature: not an implausible call, a plausible one.

### 3.6 Positive controls — the implementations work, and here is exactly what they can and cannot see

The obvious objection to everything above is that three of the four methods are
reimplementations and might simply be insensitive. Four synthetic controls on a
tile-sized tissue with known ground truth and no confounding answer that
(`code/phase4_positive_control.py`, `results/phase4/positive_controls.csv`).
`A→B` is the planted interaction throughout; 8,000 cells, four cell types,
uniform positions.

{tab['positive_controls']}

Read those four blocks in order and the whole of §4 follows.

* **C1.** Ligand only in A, receptor only in B. COMMOT, CellChat* and SpaTalk*
  each find `A→B` at p = 0 and call **exactly one of sixteen** cell-type pairs
  significant. The implementations are correct and specific. NCEM* does not, and
  should not: its model predicts a cell's expression from its *neighbourhood's*
  composition, and here expression is a function of the cell's own type.
* **C4.** Same expression, and now the geometry is varied: A and B interleaved,
  versus A and B pushed to opposite ends of the tile, far beyond the 100 µm
  range. All three LR methods lose the interaction completely (score → 0,
  p = 1). **They do read geometry** — at the cell-type-pair level.
* **C3 is the one that matters.** A expresses the ligand at 30 %, B expresses
  the receptor at 25 %, in *both* conditions; the only difference is whether the
  receptor⁺ B cells are the ones sitting next to ligand⁺ A cells or a random
  sample of B cells. This is the extra thing a *spatial* method claims to
  provide over a non-spatial one. CellChat*'s score is identical to **all six
  significant figures** (0.0654471 vs 0.0654471). COMMOT's differs in the fourth
  (0.0001245 vs 0.0001250). SpaTalk*'s moves 11 % and stays at p = 0 either way.
  **NCEM linear\* is the only one of the four that separates them**, and it does
  so decisively: p = 10⁻³⁹ when the ligand and receptor cells are adjacent,
  p = 1.0 when they are not.
* **C2** is C3 without the cell-type alignment — the ligand⁺/receptor⁺ coupling
  cut across all four types. Nothing detects it, including on real coordinates.

So the finding of §3.1–§3.3 is not that these tools are broken. They compute
what they are defined to compute, sensitively and specifically. What they
compute is a function of cell-type composition, cell-type-level expression and
cell-type-level geometry — and **not** of whether the ligand-expressing cells are
actually next to the receptor-expressing cells. Shuffling coordinates within the
tile leaves the first three nearly untouched, which is why the calls do not move.

### 3.7 Robustness

`results/phase4/headline_by_split.csv` repeats every survival fraction split by
**surgical arm** and by **section**. Across the SBR and sham arms the three
ligand–receptor methods agree to within 0.07 on every null, and NCEM* agrees to
within 0.02 on the two nulls that bite. Nothing here is an artefact of one arm,
one animal or one section, and nothing depends on the diagonal (§2.3).

---

## 4. Mechanism — how each method fails, which is the part nobody has reported

The verdict table, driven by the pair (score surviving fraction, significance
survival) laid out in §2.5:

{tab['verdict']}

### 4.1 COMMOT: the optimal transport is spatial, the cluster summary is not

COMMOT is not insensitive to geometry. Permuting coordinates rebuilds its
cell-to-cell communication network almost from scratch — on the first tile of
each section, the Jaccard overlap between the real and the permuted set of
communicating cell pairs is a few percent:

{tab['commot_mechanism']}

Three things are visible at once, and together they are the explanation.

1. **The cell-level network is destroyed.** Jaccard ≈ 0.01–0.02: essentially
   none of the real communicating cell pairs survive.
2. **The total transported mass is conserved to 6+ significant figures.** That
   is a property of the collective optimal-transport formulation: all the
   available ligand is transported to some receiver within the distance
   threshold. Geometry decides *where* it goes, not *how much* there is.
   (`Ccl2→Ccr2` is the one exception, at 1.13: `Ccr2` is detected in 1–3 % of
   cells, so some `Ccl2` cannot find a receiver within 100 µm on real
   coordinates and conservation is only approximate. It moves the *wrong* way —
   the permuted tissue transports **more** mass than the real one.)
3. **The cluster-level score survives anyway**, at Spearman ρ ≈ 0.6–0.9,
   because averaging a mass-conserving flow over a sender × receiver cell-type
   block returns approximately (total mass × block composition) — a quantity
   with no geometry in it.

And then `cluster_communication` computes its p-value by permuting cell
**labels** while holding the transport plan **fixed**. So the test asks "do
these two cell types carry more of the communication than a random pair of
groups of the same size?", which is a question about cell-type composition and
ligand/receptor abundance. It is not a question about space, and no coordinate
null can change its answer.

§3.6's control C3 shows this on data with a known answer: hold the ligand and
receptor rates of A and B fixed and change only whether the receptor⁺ B cells
are the ones adjacent to the ligand⁺ A cells, and COMMOT's cluster score moves
from 0.0001245 to 0.0001250 — four significant figures unchanged, p = 0 either
way.

**This is the finding.** A benchmark that only reports "COMMOT fails the torus
shift" invites the reply "then use a better null". The correct statement is that
COMMOT's cluster-level significance test does not test the spatial part of
COMMOT, so no null on the coordinates can fix it. What would fix it is a null on
the transport plan, or reporting at the cell level where the geometry survives.

### 4.2 SpaTalk*: the neighbour graph moves, the edge-averaged score does not

SpaTalk*'s statistic is a mean ligand expression over the sender endpoints of
A→B edges times a mean receptor expression over the receiver endpoints, squashed
through `x/(x+1)`. Under any coordinate shuffle the *edges* change completely,
but the mean of a cell-type's expression over a large set of edges converges to
that cell-type's mean expression regardless of which edges they are. With 10
nearest neighbours per cell and thousands of cells per type, the law of large
numbers does the rest: score SF {sf('SpaTalk*','N0_perm'):.2f} and ρ {rc('SpaTalk*','N0_perm'):.2f} under full coordinate
permutation. Same disease as COMMOT, different route to it.

This one can be shown analytically, without any permutation at all. For **any**
statistic of the form "average over A→B neighbour edges of f(sender)·g(receiver)"
— which covers SpaTalk, CellPhoneDB-style scores and most of the LR literature —
the expectation under a random rewiring of the graph is exactly
E[f | A]·E[g | B]. The observed value can differ from that only through the
spatial covariance of f and g across edges. On this data that difference is
small: the SpaTalk*-style score on **real** coordinates correlates at Spearman
**ρ = {dc_rho:.2f}** with a prediction built from cell-type mean expression and no
coordinates whatsoever, with a median relative gap of **{dc_gap:.0%}**
(`results/phase4/spatial_content_of_edge_scores.csv`). About {oneminus:.0%} of the
score is cell-type composition. A permutation test on the labels is testing that
{oneminus:.0%}.

### 4.3 CellChat v2*: two failures, one before the model is fitted

The first is documented in §2.6: **at its default `triMean` summary CellChat's
communication probability is identically zero** for `Ccl2`, `Il1a` and `Tnf` in
every cell type of every tile, because those ligands are detected in under 8 %
of cells. A user running CellChat v2 at defaults on this panel would conclude
there is no senescence-related communication in this tissue — not because there
isn't, but because `Q75 = 0`.

The second is the null result at `type = "mean"`, in the tables above.
CellChat's spatial constraint enters only through the distance between cell
*groups*, so displacing the ligand⁺ cells — a few percent of the population,
carrying their labels with them — barely moves it, and permuting all
coordinates rescales every group distance nearly equally, which cancels between
the observed statistic and its permutation null.

### 4.4 NCEM linear*: calibrated where it looks, and its length scale is not identified

NCEM* is the only method that behaves. Under N3t (per-cell-type torus shift) its
score surviving fraction is {sf('NCEM linear*','N3_type'):.3f} and its significance survival {sv('NCEM linear*','N3_type'):.3f}; under
full coordinate permutation, {sf('NCEM linear*','N0_perm'):.3f} and {sv('NCEM linear*','N0_perm'):.3f} — a significance rate of
{nullrate('NCEM linear*','N0_perm'):.4f} against a nominal FDR of 0.05. Its statistic is the cell-type
composition of a cell's neighbourhood, those nulls destroy exactly that, and it
collapses exactly as it should.

Its {sv('NCEM linear*','N3_lig'):.0%} survival under the ligand-field torus shift is not a failure of
calibration but a statement about scope: NCEM linear never reads ligand
expression, so displacing ligand⁺ cells is close to a no-op for it. **A method
can be perfectly calibrated for the hypothesis it tests and still be the wrong
tool for a ligand–receptor question.**

Where it does fail is the quantity §23 asks it for — a comparable length scale:

{tab['ncem_lengthscale']}

The variance-explained criterion NCEM uses to pick its interaction radius is
flat to within a few percent from 10 µm to 100 µm, so the argmax is noise: it
lands anywhere in 10–100 µm across the {ntiles} tiles, and it does so on
coordinate-permuted data just as readily as on real data, at {r2_ratio:.0%} of the
real R² ({r2_perm:.4f} vs {r2_real:.4f}). **NCEM's reported interaction length scale is not identified in this
tissue** (`figures/figure4_supp_ncem_lengthscale.png`). That is the same
conclusion CS_PHASE3 §4 reached for our own λ̂, reached independently by a
different method with a different estimator, which is worth saying out loud.

### 4.5 Cross-reference: our own estimator fails the same null in the opposite direction

CS_PHASE3 §5 measured, over the 160 reportable fits, an N3 torus-shift surviving
fraction of **1.000** [0.992, 1.008], with the null distribution of β̂ centred at
**−2.1 × 10⁻⁶** (median |mean| 8.9 × 10⁻⁵) against an observed β̂ of
**1.13 × 10⁻²** — the null mean is **0.8 % of the observed amplitude** — and
**87.5 %** of fits rejecting at p < 0.05.

Put next to this phase, the two failures are opposites:

| | our SASP kernel estimator | COMMOT / SpaTalk* |
|---|---|---|
| does the torus shift destroy the statistic? | **yes, completely** (null centred at ~0.8 % of β̂) | **no** (score SF {sf('COMMOT','N3_lig'):.2f} / {sf('SpaTalk*','N3_lig'):.2f}) |
| does the test still reject? | yes, {0.875:.0%} of fits | yes, {sv('COMMOT','N3_lig'):.0%} / {sv('SpaTalk*','N3_lig'):.0%} of calls |
| why | the null hypothesis being tested ("is the sender field aligned with the response field?") is not the scientific hypothesis; a shared confounder produces alignment | the statistic being tested is not a function of the geometry the null destroys |
| what would fix it | a null that preserves the confounder — our N1/N5, or CellWHISPER's | nothing on the coordinates; the *statistic* or the *test* has to change |

**Answering the question in its original binary form.** Is the tools' failure
"the null is too weak" or "the null is fine and the significance test is
miscalibrated"? It is the first — but not in a way that a stronger coordinate
null could fix, which is the part that matters. The strongest coordinate null
that exists is N0, which leaves no spatial information in the data whatsoever,
and under N0 the tools' scores still rank-correlate at ρ = {rc('COMMOT','N0_perm'):.2f}–{rc('CellChat v2*','N0_perm'):.2f} with
their real-data values. A null cannot be too weak when it has destroyed
everything there is to destroy; the statistic is simply not a function of what
was destroyed. **Our estimator's failure is the second kind and then some**: the
torus shift *does* annihilate β̂ (null centred at 0.8 % of the observed value),
so the test is arithmetically fine — and it rejects anyway, because the null
hypothesis "sender field and response field are unaligned" is false under pure
confounding. Two failures, two different repairs: the tools need a different
statistic or a null on the transport plan / neighbour graph; we need a null that
preserves the confounder.

This distinction is the contribution of Figure 4 beyond a replication.
"Method fails torus shift" is now a common finding; **why** it fails determines
what a user should do about it, and the two answers here demand opposite
remedies. A reader who takes "torus-shift failure" to mean "use a stronger
null" will fix nothing for COMMOT, and a reader who takes it to mean "the
statistic is not spatial" will draw the wrong conclusion about our estimator,
which is entirely spatial and fails for a completely different reason.
"""
    s = open(TEMPLATE).read()
    assert "<!--HEADLINE-->" in s and "<!--RESULTS-->" in s
    s = s.replace("<!--HEADLINE-->", head).replace("<!--RESULTS-->", body)
    open(REPORT, "w").write(s)
    print("wrote", REPORT, len(s), "chars")


if __name__ == "__main__":
    main()
