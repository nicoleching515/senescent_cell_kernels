# Phase 0 — Data Acquisition, Storage, and Environment

**Date:** 2026-08-20
**Owner:** CS lead (automated)
**Plan sections exercised:** 0.1 item 1, 7 (Rank 1), 8 (Test 1 prep), 15–17 (RunPod storage), 20 (failure modes)

---

## 1. What was requested

Download `GSM9295281_7259_liver_sbr_Male_26-U1.tar.gz` and
`GSM9295284_7250_liver_sham_Male_26-U1.tar.gz` from GSE310392, natively on the RunPod pod.

## 2. What happened

### 2.1 Do not download `GSE310392_RAW.tar`

GEO exposes every sample's supplementary file standalone at
`ftp.ncbi.nlm.nih.gov/geo/samples/GSM9295nnn/<GSM>/suppl/`. Fetching the two
samples directly avoids pulling the whole series tarball. Confirmed present:

| Sample | File | Remote size |
|---|---|---|
| GSM9295281 (sbr) | `GSM9295281_7259_liver_sbr_Male_26-U1.tar.gz` | 17,302,770,473 B (16 GB) |
| GSM9295284 (sham) | `GSM9295284_7250_liver_sham_Male_26-U1.tar.gz` | 14,212,427,791 B (13 GB) |

### 2.2 First attempt failed silently — the storage quota, not the network

Two concurrent `wget -c` jobs both stalled at ~62–63% and **exited with status 0
while leaving truncated files on disk**. Initial reading was an NCBI-side
connection drop. That was wrong.

**Actual cause: the RunPod network volume has a hard quota, hit at ~20 GB
combined.** `df -h /workspace` reports the *underlying MooseFS cluster*
(2.1 PB, 504 TB free), not the per-volume cap, so there appeared to be room.
The two transfers died at the same instant because they collectively exhausted
the quota. Confirmed by `dd`: writes failed with `EDQUOT` / "Disk quota
exceeded" even at 10 KB.

> **Operational lesson for this pod, add to Section 20's failure table:**
> `df` on `/workspace` is meaningless for capacity planning on a RunPod network
> volume. Probe with `dd` instead. And **`wget` can exit 0 on a quota-truncated
> file** — always verify against the remote `Content-Length`.

### 2.3 Archive contents — Section 17.2's advice is exactly right, and then some

Listing the members changed the whole approach:

**GSM9295284 (sham)** — standard Xenium output bundle, members in this order:

| Member | Size | Needed? |
|---|---|---|
| `cell_boundaries.parquet` | 22 MB | useful |
| `cell_feature_matrix.h5` | 100 MB | **yes** |
| `cells.parquet` | 5 MB | **yes** |
| `morphology.ome.tif` | **11.2 GB** | no (image) |

The three needed files total **127 MB and sit first in the archive**.

**GSM9295281 (sbr)** — full manifest, in archive order:

| Member | Size | Needed? |
|---|---|---|
| `7259_..._processed.rds` | **11.66 GB** | no (R object) |
| `cell_boundaries.parquet` | 11 MB | useful |
| `cell_feature_matrix.h5` | 62 MB | **yes** |
| `cells.parquet` | 2.8 MB | **yes** |
| `morphology.ome.tif` | 3.95 GB | no (image) |
| `nucleus_boundaries.parquet` | 9.5 MB | no |
| `transcripts.parquet` | **1.85 GB** | no — Section 17.2 says explicitly not to pull this |

Both samples *do* ship a standard Xenium bundle; 281 additionally ships a
pre-processed Seurat `.rds` and puts it **first**, so the 76 MB of analysis
files sit behind 11.66 GB of object we do not want. Streamed with a cut at
11.8 GB decompressed: transfers ~11.7 GB, writes 76 MB, never touches the
3.95 GB image or the 1.85 GB transcript table.

> **Note on member ordering:** GSM9295284's listing was taken from a truncated
> partial and stopped at `morphology.ome.tif`. By analogy with 281 it almost
> certainly also carries `nucleus_boundaries.parquet` and `transcripts.parquet`
> after the image. Neither is needed.

### 2.4 Resolution: stream-extract, never store the archive

```bash
curl -sL "<url>" | gzip -dc | head -c 220000000 \
  | tar -x -C data/raw/ --wildcards '*/cell_feature_matrix.h5' '*/cells.parquet' '*/cell_boundaries.parquet'
```

`head -c` closes the pipe once past the needed members, SIGPIPEs `gzip` and
`curl`, and the 11.2 GB image is never transferred. Because pods have free
egress (Section 14) but a tight disk quota, **trading bandwidth for disk is the
correct optimisation on this platform.**

**Result for GSM9295284 — complete, all three files byte-exact:**

```
22019666  cell_boundaries.parquet
100219736 cell_feature_matrix.h5
5084490   cells.parquet
```

**127 MB transferred and stored instead of 14.2 GB — a 112× reduction.**

## 3. Status

| Item | Status |
|---|---|
| GSM9295284 (sham) — analysis files | **COMPLETE**, verified byte-exact |
| GSM9295281 (sbr) — full manifest | **COMPLETE** (7 members, see 2.3) |
| GSM9295281 — analysis files | streaming, cut at 11.8 GB to skip the `.rds` |
| `morphology.ome.tif`, both samples | deliberately NOT downloaded |
| Python environment | **READY** |
| Volume usage | 131 MB of ~20 GB |

Environment: Python 3.11.10, numpy 1.26.3, scipy 1.17.1, pandas 2.3.3,
pyarrow, h5py, scikit-learn, statsmodels, matplotlib, seaborn, joblib, tqdm,
anndata, scanpy 1.11.5, pointpats 2.5.5. Hardware: 48 cores, 251 GB RAM.

Directory layout per Section 17.1: `data/{raw,interim,processed}`, `results/`,
`code/`, `figures/`, `genesets/`, `logs/`, `reports/`.

Note: packages are installed to the **container disk**, which does *not* survive
pod termination (Section 20). Section 16.2's custom Docker image, or a
`/workspace/envs` conda prefix, is still the durable fix — but `/workspace` is
quota-tight, so the image is the better option here.

## 4. Open issues for the next phase

1. **Sample 281's `.rds` route.** 11.66 GB extracted would consume >half the
   volume quota and needs an R toolchain. Decide whether 281 is worth it or
   whether the analysis proceeds on 284 plus other series samples deposited as
   Xenium bundles.
2. **Series identity.** The plan's Section 7 describes GSE310392 as a *human*
   liver aging/cancer atlas with 43 donors. The deposit we pulled is *mouse*
   (`liver_sbr` / `liver_sham`, ENSMUSG IDs, Xenium Prime **Mouse** 5K panel).
   Assigned to the Bio agent as its top-priority deliverable.
3. **Replication units.** Two samples, one male mouse each, sbr vs sham, is not
   a donor bootstrap (Section 24.1). Establish the true sample count in the
   series before committing.

---

*Both downloads ran genuinely concurrently as separate processes and both are
resolved: 284 complete, 281 pending a structural decision rather than a
transfer.*

---

# Phase 0b — Series Identity Resolved (2026-08-20, later)

## 5. GSE310392 is the mouse arm — plan Section 7 cited the right paper, wrong accession

Per the paper's data availability statement, confirmed against the GEO record:

| Arm | Where it lives |
|---|---|
| Human liver (normal aging) | **SenNet portal**, data provider group = Washington University — *not GEO* |
| Human mCRC liver metastases | **HTAN DCC**, HTAN WUSTL Atlas |
| **Mouse Xenium** | **GEO GSE310392** ← what we have |
| Mouse snRNA-seq | GEO GSE311064 |
| Mouse 2-week snRNA-seq | GEO GSE293958 |

GEO record: series title still reads "Cellular senescence in human liver under
normal aging and cancer **[Xenium]**", but organism is *Mus musculus* and the
platform is **GPL33896, Xenium In Situ Analyzer: Mus musculus**. The "human"
in the title refers to the parent study, not this subseries. That is exactly
how the plan's Rank 1 entry went wrong.

## 6. The model is IFALD via small bowel resection — not CCl4

From the Zenodo code deposit (10.5281/zenodo.17584291, 110 MB, scripts
extracted to `data/raw/zenodo_code/`): directories named `IFALD_mouse`,
scripts named `Mouse_SBR_Recluster_*`.

- **IFALD** = Intestinal Failure-Associated Liver Disease
- **sbr** = small bowel resection (the model); **sham** = surgical control
- Cell types reclustered by the authors: Hepatocytes, Cholangiocytes,
  Hepatic stellate cells, LSECs

> **Limitation of the deposit:** it contains `human_Xenium` merge scripts but
> **no mouse Xenium script**. It gives mouse marker logic and QC thresholds
> from the snRNA-seq arm, but no ready-made senescence calls or cell-type
> annotations for our sections. It does not compress Days 1–3 as much as hoped.

## 7. Series structure — 12 samples, and this fixes the replication problem

| Timepoint | SBR | Sham |
|---|---|---|
| 2 wk | 7361 | 7352 |
| 10 wk | 7448, 7450 | 7435 |
| 26 wk | **7259**, 7260 | 7248, **7250** |
| 52 wk | 7239 (+ 7239 *tumor*) | 7001 |

**11 liver sections from 11 distinct animals**, plus 1 tumor section from the
52-week SBR animal. Bold = downloaded.

Three consequences, all favourable:

1. **Replication units exist.** 11 animals supports an animal-level bootstrap.
   Section 24.1's "fewer than 3 donors = label it a case study" is cleared,
   though 11 is modest and the design is unbalanced.
2. **A progression axis replaces the human aging axis.** 2 / 10 / 26 / 52 weeks.
3. **A built-in dose-response positive control.** SBR vs sham is a high- vs
   low-sender-burden contrast. If the kernel is real, effect size should scale
   with burden. Section 29 objection 6 asks how much true signal the decoy
   control removes; this contrast tests it on real tissue, not just synthetic.

The two samples already pulled (7259 sbr / 7250 sham) are a **matched SBR/sham
pair at the same 26-week timepoint** — the right pair for a first pass.

## 8. Panel verdict — the go/no-go PASSES

Measured directly from `7250_liver_sham_Male_26-U1/cell_feature_matrix.h5`:

| Quantity | Value |
|---|---|
| Cells in section | **237,982** |
| Total features | 13,590 |
| **`Gene Expression` features** | **5,106** |
| Negative Control Probe / Genomic Control / Negative Control Codeword | 40 / 21 / 609 |
| Nonzero entries | 73,525,257 (~309 per cell) |

This is the **Xenium Prime Mouse 5K panel**. The concern that the mouse arm
might use a 250–400 gene panel is resolved: it does not. The paper's own
caveat about "a limited number of genes available in predesigned and custom
panels" does not apply to this deposit.

Per-cell nonzero count (~309) sits at the low end of Section 15.4's 300–800
estimate, so memory budgeting there holds.

### 8.1 Preliminary Tier intersections (CS lead; Bio agent to verify and extend)

| Tier | On-panel | Missing |
|---|---|---|
| A1 arrest | 5/6 | Trp53i3 |
| A2 proliferation | 8/15 | Tyms, Rrm2, Mcm3–Mcm7 |
| A3 envelope/chromatin | 2/3 | Hmgb1 |
| A4 DDR | 6/7 | H2ax |
| **Tier A total** | **~21** | **clears the ≥15 bar** |
| B1 NF-κB proximal | 13/15 | Fosb, Zfp36 |
| B2 IL6/JAK/STAT3 | 9/9 | — |
| C ligands | 11/13 | Il1b, Igfbp7 |
| C receptors | 13/13 | — |

Tier C receptor coverage includes `Ccr2`, `Ackr3`, `Cxcr4` and `Dpp4`, so the
CCL2 / CXCR7 / DPP4 axis is fully covered — good for the Section 9 Tier C
internal control (ligand range ordering).

### 8.2 Two real problems the panel creates

**Problem 1 — Tier E housekeeping controls are almost absent: 1/6 on-panel.**
Only `Tbp`. `Actb`, `Gapdh`, `Rpl13a`, `Rps18`, `Ppia` are all off-panel —
Xenium excludes high-expressors to avoid optical crowding. Section 9 Tier E
requires a flat-kernel technical control.
*Proposed substitute:* the 40 Negative Control Probe + 21 Genomic Control +
609 Negative Control Codeword features. These measure background directly and
are arguably a **better** technical null than housekeeping genes. To be written
up as a justified deviation from Section 9, not an omission.

**Problem 2 — zonation markers are thin.** Periportal 2/5 (`Hal`, `Cps1`
present; `Ass1`, `Sds`, `Alb` missing), pericentral 3/4 (`Glul`, `Cyp2e1`,
`Cyp1a2` present; `Oat` missing). `Alb` absent from a liver panel is the same
optical-crowding exclusion. Section 11 makes zonation *the* liver confound and
requires a continuous zonation covariate, so the score must be rebuilt from
additional on-panel markers. Assigned to the Bio agent as top priority.

## 9. Revised Section 7 ranking

1. **GSE310392 mouse Xenium** — in hand, 5K panel **verified**, 11 animals, 4 timepoints, burden contrast
2. **SenNet portal, human WashU arm** — better panel, access friction; register now, decide by Day 3
3. **10x public Xenium Prime 5K human** — 5K guaranteed, zero friction, no senescence annotation
4. **HTAN DCC mCRC** — some levels open access

---

# Phase 0c — Bulk Fetch Engineering (2026-08-20)

Pulling the remaining 9 liver sections needed a smarter fetcher than `wget`,
because the three files we want (~76–130 MB) sit at unpredictable offsets
inside 12–31 GB archives, and every archive ends with `morphology.ome.tif`
(4–11 GB) plus `transcripts.parquet` (1.9 GB) that Section 17.2 forbids.

**Approach** (`code/fetch_xenium_bundle.sh`): stream the archive, extract only
the three wanted members, and tear the pipeline down the moment they stop
growing. Free egress plus a hard disk quota means bandwidth is the cheap
resource — so we trade it deliberately.

## Three bugs, all worth recording

**1. GEO's accession masking is the last THREE digits.** `GSM9295276` lives
under `GSM9295nnn`, so the prefix is `${GSM%???}nnn`, not `${GSM:0:6}nnn`.
The wrong prefix produced empty filenames and nine silent no-ops.

**2. `pkill -f "$GSM"` killed the caller.** The accession appears in the
invoking script's own command line, so pattern-killing by accession matched
the parent and terminated it (exit 144). Never pattern-kill on a string that
appears in your own argv.

**3. `kill -TERM -$PGID` killed the runner.** After switching to `setsid` +
group kill, the teardown on GSM9295276 took down the whole chain: the sample
script, and the sequential runner above it. `fetch_all.log` never advanced
past sample 1 and sat dead for ~14 minutes before it was noticed. Fixed with a
guard that group-kills only when `setsid` genuinely produced a *different*
process group, falling back to killing the pipeline PID and its children:

```bash
MYPGID=$(ps -o pgid= -p $$ | tr -d ' ')
if [ -n "$PGID" ] && [ "$PGID" != "$MYPGID" ]; then
  kill -TERM -"$PGID"
else
  pkill -P "$PIPE"; kill -TERM "$PIPE"
fi
```

**Diagnostic trap worth noting:** `pgrep -f fetch_all_livers` returns its own
shell, so it reports the runner as ALIVE when the runner is dead. Confirm with
`ps -o pid,etime,cmd` and check elapsed time, or match on the parent PID.

**4. Teardown keyed on one file is fragile.** The original check watched only
`cells.parquet`, assuming it is last of the three. Member ordering is not
guaranteed constant across samples, so the check now watches the *combined*
size of all three targets.

## Fetch status

| Sample | Timepoint | Cond. | Cells | Status |
|---|---|---|---|---|
| 7250 | 26 wk | sham | 237,982 | done |
| 7259 | 26 wk | sbr | 128,030 | done |
| 7239 | 52 wk | sbr | — | done |
| 7361 | 2 wk | sbr | 194,740 | done |
| 7352, 7448, 7450, 7435, 7260, 7248, 7001 | 2–52 wk | both | — | in progress |

Cost per sample so far: ~5–15 min of streaming for ~76–130 MB stored.

---

# Phase 0d — Transcript Assignment Rate (Section 8 Test 1)

The one audit cell left blank through Phases 0–5, because it requires
`transcripts.parquet` (1.85 GB), which Section 17.2 tells you not to download.
It was worth the exception: the plan flags >30% unassigned as evidence of
bleed-through, and bleed-through *manufactures spatial autocorrelation* — the
artefact class this project accuses other methods of failing to control.
Leaving it unmeasured while attacking others' calibration would have been an
obvious asymmetry at review.

Streamed for **7259** (an *admissible* section, not a discarded one); the member
is last in the archive, so it cost a full 17.3 GB traversal.

| Quantity | Value |
|---|---|
| Transcripts total | 129,104,526 |
| **Assigned to a cell** | **113,954,657 (88.27%)** |
| **Unassigned** | **15,149,869 (11.73%)** |
| Q≥20 transcripts | 105,447,674 (81.68% of all) |
| Q≥20 and assigned | 93,086,712 (**88.28%** of Q≥20) |
| Q≥20 unassigned | 12,360,962 (**11.72%**) |
| Q≥20 assigned, nuclear | 31,479,636 (33.82% of assigned) |

**VERDICT: PASS.** 11.72% unassigned against a 30% threshold, and the rate is
unchanged by quality filtering (88.27% raw vs 88.28% at Q≥20), so the assignment
is not being propped up by low-confidence calls. The ~34% nuclear fraction is
normal for Xenium's nucleus-expansion segmentation.

**Why it matters for the argument.** Segmentation bleed-through is one of the
standard explanations for spurious spatial autocorrelation (Section 3). It is
ruled out here. The confounding this project documents — 92% of the naive
amplitude, driven by depth, cell size, local density and receiver cell-type
composition — is *not* an artefact of transcripts leaking between cells. That
strengthens the negative result rather than weakening it.

Reproduce with `python3 code/assignment_rate.py`.
