# Phase 7 §12.1 — GEO panel-first screen for a human Xenium Prime 5K arm

> **⚠ TWO ROWS IN §2 WERE CORRECTED 2026-08-27.** Both gave **sample** counts where
> **donor** counts were needed, and both were wrong — in opposite directions.
> GSE336890 has *more* donors than stated (20, not 9); GSE335963 has far *fewer*
> (4, not 39 — it is a SuperSeries and only 6 GSM are Xenium at all). Neither error
> changed the H1 selection, but both would have misled anyone using this table as a
> fallback list. See `reports/A3_FALLBACK_SCREEN.md`.


**Date:** 2026-08-27
**Scope:** the §12.1 "Search protocol, GEO only". Panel-first, tissue-agnostic.
**Status:** candidate selected and acquired; **freeze/pre-registration (§15) NOT yet committed** — see the caveat at the end.

## 1. What was searched

Platform `GPL33762` (Xenium In Situ Analyzer, *Homo sapiens*), queried via NCBI
E-utilities on 2026-08-27:

```
esearch db=gds  term="GPL33762[Accession] AND gse[Entry Type]"
```

**132 series** (the doc recorded 117 series / 1,044 samples as of August 2026;
the platform now carries 1,376 samples). Brief SOFT metadata was pulled for all
132 and screened for Prime-5K evidence in series titles, summaries, overall
design, sample titles, sample characteristics and data-processing fields.

- Title-only screening for "5K"/"Prime" returns **5** series — which is why
  §12.1 says to confirm on the data.
- Full-metadata screening returns **19** series with explicit Prime 5K evidence.

**Screen validated against a known answer.** `GSE335761`, the doc's own primary
candidate, is in the 132 and does *not* pass the 5K screen. Its own metadata
gives the reason: *"The Xenium Human Lung Gene Expression panel with a custom
add-on (Design ID 77VFDX; 386 total targets)"* — reproducing §12.1's finding
independently. Worth noting for later: that add-on is *"a custom 100-gene
senescence panel"*, making GSE335761 the only human deposit found with a
purpose-built senescence panel. It still fails A2 as a primary arm (7 Tier B
modules x 30 genes will not come out of 386 targets), but it is a defensible
secondary/validation target.

## 2. The 19 Prime-5K series, screened on donors, tissue and file format

| Series | Human n | Tissue | Files | Verdict |
|---|---|---|---|---|
| **GSE326743** | **7** | **Spleen, normal FFPE** | h5 + parquet | **SELECTED** |
| GSE336890 | 9 slides / **20 donors** | Kidney biopsy (AIN / ATI / reference) | h5 + parquet | **CORRECTED 2026-08-27:** the 9 are Region *slides*, carrying 20 patient specimens (8 AIN / 7 ATI / **only 5 reference**), 358,005 cells. See `reports/A3_FALLBACK_SCREEN.md` |
| GSE335963 | 39 GSM / **4 donors** | Bone marrow (CHIP vs non-CHIP) | h5 + parquet | **CORRECTED 2026-08-27: this is a SuperSeries.** Only **6 of 39** GSM are Xenium (= GSE335962), from **4 donors** (two are repeat sections). **Fails A3 structurally** — 0 of 88 donor x type strata clear the floors at p=2%. **Strike from the runner-up list** |
| GSE335962 | 6 | Bone marrow biopsy (CHIP) | h5 + parquet | Subset of the above |
| GSE311609 | 41 | Primary NSCLC / breast cancer | h5 + parquet | Tumour |
| GSE343063 | 6 | Lung (SCLC) | h5 + parquet + zarr | Tumour |
| GSE328422 | 6 | Lymph node, metastasis stages | h5 + parquet | Tumour |
| GSE313662 | 6 | PDAC, post-trial archival | h5 + parquet | Tumour + therapy |
| GSE326226 | 5 | Bladder cancer | h5 only | Tumour; no cells.parquet listed |
| GSE319763 | 3 | Lung, normal + AML-infiltrated | h5 + parquet + zarr | n=3, mixed |
| GSE315435 | 5 | Fetal kidney | h5 + parquet + zarr | Age axis runs the wrong way |
| GSE322974 | 1 | Bone marrow clot | h5 + parquet | n=1, donor bootstrap gone |
| GSE312420 | 32 | Ileum / rectum (IBD) | parquet | Only 4 samples are the 5K panel; mixed platforms |
| GSE282123 | 7 | Colon (UC / ICI colitis) | tar | Early-access 4,737-gene panel, not stock Prime 5K |
| GSE315411 | 5 | Lung | tar | Hybrid V1+Prime design |
| GSE299193 / GSE299207 | 22 each | Mixed | h5 + zarr | Mixed 5K and Multi-Tissue-and-Cancer panels |
| GSE290468 / GSE301112 | 157 / 51 | Cell lines, benchmarking | — | Not tissue |

## 3. The selected arm: GSE326743

*"Spatial transcriptomics using Xenium, ran over human spleens using Xenium
prime 5k standard panel plus 100 addon genes."*

| Sample | GSM | Age | Sex | Cells |
|---|---|---|---|---|
| SPLN07 | GSM9638040 | 17 | F | 249,420 |
| SPLN14 | GSM9638041 | 37 | F | 329,371 |
| SPLN21 | GSM9638042 | 32 | M | 220,435 |
| SPLN24 | GSM9638043 | 32 | F | 396,173 |
| SPLN30 | GSM9638044 | 57 | M | 366,199 |
| SPLN43 | GSM9638045 | 31 | M | 331,582 |
| SPLN44 | GSM9638046 | 59 | M | 314,413 |

**7 donors, 2,207,593 cells** (M1: 11 sections, 1,834,806 cells).

### Panel verified on the data, not the title (§12.1 step 2)

| Feature class | GSE326743 | M1 / GSE310392 |
|---|---|---|
| Gene Expression | **5,093** (all ENSG, 0 duplicate symbols) | 5,106 (5,097 ENSMUSG) |
| Negative Control Codeword | 609 | 609 |
| Negative Control Probe | 40 | 40 |
| Genomic Control | 21 | 21 |

The panel is **byte-identical across SPLN07 / SPLN30 / SPLN44** (symmetric
difference = 0), so this is one stock panel, not per-sample designs. The
609/40/21 control complement matches the pre-designed Prime 5K figures recorded
in §12.4 — independent confirmation the panel is stock Prime 5K, and the 40
negative control probes that audit test **A7** requires are present.

"5K" is a product name: the measured panel is **5,093 genes, not 5,000**. The
series describes it as the standard panel plus 100 addon genes, which is
arithmetically consistent with ~4,993 predesigned + 100 custom, but that split
has not been confirmed against 10x's published base panel and is not asserted here.

### Why this arm

- **Only normal-tissue Prime 5K human deposit with >=3 donors found in the 132.**
  Every other multi-donor candidate is tumour or disease-defined.
- **A real chronological human ageing axis** (17-59), which §12.1 argues is what
  the senescence field actually cares about and which the mouse arm cannot give.
- **Spleen has a clean, well-defined anatomical axis** — red pulp vs white pulp —
  which supplies the arm-specific covariate that test **A6** demands, playing the
  role liver zonation plays in M1. This is a better structural match than lung's
  airway-to-alveolar axis would have been.
- **Design mirrors M1 exactly**: Prime 5K + 100 custom addon, both arms.
- Per-sample cell counts (220k-396k) sit **above** the mouse range (84k-238k) and overlap it only
  at the top; every H1 sample is larger than the median M1 section. *(Corrected 2026-08-27: this
  read "bracket the mouse range", which they do not — `AUDIT_PHASE8_FACTCHECK.md` M10. Both ranges
  are right; the relation was not. It is a point in H1's favour, not against, but it must be
  stated as "comparable or larger", never as "bracketing", and per-sample depth still has to be
  matched explicitly in the two-arm comparison.)*

### Honest limitations, recorded before the freeze

- **The old end of the age axis is thin.** Ages are 17, 31, 32, 32, 37, 57, 59 —
  five donors under 40, two over 55. This is a sparse continuum, not a
  young-vs-old two-group design. Any age-stratified claim rests on n=2 at the
  top end.
- **Spleen is a third tissue.** M1 is liver, the doc's plan was lung, this is
  spleen. §15's declared deviation "tissue differs between arms, so no cross-arm
  difference can be attributed to species or tissue" stands unchanged and if
  anything is now more emphatic.
- **No disease/treatment contrast at all.** All seven are normal. The within-
  section distance-to-sender machinery does not need a group contrast, but the
  depth-vs-prevalence correlation analysis (rho = 0.94 in M1 SBR) has no direct
  analogue and will have to be reframed as across-donor rather than across-arm.
- **Sex is unbalanced** (4M / 3F) and race is unknown for two donors.

## 4. Acquisition (§12.2)

GEO deposits each Xenium output as its own supplementary file here, so the tar-
streaming machinery in `code/fetch_xenium_bundle.sh` — written because GSE310392
buried the bundle behind an 11.7 GB .rds — is not needed. The **discipline** is
identical: take counts, coordinates and segmentation polygons; skip the rest.
`code/fetch_h1_geo.sh` implements it.

**Taken (28 files, 525 MB):** `cell_feature_matrix.h5`, `cells.parquet.gz`,
`cell_boundaries.parquet.gz`, `annotations.csv.gz`.
**Skipped (21 files):** `transcripts.parquet.gz`, `morphology.ome.tif.gz`,
`nucleus_boundaries.parquet.gz`.

**Not available in this deposit:** `gene_panel.json` and `experiment.xenium`,
both listed in §12.2. Panel membership and the negative-control probe list were
recovered from the `.h5` feature table instead and archived at
`genesets/h1_candidate/`. Run chemistry and segmentation method are recoverable
from the `segmentation_method` column of `cells.parquet` rather than from
`experiment.xenium`.

### Integrity check

All 7 sections: panel = 5,093, `cell_id` sets in `cell_feature_matrix.h5` and
`cells.parquet` **match exactly**, centroids are in microns (ranges ~5-6,500 um),
and `cells.parquet` carries every Tier D nuisance column the plan asks for
(`transcript_counts`, `control_probe_counts`, `cell_area`, `nucleus_area`,
`nucleus_count`, `segmentation_method`).

**Deviation from §12.3.** The plan states cell type labels are "not provided in
usable form". For this deposit they **are**: `annotations.csv.gz` ships four
nested annotation levels (`Level_1` to `Level_4`). This does not replace Job B —
the plan's own pipeline still runs — but it gives an external label set to check
the human-marker annotation against, which is strictly better than eyeballing
cold. Note the depositors' annotations cover **fewer cells than the matrix**
(e.g. SPLN07: 239,167 of 249,420), so they applied a QC filter that must be
characterised before the two label sets are compared.

## 5. Freeze caveat

§0.4 and §15 require the pipeline frozen and `PREREG_PHASE7.md` committed
**before** the human data is downloaded. The data is now on disk ahead of that.

What has actually been looked at is **panel membership, file structure, cell
counts and coordinate ranges** — the §12.1 step-2 screen plus a structural
integrity check. **No expression values, no senescence score, no cell-type
assignment and no outcome-bearing quantity has been computed.** The pre-
registration is still writable without contamination, but it must be committed
before Job B or any A3/A5 test runs against these files.
