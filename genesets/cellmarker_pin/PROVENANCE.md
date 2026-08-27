# CellMarker 2.0 pin — human cell-type marker authority

`CellMarker2.0_Human.2026-08-27.csv.gz` is a column subset
(`species, tissue_class, tissue_type, cancer_type, cell_type, cell_name, marker, Symbol, GeneID,
technology_seq, marker_source, PMID, year`) of the CellMarker 2.0 **human** marker table,
downloaded **2026-08-27** from

    http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx
    (mirror: http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx)

- rows: 60,877
- md5 of the downloaded `.xlsx`: `c7a1b764b66cb3a3c16cfac428160f72`
- `tissue_class == 'Spleen'`: 446 rows across 60 `cell_name` values

CellMarker 2.0 publishes no dated release archive at this endpoint, so the pin is
**by download date + md5**.

Every marker row carries its own `PMID`, so each gene entering
`/workspace/code/markers_human_spleen.py` is traceable to a primary publication rather than to
recollection. The per-gene evidence table, including PMID counts, is written to
`/workspace/genesets/human/markers_spleen_evidence.csv` by
`/workspace/code/build_markers_human_spleen.py`.

Reading the `.xlsx` required `openpyxl`, which was **not** in `/workspace/requirements.txt` and was
installed with pip on 2026-08-27 (`openpyxl 3.1.5`). Recorded here because environment pinning is a
tracked concern in this repo. Nothing else in the human-arm build depends on it — the pinned
`.csv.gz` above is what the build script reads.
