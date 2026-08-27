# HGNC pin — human symbol / alias / previous-symbol / Ensembl authority

`hgnc_symbol_alias_ensembl.2026-08-27.csv.gz` is a column subset
(`hgnc_id, symbol, alias_symbol, prev_symbol, ensembl_gene_id, locus_group, status, gene_group`)
of the HGNC **complete set**, downloaded **2026-08-27** from

    https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

- rows: 45,045 (plus header)
- md5 of the downloaded source file: `2d741e796d5538cc48d3696452237781`

HGNC publishes this endpoint as a rolling "latest" file. The dated monthly archive URLs
(`.../archive/archive/monthly/tsv/hgnc_complete_set_YYYY-MM-01.txt`) returned HTTP 404 on
2026-08-27, so the pin is **by download date + md5**, not by an upstream release tag.

Used by `/workspace/code/build_genesets_human.py` for exactly one purpose: resolving a human
gene symbol to the symbol the GSE326743 Xenium panel actually uses, by going
symbol / alias_symbol / prev_symbol -> ensembl_gene_id -> the panel's own `gene_id` column.
This is why `H2AX` resolves to the panel's legacy `H2AFX` (both ENSG00000188486) rather than
being silently dropped. Every resolution actually applied is logged in
`/workspace/genesets/human/_symbol_resolutions.csv`.
