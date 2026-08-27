# H1 candidate panel — provenance

Source: GEO GSE326743, sample GSM9638040 (`GSM9638040_SPLN07_cell_feature_matrix.h5`),
downloaded 2026-08-27 from
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9638nnn/GSM9638040/suppl/

Panel verified ON THE DATA, not from the series title (Phase 7 §12.1 step 2):
  Gene Expression          5093   (all ENSG-prefixed, 0 duplicate symbols)
  Negative Control Codeword 609
  Negative Control Probe     40
  Genomic Control            21
  Unassigned Codeword       695
  Deprecated Codeword      3291

Cross-checked identical (symmetric difference = 0) against GSM9638044 (SPLN30)
and GSM9638046 (SPLN44). One stock panel across the series.

The 609/40/21 control counts match the pre-designed Xenium Prime 5K control
complement recorded in Phase7 §12.4, which is independent confirmation the
panel is stock Prime 5K rather than a mislabelled v1 design.

NOTE ON THE FREEZE: this is a pre-freeze PANEL SCREEN, explicitly sanctioned by
§12.1 ("Confirm on the data, not the title. Pull one cell_feature_matrix.h5 per
candidate"). No expression values, no cell-level data, and no outcome-bearing
quantity has been looked at. The full H1 acquisition stays gated on §15.
