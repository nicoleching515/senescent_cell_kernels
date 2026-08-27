"""Human spleen cell-type markers, Xenium Prime 5K Human panel (GSE326743, 5,093 genes).

GENERATED FILE -- do not hand-edit. Rebuild with:
    python3 /workspace/code/build_markers_human_spleen.py

TAXONOMY SOURCE: Phase 7 Job A follow-on task 1 (the compartments the coordinator specified:
  red pulp, white pulp T zone / follicle / marginal zone / germinal centre, FDC, sinusoidal vs
  other endothelium, FRC, plasma, T/NK, B, mono/DC, erythroid, megakaryocyte, capsule/trabecular
  smooth muscle). It replaces section 14's lung list, which does not apply to this arm.
MARKER SOURCE: CellMarker 2.0 human, pinned at genesets/cellmarker_pin/ (md5 recorded).
  *** No marker gene here was written from memory. *** Each gene carries its PMIDs in
  genesets/human/markers_spleen_evidence.csv. Filters: >= 1 PMID, on-panel, claimed by <= 3
  labels, capped at the 15 best-evidenced genes, and >= 4 surviving genes per label.

*** NOT VALIDATED ON DATA. *** The mouse equivalent (markers_mouse_liver.py / annotate_pipeline.py)
had its non-specific labels removed on MEASURED expression. H1 is behind the section 15 freeze, so
no H1 expression value has been read and no such check was possible. Re-gate this label set
against measured expression before trusting the fine labels.
"""
MARKERS = {
'Red pulp macrophages':'APOC1 APOE C1QA C1QB CD163 CD33 CD5L CD68 FCGR1A FCGR2A FOLR2 HMOX1 LGMN MARCO MRC1'.split(),
'Monocytes':'C5AR1 FCAR FCN1 S100A12 VCAN'.split(),
'cDC1':'CLEC9A PTTG1 THBD XCR1'.split(),
'cDC2':'C1QA C1QB CD1C CLEC10A FCER1A'.split(),
'pDC':'BCL11A CCDC50 CLEC4C IRF7 IRF8 LILRA4 SPIB TCF4'.split(),
'Follicular B cells':'CD19 CD27 CD79A CD79B FCER2 IL4R MS4A1 TCL1A TNFRSF13B'.split(),
'Marginal zone B cells':'CD180 CR2 EBF1 FCRL4 ITGAE NOTCH2 PAX5 TNFRSF13C'.split(),
'Germinal centre B cells':'AICDA BCL6 MME PTTG1'.split(),
'CD4 T cells':'CD3D CD3E CD4 CD40LG CXCR5 ICOS IL7R LEF1 MAL MYC PDCD1 SOCS3 TCF7 TRAC'.split(),
'CD8 T cells':'CD3D CD3E CD8A CD8B CXCR3 GZMA IL7R ITGA1 LEF1 MAL MYC PRF1 SOCS3 TCF7 TRAC'.split(),
'NK cells':'B3GAT1 CCL4 CCL5 CD160 CD247 CD7 CST7 CX3CR1 CXCR3 CXCR6 FCGR3A KLRD1 KLRF1 NCAM1 NCR1'.split(),
'Fibroblastic reticular cells':'CCL21 CD36 MFAP5 PDPN'.split(),
'Sinusoidal endothelium':'CLEC4M DNASE1L3 F8 FCGR2B FCN2 FCN3 MRC1 OIT3 STAB2'.split(),
'Endothelial cells':'ITGA1 PCDH17 PROX1 VWF'.split(),
'Lymphatic endothelium':'CCL21 FLT4 LYVE1 PDPN'.split(),
'Smooth muscle / capsule':'CNN1 DES MYH11 MYLK NOTCH3 PLN'.split(),
'Pericytes':'ABCC9 ANPEP CD248 CD44 CSPG4 KCNJ8 MCAM NT5E PDGFRB RGS5 SMN1 THY1'.split(),
'Fibroblasts':'COL11A1 COL4A1 CXCL13 IGFBP2 PTGDS'.split(),
'Erythroid cells':'GYPA HBG1 SLC4A1 SPI1 TFRC'.split(),
'Megakaryocytes':'CAVIN2 CCL5 CD33 CD79A CPA3 FLI1 FUT4 GP1BA GP5 IL3RA ITGA2B ITGA6 ITGAV ITGB3 SELP'.split(),
'Neutrophils':'ANPEP CD177 CEACAM8 CSF3R ELANE FCGR3A FCGR3B FPR1 FPR2 FUT4 IL1R2 IL3RA ITGAM MME MPO'.split(),
'Mesothelial cells':'ADAM28 ADAMDEC1 CDH2 CENPF GATA6 GREM2 HPGD IL18 MSLN MYH11 PDPN PRG4 SOX6 UPK3B WT1'.split(),
}
# Labels removed by the >= 4 on-panel marker gate, with what survived:
#   Plasma cells                   JCHAIN,MZB1,XBP1
#   Follicular dendritic cells     CR1,CR2,FCER2
#   Proliferating cells            CD3D,IL7R,TUBB

# Fine labels that share a compartment. annotate_pipeline.py recomputes the
# assignment over each group's union of markers and writes it as cell_type_merged.
MERGE = {
 'B cells':['Follicular B cells', 'Marginal zone B cells', 'Germinal centre B cells'],
 'T/NK cells':['CD4 T cells', 'CD8 T cells', 'NK cells'],
 'Mono/Mac/DC':['Red pulp macrophages', 'Monocytes', 'cDC1', 'cDC2', 'pDC'],
 'Endothelial':['Sinusoidal endothelium', 'Endothelial cells', 'Lymphatic endothelium'],
 'Stromal':['Fibroblastic reticular cells', 'Fibroblasts', 'Pericytes', 'Smooth muscle / capsule'],
}

# No entries: DROP_NONSPECIFIC is an EXPRESSION-based judgement and H1 expression
# is behind the section 15 freeze. Populate it after the freeze, from measured data.
DROP_NONSPECIFIC = {}
