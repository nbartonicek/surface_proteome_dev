#!/usr/bin/env Rscript

# Builds the Tier 1 + Tier 2 "broad AML biology" gene signatures requested
# for the per-sample report (22.patient_qc_report.Rmd), beyond the Van Galen
# lineage signatures already in annotation/signatures/. Saves everything into
# annotation/signatures/signatures_broad/, mirroring the van_galen convention
# (one consolidated named-list RDS + a documentation CSV of sources/confidence).
#
# Item 1 (differentiation-state composition) is intentionally NOT a gene
# signature here - it's the existing BoneMarrowMap projection
# (predicted_CellType_Broad), already wired into the report.
#
# CONFIDENCE LEVELS (see aml_broad_signature_sources.csv for the per-signature
# column): "high" = gene list confirmed against a citable, exact source.
# "medium" = extracted from an author-provided public resource (e.g. a GitHub
# GMT file) via automated fetch, not manually cross-checked against the
# original supplementary PDF. "literature-derived" = built directly from gene
# symbols given in the request/literature review, not itself a formally
# validated classifier (per the "practical deployment notes" caveat) - true
# for most Tier 1/2 items other than LSC17. "pending" = no gene/protein list
# could be sourced yet (Jayavelu 27-protein Mito-AML classifier, Lasry iScore
# adult/pediatric) - these are NOT fabricated, and are absent from the RDS.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
out_dir <- file.path(proj, "annotation", "signatures", "signatures_broad")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sigs <- list()
sources <- list()

add_sig <- function(name, genes, tier, source, confidence, notes = "") {
  sigs[[name]] <<- unique(genes)
  sources[[name]] <<- tibble(signature = name, tier = tier, n_genes = length(unique(genes)),
                             source = source, confidence = confidence, notes = notes)
}

# ------------------------------------------------------------------
# Tier 1
# ------------------------------------------------------------------

# #2 LSC17 (Ng et al 2016, Nature) - exact 17 genes + published regression
# coefficients, cross-checked against the reproduction in the ALFA study
# (Blood Advances 2023, PMC10410128). Coefficients are used for the dedicated
# weighted LSC17_score in the report (NOT just AUCell) - see
# score_lsc17_weighted() in 22.patient_qc_report.Rmd.
lsc17_genes <- c("DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B", "MMRN1",
                 "DPYSL3", "KIAA0125", "CDK6", "CPXM1", "SOCS2", "SMIM24",
                 "EMP1", "NGFRAP1", "CD34", "AKR1C3", "GPR56")
lsc17_coefs <- c(DNMT3B = 0.0874, ZBTB46 = -0.0347, NYNRIN = 0.00865, ARHGAP22 = -0.0138,
                 LAPTM4B = 0.00582, MMRN1 = 0.0258, DPYSL3 = 0.0284, KIAA0125 = 0.0196,
                 CDK6 = -0.0704, CPXM1 = -0.0258, SOCS2 = 0.0271, SMIM24 = -0.0226,
                 EMP1 = 0.0146, NGFRAP1 = 0.0465, CD34 = 0.0338, AKR1C3 = -0.0402, GPR56 = 0.0501)
add_sig("LSC17", lsc17_genes, "Tier1_2_LSC_stemness",
        "Ng et al 2016 Nature (nature20598); coefficients cross-checked vs Blood Advances 2023 ALFA study (PMC10410128)",
        "high",
        "Older gene aliases (NGFRAP1=BEX3, KIAA0125=FAM30A, SMIM24=C11orf21) kept as-published; add alias fallback if not found in a newer reference annotation.")

# #2 Quiescent LSPC (Zeng et al 2022 Nature Medicine cellular-hierarchy paper).
# Extracted from the authors' own public GMT (github.com/andygxzeng/AMLHierarchies,
# Data/AMLCellType_Genesets.gmt), LSPC-Quiescent geneset - NOT manually
# cross-checked against the original supplementary table, hence "medium".
zeng_quiescent_lspc <- c("GAS5", "RPLP0", "LRRC75A-AS1", "SPINK2", "FAM30A", "RPLP0P2",
                         "ANGPT1", "RPS18P9", "DSE", "SNHG8", "CD34", "PTPRCAP", "GUCY1A3",
                         "HOPX", "EGFL7", "CD52", "EBPL", "ITM2A", "SOX4", "SELENOP",
                         "RPL31P11", "SMIM24", "TFDP2", "SNHG7", "BCL11A", "MSI2", "H1F0",
                         "NPR3", "CD69", "TFPI", "SOCS2", "MPG", "MAP7", "AKR1C3", "SMYD3",
                         "C1orf186", "GNA15", "BCAT1", "SCHIP1", "SSBP2", "GATA2", "CD200",
                         "STON2", "CD3D", "TSC22D1", "MMRN1", "SLC35F2", "NRIP1", "MAP1A",
                         "CEP70", "GCSAML", "OAF", "NPDC1", "TPSB2", "MECOM", "HEMGN",
                         "CPA3", "TPSD1", "XIRP2", "GPM6A", "APOOL")
add_sig("Zeng2022_QuiescentLSPC", zeng_quiescent_lspc, "Tier1_2_LSC_stemness",
        "Zeng et al 2022 Nature Medicine; gene list via github.com/andygxzeng/AMLHierarchies Data/AMLCellType_Genesets.gmt (LSPC-Quiescent)",
        "medium", "Auto-extracted from public GMT, not manually verified against original PDF supplement.")

# #3 Stetson et al 2021 relapse-resistant LIC / stem-survival module + its
# metabolic companion module (same paper). "NDUFA/NDUFB/NDUFS family" and
# "COX family" resolved to the actual complex I / complex IV subunit genes
# already curated in MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION (see below),
# rather than re-deriving that family list from scratch.
stetson_stem_survival <- c("BCL2", "MCL1", "CXCR4", "IRF8", "GADD45A", "CTNNB1", "LEF1")
add_sig("Stetson2021_StemSurvival", stetson_stem_survival, "Tier1_3_relapse_LIC",
        "Stetson et al 2021", "literature-derived", "Genes as given in request; not independently re-derived from the paper's supplement.")

hallmark_oxphos <- c("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2",
                     "AIFM1","ALAS1","ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D",
                     "ATP5F1E","ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5PD","ATP5ME","ATP5PF",
                     "ATP5MF","ATP5MG","ATP5PO","ATP6AP1","ATP6V0B","ATP6V0C","ATP6V0E1","ATP6V1C1",
                     "ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA","BDH2",
                     "MPC1","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B",
                     "COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A",
                     "CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT","DLD","DLST","ECH1","ECHS1",
                     "ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4",
                     "GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A",
                     "IDH3B","IDH3G","IMMT","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2",
                     "MFN2","MGST3","MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15",
                     "MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3","NDUFA4",
                     "NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2",
                     "NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2",
                     "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1",
                     "NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX",
                     "PDK4","PDP1","PHB2","PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1",
                     "RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11","SLC25A12","SLC25A20","SLC25A3",
                     "SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1",
                     "TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70",
                     "UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ",
                     "VDAC1","VDAC2","VDAC3")
stetson_metabolic_companion <- c("CPT1C", "CPT2", "CD36", "SLC7A11",
                                 grep("^NDUF|^COX", hallmark_oxphos, value = TRUE))
add_sig("Stetson2021_MetabolicCompanion", stetson_metabolic_companion, "Tier1_3_relapse_LIC",
        "Stetson et al 2021 (named genes) + NDUF/COX family members from MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "literature-derived", "NDUF/COX 'family' resolved via the standard Hallmark gene set, not re-derived from Stetson's own supplement.")

# #4 Mito-AML / OXPHOS burden. Jayavelu's exact 27-protein classifier could
# not be sourced (paywalled supplement) - flagged "pending", NOT fabricated.
# Using MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION as the RNA-level first-pass
# proxy, exactly as the request specifies, plus a narrower MRPL/MRPS/NDUF-only
# subset for the Stratmann et al 2023 mitoribosomal/respiratory-chain angle.
add_sig("Hallmark_OXPHOS_MitoAML_proxy", hallmark_oxphos, "Tier1_4_mito_oxphos",
        "MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION (gsea-msigdb.org), used per request as first-pass RNA proxy for Jayavelu et al 2022 Mito-AML",
        "high_as_geneset_low_as_MitoAML_proxy",
        "This is the standard, exact Hallmark gene set (high confidence as a gene list) but only a proxy for Mito-AML (RNA vs protein divergence flagged in the report - see Jayavelu et al 2022). The actual 27-protein classifier is NOT sourced here (paywalled); treat this column as lower-confidence than the 4 fixed classifiers until that list is added.")

stratmann_mitoribosomal_ndufs <- grep("^MRPL|^MRPS|^NDUF", hallmark_oxphos, value = TRUE)
add_sig("Stratmann2023_MitoRibosomal_NDUF", stratmann_mitoribosomal_ndufs, "Tier1_4_mito_oxphos",
        "Stratmann et al 2023 (MRPL/MRPS mitoribosomal + NDUF respiratory-chain family); gene family resolved via MSigDB Hallmark OXPHOS membership",
        "literature-derived", "")

# #5 Interferon-response / CD32A-SNRNP200 axis (Knorr et al 2023 - this
# project's own core biology). FCGR2A and SNRNP200 are tracked as standalone
# expression columns in the report, not blended into this module score.
knorr_ifn_axis <- c("IFITM2", "IFITM3", "STAT1", "ISG15", "MX1", "OAS1", "IFI27", "IFI44L")
add_sig("Knorr2023_IFN_axis", knorr_ifn_axis, "Tier1_5_IFN_axis",
        "Knorr et al 2023 (IFITM2/IFITM3 explicitly named + standard ISG panel per request)",
        "literature-derived", "FCGR2A (CD32A) and SNRNP200 tracked as separate standalone columns, not part of this gene set - see companion_genes in the Rmd.")

# #6 Immunotherapy surface target panel (compiled across de Boer 2018,
# Bordeleau 2024, Perna 2017, Kohnke 2022, Knorr 2023). Split into the three
# sub-groups from the request; the report also cross-references this against
# the actual CITE_DSB antibody panel (annotation/adt_tags.csv) for a
# CITE-measured companion score wherever an antibody exists.
target_canonical <- c("CD33", "IL3RA", "CLEC12A", "CD70", "CD47", "CD44", "CD96", "HAVCR2",
                      "IL2RA", "CD99", "IL1RAP", "FLT3", "ADGRE2", "CCR1", "LILRB2", "ITGB7",
                      "PTPRJ", "ITGA4")
target_subtype_linked <- c("SEMA4D", "ENG", "NCAM1", "VSIR", "ADGRG1")
target_noncanonical <- c("SNRNP200", "CD180", "MRC1")
add_sig("ImmunotherapyTargets_canonical", target_canonical, "Tier1_6_surface_targets",
        "Compiled: de Boer 2018, Bordeleau 2024, Perna 2017, Kohnke 2022, Knorr 2023", "literature-derived", "")
add_sig("ImmunotherapyTargets_subtype_linked", target_subtype_linked, "Tier1_6_surface_targets",
        "Bordeleau 2024 (subtype-linked targets: RUNX1-mutant/complex karyotype-MECOM/KMT2A-r monocytic/NPM1-NK-triple-mutant)",
        "literature-derived", "")
add_sig("ImmunotherapyTargets_noncanonical", target_noncanonical, "Tier1_6_surface_targets",
        "Compiled across the same surface-target literature", "literature-derived",
        "SNRNP200 is transcript-only here (also tracked standalone per Knorr axis above); CD180/MRC1 are genes, MRC1(CD206) also has a CITE antibody ('MMR').")

# #7 Genotype-linked transcriptional programs. NPM1 180-gene classifier
# (Lilljebjorn et al) is ALREADY in this project as NPM1_classI/NPM1_classII
# inside annotation/signatures/van_galen_signatures_top50.rds (built earlier
# from the Uckelmann et al 2025 signature file) - not duplicated here.
# TP53/complex-karyotype erythroid program and RUNX1 early-lymphoid/pDC
# program are NOT new gene sets - both are scored directly against the
# existing BoneMarrowMap differentiation-state composition (item #1) in the
# report (erythroid-state fraction; Early Lymphoid + pDC fraction), per the
# request's own framing ("score against your differentiation-state output").
npm1_hox_program <- c("HOXA9", "HOXA10", "HOXB3", "MEIS1")
add_sig("NPM1_HOX_program", npm1_hox_program, "Tier1_7_genotype_programs",
        "Zeng 2025, Lilljebjorn 2025", "literature-derived", "")

# ------------------------------------------------------------------
# Tier 2
# ------------------------------------------------------------------

# #8 iScore (Lasry et al 2023) - exact 38-gene adult / 11-gene pediatric list
# could not be sourced (supplementary table not accessible). NOT fabricated -
# deliberately absent from the RDS; the report shows a "pending" placeholder.

# #9 Immune microenvironment composition (Lasry 2023, Mikami 2025). HLA-DRA/
# HLA-DRB1 tracked as standalone companion columns (the M-MDSC-like call is
# "LOW HLA-DR", which a positive-gene-set AUCell score can't represent).
add_sig("ImmuneMicroenvironment_AtypicalB", c("ITGAX", "FCRL3", "FCRL5", "IRF8", "CD72"),
        "Tier2_9_immune_microenvironment", "Lasry et al 2023, Mikami et al 2025", "literature-derived", "")
add_sig("ImmuneMicroenvironment_ExhaustedCD8T", c("GZMK", "PDCD1", "TIGIT", "TOX"),
        "Tier2_9_immune_microenvironment", "Lasry et al 2023, Mikami et al 2025", "literature-derived", "")
add_sig("ImmuneMicroenvironment_MDSClike", c("VEGFA", "ORM1", "CD82", "MSRB1"),
        "Tier2_9_immune_microenvironment", "Lasry et al 2023, Mikami et al 2025", "literature-derived",
        "Diagnostic signature is 'high on these 4 genes AND low HLA-DR' - HLA-DRA/HLA-DRB1 tracked as separate companion columns in the report, not part of this gene set.")

# #10 AP-1 stress/leukemic program (Velten 2021, Zhai 2022)
add_sig("AP1_StressProgram", c("FOS", "JUN", "JUNB", "JUND", "FOSB", "EGR1"),
        "Tier2_10_AP1_program", "Velten et al 2021, Zhai et al 2022", "literature-derived", "")

# #11 Glycocalyx / physical immune-evasion barrier (Chung et al 2026). Split
# into 3 axes per the paper's own framing (core glycocalyx, MHC-I pathway,
# LILRB evasion) for more interpretable separate scores. Sialyltransferase and
# MHC-I antigen-presentation gene lists are standard/stable HGNC gene-family
# members, not re-derived from Chung et al's own supplement (2026 paper, not
# independently retrievable here) - flagged "literature-derived", not "high".
glycocalyx_core <- c("SPN", "MUC1", "C1GALT1", "C1GALT1C1",
                     "ST3GAL1", "ST3GAL2", "ST3GAL3", "ST3GAL4", "ST3GAL5", "ST3GAL6",
                     "ST6GAL1", "ST6GAL2", "ST6GALNAC1", "ST6GALNAC2", "ST6GALNAC3",
                     "ST6GALNAC4", "ST8SIA1", "ST8SIA4", "ST8SIA6")
mhci_antigen_presentation <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9")
lilrb_evasion <- c("LILRB1", "LILRB2")
add_sig("Glycocalyx_core", glycocalyx_core, "Tier2_11_glycocalyx",
        "Chung et al 2026 (SPN/MUC1/C1GALT1/C1GALT1C1 as named) + standard sialyltransferase gene family",
        "literature-derived", "Sialyltransferase family membership is standard HGNC nomenclature, not re-derived from this specific paper's supplement.")
add_sig("MHCI_AntigenPresentation", mhci_antigen_presentation, "Tier2_11_glycocalyx",
        "Chung et al 2026 (MHC-I pathway, first evasion axis)", "literature-derived", "")
add_sig("LILRB_Evasion", lilrb_evasion, "Tier2_11_glycocalyx",
        "Chung et al 2026 (LILRB1/LILRB2, second evasion axis)", "literature-derived", "")

# #12 Cell-cycle phase composition: deliberately NOT built here. Scored
# directly in the report via Seurat::CellCycleScoring() using Seurat's own
# built-in Tirosh S/G2M gene lists (cc.genes.updated.2019) - a base package
# resource, not a project-specific signature file.

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

saveRDS(sigs, file.path(out_dir, "aml_broad_signatures.rds"))

source_tbl <- bind_rows(sources)
write.csv(source_tbl, file.path(out_dir, "aml_broad_signature_sources.csv"), row.names = FALSE)

# LSC17 coefficients saved separately - used for the weighted LSC17_score,
# not for AUCell (which treats all genes as equally-weighted positive markers).
saveRDS(lsc17_coefs, file.path(out_dir, "LSC17_coefficients.rds"))

cat("Saved", length(sigs), "gene signatures to", out_dir, "\n")
cat("\nPer-signature gene counts:\n")
print(sapply(sigs, length))
cat("\nPending (not built - no source list found):\n")
cat("  - Lasry2023_iScore_adult (38 genes)\n")
cat("  - Lasry2023_iScore_pediatric (11 genes)\n")
cat("  - Jayavelu2022_MitoAML_27protein (using Hallmark_OXPHOS_MitoAML_proxy as RNA proxy instead)\n")
cat("\nDone.\n")
