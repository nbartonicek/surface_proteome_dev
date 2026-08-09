#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Biotin report - step 02 of 4
#
# Splits HSC/MPP cells into biotin-high and biotin-low on CITE_DSB, then runs FindMarkers, fgsea on the ranking, pseudobulk DESeq2, per-protein CITE DE, glycosylation by fine cell type, surface-load-adjusted residuals, and cell-cycle proportions.
#
# Frozen for the lab archive 2026-08-04 from scripts/27a.hsc_mpp_biotin_split_analysis.R (mtime 2026-07-22).
# md5 of the original: 92093b7d4293adfbf73f48f92931b5a6
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# HSC MPP cells (predicted_CellType_Broad), split into biotin-high vs
# biotin-low (CITE_DSB), then compared:
#   1. Seurat FindMarkers (single-cell RNA, Wilcoxon)
#   1b. fgsea/MSigDB Hallmark pathway analysis on the FindMarkers ranking
#   2. Pseudobulk RNA DE (sum counts -> DESeq2, fixed dispersion)
#   3. Differential CITE-seq binding (per-protein Wilcoxon on DSB values)
#   3b. Glycosylation (biotin) by fine cell type, with cell counts
#   3c. CITE-seq residuals vs glycosylation (biotin), surface-load-adjusted
#   4. Cell-cycle phase (Seurat CellCycleScoring) proportions
#
# Same Seurat object as 27.removeme_quick_biotin_test.R, and the same
# CITE-DE/cell-cycle conventions as 22.patient_qc_report.Rmd
# (run_cite_wilcoxon_de's raw-value rank test with floored-at-0 display
# means; score_cell_cycle's default Tirosh S/G2M gene lists) - not sourced
# from that Rmd directly since it's parameterized per-donor/report, but
# deliberately kept identical in method.
#
# IMPORTANT SCOPE NOTE for sections 2 and 3c: this script is restricted to
# ONE sample (SAMPLE_NAME below). Both the pseudobulk DE (1 profile per
# biotin group) and the CITE/glyco residual regressions are therefore
# single-patient, cell-level analyses - there is no second patient level to
# estimate a real dispersion from, or to put a "(1|patient)" random effect
# on. That is a different, larger analysis (comparing multiple patients,
# e.g. AML vs the normal-02/normal-CD34-01 samples this script deliberately
# excludes) than "add a section to the WEI-only script" - see the response
# this was requested in for what would need to change to do that properly.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(patchwork)
})

# ============================================================
# Settings
# ============================================================

SEURAT_FILE  <- "../results_nf/260717_VH01624_477_222KG22NX/rds/11_final/LK3-GEX/rds/LK3_final.rds"
# LK3-GEX pools 3 samples: "normal-02", "normal-CD34-01", "WEI21_26-17_NK" -
# this analysis is about WEI21_26-17_NK specifically, not the normal
# controls (which contribute almost no HSC MPP cells anyway: 13 + 6 vs
# 1369, but should be excluded outright, not just outnumbered).
SAMPLE_NAME  <- "WEI21_26-17_NK"
CELLTYPE_COL <- "predicted_CellType_Broad"
CELLTYPE_VAL <- "HSC MPP"
BIOTIN_MARKER <- "biotin"

MIN_N_DE     <- 20    # min cells per group for any DE test, same as 22.patient_qc_report.Rmd's min_cells_de
BCV          <- 0.4   # fixed dispersion for the pseudobulk edgeR exact test, same as the report's default
FDR_CUTOFF   <- 0.05
LOGFC_CUTOFF <- 1

out_dir <- "../results/hsc_mpp_biotin_split"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load object, join Seurat v5 layers if needed (same guard as
# 22.patient_qc_report.Rmd uses on object load)
# ============================================================

message("Loading: ", SEURAT_FILE)
seu <- readRDS(SEURAT_FILE)

Seurat::DefaultAssay(seu) <- "RNA"
if (inherits(seu[["RNA"]], "Assay5")) seu <- SeuratObject::JoinLayers(seu, assay = "RNA")
if ("CITE_DSB" %in% names(seu@assays) && inherits(seu[["CITE_DSB"]], "Assay5")) {
  seu <- SeuratObject::JoinLayers(seu, assay = "CITE_DSB")
}

stopifnot(CELLTYPE_COL %in% colnames(seu@meta.data))
stopifnot("mapping_error_QC" %in% colnames(seu@meta.data))
stopifnot("CITE_DSB" %in% names(seu@assays))
stopifnot(BIOTIN_MARKER %in% rownames(seu[["CITE_DSB"]]))
stopifnot("sample_name" %in% colnames(seu@meta.data))
stopifnot(SAMPLE_NAME %in% seu$sample_name)

seu <- subset(seu, subset = sample_name == SAMPLE_NAME)
message("Restricted to sample_name == '", SAMPLE_NAME, "': ", ncol(seu), " cells")

# Cell-cycle phase, scored on the whole object (not just HSC MPP) before
# subsetting - it's a per-cell score off fixed marker genes (Seurat's own
# Tirosh S/G2M lists), so it doesn't need to be, and shouldn't be,
# recalibrated to a single cell type's expression range.
s_use   <- intersect(Seurat::cc.genes.updated.2019$s.genes, rownames(seu))
g2m_use <- intersect(Seurat::cc.genes.updated.2019$g2m.genes, rownames(seu))
message("Cell-cycle genes found: ", length(s_use), " S-phase, ", length(g2m_use), " G2M")
seu <- Seurat::CellCycleScoring(seu, s.features = s_use, g2m.features = g2m_use)

# ============================================================
# RNA+Protein (WNN) UMAP on the whole object, same construction as
# 22.patient_qc_report.Rmd's compute_cite_umap()/compute_wnn_umap(): a
# CITE_DSB-only PCA, then FindMultiModalNeighbors(pca, pca_cite) -> wsnn
# graph -> umap_wnn (no FindClusters step here - wnn_clusters themselves
# aren't needed for this script, unlike the report's cluster-flow
# alluvial). Computed on the full object, not just HSC MPP, so the
# embedding reflects real cell-type structure rather than the geometry of
# one ~1000-cell subset - subset() below carries this reduction (and every
# other existing one) through to `hsc` automatically.
# ============================================================

message("\nComputing RNA+Protein (WNN) UMAP on the full object...")

Seurat::DefaultAssay(seu) <- "CITE_DSB"
seu <- Seurat::ScaleData(seu, assay = "CITE_DSB", verbose = FALSE)
n_cite_features <- nrow(seu[["CITE_DSB"]])
npcs_cite <- max(2, min(20, n_cite_features - 1, ncol(seu) - 1))
seu <- Seurat::RunPCA(seu, assay = "CITE_DSB", reduction.name = "pca_cite",
                      npcs = npcs_cite, features = rownames(seu[["CITE_DSB"]]), verbose = FALSE)

Seurat::DefaultAssay(seu) <- "RNA"
npcs_rna <- min(30, ncol(Seurat::Embeddings(seu, "pca")))
seu <- Seurat::FindMultiModalNeighbors(
  seu, reduction.list = list("pca", "pca_cite"),
  dims.list = list(seq_len(npcs_rna), seq_len(npcs_cite)),
  modality.weight.name = c("RNA.weight", "DSB.weight"), verbose = FALSE
)
seu <- Seurat::RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "umap_wnn",
                       reduction.key = "wnnUMAP_", verbose = FALSE)

# ============================================================
# Subset to HSC MPP, split by biotin into high/low (median split within
# this cell type - "high"/"low" relative to HSC MPP's own biotin
# distribution, not the whole object's). Ties at exactly the median land
# in the low group (ifelse uses > , not >=).
#
# predicted_CellType_Broad == "HSC MPP" is confirmed equal to the union of
# the fine-grained BoneMarrowMap calls HSC (216) + MPP-MkEry (589) +
# MPP-MyLy (564) = 1369, restricted to mapping_error_QC == "Pass" (all 1369
# already pass QC here, but the condition is kept explicit rather than
# assumed). Those 3 fine states are exactly the kind of substates a real
# quiescent-vs-cycling/early-vs-late sub-split would need to consider before
# trusting any AML-vs-normal comparison - not done here (single-sample,
# single-cell-type scope), flagging it rather than silently ignoring it.
# ============================================================

hsc <- subset(seu, subset = predicted_CellType_Broad == CELLTYPE_VAL & mapping_error_QC == "Pass")
message("HSC MPP cells (Pass QC): ", ncol(hsc))

# FetchData has no assay= argument - it reads from whichever assay is
# currently DefaultAssay(), so switch to CITE_DSB just for this fetch.
Seurat::DefaultAssay(hsc) <- "CITE_DSB"
biotin_vals <- Seurat::FetchData(hsc, vars = BIOTIN_MARKER, layer = "data")[[1]]
biotin_median <- median(biotin_vals)
hsc$biotin_value <- biotin_vals
hsc$biotin_group <- ifelse(biotin_vals > biotin_median, "biotin_high", "biotin_low")

cells_high <- colnames(hsc)[hsc$biotin_group == "biotin_high"]
cells_low  <- colnames(hsc)[hsc$biotin_group == "biotin_low"]
message("Biotin median = ", round(biotin_median, 3), " -> ",
        length(cells_high), " biotin_high / ", length(cells_low), " biotin_low")

saveRDS(hsc, file.path(out_dir, "hsc_mpp_biotin_split_seurat.rds"))

# ============================================================
# 1. Seurat FindMarkers - single-cell RNA DE, biotin_high vs biotin_low
# ============================================================

message("\n=== 1. Seurat FindMarkers (RNA, Wilcoxon) ===")

if (length(cells_high) < MIN_N_DE || length(cells_low) < MIN_N_DE) {
  message("Fewer than ", MIN_N_DE, " cells in one biotin group - skipping FindMarkers.")
  findmarkers_res <- NULL
} else {
  Seurat::DefaultAssay(hsc) <- "RNA"
  Seurat::Idents(hsc) <- "biotin_group"
  findmarkers_res <- Seurat::FindMarkers(
    hsc, ident.1 = "biotin_high", ident.2 = "biotin_low",
    assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1
  ) %>%
    rownames_to_column("gene") %>% as_tibble() %>% arrange(p_val_adj)

  write_csv(findmarkers_res, file.path(out_dir, "1_findmarkers_biotin_high_vs_low.csv"))
  sig_fm <- findmarkers_res %>% filter(p_val_adj < FDR_CUTOFF, abs(avg_log2FC) >= LOGFC_CUTOFF)
  message(nrow(sig_fm), " genes at FDR < ", FDR_CUTOFF, ", |avg_log2FC| >= ", LOGFC_CUTOFF)
}

# ============================================================
# 1b. Pathway analysis (fgsea, MSigDB Hallmark + KEGG) on the FindMarkers
#     ranking. Ranked by signed significance
#     (sign(avg_log2FC) * -log10(p_val)), not avg_log2FC alone - standard
#     fgsea practice, since a gene can have a large fold-change on very few
#     cells (noisy) vs a smaller but highly consistent one; the signed
#     p-value rewards the latter appropriately.
#
#     Hallmark alone is broad/summary-level (50 curated sets spanning many
#     processes at once) - KEGG_LEGACY added for actual named metabolic
#     pathways (glycolysis, TCA cycle, OXPHOS, amino/nucleotide sugar
#     metabolism, etc.). Uses collection = "C2", subcollection =
#     "CP:KEGG_LEGACY" (186 sets, the classic KEGG pathway set, license-
#     cleared) rather than "CP:KEGG_MEDICUS" (658 sets - a newer, more
#     complex network-based reorganization that also covers drugs/disease,
#     not the recognizable named metabolic pathways this was asked for).
#     msigdbr's `category`/`subcategory` args were renamed
#     `collection`/`subcollection` in the installed version (26.1.0) -
#     using the current names directly, not the deprecated ones.
# ============================================================

message("\n=== 1b. Pathway analysis (fgsea, MSigDB Hallmark + KEGG) ===")

HAS_MSIGDBR <- requireNamespace("msigdbr", quietly = TRUE)
HAS_FGSEA   <- requireNamespace("fgsea", quietly = TRUE)

fgsea_res <- NULL
if (is.null(findmarkers_res) || nrow(findmarkers_res) == 0) {
  message("No FindMarkers result - skipping pathway analysis.")
} else if (!HAS_MSIGDBR || !HAS_FGSEA) {
  message("msigdbr (", HAS_MSIGDBR, ") and/or fgsea (", HAS_FGSEA,
          ") not installed - skipping pathway analysis.")
} else {
  gene_sets <- bind_rows(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>% mutate(source = "HALLMARK"),
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
      mutate(source = "KEGG")
  )
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  pathway_source <- gene_sets %>% distinct(gs_name, source)

  ranks <- findmarkers_res %>%
    mutate(rank_stat = sign(avg_log2FC) * -log10(pmax(p_val, .Machine$double.xmin))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(desc(rank_stat))
  rank_vec <- setNames(ranks$rank_stat, ranks$gene)

  fgsea_res <- fgsea::fgsea(pathways = pathways, stats = rank_vec, eps = 0) %>%
    as_tibble() %>%
    left_join(pathway_source, by = c("pathway" = "gs_name")) %>%
    arrange(pval)

  write_csv(
    fgsea_res %>% mutate(leadingEdge = purrr::map_chr(leadingEdge, paste, collapse = ";")),
    file.path(out_dir, "1b_fgsea_hallmark_kegg_biotin_high_vs_low.csv")
  )
  sig_path <- fgsea_res %>% filter(padj < FDR_CUTOFF)
  message(nrow(sig_path), " pathways at FDR < ", FDR_CUTOFF, " (",
          sum(sig_path$source == "HALLMARK"), " Hallmark, ", sum(sig_path$source == "KEGG"), " KEGG)")

  if (nrow(fgsea_res) > 0) {
    # Top N per source, not top N overall - Hallmark's larger effect sizes
    # would otherwise crowd every KEGG pathway out of the plot.
    top_path <- fgsea_res %>% group_by(source) %>% slice_min(padj, n = 10, with_ties = FALSE) %>%
      ungroup() %>% mutate(pathway = str_remove(pathway, "^HALLMARK_|^KEGG_"))
    p_path <- ggplot(top_path, aes(x = NES, y = reorder(pathway, NES), fill = padj < FDR_CUTOFF)) +
      geom_col() +
      facet_wrap(~ source, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c(`TRUE` = "firebrick3", `FALSE` = "grey60"), guide = "none") +
      labs(title = paste0(CELLTYPE_VAL, " - top Hallmark/KEGG pathways (fgsea), biotin high vs low"),
           x = "NES (positive = higher in biotin_high)", y = "")
    ggsave(file.path(out_dir, "1b_fgsea_hallmark_kegg_biotin_high_vs_low.pdf"), p_path, width = 8, height = 11)
  }
}

# ============================================================
# 2. Pseudobulk RNA DE - DESeq2, fixed dispersion (one pseudobulk profile
#    per biotin group, no per-patient replicates here, so a real dispersion
#    can't be estimated by ANY tool from n=1 per group - this is DESeq2's
#    own counterpart of the previous edgeR::exactTest approach
#    (skip estimateDispersions(), plug in the same fixed bcv, Wald test
#    instead of exact test), not a fix for the underlying replication
#    problem. Output columns kept as logFC/PValue/FDR (renamed from
#    DESeq2's log2FoldChange/pvalue/padj) so the volcano/violin code below
#    didn't need to change when the method did.
# ============================================================

message("\n=== 2. Pseudobulk RNA DE (DESeq2, fixed dispersion = ", BCV^2, ") ===")

make_pseudobulk_from_cells <- function(count_mat, cells) {
  cells <- intersect(cells, colnames(count_mat))
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}

run_deseq2_fixed_dispersion <- function(count_mat, group_vec, bcv = BCV) {
  group_vec <- factor(group_vec, levels = c("reference", "test"))
  keep <- edgeR::filterByExpr(edgeR::DGEList(counts = count_mat, group = group_vec))
  count_mat <- round(count_mat[keep, , drop = FALSE])

  coldata <- data.frame(group = group_vec, row.names = colnames(count_mat))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~group)
  dds <- DESeq2::estimateSizeFactors(dds)
  DESeq2::dispersions(dds) <- rep(bcv^2, nrow(dds))
  dds <- DESeq2::nbinomWaldTest(dds)
  DESeq2::results(dds, contrast = c("group", "test", "reference")) %>%
    as.data.frame() %>% rownames_to_column("gene") %>% as_tibble() %>%
    transmute(gene, logFC = log2FoldChange, logCPM = log2(baseMean + 1), PValue = pvalue, FDR = padj) %>%
    arrange(FDR)
}

if (length(cells_high) < MIN_N_DE || length(cells_low) < MIN_N_DE) {
  message("Fewer than ", MIN_N_DE, " cells in one biotin group - skipping pseudobulk DE.")
  pseudobulk_res <- NULL
} else {
  rna_counts <- Seurat::GetAssayData(hsc, assay = "RNA", layer = "counts")
  pb <- cbind(
    reference = make_pseudobulk_from_cells(rna_counts, cells_low),
    test      = make_pseudobulk_from_cells(rna_counts, cells_high)
  )
  pseudobulk_res <- tryCatch(
    run_deseq2_fixed_dispersion(pb, group_vec = c("reference", "test"), bcv = BCV) %>%
      mutate(n_high = length(cells_high), n_low = length(cells_low)),
    error = function(e) { message("DESeq2 failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(pseudobulk_res)) {
    write_csv(pseudobulk_res, file.path(out_dir, "2_pseudobulk_deseq2_biotin_high_vs_low.csv"))
    sig_pb <- pseudobulk_res %>% filter(FDR < FDR_CUTOFF, abs(logFC) >= LOGFC_CUTOFF)
    message(nrow(sig_pb), " genes at FDR < ", FDR_CUTOFF, ", |logFC| >= ", LOGFC_CUTOFF)

    # Volcano - FDR of exactly 0 is floored to the smallest representable
    # double before -log10() so it doesn't plot as Inf and blow out the
    # y-axis scale for every other gene.
    volcano_df <- pseudobulk_res %>%
      mutate(FDR_floor = pmax(FDR, .Machine$double.xmin),
             sig = FDR < FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)
    p_volcano_pb <- ggplot(volcano_df, aes(x = logFC, y = -log10(FDR_floor), color = sig)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c(`TRUE` = "firebrick3", `FALSE` = "grey70"), guide = "none") +
      geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
      ggrepel::geom_text_repel(data = ~filter(.x, sig), aes(label = gene), size = 3, max.overlaps = 20) +
      labs(title = paste0(CELLTYPE_VAL, " - pseudobulk DESeq2 volcano, biotin high vs low"),
           x = "log2FC (high vs low)", y = expression(-log[10](FDR)))
    ggsave(file.path(out_dir, "2_pseudobulk_deseq2_volcano_biotin_high_vs_low.pdf"), p_volcano_pb,
           width = 6, height = 5)
  }
}

# ============================================================
# 3. Differential CITE-seq binding (DSB) - per-protein Wilcoxon,
#    biotin_high vs biotin_low. Rank-based test on the raw (unfloored) DSB
#    values - negative-vs-negative ordering is real signal, not noise.
#    Reported means are floored at 0 for display only, same convention as
#    22.patient_qc_report.Rmd's run_cite_wilcoxon_de(). biotin itself is
#    included as a sanity check - it defines the grouping, so it will
#    trivially be the top hit.
# ============================================================

message("\n=== 3. CITE-seq DSB DE (per-protein Wilcoxon) ===")

if (length(cells_high) < MIN_N_DE || length(cells_low) < MIN_N_DE) {
  message("Fewer than ", MIN_N_DE, " cells in one biotin group - skipping CITE DSB DE.")
  cite_de_res <- NULL
} else {
  dsb_mat <- Seurat::GetAssayData(hsc, assay = "CITE_DSB", layer = "data")
  cite_de_res <- purrr::map_dfr(rownames(dsb_mat), function(p) {
    x <- dsb_mat[p, cells_high]
    y <- dsb_mat[p, cells_low]
    wt <- suppressWarnings(wilcox.test(x, y))
    mean_high <- mean(pmax(x, 0)); mean_low <- mean(pmax(y, 0))
    tibble(protein = p, n_high = length(x), n_low = length(y),
           mean_high = mean_high, mean_low = mean_low,
           mean_diff = mean_high - mean_low, pvalue = wt$p.value)
  }) %>% mutate(FDR = p.adjust(pvalue, method = "BH")) %>% arrange(FDR)

  write_csv(cite_de_res, file.path(out_dir, "3_cite_dsb_wilcoxon_biotin_high_vs_low.csv"))
  message(sum(cite_de_res$FDR < FDR_CUTOFF), " / ", nrow(cite_de_res),
          " proteins at FDR < ", FDR_CUTOFF)
}

# ============================================================
# 3b. Glycosylation (biotin) by fine cell type - raw biotin only (no
#     surface-load adjustment here; that's introduced in 3c below for the
#     protein-correlation regressions, not needed for this plot), with
#     cell counts per type. This is the sanity check for the "lumping fine
#     substates manufactures spurious signal" concern (methodology step 3,
#     not otherwise implemented): does biotin actually look similar across
#     HSC MPP's 3 fine substates (HSC/MPP-MkEry/MPP-MyLy, see the subset
#     comment above)? Run first, before 3c's pooled regression, since it's
#     exactly the assumption that regression depends on - if these look
#     similar across fine types, pooling is fine; if not, the
#     biotin_high/low split used everywhere in this script may really just
#     be a cell-type-composition split in disguise.
# ============================================================

message("\n=== 3b. Glycosylation (biotin) by fine cell type ===")

Seurat::DefaultAssay(hsc) <- "CITE_DSB"
dsb_all <- Seurat::GetAssayData(hsc, assay = "CITE_DSB", layer = "data")

if ("predicted_CellType" %in% colnames(hsc@meta.data)) {
  finetype_counts <- hsc@meta.data %>% as_tibble() %>% count(predicted_CellType, name = "n_cells")
  write_csv(finetype_counts, file.path(out_dir, "3b_fine_celltype_counts.csv"))
  message("Fine cell type counts:")
  print(finetype_counts)

  finetype_df <- tibble(fine_type = hsc$predicted_CellType, biotin = dsb_all[BIOTIN_MARKER, ]) %>%
    left_join(finetype_counts, by = c("fine_type" = "predicted_CellType")) %>%
    mutate(fine_type_label = paste0(fine_type, " (n=", n_cells, ")"))

  p_glyco_finetype <- ggplot(finetype_df, aes(x = fine_type_label, y = biotin)) +
    geom_violin(fill = "grey85") +
    geom_boxplot(width = 0.12, outlier.size = 0.5) +
    labs(title = paste0(CELLTYPE_VAL, " (", SAMPLE_NAME, ") - biotin (glycosylation) by fine cell type"),
         x = "", y = "Biotin (raw DSB)")
  ggsave(file.path(out_dir, "3b_glyco_by_fine_celltype_violin.pdf"), p_glyco_finetype, width = 7, height = 5)

  kw_raw <- kruskal.test(dsb_all[BIOTIN_MARKER, ] ~ factor(hsc$predicted_CellType))
  message("Kruskal-Wallis, raw biotin across fine cell types: p = ", signif(kw_raw$p.value, 3))
} else {
  message("predicted_CellType (fine) not available - skipping fine-cell-type glyco check.")
}

# ============================================================
# 3c. CITE-seq residuals vs glycosylation (biotin). biotin is the only
#     glyco-relevant channel in this panel - it's the readout, not a
#     comparator, so it's excluded from the "other proteins" list below.
#
#     Step 1 (surface-load regression): bigger/more granular cells bind
#     more antibody and more biotin-reagent non-specifically, independent
#     of true biology. CD45 + HLA-ABC (MHC-I, the closest available proxy
#     to CD45/beta-2-microglobulin - no B2M antibody in this panel) are the
#     two most uniformly-expressed pan-leukocyte ADTs here, so their mean
#     is used as a per-cell "surface load" proxy; biotin and every other
#     protein are each regressed against it, and everything below works
#     with residuals, not raw DSB values.
#
#     Step 2 (protein vs glyco): for each other protein, fit
#     glyco_residual ~ protein_residual. A significant positive/negative
#     coefficient means cells with more of that protein carry more/less
#     glycosylation, above what's explained by cell size.
#
#     Scope, and what this is NOT: single sample (SAMPLE_NAME), so this is
#     a cell-level model within one patient - not a patient-level/mixed-
#     effects model (there is no second patient level here to put a
#     "(1 | patient)" term on), and not a substitute for a real AML-vs-
#     normal comparison (which needs multiple patients - see the header
#     note at the top of this script). Not applicable here: a
#     MOFA+/joint-factor step, since there is exactly one glyco channel
#     (biotin), not a multi-lectin panel.
# ============================================================

message("\n=== 3c. CITE-seq residuals vs glycosylation (biotin) ===")

SURFACE_LOAD_MARKERS <- c("CD45", "HLA-ABC")
# Always shown in the summary/scatter below regardless of rank - CD47
# ("don't eat me") and CD33 (canonical AML/myeloid marker) requested
# specifically, not because they ranked top by effect size or FDR.
GLYCO_ALWAYS_INCLUDE <- c("CD47", "CD33")

stopifnot(all(SURFACE_LOAD_MARKERS %in% rownames(hsc[["CITE_DSB"]])))
surface_load <- Matrix::colMeans(dsb_all[SURFACE_LOAD_MARKERS, , drop = FALSE])

residualize <- function(x, covariate) resid(lm(x ~ covariate))

glyco_resid <- residualize(dsb_all[BIOTIN_MARKER, ], surface_load)

other_proteins <- setdiff(rownames(dsb_all), c(BIOTIN_MARKER, SURFACE_LOAD_MARKERS))

glyco_protein_cor <- purrr::map_dfr(other_proteins, function(p) {
  protein_resid <- residualize(dsb_all[p, ], surface_load)
  fit <- lm(glyco_resid ~ protein_resid)
  s <- summary(fit)
  tibble(protein = p,
         coef = coef(fit)[["protein_resid"]],
         se = s$coefficients["protein_resid", "Std. Error"],
         pvalue = s$coefficients["protein_resid", "Pr(>|t|)"],
         r_squared = s$r.squared)
}) %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  arrange(FDR)

write_csv(glyco_protein_cor, file.path(out_dir, "3c_glyco_residual_vs_protein_residual.csv"))
message(sum(glyco_protein_cor$FDR < FDR_CUTOFF), " / ", nrow(glyco_protein_cor),
        " proteins associated with the glycosylation residual at FDR < ", FDR_CUTOFF)

# Summary plot: per-protein coefficient +/- SE, ranked by effect size
# (|coefficient|), not FDR - FDR conflates effect size with how precisely
# it's estimated (a big but noisy coefficient can lose to a small but very
# precise one), which isn't what "which proteins have more glyco" is
# asking. Significance is still shown, via point color, just not used to
# pick which proteins make the plot. CD47/CD33 always included (see
# GLYCO_ALWAYS_INCLUDE above) even if they don't rank in the top 20.
glyco_plot_df <- glyco_protein_cor %>%
  mutate(sig = FDR < FDR_CUTOFF) %>%
  arrange(desc(abs(coef)))
glyco_plot_df <- bind_rows(
  glyco_plot_df %>% slice_head(n = min(20, nrow(glyco_plot_df))),
  glyco_plot_df %>% filter(protein %in% GLYCO_ALWAYS_INCLUDE)
) %>% distinct(protein, .keep_all = TRUE)

p_glyco_summary <- ggplot(glyco_plot_df, aes(x = coef, y = reorder(protein, coef), color = sig)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = coef - se, xmax = coef + se), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(`TRUE` = "firebrick3", `FALSE` = "grey50"), guide = "none") +
  labs(title = paste0(CELLTYPE_VAL, " (", SAMPLE_NAME, ") - proteins vs glycosylation residual"),
       x = "Coefficient (glyco_resid ~ protein_resid), +/- SE", y = "")
ggsave(file.path(out_dir, "3c_glyco_residual_vs_protein_summary.pdf"), p_glyco_summary, width = 7, height = 8)

# Scatter for the top hits (by effect size, same criterion as the summary
# plot above, plus CD47/CD33), so the summary coefficients have an actual
# picture behind them.
top_glyco_proteins <- glyco_protein_cor %>% arrange(desc(abs(coef))) %>% slice_head(n = 4) %>% pull(protein)
top_glyco_proteins <- union(top_glyco_proteins, intersect(GLYCO_ALWAYS_INCLUDE, glyco_protein_cor$protein))
scatter_df <- tibble(cell = colnames(hsc), glyco_resid = glyco_resid) %>%
  bind_cols(
    map_dfc(top_glyco_proteins, function(p) {
      tibble(!!p := residualize(dsb_all[p, ], surface_load))
    })
  ) %>%
  pivot_longer(-c(cell, glyco_resid), names_to = "protein", values_to = "protein_resid")

p_glyco_scatter <- ggplot(scatter_df, aes(x = protein_resid, y = glyco_resid)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick3") +
  facet_wrap(~ protein, scales = "free_x") +
  labs(title = paste0(CELLTYPE_VAL, " (", SAMPLE_NAME, ") - top proteins vs glycosylation residual"),
       x = "Protein residual (surface-load-adjusted)", y = "Glycosylation (biotin) residual")
ggsave(file.path(out_dir, "3c_glyco_residual_vs_protein_top_scatter.pdf"), p_glyco_scatter, width = 8, height = 7)

# Bonus: same logic against RNA glycosyltransferase expression, same
# single-patient/cell-level scope and caveats as above.
glyco_genes <- intersect(c("ST6GAL1", "FUT8", "MGAT5", "B4GALT1"), rownames(hsc[["RNA"]]))
if (length(glyco_genes) > 0) {
  Seurat::DefaultAssay(hsc) <- "RNA"
  rna_expr <- Seurat::FetchData(hsc, vars = glyco_genes, layer = "data")
  glyco_gene_cor <- purrr::map_dfr(glyco_genes, function(g) {
    fit <- lm(glyco_resid ~ rna_expr[[g]])
    s <- summary(fit)
    tibble(gene = g, coef = coef(fit)[[2]], pvalue = s$coefficients[2, "Pr(>|t|)"])
  }) %>% mutate(FDR = p.adjust(pvalue, method = "BH")) %>% arrange(FDR)
  write_csv(glyco_gene_cor, file.path(out_dir, "3c_glyco_residual_vs_glycosyltransferase_RNA.csv"))
  message("Glycosyltransferase RNA vs glyco residual (ST6GAL1/FUT8/MGAT5/B4GALT1):")
  print(glyco_gene_cor)
} else {
  message("No glycosyltransferase genes (ST6GAL1/FUT8/MGAT5/B4GALT1) found in RNA assay.")
}

# ============================================================
# 4. Cell-cycle phase proportions, biotin_high vs biotin_low
# ============================================================

message("\n=== 4. Cell-cycle phase proportions ===")

phase_tab <- hsc@meta.data %>%
  as_tibble() %>%
  count(biotin_group, Phase) %>%
  group_by(biotin_group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

write_csv(phase_tab, file.path(out_dir, "4_cell_cycle_phase_proportions_biotin_high_vs_low.csv"))
print(phase_tab)

p_phase <- ggplot(phase_tab, aes(x = biotin_group, y = pct, fill = Phase)) +
  geom_col(position = "stack") +
  labs(title = paste0(CELLTYPE_VAL, " - cell-cycle phase proportions, biotin high vs low"),
       x = "", y = "% of cells", fill = "Phase")
ggsave(file.path(out_dir, "4_cell_cycle_phase_proportions_biotin_high_vs_low.pdf"), p_phase, width = 6, height = 5)
print(p_phase)

phase_chisq <- tryCatch(
  suppressWarnings(chisq.test(table(hsc$biotin_group, hsc$Phase))),
  error = function(e) NULL
)
if (!is.null(phase_chisq)) {
  message("Chi-square test (phase distribution, biotin high vs low): p = ", signif(phase_chisq$p.value, 3))
}

# ============================================================
# 5. UMAP visualizations - biotin split (categorical + continuous) and
#    cell-cycle phase only. Top genes/proteins from the DE analyses above
#    are violin plots instead (section 6 below) - a spatial embedding adds
#    nothing for "is this gene/protein different between two groups I
#    already defined", it's just noise around the same signal a violin
#    shows directly. Shown on three reductions: RNA-only ("umap") and
#    RNA+Protein/WNN ("umap_wnn") - both "normal" (de novo) embeddings,
#    same as 22.patient_qc_report.Rmd's "RNA and RNA+CITE (WNN) UMAPs"
#    section - plus the BoneMarrowMap reference projection
#    ("umap_projected"). Restricted to the HSC MPP cells (both biotin
#    groups) at their real coordinates within each embedding.
# ============================================================

message("\n=== 5. UMAP visualizations ===")

umap_dir <- file.path(out_dir, "umaps")
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)

reductions_to_plot <- c(
  "RNA"          = "umap",
  "WNN"          = "umap_wnn",
  "BoneMarrowMap" = "umap_projected"
)
reductions_to_plot <- reductions_to_plot[reductions_to_plot %in% Seurat::Reductions(hsc)]
message("Reductions available on this object: ", paste(reductions_to_plot, collapse = ", "))

for (i in seq_along(reductions_to_plot)) {
  red_name  <- reductions_to_plot[i]
  red_label <- names(reductions_to_plot)[i]

  # -- biotin split itself: categorical group + the raw continuous value,
  #    both cell here, no split.by needed - the point is to see whether the
  #    two colors/gradient separate spatially within one panel --
  Seurat::DefaultAssay(hsc) <- "CITE_DSB"

  p_group <- Seurat::DimPlot(hsc, reduction = red_name, group.by = "biotin_group",
                             cols = c(biotin_high = "firebrick3", biotin_low = "navy"),
                             pt.size = 0.4, shuffle = TRUE, seed = 1) +
    ggtitle(paste0(CELLTYPE_VAL, " - biotin high/low (", red_label, ")"))
  ggsave(file.path(umap_dir, paste0("5a_biotin_group_", red_name, ".pdf")), p_group, width = 6, height = 5)

  p_value <- Seurat::FeaturePlot(hsc, features = "biotin_value", reduction = red_name,
                                 pt.size = 0.4, min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(paste0(CELLTYPE_VAL, " - biotin value (", red_label, ")"))
  ggsave(file.path(umap_dir, paste0("5b_biotin_value_", red_name, ".pdf")), p_value, width = 6, height = 5)

  # -- cell-cycle phase, split by biotin group --
  p_phase_umap <- Seurat::DimPlot(hsc, reduction = red_name, group.by = "Phase",
                                  split.by = "biotin_group", pt.size = 0.4,
                                  shuffle = TRUE, seed = 1) +
    plot_annotation(title = paste0(CELLTYPE_VAL, " - cell-cycle phase, biotin high vs low (", red_label, ")"))
  ggsave(file.path(umap_dir, paste0("5c_cell_cycle_phase_", red_name, ".pdf")), p_phase_umap, width = 10, height = 5)
}

# ============================================================
# 6. Top genes/antibodies - violin, not UMAP. Comparing an already-defined
#    two-group split doesn't need a spatial embedding: a violin of
#    biotin_high vs biotin_low directly shows the distribution shift, one
#    panel per gene/antibody, computed once rather than once per reduction
#    like section 5. biotin itself excluded from the CITE panel - it
#    trivially separates the groups by construction, not informative.
# ============================================================

message("\n=== 6. Top genes/antibodies (violin) ===")

if (!is.null(findmarkers_res) && nrow(findmarkers_res) > 0) {
  top_fm_genes <- findmarkers_res %>% slice_min(p_val_adj, n = 6, with_ties = FALSE) %>% pull(gene)
  Seurat::DefaultAssay(hsc) <- "RNA"
  p_fm_vln <- Seurat::VlnPlot(hsc, features = top_fm_genes, assay = "RNA", layer = "data",
                              group.by = "biotin_group", pt.size = 0) &
    theme(legend.position = "none")
  ggsave(file.path(umap_dir, "6a_findmarkers_top_genes_violin.pdf"), p_fm_vln,
         width = 3 * min(length(top_fm_genes), 3), height = 4 * ceiling(length(top_fm_genes) / 3))
} else {
  message("No FindMarkers result/genes to plot.")
}

if (!is.null(pseudobulk_res) && nrow(pseudobulk_res) > 0) {
  top_pb_genes <- pseudobulk_res %>% slice_min(FDR, n = 6, with_ties = FALSE) %>% pull(gene)
  Seurat::DefaultAssay(hsc) <- "RNA"
  p_pb_vln <- Seurat::VlnPlot(hsc, features = top_pb_genes, assay = "RNA", layer = "data",
                              group.by = "biotin_group", pt.size = 0) &
    theme(legend.position = "none")
  ggsave(file.path(umap_dir, "6b_pseudobulk_top_genes_violin.pdf"), p_pb_vln,
         width = 3 * min(length(top_pb_genes), 3), height = 4 * ceiling(length(top_pb_genes) / 3))
} else {
  message("No pseudobulk DE result/genes to plot.")
}

if (!is.null(cite_de_res) && nrow(cite_de_res) > 0) {
  top_cite <- cite_de_res %>% filter(protein != BIOTIN_MARKER) %>%
    slice_min(FDR, n = 6, with_ties = FALSE) %>% pull(protein)
  if (length(top_cite) > 0) {
    Seurat::DefaultAssay(hsc) <- "CITE_DSB"
    p_cite_vln <- Seurat::VlnPlot(hsc, features = top_cite, assay = "CITE_DSB", layer = "data",
                                  group.by = "biotin_group", pt.size = 0) &
      theme(legend.position = "none")
    ggsave(file.path(umap_dir, "6c_cite_dsb_top_proteins_violin.pdf"), p_cite_vln,
           width = 3 * min(length(top_cite), 3), height = 4 * ceiling(length(top_cite) / 3))
  }
} else {
  message("No CITE DSB DE result/proteins to plot.")
}

message("\nDone. Output directory: ", out_dir, " (figures under ", umap_dir, ")")
