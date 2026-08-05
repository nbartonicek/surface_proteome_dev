#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 16 of 19
#
# CopyKAT on all cells pooled. First CNV pass.
#
# Frozen for the lab archive 2026-07-31 from scripts/9.copykat.R (mtime 2026-05-11).
# md5 of the original: bf86300348b3d0fe4d1f7df4a2606eb8
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(data.table)
  library(copykat)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260423_VH01624_453_222HWMYNX"

annotation_dir <- file.path("../results/seurat_annotated", run)

seurat_projected_rds <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

if (!file.exists(seurat_projected_rds)) {
  stop("Cannot find projected Seurat object: ", seurat_projected_rds)
}

seu <- readRDS(seurat_projected_rds)

copykat_out <- file.path(annotation_dir, "copykat_annotated")
dir.create(copykat_out, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load analysis bundle
# ----------------------------

DefaultAssay(seu) <- "RNA"

stopifnot("RNA" %in% Assays(seu))
stopifnot("sampleID" %in% colnames(seu@meta.data))

cat("Loaded cells:", ncol(seu), "\n")
cat("RNA genes:", nrow(seu[["RNA"]]), "\n")

# ----------------------------
# Optional filtering
# ----------------------------

if ("mapping_error_QC" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = mapping_error_QC == "Pass")
}

#if ("HTO_classification.global" %in% colnames(seu@meta.data)) {
#  seu <- subset(seu, subset = HTO_classification.global == "Singlet")
#}

cat("Cells after filtering:", ncol(seu), "\n")

# ----------------------------
# Prepare raw count matrix
# ----------------------------
# CopyKAT expects genes x cells raw counts.
# Your object appears to use Ensembl IDs in RNA rownames from the CellRanger workflow,
# so gene IDs may need conversion to symbols if CopyKAT fails. The upstream object
# contains RNA/HTO/ADT assays from your demux pipeline. :contentReference[oaicite:0]{index=0}

raw_counts <- GetAssayData(
  seu,
  assay = "RNA",
  layer = "counts"
)

# Remove genes detected in very few cells to reduce runtime
min_cells <- 5

keep_genes <- Matrix::rowSums(raw_counts > 0) >= min_cells
raw_counts <- raw_counts[keep_genes, ]

cat("Genes retained for CopyKAT:", nrow(raw_counts), "\n")
cat("Cells retained for CopyKAT:", ncol(raw_counts), "\n")

# CopyKAT sometimes dislikes sparse matrices
raw_counts_dense <- as.matrix(raw_counts)

# ----------------------------
# Run CopyKAT
# ----------------------------
# id.type:
#   "S" = gene symbols
#   "E" = Ensembl IDs
#
# Since your RNA rownames appear Ensembl-like, use id.type = "E".
# If this fails, convert Ensembl IDs to gene symbols and rerun with id.type = "S".

setwd(copykat_out)
normal_celltypes <- c(
  "Naive T",
  "CD4 Memory T",
  "CD8 Memory T",
  "NK",
  "B"
)
norm_cells <- rownames(seu@meta.data)[
  seu$predicted_CellType_Broad %in% normal_celltypes
]

copykat_res <- copykat(
  rawmat = raw_counts_dense,
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = "LK1_copykat",
  distance = "euclidean",
  norm.cell.names = norm_cells,
  output.seg = "FALSE",
  plot.genes = TRUE,
  genome = "hg20",
  n.cores = 4
)

saveRDS(
  copykat_res,
  file.path(copykat_out, "copykat_result.rds")
)

# ----------------------------
# Extract CopyKAT prediction
# ----------------------------

copykat_res <- readRDS(paste0(copykat_out,"/LK1_copykat_copykat_clustering_results.rds"))

copykat_pred <- read.table(paste0(copykat_out,"/LK1_copykat_copykat_prediction.txt"), header=T)

copykat_pred <- copykat_pred %>%
  as_tibble() %>%
  rename(
    cell = cell.names,
    copykat_call = copykat.pred
  )

write.csv(
  copykat_pred,
  file.path(copykat_out, "copykat_prediction_raw.csv"),
  row.names = FALSE
)

# Usually columns include:
# cell.names
# copykat.pred
# sometimes prediction confidence columns depending on version

# ----------------------------
# Add CopyKAT calls back to Seurat
# ----------------------------

seu$copykat_call <- NA_character_

common_cells <- intersect(colnames(seu), copykat_pred$cell)

seu$copykat_call[common_cells] <- copykat_pred$copykat_call[
  match(common_cells, copykat_pred$cell)
]

# Optional cleaner labels
seu$copykat_malignancy <- case_when(
  seu$copykat_call == "aneuploid" ~ "CNV_aberrant",
  seu$copykat_call == "diploid" ~ "CNV_neutral",
  is.na(seu$copykat_call) ~ "Not_called",
  TRUE ~ seu$copykat_call
)

# ----------------------------
# Add metadata to CopyKAT table
# ----------------------------

copykat_annotated <- copykat_pred %>%
  mutate(
    sampleID = seu$sampleID[cell],
    predicted_CellType = if ("predicted_CellType" %in% colnames(seu@meta.data)) {
      seu$predicted_CellType[cell]
    } else {
      NA_character_
    },
    predicted_CellType_Broad = if ("predicted_CellType_Broad" %in% colnames(seu@meta.data)) {
      seu$predicted_CellType_Broad[cell]
    } else {
      NA_character_
    },
    WNN_cluster = if ("WNN_cluster" %in% colnames(seu@meta.data)) {
      as.character(seu$WNN_cluster[cell])
    } else {
      NA_character_
    }
  )

write.csv(
  copykat_annotated,
  file.path(copykat_out, "copykat_prediction_with_metadata.csv"),
  row.names = FALSE
)

# ----------------------------
# Summary tables
# ----------------------------

sample_summary <- copykat_annotated %>%
  count(sampleID, copykat_call, name = "n_cells") %>%
  group_by(sampleID) %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  sample_summary,
  file.path(copykat_out, "copykat_call_composition_by_sampleID.csv"),
  row.names = FALSE
)

celltype_summary <- copykat_annotated %>%
  count(sampleID, predicted_CellType_Broad, copykat_call, name = "n_cells") %>%
  group_by(sampleID, predicted_CellType_Broad) %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  celltype_summary,
  file.path(copykat_out, "copykat_call_by_sampleID_and_celltype.csv"),
  row.names = FALSE
)

# ----------------------------
# Plots
# ----------------------------

cols <- brewer.pal(n = 8, name = "Set1")
copykat_levels <- levels(factor(seu$copykat_call))

copykat_cols <- setNames(
  cols[seq_along(copykat_levels)],
  copykat_levels
)

seu$copykat_call <- factor(seu$copykat_call, levels = copykat_levels)
sample_summary <- sample_summary %>%
  mutate(copykat_call = factor(copykat_call, levels = copykat_levels))

celltype_summary <- celltype_summary %>%
  mutate(copykat_call = factor(copykat_call, levels = copykat_levels))

pdf(file.path(copykat_out, "01_copykat_call_composition_by_sampleID.pdf"), width = 8, height = 5)
print(
  ggplot(sample_summary, aes(sampleID, percent, fill = copykat_call)) +
    geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
    theme_bw() +
    labs(
      x = "Sample",
      y = "Cells (%)",
      fill = "CopyKAT call",
      title = "CopyKAT malignant/normal-like composition"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()
if ("umap_projected" %in% Reductions(seu)) {
  
  pdf(file.path(copykat_out, "02_copykat_calls_on_projected_RNA_UMAP.pdf"), width = 8, height = 6)
  
  print(
    DimPlot(
      seu,
      reduction = "umap_projected",
      group.by = "copykat_call",
      cols = copykat_cols
    ) +
      ggtitle("CopyKAT calls on projected RNA UMAP")
  )
  
  dev.off()
}

if ("umap_projected" %in% Reductions(seu)) {
  
  pdf(file.path(copykat_out, "02a_copykat_calls_on_projected_RNA_UMAP.pdf"), width = 8, height = 6)
  
  print(
    DimPlot(
      seu,
      reduction = "umap_projected",
      group.by = "copykat_call",
      split.by = "sampleID",
      cols = copykat_cols,
      ncol = 2
    ) +
      ggtitle("CopyKAT calls on projected RNA UMAP by sampleID")
  )
  
  dev.off()
}

if ("predicted_CellType_Broad" %in% colnames(seu@meta.data)) {
  
  pdf(file.path(copykat_out, "03_copykat_by_BoneMarrowMap_broad_celltype.pdf"), width = 10, height = 5)
  
  print(
    ggplot(
      celltype_summary,
      aes(predicted_CellType_Broad, percent, fill = copykat_call)
    ) +
      geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
      facet_wrap(~ sampleID) +
      scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
      theme_bw() +
      labs(
        x = "BoneMarrowMap broad cell type",
        y = "Cells (%)",
        fill = "CopyKAT call"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  )
  
  dev.off()
  pdf(file.path(copykat_out, "03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf"), width = 10, height = 5)
  
  print(
    seu@meta.data %>%
      as_tibble() %>%
      filter(!is.na(predicted_CellType_Broad), !is.na(copykat_call)) %>%
      count(sampleID, predicted_CellType_Broad, copykat_call) %>%
      ggplot(
        aes(predicted_CellType_Broad, n, fill = copykat_call)
      ) +
      geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
      facet_wrap(~ sampleID) +
      scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
      theme_bw() +
      labs(
        x = "BoneMarrowMap broad cell type",
        y = "Cell count",
        fill = "CopyKAT call"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  )
  
  dev.off()
}



# ----------------------------
# Save Seurat object with CopyKAT calls
# ----------------------------

saveRDS(
  seu,
  file.path(copykat_out, "query_with_copykat_calls.rds")
)

# Also update a compact bundle
copykat_bundle <- list(
  seu_copykat = seu,
  copykat_result = copykat_res,
  copykat_prediction = copykat_annotated,
  sample_summary = sample_summary,
  celltype_summary = celltype_summary
)

saveRDS(
  copykat_bundle,
  file.path(copykat_out, "copykat_bundle.rds")
)

cat("\nDone.\n")
cat("CopyKAT output:", copykat_out, "\n")
cat("Saved Seurat object:", file.path(copykat_out, "query_with_copykat_calls.rds"), "\n")
cat("Saved bundle:", file.path(copykat_out, "copykat_bundle.rds"), "\n")

