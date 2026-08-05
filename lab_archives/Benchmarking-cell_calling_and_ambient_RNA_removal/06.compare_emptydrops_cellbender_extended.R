#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Cell calling and ambient RNA removal - step 06 of 7
#
# Extends step 05 with the +30k droplet arm, a four-way Venn, and per-category QC: genes, UMIs, MT counts and percent MT by barcode category.
#
# Kept alongside step 05 rather than replacing it: the two produced different
# outputs on different days and neither reproduces the other's. Step 05 wrote
# the 25k/32k Venns on 2026-06-24; this wrote the 32k30k arm, the four-way
# Venn and the QC violins on 2026-06-29.
#
# Frozen for the lab archive 2026-08-04 from scripts/19a.compare_emptydrops_cellbender.R (mtime 2026-06-29).
# md5 of the original: a4aac3e818ee61d448e5cd2e6ee8735d
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggVennDiagram)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(DropletUtils)
  library(SingleCellExperiment)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
run <- "260528_VH01624_464_222K7VKNX"
sample <- "LK2-GEX"

out_dir <- file.path(proj, "results", "cellbender_comparison", run, sample)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------

emptydrops_barcode_file <- file.path(
  proj, "results", "emptydrops", run, sample, "emptydrops_barcodes.txt"
)

cellranger_matrix_dir <- file.path(
  proj, "results", "cellranger_withbam", run, sample, "outs", "raw_feature_bc_matrix"
)

cellbender_dir <- file.path(
  proj, "results", "cellbender", run, paste0(sample, "_e30000")
)

cellbender_barcode_file <- file.path(
  cellbender_dir,
  paste0(sample, "_cellbender_cell_barcodes.csv")
)

cellbender_h5 <- file.path(
  cellbender_dir,
  paste0(sample, "_cellbender_filtered.h5")
)

stopifnot(file.exists(emptydrops_barcode_file))
stopifnot(dir.exists(cellranger_matrix_dir))
stopifnot(file.exists(cellbender_barcode_file))
stopifnot(file.exists(cellbender_h5))

# ------------------------------------------------------------------
# 2. Load barcode lists
# ------------------------------------------------------------------

emptydrops_barcodes <- unique(readLines(emptydrops_barcode_file))
cellbender_barcodes <- unique(readLines(cellbender_barcode_file))

cat("Barcode counts from caller outputs:\n")
cat("  EmptyDrops:        ", length(emptydrops_barcodes), "\n")
cat("  CellBender 32k30k: ", length(cellbender_barcodes), "\n")

# ------------------------------------------------------------------
# 3. Load CellRanger raw matrix for EmptyDrops object
# ------------------------------------------------------------------

cat("\nLoading CellRanger raw matrix:\n", cellranger_matrix_dir, "\n")

cr_counts <- Read10X(cellranger_matrix_dir)

if (is.list(cr_counts)) {
  cr_counts <- cr_counts[["Gene Expression"]]
}

emptydrops_barcodes <- intersect(emptydrops_barcodes, colnames(cr_counts))

seu_emptydrops <- CreateSeuratObject(
  counts = cr_counts[, emptydrops_barcodes, drop = FALSE],
  project = paste0(sample, "_EmptyDrops")
)

# ------------------------------------------------------------------
# 4. Load CellBender corrected filtered matrix
# ------------------------------------------------------------------

cat("\nLoading CellBender corrected matrix:\n", cellbender_h5, "\n")

cb_sce <- DropletUtils::read10xCounts(
  cellbender_h5,
  type = "HDF5",
  col.names = TRUE
)

cb_counts <- SingleCellExperiment::counts(cb_sce)

if ("Symbol" %in% colnames(SingleCellExperiment::rowData(cb_sce))) {
  rownames(cb_counts) <- make.unique(as.character(SingleCellExperiment::rowData(cb_sce)$Symbol))
}

cellbender_barcodes <- intersect(cellbender_barcodes, colnames(cb_counts))

seu_cellbender <- CreateSeuratObject(
  counts = cb_counts[, cellbender_barcodes, drop = FALSE],
  project = paste0(sample, "_CellBender_32k30k")
)

# ------------------------------------------------------------------
# 5. Store Seurat objects in a list
# ------------------------------------------------------------------

seurat_objects <- list(
  EmptyDrops = seu_emptydrops,
  CellBender_32k30k = seu_cellbender
)

saveRDS(
  seurat_objects,
  file.path(out_dir, "seurat_objects_emptydrops_cellbender_32k30k.rds")
)

cat("\nCreated Seurat objects:\n")
cat("  EmptyDrops cells:        ", ncol(seurat_objects$EmptyDrops), "\n")
cat("  CellBender 32k30k cells: ", ncol(seurat_objects$CellBender_32k30k), "\n")

# ------------------------------------------------------------------
# 6. Define barcode overlap categories
# ------------------------------------------------------------------

emptydrops_cells <- colnames(seurat_objects$EmptyDrops)
cellbender_cells <- colnames(seurat_objects$CellBender_32k30k)

common_barcodes <- intersect(emptydrops_cells, cellbender_cells)
emptydrops_only <- setdiff(emptydrops_cells, cellbender_cells)
cellbender_only <- setdiff(cellbender_cells, emptydrops_cells)

cat("\nBarcode overlap:\n")
cat("  Common:          ", length(common_barcodes), "\n")
cat("  EmptyDrops only: ", length(emptydrops_only), "\n")
cat("  CellBender only: ", length(cellbender_only), "\n")

#Then your case_when() will work.
# ------------------------------------------------------------------
# 5. Extract QC stats
# ------------------------------------------------------------------
add_qc_metrics <- function(seu) {
  counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
  
  mt_genes <- grep("^MT-", rownames(counts_mat), value = TRUE)
  
  seu$nCount_MT <- if (length(mt_genes) > 0) {
    Matrix::colSums(counts_mat[mt_genes, , drop = FALSE])
  } else {
    rep(0, ncol(seu))
  }
  
  seu$percent.mt <- ifelse(
    seu$nCount_RNA > 0,
    100 * seu$nCount_MT / seu$nCount_RNA,
    0
  )
  
  seu
}

seurat_objects <- lapply(seurat_objects, add_qc_metrics)

get_qc_df <- function(seu, label, cells_use) {
  seu@meta.data %>%
    rownames_to_column("barcode") %>%
    filter(barcode %in% cells_use) %>%
    mutate(barcode_category = label) %>%
    select(
      barcode,
      barcode_category,
      nFeature_RNA,
      nCount_RNA,
      nCount_MT,
      percent.mt
    )
}

qc_df <- bind_rows(
  get_qc_df(
    seurat_objects$EmptyDrops,
    label = "Common",
    cells_use = common_barcodes
  ),
  get_qc_df(
    seurat_objects$CellBender_32k30k,
    label = "Common corrected",
    cells_use = common_barcodes
  ),
  get_qc_df(
    seurat_objects$EmptyDrops,
    label = "EmptyDrops only",
    cells_use = emptydrops_only
  ),
  get_qc_df(
    seurat_objects$CellBender_32k30k,
    label = "CellBender only",
    cells_use = cellbender_only
  )
)

qc_df$barcode_category <- factor(
  qc_df$barcode_category,
  levels = c(
    "Common",
    "Common corrected",
    "EmptyDrops only",
    "CellBender only"
  )
)



write.csv(
  qc_df,
  file.path(out_dir, "emptydrops_cellbender_32k30k_qc_by_barcode_category.csv"),
  row.names = FALSE
)

qc_summary <- qc_df %>%
  group_by(barcode_category) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_umis = median(nCount_RNA),
    median_mt_counts = median(nCount_MT),
    median_percent_mt = median(percent.mt),
    mean_percent_mt = mean(percent.mt),
    pct_mt_gt_20 = mean(percent.mt > 20) * 100,
    pct_mt_gt_50 = mean(percent.mt > 50) * 100,
    pct_mt_gt_80 = mean(percent.mt > 80) * 100,
    pct_mt_gt_90 = mean(percent.mt > 90) * 100,
    .groups = "drop"
  )

write.csv(
  qc_summary,
  file.path(out_dir, "emptydrops_cellbender_32k30k_qc_summary.csv"),
  row.names = FALSE
)

print(qc_summary)

# ------------------------------------------------------------------
# 6. Venn diagram
# ------------------------------------------------------------------

barcode_sets <- list(
  EmptyDrops = emptydrops_cells,
  `CellBender 32k30k` = cellbender_cells
)

p_venn <- ggVennDiagram(barcode_sets, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#4292C6") +
  ggtitle(paste0("Cell calling comparison - ", sample)) +
  coord_cartesian(clip = "off") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.margin = margin(t = 10, r = 35, b = 10, l = 35)
  )

ggsave(
  file.path(out_dir, "venn_emptydrops_vs_cellbender_32k30k.pdf"),
  p_venn,
  width = 8,
  height = 6,
  limitsize = FALSE
)

# ------------------------------------------------------------------
# 7. Violin plots
# ------------------------------------------------------------------
plot_violin <- function(df, y, ylab, filename, log_y = FALSE) {
  p <- ggplot(df, aes(x = barcode_category, y = .data[[y]], fill = barcode_category)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.65) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    ) +
    labs(
      title = paste0(sample, ": ", ylab),
      x = NULL,
      y = ylab
    )
  
  #if (log_y) {
  #  p <- p + scale_y_continuous(trans = "log1p")
  #}
  #
  ggsave(file.path(out_dir, filename), p, width = 8, height = 5)
  p
}
p_genes <- plot_violin(
  qc_df,
  "nFeature_RNA",
  "Number of detected genes",
  "violin_genes_by_barcode_category_unfaceted.pdf",
  log_y = TRUE
)

p_umis <- plot_violin(
  qc_df,
  "nCount_RNA",
  "Number of UMIs",
  "violin_umis_by_barcode_category_unfaceted.pdf",
  log_y = TRUE
)

p_mt_counts <- plot_violin(
  qc_df,
  "nCount_MT",
  "Mitochondrial UMI counts",
  "violin_mt_counts_by_barcode_category_unfaceted.pdf",
  log_y = TRUE
)

p_percent_mt <- plot_violin(
  qc_df,
  "percent.mt",
  "Percent mitochondrial",
  "violin_percent_mt_by_barcode_category_unfaceted.pdf",
  log_y = FALSE
)
p_combined <- p_genes / p_umis / p_mt_counts / p_percent_mt

ggsave(
  file.path(out_dir, "violin_qc_metrics_by_barcode_category_combined.pdf"),
  p_combined,
  width = 10,
  height = 16
)

cat("\nOutputs saved to:", out_dir, "\n")
cat("Done.\n")
