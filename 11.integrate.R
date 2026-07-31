library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)

run <- "260528_VH01624_464_222K7VKNX"
proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
sample_short <- "LK2"

annotation_dir <- file.path(
  proj,
  "results/seurat_annotated",
  run
)

seurat_file <- file.path(
  annotation_dir,
  "demux_singlets_annotated_seurat.rds"
)

out_dir <- file.path(annotation_dir, "harmony_integration_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(seurat_file)

# ----------------------------
# Find useful metadata columns
# ----------------------------

meta_cols <- colnames(seu@meta.data)
print(meta_cols)

annotation_col <- meta_cols[
  grepl("bone|atlas|celltype|cell_type|predicted|annotation|label", meta_cols, ignore.case = TRUE)
]

print(annotation_col)

# Pick one manually if needed
annotation_col <- annotation_col[1]

if (is.na(annotation_col)) {
  stop("Could not guess BoneMarrowAtlas annotation column. Check colnames(seu@meta.data).")
}

message("Using annotation column: ", annotation_col)

# Choose batch/sample column
batch_candidates <- c(
  "sample_name",
  "sample",
  "orig.ident",
  "Sample",
  "donor",
  "run",
  "library"
)

batch_col <- batch_candidates[batch_candidates %in% meta_cols][1]

if (is.na(batch_col)) {
  stop("Could not find batch column. Add one manually, e.g. batch_col <- 'sample_name'")
}

message("Using batch column: ", batch_col)

# ----------------------------
# Standard Seurat processing
# ----------------------------

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# ----------------------------
# Harmony
# ----------------------------

seu <- RunHarmony(
  object = seu,
  group.by.vars = batch_col,
  reduction = "pca",
  dims.use = 1:50,
  assay.use = "RNA",
  reduction.save = "harmony"
)

seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  reduction.key = "hUMAP_"
)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

# ----------------------------
# Plots
# ----------------------------

p_batch <- DimPlot(
  seu,
  reduction = "umap_harmony",
  group.by = batch_col
) +
  ggtitle(paste("Harmony UMAP by", batch_col))

p_annot <- DimPlot(
  seu,
  reduction = "umap_harmony",
  group.by = annotation_col,
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste("Harmony UMAP by", annotation_col))

p_cluster <- DimPlot(
  seu,
  reduction = "umap_harmony",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Harmony UMAP by Seurat clusters")

ggsave(
  file.path(out_dir, "harmony_umap_by_batch.pdf"),
  p_batch,
  width = 8,
  height = 6
)

ggsave(
  file.path(out_dir, "harmony_umap_by_bonemarrow_annotation.pdf"),
  p_annot,
  width = 10,
  height = 6
)

ggsave(
  file.path(out_dir, "harmony_umap_by_cluster.pdf"),
  p_cluster,
  width = 8,
  height = 6
)

ggsave(
  file.path(out_dir, "harmony_umap_combined.pdf"),
  p_batch + p_annot + p_cluster,
  width = 22,
  height = 7
)

# ----------------------------
# Integration diagnostics
# ----------------------------

tab_batch_cluster <- table(
  cluster = seu$seurat_clusters,
  batch = seu@meta.data[[batch_col]]
)

tab_annot_cluster <- table(
  cluster = seu$seurat_clusters,
  annotation = seu@meta.data[[annotation_col]]
)

write.csv(
  as.data.frame.matrix(tab_batch_cluster),
  file.path(out_dir, "cluster_by_batch_counts.csv")
)

write.csv(
  as.data.frame.matrix(tab_annot_cluster),
  file.path(out_dir, "cluster_by_bonemarrow_annotation_counts.csv")
)

# Proportions per cluster
prop_batch_cluster <- prop.table(tab_batch_cluster, margin = 1)
prop_annot_cluster <- prop.table(tab_annot_cluster, margin = 1)

write.csv(
  as.data.frame.matrix(prop_batch_cluster),
  file.path(out_dir, "cluster_by_batch_proportions.csv")
)

write.csv(
  as.data.frame.matrix(prop_annot_cluster),
  file.path(out_dir, "cluster_by_bonemarrow_annotation_proportions.csv")
)

# Save object
saveRDS(
  seu,
  file.path(out_dir, "demux_singlets_annotated_seurat_harmony.rds")
)

message("Done. Outputs written to: ", out_dir)