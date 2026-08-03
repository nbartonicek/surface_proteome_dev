#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 09 of 19
#
# Fold the CITE assays into one object and write the analysis_bundle the later steps read.
#
# Frozen for the lab archive 2026-07-31 from scripts/backup/11.combine_cite.R (mtime 2026-05-03).
# md5 of the original: c3bd168321fa1e8cf6eec5e1c6c1a960
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260423_VH01624_453_222HWMYNX"

out_dir <- file.path("../results/seurat_demux", run)
bundle_file <- file.path(out_dir, "analysis_bundle.rds")

cite_out <- file.path(out_dir, "CITE_refinement")
dir.create(cite_out, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load analysis bundle
# ----------------------------

bundle <- readRDS(bundle_file)

query <- bundle$query
seu <- query

stopifnot("RNA" %in% Assays(seu))
stopifnot("ADT" %in% Assays(seu))

if (!"predicted_CellType_Broad" %in% colnames(seu@meta.data)) {
  stop("Missing predicted_CellType_Broad in query object.")
}

if (!"sampleID" %in% colnames(seu@meta.data)) {
  stop("Missing sampleID in query object.")
}

celltype_cols <- bundle$celltype_cols
lineage_cols <- bundle$lineage_cols

cat("Loaded bundle:\n")
cat("Cells:", ncol(seu), "\n")
cat("RNA genes:", nrow(seu[["RNA"]]), "\n")
cat("ADT features:", nrow(seu[["ADT"]]), "\n")

# ----------------------------
# Optional filtering
# ----------------------------

if ("mapping_error_QC" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = mapping_error_QC == "Pass")
}

if ("HTO_classification.global" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = HTO_classification.global == "Singlet")
}

cat("Cells after filtering:", ncol(seu), "\n")

# ----------------------------
# Helper functions
# ----------------------------

find_adt <- function(pattern, object = seu) {
  grep(pattern, rownames(object[["ADT"]]), value = TRUE, ignore.case = TRUE)
}

pick_adt <- function(patterns, object = seu) {
  hits <- unique(unlist(lapply(patterns, find_adt, object = object)))
  hits[hits %in% rownames(object[["ADT"]])]
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(cite_out, filename),
    plot = plot,
    width = width,
    height = height
  )
}

# ----------------------------
# ADT feature lookup
# ----------------------------

adt_feature_df <- tibble(
  adt_feature = rownames(seu[["ADT"]]),
  adt_marker = gsub("-[ACGT]+$", "", rownames(seu[["ADT"]]))
)

write.csv(
  adt_feature_df,
  file.path(cite_out, "ADT_feature_lookup.csv"),
  row.names = FALSE
)

# ----------------------------
# Ensure RNA PCA exists
# ----------------------------

DefaultAssay(seu) <- "RNA"

if (!"pca" %in% Reductions(seu)) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
}

# ----------------------------
# ADT-only analysis
# ----------------------------

DefaultAssay(seu) <- "ADT"

seu <- NormalizeData(
  seu,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2,
  verbose = FALSE
)

seu <- ScaleData(
  seu,
  assay = "ADT",
  features = rownames(seu[["ADT"]]),
  verbose = FALSE
)

seu <- RunPCA(
  seu,
  assay = "ADT",
  features = rownames(seu[["ADT"]]),
  reduction.name = "apca",
  reduction.key = "APCA_",
  npcs = 30,
  verbose = FALSE
)

seu <- FindNeighbors(
  seu,
  reduction = "apca",
  dims = 1:20,
  graph.name = "adt_snn",
  verbose = FALSE
)

seu <- FindClusters(
  seu,
  graph.name = "adt_snn",
  resolution = 0.4,
  algorithm = 1
)

adt_cluster_col <- grep("^adt_snn_res", colnames(seu@meta.data), value = TRUE)
adt_cluster_col <- adt_cluster_col[length(adt_cluster_col)]
seu$ADT_cluster <- seu@meta.data[[adt_cluster_col]]

seu <- RunUMAP(
  seu,
  reduction = "apca",
  dims = 1:20,
  reduction.name = "adt_umap",
  reduction.key = "ADTUMAP_",
  verbose = FALSE
)

# ----------------------------
# WNN RNA + ADT analysis
# ----------------------------

DefaultAssay(seu) <- "RNA"

seu <- FindMultiModalNeighbors(
  seu,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:20),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)

seu <- RunUMAP(
  seu,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = FALSE
)

seu <- FindClusters(
  seu,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.4
)

wnn_cluster_col <- grep("^wsnn_res", colnames(seu@meta.data), value = TRUE)
wnn_cluster_col <- wnn_cluster_col[length(wnn_cluster_col)]
seu$WNN_cluster <- seu@meta.data[[wnn_cluster_col]]

# ----------------------------
# UMAP plots
# ----------------------------

p1 <- DimPlot(
  seu,
  reduction = "adt_umap",
  group.by = "predicted_CellType_Broad",
  label = TRUE,
  repel = TRUE,
  cols = celltype_cols
) +
  ggtitle("ADT UMAP: BoneMarrowMap broad labels")

p2 <- DimPlot(
  seu,
  reduction = "adt_umap",
  group.by = "ADT_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("ADT-only clusters")

p3 <- DimPlot(
  seu,
  reduction = "adt_umap",
  group.by = "sampleID"
) +
  ggtitle("ADT UMAP: sample")

save_plot(p1, "01_ADT_UMAP_BoneMarrowMap_broad.pdf", 10, 8)
save_plot(p2, "02_ADT_UMAP_ADT_clusters.pdf", 10, 8)
save_plot(p3, "03_ADT_UMAP_sampleID.pdf", 10, 8)

p4 <- DimPlot(
  seu,
  reduction = "wnn.umap",
  group.by = "predicted_CellType_Broad",
  label = TRUE,
  repel = TRUE,
  cols = celltype_cols
) +
  ggtitle("WNN UMAP: BoneMarrowMap broad labels")

p5 <- DimPlot(
  seu,
  reduction = "wnn.umap",
  group.by = "WNN_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("WNN clusters")

p6 <- DimPlot(
  seu,
  reduction = "wnn.umap",
  group.by = "sampleID"
) +
  ggtitle("WNN UMAP: sample")

save_plot(p4, "04_WNN_UMAP_BoneMarrowMap_broad.pdf", 10, 8)
save_plot(p5, "05_WNN_UMAP_WNN_clusters.pdf", 10, 8)
save_plot(p6, "06_WNN_UMAP_sampleID.pdf", 10, 8)

# ----------------------------
# Canonical ADT marker feature plots
# ----------------------------

canonical_adt <- pick_adt(c(
  "CD34", "c-Kit", "CD117", "Thy1",
  "CD33", "CD13", "CD14", "CD16", "CD11b",
  "HLA", "CD80", "CD86",
  "CD3", "CD4", "CD8", "CD19", "CD20", "CD56",
  "CD38", "Syndecan", "BCMA"
))

write.csv(
  tibble(selected_adt_features = canonical_adt),
  file.path(cite_out, "selected_canonical_ADT_features.csv"),
  row.names = FALSE
)

if (length(canonical_adt) > 0) {
  p_adt <- FeaturePlot(
    seu,
    reduction = "wnn.umap",
    features = canonical_adt,
    ncol = 4
  )
  
  save_plot(
    p_adt,
    "07_WNN_canonical_ADT_markers.pdf",
    width = 16,
    height = ceiling(length(canonical_adt) / 4) * 3
  )
}

# ----------------------------
# Composition: WNN clusters per sample
# ----------------------------

wnn_comp <- seu@meta.data %>%
  count(sampleID, WNN_cluster, name = "n_cells") %>%
  group_by(sampleID) %>%
  mutate(
    total_cells = sum(n_cells),
    percent = 100 * n_cells / total_cells
  ) %>%
  ungroup()

write.csv(
  wnn_comp,
  file.path(cite_out, "WNN_cluster_composition_by_sampleID.csv"),
  row.names = FALSE
)

p_wnn_comp <- ggplot(
  wnn_comp,
  aes(sampleID, percent, fill = WNN_cluster)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Composition (%)",
    fill = "WNN cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(p_wnn_comp, "08_WNN_cluster_composition_by_sampleID.pdf", 9, 5)

# ----------------------------
# Cross-tab: BoneMarrowMap labels vs WNN clusters
# ----------------------------

label_vs_wnn <- seu@meta.data %>%
  count(predicted_CellType_Broad, WNN_cluster, name = "n_cells") %>%
  group_by(predicted_CellType_Broad) %>%
  mutate(percent_within_label = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  label_vs_wnn,
  file.path(cite_out, "BoneMarrowMapBroad_vs_WNN_clusters.csv"),
  row.names = FALSE
)

label_vs_wnn_mat <- label_vs_wnn %>%
  select(predicted_CellType_Broad, WNN_cluster, percent_within_label) %>%
  pivot_wider(
    names_from = WNN_cluster,
    values_from = percent_within_label,
    values_fill = 0
  )

mat_label_wnn <- as.matrix(label_vs_wnn_mat[, -1])
rownames(mat_label_wnn) <- label_vs_wnn_mat$predicted_CellType_Broad

pdf(file.path(cite_out, "09_BoneMarrowMapBroad_vs_WNN_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  mat_label_wnn,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "BoneMarrowMap broad labels vs WNN clusters"
)
dev.off()

# ----------------------------
# ADT marker discovery for WNN clusters
# ----------------------------

DefaultAssay(seu) <- "ADT"
Idents(seu) <- "WNN_cluster"

adt_markers <- FindAllMarkers(
  seu,
  assay = "ADT",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.05
)

write.csv(
  adt_markers,
  file.path(cite_out, "ADT_markers_by_WNN_cluster.csv"),
  row.names = FALSE
)

top_adt_markers <- adt_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_adt_markers,
  file.path(cite_out, "top5_ADT_markers_by_WNN_cluster.csv"),
  row.names = FALSE
)

top_features <- unique(top_adt_markers$gene)

if (length(top_features) > 1) {
  pdf(file.path(cite_out, "10_top_ADT_markers_by_WNN_cluster_heatmap.pdf"), width = 12, height = 10)
  print(
    DoHeatmap(
      seu,
      assay = "ADT",
      features = top_features,
      group.by = "WNN_cluster"
    ) +
      NoLegend()
  )
  dev.off()
}

# ----------------------------
# Within-lineage refinement examples
# ----------------------------

lineages_to_test <- c(
  "HSC MPP", "LMPP", "GMP", "Early GMP", "Late GMP",
  "Monocyte", "cDC", "pDC", "B", "Naive T", "NK"
)

lineage_refinement <- seu@meta.data %>%
  filter(predicted_CellType_Broad %in% lineages_to_test) %>%
  count(predicted_CellType_Broad, WNN_cluster, name = "n_cells") %>%
  group_by(predicted_CellType_Broad) %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  lineage_refinement,
  file.path(cite_out, "within_BoneMarrowMap_lineage_WNN_refinement.csv"),
  row.names = FALSE
)

p_lineage_refine <- ggplot(
  lineage_refinement,
  aes(predicted_CellType_Broad, percent, fill = WNN_cluster)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  theme_bw() +
  labs(
    x = "BoneMarrowMap broad label",
    y = "WNN cluster composition (%)",
    fill = "WNN cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(p_lineage_refine, "11_within_lineage_WNN_refinement.pdf", 11, 5)

# ----------------------------
# Save updated analysis bundle
# ----------------------------

cite_analysis_bundle <- list(
  seu_cite_refined = seu,
  adt_feature_lookup = adt_feature_df,
  canonical_adt = canonical_adt,
  wnn_comp = wnn_comp,
  label_vs_wnn = label_vs_wnn,
  adt_markers = adt_markers,
  top_adt_markers = top_adt_markers,
  lineage_refinement = lineage_refinement,
  celltype_cols = celltype_cols,
  lineage_cols = lineage_cols
)

saveRDS(
  cite_analysis_bundle,
  file.path(cite_out, "CITE_refinement_bundle.rds")
)

saveRDS(
  seu,
  file.path(cite_out, "query_CITE_WNN_refined.rds")
)

cat("\nDone.\n")
cat("CITE refinement output:", cite_out, "\n")
cat("Saved object:", file.path(cite_out, "query_CITE_WNN_refined.rds"), "\n")
cat("Saved bundle:", file.path(cite_out, "CITE_refinement_bundle.rds"), "\n")
