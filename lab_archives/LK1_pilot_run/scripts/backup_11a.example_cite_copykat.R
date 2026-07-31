#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(RColorBrewer)
  library(pheatmap)
  library(patchwork)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260423_VH01624_453_222HWMYNX"

out_dir <- file.path("../results/seurat_demux", run)
copykat_out <- file.path(out_dir, "copykat")

seu_file <- file.path(copykat_out, "query_with_copykat_calls.rds")

cite_out <- file.path(out_dir, "CITE_Late_GMP_HBDN206_MNpCT")
dir.create(cite_out, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper
# ----------------------------

clean_adt_names <- function(x) {
  gsub("-[ACGT]{15}$", "", x)
}

make_unique_clean_names <- function(x) {
  make.unique(clean_adt_names(x), sep = "_")
}

# ----------------------------
# Load Seurat object with CopyKAT calls
# ----------------------------

seu <- readRDS(seu_file)

stopifnot("ADT" %in% Assays(seu))
stopifnot("sampleID" %in% colnames(seu@meta.data))
stopifnot("predicted_CellType_Broad" %in% colnames(seu@meta.data))

cat("Loaded cells:", ncol(seu), "\n")
cat("ADT markers:", nrow(seu[["ADT"]]), "\n")

# ----------------------------
# Filter target sample + Late GMP
# ----------------------------
target_sample <- "HBDN206-MNpCT"
target_cluster <- "Late GMP"

seu$sampleID_clean <- gsub("−", "-", seu$sampleID)

target_cells <- colnames(seu)[
  seu$sampleID_clean == target_sample &
    seu$predicted_CellType_Broad == target_cluster
]

cat("Target cells:", length(target_cells), "\n")

if (length(target_cells) < 20) {
  stop("Too few target cells found. Check sampleID and predicted_CellType_Broad names.")
}

late_gmp <- subset(seu, cells = target_cells)

# ----------------------------
# Shorten ADT feature names globally in this object
# ----------------------------

old_adt_names <- rownames(late_gmp[["ADT"]])
new_adt_names <- make_unique_clean_names(old_adt_names)

adt_name_lookup <- tibble(
  original_adt_name = old_adt_names,
  clean_adt_name = new_adt_names
)

write.csv(
  adt_name_lookup,
  file.path(cite_out, "ADT_name_lookup_original_to_clean.csv"),
  row.names = FALSE
)

rownames(late_gmp[["ADT"]]) <- new_adt_names

# ----------------------------
# Normalise ADT
# ----------------------------

DefaultAssay(late_gmp) <- "ADT"

late_gmp <- NormalizeData(
  late_gmp,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2,
  verbose = FALSE
)

late_gmp <- ScaleData(
  late_gmp,
  assay = "ADT",
  features = rownames(late_gmp[["ADT"]]),
  verbose = FALSE
)

# ----------------------------
# Find top 12 abundant / variable CITE markers
# ----------------------------

adt_mat <- GetAssayData(
  late_gmp,
  assay = "ADT",
  values = "data"
)

adt_summary <- tibble(
  marker = rownames(adt_mat),
  mean_CLR = Matrix::rowMeans(adt_mat),
  pct_positive = Matrix::rowMeans(adt_mat > 0) * 120,
  variance = apply(as.matrix(adt_mat), 1, var)
) %>%
  arrange(desc(mean_CLR))

write.csv(
  adt_summary,
  file.path(cite_out, "LateGMP_HBDN206_MNpCT_ADT_marker_summary.csv"),
  row.names = FALSE
)

top12_abundant <- adt_summary %>%
  slice_max(mean_CLR, n = 12, with_ties = FALSE) %>%
  pull(marker)

top12_variable <- adt_summary %>%
  slice_max(variance, n = 12, with_ties = FALSE) %>%
  pull(marker)

write.csv(
  tibble(top12_abundant = top12_abundant),
  file.path(cite_out, "top12_abundant_ADT_markers.csv"),
  row.names = FALSE
)

write.csv(
  tibble(top12_variable = top12_variable),
  file.path(cite_out, "top12_variable_ADT_markers.csv"),
  row.names = FALSE
)

# ----------------------------
# ADT-only PCA / UMAP / clustering within Late GMP
# ----------------------------

late_gmp <- RunPCA(
  late_gmp,
  assay = "ADT",
  features = rownames(late_gmp[["ADT"]]),
  reduction.name = "apca_late_gmp",
  reduction.key = "APCA_LGMP_",
  npcs = 20,
  verbose = FALSE
)

late_gmp <- FindNeighbors(
  late_gmp,
  reduction = "apca_late_gmp",
  dims = 1:15,
  graph.name = "late_gmp_adt_snn",
  verbose = FALSE
)

late_gmp <- FindClusters(
  late_gmp,
  graph.name = "late_gmp_adt_snn",
  resolution = 0.3,
  algorithm = 1
)

cluster_col <- grep("^late_gmp_adt_snn_res", colnames(late_gmp@meta.data), value = TRUE)
cluster_col <- cluster_col[length(cluster_col)]

late_gmp$LateGMP_ADT_subcluster <- late_gmp@meta.data[[cluster_col]]

late_gmp <- RunUMAP(
  late_gmp,
  reduction = "apca_late_gmp",
  dims = 1:15,
  reduction.name = "late_gmp_adt_umap",
  reduction.key = "LGMPADTUMAP_",
  verbose = FALSE
)

# ----------------------------
# Colours
# ----------------------------

subclust_levels <- levels(factor(late_gmp$LateGMP_ADT_subcluster))

subclust_cols <- setNames(
  brewer.pal(
    n = max(3, min(8, length(subclust_levels))),
    name = "Set2"
  )[seq_along(subclust_levels)],
  subclust_levels
)

copykat_levels <- levels(factor(late_gmp$copykat_call))

copykat_cols <- setNames(
  brewer.pal(
    n = max(3, min(8, length(copykat_levels))),
    name = "Set1"
  )[seq_along(copykat_levels)],
  copykat_levels
)

# ----------------------------
# Plot 1: Late GMP ADT UMAP subclusters
# ----------------------------

pdf(file.path(cite_out, "01_LateGMP_ADT_subclusters.pdf"), width = 7, height = 6)

print(
  DimPlot(
    late_gmp,
    reduction = "late_gmp_adt_umap",
    group.by = "LateGMP_ADT_subcluster",
    label = TRUE,
    repel = TRUE,
    cols = subclust_cols
  ) +
    ggtitle("HBDN206-MNpCT Late GMP: ADT-defined subclusters")
)

dev.off()

# ----------------------------
# Plot 2: CopyKAT calls in Late GMP ADT space
# ----------------------------

if ("copykat_call" %in% colnames(late_gmp@meta.data)) {
  
  pdf(file.path(cite_out, "02_LateGMP_CopyKAT_on_ADT_UMAP.pdf"), width = 7, height = 6)
  
  print(
    DimPlot(
      late_gmp,
      reduction = "late_gmp_adt_umap",
      group.by = "copykat_call",
      cols = copykat_cols
    ) +
      ggtitle("CopyKAT calls within HBDN206-MNpCT Late GMP")
  )
  
  dev.off()
}

# ----------------------------
# Plot 3: top 12 abundant CITE markers
# ----------------------------
top12_abundant_clean <- clean_adt_names(top12_abundant)

p_abundant <- FeaturePlot(
  late_gmp,
  reduction = "late_gmp_adt_umap",
  features = top12_abundant,
  ncol = 4,
  combine = FALSE
)

# Replace titles
for (i in seq_along(p_abundant)) {
  p_abundant[[i]] <- p_abundant[[i]] +
    ggtitle(top12_abundant_clean[i])
}

pdf(file.path(cite_out, "03_top12_abundant_CITE_markers_LateGMP.pdf"),
    width = 14, height = 8)

print(
  wrap_plots(p_abundant, ncol = 4) +
    plot_annotation(
      title = "Top 12 abundant CITE-seq markers: HBDN206-MNpCT Late GMP"
    )
)

dev.off()
# ----------------------------
# Plot 4: top 12 variable CITE markers
# ----------------------------
# ----------------------------
# Plot 4: top 12 variable CITE markers
# ----------------------------

top12_variable_clean <- clean_adt_names(top12_variable)

p_variable <- FeaturePlot(
  late_gmp,
  reduction = "late_gmp_adt_umap",
  features = top12_variable,
  ncol = 4,
  combine = FALSE
)

for (i in seq_along(p_variable)) {
  p_variable[[i]] <- p_variable[[i]] +
    ggtitle(top12_variable_clean[i])
}

pdf(file.path(cite_out, "04_top12_variable_CITE_markers_LateGMP.pdf"),
    width = 14, height = 8)

print(
  wrap_plots(p_variable, ncol = 4) +
    plot_annotation(
      title = "Top 12 variable CITE-seq markers: HBDN206-MNpCT Late GMP"
    )
)

dev.off()
# ----------------------------
# Find ADT markers separating Late GMP subpopulations
# ----------------------------

DefaultAssay(late_gmp) <- "ADT"
Idents(late_gmp) <- "LateGMP_ADT_subcluster"

late_gmp_adt_markers <- FindAllMarkers(
  late_gmp,
  assay = "ADT",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  logfc.threshold = 0.05,
  min.pct = 0.05
)

late_gmp_adt_markers <- late_gmp_adt_markers %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  late_gmp_adt_markers,
  file.path(cite_out, "LateGMP_ADT_markers_by_subcluster.csv"),
  row.names = FALSE
)

top_sep_markers <- late_gmp_adt_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_sep_markers,
  file.path(cite_out, "top_ADT_markers_separating_LateGMP_subclusters.csv"),
  row.names = FALSE
)

top_sep_features <- unique(top_sep_markers$gene)

# ----------------------------
# Plot 5: most separating CITE markers on UMAP
# ----------------------------
top_sep_features_clean <- clean_adt_names(top_sep_features)

p_sep <- FeaturePlot(
  late_gmp,
  reduction = "late_gmp_adt_umap",
  features = top_sep_features,
  ncol = 5,
  combine = FALSE
)

for (i in seq_along(p_sep)) {
  p_sep[[i]] <- p_sep[[i]] +
    ggtitle(top_sep_features_clean[i])
}

pdf(file.path(cite_out, "05_top_separating_CITE_markers_on_ADT_UMAP.pdf"),
    width = 14, height = 12)

print(
  wrap_plots(p_sep, ncol = 5) +
    plot_annotation(
      title = "CITE markers separating ADT-defined Late GMP subpopulations"
    )
)

dev.off()
# ----------------------------
# Plot 6: heatmap of separating CITE markers
# ----------------------------

pdf(file.path(cite_out, "06_top_separating_CITE_markers_heatmap.pdf"), width = 12, height = 8)

print(
  DoHeatmap(
    late_gmp,
    assay = "ADT",
    features = top_sep_features,
    group.by = "LateGMP_ADT_subcluster"
  ) +
    NoLegend()
)

dev.off()

# ----------------------------
# Plot 7: dot plot of separating CITE markers
# ----------------------------

pdf(file.path(cite_out, "07_top_separating_CITE_markers_dotplot.pdf"), width = 12, height = 6)

print(
  DotPlot(
    late_gmp,
    assay = "ADT",
    features = top_sep_features,
    group.by = "LateGMP_ADT_subcluster"
  ) +
    RotatedAxis() +
    ggtitle("ADT markers separating Late GMP subclusters")
)

dev.off()

# ----------------------------
# Composition of ADT subclusters by CopyKAT
# ----------------------------

if ("copykat_call" %in% colnames(late_gmp@meta.data)) {
  
  copykat_subcluster_summary <- late_gmp@meta.data %>%
    count(LateGMP_ADT_subcluster, copykat_call, name = "n_cells") %>%
    group_by(LateGMP_ADT_subcluster) %>%
    mutate(percent = 120 * n_cells / sum(n_cells)) %>%
    ungroup()
  
  write.csv(
    copykat_subcluster_summary,
    file.path(cite_out, "LateGMP_ADT_subcluster_by_CopyKAT.csv"),
    row.names = FALSE
  )
  
  pdf(file.path(cite_out, "08_LateGMP_ADT_subcluster_CopyKAT_composition.pdf"), width = 7, height = 5)
  
  print(
    ggplot(
      copykat_subcluster_summary,
      aes(LateGMP_ADT_subcluster, percent, fill = copykat_call)
    ) +
      geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
      theme_bw() +
      labs(
        x = "Late GMP ADT subcluster",
        y = "Cells (%)",
        fill = "CopyKAT call",
        title = "CopyKAT composition of Late GMP ADT subclusters"
      )
  )
  
  dev.off()
}

# ----------------------------
# Save object
# ----------------------------

saveRDS(
  late_gmp,
  file.path(cite_out, "HBDN206_MNpCT_LateGMP_CITE_subclustered.rds")
)

cat("\nDone.\n")
cat("Output:", cite_out, "\n")

