# ============================================================
# residual_cite_analysis.R
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(tibble)
library(purrr)
library(RColorBrewer)
library(caret)

# ============================================================
# Parameters
# ============================================================

sample_col <- "sample_name"
celltype_col <- "predicted_CellType_Broad"

out_dir <- "../results/residual_cite_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Colours
# ============================================================

# Assumes you already defined this previously
# If not, define manually here
#
# Example:
#
# celltype_cols <- c(
#   "HSC MPP" = "#1B9E77",
#   "LMPP" = "#66C2A5",
#   ...
# )

# ============================================================
# Basic setup
# ============================================================

DefaultAssay(seu_dsb) <- "RNA"

# Ensure RNA PCA exists
if (!"pca" %in% names(seu_dsb@reductions)) {
  
  seu_dsb <- NormalizeData(seu_dsb, verbose = FALSE)
  seu_dsb <- FindVariableFeatures(seu_dsb, verbose = FALSE)
  seu_dsb <- ScaleData(seu_dsb, verbose = FALSE)
  
  seu_dsb <- RunPCA(
    seu_dsb,
    npcs = 50,
    reduction.name = "pca",
    verbose = FALSE
  )
}

# ============================================================
# Extract metadata + matrices
# ============================================================

DefaultAssay(seu_dsb) <- "CITE_DSB"

adt <- GetAssayData(
  seu_dsb,
  assay = "CITE_DSB",
  layer = "data"
)

meta <- seu_dsb@meta.data %>%
  rownames_to_column("cell")

# ============================================================
# RNA PCs
# ============================================================

rna_pcs <- Embeddings(
  seu_dsb,
  reduction = "pca"
)[, 1:10, drop = FALSE]

rna_pcs <- as.data.frame(rna_pcs) %>%
  rownames_to_column("cell")

meta2 <- meta %>%
  left_join(rna_pcs, by = "cell")

rownames(meta2) <- meta2$cell

common_cells <- intersect(
  colnames(adt),
  rownames(meta2)
)

adt <- adt[, common_cells]
meta2 <- meta2[common_cells, ]

# ============================================================
# Residualise protein expression
# ============================================================

cat("Running residualisation...\n")

resid_mat <- matrix(
  NA_real_,
  nrow = nrow(adt),
  ncol = ncol(adt),
  dimnames = dimnames(adt)
)

for (marker in rownames(adt)) {
  
  df <- meta2
  df$cell <- rownames(df)
  df$protein <- as.numeric(adt[marker, rownames(df)])
  
  keep <- complete.cases(
    df[, c(
      "protein",
      "sample_name",
      "predicted_CellType_Broad",
      "PC_1", "PC_2", "PC_3", "PC_4", "PC_5"
    )]
  )
  
  df_fit <- df[keep, ]
  
  fit <- lm(
    protein ~
      sample_name +
      predicted_CellType_Broad +
      PC_1 + PC_2 + PC_3 + PC_4 + PC_5,
    data = df_fit
  )
  
  marker_resid <- rep(NA_real_, nrow(df))
  names(marker_resid) <- df$cell
  
  marker_resid[df_fit$cell] <- residuals(fit)
  
  resid_mat[marker, names(marker_resid)] <- marker_resid
}

# ============================================================
# Create residual assay
# ============================================================

seu_dsb[["CITE_residual"]] <- CreateAssayObject(
  data = resid_mat
)

DefaultAssay(seu_dsb) <- "CITE_residual"

# ============================================================
# Scale + PCA
# ============================================================

resid_features <- rownames(seu_dsb[["CITE_residual"]])

seu_dsb <- ScaleData(
  seu_dsb,
  assay = "CITE_residual",
  features = resid_features,
  verbose = FALSE
)

seu_dsb <- RunPCA(
  seu_dsb,
  assay = "CITE_residual",
  features = resid_features,
  reduction.name = "pca_cite_resid",
  reduction.key = "RESIDPC_",
  npcs = min(20, length(resid_features) - 1),
  verbose = FALSE
)

# ============================================================
# Clustering on residual signal
# ============================================================

seu_dsb <- FindNeighbors(
  seu_dsb,
  reduction = "pca_cite_resid",
  dims = 1:15,
  verbose = FALSE
)

seu_dsb <- FindClusters(
  seu_dsb,
  resolution = 0.4,
  cluster.name = "cite_resid_clusters",
  verbose = FALSE
)

seu_dsb <- RunUMAP(
  seu_dsb,
  reduction = "pca_cite_resid",
  dims = 1:15,
  reduction.name = "umap_cite_resid",
  reduction.key = "residUMAP_",
  verbose = FALSE
)

# ============================================================
# Cluster colours
# ============================================================

cluster_levels <- sort(unique(as.character(seu_dsb$cite_resid_clusters)))

cluster_cols <- setNames(
  colorRampPalette(
    brewer.pal(12, "Paired")
  )(length(cluster_levels)),
  cluster_levels
)

# ============================================================
# Plot residual clusters
# ============================================================

p1 <- DimPlot(
  seu_dsb,
  reduction = "umap_cite_resid",
  group.by = "cite_resid_clusters",
  cols = cluster_cols,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Residual CITE clusters")

ggsave(
  file.path(out_dir, "01_residual_cite_clusters.pdf"),
  p1,
  width = 8,
  height = 6
)

# ============================================================
# Plot residual clusters by sample
# ============================================================

p2 <- DimPlot(
  seu_dsb,
  reduction = "umap_cite_resid",
  group.by = "cite_resid_clusters",
  split.by = sample_col,
  cols = cluster_cols,
  label = TRUE,
  repel = TRUE,
  ncol = 3
) +
  ggtitle("Residual CITE clusters by sample")

ggsave(
  file.path(out_dir, "02_residual_cite_clusters_by_sample.pdf"),
  p2,
  width = 14,
  height = 8
)

# ============================================================
# Plot known lineages
# ============================================================

p3 <- DimPlot(
  seu_dsb,
  reduction = "umap_cite_resid",
  group.by = celltype_col,
  split.by = sample_col,
  cols = celltype_cols,
  ncol = 3
) +
  ggtitle("Known broad lineages")

ggsave(
  file.path(out_dir, "03_known_lineages_on_residual_umap.pdf"),
  p3,
  width = 15,
  height = 9
)

# ============================================================
# Find residual markers
# ============================================================

DefaultAssay(seu_dsb) <- "CITE_residual"

Idents(seu_dsb) <- "cite_resid_clusters"

resid_markers <- FindAllMarkers(
  seu_dsb,
  assay = "CITE_residual",
  layer = "data",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  test.use = "wilcox"
)

write.csv(
  resid_markers,
  file.path(out_dir, "04_residual_cite_markers.csv"),
  row.names = FALSE
)

# ============================================================
# Top markers per residual cluster
# ============================================================

top_markers <- resid_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 6, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_markers,
  file.path(out_dir, "05_top_residual_markers.csv"),
  row.names = FALSE
)

# ============================================================
# FeaturePlots of residual markers
# ============================================================

top_features <- unique(top_markers$gene)

p4 <- FeaturePlot(
  seu_dsb,
  reduction = "umap_cite_resid",
  features = top_features,
  order = TRUE,
  ncol = 3,
  keep.scale = "feature"
)

ggsave(
  file.path(out_dir, "06_top_residual_marker_featureplots.pdf"),
  p4,
  width = 14,
  height = 12
)

# ============================================================
# DotPlot
# ============================================================

p5 <- DotPlot(
  seu_dsb,
  assay = "CITE_residual",
  features = top_features,
  group.by = "cite_resid_clusters"
) +
  RotatedAxis()

ggsave(
  file.path(out_dir, "07_top_residual_marker_dotplot.pdf"),
  p5,
  width = 12,
  height = 6
)

# ============================================================
# Marker-high analysis
# ============================================================

cat("Running marker-high enrichment...\n")

DefaultAssay(seu_dsb) <- "CITE_DSB"

adt_raw <- GetAssayData(
  seu_dsb,
  assay = "CITE_DSB",
  layer = "data"
)

marker_high <- t(
  apply(
    adt_raw,
    1,
    function(x) {
      x > quantile(x, 0.90, na.rm = TRUE)
    }
  )
)

marker_high <- as.data.frame(t(marker_high))

colnames(marker_high) <- paste0(
  colnames(marker_high),
  "_high"
)

meta_high <- bind_cols(
  seu_dsb@meta.data,
  marker_high
)

# ============================================================
# Frequency summaries
# ============================================================

freq_list <- list()

for (marker in colnames(marker_high)) {
  
  tmp <- meta_high %>%
    group_by(
      .data[[sample_col]],
      .data[[celltype_col]]
    ) %>%
    summarise(
      frac_high = mean(.data[[marker]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(marker = marker)
  
  freq_list[[marker]] <- tmp
}

freq_df <- bind_rows(freq_list)

write.csv(
  freq_df,
  file.path(out_dir, "08_marker_high_frequencies.csv"),
  row.names = FALSE
)

# ============================================================
# Heatmap-style plot
# ============================================================

p6 <- ggplot(
  freq_df,
  aes(
    x = marker,
    y = .data[[celltype_col]],
    fill = frac_high
  )
) +
  geom_tile() +
  facet_wrap(
    as.formula(paste("~", sample_col))
  ) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Frac high"
  )

ggsave(
  file.path(out_dir, "09_marker_high_heatmap.pdf"),
  p6,
  width = 18,
  height = 10
)

# ============================================================
# RNA state correlations
# ============================================================

DefaultAssay(seu_dsb) <- "RNA"

# Example Hallmark-like simple signatures
ifn_genes <- c(
  "STAT1", "IFI6", "ISG15", "MX1", "IFIT1",
  "IFIT3", "OAS1", "IRF7"
)

myeloid_genes <- c(
  "LYZ", "S100A8", "S100A9", "FCN1",
  "CTSS", "TYMP"
)

mk_genes <- c(
  "PPBP", "PF4", "ITGA2B", "GP9",
  "NRGN", "TUBB1"
)

sig_list <- list(
  IFN = ifn_genes,
  MYELOID = myeloid_genes,
  MEGAKARYOCYTE = mk_genes
)

for (nm in names(sig_list)) {
  
  genes_use <- intersect(
    sig_list[[nm]],
    rownames(seu_dsb)
  )
  
  if (length(genes_use) < 3) next
  
  seu_dsb <- AddModuleScore(
    seu_dsb,
    features = list(genes_use),
    name = paste0(nm, "_score"),
    assay = "RNA"
  )
}

# ============================================================
# Correlations with proteins
# ============================================================

DefaultAssay(seu_dsb) <- "CITE_DSB"

adt_corr <- GetAssayData(
  seu_dsb,
  assay = "CITE_DSB",
  layer = "data"
)

corr_res <- list()

score_cols <- grep(
  "_score1$",
  colnames(seu_dsb@meta.data),
  value = TRUE
)

for (marker in rownames(adt_corr)) {
  
  for (score in score_cols) {
    
    cor_val <- suppressWarnings(
      cor(
        as.numeric(adt_corr[marker, ]),
        seu_dsb@meta.data[[score]],
        method = "spearman",
        use = "complete.obs"
      )
    )
    
    corr_res[[paste(marker, score, sep = "__")]] <- data.frame(
      marker = marker,
      score = score,
      spearman_cor = cor_val
    )
  }
}

corr_df <- bind_rows(corr_res)

write.csv(
  corr_df,
  file.path(out_dir, "10_RNA_program_CITE_correlations.csv"),
  row.names = FALSE
)

# ============================================================
# Top correlations plot
# ============================================================

top_corr <- corr_df %>%
  arrange(desc(abs(spearman_cor))) %>%
  slice_head(n = 40)

p7 <- ggplot(
  top_corr,
  aes(
    x = reorder(
      paste(marker, score, sep = " | "),
      spearman_cor
    ),
    y = spearman_cor
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spearman correlation"
  )

ggsave(
  file.path(out_dir, "11_top_RNA_CITE_correlations.pdf"),
  p7,
  width = 10,
  height = 12
)

# ============================================================
# Save object
# ============================================================

saveRDS(
  seu_dsb,
  file.path(out_dir, "seurat_residual_cite_analysis.rds")
)

cat("Done!\n")