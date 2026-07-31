#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 11 of 14
#
# Full DSB-vs-CLR evaluation - background suppression, marker specificity, leakage - plus composition by sample.
#
# Frozen for the lab archive 2026-07-31 from scripts/12.data_integration.R (mtime 2026-06-09).
# md5 of the original: d74562b8de6174aaec93ee5502c4abbd
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"

annotation_dir <- file.path("../results/seurat_annotated", run)

seurat_projected_rds <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

out_dir <- file.path("../results/data_integration", run)
projection_dir <- file.path(out_dir, "projection_umaps")

copykat_dir <- file.path("../results/seurat_demux", run, "copykat")
copykat_pred_file <- file.path(copykat_dir, "copykat_prediction_with_metadata.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Colours
# ----------------------------

celltype_order <- rev(c(
  "HSC MPP", "LMPP", "MEP", "Megakaryocyte Precursor",
  "GMP", "Early GMP", "Late GMP", "Cycling Progenitor",
  "EoBasoMast Precursor", "Pro-Monocyte", "Monocyte",
  "cDC", "pDC", "Early Lymphoid", "Pro-B", "Pre-B",
  "B", "Naive T", "CD4 Memory T", "CD8 Memory T",
  "NK", "Plasma Cell", "Early Erythroid", "Late Erythroid"
))

celltype_cols <- c(
  "HSC MPP" = "#1B9E77",
  "LMPP" = "#66C2A5",
  "MEP" = "#B2DF8A",
  "GMP" = "#33A02C",
  "Early GMP" = "#A6D854",
  "Late GMP" = "#006D2C",
  "Cycling Progenitor" = "#00441B",
  "EoBasoMast Precursor" = "#8DD3C7",
  "Megakocyte Precursor" = "#4DAF4A",
  "Megakaryocyte Precursor" = "#4DAF4A",
  "Monocyte" = "#E31A1C",
  "Pro-Monocyte" = "#FB6A4A",
  "cDC" = "#FD8D3C",
  "pDC" = "#FCBBA1",
  "Naive T" = "#2171B5",
  "CD4 Memory T" = "#6BAED6",
  "CD8 Memory T" = "#08519C",
  "NK" = "#54278F",
  "Early Lymphoid" = "#9E9AC8",
  "B" = "#3182BD",
  "Pre-B" = "#9ECAE1",
  "Pro-B" = "#C6DBEF",
  "Plasma Cell" = "#756BB1",
  "Early Erythroid" = "#F768A1",
  "Late Erythroid" = "#C51B8A",
  "Stromal" = "gray40"
)

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

copykat_cols <- c(
  "CNV_aberrant" = "#E41A1C",
  "CNV_neutral" = "#377EB8",
  "Not_called" = "grey80",
  "aneuploid" = "#E41A1C",
  "diploid" = "#377EB8",
  "not.defined" = "grey80"
)

granular_celltype_cols <- c(
  
  # ----------------------------
  # Stem / progenitor (greens)
  # ----------------------------
  
  "HSC" = "#1B9E77",
  "LMPP" = "#66C2A5",
  "MLP" = "#7BC87C",
  "MLP-II" = "#A1D99B",
  
  "MPP-MkEry" = "#74C476",
  "MPP-MyLy" = "#41AB5D",
  
  "MEP" = "#B2DF8A",
  
  "BFU-E" = "#C7E9C0",
  "CFU-E" = "#A1D99B",
  
  "Early GMP" = "#A6D854",
  "GMP-Cycle" = "#4DAF4A",
  "GMP-Mono" = "#238B45",
  "GMP-Neut" = "#006D2C",
  
  "Cycling Progenitor" = "#00441B",
  
  "EoBasoMast Precursor" = "#8DD3C7",
  "Megakaryocyte Precursor" = "#4DAF4A",
  
  # ----------------------------
  # Myeloid / DC (reds/oranges)
  # ----------------------------
  
  "CD14 Mono" = "#E31A1C",
  "CD16 Mono" = "#FB6A4A",
  
  "Early ProMono" = "#FC9272",
  "Late ProMono" = "#CB181D",
  
  "cDC1" = "#FD8D3C",
  "cDC2" = "#F16913",
  "Pre-cDC" = "#FDBB84",
  
  "pDC" = "#FCBBA1",
  "Pre-pDC" = "#FDD0A2",
  "Pre-pDC Cycling" = "#FEE6CE",
  
  "ASDC" = "#FDAE6B",
  
  # ----------------------------
  # Lymphoid (blues/purples)
  # ----------------------------
  
  "CLP" = "#9E9AC8",
  
  "CD4 Naive" = "#2171B5",
  "CD4 Central Memory" = "#4292C6",
  "CD4 Effector Memory" = "#6BAED6",
  "CD4 Regulatory" = "#9ECAE1",
  
  "CD8 Naive" = "#08519C",
  "CD8 Central Memory" = "#2171B5",
  "CD8 Effector Memory 1" = "#3182BD",
  "CD8 Effector Memory 2" = "#6BAED6",
  "CD8 Tissue Resident Memory" = "#9ECAE1",
  
  "T Proliferating" = "#6A51A3",
  
  "NK" = "#54278F",
  "NK CD56high" = "#756BB1",
  "NK Proliferating" = "#9E9AC8",
  
  "Immature B" = "#9ECAE1",
  "Large Pre-B" = "#C6DBEF",
  "Small Pre-B" = "#DEEBF7",
  
  "Pre-ProB" = "#C6DBEF",
  "Pro-B VDJ" = "#9ECAE1",
  
  "Mature B" = "#3182BD",
  
  "Plasma Cell" = "#756BB1",
  
  # ----------------------------
  # Erythroid (pinks)
  # ----------------------------
  
  "Pro-Erythroblast" = "#FBB4C4",
  "Basophilic Erythroblast" = "#F768A1",
  "Polychromatic Erythroblast" = "#DD3497",
  "Orthochromatic Erythroblast" = "#C51B8A",
  
  # ----------------------------
  # Other
  # ----------------------------
  
  "Stromal" = "gray40"
)

# ----------------------------
# Helper: broad lineage annotation
# ----------------------------

make_broad_lineage <- function(x) {
  case_when(
    x %in% c(
      "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
      "Cycling Progenitor", "EoBasoMast Precursor",
      "Megakaryocyte Precursor"
    ) ~ "Stem / progenitor",
    
    x %in% c("Monocyte", "Pro-Monocyte", "cDC", "pDC") ~ "Myeloid / DC",
    
    x %in% c(
      "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
      "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
    ) ~ "Lymphoid",
    
    x %in% c("Early Erythroid", "Late Erythroid") ~ "Erythroid",
    
    TRUE ~ "Other"
  )
}

# ----------------------------
# Load projected object
# ----------------------------

if (!file.exists(seurat_projected_rds)) {
  stop("Cannot find projected Seurat object: ", seurat_projected_rds)
}

seu <- readRDS(seurat_projected_rds)

cat("Loaded projected object\n")
cat("Cells:", ncol(seu), "\n")
cat("Assays:", paste(Assays(seu), collapse = ", "), "\n")
cat("Reductions:", paste(Reductions(seu), collapse = ", "), "\n\n")

stopifnot("predicted_CellType" %in% colnames(seu@meta.data))
stopifnot("sampleID" %in% colnames(seu@meta.data))

if (!"umap" %in% Reductions(seu)) {
  stop("No 'umap' reduction found in seu object.")
}

# ----------------------------
# Main UMAPs: detailed and broad annotations
# ----------------------------

pdf(file.path(projection_dir, "01_projected_predicted_celltypes_broad.pdf"), width = 14, height = 8)
p_broad <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "predicted_CellType_Broad",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  cols = celltype_cols
) +
  ggtitle("BoneMarrowMap projected cell types")
print(p_broad)
dev.off()

pdf(file.path(projection_dir, "02_projected_predicted_celltypes_granulated.pdf"), width = 28, height = 16)
print(
  DimPlot(
    seu,
    reduction = "umap",
    group.by = "predicted_CellType",
    label = TRUE,
    repel = TRUE,
    label.size = 5,
    cols = granular_celltype_cols
  ) +
    ggtitle("BoneMarrowMap broad annotations")
)
dev.off()

pdf(file.path(projection_dir, "03_projected_samples.pdf"), width = 10, height = 7)
print(
  DimPlot(
    seu,
    reduction = "umap",
    group.by = "sampleID",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Samples on projected UMAP")
)
dev.off()

if ("predicted_Pseudotime" %in% colnames(seu@meta.data)) {
  pdf(file.path(projection_dir, "04_projected_pseudotime.pdf"), width = 8, height = 6)
  print(
    FeaturePlot(
      seu,
      reduction = "umap",
      features = "predicted_Pseudotime",
      order = TRUE
    ) +
      ggtitle("Predicted pseudotime")
  )
  dev.off()
}

# ----------------------------
# Broad annotation by sample
# ----------------------------

pdf(file.path(projection_dir, "05_broad_annotations_split_by_sample.pdf"), width = 14, height = 10)
print(
  DimPlot(
    seu,
    reduction = "umap",
    group.by = "predicted_CellType_Broad",
    split.by = "sampleID",
    cols = celltype_cols,
    ncol = 2
  ) +
    ggtitle("Broad annotations by sample")
)
dev.off()

pdf(file.path(projection_dir, "06_detailed_annotations_split_by_sample.pdf"), width = 18, height = 12)
print(
  DimPlot(
    seu,
    reduction = "umap",
    group.by = "predicted_CellType",
    split.by = "sampleID",
    cols = granular_celltype_cols,
    ncol = 2
  ) +
    ggtitle("Detailed BoneMarrowMap annotations by sample")
)
dev.off()

####### Empty drop composition

emptydrop_barcodes_file <- paste0("../results/emptydrops/",run,"/",sample_name,"/emptydrops_barcodes.txt")
emptydrop_barcodes <- read.table(emptydrop_barcodes_file)

basic_seurat_barcodes_file <-  raw_dir <- paste0("../results/cellranger/",run,"/",sample_name,"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
seurat_barcodes <- read.table(basic_seurat_barcodes_file)

emptydrop_barcodes_unique <- emptydrop_barcodes$V1[!emptydrop_barcodes$V1 %in% seurat_barcodes$V1]

#seu_emptydrop_extra <- subset(seu, cells = emptydrop_barcodes_unique)

seu$emptydrop_unique <- colnames(seu) %in% emptydrop_barcodes_unique
pdf(file.path(projection_dir, "06a_detailed_annotations_split_by_sample_empty.pdf"), 
    width = 16, height = 6)
p<- DimPlot(
  seu,
  reduction = "umap",
  group.by = "emptydrop_unique",
  cols = c("lightgray", "black")
  #    cols = granular_celltype_cols,
) +
  ggtitle("Detailed BoneMarrowMap annotations by sample")

print(p_broad+p)
dev.off()

######## CITE QC

DefaultAssay(seu) <- "ADT_raw"

adt_counts <- GetAssayData(
  seu,
  assay = "ADT_raw",
  layer = "counts"
)

# ----------------------------
# Rank markers highest to lowest
# ----------------------------
marker_rank <- tibble(
  marker = rownames(adt_counts),
  total_counts = Matrix::rowSums(adt_counts),
  mean_count = Matrix::rowMeans(adt_counts),
  median_count = apply(as.matrix(adt_counts), 1, median)
) %>%
  arrange(desc(total_counts)) %>%
  mutate(
    rank = row_number(),
    tier = ntile(desc(total_counts), 3),
    tier = case_when(
      tier == 1 ~ "Tier_1_highest",
      tier == 2 ~ "Tier_2_middle",
      tier == 3 ~ "Tier_3_lowest"
    )
  )

write.csv(
  marker_rank,
  file.path(projection_dir, "ADT_raw_marker_expression_rank_tiers.csv"),
  row.names = FALSE
)

# ----------------------------
# Long format
# ----------------------------
adt_long <- as.data.frame(as.matrix(adt_counts)) %>%
  rownames_to_column("marker") %>%
  pivot_longer(
    cols = -marker,
    names_to = "cell",
    values_to = "count"
  ) %>%
  left_join(marker_rank %>% select(marker, tier, rank), by = "marker") %>%
  mutate(
    log_count = log10(count + 1),
    marker = factor(marker, levels = marker_rank$marker)
  )

# ----------------------------
# Plot each tier separately
# ----------------------------
for (this_tier in unique(marker_rank$tier)) {
  
  plot_df <- adt_long %>%
    filter(tier == this_tier)
  
  pdf(
    file.path(
      projection_dir,
      paste0("ADT_raw_marker_read_distributions_", this_tier, ".pdf")
    ),
    width = 14,
    height = 10
  )
  
  print(
    ggplot(plot_df, aes(x = log_count)) +
      geom_density(fill = "grey70", color = "black", linewidth = 0.2) +
      facet_wrap(~ marker, scales = "free_y", ncol = 5) +
      theme_bw() +
      labs(
        title = paste0("Raw ADT read distributions: ", gsub("_", " ", this_tier)),
        x = "log10(raw ADT counts + 1)",
        y = "Density"
      )
  )
  
  dev.off()
}

########## major cell types

dsb_assay <- "CITE_DSB"   # change if yours is named differently
group_col <- "predicted_CellType_Broad"

DefaultAssay(seu) <- dsb_assay

# ----------------------------
# Get DSB matrix
# ----------------------------
dsb_mat <- GetAssayData(
  seu,
  assay = dsb_assay,
  layer = "data"
)

# optional: remove isotype/background controls
markers_keep <- rownames(dsb_mat)[
  !grepl("IgG|isotype|Biotin|TotalSeq|Hash|HTO", rownames(dsb_mat), ignore.case = TRUE)
]

dsb_mat <- dsb_mat[markers_keep, , drop = FALSE]

# ----------------------------
# Long format with broad annotation
# ----------------------------
dsb_long <- as.data.frame(as.matrix(dsb_mat)) %>%
  rownames_to_column("marker") %>%
  pivot_longer(
    cols = -marker,
    names_to = "cell",
    values_to = "dsb"
  ) %>%
  left_join(
    seu@meta.data %>%
      rownames_to_column("cell") %>%
      select(cell, broad = all_of(group_col)),
    by = "cell"
  ) %>%
  filter(!is.na(broad))

# ----------------------------
# Top 3 CITE markers per broad annotation
# ----------------------------
top_markers <- dsb_long %>%
  group_by(broad, marker) %>%
  summarise(
    mean_dsb = mean(dsb, na.rm = TRUE),
    pct_positive = mean(dsb > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  group_by(broad) %>%
  slice_max(mean_dsb, n = 3, with_ties = FALSE) %>%
  ungroup()

top_features <- unique(top_markers$marker)

write.csv(
  top_markers,
  file.path(projection_dir, "top3_DSB_CITE_markers_per_broad_annotation.csv"),
  row.names = FALSE
)

dsb_cutoff <- 3   # try 2, 3, or 5

plot_df <- dsb_long %>%
  group_by(broad, marker) %>%
  summarise(
    mean_dsb = mean(dsb, na.rm = TRUE),
    pct_positive = mean(dsb > dsb_cutoff, na.rm = TRUE) * 100,
    .groups = "drop"
  )

top_markers <- plot_df %>%
  group_by(broad) %>%
  filter(pct_positive >= 10) %>%     # avoids markers positive in almost nobody
  slice_max(mean_dsb, n = 3, with_ties = FALSE) %>%
  ungroup()

top_features <- unique(top_markers$marker)

plot_df_top <- plot_df %>%
  filter(marker %in% top_features)
pdf(
  file.path(
    projection_dir,
    paste0("top3_DSB_CITE_markers_per_broad_annotation_cutoff_", dsb_cutoff, ".pdf")
  ),
  width = 12,
  height = 8
)

print(
  ggplot(plot_df_top, aes(x = marker, y = broad)) +
    geom_point(aes(size = pct_positive, color = mean_dsb)) +
    
    scale_color_gradient(
      low = "grey90",
      high = "darkblue",
      limits = c(0, 15),
      oob = squish
    ) +
    
    scale_size_continuous(
      range = c(0.5, 6),
      limits = c(0, 100)
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      title = paste0(
        "Top 3 DSB CITE markers per broad annotation, DSB > ",
        dsb_cutoff
      ),
      x = "CITE marker",
      y = "Broad annotation",
      color = "Mean DSB\n(capped at 15)",
      size = paste0("% DSB > ", dsb_cutoff)
    )
)

dev.off()


###### clustered annotation
cluster_mat <- plot_df_top %>%
  select(broad, marker, mean_dsb) %>%
  pivot_wider(
    names_from = marker,
    values_from = mean_dsb,
    values_fill = 0
  ) %>%
  column_to_rownames("broad") %>%
  as.matrix()

# row clustering = broad annotations
broad_order <- rownames(cluster_mat)[
  hclust(dist(cluster_mat))$order
]

# column clustering = CITE markers
marker_order <- colnames(cluster_mat)[
  hclust(dist(t(cluster_mat)))$order
]

plot_df_top <- plot_df_top %>%
  mutate(
    broad = factor(broad, levels = broad_order),
    marker = factor(marker, levels = marker_order)
  )

pdf(
  file.path(
    projection_dir,
    paste0("11.top3_DSB_CITE_markers_per_broad_annotation_cutoff_clustered_", dsb_cutoff, ".pdf")
  ),
  width = 12,
  height = 8
)
print(
  ggplot(plot_df_top, aes(x = marker, y = broad)) +
    geom_point(aes(size = pct_positive, color = mean_dsb)) +
    scale_color_gradient(
      low = "grey90",
      high = "darkblue",
      limits = c(0, 15),
      oob = squish
    ) +
    scale_size_continuous(
      range = c(0.5, 6),
      limits = c(0, 100)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      title = paste0(
        "Top 3 DSB CITE markers per broad annotation, DSB > ",
        dsb_cutoff
      ),
      x = "CITE marker",
      y = "Broad annotation",
      color = "Mean DSB\n(capped at 15)",
      size = paste0("% DSB > ", dsb_cutoff)
    )
)
dev.off()

########### granutlated clustering


group_col <- "predicted_CellType"

# ----------------------------
# Add detailed annotation
# ----------------------------
dsb_long_detail <- dsb_long %>%
  select(-broad) %>%
  left_join(
    seu@meta.data %>%
      rownames_to_column("cell") %>%
      select(cell, detail = all_of(group_col)),
    by = "cell"
  ) %>%
  filter(!is.na(detail))

# ----------------------------
# Summarise marker signal
# ----------------------------
plot_df_detail <- dsb_long_detail %>%
  group_by(detail, marker) %>%
  summarise(
    mean_dsb = mean(dsb, na.rm = TRUE),
    pct_positive = mean(dsb > dsb_cutoff, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# ----------------------------
# Top markers per detailed annotation
# ----------------------------
top_markers_detail <- plot_df_detail %>%
  group_by(detail) %>%
  filter(pct_positive >= 10) %>%
  slice_max(mean_dsb, n = 3, with_ties = FALSE) %>%
  ungroup()

top_features_detail <- unique(top_markers_detail$marker)

write.csv(
  top_markers_detail,
  file.path(
    projection_dir,
    "top3_DSB_CITE_markers_per_detailed_annotation.csv"
  ),
  row.names = FALSE
)

# ----------------------------
# Restrict to top features
# ----------------------------
plot_df_top_detail <- plot_df_detail %>%
  filter(marker %in% top_features_detail)

# =========================================================
# Cluster rows and columns
# =========================================================

cluster_mat_detail <- plot_df_top_detail %>%
  select(detail, marker, mean_dsb) %>%
  pivot_wider(
    names_from = marker,
    values_from = mean_dsb,
    values_fill = 0
  ) %>%
  column_to_rownames("detail") %>%
  as.matrix()

# row clustering = detailed annotations
detail_order <- rownames(cluster_mat_detail)[
  hclust(dist(cluster_mat_detail))$order
]

# column clustering = markers
marker_order_detail <- colnames(cluster_mat_detail)[
  hclust(dist(t(cluster_mat_detail)))$order
]

plot_df_top_detail <- plot_df_top_detail %>%
  mutate(
    detail = factor(detail, levels = rev(detail_order)),
    marker = factor(marker, levels = marker_order_detail)
  )

# =========================================================
# Plot
# =========================================================

pdf(
  file.path(
    projection_dir,
    paste0(
      "12.top3_DSB_CITE_markers_per_detailed_annotation_clustered_cutoff_",
      dsb_cutoff,
      ".pdf"
    )
  ),
  width = 14,
  height = 10
)

print(
  ggplot(plot_df_top_detail, aes(x = marker, y = detail)) +
    geom_point(aes(size = pct_positive, color = mean_dsb)) +
    
    scale_color_gradient(
      low = "grey90",
      high = "darkblue",
      limits = c(0, 15),
      oob = squish
    ) +
    
    scale_size_continuous(
      range = c(0.5, 6),
      limits = c(0, 100)
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 7
      ),
      axis.text.y = element_text(size = 7),
      panel.grid.minor = element_blank()
    ) +
    
    labs(
      title = paste0(
        "Top 3 DSB CITE markers per detailed annotation, DSB > ",
        dsb_cutoff
      ),
      x = "CITE marker",
      y = "Detailed annotation",
      color = "Mean DSB\n(capped at 15)",
      size = paste0("% DSB > ", dsb_cutoff)
    )
)

dev.off()


##########################



# ----------------------------
# Composition plots
# ----------------------------

meta_plot <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(sampleID))

broad_comp <- meta_plot %>%
  count(sampleID, predicted_CellType_Broad, name = "n") %>%
  group_by(sampleID) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  broad_comp,
  file.path(out_dir, "broad_celltype_composition_by_sample.csv"),
  row.names = FALSE
)

p_broad_comp <- ggplot(
  broad_comp,
  aes(x = sampleID, y = percent, fill = predicted_CellType_Broad)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cells (%)",
    fill = "Broad annotation",
    title = "Broad cell-type composition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(projection_dir, "07_broad_celltype_composition_by_sample.pdf"),
  p_broad_comp,
  width = 8,
  height = 4.5
)

detailed_comp <- meta_plot %>%
  count(sampleID, predicted_CellType, name = "n") %>%
  group_by(sampleID) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  detailed_comp,
  file.path(out_dir, "detailed_celltype_composition_by_sample.csv"),
  row.names = FALSE
)

p_detailed_comp <- ggplot(
  detailed_comp,
  aes(x = sampleID, y = percent, fill = predicted_CellType)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = granular_celltype_cols, drop = FALSE, na.value = "grey80") +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cells (%)",
    fill = "Detailed annotation",
    title = "Detailed BoneMarrowMap composition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "08_detailed_celltype_composition_by_sample.pdf"),
  p_detailed_comp,
  width = 12,
  height = 5
)

############ cluster by both cite and RNA

DefaultAssay(seu) <- "CITE_DSB"

# Remove isotypes / controls
prots <- rownames(seu[["CITE_DSB"]])

prots <- prots[!grepl("IgG|isotype|Control", prots, ignore.case = TRUE)]

# Optionally restrict to variable proteins
# prots <- head(prots, 30)

# ----------------------------
# Neighbors + clustering
# ----------------------------
seu_dsb <- FindNeighbors(
  object = seu,
  assay = "CITE_DSB",
  features = prots,
  dims = NULL,          # critical: no PCA
  k.param = 30,
  graph.name = "CITE_snn",
  verbose = FALSE
)

seu_dsb <- FindClusters(
  object = seu_dsb,
  graph.name = "CITE_snn",
  resolution = 0.8,
  algorithm = 3,
  verbose = FALSE
)

# ----------------------------
# Optional UMAP on proteins
# ----------------------------
seu_dsb <- RunUMAP(
  seu_dsb,
  assay = "CITE_DSB",
  features = prots,
  reduction.name = "umap_cite",
  verbose = FALSE
)

DefaultAssay(seu_dsb) <- "CITE_DSB"

p <- DimPlot(
  seu_dsb,
  reduction = "umap_cite",
  group.by = "predicted_CellType_Broad",
  split.by = "sample_name",
  cols = celltype_cols,
  ncol = 2
) +
  ggtitle("DSB-normalised CITE UMAP by sample")


ggsave(
  file.path(out_dir, "09_DSB_UMAP_by_sample.pdf"),
  p,
  width = 18,
  height = 8
)

prots <- prots[!grepl("IgG|isotype|Control", prots, ignore.case = TRUE)]

# Optionally restrict to variable proteins
# prots <- head(prots, 30)

# ----------------------------
# Neighbors + clustering
# ----------------------------
seu_dsb <- FindNeighbors(
  object = seu_dsb,
  assay = "CITE_DSB",
  features = prots,
  dims = NULL,          # critical: no PCA
  k.param = 30,
  graph.name = "CITE_snn",
  verbose = FALSE
)

seu_dsb <- FindClusters(
  object = seu_dsb,
  graph.name = "CITE_snn",
  resolution = 0.8,
  algorithm = 3,
  verbose = FALSE
)

# ----------------------------
# Optional UMAP on proteins
# ----------------------------
seu_dsb <- RunUMAP(
  seu_dsb,
  assay = "CITE_DSB",
  features = prots,
  reduction.name = "umap_cite",
  verbose = FALSE
)

p<- DimPlot(seu_dsb, reduction = "umap_cite", group.by = "seurat_clusters")

ggsave(
  file.path(out_dir, "10_DSB_UMAP_by_sample_variable_proteins.pdf"),
  p,
  width = 18,
  height = 8
)


######### WNN integration


# ----------------------------
# Parameters
# ----------------------------
sample_col <- "sample_name"  # change if needed
celltype_col <- "predicted_CellType_Broad"

wnn_out <- file.path(out_dir, "09_WNN_RNA_DSB")
dir.create(wnn_out, showWarnings = FALSE, recursive = TRUE)


# RNA processing
DefaultAssay(seu_dsb) <- "RNA"

#seu_dsb <- NormalizeData(seu_dsb, verbose = FALSE)
#seu_dsb <- FindVariableFeatures(seu_dsb, nfeatures = 3000, verbose = FALSE)
#seu_dsb <- ScaleData(seu_dsb, verbose = FALSE)
#seu_dsb <- RunPCA(
#  seu_dsb,
#  assay = "RNA",
#  reduction.name = "pca",
#  npcs = 50,
#  verbose = FALSE
#)

# CITE DSB processing
DefaultAssay(seu_dsb) <- "CITE_DSB"

# DSB values are already normalised, so just scale + PCA
cite_features <- rownames(seu_dsb[["CITE_DSB"]])

seu_dsb <- ScaleData(seu_dsb, verbose = FALSE)
seu_dsb <- RunPCA(
  seu_dsb,
  features = cite_features,
  assay = "CITE_DSB",
  reduction.name = "pca_dsb",
  reduction.key = "DSBPC_",
  npcs = 30,
  verbose = FALSE
)

# WNN integration: RNA + DSB protein
seu_dsb <- FindMultiModalNeighbors(
  seu_dsb,
  reduction.list = list("pca", "pca_dsb"),
  dims.list = list(1:30, 1:20),
  modality.weight.name = c("RNA.weight", "DSB.weight"),
  verbose = FALSE
)

seu_dsb <- FindClusters(
  seu_dsb,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.6,
  cluster.name = "wnn_clusters",
  verbose = FALSE
)

seu_dsb <- RunUMAP(
  seu_dsb,
  nn.name = "weighted.nn",
  reduction.name = "umap_wnn",
  reduction.key = "wnnUMAP_",
  verbose = FALSE
)

# ----------------------------
# 2. Cluster colours
# ----------------------------

wnn_levels <- sort(unique(as.character(seu_dsb$wnn_clusters)))

cluster_cols <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(wnn_levels)),
  wnn_levels
)

# ----------------------------
# 3. Plot WNN clusters + broad cell types
# ----------------------------

p_clusters <- DimPlot(
  seu_dsb,
  reduction = "umap_wnn",
  group.by = "wnn_clusters",
  cols = cluster_cols,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("RNA + DSB CITE WNN clusters")

ggsave(
  file.path(wnn_out, "01_WNN_RNA_DSB_clusters.pdf"),
  p_clusters,
  width = 8,
  height = 6
)

p_clusters_sample <- DimPlot(
  seu_dsb,
  reduction = "umap_wnn",
  group.by = "wnn_clusters",
  split.by = sample_col,
  cols = cluster_cols,
  label = TRUE,
  repel = TRUE,
  ncol = 2
) +
  ggtitle("RNA + DSB CITE WNN clusters by sample")

ggsave(
  file.path(wnn_out, "02_WNN_RNA_DSB_clusters_by_sample.pdf"),
  p_clusters_sample,
  width = 15,
  height = 9
)

p_celltype_sample <- DimPlot(
  seu_dsb,
  reduction = "umap_wnn",
  group.by = celltype_col,
  split.by = sample_col,
  cols = celltype_cols,
  ncol = 2
) +
  ggtitle("Broad cell types on RNA + DSB CITE WNN UMAP")

ggsave(
  file.path(wnn_out, "03_WNN_RNA_DSB_broad_celltypes_by_sample.pdf"),
  p_celltype_sample,
  width = 15,
  height = 9
)

# ----------------------------
# 4. Call new cluster labels by majority broad cell type
# ----------------------------

cluster_labels <- seu_dsb@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  count(wnn_clusters, .data[[celltype_col]], name = "n") %>%
  group_by(wnn_clusters) %>%
  mutate(freq = n / sum(n)) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    wnn_cluster_label = paste0(
      "WNN", wnn_clusters,
      "_",
      .data[[celltype_col]],
      "_",
      round(freq * 100),
      "pct"
    )
  )

label_map <- setNames(cluster_labels$wnn_cluster_label, cluster_labels$wnn_clusters)

seu_dsb$wnn_cluster_label <- as.character(
  label_map[as.character(seu_dsb$wnn_clusters)]
)

write.csv(
  cluster_labels,
  file.path(wnn_out, "04_WNN_cluster_majority_celltype_labels.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Find largest labelled cell type groups per sample
# ----------------------------

top_groups <- seu_dsb@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(!is.na(.data[[sample_col]]), !is.na(.data[[celltype_col]])) %>%
  count(.data[[sample_col]], .data[[celltype_col]], name = "n") %>%
  group_by(.data[[sample_col]]) %>%
  slice_max(n, n = 3, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_groups,
  file.path(wnn_out, "05_largest_broad_celltype_groups_per_sample.csv"),
  row.names = FALSE
)

# ----------------------------
# 6. Within each sample + broad cell type, find CITE markers
#    that split WNN clusters
# ----------------------------

DefaultAssay(seu_dsb) <- "CITE_DSB"

marker_list <- list()

for (i in seq_len(nrow(top_groups))) {
  
  smp <- top_groups[[sample_col]][i]
  ct  <- top_groups[[celltype_col]][i]
  
  cells_use <- rownames(seu_dsb@meta.data)[
    seu_dsb@meta.data[[sample_col]] == smp &
      seu_dsb@meta.data[[celltype_col]] == ct
  ]
  
  obj_sub <- subset(seu_dsb, cells = cells_use)
  
  # Need at least 2 WNN clusters represented
  if (length(unique(obj_sub$wnn_clusters)) < 2) next
  
  Idents(obj_sub) <- "wnn_clusters"
  
  markers <- FindAllMarkers(
    obj_sub,
    assay = "CITE_DSB",
    slot = "data",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.15,
    test.use = "wilcox"
  )
  
  if (nrow(markers) == 0) next
  
  markers <- markers %>%
    mutate(
      sample = smp,
      broad_celltype = ct
    ) %>%
    arrange(sample, broad_celltype, cluster, desc(avg_log2FC))
  
  marker_list[[paste(smp, ct, sep = "__")]] <- markers
}

cite_markers <- bind_rows(marker_list)

write.csv(
  cite_markers,
  file.path(wnn_out, "06_CITE_DSB_markers_splitting_largest_celltypes_by_sample.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Plot top CITE markers per sample / cell-type group
# ----------------------------

plot_dir <- file.path(wnn_out, "CITE_marker_featureplots_by_sample")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

top_marker_tbl <- cite_markers %>%
  group_by(sample, broad_celltype) %>%
  slice_max(avg_log2FC, n = 8, with_ties = FALSE) %>%
  ungroup()

for (nm in unique(paste(top_marker_tbl$sample, top_marker_tbl$broad_celltype, sep = "__"))) {
  
  smp <- sub("__.*$", "", nm)
  ct  <- sub("^.*__", "", nm)
  
  markers_to_plot <- top_marker_tbl %>%
    filter(sample == smp, broad_celltype == ct) %>%
    pull(gene) %>%
    unique()
  
  markers_to_plot <- intersect(
    markers_to_plot,
    rownames(seu_dsb[["CITE_DSB"]])
  )
  
  if (length(markers_to_plot) == 0) next
  
  cells_use <- rownames(seu_dsb@meta.data)[
    seu_dsb@meta.data[[sample_col]] == smp &
      seu_dsb@meta.data[[celltype_col]] == ct
  ]
  
  obj_sub <- subset(seu_dsb, cells = cells_use)
  
  DefaultAssay(obj_sub) <- "CITE_DSB"
  
  p_feat <- FeaturePlot(
    obj_sub,
    reduction = "umap_wnn",
    features = markers_to_plot,
    slot = "data",
    order = TRUE,
    ncol = 4,
    keep.scale = "feature"
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("CITE_markers_", make.names(smp), "_", make.names(ct), ".pdf")
    ),
    p_feat,
    width = 14,
    height = 10
  )
}



# ----------------------------
# Load CopyKAT calls and add to object
# ----------------------------


demux_out_dir <- file.path("../results/seurat_demux", run)
demux_seurat_file <- file.path(demux_out_dir, "LK1_GEX_HTO_demux_seurat.rds")

copykat_out <- file.path(demux_out_dir, "copykat")

copykat_meta <- read.csv(
  file.path(copykat_out, "copykat_prediction_with_metadata.csv"),
  stringsAsFactors = FALSE
)


# inspect columns
colnames(copykat_meta)

# ----------------------------
# Make sure barcode column exists
# ----------------------------
# replace "cell" if your barcode column has another name
copykat_meta <- copykat_meta %>%
  rename(cell = 1)

# ----------------------------
# Keep useful columns
# ----------------------------
copykat_meta_sub <- copykat_meta %>%
  select(
    cell,
    copykat_prediction,
    copykat_conf,
    everything()
  )

# ----------------------------
# Match to Seurat object
# ----------------------------
common_cells <- intersect(
  colnames(seu),
  copykat_meta_sub$cell
)

length(common_cells)

# subset/reorder
copykat_meta_sub <- copykat_meta_sub %>%
  filter(cell %in% common_cells) %>%
  arrange(match(cell, colnames(seu)))

# ----------------------------
# Add metadata to Seurat
# ----------------------------
meta_to_add <- copykat_meta_sub %>%
  column_to_rownames("cell")

seu <- AddMetaData(
  seu,
  metadata = meta_to_add
)

# ----------------------------
# Quick checks
# ----------------------------
table(seu$copykat_prediction, useNA = "ifany")

DimPlot(
  seu,
  reduction = "umap",
  group.by = "copykat_prediction",
  cols = c(
    "aneuploid" = "red",
    "diploid" = "grey70",
    "not.defined" = "black"
  )
)

seu$copykat_call <- "Not_called"
seu$copykat_malignancy <- "Not_called"

if (file.exists(copykat_pred_file)) {
  
  copykat_pred <- fread(copykat_pred_file) %>%
    as_tibble()
  
  cell_col <- intersect(c("cell", "barcode", "cell.names"), colnames(copykat_pred))[1]
  
  if (is.na(cell_col)) {
    warning("Could not find cell column in CopyKAT table. Skipping CopyKAT annotation.")
  } else {
    
    copykat_pred <- copykat_pred %>%
      rename(cell = all_of(cell_col))
    
    if (!"copykat_call" %in% colnames(copykat_pred)) {
      if ("copykat.pred" %in% colnames(copykat_pred)) {
        copykat_pred <- copykat_pred %>%
          rename(copykat_call = copykat.pred)
      } else {
        stop("CopyKAT table found but no copykat_call/copykat.pred column.")
      }
    }
    
    copykat_vec <- setNames(copykat_pred$copykat_call, copykat_pred$cell)
    
    common_copykat <- intersect(colnames(seu), names(copykat_vec))
    cat("\nCopyKAT cells matched:", length(common_copykat), "\n")
    
    seu$copykat_call[common_copykat] <- copykat_vec[common_copykat]
    
    seu$copykat_malignancy <- case_when(
      seu$copykat_call == "aneuploid" ~ "CNV_aberrant",
      seu$copykat_call == "diploid" ~ "CNV_neutral",
      seu$copykat_call %in% c("Not_called", NA) ~ "Not_called",
      TRUE ~ as.character(seu$copykat_call)
    )
  }
  
} else {
  warning("CopyKAT prediction file not found: ", copykat_pred_file)
}

seu$copykat_malignancy <- factor(
  seu$copykat_malignancy,
  levels = c("CNV_aberrant", "CNV_neutral", "Not_called")
)

cat("\nCopyKAT table:\n")
print(table(seu$copykat_malignancy, useNA = "ifany"))

seu_copykat_plot <- if ("mapping_error_QC" %in% colnames(seu@meta.data)) {
  subset(seu, subset = mapping_error_QC == "Pass")
} else {
  seu
}

# ----------------------------
# CopyKAT on joint UMAP
# ----------------------------

pdf(file.path(projection_dir, "09_copykat_on_joint_umap.pdf"), width = 8, height = 6)
print(
  DimPlot(
    seu_copykat_plot,
    reduction = "umap",
    group.by = "copykat_malignancy",
    cols = copykat_cols
  ) +
    ggtitle("CopyKAT calls on projected UMAP")
)
dev.off()

pdf(file.path(projection_dir, "10_copykat_on_umap_split_by_sample.pdf"), width = 14, height = 10)
print(
  DimPlot(
    seu_copykat_plot,
    reduction = "umap",
    group.by = "copykat_malignancy",
    split.by = "sampleID",
    cols = copykat_cols,
    ncol = 2
  ) +
    ggtitle("CopyKAT calls by sample")
)
dev.off()

# ----------------------------
# CopyKAT composition by sample and broad annotation
# ----------------------------

copykat_summary <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(sampleID)) %>%
  count(sampleID, predicted_CellType_Broad, copykat_malignancy, name = "n") %>%
  group_by(sampleID, predicted_CellType_Broad) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  copykat_summary,
  file.path(out_dir, "copykat_by_sample_and_broad_annotation.csv"),
  row.names = FALSE
)

p_copykat_comp <- copykat_summary %>%
  ggplot(aes(x = predicted_CellType_Broad, y = percent, fill = copykat_malignancy)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ sampleID) +
  scale_fill_manual(values = copykat_cols, drop = FALSE) +
  theme_bw() +
  labs(
    x = "Broad annotation",
    y = "Cells (%)",
    fill = "CopyKAT",
    title = "CopyKAT calls within broad BoneMarrowMap annotations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "11_copykat_by_sample_and_broad_annotation.pdf"),
  p_copykat_comp,
  width = 12,
  height = 5
)

# ----------------------------
# Save integrated object
# ----------------------------

saveRDS(
  seu,
  file.path(out_dir, "LK1_projected_CITE_DSB_CopyKAT_integrated.rds")
)

write.csv(
  seu@meta.data,
  file.path(out_dir, "LK1_projected_CITE_DSB_CopyKAT_integrated_metadata.csv")
)

cat("\nDone.\n")
cat("Output directory:", out_dir, "\n")
cat("Integrated object:", file.path(out_dir, "LK1_projected_CITE_DSB_CopyKAT_integrated.rds"), "\n")
