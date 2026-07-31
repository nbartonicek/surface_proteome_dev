#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260522_VH01624_461_222JLJVNX"
project <- "LK1"   # change to "LK1" or "LK2"

sample_name <- paste0(project, "-GEX")

gex_dir <- file.path(
  "../results/cellranger",
  run,
  sample_name,
  "outs",
  "filtered_feature_bc_matrix"
)

hto_dir <- file.path(
  "../results/cite_seq_count",
  run,
  paste0("hto_counts_", project, "_emptydrops"),
  "umi_count"
)

doublet_dir <- file.path("../results/scDblFinder", run, sample_name)
out_dir <- file.path("../results/seurat_demux", run)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# HTO sample annotation
# ----------------------------
hto_names <- c(
  "HTO1-GTCAACTCTTTAGCG" = "MOLM13",
  "HTO2-TGATGGCCTATTGGG" = "HBDN206_MNpCT",
  "HTO3-TTCCGCCTCTCTTTG" = "HBDN392_AML_MDS",
  "HTO4-AGTAAGTTCAGCGTA" = "HBDN501_KMT2A"
)

#hto_names <- c(
#  "HTO1-GTCAACTCTTTAGCG" = "HBDN498_TP53",
#  "HTO2-TGATGGCCTATTGGG" = "APOP576_TP53",
#  "HTO3-TTCCGCCTCTCTTTG" = "HBDN376_TP53",
#  "HTO4-AGTAAGTTCAGCGTA" = "healthy_01"
#)

# ----------------------------
# Load GEX
# ----------------------------

seu<- readRDS(
  file = file.path(doublet_dir, "seurat_emptyDrops_RNA_scDblFinder_filtered.rds")
)



# ----------------------------
# Basic RNA QC
# ----------------------------

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

# ----------------------------
# Load HTO matrix
# ----------------------------

hto_mat <- readMM(file.path(hto_dir, "matrix.mtx.gz"))
hto_barcodes <- fread(
  file.path(hto_dir, "barcodes.tsv.gz"),
  header = FALSE
)$V1
hto_barcodes<-paste0(hto_barcodes,"-1")
hto_features <- fread(
  file.path(hto_dir, "features.tsv.gz"),
  header = FALSE
)$V1

rownames(hto_mat) <- hto_features
colnames(hto_mat) <- hto_barcodes

# remove unmapped if present
hto_mat <- hto_mat[!grepl("^unmapped$", rownames(hto_mat), ignore.case = TRUE), , drop = FALSE]

# rename HTO rows to biological sample names
new_hto_names <- hto_names[rownames(hto_mat)]
rownames(hto_mat) <- ifelse(
  is.na(new_hto_names),
  rownames(hto_mat),
  new_hto_names
)

# ----------------------------
# Match GEX and HTO barcodes
# ----------------------------

common_barcodes <- intersect(colnames(seu), colnames(hto_mat))

cat("GEX cells:", ncol(seu), "\n")
cat("HTO cells:", ncol(hto_mat), "\n")
cat("Shared cells:", length(common_barcodes), "\n")

seu <- subset(seu, cells = common_barcodes)
hto_mat <- hto_mat[, common_barcodes, drop = FALSE]

seu[["HTO"]] <- CreateAssayObject(counts = hto_mat)

# ----------------------------
# HTO demultiplexing
# ----------------------------

DefaultAssay(seu) <- "HTO"

seu <- NormalizeData(
  seu,
  assay = "HTO",
  normalization.method = "CLR",
  margin = 2
)

seu <- HTODemux(
  seu,
  assay = "HTO",
  positive.quantile = 0.99
)

# Add simple top-HTO call based on raw counts
hto_counts <- GetAssayData(seu, assay = "HTO", layer = "counts")

seu$top_HTO <- apply(hto_counts, 2, function(x) {
  rownames(hto_counts)[which.max(x)]
})

seu$top_HTO_count <- apply(hto_counts, 2, max)

# ----------------------------
# RNA processing / UMAP
# ----------------------------

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(
  seu,
  vars.to.regress = c("percent.mt"),
  verbose = FALSE
)
seu <- RunPCA(seu, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

# ----------------------------
# QC tables
# ----------------------------

demux_summary <- seu@meta.data %>%
  dplyr::count(HTO_classification.global, hash.ID, name = "n_cells")

sample_qc <- seu@meta.data %>%
  group_by(hash.ID, HTO_classification.global) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    median_percent_ribo = median(percent.ribo),
    median_nCount_HTO = median(nCount_HTO),
    median_nFeature_HTO = median(nFeature_HTO),
    .groups = "drop"
  )

write.csv(
  demux_summary,
  file.path(out_dir, "demux_summary.csv"),
  row.names = FALSE
)

write.csv(
  sample_qc,
  file.path(out_dir, "sample_qc_summary.csv"),
  row.names = FALSE
)

write.csv(
  seu@meta.data,
  file.path(out_dir, "cell_metadata_with_demux.csv")
)

# ----------------------------
# Plots
# ----------------------------

pdf(file.path(out_dir, "01_initial_qc_violin.pdf"), width = 12, height = 5)
print(
  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "hash.ID",
    pt.size = 0.05,
    ncol = 3
  )
)
dev.off()

pdf(file.path(out_dir, "02_hto_heatmap.pdf"), width = 10, height = 8)
  HTOHeatmap(
    seu,
    ncells = 5000
  )
dev.off()

pdf(file.path(out_dir, "03_umap_demux.pdf"), width = 14, height = 5)
p1 <- DimPlot(seu, group.by = "HTO_classification.global", label = TRUE) +
  ggtitle("HTO classification")

p2 <- DimPlot(seu, group.by = "hash.ID", label = TRUE) +
  ggtitle("Assigned sample")

p3 <- DimPlot(seu, group.by = "top_HTO", label = TRUE) +
  ggtitle("Top HTO by raw count")

print(p1 + p2 + p3)
dev.off()

pdf(file.path(out_dir, "04_umap_qc_features.pdf"), width = 14, height = 5)
p1 <- FeaturePlot(seu, features = "nCount_RNA") +
  ggtitle("RNA UMI counts")

p2 <- FeaturePlot(seu, features = "nFeature_RNA") +
  ggtitle("Detected genes")

p3 <- FeaturePlot(seu, features = "percent.mt") +
  ggtitle("Mitochondrial %")

print(p1 + p2 + p3)
dev.off()

pdf(file.path(out_dir, "05_hto_featureplots.pdf"), width = 12, height = 8)
print(
  FeaturePlot(
    seu,
    features = rownames(seu[["HTO"]]),
    ncol = 2
  )
)
dev.off()

# ----------------------------
# Save object
# ----------------------------

saveRDS(
  seu,
  file.path(out_dir, "LK1_GEX_HTO_demux_seurat.rds")
)

cat("\nDone.\n")
cat("Output written to:", out_dir, "\n")
cat("\nDemux summary:\n")
print(demux_summary)

cat("\nSample QC summary:\n")
print(sample_qc)



######## ADT

adt_dir <- file.path(
  "../results/cite_seq_count_fromRaw",
  run,
  "adt_counts_LK2_raw/umi_count"
)

adt_mat <- readMM(file.path(adt_dir, "matrix.mtx.gz"))

adt_barcodes <- fread(
  file.path(adt_dir, "barcodes.tsv.gz"),
  header = FALSE
)$V1

adt_features <- fread(
  file.path(adt_dir, "features.tsv.gz"),
  header = FALSE
)

# usually feature name in V2, fallback to V1
adt_feature_names <- if (ncol(adt_features) >= 2) {
  adt_features$V2
} else {
  adt_features$V1
}

adt_barcodes <- paste0(adt_barcodes, "-1")

rownames(adt_mat) <- make.unique(adt_feature_names)
colnames(adt_mat) <- adt_barcodes

# Match cells
common_adt <- intersect(colnames(seu), colnames(adt_mat))

cat("ADT shared cells:", length(common_adt), "\n")

adt_mat <- adt_mat[, common_adt, drop = FALSE]
adt_names_clean <- gsub(
  "-[ACGT]+$",
  "",
  rownames(adt_mat)
)

rownames(adt_mat) <- adt_names_clean
adt_counts <- Matrix::rowSums(adt_mat)

plot_df <- tibble(
  feature = names(adt_counts),
  counts = as.numeric(adt_counts)
) %>%
  arrange(desc(counts)) %>%
  mutate(
    feature = factor(feature, levels = feature)
  )

ggplot(plot_df, aes(x = feature, y = counts)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size= 5
    )
  ) +
  labs(
    x = "ADT",
    y = "Total reads",
    title = "ADT counts by feature"
  )
# subset Seurat object to shared cells
seu <- subset(seu, cells = common_adt)

# add ADT assay
seu[["ADT"]] <- CreateAssayObject(counts = adt_mat)

DefaultAssay(seu) <- "ADT"

seu <- NormalizeData(
  seu,
  normalization.method = "CLR",
  margin = 2
)

DefaultAssay(seu) <- "ADT"

seu$ADT_total <- Matrix::colSums(GetAssayData(seu, assay = "ADT", layer = "counts"))
seu$ADT_features <- Matrix::colSums(GetAssayData(seu, assay = "ADT", layer = "counts") > 0)

#######
DefaultAssay(seu) <- "ADT"

adt_counts <- GetAssayData(seu, assay = "ADT", layer = "data")

adt_means <- rowMeans(as.matrix(adt_counts))

top5_adt <- sort(adt_means, decreasing = TRUE)[1:7]
top5_adt<-top5_adt[c(1,4,5,6,7)]
print(top5_adt)
biotin_feature <- grep(
  "biotin",
  rownames(seu[["ADT"]]),
  ignore.case = TRUE,
  value = TRUE
)

features_to_plot <- unique(c(
  c("CD34","CD19","CD11a","ILT3","CD123"),
  biotin_feature
))

pdf(file.path(out_dir, "07_top_ADT_markers.pdf"), width = 12, height = 8)
print(
  FeaturePlot(
    seu,
    features = features_to_plot,
    reduction = "umap",
    ncol = 3
  )
)
dev.off()
