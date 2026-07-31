#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(dsb)
  library(pheatmap)
})

# ============================================================
# CITE-seq QC first, then DSB normalisation and matched plots
# ============================================================
# Workflow:
#   1. Load demultiplexed Seurat object.
#   2. Import ADT counts from emptyDrops CITE-seq-Count output.
#   3. Do basic CITE QC without DSB:
#        - raw ADT totals/features
#        - CLR normalisation
#        - marker summaries
#        - marker distributions
#        - RNA UMAP marker overlays
#        - heatmap by RNA cluster / demux sample
#   4. Select background droplets from raw ADT output.
#   5. Run DSB normalisation.
#   6. Repeat equivalent plots with DSB values.

# ----------------------------
# Paths / parameters
# ----------------------------

#run <- "260423_VH01624_453_222HWMYNX"
run <- "260528_VH01624_464_222K7VKNX"
sample_name <- "LK2-GEX"
sample_short <- "LK2"

# Use your most recent demux object. Change if needed.
demux_dir <- file.path("../results/demux_comparison", run, sample_name)
seurat_file <- file.path(demux_dir, "seurat_HTO_vireo_combined_positive_quantile_0.99_singlets_HTO_named.rds")

# Raw 10x GEX folder used only for background droplet RNA QC
raw_rna_dir <- file.path(
  "../results/cellranger_withbam", run, sample_name,
  "outs/raw_feature_bc_matrix"
)

# ADT counts for raw droplets, needed for DSB background
adt_raw_dir <- file.path(
  "../results/cite_seq_count",
  run,
  paste0("adt_counts_", sample_short, "_raw/umi_count")
)

emptydrops_dir <- file.path(
  "../results/emptydrops",
  run,
  sample_name
)

out_dir <- file.path(
  "../results/cite_seq_dsb",
  run,
  sample_name
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper functions
# ----------------------------

read_citeseq_count_matrix <- function(umi_count_dir, add_suffix = TRUE) {
  mat <- readMM(file.path(umi_count_dir, "matrix.mtx.gz"))
  
  barcodes <- fread(
    file.path(umi_count_dir, "barcodes.tsv.gz"),
    header = FALSE
  )$V1
  
  features <- fread(
    file.path(umi_count_dir, "features.tsv.gz"),
    header = FALSE
  )
  
  feature_names <- features$V1
  feature_names <- gsub("-[ACGT]+$", "", feature_names)
  feature_names <- make.unique(feature_names)
  
  if (add_suffix && !any(grepl("-1$", barcodes))) {
    barcodes <- paste0(barcodes, "-1")
  }
  
  rownames(mat) <- feature_names
  colnames(mat) <- barcodes
  
  mat <- mat[
    !grepl("^unmapped$", rownames(mat), ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  
  mat
}

make_marker_summary <- function(mat, value_name = "value") {
  tibble(
    marker = rownames(mat),
    mean = rowMeans(mat),
    median = apply(mat, 1, median),
    sd = apply(mat, 1, sd),
    pct_gt_0 = rowSums(mat > 0) / ncol(mat) * 100,
    pct_gt_1 = rowSums(mat > 1) / ncol(mat) * 100,
    pct_gt_3 = rowSums(mat > 3) / ncol(mat) * 100,
    assay_value = value_name
  ) %>%
    arrange(desc(mean))
}

choose_markers_to_plot <- function(summary_df, all_markers, n = 12) {
  top_markers <- summary_df %>%
    filter(!grepl("isotype|igg|control|unmapped", marker, ignore.case = TRUE)) %>%
    slice_max(mean, n = n) %>%
    pull(marker)
  
  biotin_marker <- grep("biotin", all_markers, ignore.case = TRUE, value = TRUE)
  unique(c(top_markers, biotin_marker))
}

plot_marker_distributions <- function(mat, markers, xlab, title) {
  plot_long <- as.data.frame(t(mat[markers, , drop = FALSE])) %>%
    rownames_to_column("barcode") %>%
    pivot_longer(-barcode, names_to = "marker", values_to = "value")
  
  ggplot(plot_long, aes(x = value)) +
    geom_density() +
    facet_wrap(~ marker, scales = "free_y", ncol = 5) +
    theme_bw() +
    scale_x_log10()+
    labs(
      x = xlab,
      y = "Density",
      title = title
    )
}

plot_heatmap_by_group <- function(mat, meta, group_col, markers, main, output_file,
                                  width = 10, height = 8) {
  stopifnot(all(colnames(mat) %in% rownames(meta)))
  
  heat_df <- cbind(
    meta[colnames(mat), , drop = FALSE],
    as.data.frame(t(mat[markers, , drop = FALSE]))
  ) %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(across(all_of(markers), median), .groups = "drop") %>%
    column_to_rownames(group_col)
  
  pdf(output_file, width = width, height = height)
  pheatmap(
    t(as.matrix(heat_df)),
    fontsize_row = 8,
    border_color = NA,
    main = main
  )
  dev.off()
}

# ----------------------------
# Load demultiplexed object
# ----------------------------

if (!file.exists(seurat_file)) {
  stop("Cannot find demux Seurat file: ", seurat_file)
}

seu <- readRDS(seurat_file)

cat("Loaded demuxed Seurat object:\n")
cat("Cells:", ncol(seu), "\n")
cat("Features:", nrow(seu), "\n")
cat("Assays:", paste(Assays(seu), collapse = ", "), "\n\n")

# ----------------------------
# Assign clean sample names from Vireo donor relabelled to HTO names
# ----------------------------

# ----------------------------
# Import emptyDrops ADT counts and add raw + CLR assays
# ----------------------------

adt_empty <- read_citeseq_count_matrix(adt_raw_dir, add_suffix = TRUE)
common_cells <- intersect(colnames(seu), colnames(adt_empty))

cat("ADT emptyDrops barcodes:", ncol(adt_empty), "\n")
cat("Shared demux/ADT cells:", length(common_cells), "\n\n")

if (length(common_cells) == 0) {
  stop("No shared cells between demux object and emptyDrops ADT matrix.")
}

seu <- subset(seu, cells = common_cells)
adt_empty <- adt_empty[, colnames(seu), drop = FALSE]

if ("ADT_raw" %in% Assays(seu)) seu[["ADT_raw"]] <- NULL
if ("ADT_CLR" %in% Assays(seu)) seu[["ADT_CLR"]] <- NULL
if ("CITE_DSB" %in% Assays(seu)) seu[["CITE_DSB"]] <- NULL

seu[["ADT_raw"]] <- CreateAssayObject(counts = adt_empty)
seu[["ADT_CLR"]] <- CreateAssayObject(counts = adt_empty)

DefaultAssay(seu) <- "ADT_CLR"
seu <- NormalizeData(
  seu,
  assay = "ADT_CLR",
  normalization.method = "CLR",
  margin = 2,
  verbose = FALSE
)

seu$ADT_raw_total <- Matrix::colSums(adt_empty)
seu$ADT_raw_features <- Matrix::colSums(adt_empty > 0)

# ----------------------------
# RNA processing / UMAP for QC overlays
# ----------------------------

DefaultAssay(seu) <- "RNA"

if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
}

if (!"percent.ribo" %in% colnames(seu@meta.data)) {
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
}

# Recompute a simple RNA UMAP if missing.
#if (!"umap" %in% Reductions(seu)) {
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
#seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
dims_use <- 1:min(30, ncol(Embeddings(seu, "pca")))
seu <- FindNeighbors(seu, dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)
#}

# ----------------------------
# 1) Regular CITE QC, before DSB
# ----------------------------

raw_summary <- make_marker_summary(adt_empty, value_name = "raw_umi")
write.csv(raw_summary, file.path(out_dir, "01_raw_ADT_marker_summary.csv"), row.names = FALSE)

clr_mat <- GetAssayData(seu, assay = "ADT_CLR", layer = "data")
clr_summary <- make_marker_summary(as.matrix(clr_mat), value_name = "CLR")
write.csv(clr_summary, file.path(out_dir, "02_CLR_ADT_marker_summary.csv"), row.names = FALSE)

markers_raw <- choose_markers_to_plot(raw_summary, rownames(adt_empty), n = 15)
markers_clr <- choose_markers_to_plot(clr_summary, rownames(clr_mat), n = 15)
markers_to_plot <- unique(c(markers_raw, markers_clr))
markers_to_plot <- intersect(markers_to_plot, rownames(adt_empty))

write.csv(
  tibble(marker = markers_to_plot),
  file.path(out_dir, "markers_selected_for_CITE_QC_plots.csv"),
  row.names = FALSE
)

# Cell-level QC violin
pdf(file.path(out_dir, "03_cell_QC_RNA_ADT_raw.pdf"), width = 12, height = 5)
print(
  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_raw_total"),
    pt.size = 0.05,
    group.by = "sample_name",
    ncol = 4
  )
)
dev.off()

# Raw ADT distributions
markers_to_plot_clr <- intersect(markers_to_plot, rownames(clr_mat))
pdf(file.path(out_dir, "04_raw_ADT_marker_distributions.pdf"), width = 12, height = 6)
print(
  plot_marker_distributions(
    as.matrix(adt_empty),
    markers_to_plot_clr[1:15],
    xlab = "Raw ADT UMI count",
    title = "Raw CITE marker distributions before DSB"
  )
)
dev.off()

# CLR ADT distributions

pdf(file.path(out_dir, "05_CLR_ADT_marker_distributions.pdf"), width = 12, height = 6)
print(
  plot_marker_distributions(
    as.matrix(clr_mat),
    markers_to_plot_clr[1:15],
    xlab = "CLR-normalised ADT",
    title = "CLR-normalised CITE marker distributions before DSB"
  )
)
dev.off()

# RNA UMAP QC
pdf(file.path(out_dir, "06_RNA_umap_QC.pdf"), width = 20, height = 4)

p1 <- DimPlot(seu, reduction = "umap", group.by = "sample_name", label = TRUE) +
  ggtitle("RNA clusters")

p2 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("RNA clusters")

p3 <- FeaturePlot(seu, reduction = "umap", features = "nCount_RNA") +
  ggtitle("RNA UMIs")

p4 <- FeaturePlot(seu, reduction = "umap", features = "percent.mt") +
  ggtitle("Mitochondrial %")

print(
  (p1 + p2 + p3 + p4) +
    patchwork::plot_layout(ncol = 4)
)
dev.off()


# CLR marker overlays using RNA UMAP
DefaultAssay(seu) <- "ADT_CLR"
pdf(file.path(out_dir, "07_CLR_ADT_featureplots_top_markers.pdf"), width = 14, height = 10)
print(
  FeaturePlot(
    seu,
    features = markers_to_plot[1:17],
    reduction = "umap",
    ncol = 4,
    order = TRUE
  )
)
dev.off()

# CLR marker overlays using RNA UMAP - discriminative markers
clr_mat <- GetAssayData(seu, assay = "ADT_CLR", layer = "data")

cluster_col <- "seurat_clusters"

markers_use <- rownames(clr_mat) %>%
  setdiff(grep("isotype|igg|control|unmapped", ., ignore.case = TRUE, value = TRUE))

meta <- seu@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, cluster = all_of(cluster_col))

clr_df <- as.data.frame(t(as.matrix(clr_mat[markers_use, , drop = FALSE]))) %>%
  rownames_to_column("barcode") %>%
  left_join(meta, by = "barcode") %>%
  filter(!is.na(cluster))

marker_cluster_medians <- clr_df %>%
  pivot_longer(
    cols = all_of(markers_use),
    names_to = "marker",
    values_to = "CLR"
  ) %>%
  group_by(marker, cluster) %>%
  summarise(median_CLR = median(CLR, na.rm = TRUE), .groups = "drop")

discriminative_markers <- marker_cluster_medians %>%
  group_by(marker) %>%
  summarise(
    dynamic_range = max(median_CLR, na.rm = TRUE) - min(median_CLR, na.rm = TRUE),
    sd_cluster_median = sd(median_CLR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(dynamic_range))

write.csv(
  discriminative_markers,
  file.path(out_dir, "CLR_discriminative_ADT_markers_by_RNA_cluster.csv"),
  row.names = FALSE
)

top_discriminative_markers <- discriminative_markers %>%
  slice_head(n = 16) %>%
  pull(marker)

DefaultAssay(seu) <- "ADT_CLR"
pdf(file.path(out_dir, "08_CLR_ADT_discriminative_markers.pdf"), width = 14, height = 10)
print(
  FeaturePlot(
    seu,
    features = top_discriminative_markers,
    reduction = "umap",
    ncol = 4,
    order = TRUE
  )
)
dev.off()

# Heatmap by RNA clusters using CLR values
cluster_col <- if ("seurat_clusters" %in% colnames(seu@meta.data)) {
  "seurat_clusters"
} else {
  colnames(seu@meta.data)[1]
}

markers_for_heatmap <- rownames(clr_mat)[
  !grepl("isotype|igg|control|biotin|unmapped", rownames(clr_mat), ignore.case = TRUE)
]

plot_heatmap_by_group(
  mat = as.matrix(clr_mat),
  meta = seu@meta.data,
  group_col = cluster_col,
  markers = markers_for_heatmap,
  main = "Median CLR-normalised CITE signal by RNA cluster",
  output_file = file.path(out_dir, "09_CLR_heatmap_by_RNA_cluster.pdf"),
  width = 6,
  height = 18
)

# Optional demux/sample heatmap if sampleID/hash.ID exists
sample_col <- case_when(
  "sample_name" %in% colnames(seu@meta.data) ~ "sample_name",
  "hash.ID" %in% colnames(seu@meta.data) ~ "hash.ID",
  "hto_donor" %in% colnames(seu@meta.data) ~ "hto_donor",
  TRUE ~ NA_character_
)

seu_clean <- subset(seu, subset = !(sample_name =="unassigned"))
clr_mat_clean <- GetAssayData(seu_clean, assay = "ADT_CLR", layer = "data")

if (!is.na(sample_col)) {
  
  n <- length(markers_for_heatmap)
  marker_split <- list(
    markers_for_heatmap[seq_len(ceiling(n / 2))],
    markers_for_heatmap[(ceiling(n / 2) + 1):n]
  )
  
  for (i in seq_along(marker_split)) {
    plot_heatmap_by_group(
      mat = as.matrix(clr_mat_clean),
      meta = seu_clean@meta.data,
      group_col = sample_col,
      markers = marker_split[[i]],
      main = "",
      output_file = file.path(
        out_dir,
        paste0("10_CLR_heatmap_by_", sample_col, "_part", i, ".pdf")
      ),
      width = 5,
      height = 12
    )
  }
}

######### plot of ordered CITE-seq 

# Use CLR-normalised ADT values
clr_mat <- GetAssayData(seu_clean, assay = "ADT_CLR", values = "data")
meta <- seu_clean@meta.data

markers_use <- intersect(markers_for_heatmap, rownames(clr_mat))
clr_use <- as.matrix(clr_mat[markers_use, , drop = FALSE])

# Marker-specific thresholds
marker_cutoffs <- apply(clr_use, 1, function(x) {
  median(x, na.rm = TRUE) + 2 * mad(x, na.rm = TRUE)
})

# Long format
df_long <- as_tibble(clr_use, rownames = "marker") %>%
  pivot_longer(-marker, names_to = "cell", values_to = "clr") %>%
  left_join(
    meta %>%
      rownames_to_column("cell") %>%
      select(cell, sample_name),
    by = "cell"
  ) %>%
  filter(!is.na(sample_name)) %>%
  mutate(
    cutoff = marker_cutoffs[marker],
    positive = clr > cutoff
  )

# Count positives per marker per sample
plot_df <- df_long %>%
  group_by(marker, sample_name) %>%
  summarise(
    n_positive = sum(positive, na.rm = TRUE),
    .groups = "drop"
  )

# Order markers by total positivity
marker_order <- plot_df %>%
  group_by(marker) %>%
  summarise(total = sum(n_positive)) %>%
  arrange(desc(total)) %>%
  pull(marker)

plot_df <- plot_df %>%
  mutate(marker = factor(marker, levels = marker_order))

# Plot
p <- ggplot(plot_df, aes(x = marker, y = n_positive, fill = sample_name)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    x = "CITE marker",
    y = "Number of positive cells",
    fill = "Sample",
    title = "CITE marker positivity by sample",
    subtitle = "Positive = CLR > marker median + 2 × MAD"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "11_CITE_marker_CLR_positive_cells_by_sample.pdf"),
  p,
  width = 15,
  height = 3
)


########## use different cutoff

# Long format
df_long <- as_tibble(clr_use, rownames = "marker") %>%
  pivot_longer(-marker, names_to = "cell", values_to = "clr") %>%
  left_join(
    meta %>%
      rownames_to_column("cell") %>%
      select(cell, sample_name),
    by = "cell"
  ) %>%
  filter(!is.na(sample_name)) %>%
  mutate(
    cutoff = 1,
    positive = clr > cutoff
  )

# Count positives per marker per sample
plot_df <- df_long %>%
  group_by(marker, sample_name) %>%
  summarise(
    n_positive = sum(positive, na.rm = TRUE),
    .groups = "drop"
  )

# Order markers by total positivity
marker_order <- plot_df %>%
  group_by(marker) %>%
  summarise(total = sum(n_positive)) %>%
  arrange(desc(total)) %>%
  pull(marker)

plot_df <- plot_df %>%
  mutate(marker = factor(marker, levels = marker_order))

# Plot
p <- ggplot(plot_df, aes(x = marker, y = n_positive, fill = sample_name)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    x = "CITE marker",
    y = "Number of positive cells",
    fill = "Sample",
    title = "CITE marker positivity by sample",
    subtitle = "Positive = CLR > 1"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "11a_CITE_marker_CLR_positive_cells_by_sample_CLR1.pdf"),
  p,
  width = 15,
  height = 3
)


# ----------------------------
# 2) Background droplet QC for DSB
# ----------------------------

adt_raw_all <- read_citeseq_count_matrix(adt_raw_dir, add_suffix = TRUE)

# Match protein set between cell matrix and raw matrix
common_proteins <- intersect(rownames(adt_empty), rownames(adt_raw_all))
adt_empty_dsb <- adt_empty[common_proteins, , drop = FALSE]
adt_raw_all <- adt_raw_all[common_proteins, , drop = FALSE]

# Load raw RNA to classify background droplets
rna_raw <- Read10X(raw_rna_dir)
if (is.list(rna_raw)) {
  # Cell Ranger v3+ raw_feature_bc_matrix can return list of assays.
  # Prefer Gene Expression.
  rna_raw <- rna_raw[[grep("Gene Expression|RNA", names(rna_raw), ignore.case = TRUE)[1]]]
}

stained_cells <- readLines(file.path(emptydrops_dir, "emptydrops_barcodes.txt"))
if (!any(grepl("-1$", stained_cells)) && any(paste0(stained_cells, "-1") %in% colnames(adt_raw_all))) {
  stained_cells <- paste0(stained_cells, "-1")
}

if (!any(grepl("-1$", colnames(rna_raw))) && any(paste0(colnames(rna_raw), "-1") %in% colnames(adt_raw_all))) {
  colnames(rna_raw) <- paste0(colnames(rna_raw), "-1")
}

shared_raw_barcodes <- intersect(colnames(rna_raw), colnames(adt_raw_all))
rna_raw <- rna_raw[, shared_raw_barcodes, drop = FALSE]
adt_raw_all <- adt_raw_all[, shared_raw_barcodes, drop = FALSE]

mtgene <- grep("^MT-", rownames(rna_raw), value = TRUE)

md <- data.frame(
  barcode = colnames(rna_raw),
  rna.size = log10(Matrix::colSums(rna_raw) + 1),
  n.gene = Matrix::colSums(rna_raw > 0),
  mt.prop = Matrix::colSums(rna_raw[mtgene, , drop = FALSE]) /
    pmax(Matrix::colSums(rna_raw), 1),
  prot.size = log10(Matrix::colSums(adt_raw_all) + 1),
  drop.class = ifelse(colnames(rna_raw) %in% stained_cells, "cell", "background"),
  row.names = colnames(rna_raw)
)

write.csv(md, file.path(out_dir, "11_raw_droplet_QC_metadata.csv"), row.names = FALSE)

pdf(file.path(out_dir, "12_raw_droplet_QC_background_vs_cells.pdf"), width = 10, height = 5)
print(
  ggplot(md, aes(x = log10(n.gene + 1), y = prot.size)) +
    geom_bin2d(bins = 200) +
    facet_wrap(~ drop.class) +
    theme_bw() +
    labs(
      x = "log10 detected RNA genes + 1",
      y = "log10 ADT UMIs + 1",
      title = "Raw droplet QC: cells versus background droplets"
    )
)
dev.off()

# Deliberately broad defaults. Adjust after inspecting plot 12.
background_drops <- md %>%
  filter(
    drop.class == "background",
    prot.size > 1.5,
    prot.size < 4.0,
    rna.size < 2.5
  ) %>%
  pull(barcode)

background_drops <- intersect(background_drops, colnames(adt_raw_all))

cat("Background droplets used for DSB:", length(background_drops), "\n")

if (length(background_drops) < 1000) {
  warning("Few background droplets selected. Consider relaxing background thresholds after checking 12_raw_droplet_QC_background_vs_cells.pdf")
}

background_adt <- adt_raw_all[, background_drops, drop = FALSE]

# Remove near-empty proteins from both matrices before DSB
protein_max <- sort(apply(adt_empty_dsb, 1, max))
write.csv(
  data.frame(
    protein = names(protein_max),
    max_raw_umi = as.numeric(protein_max)
  ),
  file.path(out_dir, "13_raw_ADT_marker_max_counts_in_cells.csv"),
  row.names = FALSE
)

low_signal_proteins <- names(protein_max)[protein_max <= 5]
cat("Low-signal ADTs removed before DSB:", length(low_signal_proteins), "\n")
print(low_signal_proteins)

if (length(low_signal_proteins) > 0) {
  adt_empty_dsb <- adt_empty_dsb[
    !rownames(adt_empty_dsb) %in% low_signal_proteins,
    ,
    drop = FALSE
  ]
  background_adt <- background_adt[
    !rownames(background_adt) %in% low_signal_proteins,
    ,
    drop = FALSE
  ]
}

isotype_controls <- grep(
  "isotype|igg|control",
  rownames(adt_empty_dsb),
  ignore.case = TRUE,
  value = TRUE
)

cat("Detected isotype/control ADTs:\n")
print(isotype_controls)

use_isotype <- length(isotype_controls) > 0

# ----------------------------
# 3) DSB normalisation
# ----------------------------

set.seed(123)

cells_dsb_norm <- DSBNormalizeProtein(
  cell_protein_matrix = as.matrix(adt_empty_dsb),
  empty_drop_matrix = as.matrix(background_adt),
  denoise.counts = TRUE,
  use.isotype.control = use_isotype,
  isotype.control.name.vec = if (use_isotype) isotype_controls else NULL,
  quantile.clipping = TRUE,
  quantile.clip = c(0.001, 0.999)
)

# Add DSB assay to the same demuxed/RNA UMAP object
seu[["CITE_DSB"]] <- CreateAssayObject(data = cells_dsb_norm)

# ----------------------------
# 4) Same plots after DSB
# ----------------------------

dsb_mat <- GetAssayData(seu, assay = "CITE_DSB", layer = "data")
dsb_summary <- make_marker_summary(as.matrix(dsb_mat), value_name = "DSB")
write.csv(dsb_summary, file.path(out_dir, "14_DSB_marker_summary.csv"), row.names = FALSE)
write.csv(seu@meta.data, file.path(out_dir, "15_cell_metadata_CITE_QC_DSB.csv"))

markers_dsb <- choose_markers_to_plot(dsb_summary, rownames(dsb_mat), n = 12)
markers_dsb <- intersect(unique(c(markers_to_plot, markers_dsb)), rownames(dsb_mat))

pdf(file.path(out_dir, "16_DSB_marker_distributions.pdf"), width = 12, height = 8)
print(
  plot_marker_distributions(
    as.matrix(dsb_mat),
    markers_dsb,
    xlab = "DSB-normalised ADT",
    title = "DSB-normalised CITE marker distributions"
  )
)
dev.off()

DefaultAssay(seu) <- "CITE_DSB"
pdf(file.path(out_dir, "17_DSB_featureplots_top_markers.pdf"), width = 14, height = 10)
print(
  FeaturePlot(
    seu,
    features = markers_dsb,
    reduction = "umap",
    ncol = 4,
    order = TRUE
  )
)
dev.off()

markers_for_dsb_heatmap <- rownames(dsb_mat)[
  !grepl("isotype|igg|control|biotin|unmapped", rownames(dsb_mat), ignore.case = TRUE)
]

plot_heatmap_by_group(
  mat = as.matrix(dsb_mat),
  meta = seu@meta.data,
  group_col = cluster_col,
  markers = markers_for_dsb_heatmap,
  main = "Median DSB-normalised CITE signal by RNA cluster",
  output_file = file.path(out_dir, "18_DSB_heatmap_by_RNA_cluster.pdf"),
  width = 10,
  height = 8
)

if (!is.na(sample_col)) {
  plot_heatmap_by_group(
    mat = as.matrix(dsb_mat),
    meta = seu@meta.data,
    group_col = sample_col,
    markers = markers_for_dsb_heatmap,
    main = paste0("Median DSB-normalised CITE signal by ", sample_col),
    output_file = file.path(out_dir, paste0("19_DSB_heatmap_by_", sample_col, ".pdf")),
    width = 10,
    height = 8
  )
}

# Direct side-by-side summary of CLR vs DSB marker behaviour
compare_summary <- clr_summary %>%
  select(marker, CLR_mean = mean, CLR_median = median, CLR_pct_gt_0 = pct_gt_0) %>%
  inner_join(
    dsb_summary %>%
      select(marker, DSB_mean = mean, DSB_median = median, DSB_pct_gt_0 = pct_gt_0),
    by = "marker"
  )

write.csv(compare_summary, file.path(out_dir, "20_CLR_vs_DSB_marker_summary.csv"), row.names = FALSE)

pdf(file.path(out_dir, "21_CLR_vs_DSB_marker_mean_scatter.pdf"), width = 6, height = 5)
print(
  ggplot(compare_summary, aes(x = CLR_mean, y = DSB_mean, label = marker)) +
    geom_point() +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 25) +
    theme_bw() +
    labs(
      x = "Mean CLR-normalised ADT",
      y = "Mean DSB-normalised ADT",
      title = "Marker-level comparison: CLR before DSB versus DSB"
    )
)
dev.off()

# ----------------------------
# Save object
# ----------------------------
DefaultAssay(seu) <- "RNA"

saveRDS(
  seu,
  file.path(out_dir, "seurat_demux_ADT_CLR_DSB_QC.rds")
)

cat("\nDone.\n")
cat("Output written to:", out_dir, "\n")
cat("Saved object:", file.path(out_dir, "seurat_demux_ADT_CLR_DSB_QC.rds"), "\n\n")

cat("Raw/CLR markers selected for plots:\n")
print(markers_to_plot)

cat("\nDSB markers selected for plots:\n")
print(markers_dsb)



####### comparisons

compare_summary <- compare_summary %>%
  mutate(delta_mean = DSB_mean - CLR_mean)

p <- ggplot(compare_summary, aes(x = CLR_mean, y = DSB_mean)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = compare_summary %>%
      arrange(desc(abs(delta_mean))) %>%
      slice_head(n = 20),
    aes(label = marker),
    size = 3
  ) +
  theme_bw() +
  labs(
    x = "CLR mean",
    y = "DSB mean",
    title = "CLR vs DSB marker signal",
    subtitle = "Dashed line = no change"
  )

p

compare_summary <- compare_summary %>%
  mutate(type = case_when(
    grepl("IgG|isotype|control", marker, ignore.case = TRUE) ~ "Isotype/control",
    TRUE ~ "Protein marker"
  ))

ggplot(compare_summary, aes(CLR_mean, DSB_mean, color = type)) +
  geom_point(alpha = 0.7) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(title = "DSB suppresses isotype/background signal")

compare_summary <- compare_summary %>%
  mutate(
    rank_CLR = rank(-CLR_mean),
    rank_DSB = rank(-DSB_mean)
  )

ggplot(compare_summary, aes(rank_CLR, rank_DSB)) +
  geom_point(alpha = 0.5) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Rank (CLR)",
    y = "Rank (DSB)",
    title = "Marker ranking shift after DSB"
  )

