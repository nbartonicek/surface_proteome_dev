#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 03 of 14
#
# Add vireo donor calls to the object, cross-tabulate against HTO, and run doublet detection.
#
# Frozen for the lab archive 2026-07-31 from scripts/backup/6.demux.R (mtime 2026-05-05).
# md5 of the original: 8d35cb95e45928ba6d286bff8529b0e0
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(data.table)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(RColorBrewer)
})

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"

positive_quantile <- 0.99
pq_tag <- paste0("positive_quantile_", positive_quantile)

out_dir <- file.path("../results/scDblFinder", run, sample_name)

demux_out_dir <- file.path("../results/demux_comparison", run, sample_name)
dir.create(demux_out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load existing Seurat object
# ----------------------------

seu <- readRDS(
  file.path(out_dir, "seurat_emptyDrops_RNA_scDblFinder_nonfiltered.rds")
)

cat("Loaded Seurat object:\n")
cat("Cells:", ncol(seu), "\n")
cat("Genes:", nrow(seu), "\n\n")

# ----------------------------
# Load Vireo donor calls
# ----------------------------

vireo_dir <- file.path(
  "../results/genotype_demux",
  run,
  sample_name,
  "vireo"
)

donors <- read_tsv(
  file.path(vireo_dir, "donor_ids.tsv"),
  show_col_types = FALSE
)

cat("Vireo donor table columns:\n")
print(colnames(donors))

if ("cell" %in% colnames(donors)) {
  donors <- donors %>% rename(barcode = cell)
} else if ("cell_id" %in% colnames(donors)) {
  donors <- donors %>% rename(barcode = cell_id)
} else if (!"barcode" %in% colnames(donors)) {
  stop("Could not find barcode column in donor_ids.tsv")
}

if ("donor_id" %in% colnames(donors)) {
  donors <- donors %>% rename(vireo_donor = donor_id)
} else if ("donor" %in% colnames(donors)) {
  donors <- donors %>% rename(vireo_donor = donor)
} else if ("best_singlet" %in% colnames(donors)) {
  donors <- donors %>% rename(vireo_donor = best_singlet)
} else {
  stop("Could not find donor assignment column in donor_ids.tsv")
}

keep_cols <- intersect(
  c("barcode", "vireo_donor", "doublet_prob", "prob_doublet", "prob_max", "n_vars"),
  colnames(donors)
)

donors <- donors %>%
  select(all_of(keep_cols)) %>%
  mutate(barcode = as.character(barcode))

# ----------------------------
# Add Vireo calls to Seurat metadata
# ----------------------------

seu$vireo_donor <- donors$vireo_donor[
  match(colnames(seu), donors$barcode)
]

if ("doublet_prob" %in% colnames(donors)) {
  seu$vireo_doublet_prob <- donors$doublet_prob[
    match(colnames(seu), donors$barcode)
  ]
} else if ("prob_doublet" %in% colnames(donors)) {
  seu$vireo_doublet_prob <- donors$prob_doublet[
    match(colnames(seu), donors$barcode)
  ]
} else {
  seu$vireo_doublet_prob <- NA_real_
}

vireo_doublet_threshold <- 0.5

seu$vireo_class <- case_when(
  is.na(seu$vireo_donor) ~ "unassigned",
  !is.na(seu$vireo_doublet_prob) &
    seu$vireo_doublet_prob >= vireo_doublet_threshold ~ "doublet",
  TRUE ~ "singlet"
)

# ----------------------------
# Load HTO counts
# ----------------------------

hto_dir <- file.path(
  "../results/cite_seq_count",
  run,
  "hto_counts_LK1_emptydrops/umi_count"
)

hto_names <- c(
  "HTO1-GTCAACTCTTTAGCG" = "MOLM13",
  "HTO2-TGATGGCCTATTGGG" = "HBDN206_MNpCT",
  "HTO3-TTCCGCCTCTCTTTG" = "HBDN392_AML_MDS",
  "HTO4-AGTAAGTTCAGCGTA" = "HBDN501_AML_KMT2A"
)

hto_mat <- readMM(file.path(hto_dir, "matrix.mtx.gz"))

hto_barcodes <- fread(
  file.path(hto_dir, "barcodes.tsv.gz"),
  header = FALSE
)$V1

hto_features <- fread(
  file.path(hto_dir, "features.tsv.gz"),
  header = FALSE
)$V1

if (!any(hto_barcodes %in% colnames(seu)) &&
    any(paste0(hto_barcodes, "-1") %in% colnames(seu))) {
  hto_barcodes <- paste0(hto_barcodes, "-1")
}

rownames(hto_mat) <- hto_features
colnames(hto_mat) <- hto_barcodes

hto_mat <- hto_mat[
  !grepl("^unmapped$", rownames(hto_mat), ignore.case = TRUE),
  ,
  drop = FALSE
]

new_hto_names <- hto_names[rownames(hto_mat)]

rownames(hto_mat) <- ifelse(
  is.na(new_hto_names),
  rownames(hto_mat),
  new_hto_names
)

rownames(hto_mat) <- make.unique(rownames(hto_mat))

common_barcodes <- intersect(colnames(seu), colnames(hto_mat))

cat("Seurat cells:", ncol(seu), "\n")
cat("HTO cells:", ncol(hto_mat), "\n")
cat("Shared cells:", length(common_barcodes), "\n")

if (length(common_barcodes) == 0) {
  stop("No shared barcodes between Seurat object and HTO matrix. Check barcode suffixes, sample name, and HTO path.")
}

seu <- subset(seu, cells = common_barcodes)
hto_mat <- hto_mat[, colnames(seu), drop = FALSE]

if ("HTO" %in% Assays(seu)) {
  seu[["HTO"]] <- NULL
}

seu[["HTO"]] <- CreateAssayObject(counts = hto_mat)

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
  positive.quantile = positive_quantile
)

hto_counts <- GetAssayData(seu, assay = "HTO", layer = "counts")

seu$top_HTO <- apply(hto_counts, 2, function(x) {
  rownames(hto_counts)[which.max(x)]
})

seu$top_HTO_count <- apply(hto_counts, 2, max)

DefaultAssay(seu) <- "RNA"

# ----------------------------
# Check HTO columns
# ----------------------------

cat("\nAvailable metadata columns:\n")
print(colnames(seu@meta.data))

if (!"HTO_classification.global" %in% colnames(seu@meta.data)) {
  stop("HTO_classification.global not found in Seurat metadata. Did HTODemux fail?")
}

if (!"hash.ID" %in% colnames(seu@meta.data)) {
  stop("hash.ID not found in Seurat metadata. Did HTODemux fail?")
}

seu$hto_class <- seu$HTO_classification.global
seu$hto_donor <- seu$hash.ID

# ----------------------------
# Demux summaries
# ----------------------------

meta <- seu@meta.data %>%
  tibble::rownames_to_column("barcode")

hto_summary <- meta %>%
  count(hto_class, hto_donor, name = "n_cells") %>%
  arrange(desc(n_cells))

vireo_summary <- meta %>%
  count(vireo_class, vireo_donor, name = "n_cells") %>%
  arrange(desc(n_cells))

combined_summary <- meta %>%
  count(hto_class, hto_donor, vireo_class, vireo_donor, name = "n_cells") %>%
  arrange(desc(n_cells))

efficiency_summary <- tibble(
  method = c("HTO", "cellSNP_vireo"),
  total_cells = nrow(meta),
  singlets = c(
    sum(meta$hto_class == "Singlet", na.rm = TRUE),
    sum(meta$vireo_class == "singlet", na.rm = TRUE)
  ),
  doublets = c(
    sum(meta$hto_class == "Doublet", na.rm = TRUE),
    sum(meta$vireo_class == "doublet", na.rm = TRUE)
  ),
  negatives_or_unassigned = c(
    sum(meta$hto_class == "Negative", na.rm = TRUE),
    sum(meta$vireo_class == "unassigned", na.rm = TRUE)
  )
) %>%
  mutate(
    singlet_fraction = singlets / total_cells,
    doublet_fraction = doublets / total_cells,
    negative_or_unassigned_fraction = negatives_or_unassigned / total_cells
  )

agreement_summary <- meta %>%
  mutate(
    hto_is_singlet = hto_class == "Singlet",
    hto_is_doublet = hto_class == "Doublet",
    vireo_is_singlet = vireo_class == "singlet",
    vireo_is_doublet = vireo_class == "doublet"
  ) %>%
  summarise(
    total_cells = n(),
    both_singlet = sum(hto_is_singlet & vireo_is_singlet, na.rm = TRUE),
    HTO_doublet_only = sum(hto_is_doublet & !vireo_is_doublet, na.rm = TRUE),
    Vireo_doublet_only = sum(!hto_is_doublet & vireo_is_doublet, na.rm = TRUE),
    both_doublet = sum(hto_is_doublet & vireo_is_doublet, na.rm = TRUE),
    HTO_singlet_Vireo_unassigned = sum(hto_is_singlet & vireo_class == "unassigned", na.rm = TRUE),
    HTO_negative_Vireo_singlet = sum(hto_class == "Negative" & vireo_is_singlet, na.rm = TRUE)
  ) %>%
  mutate(
    both_singlet_fraction = both_singlet / total_cells,
    both_doublet_fraction = both_doublet / total_cells
  )

donor_cross_tab <- meta %>%
  filter(hto_class == "Singlet", vireo_class == "singlet") %>%
  count(hto_donor, vireo_donor, name = "n_cells") %>%
  arrange(hto_donor, desc(n_cells))

# ----------------------------
# Write CSV/RDS outputs
# ----------------------------

write.csv(hto_summary, file.path(demux_out_dir, "hto_summary.csv"), row.names = FALSE)
write.csv(vireo_summary, file.path(demux_out_dir, "vireo_summary.csv"), row.names = FALSE)
write.csv(combined_summary, file.path(demux_out_dir, "hto_vireo_combined_summary.csv"), row.names = FALSE)
write.csv(efficiency_summary, file.path(demux_out_dir, "demux_efficiency_summary.csv"), row.names = FALSE)
write.csv(agreement_summary, file.path(demux_out_dir, "hto_vireo_agreement_summary.csv"), row.names = FALSE)
write.csv(donor_cross_tab, file.path(demux_out_dir, "hto_vireo_donor_cross_tab.csv"), row.names = FALSE)
write.csv(meta, file.path(demux_out_dir, "cell_metadata_HTO_vireo_combined.csv"), row.names = FALSE)

saveRDS(
  seu,
  file.path(demux_out_dir, paste0("seurat_HTO_vireo_combined_", pq_tag, ".rds"))
)

# ----------------------------
# Plots
# ----------------------------

efficiency_long <- efficiency_summary %>%
  select(method, singlets, doublets, negatives_or_unassigned) %>%
  pivot_longer(
    cols = c(singlets, doublets, negatives_or_unassigned),
    names_to = "category",
    values_to = "n_cells"
  )

pdf(
  file.path(demux_out_dir, paste0("01_demux_efficiency_barplot_", pq_tag, ".pdf")),
  width = 7,
  height = 5
)

print(
  ggplot(efficiency_long, aes(x = method, y = n_cells, fill = category)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Fraction of cells",
      fill = "Classification",
      title = paste0("Demultiplexing efficiency: HTO vs cellSNP-vireo\nHTODemux positive.quantile = ", positive_quantile)
    )
)

dev.off()

doublet_df <- tibble(
  category = c("HTO doublets", "Vireo doublets", "Both doublet"),
  n_cells = c(
    sum(meta$hto_class == "Doublet", na.rm = TRUE),
    sum(meta$vireo_class == "doublet", na.rm = TRUE),
    sum(meta$hto_class == "Doublet" & meta$vireo_class == "doublet", na.rm = TRUE)
  )
)

pdf(
  file.path(demux_out_dir, paste0("02_hto_vireo_doublet_comparison_", pq_tag, ".pdf")),
  width = 6,
  height = 5
)

print(
  ggplot(doublet_df, aes(x = category, y = n_cells)) +
    geom_col() +
    theme_bw() +
    labs(
      x = NULL,
      y = "Number of cells",
      title = paste0("Doublet detection overlap\nHTODemux positive.quantile = ", positive_quantile)
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

pdf(
  file.path(demux_out_dir, paste0("03_hto_vireo_donor_heatmap_", pq_tag, ".pdf")),
  width = 8,
  height = 6
)

print(
  ggplot(donor_cross_tab, aes(x = vireo_donor, y = hto_donor, fill = n_cells)) +
    geom_tile() +
    geom_text(aes(label = n_cells), size = 3) +
    theme_bw() +
    labs(
      x = "Vireo donor",
      y = "HTO donor",
      fill = "Cells",
      title = paste0("HTO versus Vireo donor assignment\nHTODemux positive.quantile = ", positive_quantile)
    )
)

dev.off()

upset_df <- meta %>%
  mutate(
    HTO_singlet = hto_class == "Singlet",
    HTO_doublet = hto_class == "Doublet",
    HTO_negative = hto_class == "Negative",
    Vireo_singlet = vireo_class == "singlet",
    Vireo_doublet = vireo_class == "doublet",
    Vireo_unassigned = vireo_class == "unassigned"
  ) %>%
  select(
    barcode,
    HTO_singlet,
    HTO_doublet,
    HTO_negative,
    Vireo_singlet,
    Vireo_doublet,
    Vireo_unassigned
  )

pdf(
  file.path(demux_out_dir, paste0("04_upset_HTO_vireo_overlap_", pq_tag, ".pdf")),
  width = 10,
  height = 6
)

print(
  upset(
    upset_df,
    intersect = c(
      "HTO_singlet",
      "HTO_doublet",
      "HTO_negative",
      "Vireo_singlet",
      "Vireo_doublet",
      "Vireo_unassigned"
    ),
    name = "Demux call",
    width_ratio = 0.25,
    min_size = 1
  ) +
    ggtitle(paste0(
      "Overlap between HTO and Vireo demultiplexing calls\n",
      "HTODemux positive.quantile = ", positive_quantile
    ))
)

dev.off()

# ----------------------------
# HTO QC plots
# ----------------------------

DefaultAssay(seu) <- "HTO"

hto_counts <- GetAssayData(seu, assay = "HTO", layer = "counts")

seu$HTO_total <- Matrix::colSums(hto_counts)
seu$HTO_max <- apply(hto_counts, 2, max)
seu$HTO_second <- apply(hto_counts, 2, function(x) sort(x, decreasing = TRUE)[2])
seu$HTO_ratio <- seu$HTO_max / (seu$HTO_second + 1)

p1 <- ggplot(seu@meta.data, aes(x = HTO_total)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  theme_bw() +
  ggtitle("HTO total counts per cell")

p2 <- ggplot(seu@meta.data, aes(x = HTO_max)) +
  geom_histogram(bins = 50, fill = "darkorange") +
  scale_x_log10() +
  theme_bw() +
  ggtitle("Max HTO counts per cell")

p3 <- ggplot(seu@meta.data, aes(x = HTO_ratio)) +
  geom_histogram(bins = 50, fill = "purple") +
  theme_bw() +
  ggtitle("HTO max / second ratio")

pdf(
  file.path(demux_out_dir, paste0("05_HTO_QC_", pq_tag, ".pdf")),
  width = 10,
  height = 6
)

print(
  p1 + p2 + p3 +
    plot_annotation(
      title = paste0("HTO QC; HTODemux positive.quantile = ", positive_quantile)
    )
)

dev.off()

pdf(
  file.path(demux_out_dir, paste0("06_HTO_ridge_", pq_tag, ".pdf")),
  width = 12,
  height = 16
)

print(
  RidgePlot(
    seu,
    assay = "HTO",
    features = rownames(seu[["HTO"]]),
    ncol = 1
  ) +
    plot_annotation(
      title = paste0("HTO ridge plot; HTODemux positive.quantile = ", positive_quantile)
    )
)

dev.off()

p_ratio <- ggplot(seu@meta.data, aes(x = HTO_ratio)) +
  geom_histogram(bins = 50, fill = "grey40") +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 5, color = "blue", linetype = "dashed") +
  theme_bw() +
  xlim(c(0, 10)) +
  ggtitle(paste0(
    "HTO ratio thresholds; HTODemux positive.quantile = ",
    positive_quantile
  ))

pdf(
  file.path(demux_out_dir, paste0("06b_HTO_ratio_thresholds_", pq_tag, ".pdf")),
  width = 7,
  height = 5
)

print(p_ratio)

dev.off()

# ----------------------------
# UMAP plots
# ----------------------------

seu$HTO_singlet <- seu$HTO_classification.global == "Singlet"
seu$Vireo_singlet <- seu$vireo_doublet_prob < 0.5

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

dims_use <- 1:min(30, ncol(Embeddings(seu, "pca")))
seu <- RunUMAP(seu, dims = dims_use)

seu_clean <- subset(seu, subset = HTO_singlet & Vireo_singlet)
seu_clean_vireo <- subset(seu_clean, subset = vireo_donor != "unassigned")

hto_levels <- levels(factor(seu_clean$hash.ID))
vireo_levels <- levels(factor(seu_clean_vireo$vireo_donor))

hto_cols <- setNames(
  brewer.pal(max(3, length(hto_levels)), "Set1")[seq_along(hto_levels)],
  hto_levels
)

vireo_cols <- setNames(
  brewer.pal(max(3, length(vireo_levels)), "Set1")[seq_along(vireo_levels)],
  vireo_levels
)

p1 <- DimPlot(
  seu_clean,
  group.by = "hash.ID",
  reduction = "umap",
  cols = hto_cols
) +
  ggtitle("HTO assignment") +
  theme_bw()

p2 <- DimPlot(
  seu_clean_vireo,
  group.by = "vireo_donor",
  reduction = "umap",
  cols = vireo_cols
) +
  ggtitle("Vireo assignment (singlets only)") +
  theme_bw()

pdf(
  file.path(demux_out_dir, paste0("07_UMAP_HTO_cellSNP_", pq_tag, ".pdf")),
  width = 12,
  height = 8
)

print(
  p1 + p2 +
    plot_annotation(
      title = paste0("HTO and Vireo assignments; HTODemux positive.quantile = ", positive_quantile)
    )
)

dev.off()

saveRDS(
  seu,
  file.path(demux_out_dir, paste0("seurat_HTO_vireo_combined_", pq_tag, ".rds"))
)

cat("\nDone.\n")
cat("Output written to:", demux_out_dir, "\n")
cat("HTODemux positive.quantile:", positive_quantile, "\n\n")

cat("Demux efficiency:\n")
print(efficiency_summary)

cat("\nAgreement summary:\n")
print(agreement_summary)

cat("\nHTO summary:\n")
print(hto_summary)

cat("\nVireo summary:\n")
print(vireo_summary)
