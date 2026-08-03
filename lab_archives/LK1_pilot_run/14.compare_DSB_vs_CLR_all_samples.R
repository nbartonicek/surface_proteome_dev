# ------------------------------------------------------------------
# LK1 pilot run - step 14 of 19
#
# DSB-vs-CLR across all samples: background suppression, expected-marker specificity, marker leakage. Writes data_integration/ 01-06.
#
# Frozen for the lab archive 2026-07-31 from scripts/backup/10c.compare_DSB_CLR.R (mtime 2026-05-20).
# md5 of the original: f87369da7fdb179ad8aa0eb70d5f1790
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

sample_to_use <- "HBDN206-MNpCT"

seu_sub <- subset(
  seu,
  subset = sampleID == sample_to_use
)

cat("Cells:", ncol(seu_sub), "\n")
celltype_col <- "predicted_CellType_Broad"

celltype_vec <- as.character(
  seu_sub@meta.data[[celltype_col]]
)

table(celltype_vec)

#################### cleaning
#celltype_vec <- seu_sub@meta.data[[broad_col]]

seu_sub$lineage_broad <- dplyr::case_when(
  
  celltype_vec %in% c(
    "HSC MPP", "LMPP", "MEP",
    "Early GMP", "Late GMP",
    "Cycling Progenitor",
    "EoBasoMast Precursor",
    "Megakaryocyte Precursor"
  ) ~ "Stem / progenitor",
  
  celltype_vec %in% c(
    "Monocyte", "Pro-Monocyte", "pDC"
  ) ~ "Myeloid / DC",
  
  celltype_vec %in% c(
    "B", "Pro-B", "Early Lymphoid",
    "Naive T", "CD4 Memory T",
    "CD8 Memory T", "NK",
    "Plasma Cell"
  ) ~ "Lymphoid",
  
  celltype_vec %in% c(
    "Early Erythroid",
    "Late Erythroid"
  ) ~ "Erythroid",
  
  TRUE ~ "Other"
)

table(seu_sub$lineage_broad)

get_mat <- function(obj, assay) {
  GetAssayData(obj, assay = assay, layer = "data")
}

is_bg <- function(x) {
  grepl("IgG|IgM|isotype|control|Fc|Hash|HTO|TotalSeq",
        x, ignore.case = TRUE)
}

marker_exists <- function(markers, mat) {
  markers[markers %in% rownames(mat)]
}

safe_ratio <- function(a, b) {
  ifelse(is.na(b) | b == 0, NA_real_, a / b)
}

# ----------------------------
# Matrices
# ----------------------------
clr_assay <- "ADT_CLR"
dsb_assay <- "CITE_DSB"
broad_col <- "lineage_broad"

clr <- get_mat(seu_sub, clr_assay)
dsb <- get_mat(seu_sub, dsb_assay)

common_markers <- intersect(rownames(clr), rownames(dsb))

clr <- clr[common_markers, , drop = FALSE]
dsb <- dsb[common_markers, , drop = FALSE]

meta <- seu_sub@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, broad = all_of(broad_col)) %>%
  dplyr::filter(!is.na(broad))

# ============================================================
# 1. Background suppression
# ============================================================

bg_markers <- common_markers[is_bg(common_markers)]

bg_summary <- tibble(
  marker = bg_markers,
  mean_CLR = Matrix::rowMeans(clr[bg_markers, , drop = FALSE]),
  mean_DSB = Matrix::rowMeans(dsb[bg_markers, , drop = FALSE]),
  delta_DSB_minus_CLR = mean_DSB - mean_CLR
) %>%
  arrange(delta_DSB_minus_CLR)

write.csv(
  bg_summary,
  file.path(out_dir, "01_background_suppression.csv"),
  row.names = FALSE
)

p_bg <- bg_summary %>%
  mutate(marker = fct_reorder(marker, delta_DSB_minus_CLR)) %>%
  ggplot(aes(x = delta_DSB_minus_CLR, y = marker)) +
  geom_col(fill = "grey35") +
  theme_bw(base_size = 13) +
  labs(
    title = "Background / isotype suppression",
    subtitle = "Negative values mean DSB reduced background relative to CLR",
    x = "Mean DSB - mean CLR",
    y = NULL
  )

ggsave(
  file.path(out_dir, "01_background_suppression.pdf"),
  p_bg,
  width = 8,
  height = max(4, 0.25 * length(bg_markers))
)

# ============================================================
# 2. Expected marker specificity
# ============================================================

expected_markers <- tribble(
  ~expected_broad,        ~marker,
  "Myeloid / DC",         "CD163",
  "Myeloid / DC",         "CD14",
  "Myeloid / DC",         "CD33",
  "Myeloid / DC",         "CD11b",
  "Stem / progenitor",    "CD34",
  "Stem / progenitor",    "CD117",
  "Stem / progenitor",    "CD133",
  "Lymphoid",             "CD4",
  "Lymphoid",             "CD8",
  "Lymphoid",             "CD19",
  "Lymphoid",             "CD20",
  "Erythroid",            "CD71"
) %>%
  dplyr::filter(marker %in% common_markers)

make_long <- function(mat, norm_name) {
  as.data.frame(as.matrix(mat[expected_markers$marker, , drop = FALSE])) %>%
    rownames_to_column("marker") %>%
    pivot_longer(
      cols = -marker,
      names_to = "cell",
      values_to = "value"
    ) %>%
    left_join(meta, by = "cell") %>%
    left_join(expected_markers, by = "marker") %>%
    mutate(
      norm = norm_name,
      is_expected = broad == expected_broad
    )
}

long_df <- bind_rows(
  make_long(clr, "CLR"),
  make_long(dsb, "DSB")
)

specificity <- long_df %>%
  group_by(norm, marker, expected_broad) %>%
  summarise(
    expected_mean = mean(value[is_expected], na.rm = TRUE),
    other_median = median(value[!is_expected], na.rm = TRUE),
    specificity = expected_mean - other_median,
    expected_pct_positive = mean(value[is_expected] > 0, na.rm = TRUE) * 100,
    other_pct_positive = mean(value[!is_expected] > 0, na.rm = TRUE) * 100,
    pct_specificity = expected_pct_positive - other_pct_positive,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = norm,
    values_from = c(
      expected_mean,
      other_median,
      specificity,
      expected_pct_positive,
      other_pct_positive,
      pct_specificity
    )
  ) %>%
  mutate(
    delta_specificity_DSB_minus_CLR = specificity_DSB - specificity_CLR,
    delta_pct_specificity_DSB_minus_CLR = pct_specificity_DSB - pct_specificity_CLR
  ) %>%
  arrange(desc(delta_specificity_DSB_minus_CLR))

write.csv(
  specificity,
  file.path(out_dir, "02_expected_marker_specificity.csv"),
  row.names = FALSE
)

p_spec <- specificity %>%
  ggplot(aes(x = specificity_CLR, y = specificity_DSB, label = marker)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Expected marker specificity: DSB vs CLR",
    subtitle = "Above diagonal means marker is more cell-type-specific after DSB",
    x = "CLR specificity",
    y = "DSB specificity"
  )

ggsave(
  file.path(out_dir, "02_expected_marker_specificity.pdf"),
  p_spec,
  width = 8,
  height = 7
)

p_pct_spec <- specificity %>%
  ggplot(aes(
    x = pct_specificity_CLR,
    y = pct_specificity_DSB,
    label = marker
  )) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    max.overlaps = 50
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Expected marker percent-specificity: DSB vs CLR",
    subtitle = paste(
      "Percent specificity =",
      "% positive expected cells - % positive other cells",
      "\nAbove diagonal means DSB restricted positivity",
      "more strongly to expected populations"
    ),
    x = "CLR percent specificity",
    y = "DSB percent specificity"
  )

ggsave(
  file.path(out_dir, "02_expected_marker_percent_specificity.pdf"),
  p_pct_spec,
  width = 8,
  height = 7
)

# ============================================================
# 3. Marker leakage / diffusion
#    Signal outside expected population should decrease.
# ============================================================

leakage <- long_df %>%
  group_by(norm, marker, expected_broad) %>%
  summarise(
    inside_mean = mean(value[is_expected], na.rm = TRUE),
    outside_mean = mean(value[!is_expected], na.rm = TRUE),
    outside_positive = mean(value[!is_expected] > 0, na.rm = TRUE) * 100,
    leakage_ratio = safe_ratio(outside_mean, inside_mean),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = norm,
    values_from = c(
      inside_mean,
      outside_mean,
      outside_positive,
      leakage_ratio
    )
  ) %>%
  mutate(
    delta_outside_mean_DSB_minus_CLR = outside_mean_DSB - outside_mean_CLR,
    delta_outside_positive_DSB_minus_CLR = outside_positive_DSB - outside_positive_CLR,
    delta_leakage_ratio_DSB_minus_CLR = leakage_ratio_DSB - leakage_ratio_CLR
  ) %>%
  arrange(delta_leakage_ratio_DSB_minus_CLR)

write.csv(
  leakage,
  file.path(out_dir, "03_marker_leakage_diffusion.csv"),
  row.names = FALSE
)

p_leak <- leakage %>%
  ggplot(aes(x = leakage_ratio_CLR, y = leakage_ratio_DSB, label = marker)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Marker leakage outside expected population",
    subtitle = "Below diagonal means DSB reduced off-target/background spread",
    x = "CLR outside / inside signal",
    y = "DSB outside / inside signal"
  )

ggsave(
  file.path(out_dir, "03_marker_leakage_diffusion.pdf"),
  p_leak,
  width = 8,
  height = 7
)

# ============================================================
# 4. Per-marker violin plots: expected vs other cells
# ============================================================

p_violin <- long_df %>%
  mutate(
    population = if_else(is_expected, "Expected cell type", "Other cells"),
    population = factor(population, levels = c("Other cells", "Expected cell type"))
  ) %>%
  ggplot(aes(x = population, y = value, fill = population)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.6) +
  facet_grid(marker ~ norm, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Expected marker signal in expected vs other cells",
    x = NULL,
    y = "Protein value"
  )

ggsave(
  file.path(out_dir, "04_expected_vs_other_violin.pdf"),
  p_violin,
  width = 10,
  height = max(7, 2 * nrow(expected_markers))
)

# ============================================================
# 5. FeaturePlot comparisons
# ============================================================

reduction_use <- if ("umap_projected" %in% Reductions(seu_sub)) {
  "umap_projected"
} else {
  "umap"
}

feature_markers <- expected_markers$marker %>%
  unique() %>%
  marker_exists(clr)

pdf(
  file.path(out_dir, "05_featureplots_CLR_vs_DSB_expected_markers.pdf"),
  width = 13,
  height = 6
)

for (mk in feature_markers) {
  DefaultAssay(seu_sub) <- clr_assay
  p1 <- FeaturePlot(
    seu_sub,
    features = mk,
    reduction = reduction_use,
    order = TRUE
  ) +
    ggtitle(paste0(mk, " CLR"))
  
  DefaultAssay(seu_sub) <- dsb_assay
  p2 <- FeaturePlot(
    seu_sub,
    features = mk,
    reduction = reduction_use,
    order = TRUE
  ) +
    ggtitle(paste0(mk, " DSB"))
  
  print(p1 + p2)
}

dev.off()

# ============================================================
# 6. Compact final summary
# ============================================================

summary_metrics <- tibble(
  metric = c(
    "Median background marker mean CLR",
    "Median background marker mean DSB",
    "Median background delta DSB-CLR",
    "Median expected-marker specificity CLR",
    "Median expected-marker specificity DSB",
    "Median specificity delta DSB-CLR",
    "Median outside/inside leakage ratio CLR",
    "Median outside/inside leakage ratio DSB",
    "Median leakage-ratio delta DSB-CLR"
  ),
  value = c(
    median(bg_summary$mean_CLR, na.rm = TRUE),
    median(bg_summary$mean_DSB, na.rm = TRUE),
    median(bg_summary$delta_DSB_minus_CLR, na.rm = TRUE),
    median(specificity$specificity_CLR, na.rm = TRUE),
    median(specificity$specificity_DSB, na.rm = TRUE),
    median(specificity$delta_specificity_DSB_minus_CLR, na.rm = TRUE),
    median(leakage$leakage_ratio_CLR, na.rm = TRUE),
    median(leakage$leakage_ratio_DSB, na.rm = TRUE),
    median(leakage$delta_leakage_ratio_DSB_minus_CLR, na.rm = TRUE)
  )
)

write.csv(
  summary_metrics,
  file.path(out_dir, "06_summary_metrics.csv"),
  row.names = FALSE
)

print(summary_metrics)

cat("\nDone. Files written to:\n", out_dir, "\n")

