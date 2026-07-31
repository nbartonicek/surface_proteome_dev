suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(scales)
})

# ============================================================
# Compare CITE_DSB vs ADT_CLR in one sample
# ============================================================

# ----------------------------
# Inputs
# ----------------------------
run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"
sample_short <- "LK1"

projection_path <- "../annotation/"
out_dir <- file.path("../results/seurat_annotated", run)
projection_dir <- file.path(out_dir, "projectionFigures/")

cite_dsb_dir <- file.path(
  "../results/cite_qc_then_dsb",
  run,
  sample_name
)

annotation_dir <- file.path("../results/seurat_annotated", run)

seurat_input <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

seu <- readRDS(seurat_input)

sample_to_use <- "HBDN206-MNpCT"

rna_assay <- "RNA"
clr_assay <- "ADT_CLR"
dsb_assay <- "CITE_DSB"

sample_col <- "sample_name"
broad_col  <- "predicted_CellType_Broad"
fine_col   <- "predicted_CellType"


reduction_to_use <- "umap_projected"
if (!reduction_to_use %in% Reductions(seu)) {
  reduction_to_use <- "umap"
}

compare_dir <- file.path(
  projection_dir,
  paste0("DSB_vs_CLR_", sample_to_use)
)
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Checks
# ----------------------------

stopifnot(rna_assay %in% Assays(seu))
stopifnot(clr_assay %in% Assays(seu))
stopifnot(dsb_assay %in% Assays(seu))
stopifnot(sample_col %in% colnames(seu@meta.data))
stopifnot(broad_col %in% colnames(seu@meta.data))
stopifnot(reduction_to_use %in% Reductions(seu))

# ----------------------------
# Subset one sample
# ----------------------------

seu_sub <- subset(
  seu,
  cells = rownames(seu@meta.data)[seu@meta.data[[sample_col]] == sample_to_use]
)

cat("Using sample:", sample_to_use, "\n")
cat("Cells:", ncol(seu_sub), "\n")
cat("Assays:", paste(Assays(seu_sub), collapse = ", "), "\n")
cat("Reduction:", reduction_to_use, "\n\n")

# ----------------------------
# Helper functions
# ----------------------------

get_layer <- function(object, assay, layer = "data") {
  GetAssayData(object, assay = assay, layer = layer)
}

clean_marker_name <- function(x) {
  x %>%
    str_replace_all("-[ACGT]{10,}$", "") %>%
    str_replace_all("_TotalSeq.*$", "") %>%
    str_replace_all("TotalSeq.*$", "") %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "") %>%
    str_replace_all("_", "") %>%
    toupper()
}

is_background_marker <- function(x) {
  grepl(
    "IgG|isotype|Biotin|TotalSeq|Hash|HTO|control|Ctrl",
    x,
    ignore.case = TRUE
  )
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

# ----------------------------
# Pull matrices
# ----------------------------

rna_mat <- get_layer(seu_sub, rna_assay, "data")
clr_mat <- get_layer(seu_sub, clr_assay, "data")
dsb_mat <- get_layer(seu_sub, dsb_assay, "data")

shared_proteins <- intersect(rownames(clr_mat), rownames(dsb_mat))
shared_proteins <- shared_proteins[!is_background_marker(shared_proteins)]

clr_mat <- clr_mat[shared_proteins, , drop = FALSE]
dsb_mat <- dsb_mat[shared_proteins, , drop = FALSE]

meta <- seu_sub@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    broad = .data[[broad_col]],
    fine = if (fine_col %in% colnames(.)) .data[[fine_col]] else broad
  )

# ============================================================
# 1. Global DSB - CLR marker shift
# ============================================================

marker_shift <- tibble(
  marker = shared_proteins,
  mean_CLR = Matrix::rowMeans(clr_mat),
  mean_DSB = Matrix::rowMeans(dsb_mat),
  delta_DSB_minus_CLR = mean_DSB - mean_CLR
) %>%
  arrange(desc(delta_DSB_minus_CLR))

write.csv(
  marker_shift,
  file.path(compare_dir, "01_marker_mean_shift_DSB_minus_CLR.csv"),
  row.names = FALSE
)

p_shift <- marker_shift %>%
  slice_max(abs(delta_DSB_minus_CLR), n = 30) %>%
  mutate(marker = fct_reorder(marker, delta_DSB_minus_CLR)) %>%
  ggplot(aes(x = delta_DSB_minus_CLR, y = marker)) +
  geom_col(fill = "grey35") +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("Markers most changed by DSB: ", sample_to_use),
    x = "Mean DSB - mean CLR",
    y = NULL
  )

ggsave(
  file.path(compare_dir, "01_marker_mean_shift_DSB_minus_CLR.pdf"),
  p_shift,
  width = 10,
  height = 8
)

# ============================================================
# 2. RNA-protein correlation
# ============================================================

# Manual ADT-to-RNA gene map.
# Add/remove based on your panel names.
adt_gene_map <- tribble(
  ~marker_pattern, ~gene,
  "CD3", "CD3D",
  "CD4", "CD4",
  "CD8", "CD8A",
  "CD11A", "ITGAL",
  "CD11B", "ITGAM",
  "CD11C", "ITGAX",
  "CD14", "CD14",
  "CD16", "FCGR3A",
  "CD18", "ITGB2",
  "CD19", "CD19",
  "CD20", "MS4A1",
  "CD22", "CD22",
  "CD23", "FCER2",
  "CD31", "PECAM1",
  "CD33", "CD33",
  "CD34", "CD34",
  "CD38", "CD38",
  "CD41", "ITGA2B",
  "CD43", "SPN",
  "CD44", "CD44",
  "CD45", "PTPRC",
  "CD45RA", "PTPRC",
  "CD45RO", "PTPRC",
  "CD47", "CD47",
  "CD49D", "ITGA4",
  "CD54", "ICAM1",
  "CD56", "NCAM1",
  "CD61", "ITGB3",
  "CD62L", "SELL",
  "CD64", "FCGR1A",
  "CD71", "TFRC",
  "CD83", "CD83",
  "CD86", "CD86",
  "CD99", "CD99",
  "CD117", "KIT",
  "CD123", "IL3RA",
  "CD133", "PROM1",
  "CD138", "SDC1",
  "SYNDECAN1", "SDC1",
  "CD163", "CD163",
  "CD244", "CD244",
  "CD274", "CD274",
  "BCMA", "TNFRSF17",
  "FLT3", "FLT3",
  "GPR56", "ADGRG1",
  "THROMBOMODULIN", "THBD",
  "CXCR1", "CXCR1",
  "CXCR2", "CXCR2",
  "CCR3", "CCR3",
  "CCR7", "CCR7",
  "OX40", "TNFRSF4",
  "IL2RB", "IL2RB"
)

marker_key <- tibble(
  marker = shared_proteins,
  marker_clean = clean_marker_name(shared_proteins)
)

corr_input <- marker_key %>%
  crossing(adt_gene_map) %>%
  dplyr::filter(str_detect(marker_clean, marker_pattern)) %>%
  distinct(marker, gene) %>%
  dplyr::filter(gene %in% rownames(rna_mat))

rna_protein_cor <- corr_input %>%
  rowwise() %>%
  mutate(
    cor_CLR = safe_cor(
      as.numeric(rna_mat[gene, colnames(clr_mat)]),
      as.numeric(clr_mat[marker, ])
    ),
    cor_DSB = safe_cor(
      as.numeric(rna_mat[gene, colnames(dsb_mat)]),
      as.numeric(dsb_mat[marker, ])
    ),
    delta_DSB_minus_CLR = cor_DSB - cor_CLR
  ) %>%
  ungroup() %>%
  arrange(desc(delta_DSB_minus_CLR))

write.csv(
  rna_protein_cor,
  file.path(compare_dir, "02_RNA_protein_correlation_DSB_vs_CLR.csv"),
  row.names = FALSE
)

p_cor <- rna_protein_cor %>%
  dplyr::filter(is.finite(cor_CLR), is.finite(cor_DSB)) %>%
  ggplot(aes(x = cor_CLR, y = cor_DSB, label = marker)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 30) +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("RNA-protein correlation: DSB vs CLR: ", sample_to_use),
    subtitle = "Points above diagonal support better RNA-protein agreement after DSB",
    x = "Spearman RNA vs CLR protein",
    y = "Spearman RNA vs DSB protein"
  )

ggsave(
  file.path(compare_dir, "02_RNA_protein_correlation_DSB_vs_CLR.pdf"),
  p_cor,
  width = 8,
  height = 7
)

# ============================================================
# 3. Cell-type specificity score
#    For each marker:
#    max broad mean - median broad mean
# ============================================================

make_long_protein <- function(mat, norm_name) {
  as.data.frame(as.matrix(mat)) %>%
    rownames_to_column("marker") %>%
    pivot_longer(
      cols = -marker,
      names_to = "cell",
      values_to = "value"
    ) %>%
    left_join(meta %>% select(cell, broad), by = "cell") %>%
    filter(!is.na(broad)) %>%
    mutate(norm = norm_name)
}

protein_long <- bind_rows(
  make_long_protein(clr_mat, "CLR"),
  make_long_protein(dsb_mat, "DSB")
)

broad_summary <- protein_long %>%
  group_by(norm, broad, marker) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    pct_positive_0 = mean(value > 0, na.rm = TRUE) * 100,
    pct_positive_3 = mean(value > 3, na.rm = TRUE) * 100,
    .groups = "drop"
  )

specificity <- broad_summary %>%
  group_by(norm, marker) %>%
  summarise(
    max_broad = broad[which.max(mean_value)],
    max_mean = max(mean_value, na.rm = TRUE),
    median_mean = median(mean_value, na.rm = TRUE),
    specificity_score = max_mean - median_mean,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = norm,
    values_from = c(max_broad, max_mean, median_mean, specificity_score)
  ) %>%
  mutate(
    delta_specificity_DSB_minus_CLR =
      specificity_score_DSB - specificity_score_CLR
  ) %>%
  arrange(desc(delta_specificity_DSB_minus_CLR))

write.csv(
  specificity,
  file.path(compare_dir, "03_marker_celltype_specificity_DSB_vs_CLR.csv"),
  row.names = FALSE
)

p_spec <- specificity %>%
  filter(is.finite(specificity_score_CLR), is.finite(specificity_score_DSB)) %>%
  ggplot(aes(
    x = specificity_score_CLR,
    y = specificity_score_DSB,
    label = marker
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 25) +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("Cell-type specificity: DSB vs CLR: ", sample_to_use),
    subtitle = "Points above diagonal are more cell-type-specific after DSB",
    x = "CLR specificity score",
    y = "DSB specificity score"
  )

ggsave(
  file.path(compare_dir, "03_celltype_specificity_DSB_vs_CLR.pdf"),
  p_spec,
  width = 8,
  height = 7
)

# ============================================================
# 4. Dotplots: same markers, CLR vs DSB
# ============================================================

markers_for_dotplot <- specificity %>%
  slice_max(abs(delta_specificity_DSB_minus_CLR), n = 25) %>%
  pull(marker) %>%
  unique()

markers_for_dotplot <- markers_for_dotplot[markers_for_dotplot %in% shared_proteins]

DefaultAssay(seu_sub) <- clr_assay
p_dot_clr <- DotPlot(
  seu_sub,
  features = markers_for_dotplot,
  group.by = broad_col,
  assay = clr_assay
) +
  RotatedAxis() +
  theme_bw(base_size = 11) +
  labs(title = "CLR")

DefaultAssay(seu_sub) <- dsb_assay
p_dot_dsb <- DotPlot(
  seu_sub,
  features = markers_for_dotplot,
  group.by = broad_col,
  assay = dsb_assay
) +
  RotatedAxis() +
  theme_bw(base_size = 11) +
  labs(title = "DSB")

ggsave(
  file.path(compare_dir, "04_dotplot_CLR_vs_DSB_top_changed_specificity.pdf"),
  p_dot_clr / p_dot_dsb,
  width = 14,
  height = 10
)

# ============================================================
# 5. FeaturePlots on projected UMAP
# ============================================================

markers_for_feature <- c(
  "CD34", "CD33", "CD117", "CD14", "CD163",
  "CD4", "CD8", "CD19", "CD38", "BCMA", "Syndecan-1"
)

markers_for_feature <- markers_for_feature[
  markers_for_feature %in% shared_proteins
]

if (length(markers_for_feature) > 0) {
  
  pdf(
    file.path(compare_dir, "05_featureplots_CLR_vs_DSB_selected_markers.pdf"),
    width = 14,
    height = 8
  )
  
  for (mk in markers_for_feature) {
    
    DefaultAssay(seu_sub) <- clr_assay
    p1 <- FeaturePlot(
      seu_sub,
      features = mk,
      reduction = reduction_to_use,
      order = TRUE
    ) +
      labs(title = paste0(mk, " CLR"))
    
    DefaultAssay(seu_sub) <- dsb_assay
    p2 <- FeaturePlot(
      seu_sub,
      features = mk,
      reduction = reduction_to_use,
      order = TRUE
    ) +
      labs(title = paste0(mk, " DSB"))
    
    print(p1 + p2)
  }
  
  dev.off()
}

# ============================================================
# 6. Violin plots by broad cell type
# ============================================================

if (length(markers_for_feature) > 0) {
  
  pdf(
    file.path(compare_dir, "06_violin_CLR_vs_DSB_selected_markers_by_broad.pdf"),
    width = 14,
    height = 8
  )
  
  for (mk in markers_for_feature) {
    
    df_mk <- protein_long %>%
      filter(marker == mk) %>%
      mutate(
        broad = fct_reorder(broad, value, .fun = median, .desc = TRUE)
      )
    
    p <- ggplot(df_mk, aes(x = broad, y = value, fill = norm)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.5) +
      facet_wrap(~norm, ncol = 1, scales = "free_y") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      labs(
        title = paste0(mk, ": CLR vs DSB by broad annotation"),
        x = NULL,
        y = "Protein value"
      )
    
    print(p)
  }
  
  dev.off()
}

# ============================================================
# 7. Expected marker panel summary
# ============================================================

expected_markers <- tribble(
  ~expected_broad, ~marker,
  "Stem / progenitor", "CD34",
  "Stem / progenitor", "CD117",
  "Stem / progenitor", "CD133",
  "Myeloid / DC", "CD33",
  "Myeloid / DC", "CD14",
  "Myeloid / DC", "CD163",
  "Myeloid / DC", "CD11b",
  "Lymphoid", "CD4",
  "Lymphoid", "CD8",
  "Lymphoid", "CD19",
  "Lymphoid", "CD20",
  "Lymphoid", "BCMA",
  "Lymphoid", "Syndecan-1",
  "Erythroid", "CD71"
) %>%
  filter(marker %in% shared_proteins)

expected_eval <- broad_summary %>%
  inner_join(expected_markers, by = "marker") %>%
  group_by(norm, marker, expected_broad) %>%
  summarise(
    expected_mean = mean(mean_value[broad == expected_broad], na.rm = TRUE),
    other_mean = median(mean_value[broad != expected_broad], na.rm = TRUE),
    expected_specificity = expected_mean - other_mean,
    expected_pct_positive_0 =
      mean(pct_positive_0[broad == expected_broad], na.rm = TRUE),
    other_pct_positive_0 =
      median(pct_positive_0[broad != expected_broad], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = norm,
    values_from = c(
      expected_mean,
      other_mean,
      expected_specificity,
      expected_pct_positive_0,
      other_pct_positive_0
    )
  ) %>%
  mutate(
    delta_expected_specificity_DSB_minus_CLR =
      expected_specificity_DSB - expected_specificity_CLR
  ) %>%
  arrange(desc(delta_expected_specificity_DSB_minus_CLR))

write.csv(
  expected_eval,
  file.path(compare_dir, "07_expected_marker_specificity_DSB_vs_CLR.csv"),
  row.names = FALSE
)

p_expected <- expected_eval %>%
  ggplot(aes(
    x = expected_specificity_CLR,
    y = expected_specificity_DSB,
    label = marker
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("Expected marker specificity: DSB vs CLR: ", sample_to_use),
    subtitle = "Above diagonal means expected marker is more specific after DSB",
    x = "CLR expected specificity",
    y = "DSB expected specificity"
  )

ggsave(
  file.path(compare_dir, "07_expected_marker_specificity_DSB_vs_CLR.pdf"),
  p_expected,
  width = 8,
  height = 7
)

# ============================================================
# 8. Background/isotype suppression check
# ============================================================

bg_markers <- rownames(get_layer(seu_sub, clr_assay, "data"))[
  is_background_marker(rownames(get_layer(seu_sub, clr_assay, "data")))
]

bg_markers <- intersect(
  bg_markers,
  rownames(get_layer(seu_sub, dsb_assay, "data"))
)

if (length(bg_markers) > 0) {
  
  clr_bg <- get_layer(seu_sub, clr_assay, "data")[bg_markers, , drop = FALSE]
  dsb_bg <- get_layer(seu_sub, dsb_assay, "data")[bg_markers, , drop = FALSE]
  
  bg_summary <- tibble(
    marker = bg_markers,
    mean_CLR = Matrix::rowMeans(clr_bg),
    mean_DSB = Matrix::rowMeans(dsb_bg),
    delta_DSB_minus_CLR = mean_DSB - mean_CLR
  ) %>%
    arrange(delta_DSB_minus_CLR)
  
  write.csv(
    bg_summary,
    file.path(compare_dir, "08_background_marker_suppression_DSB_vs_CLR.csv"),
    row.names = FALSE
  )
  
  p_bg <- bg_summary %>%
    mutate(marker = fct_reorder(marker, delta_DSB_minus_CLR)) %>%
    ggplot(aes(x = delta_DSB_minus_CLR, y = marker)) +
    geom_col(fill = "grey35") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Background/isotype suppression by DSB: ", sample_to_use),
      x = "Mean DSB - mean CLR",
      y = NULL
    )
  
  ggsave(
    file.path(compare_dir, "08_background_marker_suppression_DSB_vs_CLR.pdf"),
    p_bg,
    width = 8,
    height = max(4, length(bg_markers) * 0.25)
  )
}

# ============================================================
# 9. Final compact summary table
# ============================================================

summary_table <- tibble(
  metric = c(
    "Median RNA-protein correlation CLR",
    "Median RNA-protein correlation DSB",
    "Median delta RNA-protein correlation DSB-CLR",
    "Median marker specificity CLR",
    "Median marker specificity DSB",
    "Median delta marker specificity DSB-CLR",
    "Median expected-marker specificity CLR",
    "Median expected-marker specificity DSB",
    "Median delta expected-marker specificity DSB-CLR"
  ),
  value = c(
    median(rna_protein_cor$cor_CLR, na.rm = TRUE),
    median(rna_protein_cor$cor_DSB, na.rm = TRUE),
    median(rna_protein_cor$delta_DSB_minus_CLR, na.rm = TRUE),
    median(specificity$specificity_score_CLR, na.rm = TRUE),
    median(specificity$specificity_score_DSB, na.rm = TRUE),
    median(specificity$delta_specificity_DSB_minus_CLR, na.rm = TRUE),
    median(expected_eval$expected_specificity_CLR, na.rm = TRUE),
    median(expected_eval$expected_specificity_DSB, na.rm = TRUE),
    median(expected_eval$delta_expected_specificity_DSB_minus_CLR, na.rm = TRUE)
  )
)

write.csv(
  summary_table,
  file.path(compare_dir, "09_summary_metrics_DSB_vs_CLR.csv"),
  row.names = FALSE
)

print(summary_table)

cat("\nDone. Outputs written to:\n", compare_dir, "\n")