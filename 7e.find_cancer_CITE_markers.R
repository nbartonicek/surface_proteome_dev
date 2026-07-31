# ============================================================
# CITE_DSB HSC compartment analysis with Numbat calls
# - HSC CITE marker enrichment
# - HSC cancer vs normal-by-Numbat DE
# - Validation in true normal-sample HSCs
# - BAFF-R RNA/protein correlation
# - Manual selected-marker dotplot
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
  library(scales)
})

# ============================================================
# Inputs
# ============================================================

run <- "260528_VH01624_464_222K7VKNX"
sample_name <- "LK2-GEX"
sample_short <- "LK2"

out_dir <- file.path("../results/seurat_annotated", run, "numbat")

seurat_input <- file.path(
  out_dir,
  paste0(sample_short, "_projected_CITE_DSB_Numbat_integrated.rds")
)

seu <- readRDS(seurat_input)

analysis_dir <- file.path(out_dir, "CITE_HSC_Numbat_analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Columns / assays
# ============================================================

dsb_assay <- "CITE_DSB"
rna_assay <- "RNA"

sample_col <- "sample_name"
broad_col  <- "predicted_CellType_Broad"
fine_col   <- "predicted_CellType"
numbat_call_col <- "numbat_call"

stopifnot(dsb_assay %in% Assays(seu))
stopifnot(rna_assay %in% Assays(seu))
stopifnot(sample_col %in% colnames(seu@meta.data))
stopifnot(broad_col %in% colnames(seu@meta.data))
stopifnot(numbat_call_col %in% colnames(seu@meta.data))

reduction_to_use <- "umap_projected"
if (!reduction_to_use %in% Reductions(seu)) {
  reduction_to_use <- "umap"
}

# ============================================================
# Helper functions
# ============================================================

get_layer <- function(object, assay, layer = "data") {
  GetAssayData(object, assay = assay, layer = layer)
}

is_background_marker <- function(x) {
  grepl(
    "IgG|isotype|Biotin|TotalSeq|Hash|HTO|control|Ctrl",
    x,
    ignore.case = TRUE
  )
}

safe_wilcox <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- droplevels(factor(group[ok]))
  
  if (length(levels(group)) != 2) return(NA_real_)
  if (sum(group == levels(group)[1]) < 5) return(NA_real_)
  if (sum(group == levels(group)[2]) < 5) return(NA_real_)
  
  suppressWarnings(wilcox.test(x ~ group)$p.value)
}

# ============================================================
# Define HSC compartment
# ============================================================

hsc_patterns <- c(
  "HSC",
  "MPP",
  "HSC MPP",
  "stem",
  "progenitor"
)

seu$hsc_compartment <- ifelse(
  grepl(
    paste(hsc_patterns, collapse = "|"),
    seu@meta.data[[broad_col]],
    ignore.case = TRUE
  ),
  "HSC_MPP",
  "Other"
)

seu$hsc_compartment <- factor(
  seu$hsc_compartment,
  levels = c("Other", "HSC_MPP")
)

message("HSC compartment counts:")
print(table(seu$hsc_compartment, useNA = "ifany"))

# ============================================================
# Define Numbat cancer / normal
# ============================================================

raw_call <- as.character(seu@meta.data[[numbat_call_col]])

seu$numbat_status_simplified <- case_when(
  grepl("normal|diploid|neutral|non.?malig|healthy", raw_call, ignore.case = TRUE) ~ "Normal_by_Numbat",
  grepl("cancer|malig|tumou|tumor|aneuploid|aberrant|clone|CNV", raw_call, ignore.case = TRUE) ~ "Cancer_by_Numbat",
  TRUE ~ NA_character_
)

# fallback for numeric/clone-style calls
seu$numbat_status_simplified <- ifelse(
  is.na(seu$numbat_status_simplified) &
    grepl("^0$|normal", raw_call, ignore.case = TRUE),
  "Normal_by_Numbat",
  seu$numbat_status_simplified
)

seu$numbat_status_simplified <- factor(
  seu$numbat_status_simplified,
  levels = c("Normal_by_Numbat", "Cancer_by_Numbat")
)

message("Numbat status counts:")
print(table(seu$numbat_status_simplified, useNA = "ifany"))

message("HSC compartment by Numbat status:")
print(table(seu$hsc_compartment, seu$numbat_status_simplified, useNA = "ifany"))

# ============================================================
# Normal sample flag
# ============================================================

seu$sample_origin <- ifelse(
  grepl("normal|healthy|control", seu@meta.data[[sample_col]], ignore.case = TRUE),
  "Normal_sample",
  "Patient_sample"
)

seu$sample_origin <- factor(
  seu$sample_origin,
  levels = c("Normal_sample", "Patient_sample")
)

message("Sample origin:")
print(table(seu@meta.data[[sample_col]], seu$sample_origin))

# ============================================================
# Protein matrix
# ============================================================

DefaultAssay(seu) <- dsb_assay

adt_mat <- get_layer(seu, dsb_assay, "data")

markers <- rownames(adt_mat)
markers <- markers[!is_background_marker(markers)]

adt_mat <- adt_mat[markers, , drop = FALSE]

meta <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    broad = .data[[broad_col]],
    fine = if (fine_col %in% colnames(.)) .data[[fine_col]] else broad,
    hsc_compartment = seu$hsc_compartment,
    numbat_status = seu$numbat_status_simplified,
    sample_origin = seu$sample_origin,
    sample = .data[[sample_col]]
  )

# ============================================================
# 1. Which CITE markers are enriched in HSC/MPP vs other cells?
# ============================================================

hsc_cells <- rownames(seu@meta.data)[seu$hsc_compartment == "HSC_MPP"]
other_cells <- rownames(seu@meta.data)[seu$hsc_compartment == "Other"]

hsc_cells <- intersect(hsc_cells, colnames(adt_mat))
other_cells <- intersect(other_cells, colnames(adt_mat))

if (length(hsc_cells) == 0) {
  stop("No HSC_MPP cells found in adt_mat.")
}

if (length(other_cells) == 0) {
  stop("No Other cells found in adt_mat.")
}

hsc_marker_summary <- tibble(
  marker = markers,
  mean_HSC = Matrix::rowMeans(adt_mat[, hsc_cells, drop = FALSE]),
  mean_Other = Matrix::rowMeans(adt_mat[, other_cells, drop = FALSE]),
  pct_HSC_pos_0 = Matrix::rowMeans(adt_mat[, hsc_cells, drop = FALSE] > 0) * 100,
  pct_Other_pos_0 = Matrix::rowMeans(adt_mat[, other_cells, drop = FALSE] > 0) * 100,
  delta_HSC_vs_Other = mean_HSC - mean_Other,
  ratio_pct_HSC_vs_Other = (pct_HSC_pos_0 + 1) / (pct_Other_pos_0 + 1)
) %>%
  arrange(desc(delta_HSC_vs_Other))

write.csv(
  hsc_marker_summary,
  file.path(analysis_dir, "01_HSC_marker_presence_vs_other_cells.csv"),
  row.names = FALSE
)

p_hsc_markers <- hsc_marker_summary %>%
  slice_max(delta_HSC_vs_Other, n = 30) %>%
  mutate(marker = fct_reorder(marker, delta_HSC_vs_Other)) %>%
  ggplot(aes(x = delta_HSC_vs_Other, y = marker)) +
  geom_col(fill = "grey35") +
  theme_bw(base_size = 12) +
  labs(
    title = "CITE_DSB markers enriched in HSC/MPP compartment",
    x = "Mean DSB: HSC/MPP - other cells",
    y = NULL
  )

ggsave(
  file.path(analysis_dir, "01_HSC_marker_presence_vs_other_cells.pdf"),
  p_hsc_markers,
  width = 9,
  height = 8
)

# ============================================================
# 2. Differential CITE in HSC: Cancer vs Normal by Numbat
# ============================================================

hsc_meta <- meta %>%
  filter(
    hsc_compartment == "HSC_MPP",
    !is.na(numbat_status)
  )

hsc_cells_use <- intersect(hsc_meta$cell, colnames(adt_mat))
hsc_meta <- hsc_meta %>% filter(cell %in% hsc_cells_use)

message("HSC cells used for Numbat DE:")
print(table(hsc_meta$numbat_status, useNA = "ifany"))

if (length(unique(na.omit(hsc_meta$numbat_status))) < 2) {
  warning("HSC compartment does not contain both Cancer_by_Numbat and Normal_by_Numbat.")
}

hsc_numbat_de <- map_dfr(markers, function(mk) {
  
  x <- as.numeric(adt_mat[mk, hsc_cells_use, drop = TRUE])
  group <- hsc_meta$numbat_status
  
  tibble(
    marker = mk,
    mean_normal = mean(x[group == "Normal_by_Numbat"], na.rm = TRUE),
    mean_cancer = mean(x[group == "Cancer_by_Numbat"], na.rm = TRUE),
    median_normal = median(x[group == "Normal_by_Numbat"], na.rm = TRUE),
    median_cancer = median(x[group == "Cancer_by_Numbat"], na.rm = TRUE),
    pct_normal_pos_0 = mean(x[group == "Normal_by_Numbat"] > 0, na.rm = TRUE) * 100,
    pct_cancer_pos_0 = mean(x[group == "Cancer_by_Numbat"] > 0, na.rm = TRUE) * 100,
    delta_cancer_minus_normal = mean_cancer - mean_normal,
    p_value = safe_wilcox(x, group)
  )
}) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(padj, desc(abs(delta_cancer_minus_normal)))

write.csv(
  hsc_numbat_de,
  file.path(analysis_dir, "02_HSC_CITE_DE_Cancer_vs_Normal_by_Numbat.csv"),
  row.names = FALSE
)

p_de <- hsc_numbat_de %>%
  mutate(
    neglog10_padj = -log10(padj),
    significant = padj < 0.05 & abs(delta_cancer_minus_normal) > 0.25
  ) %>%
  ggplot(aes(
    x = delta_cancer_minus_normal,
    y = neglog10_padj,
    label = ifelse(significant, marker, NA)
  )) +
  geom_point(aes(color = significant), size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 40) +
  theme_bw(base_size = 12) +
  labs(
    title = "HSC/MPP CITE_DSB: cancer vs normal by Numbat",
    x = "Mean DSB cancer - normal",
    y = "-log10 adjusted P"
  )

ggsave(
  file.path(analysis_dir, "02_volcano_HSC_CITE_DE_Cancer_vs_Normal_by_Numbat.pdf"),
  p_de,
  width = 8,
  height = 7
)

# ============================================================
# 3. Confirm markers in cells from normal sample
# ============================================================

candidate_markers <- hsc_numbat_de %>%
  filter(
    padj < 0.05,
    delta_cancer_minus_normal > 0
  ) %>%
  slice_max(delta_cancer_minus_normal, n = 20) %>%
  pull(marker)

if (length(candidate_markers) == 0) {
  candidate_markers <- hsc_numbat_de %>%
    slice_max(delta_cancer_minus_normal, n = 20) %>%
    pull(marker)
}

candidate_markers <- intersect(candidate_markers, rownames(adt_mat))

normal_sample_cells <- rownames(seu@meta.data)[
  seu$sample_origin == "Normal_sample" &
    seu$hsc_compartment == "HSC_MPP"
]

patient_cancer_hsc_cells <- rownames(seu@meta.data)[
  seu$sample_origin == "Patient_sample" &
    seu$hsc_compartment == "HSC_MPP" &
    seu$numbat_status_simplified == "Cancer_by_Numbat"
]

patient_normal_hsc_cells <- rownames(seu@meta.data)[
  seu$sample_origin == "Patient_sample" &
    seu$hsc_compartment == "HSC_MPP" &
    seu$numbat_status_simplified == "Normal_by_Numbat"
]

validation_groups <- c(
  setNames(rep("Normal_sample_HSC", length(normal_sample_cells)), normal_sample_cells),
  setNames(rep("Patient_normal_by_Numbat_HSC", length(patient_normal_hsc_cells)), patient_normal_hsc_cells),
  setNames(rep("Patient_cancer_by_Numbat_HSC", length(patient_cancer_hsc_cells)), patient_cancer_hsc_cells)
)

validation_cells <- names(validation_groups)
validation_cells <- intersect(validation_cells, colnames(adt_mat))
validation_groups <- validation_groups[validation_cells]

cat("Candidate markers found:", length(candidate_markers), "\n")
cat("Validation cells found:", length(validation_cells), "\n")

if (length(candidate_markers) == 0) {
  stop("No candidate markers found in adt_mat.")
}

if (length(validation_cells) == 0) {
  stop("No validation cells found in adt_mat.")
}

validation_summary <- map_dfr(candidate_markers, function(mk) {
  
  x <- as.numeric(adt_mat[mk, validation_cells, drop = TRUE])
  group <- validation_groups[validation_cells]
  
  tibble(
    marker = mk,
    group = as.character(group),
    value = x
  ) %>%
    group_by(marker, group) %>%
    summarise(
      n_cells = n(),
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      pct_positive_0 = mean(value > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
})

write.csv(
  validation_summary,
  file.path(analysis_dir, "03_validation_candidate_markers_in_normal_sample_HSC.csv"),
  row.names = FALSE
)

validation_long <- map_dfr(candidate_markers, function(mk) {
  tibble(
    marker = mk,
    cell = validation_cells,
    group = as.character(validation_groups[validation_cells]),
    value = as.numeric(adt_mat[mk, validation_cells, drop = TRUE])
  )
})

p_validation <- validation_long %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Normal_sample_HSC",
        "Patient_normal_by_Numbat_HSC",
        "Patient_cancer_by_Numbat_HSC"
      )
    )
  ) %>%
  ggplot(aes(x = group, y = value)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.2) +
  facet_wrap(~ marker, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Candidate cancer-HSC CITE markers checked in true normal sample HSCs",
    x = NULL,
    y = "CITE_DSB"
  )

ggsave(
  file.path(analysis_dir, "03_validation_candidate_markers_in_normal_sample_HSC.pdf"),
  p_validation,
  width = 14,
  height = 10
)

# ============================================================
# 4. Heatmap of HSC candidate markers
# ============================================================

heat_cells <- validation_cells
heat_markers <- intersect(candidate_markers, rownames(adt_mat))

if (length(heat_markers) > 1 && length(heat_cells) > 1) {
  
  heat_mat <- as.matrix(adt_mat[heat_markers, heat_cells, drop = FALSE])
  
  ann_col <- data.frame(
    group = as.character(validation_groups[heat_cells]),
    sample = seu@meta.data[heat_cells, sample_col],
    broad = seu@meta.data[heat_cells, broad_col]
  )
  
  rownames(ann_col) <- heat_cells
  
  pdf(
    file.path(analysis_dir, "04_heatmap_candidate_CITE_markers_HSC.pdf"),
    width = 10,
    height = 8
  )
  
  pheatmap(
    heat_mat,
    scale = "row",
    show_colnames = FALSE,
    annotation_col = ann_col,
    main = "Candidate CITE markers in HSC/MPP cells"
  )
  
  dev.off()
}

# ============================================================
# 5. Feature plots for candidate markers
# ============================================================

DefaultAssay(seu) <- dsb_assay

if (length(candidate_markers) > 0) {
  
  pdf(
    file.path(analysis_dir, "05_featureplots_candidate_HSC_CITE_markers.pdf"),
    width = 12,
    height = 8
  )
  
  for (mk in candidate_markers[seq_len(min(12, length(candidate_markers)))]) {
    
    if (!mk %in% rownames(seu[[dsb_assay]])) next
    
    p <- FeaturePlot(
      seu,
      features = mk,
      reduction = reduction_to_use,
      order = TRUE
    ) +
      ggtitle(paste0(mk, " CITE_DSB"))
    
    print(p)
  }
  
  dev.off()
}

# ============================================================
# 6. BAFF-R RNA-protein correlation by cell type
# ============================================================

rna_gene <- "TNFRSF13C"
protein_patterns <- c("BAFF", "BAFF-R", "BAFFR", "TNFRSF13C")

protein_features <- rownames(seu[[dsb_assay]])

baffr_protein <- protein_features[
  grepl(
    paste(protein_patterns, collapse = "|"),
    protein_features,
    ignore.case = TRUE
  )
]

print(baffr_protein)

if (length(baffr_protein) == 0) {
  stop("Could not find BAFF-R protein feature in CITE assay.")
}

baffr_protein <- baffr_protein[1]

if (!rna_gene %in% rownames(seu[[rna_assay]])) {
  stop("Could not find TNFRSF13C in RNA assay.")
}

rna_vals <- GetAssayData(
  seu,
  assay = rna_assay,
  layer = "data"
)[rna_gene, ]

protein_vals <- GetAssayData(
  seu,
  assay = dsb_assay,
  layer = "data"
)[baffr_protein, ]

baffr_df <- tibble(
  cell = colnames(seu),
  celltype = seu@meta.data[[broad_col]],
  sample = seu@meta.data[[sample_col]],
  RNA_TNFRSF13C = as.numeric(rna_vals[colnames(seu)]),
  protein_BAFFR = as.numeric(protein_vals[colnames(seu)])
) %>%
  filter(!is.na(celltype))

cor_by_celltype <- baffr_df %>%
  group_by(celltype) %>%
  summarise(
    n_cells = n(),
    n_rna_positive = sum(RNA_TNFRSF13C > 0, na.rm = TRUE),
    n_protein_positive = sum(protein_BAFFR > 0, na.rm = TRUE),
    pct_rna_positive = mean(RNA_TNFRSF13C > 0, na.rm = TRUE) * 100,
    pct_protein_positive = mean(protein_BAFFR > 0, na.rm = TRUE) * 100,
    mean_rna = mean(RNA_TNFRSF13C, na.rm = TRUE),
    mean_protein = mean(protein_BAFFR, na.rm = TRUE),
    spearman_cor = suppressWarnings(
      cor(RNA_TNFRSF13C, protein_BAFFR, method = "spearman")
    ),
    pearson_cor = suppressWarnings(
      cor(RNA_TNFRSF13C, protein_BAFFR, method = "pearson")
    ),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 20) %>%
  arrange(desc(spearman_cor))

write.csv(
  cor_by_celltype,
  file.path(analysis_dir, "06_BAFFR_RNA_protein_correlation_by_celltype.csv"),
  row.names = FALSE
)

p_mean_scatter <- cor_by_celltype %>%
  ggplot(aes(
    x = mean_rna,
    y = mean_protein,
    label = celltype
  )) +
  geom_point(aes(size = n_cells), alpha = 0.8) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 12) +
  labs(
    title = "BAFF-R mean RNA vs mean protein by cell type",
    subtitle = paste0("RNA = TNFRSF13C; protein = ", baffr_protein),
    x = "Mean RNA expression: TNFRSF13C",
    y = "Mean CITE_DSB protein: BAFF-R",
    size = "Cells"
  )

ggsave(
  file.path(analysis_dir, "06_BAFFR_mean_RNA_vs_mean_protein_by_celltype.pdf"),
  p_mean_scatter,
  width = 8,
  height = 7
)

# ============================================================
# 7. Manual selected-marker CITE dotplot
# ============================================================

features_use <- c("CD34", "c-Kit", "BAFF-R", "CD54", "GPR56")
features_use <- intersect(features_use, rownames(adt_mat))

if (length(features_use) > 0) {
  
  meta_dot <- seu@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      celltype_broad_clean = as.character(.data[[broad_col]]),
      celltype_broad_clean = ifelse(
        is.na(celltype_broad_clean),
        "Unknown",
        celltype_broad_clean
      )
    ) %>%
    select(cell, celltype_broad_clean)
  
  df_dot <- as.data.frame(t(as.matrix(adt_mat[features_use, , drop = FALSE]))) %>%
    rownames_to_column("cell") %>%
    pivot_longer(
      cols = all_of(features_use),
      names_to = "marker",
      values_to = "value"
    ) %>%
    left_join(meta_dot, by = "cell") %>%
    filter(!is.na(celltype_broad_clean)) %>%
    group_by(celltype_broad_clean, marker) %>%
    summarise(
      avg_exp = mean(value, na.rm = TRUE),
      pct_exp = mean(value > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  p_manual_dot <- ggplot(
    df_dot,
    aes(x = marker, y = celltype_broad_clean)
  ) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_size(range = c(0, 6)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(
      title = "Selected CITE_DSB markers by broad cell type",
      x = NULL,
      y = NULL,
      size = "% positive",
      color = "Mean CITE_DSB"
    )
  
  ggsave(
    file.path(analysis_dir, "07_selected_CITE_markers_by_celltype_manual_dotplot.pdf"),
    p_manual_dot,
    width = 8,
    height = 7
  )
  
  print(p_manual_dot)
}

# ============================================================
# 8. Final summary
# ============================================================

summary_table <- tibble(
  metric = c(
    "Total cells",
    "HSC/MPP cells",
    "HSC/MPP normal by Numbat",
    "HSC/MPP cancer by Numbat",
    "HSC/MPP cells from normal sample",
    "Candidate cancer-HSC CITE markers"
  ),
  value = c(
    ncol(seu),
    sum(seu$hsc_compartment == "HSC_MPP", na.rm = TRUE),
    sum(
      seu$hsc_compartment == "HSC_MPP" &
        seu$numbat_status_simplified == "Normal_by_Numbat",
      na.rm = TRUE
    ),
    sum(
      seu$hsc_compartment == "HSC_MPP" &
        seu$numbat_status_simplified == "Cancer_by_Numbat",
      na.rm = TRUE
    ),
    sum(
      seu$hsc_compartment == "HSC_MPP" &
        seu$sample_origin == "Normal_sample",
      na.rm = TRUE
    ),
    length(candidate_markers)
  )
)

write.csv(
  summary_table,
  file.path(analysis_dir, "08_summary_table.csv"),
  row.names = FALSE
)

print(summary_table)

message("Done. Outputs written to: ", analysis_dir)