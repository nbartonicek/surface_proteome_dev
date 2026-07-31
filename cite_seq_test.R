# ============================================================
# CITE_DSB AML-vs-normal markers per specific cell type
# LK1 + LK2 combined
# Normal = normal sample OR Numbat-normal cells
# AML    = Numbat-tumor cells
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(pheatmap)
})

# ----------------------------
# Settings
# ----------------------------

cite_assay <- "CITE_DSB"
fine_col <- "predicted_CellType"
sample_col <- "sample_name"

normal_samples <- c("normal-01", "normal", "Normal", "NORMAL")

min_cells_per_group <- 10
min_total_cells <- 30
top_n_per_celltype <- 5

out_cite_dir <- file.path(out_dir, "CITE_DSB_AML_vs_Normal_by_specific_celltype_LK1_LK2")
dir.create(out_cite_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Input RDS files
# ----------------------------

numbat_files <- tibble(
  run_id = c("LK1", "LK2"),
  rds = c(
    file.path(
      proj,
      "results/seurat_annotated",
      "260423_VH01624_453_222HWMYNX",
      "numbat",
      "LK1_projected_CITE_DSB_Numbat_integrated.rds"
    ),
    file.path(
      proj,
      "results/seurat_annotated",
      "260528_VH01624_464_222K7VKNX",
      "numbat",
      "LK2_projected_CITE_DSB_Numbat_integrated.rds"
    )
  )
)

missing_files <- numbat_files$rds[!file.exists(numbat_files$rds)]

if (length(missing_files) > 0) {
  stop(
    "Missing Numbat-integrated RDS file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# ----------------------------
# Load and check objects
# ----------------------------

required_numbat_cols <- c("numbat_compartment", "numbat_call", "numbat_clone")

seu_list <- purrr::map2(
  numbat_files$rds,
  numbat_files$run_id,
  function(rds, run_id) {
    
    message("Loading ", run_id, ": ", rds)
    
    obj <- readRDS(rds)
    
    obj$run_id <- run_id
    
    missing_cols <- setdiff(required_numbat_cols, colnames(obj@meta.data))
    
    if (length(missing_cols) > 0) {
      stop(
        "Numbat metadata is missing from ", run_id, ".\n",
        "Missing columns: ", paste(missing_cols, collapse = ", ")
      )
    }
    
    if (!cite_assay %in% Assays(obj)) {
      stop(cite_assay, " assay is missing from ", run_id)
    }
    
    if (!fine_col %in% colnames(obj@meta.data)) {
      stop(fine_col, " is missing from metadata in ", run_id)
    }
    
    if (!sample_col %in% colnames(obj@meta.data)) {
      stop(sample_col, " is missing from metadata in ", run_id)
    }
    
    cat("\n", run_id, " Numbat compartment:\n", sep = "")
    print(table(obj$numbat_compartment, useNA = "ifany"))
    
    cat("\n", run_id, " Numbat call:\n", sep = "")
    print(table(obj$numbat_call, useNA = "ifany"))
    
    cat("\n", run_id, " Numbat clone:\n", sep = "")
    print(table(obj$numbat_clone, useNA = "ifany"))
    
    obj
  }
)

names(seu_list) <- numbat_files$run_id

# ----------------------------
# Merge LK1 + LK2
# ----------------------------

seu <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "LK1_LK2_CITE_Numbat"
)

DefaultAssay(seu) <- cite_assay

cat("\nMerged object:\n")
print(seu)

cat("\nCells by run:\n")
print(table(seu$run_id, useNA = "ifany"))

cat("\nCells by sample:\n")
print(table(seu@meta.data[[sample_col]], useNA = "ifany"))

# ----------------------------
# Define AML vs normal status
# ----------------------------

seu$AML_status_numbat <- dplyr::case_when(
  as.character(seu@meta.data[[sample_col]]) %in% normal_samples ~ "Normal",
  tolower(as.character(seu$numbat_compartment)) == "normal" ~ "Normal",
  tolower(as.character(seu$numbat_compartment)) == "tumor" ~ "AML",
  TRUE ~ NA_character_
)

cat("\nAML_status_numbat table:\n")
print(table(seu$AML_status_numbat, useNA = "ifany"))

cat("\nAML_status_numbat by run:\n")
print(table(seu$run_id, seu$AML_status_numbat, useNA = "ifany"))

write.csv(
  seu@meta.data,
  file.path(out_cite_dir, "metadata_with_AML_status_numbat_LK1_LK2.csv")
)

# ----------------------------
# Get CITE_DSB matrix
# ----------------------------

cite_mat <- tryCatch(
  GetAssayData(seu, assay = cite_assay, layer = "data"),
  error = function(e) {
    GetAssayData(seu, assay = cite_assay, slot = "data")
  }
)

if (nrow(cite_mat) == 0) {
  cite_mat <- tryCatch(
    GetAssayData(seu, assay = cite_assay, layer = "counts"),
    error = function(e) {
      GetAssayData(seu, assay = cite_assay, slot = "counts")
    }
  )
}

cite_features <- rownames(cite_mat)

# ----------------------------
# Cell type QC summary
# ----------------------------

celltype_status_counts <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(
    run_id,
    specific_celltype = .data[[fine_col]],
    AML_status_numbat,
    name = "n_cells"
  ) %>%
  arrange(run_id, specific_celltype, AML_status_numbat)

write.csv(
  celltype_status_counts,
  file.path(out_cite_dir, "celltype_AML_normal_counts_by_run.csv"),
  row.names = FALSE
)

celltype_status_counts_combined <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(
    specific_celltype = .data[[fine_col]],
    AML_status_numbat,
    name = "n_cells"
  ) %>%
  arrange(specific_celltype, AML_status_numbat)

write.csv(
  celltype_status_counts_combined,
  file.path(out_cite_dir, "celltype_AML_normal_counts_combined.csv"),
  row.names = FALSE
)

print(celltype_status_counts_combined, n = 100)

# ----------------------------
# Per-specific-cell-type CITE testing
# ----------------------------

fine_types <- sort(unique(seu@meta.data[[fine_col]]))
fine_types <- fine_types[!is.na(fine_types)]

cite_de_list <- list()
skipped_list <- list()

for (ct in fine_types) {
  
  message("Testing: ", ct)
  
  cells_ct <- rownames(seu@meta.data)[
    seu@meta.data[[fine_col]] == ct
  ]
  
  status_ct <- as.character(seu@meta.data[cells_ct, "AML_status_numbat", drop = TRUE])
  
  keep_status <- !is.na(status_ct) & status_ct %in% c("AML", "Normal")
  
  cells_ct <- cells_ct[keep_status]
  status_ct <- status_ct[keep_status]
  
  n_aml <- sum(status_ct == "AML", na.rm = TRUE)
  n_norm <- sum(status_ct == "Normal", na.rm = TRUE)
  
  if (length(cells_ct) < min_total_cells) {
    skipped_list[[ct]] <- tibble(
      celltype = ct,
      reason = "too_few_total_valid_cells",
      n_AML = n_aml,
      n_Normal = n_norm,
      n_total = length(cells_ct)
    )
    next
  }
  
  if (n_aml < min_cells_per_group || n_norm < min_cells_per_group) {
    skipped_list[[ct]] <- tibble(
      celltype = ct,
      reason = "too_few_cells_in_one_group",
      n_AML = n_aml,
      n_Normal = n_norm,
      n_total = length(cells_ct)
    )
    next
  }
  
  mat_ct <- as.matrix(cite_mat[, cells_ct, drop = FALSE])
  
  res_ct <- purrr::map_dfr(cite_features, function(marker) {
    
    x_aml <- mat_ct[marker, status_ct == "AML"]
    x_norm <- mat_ct[marker, status_ct == "Normal"]
    
    wt <- tryCatch(
      wilcox.test(x_aml, x_norm),
      error = function(e) NULL
    )
    
    tibble(
      celltype = ct,
      marker = marker,
      n_AML = length(x_aml),
      n_Normal = length(x_norm),
      mean_AML = mean(x_aml, na.rm = TRUE),
      mean_Normal = mean(x_norm, na.rm = TRUE),
      median_AML = median(x_aml, na.rm = TRUE),
      median_Normal = median(x_norm, na.rm = TRUE),
      delta_mean = mean_AML - mean_Normal,
      delta_median = median_AML - median_Normal,
      p_value = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }) %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      direction = case_when(
        FDR < 0.05 & delta_mean > 0 ~ "Higher in AML",
        FDR < 0.05 & delta_mean < 0 ~ "Higher in normal",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(FDR, desc(abs(delta_mean)))
  
  cite_de_list[[ct]] <- res_ct
}

cite_de <- bind_rows(cite_de_list)
skipped_celltypes <- bind_rows(skipped_list)

write.csv(
  cite_de,
  file.path(out_cite_dir, "CITE_DSB_AML_vs_Normal_by_specific_celltype_LK1_LK2.csv"),
  row.names = FALSE
)

write.csv(
  skipped_celltypes,
  file.path(out_cite_dir, "skipped_specific_celltypes_LK1_LK2.csv"),
  row.names = FALSE
)

# ----------------------------
# Top AML-up markers per specific cell type
# ----------------------------

top_cite_by_celltype <- cite_de %>%
  filter(FDR < 0.05, delta_mean > 0) %>%
  group_by(celltype) %>%
  slice_max(
    order_by = delta_mean,
    n = top_n_per_celltype,
    with_ties = FALSE
  ) %>%
  ungroup()

if (nrow(top_cite_by_celltype) == 0) {
  
  warning("No FDR-significant AML-up CITE markers found. Using top positive delta_mean markers instead.")
  
  top_cite_by_celltype <- cite_de %>%
    filter(delta_mean > 0) %>%
    group_by(celltype) %>%
    slice_max(
      order_by = delta_mean,
      n = top_n_per_celltype,
      with_ties = FALSE
    ) %>%
    ungroup()
}

write.csv(
  top_cite_by_celltype,
  file.path(out_cite_dir, "top_AML_up_CITE_DSB_markers_by_specific_celltype_LK1_LK2.csv"),
  row.names = FALSE
)

top_markers <- unique(top_cite_by_celltype$marker)

# ----------------------------
# Heatmap: delta mean DSB, AML - Normal
# ----------------------------

heat_df <- cite_de %>%
  filter(marker %in% top_markers) %>%
  select(celltype, marker, delta_mean) %>%
  pivot_wider(
    names_from = marker,
    values_from = delta_mean,
    values_fill = 0
  )

heat_mat <- heat_df %>%
  column_to_rownames("celltype") %>%
  as.matrix()

heat_mat <- heat_mat[, colSums(abs(heat_mat)) > 0, drop = FALSE]

if (nrow(heat_mat) > 1 && ncol(heat_mat) > 1) {
  
  pdf(
    file.path(out_cite_dir, "heatmap_top_AML_up_CITE_DSB_markers_by_specific_celltype_LK1_LK2.pdf"),
    width = max(8, ncol(heat_mat) * 0.35),
    height = max(6, nrow(heat_mat) * 0.35)
  )
  
  pheatmap(
    heat_mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "AML-enriched CITE_DSB markers by specific cell type\nLK1 + LK2, AML - Normal mean DSB"
  )
  
  dev.off()
}

# ----------------------------
# Significant-only heatmap
# ----------------------------

sig_heat_df <- cite_de %>%
  filter(marker %in% top_markers) %>%
  mutate(delta_sig = if_else(FDR < 0.05, delta_mean, 0)) %>%
  select(celltype, marker, delta_sig) %>%
  pivot_wider(
    names_from = marker,
    values_from = delta_sig,
    values_fill = 0
  )

sig_heat_mat <- sig_heat_df %>%
  column_to_rownames("celltype") %>%
  as.matrix()

sig_heat_mat <- sig_heat_mat[, colSums(abs(sig_heat_mat)) > 0, drop = FALSE]

if (nrow(sig_heat_mat) > 1 && ncol(sig_heat_mat) > 1) {
  
  pdf(
    file.path(out_cite_dir, "heatmap_significant_only_AML_up_CITE_DSB_markers_by_specific_celltype_LK1_LK2.pdf"),
    width = max(8, ncol(sig_heat_mat) * 0.35),
    height = max(6, nrow(sig_heat_mat) * 0.35)
  )
  
  pheatmap(
    sig_heat_mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Significant AML-enriched CITE_DSB markers\nLK1 + LK2, non-significant values set to 0"
  )
  
  dev.off()
}

# ----------------------------
# Save merged object
# ----------------------------

saveRDS(
  seu,
  file.path(out_cite_dir, "seu_LK1_LK2_with_AML_status_numbat_for_CITE_DE.rds")
)

cat("\nDone.\n")
cat("Output directory:\n", out_cite_dir, "\n")