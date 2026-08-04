# ------------------------------------------------------------------
# LK1 re-sequenced - step 02 of 5
#
# First run vs re-sequenced: cells, genes and UMIs per cell, and the gain summary.
#
# Frozen for the lab archive 2026-08-03 from scripts/1b.investigate_resequencing.R (mtime 2026-06-01).
# md5 of the original: 261dd563cda48582536484f02b4f29d1
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

library(tidyverse)
library(rhdf5)
library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(scater)
library(patchwork)
library(Seurat)

outs <- list(
  resequenced = "../results/cellranger_withbam_reseq/260522_VH01624_461_222JLJVNX/LK1-GEX/outs",
  first_run   = "../results/cellranger/260423_VH01624_453_222HWMYNX/LK1-GEX/outs"
)

aml_h5ad <- "../results/AML_atlas/6e37b9b1-185c-4505-9842-8138157c1923.h5ad"

output_dir <- "../results/depth_benchmarking"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_obs_column <- function(col) {
  read_h5ad_col(aml_h5ad, file.path("obs", col))
}

read_h5ad_obs_col <- function(file, col) {
  
  base <- paste0("/obs/", col)
  
  # 1. Try categorical AnnData format: /obs/col/codes + /obs/col/categories
  out <- tryCatch({
    
    codes <- rhdf5::h5read(file, paste0(base, "/codes"))
    cats  <- rhdf5::h5read(file, paste0(base, "/categories"))
    
    values <- rep(NA_character_, length(codes))
    keep <- !is.na(codes) & codes >= 0
    values[keep] <- as.character(cats[codes[keep] + 1])
    values
    
  }, error = function(e) NULL)
  
  if (!is.null(out)) return(out)
  
  # 2. Otherwise try plain dataset, e.g. /obs/_index
  out <- tryCatch({
    rhdf5::h5read(file, base)
  }, error = function(e) NULL)
  
  if (!is.null(out)) return(out)
  
  stop("Could not read /obs column: ", col)
}

obs_cols <- c(
  "_index",
  "Study",
  "Sample",
  "donor_id",
  "cell_type",
  "Author Cell Type",
  "HSPC Cell Type",
  "disease",
  "assay",
  "sex",
  "age_floor",
  "Cytogenetics",
  "ELN Classification",
  "ELN Risk Group",
  "Mutations",
  "Translocations"
)

obs <- purrr::map(obs_cols, function(x) {
  message("Reading: ", x)
  read_h5ad_obs_col(aml_h5ad, x)
})

names(obs) <- make.names(obs_cols)

obs <- tibble::as_tibble(obs) %>%
  dplyr::rename(cell_id = X_index)

glimpse(obs)

write.csv(obs, file.path(output_dir, "aml_atlas_obs_metadata_selected.csv"))

message("Reading X/data and X/indptr...")

x_data   <- h5read(aml_h5ad, "raw/X/data")
x_indptr <- h5read(aml_h5ad, "raw/X/indptr")

n_cells <- length(x_indptr) - 1

stopifnot(nrow(obs) == n_cells)

detected_genes <- diff(x_indptr)

total_counts <- vapply(seq_len(n_cells), function(i) {
  start <- x_indptr[i] + 1
  end   <- x_indptr[i + 1]
  
  if (end < start) return(0)
  sum(x_data[start:end])
}, numeric(1))

obs_qc <- obs %>%
  mutate(
    n_detected_genes = detected_genes,
    total_counts = total_counts
  )

write.csv(obs_qc, file.path(output_dir, "aml_atlas_obs_metadata_with_qc.csv"))
# ============================================================
# Donor-level summary
# ============================================================

donor_summary <- obs_qc %>%
  dplyr::filter(!is.na(donor_id)) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    median_total_counts = median(total_counts, na.rm = TRUE),
    mean_total_counts = mean(total_counts, na.rm = TRUE),
    
    median_detected_genes = median(n_detected_genes, na.rm = TRUE),
    mean_detected_genes = mean(n_detected_genes, na.rm = TRUE),
    
    Study = Study[1],
    Sample = Sample[1],
    assay = assay[1],
    disease = disease[1],
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pct_cells = 100 * n_cells / sum(n_cells)
  ) %>%
  dplyr::arrange(desc(n_cells))

readr::write_csv(
  donor_summary,
  file.path(output_dir, "aml_atlas_donor_summary.csv")
)

readr::write_csv(
  donor_summary %>%
    dplyr::select(donor_id, n_cells, pct_cells, Study, Sample, assay, disease),
  file.path(output_dir, "aml_atlas_cells_per_donor.csv")
)

print(donor_summary, n = 100)

# ============================================================
# Cell counts per donor: by Study
# ============================================================

cells_per_donor_study <- obs_qc %>%
  dplyr::filter(!is.na(donor_id), !is.na(Study)) %>%
  dplyr::group_by(Study, donor_id) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  )

study_summary <- cells_per_donor_study %>%
  dplyr::group_by(Study) %>%
  dplyr::summarise(
    n_donors = dplyr::n(),
    median_cells_per_donor = median(n_cells, na.rm = TRUE),
    mean_cells_per_donor = mean(n_cells, na.rm = TRUE),
    p10_cells_per_donor = quantile(n_cells, 0.10, na.rm = TRUE),
    p90_cells_per_donor = quantile(n_cells, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(median_cells_per_donor))

readr::write_csv(
  cells_per_donor_study,
  file.path(output_dir, "aml_atlas_cells_per_donor_by_study.csv")
)

readr::write_csv(
  study_summary,
  file.path(output_dir, "aml_atlas_cell_count_summary_by_study.csv")
)

study_order <- study_summary %>%
  dplyr::arrange(median_cells_per_donor) %>%
  dplyr::pull(Study)

p_cells_by_study <- cells_per_donor_study %>%
  dplyr::mutate(
    Study = factor(Study, levels = study_order)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Study,
      y = n_cells
    )
  ) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.55,
    size = 1.8
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "AML Atlas: cells per donor by study",
    x = "Study",
    y = "Cells per donor"
  )

ggplot2::ggsave(
  file.path(output_dir, "boxplot_cells_per_donor_by_study.pdf"),
  p_cells_by_study,
  width = 6,
  height = 8
)

# ============================================================
# Cells per donor: by cell type
# ============================================================

cells_per_donor_cell_type <- obs_qc %>%
  dplyr::filter(!is.na(donor_id), !is.na(cell_type)) %>%
  dplyr::group_by(cell_type, donor_id) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  )

cell_type_summary <- cells_per_donor_cell_type %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    n_donors = dplyr::n(),
    median_cells_per_donor = median(n_cells, na.rm = TRUE),
    mean_cells_per_donor = mean(n_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(median_cells_per_donor))

readr::write_csv(
  cells_per_donor_cell_type,
  file.path(output_dir, "aml_atlas_cells_per_donor_by_cell_type.csv")
)

readr::write_csv(
  cell_type_summary,
  file.path(output_dir, "aml_atlas_cell_count_summary_by_cell_type.csv")
)

cell_type_order <- cell_type_summary %>%
  dplyr::arrange(median_cells_per_donor) %>%
  dplyr::pull(cell_type)

p_cells_by_cell_type <- cells_per_donor_cell_type %>%
  dplyr::mutate(
    cell_type = factor(cell_type, levels = cell_type_order)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = cell_type,
      y = n_cells
    )
  ) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.55,
    size = 1.8
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "AML Atlas: cells per donor by cell type",
    x = "Cell type",
    y = "Cells per donor"
  )

ggplot2::ggsave(
  file.path(output_dir, "boxplot_cells_per_donor_by_cell_type.pdf"),
  p_cells_by_cell_type,
  width = 8,
  height = 6
)
# ============================================================
# Cell counts by cell type
# ============================================================

cell_type_counts <- obs_qc %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::count(cell_type, name = "n_cells") %>%
  dplyr::mutate(
    pct_cells = 100 * n_cells / sum(n_cells)
  ) %>%
  dplyr::arrange(desc(n_cells))

readr::write_csv(
  cell_type_counts,
  file.path(output_dir, "aml_atlas_cell_counts_by_cell_type.csv")
)

p_cell_type_counts <- cell_type_counts %>%
  dplyr::mutate(
    cell_type = factor(cell_type, levels = rev(cell_type))
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = cell_type,
      y = n_cells
    )
  ) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "AML Atlas: cell counts by cell type",
    x = "Cell type",
    y = "Number of cells"
  )

ggplot2::ggsave(
  file.path(output_dir, "barplot_cell_counts_by_cell_type.pdf"),
  p_cell_type_counts,
  width = 8,
  height = 6
)


# ============================================================
# Cell counts per donor: by assay / technology
# ============================================================

cells_per_donor_assay <- obs_qc %>%
  dplyr::filter(!is.na(donor_id), !is.na(assay)) %>%
  dplyr::group_by(assay, donor_id) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  )

assay_summary <- cells_per_donor_assay %>%
  dplyr::group_by(assay) %>%
  dplyr::summarise(
    n_donors = dplyr::n(),
    median_cells_per_donor = median(n_cells, na.rm = TRUE),
    mean_cells_per_donor = mean(n_cells, na.rm = TRUE),
    p10_cells_per_donor = quantile(n_cells, 0.10, na.rm = TRUE),
    p90_cells_per_donor = quantile(n_cells, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(median_cells_per_donor))

readr::write_csv(
  cells_per_donor_assay,
  file.path(output_dir, "aml_atlas_cells_per_donor_by_assay.csv")
)

readr::write_csv(
  assay_summary,
  file.path(output_dir, "aml_atlas_cell_count_summary_by_assay.csv")
)

assay_order <- assay_summary %>%
  dplyr::arrange(median_cells_per_donor) %>%
  dplyr::pull(assay)

p_cells_by_assay <- cells_per_donor_assay %>%
  dplyr::mutate(
    assay = factor(assay, levels = assay_order)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = assay,
      y = n_cells
    )
  ) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.55,
    size = 1.8
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "Cells per donor by assay",
    x = "Technology / assay",
    y = "Cells per donor"
  )

ggplot2::ggsave(
  file.path(output_dir, "boxplot_cells_per_donor_by_assay.pdf"),
  p_cells_by_assay,
  width = 4,
  height = 5
)


# ============================================================
# Generic metadata QC summaries
# ============================================================

summarise_qc <- function(df, group_col) {
  
  df %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::group_by(group = .data[[group_col]]) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_total_counts = median(total_counts, na.rm = TRUE),
      mean_total_counts = mean(total_counts, na.rm = TRUE),
      median_detected_genes = median(n_detected_genes, na.rm = TRUE),
      mean_detected_genes = mean(n_detected_genes, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      metadata_column = group_col,
      group = as.character(group)
    ) %>%
    dplyr::select(
      metadata_column,
      group,
      n_cells,
      median_total_counts,
      mean_total_counts,
      median_detected_genes,
      mean_detected_genes
    )
}

qc_by_group <- dplyr::bind_rows(
  summarise_qc(obs_qc, "Study"),
  summarise_qc(obs_qc, "Sample"),
  summarise_qc(obs_qc, "donor_id"),
  summarise_qc(obs_qc, "cell_type"),
  summarise_qc(obs_qc, "Author.Cell.Type"),
  summarise_qc(obs_qc, "HSPC.Cell.Type"),
  summarise_qc(obs_qc, "assay")
)

readr::write_csv(
  qc_by_group,
  file.path(output_dir, "aml_atlas_qc_by_metadata_group.csv")
)

qc_by_group %>%
  dplyr::arrange(metadata_column, median_total_counts) %>%
  print(n = 100)
# =========================================================
# Seurat-style violin QC plots: genes / UMIs per cell
# =========================================================
aml_qc <- obs_qc %>%
  dplyr::transmute(
    dataset = "AML scAtlas",
    barcode = cell_id,
    fraction = 1,
    nGene = n_detected_genes,
    nUMI = total_counts
  )

comparison_full_qc <- dplyr::bind_rows(
  own_full_qc %>%
    dplyr::select(dataset, barcode, fraction, nGene, nUMI),
  aml_qc
)


plot_ylim_gene <- quantile(comparison_full_qc$nGene, 0.99, na.rm = TRUE)
plot_ylim_umi  <- quantile(comparison_full_qc$nUMI,  0.99, na.rm = TRUE)

p_vln_genes <- ggplot(
  comparison_full_qc,
  aes(x = dataset, y = nGene, fill = dataset)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
  coord_cartesian(ylim = c(0, plot_ylim_gene)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Detected genes per cell",
    x = NULL,
    y = "Genes per cell"
  )

p_vln_umis <- ggplot(
  comparison_full_qc,
  aes(x = dataset, y = nUMI, fill = dataset)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
  coord_cartesian(ylim = c(0, plot_ylim_umi)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "UMIs per cell",
    x = NULL,
    y = "UMIs per cell"
  )

p_vln_qc <- p_vln_genes + p_vln_umis

ggsave(
  file.path(output_dir, "violin_genes_UMIs_per_cell_full_depth_vs_AML.pdf"),
  p_vln_qc,
  width = 12,
  height = 4.5
)

ggsave(
  file.path(output_dir, "violin_genes_per_cell_full_depth_vs_AML.pdf"),
  p_vln_genes,
  width = 5,
  height = 4.5
)

ggsave(
  file.path(output_dir, "violin_UMIs_per_cell_full_depth_vs_AML.pdf"),
  p_vln_umis
  ,
  width = 5,
  height = 4.5
)

aml_qc <- obs_qc %>%
  dplyr::filter(
    assay %in% c("10x 3' v3", "10x 3' v3.1", "10x 3 prime v3")
  ) %>%
  dplyr::transmute(
    dataset = "AML scAtlas 10x 3' v3",
    barcode = cell_id,
    fraction = 1,
    nGene = n_detected_genes,
    nUMI = total_counts
  )

comparison_full_qc <- dplyr::bind_rows(
  own_full_qc %>%
    dplyr::select(dataset, barcode, fraction, nGene, nUMI),
  aml_qc
) %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = c("AML scAtlas 10x 3' v3", "first_run", "resequenced")
    )
  )


plot_ylim_gene <- quantile(comparison_full_qc$nGene, 0.99, na.rm = TRUE)
plot_ylim_umi  <- quantile(comparison_full_qc$nUMI,  0.99, na.rm = TRUE)

p_vln_genes <- ggplot(
  comparison_full_qc,
  aes(x = dataset, y = nGene, fill = dataset)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
  coord_cartesian(ylim = c(0, plot_ylim_gene)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Detected genes per cell",
    x = NULL,
    y = "Genes per cell"
  )

p_vln_umis <- ggplot(
  comparison_full_qc,
  aes(x = dataset, y = nUMI, fill = dataset)
) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
  coord_cartesian(ylim = c(0, plot_ylim_umi)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "UMIs per cell",
    x = NULL,
    y = "UMIs per cell"
  )

p_vln_qc <- p_vln_genes + p_vln_umis

ggsave(
  file.path(output_dir, "violin_genes_UMIs_per_cell_full_depth_vs_AML_10x3v3.pdf"),
  p_vln_qc,
  width = 12,
  height = 4.5
)

ggsave(
  file.path(output_dir, "violin_genes_per_cell_full_depth_vs_AML_10x3v3.pdf"),
  p_vln_genes,
  width = 5,
  height = 4.5
)

ggsave(
  file.path(output_dir, "violin_UMIs_per_cell_full_depth_vs_AML_10x3v3.pdf"),
  p_vln_umis
  ,
  width = 5,
  height = 4.5
)

# ============================================================
# Generic QC boxplots
# ============================================================

plot_qc_box <- function(df, group_col, y_col, file_stub) {
  
  safe_group_col <- stringr::str_replace_all(group_col, "[^A-Za-z0-9_]+", "_")
  
  p <- df %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = reorder(as.character(.data[[group_col]]), .data[[y_col]], median),
        y = .data[[y_col]]
      )
    ) +
    ggplot2::geom_boxplot(outlier.size = 0.2) +
    scale_y_log10() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = paste(y_col, "by", group_col),
      x = group_col,
      y = y_col
    )
  
  ggplot2::ggsave(
    file.path(output_dir, paste0(file_stub, "_by_", safe_group_col, ".pdf")),
    p,
    width = 8,
    height = 10
  )
}

plot_cols <- c(
  "Study",
  "Sample",
  "donor_id",
  "cell_type",
  "Author.Cell.Type",
  "HSPC.Cell.Type",
  "assay"
)

for (cc in plot_cols) {
  plot_qc_box(obs_qc, cc, "total_counts", "total_counts")
  plot_qc_box(obs_qc, cc, "n_detected_genes", "detected_genes")
}

# ============================================================
# Identify shallow metadata groups
# ============================================================

shallow_groups <- qc_by_group %>%
  dplyr::filter(n_cells >= 50) %>%
  dplyr::arrange(median_total_counts)

readr::write_csv(
  shallow_groups,
  file.path(output_dir, "aml_atlas_shallow_groups_ranked.csv")
)

print(shallow_groups, n = 50)

message("Done. Outputs written to: ", output_dir)

write.csv(
  obs,
  file.path(output_dir, "aml_atlas_obs_metadata_selected.csv")
)

write.csv(
  obs_qc,
  file.path(output_dir, "aml_atlas_obs_metadata_with_qc.csv")
)


