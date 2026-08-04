# ------------------------------------------------------------------
# LK1 re-sequenced - step 03 of 5
#
# Our depth against the AML scAtlas, plus the molecule-level saturation curves.
#
# Supersedes earlier revisions of the same analysis, which are not kept here:
#   scripts/1c.investigate_resequencing_subsampling.R  (2026-05-27, first pass)
# A LATER edit also exists, scripts/1c2....R (2026-06-16), but it produced no
# output - nothing under depth_benchmarking* postdates 2026-06-01 - so this is
# the revision that made the figures in the report. 1c2 is also shorter than
# this one, so it is not a strict superset; check it before reusing.
#
# Frozen for the lab archive 2026-08-03 from scripts/1c1.investigate_resequencing_subsampling.R (mtime 2026-06-01).
# md5 of the original: 4de3ad15758d530eab1109a6791266fe
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# ============================================================
# Cell Ranger downsampled outputs vs AML atlas benchmarking
# Use 100% run sample IDs to annotate all subsamples
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Seurat)
  library(patchwork)
  library(scales)
})

# ============================================================
# Paths
# ============================================================

downsample_base <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/results/cellranger_downsampled_reseq/260522_VH01624_461_222JLJVNX"

aml_intermediate_dir <- "../results/depth_benchmarking"

aml_qc_file <- file.path(
  aml_intermediate_dir,
  "aml_atlas_obs_metadata_with_qc.csv"
)

output_dir <- "../results/depth_benchmarking_cellranger_downsampled"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

clean_barcode <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace("-[0-9]+$", "")
}

find_cr_matrix_dir <- function(run_dir) {
  hits <- list.dirs(run_dir, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == "filtered_feature_bc_matrix"]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

find_cr_metrics_file <- function(run_dir) {
  hits <- list.files(
    run_dir,
    pattern = "^metrics_summary\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

read_cr_matrix_qc <- function(matrix_dir, dataset, fraction) {
  
  message("Reading matrix: ", matrix_dir)
  
  mat <- Seurat::Read10X(data.dir = matrix_dir)
  
  if (is.list(mat)) {
    mat <- mat[["Gene Expression"]]
  }
  
  tibble(
    dataset = dataset,
    source = "Cell Ranger downsampled",
    barcode = colnames(mat),
    barcode_clean = clean_barcode(colnames(mat)),
    fraction = fraction,
    nUMI = as.numeric(Matrix::colSums(mat)),
    nGene = as.integer(Matrix::colSums(mat > 0)),
    UMIs_per_gene = nUMI / pmax(nGene, 1)
  )
}

read_cr_metrics <- function(metrics_file, dataset, fraction) {
  
  if (is.na(metrics_file) || !file.exists(metrics_file)) {
    warning("Missing metrics_summary.csv for: ", dataset)
    return(tibble(
      dataset = dataset,
      source = "Cell Ranger downsampled",
      fraction = fraction
    ))
  }
  
  readr::read_csv(metrics_file, show_col_types = FALSE) %>%
    mutate(
      dataset = dataset,
      source = "Cell Ranger downsampled",
      fraction = fraction
    ) %>%
    relocate(dataset, source, fraction)
}

extract_metric_numeric <- function(metrics_df, pattern) {
  
  col <- names(metrics_df)[stringr::str_detect(
    names(metrics_df),
    stringr::regex(pattern, ignore_case = TRUE)
  )]
  
  if (length(col) == 0) return(NA_real_)
  
  x <- metrics_df[[col[1]]] %>%
    as.character() %>%
    stringr::str_remove_all(",") %>%
    stringr::str_remove_all("%")
  
  as.numeric(x)
}

summarise_per_cell_qc <- function(df) {
  
  df %>%
    group_by(dataset, source, fraction) %>%
    summarise(
      cells_detected = n(),
      median_genes = median(nGene, na.rm = TRUE),
      mean_genes = mean(nGene, na.rm = TRUE),
      p10_genes = quantile(nGene, 0.10, na.rm = TRUE),
      p90_genes = quantile(nGene, 0.90, na.rm = TRUE),
      median_UMI = median(nUMI, na.rm = TRUE),
      mean_UMI = mean(nUMI, na.rm = TRUE),
      p10_UMI = quantile(nUMI, 0.10, na.rm = TRUE),
      p90_UMI = quantile(nUMI, 0.90, na.rm = TRUE),
      median_UMIs_per_gene = median(UMIs_per_gene, na.rm = TRUE),
      mean_UMIs_per_gene = mean(UMIs_per_gene, na.rm = TRUE),
      p10_UMIs_per_gene = quantile(UMIs_per_gene, 0.10, na.rm = TRUE),
      p90_UMIs_per_gene = quantile(UMIs_per_gene, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
}

# ============================================================
# Discover Cell Ranger downsample folders
# ============================================================

outs <- tibble(
  run_dir = list.dirs(downsample_base, recursive = FALSE, full.names = TRUE),
  run_name = basename(run_dir)
) %>%
  filter(str_detect(run_name, "LK1-GEX_frac_")) %>%
  mutate(
    fraction_label = str_extract(run_name, "frac_[0-9]+"),
    fraction_num = as.numeric(str_remove(fraction_label, "frac_")),
    fraction = fraction_num / 100,
    matrix_dir = map_chr(run_dir, find_cr_matrix_dir),
    metrics_file = map_chr(run_dir, find_cr_metrics_file),
    has_matrix = !is.na(matrix_dir),
    has_metrics = !is.na(metrics_file)
  ) %>%
  filter(has_matrix) %>%
  arrange(fraction)

if (nrow(outs) == 0) {
  stop("No completed Cell Ranger filtered_feature_bc_matrix folders found under: ", downsample_base)
}

write_csv(
  outs,
  file.path(output_dir, "cellranger_downsampled_dirs_used.csv")
)

print(outs, n = Inf, width = Inf)

# ============================================================
# Read all Cell Ranger matrices
# ============================================================

cr_cell_qc_raw <- purrr::pmap_dfr(
  outs %>% select(matrix_dir, run_name, fraction),
  function(matrix_dir, run_name, fraction) {
    read_cr_matrix_qc(
      matrix_dir = matrix_dir,
      dataset = run_name,
      fraction = fraction
    )
  }
)

# ============================================================
# Build sample_name lookup from 100% subsample
# ============================================================

full_fraction <- max(cr_cell_qc_raw$fraction, na.rm = TRUE)

message("Using fraction ", full_fraction, " as the sample_name lookup.")

full_run <- outs %>%
  dplyr::filter(fraction == full_fraction) %>%
  dplyr::slice(1)

# This assumes the 100% output has per-sample Cell Ranger folders.
# Usually they are under something like:
# outs/per_sample_outs/<sample_name>/count/sample_filtered_feature_bc_matrix
full_sample_matrix_dirs <- list.dirs(
  full_run$run_dir,
  recursive = TRUE,
  full.names = TRUE
)

full_sample_matrix_dirs <- full_sample_matrix_dirs[
  basename(full_sample_matrix_dirs) %in% c(
    "sample_filtered_feature_bc_matrix",
    "filtered_feature_bc_matrix"
  )
]

# Avoid using the pooled matrix itself as a sample.
full_sample_matrix_dirs <- full_sample_matrix_dirs[
  stringr::str_detect(full_sample_matrix_dirs, "per_sample_outs")
]

if (length(full_sample_matrix_dirs) == 0) {
  stop(
    "Could not find per-sample matrices in the 100% run.\n",
    "Expected paths containing per_sample_outs/<sample_name>/count/sample_filtered_feature_bc_matrix.\n",
    "The pooled filtered_feature_bc_matrix alone does not contain sample_name."
  )
}

sample_lookup <- purrr::map_dfr(full_sample_matrix_dirs, function(matrix_dir) {
  
  sample_name <- stringr::str_match(
    matrix_dir,
    "per_sample_outs/([^/]+)/"
  )[, 2]
  
  mat <- Seurat::Read10X(data.dir = matrix_dir)
  
  if (is.list(mat)) {
    mat <- mat[["Gene Expression"]]
  }
  
  tibble(
    sample_name = sample_name,
    barcode = colnames(mat),
    barcode_clean = clean_barcode(colnames(mat))
  )
}) %>%
  distinct(barcode_clean, .keep_all = TRUE)

write_csv(
  sample_lookup,
  file.path(output_dir, "sample_name_lookup_from_100_percent_run.csv")
)

print(sample_lookup %>% count(sample_name), n = Inf)

# ============================================================
# Apply 100% sample_name lookup to all subsamples
# ============================================================

cr_cell_qc <- cr_cell_qc_raw %>%
  left_join(
    sample_lookup %>% select(barcode_clean, sample_name),
    by = "barcode_clean"
  )

missing_sample <- cr_cell_qc %>%
  filter(is.na(sample_name)) %>%
  count(fraction, name = "missing_cells")

if (nrow(missing_sample) > 0) {
  warning("Some cells did not match the 100% sample lookup.")
  print(missing_sample)
}

cr_cell_qc <- cr_cell_qc %>%
  mutate(
    sample_name = tidyr::replace_na(sample_name, "unmatched_sample")
  )

write_csv(
  cr_cell_qc,
  file.path(output_dir, "cellranger_downsampled_per_cell_qc_with_sample_name.csv")
)

# Check that this is now per sample_name, not pooled.
print(
  cr_cell_qc %>%
    count(fraction, sample_name) %>%
    arrange(fraction, sample_name),
  n = Inf
)

# ============================================================
# Read Cell Ranger metrics
# ============================================================

cr_metrics_raw <- purrr::pmap_dfr(
  outs %>% select(metrics_file, run_name, fraction),
  function(metrics_file, run_name, fraction) {
    read_cr_metrics(
      metrics_file = metrics_file,
      dataset = run_name,
      fraction = fraction
    )
  }
)

write_csv(
  cr_metrics_raw,
  file.path(output_dir, "cellranger_downsampled_metrics_raw.csv")
)

cr_metrics_clean <- cr_metrics_raw %>%
  rowwise() %>%
  mutate(
    estimated_cells_metric = extract_metric_numeric(cur_data(), "Estimated Number of Cells"),
    mean_reads_per_cell_metric = extract_metric_numeric(cur_data(), "Mean Reads per Cell"),
    median_genes_per_cell_metric = extract_metric_numeric(cur_data(), "Median Genes per Cell"),
    median_umi_per_cell_metric = extract_metric_numeric(cur_data(), "Median UMI Counts per Cell|Median UMI"),
    sequencing_saturation_metric = extract_metric_numeric(cur_data(), "Sequencing Saturation")
  ) %>%
  ungroup() %>%
  select(
    dataset,
    source,
    fraction,
    estimated_cells_metric,
    mean_reads_per_cell_metric,
    median_genes_per_cell_metric,
    median_umi_per_cell_metric,
    sequencing_saturation_metric
  )

write_csv(
  cr_metrics_clean,
  file.path(output_dir, "cellranger_downsampled_metrics_clean.csv")
)

# ============================================================
# Load AML atlas QC
# ============================================================

if (!file.exists(aml_qc_file)) {
  stop("Missing AML QC file: ", aml_qc_file)
}

aml_qc_raw <- readr::read_csv(
  aml_qc_file,
  show_col_types = FALSE
)

aml_qc <- aml_qc_raw %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(
    dataset = "AML scAtlas",
    source = "AML atlas",
    barcode = cell_id,
    fraction = 1,
    sample_name = donor_id,
    nGene = as.numeric(n_detected_genes),
    nUMI = as.numeric(total_counts),
    UMIs_per_gene = nUMI / pmax(nGene, 1)
  ) %>%
  select(
    dataset,
    source,
    barcode,
    fraction,
    sample_name,
    donor_id,
    nGene,
    nUMI,
    UMIs_per_gene,
    everything()
  )

# ============================================================
# Standard summaries
# ============================================================

cr_summary <- summarise_per_cell_qc(cr_cell_qc)

cr_summary <- cr_summary %>%
  left_join(
    cr_metrics_clean %>%
      select(dataset, fraction, mean_reads_per_cell_metric),
    by = c("dataset", "fraction")
  ) %>%
  mutate(
    mean_reads_per_cell = mean_reads_per_cell_metric
  ) %>%
  select(-mean_reads_per_cell_metric)

aml_summary <- summarise_per_cell_qc(aml_qc) %>%
  mutate(mean_reads_per_cell = NA_real_)

combined_summary <- bind_rows(cr_summary, aml_summary) %>%
  arrange(source, fraction)

combined_per_cell <- bind_rows(
  cr_cell_qc %>%
    select(dataset, source, barcode, fraction, sample_name, nUMI, nGene, UMIs_per_gene),
  aml_qc %>%
    select(dataset, source, barcode, fraction, sample_name, nUMI, nGene, UMIs_per_gene)
)

write_csv(cr_summary, file.path(output_dir, "cellranger_downsampled_summary.csv"))
write_csv(aml_summary, file.path(output_dir, "aml_atlas_summary.csv"))
write_csv(combined_summary, file.path(output_dir, "combined_cellranger_downsampled_vs_AML_summary.csv"))
write_csv(combined_per_cell, file.path(output_dir, "combined_cellranger_downsampled_vs_AML_per_cell_qc.csv"))

# ============================================================
# Correct cells-per-sample benchmark
# ============================================================

cr_cells_per_sample <- cr_cell_qc %>%
  filter(sample_name != "unmatched_sample") %>%
  group_by(fraction, sample_name) %>%
  summarise(
    cells = n(),
    .groups = "drop"
  )

aml_cells_per_donor <- aml_qc %>%
  filter(!is.na(donor_id)) %>%
  group_by(donor_id) %>%
  summarise(
    cells = n(),
    .groups = "drop"
  )

aml_median_cells_per_donor <- median(
  aml_cells_per_donor$cells,
  na.rm = TRUE
)

write_csv(
  cr_cells_per_sample,
  file.path(output_dir, "cellranger_cells_per_sample_name_per_subsample.csv")
)

write_csv(
  aml_cells_per_donor,
  file.path(output_dir, "aml_atlas_cells_per_donor_id.csv")
)

p_cells <- cr_cells_per_sample %>%
  mutate(
    fraction_label = paste0(fraction * 100, "%"),
    fraction_label = factor(
      fraction_label,
      levels = paste0(sort(unique(fraction)) * 100, "%")
    )
  ) %>%
  ggplot(aes(x = fraction_label, y = cells)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_point(
    aes(colour = sample_name),
    position = position_jitter(width = 0.12, height = 0),
    size = 2,
    alpha = 0.85
  ) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 1
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5
  ) +
  geom_hline(
    yintercept = aml_median_cells_per_donor,
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cell Ranger downsample fraction",
    y = "Cells detected per sample_name",
    title = "Cells detected per sample across Cell Ranger downsampled outputs",
    subtitle = paste0(
      "Boxplot = sample_name values per subsample; dashed line = AML atlas median cells per donor_id: ",
      scales::comma(round(aml_median_cells_per_donor))
    ),
    colour = "sample_name"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ============================================================
# Other reference values
# ============================================================

aml_ref <- aml_summary %>% slice(1)

aml_median_genes <- aml_ref$median_genes
aml_median_umi <- aml_ref$median_UMI
aml_median_umi_per_gene <- aml_ref$median_UMIs_per_gene

# ============================================================
# Other plots
# ============================================================

plot_curve_with_aml <- function(summary_df, y_col, aml_value, y_label, title) {
  
  summary_df %>%
    filter(source == "Cell Ranger downsampled") %>%
    ggplot(aes(x = fraction, y = .data[[y_col]])) +
    geom_hline(yintercept = aml_value, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_bw(base_size = 12) +
    labs(
      x = "Cell Ranger downsample fraction",
      y = y_label,
      title = title,
      subtitle = "Dashed line = AML atlas reference"
    )
}

plot_distribution_with_aml <- function(per_cell_df, y_col, y_label, title) {
  
  available_depths <- per_cell_df %>%
    filter(source == "Cell Ranger downsampled") %>%
    distinct(fraction) %>%
    arrange(fraction) %>%
    pull(fraction)
  
  df <- bind_rows(
    per_cell_df %>%
      filter(source == "Cell Ranger downsampled") %>%
      mutate(depth_label = paste0(fraction * 100, "%")),
    per_cell_df %>%
      filter(source == "AML atlas") %>%
      mutate(depth_label = "AML atlas")
  ) %>%
    mutate(
      depth_label = factor(
        depth_label,
        levels = c(paste0(available_depths * 100, "%"), "AML atlas")
      )
    )
  
  ymax <- quantile(df[[y_col]], 0.99, na.rm = TRUE)
  
  ggplot(df, aes(x = depth_label, y = .data[[y_col]], fill = depth_label)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_bw(base_size = 12) +
    labs(
      x = NULL,
      y = y_label,
      title = title
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

p_genes <- plot_curve_with_aml(
  combined_summary,
  "median_genes",
  aml_median_genes,
  "Median genes per cell",
  "Genes per cell across Cell Ranger downsampled outputs"
)

p_umi <- plot_curve_with_aml(
  combined_summary,
  "median_UMI",
  aml_median_umi,
  "Median UMIs per cell",
  "UMIs per cell across Cell Ranger downsampled outputs"
)

p_reads <- combined_summary %>%
  filter(source == "Cell Ranger downsampled") %>%
  ggplot(aes(x = fraction, y = mean_reads_per_cell)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cell Ranger downsample fraction",
    y = "Mean reads per cell",
    title = "Reads per cell across Cell Ranger downsampled outputs",
    subtitle = "Reads per cell unavailable for AML atlas h5ad count matrix"
  )

p_umi_per_gene <- plot_curve_with_aml(
  combined_summary,
  "median_UMIs_per_gene",
  aml_median_umi_per_gene,
  "Median UMIs per detected gene",
  "UMIs per gene across Cell Ranger downsampled outputs"
)

p_genes_dist <- plot_distribution_with_aml(
  combined_per_cell,
  "nGene",
  "Genes per cell",
  "Genes per cell: Cell Ranger downsampled outputs vs AML atlas"
)

p_umi_dist <- plot_distribution_with_aml(
  combined_per_cell,
  "nUMI",
  "UMIs per cell",
  "UMIs per cell: Cell Ranger downsampled outputs vs AML atlas"
)

p_umi_per_gene_dist <- plot_distribution_with_aml(
  combined_per_cell,
  "UMIs_per_gene",
  "UMIs per detected gene",
  "UMIs per gene: Cell Ranger downsampled outputs vs AML atlas"
)

# ============================================================
# Save plots
# ============================================================

ggsave(file.path(output_dir, "boxplot_cells_per_sample_name_per_subsample_vs_AML.pdf"), p_cells, width = 8, height = 5)
ggsave(file.path(output_dir, "curve_genes_per_cell_vs_AML.pdf"), p_genes, width = 7, height = 5)
ggsave(file.path(output_dir, "curve_UMIs_per_cell_vs_AML.pdf"), p_umi, width = 7, height = 5)
ggsave(file.path(output_dir, "curve_reads_per_cell.pdf"), p_reads, width = 7, height = 5)
ggsave(file.path(output_dir, "curve_UMIs_per_gene_vs_AML.pdf"), p_umi_per_gene, width = 7, height = 5)

ggsave(file.path(output_dir, "distribution_genes_per_cell_vs_AML.pdf"), p_genes_dist, width = 9, height = 5)
ggsave(file.path(output_dir, "distribution_UMIs_per_cell_vs_AML.pdf"), p_umi_dist, width = 9, height = 5)
ggsave(file.path(output_dir, "distribution_UMIs_per_gene_vs_AML.pdf"), p_umi_per_gene_dist, width = 9, height = 5)

ggsave(
  file.path(output_dir, "combined_curves_cellranger_downsampled_vs_AML.pdf"),
  p_cells / p_genes / p_umi / p_reads / p_umi_per_gene,
  width = 8,
  height = 18
)

ggsave(
  file.path(output_dir, "combined_distributions_cellranger_downsampled_vs_AML.pdf"),
  p_genes_dist / p_umi_dist / p_umi_per_gene_dist,
  width = 9,
  height = 14
)

print(cr_cells_per_sample, n = Inf)
print(aml_cells_per_donor)
print(combined_summary)
print(cr_metrics_clean)

message("Analysis complete. Outputs written to: ", output_dir)

