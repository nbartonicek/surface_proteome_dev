
# =========================================================
# Helper functions
# =========================================================

read_filtered_barcodes <- function(outdir) {
  
  barcodes <- read_tsv(
    file.path(outdir, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
    col_names = "barcode",
    show_col_types = FALSE
  ) %>%
    pull(barcode)
  
  barcodes
}

read_gene_names <- function(outdir) {
  
  features <- read_tsv(
    file.path(outdir, "filtered_feature_bc_matrix/features.tsv.gz"),
    col_names = c("gene_id", "gene_name", "feature_type"),
    show_col_types = FALSE
  )
  
  features
}

read_cr_qc_from_matrix <- function(outdir, run_name) {
  
  message("Reading filtered matrix for full-depth QC: ", run_name)
  
  mat <- Seurat::Read10X(
    data.dir = file.path(outdir, "filtered_feature_bc_matrix")
  )
  
  if (is.list(mat)) {
    mat <- mat[["Gene Expression"]]
  }
  
  tibble(
    dataset = run_name,
    barcode = colnames(mat),
    fraction = 1,
    nUMI = as.numeric(Matrix::colSums(mat)),
    nGene = as.integer(Matrix::colSums(mat > 0))
  )
}

read_cellranger_metrics <- function(outdir, run_name) {
  
  read_csv(
    file.path(outdir, "metrics_summary.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(dataset = run_name) %>%
    relocate(dataset)
}

read_molecule_info_for_filtered_cells <- function(outdir) {
  
  h5 <- file.path(outdir, "molecule_info.h5")
  
  message("Reading molecule_info.h5: ", h5)
  
  filtered_barcodes <- read_filtered_barcodes(outdir)
  
  h5_barcodes <- h5read(h5, "barcodes")
  
  h5_contents <- h5ls(h5)
  
  feature_path <- if (any(h5_contents$group == "/features" & h5_contents$name == "name")) {
    "features/name"
  } else if (any(h5_contents$group == "/features" & h5_contents$name == "id")) {
    "features/id"
  } else {
    stop("Could not find features/name or features/id in molecule_info.h5")
  }
  
  h5_features <- h5read(h5, feature_path)
  
  barcode_idx <- as.integer(h5read(h5, "barcode_idx")) + 1L
  feature_idx <- as.integer(h5read(h5, "feature_idx")) + 1L
  
  # This is UMI molecule count/reads-per-molecule-like field in Cell Ranger molecule info.
  # It is useful as a weight for downsampling molecules.
  mol_count <- as.numeric(h5read(h5, "count"))
  
  mol <- tibble(
    barcode = h5_barcodes[barcode_idx],
    gene = h5_features[feature_idx],
    count = mol_count
  )
  
  # Critical: barcode matching can fail if suffixes differ.
  # Try exact match first.
  n_exact <- sum(mol$barcode %in% filtered_barcodes)
  
  message("Exact barcode matches: ", n_exact)
  
  if (n_exact == 0) {
    
    message("No exact matches. Trying barcode suffix harmonisation...")
    
    mol <- mol %>%
      mutate(
        barcode_raw = barcode,
        barcode = str_replace(barcode, "-1$", "")
      )
    
    filtered_barcodes2 <- str_replace(filtered_barcodes, "-1$", "")
    
    mol <- mol %>%
      filter(barcode %in% filtered_barcodes2)
    
    filtered_barcodes_final <- filtered_barcodes2
    
  } else {
    
    mol <- mol %>%
      filter(barcode %in% filtered_barcodes)
    
    filtered_barcodes_final <- filtered_barcodes
  }
  
  message("Molecules retained after barcode filtering: ", nrow(mol))
  message("Cells represented in molecule table: ", n_distinct(mol$barcode))
  
  list(
    mol = mol,
    filtered_barcodes = filtered_barcodes_final
  )
}

simulate_depth_from_molecules <- function(
    mol,
    filtered_barcodes,
    fractions,
    seed = 1
) {
  
  set.seed(seed)
  
  if (nrow(mol) == 0) {
    stop("No molecules available after barcode filtering.")
  }
  
  map_dfr(fractions, function(f) {
    
    message("Simulating fraction: ", f)
    
    # For each original molecule with 'count' supporting reads,
    # probability the molecule remains observed after read downsampling:
    # P(observed) = 1 - (1 - f)^count
    mol_sub <- mol %>%
      mutate(
        keep_prob = 1 - (1 - f)^count,
        kept = runif(n()) < keep_prob
      ) %>%
      filter(kept)
    
    mol_sub %>%
      group_by(barcode) %>%
      summarise(
        nUMI = n(),
        nGene = n_distinct(gene),
        .groups = "drop"
      ) %>%
      right_join(
        tibble(barcode = filtered_barcodes),
        by = "barcode"
      ) %>%
      mutate(
        nUMI = replace_na(nUMI, 0L),
        nGene = replace_na(nGene, 0L),
        fraction = f
      )
  })
}

summarise_qc <- function(df) {
  
  df %>%
    group_by(dataset, fraction) %>%
    summarise(
      cells = n(),
      
      median_genes = median(nGene, na.rm = TRUE),
      mean_genes = mean(nGene, na.rm = TRUE),
      p10_genes = quantile(nGene, 0.10, na.rm = TRUE),
      p90_genes = quantile(nGene, 0.90, na.rm = TRUE),
      
      median_UMI = median(nUMI, na.rm = TRUE),
      mean_UMI = mean(nUMI, na.rm = TRUE),
      p10_UMI = quantile(nUMI, 0.10, na.rm = TRUE),
      p90_UMI = quantile(nUMI, 0.90, na.rm = TRUE),
      
      .groups = "drop"
    )
}

# =========================================================
# Run full-depth Cell Ranger QC
# =========================================================

own_full_qc <- imap_dfr(
  outs,
  ~ read_cr_qc_from_matrix(.x, .y)
)

write_csv(
  own_full_qc,
  file.path(output_dir, "own_full_depth_QC_from_filtered_matrix.csv")
)

# =========================================================
# Run molecule-level downsampling
# =========================================================

depth_results <- imap_dfr(
  outs,
  function(outdir, run_name) {
    
    message("Processing run: ", run_name)
    
    mol_list <- read_molecule_info_for_filtered_cells(outdir)
    
    simulate_depth_from_molecules(
      mol = mol_list$mol,
      filtered_barcodes = mol_list$filtered_barcodes,
      fractions = fractions
    ) %>%
      mutate(dataset = run_name)
  }
)

write_csv(
  depth_results,
  file.path(output_dir, "cellranger_molecule_downsampling_per_cell.csv")
)

depth_summary <- summarise_qc(depth_results)

write_csv(
  depth_summary,
  file.path(output_dir, "cellranger_molecule_downsampling_summary.csv")
)

# =========================================================
# Add actual full-depth QC check
# =========================================================

full_depth_check <- own_full_qc %>%
  group_by(dataset) %>%
  summarise(
    cells = n(),
    median_genes_matrix = median(nGene),
    median_UMI_matrix = median(nUMI),
    mean_genes_matrix = mean(nGene),
    mean_UMI_matrix = mean(nUMI),
    .groups = "drop"
  ) %>%
  left_join(
    depth_summary %>%
      filter(fraction == 1) %>%
      select(
        dataset,
        median_genes_molecule = median_genes,
        median_UMI_molecule = median_UMI,
        mean_genes_molecule = mean_genes,
        mean_UMI_molecule = mean_UMI
      ),
    by = "dataset"
  )

write_csv(
  full_depth_check,
  file.path(output_dir, "full_depth_matrix_vs_molecule_check.csv")
)

# =========================================================
# Cell Ranger metrics
# =========================================================

cr_metrics <- imap_dfr(
  outs,
  ~ read_cellranger_metrics(.x, .y)
)

write_csv(
  cr_metrics,
  file.path(output_dir, "cellranger_metrics_summary_combined.csv")
)

# =========================================================
# Benchmark summaries against AML atlas
# =========================================================

comparison_full_qc <- bind_rows(
  own_full_qc %>%
    select(dataset, barcode, fraction, nGene, nUMI),
  aml_qc %>%
    select(dataset, barcode, fraction, nGene, nUMI)
)

benchmark_summary <- summarise_qc(comparison_full_qc)

write_csv(
  comparison_full_qc,
  file.path(output_dir, "own_full_depth_vs_AML_scAtlas_per_cell_QC.csv")
)

write_csv(
  benchmark_summary,
  file.path(output_dir, "own_full_depth_vs_AML_scAtlas_summary.csv")
)

# =========================================================
# AML benchmark lines
# =========================================================

aml_summary <- benchmark_summary %>%
  filter(dataset == "AML scAtlas", fraction == 1)

aml_median_genes <- aml_summary$median_genes
aml_median_UMI <- aml_summary$median_UMI

# =========================================================
# Plotting
# =========================================================

depth_results <- depth_results %>%
  mutate(
    dataset = factor(dataset, levels = c("first_run", "resequenced")),
    fraction = as.numeric(fraction)
  )

depth_summary <- depth_summary %>%
  mutate(
    dataset = factor(dataset, levels = c("first_run", "resequenced")),
    fraction = as.numeric(fraction)
  )

comparison_full_qc <- comparison_full_qc %>%
  mutate(
    dataset = factor(dataset, levels = c("AML scAtlas", "first_run", "resequenced"))
  )

plot_ylim_gene <- quantile(comparison_full_qc$nGene, 0.99, na.rm = TRUE)
plot_ylim_umi <- quantile(comparison_full_qc$nUMI, 0.99, na.rm = TRUE)

# ---------------------------------------------------------
# Downsampling curves: median genes / UMIs
# ---------------------------------------------------------

p_downsample_genes <- depth_summary %>%
  ggplot(aes(x = fraction, y = median_genes, colour = dataset)) +
  geom_hline(
    yintercept = aml_median_genes,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    x = "Fraction of current sequencing depth",
    y = "Median genes per cell",
    title = "Downsampled gene detection vs AML scAtlas median"
  )

p_downsample_umi <- depth_summary %>%
  ggplot(aes(x = fraction, y = median_UMI, colour = dataset)) +
  geom_hline(
    yintercept = aml_median_UMI,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    x = "Fraction of current sequencing depth",
    y = "Median UMIs per cell",
    title = "Downsampled UMI detection vs AML scAtlas median"
  )

# ---------------------------------------------------------
# Distributions at selected depths
# ---------------------------------------------------------

selected_fractions <- c(0.10, 0.30, 0.50, 1.00)

depth_selected <- depth_results %>%
  filter(fraction %in% selected_fractions) %>%
  mutate(
    fraction_label = paste0(fraction * 100, "%")
  )

p_depth_genes_dist <- depth_selected %>%
  ggplot(aes(x = fraction_label, y = nGene, fill = dataset)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Downsampled depth",
    y = "Genes per cell",
    title = "Genes per cell at selected sequencing depths"
  )

p_depth_umi_dist <- depth_selected %>%
  ggplot(aes(x = fraction_label, y = nUMI, fill = dataset)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Downsampled depth",
    y = "UMIs per cell",
    title = "UMIs per cell at selected sequencing depths"
  )

# ---------------------------------------------------------
# Full-depth comparison to AML atlas
# ---------------------------------------------------------

p_compare_genes <- comparison_full_qc %>%
  ggplot(aes(x = dataset, y = nGene, fill = dataset)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, plot_ylim_gene)) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Genes per cell",
    title = "Full-depth genes per cell vs AML scAtlas"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_compare_umi <- comparison_full_qc %>%
  ggplot(aes(x = dataset, y = nUMI, fill = dataset)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, plot_ylim_umi)) +
  theme_bw() +
  labs(
    x = NULL,
    y = "UMIs per cell",
    title = "Full-depth UMIs per cell vs AML scAtlas"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# =========================================================
# Save plots
# =========================================================

ggsave(
  file.path(output_dir, "downsampled_median_genes_vs_AML_atlas.pdf"),
  p_downsample_genes,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_dir, "downsampled_median_UMIs_vs_AML_atlas.pdf"),
  p_downsample_umi,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_dir, "downsampled_genes_selected_depths.pdf"),
  p_depth_genes_dist,
  width = 9,
  height = 5
)

ggsave(
  file.path(output_dir, "downsampled_UMIs_selected_depths.pdf"),
  p_depth_umi_dist,
  width = 9,
  height = 5
)

ggsave(
  file.path(output_dir, "full_depth_genes_vs_AML_atlas.pdf"),
  p_compare_genes,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "full_depth_UMIs_vs_AML_atlas.pdf"),
  p_compare_umi,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "combined_downsampling_curves.pdf"),
  p_downsample_genes / p_downsample_umi,
  width = 8,
  height = 9
)

ggsave(
  file.path(output_dir, "combined_full_depth_vs_AML_atlas.pdf"),
  p_compare_genes / p_compare_umi,
  width = 8,
  height = 9
)

# =========================================================
# Print summaries
# =========================================================

print(depth_summary)
print(full_depth_check)
print(benchmark_summary)

message("Analysis complete.")