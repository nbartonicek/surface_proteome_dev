#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Demultiplexing benchmark: genotyping vs HTO - step 05 of 7
#
# Compares HTO and ADT behaviour between the LK1 and LK2 runs.
#
# The biotin-specific part of this script was later split out into
# scripts/6b1.biotin_ADT_comparison.R, which belongs to the biotin report.
#
# Frozen for the lab archive 2026-08-03 from scripts/6b.QC_HTO_comparison.R (mtime 2026-06-03).
# md5 of the original: ca9bc9dca7bb0247a98da91ecb419acb
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# Compare LK1 vs LK2 HTO/ADT runs
# ============================================================

run_tbl <- tibble::tribble(
  ~run,                              ~project, ~label,
  "260522_VH01624_461_222JLJVNX",    "LK1",    "LK1",
  "260528_VH01624_464_222K7VKNX",    "LK2",    "LK2"
)

out_dir <- file.path(
  "../results/seurat_demux_test",
  "LK1_vs_LK2_test_HTO_ADT_comparison"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hto_positive_quantile <- 0.99

hto_names <- c(
  "HTO1-GTCAACTCTTTAGCG" = "HTO1",
  "HTO2-TGATGGCCTATTGGG" = "HTO2",
  "HTO3-TTCCGCCTCTCTTTG" = "HTO3",
  "HTO4-AGTAAGTTCAGCGTA" = "HTO4"
)

clean_adt_name <- function(x) {
  x %>%
    gsub("-[ACGT]+$", "", .) %>%
    gsub("_TotalSeq.*$", "", .)
}

read_csc_matrix <- function(dir, add_suffix = TRUE) {
  mat <- Matrix::readMM(file.path(dir, "matrix.mtx.gz"))
  
  barcodes <- data.table::fread(
    file.path(dir, "barcodes.tsv.gz"),
    header = FALSE
  )$V1
  
  features <- data.table::fread(
    file.path(dir, "features.tsv.gz"),
    header = FALSE
  )
  
  feature_names <- if (ncol(features) >= 2) features$V2 else features$V1
  feature_names <- make.unique(feature_names)
  
  if (add_suffix) {
    barcodes <- paste0(barcodes, "-1")
  }
  
  rownames(mat) <- feature_names
  colnames(mat) <- barcodes
  
  mat
}

process_one_run <- function(run, project, label) {
  
  message("Processing ", label, " | ", run, " | ", project)
  
  sample_name <- paste0(project, "-GEX")
  
  doublet_dir <- file.path(
    "../results/scDblFinder",
    run,
    sample_name
  )
  
  hto_dir <- file.path(
    "../results/cite_seq_count",
    run,
    paste0("hto_counts_", project, "_emptydrops"),
    "umi_count"
  )
  
  adt_dir <- file.path(
    "../results/cite_seq_count_fromRaw",
    run,
    paste0("adt_counts_", project, "_raw"),
    "umi_count"
  )
  
  seu <- readRDS(
    file.path(doublet_dir, "seurat_emptyDrops_RNA_scDblFinder_filtered.rds")
  )
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  # ----------------------------
  # HTO
  # ----------------------------
  
  hto_mat <- read_csc_matrix(hto_dir, add_suffix = TRUE)
  
  hto_mat <- hto_mat[
    !grepl("^unmapped$", rownames(hto_mat), ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  
  hto_raw_feature <- rownames(hto_mat)
  
  new_hto_names <- hto_names[rownames(hto_mat)]
  rownames(hto_mat) <- ifelse(
    is.na(new_hto_names),
    rownames(hto_mat),
    new_hto_names
  )
  
  hto_total_reads <- sum(hto_mat)
  
  common_hto <- intersect(colnames(seu), colnames(hto_mat))
  message("Number of cells",cat(length(common_hto)))
  seu_hto <- subset(seu, cells = common_hto)
  hto_mat <- hto_mat[, common_hto, drop = FALSE]
  
  seu_hto[["HTO"]] <- CreateAssayObject(counts = hto_mat)
  hto_total <- Matrix::colSums(hto_mat)
  
  DefaultAssay(seu_hto) <- "HTO"
  seu_hto <- subset(
    seu_hto,
    cells = names(hto_total)[hto_total > 0]
  )
  seu_hto <- NormalizeData(
    seu_hto,
    assay = "HTO",
    normalization.method = "CLR",
    margin = 2,
    verbose = FALSE
  )
  
  seu_hto <- HTODemux(
    seu_hto,
    assay = "HTO",
    positive.quantile = hto_positive_quantile
  )
  
  hto_counts <- GetAssayData(seu_hto, assay = "HTO", layer = "counts")
  hto_clr <- GetAssayData(seu_hto, assay = "HTO", layer = "data")
  
  hto_feature_summary <- tibble(
    run = run,
    project = project,
    label = label,
    hto = rownames(hto_counts),
    total_reads = as.numeric(Matrix::rowSums(hto_counts)),
    rank_total_reads = rank(-as.numeric(Matrix::rowSums(hto_counts)), ties.method = "first"),
    clr_global_p99_cutoff = apply(
      as.matrix(hto_clr),
      1,
      quantile,
      probs = hto_positive_quantile,
      na.rm = TRUE
    )
  ) %>%
    arrange(rank_total_reads)
  
  hto_demux_summary <- seu_hto@meta.data %>%
    as_tibble(rownames = "barcode") %>%
    dplyr::count(HTO_classification.global, hash.ID, name = "n_cells") %>%
    mutate(
      run = run,
      project = project,
      label = label,
      .before = 1
    )
  
  hto_run_summary <- seu_hto@meta.data %>%
    as_tibble(rownames = "barcode") %>%
    summarise(
      run = run,
      project = project,
      label = label,
      n_gex_cells = ncol(seu),
      n_hto_barcodes = ncol(hto_mat),
      n_shared_gex_hto_cells = n(),
      total_hto_reads = hto_total_reads,
      n_hto_singlets = sum(HTO_classification.global == "Singlet", na.rm = TRUE),
      pct_hto_singlets = 100 * mean(HTO_classification.global == "Singlet", na.rm = TRUE),
      n_hto_doublets = sum(HTO_classification.global == "Doublet", na.rm = TRUE),
      pct_hto_doublets = 100 * mean(HTO_classification.global == "Doublet", na.rm = TRUE),
      n_hto_negative = sum(HTO_classification.global == "Negative", na.rm = TRUE),
      pct_hto_negative = 100 * mean(HTO_classification.global == "Negative", na.rm = TRUE)
    )
  
  # ----------------------------
  # ADT
  # ----------------------------
  
  adt_mat <- read_csc_matrix(adt_dir, add_suffix = TRUE)
  
  rownames(adt_mat) <- clean_adt_name(rownames(adt_mat))
  rownames(adt_mat) <- make.unique(rownames(adt_mat))
  
  adt_total_reads <- sum(adt_mat)
  
  common_adt <- intersect(colnames(seu), colnames(adt_mat))
  
  adt_mat_shared <- adt_mat[, common_adt, drop = FALSE]
  
  adt_feature_counts <- Matrix::rowSums(adt_mat_shared)
  
  adt_feature_summary <- tibble(
    run = run,
    project = project,
    label = label,
    adt = names(adt_feature_counts),
    total_reads_shared_cells = as.numeric(adt_feature_counts),
    rank_shared_cells = rank(-as.numeric(adt_feature_counts), ties.method = "first"),
    is_biotin = grepl("biotin", names(adt_feature_counts), ignore.case = TRUE)
  ) %>%
    arrange(rank_shared_cells)
  
  biotin_summary <- adt_feature_summary %>%
    filter(is_biotin) %>%
    mutate(
      biotin_match = adt
    )
  
  adt_run_summary <- tibble(
    run = run,
    project = project,
    label = label,
    n_gex_cells = ncol(seu),
    n_adt_barcodes = ncol(adt_mat),
    n_shared_gex_adt_cells = length(common_adt),
    total_adt_reads_all_barcodes = adt_total_reads,
    total_adt_reads_shared_cells = sum(adt_mat_shared)
  )
  
  list(
    hto_run_summary = hto_run_summary,
    hto_demux_summary = hto_demux_summary,
    hto_feature_summary = hto_feature_summary,
    adt_run_summary = adt_run_summary,
    adt_feature_summary = adt_feature_summary,
    biotin_summary = biotin_summary
  )
}

res <- pmap(
  run_tbl,
  process_one_run
)

hto_run_summary <- map_dfr(res, "hto_run_summary")
hto_demux_summary <- map_dfr(res, "hto_demux_summary")
hto_feature_summary <- map_dfr(res, "hto_feature_summary")
adt_run_summary <- map_dfr(res, "adt_run_summary")
adt_feature_summary <- map_dfr(res, "adt_feature_summary")
biotin_summary <- map_dfr(res, "biotin_summary")

# ----------------------------
# Save tables
# ----------------------------

write_csv(hto_run_summary, file.path(out_dir, "test_hto_run_summary.csv"))
write_csv(hto_demux_summary, file.path(out_dir, "test_hto_demux_summary.csv"))
write_csv(hto_feature_summary, file.path(out_dir, "test_hto_feature_cutoffs_and_reads.csv"))

write_csv(adt_run_summary, file.path(out_dir, "test_adt_run_summary.csv"))
write_csv(adt_feature_summary, file.path(out_dir, "test_adt_feature_read_ranks.csv"))
write_csv(biotin_summary, file.path(out_dir, "test_biotin_read_rank_summary.csv"))

# ----------------------------
# Before/after comparison tables
# ----------------------------

hto_compare <- hto_run_summary %>%
  select(label, total_hto_reads, n_shared_gex_hto_cells, pct_hto_singlets,
         pct_hto_doublets, pct_hto_negative) %>%
  pivot_longer(-label, names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = label, values_from = value)

adt_compare <- adt_run_summary %>%
  select(label, total_adt_reads_all_barcodes, total_adt_reads_shared_cells,
         n_shared_gex_adt_cells) %>%
  pivot_longer(-label, names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = label, values_from = value)

biotin_compare <- biotin_summary %>%
  select(label, adt, total_reads_shared_cells, rank_shared_cells) %>%
  pivot_wider(
    names_from = label,
    values_from = c(total_reads_shared_cells, rank_shared_cells)
  )

write_csv(hto_compare, file.path(out_dir, "test_hto_before_after_comparison.csv"))
write_csv(adt_compare, file.path(out_dir, "test_adt_before_after_comparison.csv"))
write_csv(biotin_compare, file.path(out_dir, "test_biotin_before_after_comparison.csv"))

# ----------------------------
# Plots
# ----------------------------

p_hto_reads <- hto_run_summary %>%
  ggplot(aes(x = label, y = total_hto_reads)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Total HTO reads",
    x = NULL,
    y = "Total reads"
  )

ggsave(
  file.path(out_dir, "test_total_hto_reads.pdf"),
  p_hto_reads,
  width = 6,
  height = 5
)

hto_run_summary <- hto_run_summary %>%
  mutate(
    hto_reads_per_cell = total_hto_reads / n_gex_cells
  )

p_hto_reads_per_cell <- hto_run_summary %>%
  ggplot(aes(x = label, y = hto_reads_per_cell)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "HTO reads per cell",
    x = NULL,
    y = "Reads per cell"
  )

p_hto_reads_per_cell

p_hto_singlets <- hto_run_summary %>%
  ggplot(aes(x = label, y = pct_hto_singlets)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Percentage assigned as HTO singlets",
    x = NULL,
    y = "% HTO singlets"
  )

ggsave(
  file.path(out_dir, "test_hto_percent_singlets.pdf"),
  p_hto_singlets,
  width = 6,
  height = 5
)

p_adt_reads <- adt_run_summary %>%
  ggplot(aes(x = label, y = total_adt_reads_shared_cells)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Total ADT reads in shared GEX/ADT cells",
    x = NULL,
    y = "Total reads"
  )

ggsave(
  file.path(out_dir, "test_total_adt_reads_shared_cells.pdf"),
  p_adt_reads,
  width = 6,
  height = 5
)

adt_run_summary <- adt_run_summary %>%
  mutate(
    adt_reads_per_cell = total_adt_reads_shared_cells / n_gex_cells
  )

p_adt_reads_per_cell <- adt_run_summary %>%
  ggplot(aes(x = label, y = adt_reads_per_cell)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "ADT reads per cell",
    x = NULL,
    y = "Reads per cell"
  )

p_adt_reads_per_cell

p_biotin <- biotin_summary %>%
  ggplot(aes(x = label, y = total_reads_shared_cells, fill = adt)) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 12) +
  labs(
    title = "Biotin ADT reads before vs after",
    x = NULL,
    y = "Total reads in shared cells",
    fill = "Biotin feature"
  )
p_biotin



ggsave(
  file.path(out_dir, "test_biotin_reads_before_after.pdf"),
  p_biotin,
  width = 8,
  height = 5
)

biotin_summary <- biotin_summary %>%
  left_join(
    adt_run_summary %>%
      select(label, n_gex_cells),
    by = "label"
  ) %>%
  mutate(
    reads_per_cell_shared =
      total_reads_shared_cells / n_gex_cells
  )

p_biotin_per_cell <- biotin_summary %>%
  ggplot(
    aes(
      x = label,
      y = reads_per_cell_shared,
      fill = adt
    )
  ) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 12) +
  labs(
    title = "Biotin ADT reads per cell before vs after",
    x = NULL,
    y = "Reads per shared cell",
    fill = "Biotin feature"
  )

p_biotin_per_cell

p_adt_rank <- adt_feature_summary %>%
  mutate(
    adt = fct_reorder(adt, total_reads_shared_cells, .desc = TRUE)
  ) %>%
  ggplot(aes(x = adt, y = total_reads_shared_cells, fill = is_biotin)) +
  geom_col() +
  facet_wrap(~ label, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
  ) +
  labs(
    title = "ADT feature read counts and biotin rank",
    x = "ADT feature",
    y = "Total reads in shared cells",
    fill = "Biotin"
  )

ggsave(
  file.path(out_dir, "test_adt_feature_reads_ranked.pdf"),
  p_adt_rank,
  width = 18,
  height = 7
)

cat("\nDone.\n")
cat("Output written to:\n", out_dir, "\n\n")

cat("HTO summary:\n")
print(hto_run_summary)

cat("\nADT summary:\n")
print(adt_run_summary)

cat("\nBiotin summary:\n")
print(biotin_summary)