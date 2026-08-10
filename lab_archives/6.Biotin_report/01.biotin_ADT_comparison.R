#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Biotin report - step 01 of 4
#
# Biotin ADT levels compared across the LK1 and LK2 runs.
#
# Split out of scripts/6b.QC_HTO_comparison.R so the biotin part can run on its
# own without the HTO/demux steps. Same run list and output directory as the
# original; 6b itself belongs to the demultiplexing benchmark report.
#
# Frozen for the lab archive 2026-08-04 from scripts/6b1.biotin_ADT_comparison.R (mtime 2026-07-22).
# md5 of the original: a14656892d0bf6a7c31b5c96d1a827d8
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# Biotin ADT comparison across runs - extracted from 6b.QC_HTO_comparison.R
# (which does a broader HTO+ADT comparison; this pulls out just the biotin
# piece so it can be run on its own without the HTO/demux steps). Same
# run list, same output directory/filenames as the original.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
})

# ============================================================
# Compare biotin ADT levels across LK1 vs LK2 runs
# ============================================================

run_tbl <- tibble::tribble(
  ~run,                              ~project, ~label,
  "260522_VH01624_461_222JLJVNX",    "LK1",    "LK1",
  "260528_VH01624_464_222K7VKNX",    "LK2",    "LK2",
  "260717_VH01624_477_222KG22NX",    "LK3",    "LK3"
)

out_dir <- file.path(
  "../results/seurat_demux_test",
  "LK1_vs_LK2_test_HTO_ADT_comparison"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

process_one_run_biotin <- function(run, project, label) {

  message("Processing ", label, " | ", run, " | ", project)

  sample_name <- paste0(project, "-GEX")

  doublet_dir <- file.path(
    "../results_nf", run, "rds", "05_scdblfinder", sample_name, "rds"
  )

  adt_dir <- file.path(
    "../results_nf", run, "04_cite_seq_count", sample_name,
    paste0("adt_counts_", project, "_raw"),
    "umi_count"
  )

  seu <- readRDS(
    file.path(doublet_dir, "seurat_emptyDrops_RNA_scDblFinder_filtered.rds")
  )

  adt_mat <- read_csc_matrix(adt_dir, add_suffix = TRUE)

  rownames(adt_mat) <- clean_adt_name(rownames(adt_mat))
  rownames(adt_mat) <- make.unique(rownames(adt_mat))

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
      biotin_match = adt,
      n_gex_cells = ncol(seu)
    )

  biotin_summary
}

res <- pmap(
  run_tbl,
  process_one_run_biotin
)

biotin_summary <- bind_rows(res)

# ----------------------------
# Save tables
# ----------------------------

write_csv(biotin_summary, file.path(out_dir, "test_biotin_read_rank_summary.csv"))

biotin_compare <- biotin_summary %>%
  select(label, adt, total_reads_shared_cells, rank_shared_cells) %>%
  pivot_wider(
    names_from = label,
    values_from = c(total_reads_shared_cells, rank_shared_cells)
  )

write_csv(biotin_compare, file.path(out_dir, "test_biotin_before_after_comparison.csv"))

# ----------------------------
# Plots
# ----------------------------

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

ggsave(
  file.path(out_dir, "test_biotin_reads_before_after.pdf"),
  p_biotin,
  width = 8,
  height = 5
)

biotin_summary <- biotin_summary %>%
  mutate(
    reads_per_cell_shared = total_reads_shared_cells / n_gex_cells
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

ggsave(
  file.path(out_dir, "test_biotin_reads_per_cell_before_after.pdf"),
  p_biotin_per_cell,
  width = 8,
  height = 5
)

cat("\nDone.\n")
cat("Output written to:\n", out_dir, "\n\n")

cat("Biotin summary:\n")
print(biotin_summary)
VlnPlot(seu, features = dsb_markers, assay = "CITE_DSB", slot = "data",
        group.by = "predicted_CellType_Broad", split.by = "sample_name",
        pt.size = 0, cols=c("gray50","gray70","firebrick")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
