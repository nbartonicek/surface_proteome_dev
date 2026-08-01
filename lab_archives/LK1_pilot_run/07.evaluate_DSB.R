#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 07 of 18
#
# DSB normalisation using raw background droplets, with random IgG-control subsampling.
#
# Frozen for the lab archive 2026-07-31 from scripts/7c.evaluate_DSB.R (mtime 2026-05-22).
# md5 of the original: 58a0f44be40658ffc50cc0b542e73267
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(dsb)
})

set.seed(1)

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"
sample_short <- "LK1"
patient_keep <- "HBDN206-MNpCT"

annotation_dir <- file.path("../results/seurat_annotated", run)

seurat_file <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)
adt_cells_dir <- file.path(
  "../results/cite_seq_count", run,
  paste0("adt_counts_", sample_short, "_emptydrops/umi_count")
)

adt_raw_dir <- file.path(
  "../results/cite_seq_count", run,
  paste0("adt_counts_", sample_short, "/umi_count")
)

out_dir <- file.path(
  "../results/dsb_igg_subsample", run, sample_name, patient_keep
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_cite <- function(d) {
  m <- readMM(file.path(d, "matrix.mtx.gz"))
  
  bc <- fread(file.path(d, "barcodes.tsv.gz"), header = FALSE)$V1
  ft <- fread(file.path(d, "features.tsv.gz"), header = FALSE)$V1
  
  ft <- gsub("-[ACGT]+$", "", ft)
  ft <- make.unique(ft)
  
  if (!any(grepl("-1$", bc))) {
    bc <- paste0(bc, "-1")
  }
  
  rownames(m) <- ft
  colnames(m) <- bc
  
  m[!grepl("^unmapped$", rownames(m), ignore.case = TRUE), , drop = FALSE]
}

message("Loading Seurat object...")
seu <- readRDS(seurat_file)

stopifnot("sample_name" %in% colnames(seu@meta.data))

seu <- subset(seu, subset = sample_name == patient_keep)

message("Cells in ", patient_keep, ": ", ncol(seu))

message("Loading ADT matrices...")
adt_cells <- read_cite(adt_cells_dir)
adt_raw <- read_cite(adt_raw_dir)

common_cells <- intersect(colnames(seu), colnames(adt_cells))
stopifnot(length(common_cells) > 100)

seu <- subset(seu, cells = common_cells)
adt_cells <- adt_cells[, colnames(seu), drop = FALSE]

common_proteins <- intersect(rownames(adt_cells), rownames(adt_raw))
adt_cells <- adt_cells[common_proteins, , drop = FALSE]
adt_raw <- adt_raw[common_proteins, , drop = FALSE]

message("Cell ADT matrix: ", nrow(adt_cells), " proteins x ", ncol(adt_cells), " cells")
message("Raw ADT matrix: ", nrow(adt_raw), " proteins x ", ncol(adt_raw), " droplets")

# ----------------------------
# Background droplets
# ----------------------------

cell_barcodes <- colnames(adt_cells)
bg_barcodes <- setdiff(colnames(adt_raw), cell_barcodes)

stopifnot(length(bg_barcodes) > 100)

bg_sizes <- Matrix::colSums(adt_raw[, bg_barcodes, drop = FALSE])

bg_barcodes <- names(bg_sizes)[
  bg_sizes > quantile(bg_sizes, 0.50, na.rm = TRUE) &
    bg_sizes < quantile(bg_sizes, 0.99, na.rm = TRUE)
]

bg_barcodes <- head(
  bg_barcodes[order(bg_sizes[bg_barcodes], decreasing = TRUE)],
  5000
)

bg <- adt_raw[, bg_barcodes, drop = FALSE]

message("Background matrix: ", nrow(bg), " proteins x ", ncol(bg), " droplets")

# ----------------------------
# Filter proteins
# ----------------------------

keep <- apply(adt_cells, 1, max) > 5 &
  Matrix::rowSums(adt_cells) > 0 &
  Matrix::rowSums(bg) > 0

adt_cells <- adt_cells[keep, , drop = FALSE]
bg <- bg[rownames(adt_cells), , drop = FALSE]

igg <- grep(
  "IgG|Rat-IgG|Hamster-IgG|Armenian-Hamster-IgG|isotype",
  rownames(adt_cells),
  ignore.case = TRUE,
  value = TRUE
)

igg <- setdiff(
  igg,
  grep(
    "IgM|IgD|IgE|IgG-Fc|light-chain",
    igg,
    ignore.case = TRUE,
    value = TRUE
  )
)

message("After filtering:")
message("  Cell ADT: ", nrow(adt_cells), " proteins x ", ncol(adt_cells), " cells")
message("  Background: ", nrow(bg), " proteins x ", ncol(bg), " droplets")
message("  IgG/isotype controls: ", paste(igg, collapse = ", "))

stopifnot(nrow(adt_cells) > 10)
stopifnot(ncol(adt_cells) > 100)
stopifnot(ncol(bg) > 100)
stopifnot(length(igg) >= 1)

write_csv(
  tibble(igg_control = igg),
  file.path(out_dir, "igg_controls_detected.csv")
)

# ----------------------------
# DSB helper
# ----------------------------

run_dsb <- function(ctrls) {
  use_iso <- length(ctrls) > 0
  
  DSBNormalizeProtein(
    cell_protein_matrix = as.matrix(adt_cells),
    empty_drop_matrix = as.matrix(bg),
    denoise.counts = use_iso,
    use.isotype.control = use_iso,
    isotype.control.name.vec = if (use_iso) ctrls else NULL,
    quantile.clipping = FALSE
  )
}

# ----------------------------
# Baseline: ambient correction only
# ----------------------------

message("Running baseline DSB without isotype denoising...")
dsb0 <- run_dsb(character(0))

baseline <- median(abs(dsb0[igg, , drop = FALSE]), na.rm = TRUE)

message("Baseline median absolute IgG signal: ", round(baseline, 4))

# ----------------------------
# Random IgG-control subsampling
# ----------------------------

n_rep <- 30
max_k <- length(igg)

res <- map_dfr(0:max_k, function(k) {
  map_dfr(seq_len(n_rep), function(r) {
    ctrls <- if (k == 0) character(0) else sample(igg, k)
    
    dsb_mat <- run_dsb(ctrls)
    
    median_abs <- median(abs(dsb_mat[igg, , drop = FALSE]), na.rm = TRUE)
    mean_abs <- mean(abs(dsb_mat[igg, , drop = FALSE]), na.rm = TRUE)
    
    tibble(
      n_igg_controls_used = k,
      replicate = r,
      controls_used = paste(ctrls, collapse = ";"),
      median_abs_igg_dsb = median_abs,
      mean_abs_igg_dsb = mean_abs,
      pct_reduction_vs_no_isotype = 100 * (baseline - median_abs) / baseline
    )
  })
})

write_csv(
  res,
  file.path(out_dir, "dsb_igg_control_subsampling.csv")
)

saveRDS(
  res,
  file.path(out_dir, "dsb_igg_control_subsampling.rds")
)

# ----------------------------
# Plot
# ----------------------------

p <- ggplot(
  res,
  aes(
    x = factor(n_igg_controls_used),
    y = pct_reduction_vs_no_isotype
  )
) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  theme_bw() +
  labs(
    x = "Number of IgG/isotype controls supplied to DSB",
    y = "% reduction in median absolute IgG signal",
    title = "DSB reduction of IgG/isotype control signal",
    subtitle = paste0(
      "Patient subset: ", patient_keep,
      "; relative to ambient-only DSB without isotype denoising"
    )
  )

ggsave(
  file.path(out_dir, "dsb_igg_reduction_subsampling.pdf"),
  p,
  width = 7,
  height = 5
)

ggsave(
  file.path(out_dir, "dsb_igg_reduction_subsampling.png"),
  p,
  width = 7,
  height = 5,
  dpi = 300
)

message("Done.")
message("Output written to: ", out_dir)

