#!/usr/bin/env Rscript

# AML-vs-normal-01 pseudobulk DE, generalized across every broad cell type
# present in both AML samples and normal-01, RESTRICTED to the union of all
# AML patients' candidate surface-gene lists from
# 24.leucegene_reference_and_candidate_genes.R - not genome-wide. Union
# (rather than a single patient's list) because this comparison pools
# multiple AML patients together per cell type, and different patients may
# have been assigned to different Leucegene clusters; a gene flagged as
# subtype-relevant for any one of them is still worth testing here. See
# 24.leucegene_reference_and_candidate_genes.R's header comment for why
# pre-restricting the candidate list before testing matters.
#
# Reuses 17.differential_expression.rmd's exact patterns: loading + merging
# both runs' annotated Seurat objects, make_pseudobulk_from_cells(),
# fixed-dispersion edgeR exactTest (BCV = 0.4 - still appropriate, one
# pseudobulk sample per group per cell type), the >=50-cells-per-sample
# filter, and excluding MOLM13 (cell line control).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
})

# ----------------------------
# Paths / settings
# ----------------------------

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

runs <- c(
  "260423_VH01624_453_222HWMYNX",
  "260528_VH01624_464_222K7VKNX"
)

candidate_dir <- file.path(proj, "results/leucegene_reference/candidate_genes_by_patient")

out_dir <- file.path(proj, "results/differential_expression/AML_vs_normal01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normal_sample  <- "normal-01"
min_cells      <- 50
bcv            <- 0.4
fdr_cutoff     <- 0.05
logfc_cutoff   <- 1

# ----------------------------
# Helpers (same as 17.differential_expression.rmd)
# ----------------------------

make_pseudobulk_from_cells <- function(count_mat, cells) {
  cells <- intersect(cells, colnames(count_mat))
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}

run_edgeR_exact <- function(count_mat, group_vec, bcv = 0.4) {
  group_vec <- factor(group_vec, levels = c("normal", "aml"))

  dge <- edgeR::DGEList(counts = count_mat, group = group_vec)
  keep <- edgeR::filterByExpr(dge, group = group_vec)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  et <- edgeR::exactTest(dge, pair = c("normal", "aml"), dispersion = bcv^2)

  edgeR::topTags(et, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    as_tibble()
}

# ----------------------------
# 0. Union of candidate surface genes across all AML patients (excludes
# normal-01/MOLM13 candidate files, if any happen to exist)
# ----------------------------

candidate_files <- list.files(candidate_dir, pattern = "_candidate_surface_genes\\.csv$", full.names = TRUE)
candidate_files <- candidate_files[!grepl("normal-01|MOLM13", candidate_files, ignore.case = TRUE)]

if (length(candidate_files) == 0) {
  stop(
    "No candidate gene list files found in ", candidate_dir, " - ",
    "run 24.leucegene_reference_and_candidate_genes.R first."
  )
}

candidate_genes <- candidate_files %>%
  map(~ read_csv(.x, show_col_types = FALSE)) %>%
  bind_rows() %>%
  pull(gene) %>%
  unique()

message("Union of candidate surface genes across ", length(candidate_files), " AML patients: ", length(candidate_genes))

if (length(candidate_genes) == 0) {
  stop("Candidate gene union is empty - nothing to test.")
}

# ----------------------------
# 1. Load + merge both runs (same pattern as 17.differential_expression.rmd)
# ----------------------------

seuList <- list()
for (run in runs) {
  seurat_file <- file.path(proj, "results/seurat_annotated", run, "demux_singlets_annotated_seurat.rds")
  stopifnot(file.exists(seurat_file))
  seuList[[run]] <- readRDS(seurat_file)
}

seu <- merge(x = seuList[[1]], y = seuList[[2]], add.cell.ids = c("run1", "run2"), project = "combined")
rm(seuList)

DefaultAssay(seu) <- "RNA"
if (inherits(seu[["RNA"]], "Assay5")) {
  seu <- JoinLayers(seu, assay = "RNA")
}

stopifnot(all(c("sample_name", "predicted_CellType_Broad", "scDblFinder.class") %in% colnames(seu@meta.data)))

# ----------------------------
# Restrict to singlets, drop cell-line controls
# ----------------------------

seu <- subset(
  seu,
  subset = scDblFinder.class == "singlet" & !grepl("MOLM13", sample_name, ignore.case = TRUE)
)

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")

candidate_genes_present <- intersect(candidate_genes, rownames(counts))
message("Candidate genes present in RNA assay: ", length(candidate_genes_present), " / ", length(candidate_genes))

if (length(candidate_genes_present) < 3) {
  stop("Too few candidate genes present in the RNA assay to test.")
}

celltypes <- seu@meta.data %>%
  filter(!is.na(predicted_CellType_Broad)) %>%
  pull(predicted_CellType_Broad) %>%
  unique() %>%
  sort()

message("Broad cell types to test: ", paste(celltypes, collapse = ", "))

all_de_results <- list()
kept_samples_summary <- list()

for (ct in celltypes) {

  message("\n=== ", ct, " ===")

  cell_counts <- seu@meta.data %>%
    filter(predicted_CellType_Broad == ct) %>%
    count(sample_name, name = "n_cells")

  good_samples <- cell_counts %>% filter(n_cells >= min_cells) %>% pull(sample_name)

  has_normal <- normal_sample %in% good_samples
  aml_samples <- setdiff(good_samples, normal_sample)

  kept_samples_summary[[ct]] <- cell_counts %>%
    mutate(celltype = ct, kept = sample_name %in% good_samples)

  if (!has_normal || length(aml_samples) == 0) {
    message("  Skipping ", ct, ": needs normal-01 plus >=1 AML sample with >=", min_cells, " cells (has normal-01: ",
            has_normal, ", AML samples passing filter: ", length(aml_samples), ")")
    next
  }

  meta_ct <- seu@meta.data %>%
    rownames_to_column("cell") %>%
    filter(predicted_CellType_Broad == ct, sample_name %in% good_samples) %>%
    mutate(group = if_else(sample_name == normal_sample, "normal", "aml"))

  sample_ids <- unique(meta_ct$sample_name)

  pb_list <- lapply(sample_ids, function(sid) {
    cells_i <- meta_ct %>% filter(sample_name == sid) %>% pull(cell)
    pb <- make_pseudobulk_from_cells(counts, cells_i)
    names(pb) <- rownames(counts)
    pb
  })
  names(pb_list) <- sample_ids

  pb_counts <- do.call(cbind, pb_list)
  rownames(pb_counts) <- rownames(counts)

  # Restrict to the candidate surface-gene union before testing
  pb_counts_restricted <- pb_counts[candidate_genes_present, , drop = FALSE]

  pb_meta <- tibble(sample = colnames(pb_counts_restricted)) %>%
    left_join(meta_ct %>% distinct(sample_name, group) %>% rename(sample = sample_name), by = "sample")

  de_res <- run_edgeR_exact(pb_counts_restricted, group_vec = pb_meta$group, bcv = bcv) %>%
    mutate(
      celltype = ct,
      direction = case_when(
        FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ paste0("Up in AML ", ct),
        FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ paste0("Down in AML ", ct),
        TRUE ~ "Not significant"
      )
    )

  write_csv(
    de_res,
    file.path(out_dir, paste0(gsub("[^A-Za-z0-9]+", "_", ct), "_AML_vs_normal01_DE.csv"))
  )

  all_de_results[[ct]] <- de_res

  cat("  AML samples:", paste(aml_samples, collapse = ", "), "\n")
  cat("  Up:", sum(grepl("^Up", de_res$direction)),
      " | Down:", sum(grepl("^Down", de_res$direction)), "\n")
}

write_csv(
  bind_rows(kept_samples_summary),
  file.path(out_dir, "cell_counts_by_celltype_and_sample.csv")
)

if (length(all_de_results) > 0) {
  write_csv(
    bind_rows(all_de_results),
    file.path(out_dir, "all_celltypes_AML_vs_normal01_DE_combined.csv")
  )
}

cat("\nDone. Output directory:", out_dir, "\n")
