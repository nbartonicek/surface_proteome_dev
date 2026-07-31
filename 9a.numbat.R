#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(numbat)
  library(data.table)
  library(Matrix)
  library(dplyr)
})

# ----------------------------
# Args
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 9a.numbat.R <donor> [project_path]")
}

donor <- args[1]

PROJ <- if (length(args) >= 2) {
  args[2]
} else if (dir.exists("/scratch/users/nbartonicek/projects/amgen")) {
  "/scratch/users/nbartonicek/projects/amgen"
} else {
  "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
}

RUN <- "260528_VH01624_464_222K7VKNX"
SAMPLE_SHORT <- "LK2"

BASE <- file.path(PROJ, "results/seurat_annotated", RUN, "numbat")
ALLELE_BASE <- file.path(PROJ, "results/seurat_demux", RUN, "numbat")

LABEL <- paste0(SAMPLE_SHORT, "_", donor)

input_dir <- file.path(BASE, "numbat_inputs_no_seurat", LABEL)

counts_file <- file.path(input_dir, paste0(LABEL, "_counts.mtx"))
genes_file <- file.path(input_dir, paste0(LABEL, "_genes.tsv"))
barcodes_file <- file.path(input_dir, paste0(LABEL, "_barcodes.tsv"))
cell_annot_file <- file.path(input_dir, paste0(LABEL, "_cell_annot.tsv"))

allele_file <- file.path(ALLELE_BASE, LABEL, paste0(LABEL, "_allele_counts.tsv.gz"))

out_dir <- file.path(BASE, LABEL, "numbat_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, paste0(LABEL, "_numbat_input_diagnostics.txt"))
sink(log_file, split = TRUE)

cat("Donor:", donor, "\n")
cat("Label:", LABEL, "\n")
cat("Project:", PROJ, "\n")
cat("Input dir:", input_dir, "\n")
cat("Allele file:", allele_file, "\n")
cat("Output:", out_dir, "\n\n")

# ----------------------------
# Robust function patch helpers
# ----------------------------

patch_numbat_if_guard <- function(
    fname,
    variable,
    threshold_pattern = "[0-9.]+"
) {
  ns <- asNamespace("numbat")
  
  if (!exists(fname, envir = ns, inherits = FALSE)) {
    cat(fname, "not found in numbat namespace; skipping patch.\n")
    return(invisible(FALSE))
  }
  
  old_fun <- get(fname, envir = ns)
  body_txt <- deparse(body(old_fun))
  
  cat("\nBefore patch:", fname, "lines containing", variable, "\n")
  print(grep(variable, body_txt, value = TRUE))
  
  # Match:
  # if (mse > 0.5)
  # if(mse>0.5)
  # if ( hom_rate > 0.4 )
  pattern <- paste0(
    "if[[:space:]]*\\([[:space:]]*",
    variable,
    "[[:space:]]*>[[:space:]]*(",
    threshold_pattern,
    ")[[:space:]]*\\)"
  )
  
  replacement <- paste0(
    "if (length(",
    variable,
    ") > 0 && !is.na(",
    variable,
    ") && is.finite(",
    variable,
    ") && ",
    variable,
    " > \\1)"
  )
  
  body_txt2 <- gsub(pattern, replacement, body_txt)
  
  if (identical(body_txt, body_txt2)) {
    cat("WARNING: No replacement made for", fname, "\n")
    cat("Function body may use different formatting. Relevant lines:\n")
    print(grep(variable, body_txt, value = TRUE))
    return(invisible(FALSE))
  }
  
  unlockBinding(fname, ns)
  body(old_fun) <- parse(text = paste(body_txt2, collapse = "\n"))[[1]]
  assign(fname, old_fun, envir = ns)
  lockBinding(fname, ns)
  
  patched_txt <- deparse(body(get(fname, envir = ns)))
  
  cat("\nAfter patch:", fname, "lines containing", variable, "\n")
  print(grep(variable, patched_txt, value = TRUE))
  
  # Defensive verification
  still_bad <- any(grepl(pattern, patched_txt))
  
  if (still_bad) {
    cat("WARNING: unguarded if statement may still remain in", fname, "\n")
  } else {
    cat("Patch verified for", fname, "\n")
  }
  
  invisible(TRUE)
}

patch_numbat_functions <- function() {
  ok1 <- patch_numbat_if_guard(
    fname = "check_contam",
    variable = "hom_rate",
    threshold_pattern = "[0-9.]+"
  )
  
  ok2 <- patch_numbat_if_guard(
    fname = "check_exp_noise",
    variable = "mse",
    threshold_pattern = "[0-9.]+"
  )
  
  cat("\nPatch summary:\n")
  cat("  check_contam patched:", ok1, "\n")
  cat("  check_exp_noise patched:", ok2, "\n\n")
  
  invisible(ok1 && ok2)
}

patch_numbat_functions()

# ----------------------------
# Check files
# ----------------------------

required_files <- c(
  counts_file,
  genes_file,
  barcodes_file,
  cell_annot_file,
  allele_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

# ----------------------------
# Load expression matrix
# ----------------------------

expr <- Matrix::readMM(counts_file)
expr <- as(expr, "CsparseMatrix")

genes <- fread(genes_file, header = FALSE)$V1
barcodes <- fread(barcodes_file, header = FALSE)$V1

if (length(genes) != nrow(expr)) {
  stop("Number of genes does not match matrix rows.")
}

if (length(barcodes) != ncol(expr)) {
  stop("Number of barcodes does not match matrix columns.")
}

rownames(expr) <- make.unique(as.character(genes))
colnames(expr) <- as.character(barcodes)

cat("\nExpression matrix loaded:\n")
cat("  genes:", nrow(expr), "\n")
cat("  cells:", ncol(expr), "\n")
cat("  first genes:", paste(head(rownames(expr)), collapse = ", "), "\n")
cat("  first cells:", paste(head(colnames(expr)), collapse = ", "), "\n\n")

keep_genes <- Matrix::rowSums(expr > 0) >= 10
expr <- expr[keep_genes, , drop = FALSE]

keep_cells <- Matrix::colSums(expr) > 0
expr <- expr[, keep_cells, drop = FALSE]

cat("After basic filtering:\n")
cat("  genes:", nrow(expr), "\n")
cat("  cells:", ncol(expr), "\n\n")

# ----------------------------
# Load cell annotation
# ----------------------------

cell_annot <- fread(cell_annot_file)
cell_annot <- as.data.frame(cell_annot)

required_cols <- c("cell", "sample", "clone", "cell_type")
missing_cols <- setdiff(required_cols, colnames(cell_annot))

if (length(missing_cols) > 0) {
  stop("Missing columns in cell annotation: ", paste(missing_cols, collapse = ", "))
}

cell_annot$cell <- as.character(cell_annot$cell)
rownames(cell_annot) <- cell_annot$cell

common_cells <- intersect(colnames(expr), cell_annot$cell)

cat("Cell annotation loaded:\n")
cat("  annotation rows:", nrow(cell_annot), "\n")
cat("  common cells expr/annotation:", length(common_cells), "\n\n")

if (length(common_cells) < 100) {
  stop("Too few common cells between expression matrix and annotation.")
}

expr <- expr[, common_cells, drop = FALSE]
cell_annot <- cell_annot[common_cells, , drop = FALSE]

# ----------------------------
# Load allele counts
# ----------------------------

df_allele <- fread(allele_file, fill = TRUE)
df_allele <- as.data.frame(df_allele)

if (!"cell" %in% colnames(df_allele)) {
  stop("Allele file does not contain a 'cell' column.")
}

df_allele$cell <- as.character(df_allele$cell)

common_cells2 <- intersect(colnames(expr), unique(df_allele$cell))

cat("Allele table loaded:\n")
cat("  allele rows:", nrow(df_allele), "\n")
cat("  allele cells:", length(unique(df_allele$cell)), "\n")
cat("  common cells expr/allele:", length(common_cells2), "\n\n")

if (length(common_cells2) < 100) {
  stop("Too few common cells between expression matrix and allele counts.")
}

expr <- expr[, common_cells2, drop = FALSE]
cell_annot <- cell_annot[common_cells2, , drop = FALSE]
df_allele <- df_allele[df_allele$cell %in% common_cells2, , drop = FALSE]

# ----------------------------
# Clean allele table
# ----------------------------

df_allele <- df_allele[!is.na(df_allele$cell), , drop = FALSE]

if ("AD" %in% colnames(df_allele)) {
  df_allele <- df_allele[!is.na(df_allele$AD), , drop = FALSE]
}

if ("DP" %in% colnames(df_allele)) {
  df_allele <- df_allele[!is.na(df_allele$DP) & df_allele$DP > 0, , drop = FALSE]
}

if ("GT" %in% colnames(df_allele)) {
  df_allele <- df_allele[
    !is.na(df_allele$GT) &
      df_allele$GT %in% c("0|1", "1|0", "0/1", "1/0"),
    ,
    drop = FALSE
  ]
}

allele_per_cell <- table(df_allele$cell)
good_allele_cells <- names(allele_per_cell)[allele_per_cell >= 50]

cat("Cells with >=50 allele rows:", length(good_allele_cells), "\n")

expr <- expr[, intersect(colnames(expr), good_allele_cells), drop = FALSE]
cell_annot <- cell_annot[colnames(expr), , drop = FALSE]
df_allele <- df_allele[df_allele$cell %in% colnames(expr), , drop = FALSE]

cat("After allele filtering:\n")
cat("  cells:", ncol(expr), "\n")
cat("  allele rows:", nrow(df_allele), "\n")
cat("  allele cells:", length(unique(df_allele$cell)), "\n\n")

if (ncol(expr) < 100) {
  stop("Too few cells after allele filtering.")
}

# ----------------------------
# Build reference expression
# ----------------------------

normal_ref_types <- c(
  "B",
  "CD4 Memory T",
  "CD8 Memory T",
  "Naive T",
  "NK"
)

ref_cells <- rownames(cell_annot)[
  cell_annot$cell_type %in% normal_ref_types
]

ref_cells <- intersect(ref_cells, colnames(expr))

cat("Candidate internal reference cells:", length(ref_cells), "\n")
cat("Cell type table:\n")
print(table(cell_annot$cell_type, useNA = "ifany"))
cat("\n")

if (length(ref_cells) >= 100) {
  
  ref_annot <- data.frame(
    cell = ref_cells,
    group = as.character(cell_annot[ref_cells, "cell_type"]),
    stringsAsFactors = FALSE
  )
  
  lambdas_ref <- numbat::aggregate_counts(
    expr[, ref_cells, drop = FALSE],
    ref_annot
  )
  
  cat("Using internal reference matrix from normal-like cells.\n")
  cat("Reference cell types:\n")
  print(table(ref_annot$group))
  
} else {
  
  lambdas_ref <- numbat::ref_hca
  cat("Too few internal reference cells. Using numbat::ref_hca.\n")
}

cat("Reference dimensions before gene matching:\n")
print(dim(lambdas_ref))
cat("\n")

# ----------------------------
# Match genes
# ----------------------------

common_genes <- intersect(rownames(expr), rownames(lambdas_ref))

cat("Common genes between expr and lambdas_ref:", length(common_genes), "\n\n")

if (length(common_genes) < 3000) {
  cat("First expr genes:\n")
  print(head(rownames(expr), 20))
  cat("First reference genes:\n")
  print(head(rownames(lambdas_ref), 20))
  stop("Too few common genes. Gene identifiers may be Ensembl IDs vs symbols.")
}

expr <- expr[common_genes, , drop = FALSE]
lambdas_ref <- lambdas_ref[common_genes, , drop = FALSE]

keep_genes2 <- Matrix::rowSums(expr) > 0
expr <- expr[keep_genes2, , drop = FALSE]
lambdas_ref <- lambdas_ref[rownames(expr), , drop = FALSE]

cat("Final expression/reference dimensions:\n")
cat("  expr:", paste(dim(expr), collapse = " x "), "\n")
cat("  lambdas_ref:", paste(dim(lambdas_ref), collapse = " x "), "\n\n")

# ----------------------------
# Allele diagnostics
# ----------------------------

cat("Allele diagnostics before run_numbat:\n")
print(colnames(df_allele))
print(head(df_allele))

cat("\nGT table:\n")
print(table(df_allele$GT, useNA = "ifany"))

cat("\nDP summary:\n")
print(summary(df_allele$DP))

cat("\nAD summary:\n")
print(summary(df_allele$AD))

bulk_test <- df_allele |>
  dplyr::group_by(snp_id, CHROM, POS, REF, ALT, GT) |>
  dplyr::summarise(
    AD = sum(AD, na.rm = TRUE),
    DP = sum(DP, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(AR = AD / DP)

cat("\nPseudo-bulk SNP diagnostics:\n")
cat("bulk SNPs total:", nrow(bulk_test), "\n")
cat("bulk SNPs DP >= 8:", sum(bulk_test$DP >= 8, na.rm = TRUE), "\n")
cat("bulk SNPs finite AR:", sum(is.finite(bulk_test$AR)), "\n")
cat("bulk SNPs DP>=8 and finite AR:", sum(bulk_test$DP >= 8 & is.finite(bulk_test$AR)), "\n")
cat("hom_rate test:\n")
print(mean(na.omit((bulk_test$AR == 0 | bulk_test$AR == 1)[bulk_test$DP >= 8])))

good_snps <- bulk_test |>
  dplyr::filter(
    DP >= 8,
    is.finite(AR),
    AR > 0,
    AR < 1
  ) |>
  dplyr::pull(snp_id)

cat("Good informative SNPs:", length(good_snps), "\n")

if (length(good_snps) >= 1000) {
  df_allele <- df_allele[df_allele$snp_id %in% good_snps, , drop = FALSE]
  cat("After informative SNP filtering:\n")
  cat("  allele rows:", nrow(df_allele), "\n")
  cat("  allele cells:", length(unique(df_allele$cell)), "\n")
} else {
  cat("Too few informative SNPs to filter; keeping original allele table.\n")
}

# ----------------------------
# Final checks
# ----------------------------

if (anyDuplicated(rownames(expr)) > 0) {
  stop("Duplicated gene names remain in expr.")
}

if (anyDuplicated(colnames(expr)) > 0) {
  stop("Duplicated cell names remain in expr.")
}

if (nrow(expr) < 3000) {
  stop("Too few genes for Numbat.")
}

if (nrow(df_allele) < 1000) {
  stop("Too few allele observations for Numbat.")
}

saveRDS(
  list(
    expr = expr,
    lambdas_ref = lambdas_ref,
    df_allele = df_allele,
    cell_annot = cell_annot,
    bulk_test = bulk_test
  ),
  file.path(out_dir, paste0(LABEL, "_numbat_inputs_checked.rds"))
)

# ----------------------------
# Run Numbat
# ----------------------------

cat("\nStarting Numbat...\n\n")

out <- numbat::run_numbat(
  count_mat = expr,
  lambdas_ref = lambdas_ref,
  df_allele = df_allele,
  genome = "hg38",
  
  t = 1e-3,
  gamma = 20,
  init_k = 1,
  max_iter = 1,
  min_cells = 50,
  min_LLR = 1,
  max_entropy = 1,
  
  ncores = 8,
  ncores_nni = 8,
  
  plot = FALSE,
  call_clonal_loh = TRUE,
  common_diploid = FALSE,
  multi_allelic = FALSE,
  check_convergence = FALSE,
  
  out_dir = out_dir
)

saveRDS(out, file.path(out_dir, paste0(LABEL, "_run_numbat_output.rds")))

cat("\nDone.\n")
cat("Output directory:", out_dir, "\n")
cat("Saved checked inputs:", file.path(out_dir, paste0(LABEL, "_numbat_inputs_checked.rds")), "\n")
cat("Saved Numbat output:", file.path(out_dir, paste0(LABEL, "_run_numbat_output.rds")), "\n")

sink()