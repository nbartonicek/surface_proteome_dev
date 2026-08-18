#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Calling cancer cells - step 01
#
# CopyKAT on each demultiplexed donor separately, using the lymphoid
# populations as the confident-normal reference.
#
# Collated from 9.copykat.R and 9a.copykat_perSample.R. 9.copykat.R ran CopyKAT
# over the whole GEM well at once, which is wrong here - the well is a pool of
# four donors, and CopyKAT's baseline is estimated from the cells you give it,
# so pooling donors lets one donor's aneuploidy define another's baseline. The
# per-donor version is what produced the surviving outputs in
# results/.../copykat_annotated/<donor>/, so that is what this step keeps.
#
# Four things were cleaned up on the way in:
#   - the donor loop was `unique(temp$sample_name)[c(4)]`, a hardcoded index, so
#     the script processed exactly one donor per invocation and the index was
#     edited by hand between runs. It now loops over all donors, or over the
#     ones named on the command line.
#   - CopyKAT writes its output files into the working directory and names them
#     from `sam.name`, which is why the original did a bare setwd() to an
#     absolute path. That is kept (there is no other way to steer CopyKAT's
#     output) but the original directory is now restored with on.exit, so a
#     failure part-way through does not leave the session in the wrong place.
#   - `copykat_res <- copykat(...)` was assigned and then never used - the
#     results were re-read from the files CopyKAT had just written. The
#     assignment is dropped and the files are read directly.
#   - the closing cat() claimed to have written query_with_copykat_calls.rds and
#     copykat_bundle.rds. Neither is ever written. Removed.
#
# NOTE ON MOLM13: CopyKAT is never run on it. MOLM13 is an AML cell line, so
# every cell is malignant by construction and there are no normal cells in it to
# anchor a baseline. It is assigned "aneuploid" directly. That is a manual
# override, not a call - see step 09, which reports it separately so it cannot
# be mistaken for a result.
#
# Run from the scripts/ directory - paths are relative to it.
#
#   Rscript lab_archives/8.Benchmarking-calling_cancer_cells/01.copykat_per_donor.R
#   Rscript lab_archives/8.Benchmarking-calling_cancer_cells/01.copykat_per_donor.R HBDN206-MNpCT
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(data.table)
  library(copykat)
  library(RColorBrewer)
})

# ----------------------------
# Parameters
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

RUN          <- "260423_VH01624_453_222HWMYNX"
SAMPLE_SHORT <- "LK1"

# Cell line - malignant by construction, see header
CELL_LINES <- c("MOLM13")

# The confident-normal anchor. Lymphoid cells are used because they are not part
# of the myeloid compartment the leukaemia comes from, so they are the
# populations least likely to carry the tumour's CNVs.
NORMAL_CELLTYPES <- c("Naive T", "CD4 Memory T", "CD8 Memory T", "NK", "B")

scripts_dir <- normalizePath(".")

annotation_dir <- file.path("../results/seurat_annotated", RUN)
seurat_projected_rds <- file.path(
  annotation_dir, "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

copykat_out <- file.path(annotation_dir, "copykat_annotated")
dir.create(copykat_out, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seurat_projected_rds)) {
  stop("Cannot find projected Seurat object: ", seurat_projected_rds)
}

# ----------------------------
# Load the annotated object
# ----------------------------

full <- readRDS(seurat_projected_rds)
DefaultAssay(full) <- "RNA"

stopifnot("RNA" %in% Assays(full))
stopifnot("sampleID" %in% colnames(full@meta.data))

cat("Loaded cells:", ncol(full), "\n")
cat("RNA genes:", nrow(full[["RNA"]]), "\n")

# Only cells the reference could place - a cell with no cell type cannot be used
# as a normal anchor and cannot be interpreted afterwards.
if ("mapping_error_QC" %in% colnames(full@meta.data)) {
  full <- subset(full, subset = mapping_error_QC == "Pass")
}

cat("Cells after mapping QC:", ncol(full), "\n")

donors <- sort(unique(as.character(full$sample_name)))
if (length(args)) donors <- intersect(args, donors)
donors_to_run <- setdiff(donors, CELL_LINES)

cat("Donors to run CopyKAT on:", paste(donors_to_run, collapse = ", "), "\n")
cat("Donors assigned aneuploid without running:",
    paste(intersect(donors, CELL_LINES), collapse = ", "), "\n\n")

# ----------------------------
# CopyKAT, one donor at a time
# ----------------------------

for (donor in donors_to_run) {

  cat("=====", donor, "=====\n")

  copykat_out_sample <- file.path(copykat_out, donor)
  dir.create(copykat_out_sample, recursive = TRUE, showWarnings = FALSE)

  seu <- subset(full, cells = colnames(full)[full$sample_name == donor])

  norm_cells <- rownames(seu@meta.data)[
    seu$predicted_CellType_Broad %in% NORMAL_CELLTYPES
  ]

  cat("  cells:", ncol(seu), " normal anchor cells:", length(norm_cells), "\n")

  if (length(norm_cells) < 20) {
    warning("Fewer than 20 normal anchor cells for ", donor,
            " - CopyKAT's baseline will be unstable.", call. = FALSE)
  }

  raw_counts_dense <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))

  # CopyKAT writes into the working directory and names files from sam.name,
  # so the only way to place its output is to be standing in the right folder.
  old_wd <- setwd(copykat_out_sample)
  on.exit(setwd(old_wd), add = TRUE)

  copykat(
    rawmat = raw_counts_dense,
    id.type = "S",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = paste0(SAMPLE_SHORT, "_copykat"),
    distance = "euclidean",
    norm.cell.names = norm_cells,
    output.seg = "FALSE",
    plot.genes = "FALSE",
    genome = "hg20",
    n.cores = 2
  )

  setwd(old_wd)

  # ----------------------------
  # Read back what CopyKAT wrote
  # ----------------------------

  pred_file <- file.path(copykat_out_sample,
                         paste0(SAMPLE_SHORT, "_copykat_copykat_prediction.txt"))

  if (!file.exists(pred_file)) {
    warning("CopyKAT produced no prediction file for ", donor, " - skipping.",
            call. = FALSE)
    next
  }

  copykat_pred <- read.table(pred_file, header = TRUE) %>%
    as_tibble() %>%
    rename(cell = cell.names, copykat_call = copykat.pred)

  seu$copykat_call <- NA_character_
  common_cells <- intersect(colnames(seu), copykat_pred$cell)
  seu$copykat_call[common_cells] <-
    copykat_pred$copykat_call[match(common_cells, copykat_pred$cell)]

  seu$copykat_malignancy <- case_when(
    seu$copykat_call == "aneuploid" ~ "CNV_aberrant",
    seu$copykat_call == "diploid"   ~ "CNV_neutral",
    is.na(seu$copykat_call)         ~ "Not_called",
    TRUE                            ~ seu$copykat_call
  )

  # ----------------------------
  # Tables
  # ----------------------------

  copykat_annotated <- copykat_pred %>%
    mutate(
      sampleID = seu$sampleID[cell],
      predicted_CellType = if ("predicted_CellType" %in% colnames(seu@meta.data)) {
        seu$predicted_CellType[cell]
      } else NA_character_,
      predicted_CellType_Broad = if ("predicted_CellType_Broad" %in% colnames(seu@meta.data)) {
        seu$predicted_CellType_Broad[cell]
      } else NA_character_
    )

  write.csv(copykat_annotated,
            file.path(copykat_out_sample, "copykat_prediction_with_metadata.csv"),
            row.names = FALSE)

  sample_summary <- copykat_annotated %>%
    count(sampleID, copykat_call, name = "n_cells") %>%
    group_by(sampleID) %>%
    mutate(percent = 100 * n_cells / sum(n_cells)) %>%
    ungroup()

  write.csv(sample_summary,
            file.path(copykat_out_sample, "copykat_call_composition_by_sampleID.csv"),
            row.names = FALSE)

  celltype_summary <- copykat_annotated %>%
    count(sampleID, predicted_CellType_Broad, copykat_call, name = "n_cells") %>%
    group_by(sampleID, predicted_CellType_Broad) %>%
    mutate(percent = 100 * n_cells / sum(n_cells)) %>%
    ungroup()

  write.csv(celltype_summary,
            file.path(copykat_out_sample, "copykat_call_by_sampleID_and_celltype.csv"),
            row.names = FALSE)

  # ----------------------------
  # Plots
  # ----------------------------

  cols <- brewer.pal(n = 8, name = "Set1")
  copykat_levels <- levels(factor(seu$copykat_call))
  copykat_cols <- setNames(cols[seq_along(copykat_levels)], copykat_levels)

  seu$copykat_call  <- factor(seu$copykat_call, levels = copykat_levels)
  sample_summary    <- mutate(sample_summary,   copykat_call = factor(copykat_call, levels = copykat_levels))
  celltype_summary  <- mutate(celltype_summary, copykat_call = factor(copykat_call, levels = copykat_levels))

  pdf(file.path(copykat_out_sample, "01_copykat_call_composition_by_sampleID.pdf"),
      width = 8, height = 5)
  print(
    ggplot(sample_summary, aes(sampleID, percent, fill = copykat_call)) +
      geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
      theme_bw() +
      labs(x = "Sample", y = "Cells (%)", fill = "CopyKAT call",
           title = "CopyKAT malignant/normal-like composition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()

  if ("umap_projected" %in% Reductions(seu)) {

    pdf(file.path(copykat_out_sample, "02_copykat_calls_on_projected_RNA_UMAP.pdf"),
        width = 8, height = 6)
    print(
      DimPlot(seu, reduction = "umap_projected", group.by = "copykat_call",
              cols = copykat_cols) +
        ggtitle("CopyKAT calls on projected RNA UMAP")
    )
    dev.off()

    pdf(file.path(copykat_out_sample, "02a_copykat_calls_on_projected_RNA_UMAP.pdf"),
        width = 8, height = 6)
    print(
      DimPlot(seu, reduction = "umap_projected", group.by = "copykat_call",
              split.by = "sampleID", cols = copykat_cols, ncol = 2) +
        ggtitle("CopyKAT calls on projected RNA UMAP by sampleID")
    )
    dev.off()
  }

  if ("predicted_CellType_Broad" %in% colnames(seu@meta.data)) {

    pdf(file.path(copykat_out_sample, "03_copykat_by_BoneMarrowMap_broad_celltype.pdf"),
        width = 10, height = 5)
    print(
      ggplot(celltype_summary,
             aes(predicted_CellType_Broad, percent, fill = copykat_call)) +
        geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
        facet_wrap(~ sampleID) +
        scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
        theme_bw() +
        labs(x = "BoneMarrowMap broad cell type", y = "Cells (%)",
             fill = "CopyKAT call") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank())
    )
    dev.off()

    pdf(file.path(copykat_out_sample, "03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf"),
        width = 10, height = 5)
    print(
      seu@meta.data %>%
        as_tibble() %>%
        filter(!is.na(predicted_CellType_Broad), !is.na(copykat_call)) %>%
        count(sampleID, predicted_CellType_Broad, copykat_call) %>%
        ggplot(aes(predicted_CellType_Broad, n, fill = copykat_call)) +
        geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
        facet_wrap(~ sampleID) +
        scale_fill_manual(values = copykat_cols, drop = FALSE, na.value = "grey80") +
        theme_bw() +
        labs(x = "BoneMarrowMap broad cell type", y = "Cell count",
             fill = "CopyKAT call") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank())
    )
    dev.off()
  }

  rm(seu, raw_counts_dense); gc(verbose = FALSE)
}

# ----------------------------
# Collect every donor's calls back onto the full object
# ----------------------------

read_copykat_calls <- function(d) {

  csv_file <- file.path(d, "copykat_prediction_with_metadata.csv")
  txt_file <- file.path(d, paste0(SAMPLE_SHORT, "_copykat_copykat_prediction.txt"))

  if (file.exists(csv_file)) {
    read.csv(csv_file, stringsAsFactors = FALSE) %>% dplyr::select(cell, copykat_call)
  } else if (file.exists(txt_file)) {
    read.delim(txt_file, stringsAsFactors = FALSE) %>%
      dplyr::rename(cell = cell.names, copykat_call = copykat.pred)
  } else {
    warning("No CopyKAT prediction file found in: ", d, call. = FALSE)
    NULL
  }
}

copykat_dirs <- list.dirs(copykat_out, recursive = FALSE, full.names = TRUE)
copykat_dirs <- copykat_dirs[basename(copykat_dirs) %in% donors]

copykat_all <- purrr::map_dfr(copykat_dirs, read_copykat_calls) %>%
  dplyr::distinct(cell, .keep_all = TRUE) %>%
  dplyr::mutate(
    copykat_malignancy = dplyr::case_when(
      copykat_call == "aneuploid" ~ "CNV_aberrant",
      copykat_call == "diploid"   ~ "CNV_neutral",
      is.na(copykat_call)         ~ "Not_called",
      TRUE                        ~ as.character(copykat_call)
    ),
    copykat_source = "CopyKAT"
  )

# The cell-line override, kept explicit and flagged in its own column so that
# anything downstream can exclude it from concordance statistics.
for (cl in intersect(donors, CELL_LINES)) {
  cl_cells <- colnames(full)[full$sample_name == cl]
  copykat_all <- dplyr::bind_rows(
    copykat_all,
    tibble(cell = cl_cells,
           copykat_call = "aneuploid",
           copykat_malignancy = "CNV_aberrant",
           copykat_source = "cell line, assigned not called")
  )
}

copykat_all <- dplyr::distinct(copykat_all, cell, .keep_all = TRUE)

full$copykat_call       <- NA_character_
full$copykat_malignancy <- "Not_called"
full$copykat_source     <- NA_character_

common_cells <- intersect(colnames(full), copykat_all$cell)
idx <- match(common_cells, copykat_all$cell)

full$copykat_call[common_cells]       <- copykat_all$copykat_call[idx]
full$copykat_malignancy[common_cells] <- copykat_all$copykat_malignancy[idx]
full$copykat_source[common_cells]     <- copykat_all$copykat_source[idx]

saveRDS(full, file.path(copykat_out, "temp_with_all_copykat_calls.rds"))

write.csv(full@meta.data,
          file.path(copykat_out, "temp_with_all_copykat_calls_metadata.csv"),
          row.names = TRUE)

cat("\nDone.\n")
cat("CopyKAT output:", copykat_out, "\n\n")
cat("Calls by donor:\n")
print(table(full$sample_name, full$copykat_malignancy, useNA = "ifany"))
