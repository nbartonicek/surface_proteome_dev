library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

RESULTS <- "/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX"
SAMPLE  <- "LK2-GEX"
OUT_DIR <- file.path(RESULTS, "velocity", SAMPLE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Final metadata covers all 4 donors in the LK2 run (HBDN498-TP53, APOP576-TP53,
# HBDN376-TP53, normal-01). Each barcode has donor assignment (sample_name),
# broad cell type, CopyKAT/Numbat calls, and mapping QC pass/fail.
meta <- read_csv(file.path(RESULTS, "11_final", SAMPLE, "tables", "LK2_final_metadata.csv"),
                 show_col_types = FALSE) %>%
  rename(barcode = 1) %>%
  filter(mapping_error_QC == "Pass")

message(sprintf("Donors in this run: %s", paste(sort(unique(meta$sample_name)), collapse = ", ")))
message(sprintf("Total QC-passing cells: %d", nrow(meta)))

# Projected UMAP coordinates from BoneMarrowMap (08_annotate)
umap <- read_csv(file.path(RESULTS, "08_annotate", SAMPLE, "tables", "querydata_projected_labeled.csv"),
                 show_col_types = FALSE) %>%
  rename(barcode = 1) %>%
  select(barcode, UMAP1_projected, UMAP2_projected)

meta <- meta %>% left_join(umap, by = "barcode")

# Raw RNA counts from CellRanger filtered matrix (all barcodes — subset below)
cr_mat_dir <- file.path(RESULTS, "01_cellranger", SAMPLE, SAMPLE, "outs", "filtered_feature_bc_matrix")
seu <- CreateSeuratObject(counts = Read10X(cr_mat_dir))

# Subset to QC-passing cells and attach full metadata
keep <- intersect(colnames(seu), meta$barcode)
seu  <- seu[, keep]
meta_mat <- meta %>% filter(barcode %in% keep) %>% tibble::column_to_rownames("barcode")
seu  <- AddMetaData(seu, meta_mat)

message(sprintf("Cells after intersection with CellRanger matrix: %d", ncol(seu)))

# Write outputs for scVelo
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
writeMM(counts_mat, file.path(OUT_DIR, "counts.mtx"))

write_csv(
  seu@meta.data %>% tibble::rownames_to_column("barcode"),
  file.path(OUT_DIR, "metadata.csv")
)

write_lines(rownames(counts_mat), file.path(OUT_DIR, "gene_names.csv"))

message("Done. Exported to: ", OUT_DIR)
