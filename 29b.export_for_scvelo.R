library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

RESULTS  <- "/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX"
SAMPLE   <- "LK2-GEX"
OUT_DIR  <- file.path(RESULTS, "29_velocity", SAMPLE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load final metadata (has donor, cell type, UMAP, CopyKAT/Numbat calls)
meta <- read_csv(file.path(RESULTS, "11_final", SAMPLE, "tables", "LK2_final_metadata.csv")) %>%
  rename(barcode = 1) %>%
  filter(mapping_error_QC == "Pass")

# Load projected UMAP coordinates
umap <- read_csv(file.path(RESULTS, "08_annotate", SAMPLE, "tables", "querydata_projected_labeled.csv")) %>%
  rename(barcode = 1) %>%
  select(barcode, UMAP1_projected, UMAP2_projected)

meta <- meta %>% left_join(umap, by = "barcode")

# Load raw RNA counts from CellRanger filtered matrix
cr_mat_dir <- file.path(RESULTS, "01_cellranger", SAMPLE, SAMPLE, "outs", "filtered_feature_bc_matrix")
seu <- CreateSeuratObject(counts = Read10X(cr_mat_dir))

# Subset to QC-passing cells
keep <- intersect(colnames(seu), meta$barcode)
seu  <- seu[, keep]
meta <- meta %>% filter(barcode %in% keep) %>% tibble::column_to_rownames("barcode")
seu  <- AddMetaData(seu, meta)

# Write counts, metadata, gene names, and PCA for scVelo
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
writeMM(counts_mat, file.path(OUT_DIR, "counts.mtx"))

write_csv(
  seu@meta.data %>% tibble::rownames_to_column("barcode"),
  file.path(OUT_DIR, "metadata.csv")
)

write_lines(rownames(counts_mat), file.path(OUT_DIR, "gene_names.csv"))

message("Exported ", ncol(seu), " cells to ", OUT_DIR)
