
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(dsb)
})

# ----------------------------
# Paths
# ----------------------------

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"
sample_short <- "LK1"

seurat_dir <- file.path("../results/scDblFinder", run, sample_name)

seurat_file <- file.path(
  seurat_dir,
  "seurat_emptyDrops_RNA_scDblFinder_nonfiltered.rds"
)

raw_dir <- paste0(
  "../results/cellranger/", run, "/", sample_name,
  "/outs/raw_feature_bc_matrix/"
)

adt_dir <- file.path(
  "../results/cite_seq_count",
  run,
  "adt_counts_LK1/umi_count"
)

emptydrops_dir <- file.path(
  "../results/emptydrops",
  run,
  sample_name
)

out_dir <- file.path(
  "../results/cite_dsb_qc",
  run,
  sample_name
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#rna <- Read10X(raw_dir)


adt_mat <- readMM(file.path(adt_dir, "matrix.mtx.gz"))

adt_barcodes <- fread(
  file.path(adt_dir, "barcodes.tsv.gz"),
  header = FALSE
)$V1

adt_features <- fread(
  file.path(adt_dir, "features.tsv.gz"),
  header = FALSE
)

# usually feature name in V2, fallback to V1
adt_feature_names <- if (ncol(adt_features) >= 2) {
  adt_features$V2
} else {
  adt_features$V1
}

adt_barcodes <- paste0(adt_barcodes, "-1")

rownames(adt_mat) <- make.unique(adt_feature_names)
colnames(adt_mat) <- adt_barcodes

# Match cells
common_adt <- intersect(colnames(seu), colnames(adt_mat))

cat("ADT shared cells:", length(common_adt), "\n")

adt_mat <- adt_mat[, common_adt, drop = FALSE]

# subset Seurat object to shared cells
seu <- subset(seu, cells = common_adt)

# add ADT assay
seu[["ADT"]] <- CreateAssayObject(counts = adt_mat)

DefaultAssay(seu) <- "ADT"

seu <- NormalizeData(
  seu,
  normalization.method = "CLR",
  margin = 2
)

DefaultAssay(seu) <- "ADT"

seu$ADT_total <- Matrix::colSums(GetAssayData(seu, assay = "ADT", layer = "counts"))
seu$ADT_features <- Matrix::colSums(GetAssayData(seu, assay = "ADT", layer = "counts") > 0)


vireo_to_hto_map <- meta %>%
  filter(
    hto_class == "Singlet",
    vireo_class == "singlet",
    !is.na(hto_donor),
    !is.na(vireo_donor),
    hto_donor != "Negative",
    vireo_donor != "unassigned"
  ) %>%
  count(vireo_donor, hto_donor, name = "n_cells") %>%
  group_by(vireo_donor) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup()

vireo_to_hto_vec <- setNames(
  vireo_to_hto_map$hto_donor,
  vireo_to_hto_map$vireo_donor
)

# Keep original Vireo donor
seu$vireo_donor_original <- seu$vireo_donor

# New column: Vireo donor relabelled using HTO names
seu$vireo_donor_hto_name <- unname(vireo_to_hto_vec[as.character(seu$vireo_donor)])

# Preserve unassigned/doublets/unknowns
seu$sample_name <- case_when(
  is.na(seu$vireo_donor) ~ "unassigned",
  seu$vireo_donor == "unassigned" ~ "unassigned",
  seu$vireo_class == "doublet" ~ "doublet",
  is.na(seu$vireo_donor_hto_name) ~ as.character(seu$vireo_donor),
  TRUE ~ seu$vireo_donor_hto_name
)

print(vireo_to_hto_map)

saveRDS(
  seu,
  file.path(demux_out_dir, "seurat_demultiplexed_CITE.rds")
)
