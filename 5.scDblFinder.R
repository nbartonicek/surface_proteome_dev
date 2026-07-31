library(Seurat)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)

run <- "260528_VH01624_464_222K7VKNX"
#run <- "260522_VH01624_461_222JLJVNX"
sample_name <- "LK2-GEX"
#sample_name <- "LK1-GEX"
raw_dir <- paste0(
  "../results/cellranger_withbam/", run, "/", sample_name,
  "/outs/raw_feature_bc_matrix/"
)
in_dir <- paste0(
  "../results/emptydrops/", run, "/", sample_name
)
out_dir <- paste0(
  "../results/scDblFinder/", run, "/", sample_name
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load raw matrix
# -----------------------------

gex_counts <- Read10X(raw_dir)

if (is.list(gex_counts)) {
  gex_counts <- gex_counts[["Gene Expression"]]
}

# -----------------------------
# Load EmptyDrops barcodes
# -----------------------------

emptydrops_barcodes <- readLines(
  file.path(in_dir, "emptydrops_barcodes.txt")
)

emptydrops_barcodes <- intersect(emptydrops_barcodes, colnames(gex_counts))

filtered_counts <- gex_counts[, emptydrops_barcodes, drop = FALSE]

# -----------------------------
# Create Seurat object
# -----------------------------

seu <- CreateSeuratObject(
  counts = filtered_counts,
  project = sample_name,
  min.cells = 3,
  min.features = 100
)

# -----------------------------
# Standard RNA QC
# -----------------------------

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

pdf(file.path(out_dir, "RNA_QC_before_filtering.pdf"), width = 8, height = 4)
print(
  VlnPlot(
    seu,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    pt.size = 0,
    ncol = 3
  )
)
dev.off()

# Adjust thresholds after inspecting QC
seu <- subset(
  seu,
  subset =
    nFeature_RNA > 100 
#    percent.mt < 30
)

# -----------------------------
# Run RNA-only scDblFinder
# -----------------------------

sce <- SingleCellExperiment(
  assays = list(
    counts = GetAssayData(seu, assay = "RNA", layer = "counts")
  )
)

sce$sample_id <- sample_name

set.seed(123)
sce <- scDblFinder(sce,dbr = 0.1)

#hist(sce$scDblFinder.score)
pdf(file.path(out_dir, "hist_score.pdf"), width = 8, height = 4)
hist(sce$scDblFinder.score)
dev.off()
# -----------------------------
# Add results back to Seurat
# -----------------------------

seu$scDblFinder.score <- sce$scDblFinder.score
seu$scDblFinder.class <- sce$scDblFinder.class

pdf(file.path(out_dir, "scDblFinder_QC.pdf"), width = 8, height = 4)
print(
  VlnPlot(
    seu,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder.score"),
    group.by = "scDblFinder.class",
    pt.size = 0,
    ncol = 4
  )
)
dev.off()

# -----------------------------
# Keep transcriptomic singlets
# -----------------------------

seu_filtered <- subset(
  seu,
  subset = scDblFinder.class == "singlet"
)

# -----------------------------
# Save
# -----------------------------

saveRDS(
  seu_filtered,
  file = file.path(out_dir, "seurat_emptyDrops_RNA_scDblFinder_filtered.rds")
)

saveRDS(
  seu,
  file = file.path(out_dir, "seurat_emptyDrops_RNA_scDblFinder_nonfiltered.rds")
)

write.csv(
  seu_filtered@meta.data,
  file = file.path(out_dir, "cell_metadata_after_RNA_scDblFinder.csv")
)

message("Starting cells after EmptyDrops + min.features: ", ncol(seu))
message("Final cells after RNA QC + scDblFinder: ", ncol(seu_filtered))
