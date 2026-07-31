library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scIntegrationMetrics)

# remotes::install_github("carmonalab/scIntegrationMetrics")

run <- "260528_VH01624_464_222K7VKNX"
proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

annotation_dir <- file.path(proj, "results/seurat_annotated", run)
seurat_file <- file.path(annotation_dir, "demux_singlets_annotated_seurat.rds")

out_dir <- file.path(annotation_dir, "harmony_parameter_qc_scIntegrationMetrics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(seurat_file)

FeaturePlot(
  seu,
  features = "SNRNP200",
  reduction = "umap_projected",
  split.by = "sample_name",
  order = TRUE
)
