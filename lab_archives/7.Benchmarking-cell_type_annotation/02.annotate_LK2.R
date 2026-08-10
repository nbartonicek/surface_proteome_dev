#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 02 of 05
#
# The same BoneMarrowMap projection applied to LK2, which is the version that
# went on to become the pipeline's 08_annotate process.
#
# Collated from 8.annotate.R. Three things changed on the way in:
#   - 8.annotate.R computes the composition table three separate times and
#     writes 05/06 (fine), 08 (broad) and 09 (broad again, no palette). 09 is a
#     duplicate of 08 and is dropped.
#   - 10a.heatmap_cell_type_all.pdf comes out 0 bytes because scaling a
#     zero-variance row gives NaN and pheatmap will not draw it. Same fix as
#     step 01: NaN -> 0.
#   - the numbat/barcode export loop at the tail belongs to the CNV work, not
#     to annotation. Dropped here.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
  library(patchwork)
  library(BoneMarrowMap)
  library(symphony)
  library(pheatmap)
})

# ----------------------------
# Paths
# ----------------------------

run          <- "260528_VH01624_464_222K7VKNX"
sample_name  <- "LK2-GEX"
sample_short <- "LK2"

projection_path <- "../annotation/"

in_dir <- file.path("../results/cite_seq_dsb", run, sample_name)

out_dir        <- file.path("../results/seurat_annotated", run)
projection_dir <- file.path(out_dir, "projectionFigures/")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

seurat_rds           <- file.path(in_dir, "seurat_demux_ADT_CLR_DSB_QC.rds")
seurat_projected_rds <- file.path(out_dir, "demux_singlets_annotated_seurat.rds")

if (!file.exists(seurat_rds)) stop("Cannot find input Seurat object: ", seurat_rds)

message("Loading existing Seurat object: ", seurat_rds)
seu <- readRDS(seurat_rds)

# ----------------------------
# Reference
# ----------------------------

ref <- readRDS(paste0(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
ref$save_uwot_path <- paste0(projection_path, "BoneMarrowMap_uwot_model.uwot")

ReferenceSeuratObj <- create_ReferenceObject(ref)

# ----------------------------
# Map, QC, predict
# ----------------------------

batchvar <- "sample_name"

query <- map_Query(query = seu, ref_obj = ref, vars = batchvar)

query <- calculate_MappingError(query, reference = ref, MAD_threshold = 2.5)

pdf(file.path(out_dir, "01_mapping_error_QC.pdf"), width = 8, height = 6)
print(plot_MappingErrorQC(query))
dev.off()

query <- predict_CellTypes(
  query_obj = query, ref_obj = ref, final_label = "predicted_CellType"
)

pdf(file.path(out_dir, "02_projected_predicted_celltypes_pass_only.pdf"), width = 20, height = 12)
print(
  DimPlot(subset(query, mapping_error_QC == "Pass"),
          group.by = "predicted_CellType", label = TRUE, label.size = 4)
)
dev.off()

query <- predict_Pseudotime(
  query_obj = query, ref_obj = ref, final_label = "predicted_Pseudotime"
)

pdf(file.path(out_dir, "03_projected_pseudotime_pass_only.pdf"), width = 8, height = 6)
print(
  FeaturePlot(subset(query, mapping_error_QC == "Pass"),
              features = "predicted_Pseudotime")
)
dev.off()

# ----------------------------
# Per-donor projection panels
# ----------------------------

sample_ids <- levels(factor(query$sample_name))

projection_plots <- lapply(sample_ids, function(sid) {
  q_sub  <- subset(query, subset = sample_name == sid)
  out_i  <- file.path(projection_dir, sid)
  dir.create(out_i, recursive = TRUE, showWarnings = FALSE)
  plot_Projection_byDonor(
    query_obj = q_sub, batch_key = "sample_name",
    ref_obj = ref, save_folder = out_i
  )[[1]] +
    ggtitle(sid)
})

names(projection_plots) <- sample_ids

pdf(file.path(projection_dir, "04_projection_by_sample_name.pdf"), width = 14, height = 10)
print(patchwork::wrap_plots(projection_plots, ncol = 2))
dev.off()

# ----------------------------
# Result tables
# ----------------------------

save_ProjectionResults(
  query_obj = query,
  file_name = file.path(out_dir, "querydata_projected_labeled.csv")
)

write.csv(
  query@meta.data,
  file.path(out_dir, "cell_metadata_demux_ADT_BoneMarrowMap.csv")
)

saveRDS(query, seurat_projected_rds)

# ----------------------------
# Composition, fine labels
# ----------------------------

meta_pass <- query@meta.data %>%
  filter(mapping_error_QC == "Pass") %>%
  filter(!is.na(predicted_CellType), !is.na(sample_name))

celltype_composition <- meta_pass %>%
  count(sample_name, predicted_CellType, name = "n_cells") %>%
  group_by(sample_name) %>%
  mutate(
    total_cells = sum(n_cells),
    fraction = n_cells / total_cells,
    percent = 100 * fraction
  ) %>%
  ungroup() %>%
  arrange(sample_name, desc(percent))

write.csv(celltype_composition,
          file.path(out_dir, "celltype_composition_by_sample_name.csv"),
          row.names = FALSE)

p_comp <- celltype_composition %>%
  ggplot(aes(x = sample_name, y = percent, fill = predicted_CellType)) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(x = "Sample ID", y = "Cellular composition (%)", fill = "Predicted cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "05.celltype_composition_by_sample_name_stacked_bar.pdf"),
       p_comp, width = 10, height = 6)

p_counts <- celltype_composition %>%
  ggplot(aes(x = sample_name, y = n_cells, fill = predicted_CellType)) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(x = "Sample ID", y = "Number of cells", fill = "Predicted cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "06.celltype_counts_by_sample_name_stacked_bar.pdf"),
       p_counts, width = 10, height = 6)

# ----------------------------
# Composition, broad labels and lineage
# ----------------------------
# Keyed on predicted_CellType_Broad. This is the column the lineage vocabulary
# and the bespoke palette are written against - see step 04.

# Shared palette, lineage vocabulary and the check_palette() guard.
# Defined once in step 00 so the copies cannot drift apart again.
source("lab_archives/7.Benchmarking-cell_type_annotation/00.celltype_palette.R")

composition_broad <- meta_pass %>%
  count(sample_name, predicted_CellType_Broad) %>%
  group_by(sample_name) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    lineage = case_when(
      predicted_CellType_Broad %in% c(
        "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
        "Cycling Progenitor", "EoBasoMast Precursor",
        "Megakaryocyte Precursor"
      ) ~ "Stem / progenitor",

      predicted_CellType_Broad %in% c(
        "Monocyte", "Pro-Monocyte", "cDC", "pDC"
      ) ~ "Myeloid / DC",

      predicted_CellType_Broad %in% c(
        "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
        "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
      ) ~ "Lymphoid",

      predicted_CellType_Broad %in% c(
        "Early Erythroid", "Late Erythroid"
      ) ~ "Erythroid",

      TRUE ~ "Other"
    )
  )

write.csv(composition_broad,
          file.path(out_dir, "broad_celltype_composition_by_sample.csv"),
          row.names = FALSE)

composition_lineage <- composition_broad %>%
  group_by(sample_name, lineage) %>%
  summarise(percent = sum(percent), .groups = "drop")

pdf(file.path(out_dir, "07.broad_celltype_composition.pdf"), width = 8, height = 4)
print(
  ggplot(composition_lineage, aes(sample_name, percent, fill = lineage)) +
    geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = lineage_cols) +
    theme_bw() +
    labs(x = "Sample", y = "Cellular composition (%)", fill = "Lineage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
)
dev.off()

composition_broad$predicted_CellType_Broad <- factor(
  composition_broad$predicted_CellType_Broad,
  levels = celltype_order
)

pdf(file.path(out_dir, "08.specific_celltype_composition.pdf"), width = 8, height = 4)
print(
  ggplot(composition_broad, aes(sample_name, percent, fill = predicted_CellType_Broad)) +
    geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = celltype_cols, na.value = "grey80") +
    theme_bw() +
    labs(x = "Sample", y = "Cellular composition (%)", fill = "Broad cell type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
)
dev.off()

# ----------------------------
# Composition heatmaps
# ----------------------------

heat_broad <- composition_broad %>%
  dplyr::select(sample_name, predicted_CellType_Broad, percent) %>%
  pivot_wider(names_from = predicted_CellType_Broad, values_from = percent, values_fill = 0)

mat <- as.matrix(heat_broad[, -1, drop = FALSE])
rownames(mat) <- heat_broad$sample_name
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

pdf(file.path(out_dir, "10.heatmap_cell_type_broad.pdf"), width = 8, height = 3)
pheatmap(mat_scaled, cluster_rows = FALSE, border_color = NA)
dev.off()

heat_fine <- celltype_composition %>%
  dplyr::select(sample_name, predicted_CellType, percent) %>%
  pivot_wider(names_from = predicted_CellType, values_from = percent, values_fill = 0)

mat_fine <- as.matrix(heat_fine[, -1, drop = FALSE])
rownames(mat_fine) <- heat_fine$sample_name
mat_fine_scaled <- t(scale(t(mat_fine)))
mat_fine_scaled[is.na(mat_fine_scaled)] <- 0

pdf(file.path(out_dir, "10a.heatmap_cell_type_all.pdf"), width = 8, height = 3.5)
pheatmap(mat_fine_scaled, cluster_rows = FALSE, border_color = NA)
dev.off()

# ----------------------------
# Reference UMAP in the project palette
# ----------------------------

ReferenceSeuratObj$CellType_Broad <- factor(
  ReferenceSeuratObj$CellType_Broad,
  levels = names(celltype_cols)
)

pdf(file.path(out_dir, "11.database_umap.pdf"), width = 12, height = 10)
print(
  DimPlot(ReferenceSeuratObj, reduction = "umap", group.by = "CellType_Broad",
          raster = FALSE, label = TRUE, repel = TRUE, label.size = 6,
          cols = celltype_cols) +
    theme_bw()
)
dev.off()

# ----------------------------
# Bundle
# ----------------------------

analysis_bundle <- list(
  seu = seu,
  query = query,
  ReferenceSeuratObj = ReferenceSeuratObj,
  ref = ref,
  meta = query@meta.data,
  composition_broad = composition_broad,
  composition_lineage = composition_lineage,
  celltype_composition = celltype_composition,
  celltype_order = celltype_order,
  celltype_cols = celltype_cols,
  lineage_cols = lineage_cols
)

saveRDS(analysis_bundle, file.path(out_dir, "analysis_bundle.rds"))

cat("\nDone.\n")
cat("Output directory:", out_dir, "\n")
cat("Projected object:", seurat_projected_rds, "\n\n")

cat("Mapping QC summary:\n")
print(table(query$mapping_error_QC))
