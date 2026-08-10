#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 01 of 05
#
# BoneMarrowMap/Symphony projection of the LK1 pilot run.
#
# Collated from backup/10.annotate.R, backup/10a.annotate.R and
# backup/10b.annotate_DSB.R. 10.annotate.R is superseded - none of its output
# names survive in results/. 10a and 10b are the same script run twice, once
# with the ADT_CLR assay and once with CITE_DSB, for the marker-projection
# panels only; that is the ADT_ASSAY argument here. Both sets of figures are on
# disk (06_ADT_CLR_* and 06_CITE_DSB_*), which is why the argument exists
# rather than a hardcoded choice.
#
# Two things were cleaned up on the way in:
#   - 10a had a marker_map / plot_df block sitting above the code that creates
#     adt_df and plot_df, so it referenced objects that did not exist yet. It
#     is removed; the working copy of that logic is further down.
#   - the numbat export loop at the tail of 10a belongs to the CNV work, not to
#     annotation. Dropped here.
#
# Run from the scripts/ directory - paths are relative to it.
#
#   Rscript lab_archives/7.Benchmarking-cell_type_annotation/01.annotate_LK1_pilot.R CITE_DSB
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
# Parameters
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)
ADT_ASSAY <- if (length(args) >= 1) args[1] else "CITE_DSB"

run          <- "260423_VH01624_453_222HWMYNX"
sample_name  <- "LK1-GEX"
sample_short <- "LK1"

projection_path <- "../annotation/"

cite_dsb_dir <- file.path("../results/cite_qc_then_dsb", run, sample_name)
seurat_input <- file.path(cite_dsb_dir, "seurat_demux_ADT_CLR_DSB_QC.rds")

out_dir        <- file.path("../results/seurat_annotated", run)
projection_dir <- file.path(out_dir, "projectionFigures/")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

seurat_projected_rds <- file.path(
  out_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

# ----------------------------
# Load the CITE/DSB object
# ----------------------------

if (!file.exists(seurat_input)) {
  stop("Cannot find input Seurat object: ", seurat_input)
}

seu <- readRDS(seurat_input)

cat("Loaded Seurat object:\n")
cat("Cells:", ncol(seu), "\n")
cat("Features:", nrow(seu), "\n")
cat("Assays:", paste(Assays(seu), collapse = ", "), "\n\n")

stopifnot("RNA" %in% Assays(seu))
stopifnot("sample_name" %in% colnames(seu@meta.data))

seu$sampleID <- as.character(seu$sample_name)

# Only demultiplexed singlets go into the projection
seu <- subset(
  seu,
  subset = !sampleID %in% c(
    "Doublet", "doublet",
    "Negative", "negative",
    "unassigned", "Unassigned",
    NA
  )
)

seu$sampleID <- factor(seu$sampleID)

cat("Cells after removing doublet/negative/unassigned:\n")
print(table(seu$sampleID))
cat("\n")

# ----------------------------
# BoneMarrowMap / Symphony reference
# ----------------------------
# The reference and its uwot model are downloaded once from
# https://bonemarrowmap.s3.us-east-2.amazonaws.com/ into ../annotation/.

ref_file  <- file.path(projection_path, "BoneMarrowMap_SymphonyReference.rds")
uwot_file <- file.path(projection_path, "BoneMarrowMap_uwot_model.uwot")

if (!file.exists(ref_file))  stop("Missing BoneMarrowMap reference: ", ref_file)
if (!file.exists(uwot_file)) stop("Missing BoneMarrowMap uwot model: ", uwot_file)

ref <- readRDS(ref_file)
ref$save_uwot_path <- uwot_file

ReferenceSeuratObj <- create_ReferenceObject(ref)

pdf(file.path(out_dir, "00_BoneMarrowMap_reference_celltype.pdf"), width = 12, height = 10)
print(
  DimPlot(
    ReferenceSeuratObj,
    reduction = "umap",
    group.by = "CellType_Annotation_formatted",
    raster = FALSE, label = TRUE, repel = TRUE, label.size = 4
  ) +
    ggtitle("BoneMarrowMap reference: formatted annotation")
)
dev.off()

pdf(file.path(out_dir, "00a_BoneMarrowMap_reference_broad_celltype.pdf"), width = 12, height = 10)
print(
  DimPlot(
    ReferenceSeuratObj,
    reduction = "umap",
    group.by = "CellType_Broad",
    raster = FALSE, label = TRUE, repel = TRUE, label.size = 4
  ) +
    ggtitle("BoneMarrowMap reference: broad cell type")
)
dev.off()

# ----------------------------
# Map the query into reference space
# ----------------------------
# Batch variable is the demultiplexed sample, so Symphony harmonises across
# donors within the pool rather than treating the GEM well as one batch.

DefaultAssay(seu) <- "RNA"

batchvar <- "sampleID"
cat("Running Symphony map_Query using batch variable:", batchvar, "\n")

query <- map_Query(query = seu, ref_obj = ref, vars = batchvar)

# MAD_threshold = 2.5 is the BoneMarrowMap default and is what produced the
# figures in results/. Step 03 re-runs the same call at 4 on a later object;
# step 04 compares what the two thresholds cost.
query <- calculate_MappingError(query, reference = ref, MAD_threshold = 2.5)

pdf(file.path(projection_dir, "01_mapping_error_QC.pdf"), width = 8, height = 6)
print(plot_MappingErrorQC(query))
dev.off()

# ----------------------------
# Predict cell types and pseudotime
# ----------------------------

query <- predict_CellTypes(
  query_obj = query,
  ref_obj = ref,
  final_label = "predicted_CellType"
)

query <- predict_Pseudotime(
  query_obj = query,
  ref_obj = ref,
  final_label = "predicted_Pseudotime"
)

# NOTE: predicted_CellType is the FINE label (53 levels: "HSC", "MPP-MyLy",
# "CD14 Mono", "BFU-E", ...). predicted_CellType_Broad is the 24-level
# vocabulary. The lineage grouping below has to key on the broad column - see
# step 04, which quantifies what happens if it is keyed on the fine one.

query$predicted_Lineage <- case_when(
  query$predicted_CellType_Broad %in% c(
    "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
    "Cycling Progenitor", "EoBasoMast Precursor",
    "Megakaryocyte Precursor"
  ) ~ "Stem / progenitor",

  query$predicted_CellType_Broad %in% c(
    "Monocyte", "Pro-Monocyte", "cDC", "pDC"
  ) ~ "Myeloid / DC",

  query$predicted_CellType_Broad %in% c(
    "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
    "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
  ) ~ "Lymphoid",

  query$predicted_CellType_Broad %in% c(
    "Early Erythroid", "Late Erythroid"
  ) ~ "Erythroid",

  TRUE ~ "Other"
)

# ----------------------------
# Projection plots
# ----------------------------

query_pass <- subset(query, subset = mapping_error_QC == "Pass")

pdf(file.path(projection_dir, "02_projected_predicted_celltypes_pass_only.pdf"), width = 20, height = 12)
print(
  DimPlot(
    query_pass, reduction = "umap", group.by = "predicted_CellType",
    label = TRUE, repel = TRUE, label.size = 4
  ) +
    ggtitle("BoneMarrowMap projected cell types: mapping pass only")
)
dev.off()

pdf(file.path(projection_dir, "03_projected_samples_pass_only.pdf"), width = 10, height = 7)
print(
  DimPlot(
    query_pass, reduction = "umap", group.by = "sampleID",
    label = TRUE, repel = TRUE
  ) +
    ggtitle("Samples on BoneMarrowMap projection: mapping pass only")
)
dev.off()

pdf(file.path(projection_dir, "04_projected_pseudotime_pass_only.pdf"), width = 8, height = 6)
print(
  FeaturePlot(
    query_pass, reduction = "umap",
    features = "predicted_Pseudotime", order = TRUE
  ) +
    ggtitle("Predicted pseudotime: mapping pass only")
)
dev.off()

sample_ids <- levels(query_pass$sampleID)

projection_plots <- lapply(sample_ids, function(sid) {
  q_sub <- subset(query_pass, subset = sampleID == sid)
  plot_Projection_byDonor(
    query_obj = q_sub,
    batch_key = "sampleID",
    ref_obj = ref,
    save_folder = file.path(projection_dir, sid)
  )[[1]] +
    ggtitle(sid)
})

names(projection_plots) <- sample_ids

pdf(file.path(projection_dir, "05_projection_by_sampleID_BMM_helper_manual_split.pdf"),
    width = 14, height = 10)
print(patchwork::wrap_plots(projection_plots, ncol = 2))
dev.off()

# ----------------------------
# ADT marker projection sanity check
# ----------------------------
# Does surface protein land where the projected label says it should. Run once
# per assay (ADT_CLR and CITE_DSB) so the two can be compared like for like.

cite_assay <- if (ADT_ASSAY %in% Assays(query)) {
  ADT_ASSAY
} else if ("ADT_CLR" %in% Assays(query)) {
  message("Requested assay '", ADT_ASSAY, "' not present - falling back to ADT_CLR.")
  "ADT_CLR"
} else if ("ADT" %in% Assays(query)) {
  "ADT"
} else {
  NA_character_
}

if (!is.na(cite_assay)) {

  cat("Marker projection assay:", cite_assay, "\n")
  DefaultAssay(query) <- cite_assay

  cite_mat <- GetAssayData(query, assay = cite_assay, layer = "data")

  cite_markers <- rownames(cite_mat)
  cite_markers <- cite_markers[
    !grepl("isotype|igg|control|unmapped", cite_markers, ignore.case = TRUE)
  ]

  marker_means <- Matrix::rowMeans(as.matrix(cite_mat[cite_markers, , drop = FALSE]))
  top_cite_markers <- names(sort(marker_means, decreasing = TRUE))[1:min(16, length(marker_means))]

  pdf(file.path(projection_dir, paste0("06_", cite_assay, "_top_marker_projection.pdf")),
      width = 14, height = 10)
  print(
    FeaturePlot(query, reduction = "umap_projected",
                features = top_cite_markers, ncol = 4, order = TRUE)
  )
  dev.off()

  # Markers of interest for the surfaceome work, plotted by hand rather than
  # through FeaturePlot so the colour scale can be clipped at the 5th/95th
  # percentile - otherwise one bright cell flattens the whole panel.
  extra_markers <- c("CD33", "CD133", "CD47", "Flt3-Flk2", "TACTILE", "GPR56", "CD7")
  extra_markers <- extra_markers[extra_markers %in% rownames(cite_mat)]

  if (length(extra_markers)) {

    emb <- Embeddings(query, "umap_projected") %>%
      as.data.frame() %>%
      rownames_to_column("cell")
    colnames(emb)[2:3] <- c("UMAP_1", "UMAP_2")

    adt_df <- FetchData(query, vars = extra_markers, layer = "data") %>%
      rownames_to_column("cell")

    plot_df <- emb %>%
      left_join(adt_df, by = "cell") %>%
      pivot_longer(cols = all_of(extra_markers), names_to = "marker", values_to = "value") %>%
      group_by(marker) %>%
      mutate(
        value_clip = pmin(
          pmax(value, quantile(value, 0.05, na.rm = TRUE)),
          quantile(value, 0.95, na.rm = TRUE)
        )
      ) %>%
      ungroup()

    # Display names for the panels, restricted to what the panel actually has
    marker_map <- c(
      "CD33" = "CD33", "CD133" = "CD133", "CD47" = "CD47",
      "CD135 (FLT3)" = "Flt3-Flk2", "GPR56" = "GPR56",
      "CD7" = "CD7", "CD96" = "TACTILE"
    )
    marker_map <- marker_map[marker_map %in% extra_markers]

    plot_df$marker <- factor(plot_df$marker, levels = marker_map, labels = names(marker_map))

    p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = value_clip)) +
      geom_point(size = 0.15) +
      facet_wrap(~ marker, ncol = 4) +
      scale_colour_gradientn(colours = c("grey90", "#FDB863", "#E66101", "#B2182B")) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank()) +
      labs(colour = cite_assay, x = "Projected UMAP 1", y = "Projected UMAP 2")

    ggsave(
      file.path(projection_dir, paste0("07_", cite_assay, "_selected_markers_manual_projectedUMAP.pdf")),
      p, width = 14, height = 8
    )
  }
}

DefaultAssay(query) <- "RNA"

# ----------------------------
# Projection result tables
# ----------------------------

save_ProjectionResults(
  query_obj = query,
  file_name = file.path(out_dir, "querydata_projected_labeled.csv")
)

write.csv(
  query@meta.data,
  file.path(out_dir, "cell_metadata_CITE_DSB_BoneMarrowMap.csv")
)

# ----------------------------
# Composition summaries
# ----------------------------

meta_pass <- query@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  filter(
    mapping_error_QC == "Pass",
    !is.na(predicted_CellType),
    !is.na(sampleID)
  )

celltype_composition <- meta_pass %>%
  count(sampleID, predicted_CellType, name = "n_cells") %>%
  group_by(sampleID) %>%
  mutate(
    total_cells = sum(n_cells),
    fraction = n_cells / total_cells,
    percent = 100 * fraction
  ) %>%
  ungroup() %>%
  arrange(sampleID, desc(percent))

write.csv(
  celltype_composition,
  file.path(out_dir, "07_celltype_composition_by_sampleID.csv"),
  row.names = FALSE
)

p_comp <- celltype_composition %>%
  ggplot(aes(x = sampleID, y = percent, fill = predicted_CellType)) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(x = "Sample ID", y = "Cellular composition (%)", fill = "Predicted cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "08_celltype_composition_by_sampleID_stacked_bar.pdf"),
       p_comp, width = 11, height = 6)

p_counts <- celltype_composition %>%
  ggplot(aes(x = sampleID, y = n_cells, fill = predicted_CellType)) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(x = "Sample ID", y = "Number of cells", fill = "Predicted cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "09_celltype_counts_by_sampleID_stacked_bar.pdf"),
       p_counts, width = 11, height = 6)

# ----------------------------
# Broad / lineage composition
# ----------------------------

# Shared palette, lineage vocabulary and the check_palette() guard.
# Defined once in step 00 so the copies cannot drift apart again.
source("lab_archives/7.Benchmarking-cell_type_annotation/00.celltype_palette.R")

# ----------------------------
# Reference UMAP in the project palette
# ----------------------------

if ("CellType_Broad" %in% colnames(ReferenceSeuratObj@meta.data)) {

  ReferenceSeuratObj$CellType_Broad <- factor(
    ReferenceSeuratObj$CellType_Broad,
    levels = names(celltype_cols)
  )

  pdf(file.path(out_dir, "15_BoneMarrowMap_reference_broad_coloured.pdf"), width = 12, height = 10)
  print(
    DimPlot(
      ReferenceSeuratObj, reduction = "umap", group.by = "CellType_Broad",
      raster = FALSE, label = TRUE, repel = TRUE, label.size = 6,
      cols = celltype_cols
    ) +
      theme_bw()
  )
  dev.off()
}

# ----------------------------
# Save objects and bundle
# ----------------------------

saveRDS(query, seurat_projected_rds)

analysis_bundle <- list(
  seu_input = seu,
  query = query,
  ReferenceSeuratObj = ReferenceSeuratObj,
  ref = ref,
  meta = query@meta.data,
  celltype_composition = celltype_composition,
  composition_broad = composition_broad,
  composition_lineage = composition_lineage,
  celltype_order = celltype_order,
  celltype_cols = celltype_cols,
  lineage_cols = lineage_cols,
  cite_assay_used = cite_assay
)

saveRDS(analysis_bundle, file.path(out_dir, "analysis_bundle.rds"))

cat("\nDone.\n")
cat("Input object:", seurat_input, "\n")
cat("Output directory:", out_dir, "\n")
cat("Projected object:", seurat_projected_rds, "\n\n")

cat("Sample summary after filtering:\n")
print(table(query$sampleID))

cat("\nMapping QC summary:\n")
print(table(query$mapping_error_QC))

cat("\nPredicted lineage summary:\n")
print(table(query$predicted_Lineage, useNA = "ifany"))
