#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 03 of 05
#
# What are the cells BoneMarrowMap refuses to label, and where do they sit.
#
# Collated from 8a.annotate_and_plot_unannotated.R, run on the LK2-GEX object
# taken from the results_nf pipeline checkpoint rather than from step 02.
#
# The starting observation was that the checkpoint object already carries
# predicted_CellType_Broad and mapping_error_QC, and that its NA cells in
# predicted_CellType_Broad are exactly the mapping_error_QC == "Fail" cells,
# 1:1. So the "unannotated" cells are not cells the reference missed - they are
# cells its own mapping-error QC threw out.
#
# The checkpoint does not carry a umap_projected reduction or the fine-grained
# predicted_CellType column, so map_Query() has to be re-run regardless. It is
# re-run rather than bolted onto the existing labels, so that the projected
# coordinates and the labels come out of the same single mapping call. This
# overwrites predicted_CellType_Broad with a freshly computed version.
#
# Same call sequence as steps 01 and 02 (map_Query -> calculate_MappingError ->
# predict_CellTypes), but at MAD_threshold = 4 rather than 2.5. Step 04
# compares what the two thresholds cost.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(class)
  library(tidyverse)
  library(BoneMarrowMap)
  library(symphony)
})

# ----------------------------
# Paths
# ----------------------------

projection_path <- "../annotation/"

input_rds <- "../results_nf/260528_VH01624_464_222K7VKNX/rds/07_cite_qc/LK2-GEX/rds/seurat_demux_ADT_CLR_DSB_QC.rds"

out_dir <- "../results/seurat_annotated/260528_VH01624_464_222K7VKNX/LK2-GEX_results_nf_unannotated_check"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annotated_rds <- file.path(out_dir, "seurat_demux_ADT_CLR_DSB_QC_BoneMarrowMap_annotated.rds")

MAD_THRESHOLD <- 4

# ----------------------------
# Load
# ----------------------------

stopifnot(file.exists(input_rds))
message("Loading: ", input_rds)
seu <- readRDS(input_rds)

message("Cells: ", ncol(seu))
message("Existing reductions: ", paste(Reductions(seu), collapse = ", "))
message("Has predicted_CellType_Broad: ",
        "predicted_CellType_Broad" %in% colnames(seu@meta.data))

# ----------------------------
# Re-annotate
# ----------------------------
# The input checkpoint carries predicted_CellType_Broad but no projected
# coordinates, so the projection has to run (see header). Once it has run once
# the result is on disk, and re-running only to redraw a figure would cost the
# whole mapping again - so a complete saved object is reused.

if (file.exists(annotated_rds)) {

  message("Loading previously annotated object: ", annotated_rds)
  seu_done <- readRDS(annotated_rds)

  if ("umap_projected" %in% Reductions(seu_done) &&
      "predicted_CellType" %in% colnames(seu_done@meta.data)) {
    message("Reusing it - projection is complete.")
    seu <- seu_done
  } else {
    message("Saved object is incomplete - re-running the projection.")
    seu_done <- NULL
  }
  rm(seu_done)
}

if (!("umap_projected" %in% Reductions(seu))) {

  message("Running BoneMarrowMap projection (umap_projected missing)...")

  ref <- readRDS(paste0(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
  ref$save_uwot_path <- paste0(projection_path, "BoneMarrowMap_uwot_model.uwot")

  stopifnot("sample_name" %in% colnames(seu@meta.data))
  batchvar <- "sample_name"

  seu <- map_Query(query = seu, ref_obj = ref, vars = batchvar)
  seu <- calculate_MappingError(seu, reference = ref, MAD_threshold = MAD_THRESHOLD)
  seu <- predict_CellTypes(query_obj = seu, ref_obj = ref, final_label = "predicted_CellType")

  saveRDS(seu, annotated_rds)
  message("Saved annotated object: ", annotated_rds)
}

if (!all(c("predicted_CellType_Broad", "mapping_error_QC") %in% colnames(seu@meta.data))) {
  stop("predicted_CellType_Broad / mapping_error_QC missing after annotation step.")
}

# ----------------------------
# Fold the QC-failed cells into an explicit "Unknown" label
# ----------------------------

seu$celltype_display <- as.character(seu$predicted_CellType_Broad)
seu$celltype_display[is.na(seu$celltype_display) | seu$celltype_display == ""] <- "Unknown"

n_unknown <- sum(seu$celltype_display == "Unknown")
message("Unknown/unannotated cells: ", n_unknown, " / ", ncol(seu),
        " (", round(100 * n_unknown / ncol(seu), 1), "%)")

# ----------------------------
# Which native clusters do they sit in
# ----------------------------

if ("seurat_clusters" %in% colnames(seu@meta.data)) {

  cluster_summary <- seu@meta.data %>%
    as_tibble() %>%
    mutate(is_unknown = celltype_display == "Unknown") %>%
    count(seurat_clusters, is_unknown) %>%
    group_by(seurat_clusters) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup() %>%
    arrange(desc(is_unknown), desc(percent))

  write_csv(cluster_summary, file.path(out_dir, "unknown_cells_by_seurat_cluster.csv"))
  print(cluster_summary %>% filter(is_unknown) %>% arrange(desc(percent)))

  p_cluster_bar <- ggplot(cluster_summary,
                          aes(x = factor(seurat_clusters), y = percent, fill = is_unknown)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                      labels = c("Annotated", "Unknown"), name = NULL) +
    theme_bw(base_size = 12) +
    labs(title = "Unknown/unannotated cell fraction per Seurat cluster (native clustering)",
         x = "seurat_clusters", y = "Percent of cluster")

  ggsave(file.path(out_dir, "01_unknown_fraction_by_cluster.pdf"),
         p_cluster_bar, width = 9, height = 5)

} else {
  message("No seurat_clusters column found - skipping per-cluster breakdown.")
}

# ----------------------------
# What is different about them
# ----------------------------
# map_Query() embeds every cell into reference space regardless of QC outcome;
# it is calculate_MappingError() that gates label assignment. So the Unknown
# cells still have real coordinates, and their QC covariates can be compared
# against the annotated cells to see what drives the mapping failure.

qc_cols <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "scDblFinder.score"),
  colnames(seu@meta.data)
)

qc_comparison <- seu@meta.data %>%
  as_tibble() %>%
  mutate(is_unknown = celltype_display == "Unknown") %>%
  group_by(is_unknown) %>%
  summarise(n = n(), across(all_of(qc_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")

write_csv(qc_comparison, file.path(out_dir, "unknown_vs_annotated_qc_medians.csv"))
message("\nQC covariate medians, Unknown vs annotated:")
print(qc_comparison)

if ("sample_name" %in% colnames(seu@meta.data)) {
  patient_summary <- seu@meta.data %>%
    as_tibble() %>%
    mutate(is_unknown = celltype_display == "Unknown") %>%
    count(sample_name, is_unknown) %>%
    group_by(sample_name) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()

  write_csv(patient_summary, file.path(out_dir, "unknown_fraction_by_patient.csv"))
  message("\nUnknown fraction by patient:")
  print(patient_summary %>% filter(is_unknown) %>% arrange(desc(percent)))
}

if ("HTO_classification.global" %in% colnames(seu@meta.data)) {
  hto_summary <- seu@meta.data %>%
    as_tibble() %>%
    mutate(is_unknown = celltype_display == "Unknown") %>%
    count(is_unknown, HTO_classification.global) %>%
    group_by(is_unknown) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()

  write_csv(hto_summary, file.path(out_dir, "unknown_vs_annotated_HTO_classification.csv"))
  message("\nHTO classification, Unknown vs annotated:")
  print(hto_summary)
}

# ----------------------------
# What would they be if forced to pick
# ----------------------------
# KNN-classify each Unknown cell against the annotated cells' coordinates,
# majority vote, k = 15. Run once per embedding: the object's own UMAP and the
# reference-projected one. If the two agree the Unknown cells have a real
# identity; if they disagree the cells are being placed by depth rather than by
# biology.

assign_nearest_known <- function(seu, reduction, k = 15) {

  if (!(reduction %in% Reductions(seu))) {
    message("Reduction '", reduction, "' not present - skipping.")
    return(seu)
  }

  emb         <- Embeddings(seu, reduction)
  known_idx   <- seu$celltype_display != "Unknown"
  unknown_idx <- !known_idx

  if (sum(unknown_idx) == 0) {
    message("No Unknown cells to classify for reduction '", reduction, "'.")
    return(seu)
  }

  pred <- class::knn(
    train = emb[known_idx, , drop = FALSE],
    test  = emb[unknown_idx, , drop = FALSE],
    cl    = seu$celltype_display[known_idx],
    k     = k
  )

  col_name <- paste0("nearest_known_celltype_", reduction)
  seu@meta.data[[col_name]] <- NA_character_
  seu@meta.data[[col_name]][known_idx]   <- seu$celltype_display[known_idx]
  seu@meta.data[[col_name]][unknown_idx] <- as.character(pred)

  nn_summary <- tibble(nearest_known_celltype = as.character(pred)) %>%
    count(nearest_known_celltype, sort = TRUE) %>%
    mutate(percent = round(100 * n / sum(n), 1))

  write_csv(nn_summary,
            file.path(out_dir, paste0("unknown_nearest_neighbor_identity_", reduction, ".csv")))

  message("\nNearest-neighbor identity of Unknown cells (", reduction, ", k=", k, "):")
  print(nn_summary)

  p_nn <- ggplot(nn_summary, aes(x = reorder(nearest_known_celltype, percent), y = percent)) +
    geom_col(fill = "firebrick", width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(title = paste0("Unknown cells' nearest annotated neighbor (", reduction, ", k=", k, ")"),
         x = NULL, y = "Percent of Unknown cells")

  ggsave(file.path(out_dir, paste0("04_unknown_nearest_neighbor_identity_", reduction, ".pdf")),
         p_nn, width = 8, height = 6)

  seu
}

seu <- assign_nearest_known(seu, "umap")
seu <- assign_nearest_known(seu, "umap_projected")

# ----------------------------
# Palette
# ----------------------------
# Same bespoke broad-cell-type palette as steps 01 and 02, so cell types read
# the same colour across every annotation script in the project. "Unknown" uses
# the grey70 that "Other" uses in the lineage palette. Anything in this object
# but not in the palette (e.g. "Stromal") gets a visibly distinct fallback
# rather than silently reusing a colour.

# Shared palette, lineage vocabulary and the check_palette() guard.
# Defined once in step 00 so the copies cannot drift apart again.
source("lab_archives/7.Benchmarking-cell_type_annotation/00.celltype_palette.R")

celltype_cols <- c(celltype_cols, "Unknown" = "grey70")

known_types <- unique(seu$celltype_display)
missing_from_palette <- setdiff(known_types, names(celltype_cols))

if (length(missing_from_palette) > 0) {
  message("Cell type(s) present in this object but not in the bespoke palette: ",
          paste(missing_from_palette, collapse = ", "),
          " - assigning fallback colours.")
  celltype_cols <- c(
    celltype_cols,
    setNames(scales::hue_pal()(length(missing_from_palette)), missing_from_palette)
  )
}

type_cols <- celltype_cols[known_types]

# ----------------------------
# UMAP 1: the object's own (not-projected) embedding
# ----------------------------

if ("umap" %in% Reductions(seu)) {

  p_native <- DimPlot(seu, reduction = "umap", group.by = "celltype_display",
                      cols = type_cols, raster = FALSE) +
    ggtitle("Native (not-projected) UMAP - BoneMarrowMap annotation incl. Unknown")

  ggsave(file.path(out_dir, "02_umap_native_with_unknown.pdf"), p_native, width = 10, height = 8)

  p_native_highlight <- DimPlot(
    seu, reduction = "umap",
    cells.highlight = list(Unknown = colnames(seu)[seu$celltype_display == "Unknown"]),
    cols.highlight = "firebrick", cols = "grey85",
    sizes.highlight = 0.6, raster = FALSE
  ) +
    ggtitle("Native UMAP - Unknown cells highlighted") +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "02a_umap_native_unknown_highlighted.pdf"),
         p_native_highlight, width = 9, height = 8)

} else {
  message("No native 'umap' reduction found - skipping not-projected plot.")
}

# ----------------------------
# UMAP 2: reference-projected
# ----------------------------

if ("umap_projected" %in% Reductions(seu)) {

  p_projected <- DimPlot(seu, reduction = "umap_projected", group.by = "celltype_display",
                         cols = type_cols, raster = FALSE) +
    ggtitle("BoneMarrowMap-projected UMAP - annotation incl. Unknown")

  ggsave(file.path(out_dir, "03_umap_projected_with_unknown.pdf"),
         p_projected, width = 10, height = 8)

  p_projected_highlight <- DimPlot(
    seu, reduction = "umap_projected",
    cells.highlight = list(Unknown = colnames(seu)[seu$celltype_display == "Unknown"]),
    cols.highlight = "firebrick", cols = "grey85",
    sizes.highlight = 0.6, raster = FALSE
  ) +
    ggtitle("Projected UMAP - Unknown cells highlighted") +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "03a_umap_projected_unknown_highlighted.pdf"),
         p_projected_highlight, width = 9, height = 8)

} else {
  message("No 'umap_projected' reduction found after annotation - something went wrong upstream.")
}

# ----------------------------
# Save
# ----------------------------

write_csv(
  seu@meta.data %>% as_tibble(rownames = "cell"),
  file.path(out_dir, "cell_metadata_with_unknown_flag.csv")
)

saveRDS(seu, annotated_rds)

message("\nDone.")
message("Output directory: ", out_dir)
message("Annotated object: ", annotated_rds)
