#!/usr/bin/env Rscript

# Investigates cells that BoneMarrowMap couldn't confidently annotate in the
# LK2-GEX object from the newer results_nf pipeline checkpoint, and shows
# where they fall relative to annotated cells on both the object's own
# ("not-projected") UMAP and the BoneMarrowMap reference-projected UMAP.
#
# On inspection this object already carries predicted_CellType_Broad (and
# mapping_error_QC), but the 548 NA cells in predicted_CellType_Broad are
# exactly the mapping_error_QC == "Fail" cells (1:1, confirmed directly on
# this object) - i.e. BoneMarrowMap's own mapping-error QC step is what
# already produced the "unannotated" cells, they weren't just missed.
# However this checkpoint does NOT have a umap_projected reduction or the
# fine-grained predicted_CellType column, meaning map_Query() was either
# never run to completion on this exact object or its output wasn't kept at
# this pipeline stage. Since a "projected" UMAP was explicitly asked for,
# map_Query()/predict_CellTypes() has to run regardless of the existing
# partial annotation - re-running it also keeps the annotation and the
# projected coordinates internally consistent (same single mapping call),
# rather than bolting a fresh projection onto a possibly differently-sourced
# old annotation. This will overwrite predicted_CellType_Broad with a
# freshly (re-)computed version.
#
# Same BoneMarrowMap call sequence as 8.annotate.R (map_Query ->
# calculate_MappingError -> predict_CellTypes), reused rather than guessed.

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

# ----------------------------
# Load object
# ----------------------------

stopifnot(file.exists(input_rds))
message("Loading: ", input_rds)
seu <- readRDS(input_rds)

message("Cells: ", ncol(seu))
message("Existing reductions: ", paste(Reductions(seu), collapse = ", "))
message(
  "Has predicted_CellType_Broad: ",
  "predicted_CellType_Broad" %in% colnames(seu@meta.data)
)

# ----------------------------
# Annotate with BoneMarrowMap only if projected coordinates are missing
# (see header comment for why predicted_CellType_Broad alone isn't
# sufficient to skip this step here)
# ----------------------------

needs_annotation <- !("umap_projected" %in% Reductions(seu))

#if (!needs_annotation) {
#
#  message("umap_projected already present - reusing existing annotation, skipping BoneMarrowMap.")
#
#} else if (file.exists(annotated_rds)) {
#
#  message("Loading previously-computed annotated object: ", annotated_rds)
#  #seu <- readRDS(annotated_rds)
#
#} else {

  message("Running BoneMarrowMap projection (umap_projected missing)...")

  ref <- readRDS(paste0(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
  ref$save_uwot_path <- paste0(projection_path, "BoneMarrowMap_uwot_model.uwot")

  stopifnot("sample_name" %in% colnames(seu@meta.data))
  batchvar <- "sample_name"

  seu <- map_Query(
    query = seu,
    ref_obj = ref,
    vars = batchvar
  )

  seu <- seu %>% calculate_MappingError(., reference = ref, MAD_threshold = 4)

  seu <- predict_CellTypes(
    query_obj = seu,
    ref_obj = ref,
    final_label = "predicted_CellType"
  )

  saveRDS(seu, annotated_rds)
  message("Saved annotated object: ", annotated_rds)
#}

if (!all(c("predicted_CellType_Broad", "mapping_error_QC") %in% colnames(seu@meta.data))) {
  stop("predicted_CellType_Broad / mapping_error_QC missing after annotation step - cannot continue.")
}

# ----------------------------
# Derive a display label that folds NA/failed-QC cells into "Unknown"
# ----------------------------

seu$celltype_display <- as.character(seu$predicted_CellType_Broad)
seu$celltype_display[is.na(seu$celltype_display) | seu$celltype_display == ""] <- "Unknown"

n_unknown <- sum(seu$celltype_display == "Unknown")
message("Unknown/unannotated cells: ", n_unknown, " / ", ncol(seu),
        " (", round(100 * n_unknown / ncol(seu), 1), "%)")

# ----------------------------
# Where do the Unknown cells sit relative to the object's own clustering?
# (This is the direct, non-visual answer to "which cluster do they belong to")
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

  write_csv(
    cluster_summary,
    file.path(out_dir, "unknown_cells_by_seurat_cluster.csv")
  )

  print(
    cluster_summary %>%
      filter(is_unknown) %>%
      arrange(desc(percent))
  )

  p_cluster_bar <- ggplot(
    cluster_summary,
    aes(x = factor(seurat_clusters), y = percent, fill = is_unknown)
  ) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                       labels = c("Annotated", "Unknown"), name = NULL) +
    theme_bw(base_size = 12) +
    labs(
      title = "Unknown/unannotated cell fraction per Seurat cluster (native clustering)",
      x = "seurat_clusters", y = "Percent of cluster"
    )

  ggsave(
    file.path(out_dir, "01_unknown_fraction_by_cluster.pdf"),
    p_cluster_bar, width = 9, height = 5
  )
  print(p_cluster_bar)

} else {
  message("No seurat_clusters column found - skipping per-cluster breakdown.")
}

# ----------------------------
# What's actually different about the Unknown cells? map_Query() embeds
# every cell into reference space regardless of QC outcome - it's
# calculate_MappingError()/predict_CellTypes() downstream that gate label
# assignment. So Unknown cells still have real coordinates; compare their
# QC covariates against annotated cells to see what's driving the mapping
# failure (low depth, doublets, HTO ambiguity, patient skew).
# ----------------------------

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
# Nearest-neighbor identity: for each Unknown cell, KNN-classify against
# the annotated cells' embedding coordinates (majority vote, k=15). This
# gives Unknown cells a "what would this cell be if forced to pick" label
# instead of leaving them as an undifferentiated blob - directly answers
# "group them by the cluster they belong to". Run once per available
# embedding (native umap always available; umap_projected once annotation
# has produced it).
# ----------------------------

assign_nearest_known <- function(seu, reduction, k = 15) {

  if (!(reduction %in% Reductions(seu))) {
    message("Reduction '", reduction, "' not present - skipping nearest-neighbor assignment.")
    return(seu)
  }

  emb <- Embeddings(seu, reduction)
  known_idx <- seu$celltype_display != "Unknown"
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
  seu@meta.data[[col_name]][known_idx] <- seu$celltype_display[known_idx]
  seu@meta.data[[col_name]][unknown_idx] <- as.character(pred)

  nn_summary <- tibble(
    nearest_known_celltype = as.character(pred),
    sample_name = if ("sample_name" %in% colnames(seu@meta.data)) seu$sample_name[unknown_idx] else NA_character_
  ) %>%
    count(nearest_known_celltype, sort = TRUE) %>%
    mutate(percent = round(100 * n / sum(n), 1))

  write_csv(
    nn_summary,
    file.path(out_dir, paste0("unknown_nearest_neighbor_identity_", reduction, ".csv"))
  )

  message("\nNearest-neighbor identity of Unknown cells (", reduction, ", k=", k, "):")
  print(nn_summary)

  p_nn <- ggplot(nn_summary, aes(x = reorder(nearest_known_celltype, percent), y = percent)) +
    geom_col(fill = "firebrick", width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Unknown cells' nearest annotated neighbor (", reduction, ", k=", k, ")"),
      x = NULL, y = "Percent of Unknown cells"
    )

  ggsave(
    file.path(out_dir, paste0("04_unknown_nearest_neighbor_identity_", reduction, ".pdf")),
    p_nn, width = 8, height = 6
  )
  print(p_nn)

  seu
}

seu <- assign_nearest_known(seu, "umap")
seu <- assign_nearest_known(seu, "umap_projected")

# ----------------------------
# Colour palette: same bespoke predicted_CellType_Broad palette used in
# 8.annotate.R (celltype_cols), so cell types read the same colour across
# every annotation script in this project. "Unknown" uses the same grey70
# convention as 8.annotate.R's lineage_cols "Other" entry. Any cell type
# present in this object but not in the bespoke palette (e.g. "Stromal",
# which appears here but isn't in 8.annotate.R's list) gets a visibly
# distinct fallback colour rather than silently reusing/clashing with an
# existing one.
# ----------------------------

celltype_cols <- c(
  # Stem / progenitor - greens
  "HSC MPP" = "#1B9E77",
  "LMPP" = "#66C2A5",
  "MEP" = "#B2DF8A",
  "GMP" = "#33A02C",
  "Early GMP" = "#A6D854",
  "Late GMP" = "#006D2C",
  "Cycling Progenitor" = "#00441B",
  "EoBasoMast Precursor" = "#8DD3C7",
  "Megakaryocyte Precursor" = "#4DAF4A",

  # Myeloid / DC - oranges/reds
  "Monocyte" = "#E31A1C",
  "Pro-Monocyte" = "#FB6A4A",
  "cDC" = "#FD8D3C",
  "pDC" = "#FCBBA1",

  # Lymphoid - blues/purples
  "Naive T" = "#2171B5",
  "CD4 Memory T" = "#6BAED6",
  "CD8 Memory T" = "#08519C",
  "NK" = "#54278F",
  "Early Lymphoid" = "#9E9AC8",
  "B" = "#3182BD",
  "Pre-B" = "#9ECAE1",
  "Pro-B" = "#C6DBEF",
  "Plasma Cell" = "#756BB1",

  # Erythroid - pinks
  "Early Erythroid" = "#F768A1",
  "Late Erythroid" = "#C51B8A",

  # Unknown / QC-fail - same grey70 convention as 8.annotate.R's
  # lineage_cols "Other" entry
  "Unknown" = "grey70"
)

known_types <- unique(seu$celltype_display)
missing_from_palette <- setdiff(known_types, names(celltype_cols))

if (length(missing_from_palette) > 0) {
  message(
    "Cell type(s) present in this object but not in the bespoke celltype_cols palette: ",
    paste(missing_from_palette, collapse = ", "),
    " - assigning fallback colours so they're still visibly distinct."
  )
  fallback_cols <- setNames(
    scales::hue_pal()(length(missing_from_palette)),
    missing_from_palette
  )
  celltype_cols <- c(celltype_cols, fallback_cols)
}

type_cols <- celltype_cols[known_types]

# ----------------------------
# UMAP 1: not-projected (object's own native umap/pca-derived clustering)
# ----------------------------

if ("umap" %in% Reductions(seu)) {

  p_native <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "celltype_display",
    cols = type_cols,
    raster = FALSE
  ) +
    ggtitle("Native (not-projected) UMAP - BoneMarrowMap annotation incl. Unknown")

  ggsave(file.path(out_dir, "02_umap_native_with_unknown.pdf"), p_native, width = 10, height = 8)
  print(p_native)

  p_native_highlight <- DimPlot(
    seu,
    reduction = "umap",
    cells.highlight = list(Unknown = colnames(seu)[seu$celltype_display == "Unknown"]),
    cols.highlight = "firebrick",
    cols = "grey85",
    sizes.highlight = 0.6,
    raster = FALSE
  ) +
    ggtitle("Native UMAP - Unknown cells highlighted") +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "02a_umap_native_unknown_highlighted.pdf"), p_native_highlight, width = 9, height = 8)
  print(p_native_highlight)

} else {
  message("No native 'umap' reduction found - skipping not-projected plot.")
}

# ----------------------------
# UMAP 2: BoneMarrowMap-projected
# ----------------------------

if ("umap_projected" %in% Reductions(seu)) {

  p_projected <- DimPlot(
    seu,
    reduction = "umap_projected",
    group.by = "celltype_display",
    cols = type_cols,
    raster = FALSE
  ) +
    ggtitle("BoneMarrowMap-projected UMAP - annotation incl. Unknown")

  ggsave(file.path(out_dir, "03_umap_projected_with_unknown.pdf"), p_projected, width = 10, height = 8)
  print(p_projected)

  p_projected_highlight <- DimPlot(
    seu,
    reduction = "umap_projected",
    cells.highlight = list(Unknown = colnames(seu)[seu$celltype_display == "Unknown"]),
    cols.highlight = "firebrick",
    cols = "grey85",
    sizes.highlight = 0.6,
    raster = FALSE
  ) +
    ggtitle("Projected UMAP - Unknown cells highlighted") +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "03a_umap_projected_unknown_highlighted.pdf"), p_projected_highlight, width = 9, height = 8)
  print(p_projected_highlight)

} else {
  message("No 'umap_projected' reduction found after annotation step - something went wrong upstream.")
}

# ----------------------------
# Save final metadata + object
# ----------------------------

write_csv(
  seu@meta.data %>% as_tibble(rownames = "cell"),
  file.path(out_dir, "cell_metadata_with_unknown_flag.csv")
)

saveRDS(seu, annotated_rds)

message("\nDone.")
message("Output directory: ", out_dir)
message("Annotated object: ", annotated_rds)
