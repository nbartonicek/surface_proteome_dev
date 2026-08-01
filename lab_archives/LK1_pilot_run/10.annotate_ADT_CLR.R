#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 10 of 18
#
# Project the query onto BoneMarrowMap using the ADT/CLR object. Writes seurat_annotated/.
#
# Frozen for the lab archive 2026-07-31 from scripts/backup/10a.annotate.R (mtime 2026-05-06).
# md5 of the original: 8e2bc2377d3152df370da10966e69fb9
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
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

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"
sample_short <- "LK1"

projection_path <- "../annotation/"

cite_dsb_dir <- file.path(
  "../results/cite_qc_then_dsb",
  run,
  sample_name
)

seurat_input <- file.path(
  cite_dsb_dir,
  "seurat_demux_ADT_CLR_DSB_QC.rds"
)

out_dir <- file.path("../results/seurat_annotated", run)
projection_dir <- file.path(out_dir, "projectionFigures/")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

seurat_projected_rds <- file.path(
  out_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

# ----------------------------
# Load CITE/DSB Seurat object
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

# Clean sample labels
seu$sampleID <- as.character(seu$sample_name)

# Remove unwanted demux labels if present
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
# Load BoneMarrowMap / Symphony reference
# ----------------------------

ref_file <- file.path(projection_path, "BoneMarrowMap_SymphonyReference.rds")
uwot_file <- file.path(projection_path, "BoneMarrowMap_uwot_model.uwot")

if (!file.exists(ref_file)) {
  stop("Missing BoneMarrowMap reference: ", ref_file)
}

if (!file.exists(uwot_file)) {
  stop("Missing BoneMarrowMap uwot model: ", uwot_file)
}

ref <- readRDS(ref_file)
ref$save_uwot_path <- uwot_file

ReferenceSeuratObj <- create_ReferenceObject(ref)

pdf(file.path(out_dir, "00_BoneMarrowMap_reference_celltype.pdf"), width = 12, height = 10)
print(
  DimPlot(
    ReferenceSeuratObj,
    reduction = "umap",
    group.by = "CellType_Annotation_formatted",
    raster = FALSE,
    label = TRUE,
    repel = TRUE,
    label.size = 4
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
    raster = FALSE,
    label = TRUE,
    repel = TRUE,
    label.size = 4
  ) +
    ggtitle("BoneMarrowMap reference: broad cell type")
)
dev.off()

# ----------------------------
# Map query using Symphony
# ----------------------------

DefaultAssay(seu) <- "RNA"

batchvar <- "sampleID"

cat("Running Symphony map_Query using batch variable:", batchvar, "\n")

query <- map_Query(
  query = seu,
  ref_obj = ref,
  vars = batchvar
)

query <- calculate_MappingError(
  query,
  reference = ref,
  MAD_threshold = 2.5
)

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

# ----------------------------
# Add broad cell type if possible
# ----------------------------
# Some BoneMarrowMap versions add this automatically; if not, derive from metadata where possible.

if (!"predicted_CellType_Broad" %in% colnames(query@meta.data)) {
  if ("predicted_CellType_Broad" %in% colnames(query@meta.data)) {
    message("predicted_CellType_Broad already present.")
  } else {
    message("predicted_CellType_Broad not found; deriving broad lineage manually.")
  }
}

# Manual lineage grouping based on predicted_CellType / predicted_CellType_Broad
query$predicted_Lineage <- case_when(
  query$predicted_CellType %in% c(
    "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
    "Cycling Progenitor", "EoBasoMast Precursor",
    "Megakaryocyte Precursor"
  ) ~ "Stem / progenitor",
  
  query$predicted_CellType %in% c(
    "Monocyte", "Pro-Monocyte", "cDC", "pDC"
  ) ~ "Myeloid / DC",
  
  query$predicted_CellType %in% c(
    "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
    "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
  ) ~ "Lymphoid",
  
  query$predicted_CellType %in% c(
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
    query_pass,
    reduction = "umap",
    group.by = "predicted_CellType",
    label = TRUE,
    repel = TRUE,
    label.size = 4
  ) +
    ggtitle("BoneMarrowMap projected cell types: mapping pass only")
)
dev.off()

pdf(file.path(projection_dir, "03_projected_samples_pass_only.pdf"), width = 10, height = 7)
print(
  DimPlot(
    query_pass,
    reduction = "umap",
    group.by = "sampleID",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Samples on BoneMarrowMap projection: mapping pass only")
)
dev.off()

pdf(file.path(projection_dir, "04_projected_pseudotime_pass_only.pdf"), width = 8, height = 6)
print(
  FeaturePlot(
    query_pass,
    reduction = "umap",
    features = "predicted_Pseudotime",
    order = TRUE
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
# CITE marker projection sanity checks
# ----------------------------

cite_assay <- case_when(
  #"CITE_DSB" %in% Assays(query) ~ "CITE_DSB",
  "ADT_CLR" %in% Assays(query) ~ "ADT_CLR",
  "ADT" %in% Assays(query) ~ "ADT",
  TRUE ~ NA_character_
)

if (!is.na(cite_assay)) {
  DefaultAssay(query) <- cite_assay
  
  cite_mat <- GetAssayData(query, assay = cite_assay, layer = "data")
  
  cite_markers <- rownames(cite_mat)
  cite_markers <- cite_markers[
    !grepl("isotype|igg|control|unmapped", cite_markers, ignore.case = TRUE)
  ]
  
  marker_means <- Matrix::rowMeans(as.matrix(cite_mat[cite_markers, , drop = FALSE]))
  
  top_cite_markers <- names(sort(marker_means, decreasing = TRUE))[1:min(16, length(marker_means))]
  
  pdf(file.path(projection_dir, paste0("06_", cite_assay, "_top_marker_projection.pdf")), width = 14, height = 10)
  print(
    FeaturePlot(
      query,
      reduction = "umap_projected",
      features = top_cite_markers,
      ncol = 4,
      order = TRUE
    )
  )
  dev.off()
}

extra_markers <- c(
  "CD33",
  "CD133",
  "CD47",
  "Flt3-Flk2",   # correct (not FLT3)
  "TACTILE",
  "GPR56",
  "CD7"
)

marker_map <- c(
  "CD33" = "CD33",
  "CD133" = "CD133",
  "CD47" = "CD47",
  "CD135 (FLT3)" = "Flt3-Flk2",
  "GPR56" = "GPR56",
  "CD7" = "CD7",
  "CD96" = "TACTILE"
)
marker_map <- marker_map[marker_map %in% colnames(adt_df)]

plot_df$marker <- factor(
  plot_df$marker,
  levels = marker_map,
  labels = names(marker_map)
)

extra_markers %in% rownames(cite_mat)

emb <- Embeddings(query, "umap_projected") %>%
  as.data.frame() %>%
  rownames_to_column("cell")

colnames(emb)[2:3] <- c("UMAP_1", "UMAP_2")

adt_df <- FetchData(
  query,
  vars = extra_markers,
  values = "data"
) %>%
  rownames_to_column("cell")

plot_df <- emb %>%
  left_join(adt_df, by = "cell") %>%
  pivot_longer(
    cols = all_of(extra_markers),
    names_to = "marker",
    values_to = "value"
  ) %>%
  group_by(marker) %>%
  mutate(
    value_clip = pmin(
      pmax(value, quantile(value, 0.05, na.rm = TRUE)),
      quantile(value, 0.95, na.rm = TRUE)
    )
  ) %>%
  ungroup()

p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = value_clip)) +
  geom_point(size = 0.15) +
  facet_wrap(~ marker, ncol = 4) +
  scale_colour_gradientn(
    colours = c("grey90", "#FDB863", "#E66101", "#B2182B")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(colour = "CLR", x = "Projected UMAP 1", y = "Projected UMAP 2")

ggsave(
  file.path(projection_dir, paste0("07_", cite_assay, "_selected_markers_manual_projectedUMAP.pdf")),
  p,
  width = 14,
  height = 8
)


DefaultAssay(query) <- "RNA"

# ----------------------------
# Save projection result tables
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
  labs(
    x = "Sample ID",
    y = "Cellular composition (%)",
    fill = "Predicted cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "08_celltype_composition_by_sampleID_stacked_bar.pdf"),
  p_comp,
  width = 11,
  height = 6
)

p_counts <- celltype_composition %>%
  ggplot(aes(x = sampleID, y = n_cells, fill = predicted_CellType)) +
  geom_col(width = 0.8) +
  theme_bw() +
  labs(
    x = "Sample ID",
    y = "Number of cells",
    fill = "Predicted cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "09_celltype_counts_by_sampleID_stacked_bar.pdf"),
  p_counts,
  width = 11,
  height = 6
)

# ----------------------------
# Broad / lineage composition
# ----------------------------

celltype_order <- rev(c(
  "HSC MPP",
  "LMPP",
  "MEP",
  "Megakaryocyte Precursor",
  "GMP",
  "Early GMP",
  "Late GMP",
  "Cycling Progenitor",
  "EoBasoMast Precursor",
  "Pro-Monocyte",
  "Monocyte",
  "cDC",
  "pDC",
  "Early Lymphoid",
  "Pro-B",
  "Pre-B",
  "B",
  "Naive T",
  "CD4 Memory T",
  "CD8 Memory T",
  "NK",
  "Plasma Cell",
  "Early Erythroid",
  "Late Erythroid"
))

celltype_cols <- c(
  "HSC MPP" = "#1B9E77",
  "LMPP" = "#66C2A5",
  "MEP" = "#B2DF8A",
  "GMP" = "#33A02C",
  "Early GMP" = "#A6D854",
  "Late GMP" = "#006D2C",
  "Cycling Progenitor" = "#00441B",
  "EoBasoMast Precursor" = "#8DD3C7",
  "Megakaryocyte Precursor" = "#4DAF4A",
  "Monocyte" = "#E31A1C",
  "Pro-Monocyte" = "#FB6A4A",
  "cDC" = "#FD8D3C",
  "pDC" = "#FCBBA1",
  "Naive T" = "#2171B5",
  "CD4 Memory T" = "#6BAED6",
  "CD8 Memory T" = "#08519C",
  "NK" = "#54278F",
  "Early Lymphoid" = "#9E9AC8",
  "B" = "#3182BD",
  "Pre-B" = "#9ECAE1",
  "Pro-B" = "#C6DBEF",
  "Plasma Cell" = "#756BB1",
  "Early Erythroid" = "#F768A1",
  "Late Erythroid" = "#C51B8A"
)

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

composition_lineage <- meta_pass %>%
  count(sampleID, predicted_Lineage, name = "n_cells") %>%
  group_by(sampleID) %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(
  composition_lineage,
  file.path(out_dir, "10_lineage_composition_by_sampleID.csv"),
  row.names = FALSE
)

p_lineage <- ggplot(
  composition_lineage,
  aes(sampleID, percent, fill = predicted_Lineage)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = lineage_cols) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cellular composition (%)",
    fill = "Lineage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "11_lineage_composition_by_sampleID.pdf"),
  p_lineage,
  width = 8,
  height = 4
)

composition_broad <- meta_pass %>%
  count(sampleID, predicted_CellType, name = "n_cells") %>%
  group_by(sampleID) %>%
  mutate(percent = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    predicted_CellType = factor(predicted_CellType, levels = celltype_order)
  )

write.csv(
  composition_broad,
  file.path(out_dir, "12_specific_celltype_composition_by_sampleID.csv"),
  row.names = FALSE
)

p_broad <- ggplot(
  composition_broad,
  aes(sampleID, percent, fill = predicted_CellType)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = celltype_cols, na.value = "grey80") +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cellular composition (%)",
    fill = "Predicted cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "13_specific_celltype_composition_by_sampleID.pdf"),
  p_broad,
  width = 10,
  height = 5
)

# ----------------------------
# Heatmaps
# ----------------------------

heat_df <- composition_broad %>%
  dplyr::select(sampleID, predicted_CellType, percent) %>%
  pivot_wider(
    names_from = predicted_CellType,
    values_from = percent,
    values_fill = 0
  )

mat <- as.matrix(heat_df[, -1, drop = FALSE])
rownames(mat) <- heat_df$sampleID

# Avoid NaNs from zero-variance rows
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

pdf(file.path(out_dir, "14_heatmap_cell_type_composition_scaled.pdf"), width = 10, height = 4)
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  border_color = NA,
  main = "Scaled cell-type composition by sample"
)
dev.off()

mat_unscaled <- mat

pdf(file.path(out_dir, "14a_heatmap_cell_type_composition_percent.pdf"), width = 10, height = 4)
pheatmap(
  mat_unscaled,
  cluster_rows = FALSE,
  border_color = NA,
  main = "Cell-type composition (%) by sample"
)
dev.off()

# ----------------------------
# Reference UMAP with consistent colours
# ----------------------------

if ("CellType_Broad" %in% colnames(ReferenceSeuratObj@meta.data)) {
  ReferenceSeuratObj$CellType_Broad <- factor(
    ReferenceSeuratObj$CellType_Broad,
    levels = names(celltype_cols)
  )
  
  pdf(file.path(out_dir, "15_BoneMarrowMap_reference_broad_coloured.pdf"), width = 12, height = 10)
  print(
    DimPlot(
      ReferenceSeuratObj,
      reduction = "umap",
      group.by = "CellType_Broad",
      raster = FALSE,
      label = TRUE,
      repel = TRUE,
      label.size = 6,
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

saveRDS(
  analysis_bundle,
  file.path(out_dir, "analysis_bundle.rds")
)

cat("\nDone.\n")
cat("Input object:", seurat_input, "\n")
cat("Output directory:", out_dir, "\n")
cat("Projected object:", seurat_projected_rds, "\n")
cat("Analysis bundle:", file.path(out_dir, "analysis_bundle.rds"), "\n\n")

cat("Sample summary after filtering:\n")
print(table(query$sampleID))

cat("\nMapping QC summary:\n")
print(table(query$mapping_error_QC))

cat("\nPredicted lineage summary:\n")
print(table(query$predicted_Lineage, useNA = "ifany"))


######### annotate numbat
RUN <- "260423_VH01624_453_222HWMYNX"
SAMPLE_SHORT <- "LK1"
PROJ <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

seu_file <- file.path(
  PROJ,
  "results/seurat_annotated",
  RUN,
  paste0(SAMPLE_SHORT, "_seurat_annotated.rds")
)

OUT_BASE <- file.path(
  PROJ,
  "results/seurat_annotated",
  RUN,
  "numbat",
  "numbat_inputs_no_seurat"
)

dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(seu) <- "RNA"

# ----------------------------
# Pick donor/sample column
# ----------------------------
sample_col <- if ("sample_name" %in% colnames(seu@meta.data)) {
  "sample_name"
} else if ("orig.ident" %in% colnames(seu@meta.data)) {
  "orig.ident"
} else {
  stop("No sample_name or orig.ident column found.")
}

samples <- sort(unique(na.omit(seu@meta.data[[sample_col]])))

# ----------------------------
# Export per donor/sample
# ----------------------------
for (donor in samples) {
  
  message("Exporting: ", donor)
  
  label <- paste0(SAMPLE_SHORT, "_", donor)
  out_dir <- file.path(OUT_BASE, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  seu_sub <- subset(
    seu,
    cells = colnames(seu)[seu@meta.data[[sample_col]] == donor]
  )
  
  expr <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")
  expr <- as(expr, "dgCMatrix")
  
  # ----------------------------
  # Cell type annotation
  # ----------------------------
  if (!"predicted_CellType_Broad" %in% colnames(seu_sub@meta.data)) {
    stop("predicted_CellType_Broad column not found.")
  }
  
  cell_type <- as.character(seu_sub$predicted_CellType_Broad)
  cell_type[is.na(cell_type) | cell_type == ""] <- "unknown"
  
  cell_annot <- data.table(
    cell = colnames(seu_sub),
    sample = label,
    clone = "unknown",
    cell_type = cell_type
  )
  
  # ----------------------------
  # Write files
  # ----------------------------
  Matrix::writeMM(
    expr,
    file.path(out_dir, paste0(label, "_counts.mtx"))
  )
  
  fwrite(
    data.table(gene = rownames(expr)),
    file.path(out_dir, paste0(label, "_genes.tsv")),
    sep = "\t",
    col.names = FALSE
  )
  
  fwrite(
    data.table(cell = colnames(expr)),
    file.path(out_dir, paste0(label, "_barcodes.tsv")),
    sep = "\t",
    col.names = FALSE
  )
  
  fwrite(
    cell_annot,
    file.path(out_dir, paste0(label, "_cell_annot.tsv")),
    sep = "\t"
  )
  
  message("  cells: ", ncol(expr))
  message("  genes: ", nrow(expr))
  message("  cell types:")
  print(table(cell_annot$cell_type, useNA = "ifany"))
}

message("Done. Files written to: ", OUT_BASE)



