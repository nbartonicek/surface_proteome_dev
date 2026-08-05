#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# LK1 pilot run - step 10 of 19
#
# Load and visualise the BoneMarrowMap Symphony reference, set the uwot path, first projection.
#
# Frozen for the lab archive 2026-07-31 from scripts/backup/10.annotate.R (mtime 2026-05-05).
# md5 of the original: 06d41e155729a041bd60c2e5e7abbf95
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

# Set directory to store projection reference files
projection_path = '../annotation/'

run <- "260423_VH01624_453_222HWMYNX"
sample_short <- "LK1"

projection_path <- "../annotation/"

gex_dir <- file.path(
  "../results/cellranger",
  run,
  "LK1-GEX/outs/filtered_feature_bc_matrix"
)

hto_dir <- file.path(
  "../results/cite_seq_count",
  run,
  "hto_counts_LK1/umi_count"
)

adt_dir <- file.path(
  "../results/cite_seq_count",
  run,
  "adt_counts_LK1/umi_count"
)

out_dir <- file.path("../results/seurat_demux", run)
projection_dir <- file.path(out_dir, "projectionFigures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

seurat_rds <- file.path(out_dir, "LK1_GEX_HTO_demux_seurat.rds")
seurat_projected_rds <- file.path(out_dir, "LK1_GEX_HTO_ADT_BoneMarrowMap_projected.rds")

# ----------------------------
# HTO sample annotation
# ----------------------------

hto_names <- c(
  "HTO1-GTCAACTCTTTAGCG" = "MOLM13",
  "HTO2-TGATGGCCTATTGGG" = "HBDN206_MNpCT",
  "HTO3-TTCCGCCTCTCTTTG" = "HBDN392_AML_MDS",
  "HTO4-AGTAAGTTCAGCGTA" = "HBDN501_AML_KMT2A"
)

# ----------------------------
# Load existing demux object if available,
# otherwise create it from GEX + HTO
# ----------------------------

if (file.exists(seurat_rds)) {
  
  message("Loading existing Seurat object: ", seurat_rds)
  seu <- readRDS(seurat_rds)
  
}

if (!"ADT" %in% Assays(seu)) {
  
  message("Adding ADT assay.")
  
  adt_mat <- readMM(file.path(adt_dir, "matrix.mtx.gz"))
  
  adt_barcodes <- fread(
    file.path(adt_dir, "barcodes.tsv.gz"),
    header = FALSE
  )$V1
  
  adt_features <- fread(
    file.path(adt_dir, "features.tsv.gz"),
    header = FALSE
  )
  
  adt_feature_names <- if (ncol(adt_features) >= 2) {
    adt_features$V2
  } else {
    adt_features$V1
  }
  
  adt_barcodes <- paste0(adt_barcodes, "-1")
  
  rownames(adt_mat) <- make.unique(adt_feature_names)
  colnames(adt_mat) <- adt_barcodes
  
  common_adt <- intersect(colnames(seu), colnames(adt_mat))
  
  message("ADT shared cells: ", length(common_adt))
  
  seu <- subset(seu, cells = common_adt)
  adt_mat <- adt_mat[, common_adt, drop = FALSE]
  
  seu[["ADT"]] <- CreateAssayObject(counts = adt_mat)
  
  DefaultAssay(seu) <- "ADT"
  
  seu <- NormalizeData(
    seu,
    normalization.method = "CLR",
    margin = 2
  )
  
  seu$ADT_total <- Matrix::colSums(
    GetAssayData(seu, assay = "ADT", layer = "counts")
  )
  
  seu$ADT_features <- Matrix::colSums(
    GetAssayData(seu, assay = "ADT", layer = "counts") > 0
  )
  
} else {
  message("ADT assay already present. Skipping ADT import.")
}

seu$sampleID <- seu$sample_name

#curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds', 
#                    destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Download uwot model file - 221 Mb
#curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_uwot_model.uwot', 
#                    destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot'))

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

# Visualize Bone Marrow reference
ReferenceSeuratObj <- create_ReferenceObject(ref)
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', raster=FALSE, label=TRUE, label.size = 4)
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Broad', raster=FALSE, label=TRUE, label.size = 4)

######### load in demultiplexed seurat file

batchvar <- 'sampleID'

# Map query dataset using Symphony 
query <- map_Query(
  query = seu_filter,  # load in query seurat object. Can also load counts and meta-data separately. 
  ref_obj = ref,
  vars = batchvar
)

query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

pdf(file.path(out_dir, "09_mapping_error_QC.pdf"), width = 8, height = 6)
print(plot_MappingErrorQC(query))
dev.off()

query <- predict_CellTypes(
  query_obj = query,
  ref_obj = ref,
  final_label = "predicted_CellType"
)

pdf(file.path(out_dir, "10_projected_predicted_celltypes_pass_only.pdf"), width = 20, height = 12)
print(
  DimPlot(
    subset(query, mapping_error_QC == "Pass"),
    group.by = "predicted_CellType",
    label = TRUE,
    label.size = 4
  )
)
dev.off()

query <- predict_Pseudotime(
  query_obj = query,
  ref_obj = ref,
  final_label = "predicted_Pseudotime"
)

pdf(file.path(out_dir, "11_projected_pseudotime_pass_only.pdf"), width = 8, height = 6)
print(
  FeaturePlot(
    subset(query, mapping_error_QC == "Pass"),
    features = "predicted_Pseudotime"
  )
)
dev.off()

projection_plots <- plot_Projection_byDonor(
  query_obj = query,
  batch_key = "sampleID",
  ref_obj = ref,
  save_folder = projection_dir
)

pdf(file.path(out_dir, "12_projection_by_sampleID.pdf"), width = 14, height = 10)
print(patchwork::wrap_plots(projection_plots))
dev.off()

save_ProjectionResults(
  query_obj = query,
  file_name = file.path(out_dir, "querydata_projected_labeled.csv")
)

write.csv(
  query@meta.data,
  file.path(out_dir, "cell_metadata_demux_ADT_BoneMarrowMap.csv")
)

saveRDS(query, seurat_projected_rds)

cat("\nDone.\n")
cat("Output directory:", out_dir, "\n")
cat("Projected object:", seurat_projected_rds, "\n")
cat("\nDemux summary:\n")
print(demux_summary)
cat("\nSample QC summary:\n")
print(sample_qc)

FeaturePlot(subset(query, mapping_error_QC == 'Pass'), features = c('predicted_Pseudotime'))

# Set batch/condition to be visualized individually
batch_key <- 'sampleID'

# returns a list of plots for each donor from a pre-specified batch variable
projection_plots <- plot_Projection_byDonor(
  query_obj = query, 
  batch_key = batch_key, 
  ref_obj = ref, 
  save_folder = 'projectionFigures/'
)

# show plots together with patchwork
patchwork::wrap_plots(projection_plots)

save_ProjectionResults(
  query_obj = query, 
  file_name = 'querydata_projected_labeled.csv')

meta_comp <- query@meta.data %>%
  filter(mapping_error_QC == "Pass") %>%
  filter(!is.na(predicted_CellType), !is.na(sampleID))

celltype_composition <- meta_comp %>%
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
  file.path(out_dir, "celltype_composition_by_sampleID.csv"),
  row.names = FALSE
)

celltype_composition

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
  file.path(out_dir, "13.celltype_composition_by_sampleID_stacked_bar.pdf"),
  p_comp,
  width = 10,
  height = 6
)

p_comp

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
  file.path(out_dir, "14.celltype_counts_by_sampleID_stacked_bar.pdf"),
  p_counts,
  width = 10,
  height = 6
)

p_counts

######### annotate per cell type
meta <- query@meta.data

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  count(sampleID, predicted_CellType_Broad) %>%
  group_by(sampleID) %>%
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

composition_lineage <- composition_broad %>%
  group_by(sampleID, lineage) %>%
  summarise(percent = sum(percent), .groups = "drop")

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

pdf(paste0(out_dir, "/13.broad_celltype_composition.pdf"),width=8,height=4)

ggplot(
  composition_lineage,
  aes(sampleID, percent, fill = lineage)
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
dev.off()

#######
celltype_order <- rev(c(
  
  # Stem / progenitor
  "HSC MPP",
  "LMPP",
  "MEP",
  "Megakaryocyte Precursor",
  "Early GMP",
  "Late GMP",
  "Cycling Progenitor",
  "EoBasoMast Precursor",
  
  # Myeloid
  "Pro-Monocyte",
  "Monocyte",
  "cDC",
  "pDC",
  
  # Lymphoid
  "Early Lymphoid",
  "Pro-B",
  "Pre-B",
  "B",
  "Naive T",
  "CD4 Memory T",
  "CD8 Memory T",
  "NK",
  "Plasma Cell",
  
  # Erythroid
  "Early Erythroid",
  "Late Erythroid"
))

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
  "Late Erythroid" = "#C51B8A"
)

composition_broad$predicted_CellType_Broad <- factor(
  composition_broad$predicted_CellType_Broad,
  levels = celltype_order
)
pdf(paste0(out_dir, "/14.specific_celltype_composition.pdf"),width=8,height=4)

ggplot(
  composition_broad,
  aes(sampleID, percent, fill = predicted_CellType_Broad)
) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = celltype_cols) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cellular composition (%)",
    fill = "Broad cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
dev.off()
########## annotate per cell type
meta <- query@meta.data

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  count(sampleID, predicted_CellType_Broad) %>%
  group_by(sampleID) %>%
  mutate(percent = n / sum(n) * 100)

composition_plot <- ggplot(
  composition_broad,
  aes(sampleID, percent, fill = predicted_CellType_Broad)
) +
  geom_col(width = 0.85) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cellular composition (%)",
    fill = "Broad cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(out_dir, "14.celltype_composition_by_sampleID_stacked_bar.pdf"),
  composition_plot,
  width = 10,
  height = 6
)



heat_df <- composition_broad %>%
  select(sampleID, predicted_CellType_Broad, percent) %>%
  pivot_wider(
    names_from = predicted_CellType_Broad,
    values_from = percent,
    values_fill = 0
  )

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$sampleID

mat_scaled <- t(scale(t(mat)))

pdf(paste0(out_dir, "/15.heatmap_cell_type_broad.pdf"),width=8,height=3)
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  border_color = NA
)
dev.off()

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  count(sampleID, predicted_CellType) %>%
  group_by(sampleID) %>%
  mutate(percent = n / sum(n) * 100)

heat_df <- composition_broad %>%
  select(sampleID, predicted_CellType, percent) %>%
  pivot_wider(
    names_from = predicted_CellType,
    values_from = percent,
    values_fill = 0
  )

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$sampleID

mat_scaled <- t(scale(t(mat)))

pdf(paste0(out_dir, "/15a.heatmap_cell_type_all.pdf"),width=8,height=3.5)
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  border_color = NA
)
dev.off()


#### plot first heatmap based on groups
ReferenceSeuratObj$CellType_Broad <- factor(
  ReferenceSeuratObj$CellType_Broad,
  levels = names(celltype_cols)
)
pdf(paste0(out_dir, "/10a.database_umap.pdf"),width=12,height=10)
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
dev.off()



#save the objects
analysis_bundle <- list(
  
  # Core objects
  seu = seu,
  query = query,
  ReferenceSeuratObj = ReferenceSeuratObj,
  
  # Reference
  ref = ref,
  
  # Metadata tables
  meta = query@meta.data,
  composition_broad = composition_broad,
  composition_lineage = composition_lineage,
  celltype_composition = celltype_composition,
  
  # Plot ordering / palettes
  celltype_order = celltype_order,
  celltype_cols = celltype_cols,
  lineage_cols = lineage_cols
  
)

saveRDS(
  analysis_bundle,
  file.path(out_dir, "analysis_bundle.rds")
)

