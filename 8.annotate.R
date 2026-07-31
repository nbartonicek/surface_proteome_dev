#!/usr/bin/env Rscript

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

run <- "260528_VH01624_464_222K7VKNX"
sample_name <- "LK2-GEX"
sample_short <- "LK2"

projection_path <- "../annotation/"

in_dir <- file.path(
  "../results/cite_seq_dsb",
  run,
  sample_name
)

out_dir <- file.path("../results/seurat_annotated", run)
projection_dir <- file.path(out_dir, "projectionFigures/")
barcode_dir=file.path("../results/seurat_demux",run,"/numbat/barcodes_by_sample_name/")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barcode_dir,recursive = TRUE)
seurat_rds <- file.path(in_dir, paste0("seurat_demux_ADT_CLR_DSB_QC.rds"))
seurat_projected_rds <- file.path(out_dir, "demux_singlets_annotated_seurat.rds")

out_numbat <- file.path(
  "../results/seurat_annotated",
  run,
  "numbat",
  "numbat_inputs_no_seurat"
)

# ----------------------------
# Load existing demux object if available,
# otherwise create it from GEX + HTO
# ----------------------------

if (file.exists(seurat_rds)) {
  
  message("Loading existing Seurat object: ", seurat_rds)
  seu <- readRDS(seurat_rds)
  
}

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

batchvar <- 'sample_name'

# Map query dataset using Symphony 
query <- map_Query(
  query = seu,  # load in query seurat object. Can also load counts and meta-data separately. 
  ref_obj = ref,
  vars = batchvar
)

query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

pdf(file.path(out_dir, "01_mapping_error_QC.pdf"), width = 8, height = 6)
print(plot_MappingErrorQC(query))
dev.off()

query <- predict_CellTypes(
  query_obj = query,
  ref_obj = ref,
  final_label = "predicted_CellType"
)

pdf(file.path(out_dir, "02_projected_predicted_celltypes_pass_only.pdf"), width = 20, height = 12)
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

pdf(file.path(out_dir, "03_projected_pseudotime_pass_only.pdf"), width = 8, height = 6)
print(
  FeaturePlot(
    subset(query, mapping_error_QC == "Pass"),
    features = "predicted_Pseudotime"
  )
)
dev.off()

sample_ids <- levels(factor(query$sample_name))

projection_plots <- lapply(sample_ids, function(sid) {
  
  q_sub <- subset(query, subset = sample_name == sid)
  
  out_i <- file.path(projection_dir, sid)
  dir.create(out_i, recursive = TRUE, showWarnings = FALSE)
  
  plot_Projection_byDonor(
    query_obj = q_sub,
    batch_key = "sample_name",
    ref_obj = ref,
    save_folder = out_i
  )[[1]] +
    ggtitle(sid)
})

names(projection_plots) <- sample_ids

pdf(
  file.path(projection_dir, "04_projection_by_sample_name.pdf"),
  width = 14,
  height = 10
)

print(
  patchwork::wrap_plots(projection_plots, ncol = 2)
)

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
batch_key <- 'sample_name'

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
  filter(!is.na(predicted_CellType), !is.na(sample_name))

celltype_composition <- meta_comp %>%
  count(sample_name, predicted_CellType, name = "n_cells") %>%
  group_by(sample_name) %>%
  mutate(
    total_cells = sum(n_cells),
    fraction = n_cells / total_cells,
    percent = 100 * fraction
  ) %>%
  ungroup() %>%
  arrange(sample_name, desc(percent))

write.csv(
  celltype_composition,
  file.path(out_dir, "celltype_composition_by_sample_name.csv"),
  row.names = FALSE
)

celltype_composition

p_comp <- celltype_composition %>%
  ggplot(aes(x = sample_name, y = percent, fill = predicted_CellType)) +
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
  file.path(out_dir, "05.celltype_composition_by_sample_name_stacked_bar.pdf"),
  p_comp,
  width = 10,
  height = 6
)

p_comp

p_counts <- celltype_composition %>%
  ggplot(aes(x = sample_name, y = n_cells, fill = predicted_CellType)) +
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
  file.path(out_dir, "06.celltype_counts_by_sample_name_stacked_bar.pdf"),
  p_counts,
  width = 10,
  height = 6
)

p_counts

######### annotate per cell type
meta <- query@meta.data

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
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

composition_lineage <- composition_broad %>%
  group_by(sample_name, lineage) %>%
  summarise(percent = sum(percent), .groups = "drop")

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

pdf(paste0(out_dir, "/07.broad_celltype_composition.pdf"),width=8,height=4)

ggplot(
  composition_lineage,
  aes(sample_name, percent, fill = lineage)
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
pdf(paste0(out_dir, "/08.specific_celltype_composition.pdf"),width=8,height=4)

ggplot(
  composition_broad,
  aes(sample_name, percent, fill = predicted_CellType_Broad)
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
  count(sample_name, predicted_CellType_Broad) %>%
  group_by(sample_name) %>%
  mutate(percent = n / sum(n) * 100)

composition_plot <- ggplot(
  composition_broad,
  aes(sample_name, percent, fill = predicted_CellType_Broad)
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
  file.path(out_dir, "09.celltype_composition_by_sample_name_stacked_bar.pdf"),
  composition_plot,
  width = 10,
  height = 6
)



heat_df <- composition_broad %>%
  select(sample_name, predicted_CellType_Broad, percent) %>%
  pivot_wider(
    names_from = predicted_CellType_Broad,
    values_from = percent,
    values_fill = 0
  )

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$sample_name

mat_scaled <- t(scale(t(mat)))

pdf(paste0(out_dir, "/10.heatmap_cell_type_broad.pdf"),width=8,height=3)
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  border_color = NA
)
dev.off()

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  count(sample_name, predicted_CellType) %>%
  group_by(sample_name) %>%
  mutate(percent = n / sum(n) * 100)

heat_df <- composition_broad %>%
  select(sample_name, predicted_CellType, percent) %>%
  pivot_wider(
    names_from = predicted_CellType,
    values_from = percent,
    values_fill = 0
  )

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$sample_name

mat_scaled <- t(scale(t(mat)))

pdf(paste0(out_dir, "/10a.heatmap_cell_type_all.pdf"),width=8,height=3.5)
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
pdf(paste0(out_dir, "/11.database_umap.pdf"),width=12,height=10)
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

for(sample_name_i in unique(seu$sample_name)){
  cat(sample_name_i)
  barcodes<-colnames(seu)[seu$sample_name==sample_name_i]
  write.table(as.data.frame(barcodes),file=paste0(barcode_dir,sample_name_i,"_barcodes.tsv"),col.names = F,quote=F,row.names=F)
}


###### annotate numbat
# ----------------------------
# Export per donor/sample
# ----------------------------
seu <- query
for (donor in unique(seu$sample_name)) {
  
  message("Exporting: ", donor)
  
  label <- paste0(sample_short, "_", donor)
  out_dir <- file.path(out_numbat, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  seu_sub <- subset(
    seu,
    cells = colnames(seu)[seu@meta.data$sample_name == donor]
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

message("Done. Files written to: ", out_numbat)
