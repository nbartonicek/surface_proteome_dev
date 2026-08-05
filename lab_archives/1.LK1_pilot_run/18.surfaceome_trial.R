# ------------------------------------------------------------------
# LK1 pilot run - step 18 of 19
#
# First surfaceome pass - surfaceome-only UMAP on the Numbat-annotated object. Exploratory, saves nothing.
#
# Frozen for the lab archive 2026-07-31 from scripts/14.surfaceome_trial.R (mtime 2026-05-25).
# md5 of the original: e19c220260e1ad25dece280330c450f6
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

run <- "260423_VH01624_453_222HWMYNX"
sample_name <- "LK1-GEX"
sample_short <- "LK1"

sample_to_use <- "HBDN206-MNpCT"
out_dir <- file.path("../results/data_integration", run)
projection_dir <- file.path(out_dir, "projection_umaps")

numbat_dir <- file.path(
  "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/results/seurat_annotated/260423_VH01624_453_222HWMYNX/numbat",
  paste0("LK1_", sample_to_use),
  "numbat_final"
)
annotation_dir <- file.path("../results/seurat_annotated", run)

copykat_file <- file.path(
  annotation_dir,
  "copykat_annotated",
  sample_to_use,
  "copykat_prediction_with_metadata.csv"
)

out_dir <- file.path(
  projection_dir,
  paste0("numbat_copykat_comparison_", sample_to_use)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD SAMPLE SEURAT OBJECT
# ============================================================

seurat_file <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)
seu <- readRDS(seurat_file)
seu_sub <- subset(seu, subset = sampleID == sample_to_use)

cat("Sample:", sample_to_use, "\n")
cat("Seurat cells:", ncol(seu_sub), "\n")

# ============================================================
# LOAD NUMBAT
# ============================================================

numbat_files <- list.files(
  numbat_dir,
  pattern = "clone_post|cell_anno|numbat.*cell|posterior|subclone",
  full.names = TRUE,
  recursive = TRUE
)

clone_file <- numbat_files[grepl("clone_post", basename(numbat_files))][1]

if (is.na(clone_file)) {
  stop("Could not find clone_post file in: ", numbat_dir)
}

numbat_meta <- fread(clone_file)

clone_col <- dplyr::case_when(
  "clone" %in% colnames(numbat_meta) ~ "clone",
  "clone_opt" %in% colnames(numbat_meta) ~ "clone_opt",
  "clone_post" %in% colnames(numbat_meta) ~ "clone_post",
  "subclone" %in% colnames(numbat_meta) ~ "subclone",
  TRUE ~ NA_character_
)

if (is.na(clone_col)) {
  warning("No clone-like column found in Numbat output. Using Unknown.")
  numbat_meta$numbat_clone <- "Unknown"
} else {
  message("Using Numbat clone column: ", clone_col)
  numbat_meta$numbat_clone <- as.character(numbat_meta[[clone_col]])
}

numbat_meta <- numbat_meta %>%
  mutate(
    numbat_malignancy = case_when(
      numbat_clone %in% c("normal", "diploid", "neutral", "0") ~ "CNV_neutral",
      is.na(numbat_clone) | numbat_clone == "Unknown" ~ "Not_called",
      TRUE ~ "CNV_aberrant"
    )
  ) %>%
  dplyr::select(cell, numbat_clone, numbat_malignancy, everything())
cat("Numbat cells:", nrow(numbat_meta), "\n")
print(table(numbat_meta$numbat_malignancy, useNA = "ifany"))

# match metadata to Seurat cells
numbat_match <- numbat_meta %>%
  dplyr::distinct(cell, .keep_all = TRUE)

# add metadata
seu_sub$numbat_clone <- numbat_match$numbat_clone[
  match(colnames(seu_sub), numbat_match$cell)
]

seu_sub$numbat_malignancy <- numbat_match$numbat_malignancy[
  match(colnames(seu_sub), numbat_match$cell)
]

# replace missing
seu_sub$numbat_clone <- tidyr::replace_na(
  seu_sub$numbat_clone,
  "Not_called"
)

seu_sub$numbat_malignancy <- tidyr::replace_na(
  seu_sub$numbat_malignancy,
  "Not_called"
)

# QC
cat("Matched cells:",
    sum(seu_sub$numbat_malignancy != "Not_called"),
    "\n")

table(seu_sub$numbat_malignancy, useNA = "ifany")

table(seu_sub$numbat_clone, useNA = "ifany")

DefaultAssay(seu_sub) <- "RNA"


###### import surfaceome

surfaceome_file <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/proteome/S2_File.csv"

# seu_sub should already exist and contain:
# - RNA assay
# - umap_projected reduction
# - predicted_CellType_Broad or equivalent annotation
# - Numbat metadata if already added

DefaultAssay(seu_sub) <- "RNA"

surfaceome <- data.table::fread(surfaceome_file) %>%
  as_tibble()

colnames(surfaceome)

surfaceome_genes <- surfaceome %>%
  dplyr::filter(organism == "Human") %>%
  dplyr::pull(`ENTREZ gene symbol`) %>%
  unique() %>%
  na.omit() %>%
  as.character()


surfaceome <- data.table::fread(surfaceome_file) %>%
  as_tibble()

colnames(surfaceome)

surfaceome_genes <- surfaceome %>%
  dplyr::filter(organism == "Human") %>%
  dplyr::pull(`ENTREZ gene symbol`) %>%
  unique() %>%
  na.omit() %>%
  as.character()

surfaceome_genes <- surfaceome_genes[surfaceome_genes != ""]

surfaceome_genes_use <- intersect(surfaceome_genes, rownames(seu_sub))

cat("Surfaceome genes in file:", length(surfaceome_genes), "\n")
cat("Surfaceome genes found in Seurat object:", length(surfaceome_genes_use), "\n")

# Optional: keep only high-confidence CSPA proteins
surfaceome_highconf <- surfaceome %>%
  dplyr::filter(
    organism == "Human",
    stringr::str_detect(`CSPA category`, "high confidence")
  ) %>%
  dplyr::pull(`ENTREZ gene symbol`) %>%
  unique() %>%
  na.omit() %>%
  as.character()

seu_surface <- seu_sub

VariableFeatures(seu_surface) <- surface_features

seu_surface <- ScaleData(
  seu_surface,
  features = surface_features,
  verbose = FALSE
)

seu_surface <- RunPCA(
  seu_surface,
  features = surface_features,
  npcs = 30,
  reduction.name = "surface_pca",
  reduction.key = "surfacePC_",
  verbose = FALSE
)

ElbowPlot(seu_surface, reduction = "surface_pca", ndims = 30)

seu_surface <- RunUMAP(
  seu_surface,
  reduction = "surface_pca",
  dims = 1:20,
  reduction.name = "surface_umap",
  reduction.key = "surfaceUMAP_",
  verbose = FALSE
)

seu_surface <- FindNeighbors(
  seu_surface,
  reduction = "surface_pca",
  dims = 1:20,
  graph.name = "surface_snn",
  verbose = FALSE
)

seu_surface <- FindClusters(
  seu_surface,
  graph.name = "surface_snn",
  resolution = 0.4,
  cluster.name = "surface_clusters",
  verbose = FALSE
)



# ------------------------------------------------------------
# Transfer surfaceome reductions/clusters back to seu_sub
# ------------------------------------------------------------

seu_sub[["surface_pca"]] <- seu_surface[["surface_pca"]]
seu_sub[["surface_umap"]] <- seu_surface[["surface_umap"]]

seu_sub$surface_clusters <- seu_surface$surface_clusters

# ------------------------------------------------------------
# Plot surfaceome-only UMAP
# ------------------------------------------------------------

p_surface_clusters <- DimPlot(
  seu_sub,
  reduction = "surface_umap",
  group.by = "surface_clusters",
  label = TRUE
) +
  ggtitle("Surfaceome-only RNA UMAP / surfaceome clusters")

p_surface_clusters

p_surface_sample <- DimPlot(
  seu_sub,
  reduction = "surface_umap",
  group.by = "sample_name"
) +
  ggtitle("Surfaceome-only RNA UMAP / sample")

p_surface_sample

p_surface_numbat <- DimPlot(
  seu_sub,
  reduction = "surface_umap",
  group.by = "numbat_malignancy",
  cols = c(
    CNV_aberrant = "#E41A1C",
    CNV_neutral = "#377EB8",
    Not_called = "grey80"
  )
) +
  ggtitle("Surfaceome-only RNA UMAP / Numbat malignancy")

p_surface_numbat



# ------------------------------------------------------------
# Use colours from existing umap_projected annotation
# ------------------------------------------------------------

celltype_col <- "predicted_CellType_Broad"

if (!celltype_col %in% colnames(seu_sub@meta.data)) {
  stop("Column not found in metadata: ", celltype_col)
}

celltypes <- sort(unique(as.character(seu_sub@meta.data[[celltype_col]])))

celltype_cols <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(celltypes)),
  celltypes
)

p_projected_celltypes <- DimPlot(
  seu_sub,
  reduction = "umap_projected",
  group.by = celltype_col,
  cols = celltype_cols,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Original projected UMAP / cell type annotation")

p_surface_celltypes <- DimPlot(
  seu_sub,
  reduction = "surface_umap",
  group.by = celltype_col,
  cols = celltype_cols,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Surfaceome-only UMAP / same cell type colours")

p_projected_celltypes
p_surface_celltypes

# ------------------------------------------------------------
# Show how original projected annotations split by surfaceome clusters
# ------------------------------------------------------------

p_split_table <- seu_sub@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::count(
    .data[[celltype_col]],
    surface_clusters,
    name = "n"
  ) %>%
  dplyr::group_by(.data[[celltype_col]]) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

p_split_bar <- ggplot(
  p_split_table,
  aes(
    x = .data[[celltype_col]],
    y = percent,
    fill = surface_clusters
  )
) +
  geom_col(color = "black", linewidth = 0.15) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "How projected cell types split across surfaceome clusters",
    x = "Projected cell type",
    y = "Percent of cells",
    fill = "Surfaceome cluster"
  )

p_split_bar

Idents(seu_sub) <- "surface_clusters"

surface_cluster_markers <- FindAllMarkers(
  seu_sub,
  assay = "RNA",
  features = surface_features,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

surface_cluster_markers <- surface_cluster_markers %>%
  dplyr::left_join(
    surfaceome %>%
      dplyr::filter(organism == "Human") %>%
      dplyr::select(
        gene = `ENTREZ.gene.symbol`,
        CSPA_category = `CSPA.category`,
        UP_Protein_name,
        CD,
        count.detection.in.different.cell.types,
        UniProt.Cell.surface
      ) %>%
      dplyr::distinct(gene, .keep_all = TRUE),
    by = "gene"
  )

surface_cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_min(order_by = p_val_adj, n = 20) %>%
  dplyr::ungroup() %>%
  print(n = 100)

write.csv(
  surface_cluster_markers,
  "surfaceome_cluster_markers.csv",
  row.names = FALSE
)







