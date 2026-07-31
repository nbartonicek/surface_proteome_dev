# ============================================================
# LK1 + LK2 CD14 monocyte analysis
# Uses Numbat-integrated Seurat objects directly
# Outputs:
#   1. merged CD14 object
#   2. RNA UMAP / Harmony UMAP
#   3. CITE_DSB UMAP
#   4. WNN RNA+CITE UMAP
#   5. CD14 CITE_DSB AML-vs-normal marker table
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(knitr)
})

set.seed(123)

# ----------------------------
# Settings
# ----------------------------

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

out_dir <- file.path("../merged_run_1_2", "DE", "CD14_monocyte_DE_Rmd_numbat")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

assay_use <- "RNA"
cite_assay <- "CITE_DSB"

sample_col <- "sample_name"
fine_col <- "predicted_CellType"
broad_col <- "predicted_CellType_Broad"
doublet_col <- "scDblFinder.class"

normal_sample <- "normal-01"

min_cd14_cells <- 50
min_cells_per_group <- 10

fdr_cutoff <- 0.05
logfc_cutoff <- 1
bcv <- 0.4

sample_cols <- c(
  "APOP576-TP53"      = "#66C2A5",
  "HBDN206-MNpCT"     = "#FC8D62",
  "HBDN376-TP53"      = "#8DA0CB",
  "HBDN392-AML-MDS"   = "#E78AC3",
  "HBDN498-TP53"      = "#FFD92F",
  "HBDN501-AML-KMT2A" = "#A6D854",
  "normal-01"         = "darkred"
)

aml_status_cols <- c(
  "AML"    = "#D95F02",
  "Normal" = "#1F78B4"
)

# ----------------------------
# Input Numbat-integrated objects
# ----------------------------

numbat_files <- tibble(
  run_id = c("LK1", "LK2"),
  rds = c(
    file.path(
      proj,
      "results/seurat_annotated",
      "260423_VH01624_453_222HWMYNX",
      "numbat",
      "LK1_projected_CITE_DSB_Numbat_integrated.rds"
    ),
    file.path(
      proj,
      "results/seurat_annotated",
      "260528_VH01624_464_222K7VKNX",
      "numbat",
      "LK2_projected_CITE_DSB_Numbat_integrated.rds"
    )
  )
)

missing_files <- numbat_files$rds[!file.exists(numbat_files$rds)]

if (length(missing_files) > 0) {
  stop(
    "Missing Numbat RDS file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

required_cols <- c(
  sample_col,
  fine_col,
  broad_col,
  doublet_col,
  "numbat_compartment",
  "numbat_call",
  "numbat_clone"
)

# ----------------------------
# Helper functions
# ----------------------------

get_assay_matrix <- function(obj, assay_name) {
  
  mat <- tryCatch(
    GetAssayData(obj, assay = assay_name, layer = "data"),
    error = function(e) {
      GetAssayData(obj, assay = assay_name, slot = "data")
    }
  )
  
  if (nrow(mat) == 0) {
    mat <- tryCatch(
      GetAssayData(obj, assay = assay_name, layer = "counts"),
      error = function(e) {
        GetAssayData(obj, assay = assay_name, slot = "counts")
      }
    )
  }
  
  mat
}

make_aml_status <- function(obj) {
  dplyr::case_when(
    obj@meta.data[[sample_col]] == normal_sample ~ "Normal",
    tolower(as.character(obj$numbat_compartment)) == "normal" ~ "Normal",
    tolower(as.character(obj$numbat_compartment)) == "tumor" ~ "AML",
    TRUE ~ NA_character_
  )
}

plot_volcano <- function(
    top,
    FDR_cutoff = 0.05,
    logFC_cutoff = 1,
    max.overlaps = 40,
    title = NULL,
    subtitle = NULL
) {
  
  top <- as.data.frame(top)
  
  if (!"gene" %in% colnames(top)) {
    top$gene <- rownames(top)
  }
  
  top <- top %>%
    mutate(
      diff_expression = case_when(
        logFC >= logFC_cutoff & FDR <= FDR_cutoff ~ "UP",
        logFC <= -logFC_cutoff & FDR <= FDR_cutoff ~ "DOWN",
        TRUE ~ "NO"
      ),
      labelgenes = if_else(diff_expression != "NO", gene, "")
    )
  
  ggplot(
    top,
    aes(
      x = logFC,
      y = -log10(FDR),
      colour = diff_expression,
      label = labelgenes
    )
  ) +
    geom_point(alpha = 0.8, size = 1.3) +
    ggrepel::geom_text_repel(
      max.overlaps = max.overlaps,
      min.segment.length = 5,
      size = 3
    ) +
    scale_color_manual(
      values = c(
        "DOWN" = "darkblue",
        "NO" = "black",
        "UP" = "darkred"
      )
    ) +
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      colour = "darkred"
    ) +
    geom_hline(
      yintercept = -log10(FDR_cutoff),
      linetype = "dashed",
      colour = "darkred"
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "log2FC AML / Normal",
      y = "-log10(FDR)",
      colour = "DE"
    )
}

run_edgeR_exact <- function(count_mat, group_vec, bcv = 0.4) {
  
  group_vec <- factor(group_vec, levels = c("Normal", "AML"))
  
  dge <- edgeR::DGEList(counts = count_mat, group = group_vec)
  
  keep <- edgeR::filterByExpr(dge, group = group_vec)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  dge <- edgeR::calcNormFactors(dge)
  
  et <- edgeR::exactTest(
    dge,
    pair = c("Normal", "AML"),
    dispersion = bcv^2
  )
  
  edgeR::topTags(et, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    as_tibble()
}

# ----------------------------
# Load objects
# ----------------------------

seu_list <- purrr::map2(
  numbat_files$rds,
  numbat_files$run_id,
  function(rds, run_id) {
    
    message("Loading ", run_id, ": ", rds)
    
    obj <- readRDS(rds)
    obj$run_id <- run_id
    
    missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
    
    if (length(missing_cols) > 0) {
      stop(
        "Missing metadata in ", run_id, ":\n",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    obj$AML_status_numbat <- make_aml_status(obj)
    
    cat("\n", run_id, " Numbat compartment:\n", sep = "")
    print(table(obj$numbat_compartment, useNA = "ifany"))
    
    cat("\n", run_id, " AML_status_numbat:\n", sep = "")
    print(table(obj$AML_status_numbat, useNA = "ifany"))
    
    obj
  }
)

names(seu_list) <- numbat_files$run_id

# ----------------------------
# Merge LK1 + LK2
# ----------------------------

seu <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "LK1_LK2_Numbat"
)

if (inherits(seu[[assay_use]], "Assay5")) {
  seu <- JoinLayers(seu, assay = assay_use)
}

if (cite_assay %in% Assays(seu) && inherits(seu[[cite_assay]], "Assay5")) {
  seu <- JoinLayers(seu, assay = cite_assay)
}

cat("\nMerged object:\n")
print(seu)

cat("\nCells by run and AML status:\n")
print(table(seu$run_id, seu$AML_status_numbat, useNA = "ifany"))

# ----------------------------
# Subset CD14 monocytes
# ----------------------------

seu_cd14_all <- subset(
  seu,
  subset =
    scDblFinder.class == "singlet" &
    predicted_CellType == "CD14 Mono" &
    !grepl("MOLM13", sample_name, ignore.case = TRUE) &
    !is.na(AML_status_numbat)
)

cd14_counts <- seu_cd14_all@meta.data %>%
  rownames_to_column("cell") %>%
  count(run_id, sample_name, AML_status_numbat, name = "n_cd14_cells") %>%
  arrange(run_id, sample_name)

write.csv(
  cd14_counts,
  file.path(out_dir, "CD14_counts_before_filtering.csv"),
  row.names = FALSE
)

print(cd14_counts)

good_samples <- cd14_counts %>%
  group_by(sample_name) %>%
  summarise(n_total = sum(n_cd14_cells), .groups = "drop") %>%
  filter(n_total >= min_cd14_cells) %>%
  pull(sample_name)

seu_cd14 <- subset(
  seu_cd14_all,
  subset = sample_name %in% good_samples
)

cat("\nCD14 cells after filtering:\n")
print(table(seu_cd14$sample_name, seu_cd14$AML_status_numbat))

write.csv(
  seu_cd14@meta.data,
  file.path(out_dir, "CD14_metadata_with_numbat_status.csv")
)

# ----------------------------
# RNA PCA / UMAP / Harmony
# ----------------------------

DefaultAssay(seu_cd14) <- assay_use

seu_cd14 <- NormalizeData(seu_cd14)
seu_cd14 <- FindVariableFeatures(seu_cd14, nfeatures = 3000)
seu_cd14 <- ScaleData(seu_cd14, verbose = FALSE)
seu_cd14 <- RunPCA(seu_cd14, npcs = 50, verbose = FALSE)

seu_cd14 <- RunUMAP(
  seu_cd14,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap_pca",
  reduction.key = "pcaUMAP_"
)

p_raw_sample <- DimPlot(
  seu_cd14,
  reduction = "umap_pca",
  group.by = sample_col,
  raster = FALSE
) +
  scale_color_manual(values = sample_cols, na.value = "grey80") +
  ggtitle("CD14 monocytes: non-integrated RNA UMAP by sample")

p_raw_status <- DimPlot(
  seu_cd14,
  reduction = "umap_pca",
  group.by = "AML_status_numbat",
  raster = FALSE
) +
  scale_color_manual(values = aml_status_cols, na.value = "grey80") +
  ggtitle("CD14 monocytes: non-integrated RNA UMAP by AML/Numbat status")

ggsave(
  file.path(plot_dir, "CD14_RNA_nonintegrated_UMAP_sample_status.pdf"),
  p_raw_sample + p_raw_status,
  width = 14,
  height = 6
)

seu_cd14 <- RunHarmony(
  object = seu_cd14,
  group.by.vars = sample_col,
  reduction = "pca",
  dims.use = 1:30,
  theta = 1,
  lambda = 2,
  reduction.save = "harmony",
  verbose = FALSE
)

seu_cd14 <- RunUMAP(
  seu_cd14,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap_harmony",
  reduction.key = "harmonyUMAP_"
)

p_harmony_sample <- DimPlot(
  seu_cd14,
  reduction = "umap_harmony",
  group.by = sample_col,
  raster = FALSE
) +
  scale_color_manual(values = sample_cols, na.value = "grey80") +
  ggtitle("CD14 monocytes: Harmony RNA UMAP by sample")

p_harmony_status <- DimPlot(
  seu_cd14,
  reduction = "umap_harmony",
  group.by = "AML_status_numbat",
  raster = FALSE
) +
  scale_color_manual(values = aml_status_cols, na.value = "grey80") +
  ggtitle("CD14 monocytes: Harmony RNA UMAP by AML/Numbat status")

ggsave(
  file.path(plot_dir, "CD14_RNA_Harmony_UMAP_sample_status.pdf"),
  p_harmony_sample + p_harmony_status,
  width = 14,
  height = 6
)

# ----------------------------
# CITE_DSB-only UMAP
# ----------------------------

if (cite_assay %in% Assays(seu_cd14)) {
  
  DefaultAssay(seu_cd14) <- cite_assay
  
  cite_features <- rownames(seu_cd14[[cite_assay]])
  
  seu_cd14 <- ScaleData(
    seu_cd14,
    assay = cite_assay,
    features = cite_features,
    verbose = FALSE
  )
  
  seu_cd14 <- RunPCA(
    seu_cd14,
    assay = cite_assay,
    features = cite_features,
    reduction.name = "apca",
    reduction.key = "apca_",
    npcs = min(30, length(cite_features)),
    verbose = FALSE
  )
  
  cite_dims <- 1:min(20, ncol(Embeddings(seu_cd14, "apca")))
  
  seu_cd14 <- FindNeighbors(
    seu_cd14,
    reduction = "apca",
    dims = cite_dims,
    graph.name = "CITE_DSB_snn"
  )
  
  seu_cd14 <- FindClusters(
    seu_cd14,
    graph.name = "CITE_DSB_snn",
    resolution = 0.3,
    cluster.name = "CITE_DSB_cluster"
  )
  
  seu_cd14 <- RunUMAP(
    seu_cd14,
    reduction = "apca",
    dims = cite_dims,
    reduction.name = "umap_CITE_DSB",
    reduction.key = "citeUMAP_"
  )
  
  p_cite_cluster <- DimPlot(
    seu_cd14,
    reduction = "umap_CITE_DSB",
    group.by = "CITE_DSB_cluster",
    label = TRUE,
    raster = FALSE
  ) +
    ggtitle("CD14 monocytes: CITE_DSB-only clusters")
  
  p_cite_sample <- DimPlot(
    seu_cd14,
    reduction = "umap_CITE_DSB",
    group.by = sample_col,
    raster = FALSE
  ) +
    scale_color_manual(values = sample_cols, na.value = "grey80") +
    ggtitle("CITE_DSB UMAP by sample")
  
  p_cite_status <- DimPlot(
    seu_cd14,
    reduction = "umap_CITE_DSB",
    group.by = "AML_status_numbat",
    raster = FALSE
  ) +
    scale_color_manual(values = aml_status_cols, na.value = "grey80") +
    ggtitle("CITE_DSB UMAP by AML/Numbat status")
  
  ggsave(
    file.path(plot_dir, "CD14_CITE_DSB_only_clustering.pdf"),
    p_cite_cluster + p_cite_sample + p_cite_status,
    width = 15,
    height = 5
  )
}

# ----------------------------
# WNN RNA + CITE_DSB UMAP
# ----------------------------

if (cite_assay %in% Assays(seu_cd14) && "apca" %in% Reductions(seu_cd14)) {
  
  DefaultAssay(seu_cd14) <- assay_use
  
  rna_dims <- 1:30
  cite_dims <- 1:min(20, ncol(Embeddings(seu_cd14, "apca")))
  
  seu_cd14 <- FindMultiModalNeighbors(
    seu_cd14,
    reduction.list = list("pca", "apca"),
    dims.list = list(rna_dims, cite_dims),
    modality.weight.name = c("RNA.weight", "CITE_DSB.weight")
  )
  
  seu_cd14 <- RunUMAP(
    seu_cd14,
    nn.name = "weighted.nn",
    reduction.name = "umap_WNN_RNA_CITE",
    reduction.key = "wnnUMAP_"
  )
  
  seu_cd14 <- FindClusters(
    seu_cd14,
    graph.name = "wsnn",
    resolution = 0.3,
    cluster.name = "WNN_RNA_CITE_cluster"
  )
  
  p_wnn_cluster <- DimPlot(
    seu_cd14,
    reduction = "umap_WNN_RNA_CITE",
    group.by = "WNN_RNA_CITE_cluster",
    label = TRUE,
    raster = FALSE
  ) +
    ggtitle("CD14 monocytes: RNA + CITE_DSB WNN clusters")
  
  p_wnn_sample <- DimPlot(
    seu_cd14,
    reduction = "umap_WNN_RNA_CITE",
    group.by = sample_col,
    raster = FALSE
  ) +
    scale_color_manual(values = sample_cols, na.value = "grey80") +
    ggtitle("WNN UMAP by sample")
  
  p_wnn_status <- DimPlot(
    seu_cd14,
    reduction = "umap_WNN_RNA_CITE",
    group.by = "AML_status_numbat",
    raster = FALSE
  ) +
    scale_color_manual(values = aml_status_cols, na.value = "grey80") +
    ggtitle("WNN UMAP by AML/Numbat status")
  
  ggsave(
    file.path(plot_dir, "CD14_RNA_CITE_DSB_WNN_clustering.pdf"),
    p_wnn_cluster + p_wnn_sample + p_wnn_status,
    width = 15,
    height = 5
  )
}

# ----------------------------
# RNA marker visualisation
# ----------------------------

DefaultAssay(seu_cd14) <- assay_use

rna_markers <- c(
  "LTF", "MPO", "ELANE", "AZU1", "CTSG", "LCN2", "FCGR3B",
  "S100A8", "S100A9", "FCN1", "LYZ", "CD14", "MS4A7", "VCAN"
)

rna_markers <- intersect(rna_markers, rownames(seu_cd14))

if (length(rna_markers) > 0) {
  
  pdf(
    file.path(plot_dir, "CD14_RNA_neutrophil_monocyte_marker_UMAPs.pdf"),
    width = 10,
    height = 8
  )
  
  for (gene in rna_markers) {
    print(
      FeaturePlot(
        seu_cd14,
        reduction = "umap_harmony",
        features = gene,
        order = TRUE,
        pt.size = 0.25
      ) +
        ggtitle(paste0("RNA: ", gene))
    )
  }
  
  dev.off()
}

# ----------------------------
# CITE marker visualisation
# ----------------------------

if (cite_assay %in% Assays(seu_cd14)) {
  
  DefaultAssay(seu_cd14) <- cite_assay
  
  cite_markers <- c(
    "BAFF-R", "CD54", "CD69", "CD9", "CD35",
    "NCAM", "CD14", "CD16", "CD11c", "CD163",
    "CD33", "CD64", "HLA-DR", "HLA-DR-DP-DQ"
  )
  
  cite_markers <- intersect(cite_markers, rownames(seu_cd14[[cite_assay]]))
  
  if (length(cite_markers) > 0) {
    
    pdf(
      file.path(plot_dir, "CD14_CITE_DSB_key_marker_UMAPs.pdf"),
      width = 10,
      height = 8
    )
    
    for (mk in cite_markers) {
      print(
        FeaturePlot(
          seu_cd14,
          reduction = "umap_CITE_DSB",
          features = mk,
          order = TRUE,
          pt.size = 0.25
        ) +
          ggtitle(paste0("CITE_DSB: ", mk))
      )
    }
    
    dev.off()
    
    p_dot <- DotPlot(
      seu_cd14,
      features = cite_markers,
      group.by = "AML_status_numbat"
    ) +
      RotatedAxis() +
      ggtitle("CD14 monocyte CITE_DSB markers by AML/Numbat status")
    
    ggsave(
      file.path(plot_dir, "CD14_CITE_DSB_marker_DotPlot_AML_status.pdf"),
      p_dot,
      width = 10,
      height = 5
    )
  }
}

# ----------------------------
# Pseudobulk RNA DE: AML vs Normal CD14 monocytes
# ----------------------------

DefaultAssay(seu_cd14) <- assay_use

count_mat <- tryCatch(
  GetAssayData(seu_cd14, assay = assay_use, layer = "counts"),
  error = function(e) {
    GetAssayData(seu_cd14, assay = assay_use, slot = "counts")
  }
)

meta_cd14 <- seu_cd14@meta.data %>%
  rownames_to_column("cell") %>%
  filter(AML_status_numbat %in% c("AML", "Normal"))

pb_samples <- meta_cd14 %>%
  distinct(sample_name, AML_status_numbat)

pb_counts <- purrr::map_dfc(pb_samples$sample_name, function(smp) {
  
  cells <- meta_cd14 %>%
    filter(sample_name == smp) %>%
    pull(cell)
  
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
})

pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- pb_samples$sample_name
rownames(pb_counts) <- rownames(count_mat)

group_vec <- pb_samples$AML_status_numbat

if (sum(group_vec == "AML") >= 1 && sum(group_vec == "Normal") >= 1) {
  
  de_res <- run_edgeR_exact(
    count_mat = pb_counts,
    group_vec = group_vec,
    bcv = bcv
  )
  
  write.csv(
    de_res,
    file.path(out_dir, "CD14_monocyte_pseudobulk_edgeR_AML_vs_Normal.csv"),
    row.names = FALSE
  )
  
  p_volcano <- plot_volcano(
    de_res,
    FDR_cutoff = fdr_cutoff,
    logFC_cutoff = logfc_cutoff,
    max.overlaps = 50,
    title = "CD14 monocytes: AML vs Normal",
    subtitle = "Pseudobulk by sample; edgeR exactTest"
  )
  
  ggsave(
    file.path(plot_dir, "CD14_monocyte_pseudobulk_edgeR_volcano.pdf"),
    p_volcano,
    width = 10,
    height = 8
  )
}

# ----------------------------
# CITE_DSB AML-vs-normal within CD14 monocytes
# ----------------------------

if (cite_assay %in% Assays(seu_cd14)) {
  
  DefaultAssay(seu_cd14) <- cite_assay
  
  cite_mat <- get_assay_matrix(seu_cd14, cite_assay)
  cite_features <- rownames(cite_mat)
  
  status_ct <- as.character(seu_cd14$AML_status_numbat)
  keep <- !is.na(status_ct) & status_ct %in% c("AML", "Normal")
  
  cells_use <- colnames(seu_cd14)[keep]
  status_use <- status_ct[keep]
  
  n_aml <- sum(status_use == "AML", na.rm = TRUE)
  n_norm <- sum(status_use == "Normal", na.rm = TRUE)
  
  if (n_aml >= min_cells_per_group && n_norm >= min_cells_per_group) {
    
    mat_use <- as.matrix(cite_mat[, cells_use, drop = FALSE])
    
    cite_de_cd14 <- purrr::map_dfr(cite_features, function(marker) {
      
      x_aml <- mat_use[marker, status_use == "AML"]
      x_norm <- mat_use[marker, status_use == "Normal"]
      
      wt <- tryCatch(
        wilcox.test(x_aml, x_norm),
        error = function(e) NULL
      )
      
      tibble(
        marker = marker,
        n_AML = length(x_aml),
        n_Normal = length(x_norm),
        mean_AML = mean(x_aml, na.rm = TRUE),
        mean_Normal = mean(x_norm, na.rm = TRUE),
        median_AML = median(x_aml, na.rm = TRUE),
        median_Normal = median(x_norm, na.rm = TRUE),
        delta_mean = mean_AML - mean_Normal,
        delta_median = median_AML - median_Normal,
        p_value = ifelse(is.null(wt), NA_real_, wt$p.value)
      )
    }) %>%
      mutate(
        FDR = p.adjust(p_value, method = "BH"),
        direction = case_when(
          FDR < 0.05 & delta_mean > 0 ~ "Higher in AML",
          FDR < 0.05 & delta_mean < 0 ~ "Higher in normal",
          TRUE ~ "NS"
        )
      ) %>%
      arrange(FDR, desc(abs(delta_mean)))
    
    write.csv(
      cite_de_cd14,
      file.path(out_dir, "CD14_CITE_DSB_AML_vs_Normal_numbat.csv"),
      row.names = FALSE
    )
  }
}

# ----------------------------
# Save object
# ----------------------------

saveRDS(
  seu_cd14,
  file.path(out_dir, "seu_cd14_LK1_LK2_numbat_RNA_CITE_WNN.rds")
)

cat("\nDone.\n")
cat("Output directory:\n", out_dir, "\n")
cat("Plot directory:\n", plot_dir, "\n")