# ------------------------------------------------------------------
# LK1 re-sequenced - step 04 of 5
#
# Pseudobulk DE at each depth, then cluster HBDN501-AML-KMT2A and score separation.
#
# Supersedes earlier revisions of the same analysis, which are not kept here:
#   scripts/1d.DE_subsampling.R  (2026-05-28, first pass)
#   scripts/1d3.DE_subsampling_LK2.R  (byte-identical to this despite the name)
#
# Frozen for the lab archive 2026-08-03 from scripts/1d2.DE_subsampling.R (mtime 2026-06-01).
# md5 of the original: 0ec5c5b879aa541b853d0ae561f41c00
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# ============================================================
# Cell Ranger downsampled matrices:
# 1) load real downsampled count matrices
# 2) map cells to annotated Seurat metadata
# 3) run pseudobulk DE at each depth
# 4) cluster HBDN501-AML-KMT2A at each depth and assess separation
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(edgeR)
  library(patchwork)
})

set.seed(123)

# Optional metrics
has_mclust <- requireNamespace("mclust", quietly = TRUE)
has_cluster <- requireNamespace("cluster", quietly = TRUE)

# ============================================================
# Paths
# ============================================================

downsample_base <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/results/cellranger_downsampled_reseq/260522_VH01624_461_222JLJVNX"

run_annotation <- "260423_VH01624_453_222HWMYNX"
annotation_dir <- file.path("../results/seurat_annotated", run_annotation)

seurat_file <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

out_dir <- file.path(
  "../results/data_integration",
  run_annotation,
  "cellranger_downsampled_DE_and_clustering"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

de_dir <- file.path(out_dir, "DE_results")
umap_dir <- file.path(out_dir, "HBDN501_downsampled_umaps")
qc_dir <- file.path(out_dir, "QC")

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Parameters
# ============================================================

assay_use <- "RNA"
sample_col <- "sampleID"
celltype_col <- "predicted_CellType_Broad"

fdr_cutoff <- 0.05
logfc_cutoff <- 1
bcv <- 0.4

cluster_sample <- "HBDN501-AML-KMT2A"
cluster_label_col <- "predicted_CellType_Broad"

comparisons <- tibble::tribble(
  ~sampleID,          ~comparison_name,        ~group_a,        ~group_b,
  "HBDN392-AML-MDS",  "pDC_vs_cDC",            "cDC",          "pDC",
  "HBDN392-AML-MDS",  "NaiveT_vs_CD8MemoryT",  "Naive T",      "CD8 Memory T",
  "HBDN206-MNpCT",    "HSCMPP_vs_LateGMP",     "HSC MPP",      "Late GMP",
  "HBDN206-MNpCT",    "NaiveT_vs_CD8MemoryT",  "Naive T",      "CD8 Memory T"
)

# ============================================================
# Helpers
# ============================================================

clean_barcode <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace(".*_", "") %>%
    stringr::str_replace("-1$", "")
}

find_cr_matrix_dir <- function(run_dir) {
  hits <- list.dirs(run_dir, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == "filtered_feature_bc_matrix"]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

find_cr_metrics_file <- function(run_dir) {
  hits <- list.files(
    run_dir,
    pattern = "^metrics_summary\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

read_cr_matrix <- function(matrix_dir) {
  mat <- Seurat::Read10X(data.dir = matrix_dir)
  if (is.list(mat)) {
    mat <- mat[["Gene Expression"]]
  }
  mat
}

make_pseudobulk_from_cells <- function(count_mat, cells) {
  cells <- intersect(cells, colnames(count_mat))
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}

run_edgeR_exact <- function(count_mat, group_vec, bcv = 0.4) {
  
  group_vec <- factor(group_vec, levels = c("group_a", "group_b"))
  
  dge <- edgeR::DGEList(counts = count_mat, group = group_vec)
  
  keep <- edgeR::filterByExpr(dge, group = group_vec)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) < 10) return(NULL)
  
  dge <- edgeR::calcNormFactors(dge)
  
  et <- edgeR::exactTest(dge, dispersion = bcv^2)
  
  edgeR::topTags(et, n = Inf)$table %>%
    tibble::rownames_to_column("gene") %>%
    tibble::as_tibble()
}

read_cr_qc <- function(mat, dataset, fraction) {
  tibble(
    dataset = dataset,
    fraction = fraction,
    barcode = colnames(mat),
    barcode_clean = clean_barcode(colnames(mat)),
    nUMI = as.numeric(Matrix::colSums(mat)),
    nGene = as.integer(Matrix::colSums(mat > 0))
  )
}

compute_cluster_metrics <- function(seu_obj, label_col, cluster_col = "seurat_clusters") {
  
  meta_i <- seu_obj@meta.data
  
  out <- tibble(
    n_cells = nrow(meta_i),
    n_clusters = dplyr::n_distinct(meta_i[[cluster_col]]),
    n_labels = dplyr::n_distinct(meta_i[[label_col]]),
    ari_cluster_vs_label = NA_real_,
    mean_silhouette_label_pca = NA_real_
  )
  
  if (has_mclust) {
    out$ari_cluster_vs_label <- mclust::adjustedRandIndex(
      as.character(meta_i[[cluster_col]]),
      as.character(meta_i[[label_col]])
    )
  }
  
  if (has_cluster && "pca" %in% names(seu_obj@reductions)) {
    emb <- Seurat::Embeddings(seu_obj, "pca")[, 1:min(20, ncol(Seurat::Embeddings(seu_obj, "pca"))), drop = FALSE]
    labels <- as.factor(meta_i[[label_col]])
    
    if (length(unique(labels)) > 1 && nrow(emb) > length(unique(labels))) {
      d <- dist(emb)
      sil <- cluster::silhouette(as.integer(labels), d)
      out$mean_silhouette_label_pca <- mean(sil[, "sil_width"], na.rm = TRUE)
    }
  }
  
  out
}

# ============================================================
# Discover completed downsample Cell Ranger outputs
# ============================================================

outs <- tibble(
  run_dir = list.dirs(downsample_base, recursive = FALSE, full.names = TRUE),
  run_name = basename(run_dir)
) %>%
  filter(stringr::str_detect(run_name, "LK1-GEX_frac_")) %>%
  mutate(
    fraction_label = stringr::str_extract(run_name, "frac_[0-9]+"),
    fraction_num = as.numeric(stringr::str_remove(fraction_label, "frac_")),
    fraction = fraction_num / 100,
    matrix_dir = purrr::map_chr(run_dir, find_cr_matrix_dir),
    metrics_file = purrr::map_chr(run_dir, find_cr_metrics_file),
    has_matrix = !is.na(matrix_dir),
    has_metrics = !is.na(metrics_file)
  ) %>%
  filter(has_matrix) %>%
  arrange(fraction)

if (nrow(outs) == 0) {
  stop("No completed Cell Ranger filtered_feature_bc_matrix folders found.")
}

print(outs, n = Inf, width = Inf)

readr::write_csv(
  outs,
  file.path(qc_dir, "cellranger_downsampled_dirs_used.csv")
)

# ============================================================
# Load annotated Seurat object metadata
# ============================================================

message("Loading annotated Seurat object: ", seurat_file)

seu_annot <- readRDS(seurat_file)

annot_meta <- seu_annot@meta.data %>%
  tibble::rownames_to_column("annot_cell") %>%
  mutate(
    barcode_clean = clean_barcode(annot_cell)
  )

stopifnot(sample_col %in% colnames(annot_meta))
stopifnot(celltype_col %in% colnames(annot_meta))

readr::write_csv(
  annot_meta,
  file.path(qc_dir, "annotated_seurat_metadata_used.csv")
)

# ============================================================
# Load all downsample matrices and attach annotation
# ============================================================

downsample_objects <- purrr::pmap(
  outs %>% select(run_name, fraction, matrix_dir),
  function(run_name, fraction, matrix_dir) {
    
    message("Reading downsample: ", run_name)
    
    mat <- read_cr_matrix(matrix_dir)
    
    barcode_map <- tibble(
      barcode = colnames(mat),
      barcode_clean = clean_barcode(colnames(mat))
    ) %>%
      left_join(annot_meta, by = "barcode_clean")
    
    matched <- !is.na(barcode_map[[sample_col]])
    
    message("  Cells in matrix: ", ncol(mat))
    message("  Matched to annotation: ", sum(matched))
    
    mat_matched <- mat[, barcode_map$barcode[matched], drop = FALSE]
    
    meta_matched <- barcode_map %>%
      filter(matched) %>%
      as.data.frame()
    
    rownames(meta_matched) <- meta_matched$barcode
    
    list(
      run_name = run_name,
      fraction = fraction,
      matrix = mat_matched,
      meta = meta_matched,
      qc = read_cr_qc(mat, run_name, fraction) %>%
        left_join(
          barcode_map %>%
            select(barcode, barcode_clean, all_of(sample_col), all_of(celltype_col)),
          by = c("barcode", "barcode_clean")
        )
    )
  }
)

names(downsample_objects) <- outs$run_name

cell_qc <- purrr::map_dfr(downsample_objects, "qc")

readr::write_csv(
  cell_qc,
  file.path(qc_dir, "downsampled_cell_qc_with_annotation.csv")
)

cell_counts_by_depth <- cell_qc %>%
  dplyr::filter(!is.na(.data[[sample_col]]), !is.na(.data[[celltype_col]])) %>%
  dplyr::count(dataset, fraction, .data[[sample_col]], .data[[celltype_col]], name = "n_cells") %>%
  arrange(fraction, .data[[sample_col]], desc(n_cells))

readr::write_csv(
  cell_counts_by_depth,
  file.path(qc_dir, "cell_counts_by_depth_sample_celltype.csv")
)

# ============================================================
# Pseudobulk DE at each Cell Ranger downsample depth
# ============================================================

run_de_for_depth <- function(ds_obj, comparison_row) {
  
  sample_i <- comparison_row$sampleID
  comparison_name <- comparison_row$comparison_name
  group_a_i <- comparison_row$group_a
  group_b_i <- comparison_row$group_b
  
  mat <- ds_obj$matrix
  meta <- ds_obj$meta
  fraction <- ds_obj$fraction
  run_name <- ds_obj$run_name
  
  message("DE: ", run_name, " | ", sample_i, " | ", comparison_name)
  
  meta_i <- meta %>%
    filter(
      .data[[sample_col]] == sample_i,
      .data[[celltype_col]] %in% c(group_a_i, group_b_i)
    )
  
  cells_a <- meta_i %>%
    filter(.data[[celltype_col]] == group_a_i) %>%
    pull(barcode)
  
  cells_b <- meta_i %>%
    filter(.data[[celltype_col]] == group_b_i) %>%
    pull(barcode)
  
  if (length(cells_a) == 0 || length(cells_b) == 0) {
    warning("Missing cells: ", run_name, " | ", comparison_name)
    return(tibble(
      dataset = run_name,
      fraction = fraction,
      sampleID = sample_i,
      comparison = comparison_name,
      group_a = group_a_i,
      group_b = group_b_i,
      n_cells_group_a = length(cells_a),
      n_cells_group_b = length(cells_b),
      n_up = NA_integer_,
      n_down = NA_integer_,
      n_de_total = NA_integer_,
      status = "missing_cells"
    ))
  }
  
  pb_a <- make_pseudobulk_from_cells(mat, cells_a)
  pb_b <- make_pseudobulk_from_cells(mat, cells_b)
  
  names(pb_a) <- rownames(mat)
  names(pb_b) <- rownames(mat)
  
  count_mat <- cbind(
    group_a = pb_a,
    group_b = pb_b
  )
  
  rownames(count_mat) <- rownames(mat)
  
  de_res <- run_edgeR_exact(
    count_mat = count_mat,
    group_vec = c("group_a", "group_b"),
    bcv = bcv
  )
  
  if (is.null(de_res)) {
    return(tibble(
      dataset = run_name,
      fraction = fraction,
      sampleID = sample_i,
      comparison = comparison_name,
      group_a = group_a_i,
      group_b = group_b_i,
      n_cells_group_a = length(cells_a),
      n_cells_group_b = length(cells_b),
      n_up = NA_integer_,
      n_down = NA_integer_,
      n_de_total = NA_integer_,
      status = "edgeR_failed"
    ))
  }
  
  de_res <- de_res %>%
    mutate(
      dataset = run_name,
      fraction = fraction,
      sampleID = sample_i,
      comparison = comparison_name,
      group_a = group_a_i,
      group_b = group_b_i,
      n_cells_group_a = length(cells_a),
      n_cells_group_b = length(cells_b),
      lib_size_group_a = sum(pb_a),
      lib_size_group_b = sum(pb_b),
      direction = case_when(
        FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Up in group_b",
        FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Down in group_b",
        TRUE ~ "Not significant"
      )
    )
  
  safe_name <- paste(
    run_name,
    sample_i,
    comparison_name,
    sep = "__"
  ) %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_")
  
  readr::write_csv(
    de_res,
    file.path(de_dir, paste0(safe_name, "_DE_results.csv"))
  )
  
  tibble(
    dataset = run_name,
    fraction = fraction,
    sampleID = sample_i,
    comparison = comparison_name,
    group_a = group_a_i,
    group_b = group_b_i,
    n_cells_group_a = length(cells_a),
    n_cells_group_b = length(cells_b),
    lib_size_group_a = sum(pb_a),
    lib_size_group_b = sum(pb_b),
    n_up = sum(de_res$direction == "Up in group_b", na.rm = TRUE),
    n_down = sum(de_res$direction == "Down in group_b", na.rm = TRUE),
    n_de_total = n_up + n_down,
    status = "ok"
  )
}

de_summary <- purrr::map_dfr(downsample_objects, function(ds_obj) {
  purrr::pmap_dfr(
    comparisons,
    function(sampleID, comparison_name, group_a, group_b) {
      run_de_for_depth(
        ds_obj,
        tibble(
          sampleID = sampleID,
          comparison_name = comparison_name,
          group_a = group_a,
          group_b = group_b
        )
      )
    }
  )
})

readr::write_csv(
  de_summary,
  file.path(de_dir, "cellranger_downsampled_pseudobulk_DE_summary.csv")
)

saveRDS(
  de_summary,
  file.path(de_dir, "cellranger_downsampled_pseudobulk_DE_summary.rds")
)

# ============================================================
# DE plots
# ============================================================

plot_df <- de_summary %>%
  filter(status == "ok") %>%
  pivot_longer(
    cols = c(n_up, n_down),
    names_to = "direction",
    values_to = "n_genes"
  ) %>%
  mutate(
    direction = recode(
      direction,
      n_up = "Up in group_b",
      n_down = "Down in group_b"
    ),
    comparison_label = paste0(
      sampleID,
      "\n",
      group_b,
      " vs ",
      group_a
    )
  )

p_de <- ggplot(
  plot_df,
  aes(
    x = fraction * 100,
    y = n_genes,
    colour = direction,
    group = direction
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ comparison_label, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Pseudobulk DE across Cell Ranger downsampled depths",
    subtitle = paste0("edgeR exactTest; BCV = ", bcv, "; FDR < ", fdr_cutoff, ", |log2FC| > ", logfc_cutoff),
    x = "Cell Ranger downsample fraction (%)",
    y = "Number of DE genes",
    colour = "Direction"
  )

ggsave(
  file.path(de_dir, "cellranger_downsampled_DE_up_down.pdf"),
  p_de,
  width = 12,
  height = 8
)

p_total <- de_summary %>%
  filter(status == "ok") %>%
  mutate(
    comparison_label = paste0(
      sampleID,
      "\n",
      group_b,
      " vs ",
      group_a
    )
  ) %>%
  ggplot(aes(x = fraction * 100, y = n_de_total)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ comparison_label, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Total DE genes across Cell Ranger downsampled depths",
    x = "Cell Ranger downsample fraction (%)",
    y = "Total DE genes"
  )

ggsave(
  file.path(de_dir, "cellranger_downsampled_DE_total.pdf"),
  p_total,
  width = 12,
  height = 8
)

# ============================================================
# Clustering / UMAP QC for HBDN501-AML-KMT2A at each depth
# ============================================================

cluster_metrics <- purrr::map_dfr(downsample_objects, function(ds_obj) {
  
  run_name <- ds_obj$run_name
  fraction <- ds_obj$fraction
  mat <- ds_obj$matrix
  meta <- ds_obj$meta
  
  message("Clustering sample ", cluster_sample, " at depth: ", run_name)
  
  meta_i <- meta %>%
    filter(
      .data[[sample_col]] == cluster_sample,
      !is.na(.data[[cluster_label_col]])
    )
  
  cells_i <- intersect(meta_i$barcode, colnames(mat))
  
  if (length(cells_i) < 50) {
    warning("Too few cells for clustering: ", run_name, " | n=", length(cells_i))
    return(tibble(
      dataset = run_name,
      fraction = fraction,
      n_cells = length(cells_i),
      n_clusters = NA_integer_,
      n_labels = NA_integer_,
      ari_cluster_vs_label = NA_real_,
      mean_silhouette_label_pca = NA_real_,
      status = "too_few_cells"
    ))
  }
  
  seu_i <- CreateSeuratObject(
    counts = mat[, cells_i, drop = FALSE],
    meta.data = meta_i %>%
      filter(barcode %in% cells_i) %>%
      as.data.frame()
  )
  
  seu_i <- NormalizeData(seu_i, verbose = FALSE)
  seu_i <- FindVariableFeatures(seu_i, nfeatures = 2000, verbose = FALSE)
  seu_i <- ScaleData(seu_i, verbose = FALSE)
  seu_i <- RunPCA(seu_i, npcs = 30, verbose = FALSE)
  seu_i <- FindNeighbors(seu_i, dims = 1:20, verbose = FALSE)
  seu_i <- FindClusters(seu_i, resolution = 0.5, verbose = FALSE)
  seu_i <- RunUMAP(seu_i, dims = 1:20, verbose = FALSE)
  
  metrics_i <- compute_cluster_metrics(
    seu_i,
    label_col = cluster_label_col,
    cluster_col = "seurat_clusters"
  ) %>%
    mutate(
      dataset = run_name,
      fraction = fraction,
      status = "ok",
      .before = 1
    )
  
  p_label <- DimPlot(
    seu_i,
    reduction = "umap",
    group.by = cluster_label_col,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0(cluster_sample, " | ", run_name, " | broad labels"))
  
  p_cluster <- DimPlot(
    seu_i,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0(cluster_sample, " | ", run_name, " | Seurat clusters"))
  
  ggsave(
    file.path(
      umap_dir,
      paste0(run_name, "_", cluster_sample, "_UMAP_broad_celltypes.pdf")
    ),
    p_label,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(
      umap_dir,
      paste0(run_name, "_", cluster_sample, "_UMAP_clusters.pdf")
    ),
    p_cluster,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(
      umap_dir,
      paste0(run_name, "_", cluster_sample, "_UMAP_combined.pdf")
    ),
    p_label + p_cluster,
    width = 14,
    height = 6
  )
  
  saveRDS(
    seu_i,
    file.path(
      umap_dir,
      paste0(run_name, "_", cluster_sample, "_clustered_seurat.rds")
    )
  )
  
  metrics_i
})

readr::write_csv(
  cluster_metrics,
  file.path(umap_dir, "HBDN501_cluster_separation_metrics.csv")
)

p_cluster_metrics <- cluster_metrics %>%
  filter(status == "ok") %>%
  pivot_longer(
    cols = c(ari_cluster_vs_label, mean_silhouette_label_pca),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = fraction * 100, y = value)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = paste0(cluster_sample, ": cluster/label separation across downsample depths"),
    x = "Cell Ranger downsample fraction (%)",
    y = "Metric value"
  )

ggsave(
  file.path(umap_dir, "HBDN501_cluster_separation_metrics.pdf"),
  p_cluster_metrics,
  width = 8,
  height = 5
)

message("Done. Outputs written to: ", out_dir)

# ============================================================
# Cell annotation robustness across downsampled depths
# ============================================================

suppressPackageStartupMessages({
  library(BoneMarrowMap)
  library(symphony)
})

annotation_eval_dir <- file.path(out_dir, "annotation_robustness")
dir.create(annotation_eval_dir, recursive = TRUE, showWarnings = FALSE)

projection_path <- "../annotation/"

ref <- readRDS(file.path(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
ref$save_uwot_path <- file.path(projection_path, "BoneMarrowMap_uwot_model.uwot")

annotation_truth_cols <- c(
  "predicted_CellType_Broad",
  "predicted_CellType"
)

annotation_truth_cols <- annotation_truth_cols[
  annotation_truth_cols %in% colnames(annot_meta)
]

stopifnot(length(annotation_truth_cols) > 0)

run_annotation_for_depth <- function(ds_obj) {
  
  run_name <- ds_obj$run_name
  fraction <- ds_obj$fraction
  mat <- ds_obj$matrix
  meta <- ds_obj$meta
  
  message("Annotating downsampled object: ", run_name)
  
  if (ncol(mat) < 50) {
    warning("Too few cells for annotation: ", run_name)
    return(NULL)
  }
  
  seu_i <- CreateSeuratObject(
    counts = mat,
    meta.data = meta
  )
  
  seu_i$sampleID <- seu_i[[sample_col]][, 1]
  
  seu_i <- NormalizeData(seu_i, verbose = FALSE)
  seu_i <- FindVariableFeatures(seu_i, nfeatures = 2000, verbose = FALSE)
  
  query_i <- map_Query(
    query = seu_i,
    ref_obj = ref,
    vars = "sampleID"
  )
  
  query_i <- calculate_MappingError(
    query_i,
    reference = ref,
    MAD_threshold = 2.5
  )
  
  query_i <- predict_CellTypes(
    query_obj = query_i,
    ref_obj = ref,
    final_label = "predicted_CellType_downsampled"
  )
  
  # Save Symphony mapping QC plot
  pdf(
    file.path(annotation_eval_dir, paste0(run_name, "_mapping_error_QC.pdf")),
    width = 8,
    height = 6
  )
  print(plot_MappingErrorQC(query_i))
  dev.off()
  
  # Save projected annotation plot
  pdf(
    file.path(annotation_eval_dir, paste0(run_name, "_projected_predicted_celltypes.pdf")),
    width = 14,
    height = 10
  )
  print(
    DimPlot(
      subset(query_i, mapping_error_QC == "Pass"),
      group.by = "predicted_CellType_downsampled",
      label = TRUE,
      repel = TRUE,
      raster = FALSE
    ) +
      ggtitle(paste0(run_name, " | Symphony predicted cell types"))
  )
  dev.off()
  
  pred_meta <- query_i@meta.data %>%
    tibble::as_tibble(rownames = "cell_barcode") %>%
    dplyr::mutate(
      barcode_clean = clean_barcode(cell_barcode),
      dataset = run_name,
      fraction = fraction
    )
  
  saveRDS(
    query_i,
    file.path(annotation_eval_dir, paste0(run_name, "_symphony_annotated.rds"))
  )
  
  pred_meta
}

annotation_meta <- purrr::map_dfr(
  downsample_objects,
  run_annotation_for_depth
)

readr::write_csv(
  annotation_meta,
  file.path(annotation_eval_dir, "downsampled_symphony_annotation_metadata.csv")
)

# ============================================================
# Compare downsampled annotation to 100% annotation
# ============================================================

annotation_eval <- annotation_meta %>%
  left_join(
    annot_meta %>%
      select(
        barcode_clean,
        sampleID_truth = all_of(sample_col),
        true_CellType_Broad = predicted_CellType_Broad,
        true_CellType = predicted_CellType
      ),
    by = "barcode_clean"
  ) %>%
  mutate(
    pass_mapping = mapping_error_QC == "Pass",
    
    broad_match = !is.na(predicted_CellType_downsampled_Broad) &
      !is.na(true_CellType_Broad) &
      predicted_CellType_downsampled_Broad == true_CellType_Broad,
    
    fine_match = !is.na(predicted_CellType_downsampled) &
      !is.na(true_CellType) &
      predicted_CellType_downsampled == true_CellType
  )

readr::write_csv(
  annotation_eval,
  file.path(annotation_eval_dir, "downsampled_annotation_vs_full_depth_per_cell.csv")
)

annotation_summary <- annotation_eval %>%
  group_by(dataset, fraction) %>%
  summarise(
    n_cells = n(),
    n_mapping_pass = sum(pass_mapping, na.rm = TRUE),
    mapping_pass_rate = mean(pass_mapping, na.rm = TRUE),
    
    broad_accuracy_all_cells = mean(broad_match, na.rm = TRUE),
    broad_accuracy_pass_only = mean(broad_match[pass_mapping], na.rm = TRUE),
    
    fine_accuracy_all_cells = mean(fine_match, na.rm = TRUE),
    fine_accuracy_pass_only = mean(fine_match[pass_mapping], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(fraction)

readr::write_csv(
  annotation_summary,
  file.path(annotation_eval_dir, "downsampled_annotation_accuracy_summary.csv")
)

# ============================================================
# Confusion matrices
# ============================================================

confusion_broad <- annotation_eval %>%
  filter(
    pass_mapping,
    !is.na(true_CellType_Broad),
    !is.na(predicted_CellType_downsampled_Broad)
  ) %>%
  count(
    dataset,
    fraction,
    true_CellType_Broad,
    predicted_CellType_downsampled_Broad,
    name = "n_cells"
  ) %>%
  group_by(dataset, fraction, true_CellType_Broad) %>%
  mutate(
    percent_of_truth = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

readr::write_csv(
  confusion_broad,
  file.path(annotation_eval_dir, "downsampled_annotation_confusion_broad.csv")
)

confusion_fine <- annotation_eval %>%
  filter(
    pass_mapping,
    !is.na(true_CellType),
    !is.na(predicted_CellType_downsampled)
  ) %>%
  count(
    dataset,
    fraction,
    true_CellType,
    predicted_CellType_downsampled,
    name = "n_cells"
  ) %>%
  group_by(dataset, fraction, true_CellType) %>%
  mutate(
    percent_of_truth = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

readr::write_csv(
  confusion_fine,
  file.path(annotation_eval_dir, "downsampled_annotation_confusion_fine.csv")
)

# ============================================================
# Annotation robustness plots
# ============================================================

p_annotation_accuracy <- annotation_summary %>%
  select(
    dataset,
    fraction,
    mapping_pass_rate,
    broad_accuracy_pass_only,
    fine_accuracy_pass_only
  ) %>%
  pivot_longer(
    cols = c(
      mapping_pass_rate,
      broad_accuracy_pass_only,
      fine_accuracy_pass_only
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = fraction * 100, y = value)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Cell annotation robustness across downsampled sequencing depths",
    x = "Cell Ranger downsample fraction (%)",
    y = "Proportion"
  )

ggsave(
  file.path(annotation_eval_dir, "annotation_accuracy_across_downsample_depths.pdf"),
  p_annotation_accuracy,
  width = 9,
  height = 3.5
)
library(gtools)
confusion_fine_plot <- confusion_fine %>%
  dplyr::mutate(
    fraction = factor(
      fraction,
      levels = mixedsort(unique(as.character(fraction)))
    ),
    true_CellType = factor(
      true_CellType,
      levels = mixedsort(unique(as.character(true_CellType)))
    ),
    predicted_CellType_downsampled = factor(
      predicted_CellType_downsampled,
      levels = levels(true_CellType)
    )
  )

p_confusion_fine <- confusion_fine_plot %>%
  ggplot(
    aes(
      x = predicted_CellType_downsampled,
      y = true_CellType,
      fill = percent_of_truth
    )
  ) +
  geom_tile() +
  facet_wrap(~ fraction) +
  theme_bw(base_size = 9) +
  labs(
    title = "Fine cell-type annotation agreement with full-depth annotation",
    x = "Downsampled Symphony fine annotation",
    y = "Full-depth fine annotation",
    fill = "% of true label"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(annotation_eval_dir, "annotation_confusion_fine_across_depths.pdf"),
  p_confusion_fine,
  width = 25,
  height = 14
)


# ----------------------------
# Broad confusion plot
# ----------------------------

confusion_broad_plot <- confusion_broad %>%
  dplyr::mutate(
    fraction = factor(
      fraction,
      levels = mixedsort(unique(as.character(fraction)))
    ),
    true_CellType_Broad = factor(
      true_CellType_Broad,
      levels = mixedsort(unique(as.character(true_CellType_Broad)))
    ),
    predicted_CellType_downsampled_Broad = factor(
      predicted_CellType_downsampled_Broad,
      levels = levels(true_CellType_Broad)
    )
  )

p_confusion_broad <- confusion_broad_plot %>%
  ggplot(
    aes(
      x = predicted_CellType_downsampled_Broad,
      y = true_CellType_Broad,
      fill = percent_of_truth
    )
  ) +
  geom_tile() +
  facet_wrap(~ fraction) +
  theme_bw(base_size = 9) +
  labs(
    title = "Broad cell-type annotation agreement with full-depth annotation",
    x = "Downsampled Symphony broad annotation",
    y = "Full-depth broad annotation",
    fill = "% of true label"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(annotation_eval_dir, "annotation_confusion_broad_across_depths.pdf"),
  p_confusion_broad,
  width = 16,
  height = 10
)

message("Annotation robustness outputs written to: ", annotation_eval_dir)