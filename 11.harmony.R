library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scIntegrationMetrics)

# remotes::install_github("carmonalab/scIntegrationMetrics")

run <- "260528_VH01624_464_222K7VKNX"
proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

annotation_dir <- file.path(proj, "results/seurat_annotated", run)
seurat_file <- file.path(annotation_dir, "demux_singlets_annotated_seurat.rds")

out_dir <- file.path(annotation_dir, "harmony_parameter_qc_scIntegrationMetrics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(seurat_file)

# ----------------------------
# Metadata columns
# ----------------------------
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

meta_cols <- colnames(seu@meta.data)
print(meta_cols)

annotation_col <- "predicted_CellType_Broad"

batch_col <- "sample_name"

message("Using annotation column: ", annotation_col)
message("Using batch column: ", batch_col)

# Remove cells without batch / annotation
keep <- !is.na(seu@meta.data[[batch_col]]) &
  !is.na(seu@meta.data[[annotation_col]])

seu <- subset(seu, cells = colnames(seu)[keep])

# ----------------------------
# Standard Seurat processing
# ----------------------------

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# ----------------------------
# Baseline non-Harmony UMAP
# ----------------------------

seu <- RunUMAP(
  seu,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap_pca",
  reduction.key = "pcaUMAP_"
)

# ----------------------------
# Helper function
# ----------------------------

score_reduction <- function(
    seu,
    reduction,
    method_name,
    batch_col,
    annotation_col,
    ndim = 30,
    iLISI_perplexity = 30
) {
  
  metrics <- scIntegrationMetrics::getIntegrationMetrics(
    object = seu,
    meta.label = annotation_col,
    meta.batch = batch_col,
    method.reduction = reduction,
    ndim = ndim,
    iLISI_perplexity = iLISI_perplexity
  )
  
  data.frame(
    method = method_name,
    reduction = reduction,
    metric = names(unlist(metrics)),
    value = as.numeric(unlist(metrics))
  )
}

plot_umap_set <- function(seu, reduction, method_name) {
  
  p_batch <- DimPlot(
    seu,
    reduction = reduction,
    group.by = batch_col
  ) +
    ggtitle(paste(method_name, "by", batch_col))
  
  p_annot <- DimPlot(
    seu,
    reduction = reduction,
    group.by = annotation_col,
    label = TRUE,
    repel = TRUE,
    cols = celltype_cols
  ) +
    ggtitle(paste(method_name, "by", annotation_col))
  
  ggsave(
    file.path(out_dir, paste0(method_name, "_batch.pdf")),
    p_batch,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(out_dir, paste0(method_name, "_annotation.pdf")),
    p_annot,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(out_dir, paste0(method_name, "_combined.pdf")),
    p_batch + p_annot,
    width = 16,
    height = 6
  )
}

# ----------------------------
# Score PCA baseline
# ----------------------------

all_scores <- list()

all_scores[["pca_no_harmony"]] <- score_reduction(
  seu = seu,
  reduction = "pca",
  method_name = "pca_no_harmony",
  batch_col = batch_col,
  annotation_col = annotation_col,
  ndim = 30
)

plot_umap_set(
  seu,
  reduction = "umap_pca",
  method_name = "pca_no_harmony_umap"
)

# ----------------------------
# Harmony parameter sweep
# ----------------------------

param_grid <- expand.grid(
  theta = c(0, 1, 2, 5, 10),
  lambda = c(0.1, 0.5, 1, 2,5),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(param_grid))) {
  
  theta_i <- param_grid$theta[i]
  lambda_i <- param_grid$lambda[i]
  
  harmony_reduction <- paste0(
    "harmony_theta", theta_i,
    "_lambda", lambda_i
  )
  
  umap_reduction <- paste0("umap_", harmony_reduction)
  
  message("Running ", harmony_reduction)
  
  seu <- RunHarmony(
    object = seu,
    group.by.vars = batch_col,
    reduction = "pca",
    dims.use = 1:50,
    assay.use = "RNA",
    theta = theta_i,
    lambda = lambda_i,
    reduction.save = harmony_reduction,
    verbose = FALSE
  )
  
  seu <- RunUMAP(
    seu,
    reduction = harmony_reduction,
    dims = 1:30,
    reduction.name = umap_reduction,
    reduction.key = paste0("hUMAP", i, "_"),
    verbose = FALSE
  )
  
  all_scores[[harmony_reduction]] <- score_reduction(
    seu = seu,
    reduction = harmony_reduction,
    method_name = harmony_reduction,
    batch_col = batch_col,
    annotation_col = annotation_col
  )
  
  plot_umap_set(
    seu,
    reduction = umap_reduction,
    method_name = paste0(harmony_reduction, "_umap")
  )
}

# ----------------------------
# Save metrics
# ----------------------------

scores_long <- bind_rows(all_scores)

write.csv(
  scores_long,
  file.path(out_dir, "integration_metrics_long.csv"),
  row.names = FALSE
)

scores_wide <- scores_long |>
  tidyr::pivot_wider(
    names_from = metric,
    values_from = value
  )

write.csv(
  scores_wide,
  file.path(out_dir, "integration_metrics_wide.csv"),
  row.names = FALSE
)

# ----------------------------
# Ranking
# ----------------------------

scores_ranked <- scores_wide |>
  mutate(
    rank_norm_iLISI = rank(-norm_iLISI, ties.method = "average"),
    rank_norm_cLISI = rank(-norm_cLISI, ties.method = "average"),
    rank_celltype_ASW = rank(-celltype_ASW, ties.method = "average"),
    
    overall_rank =
      rank_norm_iLISI +
      rank_norm_cLISI +
      rank_celltype_ASW
  ) |>
  arrange(overall_rank)

write.csv(
  scores_ranked,
  file.path(out_dir, "integration_metrics_ranked.csv"),
  row.names = FALSE
)

# ----------------------------
# Plot metrics
# ----------------------------

p_metrics <- scores_long |>
  filter(metric %in% c(
    "norm_iLISI",
    "norm_cLISI",
    "celltype_ASW_means"
  )) |>
  ggplot(aes(x = reorder(method, value), y = value)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ metric, scales = "free_x") +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = NULL)

ggsave(
  file.path(out_dir, "integration_metrics_overview.pdf"),
  p_metrics,
  width = 14,
  height = 5
)

# ----------------------------
# Save object
# ----------------------------

saveRDS(
  seu,
  file.path(out_dir, "demux_singlets_annotated_seurat_harmony_parameter_sweep.rds")
)

message("Done. Outputs written to: ", out_dir)


metrics_to_plot <- c(
  "norm_iLISI",
  "norm_cLISI",
  "celltype_ASW_means"
)

scores_plot <- scores_long |>
  filter(metric %in% metrics_to_plot) |>
  mutate(
    is_pca = method == "pca_no_harmony",
    theta = as.numeric(stringr::str_match(method, "theta([0-9.]+)")[, 2]),
    lambda = as.numeric(stringr::str_match(method, "lambda([0-9.]+)")[, 2])
  )

pca_baseline <- scores_plot |>
  filter(is_pca) |>
  select(metric, pca_value = value)

harmony_plot <- scores_plot |>
  filter(!is_pca) |>
  left_join(pca_baseline, by = "metric")

# ----------------------------
# 1. Metric vs theta, coloured by lambda
# ----------------------------

p_theta <- harmony_plot |>
  ggplot(aes(
    x = theta,
    y = value,
    group = factor(lambda),
    color = factor(lambda)
  )) +
  geom_hline(
    aes(yintercept = pca_value),
    linetype = 2,
    color = "black"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Harmony parameter sweep: metric vs theta",
    subtitle = "Dashed black line = PCA without Harmony",
    x = "Harmony theta",
    y = "Metric value",
    color = "lambda"
  )

ggsave(
  file.path(out_dir, "integration_metrics_vs_theta.pdf"),
  p_theta,
  width = 13,
  height = 4
)

# ----------------------------
# 2. Metric vs lambda, coloured by theta
# ----------------------------

p_lambda <- harmony_plot |>
  ggplot(aes(
    x = lambda,
    y = value,
    group = factor(theta),
    color = factor(theta)
  )) +
  geom_hline(
    aes(yintercept = pca_value),
    linetype = 2,
    color = "black"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Harmony parameter sweep: metric vs lambda",
    subtitle = "Dashed black line = PCA without Harmony",
    x = "Harmony lambda",
    y = "Metric value",
    color = "theta"
  )

ggsave(
  file.path(out_dir, "integration_metrics_vs_lambda.pdf"),
  p_lambda,
  width = 13,
  height = 4
)

# ----------------------------
# 3. Heatmap-style plot: theta x lambda per metric
# ----------------------------

p_heat <- harmony_plot |>
  ggplot(aes(
    x = factor(theta),
    y = factor(lambda),
    fill = value
  )) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), size = 3) +
  facet_wrap(~ metric, scales = "free") +
  theme_bw(base_size = 11) +
  labs(
    title = "Harmony parameter sweep heatmap",
    x = "theta",
    y = "lambda",
    fill = "value"
  )

ggsave(
  file.path(out_dir, "integration_metrics_theta_lambda_heatmap.pdf"),
  p_heat,
  width = 13,
  height = 8
)

# ----------------------------
# 4. Ranked objective score
#    Higher score = better.
#    Adjust weights if you care more about biology than batch mixing.
# ----------------------------

scores_objective <- scores_wide |>
  mutate(
    is_pca = method == "pca_no_harmony",
    theta = as.numeric(stringr::str_match(method, "theta([0-9.]+)")[, 2]),
    lambda = as.numeric(stringr::str_match(method, "lambda([0-9.]+)")[, 2])
  ) |>
  mutate(
    z_norm_iLISI = as.numeric(scale(norm_iLISI)),
    z_norm_cLISI = as.numeric(scale(-norm_cLISI)),
    z_celltype_ASW_means = as.numeric(scale(celltype_ASW_means)),
    
    # Conservative score:
    # 40% batch mixing, 60% biology preservation.
    objective_score =
      0.20 * z_norm_iLISI +
      0.30 * z_norm_cLISI +
      0.50 * z_celltype_ASW_means
  ) |>
  arrange(desc(objective_score))

write.csv(
  scores_objective,
  file.path(out_dir, "integration_metrics_objective_score.csv"),
  row.names = FALSE
)

p_objective <- scores_objective |>
  mutate(method = forcats::fct_reorder(method, objective_score)) |>
  ggplot(aes(x = method, y = objective_score, fill = is_pca)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Overall Harmony parameter score",
    subtitle = "Higher = better; score uses scaled metrics with \n20% batch / 80% biology weighting",
    x = NULL,
    y = "Objective score",
    fill = "PCA baseline"
  )

ggsave(
  file.path(out_dir, "integration_metrics_objective_score.pdf"),
  p_objective,
  width = 5,
  height = 4
)

# ----------------------------
# 5. Objective score as theta/lambda heatmap
# ----------------------------

p_objective_heat <- scores_objective |>
  filter(!is_pca) |>
  ggplot(aes(
    x = factor(theta),
    y = factor(lambda),
    fill = objective_score
  )) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(objective_score, 2)), size = 4) +
  theme_bw(base_size = 12) +
  labs(
    title = "Overall objective score by Harmony parameters",
    subtitle = "Higher = better",
    x = "theta",
    y = "lambda",
    fill = "score"
  )

ggsave(
  file.path(out_dir, "integration_metrics_objective_score_heatmap.pdf"),
  p_objective_heat,
  width = 8,
  height = 5
)

p_batch <- DimPlot(
  seu,
  reduction = "umap_harmony_theta2_lambda1",
  group.by = batch_col
) +
  ggtitle("Harmony theta=2 lambda=1 by sample")

p_annot <- DimPlot(
  seu,
  reduction = "umap_harmony_theta2_lambda1",
  group.by = annotation_col,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Harmony theta=2 lambda=1 by cell type")

ggsave(
  file.path(
    out_dir,
    "harmony_theta2_lambda1_umap_combined.pdf"
  ),
  p_batch + p_annot,
  width = 16,
  height = 6
)
####### alluvial

# ----------------------------
# Cluster cells: before (PCA) vs after (Harmony theta=2, lambda=1)
# ----------------------------

library(ggalluvial)

seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30, graph.name = c("pca_nn", "pca_snn"))
seu <- FindClusters(seu, graph.name = "pca_snn", resolution = 0.8)
seu$cluster_pca <- Idents(seu)

seu <- FindNeighbors(seu, reduction = "harmony_theta2_lambda1", dims = 1:30, graph.name = c("harmony_nn", "harmony_snn"))
seu <- FindClusters(seu, graph.name = "harmony_snn", resolution = 0.8)
seu$cluster_harmony <- Idents(seu)

alluvial_df <- seu@meta.data %>%
  dplyr::select(
    cluster_pca, cluster_harmony,
    sample = .data[[batch_col]],
    annotation = .data[[annotation_col]]
  ) %>%
  dplyr::mutate(
    cluster_pca = paste0("PCA_", cluster_pca),
    cluster_harmony = paste0("Harmony_", cluster_harmony)
  )

# ----------------------------
# Alluvial 1: PCA cluster -> Harmony cluster, coloured by sample (batch)
# ----------------------------

alluvial_sample <- alluvial_df %>%
  dplyr::count(cluster_pca, cluster_harmony, sample, name = "n")

p_alluvial_sample <- ggplot(
  alluvial_sample,
  aes(axis1 = cluster_pca, axis2 = cluster_harmony, y = n)
) +
  geom_alluvium(aes(fill = sample), width = 0.3) +
  geom_stratum(width = 0.3, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("PCA cluster", "Harmony cluster"), expand = c(.1, .1)) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Cell flow between PCA and Harmony (theta=2, lambda=1) clusters, by sample",
    subtitle = "Each ribbon is a group of cells; ribbon width = number of cells",
    y = "Number of cells",
    fill = "Sample"
  )

ggsave(
  file.path(out_dir, "alluvial_pca_to_harmony_by_sample.pdf"),
  p_alluvial_sample,
  width = 12,
  height = 8
)

# ----------------------------
# Alluvial 2: PCA cluster -> Harmony cluster, coloured by broad cell type
# ----------------------------

alluvial_annot <- alluvial_df %>%
  dplyr::count(cluster_pca, cluster_harmony, annotation, name = "n")

p_alluvial_annot <- ggplot(
  alluvial_annot,
  aes(axis1 = cluster_pca, axis2 = cluster_harmony, y = n)
) +
  geom_alluvium(aes(fill = annotation), width = 0.3) +
  geom_stratum(width = 0.3, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("PCA cluster", "Harmony cluster"), expand = c(.1, .1)) +
  scale_fill_manual(values = celltype_cols) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Cell flow between PCA and Harmony (theta=2, lambda=1) clusters, by cell type",
    subtitle = "Each ribbon is a group of cells; ribbon width = number of cells",
    y = "Number of cells",
    fill = "Cell type"
  )

ggsave(
  file.path(out_dir, "alluvial_pca_to_harmony_by_annotation.pdf"),
  p_alluvial_annot,
  width = 12,
  height = 8
)

# ----------------------------
# Re-save object with new cluster assignments
# ----------------------------

saveRDS(
  seu,
  file.path(out_dir, "demux_singlets_annotated_seurat_harmony_parameter_sweep.rds")
)
