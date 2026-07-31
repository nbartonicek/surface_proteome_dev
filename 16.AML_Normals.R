# ============================================================
# Project NORMAL AML atlas cells onto BoneMarrowMap
# Direct h5ad reader: no zellkonverter
# Reads only selected normal cells from AnnData CSR /X
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(rhdf5)
  library(Seurat)
  library(BoneMarrowMap)
  library(symphony)
  library(patchwork)
  library(RColorBrewer)
  library(curl)
})

# ============================================================
# Paths
# ============================================================

aml_h5ad <- "../results/AML_atlas/6e37b9b1-185c-4505-9842-8138157c1923.h5ad"

output_dir <- "../results/AML_atlas_bonemarrowmap_normal_projection"
metadata_dir <- "../results/depth_benchmarking"
figure_dir <- file.path(output_dir, "figures")
reference_dir <- file.path(output_dir, "reference")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper readers
# ============================================================

read_h5ad_obs_col <- function(file, col) {
  
  base <- paste0("/obs/", col)
  
  out <- tryCatch({
    codes <- rhdf5::h5read(file, paste0(base, "/codes"))
    cats <- rhdf5::h5read(file, paste0(base, "/categories"))
    
    values <- rep(NA_character_, length(codes))
    keep <- !is.na(codes) & codes >= 0
    values[keep] <- as.character(cats[codes[keep] + 1])
    values
  }, error = function(e) NULL)
  
  if (!is.null(out)) return(out)
  
  out <- tryCatch({
    rhdf5::h5read(file, base)
  }, error = function(e) NULL)
  
  if (!is.null(out)) return(as.character(out))
  
  stop("Could not read /obs column: ", col)
}

read_h5ad_var <- function(file) {
  
  var <- rhdf5::h5read(file, "/var")
  
  if (!is.list(var)) {
    stop("/var was not read as a list-like object.")
  }
  
  message("Available /var fields:")
  print(names(var))
  print(sapply(var, length))
  
  var
}

get_gene_vector <- function(var) {
  
  # Prefer gene symbols if present, otherwise Ensembl IDs.
  possible_symbol_cols <- c("gene_name", "feature_name", "symbol", "gene_symbols")
  possible_id_cols <- c("_index", "gene_id", "ensembl_id")
  
  symbol_col <- intersect(possible_symbol_cols, names(var))
  id_col <- intersect(possible_id_cols, names(var))
  
  if (length(symbol_col) > 0) {
    gene_names <- as.character(var[[symbol_col[1]]])
  } else if (length(id_col) > 0) {
    gene_names <- as.character(var[[id_col[1]]])
  } else {
    gene_names <- as.character(var[[1]])
  }
  
  gene_names <- make.unique(gene_names)
  
  if (length(gene_names) < 1000) {
    stop("gene_names looks too short: ", length(gene_names))
  }
  
  gene_names
}

read_h5ad_X_selected_as_genes_by_cells <- function(
    file,
    selected_cell_ids,
    all_cell_ids,
    gene_names
) {
  
  data <- rhdf5::h5read(file, "/X/data")
  indices <- rhdf5::h5read(file, "/X/indices")
  indptr <- rhdf5::h5read(file, "/X/indptr")
  
  data <- as.numeric(drop(data))
  indices <- as.integer(drop(indices))
  indptr <- as.integer(drop(indptr))
  
  dim(data) <- NULL
  dim(indices) <- NULL
  dim(indptr) <- NULL
  
  n_cells_all <- length(all_cell_ids)
  n_genes <- length(gene_names)
  
  message("Detected sparse AnnData CSR matrix.")
  message("Full cells from obs: ", n_cells_all)
  message("Selected cells requested: ", length(selected_cell_ids))
  message("Genes from var: ", n_genes)
  message("length(indptr): ", length(indptr))
  message("length(data): ", length(data))
  message("length(indices): ", length(indices))
  
  if (length(indptr) != n_cells_all + 1L) {
    stop(
      "indptr length does not equal full obs cells + 1. ",
      "Expected ", n_cells_all + 1L,
      ", got ", length(indptr), "."
    )
  }
  
  if (length(data) != length(indices)) {
    stop("Length of /X/data does not match length of /X/indices.")
  }
  
  if (max(indices, na.rm = TRUE) + 1L > n_genes) {
    stop(
      "Max gene index in /X/indices exceeds n_genes. ",
      "Check /var gene vector."
    )
  }
  
  cell_match <- match(selected_cell_ids, all_cell_ids)
  keep <- !is.na(cell_match)
  
  if (!any(keep)) {
    stop("None of selected_cell_ids matched all_cell_ids.")
  }
  
  selected_cell_ids <- selected_cell_ids[keep]
  cell_idx <- cell_match[keep]
  
  message("Selected cells matched: ", length(cell_idx))
  
  out_i <- integer()
  out_j <- integer()
  out_x <- numeric()
  
  for (k in seq_along(cell_idx)) {
    
    r <- cell_idx[k]
    
    start <- indptr[r] + 1L
    end <- indptr[r + 1L]
    
    if (end >= start) {
      pos <- start:end
      
      out_i <- c(out_i, rep(k, length(pos)))
      out_j <- c(out_j, indices[pos] + 1L)
      out_x <- c(out_x, data[pos])
    }
  }
  
  X_cells_by_genes <- Matrix::sparseMatrix(
    i = out_i,
    j = out_j,
    x = out_x,
    dims = c(length(cell_idx), n_genes)
  )
  
  rownames(X_cells_by_genes) <- selected_cell_ids
  colnames(X_cells_by_genes) <- gene_names
  
  X_genes_by_cells <- Matrix::t(X_cells_by_genes)
  X_genes_by_cells <- as(X_genes_by_cells, "dgCMatrix")
  
  X_genes_by_cells
}

# ============================================================
# Read obs metadata
# ============================================================

obs_cols <- c(
  "_index",
  "Study",
  "Sample",
  "donor_id",
  "cell_type",
  "Author Cell Type",
  "HSPC Cell Type",
  "disease",
  "assay",
  "sex",
  "age_floor",
  "Cytogenetics",
  "ELN Classification",
  "ELN Risk Group",
  "Mutations",
  "Translocations"
)

obs_list <- purrr::map(obs_cols, function(x) {
  message("Reading obs column: ", x)
  read_h5ad_obs_col(aml_h5ad, x)
})

names(obs_list) <- make.names(obs_cols)

obs <- tibble::as_tibble(obs_list) %>%
  dplyr::rename(cell_id = X_index)

write.csv(
  obs,
  file.path(metadata_dir, "aml_atlas_obs_metadata_selected.csv"),
  row.names = FALSE
)

message("Obs dimensions:")
print(dim(obs))

message("Disease table:")
print(table(obs$disease, useNA = "ifany"))

# ============================================================
# Select normal / healthy / control cells
# ============================================================

normal_obs <- obs %>%
  dplyr::filter(
    stringr::str_detect(
      stringr::str_to_lower(disease),
      "normal|healthy|control"
    )
  ) %>%
  dplyr::distinct(cell_id, .keep_all = TRUE)

if (nrow(normal_obs) == 0) {
  stop("No normal cells found. Inspect table(obs$disease) and adjust filter.")
}

normal_cells <- normal_obs$cell_id

message("Normal cells selected: ", length(normal_cells))
message("Normal donors selected: ", dplyr::n_distinct(normal_obs$donor_id))
print(table(normal_obs$disease, useNA = "ifany"))

write.csv(
  normal_obs,
  file.path(output_dir, "aml_atlas_normal_obs_metadata.csv"),
  row.names = FALSE
)

# ============================================================
# Read genes
# ============================================================

#var <- read_h5ad_var(aml_h5ad)
#gene_names <- get_gene_vector(var)



gene_ids <- as.character(var$`_index`)
gene_names <- as.character(var$gene_symbol)

gene_names[is.na(gene_names) | gene_names == ""] <-
  gene_ids[is.na(gene_names) | gene_names == ""]

gene_names <- make.unique(gene_names)

length(gene_names)
head(gene_names)

message("Genes read: ", length(gene_names))
message("First genes:")
print(head(gene_names))

# ============================================================
# Read only normal cells from /X
# ============================================================
use_python(
  "/Users/bartoniceknenad/miniconda3/bin/python",
  required = TRUE
)

py_config()

py_module_available("anndata")
py_module_available("scipy")

expr <- read_h5ad_X_selected_as_genes_by_cells(
  file = aml_h5ad,
  selected_cell_ids = normal_cells,
  all_cell_ids = obs$cell_id,
  gene_names = gene_names
)

saveRDS(expr,file="../results/normal_cells.rds")
message("Expression matrix genes x normal cells:")
print(dim(expr))

if (nrow(expr) != length(gene_names)) {
  stop("Number of rows in expression matrix does not match number of genes.")
}

normal_cells_present <- intersect(normal_cells, colnames(expr))

message("Normal cells present in X: ", length(normal_cells_present))

if (length(normal_cells_present) == 0) {
  stop("No normal cells matched colnames of expression matrix.")
}

expr <- expr[, normal_cells_present, drop = FALSE]
expr <- as(expr, "dgCMatrix")

# ============================================================
# Metadata aligned to expression matrix
# ============================================================

meta_query <- normal_obs %>%
  dplyr::filter(cell_id %in% colnames(expr)) %>%
  dplyr::distinct(cell_id, .keep_all = TRUE)

meta_query <- as.data.frame(meta_query)
rownames(meta_query) <- meta_query$cell_id
meta_query <- meta_query[colnames(expr), , drop = FALSE]

if (!identical(rownames(meta_query), colnames(expr))) {
  stop("Metadata rows do not match expression matrix columns.")
}

meta_query <- meta_query %>%
  dplyr::mutate(
    projection_donor = donor_id,
    projection_sample = Sample,
    projection_disease = disease
  )

message("Final normal expression matrix genes x cells:")
print(dim(expr))

message("Final normal metadata:")
print(dim(meta_query))

# ============================================================
# Optional: cap cells per donor
# ============================================================

max_cells_per_donor <- Inf

if (is.finite(max_cells_per_donor)) {
  
  set.seed(1)
  
  keep_cells <- meta_query %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::group_by(projection_donor) %>%
    dplyr::slice_sample(n = min(dplyr::n(), max_cells_per_donor)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell)
  
  expr <- expr[, keep_cells, drop = FALSE]
  meta_query <- meta_query[keep_cells, , drop = FALSE]
}

# ============================================================
# Load BoneMarrowMap reference
# ============================================================

bm_ref_rds <- file.path(reference_dir, "BoneMarrowMap_SymphonyReference.rds")
bm_uwot <- file.path(reference_dir, "BoneMarrowMap_uwot_model.uwot")

if (!file.exists(bm_ref_rds)) {
  curl::curl_download(
    "https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds",
    bm_ref_rds
  )
}

if (!file.exists(bm_uwot)) {
  curl::curl_download(
    "https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_uwot_model.uwot",
    bm_uwot
  )
}

BM_ref <- readRDS(bm_ref_rds)
BM_ref$save_uwot_path <- bm_uwot

# ============================================================
# Project onto BoneMarrowMap
# ============================================================

query <- map_Query(
  exp_query = expr,
  metadata_query = meta_query,
  ref_obj = BM_ref,
  vars = "projection_donor"
)

query <- calculate_MappingError(
  query,
  reference = BM_ref,
  MAD_threshold = 2.5,
  threshold_by_donor = TRUE,
  donor_key = "projection_donor"
)

query <- predict_CellTypes(
  query_obj = query,
  ref_obj = BM_ref,
  initial_label = "initial_CellType_BoneMarrowMap",
  final_label = "predicted_CellType_BoneMarrowMap"
)

query <- predict_Pseudotime(
  query_obj = query,
  ref_obj = BM_ref,
  initial_label = "initial_Pseudotime",
  final_label = "predicted_Pseudotime"
)

saveRDS(
  query,
  file.path(output_dir, "aml_atlas_normal_BoneMarrowMap_projected.rds")
)

# ============================================================
# Save metadata / labels
# ============================================================

projection_meta <- query@meta.data %>%
  tibble::rownames_to_column("cell_id")

write.csv(
  projection_meta,
  file.path(output_dir, "aml_atlas_normal_BoneMarrowMap_projected_metadata.csv"),
  row.names = FALSE
)

save_ProjectionResults(
  query_obj = query,
  file_name = file.path(output_dir, "aml_atlas_normal_BoneMarrowMap_projected_labels.csv")
)

# ============================================================
# Composition summaries
# ============================================================

normal_composition <- projection_meta %>%
  dplyr::count(
    projection_donor,
    predicted_CellType_BoneMarrowMap,
    name = "n"
  ) %>%
  dplyr::group_by(projection_donor) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

write.csv(
  normal_composition,
  file.path(output_dir, "normal_donor_celltype_composition.csv"),
  row.names = FALSE
)

normal_broad_composition <- projection_meta %>%
  dplyr::count(
    projection_donor,
    predicted_CellType_BoneMarrowMap_Broad,
    name = "n"
  ) %>%
  dplyr::group_by(projection_donor) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

write.csv(
  normal_broad_composition,
  file.path(output_dir, "normal_donor_broad_celltype_composition.csv"),
  row.names = FALSE
)

# ============================================================
# Plots
# ============================================================

query_pass <- subset(query, mapping_error_QC == "Pass")

p_donor <- DimPlot(
  query_pass,
  reduction = "umap_projected",
  group.by = "projection_donor",
  raster = FALSE
)

p_celltype <- DimPlot(
  query_pass,
  reduction = "umap_projected",
  group.by = "predicted_CellType_BoneMarrowMap",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + NoLegend()

p_broad <- DimPlot(
  query_pass,
  reduction = "umap_projected",
  group.by = "predicted_CellType_BoneMarrowMap_Broad",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + NoLegend()

p_pseudotime <- FeaturePlot(
  query_pass,
  reduction = "umap_projected",
  features = "predicted_Pseudotime",
  raster = FALSE
) +
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )

ggsave(file.path(figure_dir, "normal_projection_by_donor.pdf"), p_donor, width = 9, height = 7)
ggsave(file.path(figure_dir, "normal_projection_predicted_celltype.pdf"), p_celltype, width = 10, height = 8)
ggsave(file.path(figure_dir, "normal_projection_predicted_broad_celltype.pdf"), p_broad, width = 9, height = 7)
ggsave(file.path(figure_dir, "normal_projection_pseudotime.pdf"), p_pseudotime, width = 8, height = 6)

p_mapping_error <- query@meta.data %>%
  ggplot(aes(x = mapping_error_score, fill = mapping_error_QC)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ projection_donor, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "BoneMarrowMap mapping error QC",
    x = "Mapping error score",
    y = "Cells"
  )

ggsave(
  file.path(figure_dir, "normal_projection_mapping_error_by_donor.pdf"),
  p_mapping_error,
  width = 12,
  height = 8
)

p_broad_comp <- normal_broad_composition %>%
  ggplot(aes(
    x = projection_donor,
    y = prop,
    fill = predicted_CellType_BoneMarrowMap_Broad
  )) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    title = "Normal AML atlas donor composition after BoneMarrowMap projection",
    x = "Donor",
    y = "Proportion",
    fill = "Broad cell type"
  )

ggsave(
  file.path(figure_dir, "normal_donor_broad_celltype_composition.pdf"),
  p_broad_comp,
  width = 11,
  height = 9
)

message("Done. Outputs written to: ", output_dir)