# ============================================================
# Project subsetted normal AML atlas onto BoneMarrowMap
# 5000 cells per donor
# Uses rhdf5 instead of zellkonverter/readH5AD
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(Seurat)
  library(BoneMarrowMap)
  library(symphony)
  library(curl)
  library(RColorBrewer)
  library(rhdf5)
})

output_dir <- "../results/AML_atlas_bonemarrowmap_normal_projection"
figure_dir <- file.path(output_dir, "figures_5000")
reference_dir <- file.path(output_dir, "reference")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

h5ad_file <- file.path(
  output_dir,
  "aml_atlas_normal_5000_per_donor.h5ad"
)

if (!file.exists(h5ad_file)) {
  stop("Could not find h5ad file: ", h5ad_file)
}

# ============================================================
# rhdf5 h5ad reader helpers
# ============================================================

# h5ls() has NO "group path" argument. Its signature is:
#   h5ls(file, recursive = TRUE, all = FALSE, datasetinfo = TRUE, ...)
# Calling h5ls(file, "/var", recursive = FALSE) silently matches "/var"
# positionally to the `all` argument (since `recursive` was already
# matched by name) -- which is how you get "Error in !all : invalid
# argument type" deep inside rhdf5. The fix is to always do a full
# listing and filter by the `group` column ourselves.
h5ls_at <- function(h5ad_file, group_path) {
  full <- h5ls(h5ad_file, recursive = TRUE)
  full[full$group == group_path, , drop = FALSE]
}

h5_path_exists <- function(h5ad_file, path) {
  all_paths <- paste0(h5ls(h5ad_file, recursive = TRUE)$group, "/", h5ls(h5ad_file, recursive = TRUE)$name)
  path %in% all_paths
}

read_h5ad_var_names <- function(h5ad_file) {
  
  var_items <- h5ls_at(h5ad_file, "/var")
  
  if ("gene_symbol" %in% var_items$name) {
    genes <- as.character(as.vector(h5read(h5ad_file, "/var/gene_symbol")))
  } else if ("features" %in% var_items$name) {
    genes <- as.character(as.vector(h5read(h5ad_file, "/var/features")))
  } else if ("_index" %in% var_items$name) {
    genes <- as.character(as.vector(h5read(h5ad_file, "/var/_index")))
  } else if ("index" %in% var_items$name) {
    genes <- as.character(as.vector(h5read(h5ad_file, "/var/index")))
  } else {
    stop(
      "Could not find gene names in /var. Available /var fields: ",
      paste(var_items$name, collapse = ", ")
    )
  }
  
  make.unique(genes)
}

read_h5ad_obs_names <- function(h5ad_file) {
  
  obs_items <- h5ls_at(h5ad_file, "/obs")
  
  if ("_index" %in% obs_items$name) {
    return(as.character(as.vector(h5read(h5ad_file, "/obs/_index"))))
  }
  
  if ("index" %in% obs_items$name) {
    return(as.character(as.vector(h5read(h5ad_file, "/obs/index"))))
  }
  
  stop(
    "Could not find cell names in /obs. Available /obs fields: ",
    paste(obs_items$name, collapse = ", ")
  )
}

read_h5ad_obs <- function(h5ad_file) {
  
  message("Reading /obs metadata from h5ad...")
  
  obs_names <- read_h5ad_obs_names(h5ad_file)
  n_cells <- length(obs_names)
  
  obs_items <- h5ls_at(h5ad_file, "/obs")
  
  meta <- data.frame(
    cell_id = obs_names,
    row.names = obs_names,
    stringsAsFactors = FALSE
  )
  
  for (nm in obs_items$name) {
    
    if (nm %in% c("_index", "index")) next
    
    path <- paste0("/obs/", nm)
    
    # AnnData categorical: /obs/field/categories + /obs/field/codes
    cat_try <- tryCatch({
      
      sub_items <- h5ls_at(h5ad_file, path)
      
      if (all(c("categories", "codes") %in% sub_items$name)) {
        
        categories <- as.character(as.vector(h5read(h5ad_file, paste0(path, "/categories"))))
        codes <- as.integer(as.vector(h5read(h5ad_file, paste0(path, "/codes"))))
        
        out <- rep(NA_character_, length(codes))
        ok <- codes >= 0
        out[ok] <- categories[codes[ok] + 1]
        
        if (length(out) == n_cells) out else NULL
        
      } else {
        NULL
      }
      
    }, error = function(e) NULL)
    
    if (!is.null(cat_try)) {
      meta[[nm]] <- cat_try
      next
    }
    
    # Plain vector column
    vec_try <- tryCatch({
      val <- h5read(h5ad_file, path)
      
      if (is.list(val)) return(NULL)
      
      val <- as.vector(val)
      
      if (length(val) != n_cells) return(NULL)
      
      if (is.raw(val)) {
        val <- as.character(val)
      }
      
      val
      
    }, error = function(e) NULL)
    
    if (!is.null(vec_try)) {
      meta[[nm]] <- vec_try
    }
  }
  
  meta
}

read_h5ad_X_as_genes_by_cells <- function(h5ad_file) {
  
  message("Reading /X sparse matrix from h5ad using rhdf5...")
  
  x_items <- h5ls_at(h5ad_file, "/X")
  
  message("Contents of /X:")
  print(x_items)
  
  needed <- c("data", "indices", "indptr")
  missing_x <- setdiff(needed, x_items$name)
  
  if (length(missing_x) > 0) {
    stop(
      "/X is not stored as sparse data/indices/indptr. Missing: ",
      paste(missing_x, collapse = ", "),
      "\nFound under /X: ",
      paste(x_items$name, collapse = ", ")
    )
  }
  
  data <- as.numeric(as.vector(h5read(h5ad_file, "/X/data")))
  indices <- as.integer(as.vector(h5read(h5ad_file, "/X/indices")))
  indptr <- as.integer(as.vector(h5read(h5ad_file, "/X/indptr")))
  
  obs_names <- read_h5ad_obs_names(h5ad_file)
  gene_names <- read_h5ad_var_names(h5ad_file)
  
  n_cells <- length(obs_names)
  n_genes <- length(gene_names)
  
  message("Detected cells: ", n_cells)
  message("Detected genes: ", n_genes)
  message("length(data): ", length(data))
  message("length(indices): ", length(indices))
  message("length(indptr): ", length(indptr))
  
  if (length(data) != length(indices)) {
    stop("Length mismatch: /X/data and /X/indices are different lengths.")
  }
  
  if (length(indptr) == n_cells + 1) {
    
    message("Detected CSR layout: cells x genes")
    
    X_cells_by_genes <- Matrix::sparseMatrix(
      j = indices + 1L,
      p = indptr,
      x = data,
      dims = c(n_cells, n_genes),
      index1 = TRUE
    )
    
    rownames(X_cells_by_genes) <- obs_names
    colnames(X_cells_by_genes) <- gene_names
    
    expr <- Matrix::t(X_cells_by_genes)
    
  } else if (length(indptr) == n_genes + 1) {
    
    message("Detected CSC layout: genes x cells")
    
    expr <- Matrix::sparseMatrix(
      i = indices + 1L,
      p = indptr,
      x = data,
      dims = c(n_genes, n_cells),
      index1 = TRUE
    )
    
    rownames(expr) <- gene_names
    colnames(expr) <- obs_names
    
  } else {
    
    stop(
      "Cannot infer sparse layout.\n",
      "length(indptr) = ", length(indptr), "\n",
      "n_cells + 1 = ", n_cells + 1, "\n",
      "n_genes + 1 = ", n_genes + 1
    )
  }
  
  expr <- as(expr, "dgCMatrix")
  expr
}

# ============================================================
# Expression matrix: genes x cells
# ============================================================

expr <- read_h5ad_X_as_genes_by_cells(h5ad_file)

stopifnot(!is.null(rownames(expr)))
stopifnot(!is.null(colnames(expr)))

message("Expression matrix genes x cells:")
print(dim(expr))

message("First genes:")
print(head(rownames(expr)))

keep_genes <- Matrix::rowSums(expr != 0) > 0
expr <- expr[keep_genes, , drop = FALSE]

message("Expression matrix after dropping zero genes:")
print(dim(expr))

# ============================================================
# Metadata
# ============================================================

meta_query <- read_h5ad_obs(h5ad_file)

meta_query <- meta_query[colnames(expr), , drop = FALSE]

if (!identical(rownames(meta_query), colnames(expr))) {
  stop("Metadata rows do not match expression matrix columns.")
}

required_cols <- c("donor_id", "Sample", "disease")
missing_cols <- setdiff(required_cols, colnames(meta_query))

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable metadata columns: ",
    paste(colnames(meta_query), collapse = ", ")
  )
}

meta_query <- meta_query %>%
  dplyr::mutate(
    projection_donor = donor_id,
    projection_sample = Sample,
    projection_disease = disease
  )

message("Metadata:")
print(dim(meta_query))

# count checks
donor_counts <- meta_query %>%
  dplyr::count(projection_donor, name = "n_cells") %>%
  dplyr::arrange(projection_donor)

sample_counts <- meta_query %>%
  dplyr::count(projection_sample, projection_donor, name = "n_cells") %>%
  dplyr::arrange(projection_donor, projection_sample)

message("Cells per donor:")
print(donor_counts)

message("Cells per sample/library:")
print(sample_counts)

write.csv(
  donor_counts,
  file.path(output_dir, "aml_atlas_normal_5000_donor_counts.csv"),
  row.names = FALSE
)

write.csv(
  sample_counts,
  file.path(output_dir, "aml_atlas_normal_5000_sample_counts.csv"),
  row.names = FALSE
)

saveRDS(
  expr,
  file.path(output_dir, "aml_atlas_normal_5000_expr_genes_by_cells.rds")
)

saveRDS(
  meta_query,
  file.path(output_dir, "aml_atlas_normal_5000_metadata.rds")
)

# ============================================================
# Optional Seurat object
# ============================================================

seu_normal <- CreateSeuratObject(
  counts = expr,
  meta.data = meta_query,
  project = "AML_atlas_normal_5000"
)

saveRDS(
  seu_normal,
  file.path(output_dir, "aml_atlas_normal_5000_seurat.rds")
)

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
  file.path(output_dir, "aml_atlas_normal_5000_BoneMarrowMap_projected.rds")
)

# ============================================================
# Save projection metadata
# ============================================================

projection_meta <- query@meta.data %>%
  tibble::rownames_to_column("cell_id")

write.csv(
  projection_meta,
  file.path(output_dir, "aml_atlas_normal_5000_BoneMarrowMap_projected_metadata.csv"),
  row.names = FALSE
)

save_ProjectionResults(
  query_obj = query,
  file_name = file.path(output_dir, "aml_atlas_normal_5000_BoneMarrowMap_projected_labels.csv")
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

p_sample <- DimPlot(
  query_pass,
  reduction = "umap_projected",
  group.by = "projection_sample",
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

ggsave(
  file.path(figure_dir, "normal_5000_projection_by_donor.pdf"),
  p_donor,
  width = 9,
  height = 7
)

ggsave(
  file.path(figure_dir, "normal_5000_projection_by_sample.pdf"),
  p_sample,
  width = 9,
  height = 7
)

ggsave(
  file.path(figure_dir, "normal_5000_projection_predicted_celltype.pdf"),
  p_celltype,
  width = 10,
  height = 8
)

ggsave(
  file.path(figure_dir, "normal_5000_projection_predicted_broad_celltype.pdf"),
  p_broad,
  width = 9,
  height = 7
)

ggsave(
  file.path(figure_dir, "normal_5000_projection_pseudotime.pdf"),
  p_pseudotime,
  width = 8,
  height = 6
)

message("Done. Outputs written to: ", output_dir)