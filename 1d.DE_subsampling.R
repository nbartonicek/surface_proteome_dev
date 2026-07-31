# ============================================================
# Read/UMI downsampling pseudobulk DE analysis
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(edgeR)
  library(Matrix)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(DropletUtils)
})

set.seed(123)

# ============================================================
# PATHS
# ============================================================

run <- "260423_VH01624_453_222HWMYNX"

annotation_dir <- file.path("../results/seurat_annotated", run)

seurat_file <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

out_dir <- file.path(
  "../results/data_integration",
  run,
  "read_downsampled_pseudobulk_DE"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PARAMETERS
# ============================================================

assay_use <- "RNA"
sample_col <- "sampleID"
celltype_col <- "predicted_CellType_Broad"

downsample_props <- seq(0.1, 1.0, by = 0.1)
n_downsample_reps <- 5

fdr_cutoff <- 0.05
logfc_cutoff <- 1

comparisons <- tibble::tribble(
  ~sampleID,          ~comparison_name,        ~group_a,        ~group_b,
  "HBDN392-AML-MDS",  "pDC_vs_cDC",            "cDC",          "pDC",
  "HBDN392-AML-MDS",  "NaiveT_vs_CD4MemoryT",  "Naive T",      "CD8 Memory T",
  "HBDN206-MNpCT",    "HSCMPP_vs_LateGMP",     "HSC MPP",      "Late GMP",
  "HBDN206-MNpCT",    "NaiveT_vs_CD4MemoryT",  "Naive T",      "CD8 Memory T"
)

# ============================================================
# LOAD DATA
# ============================================================

seu <- readRDS(seurat_file)
DefaultAssay(seu) <- assay_use

counts <- GetAssayData(seu, assay = assay_use, layer = "counts")

if (is.null(counts) || ncol(counts) == 0) {
  counts <- GetAssayData(seu, assay = assay_use, slot = "counts")
}

meta <- seu@meta.data %>%
  rownames_to_column("cell")

stopifnot(sample_col %in% colnames(meta))
stopifnot(celltype_col %in% colnames(meta))

# ============================================================
# CHECK CELL NUMBERS
# ============================================================

cell_counts <- meta %>%
  filter(.data[[sample_col]] %in% unique(comparisons$sampleID)) %>%
  mutate(
    sampleID = .data[[sample_col]],
    celltype = .data[[celltype_col]]
  ) %>%
  dplyr::count(sampleID, celltype, name = "n_cells") %>%
  arrange(sampleID, desc(n_cells))

write_csv(cell_counts, file.path(out_dir, "cell_counts_by_sample_celltype.csv"))
print(cell_counts)

# ============================================================
# FUNCTIONS
# ============================================================

make_pseudobulk_from_cells <- function(count_mat, cells) {
  
  cells <- intersect(cells, colnames(count_mat))
  
  if (length(cells) == 0) {
    stop("No matching cells found in count matrix.")
  }
  
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}


downsample_pseudobulk <- function(pb_counts, prop) {
  
  pb_mat <- Matrix::Matrix(
    pb_counts,
    ncol = 1,
    sparse = TRUE
  )
  
  rownames(pb_mat) <- names(pb_counts)
  colnames(pb_mat) <- "sample"
  
  downsampled <- DropletUtils::downsampleMatrix(
    pb_mat,
    prop = prop
  )
  
  as.numeric(downsampled[, 1])
}


run_edgeR_de <- function(count_mat, group_vec) {
  
  group_vec <- factor(group_vec, levels = c("group_a", "group_b"))
  
  dge <- DGEList(counts = count_mat, group = group_vec)
  
  keep <- filterByExpr(dge, group = group_vec)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) < 10) {
    return(NULL)
  }
  
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ 0 + group_vec)
  colnames(design) <- levels(group_vec)
  
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  contrast <- makeContrasts(
    group_b_vs_group_a = group_b - group_a,
    levels = design
  )
  
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  topTags(qlf, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    as_tibble()
}


run_read_downsampled_comparison <- function(sample_i, comparison_name, group_a_i, group_b_i) {
  
  message("Running: ", sample_i, " | ", comparison_name)
  
  meta_i <- meta %>%
    filter(
      .data[[sample_col]] == sample_i,
      .data[[celltype_col]] %in% c(group_a_i, group_b_i)
    )
  
  cells_a <- meta_i %>%
    filter(.data[[celltype_col]] == group_a_i) %>%
    pull(cell)
  
  cells_b <- meta_i %>%
    filter(.data[[celltype_col]] == group_b_i) %>%
    pull(cell)
  
  message("  ", group_a_i, ": ", length(cells_a), " cells")
  message("  ", group_b_i, ": ", length(cells_b), " cells")
  
  if (length(cells_a) == 0 || length(cells_b) == 0) {
    warning("Missing cells for comparison: ", comparison_name)
    return(tibble())
  }
  
  pb_a_full <- make_pseudobulk_from_cells(counts, cells_a)
  pb_b_full <- make_pseudobulk_from_cells(counts, cells_b)
  
  names(pb_a_full) <- rownames(counts)
  names(pb_b_full) <- rownames(counts)
  
  full_lib_a <- sum(pb_a_full)
  full_lib_b <- sum(pb_b_full)
  
  message("  Full pseudobulk library sizes:")
  message("    ", group_a_i, ": ", full_lib_a)
  message("    ", group_b_i, ": ", full_lib_b)
  
  summary_res <- purrr::map_dfr(downsample_props, function(prop_i) {
    
    purrr::map_dfr(seq_len(n_downsample_reps), function(rep_i) {
      
      set.seed(1000 + rep_i + round(prop_i * 1000))
      
      pb_a_ds <- downsample_pseudobulk(pb_a_full, prop_i)
      pb_b_ds <- downsample_pseudobulk(pb_b_full, prop_i)
      
      names(pb_a_ds) <- rownames(counts)
      names(pb_b_ds) <- rownames(counts)
      
      count_mat <- cbind(
        group_a = pb_a_ds,
        group_b = pb_b_ds
      )
      
      rownames(count_mat) <- rownames(counts)
      
      group_vec <- c("group_a", "group_b")
      
      # edgeR needs replication for dispersion estimation.
      # With only one pseudobulk per group, we use exactTest with a fixed dispersion.
      dge <- DGEList(counts = count_mat, group = group_vec)
      
      keep <- filterByExpr(dge, group = group_vec)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- calcNormFactors(dge)
      
      bcv <- 0.4
      et <- exactTest(dge, dispersion = bcv^2)
      
      de_res <- topTags(et, n = Inf)$table %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        mutate(
          sampleID = sample_i,
          comparison = comparison_name,
          group_a = group_a_i,
          group_b = group_b_i,
          downsample_prop = prop_i,
          downsample_percent = prop_i * 100,
          downsample_rep = rep_i,
          full_lib_size_group_a = full_lib_a,
          full_lib_size_group_b = full_lib_b,
          expected_lib_size_group_a = full_lib_a * prop_i,
          expected_lib_size_group_b = full_lib_b * prop_i,
          direction = case_when(
            FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Up in group_b",
            FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Down in group_b",
            TRUE ~ "Not significant"
          )
        )
      
      write_csv(
        de_res,
        file.path(
          out_dir,
          paste0(
            sample_i, "_",
            comparison_name, "_",
            "downsample_",
            prop_i * 100,
            "pct_rep",
            rep_i,
            "_DE_results.csv"
          )
        )
      )
      
      tibble(
        sampleID = sample_i,
        comparison = comparison_name,
        group_a = group_a_i,
        group_b = group_b_i,
        downsample_prop = prop_i,
        downsample_percent = prop_i * 100,
        downsample_rep = rep_i,
        full_lib_size_group_a = full_lib_a,
        full_lib_size_group_b = full_lib_b,
        expected_lib_size_group_a = full_lib_a * prop_i,
        expected_lib_size_group_b = full_lib_b * prop_i,
        n_up = sum(de_res$direction == "Up in group_b", na.rm = TRUE),
        n_down = sum(de_res$direction == "Down in group_b", na.rm = TRUE),
        n_de_total = n_up + n_down
      )
    })
  })
  
  summary_res
}

# ============================================================
# RUN ALL COMPARISONS
# ============================================================

summary_results <- pmap_dfr(
  comparisons,
  function(sampleID, comparison_name, group_a, group_b) {
    run_read_downsampled_comparison(
      sample_i = sampleID,
      comparison_name = comparison_name,
      group_a_i = group_a,
      group_b_i = group_b
    )
  }
)

write_csv(
  summary_results,
  file.path(out_dir, "read_downsampled_pseudobulk_DE_summary.csv")
)

saveRDS(
  summary_results,
  file.path(out_dir, "read_downsampled_pseudobulk_DE_summary.rds")
)

# ============================================================
# PLOT UP / DOWN DE GENES
# ============================================================

plot_df <- summary_results %>%
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
    x = downsample_percent,
    y = n_genes,
    colour = direction,
    group = interaction(direction, downsample_rep)
  )
) +
  geom_line(alpha = 0.35) +
  geom_point(alpha = 0.5, size = 1.8) +
  stat_summary(
    aes(group = direction),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(group = direction),
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  facet_wrap(~ comparison_label, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Read/UMI downsampling pseudobulk DE saturation",
    subtitle = paste0(
      "DropletUtils::downsampleMatrix; edgeR exactTest; FDR < ",
      fdr_cutoff,
      ", |log2FC| > ",
      logfc_cutoff
    ),
    x = "Retained reads / UMIs (%)",
    y = "Number of DE genes",
    colour = "Direction"
  )

ggsave(
  file.path(out_dir, "read_downsampled_pseudobulk_DE_up_down.pdf"),
  p_de,
  width = 12,
  height = 8
)

ggsave(
  file.path(out_dir, "read_downsampled_pseudobulk_DE_up_down.png"),
  p_de,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# PLOT TOTAL DE GENES
# ============================================================

p_total <- summary_results %>%
  mutate(
    comparison_label = paste0(
      sampleID,
      "\n",
      group_b,
      " vs ",
      group_a
    )
  ) %>%
  ggplot(
    aes(
      x = downsample_percent,
      y = n_de_total,
      group = downsample_rep
    )
  ) +
  geom_line(alpha = 0.35) +
  geom_point(alpha = 0.5, size = 1.8) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  facet_wrap(~ comparison_label, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Total DE genes after read/UMI downsampling",
    x = "Retained reads / UMIs (%)",
    y = "Total DE genes"
  )

ggsave(
  file.path(out_dir, "read_downsampled_pseudobulk_DE_total.pdf"),
  p_total,
  width = 12,
  height = 8
)

ggsave(
  file.path(out_dir, "read_downsampled_pseudobulk_DE_total.png"),
  p_total,
  width = 12,
  height = 8,
  dpi = 300
)

message("Done. Outputs written to: ", out_dir)
