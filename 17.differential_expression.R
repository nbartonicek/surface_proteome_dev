suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(msigdbr)
  library(fgsea)
})

set.seed(123)

# ============================================================
# Paths
# ============================================================

runs <- c(
  "260423_VH01624_453_222HWMYNX",
  "260528_VH01624_464_222K7VKNX"
)

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

out_dir <- file.path("../merged_run_1_2", "DE", "monocyte_pseudobulk_vs_normal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Parameters
# ============================================================

assay_use <- "RNA"

sample_col <- "sample_name"
celltype_col <- "predicted_CellType_Broad"
doublet_col <- "scDblFinder.class"

normal_patterns <- c("^normal", "normal-01", "normal")
exclude_patterns <- c("MOLM13")

fdr_cutoff <- 0.05
logfc_cutoff <- 1
bcv <- 0.4

# ============================================================
# Helpers
# ============================================================

guess_sample_column <- function(meta) {
  candidates <- c("sample_name", "sampleID", "orig.ident", "sample", "donor")
  hit <- candidates[candidates %in% colnames(meta)][1]
  if (is.na(hit)) stop("Could not find a sample column.")
  hit
}

make_pseudobulk_from_cells <- function(count_mat, cells) {
  cells <- intersect(cells, colnames(count_mat))
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}

run_edgeR_exact <- function(count_mat, group_vec, bcv = 0.4) {
  
  group_vec <- factor(group_vec, levels = c("normal", "aml"))
  
  dge <- edgeR::DGEList(counts = count_mat, group = group_vec)
  
  keep <- edgeR::filterByExpr(dge, group = group_vec)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) < 10) return(NULL)
  
  dge <- edgeR::calcNormFactors(dge)
  
  et <- edgeR::exactTest(
    dge,
    pair = c("normal", "aml"),
    dispersion = bcv^2
  )
  
  edgeR::topTags(et, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    as_tibble()
}

# ============================================================
# Load and merge Seurat objects
# ============================================================

seuList <- list()

for (run in runs) {
  annotation_dir <- file.path(proj, "results/seurat_annotated", run)
  seurat_file <- file.path(annotation_dir, "demux_singlets_annotated_seurat.rds")
  
  message("Reading: ", seurat_file)
  seuList[[run]] <- readRDS(seurat_file)
}

seu <- merge(
  x = seuList[[1]],
  y = seuList[[2]],
  add.cell.ids = c("run1", "run2"),
  project = "combined"
)

DefaultAssay(seu) <- assay_use

meta <- seu@meta.data

if (!sample_col %in% colnames(meta)) {
  sample_col <- guess_sample_column(meta)
  message("Using guessed sample column: ", sample_col)
}

stopifnot(celltype_col %in% colnames(meta))
stopifnot(doublet_col %in% colnames(meta))

# ============================================================
# Filter cells
# ============================================================

meta <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    sample_value = as.character(.data[[sample_col]]),
    is_normal = str_detect(
      sample_value,
      regex(paste(normal_patterns, collapse = "|"), ignore_case = TRUE)
    ),
    is_excluded = str_detect(
      sample_value,
      regex(paste(exclude_patterns, collapse = "|"), ignore_case = TRUE)
    )
  )

cells_keep <- meta %>%
  filter(
    !is_excluded,
    .data[[doublet_col]] == "singlet",
    .data[[celltype_col]] == "Monocyte",
    !is.na(sample_value)
  ) %>%
  pull(cell)

seu_mono <- subset(seu, cells = cells_keep)

message("Cells after filtering: ", ncol(seu_mono))

cell_summary <- seu_mono@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    sample_value = as.character(.data[[sample_col]]),
    is_normal = str_detect(
      sample_value,
      regex(paste(normal_patterns, collapse = "|"), ignore_case = TRUE)
    ),
    group = if_else(is_normal, "normal", "aml")
  ) %>%
  count(group, sample_value, name = "n_cells") %>%
  arrange(group, sample_value)

write_csv(cell_summary, file.path(out_dir, "monocyte_cell_counts_by_sample.csv"))
print(cell_summary)

# ============================================================
# Pseudobulk by sample
# ============================================================

# Seurat v5: merge creates multiple layers, e.g. counts.1/counts.2
if (inherits(seu_mono[[assay_use]], "Assay5")) {
  seu_mono[[assay_use]] <- JoinLayers(seu_mono[[assay_use]])
}

counts <- GetAssayData(
  seu_mono,
  assay = assay_use,
  layer = "counts"
)

meta_mono <- seu_mono@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    sample_value = as.character(.data[[sample_col]]),
    is_normal = str_detect(
      sample_value,
      regex(paste(normal_patterns, collapse = "|"), ignore_case = TRUE)
    ),
    group = if_else(is_normal, "normal", "aml")
  )

sample_ids <- unique(meta_mono$sample_value)

pb_list <- lapply(sample_ids, function(sid) {
  cells_i <- meta_mono %>%
    filter(sample_value == sid) %>%
    pull(cell)
  
  pb <- make_pseudobulk_from_cells(counts, cells_i)
  names(pb) <- rownames(counts)
  pb
})

names(pb_list) <- sample_ids

pb_counts <- do.call(cbind, pb_list)
rownames(pb_counts) <- rownames(counts)

pb_meta <- tibble(
  sample = colnames(pb_counts)
) %>%
  left_join(
    meta_mono %>%
      distinct(sample_value, group) %>%
      rename(sample = sample_value),
    by = "sample"
  )

write_csv(pb_meta, file.path(out_dir, "pseudobulk_sample_metadata.csv"))

# ============================================================
# DE: AML monocytes vs normal monocytes
# ============================================================

if (sum(pb_meta$group == "normal") < 1) {
  stop("No normal pseudobulk sample found.")
}

if (sum(pb_meta$group == "aml") < 1) {
  stop("No AML pseudobulk samples found.")
}

message("Normal pseudobulk samples: ", sum(pb_meta$group == "normal"))
message("AML pseudobulk samples: ", sum(pb_meta$group == "aml"))

de_res <- run_edgeR_exact(
  count_mat = pb_counts,
  group_vec = pb_meta$group,
  bcv = bcv
)

if (is.null(de_res)) {
  stop("edgeR failed or too few genes passed filtering.")
}

de_res <- de_res %>%
  mutate(
    direction = case_when(
      FDR < fdr_cutoff & logFC >= logfc_cutoff  ~ "Up in AML monocytes",
      FDR < fdr_cutoff & logFC <= -logfc_cutoff ~ "Down in AML monocytes",
      TRUE ~ "Not significant"
    ),
    neg_log10_fdr = -log10(pmax(FDR, 1e-300))
  )

write_csv(de_res, file.path(out_dir, "AML_monocyte_vs_normal_monocyte_edgeR_exactTest.csv"))

# ============================================================
# Volcano plot
# ============================================================

label_genes <- de_res %>%
  filter(FDR < fdr_cutoff, abs(logFC) >= logfc_cutoff) %>%
  arrange(FDR) %>%
  slice_head(n = 30)

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
  
  if ("adj.P.Val" %in% colnames(top) && !"FDR" %in% colnames(top)) {
    top <- top %>%
      dplyr::rename(FDR = adj.P.Val)
  }
  
  top <- top %>%
    dplyr::mutate(
      diff_expression = dplyr::case_when(
        logFC >= logFC_cutoff & FDR <= FDR_cutoff ~ "UP",
        logFC <= -logFC_cutoff & FDR <= FDR_cutoff ~ "DOWN",
        TRUE ~ "NO"
      ),
      labelgenes = dplyr::if_else(diff_expression != "NO", gene, "")
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
      x = "log2FC (AML / Normal)",
      y = "-log10(FDR)",
      colour = "DE"
    )
}

n_normal <- sum(pb_meta$group == "normal")
n_aml <- sum(pb_meta$group == "aml")

p_volcano <- plot_volcano(
  de_res,
  FDR_cutoff = 0.05,
  logFC_cutoff = 1,
  title = "AML Monocytes vs Normal Monocytes",
  subtitle = paste(
    paste0("AML samples = ", n_aml),
    paste0("Normal samples = ", n_normal),
    paste0("Monocyte cells = ", ncol(seu_mono)),
    "Pseudobulk edgeR exactTest",
    paste0("BCV = ", bcv),
    sep = " | "
  )
)

ggsave(
  file.path(plot_dir, "AML_monocyte_vs_normal_volcano.pdf"),
  p_volcano,
  width = 8,
  height = 7
)

ggsave(
  file.path(plot_dir, "AML_monocyte_vs_normal_volcano.png"),
  p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# Hallmark enrichment using fgsea
# ============================================================

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- hallmark_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

ranks <- de_res %>%
  filter(!is.na(logFC), !is.na(PValue)) %>%
  mutate(
    rank_stat = sign(logFC) * -log10(pmax(PValue, 1e-300))
  ) %>%
  arrange(desc(rank_stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(gene, rank_stat) %>%
  deframe()

fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = ranks,
  minSize = 10,
  maxSize = 500,
  nperm = 10000
) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    direction = if_else(NES > 0, "Enriched in AML monocytes", "Enriched in normal monocytes")
  )

write_csv(
  fgsea_res,
  file.path(out_dir, "AML_monocyte_vs_normal_Hallmark_fgsea.csv")
)

top_fgsea <- fgsea_res %>%
  filter(padj < 0.25) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  mutate(
    pathway = str_replace(pathway, "^HALLMARK_", ""),
    pathway = str_replace_all(pathway, "_", " "),
    pathway = factor(pathway, levels = rev(pathway))
  )

p_fgsea <- ggplot(
  top_fgsea,
  aes(x = NES, y = pathway, fill = direction)
) +
  geom_col() +
  theme_bw(base_size = 11) +
  labs(
    title = "MSigDB Hallmark enrichment",
    subtitle = "Ranked by signed -log10(P value) from edgeR",
    x = "Normalized enrichment score",
    y = NULL,
    fill = NULL
  )

ggsave(
  file.path(plot_dir, "AML_monocyte_vs_normal_Hallmark_fgsea_top20.pdf"),
  p_fgsea,
  width = 8,
  height = 6
)

ggsave(
  file.path(plot_dir, "AML_monocyte_vs_normal_Hallmark_fgsea_top20.png"),
  p_fgsea,
  width = 8,
  height = 6,
  dpi = 300
)

message("Done. Outputs written to: ", out_dir)
