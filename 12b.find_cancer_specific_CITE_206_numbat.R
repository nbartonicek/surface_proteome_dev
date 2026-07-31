suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

sample_to_use <- "HBDN206-MNpCT"

numbat_dir <- file.path(
  "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/results/seurat_annotated/260423_VH01624_453_222HWMYNX/numbat",
  paste0("LK1_", sample_to_use),
  "numbat_final"
)

out_dir <- file.path(
  projection_dir,
  paste0("HSC_MPP_numbat_DSB_", sample_to_use)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD NUMBAT CALLS
# ============================================================

numbat_files <- list.files(
  numbat_dir,
  pattern = "clone_post|cell_anno|numbat.*cell|posterior|subclone",
  full.names = TRUE,
  recursive = TRUE
)

print(numbat_files)

# Most useful file is usually clone_post_*.tsv
clone_file <- numbat_files[grepl("clone_post", basename(numbat_files))][1]

if (is.na(clone_file)) {
  stop("Could not find clone_post file in: ", numbat_dir)
}

numbat_meta <- fread(clone_file)

cat("Loaded Numbat file:\n", clone_file, "\n")
print(colnames(numbat_meta))
head(numbat_meta)

# ============================================================
# STANDARDISE NUMBAT COLUMNS
# ============================================================

# cell column
if (!"cell" %in% colnames(numbat_meta)) {
  colnames(numbat_meta)[1] <- "cell"
}

# Try to find malignancy / clone-like columns
possible_malignancy_cols <- c(
  "malignant", "malignancy", "numbat_malignancy",
  "tumor", "clone", "clone_opt", "clone_post", "cnv_state"
)

available_malignancy_cols <- intersect(
  possible_malignancy_cols,
  colnames(numbat_meta)
)

cat("Candidate Numbat state columns:\n")
print(available_malignancy_cols)

# If clone column exists, use it.
# Numbat often has clone/subclone posterior info rather than simple tumor/normal labels.
if ("clone" %in% colnames(numbat_meta)) {
  numbat_meta <- numbat_meta %>%
    mutate(numbat_clone = as.character(clone))
} else if ("clone_opt" %in% colnames(numbat_meta)) {
  numbat_meta <- numbat_meta %>%
    mutate(numbat_clone = as.character(clone_opt))
} else {
  numbat_meta <- numbat_meta %>%
    mutate(numbat_clone = "Unknown")
}

# Define malignant/normal.
# Adjust this if your Numbat output uses different labels.
numbat_meta <- numbat_meta %>%
  mutate(
    numbat_malignancy = case_when(
      numbat_clone %in% c("normal", "diploid", "neutral", "0") ~ "CNV_neutral",
      is.na(numbat_clone) | numbat_clone == "Unknown" ~ "Not_called",
      TRUE ~ "CNV_aberrant"
    )
  )

table(numbat_meta$numbat_malignancy, useNA = "ifany")
table(numbat_meta$numbat_clone, useNA = "ifany")

# ============================================================
# SUBSET SAMPLE
# ============================================================

seu_sub <- subset(
  seu,
  subset = sampleID == sample_to_use
)

cat("Sample cells:", ncol(seu_sub), "\n")

# ============================================================
# ADD NUMBAT METADATA
# ============================================================

seu_sub$numbat_clone <- NA_character_
seu_sub$numbat_malignancy <- "Not_called"

common_cells <- intersect(colnames(seu_sub), numbat_meta$cell)

idx_seu <- match(common_cells, colnames(seu_sub))
idx_num <- match(common_cells, numbat_meta$cell)

seu_sub$numbat_clone[idx_seu] <- numbat_meta$numbat_clone[idx_num]
seu_sub$numbat_malignancy[idx_seu] <- numbat_meta$numbat_malignancy[idx_num]

cat("Numbat annotation overlap:", length(common_cells), "\n")
table(seu_sub$numbat_malignancy, useNA = "ifany")
table(seu_sub$numbat_clone, useNA = "ifany")

# ============================================================
# HSC-MPP ONLY
# ============================================================

hsc <- subset(
  seu_sub,
  subset = predicted_CellType_Broad == "HSC MPP"
)

cat("HSC-MPP cells:", ncol(hsc), "\n")
table(hsc$numbat_malignancy, useNA = "ifany")
table(hsc$numbat_clone, useNA = "ifany")

hsc <- subset(
  hsc,
  subset = numbat_malignancy %in% c("CNV_aberrant", "CNV_neutral")
)

table(hsc$numbat_malignancy)

# Stop if groups are unusable
tab <- table(hsc$numbat_malignancy)

if (any(tab[c("CNV_aberrant", "CNV_neutral")] < 5, na.rm = TRUE)) {
  warning(
    "One Numbat group has fewer than 5 cells. Differential marker testing may be unstable."
  )
}

# ============================================================
# DSB MATRIX
# ============================================================

DefaultAssay(hsc) <- "CITE_DSB"

dsb <- GetAssayData(
  hsc,
  assay = "CITE_DSB",
  layer = "data"
)

keep_markers <- rownames(dsb)[
  !grepl(
    "IgG|IgM|isotype|Hash|HTO|TotalSeq|control|Fc",
    rownames(dsb),
    ignore.case = TRUE
  )
]

dsb <- dsb[keep_markers, , drop = FALSE]

# ============================================================
# DIFFERENTIAL DSB MARKERS: NUMBAT ABERRANT VS NEUTRAL
# ============================================================

group_vec <- hsc$numbat_malignancy

marker_stats <- purrr::map_dfr(
  rownames(dsb),
  function(marker) {
    
    vals <- as.numeric(dsb[marker, ])
    
    cancer_vals <- vals[group_vec == "CNV_aberrant"]
    normal_vals <- vals[group_vec == "CNV_neutral"]
    
    tibble(
      marker = marker,
      
      n_cancer = length(cancer_vals),
      n_normal = length(normal_vals),
      
      mean_cancer = mean(cancer_vals, na.rm = TRUE),
      mean_normal = mean(normal_vals, na.rm = TRUE),
      
      median_cancer = median(cancer_vals, na.rm = TRUE),
      median_normal = median(normal_vals, na.rm = TRUE),
      
      pct_cancer_positive = mean(cancer_vals > 0, na.rm = TRUE) * 100,
      pct_normal_positive = mean(normal_vals > 0, na.rm = TRUE) * 100,
      
      delta_mean = mean_cancer - mean_normal,
      delta_median = median_cancer - median_normal,
      delta_pct_positive = pct_cancer_positive - pct_normal_positive,
      
      log2FC = log2(
        (mean(cancer_vals, na.rm = TRUE) + 1e-6) /
          (mean(normal_vals, na.rm = TRUE) + 1e-6)
      ),
      
      p_value = tryCatch(
        wilcox.test(cancer_vals, normal_vals)$p.value,
        error = function(e) NA_real_
      )
    )
  }
) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

write.csv(
  marker_stats,
  file.path(out_dir, "HSC_MPP_DSB_numbat_aberrant_vs_neutral_markers.csv"),
  row.names = FALSE
)

# ============================================================
# TOP MARKERS
# ============================================================

top_up <- marker_stats %>%
  dplyr::filter(is.finite(log2FC), !is.na(FDR)) %>%
  dplyr::arrange(FDR, dplyr::desc(log2FC)) %>%
  dplyr::slice_head(n = 10)

top_down <- marker_stats %>%
  dplyr::filter(is.finite(log2FC), !is.na(FDR)) %>%
  dplyr::arrange(FDR, log2FC) %>%
  dplyr::slice_head(n = 10)

top_markers <- unique(c(top_up$marker, top_down$marker))

# ============================================================
# VOLCANO
# ============================================================

volcano_df <- marker_stats %>%
  mutate(
    significant = FDR < 0.05,
    neglog10FDR = -log10(FDR)
  )

p_volcano <- ggplot(
  volcano_df,
  aes(x = log2FC, y = neglog10FDR)
) +
  geom_point(aes(color = significant), alpha = 0.8) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "red")
  ) +
  ggrepel::geom_text_repel(
    data = volcano_df %>%
      dplyr::filter(marker %in% top_markers),
    aes(label = marker),
    max.overlaps = 50,
    size = 3
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("HSC-MPP: Numbat CNV-aberrant vs CNV-neutral DSB markers\n", sample_to_use),
    x = "log2FC CNV_aberrant / CNV_neutral",
    y = "-log10(FDR)"
  )

ggsave(
  file.path(out_dir, "01_HSC_MPP_DSB_numbat_volcano.pdf"),
  p_volcano,
  width = 8,
  height = 7
)

# ============================================================
# UMAP
# ============================================================

reduction_use <- if ("umap_projected" %in% Reductions(seu_sub)) {
  "umap_projected"
} else {
  "umap"
}

numbat_cols <- c(
  "CNV_aberrant" = "#E41A1C",
  "CNV_neutral" = "#377EB8",
  "Not_called" = "grey80"
)

p_umap_all <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "numbat_malignancy",
  cols = numbat_cols
) +
  theme_bw(base_size = 12) +
  ggtitle(paste0("Numbat malignancy: ", sample_to_use))

ggsave(
  file.path(out_dir, "02_all_cells_numbat_malignancy_umap.pdf"),
  p_umap_all,
  width = 6,
  height = 5
)

p_umap_hsc <- DimPlot(
  hsc,
  reduction = reduction_use,
  group.by = "numbat_malignancy",
  cols = numbat_cols
) +
  theme_bw(base_size = 12) +
  ggtitle(paste0("HSC-MPP Numbat states: ", sample_to_use))

ggsave(
  file.path(out_dir, "03_HSC_MPP_numbat_malignancy_umap.pdf"),
  p_umap_hsc,
  width = 6,
  height = 5
)

# ============================================================
# FEATURE PLOTS
# ============================================================

pdf(
  file.path(out_dir, "04_HSC_MPP_top_DSB_featureplots_numbat.pdf"),
  width = 12,
  height = 6
)

for (mk in top_markers) {
  
  if (!mk %in% rownames(hsc[["CITE_DSB"]])) next
  
  p1 <- FeaturePlot(
    hsc,
    features = mk,
    reduction = reduction_use,
    order = TRUE
  ) +
    ggtitle(mk)
  
  p2 <- DimPlot(
    hsc,
    reduction = reduction_use,
    group.by = "numbat_malignancy",
    cols = numbat_cols
  ) +
    ggtitle("Numbat")
  
  print(p1 + p2)
}

dev.off()

# ============================================================
# VIOLIN PLOTS
# ============================================================

pdf(
  file.path(out_dir, "05_HSC_MPP_top_DSB_violins_numbat.pdf"),
  width = 10,
  height = 5
)

for (mk in top_markers) {
  
  if (!mk %in% rownames(hsc[["CITE_DSB"]])) next
  
  p <- VlnPlot(
    hsc,
    features = mk,
    group.by = "numbat_malignancy",
    assay = "CITE_DSB",
    pt.size = 0
  ) +
    theme_bw(base_size = 12) +
    ggtitle(mk)
  
  print(p)
}

dev.off()

# ============================================================
# DOTPLOT
# ============================================================

p_dot <- DotPlot(
  hsc,
  features = rev(top_markers),
  group.by = "numbat_malignancy",
  assay = "CITE_DSB"
) +
  RotatedAxis() +
  theme_bw(base_size = 11)

ggsave(
  file.path(out_dir, "06_HSC_MPP_DSB_dotplot_numbat.pdf"),
  p_dot,
  width = 10,
  height = 7
)

# ============================================================
# OPTIONAL: DOTPLOT BY NUMBAT CLONE
# ============================================================

if (length(unique(na.omit(hsc$numbat_clone))) > 1) {
  
  p_clone_dot <- DotPlot(
    hsc,
    features = rev(top_markers),
    group.by = "numbat_clone",
    assay = "CITE_DSB"
  ) +
    RotatedAxis() +
    theme_bw(base_size = 11) +
    ggtitle("Top DSB markers by Numbat clone")
  
  ggsave(
    file.path(out_dir, "07_HSC_MPP_DSB_dotplot_by_numbat_clone.pdf"),
    p_clone_dot,
    width = 10,
    height = 7
  )
}

# ============================================================
# SAVE
# ============================================================

saveRDS(
  hsc,
  file.path(out_dir, "HSC_MPP_numbat_annotated.rds")
)

saveRDS(
  seu_sub,
  file.path(out_dir, "sample_numbat_annotated.rds")
)

cat("\nDone.\n")
cat("Output:", out_dir, "\n")