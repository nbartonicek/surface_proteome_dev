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

copykat_file <- file.path(
  annotation_dir,
  "copykat_annotated",
  sample_to_use,
  "copykat_prediction_with_metadata.csv"
)

out_dir <- file.path(
  projection_dir,
  paste0("HSC_MPP_copykat_DSB_", sample_to_use)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD COPYKAT
# ============================================================

copykat_meta <- fread(copykat_file)
copykat_meta <- copykat_meta %>%
  dplyr::mutate(
    copykat_malignancy = dplyr::case_when(
      copykat_call == "aneuploid" ~ "CNV_aberrant",
      copykat_call == "diploid" ~ "CNV_neutral",
      is.na(copykat_call) ~ "Not_called",
      TRUE ~ as.character(copykat_call)
    )
  )

head(copykat_meta)

# expected columns:
# cell
# copykat_call
# copykat_malignancy

# ============================================================
# SUBSET SAMPLE
# ============================================================

seu_sub <- subset(
  seu,
  subset = sampleID == sample_to_use
)

cat("Cells:", ncol(seu_sub), "\n")

# ============================================================
# ADD COPYKAT METADATA
# ============================================================

seu_sub$copykat_call <- NA_character_
seu_sub$copykat_malignancy <- "Not_called"

common_cells <- intersect(
  colnames(seu_sub),
  copykat_meta$cell
)

seu_sub$copykat_call[common_cells] <-
  copykat_meta$copykat_call[
    match(common_cells, copykat_meta$cell)
  ]

seu_sub$copykat_malignancy[common_cells] <-
  copykat_meta$copykat_malignancy[
    match(common_cells, copykat_meta$cell)
  ]

table(seu_sub$copykat_malignancy)

# ============================================================
# LOOK AT HSC-MPP ONLY
# ============================================================

hsc <- subset(
  seu_sub,
  subset = predicted_CellType_Broad == "HSC MPP"
)

cat("HSC-MPP cells:", ncol(hsc), "\n")

table(hsc$copykat_malignancy)

# ============================================================
# REMOVE NOT CALLED
# ============================================================

hsc <- subset(
  hsc,
  subset = copykat_malignancy %in% c(
    "CNV_aberrant",
    "CNV_neutral"
  )
)

table(hsc$copykat_malignancy)

# ============================================================
# DSB MATRIX
# ============================================================

DefaultAssay(hsc) <- "CITE_DSB"

dsb <- GetAssayData(
  hsc,
  assay = "CITE_DSB",
  layer = "data"
)

# remove obvious background controls
keep_markers <- rownames(dsb)[
  !grepl(
    "IgG|IgM|isotype|Hash|HTO|TotalSeq|control|Fc",
    rownames(dsb),
    ignore.case = TRUE
  )
]

dsb <- dsb[keep_markers, , drop = FALSE]

# ============================================================
# DIFFERENTIAL DSB SIGNAL
# ============================================================

group_vec <- hsc$copykat_malignancy

marker_stats <- map_dfr(
  rownames(dsb),
  function(marker) {
    
    vals <- as.numeric(dsb[marker, ])
    
    cancer_vals <- vals[group_vec == "CNV_aberrant"]
    normal_vals <- vals[group_vec == "CNV_neutral"]
    
    tibble(
      marker = marker,
      
      mean_cancer = mean(cancer_vals, na.rm = TRUE),
      mean_normal = mean(normal_vals, na.rm = TRUE),
      
      median_cancer = median(cancer_vals, na.rm = TRUE),
      median_normal = median(normal_vals, na.rm = TRUE),
      
      pct_cancer_positive =
        mean(cancer_vals > 0, na.rm = TRUE) * 100,
      
      pct_normal_positive =
        mean(normal_vals > 0, na.rm = TRUE) * 100,
      
      log2FC =
        log2(
          (mean(cancer_vals, na.rm = TRUE) + 1e-6) /
            (mean(normal_vals, na.rm = TRUE) + 1e-6)
        ),
      
      delta_pct_positive =
        pct_cancer_positive -
        pct_normal_positive,
      
      p_value = tryCatch(
        wilcox.test(cancer_vals, normal_vals)$p.value,
        error = function(e) NA_real_
      )
    )
  }
)

marker_stats <- marker_stats %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(FDR)

write.csv(
  marker_stats,
  file.path(out_dir, "HSC_MPP_DSB_copykat_markers.csv"),
  row.names = FALSE
)

head(marker_stats)

# ============================================================
# TOP MARKERS
# ============================================================

top_up <- marker_stats %>%
  dplyr::filter(
    is.finite(log2FC),
    !is.na(FDR)
  ) %>%
  arrange(FDR, desc(log2FC)) %>%
  slice_head(n = 10)

top_down <- marker_stats %>%
  dplyr::filter(
    is.finite(log2FC),
    !is.na(FDR)
  ) %>%
  arrange(FDR, log2FC) %>%
  slice_head(n = 10)

top_markers <- unique(c(
  top_up$marker,
  top_down$marker
))

# ============================================================
# VOLCANO PLOT
# ============================================================

volcano_df <- marker_stats %>%
  mutate(
    significant = FDR < 0.05
  )

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FC,
    y = -log10(FDR)
  )
) +
  geom_point(
    aes(color = significant),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE" = "red"
    )
  ) +
  ggrepel::geom_text_repel(
    data = volcano_df %>%
      filter(marker %in% top_markers),
    aes(label = marker),
    max.overlaps = 50,
    size = 3
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = paste0(
      "HSC-MPP: malignant vs normal DSB markers\n",
      sample_to_use
    ),
    x = "log2FC (cancer / normal)",
    y = "-log10(FDR)"
  )

ggsave(
  file.path(out_dir, "01_HSC_MPP_DSB_volcano.pdf"),
  p_volcano,
  width = 8,
  height = 7
)

# ============================================================
# FEATURE PLOTS
# ============================================================

reduction_use <- if (
  "umap_projected" %in% Reductions(seu_sub)
) {
  "umap_projected"
} else {
  "umap"
}

pdf(
  file.path(out_dir, "02_HSC_MPP_top_DSB_featureplots.pdf"),
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
    group.by = "copykat_malignancy",
    cols = c(
      "CNV_aberrant" = "#E41A1C",
      "CNV_neutral" = "#377EB8"
    )
  ) +
    ggtitle("CopyKAT")
  
  print(p1 + p2)
}

dev.off()

# ============================================================
# VIOLIN PLOTS
# ============================================================

pdf(
  file.path(out_dir, "03_HSC_MPP_top_DSB_violins.pdf"),
  width = 10,
  height = 5
)

for (mk in top_markers) {
  
  if (!mk %in% rownames(hsc[["CITE_DSB"]])) next
  
  p <- VlnPlot(
    hsc,
    features = mk,
    group.by = "copykat_malignancy",
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
  group.by = "copykat_malignancy",
  assay = "CITE_DSB"
) +
  RotatedAxis() +
  theme_bw(base_size = 11)

ggsave(
  file.path(out_dir, "04_HSC_MPP_DSB_dotplot.pdf"),
  p_dot,
  width = 10,
  height = 7
)

# ============================================================
# UMAP OF HSC-MPP
# ============================================================

p_umap <- DimPlot(
  hsc,
  reduction = reduction_use,
  group.by = "copykat_malignancy",
  cols = c(
    "CNV_aberrant" = "#E41A1C",
    "CNV_neutral" = "#377EB8"
  )
) +
  theme_bw(base_size = 12) +
  ggtitle(
    paste0(
      "HSC-MPP CopyKAT states: ",
      sample_to_use
    )
  )

ggsave(
  file.path(out_dir, "05_HSC_MPP_copykat_umap.pdf"),
  p_umap,
  width = 6,
  height = 5
)

# ============================================================
# SAVE OBJECT
# ============================================================

saveRDS(
  hsc,
  file.path(out_dir, "HSC_MPP_copykat_annotated.rds")
)

cat("\nDone.\n")
cat("Output:", out_dir, "\n")