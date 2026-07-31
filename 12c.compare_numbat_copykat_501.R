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

sample_to_use <- "HBDN501-AML-KMT2A"
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
  select(cell, numbat_clone, numbat_malignancy, everything())
cat("Numbat cells:", nrow(numbat_meta), "\n")
print(table(numbat_meta$numbat_malignancy, useNA = "ifany"))

# ============================================================
# LOAD COPYKAT
# ============================================================

copykat_meta <- fread(copykat_file) %>%
  mutate(
    copykat_malignancy = case_when(
      copykat_call == "aneuploid" ~ "CNV_aberrant",
      copykat_call == "diploid" ~ "CNV_neutral",
      is.na(copykat_call) ~ "Not_called",
      TRUE ~ as.character(copykat_call)
    )
  ) %>%
  select(cell, copykat_call, copykat_malignancy, everything())

cat("CopyKAT cells:", nrow(copykat_meta), "\n")
print(table(copykat_meta$copykat_malignancy, useNA = "ifany"))

# ============================================================
# MERGE CALLS ONTO SEURAT METADATA
# ============================================================

meta <- seu_sub@meta.data %>%
  as_tibble(rownames = "cell") %>%
  left_join(
    numbat_meta %>% select(cell, numbat_clone, numbat_malignancy),
    by = "cell"
  ) %>%
  left_join(
    copykat_meta %>% select(cell, copykat_call, copykat_malignancy),
    by = "cell"
  ) %>%
  mutate(
    numbat_malignancy = replace_na(numbat_malignancy, "Not_called"),
    copykat_malignancy = replace_na(copykat_malignancy, "Not_called"),
    
    both_called = numbat_malignancy %in% c("CNV_aberrant", "CNV_neutral") &
      copykat_malignancy %in% c("CNV_aberrant", "CNV_neutral"),
    
    cnv_overlap_class = case_when(
      !both_called ~ "Not_called_by_one_or_both",
      numbat_malignancy == "CNV_aberrant" &
        copykat_malignancy == "CNV_aberrant" ~ "Both_CNV_aberrant",
      numbat_malignancy == "CNV_neutral" &
        copykat_malignancy == "CNV_neutral" ~ "Both_CNV_neutral",
      numbat_malignancy == "CNV_aberrant" &
        copykat_malignancy == "CNV_neutral" ~ "Numbat_aberrant_CopyKAT_neutral",
      numbat_malignancy == "CNV_neutral" &
        copykat_malignancy == "CNV_aberrant" ~ "Numbat_neutral_CopyKAT_aberrant",
      TRUE ~ "Other"
    ),
    
    cnv_binary_agreement = case_when(
      !both_called ~ "Not_called_by_one_or_both",
      numbat_malignancy == copykat_malignancy ~ "Agree",
      TRUE ~ "Disagree"
    )
  )

# Add back to Seurat
seu_sub$numbat_malignancy <- meta$numbat_malignancy[match(colnames(seu_sub), meta$cell)]
seu_sub$numbat_clone <- meta$numbat_clone[match(colnames(seu_sub), meta$cell)]
seu_sub$copykat_malignancy <- meta$copykat_malignancy[match(colnames(seu_sub), meta$cell)]
seu_sub$copykat_call <- meta$copykat_call[match(colnames(seu_sub), meta$cell)]
seu_sub$cnv_overlap_class <- meta$cnv_overlap_class[match(colnames(seu_sub), meta$cell)]
seu_sub$cnv_binary_agreement <- meta$cnv_binary_agreement[match(colnames(seu_sub), meta$cell)]

# ============================================================
# CORE OVERLAP REPORT
# ============================================================
# ============================================================
# CORE OVERLAP REPORT — robust to missing classes
# ============================================================

celltype_col <- "predicted_CellType_Broad"
# ============================================================
# BINARY CNV ABERRANT COMPARISON
# CNV_aberrant vs everything else
# ============================================================

binary_levels <- c("CNV_aberrant", "Not_aberrant_or_not_called")

meta <- meta %>%
  mutate(
    numbat_binary = ifelse(
      numbat_malignancy == "CNV_aberrant",
      "CNV_aberrant",
      "Not_aberrant_or_not_called"
    ),
    copykat_binary = ifelse(
      copykat_malignancy == "CNV_aberrant",
      "CNV_aberrant",
      "Not_aberrant_or_not_called"
    ),
    numbat_binary = factor(numbat_binary, levels = binary_levels),
    copykat_binary = factor(copykat_binary, levels = binary_levels),
    
    cnv_binary_overlap = case_when(
      numbat_binary == "CNV_aberrant" &
        copykat_binary == "CNV_aberrant" ~
        "Both_CNV_aberrant",
      
      numbat_binary == "CNV_aberrant" &
        copykat_binary == "Not_aberrant_or_not_called" ~
        "Numbat_only_CNV_aberrant",
      
      numbat_binary == "Not_aberrant_or_not_called" &
        copykat_binary == "CNV_aberrant" ~
        "CopyKAT_only_CNV_aberrant",
      
      TRUE ~ "Neither_CNV_aberrant"
    ),
    cnv_binary_overlap = factor(
      cnv_binary_overlap,
      levels = c(
        "Both_CNV_aberrant",
        "Numbat_only_CNV_aberrant",
        "CopyKAT_only_CNV_aberrant",
        "Neither_CNV_aberrant"
      )
    )
  )

seu_sub$numbat_binary <- meta$numbat_binary[match(colnames(seu_sub), meta$cell)]
seu_sub$copykat_binary <- meta$copykat_binary[match(colnames(seu_sub), meta$cell)]
seu_sub$cnv_binary_overlap <- meta$cnv_binary_overlap[match(colnames(seu_sub), meta$cell)]

# ============================================================
# TABLES
# ============================================================

binary_overlap_table <- meta %>%
  count(numbat_binary, copykat_binary, name = "n_cells") %>%
  tidyr::complete(
    numbat_binary = binary_levels,
    copykat_binary = binary_levels,
    fill = list(n_cells = 0)
  ) %>%
  group_by(numbat_binary) %>%
  mutate(
    total_numbat = sum(n_cells),
    percent_within_numbat = ifelse(
      total_numbat > 0,
      100 * n_cells / total_numbat,
      0
    )
  ) %>%
  ungroup()

binary_overlap_summary <- meta %>%
  count(cnv_binary_overlap, name = "n_cells") %>%
  tidyr::complete(
    cnv_binary_overlap = levels(meta$cnv_binary_overlap),
    fill = list(n_cells = 0)
  ) %>%
  mutate(percent = 100 * n_cells / sum(n_cells))

write.csv(
  binary_overlap_table,
  file.path(out_dir, "01_binary_numbat_copykat_overlap_table.csv"),
  row.names = FALSE
)

write.csv(
  binary_overlap_summary,
  file.path(out_dir, "02_binary_overlap_summary.csv"),
  row.names = FALSE
)

# ============================================================
# CELL TYPE TABLE
# ============================================================

celltype_col <- "predicted_CellType_Broad"

celltype_by_binary_overlap <- meta %>%
  count(cnv_binary_overlap, .data[[celltype_col]], name = "n_cells") %>%
  tidyr::complete(
    cnv_binary_overlap = levels(meta$cnv_binary_overlap),
    !!rlang::sym(celltype_col) := sort(unique(meta[[celltype_col]])),
    fill = list(n_cells = 0)
  ) %>%
  group_by(cnv_binary_overlap) %>%
  mutate(
    total_overlap_class = sum(n_cells),
    percent_within_overlap_class = ifelse(
      total_overlap_class > 0,
      100 * n_cells / total_overlap_class,
      0
    )
  ) %>%
  ungroup()

write.csv(
  celltype_by_binary_overlap,
  file.path(out_dir, "03_celltype_by_binary_overlap.csv"),
  row.names = FALSE
)

# ============================================================
# HEATMAP
# ============================================================

p_binary_overlap <- ggplot(
  binary_overlap_table,
  aes(
    x = copykat_binary,
    y = numbat_binary,
    fill = n_cells
  )
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(
      label = n_cells,
      color = n_cells > median(n_cells)
    ),
    size = 5,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"),
    guide = "none"
  ) +
  scale_fill_gradient(
    low = "#DCEAF7",
    high = "#08519C"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = paste0("Binary CNV-aberrant overlap: ", sample_to_use),
    x = "CopyKAT",
    y = "Numbat",
    fill = "Cells"
  )

ggsave(
  file.path(out_dir, "04_binary_overlap_heatmap.pdf"),
  p_binary_overlap,
  width = 7,
  height = 5
)

# ============================================================
# BARPLOT
# ============================================================

p_binary_bar <- binary_overlap_summary %>%
  ggplot(
    aes(
      x = cnv_binary_overlap,
      y = n_cells,
      fill = cnv_binary_overlap
    )
  ) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Binary CNV-aberrant overlap classes: ", sample_to_use),
    x = NULL,
    y = "Cells"
  )

ggsave(
  file.path(out_dir, "05_binary_overlap_classes.pdf"),
  p_binary_bar,
  width = 8,
  height = 5
)

# ============================================================
# CELL TYPE COMPOSITION
# ============================================================
celltype_order <- rev(c(
  "HSC MPP", "LMPP", "MEP", "Megakaryocyte Precursor",
  "GMP", "Early GMP", "Late GMP", "Cycling Progenitor",
  "EoBasoMast Precursor", "Pro-Monocyte", "Monocyte",
  "cDC", "pDC", "Early Lymphoid", "Pro-B", "Pre-B",
  "B", "Naive T", "CD4 Memory T", "CD8 Memory T",
  "NK", "Plasma Cell", "Early Erythroid", "Late Erythroid",
  "Stromal"
))

# ------------------------------------------------------------
# preserve order + add unexpected cell types
# ------------------------------------------------------------

observed_celltypes <- unique(
  celltype_by_binary_overlap[[celltype_col]]
)

celltype_levels <- c(
  celltype_order[celltype_order %in% observed_celltypes],
  setdiff(observed_celltypes, celltype_order)
)

# ------------------------------------------------------------
# ensure all plotted cell types have colours
# ------------------------------------------------------------

missing_cols <- setdiff(
  celltype_levels,
  names(celltype_cols)
)

if (length(missing_cols) > 0) {
  
  warning(
    "Missing colours for: ",
    paste(missing_cols, collapse = ", ")
  )
  
  extra_cols <- setNames(
    rep("grey70", length(missing_cols)),
    missing_cols
  )
  
  celltype_cols <- c(celltype_cols, extra_cols)
}

celltype_cols_use <- celltype_cols[celltype_levels]

# ------------------------------------------------------------
# plotting dataframe
# ------------------------------------------------------------

celltype_by_binary_overlap_plot <- celltype_by_binary_overlap %>%
  mutate(
    celltype_plot = factor(
      .data[[celltype_col]],
      levels = celltype_levels
    ),
    cnv_binary_overlap = factor(
      cnv_binary_overlap,
      levels = c(
        "Both_CNV_aberrant",
        "Numbat_only_CNV_aberrant",
        "CopyKAT_only_CNV_aberrant",
        "Neither_CNV_aberrant"
      )
    )
  )

# ------------------------------------------------------------
# plot
# ------------------------------------------------------------

p_celltype_binary <- ggplot(
  celltype_by_binary_overlap_plot,
  aes(
    x = cnv_binary_overlap,
    y = percent_within_overlap_class,
    fill = celltype_plot
  )
) +
  geom_col(
    color = "white",
    linewidth = 0.15
  ) +
  
  scale_fill_manual(
    values = celltype_cols_use,
    drop = FALSE,
    na.value = "grey80"
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm")
  ) +
  
  labs(
    title = paste0(
      "Cell type composition by binary CNV-aberrant overlap: ",
      sample_to_use
    ),
    x = NULL,
    y = "% within overlap class",
    fill = "Cell type"
  )

ggsave(
  file.path(out_dir, "06_celltype_by_binary_overlap.pdf"),
  p_celltype_binary,
  width = 11,
  height = 7
)
# ============================================================
# UMAP
# ============================================================

p_binary_umap <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "cnv_binary_overlap"
) +
  ggtitle(paste0("Binary CNV-aberrant overlap: ", sample_to_use))

ggsave(
  file.path(out_dir, "07_binary_overlap_umap.pdf"),
  p_binary_umap,
  width = 7,
  height = 6
)

cat("\nBinary overlap summary:\n")
print(binary_overlap_summary) 

DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "cnv_binary_overlap",
  cols = c(
    "Both_CNV_aberrant" = "red",
    "Numbat_only_CNV_aberrant" = "orange",
    "CopyKAT_only_CNV_aberrant" = "blue",
    "Neither_CNV_aberrant" = "grey80"
  )
)

DimPlot(
  seu_sub,
  reduction = reduction_use,
  split.by = "cnv_binary_overlap",
  group.by = celltype_col
)  

seu_sub$numbat_clone <- factor(
  seu_sub$numbat_clone,
  levels = sort(unique(seu_sub$numbat_clone))
)

# nice colours
clone_cols <- setNames(
  c(
    "#D73027",
    "#4575B4",
    "#1A9850",
    "#984EA3",
    "#FF7F00",
    "#A65628",
    "#F781BF",
    "#999999"
  )[seq_len(length(levels(seu_sub$numbat_clone)))],
  levels(seu_sub$numbat_clone)
)

# ------------------------------------------------------------
# clone-only UMAP
# ------------------------------------------------------------

p_clone_umap <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "numbat_clone",
  cols = clone_cols,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle(
    paste0(
      "Numbat clones: ",
      sample_to_use
    )
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "16_numbat_clone_umap.pdf"),
  p_clone_umap,
  width = 8,
  height = 7
)

# ------------------------------------------------------------
# split by overlap class
# very useful biologically
# ------------------------------------------------------------

p_clone_split <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "numbat_clone",
  split.by = "cnv_binary_overlap",
  cols = clone_cols,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  ncol = 2
) +
  ggtitle(
    paste0(
      "Numbat clones split by binary overlap: ",
      sample_to_use
    )
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(out_dir, "17_numbat_clone_split_by_overlap.pdf"),
  p_clone_split,
  width = 8,
  height = 10
)

# ------------------------------------------------------------
# optional:
# highlight Numbat-only cells
# ------------------------------------------------------------

numbat_only_cells <- meta$cell[
  meta$cnv_binary_overlap ==
    "Numbat_only_CNV_aberrant"
]

p_numbat_only <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  cells.highlight = numbat_only_cells,
  cols.highlight = "#FF7F00",
  cols = "grey85",
  sizes.highlight = 0.6,
  raster = FALSE
) +
  ggtitle(
    paste0(
      "Numbat-only CNV-aberrant cells: ",
      sample_to_use
    )
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "18_numbat_only_cells_umap.pdf"),
  p_numbat_only,
  width = 7,
  height = 6
)

# ------------------------------------------------------------
# clone x overlap contingency table
# useful for manuscript figures
# ------------------------------------------------------------

clone_overlap_table <- meta %>%
  count(
    numbat_clone,
    cnv_binary_overlap,
    name = "n_cells"
  ) %>%
  group_by(numbat_clone) %>%
  mutate(
    percent_within_clone =
      100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

write.csv(
  clone_overlap_table,
  file.path(out_dir, "19_clone_by_overlap_table.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# stacked barplot
# ------------------------------------------------------------

p_clone_bar <- ggplot(
  clone_overlap_table,
  aes(
    x = factor(numbat_clone),
    y = percent_within_clone,
    fill = cnv_binary_overlap
  )
) +
  geom_col(color = "white") +
  theme_bw(base_size = 12) +
  labs(
    title = paste0(
      "Binary overlap composition within Numbat clones: ",
      sample_to_use
    ),
    x = "Numbat clone",
    y = "% within clone",
    fill = "Binary overlap"
  )

ggsave(
  file.path(out_dir, "20_clone_overlap_barplot.pdf"),
  p_clone_bar,
  width = 8,
  height = 5
)


