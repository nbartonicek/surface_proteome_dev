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
reduction_use <- "umap_projected"
 
sample_to_use <- "HBDN206-MNpCT"
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
  paste0("numbat_copykat_comparison_noDoublet_", sample_to_use)
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
seu_sub <- subset(seu_sub, subset = scDblFinder.class == "singlet")
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

numbat_keep <- numbat_meta %>%
  mutate(cell = as.character(cell)) %>%
  select(any_of(c("cell", "numbat_clone", "numbat_malignancy"))) %>%
  distinct(cell, .keep_all = TRUE)

copykat_keep <- copykat_meta %>%
  mutate(cell = as.character(cell)) %>%
  select(any_of(c("cell", "copykat_call", "copykat_malignancy"))) %>%
  distinct(cell, .keep_all = TRUE)

# ============================================================
# MERGE NUMBAT + COPYKAT CALLS
# Binary comparison: CNV_aberrant vs everything else
# ============================================================

# Start clean from Seurat metadata
meta <- seu_sub@meta.data %>%
  as_tibble(rownames = "cell") %>%
  dplyr::select(
    -matches("numbat", ignore.case = TRUE),
    -matches("copykat", ignore.case = TRUE),
    -matches("^cnv_", ignore.case = TRUE)
  ) %>%
  mutate(cell = as.character(cell)) %>%
  left_join(numbat_keep, by = "cell") %>%
  left_join(copykat_keep, by = "cell")

# Add missing columns if needed
needed_cols <- c(
  "numbat_clone",
  "numbat_malignancy",
  "copykat_call",
  "copykat_malignancy"
)

for (cc in needed_cols) {
  if (!cc %in% colnames(meta)) {
    meta[[cc]] <- NA_character_
  }
}

# Binary CNV-aberrant classification
meta <- meta %>%
  mutate(
    numbat_clone = as.character(numbat_clone),
    
    numbat_malignancy = replace_na(
      as.character(numbat_malignancy),
      "Not_called"
    ),
    
    copykat_call = replace_na(
      as.character(copykat_call),
      "Not_called"
    ),
    
    copykat_malignancy = replace_na(
      as.character(copykat_malignancy),
      "Not_called"
    ),
    
    numbat_binary = if_else(
      numbat_malignancy == "CNV_aberrant",
      "CNV_aberrant",
      "Not_CNV_aberrant"
    ),
    
    copykat_binary = if_else(
      copykat_malignancy == "CNV_aberrant",
      "CNV_aberrant",
      "Not_CNV_aberrant"
    ),
    
    cnv_binary_overlap = case_when(
      numbat_binary == "CNV_aberrant" &
        copykat_binary == "CNV_aberrant" ~
        "Both_CNV_aberrant",
      
      numbat_binary == "CNV_aberrant" &
        copykat_binary != "CNV_aberrant" ~
        "Numbat_only_CNV_aberrant",
      
      numbat_binary != "CNV_aberrant" &
        copykat_binary == "CNV_aberrant" ~
        "CopyKAT_only_CNV_aberrant",
      
      TRUE ~ "Neither_CNV_aberrant"
    ),
    
    cnv_binary_agreement = case_when(
      cnv_binary_overlap %in% c(
        "Both_CNV_aberrant",
        "Neither_CNV_aberrant"
      ) ~ "Agree",
      
      cnv_binary_overlap %in% c(
        "Numbat_only_CNV_aberrant",
        "CopyKAT_only_CNV_aberrant"
      ) ~ "Disagree",
      
      TRUE ~ "Other"
    ),
    
    numbat_binary = factor(
      numbat_binary,
      levels = c("CNV_aberrant", "Not_CNV_aberrant")
    ),
    
    copykat_binary = factor(
      copykat_binary,
      levels = c("CNV_aberrant", "Not_CNV_aberrant")
    ),
    
    cnv_binary_overlap = factor(
      cnv_binary_overlap,
      levels = c(
        "Both_CNV_aberrant",
        "Numbat_only_CNV_aberrant",
        "CopyKAT_only_CNV_aberrant",
        "Neither_CNV_aberrant"
      )
    ),
    
    cnv_binary_agreement = factor(
      cnv_binary_agreement,
      levels = c("Agree", "Disagree", "Other")
    )
  )

# ============================================================
# ADD BACK TO SEURAT
# ============================================================

add_meta_cols <- c(
  "numbat_clone",
  "numbat_malignancy",
  "copykat_call",
  "copykat_malignancy",
  "numbat_binary",
  "copykat_binary",
  "cnv_binary_overlap",
  "cnv_binary_agreement"
)

for (cc in add_meta_cols) {
  seu_sub[[cc]] <- meta[[cc]][match(colnames(seu_sub), meta$cell)]
}

# ============================================================
# QUICK CHECKS
# ============================================================

cat("\nNumbat malignancy:\n")
print(table(meta$numbat_malignancy, useNA = "ifany"))

cat("\nCopyKAT malignancy:\n")
print(table(meta$copykat_malignancy, useNA = "ifany"))

cat("\nBinary overlap:\n")
print(table(meta$cnv_binary_overlap, useNA = "ifany"))

cat("\nBinary agreement:\n")
print(table(meta$cnv_binary_agreement, useNA = "ifany"))

cat("\nNumbat clone x binary overlap:\n")
print(table(meta$numbat_clone, meta$cnv_binary_overlap, useNA = "ifany"))
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

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

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
  width = 14,
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

##########
# ============================================================
# INVESTIGATE NUMBAT CLONE 1 — especially lymphoid cells
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Seurat)
  library(patchwork)
})

clone_to_check <- "1"

clone_out_dir <- file.path(out_dir, paste0("clone_", clone_to_check, "_inspection"))
dir.create(clone_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Make sure clone is character
# ------------------------------------------------------------

meta <- meta %>%
  mutate(
    numbat_clone = as.character(numbat_clone),
    clone_focus = case_when(
      numbat_clone == clone_to_check ~ paste0("Clone_", clone_to_check),
      is.na(numbat_clone) ~ "No_clone",
      TRUE ~ "Other_clones"
    )
  )

seu_sub$numbat_clone <- as.character(seu_sub$numbat_clone)
seu_sub$clone_focus <- meta$clone_focus[match(colnames(seu_sub), meta$cell)]

# ------------------------------------------------------------
# Define broad lymphoid labels
# ------------------------------------------------------------

lymphoid_celltypes <- c(
  "Early Lymphoid",
  "Pro-B", "Pre-B", "B",
  "Naive T", "CD4 Memory T", "CD8 Memory T",
  "NK",
  "Plasma Cell"
)

meta <- meta %>%
  mutate(
    is_clone_focus = numbat_clone == clone_to_check,
    is_lymphoid = .data[[celltype_col]] %in% lymphoid_celltypes,
    clone1_lymphoid_status = case_when(
      is_clone_focus & is_lymphoid ~ paste0("Clone_", clone_to_check, "_lymphoid"),
      is_clone_focus & !is_lymphoid ~ paste0("Clone_", clone_to_check, "_non_lymphoid"),
      !is_clone_focus & is_lymphoid ~ "Other_lymphoid",
      TRUE ~ "Other_non_lymphoid"
    )
  )

seu_sub$clone1_lymphoid_status <- meta$clone1_lymphoid_status[
  match(colnames(seu_sub), meta$cell)
]

# ============================================================
# 1. BASIC TABLES
# ============================================================

clone_celltype_table <- meta %>%
  count(numbat_clone, .data[[celltype_col]], name = "n_cells") %>%
  group_by(numbat_clone) %>%
  mutate(percent_within_clone = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(numbat_clone, desc(n_cells))

clone_overlap_table <- meta %>%
  count(numbat_clone, cnv_binary_overlap, name = "n_cells") %>%
  group_by(numbat_clone) %>%
  mutate(percent_within_clone = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(numbat_clone, desc(n_cells))

clone1_summary <- meta %>%
  filter(numbat_clone == clone_to_check) %>%
  count(.data[[celltype_col]], cnv_binary_overlap, name = "n_cells") %>%
  group_by(.data[[celltype_col]]) %>%
  mutate(percent_within_celltype = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(desc(n_cells))

write.csv(
  clone_celltype_table,
  file.path(clone_out_dir, "01_numbat_clone_by_celltype.csv"),
  row.names = FALSE
)

write.csv(
  clone_overlap_table,
  file.path(clone_out_dir, "02_numbat_clone_by_binary_overlap.csv"),
  row.names = FALSE
)

write.csv(
  clone1_summary,
  file.path(clone_out_dir, paste0("03_clone_", clone_to_check, "_celltype_overlap_summary.csv")),
  row.names = FALSE
)

cat("\nClone by cell type:\n")
print(clone_celltype_table %>% filter(numbat_clone == clone_to_check))

cat("\nClone by binary overlap:\n")
print(clone_overlap_table %>% filter(numbat_clone == clone_to_check))

# ============================================================
# 2. UMAPS: clone 1, lymphoid status, and overlap
# ============================================================

clone_focus_cols <- c(
  "#E41A1C",
  "grey75",
  "grey90"
)

names(clone_focus_cols) <- c(
  paste0("Clone_", clone_to_check),
  "Other_clones",
  "No_clone"
)

status_cols <- c(
  "#E41A1C",
  "#FF7F00",
  "#377EB8",
  "grey85"
)

names(status_cols) <- c(
  paste0("Clone_", clone_to_check, "_lymphoid"),
  paste0("Clone_", clone_to_check, "_non_lymphoid"),
  "Other_lymphoid",
  "Other_non_lymphoid"
)

p_clone_focus <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "clone_focus",
  cols = clone_focus_cols,
  raster = FALSE
) +
  ggtitle(paste0("Numbat clone ", clone_to_check, " vs other cells"))

p_clone_status <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "clone1_lymphoid_status",
  cols = status_cols,
  raster = FALSE
) +
  ggtitle(paste0("Clone ", clone_to_check, ": lymphoid vs non-lymphoid"))

p_overlap <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = "cnv_binary_overlap",
  raster = FALSE
) +
  ggtitle("Binary Numbat/CopyKAT overlap")

p_celltype <- DimPlot(
  seu_sub,
  reduction = reduction_use,
  group.by = celltype_col,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Cell type")

ggsave(
  file.path(clone_out_dir, paste0("04_clone_", clone_to_check, "_inspection_umaps.pdf")),
  (p_clone_focus | p_clone_status) / (p_overlap | p_celltype),
  width = 14,
  height = 11
)

# ============================================================
# 3. QC: RNA depth and features
# ============================================================

qc_features <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  colnames(seu_sub@meta.data)
)

if (length(qc_features) > 0) {
  
  p_qc_clone <- VlnPlot(
    seu_sub,
    features = qc_features,
    group.by = "clone1_lymphoid_status",
    pt.size = 0,
    ncol = length(qc_features)
  ) +
    ggtitle(paste0("QC metrics for clone ", clone_to_check, " lymphoid calls"))
  
  ggsave(
    file.path(clone_out_dir, paste0("05_clone_", clone_to_check, "_qc_violin.pdf")),
    p_qc_clone,
    width = 14,
    height = 5
  )
}

qc_table <- meta %>%
  group_by(clone1_lymphoid_status) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = if ("percent.mt" %in% colnames(meta)) median(percent.mt, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

write.csv(
  qc_table,
  file.path(clone_out_dir, paste0("06_clone_", clone_to_check, "_qc_summary.csv")),
  row.names = FALSE
)

# ============================================================
# 4. Look for Numbat confidence / posterior columns
# ============================================================

possible_conf_cols <- colnames(meta)[
  grepl(
    "prob|post|posterior|conf|lik|score|entropy|p_",
    colnames(meta),
    ignore.case = TRUE
  )
]

writeLines(
  possible_conf_cols,
  file.path(clone_out_dir, "07_possible_numbat_confidence_columns.txt")
)

cat("\nPossible Numbat confidence/posterior columns:\n")
print(possible_conf_cols)

# plot numeric confidence-like columns if present
numeric_conf_cols <- possible_conf_cols[
  possible_conf_cols %in% colnames(meta) &
    sapply(meta[possible_conf_cols], is.numeric)
]

if (length(numeric_conf_cols) > 0) {
  
  for (cc in numeric_conf_cols) {
    seu_sub[[cc]] <- meta[[cc]][match(colnames(seu_sub), meta$cell)]
  }
  
  pdf(
    file.path(clone_out_dir, "08_possible_confidence_featureplots.pdf"),
    width = 7,
    height = 6
  )
  
  for (cc in numeric_conf_cols) {
    print(
      FeaturePlot(
        seu_sub,
        features = cc,
        reduction = reduction_use,
        order = TRUE,
        raster = FALSE
      ) +
        ggtitle(cc)
    )
  }
  
  dev.off()
  
  conf_summary <- meta %>%
    group_by(clone1_lymphoid_status) %>%
    summarise(
      across(
        all_of(numeric_conf_cols),
        list(
          median = ~ median(.x, na.rm = TRUE),
          mean = ~ mean(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_cells = n(),
      .groups = "drop"
    )
  
  write.csv(
    conf_summary,
    file.path(clone_out_dir, "09_possible_confidence_summary_by_status.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 5. Try to load useful Numbat output files for segment/posterior inspection
# ============================================================

numbat_all_files <- list.files(
  numbat_dir,
  full.names = TRUE,
  recursive = TRUE
)

numbat_candidate_files <- numbat_all_files[
  grepl(
    "post|posterior|cnv|seg|bulk|clone|loh|allele|haplo",
    basename(numbat_all_files),
    ignore.case = TRUE
  )
]

candidate_file_table <- tibble(
  file = numbat_candidate_files,
  basename = basename(numbat_candidate_files)
)

write.csv(
  candidate_file_table,
  file.path(clone_out_dir, "10_numbat_candidate_files.csv"),
  row.names = FALSE
)

cat("\nCandidate Numbat files:\n")
print(candidate_file_table)

# ============================================================
# 6. If posterior-like files are readable, summarize clone 1 rows
# ============================================================

readable_tables <- list()

for (ff in numbat_candidate_files) {
  
  dat <- tryCatch(
    fread(ff, nThread = 1),
    error = function(e) NULL
  )
  
  if (is.null(dat)) next
  if (nrow(dat) == 0) next
  
  readable_tables[[basename(ff)]] <- dat
  
  out_preview <- as.data.frame(head(dat, 20))
  
  write.csv(
    out_preview,
    file.path(
      clone_out_dir,
      paste0("preview_", gsub("[^A-Za-z0-9_]", "_", basename(ff)), ".csv")
    ),
    row.names = FALSE
  )
}

# save column names from readable tables
if (length(readable_tables) > 0) {
  
  readable_cols <- imap_dfr(
    readable_tables,
    ~ tibble(
      file = .y,
      column = colnames(.x)
    )
  )
  
  write.csv(
    readable_cols,
    file.path(clone_out_dir, "11_readable_numbat_table_columns.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 7. Try to find per-cell posterior table and merge if possible
# ============================================================

posterior_tables <- readable_tables[
  sapply(readable_tables, function(x) {
    any(colnames(x) %in% c("cell", "barcode")) &&
      any(grepl("prob|post|posterior|clone|p_", colnames(x), ignore.case = TRUE))
  })
]

if (length(posterior_tables) > 0) {
  
  for (nm in names(posterior_tables)) {
    
    tab <- posterior_tables[[nm]]
    
    if (!"cell" %in% colnames(tab) && "barcode" %in% colnames(tab)) {
      tab <- tab %>% rename(cell = barcode)
    }
    
    clone1_cells <- meta %>%
      filter(numbat_clone == clone_to_check) %>%
      pull(cell)
    
    tab_clone1 <- tab %>%
      filter(cell %in% clone1_cells)
    
    write.csv(
      tab_clone1,
      file.path(
        clone_out_dir,
        paste0("12_clone_", clone_to_check, "_rows_from_", gsub("[^A-Za-z0-9_]", "_", nm), ".csv")
      ),
      row.names = FALSE
    )
  }
}

# ============================================================
# 8. Compare marker/cell-type identity of clone 1 lymphoid cells
# ============================================================

Idents(seu_sub) <- "clone1_lymphoid_status"

groups_present <- unique(seu_sub$clone1_lymphoid_status)

if (
  paste0("Clone_", clone_to_check, "_lymphoid") %in% groups_present &&
  "Other_lymphoid" %in% groups_present
) {
  
  de_clone1_lymphoid_vs_other_lymphoid <- FindMarkers(
    seu_sub,
    ident.1 = paste0("Clone_", clone_to_check, "_lymphoid"),
    ident.2 = "Other_lymphoid",
    assay = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05
  ) %>%
    rownames_to_column("gene") %>%
    arrange(p_val_adj)
  
  write.csv(
    de_clone1_lymphoid_vs_other_lymphoid,
    file.path(clone_out_dir, paste0("13_clone_", clone_to_check, "_lymphoid_vs_other_lymphoid_DE.csv")),
    row.names = FALSE
  )
}

# ============================================================
# 9. Final summary to console
# ============================================================

cat("\nDone clone inspection.\n")
cat("Output directory:\n", clone_out_dir, "\n\n")

cat("Clone 1 lymphoid status counts:\n")
print(table(meta$clone1_lymphoid_status, useNA = "ifany"))

cat("\nQC summary:\n")
print(qc_table)

