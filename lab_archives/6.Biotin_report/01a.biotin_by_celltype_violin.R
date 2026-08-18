#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Biotin report - step 01a
#
# Biotin by broad cell type, split by sample, with the cell types ordered by
# median biotin in the AML sample rather than alphabetically.
#
# WHERE THIS CAME FROM
#
# This figure existed only as a slide PNG. The code that made it was four
# unassigned lines at the very end of 6b1.biotin_ADT_comparison.R:
#
#   VlnPlot(seu, features = dsb_markers, assay = "CITE_DSB", slot = "data",
#           group.by = "predicted_CellType_Broad", split.by = "sample_name",
#           pt.size = 0, cols = c("gray50","gray70","firebrick")) + ...
#
# It was never assigned, never saved, and `seu` at that point was whatever the
# loop above had left in the session - which is not the LK3 object the figure
# actually shows. So the figure could not be regenerated from the script. This
# step loads the LK3 object explicitly and writes the figure to disk.
#
# WHAT CHANGED
#
# Seurat orders the x axis by factor level, which for a character column is
# alphabetical - that is why the original ran B, CD4 Memory T, CD8 Memory T,
# cDC, ... The cell types are now ordered by median biotin in the AML sample
# (descending), so the panel reads as a ranking. The ordering sample is an
# argument, and the ordering itself is written out as a table so the figure and
# the numbers cannot drift apart.
#
# NA cells - those with no predicted_CellType_Broad - are dropped rather than
# plotted as an "NA" category, which is what the original showed at the right
# hand end.
#
# Run from the scripts/ directory - paths are relative to it.
#
#   Rscript lab_archives/6.Biotin_report/01a.biotin_by_celltype_violin.R
#   Rscript lab_archives/6.Biotin_report/01a.biotin_by_celltype_violin.R WEI21_26-17_NK
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

ORDER_BY_SAMPLE <- if (length(args) >= 1) args[1] else "WEI21_26-17_NK"
MARKER          <- "biotin"
ASSAY           <- "CITE_DSB"

SEURAT_FILE <- "../results_nf/260717_VH01624_477_222KG22NX/rds/11_final/LK3-GEX/rds/LK3_final.rds"
OUT_DIR     <- "../results/biotin_by_celltype"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SEURAT_FILE)) stop("Cannot find LK3 object: ", SEURAT_FILE)

message("Loading: ", SEURAT_FILE)
seu <- readRDS(SEURAT_FILE)

stopifnot(ASSAY %in% Assays(seu))
stopifnot("predicted_CellType_Broad" %in% colnames(seu@meta.data))
stopifnot("sample_name" %in% colnames(seu@meta.data))

# ----------------------------
# Pull the biotin values
# ----------------------------
# The ADT feature name is matched case-insensitively - it has been "biotin" and
# "Biotin" in different panel versions.

adt_features <- rownames(GetAssayData(seu, assay = ASSAY, layer = "data"))
marker_hit <- adt_features[grepl(paste0("^", MARKER, "$"), adt_features, ignore.case = TRUE)]

if (!length(marker_hit)) {
  marker_hit <- adt_features[grepl(MARKER, adt_features, ignore.case = TRUE)]
}
if (!length(marker_hit)) {
  stop("No ADT feature matching '", MARKER, "' in assay ", ASSAY,
       ". Available: ", paste(head(adt_features, 40), collapse = ", "))
}
marker_hit <- marker_hit[1]
message("Using ADT feature: ", marker_hit)

df <- tibble(
  cell      = colnames(seu),
  value     = as.numeric(GetAssayData(seu, assay = ASSAY, layer = "data")[marker_hit, ]),
  celltype  = as.character(seu$predicted_CellType_Broad),
  sample    = as.character(seu$sample_name)
) %>%
  filter(!is.na(celltype), celltype != "", !is.na(sample))

message("Cells plotted: ", nrow(df))
message("Samples: ", paste(sort(unique(df$sample)), collapse = ", "))

if (!ORDER_BY_SAMPLE %in% df$sample) {
  stop("Ordering sample '", ORDER_BY_SAMPLE, "' not present. Available: ",
       paste(sort(unique(df$sample)), collapse = ", "))
}

# ----------------------------
# Order the cell types by median biotin in the ordering sample
# ----------------------------
# Cell types absent from the ordering sample cannot be ranked by it. They are
# kept, ranked among themselves by their overall median, and placed at the end
# so it is obvious they are not part of the ranking.

MIN_CELLS <- 10

rank_tbl <- df %>%
  filter(sample == ORDER_BY_SAMPLE) %>%
  group_by(celltype) %>%
  summarise(n_cells = n(), median_biotin = median(value), .groups = "drop") %>%
  filter(n_cells >= MIN_CELLS) %>%
  arrange(desc(median_biotin))

unranked <- df %>%
  filter(!celltype %in% rank_tbl$celltype) %>%
  group_by(celltype) %>%
  summarise(n_cells = n(), median_biotin = median(value), .groups = "drop") %>%
  arrange(desc(median_biotin))

if (nrow(unranked)) {
  message("Cell types with fewer than ", MIN_CELLS, " cells in ", ORDER_BY_SAMPLE,
          " - placed at the end, not ranked: ",
          paste(unranked$celltype, collapse = ", "))
}

celltype_order <- c(rank_tbl$celltype, unranked$celltype)

write_csv(
  bind_rows(
    rank_tbl %>% mutate(ranked_by = ORDER_BY_SAMPLE),
    unranked %>% mutate(ranked_by = "not ranked - too few cells in ordering sample")
  ),
  file.path(OUT_DIR, "biotin_median_by_celltype_ordering.csv")
)

df$celltype <- factor(df$celltype, levels = celltype_order)

# ----------------------------
# Sample colours
# ----------------------------
# Same convention as the original: normals grey, the AML sample red.

samples <- sort(unique(df$sample))
normals <- setdiff(samples, ORDER_BY_SAMPLE)

sample_cols <- setNames(
  c(rep(c("gray50", "gray70"), length.out = length(normals)), "firebrick"),
  c(normals, ORDER_BY_SAMPLE)
)

df$sample <- factor(df$sample, levels = c(normals, ORDER_BY_SAMPLE))

# ----------------------------
# Plot
# ----------------------------

p <- ggplot(df, aes(x = celltype, y = value, fill = sample)) +
  geom_violin(scale = "width", linewidth = 0.25, colour = "black",
              position = position_dodge(width = 0.85), trim = TRUE) +
  scale_fill_manual(values = sample_cols, name = NULL) +
  theme_classic(base_size = 12) +
  labs(
    title = marker_hit,
    subtitle = paste0("broad cell types ordered by median ", marker_hit,
                      " in ", ORDER_BY_SAMPLE),
    x = "Identity", y = "Expression Level"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "biotin_by_celltype_violin_ordered.pdf"),
       p, width = 14, height = 5.5)

# A median-only version - with 20-odd cell types the violins are narrow, and
# the ranking is the point
p_med <- rank_tbl %>%
  mutate(celltype = factor(celltype, levels = rank_tbl$celltype)) %>%
  ggplot(aes(celltype, median_biotin)) +
  geom_col(fill = "firebrick", width = 0.75) +
  geom_text(aes(label = n_cells), vjust = -0.4, size = 3, colour = "grey30") +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Median ", marker_hit, " by broad cell type in ", ORDER_BY_SAMPLE),
    subtitle = "cell counts above each bar",
    x = NULL, y = paste0("Median ", marker_hit, " (", ASSAY, ")")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "biotin_median_by_celltype_ranked.pdf"),
       p_med, width = 11, height = 5)

cat("\nDone. Output:", OUT_DIR, "\n\n")
cat("Cell type ranking by median", marker_hit, "in", ORDER_BY_SAMPLE, ":\n")
print(as.data.frame(rank_tbl))
