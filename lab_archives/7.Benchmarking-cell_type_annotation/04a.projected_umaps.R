#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 04a of 05
#
# The projected UMAP for both runs, coloured by predicted_CellType_Broad in
# the project palette from step 00.
#
# The originals in results/seurat_annotated/ are coloured by the FINE label
# (53 levels, ggplot default colours), which is unreadable and does not match
# any other figure in the project. They also only exist for one run at a time.
# This step draws the same thing for LK1 and LK2 side by side in the bespoke
# palette, so the two are comparable and so Stromal is visible rather than
# silently dropped.
#
# Loads the two projected Seurat objects, which are ~450 MB each - this is the
# only step in the folder that does, and it is why it is separate from step 04.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

source("lab_archives/7.Benchmarking-cell_type_annotation/00.celltype_palette.R")

res     <- "../results"
out_dir <- file.path(res, "benchmarking", "celltype_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PILOT <- "260423_VH01624_453_222HWMYNX"
LK2   <- "260528_VH01624_464_222K7VKNX"

objects <- list(
  "LK1 pilot" = file.path(res, "seurat_annotated", PILOT,
                          "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"),
  "LK2"       = file.path(res, "seurat_annotated", LK2,
                          "demux_singlets_annotated_seurat.rds")
)

plot_one <- function(path, label) {

  if (!file.exists(path)) {
    message("Missing object, skipping: ", path)
    return(NULL)
  }

  message("Loading ", label, ": ", path)
  seu <- readRDS(path)

  red <- if ("umap_projected" %in% Reductions(seu)) "umap_projected" else "umap"
  message("  reduction: ", red, "; cells: ", ncol(seu))

  seu <- subset(seu, subset = mapping_error_QC == "Pass")

  cols <- check_palette(seu$predicted_CellType_Broad)

  seu$predicted_CellType_Broad <- factor(
    seu$predicted_CellType_Broad,
    levels = intersect(names(cols), unique(as.character(seu$predicted_CellType_Broad)))
  )

  present <- table(seu$predicted_CellType_Broad)
  message("  broad types present: ", sum(present > 0))
  print(present[present > 0])

  p <- DimPlot(
    seu, reduction = red, group.by = "predicted_CellType_Broad",
    cols = cols, raster = FALSE, label = TRUE, repel = TRUE, label.size = 3.5
  ) +
    ggtitle(paste0(label, " - projected UMAP, broad cell type")) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, paste0("07_projected_umap_broad_",
                                   gsub("[^A-Za-z0-9]+", "_", label), ".pdf")),
         p, width = 10, height = 8)

  rm(seu); gc(verbose = FALSE)
  p
}

plots <- imap(objects, ~ plot_one(.x, .y))
plots <- compact(plots)

if (length(plots) == 2) {
  combined <- plots[[1]] + plots[[2]] + plot_layout(ncol = 2, guides = "collect")
  ggsave(file.path(out_dir, "07_projected_umap_broad_both_runs.pdf"),
         combined, width = 18, height = 8)
  message("Wrote the combined two-run panel.")
} else {
  message("Only ", length(plots), " object(s) available - no combined panel written.")
}

cat("\nDone. Output directory:", out_dir, "\n")
