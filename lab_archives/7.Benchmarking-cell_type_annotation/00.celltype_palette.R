#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 00 of 05
#
# The project's bespoke broad-cell-type palette, in one place, plus the
# reference UMAPs drawn with it. Steps 01-03 source this file rather than
# each carrying their own copy, which is how the palettes drifted apart in
# the first place.
#
# WHY THIS EXISTS AS ITS OWN STEP
#
# The palette used to list "GMP" - a label the BoneMarrowMap reference never
# emits, it splits that population into Early GMP and Late GMP - and it had no
# entry for "Stromal", which the reference does emit. Every plot then did
#
#   factor(CellType_Broad, levels = names(celltype_cols))
#
# which silently sets every Stromal cell to NA. That is why stromal cells
# looked unannotated: they were annotated, they were being dropped at the
# plotting step. The reference carries a whole Stromal population and our own
# data carries 11 cells of it.
#
# Stromal is grey because it is not haematopoietic and should read as
# background rather than as another lineage. It is grey25 against the grey70
# used for "Unknown" in step 03 - dark against light, so that in a plot
# carrying both it is obvious which grey is which, and so the Unknown cells
# can be told apart from the stromal island they sit near.
#
# Run directly to regenerate the reference UMAPs:
#   Rscript lab_archives/7.Benchmarking-cell_type_annotation/00.celltype_palette.R
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

celltype_cols <- c(
  # Stem / progenitor - greens
  "HSC MPP" = "#1B9E77",
  "LMPP" = "#66C2A5",
  "MEP" = "#B2DF8A",
  "GMP" = "#33A02C",
  "Early GMP" = "#A6D854",
  "Late GMP" = "#006D2C",
  "Cycling Progenitor" = "#00441B",
  "EoBasoMast Precursor" = "#8DD3C7",
  "Megakaryocyte Precursor" = "#4DAF4A",

  # Myeloid / DC - oranges/reds
  "Monocyte" = "#E31A1C",
  "Pro-Monocyte" = "#FB6A4A",
  "cDC" = "#FD8D3C",
  "pDC" = "#FCBBA1",

  # Lymphoid - blues/purples
  "Naive T" = "#2171B5",
  "CD4 Memory T" = "#6BAED6",
  "CD8 Memory T" = "#08519C",
  "NK" = "#54278F",
  "Early Lymphoid" = "#9E9AC8",
  "B" = "#3182BD",
  "Pre-B" = "#9ECAE1",
  "Pro-B" = "#C6DBEF",
  "Plasma Cell" = "#756BB1",

  # Erythroid - pinks
  "Early Erythroid" = "#F768A1",
  "Late Erythroid" = "#C51B8A",

  # Non-haematopoietic - grey, so it reads as background
  "Stromal" = "grey25"
)

celltype_order <- rev(c(
  "HSC MPP", "LMPP", "MEP", "Megakaryocyte Precursor", "GMP", "Early GMP",
  "Late GMP", "Cycling Progenitor", "EoBasoMast Precursor",
  "Pro-Monocyte", "Monocyte", "cDC", "pDC",
  "Early Lymphoid", "Pro-B", "Pre-B", "B", "Naive T", "CD4 Memory T",
  "CD8 Memory T", "NK", "Plasma Cell",
  "Early Erythroid", "Late Erythroid",
  "Stromal"
))

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Stromal" = "grey25",
  "Other" = "grey70"
)

# Lineage grouping. Keyed on predicted_CellType_Broad - see step 04 for what
# happens if it is keyed on the fine label instead.
lineage_of <- function(x) {
  dplyr::case_when(
    x %in% c("HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
             "Cycling Progenitor", "EoBasoMast Precursor",
             "Megakaryocyte Precursor") ~ "Stem / progenitor",
    x %in% c("Monocyte", "Pro-Monocyte", "cDC", "pDC") ~ "Myeloid / DC",
    x %in% c("Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
             "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell") ~ "Lymphoid",
    x %in% c("Early Erythroid", "Late Erythroid") ~ "Erythroid",
    x %in% c("Stromal") ~ "Stromal",
    TRUE ~ "Other"
  )
}

# Guard: anything the data carries that the palette does not know about would
# be silently dropped by factor(levels = names(celltype_cols)). Warn loudly
# instead, and hand back a palette extended with visibly distinct fallbacks.
check_palette <- function(labels, cols = celltype_cols) {
  labels <- unique(as.character(labels[!is.na(labels)]))
  missing <- setdiff(labels, names(cols))
  if (length(missing)) {
    warning("Cell type(s) with no palette entry - they would be dropped from ",
            "every plot: ", paste(missing, collapse = ", "), call. = FALSE)
    cols <- c(cols, setNames(scales::hue_pal()(length(missing)), missing))
  }
  cols
}

# ----------------------------
# Reference UMAPs
# ----------------------------
# Only runs when this file is executed directly, not when it is sourced.

if (sys.nframe() == 0L) {

  suppressPackageStartupMessages({
    library(Seurat); library(ggplot2); library(BoneMarrowMap); library(symphony)
  })

  projection_path <- "../annotation/"
  ref <- readRDS(paste0(projection_path, "BoneMarrowMap_SymphonyReference.rds"))
  ref$save_uwot_path <- paste0(projection_path, "BoneMarrowMap_uwot_model.uwot")

  ReferenceSeuratObj <- create_ReferenceObject(ref)

  cols <- check_palette(ReferenceSeuratObj$CellType_Broad)

  ReferenceSeuratObj$CellType_Broad <- factor(
    ReferenceSeuratObj$CellType_Broad,
    levels = intersect(names(cols), unique(as.character(ReferenceSeuratObj$CellType_Broad)))
  )

  stopifnot(!any(is.na(ReferenceSeuratObj$CellType_Broad)))

  p_ref <- DimPlot(
    ReferenceSeuratObj, reduction = "umap", group.by = "CellType_Broad",
    raster = FALSE, label = TRUE, repel = TRUE, label.size = 5, cols = cols
  ) +
    ggtitle("BoneMarrowMap reference: broad cell type") +
    theme_bw()

  targets <- c(
    "../results/seurat_annotated/260423_VH01624_453_222HWMYNX/00a_BoneMarrowMap_reference_broad_celltype.pdf",
    "../results/seurat_annotated/260528_VH01624_464_222K7VKNX/11.database_umap.pdf"
  )

  for (f in targets) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    pdf(f, width = 12, height = 10)
    print(p_ref)
    dev.off()
    message("Wrote ", f)
  }

  # The fine-label version keeps default colours - there are 53 of them and no
  # bespoke palette exists at that level.
  p_fine <- DimPlot(
    ReferenceSeuratObj, reduction = "umap",
    group.by = "CellType_Annotation_formatted",
    raster = FALSE, label = TRUE, repel = TRUE, label.size = 4
  ) +
    ggtitle("BoneMarrowMap reference: formatted annotation")

  f_fine <- "../results/seurat_annotated/260423_VH01624_453_222HWMYNX/00_BoneMarrowMap_reference_celltype.pdf"
  pdf(f_fine, width = 12, height = 10)
  print(p_fine)
  dev.off()
  message("Wrote ", f_fine)

  cat("\nBroad cell types drawn:", nlevels(ReferenceSeuratObj$CellType_Broad), "\n")
  print(table(ReferenceSeuratObj$CellType_Broad))
}
