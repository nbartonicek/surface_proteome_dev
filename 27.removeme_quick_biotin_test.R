suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(patchwork)
})

seu <- readRDS("../results_nf/260717_VH01624_477_222KG22NX/rds/11_final/LK3-GEX/rds/LK3_final.rds")

# DSB-normalized CITE-seq assay as the active assay
DefaultAssay(seu) <- "CITE_DSB"

# all protein markers in the DSB assay
dsb_markers <- "CD34"

# one marker, split into one panel per patient (sample_name) on the same
# umap_projected layout - keep.scale = "feature" forces every patient's
# panel to share the same color scale, which is what actually makes them
# comparable (Seurat's default per-panel scaling would otherwise rescale
# min/max independently per patient, making side-by-side comparison
# meaningless). plot_annotation() adds the marker name as a figure-level
# title without overwriting each panel's own per-patient title.
FeaturePlot(seu, features = dsb_markers[1], reduction = "umap_projected",
            min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.3,
            split.by = "sample_name", keep.scale = "feature") +
  plot_annotation(title = dsb_markers[1])

# grid of every DSB marker (rows) x patient (columns) on the same projected
# UMAP, same shared-scale-per-marker logic as above. `ncol` dropped since
# Seurat ignores it once split.by is set - layout is fixed to markers x
# patients.
FeaturePlot(seu, features = dsb_markers, reduction = "umap_projected",
            min.cutoff = "q05", max.cutoff = "q95", pt.size = 0.1,
            split.by = "sample_name", keep.scale = "feature") &
  theme(legend.position = "right")

# violin of DSB values by cell type, one violin per patient placed side by
# side within each cell type - split.by (not split.plot, which instead
# overlays each split as a half-violin and only really works for exactly 2
# groups) draws adjacent violins per patient within every
# predicted_CellType_Broad group, same "compare across patients" intent as
# the FeaturePlots above. pt.size = 0 since jittering every cell as a point
# is unreadable/slow at this cell count.
VlnPlot(seu, features = dsb_markers, assay = "CITE_DSB", slot = "data",
        group.by = "predicted_CellType_Broad", split.by = "sample_name",
        pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))