#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

run <- "260528_VH01624_464_222K7VKNX"

annotation_dir <- file.path(
  proj,
  "results",
  "seurat_annotated",
  run
)

seurat_file <- file.path(
  annotation_dir,
  "demux_singlets_annotated_seurat.rds"
)

sig_file <- file.path(
  proj,
  "annotation",
  "van_galen_signatures_compartment_matched",
  "van_galen_AUC_signatures_top_compartment_matched.rds"
)

out_dir <- file.path(
  proj,
  "results",
  "van_galen_signature_scores"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------

celltype_col <- "predicted_CellType"
status_col   <- "AML_status_numbat"

# ----------------------------
# Load data
# ----------------------------

seu <- readRDS(seurat_file)
sigs <- readRDS(sig_file)

DefaultAssay(seu) <- "RNA"

# ----------------------------
# Make broad compartments
# ----------------------------

meta <- seu@meta.data

seu$major_compartment <- case_when(
  meta[[celltype_col]] %in% c(
    "HSC", "HSC MPP", "MPP", "LMPP", "Early Lymphoid"
  ) ~ "HSC_Prog",
  
  meta[[celltype_col]] %in% c(
    "GMP", "GMP Early", "GMP Late", "Pro-Monocyte",
    "Cycling Progenitor"
  ) ~ "GMP",
  
  meta[[celltype_col]] %in% c(
    "CD14 Mono", "CD16 Mono", "Monocyte", "cDC", "pDC"
  ) ~ "Myeloid",
  
  TRUE ~ NA_character_
)

# ----------------------------
# Score signatures
# ----------------------------

for (nm in names(sigs)) {
  genes_use <- intersect(sigs[[nm]], rownames(seu))
  
  message(nm, ": ", length(genes_use), " genes found in object")
  
  if (length(genes_use) < 5) {
    warning("Skipping ", nm, " because fewer than 5 genes were found.")
    next
  }
  
  seu <- AddModuleScore(
    object = seu,
    features = list(genes_use),
    name = paste0(nm, "_")
  )
}

# AddModuleScore creates columns ending in 1
score_cols <- grep("^VanGalen_.*_compartment_matched_1$", colnames(seu@meta.data), value = TRUE)

print(score_cols)

# Rename to cleaner names
rename_map <- c(
  "VanGalen_HSC_Prog_compartment_matched_1" = "VG_HSC_Prog_score",
  "VanGalen_GMP_compartment_matched_1"      = "VG_GMP_score",
  "VanGalen_Myeloid_compartment_matched_1"  = "VG_Myeloid_score"
)

for (old in names(rename_map)) {
  if (old %in% colnames(seu@meta.data)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == old] <- rename_map[[old]]
  }
}

# ----------------------------
# UMAP plots
# ----------------------------

score_cols_clean <- intersect(
  c("VG_HSC_Prog_score", "VG_GMP_score", "VG_Myeloid_score"),
  colnames(seu@meta.data)
)

for (score in score_cols_clean) {
  
  p <- FeaturePlot(
    seu,
    features = score,
    reduction = "umap",
    order = TRUE
  ) +
    ggtitle(score)
  
  ggsave(
    file.path(out_dir, paste0(score, "_UMAP.png")),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# ----------------------------
# Matched compartment violin plots
# ----------------------------

plot_matched_score <- function(seu, compartment, score_col) {
  
  cells_use <- rownames(seu@meta.data)[
    seu$major_compartment == compartment &
      !is.na(seu@meta.data[[status_col]]) &
      !is.na(seu@meta.data[[score_col]])
  ]
  
  if (length(cells_use) == 0) {
    warning("No cells found for ", compartment)
    return(NULL)
  }
  
  df <- seu@meta.data[cells_use, , drop = FALSE] %>%
    tibble::rownames_to_column("cell") %>%
    mutate(score = .data[[score_col]])
  
  ggplot(df, aes(x = .data[[status_col]], y = score, fill = .data[[status_col]])) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste0(score_col, " in ", compartment),
      x = status_col,
      y = "Module score"
    )
}

matched_plots <- list(
  HSC_Prog = plot_matched_score(seu, "HSC_Prog", "VG_HSC_Prog_score"),
  GMP      = plot_matched_score(seu, "GMP", "VG_GMP_score"),
  Myeloid  = plot_matched_score(seu, "Myeloid", "VG_Myeloid_score")
)

for (nm in names(matched_plots)) {
  
  p <- matched_plots[[nm]]
  
  if (!is.null(p)) {
    ggsave(
      file.path(out_dir, paste0("matched_", nm, "_signature_score_violin.png")),
      p,
      width = 5,
      height = 4,
      dpi = 300
    )
  }
}

# ----------------------------
# Combined plot
# ----------------------------

combined <- wrap_plots(matched_plots[!vapply(matched_plots, is.null, logical(1))], ncol = 3)

ggsave(
  file.path(out_dir, "matched_compartment_signature_scores_combined.png"),
  combined,
  width = 14,
  height = 4,
  dpi = 300
)

# ----------------------------
# Save scored object
# ----------------------------

saveRDS(
  seu,
  file.path(out_dir, "seurat_with_van_galen_compartment_scores.rds")
)

cat("\nSaved plots and scored Seurat object to:\n")
cat(out_dir, "\n")