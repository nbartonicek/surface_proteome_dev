#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(SingleCellExperiment)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

aml_runs <- c(
  "260528_VH01624_464_222K7VKNX"
)

sig_dir <- file.path(proj, "annotation", "van_galen_signatures")
out_dir <- file.path(proj, "results", "van_galen_analysis")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")

for (d in c(fig_dir, tab_dir, rds_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Load signatures
# ------------------------------------------------------------------

signatures_top <- readRDS(file.path(sig_dir, "van_galen_signatures_top50.rds"))
signatures_all <- readRDS(file.path(sig_dir, "van_galen_signatures_all_positive.rds"))

cat("Signatures loaded:\n")
for (nm in names(signatures_top)) {
  cat(sprintf("  %s: %d (top50), %d (all positive)\n", nm, length(signatures_top[[nm]]), length(signatures_all[[nm]])))
}

# ------------------------------------------------------------------
# 2. Load annotated Seurat objects (with numbat calls)
# ------------------------------------------------------------------

seu_list <- list()
for (run in aml_runs) {
  rds_file <- file.path(
    proj, "results", "seurat_annotated", run,
    "numbat", list.files(
      file.path(proj, "results", "seurat_annotated", run, "numbat"),
      pattern = "_projected_CITE_DSB_Numbat_integrated\\.rds$"
    )[1]
  )

  if (!is.null(rds_file) && file.exists(rds_file)) {
    cat("Loading:", rds_file, "\n")
    seu_list[[run]] <- readRDS(rds_file)
  } else {
    cat("WARNING: No numbat-annotated RDS found for", run, "\n")
  }
}

if (length(seu_list) == 0) stop("No Seurat objects loaded")

if (length(seu_list) == 1) {
  seu <- seu_list[[1]]
} else {
  seu <- merge(seu_list[[1]], seu_list[-1], add.cell.ids = paste0("run", seq_along(seu_list)))
}

DefaultAssay(seu) <- "RNA"
if (inherits(seu[["RNA"]], "Assay5")) {
  seu <- JoinLayers(seu, assay = "RNA")
}

cat("Total cells:", ncol(seu), "\n")

# ------------------------------------------------------------------
# 3. Run AUCell
# ------------------------------------------------------------------

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")

cells_rankings <- AUCell_buildRankings(counts, plotStats = FALSE)

gene_sets <- c(signatures_top, setNames(signatures_all, paste0(names(signatures_all), "_all")))

cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(counts)))

auc_mat <- getAUC(cells_AUC)

for (sig_name in rownames(auc_mat)) {
  safe_name <- gsub("/", "_", sig_name)
  seu@meta.data[[paste0("AUCell_", safe_name)]] <- auc_mat[sig_name, colnames(seu)]
}

cat("AUCell scores added for", nrow(auc_mat), "signatures\n")

# ------------------------------------------------------------------
# 4. Define cancer vs non-cancer based on numbat
# ------------------------------------------------------------------

if ("numbat_compartment" %in% colnames(seu@meta.data)) {
  seu$cancer_status <- case_when(
    seu$numbat_compartment == "tumor" ~ "Cancer",
    seu$numbat_compartment == "normal" ~ "Normal",
    TRUE ~ "Not_called"
  )
} else if ("numbat_call" %in% colnames(seu@meta.data)) {
  seu$cancer_status <- case_when(
    grepl("aneuploid|tumor", seu$numbat_call, ignore.case = TRUE) ~ "Cancer",
    grepl("diploid|normal|reference", seu$numbat_call, ignore.case = TRUE) ~ "Normal",
    TRUE ~ "Not_called"
  )
} else {
  seu$cancer_status <- "Not_called"
  message("WARNING: No numbat compartment/call column found")
}

cat("\nCancer status breakdown:\n")
print(table(seu$cancer_status, useNA = "ifany"))

# ------------------------------------------------------------------
# 5. Violin plots: AUCell scores by cancer status
# ------------------------------------------------------------------

aucell_cols <- grep("^AUCell_", colnames(seu@meta.data), value = TRUE)
top50_cols <- aucell_cols[!grepl("_all$", aucell_cols)]

for (col in top50_cols) {
  sig_label <- gsub("AUCell_", "", col)

  p <- VlnPlot(seu, features = col, group.by = "cancer_status", pt.size = 0) +
    ggtitle(paste("Van Galen", sig_label, "signature")) +
    ylab("AUCell score") +
    theme(legend.position = "none")

  ggsave(
    file.path(fig_dir, paste0("01_violin_", sig_label, "_by_cancer_status.pdf")),
    p, width = 6, height = 5
  )
}

# ------------------------------------------------------------------
# 6. Violin plots: AUCell scores by cancer status, split by sample
# ------------------------------------------------------------------

if ("sample_name" %in% colnames(seu@meta.data)) {
  for (col in top50_cols) {
    sig_label <- gsub("AUCell_", "", col)

    p <- VlnPlot(
      seu,
      features = col,
      group.by = "cancer_status",
      split.by = "sample_name",
      pt.size = 0
    ) +
      ggtitle(paste("Van Galen", sig_label, "by sample")) +
      ylab("AUCell score")

    ggsave(
      file.path(fig_dir, paste0("02_violin_", sig_label, "_by_cancer_status_per_sample.pdf")),
      p, width = 10, height = 5
    )
  }
}

# ------------------------------------------------------------------
# 7. AUCell heatmap: mean score per sample x signature
# ------------------------------------------------------------------

sample_col <- if ("sample_name" %in% colnames(seu@meta.data)) "sample_name" else "orig.ident"

aucell_summary <- seu@meta.data %>%
  dplyr::select(all_of(c(sample_col, "cancer_status", aucell_cols))) %>%
  pivot_longer(cols = all_of(aucell_cols), names_to = "signature", values_to = "score") %>%
  group_by(.data[[sample_col]], cancer_status, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

write.csv(aucell_summary, file.path(tab_dir, "aucell_mean_scores_by_sample_status.csv"), row.names = FALSE)

aucell_per_sample <- seu@meta.data %>%
  dplyr::select(all_of(c(sample_col, aucell_cols))) %>%
  group_by(.data[[sample_col]]) %>%
  summarise(across(all_of(aucell_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames(sample_col) %>%
  as.matrix()

pdf(file.path(fig_dir, "03_aucell_heatmap_per_sample.pdf"), width = 10, height = 8)
pheatmap(
  t(aucell_per_sample),
  scale = "row",
  fontsize = 12,
  fontsize_row = 10,
  angle_col = 45,
  main = "Van Galen AUCell scores per sample",
  color = colorRampPalette(c("Darkblue", "white", "red"))(100),
  border_color = NA
)
dev.off()

# ------------------------------------------------------------------
# 8. Statistical test: cancer vs normal per signature
# ------------------------------------------------------------------

stat_results <- list()
for (col in top50_cols) {
  sig_label <- gsub("AUCell_", "", col)

  cancer_scores <- seu@meta.data[[col]][seu$cancer_status == "Cancer"]
  normal_scores <- seu@meta.data[[col]][seu$cancer_status == "Normal"]

  if (length(cancer_scores) >= 10 && length(normal_scores) >= 10) {
    wt <- wilcox.test(cancer_scores, normal_scores, alternative = "greater")
    stat_results[[sig_label]] <- tibble(
      signature     = sig_label,
      n_cancer      = length(cancer_scores),
      n_normal      = length(normal_scores),
      mean_cancer   = mean(cancer_scores, na.rm = TRUE),
      mean_normal   = mean(normal_scores, na.rm = TRUE),
      fold_change   = mean(cancer_scores, na.rm = TRUE) / mean(normal_scores, na.rm = TRUE),
      wilcox_pvalue = wt$p.value
    )
  }
}

if (length(stat_results) > 0) {
  stat_df <- bind_rows(stat_results) %>%
    mutate(padj = p.adjust(wilcox_pvalue, method = "BH"))

  write.csv(stat_df, file.path(tab_dir, "aucell_cancer_vs_normal_wilcox.csv"), row.names = FALSE)

  cat("\nCancer vs Normal enrichment (Wilcoxon, one-sided):\n")
  print(as.data.frame(stat_df), row.names = FALSE)
} else {
  cat("\nNot enough cancer/normal cells for statistical testing\n")
}

# ------------------------------------------------------------------
# 9. Feature plots on UMAP (all cells)
# ------------------------------------------------------------------

aucell_cols <- grep("^AUCell_", colnames(seu@meta.data), value = TRUE)

for (reduction in intersect(c("umap", "umap_projected", "umap_harmony"), Reductions(seu))) {
  plots <- list()
  for (col in aucell_cols) {
    sig_label <- gsub("AUCell_", "", col)
    plots[[sig_label]] <- FeaturePlot(
      seu,
      features = col,
      reduction = reduction,
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) + ggtitle(sig_label)
  }

  p_combined <- wrap_plots(plots, ncol = 2)
  ggsave(
    file.path(fig_dir, paste0("04_featureplot_aucell_", reduction, ".pdf")),
    p_combined, width = 12, height = 5 * ceiling(length(plots) / 2)
  )
}

# ------------------------------------------------------------------
# 9b. Per-patient UMAPs with all signatures
# ------------------------------------------------------------------

sample_col <- if ("sample_name" %in% colnames(seu@meta.data)) "sample_name" else "orig.ident"
patients <- sort(unique(seu@meta.data[[sample_col]]))

for (reduction in intersect(c("umap", "umap_projected", "umap_harmony"), Reductions(seu))) {
  for (patient in patients) {
    safe_patient <- gsub("[^A-Za-z0-9_-]", "_", patient)
    patient_cells <- colnames(seu)[seu@meta.data[[sample_col]] == patient]

    if (length(patient_cells) < 50) {
      message("Skipping ", patient, ": too few cells (", length(patient_cells), ")")
      next
    }

    seu_sub <- subset(seu, cells = patient_cells)

    plots <- list()

    p_status <- DimPlot(
      seu_sub,
      reduction = reduction,
      group.by = "cancer_status",
      raster = FALSE
    ) +
      scale_color_manual(values = c("Cancer" = "#D95F02", "Normal" = "#1F78B4", "Not_called" = "grey70")) +
      ggtitle(paste(patient, "- Cancer status"))
    plots[["cancer_status"]] <- p_status

    celltype_col_use <- if ("predicted_CellType_Broad" %in% colnames(seu_sub@meta.data)) {
      "predicted_CellType_Broad"
    } else if ("predicted_CellType" %in% colnames(seu_sub@meta.data)) {
      "predicted_CellType"
    } else {
      NULL
    }

    if (!is.null(celltype_col_use)) {
      p_ct <- DimPlot(
        seu_sub,
        reduction = reduction,
        group.by = celltype_col_use,
        raster = FALSE
      ) + ggtitle(paste(patient, "- Cell type"))
      plots[["celltype"]] <- p_ct
    }

    for (col in aucell_cols) {
      sig_label <- gsub("AUCell_", "", col)
      plots[[sig_label]] <- FeaturePlot(
        seu_sub,
        features = col,
        reduction = reduction,
        order = TRUE,
        min.cutoff = "q05",
        max.cutoff = "q95"
      ) + ggtitle(sig_label)
    }

    n_plots <- length(plots)
    p_combined <- wrap_plots(plots, ncol = 2)
    ggsave(
      file.path(fig_dir, paste0("05_per_patient_", safe_patient, "_", reduction, ".pdf")),
      p_combined, width = 12, height = 5 * ceiling(n_plots / 2)
    )

    message("Per-patient UMAP done for ", patient, " (", length(patient_cells), " cells, ", n_plots, " panels)")
  }
}

# ------------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------------

saveRDS(seu, file.path(rds_dir, "seurat_van_galen_aucell.rds"))
saveRDS(cells_AUC, file.path(rds_dir, "aucell_results.rds"))

cat("\nOutputs saved to:", out_dir, "\n")
cat("Done.\n")
