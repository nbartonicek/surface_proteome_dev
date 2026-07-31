# ============================================================
# Integrate per-donor Numbat calls into Seurat and plot
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(stringr)

# ----------------------------
# Paths
# ----------------------------
run <- "260423_VH01624_453_222HWMYNX"
proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

out_dir <- file.path(
  proj,
  "results/seurat_annotated",
  run,
  "numbat"
)

projection_dir <- out_dir

annotation_dir <- file.path(
  proj,
  "results/seurat_annotated",
  run
)

seurat_file <- file.path(
  annotation_dir,
  "LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds"
)

# Uncomment if needed
# seu <- readRDS(seurat_file)

# ----------------------------
# Donor Numbat directories
# ----------------------------
numbat_dirs <- c(
  "HBDN206-MNpCT" = file.path(out_dir, "LK1_HBDN206-MNpCT", "numbat_final"),
  "HBDN501-AML-KMT2A" = file.path(out_dir, "LK1_HBDN501-AML-KMT2A", "numbat_final")
  # "HBDN392-AML-MDS" = file.path(out_dir, "LK1_HBDN392-AML-MDS", "numbat_final")
)

# ----------------------------
# Read Numbat clone_post_1.tsv
# ----------------------------
read_numbat_calls <- function(numbat_dir, donor_name) {
  
  message("\nReading Numbat output for: ", donor_name)
  
  cell_file <- file.path(numbat_dir, "clone_post_1.tsv")
  
  if (!file.exists(cell_file)) {
    stop("Missing file: ", cell_file)
  }
  
  x <- fread(cell_file) %>%
    as_tibble()
  
  required_cols <- c("cell", "clone_opt", "GT_opt", "p_opt", "p_cnv", "compartment_opt")
  missing_cols <- setdiff(required_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing expected columns in clone_post_1.tsv for ",
      donor_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  x %>%
    transmute(
      cell = as.character(cell),
      numbat_donor = donor_name,
      numbat_clone = as.character(clone_opt),
      numbat_GT = as.character(GT_opt),
      numbat_p_opt = as.numeric(p_opt),
      numbat_p_cnv = as.numeric(p_cnv),
      numbat_compartment = case_when(
        compartment_opt == "tumor" ~ "tumor",
        compartment_opt == "normal" ~ "normal",
        TRUE ~ "not_called"
      ),
      numbat_call = case_when(
        compartment_opt == "tumor" ~ "tumor",
        compartment_opt == "normal" ~ "normal",
        TRUE ~ "not_called"
      )
    )
}

numbat_meta <- bind_rows(
  lapply(names(numbat_dirs), function(donor) {
    read_numbat_calls(numbat_dirs[[donor]], donor)
  })
)

# ----------------------------
# Check barcode overlap
# ----------------------------
cat("\nNumbat rows:", nrow(numbat_meta), "\n")
cat("Seurat cells:", ncol(seu), "\n")

direct_overlap <- sum(numbat_meta$cell %in% colnames(seu))
cat("Direct barcode overlap:", direct_overlap, "\n")

# If direct overlap is poor, try mapping by stripped barcode
if (direct_overlap < 100) {
  
  message("Low direct overlap. Trying stripped barcode matching...")
  
  seu_barcode_stripped <- gsub("-1$", "", colnames(seu))
  numbat_barcode_stripped <- gsub("-1$", "", numbat_meta$cell)
  
  matched_idx <- match(numbat_barcode_stripped, seu_barcode_stripped)
  
  numbat_meta$cell_original_numbat <- numbat_meta$cell
  numbat_meta$cell <- colnames(seu)[matched_idx]
  
  numbat_meta <- numbat_meta %>%
    filter(!is.na(cell))
  
  cat("Overlap after stripped matching:", sum(numbat_meta$cell %in% colnames(seu)), "\n")
}

common_cells <- intersect(colnames(seu), numbat_meta$cell)

cat("\nNumbat cells matched to Seurat:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  stop("No Numbat cells matched to Seurat. Check barcode format/sample filtering.")
}

# ----------------------------
# Add metadata to Seurat
# ----------------------------
numbat_meta_sub <- numbat_meta %>%
  filter(cell %in% common_cells) %>%
  distinct(cell, .keep_all = TRUE) %>%
  arrange(match(cell, colnames(seu)))

meta_to_add <- numbat_meta_sub %>%
  column_to_rownames("cell")

# Remove old Numbat columns if re-running
old_numbat_cols <- grep("^numbat_", colnames(seu@meta.data), value = TRUE)

if (length(old_numbat_cols) > 0) {
  seu@meta.data[, old_numbat_cols] <- NULL
}

seu <- AddMetaData(seu, metadata = meta_to_add)

# Fill missing cells
for (cc in c("numbat_compartment", "numbat_call", "numbat_clone")) {
  if (!cc %in% colnames(seu@meta.data)) {
    seu[[cc]] <- "not_called"
  }
  seu@meta.data[[cc]][is.na(seu@meta.data[[cc]])] <- "not_called"
}

seu$numbat_compartment <- factor(
  seu$numbat_compartment,
  levels = c("tumor", "normal", "not_called")
)

seu$numbat_call <- factor(
  seu$numbat_call,
  levels = c("tumor", "normal", "not_called")
)

cat("\nNumbat compartment table:\n")
print(table(seu$numbat_compartment, useNA = "ifany"))

cat("\nNumbat clone table:\n")
print(table(seu$numbat_clone, useNA = "ifany"))

# ----------------------------
# Plot object
# ----------------------------
seu_numbat_plot <- if ("mapping_error_QC" %in% colnames(seu@meta.data)) {
  subset(seu, subset = mapping_error_QC == "Pass")
} else {
  seu
}

numbat_cols <- c(
  "tumor" = "red",
  "normal" = "grey70",
  "not_called" = "black"
)

split_col <- if ("sampleID" %in% colnames(seu@meta.data)) {
  "sampleID"
} else if ("sample_name" %in% colnames(seu@meta.data)) {
  "sample_name"
} else {
  stop("No sampleID or sample_name column found.")
}

# ----------------------------
# Detect reductions
# ----------------------------
available_reductions <- Reductions(seu_numbat_plot)

cat("\nAvailable reductions:\n")
print(available_reductions)

standard_reduction <- if ("umap" %in% available_reductions) "umap" else NA_character_

projected_reduction <- case_when(
  #"harmony_projected" %in% available_reductions ~ "harmony_projected",
  "umap_projected" %in% available_reductions ~ "umap_projected",
  TRUE ~ NA_character_
)

plot_reductions <- c(
  standard = standard_reduction,
  projected = projected_reduction
)

plot_reductions <- plot_reductions[!is.na(plot_reductions)]

if (length(plot_reductions) == 0) {
  stop(
    "No suitable UMAP reduction found. Available reductions: ",
    paste(available_reductions, collapse = ", ")
  )
}

cat("\nPlotting reductions:\n")
print(plot_reductions)

# ----------------------------
# Helper for safe DimPlots
# ----------------------------
save_dimplot <- function(filename, reduction, group_by, title, split_by = NULL,
                         cols = NULL, width = 8, height = 6, ncol = 2) {
  
  pdf(file.path(projection_dir, filename), width = width, height = height)
  
  p <- DimPlot(
    seu_numbat_plot,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    cols = cols,
    ncol = ncol
  ) +
    ggtitle(title)
  
  print(p)
  dev.off()
}

# ----------------------------
# UMAP plots: compartment
# ----------------------------
for (red_name in names(plot_reductions)) {
  
  red <- plot_reductions[[red_name]]
  
  save_dimplot(
    filename = paste0("09_numbat_compartment_on_", red, ".pdf"),
    reduction = red,
    group_by = "numbat_compartment",
    cols = numbat_cols,
    title = paste0("Numbat compartment on ", red),
    width = 8,
    height = 6
  )
  
  save_dimplot(
    filename = paste0("10_numbat_compartment_on_", red, "_split_by_sample.pdf"),
    reduction = red,
    group_by = "numbat_compartment",
    split_by = split_col,
    cols = numbat_cols,
    title = paste0("Numbat compartment by sample on ", red),
    width = 14,
    height = 10,
    ncol = 2
  )
  
  save_dimplot(
    filename = paste0("10b_numbat_clone_on_", red, "_split_by_sample.pdf"),
    reduction = red,
    group_by = "numbat_clone",
    split_by = split_col,
    title = paste0("Numbat clone assignments by sample on ", red),
    width = 14,
    height = 10,
    ncol = 2
  )
}

# ----------------------------
# Composition by broad annotation
# ----------------------------
numbat_summary <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(.data[[split_col]])) %>%
  count(
    sample = .data[[split_col]],
    predicted_CellType_Broad,
    numbat_compartment,
    name = "n"
  ) %>%
  group_by(sample, predicted_CellType_Broad) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  numbat_summary,
  file.path(out_dir, "numbat_compartment_by_sample_and_broad_annotation.csv"),
  row.names = FALSE
)

p_numbat_comp <- numbat_summary %>%
  ggplot(aes(x = predicted_CellType_Broad, y = percent, fill = numbat_compartment)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ sample) +
  scale_fill_manual(values = numbat_cols, drop = FALSE) +
  theme_bw() +
  labs(
    x = "Broad annotation",
    y = "Cells (%)",
    fill = "Numbat compartment",
    title = "Numbat compartments within broad BoneMarrowMap annotations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "11_numbat_compartment_by_sample_and_broad_annotation.pdf"),
  p_numbat_comp,
  width = 12,
  height = 5
)

# ----------------------------
# Clone composition
# ----------------------------
numbat_clone_summary <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(.data[[split_col]])) %>%
  count(
    sample = .data[[split_col]],
    predicted_CellType_Broad,
    numbat_clone,
    name = "n"
  ) %>%
  group_by(sample, predicted_CellType_Broad) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  numbat_clone_summary,
  file.path(out_dir, "numbat_clone_by_sample_and_broad_annotation.csv"),
  row.names = FALSE
)

p_numbat_clone_comp <- numbat_clone_summary %>%
  ggplot(aes(x = predicted_CellType_Broad, y = percent, fill = numbat_clone)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ sample) +
  theme_bw() +
  labs(
    x = "Broad annotation",
    y = "Cells (%)",
    fill = "Numbat clone",
    title = "Numbat clone assignments within broad BoneMarrowMap annotations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "12_numbat_clone_by_sample_and_broad_annotation.pdf"),
  p_numbat_clone_comp,
  width = 12,
  height = 5
)

# ----------------------------
# Save integrated object
# ----------------------------
saveRDS(
  seu,
  file.path(out_dir, "LK1_projected_CITE_DSB_Numbat_integrated.rds")
)

write.csv(
  seu@meta.data,
  file.path(out_dir, "LK1_projected_CITE_DSB_Numbat_integrated_metadata.csv")
)

cat("\nDone. Integrated Numbat metadata and saved plots to:\n")
cat(out_dir, "\n")