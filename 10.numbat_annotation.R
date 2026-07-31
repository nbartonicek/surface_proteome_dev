# ============================================================
# Integrate available per-donor Numbat outputs into Seurat
# Handles incomplete Numbat runs gracefully
# ============================================================

library(Seurat)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(stringr)
library(purrr)

# ----------------------------
# Paths
# ----------------------------
run <- "260528_VH01624_464_222K7VKNX"
proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
sample_short <- "LK2"

out_dir <- file.path(
  proj,
  "results/seurat_annotated",
  run,
  "numbat"
)

annotation_dir <- file.path(
  proj,
  "results/seurat_annotated",
  run
)

seurat_file <- file.path(
  annotation_dir,
  "demux_singlets_annotated_seurat.rds"
)

seu <- readRDS(seurat_file)

# ----------------------------
# Donor Numbat directories
# ----------------------------
samples <- unique(seu$sample_name)

numbat_dirs <- setNames(
  file.path(
    out_dir,
    paste0(sample_short, "_", samples),
    "numbat_final"
  ),
  samples
)

numbat_dirs <- numbat_dirs[dir.exists(numbat_dirs)]

cat("\nDetected Numbat dirs:\n")
print(numbat_dirs)

# ----------------------------
# Helper: inspect Numbat output status
# ----------------------------
inspect_numbat_dir <- function(numbat_dir, donor_name) {
  
  tibble(
    donor = donor_name,
    numbat_dir = numbat_dir,
    has_clone_post = file.exists(file.path(numbat_dir, "clone_post_1.tsv")),
    has_sc_refs = file.exists(file.path(numbat_dir, "sc_refs.rds")),
    has_segs_consensus = file.exists(file.path(numbat_dir, "segs_consensus_1.tsv")),
    has_segs_loh = file.exists(file.path(numbat_dir, "segs_loh.tsv")),
    has_bulk_clones = file.exists(file.path(numbat_dir, "bulk_clones_1.tsv.gz")),
    has_log = file.exists(file.path(numbat_dir, "log.txt"))
  )
}

numbat_status <- bind_rows(
  imap(numbat_dirs, inspect_numbat_dir)
)

write.csv(
  numbat_status,
  file.path(out_dir, paste0(sample_short, "_numbat_output_status.csv")),
  row.names = FALSE
)

cat("\nNumbat output status:\n")
print(numbat_status)

# ----------------------------
# Helper: safely extract cell names from sc_refs.rds
# ----------------------------
extract_sc_ref_cells <- function(sc_refs_file) {
  
  obj <- readRDS(sc_refs_file)
  
  possible_cells <- character()
  
  if (is.character(obj)) {
    possible_cells <- obj
  } else if (is.data.frame(obj) || is.matrix(obj)) {
    possible_cells <- rownames(obj)
  } else if (is.list(obj)) {
    
    list_cells <- unlist(
      lapply(obj, function(x) {
        if (is.character(x)) return(x)
        if (is.data.frame(x) || is.matrix(x)) return(rownames(x))
        return(NULL)
      }),
      use.names = FALSE
    )
    
    possible_cells <- list_cells
  }
  
  possible_cells <- unique(as.character(possible_cells))
  possible_cells <- possible_cells[!is.na(possible_cells) & possible_cells != ""]
  
  possible_cells
}

# ----------------------------
# Read best available Numbat calls
# ----------------------------
read_numbat_available <- function(numbat_dir, donor_name) {
  
  message("\nReading available Numbat output for: ", donor_name)
  
  clone_file <- file.path(numbat_dir, "clone_post_1.tsv")
  sc_refs_file <- file.path(numbat_dir, "sc_refs.rds")
  seg_file <- file.path(numbat_dir, "segs_consensus_1.tsv")
  loh_file <- file.path(numbat_dir, "segs_loh.tsv")
  
  # Case 1: complete per-cell Numbat clone calls
  if (file.exists(clone_file)) {
    
    x <- fread(clone_file) %>%
      as_tibble()
    
    required_cols <- c(
      "cell", "clone_opt", "GT_opt",
      "p_opt", "p_cnv", "compartment_opt"
    )
    
    missing_cols <- setdiff(required_cols, colnames(x))
    
    if (length(missing_cols) > 0) {
      warning(
        "clone_post_1.tsv exists for ", donor_name,
        " but is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
      return(NULL)
    }
    
    return(
      x %>%
        transmute(
          cell = as.character(cell),
          numbat_donor = donor_name,
          numbat_source = "clone_post_1.tsv",
          numbat_status = "complete_per_cell_calls",
          numbat_clone = as.character(clone_opt),
          numbat_GT = as.character(GT_opt),
          numbat_p_opt = as.numeric(p_opt),
          numbat_p_cnv = as.numeric(p_cnv),
          numbat_compartment = case_when(
            compartment_opt == "tumor" ~ "tumor",
            compartment_opt == "normal" ~ "normal",
            TRUE ~ "not_called"
          ),
          numbat_call = numbat_compartment
        )
    )
  }
  
  # Case 2: no clone_post, but reference cells exist
  if (file.exists(sc_refs_file)) {
    
    ref_cells <- extract_sc_ref_cells(sc_refs_file)
    
    if (length(ref_cells) > 0) {
      
      return(
        tibble(
          cell = ref_cells,
          numbat_donor = donor_name,
          numbat_source = "sc_refs.rds",
          numbat_status = "reference_cells_only",
          numbat_clone = "reference",
          numbat_GT = NA_character_,
          numbat_p_opt = NA_real_,
          numbat_p_cnv = NA_real_,
          numbat_compartment = "normal",
          numbat_call = "normal"
        )
      )
    }
  }
  
  # Case 3: only sample-level segment outputs exist
  if (file.exists(seg_file) || file.exists(loh_file)) {
    
    warning(
      donor_name,
      " has segment-level Numbat output but no per-cell clone_post_1.tsv. ",
      "Will not add per-cell calls for this donor."
    )
    
    return(
      tibble(
        cell = character(),
        numbat_donor = donor_name,
        numbat_source = "segments_only",
        numbat_status = "segments_only_no_per_cell_calls",
        numbat_clone = character(),
        numbat_GT = character(),
        numbat_p_opt = numeric(),
        numbat_p_cnv = numeric(),
        numbat_compartment = character(),
        numbat_call = character()
      )
    )
  }
  
  warning("No usable Numbat output found for ", donor_name)
  NULL
}

numbat_meta_list <- imap(numbat_dirs, read_numbat_available)
numbat_meta_list <- Filter(Negate(is.null), numbat_meta_list)

if (length(numbat_meta_list) == 0) {
  stop("No usable Numbat metadata found.")
}

numbat_meta <- bind_rows(numbat_meta_list)

# Remove empty segment-only rows if present
numbat_meta <- numbat_meta %>%
  filter(!is.na(cell), cell != "")

cat("\nNumbat metadata rows:", nrow(numbat_meta), "\n")
cat("Seurat cells:", ncol(seu), "\n")

# ----------------------------
# Barcode matching
# ----------------------------
direct_overlap <- sum(numbat_meta$cell %in% colnames(seu))
cat("Direct barcode overlap:", direct_overlap, "\n")

if (direct_overlap < 100 && nrow(numbat_meta) > 0) {
  
  message("Low direct overlap. Trying stripped barcode matching...")
  
  seu_barcode_stripped <- gsub("-1$", "", colnames(seu))
  numbat_barcode_stripped <- gsub("-1$", "", numbat_meta$cell)
  
  matched_idx <- match(numbat_barcode_stripped, seu_barcode_stripped)
  
  numbat_meta$cell_original_numbat <- numbat_meta$cell
  numbat_meta$cell <- colnames(seu)[matched_idx]
  
  numbat_meta <- numbat_meta %>%
    filter(!is.na(cell))
  
  cat(
    "Overlap after stripped matching:",
    sum(numbat_meta$cell %in% colnames(seu)),
    "\n"
  )
}

common_cells <- intersect(colnames(seu), numbat_meta$cell)

cat("\nNumbat cells matched to Seurat:", length(common_cells), "\n")

# ----------------------------
# Add metadata to Seurat
# ----------------------------
numbat_meta_sub <- numbat_meta %>%
  filter(cell %in% common_cells) %>%
  arrange(match(cell, colnames(seu))) %>%
  distinct(cell, .keep_all = TRUE)

meta_to_add <- numbat_meta_sub %>%
  column_to_rownames("cell")

old_numbat_cols <- grep("^numbat_", colnames(seu@meta.data), value = TRUE)

if (length(old_numbat_cols) > 0) {
  seu@meta.data[, old_numbat_cols] <- NULL
}

if (nrow(meta_to_add) > 0) {
  seu <- AddMetaData(seu, metadata = meta_to_add)
}

# Fill missing values
default_numbat_cols <- list(
  numbat_donor = "not_available",
  numbat_source = "not_available",
  numbat_status = "not_available",
  numbat_clone = "not_called",
  numbat_GT = NA_character_,
  numbat_p_opt = NA_real_,
  numbat_p_cnv = NA_real_,
  numbat_compartment = "not_called",
  numbat_call = "not_called"
)

for (cc in names(default_numbat_cols)) {
  
  if (!cc %in% colnames(seu@meta.data)) {
    seu@meta.data[[cc]] <- default_numbat_cols[[cc]]
  }
  
  if (is.character(default_numbat_cols[[cc]])) {
    seu@meta.data[[cc]][is.na(seu@meta.data[[cc]])] <- default_numbat_cols[[cc]]
  }
}

seu$numbat_compartment <- factor(
  seu$numbat_compartment,
  levels = c("tumor", "normal", "not_called")
)

seu$numbat_call <- factor(
  seu$numbat_call,
  levels = c("tumor", "normal", "not_called")
)

cat("\nNumbat status table:\n")
print(table(seu$numbat_status, useNA = "ifany"))

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

# ----------------------------
# Helper for DimPlots
# ----------------------------
save_dimplot <- function(filename, reduction, group_by, title, split_by = NULL,
                         cols = NULL, width = 8, height = 6, ncol = 2) {
  
  pdf(file.path(out_dir, filename), width = width, height = height)
  
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
# UMAP plots
# ----------------------------
for (red_name in names(plot_reductions)) {
  
  red <- plot_reductions[[red_name]]
  
  save_dimplot(
    filename = paste0(sample_short, "_09_numbat_compartment_on_", red, ".pdf"),
    reduction = red,
    group_by = "numbat_compartment",
    cols = numbat_cols,
    title = paste0("Numbat compartment on ", red),
    width = 8,
    height = 6
  )
  
  save_dimplot(
    filename = paste0(sample_short, "_10_numbat_compartment_on_", red, "_split_by_sample.pdf"),
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
    filename = paste0(sample_short, "_10b_numbat_clone_on_", red, "_split_by_sample.pdf"),
    reduction = red,
    group_by = "numbat_clone",
    split_by = split_col,
    title = paste0("Numbat clone assignments by sample on ", red),
    width = 14,
    height = 10,
    ncol = 2
  )
  
  save_dimplot(
    filename = paste0(sample_short, "_10c_numbat_status_on_", red, "_split_by_sample.pdf"),
    reduction = red,
    group_by = "numbat_status",
    split_by = split_col,
    title = paste0("Numbat output status by sample on ", red),
    width = 14,
    height = 10,
    ncol = 2
  )
}

# ----------------------------
# Composition by broad annotation
# ----------------------------
annotation_col <- if ("predicted_CellType_Broad" %in% colnames(seu@meta.data)) {
  "predicted_CellType_Broad"
} else if ("BoneMarrowMap_cell_type" %in% colnames(seu@meta.data)) {
  "BoneMarrowMap_cell_type"
} else {
  stop("No broad annotation column found.")
}

numbat_summary <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!is.na(.data[[split_col]])) %>%
  count(
    sample = .data[[split_col]],
    broad_annotation = .data[[annotation_col]],
    numbat_compartment,
    name = "n"
  ) %>%
  group_by(sample, broad_annotation) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  numbat_summary,
  file.path(out_dir, paste0(sample_short, "_numbat_compartment_by_sample_and_broad_annotation.csv")),
  row.names = FALSE
)

p_numbat_comp <- numbat_summary %>%
  ggplot(aes(x = broad_annotation, y = percent, fill = numbat_compartment)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ sample) +
  scale_fill_manual(values = numbat_cols, drop = FALSE) +
  theme_bw() +
  labs(
    x = "Broad annotation",
    y = "Cells (%)",
    fill = "Numbat compartment",
    title = "Numbat compartments within broad annotations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, paste0(sample_short, "_11_numbat_compartment_by_sample_and_broad_annotation.pdf")),
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
    broad_annotation = .data[[annotation_col]],
    numbat_clone,
    name = "n"
  ) %>%
  group_by(sample, broad_annotation) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(
  numbat_clone_summary,
  file.path(out_dir, paste0(sample_short, "_numbat_clone_by_sample_and_broad_annotation.csv")),
  row.names = FALSE
)

p_numbat_clone_comp <- numbat_clone_summary %>%
  ggplot(aes(x = broad_annotation, y = percent, fill = numbat_clone)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ sample) +
  theme_bw() +
  labs(
    x = "Broad annotation",
    y = "Cells (%)",
    fill = "Numbat clone",
    title = "Numbat clone/reference assignments within broad annotations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, paste0(sample_short, "_12_numbat_clone_by_sample_and_broad_annotation.pdf")),
  p_numbat_clone_comp,
  width = 12,
  height = 5
)

# ----------------------------
# Save integrated object
# ----------------------------
saveRDS(
  seu,
  file.path(out_dir, paste0(sample_short, "_projected_CITE_DSB_Numbat_integrated.rds"))
)

write.csv(
  seu@meta.data,
  file.path(out_dir, paste0(sample_short, "_projected_CITE_DSB_Numbat_integrated_metadata.csv"))
)

cat("\nDone. Integrated available Numbat metadata and saved plots to:\n")
cat(out_dir, "\n")