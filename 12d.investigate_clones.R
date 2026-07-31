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
  paste0("Clone_", clone_to_check) = "#E41A1C",
  "Other_clones" = "grey75",
  "No_clone" = "grey90"
)

status_cols <- c(
  paste0("Clone_", clone_to_check, "_lymphoid") = "#E41A1C",
  paste0("Clone_", clone_to_check, "_non_lymphoid") = "#FF7F00",
  "Other_lymphoid" = "#377EB8",
  "Other_non_lymphoid" = "grey85"
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