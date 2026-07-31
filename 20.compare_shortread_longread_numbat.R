#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(purrr)
  library(readr)
  library(stringr)
})

PROJECT_DIR <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
RUN <- "260528_VH01624_464_222K7VKNX"

SR_BASE <- file.path(PROJECT_DIR, "results/seurat_annotated", RUN, "numbat")
LR_BASE <- file.path(PROJECT_DIR, "results/longread_numbat")

SEURAT_FILE <- file.path(
  PROJECT_DIR,
  "results/seurat_annotated",
  RUN,
  "demux_singlets_annotated_seurat.rds"
)

OUT_FIG <- file.path(PROJECT_DIR, "results/longread_numbat/comparison_figures")
OUT_TBL <- file.path(PROJECT_DIR, "results/longread_numbat/comparison_tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)

donor_map <- c(
  "donor0" = "LK2_APOP576-TP53",
  "donor1" = "LK2_HBDN498-TP53",
  "donor2" = "LK2_HBDN376-TP53",
  "donor3" = "LK2_normal-01"
)

load_clone_post <- function(numbat_dir) {
  candidates <- list.files(numbat_dir, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) return(NULL)
  iterations <- as.integer(str_extract(basename(candidates), "\\d+"))
  latest <- candidates[which.max(iterations)]
  message("  Loading: ", latest)
  read_tsv(latest, show_col_types = FALSE)
}

load_segs <- function(numbat_dir) {
  candidates <- list.files(numbat_dir, pattern = "^segs_consensus_\\d+\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) return(NULL)
  iterations <- as.integer(str_extract(basename(candidates), "\\d+"))
  latest <- candidates[which.max(iterations)]
  read_tsv(latest, show_col_types = FALSE)
}

load_geno <- function(numbat_dir) {
  candidates <- list.files(numbat_dir, pattern = "^geno_\\d+\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) return(NULL)
  iterations <- as.integer(str_extract(basename(candidates), "\\d+"))
  latest <- candidates[which.max(iterations)]
  read_tsv(latest, show_col_types = FALSE)
}

# ---------------------------------------------------------------
# Load Seurat metadata for broad cell type composition
# ---------------------------------------------------------------

seu <- readRDS(SEURAT_FILE)
meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(cell = sub("-1$", "", cell))

composition_broad <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  dplyr::count(sample_name, predicted_CellType_Broad) %>%
  group_by(sample_name) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    lineage = case_when(
      predicted_CellType_Broad %in% c(
        "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
        "Cycling Progenitor", "EoBasoMast Precursor",
        "Megakaryocyte Precursor"
      ) ~ "Stem / progenitor",
      
      predicted_CellType_Broad %in% c(
        "Monocyte", "Pro-Monocyte", "cDC", "pDC"
      ) ~ "Myeloid / DC",
      
      predicted_CellType_Broad %in% c(
        "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
        "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
      ) ~ "Lymphoid",
      
      predicted_CellType_Broad %in% c(
        "Early Erythroid", "Late Erythroid"
      ) ~ "Erythroid",
      
      TRUE ~ "Other"
    )
  )

write_csv(
  composition_broad,
  file.path(OUT_TBL, "seurat_broad_celltype_composition_by_sample.csv")
)

# ---------------------------------------------------------------
# Compare donors
# ---------------------------------------------------------------

comparison_results <- list()
composition_results <- list()

for (lr_donor in names(donor_map)) {
  
  sr_donor <- donor_map[[lr_donor]]
  message("\n=== Comparing ", lr_donor, " vs ", sr_donor, " ===")
  
  lr_dir <- file.path(LR_BASE, lr_donor)
  sr_dir <- file.path(SR_BASE, sr_donor, "numbat_final")
  
  if (!dir.exists(lr_dir)) { message("  Long-read dir missing: ", lr_dir); next }
  if (!dir.exists(sr_dir)) { message("  Short-read dir missing: ", sr_dir); next }
  
  lr_clones <- load_clone_post(lr_dir)
  sr_clones <- load_clone_post(sr_dir)
  
  if (is.null(lr_clones) || is.null(sr_clones)) {
    message("  Skipping - missing clone_post files")
    next
  }
  
  lr_clones <- lr_clones %>%
    select(cell, any_of(c("clone_opt", "compartment", "compartment_opt", "p_cnv", "clone_post"))) %>%
    rename_with(~ paste0(.x, "_lr"), -cell)
  
  sr_clones <- sr_clones %>%
    mutate(cell = sub("-1$", "", cell)) %>%
    select(cell, any_of(c("clone_opt", "compartment", "compartment_opt", "p_cnv", "clone_post"))) %>%
    rename_with(~ paste0(.x, "_sr"), -cell)
  
  merged <- inner_join(lr_clones, sr_clones, by = "cell") %>%
    left_join(
      meta %>%
        select(
          cell,
          sample_name,
          mapping_error_QC,
          predicted_CellType_Broad
        ),
      by = "cell"
    ) %>%
    mutate(
      lineage = case_when(
        predicted_CellType_Broad %in% c(
          "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
          "Cycling Progenitor", "EoBasoMast Precursor",
          "Megakaryocyte Precursor"
        ) ~ "Stem / progenitor",
        
        predicted_CellType_Broad %in% c(
          "Monocyte", "Pro-Monocyte", "cDC", "pDC"
        ) ~ "Myeloid / DC",
        
        predicted_CellType_Broad %in% c(
          "Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
          "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell"
        ) ~ "Lymphoid",
        
        predicted_CellType_Broad %in% c(
          "Early Erythroid", "Late Erythroid"
        ) ~ "Erythroid",
        
        TRUE ~ "Other"
      )
    )
  
  message("  Cells in common: ", nrow(merged),
          " (LR: ", nrow(lr_clones), ", SR: ", nrow(sr_clones), ")")
  
  if (nrow(merged) == 0) next
  
  comparison_results[[sr_donor]] <- merged
  
  compartment_lr_col <- intersect(c("compartment_lr", "compartment_opt_lr"), names(merged))[1]
  compartment_sr_col <- intersect(c("compartment_sr", "compartment_opt_sr"), names(merged))[1]
  
  # -------------------------------------------------------------
  # Broad composition among overlapping cells
  # -------------------------------------------------------------
  
  if (!is.na(compartment_lr_col) && !is.na(compartment_sr_col)) {
    
    donor_comp <- merged %>%
      filter(mapping_error_QC == "Pass") %>%
      mutate(
        donor = sr_donor,
        sr_compartment = .data[[compartment_sr_col]],
        lr_compartment = .data[[compartment_lr_col]]
      )
    
    comp_broad_counts <- donor_comp %>%
      dplyr::count(
        donor,
        predicted_CellType_Broad,
        lineage,
        sr_compartment,
        lr_compartment,
        name = "n"
      ) %>%
      group_by(donor) %>%
      mutate(percent_total_overlap = n / sum(n) * 100) %>%
      ungroup()
    
    composition_results[[sr_donor]] <- comp_broad_counts
    
    write_csv(
      comp_broad_counts,
      file.path(OUT_TBL, paste0(sr_donor, "_broad_composition_by_SR_LR_compartment.csv"))
    )
    
    # Broad cell type composition by Numbat compartment and method
    comp_broad_long <- donor_comp %>%
      select(cell, donor, predicted_CellType_Broad, lineage, sr_compartment, lr_compartment) %>%
      pivot_longer(
        cols = c(sr_compartment, lr_compartment),
        names_to = "method",
        values_to = "compartment"
      ) %>%
      mutate(
        method = recode(
          method,
          sr_compartment = "Short-read",
          lr_compartment = "Long-read"
        )
      ) %>%
      dplyr::count(donor, method, compartment, lineage, predicted_CellType_Broad, name = "n") %>%
      group_by(donor, method, compartment) %>%
      mutate(percent = n / sum(n) * 100) %>%
      ungroup()
    
    write_csv(
      comp_broad_long,
      file.path(OUT_TBL, paste0(sr_donor, "_broad_composition_long.csv"))
    )
    
    p_broad <- ggplot(
      comp_broad_long,
      aes(x = compartment, y = percent, fill = lineage)
    ) +
      geom_col(width = 0.7) +
      facet_wrap(~ method, nrow = 1) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(sr_donor, " — Broad lineage composition by Numbat compartment"),
        x = "Numbat compartment",
        y = "Percent of cells",
        fill = "Lineage"
      )
    
    ggsave(
      file.path(OUT_FIG, paste0(sr_donor, "_broad_lineage_composition_by_compartment.pdf")),
      p_broad,
      width = 9,
      height = 5
    )
    print(p_broad)
  }
  
  # -------------------------------------------------------------
  # Compartment alluvial
  # -------------------------------------------------------------
  
  if (!is.na(compartment_lr_col) && !is.na(compartment_sr_col)) {
    
    comp_counts <- merged %>%
      dplyr::count(
        short_read = !!sym(compartment_sr_col),
        long_read = !!sym(compartment_lr_col),
        name = "count"
      )
    
    p_comp <- ggplot(comp_counts, aes(y = count, axis1 = short_read, axis2 = long_read)) +
      geom_alluvium(aes(fill = short_read), width = 1/6, alpha = 0.7) +
      geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
      scale_x_discrete(limits = c("Short-read", "Long-read"), expand = c(0.15, 0.05)) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none") +
      labs(
        title = paste0(sr_donor, " — Compartment: short-read vs long-read"),
        y = "Number of cells"
      )
    
    ggsave(file.path(OUT_FIG, paste0(sr_donor, "_compartment_alluvial.pdf")),
           p_comp, width = 6, height = 6)
    print(p_comp)
    
    write_csv(
      comp_counts,
      file.path(OUT_TBL, paste0(sr_donor, "_compartment_concordance.csv"))
    )
    
    concordance <- mean(
      merged[[compartment_sr_col]] == merged[[compartment_lr_col]],
      na.rm = TRUE
    )
    message("  Compartment concordance: ", round(concordance * 100, 1), "%")
  }
  
  # -------------------------------------------------------------
  # Clone assignment alluvial
  # -------------------------------------------------------------
  
  clone_lr_col <- intersect(c("clone_opt_lr"), names(merged))[1]
  clone_sr_col <- intersect(c("clone_opt_sr"), names(merged))[1]
  
  if (!is.na(clone_lr_col) && !is.na(clone_sr_col)) {
    
    clone_counts <- merged %>%
      mutate(
        clone_sr = paste0("SR_", !!sym(clone_sr_col)),
        clone_lr = paste0("LR_", !!sym(clone_lr_col))
      ) %>%
      dplyr::count(clone_sr, clone_lr, name = "count")
    
    p_clone <- ggplot(clone_counts, aes(y = count, axis1 = clone_sr, axis2 = clone_lr)) +
      geom_alluvium(aes(fill = clone_sr), width = 1/6, alpha = 0.7) +
      geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(limits = c("Short-read", "Long-read"), expand = c(0.15, 0.05)) +
      scale_fill_brewer(palette = "Set3") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none") +
      labs(
        title = paste0(sr_donor, " — Clone: short-read vs long-read"),
        y = "Number of cells"
      )
    
    ggsave(file.path(OUT_FIG, paste0(sr_donor, "_clone_alluvial.pdf")),
           p_clone, width = 7, height = 6)
    print(p_clone)
    
    write_csv(
      clone_counts,
      file.path(OUT_TBL, paste0(sr_donor, "_clone_concordance.csv"))
    )
  }
  
  # -------------------------------------------------------------
  # p(CNV) scatter
  # -------------------------------------------------------------
  
  pcnv_lr_col <- intersect(c("p_cnv_lr"), names(merged))[1]
  pcnv_sr_col <- intersect(c("p_cnv_sr"), names(merged))[1]
  
  if (!is.na(pcnv_lr_col) && !is.na(pcnv_sr_col)) {
    
    cor_val <- cor(merged[[pcnv_sr_col]], merged[[pcnv_lr_col]], use = "complete.obs")
    
    p_scatter <- ggplot(merged, aes(x = !!sym(pcnv_sr_col), y = !!sym(pcnv_lr_col))) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      annotate("text", x = 0.1, y = 0.9, label = paste0("r = ", round(cor_val, 3)), size = 4) +
      theme_bw(base_size = 12) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        title = paste0(sr_donor, " — p(CNV) correlation"),
        x = "p(CNV) short-read",
        y = "p(CNV) long-read"
      )
    
    ggsave(file.path(OUT_FIG, paste0(sr_donor, "_pCNV_scatter.pdf")),
           p_scatter, width = 6, height = 6)
    print(p_scatter)
  }
}

# ---------------------------------------------------------------
# Summary table across donors
# ---------------------------------------------------------------

if (length(comparison_results) > 0) {
  
  summary_df <- imap_dfr(comparison_results, function(merged, donor) {
    
    compartment_lr_col <- intersect(c("compartment_lr", "compartment_opt_lr"), names(merged))[1]
    compartment_sr_col <- intersect(c("compartment_sr", "compartment_opt_sr"), names(merged))[1]
    pcnv_lr_col <- intersect(c("p_cnv_lr"), names(merged))[1]
    pcnv_sr_col <- intersect(c("p_cnv_sr"), names(merged))[1]
    
    n_total <- nrow(merged)
    
    comp_conc <- NA_real_
    n_sr_tumor <- n_sr_normal <- n_lr_tumor <- n_lr_normal <- NA_integer_
    n_agree_tumor <- n_agree_normal <- NA_integer_
    n_sr_tumor_lr_normal <- n_sr_normal_lr_tumor <- NA_integer_
    
    if (!is.na(compartment_lr_col) && !is.na(compartment_sr_col)) {
      
      sr <- merged[[compartment_sr_col]]
      lr <- merged[[compartment_lr_col]]
      
      comp_conc <- mean(sr == lr, na.rm = TRUE)
      
      n_sr_tumor  <- sum(sr == "tumor", na.rm = TRUE)
      n_sr_normal <- sum(sr == "normal", na.rm = TRUE)
      
      n_lr_tumor  <- sum(lr == "tumor", na.rm = TRUE)
      n_lr_normal <- sum(lr == "normal", na.rm = TRUE)
      
      n_agree_tumor <- sum(sr == "tumor" & lr == "tumor", na.rm = TRUE)
      n_agree_normal <- sum(sr == "normal" & lr == "normal", na.rm = TRUE)
      
      n_sr_tumor_lr_normal <- sum(sr == "tumor" & lr == "normal", na.rm = TRUE)
      n_sr_normal_lr_tumor <- sum(sr == "normal" & lr == "tumor", na.rm = TRUE)
    }
    
    pcnv_cor <- NA_real_
    if (!is.na(pcnv_lr_col) && !is.na(pcnv_sr_col)) {
      pcnv_cor <- cor(merged[[pcnv_sr_col]], merged[[pcnv_lr_col]], use = "complete.obs")
    }
    
    tibble(
      donor = donor,
      n_cells_overlap = n_total,
      
      n_sr_total = n_total,
      n_lr_total = n_total,
      
      n_sr_tumor = n_sr_tumor,
      n_sr_normal = n_sr_normal,
      n_lr_tumor = n_lr_tumor,
      n_lr_normal = n_lr_normal,
      
      n_agree_tumor = n_agree_tumor,
      n_agree_normal = n_agree_normal,
      n_sr_tumor_lr_normal = n_sr_tumor_lr_normal,
      n_sr_normal_lr_tumor = n_sr_normal_lr_tumor,
      
      compartment_concordance = comp_conc,
      pcnv_correlation = pcnv_cor
    )
  })
  
  write_csv(summary_df, file.path(OUT_TBL, "comparison_summary.csv"))
  print(summary_df)
  
  # -------------------------------------------------------------
  # Summary plot: total, tumour, normal counts
  # -------------------------------------------------------------
  
  cell_counts_long <- summary_df %>%
    select(
      donor,
      n_sr_total, n_lr_total,
      n_sr_tumor, n_lr_tumor,
      n_sr_normal, n_lr_normal
    ) %>%
    pivot_longer(
      cols = -donor,
      names_to = "metric",
      values_to = "n_cells"
    ) %>%
    mutate(
      method = case_when(
        str_detect(metric, "_sr_") ~ "Short-read",
        str_detect(metric, "_lr_") ~ "Long-read",
        TRUE ~ NA_character_
      ),
      compartment = case_when(
        str_detect(metric, "total") ~ "Total",
        str_detect(metric, "tumor") ~ "Tumor",
        str_detect(metric, "normal") ~ "Normal",
        TRUE ~ NA_character_
      ),
      compartment = factor(compartment, levels = c("Total", "Tumor", "Normal"))
    )
  
  p_cell_counts <- ggplot(
    cell_counts_long,
    aes(x = donor, y = n_cells, fill = method)
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(
      aes(label = n_cells),
      position = position_dodge(width = 0.7),
      vjust = -0.4,
      size = 3
    ) +
    facet_wrap(~ compartment, scales = "free_y", nrow = 1) +
    theme_bw(base_size = 12) +
    labs(
      title = "Total, tumour, and normal cell counts: short-read vs long-read",
      x = NULL,
      y = "Number of overlapping cells",
      fill = NULL
    )
  
  ggsave(
    file.path(OUT_FIG, "summary_total_tumor_normal_counts.pdf"),
    p_cell_counts,
    width = 12,
    height = 5
  )
  print(p_cell_counts)
  
  # -------------------------------------------------------------
  # Summary broad composition across donors
  # -------------------------------------------------------------
  
  all_broad <- bind_rows(composition_results)
  
  if (nrow(all_broad) > 0) {
    
    write_csv(
      all_broad,
      file.path(OUT_TBL, "summary_broad_composition_by_SR_LR_compartment.csv")
    )
    
    broad_lineage_summary <- all_broad %>%
      group_by(donor, lineage) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      group_by(donor) %>%
      mutate(percent = n / sum(n) * 100) %>%
      ungroup()
    
    p_lineage_summary <- ggplot(
      broad_lineage_summary,
      aes(x = donor, y = percent, fill = lineage)
    ) +
      geom_col(width = 0.7) +
      theme_bw(base_size = 12) +
      labs(
        title = "Broad lineage composition of overlapping cells",
        x = NULL,
        y = "Percent of overlapping cells",
        fill = "Lineage"
      )
    
    ggsave(
      file.path(OUT_FIG, "summary_broad_lineage_composition_overlap.pdf"),
      p_lineage_summary,
      width = 8,
      height = 5
    )
    print(p_lineage_summary)
  }
}

message("\nDone. Figures: ", OUT_FIG, " | Tables: ", OUT_TBL)

