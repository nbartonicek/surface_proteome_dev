#!/usr/bin/env Rscript

# Cross-method, cross-modality clone comparison.
#
# Brings together the four clone/malignancy callers used in this project:
#   - numbat      (CNV-based, run separately per donor, short-read + long-read)
#   - copykat     (CNV-based binary aneuploid/diploid call, short-read only -
#                  no long-read copykat script exists in this project, since
#                  FLAMES long-read data hasn't been run through it)
#   - mgatk       (mitochondrial variant clusters, via 15.mgatk.R / 15a.mgatk_longread.sh)
#   - mitoClone2  (mitochondrial variant clusters, via 15b.mitoclone2_calling.R)
#
# Two comparisons are produced per donor, mirroring the existing
# 20.compare_shortread_longread_numbat.R pattern (donor_map, per-donor
# directories, comparison_figures/comparison_tables outputs):
#
#   A) Within-modality, cross-method concordance: for short-read and for
#      long-read separately, how well do numbat/copykat/mgatk/mitoClone2
#      agree on clone/clustering structure for the same cells? Measured
#      with the adjusted Rand index (ARI), since these methods don't share
#      a common label naming scheme.
#
#   B) Cross-modality concordance for the two genetic methods (mgatk,
#      mitoClone2): do short-read and long-read data give the same clone
#      calls for the same cells? Same treatment numbat already gets in
#      20.compare_shortread_longread_numbat.R.
#
# Expects the outputs of 15.mgatk.R, 15a.mgatk_longread.sh + 15.mgatk.R,
# 15b.mitoclone2_calling.R, and the existing numbat/copykat pipelines to
# already exist on disk. Any missing input is skipped with a message rather
# than erroring, since not everything may have been run for every
# donor/modality yet.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

PROJECT_DIR <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
RUN <- "260528_VH01624_464_222K7VKNX"

SR_SAMPLE <- "LK2-GEX"   # short-read pooled GEM well (cellranger sample name)
LR_SAMPLE <- "LK2"       # long-read pooled GEM well (results/long_read/<LR_SAMPLE>/bam)

# Pipeline-native numbat output (results_nf/<run>/09c_numbat_run/) - confirmed
# complete (has clone_post_1.tsv) for all 4 donors including normal-01, unlike
# the results/seurat_annotated/<run>/numbat/ copy this used to point to, which
# is missing clone_post and most other converged-run files for normal-01
# specifically (a stale/partial copy, not the source of truth).
SR_NUMBAT_BASE <- file.path(PROJECT_DIR, "results_nf", RUN, "09c_numbat_run")
LR_NUMBAT_BASE <- file.path(PROJECT_DIR, "results/longread_numbat")

# Per-donor copykat prediction tables (the "temp_with_all_copykat_calls_
# metadata.csv" pre-merged file this used to point to doesn't exist for
# this run - copykat is run and written out per-donor, not pooled).
SR_COPYKAT_BASE <- file.path(PROJECT_DIR, "results_nf", RUN, "09_copykat")

SR_MGATK_CSV <- file.path(
  PROJECT_DIR, "results/mitochondrial_clones", RUN, SR_SAMPLE,
  paste0(SR_SAMPLE, ".mgatk_clones.csv")
)
LR_MGATK_CSV <- file.path(
  PROJECT_DIR, "results/mitochondrial_clones_longread", LR_SAMPLE,
  paste0(LR_SAMPLE, ".mgatk_clones.csv")
)

SR_MITOCLONE2_CSV <- file.path(
  PROJECT_DIR, "results/mitochondrial_clones", RUN, SR_SAMPLE, "mitoclone2",
  paste0(SR_SAMPLE, ".mitoclone2_clones.csv")
)
LR_MITOCLONE2_CSV <- file.path(
  PROJECT_DIR, "results/mitochondrial_clones_longread", LR_SAMPLE, "mitoclone2",
  paste0(LR_SAMPLE, ".mitoclone2_clones.csv")
)

SEURAT_FILE <- file.path(
  PROJECT_DIR, "results/seurat_annotated", RUN,
  "demux_singlets_annotated_seurat.rds"
)

OUT_FIG <- file.path(PROJECT_DIR, "results/mitochondrial_clones/comparison_all_methods_figures")
OUT_TBL <- file.path(PROJECT_DIR, "results/mitochondrial_clones/comparison_all_methods_tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)

donor_map <- c(
  "donor0" = "LK2_APOP576-TP53",
  "donor1" = "LK2_HBDN498-TP53",
  "donor2" = "LK2_HBDN376-TP53",
  "donor3" = "LK2_normal-01"
)

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------

# Standard Hubert & Arabie (1985) adjusted Rand index. Used instead of a
# simple percent-agreement because the clone labels across methods (e.g.
# numbat's clone_opt vs mgatk's kmeans cluster ids) have no shared naming
# scheme - only their partition structure is comparable.
adjusted_rand_index <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  tab <- table(x, y)
  n <- sum(tab)
  a <- sum(choose(rowSums(tab), 2))
  b <- sum(choose(colSums(tab), 2))
  c_sum <- sum(choose(tab, 2))
  expected <- a * b / choose(n, 2)
  max_index <- (a + b) / 2
  if (max_index == expected) return(NA_real_)
  (c_sum - expected) / (max_index - expected)
}

load_clone_post <- function(numbat_dir) {
  candidates <- list.files(numbat_dir, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) return(NULL)
  iterations <- as.integer(str_extract(basename(candidates), "\\d+"))
  latest <- candidates[which.max(iterations)]
  message("  Loading: ", latest)
  read_tsv(latest, show_col_types = FALSE)
}

load_csv_if_exists <- function(path, label) {
  if (!file.exists(path)) {
    message("  ", label, " not found, skipping: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

# bare_donor: donor_id with no sample prefix (e.g. "APOP576-TP53"), matching
# both this file's naming (copykat_prediction_with_metadata.csv lives under
# 09_copykat/<bare_donor>/tables/) and the sample_name values used in the
# Seurat metadata / composition tables throughout the rest of this project.
load_sr_copykat_donor <- function(bare_donor) {
  path <- file.path(SR_COPYKAT_BASE, bare_donor, "tables", "copykat_prediction_with_metadata.csv")
  if (!file.exists(path)) {
    message("  SR copykat not found for ", bare_donor, ", skipping: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) %>%
    mutate(cell = sub("-1$", "", cell)) %>%
    dplyr::select(cell, copykat = copykat_call)
}

# ---------------------------------------------------------------
# Load pooled per-well method outputs once (filtered per-donor below)
# ---------------------------------------------------------------

sr_mgatk_all      <- load_csv_if_exists(SR_MGATK_CSV, "SR mgatk")
lr_mgatk_all      <- load_csv_if_exists(LR_MGATK_CSV, "LR mgatk")
sr_mitoclone2_all <- load_csv_if_exists(SR_MITOCLONE2_CSV, "SR mitoClone2")
lr_mitoclone2_all <- load_csv_if_exists(LR_MITOCLONE2_CSV, "LR mitoClone2")
# copykat is loaded per-donor inside the main loop below (see load_sr_copykat_donor)

sr_meta <- NULL
if (file.exists(SEURAT_FILE)) {
  seu <- readRDS(SEURAT_FILE)
  sr_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(cell = sub("-1$", "", cell)) %>%
    dplyr::select(cell, sample_name, mapping_error_QC)
} else {
  message("Seurat file not found, cannot derive SR per-donor cell membership: ", SEURAT_FILE)
}

# ---------------------------------------------------------------
# Per-donor, per-modality method table builder
# ---------------------------------------------------------------

# donor_cells: character vector of cell barcodes belonging to this donor
# in this modality. method_tables: named list of data.frames each with a
# "cell" column plus exactly one clone/call column.
build_method_wide_table <- function(donor_cells, method_tables) {
  base <- tibble(cell = donor_cells)
  for (nm in names(method_tables)) {
    tbl <- method_tables[[nm]]
    if (is.null(tbl)) next
    value_col <- setdiff(names(tbl), "cell")[1]
    tbl_sub <- tbl %>%
      filter(cell %in% donor_cells) %>%
      dplyr::select(cell, !!nm := all_of(value_col))
    base <- base %>% left_join(tbl_sub, by = "cell")
  }
  base
}

pairwise_ari_table <- function(wide_tbl, method_cols) {
  present <- method_cols[method_cols %in% names(wide_tbl)]
  if (length(present) < 2) return(tibble())
  combn(present, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      x <- wide_tbl[[pair[1]]]
      y <- wide_tbl[[pair[2]]]
      n_shared <- sum(!is.na(x) & !is.na(y))
      tibble(
        method_1 = pair[1],
        method_2 = pair[2],
        n_cells_compared = n_shared,
        ari = adjusted_rand_index(x, y)
      )
    })
}

plot_alluvial <- function(wide_tbl, col1, col2, title, out_path) {
  if (!all(c(col1, col2) %in% names(wide_tbl))) return(invisible(NULL))
  d <- wide_tbl %>%
    filter(!is.na(.data[[col1]]), !is.na(.data[[col2]])) %>%
    # numbat's clone_opt is a plain integer (1, 2, 3...), unlike the other
    # methods' string labels - without coercing to character here, "a" stays
    # numeric and scale_fill_brewer() (discrete-only) errors on it.
    mutate(across(all_of(c(col1, col2)), as.character)) %>%
    count(.data[[col1]], .data[[col2]], name = "count")
  if (nrow(d) == 0) return(invisible(NULL))
  names(d)[1:2] <- c("a", "b")
  p <- ggplot(d, aes(y = count, axis1 = a, axis2 = b)) +
    geom_alluvium(aes(fill = a), width = 1/6, alpha = 0.7) +
    geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c(col1, col2), expand = c(0.15, 0.05)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = title, y = "Number of cells")
  ggsave(out_path, p, width = 7, height = 6)
  print(p)
}

# ---------------------------------------------------------------
# A) Within-modality, cross-method concordance
# ---------------------------------------------------------------

within_modality_results <- list()

for (lr_donor in names(donor_map)) {

  sr_donor <- donor_map[[lr_donor]]
  # donor_map values carry the "LK2_" numbat-directory prefix (real, confirmed
  # against results/seurat_annotated/.../numbat/LK2_<donor>/), but sample_name
  # in the Seurat metadata / composition / copykat tables never has that
  # prefix (confirmed against the real composition CSVs) - bare_donor is the
  # form to use for every sample_name-based filter below.
  sr_donor_bare <- sub(paste0("^", LR_SAMPLE, "_"), "", sr_donor)
  message("\n=== Within-modality comparison: ", sr_donor, " ===")

  # ---- Short-read ----
  if (!is.null(sr_meta)) {
    sr_donor_cells <- sr_meta %>%
      filter(sample_name == sr_donor_bare, mapping_error_QC == "Pass") %>%
      pull(cell)

    sr_numbat_dir <- file.path(SR_NUMBAT_BASE, sr_donor, "numbat_final")
    sr_numbat <- if (dir.exists(sr_numbat_dir)) load_clone_post(sr_numbat_dir) else NULL
    if (!is.null(sr_numbat)) {
      sr_numbat <- sr_numbat %>%
        mutate(cell = sub("-1$", "", cell)) %>%
        dplyr::select(cell, numbat = any_of("clone_opt"))
    }

    sr_copykat_donor <- load_sr_copykat_donor(sr_donor_bare)

    sr_mgatk_donor <- if (!is.null(sr_mgatk_all)) {
      sr_mgatk_all %>% dplyr::select(cell, mgatk = mgatk_clone)
    } else NULL

    sr_mitoclone2_donor <- if (!is.null(sr_mitoclone2_all)) {
      sr_mitoclone2_all %>% dplyr::select(cell, mitoclone2 = clone)
    } else NULL

    sr_wide <- build_method_wide_table(
      sr_donor_cells,
      list(numbat = sr_numbat, copykat = sr_copykat_donor,
           mgatk = sr_mgatk_donor, mitoclone2 = sr_mitoclone2_donor)
    )

    sr_ari <- pairwise_ari_table(sr_wide, c("numbat", "copykat", "mgatk", "mitoclone2")) %>%
      mutate(donor = sr_donor, modality = "Short-read")

    within_modality_results[[paste0(sr_donor, "_SR")]] <- sr_ari

    write_csv(sr_wide, file.path(OUT_TBL, paste0(sr_donor, "_SR_all_methods_wide.csv")))

    plot_alluvial(sr_wide, "numbat", "mgatk",
                  paste0(sr_donor, " (SR) — numbat vs mgatk"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_numbat_vs_mgatk_alluvial.pdf")))
    plot_alluvial(sr_wide, "mgatk", "mitoclone2",
                  paste0(sr_donor, " (SR) — mgatk vs mitoClone2"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_mgatk_vs_mitoclone2_alluvial.pdf")))
    # Direct comparisons against the two CNV-based callers, not just via mgatk -
    # ARI for these pairs was already computed above (pairwise_ari_table runs
    # combn() over all present methods), this just adds the matching alluvial plots.
    plot_alluvial(sr_wide, "numbat", "mitoclone2",
                  paste0(sr_donor, " (SR) — numbat vs mitoClone2"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_numbat_vs_mitoclone2_alluvial.pdf")))
    plot_alluvial(sr_wide, "copykat", "mitoclone2",
                  paste0(sr_donor, " (SR) — copykat vs mitoClone2"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_copykat_vs_mitoclone2_alluvial.pdf")))
    # The two CNV-based callers directly against each other, and copykat vs
    # mgatk - completes every pairwise combination of the four methods.
    plot_alluvial(sr_wide, "numbat", "copykat",
                  paste0(sr_donor, " (SR) — numbat vs copykat"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_numbat_vs_copykat_alluvial.pdf")))
    plot_alluvial(sr_wide, "copykat", "mgatk",
                  paste0(sr_donor, " (SR) — copykat vs mgatk"),
                  file.path(OUT_FIG, paste0(sr_donor, "_SR_copykat_vs_mgatk_alluvial.pdf")))
  }

  # ---- Long-read ----
  lr_numbat_dir <- file.path(LR_NUMBAT_BASE, lr_donor)
  lr_numbat <- if (dir.exists(lr_numbat_dir)) load_clone_post(lr_numbat_dir) else NULL

  if (!is.null(lr_numbat)) {
    lr_numbat <- lr_numbat %>%
      mutate(cell = sub("-1$", "", cell)) %>%
      dplyr::select(cell, numbat = any_of("clone_opt"))

    lr_donor_cells <- lr_numbat$cell

    lr_mgatk_donor <- if (!is.null(lr_mgatk_all)) {
      lr_mgatk_all %>% dplyr::select(cell, mgatk = mgatk_clone)
    } else NULL

    lr_mitoclone2_donor <- if (!is.null(lr_mitoclone2_all)) {
      lr_mitoclone2_all %>% dplyr::select(cell, mitoclone2 = clone)
    } else NULL

    lr_wide <- build_method_wide_table(
      lr_donor_cells,
      list(numbat = lr_numbat,
           mgatk = lr_mgatk_donor, mitoclone2 = lr_mitoclone2_donor)
    )

    lr_ari <- pairwise_ari_table(lr_wide, c("numbat", "mgatk", "mitoclone2")) %>%
      mutate(donor = sr_donor, modality = "Long-read")

    within_modality_results[[paste0(sr_donor, "_LR")]] <- lr_ari

    write_csv(lr_wide, file.path(OUT_TBL, paste0(sr_donor, "_LR_all_methods_wide.csv")))

    plot_alluvial(lr_wide, "numbat", "mgatk",
                  paste0(sr_donor, " (LR) — numbat vs mgatk"),
                  file.path(OUT_FIG, paste0(sr_donor, "_LR_numbat_vs_mgatk_alluvial.pdf")))
    plot_alluvial(lr_wide, "mgatk", "mitoclone2",
                  paste0(sr_donor, " (LR) — mgatk vs mitoClone2"),
                  file.path(OUT_FIG, paste0(sr_donor, "_LR_mgatk_vs_mitoclone2_alluvial.pdf")))
    # Direct numbat-vs-mitoClone2 comparison (no copykat here - no long-read
    # copykat script exists in this project, per the header comment above).
    plot_alluvial(lr_wide, "numbat", "mitoclone2",
                  paste0(sr_donor, " (LR) — numbat vs mitoClone2"),
                  file.path(OUT_FIG, paste0(sr_donor, "_LR_numbat_vs_mitoclone2_alluvial.pdf")))
  }
}

within_modality_summary <- bind_rows(within_modality_results)

if (nrow(within_modality_summary) > 0) {
  write_csv(within_modality_summary, file.path(OUT_TBL, "within_modality_method_ari_summary.csv"))

  p_ari <- ggplot(
    within_modality_summary,
    aes(x = paste(method_1, method_2, sep = " vs "), y = ari, fill = modality)
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ donor, scales = "free_x") +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = "Cross-method clustering concordance (adjusted Rand index)",
      x = NULL, y = "Adjusted Rand index", fill = NULL
    )

  ggsave(file.path(OUT_FIG, "within_modality_method_ari_summary.pdf"), p_ari, width = 10, height = 7)
  print(p_ari)
}

# ---------------------------------------------------------------
# B) Cross-modality (SR vs LR) concordance for genetic clone callers
# ---------------------------------------------------------------

cross_modality_results <- list()

for (lr_donor in names(donor_map)) {

  sr_donor <- donor_map[[lr_donor]]
  sr_donor_bare <- sub(paste0("^", LR_SAMPLE, "_"), "", sr_donor)
  message("\n=== Cross-modality (SR vs LR) genetic clone comparison: ", sr_donor, " ===")

  if (is.null(sr_mgatk_all) || is.null(lr_mgatk_all)) {
    message("  mgatk SR or LR clone table missing, skipping mgatk SR-vs-LR for this donor")
  } else if (is.null(sr_meta)) {
    message("  Seurat metadata missing, cannot restrict SR mgatk cells to this donor")
  } else {
    sr_donor_cells <- sr_meta %>%
      filter(sample_name == sr_donor_bare, mapping_error_QC == "Pass") %>%
      pull(cell)

    lr_numbat_dir <- file.path(LR_NUMBAT_BASE, lr_donor)
    lr_numbat <- if (dir.exists(lr_numbat_dir)) load_clone_post(lr_numbat_dir) else NULL
    lr_donor_cells <- if (!is.null(lr_numbat)) sub("-1$", "", lr_numbat$cell) else character(0)

    sr_mgatk_donor <- sr_mgatk_all %>% filter(cell %in% sr_donor_cells) %>%
      dplyr::select(cell, mgatk_sr = mgatk_clone)
    lr_mgatk_donor <- lr_mgatk_all %>% filter(cell %in% lr_donor_cells) %>%
      dplyr::select(cell, mgatk_lr = mgatk_clone)

    merged_mgatk <- inner_join(sr_mgatk_donor, lr_mgatk_donor, by = "cell")

    if (nrow(merged_mgatk) > 0) {
      ari_val <- adjusted_rand_index(merged_mgatk$mgatk_sr, merged_mgatk$mgatk_lr)
      message("  mgatk SR-vs-LR: ", nrow(merged_mgatk), " shared cells, ARI = ", round(ari_val, 3))

      cross_modality_results[[paste0(sr_donor, "_mgatk")]] <- tibble(
        donor = sr_donor, method = "mgatk",
        n_cells_shared = nrow(merged_mgatk), ari = ari_val
      )

      plot_alluvial(
        merged_mgatk %>% dplyr::rename(mgatk_SR = mgatk_sr, mgatk_LR = mgatk_lr),
        "mgatk_SR", "mgatk_LR",
        paste0(sr_donor, " — mgatk clones: short-read vs long-read"),
        file.path(OUT_FIG, paste0(sr_donor, "_mgatk_SR_vs_LR_alluvial.pdf"))
      )
    }
  }

  if (is.null(sr_mitoclone2_all) || is.null(lr_mitoclone2_all)) {
    message("  mitoClone2 SR or LR clone table missing, skipping mitoClone2 SR-vs-LR for this donor")
  } else if (is.null(sr_meta)) {
    message("  Seurat metadata missing, cannot restrict SR mitoClone2 cells to this donor")
  } else {
    sr_donor_cells <- sr_meta %>%
      filter(sample_name == sr_donor_bare, mapping_error_QC == "Pass") %>%
      pull(cell)

    lr_numbat_dir <- file.path(LR_NUMBAT_BASE, lr_donor)
    lr_numbat <- if (dir.exists(lr_numbat_dir)) load_clone_post(lr_numbat_dir) else NULL
    lr_donor_cells <- if (!is.null(lr_numbat)) sub("-1$", "", lr_numbat$cell) else character(0)

    sr_mc2_donor <- sr_mitoclone2_all %>% filter(cell %in% sr_donor_cells) %>%
      dplyr::select(cell, mitoclone2_sr = clone)
    lr_mc2_donor <- lr_mitoclone2_all %>% filter(cell %in% lr_donor_cells) %>%
      dplyr::select(cell, mitoclone2_lr = clone)

    merged_mc2 <- inner_join(sr_mc2_donor, lr_mc2_donor, by = "cell")

    if (nrow(merged_mc2) > 0) {
      ari_val <- adjusted_rand_index(merged_mc2$mitoclone2_sr, merged_mc2$mitoclone2_lr)
      message("  mitoClone2 SR-vs-LR: ", nrow(merged_mc2), " shared cells, ARI = ", round(ari_val, 3))

      cross_modality_results[[paste0(sr_donor, "_mitoclone2")]] <- tibble(
        donor = sr_donor, method = "mitoClone2",
        n_cells_shared = nrow(merged_mc2), ari = ari_val
      )

      plot_alluvial(
        merged_mc2 %>% dplyr::rename(mitoclone2_SR = mitoclone2_sr, mitoclone2_LR = mitoclone2_lr),
        "mitoclone2_SR", "mitoclone2_LR",
        paste0(sr_donor, " — mitoClone2 clones: short-read vs long-read"),
        file.path(OUT_FIG, paste0(sr_donor, "_mitoclone2_SR_vs_LR_alluvial.pdf"))
      )
    }
  }
}

cross_modality_summary <- bind_rows(cross_modality_results)

if (nrow(cross_modality_summary) > 0) {
  write_csv(cross_modality_summary, file.path(OUT_TBL, "cross_modality_genetic_method_ari_summary.csv"))

  p_cross <- ggplot(cross_modality_summary, aes(x = donor, y = ari, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_bw(base_size = 12) +
    labs(
      title = "Short-read vs long-read genetic clone concordance (adjusted Rand index)",
      x = NULL, y = "Adjusted Rand index", fill = NULL
    )

  ggsave(file.path(OUT_FIG, "cross_modality_genetic_method_ari_summary.pdf"), p_cross, width = 8, height = 5)
  print(p_cross)
}

message("\nDone. Figures: ", OUT_FIG, " | Tables: ", OUT_TBL)
