#!/usr/bin/env Rscript

# Shared helpers sourced by both 15a.mito_mgatk_stage.R (run under the
# mgatk env) and 15b.mito_mitoclone2_and_compare.R (run under r_env).
# Only needs base tidyverse packages (dplyr/tidyr/tibble/stringr/readr/
# purrr/ggplot2/ggalluvial) - no mgatk CLI, no Signac, no mitoClone2 - so
# it's safe to source from either environment.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(ggalluvial)
})

# =================================================================
# Barcode helpers
# =================================================================

extract_10x <- function(x) {
  x <- as.character(x)
  bc <- stringr::str_extract(x, "[ACGT]{16}(-[0-9]+)?")
  bc[is.na(bc)] <- x[is.na(bc)]
  bc
}

add_dash_one <- function(x) {
  x <- as.character(x)
  ifelse(grepl("-[0-9]+$", x), x, paste0(x, "-1"))
}

# Writes a suffix-stripped copy of a "-1"-suffixed (CellRanger-convention)
# barcode file, for use against long-read BAMs - confirmed directly that
# the FLAMES/minimap2 long-read BAMs' CB tags carry no GEM-well suffix at
# all (0 of 7367 whitelist barcodes matched any CB tag sampled from the
# real BAM before this fix), unlike the short-read CellRanger BAM.
strip_barcode_suffix <- function(barcodes_in, barcodes_out) {
  bc <- readLines(barcodes_in)
  writeLines(sub("-1$", "", bc), barcodes_out)
  invisible(barcodes_out)
}

# =================================================================
# Per-donor discovery - replaces the hand-maintained donor_map in the old
# 20a.compare_all_clone_methods.R. Returns one row per donor for the given
# SR sample, with everything needed downstream: the donor's SR barcode
# file (already "-1"-suffixed, from the pipeline's own per-donor split
# under numbat/barcodes_by_sample_name/), and whether a matching
# long-read run exists.
# =================================================================

discover_donors <- function(run_id, results_dir, sr_sample) {
  lr_sample <- sub("-GEX$", "", sr_sample)

  comp_path <- file.path(results_dir, run_id, "08_annotate", sr_sample, "tables",
                          "composition_broad_by_sample_name.csv")
  if (!file.exists(comp_path)) stop("Composition table not found: ", comp_path)

  donors <- sort(unique(read_csv(comp_path, show_col_types = FALSE)$sample_name))

  barcodes_dir <- file.path(results_dir, run_id, "numbat", "barcodes_by_sample_name", "barcodes")
  # results_dir is "<proj>/results_nf"; long_read lives at "<proj>/results/long_read"
  lr_bam <- file.path(dirname(results_dir), "results", "long_read", lr_sample, "bam", "align2genome.bam")
  has_lr <- file.exists(lr_bam)

  tibble(
    donor = donors,
    sr_sample = sr_sample,
    lr_sample = lr_sample,
    sr_barcodes_file = file.path(barcodes_dir, paste0(donors, "_barcodes.tsv")),
    has_lr = has_lr
  ) %>%
    filter(file.exists(sr_barcodes_file))
}

# =================================================================
# Cross-method comparison (ARI + alluvial) - from
# 20a.compare_all_clone_methods.R, unchanged. Pure tidyverse, no
# mgatk/mitoClone2/Signac dependency, so it can run in either
# environment - lives here and gets called from 15b (the stage that
# runs after both mgatk and mitoClone2 output already exist on disk).
# =================================================================

adjusted_rand_index <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  tab <- table(x, y)
  n <- sum(tab)
  a <- sum(choose(rowSums(tab), 2)); b <- sum(choose(colSums(tab), 2)); c_sum <- sum(choose(tab, 2))
  expected <- a * b / choose(n, 2); max_index <- (a + b) / 2
  if (max_index == expected) return(NA_real_)
  (c_sum - expected) / (max_index - expected)
}

build_method_wide_table <- function(donor_cells, method_tables) {
  base <- tibble(cell = donor_cells)
  for (nm in names(method_tables)) {
    tbl <- method_tables[[nm]]
    if (is.null(tbl)) next
    value_col <- setdiff(names(tbl), "cell")[1]
    tbl_sub <- tbl %>% filter(cell %in% donor_cells) %>% dplyr::select(cell, !!nm := all_of(value_col))
    base <- base %>% left_join(tbl_sub, by = "cell")
  }
  base
}

pairwise_ari_table <- function(wide_tbl, method_cols) {
  present <- method_cols[method_cols %in% names(wide_tbl)]
  if (length(present) < 2) return(tibble())
  combn(present, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      x <- wide_tbl[[pair[1]]]; y <- wide_tbl[[pair[2]]]
      tibble(method_1 = pair[1], method_2 = pair[2],
             n_cells_compared = sum(!is.na(x) & !is.na(y)), ari = adjusted_rand_index(x, y))
    })
}

plot_alluvial <- function(wide_tbl, col1, col2, title, out_path) {
  if (!all(c(col1, col2) %in% names(wide_tbl))) return(invisible(NULL))
  d <- wide_tbl %>%
    filter(!is.na(.data[[col1]]), !is.na(.data[[col2]])) %>%
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
    theme_bw(base_size = 12) + theme(legend.position = "none") +
    labs(title = title, y = "Number of cells")
  ggsave(out_path, p, width = 7, height = 6)
}

load_clone_post <- function(numbat_dir) {
  candidates <- list.files(numbat_dir, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) return(NULL)
  iterations <- as.integer(str_extract(basename(candidates), "\\d+"))
  read_tsv(candidates[which.max(iterations)], show_col_types = FALSE)
}

load_copykat_donor <- function(copykat_base, bare_donor) {
  path <- file.path(copykat_base, bare_donor, "tables", "copykat_prediction_with_metadata.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(cell = sub("-1$", "", cell)) %>%
    dplyr::select(cell, copykat = copykat_call)
}

# Output layout shared by both stages, so 15a writes exactly where 15b
# expects to find it.
donor_out_dir <- function(project_dir, run_id, donor) {
  file.path(project_dir, "results/mitochondrial_clones/per_patient", run_id, donor)
}
