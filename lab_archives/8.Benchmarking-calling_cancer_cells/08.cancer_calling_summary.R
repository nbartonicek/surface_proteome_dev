#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Calling cancer cells - step 08
#
# Scores the callers against each other rather than producing more calls.
# Reads only the tables written by steps 01-07, so it runs in seconds and does
# not need the Seurat objects, the BAMs or the container.
#
# Four questions:
#   1. what fraction of each donor does each caller call malignant
#   2. where CopyKAT and Numbat are both available, how far do they agree
#   3. what does the normal donor get called - the false positive floor
#   4. short read versus long read: compartment concordance and pCNV agreement
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

options(scipen = 999)

res     <- "../results"
out_dir <- file.path(res, "benchmarking", "cancer_calling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PILOT <- "260423_VH01624_453_222HWMYNX"
LK2   <- "260528_VH01624_464_222K7VKNX"

read_if <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  out <- try(suppressMessages(read_csv(path, show_col_types = FALSE, ...)), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

# ----------------------------
# 1. CopyKAT calls per donor
# ----------------------------

ck_root <- file.path(res, "seurat_annotated", PILOT, "copykat_annotated")
ck_donors <- list.dirs(ck_root, recursive = FALSE, full.names = TRUE)
ck_donors <- ck_donors[!grepl("(^|/)(results|\\._)", basename(ck_donors))]

copykat_calls <- map_dfr(ck_donors, function(d) {
  f <- file.path(d, "copykat_prediction_with_metadata.csv")
  x <- read_if(f)
  if (is.null(x)) return(NULL)
  x$donor <- basename(d)
  x
})

if (!is.null(copykat_calls) && nrow(copykat_calls)) {

  copykat_summary <- copykat_calls %>%
    count(donor, copykat_call, name = "n_cells") %>%
    group_by(donor) %>%
    mutate(total = sum(n_cells), percent = round(100 * n_cells / total, 2)) %>%
    ungroup() %>%
    arrange(donor, desc(percent))

  write_csv(copykat_summary, file.path(out_dir, "01_copykat_calls_by_donor.csv"))

  p_ck <- copykat_summary %>%
    filter(copykat_call == "aneuploid") %>%
    ggplot(aes(x = reorder(donor, percent), y = percent)) +
    geom_col(fill = "firebrick", width = 0.7) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(title = "CopyKAT: percent of cells called aneuploid",
         subtitle = "MOLM13 is not shown - it is a cell line and was assigned, not called",
         x = NULL, y = "Percent aneuploid")

  ggsave(file.path(out_dir, "01_copykat_aneuploid_by_donor.pdf"), p_ck,
         width = 8, height = 4.5)
} else {
  message("No CopyKAT prediction tables found under ", ck_root)
  copykat_summary <- NULL
}

# ----------------------------
# 2. Numbat compartment per donor, both runs
# ----------------------------

numbat_summary <- map_dfr(c(LK1 = PILOT, LK2 = LK2), function(r) {
  f <- file.path(res, "seurat_annotated", r, "numbat",
                 "numbat_compartment_by_sample_and_broad_annotation.csv")
  x <- read_if(f)
  if (is.null(x)) return(NULL)
  x %>%
    group_by(sample, numbat_compartment) %>%
    summarise(n_cells = sum(n), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(total = sum(n_cells), percent = round(100 * n_cells / total, 2)) %>%
    ungroup() %>%
    mutate(run = r)
}, .id = "run_label")

if (nrow(numbat_summary)) {
  write_csv(numbat_summary, file.path(out_dir, "02_numbat_compartment_by_donor.csv"))
} else {
  message("No Numbat compartment tables found - run step 04 first.")
}

# ----------------------------
# 3. CopyKAT against Numbat, cell by cell
# ----------------------------
# Only possible where both callers ran on the same donor. The comparison is on
# the malignant/normal call, not on clones - the two callers do not share a
# clone vocabulary.

numbat_meta <- read_if(file.path(res, "seurat_annotated", PILOT, "numbat",
                                 "LK1_projected_CITE_DSB_Numbat_integrated_metadata.csv"))

if (!is.null(numbat_meta) && !is.null(copykat_calls)) {

  names(numbat_meta)[1] <- "cell"

  joint <- numbat_meta %>%
    dplyr::select(cell, sample_name, numbat_compartment,
                  predicted_CellType_Broad) %>%
    inner_join(
      copykat_calls %>% dplyr::select(cell, donor, copykat_call),
      by = "cell"
    ) %>%
    filter(numbat_compartment %in% c("tumor", "normal"),
           copykat_call %in% c("aneuploid", "diploid")) %>%
    mutate(
      numbat  = ifelse(numbat_compartment == "tumor", "malignant", "normal"),
      copykat = ifelse(copykat_call == "aneuploid", "malignant", "normal")
    )

  if (nrow(joint)) {

    concordance <- joint %>%
      count(donor, numbat, copykat, name = "n_cells") %>%
      group_by(donor) %>%
      mutate(percent = round(100 * n_cells / sum(n_cells), 2)) %>%
      ungroup()

    write_csv(concordance, file.path(out_dir, "03_copykat_vs_numbat_confusion.csv"))

    agreement <- joint %>%
      group_by(donor) %>%
      summarise(
        n_cells = n(),
        percent_agree = round(100 * mean(numbat == copykat), 2),
        n_numbat_only = sum(numbat == "malignant" & copykat == "normal"),
        n_copykat_only = sum(copykat == "malignant" & numbat == "normal"),
        .groups = "drop"
      ) %>%
      arrange(desc(percent_agree))

    write_csv(agreement, file.path(out_dir, "04_copykat_vs_numbat_agreement.csv"))

    p_conf <- ggplot(concordance, aes(copykat, numbat, fill = percent)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = paste0(n_cells, "\n", percent, "%")), size = 3) +
      facet_wrap(~ donor) +
      scale_fill_gradient(low = "grey95", high = "firebrick") +
      theme_bw(base_size = 12) +
      labs(title = "CopyKAT against Numbat, same cells",
           x = "CopyKAT", y = "Numbat", fill = "% of donor")

    ggsave(file.path(out_dir, "02_copykat_vs_numbat_confusion.pdf"), p_conf,
           width = 10, height = 4.5)

    # Which cell types the two disagree on - a caller that is wrong tends to be
    # wrong in one compartment
    disagree_by_type <- joint %>%
      mutate(agree = numbat == copykat) %>%
      group_by(predicted_CellType_Broad) %>%
      summarise(n_cells = n(), percent_agree = round(100 * mean(agree), 1),
                .groups = "drop") %>%
      filter(n_cells >= 20) %>%
      arrange(percent_agree)

    write_csv(disagree_by_type,
              file.path(out_dir, "05_copykat_vs_numbat_agreement_by_celltype.csv"))
  } else {
    message("No cells with both a CopyKAT and a Numbat call.")
    agreement <- NULL
  }
} else {
  message("Need both step 01 and step 04 outputs for the caller comparison.")
  agreement <- NULL
}

# ----------------------------
# 4. Short read against long read
# ----------------------------

lr_summary <- read_if(file.path(res, "longread_numbat", "comparison_tables",
                                "comparison_summary.csv"))

if (!is.null(lr_summary)) {

  lr_tidy <- lr_summary %>%
    transmute(
      donor,
      n_cells = n_cells_overlap,
      sr_percent_tumor = round(100 * n_sr_tumor / n_sr_total, 2),
      lr_percent_tumor = round(100 * n_lr_tumor / n_lr_total, 2),
      compartment_concordance = round(100 * compartment_concordance, 2),
      pcnv_correlation = round(pcnv_correlation, 3),
      sr_tumor_lr_normal = n_sr_tumor_lr_normal,
      sr_normal_lr_tumor = n_sr_normal_lr_tumor
    ) %>%
    arrange(desc(compartment_concordance))

  write_csv(lr_tidy, file.path(out_dir, "06_shortread_vs_longread.csv"))

  p_lr <- lr_tidy %>%
    dplyr::select(donor, sr_percent_tumor, lr_percent_tumor) %>%
    pivot_longer(-donor, names_to = "platform", values_to = "percent_tumor") %>%
    mutate(platform = recode(platform,
                             sr_percent_tumor = "short read",
                             lr_percent_tumor = "long read")) %>%
    ggplot(aes(donor, percent_tumor, fill = platform)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("short read" = "#1B9E77", "long read" = "#E31A1C")) +
    theme_bw(base_size = 12) +
    labs(title = "Percent of cells called tumour, short read against long read",
         subtitle = "same cells, same caller, different sequencing chemistry",
         x = NULL, y = "Percent tumour", fill = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(out_dir, "03_shortread_vs_longread_tumor_fraction.pdf"), p_lr,
         width = 9, height = 4.5)

  # The disagreement is one-directional; show it as such
  p_dir <- lr_tidy %>%
    dplyr::select(donor, sr_tumor_lr_normal, sr_normal_lr_tumor) %>%
    pivot_longer(-donor, names_to = "direction", values_to = "n_cells") %>%
    mutate(direction = recode(direction,
                              sr_tumor_lr_normal = "short read tumour, long read normal",
                              sr_normal_lr_tumor = "short read normal, long read tumour")) %>%
    ggplot(aes(donor, n_cells, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("short read tumour, long read normal" = "#E31A1C",
                                 "short read normal, long read tumour" = "#2171B5")) +
    theme_bw(base_size = 12) +
    labs(title = "Which way the two platforms disagree",
         x = NULL, y = "Cells", fill = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "bottom")

  ggsave(file.path(out_dir, "04_shortread_vs_longread_direction.pdf"), p_dir,
         width = 9, height = 5)
} else {
  message("No long-read comparison summary found - run step 07 first.")
  lr_tidy <- NULL
}

# ----------------------------
# Report
# ----------------------------

cat("\n=== CopyKAT calls by donor ===\n")
if (!is.null(copykat_summary)) print(as.data.frame(copykat_summary))

cat("\n=== Numbat compartment by donor ===\n")
if (nrow(numbat_summary)) print(as.data.frame(numbat_summary))

cat("\n=== CopyKAT vs Numbat agreement ===\n")
if (!is.null(agreement)) print(as.data.frame(agreement))

cat("\n=== Short read vs long read ===\n")
if (!is.null(lr_tidy)) print(as.data.frame(lr_tidy))

cat("\nDone. Output directory:", out_dir, "\n")
