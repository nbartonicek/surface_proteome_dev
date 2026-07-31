#!/usr/bin/env Rscript

# Recovery/refresh utility - re-runs ONLY the cross-method comparison step
# of 15b.mito_mitoclone2_and_compare.R (build_method_wide_table +
# pairwise_ari_table + plot_alluvial), reading whatever per-method clone
# CSVs already exist on disk (mgatk_sr/lr/merged from 15a, mitoclone2_sr/
# lr/merged from 15b, numbat/copykat from the main pipeline). Does NOT
# rerun mgatk tenx or mitoClone2 at all - useful after backfilling one
# method's output (e.g. 15a2.mito_mgatk_remerge.R for mgatk_merged) without
# re-paying for mitoClone2's samtools+kmeans steps, which already succeeded
# for donors whose CSVs are already there. A donor missing some methods
# (e.g. normal-01, where mitoClone2 itself fails on this donor's data) just
# gets a comparison over whatever methods ARE available.
#
# Run under r_env (sources 15b.mito_mitoclone2_and_compare.R for its
# function defs, which needs the mitoClone2 library loadable even though
# this script never calls into it).
#
# Usage:
#   Rscript 15b2.mito_compare_only.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv]

source("15b.mito_mitoclone2_and_compare.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 15b2.mito_compare_only.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv]")
}

RUN_ID <- args[1]
SR_SAMPLE <- if (length(args) >= 2) args[2] else "LK2-GEX"
RESULTS_DIR <- if (length(args) >= 3) args[3] else "/scratch/users/nbartonicek/projects/amgen/results_nf"
PROJECT_DIR <- dirname(RESULTS_DIR)
EXCLUDE_DONORS <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else character(0)

read_clone_csv <- function(path, colname) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% dplyr::select(cell, !!colname := clone)
}
read_mgatk_csv <- function(path, colname) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% dplyr::select(cell, !!colname := mgatk_clone)
}

donors <- discover_donors(RUN_ID, RESULTS_DIR, SR_SAMPLE) %>% filter(!donor %in% EXCLUDE_DONORS)
message("Donors found: ", paste(donors$donor, collapse = ", "))

all_methods <- c("numbat", "copykat", "mgatk_sr", "mgatk_lr", "mgatk_merged",
                  "mitoclone2_sr", "mitoclone2_lr", "mitoclone2_merged")

for (i in seq_len(nrow(donors))) {
  donor <- donors$donor[i]
  lr_sample <- donors$lr_sample[i]
  message("\n========== [compare-only] Donor: ", donor, " ==========")
  out_base <- donor_out_dir(PROJECT_DIR, RUN_ID, donor)

  method_tables <- list(
    mgatk_sr = read_mgatk_csv(file.path(out_base, paste0(donor, ".sr_mgatk_clones.csv")), "mgatk_sr"),
    mgatk_lr = read_mgatk_csv(file.path(out_base, paste0(donor, ".lr_mgatk_clones.csv")), "mgatk_lr"),
    mgatk_merged = read_mgatk_csv(file.path(out_base, paste0(donor, ".merged_mgatk_clones.csv")), "mgatk_merged"),
    mitoclone2_sr = read_clone_csv(file.path(out_base, "sr_mitoclone2", paste0(donor, ".mitoclone2_clones.csv")), "mitoclone2_sr"),
    mitoclone2_lr = read_clone_csv(file.path(out_base, "lr_mitoclone2", paste0(donor, ".mitoclone2_clones.csv")), "mitoclone2_lr"),
    mitoclone2_merged = read_clone_csv(file.path(out_base, "merged_mitoclone2", paste0(donor, ".mitoclone2_merged_clones.csv")), "mitoclone2_merged")
  )
  method_tables <- method_tables[!sapply(method_tables, is.null)]

  numbat_dir <- file.path(RESULTS_DIR, RUN_ID, "09c_numbat_run", paste0(lr_sample, "_", donor), "numbat_final")
  numbat_tbl <- if (dir.exists(numbat_dir)) load_clone_post(numbat_dir) else NULL
  if (!is.null(numbat_tbl)) {
    method_tables$numbat <- numbat_tbl %>% mutate(cell = sub("-1$", "", cell)) %>%
      dplyr::select(cell, numbat = any_of("clone_opt"))
  }
  copykat_tbl <- load_copykat_donor(file.path(RESULTS_DIR, RUN_ID, "09_copykat"), donor)
  if (!is.null(copykat_tbl)) method_tables$copykat <- copykat_tbl

  message("Methods available: ", paste(names(method_tables), collapse = ", "))
  if (length(method_tables) == 0) { message("  no method output found for ", donor, " - skipping"); next }

  donor_cells <- readLines(donors$sr_barcodes_file[i]) %>% sub("-1$", "", .)
  wide <- build_method_wide_table(donor_cells, method_tables)
  write_csv(wide, file.path(out_base, paste0(donor, "_all_methods_wide.csv")))

  ari <- pairwise_ari_table(wide, all_methods) %>% mutate(donor = donor)
  write_csv(ari, file.path(out_base, paste0(donor, "_ari_summary.csv")))

  present <- all_methods[all_methods %in% names(wide)]
  if (length(present) >= 2) {
    for (pair in combn(present, 2, simplify = FALSE)) {
      plot_alluvial(wide, pair[1], pair[2], paste0(donor, ": ", pair[1], " vs ", pair[2]),
                    file.path(out_base, paste0(donor, "_", pair[1], "_vs_", pair[2], "_alluvial.pdf")))
    }
  }
}

message("\nCompare-only pass complete.")
