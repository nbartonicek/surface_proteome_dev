#!/usr/bin/env Rscript

# Stage 2 of 2 - runs under r_env (needs: mitoClone2, samtools, R with the
# shared tidyverse packages in 15.mito_clone_common.R). Does NOT need the
# mgatk CLI or Signac at all - it only reads mgatk's CSV output from disk
# (written by 15a.mito_mgatk_stage.R, run separately under the mgatk env).
#
# For every donor: runs mitoClone2 clustering (SR, and LR if available),
# merges SR+LR, then builds the full cross-method comparison (numbat +
# copykat + mgatk-SR/LR/merged + mitoClone2-SR/LR/merged) - ARI table +
# alluvial plots - per donor. Run 15a first; if its output isn't there yet
# for a donor, the mgatk columns are just absent from that donor's
# comparison rather than erroring (same "skip missing input" pattern used
# throughout this project's scripts).
#
# Usage:
#   Rscript 15b.mito_mitoclone2_and_compare.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv] [n_parallel_donors] [ncores_per_donor]
#
# Donors are independent (each writes to its own donor_out_dir()) so
# n_parallel_donors donors are run concurrently via parallel::mclapply
# (forked, Unix-only), each mitoClone2/samtools call using ncores_per_donor
# cores. Keep n_parallel_donors * ncores_per_donor <= the cores actually
# allocated to the job (see --cpus-per-task in
# 15b.mito_mitoclone2_and_compare.sh), and remember memory scales with
# n_parallel_donors too since donors' runs then overlap.

suppressPackageStartupMessages({
  library(mitoClone2)
})

# Run from the scripts/ directory (same convention as the old
# 15c.mgatk_merged_clones.R, which sourced 15.mgatk.R the same way).
source("15.mito_clone_common.R")

# =================================================================
# mitoClone2 clone calling (from 15b.mitoclone2_calling.R, unchanged)
# =================================================================

call_mitoclone2_clones <- function(mt_bam, barcodes, out_dir, sample, n_clones = 4, ncores = 8) {
  stopifnot(file.exists(mt_bam), file.exists(barcodes))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (Sys.which("samtools") == "") stop("samtools not found on PATH")

  # bam2R_10x() has no barcode-whitelist filter of its own beyond
  # min_reads_per_barcode, and builds a full-length mtDNA matrix per
  # barcode it sees - pre-filter to the (donor-specific) whitelist first.
  filtered_bam <- file.path(out_dir, paste0(sample, ".mt.whitelisted.bam"))
  if (!file.exists(filtered_bam)) {
    filter_cmd <- sprintf(
      "samtools view -@ %d -b -D CB:%s %s > %s && samtools index -@ %d %s",
      ncores, shQuote(barcodes), shQuote(mt_bam), shQuote(filtered_bam), ncores, shQuote(filtered_bam)
    )
    status <- system(filter_cmd)
    if (status != 0 || !file.exists(filtered_bam)) { warning("BAM whitelist filtering failed for ", sample); return(NULL) }
  }

  baseCounts <- bam2R_10x(file = filtered_bam, sites = "chrM:1-16569", ncores = ncores)
  message(sample, ": cells with mitochondrial coverage: ", length(baseCounts))
  saveRDS(baseCounts, file.path(out_dir, paste0(sample, ".baseCounts.rds")))
  if (length(baseCounts) < 10) { warning("Too few cells with mt coverage for ", sample); return(NULL) }

  mutCalls <- mutationCallsFromExclusionlist(
    baseCounts, min.af = 0.05, min.num.samples = 5,
    universal.var.cells = 0.5 * length(baseCounts), binarize = 0.1
  )
  saveRDS(mutCalls, file.path(out_dir, paste0(sample, ".mutCalls.rds")))

  pdf(file.path(out_dir, paste0(sample, ".mitoclone2_heatmap.pdf")), width = 8, height = 10)
  clustered <- tryCatch(
    quick_cluster(mutCalls, binarize = TRUE, drop_empty = TRUE,
                  clustering.method = "ward.D2", show_colnames = FALSE, fontsize_row = 7),
    error = function(e) { warning("quick_cluster() failed for ", sample, ": ", conditionMessage(e)); NULL }
  )
  dev.off()
  if (is.null(clustered)) return(NULL)
  saveRDS(clustered, file.path(out_dir, paste0(sample, ".clustered.rds")))

  # quick_cluster()'s object is pheatmap-style with variants on rows, cells
  # on columns - confirmed from real output (an earlier version of this
  # logic cut $tree_row and silently produced a "clone" table keyed by
  # mtDNA variant id instead of cell barcode). Cut $tree_col for cells.
  if (is.null(clustered$tree_col)) {
    warning("quick_cluster() output for ", sample, " had no $tree_col"); return(NULL)
  }
  clone_labels <- cutree(clustered$tree_col, k = min(n_clones, length(clustered$tree_col$order)))

  clone_df <- data.frame(cell = sub("-1$", "", names(clone_labels)),
                          clone = paste0("mitoclone2_", clone_labels), row.names = NULL)
  write_csv(clone_df, file.path(out_dir, paste0(sample, ".mitoclone2_clones.csv")))
  message("Clone sizes:"); print(table(clone_df$clone))
  clone_df
}

call_mitoclone2_clones_merged <- function(sr_basecounts_rds, lr_basecounts_rds, out_dir, sample, n_clones = 4) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sr_counts <- readRDS(sr_basecounts_rds); lr_counts <- readRDS(lr_basecounts_rds)
  norm_bc <- function(x) sub("-1$", "", x)
  names(sr_counts) <- norm_bc(names(sr_counts)); names(lr_counts) <- norm_bc(names(lr_counts))

  shared_barcodes <- intersect(names(sr_counts), names(lr_counts))
  message("Cells shared between SR and LR mitoClone2 base counts: ", length(shared_barcodes))
  if (length(shared_barcodes) == 0) { warning("No cells shared between SR and LR baseCounts for ", sample); return(NULL) }

  merged_counts <- lapply(shared_barcodes, function(bc) {
    sr_mat <- sr_counts[[bc]]; lr_mat <- lr_counts[[bc]]
    if (!identical(dim(sr_mat), dim(lr_mat)) || !identical(colnames(sr_mat), colnames(lr_mat))) {
      stop("SR/LR base-count matrices for cell ", bc, " have mismatched dimensions")
    }
    sr_mat + lr_mat
  })
  names(merged_counts) <- shared_barcodes

  mutCalls <- mutationCallsFromExclusionlist(
    merged_counts, min.af = 0.05, min.num.samples = 5,
    universal.var.cells = 0.5 * length(merged_counts), binarize = 0.1
  )
  pdf(file.path(out_dir, paste0(sample, ".merged_mitoclone2_heatmap.pdf")), width = 8, height = 10)
  clustered <- tryCatch(
    quick_cluster(mutCalls, binarize = TRUE, drop_empty = TRUE,
                  clustering.method = "ward.D2", show_colnames = FALSE, fontsize_row = 7),
    error = function(e) { warning("quick_cluster() failed for merged ", sample, ": ", conditionMessage(e)); NULL }
  )
  dev.off()
  if (is.null(clustered) || is.null(clustered$tree_col)) return(NULL)

  clone_labels <- cutree(clustered$tree_col, k = min(n_clones, length(clustered$tree_col$order)))
  clone_df <- data.frame(cell = sub("-1$", "", names(clone_labels)),
                          clone = paste0("mitoclone2_merged_", clone_labels), row.names = NULL)
  write_csv(clone_df, file.path(out_dir, paste0(sample, ".mitoclone2_merged_clones.csv")))
  clone_df
}

# =================================================================
# Per-donor driver: mitoClone2 stage + full comparison
# =================================================================

run_donor_mitoclone2_and_compare <- function(run_id, results_dir, project_dir, donor_row, n_clones = 4, ncores = 8) {
  donor <- donor_row$donor
  sr_sample <- donor_row$sr_sample
  lr_sample <- donor_row$lr_sample
  message("\n========== [mitoClone2+compare] Donor: ", donor, " ==========")

  out_base <- donor_out_dir(project_dir, run_id, donor)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  sr_mt_bam <- file.path(project_dir, "results/mitochondrial_clones", run_id, sr_sample, paste0(sr_sample, ".mt.bam"))
  if (!file.exists(sr_mt_bam)) {
    warning("SR .mt.bam not found for ", sr_sample); return(invisible(NULL))
  }

  method_tables <- list()

  # ---- mitoClone2: SR ----
  # tryCatch'd like the merge call below - mitoClone2's own internals (e.g.
  # mutationCallsFromExclusionlist()'s GRanges construction) can throw for a
  # particular donor's data (observed for normal-01: "must contain strings
  # of the form chr:start-end..." - a mitoClone2-internal issue, not ours).
  # One method failing for one donor must not lose numbat/copykat/mgatk too.
  sr_mc2_dir <- file.path(out_base, "sr_mitoclone2")
  sr_mc2_clones <- tryCatch(
    call_mitoclone2_clones(sr_mt_bam, donor_row$sr_barcodes_file, sr_mc2_dir, donor, n_clones, ncores),
    error = function(e) { warning("mitoClone2 SR failed for ", donor, ": ", conditionMessage(e)); NULL }
  )
  if (!is.null(sr_mc2_clones)) method_tables$mitoclone2_sr <- sr_mc2_clones %>% dplyr::select(cell, mitoclone2_sr = clone)

  # ---- mitoClone2: LR + merge ----
  if (isTRUE(donor_row$has_lr)) {
    lr_mt_bam <- file.path(project_dir, "results/mitochondrial_clones_longread", lr_sample, paste0(lr_sample, ".mt.bam"))
    if (file.exists(lr_mt_bam)) {
      lr_barcodes_file <- file.path(out_base, paste0(donor, "_barcodes_nosuffix.txt"))
      if (!file.exists(lr_barcodes_file)) strip_barcode_suffix(donor_row$sr_barcodes_file, lr_barcodes_file)

      lr_mc2_dir <- file.path(out_base, "lr_mitoclone2")
      lr_mc2_clones <- tryCatch(
        call_mitoclone2_clones(lr_mt_bam, lr_barcodes_file, lr_mc2_dir, donor, n_clones, ncores),
        error = function(e) { warning("mitoClone2 LR failed for ", donor, ": ", conditionMessage(e)); NULL }
      )
      if (!is.null(lr_mc2_clones)) {
        method_tables$mitoclone2_lr <- lr_mc2_clones %>% dplyr::select(cell, mitoclone2_lr = clone)

        sr_basecounts <- file.path(sr_mc2_dir, paste0(donor, ".baseCounts.rds"))
        lr_basecounts <- file.path(lr_mc2_dir, paste0(donor, ".baseCounts.rds"))
        if (file.exists(sr_basecounts) && file.exists(lr_basecounts)) {
          merged_mc2 <- tryCatch(
            call_mitoclone2_clones_merged(sr_basecounts, lr_basecounts, file.path(out_base, "merged_mitoclone2"), donor, n_clones),
            error = function(e) { warning("mitoClone2 merge failed for ", donor, ": ", conditionMessage(e)); NULL }
          )
          if (!is.null(merged_mc2)) method_tables$mitoclone2_merged <- merged_mc2 %>% dplyr::select(cell, mitoclone2_merged = clone)
        }
      }
    } else {
      message("No long-read .mt.bam found at ", lr_mt_bam, " - skipping long-read mitoClone2 for this donor")
    }
  }

  # ---- mgatk output from stage 1, if it exists ----
  read_mgatk_csv <- function(suffix, colname) {
    path <- file.path(out_base, paste0(donor, suffix))
    if (!file.exists(path)) { message("  ", path, " not found - run 15a.mito_mgatk_stage.R first for this to be included"); return(NULL) }
    read_csv(path, show_col_types = FALSE) %>% dplyr::select(cell, !!colname := mgatk_clone)
  }
  mgatk_sr <- read_mgatk_csv(".sr_mgatk_clones.csv", "mgatk_sr")
  if (!is.null(mgatk_sr)) method_tables$mgatk_sr <- mgatk_sr
  mgatk_lr <- read_mgatk_csv(".lr_mgatk_clones.csv", "mgatk_lr")
  if (!is.null(mgatk_lr)) method_tables$mgatk_lr <- mgatk_lr
  mgatk_merged <- read_mgatk_csv(".merged_mgatk_clones.csv", "mgatk_merged")
  if (!is.null(mgatk_merged)) method_tables$mgatk_merged <- mgatk_merged

  # ---- Numbat + CopyKAT (already per-donor - just load) ----
  numbat_dir <- file.path(results_dir, run_id, "09c_numbat_run", paste0(lr_sample, "_", donor), "numbat_final")
  numbat_tbl <- if (dir.exists(numbat_dir)) load_clone_post(numbat_dir) else NULL
  if (!is.null(numbat_tbl)) {
    method_tables$numbat <- numbat_tbl %>% mutate(cell = sub("-1$", "", cell)) %>%
      dplyr::select(cell, numbat = any_of("clone_opt"))
  }
  copykat_tbl <- load_copykat_donor(file.path(results_dir, run_id, "09_copykat"), donor)
  if (!is.null(copykat_tbl)) method_tables$copykat <- copykat_tbl

  # ---- Consolidated comparison for this donor ----
  donor_cells <- readLines(donor_row$sr_barcodes_file) %>% sub("-1$", "", .)
  wide <- build_method_wide_table(donor_cells, method_tables)
  write_csv(wide, file.path(out_base, paste0(donor, "_all_methods_wide.csv")))

  all_methods <- c("numbat", "copykat", "mgatk_sr", "mgatk_lr", "mgatk_merged",
                    "mitoclone2_sr", "mitoclone2_lr", "mitoclone2_merged")
  ari <- pairwise_ari_table(wide, all_methods) %>% mutate(donor = donor)
  write_csv(ari, file.path(out_base, paste0(donor, "_ari_summary.csv")))

  present <- all_methods[all_methods %in% names(wide)]
  if (length(present) >= 2) {
    for (pair in combn(present, 2, simplify = FALSE)) {
      plot_alluvial(wide, pair[1], pair[2], paste0(donor, ": ", pair[1], " vs ", pair[2]),
                    file.path(out_base, paste0(donor, "_", pair[1], "_vs_", pair[2], "_alluvial.pdf")))
    }
  }

  message("Donor ", donor, " done. Methods available: ", paste(present, collapse = ", "))
  invisible(list(wide = wide, ari = ari))
}

# =================================================================
# CLI entry point
# =================================================================

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript 15b.mito_mitoclone2_and_compare.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv] [n_parallel_donors] [ncores_per_donor]")
  }

  RUN_ID <- args[1]
  SR_SAMPLE <- if (length(args) >= 2) args[2] else "LK2-GEX"
  RESULTS_DIR <- if (length(args) >= 3) args[3] else "/scratch/users/nbartonicek/projects/amgen/results_nf"
  PROJECT_DIR <- dirname(RESULTS_DIR)
  EXCLUDE_DONORS <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else character(0)
  N_PARALLEL_DONORS <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1
  NCORES_PER_DONOR <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else 8

  donors <- discover_donors(RUN_ID, RESULTS_DIR, SR_SAMPLE) %>% filter(!donor %in% EXCLUDE_DONORS)
  message("Donors found: ", paste(donors$donor, collapse = ", "))
  message("Running ", nrow(donors), " donor(s), ", N_PARALLEL_DONORS, " in parallel, ",
          NCORES_PER_DONOR, " core(s) each.")

  results <- parallel::mclapply(seq_len(nrow(donors)), function(i) {
    tryCatch(
      run_donor_mitoclone2_and_compare(RUN_ID, RESULTS_DIR, PROJECT_DIR, donors[i, ], ncores = NCORES_PER_DONOR),
      error = function(e) { message("FAILED for donor ", donors$donor[i], ": ", conditionMessage(e)); NULL }
    )
  }, mc.cores = N_PARALLEL_DONORS, mc.preschedule = FALSE)
  names(results) <- donors$donor

  message("\nAll donors processed. Failed: ", paste(names(results)[sapply(results, is.null)], collapse = ", "))
}
