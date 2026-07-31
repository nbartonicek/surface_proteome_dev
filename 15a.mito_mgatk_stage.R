#!/usr/bin/env Rscript

# Stage 1 of 2 - runs under the `mgatk` conda env (needs: mgatk CLI,
# samtools, R with Matrix/Signac + the shared tidyverse packages in
# 15.mito_clone_common.R). Does NOT need mitoClone2 at all.
#
# For every donor in a run: runs `mgatk tenx` (SR, and LR if that GEM well
# has matching long-read data) scoped to just that donor's cells, converts
# the output to a consolidated signac.rds, clusters into clones, and - if
# both modalities succeeded - merges them. Writes one CSV per stage to
# donor_out_dir() so stage 2 (15b.mito_mitoclone2_and_compare.R, run under
# r_env) can pick them up later.
#
# Per-donor scoping rationale, and why this used to be one combined script:
# see the header comment in the (now-split) original 15.mito_clone_analysis.R
# history / 15.mito_clone_common.R's discover_donors() - short version: mgatk
# and mitoClone2 used to run on the WHOLE pooled GEM well, which was checked
# directly and found to silently drop 3 of 4 donors (100% of surviving cells
# belonged to one donor). Numbat/CopyKAT already ran per-donor; this brings
# mgatk in line using the same per-donor barcode files they already use.
#
# Usage:
#   Rscript 15a.mito_mgatk_stage.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv] [n_parallel_donors] [ncores_per_donor]
#
# Donors are independent (each writes to its own donor_out_dir(), reads its
# own barcode file/mt.bam) so n_parallel_donors donors are run concurrently
# via parallel::mclapply (forked, Unix-only - fine on the cluster), each
# mgatk tenx call using ncores_per_donor cores. Keep n_parallel_donors *
# ncores_per_donor <= the cores actually allocated to the job (see
# --cpus-per-task in 15a.mito_mgatk_stage.sh), and remember memory scales
# with n_parallel_donors too since donors' mgatk runs then overlap.

suppressPackageStartupMessages({
  library(Matrix)
  library(Signac)
})

# Run from the scripts/ directory (same convention as the old
# 15c.mgatk_merged_clones.R, which sourced 15.mgatk.R the same way).
source("15.mito_clone_common.R")

# =================================================================
# mgatk CLI invocation
# =================================================================

# mt_bam:    whole-well (not donor-subsetted) *.mt.bam - already extracted
#            once per modality; mgatk's own -b filters to the barcode list
# barcodes:  donor-specific barcode file (SR: "-1" suffix; LR: stripped)
# is_longread: loosens --NMmax/--NHmax (defaults 4/1, tuned for near-perfect
#            Illumina reads) - confirmed against the real long-read BAM:
#            median NM is 5 (already above the default cutoff), only 42% of
#            reads pass at NM<=4 vs 99.6% at NM<=100; separately, 0 of 200k
#            sampled long-read alignments carry an NH tag at all.
run_mgatk_tenx <- function(mt_bam, barcodes, out_dir, sample, ncores = 8, is_longread = FALSE) {
  stopifnot(file.exists(mt_bam), file.exists(barcodes))

  # mgatk's internal Snakemake workflow triggers rules off output-file
  # existence, not correctness - if out_dir already has files in it (e.g.
  # from a previous attempt made before a fix), it silently reuses them
  # instead of regenerating. Always start clean.
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

  nmax <- if (is_longread) 100 else 4
  hmax <- if (is_longread) 100 else 1

  cmd <- sprintf(
    "mgatk tenx -i %s -o %s -n %s -c %d -b %s --barcode-tag CB --umi-barcode UB --NMmax %d --NHmax %d",
    shQuote(mt_bam), shQuote(out_dir), shQuote(sample), ncores, shQuote(barcodes), nmax, hmax
  )
  message("Running: ", cmd)
  status <- system(cmd)

  final_dir <- file.path(out_dir, "final")
  depth_file <- list.files(final_dir, pattern = "\\.depthTable\\.txt$", full.names = TRUE)

  if (status != 0 || length(depth_file) != 1) {
    warning("mgatk tenx did not produce the expected output for ", sample, " in ", out_dir)
    return(NA_character_)
  }
  final_dir
}

mgatk_to_signac <- function(mgatk_final_dir) {
  depth_files <- list.files(mgatk_final_dir, pattern = "\\.depthTable\\.txt$", full.names = FALSE)
  if (length(depth_files) != 1) {
    warning("Expected exactly one *.depthTable.txt in ", mgatk_final_dir, ", found ", length(depth_files))
    return(NA_character_)
  }
  sample_prefix <- sub("\\.depthTable\\.txt$", "", depth_files)
  out_path <- file.path(mgatk_final_dir, paste0(sample_prefix, ".signac.rds"))
  if (file.exists(out_path)) return(out_path)

  mgatk_obj <- tryCatch(ReadMGATK(dir = mgatk_final_dir), error = function(e) {
    warning("ReadMGATK() failed for ", mgatk_final_dir, ": ", conditionMessage(e)); NULL
  })
  if (is.null(mgatk_obj) || !all(c("counts", "depth") %in% names(mgatk_obj))) return(NA_character_)

  saveRDS(mgatk_obj, out_path)
  out_path
}

# =================================================================
# mgatk clone calling (from 15.mgatk.R, unchanged)
# =================================================================

load_mgatk_matrix <- function(mgatk_rds) {
  stopifnot(file.exists(mgatk_rds))
  mgatk_obj <- readRDS(mgatk_rds)
  if (!all(c("counts", "depth") %in% names(mgatk_obj))) stop("mgatk_obj must contain counts and depth: ", mgatk_rds)

  mgatk_depth <- mgatk_obj$depth %>%
    as.data.frame() %>%
    rownames_to_column("mgatk_depth_row") %>%
    mutate(barcode = add_dash_one(extract_10x(sample)), mgatk_mt_depth = as.numeric(depth))

  mt_counts <- mgatk_obj$counts
  if (ncol(mt_counts) != nrow(mgatk_depth)) stop("ncol(counts) != nrow(depth) in ", mgatk_rds)
  colnames(mt_counts) <- mgatk_depth$barcode

  list(counts = mt_counts, depth = setNames(mgatk_depth$mgatk_mt_depth, mgatk_depth$barcode))
}

cluster_mgatk_counts <- function(mt_counts, mt_depth, n_clones = 4, min_frac = 0.01, max_frac = 0.80,
                                  seed = 1, label = "mgatk") {
  variant_detected_fraction <- Matrix::rowMeans(mt_counts > 0)
  keep_variants <- which(variant_detected_fraction >= min_frac & variant_detected_fraction <= max_frac)
  mt_counts_filt <- mt_counts[keep_variants, , drop = FALSE]

  message(label, ": informative mt variants retained: ", nrow(mt_counts_filt))
  if (nrow(mt_counts_filt) < 2) { warning("Too few informative mt variants retained for ", label); return(NULL) }

  mt_binary <- as.matrix(t(mt_counts_filt > 0))
  set.seed(seed)
  n_clones_use <- min(n_clones, nrow(mt_binary))
  km <- kmeans(mt_binary, centers = n_clones_use, nstart = 50)

  data.frame(
    cell = sub("-1$", "", rownames(mt_binary)),
    mgatk_clone = paste0("mgatk_clone_", km$cluster),
    mgatk_mt_depth = mt_depth[rownames(mt_binary)],
    row.names = NULL
  )
}

call_mgatk_clones <- function(mgatk_rds, n_clones = 4, min_frac = 0.01, max_frac = 0.80, seed = 1) {
  m <- load_mgatk_matrix(mgatk_rds)
  cluster_mgatk_counts(m$counts, m$depth, n_clones, min_frac, max_frac, seed, label = basename(mgatk_rds))
}

call_mgatk_clones_merged <- function(sr_rds, lr_rds, n_clones = 4, min_frac = 0.01, max_frac = 0.80, seed = 1) {
  sr <- load_mgatk_matrix(sr_rds); lr <- load_mgatk_matrix(lr_rds)
  shared_cells <- intersect(colnames(sr$counts), colnames(lr$counts))
  message("Cells shared between SR and LR mgatk runs: ", length(shared_cells))
  if (length(shared_cells) == 0) { warning("No cells shared between SR and LR mgatk outputs"); return(NULL) }

  # load_mgatk_matrix() never sets rownames(counts) - each row is one of
  # mgatk's fixed (position x base x strand) cells over the mtDNA genome
  # (16569bp x 4 bases x 2 strands = 132552 rows, always in the same order
  # regardless of input data), so SR and LR are already row-aligned by
  # construction - confirmed nrow(sr$counts) == nrow(lr$counts) == 132552
  # directly against real output. No variant-name reindexing is needed (an
  # earlier version tried to union/intersect by rownames(), which were
  # always NULL - that silently collapsed every donor's merge to 0 rows).
  # Only the cell (column) axis actually needs aligning, via shared_cells.
  stopifnot(nrow(sr$counts) == nrow(lr$counts))
  merged_counts <- sr$counts[, shared_cells, drop = FALSE] + lr$counts[, shared_cells, drop = FALSE]
  merged_depth <- sr$depth[shared_cells] + lr$depth[shared_cells]
  cluster_mgatk_counts(merged_counts, merged_depth, n_clones, min_frac, max_frac, seed, label = "merged SR+LR")
}

# =================================================================
# Per-donor driver
# =================================================================

run_donor_mgatk <- function(run_id, results_dir, project_dir, donor_row, n_clones = 4, ncores = 8) {
  donor <- donor_row$donor
  sr_sample <- donor_row$sr_sample
  lr_sample <- donor_row$lr_sample
  message("\n========== [mgatk] Donor: ", donor, " ==========")

  out_base <- donor_out_dir(project_dir, run_id, donor)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  sr_mt_bam <- file.path(project_dir, "results/mitochondrial_clones", run_id, sr_sample, paste0(sr_sample, ".mt.bam"))
  if (!file.exists(sr_mt_bam)) {
    warning("SR .mt.bam not found for ", sr_sample, " - run the whole-well mitochondrial extraction first")
    return(invisible(NULL))
  }

  sr_mgatk_dir <- file.path(out_base, "sr_mgatk")
  sr_mgatk_final <- run_mgatk_tenx(sr_mt_bam, donor_row$sr_barcodes_file, sr_mgatk_dir, donor, ncores, is_longread = FALSE)
  sr_signac <- if (!is.na(sr_mgatk_final)) mgatk_to_signac(sr_mgatk_final) else NA_character_
  if (!is.na(sr_signac)) {
    sr_clones <- tryCatch(call_mgatk_clones(sr_signac, n_clones = n_clones), error = function(e) {
      warning("SR call_mgatk_clones() failed for ", donor, ": ", conditionMessage(e)); NULL
    })
    if (!is.null(sr_clones)) write_csv(sr_clones, file.path(out_base, paste0(donor, ".sr_mgatk_clones.csv")))
  }

  if (isTRUE(donor_row$has_lr)) {
    lr_mt_bam <- file.path(project_dir, "results/mitochondrial_clones_longread", lr_sample, paste0(lr_sample, ".mt.bam"))
    if (file.exists(lr_mt_bam)) {
      lr_barcodes_file <- file.path(out_base, paste0(donor, "_barcodes_nosuffix.txt"))
      strip_barcode_suffix(donor_row$sr_barcodes_file, lr_barcodes_file)

      lr_mgatk_dir <- file.path(out_base, "lr_mgatk")
      lr_mgatk_final <- run_mgatk_tenx(lr_mt_bam, lr_barcodes_file, lr_mgatk_dir, donor, ncores, is_longread = TRUE)
      lr_signac <- if (!is.na(lr_mgatk_final)) mgatk_to_signac(lr_mgatk_final) else NA_character_
      if (!is.na(lr_signac)) {
        lr_clones <- tryCatch(call_mgatk_clones(lr_signac, n_clones = n_clones), error = function(e) {
          warning("LR call_mgatk_clones() failed for ", donor, ": ", conditionMessage(e)); NULL
        })
        if (!is.null(lr_clones)) write_csv(lr_clones, file.path(out_base, paste0(donor, ".lr_mgatk_clones.csv")))

        if (!is.na(sr_signac)) {
          merged <- tryCatch(call_mgatk_clones_merged(sr_signac, lr_signac, n_clones = n_clones), error = function(e) {
            warning("mgatk merge failed for ", donor, ": ", conditionMessage(e)); NULL
          })
          if (!is.null(merged)) write_csv(merged, file.path(out_base, paste0(donor, ".merged_mgatk_clones.csv")))
        }
      }
    } else {
      message("No long-read .mt.bam found at ", lr_mt_bam, " - skipping long-read mgatk for this donor")
    }
  }

  message("Donor ", donor, " mgatk stage done.")
}

# =================================================================
# CLI entry point
# =================================================================

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript 15a.mito_mgatk_stage.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv] [n_parallel_donors] [ncores_per_donor]")
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
          NCORES_PER_DONOR, " mgatk core(s) each.")

  invisible(parallel::mclapply(seq_len(nrow(donors)), function(i) {
    tryCatch(
      run_donor_mgatk(RUN_ID, RESULTS_DIR, PROJECT_DIR, donors[i, ], ncores = NCORES_PER_DONOR),
      error = function(e) message("FAILED for donor ", donors$donor[i], ": ", conditionMessage(e))
    )
  }, mc.cores = N_PARALLEL_DONORS, mc.preschedule = FALSE))

  message("\nmgatk stage complete for all donors.")
}
