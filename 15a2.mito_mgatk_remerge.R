#!/usr/bin/env Rscript

# Recovery utility - re-runs ONLY the SR+LR merge step of
# 15a.mito_mgatk_stage.R, for donors whose sr_mgatk/final/*.signac.rds and
# lr_mgatk/final/*.signac.rds already exist on disk (the expensive mgatk
# tenx step - the part that took >12h for the run that finished
# 2026-07-28). Does NOT re-run mgatk tenx at all - just rebuilds each
# donor's merged clone calls in memory from the existing .signac.rds files,
# using the fixed align_to() in 15a.mito_mgatk_stage.R (the old version hit
# a Matrix-package `[<-` class-dispatch bug and failed identically for
# every donor in that run).
#
# Run under the same `mgatk` env as 15a.mito_mgatk_stage.sh (this sources
# that file for its function defs, which needs Signac/Matrix - even though
# the merge step itself doesn't call Signac's ReadMGATK()).
#
# Usage:
#   Rscript 15a2.mito_mgatk_remerge.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv]

source("15a.mito_mgatk_stage.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 15a2.mito_mgatk_remerge.R <run_id> [sr_sample] [results_dir] [exclude_donors_csv]")
}

RUN_ID <- args[1]
SR_SAMPLE <- if (length(args) >= 2) args[2] else "LK2-GEX"
RESULTS_DIR <- if (length(args) >= 3) args[3] else "/scratch/users/nbartonicek/projects/amgen/results_nf"
PROJECT_DIR <- dirname(RESULTS_DIR)
EXCLUDE_DONORS <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else character(0)

donors <- discover_donors(RUN_ID, RESULTS_DIR, SR_SAMPLE) %>% filter(!donor %in% EXCLUDE_DONORS)
message("Donors found: ", paste(donors$donor, collapse = ", "))

for (i in seq_len(nrow(donors))) {
  donor <- donors$donor[i]
  out_base <- donor_out_dir(PROJECT_DIR, RUN_ID, donor)
  sr_signac <- file.path(out_base, "sr_mgatk", "final", paste0(donor, ".signac.rds"))
  lr_signac <- file.path(out_base, "lr_mgatk", "final", paste0(donor, ".signac.rds"))

  if (!file.exists(sr_signac) || !file.exists(lr_signac)) {
    message(donor, ": missing sr/lr signac.rds under ", out_base, " - skipping (rerun 15a for this donor first)")
    next
  }

  message("\n========== [remerge] Donor: ", donor, " ==========")
  merged <- tryCatch(
    call_mgatk_clones_merged(sr_signac, lr_signac),
    error = function(e) { message("FAILED for donor ", donor, ": ", conditionMessage(e)); NULL }
  )
  if (!is.null(merged)) {
    write_csv(merged, file.path(out_base, paste0(donor, ".merged_mgatk_clones.csv")))
    message(donor, ": wrote merged_mgatk_clones.csv")
  }
}

message("\nRemerge complete.")
