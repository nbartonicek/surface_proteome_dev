#!/usr/bin/env Rscript

# Batch driver for 22.patient_qc_report.Rmd - renders one HTML report per
# donor for every GEX sample found in a given run, instead of one at a time
# by hand. Does not modify 22.patient_qc_report.Rmd itself; this only loops
# rmarkdown::render() over it, matching the approach the Rmd's own trailing
# notes describe ("loop over the donor list ... call rmarkdown::render() in
# a for loop with a different output_file each time").
#
# Sample/donor discovery: GEX sample names come from the subdirectories of
# results_nf/<run_id>/08_annotate/ (one per pooled GEM well). Within each
# sample, donor_ids come from the "sample_name" column of that sample's
# composition_broad_by_sample_name.csv - confirmed against real data to be
# exactly the same string values used as `donor_id` throughout this project
# (e.g. "HBDN498-TP53", "APOP576-TP53", "normal-01"), not derived/guessed.
#
# Each donor's render is wrapped in tryCatch so one failure (e.g. a donor
# with no usable Seurat object, same failure mode diagnosed earlier for
# LK1-GEX/HBDN501-KMT2A) does not abort the whole batch - failures are
# logged to the summary CSV instead, with the actual error message, so they
# can be triaged after the fact rather than losing every other donor's
# report because of one bad one.

suppressPackageStartupMessages({
  library(rmarkdown)
  library(dplyr)
  library(readr)
})

# ----------------------------
# Parameters
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript 22a.patient_qc_report_batch.R <run_id> [results_dir] [output_dir] [exclude_donors_csv] [force]\n",
    "  run_id:              e.g. 260528_VH01624_464_222K7VKNX\n",
    "  results_dir:         default '../results_nf'\n",
    "  output_dir:          where HTMLs land, default '../results/patient_reports/<run_id>'\n",
    "  exclude_donors_csv:  comma-separated donor_ids to skip, e.g. 'normal-01' (default: none)\n",
    "  force:               'true' to re-render donors whose HTML already exists (default: false, skip existing)\n",
    "Example:\n",
    "  Rscript 22a.patient_qc_report_batch.R 260528_VH01624_464_222K7VKNX"
  )
}

RUN_ID       <- args[1]
RESULTS_DIR  <- if (length(args) >= 2) args[2] else "../results_nf"
OUTPUT_DIR   <- if (length(args) >= 3) args[3] else file.path("../results/patient_reports", RUN_ID)
EXCLUDE_DONORS <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else character(0)
FORCE_RERUN  <- length(args) >= 5 && tolower(args[5]) == "true"

RMD_FILE <- "22.patient_qc_report.Rmd"
DE_OUT_DIR <- "../results/patient_reports"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(RMD_FILE))

cat("Run:", RUN_ID, "\n")
cat("Results dir:", RESULTS_DIR, "\n")
cat("Output dir:", OUTPUT_DIR, "\n")
if (length(EXCLUDE_DONORS) > 0) cat("Excluding donors:", paste(EXCLUDE_DONORS, collapse = ", "), "\n")

# ----------------------------
# 1. Discover GEX samples in this run
# ----------------------------

annotate_dir <- file.path(RESULTS_DIR, RUN_ID, "08_annotate")

if (!dir.exists(annotate_dir)) {
  stop("No 08_annotate directory found for this run: ", annotate_dir)
}

samples <- list.dirs(annotate_dir, full.names = FALSE, recursive = FALSE)
samples <- samples[nzchar(samples)]

if (length(samples) == 0) {
  stop("No GEX sample subdirectories found under: ", annotate_dir)
}

cat("Samples found:", paste(samples, collapse = ", "), "\n")

# ----------------------------
# 2. Discover donors per sample, from that sample's composition table
# ----------------------------

find_donors_for_sample <- function(sample_name) {
  comp_path <- file.path(annotate_dir, sample_name, "tables", "composition_broad_by_sample_name.csv")
  if (!file.exists(comp_path)) {
    warning("No composition_broad_by_sample_name.csv for sample ", sample_name, " - skipping this sample.")
    return(character(0))
  }
  comp <- read_csv(comp_path, show_col_types = FALSE)
  if (!"sample_name" %in% names(comp)) {
    warning("composition table for ", sample_name, " has no sample_name column - skipping this sample.")
    return(character(0))
  }
  sort(unique(comp$sample_name))
}

jobs <- bind_rows(lapply(samples, function(s) {
  donors <- find_donors_for_sample(s)
  donors <- setdiff(donors, EXCLUDE_DONORS)
  if (length(donors) == 0) return(NULL)
  tibble(sample_name = s, donor_id = donors)
}))

if (is.null(jobs) || nrow(jobs) == 0) {
  stop("No (sample, donor) pairs discovered for run ", RUN_ID, " - nothing to render.")
}

cat("\nDonor/sample pairs to render (", nrow(jobs), "):\n", sep = "")
print(jobs)

# ----------------------------
# 3. Render one HTML per donor, tolerating individual failures
# ----------------------------

render_one <- function(sample_name, donor_id) {
  out_file <- paste0(donor_id, "_patient_qc_report.html")
  out_path <- file.path(OUTPUT_DIR, out_file)

  if (file.exists(out_path) && !FORCE_RERUN) {
    cat("\n[SKIP] ", donor_id, " - output already exists: ", out_path, "\n", sep = "")
    return(list(sample_name = sample_name, donor_id = donor_id, status = "skipped", error = NA_character_))
  }

  cat("\n[RENDER] ", sample_name, " / ", donor_id, " -> ", out_path, "\n", sep = "")

  result <- tryCatch({
    rmarkdown::render(
      RMD_FILE,
      params = list(
        run_id        = RUN_ID,
        sample_name   = sample_name,
        donor_id      = donor_id,
        results_dir   = RESULTS_DIR,
        seurat_path   = NA,
        min_cells_group = 10,
        min_cells_de  = 20,
        bcv           = 0.4,
        fdr_cutoff    = 0.05,
        logfc_cutoff  = 1,
        de_out_dir    = DE_OUT_DIR
      ),
      output_file = out_file,
      output_dir  = OUTPUT_DIR,
      envir       = new.env()   # isolate each render's globals from the next
    )
    list(sample_name = sample_name, donor_id = donor_id, status = "success", error = NA_character_)
  }, error = function(e) {
    cat("[FAILED] ", donor_id, ": ", conditionMessage(e), "\n", sep = "")
    list(sample_name = sample_name, donor_id = donor_id, status = "failed", error = conditionMessage(e))
  })

  result
}

summary_rows <- Map(render_one, jobs$sample_name, jobs$donor_id)
summary_df <- bind_rows(summary_rows)

summary_path <- file.path(OUTPUT_DIR, paste0(RUN_ID, "_batch_render_summary.csv"))
write_csv(summary_df, summary_path)

cat("\n========================================\n")
cat("Batch render complete for run", RUN_ID, "\n")
print(summary_df)
cat("Summary written to:", summary_path, "\n")

if (any(summary_df$status == "failed")) {
  cat("\nWARNING:", sum(summary_df$status == "failed"), "donor(s) failed - see summary CSV for error messages.\n")
}
