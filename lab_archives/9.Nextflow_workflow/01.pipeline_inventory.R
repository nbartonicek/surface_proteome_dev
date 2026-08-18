#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Nextflow workflow - step 01
#
# Reads the pipeline and writes down what it actually is: the processes, the
# script each one runs, where it publishes, the resources it asks for, and the
# parameters it is configured with.
#
# This report is a different shape from reports 1-8. There are no analysis
# scripts of mine to collate - the pipeline is the artefact, and it lives in a
# DIFFERENT repository (nbartonicek/surface_proteome) from the lab archives
# (nbartonicek/surface_proteome_dev). So this step reads the pipeline rather
# than rewriting it, and nothing here modifies anything under nextflow/.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

options(scipen = 999)

NF      <- "../nextflow"
out_dir <- "../results/benchmarking/nextflow_workflow"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(NF)) stop("Cannot find the pipeline at ", NF)

read_lines_safe <- function(p) if (file.exists(p)) readLines(p, warn = FALSE) else character(0)

# ----------------------------
# 1. Processes, from modules/
# ----------------------------

module_files <- list.files(file.path(NF, "modules"), pattern = "\\.nf$", full.names = TRUE)
module_files <- module_files[!grepl("(^|/)\\._", module_files)]

parse_module <- function(f) {

  txt <- read_lines_safe(f)
  if (!length(txt)) return(NULL)

  process <- str_match(txt, "^\\s*process\\s+([A-Za-z0-9_]+)\\s*\\{")[, 2]
  process <- process[!is.na(process)]

  # the command actually invoked: first non-comment line inside script:"""
  script_start <- grep('script:', txt)
  cmd <- NA_character_
  if (length(script_start)) {
    tail_txt <- txt[(script_start[1] + 1):length(txt)]
    cand <- tail_txt[!grepl('^\\s*("""|\\}|//|$)', tail_txt)]
    if (length(cand)) cmd <- str_trim(str_remove(cand[1], "\\\\\\s*$"))
  }

  publish <- str_match_all(txt, 'publishDir\\s*\\{?\\s*"([^"]+)"')
  publish <- unlist(lapply(publish, function(m) if (nrow(m)) m[, 2] else character(0)))
  publish <- unique(str_replace_all(publish, "\\$\\{params\\.outdir\\}/", ""))

  emits <- str_match_all(txt, "emit:\\s*([A-Za-z0-9_]+)")
  emits <- unlist(lapply(emits, function(m) if (nrow(m)) m[, 2] else character(0)))

  tibble(
    module   = basename(f),
    process  = paste(process, collapse = ", "),
    command  = cmd %||% NA_character_,
    publishes = paste(unique(publish), collapse = "; "),
    emits    = paste(unique(emits), collapse = ", ")
  )
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

processes <- map_dfr(module_files, parse_module)

write_csv(processes, file.path(out_dir, "01_pipeline_processes.csv"))

# ----------------------------
# 2. Workflow order, from main.nf
# ----------------------------
# The numbered section comments in main.nf are the author's own ordering, so
# they are used rather than a topological sort of the channel graph.

main <- read_lines_safe(file.path(NF, "main.nf"))

sections <- tibble(line = seq_along(main), txt = main) %>%
  filter(str_detect(txt, "^\\s*//\\s*[0-9]+[a-z]?\\.\\s")) %>%
  mutate(
    stage = str_trim(str_remove(txt, "^\\s*//\\s*")),
    order = row_number()
  ) %>%
  dplyr::select(order, line, stage)

calls <- tibble(line = seq_along(main), txt = main) %>%
  filter(str_detect(txt, "^\\s*[A-Z][A-Z0-9_]+\\s*\\(")) %>%
  mutate(process = str_match(txt, "^\\s*([A-Z][A-Z0-9_]+)\\s*\\(")[, 2]) %>%
  dplyr::select(line, process)

stage_of <- function(l) {
  s <- sections$stage[sections$line < l]
  if (!length(s)) NA_character_ else tail(s, 1)
}

workflow_order <- calls %>%
  mutate(stage = map_chr(line, stage_of)) %>%
  distinct(process, .keep_all = TRUE) %>%
  dplyr::select(stage, process)

write_csv(workflow_order, file.path(out_dir, "02_workflow_order.csv"))

# ----------------------------
# 3. Parameters
# ----------------------------

cfg <- read_lines_safe(file.path(NF, "nextflow.config"))

param_block <- {
  s <- grep("^\\s*params\\s*\\{", cfg)[1]
  if (is.na(s)) character(0) else {
    depth <- 0; out <- character(0)
    for (i in s:length(cfg)) {
      depth <- depth + str_count(cfg[i], "\\{") - str_count(cfg[i], "\\}")
      out <- c(out, cfg[i])
      if (depth == 0 && i > s) break
    }
    out
  }
}

params_tbl <- tibble(txt = param_block) %>%
  filter(str_detect(txt, "^\\s*[A-Za-z_][A-Za-z0-9_]*\\s*=")) %>%
  mutate(
    param = str_trim(str_extract(txt, "^\\s*[A-Za-z_][A-Za-z0-9_]*")),
    value = str_trim(str_remove(txt, "^\\s*[A-Za-z_][A-Za-z0-9_]*\\s*=\\s*"))
  ) %>%
  dplyr::select(param, value)

write_csv(params_tbl, file.path(out_dir, "03_pipeline_params.csv"))

# ----------------------------
# 4. Resources per process
# ----------------------------

res_tbl <- {
  idx <- grep("withName:", cfg)
  map_dfr(idx, function(i) {
    name <- str_match(cfg[i], "withName:\\s*'?\"?([A-Za-z0-9_|]+)'?\"?")[, 2]
    depth <- 0; blk <- character(0)
    for (j in i:length(cfg)) {
      depth <- depth + str_count(cfg[j], "\\{") - str_count(cfg[j], "\\}")
      blk <- c(blk, cfg[j])
      if (depth == 0 && j > i) break
    }
    grab <- function(k) {
      hit <- str_match(blk, paste0("^\\s*", k, "\\s*=\\s*'?\"?([^'\"\n]+)'?\"?\\s*$"))[, 2]
      hit <- hit[!is.na(hit)]
      if (!length(hit)) NA_character_ else str_trim(hit[1])
    }
    tibble(process = name, cpus = grab("cpus"), memory = grab("memory"), time = grab("time"))
  })
}

write_csv(res_tbl, file.path(out_dir, "04_process_resources.csv"))

# ----------------------------
# 5. Samplesheet
# ----------------------------

samples <- suppressMessages(read_csv(file.path(NF, "conf", "samples.csv"),
                                     show_col_types = FALSE))
write_csv(samples, file.path(out_dir, "05_samplesheet.csv"))

# ----------------------------
# 6. Things that look wrong
# ----------------------------
# Recorded, not fixed - the pipeline is a separate repository and is not this
# report's to edit.

issues <- tibble::tribble(
  ~issue, ~where, ~detail,

  "params.cell_caller is declared and never used",
  "nextflow.config",
  paste("It is set to \"emptydrops\", but the workflow assigns",
        "cell_barcodes = CELLBENDER.out.cellbender_barcodes unconditionally.",
        "The parameter says the opposite of what the DAG does, and changing it",
        "would have no effect."),

  "MAD_threshold differs between the two annotation code paths",
  "bin/annotate.R vs bin/lib_early_qc_annotation.R",
  paste("annotate.R, the live ANNOTATE process, uses 4.",
        "lib_early_qc_annotation.R uses the BoneMarrowMap default of 2.5.",
        "The standalone benchmark scripts also used 2.5, so results_nf and",
        "results/seurat_annotated are not directly comparable."),

  "leftover debug logging in the workflow",
  "main.nf",
  paste("A log.info marked \"TEMP DEBUG - remove once NUMBAT_PILEUP null-path",
        "bug is confirmed fixed\" still fires for every donor of every run."),

  "donor names differ between the pipeline and the standalone analyses",
  "conf/samples.csv",
  paste("samples.csv has HBDN206-MN-pCT, HBDN392-MDS, HBDN501-KMT2A;",
        "the standalone scripts and results use HBDN206-MNpCT,",
        "HBDN392-AML-MDS, HBDN501-AML-KMT2A. Same donors, different strings,",
        "so joins between the two output trees will not match on donor."),

  "the Numbat container is pinned to `latest`",
  "nextflow.config",
  paste("numbat_sif points at numbat-rbase_latest.sif. The version that",
        "produced the published results is not recoverable from the image name."),

  "conda environments are absolute paths in one home directory",
  "nextflow.config",
  paste("r_env, cite_seq_count_env, cellsnp_env and cellranger_bin are all",
        "under /home/nbartonicek, so the pipeline only runs as this user."),

  "DOROTHEA_VIPER is a dead module",
  "modules/dorothea_viper.nf, bin/dorothea_viper.R",
  paste("The module and its R script exist, but the process is neither included",
        "nor called anywhere in main.nf, so it never runs. It is the only module",
        "of the 18 that is not wired into the workflow.")
)

# Cross-check the dead-module claim rather than asserting it
included <- str_match_all(paste(main, collapse = "\n"),
                          "include\\s*\\{([^}]*)\\}")[[1]][, 2]
included <- str_trim(unlist(str_split(paste(included, collapse = ";"), "[;\n]")))
included <- included[nzchar(included)]

declared <- unlist(str_split(processes$process, ",\\s*"))
declared <- declared[nzchar(declared)]

never_included <- setdiff(declared, included)
if (length(never_included)) {
  message("Processes declared in modules/ but never included in main.nf: ",
          paste(never_included, collapse = ", "))
}

write_csv(tibble(process = declared,
                 included_in_main = declared %in% included),
          file.path(out_dir, "07_process_wired_into_workflow.csv"))

write_csv(issues, file.path(out_dir, "06_pipeline_observations.csv"))

# ----------------------------
# Report
# ----------------------------

cat("\n=== Processes (", nrow(processes), ") ===\n", sep = "")
print(as.data.frame(processes[, c("process", "command")]))

cat("\n=== Workflow order ===\n")
print(as.data.frame(workflow_order))

cat("\n=== Params ===\n")
print(as.data.frame(params_tbl))

cat("\n=== Samplesheet ===\n")
print(as.data.frame(samples))

cat("\n=== Observations ===\n")
print(as.data.frame(issues[, c("issue", "where")]))

cat("\nDone. Output directory:", out_dir, "\n")
