#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Nextflow workflow - step 02
#
# What the pipeline actually produced, and whether the parameters it runs with
# are the ones the benchmark reports settled on.
#
# Two halves:
#   1. the published output tree per run - which stages exist, how many cells
#      survive each one, and how the two cell callers overlap
#   2. a side-by-side of every parameter that a benchmark report has an opinion
#      about, pipeline value against benchmark conclusion
#
# The second half is the point of this report. The pipeline is where the
# benchmarks are supposed to have landed, so any row where the two disagree is
# either a deliberate change nobody wrote down or a benchmark that never made
# it in.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

options(scipen = 999)

NF      <- "../nextflow"
RES_NF  <- "../results_nf"
out_dir <- "../results/benchmarking/nextflow_workflow"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- list.dirs(RES_NF, recursive = FALSE, full.names = FALSE)
runs <- runs[grepl("^[0-9]{6}_", runs)]

message("Runs found: ", paste(runs, collapse = ", "))

read_if <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- try(suppressMessages(read_csv(path, show_col_types = FALSE)), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

# ----------------------------
# 1. Stage inventory
# ----------------------------
# Shallow listing per stage - the Cell Ranger directories contain thousands of
# files and this is over a network mount.

stage_inventory <- map_dfr(runs, function(r) {
  stages <- list.dirs(file.path(RES_NF, r), recursive = FALSE, full.names = FALSE)
  stages <- stages[!grepl("^\\._", stages)]
  map_dfr(stages, function(s) {
    p <- file.path(RES_NF, r, s)
    f <- list.files(p, recursive = FALSE, full.names = TRUE)
    f <- f[!grepl("(^|/)\\._", f)]
    tibble(
      run = r, stage = s,
      n_entries = length(f),
      last_modified = if (length(f)) as.character(max(as.Date(file.info(f)$mtime))) else NA_character_
    )
  })
})

write_csv(stage_inventory, file.path(out_dir, "08_stage_inventory.csv"))

stage_matrix <- stage_inventory %>%
  filter(grepl("^[0-9]", stage)) %>%
  dplyr::select(run, stage, n_entries) %>%
  pivot_wider(names_from = run, values_from = n_entries)

write_csv(stage_matrix, file.path(out_dir, "09_stage_presence_by_run.csv"))

# ----------------------------
# 2. Cell count funnel
# ----------------------------

funnel <- map_dfr(runs, function(r) {
  d <- read_if(file.path(RES_NF, r, paste0(r, "_cell_count_funnel.csv")))
  if (is.null(d)) return(NULL)
  d$run <- r
  d
})

if (nrow(funnel)) {

  write_csv(funnel, file.path(out_dir, "10_cell_count_funnel.csv"))

  step_levels <- c("CellRanger", "EmptyDrops", "CellBender", "scDblFinder",
                   "HTO+Vireo demux", "Annotated (BoneMarrowMap)",
                   "Annotated & mapping QC pass")

  funnel_p <- funnel %>%
    mutate(step = factor(step, levels = step_levels)) %>%
    filter(!is.na(step))

  p_funnel <- ggplot(funnel_p, aes(step, n_cells, group = sample, colour = sample)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    facet_wrap(~ run, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 11) +
    labs(
      title = "Cells surviving each pipeline stage",
      subtitle = "EmptyDrops and CellBender are alternatives at the same point, not sequential",
      x = NULL, y = "Cells", colour = NULL
    ) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom")

  ggsave(file.path(out_dir, "01_cell_count_funnel.pdf"), p_funnel,
         width = 10, height = 9)
}

lost <- map_dfr(runs, function(r) {
  d <- read_if(file.path(RES_NF, r, paste0(r, "_cells_lost_by_step.csv")))
  if (is.null(d)) return(NULL)
  d$run <- r
  d
})

if (nrow(lost)) write_csv(lost, file.path(out_dir, "11_cells_lost_by_step.csv"))

# ----------------------------
# 3. The two cell callers
# ----------------------------

overlap <- map_dfr(runs, function(r) {
  d <- read_if(file.path(RES_NF, r, paste0(r, "_emptydrops_cellbender_overlap.csv")))
  if (is.null(d)) return(NULL)
  d$run <- r
  d
})

if (nrow(overlap)) {

  write_csv(overlap, file.path(out_dir, "12_emptydrops_cellbender_overlap.csv"))

  p_overlap <- overlap %>%
    dplyr::select(run, sample, n_emptydrops_only, n_common, n_cellbender_only) %>%
    pivot_longer(-c(run, sample), names_to = "set", values_to = "n_cells") %>%
    mutate(set = recode(set,
                        n_emptydrops_only = "EmptyDrops only",
                        n_common = "both",
                        n_cellbender_only = "CellBender only"),
           set = factor(set, levels = c("EmptyDrops only", "both", "CellBender only"))) %>%
    ggplot(aes(sample, n_cells, fill = set)) +
    geom_col(width = 0.7) +
    facet_wrap(~ run, scales = "free_x") +
    scale_fill_manual(values = c("EmptyDrops only" = "#1B9E77",
                                 "both" = "grey70",
                                 "CellBender only" = "#E31A1C")) +
    theme_bw(base_size = 11) +
    labs(title = "Barcodes called by each cell caller",
         subtitle = "CellBender barcodes are the pipeline's primary cell calls; EmptyDrops feeds the DSB background",
         x = NULL, y = "Barcodes", fill = NULL) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom")

  ggsave(file.path(out_dir, "02_cell_caller_overlap.pdf"), p_overlap,
         width = 10, height = 5)
}

# ----------------------------
# 4. Pipeline parameters against the benchmark conclusions
# ----------------------------
# The benchmark column is what the corresponding lab archive report concluded.
# Rows are marked "matches", "differs" or "not benchmarked".

params_tbl <- read_if(file.path(out_dir, "03_pipeline_params.csv"))

get_param <- function(p) {
  if (is.null(params_tbl)) return(NA_character_)
  v <- params_tbl$value[params_tbl$param == p]
  if (!length(v)) NA_character_ else v[1]
}

comparison <- tibble::tribble(
  ~setting, ~pipeline_value, ~benchmark_conclusion, ~report, ~verdict,

  "emptyDrops FDR", get_param("emptydrops_fdr"),
  "FDR <= 0.01", "04 cell calling", "matches",

  "emptyDrops minimum UMI", get_param("emptydrops_umi_min"),
  "100", "04 cell calling", "matches",

  "primary cell caller", "CellBender (hardcoded in main.nf)",
  "CellBender for cell calls, emptyDrops retained for DSB background",
  "04 cell calling", "matches the DAG, contradicts params.cell_caller",

  "params.cell_caller", get_param("cell_caller"),
  "n/a - the parameter is never read", "04 cell calling", "dead parameter",

  "HTO positive.quantile", get_param("hto_positive_quantile"),
  "0.99", "03 demultiplexing", "matches",

  "demultiplexing method", "cellSNP-lite + vireo, compared against HTO",
  "genotyping outperforms HTO hashing", "03 demultiplexing", "matches",

  "cellSNP min count", get_param("cellsnp_min_count"),
  "20", "03 demultiplexing", "matches",

  "cellSNP min MAF", get_param("cellsnp_min_maf"),
  "0.10", "03 demultiplexing", "matches",

  "ADT normalisation", "DSB, in CITE_QC",
  "DSB outperforms CLR", "05 CITE-seq", "matches",

  "BoneMarrowMap MAD_threshold", "4 (bin/annotate.R)",
  "2.5 - the default, and what the standalone scripts used",
  "07 cell type annotation", "differs",

  "Numbat container", get_param("numbat_sif"),
  "version not recorded - `latest` tag", "08 calling cancer cells", "unpinned",

  "CellBender container", get_param("cellbender_sif"),
  "CellBender 0.3.0", "04 cell calling", "matches",

  "transcriptome", get_param("transcriptome"),
  "refdata-gex-GRCh38-2024-A", "01 LK1 pilot", "matches"
)

write_csv(comparison, file.path(out_dir, "13_params_vs_benchmarks.csv"))

# ----------------------------
# 5. Which run is which
# ----------------------------
# The pipeline runs the RE-SEQUENCED LK1 (260522), not the pilot (260423) that
# reports 01, 07 and 08 analysed. Worth stating explicitly, because the two
# output trees look interchangeable and are not.

samples <- read_if(file.path(NF, "conf", "samples.csv"))

run_map <- tibble::tribble(
  ~run, ~short, ~note,
  "260423_VH01624_453_222HWMYNX", "LK1 pilot",
  "NOT in the pipeline samplesheet. Superseded by the re-sequencing; analysed only by the standalone scripts in reports 01, 07 and 08.",
  "260522_VH01624_461_222JLJVNX", "LK1 re-sequenced",
  "What the pipeline runs for LK1.",
  "260528_VH01624_464_222K7VKNX", "LK2",
  "Run by both the pipeline and the standalone scripts.",
  "260717_VH01624_477_222KG22NX", "LK3",
  "Pipeline only. Carries the unhashed donor WEI21_26-17_NK."
)

write_csv(run_map, file.path(out_dir, "14_run_identity.csv"))

# ----------------------------
# Report
# ----------------------------

cat("\n=== Stages present per run ===\n")
print(as.data.frame(stage_matrix))

cat("\n=== Cell count funnel ===\n")
if (nrow(funnel)) print(as.data.frame(funnel))

cat("\n=== Cell caller overlap ===\n")
if (nrow(overlap)) print(as.data.frame(overlap))

cat("\n=== Params against benchmarks ===\n")
print(as.data.frame(comparison[, c("setting", "pipeline_value", "verdict")]))

cat("\n=== Samplesheet ===\n")
if (!is.null(samples)) print(as.data.frame(samples))

cat("\nDone. Output directory:", out_dir, "\n")
