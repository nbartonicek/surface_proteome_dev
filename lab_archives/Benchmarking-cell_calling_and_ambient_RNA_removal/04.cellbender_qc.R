#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Cell calling and ambient RNA removal - step 04 of 7
#
# QC after CellBender: knee plot, UMI before/after (density and ECDF), percent background removed, raw-vs-corrected scatter, summary CSV.
#
# Frozen for the lab archive 2026-08-04 from scripts/2a.cellbender.R (mtime 2026-06-26).
# md5 of the original: 3a705cec70dad29382b616da2545dbd2
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# =============================================================================
# cellbender_qc_plots.R
# QC plots after running CellBender on Cell Ranger raw output.
#
# Produces:
#   1. Knee / barcode rank plot (raw)
#   2. UMI distribution: before vs after CellBender (density + ECDF)
#   3. % background RNA removed per cell (histogram + summary)
#   4. Cell-level scatter: raw UMI vs CellBender UMI
#   5. Summary stats table (CSV)
#
# Usage (called from run_cellbender.sh):
#   Rscript cellbender_qc_plots.R \
#       --sample SAMPLE_ID \
#       --raw    /path/to/raw_feature_bc_matrix \
#       --cb-h5  /path/to/cellbender_output.h5 \
#       --outdir /path/to/plots \
#       --fdr    0.01
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Prefer DropletUtils for reading 10x; fall back gracefully
has_dropletutils <- requireNamespace("DropletUtils", quietly = TRUE)
has_rhdf5        <- requireNamespace("rhdf5",        quietly = TRUE)

# -------------------------------------------------------
# Parse arguments
# -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}

sample_id <- parse_arg("--sample")
raw_dir   <- parse_arg("--raw")
cb_h5     <- parse_arg("--cb-h5")
out_dir   <- parse_arg("--outdir")
fdr_thr   <- as.numeric(parse_arg("--fdr", "0.01"))

stopifnot(
  !is.null(sample_id),
  !is.null(raw_dir),
  !is.null(cb_h5),
  !is.null(out_dir)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:  ", sample_id)
message("Raw dir: ", raw_dir)
message("CB H5:   ", cb_h5)
message("Out dir: ", out_dir)

# -------------------------------------------------------
# Helper: read raw Cell Ranger matrix
# -------------------------------------------------------
read_raw_matrix <- function(raw_dir) {
  
  # Support both flat and nested (Gene Expression) layouts
  candidates <- c(
    raw_dir,
    file.path(raw_dir, "Gene Expression")
  )
  
  valid <- candidates[sapply(candidates, function(d) {
    file.exists(file.path(d, "matrix.mtx.gz")) ||
      file.exists(file.path(d, "matrix.mtx"))
  })]
  
  if (length(valid) == 0) stop("Cannot find matrix.mtx[.gz] under: ", raw_dir)
  
  mat_dir <- valid[1]
  
  if (has_dropletutils) {
    DropletUtils::read10xCounts(mat_dir, col.names = TRUE)
  } else {
    # Manual fallback
    mtx_file <- if (file.exists(file.path(mat_dir, "matrix.mtx.gz"))) {
      file.path(mat_dir, "matrix.mtx.gz")
    } else {
      file.path(mat_dir, "matrix.mtx")
    }
    bc_file <- if (file.exists(file.path(mat_dir, "barcodes.tsv.gz"))) {
      file.path(mat_dir, "barcodes.tsv.gz")
    } else {
      file.path(mat_dir, "barcodes.tsv")
    }
    mat <- readMM(mtx_file)
    barcodes <- readLines(bc_file)
    colnames(mat) <- barcodes
    mat
  }
}

# -------------------------------------------------------
# Helper: read CellBender H5 output
# -------------------------------------------------------
read_cellbender_h5 <- function(h5_path) {
  
  if (!file.exists(h5_path)) {
    stop("CellBender H5 not found: ", h5_path)
  }
  
  if (!has_rhdf5) stop("rhdf5 package required to read CellBender H5 output.")
  
  message("Reading CellBender H5: ", h5_path)
  
  h5 <- rhdf5::h5dump(h5_path, load = FALSE)
  
  # CellBender writes under /matrix or /background_fraction
  counts_path <- "/matrix/data"
  indices_path <- "/matrix/indices"
  indptr_path  <- "/matrix/indptr"
  shape_path   <- "/matrix/shape"
  barcodes_path <- "/matrix/barcodes"
  
  data    <- rhdf5::h5read(h5_path, counts_path)
  indices <- rhdf5::h5read(h5_path, indices_path)
  indptr  <- rhdf5::h5read(h5_path, indptr_path)
  shape   <- rhdf5::h5read(h5_path, shape_path)
  barcodes <- rhdf5::h5read(h5_path, barcodes_path)
  
  mat <- sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = as.numeric(data),
    dims = as.integer(shape),
    dimnames = list(NULL, as.character(barcodes))
  )
  
  # Also try to read background_fraction if available
  bg_frac <- tryCatch(
    rhdf5::h5read(h5_path, "/background_fraction"),
    error = function(e) NULL
  )
  
  list(mat = mat, bg_frac = bg_frac, barcodes = barcodes)
}

# -------------------------------------------------------
# Helper: match barcodes between raw and CB matrices
# -------------------------------------------------------
strip_suffix <- function(barcodes) sub("-[0-9]+$", "", barcodes)

# -------------------------------------------------------
# Load data
# -------------------------------------------------------
message("\nLoading raw Cell Ranger matrix...")

raw_obj <- read_raw_matrix(raw_dir)

if (has_dropletutils && is(raw_obj, "SingleCellExperiment")) {
  raw_mat  <- counts(raw_obj)
  raw_bcs  <- colnames(raw_mat)
} else {
  raw_mat  <- raw_obj
  raw_bcs  <- colnames(raw_mat)
}

message("  Droplets: ", ncol(raw_mat), "  Features: ", nrow(raw_mat))

message("\nLoading CellBender output...")
cb <- read_cellbender_h5(cb_h5)
cb_mat  <- cb$mat
cb_bcs  <- strip_suffix(cb$barcodes)
bg_frac <- cb$bg_frac

# Also try filtered H5
cb_h5_filtered <- sub("\\.h5$", "_filtered.h5", cb_h5)
cb_filtered <- if (file.exists(cb_h5_filtered)) {
  message("Reading filtered H5: ", cb_h5_filtered)
  read_cellbender_h5(cb_h5_filtered)
} else {
  NULL
}

message("  CB droplets: ", ncol(cb_mat))

# -------------------------------------------------------
# UMI totals
# -------------------------------------------------------
raw_umi   <- colSums(raw_mat)
raw_umi   <- sort(raw_umi, decreasing = TRUE)

# Intersect barcodes for before/after comparison
raw_bcs_stripped  <- strip_suffix(raw_bcs)
common_bcs_raw_idx <- match(cb_bcs, raw_bcs_stripped)
valid_idx <- !is.na(common_bcs_raw_idx)

cb_umi_matched  <- colSums(cb_mat[, valid_idx, drop = FALSE])
raw_umi_matched <- raw_umi[common_bcs_raw_idx[valid_idx]]

# -------------------------------------------------------
# % background RNA removed
# -------------------------------------------------------
pct_bg_removed <- 100 * (raw_umi_matched - cb_umi_matched) / (raw_umi_matched + 1e-6)
pct_bg_removed[pct_bg_removed < 0] <- 0  # clamp floating point negatives

# If CellBender wrote background_fraction field, use it directly for cells
if (!is.null(bg_frac)) {
  message("  background_fraction field found in H5 — will use for per-cell plot")
  bg_frac_pct <- bg_frac * 100
} else {
  bg_frac_pct <- NULL
}

# -------------------------------------------------------
# Summary stats
# -------------------------------------------------------
summary_stats <- tibble(
  metric = c(
    "Total raw droplets",
    "Total features",
    "CB output droplets",
    "Median raw UMI (CB barcodes)",
    "Median CB UMI (CB barcodes)",
    "Median % background removed",
    "Mean % background removed",
    "Cells with > 20% background removed"
  ),
  value = c(
    ncol(raw_mat),
    nrow(raw_mat),
    ncol(cb_mat),
    round(median(raw_umi_matched)),
    round(median(cb_umi_matched)),
    round(median(pct_bg_removed), 2),
    round(mean(pct_bg_removed), 2),
    sum(pct_bg_removed > 20)
  )
)

write.csv(
  summary_stats,
  file.path(out_dir, paste0(sample_id, "_cellbender_qc_summary.csv")),
  row.names = FALSE
)

message("\nSummary stats:")
print(summary_stats, n = Inf)

# -------------------------------------------------------
# Plot theme
# -------------------------------------------------------
theme_qc <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92")
    )
}

# -------------------------------------------------------
# Plot 1: Knee / barcode rank plot (raw)
# -------------------------------------------------------
message("\nPlot 1: Knee plot...")

knee_df <- tibble(
  rank = seq_along(raw_umi),
  umi  = raw_umi
)

# Knee point: largest second derivative of log-log curve
log_rank <- log10(knee_df$rank)
log_umi  <- log10(knee_df$umi + 1)
d2       <- diff(diff(log_umi) / diff(log_rank)) / diff(log_rank[-1])
knee_rank <- which.min(d2)

p_knee <- ggplot(knee_df, aes(rank, umi)) +
  geom_line(colour = "#2166ac", linewidth = 0.7) +
  geom_vline(
    xintercept = knee_rank,
    linetype = "dashed", colour = "#d6604d", linewidth = 0.7
  ) +
  annotate(
    "text",
    x = knee_rank * 1.5,
    y = max(raw_umi) * 0.5,
    label = paste0("Knee ~rank ", comma(knee_rank)),
    colour = "#d6604d", size = 3.2, hjust = 0
  ) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(
    title    = paste0(sample_id, " — Barcode rank plot (raw)"),
    subtitle = paste0(comma(ncol(raw_mat)), " total droplets"),
    x = "Barcode rank",
    y = "Total UMI count"
  ) +
  theme_qc()

# -------------------------------------------------------
# Plot 2a: UMI distribution before vs after (density)
# -------------------------------------------------------
message("Plot 2: UMI distributions...")

umi_df <- bind_rows(
  tibble(umi = raw_umi_matched, source = "Raw (pre-CellBender)"),
  tibble(umi = cb_umi_matched,  source = "CellBender output")
) %>%
  filter(umi > 0)

p_umi_density <- ggplot(umi_df, aes(umi, colour = source, fill = source)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = c("Raw (pre-CellBender)" = "#4393c3",
                                 "CellBender output"    = "#d6604d")) +
  scale_fill_manual(values  = c("Raw (pre-CellBender)" = "#4393c3",
                                "CellBender output"    = "#d6604d")) +
  labs(
    title    = "UMI distribution: raw vs CellBender",
    subtitle = paste0(comma(length(cb_umi_matched)), " matched barcodes"),
    x = "Total UMI (log10)",
    y = "Density",
    colour = NULL, fill = NULL
  ) +
  theme_qc() +
  theme(legend.position = "top")

# Plot 2b: ECDF
p_umi_ecdf <- ggplot(umi_df, aes(umi, colour = source)) +
  stat_ecdf(linewidth = 0.7) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = c("Raw (pre-CellBender)" = "#4393c3",
                                 "CellBender output"    = "#d6604d")) +
  labs(
    title  = "ECDF: raw vs CellBender UMI",
    x = "Total UMI (log10)",
    y = "Cumulative fraction",
    colour = NULL
  ) +
  theme_qc() +
  theme(legend.position = "top")

# -------------------------------------------------------
# Plot 3: % background RNA removed
# -------------------------------------------------------
message("Plot 3: % background removed...")

bg_df <- tibble(pct_bg = pct_bg_removed)

med_bg <- median(pct_bg_removed)
mn_bg  <- mean(pct_bg_removed)

p_bg_hist <- ggplot(bg_df, aes(pct_bg)) +
  geom_histogram(
    binwidth = 2, fill = "#74add1", colour = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = med_bg, colour = "#d73027", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = mn_bg,  colour = "#fc8d59", linetype = "dotted",  linewidth = 0.8) +
  annotate("text", x = med_bg + 0.5, y = Inf,
           label = paste0("Median: ", round(med_bg, 1), "%"),
           vjust = 2, hjust = 0, colour = "#d73027", size = 3.2) +
  annotate("text", x = mn_bg + 0.5, y = Inf,
           label = paste0("Mean: ", round(mn_bg, 1), "%"),
           vjust = 4, hjust = 0, colour = "#fc8d59", size = 3.2) +
  labs(
    title    = "% background RNA removed per barcode",
    subtitle = paste0(
      comma(sum(pct_bg_removed > 20)), " barcodes (>20% removed); ",
      comma(length(pct_bg_removed)), " total"
    ),
    x = "% UMI attributed to background (removed)",
    y = "Number of barcodes"
  ) +
  theme_qc()

# -------------------------------------------------------
# Plot 4: raw vs CB UMI scatter (hex)
# -------------------------------------------------------
message("Plot 4: raw vs CB scatter...")

scatter_df <- tibble(
  raw_umi = raw_umi_matched,
  cb_umi  = cb_umi_matched,
  pct_bg  = pct_bg_removed
)

p_scatter <- ggplot(scatter_df, aes(raw_umi + 1, cb_umi + 1, colour = pct_bg)) +
  geom_hex(bins = 80) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30", linewidth = 0.6) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_fill_viridis_c(option = "plasma", name = "Cell count") +
  scale_colour_viridis_c(option = "inferno", name = "% bg removed") +
  labs(
    title    = "Raw UMI vs CellBender UMI per barcode",
    subtitle = "Dashed line = no change; colour = % background removed",
    x = "Raw total UMI (log10)",
    y = "CellBender UMI (log10)"
  ) +
  theme_qc()

# -------------------------------------------------------
# Assemble + save combined PDF
# -------------------------------------------------------
message("\nSaving plots...")

pdf(
  file.path(out_dir, paste0(sample_id, "_cellbender_qc.pdf")),
  width = 12, height = 16
)

print(
  (p_knee | plot_spacer()) /
    (p_umi_density | p_umi_ecdf) /
    (p_bg_hist | p_scatter) +
    plot_annotation(
      title    = paste0("CellBender QC — ", sample_id),
      subtitle = paste0(
        "Raw: ", comma(ncol(raw_mat)), " droplets  |  ",
        "CB output: ", comma(ncol(cb_mat)), " droplets  |  ",
        "Median background removed: ", round(med_bg, 1), "%"
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, colour = "grey40")
      )
    )
)

dev.off()

# Also save individual PNGs for quick inspection
ggsave(file.path(out_dir, paste0(sample_id, "_01_knee_plot.png")),
       p_knee,       width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, paste0(sample_id, "_02_umi_density.png")),
       p_umi_density, width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, paste0(sample_id, "_03_umi_ecdf.png")),
       p_umi_ecdf,    width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, paste0(sample_id, "_04_pct_bg_removed.png")),
       p_bg_hist,     width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, paste0(sample_id, "_05_raw_vs_cb_scatter.png")),
       p_scatter,     width = 7, height = 5, dpi = 150)

message("\nDone. Plots saved to: ", out_dir)
