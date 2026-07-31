#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

root <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX/09c_numbat_run"
out_dir <- file.path(root, "numbat_interpretable_cnv_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- tibble(
  sample = list.dirs(root, recursive = FALSE, full.names = FALSE)
) |>
  filter(str_detect(sample, "^LK2_")) |>
  mutate(
    dir = file.path(root, sample, "numbat_final"),
    group = if_else(str_detect(sample, regex("normal", ignore_case = TRUE)), "Normal", "AML")
  )

safe_fread <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  tryCatch(fread(path), error = function(e) NULL)
}

guess_col <- function(dt, patterns) {
  hits <- names(dt)[str_detect(names(dt), regex(paste(patterns, collapse = "|"), ignore_case = TRUE))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

summarise_joint <- function(sample, dir, group) {
  x <- safe_fread(file.path(dir, "joint_post_1.tsv"))
  if (is.null(x) || nrow(x) == 0) return(NULL)
  
  message("\n--- ", sample, " joint_post columns ---")
  print(names(x))
  
  cell_col <- guess_col(x, c("^cell$", "barcode", "cell_id"))
  if (is.na(cell_col)) cell_col <- names(x)[1]
  
  seg_col <- guess_col(x, c("seg", "region", "bin", "locus"))
  cnv_col <- guess_col(x, c("cnv", "copy", "state", "z_cnv", "Z_cnv"))
  prob_col <- guess_col(x, c("p_cnv", "prob", "post", "p_"))
  
  num_cols <- names(x)[map_lgl(x, is.numeric)]
  
  z_cols <- num_cols[str_detect(num_cols, regex("z_cnv|Z_cnv", ignore_case = TRUE))]
  p_cols <- num_cols[str_detect(num_cols, regex("p_cnv|prob|post|p_", ignore_case = TRUE))]
  
  z_use <- if (length(z_cols) > 0) z_cols[1] else NA_character_
  p_use <- if (length(p_cols) > 0) p_cols[1] else NA_character_
  
  x <- as_tibble(x)
  
  out <- x |>
    group_by(cell = .data[[cell_col]]) |>
    summarise(
      n_segments = n(),
      
      cnv_z_mean = if (!is.na(z_use)) mean(.data[[z_use]], na.rm = TRUE) else NA_real_,
      cnv_z_abs_mean = if (!is.na(z_use)) mean(abs(.data[[z_use]]), na.rm = TRUE) else NA_real_,
      cnv_z_min = if (!is.na(z_use)) min(.data[[z_use]], na.rm = TRUE) else NA_real_,
      cnv_z_max = if (!is.na(z_use)) max(.data[[z_use]], na.rm = TRUE) else NA_real_,
      cnv_z_burden_25 = if (!is.na(z_use)) mean(abs(.data[[z_use]]) > 25, na.rm = TRUE) else NA_real_,
      cnv_z_burden_50 = if (!is.na(z_use)) mean(abs(.data[[z_use]]) > 50, na.rm = TRUE) else NA_real_,
      cnv_z_burden_100 = if (!is.na(z_use)) mean(abs(.data[[z_use]]) > 100, na.rm = TRUE) else NA_real_,
      
      cnv_p_mean = if (!is.na(p_use)) mean(.data[[p_use]], na.rm = TRUE) else NA_real_,
      cnv_p_max = if (!is.na(p_use)) max(.data[[p_use]], na.rm = TRUE) else NA_real_,
      cnv_p_burden_50 = if (!is.na(p_use)) mean(.data[[p_use]] > 0.50, na.rm = TRUE) else NA_real_,
      cnv_p_burden_90 = if (!is.na(p_use)) mean(.data[[p_use]] > 0.90, na.rm = TRUE) else NA_real_,
      
      .groups = "drop"
    ) |>
    mutate(sample = sample, group = group, source = "joint_post")
  
  out
}

summarise_exp <- function(sample, dir, group) {
  x <- safe_fread(file.path(dir, "exp_post_1.tsv"))
  if (is.null(x) || nrow(x) == 0) return(NULL)
  
  message("\n--- ", sample, " exp_post columns ---")
  print(names(x))
  
  cell_col <- guess_col(x, c("^cell$", "barcode", "cell_id"))
  if (is.na(cell_col)) cell_col <- names(x)[1]
  
  num_cols <- names(x)[map_lgl(x, is.numeric)]
  z_cols <- num_cols[str_detect(num_cols, regex("z_cnv|Z_cnv", ignore_case = TRUE))]
  p_cols <- num_cols[str_detect(num_cols, regex("p_cnv|prob|post|p_", ignore_case = TRUE))]
  
  z_use <- if (length(z_cols) > 0) z_cols[1] else NA_character_
  p_use <- if (length(p_cols) > 0) p_cols[1] else NA_character_
  
  as_tibble(x) |>
    group_by(cell = .data[[cell_col]]) |>
    summarise(
      n_exp_records = n(),
      exp_z_mean = if (!is.na(z_use)) mean(.data[[z_use]], na.rm = TRUE) else NA_real_,
      exp_z_abs_mean = if (!is.na(z_use)) mean(abs(.data[[z_use]]), na.rm = TRUE) else NA_real_,
      exp_z_burden_25 = if (!is.na(z_use)) mean(abs(.data[[z_use]]) > 25, na.rm = TRUE) else NA_real_,
      exp_z_burden_50 = if (!is.na(z_use)) mean(abs(.data[[z_use]]) > 50, na.rm = TRUE) else NA_real_,
      exp_p_mean = if (!is.na(p_use)) mean(.data[[p_use]], na.rm = TRUE) else NA_real_,
      exp_p_max = if (!is.na(p_use)) max(.data[[p_use]], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    mutate(sample = sample, group = group, source = "exp_post")
}

joint_metrics <- pmap_dfr(samples, \(sample, dir, group) summarise_joint(sample, dir, group))
exp_metrics   <- pmap_dfr(samples, \(sample, dir, group) summarise_exp(sample, dir, group))

cell_metrics <- full_join(
  joint_metrics |> select(-source),
  exp_metrics |> select(-source, -group),
  by = c("cell", "sample")
) |>
  left_join(samples |> select(sample, group), by = "sample", suffix = c("", ".y")) |>
  mutate(group = coalesce(group, group.y)) |>
  select(-any_of("group.y"))

write_csv(cell_metrics, file.path(out_dir, "per_cell_interpretable_cnv_metrics.csv"))

metrics_to_plot <- c(
  "cnv_z_abs_mean",
  "cnv_z_burden_25",
  "cnv_z_burden_50",
  "cnv_z_burden_100",
  "cnv_p_mean",
  "cnv_p_max",
  "cnv_p_burden_50",
  "cnv_p_burden_90",
  "exp_z_abs_mean",
  "exp_z_burden_25",
  "exp_z_burden_50",
  "exp_p_mean",
  "exp_p_max"
)

plot_df <- cell_metrics |>
  select(cell, sample, group, any_of(metrics_to_plot)) |>
  pivot_longer(any_of(metrics_to_plot), names_to = "metric", values_to = "value") |>
  filter(is.finite(value))

p1 <- ggplot(plot_df, aes(x = sample, y = value, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.65) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, alpha = 0.85) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Interpretable Numbat CNV evidence per cell",
    subtitle = "Use CNV burden / altered-segment burden, not clone posterior, to judge whether normal cells are truly tumour-like",
    x = NULL,
    y = "Value"
  )

ggsave(file.path(out_dir, "01_interpretable_per_cell_cnv_metrics.pdf"),
       p1, width = 15, height = 10)

normal_cutoffs <- plot_df |>
  filter(group == "Normal") |>
  group_by(metric) |>
  summarise(
    normal_median = median(value, na.rm = TRUE),
    normal_p95 = quantile(value, 0.95, na.rm = TRUE),
    normal_p99 = quantile(value, 0.99, na.rm = TRUE),
    normal_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(normal_cutoffs, file.path(out_dir, "normal_empirical_cutoffs.csv"))

exceedance <- plot_df |>
  left_join(normal_cutoffs, by = "metric") |>
  group_by(sample, group, metric) |>
  summarise(
    n_cells = n(),
    median = median(value, na.rm = TRUE),
    frac_above_normal_p95 = mean(value > normal_p95, na.rm = TRUE),
    frac_above_normal_p99 = mean(value > normal_p99, na.rm = TRUE),
    frac_above_normal_max = mean(value > normal_max, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(exceedance, file.path(out_dir, "sample_exceedance_above_normal.csv"))

p2 <- exceedance |>
  filter(str_detect(metric, "burden|abs_mean|p_max")) |>
  pivot_longer(starts_with("frac_above"), names_to = "cutoff", values_to = "fraction") |>
  ggplot(aes(x = sample, y = fraction, fill = group)) +
  geom_col() +
  facet_grid(cutoff ~ metric, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "How many cells exceed the normal empirical CNV-evidence distribution?",
    x = NULL,
    y = "Fraction of cells"
  )

ggsave(file.path(out_dir, "02_fraction_above_normal_cnv_cutoffs.pdf"),
       p2, width = 16, height = 9)

sample_summary <- samples |>
  rowwise() |>
  mutate(
    allele_post_size = file.info(file.path(dir, "allele_post_1.tsv"))$size,
    joint_post_size = file.info(file.path(dir, "joint_post_1.tsv"))$size,
    exp_post_size = file.info(file.path(dir, "exp_post_1.tsv"))$size,
    clone_post_size = file.info(file.path(dir, "clone_post_1.tsv"))$size,
    n_consensus_segments = {
      x <- safe_fread(file.path(dir, "segs_consensus_1.tsv"))
      if (is.null(x)) NA_integer_ else nrow(x)
    },
    n_loh_segments = {
      x <- safe_fread(file.path(dir, "segs_loh.tsv"))
      if (is.null(x)) NA_integer_ else nrow(x)
    }
  ) |>
  ungroup()

write_csv(sample_summary, file.path(out_dir, "sample_level_file_segment_summary.csv"))

p3 <- sample_summary |>
  select(sample, group, allele_post_size, joint_post_size, exp_post_size,
         clone_post_size, n_consensus_segments, n_loh_segments) |>
  pivot_longer(-c(sample, group), names_to = "metric", values_to = "value") |>
  ggplot(aes(x = sample, y = value, fill = group)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Sample-level Numbat evidence",
    subtitle = "Normal may have consensus segments, but should have much weaker allele/joint evidence",
    x = NULL,
    y = "Value"
  )

ggsave(file.path(out_dir, "03_sample_level_evidence.pdf"),
       p3, width = 13, height = 7)

# Simple call table using normal p99 on the best interpretable metrics
call_metrics <- c(
  "cnv_z_abs_mean",
  "cnv_z_burden_50",
  "cnv_z_burden_100",
  "cnv_p_burden_90",
  "exp_z_abs_mean",
  "exp_z_burden_50"
)

call_table <- exceedance |>
  filter(metric %in% call_metrics) |>
  select(sample, group, metric, median, frac_above_normal_p99, frac_above_normal_max) |>
  arrange(metric, group, sample)

write_csv(call_table, file.path(out_dir, "tumour_like_fraction_by_metric.csv"))

message("\nDone.")
message("Output directory: ", out_dir)
message("\nMost useful files:")
message("  01_interpretable_per_cell_cnv_metrics.pdf")
message("  02_fraction_above_normal_cnv_cutoffs.pdf")
message("  tumour_like_fraction_by_metric.csv")
message("  normal_empirical_cutoffs.csv")