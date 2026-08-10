#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Cell type annotation benchmark - step 04 of 05
#
# Scores the annotation itself, rather than producing more of it. Reads only
# the per-cell metadata tables written by steps 01-03, so it runs in seconds
# and does not need the Seurat objects or the reference.
#
# Six questions:
#   1. how many cells does BoneMarrowMap's own mapping-error QC reject, per run
#      and per sample
#   2. what are those cells like - depth, mitochondrial and ribosomal content
#   3. how confident are the labels that survive
#   4. does predict_CellTypes() change the label relative to initial_CellType,
#      or does it only gate on QC
#   5. fine (predicted_CellType) versus broad (predicted_CellType_Broad)
#      vocabulary, and what it costs to key downstream grouping on the wrong one
#   6. what MAD_threshold 2.5 versus 4 costs
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

options(scipen = 999)

res     <- "../results"
out_dir <- file.path(res, "benchmarking", "celltype_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PILOT <- "260423_VH01624_453_222HWMYNX"
LK2   <- "260528_VH01624_464_222K7VKNX"

meta_files <- c(
  "LK1 pilot" = file.path(res, "seurat_annotated", PILOT, "cell_metadata_CITE_DSB_BoneMarrowMap.csv"),
  "LK2"       = file.path(res, "seurat_annotated", LK2,   "cell_metadata_demux_ADT_BoneMarrowMap.csv")
)

for (f in meta_files) if (!file.exists(f)) stop("Missing metadata table: ", f)

read_meta <- function(path, label) {
  d <- suppressMessages(read_csv(path, show_col_types = FALSE))
  names(d)[1] <- "barcode"
  # the pilot calls the demultiplexed sample sampleID, LK2 calls it sample_name
  if (!"sample" %in% names(d)) {
    d$sample <- if ("sampleID" %in% names(d)) d$sampleID else d$sample_name
  }
  d$run <- label
  d
}

meta <- imap_dfr(meta_files, ~ read_meta(.x, .y))

# ----------------------------
# 1. Mapping-error QC
# ----------------------------

qc_by_run <- meta %>%
  count(run, mapping_error_QC) %>%
  group_by(run) %>%
  mutate(total = sum(n), percent = round(100 * n / total, 2)) %>%
  ungroup() %>%
  arrange(run, mapping_error_QC)

write_csv(qc_by_run, file.path(out_dir, "01_mapping_QC_by_run.csv"))

qc_by_sample <- meta %>%
  count(run, sample, mapping_error_QC) %>%
  group_by(run, sample) %>%
  mutate(total = sum(n), percent = round(100 * n / total, 2)) %>%
  ungroup() %>%
  filter(mapping_error_QC == "Fail") %>%
  arrange(desc(percent))

write_csv(qc_by_sample, file.path(out_dir, "02_mapping_QC_fail_by_sample.csv"))

p_fail <- ggplot(qc_by_sample, aes(x = reorder(sample, percent), y = percent, fill = run)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("LK1 pilot" = "#1B9E77", "LK2" = "#E31A1C")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Cells rejected by BoneMarrowMap mapping-error QC",
    subtitle = "MAD_threshold = 2.5",
    x = NULL, y = "Percent of sample rejected", fill = "Run"
  )

ggsave(file.path(out_dir, "01_mapping_QC_fail_by_sample.pdf"), p_fail, width = 9, height = 5)

# ----------------------------
# 2. What the rejected cells look like
# ----------------------------

qc_cols <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
    "scDblFinder.score", "mapping_error_score"),
  names(meta)
)

fail_profile <- meta %>%
  group_by(run, mapping_error_QC) %>%
  summarise(n = n(), across(all_of(qc_cols), ~ round(median(.x, na.rm = TRUE), 3)),
            .groups = "drop")

write_csv(fail_profile, file.path(out_dir, "03_rejected_cell_profile.csv"))

p_profile <- meta %>%
  dplyr::select(run, mapping_error_QC, nFeature_RNA, percent.mt, percent.ribo) %>%
  pivot_longer(c(nFeature_RNA, percent.mt, percent.ribo),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = mapping_error_QC, y = value, fill = mapping_error_QC)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.2) +
  facet_grid(metric ~ run, scales = "free_y") +
  scale_fill_manual(values = c(Pass = "grey70", Fail = "firebrick"), guide = "none") +
  theme_bw(base_size = 12) +
  labs(
    title = "QC covariates of cells the reference rejects",
    x = NULL, y = NULL
  )

ggsave(file.path(out_dir, "02_rejected_cell_profile.pdf"), p_profile, width = 9, height = 7)

# ----------------------------
# 3. Label confidence
# ----------------------------

conf_summary <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  group_by(run) %>%
  summarise(
    n = n(),
    prob_p05 = round(quantile(predicted_CellType_prob, 0.05, na.rm = TRUE), 3),
    prob_q25 = round(quantile(predicted_CellType_prob, 0.25, na.rm = TRUE), 3),
    prob_median = round(median(predicted_CellType_prob, na.rm = TRUE), 3),
    prob_q75 = round(quantile(predicted_CellType_prob, 0.75, na.rm = TRUE), 3),
    percent_above_0.5 = round(100 * mean(predicted_CellType_prob > 0.5, na.rm = TRUE), 1),
    percent_above_0.8 = round(100 * mean(predicted_CellType_prob > 0.8, na.rm = TRUE), 1),
    .groups = "drop"
  )

write_csv(conf_summary, file.path(out_dir, "04_label_confidence.csv"))

# predict_CellTypes() calls knnPredict_Seurat() with k = 30 and applies NO
# probability cutoff - the only thing it gates on is mapping_error_QC. So
# predicted_CellType_prob is reported, never acted on, and every value is a
# multiple of 1/30. The reference lines below are ours, for reading the plot;
# they are not thresholds the tool applies.
K_NN <- 30

conf_lines <- tibble(
  x = c(1 / K_NN, 0.5, 1),
  label = c("1/30 - one neighbour", "0.5 - majority of the 30", "1.0 - unanimous")
)

p_conf <- meta %>%
  filter(mapping_error_QC == "Pass", !is.na(predicted_CellType_prob)) %>%
  ggplot(aes(x = predicted_CellType_prob, fill = run)) +
  geom_histogram(binwidth = 1 / K_NN, alpha = 0.75, position = "identity") +
  geom_vline(data = conf_lines, aes(xintercept = x),
             linetype = "dashed", linewidth = 0.4, colour = "grey20") +
  geom_text(data = conf_lines, aes(x = x, y = Inf, label = label),
            inherit.aes = FALSE, angle = 90, vjust = -0.4, hjust = 1.05,
            size = 3, colour = "grey20") +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1.02)) +
  scale_fill_manual(values = c("LK1 pilot" = "#1B9E77", "LK2" = "#E31A1C")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Fraction of the k = 30 nearest reference neighbours voting for the assigned label",
    subtitle = paste("BoneMarrowMap applies no cutoff to this value - the dashed lines are",
                     "reading aids, not thresholds.\nOne bar = one neighbour changing its vote."),
    x = "predicted_CellType_prob", y = "Cells", fill = "Run"
  )

ggsave(file.path(out_dir, "03_label_confidence.pdf"), p_conf, width = 10, height = 5.5)

# Confidence per broad cell type - which populations the reference is unsure of
conf_by_type <- meta %>%
  filter(mapping_error_QC == "Pass", !is.na(predicted_CellType_Broad)) %>%
  group_by(run, predicted_CellType_Broad) %>%
  summarise(n = n(), median_prob = round(median(predicted_CellType_prob, na.rm = TRUE), 3),
            .groups = "drop") %>%
  filter(n >= 20) %>%
  arrange(median_prob)

write_csv(conf_by_type, file.path(out_dir, "05_label_confidence_by_celltype.csv"))

p_conf_type <- conf_by_type %>%
  ggplot(aes(x = reorder(predicted_CellType_Broad, median_prob), y = median_prob, fill = run)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("LK1 pilot" = "#1B9E77", "LK2" = "#E31A1C")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Median label confidence by broad cell type",
    subtitle = "populations with at least 20 cells",
    x = NULL, y = "Median predicted_CellType_prob", fill = "Run"
  )

ggsave(file.path(out_dir, "04_label_confidence_by_celltype.pdf"), p_conf_type,
       width = 9, height = 7)

# ----------------------------
# 4. Does predict_CellTypes() re-label, or only gate
# ----------------------------

relabel <- meta %>%
  group_by(run) %>%
  summarise(
    n = n(),
    n_initial_NA = sum(is.na(initial_CellType)),
    n_final_NA = sum(is.na(predicted_CellType)),
    n_QC_fail = sum(mapping_error_QC == "Fail"),
    final_NA_equals_QC_fail = identical(
      which(is.na(predicted_CellType)), which(mapping_error_QC == "Fail")
    ),
    percent_label_unchanged = round(
      100 * mean(initial_CellType == predicted_CellType, na.rm = TRUE), 2
    ),
    .groups = "drop"
  )

write_csv(relabel, file.path(out_dir, "06_initial_vs_final_label.csv"))

# ----------------------------
# 5. Fine versus broad vocabulary
# ----------------------------

broad_vocab <- c(
  "HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
  "Cycling Progenitor", "EoBasoMast Precursor", "Megakaryocyte Precursor",
  "Monocyte", "Pro-Monocyte", "cDC", "pDC",
  "Naive T", "CD4 Memory T", "CD8 Memory T", "NK", "Early Lymphoid",
  "B", "Pre-B", "Pro-B", "Plasma Cell",
  "Early Erythroid", "Late Erythroid"
)

vocab <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  group_by(run) %>%
  summarise(
    n_fine_labels = n_distinct(predicted_CellType),
    n_broad_labels = n_distinct(predicted_CellType_Broad),
    fine_labels_in_broad_vocab = sum(unique(predicted_CellType) %in% broad_vocab),
    percent_cells_Other_if_keyed_on_fine =
      round(100 * mean(!predicted_CellType %in% broad_vocab), 1),
    percent_cells_Other_if_keyed_on_broad =
      round(100 * mean(!predicted_CellType_Broad %in% broad_vocab), 1),
    .groups = "drop"
  )

write_csv(vocab, file.path(out_dir, "07_fine_vs_broad_vocabulary.csv"))

# Same thing per sample, which is what the lineage bar chart actually plots
lineage_of <- function(x) {
  case_when(
    x %in% c("HSC MPP", "LMPP", "MEP", "GMP", "Early GMP", "Late GMP",
             "Cycling Progenitor", "EoBasoMast Precursor",
             "Megakaryocyte Precursor") ~ "Stem / progenitor",
    x %in% c("Monocyte", "Pro-Monocyte", "cDC", "pDC") ~ "Myeloid / DC",
    x %in% c("Naive T", "CD4 Memory T", "CD8 Memory T", "NK",
             "Early Lymphoid", "B", "Pre-B", "Pro-B", "Plasma Cell") ~ "Lymphoid",
    x %in% c("Early Erythroid", "Late Erythroid") ~ "Erythroid",
    TRUE ~ "Other"
  )
}

lineage_cmp <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  mutate(
    lineage_from_fine  = lineage_of(predicted_CellType),
    lineage_from_broad = lineage_of(predicted_CellType_Broad)
  ) %>%
  dplyr::select(run, sample, lineage_from_fine, lineage_from_broad) %>%
  pivot_longer(starts_with("lineage_from_"), names_to = "keyed_on", values_to = "lineage") %>%
  mutate(keyed_on = recode(keyed_on,
                           lineage_from_fine = "predicted_CellType (fine)",
                           lineage_from_broad = "predicted_CellType_Broad")) %>%
  count(run, sample, keyed_on, lineage) %>%
  group_by(run, sample, keyed_on) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup()

write_csv(lineage_cmp, file.path(out_dir, "08_lineage_keyed_on_fine_vs_broad.csv"))

lineage_cols <- c(
  "Stem / progenitor" = "#1B9E77",
  "Myeloid / DC" = "#E31A1C",
  "Lymphoid" = "#2171B5",
  "Erythroid" = "#C51B8A",
  "Other" = "grey70"
)

p_lineage_cmp <- ggplot(lineage_cmp, aes(sample, percent, fill = lineage)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ keyed_on, ncol = 2) +
  scale_fill_manual(values = lineage_cols) +
  theme_bw(base_size = 11) +
  labs(
    title = "The same lineage grouping, keyed on the fine label and on the broad label",
    subtitle = "the fine vocabulary does not overlap the lineage lists, so almost everything falls to Other",
    x = NULL, y = "Cellular composition (%)", fill = "Lineage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "05_lineage_keyed_on_fine_vs_broad.pdf"), p_lineage_cmp,
       width = 11, height = 5.5)

# Which fine labels roll up into which broad label - the mapping itself
vocab_map <- meta %>%
  filter(mapping_error_QC == "Pass") %>%
  count(predicted_CellType_Broad, predicted_CellType, name = "n_cells") %>%
  arrange(predicted_CellType_Broad, desc(n_cells))

write_csv(vocab_map, file.path(out_dir, "09_fine_to_broad_label_map.csv"))

# ----------------------------
# 6. What the MAD threshold costs
# ----------------------------
# calculate_MappingError() stores the per-cell mapping_error_score and then
# calls Pass/Fail on
#
#   score < median(score) + MAD_threshold * mad(score)
#
# with median and mad taken over the whole object. The score does not depend on
# the threshold, so the sweep can be done directly from the stored score - no
# need to re-project. The 2.5 column is checked against the stored
# mapping_error_QC as a correctness test before anything else is reported.

thresholds <- c(2.5, 3, 4)

sweep <- meta %>%
  group_by(run) %>%
  group_modify(~ {
    s   <- .x$mapping_error_score
    cut <- stats::median(s) + thresholds * stats::mad(s)
    map_dfr(seq_along(thresholds), function(i) {
      .x %>%
        mutate(fail = s >= cut[i]) %>%
        count(sample, fail) %>%
        group_by(sample) %>%
        mutate(total = sum(n)) %>%
        ungroup() %>%
        filter(fail) %>%
        transmute(MAD_threshold = thresholds[i], sample,
                  n_rejected = n, n_cells = total,
                  percent_rejected = round(100 * n / total, 2))
    })
  }) %>%
  ungroup()

# Correctness check: the recomputed 2.5 call must reproduce the stored one
check <- meta %>%
  group_by(run) %>%
  mutate(recomputed = ifelse(
    mapping_error_score < stats::median(mapping_error_score) +
      2.5 * stats::mad(mapping_error_score), "Pass", "Fail")) %>%
  summarise(agreement = round(100 * mean(recomputed == mapping_error_QC), 4),
            .groups = "drop")

print(as.data.frame(check))
if (any(check$agreement < 100)) {
  warning("Recomputed MAD 2.5 call does not reproduce the stored mapping_error_QC - ",
          "the sweep below cannot be trusted.")
}

write_csv(sweep, file.path(out_dir, "11_MAD_threshold_sweep_by_sample.csv"))

sweep_run <- meta %>%
  group_by(run) %>%
  group_modify(~ {
    s   <- .x$mapping_error_score
    cut <- stats::median(s) + thresholds * stats::mad(s)
    tibble(
      MAD_threshold = thresholds,
      n_cells = nrow(.x),
      n_rejected = map_int(cut, ~ sum(s >= .x)),
      percent_rejected = round(100 * map_dbl(cut, ~ mean(s >= .x)), 2)
    )
  }) %>%
  ungroup()

write_csv(sweep_run, file.path(out_dir, "12_MAD_threshold_sweep_by_run.csv"))

p_sweep <- ggplot(sweep,
                  aes(x = factor(MAD_threshold), y = percent_rejected,
                      group = sample, colour = sample)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ run, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Per-sample rejection rate against MAD_threshold",
    subtitle = "pipeline annotate.R runs at 4; the standalone scripts ran at 2.5",
    x = "MAD_threshold", y = "Percent of sample rejected", colour = NULL
  )

ggsave(file.path(out_dir, "06_MAD_threshold_sweep.pdf"), p_sweep, width = 10, height = 5)

# The original two-object comparison, kept because step 03 genuinely ran on a
# different object at 4 rather than being a recomputation of the same cells.

nf_qc <- file.path(res, "seurat_annotated", LK2,
                   "LK2-GEX_results_nf_unannotated_check",
                   "unknown_vs_annotated_qc_medians.csv")

threshold_rows <- qc_by_run %>%
  filter(mapping_error_QC == "Fail") %>%
  transmute(
    object = paste0(run, " (step ", ifelse(run == "LK2", "02", "01"), ")"),
    MAD_threshold = 2.5,
    n_cells = total,
    n_rejected = n,
    percent_rejected = percent
  )

if (file.exists(nf_qc)) {
  nfd <- suppressMessages(read_csv(nf_qc, show_col_types = FALSE))
  threshold_rows <- bind_rows(
    threshold_rows,
    tibble(
      object = "LK2 results_nf checkpoint (step 03)",
      MAD_threshold = 4,
      n_cells = sum(nfd$n),
      n_rejected = nfd$n[nfd$is_unknown],
      percent_rejected = round(100 * nfd$n[nfd$is_unknown] / sum(nfd$n), 2)
    )
  )
} else {
  message("Step 03 output not found - MAD threshold table will cover 2.5 only.")
}

write_csv(threshold_rows, file.path(out_dir, "10_MAD_threshold_cost.csv"))

# ----------------------------
# Report
# ----------------------------

cat("\n=== Mapping QC by run ===\n");            print(as.data.frame(qc_by_run))
cat("\n=== Fail rate by sample ===\n");          print(as.data.frame(qc_by_sample))
cat("\n=== Rejected cell profile ===\n");        print(as.data.frame(fail_profile))
cat("\n=== Label confidence ===\n");             print(as.data.frame(conf_summary))
cat("\n=== initial vs final label ===\n");       print(as.data.frame(relabel))
cat("\n=== Fine vs broad vocabulary ===\n");     print(as.data.frame(vocab))
cat("\n=== MAD threshold, two objects ===\n");   print(as.data.frame(threshold_rows))
cat("\n=== MAD threshold sweep, by run ===\n");  print(as.data.frame(sweep_run))
cat("\n=== MAD threshold sweep, by sample ===\n")
print(as.data.frame(sweep %>% arrange(run, sample, MAD_threshold)))

cat("\nDone. Output directory:", out_dir, "\n")
