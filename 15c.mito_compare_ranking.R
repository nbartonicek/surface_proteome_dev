#!/usr/bin/env Rscript

# Aggregates every <donor>_ari_summary.csv written by
# 15b.mito_mitoclone2_and_compare.R / 15b2.mito_compare_only.R across all
# donors in a run, and visualises cross-method clone-calling concordance:
# a method x method heatmap (joint mean-across-donors panel + one panel
# per donor), a per-method concordance ranking (average ARI vs every other
# method), and a category breakdown.
#
# Excludes mgatk_merged/mitoclone2_merged (SR+LR combined - dropped per
# request, not yet fully understood/trusted) and copykat (a CNV-subclone
# caller, not a clone-calling method comparable to the rest). Only
# numbat, mgatk_sr, mgatk_lr, mitoclone2_sr, mitoclone2_lr remain.
#
# The category split is the actual diagnostic tool: "same algorithm,
# cross-platform" (e.g. mgatk_sr vs mgatk_lr) tests whether ONE method
# reproduces across SR/LR of the SAME cells, while "cross-algorithm,
# mtDNA-only" (mgatk vs mitoClone2, any platform) tests whether TWO
# different algorithms find the same lineage structure from the SAME
# reads. If even the same-algorithm cross-platform numbers are weak, that
# points to a signal/depth problem (e.g. no mtDNA enrichment) rather than
# the algorithms just disagreeing with each other.
#
# Needs only base tidyverse (dplyr/tidyr/readr/purrr/ggplot2) - no mgatk,
# Signac, or mitoClone2 - so it runs in any env, including outside the
# cluster (e.g. against an SMB-mounted copy of results/).
#
# Usage:
#   Rscript 15c.mito_compare_ranking.R [run_id] [project_dir]
#   Both optional - default to the run this pipeline has been analyzing
#   (260528_VH01624_464_222K7VKNX) and ".." (project_dir relative to this
#   script's own location, same convention as compare_isoform_annotations.R -
#   run from the scripts/ directory), not a hardcoded absolute cluster path,
#   so `Rscript 15c.mito_compare_ranking.R` with no args just works regardless
#   of the actual absolute mount path on whatever machine it's run from.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else "260528_VH01624_464_222K7VKNX"
PROJECT_DIR <- if (length(args) >= 2) args[2] else ".."
PER_PATIENT_DIR <- file.path(PROJECT_DIR, "results/mitochondrial_clones/per_patient", RUN_ID)
OUT_DIR <- file.path(PROJECT_DIR, "results/mitochondrial_clones/comparison_ranking", RUN_ID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(PER_PATIENT_DIR, pattern = "_ari_summary\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No *_ari_summary.csv files found under ", PER_PATIENT_DIR)
message("Found ", length(files), " donor ARI summary file(s): ", paste(basename(files), collapse = ", "))

ari_all_raw <- map_dfr(files, read_csv, show_col_types = FALSE)

# Dropped per request: "merged" (SR+LR combined) categories, since their
# construction isn't fully trusted/understood yet, and copykat, since it's
# a CNV-subclone caller rather than a clone-calling method comparable to
# the others here. Everything below (heatmaps, ranking, categories) only
# ever sees the 5 remaining methods: numbat, mgatk_sr, mgatk_lr,
# mitoclone2_sr, mitoclone2_lr.
DROPPED_METHODS <- c("mgatk_merged", "mitoclone2_merged", "copykat")
ari_all <- ari_all_raw %>% filter(!method_1 %in% DROPPED_METHODS, !method_2 %in% DROPPED_METHODS)

classify_pair <- function(m1, m2) {
  base <- function(x) sub("_(sr|lr|merged)$", "", x)
  plat <- function(x) ifelse(grepl("_(sr|lr|merged)$", x), sub(".*_", "", x), NA)
  b1 <- base(m1); b2 <- base(m2); p1 <- plat(m1); p2 <- plat(m2)
  dplyr::case_when(
    b1 == "numbat" & b2 == "numbat" ~ "CNV vs CNV",
    b1 == "numbat" | b2 == "numbat" ~ "CNV vs mtDNA",
    b1 == b2 & p1 != p2 ~ "Same algorithm, cross-platform",
    b1 != b2 ~ "Cross-algorithm (mgatk vs mitoClone2)",
    TRUE ~ "other"
  )
}
ari_all <- ari_all %>% mutate(category = classify_pair(method_1, method_2))
write_csv(ari_all, file.path(OUT_DIR, "ari_all_donors_categorized.csv"))

method_order <- c("numbat", "mgatk_sr", "mgatk_lr", "mitoclone2_sr", "mitoclone2_lr")
methods_present <- intersect(method_order, unique(c(ari_all$method_1, ari_all$method_2)))
donors_present <- sort(unique(ari_all$donor))

# ---- Final heatmap: joint (mean across donors) + one panel per donor ----
joint_data <- ari_all %>%
  group_by(method_1, method_2) %>%
  summarise(ari = mean(ari), .groups = "drop") %>%
  mutate(panel = "Joint (mean across donors)")
per_patient_data <- ari_all %>% transmute(method_1, method_2, ari, panel = donor)

heat_data <- bind_rows(joint_data, per_patient_data)
heat_data <- bind_rows(
  heat_data %>% rename(row = method_1, col = method_2),
  heat_data %>% rename(row = method_2, col = method_1)
) %>%
  mutate(row = factor(row, levels = methods_present),
         col = factor(col, levels = methods_present),
         panel = factor(panel, levels = c("Joint (mean across donors)", donors_present)),
         label_color = if_else(abs(ari) > 0.4, "white", "grey20"))

max_abs_ari <- max(abs(heat_data$ari), na.rm = TRUE)

p_heat <- ggplot(heat_data, aes(x = col, y = row, fill = ari)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", ari), color = label_color), size = 2.8, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(low = "#e34948", mid = "#f0efec", high = "#2a78d6", midpoint = 0,
                        limits = c(-max_abs_ari, max_abs_ari), name = "ARI") +
  coord_fixed() +
  facet_wrap(~panel, ncol = 2) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10)) +
  labs(title = paste0("Cross-method clone-calling concordance - ", RUN_ID),
       subtitle = "copykat and merged (SR+LR) methods excluded; joint panel = mean ARI across donors",
       x = NULL, y = NULL)
ggsave(file.path(OUT_DIR, "ari_heatmap.pdf"), p_heat, width = 10, height = 9)
ggsave(file.path(OUT_DIR, "ari_heatmap.png"), p_heat, width = 10, height = 9, dpi = 150)

# ---- Ranking: per-method average concordance with all other methods ----
long <- bind_rows(
  ari_all %>% transmute(method = method_1, other = method_2, ari, donor),
  ari_all %>% transmute(method = method_2, other = method_1, ari, donor)
)
ranking <- long %>%
  group_by(method) %>%
  summarise(mean_ari_vs_others = mean(ari), n = n(), .groups = "drop") %>%
  arrange(desc(mean_ari_vs_others)) %>%
  mutate(method = factor(method, levels = rev(method)),
         label_hjust = if_else(mean_ari_vs_others >= 0, -0.15, 1.15))
write_csv(ranking, file.path(OUT_DIR, "method_ranking.csv"))

p_rank <- ggplot(ranking, aes(x = mean_ari_vs_others, y = method)) +
  geom_col(fill = "#2a78d6", width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", mean_ari_vs_others), hjust = label_hjust), size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "Method ranking: average ARI vs every other method",
       subtitle = "Higher = this method's clone calls agree more with the rest of the panel, on average",
       x = "Mean ARI vs all other methods", y = NULL)
ggsave(file.path(OUT_DIR, "method_ranking.pdf"), p_rank, width = 8, height = 5)

# ---- Category summary: where does concordance actually come from? ----
category_summary <- ari_all %>%
  group_by(category) %>%
  summarise(mean_ari = mean(ari), median_ari = median(ari), n_pairs = n(), .groups = "drop") %>%
  arrange(desc(mean_ari)) %>%
  mutate(category = factor(category, levels = rev(category)),
         label = sprintf("%.3f (n=%d)", mean_ari, n_pairs))
write_csv(category_summary, file.path(OUT_DIR, "category_summary.csv"))

p_cat <- ggplot(category_summary, aes(x = mean_ari, y = category)) +
  geom_col(fill = "#1baf7a", width = 0.6) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.5) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "Where does concordance come from?",
       subtitle = "Mean ARI grouped by what kind of comparison the method pair represents",
       x = "Mean ARI", y = NULL) +
  xlim(min(category_summary$mean_ari) - 0.02, max(category_summary$mean_ari) + 0.12)
ggsave(file.path(OUT_DIR, "category_summary.pdf"), p_cat, width = 8, height = 4)

message("\nWrote outputs to: ", OUT_DIR)
message("  - ari_heatmap.pdf / ari_heatmap.png (joint + per-donor panels)")
message("  - method_ranking.pdf / method_ranking.csv")
message("  - category_summary.pdf / category_summary.csv")
message("  - ari_all_donors_categorized.csv (raw, for further slicing; copykat/merged already excluded)")
