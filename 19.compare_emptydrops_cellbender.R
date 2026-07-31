#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(dplyr)
  library(tibble)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
run <- "260528_VH01624_464_222K7VKNX"
sample <- "LK2-GEX"

out_dir <- file.path(proj, "results", "cellbender_comparison", run, sample)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Load barcode lists
# ------------------------------------------------------------------

emptydrops_barcodes <- readLines(
  file.path(proj, "results", "emptydrops", run, sample, "emptydrops_barcodes.txt")
)

cellbender_32k15k_barcodes <- readLines(
  file.path(proj, "results", "cellbender", run, sample, paste0(sample, "_cellbender_cell_barcodes.csv"))
)

cellbender_25k15k_barcodes <- readLines(
  file.path(proj, "results", "cellbender", run, paste0(sample, "_25000"), paste0(sample, "_cellbender_cell_barcodes.csv"))
)

cellbender_32k30k_barcodes <- readLines(
  file.path(proj, "results", "cellbender", run, paste0(sample, "_e30000"), paste0(sample, "_cellbender_cell_barcodes.csv"))
)

cat("Barcode counts:\n")
cat("  EmptyDrops:           ", length(emptydrops_barcodes), "\n")
cat("  CellBender (32k exp): ", length(cellbender_32k_barcodes), "\n")
cat("  CellBender (25k exp): ", length(cellbender_25k_barcodes), "\n")

# ------------------------------------------------------------------
# 2. Venn diagram: all three methods
# ------------------------------------------------------------------

barcode_sets <- list(
  EmptyDrops           = emptydrops_barcodes,
  `CellBender (32k15k)`   = cellbender_32k15k_barcodes,
  `CellBender (25k15k)`   = cellbender_25k15k_barcodes,
  `CellBender (32k30k)`   = cellbender_32k30k_barcodes
)

p_venn4 <- ggVennDiagram(barcode_sets, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#4292C6") +
  ggtitle(paste0("Cell calling comparison - ", sample)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave(
  file.path(out_dir, "venn_emptydrops_cellbender_4way.pdf"),
  p_venn4, width = 8, height = 7
)

# ------------------------------------------------------------------
# 3. Pairwise Venn diagrams
# ------------------------------------------------------------------

pairs <- list(
  list(
    sets = list(EmptyDrops = emptydrops_barcodes, `CellBender (32k)` = cellbender_32k_barcodes),
    name = "venn_emptydrops_vs_cellbender_32k.pdf"
  ),
  list(
    sets = list(EmptyDrops = emptydrops_barcodes, `CellBender (25k)` = cellbender_25k_barcodes),
    name = "venn_emptydrops_vs_cellbender_25k.pdf"
  ),
  list(
    sets = list(`CellBender (32k)` = cellbender_32k_barcodes, `CellBender (25k)` = cellbender_25k_barcodes),
    name = "venn_cellbender_32k_vs_25k.pdf"
  )
)

for (pair in pairs) {
  p <- ggVennDiagram(pair$sets, label = "count", label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "#4292C6") +
    ggtitle(paste0(sample, ": ", paste(names(pair$sets), collapse = " vs "))) +
    theme(plot.title = element_text(hjust = 0.5, size = 13))

  ggsave(file.path(out_dir, pair$name), p, width = 7, height = 6)
}

# ------------------------------------------------------------------
# 4. Summary table
# ------------------------------------------------------------------

all_barcodes <- unique(c(emptydrops_barcodes, cellbender_32k_barcodes, cellbender_25k_barcodes))

summary_df <- tibble(
  barcode          = all_barcodes,
  EmptyDrops       = all_barcodes %in% emptydrops_barcodes,
  CellBender_32k   = all_barcodes %in% cellbender_32k_barcodes,
  CellBender_25k   = all_barcodes %in% cellbender_25k_barcodes
) %>%
  mutate(
    n_methods = EmptyDrops + CellBender_32k + CellBender_25k,
    category = case_when(
      EmptyDrops & CellBender_32k & CellBender_25k ~ "All three",
      EmptyDrops & CellBender_32k                  ~ "EmptyDrops + CB_32k",
      EmptyDrops & CellBender_25k                  ~ "EmptyDrops + CB_25k",
      CellBender_32k & CellBender_25k              ~ "CB_32k + CB_25k",
      EmptyDrops                                    ~ "EmptyDrops only",
      CellBender_32k                                ~ "CellBender_32k only",
      CellBender_25k                                ~ "CellBender_25k only"
    )
  )

cat("\nCategory breakdown:\n")
print(summary_df %>% dplyr::count(category, name = "n_barcodes") %>% arrange(desc(n_barcodes)))

write.csv(summary_df, file.path(out_dir, "barcode_comparison_full.csv"), row.names = FALSE)

category_summary <- summary_df %>%
  dplyr::count(category, name = "n_barcodes") %>%
  arrange(desc(n_barcodes))
write.csv(category_summary, file.path(out_dir, "barcode_comparison_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 5. Bar plot of category breakdown
# ------------------------------------------------------------------

p_bar <- ggplot(category_summary, aes(x = reorder(category, n_barcodes), y = n_barcodes, fill = category)) +
  geom_col() +
  geom_text(aes(label = n_barcodes), hjust = -0.1, size = 3.5) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("Cell calling agreement - ", sample),
    x = NULL, y = "Number of barcodes"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(out_dir, "barplot_category_breakdown.pdf"), p_bar, width = 9, height = 5)

cat("\nOutputs saved to:", out_dir, "\n")
cat("Done.\n")
