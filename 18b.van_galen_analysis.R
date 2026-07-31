#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
xlsx_file <- file.path(
  proj,
  "annotation",
  "van_galen_1-s2.0-S0092867419300947-mmc4.xlsx"
)

out_dir <- file.path(proj, "annotation", "van_galen_signatures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------

cor_threshold <- 0.25
top_n <- 50

# ----------------------------
# Read Table S4A
# ----------------------------

s4a <- read_excel(xlsx_file, sheet = "Table S4A", skip = 1)

colnames(s4a) <- c(
  "gene",
  "expr_normal_HSC_Prog", "expr_normal_GMP", "expr_normal_Myeloid",
  "expr_malignant_HSC_Prog", "expr_malignant_GMP", "expr_malignant_Myeloid",
  "cor_HSC_Prog", "cor_GMP", "cor_Myeloid"
)

s4a <- s4a %>%
  mutate(
    gene = as.character(gene),
    across(starts_with("expr_") | starts_with("cor_"), as.numeric)
  ) %>%
  filter(!is.na(gene), gene != "")

write.csv(
  s4a,
  file.path(out_dir, "van_galen_S4A_cleaned.csv"),
  row.names = FALSE
)

cat("Table S4A:", nrow(s4a), "genes\n")
cat("Correlation ranges:\n")
cat("  HSC/Prog:", range(s4a$cor_HSC_Prog, na.rm = TRUE), "\n")
cat("  GMP:     ", range(s4a$cor_GMP, na.rm = TRUE), "\n")
cat("  Myeloid: ", range(s4a$cor_Myeloid, na.rm = TRUE), "\n")

# ----------------------------
# Build correlation-only signatures
# ----------------------------

build_cor_signature <- function(df, cor_col, label, threshold, top_n) {
  sig_all <- df %>%
    filter(!is.na(.data[[cor_col]])) %>%
    filter(.data[[cor_col]] > threshold) %>%
    arrange(desc(.data[[cor_col]])) %>%
    dplyr::select(gene, correlation = all_of(cor_col))
  
  sig_top <- sig_all %>%
    slice_head(n = top_n)
  
  cat(sprintf(
    "%s: %d genes above cor > %.2f; using top %d for AUC\n",
    label, nrow(sig_all), threshold, nrow(sig_top)
  ))
  
  list(top = sig_top, all = sig_all)
}

hsc  <- build_cor_signature(s4a, "cor_HSC_Prog", "HSC/Prog", cor_threshold, top_n)
gmp  <- build_cor_signature(s4a, "cor_GMP", "GMP", cor_threshold, top_n)
myel <- build_cor_signature(s4a, "cor_Myeloid", "Myeloid", cor_threshold, top_n)

# ----------------------------
# AUCell / AUC gene sets
# ----------------------------

auc_signatures_top <- list(
  VanGalen_HSC_Prog_cor = hsc$top$gene,
  VanGalen_GMP_cor      = gmp$top$gene,
  VanGalen_Myeloid_cor  = myel$top$gene
)

auc_signatures_all <- list(
  VanGalen_HSC_Prog_cor = hsc$all$gene,
  VanGalen_GMP_cor      = gmp$all$gene,
  VanGalen_Myeloid_cor  = myel$all$gene
)

saveRDS(
  auc_signatures_top,
  file.path(out_dir, "van_galen_AUC_signatures_top_correlated.rds")
)

saveRDS(
  auc_signatures_all,
  file.path(out_dir, "van_galen_AUC_signatures_all_correlated.rds")
)

# Save CSVs with correlation values
write.csv(hsc$top,  file.path(out_dir, "van_galen_HSC_Prog_top_correlated.csv"), row.names = FALSE)
write.csv(gmp$top,  file.path(out_dir, "van_galen_GMP_top_correlated.csv"), row.names = FALSE)
write.csv(myel$top, file.path(out_dir, "van_galen_Myeloid_top_correlated.csv"), row.names = FALSE)

write.csv(hsc$all,  file.path(out_dir, "van_galen_HSC_Prog_all_correlated.csv"), row.names = FALSE)
write.csv(gmp$all,  file.path(out_dir, "van_galen_GMP_all_correlated.csv"), row.names = FALSE)
write.csv(myel$all, file.path(out_dir, "van_galen_Myeloid_all_correlated.csv"), row.names = FALSE)

# Save plain text gene lists
for (nm in names(auc_signatures_top)) {
  writeLines(
    auc_signatures_top[[nm]],
    file.path(out_dir, paste0(nm, "_top", top_n, ".txt"))
  )
  
  writeLines(
    auc_signatures_all[[nm]],
    file.path(out_dir, paste0(nm, "_all_cor_gt_", cor_threshold, ".txt"))
  )
}

cat("\nSaved AUC correlation-only signatures to:\n")
cat(out_dir, "\n")
cat("Done.\n")

