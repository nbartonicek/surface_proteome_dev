#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
})

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
xlsx_file <- file.path(proj, "annotation", "van_galen_1-s2.0-S0092867419300947-mmc4.xlsx")
out_dir <- file.path(proj, "annotation", "van_galen_signatures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# 1. Read Table S4A (gene-level expression + correlation scores)
# ------------------------------------------------------------------

s4a <- read_excel(xlsx_file, sheet = "Table S4A", skip = 1)
colnames(s4a) <- c(
  "gene",
  "expr_normal_HSC_Prog", "expr_normal_GMP", "expr_normal_Myeloid",
  "expr_malignant_HSC_Prog", "expr_malignant_GMP", "expr_malignant_Myeloid",
  "cor_HSC_Prog", "cor_GMP", "cor_Myeloid"
)

s4a <- s4a %>%
  mutate(across(starts_with("expr_") | starts_with("cor_"), as.numeric))

write.csv(s4a, file.path(out_dir, "van_galen_S4A_cleaned.csv"), row.names = FALSE)

cat("Table S4A:", nrow(s4a), "genes\n")
cat("Correlation columns:\n")
cat("  HSC/Prog range:", range(s4a$cor_HSC_Prog, na.rm = TRUE), "\n")
cat("  GMP range:     ", range(s4a$cor_GMP, na.rm = TRUE), "\n")
cat("  Myeloid range: ", range(s4a$cor_Myeloid, na.rm = TRUE), "\n")

# ------------------------------------------------------------------
# 2. Build signatures from correlation scores
#    Positive correlation = gene associated with that lineage in
#    malignant cells. Use top genes by correlation as signature.
# ------------------------------------------------------------------

cor_threshold <- 0.25
top_n <- 50

build_signature <- function(df, cor_col, label, threshold = cor_threshold, n = top_n) {
  sig <- df %>%
    filter(.data[[cor_col]] > threshold) %>%
    arrange(desc(.data[[cor_col]])) %>%
    slice_head(n = n) %>%
    dplyr::select(gene, correlation = !!sym(cor_col))

  cat(sprintf("  %s: %d genes (cor > %.2f, top %d)\n", label, nrow(sig), threshold, n))
  sig
}

cat("\nBuilding signatures (correlation-based, top", top_n, "genes with cor >", cor_threshold, "):\n")

sig_hsc  <- build_signature(s4a, "cor_HSC_Prog", "HSC/Prog")
sig_gmp  <- build_signature(s4a, "cor_GMP", "GMP")
sig_myel <- build_signature(s4a, "cor_Myeloid", "Myeloid")

# ------------------------------------------------------------------
# 3. Also build signatures with all positively correlated genes
# ------------------------------------------------------------------

cat("\nAll positively correlated genes (cor >", cor_threshold, "):\n")

sig_hsc_all  <- s4a %>% filter(cor_HSC_Prog > cor_threshold) %>%
  arrange(desc(cor_HSC_Prog)) %>% dplyr::select(gene, correlation = cor_HSC_Prog)
sig_gmp_all  <- s4a %>% filter(cor_GMP > cor_threshold) %>%
  arrange(desc(cor_GMP)) %>% dplyr::select(gene, correlation = cor_GMP)
sig_myel_all <- s4a %>% filter(cor_Myeloid > cor_threshold) %>%
  arrange(desc(cor_Myeloid)) %>% dplyr::select(gene, correlation = cor_Myeloid)

cat("  HSC/Prog:", nrow(sig_hsc_all), "genes\n")
cat("  GMP:     ", nrow(sig_gmp_all), "genes\n")
cat("  Myeloid: ", nrow(sig_myel_all), "genes\n")

# ------------------------------------------------------------------
# 4. Save signatures as gene lists and as a named list object
# ------------------------------------------------------------------

signatures_top <- list(
  HSC_Prog = sig_hsc$gene,
  GMP      = sig_gmp$gene,
  Myeloid  = sig_myel$gene
)

signatures_all <- list(
  HSC_Prog = sig_hsc_all$gene,
  GMP      = sig_gmp_all$gene,
  Myeloid  = sig_myel_all$gene
)

saveRDS(signatures_top, file.path(out_dir, "van_galen_signatures_top50.rds"))
saveRDS(signatures_all, file.path(out_dir, "van_galen_signatures_all_positive.rds"))

for (nm in names(signatures_top)) {
  writeLines(signatures_top[[nm]], file.path(out_dir, paste0("van_galen_", nm, "_top50.txt")))
  writeLines(signatures_all[[nm]], file.path(out_dir, paste0("van_galen_", nm, "_all_positive.txt")))
}

write.csv(sig_hsc, file.path(out_dir, "van_galen_HSC_Prog_top50.csv"), row.names = FALSE)
write.csv(sig_gmp, file.path(out_dir, "van_galen_GMP_top50.csv"), row.names = FALSE)
write.csv(sig_myel, file.path(out_dir, "van_galen_Myeloid_top50.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 5. Also save Table S4B (monocyte expression across AML patients)
# ------------------------------------------------------------------

s4b <- read_excel(xlsx_file, sheet = "Table S4B", skip = 1)
colnames(s4b)[1] <- "gene"

write.csv(s4b, file.path(out_dir, "van_galen_S4B_monocyte_expression.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 6. NPM1 class I / class II signatures (Uckelmann et al. 2025)
# ------------------------------------------------------------------

npm1_file <- file.path(proj, "annotation", "41467_2025_66546_MOESM3_ESM_signature.csv")
npm1 <- read.csv(npm1_file)
colnames(npm1) <- c("gene", "class", "t_stat", "pvalue", "qvalue", "fold_change")

npm1 <- npm1 %>%
  mutate(across(c(t_stat, pvalue, qvalue, fold_change), as.numeric))

sig_npm1_I  <- npm1 %>% filter(class == "NPM1 class I") %>% arrange(qvalue)
sig_npm1_II <- npm1 %>% filter(class == "NPM1 class II") %>% arrange(qvalue)

cat("\nNPM1 signatures:\n")
cat("  Class I: ", nrow(sig_npm1_I), "genes\n")
cat("  Class II:", nrow(sig_npm1_II), "genes\n")

signatures_top[["NPM1_classI"]]  <- sig_npm1_I$gene
signatures_top[["NPM1_classII"]] <- sig_npm1_II$gene

signatures_all[["NPM1_classI"]]  <- sig_npm1_I$gene
signatures_all[["NPM1_classII"]] <- sig_npm1_II$gene

saveRDS(signatures_top, file.path(out_dir, "van_galen_signatures_top50.rds"))
saveRDS(signatures_all, file.path(out_dir, "van_galen_signatures_all_positive.rds"))

writeLines(sig_npm1_I$gene, file.path(out_dir, "NPM1_classI.txt"))
writeLines(sig_npm1_II$gene, file.path(out_dir, "NPM1_classII.txt"))
write.csv(sig_npm1_I, file.path(out_dir, "NPM1_classI.csv"), row.names = FALSE)
write.csv(sig_npm1_II, file.path(out_dir, "NPM1_classII.csv"), row.names = FALSE)

cat("\nAll signatures saved to:", out_dir, "\n")
cat("Done.\n")
