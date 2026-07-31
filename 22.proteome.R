#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# ============================================================
# Inputs
# ============================================================

knorr_file     <- "../proteome/2119KK_S1_C_ProteinQuant_Rerun.txt"
bordeleau_file <- "../proteome/evidence.txt"

out_dir <- "../results/cryptic_surfaceome"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target   <- "SNRNP200"
receptor <- "FCGR2A"

# ============================================================
# Helpers
# ============================================================

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

clean_base_names <- function(df) {
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  df
}

standardise_gene_names <- function(x) {
  x <- as.character(x)
  x <- gsub("^;+|;+$", "", x)
  x <- sub(";.*$", "", x)
  trimws(x)
}

# ============================================================
# Read Knorr: wide protein-level table
# ============================================================

read_knorr <- function(file) {
  
  message("Reading Knorr file: ", file)
  
  knorr <- fread(file, header = TRUE, sep = "\t", data.table = FALSE)
  knorr <- clean_base_names(knorr)
  
  quantity_cols <- grep("PG\\.Quantity$", colnames(knorr), value = TRUE)
  
  if (length(quantity_cols) == 0) {
    stop("No PG.Quantity columns found in Knorr file.")
  }
  
  required <- c("PG.Genes", "PG.ProteinDescriptions")
  missing <- setdiff(required, colnames(knorr))
  
  if (length(missing) > 0) {
    stop("Knorr file missing columns: ", paste(missing, collapse = ", "))
  }
  
  knorr %>%
    as_tibble() %>%
    select(
      gene_names = PG.Genes,
      protein_names = PG.ProteinDescriptions,
      all_of(quantity_cols)
    ) %>%
    pivot_longer(
      cols = all_of(quantity_cols),
      names_to = "sample_id",
      values_to = "intensity"
    ) %>%
    mutate(
      dataset = "Knorr",
      gene_names = standardise_gene_names(gene_names),
      intensity = as.numeric(intensity)
    ) %>%
    filter(!is.na(gene_names), gene_names != "") %>%
    select(dataset, sample_id, gene_names, protein_names, intensity)
}

# ============================================================
# Read Bordeleau: long peptide-level evidence table
# ============================================================

read_bordeleau <- function(file) {
  
  message("Reading Bordeleau file: ", file)
  
  bordeleau <- fread(file, header = TRUE, sep = "\t", data.table = FALSE)
  bordeleau <- clean_base_names(bordeleau)
  
  required <- c("Gene.names", "Protein.names", "Raw.file", "Intensity")
  missing <- setdiff(required, colnames(bordeleau))
  
  if (length(missing) > 0) {
    stop("Bordeleau file missing columns: ", paste(missing, collapse = ", "))
  }
  
  bordeleau %>%
    as_tibble() %>%
    transmute(
      dataset = "Bordeleau",
      sample_id = Raw.file,
      gene_names = standardise_gene_names(Gene.names),
      protein_names = Protein.names,
      intensity = as.numeric(Intensity)
    ) %>%
    filter(!is.na(gene_names), gene_names != "")
}

# ============================================================
# Protein-level summary
# ============================================================

make_protein_summary <- function(df) {
  
  df %>%
    group_by(dataset, gene_names) %>%
    summarise(
      protein_names = paste(unique(na.omit(protein_names)), collapse = "; "),
      mean_intensity = safe_mean(intensity),
      count = n(),
      n_detected = sum(!is.na(intensity)),
      n_samples_detected = n_distinct(sample_id[!is.na(intensity)]),
      .groups = "drop"
    ) %>%
    filter(n_detected > 0, !is.na(mean_intensity)) %>%
    mutate(
      detection_score = log10(mean_intensity + 1) * log2(n_samples_detected + 1)
    ) %>%
    arrange(dataset, desc(detection_score))
}

# ============================================================
# Sample x protein matrix
# ============================================================

make_sample_matrix <- function(df) {
  
  df %>%
    filter(!is.na(gene_names), gene_names != "", !is.na(intensity)) %>%
    group_by(dataset, sample_id, gene_names) %>%
    summarise(
      intensity = safe_mean(intensity),
      .groups = "drop"
    ) %>%
    mutate(sample_id = paste(dataset, sample_id, sep = "__")) %>%
    select(sample_id, gene_names, intensity) %>%
    pivot_wider(
      names_from = gene_names,
      values_from = intensity
    )
}

# ============================================================
# Correlations
# ============================================================

make_correlations <- function(protein_mat,
                              target = "SNRNP200",
                              receptor = "FCGR2A") {
  
  if (is.null(protein_mat) || nrow(protein_mat) < 3) {
    return(tibble(
      gene_names = character(),
      cor_with_snrnp200 = numeric(),
      cor_with_fcgr2a = numeric()
    ))
  }
  
  mat <- protein_mat %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
  
  mat_log <- log2(mat + 1)
  
  if (!target %in% colnames(mat_log)) {
    warning(target, " not found in matrix.")
    return(tibble(
      gene_names = colnames(mat_log),
      cor_with_snrnp200 = NA_real_,
      cor_with_fcgr2a = NA_real_
    ))
  }
  
  snrnp200_cor <- tibble(
    gene_names = colnames(mat_log),
    cor_with_snrnp200 = apply(
      mat_log,
      2,
      function(x) suppressWarnings(
        cor(
          x,
          mat_log[, target],
          use = "pairwise.complete.obs",
          method = "spearman"
        )
      )
    )
  )
  
  if (receptor %in% colnames(mat_log)) {
    fcgr2a_cor <- tibble(
      gene_names = colnames(mat_log),
      cor_with_fcgr2a = apply(
        mat_log,
        2,
        function(x) suppressWarnings(
          cor(
            x,
            mat_log[, receptor],
            use = "pairwise.complete.obs",
            method = "spearman"
          )
        )
      )
    )
  } else {
    warning(receptor, " not found in matrix.")
    
    fcgr2a_cor <- tibble(
      gene_names = colnames(mat_log),
      cor_with_fcgr2a = NA_real_
    )
  }
  
  snrnp200_cor %>%
    left_join(fcgr2a_cor, by = "gene_names") %>%
    filter(gene_names != target)
}

# ============================================================
# 1. Load and standardise both datasets
# ============================================================

knorr_std <- read_knorr(knorr_file)
bordeleau_std <- read_bordeleau(bordeleau_file)

combined <- bind_rows(knorr_std, bordeleau_std)

write.csv(
  combined,
  file.path(out_dir, "combined_standardised_long.csv"),
  row.names = FALSE
)

cat("\nCombined standardised data:\n")
print(table(combined$dataset))

# ============================================================
# 2. Protein summary
# ============================================================

protein_summary <- make_protein_summary(combined)

write.csv(
  protein_summary,
  file.path(out_dir, "protein_summary_by_dataset.csv"),
  row.names = FALSE
)

# ============================================================
# 3. SNRNP200 summary
# ============================================================

snrnp200_summary <- protein_summary %>%
  filter(gene_names == target)

write.csv(
  snrnp200_summary,
  file.path(out_dir, "SNRNP200_summary_by_dataset.csv"),
  row.names = FALSE
)

cat("\nSNRNP200 summary:\n")
print(snrnp200_summary)

# ============================================================
# 4. Within-dataset normalisation
# ============================================================

protein_summary_norm <- protein_summary %>%
  group_by(dataset) %>%
  mutate(
    intensity_rank = rank(-mean_intensity, ties.method = "average"),
    intensity_percentile = percent_rank(mean_intensity),
    log10_intensity = log10(mean_intensity + 1),
    z_log10_intensity = as.numeric(scale(log10_intensity)),
    detection_percent = n_samples_detected / max(n_samples_detected, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv(
  protein_summary_norm,
  file.path(out_dir, "protein_summary_by_dataset_normalised.csv"),
  row.names = FALSE
)

# ============================================================
# 5. Similar to SNRNP200 within each dataset
# ============================================================

snrnp200_ref <- protein_summary_norm %>%
  filter(gene_names == target) %>%
  select(
    dataset,
    snrnp200_z = z_log10_intensity,
    snrnp200_percentile = intensity_percentile
  )

similar_to_snrnp200 <- protein_summary_norm %>%
  left_join(snrnp200_ref, by = "dataset") %>%
  filter(!is.na(snrnp200_z)) %>%
  mutate(
    delta_z = abs(z_log10_intensity - snrnp200_z),
    delta_percentile = abs(intensity_percentile - snrnp200_percentile)
  ) %>%
  filter(delta_z < 0.5) %>%
  arrange(dataset, delta_z)

write.csv(
  similar_to_snrnp200,
  file.path(out_dir, "proteins_with_intensity_similar_to_SNRNP200.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Overlap across datasets
# ============================================================

overlap_summary <- protein_summary_norm %>%
  group_by(gene_names) %>%
  summarise(
    present_in = paste(sort(unique(dataset)), collapse = ";"),
    n_datasets = n_distinct(dataset),
    mean_detection_score = mean(detection_score, na.rm = TRUE),
    mean_detection_percent = mean(detection_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_datasets), desc(mean_detection_score))

write.csv(
  overlap_summary,
  file.path(out_dir, "protein_overlap_knorr_bordeleau.csv"),
  row.names = FALSE
)

# ============================================================
# 7. Correlations within each dataset
# ============================================================

correlations <- combined %>%
  split(.$dataset) %>%
  imap_dfr(function(df, ds) {
    mat <- make_sample_matrix(df)
    
    make_correlations(
      mat,
      target = target,
      receptor = receptor
    ) %>%
      mutate(dataset = ds)
  })

write.csv(
  correlations,
  file.path(out_dir, "correlations_with_SNRNP200_and_FCGR2A.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Rank cryptic surface candidates
# ============================================================

interesting_intracellular_keywords <- c(
  "splice", "snrnp", "ribonucleoprotein", "rna.binding",
  "helicase", "nucleolar", "nuclear", "chromatin",
  "hnrnp", "splicing", "mrna", "spliceosome"
)

likely_contaminant_keywords <- c(
  "ribosomal", "ribosome",
  "histone",
  "heat shock", "chaperone",
  "elongation factor", "translation elongation",
  "actin", "tubulin", "keratin",
  "mitochondrial", "mitochondrion",
  "atp synthase",
  "glyceraldehyde", "gapdh",
  "enolase", "aldolase",
  "peroxiredoxin",
  "proteasome"
)

surface_keywords <- c(
  "cd", "integrin", "receptor", "transmembrane",
  "membrane", "surface", "secreted", "extracellular",
  "adhesion", "glycoprotein"
)

ranked_candidates <- protein_summary_norm %>%
  left_join(correlations, by = c("dataset", "gene_names")) %>%
  left_join(overlap_summary, by = "gene_names") %>%
  mutate(
    protein_text = str_to_lower(paste(gene_names, protein_names)),
    
    interesting_intracellular_like = str_detect(
      protein_text,
      paste(interesting_intracellular_keywords, collapse = "|")
    ),
    
    likely_contaminant_like =
      str_detect(
        protein_text,
        paste(likely_contaminant_keywords, collapse = "|")
      ) |
      str_detect(
        gene_names,
        "^(RPL|RPS|MRPL|MRPS|HIST|HSP|HSPA|HSPD|EEF|TUB|ACT|KRT)"
      ),
    
    canonical_surface_like = str_detect(
      protein_text,
      paste(surface_keywords, collapse = "|")
    ),
    
    weird_surface_score =
      z_log10_intensity +
      log2(n_samples_detected + 1) +
      5 * coalesce(cor_with_snrnp200, 0) +
      3 * coalesce(cor_with_fcgr2a, 0) +
      2 * coalesce(n_datasets, 1) +
      if_else(interesting_intracellular_like, 4, 0) -
      if_else(canonical_surface_like, 3, 0) -
      if_else(likely_contaminant_like, 8, 0)
  ) %>%
  arrange(desc(weird_surface_score))

write.csv(
  ranked_candidates,
  file.path(out_dir, "ranked_cryptic_surface_candidates_by_dataset.csv"),
  row.names = FALSE
)

# ============================================================
# 9. Collapse to gene level
# ============================================================

ranked_candidates_gene_level <- ranked_candidates %>%
  group_by(gene_names) %>%
  summarise(
    protein_names = paste(unique(na.omit(protein_names)), collapse = "; "),
    present_in = paste(sort(unique(dataset)), collapse = ";"),
    n_datasets = n_distinct(dataset),
    best_weird_surface_score = max(weird_surface_score, na.rm = TRUE),
    mean_weird_surface_score = mean(weird_surface_score, na.rm = TRUE),
    mean_detection_score = mean(detection_score, na.rm = TRUE),
    max_z_log10_intensity = max(z_log10_intensity, na.rm = TRUE),
    mean_cor_with_snrnp200 = mean(cor_with_snrnp200, na.rm = TRUE),
    mean_cor_with_fcgr2a = mean(cor_with_fcgr2a, na.rm = TRUE),
    interesting_intracellular_like = any(interesting_intracellular_like, na.rm = TRUE),
    likely_contaminant_like = any(likely_contaminant_like, na.rm = TRUE),
    canonical_surface_like = any(canonical_surface_like, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(best_weird_surface_score))

write.csv(
  ranked_candidates_gene_level,
  file.path(out_dir, "ranked_cryptic_surface_candidates_gene_level.csv"),
  row.names = FALSE
)

# ============================================================
# 10. Clean cryptic shortlist
# ============================================================

cryptic_surface_shortlist <- ranked_candidates_gene_level %>%
  filter(
    interesting_intracellular_like,
    !likely_contaminant_like,
    !canonical_surface_like
  ) %>%
  arrange(desc(best_weird_surface_score))

write.csv(
  cryptic_surface_shortlist,
  file.path(out_dir, "cryptic_surface_shortlist_no_ribosomal_histone_hsp.csv"),
  row.names = FALSE
)

# ============================================================
# 11. Contaminant-like high scoring proteins
# ============================================================

contaminant_like_hits <- ranked_candidates_gene_level %>%
  filter(likely_contaminant_like) %>%
  arrange(desc(best_weird_surface_score))

write.csv(
  contaminant_like_hits,
  file.path(out_dir, "likely_contaminant_like_hits.csv"),
  row.names = FALSE
)

# ============================================================
# 12. Print outputs
# ============================================================

cat("\nTop cryptic surface candidates, excluding ribosomal/histone/HSP/etc:\n")

print(
  cryptic_surface_shortlist %>%
    select(
      gene_names,
      present_in,
      n_datasets,
      best_weird_surface_score,
      mean_detection_score,
      mean_cor_with_snrnp200,
      mean_cor_with_fcgr2a,
      interesting_intracellular_like,
      likely_contaminant_like,
      canonical_surface_like
    ) %>%
    head(50)
)

cat("\nTop likely contaminant-like hits:\n")

print(
  contaminant_like_hits %>%
    select(
      gene_names,
      present_in,
      n_datasets,
      best_weird_surface_score,
      mean_detection_score,
      mean_cor_with_snrnp200,
      mean_cor_with_fcgr2a,
      likely_contaminant_like
    ) %>%
    head(30)
)

cat("\nDone. Outputs written to:\n")
cat(out_dir, "\n")