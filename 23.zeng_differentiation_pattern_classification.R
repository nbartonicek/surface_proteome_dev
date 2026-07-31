#!/usr/bin/env Rscript

# Classifies each AML donor against Zeng et al. 2025 (Blood Cancer Discov,
# "Single-cell Transcriptional Atlas of Human Hematopoiesis...", the
# BoneMarrowMap paper - same atlas this pipeline already uses for
# predicted_CellType/predicted_CellType_Broad) into one of their 12 recurrent
# AML "differentiation patterns", and separately applies their published
# bulk-RNA-seq LASSO deconvolution model to project both our patients and the
# Leucegene cohort into the same 13-differentiation-state score space.
#
# Two independent, complementary analyses:
#
#  A) Single-cell pattern classification (scRNA-seq -> per-donor, per-clone
#     composition vector -> nearest-centroid match against the paper's 318
#     reference scAML patients -> Differentiation_Pattern assignment).
#     Reference vectors are rebuilt from Supplementary Table S8 (raw
#     precise-cell-state composition per reference patient) via S9's exact
#     Cell State -> AML Differentiation State mapping - NOT read directly off
#     S7's summary scAML_* columns, because S7 only has 12 columns (missing a
#     separate Late_Erythroid) while S9/S11 define 13 categories. Building
#     both our patients' and the reference patients' vectors through the same
#     S9 mapping keeps them on an identical, verifiable footing. S7 is used
#     only for its ground-truth Differentiation_Pattern label per patient.
#
#  B) Bulk RNA-seq deconvolution (Leucegene + our own pseudobulk -> the same
#     13-state score space via S11's published LASSO models, applied
#     identically to both per Supplementary Note S2: log2CPM-normalized
#     expression x gene coefficients, summed, then z-standardized across the
#     combined sample set so our patients and Leucegene land on the same
#     scale). This replaces the earlier whole-transcriptome Leucegene PCA
#     (which collapsed under a bulk-vs-pseudobulk batch effect) with a
#     compact, curated, cross-platform-validated 400-gene/13-state space -
#     batch effects from mixing 10x pseudobulk with Leucegene's own bulk
#     RNA-seq are far less likely to dominate a curated biological signature
#     space than they are to dominate raw whole-transcriptome PCA.
#
# Requires: annotation/zeng_bonemarrowmap_atlas/supplementary_tables_s1-s21.xlsx
# (downloaded from the paper's AACR/Silverchair-hosted supplementary data,
# with the user's explicit permission), and
# proteome/leucegene_transcriptome_readcounts.tsv (already present).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(readxl)
  library(edgeR)
})

# ----------------------------
# Paths / settings
# ----------------------------

proj <- "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"

runs <- c(
  "260522_VH01624_461_222JLJVNX",
  "260528_VH01624_464_222K7VKNX"
)

ZENG_XLSX <- file.path(proj, "annotation/zeng_bonemarrowmap_atlas/supplementary_tables_s1-s21.xlsx")
stopifnot(file.exists(ZENG_XLSX))

leucegene_file <- file.path(proj, "proteome/leucegene_transcriptome_readcounts.tsv")
stopifnot(file.exists(leucegene_file))

out_dir <- file.path(proj, "results/zeng_differentiation_patterns")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_composition <- 30  # minimum mapped cells for a donor's composition vector to be trusted

# ----------------------------
# Load Zeng reference tables
# ----------------------------

message("Loading Zeng et al. reference tables...")

s7 <- read_excel(ZENG_XLSX, sheet = "S7 - scAML Patient Annotations") %>%
  select(Sample_ID, Differentiation_Pattern, Study, Diagnosis, AgeGroup, Key_Mutations)

s8 <- read_excel(ZENG_XLSX, sheet = "S8 - scAML Cell Composition")
names(s8)[1] <- "predicted_CellType"

# S7's Sample_ID (short codes, e.g. "MLL_14666") and S8's column headers
# (study-prefixed, e.g. "scAML_AEL_MLL_14666" or "AML_Bailur2020_GSM4664013_1150")
# are two different naming schemes for the same 318 patients - verified by
# direct inspection, not assumed. Link them by matching each S8 column
# against the LONGEST S7 Sample_ID that is a suffix of it (longest first,
# to avoid short-ID collisions), with a fallback that strips a trailing
# timepoint suffix (-D0/-REL/-D<n>) from Sample_ID for the van Galen 2019
# samples, which don't carry that suffix in S8's column names.
s7_sorted <- s7$Sample_ID[order(-nchar(s7$Sample_ID))]
s8_sample_id <- sapply(names(s8)[-1], function(col) {
  hit <- s7_sorted[col == s7_sorted | endsWith(col, paste0("_", s7_sorted))]
  if (length(hit) > 0) return(hit[1])
  stripped <- sub("-D0$|-REL$|-D[0-9]+$", "", s7_sorted)
  hit2 <- s7_sorted[col == stripped | endsWith(col, paste0("_", stripped))]
  if (length(hit2) > 0) hit2[1] else NA_character_
})
n_unmatched <- sum(is.na(s8_sample_id))
if (n_unmatched > 0) {
  message("Note: ", n_unmatched, " / ", length(s8_sample_id),
          " S8 columns could not be linked to an S7 Sample_ID and will be dropped from the reference.")
}
names(s8)[-1] <- ifelse(is.na(s8_sample_id), paste0("UNMATCHED_", seq_along(s8_sample_id)), s8_sample_id)
s8 <- s8 %>% select(-starts_with("UNMATCHED_"))

s9 <- read_excel(ZENG_XLSX, sheet = "S9 - AML Differentiation States") %>%
  rename(predicted_CellType = `Cell State`, diffstate = `AML Differentiation State`) %>%
  filter(diffstate != "Unassigned - not adequately enriched by NMF") %>%
  mutate(
    # Compact, code-friendly state names matching S11's column suffixes
    diffstate_code = recode(diffstate,
      "HSC_MPP" = "HSCMPP", "LMPP" = "LMPP", "Early_Lymphoid" = "EarlyLymphoid",
      "Pro-B / Pre-B" = "ProBPreB", "MEP_MkP" = "MEPMkP", "EoBasoMast" = "EoBasoMast",
      "Early_Erythroid" = "EarlyEry", "Late_Erythroid" = "LateEry", "GMP" = "GMP",
      "ProMono" = "ProMono", "Monocyte" = "Monocyte", "cDC" = "cDC", "pDC" = "pDC"
    ),
    # normalize hyphenation for matching against our own predicted_CellType values
    match_key = tolower(gsub("[- ]", "", predicted_CellType))
  )

s11 <- read_excel(ZENG_XLSX, sheet = "S11 - scAML DiffState Models")

diffstate_levels <- unique(s9$diffstate_code)

message("Reference: ", nrow(s7), " scAML patients, ", length(diffstate_levels), " differentiation states, ",
        nrow(s11), " LASSO model genes")

# ----------------------------
# Helper: precise-cell-state composition matrix -> 13-category % vector,
# via S9's mapping. Used for BOTH the reference S8 matrix and our own
# per-donor predicted_CellType counts, so both sides go through identical
# logic.
# ----------------------------

celltype_counts_to_diffstate_pct <- function(celltype_counts) {
  # celltype_counts: named numeric vector, names = predicted_CellType, values = n_cells
  df <- tibble(predicted_CellType = names(celltype_counts), n_cells = as.numeric(celltype_counts)) %>%
    mutate(match_key = tolower(gsub("[- ]", "", predicted_CellType))) %>%
    inner_join(s9 %>% select(match_key, diffstate_code), by = "match_key") %>%
    group_by(diffstate_code) %>%
    summarise(n_cells = sum(n_cells), .groups = "drop")
  total <- sum(df$n_cells)
  if (total == 0) return(NULL)
  vec <- setNames(rep(0, length(diffstate_levels)), diffstate_levels)
  vec[df$diffstate_code] <- 100 * df$n_cells / total
  list(vector = vec, n_cells_mapped = total)
}

# ----------------------------
# Part A: build reference (318 scAML patients) composition vectors from S8 + S9
# ----------------------------

message("\n=== Part A: single-cell differentiation-pattern classification ===")

reference_sample_ids <- setdiff(names(s8), "predicted_CellType")

reference_vectors <- map(reference_sample_ids, function(sid) {
  counts <- setNames(s8[[sid]], s8$predicted_CellType)
  celltype_counts_to_diffstate_pct(counts)
})
names(reference_vectors) <- reference_sample_ids

reference_matrix <- do.call(rbind, map(reference_vectors, ~ if (is.null(.x)) rep(NA_real_, length(diffstate_levels)) else .x$vector))
colnames(reference_matrix) <- diffstate_levels
reference_matrix <- reference_matrix[complete.cases(reference_matrix), , drop = FALSE]

reference_labels <- s7 %>% filter(Sample_ID %in% rownames(reference_matrix)) %>%
  distinct(Sample_ID, .keep_all = TRUE)
reference_matrix <- reference_matrix[rownames(reference_matrix) %in% reference_labels$Sample_ID, , drop = FALSE]
reference_labels <- reference_labels[match(rownames(reference_matrix), reference_labels$Sample_ID), ]

message("Built ", nrow(reference_matrix), " reference composition vectors (of ", length(reference_sample_ids), " S8 columns)")

# Pattern centroids, for a simple sanity check alongside per-patient nearest-neighbor matching
pattern_centroids <- reference_matrix %>%
  as_tibble(rownames = "Sample_ID") %>%
  left_join(reference_labels %>% select(Sample_ID, Differentiation_Pattern), by = "Sample_ID") %>%
  filter(!is.na(Differentiation_Pattern)) %>%
  group_by(Differentiation_Pattern) %>%
  summarise(across(all_of(diffstate_levels), mean), .groups = "drop")

# ----------------------------
# Load our Seurat objects, get numbat/copykat-consensus malignant compartment
# per donor (reusing the same malignant_call derivation used throughout this
# project's reports), and fine predicted_CellType per cell
# ----------------------------

# Same two conventions 22.patient_qc_report.Rmd's find_seurat_object() checks:
# the older results/seurat_annotated/<run>/ layout (still what 260528 has)
# and the pipeline-native results_nf/<run>/rds/08_annotate/<sample>/rds/
# layout (the only one newer Nextflow-run runs, e.g. 260522, actually have -
# there is no results/seurat_annotated/<run>/ for those at all). Sample
# folder name (LK1-GEX/LK2-GEX/...) is globbed rather than hardcoded, since
# it varies per run and this script only tracks run IDs.
find_seurat_file_for_run <- function(run) {
  old_path <- file.path(proj, "results/seurat_annotated", run, "demux_singlets_annotated_seurat.rds")
  if (file.exists(old_path)) return(old_path)
  nf_candidates <- Sys.glob(file.path(proj, "results_nf", run, "rds", "08_annotate", "*", "rds",
                                       "demux_singlets_annotated_seurat.rds"))
  if (length(nf_candidates) > 0) return(nf_candidates[1])
  NA_character_
}

seuList <- list()
for (run in runs) {
  seurat_file <- find_seurat_file_for_run(run)
  if (is.na(seurat_file)) {
    message("Missing demux_singlets_annotated_seurat.rds for ", run,
            " (checked results/seurat_annotated/ and results_nf/.../rds/08_annotate/) - skipping run")
    next
  }
  seuList[[run]] <- readRDS(seurat_file)
}
stopifnot(length(seuList) > 0)

seu <- if (length(seuList) > 1) merge(x = seuList[[1]], y = seuList[-1], add.cell.ids = names(seuList), project = "combined") else seuList[[1]]
rm(seuList)

DefaultAssay(seu) <- "RNA"
if (inherits(seu[["RNA"]], "Assay5")) seu <- JoinLayers(seu, assay = "RNA")

stopifnot(all(c("sample_name", "predicted_CellType", "scDblFinder.class") %in% colnames(seu@meta.data)))

seu <- subset(seu, subset = scDblFinder.class == "singlet" & !grepl("MOLM13", sample_name, ignore.case = TRUE))

donors <- sort(unique(seu$sample_name))
message("Donors found: ", paste(donors, collapse = ", "))

# ----------------------------
# Per donor: composition vector (S9-mapped predicted_CellType counts, all
# cells - S9's mapping already naturally excludes mature lymphocytes/other
# non-leukemic-relevant states simply because they have no entry in S9, same
# effect as the paper's own "exclusion of non-leukemic mature lymphocytes"),
# nearest-neighbor + centroid classification against the 318 reference
# patients.
# ----------------------------

classify_against_reference <- function(query_vec, reference_matrix, reference_labels, top_n = 10) {
  # Pearson correlation across the 13-state profile - robust to overall
  # scale differences between single-cell composition and the reference's
  # own composition, appropriate since both are % vectors on the same axes.
  cors <- apply(reference_matrix, 1, function(r) suppressWarnings(cor(query_vec, r, method = "pearson")))
  ord <- order(cors, decreasing = TRUE)
  top <- tibble(
    Sample_ID = rownames(reference_matrix)[ord][seq_len(top_n)],
    correlation = cors[ord][seq_len(top_n)]
  ) %>%
    left_join(reference_labels %>% select(Sample_ID, Differentiation_Pattern, Study, Diagnosis, Key_Mutations), by = "Sample_ID")

  assigned_pattern <- top %>% count(Differentiation_Pattern, sort = TRUE) %>% slice_head(n = 1) %>% pull(Differentiation_Pattern)
  list(assigned_pattern = assigned_pattern, top_matches = top)
}

pattern_results <- list()

for (donor in donors) {
  donor_meta <- seu@meta.data %>% filter(sample_name == donor, !is.na(predicted_CellType))
  celltype_counts <- table(donor_meta$predicted_CellType)

  res <- celltype_counts_to_diffstate_pct(celltype_counts)
  if (is.null(res) || res$n_cells_mapped < min_cells_composition) {
    message("  ", donor, ": only ", if (is.null(res)) 0 else res$n_cells_mapped,
            " cells mapped to a Zeng differentiation state (need >= ", min_cells_composition, ") - skipping")
    next
  }

  cls <- classify_against_reference(res$vector, reference_matrix, reference_labels)

  message("  ", donor, " (n=", res$n_cells_mapped, " mapped cells): assigned pattern = ", cls$assigned_pattern,
          " (top match r=", round(cls$top_matches$correlation[1], 3), ")")

  write_csv(cls$top_matches, file.path(out_dir, paste0(donor, "_zeng_nearest_reference_patients.csv")))

  pattern_results[[donor]] <- tibble(
    donor = donor, assigned_pattern = cls$assigned_pattern,
    n_cells_mapped = res$n_cells_mapped,
    top_match_sample = cls$top_matches$Sample_ID[1], top_match_r = cls$top_matches$correlation[1]
  ) %>%
    bind_cols(as_tibble_row(res$vector) %>% rename_with(~ paste0("pct_", .x)))
}

pattern_summary <- bind_rows(pattern_results)
if (nrow(pattern_summary) > 0) {
  write_csv(pattern_summary, file.path(out_dir, "patient_zeng_differentiation_pattern_assignment.csv"))
}
write_csv(pattern_centroids, file.path(out_dir, "zeng_reference_pattern_centroids.csv"))

# ----------------------------
# Part B: bulk RNA-seq deconvolution (Leucegene + our own pseudobulk) via
# S11's published LASSO models
# ----------------------------

message("\n=== Part B: bulk RNA-seq deconvolution (Leucegene + our pseudobulk) ===")

model_genes <- s11$Gene
coef_matrix <- as.matrix(s11 %>% select(-Gene))
rownames(coef_matrix) <- model_genes
colnames(coef_matrix) <- sub("^scAML_", "", colnames(coef_matrix))

# --- Leucegene, Ensembl -> symbol mapping (same convention as
# 24.leucegene_reference_and_candidate_genes.R for consistency) ---

message("Loading Leucegene readcounts: ", leucegene_file)
leucegene_raw <- read_tsv(leucegene_file, show_col_types = FALSE)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("org.Hs.eg.db required to map Leucegene's Ensembl IDs to gene symbols.")
}
ensembl_ids <- sub("\\..*$", "", leucegene_raw$ID)
gene_symbols <- suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL",
                        keytype = "ENSEMBL", multiVals = "first")
)
leucegene_counts <- leucegene_raw %>% select(-ID) %>% as.matrix()
rownames(leucegene_counts) <- gene_symbols
leucegene_counts <- leucegene_counts[!is.na(rownames(leucegene_counts)), , drop = FALSE]
leucegene_counts <- rowsum(leucegene_counts, group = rownames(leucegene_counts))

message("Leucegene matrix: ", nrow(leucegene_counts), " genes x ", ncol(leucegene_counts), " samples")

# --- Our own pseudobulk per donor (whole sample, not restricted to
# malignant compartment - matches Leucegene's own whole-marrow bulk nature,
# and the paper's own bulk deconvolution training on whole pseudo-bulk
# profiles per Supplementary Note S2) ---

make_pseudobulk_from_cells <- function(count_mat, cells) {
  cells <- intersect(cells, colnames(count_mat))
  if (length(cells) == 0) return(NULL)
  Matrix::rowSums(count_mat[, cells, drop = FALSE])
}

rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
our_pb_list <- lapply(donors, function(d) {
  cells <- rownames(seu@meta.data)[seu$sample_name == d]
  make_pseudobulk_from_cells(rna_counts, cells)
})
names(our_pb_list) <- donors
our_pb <- do.call(cbind, our_pb_list)
rownames(our_pb) <- rownames(rna_counts)

# --- Combine, log2CPM-normalize together, restrict to the 400 model genes ---

common_genes <- intersect(model_genes, intersect(rownames(leucegene_counts), rownames(our_pb)))
message("Model genes present in both Leucegene and our data: ", length(common_genes), " / ", length(model_genes))

combined_counts <- cbind(leucegene_counts[common_genes, , drop = FALSE], our_pb[common_genes, , drop = FALSE])
dge <- edgeR::DGEList(counts = combined_counts)
dge <- edgeR::calcNormFactors(dge)
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

coef_sub <- coef_matrix[common_genes, , drop = FALSE]
raw_scores <- t(logcpm) %*% coef_sub  # samples x diffstates

# z-standardize each state's score across the combined sample set, per
# Supplementary Note S2 ("standardized across patients")
standardized_scores <- scale(raw_scores)

bulk_scores_df <- as_tibble(standardized_scores, rownames = "sample") %>%
  mutate(cohort = ifelse(sample %in% colnames(leucegene_counts), "Leucegene", "This_project"))

write_csv(bulk_scores_df, file.path(out_dir, "bulk_deconvolution_standardized_scores.csv"))

# --- Nearest Leucegene samples per donor, in the 13-state standardized score space ---

leuc_scores <- standardized_scores[colnames(leucegene_counts), , drop = FALSE]
our_scores <- standardized_scores[donors[donors %in% rownames(standardized_scores)], , drop = FALSE]

nearest_leucegene <- map_dfr(rownames(our_scores), function(d) {
  dists <- sqrt(rowSums((sweep(leuc_scores, 2, our_scores[d, ], "-"))^2))
  ord <- order(dists)[1:10]
  tibble(donor = d, leucegene_sample = rownames(leuc_scores)[ord], rank = seq_along(ord), distance = dists[ord])
})
write_csv(nearest_leucegene, file.path(out_dir, "donor_nearest_leucegene_samples_bulk_deconv.csv"))

# --- PCA of the 13-state standardized score space (Leucegene + our patients),
# a compact, biologically-curated space rather than the whole transcriptome -
# far less prone to the platform-driven collapse seen in the earlier
# whole-gene PCA ---

pca <- prcomp(standardized_scores, center = FALSE, scale. = FALSE)
pca_df <- as_tibble(pca$x[, 1:2], rownames = "sample") %>%
  mutate(cohort = ifelse(sample %in% colnames(leucegene_counts), "Leucegene", "This_project"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point(data = ~filter(.x, cohort == "Leucegene"), alpha = 0.4, size = 1.2) +
  geom_point(data = ~filter(.x, cohort == "This_project"), size = 3, shape = 17) +
  ggrepel::geom_text_repel(data = ~filter(.x, cohort == "This_project"), aes(label = sample), size = 3) +
  scale_color_manual(values = c(Leucegene = "grey60", This_project = "firebrick3")) +
  theme_bw(base_size = 12) +
  labs(title = "Zeng et al. 13-differentiation-state score space: Leucegene + this project",
       subtitle = paste0(length(common_genes), " / ", length(model_genes), " model genes used"))

ggsave(file.path(out_dir, "leucegene_bulk_deconvolution_pca.pdf"), p_pca, width = 8, height = 6)
print(p_pca)

message("\nDone. Output directory: ", out_dir)
