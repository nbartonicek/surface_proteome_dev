#!/usr/bin/env Rscript

# Compares three isoform annotations structurally (by genomic coordinates,
# not ID): AMLisoDB.1.0 (Shi et al., FLAIR+TALON+StringTie novel calls),
# the novel long-read isoforms from aml_data.gtf.gz (Miller et al.,
# ESPRESSO calls, source="annotated_isoform"), and GENCODE/Ensembl (the
# reference annotation that's ALSO concatenated into aml_data.gtf.gz under
# source %in% ensembl/ensembl_havana/havana/mirbase - confirmed this
# matches the paper's "206,601 known Ensembl transcripts" figure exactly
# when unfiltered by chromosome).
#
# Restricted to primary chromosomes (chr1-22, chrX, chrY, chrM) - all
# three carry scaffold/decoy contigs under different naming conventions,
# which would only add naming-mismatch noise.
#
# Parsing/grouping uses data.table (fread + grouped shift()), not
# rtracklayer::import() + per-transcript lapply() over GRanges (the first
# version of this script took >20 min and was killed still running).
# GenomicRanges is still used for the interval ops (reduce/intersect/
# findOverlaps), which are C-level and fast regardless of transcript count.
#
# Outputs (in out_dir):
#   - summary_all_pairs.csv: the full pairwise metric table (transcript/
#     gene counts, exonic bp coverage + Jaccard, splice junction/site
#     sharing, any-overlap transcripts, exact structural matches) for all
#     3 pairs (AMLisoDB-Miller, AMLisoDB-GENCODE, Miller-GENCODE)
#   - venn_splice_junctions.pdf: 3-way Venn of distinct (chr,strand,intron
#     start,intron end) junctions
#   - venn_transcribed_loci.pdf: 3-way Venn of "transcribed loci" - all 3
#     files' transcript spans pooled and merged into loci, then each file
#     counted as touching a locus if >=1 of its transcripts overlaps it
#     (this is the natural 3-way generalization of the pairwise
#     "n_transcripts_any_overlap" metric, since transcript IDs aren't
#     comparable across differently-sourced annotations)
#   - venn_exonic_footprint_bp.pdf: 3-way Venn of exonic genome coverage
#     (bp), i.e. "percentage of genome covered in common"
#   These three were picked as the most informative for a 3-way view;
#   splice sites and exact-structural-match are still in the pairwise CSV
#   but not Venn'd (site-level Venns are largely redundant with junctions,
#   and exact-match counts are too small/noisy to read well as a Venn).
#
# Needs data.table + GenomicRanges + ggVennDiagram + VennDiagram. Run with
# an R that has these installed (not necessarily the cluster's mgatk/r_env
# conda envs).
#
# Usage:
#   Rscript compare_isoform_annotations.R [gtf_a] [gtf_b] [out_dir]
#   Defaults: AMLisoDB.1.0.gtf.gz (Shi et al.) vs aml_data.gtf.gz (Miller et al.,
#   novel+GENCODE both extracted from this one file) in ../annotation/isoforms
#   relative to this script's directory.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(ggvenn)
  library(ggforce)
  library(ggplot2)
})

PALETTE3 <- c("#2a78d6", "#eb6834", "#1baf7a")

args <- commandArgs(trailingOnly = TRUE)
default_dir <- "../annotation/isoforms"
GTF_A <- if (length(args) >= 1) args[1] else file.path(default_dir, "AMLisoDB.1.0.gtf.gz")
GTF_B <- if (length(args) >= 2) args[2] else file.path(default_dir, "aml_data.gtf.gz")
OUT_DIR <- if (length(args) >= 3) args[3] else file.path(default_dir, "comparison_AMLisoDB_vs_aml_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LABEL_AML <- "AMLisoDB_Shi"
LABEL_MILLER <- "aml_data_novel_Miller"
LABEL_GENCODE <- "GENCODE_ref"

PRIMARY_CHROMS <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
REFERENCE_SOURCES <- c("ensembl", "ensembl_havana", "havana", "mirbase",
                        "ensembl_havana_tagene", "havana_tagene")

# ---- fast GTF -> exon-only data.table, primary chroms only ----
read_exons_raw <- function(path) {
  message("Reading ", path, " ...")
  dt <- fread(cmd = paste("zcat <", shQuote(path), "| grep -v '^#'"),
              sep = "\t", header = FALSE, quote = "",
              col.names = c("seqnames", "source", "feature", "start", "end",
                            "score", "strand", "frame", "attribute"))
  dt <- dt[feature == "exon" & seqnames %in% PRIMARY_CHROMS]
  dt[, transcript_id := sub('.*transcript_id "([^"]*)".*', "\\1", attribute)]
  dt[, gene_id := sub('.*gene_id "([^"]*)".*', "\\1", attribute)]
  # strip Ensembl version suffix (AMLisoDB embeds it directly, e.g.
  # "ENSG00000175756.13"; aml_data's GENCODE rows keep gene_id unversioned
  # already and store the version separately in "gene_version", so this is
  # a no-op there) - safe against non-Ensembl gene_ids too (AMLisoDB's
  # "novelGene_733"/"AS:CEP104" labels have no trailing ".<digits>").
  dt[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
  dt[, attribute := NULL]
  setkey(dt, transcript_id, start)
  dt
}

# ---- splice junctions (grouped shift(), no per-transcript loop) ----
junctions_from_exons <- function(dt) {
  d <- copy(dt)[order(transcript_id, start)]
  d[, next_start := shift(start, type = "lead"), by = transcript_id]
  d <- d[!is.na(next_start)]
  unique(d[, .(seqnames, strand, istart = end + 1L, iend = next_start - 1L)])
}
sites_from_junctions_dt <- function(jdt) {
  unique(rbind(jdt[, .(seqnames, strand, pos = istart)],
               jdt[, .(seqnames, strand, pos = iend)]))
}
build_chains <- function(exons_dt) {
  d <- copy(exons_dt)[order(transcript_id, start)]
  d[, next_start := shift(start, type = "lead"), by = transcript_id]
  d <- d[!is.na(next_start), .(transcript_id, seqnames, strand, istart = end + 1L, iend = next_start - 1L)]
  d[order(transcript_id, istart), .(sig = paste0(seqnames[1], ":", strand[1], ":",
                                                  paste(istart, iend, sep = "-", collapse = ","))),
    by = transcript_id]
}

# ---- bundles all derived structures for one dataset ----
process_dataset <- function(dt, label) {
  span_dt <- dt[, .(seqnames = seqnames[1], strand = strand[1],
                     start = min(start), end = max(end)), by = transcript_id]
  span <- GRanges(span_dt$seqnames, IRanges(span_dt$start, span_dt$end), span_dt$strand)
  gr <- GRanges(dt$seqnames, IRanges(dt$start, dt$end), dt$strand)
  footprint <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
  junctions_dt <- junctions_from_exons(dt)
  junctions <- GRanges(junctions_dt$seqnames, IRanges(junctions_dt$istart, junctions_dt$iend), junctions_dt$strand)
  sites_dt <- sites_from_junctions_dt(junctions_dt)
  sites <- GRanges(sites_dt$seqnames, IRanges(sites_dt$pos, sites_dt$pos), sites_dt$strand)
  chains <- build_chains(dt)
  n_tx <- uniqueN(dt$transcript_id); n_gene <- uniqueN(dt$gene_id)
  message(label, ": ", n_tx, " transcripts, ", n_gene, " genes (", nrow(dt), " exon rows)")
  list(label = label, exons = dt, span_dt = span_dt, span = span, footprint = footprint,
       junctions_dt = junctions_dt, junctions = junctions, sites = sites,
       chains = chains, n_tx = n_tx, n_gene = n_gene)
}

message("Loading ", GTF_A, " (AMLisoDB - all novel, no reference contamination) ...")
exons_aml <- read_exons_raw(GTF_A)
message("Loading ", GTF_B, " (aml_data - splitting novel ESPRESSO calls from embedded GENCODE) ...")
exons_b_raw <- read_exons_raw(GTF_B)
exons_miller <- exons_b_raw[!(source %in% REFERENCE_SOURCES)]
exons_gencode <- exons_b_raw[source %in% REFERENCE_SOURCES]

P_AML <- process_dataset(exons_aml, LABEL_AML)
P_MILLER <- process_dataset(exons_miller, LABEL_MILLER)
P_GENCODE <- process_dataset(exons_gencode, LABEL_GENCODE)

# ---- pairwise summary metrics (same definitions as before) ----
n_shared <- function(x, y) length(GenomicRanges::intersect(x, y, ignore.strand = FALSE))

reciprocal_overlap_match <- function(x, y, min_frac = 0.9) {
  if (length(x) == 0 || length(y) == 0) return(0)
  hits <- findOverlaps(x, y)
  if (length(hits) == 0) return(0)
  ov <- width(pintersect(x[queryHits(hits)], y[subjectHits(hits)]))
  frac_x <- ov / width(x[queryHits(hits)]); frac_y <- ov / width(y[subjectHits(hits)])
  length(unique(queryHits(hits)[frac_x >= min_frac & frac_y >= min_frac]))
}

pairwise_summary <- function(P1, P2) {
  bp1 <- sum(width(P1$footprint)); bp2 <- sum(width(P2$footprint))
  bp_shared <- sum(width(GenomicRanges::intersect(P1$footprint, P2$footprint, ignore.strand = TRUE)))
  bp_union <- sum(width(GenomicRanges::reduce(c(P1$footprint, P2$footprint), ignore.strand = TRUE)))

  junc_shared <- n_shared(P1$junctions, P2$junctions)
  sites_shared <- n_shared(P1$sites, P2$sites)

  any_ov_1 <- length(unique(queryHits(findOverlaps(P1$span, P2$span))))
  any_ov_2 <- length(unique(queryHits(findOverlaps(P2$span, P1$span))))

  exact_multi_1 <- sum(P1$chains$sig %in% P2$chains$sig)
  exact_multi_2 <- sum(P2$chains$sig %in% P1$chains$sig)
  single_1 <- P1$span[!(P1$span_dt$transcript_id %in% P1$chains$transcript_id)]
  single_2 <- P2$span[!(P2$span_dt$transcript_id %in% P2$chains$transcript_id)]
  exact_single_1 <- reciprocal_overlap_match(single_1, single_2)
  exact_single_2 <- reciprocal_overlap_match(single_2, single_1)
  exact_1 <- exact_multi_1 + exact_single_1; exact_2 <- exact_multi_2 + exact_single_2

  tibble(
    comparison = paste0(P1$label, " vs ", P2$label),
    metric = c("n_transcripts", "n_genes", "exonic_bp_covered", "bp_shared_between_files",
               "pct_of_this_file_covered_by_other", "jaccard_exonic_footprint",
               "n_splice_junctions", "n_splice_junctions_shared", "pct_junctions_shared",
               "n_splice_sites", "n_splice_sites_shared", "pct_splice_sites_shared",
               "n_transcripts_any_overlap", "pct_transcripts_any_overlap",
               "n_transcripts_exact_structural_match", "pct_transcripts_exact_structural_match"),
    value_1 = c(P1$n_tx, P1$n_gene, bp1, bp_shared, round(100 * bp_shared / bp1, 2), round(bp_shared / bp_union, 4),
                length(P1$junctions), junc_shared, round(100 * junc_shared / length(P1$junctions), 2),
                length(P1$sites), sites_shared, round(100 * sites_shared / length(P1$sites), 2),
                any_ov_1, round(100 * any_ov_1 / P1$n_tx, 2), exact_1, round(100 * exact_1 / P1$n_tx, 2)),
    value_2 = c(P2$n_tx, P2$n_gene, bp2, bp_shared, round(100 * bp_shared / bp2, 2), round(bp_shared / bp_union, 4),
                length(P2$junctions), junc_shared, round(100 * junc_shared / length(P2$junctions), 2),
                length(P2$sites), sites_shared, round(100 * sites_shared / length(P2$sites), 2),
                any_ov_2, round(100 * any_ov_2 / P2$n_tx, 2), exact_2, round(100 * exact_2 / P2$n_tx, 2))
  ) %>% rename(!!P1$label := value_1, !!P2$label := value_2)
}

summary_all <- bind_rows(
  pairwise_summary(P_AML, P_MILLER),
  pairwise_summary(P_AML, P_GENCODE),
  pairwise_summary(P_MILLER, P_GENCODE)
)
write_csv(summary_all, file.path(OUT_DIR, "summary_all_pairs.csv"))
print(summary_all, n = Inf)

# ---- 3-way Venn 1: splice junctions (exact set membership) ----
junction_key <- function(jdt) paste(jdt$seqnames, jdt$strand, jdt$istart, jdt$iend, sep = ":")
venn_junctions <- list(junction_key(P_AML$junctions_dt), junction_key(P_MILLER$junctions_dt), junction_key(P_GENCODE$junctions_dt))
names(venn_junctions) <- c(LABEL_AML, LABEL_MILLER, LABEL_GENCODE)
p1 <- ggvenn(venn_junctions, fill_color = PALETTE3, fill_alpha = 0.55, stroke_color = "white",
             stroke_size = 1, set_name_size = 5, text_size = 4, show_percentage = TRUE) +
  labs(title = "Distinct splice junctions shared across annotations") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
ggsave(file.path(OUT_DIR, "venn_splice_junctions.pdf"), p1, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "venn_splice_junctions.png"), p1, width = 8, height = 7, dpi = 150)

# ---- 3-way Venn 2: transcribed loci (generalizes "any overlap" transcripts) ----
all_spans <- c(P_AML$span, P_MILLER$span, P_GENCODE$span)
combined_loci <- GenomicRanges::reduce(all_spans, ignore.strand = FALSE)
loci_touched <- function(span_x) unique(subjectHits(findOverlaps(span_x, combined_loci)))
venn_loci <- list(loci_touched(P_AML$span), loci_touched(P_MILLER$span), loci_touched(P_GENCODE$span))
names(venn_loci) <- c(LABEL_AML, LABEL_MILLER, LABEL_GENCODE)
p2 <- ggvenn(venn_loci, fill_color = PALETTE3, fill_alpha = 0.55, stroke_color = "white",
             stroke_size = 1, set_name_size = 5, text_size = 4, show_percentage = TRUE) +
  labs(title = "Genomic loci with >=1 overlapping transcript, by annotation",
       subtitle = paste0(length(combined_loci), " merged loci total (union of all 3 annotations' transcript spans, reduced) - ",
                          "NOTE: every Miller-novel locus is also touched by GENCODE (0 in that exclusive region); this reflects ",
                          "GENCODE's near-total genome coverage plus locus-merging being transitive, not that Miller's calls ",
                          "lack unannotated regions (the paper reports 13.3% of its novel transcripts as entirely intergenic)")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 8, margin = margin(t = 4)))
ggsave(file.path(OUT_DIR, "venn_transcribed_loci.pdf"), p2, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "venn_transcribed_loci.png"), p2, width = 8, height = 7, dpi = 150)

# ---- 3-way Venn 3: exonic genome coverage (bp) - a continuous quantity, not
# a list of discrete items, so it doesn't fit ggvenn's/ggVennDiagram's item-
# list API. Built directly as 3 circles (ggforce::geom_circle) at fixed,
# standard symmetric-Venn coordinates, labelled with the exact precomputed
# region sizes (Mb) - same visual language as the two ggvenn plots above,
# without needing to approximate bp totals as enumerable "items".
bp <- function(gr) sum(width(gr))
bp_a <- bp(P_AML$footprint); bp_m <- bp(P_MILLER$footprint); bp_g <- bp(P_GENCODE$footprint)
bp_am <- bp(GenomicRanges::intersect(P_AML$footprint, P_MILLER$footprint, ignore.strand = TRUE))
bp_ag <- bp(GenomicRanges::intersect(P_AML$footprint, P_GENCODE$footprint, ignore.strand = TRUE))
bp_mg <- bp(GenomicRanges::intersect(P_MILLER$footprint, P_GENCODE$footprint, ignore.strand = TRUE))
bp_amg <- bp(GenomicRanges::intersect(GenomicRanges::intersect(P_AML$footprint, P_MILLER$footprint, ignore.strand = TRUE),
                                       P_GENCODE$footprint, ignore.strand = TRUE))
mb <- function(x) sprintf("%.1f Mb", x / 1e6)

circles <- tibble(
  x = c(-0.8, 0.8, 0), y = c(0.55, 0.55, -0.7), r = 1.5,
  label = c(LABEL_AML, LABEL_MILLER, LABEL_GENCODE), fill = PALETTE3,
  name_x = c(-1.35, 1.35, 0), name_y = c(2.35, 2.35, -2.55)
)
region_labels <- tibble(
  x = c(-1.35, 1.35, 0,    0,    -0.95, 0.95, 0),
  y = c(1.05,  1.05, -1.9, 1.35, -0.55, -0.55, 0.0),
  label = c(mb(bp_a - bp_am - bp_ag + bp_amg), mb(bp_m - bp_am - bp_mg + bp_amg), mb(bp_g - bp_ag - bp_mg + bp_amg),
            mb(bp_am - bp_amg), mb(bp_ag - bp_amg), mb(bp_mg - bp_amg), mb(bp_amg))
)
p3 <- ggplot() +
  geom_circle(data = circles, aes(x0 = x, y0 = y, r = r, fill = label), color = "white", linewidth = 1, alpha = 0.55) +
  geom_text(data = circles, aes(x = name_x, y = name_y, label = label), size = 5) +
  geom_text(data = region_labels, aes(x = x, y = y, label = label), size = 4) +
  scale_fill_manual(values = setNames(PALETTE3, circles$label), guide = "none") +
  coord_fixed(clip = "off") + theme_void() +
  labs(title = "Exonic genome coverage shared across annotations") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)),
        plot.margin = margin(20, 20, 20, 20))
ggsave(file.path(OUT_DIR, "venn_exonic_footprint_bp.pdf"), p3, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "venn_exonic_footprint_bp.png"), p3, width = 8, height = 7, dpi = 150)

message("\nWrote outputs to: ", OUT_DIR)
message("  - summary_all_pairs.csv")
message("  - venn_splice_junctions.pdf")
message("  - venn_transcribed_loci.pdf")
message("  - venn_exonic_footprint_bp.pdf")
