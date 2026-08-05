#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# CITE-seq benchmark - step 06 of 07
#
# Diagnoses why ADT reads come back unmapped: observed barcode counts, edit distance to the known barcodes, an offset scan, and a reverse-complement check.
#
# Frozen for the lab archive 2026-08-04 from scripts/7g.CITE_unmapped_QC.R (mtime 2026-07-07).
# md5 of the original: b50c60ed5b28b9af56ccb70d3e0ad0d2
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

# ============================================================================
# investigate_adt_unmapped_reads.R
#
# Purpose:
#   Diagnose why ADT / CITE-seq reads are reported as "unmapped".
#
# Key point:
#   If you provide the full ADT FASTQ, this diagnoses barcode structure/reference
#   globally. If you provide an unmapped-only FASTQ, this directly diagnoses the
#   failed reads.
#
# Outputs:
#   - observed_barcode_counts.csv
#   - observed_vs_known_distance.csv
#   - offset_scan_summary.csv
#   - reverse_complement_distance.csv
#   - base_composition_by_position.csv
#   - diagnostic_summary.txt
#   - plots/*.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(stringdist)
})

# ----------------------------
# User parameters
# ----------------------------

input_mode <- "full_fastq"
# input_mode <- "unmapped_fastq"

fastq_r2_files <- c(
  "../raw/260528_VH01624_464_222K7VKNX/LK2-TotalSeqA-ADT_S3_R2_001.fastq.gz"
)

feature_reference_file <- "../annotation/adt_tags.csv"

out_dir <- "../results/adt_unmapped_investigator"

barcode_length <- 15L
max_reads <- 3.32e7       # set Inf for full file
chunk_size <- 1e6
top_n_observed <- 1000L

near_match_distance <- 1L
suspicious_distance <- 3L

offset_starts <- 1:30
checkpoint_every <- 5e6

# ----------------------------
# Helpers
# ----------------------------

stop_if_missing <- function(x, label) {
  missing <- x[!file.exists(x)]
  if (length(missing) > 0) {
    stop(label, " not found:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
}

revcomp <- function(x) {
  chartr("ACGTN", "TGCAN", stringi::stri_reverse(x))
}

read_known_barcodes <- function(feature_reference_file) {
  message("Reading feature reference: ", feature_reference_file)
  
  ref <- fread(feature_reference_file, header = FALSE)
  
  if (ncol(ref) < 2) {
    stop("Expected at least two columns: barcode sequence and antibody name.")
  }
  
  tibble(
    sequence = toupper(as.character(ref[[1]])),
    antibody = as.character(ref[[2]])
  ) %>%
    filter(!is.na(sequence), !is.na(antibody), sequence != "", antibody != "") %>%
    mutate(antibody = make.unique(antibody))
}

count_barcodes_from_fastq <- function(
    fastq_files,
    barcode_length,
    barcode_start = 1L,
    max_reads = Inf,
    chunk_size = 1e6,
    checkpoint_file = NULL,
    checkpoint_every = 5e6
) {
  all_counts <- integer()
  total_reads <- 0
  
  for (fq in fastq_files) {
    message("Reading FASTQ: ", fq)
    con <- gzfile(fq, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    
    repeat {
      remaining <- max_reads - total_reads
      if (remaining <= 0) break
      
      reads_this_chunk <- min(chunk_size, remaining)
      lines <- readLines(con, n = reads_this_chunk * 4L)
      if (length(lines) == 0L) break
      
      n_complete <- floor(length(lines) / 4L)
      lines <- lines[seq_len(n_complete * 4L)]
      
      seqs <- toupper(lines[seq(2, length(lines), by = 4)])
      bcs <- substr(seqs, barcode_start, barcode_start + barcode_length - 1L)
      
      bcs <- bcs[nchar(bcs) == barcode_length]
      
      tab <- table(bcs)
      idx <- names(tab)
      
      old <- all_counts[idx]
      old[is.na(old)] <- 0L
      all_counts[idx] <- old + as.integer(tab)
      
      total_reads <- total_reads + length(seqs)
      message("  processed reads: ", format(total_reads, scientific = TRUE))
      
      if (!is.null(checkpoint_file) && total_reads %% checkpoint_every < chunk_size) {
        saveRDS(
          list(
            counts = all_counts,
            total_reads = total_reads,
            barcode_start = barcode_start,
            barcode_length = barcode_length
          ),
          checkpoint_file
        )
      }
    }
    
    try(close(con), silent = TRUE)
    if (total_reads >= max_reads) break
  }
  
  observed <- tibble(
    observed_sequence = names(all_counts),
    read_count = as.integer(all_counts)
  ) %>%
    arrange(desc(read_count)) %>%
    mutate(
      rank = row_number(),
      pct_reads = 100 * read_count / sum(read_count)
    )
  
  attr(observed, "total_fastq_reads_processed") <- total_reads
  observed
}

nearest_known_distance <- function(observed, known, seq_col = "observed_sequence") {
  obs_seq <- observed[[seq_col]]
  known_seq <- known$sequence
  
  dmat <- stringdistmatrix(obs_seq, known_seq, method = "lv")
  min_idx <- apply(dmat, 1, which.min)
  min_dist <- dmat[cbind(seq_along(obs_seq), min_idx)]
  
  observed %>%
    mutate(
      nearest_known_sequence = known_seq[min_idx],
      nearest_antibody = known$antibody[min_idx],
      edit_distance = as.integer(min_dist),
      exact_match = edit_distance == 0,
      near_match = edit_distance <= near_match_distance
    )
}

offset_scan_fastq <- function(
    fastq_files,
    known,
    barcode_length,
    offset_starts,
    max_reads = 1e6,
    chunk_size = 1e6
) {
  known_set <- known$sequence
  results <- list()
  
  for (start in offset_starts) {
    message("Offset scan: R2 bases ", start, "-", start + barcode_length - 1L)
    
    obs <- count_barcodes_from_fastq(
      fastq_files = fastq_files,
      barcode_length = barcode_length,
      barcode_start = start,
      max_reads = max_reads,
      chunk_size = chunk_size
    )
    
    exact_reads <- sum(obs$read_count[obs$observed_sequence %in% known_set])
    total_reads <- sum(obs$read_count)
    
    top_obs <- obs %>% slice_head(n = min(500, nrow(obs)))
    top_dist <- nearest_known_distance(top_obs, known)
    
    results[[as.character(start)]] <- tibble(
      barcode_start = start,
      barcode_end = start + barcode_length - 1L,
      reads_checked = total_reads,
      exact_match_reads = exact_reads,
      exact_match_pct = 100 * exact_reads / total_reads,
      top500_near_match_pct = 100 * sum(top_dist$read_count[top_dist$edit_distance <= 1]) / sum(top_dist$read_count),
      median_top500_distance = weighted_median(top_dist$edit_distance, top_dist$read_count)
    )
  }
  
  bind_rows(results)
}

weighted_median <- function(x, w) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

base_composition <- function(observed, max_rank = 1000L) {
  x <- observed %>%
    slice_head(n = min(max_rank, nrow(.))) %>%
    select(observed_sequence, read_count)
  
  seqs <- x$observed_sequence
  weights <- x$read_count
  len <- nchar(seqs[1])
  
  out <- list()
  
  for (pos in seq_len(len)) {
    bases <- substr(seqs, pos, pos)
    tab <- tapply(weights, bases, sum)
    
    out[[pos]] <- tibble(
      position = pos,
      base = c("A", "C", "G", "T", "N"),
      read_count = as.numeric(tab[c("A", "C", "G", "T", "N")])
    ) %>%
      mutate(read_count = replace_na(read_count, 0),
             pct = 100 * read_count / sum(read_count))
  }
  
  bind_rows(out)
}

make_plots <- function(matches, offset_scan, base_comp, out_dir) {
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  p1 <- ggplot(matches, aes(x = edit_distance, weight = read_count)) +
    geom_bar() +
    theme_bw() +
    labs(
      title = "Distance of observed ADT barcode sequences to nearest known barcode",
      subtitle = "Weighted by read count",
      x = "Levenshtein edit distance",
      y = "Read count"
    )
  ggsave(file.path(plot_dir, "edit_distance_histogram.pdf"), p1, width = 7, height = 4)
  
  top_plot <- matches %>%
    slice_head(n = min(50, nrow(.))) %>%
    mutate(label = paste0(rank, ": ", observed_sequence),
           label = factor(label, levels = rev(label)))
  
  p2 <- ggplot(top_plot, aes(x = label, y = read_count, fill = factor(edit_distance))) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Top observed ADT barcode sequences",
      x = NULL,
      y = "Read count",
      fill = "Edit distance"
    )
  ggsave(file.path(plot_dir, "top_observed_barcodes.pdf"), p2, width = 8, height = 10)
  
  p3 <- ggplot(offset_scan, aes(x = barcode_start, y = exact_match_pct)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(
      title = "R2 offset scan",
      subtitle = "Where does the known ADT barcode best match?",
      x = "R2 barcode start position",
      y = "% reads exactly matching known ADT barcodes"
    )
  ggsave(file.path(plot_dir, "offset_scan_exact_match_pct.pdf"), p3, width = 7, height = 4)
  
  p4 <- ggplot(base_comp, aes(x = position, y = pct, fill = base)) +
    geom_col(position = "stack") +
    theme_bw() +
    labs(
      title = "Base composition of observed ADT barcode sequences",
      subtitle = "Weighted by read count among top observed sequences",
      x = "Barcode position",
      y = "% reads",
      fill = "Base"
    )
  ggsave(file.path(plot_dir, "base_composition_by_position.pdf"), p4, width = 7, height = 4)
}

write_summary <- function(
    input_mode,
    observed,
    known,
    matches,
    rc_matches,
    offset_scan,
    suspicious,
    out_dir
) {
  total_reads <- attr(observed, "total_fastq_reads_processed")
  
  total_top_reads <- sum(matches$read_count)
  exact_reads <- sum(matches$read_count[matches$edit_distance == 0])
  near_reads <- sum(matches$read_count[matches$edit_distance <= near_match_distance])
  suspicious_reads <- sum(suspicious$read_count)
  
  rc_exact_reads <- sum(rc_matches$read_count[rc_matches$edit_distance == 0])
  rc_near_reads <- sum(rc_matches$read_count[rc_matches$edit_distance <= near_match_distance])
  
  best_offset <- offset_scan %>% arrange(desc(exact_match_pct)) %>% slice_head(n = 1)
  top1 <- matches %>% slice_head(n = 1)
  
  likely <- c()
  
  if (best_offset$barcode_start != 1 && best_offset$exact_match_pct > 50) {
    likely <- c(likely, paste0(
      "Barcode appears shifted: best R2 start is ",
      best_offset$barcode_start,
      " with ",
      round(best_offset$exact_match_pct, 2),
      "% exact matches."
    ))
  }
  
  if ((100 * exact_reads / total_top_reads) > 90) {
    likely <- c(likely, "Most abundant observed sequences match the known ADT reference. Wrong antibody barcode reference is unlikely.")
  }
  
  if ((100 * suspicious_reads / total_top_reads) > 20) {
    likely <- c(likely, "Many high-abundance sequences are far from known barcodes. Possible wrong panel, wrong barcode length, wrong R2 structure, or non-ADT contamination.")
  }
  
  if ((100 * rc_exact_reads / total_top_reads) > (100 * exact_reads / total_top_reads) + 20) {
    likely <- c(likely, "Reverse-complement orientation fits much better than forward orientation.")
  }
  
  if (length(likely) == 0) {
    likely <- "No single dominant failure mode detected from barcode sequence structure."
  }
  
  txt <- c(
    "ADT unmapped-read investigator summary",
    "======================================",
    paste0("Input mode: ", input_mode),
    "",
    "Important interpretation:",
    ifelse(
      input_mode == "unmapped_fastq",
      "  These results describe reads already classified as unmapped.",
      "  These results describe all reads in the supplied FASTQ, not only reads classified as unmapped."
    ),
    "",
    paste0("FASTQ reads processed: ", format(total_reads, big.mark = ",")),
    paste0("Unique observed ", barcode_length, "-bp sequences: ", format(nrow(observed), big.mark = ",")),
    paste0("Known ADT barcodes: ", nrow(known)),
    paste0("Top observed sequences analysed: ", nrow(matches)),
    "",
    "Forward orientation, R2 bases 1-15:",
    paste0("  Exact-match reads: ", round(100 * exact_reads / total_top_reads, 2), "%"),
    paste0("  Near-match reads, edit distance <= ", near_match_distance, ": ", round(100 * near_reads / total_top_reads, 2), "%"),
    paste0("  Suspicious reads, edit distance >= ", suspicious_distance, ": ", round(100 * suspicious_reads / total_top_reads, 2), "%"),
    "",
    "Reverse-complement test:",
    paste0("  Exact-match reads after reverse-complement: ", round(100 * rc_exact_reads / total_top_reads, 2), "%"),
    paste0("  Near-match reads after reverse-complement: ", round(100 * rc_near_reads / total_top_reads, 2), "%"),
    "",
    "Best R2 offset scan:",
    paste0("  Best barcode start: ", best_offset$barcode_start),
    paste0("  Best barcode interval: R2 bases ", best_offset$barcode_start, "-", best_offset$barcode_end),
    paste0("  Exact-match reads at best offset: ", round(best_offset$exact_match_pct, 2), "%"),
    "",
    "Most abundant observed sequence:",
    paste0("  Sequence: ", top1$observed_sequence),
    paste0("  Read count: ", format(top1$read_count, big.mark = ",")),
    paste0("  % of observed barcode reads: ", round(top1$pct_reads, 2), "%"),
    paste0("  Nearest antibody: ", top1$nearest_antibody),
    paste0("  Nearest known sequence: ", top1$nearest_known_sequence),
    paste0("  Edit distance: ", top1$edit_distance),
    "",
    "Likely interpretation:",
    paste0("  - ", likely),
    "",
    "How to interpret true unmapped-only results:",
    "  - Mostly edit distance 0: reads are not unmapped because of ADT barcode sequence; check cell barcode/UMI parsing or software settings.",
    "  - Mostly edit distance 1: sequencing errors or too-strict mismatch threshold.",
    "  - Mostly high edit distance: wrong ADT reference, wrong barcode length, shifted R2 structure, or wrong read file.",
    "  - Reverse-complement much better: orientation problem.",
    "  - Offset other than 1 much better: barcode start position/adapter trimming problem."
  )
  
  writeLines(txt, file.path(out_dir, "diagnostic_summary.txt"))
  cat(paste(txt, collapse = "\n"), "\n")
}

# ----------------------------
# Run
# ----------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing(fastq_r2_files, "FASTQ R2 file")
stop_if_missing(feature_reference_file, "Feature reference file")

known <- read_known_barcodes(feature_reference_file)
write.csv(known, file.path(out_dir, "known_adt_barcodes.csv"), row.names = FALSE)

checkpoint_file <- file.path(out_dir, "observed_barcode_counts_checkpoint.rds")

observed <- count_barcodes_from_fastq(
  fastq_files = fastq_r2_files,
  barcode_length = barcode_length,
  barcode_start = 1L,
  max_reads = max_reads,
  chunk_size = chunk_size,
  checkpoint_file = checkpoint_file,
  checkpoint_every = checkpoint_every
)

saveRDS(observed, file.path(out_dir, "observed_barcode_counts.rds"))
write.csv(observed, file.path(out_dir, "observed_barcode_counts.csv"), row.names = FALSE)

top_observed <- observed %>% slice_head(n = top_n_observed)

matches <- nearest_known_distance(top_observed, known)
write.csv(matches, file.path(out_dir, "observed_vs_known_distance.csv"), row.names = FALSE)

rc_observed <- top_observed %>%
  mutate(observed_sequence_rc = revcomp(observed_sequence)) %>%
  select(-observed_sequence) %>%
  rename(observed_sequence = observed_sequence_rc)

rc_matches <- nearest_known_distance(rc_observed, known)
write.csv(rc_matches, file.path(out_dir, "reverse_complement_distance.csv"), row.names = FALSE)

offset_scan <- offset_scan_fastq(
  fastq_files = fastq_r2_files,
  known = known,
  barcode_length = barcode_length,
  offset_starts = offset_starts,
  max_reads = min(max_reads, 1e6),
  chunk_size = min(chunk_size, 1e6)
)
write.csv(offset_scan, file.path(out_dir, "offset_scan_summary.csv"), row.names = FALSE)

base_comp <- base_composition(observed, max_rank = top_n_observed)
write.csv(base_comp, file.path(out_dir, "base_composition_by_position.csv"), row.names = FALSE)

suspicious <- matches %>%
  filter(edit_distance >= suspicious_distance) %>%
  arrange(desc(read_count))
write.csv(suspicious, file.path(out_dir, "suspected_missing_or_unmapped_sequences.csv"), row.names = FALSE)

make_plots(matches, offset_scan, base_comp, out_dir)

write_summary(
  input_mode = input_mode,
  observed = observed,
  known = known,
  matches = matches,
  rc_matches = rc_matches,
  offset_scan = offset_scan,
  suspicious = suspicious,
  out_dir = out_dir
)

message("\nDone. Outputs written to: ", out_dir)