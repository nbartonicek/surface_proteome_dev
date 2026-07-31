#!/usr/bin/env Rscript
# ============================================================================
# aml_scaffold_interactors.R
#
# R port of aml_scaffold_interactors.py, integrated with:
#   (A) IntAct interactor retrieval + intracellular-localisation flagging   [port]
#   (B) AlphaFold / ColabFold iPTM analysis + plotting  (your FCGR2A screen)
#   (C) three-way overlap with your single-cell surfaceome pipeline
#   (D) ColabFold multimer input generation, with your validated controls
#       (CRP positive; ACTB / BRCA1 / HSP90AA1 negatives) and SNRNP200-style
#       chunking baked in.
#
# The discovery logic is identical to the validated Python version; only the
# language changed. No R was available to execute this at authoring time, so
# sanity-check the first live run (especially the httr2 calls).
#
# Deps: tidyverse, data.table, httr2, digest, ggplot2  (Seurat only for part C)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(httr2)
  library(digest)
})

# ---------------------------------------------------------------------------
# SECTION 0 — CONFIG
# ---------------------------------------------------------------------------

# Bait scaffolds. FCGR2A included so the pipeline recovers its own positive
# control (FCGR2A <-> SNRNP200) end-to-end.
BAIT_GENES <- c("FCGR2A", "CD74", "CD47", "CD37", "ENG",
                "IL1RAP", "SEMA4D", "FLT3", "ITGA4")

TAXID_HUMAN  <- "9606"
UNIPROT_REST <- "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ACC  <- "https://rest.uniprot.org/uniprotkb"   # for .fasta

PSICQUIC_SERVICES <- c(
  IntAct = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/",
  MINT   = "https://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/"
)

ACCEPTED_INTERACTION_TERMS <- c("physical association", "direct interaction",
                                "association", "covalent binding")

INTRACELLULAR_TERMS <- c("nucleus", "nucleoplasm", "nucleolus", "cytoplasm",
                         "cytosol", "mitochond", "endoplasmic reticulum",
                         "spliceosome", "ribosom", "chromosome", "cytoskeleton",
                         "golgi", "peroxisome", "lysosome lumen")
SURFACE_TERMS <- c("cell membrane", "plasma membrane", "cell surface",
                   "secreted", "extracellular", "membrane raft", "cell junction")

# Controls for the AlphaFold step (edit if you swap them).
AF_POSITIVE_CONTROLS <- c("CRP")
AF_NEGATIVE_CONTROLS <- c("ACTB", "BRCA1", "HSP90AA1")

CACHE_DIR      <- ".scaffold_cache_R"
REQUEST_PAUSE  <- 0.34
PSICQUIC_PAGE  <- 500L
HTTP_TIMEOUT   <- 60
HTTP_RETRIES   <- 4L

# chunking (matches your SNRNP200 windows: width 600, step 450 -> 150 overlap)
CHUNK_WIDTH <- 600L
CHUNK_STEP  <- 450L

# ---------------------------------------------------------------------------
# SECTION 1 — HTTP with on-disk cache + retry
# ---------------------------------------------------------------------------

.cache_path <- function(url, query) {
  if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, showWarnings = FALSE)
  key <- paste0(url, "|", jsonlite::toJSON(query, auto_unbox = TRUE))
  file.path(CACHE_DIR, paste0(digest(key, algo = "sha1"), ".txt"))
}

http_get_text <- function(url, query = list(), use_cache = TRUE) {
  cp <- .cache_path(url, query)
  if (use_cache && file.exists(cp)) return(readr::read_file(cp))
  resp <- request(url) |>
    req_url_query(!!!query) |>
    req_timeout(HTTP_TIMEOUT) |>
    req_retry(max_tries = HTTP_RETRIES, backoff = \(i) 1.5 * i) |>
    req_error(is_error = \(resp) FALSE) |>     # handle status manually
    req_perform()
  if (resp_status(resp) != 200)
    stop("GET ", url, " -> HTTP ", resp_status(resp))
  txt <- resp_body_string(resp)
  writeLines(txt, cp)
  Sys.sleep(REQUEST_PAUSE)
  txt
}

# ---------------------------------------------------------------------------
# SECTION 2 — UniProt: gene -> accession, subcellular location, sequence
# ---------------------------------------------------------------------------

gene_to_accession <- function(gene) {
  txt <- http_get_text(UNIPROT_REST, list(
    query  = sprintf("gene_exact:%s AND organism_id:%s AND reviewed:true",
                     gene, TAXID_HUMAN),
    fields = "accession,gene_names", format = "tsv", size = 1))
  lines <- str_split(txt, "\n")[[1]] |> discard(~ .x == "")
  if (length(lines) < 2) return(NA_character_)
  str_split(lines[2], "\t")[[1]][1] |> str_trim()
}

fetch_localizations <- function(accessions) {
  out <- list()
  for (i in seq(1, length(accessions), by = 40)) {
    chunk <- accessions[i:min(i + 39, length(accessions))]
    q <- paste(sprintf("accession:%s", chunk), collapse = " OR ")
    txt <- http_get_text(UNIPROT_REST, list(
      query = q, fields = "accession,cc_subcellular_location",
      format = "tsv", size = length(chunk)))
    lines <- str_split(txt, "\n")[[1]][-1] |> discard(~ .x == "")
    for (ln in lines) {
      p <- str_split(ln, "\t")[[1]]
      out[[str_trim(p[1])]] <- if (length(p) > 1) str_trim(p[2]) else ""
    }
  }
  out
}

fetch_sequence <- function(acc) {
  txt <- http_get_text(sprintf("%s/%s.fasta", UNIPROT_ACC, acc), list())
  lines <- str_split(txt, "\n")[[1]]
  paste(lines[!str_starts(lines, ">")], collapse = "") |> str_remove_all("\\s")
}

classify_localization <- function(loc_text) {
  t <- str_to_lower(loc_text %||% "")
  c(has_intra   = any(str_detect(t, fixed(INTRACELLULAR_TERMS))),
    has_surface = any(str_detect(t, fixed(SURFACE_TERMS))))
}

# ---------------------------------------------------------------------------
# SECTION 3 — PSICQUIC / IntAct MITAB (2.5) fetch + parse
# ---------------------------------------------------------------------------

.first_uniprot <- function(field) {
  m <- str_match(field %||% "", "uniprotkb:([A-Z0-9]+(?:-\\d+)?)")[, 2]
  if (is.na(m)) return(NA_character_)
  str_split(m, "-")[[1]][1]
}
.taxid    <- function(field) str_match(field %||% "", "taxid:(-?\\d+)")[, 2]
.gene_nm  <- function(a, b) {
  for (f in c(a, b)) {
    m <- str_match(f %||% "", "uniprotkb:([^()|]+)\\(gene name\\)")[, 2]
    if (!is.na(m)) return(str_trim(m))
  }
  NA_character_
}
.clean_mi <- function(field) {
  m <- str_match(str_trim(field %||% ""), "\\(([^()]+)\\)\\s*$")[, 2]
  if (is.na(m)) str_trim(field %||% "") else m
}

parse_mitab <- function(text, bait_acc) {
  lines <- str_split(text, "\n")[[1]]
  lines <- lines[lines != "" & !str_starts(lines, "#")]
  recs <- list()
  for (ln in lines) {
    c <- str_split(ln, "\t")[[1]]
    if (length(c) < 15) next
    accA <- .first_uniprot(c[1]); accB <- .first_uniprot(c[2])
    if (is.na(accA) || is.na(accB)) next
    if (.taxid(c[10]) != TAXID_HUMAN || .taxid(c[11]) != TAXID_HUMAN) next
    itype <- str_to_lower(c[12] %||% "")
    if (!any(str_detect(itype, fixed(ACCEPTED_INTERACTION_TERMS)))) next
    if (accA == bait_acc && accB != bait_acc) {
      p_acc <- accB; p_al <- c[6]; p_alt <- c[4]
    } else if (accB == bait_acc && accA != bait_acc) {
      p_acc <- accA; p_al <- c[5]; p_alt <- c[3]
    } else next
    recs[[length(recs) + 1]] <- tibble(
      partner_acc  = p_acc,
      partner_gene = .gene_nm(p_al, p_alt),
      method       = .clean_mi(c[7]),
      pub          = str_trim(c[9] %||% ""),
      source       = .clean_mi(c[13]))
  }
  if (length(recs) == 0) return(tibble())
  bind_rows(recs)
}

query_intact <- function(bait_acc, services) {
  all <- list()
  for (nm in names(services)) {
    first <- 0L
    repeat {
      txt <- http_get_text(paste0(services[[nm]], "id:", bait_acc),
                           list(format = "tab25", firstResult = first,
                                maxResults = PSICQUIC_PAGE))
      rows <- str_split(txt, "\n")[[1]]
      rows <- rows[rows != "" & !str_starts(rows, "#")]
      recs <- parse_mitab(txt, bait_acc)
      if (nrow(recs) > 0) { recs$db <- nm; all[[length(all) + 1]] <- recs }
      if (length(rows) < PSICQUIC_PAGE) break
      first <- first + PSICQUIC_PAGE
    }
  }
  if (length(all) == 0) return(tibble())
  bind_rows(all)
}

# ---------------------------------------------------------------------------
# SECTION 4 — aggregate + candidate flag
# ---------------------------------------------------------------------------

aggregate_partners <- function(bait_gene, bait_acc, recs) {
  if (nrow(recs) == 0) return(tibble())
  recs |>
    group_by(partner_acc) |>
    summarise(
      partner_gene = first(na.omit(partner_gene)) %||% "",
      n_pubs     = n_distinct(pub[pub != ""]),
      n_methods  = n_distinct(method[method != ""]),
      source_dbs = paste(sort(unique(db)), collapse = ","),
      .groups = "drop") |>
    mutate(bait_gene = bait_gene, bait_acc = bait_acc, .before = 1)
}

# ---------------------------------------------------------------------------
# SECTION 5 — run interactor discovery over baits
# ---------------------------------------------------------------------------

run_interactor_discovery <- function(genes = BAIT_GENES,
                                     extra_dbs = FALSE,
                                     outdir = "scaffold_out") {
  dir.create(outdir, showWarnings = FALSE)
  services <- PSICQUIC_SERVICES["IntAct"]
  if (extra_dbs) services <- PSICQUIC_SERVICES[c("IntAct", "MINT")]
  
  frames <- list()
  for (g in genes) {
    message("[bait] ", g, " ...")
    acc <- gene_to_accession(g)
    if (is.na(acc)) { message("  ! no accession; skipping"); next }
    df <- aggregate_partners(g, acc, query_intact(acc, services))
    message("  ", acc, ": ", nrow(df), " unique human partners")
    if (nrow(df) > 0) frames[[length(frames) + 1]] <- df
  }
  combined <- bind_rows(frames)
  if (nrow(combined) == 0) stop("No interactors found.")
  
  loc <- fetch_localizations(sort(unique(combined$partner_acc)))
  fl  <- map(combined$partner_acc, ~ classify_localization(loc[[.x]] %||% ""))
  combined <- combined |>
    mutate(
      subcellular_location  = map_chr(partner_acc, ~ loc[[.x]] %||% ""),
      has_intracellular     = map_lgl(fl, "has_intra"),
      has_surface_secreted  = map_lgl(fl, "has_surface"),
      intracellular_only_candidate = has_intracellular & !has_surface_secreted) |>
    arrange(desc(intracellular_only_candidate), desc(n_pubs), desc(n_methods))
  
  write_csv(combined, file.path(outdir, "all_bait_interactors_annotated.csv"))
  combined |> group_by(bait_gene) |>
    group_walk(~ write_csv(.x, file.path(outdir,
                                         paste0("interactors_", .y$bait_gene, ".csv"))))
  cand <- filter(combined, intracellular_only_candidate)
  write_csv(cand, file.path(outdir, "CANDIDATES_intracellular_only.csv"))
  
  message("\nTop intracellular-only candidates:")
  print(head(cand[c("bait_gene", "partner_gene", "partner_acc",
                    "n_pubs", "n_methods", "subcellular_location")], 30))
  invisible(list(all = combined, candidates = cand))
}

# ---------------------------------------------------------------------------
# SECTION 6 — AlphaFold / ColabFold iPTM analysis + plot
# ---------------------------------------------------------------------------

# Parse target name: "O75643_SNRNP200_HUMAN_chunk4_1351-1950"
parse_af_target <- function(subfolder) {
  s <- str_remove_all(subfolder, '"')
  acc  <- str_match(s, "^([A-Z0-9]+)_")[, 2]
  gene <- str_match(s, "^[A-Z0-9]+_([A-Za-z0-9]+)_")[, 2]
  ch   <- str_match(s, "chunk\\d+_(\\d+)-(\\d+)")
  tibble(subfolder = subfolder, acc = acc, gene = gene,
         is_chunk = !is.na(ch[, 1]),
         chunk_start = as.integer(ch[, 2]),
         chunk_end   = as.integer(ch[, 3]))
}

analyse_alphafold <- function(csv_path,
                              positive = AF_POSITIVE_CONTROLS,
                              negative = AF_NEGATIVE_CONTROLS,
                              plot_file = "alphafold_iptm.pdf") {
  raw <- fread(csv_path) |> as_tibble()
  meta <- distinct(raw, subfolder) |>
    mutate(parse_af_target(subfolder)) |>
    tidyr::unnest(cols = everything())   # in case of list-cols; harmless otherwise
  
  ranked <- raw |>
    group_by(subfolder) |>
    summarise(n_models   = n(),
              best_iptm  = max(iPTM, na.rm = TRUE),
              mean_iptm  = mean(iPTM, na.rm = TRUE),
              sd_iptm    = sd(iPTM,  na.rm = TRUE),
              best_ptm   = max(pTM,  na.rm = TRUE),
              best_plddt = max(pLDDT, na.rm = TRUE), .groups = "drop") |>
    left_join(distinct(meta), by = "subfolder") |>
    mutate(control_type = case_when(
      gene %in% positive ~ "positive",
      gene %in% negative ~ "negative",
      TRUE               ~ "candidate"),
      label = ifelse(is_chunk,
                     paste0(gene, "_", chunk_start, "-", chunk_end), gene)) |>
    arrange(best_iptm)
  
  # interface localisation: best-scoring chunk per candidate gene
  interface <- ranked |>
    filter(is_chunk, control_type == "candidate") |>
    group_by(gene) |>
    slice_max(best_iptm, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(gene, interface_region = paste0(chunk_start, "-", chunk_end),
              interface_iptm = best_iptm)
  
  # NB: max-of-N-models inflates iPTM; mean_iptm is the conservative read.
  p <- ranked |>
    mutate(label = factor(label, levels = label[order(best_iptm)])) |>
    ggplot(aes(label, best_iptm, colour = control_type)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey70") +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(positive = "#D62728",
                                   negative = "grey55",
                                   candidate = "#1F3B99")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "AlphaFold interface confidence (best iPTM per target)",
         subtitle = "dashed line = 0.5; positive control should top the ranking",
         x = NULL, y = "iPTM", colour = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  ggsave(plot_file, p, width = 7, height = 4.5)
  message("Wrote ", plot_file)
  if (nrow(interface) > 0) {
    message("Interface localisation (best chunk):")
    print(interface)
  }
  list(ranked = ranked, interface = interface, plot = p)
}

# ---------------------------------------------------------------------------
# SECTION 7 — integration with the single-cell surfaceome pipeline
# ---------------------------------------------------------------------------

# (7a) RNA evidence from YOUR object: is a candidate expressed in the malignant
#      (Numbat CNV_aberrant) compartment, and enriched vs CNV_neutral?
#      Uses the numbat_malignancy column your script already creates.
rna_malignancy_evidence <- function(seu, genes,
                                    malig_col = "numbat_malignancy",
                                    assay = "RNA") {
  genes <- intersect(genes, rownames(seu))
  if (length(genes) == 0) { message("none of the genes are in the object"); return(tibble()) }
  m   <- Seurat::GetAssayData(seu, assay = assay, slot = "data")[genes, , drop = FALSE]
  grp <- seu[[malig_col]][, 1]
  keep <- grp %in% c("CNV_aberrant", "CNV_neutral")
  m <- m[, keep, drop = FALSE]; grp <- grp[keep]
  aber <- grp == "CNV_aberrant"; norm <- grp == "CNV_neutral"
  tibble(
    gene        = genes,
    pct_malig   = Matrix::rowMeans(m[, aber, drop = FALSE] > 0),
    pct_normal  = Matrix::rowMeans(m[, norm, drop = FALSE] > 0),
    mean_malig  = Matrix::rowMeans(m[, aber, drop = FALSE]),
    mean_normal = Matrix::rowMeans(m[, norm, drop = FALSE])) |>
    mutate(log2fc_malig_vs_normal = log2((mean_malig + 1e-9) / (mean_normal + 1e-9))) |>
    arrange(desc(log2fc_malig_vs_normal))
}

# (7b) The three-way overlap.
#   circle 1: intracellular interactors of baits            (candidates$partner_gene)
#   circle 2: proteins empirically detected on the AML surface
#             (Bordeleau/Knorr surface proteomics; supply a character vector)
#   circle 3: expressed/enriched in malignant cells in YOUR data (7a)
triangulate <- function(candidates, surface_detected_genes, rna_evidence,
                        min_log2fc = 0) {
  candidates |>
    distinct(bait_gene, partner_gene, partner_acc, n_pubs, n_methods) |>
    filter(partner_gene %in% surface_detected_genes) |>            # circle 2
    left_join(rna_evidence, by = c("partner_gene" = "gene")) |>    # circle 3
    filter(!is.na(log2fc_malig_vs_normal),
           log2fc_malig_vs_normal >= min_log2fc) |>
    arrange(desc(log2fc_malig_vs_normal), desc(n_pubs))
}

# ---------------------------------------------------------------------------
# SECTION 8 — generate ColabFold multimer inputs (with your controls)
# ---------------------------------------------------------------------------

.chunk_ranges <- function(len, width = CHUNK_WIDTH, step = CHUNK_STEP) {
  if (len <= width) return(list(c(1L, len)))
  starts <- seq(1L, len, by = step)
  map(starts, ~ c(.x, min(.x + width - 1L, len))) |>
    keep(~ .x[1] <= len) |> unique()
}

# Writes one FASTA per bait-candidate pair: ">bait__cand\nBAITSEQ:CANDSEQ"
# Large candidates are chunked; positive/negative controls are appended so every
# run is internally calibrated exactly like your FCGR2A screen.
build_colabfold_inputs <- function(bait_gene, candidate_genes,
                                   positive = AF_POSITIVE_CONTROLS,
                                   negative = AF_NEGATIVE_CONTROLS,
                                   outdir = "colabfold_inputs") {
  dir.create(outdir, showWarnings = FALSE)
  bait_acc <- gene_to_accession(bait_gene); bait_seq <- fetch_sequence(bait_acc)
  targets  <- unique(c(candidate_genes, positive, negative))
  
  manifest <- list()
  for (g in targets) {
    acc <- gene_to_accession(g); if (is.na(acc)) { message("skip ", g); next }
    seq <- fetch_sequence(acc)
    role <- if (g %in% positive) "positive" else
      if (g %in% negative) "negative" else "candidate"
    for (rg in .chunk_ranges(nchar(seq))) {
      s <- substr(seq, rg[1], rg[2])
      tag <- if (nchar(seq) > CHUNK_WIDTH)
        sprintf("%s_%s_%s_%d-%d", bait_gene, g, acc, rg[1], rg[2]) else
          sprintf("%s_%s_%s", bait_gene, g, acc)
      writeLines(c(paste0(">", tag), paste0(bait_seq, ":", s)),
                 file.path(outdir, paste0(tag, ".fasta")))
      manifest[[length(manifest) + 1]] <-
        tibble(fasta = paste0(tag, ".fasta"), bait = bait_gene, target = g,
               role = role, start = rg[1], end = rg[2])
    }
  }
  man <- bind_rows(manifest)
  write_csv(man, file.path(outdir, "manifest.csv"))
  message("Wrote ", nrow(man), " FASTA pairs -> ", outdir,
          "  (run: colabfold_batch ", outdir, "/ <results>/)")
  invisible(man)
}

# ---------------------------------------------------------------------------
# EXAMPLE WORKFLOW (uncomment to run)
# ---------------------------------------------------------------------------
res  <- run_interactor_discovery(extra_dbs = FALSE)
#
# af   <- analyse_alphafold("colabfold_results_summary_FCGR2A.csv")
#         # reproduces your plot + reports SNRNP200 interface region (~1351-1950)
#
# # circle 2: your surface-proteomics hits (Bordeleau/Knorr), as a vector:
# surf <- fread("surface_proteomics_hits.csv")$gene
# # circle 3: RNA evidence from your projected object (from script 3):
# rna  <- rna_malignancy_evidence(seu_sub, res$candidates$partner_gene)
# hits <- triangulate(res$candidates, surf, rna)
# print(hits)
#
# # feed the winners back into ColabFold, controls + chunking auto-included:
# build_colabfold_inputs("FCGR2A", head(hits$partner_gene, 10))
# ============================================================================