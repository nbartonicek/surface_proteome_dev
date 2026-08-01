# ------------------------------------------------------------------
# LK1 pilot run - step 03 of 18
#
# emptyDrops on the raw Cell Ranger matrix, FDR <= 0.01. Its barcode list is the whitelist CITE-seq-Count and cellSNP were run against.
#
# Frozen for the lab archive 2026-07-31 from scripts/2.emptydrops.R (mtime 2026-06-03).
# md5 of the original: 6e07834dc7d64ade1c65b49f577d8dac
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
#
# NOTE: `run` in this copy is the RE-SEQUENCED run (260522), because the
# script was re-pointed at it later. It ran on the pilot first - the
# 260423 outputs predate 260522 existing. Set `run` back to the pilot run
# to reproduce.
# ------------------------------------------------------------------

library(DropletUtils)
library(Matrix)
library(Seurat)

run <- "260522_VH01624_461_222JLJVNX"
sample_name <- "LK1-GEX"
raw_dir <- paste0("../results/cellranger_withbam/",run,"/",sample_name,"/outs/raw_feature_bc_matrix/")
out_dir <- paste0("../results/emptydrops/",run,"/",sample_name)

system(paste0("mkdir -p ", out_dir))

gex_counts <- Read10X(raw_dir)

set.seed(123)

e.out <- emptyDrops(gex_counts)

keep <- which(!is.na(e.out$FDR) & e.out$FDR <= 0.01)

filtered_counts <- gex_counts[, keep]

barcodes <- colnames(filtered_counts)

write.table(
  barcodes,
  file = paste0("../results/emptydrops/",run,"/",sample_name,"/emptydrops_barcodes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

plot_df <- data.frame(
  total = e.out$Total,
  fdr = e.out$FDR
)

plot_df <- plot_df %>%
  filter(!is.na(fdr), fdr > 0)

plot(
  log10(plot_df$total + 1),
  -log10(plot_df$fdr),
  pch = 16,
  cex = 0.3,
  col = rgb(0,0,0,0.3),
  xlab = "Log10 total UMIs",
  ylab = "-log10 FDR"
)

abline(h = -log10(0.01), col = "red", lty = 2)

#gex_counts <- Read10X(raw_rna_dir)

# If Read10X returns a list, take Gene Expression
if (is.list(gex_counts)) {
  gex_counts <- gex_counts[["Gene Expression"]]
}

raw_umi <- Matrix::colSums(gex_counts)

barcode_keep <- names(raw_umi)[raw_umi > 100]

length(barcode_keep)

writeLines(
  barcode_keep,
  file.path(out_dir, paste0("barcodes_raw_", sample_name, "_RNAumi_gt100.tsv"))
)

summary(raw_umi)

