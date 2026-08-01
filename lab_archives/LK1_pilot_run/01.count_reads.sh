#!/bin/bash
# ------------------------------------------------------------------
# LK1 pilot run - step 01 of 18
#
# Count R1 reads per library straight off the FASTQ. The denominator for everything else.
#
# Frozen for the lab archive 2026-07-31 from scripts/0.count_reads.sh (mtime 2026-05-05).
# md5 of the original: e2255232a77cca01925680c556f56fce
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

#SBATCH -J fastq_count
#SBATCH --partition=rhel_short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-2
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err

set -euo pipefail

# ----------------------------
# Input directory
# ----------------------------
BASE_DIR="/scratch/users/nbartonicek/projects/amgen/raw/260423_VH01624_453_222HWMYNX"

# Find all R1 files
FILES=(
  "${BASE_DIR}/Sample_LK1-GEX/"*R1*.fastq.gz
  "${BASE_DIR}/Sample_LK1-TotalSeqA-ADT/"*R1*.fastq.gz
  "${BASE_DIR}/Sample_LK1-TotalSeqA-HTO/"*R1*.fastq.gz
)

# Select file for this task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

OUT="${BASE_DIR}/fastq_read_counts.tsv"

echo "Processing: $FILE"

# Count lines → divide by 4 = reads
READS=$(zcat "$FILE" | wc -l)
READS=$((READS / 4))

# Write output (append safely)
echo -e "$(basename "$FILE")\t$READS" >> "$OUT"
