#!/bin/bash
# ------------------------------------------------------------------
# CITE-seq benchmark - step 01 of 07
#
# CITE-seq-Count on the ADT library, whitelisted to the emptyDrops barcodes. The production counting run.
#
# Frozen for the lab archive 2026-08-04 from scripts/4.cite_seq_count.sbatch (mtime 2026-06-03).
# md5 of the original: 3653e7d54eada6785f62c77a22b9aa1d
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

#SBATCH -J cite_seq_raw
#SBATCH --partition=rhel_long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

RUN="260522_VH01624_461_222JLJVNX"
SAMPLE_SHORT="LK1"
SAMPLE_GEX="${SAMPLE_SHORT}-GEX"
PROJ="/scratch/users/nbartonicek/projects/amgen"

ADT_FASTQ_DIR="${PROJ}/raw/${RUN}/Sample_${SAMPLE_SHORT}-TotalSeqA-ADT"
OUT_DIR="${PROJ}/results/cite_seq_count/${RUN}"
ADT_TAGS="${PROJ}/annotation/adt_tags.csv"

CORES=16

mkdir -p "${OUT_DIR}" logs

# ----------------------------
# Barcode whitelist to use
# ----------------------------

RAW_WL="${PROJ}/results/emptydrops/${RUN}/${SAMPLE_GEX}/barcodes_raw_${SAMPLE_GEX}_RNAumi_gt100.tsv"

echo "Using barcode whitelist:"
echo "${RAW_WL}"

if [[ ! -s "${RAW_WL}" ]]; then
  echo "ERROR: RAW_WL does not exist or is empty:"
  echo "${RAW_WL}"
  exit 1
fi

CELLS=$(wc -l < "${RAW_WL}")

echo "Number of barcodes/cells used for -cells:"
echo "${CELLS}"

# ----------------------------
# ADT FASTQs
# ----------------------------

R1=$(ls "${ADT_FASTQ_DIR}/${SAMPLE_SHORT}-TotalSeqA-ADT"_S*_R1_001.fastq.gz)
R2=$(ls "${ADT_FASTQ_DIR}/${SAMPLE_SHORT}-TotalSeqA-ADT"_S*_R2_001.fastq.gz)

echo "R1:"
echo "${R1}"
echo "R2:"
echo "${R2}"

# ----------------------------
# ADT counting
# ----------------------------

CITE-seq-Count \
  -R1 ${R1} \
  -R2 ${R2} \
  -t "${ADT_TAGS}" \
  -wl "${RAW_WL}" \
  -cbf 1 -cbl 16 \
  -umif 17 -umil 28 \
  -o "${OUT_DIR}/adt_counts_${SAMPLE_SHORT}_raw" \
  -T "${CORES}" \
  -cells "${CELLS}"

echo "Done."
echo "ADT raw output: ${OUT_DIR}/adt_counts_${SAMPLE_SHORT}_raw"