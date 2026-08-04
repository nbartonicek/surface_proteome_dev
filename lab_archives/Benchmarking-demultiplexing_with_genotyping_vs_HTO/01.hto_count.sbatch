#!/bin/bash
# ------------------------------------------------------------------
# Demultiplexing benchmark: genotyping vs HTO - step 01 of 7
#
# CITE-seq-Count on the HTO library, whitelisted to the emptyDrops barcodes.
#
# Frozen for the lab archive 2026-08-03 from scripts/3.hto_count.sbatch (mtime 2026-06-03).
# md5 of the original: 59aafb846a433612bdb08433f1433363
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

#SBATCH -J hto
#SBATCH --partition=rhel_long
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# ----------------------------
# Paths & parameters
# ----------------------------
#mamba activate cite_seq_count

RUN="260528_VH01624_464_222K7VKNX"
SAMPLE_SHORT="LK2"
SAMPLE_GEX="${SAMPLE_SHORT}-GEX"

PROJ="/scratch/users/nbartonicek/projects/amgen"

HTO_FASTQ_DIR="${PROJ}/raw/${RUN}/Sample_${SAMPLE_SHORT}-TotalSeqA-HTO"
OUT_DIR="${PROJ}/results/cite_seq_count/${RUN}"
EMPTYDROPS_DIR="${PROJ}/results/emptydrops/${RUN}/${SAMPLE_GEX}"

HTO_TAGS="${PROJ}/annotation/hto_tags_${SAMPLE_SHORT}.csv"

CORES=16
MEM_GB=128

# ----------------------------
# Input whitelist from EmptyDrops
# ----------------------------

WL="${EMPTYDROPS_DIR}/emptydrops_barcodes.txt"

if [[ ! -s "${WL}" ]]; then
  echo "ERROR: whitelist not found or empty:"
  echo "${WL}"
  exit 1
fi

CELLS=$(wc -l < "${WL}" | tr -d ' ')

if [[ "${CELLS}" -eq 0 ]]; then
  echo "ERROR: zero cells detected in whitelist:"
  echo "${WL}"
  exit 1
fi

echo "Using whitelist:"
echo "${WL}"
echo "Number of cells:"
echo "${CELLS}"

# ----------------------------
# Setup folders
# ----------------------------

mkdir -p "${OUT_DIR}" logs "${OUT_DIR}/hto_counts_${SAMPLE_SHORT}_emptydrops"

# ----------------------------
# FASTQs
# ----------------------------

R1=$(ls "${HTO_FASTQ_DIR}/${SAMPLE_SHORT}-TotalSeqA-HTO"_S*_R1_001.fastq.gz)
R2=$(ls "${HTO_FASTQ_DIR}/${SAMPLE_SHORT}-TotalSeqA-HTO"_S*_R2_001.fastq.gz)

echo "R1: ${R1}"
echo "R2: ${R2}"
echo "HTO tags: ${HTO_TAGS}"

# ----------------------------
# Run CITE-seq-Count: HTO
# ----------------------------

CITE-seq-Count \
  -R1 "${R1}" \
  -R2 "${R2}" \
  -t "${HTO_TAGS}" \
  -wl "${WL}" \
  -cbf 1 -cbl 16 \
  -umif 17 -umil 28 \
  -o "${OUT_DIR}/hto_counts_${SAMPLE_SHORT}_emptydrops" \
  -T "${CORES}" \
  -cells "${CELLS}"

echo "Done."
echo "Output written to:"
echo "${OUT_DIR}/hto_counts_${SAMPLE_SHORT}_emptydrops"