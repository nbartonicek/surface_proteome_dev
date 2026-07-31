#!/bin/bash
#SBATCH --job-name=cellsnp_vireo_LK1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=rhel_long                # use long queue (CellRanger is heavy)
#SBATCH --output=logs/cellsnp_vireo_LK2.%j.out
#SBATCH --error=logs/cellsnp_vireo_LK2.%j.err
#SBATCH --mail-type=ALL
#SBATCH -o logs/%j.out                       # stdout
#SBATCH -e logs/%j.err                       # stderr

#set -euo pipefail

module load samtools
# Or activate conda env containing cellsnp-lite and vireo
#mamba activate cellsnp_env

# ----------------------------
# Paths & parameters
# ----------------------------

RUN="260522_VH01624_461_222JLJVNX"
SAMPLE_SHORT="LK1"
SAMPLE="${SAMPLE_SHORT}-GEX"

PROJ="/scratch/users/nbartonicek/projects/amgen"

CELLRANGER_DIR="${PROJ}/results/cellranger_withbam/${RUN}/${SAMPLE}/outs"
BAM="${CELLRANGER_DIR}/possorted_genome_bam.bam"

BARCODES="${PROJ}/results/emptydrops/${RUN}/${SAMPLE}/emptydrops_barcodes.txt"

OUT_BASE="${PROJ}/results/genotype_demux/${RUN}/${SAMPLE}"
CELLSNP_OUT="${OUT_BASE}/cellsnp"
VIREO_OUT="${OUT_BASE}/vireo"

COMMON_SNP_VCF="${PROJ}/annotation/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.sorted.vcf.gz"

N_DONORS=4
N_THREADS=8

mkdir -p "${CELLSNP_OUT}" "${VIREO_OUT}" "${PROJ}/scripts/logs"

# ----------------------------
# Sanity checks
# ----------------------------

[[ -f "${BAM}" ]] || { echo "Missing BAM: ${BAM}"; }
[[ -f "${BAM}.bai" ]] || { echo "Missing BAM index: ${BAM}.bai"; }
[[ -f "${BARCODES}" ]] || { echo "Missing barcodes: ${BARCODES}"; }
[[ -f "${COMMON_SNP_VCF}" ]] || { echo "Missing SNP VCF: ${COMMON_SNP_VCF}"; }
[[ -f "${COMMON_SNP_VCF}.tbi" ]] || { echo "Missing SNP VCF index: ${COMMON_SNP_VCF}.tbi"; }

# ----------------------------
# cellSNP-lite
# ----------------------------

cellsnp-lite \
  -s "${BAM}" \
  -b "${BARCODES}" \
  -O "${CELLSNP_OUT}" \
  -R "${COMMON_SNP_VCF}" \
  -p "${N_THREADS}" \
  --minCOUNT 20 \
  --minMAF 0.10 \
  --cellTAG CB \
  --UMItag UB \
  --gzip

# ----------------------------
# vireo
# ----------------------------

vireo \
  -c "${CELLSNP_OUT}" \
  -N "${N_DONORS}" \
  -o "${VIREO_OUT}"

echo "Done."
echo "cellSNP-lite output: ${CELLSNP_OUT}"
echo "vireo output: ${VIREO_OUT}"