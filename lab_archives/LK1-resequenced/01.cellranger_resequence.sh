#!/bin/bash
# ------------------------------------------------------------------
# LK1 re-sequenced - step 01 of 5
#
# Cell Ranger 10.0.0 count on the re-sequenced LK1 GEX library.
#
# Frozen for the lab archive 2026-08-03 from scripts/1a.cellranger_resequence.sh (mtime 2026-05-25).
# md5 of the original: 7e14f4a563dcc673ae8b608e79caeb13
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
# ------------------------------------------------------------------

#SBATCH -J cellranger_LN         # job name
#SBATCH --partition=rhel_long                # use long queue (CellRanger is heavy)
#SBATCH --time=2-00:00:00                    # 2 days walltime (adjust if needed)
#SBATCH --cpus-per-task=16                   # must match Cell Ranger threads
#SBATCH --mem=128G                            # total memory
#SBATCH -o logs/%j.out                       # stdout
#SBATCH -e logs/%j.err                       # stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nenad.bartonicek@petermac.org   # change if needed

set -euo pipefail

#warning: cellranger-10.0.0 required 

# ----------------------------
# Paths & parameters
# ----------------------------

RUN="260522_VH01624_461_222JLJVNX"
SAMPLE_SHORT="LK1-GEX"
PROJ="/scratch/users/nbartonicek/projects/amgen"
RAW_FASTQ_DIR="${PROJ}/raw_reseq/$RUN/Sample_$SAMPLE_SHORT"
WORK_FASTQ_DIR="${PROJ}/cellranger_reseq/fastq/$RUN"
OUT_DIR="${PROJ}/results/cellranger_withbam_reseq/$RUN"
SAMPLE="LK1-GEX_S5"


TRANSCRIPTOME="${PROJ}/annotation/refdata-gex-GRCh38-2024-A"

CORES=16
MEM_GB=128

# ----------------------------
# Setup folders
# ----------------------------
mkdir -p "${WORK_FASTQ_DIR}" "${OUT_DIR}" logs

# ----------------------------
# Symlink FASTQs into Illumina naming
# ----------------------------
#ln -sf "${RAW_FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz" \
#       "${WORK_FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz"
#
#ln -sf "${RAW_FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz" \
#       "${WORK_FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

# ----------------------------
# Run Cell Ranger
# ----------------------------
cd "${OUT_DIR}"
#   --chemistry=SC3Pv3 \

cellranger count \
   --id="${SAMPLE_SHORT}" \
   --create-bam=true \
   --fastqs="${RAW_FASTQ_DIR}" \
   --sample="${SAMPLE_SHORT}" \
   --transcriptome="${TRANSCRIPTOME}" \
   --localcores="${CORES}" \
   --localmem="${MEM_GB}"