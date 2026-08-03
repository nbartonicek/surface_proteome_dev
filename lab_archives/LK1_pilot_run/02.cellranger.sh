#!/bin/bash
# ------------------------------------------------------------------
# LK1 pilot run - step 02 of 19
#
# Cell Ranger count on the GEX library. Produces the filtered_feature_bc_matrix
# that step 03 reads.
#
# Frozen for the lab archive 2026-08-01 from scripts/1.cellranger.sh (mtime 2026-06-01).
# md5 of the original: 3357bae744004d03ef9be94ff4d16758
# Body is unmodified - only this header was added, so the paths inside
# are still the ones that ran (relative to scripts/, i.e. ../results/...).
#
# NOTE: `RUN` and `SAMPLE_SHORT` in this copy point at LK2 (260528), because the
# script was re-pointed later. It ran on the pilot first, and the re-pointing was
# incomplete - `SAMPLE` is still "LK1-GEX_S5". Set RUN/SAMPLE_SHORT back to the
# pilot run to reproduce. Note also that OUT_DIR here is results/cellranger_withbam,
# whereas the pilot wrote to results/cellranger, which step 03 reads.
#
# Cell Ranger 10.0.0 from ~/tools/cellranger-10.0.0/cellranger. The commented
# `module load cellranger/7.2.0-gcc-13.2.0` line is stale and was not used.
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

#module load cellranger/7.2.0-gcc-13.2.0 

# ----------------------------
# Paths & parameters
# ----------------------------

RUN="260528_VH01624_464_222K7VKNX"
SAMPLE_SHORT="LK2-GEX"
PROJ="/scratch/users/nbartonicek/projects/amgen"
RAW_FASTQ_DIR="${PROJ}/raw/$RUN/Sample_$SAMPLE_SHORT"
WORK_FASTQ_DIR="${PROJ}/cellranger/fastq/$RUN"
OUT_DIR="${PROJ}/results/cellranger_withbam/$RUN"
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