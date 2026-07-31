#!/bin/bash
#SBATCH --job-name=numbat_final_LK1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --partition=rhel_long
#SBATCH --output=logs/numbat_final_LK1_%A_%a.out
#SBATCH --error=logs/numbat_final_LK1_%A_%a.err
#SBATCH --array=0-3

#set -euo pipefail

PROJ="/scratch/users/nbartonicek/projects/amgen"
SIF="${PROJ}/containers/numbat-rbase_latest.sif"

SAMPLES=(
#  "HBDN206-MNpCT"
  "HBDN392-AML-MDS"
#  "HBDN501-AML-KMT2A"
)
SAMPLES=(
  "HBDN498-TP53"
  "APOP576-TP53"
  "HBDN376-TP53"
  "normal-01"
)

DONOR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

mkdir -p logs

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

apptainer exec \
  --cleanenv \
  --containall \
  --bind "${PROJ}:${PROJ}" \
  "${SIF}" \
  bash -lc "ulimit -s unlimited || true;
  Rscript "${PROJ}/scripts/9a.numbat.R" "${DONOR}"
  "
