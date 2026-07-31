#!/bin/bash
#SBATCH --job-name=cellbender
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=rhel_gpu
#SBATCH --gres=gpu:1
#SBATCH --output=logs/cellbender_%A_%a.out
#SBATCH --error=logs/cellbender_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail

PROJ="/scratch/users/nbartonicek/projects/amgen"
CELLRANGER_DIR="${PROJ}/results/cellranger_withbam"
OUT_BASE="${PROJ}/results/cellbender"

RUNS=(
  "260522_VH01624_461_222JLJVNX/LK1-GEX"
  "260528_VH01624_464_222K7VKNX/LK2-GEX"
)

RUN_SAMPLE="${RUNS[$SLURM_ARRAY_TASK_ID]}"
RUN=$(dirname "${RUN_SAMPLE}")
SAMPLE=$(basename "${RUN_SAMPLE}")

INPUT="${CELLRANGER_DIR}/${RUN_SAMPLE}/outs/raw_feature_bc_matrix.h5"
METRICS="${CELLRANGER_DIR}/${RUN_SAMPLE}/outs/metrics_summary.csv"
OUTDIR="${OUT_BASE}/${RUN}/${SAMPLE}"

mkdir -p logs "${OUTDIR}"

[[ -f "${INPUT}" ]] || { echo "Missing input: ${INPUT}"; exit 1; }
[[ -f "${METRICS}" ]] || { echo "Missing metrics: ${METRICS}"; exit 1; }

EXPECTED_CELLS=$(python3 -c "
import csv, sys
with open(sys.argv[1]) as f:
    r = csv.reader(f)
    next(r)
    row = next(r)
    print(row[0].replace(',', ''))
" "${METRICS}")
TOTAL_DROPLETS=$(( EXPECTED_CELLS + 15000 ))

echo "Run: ${RUN}"
echo "Sample: ${SAMPLE}"
echo "Input: ${INPUT}"
echo "Expected cells: ${EXPECTED_CELLS}"
echo "Total droplets: ${TOTAL_DROPLETS}"
echo "Output dir: ${OUTDIR}"

export PATH=/home/nbartonicek/miniconda3/bin:$PATH
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate cellbender

cellbender remove-background \
    --input "${INPUT}" \
    --output "${OUTDIR}/${SAMPLE}_cellbender.h5" \
    --cuda \
    --expected-cells "${EXPECTED_CELLS}" \
    --total-droplets-included "${TOTAL_DROPLETS}" \
    --epochs 150 \
    --fpr 0.01 \
    --learning-rate 0.00005
