#!/bin/bash
#SBATCH -J patient_qc_batch
#SBATCH --partition=rhel_long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

# sbatch wrapper around 22a.patient_qc_report_batch.R - renders all donor
# reports for one run unattended. Needs pandoc (via rmarkdown::render()),
# which is why this runs on the cluster rather than being knit locally.
#
# Edit RUN_ID (and EXCLUDE_DONORS if wanted) below, or override at submit
# time with:
#   sbatch --export=RUN_ID=<run>,EXCLUDE_DONORS=<csv> 22a.patient_qc_report_batch.sh

RUN_ID="${RUN_ID:-260528_VH01624_464_222K7VKNX}"
RESULTS_DIR="${RESULTS_DIR:-../results_nf}"
OUTPUT_DIR="${OUTPUT_DIR:-../results/patient_reports/${RUN_ID}}"
EXCLUDE_DONORS="${EXCLUDE_DONORS:-}"
FORCE_RERUN="${FORCE_RERUN:-false}"

mkdir -p logs

# ----------------------------
# Activate environment
# ----------------------------

# Example:
# mamba activate <your R env with rmarkdown, Seurat, etc. + pandoc on PATH>

echo "Run: ${RUN_ID}"
echo "Results dir: ${RESULTS_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Exclude donors: ${EXCLUDE_DONORS:-none}"
echo "Force rerun: ${FORCE_RERUN}"

Rscript 22a.patient_qc_report_batch.R \
  "${RUN_ID}" \
  "${RESULTS_DIR}" \
  "${OUTPUT_DIR}" \
  "${EXCLUDE_DONORS}" \
  "${FORCE_RERUN}"

echo "Done."
