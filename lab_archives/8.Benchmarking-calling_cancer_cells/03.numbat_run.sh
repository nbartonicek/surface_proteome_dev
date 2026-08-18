#!/bin/bash

# ------------------------------------------------------------------
# Calling cancer cells - step 03
#
# Numbat stage 2: the actual CNV and clone inference, one array task per donor,
# inside the container. Takes the allele counts from step 02 and the expression
# matrix exported by the annotation step.
#
# Collated from 9b.numbat.sh. Three things were wrong with the original:
#   - it declared SAMPLES twice, LK1 first and LK2 second, so the second
#     assignment silently overrode the first and the LK1 block was dead code
#     that looked live. RUN is now an argument.
#   - `set -euo pipefail` was commented out, so a failure inside the container
#     did not fail the array task and a broken donor looked like a finished one.
#     It is back on.
#   - the Rscript call was nested inside a double-quoted `bash -lc` string while
#     itself using double quotes, so the inner quotes closed the outer string.
#     It happened to work through word splitting, but any path with a space
#     would have broken it. Now single-quoted with the arguments passed
#     positionally.
#
# Submit from the scripts/ directory on the cluster:
#
#   sbatch lab_archives/8.Benchmarking-calling_cancer_cells/03.numbat_run.sh LK2
# ------------------------------------------------------------------

#SBATCH --job-name=numbat_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --partition=rhel_long
#SBATCH --output=logs/numbat_final_%A_%a.out
#SBATCH --error=logs/numbat_final_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail

SAMPLE_SHORT="${1:-LK2}"

case "${SAMPLE_SHORT}" in
  LK1) SAMPLES=("HBDN206-MNpCT" "HBDN392-AML-MDS" "HBDN501-AML-KMT2A" "MOLM13") ;;
  LK2) SAMPLES=("HBDN498-TP53" "APOP576-TP53" "HBDN376-TP53" "normal-01") ;;
  *)   echo "Unknown run: ${SAMPLE_SHORT} (expected LK1 or LK2)"; exit 1 ;;
esac

PROJ="/scratch/users/nbartonicek/projects/amgen"
SIF="${PROJ}/containers/numbat-rbase_latest.sif"
RSCRIPT="${PROJ}/scripts/lab_archives/8.Benchmarking-calling_cancer_cells/03a.numbat_run.R"

DONOR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

mkdir -p logs

echo "RUN:   ${SAMPLE_SHORT}"
echo "DONOR: ${DONOR}"

[[ -f "${SIF}" ]]     || { echo "Missing SIF: ${SIF}"; exit 1; }
[[ -f "${RSCRIPT}" ]] || { echo "Missing R script: ${RSCRIPT}"; exit 1; }

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

apptainer exec \
  --cleanenv \
  --containall \
  --bind "${PROJ}:${PROJ}" \
  "${SIF}" \
  bash -lc 'ulimit -s unlimited || true; Rscript "$1" "$2" "$3"' _ \
    "${RSCRIPT}" "${DONOR}" "${PROJ}"

echo "Done."
