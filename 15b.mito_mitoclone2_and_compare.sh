#!/bin/bash
#SBATCH -J mito_mitoclone2
#SBATCH --partition=rhel_long
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

# Stage 2 of 2 - runs 15b.mito_mitoclone2_and_compare.R under r_env.
# Needs: mitoClone2, samtools, R with the shared tidyverse packages in
# 15.mito_clone_common.R. Does NOT need the mgatk CLI or Signac - it only
# reads mgatk's CSV output from disk. Run 15a.mito_mgatk_stage.sh (under
# the mgatk env) first, or this donor's comparison just won't have mgatk
# columns (skipped, not an error).
#
# Donors run PARALLEL_DONORS-at-a-time (parallel::mclapply), each
# mitoClone2/samtools call using NCORES_PER_DONOR cores - keep
# --cpus-per-task >= PARALLEL_DONORS * NCORES_PER_DONOR (default 16 = 4*4),
# and bump --mem if you raise PARALLEL_DONORS since donors' runs then
# overlap in memory.

RUN_ID="${RUN_ID:-260528_VH01624_464_222K7VKNX}"
SR_SAMPLE="${SR_SAMPLE:-LK2-GEX}"
RESULTS_DIR="${RESULTS_DIR:-/scratch/users/nbartonicek/projects/amgen/results_nf}"
EXCLUDE_DONORS="${EXCLUDE_DONORS:-}"
PARALLEL_DONORS="${PARALLEL_DONORS:-4}"
NCORES_PER_DONOR="${NCORES_PER_DONOR:-4}"

mkdir -p logs

#mamba activate r_env

# mitoClone2's own BAM-filtering step shells out to samtools - make sure
# it's on PATH regardless of what r_env itself bundles.
module load samtools 2>/dev/null || true

echo "Run: ${RUN_ID}  |  SR sample: ${SR_SAMPLE}  |  Results dir: ${RESULTS_DIR}"
echo "Parallel donors: ${PARALLEL_DONORS}  |  Cores per donor: ${NCORES_PER_DONOR}"

Rscript 15b.mito_mitoclone2_and_compare.R "${RUN_ID}" "${SR_SAMPLE}" "${RESULTS_DIR}" "${EXCLUDE_DONORS}" "${PARALLEL_DONORS}" "${NCORES_PER_DONOR}"

echo "Done."
