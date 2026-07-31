#!/bin/bash
#SBATCH -J mito_mgatk
#SBATCH --partition=rhel_long
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

# Stage 1 of 2 - runs 15a.mito_mgatk_stage.R under the `mgatk` env.
# Needs: mgatk CLI, samtools, R with Matrix/Signac (+ shared tidyverse
# packages in 15.mito_clone_common.R). Does NOT need mitoClone2.
# Run 15b.mito_mitoclone2_and_compare.sh (under r_env) after this finishes.
#
# Donors run PARALLEL_DONORS-at-a-time (parallel::mclapply), each mgatk tenx
# call using NCORES_PER_DONOR cores - keep --cpus-per-task >=
# PARALLEL_DONORS * NCORES_PER_DONOR (default 16 = 4*4), and bump --mem if
# you raise PARALLEL_DONORS since donors' mgatk runs then overlap in memory.

RUN_ID="${RUN_ID:-260528_VH01624_464_222K7VKNX}"
SR_SAMPLE="${SR_SAMPLE:-LK2-GEX}"
RESULTS_DIR="${RESULTS_DIR:-/scratch/users/nbartonicek/projects/amgen/results_nf}"
EXCLUDE_DONORS="${EXCLUDE_DONORS:-}"
PARALLEL_DONORS="${PARALLEL_DONORS:-4}"
NCORES_PER_DONOR="${NCORES_PER_DONOR:-4}"

mkdir -p logs

#mamba activate mgatk

# See old 15.mgatk.sh for why this matters: mgatk's internal Snakemake
# workflow shells out via a bare `python` subprocess call resolved via
# PATH. Spack's lmod-managed python can end up ahead of the conda env's
# bin/ in PATH even after `mamba activate`, breaking pysam imports in
# those subprocesses. Forcing the active env's bin/ to the front fixes it.
if [[ -n "${CONDA_PREFIX:-}" ]]; then
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi

echo "Run: ${RUN_ID}  |  SR sample: ${SR_SAMPLE}  |  Results dir: ${RESULTS_DIR}"
echo "Parallel donors: ${PARALLEL_DONORS}  |  Cores per donor: ${NCORES_PER_DONOR}"

Rscript 15a.mito_mgatk_stage.R "${RUN_ID}" "${SR_SAMPLE}" "${RESULTS_DIR}" "${EXCLUDE_DONORS}" "${PARALLEL_DONORS}" "${NCORES_PER_DONOR}"

echo "Done."
