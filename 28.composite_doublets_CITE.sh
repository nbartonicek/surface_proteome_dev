#!/bin/bash
#SBATCH -J composite_doublets
#SBATCH --partition=rhel_long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=FAIL

set -euo pipefail

# sbatch wrapper around 28.composite_doublets_CITE.py - runs COMPOSITE
# (sccomposite) doublet calling on the CITE-seq ADT modality plus RNA for one
# sample, and benchmarks it against the HTO / vireo / scDblFinder calls already
# in results_nf.
#
# Needs a python env with: sccomposite==1.0.0, torch, scipy, scikit-learn,
# pandas, matplotlib. COMPOSITE will use a GPU if it finds one but does not
# need one; everything below is sized for CPU.
#
# Memory/time notes:
#  * RNA counts are streamed out of the 456 MB gzipped CellRanger raw matrix in
#    pandas chunks rather than mmread'ing the lot, so the read is memory-light
#    but slow - budget ~10-20 min for that step alone.
#  * COMPOSITE's select_stable() densifies whatever matrix it is given. The
#    python script pre-filters RNA genes to those detected in >50% of cells
#    (exactly select_stable()'s own first step) so the dense array stays around
#    1 GB rather than ~10 GB, but the autograd graph over 33k cells still wants
#    room - hence 64G. Drop to 32G if you use --subsample.
#
# ADT_MODE:
#   raw                 raw ADT UMI counts - the paper's validated input
#   ambient_subtracted  crude per-marker background subtraction first
# DSB-normalised ADT is deliberately not an option; see the docstring in
# 28.composite_doublets_CITE.py for why it cannot work.
#
# Override at submit time, e.g.:
#   sbatch --export=ALL,ADT_MODE=ambient_subtracted 28.composite_doublets_CITE.sh
#   sbatch --export=ALL,SUBSAMPLE=3000 28.composite_doublets_CITE.sh   # smoke test

RUN_ID="${RUN_ID:-260528_VH01624_464_222K7VKNX}"
SAMPLE="${SAMPLE:-LK2-GEX}"
SAMPLE_SHORT="${SAMPLE_SHORT:-LK2}"
RESULTS_NF="${RESULTS_NF:-/scratch/users/nbartonicek/projects/amgen/results_nf}"
OUTROOT="${OUTROOT:-/scratch/users/nbartonicek/projects/amgen/results/composite_doublets}"
ADT_MODE="${ADT_MODE:-raw}"
# scdblfinder = the same droplets scDblFinder scored, from the same CellRanger
# raw counts, so the two are directly comparable. Do not change this unless you
# specifically want the post-scDblFinder demux subset.
CELL_SET="${CELL_SET:-scdblfinder}"
EXCLUDE_FEATURES="${EXCLUDE_FEATURES:-}"
SUBSAMPLE="${SUBSAMPLE:-0}"

mkdir -p logs

# ----------------------------
# Activate environment
# ----------------------------

# Example:
# mamba activate <python env with sccomposite + torch>

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "Run:          ${RUN_ID}"
echo "Sample:       ${SAMPLE} (${SAMPLE_SHORT})"
echo "results_nf:   ${RESULTS_NF}"
echo "Output root:  ${OUTROOT}"
echo "ADT mode:     ${ADT_MODE}"
echo "Cell set:     ${CELL_SET}"
echo "Exclude:      ${EXCLUDE_FEATURES:-none}"
echo "Subsample:    ${SUBSAMPLE}"
echo "Threads:      ${OMP_NUM_THREADS}"

if ! python -c "import sccomposite" 2>/dev/null; then
  echo "ERROR: sccomposite not importable in the active env." >&2
  echo "       pip install sccomposite==1.0.0" >&2
  exit 1
fi

ARGS=(
  --run-id "${RUN_ID}"
  --sample "${SAMPLE}"
  --sample-short "${SAMPLE_SHORT}"
  --results-nf "${RESULTS_NF}"
  --outroot "${OUTROOT}"
  --adt-mode "${ADT_MODE}"
  --cell-set "${CELL_SET}"
)

if [[ -n "${EXCLUDE_FEATURES}" ]]; then
  # shellcheck disable=SC2206
  ARGS+=(--exclude-features ${EXCLUDE_FEATURES//,/ })
fi

if [[ "${SUBSAMPLE}" -gt 0 ]]; then
  ARGS+=(--subsample "${SUBSAMPLE}")
fi

python 28.composite_doublets_CITE.py "${ARGS[@]}"

echo "Done."
