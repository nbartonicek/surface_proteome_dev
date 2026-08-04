#!/bin/bash
# ------------------------------------------------------------------
# Demultiplexing benchmark: genotyping vs HTO - step 06 of 7
#
# Verifies the unhashed-donor HTO/vireo matching fix by running the production nextflow script against LK3's completed cellSNP/vireo and HTO output.
#
# Explicitly a trial/diagnostic. It runs the REAL production script
# (nextflow/bin/process_cellsnp_hto.R) with the arguments the pipeline would
# pass, rather than a reimplementation, so it exercises the code that ships.
#
# Frozen for the lab archive 2026-08-03 from scripts/25.trial_test_unhashed_donor_demux.sh (mtime 2026-07-21).
# md5 of the original: 0d72bf91e133e51b1fb9b5d66a0fa8c4
# Body is unmodified - only this header was added.
# ------------------------------------------------------------------

#
# *** TRIAL / DIAGNOSTIC SCRIPT - NOT PART OF THE NEXTFLOW PIPELINE ***
#
# Lives here in scripts/ (not in nextflow/) deliberately - testing/trial
# scripts should never be placed inside the pipeline directory itself.
#
# Purpose: verify the unhashed-donor HTO/vireo matching fix (in
# nextflow/bin/process_cellsnp_hto.R) actually works, using LK3's real,
# already-completed cellSNP/vireo + HTO output - without waiting hours for
# another full nextflow run through CELLRANGER/CellBender/cellSNP-Vireo
# again.
#
# Runs the REAL production script (nextflow/bin/process_cellsnp_hto.R)
# directly with the exact same arguments the nextflow module would pass for
# LK3-GEX - this is not a reimplementation, it exercises the actual code
# that will run in the pipeline.
#
# Writes to scripts/TRIAL_output_unhashed_donor_demux/ only - never touches
# results_nf/ or nextflow/work/. Safe to re-run repeatedly; delete that
# directory when done experimenting.
#
# Usage - self-contained, activates the correct R env (r_env, same one the
# real PROCESS_CELLSNP_HTO uses) internally, no pre-activation needed:
#   cd /scratch/users/nbartonicek/projects/amgen/scripts
#   ./25.trial_test_unhashed_donor_demux.sh

set -euo pipefail

PROJ="/scratch/users/nbartonicek/projects/amgen"
NF_DIR="${PROJ}/nextflow"
SCRIPTS_DIR="${PROJ}/scripts"
RUN="260717_VH01624_477_222KG22NX"
RESULTS="${PROJ}/results_nf/${RUN}"

# Real files already on disk from gigantic_keller's run (verified present
# before writing this script - not guessed paths):
SCDBLFINDER_RDS="${RESULTS}/rds/05_scdblfinder/LK3-GEX/rds/seurat_emptyDrops_RNA_scDblFinder_filtered.rds"
VIREO_DONOR_IDS="${RESULTS}/03_genotype_demux/LK3-GEX/donor_ids.tsv"
HTO_DIR="${RESULTS}/04_cite_seq_count/LK3-GEX/hto_counts_LK3_emptydrops"

SAMPLE_SHORT="LK3"
SAMPLE_NAME="LK3-GEX"
POSITIVE_QUANTILE="0.99"
DONOR_IDS="normal_CD34-01;;;normal-02"
UNHASHED="WEI21_26-17_NK"

TRIAL_DIR="${SCRIPTS_DIR}/TRIAL_output_unhashed_donor_demux"
FIGURES_DIR="${TRIAL_DIR}/figures"
TABLES_DIR="${TRIAL_DIR}/tables"
RDS_DIR="${TRIAL_DIR}/rds"

echo "=============================================================="
echo " TRIAL RUN - unhashed-donor HTO/vireo matching fix"
echo " This is a diagnostic script, NOT the real pipeline."
echo " Testing the real script at: ${NF_DIR}/bin/process_cellsnp_hto.R"
echo " Output directory: ${TRIAL_DIR}"
echo "=============================================================="
echo ""

echo "--- Checking required inputs exist ---"
missing=0
for f in "$SCDBLFINDER_RDS" "$VIREO_DONOR_IDS" "${HTO_DIR}/umi_count/matrix.mtx.gz"; do
  if [[ ! -f "$f" ]]; then
    echo "  MISSING: $f"
    missing=1
  else
    echo "  found:   $f"
  fi
done
if [[ "$missing" -eq 1 ]]; then
  echo ""
  echo "TRIAL ABORTED: one or more required inputs are missing. Check that"
  echo "gigantic_keller's outputs (or a later successful run's) are still"
  echo "present at the paths above before re-running this trial."
  exit 1
fi

mkdir -p "$FIGURES_DIR" "$TABLES_DIR" "$RDS_DIR"

# Activate the exact same R environment PROCESS_CELLSNP_HTO uses in the real
# pipeline (see nextflow.config's withName: 'PROCESS_CELLSNP_HTO|CITE_QC'
# beforeScript) - not just whatever env happens to be active in the calling
# shell, so this trial genuinely mirrors the real pipeline's execution.
R_ENV="/home/nbartonicek/miniconda3/envs/r_env"
export PATH="/home/nbartonicek/miniconda3/bin:$PATH"
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "$R_ENV"

echo ""
echo "--- Running nextflow/bin/process_cellsnp_hto.R (the real script) ---"
echo "    (using R env: ${R_ENV})"
Rscript "${NF_DIR}/bin/process_cellsnp_hto.R" \
    "$SCDBLFINDER_RDS" \
    "$VIREO_DONOR_IDS" \
    "${HTO_DIR}/umi_count" \
    "$SAMPLE_SHORT" \
    "$SAMPLE_NAME" \
    "$POSITIVE_QUANTILE" \
    "$FIGURES_DIR" \
    "$TABLES_DIR" \
    "$RDS_DIR" \
    "$DONOR_IDS" \
    "$UNHASHED"

echo ""
echo "=============================================================="
echo " TRIAL RESULTS"
echo "=============================================================="

MAP_CSV=$(ls "${TABLES_DIR}"/vireo_to_hto_sample_name_map_*.csv 2>/dev/null | head -1 || true)
if [[ -z "$MAP_CSV" ]]; then
  echo "FAIL: no vireo_to_hto_sample_name_map_*.csv was produced at all."
  echo "The script likely crashed before reaching the matching step - check"
  echo "the R output above for the actual error."
  exit 1
fi

echo ""
echo "vireo -> donor mapping (${MAP_CSV}):"
column -s, -t "$MAP_CSV"

echo ""
if grep -q "$UNHASHED" "$MAP_CSV"; then
  echo "PASS: unhashed donor '${UNHASHED}' was assigned to a vireo cluster by elimination."
else
  echo "FAIL: unhashed donor '${UNHASHED}' does NOT appear in the vireo->donor mapping."
  echo "The elimination logic did not resolve to exactly one unmatched vireo"
  echo "cluster - check the R script's warning messages above for details."
  exit 1
fi

SINGLETS_CSV=$(ls "${TABLES_DIR}"/cell_metadata_HTO_vireo_combined_*_singlets_HTO_named.csv 2>/dev/null | head -1 || true)
if [[ -n "$SINGLETS_CSV" ]]; then
  echo ""
  echo "Final sample_name counts (HTO+vireo singlets only):"
  awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_name") c=i} NR>1 && c{print $c}' "$SINGLETS_CSV" \
    | sort | uniq -c | sort -rn
  echo ""
  echo "Sanity check: expect roughly 3 groups here (normal_CD34-01,"
  echo "normal-02, and WEI21_26-17_NK), each with a plausible non-trivial"
  echo "cell count - not one group dominating everything, and not a 4th"
  echo "or missing group."
fi

echo ""
echo "=============================================================="
echo " TRIAL complete. This directory (${TRIAL_DIR}) is scratch output"
echo " for this diagnostic only - delete it whenever you're done. It does"
echo " NOT feed into or affect the real nextflow pipeline in any way."
echo "=============================================================="
