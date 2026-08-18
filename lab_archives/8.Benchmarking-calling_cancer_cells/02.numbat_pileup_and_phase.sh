#!/bin/bash

# ------------------------------------------------------------------
# Calling cancer cells - step 02
#
# Numbat stage 1: pileup and phasing. Runs on the cluster, one array task per
# donor, inside the Numbat container.
#
# Collated from 9a.numbat.sh. The original had both runs' donor lists in the
# file with one of them commented out, so which run it processed depended on
# which block had been uncommented last. RUN and the donor list are now
# arguments, and the two runs' lists are kept side by side rather than one
# being commented out.
#
# This reads the Cell Ranger BAM, so it needs the cellranger_withbam output,
# and the per-donor barcode lists written by the annotation step.
#
# Submit from the scripts/ directory on the cluster:
#
#   sbatch lab_archives/8.Benchmarking-calling_cancer_cells/02.numbat_pileup_and_phase.sh LK2
#   sbatch lab_archives/8.Benchmarking-calling_cancer_cells/02.numbat_pileup_and_phase.sh LK1
#
# Stage 2 is step 03, which needs this step's allele counts.
# ------------------------------------------------------------------

#SBATCH --job-name=numbat_pileup
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --partition=rhel_long
#SBATCH --output=logs/numbat_pileup_%A_%a.out
#SBATCH --error=logs/numbat_pileup_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail

# ----------------------------
# Which run
# ----------------------------

SAMPLE_SHORT="${1:-LK2}"

case "${SAMPLE_SHORT}" in
  LK1)
    RUN="260423_VH01624_453_222HWMYNX"
    SAMPLES=("HBDN206-MNpCT" "HBDN392-AML-MDS" "HBDN501-AML-KMT2A" "MOLM13")
    ;;
  LK2)
    RUN="260528_VH01624_464_222K7VKNX"
    SAMPLES=("HBDN498-TP53" "APOP576-TP53" "HBDN376-TP53" "normal-01")
    ;;
  *)
    echo "Unknown run: ${SAMPLE_SHORT} (expected LK1 or LK2)"; exit 1
    ;;
esac

SAMPLE="${SAMPLE_SHORT}-GEX"

# ----------------------------
# Paths
# ----------------------------

PROJ="/scratch/users/nbartonicek/projects/amgen"

SIF="${PROJ}/containers/numbat-rbase_latest.sif"

CELLRANGER_DIR="${PROJ}/results/cellranger_withbam/${RUN}/${SAMPLE}/outs"
BAM="${CELLRANGER_DIR}/possorted_genome_bam.bam"

BARCODE_DIR="${PROJ}/results/seurat_demux/${RUN}/numbat/barcodes_by_sample_name"
OUT_BASE="${PROJ}/results/seurat_demux/${RUN}/numbat"

SNPVCF="${PROJ}/annotation/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.sorted.vcf.gz"

# These live inside the official Numbat container
GMAP="/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
PANELDIR="/data/1000G_hg38"
EAGLE="eagle"

DONOR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
LABEL="${SAMPLE_SHORT}_${DONOR}"

BARCODES="${BARCODE_DIR}/${DONOR}_barcodes.tsv"
OUTDIR="${OUT_BASE}/${LABEL}"

mkdir -p logs "${OUTDIR}"

# ----------------------------
# Checks
# ----------------------------

echo "RUN:      ${RUN}"
echo "SAMPLE:   ${SAMPLE}"
echo "DONOR:    ${DONOR}"
echo "LABEL:    ${LABEL}"
echo "BAM:      ${BAM}"
echo "BARCODES: ${BARCODES}"
echo "SNPVCF:   ${SNPVCF}"
echo "OUTDIR:   ${OUTDIR}"

[[ -f "${SIF}" ]]        || { echo "Missing SIF: ${SIF}"; exit 1; }
[[ -f "${BAM}" ]]        || { echo "Missing BAM: ${BAM}"; exit 1; }
[[ -f "${BAM}.bai" ]]    || { echo "Missing BAM index: ${BAM}.bai"; exit 1; }
[[ -f "${BARCODES}" ]]   || { echo "Missing barcodes: ${BARCODES}"; exit 1; }
[[ -f "${SNPVCF}" ]]     || { echo "Missing SNP VCF: ${SNPVCF}"; exit 1; }
[[ -f "${SNPVCF}.tbi" ]] || { echo "Missing SNP VCF index: ${SNPVCF}.tbi"; exit 1; }

# ----------------------------
# Pileup and phasing
# ----------------------------
# R_LIBS is unset so the container does not pick up the host library, which
# has a different Numbat build.

unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE
ulimit -s unlimited

apptainer exec \
  --cleanenv \
  --containall \
  --bind "${PROJ}:${PROJ}" \
  "${SIF}" \
  Rscript /numbat/inst/bin/pileup_and_phase.R \
    --label "${LABEL}" \
    --samples "${LABEL}" \
    --bams "${BAM}" \
    --barcodes "${BARCODES}" \
    --gmap "${GMAP}" \
    --eagle "${EAGLE}" \
    --snpvcf "${SNPVCF}" \
    --paneldir "${PANELDIR}" \
    --outdir "${OUTDIR}" \
    --ncores "${SLURM_CPUS_PER_TASK}" \
    --UMItag UB \
    --cellTAG CB

echo "Done."
echo "Expected allele count file:"
echo "${OUTDIR}/${LABEL}_allele_counts.tsv.gz"
