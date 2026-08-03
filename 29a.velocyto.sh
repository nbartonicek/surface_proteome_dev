#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH --partition=rhel_long
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX/velocity/logs/velocyto_%a_%j.log
#SBATCH --array=0-0  # one element per GEX sample; extend range if more samples added

RESULTS="/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX"
GTF="/scratch/users/nbartonicek/projects/amgen/annotation/refdata-gex-GRCh38-2024-A/genes/genes.gtf"

# List all GEX sample directories under 01_cellranger
SAMPLES=($(ls -d "${RESULTS}/01_cellranger"/*/  | xargs -n1 basename))
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

CR_OUTS="${RESULTS}/01_cellranger/${SAMPLE}/${SAMPLE}/outs"
OUT_DIR="${RESULTS}/velocity/${SAMPLE}"

mkdir -p "${OUT_DIR}" "${RESULTS}/velocity/logs"

echo "Processing sample: ${SAMPLE}"

module load samtools
source activate velocity

# Extract QC-passing barcodes from pipeline final metadata.
# The final metadata covers all donors in this GEX run — velocyto will only
# count spliced/unspliced for these barcodes, discarding empty droplets,
# doublets, and mapping-QC failures. Cell identity (donor, cell type, UMAP)
# is stored in the metadata itself and joined in 29b/29c.
FINAL_META="${RESULTS}/11_final/${SAMPLE}/tables"
BARCODES_FILE="${OUT_DIR}/qc_pass_barcodes.txt"

# Extract barcodes from the first column, skip header, filter mapping_error_QC == Pass
awk -F',' 'NR==1 { for(i=1;i<=NF;i++) { if($i=="\"barcode\"" || $i=="barcode" || i==1) bc=i; if($i=="\"mapping_error_QC\"" || $i=="mapping_error_QC") qc=i } next }
           $qc=="\"Pass\"" || $qc=="Pass" { gsub(/"/, "", $bc); print $bc }' \
    "${FINAL_META}"/*_final_metadata.csv > "${BARCODES_FILE}"

N_BARCODES=$(wc -l < "${BARCODES_FILE}")
echo "Using ${N_BARCODES} QC-passing barcodes from $(ls ${FINAL_META}/*_final_metadata.csv)"

BAM="${CR_OUTS}/possorted_genome_bam.bam"
SORTED_BAM="${OUT_DIR}/cellsorted_possorted_genome_bam.bam"

# Sort BAM by cell barcode (required by velocyto)
samtools sort -l 7 -m 48G -t CB -O BAM -@ "${SLURM_CPUS_PER_TASK}" \
    -o "${SORTED_BAM}" \
    "${BAM}"

# Run velocyto with filtered barcodes — only counts cells in our QC-passing set
velocyto run10x \
    --samtools-threads "${SLURM_CPUS_PER_TASK}" \
    --samtools-memory 8000 \
    -b "${BARCODES_FILE}" \
    "${CR_OUTS}/.." \
    "${GTF}"

echo "Done. Loom: ${CR_OUTS}/velocyto/${SAMPLE}.loom"
