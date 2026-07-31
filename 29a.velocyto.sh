#!/bin/bash
# 29a.velocyto.sh
# Run velocyto on the LK2 CellRanger output to generate spliced/unspliced loom file.
# Run inside sinteractive: sinteractive --time 0-12:00 -c 4 --mem 64G

module load samtools

mamba activate velocyto

RESULTS="/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX"
SAMPLE="LK2-GEX"
CR_OUTS="${RESULTS}/01_cellranger/${SAMPLE}/${SAMPLE}/outs"
GTF="/scratch/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
OUT_DIR="${RESULTS}/29_velocity/${SAMPLE}"

mkdir -p "${OUT_DIR}"

BAM="${CR_OUTS}/possorted_genome_bam.bam"
SORTED_BAM="${OUT_DIR}/cellsorted_possorted_genome_bam.bam"

# Sort BAM by cell barcode (required by velocyto)
samtools sort -l 7 -m 48G -t CB -O BAM -@ 4 \
    -o "${SORTED_BAM}" \
    "${BAM}"

# Run velocyto — outputs a .loom file into ${CR_OUTS}/velocyto/
velocyto run10x \
    --samtools-threads 4 \
    --samtools-memory 8000 \
    "${CR_OUTS}/.." \
    "${GTF}"

echo "Done. Loom file: ${CR_OUTS}/velocyto/${SAMPLE}.loom"
