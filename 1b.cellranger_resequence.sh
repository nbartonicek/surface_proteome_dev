#!/bin/bash
#SBATCH -J cr_downsample_LK1
#SBATCH --partition=rhel_long
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=1-7
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nenad.bartonicek@petermac.org

#module load cellranger/10.0.0
module load seqtk

RUN="260528_VH01624_464_222K7VKNX"
SAMPLE_SHORT="LK2-GEX"
PROJ="/scratch/users/nbartonicek/projects/amgen"

RAW_FASTQ_DIR="${PROJ}/raw/${RUN}/Sample_${SAMPLE_SHORT}"
BASE_FASTQ_DIR="${PROJ}/cellranger_reseq/downsampled_fastq/${RUN}"
OUT_DIR="${PROJ}/results/cellranger_downsampled_reseq/${RUN}"
TRANSCRIPTOME="${PROJ}/annotation/refdata-gex-GRCh38-2024-A"

CORES=16
MEM_GB=128
SEED=123

#PROPS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
PROPS=(0.9 0.8 0.7 0.6 0.5 0.4 0.3)
PROP=${PROPS[$SLURM_ARRAY_TASK_ID-1]}
PROP_LABEL=$(printf "%03d" $(echo "$PROP * 100" | bc | cut -d. -f1))

FASTQ_DIR="${BASE_FASTQ_DIR}/frac_${PROP_LABEL}"
RUN_ID="${SAMPLE_SHORT}_frac_${PROP_LABEL}"

mkdir -p "${FASTQ_DIR}" "${OUT_DIR}" logs

R1_IN="${RAW_FASTQ_DIR}/${SAMPLE_SHORT}_S2_R1_001.fastq.gz"
R2_IN="${RAW_FASTQ_DIR}/${SAMPLE_SHORT}_S2_R2_001.fastq.gz"

R1_OUT="${FASTQ_DIR}/${SAMPLE_SHORT}_S00_R1_001.fastq.gz"
R2_OUT="${FASTQ_DIR}/${SAMPLE_SHORT}_S00_R2_001.fastq.gz"

if [[ "$PROP" == "1.0" ]]; then
  ln -sf "$R1_IN" "$R1_OUT"
  ln -sf "$R2_IN" "$R2_OUT"
else
  seqtk sample -s${SEED} "$R1_IN" "$PROP" | gzip > "$R1_OUT"
  seqtk sample -s${SEED} "$R2_IN" "$PROP" | gzip > "$R2_OUT"
fi

cd "${OUT_DIR}"

cellranger count \
  --id="${RUN_ID}" \
  --create-bam=true \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE_SHORT}" \
  --transcriptome="${TRANSCRIPTOME}" \
  --localcores="${CORES}" \
  --localmem="${MEM_GB}"