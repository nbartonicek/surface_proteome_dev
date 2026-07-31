#!/bin/bash
#SBATCH --job-name=normal_aml
#SBATCH --partition=rhel_long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err
#SBATCH --mail-type=ALL

#set -euo pipefail

#source ~/.bashrc

# activate your conda environment
#mamba activate scanpy
# or:
# conda activate scanpy

#cd /scratch/users/nbartonicek/projects/amgen/scripts

python normal_aml.py