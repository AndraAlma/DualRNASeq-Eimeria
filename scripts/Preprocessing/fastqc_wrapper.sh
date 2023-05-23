#!/bin/bash -l

#SBATCH -A snic2022-22-1221
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J fastqc_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user alma.hansen.3759@student.uu.se

module load bioinfo-tools
module load FastQC/0.11.8
module load MultiQC/1.8

# Your commands
bash scripts/fastqc_run.sh results/trimmomatic/VJ3440/221111


