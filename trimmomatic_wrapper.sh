#!/bin/bash -l

#SBATCH -A snic2022-22-1221
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user alma.hansen.3759@student.uu.se

module load bioinfo-tools
module load trimmomatic

# Your commands
bash scripts/trimmomatic_run.sh data/VJ3440/221111/concatenated
